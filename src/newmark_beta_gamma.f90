#include "scalar.fpp"
!=====================================================================!
! Newmark-Beta-Gamma Integration Module for first and second order
! systems with adjoint derivative capabilities.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module nbg_integrator

  use integrator_class, only : integrator
  use physics_class,    only : physics

  implicit none

  private
  public :: NBG
  
  type(scalar), parameter :: ONE  = 1.0d0
  type(scalar), parameter :: ZERO = 0.0d0

  !===================================================================! 
  ! NBG Integrator type
  !===================================================================! 

  type, extends(integrator) :: NBG
     
     type(scalar) :: BETA   = 0.25d0
     type(scalar) :: GAMMA  = 0.50d0

   contains

     ! Destructor
     procedure, public  :: finalize
     
     ! Routines for integration
     procedure, public  :: integrate
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff

     ! Routines for adjoint gradient
     procedure, public  :: marchBackwards
     procedure, private :: assembleRHS
     procedure :: computeTotalDerivative

  end type NBG

  interface NBG
     module procedure initialize
  end interface NBG

contains


  !===================================================================!
  ! Initialize the NBG datatype and allocate required variables
  !===================================================================!
  
  type(nbg) function initialize( system, tinit, tfinal, h, second_order )  result (this)
   
    class(physics), target          :: system
    real(dp) , OPTIONAL, intent(in) :: tinit, tfinal
    real(dp) , OPTIONAL, intent(in) :: h
    logical  , OPTIONAL, intent(in) :: second_order

    print *, "======================================"
    print *, ">>   Adams Bashforth Moulton       << "
    print *, "======================================"

    !-----------------------------------------------------------------!
    ! Set the physical system in to the integrator                    !
    !-----------------------------------------------------------------!
    
    call this % setPhysicalSystem(system)

    !-----------------------------------------------------------------!
    ! Fetch the number of state variables from the system object
    !-----------------------------------------------------------------!

    this % nsvars = system % getNumStateVars()
    print '("  >> Number of variables    : ",i4)', this % nsvars
    
    if ( .not. (this % nsvars .gt. 0) ) then
       stop ">> Error: No state variable. Stopping."
    end if

    !-----------------------------------------------------------------!
    ! Set the order of the governing equations
    !-----------------------------------------------------------------!
    
    if (present(second_order)) then
       this % second_order = second_order
    end if
    print '("  >> Second order           : ",L1)', this % second_order

    !-----------------------------------------------------------------!
    ! Set the initial and final time
    !-----------------------------------------------------------------!

    if (present(tinit)) then
       this % tinit = tinit
    end if
    print '("  >> Start time             : ",F8.3)', this % tinit

    if (present(tfinal)) then
       this % tfinal = tfinal
    end if
    print '("  >> End time               : ",F8.3)', this % tfinal

    !-----------------------------------------------------------------!
    ! Set the user supplied initial step size
    !-----------------------------------------------------------------!

    if (present(h)) then
       this % h = h 
    end if
    print '("  >> Step size              : ",E9.3)', this % h
    
    !-----------------------------------------------------------------!
    ! Find the number of time steps required during integration
    !-----------------------------------------------------------------!

    this % num_steps = int((this % tfinal - this % tinit)/this % h) + 1 
    print '("  >> Number of steps        : ",i6)', this % num_steps

    !-----------------------------------------------------------------!
    ! Allocate space for the RHS of adjoint equations
    !-----------------------------------------------------------------!

    allocate(this % psi(this % num_steps, this % nsvars))
    this % psi = 0.0d0

    !-----------------------------------------------------------------!
    ! Allocate space for the global states and time
    !-----------------------------------------------------------------!

    allocate(this % time(this % num_steps))
    this % time = 0.0d0

    allocate(this % U(this % num_steps, this % nsvars))
    this % U = 0.0d0

    allocate(this % UDOT(this % num_steps, this % nsvars))
    this % UDOT = 0.0d0

    allocate(this % UDDOT(this % num_steps, this % nsvars))
    this % UDDOT = 0.0d0
    
    ! Set the start time
    this % time(1) = this % tinit

  end function initialize

  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!
  
  subroutine finalize( this )

    class(NBG) :: this

    ! Clear global states and time
    if(allocated(this % UDDOT)) deallocate(this % UDDOT)
    if(allocated(this % UDOT)) deallocate(this % UDOT)
    if(allocated(this % U)) deallocate(this % U)
    if(allocated(this % time)) deallocate(this % time)

    if(allocated(this % psi)) deallocate(this % psi)

  end subroutine finalize
  
  !===================================================================!
  ! Approximate the state variables at each step using NBG formula
  !===================================================================!
  
  subroutine approximateStates( this )

    class(NBG), intent(inout) :: this
    integer                   :: k
    type(scalar)              :: scale

    k = this % current_step

    !-----------------------------------------------------------------!
    ! Assume a UDDOT for the next time step
    !-----------------------------------------------------------------!

    if ( k .eq. 2 ) then
       this % uddot(k,:) = 1.0d0
    else
       this % uddot(k,:) = this % uddot(k-1,:) 
    end if

    !-----------------------------------------------------------------!
    ! Approximate UDOT using NBG
    !-----------------------------------------------------------------!

    this % udot(k,:) = this % udot(k-1,:) 

    scale = this % h * (1.0d0 - this % GAMMA)
    this % udot(k,:) = scale*this % uddot(k-1,:) 

    scale = this % h * this % GAMMA
    this % udot(k,:) = scale*this % uddot(k,:) 
    
    !-----------------------------------------------------------------!
    ! Approximate U using NBG
    !-----------------------------------------------------------------!

    this % u(k,:) = this % u(k-1,:) 

    scale = this % h
    this % u(k,:) = scale*this % udot(k-1,:) 

    scale = this % h * this % h * (1.0d0 - 2.0d0 * this % BETA)/2.0d0
    this % u(k,:) = scale*this % uddot(k-1,:) 

    scale = this % h * this % h * this % BETA
    this % u(k,:) = scale*this % uddot(k,:) 
    
  end subroutine approximateStates
  
  !===================================================================!
  ! Returns the linearization scalar coefficients: alpha, beta, gamma
  !===================================================================!
  
  subroutine getLinearCoeff( this, k, alpha, beta, gamma )

    class(NBG), intent(inout)   :: this
    integer, intent(in)         :: k
    type(scalar), intent(inout) :: alpha, beta, gamma
    integer :: m

    gamma = 1.0d0
    beta  = this % h * this % GAMMA
    alpha = this % h * this % h * this % BETA

  end subroutine getLinearCoeff

  !===================================================================!
  ! Time integration logic
  !===================================================================!

  subroutine Integrate( this )

    class(NBG)   :: this
    type(scalar) :: alpha, beta, gamma
    integer      :: k

    ! Set states to zeror
    this % U     = 0.0d0
    this % UDOT  = 0.0d0
    this % UDDOT = 0.0d0
    this % time  = 0.0d0

    ! Set the initial condition
    call this % system % getInitialStates(this % time(1), &
         & this % u(1,:), this % udot(1,:))

    this % current_step = 1

    ! March in time
    time: do k = 2, this % num_steps

       this % current_step =  k
       
       ! Increment the time (states are already advanced after the
       ! Newton solve)
       this % time(k) = this % time(k-1) + this % h
       
       ! Approximate the states u, udot and uddot using NBG stencil
       call this % approximateStates()

       ! Determine the coefficients for linearing the Residual
       call this % getLinearCoeff(k, alpha, beta, gamma)

       ! Solve the nonlinear system at each step by driving the
       ! residual to zero
       call this % newtonSolve(alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))

    end do time

  end subroutine Integrate
  
  !===================================================================!
  ! Subroutine that marches backwards in time to compute the lagrange
  ! multipliers (adjoint variables for the function)
  ! ===================================================================!
  
  subroutine marchBackwards( this )

    class(NBG)                :: this
    integer                   :: k
    type(scalar)              :: alpha, beta, gamma

    time: do k = this % num_steps, 2, -1

       this % current_step = k 

       !--------------------------------------------------------------!
       ! Determine the linearization coefficients for the Jacobian
       !--------------------------------------------------------------!

       call this % getLinearCoeff(k, alpha, beta, gamma)

       !--------------------------------------------------------------!
       ! Solve the adjoint equation at each step
       !--------------------------------------------------------------!

       call this % adjointSolve(this % psi(k,:), alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))

       print*, "k,psi=", k, this % psi(k,:)

    end do time

  end subroutine marchBackwards
  
  !===================================================================!
  ! Function that puts together the right hand side of the adjoint
  ! equation into the supplied rhs vector.
  !===================================================================!
  
  subroutine assembleRHS( this, rhs )

    class(NBG)                                :: this
    type(scalar), dimension(:), intent(inout) :: rhs
    type(scalar), dimension(:,:), allocatable :: jac
    type(scalar)                              :: scale
    integer                                   :: k, i, m1, m2, idx, m

    allocate( jac(this % nSVars, this % nSVars) )

    ! Zero the RHS first
    rhs = 0.0d0

    k = this % current_step

    ! Negate the RHS
    rhs = - rhs

    if(allocated(jac)) deallocate(jac)

  end subroutine assembleRHS


  !===================================================================!
  ! Compute the total derivative of the function with respect to the
  ! design variables and return the gradient 
  !===================================================================!

  subroutine computeTotalDerivative( this, dfdx )

    class(NBG)                                         :: this
    type(scalar) , dimension(:), intent(inout)             :: dfdx
    type(scalar) , dimension(this % nSVars, this % nDVars) :: dRdX
    type(scalar)                                           :: scale = 1.0d0
    integer                                            :: k

    !scale = this % h

    dfdx = 0.0d0

    !-----------------------------------------------------------------!
    ! Compute dfdx
    !-----------------------------------------------------------------!

    do k = 2, this % num_steps
       call this % system % func % addFuncDVSens(dfdx, scale, this % time(k), &
            & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:) )
    end do

    ! Initial condition
!!$    call this % system % func % addFuncDVSens(dfdx, scale, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(1,:) )

    !-----------------------------------------------------------------!
    ! Compute the total derivative
    !-----------------------------------------------------------------!

    do k = 2, this % num_steps
       call this % system % getResidualDVSens(dRdX, scale, this % time(k), &
            & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:))
       dfdx = dfdx + matmul(this % psi(k,:), dRdX) ! check order
    end do

!!$    ! Add constraint contribution
!!$    call this % system % getResidualDVSens(dRdX, scale, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(1,:))
!!$    dfdx = dfdx + matmul(this % psi(2,:), dRdX)

    ! Finally multiply by the scalar
    dfdx = this % h * dfdx

    print*, "Check scaling of dfdx, and transpose"

  end subroutine computeTotalDerivative

end module nbg_integrator
