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
     
     ! Average Constant Accelearation (second order unconditionally stable)
     type(scalar) :: BETA   = 0.25d0
     type(scalar) :: GAMMA  = 0.50d0

!!$     ! Fox & Goodwin  (third order & conditionally stable wh=2.45)
!!$     type(scalar) :: BETA   = 1.0d0/12.0d0
!!$     type(scalar) :: GAMMA  = 0.50d0
!!$
!!$     ! Linear Acceleration (second order & conditionally stable wh=3.46)
!!$     type(scalar) :: BETA   = 1.0d0/6.0d0
!!$     type(scalar) :: GAMMA  = 0.50d0
!!$
!!$     ! Central Difference (second order & conditionally stable wh=2)
!!$     type(scalar) :: BETA   = 1.0d0/2.0d0
!!$     type(scalar) :: GAMMA  = 0.50d0
!!$
!!$     ! Purely Explicit
!!$     type(scalar) :: BETA   = 0.0d0
!!$     type(scalar) :: GAMMA  = 0.0d0
     
     type(scalar), allocatable, dimension(:) :: rho
     type(scalar), allocatable, dimension(:) :: sigma
     
   contains

     ! Destructor
     procedure, public  :: finalize
     
     ! Routines for integration
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff

     ! Routines for adjoint gradient
     procedure, public  :: marchBackwards
     procedure, private :: assembleRHS
     procedure          :: computeTotalDerivative

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
    print *, ">>   Newmark Beta Gamma (NBG)      << "
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

    allocate(this % rho(this % nsvars))
    this % rho = 0.0d0
    
    allocate(this % sigma(this % nsvars))
    this % sigma = 0.0d0

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
    if(allocated(this % rho)) deallocate(this % rho)
    if(allocated(this % sigma)) deallocate(this % sigma)
   
  end subroutine finalize
  
  !===================================================================!
  ! Approximate the state variables at each step using NBG formula
  !===================================================================!
  
  subroutine approximateStates( this )

    class(NBG) :: this
    integer                   :: k
    type(scalar)              :: scale

    k = this % current_step

    !-----------------------------------------------------------------!
    ! Assume a UDDOT for the next time step
    !-----------------------------------------------------------------!

    this % uddot(k,:) = this % uddot(k-1,:)

    !-----------------------------------------------------------------!
    ! Approximate UDOT using NBG
    !-----------------------------------------------------------------!

    this % udot(k,:) = this % udot(k-1,:) 

    scale = this % h * (1.0d0 - this % GAMMA)
    this % udot(k,:) = this % udot(k,:) + scale*this % uddot(k-1,:) 

    scale = this % h * this % GAMMA
    this % udot(k,:) = this % udot(k,:) + scale*this % uddot(k,:) 

    !-----------------------------------------------------------------!
    ! Approximate U using NBG
    !-----------------------------------------------------------------!

    this % u(k,:) = this % u(k-1,:) 

    scale = this % h
    this % u(k,:) = this % u(k,:) + scale*this % udot(k-1,:) 

    scale = this % h * this % h * (1.0d0 - 2.0d0 * this % BETA)/2.0d0
    this % u(k,:) = this % u(k,:) + scale*this % uddot(k-1,:) 

    scale = this % h * this % h * this % BETA
    this % u(k,:) = this % u(k,:) + scale*this % uddot(k,:) 

  end subroutine approximateStates
  
  !===================================================================!
  ! Returns the linearization scalar coefficients: alpha, beta, gamma
  !===================================================================!
  
  subroutine getLinearCoeff( this, alpha, beta, gamma )

    class(NBG) :: this
    type(scalar), intent(out) :: alpha, beta, gamma
    integer                     :: m
    
    if ( this % second_order ) then
       gamma = 1.0d0/this % h/ this % h
       beta  = this % GAMMA/ this % h
       alpha = this % BETA
    else
       stop"Error: Newmark-Beta-Gamma method works for second order systems in current form..."
    end if

!!$    gamma = 1.0d0
!!$    beta  = this % h * this % GAMMA
!!$    alpha = this % h * this % h * this % BETA

  end subroutine getLinearCoeff

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
       call this % getLinearCoeff(alpha, beta, gamma)

       !--------------------------------------------------------------!
       ! Solve the adjoint equation at each step
       !--------------------------------------------------------------!

       call this % adjointSolve(this % psi(k,:), alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))

       print*, "k, psi=", k, this % psi(k,:)       

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
    type(scalar)                              :: alpha, beta, gamma
    type(integer)                             :: k

    allocate( jac(this % nSVars, this % nSVars) )
    
    k = this % current_step

    ! Zero the RHS first
    rhs = 0.0d0

    if ( k .ne. this % num_steps) then

       !-----------------------------------------------------!
       ! Add previous residual contributions to RHS
       !-----------------------------------------------------!

       gamma = 0.0d0
       beta  = 1.0d0/this % h
       alpha = 0.5d0 + this % GAMMA

       ! Add the state variable sensitivity from the previous step
       call this % system % func % addFuncSVSens(rhs, &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % system % X, &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       ! Add the residual adjoint product from the previous step
       call this % system % assembleJacobian(jac, &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       rhs = rhs + matmul(transpose(jac(:,:)), this % psi(k+1,:))

       !-------------------------------------------------------!
       ! Add similar contributions to RHS
       !-------------------------------------------------------!

       rhs = rhs + beta  * this % sigma/this % h
       rhs = rhs + alpha * this % rho/this % h
       
       !-----------------------------------------------------------!
       ! Compute NEW THIS % SIGMA (k)
       !-----------------------------------------------------------!

       this % sigma = this % sigma + this % h * this % rho

       gamma = 0.0d0
       beta  = this % h
       alpha = this % h * this % h

       ! Add the state variable sensitivity from the previous step
       call this % system % func % addFuncSVSens(this % sigma, &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % system % X, &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       ! Add the residual adjoint product from the previous step
       call this % system % assembleJacobian(jac, &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       this % sigma = this % sigma + matmul(transpose(jac(:,:)), this % psi(k+1,:))

       !-----------------------------------------------------------!
       ! Compute THIS % RHO (K)
       !-----------------------------------------------------------!

       gamma = 0.0d0
       beta  = 0.0d0
       alpha = this % h

       ! Add the state variable sensitivity from the previous step
       call this % system % func % addFuncSVSens(this % rho, &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % system % X, &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       ! Add the residual adjoint product from the previous step
       call this % system % assembleJacobian(jac, &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       this % rho = this % rho + matmul(transpose(jac(:,:)), this % psi(k+1,:))
      
    end if

    ! Get the coefficients
    call this % getLinearCoeff( alpha, beta, gamma)
    
    ! Add the state variable sensitivity
    call this % system % func % addFuncSVSens(rhs, &
         & alpha, beta, gamma, &
         & this % time(k), &
         & this % system % X, &
         & this % u(k,:), &
         & this % udot(k,:), &
         & this % uddot(k,:))

    ! Negate the RHS
    rhs = -rhs

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
