#include "scalar.fpp"
!=====================================================================!
! Adams Bashworth Moulton Integration Module for first and second
! order systems with adjoint derivative capabilities.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module abm_integrator

  use integrator_class, only : integrator
  use physics_class,    only : physics

  implicit none

  private
  public :: ABM
  
  type(scalar), parameter :: ONE  = 1.0d0
  type(scalar), parameter :: ZERO = 0.0d0

  !===================================================================! 
  ! ABM Integrator type
  !===================================================================! 

  type, extends(integrator) :: ABM
     
     ! ABM variables
     integer :: max_abm_order = 3
     type(scalar), allocatable, dimension(:,:) :: A

   contains

     ! Routines for integration
     procedure, public  :: finalize, integrate     
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff
     procedure, private :: getOrder

     ! Routines for adjoint gradient
     procedure, public  :: marchBackwards
     procedure, private :: assembleRHS
     procedure :: computeTotalDerivative

  end type ABM

  interface ABM
     module procedure initialize
  end interface ABM

contains


  !===================================================================!
  ! Initialize the ABM datatype and allocate required variables
  !===================================================================!
  
  type(abm) function initialize( system, tinit, tfinal, h, second_order, max_abm_order )  result (this)
   
    class(physics), target          :: system
    integer  , OPTIONAL, intent(in) :: max_abm_order
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
    ! Set the order of integration
    !-----------------------------------------------------------------!

    if (present(max_abm_order)) then
       this % max_abm_order = max_abm_order
    end if
    print '("  >> Max ABM Order          : ",i4)', this % max_abm_order

    allocate( this % A (this % max_abm_order, this % max_abm_order) )
    this % A = 0.0d0 
 
    ! Set the coefficients
    if ( this % max_abm_order .eq. 1 ) then       
       this % A(1,1:1) = (/ 1.0d0 /)
    else if ( this % max_abm_order .eq. 2 ) then
       this % A(1,1:1) = (/ 1.0d0 /)
       this % A(2,1:2) = (/ 1.0d0/2.0d0, 1.0d0/2.0d0 /)
    else if ( this % max_abm_order .eq. 3 ) then
       this % A(1,1:1) = (/ 1.0d0 /)
       this % A(2,1:2) = (/ 1.0d0/2.0d0, 1.0d0/2.0d0 /)
       this % A(3,1:3) = (/ 5.0d0/12.0d0, 8.0d0/12.0d0, -1.0d0/12.0d0 /)
    else 
       print *,  "Wrong max_abm_order:", this % max_abm_order
       stop
    end if

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

    class(ABM) :: this

    ! Clear global states and time
    if(allocated(this % UDDOT)) deallocate(this % UDDOT)
    if(allocated(this % UDOT)) deallocate(this % UDOT)
    if(allocated(this % U)) deallocate(this % U)
    if(allocated(this % time)) deallocate(this % time)

    if(allocated(this % psi)) deallocate(this % psi)

    ! call the ABM coeff destructor
    if(allocated(this % A)) deallocate(this % A)

  end subroutine finalize

  !===================================================================!
  ! Time integration logic
  !===================================================================!

  subroutine Integrate( this )

    class(ABM)   :: this
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
       
       ! Approximate the states u, udot and uddot using ABM stencil
       call this % approximateStates()

       ! Determine the coefficients for linearing the Residual
       call this % getLinearCoeff(k, alpha, beta, gamma)

       ! Solve the nonlinear system at each step by driving the
       ! residual to zero

!       print *, k, alpha, beta, gamma
       call this % newtonSolve(alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))

    end do time

  end subroutine Integrate
  
  !===================================================================!
  ! Returns the linearization scalar coefficients: alpha, beta, gamma
  !===================================================================!
  
  subroutine getLinearCoeff( this, k, alpha, beta, gamma )
   
    class(ABM), intent(inout)   :: this
    integer, intent(in)         :: k
    type(scalar), intent(inout) :: alpha, beta, gamma
    integer :: m

    m = this % getOrder(k)

    if ( this % second_order ) then
       gamma = 1.0d0
       beta  = this  % A(m, 1)*this % h
       alpha = (this % A(m, 1)*this % h)**2
    else 
       gamma = 0.0d0
       beta  = 1.0d0
       alpha = this  % A(m, 1)*this % h
    end if
    
  end subroutine getLinearCoeff
  
  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!
  
  pure integer function getOrder(this, k)

    class(ABM), intent(in)    :: this
    integer, intent(in)                  :: k

    ! Find the order of approximation
    getOrder = k -1

    ! Do not let exceed the max order sought
    if ( getOrder .gt. this % max_abm_order ) getOrder = this % max_abm_order

  end function getOrder
 
  !===================================================================!
  ! Compute the total derivative of the function with respect to the
  ! design variables and return the gradient 
  !===================================================================!
  
  subroutine computeTotalDerivative( this, dfdx )
    
    class(ABM)                                         :: this
    type(scalar) , dimension(:), intent(inout)             :: dfdx
    type(scalar) , dimension(this % nSVars, this % nDVars) :: dRdX
    type(scalar)                                           :: scale = 1.0d0
    integer                                            :: k
    
!    scale = this % h
    
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

  !===================================================================!
  ! Subroutine that marches backwards in time to compute the lagrange
  ! multipliers (adjoint variables for the function)
  ! ===================================================================!
  
  subroutine marchBackwards( this )

    class(ABM)                :: this
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
  ! Approximate the state variables at each step using ABM formulae
  !===================================================================!
  
  subroutine approximateStates( this )

    class(ABM), intent(inout) :: this
    integer                   :: k, m, i

    k = this % current_step
    
    !-----------------------------------------------------------------!
    ! Assume a UDDOT for the next time step
    !-----------------------------------------------------------------!
    
    if ( k .eq. 2 ) then
       this % uddot(k,:) = 1.0d0
    else
       this % uddot(k,:) = this % uddot(k-1,:) 
    end if
    
    m = this % getOrder(k)

    !-----------------------------------------------------------------!
    ! Approximate UDOT using ABM
    !-----------------------------------------------------------------!
      
    do i = 1, m
       this % udot(k,:) = this % udot(k,:) &
            & + this % h*this % A(m,i) * this % uddot(k-i+1,:)
    end do
    this % udot(k,:) = this % udot(k,:) + this % udot(k-1,:) 

    !-----------------------------------------------------------------!
    ! Approximate U using ABM
    !-----------------------------------------------------------------!
   
    do i = 1, m
       this % u(k,:) = this % u(k,:) &
            & + (this % h*this % A(m, i)**2)*this % udot(k-i+1,:)
    end do
    this % u(k,:) = this % u(k,:) + this % u(k-1,:) 
    
  end subroutine approximateStates

  !===================================================================!
  ! Function that puts together the right hand side of the adjoint
  ! equation into the supplied rhs vector.
  !===================================================================!
  
  subroutine assembleRHS( this, rhs )

    class(ABM)                                :: this
    type(scalar), dimension(:), intent(inout) :: rhs
    type(scalar), dimension(:,:), allocatable :: jac
    type(scalar)                              :: scale
    integer                                   :: k, i, m1, m2, idx
    
    allocate( jac(this % nSVars, this % nSVars) )
    
    k = this % current_step 
    m1 = this % getOrder(k)
    m2 = this % getOrder(k)

    if (m2 .eq. 0) m2 = 1

    ! Zero the RHS first
    rhs = 0.0d0
    
    !-----------------------------------------------------------------!
    ! Add function contribution (dfdu)
    !-----------------------------------------------------------------!
    
    call this % system % func % addDFdU(rhs, ONE, this % time(k), &
         & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:))

    do i = 0, m1
       idx = k + i
       if ( idx .le. this % num_steps) then
!          scale = this % coeff % beta(m1, i+1)/this % h
          call this % system % func % addDFdUDot(rhs, scale, this % time(idx), &
               & this % system % x, this % u(idx,:), this % udot(idx,:), this % uddot(idx,:))
       end if
    end do

    do i = 0, 2*m2
       idx = k + i
       if ( idx .le. this % num_steps) then
 !         scale = this % coeff % gamma(m2, i+1)/this % h/this % h
          call this % system % func % addDFdUDDot(rhs, scale, this % time(idx), &
               & this % system % x, this % u(idx,:), this % udot(idx,:), this % uddot(idx,:))
       end if
    end do

    !-----------------------------------------------------------------!
    ! Add residual contribution
    !-----------------------------------------------------------------!
    
    do i = 1, m1 
       idx = k + i
       if ( idx .le. this % num_steps) then
  !        scale = this % coeff % beta(m1, i+1)/this % h
          call this % system % assembleJacobian(jac, ZERO, ONE, ZERO, &
               & this % time(idx), this % u(idx,:), this % udot(idx,:), this % uddot(idx,:))
          rhs = rhs + scale*matmul( transpose(jac), this % psi(idx,:) )
       end if
    end do

    do i = 1, 2*m2
       idx = k + i
       if ( idx .le. this % num_steps) then
   !       scale = this % coeff % gamma(m2, i+1)/this % h/this % h
          call this % system % assembleJacobian(jac, ZERO, ZERO, ONE, &
               & this % time(idx), this % u(idx,:), this % udot(idx,:), this % uddot(idx,:))
          rhs = rhs + scale*matmul( transpose(jac), this % psi(idx,:) )
       end if
    end do
    
    ! Negate the RHS
    rhs = - rhs
    
    if(allocated(jac)) deallocate(jac)
    
  end subroutine assembleRHS

end module abm_integrator
