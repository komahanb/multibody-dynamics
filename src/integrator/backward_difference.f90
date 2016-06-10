!=====================================================================!
! Backward Difference Formula Integration Module for first and second
! order systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module bdf_integrator

  use iso_fortran_env , only : dp => real64
  use integrator_class, only : integrator
  use physics_class,    only : physics

  implicit none

  private
  public :: BDF

  !===================================================================! 
  ! BDF Integrator type
  !===================================================================! 

  type, extends(integrator) :: BDF

     integer                                :: max_bdf_order
     real(dp) , dimension(:), allocatable   :: beta, gamm
     real(dp) , dimension(:,:), allocatable :: psi
     integer :: ndvars

   contains
     
     private 

     ! Routines for integration
     procedure, public :: initialize, finalize, integrate     
     procedure         :: approximateStates
     procedure         :: getOrderCoeff

     ! Routines for adjoint gradient
     procedure         :: computeTotalDerivative
     procedure, public :: marchBackwards
     procedure, public :: getAdjointGradient
     procedure, public :: assembleRHS

  end type BDF

contains

  !===================================================================!
  ! Subroutine that marches backwards in time to compute the lagrange
  ! multipliers (adjoint variables for the function)
  ! ===================================================================!
  
  subroutine marchBackwards(this)

    class(BDF) :: this
    integer    :: k, dummy
    real(dp)   :: alpha, beta, gamma
    
    do k = this % num_steps, 2, -1
       
       this % current_step = k
       
       !--------------------------------------------------------------!
       ! Determine the linearization coefficients for the Jacobian
       !--------------------------------------------------------------!
       
       call this % getOrderCoeff(dummy, this % beta, 1, k)
       call this % getOrderCoeff(dummy, this % gamm, 2, k)
       
       alpha = 1.0d0       
       beta  = this % beta(1)/this % h
       
       if ( this % second_order ) then
          gamma = this % gamm(1)/this % h/this % h
       else 
          gamma = 0.0d0
       end if

       !--------------------------------------------------------------!
       ! Solve the adjoint equation at each step
       !--------------------------------------------------------------!

       call this % adjointSolve(this % psi(k,:), alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))
       
       print *, this % psi(k,:), k

    end do

  end subroutine marchBackwards

  !===================================================================!
  ! Compute the total derivative of the function with respect to the
  ! design variables and return the gradient 
  !===================================================================!

  subroutine computeTotalDerivative(this, dfdx)
    
    class(BDF)                               :: this
    real(dp) , dimension(:), intent(inout)   :: dfdx
    real(dp) , dimension(:,:), allocatable   :: dRdX
    integer                                  :: k
    
    allocate(dRdX(this % nsvars, this % ndvars))
    dfdx = 0.0d0

    ! compute dfdx
    do k = 2, this % num_steps
       call this % system % func % addDfdx(dfdx, 1.0d0, this % time(k), &
            & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:) )
    end do

    ! Compute the total derivative
    do k = 2, this % num_steps
       call this % system % getResidualDVSens(dRdX, 1.0d0, this % time(k), &
            & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:))
       dfdx = dfdx + matmul(this % psi(k,:), dRdX) ! check order
    end do

    deallocate(dRdX)

  end subroutine computeTotalDerivative

  !===================================================================!
  ! Public wrapper for all the adjoint gradient related sequence of
  ! calls
  !===================================================================!
  
  subroutine getAdjointGradient(this, dfdx)

    class(BDF)                             :: this
    real(dp) , dimension(:), intent(inout) :: dfdx

    !-----------------------------------------------------------------!
    ! Integrate backwards to solve for lagrange multipliers for the
    ! set design variable
    !-----------------------------------------------------------------!

    call this % marchBackwards()
    
    !-----------------------------------------------------------------!
    ! Compute the total derivative of the function with respect to the
    ! design variables
    !-----------------------------------------------------------------!

    call this % computeTotalDerivative(dfdx)

  end subroutine getAdjointGradient

  !===================================================================!
  ! Initialize the BDF datatype and allocate required variables
  !===================================================================!
  
  subroutine initialize(this, system, tinit, tfinal, h, second_order, max_bdf_order)
    
    class(BDF)                      :: this
    class(physics), target          :: system
    integer  , OPTIONAL, intent(in) :: max_bdf_order
    real(dp) , OPTIONAL, intent(in) :: tinit, tfinal
    real(dp) , OPTIONAL, intent(in) :: h
    logical  , OPTIONAL, intent(in) :: second_order

    print *, "======================================"
    print *, ">>   Backward Difference Formula    << "
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
    if (.not.(this % nsvars .gt. 0)) stop ">> Error: Zero state variable. Stopping."

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

    if (present(max_bdf_order)) then
       this % max_bdf_order = max_bdf_order
    end if
    print '("  >> Max BDF Order          : ",i4)', this % max_bdf_order

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
    
    !-----------------------------------------------------------------!
    ! Allocate space for BDF coefficients
    !-----------------------------------------------------------------!
    
    allocate(this % beta(this % max_bdf_order + 1))
    this % beta = 0.0d0

    allocate(this % gamm(2*this % max_bdf_order + 1))
    this % gamm = 0.0d0
    
    ! Set the start time
    this % time(1) = this % tinit

  end subroutine initialize

  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!
  
  subroutine finalize(this)

    class(BDF) :: this

    ! Clear global states and time
    if(allocated(this % UDDOT)) deallocate(this % UDDOT)
    if(allocated(this % UDOT)) deallocate(this % UDOT)
    if(allocated(this % U)) deallocate(this % U)
    if(allocated(this % time)) deallocate(this % time)

    if(allocated(this % psi)) deallocate(this % psi)

    if(allocated(this % beta)) deallocate(this % beta)
    if(allocated(this % gamm)) deallocate(this % gamm)

  end subroutine finalize

  !===================================================================!
  ! Time integration logic
  !===================================================================!

  subroutine Integrate(this)

    class(BDF) :: this
    real(dp)   :: alpha, beta, gamma
    integer    :: k

    ! Set the initial condition
    call this % system % getInitialStates(this % time(1), &
         & this % u(1,:), this % udot(1,:))

    this % current_step = 1

    ! March in time
    march: do k = 2, this % num_steps
       
       this % current_step = this % current_step + 1

       ! Increment the time (states are already advanced after the
       ! Newton solve)
       this % time(k) = this % time(k-1) + this % h
       
       ! Approximate the states u, udot and uddot using BDF stencil
       call this % approximateStates()
       
       ! Determine the coefficients for linearing the Residual
       alpha = 1.0d0
       beta  = this % beta(1)/this % h
       if ( this % second_order ) then
          gamma = this % gamm(1)/this % h/this % h
       else 
          gamma = 0.0d0
       end if

       ! Solve the nonlinear system at each step by driving the
       ! residual to zero
       call this % newtonSolve(alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))
       
    end do march

  end subroutine Integrate
  
  !===================================================================!
  ! Approximate the state variables at each step using BDF formulae
  !===================================================================!
  
  subroutine approximateStates(this)

    class(BDF) :: this
    integer :: k, m, i

    k = this % current_step
    
    !-----------------------------------------------------------------!
    ! Extrapolate U to next time step
    !-----------------------------------------------------------------!

    this % u(k,:) = this % u(k-1,:) + this % udot(k-1,:)*this % h &
         & + this % uddot(k-1,:)*this % h*this % h/2.0d0

    !-----------------------------------------------------------------!
    ! Approximate UDOT using BDF
    !-----------------------------------------------------------------!
   
    call this % getOrderCoeff(m, this % beta, 1, k)

    do i = 1, m + 1
       this % udot(k,:) = this % udot(k,:) &
            & + this % beta(i)*this % u(k-i+1,:)/this % h
    end do
    
    !-----------------------------------------------------------------!
    ! Approximate UDDOT using BDF
    !-----------------------------------------------------------------!
    
    call this % getOrderCoeff(m, this % gamm, 2, k)

    if ( m .eq. 0) then

       ! We dont have enought points yet
       this % uddot(k,:) = (this % udot(k,:) - this % udot(k-1,:))/this % h
       
    else
       
       do i = 1, 2*m + 1
          this % uddot(k,:) = this % uddot(k,:) &
               & + this % gamm(i)*this % u(k-i+1,:)/this % h/this % h
       end do
       
    end if
    
  end subroutine approximateStates

  !-------------------------------------------------------------------!
  ! Returns the bdf coeffs (unscaled with respect to the step size h)
  ! and the order
  !-------------------------------------------------------------------!
  ! Input:
  !-------------------------------------------------------------------!
  ! d : d-th derivative
  !-------------------------------------------------------------------!
  ! Output:
  ! m : order of accuracy
  ! coeff: the array of coefficients
  !-------------------------------------------------------------------!
  
  subroutine getOrderCoeff(this, m, coeff, d, k)
    
    class(BDF)              :: this
    integer  , intent(in)   :: d
    integer  , intent(out)  :: m
    real(dp) , intent(out)  :: coeff(:)
    integer  , intent(in)   :: k

    m = (k-1)/d ! k = md + 1
    if ( m .gt. this % max_bdf_order ) m = this % max_bdf_order

    ! set the BDF coefficient for first derivative    
    if (d.eq.1) then

       if (m .eq. 1) then
          coeff(1:m+1) = (/ 1.0, -1.0 /)
       else if (m .eq. 2) then
          coeff(1:m+1) = (/ 1.5_dp, -2.0_dp, 0.5_dp /)
       else  if (m .eq. 3) then
          coeff(1:m+1) = (/ 35.d0/24.0d0, -15.d0/8.0d0, 3.0d0/8.0d0, 1.0d0/24.0d0 /)
       else 
          print *, "wrong order", m
          stop
       end if

    else if (d.eq.2) then

       ! set the BDF coefficient for second derivative   
       if (m.eq.0) then
          coeff(1:2*m+1) = (/ 1.0_dp /) ! used just during linearization
       else if (m.eq.1) then
          ! coeff(1:2*m+1) = (/ 2.25d0, -5.0d0, 2.75d0 /) ! TACS
          coeff(1:2*m+1) = (/ 1.0_dp, -2.0_dp, 1.0_dp /)
       else if (m.eq.2) then
          coeff(1:2*m+1) = (/ 2.25_dp, -6.0_dp, 5.5_dp, -2.0_dp, 0.25_dp /)
       else if (m.eq.3) then
          coeff(1:2*m+1) = (/ 2.126736d0, -5.468750d0, 4.609375d0, &
               & -1.284722d0, -0.015625d0, 0.031250d0, 0.001736d0 /)
       else
          print *, "wrong order", m
          stop
       end if

    else
       print*, "wrong degree", d
       stop
    end if

  end subroutine getOrderCoeff

  !===================================================================!
  ! Function that puts together the right hand side of the adjoint
  ! equation into the supplied rhs vector.
  !===================================================================!
  
  subroutine assembleRHS(this, rhs)

    class(BDF)                            :: this
    real(dp), dimension(:), intent(inout) :: rhs
    integer                               :: k, i, m1, m2
    real(dp)                              :: scale = 0.0d0
    real(dp), allocatable, dimension(:,:) :: jac

    allocate(jac(this%nsvars,this%ndvars))

    k = this % current_step
    
    ! Zero the RHS first
    rhs = 0.0d0
    
    !-----------------------------------------------------------------!
    ! Add function contribution
    !-----------------------------------------------------------------!
    
    scale = 1.0d0
    call this % system % func % addDFdU(rhs, scale, this % time(k), &
         & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:))

    call this % getOrderCoeff(m1, this % beta, 1, k )
    do i = 1, m1 + 1
       if ( k+i-1 .le. this % num_steps) then
          call this % getOrderCoeff(m1, this % beta, 1, k+i-1 )
          scale = this % beta(i)/this % h
          call this % system % func % addDFdUDot(rhs, scale, this % time(k+i-1), &
               & this % system % x, this % u(k+i-1,:), this % udot(k+i-1,:), this % uddot(k+i-1,:))
       end if
    end do

    call this % getOrderCoeff(m2, this % gamm, 2, k)
    do i = 1, 2*m2 + 1
       if ( k+i-1 .le. this % num_steps) then
          call this % getOrderCoeff(m2, this % gamm, 2, k+i-1)
          scale = this % gamm(i)/this % h/this % h
          call this % system % func % addDFdUDDot(rhs, scale, this % time(k+i-1), &
               & this % system % x, this % u(k+i-1,:), this % udot(k+i-1,:), this % uddot(k+i-1,:))
       end if
    end do

    !-----------------------------------------------------------------!
    ! Add residual contribution
    !-----------------------------------------------------------------!

    call this % getOrderCoeff(m1, this % beta, 1, k)    
    do i = 2, m1 + 1
       if ( k+i-1 .le. this % num_steps) then
          call this % getOrderCoeff(m1, this % beta, 1, k+i-1)
          scale = this % beta(i)/this % h
          call this % system % assembleJacobian(jac, 0.0d0, 1.0d0, 0.0d0, &
               & this % time(k+i-1), this % u(k+i-1,:), this % udot(k+i-1,:), this % uddot(k+i-1,:))
          rhs = rhs + scale*matmul( transpose(jac), this % psi(k+i-1,:) )
       end if
    end do

    call this % getOrderCoeff(m2, this % gamm, 2, k)
    do i = 2, 2*m2 + 1
        if ( k+i-1 .le. this % num_steps) then
          call this % getOrderCoeff(m2, this % gamm, 2, k+i-1)
          scale = this % gamm(i)/this % h/this % h
          call this % system % assembleJacobian(jac, 0.0d0, 0.0d0, 1.0d0, &
               & this % time(k+i-1), this % u(k+i-1,:), this % udot(k+i-1,:), this % uddot(k+i-1,:))
          rhs = rhs + scale*matmul( transpose(jac), this % psi(k+i-1,:) )
       end if
    end do

    ! Negate the RHS
    rhs = -rhs
    
    if(allocated(jac)) deallocate(jac)

  end subroutine assembleRHS

end module bdf_integrator
