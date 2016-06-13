!=====================================================================!
! Backward Difference Formula Integration Module for first and second
! order systems with adjoint derivative capabilities.
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
  ! A derived type for the bdf coefficients
  !===================================================================!
  
  type :: bdf_coeff
     
     ! Information
     integer                               :: max_order = 3
     
     ! Coeff values
     real(dp)                              :: alpha
     real(dp), dimension(:,:), allocatable :: beta
     real(dp), dimension(:,:), allocatable :: gamma

   contains

     private
     
     procedure :: destruct
     procedure :: getOrder
     
  end type bdf_coeff
  
  ! Interface for the constructor of bdf_coeff type
  interface bdf_coeff
     module procedure construct_bdf_coeff
  end interface bdf_coeff
  
  !===================================================================! 
  ! BDF Integrator type
  !===================================================================! 

  type, extends(integrator) :: BDF

! BDF variables
     integer                                :: max_bdf_order = 3
     type(bdf_coeff)                        :: coeff

     real(dp) , dimension(:,:), allocatable :: psi

   contains

 ! Routines for integration
     procedure, public  :: finalize, integrate     
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff

 ! Routines for adjoint gradient

     procedure, public  :: marchBackwards
     procedure, private :: assembleRHS
     procedure, private :: computeTotalDerivative
     procedure, public  :: evalFunc

! Overridden procedure
     procedure :: writeSolution => writeSolutionAdjoint

  end type BDF
  
  interface BDF
     module procedure initialize
  end interface BDF

contains
  
  !===================================================================!
  ! Returns the linearization scalar coefficients: alpha, beta, gamma
  !===================================================================!
  
  subroutine getLinearCoeff( this, k, alpha, beta, gamma )
    
    class(BDF), intent(inout) :: this
    integer, intent(in)       :: k
    real(dp), intent(inout)   :: alpha, beta, gamma
    
    alpha = this % coeff % alpha

    beta  = this % coeff % beta(this % coeff % getOrder(k,1), 1)/this % h
    
    if ( this % second_order ) then
       gamma = this % coeff % gamma(this % coeff % getOrder(k, 2), 1)/this % h/this % h
    else 
       gamma = 0.0d0
    end if
    
  end subroutine getLinearCoeff

  !===================================================================!
  ! Constructor for BDF coeff object
  !===================================================================!
  
  type(bdf_coeff) function construct_bdf_coeff( max_order ) result( this )

    integer, intent(in), OPTIONAL :: max_order
    
    if( present(max_order) ) then
       this % max_order = max_order
    end if
    
    allocate( this % beta (0:this % max_order, this % max_order + 1) )
    this % beta = 0.0d0 

    allocate( this % gamma (0:this % max_order, 2*this % max_order + 1) ) ! 2m + 1
    this % gamma = 0.0d0

    ! Set the coefficients

    this % alpha = 1.0d0
    this % gamma(0, 1:1) = 1.0d0

    if (this % max_order .eq. 3) then

       this % beta (1, 1:2) = (/ 1.0, -1.0 /)
       this % beta (2, 1:3) = (/ 1.5_dp, -2.0_dp, 0.5_dp /)
       this % beta (3, 1:4) = (/ 35.d0/24.0d0, -15.d0/8.0d0, 3.0d0/8.0d0, 1.0d0/24.0d0 /)

       this % gamma(1, 1:3) = (/ 1.0_dp, -2.0_dp, 1.0_dp /)
       this % gamma(2, 1:5) = (/ 2.25_dp, -6.0_dp, 5.5_dp, -2.0_dp, 0.25_dp /)
       this % gamma(3, 1:7) =  (/ 2.126736d0, -5.468750d0, 4.609375d0, &
            & -1.284722d0, -0.015625d0, 0.031250d0, 0.001736d0 /)

    else if (this % max_order .eq. 2) then

       this % beta (1, 1:2) = (/ 1.0, -1.0 /)
       this % beta (2, 1:3) = (/ 1.5_dp, -2.0_dp, 0.5_dp /)

       this % gamma(1, 1:3) = (/ 1.0_dp, -2.0_dp, 1.0_dp /)
       this % gamma(2, 1:5) = (/ 2.25_dp, -6.0_dp, 5.5_dp, -2.0_dp, 0.25_dp /)

    else if (this % max_order .eq. 1) then

       this % beta (1, 1:2) = (/ 1.0, -1.0 /)
       this % gamma(1, 1:3) = (/ 1.0_dp, -2.0_dp, 1.0_dp /)
    else 
       print *,  "Wrong max_bdf_order:", this % max_order
       stop
    end if

  end function construct_bdf_coeff
  
  !===================================================================!
  ! Destructor for BDF coeff object
  !===================================================================!
  
  subroutine destruct( this ) 
    
    class(bdf_coeff) :: this
    
    if ( allocated(this % beta) ) deallocate(this % beta)
    if ( allocated(this % gamma) ) deallocate(this % gamma)
    
  end subroutine destruct
  
  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!
  
  pure integer function getOrder(this, k, d )

    class(bdf_coeff), intent(in)    :: this
    integer, intent(in)             :: k, d

    ! find the order of approximation
    getOrder = (k-1)/d ! k = md + 1

    ! Do not let exceed the max order sought
    if ( getOrder .gt. this % max_order ) getOrder = this % max_order

  end function getOrder
  
  !===================================================================!
  ! Write solution to file
  !===================================================================!
  
  subroutine writeSolutionAdjoint(this, filename)

    class(BDF)                             :: this
    character(len=*), OPTIONAL, intent(in) :: filename
    character(len=7), parameter            :: directory = "output/"
    character(len=32)                      :: path = ""
    integer                                :: k, j, ierr

    path = trim(path)

    if (present(filename)) then
       path = directory//filename
    else
       path = directory//"solution.dat"
    end if

    open(unit=90, file=trim(path), iostat= ierr)

    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') path
       return
    end if

    do k = 1, this % num_steps
       write(90, *)  this % time(k), (this % u(k,j), j=1,this % nsvars ), &
            & (this % udot(k,j), j=1,this % nsvars ), &
            & (this % uddot(k,j), j=1,this % nsvars ), &
            & (this % psi(k,j), j=1,this % nsvars)
    end do

    close(90)

  end subroutine writeSolutionAdjoint

  !===================================================================!
  ! Compute the total derivative of the function with respect to the
  ! design variables and return the gradient 
  !===================================================================!
  
  subroutine computeTotalDerivative( this, dfdx )
    
    class(BDF)                                         :: this
    real(dp) , dimension(:), intent(inout)             :: dfdx
    real(dp) , dimension(this % nSVars, this % nDVars) :: dRdX
    real(dp)                                           :: scale
    integer                                            :: k
    
    scale = this % h
    
    dfdx = 0.0d0
    
    !-----------------------------------------------------------------!
    ! Compute dfdx
    !-----------------------------------------------------------------!

    do k = 2, this % num_steps
       call this % system % func % addDfdx(dfdx, scale, this % time(k), &
            & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:) )
    end do

    ! Initial condition
    call this % system % func % addDfdx(dfdx, scale, this % time(1), &
         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(2,:) )
    
    !-----------------------------------------------------------------!
    ! Compute the total derivative
    !-----------------------------------------------------------------!

    do k = 2, this % num_steps
       call this % system % getResidualDVSens(dRdX, scale, this % time(k), &
            & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:))
       dfdx = dfdx + matmul(this % psi(k,:), dRdX) ! check order
    end do

    ! Add constraint contribution
    call this % system % getResidualDVSens(dRdX, scale, this % time(1), &
         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(2,:))
    dfdx = dfdx + matmul(this % psi(2,:), dRdX)

    ! Finally multiply by the scalar
    !   dfdx = this %  * dfdx
    print*, "Check scaling of dfdx"

  end subroutine computeTotalDerivative

  !===================================================================!
  ! Initialize the BDF datatype and allocate required variables
  !===================================================================!
  
  type(bdf) function initialize( system, tinit, tfinal, h, second_order, max_bdf_order )  result (this)
   
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
    
    if ( .not. (this % nsvars .gt. 0) ) then
       stop ">> Error: Zero state variable. Stopping."
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
    ! Create the BDF-coeff object
    !-----------------------------------------------------------------!
    
    this % coeff = bdf_coeff(this % max_bdf_order)
    
    ! Set the start time
    this % time(1) = this % tinit

  end function initialize

  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!
  
  subroutine finalize( this )

    class(BDF) :: this

    ! Clear global states and time
    if(allocated(this % UDDOT)) deallocate(this % UDDOT)
    if(allocated(this % UDOT)) deallocate(this % UDOT)
    if(allocated(this % U)) deallocate(this % U)
    if(allocated(this % time)) deallocate(this % time)

    if(allocated(this % psi)) deallocate(this % psi)

    ! call the BDF coeff destructor
    call this % coeff % destruct()

  end subroutine finalize

  !===================================================================!
  ! Time integration logic
  !===================================================================!

  subroutine Integrate( this )

    class(BDF) :: this
    real(dp)   :: alpha, beta, gamma
    integer    :: k

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
       
       ! Approximate the states u, udot and uddot using BDF stencil
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

    class(BDF)                :: this
    integer                   :: k
    real(dp)                  :: alpha, beta, gamma
    
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
       
       ! print *, this % psi(k,:), k

    end do time

  end subroutine marchBackwards
  
  !===================================================================!
  ! Approximate the state variables at each step using BDF formulae
  !===================================================================!
  
  subroutine approximateStates( this )

    class(BDF), intent(inout) :: this
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
    
    m = this % coeff % getOrder(k, 1)

    do i = 1, m + 1
       this % udot(k,:) = this % udot(k,:) &
            & + this % coeff % beta(m, i)*this % u(k-i+1,:)/this % h
    end do
    
    !-----------------------------------------------------------------!
    ! Approximate UDDOT using BDF
    !-----------------------------------------------------------------!
    
    m = this % coeff % getOrder(k, 2)

    if ( m .eq. 0) then

       ! We dont have enought points yet
       this % uddot(k,:) = (this % udot(k,:) - this % udot(k-1,:))/this % h
       
    else
       
       do i = 1, 2*m + 1
          this % uddot(k,:) = this % uddot(k,:) &
               & + this % coeff % gamma(m, i)*this % u(k-i+1,:)/this % h/this % h
       end do
       
    end if
    
  end subroutine approximateStates

  !===================================================================!
  ! Function that puts together the right hand side of the adjoint
  ! equation into the supplied rhs vector.
  !===================================================================!
  
  subroutine assembleRHS( this, rhs )

    class(BDF)                            :: this
    real(dp), dimension(:), intent(inout) :: rhs
    real(dp), dimension(:,:), allocatable :: jac
    real(dp)                              :: scale
    integer                               :: k, i, m1, m2
    
    allocate( jac(this % nSVars, this % nSVars) )
    
    k = this % current_step
    
    ! Zero the RHS first
    rhs = 0.0d0
    
    !-----------------------------------------------------------------!
    ! Add function contribution (dfdu)
    !-----------------------------------------------------------------!
    
    call this % system % func % addDFdU(rhs, 1.0d0, this % time(k), &
         & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:))

    m1 = this % coeff % getOrder(k, 1)
    do i = 1, m1 + 1
       if ( k+i-1 .le. this % num_steps) then
          scale = this % coeff % beta(m1, 1)/this % h
          call this % system % func % addDFdUDot(rhs, scale, this % time(k+i-1), &
               & this % system % x, this % u(k+i-1,:), this % udot(k+i-1,:), this % uddot(k+i-1,:))
       end if
    end do

    m2 = this % coeff % getOrder(k, 2)
    do i = 1, 2*m2 + 1
       if ( k+i-1 .le. this % num_steps) then
          scale = this % coeff % gamma(m2, 1)/this % h/this % h
          call this % system % func % addDFdUDDot(rhs, scale, this % time(k+i-1), &
               & this % system % x, this % u(k+i-1,:), this % udot(k+i-1,:), this % uddot(k+i-1,:))
       end if
    end do

    !-----------------------------------------------------------------!
    ! Add residual contribution
    !-----------------------------------------------------------------!

    do i = 2, m1 + 1
       if ( k+i-1 .le. this % num_steps) then
          scale = this % coeff % beta(m1, 1)/this % h
          call this % system % assembleJacobian(jac, 0.0d0, 1.0d0, 0.0d0, &
               & this % time(k+i-1), this % u(k+i-1,:), this % udot(k+i-1,:), this % uddot(k+i-1,:))
          rhs = rhs + scale*matmul( transpose(jac), this % psi(k+i-1,:) )
       end if
    end do

    do i = 2, 2*m2 + 1
       if ( k+i-1 .le. this % num_steps) then
          scale = this % coeff % gamma(m2, 1)/this % h/this % h
          call this % system % assembleJacobian(jac, 0.0d0, 0.0d0, 1.0d0, &
               & this % time(k+i-1), this % u(k+i-1,:), this % udot(k+i-1,:), this % uddot(k+i-1,:))
          rhs = rhs + scale*matmul( transpose(jac), this % psi(k+i-1,:) )
       end if
    end do
    
    ! Negate the RHS
    rhs = - rhs
    
    if(allocated(jac)) deallocate(jac)
    
  end subroutine assembleRHS

  !===================================================================!
  ! Evaluating the function of interest
  !===================================================================!

  subroutine evalFunc(this, x, fval)

    class(BDF)                            :: this
    real(dp), dimension(:), intent(in)    :: x
    real(dp), intent(inout)               :: fval
    integer                               :: k
    real(dp), dimension(this % num_steps) :: ftmp
    
!    print*, "Evaluating function of interest"
    
    do concurrent(k = 1 : this % num_steps)
       call this % system % func % getFunctionValue(ftmp(k), this % time(k), &
            & x, this % U(k,:), this % UDOT(k,:), this % UDDOT(k,:))
       ftmp(k) = this % h * ftmp(k)
    end do
    
    ! fval = sum(ftmp)/dble(this % num_steps)
    fval = sum(ftmp)
   
  end subroutine evalFunc

end module bdf_integrator
