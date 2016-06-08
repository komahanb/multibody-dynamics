!=====================================================================!
! Module to handle precision of variables
!=====================================================================!

module precision
  
  use, intrinsic :: iso_fortran_env

  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  integer, parameter :: qp = REAL128
  
end module precision

!=====================================================================!
! Parent class for integration schemes
!=====================================================================!

module integrator_class
  
  use precision, only : dp
  use physics_class

  implicit none

  private
  public ::  integrator

  !===================================================================! 
  ! Abstract Integrator type
  !===================================================================! 

  type, abstract :: integrator
     
     !----------------------------------------------------------------!
     ! Contains the actual physical system
     !----------------------------------------------------------------!
     
     class(physics), pointer :: system => null()
     integer                 :: nsvars = 1 ! number of states/equations

     !----------------------------------------------------------------!
     ! Variables for managing time marching
     !----------------------------------------------------------------!

     integer  :: num_steps                     ! number of time steps
     real(dp) :: tinit = 0.0d0, tfinal = 1.0d0 ! initial and final times
     real(dp) :: h = 0.1d0                     ! default step size

     !----------------------------------------------------------------!
     ! Nonlinear solution at each stage
     !----------------------------------------------------------------!

     integer  :: max_newton = 25
     real(dp) :: atol = 1.0d-12, rtol = 1.0d-8

     !----------------------------------------------------------------!
     ! Track global time and states
     !----------------------------------------------------------------!

     real(dp), dimension(:), allocatable   :: time
     real(dp), dimension(:,:), allocatable :: U
     real(dp), dimension(:,:), allocatable :: UDOT
     real(dp), dimension(:,:), allocatable :: UDDOT

     !----------------------------------------------------------------!
     ! Miscellaneous variables
     !----------------------------------------------------------------!

     logical :: second_order = .false.
     logical :: forward = .true.
     integer :: print_level = 0
     integer :: current_step
     logical :: approximate_jacobian = .false.

   contains
     
     !----------------------------------------------------------------!
     ! Procedures                                                     !
     !----------------------------------------------------------------!

     procedure :: writeSolution
     procedure :: setPhysicalSystem
     procedure :: newtonSolve   
     procedure :: setPrintLevel
     procedure :: approximateJacobian

     !----------------------------------------------------------------!
     ! Important setters
     !----------------------------------------------------------------!
     
     procedure :: setApproximateJacobian

  end type integrator

contains
  
  !===================================================================!
  ! Setter that can be used to set the method in which jacobian needs
  ! to be computed. Setting this to .true. would make the code use
  ! finite differences, this is enabled by default too. If set to
  ! .false. the expects to provide implementation in assembleJacobian
  ! in a type that extends PHYSICS.
  !===================================================================!
  
  subroutine setApproximateJacobian(this, approximateJacobian)

    class(integrator) :: this
    logical :: approximateJacobian

    this % approximate_jacobian = approximateJacobian

  end subroutine setApproximateJacobian

  !===================================================================!
  ! Set ANY physical system that extends the type PHYSICS and provides
  ! implementation to the mandatory functions assembleResidual and
  ! getInitialStates
  !===================================================================!
  
  subroutine setPhysicalSystem(this, physical_system)

    class(integrator)      :: this
    class(physics), target :: physical_system

    this % system => physical_system

  end subroutine setPhysicalSystem
  
  !===================================================================!
  ! Manages the amount of print
  !===================================================================!
  
  subroutine setPrintLevel(this,print_level)

    class(integrator) :: this
    integer           :: print_level

    this % print_level = print_level

  end subroutine setPrintLevel

  !===================================================================!
  ! Write solution to file
  !===================================================================!

  subroutine writeSolution(this, filename)
    
    class(integrator)                       :: this
    character(len=32), OPTIONAL, intent(in) :: filename
    character(len=7), parameter             :: directory = "output/"
    character(len=39)                       :: path=""
    integer                                 :: k, j, ierr
    
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
       write(90, *)  this % time(k), (this % u(k,j), j=1,this%nsvars ), &
            & (this % udot(k,j), j=1,this%nsvars ), &
            & (this % uddot(k,j), j=1,this%nsvars )
    end do

    close(90)

  end subroutine writeSolution
  
  !===================================================================!
  ! Solve the nonlinear system at each step by driving the
  ! residual to zero
  !
  ! Input: 
  ! The guessed (initial) state variable values q, qdot, qddot are
  ! supplied
  !
  ! Output: q, qdot, qddot updated iteratively until the corresponding
  ! residual R = 0
  !
  ! alpha: multiplier for derivative of Residual wrt to q
  ! beta : multiplier for derivative of Residual wrt to qdot
  ! gamma: multiplier for derivative of Residual wrt to qddot
  !===================================================================!

  subroutine newtonSolve( this, alpha, beta, gamma, t, q, qdot, qddot )
    
    class(integrator)                     :: this

 ! Arguments
    real(dp), intent(inout)                  :: alpha, beta, gamma
    real(dp), intent(in)                  :: t
    real(dp), intent(inout), dimension(:) :: q, qdot, qddot
    
 ! Lapack variables
    integer, allocatable, dimension(:)    :: ipiv
    integer                               ::  info, size
   
 ! Norms to tracking progress
    real(dp)                              :: abs_res_norm
    real(dp)                              :: rel_res_norm
    real(dp)                              :: init_norm
    
 ! Other Local variables
    real(dp), allocatable, dimension(:)   :: res, dq
    real(dp), allocatable, dimension(:,:) :: jac, fd_jac

    integer                               :: n, k
    logical                               :: conv = .false.

    real(dp)                              :: jac_err

    ! find the size of the linear system based on the calling object
    size = this % nsvars; k = this % current_step

    if ( .not. allocated(ipiv)   ) allocate( ipiv(size)        )
    if ( .not. allocated(res)    ) allocate( res(size)         )
    if ( .not. allocated(dq)     ) allocate( dq(size)          )
    if ( .not. allocated(jac)    ) allocate( jac(size,size)    )
    if ( .not. allocated(fd_jac) ) allocate( fd_jac(size,size) )

    if ( this % print_level .ge. 1 .and. k .eq. 2) then
       write(*,'(/2A5, 2A12/)') "Step", "Iter", "|R|", "|R|/|R1|"
    end if
    
    newton: do n = 1, this % max_newton

       ! Get the residual of the function
       call this % system % assembleResidual(res, t, q, qdot, qddot)

       ! Get the jacobian matrix
       if ( this % approximate_jacobian ) then
                    
          ! Compute an approximate Jacobian using finite differences
          call this % approximateJacobian(jac, alpha, beta, gamma, t, q, qdot, qddot)

       else
          
          ! Use the user supplied Jacobian implementation
          call this % system % assembleJacobian(jac, alpha, beta, gamma, t, q, qdot, qddot)
          
          ! Check the Jacobian implementation once at the beginning of integration
          if ( k .eq. 2 .and. n .eq. 1 ) then
             
             ! Compute an approximate Jacobian using finite differences
             call this % approximateJacobian(fd_jac, alpha, beta, gamma, t, q, qdot, qddot)

             ! Compare the exact and approximate Jacobians and
             ! complain about the error in Jacobian if there is any
             jac_err = maxval(abs(fd_jac - jac))
             if ( jac_err .gt. 1.0d-6) then
                print *, "WARNING: Possible error in jacobian", jac_err
             end if
             
          end if

       end if
       
       ! Find norm of the residual
       abs_res_norm = norm2(res)
       if ( n .eq. 1) init_norm = abs_res_norm
       rel_res_norm = abs_res_norm/init_norm

       if ( this % print_level .eq. 2) then
          write(*, "(2I5,2ES12.2)") k, n, abs_res_norm, rel_res_norm
       end if

       ! Check stopping
       if ((abs_res_norm .le. this % atol) .or. (rel_res_norm .le. this % rtol)) then
          conv = .true.
          exit newton
       else if (abs_res_norm .ne. abs_res_norm .or. rel_res_norm .ne. rel_res_norm ) then
          conv = .false.
          exit newton
       end if

       ! Call LAPACK to solve the stage values system
       dq = -res
       call DGESV(size, 1, jac, size, IPIV, dq, size, INFO)
       
       ! Update the solution
       qddot = qddot + gamma * dq
       qdot  = qdot  + beta  * dq
       q     = q     + alpha * dq
       
    end do newton

    if (this % print_level .eq. 1) then 
       write(*, "(2I5,2ES12.2)") k, n, abs_res_norm, rel_res_norm
    end if

    ! Print warning message if not converged
    if (.not. conv) then
       write(*,'(/2A5, 2A12)') "STEP", "ITER", "|R|", "|R|/|R1|"
       write(*, "(2I5,2ES12.2)") k, n, abs_res_norm, rel_res_norm
       stop "Newton Solve Failed"
    end if

    if (allocated(ipiv)) deallocate(ipiv)
    if (allocated(res)) deallocate(res)
    if (allocated(dq)) deallocate(dq)
    if (allocated(jac)) deallocate(jac)
    if (allocated(fd_jac)) deallocate(fd_jac)
            
  end subroutine newtonSolve

  !===================================================================! 
  ! Routine that approximates the Jacobian based on finite differences
  ! [d{R}/d{q}] = alpha*[dR/dq] + beta*[dR/dqdot] + gamma*[dR/dqddot]
  !===================================================================!
  
  subroutine approximateJacobian( this, jac, alpha, beta, gamma, t, q, qdot, qddot )

    class(integrator)                       :: this
    
    ! Matrices
    real(dp), intent(inout), dimension(:,:) :: jac
    
    ! Arrays
    real(dp), intent(in)                    :: t
    real(dp), intent(in), dimension(:)      :: q, qdot, qddot     ! states

    real(dp), allocatable, dimension(:)     :: pstate             ! perturbed states
    real(dp), allocatable, dimension(:)     :: R, Rtmp            ! original residual and perturbed residual

    ! Scalars
    real(dp)                                :: dh = 1.0d-5        ! finite-diff step size
    real(dp), intent(in)                    :: alpha, beta, gamma ! linearization coefficients
    integer                                 :: m                  ! loop variables

    !  Zero the supplied jacobian matrix for safety (as we are
    !  computing everything newly here)
    jac = 0.0d0
    
    ! Allocate required arrays
    allocate(pstate(this % nsvars)); pstate = 0.0d0;
    allocate(R(this % nsvars));      R = 0.0d0;
    allocate(Rtmp(this % nsvars));   Rtmp = 0.0d0;

    ! Make a residual call with original variables
    call this % system % assembleResidual(R, t, q, qdot, qddot)

    !-----------------------------------------------------------!
    ! Derivative of R WRT Q: dR/dQ
    !-----------------------------------------------------------!

    pstate = q

    loop_vars: do m = 1, this % nsvars

       ! Perturb the k-th variable
       pstate(m) = pstate(m) + dh

       ! Make a residual call with the perturbed variable
       call this % system % assembleResidual(Rtmp, t, pstate, qdot, qddot)

       ! Unperturb (restore) the k-th variable
       pstate(m) =  q(m)

       ! Approximate the jacobian with respect to the k-th variable
       jac(:,m) = jac(:,m) + alpha*(Rtmp-R)/dh

    end do loop_vars

    !-----------------------------------------------------------!
    ! Derivative of R WRT QDOT: dR/dQDOT
    !-----------------------------------------------------------!

    pstate = qdot

    do m = 1, this % nsvars

       ! Perturb the k-th variable
       pstate(m) = pstate(m) + dh

       ! Make a residual call with the perturbed variable
       call this % system % assembleResidual(Rtmp, t, q, pstate, qddot)

       ! Unperturb (restore) the k-th variable
       pstate(m) =  qdot(m)

       ! Approximate the jacobian with respect to the k-th variable
       Jac(:,m) = Jac(:,m) + beta*(Rtmp-R)/dh

    end do

    ! Second order equations have an extra block to add
    if (this % second_order) then

       !-----------------------------------------------------------!
       ! Derivative of R WRT QDDOT: dR/dQDDOT
       !-----------------------------------------------------------!     

       pstate = qddot

       do m = 1, this % nsvars

          ! Perturb the k-th variable
          pstate(m) = pstate(m) + dh

          ! Make a residual call with the perturbed variable
          call this % system % assembleResidual(Rtmp, t, q, qdot, pstate)

          ! Unperturb (restore) the k-th variable
          pstate(m) =  qddot(m)

          ! Approximate the jacobian with respect to the k-th variable
          Jac(:,m) = Jac(:,m) + gamma*(Rtmp-R)/dh

       end do

    end if ! first or second order

    deallocate(pstate)
    deallocate(R,Rtmp)

  end subroutine approximateJacobian

end module integrator_class

!=====================================================================! 
! Backward Difference Formula Integration Module
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module bdf_integrator

  use precision, only : dp
  use integrator_class, only: integrator

  implicit none

  private
  public :: BDF

  !===================================================================! 
  ! BDF Integrator type
  !===================================================================! 

  type, extends(integrator) :: BDF

     integer                                :: max_bdf_order
     real(dp) , dimension(:), allocatable   :: beta, gamm
     real(dp) , dimension(:,:), allocatable :: psi, rhs

   contains
     
     private

     ! Routines for integration
     procedure, public :: initialize, finalize, integrate     
     procedure         :: approximateStates
     procedure         :: getOrderCoeff

     ! Routines for adjoint gradient
     procedure         :: adjointSolve, computeTotalDerivative
     procedure, public :: getAdjointGradient

  end type BDF

contains
  
  !===================================================================!
  ! Subroutine that integrates backwards in time to compute the
  ! lagrange multipliers (adjoint variables for the function)
  !===================================================================!
  
  subroutine adjointSolve(this)

    class(BDF)               :: this

  end subroutine adjointSolve

  !===================================================================!
  ! Compute the total derivative of the function with respect to the
  ! design variables and return the gradient 
  !===================================================================!

  subroutine computeTotalDerivative(this, x, dfdx)

    class(BDF)               :: this
    real(dp) , intent(in)    :: x
    real(dp) , intent(inout) :: dfdx

  end subroutine computeTotalDerivative

  !===================================================================!
  ! Public wrapper for all the adjoint gradient related sequence of
  ! calls
  !===================================================================!
  
  subroutine getAdjointGradient(this, x, dfdx)

    class(BDF)               :: this
    real(dp) , intent(in)    :: x
    real(dp) , intent(inout) :: dfdx

    !-----------------------------------------------------------------!
    ! Set the design variable into the system
    !-----------------------------------------------------------------!
    
    
    !-----------------------------------------------------------------!
    ! First integrate forward in time for the set design variable call
    !-----------------------------------------------------------------!
    
    call this % integrate()

    !-----------------------------------------------------------------!
    ! Integrate backwards to solve for lagrange multipliers for the
    ! set design variable
    !-----------------------------------------------------------------!

    call this % adjointSolve()

    !-----------------------------------------------------------------!
    ! Compute the total derivative of the function with respect to the
    ! design variables
    !-----------------------------------------------------------------!

    call this % computeTotalDerivative(x, dfdx)

  end subroutine getAdjointGradient

  !===================================================================!
  ! Initialize the BDF datatype and allocate required variables
  !===================================================================!
  
  subroutine initialize(this, nsvars, max_bdf_order, tinit, tfinal, h, second_order)
    
    class(BDF)                      :: this
    integer  , OPTIONAL, intent(in) :: max_bdf_order
    integer  , OPTIONAL, intent(in) :: nsvars
    real(dp) , OPTIONAL, intent(in) :: tinit, tfinal
    real(dp) , OPTIONAL, intent(in) :: h
    logical  , OPTIONAL, intent(in) :: second_order

    print *, "======================================"
    print *, ">>   Backward Difference Formula    << "
    print *, "======================================"

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
    ! Set the user supplied number of variables
    !-----------------------------------------------------------------!

    if (present(nsvars)) then
       this % nsvars = nsvars 
    end if
    print '("  >> Number of variables    : ",i4)', this % nsvars

    !-----------------------------------------------------------------!
    ! Find the number of time steps required during integration
    !-----------------------------------------------------------------!

    this % num_steps = int((this % tfinal - this % tinit)/this % h) + 1 
    print '("  >> Number of steps        : ",i6)', this % num_steps

    !-----------------------------------------------------------------!
    ! Allocate space for the RHS of adjoint equations
    !-----------------------------------------------------------------!

    allocate(this % RHS(this % num_steps, this % nsvars))
    this % RHS = 0.0d0

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

    if(allocated(this % RHS)) deallocate(this % RHS)
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
   
    call this % getOrderCoeff(m, this % beta, 1)

    do i = 1, m + 1
       this % udot(k,:) = this % udot(k,:) &
            & + this % beta(i)*this % u(k-i+1,:)/this % h
    end do
    
    !-----------------------------------------------------------------!
    ! Approximate UDDOT using BDF
    !-----------------------------------------------------------------!
    
    call this % getOrderCoeff(m, this % gamm, 2)

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
  
  subroutine getOrderCoeff(this, m, coeff, d)
    
    class(BDF)              :: this
    integer  , intent(in)   :: d
    integer  , intent(out)  :: m
    real(dp) , intent(out)  :: coeff(:)
    integer                 :: k

    k = this % current_step

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

end module bdf_integrator

!=====================================================================!
! A Diagonally Implicit Runge Kutta integrator module for first and
! second order systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module runge_kutta_integrator

  use precision, only : dp
  use integrator_class, only: integrator

  implicit none

  private
  public :: DIRK

  !===================================================================! 
  ! Abstract Runge-Kutta type
  !===================================================================! 

  type, abstract, extends(integrator) :: RK

     integer :: num_stages = 1 ! default number of stages
     integer :: order = 2      ! order of accuracy, only for informatory purposes
     integer :: current_stage = 0

     !----------------------------------------------------------------!
     ! The Butcher Tableau 
     !----------------------------------------------------------------!

     real(dp), dimension(:,:), allocatable :: A ! forms the coeff matrix
     real(dp), dimension(:), allocatable   :: B ! multiplies the state derivatives
     real(dp), dimension(:), allocatable   :: C ! multiplies the time

     !----------------------------------------------------------------!
     ! The stage time and its corresponding derivatives
     !----------------------------------------------------------------!

     real(dp), dimension(:), allocatable     :: T
     real(dp), dimension(:,:,:), allocatable :: Q
     real(dp), dimension(:,:,:), allocatable :: QDOT
     real(dp), dimension(:,:,:), allocatable :: QDDOT

     !----------------------------------------------------------------!
     ! The stage residual and jacobian
     !----------------------------------------------------------------!

     real(dp), dimension(:,:), allocatable     :: R ! stage residual
     real(dp), dimension(:,:,:,:), allocatable :: J ! stage jacobian

     real(dp), dimension(:,:,:), allocatable   :: psi, rhs

   contains

     !----------------------------------------------------------------!
     ! Implemented common procedures (visible to the user)
     !----------------------------------------------------------------!

     procedure :: initialize, finalize, integrate

     !----------------------------------------------------------------!
     ! Implemented procedures (not callable by the user)
     !----------------------------------------------------------------!

     procedure, private :: TimeMarch
     procedure, private :: CheckButcherTableau

     !----------------------------------------------------------------!
     ! Deferred common procedures
     !----------------------------------------------------------------!

     procedure(computeStageStateValues_interface), private, deferred :: ComputeStageStateValues
     procedure(buthcher_interface), private, deferred                :: SetupButcherTableau

  end type RK

  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(RK) :: DIRK

   contains

     private

     procedure :: setupButcherTableau => ButcherDIRK
     procedure :: computeStageStateValues

     !----------------------------------------------------------------!
     ! Adjoint procedures
     !----------------------------------------------------------------!

     procedure, public :: IntegrateBackward
     procedure, private :: AdjointSolve
     procedure, private :: AssembleRHS
     procedure, private :: AddFunctionDependency
     procedure, private :: AddTimeDependency
     procedure, private :: AddStageDependency

  end type DIRK

  !===================================================================!
  ! Interfaces for deferred specialized procedures 
  !===================================================================!

  interface

     !================================================================!
     ! Interface for finding the stage derivatives at each time step
     !================================================================!

     subroutine computeStageStateValues_interface(this, q, qdot)
       use precision, only : dp
       import RK
       class(RK) :: this
       real(dp), intent(in), dimension(:,:)           :: q
       real(dp), OPTIONAL, intent(in), dimension(:,:) :: qdot
     end subroutine computeStageStateValues_interface

     !================================================================!
     ! Interface for setting the Butcher tableau for each type of RK
     ! scheme
     !================================================================!

     subroutine buthcher_interface(this)
       use precision, only : dp
       import RK
       class(RK) :: this
     end subroutine buthcher_interface

  end interface

contains

  !===================================================================!
  ! Initialize the dirk datatype and construct the tableau
  !===================================================================!
  
  subroutine initialize(this, nsvars, num_stages, tinit, tfinal, h, second_order)
    
    class(RK)                       :: this
    integer  , OPTIONAL, intent(in) :: num_stages
    integer  , OPTIONAL, intent(in) :: nsvars
    real(dp) , OPTIONAL, intent(in) :: tinit, tfinal
    real(dp) , OPTIONAL, intent(in) :: h
    logical  , OPTIONAL, intent(in) :: second_order

    print *, "======================================"
    print *, ">> Diagonally-Implicit-Runge-Kutta  <<"
    print *, "======================================"

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

    if (present(num_stages)) then
       this % num_stages = num_stages
    end if
    print '("  >> Number of stages       : ",i4)', this % num_stages

    !-----------------------------------------------------------------!
    ! Set the user supplied initial step size
    !-----------------------------------------------------------------!

    if (present(h)) then
       this % h = h 
    end if
    print '("  >> Step size              : ",E9.3)', this % h
    
    !-----------------------------------------------------------------!
    ! Set the user supplied number of variables
    !-----------------------------------------------------------------!

    if (present(nsvars)) then
       this % nsvars = nsvars 
    end if
    print '("  >> Number of variables    : ",i4)', this % nsvars

    !-----------------------------------------------------------------!
    ! Find the number of time steps required during integration
    !-----------------------------------------------------------------!

    this % num_steps = int((this % tfinal - this % tinit)/this % h) + 1 
    print '("  >> Number of steps        : ",i6)', this % num_steps
    
    !-----------------------------------------------------------------!
    ! Allocate space for the tableau
    !-----------------------------------------------------------------!

    allocate(this % A(this % num_stages, this % num_stages))
    this % A = 0.0d0

    allocate(this % B(this % num_stages))    
    this % B = 0.0d0

    allocate(this % C(this % num_stages))
    this % C = 0.0d0

    !-----------------------------------------------------------------!
    ! Allocate space for the stage states and time
    !-----------------------------------------------------------------!

    allocate(this % T(this % num_stages))
    this % T = 0.0d0

    allocate(this % Q(this % num_steps, this % num_stages, this % nsvars))
    this % Q = 0.0d0

    allocate(this % QDOT(this % num_steps, this % num_stages, this % nsvars))
    this % QDOT = 0.0d0

    allocate(this % QDDOT(this % num_steps, this % num_stages, this % nsvars))
    this % QDDOT = 0.0d0

    !-----------------------------------------------------------------!
    ! Allocate space for the stage residual and jacobian
    !-----------------------------------------------------------------!

    allocate(this % R(this % num_stages, this % nsvars))
    this % R = 0.0d0

    allocate(this % J(this % num_stages, this % num_stages, this % nsvars, this % nsvars))
    this % J = 0.0d0

    !-----------------------------------------------------------------!
    ! Allocate space for the RHS of adjoint equations
    !-----------------------------------------------------------------!

    allocate(this % RHS(this % num_steps, this % num_stages, this % nsvars))
    this % RHS = 0.0d0

    allocate(this % psi(this % num_steps, this % num_stages, this % nsvars))
    this % psi = 0.0d0

    !-----------------------------------------------------------------!
    ! Allocate space for the global states and time
    !-----------------------------------------------------------------!

    allocate(this % time(this % num_steps))
    this % time = 0.0d0

    allocate(this % U(this % num_steps, this % nsvars))
    this % U = 0.0d0

    allocate(this % uDOT(this % num_steps, this % nsvars))
    this % UDOT = 0.0d0

    allocate(this % UDDOT(this % num_steps, this % nsvars))
    this % UDDOT = 0.0d0

    !-----------------------------------------------------------------!
    ! Put values into the Butcher tableau
    !-----------------------------------------------------------------!

    call this % setupButcherTableau()

    !-----------------------------------------------------------------!
    ! Sanity check for consistency of Butcher Tableau
    !-----------------------------------------------------------------!

    call this % checkButcherTableau()

    ! set the start time
    this % time(1) = this % tinit

  end subroutine initialize

  !===================================================================!
  ! Routine that checks if the Butcher Tableau entries are valid for
  ! the chosen number of stages/order
  !===================================================================!

  subroutine checkButcherTableau(this)

    class(RK) :: this
    integer :: i

    do i = 1, this  % num_stages
       if (abs(this % C(i) - sum(this % A(i,:))) .gt. 5.0d-16) then
          print *, "WARNING: sum(A(i,j)) != C(i)", i, this % num_stages
       end if
    end do

    if ((sum(this % B) - 1.0d0) .gt. 5.0d-16) then
       print *, "WARNING: sum(B) != 1", this % num_stages
    end if

  end subroutine checkButcherTableau

  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!

  subroutine finalize(this)

    class(RK) :: this

    ! Clear butcher's tableau
    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

    ! Clear stage value
    if(allocated(this % QDDOT)) deallocate(this % QDDOT)
    if(allocated(this % QDOT)) deallocate(this % QDOT)
    if(allocated(this % Q)) deallocate(this % Q)
    if(allocated(this % T)) deallocate(this % T)

    ! Clear the stage residual and jacobian
    if(allocated(this % R)) deallocate(this % R)
    if(allocated(this % J)) deallocate(this % J)

    ! Clear global states and time
    if(allocated(this % UDDOT)) deallocate(this % UDDOT)
    if(allocated(this % UDOT)) deallocate(this % UDOT)
    if(allocated(this % U)) deallocate(this % U)
    if(allocated(this % time)) deallocate(this % time)

    if(allocated(this % RHS)) deallocate(this % RHS)
    if(allocated(this % psi)) deallocate(this % psi)

  end subroutine finalize

  !===================================================================!
  ! Time integration logic
  !===================================================================!

  subroutine Integrate(this)

    class(RK) :: this
    integer :: k

    call this % system % getInitialStates(this % time(1), &
         & this % u(1,:), this % udot(1,:))
    
    this % current_step = 1

    ! March in time
    march: do k = 2, this % num_steps

       this % current_step = this % current_step + 1

       ! Find the stage derivatives at the current step
       call this % computeStageStateValues(this % u, this % udot)

       ! Advance the state to the current step
       call this % timeMarch(this % u, this % udot, this % uddot)

    end do march

  end subroutine Integrate
  
  !===================================================================!
  ! Time backwards in stage and backwards in time to solve for the
  ! adjoint variables
  !===================================================================!
  
  subroutine IntegrateBackward(this)
    
    class(DIRK) :: this
    integer :: k, i
    integer :: ndvars

    real(dp), allocatable, dimension(:) :: dfdx, dLdx, tmp
    real(dp), allocatable, dimension(:,:) :: dRdx

    allocate(dfdx(ndvars))
    allocate(dLdx(ndvars))
    allocate(dRdx(this%nsvars,ndvars))
    allocate(tmp(ndvars))

    do k = this % num_steps, 1, -1

       this % current_step = k

       do i = this % num_stages, 1, -1

          this % current_stage = i

          call this % AdjointSolve()

       end do

    end do

    ! Compute the total derivative
    tmp = 0.0d0
    do k = 1, this % num_steps
       tmp = tmp + matmul(this%psi(k,i,:), dRdx)
    end do

    ! call addDVSens
    dLdx = dfdx + tmp

    ! Write the adjoint variables 
    open(unit=90, file='output/adjoint.dat')
    do k = 1, this % num_steps
       write(90, *)  this % time(k), this % psi(k,1,:), this % psi(k,2,:), this % psi(k,3,:)
    end do
    close(90)
    
    deallocate(dfdx)
    deallocate(dLdx)
    deallocate(dRdx)
    deallocate(tmp)

  end subroutine IntegrateBackward

  !===================================================================!
  ! Update the states based on RK Formulae
  !===================================================================!

  subroutine timeMarch(this, q, qdot, qddot)

    implicit none

    class(RK) :: this
    real(dp),  dimension(:,:) :: q, qdot, qddot ! current state
    integer :: m, k

    ! Store the current time step
    k = this % current_step

    ! Increment the time
    this % time(k) = this % time(k-1) + this % h

    ! March q to next time step
    forall(m=1:this%nsvars)
       q(k,m) = q(k-1,m) + this % h*sum(this % B(1:this%num_stages) &
            &* this % QDOT(k, 1:this%num_stages, m))
    end forall

    if (this % second_order) then

       ! March qdot
       forall(m=1:this%nsvars)
          qdot(k,m) = qdot(k-1,m) + this % h*sum(this % B(1:this%num_stages) &
               &* this % QDDOT(k, 1:this%num_stages, m))
       end forall

       ! March qddot
       forall(m=1:this%nsvars)
          qddot(k,m) = sum(this % B(1:this%num_stages) &
               &* this % QDDOT(k,1:this%num_stages,m))
       end forall

    else

       ! March qdot
       forall(m=1:this%nsvars)
          qdot(k,m) = sum(this % B(1:this%num_stages) &
               & * this % QDOT(k,1:this%num_stages,m))
       end forall

    end if

  end subroutine timeMarch

  !===================================================================!
  ! Butcher's tableau for DIRK 
  !===================================================================!

  subroutine ButcherDIRK(this)

    class(DIRK) :: this
    real(dp), parameter :: PI = 22.0d0/7.0d0
    real(dp), parameter :: tmp  = 1.0d0/(2.0d0*dsqrt(3.0d0))
    real(dp), parameter :: half = 1.0d0/2.0d0
    real(dp), parameter :: one  = 1.0d0
    real(dp), parameter :: alpha = 2.0d0*cos(PI/18.0d0)/dsqrt(3.0d0)

    ! Put the entries into the tableau (ROGER ALEXANDER 1977)
    if (this % num_stages .eq. 1) then 

       ! Implicit mid-point rule (A-stable)

       this % A(1,1)    = half
       this % B(1)      = one
       this % C(1)      = half

       this % order     = 2

!!$       ! Implicit Euler (Backward Euler) but first order accurate
!!$       this % A(1,1) = one
!!$       this % B(1)   = one
!!$       this % C(1)   = one
!!$       this % order  = 1


    else if (this % num_stages .eq. 2) then

       ! Crouzeix formula (A-stable)

       this % A(1,1)    = half + tmp
       this % A(2,1)    = -one/dsqrt(3.0d0)
       this % A(2,2)    = this % A(1,1)

       this % B(1)      = half
       this % B(2)      = half

       this % C(1)      = half + tmp
       this % C(2)      = half - tmp

       this % order     = 3

    else if (this % num_stages .eq. 3) then

       ! Crouzeix formula (A-stable)

       this % A(1,1)    = (one+alpha)*half
       this % A(2,1)    = -half*alpha
       this % A(3,1)    =  one + alpha

       this % A(2,2)    = this % A(1,1)
       this % A(3,2)    = -(one + 2.0d0*alpha)
       this % A(3,3)    = this % A(1,1)

       this % B(1)      = one/(6.0d0*alpha*alpha)
       this % B(2)      = one - one/(3.0d0*alpha*alpha)
       this % B(3)      = this % B(1)

       this % C(1)      = (one + alpha)*half
       this % C(2)      = half
       this % C(3)      = (one - alpha)*half

       this % order     = 4

    else if (this % num_stages .eq. 4) then

       stop "Four stage DIRK formula does not exist"

    else
       
       print *, this % num_stages
       stop "DIRK Butcher tableau is not implemented for the requested&
            & order/stages"

    end if

  end subroutine ButcherDIRK

  !===================================================================!
  ! Get the stage derivative array for the current step and states for
  ! DIRK
  !===================================================================!

  subroutine computeStageStateValues(this, q, qdot)

    class(DIRK)                                    :: this
    real(dp), intent(in), dimension(:,:)           :: q
    real(dp), OPTIONAL, intent(in), dimension(:,:) :: qdot
    integer                                        :: k, j, m
    real(dp)                                       :: alpha, beta, gamma

    k = this % current_step

    this % current_stage = 0

    do j = 1, this % num_stages

       this % current_stage = this % current_stage + 1

       ! Find the stage times
       this % T(j) = this % time(k-1) + this % C(j)*this % h

       ! Guess the solution for stage states

       if (this % second_order) then
          
          ! guess qddot
          if ( k .eq. 2 .and. j .eq. 1) then
             this % QDDOT(k,j,:) = 1.0d0
          else if ( k .eq. 2 .and. j .gt. 1) then
             this % QDDOT(k,j,:) = this % UDDOT(k-1,:)
          else
             this % QDDOT(k,j,:) = this % QDDOT(k,j-1,:)
          end if

          ! compute the stage velocities for the guessed QDDOT
          forall(m = 1 : this % nsvars)
             this % QDOT(k,j,m) = qdot(k-1,m) &
                  & + this % h*sum(this % A(j,:)&
                  & * this % QDDOT(k,:, m))
          end forall

          ! compute the stage states for the guessed QDDOT
          forall(m = 1 : this % nsvars)
             this % Q(k,j,m) = q(k-1,m) &
                  & + this % h*sum(this % A(j,:)*this % QDOT(k,:, m))
          end forall

       else

          ! guess qdot
          if ( k .eq. 2 .and. j .eq. 1) then
             this % QDOT(k,j,:) = 1.0d0
          else if ( k .eq. 2 .and. j .gt. 1) then
             this % QDOT(k,j,:) = this % UDOT(k-1,:)
          else
             this % QDOT(k,j,:) = this % QDOT(k,j-1,:)
          end if

          ! compute the stage states for the guessed 
          forall(m = 1 : this % nsvars)
             this % Q(k,j,m) = q(k-1,m) &
                  & + this % h*sum(this % A(j,:)*this % QDOT(k,:, m))
          end forall

       end if
       
       ! solve the non linear stage equations using Newton's method for
       ! the actual stage states 
       if (this % second_order) then
          gamma = 1.0d0
          beta  = this % h * this%A(j,j)
          alpha = this % h * this%A(j,j)* this % h * this%A(j,j)
       else
          gamma = 0.0d0
          beta  = 1.0d0
          alpha = this % h * this%A(j,j)
       end if

       call this % newtonSolve(alpha, beta, gamma, &
            & this % time(k), this % q(k,j,:), this % qdot(k,j,:), this % qddot(k,j,:))
       
    end do

  end subroutine computeStageStateValues
  
  !===================================================================!
  ! Solve the linear adjoint equation at each stage and time step
  !===================================================================!
  
  subroutine AdjointSolve(this)

    class(DIRK)                           :: this
    real(dp), allocatable, dimension(:)   :: res, dq
    real(dp), allocatable, dimension(:,:) :: jac
    integer, allocatable, dimension(:)    :: ipiv
    integer                               :: info, size, k, i
    logical                               :: conv = .false.
    real(dp)                              :: alpha, beta, gamma

    k = this % current_step
    i = this % current_stage

    ! find the size of the linear system based on the calling object
    size = this % nsvars

    if (.not.allocated(ipiv)) allocate(ipiv(size))
    if (.not.allocated(res)) allocate(res(size))
    if (.not.allocated(dq)) allocate(dq(size))
    if (.not.allocated(jac)) allocate(jac(size,size))

    ! Get the residual of the function
    call this % assembleRHS(this % RHS(k,i,:))

    ! Get the jacobian matrix
    alpha = this % h * this % A(i,i)* this % h * this % A(i,i)
    beta  = this % h * this % A(i,i)
    gamma = 1.0d0
    call this % system % assembleJacobian(this % J(i, i,:,:), alpha, beta, gamma, &
         & this % T(i), this % Q(k,i,:), this % QDOT(k,i,:), this % QDDOT(k,i,:))

    ! Setup linear system in lapack format
    res = this % rhs(k,i,:)
    jac = transpose(this % J(i,i,:,:))

    ! Call lapack to solve the stage values system
    dq = -res
    call DGESV(size, 1, jac, size, IPIV, dq, size, INFO)

    this % psi(k,i,:) = dq

    if (allocated(ipiv)) deallocate(ipiv)
    if (allocated(res)) deallocate(res)
    if (allocated(dq)) deallocate(dq)
    if (allocated(jac)) deallocate(jac)

  end subroutine AdjointSolve

  !===================================================================!
  ! The derivative of the objective function with respect to the
  ! state.
  !===================================================================!
  
  subroutine AddFunctionDependency(this, rhs)

    class(DIRK) :: this
    real(dp), dimension(:) :: rhs
    real(dp) :: scale
    integer :: k, i, j

    i = this % current_stage
    k = this % current_step
    
    ! if (.not. k .eq. this % num_steps) return

    scale = 0.0d0
    do j = i, this % num_stages
       scale = scale + this % h * this % B(j) *  this % h * this %A(j,i)
    end do

    rhs(:) = scale*2.0d0*this % q(k,i,:)
    
    ! rhs(:) =  2.0d0 * this % u(k,:)
    ! rhs(:) =  2.0d0 * this % q(k,i,:)

  end subroutine AddFunctionDependency

  !===================================================================!
  ! Add the contribution to the RHS from the time-dependent terms
  !===================================================================!

  subroutine AddTimeDependency(this, rhs)
    
    class(DIRK) :: this
    real(dp), dimension(:) :: rhs !nsvars length
    integer  :: i, j, k, p
    real(dp) :: scal1, scal2
    real(dp), allocatable, dimension(:,:) :: jac1, jac2

    ! if this is the last step skip this term
    if (this % current_step .eq. this % num_steps) then
       ! print*, "no time dependent terms", this % current_step, this % num_steps
       return
    end if

    i = this % current_stage
    k = this % current_step + 1

    allocate(jac1(this%nsvars,this%nsvars))
    allocate(jac2(this%nsvars,this%nsvars))

    do j = 1, this % num_stages

       ! first term
       scal1 = this % h * this % B(i)
       jac1 = 0.0d0
       call this % system % assembleJacobian(jac1, 0.0d0, scal1, 0.0d0, &
            & this % T(j), this % Q(k,j,:), this % QDOT(k,j,:), this % QDDOT(k,j,:))

       ! second term
       scal2 = 0.0d0

       do p = i, this % num_stages
          scal2 = scal2 + this % h * this % h * this % A(p,i) * this % B(p)
       end do

       do p = 1, j
          scal2 = scal2 + this % h * this % h * this % B(i) * this % A(j,p)
       end do
       jac2 = 0.0d0
       call this % system % assembleJacobian(jac2, scal2, 0.0d0, 0.0d0, &
            & this % T(j), this % Q(k,j,:), this % QDOT(k,j,:), this % QDDOT(k,j,:))

       ! multiply the jacobian with the adjoint vector from other stages
       rhs = rhs + matmul(transpose(jac1+jac2), this % psi(k,j,:))

    end do

    if(allocated(jac1)) deallocate(jac1)
    if(allocated(jac2)) deallocate(jac2)

  end subroutine AddTimeDependency

  !===================================================================!
  ! Add the contribution to the RHS from the stage-dependent terms
  !===================================================================!

  subroutine AddStageDependency(this, rhs)

    class(DIRK) :: this
    real(dp), dimension(:) :: rhs !nsvars length
    integer  :: i, j, k, p
    real(dp) :: scal1, scal2
    real(dp), allocatable, dimension(:,:) :: jac1, jac2

    allocate(jac1(this%nsvars,this%nsvars))
    allocate(jac2(this%nsvars,this%nsvars))

    i = this % current_stage
    k = this % current_step

    do j = i+1, this % num_stages

       ! first term
       scal1 = this % h * this % A(j,i)

       call this % system % assembleJacobian(jac1, 0.0d0, scal1, 0.0d0, &
            & this % T(j), this % Q(k,j,:), this % QDOT(k,j,:), this % QDDOT(k,j,:))

       ! second term
       scal2 = 0.0d0
       do p = 1, j-i+1
          scal2 = scal2 + this % h * this % h * this % A(j,p+i-1) * this % A(p+i-1,i)
       end do

       call this % system % assembleJacobian(jac2, scal2, 0.0d0, 0.0d0, &
            & this % T(j), this % Q(k,j,:), this % QDOT(k,j,:), this % QDDOT(k,j,:))

       ! multiply the jacobian with the adjoint vector from other stages
       rhs = rhs + matmul(transpose(jac1+jac2), this % psi(k,j,:))

    end do

    if(allocated(jac1)) deallocate(jac1)
    if(allocated(jac2)) deallocate(jac2)

  end subroutine AddStageDependency
  
  !===================================================================!
  ! Function that puts together the right hand side of the adjoint
  ! equation into the supplied rhs vector.
  !===================================================================!

  subroutine assembleRHS(this, rhs)
 
    class(DIRK) :: this
    integer     :: k , i
    real(dp)    :: rhs(:) ! of length nsvars

    k = this % current_step
    i = this % current_stage
    
    ! Add the contributions from the objective function
    call this % AddFunctionDependency(rhs)
    call this % AddTimeDependency(rhs)
    call this % AddStageDependency(rhs)

  end subroutine assembleRHS

end module runge_kutta_integrator
