#include "scalar.fpp"

!=====================================================================!
! Parent class for integration schemes to extend. This has some common
! logic such as:
!
! (1) nonlinear solution process
! (2) approximating derivatives 
! (3) writing solution to files
! (4) adjoint system solution
!
!=====================================================================!

module integrator_class

  use physics_class, only : physics

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
     integer                 :: nsvars = 0 ! number of states/equations

     !----------------------------------------------------------------!
     ! Variables for managing time marching
     !----------------------------------------------------------------!

     integer  :: num_steps = 0     ! number of time steps
     real(dp) :: tinit     = 0.0d0 ! initial time
     real(dp) :: tfinal    = 1.0d0 ! final time
     real(dp) :: h         = 0.1d0 ! default step size

     !----------------------------------------------------------------!
     ! Track global time and states
     !----------------------------------------------------------------!

     real(dp), dimension(:), allocatable       :: time
     type(scalar), dimension(:,:), allocatable :: U
     type(scalar), dimension(:,:), allocatable :: UDOT
     type(scalar), dimension(:,:), allocatable :: UDDOT

     !----------------------------------------------------------------!
     ! Miscellaneous variables
     !----------------------------------------------------------------!

     logical :: second_order         = .true.
     logical :: forward              = .true.
     integer :: print_level          = 0
     integer :: current_step         = 0
     logical :: approximate_jacobian = .false.

   contains

     !----------------------------------------------------------------!
     ! Deferred procedures for subtypes to implement                  !
     !----------------------------------------------------------------!

     procedure(InterfaceDefault)        , private, deferred :: approximateStates
     procedure(InterfaceGetLinearCoeff) , private, deferred :: getLinearCoeff

     !----------------------------------------------------------------!
     ! Procedures                                                     !
     !----------------------------------------------------------------!

     procedure :: writeSolution
     procedure :: setPhysicalSystem 
     procedure :: setPrintLevel
     procedure :: setApproximateJacobian
     procedure :: construct, destruct
     procedure :: integrate
     
  end type integrator

  interface

     !===================================================================!
     ! Default interface without any arguments
     !===================================================================!

     subroutine InterfaceDefault(this)

       import integrator

       class(integrator) :: this

     end subroutine InterfaceDefault

     !===================================================================!
     ! Interface for getting the coefficients for Residual linearization
     !===================================================================!

     subroutine InterfaceGetLinearCoeff(this, alpha, beta, gamma)

       import integrator

       class(integrator)         :: this
       type(scalar), intent(out) :: alpha, beta, gamma

     end subroutine InterfaceGetLinearCoeff

  end interface

contains

  !======================================================================!
  ! Base class constructor logic
  !======================================================================!

  subroutine construct(this, system, tinit, tfinal, h, second_order)

    class(integrator)               :: this
    class(physics), target          :: system
    real(dp), OPTIONAL, intent(in)  :: tinit, tfinal
    real(dp), OPTIONAL, intent(in)  :: h
    logical , OPTIONAL, intent(in)  :: second_order

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
    print '("  >> Number of steps        : ",i10)', this % num_steps

    print '("  >> Physical System        : ",A10)', this % system % name

    !-----------------------------------------------------------------!
    ! Allocate space for the global states and time
    !-----------------------------------------------------------------!

    allocate(this % time(this % num_steps))
    this % time = 0.0d0

    ! Set the start time
    this % time(1) = this % tinit

    allocate(this % U(this % num_steps, this % nsvars))
    this % U = 0.0d0

    allocate(this % UDOT(this % num_steps, this % nsvars))
    this % UDOT = 0.0d0

    allocate(this % UDDOT(this % num_steps, this % nsvars))
    this % UDDOT = 0.0d0

  end subroutine construct

  !======================================================================!
  ! Base class destructor
  !======================================================================!

  subroutine destruct(this)

    class(integrator) :: this

    ! Clear global states and time
    if(allocated(this % UDDOT)) deallocate(this % UDDOT)
    if(allocated(this % UDOT)) deallocate(this % UDOT)
    if(allocated(this % U)) deallocate(this % U)
    if(allocated(this % time)) deallocate(this % time)

  end subroutine destruct

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

    class(integrator)                      :: this
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
       write(90, *)  this % time(k), &
            & (dble(this % u     (k,j) ), j=1,this%nsvars ), &
            & (dble(this % udot  (k,j) ), j=1,this%nsvars ), &
            & (dble(this % uddot (k,j) ), j=1,this%nsvars )
    end do

    close(90)

  end subroutine writeSolution

  !===================================================================!
  ! Time integration logic
  !===================================================================!

  subroutine Integrate( this )
    
    use nonlinear_algebra, only: nonlinear_solve

    class(integrator) :: this
    type(scalar)      :: alpha, beta, gamma
    integer           :: k
    type(scalar)      :: coeff(3)

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

       ! Increment the time
       this % time(k) = this % time(k-1) + this % h

       ! Approximate the states u, udot and uddot using ABM stencil
       call this % approximateStates()

       ! Determine the coefficients for linearing the Residual
       call this % getLinearCoeff(alpha, beta, gamma)

       ! Solve the nonlinear system at each step by driving the
       ! residual to zero
       call nonlinear_solve(this % system, &
            & alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))
       
    end do time

  end subroutine Integrate
  
end module integrator_class
