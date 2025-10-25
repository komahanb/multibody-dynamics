!=====================================================================!
! A Diagonally Implicit Runge Kutta integrator module for first and
! second order systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!
module runge_kutta_integrator

  use iso_fortran_env, only: dp => real64
  use physics_class

  implicit none

  private

  public :: DIRK

  !-------------------------------------------------------------------!
  ! Abstract Runge-Kutta type
  !-------------------------------------------------------------------!

  type, abstract :: RK

     !----------------------------------------------------------------!
     ! Contains the actual physical system
     !----------------------------------------------------------------!

     class(physics), pointer :: system => null()

     integer :: num_stages = 1  ! default number of stages
     integer :: nvars = 1       ! number of states/equations
     integer :: order           ! order of accuracy
     integer :: num_steps       ! number of time steps

     real(8) :: tinit = 0.0d0, tfinal = 1.0d0
     real(8) :: h = 0.1d0       ! default step size

     logical :: second_order = .false.

     integer :: current_stage = 0
     integer :: current_step  = 1

     logical :: approximate_jacobian = .true. ! by default approximate
                                              ! Jacobian using finite
                                              ! differences

     !----------------------------------------------------------------!
     ! Track global time and states
     !----------------------------------------------------------------!

     real(8), dimension(:), allocatable   :: time
     real(8), dimension(:,:), allocatable :: U
     real(8), dimension(:,:), allocatable :: UDOT
     real(8), dimension(:,:), allocatable :: UDDOT

     !----------------------------------------------------------------!
     ! The Butcher Tableau 
     !----------------------------------------------------------------!

     real(8), dimension(:,:), allocatable :: A ! forms the coeff matrix
     real(8), dimension(:), allocatable   :: B ! multiplies the state derivatives
     real(8), dimension(:), allocatable   :: C ! multiplies the time

     !----------------------------------------------------------------!
     ! The stage time and its corresponding derivatives
     !----------------------------------------------------------------!

     real(8), dimension(:), allocatable   :: T
     real(8), dimension(:,:), allocatable :: Q
     real(8), dimension(:,:), allocatable :: QDOT
     real(8), dimension(:,:), allocatable :: QDDOT

     !----------------------------------------------------------------!
     ! The stage residual and jacobian
     !----------------------------------------------------------------!

     real(8), dimension(:,:), allocatable     :: R ! stage residual
     real(8), dimension(:,:,:,:), allocatable :: J ! stage jacobian

   contains

     !----------------------------------------------------------------!
     ! Implemented common procedures (visible to the user)
     !----------------------------------------------------------------!

     procedure :: initialize, finalize, integrate, write_solution

     !----------------------------------------------------------------!
     ! Implemented procedures (not callable by the user)
     !----------------------------------------------------------------!

     procedure, private :: time_march
     procedure, private :: reset_stage_values
     procedure, private :: check_butcher_tableau

     !----------------------------------------------------------------!
     ! Deferred common procedures
     !----------------------------------------------------------------!

     procedure(compute_stage_values_interface), private, deferred :: compute_stage_values
     procedure(buthcher_interface), private, deferred :: setup_butcher_tableau

     !----------------------------------------------------------------!
     ! Important setters
     !----------------------------------------------------------------!

     procedure :: setPhysicalSystem
     procedure :: setApproximateJacobian
     
  end type RK

  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(RK) :: DIRK

     !----------------------------------------------------------------!
     ! Nonlinear solution at each stage
     !----------------------------------------------------------------!
     
     integer :: max_newton = 25
     real(8) :: tol = 1.0d-12

   contains

     private

     !----------------------------------------------------------------!
     ! Implement/override the abstract class routines
     !----------------------------------------------------------------!

     procedure :: setup_butcher_tableau => ButcherDIRK

     !----------------------------------------------------------------!
     ! More specialized procedures
     !----------------------------------------------------------------!

     procedure :: compute_stage_values

     procedure :: newton_solve
     procedure :: state_update
     procedure :: setup_linear_system

     procedure :: get_residual
     procedure :: approx_jacobian

     procedure :: check_jacobian

  end type DIRK

  !-------------------------------------------------------------------!
  ! Interfaces for deferred specialized procedures 
  !-------------------------------------------------------------------!
  
  interface

     !----------------------------------------------------------------!
     ! Interface for finding the stage derivatives at each time step
     !----------------------------------------------------------------!

     subroutine compute_stage_values_interface(this, q, qdot)
       import RK
       class(RK) :: this
       real(8), intent(in), dimension(:,:) :: q
       real(8), OPTIONAL, intent(in), dimension(:,:) :: qdot
     end subroutine compute_stage_values_interface

     !----------------------------------------------------------------!
     ! Interface for setting the Butcher tableau for each type of RK
     ! scheme
     !----------------------------------------------------------------!

     subroutine buthcher_interface(this)
       import RK
       class(RK) :: this
     end subroutine buthcher_interface

  end interface

contains
  
  
  !--------------------------------------------------------------------!
  ! Setter that can be used to set the method in which jacobian needs
  ! to be computed. Setting this to .true. would make the code use
  ! finite differences, this is enabled by default too. If set to
  ! .false. the expects to provide implementation in assembleJacobian
  ! in a type that extends PHYSICS.
  !--------------------------------------------------------------------!

  subroutine setApproximateJacobian(this, approx_jacobian)

    class(RK) :: this
    logical :: approx_jacobian

    this % approximate_jacobian = approx_jacobian

  end subroutine setApproximateJacobian
  
  !--------------------------------------------------------------------!
  ! Set ANY physical system that extends the type PHYSICS and provides
  ! implementation to the mandatory functions assembleResidual and
  ! getInitialStates
  ! --------------------------------------------------------------------!

  subroutine setPhysicalSystem(this, physical_system)
    
    class(RK) :: this
    class(physics), target :: physical_system

    this % system => physical_system

  end subroutine setPhysicalSystem

  !-------------------------------------------------------------------!
  ! Initialize the dirk datatype and construct the tableau
  !-------------------------------------------------------------------!
  
  subroutine initialize(this, nvars, num_stages, tinit, tfinal, h, second_order)
    
    class(RK) :: this
    integer, OPTIONAL, intent(in) :: num_stages
    integer, OPTIONAL, intent(in) :: nvars
    real(8), OPTIONAL, intent(in) :: tinit, tfinal
    real(8), OPTIONAL, intent(in) :: h
    logical, OPTIONAL, intent(in) :: second_order

    !-----------------------------------------------------------------!
    ! Set the order of the governing equations
    !-----------------------------------------------------------------!
    
    if (present(second_order)) then
       this % second_order = second_order
    end if
    print '(">> Second order           : ",L1)', this % second_order

    !-----------------------------------------------------------------!
    ! Set the initial and final time
    !-----------------------------------------------------------------!

    if (present(tinit)) then
       this % tinit = tinit
    end if
    print '(">> Start time             : ",F8.3)', this % tinit

    if (present(tfinal)) then
       this % tfinal = tfinal
    end if
    print '(">> End time               : ",F8.3)', this % tfinal

    !-----------------------------------------------------------------!
    ! Set the order of integration
    !-----------------------------------------------------------------!

    if (present(num_stages)) then
       this % num_stages = num_stages
    end if
    print '(">> Number of stages       : ",i4)', this % num_stages

    !-----------------------------------------------------------------!
    ! Set the user supplied initial step size
    !-----------------------------------------------------------------!

    if (present(h)) then
       this % h = h 
    end if
    print '(">> Step size              : ",E9.3)', this % h
    
    !-----------------------------------------------------------------!
    ! Set the user supplied number of variables
    !-----------------------------------------------------------------!

    if (present(nvars)) then
       this % nvars = nvars 
    end if
    print '(">> Number of variables    : ",i4)', this % nvars

    !-----------------------------------------------------------------!
    ! Find the number of time steps required during integration
    !-----------------------------------------------------------------!

    this % num_steps = int((this % tfinal - this % tinit)/this % h) + 1 
    print '(">> Number of steps        : ",i6)', this % num_steps

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

    allocate(this % Q(this % num_stages, this % nvars))
    this % Q = 0.0d0

    allocate(this % QDOT(this % num_stages, this % nvars))
    this % QDOT = 0.0d0

    allocate(this % QDDOT(this % num_stages, this % nvars))
    this % QDDOT = 0.0d0

    !-----------------------------------------------------------------!
    ! Allocate space for the stage residual and jacobian
    !-----------------------------------------------------------------!

    allocate(this % R(this % num_stages, this % nvars))
    this % R = 0.0d0

    allocate(this % J(this % num_stages, this % num_stages, this % nvars, this % nvars))
    this % J = 0.0d0

    !-----------------------------------------------------------------!
    ! Allocate space for the global states and time
    !-----------------------------------------------------------------!

    allocate(this % time(this % num_steps))
    this % time = 0.0d0

    allocate(this % U(this % num_steps, this % nvars))
    this % U = 0.0d0

    allocate(this % uDOT(this % num_steps, this % nvars))
    this % UDOT = 0.0d0

    allocate(this % UDDOT(this % num_steps, this % nvars))
    this % UDDOT = 0.0d0

    !-----------------------------------------------------------------!
    ! Put values into the Butcher tableau
    !-----------------------------------------------------------------!

    call this % setup_butcher_tableau()

    !-----------------------------------------------------------------!
    ! Sanity check for consistency of Butcher Tableau
    !-----------------------------------------------------------------!

    call this % check_butcher_tableau()

    ! set the start time
    this % time(1) = this % tinit

  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Routine that checks if the Butcher Tableau entries are valid for
  ! the chosen number of stages/order
  !--------------------------------------------------------------------!

  subroutine check_butcher_tableau(this)

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

  end subroutine check_butcher_tableau

  !-------------------------------------------------------------------!
  ! Deallocate the tableau entries
  !-------------------------------------------------------------------!

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

  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Time integration logic
  !-------------------------------------------------------------------!
  ! Input: 
  ! o state arrays q and qdot with initial conditions set at q (1)
  !-------------------------------------------------------------------!
  ! Output:
  ! o q, qdot arrays are modified by the routine
  !-------------------------------------------------------------------!

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
       call this % compute_stage_values(this % u, this % udot)

       ! Advance the state to the current step
       call this % time_march(this % u, this % udot, this % uddot)

    end do march

  end subroutine Integrate

  !-------------------------------------------------------------------!
  ! Write solution to file
  !-------------------------------------------------------------------!

  subroutine write_solution(this)

    class(RK) :: this
    integer   :: k, j

    open(unit=90, file='solution.dat')
    
    do k = 1, this % num_steps
       write(90, *)  this % time(k), (this % u(k,j), j=1,this%nvars ), &
            & (this % udot(k,j), j=1,this%nvars ), &
            & (this % uddot(k,j), j=1,this%nvars )
    end do

    close(90)

    ! exact_solution(this % time(k),1.0d0,0.0d0)
    
  end subroutine write_solution

  !-------------------------------------------------------------------!
  ! Update the states based on RK Formulae
  !-------------------------------------------------------------------!

  subroutine time_march(this, q, qdot, qddot)

    implicit none

    class(RK) :: this
    real(8),  dimension(:,:) :: q, qdot, qddot ! current state
    integer :: m, k

    ! Store the current time step
    k = this % current_step

    ! Increment the time
    this % time(k) = this % time(k-1) + this % h

    ! March q to next time step
    forall(m=1:this%nvars)
       q(k,m) = q(k-1,m) + this % h*sum(this % B(1:this%num_stages) &
            &* this % QDOT(1:this%num_stages,m))
    end forall

    if (this % second_order) then

       ! March qdot
       forall(m=1:this%nvars)
          qdot(k,m) = qdot(k-1,m) + this % h*sum(this % B(1:this%num_stages) &
               &* this % QDDOT(1:this%num_stages,m))
       end forall

       ! March qddot
       forall(m=1:this%nvars)
          qddot(k,m) = sum(this % B(1:this%num_stages) &
               &* this % QDDOT(1:this%num_stages,m))
       end forall

    else

       ! March qdot
       forall(m=1:this%nvars)
          qdot(k,m) = sum(this % B(1:this%num_stages) &
               & * this % QDOT(1:this%num_stages,m))
       end forall

    end if

  end subroutine time_march

  !-------------------------------------------------------------------!
  ! Reset the array to store new stage values at each time step
  !-------------------------------------------------------------------!

  subroutine reset_stage_values(this)

    class(RK) :: this

    ! reset the variables that are computed during each time step    
    this % QDDOT = 0.0d0
    this % QDOT = 0.0d0
    this % Q = 0.0d0
    this % T = 0.0d0

    this % R = 0.0d0
    this % J = 0.0d0

  end subroutine reset_stage_values

  !-------------------------------------------------------------------!
  ! Butcher's tableau for DIRK 
  !-------------------------------------------------------------------!

  subroutine ButcherDIRK(this)

    class(DIRK) :: this
    real(8), parameter :: PI = 22.0d0/7.0d0
    real(8), parameter :: tmp  = 1.0d0/(2.0d0*dsqrt(3.0d0))
    real(8), parameter :: half = 1.0d0/2.0d0
    real(8), parameter :: one  = 1.0d0
    real(8), parameter :: alpha = 2.0d0*cos(PI/18.0d0)/dsqrt(3.0d0)

    ! put the entries into the tableau (ROGER ALEXANDER 1977)
    if (this % num_stages .eq. 1) then 

       ! Implicit mid-point rule (A-stable)

       this % A(1,1) = half
       this % B(1)   = one
       this % C(1)   = half

       this % order = 2

!!$       ! Implicit Euler (Backward Euler) but first order accurate
!!$       this % A(1,1) = one
!!$       this % B(1)   = one
!!$       this % C(1)   = one
!!$       this % order = 1


    else if (this % num_stages .eq. 2) then

       ! Crouzeix formula (A-stable)

       this % A(1,1) = half + tmp
       this % A(2,1) = -one/dsqrt(3.0d0)
       this % A(2,2) = this % A(1,1)

       this % B(1)   = half
       this % B(2)   = half

       this % C(1)   = half + tmp
       this % C(2)   = half - tmp

       this % order = 3

    else if (this % num_stages .eq. 3) then

       ! Crouzeix formula (A-stable)

       this % A(1,1) = (one+alpha)*half
       this % A(2,1) = -half*alpha
       this % A(3,1) =  one + alpha

       this % A(2,2) = this % A(1,1)
       this % A(3,2) = -(one + 2.0d0*alpha)
       this % A(3,3) = this % A(1,1)

       this % B(1)   = one/(6.0d0*alpha*alpha)
       this % B(2)   = one - one/(3.0d0*alpha*alpha)
       this % B(3)   = this % B(1)

       this % C(1) = (one + alpha)*half
       this % C(2) = half
       this % C(3) = (one - alpha)*half

       this % order = 4

    else if (this % num_stages .eq. 4) then

       stop "Four stage DIRK formula does not exist"

    else

       print *, this % num_stages
       stop "DIRK Butcher tableau is not implemented for the requested&
            & order/stages"

    end if

  end subroutine ButcherDIRK

  !-------------------------------------------------------------------!
  ! Get the stage derivative array for the current step and states for
  ! DIRK
  !-------------------------------------------------------------------!

  subroutine compute_stage_values(this, q, qdot)

    class(DIRK) :: this
    real(8), intent(in), dimension(:,:) :: q
    real(8), OPTIONAL, intent(in), dimension(:,:) :: qdot
    integer :: k, j, m

    ! set the stage values to zero
    call this % reset_stage_values()
    
    k = this % current_step

    this % current_stage = 0

    do j = 1, this % num_stages

       this % current_stage = this % current_stage + 1

       ! Find the stage times
       this % T(j) = this % time(k-1) + this % C(j)*this % h

       ! Guess the solution for stage states

       if (this % second_order) then

          ! guess qddot
          this % QDDOT(j,:) = 1.0d0 

          ! compute the stage velocities for the guessed QDDOT
          forall(m = 1 : this % nvars)
             this % QDOT(j,m) = qdot(k-1,m) &
                  & + this % h*sum(this % A(j,:)&
                  & * this % QDDOT(:, m))
          end forall

          ! compute the stage states for the guessed QDDOT
          forall(m = 1 : this % nvars)
             this % Q(j,m) = q(k-1,m) &
                  & + this % h*sum(this % A(j,:)*this % QDOT(:, m))
          end forall

       else

          ! guess qdot
          this % QDOT(j,:) = 1.0d0 

          ! compute the stage states for the guessed 
          forall(m = 1 : this % nvars)
             this % Q(j,m) = q(k-1,m) &
                  & + this % h*sum(this % A(j,:)*this % QDOT(:, m))
          end forall

       end if

       ! solve the non linear stage equations using Newton's method for
       ! the actual stage states 
       call this % newton_solve()

    end do

  end subroutine compute_stage_values

  !-------------------------------------------------------------------!
  ! Solve nonlinear stage equations using Newton's method at each time
  ! step.
  !
  ! q_{k,i} = q_{k} + h \sum_{j=1}^s {a_{i}j f(t_{k,j}, q_{k,j})
  ! i = 1,\ldots,s 
  !
  ! This yields $s$ equations and $s$ unknown stage values, $q_{k,i}$,
  ! that are solved using Newton's method at each time step
  ! -------------------------------------------------------------------!

  subroutine newton_solve(this)

    class(DIRK) :: this
    real(8), allocatable, dimension(:)   :: res, dq
    real(8), allocatable, dimension(:,:) :: jac
    integer, allocatable, dimension(:)   :: ipiv
    integer :: n, info, size, k, j
    logical :: conv = .false.
    real(8) :: alpha, beta, gamma

    k = this % current_step
    j = this % current_stage

    ! find the size of the linear system based on the calling object
    size = this % nvars

    if (.not.allocated(ipiv)) allocate(ipiv(size))
    if (.not.allocated(res)) allocate(res(size))
    if (.not.allocated(dq)) allocate(dq(size))
    if (.not.allocated(jac)) allocate(jac(size,size))

    newton: do n = 1, this % max_newton
       
       ! Get the residual of the function
       this % R = 0.0d0
       call this % system % assembleResidual( this % R(j,:), this % T(j), &
            & this % Q(j,:), &
            & this % QDOT(j,:), &
            & this % QDDOT(j,:))

       ! Get the jacobian matrix
       this % J = 0.0d0

       alpha = this % h * this % A(j,j)* this % h * this % A(j,j)
       beta  = this % h * this % A(j,j)
       gamma = 1.0d0
       
       if (this % approximate_jacobian) then

          ! Use finite difference to approximate the Jacobian
          call this % approx_jacobian()

       else

          ! Use the analytical Jacobian the user provided
          call this % system % assembleJacobian(this % J(j, j,:,:),&
               & alpha, beta, gamma, &
               & this % T(j), &
               & this % Q(j,:), &
               & this % QDOT(j,:), &
               & this % QDDOT(j,:))

          ! Check the jacobian implementation once at the beginning of integration
          if (this % current_step .eq. 2 .and. this % current_stage .eq. 1 .and. n .eq. 1 ) then
           
             ! print *, ">> Checking Jacobian Implementation..."

             call this % check_jacobian(this % current_stage, &
                  & this % J(this % current_stage,this % current_stage,:,:) )

          end if

       end if

       ! setup linear system in lapack format
       call this % setup_linear_system(res, jac)

       ! check stopping
       if (norm2(res) .le. this % tol) then
          conv = .true.
          exit newton
       end if

       ! call lapack to solve the stage values system
       dq = -res
       call DGESV(size, 1, jac, size, IPIV, dq, size, INFO)

       ! check stopping
       if (norm2(dq) .le. this % tol) then
          conv = .true.
          exit newton
       end if

       ! update the solution
       call this % state_update(dq)

    end do newton

    ! print warning message if not converged
    if (.not. conv) then
       
       print '("Newton solve: step = ", i3 , " iters = ", i3,&
            & " |R| = ",E10.3," |dq| = ",E10.3)',&
            & k, n, norm2(res), norm2(dq)

       stop "Newton Solve Failed"
    end if

    if (allocated(ipiv)) deallocate(ipiv)
    if (allocated(res)) deallocate(res)
    if (allocated(dq)) deallocate(dq)
    if (allocated(jac)) deallocate(jac)

  end subroutine newton_solve

  !-------------------------------------------------------------------!
  ! Routine that packs the matrix in a form that is used in lapack
  !-------------------------------------------------------------------!

  subroutine setup_linear_system(this, res, jac)

    implicit none

    class(dirk) :: this

    real(8), intent(inout), dimension(:) :: res
    real(8), intent(inout), dimension(:,:) :: jac

    res = this % R(this % current_stage,:)
    jac = this % J(this % current_stage, this % current_stage, :, :)

  end subroutine setup_linear_system

  !-------------------------------------------------------------------!
  ! After the solution of stage equations, we update the states using
  ! this call
  ! -------------------------------------------------------------------!

  subroutine state_update(this, sol)

    implicit none

    class(dirk) :: this
    real(8) :: sol(:)
    integer :: i

    i = this % current_stage

    if (this % second_order) then

       ! update qddot          
       this % QDDOT(i,:) = this % QDDOT(i,:) + sol(:)

       ! update qdot
       this % QDOT(i,:) = this % QDOT(i,:) &
            & + this % h * this % A(i,i) * sol(:)

       ! update q
       this % Q(i,:) = this % Q(i,:) &
            & + this % h * this % A(i,i) * this % h * this % A(i,i) &
            & * sol(:)

    else

       ! update qdot(k,i) for i-th stage
       this % QDOT(i,:) = this % QDOT(i,:) + sol(:)

       ! update q for i-th stage
       this % Q(i,:) = this % Q(i,:) &
            & + this % h * this % A(i,i)*sol(:)

    end if


  end subroutine state_update

  !-------------------------------------------------------------------!
  ! Computes the stage residual for the set stage state Y (comes from
  ! Newton's iteration) and sets into the same instance
  !
  ! R_{i}= q_{k,i} - q_{k} - h \sum_{j=1}^s {a_{ij} f(t_{k,j}, q_{k,j})
  ! i = 1,\ldots,s 
  !-------------------------------------------------------------------!
  
  subroutine get_residual(this)

    class(DIRK) :: this
    integer ::  m
    external :: R

    if (this % second_order) then

       ! compute the stage residuals for Q, QDOT, QDDOT
       call R(this % R(this % current_stage,:), this % nvars, &
            & this % T(this % current_stage), this % Q(this % current_stage,:), &
            & this % QDOT(this % current_stage,:), this % QDDOT(this % current_stage,:))

    else 

       ! compute the stage residuals
       call R(this % R(this % current_stage,:), this % nvars, &
            & this % T(this % current_stage), this % Q(this % current_stage,:), &
            & this % QDOT(this % current_stage,:))

    end if

  end subroutine get_residual

  !-------------------------------------------------------------------!
  ! Computes the stage jacobian and sets into the same instance
  !          J(i,j) = [ 1 - h A(i,j) DFDQ(T(j),Y(j))]
  !-------------------------------------------------------------------!

  subroutine approx_jacobian(this)

    class(DIRK) :: this
    integer :: i
    real(8) :: alpha
    external :: DRDQ, DRDQDOT, DRDQDDOT
    
    i = this % current_stage
!!$    
!!$    if (this % second_order) then
!!$
!!$       ! get the q block
!!$       alpha =  this % h * this % A(i,i)* this % h * this % A(i,i)
!!$       call DRDQ(this % J(i,i,:,:), alpha, this % nvars, this % T(i), &
!!$            & this % Q(i,:), this % QDOT(i,:), this % QDDOT(i,:))
!!$
!!$       ! get the qdot block
!!$       alpha =  this % h * this % A(i,i)
!!$       call DRDQDOT(this % J(i,i,:,:), alpha, this % nvars, this % T(i), &
!!$            & this % Q(i,:), this % QDOT(i,:), this % QDDOT(i,:))
!!$
!!$       ! get the qddot block
!!$       alpha = 1.0d0
!!$       call DRDQDDOT(this % J(i,i,:,:), alpha, this % nvars, this % T(i), &
!!$            & this % Q(i,:), this % QDOT(i,:), this % QDDOT(i,:))
!!$
!!$    else
!!$
!!$       ! get the q block
!!$       alpha = this % h * this % A(i,i)
!!$       call DRDQ(this % J(i,i,:,:), alpha, this % nvars, this % T(i), &
!!$            & this % Q(i,:), this % QDOT(i,:))
!!$
!!$       ! get the qdot block
!!$       alpha = 1.0d0
!!$       call DRDQDOT(this % J(i,i,:,:), alpha, this % nvars, this % T(i), &
!!$            & this % Q(i,:), this % QDOT(i,:))
!!$
!!$    end if

    ! check with FD
    call this % check_jacobian(i, this % J(i,i,:,:))

  end subroutine approx_jacobian

  !-------------------------------------------------------------------!  
  ! Routine to sanity check the jacobian of the governing equations
  !-------------------------------------------------------------------!

  subroutine check_jacobian(this, i, exact_jac)

    class(dirk) :: this

    integer, intent(in) :: i
    real(8), intent(inout) :: exact_jac(:,:)

    real(8), allocatable, dimension(:) :: tmp1, tmp2, qtmp, qdottmp, qddottmp
    real(8), allocatable, dimension(:,:) :: jtmp1, jtmp2, jtmp, jtmp3
    real(8) :: small = 1.0d-6
    integer :: k

    allocate(qtmp(this % nvars)); qtmp = 0.0d0;
    allocate(qdottmp(this % nvars)); qdottmp = 0.0d0;
    allocate(qddottmp(this % nvars)); qddottmp = 0.0d0;

    allocate(tmp1(this % nvars)); tmp1 = 0.0d0;
    allocate(tmp2(this % nvars)); tmp2 = 0.0d0;

    allocate(jtmp (this % nvars, this % nvars)); jtmp = 0.0d0;
    allocate(jtmp1(this % nvars, this % nvars)); jtmp1 = 0.0d0;
    allocate(jtmp2(this % nvars, this % nvars)); jtmp2 = 0.0d0;
    allocate(jtmp3(this % nvars, this % nvars)); jtmp3 = 0.0d0;

    if (this % second_order) then

       ! original function call

       call this % system % assembleResidual(tmp2, &
            & this % T(i), this % Q(i,:), &
            & this % QDOT(i,:), this % QDDOT(i,:))

       !-----------------------------------------------------------!
       ! Derivative of R WRT Q
       !-----------------------------------------------------------!

       qtmp(:) = this % Q(i,:)

       loop_vars: do k = 1, this % nvars

          ! perturb the k-th variable
          qtmp(k) = this % Q(i,k) + small

          call this % system % assembleResidual(tmp1, &
               & this % T(i), qtmp, &
               & this % QDOT(i,:), this % QDDOT(i,:))

          ! unperturb the k-th variable
          qtmp(k) = this % Q(i,k)

          ! approximate the jacobian with respect to the k-th variable
          jtmp1(:,k) = (tmp1-tmp2)/small

       end do loop_vars

       ! scale/multiply the block with the corresponding coefficient
       jtmp1 =  this % h * this % A(i,i) * this % h * this % A(i,i) * jtmp1

       !-----------------------------------------------------------!
       ! Derivative of R WRT QDOT
       !-----------------------------------------------------------!

       qdottmp(:) = this % QDOT(i,:)

       do k = 1, this % nvars

          ! perturb the k-th variable
          qdottmp(k) = this % QDOT(i,k) + small

          call this % system % assembleResidual(tmp1, &
               & this % T(i), this % Q(i,:), &
               & qdottmp, this % QDDOT(i,:))

          ! unperturb the k-th variable
          qdottmp(k) = this % QDOT(i,k)

          jtmp2(:,k) = (tmp1-tmp2)/small

       end do

       ! scale/multiply the block with the corresponding coefficient
       jtmp2 =  this % h * this % A(i,i) * jtmp2

       !-----------------------------------------------------------!
       ! Derivative of R WRT QDDOT
       !-----------------------------------------------------------!

       qddottmp(:) = this % QDDOT(i,:)

       do k = 1, this % nvars

          ! perturb the k-th variable
          qddottmp(k) = this % QDDOT(i,k) + small

          call this % system % assembleResidual(tmp1,&
               & this % T(i), this % Q(i,:), &
               & this % QDOT(i,:), qddottmp)

          ! unperturb the k-th variable
          qddottmp(k) = this % QDDOT(i,k)

          jtmp3(:,k) = (tmp1-tmp2)/small

       end do

       ! Scale/multiply the block with the corresponding coefficient
       jtmp3 =  1.0d0* jtmp3
       
       ! Add the blocks together
       jtmp = jtmp3 + jtmp2 + jtmp1 

    else

       ! original function call
       
       call this % system % assembleResidual(tmp2, &
            & this % T(i), this % Q(i,:), &
            & this % QDOT(i,:), this % QDDOT(i,:))

       !--------------------------------------------------------------!
       ! Derivative of R WRT Q
       !--------------------------------------------------------------!

       qtmp(:) = this % Q(i,:)

       loopvars: do k = 1, this % nvars

          ! perturb the k-th variable
          qtmp(k) = this % Q(i,k) + small

          call this % system % assembleResidual(tmp1, &
               & this % T(i), qtmp, &
               & this % QDOT(i,:), this % QDDOT(i,:))

          !call R(tmp1, this % nvars, this % T(i), qtmp, this % QDOT(i,:))

          ! unperturb the k-th variable
          qtmp(k) = this % Q(i,k)

          ! approximate the jacobian with respect to the k-th variable
          jtmp1(:,k) = (tmp1-tmp2)/small

       end do loopvars

       ! Scale the block with the coefficient
       jtmp1 =  this % h * this % A(i,i) * jtmp1

       !--------------------------------------------------------------!
       ! Derivative of R WRT QDOT
       !--------------------------------------------------------------!

       qdottmp(:) = this % QDOT(i,:)

       do k = 1, this % nvars

          ! perturb the k-th variable
          qdottmp(k) = this % qdot(i,k) + small

          call this % system % assembleResidual(tmp1, &
               & this % T(i), this % Q(i,:), &
               & qdottmp,this % QDDOT(i,:))

!!$          call R(tmp1, this % nvars, this % T(i), this % Q(i,:), &
!!$               & qdottmp)

          ! unperturb the k-th variable
          qdottmp(k) = this % qdot(i,k)

          jtmp2(:,k) = (tmp1-tmp2)/small

       end do

       ! sum the jacobian components to get the total derivative
       jtmp = jtmp2 + jtmp1 
       
    end if ! first or second order
    
    ! Store the jacobian for return
    if (this % approximate_jacobian) then

       exact_jac = jtmp

    else

       ! Complain about the error in Jacobian if there is any
       if (maxval(abs(exact_jac - jtmp)) .gt. small) then
          print *, "WARNING: Possible error in jacobian", &
               & maxval(abs(exact_jac - jtmp))
       end if

    end if

    deallocate(qtmp,qdottmp)
    deallocate(tmp1,tmp2)
    deallocate(jtmp1,jtmp2,jtmp)

  end subroutine check_jacobian


  !===================================================================!
  ! Exact solution to the spring mass damper system
  !===================================================================!

  function exact_solution(t, x0, v0) result (x)

    real(8) :: t, x, x0, v0
    complex(8) :: mul, a, b, term1, term2, term3

    a = 0.020d0
    b = 5.00d0

    mul = exp(-a*t*0.50d0)/sqrt(a*a - 4.00d0*b)
    term1 = a*sinh(0.50d0*t*sqrt(a*a - 4.00d0*b))
    term2 = sqrt(a*a - 4.00d0*b)*cosh(0.50d0*t*sqrt(a*a - 4.00d0*b))
    term3 = 2.00d0*v0*sinh(0.50d0*t*sqrt(a*a - 4.00d0*b))

    x = real(mul*((term1 + term2)*x0 + term3))

  end function exact_solution

end module runge_kutta_integrator
