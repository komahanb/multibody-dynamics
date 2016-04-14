!=====================================================================!
! A Diagonally Implicit Runge Kutta integrator module for first and
! second order systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!
module runge_kutta_integrator

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: RK

  !-------------------------------------------------------------------!
  ! Abstract Runge-Kutta type
  !-------------------------------------------------------------------!

  type, abstract :: RK

     integer :: num_stages = 1  ! default number of stages
     integer :: order           ! order of accuracy
     integer :: nvars = 1       ! number of states/equations

     real(8) :: h = 0.1d0       ! default step size (will reconsider when implementing adaptive step size)
     real(8) :: time            ! scalar to track integration time

     logical :: second_order = .false.

     ! The Butcher Tableau 
     real(8), dimension(:,:), allocatable :: A ! forms the coeff matrix
     real(8), dimension(:), allocatable :: B ! multiplies the state derivatives
     real(8), dimension(:), allocatable :: C ! multiplies the time

     ! The stage time and its corresponding derivatives
     real(8), dimension(:), allocatable :: T ! the corresponding stage time
     real(8), dimension(:,:), allocatable :: Q ! the corresponding state
     real(8), dimension(:,:), allocatable :: QDOT ! the stage derivatives K = F(T,Q)
     real(8), dimension(:,:), allocatable :: QDDOT ! the stage derivatives K = F(T,Q)

     real(8), dimension(:,:), allocatable :: R ! stage residual
     real(8), dimension(:,:,:,:), allocatable :: J ! stage jacobian

     ! The form of the governing equation
     logical :: descriptor_form = .true.

     ! number of function and gradient calls
     integer :: fcnt=0, fgcnt = 0

   contains

     ! Implemented common procedures (visible to the user)
     procedure :: initialize, finalize, integrate

     ! Implemented procedures (not callable by the user)
     procedure, private :: time_march
     procedure, private :: reset_stage_values
     procedure, private :: check_butcher_tableau

     ! Deferred common procedures
     procedure(compute_stage_values_interface), private, deferred :: compute_stage_values
     procedure(buthcher_interface), private, deferred :: setup_butcher_tableau

  end type RK
  
  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  
  
  type, extends(RK) :: DIRK

     ! For the nonlinear solution at each stage
     integer :: max_newton = 25
     real(8) :: tol = 1.0d-12

   contains

     private

     ! Implement/override the abstract class routines
     procedure :: setup_butcher_tableau => ButcherDIRK

     ! More specialized procedures
     procedure :: compute_stage_values

     procedure :: newton_solve
     procedure :: state_update
     procedure :: setup_linear_system

     procedure :: get_residual
     procedure :: get_jacobian

     procedure :: find_indices
     procedure :: check_jacobian

  end type DIRK
  
  !-------------------------------------------------------------------!
  ! Interfaces for deferred specialized procedures 
  !-------------------------------------------------------------------!
  
  interface

     !----------------------------------------------------------------!
     ! Interface for finding the stage derivatives at each time step
     !----------------------------------------------------------------!

     subroutine compute_stage_values_interface(this, k, q, qdot)
       import RK
       class(RK) :: this
       integer, intent(in) :: k 
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

  !-------------------------------------------------------------------!
  ! Initialize the dirk datatype and construct the tableau
  !-------------------------------------------------------------------!
  
  subroutine initialize(this, nvars, tinit, num_stages, h)

    class(RK) :: this
    integer, OPTIONAL, intent(in) :: num_stages
    integer, OPTIONAL, intent(in) :: nvars
    real(8), OPTIONAL, intent(in) :: tinit
    real(8), OPTIONAL, intent(in) :: h

    !-----------------------------------------------------------------!
    ! set the initial time
    !-----------------------------------------------------------------!

    if (present(tinit)) then
       this % time = tinit
    else
       print '("Using default start time : ",F8.3)', this % time
    end if

    !-----------------------------------------------------------------!
    ! set the order of integration
    !-----------------------------------------------------------------!

    if (present(num_stages)) then
       this % num_stages = num_stages
    else
       print '("Using default number of stages : ",i4)', this % num_stages
    end if

    !-----------------------------------------------------------------!
    ! set the user supplied initial step size
    !-----------------------------------------------------------------!

    if (present(h)) then
       this % h = h 
    else
       print '("Using default step size h : ", E9.3)', this % h
    end if

    !-----------------------------------------------------------------!
    ! set the user supplied number of variables
    !-----------------------------------------------------------------!

    if (present(nvars)) then
       this % nvars = nvars 
    else
       print '("Using default nvars : ", i4)', this % nvars
    end if

    !-----------------------------------------------------------------!
    ! allocate space for the tableau
    !-----------------------------------------------------------------!

    allocate(this % A(this % num_stages, this % num_stages))
    this % A = 0.0d0

    allocate(this % B(this % num_stages))    
    this % B = 0.0d0

    allocate(this % C(this % num_stages))
    this % C = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the stage state
    !-----------------------------------------------------------------!

    allocate(this % Q(this % num_stages, this % nvars))
    this % Q = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the stage derivatives
    !-----------------------------------------------------------------!

    allocate(this % QDOT(this % num_stages, this % nvars))
    this % QDOT = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the second stage derivatives
    !-----------------------------------------------------------------!

    allocate(this % QDDOT(this % num_stages, this % nvars))
    this % QDDOT = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the stage time
    !-----------------------------------------------------------------!

    allocate(this % T(this % num_stages))
    this % T = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the stage time
    !-----------------------------------------------------------------!

    allocate(this % R(this % num_stages, this % nvars))
    this % R = 0.0d0

    !-----------------------------------------------------------------!
    ! allocate space for the stage time
    !-----------------------------------------------------------------!

    allocate(this % J(this % num_stages,&
         & this % num_stages, this % nvars, this % nvars))
    this % J = 0.0d0

    !-----------------------------------------------------------------!
    ! this subroutine puts values into the Butcher tableau
    !-----------------------------------------------------------------!

    call this % setup_butcher_tableau()

    !-----------------------------------------------------------------!
    ! sanity check
    !-----------------------------------------------------------------!

    call this % check_butcher_tableau()

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

    ! clear butcher's tableau
    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

    ! clear stage value
    if(allocated(this % QDDOT)) deallocate(this % QDDOT)
    if(allocated(this % QDOT)) deallocate(this % QDOT)
    if(allocated(this % Q)) deallocate(this % Q)
    if(allocated(this % T)) deallocate(this % T)

    ! clear the stage residual and jacobian
    if(allocated(this % R)) deallocate(this % R)
    if(allocated(this % J)) deallocate(this % J)

  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Time integration logic
  !-------------------------------------------------------------------!
  ! Input: 
  ! o state arrays q and qdot with initial conditions set at q(1)
  ! o number of steps N
  ! o step size h
  !-------------------------------------------------------------------!
  ! Output:
  ! o q, qdot arrays are modified by the routine
  !-------------------------------------------------------------------!
  
  subroutine Integrate(this, N, q, qdot)

    class(RK) :: this
    real(8), intent(inout), dimension(:,:) :: q, qdot
    integer, intent(in) :: N 
    integer :: k

    ! March in time
    march: do k = 2, N + 1

       ! find the stage derivatives at the current step
       call this % compute_stage_values(k, q, qdot)

       ! advance the state to the current step
       call this % time_march(k, q, qdot)

       ! set the stage values to zero
       call this % reset_stage_values()

    end do march

  end subroutine Integrate
  
  !-------------------------------------------------------------------!
  ! Update the states based on RK Formulae
  !-------------------------------------------------------------------!
  
  subroutine time_march(this, k, q, qdot)

    implicit none

    class(RK) :: this
    integer, intent(in) :: k ! current time step
    real(8),  dimension(:,:) :: q, qdot ! current state
    integer :: m

    ! march q to next time step
    forall(m=1:this%nvars)
       q(k,m) = q(k-1,m) + this % h*sum(this % B(1:this%num_stages) &
            &* this % QDOT(1:this%num_stages,m))
    end forall

    ! march qdot to next time step for second order system
    if (this % second_order) then
       forall(m=1:this%nvars)
          qdot(k,m) = qdot(k-1,m) + this % h*sum(this % B(1:this%num_stages) &
               &* this % QDDOT(1:this%num_stages,m))
       end forall
    end if

    ! increment the time
    this % time = this % time + this % h

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
  
  subroutine compute_stage_values(this, k, q, qdot)

    class(DIRK) :: this
    integer, intent(in) :: k 
    real(8), intent(in), dimension(:,:) :: q
    real(8), OPTIONAL, intent(in), dimension(:,:) :: qdot
    integer :: j

    do j = 1, this % num_stages

       ! Find the stage times

       this % T(j) = this % time + this % C(j)*this % h

       ! Guess the solution for stage states

       if (.not. this % descriptor_form) then
          ! guess q
          this % Q(j,:) = 1.0d0
       else 
          if (this % second_order) then
             ! guess qddot
             this % QDDOT(j,:) = 1.0d0 
          else
             ! guess qdot
             this % QDOT(j,:) = 1.0d0 
          end if
       end if

       ! solve the non linear stage equations using Newton's method for
       ! the actual stage states 

       call this % newton_solve(q(k-1,:), qdot(k-1,:))

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
  
  subroutine newton_solve(this, qk, qdotk)

    class(DIRK) :: this
    real(8), intent(in), dimension(:)    :: qk, qdotk
    real(8), allocatable, dimension(:)   :: res, dq
    real(8), allocatable, dimension(:,:) :: jac
    integer, allocatable, dimension(:)   :: ipiv
    integer :: n, info, size
    logical :: conv = .false.

    ! find the size of the linear system based on the calling object
    size = this % nvars

    if (.not.allocated(ipiv)) allocate(ipiv(size))
    if (.not.allocated(res)) allocate(res(size))
    if (.not.allocated(dq)) allocate(dq(size))
    if (.not.allocated(jac)) allocate(jac(size,size))

    newton: do n = 1, this % max_newton

       ! Get the residual of the function
       call this % get_residual(qk, qdotk)
       this % fcnt = this % fcnt + 1

       ! Get the jacobian matrix
       call this % get_jacobian()
       this % fgcnt = this % fgcnt + 1

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
       print '("Newton solve failed : time = ", E10.3, " iters = ", i3,&
            & " |R| = ",E10.3," |dq| = ",E10.3)',&
            & this % time, n, norm2(res), norm2(dq)
       stop
    else
       print '("Newton solve: time = ", E10.3, " iters = ", i3,&
            & " |R| = ",E10.3," |dq| = ",E10.3)',&
            & this % time, n, norm2(res), norm2(dq)
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

    integer :: stage_num

    call this % find_indices(stage_num, stage_num)

    res = this % R(stage_num,:)

    jac = this % J(stage_num, stage_num, :, :)

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

    call this % find_indices(i, i)

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
  
  !-----------------------------------------------------------------!    
  ! Select type and set appropriate indices for looping  
  !-----------------------------------------------------------------!
  
  subroutine find_indices(this, istart, iend)

    class(DIRK) :: this
    integer, intent(inout) :: istart, iend
    logical :: found = .false.
    integer :: i

    findstagenum: do i = this % num_stages, 1, -1

       ! we hope to find the last non-zero time state
       if (this % T(i) .ne. 0.0d0) then
          istart = i
          iend = i
          found = .true.
          exit findstagenum
       end if

       if (.not. found ) stop "index finding failed!"

    end do findstagenum

  end subroutine find_indices
  
  !-------------------------------------------------------------------!
  ! Computes the stage residual for the set stage state Y (comes from
  ! Newton's iteration) and sets into the same instance
  !
  ! R_{i}= q_{k,i} - q_{k} - h \sum_{j=1}^s {a_{ij} f(t_{k,j}, q_{k,j})
  ! i = 1,\ldots,s 
  !-------------------------------------------------------------------!
  
  subroutine get_residual(this, qk, qdotk)

    class(DIRK) :: this
    real(8), intent(in), dimension(:) :: qk, qdotk
    integer :: i, m
    integer :: istart, iend

    external :: R

    ! get the appropriate indices based on type and stage number
    call this % find_indices(istart, iend)

    if (this % second_order) then

       ! compute the stage velocities for the guessed QDDOT
       do i = istart, iend
          forall(m = 1 : this % nvars)
             this % QDOT(i,m) = qdotk(m) &
                  & + this % h*sum(this % A(i,:)*this % QDDOT(:, m))
          end forall
       end do

       ! compute the stage states for the guessed QDDOT
       do i = istart, iend
          forall(m = 1 : this % nvars)
             this % Q(i,m) = qk(m) &
                  & + this % h*sum(this % A(i,:)*this % QDOT(:, m))
          end forall
       end do

       ! compute the stage residuals for Q, QDOT, QDDOT
       do i = istart, iend
          call R(this % R(i,:), this % nvars, this % T(i), this % Q(i,:), &
               & this % QDOT(i,:), this % QDDOT(i,:))
       end do

    else 

       ! compute the stage states for the guessed QDOT
       do i = istart, iend
          forall(m = 1 : this % nvars)
             this % Q(i,m) = qk(m) &
                  & + this % h*sum(this % A(i,:)*this % QDOT(:, m))
          end forall
       end do

       ! compute the stage residuals
       do i = istart, iend
          call R(this % R(i,:), this % nvars, this % T(i), this % Q(i,:), &
               & this % QDOT(i,:))
       end do

    end if

  end subroutine get_residual

  !-------------------------------------------------------------------!
  ! Computes the stage jacobian and sets into the same instance
  !          J(i,j) = [ 1 - h A(i,j) DFDQ(T(j),Y(j))]
  !-------------------------------------------------------------------!
  
  subroutine get_jacobian(this)

    class(DIRK) :: this
    integer :: i, j
    real(8) :: alpha
    external :: DRDQ, DRDQDOT, DRDQDDOT

    ! get the appropriate indices based on type and stage number
    call this % find_indices(j, i)

    this % J(i,j,:,:) = 0.0d0


    if (this % second_order) then

       ! get the q block
       alpha =  this % h * this % A(i,i)* this % h * this % A(i,i)
       call DRDQ(this % J(i,j,:,:), alpha, this % nvars, this % T(j), &
            & this % Q(j,:), this % QDOT(j,:), this % QDDOT(j,:))

       ! get the qdot block
       alpha =  this % h * this % A(i,i)
       call DRDQDOT(this % J(i,j,:,:), alpha, this % nvars, this % T(j), &
            & this % Q(j,:), this % QDOT(j,:), this % QDDOT(j,:))

       ! get the qddot block
       alpha = 1.0d0
       call DRDQDDOT(this % J(i,j,:,:), alpha, this % nvars, this % T(j), &
            & this % Q(j,:), this % QDOT(j,:), this % QDDOT(j,:))

    else

       ! get the q block
       alpha = this % h * this % A(i,i)
       call DRDQ(this % J(i,j,:,:), alpha, this % nvars, this % T(j), &
            & this % Q(j,:), this % QDOT(j,:))

       ! get the qdot block
       alpha = 1.0d0
       call DRDQDOT(this % J(i,j,:,:), alpha, this % nvars, this % T(j), &
            & this % Q(j,:), this % QDOT(j,:))

    end if

    ! check with FD
    call this % check_jacobian(i, this % J(i,j,:,:))

  end subroutine get_jacobian

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

       call R(tmp2, this % nvars, this % T(i), this % Q(i,:), &
            & this % QDOT(i,:), this % QDDOT(i,:))

       !-----------------------------------------------------------!
       ! Derivative of R WRT Q
       !-----------------------------------------------------------!

       qtmp(:) = this % Q(i,:)

       loop_vars: do k = 1, this % nvars

          ! perturb the k-th variable
          qtmp(k) = this % Q(i,k) + small

          call R(tmp1, this % nvars, this % T(i), qtmp, &
               & this % QDOT(i,:), this % QDDOT(i,:))

          ! unperturb the k-th variable
          qtmp(k) = this % Q(i,k)

          ! approximate the jacobian with respect to the k-th variable
          jtmp1(:,k) = (tmp1-tmp2)/small

       end do loop_vars

       jtmp1 =  this % h * this % A(i,i) * this % h * this % A(i,i) * jtmp1

       !-----------------------------------------------------------!
       ! Derivative of R WRT QDOT
       !-----------------------------------------------------------!

       qdottmp(:) = this % QDOT(i,:)

       do k = 1, this % nvars

          ! perturb the k-th variable
          qdottmp(k) = this % QDOT(i,k) + small

          call R(tmp1, this % nvars, this % T(i), this % Q(i,:), &
               & qdottmp, this % QDDOT(i,:))

          ! unperturb the k-th variable
          qdottmp(k) = this % QDOT(i,k)

          jtmp2(:,k) = (tmp1-tmp2)/small

       end do

       jtmp2 =  this % h * this % A(i,i) * jtmp2

       !-----------------------------------------------------------!
       ! Derivative of R WRT QDDOT
       !-----------------------------------------------------------!

       qddottmp(:) = this % QDDOT(i,:)

       do k = 1, this % nvars

          ! perturb the k-th variable
          qddottmp(k) = this % QDDOT(i,k) + small

          call R(tmp1, this % nvars, this % T(i), this % Q(i,:), &
               & this % QDOT(i,:), qddottmp)

          ! unperturb the k-th variable
          qddottmp(k) = this % QDDOT(i,k)

          jtmp3(:,k) = (tmp1-tmp2)/small

       end do

       jtmp = jtmp3 + jtmp2 + jtmp1 

    else

       ! original function call

       call R(tmp2, this % nvars, this % T(i), this % Q(i,:), &
            & this % QDOT(i,:))

       !--------------------------------------------------------------!
       ! Derivative of R WRT Q
       !--------------------------------------------------------------!

       qtmp(:) = this % Q(i,:)

       loopvars: do k = 1, this % nvars

          ! perturb the k-th variable
          qtmp(k) = this % Q(i,k) + small

          call R(tmp1, this % nvars, this % T(i), qtmp, this % QDOT(i,:))

          ! unperturb the k-th variable
          qtmp(k) = this % Q(i,k)

          ! approximate the jacobian with respect to the k-th variable
          jtmp1(:,k) = (tmp1-tmp2)/small

       end do loopvars

       jtmp1 =  this % h * this % A(i,i) * jtmp1

       !-----------------------------------------------------------------!
       ! Derivative of R WRT QDOT
       !-----------------------------------------------------------------!

       qdottmp(:) = this % QDOT(i,:)

       do k = 1, this % nvars

          ! perturb the k-th variable
          qdottmp(k) = this % Qdot(i,k) + small

          call R(tmp1, this % nvars, this % T(i), this % Q(i,:), &
               & qdottmp)

          ! unperturb the k-th variable
          qdottmp(k) = this % Qdot(i,k)

          jtmp2(:,k) = (tmp1-tmp2)/small

       end do

       ! sum the jacobian components to get the total derivative
       jtmp = jtmp2 + jtmp1 

    end if ! first or second order

    if (maxval(abs(exact_jac - jtmp)) .gt. small) then
       print *, "WARNING: Possible error in jacobian", &
            & maxval(abs(exact_jac - jtmp))
    end if

!!$    exact_jac = jtmp

    deallocate(qtmp,qdottmp)
    deallocate(tmp1,tmp2)
    deallocate(jtmp1,jtmp2,jtmp)

  end subroutine check_jacobian

end module runge_kutta_integrator
