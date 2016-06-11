!=====================================================================!
! A Diagonally Implicit Runge Kutta integrator module for first and
! second order systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module runge_kutta_integrator

  use iso_fortran_env , only : dp => real64
  use integrator_class, only : integrator
  use physics_class,    only : physics

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
     ! The lagrange multipliers
     !----------------------------------------------------------------!
     
     real(dp), dimension(:,:,:), allocatable   :: psi

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
     procedure, public  :: marchBackwards
     procedure, private :: assembleRHS
     procedure, private :: computeTotalDerivative

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
       use iso_fortran_env , only : dp => real64
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
       use iso_fortran_env , only : dp => real64
       import RK
       class(RK) :: this
     end subroutine buthcher_interface

  end interface

contains

  !===================================================================!
  ! Initialize the dirk datatype and construct the tableau
  !===================================================================!
  
  subroutine initialize(this, system, tinit, tfinal, h, second_order, num_stages)
    
    class(RK)                       :: this
    class(physics), target :: system
    integer  , OPTIONAL, intent(in) :: num_stages
    real(dp) , OPTIONAL, intent(in) :: tinit, tfinal
    real(dp) , OPTIONAL, intent(in) :: h
    logical  , OPTIONAL, intent(in) :: second_order

    print *, "======================================"
    print *, ">> Diagonally-Implicit-Runge-Kutta  <<"
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
    ! Allocate space for the lagrange multipliers
    !-----------------------------------------------------------------!
    
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

    ! Clear global states and time
    if(allocated(this % UDDOT)) deallocate(this % UDDOT)
    if(allocated(this % UDOT)) deallocate(this % UDOT)
    if(allocated(this % U)) deallocate(this % U)
    if(allocated(this % time)) deallocate(this % time)

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
  
  subroutine marchBackwards(this)

    class(DIRK) :: this
    integer     :: k, i
    real(dp)    :: alpha, beta, gamma

    do k = this % num_steps, 2, -1
       
       this % current_step = k
       
       do i = this % num_stages, 1, -1

          this % current_stage = i
          
          !--------------------------------------------------------------!
          ! Determine the linearization coefficients for the Jacobian
          !--------------------------------------------------------------!
          
          if (this % second_order) then
             gamma = 1.0d0
             beta  = this % h * this%A(i,i)
             alpha = this % h * this%A(i,i)* this % h * this%A(i,i)
          else
             gamma = 0.0d0
             beta  = 1.0d0
             alpha = this % h * this%A(i,i)
          end if
          
          !--------------------------------------------------------------!
          ! Solve the adjoint equation at each step
          !--------------------------------------------------------------!

          call this % adjointSolve(this % psi(k,i,:), alpha, beta, gamma,&
               & this % t(i), this % q(k,i,:), this % qdot(k,i,:), this % qddot(k,i,:))
          
       end do

    end do

!!$
!!$    ! Compute the total derivative
!!$    tmp = 0.0d0
!!$    do k = 1, this % num_steps
!!$       tmp = tmp + matmul(this%psi(k,i,:), dRdx)
!!$    end do
!!$
!!$    ! call addDVSens
!!$    dLdx = dfdx + tmp
!!$
!!$    ! Write the adjoint variables 
!!$    open(unit=90, file='output/adjoint.dat')
!!$    do k = 1, this % num_steps
!!$       write(90, *)  this % time(k), this % psi(k,1,:), this % psi(k,2,:), this % psi(k,3,:)
!!$    end do
!!$    close(90)
    
  end subroutine marchBackwards

  !===================================================================!
  ! Compute the total derivative of the function with respect to the
  ! design variables and return the gradient 
  !===================================================================!
  
  subroutine computeTotalDerivative( this, dfdx )
    
    class(DIRK)                              :: this
    real(dp) , dimension(:), intent(inout)   :: dfdx
    real(dp) , dimension(:,:), allocatable   :: dRdX
    integer                                  :: k
    
!!$    allocate(dRdX(this % nsvars, this % ndvars))
!!$    dfdx = 0.0d0
!!$    
!!$    !-----------------------------------------------------------------!
!!$    ! Compute dfdx
!!$    !-----------------------------------------------------------------!
!!$
!!$    do k = 2, this % num_steps
!!$       call this % system % func % addDfdx(dfdx, 1.0d0, this % time(k), &
!!$            & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:) )
!!$    end do
!!$
!!$    !-----------------------------------------------------------------!
!!$    ! Compute the total derivative
!!$    !-----------------------------------------------------------------!
!!$
!!$    do k = 2, this % num_steps
!!$       call this % system % getResidualDVSens(dRdX, 1.0d0, this % time(k), &
!!$            & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:))
!!$       dfdx = dfdx + matmul(this % psi(k,:), dRdX) ! check order
!!$    end do
!!$    
!!$    !-----------------------------------------------------------------!
!!$    ! Special logic for initial condition (use the adjoint variable
!!$    ! for the second time step)
!!$    !-----------------------------------------------------------------!
!!$    
!!$    call this % system % func % addDfdx(dfdx, 1.0d0, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(2,:) )
!!$    
!!$    call this % system % getResidualDVSens(dRdX, 1.0d0, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(2,:))
!!$    dfdx = dfdx + matmul(this % psi(2,:), dRdX)
!!$    
!!$    deallocate(dRdX)

  end subroutine computeTotalDerivative

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

  subroutine computeStageStateValues( this, q, qdot )

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
          if (k .eq. 2) then ! initialize with a starting value
             this % QDDOT(k,j,:) = 1.0d0
          else 
             if (j .eq. 1) then ! copy previous global state
                this % QDDOT(k,j,:) = this % UDDOT(k-1,:)
             else ! copy previous local state
                this % QDDOT(k,j,:) = this % QDDOT(k,j-1,:)
             end if
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
          if (k .eq. 2) then ! initialize with a starting value
             this % QDOT(k,j,:) = 1.0d0
          else 
             if (j .eq. 1) then ! copy previous global state
                this % QDOT(k,j,:) = this % UDOT(k-1,:)
             else ! copy previous local state
                this % QDOT(k,j,:) = this % QDOT(k,j-1,:)
             end if
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
 
    class(DIRK)                           :: this
    real(dp), dimension(:), intent(inout) :: rhs
   
    ! Add the contributions from the objective function
    call this % AddFunctionDependency(rhs)
    call this % AddTimeDependency(rhs)
    call this % AddStageDependency(rhs)

  end subroutine assembleRHS

end module runge_kutta_integrator
