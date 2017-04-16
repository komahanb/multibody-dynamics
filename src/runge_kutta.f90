#include "scalar.fpp"
!=====================================================================!
! A Diagonally Implicit Runge Kutta integrator module for first and
! second order systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module runge_kutta_integrator

  use integrator_class, only : integrator
  use physics_class,    only : physics
  use linear_algebra,   only : solve
  implicit none

  private
  public :: DIRK

  ! Useful constants
  type(scalar), parameter :: ONE  = 1.0d0
  type(scalar), parameter :: ZERO = 0.0d0

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

     type(scalar), dimension(:,:), allocatable :: A ! forms the coeff matrix
     type(scalar), dimension(:), allocatable   :: B ! multiplies the state derivatives
     type(scalar), dimension(:), allocatable   :: C ! multiplies the time

     !----------------------------------------------------------------!
     ! The stage time and its corresponding derivatives
     !----------------------------------------------------------------!

     real(dp), dimension(:), allocatable :: T
     type(scalar), dimension(:,:,:), allocatable :: Q
     type(scalar), dimension(:,:,:), allocatable :: QDOT
     type(scalar), dimension(:,:,:), allocatable :: QDDOT

     !----------------------------------------------------------------!
     ! The lagrange multipliers
     !----------------------------------------------------------------!
     
     type(scalar), dimension(:,:,:), allocatable :: LAM

   contains

     ! Routines for integration
     
     procedure          :: integrate => integrate
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff

     procedure, private :: TimeMarch
     procedure, private :: CheckButcherTableau
     
     procedure(buthcher_interface), private, deferred :: SetupButcherTableau
     
     procedure :: testAdjoint6

  end type RK

  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(RK) :: DIRK

     type(scalar), dimension(:,:,:), allocatable :: rhsbin
     type(scalar), dimension(:,:), allocatable :: psibin
     type(scalar), dimension(:,:), allocatable :: phibin

   contains

     ! Destructor
     procedure :: finalize

     procedure :: setupButcherTableau => ButcherDIRK

     !----------------------------------------------------------------!
     ! Adjoint procedures
     !----------------------------------------------------------------!
     procedure, public  :: marchBackwards
     procedure, private :: assembleRHS
     procedure, public :: computeTotalDerivative
     procedure, public :: evalFunc => evalFuncDIRK

     procedure, private :: evaluate_adjoint
     procedure, private :: distribute_contributions
    
  end type DIRK

  !===================================================================!
  ! Interfaces for deferred specialized procedures 
  !===================================================================!

  interface

     !================================================================!
     ! Interface for setting the Butcher tableau for each type of RK
     ! scheme
     !================================================================!

     subroutine buthcher_interface(this)
       import RK
       class(RK) :: this
     end subroutine buthcher_interface

  end interface
  
  interface DIRK
     module procedure initialize
  end interface DIRK

contains

  !===================================================================!
  ! Initialize the dirk datatype and construct the tableau
  !===================================================================!
  
  type(DIRK) function initialize( system, tinit, tfinal, h, second_order, num_stages ) result(this)
    
    class(physics), target :: system
    integer  , OPTIONAL, intent(in) :: num_stages
    real(dp) , OPTIONAL, intent(in) :: tinit, tfinal
    real(dp), OPTIONAL, intent(in) :: h
    logical  , OPTIONAL, intent(in) :: second_order

    print *, "======================================"
    print *, ">> Diagonally-Implicit-Runge-Kutta  <<"
    print *, "======================================"

    call this % construct(system, tinit, tfinal, h, second_order)
    
    !-----------------------------------------------------------------!
    ! Set the order/stages of integration
    !-----------------------------------------------------------------!

    if (present(num_stages)) then
       this % num_stages = num_stages
    end if
    print '("  >> Number of stages       : ",i4)', this % num_stages
    
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
    
    allocate(this % lam(this % num_steps, this % num_stages, this % nsvars))
    this % lam = 0.0d0

    !-----------------------------------------------------------------!
    ! Allocate space for the tableau
    !-----------------------------------------------------------------!

    allocate(this % A(this % num_stages, this % num_stages))
    this % A = 0.0d0

    allocate(this % B(this % num_stages))    
    this % B = 0.0d0

    allocate(this % C(this % num_stages))
    this % C = 0.0d0

    call this % setupButcherTableau()
    call this % checkButcherTableau()

    !-----------------------------------------------------------------!
    ! Setup adjoint RHS
    !-----------------------------------------------------------------!
    
    this % num_rhs_bins = this % num_stages

    allocate(this % rhs(this % nsvars))
    this % rhs = 0.0d0
    
    !-----------------------------------------------------------------!
    ! Setup adjoint RHS
    !-----------------------------------------------------------------!   

    this % num_rhs_bins = this % num_stages

    allocate(this % rhsbin(this % num_steps, this % num_stages, this % nsvars))
    this % rhsbin = 0.0d0

    allocate(this % phibin(this % num_steps, this % nsvars))
    this % phibin = 0.0d0

    allocate(this % psibin(this % num_steps, this % nsvars))
    this % psibin = 0.0d0

  end function initialize

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

    if (abs(sum(this % B) - 1.0d0) .gt. 5.0d-16) then
       print *, "WARNING: sum(B) != 1", this % num_stages
    end if

  end subroutine checkButcherTableau
  
  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!
  
  subroutine finalize(this)

    class(DIRK) :: this

    ! Clear butcher's tableau
    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

    ! Clear stage state values
    if(allocated(this % QDDOT)) deallocate(this % QDDOT)
    if(allocated(this % QDOT)) deallocate(this % QDOT)
    if(allocated(this % Q)) deallocate(this % Q)
    if(allocated(this % T)) deallocate(this % T)

    ! Clear stage lagrange multipliers
    if(allocated(this % lam)) deallocate(this % lam)

    ! Adjoint variables and RHS
    if(allocated(this % rhs)) deallocate(this % rhs)
    if(allocated(this % rhsbin)) deallocate(this % rhsbin)
    if(allocated(this % psibin)) deallocate(this % psibin)
    if(allocated(this % phibin)) deallocate(this % phibin)

  end subroutine finalize
  
  !===================================================================!
  ! Overridden Time integration logic
  !===================================================================!
  
  subroutine Integrate(this)

    use nonlinear_algebra, only: nonlinear_solve

    class(RK) :: this
    type(scalar)  :: alpha, beta, gamma
    integer :: k, i

    ! Set states to zero

    this % U     = 0.0d0
    this % UDOT  = 0.0d0
    this % UDDOT = 0.0d0
    this % time  = 0.0d0

    this % Q     = 0.0d0
    this % QDOT  = 0.0d0
    this % QDDOT = 0.0d0
    this % T     = 0.0d0

    call this % system % getInitialStates(this % time(1), this % u(1,:), this % udot(1,:))

    !-----------------------------------------------------------------!
    ! March in time
    !-----------------------------------------------------------------!
    
    this % current_step = 1
    
    time: do k = 2, this % num_steps
       
       this % current_step = k

       !-----------------------------------------------------------------!
       ! Loop over stages
       !-----------------------------------------------------------------!
       
       do i = 1, this % num_stages

          this % current_stage = i

          ! Find the stage times
          this % T(i) = this % time(k-1) + this % C(i)*this % h

          ! Guess the solution for stage states
          call this % approximateStates()

          call this % getLinearCoeff(alpha, beta, gamma)

          ! solve the non linear stage equations using Newton's method
          call nonlinear_solve(this % system, &
               & alpha, beta, gamma, &
               & this % time(k), this % q(k,i,:), this % qdot(k,i,:), this % qddot(k,i,:))

       end do

       ! Advance the state to the current step
       call this % timeMarch(this % u, this % udot, this % uddot)

    end do time

  end subroutine Integrate


  !===================================================================!
  ! Time backwards in stage and backwards in time to solve for the
  ! adjoint variables
  !===================================================================!
  
  subroutine marchBackwards( this )

    class(DIRK)    :: this
    type(integer) :: k,i
    
    time: do k = this % num_steps, 2, -1

       this % current_step = k
       
       stage: do i = this % num_stages, 1, -1
          
          this % current_stage = i

          !--------------------------------------------------------------!
          ! Solve the adjoint equation at each step
          !--------------------------------------------------------------!

          ! print *, k, i, this % lam(k,i,:), this % psi(k,:), this % phi(k,:)
          call this % evaluate_adjoint(this % lam(k,i,:), this % psi(k,:), this % phi(k,:))
          !print *,  k, i, this % lam(k,i,:), this % psi(k,:), this % phi(k,:)

          !--------------------------------------------------------------!
          ! Drop the contributions from this step to corresponding bins
          !--------------------------------------------------------------!

          call this % distribute_contributions(this % rhsbin, this % psibin, this % phibin)

       end do stage
       
    end do time
    
  end subroutine marchBackwards

  !===================================================================!
  ! Evaluates the adjoint variable values at the current step
  !===================================================================!
  
  subroutine evaluate_adjoint(this, mu, psi, phi)

    class(DIRK)                  :: this
    type(integer)                :: k, i
    type(scalar)                 :: alpha, beta, gamma
    type(scalar), dimension(:) :: mu, psi, phi

    ! Retrive the current step number
    k = this % current_step
    i = this % current_stage

    ! Evaluate PHI (no linear solution necessary)
    phi = this % phibin(k,:)

    ! Evaluate PSI (no linear solution necessary)
    psi = this % psibin(k,:)

    ! Evaluate MU
    call this % getLinearCoeff(alpha, beta, gamma)

    call this % adjointSolve(mu, &
         & alpha, beta, gamma, &
         & this % T(i), &
         & this % q(k,i,:), this % qdot(k,i,:), this % qddot(k,i,:))

  end subroutine evaluate_adjoint
  
  !===================================================================!
  ! Add contributions from kth step to k- steps
  !===================================================================!
  
  subroutine distribute_contributions(this, rhsbin, psibin, phibin)

    class(DIRK)                  :: this
    type(integer)                :: k, s, i, j, p
    type(scalar)                 :: alpha, beta, gamma
    type(scalar), dimension(:,:,:) :: rhsbin
    type(scalar), dimension(:,:) :: phibin, psibin

    ! Retrive the current step number
    k = this % current_step
    i = this % current_stage
    s = this % num_stages

    !-----------------------------------------------------------------!
    ! Add contributions from k to k- PHIBIN
    !-----------------------------------------------------------------!

    phibin(k-1,:) = phibin(k-1,:) + this % phi(k,:)

    gamma = 0.0d0
    beta  = 0.0d0
    alpha = this%h*this%B(i)

    call this % addFuncResAdjPdt(phibin(k-1,:), &
         & alpha, beta, gamma, this % T(i), &
         & this % Q(k,i,:), this % Qdot(k,i,:), this % Qddot(k,i,:), &
         & this % lam(k,i,:))

    !-----------------------------------------------------------------!
    ! Add contributions from k to k- PSIBIN    
    !-----------------------------------------------------------------!

    psibin(k-1,:) = psibin(k-1,:) + this % psi(k,:)

    gamma = 0.0d0
    beta  = this % h * this % B(i)
    alpha = this % h * this % B(i) * this % h * sum(this % A(i,1:i))
    
    call this % addFuncResAdjPdt(psibin(k-1,:), &
         & alpha, beta, gamma, this % T(i), &
         & this % Q(k,i,:), this % Qdot(k,i,:), this % Qddot(k,i,:), &
         & this % lam(k,i,:))

    psibin(k-1,:) = psibin(k-1,:) + beta*this%phi(k,:)

  end subroutine distribute_contributions

  !===================================================================!
  ! Update the states based on RK Formulae
  !===================================================================!
  
  subroutine timeMarch(this, q, qdot, qddot)

    implicit none

    class(RK) :: this
    type(scalar),  dimension(:,:) :: q, qdot, qddot ! current state
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
    type(scalar), parameter :: PI = 22.0d0/7.0d0
    type(scalar), parameter :: tmp  = 1.0d0/(2.0d0*dsqrt(3.0d0))
    type(scalar), parameter :: half = 1.0d0/2.0d0
    type(scalar), parameter :: one  = 1.0d0
    type(scalar), parameter :: alpha = 2.0d0*cos(PI/18.0d0)/dsqrt(3.0d0)

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
  ! Returns the linearization scalar coefficients: alpha, beta, gamma
  !===================================================================!
  
  subroutine getLinearCoeff( this, alpha, beta, gamma )
    
    class(RK) :: this
    type(scalar), intent(out)   :: alpha, beta, gamma
    integer :: k, i

    k = this % current_step
    i = this % current_stage

    if (this % second_order) then
       gamma = 1.0d0
       beta  = this % h * this%A(i,i)
       alpha = this % h * this%A(i,i)* this % h * this%A(i,i)
    else
       gamma = 0.0d0
       beta  = 1.0d0
       alpha = this % h * this % A(i,i)
    end if

  end subroutine getLinearCoeff

  !===================================================================!
  ! Approximate the state variables at each step using BDF formulae
  !===================================================================!
  
  subroutine approximateStates( this )

    class(RK) :: this
    integer    :: k, m, i

    k = this % current_step
    i = this % current_stage
    
    if (this % second_order) then

       ! guess qddot
       if (i .eq. 1) then ! copy previous global state
          this % QDDOT(k,i,:) = this % UDDOT(k-1,:)
       else ! copy previous local state
          this % QDDOT(k,i,:) = this % QDDOT(k,i-1,:)
       end if

       ! compute the stage velocity states for the guessed QDDOT
       forall(m = 1 : this % nsvars)
          this % QDOT(k,i,m) = this % udot(k-1,m) &
               & + this % h*sum(this % A(i,:)&
               & * this % QDDOT(k,:, m))
       end forall

       ! compute the stage states for the guessed QDDOT
       forall(m = 1 : this % nsvars)
          this % Q(k,i,m) = this % u(k-1,m) &
               & + this % h*sum(this % A(i,:)*this % QDOT(k,:, m))
       end forall

    else

       ! guess qdot
       if (i .eq. 1) then ! copy previous global state
          this % QDOT(k,i,:) = this % UDOT(k-1,:)
       else ! copy previous local state
          this % QDOT(k,i,:) = this % QDOT(k,i-1,:)
       end if

       ! compute the stage states for the guessed 
       forall(m = 1 : this % nsvars)
          this % Q(k,i,m) = this % u(k-1,m) &
               & + this % h*sum(this % A(i,:)*this % QDOT(k,:, m))
       end forall

    end if

  end subroutine approximateStates
  
  !===================================================================!
  ! Function that puts together the right hand side of the adjoint
  ! equation into the supplied rhs vector.
  !===================================================================!

  subroutine assembleRHS(this, rhs)
 
    class(DIRK)                                :: this
    type(scalar), dimension(:), intent(inout) :: rhs
    type(scalar)                              :: scale
    integer                                   :: k, j, i, m, s, p
    type(scalar)                              :: alpha, beta, gamma

    k = this % current_step
    i = this % current_stage
    s = this % num_stages

    ! Replace from rhsbin for this step (contains k+ terms)
    ! rhs = this % rhsbin(k,i,:)
    !print *, k , i, rhs

    ! Calculate and add terms from k-th step
    call this % getLinearCoeff(alpha, beta, gamma)

    alpha = alpha* this % B(i)
    beta  = beta* this % B(i)
    gamma = gamma* this % B(i)

    ! Add the state variable sensitivity from the previous step
    call this % system % func % addFuncSVSens(rhs, &
         & alpha, beta, gamma, &
         & this % T(i), &
         & this % system % X, &
         & this % Q(k,i,:), &
         & this % Qdot(k,i,:), &
         & this % Qddot(k,i,:))

    ! Add contributions from psi
    rhs = rhs + this% B(i)*this % psi(k,:)

    ! Add contributions from PHI
    scale = 0.0d0
    do j = i, s
       scale = scale + this % h * this % B(j) * this % A(j,i)
    end do
    rhs = rhs + scale * this % phi(k,:)

    !-----------------------------------------------------------------!
    ! Add contributions from k to k- RHSBIN    
    !-----------------------------------------------------------------!
    
    do j = i+1, s
!!$       gamma = 0.0d0
!!$       beta  = this % B(j) * this % h * this % A(j,i)
!!$       alpha = 0.0d0
!!$       do p = i, j
!!$          alpha = alpha + this % A(j,p) * this % A(p,i)
!!$       end do
!!$       alpha = alpha * this % h**2 * this % (j)

       gamma = 0.0d0
       beta  = this % B(j) * this % h * this % A(j,i)
       do p = i, j
          alpha = this % B(j) * this % h * this % h * this % A(j,p) * this % A(p,i)
       end do
       
       !! Need to fix this
       call this % addFuncResAdjPdt(rhs, &
            & alpha, beta, gamma, this % T(j), &
            & this % Q(k,j,:), this % Qdot(k,j,:), this % qddot(k,j,:), &
            & this % lam(k,j,:))

    end do

    ! Negate the RHS
    rhs = -rhs

  end subroutine assembleRHS

  !===================================================================!
  ! Compute the total derivative of the function with respect to the
  ! design variables and return the gradient 
  !===================================================================!
  
  subroutine computeTotalDerivative( this, dfdx )

    class(DIRK)                                  :: this
    type(scalar) , dimension(:), intent(inout)   :: dfdx
    type(scalar) , dimension(:,:), allocatable   :: dRdX
    integer                                      :: k, j
  
    if (.not.allocated(dRdX)) allocate(dRdX(this % nSVars, this % nDVars))

    dRdX = 0.0d0
    dfdx = 0.0d0

    !-----------------------------------------------------------------!
    ! Compute dfdx
    !-----------------------------------------------------------------!

    ! always use scale
    do k = 2, this % num_steps
       do j = 1, this % num_stages
          call this % system % func % addFuncDVSens(dfdx, this % h * this % B(j),&
               & this % T(j),  &
               & this % system % x, this % Q(k,j,:), this % QDOT(k,j,:), &
               & this % QDDOT(k,j,:) )
       end do
    end do

    !-----------------------------------------------------------------!
    ! Compute the total derivative
    !-----------------------------------------------------------------!

    do k = 2, this % num_steps
       do j = 1, this % num_stages

          call this % system % getResidualDVSens(dRdX, this % h * this % B(j), this % T(j), &
               & this % system % x, this % Q(k,j,:), this % QDOT(k,j,:), this % QDDOT(k,j,:))

          dfdx = dfdx + matmul(this % lam(k,j,:), drdx) ! check order
          
       end do
    end do

    !-----------------------------------------------------------------!
    ! Special logic for initial condition (use the adjoint variable
    ! for the second time step)
    !-----------------------------------------------------------------!

!!$    call this % system % func % addFuncDVSens(dfdx, 1.0d0, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(2,:) )
!!$    
!!$    call this % system % getResidualDVSens(dRdX, 1.0d0, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(2,:))
!!$    dfdx = dfdx + matmul(this % lam(2,:), dRdX)

    deallocate(dRdX)
 
  end subroutine computeTotalDerivative

  ! loop second order
  subroutine testAdjoint6(this, num_func, func, num_dv, x, dfdx, dfdxTmp)

    use function_class , only : abstract_function

    class(RK)                                  :: this
    class(abstract_function)                   :: func
    type(scalar), dimension(:), intent(inout)  :: x
    integer, intent(in)                        :: num_func, num_dv
    type(scalar), allocatable                  :: mat(:,:), tmpmat(:,:)
    type(scalar), allocatable                  :: rhs(:)
    type(scalar)                               :: alpha, beta, gamma
    type(scalar), allocatable                  :: dfdx(:), dfdxTmp(:)
    integer                                    :: size, info, ipiv
    type(scalar) , dimension(:,:), allocatable :: dRdX
    type(scalar)                               :: fvals, fvalstmp
    integer                                    :: k, ii, j, p

    allocate(rhs(this%nsvars))
    allocate(mat(this%nsvars, this%nsvars))
    allocate(tmpmat(this%nsvars, this%nsvars))

    !-----------------------------------------------------------------!
    ! Set the objective function into the system
    !-----------------------------------------------------------------!

    call this % system % setFunction(func)

    !-----------------------------------------------------------------!
    ! Set the number of variables, design variables into the system
    !-----------------------------------------------------------------!

    if (num_dv .ne. this % system % num_design_vars) stop "NDV mismatch"

    call this % system % setDesignVars(num_dv, x)

    this % nDVars = num_dv

    !-----------------------------------------------------------------!
    ! Integrate forward in time to solve for the state variables
    !-----------------------------------------------------------------!

    call this % integrate()

    do k = this % num_steps, 2, -1
       
       if (k .lt. this % num_steps) then

          rhs  = this % psi(k+1,:)
          
          do ii = 1, this % num_stages

             ! Add the psi terms
             alpha    = this % h * this % B(ii)
             beta     = 0.0d0
             gamma    = 0.0d0

             call this % system % assembleJacobian(mat, alpha, beta, gamma, &
                  & this % T(ii), this % Q(k+1,ii,:), this % qdot(k+1,ii,:), this % qddot(k+1,ii,:))
             
             rhs = rhs + matmul(transpose(mat),this % lam(k+1,ii,:))

             ! Add function contributions too
             call this % system % func % addFuncSVSens(rhs(1:1), alpha, beta, gamma, &
                  & this % T(ii), this % system % x, &
                  & this % Q(k+1,ii,:), this % qdot(k+1,ii,:), this % qddot(k+1,ii,:))

          end do
          
          this % psi(k,:) = rhs/1.0d0 ! note the positive sign
          
          !-----------------------------------------------------------!
          ! Add the phi terms
          !-----------------------------------------------------------!

          rhs  = this % phi(k+1,:)
          
          rhs = rhs + this % h * this % psi(k+1,:)

          do ii = 1, this % num_stages
             
             alpha = this % h * this % B(ii) * this % h * this % C(ii)
             beta  = this % h * this % B(ii)
             gamma = 0.0d0

             call this % system % assembleJacobian(mat, alpha, beta, gamma, &
                  & this % T(ii), this % Q(k+1,ii,:), this % qdot(k+1,ii,:), this % qddot(k+1,ii,:))

             rhs = rhs + matmul(transpose(mat), this % lam(k+1,ii,:))

             ! Add function contributions too
             call this % system % func % addFuncSVSens(rhs(1:1), alpha, beta, gamma, &
                  & this % T(ii), this % system % x, &
                  & this % Q(k+1,ii,:), this % qdot(k+1,ii,:), this % qddot(k+1,ii,:))

          end do
          
          this % phi(k,:) = rhs/1.0d0 ! note the positive sign
          
       else

          this % psi (k,:) = 0.0d0
          this % phi (k,:) = 0.0d0

       end if

!       print *, "PSI: ", k, this % psi(k,:)
!       print *, "PHI: ", k, this % phi(k,:)
       
       !-----------------------------------------------------------!
       ! Compute the stage adjoint variables
       !-----------------------------------------------------------!

       do ii = this % num_stages, 1, -1

          !-----------------------------------------------------------------!
          ! LAMBDA 22
          !-----------------------------------------------------------------!

          ! Assemble Jacobian
          mat      = 0.0d0

          alpha    = this % B(ii) * this % h * this % A(ii,ii) * this % h * this % A(ii,ii)
          beta     = this % B(ii) * this % h * this % A(ii,ii)
          gamma    = this % B(ii) * 1.0d0

          call this % system % assembleJacobian(mat, alpha, beta, gamma, &
               & this % T(ii), this % Q(k,ii,:), this % QDOT(k,ii,:), this % QDDOT(k,ii,:))
          
          ! Assemble RHS
          rhs = 0.0d0

          call this % system % func % addFuncSVSens(rhs, alpha, beta, gamma,  &
               & this % T(ii), this % system % X, &
               & this % Q(k,ii,:), this % QDOT(k,ii,:), this % QDDOT(k,ii,:))
          
          do j = ii+1, this % num_stages

             gamma = 0.0d0
             beta  = this % B(j) * this % h * this % A(j,ii)
             alpha = 0.0d0
             do p = ii, j
                alpha = alpha + this % A(j,p) * this % A(p,ii)
             end do
             alpha = alpha * this % h**2 * this % b(j)
             
             call this % system % assembleJacobian(tmpmat, alpha, beta, gamma, &
                  & this % T(j), this % Q(k,j,:), this % qdot(k,j,:), this % qddot(k,j,:))

             rhs(:) = rhs(:) + matmul(transpose(tmpmat(:,:)),this % lam(k,j,:))

             ! Add function contribution from next stage
             call this % system % func % addFuncSVSens(rhs, alpha, beta, gamma,  &
                  & this % T(j), this % system % x, &
                  & this % Q(k,j,:), this % qdot(k,j,:), this % qddot(k,j,:))

          end do

          rhs = rhs + this % B(ii) * this % phi(k,:)
          
          do j = ii, this % num_stages
             rhs = rhs + this % h * this % B(j) * this % A(j,ii) * this % psi(k,:)
          end do

          ! Solve for mu22
          this % lam(k,ii,:) = solve(mat, -rhs)
          
       end do

    end do

    ! Compute the adjoint total derivative
    call this % computeTotalDerivative(dfdx)

    deallocate(mat, tmpmat, rhs)
    
  end subroutine testAdjoint6

  subroutine evalFuncDIRK(this, x, fval)

    class(DIRK)                               :: this
    type(scalar), dimension(:), intent(in)    :: x
    type(scalar), intent(inout)               :: fval
    type(scalar)                              :: ftmp
    integer                                   :: k, j
    
    fval = 0.0d0
    ftmp = 0.0d0

!!$    do concurrent(k = 2 : this % num_steps)
!!$       call this % system % func % getFunctionValue(ftmp(k), this % time(k), &
!!$            & x, this % U(k,:), this % UDOT(k,:), this % UDDOT(k,:))
!!$    end do
!!$    
!!$    ! fval = sum(ftmp)/dble(this % num_steps)
!!$    fval = this % h*sum(ftmp)
!!$   
    do concurrent(k = 2 : this % num_steps)
       do concurrent(j = 1 : this % num_stages)
          call this % system % func % getFunctionValue(ftmp, this % T(j), x, &
               & this % Q(k,j,:), this % QDOT(k,j,:), this % QDDOT(k,j,:))
          fval = fval +  this % h * this % B(j) * ftmp
       end do
    end do

  end subroutine evalFuncDIRK

end module runge_kutta_integrator

