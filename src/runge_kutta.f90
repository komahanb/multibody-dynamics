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

     integer :: num_stages    = 1 ! default number of stages
     integer :: order         = 2      ! order of accuracy, only for informatory purposes
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

   contains

     ! Routines for integration     
     procedure :: integrate
     
     ! Destructor
     procedure :: finalize

     ! Helper routines
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff
     procedure, private :: TimeMarch
     procedure, private :: CheckButcherTableau
     
     procedure(buthcher_interface), private, deferred :: SetupButcherTableau
     
  end type RK

  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(RK) :: DIRK

   contains

     private

     procedure :: setupButcherTableau => ButcherDIRK

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

    class(RK) :: this

    ! Clear butcher's tableau
    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

    ! Clear stage state values
    if(allocated(this % QDDOT)) deallocate(this % QDDOT)
    if(allocated(this % QDOT)) deallocate(this % QDOT)
    if(allocated(this % Q)) deallocate(this % Q)
    if(allocated(this % T)) deallocate(this % T)

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
  
end module runge_kutta_integrator
