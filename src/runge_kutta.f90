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

     !----------------------------------------------------------------!
     ! Implemented common procedures (visible to the user)
     !----------------------------------------------------------------!

     procedure :: finalize, integrate

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

     procedure :: testAdjoint
     procedure :: testAdjoint2
     procedure :: testAdjoint3
     procedure :: testAdjoint4
     procedure :: testAdjoint5

  end type RK

  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(RK) :: DIRK

   contains

     private

     procedure :: setupButcherTableau => ButcherDIRK
     procedure, private :: computeStageStateValues

     !----------------------------------------------------------------!
     ! Adjoint procedures
     !----------------------------------------------------------------!
     procedure, public  :: marchBackwards
     procedure, private :: assembleRHS
     procedure, public :: computeTotalDerivative
     procedure, public :: evalFunc => evalFuncDIRK

  end type DIRK

  !===================================================================!
  ! Interfaces for deferred specialized procedures 
  !===================================================================!

  interface

     !================================================================!
     ! Interface for finding the stage derivatives at each time step
     !================================================================!

     subroutine computeStageStateValues_interface(this, q, qdot)
       import RK
       class(RK) :: this
       type(scalar), intent(in), dimension(:,:)           :: q
       type(scalar), OPTIONAL, intent(in), dimension(:,:) :: qdot
     end subroutine computeStageStateValues_interface

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
    
    allocate(this % lam(this % num_steps, this % num_stages, this % nsvars))
    this % lam = 0.0d0

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
    ! Put values into the Butcher tableau
    !-----------------------------------------------------------------!

    call this % setupButcherTableau()

    !-----------------------------------------------------------------!
    ! Sanity check for consistency of Butcher Tableau
    !-----------------------------------------------------------------!

    call this % checkButcherTableau()

    ! set the start time
    this % time(1) = this % tinit

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

    if(allocated(this % lam)) deallocate(this % lam)
    if(allocated(this % psi)) deallocate(this % psi)

  end subroutine finalize

  !===================================================================!
  ! Time integration logic
  !===================================================================!

  subroutine Integrate(this)

    class(RK) :: this
    integer :: k

    ! Set states to zero

    this % U     = 0.0d0
    this % UDOT  = 0.0d0
    this % UDDOT = 0.0d0
    this % time  = 0.0d0

    this % Q     = 0.0d0
    this % QDOT  = 0.0d0
    this % QDDOT = 0.0d0
    this % T     = 0.0d0

    call this % system % getInitialStates(this % time(1), &
         & this % u(1,:), this % udot(1,:))
    
    this % current_step = 1

    ! March in time
    time: do k = 2, this % num_steps

       this % current_step = k

       ! Find the stage derivatives at the current step
       call this % computeStageStateValues(this % u, this % udot)

       ! Advance the state to the current step
       call this % timeMarch(this % u, this % udot, this % uddot)
       
    end do time

  end subroutine Integrate
  
  !===================================================================!
  ! Time backwards in stage and backwards in time to solve for the
  ! adjoint variables
  !===================================================================!
  
  subroutine marchBackwards(this)

    class(DIRK) :: this
    integer     :: k, i
    type(scalar)    :: alpha, beta, gamma

    time: do k = this % num_steps, 2, -1
       
       this % current_step = k
       
       stage: do i = this % num_stages, 1, -1

          this % current_stage = i
          
          !--------------------------------------------------------------!
          ! Determine the linearization coefficients for the Jacobian
          !--------------------------------------------------------------!
          
          if (this % second_order) then
             gamma = this % B(i) * 1.0d0
             beta  = this % B(i) * this % h * this % A(i,i)
             alpha = this % B(i) * this % h * this % h * this % A(i,i) * this % A(i,i) 
          else
             stop "Reformualte first order"
             gamma = 0
             beta  = this % B(i) / this % h
             alpha = this % B(i) * this % A(i,i)
          end if
          
          !--------------------------------------------------------------!
          ! Solve the adjoint equation at each step
          !--------------------------------------------------------------!

          call this % adjointSolve(this % lam(k,i,:), alpha, beta, gamma, &
               & this % T(i), this % Q(k,i,:), this % QDOT(k,i,:), this % QDDOT(k,i,:))
          
          ! Find the adjoint variable for each time step          
          this % psi(k,:) = this % psi(k,:) + this % B(i) * this % lam (k,i,:)
          
       end do stage
       
       if ( k .lt. this % num_steps) then

!          this % psi(k,:) = this % psi(k,:) + this % h * this % psi (k+1,:) 

       end if

    end do time
    
  end subroutine marchBackwards

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
  ! Get the stage derivative array for the current step and states for
  ! DIRK
  !===================================================================!

  subroutine computeStageStateValues( this, q, qdot )

    class(DIRK)                                    :: this
    type(scalar), intent(in), dimension(:,:)           :: q
    type(scalar), OPTIONAL, intent(in), dimension(:,:) :: qdot
    integer                                        :: k, j, m
    type(scalar)                                       :: alpha, beta, gamma

    k = this % current_step

    do j = 1, this % num_stages

       this % current_stage = j

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
          alpha = this % h * this % A(j,j)
       end if

       call this % newtonSolve(alpha, beta, gamma, &
            & this % time(k), this % q(k,j,:), this % qdot(k,j,:), this % qddot(k,j,:))
       
    end do

  end subroutine computeStageStateValues
 
  !===================================================================!
  ! Function that puts together the right hand side of the adjoint
  ! equation into the supplied rhs vector.
  !===================================================================!

  subroutine assembleRHS(this, rhs)
 
    class(DIRK)                           :: this
    type(scalar), dimension(:), intent(inout) :: rhs
    type(scalar)                              :: scale1=0.0d0, scale2=0.0d0
    type(scalar), dimension(:,:), allocatable :: jac1, jac2
    integer :: k, j, i, p, s

    k = this % current_step
    i = this % current_stage
    s = this % num_stages

    allocate( jac1(this % nSVars, this % nSVars)  )
    allocate( jac2(this % nSVars, this % nSVars) )

    rhs = 0.0d0
    
    !-----------------------------------------------------------------!
    ! Add all the residual contributions first
    !-----------------------------------------------------------------!

    current_r: do j = i + 1, s

       scale1 = this % A(j,i) * this % h
       call this % system % assembleJacobian(jac1, ZERO, scale1, ZERO, &
            & this % T(j), this % Q(k,j,:), this % QDOT(k,j,:), this % QDDOT(k,j,:))

       scale2 = 0.0d0
       do p = i, j
          scale2 = scale2 +  this % A(p,i) * this % h
       end do
       call this % system % assembleJacobian(jac2,  scale1*scale2, ZERO, ZERO, &
            & this % T(j), this % Q(k,j,:), this % QDOT(k,j,:), this % QDDOT(k,j,:))

       rhs = rhs + matmul(transpose(jac1+jac2), this % lam(k,j,:))
       
    end do current_r

    
    if ( k+1 .le. this % num_steps ) then 

!!$       future_r: do j = i , s
!!$
!!$          scale1 = this % B(j) * this % B(i) / this % h
!!$
!!$          call this % system % assembleJacobian(jac1, ZERO, scale1, ZERO, &
!!$               & this % T(j), this % Q(k+1,j,:), this % QDOT(k+1,j,:), this % QDDOT(k+1,j,:))
!!$
!!$          scale2 = 0.0d0
!!$          do p = i , s
!!$             scale2 = scale2 + this % A(p,i) * this % B(p)
!!$          end do
!!$          do p = 1 , j
!!$             scale2 = scale2 + this % B(i) * this % A(j,p)
!!$          end do
!!$
!!$          scale2 = scale2 * this % B(j) 
!!$
!!$          call this % system % assembleJacobian(jac2, scale2, ZERO, ZERO, &
!!$               & this % T(j), this % Q(k+1,j,:), this % QDOT(k+1,j,:), this % QDDOT(k+1,j,:))
!!$
!!$          rhs = rhs + matmul( transpose(jac1+jac2), this % lam(k+1,j,:) )
!!$
!!$       end do future_r

    end if

    !-----------------------------------------------------------------!
    ! Now add function contributions
    !-----------------------------------------------------------------!
    
    ! Add contribution from second derivative of state
    
    scale1 = 1.0d0
    call this % system % func % addDFdUDDot(rhs, scale1, this % T(i), &
         & this % system % x, this % Q(k,i,:), this % qdot(k,i,:), this % qddot(k,i,:))
    
    current_f: do j = i, s

       scale1 = this % A(j,i) * this % h
       call this % system % func % addDFdUDot(rhs, scale1, this % T(j), &
            & this % system % x, this % Q(k,j,:), this % qdot(k,j,:), this % qddot(k,j,:))
       
       scale2 = 0.0d0
       do p = i, j
          scale2 = scale2 + this % A(p,i) * this % h
       end do
       call this % system % func % addDFdU(rhs, scale1*scale2, this % T(j), &
            & this % system % x, this % Q(k,j,:), this % Qdot(k,j,:), this % Qddot(k,j,:))
       
    end do current_f
    
    if ( k+1 .le. this % num_steps ) then 
!!$
!!$       future_f: do j = i , s
!!$
!!$          scale1 = this % B(j) * this % B(i) / this % h
!!$          call this % system % func % addDFdUDot(rhs, scale1, this % T(j), &
!!$               & this % system % x, this % Q(k+1,j,:), this % qdot(k+1,j,:), this % qddot(k+1,j,:))
!!$          
!!$          scale2 = 0.0d0
!!$          do p = i , s
!!$             scale2 = scale2 + this % A(p,i) * this % B(p)
!!$          end do
!!$          do p = 1 , j
!!$             scale2 = scale2 + this % B(i) * this % A(j,p)
!!$          end do
!!$
!!$          scale2 = scale2 * this % B(j)
!!$          call this % system % func % addDFdU(rhs, scale2, this % T(j), &
!!$               & this % system % x, this % Q(k+1,j,:), this % Qdot(k+1,j,:), this % Qddot(k+1,j,:))
!!$
!!$       end do future_f

    end if    

    rhs = -rhs

    if(allocated(jac1)) deallocate(jac1)
    if(allocated(jac2)) deallocate(jac2)

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
  
    allocate(dRdX(this % nSVars, this % nDVars))
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

!          print*, "Adding:", drdx, this % lam (k,j,:)

          dfdx = dfdx + matmul(this % lam(k,j,:), dRdX) ! check order
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
  
  subroutine testAdjoint5(this, num_func, func, num_dv, x, dfdx, dfdxTmp)

    use function_class , only : abstract_function

    class(RK)                                  :: this
    class(abstract_function)                   :: func
    type(scalar), dimension(:), intent(inout)  :: x
    integer, intent(in)                        :: num_func, num_dv
    type(scalar)                               :: mat(1,1) = 0.0d0, tmpmat(1,1) = 0.0d0
    type(scalar)                               :: rhs(1) = 0.0d0
    type(scalar)                               :: alpha, beta, gamma, dfdx(3), dfdxTmp(3)
    integer                                    :: size = 2 ,info, ipiv
    type(scalar) , dimension(:,:), allocatable :: dRdX
    type(scalar)                               :: fvals, fvalstmp, scale
    integer                                    :: k, ii
  

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

    !-----------------------------------------------------------------!
    ! PSI4
    !-----------------------------------------------------------------!

    this % psi (4,:) = 0.0d0

    !-----------------------------------------------------------------!
    ! LAMBDA 42
    !-----------------------------------------------------------------!

    mat      = 0.0d0

    ! Assemble Jacobian
    alpha    = this % h * this % A(2,2)
    beta     = 1.0d0
    gamma    = 0.0d0

    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(2), this % Q(4,2,:), this % qdot(4,2,:), this % qddot(4,2,:))

    rhs = 0.0d0

    ! Assemble RHS
    scale  = 1.0d0 !this % h * this % B(2) ! use only if the LHS is scaled too
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % X, this % Q(4,2,:), this % qdot(4,2,:), this % qddot(4,2,:))

    rhs = rhs + this % psi(4,:)

    ! Solve for lambda22
    this % lam(4,2,:) = -rhs(1)/mat(1,1)

    print *, "LAMBDA 42: ", this % lam(4,2,:)

    !-----------------------------------------------------------------!
    ! LAMBDA 41
    !-----------------------------------------------------------------!

    mat = 0.0d0

    ! Assemble Jacobian
    alpha    = this % h * this % A(1,1) 
    beta     = 1.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(1), this % Q(4,1,:), this % qdot(4,1,:), this % qddot(4,1,:))

    rhs = 0.0d0

    ! Add function contribution from this stage
    scale = 1.0d0 !this % h * this % B(1) 
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(1), this % system % x, &
         & this % Q(4,1,:), this % qdot(4,1,:), this % qddot(4,1,:))

    ! Add function contribution from next stage
    alpha    = this % h * this % A(2,1) 
    beta     = 0.0d0
    gamma    = 0.0d0
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % x, &
         & this % Q(4,2,:), this % qdot(4,2,:), this % qddot(4,2,:))
 
    ! Add RHS contribution from the next stage
    alpha    = this % h * this % A(2,1)
    beta     = 0.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(tmpmat(1:1,1:1), alpha, beta, gamma, &
         & this % T(2), this % Q(4,2,:), this % qdot(4,2,:), this % qddot(4,2,:))

    rhs = rhs + this % lam(4,2,:)*tmpmat(1,1) 

    rhs = rhs + this % psi(4,:)
    
    ! Solve for lambda21
    this % lam(4,1,:) = - rhs(1) / mat(1,1)

    print *, "LAMBDA 41: ", this % lam(4,1,:)

    !-----------------------------------------------------------------!
    ! psi 3
    !-----------------------------------------------------------------!

    mat = 0.0d0
    rhs = 0.0d0

    alpha    = 1.0d0
    beta     = 0.0d0
    gamma    = 0.0d0

    scale    = this % h * this % B(1)

    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(1), this % system % x, &
         & this % Q(4,1,:), this % qdot(4,1,:), this % qddot(4,1,:))

    call this % system % assembleJacobian(mat(1:1,1:1), alpha*scale, beta*scale, gamma*scale, &
         & this % T(1), this % Q(4,1,:), this % qdot(4,1,:), this % qddot(4,1,:))

    rhs = rhs + mat(1,1) * this % lam(4,1,:)

    !-----------------------------------------------------------------!

    alpha    = 1.0d0
    beta     = 0.0d0
    gamma    = 0.0d0

    scale    = this % h * this % B(2)

    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % x, &
         & this % Q(4,2,:), this % qdot(4,2,:), this % qddot(4,2,:))

    call this % system % assembleJacobian(mat(1:1,1:1), alpha*scale, beta*scale, gamma*scale, &
         & this % T(2), this % Q(4,2,:), this % qdot(4,2,:), this % qddot(4,2,:))

    rhs = rhs + mat(1,1) * this % lam(4,2,:)
    
    this % psi (3,:) = rhs/1.0d0 ! negative sign cancels

    print *, "PSI3=", this % psi(3,:)

    !-----------------------------------------------------------------!
    ! LAMBDA 32
    !-----------------------------------------------------------------!

    mat      = 0.0d0

    ! Assemble Jacobian
    alpha    = this % h * this % A(2,2)
    beta     = 1.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(2), this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))

    ! Assemble RHS
    rhs = 0.0d0
    scale  = 1.0d0 !this % h * this % B(2) ! use only if the LHS is scaled too
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % X, &
         & this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))

    rhs = rhs + this % psi(3,:)

    ! Solve for lambda22
    this % lam(3,2,:) = -rhs(1)/mat(1,1)

    print *, "LAMBDA 32: ", this % lam(3,2,:)

    !-----------------------------------------------------------------!
    ! LAMBDA 31
    !-----------------------------------------------------------------!
    
    mat = 0.0d0

    ! Assemble Jacobian
    alpha    = this % h * this % A(1,1) 
    beta     = 1.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(1), this % Q(3,1,:), this % qdot(3,1,:), this % qddot(3,1,:))

    !-----------------------------------------------------------------!
    ! Assemble RHS
    !-----------------------------------------------------------------!

    ! Add function contribution from this stage
    rhs = 0.0d0

    scale = 1.0d0 !this % h * this % B(1) 
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(1), this % system % x, &
         & this % Q(3,1,:), this % qdot(3,1,:), this % qddot(3,1,:))
    
    ! Add function contribution from next stage
    alpha    = this % h * this % A(2,1) 
    beta     = 0.0d0
    gamma    = 0.0d0
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % x, &
         & this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))
 
    ! Add RHS contribution from the next stage
    alpha    = this % h * this % A(2,1)
    beta     = 0.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(tmpmat(1:1,1:1), alpha, beta, gamma, &
         & this % T(2), this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))

    rhs = rhs + this % lam(3,2,:)*tmpmat(1,1) 

    ! Add contribution from time multiplier
    rhs = rhs + this % psi(3,:)
    
    ! Solve for lambda21
    this % lam(3,1,:) = - rhs(1) / mat(1,1)

    print *, "LAMBDA 31: ", this % lam(3,1,:)

    !-----------------------------------------------------------------!
    ! PSI 2
    !-----------------------------------------------------------------!

    mat = 0.0d0
    rhs = 0.0d0

    alpha    = 1.0d0
    beta     = 0.0d0
    gamma    = 0.0d0

    scale    = this % h * this % B(1)

    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(1), this % system % x, &
         & this % Q(3,1,:), this % qdot(3,1,:), this % qddot(3,1,:))

    call this % system % assembleJacobian(mat(1:1,1:1), alpha*scale, beta*scale, gamma*scale, &
         & this % T(1), this % Q(3,1,:), this % qdot(3,1,:), this % qddot(3,1,:))

    rhs = rhs + mat(1,1) * this % lam(3,1,:)

    !-----------------------------------------------------------------!

    alpha    = 1.0d0
    beta     = 0.0d0
    gamma    = 0.0d0

    scale    = this % h * this % B(2)

    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % x, &
         & this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))

    call this % system % assembleJacobian(mat(1:1,1:1), alpha*scale, beta*scale, gamma*scale, &
         & this % T(2), this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))

    rhs = rhs + mat(1,1) * this % lam(3,2,:)

    rhs = rhs + this % psi(3,:) !# missed it

    this % psi (2,:) = rhs/1.0d0 ! negative sign cancels

    print *, "PSI2=", this % psi(2,:)

    !-----------------------------------------------------------------!
    ! LAMBDA 22
    !-----------------------------------------------------------------!

    mat      = 0.0d0

    ! Assemble Jacobian
    alpha    = this % h * this % A(2,2)
    beta     = 1.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(2), this % Q(2,2,:), this % qdot(2,2,:), this % qddot(2,2,:))

    ! Assemble RHS
    rhs = 0.0d0
    scale  = 1.0d0 !this % h * this % B(2) ! use only if the LHS is scaled too
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % X, &
         & this % Q(2,2,:), this % qdot(2,2,:), this % qddot(2,2,:))

    rhs = rhs + this % psi(2,:)

    ! Solve for lambda22
    this % lam(2,2,:) = -rhs(1)/mat(1,1)

    print *, "LAMBDA 22: ", this % lam(2,2,:)

    !-----------------------------------------------------------------!
    ! LAMBDA 21 
    !-----------------------------------------------------------------!

    mat = 0.0d0

    ! Assemble Jacobian
    alpha    = this % h * this % A(1,1) 
    beta     = 1.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(1), this % Q(2,1,:), this % qdot(2,1,:), this % qddot(2,1,:))

    !-----------------------------------------------------------------!
    ! Assemble RHS
    !-----------------------------------------------------------------!

    ! Add function contribution from this stage
    rhs = 0.0d0

    scale = 1.0d0 !this % h * this % B(1) 
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(1), this % system % x, &
         & this % Q(2,1,:), this % qdot(2,1,:), this % qddot(2,1,:))
    
    ! Add function contribution from next stage
    alpha    = this % h * this % A(2,1) 
    beta     = 0.0d0
    gamma    = 0.0d0
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % x, &
         & this % Q(2,2,:), this % qdot(2,2,:), this % qddot(2,2,:))
 
    ! Add RHS contribution from the next stage
    alpha    = this % h * this % A(2,1)
    beta     = 0.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(tmpmat(1:1,1:1), alpha, beta, gamma, &
         & this % T(2), this % Q(2,2,:), this % qdot(2,2,:), this % qddot(2,2,:))

    rhs = rhs + this % lam(2,2,:)*tmpmat(1,1) 

    ! Add contribution from time multiplier
    rhs = rhs + this % psi(2,:)
    
    ! Solve for lambda21
    this % lam(2,1,:) = - rhs(1) / mat(1,1)

    print *, "LAMBDA 21: ", this % lam(2,1,:)

    ! Compute the adjoint total derivative
    call this % computeTotalDerivative(dfdx)

    !-----------------------------------------------------------------!
    ! CSD check
    !-----------------------------------------------------------------!

    call this % evalFunc(x, fvals)
    print *, "FV=", fvals

    x(1) = cmplx(dble(x(1)), 1.0d-16)
    call this % system % setDesignVars(num_dv, x)
    call this % integrate()
    call this % evalFunc(x, fvalstmp)
    dfdxtmp(1) = aimag(fvalstmp)/1.0d-16
    x(1) = cmplx(dble(x(1)), 0.0d0)

    x(2) = cmplx(dble(x(2)), 1.0d-16)
    call this % system % setDesignVars(num_dv, x)
    call this % integrate()
    call this % evalFunc(x, fvalstmp)
    dfdxtmp(2) = aimag(fvalstmp)/1.0d-16
    x(2) = cmplx(dble(x(2)), 0.0d0)

    x(3) = cmplx(dble(x(3)), 1.0d-16)
    call this % system % setDesignVars(num_dv, x)
    call this % integrate()
    call this % evalFunc(x, fvalstmp)
    dfdxtmp(3) = aimag(fvalstmp)/1.0d-16
    x(3) = cmplx(dble(x(3)), 0.0d0)
    call this % system % setDesignVars(num_dv, x)

  end subroutine testAdjoint5

  subroutine testAdjoint4(this, num_func, func, num_dv, x, dfdx, dfdxTmp)

    use function_class , only : abstract_function

    class(RK)                                  :: this
    class(abstract_function)                   :: func
    type(scalar), dimension(:), intent(inout)  :: x
    integer, intent(in)                        :: num_func, num_dv
    type(scalar)                               :: mat(1,1) = 0.0d0, tmpmat(1,1) = 0.0d0
    type(scalar)                               :: rhs(1) = 0.0d0
    type(scalar)                               :: alpha, beta, gamma, dfdx(3), dfdxTmp(3)
    integer                                    :: size = 2 ,info, ipiv
    type(scalar) , dimension(:,:), allocatable :: dRdX
    type(scalar)                               :: fvals, fvalstmp, scale
    integer                                    :: k, ii

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

          rhs = 0.0d0

          rhs  = this % psi(k+1,:) ! note the positive sign

          do ii = this % num_stages, 1, -1

             scale = this % h * this % B(ii)

             mat = 0.0d0

             alpha    = 1.0d0
             beta     = 0.0d0
             gamma    = 0.0d0

             call this % system % assembleJacobian(mat(1:1,1:1), alpha*scale, beta*scale, gamma*scale, &
                  & this % T(ii), this % Q(k+1,ii,:), this % qdot(k+1,ii,:), this % qddot(k+1,ii,:))
             
             rhs = rhs + mat(1,1)*this % lam(k+1,ii,:)

             ! Add function contributions too
             call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale, &
                  & this % T(ii), this % system % x, &
                  & this % Q(k+1,ii,:), this % qdot(k+1,ii,:), this % qddot(k+1,ii,:))

          end do

          this % psi(k,:) = rhs/1.0d0 ! note the positive sign

       else

          this % psi (k,:) = 0.0d0

       end if

       print *, "PSI: ", k, this % psi(k,:)

       do ii = this % num_stages, 1, -1
          
          !-----------------------------------------------------------------!
          ! LAMBDA 22
          !-----------------------------------------------------------------!

          ! Assemble Jacobian
          mat      = 0.0d0

          alpha    = this % h * this % A(ii,ii)
          beta     = 1.0d0
          gamma    = 0.0d0
          call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
               & this % T(ii), this % Q(k,ii,:), this % qdot(k,ii,:), this % qddot(k,ii,:))

          ! Assemble RHS
          rhs = 0.0d0

          scale  = 1.0d0 !this % h * this % B(2) ! use only if the LHS is scaled too
          call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
               & this % T(ii), this % system % X, &
               & this % Q(k,ii,:), this % qdot(k,ii,:), this % qddot(k,ii,:))

          if ( ii .lt. this % num_stages) then

             tmpmat = 0.0d0

             ! Add RHS contribution from the next stage
             alpha    = this % h * this % A(ii+1,ii)
             beta     = 0.0d0
             gamma    = 0.0d0
             call this % system % assembleJacobian(tmpmat(1:1,1:1), alpha, beta, gamma, &
                  & this % T(ii+1), this % Q(k,ii+1,:), this % qdot(k,ii+1,:), this % qddot(k,ii+1,:))
             
             rhs = rhs + this % lam(k,ii+1,:)*tmpmat(1,1)

             ! Add function contribution from next stage
             alpha    = this % h * this % A(ii+1,ii) 
             beta     = 0.0d0
             gamma    = 0.0d0
             call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
                  & this % T(ii+1), this % system % x, &
                  & this % Q(k,ii+1,:), this % qdot(k,ii+1,:), this % qddot(k,ii+1,:))
             
          end if
          
          rhs = rhs + this % psi(k,:)

          ! Solve for lambda22
          this % lam(k,ii,:) = -rhs(1)/mat(1,1)

          print *, "LAMBDA: ", k,ii, this % lam(k,ii,:)
          
       end do

    end do

    ! Compute the adjoint total derivative
    call this % computeTotalDerivative(dfdx)

    !-----------------------------------------------------------------!
    ! CSD check
    !-----------------------------------------------------------------!

    call this % evalFunc(x, fvals)
    print *, "FV=", fvals

    x(1) = cmplx(dble(x(1)), 1.0d-16)
    call this % system % setDesignVars(num_dv, x)
    call this % integrate()
    call this % evalFunc(x, fvalstmp)
    dfdxtmp(1) = aimag(fvalstmp)/1.0d-16
    x(1) = cmplx(dble(x(1)), 0.0d0)

    x(2) = cmplx(dble(x(2)), 1.0d-16)
    call this % system % setDesignVars(num_dv, x)
    call this % integrate()
    call this % evalFunc(x, fvalstmp)
    dfdxtmp(2) = aimag(fvalstmp)/1.0d-16
    x(2) = cmplx(dble(x(2)), 0.0d0)

    x(3) = cmplx(dble(x(3)), 1.0d-16)
    call this % system % setDesignVars(num_dv, x)
    call this % integrate()
    call this % evalFunc(x, fvalstmp)
    dfdxtmp(3) = aimag(fvalstmp)/1.0d-16
    x(3) = cmplx(dble(x(3)), 0.0d0)
    call this % system % setDesignVars(num_dv, x)

  end subroutine testAdjoint4

  ! First order test code
  subroutine testAdjoint3(this, num_func, func, num_dv, x, dfdx, dfdxTmp)

    use function_class , only : abstract_function

    class(RK)                                  :: this
    class(abstract_function)                   :: func
    type(scalar), dimension(:), intent(inout)  :: x
    integer, intent(in)                        :: num_func, num_dv
    type(scalar)                               :: mat(1,1) = 0.0d0, tmpmat(1,1) = 0.0d0
    type(scalar)                               :: rhs(1) = 0.0d0
    type(scalar)                               :: alpha, beta, gamma, dfdx(3), dfdxTmp(3)
    integer                                    :: size = 2 ,info, ipiv
    type(scalar) , dimension(:,:), allocatable :: dRdX
    type(scalar)                               :: fvals, fvalstmp, scale
    integer                                    :: k, ii

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

    !-----------------------------------------------------------------!
    ! Solve for psi3
    !-----------------------------------------------------------------!

    this % psi (3,:) = 0.0d0

    !-----------------------------------------------------------------!
    ! LAMBDA 32
    !-----------------------------------------------------------------!

    mat      = 0.0d0

    ! Assemble Jacobian
    alpha    = this % h * this % A(2,2)
    beta     = 1.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(2), this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))

    ! Assemble RHS
    rhs = 0.0d0
    scale  = 1.0d0 !this % h * this % B(2) ! use only if the LHS is scaled too
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % X, &
         & this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))

    ! Solve for lambda22
    this % lam(3,2,:) = -rhs(1)/mat(1,1)

    print *, "LAMBDA 32: ", this % lam(3,2,:)

    !-----------------------------------------------------------------!
    ! [dR/dQdot31] LAMBDA 31 = -dF/dQdot21 -lam22 dR21/dqdot22
    !-----------------------------------------------------------------!

    mat = 0.0d0

    ! Assemble Jacobian
    alpha    = this % h * this % A(1,1) 
    beta     = 1.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(1), this % Q(3,1,:), this % qdot(3,1,:), this % qddot(3,1,:))

    !-----------------------------------------------------------------!
    ! Assemble RHS
    !-----------------------------------------------------------------!

    rhs = 0.0d0

    ! Add function contribution from this stage
    scale = 1.0d0 !this % h * this % B(1) 
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(1), this % system % x, &
         & this % Q(3,1,:), this % qdot(3,1,:), this % qddot(3,1,:))

    ! Add function contribution from next stage
    alpha    = this % h * this % A(2,1) 
    beta     = 0.0d0
    gamma    = 0.0d0
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % x, &
         & this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))
 
    ! Add RHS contribution from the next stage
    alpha    = this % h * this % A(2,1)
    beta     = 0.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(tmpmat(1:1,1:1), alpha, beta, gamma, &
         & this % T(2), this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))

    rhs = rhs + this % lam(3,2,:)*tmpmat(1,1) 

    rhs = rhs + this % psi(3,:)
    
    ! Solve for lambda21
    this % lam(3,1,:) = - rhs(1) / mat(1,1)

    print *, "LAMBDA 31: ", this % lam(3,1,:)

    !-----------------------------------------------------------------!
    ! PSI 2
    !-----------------------------------------------------------------!

    mat = 0.0d0
    rhs = 0.0d0

    alpha    = 1.0d0
    beta     = 0.0d0
    gamma    = 0.0d0

    scale    = this % h * this % B(1)

    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(1), this % system % x, &
         & this % Q(3,1,:), this % qdot(3,1,:), this % qddot(3,1,:))

    call this % system % assembleJacobian(mat(1:1,1:1), alpha*scale, beta*scale, gamma*scale, &
         & this % T(1), this % Q(3,1,:), this % qdot(3,1,:), this % qddot(3,1,:))

    rhs = rhs + mat(1,1) * this % lam(3,1,:)

    !-----------------------------------------------------------------!

    alpha    = 1.0d0
    beta     = 0.0d0
    gamma    = 0.0d0

    scale    = this % h * this % B(2)

    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % x, &
         & this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))

    call this % system % assembleJacobian(mat(1:1,1:1), alpha*scale, beta*scale, gamma*scale, &
         & this % T(2), this % Q(3,2,:), this % qdot(3,2,:), this % qddot(3,2,:))

    rhs = rhs + mat(1,1) * this % lam(3,2,:)

    this % psi (2,:) = rhs/1.0d0 ! negative sign cancels

    print *, "PSI2=", this % psi(2,:)

    !-----------------------------------------------------------------!
    ! LAMBDA 22
    !-----------------------------------------------------------------!

    mat      = 0.0d0

    ! Assemble Jacobian
    alpha    = this % h * this % A(2,2)
    beta     = 1.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(2), this % Q(2,2,:), this % qdot(2,2,:), this % qddot(2,2,:))

    ! Assemble RHS
    rhs = 0.0d0
    scale  = 1.0d0 !this % h * this % B(2) ! use only if the LHS is scaled too
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % X, &
         & this % Q(2,2,:), this % qdot(2,2,:), this % qddot(2,2,:))

    rhs = rhs + this % psi(2,:)

    ! Solve for lambda22
    this % lam(2,2,:) = -rhs(1)/mat(1,1)

    print *, "LAMBDA 22: ", this % lam(2,2,:)

    !-----------------------------------------------------------------!
    ! [dR/dQdot21] LAMBDA 21 = -dF/dQdot21 -lam22 dR21/dqdot22
    !-----------------------------------------------------------------!
    mat = 0.0d0

    ! Assemble Jacobian
    alpha    = this % h * this % A(1,1) 
    beta     = 1.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(1), this % Q(2,1,:), this % qdot(2,1,:), this % qddot(2,1,:))

    !-----------------------------------------------------------------!
    ! Assemble RHS
    !-----------------------------------------------------------------!

    ! Add function contribution from this stage
    rhs = 0.0d0

    scale = 1.0d0 !this % h * this % B(1) 
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(1), this % system % x, &
         & this % Q(2,1,:), this % qdot(2,1,:), this % qddot(2,1,:))
    
    ! Add function contribution from next stage
    alpha    = this % h * this % A(2,1) 
    beta     = 0.0d0
    gamma    = 0.0d0
    call this % system % func % addFuncSVSens(rhs(1:1), alpha*scale, beta*scale, gamma*scale,  &
         & this % T(2), this % system % x, &
         & this % Q(2,2,:), this % qdot(2,2,:), this % qddot(2,2,:))
 
    ! Add RHS contribution from the next stage
    alpha    = this % h * this % A(2,1)
    beta     = 0.0d0
    gamma    = 0.0d0
    call this % system % assembleJacobian(tmpmat(1:1,1:1), alpha, beta, gamma, &
         & this % T(2), this % Q(2,2,:), this % qdot(2,2,:), this % qddot(2,2,:))

    rhs = rhs + this % lam(2,2,:)*tmpmat(1,1) 

    ! Add contribution from time multiplier
    rhs = rhs + this % psi(2,:)
    
    ! Solve for lambda21
    this % lam(2,1,:) = - rhs(1) / mat(1,1)

    print *, "LAMBDA 21: ", this % lam(2,1,:)

    !-----------------------------------------------------------------!

    ! Compute the adjoint total derivative
    call this % computeTotalDerivative(dfdx)

    !-----------------------------------------------------------------!
    ! CSD check
    !-----------------------------------------------------------------!

    call this % evalFunc(x, fvals)
    print *, "FV=", fvals

    x(1) = cmplx(dble(x(1)), 1.0d-25)
    call this % system % setDesignVars(num_dv, x)
    call this % integrate()
    call this % evalFunc(x, fvalstmp)
    dfdxtmp(1) = aimag(fvalstmp)/1.0d-25
    x(1) = cmplx(dble(x(1)), 0.0d0)

    x(2) = cmplx(dble(x(2)), 1.0d-25)
    call this % system % setDesignVars(num_dv, x)
    call this % integrate()
    call this % evalFunc(x, fvalstmp)
    dfdxtmp(2) = aimag(fvalstmp)/1.0d-25
    x(2) = cmplx(dble(x(2)), 0.0d0)

    x(3) = cmplx(dble(x(3)), 1.0d-25)
    call this % system % setDesignVars(num_dv, x)
    call this % integrate()
    call this % evalFunc(x, fvalstmp)
    dfdxtmp(3) = aimag(fvalstmp)/1.0d-25
    x(3) = cmplx(dble(x(3)), 0.0d0)
    call this % system % setDesignVars(num_dv, x)

  end subroutine testAdjoint3

  subroutine testAdjoint2(this, num_func, func, num_dv, x, dfdx)

    use function_class , only : abstract_function

    class(RK)                                  :: this
    class(abstract_function)                   :: func
    type(scalar), dimension(:), intent(in)     :: x
    integer, intent(in)                        :: num_func, num_dv
    type(scalar)                               :: mat(1,1) = 0.0d0
    type(scalar)                               :: rhs(1) = 0.0d0
    type(scalar)                               :: alpha, beta, gamma, dfdx(3), dfdxTmp(3)
    integer                                    :: size = 2 ,info, ipiv, k
    type(scalar) , dimension(:,:), allocatable :: dRdX
    type(scalar)                               :: scale = 0.0d0

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

    !-----------------------------------------------------------------!
    ! Matrix
    !-----------------------------------------------------------------!

    do k = this % num_steps, 2, -1

       mat = 0.0d0
       
       alpha    = this % h * this % h * this % A(1,1)
       beta     = this % h
       gamma    = 1.0d0
       call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))

       rhs = 0.0d0

       alpha = this % h * this % A(1,1) * this % h 
       call this % system % func % addDFdU(rhs(1:1), alpha, this % time(k), &
            & this % system % x, this % U(k,:), this % UDOT(k,:), this % UDDOT(k,:))

       beta  = this % h
       call this % system % func % addDFdUDOT(rhs(1:1), beta, this % time(k), &
            & this % system % x, this % U(k,:), this % UDOT(k,:), this % UDDOT(k,:))

       gamma = 1.0d0
       call this % system % func % addDFdUDDOT(rhs(1:1), gamma, this % time(k), &
            & this % system % x, this % U(k,:), this % UDOT(k,:), this % UDDOT(k,:))

       if ( k + 1 .le. this % num_steps ) then
          alpha    = 1.0d0 + this % h * this % B(1)
          beta     = 1.0d0
          gamma    = 0.0d0       
          call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
               & this % time(k+1), this % u(k+1,:), this % udot(k+1,:), this % uddot(k+1,:))
          rhs = rhs + mat(1,1) * this % psi(k+1,:)
       end if

       this % psi(k,:) = -rhs(1)/mat(1,1)

!       print*, "psi=", k, this % psi(k,:)
       
    end do

    ! Compute total derivative    

    dfdx = 0.0d0

    allocate(dRdX(this % nSVars, this % nDVars))

    scale = this % h * this % B(1)

    do k = this % num_steps,  2, -1

       dfdxTmp = 0.0d0

       call this % system % func % addFuncDVSens(dfdxtmp, scale, &
            & this % Time(k),  &
            & this % system % x, this % U(k,:), this % UDOT(k,:), &
            & this % UDDOT(k,:) )

       dRdX = 0.0d0

       call this % system % getResidualDVSens(dRdX, scale, &
            & this % Time(k), &
            & this % system % x, this % U(k,:), this % UDOT(k,:), &
            & this % UDDOT(k,:))

       dfdx(:) = dfdx + dfdxTmp(:) + this % psi(k,1) * dRdX(1,:)

    end do

    deallocate(drdx)

  end subroutine testAdjoint2

  subroutine testAdjoint(this, num_func, func, num_dv, x, dfdx)

    use function_class , only : abstract_function

    class(RK)    :: this
    class(abstract_function) ::func
    type(scalar), dimension(:), intent(in)           :: x
    integer, intent(in)                              :: num_func, num_dv
    type(scalar) :: mat(2,2) = 0.0d0
    type(scalar) :: rhs(2) = 0.0d0
    type(scalar) :: alpha, beta, gamma, dfdx(3)
    integer      :: size = 2 ,info, ipiv(2)
    type(scalar) , dimension(:,:), allocatable   :: dRdX
    type(scalar) :: psi, lam

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

    !-----------------------------------------------------------------!
    ! Matrix
    !-----------------------------------------------------------------!

    mat = 0.0d0
    rhs = 0.0d0

    alpha    = this % h * this % h * this % A(1,1) * this % A(1,1)
    beta     = this % h * this % A(1,1)
    gamma    = 1.0d0
    call this % system % assembleJacobian(mat(1:1,1:1), alpha, beta, gamma, &
         & this % T(1), this % q(2,1,:), this % qdot(2,1,:), this % qddot(2,1,:))

    alpha    = this % h * this % h * this % A(1,1) * this % B(1)
    beta     = this % h * this % B(1)
    gamma    = this % B(1)
    call this % system % assembleJacobian(mat(2:2,1:1), alpha, beta, gamma, &
         & this % T(1), this % q(2,1,:), this % qdot(2,1,:), this % qddot(2,1,:))

!!$    alpha    = this % h * this % h * this % A(1,1) * this % A(1,1) / this % B(1)
!!$    beta     = this % h * this % A(1,1) / this % B(1)
!!$    gamma    = 1.0d0 / this % B(1)
!!$    call this % system % assembleJacobian(mat(1:1,2:2), alpha, beta, gamma, &
!!$         & this % time(2), this % u(2,:), this % udot(2,:), this % uddot(2,:))

    alpha    = this % h * this % h * this % A(1,1)
    beta     = this % h
    gamma    = 1.0d0
    call this % system % assembleJacobian(mat(2:2,2:2), alpha, beta, gamma, &
         & this % time(2), this % u(2,:), this % udot(2,:), this % uddot(2,:))

    mat = transpose(mat)

    print *, mat
    !-----------------------------------------------------------------!
    ! Assemble the RHS
    !-----------------------------------------------------------------!

    ! First term
    gamma = 1.0d0 
    call this % system % func % addDFdUDDot(rhs(1:1), gamma , this % T(1), &
         & this % system % x, this % Q(2,1,:), this % QDot(2,1,:), this % Qddot(2,1,:))

    beta = this % h * this % A(1,1)
    call this % system % func % addDFdUDot(rhs(1:1), beta , this % T(1), &
         & this % system % x, this % Q(2,1,:), this % QDot(2,1,:), this % Qddot(2,1,:))

    alpha = this % h * this % A(1,1)* this % h * this % A(1,1)
    call this % system % func % addDFdU(rhs(1:1), alpha , this % T(1), &
         & this % system % x, this % Q(2,1,:), this % QDot(2,1,:), this % Qddot(2,1,:))

    ! Second term
    alpha = this % h * this % A(1,1) * this % h 
    call this % system % func % addDFdU(rhs(2:2), alpha, this % time(2), &
         & this % system % x, this % U(2,:), this % UDOT(2,:), this % UDDOT(2,:))

    beta  = this % h
    call this % system % func % addDFdUDOT(rhs(2:2), beta, this % time(2), &
         & this % system % x, this % U(2,:), this % UDOT(2,:), this % UDDOT(2,:))

    gamma = 1.0d0
    call this % system % func % addDFdUDDOT(rhs(2:2), gamma, this % time(2), &
         & this % system % x, this % U(2,:), this % UDOT(2,:), this % UDDOT(2,:))

    rhs = rhs
    ! print*, rhs(1)
    ! print*, rhs(2)

    ! Call linear solve
#if defined USE_COMPLEX
    call ZGESV(size, 1, mat, size, IPIV, rhs, size, INFO)
#else
    call DGESV(size, 1, mat, size, IPIV, rhs, size, INFO)
#endif
    if (info.ne.0) then
       print *, info
       stop"lapack error"
    end if

    this % lam(2,1,1) = rhs(1)

    print *, "lam1=", rhs(1)
    print *, "lam2=", rhs(2)

    !-----------------------------------------------------------------!
    ! Evaluate total derivative
    !-----------------------------------------------------------------!
    dfdx = 0.0d0

    call this % system % func % addFuncDVSens(dfdx, this % h * this % B(1),&
         &  this % T(1),  &
         & this % system % x, this % Q(2,1,:), this % QDOT(2,1,:), &
         & this % QDDOT(2,1,:) )

    allocate(dRdX(this % nSVars, this % nDVars))
    dRdX = 0.0d0

    call this % system % getResidualDVSens(dRdX, this % h * this % B(1), this % T(1), &
         & this % system % x, this % Q(2,1,:), this % QDOT(2,1,:), this % QDDOT(2,1,:))

    dfdx(:) = dfdx(:) + this % lam(2,1,1)*dRdX(1,:)

    deallocate(drdx)

  end subroutine testAdjoint


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

