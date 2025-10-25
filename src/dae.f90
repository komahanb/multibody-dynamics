#include "scalar.fpp"
!=====================================================================!
! Module that provides the ability to the user to implement random
! ODEs and use with the DIRK scheme
!=====================================================================!

module dae_class

  use physics_class,  only : physics
  use function_class, only : abstract_function

  implicit none
  
  private

  public :: DAE
  
  !-------------------------------------------------------------------!
  ! Type that models a spring mass damper system (single mass)
  !-------------------------------------------------------------------!
  
  type, extends(physics) :: DAE

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     type(scalar) :: m = 1.0d0
     type(scalar) :: c = 0.2d0
     type(scalar) :: k = 5.0d0

   contains

     procedure :: mapDesignVars
     procedure :: assembleResidual
     procedure :: assembleJacobian
     procedure :: getInitialStates
     procedure :: getResidualDVSens

  end type DAE

contains
  
  !-------------------------------------------------------------------!
  ! Map the the design variables into the class variables
  !-------------------------------------------------------------------!
  
  subroutine mapDesignVars(this)

    class(dae) :: this
        
    this % M = this % x(1)
    this % C = this % x(2)
    this % K = this % x(3)

  end subroutine mapDesignVars

  !---------------------------------------------------------------------!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and udot
  !---------------------------------------------------------------------!

  subroutine getInitialStates(this, time, u, udot)

    class(dae) :: this

    real(dp), intent(in) :: time
    type(scalar), intent(inout), dimension(:) :: u, udot

    u(1) = 1.5d0 ! Actual DOF
    u(2) = 0.0d0 ! DUMMY DOF

  end subroutine getInitialStates

  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! -------------------------------------------------------------------!

  subroutine assembleResidual(this, res, time, u, udot, uddot)

    class(dae)                               :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in)                      :: time
    type(scalar), intent(in), dimension(:)    :: u, udot, uddot
    type(scalar), parameter :: omega = 1.0d0

    res(1) = this % m * uddot(1) + this % c * udot(1) + this % k * u(1)
    res(2) = exp(omega*u(2))*(u(1)-1.0d0) ! actual constraint is algebraic: (u(1)-1.0d0)

  end subroutine assembleResidual

  !-------------------------------------------------------------------!
  ! Jacobian assembly at each time step. If you don't provide the
  ! analytical jacobian, set setApproximateJacobian(.true.) into the
  ! integrator object. We use finite-difference method to approximate
  ! the Jacobian.
  !
  ! Jacobian is the matrix of partial derivatives. Each row in the
  ! Jacobian matrix arises from differntiating a single equation. Each
  ! column in the Jacobian comes from a variable in the problem. Note
  ! the NEQN should be equal to NVARS for the system to be solved.
  !
  ! Note: alpha, beta and gamma are scalars that need to be multiplied
  ! with the partial derivatives DRDQ, DRDudot and DRDuddot
  ! respectively.
  ! -------------------------------------------------------------------!

  subroutine assembleJacobian( this, jac, alpha, beta, gamma, &
       & time, u, udot, uddot )

    class(dae)                                 :: this
    type(scalar), intent(inout), dimension(:,:) :: jac
    type(scalar), intent(in)                    :: alpha, beta, gamma
    real(dp), intent(in)                        :: time
    type(scalar), intent(in), dimension(:)      :: u, udot, uddot
    type(scalar), parameter :: omega = 1.0d0

    ! Zero all  entries first
    jac = 0.0d0

    ! First row
    jac(1,1) = gamma * this % m + beta * this % c + alpha * this % k

    ! Second row
    jac(2,1) = alpha*exp(omega*u(2))
    jac(2,2) = alpha*exp(u(2))*(u(1)-1.0d0)*omega
    
  end subroutine assembleJacobian

  !----------------------------------------------------------------!
  ! Routine for evaluating the gradient of Residual with respect
  ! to the design X
  !----------------------------------------------------------------!
  
  subroutine getResidualDVSens(this, jac, scale, time, x, u, udot, uddot)

    class(dae)                            :: this
    type(scalar), intent(inout), dimension(:,:) :: jac
    real(dp), intent(in) :: time
    type(scalar), intent(in), dimension(:)      :: x, u, udot, uddot
    type(scalar)                                :: scale

    jac = 0.0d0 

    jac(1,1) = scale*uddot(1)
    jac(1,2) = scale*udot(1)
    jac(1,3) = scale*u(1)

    stop"Unimplemented"

  end subroutine getResidualDVSens

end module dae_class
