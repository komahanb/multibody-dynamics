!=====================================================================!
! Module that provides the ability to the user to implement random
! ODEs and use with the DIRK scheme
!=====================================================================!

module myode_class

  use physics_class, only : physics
  
  implicit none

  private

  public :: ODE

  !-------------------------------------------------------------------!
  ! Example type that models the physics
  !-------------------------------------------------------------------!

  type, extends(physics) :: ODE

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     real(8) :: m = 1.0d0
     real(8) :: c = 0.02d0
     real(8) :: k = 5.0d0

   contains

     procedure :: assembleResidual
     procedure :: assembleJacobian
     procedure :: getInitialStates

  end type ODE

contains
  
  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! -------------------------------------------------------------------!
  
  subroutine assembleResidual(this, res, time, u, udot, uddot)

    class(ODE) :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in) :: time
    real(8), intent(in), dimension(:) :: u, udot, uddot

    res(1) = this % m * uddot(1) + this % c * udot(1) + this % k * u(1)
    
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

    class(ODE) :: this
    real(8), intent(inout), dimension(:,:) :: jac
    real(8), intent(in) :: alpha, beta, gamma
    real(8), intent(in) :: time
    real(8), intent(in), dimension(:) :: u, udot, uddot

    ! Zero all  entries first
    jac = 0.0d0

    jac(1,1) = gamma * this % m + beta * this % c + alpha * this % k

  end subroutine assembleJacobian

  !---------------------------------------------------------------------!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and udot
  !---------------------------------------------------------------------!

  subroutine getInitialStates(this, time, u, udot)

    class(ODE) :: this

    real(8), intent(in) :: time
    real(8), intent(inout), dimension(:) :: u, udot

    u(1) = 1.0d0

  end subroutine getInitialStates

end module myode_class

