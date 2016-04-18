!=====================================================================!
! Module that provides the ability to the user to implement random
! ODEs and use with the DIRK scheme
!=====================================================================!

module myode_class

  ! parent class
  use physics_class, only : physics

  ! any other classes
  use utils, only : vector, matrix, skew, unskew, &
       & operator(*), operator(+), operator(-), &
       & norm, array, dp

  implicit none

  private

  public :: ODE
  
  !-------------------------------------------------------------------!
  ! Type that models rigid body dynamics
  !-------------------------------------------------------------------!
  
  type, extends(physics) :: ODE

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     real(dp) :: mass = 1.0d0
     real(dp) :: K    = 5.0d0
     real(dp) :: C    = 0.01d0

   contains

     procedure :: assembleResidual
     procedure :: assembleJacobian
     procedure :: getInitialStates

  end type ODE

contains

  !-------------------------------------------------------------------!
  ! Residual assembly at each time step
  !-------------------------------------------------------------------!

subroutine assembleResidual(this, res, time, u, udot, uddot)

 class(ODE) :: this

 real(8), intent(inout), dimension(:) :: res
 real(8), intent(in) :: time
 real(8), intent(in), dimension(:) :: u, udot, uddot

 res(1) = this % mass * uddot(1) + this % c * udot(1) + this % K * u(1)

end subroutine assembleResidual

!-------------------------------------------------------------------!
! Jacobian assembly at each time step
!-------------------------------------------------------------------!

subroutine assembleJacobian(this, jac, alpha, beta, gamma, &
    & time, u, udot, uddot)

 class(ODE) :: this
 real(8), intent(inout), dimension(:,:) :: jac
 real(8), intent(in) :: alpha, beta, gamma
 real(8), intent(in) :: time
 real(8), intent(in), dimension(:) :: u, udot, uddot

end subroutine assembleJacobian

!---------------------------------------------------------------------!
! Sets the initial condition for use in the integator. If first order
! system just set initial Q, if a second order system set initial Q
! and QDOT
!---------------------------------------------------------------------!

subroutine getInitialStates(this, time, u, udot)

 class(ODE) :: this

 real(8), intent(in) :: time
 real(8), intent(inout), dimension(:) :: u, udot

 ! initial position
 u(1)    = 1.0d0

 ! initial velocity
 udot(1) = 2.0d0

end subroutine getInitialStates

end module myode_class
