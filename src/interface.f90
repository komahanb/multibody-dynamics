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
     real(dp) :: C    = 0.02d0

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

 res(1) = this % mass * uddot(1) + this % c * udot(1) + this % K * u(1)

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
! with the partial derivatives DRDQ, DRDQDOT and DRDQDDOT
! respectively.
! -------------------------------------------------------------------!

subroutine assembleJacobian(this, jac, alpha, beta, gamma, &
    & time, u, udot, uddot)

 class(ODE) :: this
 real(8), intent(inout), dimension(:,:) :: jac
 real(8), intent(in) :: alpha, beta, gamma
 real(8), intent(in) :: time
 real(8), intent(in), dimension(:) :: u, udot, uddot

 ! Derivative of the first equation with respect to the states 
 ! (alpha*DRDQ + beta*DRDQDOT + gamma*DRDQDDOT)

 ! first equation with repect to the first variable q(1), qdot(1) and qddot(1)
 JAC(1,1) =  alpha*this % K +  beta*this % C+ gamma*this % mass

 ! first equation with repect to the first variable q(1), qdot(1) and qddot(1)
 ! JAC(1,2) =  alpha*this % K +  beta*this % C+ gamma*this % mass

 ! Derivative of the second equation with respect to the states 
 ! (alpha*DRDQ + beta*DRDQDOT + gamma*DRDQDDOT)

 ! second equation with repect to the first variable q(1), qdot(1) and qddot(1)
 ! J(2,1) = - alpha*0.05d0*qdot(2) + beta*0.05d0*q(1)
 
 

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
