!=====================================================================!
! Module that provides the ability to the user to implement random
! ODEs and use with the DIRK scheme
!=====================================================================!

module spring_mass_damper_class

  use physics_class, only : physics
  
  implicit none
  
  private

  public :: smd1, smd2
  
  !-------------------------------------------------------------------!
  ! Type that models a spring mass damper system (single mass)
  !-------------------------------------------------------------------!
  
  type, extends(physics) :: smd1

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     real(8) :: m = 1.0d0
     real(8) :: c = 0.02d0
     real(8) :: k = 5.0d0

   contains

     procedure :: assembleResidual
     procedure :: assembleJacobian
     procedure :: getInitialStates

  end type smd1

  !-------------------------------------------------------------------!
  ! Type that models a spring mass damper system (two mass)
  !-------------------------------------------------------------------!
  
  type, extends(physics) :: smd2

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

   contains

     procedure :: assembleResidual => assembleResidual2
     procedure :: assembleJacobian => assembleJacobian2
     procedure :: getInitialStates => getInitialStates2

  end type smd2

contains

  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! -------------------------------------------------------------------!

  subroutine assembleResidual(this, res, time, u, udot, uddot)

    class(smd1) :: this
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

    class(smd1) :: this
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

    class(smd1) :: this

    real(8), intent(in) :: time
    real(8), intent(inout), dimension(:) :: u, udot

    u(1) = 1.0d0

  end subroutine getInitialStates

  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! -------------------------------------------------------------------!

  subroutine assembleResidual2(this, res, time, u, udot, uddot)

    class(smd2) :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in) :: time
    real(8), intent(in), dimension(:) :: u, udot, uddot

    res(1) = uddot(1) + 0.02d0*udot(1)*udot(2) + 5.0d0*u(1)
    res(2) = uddot(2) - 0.05d0*udot(2)*udot(1) + 1.0d0*u(2)*u(1)

  end subroutine assembleResidual2

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

  subroutine assembleJacobian2( this, jac, alpha, beta, gamma, &
       & time, u, udot, uddot )

    class(smd2) :: this
    real(8), intent(inout), dimension(:,:) :: jac
    real(8), intent(in) :: alpha, beta, gamma
    real(8), intent(in) :: time
    real(8), intent(in), dimension(:) :: u, udot, uddot

    ! Zero all  entries first
    jac = 0.0d0

    !-----------------------------------------------------------------!
    ! Get dR/dQ
    !-----------------------------------------------------------------!

    ! derivative of first equation     
    JAC(1,1) = JAC(1,1) + alpha*5.0d0
    JAC(1,2) = JAC(1,2) + alpha*0.0d0

    ! derivative of second equation     
    JAC(2,1) = JAC(2,1) + alpha*1.0d0*u(2)
    JAC(2,2) = JAC(2,2) + alpha*1.0d0*u(1)

    !-----------------------------------------------------------------!
    ! Get dR/dQDOT
    !-----------------------------------------------------------------!

    ! derivative of first equation        
    jac(1,1) = jac(1,1) + beta*0.02d0*udot(2)
    jac(1,2) = jac(1,2) + beta*0.02d0*udot(1)

    ! derivative of second equation
    jac(2,1) = jac(2,1) - beta*0.05d0*udot(2)
    jac(2,2) = jac(2,2) - beta*0.05d0*udot(1)

    !-----------------------------------------------------------------!
    ! Get dR/dQDDOT
    !-----------------------------------------------------------------!

    ! derivative of first equation
    JAC(1,1) = JAC(1,1) + gamma*1.0d0
    JAC(1,2) = JAC(1,2) + gamma*0.0d0

    ! derivative of second equation
    JAC(2,1) = JAC(2,1) + gamma*0.0d0
    JAC(2,2) = JAC(2,2) + gamma*1.0d0

  end subroutine assembleJacobian2

  !---------------------------------------------------------------------!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and udot
  !---------------------------------------------------------------------!

  subroutine getInitialStates2(this, time, u, udot)

    class(smd2) :: this

    real(8), intent(in) :: time
    real(8), intent(inout), dimension(:) :: u, udot

    u(1) = 1.0d0
    u(2) = 2.0d0

  end subroutine getInitialStates2
  
end module spring_mass_damper_class

