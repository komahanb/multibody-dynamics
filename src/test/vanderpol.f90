!=====================================================================!
! Module that provides the ability to the user to implement random
! ODEs and use with the DIRK scheme
!=====================================================================!

Module vanderpol_class

  use physics_class, only : physics
  
  implicit none

  private

  public :: vanderpol

  !-------------------------------------------------------------------!
  ! Type that implements vanderpol equations in first order form
  !-------------------------------------------------------------------!

  type, extends(physics) :: vanderpol

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here
     
   contains

     procedure :: assembleResidual
     procedure :: assembleJacobian
     procedure :: getInitialStates
     
  end type vanderpol

contains
  
  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! -------------------------------------------------------------------!
  
  subroutine assembleResidual( this, res, time, u, udot, uddot )

    class(vanderpol) :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in) :: time
    real(8), intent(in), dimension(:) :: u, udot, uddot

    res(1) = udot(1) - u(2)
    res(2) = udot(2) - ( 1.0d0 - u(1)*u(1) )*u(2) + u(1)

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

    class(vanderpol) :: this
    real(8), intent(inout), dimension(:,:) :: jac
    real(8), intent(in) :: alpha, beta, gamma
    real(8), intent(in) :: time
    real(8), intent(in), dimension(:) :: u, udot, uddot
    
    ! Zero all entries first
    jac = 0.0d0

    !-----------------------------------------------------------------!
    ! Add dR/dQ
    !-----------------------------------------------------------------!

    ! derivative of first equation
    
    jac(1,1) = jac(1,1) + alpha*0.0d0
    jac(1,2) = jac(1,2) - alpha*1.0d0

    ! derivative of second equation
    
    jac(2,1) = jac(2,1) + alpha*(1.0d0 + 2.0d0*u(1)*u(2))
    jac(2,2) = jac(2,2) + alpha*(u(1)*u(1)-1.0d0)
    
    !-----------------------------------------------------------------!
    ! Add dR/dQDOT
    !-----------------------------------------------------------------!
    
    ! derivative of first equation
    
    jac(1,1) = jac(1,1) + beta*1.0d0
    jac(1,2) = jac(1,2) + beta*0.0d0

    ! derivative of second equation

    jac(2,1) = jac(2,1) + beta*0.0d0
    jac(2,2) = jac(2,2) + beta*1.0d0

  end subroutine assembleJacobian

  !---------------------------------------------------------------------!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and udot
  !---------------------------------------------------------------------!

  subroutine getInitialStates(this, time, u, udot)

    class(vanderpol) :: this

    real(8), intent(in) :: time
    real(8), intent(inout), dimension(:) :: u, udot
    
    u(1) = 2.0d0
    u(2) = 0.0d0

  end subroutine getInitialStates

End Module vanderpol_class

