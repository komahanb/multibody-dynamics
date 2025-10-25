#include "scalar.fpp"
!=====================================================================!
! Module that implements an Aero-Elastic Oscillator reported in:
! Zhao, L. and Yang, Z., ``Chaotic motions of an airfoil with
! non-linear stiffness in incompressible flow,'' Journal of Sound and
! Vibration, Vol. 138, No. 2, 1990, pp. 245–254.
! =====================================================================!

module aero_elastic_oscillator_class
  
  use physics_class,  only : physics
  use function_class, only : abstract_function

  implicit none

  private

  public :: aero_elastic_oscillator

  !-------------------------------------------------------------------!
  ! Type that implements the aero-elastic oscillator (2 state vars)
  !-------------------------------------------------------------------!
  
  type, extends(physics) :: aero_elastic_oscillator

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     type(scalar) :: Q = 8.0d0  ! reduced dynamic pressure
     type(scalar) :: E = 20.0d0 ! non-linear stiffness factor

   contains    

     procedure :: mapDesignVars
     procedure :: assembleResidual
     procedure :: assembleJacobian
     procedure :: getInitialStates
     procedure :: getResidualDVSens

  end type aero_elastic_oscillator

contains
  
  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. 
  ! u(1) is the plunging state of the oscillator
  ! u(2) is the pitching state of the oscillator
  ! -------------------------------------------------------------------!
  
  subroutine assembleResidual( this, res, time, u, udot, uddot )

    class(aero_elastic_oscillator) :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in) :: time
    type(scalar), intent(in), dimension(:) :: u, udot, uddot

    res(1) = uddot(1) + 0.25d0*uddot(2) + 0.1d0*udot(1) &
         & + 0.2d0*u(1) + 0.1d0*this % Q*u(2)

    res(2) = 0.25d0*uddot(1) + 0.5d0*uddot(2) + 0.1d0*udot(2) &
         & + 0.5d0*u(2) + this % E *u(2)*u(2)*u(2) - 0.1d0*this % Q*u(2)

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
  ! with the partial derivatives DRDQ, DRDQdot and DRDQddot
  ! respectively.
  !-------------------------------------------------------------------!

  subroutine assembleJacobian( this, jac, alpha, beta, gamma, &
       & time, u, udot, uddot )

    class(aero_elastic_oscillator) :: this
    type(scalar), intent(inout), dimension(:,:) :: jac
    type(scalar), intent(in) :: alpha, beta, gamma
    real(dp), intent(in) :: time
    type(scalar), intent(in), dimension(:) :: u, udot, uddot
    
    ! Zero all entries first
    jac = 0.0d0

    !-----------------------------------------------------------------!
    ! Add dR/dQ
    !-----------------------------------------------------------------!

    ! derivative of first equation
    
    jac(1,1) = jac(1,1) + alpha*0.2d0
    jac(1,2) = jac(1,2) + alpha*0.1d0 * this % Q

    ! derivative of second equation
    
    jac(2,1) = jac(2,1) + alpha*0.0d0
    jac(2,2) = jac(2,2) + alpha*(0.5d0 + this % E*3.0d0*u(2)*u(2) - 0.1d0 * this % Q)
    
    !-----------------------------------------------------------------!
    ! Add dR/dQDOT
    !-----------------------------------------------------------------!
    
    ! derivative of first equation
    
    jac(1,1) = jac(1,1) + beta*0.1d0
    jac(1,2) = jac(1,2) + beta*0.0d0

    ! derivative of second equation

    jac(2,1) = jac(2,1) + beta*0.0d0
    jac(2,2) = jac(2,2) + beta*0.1d0

    !-----------------------------------------------------------------!
    ! Add dR/dQDDOT
    !-----------------------------------------------------------------!
    
    ! derivative of first equation
    
    jac(1,1) = jac(1,1) + gamma*1.0d0
    jac(1,2) = jac(1,2) + gamma*0.25d0

    ! derivative of second equation

    jac(2,1) = jac(2,1) + gamma*0.25d0
    jac(2,2) = jac(2,2) + gamma*0.50d0

  end subroutine assembleJacobian

  !---------------------------------------------------------------------!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and udot
  !---------------------------------------------------------------------!

  subroutine getInitialStates(this, time, u, udot)

    class(aero_elastic_oscillator) :: this

    real(dp), intent(in) :: time
    type(scalar), intent(inout), dimension(:) :: u, udot
    
    u(1) = 0.00d0
    u(2) = 0.00d0
    
    udot(1) = 0.01d0
    udot(2) = 0.00d0

  end subroutine getInitialStates
 
  !===================================================================!
  ! Return the number of state variables
  !===================================================================!
  
  function getNumStateVars(this)

    class(aero_elastic_oscillator) :: this
    integer          :: getNumStateVars

    getNumStateVars = this % num_state_vars

  end function getNumStateVars

  !-------------------------------------------------------------------!
  ! Routine for evaluating the gradient of Residual with respect
  ! to the design X
  !-------------------------------------------------------------------!
  
  subroutine getResidualDVSens(this, jac, scale, time, x, u, udot, uddot)

    class(aero_elastic_oscillator)         :: this
    type(scalar), intent(inout), dimension(:,:) :: jac
    real(dp), intent(in) :: time
    type(scalar), intent(in), dimension(:)      :: x, u, udot, uddot
    type(scalar)                                :: scale

    jac = 0.0d0
    
    jac(1,1) =   0.1d0*u(2)
    jac(1,2) =   0.0d0

    jac(2,1) = - 0.1d0*u(2)
    jac(2,2) =   u(2)*u(2)*u(2)

  end subroutine getResidualDVSens
  
  !-------------------------------------------------------------------!
  ! Map the the design variables into the class variables
  !-------------------------------------------------------------------!
  
  subroutine mapDesignVars(this)

    class(aero_elastic_oscillator) :: this

    this % Q = this % x(1)
    this % E = this % x(2)

  end subroutine mapDesignVars

end Module aero_elastic_oscillator_class
