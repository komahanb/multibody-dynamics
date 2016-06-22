#include "scalar.fpp"
!=====================================================================!
! Module that provides the ability to the user to implement random
! ODEs and use with the DIRK scheme
!=====================================================================!

module spring_mass_damper_class

  use physics_class,  only : physics
  use function_class, only : abstract_function

  implicit none
  
  private

  public :: smd1, smd2
  
  !-------------------------------------------------------------------!
  ! Type that models a spring mass damper system (single mass)
  !-------------------------------------------------------------------!
  
  type, extends(physics) :: smd1

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     type(scalar) :: m = 1.0d0
     type(scalar) :: c = 0.02d0
     type(scalar) :: k = 5.0d0

   contains

     procedure :: mapDesignVars
     procedure :: assembleResidual
     procedure :: assembleJacobian
     procedure :: getInitialStates
     procedure :: getResidualDVSens

  end type smd1

  !-------------------------------------------------------------------!
  ! Type that models a spring mass damper system (two mass)
  !-------------------------------------------------------------------!
  
  type, extends(physics) :: smd2

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

   contains

     procedure :: mapDesignVars    => mapDesignVars2
     procedure :: assembleResidual => assembleResidual2
     procedure :: assembleJacobian => assembleJacobian2
     procedure :: getInitialStates => getInitialStates2
     procedure :: getResidualDVSens => getResidualDVSens2
     
  end type smd2

contains
  
  !-------------------------------------------------------------------!
  ! Map the the design variables into the class variables
  !-------------------------------------------------------------------!
  
  subroutine mapDesignVars(this)

    class(smd1) :: this
        
    this % M = this % x(1)
    this % C = this % x(2)
    this % K = this % x(3)

  end subroutine mapDesignVars

!!$
!!$  !-------------------------------------------------------------------!
!!$  ! Constructor for the spring mass damper system 
!!$  !-------------------------------------------------------------------!
!!$  
!!$  subroutine initialize(this, x, function)
!!$
!!$    class(smd1)                                 :: this
!!$    class(abstract_function), target, OPTIONAL  :: function
!!$    type(scalar), intent(in), dimension(:), OPTIONAl :: x
!!$
!!$    ! Set the number of state variables
!!$    this % num_state_vars = 1
!!$
!!$    if (present(function)) call this % setFunction(function)
!!$
!!$    if (present(x)) then
!!$
!!$       this % num_design_vars = size(x)
!!$
!!$       allocate(this % x(this % num_design_vars))
!!$       this % x = x
!!$       
!!$       call this % setDesignVars(x)
!!$       
!!$    end if
!!$
!!$  end subroutine initialize

!!$  !===================================================================!
!!$  ! Sets the design variables into the system
!!$  !===================================================================!
!!$  
!!$  subroutine setDesignVars(this, num_dvs, x)
!!$
!!$    class(smd1)                        :: this
!!$    type(scalar), intent(in), dimension(:)  :: x
!!$   
!!$    ! Overwrite the values to supplied ones
!!$    if (this % num_design_vars .eq. 1) then 
!!$       this % m = x(1)
!!$    else if (this % num_design_vars .eq. 2) then
!!$       this % m = x(1); this % c = x(2);
!!$    else if (this % num_design_vars .eq. 3) then
!!$       this % m = x(1); this % c = x(2); this % K = x(3);          
!!$    end if
!!$
!!$  end subroutine setDesignVars
 
  !---------------------------------------------------------------------!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and udot
  !---------------------------------------------------------------------!

  subroutine getInitialStates(this, time, u, udot)

    class(smd1) :: this

    type(scalar), intent(in) :: time
    type(scalar), intent(inout), dimension(:) :: u, udot

    u(1) = 1.0d0

  end subroutine getInitialStates

  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! -------------------------------------------------------------------!

  subroutine assembleResidual(this, res, time, u, udot, uddot)

    class(smd1) :: this
    type(scalar), intent(inout), dimension(:) :: res
    type(scalar), intent(in) :: time
    type(scalar), intent(in), dimension(:) :: u, udot, uddot

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
    type(scalar), intent(inout), dimension(:,:) :: jac
    type(scalar), intent(in) :: alpha, beta, gamma
    type(scalar), intent(in) :: time
    type(scalar), intent(in), dimension(:) :: u, udot, uddot

    ! Zero all  entries first
    jac = 0.0d0

    jac(1,1) = gamma * this % m + beta * this % c + alpha * this % k

  end subroutine assembleJacobian

  !----------------------------------------------------------------!
  ! Routine for evaluating the gradient of Residual with respect
  ! to the design X
  !----------------------------------------------------------------!
  
  subroutine getResidualDVSens(this, jac, scale, time, x, u, udot, uddot)

    class(smd1)                            :: this
    type(scalar), intent(inout), dimension(:,:) :: jac
    type(scalar), intent(in)                    :: time
    type(scalar), intent(in), dimension(:)      :: x, u, udot, uddot
    type(scalar)                                :: scale

    jac = 0.0d0 

    jac(1,1) = uddot(1)
    jac(1,2) = udot(1)
    jac(1,3) = u(1)

  end subroutine getResidualDVSens

  !---------------------------------------------------------------------!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and udot
  !---------------------------------------------------------------------!

  subroutine getInitialStates2(this, time, u, udot)

    class(smd2) :: this

    type(scalar), intent(in) :: time
    type(scalar), intent(inout), dimension(:) :: u, udot

    u(1) = 1.0d0
    u(2) = 2.0d0

  end subroutine getInitialStates2

  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! -------------------------------------------------------------------!

  subroutine assembleResidual2(this, res, time, u, udot, uddot)

    class(smd2) :: this
    type(scalar), intent(inout), dimension(:) :: res
    type(scalar), intent(in) :: time
    type(scalar), intent(in), dimension(:) :: u, udot, uddot

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
    type(scalar), intent(inout), dimension(:,:) :: jac
    type(scalar), intent(in) :: alpha, beta, gamma
    type(scalar), intent(in) :: time
    type(scalar), intent(in), dimension(:) :: u, udot, uddot

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
    
  !-------------------------------------------------------------------!
  ! Routine for evaluating the gradient of Residual with respect
  ! to the design X
  !-------------------------------------------------------------------!
  
  subroutine getResidualDVSens2(this, jac, scale, time, x, u, udot, uddot)

    class(smd2)                            :: this
    type(scalar), intent(inout), dimension(:,:) :: jac
    type(scalar), intent(in)                    :: time
    type(scalar), intent(in), dimension(:)      :: x, u, udot, uddot
    type(scalar)                                :: scale
    
    stop"Not implemented"

  end subroutine getResidualDVSens2

  
  !-------------------------------------------------------------------!
  ! Map the the design variables into the class variables
  !-------------------------------------------------------------------!
  
  subroutine mapDesignVars2(this)

    class(smd2) :: this
        
    print *, "dummy impl"

  end subroutine mapDesignVars2

end module spring_mass_damper_class

