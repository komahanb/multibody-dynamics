!=====================================================================!
! Module that contains common procedures for any physical system
! subject to governing equations
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module physics_class

  use iso_fortran_env, only : dp => real64
  use function_class, only  : abstract_function

  implicit none

  private
  public :: physics
  
  !-------------------------------------------------------------------!
  ! Type that models any physical phenomenon
  !-------------------------------------------------------------------!

  type, abstract :: physics
     
     class(abstract_function) , pointer  :: func ! function of interest
     real(dp), dimension(:), allocatable :: x
    
   contains  
     
     procedure(residual_assembly_interface), deferred :: assembleResidual
     procedure(jacobian_assembly_interface), deferred :: assembleJacobian
     procedure(initial_condition_interface), deferred :: getInitialStates

     procedure(InterfaceInitialize), deferred      :: initialize
     procedure(InterfaceSetDesignVars), deferred   :: setDesignVars
     procedure(InterfaceGetNumStateVars), deferred :: getNumStateVars

     procedure :: setFunction
     
  end type physics
  
  interface
     
     !----------------------------------------------------------------!
     ! Interface for initialization tasks
     !----------------------------------------------------------------!
     
     subroutine InterfaceInitialize(this,  x, function)
       import physics
       import abstract_function
       class(physics) :: this
       class(abstract_function), target, OPTIONAL  :: function
       real(8), intent(in), dimension(:), OPTIONAl :: x
     end subroutine InterfaceInitialize

     !----------------------------------------------------------------!
     ! User implementation of how design variables are mapped to the
     ! local paramters
     ! ----------------------------------------------------------------!
     
     subroutine InterfaceSetDesignVars(this, x)
       import physics
       class(physics) :: this
       real(8), intent(in), dimension(:) :: x
     end subroutine InterfaceSetDesignVars
     
     !----------------------------------------------------------------!
     ! Interface for residual assembly at each time step
     !----------------------------------------------------------------!

     subroutine residual_assembly_interface(this, res, time, u, udot, uddot)

       import physics

       class(physics) :: this
       real(8), intent(inout), dimension(:) :: res
       real(8), intent(in) :: time
       real(8), intent(in), dimension(:) :: u, udot, uddot

     end subroutine residual_assembly_interface

     !----------------------------------------------------------------!
     ! Interface for jacobian assembly at each time step
     !----------------------------------------------------------------!

     subroutine jacobian_assembly_interface(this, jac, alpha, beta, gamma, &
          & time, u, udot, uddot)

       import physics

       class(physics) :: this
       real(8), intent(inout), dimension(:,:) :: jac
       real(8), intent(in) :: alpha, beta, gamma
       real(8), intent(in) :: time
       real(8), intent(in), dimension(:) :: u, udot, uddot

     end subroutine jacobian_assembly_interface
     
     !----------------------------------------------------------------!
     ! Interface for supplying the initial condition to the integrator!
     !----------------------------------------------------------------!
     
     subroutine initial_condition_interface(this, time, u, udot)

       import physics

       class(physics) :: this
       real(8), intent(in) :: time
       real(8), intent(inout), dimension(:) :: u, udot

     end subroutine initial_condition_interface

     !----------------------------------------------------------------!
     ! Return the number of state variables
     !----------------------------------------------------------------!
     
     function InterfaceGetNumStateVars(this)
       import physics
       class(physics) :: this
       integer :: InterfaceGetNumStateVars
     end function InterfaceGetNumStateVars

  end interface

contains
  
  !===================================================================!
  ! Set the function created into the system                          !
  !===================================================================!
  
  subroutine setFunction(this, func)
    
    class(physics)                   :: this
    class(abstract_function), target :: func

    this % func => func

  end subroutine setFunction

end module physics_class