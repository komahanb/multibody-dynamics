!=====================================================================!
! Module that contains common procedures for any physical system
! subject to governing equations
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module physics_class
  
  use function_class, only : abstract_function

  implicit none

  private
  public :: physics
  
  !-------------------------------------------------------------------!
  ! Type that models any physical phenomenon
  !-------------------------------------------------------------------!

  type, abstract :: physics
     
     class(abstract_function), pointer :: function

   contains

     procedure(residual_assembly_interface), deferred :: assembleResidual
     procedure(jacobian_assembly_interface), deferred :: assembleJacobian
     procedure(initial_condition_interface), deferred :: getInitialStates

  end type physics

  interface

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
     ! Interface for supplying the initial condition to the integrator
     !----------------------------------------------------------------!
     
     subroutine initial_condition_interface(this, time, u, udot)

       import physics

       class(physics) :: this
       real(8), intent(in) :: time
       real(8), intent(inout), dimension(:) :: u, udot

     end subroutine initial_condition_interface

  end interface

end module physics_class
