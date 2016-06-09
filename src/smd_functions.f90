!=====================================================================!
! Module that contains functions related to spring mass damper system
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module smd_functions_class

  use iso_fortran_env, only : dp => real64
  use function_class, only  : abstract_function


  implicit none

  private
  public :: kinetic_energy

  !-------------------------------------------------------------------!
  ! Type that models any physical phenomenon
  !-------------------------------------------------------------------!
  
  type, extends(abstract_function) :: kinetic_energy

   contains  

     procedure :: getFunctionValue ! function value at t, X, U, Udot, Uddot
     procedure :: getdFdX          ! partial derivative
     procedure :: getdFdU          ! partial derivative
     procedure :: getdFdUDot       ! partial derivative
     procedure :: getdFdUDDot      ! partial derivative

  end type kinetic_energy

contains

  !-------------------------------------------------------------------!
  ! Evaluate the kinetic energy for the supplied state and design
  ! variables
  ! -------------------------------------------------------------------!
  
  subroutine getFunctionValue(this, res, time, x, u, udot, uddot)
    
    class(kinetic_energy)                :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in)                  :: time
    real(8), intent(in), dimension(:)    :: x, u, udot, uddot

  end subroutine getFunctionValue

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dX
  ! -------------------------------------------------------------------!

  subroutine getdFdX(this, res, time, x, u, udot, uddot)

    class(kinetic_energy)                :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in)                  :: time
    real(8), intent(in), dimension(:)    :: x, u, udot, uddot

  end subroutine getdFdX

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dU
  ! -------------------------------------------------------------------!

  subroutine getdFdU(this, res, time, x, u, udot, uddot)

    class(kinetic_energy)                :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in)                  :: time
    real(8), intent(in), dimension(:)    :: x, u, udot, uddot

  end subroutine getdFdU

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dUdot
  ! -------------------------------------------------------------------!

  subroutine getdFdUDot(this, res, time, x, u, udot, uddot)

    class(kinetic_energy)                :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in)                  :: time
    real(8), intent(in), dimension(:)    :: x, u, udot, uddot
    
  end subroutine getdFdUDot
  
  !-------------------------------------------------------------------!
  ! Evaluate  dF/dUDDot
  ! -------------------------------------------------------------------!

  subroutine getdFdUDDot(this, res, time, x, u, udot, uddot)

    class(kinetic_energy)                :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in)                  :: time
    real(8), intent(in), dimension(:)    :: x, u, udot, uddot

  end subroutine getdFdUDDot

end module smd_functions_class
