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
     procedure :: addDFdX          ! partial derivative
     procedure :: addDFdU          ! partial derivative
     procedure :: addDFdUDot       ! partial derivative
     procedure :: addDFdUDDot      ! partial derivative

  end type kinetic_energy

contains

  !-------------------------------------------------------------------!
  ! Evaluate the kinetic energy for the supplied state and design
  ! variables KE = 0.5 m v^2
  ! -------------------------------------------------------------------!
  
  subroutine getFunctionValue(this, res, time, x, u, udot, uddot)
    
    class(kinetic_energy)                :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in)                  :: time
    real(8), intent(in), dimension(:)    :: x, u, udot, uddot

    res = 0.5d0*x(1)*udot(1)**2
    
  end subroutine getFunctionValue

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dX
  ! -------------------------------------------------------------------!

  subroutine addDfdX(this, res, scale, time, x, u, udot, uddot)

    class(kinetic_energy)                :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in)                  :: time
    real(8), intent(in), dimension(:)    :: x, u, udot, uddot
    real(8)                              :: scale

    res(1) = res(1) + scale*0.5d0*udot(1)**2 ! wrt to m
    res(2) = res(2) + scale*0.0d0            ! wrt to c
    res(3) = res(3) + scale*0.0d0            ! wrt to k

  end subroutine addDfdX
  
  !-------------------------------------------------------------------!
  ! Evaluate  dF/dU
  ! -------------------------------------------------------------------!

  subroutine addDfdU(this, res, scale, time, x, u, udot, uddot)

    class(kinetic_energy)                :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in)                  :: time
    real(8), intent(in), dimension(:)    :: x, u, udot, uddot
    real(8)                              :: scale

    res(1) = res(1) + scale*0.0d0            ! wrt to u(1)
    
  end subroutine addDfdU

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dUdot
  ! -------------------------------------------------------------------!

  subroutine addDfdUDot(this, res, scale, time, x, u, udot, uddot)

    class(kinetic_energy)                :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in)                  :: time
    real(8), intent(in), dimension(:)    :: x, u, udot, uddot
    real(8)                              :: scale

    res(1) = res(1) + scale*x(1)*udot(1) ! wrt to udot(1)

  end subroutine addDFdUDot
  
  !-------------------------------------------------------------------!
  ! Evaluate  dF/dUDDot
  ! -------------------------------------------------------------------!

  subroutine addDfdUDDot(this, res, scale, time, x, u, udot, uddot)

    class(kinetic_energy)                :: this
    real(8), intent(inout), dimension(:) :: res
    real(8), intent(in)                  :: time
    real(8), intent(in), dimension(:)    :: x, u, udot, uddot
    real(8)                              :: scale

    res(1) = res(1) + scale*0.0d0            ! wrt to uddot(1)

  end subroutine addDfdUDDot

end module smd_functions_class
