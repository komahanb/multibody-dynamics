#include "scalar.fpp"
!=====================================================================!
! Module that contains functions related to spring mass damper system
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module pendulum_functions_class

  use function_class, only  : abstract_function

  implicit none

  private
  public :: time_period

  !-------------------------------------------------------------------!
  ! Type that models any physical phenomenon
  !-------------------------------------------------------------------!
  
  type, extends(abstract_function) :: time_period

     type(scalar) :: m = 1.0d0
     type(scalar) :: g = 9.81d0
     type(scalar) :: pi = 3.1415926535897932384626433832d0
     
   contains  

     procedure :: getFunctionValue
     procedure :: addFuncDVSens
     procedure :: addFuncSVSens 
     procedure :: addDFdU
     procedure :: addDFdUDot
     procedure :: addDFdUDDot

  end type time_period

contains

  !-------------------------------------------------------------------!
  ! Evaluate the kinetic energy for the supplied state and design
  ! variables KE = 0.5 m v^2
  ! -------------------------------------------------------------------!
  
  pure subroutine getFunctionValue(this, f, time, x, u, udot, uddot)
    
    class(time_period), intent(inout) :: this
    type(scalar), intent(inout) :: f
    type(scalar), intent(in), dimension(:) :: x, u, udot, uddot
    real(dp), intent(in) :: time

    f = sqrt(this%g/x(1))/(2.0d0*this %pi)
    
  end subroutine getFunctionValue

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dX
  ! -------------------------------------------------------------------!

  subroutine addFuncDVSens(this, res, scale, time, x, u, udot, uddot)

    class(time_period) :: this
    type(scalar), intent(inout), dimension(:) :: res
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    real(dp), intent(in) :: time
    type(scalar) :: scale

    res(1) = res(1) - scale*sqrt(this%g/(x(1)*x(1)*x(1)))/(4.0d0*this%pi)

  end subroutine addFuncDVSens
  
  !-------------------------------------------------------------------!
  ! Evaluate alpha dF/dU + beta dF/dUDot + gamma dF/dUDDOT
  ! -------------------------------------------------------------------!

  subroutine addFuncSVSens(this, res, alpha, beta, gamma, &
       & time, x, u, udot, uddot)

    class(time_period) :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in) :: time
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    type(scalar), intent(in) :: alpha, beta, gamma

  end subroutine addFuncSVSens

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dU
  ! -------------------------------------------------------------------!

  subroutine addDfdU(this, res, scale, time, x, u, udot, uddot)

    class(time_period) :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in) :: time
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    type(scalar) :: scale

    stop
    
  end subroutine addDfdU

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dUdot
  ! -------------------------------------------------------------------!

  subroutine addDfdUDot(this, res, scale, time, x, u, udot, uddot)

    class(time_period) :: this
    type(scalar), intent(inout), dimension(:) :: res
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    real(dp), intent(in) :: time
    type(scalar) :: scale

    stop
    
  end subroutine addDFdUDot
  
  !-------------------------------------------------------------------!
  ! Evaluate  dF/dUDDot
  ! -------------------------------------------------------------------!

  subroutine addDfdUDDot(this, res, scale, time, x, u, udot, uddot)

    class(time_period) :: this
    type(scalar), intent(inout), dimension(:) :: res
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    real(dp), intent(in)                      :: time
    type(scalar) :: scale

    stop
    
  end subroutine addDfdUDDot

end module pendulum_functions_class
