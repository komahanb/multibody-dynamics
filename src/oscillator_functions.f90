#include "scalar.fpp"
!=====================================================================!
! Module that contains functions related to spring mass damper system
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module oscillator_functions_class

  use function_class, only  : abstract_function

  implicit none

  private
  public :: pitch

  interface sign
     module procedure sign_cc
     module procedure sign_cr
     module procedure sign_rc
  end interface

  !-------------------------------------------------------------------!
  ! Type that models any physical phenomenon
  !-------------------------------------------------------------------!
  
  type, extends(abstract_function) :: pitch

   contains  

     procedure :: getFunctionValue ! function value at t, X, U, Udot, Uddot
     procedure :: addFuncDVSens          ! partial derivative
     procedure :: addDFdU          ! partial derivative
     procedure :: addDFdUDot       ! partial derivative
     procedure :: addDFdUDDot      ! partial derivative
     procedure :: addFuncSVSens

  end type pitch

contains

  !-------------------------------------------------------------------!
  ! Evaluate the pitch of the airfoil
  ! -------------------------------------------------------------------!
  
  pure subroutine getFunctionValue(this, f, time, x, u, udot, uddot)
    
    class(pitch), intent(inout)       :: this
    type(scalar), intent(inout)            :: f
    real(dp), intent(in)                  :: time
    type(scalar), intent(in), dimension(:) :: x, u, udot, uddot
    
    f = abs(u(2))
    
  end subroutine getFunctionValue

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dX
  !-------------------------------------------------------------------!

  subroutine addFuncDVSens(this, res, scale, time, x, u, udot, uddot)

    class(pitch)                         :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in)                  :: time
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    type(scalar)                              :: scale
    
 !  res = 0.0d0

  end subroutine addFuncDVSens

    !-------------------------------------------------------------------!
  ! Evaluate alpha dF/dU + beta dF/dUDot + gamma dF/dUDDOT
  ! -------------------------------------------------------------------!

  subroutine addFuncSVSens(this, res, alpha, beta, gamma, &
       & time, x, u, udot, uddot)

    class(pitch)                     :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in)                      :: time
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    type(scalar), intent(in)                  :: alpha, beta, gamma

    stop
    
  end subroutine addFuncSVSens

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dU
  ! -------------------------------------------------------------------!

  subroutine addDfdU(this, res, scale, time, x, u, udot, uddot)

    class(pitch)                :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in)                  :: time
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    type(scalar)                              :: scale
    type(scalar) :: ONE = 1.0d0
!    res(1) = res(1) + scale*0.0d0
    res(2) = res(2) + scale*sign(ONE,u(2))
    
  end subroutine addDfdU

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dUdot
  ! -------------------------------------------------------------------!

  subroutine addDfdUDot(this, res, scale, time, x, u, udot, uddot)

    class(pitch)                :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in)                  :: time
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    type(scalar)                              :: scale
    
    res(1) = res(1) + scale*0.0d0 ! wrt to udot(1)

  end subroutine addDFdUDot
  
  !-------------------------------------------------------------------!
  ! Evaluate  dF/dUDDot
  ! -------------------------------------------------------------------!

  subroutine addDfdUDDot(this, res, scale, time, x, u, udot, uddot)

    class(pitch)                :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in)                  :: time
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    type(scalar)                              :: scale

  !  res(1) = res(1) + scale*0.0d0 ! wrt to uddot(1)

  end subroutine addDfdUDDot

! SIGN, intrinsic, assume that val1 is always a complex*16
!                  in reality could be int
  complex*16 function sign_cc(val1, val2)
    complex*16, intent(in) :: val1, val2
    real*8  sign
    if (real(val2) < 0.) then
      sign = -1.
    else
      sign = 1.
    endif
    sign_cc = sign * val1
    return
  end function sign_cc
  complex*16 function sign_cr(val1, val2)
    complex*16, intent(in) :: val1
    real*8, intent(in) :: val2
    real*8 sign
    if (real(val2) < 0.) then
      sign = -1.
    else
      sign = 1.
    endif
    sign_cr = sign * val1
    return
  end function sign_cr
  complex*16 function sign_rc(val1, val2)
    real*8, intent(in) :: val1
    complex*16, intent(in) :: val2
    real*8 sign
    if (real(val2) < 0.) then
      sign = -1.
    else
      sign = 1.
    endif
    sign_rc = sign * val1
    return
  end function sign_rc

end module oscillator_functions_class
