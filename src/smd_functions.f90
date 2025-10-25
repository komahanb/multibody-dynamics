#include "scalar.fpp"
!=====================================================================!
! Module that contains functions related to spring mass damper system
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module smd_functions_class

  use function_class, only  : abstract_function

  implicit none

  private
  public :: kinetic_energy

  !-------------------------------------------------------------------!
  ! Type that models any physical phenomenon
  !-------------------------------------------------------------------!
  
  type, extends(abstract_function) :: kinetic_energy

   contains  

     procedure :: getFunctionValue
     procedure :: addFuncDVSens
     procedure :: addFuncSVSens 
     procedure :: addDFdU
     procedure :: addDFdUDot
     procedure :: addDFdUDDot

  end type kinetic_energy

contains

  !-------------------------------------------------------------------!
  ! Evaluate the kinetic energy for the supplied state and design
  ! variables KE = 0.5 m v^2
  ! -------------------------------------------------------------------!
  
  pure subroutine getFunctionValue(this, f, time, x, u, udot, uddot)
    
    class(kinetic_energy), intent(inout)   :: this
    type(scalar), intent(inout)            :: f
    type(scalar), intent(in), dimension(:) :: x, u, udot, uddot
    real(dp), intent(in)                   :: time

    f = 0.5d0*x(1)*u(1)**2

!!$    if (size(x) .eq. 3) then
!!$       f = 0.5d0*x(1)*uddot(1)**2 + 0.5d0*x(2)*udot(1)**2 + 0.5d0*x(3)*u(1)**2
!!$    else if (size(x) .eq. 2) then
!!$       f = 0.0d0
!!$    else if (size(x) .eq. 1) then
!!$       f = 0.5d0*x(1)*u(1)**2
!!$       !f = uddot(1)**2 + udot(1)**2 + u(1)**2 + x(1)**2 ! Optimal control problem
!!$    end if
    !f = 0.5d0*x(3)*u(1)**2
    
  end subroutine getFunctionValue

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dX
  ! -------------------------------------------------------------------!

  subroutine addFuncDVSens(this, res, scale, time, x, u, udot, uddot)

    class(kinetic_energy)                     :: this
    type(scalar), intent(inout), dimension(:) :: res
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    real(dp), intent(in)                      :: time
    type(scalar)                              :: scale

!!$    if (size(x).eq.3) then
!!$       res(1) = res(1) + scale*0.5d0*uddot(1)**2 ! wrt to m
!!$       res(2) = res(2) + scale*0.5d0*udot(1)**2 ! wrt to c
!!$       res(3) = res(3) + scale*0.5d0*u(1)**2 ! wrt to k
!!$    else if (size(x) .eq. 1) then
    res(1) = res(1) + scale* 0.5d0*u(1)**2
!!$    end if

  end subroutine addFuncDVSens
  
  !-------------------------------------------------------------------!
  ! Evaluate alpha dF/dU + beta dF/dUDot + gamma dF/dUDDOT
  ! -------------------------------------------------------------------!

  subroutine addFuncSVSens(this, res, alpha, beta, gamma, &
       & time, x, u, udot, uddot)

    class(kinetic_energy)                     :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in)                      :: time
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    type(scalar), intent(in)                  :: alpha, beta, gamma

    res(1) = res(1) + alpha*x(1)*u(1)
!!$
!!$    if (size(x) .eq. 3) then
!!$       res(1) = res(1) + alpha*x(3)*u(1) + beta*x(2)*udot(1) + gamma*x(1)*uddot(1)
!!$    else if (size(x) .eq. 1) then
!!$       !res(1) = res(1) + alpha*2.0d0*u(1) + beta*2.0d0*udot(1) + gamma*2.0d0*uddot(1)
!!$       res(1) = res(1) +  alpha*x(1)*u(1) 
!!$    end  if
    
  end subroutine addFuncSVSens

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dU
  ! -------------------------------------------------------------------!

  subroutine addDfdU(this, res, scale, time, x, u, udot, uddot)

    class(kinetic_energy)                     :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in)                      :: time
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    type(scalar)                              :: scale

    stop

    if (size(x) .eq. 3) then
       res(1) = res(1) + scale*x(3)*u(1) ! to u(1)
    else if (size(x) .eq. 1) then
       res(1) = res(1) + scale*x(1)*u(1)
    end if
   
  end subroutine addDfdU

  !-------------------------------------------------------------------!
  ! Evaluate  dF/dUdot
  ! -------------------------------------------------------------------!

  subroutine addDfdUDot(this, res, scale, time, x, u, udot, uddot)

    class(kinetic_energy)                     :: this
    type(scalar), intent(inout), dimension(:) :: res
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    real(dp), intent(in)                      :: time
    type(scalar)                              :: scale
!!$
!!$    if (size(x) .eq. 3) then
!!$       res(1) = res(1) + scale*x(2)*udot(1) ! wrt to udot(1)
!!$    else if (size(x) .eq. 1) then
!!$       !res(1) = res(1) + scale*2.0d0*udot(1)
!!$    end if
         
  end subroutine addDFdUDot
  
  !-------------------------------------------------------------------!
  ! Evaluate  dF/dUDDot
  ! -------------------------------------------------------------------!

  subroutine addDfdUDDot(this, res, scale, time, x, u, udot, uddot)

    class(kinetic_energy)                     :: this
    type(scalar), intent(inout), dimension(:) :: res
    type(scalar), intent(in), dimension(:)    :: x, u, udot, uddot
    real(dp), intent(in)                      :: time
    type(scalar)                              :: scale
    
    if (size(x) .eq. 3) then
       res(1) = res(1) + scale*x(1)*uddot(1) ! wrt to uddot(1)
    else if (size(x) .eq. 1) then
       !res(1) = res(1) + scale*2.0d0*uddot(1)
    end if

  end subroutine addDfdUDDot

end module smd_functions_class
