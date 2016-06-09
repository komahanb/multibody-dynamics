!=====================================================================!
! Module that contains common procedures for any function of interest
! that the user wishes to implement
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module function_class
  
  !-------------------------------------------------------------------!
  ! Type that models any funtion of interest e.g. kinetic energy, lift
  !-------------------------------------------------------------------!
  
  type, abstract :: function

   contains

     procedure(interface_evaluate), deferred :: getFunctionValue ! function value at t, X, U, Udot, Uddot
     procedure(interface_evaluate), deferred :: getdFdX          ! partial derivative
     procedure(interface_evaluate), deferred :: getdFdU          ! partial derivative
     procedure(interface_evaluate), deferred :: getdFdUDot       ! partial derivative
     procedure(interface_evaluate), deferred :: getdFdUDDot      ! partial derivative

  end type function

  interface

     !----------------------------------------------------------------!
     ! Interface for evaluting the function for t, X, U, Udot, Uddot
     !----------------------------------------------------------------!

     subroutine interface_evaluate(this, res, time, x, u, udot, uddot)

       import function

       class(function) :: this
       real(8), intent(inout), dimension(:) :: res
       real(8), intent(in) :: time
       real(8), intent(in), dimension(:) :: x, u, udot, uddot
       
     end subroutine interface_evaluate

  end interface

end module physics_class
