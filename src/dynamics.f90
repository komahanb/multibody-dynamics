!=====================================================================!
! Module that contains all the implementation of rigid body dynamics
! related physics. 
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module dynamics

  use utils

  implicit none

  private

  public :: rigid_body

  !-------------------------------------------------------------------!
  ! Rigid body type that contains pertaining data and routines that
  ! operate over the type variables
  !-------------------------------------------------------------------!

  type :: rigid_body

   contains

     procedure :: get_residual
     procedure :: get_jacobian

  end type rigid_body

contains

  !-------------------------------------------------------------------!
  ! Residual of the kinematic and dynamic equations in state-space
  ! representation
  !-------------------------------------------------------------------!

  subroutine get_residual(this)

    class(rigid_body) :: this

  end subroutine get_residual

  !-------------------------------------------------------------------!
  ! Jacobian of the kinematic and dynamic equations in state-space
  ! representation
  !-------------------------------------------------------------------!

  subroutine get_jacobian(this)

    class(rigid_body) :: this

  end subroutine get_jacobian

end module dynamics
