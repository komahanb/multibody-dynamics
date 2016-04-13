!=====================================================================!
! Program to test rigid body dynamics implemetation
!=====================================================================!

program test

  use utils
  use rigid_body_dynamics

  real(dp) :: q(12), qdot(12), qddot(12)
  
  type(rigid_body) :: body

  call body % set (1.0d0, q, qdot, qddot)
  
  print *, body % get_residual()
  
end program test
