!=====================================================================!
! Program to test rigid body dynamics implemetation
!=====================================================================!

program test

  use utils
  use rigid_body_class, only: rigid_body
  use rigid_body_dynamics_class, only: rigid_body_dynamics
  use runge_kutta_integrator, only : DIRK
  
  type(rigid_body) :: body
  type(DIRK) :: integrator
  type(rigid_body_dynamics), target :: falcon
  real(dp) :: q(12)=0.0d0, qdot(12)=0.0d0, qddot(12)=0.0d0
  
!!$
!!$  v % x = 2.0d0
!!$  w % x = 3.0d0
!!$
!!$  A = matrix((/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0 /))
!!$
!!$  ! Initial position in inertial frame (r)
!!$  q(1:3) = (/ 0.0d0, 0.0d0, 25.0d0 /)
!!$
!!$  ! Initial velocity in inertial frame (v)
!!$  q(7:9) = (/ 1.0d0, 2.0d0, 3.0d0 /)
!!$
!!$  ! Initial angular velocity in inertial frame (omega)
!!$  q(10:12) = (/ -1.0d0, -2.0d0, -3.0d0 /)
!!$  
!!$  call body % set (1.0d0, q, qdot, qddot)
!!$  
!!$  print *, body % get_residual()
  
  ! Set the physics into the integrator
  integrator % system => falcon

  call integrator % initialize(tfinal = 10.0d0, num_stages=3, &
       & h=0.01d0, nvars=12, second_order=.false.)
  call integrator % integrate()
  call integrator % write_solution()
  call integrator % finalize()

end program test


