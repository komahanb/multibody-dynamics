!=====================================================================!
! Program to test rigid body dynamics implemetation
!=====================================================================!

program test

  use utils
  use rigid_body_class, only: rigid_body
  use multibody_dynamics_class, only: multibody_dynamics
  use runge_kutta_integrator, only : DIRK
  use myode_class, only : ODE

  type(DIRK) :: integrator

  type(multibody_dynamics), target :: falcon ! has equations from
                                             ! multibody dynamics

  type(ODE), target :: myode    ! has user supplied ODEs,
                                ! demonstrated with a one
                                !  DOF spring mass system
  
  !-------------------------------------------------------------------!
  ! Set the physics into the integrator
  !-------------------------------------------------------------------!

  integrator % system => falcon

  call integrator % initialize(tfinal = 10.0d0, num_stages=3, &
       & h=0.01d0, nvars=12, second_order=.true.)
  call integrator % integrate()
  call integrator % write_solution()
  call integrator % finalize()

  !-------------------------------------------------------------------!
  ! Solve an example ODE using the same DIRK scheme
  !-------------------------------------------------------------------!
  
  integrator % system => myode

  call integrator % initialize(tfinal = 10.0d0, num_stages=2, &
       & h=0.01d0, nvars=1, second_order=.true.)
  call integrator % integrate()
  call integrator % write_solution()
  call integrator % finalize()

end program test


