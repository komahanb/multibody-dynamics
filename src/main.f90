!=====================================================================!
! Main Program for the landing event simulation
!=====================================================================!

program main

  ! Import essentials
  use iso_fortran_env, only: dp => real64

  ! Import Integrators
  use runge_kutta_integrator, only : DIRK
  use bdf_integrator, only : BDF
  
  ! Import Physics
  use rigid_body_class, only: rigid_body
  use multibody_dynamics_class, only: multibody_dynamics
  use myode_class, only : ODE

  implicit none

  type(DIRK)                       :: dirkobj  ! DIRK Integrator object
  type(BDF)                        :: bdfobj   ! BDF Integrator object

  type(ODE), target                :: smd      ! Spring-mass-damper test ODE
  type(multibody_dynamics), target :: freefall ! Rigid body dynamics system
 
  call dirkobj % setPhysicalSystem(smd)
  call dirkobj % initialize(tfinal = 25.0d0, num_stages=1, &
       & h=0.01d0, nsvars=1, second_order=.true.) ! all are optional except nvars
  call dirkobj % setPrintLevel(0)
  call dirkobj % integrate()
  !  call dirkobj % integrateBackward()
  call dirkobj % writeSolution('smd-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % setPhysicalSystem(smd)
  call bdfobj % initialize(tfinal = 25.0d0, max_bdf_order=3, &
       & h=0.01d0, nsvars=1, second_order=.true.) ! all are optional except nvars
  call bdfobj % setPrintLevel(0)
  call bdfobj % integrate()
  ! call bdfobj % integrateBackward()
  call bdfobj % writeSolution('smd-bdf.dat')
  call bdfobj % finalize()


  call dirkobj % setPhysicalSystem(freefall)
  call dirkobj % initialize(tfinal = 25.0d0, num_stages=1, &
       & h=0.01d0, nsvars=1, second_order=.true.) ! all are optional except nvars
  call dirkobj % setPrintLevel(0)
  call dirkobj % integrate()
  call dirkobj % writeSolution('freefall-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % setPhysicalSystem(freefall)
  call bdfobj % initialize(tfinal = 25.0d0, max_bdf_order=3, &
       & h=0.01d0, nsvars=1, second_order=.true.) ! all are optional except nvars
  call bdfobj % setPrintLevel(0)
  call bdfobj % integrate()
  call bdfobj % writeSolution('freefall-bdf.dat')
  call bdfobj % finalize()

end program





