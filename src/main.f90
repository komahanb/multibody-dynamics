!=====================================================================!
! Main Program for testing the integrators on different test problems
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
  use spring_mass_damper_class, only: smd1, smd2
  use vanderpol_class, only : vanderpol

  implicit none

  ! Declare Integrators
  type(DIRK)                       :: dirkobj  ! DIRK Integrator object
  type(BDF)                        :: bdfobj   ! BDF Integrator object

  ! Declare Physics for testing
  type(smd1), target               :: smd1obj  ! Spring-mass-damper test ODE
  type(smd2), target               :: smd2obj  ! Spring-mass-damper test ODE
  type(vanderpol), target          :: vpl      ! Vanderpol equation
  type(multibody_dynamics), target :: freefall ! Rigid body dynamics system
 
  !-------------------------------------------------------------------!
  !                 Spring Mass Damper system                         !
  !-------------------------------------------------------------------!
  
  call dirkobj % setPhysicalSystem(smd1obj)
  call dirkobj % initialize(tfinal = 1.0d0, num_stages=1, &
       & h=0.01d0, nsvars=1, second_order=.true.) ! all are optional except nvars
  call dirkobj % setPrintLevel(0)
  call dirkobj % integrate()
  call dirkobj % writeSolution('smd-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % setPhysicalSystem(smd1obj)
  call bdfobj % initialize(tfinal = 1.0d0, max_bdf_order=3, &
       & h=0.01d0, nsvars=1, second_order=.true.) ! all are optional except nvars
  call bdfobj % setPrintLevel(0)
  call bdfobj % integrate()
  call bdfobj % writeSolution('smd-bdf.dat')
  call bdfobj % finalize()

  !-------------------------------------------------------------------!
  !        Spring Mass Damper system (2 var second order)             !
  !-------------------------------------------------------------------!
  
  call dirkobj % setPhysicalSystem(smd2obj)
  call dirkobj % initialize(tfinal = 1.0d0, num_stages=1, &
       & h=0.01d0, nsvars=2, second_order=.true.) ! all are optional except nvars
  call dirkobj % setPrintLevel(0)
  call dirkobj % integrate()
  call dirkobj % writeSolution('smd2-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % setPhysicalSystem(smd2obj)
  call bdfobj % initialize(tfinal = 1.0d0, max_bdf_order=3, &
       & h=0.01d0, nsvars=2, second_order=.true.) ! all are optional except nvars
  call bdfobj % setPrintLevel(0)
  call bdfobj % integrate()
  call bdfobj % writeSolution('smd2-bdf.dat')
  call bdfobj % finalize()
  
  !-------------------------------------------------------------------!
  !                 Vanderpol Equation ( 3 variables)
  !-------------------------------------------------------------------!
  
  call dirkobj % setPhysicalSystem(vpl)
  call dirkobj % setApproximateJacobian(.false.)
  call dirkobj % initialize(tfinal = 20.0d0, num_stages=3, &
       & h=0.01d0, nsvars=2, second_order=.false.) ! all are optional except nvars
  call dirkobj % setPrintLevel(0)
  call dirkobj % integrate()
  call dirkobj % writeSolution('vpl-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % setPhysicalSystem(vpl)
  call bdfobj % setApproximateJacobian(.false.)
  call bdfobj % initialize(tfinal = 20.0d0, max_bdf_order=3, &
       & h=0.1d0, nsvars=2, second_order=.true.) ! all are optional except nvars
  call bdfobj % setPrintLevel(0)
  call bdfobj % integrate()
  call bdfobj % writeSolution('vpl-bdf.dat')
  call bdfobj % finalize()

  !-------------------------------------------------------------------!
  !                 Rigidbody Dynamics (12 variables)                 !
  !-------------------------------------------------------------------!

  call dirkobj % setPhysicalSystem(freefall)
  call dirkobj % setApproximateJacobian(.true.)
  call dirkobj % initialize(tinit = 0.0d0, tfinal = 25.0d0, num_stages=3, &
       & h=0.01d0, nsvars=12, second_order=.true.) ! all are optional except nvars
  call dirkobj % setPrintLevel(0)
  call dirkobj % integrate()
  call dirkobj % writeSolution('freefall-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % setPhysicalSystem(freefall)
  call bdfobj % setApproximateJacobian(.true.)
  call bdfobj % initialize(tfinal = 5.0d0, max_bdf_order=3, &
       & h=0.01d0, nsvars=12, second_order=.true.) ! all are optional except nvars
  call bdfobj % setPrintLevel(0)
  call bdfobj % integrate()
  call bdfobj % writeSolution('freefall-bdf.dat')
  call bdfobj % finalize()

end program main
