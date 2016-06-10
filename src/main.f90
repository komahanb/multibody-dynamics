!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

program main

  ! Import essentials  
  use iso_fortran_env               , only : dp => real64

  ! Import Integrators
  use runge_kutta_integrator        , only : DIRK
  use bdf_integrator                , only : BDF
  
  ! Import Physics
  use rigid_body_class              , only : rigid_body
  use multibody_dynamics_class      , only : multibody_dynamics
  use spring_mass_damper_class      , only : smd1, smd2
  use vanderpol_class               , only : vanderpol
  use aero_elastic_oscillator_class , only : aero_elastic_oscillator

  ! Import functions for derivative calculation
  use smd_functions_class           , only : kinetic_energy
  implicit none

  ! Declare Integrators
  type(DIRK)                            :: dirkobj    ! DIRK Integrator object
  type(BDF)                             :: bdfobj     ! BDF Integrator object
  
  ! Declare Physics for testing
  type(smd1)                   , target :: smd1obj    ! Spring-mass-damper test ODE (1 var)
  type(smd2)                   , target :: smd2obj    ! Spring-mass-damper test ODE (2 var)
  type(vanderpol)              , target :: vpl        ! Vanderpol equation (2 var)
  type(multibody_dynamics)     , target :: freefall   ! Rigid body dynamics system (12 vars)
  type(aero_elastic_oscillator), target :: oscillator ! Aeroelastic oscillator (2 vars)

  ! Declare functions that are used
  type(kinetic_energy)         , target :: KE

  ! Design variable array
  real(dp), dimension(:), allocatable   :: X
  
  !-------------------------------------------------------------------!
  !                 Spring Mass Damper system                         !
  !-------------------------------------------------------------------!

  allocate(X(3))

  x(1) = 1.00d0 ! mass
  x(2) = 0.02d0 ! damping coeff
  x(3) = 5.00d0 ! stiffness coef

  ! Set the design variable and the function of interest into the
  ! system object
  call smd1obj % initialize(x, KE)
   
  ! Solve the system from tinit to tfinal using the integrator
  
  call dirkobj % initialize(system = smd1obj, tfinal = 1.0d0, h=0.01d0, second_order=.true., num_stages=1)
  call dirkobj % integrate() 
  call dirkobj % writeSolution('smd-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % initialize(system = smd1obj, tfinal = 1.0d0, h=0.01d0, second_order=.true., max_bdf_order=3)
  call bdfobj % integrate()
  call bdfobj % marchBackwards()
  call bdfobj % writeSolution('smd-bdf.dat')
  call bdfobj % finalize()
  
  deallocate(X)

  !-------------------------------------------------------------------!
  !        Spring Mass Damper system (2 var second order)             !
  !-------------------------------------------------------------------!

  call smd2obj % initialize()

  call dirkobj % initialize(system = smd2obj, tfinal = 1.0d0, h=0.01d0, second_order=.true., num_stages=1)
  call dirkobj % setPrintLevel(0)
  call dirkobj % integrate()
  call dirkobj % writeSolution('smd2-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % initialize(system = smd2obj, tfinal = 1.0d0, h=0.01d0, second_order=.true., max_bdf_order=2)
  call bdfobj % setPrintLevel(0)
  call bdfobj % integrate()
  call bdfobj % writeSolution('smd2-bdf.dat')
  call bdfobj % finalize()

  !-------------------------------------------------------------------!
  !                 Vanderpol Equation ( 3 variables)
  !-------------------------------------------------------------------!
  
  call vpl % initialize()

  call dirkobj % initialize(system = vpl, tfinal = 20.0d0, h=0.01d0, second_order=.true., num_stages=3)
  call dirkobj % integrate()
  call dirkobj % writeSolution('vpl-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % initialize(system = vpl, tfinal = 20.0d0, h=0.01d0, second_order=.true., max_bdf_order=3)
  call bdfobj % integrate()
  call bdfobj % writeSolution('vpl-bdf.dat')
  call bdfobj % finalize()

  !-------------------------------------------------------------------!
  !                 Rigidbody Dynamics (12 variables)                 !
  !-------------------------------------------------------------------!

  call freefall % initialize()

  call dirkobj % initialize(system = freefall, tfinal = 5.0d0, h=0.01d0, second_order=.true., num_stages=3)
  call dirkobj % setApproximateJacobian(.true.)
  call dirkobj % integrate()
  call dirkobj % writeSolution('freefall-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % initialize(system = freefall, tfinal = 5.0d0, h=0.01d0, second_order=.true., max_bdf_order=3)
  call bdfobj % setApproximateJacobian(.true.)
  call bdfobj % integrate()
  call bdfobj % writeSolution('freefall-bdf.dat')
  call bdfobj % finalize()

  !-------------------------------------------------------------------!
  !                 Aeroelastic Oscillator (2 variables)              !
  !-------------------------------------------------------------------!

  call oscillator % initialize()

  call dirkobj % initialize(system = oscillator, tinit = 0.0d0, tfinal = 50.0d0,  h=0.001d0, second_order=.true., num_stages=3)
  call dirkobj % integrate()
  call dirkobj % writeSolution('oscillator-dirk.dat')
  call dirkobj % finalize()

  call bdfobj % initialize(system = oscillator, tfinal = 50.0d0, h=0.001d0, second_order=.true., max_bdf_order=2)
  call bdfobj % integrate()
  call bdfobj % writeSolution('oscillator-bdf.dat')
  call bdfobj % finalize()

end program main
