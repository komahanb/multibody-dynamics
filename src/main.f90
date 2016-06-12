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
!!$  use rigid_body_class              , only : rigid_body
!!$  use multibody_dynamics_class      , only : multibody_dynamics
  use spring_mass_damper_class      , only : smd1, smd2
!!$  use vanderpol_class               , only : vanderpol
!!$  use aero_elastic_oscillator_class , only : aero_elastic_oscillator

  ! Import functions for derivative calculation
  use smd_functions_class           , only : kinetic_energy
  implicit none

  ! Declare Integrators
  type(DIRK)                            :: dirkobj    ! DIRK Integrator object
  type(BDF)                             :: bdfobj     ! BDF Integrator object
  
  ! Declare Physics for testing
  type(smd1)                   , target :: smd1obj    ! Spring-mass-damper test ODE (1 var)
  type(smd2)                   , target :: smd2obj    ! Spring-mass-damper test ODE (2 var)
!!$  type(vanderpol)              , target :: vpl        ! Vanderpol equation (2 var)
!!$  type(multibody_dynamics)     , target :: freefall   ! Rigid body dynamics system (12 vars)
!!$  type(aero_elastic_oscillator), target :: oscillator ! Aeroelastic oscillator (2 vars)

  ! Declare functions that are used
  type(kinetic_energy)         , target :: KE

  ! Design variable array
  real(dp), dimension(:), allocatable :: x, dfdx, dfdxtmp
  real(dp)                            :: fval, ftmp, dh = 1.0d-6
  
  !-------------------------------------------------------------------!
  !                 Spring Mass Damper system                         !
  !-------------------------------------------------------------------!
  
  allocate(X(3), dfdx(3), dfdxtmp(3))

  x(1) = 1.00d0 + dh    ! mass
  x(2) = 0.02d0    ! damping coeff
  x(3) = 5.00d0    ! stiffness coef

  ! Initialize the system
  call smd1obj % initialize(num_state_vars = 1, num_design_vars = 3)

  bdfobj =  BDF(system = smd1obj, tfinal = 1.0d0, h=0.01d0) 

  call bdfobj % evalFuncGrad(num_func=1, func = KE,  num_dv = 3, x = x, &
       & fvals = fval, dfdx= dfdx)

  call bdfobj % evalFDFuncGrad(num_func=1, func = KE,  num_dv = 3, x = x, &
       & fvals = fval, dfdx= dfdxtmp, dh=1.0d-6)

  call bdfobj % finalize()

  print*, "fval         =", fval
  print*, "Adjoint dfdx =", dfdx
  print*, "FD      dfdx =", dfdxtmp

  ! Finalize the system
  call smd1obj % finalize()

  deallocate(X, dfdx, dfdxtmp)

end program
