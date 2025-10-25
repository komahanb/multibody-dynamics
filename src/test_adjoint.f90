!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

#include "scalar.fpp"

program main

  use utils                         , only : real_part
  ! Import Integrators
  use runge_kutta_integrator        , only : DIRK
  use bdf_integrator                , only : BDF
  use abm_integrator                , only : ABM
  use nbg_integrator                , only : NBG

  ! Import physics
  use spring_mass_damper_class      , only : smd1

  ! Import functions
  use smd_functions_class           , only : kinetic_energy

  implicit none

  ! Declare integrators
  type(DIRK)                            :: dirkobj    ! DIRK Integrator object
  type(BDF)                             :: bdfobj     ! BDF Integrator object
  type(ABM)                             :: abmobj     ! ABM Integrator object
  type(NBG)                             :: nbgobj     ! NBM Integrator object
    
  ! Declare physics for testing
  type(smd1)                   , target :: smd1obj    ! Spring-mass-damper test ODE (1 var)

  ! Declare functions
  type(kinetic_energy)         , target :: KE

  ! Design variable array
  type(scalar), dimension(:), allocatable :: x, dfdx, dfdxtmp
  type(scalar)                            :: fval, ftmp

  ! Use different step sizes for real and complex mode arithmetic
#if defined USE_COMPLEX
  real(dp)                                :: dh = 1.0d-16
#else
  real(dp)                                :: dh = 1.0d-11
#endif
  
  allocate(X(1), dfdx(1), dfdxtmp(1))
  !allocate(X(1))

  !===================================================================!
  !                          TEST BDF                                 !
  !===================================================================!

  X       = 0.0d0
  dfdx    = 0.0d0
  dfdxtmp = 0.0d0

  x(1) = 5.0d0 
  !x(2) = 0.20d0  ! damping coeff
  !x(3) = 5.0d0   ! stiffness coeff
  
  ! Initialize the system
  call smd1obj % initialize("SMD", num_state_vars = 1, num_design_vars = 1)

  bdfobj = BDF(system = smd1obj, tfinal = 10.0d0, h=1.0d-3, max_bdf_order = 1)
  call bdfobj % evalFuncGrad(num_func=1, func = KE,  num_dv = 1, x = x, fvals = fval, dfdx= dfdx)
  call bdfobj % writeSolution("bdf_forward_1.dat")
  call bdfobj % evalFDFuncGrad(num_func=1, func = KE,  num_dv = 1, x = x, fvals = fval, dfdx= dfdxtmp, dh=dh)
  call bdfobj % finalize()

  call smd1obj % finalize()

  print *, "fval         =", real_part(fval)
  print *, "Adjoint dfdx =", real_part(dfdx)
  print *, "FD      dfdx =", real_part(dfdxtmp)
  print *, "Error        =", abs(dfdxtmp-dfdx)
  print *, "Rel. Error   =", abs(real_part(dfdxtmp)-real_part(dfdx))/real_part(dfdxtmp)

  deallocate(X, dfdx, dfdxtmp)
  
end program main
