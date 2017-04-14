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
  use vanderpol_class      , only : vanderpol

  ! Import functions
  use smd_functions_class  , only : kinetic_energy

  implicit none

  ! Declare integrators
  type(DIRK)                            :: dirkobj    ! DIRK Integrator object
  type(BDF)                             :: bdfobj     ! BDF Integrator object
  type(ABM)                             :: abmobj     ! ABM Integrator object
  type(NBG)                             :: nbgobj     ! NBM Integrator object

  ! Declare physics for testing
  type(vanderpol) , target :: vpl    ! Spring-mass-damper test ODE (1 var)

  ! Declare functions
  type(kinetic_energy) , target :: KE

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

  !===================================================================!
  !                       TEST ABM                                    !
  !===================================================================!

  fval    = 0.0d0
  X       = 0.0d0
  dfdx    = 0.0d0
  dfdxtmp = 0.0d0
  x(1)    = 1.0d0 

  ! Initialize the system
  call vpl % initialize("vanderpol", num_state_vars = 1, num_design_vars = size(x))
  abmobj = ABM(system = vpl, tfinal = 2.0d0, h=1.0d-3, max_abm_order = 6, second_order=.true.)
  call abmobj % evalFuncGrad(num_func=1, func = KE,  num_dv = size(x), x = x, fvals = fval, dfdx= dfdx)
  call abmobj % writeSolution("abm.dat")
  call abmobj % evalFDFuncGrad(num_func=1, func = KE,  num_dv = size(x), x = x, fvals = fval, dfdx= dfdxtmp, dh=dh)
  call abmobj % finalize()
  call vpl % finalize()

  print *, "fval         =", real_part(fval)
  print *, "Adjoint dfdx =", real_part(dfdx)
  print *, "FD      dfdx =", real_part(dfdxtmp)
  print *, "Error        =", abs(dfdxtmp-dfdx)
  print *, "Rel. Error   =", abs(real_part(dfdxtmp)-real_part(dfdx))/real_part(dfdxtmp)

  !===================================================================!
  !                       TEST NBG                                    !
  !===================================================================!

  fval    = 0.0d0
  X       = 0.0d0
  dfdx    = 0.0d0
  dfdxtmp = 0.0d0
  x(1)    = 1.0d0 

  ! Initialize the system
  call vpl % initialize("vanderpol", num_state_vars = 1, num_design_vars = size(x))
  nbgobj = NBG(system = vpl, tfinal = 2.0d0, h=1.0d-3, second_order=.true.)
  call nbgobj % evalFuncGrad(num_func=1, func = KE,  num_dv = size(x), x = x, fvals = fval, dfdx= dfdx)
  call nbgobj % writeSolution("nbg.dat")
  call nbgobj % evalFDFuncGrad(num_func=1, func = KE,  num_dv = size(x), x = x, fvals = fval, dfdx= dfdxtmp, dh=dh)
  call nbgobj % finalize()
  call vpl % finalize()

  print *, "fval         =", real_part(fval)
  print *, "Adjoint dfdx =", real_part(dfdx)
  print *, "FD      dfdx =", real_part(dfdxtmp)
  print *, "Error        =", abs(dfdxtmp-dfdx)
  print *, "Rel. Error   =", abs(real_part(dfdxtmp)-real_part(dfdx))/real_part(dfdxtmp)

  !===================================================================!
  !                          TEST BDF                                 !
  !===================================================================!

  fval    = 0.0d0
  X       = 0.0d0
  dfdx    = 0.0d0
  dfdxtmp = 0.0d0
  x(1)    = 1.0d0 

  ! Initialize the system
  call vpl % initialize("vanderpol", num_state_vars = 1, num_design_vars = size(x))
  bdfobj = BDF(system = vpl, tfinal = 2.0d0, h=1.0d-3, max_bdf_order = 1,second_order=.true.)
  call bdfobj % evalFuncGrad(num_func=1, func = KE,  num_dv = size(x), x = x, fvals = fval, dfdx= dfdx)
  call dirkobj % writeSolution("bdf.dat")
  call bdfobj % evalFDFuncGrad(num_func=1, func = KE,  num_dv = size(x), x = x, fvals = fval, dfdx= dfdxtmp, dh=dh)
  call bdfobj % finalize()
  call vpl % finalize()

  print *, "fval         =", real_part(fval)
  print *, "Adjoint dfdx =", real_part(dfdx)
  print *, "FD      dfdx =", real_part(dfdxtmp)
  print *, "Error        =", abs(dfdxtmp-dfdx)
  print *, "Rel. Error   =", abs(real_part(dfdxtmp)-real_part(dfdx))/real_part(dfdxtmp)

  !===================================================================!
  !                          TEST DIRK                                !
  !===================================================================!

  fval    = 0.0d0
  X       = 0.0d0
  dfdx    = 0.0d0
  dfdxtmp = 0.0d0
  x(1)    = 1.0d0 

  ! Initialize the system
  call vpl % initialize("vanderpol", num_state_vars = 1, num_design_vars = size(x))
  dirkobj = DIRK(system = vpl, tfinal = 2.0d0, h=1.0d-3, num_stages=3, second_order=.true.) 
  call dirkobj % testAdjoint6 ( num_func = 1, func = KE, num_dv = size(x), x = x, dfdx= dfdx, dfdxtmp=dfdxtmp )
  call dirkobj % writeSolution("dirk.dat")
  call dirkobj % evalFDFuncGrad(num_func=1, func = KE, num_dv = size(x), x = x, fvals = fval, dfdx= dfdxtmp, dh=dh)
  call dirkobj % finalize()  
  call vpl % finalize()

  print *, "fval         =", real_part(fval)
  print *, "Adjoint dfdx =", real_part(dfdx)
  print *, "FD      dfdx =", real_part(dfdxtmp)
  print *, "Error        =", abs(dfdxtmp-dfdx)
  print *, "Rel. Error   =", abs(real_part(dfdxtmp)-real_part(dfdx))/real_part(dfdxtmp)

  deallocate(X, dfdx, dfdxtmp)

end program main