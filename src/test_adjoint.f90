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
  
  allocate(X(3), dfdx(3), dfdxtmp(3))

  !===================================================================!
  !                       TEST ABM                                    !
  !===================================================================!
  
  X       = 0.0d0
  dfdx    = 0.0d0
  dfdxtmp = 0.0d0

  x(1) = 2.50d0   ! mass
  x(2) = 0.20d0    ! damping coeff
  x(3) = 5.0d0   ! stiffness coeff

  ! Initialize the system
  call smd1obj % initialize("SMD", num_state_vars = 1, num_design_vars = 3)
  abmobj = ABM(system = smd1obj, tfinal = 2.0d0, h=1.0d-3, max_abm_order = 3)
  call abmobj % evalFuncGrad(num_func=1, func = KE,  num_dv = 3, x = x, fvals = fval, dfdx= dfdx)
  call abmobj % writeSolution("abm.dat")
  call abmobj % evalFDFuncGrad(num_func=1, func = KE,  num_dv = 3, x = x, fvals = fval, dfdx= dfdxtmp, dh=dh)
  call abmobj % finalize()
  call smd1obj % finalize()

  print *, "fval         =", real_part(fval)
  print *, "Adjoint dfdx =", real_part(dfdx)
  print *, "FD      dfdx =", real_part(dfdxtmp)
  print *, "Error        =", abs(dfdxtmp-dfdx)
  print *, "Rel. Error   =", abs(real_part(dfdxtmp)-real_part(dfdx))/real_part(dfdxtmp)

  !===================================================================!
  !                       TEST NBG                                    !
  !===================================================================!

  fval = 0.0d0
  X       = 0.0d0
  dfdx    = 0.0d0
  dfdxtmp = 0.0d0

  x(1) = 2.50d0    ! mass
  x(2) = 0.20d0    ! damping coeff
  x(3) = 5.0d0    ! stiffness coeff

  ! Initialize the system
  call smd1obj % initialize("SMD", num_state_vars = 1, num_design_vars = 3)
  nbgobj = NBG(system = smd1obj, tfinal = 2.0d0, h=1.0d-3)
  call nbgobj % evalFuncGrad(num_func=1, func = KE,  num_dv = 3, x = x, fvals = fval, dfdx= dfdx)  
  call nbgobj % writeSolution("nbg.dat")
  call nbgobj % evalFDFuncGrad(num_func=1, func = KE,  num_dv = 3, x = x, fvals = fval, dfdx= dfdxtmp, dh=dh)
  call nbgobj % finalize()
  call smd1obj % finalize()

  print *, "fval         =", real_part(fval)
  print *, "Adjoint dfdx =", real_part(dfdx)
  print *, "FD      dfdx =", real_part(dfdxtmp)
  print *, "Error        =", abs(dfdxtmp-dfdx)
  print *, "Rel. Error   =", abs(real_part(dfdxtmp)-real_part(dfdx))/real_part(dfdxtmp)

  !===================================================================!
  !                          TEST BDF                                 !
  !===================================================================!

  X       = 0.0d0
  dfdx    = 0.0d0
  dfdxtmp = 0.0d0

  x(1) = 2.50d0  ! mass
  x(2) = 0.20d0  ! damping coeff
  x(3) = 5.0d0   ! stiffness coeff
  
  ! Initialize the system
  call smd1obj % initialize("SMD", num_state_vars = 1, num_design_vars = 3)
  bdfobj = BDF(system = smd1obj, tfinal = 2.0d0, h=1.0d-3, max_bdf_order = 3)
  call bdfobj % evalFuncGrad(num_func=1, func = KE,  num_dv = 3, x = x, fvals = fval, dfdx= dfdx)
  call bdfobj % evalFDFuncGrad(num_func=1, func = KE,  num_dv = 3, x = x, fvals = fval, dfdx= dfdxtmp, dh=dh)
  call bdfobj % finalize()
  call smd1obj % finalize()

  print *, "fval         =", real_part(fval)
  print *, "Adjoint dfdx =", real_part(dfdx)
  print *, "FD      dfdx =", real_part(dfdxtmp)
  print *, "Error        =", abs(dfdxtmp-dfdx)
  print *, "Rel. Error   =", abs(real_part(dfdxtmp)-real_part(dfdx))/real_part(dfdxtmp)

  !===================================================================!
  !                          TEST DIRK                                !
  !===================================================================!

  X       = 0.0d0
  dfdx    = 0.0d0
  dfdxtmp = 0.0d0
  
  x(1) = 2.50d0  ! mass
  x(2) = 0.20d0  ! damping coeff
  x(3) = 5.0d0   ! stiffness coeff

  ! Initialize the system
  call smd1obj % initialize("SMD", num_state_vars = 1, num_design_vars = 3)
  dirkobj = DIRK(system = smd1obj, tfinal = 2.0d0, h=1.0d-3, num_stages=3, second_order=.true.) 
  !  call dirkobj % testAdjointr( num_func = 1, func = KE, num_dv = 3, x = x,dfdx= dfdx)
  call dirkobj % testAdjoint6 ( num_func = 1, func = KE, num_dv = 3, &
       & x = x, dfdx= dfdx, dfdxtmp=dfdxtmp )

  !  call dirkobj % testAdjoint2( num_func = 1, func = KE, num_dv = 3, x = x,dfdx= dfdx)
  call dirkobj % writeSolution("dirksol.dat")

  !   call dirkobj % evalFuncGrad( num_func=1, func = KE, num_dv = 3, x = x, &
  !        & fvals = fval, dfdx= dfdx )

  call dirkobj % evalFDFuncGrad(num_func=1, func = KE, num_dv = 3, x = x, &
       & fvals = fval, dfdx= dfdxtmp, dh=dh)

  print *, "fval         =", real_part(fval)
  print *, "Adjoint dfdx =", real_part(dfdx)
  print *, "FD      dfdx =", real_part(dfdxtmp)
  print *, "Error        =", abs(dfdxtmp-dfdx)
  print *, "Rel. Error   =", abs(real_part(dfdxtmp)-real_part(dfdx))/real_part(dfdxtmp)

!!$  call dirkobj % evalFuncGrad(num_func=1, func = KE, num_dv = 3, x = x, &
!!$       & fvals = fval, dfdx= dfdx)
!!$  
!!$  call dirkobj % writeSolution("dirksol.dat")
!!$
!!$  call dirkobj % evalFDFuncGrad(num_func=1, func = KE, num_dv = 3, x = x, &
!!$       & fvals = fval, dfdx= dfdxtmp, dh=dh)

  call dirkobj % finalize()  
  call smd1obj % finalize()

  deallocate(X, dfdx, dfdxtmp)
  
end program main
