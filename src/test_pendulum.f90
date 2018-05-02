!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

#include "scalar.fpp"

program main

  use utils , only : real_part
  
  use runge_kutta_integrator , only : DIRK
  use bdf_integrator , only : BDF
  use abm_integrator , only : ABM
  use nbg_integrator , only : NBG
  
  !use vanderpol_class , only : vanderpol
  !use smd_functions_class , only : kinetic_energy

  use pendulum_class, only : pendulum
  use pendulum_functions_class, only : time_period
  
  implicit none

  ! Declare integrators
  type(DIRK) :: dirkobj
  type(BDF)  :: bdfobj
  type(ABM)  :: abmobj
  type(NBG)  :: nbgobj

  ! Physics and function
  type(pendulum) , target :: pend
  type(time_period) , target :: time
 
  ! Design variable array
  type(scalar), allocatable :: x(:), dfdx(:), dfdxtmp(:)
  type(scalar)              :: fval, ftmp
  type(scalar)              :: tend
  
  type(integer) :: k
  ! Use different step sizes for real and complex mode arithmetic
#if defined USE_COMPLEX
  real(dp) :: dh = 1.0d-25
#else
  real(dp) :: dh = 1.0d-11
#endif

  integer, parameter :: maxk = 6, max_order=3
  integer      :: num_steps(maxk+1), i, p
  type(real(8)) :: abserror(maxk+1,1)
  type(real(8)) :: relerror(maxk+1,1)

  allocate(X(1), dfdx(1), dfdxtmp(1))
  
  adjoint: block

    fval    = 0.0d0
    X       = 0.0d0
    dfdx    = 0.0d0
    dfdxtmp = 0.0d0
    x(1)    = 1.0d0 

    ! Initialize the system
    call pend % initialize("pendulum", num_state_vars = 1, num_design_vars = size(x))
    bdfobj = BDF(system = pend, tfinal = 2.0d0, h=1.0d-3, max_bdf_order = 1 ,second_order=.true.)    
    call bdfobj % evalFuncGrad(num_func=1, func = time,  num_dv = size(x), x = x, fvals = fval, dfdx= dfdx)
    call bdfobj % writeSolution("bdf.dat")
    call bdfobj % evalFDFuncGrad(num_func=1, func = time,  num_dv = size(x), x = x, fvals = fval, dfdx= dfdxtmp, dh=dh)
    call bdfobj % finalize()
    call pend % finalize()

    print *, "fval         = ", real_part(fval)
    print *, "Adjoint dfdx = ", real_part(dfdx)
    print *, "FD      dfdx = ", real_part(dfdxtmp)
    print *, "Error        = ", abs(dfdxtmp-dfdx)
    print *, "Rel. Error   = ", abs(real_part(dfdxtmp)-real_part(dfdx))/real_part(dfdxtmp)

  end block adjoint

end program main
