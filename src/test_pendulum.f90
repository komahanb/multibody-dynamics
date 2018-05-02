!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

#include "scalar.fpp"

program main

  use utils , only : real_part

  ! Integrators
  use bdf_integrator , only : BDF
  use abm_integrator , only : ABM

  ! Residual and functions
  use pendulum_class, only : pendulum
  use pendulum_functions_class, only : time_period

  implicit none

  ! Declare integrators
  type(BDF) :: bdfobj
  type(ABM) :: abmobj

  ! Physics and function
  type(pendulum), target :: pend
  type(time_period), target :: time

  ! Design variable array
  type(scalar), allocatable :: x(:), dfdx(:), dfdxtmp(:)
  type(scalar) :: fval, ftmp
  type(scalar) :: tend
  type(integer) :: k

  ! Use different step sizes for real and complex mode arithmetic
#if defined USE_COMPLEX
  real(dp) :: dh = 1.0d-16
#else
  real(dp) :: dh = 1.0d-7
#endif

  integer, parameter :: max_order = 3
  integer, parameter :: num_vars = 1

  allocate(X(num_vars), dfdx(num_vars), dfdxtmp(num_vars))

  adjoint: block

    fval    = 0.0d0
    X       = 0.0d0
    dfdx    = 0.0d0
    dfdxtmp = 0.0d0
    x(1)    = 1.0d0

    ! Initialize the system
    call pend % initialize("pendulum", num_state_vars = 1, num_design_vars = size(x))
    bdfobj = BDF(system = pend, tfinal = 3.14d0, h=1.0d-3, max_bdf_order = 2, second_order = .true.)

    call bdfobj % evalFuncGrad(num_func=1, func = time,  num_dv = size(x), x = x, fvals = fval, dfdx = dfdx)
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

  deallocate(X,dfdx,dfdxtmp)

end program main
