!=====================================================================!
! Main Program for the landing event simulation
!=====================================================================!

program main

  use iso_fortran_env, only: dp => real64

  use runge_kutta_integrator, only : DIRK
  use bdf_integrator, only : BDF

  use myode_class, only : ODE

  implicit none

  type(DIRK) :: dirkobj
  type(BDF)  :: bdfobj
  type(ODE), target :: myode

  call dirkobj % setPhysicalSystem(myode)
  call dirkobj % initialize(tfinal = 25.0d0, num_stages=1, &
       & h=0.01d0, nsvars=1, second_order=.true.) ! all are optional except nvars
  call dirkobj % setPrintLevel(0)
  call dirkobj % integrate()
  !  call dirkobj % integrateBackward()
  call dirkobj % writeSolution('dirk.dat')
  call dirkobj % finalize()

  call bdfobj % setPhysicalSystem(myode)
  call bdfobj % initialize(tfinal = 25.0d0, max_bdf_order=3, &
       & h=0.01d0, nsvars=1, second_order=.true.) ! all are optional except nvars
  call bdfobj % setPrintLevel(0)
  call bdfobj % integrate()
  ! call bdfobj % integrateBackward()
  call bdfobj % writeSolution('bdf.dat')
  call bdfobj % finalize()

end program





