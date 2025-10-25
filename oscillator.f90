
  !===================================================================!
  !  Aeroelastic Oscillator
  !===================================================================!
  
  allocate(X(2), dfdx(2), dfdxtmp(2))

  dfdx    = 0.0d0
  dfdxtmp = 0.0d0

  x(1) = 9.0d0  ! dynamic pressure
  x(2) = 20.0d0 ! nonlinear stiffness coeff
  
  ! Initialize the system
  call oscillator % initialize("AE-Oscillator", num_state_vars = 2, num_design_vars = 2)
  
  bdfobj = BDF(system = oscillator, tfinal = 1.0d0, h=1.0d-3, max_bdf_order = 3) 

  call bdfobj % evalFuncGrad(num_func=1, func = pitch1,  num_dv = 2, x = x, &
       & fvals = fval, dfdx= dfdx)
  
!  call bdfobj % integrate()
  call bdfobj % writeSolution()
  
  call bdfobj % evalFDFuncGrad(num_func=1, func = pitch1,  num_dv = 2, x = x, &
       & fvals = fval, dfdx= dfdxtmp, dh=1.0d-6)

  call bdfobj % finalize()

  print*, "fval         =", fval
  print*, "Adjoint dfdx =", dfdx
  print*, "FD      dfdx =", dfdxtmp
  print*, "Error        =", abs(dfdxtmp-dfdx)

  ! Finalize the system
  call oscillator % finalize()

  dfdx = 0.0d0
  dfdxtmp = 0.0d0
  
  ! Initialize the system
  call oscillator % initialize("AE-Oscillator", num_state_vars = 2, num_design_vars = 2)
  
  dirkobj = DIRK(system = oscillator, tfinal = 1.d0, h=1.0d-3, num_stages=3) 
  
  call dirkobj % evalFuncGrad(num_func = 1, func = pitch1, num_dv = 2, x = x, &
       & fvals = fval, dfdx= dfdx)
  
  call dirkobj % writeSolution("dirksol.dat")

  call dirkobj % evalFDFuncGrad(num_func=1, func = pitch1, num_dv = 2, x = x, &
       & fvals = fval, dfdx= dfdxtmp, dh=dh)

  call dirkobj % finalize()  

  print *, "fval         =", fval
  print *, "Adjoint dfdx =", dfdx
  print *, "FD      dfdx =", dfdxtmp
  print *, "Error        =", abs(dfdxtmp-dfdx)

  ! Finalize the system
  call oscillator % finalize()
