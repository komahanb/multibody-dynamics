#include "scalar.fpp"

!=====================================================================!
! Module that implements solving a nonlinear problem using different
! methods such as Newton's method.
!=====================================================================!

module nonlinear_algebra

  ! import dependencies
  use iso_fortran_env , only : wp => REAL64
  use linear_algebra  , only : solve
  use physics_class   , only : physics
  use utils           , only : norm

  ! disable implicit datatypes
  implicit none

  ! Define constants used  
  real(dp) :: abs_res_tol          = 1.0d-14
  real(dp) :: rel_res_tol          = 1.0d-10

  integer  :: max_newton_iters     = 5
  logical  :: approximate_jacobian = .false.
  logical  :: jacobian_check       = .true.
  logical  :: second_order         = .true.
  integer  :: print_level          = 0

  public :: nonlinear_solve

  !-------------------------------------------------------------------!
  ! Interface for nonlinear solution problems
  !-------------------------------------------------------------------!

  interface nonlinear_solve
     !module procedure newton_solve_first_order
     module procedure solve_second_order
     !module procedure solve_first_order
     !module procedure newton_solve_condensed
  end interface nonlinear_solve

contains

  !===================================================================!
  ! Solve the nonlinear system at each step by driving the
  ! residual to zero
  !
  ! Input: 
  ! The guessed (initial) state variable values q, qdot, qddot are
  ! supplied
  !
  ! Output:
  ! q, qdot, qddot updated iteratively until the corresponding
  ! residual R = 0
  !
  ! alpha: multiplier for derivative of Residual wrt to q
  ! beta : multiplier for derivative of Residual wrt to qdot
  ! gamma: multiplier for derivative of Residual wrt to qddot
  !===================================================================!

  subroutine solve_second_order( system, alpha, beta, gamma, t, q, qdot, qddot )

    class(physics) :: system

    ! Arguments
    type(scalar), intent(in)                  :: alpha, beta, gamma
    real(dp), intent(in)                      :: t
    type(scalar), intent(inout), dimension(:) :: q, qdot, qddot

    ! Norms for tracking progress
    real(dp)                                  :: abs_res_norm = 0.0d0
    real(dp)                                  :: rel_res_norm = 0.0d0
    real(dp)                                  :: init_norm    = 0.0d0

    ! Other Local variables
    type(scalar), allocatable, dimension(:)   :: res, dq
    type(scalar), allocatable, dimension(:,:) :: jac, fd_jac

    integer                                   :: n, nvars
    logical                                   :: conv = .false.

    type(scalar)                              :: jac_err

    ! find the size of the linear system based on the calling object
    nvars = size(q)

    if ( .not. allocated(res)    ) allocate( res(nvars)         )
    if ( .not. allocated(dq)     ) allocate( dq(nvars)          )
    if ( .not. allocated(jac)    ) allocate( jac(nvars,nvars)    )
    if ( .not. allocated(fd_jac) ) allocate( fd_jac(nvars,nvars) )

    if ( print_level .ge. 1) then
       write(*,'(A10, 2A12)') "NewtonIter", "|R|", "|R|/|R1|"
    end if

    newton: do n = 1, max_newton_iters
       
       call system % assembleResidual(res, t, q, qdot, qddot)

       ! Get the jacobian matrix
       if ( approximate_jacobian ) then

          ! Compute an approximate Jacobian using finite differences
          call approximateJacobian(system, jac, alpha, beta, gamma, t, q, qdot, qddot)

       else

          ! Use the user supplied Jacobian implementation
          call system % assembleJacobian(jac, alpha, beta, gamma, t, q, qdot, qddot)

          ! Check the Jacobian implementation once at the beginning of integration
          if ( jacobian_check .and. n .eq. 1 ) then

             ! Compute an approximate Jacobian using finite differences
             call approximateJacobian(system, fd_jac, alpha, beta, gamma, t, q, qdot, qddot)

             ! Compare the exact and approximate Jacobians and
             ! complain about the error in Jacobian if there is any
             jac_err = maxval(abs(fd_jac - jac))
             if ( abs(jac_err) .gt. 1.0d-3 ) then
                print *, "q     =", q
                print *, "qdot  =", qdot
                print *, "qddot =", qddot
                print *, "a,b,c =", alpha, beta, gamma
                print *, "J     =", jac
                print *, "Jhat  =", fd_jac
                print *, "WARNING: Possible error in jacobian", jac_err
             end if

             ! Set that the jacobian is checked
             jacobian_check = .false.

          end if

       end if

       ! Find norm of the residual
       abs_res_norm = norm(res)
       if ( n .eq. 1) init_norm = abs_res_norm
       rel_res_norm = abs_res_norm/init_norm

       if ( print_level .eq. 2) then
          write(*, "(I10,2ES12.2)") n, abs_res_norm, rel_res_norm
       end if

       ! Check stopping
       if ((abs_res_norm .le. abs_res_tol) .or. (rel_res_norm .le. rel_res_tol)) then
          conv = .true.
          exit newton
       else if ((abs_res_norm .ne. abs_res_norm) .or. (rel_res_norm .ne. rel_res_norm) ) then
          conv = .false.
          exit newton
       end if

       ! Call LAPACK to solve the linear system
       dq = solve(jac, -res)

       ! Update the solution
       qddot = qddot + gamma * dq
       qdot  = qdot  + beta  * dq
       q     = q     + alpha * dq

    end do newton

    if (print_level .eq. 1) then 
       write(*, "(I10,2ES12.2)") n, abs_res_norm, rel_res_norm
    end if

    ! Print warning message if not converged
    if (.not. conv) then
       write(*,'(A5, 2A12)') "ITER", "|R|", "|R|/|R1|"
       write(*, "(I5,2ES12.2)") n, abs_res_norm, rel_res_norm
       stop "Newton Solve Failed"
    end if

    if (allocated(res)) deallocate(res)
    if (allocated(dq)) deallocate(dq)
    if (allocated(jac)) deallocate(jac)
    if (allocated(fd_jac)) deallocate(fd_jac)

  end subroutine solve_second_order

!!$  subroutine newton_solve_condensed(system, coeff, U)
!!$
!!$    use physics_class, only  : physics
!!$
!!$    class(physics)                            :: system
!!$
!!$    ! Arguments
!!$    type(scalar), intent(in)                  :: coeff(:)
!!$    type(scalar), intent(inout)               :: U(:,:)
!!$
!!$    ! Other Local variables
!!$    type(scalar), allocatable, dimension(:)   :: res, dq
!!$    type(scalar), allocatable, dimension(:,:) :: jac, fd_jac
!!$
!!$    integer                                   :: n, k, nvars
!!$    logical                                   :: conv = .false.
!!$

  !===================================================================! 
  ! Routine that approximates the Jacobian based on finite differences
  ! [d{R}/d{q}] = alpha*[dR/dq] + beta*[dR/dqdot] + gamma*[dR/dqddot]
  !===================================================================!

  subroutine approximateJacobian( system, jac, alpha, beta, gamma, t, q, qdot, qddot )

    class(physics)                            :: system

    ! Matrices
    type(scalar) , intent(inout), dimension(:,:) :: jac

    ! Arrays
    real(dp)     , intent(in)                    :: t
    type(scalar) , intent(in), dimension(:)      :: q, qdot, qddot     ! states

    type(scalar) , allocatable, dimension(:)     :: pstate             ! perturbed states
    type(scalar) , allocatable, dimension(:)     :: R, Rtmp            ! original residual and perturbed residual

    ! Scalars
    type(scalar)                                 :: dh = 1.0d-6        ! finite-diff step size
    type(scalar) , intent(in)                    :: alpha, beta, gamma ! linearization coefficients
    integer                                      :: m                  ! loop variables
    integer :: nvars

    !  Zero the supplied jacobian matrix for safety (as we are
    !  computing everything newly here)
    jac = 0.0d0

    nvars = size(q)

    ! Allocate required arrays
    allocate(pstate(nvars)); pstate = 0.0d0;
    allocate(R(nvars));      R = 0.0d0;
    allocate(Rtmp(nvars));   Rtmp = 0.0d0;

    ! Make a residual call with original variables
    call system % assembleResidual(R, t, q, qdot, qddot)

    !-----------------------------------------------------------!
    ! Derivative of R WRT Q: dR/dQ
    !-----------------------------------------------------------!

    pstate = q

    loop_vars: do m = 1, nvars

       ! Perturb the k-th variable
       pstate(m) = pstate(m) + dh

       ! Make a residual call with the perturbed variable
       call system % assembleResidual(Rtmp, t, pstate, qdot, qddot)

       ! Unperturb (restore) the k-th variable
       pstate(m) =  q(m)

       ! Approximate the jacobian with respect to the k-th variable
       jac(:,m) = jac(:,m) + alpha*(Rtmp-R)/dh

    end do loop_vars

    !-----------------------------------------------------------!
    ! Derivative of R WRT QDOT: dR/dQDOT
    !-----------------------------------------------------------!

    pstate = qdot

    do m = 1, nvars

       ! Perturb the k-th variable
       pstate(m) = pstate(m) + dh

       ! Make a residual call with the perturbed variable
       call system % assembleResidual(Rtmp, t, q, pstate, qddot)

       ! Unperturb (restore) the k-th variable
       pstate(m) =  qdot(m)

       ! Approximate the jacobian with respect to the k-th variable
       Jac(:,m) = Jac(:,m) + beta*(Rtmp-R)/dh

    end do

    ! Second order equations have an extra block to add
    if (second_order) then

       !-----------------------------------------------------------!
       ! Derivative of R WRT QDDOT: dR/dQDDOT
       !-----------------------------------------------------------!     

       pstate = qddot

       do m = 1, nvars

          ! Perturb the k-th variable
          pstate(m) = pstate(m) + dh

          ! Make a residual call with the perturbed variable
          call system % assembleResidual(Rtmp, t, q, qdot, pstate)

          ! Unperturb (restore) the k-th variable
          pstate(m) =  qddot(m)

          ! Approximate the jacobian with respect to the k-th variable
          Jac(:,m) = Jac(:,m) + gamma*(Rtmp-R)/dh

       end do

    end if ! first or second order

    deallocate(pstate)
    deallocate(R,Rtmp)

  end subroutine approximateJacobian

end module nonlinear_algebra
