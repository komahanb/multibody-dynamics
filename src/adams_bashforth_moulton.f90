#include "scalar.fpp"
!=====================================================================!
! Adams Bashworth Moulton Integration Module for first and second
! order systems with adjoint derivative capabilities.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module abm_integrator

  use integrator_class, only : integrator
  use physics_class,    only : physics

  implicit none

  private
  public :: ABM
  
  type(scalar), parameter :: ONE  = 1.0d0
  type(scalar), parameter :: ZERO = 0.0d0

  !===================================================================! 
  ! ABM Integrator type
  !===================================================================! 

  type, extends(integrator) :: ABM
     
     ! ABM variables
     integer :: max_abm_order = 3
     type(scalar), allocatable, dimension(:,:) :: A

   contains
     
     ! Destructor
     procedure, public  :: finalize

     ! Routines for integration
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff
     procedure, private :: getOrder

     ! Routines for adjoint gradient
     procedure, public  :: marchBackwards
     procedure, private :: assembleRHS
     procedure          :: computeTotalDerivative

  end type ABM

  interface ABM
     module procedure initialize
  end interface ABM

contains


  !===================================================================!
  ! Initialize the ABM datatype and allocate required variables
  !===================================================================!
  
  type(abm) function initialize( system, tinit, tfinal, h, second_order, max_abm_order )  result (this)
   
    class(physics), target          :: system
    integer  , OPTIONAL, intent(in) :: max_abm_order
    real(dp) , OPTIONAL, intent(in) :: tinit, tfinal
    real(dp) , OPTIONAL, intent(in) :: h
    logical  , OPTIONAL, intent(in) :: second_order

    print *, "======================================"
    print *, ">>   Adams Bashforth Moulton       << "
    print *, "======================================"

    call this % construct(system, tinit, tfinal, h, second_order)

    !-----------------------------------------------------------------!
    ! Set the order of integration
    !-----------------------------------------------------------------!

    if (present(max_abm_order)) then
       this % max_abm_order = max_abm_order
    end if
    print '("  >> Max ABM Order          : ",i4)', this % max_abm_order

    allocate( this % A (this % max_abm_order, this % max_abm_order) )
    this % A = 0.0d0 

    ! Set the coefficients
    if ( this % max_abm_order .eq. 1 ) then       
       this % A(1,1:1) = (/ 1.0d0 /)
    else if ( this % max_abm_order .eq. 2 ) then
       this % A(1,1:1) = (/ 1.0d0 /)
       this % A(2,1:2) = (/ 1.0d0/2.0d0, 1.0d0/2.0d0 /)
    else if ( this % max_abm_order .eq. 3 ) then
       this % A(1,1:1) = (/ 1.0d0 /)
       this % A(2,1:2) = (/ 1.0d0/2.0d0, 1.0d0/2.0d0 /)
       this % A(3,1:3) = (/ 5.0d0/12.0d0, 8.0d0/12.0d0, -1.0d0/12.0d0 /)
    else 
       print *,  "Wrong max_abm_order:", this % max_abm_order
       stop
    end if
    
    !-----------------------------------------------------------------!
    ! Setup adjoint RHS
    !-----------------------------------------------------------------!
    
    this % num_rhs_bins = this % max_abm_order

    allocate(this % rhs(this % num_rhs_bins, this % nsvars))
    this % rhs = 0.0d0

  end function initialize

  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!
  
  subroutine finalize( this )

    class(ABM) :: this

    ! Deallocate ABM coefficient
    if(allocated(this % A)) deallocate(this % A)
    
    if ( allocated(this % rhs) ) deallocate(this % rhs)

  end subroutine finalize

  !===================================================================!
  ! Returns the linearization scalar coefficients: alpha, beta, gamma
  !===================================================================!
  
  subroutine getLinearCoeff( this, alpha, beta, gamma )
   
    class(ABM)   :: this
    type(scalar), intent(out) :: alpha, beta, gamma
    integer :: k, m

    k = this % current_step
    m = this % getOrder(k)

    if ( this % second_order ) then
       gamma = 1.0d0
       beta  = this % A(m,1) * this % h
       alpha = beta * this % A(m,1) * this % h
    else 
       gamma = 0.0d0
       beta  = 1.0d0
       alpha = this  % A(m,1) * this % h
    end if
    
  end subroutine getLinearCoeff
  
  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!
  
  pure integer function getOrder(this, k)

    class(ABM), intent(in) :: this
    integer, intent(in)    :: k

    ! Find the order of approximation
    getOrder = k - 1

    ! Do not let exceed the max order sought
    if ( getOrder .gt. this % max_abm_order ) getOrder = this % max_abm_order

  end function getOrder
 
  !===================================================================!
  ! Compute the total derivative of the function with respect to the
  ! design variables and return the gradient 
  !===================================================================!
  
  subroutine computeTotalDerivative( this, dfdx )
    
    class(ABM)                                 :: this
    type(scalar) , dimension(:), intent(inout) :: dfdx
    type(scalar) , allocatable, dimension(:,:) :: dRdX
    type(scalar)                               :: scale = 1.0d0
    integer                                    :: k

    if (.not.allocated(dRdX)) allocate(dRdX(this % nSVars, this % nDVars))
    
    dRdX = 0.0d0
    dfdx = 0.0d0
    
    !-----------------------------------------------------------------!
    ! Compute dfdx
    !-----------------------------------------------------------------!

    do k = 2, this % num_steps
       call this % system % func % addFuncDVSens(dfdx, scale, this % time(k), &
            & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:) )
    end do
    
    ! Initial condition
!!$    call this % system % func % addFuncDVSens(dfdx, scale, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(1,:) )
    
    !-----------------------------------------------------------------!
    ! Compute the total derivative
    !-----------------------------------------------------------------!
    
    do k = 2, this % num_steps
       call this % system % getResidualDVSens(dRdX, scale, this % time(k), &
            & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:))
       dfdx = dfdx + matmul( transpose(dRdX), this % lambda(k,:)) ! check order
    end do

!!$    ! Add constraint contribution
!!$    call this % system % getResidualDVSens(dRdX, scale, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(1,:))
!!$    dfdx = dfdx + matmul(this % lambda(2,:), dRdX)

    ! Finally multiply by the scalar
    dfdx = this % h * dfdx

    if (allocated(dRdX)) deallocate(dRdX)

  end subroutine computeTotalDerivative

  !===================================================================!
  ! Subroutine that marches backwards in time to compute the lagrange
  ! multipliers (adjoint variables for the function)
  ! ===================================================================!
  
  subroutine marchBackwards( this )

    class(ABM)                :: this
    integer                   :: k
    type(scalar)              :: alpha, beta, gamma
    
    time: do k = this % num_steps, 2, -1
       
       this % current_step = k 
       
       !--------------------------------------------------------------!
       ! Determine the linearization coefficients for the Jacobian
       !--------------------------------------------------------------!
              
       call this % getLinearCoeff(alpha, beta, gamma)

       !--------------------------------------------------------------!
       ! Solve the adjoint equation at each step
       !--------------------------------------------------------------!

       call this % adjointSolve(this % lambda(k,:), alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))
       
       print*, "k,lambda=", k, this % lambda(k,:)

    end do time
    
  end subroutine marchBackwards
  
  !===================================================================!
  ! Approximate the state variables at each step using ABM formulae
  !===================================================================!
  
  subroutine approximateStates( this )

    class(ABM)   :: this
    integer      :: k, m, i
    type(scalar) :: scale

    k = this % current_step
    
    m = this % getOrder(k)

    ! Approximate UDDOT
    this % uddot(k,:) = this % uddot(k-1,:)
    
    ! Approximate UDOT
    this % udot(k,:) = this % udot(k-1,:)
    
    do i = 0, m-1
       scale = this % h * this % A(m,i+1)
       this % udot(k,:) = this % udot(k,:) + scale * this % uddot(k-i,:)
    end do

    ! Approximate U
    this % u(k,:) = this % u(k-1,:)
    
    do i = 0, m-1
       scale = this % h * this % A(m,i+1)
       this % u(k,:) = this % u(k,:) + scale * this % udot(k-i,:)
    end do

    
  end subroutine approximateStates

  !===================================================================!
  ! Function that puts together the right hand side of the adjoint
  ! equation into the supplied rhs vector.
  !===================================================================!
  
  subroutine assembleRHS( this, rhs )

    class(ABM)                                :: this
    type(scalar), dimension(:), intent(inout) :: rhs
    type(scalar), dimension(:,:), allocatable :: jac
    type(scalar)                              :: scale
    integer                                   :: k, m
    type(scalar)                              :: alpha, beta, gamma
    
    allocate( jac(this % nSVars, this % nSVars) )

    ! Zero the RHS first
    rhs = 0.0d0

    ! Retrieve the current time index
    k = this % current_step
    m = this % getOrder(k)

    !-------------------------------------------------------------------!
    ! Find PHI at the new step
    !-------------------------------------------------------------------!

    if ( k+1 .le. this % num_steps ) then

       gamma = 0.0d0
       beta  = 0.0d0
       alpha = this % h

       ! Add the state variable sensitivity from the previous step
       call this % system % func % addFuncSVSens(this % phi(k,:), &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % system % X, &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       ! Add the residual adjoint product from the previous step
       call this % system % assembleJacobian(jac, &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       this % phi(k,:) = this % phi(k,:) + matmul(transpose(jac(:,:)), this % lambda(k+1,:))
       
    else
       
       this % phi(k,:) = 0.0d0
       
    end if
    
    !-------------------------------------------------------------------!
    ! Find PSI at the new step
    !-------------------------------------------------------------------!
    
    if ( k+1 .le. this % num_steps ) then

       !    do i = 1, m + 1 ! 0, m

       !      if ( k+i .le. this % num_steps ) then

       ! Add contributions from previous PSI
       !        if (i .eq. 1) then
       this % psi(k,:) = this % psi(k,:) + this % psi(k+1,:)
       !       end if

       ! Add contributions from PHI
       this % psi(k,:) = this % psi(k,:) + this % h * this % A(m,1) * this % phi(k,:)
       this % psi(k,:) = this % psi(k,:) + this % h * this % A(m,1) * this % phi(k+1,:)

       !   end if

       ! end do

!!$
!!$    do i = 2, m
!!$
!!$       if ( k+i .le. this % num_steps ) then

       gamma = 0.0d0
       beta  = this % h
       alpha = this % h * this % h * this % A(m,1)

       ! Add the state variable sensitivity from the previous step
       call this % system % func % addFuncSVSens(this % psi(k,:), &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % system % X, &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       ! Add the residual adjoint product from the previous step
       call this % system % assembleJacobian(jac, &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       this % psi(k,:) = this % psi(k,:) + matmul(transpose(jac(:,:)), this % lambda(k+1,:))
!!$
!!$       end if
!!$
!!$    end do

    else

       this % psi(k,:) = 0.0d0

    end if

    !--------------------------------------------------------------------------!
    ! Add up contribution
    !--------------------------------------------------------------------------!
    
    gamma = 1.0d0
    beta  = this % h * this % A(m,1)
    alpha = this % h * this % A(m,1) * this % h * this % A(m,1)

    call this % system % func % addFuncSVSens(rhs, &
         & alpha, beta, gamma, &
         & this % time(k), &
         & this % system % X, &
         & this % u(k,:), &
         & this % udot(k,:), &
         & this % uddot(k,:))

    if ( k+1 .le. this % num_steps ) then

       gamma = 0.0d0
       beta  = this % h * this % A(m,1)
       alpha = this % h * this % A(m,1) * this % h * this % A(m,1)

       call this % system % func % addFuncSVSens(rhs, &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % system % X, &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       ! Add the residual adjoint product from the previous step
       call this % system % assembleJacobian(jac, &
            & alpha, beta, gamma, &
            & this % time(k+1), &
            & this % u(k+1,:), &
            & this % udot(k+1,:), &
            & this % uddot(k+1,:))

       rhs = rhs + matmul(transpose(jac(:,:)), this % lambda(k+1,:))

       ! Add previous contributions
       rhs = rhs + this % A(m,1) * this % psi(k+1,:)

       rhs = rhs + this % A(m,1) * this % h * this % A(m,1) * this % phi(k+1,:)
       rhs = rhs + this % A(m,1) * this % h * this % A(m,1) * this % phi(k+1,:)
       
    end if

    ! current
    rhs = rhs + this % A(m,1) * this % psi(k,:)
    rhs = rhs + this % A(m,1) * this % h * this % A(m,1) * this % phi(k,:)

     ! Negate the RHS
    rhs = - rhs

    if(allocated(jac)) deallocate(jac)

  end subroutine assembleRHS

end module abm_integrator
