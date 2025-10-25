#include "scalar.fpp"

!=====================================================================!
! Backward Difference Formula Integration Module for first and second
! order systems with adjoint derivative capabilities.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module bdf_integrator

  use integrator_class, only : integrator
  use physics_class,    only : physics
  use utils, only : real_part

  implicit none

  private
  public :: BDF

  ! Useful constants
  type(scalar), parameter :: ONE  = 1.0d0
  type(scalar), parameter :: ZERO = 0.0d0
  
  !===================================================================! 
  ! BDF Integrator type
  !===================================================================! 

  type, extends(integrator) :: BDF
     
     ! BDF variables
     integer            :: max_bdf_order = 6
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

  end type BDF  

  ! Constructor for BDF type
  interface BDF
     module procedure initialize
  end interface BDF

contains
    
  !===================================================================!
  ! Returns the linearization scalar coefficients: alpha, beta, gamma
  !===================================================================!
  
  subroutine getLinearCoeff( this, alpha, beta, gamma )
    
    class(BDF) :: this
    type(scalar), intent(out)   :: alpha, beta, gamma
    integer :: k, m

    k = this % current_step
    m = this % getOrder(k)

    if ( this % second_order ) then
       alpha = 1.0d0
       beta  = this % A(m,1) / this % h
       gamma = beta * this % A(m,1) / this % h
    else 
       alpha = 1.0d0
       beta  = this  % A(m,1) / this % h
       gamma = 0.0d0
    end if
    
  end subroutine getLinearCoeff

  !===================================================================!
  ! Approximate the state variables at each step using BDF formulae
  !===================================================================!
  
  subroutine approximateStates( this )

    class(BDF) :: this
    integer    :: k, m, i
    type(scalar) :: scale

    k = this % current_step
    m = this % getOrder(k)

    !-----------------------------------------------------------------!
    ! Extrapolate U to next time step
    !-----------------------------------------------------------------!
    
    this % u(k,:) = this % u(k-1,:) + this % udot(k-1,:)*this % h &
         & + this % uddot(k-1,:)*this % h*this % h/2.0d0

    !-----------------------------------------------------------------!
    ! Approximate UDOT using BDF
    !-----------------------------------------------------------------!
    
    do i = 1, m + 1
       scale = this % A(m,i)/this %h
       this % udot(k,:) = this % udot(k,:) + scale*this % u(k-i+1,:)
    end do
    
    !-----------------------------------------------------------------!
    ! Approximate UDDOT using BDF
    !-----------------------------------------------------------------!

    do i = 1, m + 1
       scale = this % A(m,i)/this%h*this % A(m,i)/this%h
       this % uddot(k,:) = this % uddot(k,:) + scale*this % udot(k-i+1,:)
    end do
   
!!$    m = this % getOrder(k)
!!$
!!$    if ( m .eq. 0) then
!!$
!!$       ! We dont have enought points yet
!!$       this % uddot(k,:) = (this % udot(k,:) - this % udot(k-1,:))/this % h
!!$       
!!$    else
!!$       
!!$       do i = 1, 2*m + 1
!!$          this % uddot(k,:) = this % uddot(k,:) &
!!$               & + this % A(m, i)*this % u(k-i+1,:)/this % h/this % h
!!$       end do
!!$       
!!$    end if
    
  end subroutine approximateStates

  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!
  
  pure integer function getOrder(this, k)

    class(bdf), intent(in) :: this
    integer, intent(in)    :: k

    ! find the order of approximation
    getOrder = k-1 ! k = md + 1

    ! Do not let exceed the max order sought
    if ( getOrder .gt. this % max_bdf_order ) getOrder = this % max_bdf_order

  end function getOrder
 
  !===================================================================!
  ! Compute the total derivative of the function with respect to the
  ! design variables and return the gradient 
  !===================================================================!
  
  subroutine computeTotalDerivative( this, dfdx )
    
    class(BDF)                                         :: this
    type(scalar) , dimension(:), intent(inout)             :: dfdx
    type(scalar) , dimension(this % nSVars, this % nDVars) :: dRdX
    type(scalar)                                           :: scale = 1.0d0
    integer                                            :: k
    
!    scale = this % h
    
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
       dfdx = dfdx + matmul(this % psi(k,:), dRdX) ! check order
    end do

!!$    ! Add constraint contribution
!!$    call this % system % getResidualDVSens(dRdX, scale, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(1,:))
!!$    dfdx = dfdx + matmul(this % psi(2,:), dRdX)

    ! Finally multiply by the scalar
    dfdx = this % h * dfdx

    print*, "Check scaling of dfdx, and transpose"

  end subroutine computeTotalDerivative

  !===================================================================!
  ! Initialize the BDF datatype and allocate required variables
  !===================================================================!
  
  type(bdf) function initialize( system, tinit, tfinal, h, second_order, max_bdf_order )  result (this)
   
    class(physics), target          :: system
    integer  , OPTIONAL, intent(in) :: max_bdf_order
    real(dp) , OPTIONAL, intent(in) :: tinit, tfinal
    real(dp) , OPTIONAL, intent(in) :: h
    logical  , OPTIONAL, intent(in) :: second_order

    print *, "======================================"
    print *, ">>   Backward Difference Formula    << "
    print *, "======================================"

    call this % construct(system, tinit, tfinal, h, second_order)
    
    !-----------------------------------------------------------------!
    ! Set the order of integration
    !-----------------------------------------------------------------!

    if (max_bdf_order .le. this % max_bdf_order) this % max_bdf_order = max_bdf_order
    print '("  >> Max BDF Order          : ",i4)', this % max_bdf_order

    allocate( this % A (this % max_bdf_order, this % max_bdf_order+1) )
    this % A = 0.0d0 

    ! Set the coefficients
    ! http://www.scholarpedia.org/article/Backward_differentiation_formulas
    if ( this % max_bdf_order .ge. 1 ) this % A(1,1:2) = [1.0d0, -1.0d0]
    if ( this % max_bdf_order .ge. 2 ) this % A(2,1:3) = [3.0d0, -4.0d0, 1.0d0]/2.0d0
    if ( this % max_bdf_order .ge. 3 ) this % A(3,1:4) = [11.0d0, -18.0d0, 9.0d0, -2.0d0]/6.0d0
    if ( this % max_bdf_order .ge. 4 ) this % A(4,1:5) = [25.0d0, -48.0d0, 36.0d0, -16.0d0, 3.0d0]/12.0d0
    if ( this % max_bdf_order .ge. 5 ) this % A(5,1:6) = [137.0d0, -300.0d0, 300.0d0, -200.0d0, 75.0d0, -12.0d0]/60.0d0
    if ( this % max_bdf_order .ge. 6 ) this % A(6,1:7) = [147.0d0, -360.0d0, 450.0d0, -400.0d0, 225.0d0, -72.0d0, 10.0d0]/60.0d0

    ! Sanity check on BDF coeffs
    sanity_check: block
      type(integer) :: j
      do j = 1, this % max_bdf_order
         if (abs(sum(real_part(this % A(j,1:j+1)))) .gt. 1.0d-15 ) then
            print *, "Error in BDF Coeff for order ", abs(sum(real_part(this % A(j,1:j+1)))), j
            stop
         end if
      end do
    end block sanity_check

    !-----------------------------------------------------------------!
    ! Setup adjoint RHS
    !-----------------------------------------------------------------!
    
    this % num_rhs_bins = 2*this % max_bdf_order + 1

    allocate(this % rhs(this % nsvars))
    this % rhs = 0.0d0
    
  end function initialize

  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!
  
  subroutine finalize( this )

    class(BDF) :: this
    
    if ( allocated(this % A) ) deallocate(this % A)
    if ( allocated(this % rhs) ) deallocate(this % rhs)

  end subroutine finalize

  !===================================================================!
  ! Subroutine that marches backwards in time to compute the lagrange
  ! multipliers (adjoint variables for the function)
  ! ===================================================================!
  
  subroutine marchBackwards( this )

    class(BDF)                :: this
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

       call this % adjointSolve(this % psi(k,:), alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))
       
       !print*, "k,psi=", k, this % psi(k,:)

    end do time
    
  end subroutine marchBackwards
  
  !===================================================================!
  ! Function that puts together the right hand side of the adjoint
  ! equation into the supplied rhs vector.
  !===================================================================!
  
  subroutine assembleRHS( this, rhs )

    class(BDF)                                :: this
    type(scalar), dimension(:), intent(inout) :: rhs
    type(scalar), dimension(:,:), allocatable :: jac
    type(scalar)                              :: scale
    integer                                   :: k, i, m1, m2, idx
    
    stop"Broke"
!!$
!!$    allocate( jac(this % nSVars, this % nSVars) )
!!$    
!!$    k = this % current_step 
!!$    m1 = this % coeff % (k, 1)
!!$    m2 = this % coeff % getOrder( k, 2)
!!$
!!$    if (m2 .eq. 0) m2 = 1
!!$
!!$    ! Zero the RHS first
!!$    rhs = 0.0d0
!!$    
!!$    !-----------------------------------------------------------------!
!!$    ! Add function contribution (dfdu)
!!$    !-----------------------------------------------------------------!
!!$    
!!$    call this % system % func % addDFdU(rhs, ONE, this % time(k), &
!!$         & this % system % x, this % u(k,:), this % udot(k,:), this % uddot(k,:))
!!$
!!$    do i = 0, m1
!!$       idx = k + i
!!$       if ( idx .le. this % num_steps) then
!!$          scale = this % coeff % beta(m1, i+1)/this % h
!!$          call this % system % func % addDFdUDot(rhs, scale, this % time(idx), &
!!$               & this % system % x, this % u(idx,:), this % udot(idx,:), this % uddot(idx,:))
!!$       end if
!!$    end do
!!$
!!$    do i = 0, 2*m2
!!$       idx = k + i
!!$       if ( idx .le. this % num_steps) then
!!$          scale = this % coeff % gamma(m2, i+1)/this % h/this % h
!!$          call this % system % func % addDFdUDDot(rhs, scale, this % time(idx), &
!!$               & this % system % x, this % u(idx,:), this % udot(idx,:), this % uddot(idx,:))
!!$       end if
!!$    end do
!!$
!!$    !-----------------------------------------------------------------!
!!$    ! Add residual contribution
!!$    !-----------------------------------------------------------------!
!!$    
!!$    do i = 1, m1 
!!$       idx = k + i
!!$       if ( idx .le. this % num_steps) then
!!$          scale = this % coeff % beta(m1, i+1)/this % h
!!$          call this % system % assembleJacobian(jac, ZERO, ONE, ZERO, &
!!$               & this % time(idx), this % u(idx,:), this % udot(idx,:), this % uddot(idx,:))
!!$          rhs = rhs + scale*matmul( transpose(jac), this % psi(idx,:) )
!!$       end if
!!$    end do
!!$
!!$    do i = 1, 2*m2
!!$       idx = k + i
!!$       if ( idx .le. this % num_steps) then
!!$          scale = this % coeff % gamma(m2, i+1)/this % h/this % h
!!$          call this % system % assembleJacobian(jac, ZERO, ZERO, ONE, &
!!$               & this % time(idx), this % u(idx,:), this % udot(idx,:), this % uddot(idx,:))
!!$          rhs = rhs + scale*matmul( transpose(jac), this % psi(idx,:) )
!!$       end if
!!$    end do
!!$    
    ! Negate the RHS
    rhs = - rhs
    
    if(allocated(jac)) deallocate(jac)
    
  end subroutine assembleRHS

end module bdf_integrator
