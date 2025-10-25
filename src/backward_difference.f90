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

  implicit none

  private
  public :: BDF

  !===================================================================!
  ! A derived type for the bdf coefficients
  !===================================================================!
  
  type :: bdf_coeff
     
     ! Information
     integer                               :: max_order = 3

     ! Coeff values
     type(scalar)                              :: alpha
     type(scalar), dimension(:,:), allocatable :: beta
     type(scalar), dimension(:,:), allocatable :: gamma

   contains

     private
     
     procedure :: destruct
     procedure :: getOrder
     
  end type bdf_coeff
  
  ! Interface for the constructor of bdf_coeff type
  interface bdf_coeff
     module procedure construct_bdf_coeff
  end interface bdf_coeff
  
  !===================================================================! 
  ! BDF Integrator type
  !===================================================================! 

  type, extends(integrator) :: BDF
     
     ! BDF variables
     integer            :: max_bdf_order = 3
     type(bdf_coeff)    :: coeff

     type(scalar), dimension(:,:), allocatable :: rhsbin

   contains

     ! Destructor
     procedure, public  :: finalize

     ! Routines for integration
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff

     ! Routines for adjoint gradient

     procedure, public  :: marchBackwards
     procedure, private :: evaluate_adjoint
     procedure, private :: distribute_contributions

     procedure, private :: assembleRHS
     procedure          :: computeTotalDerivative

  end type BDF  

  ! Constructor for BDF type
  interface BDF
     module procedure initialize
  end interface BDF

contains
  
  !===================================================================!
  ! Constructor for BDF coeff object
  !===================================================================!
  
  type(bdf_coeff) function construct_bdf_coeff( max_order ) result( this )

    integer, intent(in), OPTIONAL :: max_order
    
    if( present(max_order) ) then
       this % max_order = max_order
    end if

    allocate( this % beta (0:this % max_order, this % max_order + 1) )
    this % beta = 0.0d0 

    allocate( this % gamma (0:this % max_order, 2*this % max_order + 1) ) ! 2m + 1
    this % gamma = 0.0d0

    ! Set the coefficients

    this % alpha = 1.0d0
    this % gamma(0, 1:1) = 1.0d0

    if (this % max_order .eq. 3) then

       this % beta (1, 1:2) = (/ 1.0, -1.0 /)
       this % beta (2, 1:3) = (/ 1.5d0, -2.0d0, 0.5d0 /)
       this % beta (3, 1:4) = (/ 35.d0/24.0d0, -15.d0/8.0d0, 3.0d0/8.0d0, 1.0d0/24.0d0 /)

       this % gamma(1, 1:3) = (/ 1.0d0, -2.0d0, 1.0d0 /)
       this % gamma(2, 1:5) = (/ 2.25d0, -6.0d0, 5.5d0, -2.0d0, 0.25d0 /)
       this % gamma(3, 1:7) = (/ 2.126736d0, -5.468750d0, 4.609375d0, &
            & -1.284722d0, -0.015625d0, 0.031250d0, 0.001736d0 /)

    else if (this % max_order .eq. 2) then

       this % beta (1, 1:2) = (/ 1.0, -1.0 /)
       this % beta (2, 1:3) = (/ 1.5d0, -2.0d0, 0.5d0 /)

       this % gamma(1, 1:3) = (/ 1.0d0, -2.0d0, 1.0d0 /)
       this % gamma(2, 1:5) = (/ 2.25d0, -6.0d0, 5.5d0, -2.0d0, 0.25d0 /)

    else if (this % max_order .eq. 1) then

       this % beta (1, 1:2) = (/ 1.0, -1.0 /)
       this % gamma(1, 1:3) = (/ 1.0d0, -2.0d0, 1.0d0 /)
    else 
       print *,  "Wrong max_bdf_order:", this % max_order
       stop
    end if
    
  end function construct_bdf_coeff
  
  !===================================================================!
  ! Destructor for BDF coeff object
  !===================================================================!
  
  subroutine destruct( this ) 
    
    class(bdf_coeff) :: this
    
    if ( allocated(this % beta) ) deallocate(this % beta)
    if ( allocated(this % gamma) ) deallocate(this % gamma)

  end subroutine destruct
  
  !===================================================================!
  ! Returns the linearization scalar coefficients: alpha, beta, gamma
  !===================================================================!
  
  subroutine getLinearCoeff( this, alpha, beta, gamma )
    
    class(BDF) :: this
    type(scalar), intent(out)   :: alpha, beta, gamma
    integer :: k

    k = this % current_step
   
    if ( this % second_order ) then
       gamma = this % coeff % gamma(this % coeff % getOrder(k, 2), 1)/this % h/this % h
    else 
       gamma = 0.0d0
    end if
    beta  = this % coeff % beta(this % coeff % getOrder(k,1), 1)/this % h
    alpha = this % coeff % alpha
    
  end subroutine getLinearCoeff

  !===================================================================!
  ! Approximate the state variables at each step using BDF formulae
  !===================================================================!
  
  subroutine approximateStates( this )

    class(BDF) :: this
    integer    :: k, m, i

    k = this % current_step
    
    !-----------------------------------------------------------------!
    ! Extrapolate U to next time step
    !-----------------------------------------------------------------!

    this % u(k,:) = this % u(k-1,:) + this % udot(k-1,:)*this % h &
         & + this % uddot(k-1,:)*this % h*this % h/2.0d0

    !-----------------------------------------------------------------!
    ! Approximate UDOT using BDF
    !-----------------------------------------------------------------!
    
    m = this % coeff % getOrder(k, 1)

    do i = 1, m + 1
       this % udot(k,:) = this % udot(k,:) &
            & + this % coeff % beta(m, i)*this % u(k-i+1,:)/this % h
    end do
    
    !-----------------------------------------------------------------!
    ! Approximate UDDOT using BDF
    !-----------------------------------------------------------------!
    
    m = this % coeff % getOrder(k, 2)

    if ( m .eq. 0) then

       ! We dont have enought points yet
       this % uddot(k,:) = (this % udot(k,:) - this % udot(k-1,:))/this % h
       
    else
       
       do i = 1, 2*m + 1
          this % uddot(k,:) = this % uddot(k,:) &
               & + this % coeff % gamma(m, i)*this % u(k-i+1,:)/this % h/this % h
       end do
       
    end if
    
  end subroutine approximateStates

  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!
  
  pure integer function getOrder(this, k, d )

    class(bdf_coeff), intent(in)    :: this
    integer, intent(in)             :: k, d

    ! find the order of approximation
    getOrder = (k-1)/d  ! k = md + 1

    ! Do not let exceed the max order sought
    if ( getOrder .gt. this % max_order ) getOrder = this % max_order

  end function getOrder
 
  !===================================================================!
  ! Compute the total derivative of the function with respect to the
  ! design variables and return the gradient 
  !===================================================================!
  
  subroutine computeTotalDerivative( this, dfdx )
    
    class(BDF)                                 :: this
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
       dfdx = dfdx + matmul(transpose(dRdX),this % mu(k,:)) ! same as {mu}, [dRdX]
    end do

!!$    ! Add constraint contribution
!!$    call this % system % getResidualDVSens(dRdX, scale, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(1,:))
!!$    dfdx = dfdx + matmul(this % psi(2,:), dRdX)

    ! Finally multiply by the scalar
    dfdx = this % h * dfdx

    if (allocated(dRdX)) deallocate(dRdX)

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

    if (present(max_bdf_order)) then
       this % max_bdf_order = max_bdf_order
    end if
    print '("  >> Max BDF Order          : ",i4)', this % max_bdf_order

    !-----------------------------------------------------------------!
    ! Create the BDF-coeff object
    !-----------------------------------------------------------------!
    
    this % coeff = bdf_coeff(this % max_bdf_order)

    !-----------------------------------------------------------------!
    ! Setup adjoint RHS
    !-----------------------------------------------------------------!

    allocate(this % rhs(this % nsvars))
    this % rhs = 0.0d0

    allocate(this % rhsbin(this % num_steps, this % nsvars))
    this % rhsbin = 0.0d0
    
  end function initialize

  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!
  
  subroutine finalize( this )

    class(BDF) :: this

    call this % destruct()

    ! call the BDF coeff destructor
    call this % coeff % destruct()
    
    if (allocated(this % rhs)) deallocate(this % rhs)
    if (allocated(this % rhsbin)) deallocate(this % rhsbin)

  end subroutine finalize

  !===================================================================!
  ! Subroutine that marches backwards in time to compute the Lagrange
  ! multipliers (adjoint variables for the function)
  !===================================================================!
  
  subroutine marchBackwards( this )

    class(BDF)    :: this
    type(integer) :: k
    integer, parameter :: unit = 88

    open(unit=88, file=trim("bdf_adjoint.dat"))

    time: do k = this % num_steps, 2, -1

       this % current_step = k

       !--------------------------------------------------------------!
       ! Solve the adjoint equation at each step
       !--------------------------------------------------------------!
       
       call this % evaluate_adjoint(this % mu, this % psi, this % phi)

       !--------------------------------------------------------------!
       ! Drop the contributions from this step to corresponding bins
       !--------------------------------------------------------------!

       call this % distribute_contributions(this % rhsbin)


       write_output: block

         type(integer) :: j

         ! Write the solution as output
         write(unit, *)  this % time(k), &
              & (dble(this % mu(k,j)), j=1,this%nsvars ), &
              & (dble(this % psi(k,j)), j=1,this%nsvars ), &
              & (dble(this % phi(k,j)), j=1,this%nsvars )

       end block write_output

    end do time

    close(88)

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
    type(scalar)                              :: alpha, beta, gamma

    ! Retrieve the current time index
    k = this % current_step
    
    ! Replace from rhsbin for this step (contains k+ terms)
    rhs = this % rhsbin(k,:)
    
    ! Coefficients
    call this % getLinearCoeff(alpha, beta, gamma)

    ! Add the state variable sensitivity from the previous step
    call this % system % func % addFuncSVSens(rhs, &
         & alpha, beta, gamma, &
         & this % time(k), &
         & this % system % X, &
         & this % u(k,:), &
         & this % udot(k,:), &
         & this % uddot(k,:))

    ! Negate the RHS
    rhs = - rhs

  end subroutine assembleRHS

  !===================================================================!
  ! Evaluates the adjoint variable values at the current step
  !===================================================================!
  
  subroutine evaluate_adjoint(this, mu, psi, phi)

    class(BDF)                   :: this
    type(integer)                :: k
    type(scalar)                 :: alpha, beta, gamma
    type(scalar), dimension(:,:) :: mu, psi, phi

    ! Retrive the current step number
    k = this % current_step

    ! Evaluate PHI (no linear solution necessary)
    phi(k,:) = 0.0d0

    ! Evaluate PSI (no linear solution necessary)
    psi(k,:) = 0.0d0

    ! Evaluate MU
    call this % getLinearCoeff(alpha, beta, gamma)

    call this % adjointSolve(mu(k,:), &
         & alpha, beta, gamma, &
         & this % time(k), &
         & this % u(k,:), this % udot(k,:), this % uddot(k,:))

  end subroutine evaluate_adjoint

  !===================================================================!
  ! Add contributions from kth step to k- steps
  !===================================================================!

  subroutine distribute_contributions(this, rhsbin)

    class(BDF)                   :: this
    type(integer)                :: k, p1, p2, i
    type(scalar)                 :: alpha, beta, gamma
    type(scalar), dimension(:,:) :: rhsbin

    ! Retrive the current step number
    k = this % current_step

    ! Determine the order of approximation of differential states
    p1 = this % coeff % getOrder(k, 1)
    p2 = this % coeff % getOrder(k, 2)

    !-----------------------------------------------------------------!
    ! Add contributions from k to k- RHSBIN    
    !-----------------------------------------------------------------!

    do i = 2, p1+1

       ! Coefficients
       alpha = 0.0d0
       beta  = this % coeff % beta(p1, i)/this % h
       gamma = 0.0d0

       call this % addFuncResAdjPdt(rhsbin(k-i+1,:), &
            & alpha, beta, gamma, this % time(k), &
            & this % u(k,:), this % udot(k,:), this % uddot(k,:), &
            & this % mu(k,:))

    end do
    
    do i = 2, 2*p2+1

       ! Coefficients
       alpha = 0.0d0
       beta  = 0.0d0
       gamma = this % coeff % gamma(p2, i)/this % h/this % h

       call this % addFuncResAdjPdt(rhsbin(k-i+1,:), &
            & alpha, beta, gamma, this % time(k), &
            & this % u(k,:), this % udot(k,:), this % uddot(k,:), &
            & this % mu(k,:))

    end do

  end subroutine distribute_contributions

end module bdf_integrator
