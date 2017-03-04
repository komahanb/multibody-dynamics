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

  ! Useful constants
  type(scalar), parameter :: ONE  = 1.0d0
  type(scalar), parameter :: ZERO = 0.0d0

  !===================================================================!
  ! A derived type for the bdf coefficients
  !===================================================================!
  
  type :: bdf_coeff
     
     ! Information
     integer :: max_order = 3

     ! Coeff values
     type(scalar)                              :: alpha
     type(scalar), dimension(:,:), allocatable :: beta
     type(scalar), dimension(:,:), allocatable :: gamma

   contains

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

   contains

     ! Destructor
     procedure, public  :: finalize

     ! Routines for integration
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff

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
   
  end function initialize

  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!
  
  subroutine finalize( this )

    class(BDF) :: this

    ! call the BDF coeff destructor
    call this % coeff % destruct()
   
  end subroutine finalize

end module bdf_integrator
