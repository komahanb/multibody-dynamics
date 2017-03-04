#include "scalar.fpp"
!=====================================================================!
! Newmark-Beta-Gamma Integration Module for first and second order
! systems with adjoint derivative capabilities.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module nbg_integrator

  use integrator_class, only : integrator
  use physics_class,    only : physics

  implicit none

  private
  public :: NBG
  
  type(scalar), parameter :: ONE  = 1.0d0
  type(scalar), parameter :: ZERO = 0.0d0

  !===================================================================! 
  ! NBG Integrator type
  !===================================================================! 

  type, extends(integrator) :: NBG
     
     private

     ! Average Constant Accelearation (second order unconditionally stable)
     type(scalar) :: BETA   = 0.25d0
     type(scalar) :: GAMMA  = 0.50d0
     
   contains

     ! Destructor
     procedure, public  :: finalize
     
     ! Routines for integration
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff

  end type NBG

  interface NBG
     module procedure initialize
  end interface NBG

contains
  
  !===================================================================!
  ! Initialize the NBG datatype and allocate required variables
  !===================================================================!
  
  type(nbg) function initialize( system, tinit, tfinal, h, second_order )  result (this)

    class(physics), target          :: system
    real(dp) , OPTIONAL, intent(in) :: tinit, tfinal
    real(dp) , OPTIONAL, intent(in) :: h
    logical  , OPTIONAL, intent(in) :: second_order

    print *, "======================================"
    print *, ">>   Newmark Beta Gamma (NBG)      << "
    print *, "======================================"

    call this % construct(system, tinit, tfinal, h, second_order)

    if ( .not. this % second_order ) then

       print *, " Warning: Newmark-Beta-Gamma method works for second" &
            & // " order systems in current form..."

    end if

  end function initialize

  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!

  subroutine finalize( this )

    class(NBG) :: this

  end subroutine finalize

  !===================================================================!
  ! Approximate the state variables at each step using NBG formula
  !===================================================================!

  subroutine approximateStates( this )

    class(NBG) :: this
    integer                   :: k
    type(scalar)              :: scale

    k = this % current_step

    !-----------------------------------------------------------------!
    ! Assume a UDDOT for the next time step
    !-----------------------------------------------------------------!

    this % uddot(k,:) = this % uddot(k-1,:)

    !-----------------------------------------------------------------!
    ! Approximate UDOT using NBG
    !-----------------------------------------------------------------!

    this % udot(k,:) = this % udot(k-1,:) 

    scale = this % h * (1.0d0 - this % GAMMA)
    this % udot(k,:) = this % udot(k,:) + scale*this % uddot(k-1,:) 

    scale = this % h * this % GAMMA
    this % udot(k,:) = this % udot(k,:) + scale*this % uddot(k,:) 

    !-----------------------------------------------------------------!
    ! Approximate U using NBG
    !-----------------------------------------------------------------!

    this % u(k,:) = this % u(k-1,:) 

    scale = this % h
    this % u(k,:) = this % u(k,:) + scale*this % udot(k-1,:) 

    scale = this % h * this % h * (1.0d0 - 2.0d0 * this % BETA)/2.0d0
    this % u(k,:) = this % u(k,:) + scale*this % uddot(k-1,:) 

    scale = this % h * this % h * this % BETA
    this % u(k,:) = this % u(k,:) + scale*this % uddot(k,:) 

  end subroutine approximateStates

  !===================================================================!
  ! Returns the linearization scalar coefficients: alpha, beta, gamma
  !===================================================================!

  subroutine getLinearCoeff( this, alpha, beta, gamma )

    class(NBG) :: this
    type(scalar), intent(out) :: alpha, beta, gamma
    integer                     :: m

    if ( this % second_order ) then
       gamma = 1.0d0/this % h/ this % h
       beta  = this % GAMMA/ this % h
       alpha = this % BETA
    else
       gamma = 1.0d0/this % h/ this % h
       beta  = this % GAMMA/ this % h
       alpha = this % BETA
    end if

  end subroutine getLinearCoeff

end module nbg_integrator
