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
     procedure :: writeSolution => writeSolution
     
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
       dfdx = dfdx + matmul( transpose(dRdX), this % mu(k,:)) ! check order
    end do

!!$    ! Add constraint contribution
!!$    call this % system % getResidualDVSens(dRdX, scale, this % time(1), &
!!$         & this % system % x, this % u(1,:), this % udot(1,:), this % uddot(1,:))
!!$    dfdx = dfdx + matmul(this % mu(2,:), dRdX)

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

       call this % adjointSolve(this % mu(k,:), alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))
       
    end do time

    do k = 2, this % num_steps
       print *, real(this % mu(k,:)), real(this % psi(k,:)), real(this % phi(k,:))
    end do 
           
    
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
    type(scalar)                              :: scale
    integer                                   :: k, m
    type(scalar)                              :: alpha, beta, gamma
    
    ! Zero the RHS first
    rhs = 0.0d0

    ! Retrieve the current time index
    k = this % current_step
    m = this % getOrder(k)
    
    if ( k+1 .le. this % num_steps ) then
       
       ! Find PHI
       this % phi(k,:) = this % phi(k+1,:)

       gamma = 0.0d0
       beta  = 0.0d0
       alpha = this % h

       call this % addFuncResAdjPdt(this % phi(k,:), &
            & alpha, beta, gamma, this % time(k+1), &
            & this % u(k+1,:), this % udot(k+1,:), this % uddot(k+1,:), &
            & this % mu(k+1,:))

    end if

    if ( k+1 .le. this % num_steps ) then
       
       ! Find PSI
       this % psi(k,:) = this % psi(k+1,:)
       
       this % psi(k,:) = this % psi(k,:) &
            & + this % h * this % A(1, 1) * this % phi(k+1,:) ! check coeff

       gamma = 0.0d0
       beta  = this % h
       alpha = this % h * this % h * this % A(1, 1)
       
       call this % addFuncResAdjPdt(this % psi(k,:), alpha, beta, gamma, this % time(k+1), &
            & this % u(k+1,:), this % udot(k+1,:), this % uddot(k+1,:), &
            & this % mu(k+1,:))

    end if
    
    gamma = 1.0d0
    beta  = this % h * this % A(1,1)
    alpha = this % h * this % A(1,1) * this % h * this % A(1,1)

    ! Add the state variable sensitivity from the previous step
    call this % system % func % addFuncSVSens(rhs, &
         & alpha, beta, gamma, &
         & this % time(k), &
         & this % system % X, &
         & this % u(k,:), &
         & this % udot(k,:), &
         & this % uddot(k,:))

    rhs = rhs + this % A(1,1) * this % psi(k,:)
    rhs = rhs + this % A(1,1) * this % h * this % A(1,1) * this % phi(k,:)
    
    if ( k + 1 .le. this % num_steps ) then
!!$      
!!$       gamma = 0.0d0
!!$       if ( k .eq. 2) then
!!$          beta  = this % h * this % A(2,2)
!!$          alpha = this % h * this % A(2,2) * this % h * this % A(2,1) + &
!!$               & this % h * this % A(2,2) * this % h * this % A(1,1)
!!$          rhs = rhs + beta/ this % h * this % psi(k+1,:)
!!$          rhs = rhs + alpha/ this % h * this % phi(k+1,:) ! this term corrects the error from third step onwards
!!$       else
!!$          beta  = this % h * this % A(2,2)
!!$          alpha = this % h * this % A(2,2) * this % h * this % A(2,1) + &
!!$               & this % h * this % A(2,1) * this % h * this % A(2,2)
!!$          rhs = rhs + beta/ this % h * this % psi(k+1,:)
!!$          rhs = rhs + alpha/ this % h * this % phi(k+1,:) ! this term corrects the error from third step onwards
!!$       end if
!!$
!!$       call this % addFuncResAdjPdt(alpha, beta, gamma, this % time(k+1), &
!!$            & this % u(k+1,:), this % udot(k+1,:), this % uddot(k+1,:), &
!!$            & this % mu(k+1,:), rhs)
       
    end if
    
    ! Negate the RHS
    rhs = - rhs

  end subroutine assembleRHS

  
  !===================================================================!
  ! Write solution to file
  !===================================================================!

  subroutine writeSolution(this, filename)

    class(ABM)                      :: this
    character(len=*), OPTIONAL, intent(in) :: filename
    character(len=7), parameter            :: directory = "output/"
    character(len=32)                      :: path = ""
    integer                                :: k, j, ierr

    path = trim(path)

    if (present(filename)) then
       path = directory//filename
    else
       path = directory//"solution.dat"
    end if

    open(unit=90, file=trim(path), iostat= ierr)

    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') path
       return
    end if

    do k = 1, this % num_steps
       write(90, *)  this % time(k), &
            & (dble(this % u     (k,j) ), j=1,this%nsvars ), &
            & (dble(this % udot  (k,j) ), j=1,this%nsvars ), &
            & (dble(this % uddot (k,j) ), j=1,this%nsvars ), & 
            & (dble(this % phi   (k,j) ), j=1,this%nsvars ), &
            & (dble(this % psi   (k,j) ), j=1,this%nsvars ), &
            & (dble(this % mu   (k,j) ), j=1,this%nsvars )
    end do

    close(90)

  end subroutine writeSolution

  
end module abm_integrator
