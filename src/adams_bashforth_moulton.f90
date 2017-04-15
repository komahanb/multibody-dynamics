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
  use utils,            only : real_part

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
     integer :: max_abm_order = 6
     type(scalar), allocatable, dimension(:,:) :: A

     type(scalar), dimension(:,:), allocatable :: rhsbin
     type(scalar), dimension(:,:), allocatable :: psibin
     type(scalar), dimension(:,:), allocatable :: phibin

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

     procedure, private :: evaluate_adjoint
     procedure, private :: distribute_contributions
    
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
    integer :: j

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
    if ( this % max_abm_order .ge. 1 ) this % A(1,1:1) = [1.0d0]
    if ( this % max_abm_order .ge. 2 ) this % A(2,1:2) = [1.0d0, 1.0d0]/2.0d0
    if ( this % max_abm_order .ge. 3 ) this % A(3,1:3) = [5.0d0, 8.0d0, -1.0d0]/12.0d0
    if ( this % max_abm_order .ge. 4 ) this % A(4,1:4) = [9.0d0, 19.0d0, -5.0d0, 1.0d0]/24.0d0
    if ( this % max_abm_order .ge. 5 ) this % A(5,1:5) = [251.0d0, 646.0d0, -264.0d0, 106.0d0, -19.0d0]/720.0d0
    if ( this % max_abm_order .ge. 6 ) this % A(6,1:6) = [475.0d0, 1427.0d0, -798.0d0, 482.0d0, -173.0d0, 27.0d0]/1440.0d0

    ! Sanity check on ABM coeffs
    do j = 1, this % max_abm_order
       if ( sum(abs(real_part(this % A(j,1:j))) - 1.0d0) .gt. 0.00001 ) then
          stop "Error in ABM Coeff"
       end if
    end do

    !-----------------------------------------------------------------!
    ! Setup adjoint RHS
    !-----------------------------------------------------------------!   

    this % num_rhs_bins = this % max_abm_order

    allocate(this % rhs(this % nsvars))
    this % rhs = 0.0d0

    allocate(this % rhsbin(this % num_steps, this % nsvars))
    this % rhsbin = 0.0d0

    allocate(this % phibin(this % num_steps, this % nsvars))
    this % phibin = 0.0d0

    allocate(this % psibin(this % num_steps, this % nsvars))
    this % psibin = 0.0d0

  end function initialize

  !===================================================================!
  ! Deallocate the allocated variables
  !===================================================================!
  
  subroutine finalize( this )

    class(ABM) :: this

    ! Deallocate ABM coefficient
    if(allocated(this % A)) deallocate(this % A)
    
    ! Adjoint variables and RHS
    if(allocated(this % rhs)) deallocate(this % rhs)
    if(allocated(this % rhsbin)) deallocate(this % rhsbin)
    if(allocated(this % psibin)) deallocate(this % psibin)
    if(allocated(this % phibin)) deallocate(this % phibin)

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

    class(ABM)    :: this
    type(integer) :: k

    time: do k = this % num_steps, 2, -1

       this % current_step = k

       !--------------------------------------------------------------!
       ! Solve the adjoint equation at each step
       !--------------------------------------------------------------!
       
       call this % evaluate_adjoint(this % mu, this % psi, this % phi)
       
       !--------------------------------------------------------------!
       ! Drop the contributions from this step to corresponding bins
       !--------------------------------------------------------------!

       call this % distribute_contributions(this % rhsbin, this % psibin, this % phibin)

    end do time

  end subroutine marchBackwards
  
  !===================================================================!
  ! Evaluates the adjoint variable values at the current step
  !===================================================================!
  
  subroutine evaluate_adjoint(this, mu, psi, phi)

    class(ABM)                   :: this
    type(integer)                :: k
    type(scalar)                 :: alpha, beta, gamma
    type(scalar), dimension(:,:) :: mu, psi, phi

    ! Retrive the current step number
    k = this % current_step

    ! Evaluate PHI (no linear solution necessary)
    phi(k,:) = this % phibin(k,:)

    ! Evaluate PSI (no linear solution necessary)
    psi(k,:) = this % psibin(k,:)

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

  subroutine distribute_contributions(this, rhsbin, psibin, phibin)

    class(ABM)                   :: this
    type(integer)                :: k, p, i
    type(scalar)                 :: alpha, beta, gamma
    type(scalar), dimension(:,:) :: rhsbin, phibin, psibin

    ! Retrive the current step number
    k = this % current_step
    p = this % getOrder(k)

    !-----------------------------------------------------------------!
    ! Add contributions from k to k- PHIBIN
    !-----------------------------------------------------------------!
    
    phibin(k-1,:) = phibin(k-1,:) + this % phi(k,:)

    gamma = 0.0d0
    beta  = 0.0d0
    alpha = this % h

    call this % addFuncResAdjPdt(phibin(k-1,:), &
         & alpha, beta, gamma, this % time(k), &
         & this % u(k,:), this % udot(k,:), this % uddot(k,:), &
         & this % mu(k,:))

    !-----------------------------------------------------------------!
    ! Add contributions from k to k- PSIBIN    
    !-----------------------------------------------------------------!

    psibin(k-1,:) = psibin(k-1,:) + this % psi(k,:)
    
    gamma = 0.0d0
    beta  = this % h 
    alpha = this % h * this % h * this % A(p,1)
    
    call this % addFuncResAdjPdt(psibin(k-1,:), &
         & alpha, beta, gamma, this % time(k), &
         & this % u(k,:), this % udot(k,:), this % uddot(k,:), &
         & this % mu(k,:))

    psibin(k-1,:) = psibin(k-1,:) + alpha/this%h*this%phi(k,:)

    do i = 2, p

       gamma = 0.0d0
       beta  = 0.0d0
       alpha = this % h * this % h * this % A(p,i)
       
       call this % addFuncResAdjPdt(psibin(k-i+1,:), &
            & alpha, beta, gamma, this % time(k), &
            & this % u(k,:), this % udot(k,:), this % uddot(k,:), &
            & this % mu(k,:))

       psibin(k-i+1,:) = psibin(k-i+1,:) + alpha/this%h*this % phi(k,:)

    end do

    !-----------------------------------------------------------------!
    ! Add contributions from k to k- RHSBIN    
    !-----------------------------------------------------------------!
    
    do i = 2, p

       gamma = 0.0d0
       beta  = this % h * this % A(p,i)
       alpha = this % h * this % A(p,1) * this % h * this % A(p,i)

       call this % addFuncResAdjPdt(rhsbin(k-i+1,:), &
            & alpha, beta, gamma, this % time(k), &
            & this % u(k,:), this % udot(k,:), this % uddot(k,:), &
            & this % mu(k,:))

       rhsbin(k-i+1,:) = rhsbin(k-i+1,:) + beta/this % h*this % psi(k,:)
       rhsbin(k-i+1,:) = rhsbin(k-i+1,:) + alpha/this % h*this % phi(k,:)

    end do

  end subroutine distribute_contributions

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
  
  subroutine assembleRHS(this, rhs)

    class(ABM)                                :: this
    type(scalar), dimension(:), intent(inout) :: rhs
    type(scalar)                              :: scale
    integer                                   :: k, m
    type(scalar)                              :: alpha, beta, gamma

    ! Retrieve the current time index
    k = this % current_step
    m = this % getOrder(k)

    ! Replace from rhsbin for this step (contains k+ terms)
    rhs = this % rhsbin(k,:)

    ! Calculate and add terms from k-th step
    call this % getLinearCoeff(alpha, beta, gamma)

    ! Add the state variable sensitivity from the previous step
    call this % system % func % addFuncSVSens(rhs, &
         & alpha, beta, gamma, &
         & this % time(k), &
         & this % system % X, &
         & this % u(k,:), &
         & this % udot(k,:), &
         & this % uddot(k,:))

    rhs = rhs + beta/this % h * this % psi(k,:)
    rhs = rhs + alpha/this % h * this % phi(k,:)

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
