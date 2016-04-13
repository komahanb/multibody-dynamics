!=====================================================================!
! Module handling Euler angle rotation for rigid body dynamics
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module rotation

  use utils

  implicit none

  private

  public :: get_rotation, get_angrate
  public :: get_angrate_dot, get_angrate_inv

  !-------------------------------------------------------------------!
  ! A common interface for different ways of getting rotation matrix
  !-------------------------------------------------------------------!
  ! (a) get_rotation_from_angles_vec   -- > input theta as VECTOR
  ! (b) get_rotation_from_angles_array -- > input theta(3) as array
  ! (c) get_rotation_from_cosines -- > input dir cosines and sines
  !-------------------------------------------------------------------!

  interface get_rotation
     module procedure get_rotation_from_angles_array, &
          & get_rotation_from_angles_vec, get_rotation_from_cosines
  end interface get_rotation

  !-------------------------------------------------------------------!
  ! A common interface for different ways of getting ang rate  matrix
  !-------------------------------------------------------------------!
  ! (a) get_angrate_from_angles_vec    -- > input VECTOR theta
  ! (b) get_angrate_from_angles_array  -- > input theta(3)
  ! (c) get_angrate_from_cosines       -- > input dir cosines and sines
  !-------------------------------------------------------------------!

  interface get_angrate
     module procedure get_angrate_from_angles_vec, &
          & get_angrate_from_angles_array, get_angrate_from_cosines
  end interface get_angrate

  !-------------------------------------------------------------------!
  ! A common interface for different ways of getting the time 
  ! derivative of the ang rate matrix
  !-------------------------------------------------------------------!
  ! (a) get_angrate_dot_vec -- > input VECTOR theta, theta_dot
  ! (b) get_angrate_dot_array -- > input theta(3), theta_dot(3)
  ! (c) get_angrate_dot_cosines  -- > input dir cosines and sines
  !-------------------------------------------------------------------!

  interface get_angrate_dot
     module procedure  get_angrate_dot_array, &
          & get_angrate_dot_vec, get_angrate_dot_cosines
  end interface get_angrate_dot

  !-------------------------------------------------------------------!
  ! A common interface for different ways of getting the inverse
  ! of the ang rate matrix
  !-------------------------------------------------------------------!
  ! (a) get_angrate_inv_vec -- > input VECTOR theta, theta_dot
  ! (b) get_angrate_inv_array -- > input theta(3), theta_dot(3)
  ! (c) get_angrate_inv_cosines  -- > input dir cosines and sines
  !-------------------------------------------------------------------!

  interface get_angrate_inv
     module procedure get_angrate_inv_vec, &
          & get_angrate_inv_array, &
          & get_angrate_inv_cosines
  end interface get_angrate_inv

contains

  !--------------------------------------------------------------------!
  ! Returns the rotation matrix based on the euler angles
  ! Compute the 3-2-1 Euler angle rotation

  ! C = C1(theta_1)*C2(theta_2)*C3(theta_3)
  !
  ! Input: theta as an array
  ! Output: CMAT of type MATRIX
  ! 
  ! Ref: Section 2.2 Eq. 18/19 Hughes
  !--------------------------------------------------------------------!

  pure function get_rotation_from_angles_array(theta) result(CMAT)

    real(dp), intent(in) :: theta(NUM_SPAT_DIM)    
    type(matrix)         :: CMAT
    real(dp)             :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles
    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    CMAT = get_rotation_from_cosines(c1, c2, c3, s1, s2, s3)

  end function get_rotation_from_angles_array

  !--------------------------------------------------------------------!
  ! Returns the rotation matrix based on the euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation
  !
  ! CMAT = C1(theta_1)*C2(theta_2)*C3(theta_3)
  !
  ! Input: theta of type VECTOR
  ! Output: CMAT of type MATRIX
  ! 
  ! Ref: Section 2.2 Eq. 18/19 Hughes
  !--------------------------------------------------------------------!

  pure function get_rotation_from_angles_vec(thetain) result(CMAT)

    type(vector), intent(in) :: thetain
    real(dp)                 :: theta(NUM_SPAT_DIM)
    type(matrix)             :: CMAT

    ! covert to array form
    theta = array(thetain)

    ! call the method that takes angles array
    CMAT  =  get_rotation_from_angles_array(theta)

  end function get_rotation_from_angles_vec

  !--------------------------------------------------------------------!
  ! Returns the rotation matrix (euler angles) based on the sines and 
  ! cosines of the euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation
  !
  ! Input: sines and cosines of the euler angles
  ! Output: CMAT of type MATRIX
  !
  ! Ref: Section 2.2 Eq. 18/19 Hughes
  !--------------------------------------------------------------------!

  pure function get_rotation_from_cosines(c1,c2,c3,s1,s2,s3) result(CMAT)

    real(dp), intent(in) :: c1, c2, c3, s1, s2, s3
    type(matrix)         :: CMAT

    CMAT = matrix((/ c2*c3, c2*s3, -s2,&
         & s1*s2*c3 - c1*s3, s1*s2*s3 + c1*c3, s1*c2,&
         & c1*s2*c3 + s1*s3, c1*s2*s3 - s1*c3, c1*c2 /))

  end function get_rotation_from_cosines

  !--------------------------------------------------------------------!
  ! Returns the ang rate matrix from the supplied euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation. The S matrix does not 
  ! depend on c3/s3, however we keep these as inputs in case we ever
  ! want to change the Euler sequence.
  !
  ! Input : theta of type VECTOR
  ! Output: SMAT  of type MATRIX
  ! 
  ! Ref: Section 2.3 Eq. 24/25 Hughes
  !--------------------------------------------------------------------!

  pure function get_angrate_from_angles_vec(thetain) result(SMAT)

    type(vector), intent(in) :: thetain
    real(dp)     :: theta(NUM_SPAT_DIM)    
    type(matrix) :: SMAT

    ! convert to array
    theta = array(thetain)

    ! call the function with array signature
    SMAT = get_angrate_from_angles_array(theta)

  end function get_angrate_from_angles_vec

  !--------------------------------------------------------------------!
  ! Returns the ang rate matrix from the supplied euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation. The S matrix does not 
  ! depend on c3/s3, however we keep these as inputs in case we ever
  ! want to change the Euler sequence.
  !
  ! Input : theta as an array
  ! Output: SMAT  of type MATRIX
  ! 
  ! Ref: Section 2.3 Eq. 24/25 Hughes
  !--------------------------------------------------------------------!

  pure function get_angrate_from_angles_array(theta) result(SMAT)

    real(dp), intent(in) :: theta(NUM_SPAT_DIM)    
    type(matrix) :: SMAT
    real(dp)     :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles
    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    SMAT = get_angrate_from_cosines(c1, c2, c3, s1, s2, s3)

  end function get_angrate_from_angles_array

  !--------------------------------------------------------------------!
  ! Returns the rotation matrix (euler angles) based on the sines and 
  ! cosines of the euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation. The S matrix does not 
  ! depend on c3/s3, however we keep these as inputs in case we ever
  ! want to change the Euler sequence.
  !
  ! Input : sines and cosines of euler angles
  ! Output: SMAT of type MATRIX
  ! 
  ! Ref: Section 2.3 Eq. 24/25 Hughes
  !--------------------------------------------------------------------!

  pure function get_angrate_from_cosines(c1,c2,c3,s1,s2,s3) result(SMAT)

    real(dp), intent(in) :: c1, c2, c3, s1, s2, s3
    type(matrix)         :: SMAT

    SMAT = matrix((/ 1.0_dp, 0.0_dp, -s2, &
         & 0.0_dp,  c1,  s1*c2, &
         & 0.0_dp,  -s1,  c1*c2 /))

  end function get_angrate_from_cosines

  !--------------------------------------------------------------------!
  ! Returns the time derivative of angular rate matrix
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : the euler angles and time derivatives as VECTOR
  ! Output: SMAT_DOT
  !--------------------------------------------------------------------!

  pure function get_angrate_dot_vec(thetain, dthetain) result(SMAT_DOT)

    type(vector), intent(in)  :: thetain, dthetain
    type(matrix) :: SMAT_DOT

    real(dp)     :: theta(NUM_SPAT_DIM), dtheta(NUM_SPAT_DIM)

    ! convert vec to array
    theta = array(thetain); dtheta=array(dthetain);   

    ! call the function matching array signature
    SMAT_DOT = get_angrate_dot_array(theta,dtheta)

  end function get_angrate_dot_vec

  !--------------------------------------------------------------------!
  ! Returns the time derivative of angular rate matrix
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : the euler angles and time derivatives as arrays
  ! Output: SMAT_DOT
  !--------------------------------------------------------------------!

  pure function get_angrate_dot_array(theta, dtheta) result(SMAT_DOT)

    real(dp), intent(in) :: theta(NUM_SPAT_DIM), dtheta(NUM_SPAT_DIM)
    type(matrix) :: SMAT_DOT
    real(dp)     :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles

    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    SMAT_DOT = get_angrate_dot_cosines( c1, c2, c3, s1, s2, s3, dtheta)

  end function get_angrate_dot_array

  !--------------------------------------------------------------------!
  ! Returns the time derivative of angular rate matrix
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : The sines and cosines of euler angles
  ! Output: SMAT_DOT
  !--------------------------------------------------------------------!

  pure function get_angrate_dot_cosines( c1, c2, c3, s1, s2, s3, dtheta) &
       & result(SMAT_DOT)

    real(dp), intent(in) :: dtheta(NUM_SPAT_DIM) 
    real(dp), intent(in) :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SMAT_DOT

    SMAT_DOT= matrix( (/ 0.0_dp,  0.0_dp,  -c2*dtheta(2), &
         & 0.0_dp, -s1*dtheta(1), c1*c2*dtheta(1)-s1*s2*dtheta(2),&
         & 0.0_dp, -c1*dtheta(1), -s1*c2*dtheta(1)-c1*s2*dtheta(2)/))

  end function get_angrate_dot_cosines


  !--------------------------------------------------------------------!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! theta vector
  !--------------------------------------------------------------------!

  pure function get_angrate_inv_vec(thetain) result(SMAT_INV)

    type(vector), intent(in) :: thetain
    type(matrix) :: SMAT_INV
    real(dp)     :: theta(NUM_SPAT_DIM)

    ! decompose the vector into array
    theta = array(thetain)

    ! call the method that takes array as input
    SMAT_INV =  get_angrate_inv_array(theta)

  end function get_angrate_inv_vec

  !--------------------------------------------------------------------!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! theta array
  !--------------------------------------------------------------------!

  pure function get_angrate_inv_array(theta) result(SMAT_INV)

    real(dp), intent(in) :: theta(NUM_SPAT_DIM)
    real(dp)     :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SMAT_INV

    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    SMAT_INV =  get_angrate_inv_cosines(  c1, c2, c3, s1, s2, s3)

  end function get_angrate_inv_array

  !--------------------------------------------------------------------!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! direction cosines
  !--------------------------------------------------------------------!

  pure function get_angrate_inv_cosines( c1, c2, c3, s1, s2, s3) result(SMAT_INV)

    real(dp), intent(in) :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SMAT_INV

    SMAT_INV= matrix( (/  1.0_dp, s1*s2/c2, c1*s2/c2, &
         &0.0_dp, c1, -s1, 0.0_dp,&
         & s1/c2, c1/c2 /))

  end function get_angrate_inv_cosines

end module rotation

!=====================================================================!
! Module that contains all the implementation of rigid body dynamics
! related physics. 
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module rigid_body_dynamics

  use utils
  use rotation

  implicit none

  private

  public :: rigid_body

  !-------------------------------------------------------------------!
  ! Rigid body type that contains pertaining data and routines that
  ! operate over the type variables
  !-------------------------------------------------------------------!

  type :: rigid_body

   contains

     procedure :: get_residual
     procedure :: get_jacobian

  end type rigid_body

contains

  !-------------------------------------------------------------------!
  ! Residual of the kinematic and dynamic equations in state-space
  ! representation
  !-------------------------------------------------------------------!

  subroutine get_residual(this)

    class(rigid_body) :: this

  end subroutine get_residual

  !-------------------------------------------------------------------!
  ! Jacobian of the kinematic and dynamic equations in state-space
  ! representation
  !-------------------------------------------------------------------!

  subroutine get_jacobian(this)

    class(rigid_body) :: this

  end subroutine get_jacobian

end module rigid_body_dynamics
