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
  ! A common interface for different ways of getting ROTATION matrix
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
  ! A common interface for different ways of getting DERIVATIVES of
  ! ROTATION matrix with respect to the  rotation parameters (thetas)
  ! -------------------------------------------------------------------!
  ! (a) get_rotation_dot_from_angles_vec -- > input theta as VECTOR 
  ! (b) get_rotation_dot_from_angles_array -- > input theta(3) as array 
  ! (c) get_rotation_dot_from_cosines -- > input dir cosines and sines
  ! -------------------------------------------------------------------!

  interface get_rotation_dot
     module procedure  get_rotation_dot_array, &
          & get_rotation_dot_vec, get_rotation_dot_cosines
  end interface get_rotation_dot

  !-------------------------------------------------------------------!
  ! A common interface for different ways of getting ANGULAR RATE matrix
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
  ! DERIVATIVE of the ANGULAR RATE matrix
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
  ! Returns the rotation transformation matrix based on the euler
  ! angles for 3-2-1 Euler sequence of rotation
  
  ! C = C1(theta_1)*C2(theta_2)*C3(theta_3)
  !
  ! Input: theta as an array
  ! Output: TMAT of type MATRIX
  ! 
  ! Ref: Section 2.2 Eq. 18/19 Hughes
  !--------------------------------------------------------------------!

  pure function get_rotation_from_angles_array(theta) result(TMAT)

    real(dp), intent(in) :: theta(NUM_SPAT_DIM)    
    type(matrix)         :: TMAT
    real(dp)             :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles
    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    TMAT = get_rotation_from_cosines(c1, c2, c3, s1, s2, s3)

  end function get_rotation_from_angles_array

  !--------------------------------------------------------------------!
  ! Returns the rotation matrix based on the euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation
  !
  ! TMAT = C1(theta_1)*C2(theta_2)*C3(theta_3)
  !
  ! Input: theta of type VECTOR
  ! Output: TMAT of type MATRIX
  ! 
  ! Ref: Section 2.2 Eq. 18/19 Hughes
  !--------------------------------------------------------------------!

  pure function get_rotation_from_angles_vec(thetain) result(TMAT)

    type(vector), intent(in) :: thetain
    real(dp)                 :: theta(NUM_SPAT_DIM)
    type(matrix)             :: TMAT

    ! covert to array form
    theta = array(thetain)

    ! call the method that takes angles array
    TMAT  =  get_rotation_from_angles_array(theta)

  end function get_rotation_from_angles_vec

  !--------------------------------------------------------------------!
  ! Returns the rotation matrix (euler angles) based on the sines and 
  ! cosines of the euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation
  !
  ! Input: sines and cosines of the euler angles
  ! Output: TMAT of type MATRIX
  !
  ! Ref: Section 2.2 Eq. 18/19 Hughes
  !--------------------------------------------------------------------!

  pure function get_rotation_from_cosines(c1,c2,c3,s1,s2,s3) result(TMAT)

    real(dp), intent(in) :: c1, c2, c3, s1, s2, s3
    type(matrix)         :: TMAT

    TMAT = matrix((/ c2*c3, c2*s3, -s2,&
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
  ! Output: SDOT
  !--------------------------------------------------------------------!

  pure function get_angrate_dot_vec(thetain, dthetain) result(SDOT)

    type(vector), intent(in)  :: thetain, dthetain
    type(matrix) :: SDOT

    real(dp)     :: theta(NUM_SPAT_DIM), dtheta(NUM_SPAT_DIM)

    ! convert vec to array
    theta = array(thetain); dtheta=array(dthetain);   

    ! call the function matching array signature
    SDOT = get_angrate_dot_array(theta,dtheta)

  end function get_angrate_dot_vec

  !--------------------------------------------------------------------!
  ! Returns the time derivative of angular rate matrix
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : the euler angles and time derivatives as arrays
  ! Output: SDOT
  !--------------------------------------------------------------------!

  pure function get_angrate_dot_array(theta, dtheta) result(SDOT)

    real(dp), intent(in) :: theta(NUM_SPAT_DIM), dtheta(NUM_SPAT_DIM)
    type(matrix) :: SDOT
    real(dp)     :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles

    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    SDOT = get_angrate_dot_cosines( c1, c2, c3, s1, s2, s3, dtheta)

  end function get_angrate_dot_array

  !--------------------------------------------------------------------!
  ! Returns the time derivative of angular rate matrix
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : The sines and cosines of euler angles
  ! Output: SDOT
  !--------------------------------------------------------------------!

  pure function get_angrate_dot_cosines( c1, c2, c3, s1, s2, s3, dtheta) &
       & result(SDOT)

    real(dp), intent(in) :: dtheta(NUM_SPAT_DIM) 
    real(dp), intent(in) :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SDOT
    type(matrix) :: TMP(3)
    
!!$    SDOT= matrix( (/ 0.0_dp,  0.0_dp,  -c2*dtheta(2), &
!!$         & 0.0_dp, -s1*dtheta(1), c1*c2*dtheta(1) - s1*s2*dtheta(2), &
!!$         & 0.0_dp, -c1*dtheta(1), -s1*c2*dtheta(1) - c1*s2*dtheta(2) /))
    
    TMP(1) = matrix( (/ 0.0_dp,  0.0_dp,  0.0_dp, &
         & 0.0_dp, -s1, c1*c2,&
         & 0.0_dp, -c1, -s1*c2/) )

    TMP(2) = matrix( (/ 0.0_dp,  0.0_dp,  -c2, &
         & 0.0_dp, 0.0_dp, -s1*s2,&
         & 0.0_dp, 0.0_dp, -c1*s2/))
    
    SDOT = dtheta(1)*TMP(1) + dtheta(2)*TMP(2) + dtheta(3)*TMP(3)
      
  end function get_angrate_dot_cosines

  !--------------------------------------------------------------------!
  ! Returns the time derivative of the transformation matrix with
  ! respect to the euler angles
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : The sines and cosines of euler angles
  ! Output: TDOT
  !--------------------------------------------------------------------!

  pure function get_rotation_dot_cosines( c1, c2, c3, s1, s2, s3, dtheta) &
       & result(TDOT)
    
    real(dp), intent(in) :: dtheta(NUM_SPAT_DIM) 
    real(dp), intent(in) :: c1, c2, c3, s1, s2, s3
    type(matrix) :: TDOT, TMP(3)

    TMP(1) = matrix( (/ 0.0_dp, 0.0_dp, 0.0_dp, &
         &  c1*s2*c3 + s1*s3, c1*s2*s3 - s1*c3, c1*c2, &
         & -s1*s2*c3 + c1*s3, -s1*s2*s3 - c1*c3, -s1*c2 /) )

    TMP(2) = matrix( (/ -s2*c3, -s2*s3, -c2, &
         & s1*c2*c3, s1*c2*s3, -s1*s2, &
         & c1*c2*c3, c1*c2*s3, -c1*s2 /))

    TMP(3) = matrix( (/ -c2*s3, c2*c3, 0.0_dp, &
         & -s1*s2*s3 - c1*c3, s1*s2*c3 - c1*s3, 0.0_dp, &
         & -c1*s2*s3 + s1*c3, c1*s2*c3 + s1*s3, 0.0_dp /) ) 

    TDOT = dtheta(1)*TMP(1) + dtheta(2)*TMP(2) + dtheta(3)*TMP(3)
    
  end function get_rotation_dot_cosines

  !--------------------------------------------------------------------!
  ! Returns the time derivative of transformation matrix
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : the euler angles and time derivatives as VECTOR
  ! Output: TDOT
  !--------------------------------------------------------------------!

  pure function get_rotation_dot_vec(thetain, dthetain) result(TDOT)

    type(vector), intent(in)  :: thetain, dthetain
    type(matrix) :: TDOT

    real(dp)     :: theta(NUM_SPAT_DIM), dtheta(NUM_SPAT_DIM)

    ! convert vec to array
    theta = array(thetain); dtheta=array(dthetain);   

    ! call the function matching array signature
    TDOT = get_rotation_dot_array(theta,dtheta)

  end function get_rotation_dot_vec

  !--------------------------------------------------------------------!
  ! Returns the time derivative of transformation matrix
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : the euler angles and time derivatives as arrays
  ! Output: TDOT
  !--------------------------------------------------------------------!

  pure function get_rotation_dot_array(theta, dtheta) result(TDOT)

    real(dp), intent(in) :: theta(NUM_SPAT_DIM), dtheta(NUM_SPAT_DIM)
    type(matrix) :: TDOT
    real(dp)     :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles

    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    TDOT = get_rotation_dot_cosines( c1, c2, c3, s1, s2, s3, dtheta)

  end function get_rotation_dot_array

  !--------------------------------------------------------------------!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! theta vector
  !--------------------------------------------------------------------!

  pure function get_angrate_inv_vec(thetain) result(SINV)

    type(vector), intent(in) :: thetain
    type(matrix) :: SINV
    real(dp)     :: theta(NUM_SPAT_DIM)

    ! decompose the vector into array
    theta = array(thetain)

    ! call the method that takes array as input
    SINV =  get_angrate_inv_array(theta)

  end function get_angrate_inv_vec

  !--------------------------------------------------------------------!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! theta array
  !--------------------------------------------------------------------!

  pure function get_angrate_inv_array(theta) result(SINV)

    real(dp), intent(in) :: theta(NUM_SPAT_DIM)
    real(dp)     :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SINV

    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    SINV =  get_angrate_inv_cosines(  c1, c2, c3, s1, s2, s3)

  end function get_angrate_inv_array

  !--------------------------------------------------------------------!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! direction cosines
  !--------------------------------------------------------------------!

  pure function get_angrate_inv_cosines( c1, c2, c3, s1, s2, s3) result(SINV)

    real(dp), intent(in) :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SINV

    SINV = matrix( (/  1.0_dp, s1*s2/c2, c1*s2/c2, &
         & 0.0_dp, c1, -s1, 0.0_dp,&
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

  use utils, only : vector, matrix, skew, unskew, &
       & operator(*), operator(+), operator(-), &
       & norm, array, dp

  use rotation

  implicit none

  private

  public :: rigid_body

  !-------------------------------------------------------------------!
  ! Rigid body type that contains pertaining data and routines that
  ! operate over the type variables
  ! 
  ! Usage: After creating this body, be sure to call set on
  ! this body with mandatory parameters
  ! 
  ! type(rigid_body) :: bodyA
  ! .
  ! .
  ! bodyA % set(...)
  ! -------------------------------------------------------------------!

  type :: rigid_body

     !----------------------------------------------------------------!
     ! Position state variables (all in inertial frame)
     !----------------------------------------------------------------!

     type(vector) :: r              
     type(vector) :: theta          
     type(vector) :: v              
     type(vector) :: omega          

     !----------------------------------------------------------------!
     ! Velocity state variables (all in inertial frame)
     !----------------------------------------------------------------!

     type(vector) :: rdot
     type(vector) :: thetadot
     type(vector) :: vdot
     type(vector) :: omegadot

     !----------------------------------------------------------------!
     ! Inertial Properties (all in inertial frame)
     !----------------------------------------------------------------!
     ! First moment of mass
     ! C = [ cx,  cy,  cz ]
     !----------------------------------------------------------------!
     ! Second moment of mass
     ! J = [ Jxx,  Jxy,  Jxz ] = [ J[0],  J[1],  J[2] ]
     ! . = [    ,  Jyy,  Jyz ] = [     ,  J[3],  J[4] ]
     ! . = [    ,     ,  Jzz ] = [     ,      ,  J[5] ]
     !----------------------------------------------------------------!

     real(dp)     :: mass    ! mass of the body
     type(vector) :: c       ! first moment of inertia
     type(matrix) :: J       ! second moment of inertia

     !----------------------------------------------------------------!
     ! Rotation matrices (rotates from body to inertial frame)
     !----------------------------------------------------------------!

     type(matrix) :: TIB, TIBdot ! Transformation matrix and its derivatives
     type(matrix) :: S, SDOT     ! angular velocity matrix and its derivative

     !----------------------------------------------------------------!
     ! Handy force and torque vectors
     !----------------------------------------------------------------!

     type(vector) :: rforce   ! external/reaction force
     type(vector) :: rtorque  ! external/reaction torque
     type(vector) :: grav     ! gravity vector in local frame

     !----------------------------------------------------------------!
     ! Used for energy conservation check
     !----------------------------------------------------------------!

     real(dp) :: KE          ! kinetic energy of the body
     real(dp) :: PE          ! potential energy of the body
     real(dp) :: TE          ! total energy of the body

   contains

     procedure :: set
     procedure :: get_residual 
     procedure :: get_jacobian

  end type rigid_body

contains

  !-------------------------------------------------------------------!
  ! Routine to set the states into the body and compute rotation
  ! matrices for the set state. This function is to be called every
  ! time the state of the body changes over time
  !-------------------------------------------------------------------!

  subroutine set(body, mass, q, qdot, qddot)

    class(rigid_body)    :: body
    real(dp), intent(in) :: mass
    real(dp), intent(in) :: q(12), qdot(12), qddot(12)

    ! set inertial properties
    body % mass    = mass
    body % c % x   = 0.0d0      ! Assuming body frame is located at CG
    body % J % PSI  = mass/6.0d0 ! Assuming atleast two planes of symmetry and a cube side = 1

    ! set the state into the body
    body % r        = vector(q(1:3))
    body % theta    = vector(q(4:6))
    body % v        = vector(q(7:9))
    body % omega    = vector(q(10:12))

    ! set the time derivatives of state into the body
    body % rdot     = vector(qdot(1:3))
    body % thetadot = vector(qdot(4:6))
    body % vdot     = vector(qdot(7:9))
    body % omegadot = vector(qdot(10:12))

    ! get rotation and angular rate matrices based on theta    
    body % TIB     = get_rotation(body % theta)
    body % S       = get_angrate(body % theta)
    body % SDOT    = get_angrate_dot(body % theta, body % thetadot)

  end subroutine set

  !-------------------------------------------------------------------!
  ! Residual of the kinematic and dynamic equations in state-space
  ! representation. The vector form of the residual 'R' can be
  ! converted into scalar form of 12 equations just by using
  ! 'array(R)'
  !-------------------------------------------------------------------!

  function get_residual(body) result(R)

    class(rigid_body) :: body
    type(vector) :: R(4)

    !-----------------------------------------------------------------!
    ! Kinematics eqn-1 (2 terms)
    !-----------------------------------------------------------------!
    ! [T] r_dot - v = 0
    !-----------------------------------------------------------------!

    !! May have to transform things into a frame
    R(1)  = body % TIB * body % rdot - body % v

    !-----------------------------------------------------------------!
    ! Kinematics eqn-2 (2 terms)
    !-----------------------------------------------------------------!
    ! [S] theta_dot - omega = 0
    !-----------------------------------------------------------------!

    !! May have to transform things into a frame
    R(2)  = body % S*body % thetadot - body % omega

    !-----------------------------------------------------------------!
    ! Dynamics eqn-1 (8 terms) (Force Equation)
    !-----------------------------------------------------------------!
    ! m (vdot - TIB*g) - c x omegadot + omega x (m v - c x omega) - fr = 0
    !-----------------------------------------------------------------!

    !! May have to transform things into a frame
    R(3)  =  body % mass*(body % vdot - body % grav) &
         & - skew(body % c)*body % omegadot &
         & + skew(body % omega)*(body % mass * body % v &
         & - skew(body % c)*body % omega) &
         & - body % rforce

    !-----------------------------------------------------------------!
    ! Dynamics eqn 2 (9-terms) (Moment Equation)
    !-----------------------------------------------------------------!
    ! c x vdot + J omegadot  + c x omega x v + omega x J - c x TIB*g - gr = 0
    !-----------------------------------------------------------------!

    !! May have to transform things into a frame
    R(4)  = skew(body % c) * body % vdot &
         & + body % J * body % omegadot &
         & + skew(body % c) * skew(body % omega)*body % v &
         & + skew(body % omega) * body % J * body % omega &
         & - skew(body % c) * body % grav &
         & - body % rtorque

  end function get_residual

  !-------------------------------------------------------------------!
  ! Jacobian of the kinematic and dynamic equations in state-space
  ! representation
  !-------------------------------------------------------------------!

  subroutine get_jacobian(this)

    class(rigid_body) :: this

  end subroutine get_jacobian

end module rigid_body_dynamics
