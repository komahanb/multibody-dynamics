!-------------------------------------------------------------------!
! Function in Implicit form R(t, q, qot) = 0 
!-------------------------------------------------------------------!

subroutine R(res, nvars, time, q, qdot, qddot)

  implicit none

  integer, intent(in) :: nvars
  real(8), intent(in) :: time
  real(8), intent(in) :: q(nvars), qdot(nvars)
  real(8), OPTIONAL, intent(in) :: qddot(nvars)
  real(8), intent(inout) :: res(nvars)

  if (nvars .eq. 1) then

     ! res(1) = qdot(1) - sin(q(1))

     res(1) = qddot(1) + 0.02d0*qdot(1) + 5.0d0*q(1)

  else if (nvars .eq. 2) then 

     !res(1) = qddot(1) + 0.02d0*qdot(1)*qdot(2) + 5.0d0*q(1)
     !res(2) = qddot(2) - 0.05d0*qdot(2)*qdot(1) + 1.0d0*q(2)*q(1)
     
     ! Vanderpol equation       
     res(1) = qdot(1) - q(2)
     res(2) = qdot(2) - ( 1.0d0 - q(1)*q(1) )*q(2) + q(1)

  else

     res(1) =  qdot(1) - exp(q(2)) - sin(q(3))
     res(2) =  qdot(2) + 0.5*q(1) - 0.5*exp(q(3))
     res(3) =  qdot(3) + 0.5*exp(q(2)) - 2.5*q(3)

  end if

end subroutine R

!---------------------------------------------------------------------!
! DRDQ of the function
!---------------------------------------------------------------------!

subroutine DRDQ(J, alpha, nvars, time, q, qdot, qddot)

  implicit none

  integer, intent(in) :: nvars  
  real(8), intent(in) :: time, alpha
  real(8), intent(in) :: q(nvars), qdot(nvars)
  real(8), OPTIONAL, intent(in) :: qddot(nvars)
  real(8), intent(inout) :: J(nvars,nvars)  

  if (nvars .eq. 1) then

     J(1,1) = J(1,1) + alpha*5.0d0

  else if (nvars .eq. 2) then

!!$     ! derivative of first equation
!!$     
!!$     J(1,1) = J(1,1) + alpha*5.0d0
!!$     J(1,2) = J(1,2) + alpha*0.0d0
!!$     
!!$     ! derivative of second equation
!!$     
!!$     J(2,1) = J(2,1) + alpha*1.0d0*q(2)
!!$     J(2,2) = J(2,2) + alpha*1.0d0*q(1)

     ! Vanderpol equation

     ! derivative of first equation

     J(1,1) = J(1,1) + alpha*0.0d0
     J(1,2) = J(1,2) - alpha*1.0d0

     ! derivative of second equation

     J(2,1) = J(2,1) + alpha*(1.0d0 + 2.0d0*q(1)*q(2))
     J(2,2) = J(2,2) + alpha*(q(1)*q(1)-1.0d0)

!!$     end if

  else if (nvars .eq. 3) then

     ! res(1) =  qdot(1) - exp(q(2)) - sin(q(3))
     ! res(2) =  qdot(2) + 0.5*q(1) - 0.5*exp(q(3))
     ! res(3) =  qdot(3) + 0.5*exp(q(2)) - 2.5*q(3)

     ! derivative of the first equation

     J(1,1) = J(1,1) + alpha*0.0d0
     J(1,2) = J(1,2) - alpha*exp(q(2))
     J(1,3) = J(1,3) - alpha*cos(q(3))

     ! derivative of the second equation

     J(2,1) = J(2,1) + alpha*0.5d0
     J(2,2) = J(2,2) - alpha*0.0d0
     J(2,3) = J(2,3) - alpha*0.5*exp(q(3))

     ! derivative of the third equation

     J(3,1) = J(3,1) + alpha*0.0d0
     J(3,2) = J(3,2) + alpha*0.5*exp(q(2))
     J(3,3) = J(3,3) - alpha*2.5d0

     ! qdot(1) = exp(q(2)) + sin(q(3))
     ! qdot(2) = -0.5*q(1) + 0.5*exp(q(3))
     ! qdot(3) = -0.5*exp(q(2)) + 2.5*q(3)

  end if

end subroutine DRDQ

!---------------------------------------------------------------------!
! DRDQDOT of the function
!---------------------------------------------------------------------!

subroutine DRDQDOT(J, alpha, nvars, time, q, qdot, qddot)

  implicit none

  integer, intent(in) :: nvars  
  real(8), intent(in) :: time, alpha
  real(8), intent(in) :: q(nvars), qdot(nvars)
  real(8), OPTIONAL, intent(in) :: qddot(nvars)
  real(8), intent(inout) :: J(nvars,nvars)  

  if (nvars .eq. 1) then

     ! res(1) = qdot(1) - sin(q(1))

     J(1,1) = J(1,1) + alpha*0.02d0

  else if (nvars .eq. 2) then

!!$     ! derivative of first equation
!!$     
!!$     J(1,1) = J(1,1) + alpha*0.02d0*qdot(2)
!!$     J(1,2) = J(1,2) + alpha*0.02d0*qdot(1)
!!$
!!$     ! derivative of second equation
!!$
!!$     J(2,1) = J(2,1) - alpha*0.05d0*qdot(2)
!!$     J(2,2) = J(2,2) - alpha*0.05d0*qdot(1)

     ! Vanderpol equation

     ! derivative of first equation

     J(1,1) = J(1,1) + alpha*1.0d0
     J(1,2) = J(1,2) + alpha*0.0d0

     ! derivative of second equation

     J(2,1) = J(2,1) + alpha*0.0d0
     J(2,2) = J(2,2) + alpha*1.0d0

  else if (nvars .eq. 3) then

     ! derivative of first equation

     J(1,1) = J(1,1) + alpha*1.0d0
     J(1,2) = J(1,2) + alpha*0.0d0
     J(1,3) = J(1,3) + alpha*0.0d0

     ! derivative of second equation

     J(2,1) = J(2,1) + alpha*0.0d0
     J(2,2) = J(2,2) + alpha*1.0d0
     J(2,3) = J(2,3) + alpha*0.0d0

     ! derivative of third equation

     J(3,1) = J(3,1) + alpha*0.0d0
     J(3,2) = J(3,2) + alpha*0.0d0
     J(3,3) = J(3,3) + alpha*1.0d0

  end if

end subroutine DRDQDOT

!---------------------------------------------------------------------!
! DRDQDDOT of the function
!---------------------------------------------------------------------!

subroutine DRDQDDOT(J, alpha, nvars, time, q, qdot, qddot)

  implicit none

  integer, intent(in) :: nvars  
  real(8), intent(in) :: time, alpha
  real(8), intent(in) :: q(nvars), qdot(nvars)
  real(8), OPTIONAL, intent(in) :: qddot(nvars)
  real(8), intent(inout) :: J(nvars,nvars)  

  if (nvars .eq. 1) then

     ! res(1) = qdot(1) - sin(q(1))

     J(1,1) = J(1,1) + alpha*1.0d0

  else if (nvars .eq. 2) then

     !res(1) = qdot(1) - q(2)
     !res(2) = qdot(2) + 0.5d0*q(1) - 2.5d0*q(2)

     ! derivative of first equation

     J(1,1) = J(1,1) + alpha*0.0d0
     J(1,2) = J(1,2) + alpha*0.0d0

     ! derivative of second equation

     J(2,1) = J(2,1) + alpha*0.0d0
     J(2,2) = J(2,2) + alpha*0.0d0

  else if (nvars .eq. 3) then

     ! derivative of first equation

     J(1,1) = J(1,1) + alpha*1.0d0
     J(1,2) = J(1,2) + alpha*0.0d0
     J(1,3) = J(1,3) + alpha*0.0d0

     ! derivative of second equation

     J(2,1) = J(2,1) + alpha*0.0d0
     J(2,2) = J(2,2) + alpha*1.0d0
     J(2,3) = J(2,3) + alpha*0.0d0

     ! derivative of third equation

     J(3,1) = J(3,1) + alpha*0.0d0
     J(3,2) = J(3,2) + alpha*0.0d0
     J(3,3) = J(3,3) + alpha*1.0d0

  end if

end subroutine DRDQDDOT
