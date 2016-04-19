!=====================================================================!
! Module that provides the ability to the user to implement random
! ODEs and use with the DIRK scheme
!=====================================================================!

module myode_class

  ! parent class
  use physics_class, only : physics

  implicit none

  private

  public :: ODE
  
  !-------------------------------------------------------------------!
  ! Type that models rigid body dynamics
  !-------------------------------------------------------------------!
  
  type, extends(physics) :: ODE

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

   contains

     procedure :: assembleResidual
     procedure :: assembleJacobian
     procedure :: getInitialStates

  end type ODE

contains

  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! -------------------------------------------------------------------!

subroutine assembleResidual(this, res, time, u, udot, uddot)

 class(ODE) :: this

 real(8), intent(inout), dimension(:) :: res
 real(8), intent(in) :: time
 real(8), intent(in), dimension(:) :: u, udot, uddot
 
 ! Temporary variables used to set the residuals (Some may not be used)
 real(8):: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, &
      t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, &
      t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, &
      t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, &
      t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, &
      t66, t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, &
      t79, t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, &
      t92, t93, t94, t95, t96, t97, t98, t99, t100, t101, t102, t103, &
      t104, t105, t106, t107, t108, t109, t110, t111, t112, t113, t114, &
      t115, t116, t117, t118, t119, t120, t121, t122, t123, t124, t125, &
      t126, t127, t128, t129, t130, t131, t132, t133, t134, t135, t136, &
      t137, t138, t139, t140, t141, t142, t143, t144, t145, t146, t147, &
      t148, t149, t150, t151, t152, t153, t154, t155, t156, t157, t158, &
      t159, t160, t161, t162, t163, t164, t165, t166, t167, t168, t169, &
      t170, t171, t172, t173, t174, t175, t176, t177, t178, t179, t180, &
      t181, t182, t183, t184, t185, t186, t187, t188, t189, t190, t191, &
      t192, t193, t194, t195, t196, t197, t198, t199, t200

 ! Variables for the project
 real(8):: mB, mL1, mL2, mL3, mL4, IByy, IL1yy, IL2yy, IL3yy
 real(8):: L1, L2, L3, L4
 real(8):: IL4yy, M1, M2, M3, M4, R3, R4, Radius, s
 real(8):: g

 ! Set all of these to constant values for now
 mB = 1000.0d0
 mL1 = 50.0d0
 mL2 = 50.0d0
 mL3 = 50.0d0
 mL4 = 50.0d0
 IByy = 1000.0d0
 IL1yy = 50.0d0
 IL2yy = 50.0d0
 IL3yy = 50.0d0
 IL4yy = 50.0d0
 L1 = 2.0d0
 L2 = 2.0d0
 L3 = 2.0d0
 L4 = 2.0d0
 M1 = 0.0d0
 M2 = 0.0d0
 M3 = 0.0d0
 M4 = 0.0d0
 R3 = 0.0d0
 R4 = 0.0d0
 Radius= 1.0d0
 s = 0.0d0
 g = 9.81d0


 t2 = u(3)
 t3 = sin(t2)
 t10 = cos(t2)
 t11 = uddot(1)
 t12 = t10 * t11
 t13 = uddot(2)
 t14 = t3 * t13
 t15 = uddot(3)
 t16 = t15 * s
 t17 = udot(3)
 t18 = t17 ** 2
 t19 = t18 * Radius
 t20 = uddot(4)
 t23 = u(4)
 t24 = cos(t23)
 t25 = (t15 - t20) * L1 * t24
 t27 = udot(4)
 t29 = (t17 - t27) ** 2
 t31 = sin(t23)
 t32 = t29 * L1 * t31
 t36 = uddot(7)
 t39 = u(7)
 t40 = cos(t39)
 t42 = sin(t39)
 t47 = udot(7)
 t49 = (t17 - t27 + t47) ** 2
 t63 = uddot(5)
 t66 = u(5)
 t67 = cos(t66)
 t68 = (t15 + t63) * L2 * t67
 t70 = udot(5)
 t72 = (t17 + t70) ** 2
 t74 = sin(t66)
 t75 = t72 * L2 * t74
 t79 = uddot(6)
 t82 = u(6)
 t83 = cos(t82)
 t85 = sin(t82)
 t90 = udot(6)
 t92 = (t17 + t70 - t90) ** 2
 res(1) = mB * g * t3 + mL1 * g * t3 + mL4 * g * t3 - R4 * t3 - mL1 * &
      (t12 - t14 + t16 - t19 - t25 / 0.2d0 - t32 / 0.2d0) - mL4 * (t12 &
      - t14 + t16 - t19 - t25 - t32 - (t15 - t20 + t36) * L4 * (t40 * t23 &
      + t42 * t31) / 0.2d0 + t49 * L4 * (t42 * t24 - t40 * t31) / 0.2d0 &
      ) + mL2 * g * t3 + mL3 * g * t3 - R3 * t3 - mL2 * (t12 - t14 + t16 &
      + t19 - t68 / 0.2d0 + t75 / 0.2d0) - mL3 * (t12 - t14 + t16 + t19 &
      - t68 + t75 - (t15 + t63 - t79) * L3 * (t83 * t24 - t85 * t31) / &
      0.2d0 + t92 * L3 * (-t85 * t24 - t83 * t31) / 0.2d0) - mB * (t12 &
      - t14)

 t2 = u(3)
 t3 = cos(t2)
 t10 = sin(t2)
 t11 = uddot(1)
 t12 = t10 * t11
 t13 = uddot(2)
 t14 = t3 * t13
 t15 = uddot(3)
 t16 = t15 * Radius
 t17 = udot(3)
 t18 = t17 ** 2
 t19 = t18 * s
 t20 = uddot(4)
 t23 = u(4)
 t24 = sin(t23)
 t25 = (t15 - t20) * L1 * t24
 t27 = udot(4)
 t29 = (t17 - t27) ** 2
 t31 = cos(t23)
 t32 = t29 * L1 * t31
 t36 = uddot(7)
 t39 = u(7)
 t40 = sin(t39)
 t42 = cos(t39)
 t47 = udot(7)
 t49 = (t17 - t27 + t47) ** 2
 t63 = uddot(5)
 t66 = u(5)
 t67 = sin(t66)
 t68 = (t15 + t63) * L2 * t67
 t70 = udot(5)
 t72 = (t17 + t70) ** 2
 t74 = cos(t66)
 t75 = t72 * L2 * t74
 t79 = uddot(6)
 t82 = u(6)
 t83 = sin(t82)
 t85 = cos(t82)
 t90 = udot(6)
 t92 = (t17 + t70 - t90) ** 2
 res(2) = -mB * g * t3 - mL1 * g * t3 - mL4 * g * t3 + R4 * t3 - mL1 &
      * (t12 + t14 - t16 - t19 - t25 / 0.2d0 + t32 / 0.2d0) - mL4 * (t12 &
      + t14 - t16 - t19 - t25 + t32 + (t15 - t20 + t36) * L4 * (-t42 * &
      t24 + t40 * t31) / 0.2d0 + t49 * L4 * (t40 * t24 + t42 * t31) / 0.2d0 &
      ) - mL2 * g * t3 - mL3 * g * t3 + R3 * t3 - mL2 * (t12 + t14 + &
      t16 - t19 + t68 / 0.2d0 + t75 / 0.2d0) - mL3 * (t12 + t14 + t16 - &
      t19 + t68 + t75 + (t15 + t63 - t79) * L3 * (-t85 * t24 - t83 * t31 &
      ) / 0.2d0 + t92 * L3 * (-t83 * t24 + t85 * t31) / 0.2d0) - mB * (t12 &
      + t14)

 t1 = mL1 * g
 t2 = u(3)
 t3 = sin(t2)
 t5 = mL4 * g
 t8 = cos(t2)
 t9 = uddot(1)
 t10 = t8 * t9
 t11 = uddot(2)
 t12 = t3 * t11
 t13 = uddot(3)
 t14 = t13 * s
 t15 = udot(3)
 t16 = t15 ** 2
 t17 = t16 * Radius
 t18 = uddot(4)
 t20 = (t13 - t18) * L1
 t21 = u(4)
 t22 = cos(t21)
 t23 = t20 * t22
 t25 = udot(4)
 t27 = (t15 - t25) ** 2
 t28 = t27 * L1
 t29 = sin(t21)
 t30 = t28 * t29
 t34 = uddot(7)
 t36 = (t13 - t18 + t34) * L4
 t37 = u(7)
 t38 = cos(t37)
 t40 = sin(t37)
 t42 = t38 * t22 + t40 * t29
 t45 = udot(7)
 t47 = (t15 - t25 + t45) ** 2
 t48 = t47 * L4
 t51 = t40 * t22 - t38 * t29
 t61 = t3 * t9
 t62 = t8 * t11
 t63 = t13 * Radius
 t64 = t16 * s
 t65 = t20 * t29
 t67 = t28 * t22
 t79 = mL2 * g
 t81 = mL3 * g
 t84 = uddot(5)
 t86 = (t13 + t84) * L2
 t87 = u(5)
 t88 = cos(t87)
 t89 = t86 * t88
 t91 = udot(5)
 t93 = (t15 + t91) ** 2
 t94 = t93 * L2
 t95 = sin(t87)
 t96 = t94 * t95
 t100 = uddot(6)
 t102 = (t13 + t84 - t100) * L3
 t103 = u(6)
 t104 = cos(t103)
 t106 = sin(t103)
 t108 = t104 * t22 - t106 * t29
 t111 = udot(6)
 t113 = (t15 + t91 - t111) ** 2
 t114 = t113 * L3
 t117 = -t104 * t29 - t106 * t22
 t127 = t86 * t95
 t129 = t94 * t88
 res(3) = (t1 * t3 + t5 * t3 - R4 * t3 - mL1 * (t10 - t12 + t14 - t17 &
      - t23 / 0.2d0 - t30 / 0.2d0) - mL4 * (t10 - t12 + t14 - t17 - t23 &
      - t30 - t36 * t42 / 0.2d0 + t48 * t51 / 0.2d0)) * s - (-t1 * t8 - &
      t5 * t8 + R4 * t8 - mL1 * (t61 + t62 - t63 - t64 - t65 / 0.2d0 + &
      t67 / 0.2d0) - mL4 * (t61 + t62 - t63 - t64 - t65 + t67 + t36 * t51 &
      / 0.2d0 + t48 * t42 / 0.2d0)) * Radius+ (t79 * t3 + t81 * t3 - R3 * &
      t3 - mL2 * (t10 - t12 + t14 + t17 - t89 / 0.2d0 + t96 / 0.2d0) - mL3 &
      * (t10 - t12 + t14 + t17 - t89 + t96 - t102 * t108 / 0.2d0 + t114 &
      * t117 / 0.2d0)) * s + (-t79 * t8 - t81 * t8 + R3 * t8 - mL2 * &
      (t61 + t62 + t63 - t64 + t127 / 0.2d0 + t129 / 0.2d0) - mL3 * (t61 &
      + t62 + t63 - t64 + t127 + t129 + t102 * t117 / 0.2d0 + t114 * t108 &
      / 0.2d0)) * Radius- M1 - M2 - IByy * t13

 t1 = u(4)
 t2 = sin(t1)
 t3 = L1 * t2
 t4 = mL1 * g
 t5 = u(3)
 t6 = cos(t5)
 t8 = mL4 * g
 t9 = t8 * t6
 t10 = R4 * t6
 t11 = sin(t5)
 t12 = uddot(1)
 t13 = t11 * t12
 t14 = uddot(2)
 t15 = t6 * t14
 t16 = uddot(3)
 t17 = t16 * Radius
 t18 = udot(3)
 t19 = t18 ** 2
 t20 = t19 * s
 t21 = uddot(4)
 t22 = t16 - t21
 t23 = t22 * L1
 t24 = t23 * t2
 t26 = udot(4)
 t28 = (t18 - t26) ** 2
 t29 = t28 * L1
 t30 = cos(t1)
 t31 = t29 * t30
 t35 = uddot(7)
 t37 = (t16 - t21 + t35) * L4
 t38 = u(7)
 t39 = sin(t38)
 t41 = cos(t38)
 t43 = -t41 * t2 + t39 * t30
 t46 = udot(7)
 t48 = (t18 - t26 + t46) ** 2
 t49 = t48 * L4
 t52 = t39 * t2 + t41 * t30
 t56 = mL4 * (t13 + t15 - t17 - t20 - t24 + t31 + t37 * t43 / 0.2d0 &
      + t49 * t52 / 0.2d0)
 t60 = L1 * t30
 t62 = t8 * t11
 t63 = R4 * t11
 t64 = t6 * t12
 t65 = t11 * t14
 t66 = t16 * s
 t67 = t19 * Radius
 t68 = t23 * t30
 t70 = t29 * t2
 t79 = mL4 * (t64 - t65 + t66 - t67 - t68 - t70 - t37 * t52 / 0.2d0 &
      + t49 * t43 / 0.2d0)
 res(4) = M1 - M4 - t3 * (-t4 * t6 - t9 + t10 - mL1 * (t13 + t15 - t17 &
      - t20 - t24 / 0.2d0 + t31 / 0.2d0) - t56) / 0.2d0 - t60 * (t4 * t11 &
      + t62 - t63 - mL1 * (t64 - t65 + t66 - t67 - t68 / 0.2d0 - t70 &
      / 0.2d0) - t79) / 0.2d0 - t3 * (-t9 + t10 - t56) / 0.2d0 - t60 * ( &
      t62 - t63 - t79) / 0.2d0 - IL1yy * t22

 t1 = u(5)
 t2 = sin(t1)
 t3 = L2 * t2
 t4 = mL2 * g
 t5 = u(3)
 t6 = cos(t5)
 t8 = mL3 * g
 t9 = t8 * t6
 t10 = R3 * t6
 t11 = sin(t5)
 t12 = uddot(1)
 t13 = t11 * t12
 t14 = uddot(2)
 t15 = t6 * t14
 t16 = uddot(3)
 t17 = t16 * Radius
 t18 = udot(3)
 t19 = t18 ** 2
 t20 = t19 * s
 t21 = uddot(5)
 t22 = t16 + t21
 t23 = t22 * L2
 t24 = t23 * t2
 t26 = udot(5)
 t28 = (t18 + t26) ** 2
 t29 = t28 * L2
 t30 = cos(t1)
 t31 = t29 * t30
 t35 = uddot(6)
 t37 = (t16 + t21 - t35) * L3
 t38 = u(6)
 t39 = sin(t38)
 t40 = u(4)
 t41 = cos(t40)
 t43 = cos(t38)
 t44 = sin(t40)
 t46 = -t39 * t41 - t43 * t44
 t49 = udot(6)
 t51 = (t18 + t26 - t49) ** 2
 t52 = t51 * L3
 t55 = -t39 * t44 + t43 * t41
 t59 = mL3 * (t13 + t15 + t17 - t20 + t24 + t31 + t37 * t46 / 0.2d0 &
      + t52 * t55 / 0.2d0)
 t63 = L2 * t30
 t65 = t8 * t11
 t66 = R3 * t11
 t67 = t6 * t12
 t68 = t11 * t14
 t69 = t16 * s
 t70 = t19 * Radius
 t71 = t23 * t30
 t73 = t29 * t2
 t82 = mL3 * (t67 - t68 + t69 + t70 - t71 + t73 - t37 * t55 / 0.2d0 &
      + t52 * t46 / 0.2d0)
 res(5) = M2 - M3 + (-t4 * t6 - t9 + t10 - mL2 * (t13 + t15 + t17 - t20 &
      + t24 / 0.2d0 + t31 / 0.2d0) - t59) * t3 / 0.2d0 - t63 * (t4 * t11 &
      + t65 - t66 - mL2 * (t67 - t68 + t69 + t70 - t71 / 0.2d0 + t73 &
      / 0.2d0) - t82) / 0.2d0 + t3 * (-t9 + t10 - t59) / 0.2d0 - t63 * ( &
      t65 - t66 - t82) / 0.2d0 - IL2yy * t22

 t1 = u(6)
 t2 = sin(t1)
 t3 = u(4)
 t4 = cos(t3)
 t6 = cos(t1)
 t7 = sin(t3)
 t9 = -t2 * t4 - t6 * t7
 t10 = L3 * t9
 t11 = mL3 * g
 t12 = u(3)
 t13 = cos(t12)
 t16 = sin(t12)
 t17 = uddot(1)
 t19 = uddot(2)
 t21 = uddot(3)
 t23 = udot(3)
 t24 = t23 ** 2
 t26 = uddot(5)
 t28 = (t21 + t26) * L2
 t29 = u(5)
 t30 = sin(t29)
 t32 = udot(5)
 t34 = (t23 + t32) ** 2
 t35 = t34 * L2
 t36 = cos(t29)
 t38 = uddot(6)
 t39 = t21 + t26 - t38
 t40 = t39 * L3
 t43 = udot(6)
 t45 = (t23 + t32 - t43) ** 2
 t46 = t45 * L3
 t49 = -t2 * t7 + t4 * t6
 t57 = L3 * t49
 res(6) = M3 - t10 * (-t11 * t13 + R3 * t13 - mL3 * (t16 * t17 + t13 * &
      t19 + t21 * Radius- t24 * s + t28 * t30 + t35 * t36 + t40 * t9 / 0.2d0 &
      + t46 * t49 / 0.2d0)) / 0.2d0 + t57 * (t11 * t16 - R3 * t16 - mL3 &
      * (t13 * t17 - t16 * t19 + t21 * s + t24 * Radius- t28 * t36 + t35 * &
      t30 - t40 * t49 / 0.2d0 + t46 * t9 / 0.2d0)) / 0.2d0 + R3 * (t10 &
      * t13 + t57 * t16) / 0.2d0 - IL3yy * t39

 t1 = u(7)
 t2 = sin(t1)
 t3 = u(4)
 t4 = cos(t3)
 t6 = cos(t1)
 t7 = sin(t3)
 t9 = t2 * t4 - t6 * t7
 t10 = L4 * t9
 t11 = mL4 * g
 t12 = u(3)
 t13 = cos(t12)
 t16 = sin(t12)
 t17 = uddot(1)
 t19 = uddot(2)
 t21 = uddot(3)
 t23 = udot(3)
 t24 = t23 ** 2
 t26 = uddot(4)
 t28 = (t21 - t26) * L1
 t30 = udot(4)
 t32 = (t23 - t30) ** 2
 t33 = t32 * L1
 t35 = uddot(7)
 t36 = t21 - t26 + t35
 t37 = t36 * L4
 t40 = udot(7)
 t42 = (t23 - t30 + t40) ** 2
 t43 = t42 * L4
 t46 = t2 * t7 + t4 * t6
 t54 = L4 * t46

 res(7) = M4 - t10 * (-t11 * t13 + R4 * t13 - mL4 * (t16 * t17 + t13 * &
      t19 - t21 * Radius- t24 * s - t28 * t7 + t33 * t4 + t37 * t9 / 0.2d0 &
      + t43 * t46 / 0.2d0)) / 0.2d0 + t54 * (t11 * t16 - R4 * t16 - mL4 &
      * (t13 * t17 - t16 * t19 + t21 * s - t24 * Radius - t28 * t4 - t33 * t7 &
      - t37 * t46 / 0.2d0 + t43 * t9 / 0.2d0)) / 0.2d0 + R4 * (t10 * t13 &
      + t54 * t16) / 0.2d0 - IL4yy * t36

end subroutine assembleResidual

!-------------------------------------------------------------------!
! Jacobian assembly at each time step. If you don't provide the
! analytical jacobian, set setApproximateJacobian(.true.) into the
! integrator object. We use finite-difference method to approximate
! the Jacobian.
! 
! Jacobian is the matrix of partial derivatives. Each row in the
! Jacobian matrix arises from differntiating a single equation. Each
! column in the Jacobian comes from a variable in the problem. Note
! the NEQN should be equal to NVARS for the system to be solved.
!
! Note: alpha, beta and gamma are scalars that need to be multiplied
! with the partial derivatives DRDQ, DRDQDOT and DRDQDDOT
! respectively.
! -------------------------------------------------------------------!

subroutine assembleJacobian(this, jac, alpha, beta, gamma, &
    & time, u, udot, uddot)

 class(ODE) :: this
 real(8), intent(inout), dimension(:,:) :: jac
 real(8), intent(in) :: alpha, beta, gamma
 real(8), intent(in) :: time
 real(8), intent(in), dimension(:) :: u, udot, uddot

end subroutine assembleJacobian

!---------------------------------------------------------------------!
! Sets the initial condition for use in the integator. If first order
! system just set initial Q, if a second order system set initial Q
! and QDOT
!---------------------------------------------------------------------!

subroutine getInitialStates(this, time, u, udot)

 class(ODE) :: this

 real(8), intent(in) :: time
 real(8), intent(inout), dimension(:) :: u, udot

 ! Set the initial conditions for all 7 u's
 u(1)    = 1.0d0
 u(2)    = 1.0d0
 u(3)    = 1.0d0
 u(4)    = 1.0d0
 u(5)    = 1.0d0
 u(6)    = 1.0d0
 u(7)    = 1.0d0

 ! Set the initial conditions for all 7 udot's
 udot(1) = 1.0d0
 udot(2) = 2.0d0
 udot(3) = 1.0d0
 udot(4) = 2.0d0
 udot(5) = 1.0d0
 udot(6) = 2.0d0
 udot(7) = 1.0d0

end subroutine getInitialStates

end module myode_class








