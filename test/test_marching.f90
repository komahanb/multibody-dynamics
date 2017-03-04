!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

#include "scalar.fpp"

program test_time_integration
  
  ! Import all test physics
  use spring_mass_damper_class      , only : smd1, smd2
  use vanderpol_class               , only : vanderpol
  use aero_elastic_oscillator_class , only : aero_elastic_oscillator

  implicit none

  ! Declare Physics for testing
  type(smd1)      , target :: smd1obj                 ! Spring-mass-damper test ODE (1 var)
  type(smd2)      , target :: smd2obj                 ! Spring-mass-damper test ODE (2 var)
  type(vanderpol) , target :: vpl                     ! Vanderpol equation (2 var)
  type(aero_elastic_oscillator), target :: aeosc      ! Aeroelastic oscillator (2 vars)

  ! Test the integrators in vanderpol oscillator system  
  test_vanderpol: block   
    call vpl % initialize("Vanderpol", num_state_vars = 2, num_design_vars = 1)
    call test_integrators(vpl, "vanderpol", .false.)
    call vpl % finalize()    
  end block test_vanderpol
  
  ! Test the integrators on 2 dof spring mass damper system
  test_smd2: block
    call smd2obj % initialize("SMD2", num_state_vars = 2)
    call test_integrators(smd2obj, "smd2", .true.)
    call smd2obj % finalize()
  end block test_smd2

  ! Test the integrators on aeroelastic oscillator problem with
  ! pitching and plunging degree of freedom
  test_aeosc: block
    call aeosc % initialize("AeroElasticOscillator", num_state_vars = 2)
    call test_integrators(aeosc, "aeosc", .true.)
    call aeosc % finalize()
  end block test_aeosc

contains

  subroutine test_integrators( test_system, name, second_order)

    ! Import all time integrators
    use runge_kutta_integrator        , only : DIRK
    use bdf_integrator                , only : BDF
    use abm_integrator                , only : ABM
    use nbg_integrator                , only : NBG

    ! Import physics base class
    use physics_class                 , only : physics

    class(physics)  , intent(in) :: test_system    
    character(len=*), intent(in) :: name
    logical                      :: second_order
    ! Declare integrators
    type(DIRK) :: dirkobj   ! DIRK Integrator object
    type(BDF)  :: bdfobj    ! BDF Integrator object
    type(ABM)  :: abmobj    ! ABM Integrator object
    type(NBG)  :: nbgobj    ! NBM Integrator object
    
    !=================================================================!
    !                        TEST BDF                                 !
    !=================================================================!

    bdfobj = BDF(system = test_system, tfinal = 20.0d0, h=1.0d-2, &
         & max_bdf_order = 3, second_order=second_order)
    call bdfobj % setPrintLevel(2)
    call bdfobj % integrate()
    call bdfobj % writeSolution(name//"_bdf.dat")
    call bdfobj % finalize()

    !=================================================================!
    !                     TEST NBG                                    !
    !=================================================================!

    nbgobj = NBG(system = test_system, tfinal = 20.0d0, h=1.0d-2, &
         & second_order=second_order)
    call nbgobj % setPrintLevel(2)
    call nbgobj % integrate()
    call nbgobj % writeSolution(name//"_nbg.dat")
    call nbgobj % finalize()

    !=================================================================!
    !                     TEST ABM                                    !
    !==================================================================!

    abmobj = ABM(system = test_system, tfinal = 20.0d0, h=1.0d-2, &
         & max_abm_order = 3, second_order=second_order)
    call abmobj % setPrintLevel(2)
    call abmobj % integrate()
    call abmobj % writeSolution(name//"_abm.dat")
    call abmobj % finalize()

    !=================================================================!
    !                        TEST DIRK                                !
    !=================================================================!

    dirkobj = DIRK(system = test_system, tfinal = 20.0d0, h=1.0d-2, &
         & num_stages=2, second_order=second_order) 
    call dirkobj % setPrintLevel(2)
    call dirkobj % integrate()
    call dirkobj % writeSolution(name//"_dirk.dat")
    call dirkobj % finalize()

  end subroutine test_integrators
  
end program test_time_integration

