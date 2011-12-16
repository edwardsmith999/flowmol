!=============================================================================
!
!                                  MAIN PROGRAM
!
! Viscous Laminar Flow solver
! Written by Edward Smith (unless otherwise commented), 11/01/10
! Updated 16/12/09
!
! program: continuum()
! subroutine: setup()
! subroutine: simulation()
! subroutine: finish()
!
!=============================================================================

!=============================================================================
!
!                                    SETUP  
!  
!-----------------------------------------------------------------------------
	
subroutine setup_continuum
    use continuum_coupler_socket, only : continuum_coupler_init, send_CFDvel
    implicit none

	call continuum_read_inputs
	call continuum_set_parameters
	call continuum_mesh_gen

    call messenger_init
    call continuum_coupler_init
    call continuum_setup_macrostate

    call send_CFDvel
	call continuum_initial_record
 	call continuum_set_BC

end subroutine setup_continuum

!=============================================================================
!
!                                SIMULATION
!
!-----------------------------------------------------------------------------

subroutine simulation_continuum
	use computational_constants
	use continuum_coupler_socket, only : send_CFDvel
	implicit none

	integer :: t	      !Time t used to determine frequency of plots

	call send_CFDvel	
	call continuum_calculate_flux
	!call continuum_advance_time_RK
	call continuum_advance_time
	call continuum_set_BC

end subroutine simulation_continuum

!=============================================================================
!
!                                FINISH
!
!-----------------------------------------------------------------------------

subroutine finish_continuum
	implicit none

	
end subroutine finish_continuum
