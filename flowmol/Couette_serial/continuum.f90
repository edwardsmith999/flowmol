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
#if USE_COUPLER
	use continuum_coupler_socket, only : socket_coupler_invoke, & 
                                         socket_coupler_send_CFD_to_MD, &
                                         socket_coupler_init
#endif
	implicit none

#if USE_COUPLER
	call socket_coupler_invoke	
#endif
	call continuum_read_inputs				!Read from continuum input file parameters
	call continuum_set_parameters
	call continuum_mesh_gen

	call messenger_init
#if USE_COUPLER
	call socket_coupler_init			!Initialise coupler data structure and transfer coupled/CFD input to MD
#endif
	call continuum_setup_macrostate

#if USE_COUPLER
	call socket_coupler_send_CFD_to_MD
#endif

	call continuum_initial_record

#if USE_COUPLER
 	call continuum_set_BC
#endif

end subroutine setup_continuum

!=============================================================================
!
!                                SIMULATION
!
!-----------------------------------------------------------------------------

subroutine simulation_continuum
	use computational_constants
#if USE_COUPLER
	use continuum_coupler_socket, only : socket_coupler_send_CFD_to_MD
#endif
	implicit none

	integer :: t	      !Time t used to determine frequency of plots

	call continuum_calculate_flux
	!call continuum_advance_time_RK
	call continuum_advance_time
#if USE_COUPLER
	call socket_coupler_send_CFD_to_MD
#endif
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