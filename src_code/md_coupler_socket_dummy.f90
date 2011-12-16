!=============================================================================
!				   
! Routine which interface with the coupler to the CFD code
! Corresponding empty shell dummy routine used for uncoupled calculation
!
!  Lucian Anton, November 2011
!
!-----------------------------------------------------------------------------
module md_coupler_socket
	implicit none

contains

!=============================================================================
! Setup initial times based on coupled calculation
!-----------------------------------------------------------------------------
subroutine socket_coupler_init
	implicit none

end subroutine socket_coupler_init

!=============================================================================
!Apply coupling forces so MD => CFD
!-----------------------------------------------------------------------------

subroutine socket_coupler_apply_continuum_forces(iter)
	implicit none
	
	integer, intent(in) :: iter
	
end subroutine socket_coupler_apply_continuum_forces

!=============================================================================
! Calculate averages of MD to pass to CFD
!-----------------------------------------------------------------------------

subroutine socket_coupler_average(iter)
	implicit none

	integer, intent(in) :: iter
	
end subroutine socket_coupler_average

end module md_coupler_socket
