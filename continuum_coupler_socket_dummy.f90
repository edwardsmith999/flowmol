module continuum_coupler_socket
	implicit none

	logical, parameter :: use_coupling = .true.

contains

!=============================================================================
!  call coupler routine to build CFD communicator                                  
!-----------------------------------------------------------------------------
subroutine create_communicators(comm)
	implicit none
        
	integer, intent(out) :: comm

end subroutine create_communicators

!=============================================================================
! Call coupler routine to create map from CFD to MD                          
!-----------------------------------------------------------------------------
subroutine continuum_coupler_init 
	implicit none

end subroutine continuum_coupler_init

!=============================================================================
!  Package velocities and pass to coupler to send to MD                                
!-----------------------------------------------------------------------------
subroutine send_CFDvel
	implicit none

end subroutine send_CFDvel

!=============================================================================
!  Get MD velocities from coupler and unpack for CFD boundary conditions                                   
!-----------------------------------------------------------------------------
subroutine MD_continuum_BC(u,v)
	implicit none
	real(kind(0.d0)), intent(out) :: u(:), v(:)

 	u(:) =  0.d0
	v(:) =  0.d0

end subroutine MD_continuum_BC

end module continuum_coupler_socket

