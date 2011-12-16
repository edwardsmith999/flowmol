!=============================================================================
!                                   
! Data structures needed by coupler internal subroutines
!
!-----------------------------------------------------------------------------
module coupler_internal_common
	implicit none
        
	integer COUPLER_COMM ! internal communicator, duplicate of realm_comm
	integer CFD_MD_ICOMM ! CFD - MD intracommunicator 

end module coupler_internal_common
