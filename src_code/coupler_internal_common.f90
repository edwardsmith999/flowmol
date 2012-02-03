!=============================================================================
!                                   
! Data structures needed by coupler internal subroutines
!
!-----------------------------------------------------------------------------
module coupler_internal_common
    implicit none    

    integer COUPLER_GLOBAL_COMM ! duplicate of MPI_COMM_WORLD, useful for input transfers
    integer COUPLER_REALM_COMM	! internal communicator inside the realm, split of COUPLER_GLOBAL_COMM
    integer COUPLER_ICOMM	! CFD - MD intracommunicator between COUPLER_REALM_COMM



end module coupler_internal_common
