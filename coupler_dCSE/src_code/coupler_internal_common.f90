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
    integer COUPLER_GRID_COMM ! duplicate of CFD or MD topology communicator

    type overlap_map
        integer n ! number of ranks that overlap with this domain
        integer, allocatable :: rank_list(:) ! rank list of overelapping  bins
        integer, allocatable :: domains(:,:) ! range of overlapping indices between this domain and overlapping boxes 
    end type overlap_map
    type(overlap_map) :: map

    ! flag marking 2d CFD solver
    logical :: cfd_is_2d = .false. ! set true if dz<=0 or kmax_cfd=kmin_cfd

    logical :: stop_request_activated = .false. ! request_abort is active or not (optimisation) 
    
    character(len=64) :: stop_request_name="none"
    integer, target   :: stop_request_tag

    logical           :: staggered_averages(3) = (/ .false., .false., .false. /)
    integer, target   :: staggered_averages_tag

contains

!===========================================================================
! Subroutine that can be used to stop the code when reaching a given 
! point in coupler
! useful when coupling new codes
!---------------------------------------------------------------------------
subroutine request_stop(tag)
    use mpi
    use coupler_parameters
    implicit none

    character(len=*),intent(in) ::tag
    integer myid, ierr

    ! do nothing, get out quick 
    if(.not. stop_request_activated ) return

    if (tag /= stop_request_name) return

    select case(stop_request_name)
    case("create_comm","CREATE_COMM")
        call mpi_comm_rank(COUPLER_REALM_COMM, myid,ierr)
        write(0,*) 'stop as requested at ', trim(stop_request_name), ', realm',COUPLER_REALM, 'rank', myid
        call MPI_Finalize(ierr)
        stop
    case("create_map","CREATE_MAP")
        call mpi_comm_rank(COUPLER_REALM_COMM, myid,ierr)
        write(0,*) 'stop as requested at ', trim(stop_request_name), ', realm',COUPLER_REALM, 'rank', myid
        call MPI_Finalize(ierr)
        stop    
    case default
        write(0,*) "WARNING: request abort activated, but the tag is unrecognized, check COUPLER.in"
        write(0,*) "         accepted stop tags are: create_comm"
    end select

end subroutine request_stop


end module coupler_internal_common
