! A single coupler module for both codes - this contains the same information 
! on both md and cfd side 
module coupler_module
	use coupler_parameters

	! <=><=><=><=> Simulation Time & Averaging data <=><=><=><=> 
	!Total number of steps for coupled simulation
	integer	:: nsteps_cfd, nsteps_md, nsteps_coupled
	! average period for averages ( it must come from CFD !!!)      
	integer :: average_period = 1
    ! save period ( corresponts to tplot in CFD, revise please !!!)
    integer :: save_period = 10
	!Timesteps
	real(kind(0.d0)) :: dt_md , dt_cfd

	! <=><=><=><=> Grid and domain data <=><=><=><=> 
	! Density
	real(kind(0.d0))	:: density_cfd, density_md
	! CFD number of cells
    integer	:: ngx,ngy,ngz
	integer	:: nlgx_cfd,nlgy_cfd,nlgz_cfd,nlgx_md,nlgy_md,nlgz_md
	! CFD grid indices
    integer	:: imin,imax,jmin,jmax,kmin,kmax 
	! CFD local grid indices (start and end of grid per CFD process)
    integer,dimension(:),allocatable :: iTmin,iTmax,jTmin,jTmax,kTmin,kTmax
	! Domain sizes
	real(kind(0.d0)) ::	xL_md,yL_md,zL_md,xL_cfd,yL_cfd,zL_cfd
	!CFD cells sizes 
	real(kind(0.d0)) 								   :: dx,dymin,dy,dymax,dz
    real(kind(0.d0)),dimension(:),  allocatable,target :: zpg
    real(kind(0.d0)),dimension(:,:),allocatable,target :: xpg,ypg

	! <=><=><=><=> Processor Topology & MPI <=><=><=><=> 
	!Processor id in grid
	integer	:: myid_grid
	! Number of processor in CFD grid
	integer :: npx_cfd, npy_cfd, npz_cfd, nproc_cfd	
    ! Number of processor in MD grid
    integer :: npx_md,  npy_md,  npz_md,  nproc_md
    ! Coordinates of MD/CFD topologies
    ! ATTENTION the values are shifted with +1, FORTRAN style
    ! if this array is passed by an MPI function, remove the shift
    integer,dimension(:,:),allocatable :: icoord_cfd, icoord_md
	! contains COMMS and other such things
    integer	:: COUPLER_GLOBAL_COMM ! duplicate of MPI_COMM_WORLD, useful for input transfers
    integer	:: COUPLER_REALM_COMM	! internal communicator inside the realm, split of COUPLER_GLOBAL_COMM
    integer	:: COUPLER_ICOMM		! CFD - MD INTER-communicator between COUPLER_REALM_COMM
    integer	:: COUPLER_GRID_COMM 	! duplicate of CFD or MD topology communicator
    integer	:: CFD_COMM_OVERLAP 	! communicator for tasks that overlap  MD region 


	! <=><=><=><=> COUPLER VARIABLES <=><=><=><=> 
 	integer	:: coupler_realm    	! NOT SURE YET???
    integer :: jmax_overlap = 5 ! maximum j index ( in y direction) which MD 
	real(kind(0.d0)) :: MD_initial_cellsize
    ! CFD code id
    integer :: cfd_code_id = couette_parallel

	! <=><=><=><=> OLD INTERNAL COMMON <=><=><=><=> 
    type overlap_map
        integer				 				:: n 		 ! number of ranks that overlap with this domain
        integer,dimension(:), allocatable 	:: rank_list ! rank list of overlapping bins
        integer,dimension(:,:),allocatable 	:: domains   ! range of overlapping indices between domain and overlapping boxes 
    end type overlap_map
    type(overlap_map) 	:: map

    ! flag marking 2d CFD solver
	logical 			:: cfd_is_2d = .false. ! set true if dz<=0 or kmax_cfd=kmin_cfd
    logical 			:: stop_request_activated = .false. ! request_abort is active or not (optimisation) 
    logical				:: staggered_averages(3) = (/ .false., .false., .false. /)
    integer, target		:: stop_request_tag
    integer, target		:: staggered_averages_tag
    character(len=64)	:: stop_request_name="none"

	interface error_abort
		module procedure error_abort_s, error_abort_si
	end interface error_abort

    private error_abort_si, error_abort_s

contains

!===========================================================================
!Error handling subroutines

subroutine error_abort_s(msg)
    use mpi
    implicit none

    character(len=*), intent(in), optional :: msg
   
	integer errcode,ierr

    if (present(msg)) then 
        write(*,*) msg
    endif

    call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)

end subroutine error_abort_s


subroutine error_abort_si(msg,i)
    use mpi
    implicit none

    character(len=*), intent(in) :: msg
    integer, intent(in) :: i

    integer errcode,ierr

    write(*,*) msg,i

    call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)

end subroutine error_abort_si

!===========================================================================
! Subroutine that can be used to stop the code when reaching a given 
! point in coupler -- useful when coupling new codes
!---------------------------------------------------------------------------
subroutine request_stop(tag)
    use mpi
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

end module coupler_module
