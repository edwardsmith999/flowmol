!=============================================================================
! COUPLER MODULE: 
! A single coupler module for both codes - this contains the same information 
! on both md and cfd side 
!
!	- Error handling
!	- MPI communicators
!	- Simulation realms
!	- MPI processor IDs
!	- Processor topologies
!	- Processor cartesian coords
!	- Global cell grid parameters
!	- Processor cell ranges
!	- Domain and cell dimensions
!	- Positions of CFD grid lines
!	- CFD to MD processor mapping
!	- Simulation parameters
!   - PROBABLY OBSOLETE STUFF
!
!=============================================================================

module coupler_module
	use coupler_parameters

	! MPI error flag
	integer		:: ierr	

	! MPI Communicators
	integer :: &

		CPL_WORLD_COMM,   & ! Copy of MPI_COMM_WORLD, both CFD and MD realms
		CPL_REALM_COMM,   & ! INTRA communicators within MD/CFD realms
		CPL_INTER_COMM,   & ! CFD/MD INTER communicator between realm comms
		CPL_CART_COMM,    & ! Comm w/cartesian topology for each realm
		CPL_OLAP_COMM,    & ! Local comm between only overlapping MD/CFD procs
		CPL_GRAPH_COMM,	  & ! Comm w/ graph topolgy between locally olapg procs
		CPL_REALM_INTERSECTION_COMM ! Intersecting MD/CFD procs in world

	! Simulation realms
	integer :: &
		realm

	! MPI processor IDs
	integer :: &
		myid_world,       &	!Processor ID from 0 to nproc_world-1
		rank_world,       & !Processor rank from 1 to nproc_world
		rootid_world,     &	!Root processor in world
		myid_realm,       & !Processor ID from 0 to nproc_realm-1
		rank_realm,       & !Processor rank from 1 to nproc_realm
		rootid_realm,     &	!Root processor in each realm
		myid_cart,        & !Processor ID from 0 to nproc_cart-1
		rank_cart,        & !Processor rank from 1 to nproc_cart
		rootid_cart,      &	!Root processor in each cart topology
		myid_olap,        & !Processor ID from 0 to nproc_olap-1
		rank_olap,        & !Processor rank from 1 to nproc_olap
		CFDid_olap,		  &	!Root processor in overlap is the CFD processor
		myid_graph,		  &	!Processor ID from 0 to nproc_graph-1
		rank_graph			!Processor rank from 1 to nproc_graph

	! Get rank in CPL_world_COMM from rank in local COMM
	integer, dimension(:), allocatable	:: &
		rank_world2rank_mdrealm,    &
		rank_world2rank_mdcart,     &
		rank_world2rank_cfdrealm,   &
		rank_world2rank_cfdcart,    &
		rank_world2rank_olap,       &
		rank_world2rank_graph,      &
		rank_world2rank_inter

	! Get rank in local COMM from rank in CPL_world_COMM
	integer, dimension(:), allocatable	:: &
		 rank_mdrealm2rank_world,    &
		  rank_mdcart2rank_world,    &
		rank_cfdrealm2rank_world,    &
		 rank_cfdcart2rank_world,    &
		    rank_olap2rank_world,    &
		   rank_graph2rank_world,    &
		   rank_inter2rank_world,    &
			rank_olap2rank_realm


	! Processor topologies
	integer :: &
		nproc_md,         &
		nproc_cfd,        &
		nproc_olap,		  &
		nproc_world,	  &
		npx_md,           &
		npy_md,           &
		npz_md,           &
		npx_cfd,          &
		npy_cfd,          &
		npz_cfd
	integer, dimension(:), allocatable :: &
		olap_mask
	integer, dimension(:,:), allocatable :: &
		rank2coord_cfd,   &
		rank2coord_md
	integer, dimension(:,:,:), allocatable :: &
		coord2rank_cfd,   &
		coord2rank_md

	! Processor cartesian coords	
	integer :: &
		iblock_realm,     &
		jblock_realm,     &
		kblock_realm

	! Global cell grid parameters
	integer :: &
		ncx,              &
		ncy,              &
		ncz,              &
		icmin,            &
		icmax,            &
		jcmin,            &
		jcmax,            &
		kcmin,            &
		kcmax,            &
		icmin_olap,       &
		icmax_olap,       &
		jcmin_olap,       &
		jcmax_olap,       &
		kcmin_olap,       &
		kcmax_olap,		  &
		ncx_olap,         &
		ncy_olap,         &
		ncz_olap
	
	! Processor cell ranges 
	integer, dimension(:), allocatable :: &
		icPmin_md,        &
		icPmax_md,        &
		jcPmin_md,        &
		jcPmax_md,        &
		kcPmin_md,        &
		kcPmax_md,        &
		icPmin_cfd,       &
		icPmax_cfd,       &
		jcPmin_cfd,       &
		jcPmax_cfd,       &
		kcPmin_cfd,       &
		kcPmax_cfd
	
	! Domain and cell dimensions
	real(kind(0.d0)) :: &
		xL_md,            &
		yL_md,            &
		zL_md,            &
		xL_cfd,           &
		yL_cfd,           &
		zL_cfd,           &
		xL_olap,          &
		yL_olap,          &
		zL_olap,          &
		xLl,              &
		yLl,              &
		zLl,              &
		dx,               &
		dy,               &
		dz,               &
		dymin,            &
		dymax

	! Positions of CFD grid lines
	real(kind(0.d0)), dimension(:,:), allocatable, target :: &
		xg,               &
		yg
	real(kind(0.d0)), dimension(:)  , allocatable, target :: &
		zg

	! CFD to MD processor mapping
	integer, dimension(:,:), allocatable :: &
		cfd_icoord2olap_md_icoords, &
		cfd_jcoord2olap_md_jcoords, &
		cfd_kcoord2olap_md_kcoords

	! Simulation parameters
	integer :: &
		nsteps_md,        & !MD input steps
		nsteps_cfd,       & !CFD input steps
		nsteps_coupled,   & !Total number of steps for coupled simulation
		average_period=1, & ! average period for averages ( it must come from CFD !!!)
		save_period=10      ! save period (corresponts to tplot in CFD, revise please !!!)
	real(kind(0.d0)) :: &
		dt_md,            &
		dt_cfd,           &
		density_md,       &
		density_cfd

	! PROBABLY OBSOLETE STUFF ------------------------------------------------!	
	real(kind(0.d0)) :: MD_initial_cellsize                                   !
    type overlap_map                                                          !
        integer                             :: n                              ! 
        integer,dimension(:), allocatable   :: rank_list                      !
        integer,dimension(:,:),allocatable  :: domains                        ! 
    end type overlap_map                                                      !
    type(overlap_map)   :: map                                                !
    integer :: cfd_code_id = couette_parallel     ! CFD code id  			  !
    ! flag marking 2d CFD solver
	logical 			:: cfd_is_2d = .false. 				! set true if dz<=0 or kmax_cfd=kmin_cfd
    logical 			:: stop_request_activated = .false. ! request_abort is active or not (optimisation) 
    logical				:: staggered_averages(3) = (/ .false., .false., .false. /)
    integer, target		:: stop_request_tag
    integer, target		:: staggered_averages_tag
    character(len=64)	:: stop_request_name="none"
	! CFD/MD number of cells
	!integer	:: nlgx_cfd,nlgy_cfd,nlgz_cfd, & 	!Local CFD
	!		   nlgx_md ,nlgy_md ,nlgz_md 		!Local MD
	! MD grid indices
	integer	:: ncy_puremd
	! Domain sizes
	real(kind(0.d0)) ::	yL_puremd
	
	integer :: testval
	!-------------------------------------------------------------------------!

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


subroutine messenger_lasterrorcheck
    use mpi
    implicit none

	integer resultlen
	character*12 err_buffer

	call MPI_Error_string(ierr,err_buffer,resultlen,ierr)
	print*, err_buffer

end subroutine messenger_lasterrorcheck


!--------------------------------------------------------------------------------------
! Prints formatted debug statements
subroutine printf(buf,dplaces_in)
	implicit none

	double precision,dimension(:),intent(in):: buf
	integer, intent(in), optional			:: dplaces_in

	integer				:: n,dplaces
	double precision	:: maxbuf,minbuf,order
	character*13	 	:: string
	character*42	 	:: buf_precision

	!Default number of decimal places if not supplied
	if (present(dplaces_in)) then
		if (dplaces_in .le. 9) then
			dplaces = dplaces_in
		else
			print*, 'Number of decimal places in printf if limited to 9'
			dplaces = 9 !Maximum
		endif
	else
		dplaces = 4
	endif

	!Find out required format to display maximum element in buffer
	maxbuf = maxval(buf); minbuf = minval(buf)
	maxbuf = max(maxbuf,10*abs(minbuf))	!10*Ensures extra space for minus sign
	order = 1.d0; n =1
	do while (max(maxbuf,order) .ne. order)
		order = order*10.d0
		n = n + 1
	enddo
	if (n+dplaces+2 .le. 9) then
		write(buf_precision,'(a,i1,a,i1)'), 'f',n+dplaces+2,'.', dplaces
	else
		write(buf_precision,'(a,i2,a,i1)'), 'f',n+dplaces+2,'.', dplaces
	endif

	! Build up format specifier string based on size of passed array
	string='(i3,   ' // trim(buf_precision) // ')'
	write(string(5:7),'(i3)'), size(buf) 

	!Write formatted data 
	print(string), rank_world,buf

end subroutine printf


!--------------------------------------------------------------------------------------
!Write matrix in correct format

subroutine write_matrix_int(a,varname,fh)
	implicit none

	integer					:: i,j,fh
	character(*)			:: varname
	integer, dimension(:,:) :: a

	write(fh,*) varname
	do i = lbound(a,1), ubound(a,1)
	    write(fh,*) (a(i,j), j = lbound(a,2), ubound(a,2))
	end do

end subroutine write_matrix_int

subroutine write_matrix(a,varname,fh)
	implicit none

	integer							 :: i,j,fh
	character(*)					 :: varname
	double precision, dimension(:,:) :: a

	write(fh,*) varname
	do i = lbound(a,1), ubound(a,1)
	    write(fh,*) (a(i,j), j = lbound(a,2), ubound(a,2))
	end do

end subroutine write_matrix

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
        call mpi_comm_rank(CPL_REALM_COMM, myid,ierr)
        write(0,*) 'stop as requested at ', trim(stop_request_name), ', realm',realm, 'rank', myid
        call MPI_Finalize(ierr)
        stop
    case("create_map","CREATE_MAP")
        call mpi_comm_rank(CPL_REALM_COMM, myid,ierr)
        write(0,*) 'stop as requested at ', trim(stop_request_name), ', realm',realm, 'rank', myid
        call MPI_Finalize(ierr)
        stop    
    case default
        write(0,*) "WARNING: request abort activated, but the tag is unrecognized, check COUPLER.in"
        write(0,*) "         accepted stop tags are: create_comm"
    end select

end subroutine request_stop


end module coupler_module
