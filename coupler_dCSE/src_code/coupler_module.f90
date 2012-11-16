!=============================================================================
!! COUPLER MODULE: 
!! A single coupler module for both codes - this contains the same information 
!! on both md and cfd side 
!!
!!  - Error handling
!!  - MPI communicators
!!  - Simulation realms
!!  - MPI processor IDs
!!  - Processor topologies
!!  - Processor cartesian coords
!!  - Global cell grid parameters
!!  - Processor cell ranges
!!  - Domain and cell dimensions
!!  - Positions of CFD grid lines
!!  - CFD to MD processor mapping
!!  - Simulation parameters
!!
!! The data is protected so only setup routines in this module can change it
! SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP
!! Setup routine which have access to coupler parameters
!!
!! - CPL_create_comm	      (cfd+md)   splits MPI_COMM_WORLD, create inter - 
!!                                   communicator between CFD and MD
!!
!! - CPL_create_map	      	  (cfd+md)   creates correspondence maps between 
!!                                      the CFD grid and MD domains
!!
!! - CPL_cfd_adjust_domain     (cfd)    adjust CFD tomain to an integer number 
!!                                      FCC or similar MD initial layout
!! @author David Trevelyan, Edward Smith
!! @see coupler
!=============================================================================

module coupler_module
	USE ISO_C_BINDING
	implicit none

 	integer, parameter :: VOID=-666			!!VOID value for data initialisation
	integer, parameter :: cfd_realm = 1		!! CFD realm identifier
	integer, parameter :: md_realm  = 2		!! MD realm identifier
	character(len=*),parameter :: &
		realm_name(2) = (/ "CFD", "MD " /) 	!! Used with realm identifier to get name

	!! error codes
	integer, parameter :: & 
		COUPLER_ERROR_REALM  = 1  		!! wrong realm value
	integer, parameter :: & 
		COUPLER_ERROR_ONE_REALM  = 2	!! one realm missing
	integer, parameter :: & 
		COUPLER_ERROR_INIT       = 3	!! initialisation error
	integer, parameter :: & 
		COUPLER_ERROR_INPUT_FILE = 4	!! wrong value in input file
	integer, parameter :: & 
		COUPLER_ERROR_READ_INPUT = 5	!! error in processing input file or data transfers
	integer, parameter :: & 
		COUPLER_ERROR_CONTINUUM_FORCE = 6	!!the region in which the continuum constrain force is apply spans over two MD domains
	integer, parameter :: & 
		COUPLER_ABORT_ON_REQUEST = 7	!! used in request_abort 
	integer, parameter :: & 
		COUPLER_ABORT_SEND_CFD   = 8	!! error in coupler_cfd_send
	integer, parameter :: & 
		COUPLER_ERROR_CART_COMM   = 9 	!! Wrong comm value in CPL_Cart_coords

	!! MPI error flag
	integer		:: ierr	

	! MPI Communicators
	integer,protected :: &
		CPL_WORLD_COMM !! Copy of MPI_COMM_WORLD, both CFD and MD realms;
	integer,protected :: &
		CPL_REALM_COMM !! INTRA communicators within MD/CFD realms;
	integer,protected :: &
		CPL_INTER_COMM !!  CFD/MD INTER communicator between realm comms;
	integer,protected :: &
		CPL_CART_COMM !!  Comm w/cartesian topology for each realm;
	integer,protected :: &
		CPL_OLAP_COMM !!  Local comm between only overlapping MD/CFD procs;
	integer,protected :: &
		CPL_GRAPH_COMM !!  Comm w/ graph topolgy between locally olapg procs;
	integer,protected :: &
		CPL_REALM_INTERSECTION_COMM !!  Intersecting MD/CFD procs in world;

	!! Simulation realms
	integer,protected :: &
		realm

	! MPI processor IDs
	integer,protected :: &
		myid_world   !! Processor ID from 0 to nproc_world-1;
	integer,protected :: &
		rank_world   !! Processor rank from 1 to nproc_world;
	integer,protected :: &
		rootid_world !! Root processor in world;
	integer,protected :: &
		myid_realm   !! Processor ID from 0 to nproc_realm-1;
	integer,protected :: &
		rank_realm   !! Processor rank from 1 to nproc_realm;
	integer,protected :: &
		rootid_realm !! Root processor in each realm;
	integer,protected :: &
		myid_cart   !! Processor ID from 0 to nproc_cart-1;
	integer,protected :: &
		rank_cart   !! Processor rank from 1 to nproc_cart;
	integer,protected :: &
		rootid_cart !! Root processor in each cart topology;
	integer,protected :: &
		myid_olap   !! Processor ID from 0 to nproc_olap-1;
	integer,protected :: &
		rank_olap   !! Processor rank from 1 to nproc_olap;
	integer,protected :: &
		CFDid_olap  !! Root processor in overlap is the CFD processor;
	integer,protected :: &
		myid_graph  !! Processor ID from 0 to nproc_graph-1;
	integer,protected :: &
		rank_graph  !! Processor rank from 1 to nproc_graph;

	!! Get rank in CPL_world_COMM from rank in local COMM
	integer,protected, dimension(:), allocatable	:: &
		rank_world2rank_mdrealm,    &
		rank_world2rank_mdcart,     &
		rank_world2rank_cfdrealm,   &
		rank_world2rank_cfdcart,    &
		rank_world2rank_olap,       &
		rank_world2rank_graph,      &
		rank_world2rank_inter

	!! Get rank in local COMM from rank in CPL_world_COMM
	integer,protected, dimension(:), allocatable	:: &
		 rank_mdrealm2rank_world,    &
		  rank_mdcart2rank_world,    &
		rank_cfdrealm2rank_world,    &
		 rank_cfdcart2rank_world,    &
		    rank_olap2rank_world,    &
		   rank_graph2rank_world,    &
		   rank_inter2rank_world,    &
			rank_olap2rank_realm


	! Processor topologies
	integer,protected :: &
		nproc_md     !! Total number of processor in md
	integer,protected :: &
		nproc_cfd    !! Total number of processor in cfd
	integer,protected :: &
		nproc_olap   !! Total number of processor in overlap region
	integer,protected :: &
		nproc_world  !! Total number of processor in world
	integer,protected :: &
		npx_md       !! Number of processor in x in the md
	integer,protected :: &
		npy_md       !! Number of processor in y in the md
	integer,protected :: &
		npz_md       !! Number of processor in z in the md
	integer,protected :: &
		npx_cfd      !! Number of processor in x in the cfd
	integer,protected :: &
		npy_cfd      !! Number of processor in y in the cfd
	integer,protected :: &
		npz_cfd      !! Number of processor in z in the cfd

	logical,protected, dimension(:), allocatable :: &
		olap_mask           !! Overlap mask specifying which processors overlap using world ranks
	integer,protected, dimension(:,:), allocatable :: &
		rank2coord_cfd,   &   !! Array containing coordinates for each cartesian rank 
		rank2coord_md
	integer,protected, dimension(:,:,:), allocatable :: &
		coord2rank_cfd,   &
		coord2rank_md

	!! Processor cartesian coords	
	integer,protected :: &
		iblock_realm,     &
		jblock_realm,     &
		kblock_realm

	! Global cell grid parameters
	integer,protected :: &
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
	integer,protected, dimension(:), allocatable :: &
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
	real(kind(0.d0)),protected :: &
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
	real(kind(0.d0)),protected, dimension(:,:), allocatable, target :: &
		xg,               &
		yg
	real(kind(0.d0)),protected, dimension(:)  , allocatable, target :: &
		zg

	! CFD to MD processor mapping
	integer,protected, dimension(:,:), allocatable :: &
		cfd_icoord2olap_md_icoords, &
		cfd_jcoord2olap_md_jcoords, &
		cfd_kcoord2olap_md_kcoords

	! Simulation parameters
	integer,protected :: &
		nsteps_md,        & !MD input steps
		nsteps_cfd,       & !CFD input steps
		nsteps_coupled,   & !Total number of steps for coupled simulation
		average_period=1, & ! average period for averages ( it must come from CFD !!!)
		save_period=10      ! save period (corresponts to tplot in CFD, revise please !!!)
	real(kind(0.d0)),protected :: &
		dt_md,            &
		dt_cfd,           &
		density_md,       &
		density_cfd
	integer,protected :: &
		timestep_ratio,        &
		md_cfd_match_cellsize, &
		testval
	logical,protected :: &
		staggered_averages(3) = (/.false.,.false.,.false./)

	!Interface for error handling functions
	interface error_abort
		module procedure error_abort_s, error_abort_si
	end interface error_abort

    private error_abort_si, error_abort_s

contains


!=============================================================================
!					 _____      _               
!					/  ___|    | |              
!					\ `--.  ___| |_ _   _ _ __  
!					 `--. \/ _ \ __| | | | '_ \ 
!					/\__/ /  __/ |_| |_| | |_) |
!					\____/ \___|\__|\__,_| .__/ 
!					                     | |    
!					                     |_|    
!=============================================================================

!=============================================================================
! 							coupler_create_comm	      	
!! (cfd+md) Splits MPI_COMM_WORLD in both the CFD and MD code respectively
!! 		   and create intercommunicator between CFD and MD
!-----------------------------------------------------------------------------

subroutine CPL_create_comm(callingrealm, RETURNED_REALM_COMM, ierror)
	use mpi
	!use coupler_module, only : rank_world,myid_world,rootid_world,nproc_world,&
	!                           realm, rank_realm,myid_realm,rootid_realm,ierr,& 
	!                           CPL_WORLD_COMM, CPL_REALM_COMM, CPL_INTER_COMM
	implicit none

	integer, intent(in) :: callingrealm ! CFD or MD
	integer, intent(out):: RETURNED_REALM_COMM, ierror

	!Get processor id in world across both realms
	call MPI_comm_rank(MPI_COMM_WORLD,myid_world,ierr)
	rank_world = myid_world + 1; rootid_world = 0
	call MPI_comm_size(MPI_COMM_WORLD,nproc_world,ierr)

	! test if we have a CFD and a MD realm
	ierror=0
	! Test realms are assigned correctly
	call test_realms	

	! Create intercommunicator between realms		
	realm = callingrealm
	call create_comm		

contains

!-----------------------------------------------------------------------------
!	Test if CFD and MD realms are assigned correctly
!-----------------------------------------------------------------------------

subroutine test_realms
	implicit none

	integer 			 :: i, root, nproc, ncfd, nmd
	integer, allocatable :: realm_list(:)

	! Allocate and gather array with realm (MD or CFD) of each 
	! processor in the coupled domain on the root processor
	root = 1
	if (rank_world .eq. root) then
		call MPI_comm_size(MPI_comm_world, nproc, ierr)
		allocate(realm_list(nproc))
	endif
	call MPI_gather(callingrealm,1,MPI_INTEGER,realm_list,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	!Check through array of processors on both realms
	!and return error if wrong values or either is missing
	if (rank_world .eq. root) then
		ncfd = 0; nmd = 0
		do i =1, nproc
			if ( realm_list(i) .eq. cfd_realm ) then 
				ncfd = ncfd + 1
			else if ( realm_list(i) .eq. md_realm ) then
				nmd = nmd +1
			else
				ierror = COUPLER_ERROR_REALM
				write(*,*) "wrong realm value in coupler_create_comm"
				call MPI_abort(MPI_COMM_WORLD,ierror,ierr)
			endif
		enddo

		if ( ncfd .eq. 0 .or. nmd .eq. 0) then 
			ierror = COUPLER_ERROR_ONE_REALM
			write(*,*) "CFD or MD realm is missing in MPI_COMM_WORLD"
			call MPI_abort(MPI_COMM_WORLD,ierror,ierr)
		endif

	endif

end subroutine test_realms

!-----------------------------------------------------------------------------
! Create communicators for each realm and inter-communicator
!-----------------------------------------------------------------------------

subroutine create_comm
	implicit none

	integer	::  callingrealm,ibuf(2),jbuf(2),remote_leader,comm,comm_size

	callingrealm = realm
	! Split MPI COMM WORLD ready to establish two communicators
	! 1) A global intra-communicator in each realm for communication
	! internally between CFD processes or between MD processes
	! 2) An inter-communicator which allows communication between  
	! the 'groups' of processors in MD and the group in the CFD 
	call MPI_comm_dup(MPI_COMM_WORLD,CPL_WORLD_COMM,ierr)
	RETURNED_REALM_COMM	= MPI_COMM_NULL
	CPL_REALM_COMM 		= MPI_COMM_NULL

	!------------ create realm intra-communicators -----------------------
	! Split MPI_COMM_WORLD into an intra-communicator for each realm 
	! (used for any communication within each realm - e.g. broadcast from 
	!  an md process to all other md processes) 
	call MPI_comm_split(CPL_WORLD_COMM,callingrealm,myid_world,RETURNED_REALM_COMM,ierr)

	!------------ create realm inter-communicators -----------------------
	! Create intercommunicator between the group of processor on each realm
	! (used for any communication between realms - e.g. md group rank 2 sends
	! to cfd group rank 5). inter-communication is by a single processor on each group
	! Split duplicate of MPI_COMM_WORLD
	call MPI_comm_split(CPL_WORLD_COMM,callingrealm,myid_world,CPL_REALM_COMM,ierr)
	call MPI_comm_rank(CPL_REALM_COMM,myid_realm,ierr)
	rank_realm = myid_realm + 1; rootid_realm = 0

	! Get the MPI_comm_world ranks that hold the largest ranks in cfd_comm and md_comm
	call MPI_comm_size(CPL_REALM_COMM,comm_size,ierr)
	ibuf(:) = -1
	jbuf(:) = -1
	if ( myid_realm .eq. comm_size - 1) then
		ibuf(realm) = myid_world
	endif

	call MPI_allreduce( ibuf ,jbuf, 2, MPI_INTEGER, MPI_MAX, &
						CPL_WORLD_COMM, ierr)

	!Set this largest rank on each process to be the inter-communicators (WHY NOT 0??)
	select case (realm)
	case (cfd_realm)
		remote_leader = jbuf(md_realm)
	case (md_realm)
		remote_leader = jbuf(cfd_realm)
	end select

	!print*,color, jbuf, remote_leader

	call MPI_intercomm_create(CPL_REALM_COMM, comm_size - 1, CPL_WORLD_COMM,&
									remote_leader, 1, CPL_INTER_COMM, ierr)
	print*, 'did (inter)communicators ', realm_name(realm), myid_world

end subroutine create_comm

end subroutine CPL_create_comm

!=============================================================================
!! Read Coupler input file
!-----------------------------------------------------------------------------

subroutine read_coupler_input
	!use coupler_module
	implicit none

	integer :: infileid
	logical :: found

	!Check all file ids until an unused one is found
	infileid = 100000
	do 
		inquire(unit=infileid,opened=found)
		if (.not.(found)) exit
		infileid = infileid + 1
	enddo

	!Open and read input file on all processes
	open(infileid,file='COUPLER.in',status="old",action="read", &
				  form="formatted")

	call locate(infileid,'DENSITY_CFD',found)
	if (found) then
		read(infileid,*) density_cfd
	else
		call error_abort("Density not specified in coupler input file.")
	end if
	
	call locate(infileid,'OVERLAP_EXTENTS',found)
	if (found) then
		read(infileid,*) icmin_olap
		read(infileid,*) icmax_olap
		read(infileid,*) jcmin_olap
		read(infileid,*) jcmax_olap
		read(infileid,*) kcmin_olap
		read(infileid,*) kcmax_olap
	else
		call error_abort("Ovelap extents unspecified in coupler input file.")
	end if

	call locate(infileid,'TIMESTEP_RATIO',found)
	if (found) then
		read(infileid,*) timestep_ratio !TODO name change
	else
		timestep_ratio = VOID
	end if
	
	call locate(infileid,'MATCH_CELLSIZE',found)
	if (found) then
		read(infileid,*) md_cfd_match_cellsize
	else
		md_cfd_match_cellsize = 0
	end if

	close(infileid,status="keep")

end subroutine read_coupler_input

!------------------------------------------------------------------------------
!                              coupler_cfd_init                               -
!------------------------------------------------------------------------------
!!
!! Initialisation routine for coupler module - Every variable is sent and stored
!! to ensure both md and cfd region have an identical list of parameters
!!
!! - Synopsis
!!
!!  - coupler_cfd_init(nsteps,dt_cfd,icomm_grid,icoord,npxyz_cfd,xyzL,ncxyz,
!!							   density,ijkcmax,ijkcmin,iTmin,iTmax,jTmin,
!!							   jTmax,kTmin,kTmax,xg,yg,zg)
!!
!! - Input
!!
!!  - nsteps
!!   - Number of time steps the CFD code is expected to run for (integer)
!!  - dt_cfd
!!   - CFD timestep (dp real)
!!  - icomm_grid
!!   - The MPI communicator setup by the MPI_CART_CREATE command in the 
!!     CFD region (integer)
!!  - icoord
!!   - The three coordinate for each rank in the domain (integer array nproc by 3)
!!  - npxyz_cfd
!!   - Number of processors in each cartesian dimension (integer array 3)
!!  - xyzL
!!   - Size of domain in each cartesian dimension (dp real array 3)
!!  - ncxyz
!!   - Global number of cells in each cartesian dimension (integer array 3)
!!  - density
!!   - Density of the CFD simulation (dp_real)
!!  - ijkcmax
!!   - Global minimum cell in each cartesian dimension (integer array 3)
!!  - ijkcmin
!!   - Global maximum cell in each cartesian dimension (integer array 3)
!!  - iTmin
!!   - Local minimum cell for each rank (integer array no. procs in x)
!!  - iTmax
!!   - Local maximum cell for each rank (integer array no. procs in x)
!!  - jTmin
!!   - Local minimum cell for each rank (integer array no. procs in y)
!!  - jTmax
!!   - Local maximum cell for each rank (integer array no. procs in y)
!!  - kTmin
!!   - Local minimum cell for each rank (integer array no. procs in z)
!!  - kTmax
!!   - Local maximum cell for each rank (integer array no. procs in z)
!!  - xg
!!   - Array of cell vertices in the x direction (no. cells in x by 
!!     no. cells in y)
!!  - yg
!!   - Array of cell vertices in the y direction (no. cells in x by 
!!     no. cells in y)
!!  - zg
!!   - Array of cell vertices in the z direction (no. cells in z)
!!
!! - Input/Output
!!  - NONE
!!
!! - Output
!!  - NONE
!! 
!! @author Edward Smith
!
! ----------------------------------------------------------------------------


subroutine coupler_cfd_init(nsteps,dt,icomm_grid,icoord,npxyz_cfd,xyzL,ncxyz, & 
							   density,ijkcmax,ijkcmin,iTmin,iTmax,jTmin, & 
							   jTmax,kTmin,kTmax,xgrid,ygrid,zgrid)
    use mpi
 	implicit none			

    integer,					    intent(in)	:: nsteps,icomm_grid 
    integer,dimension(3),		    intent(in)	:: ijkcmin,ijkcmax,npxyz_cfd,ncxyz
	integer,dimension(:),		    intent(in)	:: iTmin,iTmax,jTmin,jTmax,kTmin,kTmax
    integer,dimension(:,:),		    intent(in)	:: icoord
    real(kind(0.d0)),			    intent(in)	:: dt,density
    real(kind(0.d0)),dimension(3),  intent(in)	:: xyzL
    real(kind(0.d0)),dimension(:  ),intent(in)	:: zgrid
    real(kind(0.d0)),dimension(:,:),intent(in)	:: xgrid,ygrid


    integer											:: i,ib,jb,kb,pcoords(3),source,nproc
    integer,dimension(:),allocatable				:: buf,rank_world2rank_realm,rank_world2rank_cart
	real(kind=kind(0.d0))							:: dxmin,dxmax,dzmin,dzmax
    real(kind=kind(0.d0)),dimension(:),allocatable 	:: rbuf

    ! Duplicate grid communicator for coupler use
    call MPI_comm_dup(icomm_grid,CPL_CART_COMM,ierr)
    call MPI_comm_rank(CPL_CART_COMM,myid_cart,ierr) 
    rank_cart = myid_cart + 1; rootid_cart = 0
	!Send only from root processor
    if (myid_realm .eq. rootid_realm ) then
        source=MPI_ROOT
    else
        source=MPI_PROC_NULL
    endif

	! ================ Exchange and store Data ==============================
	! Data is stored to the coupler module with the same name in both realms
	! Note - MPI Broadcast between intercommunicators is only supported by MPI-2

	! ------------------------ Processor Topology ---------------------------
	! Store & Send CFD number of processors
	npx_cfd = npxyz_cfd(1)
	npy_cfd = npxyz_cfd(2)
	npz_cfd = npxyz_cfd(3)
	nproc_cfd = npx_cfd * npy_cfd * npz_cfd
    call MPI_bcast(npxyz_cfd,3,MPI_INTEGER,source,CPL_INTER_COMM,ierr)	!Send

	! Receive & Store MD number of processors
	allocate(buf(3))
    call MPI_bcast(   buf   ,3,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)	!Receive
    npx_md = buf(1)
    npy_md = buf(2)
    npz_md = buf(3)
    nproc_md = npx_md * npy_md * npz_md
	deallocate(buf)

	! Store & Send CFD processor rank to coord
    allocate(rank2coord_cfd(3,nproc_cfd),stat=ierr); rank2coord_cfd = icoord
	iblock_realm=icoord(1,rank_realm); jblock_realm=icoord(2,rank_realm); kblock_realm=icoord(3,rank_realm)
	allocate(buf(3*nproc_cfd)); buf = reshape(icoord, (/ 3*nproc_cfd /) )
    call MPI_bcast(buf,3*nproc_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr)	!Send
	deallocate(buf)

	! Receive & Store MD processor rank to coord
	allocate(buf(3*nproc_md))
    call MPI_bcast(buf,3*nproc_md ,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)	!Receive
    allocate(rank2coord_md (3,nproc_md),stat=ierr); rank2coord_md = reshape(buf,(/ 3,nproc_md /))
	deallocate(buf)

	! Setup CFD mapping from coordinate to rank
	! Store & Send CFD mapping from coordinate to rank to MD
	allocate(coord2rank_cfd(npx_cfd,npy_cfd,npz_cfd))
	do ib = 1,npx_cfd
	do jb = 1,npy_cfd
	do kb = 1,npz_cfd
		pcoords = (/ ib, jb, kb /)-1
		call MPI_Cart_rank(CPL_CART_COMM,pcoords,i,ierr)
		coord2rank_cfd(ib,jb,kb) = i + 1
	enddo
	enddo
	enddo

	allocate(buf(nproc_cfd)); buf = reshape(coord2rank_cfd, (/ nproc_cfd /) )
    call MPI_bcast(coord2rank_cfd,nproc_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr)	!Send
	deallocate(buf)

	! Receive & Store MD coordinate to rank mapping
	allocate(buf(nproc_md))
    call MPI_bcast(buf,nproc_md ,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)	!Receive
    allocate(coord2rank_md (npx_md,npy_md,npz_md)) 
	coord2rank_md = reshape(buf,(/ npx_md,npy_md,npz_md /))
	deallocate(buf)

	! Setup CFD mapping between realm & world rank 
	allocate(rank_cfdrealm2rank_world(nproc_cfd))
	allocate(rank_world2rank_realm(nproc_world))
	call CPL_rank_map(CPL_REALM_COMM,rank_realm,nproc, & 
						rank_cfdrealm2rank_world,rank_world2rank_realm,ierr)

	!World to rank is the same on both realms
	allocate(rank_world2rank_cfdrealm(nproc_world))
	allocate(rank_world2rank_mdrealm(nproc_world))
	rank_world2rank_cfdrealm = rank_world2rank_realm
	rank_world2rank_mdrealm  = rank_world2rank_realm

	! Send CFD mapping to MD
	call MPI_bcast(rank_cfdrealm2rank_world,nproc_cfd,MPI_integer,source,CPL_INTER_COMM,ierr)	 !send

	! Receive & Store MD mapping from realm to local rank from MD
	allocate(rank_mdrealm2rank_world(nproc_md))
	call MPI_bcast(rank_mdrealm2rank_world,nproc_md,MPI_integer,0,CPL_INTER_COMM,ierr)	!Receive

	! Setup CFD mapping between cartesian topology & world rank 
	allocate(rank_cfdcart2rank_world(nproc_cfd))
	allocate(rank_world2rank_cart(nproc_world))
	call CPL_rank_map(CPL_CART_COMM,rank_cart,nproc, & 
						rank_cfdcart2rank_world,rank_world2rank_cart,ierr)

	!World to rank is the same on both realms cart
	allocate(rank_world2rank_cfdcart(nproc_world))
	allocate(rank_world2rank_mdcart(nproc_world))
	rank_world2rank_cfdcart = rank_world2rank_cart
	rank_world2rank_mdcart  = rank_world2rank_cart

	! Send CFD mapping to MD
	call MPI_bcast(rank_cfdcart2rank_world,nproc_cfd,MPI_integer,source,CPL_INTER_COMM,ierr)	 !send

	! Receive & Store MD mapping from cart to local rank from MD
	allocate(rank_mdcart2rank_world(nproc_md))
	call MPI_bcast(rank_mdcart2rank_world,nproc_md,MPI_integer,0,CPL_INTER_COMM,ierr)	!Receive

	! ------------------ Timesteps and iterations ------------------------------
	! Store & send CFD nsteps and dt_cfd
	nsteps_cfd = nsteps
    call MPI_bcast(nsteps,1,MPI_integer,source,CPL_INTER_COMM,ierr)			!Send
	dt_cfd = dt
    call MPI_bcast(dt,1,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	! Receive & store MD timestep dt_md
    call MPI_bcast(dt_md,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive
    call MPI_bcast(nsteps_md,1,MPI_integer,     0,CPL_INTER_COMM,ierr)		!Receive

	! ------------------ Send CFD grid extents ------------------------------

	! Store & send CFD density
	density_cfd = density
	call MPI_bcast(density_cfd,1,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	! Receive & store MD density
	call MPI_bcast(density_md,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive

	! Store & send CFD domain size
	xL_cfd = xyzL(1); yL_cfd = xyzL(2); zL_cfd = xyzL(3)
	call MPI_bcast(xyzL,3,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	! Receive & store MD domain size
	allocate(rbuf(3))
	call MPI_bcast(rbuf,3,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive
	xL_md = rbuf(1); yL_md = rbuf(2); zL_md = rbuf(3);
	deallocate(rbuf)

	! Store & send CFD grid extents
	icmin = ijkcmin(1); jcmin = ijkcmin(2); kcmin = ijkcmin(3)
	icmax = ijkcmax(1); jcmax = ijkcmax(2); kcmax = ijkcmax(3)
    call MPI_bcast((/ icmin,icmax,jcmin,jcmax,kcmin,kcmax /),6,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send

	! Store & send global number of cells in CFD
	ncx = ncxyz(1); ncy = ncxyz(2); ncz = ncxyz(3)
    call MPI_bcast(ncxyz,3,MPI_INTEGER,source,CPL_INTER_COMM,ierr)				!Send

	! Store & send array of global grid points
    allocate(xg(size(xgrid+1,1)+1,size(xgrid,2)+1),stat=ierr); xg = xgrid
    allocate(yg(size(ygrid+1,1)+1,size(ygrid,2)+1),stat=ierr); yg = ygrid
    allocate(zg(size(zgrid+1,1)+1			   ),stat=ierr)  ; zg = zgrid
    call MPI_bcast(xgrid,size(xgrid),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(ygrid,size(ygrid),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(zgrid,size(zgrid),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send

	!call write_matrix(xg,'cfd side, xg=',50+rank_realm)
	!call write_matrix(yg,'cfd side, yg=',50+rank_realm)
	!write(50+rank_realm,*), 'CFD side',rank_realm,'zg',zg

    ! Store & Send local (processor) CFD grid extents
    allocate(icPmin_cfd(npx_cfd),stat=ierr); icPmin_cfd(:) = iTmin(:)
    allocate(icPmax_cfd(npx_cfd),stat=ierr); icPmax_cfd(:) = iTmax(:)
    allocate(jcPmin_cfd(npy_cfd),stat=ierr); jcPmin_cfd(:) = jTmin(:)
    allocate(jcPmax_cfd(npy_cfd),stat=ierr); jcPmax_cfd(:) = jTmax(:)
    allocate(kcPmin_cfd(npz_cfd),stat=ierr); kcPmin_cfd(:) = kTmin(:)
    allocate(kcPmax_cfd(npz_cfd),stat=ierr); kcPmax_cfd(:) = kTmax(:)
    call MPI_bcast(icPmin_cfd,npx_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(icPmax_cfd,npx_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(jcPmin_cfd,npy_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(jcPmax_cfd,npy_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(kcPmin_cfd,npz_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(kcPmax_cfd,npz_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send

	!Calculate the cell sizes dx,dy & dz
	dx = xL_cfd/ncx	  !xg(2,1)-xg(1,1)
	dy = yg(1,2)-yg(1,1) ! yL_cfd/ncy
	dz = zL_cfd/ncz	  !zg(2  )-zg(1  )

	ncx_olap = icmax_olap - icmin_olap + 1
	ncy_olap = jcmax_olap - jcmin_olap + 1
	ncz_olap = kcmax_olap - kcmin_olap + 1

	!Broadcast the overlap to CFD on intracommunicator
	call MPI_bcast(ncy_olap,1,MPI_INTEGER,rootid_realm,CPL_REALM_COMM,ierr)
	!Broadcast the overlap to MD over intercommunicator
	call MPI_bcast(ncy_olap,1,MPI_INTEGER,source,CPL_INTER_COMM,ierr)

	!Check for grid strectching and terminate process if found
	call check_mesh

contains

	subroutine check_mesh
		implicit none

		!Define cell sizes dx,dy & dz and check for grid stretching
		! - - x - -
		dx = xg(2,1)-xg(1,1)
		dxmax = maxval(xg(2:ncx+1,2:ncy+1)-xg(1:ncx,1:ncy))
		dxmin = minval(xg(2:ncx+1,2:ncy+1)-xg(1:ncx,1:ncy))
		if (dxmax-dx.gt.0.00001d0) call error_abort("ERROR - Grid stretching in x not supported")
		if (dx-dxmin.gt.0.00001d0) call error_abort("ERROR - Grid stretching in x not supported")
		! - - y - -
		dy = yg(1,2)-yg(1,1)
		dymax = maxval(yg(2:ncx+1,2:ncy+1)-yg(1:ncx,1:ncy))
		dymin = minval(yg(2:ncx+1,2:ncy+1)-yg(1:ncx,1:ncy))
		if (dymax-dy.gt.0.00001d0) call error_abort("ERROR - Grid stretching in y not supported")
		if (dy-dymin.gt.0.00001d0) call error_abort("ERROR - Grid stretching in y not supported")
		!if (dymax-dy.gt.0.0001 .or. dy-dymin.gt.0.0001) then
	    !    write(*,*) "********************************************************************"
	    !    write(*,*) " Grid stretching employed in CFD domain - range of dy sizes:        "
		!	write(*,*) "dymin = ", dymin, " dy = ",dy, " dymax = ", dymax
	    !    write(*,*) "********************************************************************"
	    !    write(*,*)
		!endif
		! - - z - -

		dzmax = maxval(zg(2:ncz)-zg(1:ncz))
		dzmin = minval(zg(2:ncz)-zg(1:ncz))
		if (dzmax-dz.gt.0.00001d0) call error_abort("ERROR - Grid stretching in z not supported")
		if (dz-dzmin.gt.0.00001d0) call error_abort("ERROR - Grid stretching in z not supported")


	end subroutine check_mesh

end subroutine coupler_cfd_init



!------------------------------------------------------------------------------
!                              coupler_md_init                               -
!------------------------------------------------------------------------------
!
!! Initialisation routine for coupler module - Every variable is sent and stored
!! to ensure both md and cfd region have an identical list of parameters
!!
!! - Synopsis
!!
!!  - coupler_mf_init(nsteps,dt_md,icomm_grid,icoord,npxyz_md,globaldomain,density)
!!
!! - Input
!!
!!  - nsteps
!!   - Number of time steps the MD code is expected to run for (integer)
!!  - dt_md
!!   - MD timestep (dp real)
!!  - icomm_grid
!!   - The MPI communicator setup by the MPI_CART_CREATE command in the 
!!     CFD region (integer)
!!  - icoord
!!   - The three coordinate for each rank in the domain (integer array nproc by 3)
!!  - npxyz_md
!!   - Number of processors in each cartesian dimension (integer array 3)
!!  - globaldomain
!!   - Size of domain in each cartesian dimension (dp real array 3)
!!  - density
!!   - Density of the CFD simulation (dp_real)
!!
!! - Input/Output
!!  - NONE
!!
!! - Output
!!  - NONE
!! 
!! @author Edward Smith
!
! ----------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! Initialisation routine for coupler - Every variable is sent and stored
! to ensure both md and cfd region have an identical list of parameters

subroutine coupler_md_init(nsteps,dt,icomm_grid,icoord,npxyz_md,globaldomain,density)
	use mpi
	implicit none

	integer, intent(in)	  							:: nsteps, icomm_grid
	integer,dimension(3), intent(in)	  			:: npxyz_md	
	integer,dimension(:,:),allocatable,intent(in)	:: icoord
	real(kind(0.d0)),intent(in) 					:: dt,density
    real(kind=kind(0.d0)),dimension(3),intent(in) 	:: globaldomain

    integer											:: i,ib,jb,kb,pcoords(3),source,nproc
    integer,dimension(:),allocatable  				:: buf,rank_world2rank_realm,rank_world2rank_cart
    real(kind=kind(0.d0)),dimension(:),allocatable 	:: rbuf

    ! Duplicate grid communicator for coupler use
    call MPI_comm_dup(icomm_grid,CPL_CART_COMM,ierr)
    call MPI_comm_rank(CPL_CART_COMM,myid_cart,ierr) 
    rank_cart = myid_cart + 1; rootid_cart = 0	
	!Send only from root processor
	if ( myid_realm .eq. rootid_realm ) then
		source=MPI_ROOT
	else
		source=MPI_PROC_NULL
	endif

	! ================ Exchange and store Data ==============================
	! Data is stored to the coupler module with the same name in both realms
	! Note - MPI Broadcast between intercommunicators is only supported by MPI-2

	! ------------------------ Processor Topology ---------------------------
	! Receive & Store CFD number of processors
    allocate(buf(3))
	call MPI_bcast(  buf   ,3,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr) !Receive
	npx_cfd = buf(1); npy_cfd = buf(2); npz_cfd = buf(3)
	nproc_cfd = npx_cfd * npy_cfd * npz_cfd
	deallocate(buf)

	! Store & Send MD number of processors
	npx_md = npxyz_md(1);	npy_md = npxyz_md(2);	npz_md = npxyz_md(3)	
	nproc_md = npx_md * npy_md * npz_md
	call MPI_bcast(npxyz_md,3,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send

	! Receive & Store CFD processor rank to coord
	allocate(buf(3*nproc_cfd))
    call MPI_bcast(buf,3*nproc_cfd,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr) !Receive
    allocate(rank2coord_cfd(3,nproc_cfd),stat=ierr); rank2coord_cfd = reshape(buf,(/ 3,nproc_cfd /))
	deallocate(buf)

	! Store & Send MD processor rank to coord
    allocate(rank2coord_md(3,nproc_md),stat=ierr); rank2coord_md = icoord
	iblock_realm=icoord(1,rank_realm); jblock_realm=icoord(2,rank_realm); kblock_realm=icoord(3,rank_realm)
	allocate(buf(3*nproc_md)); buf = reshape(icoord,(/ 3*nproc_md /))
    call MPI_bcast(buf ,3*nproc_md,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
	deallocate(buf)

	! Receive & Store CFD coordinate to rank mapping
	allocate(buf(nproc_cfd))
    call MPI_bcast(buf,nproc_cfd ,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)	!Receive
    allocate(coord2rank_cfd (npx_cfd,npy_cfd,npz_cfd))
	coord2rank_cfd = reshape(buf,(/ npx_cfd,npy_cfd,npz_cfd /))
	deallocate(buf)

	! Setup MD mapping from coordinate to rank, 
	! Store & Send MD mapping from coordinate to rank to CFD
	allocate(coord2rank_md(npx_md,npy_md,npz_md))
	do ib = 1,npx_md
	do jb = 1,npy_md
	do kb = 1,npz_md
		pcoords = (/ ib, jb, kb /)-1
		call MPI_Cart_rank(CPL_CART_COMM,pcoords,i,ierr)
		coord2rank_md(ib,jb,kb) = i + 1
	enddo
	enddo
	enddo
	allocate(buf(nproc_md)); buf = reshape(coord2rank_md, (/ nproc_md /) )
    call MPI_bcast(coord2rank_md,nproc_md,MPI_INTEGER,source,CPL_INTER_COMM,ierr)	!Send
	deallocate(buf)

	! Setup MD mapping between realm & world rank 
	allocate(rank_mdrealm2rank_world(nproc_md))
	allocate(rank_world2rank_realm(nproc_world))
	call CPL_rank_map(CPL_REALM_COMM,rank_realm,nproc, & 
						rank_mdrealm2rank_world,rank_world2rank_realm,ierr)

	!World to rank is the same on both realms
	allocate(rank_world2rank_cfdrealm(nproc_world))
	allocate(rank_world2rank_mdrealm(nproc_world))
	rank_world2rank_cfdrealm = rank_world2rank_realm
	rank_world2rank_mdrealm  = rank_world2rank_realm

	! Receive & Store CFD_mapping from realm to local rank
	allocate(rank_cfdrealm2rank_world(nproc_cfd))
	call MPI_bcast(rank_cfdrealm2rank_world,nproc_cfd,MPI_integer,0,CPL_INTER_COMM,ierr)	!Receive

	! Send MD mapping from realm to local rank to CFD
	call MPI_bcast(rank_mdrealm2rank_world,nproc_md,MPI_integer,source,CPL_INTER_COMM,ierr)	 !send

	! Setup MD mapping between cartesian topology & world rank 
	allocate(rank_mdcart2rank_world(nproc_md))
	allocate(rank_world2rank_cart(nproc_world))
	call CPL_rank_map(CPL_CART_COMM,rank_cart,nproc, & 
						rank_mdcart2rank_world,rank_world2rank_cart,ierr)

	!World to rank is the same on both realms cart
	allocate(rank_world2rank_cfdcart(nproc_world))
	allocate(rank_world2rank_mdcart(nproc_world))
	rank_world2rank_cfdcart = rank_world2rank_cart
	rank_world2rank_mdcart  = rank_world2rank_cart

	! Receive & Store CFD_mapping from cart to local rank
	allocate(rank_cfdcart2rank_world(nproc_cfd))
	call MPI_bcast(rank_cfdcart2rank_world,nproc_cfd,MPI_integer,0,CPL_INTER_COMM,ierr)	!Receive

	! Send MD mapping from cart to local rank to CFD
	call MPI_bcast(rank_mdcart2rank_world,nproc_md,MPI_integer,source,CPL_INTER_COMM,ierr)	 !send

	! ------------------ Timesteps and iterations ------------------------------
	! Receive & store CFD nsteps and dt_cfd
	call MPI_bcast(nsteps_cfd,1,MPI_integer,0,CPL_INTER_COMM,ierr)				!Receive
	call MPI_bcast(dt_cfd,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive

	! Store & send MD timestep to dt_md
	dt_MD = dt
    call MPI_bcast(dt,1,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send
	nsteps_MD = nsteps
    call MPI_bcast(nsteps,1,MPI_integer,source,CPL_INTER_COMM,ierr)	!Send

	! ------------------ Receive CFD grid extents ------------------------------
	! Receive & store CFD density
	call MPI_bcast(density_cfd,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive

	! Store & send MD density
	density_md = density 
	call MPI_bcast(density,1,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	! Receive & store CFD domain size
	allocate(rbuf(3))
	call MPI_bcast(rbuf,3,MPI_double_precision,0,CPL_INTER_COMM,ierr)				!Receive
	xL_cfd = rbuf(1); yL_cfd = rbuf(2); zL_cfd = rbuf(3)
	deallocate(rbuf)

	! Store & send MD domain size
	xL_md = globaldomain(1); yL_md = globaldomain(2); zL_md = globaldomain(3) 
	call MPI_bcast(globaldomain,3,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	! Receive & Store global CFD grid extents
	allocate(buf(6))
	call MPI_bcast(buf, 6, MPI_INTEGER, 0, CPL_INTER_COMM,ierr) !Send
	icmin = buf(1); icmax = buf(2)
	jcmin = buf(3); jcmax = buf(4)
	kcmin = buf(5); kcmax = buf(6)
	deallocate(buf)

	! Receive & Store array of global number of cells in CFD
	allocate(buf(3))
	call MPI_bcast(buf, 3, MPI_INTEGER, 0, CPL_INTER_COMM,ierr) !Receive
	ncx = buf(1); ncy = buf(2); ncz = buf(3)
	deallocate(buf)

	! Receive & Store array of global grid points
	allocate(xg(ncx+1,ncy+1),yg(ncx+1,ncy+1),zg(ncz+1))
	call MPI_bcast(xg,size(xg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive
	call MPI_bcast(yg,size(yg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive
	call MPI_bcast(zg,size(zg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive

	! Receive & Store local (processor) CFD grid extents
    allocate(icPmin_cfd(npx_cfd)); allocate(icPmax_cfd(npx_cfd));  
	allocate(jcPmin_cfd(npy_cfd)); allocate(jcPmax_cfd(npy_cfd));
    allocate(kcPmin_cfd(npz_cfd)); allocate(kcPmax_cfd(npz_cfd));
    call MPI_bcast(icPmin_cfd,npx_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(icPmax_cfd,npx_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(jcPmin_cfd,npy_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(jcPmax_cfd,npy_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(kcPmin_cfd,npz_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(kcPmax_cfd,npz_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive

	!Calculate the cell sizes dx,dy & dz
	dx = xL_cfd/ncx	  !xg(2,1)-xg(1,1)
	dy = yg(1,2)-yg(1,1) ! yL_cfd/ncy
	dz = zL_cfd/ncz	  !zg(2  )-zg(1  )

	ncx_olap = icmax_olap - icmin_olap + 1
	ncy_olap = jcmax_olap - jcmin_olap + 1
	ncz_olap = kcmax_olap - kcmin_olap + 1

    ! Initialise other md module variables if data is provided in coupler.in
!	if (md_average_period_tag == CPL) then 
!		average_period = md_average_period
!	endif
!	if (md_save_period_tag == CPL) then
!		save_period    = md_save_period
!	endif

	!if (cfd_code_id_tag == CPL) then
	!	cfd_code_id  = cfd_code_id_in
	!endif
    
	if ( nsteps_md <= 0 ) then 
		write(0,*) "Number of MD steps per dt interval <= 0"
		write(0,*) "Coupler will not work, quitting ..."
		call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_INIT,ierr)
	endif

	!Receive overlap from the CFD
    !call MPI_bcast(ncy_olap,1,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive

	! ================ Apply domain setup  ==============================
	! --- set the sizes of the MD domain ---
	!Fix xL_md domain size to continuum
	!xL_md = x(icmax_cfd) - x(icmin_cfd)

    ! yL_md is adjusted to an integer number of initialisation cells in the following steps
   ! if (md_ly_extension_tag == CPL) then 
    !    DY_PURE_MD = md_ly_extension
    !else 
    !    DY_PURE_MD = y(jcmin) - y(jcmino)
    !end if
    !yL_md = y(jcmax_overlap) - y(jcmino) + DY_PURE_MD
    !yL_md = real(floor(yL_md/b),kind(0.d0))*b
    !DY_PURE_MD = yL_md - (y(jcmax_overlap) - y(jcmino))

	!Fix zL_md domain size to continuum
	!zL_md = z(kcmax_cfd) - z(kcmin)

end subroutine coupler_md_init



!-----------------------------------------------------------------------------
!	THIS ROUTINE SHOULD BE PART OF THE INITIALISATION OF THE COUPLER
!	WHERE THE INITIAL TIMES (STIME FOR CFD & ELAPSEDTIME FOR MD) ARE
!	COMPARED TO SEE IF RESTART IS CONSISTENT AND NUMBER OF STEPS ON
!	BOTH SIDES OF THE COUPLER ARE CALCULATED AND STORED
! Setup ratio of CFD to MD timing and total number of timesteps

subroutine set_coupled_timing(Nsteps,initialstep,elapsedtime)
	!use coupler_module, only : timestep_ratio,dt_cfd,dt_MD, &
	!                           Nsteps_cfd,Nsteps_md_old=>Nsteps_md, & 
	!						   rank_realm,Nsteps_coupled
	implicit none

	integer,intent(in)					:: initialstep
	integer,intent(out)					:: Nsteps
	real(kind=kind(0.d0)),intent(out)	:: elapsedtime

	integer								:: Nsteps_MDperCFD

	!Set number of MD timesteps per CFD using ratio of timestep or coupler value
	if(timestep_ratio .eq. VOID) then
		Nsteps_MDperCFD = int(dt_cfd/dt_MD)
	else 
		Nsteps_MDperCFD = timestep_ratio
	endif
	Nsteps_coupled = Nsteps_cfd

 	!Set number of steps in MD simulation
	Nsteps_md   = initialstep + Nsteps_cfd * Nsteps_MDperCFD
	elapsedtime = elapsedtime + Nsteps_cfd * Nsteps_MD * dt_MD

	if (rank_realm .eq. 1) then 
		write(*,'(2(a,/),a,i7,a,i7,/a,i7,a,i7,/a,f12.4,a,/a)'), &
			"*********************************************************************", 		&
			" WARNING - WARNING - WARNING - WARNING - WARNING - WARNING - WARNING  ", 		&
			" Input number of timesteps from MD: ",Nsteps," & CFD: ", Nsteps_cfd,	&
			" is set in this coupled run to MD: ", Nsteps_md, ",CFD/Coupled: ", Nsteps_cfd,	&
			" At the end of the run, elapsed time will be: ", elapsedtime, " LJ time units ", 		&
			"*********************************************************************"   
	endif 

	!Set corrected nsteps returned to MD
	Nsteps = Nsteps_md

end subroutine set_coupled_timing


!=============================================================================
!! Establish for all MD processors the mapping (if any) 
!! to coupled CFD processors
!-----------------------------------------------------------------------------

subroutine CPL_create_map
	use mpi
	!use coupler_module
	implicit none

	! Check (as best as one can) that the inputs will work
	call check_config_feasibility

	! Get ranges of cells on each MD processor
	call get_md_cell_ranges

	! Get overlapping mapping for MD to CFD
	call get_overlap_blocks

	! Setup overlap communicators
	call prepare_overlap_comms

	! Setup graph topology
	call CPL_overlap_topology

contains

subroutine check_config_feasibility
	implicit none

	integer :: ival
	character(len=256) :: string

	! Check there is only one overlap CFD proc in y
	ival = nint( dble(ncy) / dble(npy_cfd) )
	if (ncy_olap .gt. ival) then

	    string = "This coupler will not work if there is more than one "// &
		         "CFD (y-coordinate) in the overlapping region. "       // &
		         "Aborting simulation."

		call error_abort(string)

	end if


	! Check whether ncx,ncy,ncz are an integer multiple of npx_md, etc.
	! - N.B. no need to check ncy/npy_md.
	ival = 0
	ival = ival + mod(ncx,npx_cfd)
	ival = ival + mod(ncy,npy_cfd)
	ival = ival + mod(ncz,npz_cfd)
	ival = ival + mod(ncx,npx_md)
	ival = ival + mod(ncz,npz_md)
	if (ival.ne.0) then	
		string = "The number of cells in the cfd domain is not an "    // &
		         "integer multiple of the number of processors in "    // &
		         "the x and z directions. Aborting simulation."	
		call error_abort(string)
	end if

	! Check that the MD region is large enough to cover overlap
	ival = 0
	ival = ival + floor(xL_olap/xL_md)
	ival = ival + floor(yL_olap/yL_md)
	ival = ival + floor(zL_olap/zL_md)
	if (ival.ne.0) then	
		string = "Overlap region is larger than the MD region. "       // &
		         "Aborting simulation."
		call error_abort(string)
	end if


	! Check overlap cells are within CFD extents
	ival = 0
	if (icmin_olap.lt.icmin) ival = ival + 1		
	if (icmax_olap.gt.icmax) ival = ival + 1		
	if (jcmin_olap.lt.jcmin) ival = ival + 1		
	if (jcmax_olap.gt.jcmax) ival = ival + 1		
	if (kcmin_olap.lt.kcmin) ival = ival + 1		
	if (kcmax_olap.gt.kcmax) ival = ival + 1		
	if (ival.ne.0) then
		string = "Overlap region has been specified outisde of the "  // &
		         "CFD region. Aborting simulation."
		call error_abort(string)
	end if
	
		
	! Check MD/CFD ratios are integers in x and z
	if (mod(npx_md,npx_cfd) .ne. 0) then

		print'(a,i8,a,i8)', ' number of MD processors in x ', npx_md,     & 
							' number of CFD processors in x ', npx_cfd

		call error_abort("get_overlap_blocks error - number of MD "    // & 
						 "processors in x must be an integer multiple "// &
						 "of number of CFD processors in x")

	elseif (mod(npz_md,npz_cfd) .ne. 0) then

		print'(a,i8,a,i8)', ' number of MD processors in z ', npz_md,     & 
							' number of CFD processors in z ', npz_cfd

		call error_abort("get_overlap_blocks error - number of MD "    // &
						 "processors in z must be an integer multiple "// &
						 "of number of CFD processors in z")

	endif

end subroutine check_config_feasibility


!------------------------------------------------------------
!Calculate processor cell ranges of MD code on all processors
	
subroutine get_md_cell_ranges
	!use coupler_module
	implicit none

	integer :: n
	integer :: olap_jmin_mdcoord
	integer :: ncxl, ncyl, nczl
	integer :: ncy_mdonly, ncy_md, ncyP_md

	allocate(icPmin_md(npx_md)); icPmin_md = VOID
	allocate(jcPmin_md(npy_md)); jcPmin_md = VOID
	allocate(kcPmin_md(npz_md)); kcPmin_md = VOID
	allocate(icPmax_md(npx_md)); icPmax_md = VOID
	allocate(jcPmax_md(npy_md)); jcPmax_md = VOID
	allocate(kcPmax_md(npz_md)); kcPmax_md = VOID

	! - - x - -
	ncxl = ceiling(dble(ncx)/dble(npx_md))
	do n=1,npx_md
		icPmax_md(n) = n * ncxl
		icPmin_md(n) = icPmax_md(n) - ncxl + 1
	end do	

	! - - y - -
	ncy_md   = nint(yL_md/dy)
	ncy_mdonly = ncy_md - ncy_olap
	ncyP_md = ncy_md / npy_md
	olap_jmin_mdcoord = npy_md - floor(dble(ncy_olap)/dble(ncyP_md))	 
	do n = olap_jmin_mdcoord,npy_md
		jcPmax_md(n) = n * ncyP_md - ncy_mdonly
		jcPmin_md(n) = jcPmax_md(n) - ncyP_md + 1
		if (jcPmin_md(n).le.0) jcPmin_md(n) = 1
	end do  

	! - - z - -
	nczl = ceiling(dble(ncz)/dble(npz_md))
	do n=1,npz_md
		kcPmax_md(n) = n * nczl
		kcPmin_md(n) = kcPmax_md(n) - nczl + 1
	end do

end subroutine get_md_cell_ranges

!------------------------------------------------------------
!Calculate processor overlap between CFD/MD on all processors

subroutine get_overlap_blocks
	!use coupler_module
	implicit none

	integer 			:: i,n,endproc,nolapsx,nolapsy,nolapsz
	integer             :: xLl_md, yLl_md, zLl_md
	integer				:: yLl_cfd
	integer,dimension(3):: pcoords

	xL_olap = ncx_olap * dx 
	yL_olap = ncy_olap * dy 
	zL_olap = ncz_olap * dz 

	xLl_md  = xL_md / npx_md
	yLl_md  = yL_md / npy_md
	zLl_md  = zL_md / npz_md

	if (realm .eq. md_realm) then
		xLl = xLl_md; yLl = yLl_md ; zLl = zLl_md 
	endif

	nolapsx = nint( dble( npx_md ) / dble( npx_cfd ) )
	nolapsy = ceiling ( yL_olap / yLl_md ) 
	nolapsz = nint( dble( npz_md ) / dble( npz_cfd ) )

	!Get cartesian coordinate of overlapping md cells & cfd cells
	allocate(cfd_icoord2olap_md_icoords(npx_cfd,nolapsx)) 
	allocate(cfd_jcoord2olap_md_jcoords(npy_cfd,nolapsy)) 
	allocate(cfd_kcoord2olap_md_kcoords(npz_cfd,nolapsz)) 
	cfd_icoord2olap_md_icoords = VOID
	cfd_jcoord2olap_md_jcoords = VOID
	cfd_kcoord2olap_md_kcoords = VOID

	! - - x - -
	do n = 1,npx_cfd
	do i = 1,nolapsx	
		cfd_icoord2olap_md_icoords(n,i) = (n-1)*nolapsx + i
	end do
	end do

	! - - y - -
	yLl_cfd = yL_cfd/npy_cfd
	endproc = ceiling(yL_olap/yLl_cfd)
	do n = 1,endproc
	do i = 1,nolapsy
		cfd_jcoord2olap_md_jcoords(n,i) =   (n-1)*nolapsy + i &
										  + (npy_md - nolapsy)
	end do
	end do

	! - - z - -
	do n = 1,npz_cfd
	do i = 1,nolapsz	
		cfd_kcoord2olap_md_kcoords(n,i) = (n-1)*nolapsz + i
	end do
	end do

end subroutine get_overlap_blocks

!subroutine intersect_comm
!	use coupler_module
!	use mpi
!	implicit none

!	integer :: n
!	integer,dimension(3)   :: pcoords

	!Create communicator for all intersecting processors
!	if (realm .eq. cfd_realm) then
!		map%n = npx_md/npx_cfd
!		allocate(map%rank_list(map%n))
		!Get rank(s) of overlapping MD processor(s)
!		do n = 1,map%n
!			pcoords(1)=cfd_icoord2olap_md_icoords(rank2coord_cfd(1,rank_realm),n)
!			pcoords(2)=cfd_jcoord2olap_md_jcoords(rank2coord_cfd(2,rank_realm),n)
!			pcoords(3)=cfd_kcoord2olap_md_kcoords(rank2coord_cfd(3,rank_realm),n)
!			if (any(pcoords(:).eq.VOID)) then
!				map%n = 0; map%rank_list(:) = VOID
!			else
!				map%rank_list(n) = coord2rank_md(pcoords(1),pcoords(2),pcoords(3))
!			endif
!			write(250+rank_realm,'(2a,6i5)'), 'overlap',realm_name(realm),rank_realm,map%n,map%rank_list(n),pcoords
!		enddo
!	else if (realm .eq. md_realm) then
!		map%n = 1
!		allocate(map%rank_list(map%n))	
		!Get rank of overlapping CFD processor
!		pcoords(1) = rank2coord_md(1,rank_realm)*(dble(npx_cfd)/dble(npx_md))
!		pcoords(2) = npy_cfd !rank2coord_md(2,rank_realm)*(dble(npy_cfd)/dble(npy_md))
!		pcoords(3) = rank2coord_md(3,rank_realm)*(dble(npz_cfd)/dble(npz_md))
!!		map%rank_list(1) = coord2rank_cfd(pcoords(1),pcoords(2),pcoords(3))
!		write(300+rank_realm,'(2a,6i5)'), 'overlap',realm_name(realm),rank_realm,map%n,map%rank_list(1),pcoords
!	endif

!end subroutine intersect_comm

!=========================================================================

subroutine prepare_overlap_comms
	!use coupler_module
	use mpi
	implicit none

	!loop over cfd cart ranks
	! find cfd cart coords from cfd cart rank
	!  find md cart coords (from cfd_icoord2olap_md_jcoords)
	!   find md cart rank from md cart coords (coord2rank_md) 
	!    find md world rank from md cart rank (rank_mdcart2rank_world)
	!     set group(md_world_rank) to cfd cart rank
	!split comm to groups
	!if group(world_rank) == 0, set olap_comm to null 

	integer :: i,j,k,ic,jc,kc
	integer :: trank_md, trank_cfd, trank_world, nolap
	integer, dimension(:), allocatable :: mdicoords, mdjcoords, mdkcoords
	integer, parameter :: olap_null = -666
	integer :: group(nproc_world), rank_olap2rank_realm_temp(nproc_world)
	integer :: cfdcoord(3)
	integer :: tempsize

	tempsize = size(cfd_icoord2olap_md_icoords,2)
	allocate(mdicoords(tempsize))
	tempsize = size(cfd_jcoord2olap_md_jcoords,2)
	allocate(mdjcoords(tempsize))
	tempsize = size(cfd_kcoord2olap_md_kcoords,2)
	allocate(mdkcoords(tempsize))
	
	allocate(olap_mask(nproc_world))

	!Set default values, must be done because coord2rank_md cannot
	!take "null" coordinates.
	group(:) = olap_null
	olap_mask(:) = .false.
	nolap = 0

	! Every process loop over all cfd ranks
	do trank_cfd = 1,nproc_cfd

		! Get cart coords of cfd rank
		cfdcoord(:)  = rank2coord_cfd(:,trank_cfd)

		! Get md cart coords overlapping cfd proc
		mdicoords(:) = cfd_icoord2olap_md_icoords(cfdcoord(1),:)
		mdjcoords(:) = cfd_jcoord2olap_md_jcoords(cfdcoord(2),:)
		mdkcoords(:) = cfd_kcoord2olap_md_kcoords(cfdcoord(3),:)

		! Set group and olap_mask for CFD processor if it overlaps
		if (any(mdicoords.ne.olap_null) .and. &
			any(mdjcoords.ne.olap_null) .and. &
			any(mdkcoords.ne.olap_null)) then

			trank_world = rank_cfdcart2rank_world(trank_cfd)
			olap_mask(trank_world) = .true.
			group    (trank_world) = trank_cfd

		end if

		! Set group and olap_mask for MD processors
		do i = 1,size(mdicoords)
		do j = 1,size(mdjcoords)
		do k = 1,size(mdkcoords)

			ic = mdicoords(i)
			jc = mdjcoords(j)
			kc = mdkcoords(k)

			if (any((/ic,jc,kc/).eq.olap_null)) cycle

			trank_md = coord2rank_md(ic,jc,kc)
			trank_world = rank_mdcart2rank_world(trank_md)

			olap_mask(trank_world) = .true.
			group    (trank_world) = trank_cfd
			
		end do
		end do	
		end do

	end do

	! Split world Comm into a set of comms for overlapping processors
	call MPI_comm_split(CPL_WORLD_COMM,group(rank_world),realm, &
	                    CPL_OLAP_COMM,ierr)

	!Setup Overlap comm sizes and id
	if (realm.eq.cfd_realm) CFDid_olap = myid_olap
	call MPI_bcast(CFDid_olap,1,MPI_INTEGER,CFDid_olap,CPL_OLAP_COMM,ierr)

	! USED ONLY FOR OUTPUT/TESTING??
	!if (myid_olap .eq. CFDid_olap) testval = group(rank_world)
	!call MPI_bcast(testval,1,MPI_INTEGER,CFDid_olap,CPL_OLAP_COMM,ierr)

	! Set all non-overlapping processors to MPI_COMM_NULL
	if (olap_mask(rank_world).eq..false.) then
		myid_olap = olap_null
		rank_olap = olap_null
		CPL_OLAP_COMM = MPI_COMM_NULL
	end if

	!Setup overlap map
	call CPL_rank_map(CPL_OLAP_COMM,rank_olap,nproc_olap, & 
	                  rank_olap2rank_world,rank_world2rank_olap,ierr)
	myid_olap = rank_olap - 1

	deallocate(mdicoords)
	deallocate(mdjcoords)
	deallocate(mdkcoords)
	
	!if (realm.eq.md_realm) call write_overlap_comms_md

end subroutine prepare_overlap_comms

!=========================================================================
!Setup topology graph of overlaps between CFD & MD processors

subroutine CPL_overlap_topology
	use mpi
	implicit none

	integer								:: i, n, nneighbors, nconnections
	integer, dimension(:),allocatable	:: index, edges, id_neighbors
	logical								:: reorder

	!Allow optimisations of ordering
	reorder = .true.

	!Get number of processors in communicating overlap region 
	if (olap_mask(rank_world).eq..true.) then

		!CFD processor is root and has mapping to all MD processors
		allocate(index(nproc_olap))			! Index for each processor
		allocate(edges(2*(nproc_olap)-1))	! nproc_olap-1 for CFD and one for
		                                    ! each of nproc_olap MD processors
		index = 0; 	edges = 0

		!CFD processor has connections to nproc_olap MD processors
		nconnections = nproc_olap-1
		index(1) = nconnections
		do n = 1,nconnections
			edges(n) = n !olap_list(n+1) !CFD connected to all MD processors 1 to nconnections
		enddo

		!MD processor has a single connection to CFD
		nconnections = 1; i = 2
		do n = nproc_olap+1,2*(nproc_olap)-1
			index(i) = index(i-1) + nconnections !Each successive index incremented by one
			edges(n) = CFDid_olap	!Connected to CFD processor
			i = i + 1
		enddo

		!Create graph topology for overlap region
		call MPI_Graph_create(CPL_OLAP_COMM, nproc_olap,index,edges,reorder,CPL_GRAPH_COMM,ierr)
	else
		CPL_GRAPH_COMM = MPI_COMM_NULL
	endif

	! Setup graph map
	call CPL_rank_map(CPL_GRAPH_COMM,rank_graph,nproc_olap, & 
					 rank_graph2rank_world,rank_world2rank_graph,ierr)
	myid_graph = rank_graph - 1

end subroutine CPL_overlap_topology

end subroutine CPL_create_map

!=============================================================================
!	Adjust CFD domain size to an integer number of lattice units used by  
!	MD if sizes are given in sigma units
!-----------------------------------------------------------------------------

subroutine CPL_cfd_adjust_domain(xL, yL, zL, nx, ny, nz, density_output)
    use mpi
    !use coupler_module, only : density_cfd,CPL_REALM_COMM, rank_realm, ierr
    implicit none

    integer, optional, intent(inout) 			:: nx, ny, nz
    real(kind(0.d0)), optional, intent(inout) 	:: xL,yL,zL
    real(kind(0.d0)), optional, intent(inout)  	:: density_output

    ! Internal variables
    integer										:: ierror, root
	character(1)				   				:: direction
	logical										:: changed

	!Define root processes
	root = 1

    density_output = density_cfd

	! Check CFD domain and MD domain are compatible sizes to allow a
	! stable initial MD lattice structure - resize if possible or
	! stop code and demand a regeneration of grid if vary by more than 0.01
	changed = .false.
    if (present(xL)) then
        call init_length(xL,resize=.true.,direction='x', &
		                 print_warning=changed)
    endif

	! No need to adjust y because we can adjust DY in MD to
	! have an integer number of FCC units. ??????? What

    if (present(zL)) then
        call init_length(zL,resize=.true.,direction='z', &
		                 print_warning=changed)
    endif

	if ( changed ) then
		print*, "Regenerate Grid with corrected sizes as above"
		call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
	endif

    ! check id CFD cell sizes are larger than 2*sigma 
    call test_cfd_cell_sizes

contains

!-----------------------------------------------------------------------------

subroutine init_length(rout,resize,direction,print_warning)
	!use coupler_module, only: dx,dy,dz,error_abort
	implicit none
            
    real(kind=kind(0.d0)), intent(inout) :: rout
    logical, intent(in)                  :: resize
	character(*),intent(in)              :: direction
    logical,intent(out)                  :: print_warning

	real(kind(0.d0)) :: dxyz  ! dx, dy or dz
    real(kind(0.d0)) :: rinit ! initial val of rout or rin for print

    print_warning=.false.

	select case (direction)
	case('x','X')
		dxyz = dx
	case('y','Y')
		dxyz = dy
	case('z','Z')
		dxyz = dz
	case default
		call error_abort('Wrong direction specified in init_length')
	end select

	if ( resize ) then

		rinit = rout
		rout = real(nint(rout/dxyz),kind(0.d0))*dxyz
		print_warning = .true.
		print*, direction, 'dxyz = ', dxyz 

	endif

    if (print_warning) then 

        !if (rank_realm .eq. root) then 

            write(*,'(3(a,/),3a,/,2(a,f20.10),/a,/,a)') &
                    "*********************************************************************",	&
                    "WARNING - this is a coupled run which resets CFD domain size         ",	&
                    " to an integer number of MD initial cells:		                      ", 	&
                    "	Domain resized in the ", direction, " direction			          ",	&
                    " inital size =", rinit, " resized ", rout,									&
                    "								                                      ",	& 
                    "*********************************************************************"   
        !endif

		!If resize is insignificant then return flag print_warning as false
		if (abs(rinit-rout) .lt. 0.01) print_warning = .false.

    end if
    
end subroutine init_length

!-----------------------------------------------------------------------------

subroutine test_cfd_cell_sizes
    implicit none

    if (rank_realm .eq. root) then
        if (present(xL) .and. present(nx)) then
            if (xL/nx < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in x direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          xL=",xL,"nx=",nx
            endif
        endif

        if (present(yL) .and. present(ny)) then
            if (yL/ny < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in y direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          yL=",yL,"nx=",ny
            endif
        endif

        if (present(zL) .and. present(nz)) then
            if (zL/nz < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in z direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          zL=",zL,"nx=",nz
            endif
        endif
    end if
        
end subroutine test_cfd_cell_sizes
            
end subroutine CPL_cfd_adjust_domain



!-------------------------------------------------------------------
! 					CPL_rank_map								   -
!-------------------------------------------------------------------

! Get COMM map for current communicator and relationship to 
! world rank used to link to others in the coupler hierachy

! - - - Synopsis - - -

! CPL_rank_map(COMM, rank, comm2world, world2comm, ierr)

! - - - Input Parameters - - -

!comm
!    communicator with cartesian structure (handle) 

! - - - Output Parameter - - -

!rank
!    rank of a process within group of comm (integer)
!    NOTE - fortran convention rank=1 to nproc  
!nproc
!    number of processes within group of comm (integer) 
!comm2world
!	Array of size nproc_world which for element at 
!	world_rank has local rank in COMM
!world2comm
!	Array of size nproc_COMM which for element at 
!	for local ranks in COMM has world rank 
!ierr
!    error flag


subroutine CPL_rank_map(COMM,rank,nproc,comm2world,world2comm,ierr)
	!use coupler_module, only : rank_world, nproc_world, CPL_WORLD_COMM, VOID
	use mpi
	implicit none

	integer, intent(in)								:: COMM
	integer, intent(out)							:: rank,nproc,ierr
	integer, dimension(:),allocatable,intent(out)	:: comm2world,world2comm

	allocate(world2comm( nproc_world))
	world2comm( nproc_world) = VOID

	if (COMM .ne. MPI_COMM_NULL) then

		!Mapping from comm rank to world rank
		call MPI_comm_rank(COMM,rank,ierr)
		rank = rank + 1
		call MPI_comm_size(COMM,nproc,ierr)
		allocate(comm2world(nproc))
		call MPI_allgather(rank_world,1,MPI_INTEGER, & 
						   comm2world,1,MPI_INTEGER,COMM,ierr)
	else
		rank = VOID
		allocate(comm2world(0))
	endif

	!Mapping from world rank to comm rank
	call MPI_allgather(rank      ,1,MPI_INTEGER, & 
					   world2comm,1,MPI_INTEGER,CPL_WORLD_COMM,ierr)

end subroutine CPL_rank_map

!---------------------------------------------------
! Locate file in input

subroutine locate(fileid,keyword,have_data)
	implicit none
	
	integer,intent(in)			:: fileid               ! File unit number
	character(len=*),intent(in)	:: keyword              ! Input keyword	
	logical,intent(out)			:: have_data            ! Flag: input found

	character*(100)				:: linestring           ! First 100 chars
	integer						:: keyword_length       ! Length of keyword
	integer						:: io                   ! File status flag

	keyword_length = len(keyword)
	rewind(fileid)
	
	! Loop until end of file or keyword found
	do
		! Read first 100 characters of line
		read (fileid,'(a)',iostat=io) linestring

		! If end of file is reached, exit
		if (io.ne.0) then 
			have_data = .false.
			exit
		end if
		
		! If the first characters match keyword, exit
		if (linestring(1:keyword_length).eq.keyword) then
			have_data = .true.
			exit
		endif

	end do

end subroutine locate

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
!subroutine request_stop(tag)
!    use mpi
!    implicit none
!
!    character(len=*),intent(in) ::tag
!    integer myid, ierr
!
!    ! do nothing, get out quick 
!    if(.not. stop_request_activated ) return
!
!    if (tag /= stop_request_name) return
!
!    select case(stop_request_name)
!    case("create_comm","CREATE_COMM")
!        call mpi_comm_rank(CPL_REALM_COMM, myid,ierr)
!        write(0,*) 'stop as requested at ', trim(stop_request_name), ', realm',realm, 'rank', myid
!        call MPI_Finalize(ierr)
!        stop
!    case("create_map","CREATE_MAP")
!        call mpi_comm_rank(CPL_REALM_COMM, myid,ierr)
!        write(0,*) 'stop as requested at ', trim(stop_request_name), ', realm',realm, 'rank', myid
!        call MPI_Finalize(ierr)
!        stop    
!    case default
!        write(0,*) "WARNING: request abort activated, but the tag is unrecognized, check COUPLER.in"
!        write(0,*) "         accepted stop tags are: create_comm"
!    end select
!
!end subroutine request_stop


end module coupler_module
