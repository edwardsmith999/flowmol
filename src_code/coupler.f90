!=============================================================================
!
!				  Coupler 
!
! Routines accessible from application ( molecular or continuum ) after 
! the name, in parenthesis, is the realm in which each routine must be called
!
! SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP
!
! coupler_create_comm	      (cfd+md)   splits MPI_COMM_WORLD, create inter - 
!										 communicator between CFD and MD
! coupler_create_map	      (cfd+md)   creates correspondence maps between 
!										 the CFD grid and MD domains
! coupler_cfd_init              (cfd)    initialises coupler with CFD data
! coupler_md_init               (cfd)    initialises coupler and set MD 
!										 parameters with using data from CFD 
!										 or COUPLER.in
! coupler_cfd_adjust_domain     (cfd)    adjust CFD tomain to an integer number 
!										 FCC or similar MD initial layout
!
! SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION
!
! coupler_send_data        	  (cfd+md)   sends grid data exchanged between 
!										 realms ( generic interface)
! coupler_recv_data        	  (cfd+md)   receives data exchanged between realms 
!										 ( generic interface)
! coupler_cfd_get               (cfd)    returns coupler internal parameters 
!										 for CFD realm
! coupler_md_get                 (md)    returns coupler internal parameters 
!										 for MD realm
! coupler_md_get_save_period     (md)    auxiliary used for testing
! coupler_md_get_average_period  (md)    returns average period of BC
! coupler_md_get_md_per_cfd_dt   (md) 	 returns the number of step MD does for 
										 !each CFD step
! coupler_md_get_nsteps          (md)    returm CFD nsteps  
! coupler_md_get_dt_cfd          (md)    returns MD dt
! coupler_md_set                 (md)    sets zL if CFD is 2D
! coupler_md_get_density         (md)    gets CFD density
! coupler_md_get_cfd_id          (md)    id for CFD code, possible values set 
!										 in coupler_parameters
!
!  Lucian Anton, November 2011
!  Revised April 2012
!
!=============================================================================

module coupler
	use coupler_parameters
    implicit none
    save

    interface coupler_send_data
        module procedure coupler_send_data_3d, coupler_send_data_4d
    end interface

    interface coupler_recv_data
        module procedure coupler_recv_data_3d, coupler_recv_data_4d
    end interface

    private coupler_send_data_3d, coupler_send_data_4d, &
        coupler_send_data_xd, coupler_recv_data_3d, coupler_recv_data_4d,&
        coupler_recv_data_xd

contains

!=============================================================================
! Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup
!
!								SETUP
!
! Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup
!=============================================================================

!=============================================================================
! 							coupler_create_comm	      	
! (cfd+md) Splits MPI_COMM_WORLD in both the CFD and MD code respectively
! 		   and create intercommunicator between CFD and MD
!-----------------------------------------------------------------------------

subroutine coupler_create_comm(realm, REALM_COMM, ierror)
	use mpi
	use coupler_internal_common
	use coupler_input_data
	implicit none

	integer, intent(in) :: realm ! CFD or MD
	integer, intent(out):: REALM_COMM, ierror

	! test if we have a CFD and a MD realm
	ierror=0

	call test_realms			! Test realms are assigned correctly
	COUPLER_REALM = realm
	call create_comm			! Create intercommunicator between realms

contains

!-----------------------------------------------------------------------------
!	Test if CFD and MD realms are assigned correctly
!-----------------------------------------------------------------------------

subroutine test_realms
	use mpi
	implicit none

	integer 			 :: i, myid, root, nproc, ncfd, nmd, ierr
	integer, allocatable :: realm_list(:)

	!Get processor id of all processor across both codes
	call MPI_comm_rank(MPI_COMM_WORLD,myid,ierr)

	! Allocate and gather array with realm (MD or CFD) of each 
	! processor in the coupled domain on the root processor
	root = 0
	if (myid .eq. root) then
		call MPI_comm_size(MPI_comm_world, nproc, ierr)
		allocate(realm_list(nproc))
	endif
	call MPI_gather(realm,1,MPI_INTEGER,realm_list,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	!Check through array of processors on both realms
	!and return error if wrong values or either is missing
	if (myid .eq. root) then
		ncfd = 0; nmd = 0
		do i =1, nproc
			if ( realm_list(i) .eq. COUPLER_CFD ) then 
				ncfd = ncfd + 1
			else if ( realm_list(i) .eq. COUPLER_MD ) then
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
	use mpi
	use coupler_internal_common, only : COUPLER_REALM_COMM, COUPLER_GLOBAL_COMM, COUPLER_ICOMM
	implicit none

	integer ierr, myid, myid_comm, myid_comm_max, realm, &
		ibuf(2), jbuf(2), remote_leader, comm, comm_size

	realm = COUPLER_REALM
	! Split MPI COMM WORLD ready to establish two communicators
	! 1) A global intra-communicator in each realm for communication
	! internally between CFD processes or between MD processes
	! 2) An inter-communicator which allows communication between  
	! the 'groups' of processors in MD and the group in the CFD 
	call MPI_comm_dup(MPI_COMM_WORLD,COUPLER_GLOBAL_COMM,ierr)
	call MPI_comm_rank(COUPLER_GLOBAL_COMM,myid,ierr)
	REALM_COMM         = MPI_COMM_NULL
	COUPLER_REALM_COMM = MPI_COMM_NULL

	!------------ create realm intra-communicators -----------------------
	! Split MPI_COMM_WORLD into an intra-communicator for each realm 
	! (used for any communication within each realm - e.g. broadcast from 
	!  an md process to all other md processes)
	call MPI_comm_split(MPI_COMM_WORLD,realm,myid,REALM_COMM,ierr)

	!------------ create realm inter-communicators -----------------------
	! Create intercommunicator between the group of processor on each realm
	! (used for any communication between realms - e.g. md group rank 2 sends
	! to cfd group rank 5). inter-communication is by a single processor on each group
	! Split duplicate of MPI_COMM_WORLD
	call MPI_comm_split(COUPLER_GLOBAL_COMM,realm,myid,COUPLER_REALM_COMM,ierr)

	! Get the MPI_comm_world ranks that hold the largest ranks in cfd_comm and md_comm
	call MPI_comm_rank(COUPLER_REALM_COMM,myid_comm,ierr)
	call MPI_comm_size(COUPLER_REALM_COMM,comm_size,ierr)
	ibuf(:) = -1
	jbuf(:) = -1
	if ( myid_comm .eq. comm_size - 1) then
		ibuf(COUPLER_REALM) = myid
	endif
	call MPI_allreduce( ibuf ,jbuf, 2, MPI_INTEGER, MPI_MAX, &
						COUPLER_GLOBAL_COMM, ierr)

	!Set this largest rank on each process to be the inter-communicators
	select case (COUPLER_REALM)
	case (COUPLER_CFD)
		remote_leader = jbuf(COUPLER_MD)
	case (COUPLER_MD)
		remote_leader = jbuf(COUPLER_CFD)
	end select

	!print*,color, jbuf, remote_leader

	call MPI_intercomm_create(COUPLER_REALM_COMM, comm_size - 1, COUPLER_GLOBAL_COMM,&
									remote_leader, 1, COUPLER_ICOMM, ierr)

	write(0,*) 'did (inter)communicators ', code_name(COUPLER_REALM), myid

end subroutine create_comm

end subroutine coupler_create_comm


!=============================================================================
!	Setup within each of the CFD/MD the mapping to the other realm
! --- CFD ---
! create_map_cfd is an internal CFD only routine to build the mapping
! only on CFD processors to the MD regions
! ---  MD  ---
! create_map_md is an internal MD only routine to build the mapping
! only on MD processors to the CFD regions
!-----------------------------------------------------------------------------

subroutine coupler_create_map
	use mpi
	use coupler_internal_cfd
	use coupler_internal_md
    use coupler_internal_common, only : request_stop
	implicit none

	integer ierr

	if (COUPLER_REALM .eq. COUPLER_CFD) then
		call create_map_cfd
	else if (COUPLER_REALM .eq. COUPLER_MD) then
		call create_map_md
	else
		write(*,*) "Wrong COUPLER_REALM in coupler_create_map"
		call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_REALM,ierr)
	end if

    call request_stop("create_map")
 
end subroutine coupler_create_map

!=============================================================================
! Get MD processor topology and timestep details on all CFD processors 
! and send mesh details
!-----------------------------------------------------------------------------

subroutine coupler_cfd_init(icomm_grid, imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,&
							kmax,kmaxo,nsteps,x,y,z,dx,dz,npx,npy,npz,icoord,dt)
	use mpi
	use coupler_internal_common
	use coupler_input_data
	use coupler_internal_cfd, only : imino_ => imino, imin_ => imin, imax_ => imax, imaxo_ => imaxo, &
        jmino_ => jmino, jmin_ => jmin, jmax_ => jmax, jmaxo_ => jmaxo, kmino_ => kmino, kmin_ => kmin, &
        kmax_ => kmax, kmaxo_ => kmaxo, nsteps_ => nsteps, x_ => x, y_ => y, z_ => z, dx_ => dx, dz_ => dz, &
		npx_ => npx, npy_ => npy, npz_ => npz, icoord_ => icoord, dt_ => dt, jmax_overlap, &
		npx_md, npy_md, npz_md, nproc_md, MD_initial_cellsize
	implicit none

	integer, intent(in) :: icomm_grid,imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,kmax,kmaxo
	integer, intent(in) :: nsteps,npx,npy,npz,icoord(:,:)
	real(kind(0.d0)), intent(in) :: x(:),y(:),z(:), dx, dz, dt

	integer i, myid, source, buf(6),ierr

    ! duplicate grid communicator for coupler use
    call MPI_comm_dup(icomm_grid,coupler_grid_comm,ierr)


	!Store copies of passed variables in CFD internal
	imino_ = imino; imin_ = imin; imax_ = imax; imaxo_ = imaxo;  
	jmino_ = jmino; jmin_ = jmin; jmax_ = jmax; jmaxo_ = jmaxo; 
	kmino_ = kmino; kmin_ = kmin; kmax_ = kmax; kmaxo_ = kmaxo;
	nsteps_ = nsteps;  dx_ = dx; dz_ = dz
	npx_ = npx; npy_ = npy; npz_ = npz 

    if(kmax == kmin) then
        cfd_is_2d = .true.
    endif

	allocate(x_(size(x)),stat=ierr); x_ = x
	allocate(y_(size(y)),stat=ierr); y_ = y
	allocate(z_(size(z)),stat=ierr); z_ = z
	allocate(icoord_(3,npx*npy*npz),stat=ierr)
	icoord_ = icoord

    call MPI_comm_rank(COUPLER_REALM_COMM,myid,ierr)

    ! Test if MD_init_cell size is larger than CFD cell size
    if( (MD_initial_cellsize .ge. x(2)-x(1) .or. MD_initial_cellsize .ge. y(2)-y(1)) .and. myid .eq. 0 ) then
        write(*,*) 
        write(*,*) "********************************************************************"
        write(*,*) " WARNING ...WARNING ...WARNING ...WARNING ...WARNING ...WARNING ... "
        write(*,*) " MD initialisation cell size larger than CFD x,y cell sizes         "
        write(*,*) " MD_init_cellsize=",MD_initial_cellsize
        write(*,*) " x(2)-x(1)=",x(2)-x(1), " y(2)-y(1)=",y(2)-y(1) 
        write(*,*) "********************************************************************"
        write(*,*)
    endif
  
	! send CFD processor grid and overlap parameter

	! Note: jmax_overlap default is provided in coupler_internal_cfd
    if (cfd_coupler_input%overlap%tag == CPL) then
        jmax_overlap =  jmin + cfd_coupler_input%overlap%y_overlap
    endif

	if ( myid .eq. 0 ) then
		source=MPI_ROOT
	else
		source=MPI_PROC_NULL
	endif


	!Note - MPI Broadcast between intercommunicators is only supported by MPI-2
	call MPI_bcast((/ npx, npy, npz, jmax_overlap /), 4, MPI_INTEGER,&
						source, COUPLER_ICOMM,ierr)

	! receive MD processor grid 
	call MPI_bcast(buf, 3, MPI_INTEGER,0, COUPLER_ICOMM,ierr)

	npx_md = buf(1)
	npy_md = buf(2)
	npz_md = buf(3)
	nproc_md = npx_md * npy_md * npz_md

	! send CFD mesh data
	call MPI_bcast((/ imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,kmax,kmaxo /), 12,&
						MPI_INTEGER, source, COUPLER_ICOMM,ierr)
	call MPI_bcast(x,size(x),MPI_double_precision,source,COUPLER_ICOMM,ierr)
	call MPI_bcast(y,size(y),MPI_double_precision,source,COUPLER_ICOMM,ierr)
	call MPI_bcast(z,size(z),MPI_double_precision,source,COUPLER_ICOMM,ierr)
	call MPI_bcast((/ dx, dz /),2,MPI_double_precision,source,COUPLER_ICOMM,ierr)

	! send CFD nsteps and dt
	call MPI_bcast(nsteps,1,MPI_integer,source,COUPLER_ICOMM,ierr)
	call MPI_bcast( (/ dt, density, MD_initial_cellsize /),3,MPI_double_precision,source,COUPLER_ICOMM,ierr)

	! write(0,*)' CFD: did exchange grid data'

end subroutine coupler_cfd_init




! Initialisation routine for coupler - Every variable is sent and stored
! to ensure both md and cfd region have an identical list of parameters

subroutine coupler_cfd_init_es(nsteps,dt_cfd,icomm_grid,npxyz_cfd,icoord,xyzL,ngxyz, & 
							   ijkmax,ijkmin,ijkTmax,ijkTmin,xpg,ypg,zpg)
    use mpi
    use coupler_internal_common
    use coupler_input_data
    use coupler_module,	dt_cfd_=>dt_cfd,nsteps_=>nsteps,					&	!Simulation lengths	
									xpg_=>xpg,ypg_=>ypg,zpg_=>zpg				!CFD grid arrays
	implicit none			
    !use coupler_module, only: 		dt_md,dt_cfd_=>dt_cfd,			&	!Timesteps
	!								nsteps_=>nsteps,						&	!Simulation lengths	
	!								cfd_xL,cfd_yL,cfd_zL,					&	!CFD domain length
	!								md_xL,md_yL,md_zL,						&	!MD domain length
	!								ngx, ngy, ngz,							&	!CFD global no. of cells
	!								imin,imax,jmin,jmax,kmin,kmax, 			& 	!CFD global grid limits
	!								iTmin,iTmax,jTmin,jTmax,kTmin,kTmax, 	& 	!CFD local grid limits
	!								xpg_=>xpg,ypg_=>ypg,zpg_=>zpg, 			& 	!CFD grid arrays
	!								npx_cfd,npy_cfd,npz_cfd,nproc_cfd, 		&	!CFD processors
     !   							npx_md ,npy_md ,npz_md ,nproc_md, 		& 	!MD processors
	!								icoord_cfd , icoord_md, 				&	!Processor topologies
	!								jmax_overlap, MD_initial_cellsize


    integer,					   intent(in):: icomm_grid, nsteps
    integer,dimension(3),		   intent(in):: ijkmin,ijkmax,npxyz_cfd,ngxyz
    integer,dimension(:,:),		   intent(in):: ijkTmin, ijkTmax, icoord
    real(kind(0.d0)),			   intent(in):: dt_cfd
    real(kind(0.d0)),dimension(3), intent(in):: xyzL
    real(kind(0.d0)),dimension(:  ),allocatable,intent(in):: zpg
    real(kind(0.d0)),dimension(:,:),allocatable,intent(in):: xpg,ypg

    integer											:: i, myid, source, ierr
    integer,dimension(:),allocatable				:: buf
    real(kind=kind(0.d0)),dimension(:),allocatable 	:: rbuf

    ! Duplicate grid communicator for coupler use
    call MPI_comm_dup(icomm_grid,coupler_grid_comm,ierr)
    call MPI_comm_rank(COUPLER_REALM_COMM,myid,ierr)
	!Send only from root processor
    if ( myid .eq. 0 ) then
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
    call MPI_bcast(npxyz_cfd,3,MPI_INTEGER,source,COUPLER_ICOMM,ierr)	!Send

	! Receive & Store MD number of processors
	allocate(buf(3))
    call MPI_bcast(   buf   ,3,MPI_INTEGER,  0   ,COUPLER_ICOMM,ierr)	!Receive
    npx_md = buf(1)
    npy_md = buf(2)
    npz_md = buf(3)
    nproc_md = npx_md * npy_md * npz_md
	deallocate(buf)

	! Store & Send CFD processor topology
    allocate(icoord_cfd(3,nproc_cfd),stat=ierr); icoord_cfd = icoord
	allocate(buf(3*nproc_cfd)); buf = reshape(icoord, (/ 3*nproc_cfd /) )
    call MPI_bcast(icoord_cfd,3*nproc_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr)	!Send

	! Receive & Store MD processor topology
    allocate(icoord_md (3,nproc_md),stat=ierr)
    call MPI_bcast(buf,3*nproc_md ,MPI_INTEGER,  0   ,COUPLER_ICOMM,ierr)	!Receive
	icoord_md = reshape(buf,(/ 3,nproc_md /))
	deallocate(buf)

	! ------------------ Timesteps and iterations ------------------------------
	! Store & send CFD nsteps and dt_cfd
	nsteps_ = nsteps
    call MPI_bcast(nsteps,1,MPI_integer,source,COUPLER_ICOMM,ierr)			!Send
	dt_cfd_ = dt_cfd
    call MPI_bcast(dt_cfd,1,MPI_double_precision,source,COUPLER_ICOMM,ierr)	!Send

	! Receive & store MD timestep dt_md
    call MPI_bcast(dt_md,1,MPI_double_precision,0,COUPLER_ICOMM,ierr)		!Receive

	! ------------------ Send CFD grid extents ------------------------------
	! Store & send CFD domain size
	cfd_xL = xyzL(1); cfd_yL = xyzL(2); cfd_zL = xyzL(3)
	call MPI_bcast(xyzL,3,MPI_double_precision,source,COUPLER_ICOMM,ierr)	!Send

	! Receive & store MD domain size
	allocate(rbuf(3))
	call MPI_bcast(rbuf,3,MPI_double_precision,0,COUPLER_ICOMM,ierr)		!Receive
	md_xL = rbuf(1); md_yL = rbuf(2); md_zL = rbuf(3);
	deallocate(rbuf)

	! Store & send CFD grid extents
	imin = ijkmin(1); jmin = ijkmin(2); kmin = ijkmin(3)
	imax = ijkmax(1); jmax = ijkmax(2); kmin = ijkmax(3)
    call MPI_bcast((/ imin,imax,jmin,jmax,kmin,kmax /),6,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send

	! Store & send global number of cells in CFD
	ngx = ngxyz(1); ngy = ngxyz(2); ngz = ngxyz(3)
    call MPI_bcast(ngxyz,3,MPI_INTEGER,source,COUPLER_ICOMM,ierr)				!Send

	! Store & send array of global grid points
    allocate(xpg_(size(xpg,1),size(xpg,2)),stat=ierr); xpg_ = xpg
    allocate(ypg_(size(ypg,1),size(ypg,2)),stat=ierr); ypg_ = ypg
    allocate(zpg_(size(zpg,1)			 ),stat=ierr); zpg_ = zpg
    call MPI_bcast(xpg,size(xpg),MPI_double_precision,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(ypg,size(ypg),MPI_double_precision,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(zpg,size(zpg),MPI_double_precision,source,COUPLER_ICOMM,ierr) !Send

    ! Store & Send local (processor) CFD grid extents
    allocate(iTmin(nproc_cfd),stat=ierr); iTmin(:) = ijkTmin(1,:)
    allocate(iTmax(nproc_cfd),stat=ierr); iTmax(:) = ijkTmax(1,:)
    allocate(jTmin(nproc_cfd),stat=ierr); jTmin(:) = ijkTmin(2,:)
    allocate(jTmax(nproc_cfd),stat=ierr); jTmax(:) = ijkTmax(2,:)
    allocate(kTmin(nproc_cfd),stat=ierr); kTmin(:) = ijkTmin(3,:)
    allocate(kTmax(nproc_cfd),stat=ierr); kTmax(:) = ijkTmax(3,:) 
    call MPI_bcast(iTmin,nproc_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(iTmax,nproc_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(jTmin,nproc_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(jTmax,nproc_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(kTmin,nproc_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(kTmax,nproc_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send

    ! send CFD processor grid and overlap parameter
    ! Note: jmax_overlap default is provided in coupler_internal_cfd
    if (cfd_coupler_input%overlap%tag == CPL) then
        jmax_overlap =  jmin + cfd_coupler_input%overlap%y_overlap
    endif

    ! test if MD_init_cell size is larger than CFD cell size
    if( (MD_initial_cellsize >= xpg(2,1)-xpg(1,1) .or. & 
		 MD_initial_cellsize >= ypg(1,2)-ypg(1,1) .or. & 
		 MD_initial_cellsize >= zpg( 2 )-zpg( 1 )).and. myid == 0 ) then
        write(*,*)
        write(*,*) "********************************************************************"
        write(*,*) " WARNING ...WARNING ...WARNING ...WARNING ...WARNING ...WARNING ... "
        write(*,*) " MD initialisation cell size larger than CFD x,y cell sizes         "
        write(*,*) " MD_init_cellsize=",MD_initial_cellsize
        write(*,*) " dx=",xpg(2,1)-xpg(1,1), " dy=",ypg(1,2)-ypg(1,1)," dz=",zpg(2)-zpg(1 )
        write(*,*) "********************************************************************"
        write(*,*)
    endif

end subroutine coupler_cfd_init_es



subroutine coupler_md_init_es(dt_md,icomm_grid,icoord,globaldomain,npxyz_md)
	use mpi
	use coupler_internal_common
	use coupler_input_data, cfd_code_id_in => cfd_code_id
    use coupler_module, dt_md_=>dt_md
									
   ! use coupler_module, only: 		dt_md_=>dt_md,dt_cfd,nsteps,					&	!Timesteps
	!								cfd_xL,cfd_yL,cfd_zL,					&	!CFD domain length
	!								md_xL,md_yL,md_zL,						&	!MD domain length
	!								ngx, ngy, ngz,							&	!CFD global no. of cells
	!								imin,imax,jmin,jmax,kmin,kmax, 			& 	!CFD global grid limits
	!								iTmin,iTmax,jTmin,jTmax,kTmin,kTmax, 	& 	!CFD local grid limits
	!								xpg_=>xpg,ypg_=>ypg,zpg_=>zpg, 			& 	!CFD grid arrays
	!								npx_cfd,npy_cfd,npz_cfd,nproc_cfd, 		&	!CFD processors
    !    							npx_md ,npy_md ,npz_md ,nproc_md, 		& 	!MD processors
	!								icoord_cfd , icoord_md, 				&	!Processor topologies
	!								jmax_overlap, MD_initial_cellsize
	implicit none

	integer, intent(in)	  							:: icomm_grid
	integer,dimension(3), intent(in)	  			:: npxyz_md	
	integer,dimension(:,:),allocatable,intent(in)	:: icoord
	real(kind(0.d0)),intent(in) 					:: dt_md
    real(kind=kind(0.d0)),dimension(3),intent(in) 	:: globaldomain

    integer											:: i, myid,myid_grid, source, ierr
    integer,dimension(:),allocatable  				:: buf
    real(kind=kind(0.d0)),dimension(:),allocatable 	:: rbuf

    ! Duplicate grid communicator for coupler use
	call MPI_comm_rank(COUPLER_REALM_COMM,myid,ierr)
    call MPI_comm_dup(icomm_grid,coupler_grid_comm,ierr)
    call MPI_comm_rank(coupler_grid_comm,myid_grid,ierr) 
    myid_grid = myid_grid + 1
	!Send only from root processor
	if ( myid .eq. 0 ) then
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
	call MPI_bcast(  buf   ,3,MPI_INTEGER,  0   ,COUPLER_ICOMM,ierr) !Receive
	npx_cfd = buf(1); npy_cfd = buf(2); npz_cfd = buf(3)
	nproc_cfd = npx_cfd * npy_cfd * npz_cfd
	deallocate(buf)

	! Store & Send MD number of processors
	npx_md = npxyz_md(1);	npy_md = npxyz_md(2);	npz_md = npxyz_md(3)	
	nproc_md = npx_md * npy_md * npz_md
	call MPI_bcast(npxyz_md,3,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send

	! Receive & Store CFD processor topology
	allocate(buf(3*nproc_cfd))
    allocate(icoord_cfd(3,nproc_cfd),stat=ierr)
    call MPI_bcast(buf,3*nproc_md,MPI_INTEGER,  0   ,COUPLER_ICOMM,ierr) !Receive
	icoord_cfd = reshape(buf,(/ 3,nproc_cfd /))

	! Store & Send MD processor topology
    allocate(icoord_md(3,nproc_md),stat=ierr); icoord_md = icoord
	buf = reshape(icoord,(/ 3*nproc_md /))
    call MPI_bcast(icoord_md ,3*nproc_md,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
	deallocate(buf)

	! ------------------ Timesteps and iterations ------------------------------
	! Receive & store CFD nsteps and dt_cfd
	call MPI_bcast(nsteps,1,MPI_integer,0,COUPLER_ICOMM,ierr)				!Receive
	call MPI_bcast(dt_cfd,1,MPI_double_precision,0,COUPLER_ICOMM,ierr)		!Receive

	! Store & send MD timestep to dt_md
	dt_md_ = dt_md
    call MPI_bcast(dt_md,1,MPI_double_precision,source,COUPLER_ICOMM,ierr)	!Send

	! ------------------ Receive CFD grid extents ------------------------------
	! Receive & store CFD domain size
	allocate(rbuf(3))
	call MPI_bcast(rbuf,3,MPI_double_precision,0,COUPLER_ICOMM,ierr)				!Receive
	cfd_xL = rbuf(1); cfd_yL = rbuf(2); cfd_zL = rbuf(3)
	deallocate(rbuf)

	! Store & send MD domain size
	md_xL = globaldomain(1); md_yL = globaldomain(2); md_zL = globaldomain(3) 
	call MPI_bcast(globaldomain,3,MPI_double_precision,source,COUPLER_ICOMM,ierr)	!Send

	! Receive & Store global CFD grid extents
	allocate(buf(6))
	call MPI_bcast(buf, 6, MPI_INTEGER, 0, COUPLER_ICOMM,ierr) !Send
	imin = buf(1); imax = buf(2)
	jmin = buf(3); jmax = buf(4)
	kmin = buf(5); kmax = buf(6)
	deallocate(buf)

	! Receive & Store array of global number of cells in CFD
	allocate(buf(3))
	call MPI_bcast(buf, 3, MPI_INTEGER, 0, COUPLER_ICOMM,ierr) !Receive
	ngx = buf(1); ngy = buf(2); ngz = buf(3)
	deallocate(buf)		

	! Receive & Store array of global grid points
	allocate(xpg(ngx,ngy),ypg(ngx,ngy),zpg(ngz))
	call MPI_bcast(xpg,size(xpg),MPI_double_precision,0,COUPLER_ICOMM,ierr) !Receive
	call MPI_bcast(ypg,size(ypg),MPI_double_precision,0,COUPLER_ICOMM,ierr) !Receive
	call MPI_bcast(zpg,size(zpg),MPI_double_precision,0,COUPLER_ICOMM,ierr) !Receive

	! Receive & Store local (processor) CFD grid extents
    allocate(iTmin(nproc_cfd)); allocate(jTmin(nproc_cfd));
    allocate(kTmin(nproc_cfd)); allocate(iTmax(nproc_cfd)); 
    allocate(jTmax(nproc_cfd)); allocate(kTmax(nproc_cfd)); 
    call MPI_bcast(iTmin,nproc_cfd,MPI_INTEGER,0,COUPLER_ICOMM,ierr) !Receive
    call MPI_bcast(iTmax,nproc_cfd,MPI_INTEGER,0,COUPLER_ICOMM,ierr) !Receive
    call MPI_bcast(jTmin,nproc_cfd,MPI_INTEGER,0,COUPLER_ICOMM,ierr) !Receive
    call MPI_bcast(jTmax,nproc_cfd,MPI_INTEGER,0,COUPLER_ICOMM,ierr) !Receive
    call MPI_bcast(kTmin,nproc_cfd,MPI_INTEGER,0,COUPLER_ICOMM,ierr) !Receive
    call MPI_bcast(kTmax,nproc_cfd,MPI_INTEGER,0,COUPLER_ICOMM,ierr) !Receive

	! ------------------ Apply domain setup etc -------------------
	! --- set the sizes of the MD domain ---
	!Fix xL_md domain size to continuum
	!xL_md = x(imax_cfd) - x(imin_cfd)

    ! yL_md is adjusted to an integer number of initialisation cells in the following steps
   ! if (md_ly_extension_tag == CPL) then 
    !    DY_PURE_MD = md_ly_extension
    !else 
    !    DY_PURE_MD = y(jmin_cfd) - y(jmino)
    !end if
    !yL_md = y(jmax_overlap_cfd) - y(jmino) + DY_PURE_MD
    !yL_md = real(floor(yL_md/b),kind(0.d0))*b
    !DY_PURE_MD = yL_md - (y(jmax_overlap_cfd) - y(jmino))

	!Fix zL_md domain size to continuum
	!zL_md = z(kmax_cfd) - z(kmin_cfd)

    !if( kmin_cfd == kmax_cfd) then
    !    cfd_is_2d = .true.
    !endif

    ! initialise other md module variables if data is provided in coupler.in
    !if (md_average_period_tag == CPL) then 
    !    average_period = md_average_period
   ! endif
   ! if (md_save_period_tag == CPL) then
    !    save_period    = md_save_period
    !endif

    !if (cfd_code_id_tag == CPL) then
   !     cfd_code_id  = cfd_code_id_in
    !endif
    
    !if(md_steps_per_dt_cfd_tag == CPL) then
    !    md_steps = md_steps_per_dt_cfd
    !else 
    !    md_steps = int(dt_cfd/dt_MD)
    !endif
    !if ( md_steps <= 0 ) then 
    !    write(0,*) "Number of MD steps per dt interval <= 0"
    !    write(0,*) "Coupler will not work, quitting ..."
    !    call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_INIT,ierr)
    !endif

end subroutine coupler_md_init_es


subroutine coupler_md_init(npxin,npyin,npzin,icoordin,icomm_grid,dtin)
	use mpi
	use coupler_internal_common
	use coupler_input_data, cfd_code_id_in => cfd_code_id
	use coupler_internal_md, b => MD_initial_cellsize, md_density => density, &
        md_steps => md_steps_per_dt_cfd
	implicit none

	integer, intent(in)	  :: npxin, npyin, npzin, icoordin(:,:),icomm_grid
	real(kind(0.d0)), intent(in) :: dtin

	integer i, ierr, source, buf(12)
	real(kind=kind(0.d0)) ra(3)

    ! Duplicate grid communicator for coupler use
	call MPI_comm_rank(COUPLER_REALM_COMM,myid,ierr)
    call MPI_comm_dup(icomm_grid,coupler_grid_comm,ierr)
    call MPI_comm_rank(coupler_grid_comm,myid_grid,ierr) 
    myid_grid = myid_grid + 1
	!Send only from root processor
	if ( myid .eq. 0 ) then
		source=MPI_ROOT
	else
		source=MPI_PROC_NULL
	endif

	! ------------------ Exchange and store Data -------------------
	npx = npxin; npy = npyin; npz = npzin
	nproc = npx * npy * npz

	dt_MD = dtin

	allocate(icoord(3,nproc),stat=ierr)
	icoord = icoordin

	! write(0,*) 'MD exchange grid data'

	! get CFD processor grid and the number of block in j direction
	call MPI_bcast( buf, 4, MPI_INTEGER,0, COUPLER_ICOMM,ierr)
	npx_cfd = buf(1)
	npy_cfd = buf(2)
	npz_cfd = buf(3)
	jmax_overlap_cfd = buf(4)
	nproc_cfd = npx_cfd * npy_cfd * npz_cfd



	call MPI_bcast((/ npx, npy, npz /), 3, MPI_INTEGER,&
								source, COUPLER_ICOMM,ierr)

	!! Test
	!	call MPI_comm_rank(MD_COMM,myid,ierr)!
	!
	!	write(0,*) 'exchange_grid_data: MD side', myid, bbox_cfd%xbb(1:2,1:npx_cfd),bbox_cfd%zbb(1:2,1:npz_cfd),&
	!	 bbox_md%xbb(1:2,1:npx_md),bbox_md%zbb(1:2,1:npz_md)    

	! Receive and unpack CFD mesh data 
	call MPI_bcast(buf, 12, MPI_INTEGER, 0, COUPLER_ICOMM,ierr) 
	imino = buf(1); imin_cfd = buf(2);  imax_cfd = buf(3);  imaxo = buf(4)
	jmino = buf(5); jmin_cfd = buf(6);  jmax_cfd = buf(7);  jmaxo = buf(8)
	kmino = buf(9); kmin_cfd = buf(10); kmax_cfd = buf(11); kmaxo = buf(12)
	allocate (x(imino:imaxo),y(jmino:jmaxo),z(kmino:kmaxo))

	call MPI_bcast(x,size(x),MPI_double_precision,0,COUPLER_ICOMM,ierr)
	call MPI_bcast(y,size(y),MPI_double_precision,0,COUPLER_ICOMM,ierr)
	call MPI_bcast(z,size(z),MPI_double_precision,0,COUPLER_ICOMM,ierr)
	call MPI_bcast(ra,2,MPI_double_precision,0,COUPLER_ICOMM,ierr)

	! rescale all lengths to MD units
	x = fsig * x; y = fsig * y; z = fsig * z 
	dx = fsig * ra(1); dz = fsig * ra(2)

	!write(0,*) 'MD exchange grid data: imin0, imin ...', imino, imin_cfd, imax_cfd,imaxo, &
	! x(imino),x(imin_cfd), x(imax_cfd), x(imaxo)
	!write(0,*) 'MD exchage grid data: recv dx, dz ', dx, dz, jmax_overlap_cfd

	! Get CFD nsteps
	call MPI_bcast(nsteps,1,MPI_integer,0,COUPLER_ICOMM,ierr)

	! Get CFD dt
	call MPI_bcast(ra,3,MPI_double_precision,0,COUPLER_ICOMM,ierr)


	! ------------------ Apply domain setup etc -------------------

	! should dt_CFD be scaled ?  see to it later
	dt_CFD     = ra(1) * FoP_time_ratio
	md_density = ra(2)
	b          = ra(3)

	! --- set the sizes of the MD domain ---
	!Fix xL_md domain size to continuum
	xL_md = x(imax_cfd) - x(imin_cfd)

    ! yL_md is adjusted to an integer number of initialisation cells in the following steps
    if (md_ly_extension_tag == CPL) then 
        DY_PURE_MD = md_ly_extension
    else 
        DY_PURE_MD = y(jmin_cfd) - y(jmino)
    end if
    yL_md = y(jmax_overlap_cfd) - y(jmino) + DY_PURE_MD
    yL_md = real(floor(yL_md/b),kind(0.d0))*b
    DY_PURE_MD = yL_md - (y(jmax_overlap_cfd) - y(jmino))

	!write(0,*) 'MD domain etc', DY_PURE_MD,y(jmin_cfd),yL_md,y(jmax_overlap_cfd),y(jmino)

	!Fix zL_md domain size to continuum
	zL_md = z(kmax_cfd) - z(kmin_cfd)

    if( kmin_cfd == kmax_cfd) then
        cfd_is_2d = .true.
    endif

	!write(0,*) 'MD: exchange_grid... xL_md, yL_md, zL_md',&
	! 				& myid, xL_md, yL_md, zL_md
	!write(0,*) 'MD: nsteps (CFD) ', nsteps

    ! initialise other md module variables if data is provided in coupler.in
    if (md_average_period_tag == CPL) then 
        average_period = md_average_period
    endif
    if (md_save_period_tag == CPL) then
        save_period    = md_save_period
    endif

    if (cfd_code_id_tag == CPL) then
        cfd_code_id  = cfd_code_id_in
    endif
    
    if(md_steps_per_dt_cfd_tag == CPL) then
        md_steps = md_steps_per_dt_cfd
    else 
        md_steps = int(dt_cfd/dt_MD)
    endif
    if ( md_steps <= 0 ) then 
        write(0,*) "Number of MD steps per dt interval <= 0"
        write(0,*) "Coupler will not work, quitting ..."
        call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_INIT,ierr)
    endif

end subroutine coupler_md_init

!=============================================================================
!	Adjust CFD domain size to an integer number of lattice units used by  
!	MD if sizes are given in sigma units
!-----------------------------------------------------------------------------

subroutine coupler_cfd_adjust_domain(xL, yL, zL, nx, ny, nz, density_cfd)
    use mpi
    use coupler_internal_common
    use coupler_input_data
    use coupler_internal_cfd, only : b => MD_initial_cellsize
    implicit none

    integer, optional, intent(inout) 			:: nx, ny, nz
    real(kind(0.d0)),optional, intent(inout) 	:: xL,yL,zL
    real(kind(0.d0)), optional, intent(inout) 	:: density_cfd

    ! Internal variables
    integer										:: myid, ierror, ierr
    real(kind=kind(0.d0)), pointer 				:: xyz_ptr => null()
	character(1)				   				:: direction
	logical										:: changed

    ! Local rank, useful for messeges
    call MPI_comm_rank(COUPLER_REALM_COMM,myid,ierr)

    if (density_tag == CPL) then 
        density_cfd = density
    else
        density     = density_cfd
        density_tag = CFD
    endif

	! Coupler input parameters are set in CFD 
    if ( cfd_coupler_input%domain%tag == VOID) then
        cfd_coupler_input%domain%tag = CFD
	end if

	print*, 'cell type =', cfd_coupler_input%domain%cell_type

    select case (cfd_coupler_input%domain%cell_type)
    case ("FCC","Fcc","fcc")
		b = (4.d0/density)**(1.d0/3.d0)
    case default
		write(*,*) "Wrong unit cell type in coupler_cfd_adjust_domain. Stopping ... "
		ierror = COUPLER_ERROR_INIT
 		call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
    end select

	! Check CFD domain and MD domain are compatible sizes to allow a
	! stable initial MD lattice structure - resize if possible or
	! stop code and demand a regeneration of grid if vary by more than 0.01
	changed = .false.
    if (present(xL)) then
 		xyz_ptr => cfd_coupler_input%domain%x
        call init_length(xyz_ptr,xL,resize=.true.,direction='x',print_warning=changed)
    endif

    if (present(yL)) then
		! No need to adjust y because we can adjust DY in MD to
		! have an integer number of FCC units.
		xyz_ptr => cfd_coupler_input%domain%y
		call init_length(xyz_ptr,yL,resize=.false.,direction='y',print_warning=changed)
    endif

    if (present(zL)) then
        xyz_ptr => cfd_coupler_input%domain%z
        call init_length(xyz_ptr,zL,resize=.true.,direction='z',print_warning=changed)
    endif

	if (changed .and. cfd_code_id .ne. couette_serial) then
		print*, "Regenerate Grid with corrected sizes as above"
		call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
	endif

    ! set CFD number of cells
    if (present(nx)) then
        if (cfd_coupler_input%ncells%tag == CPL ) then
            nx = cfd_coupler_input%ncells%x
        else
            if(myid == 0) then
                write(0,*)"WARNING: nx is present in coupler_cfd_adjust_domain argument list"
                write(0,*)"         but coupler_input%ncells tag is void."
                write(0,*)"         Using CFD input value!"
            endif
        endif

    endif

    if (present(ny)) then
        if (cfd_coupler_input%ncells%tag == CPL ) then
            ny = cfd_coupler_input%ncells%y
        else
            if(myid == 0) then
                write(0,*)"WARNING: ny is present in coupler_cfd_adjust_domain argument list"
                write(0,*)"         but coupler_input%ncells tag is void."
                write(0,*)"         Using CFD input value!"
            endif
        endif
    endif
        

    if (present(nz)) then
        if (cfd_coupler_input%ncells%tag == CPL ) then
            nz = cfd_coupler_input%ncells%z
        else
            if(myid == 0) then
                write(0,*)"WARNING: nz is present in coupler_cfd_adjust_domain argument list"
                write(0,*)"         but coupler_input%ncells tag is void."
                write(0,*)"         Using CFD input value!"
            endif
        endif
    endif

    ! check id CFD cell sizes are larger than 2*sigma 
    call test_cfd_cell_sizes

contains

!-----------------------------------------------------------------------------

subroutine init_length(rin,rout,resize,direction,print_warning)
	implicit none
            
    real(kind=kind(0.d0)), intent(in)    :: rin
    real(kind=kind(0.d0)), intent(inout) :: rout
    logical, intent(in)                  :: resize
	character(*),intent(in) 			 :: direction
    logical,intent(out)					 :: print_warning

    real(kind=kind(0.d0)) :: rinit  ! store the initial value of rout or rin needed for print


    print_warning=.false.

    select case (cfd_coupler_input%domain%tag)
    case (CPL)
        select case (cfd_coupler_input%domain%units ) 
        case ("CELLSIDE","CellSide","Cellside","cellside")
            rout = b * rin
        case("SIGMA", "Sigma", "sigma")
            select case (cfd_coupler_input%domain%cell_type)
            case("FCC","Fcc","fcc")
                if (resize) then 
                rinit = rin 
                rout = real(floor(rin/b),kind(0.d0))*b
                print_warning = .true.
                endif
            case default
                write(*,*) "wrong unit cell type in coupler_cfd_adjust_domain. Stopping ... "
                ierror = COUPLER_ERROR_INIT
                call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
            end select
        end select
    case (CFD) 
        if(resize) then
            rinit = rout
            rout = real(nint(rout/b),kind(0.d0))*b
            print_warning = .true.
        endif
    case default
        write(*,*) "Wrong domain tag in coupler_cfd_adjust_domain. Stopping ... "
        ierror = COUPLER_ERROR_INIT
        call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
    end select

    if (print_warning) then 
        if (myid == 0) then 
            write(*,'(3(a,/),3a,/,2(a,f20.10),/a,/,a)') &
                    "*********************************************************************",	&
                    "WARNING - this is a coupled run which resets CFD domain size         ",	&
                    " to an integer number of MD initial cells:		                      ", 	&
                    "	Domain resized in the ", direction, " direction			          ",	&
                    " inital size =", rinit, " resized ", rout,									&
                    "								                                      ",	& 
                    "*********************************************************************"   
        endif
		!If resize is insignificant then return flag print_warning as false
		if (abs(rinit-rout) .lt. 0.01) print_warning = .false.
    end if
    
end subroutine init_length

!-----------------------------------------------------------------------------


subroutine test_cfd_cell_sizes
    implicit none

    integer, pointer :: ndim => null()

    ndim => cfd_coupler_input%domain%ndim

    if ( myid == 0) then
        if ((present(xL) .and. present(nx)) .or. &
            (cfd_coupler_input%domain%tag == CPL .and. &
            cfd_coupler_input%ncells%tag == CPL)) then
            if (xL/nx < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in x direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          xL=",xL,"nx=",nx
            endif
        endif

        if ((present(yL) .and. present(ny)) .or. &
            (cfd_coupler_input%domain%tag == CPL .and. & 
            cfd_coupler_input%ncells%tag == CPL .and. ndim > 1)) then
            if (yL/ny < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in y direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          yL=",yL,"nx=",ny
            endif
        endif

        if ((present(zL) .and. present(nz)) .or. &
            (cfd_coupler_input%domain%tag == CPL .and. &
            cfd_coupler_input%ncells%tag == CPL .and. ndim > 2 )) then
            if (zL/nz < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in z direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          zL=",zL,"nx=",nz
            endif
        endif
    end if
        
end subroutine test_cfd_cell_sizes
            
end subroutine coupler_cfd_adjust_domain


!=============================================================================
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!
!							SIMULATION
!
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!=============================================================================


!=============================================================================
! coupler_send_data wrapper for 3d arrays
! see coupler_send_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_send_data_3d(temp,index_transpose,asend_lbound,&
    asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)
    implicit none

    integer , optional, intent(in) :: glower(3), gupper(3) 
    integer, optional, intent(in) :: asend_lbound(3)       
    integer, optional, intent(in) :: asend_grid_start(3)   
    integer, optional, intent(in) :: asend_grid_end(3)     
    integer, optional, intent(in) :: index_transpose(3)    
    logical, optional, intent(in) :: use_overlap_box(3)
    real(kind=kind(0.d0)),dimension(:,:,:), intent(in) :: temp
    
    integer n1,n2,n3,n4
    real(kind=kind(0.d0)),dimension(:,:,:,:),allocatable :: asend

	!print*,'3DTemp size', size(temp,1),size(temp,2),size(temp,3)
    
    n1 = 1
    n2 = size(temp,1)
    n3 = size(temp,2)
    n4 = size(temp,3)

	!Add padding column to 3D array to make it 4D
	allocate(asend(n1,n2,n3,n4))
	asend(1,:,:,:) = temp(:,:,:)
  
    call coupler_send_data_xd(asend,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)

end subroutine coupler_send_data_3d

!=============================================================================
! coupler_send_data wrapper for 4d arrays
! see coupler_send_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_send_data_4d(asend,index_transpose,asend_lbound,asend_grid_start,&
    asend_grid_end,glower,gupper,use_overlap_box)
    implicit none

    integer , optional, intent(in) :: glower(3), gupper(3)
    integer, optional, intent(in) :: asend_lbound(3)      
    integer, optional, intent(in) :: asend_grid_start(3)  
    integer, optional, intent(in) :: asend_grid_end(3)    
    integer, optional, intent(in) :: index_transpose(3)   
    logical, optional, intent(in) :: use_overlap_box(3)
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in) :: asend
    
    integer n1,n2,n3,n4

	!print*,'4DTemp size', size(asend,1),size(asend,2),size(asend,3),size(asend,4)

    !n1 = size(asend,1)
    !n2 = size(asend,2)
    !n3 = size(asend,3)
    !n4 = size(asend,4)
   
    call coupler_send_data_xd(asend,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)

end subroutine coupler_send_data_4d


!=============================================================================
! Send data from the local grid to the associated ranks from the other 
! realm
! The entities of this operation are as follow:
!  1) the local box of the global grid held by this rank
!  2) the data associate with a subset of local box which are held in some region of asend
!  3) the associated ranks from the other realm which overlap the same region of
!     the global grid (described in map)
!
!  In index range we have as follow:
!  the input array has the limits
!  a(n,ias:iae,jas:jae,kas:kae)
!  contains local grid data in the sector
!  a(n,iags:iage,jags:jage,kags:kage)
! 
!  which corresponds to the global index region
!  igs:ige jgs:jge, kgs:kge
!
!  which must be contain in the local grid box
!  ibs:ibs, jbs:jbe, kb:kbe 
!
!  so far the transfers done correspond to ias=iags=igs=ibs, etc 
!
! Main steps of the algorithm:
!  1) Using the arguments one must find the range of the grid data in asend
!     that is exclude halos, etc
!     This can be done by passing array sector of asend 
!     e.g. asend(igrid_start:igrid_end,jgrid_start_j:grid_end,kgrid_start:kgrid_end)
!     So far only this method was tested
!  2) Use the coupler map to find which sector of the grid go to which rank from the other
!     realm
!  3) send a quick description of the data ( an array of 8 integers )
!  4) send the data if there is anything to send
!-----------------------------------------------------------------------------

subroutine coupler_send_data_xd(asend,index_transpose,asend_lbound,&
    asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)
	use mpi
	use coupler_internal_cfd, only :  CFD_COMM_OVERLAP, &
        bbox_cfd, jmax_overlap, dx, dz, icoord, nlz
    use coupler_internal_md, only : bbox_md => bbox
	use coupler_internal_common
	implicit none
 
    ! asend sizes
    !integer, intent(in) :: n1,n2,n3,n4
    ! array containing data distributed on the grid
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in) :: asend
    ! lower and upper corners of the grid to which the data of asend
    ! must be associated
    ! if not present the local grid box limits are used
    ! any negative value is discared for the default
    integer , optional, intent(in) :: glower(3), gupper(3)
    ! starting index of asend if is not (1,1,1) 
    integer, optional, intent(in) :: asend_lbound(3)
    ! starting index coordinated of grid data if is not (1,1,1)      
    integer, optional, intent(in) :: asend_grid_start(3)   
    ! upper corner of grid data stored in asend
    ! if not present the the local grid box limits are  used
    integer, optional, intent(in) :: asend_grid_end(3)
    ! specifyle the order of dimensions in asend
    ! default is ( x, y, z) but e.g. cartesian velocity
    ! in CFD solver uses (z,x,y)     
    integer, optional, intent(in) :: index_transpose(3)    ! rule of transposition for the coordinates 
    ! whether to use on MD side the extended grid of the map ( which
    ! has a halo of CFD cells around MD local box)
    ! this is useful to trap molecules that left 
    ! the local domain or if the CFD cells sixe is not
    ! conmensurate with sizes of MD local box
    logical, optional, intent(in) :: use_overlap_box(3)
    

    ! local indices 
    integer is,ie,js,je,ks,ke,bis,bie,bjs,bje,bks,bke,ix,iy,iz,ais,aie,ajs,aje,aks,ake,&
        mis,mie,mjs,mje,mks,mke,at(2,3),ig(2,3)

    ! auxiliaries 
    integer i,ndata,itag,myid, dest, req(map%n), vel_indx(8,map%n), &
        start_address(map%n+1), a_components, ierr

	real(kind=kind(0.d0)), allocatable :: vbuf(:)

	integer, save :: ncalls = 0

	! This local CFD domain is outside MD overlap zone 
	if ( map%n .eq. 0 ) return 

	call mpi_comm_rank(coupler_grid_comm,myid,ierr)

	print('(a,2i5,11f10.4)'),'sent array',myid,1,asend(1,:,1,1)
	print('(a,2i5,11f10.4)'),'sent array',myid,2,asend(1,:,2,1)
	print('(a,2i5,11f10.4)'),'sent array',myid,63,asend(1,:,63,1)
	print('(a,2i5,11f10.4)'),'sent array',myid,64,asend(1,:,64,1)

	ncalls = ncalls + 1

    ! local grid box ranges seen by this rank
    if (COUPLER_REALM .eq. COUPLER_CFD) then 
        bis = bbox_cfd%xbb(1,icoord(1,myid+1))
        bie = bbox_cfd%xbb(2,icoord(1,myid+1))
        bjs = bbox_cfd%ybb(1,icoord(2,myid+1))
        bje = bbox_cfd%ybb(2,icoord(2,myid+1))
        bks = bbox_cfd%zbb(1,icoord(3,myid+1))
        bke = bbox_cfd%zbb(2,icoord(3,myid+1))
    else
        ! using by default the indices of CFD grid inside the MD domain
        bis = bbox_md%is
        bie = bbox_md%ie
        bjs = bbox_md%js
        bje = bbox_md%je
        bks = bbox_md%ks
        bke = bbox_md%ke
        if(present(use_overlap_box)) then 
            if (use_overlap_box(1)) then
                bis = bbox_md%iso
                bie = bbox_md%ieo
            endif
            if(use_overlap_box(2)) then
                bjs = bbox_md%jso
                bje = bbox_md%jeo
            endif
            if(use_overlap_box(3)) then
                bks = bbox_md%kso
                bke = bbox_md%keo
            endif
       endif
    endif

    ! get the index order of x,y,z coordinates (handles transpose arrays)
    ix = 1 ; iy = 2; iz = 3
    if( present(index_transpose)) then
        ix = index_transpose(1)
        iy = index_transpose(2)
        iz = index_transpose(3)
    endif


    ! default start/end points for grid data in asend
    ! +1 shift in size is because asend grid indices start from 2
    is = bis
    ie = min(is + size(asend,ix+1) - 1,bie)
    js = bjs
    je = min(js + size(asend,iy+1) - 1,bje)
    ks = bks
    ke = min(ks + size(asend,iz+1) - 1,bke)

    ! Warning if asend goes over local grid box
	if (ie-is .ne. size(asend,ix+1)) & 
		write(0,'(3(a,i8),a)') "Proc=",myid, " Sent datasize ",size(asend,ix+1)," not equal to gridsize ",ie-is," in x" 
	if (je-js .ne. size(asend,iy+1)) & 
		write(0,'(3(a,i8),a)') "Proc=",myid, " Sent datasize ",size(asend,iy+1)," not equal to gridsize ",je-js," in y" 
	if (ke-ks .ne. size(asend,iz+1)) & 
		write(0,'(3(a,i8),a)') "Proc=",myid, " Sent datasize ",size(asend,iz+1)," not equal to gridsize ",ke-ks," in z"  

    ! grid data in asend data can be mapped to a specific region
    ! of the local grid 
    ! negative values are discarded in favor of defaults
    if (present(glower))then 
        if (glower(1) > 0) is = glower(1)
        if (glower(2) > 0) js = glower(2)
        if (glower(3) > 0) ks = glower(3)
    endif

    ! sanity check is needed here

    !  upper limit exteed for grid data point stored in asend
    if (present(gupper))then 
        if (gupper(1) > 0) ie = gupper(1)
        if (gupper(2) > 0) je = gupper(2)
        if (gupper(3) > 0) ke = gupper(3)
    endif

    ! sanity check is needed here

    ! Array indices in asend for the data mapped on grid
    ! +1 shift is because asend grid indices start from 2
    ais = 1
    ajs = 1
    aks = 1
    if (present(asend_lbound)) then 
        if (.not. present(asend_grid_start))then 
            write(0,*) "because the asend lower bound is not default asend_grid_start argument must be provided"
            call MPI_Abort(MPI_COMM_WORLD,COUPLER_ABORT_SEND_CFD,ierr)
        endif
        ais = asend_grid_start(ix) - asend_lbound(ix) + 1
        ajs = asend_grid_start(iy) - asend_lbound(iy) + 1
        aks = asend_grid_start(iz) - asend_lbound(iz) + 1
    endif

    aie = ais + (ie-is) 
    aje = ajs + (je-js)
    ake = aks + (ke-ks)

    ! sanity checks are needed 
    
    ! number of components at each grid point
    a_components = size(asend,dim=1)

    ! loop through the maps and send the corresponding sections of asend
    do i = 1, map%n
        dest = map%rank_list(i)

        ! indices for the grid sector to be sent to this dest
        mis = max(is, map%domains(1,i)) 
        mie = min(ie, map%domains(2,i))
        mjs = max(js, map%domains(3,i))
        mje = min(je, map%domains(4,i))
        mks = max(ks, map%domains(5,i)) 
        mke = min(ke, map%domains(6,i))	
        ndata = a_components * (mie - mis + 1) * (mje - mjs + 1) * (mke - mks + 1)

        write(0,'(a,13i8)') 'coupler indices', 	myid, is ,ie ,js ,je ,ks ,ke,map%domains(:,i)
        write(0,'(a,9i8)') 	'coupler send', 	myid,mis,mie,mjs,mje,mks,mke,a_components,ndata

        if ( ndata > 0) then 
            if(allocated(vbuf)) deallocate(vbuf)
            
            allocate(vbuf(ndata))
            
            ! Location in asend of the domain to be sent
            ig(1,ix) = mis - is + ais 
            ig(2,ix) = mie - is + ais
            ig(1,iy) = mjs - js + ajs
            ig(2,iy) = mje - js + ajs
            ig(1,iz) = mks - ks + aks
            ig(2,iz) = mke - ks + aks
 
            write(0,'(a,14i5)') ' coupler send a*, ig ...', ncalls, myid,ais,aie,ajs,aje,aks,ake,ig 

            vbuf(1:ndata) = reshape(asend(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)), (/ ndata /) )
        endif
        ! Attention ncall could go over max tag value for long runs!!
        itag = mod( ncalls, MPI_TAG_UB)
        ! send the info about the data to come
        call MPI_send((/ndata,a_components,mis,mie,mjs,mje,mks,mke/),8,MPI_INTEGER,&
            			dest, itag, COUPLER_ICOMM, ierr)
        ! send data only if there is anything to send
        if (ndata > 0) then 
			!vbuf = 1.234567891011121314151617d0
			!print'(2a,2i8,4f25.16)', 'ICP send data',code_name(COUPLER_REALM), myid, & 
			!							 size(vbuf), maxval(vbuf),minval(vbuf),sum(vbuf),vbuf(10)
            call MPI_send(vbuf, ndata, MPI_DOUBLE_PRECISION, dest, itag, COUPLER_ICOMM, ierr)
        endif
    enddo

end subroutine coupler_send_data_xd

!=============================================================================
! coupler_recv_data wrapper for 3d arrays
! see coupler_recv_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_recv_data_3d(temp,index_transpose,asend_lbound,asend_grid_start,&
    asend_grid_end,glower,gupper,accumulate,pbc)
    implicit none

    integer, optional, intent(in) :: asend_lbound(3)
    integer, optional, intent(in) :: asend_grid_start(3) 
    integer, optional, intent(in) :: asend_grid_end(3)   
    integer, optional, intent(in) :: index_transpose(3)  
    integer, optional, intent(in) :: glower(3), gupper(3)
    logical, optional, intent(in) :: accumulate          
    integer, optional, intent(in) :: pbc  
    real(kind(0.d0)),dimension(:,:,:),intent(inout) :: temp 
                                                        
	integer n1, n2, n3, n4
	real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: arecv    

    n1 = 1 
    n2 = size(temp,1)
    n3 = size(temp,2)
    n4 = size(temp,3)

	!Add padding column to 3D array to make it 4D
	allocate(arecv(n1,n2,n3,n4))
    call coupler_recv_data_xd(arecv,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,accumulate,pbc)
	temp(:,:,:) = 	arecv(1,:,:,:) 

end subroutine coupler_recv_data_3d

!=============================================================================
! coupler_recv_data wrapper for 4d arrays
! see coupler_recv_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_recv_data_4d(arecv,index_transpose,asend_lbound,&
    asend_grid_start,asend_grid_end,glower,gupper,accumulate,pbc)

    implicit none

    integer, optional, intent(in) :: asend_lbound(3)     
    integer, optional, intent(in) :: asend_grid_start(3) 
    integer, optional, intent(in) :: asend_grid_end(3)   
    integer, optional, intent(in) :: index_transpose(3)  
    integer, optional, intent(in) :: glower(3), gupper(3)
    logical, optional, intent(in) :: accumulate          
    integer, optional, intent(in) :: pbc  
    real(kind(0.d0)),dimension(:,:,:,:),intent(inout) :: arecv
                                                         

    !integer n1, n2, n3, n4

    !n1 = size(arecv,1)
    !n2 = size(arecv,2)
    !n3 = size(arecv,3)
    !n4 = size(arecv,4)

    call coupler_recv_data_xd(arecv,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,accumulate,pbc)

end subroutine coupler_recv_data_4d

!=============================================================================
! Receive data from to local grid from the associated ranks from the other 
! realm
! The entities of this operation are as follow:
!  1) the local box of the global grid held by this rank
!  2) the grid data region to be held in some region of arecv
!  3) the associated ranks from the other realm which overlap the same region of
!     the global grid
!
! Main steps of the algorithm:
!  1) Using the arguments find the range of the grid data in arecv
!     that is: exclude halos, etc
!     This can be done by passing array sector of arecv that is meant to hold grid data 
!     e.g. arecv(igrid_start:igrid_end,jgrid_start_j:grid_end,kgrid_start:kgrid_end)
!     So far only this method was tested
!  2) Use the coupler map to find where the arriving data should land in a temporary array
!  3) copy the whole lot, after the receives completes, to the corresponding location
!     of arecv
!
!  Notes: 
!  1) This subroutine can perform accumulation of data, if one chooses to 
!     send sum of cell quantities and respective numbers. Could be useful to 
!     collect contribution from molecule leaving local MD box.
!  2) Periodic boundary condition can be applied, useful for averages with
!     staggered boxes in which the end boxes collect only form half box
!     at the moment this works only with (z,x,y) array in z and x directions
!   
!    
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine coupler_recv_data_xd(arecv,index_transpose,a_lbound,&
    									a_grid_start,a_grid_end,glower,gupper,accumulate,pbc)
    use mpi
	use coupler_internal_cfd, only : bbox_cfd, jmax_overlap, dx, dz, icoord, nlz
    use coupler_internal_md, only :  bbox_md => bbox
	use coupler_internal_common
    implicit none

    ! arecv sizes
    !integer, intent(in) :: n1, n2, n3, n4
    ! array that recieves grid distributed data 
    real(kind(0.d0)), dimension(:,:,:,:),intent(inout) :: arecv
    ! lower and upper corners of the grid to which the data of asend
    ! must be associated
    ! if not present the local grid box limits are used
    ! any negative value is discared for the default
    integer, optional, intent(in) :: glower(3), gupper(3)
    ! starting index of asend if is not (1,1,1)
    integer, optional, intent(in) :: a_lbound(3) 
    ! starting index coordinated of grid data if is not (1,1,1)    
    integer, optional, intent(in) :: a_grid_start(3)
    ! upper corner of grid data stored in asend
    ! if not present the the local grid box limits are used
    integer, optional, intent(in) :: a_grid_end(3)  
    ! specifyle the order of dimensions in arecv
    ! default is ( x, y, z) but e.g. cartesian velocity
    ! in CFD solver uses (z,x,y)
    integer, optional, intent(in) :: index_transpose(3) 
    ! if present add the receiving value if they land in the same
    ! grid position
    logical, optional, intent(in) :: accumulate   
    ! if present applies periodic boundary condition
    ! posible values 1 for x direction, 3 for z direction
    integer, optional, intent(in) :: pbc                 
                                                         
    ! local indices 
    integer is,ie,js,je,ks,ke,bis,bie,bjs,bje,bks,bke,ix,iy,iz,ais,aie,ajs,aje,aks,ake,&
        at(2,3),ig(2,3)
    ! auxiliaries 
    integer i,j,k,ii,jj,kk,ndata,itag,myid, source,req(map%n), vel_indx(8,map%n), &
        p1s,p1e,p2s,p2e,p3s,p3e,pt1s,pt1e,pt2s,pt2e,pt3s,pt3e, bgt(2,3),& 
        start_address(map%n+1), status(MPI_STATUS_SIZE,map%n),ierr
    real(kind(0.d0)), allocatable ::  vbuf(:), atmp(:,:,:,:)
    integer, save :: ncalls = 0

	! This local CFD domain is outside MD overlap zone 
	if ( map%n .eq. 0 ) return 
 
	ncalls = ncalls + 1
	call MPI_comm_rank(coupler_grid_comm,myid,ierr)

    ! Local grid box
    if (COUPLER_REALM == COUPLER_CFD) then 
         bis = bbox_cfd%xbb(1,icoord(1,myid+1))
         bie = bbox_cfd%xbb(2,icoord(1,myid+1))
         bjs = bbox_cfd%ybb(1,icoord(2,myid+1))
         bje = bbox_cfd%ybb(2,icoord(2,myid+1))
         bks = bbox_cfd%zbb(1,icoord(3,myid+1))
         bke = bbox_cfd%zbb(2,icoord(3,myid+1))
     else 
		! Use default indices of CFD grid inside the MD domain
        bis = bbox_md%is
        bie = bbox_md%ie
        bjs = bbox_md%js
        bje = bbox_md%je
        bks = bbox_md%ks
        bke = bbox_md%ke
    endif

    ! Get the indices in x,y,z direction from transpose array
    ix = 1 ; iy = 2; iz = 3
    if( present(index_transpose)) then
        ix = index_transpose(1)
        iy = index_transpose(2)
        iz = index_transpose(3)
    endif

    ! default start/end points on grid for data in asend
	is = bis
    ie = bie
	js = bjs
    je = bje
	ks = bks
    ke = bke 

    ! grid data in asend data can be mapped to a specific region
    ! of the local grid -- negative values are discarded in favor of defaults
    if (present(glower))then 
        if (glower(1) > 0) is = glower(1)
        if (glower(2) > 0) js = glower(2)
        if (glower(3) > 0) ks = glower(3)
    endif

    ! sanity check is needed here

    !  upper limit exteed for grid data point stored in asend
    if (present(gupper))then 
        if (gupper(1) > 0) ie = gupper(1)
        if (gupper(2) > 0) je = gupper(2)
        if (gupper(3) > 0) ke = gupper(3)
    endif

    ! sanity check is needed here

    ! put the mapped block limits in a transposed array
    bgt(1,ix) = is
    bgt(2,ix) = ie
    bgt(1,iy) = js
    bgt(2,iy) = je
    bgt(1,iz) = ks
    bgt(2,iz) = ke

    ! Array indices in arecv for the data mapped on grid
    ! +1 shift is because asend grid indices start from 2
    ais = 1
    ajs = 1
    aks = 1
    if (present(a_lbound)) then 
        if (.not. present(a_grid_start))then 
            write(0,*) "because the arecv lower bound is not default as_grid_start argument must be provided"
            call MPI_Abort(MPI_COMM_WORLD,COUPLER_ABORT_SEND_CFD,ierr)
        endif
        ais = a_grid_start(1) - a_lbound(1) + 1
        ajs = a_grid_start(2) - a_lbound(2) + 1
        aks = a_grid_start(3) - a_lbound(3) + 1
    endif

    aie = min(ais + (ie-is),size(arecv,ix+1)) 
    aje = min(ajs + (je-js),size(arecv,iy+1))
    ake = min(aks + (ke-ks),size(arecv,iz+1))

    ! sanity checks are needed 
    
    ! Store the transposition for grid boundaries limits
    at(1,ix) = ais
    at(2,ix) = aie
    at(1,iy) = ajs
    at(2,iy) = aje
    at(1,iz) = aks
    at(2,iz) = ake

    ! First get info about what's comming
    do i = 1, map%n
        source =  map%rank_list(i)
        itag = mod( ncalls, MPI_TAG_UB)
        call MPI_irecv(vel_indx(1,i),8,MPI_Integer,source,itag,&
            			COUPLER_ICOMM,req(i),ierr)
    enddo

    call MPI_waitall(map%n, req, status, ierr)

    write(0,'(a,16i5)') 'coupler recv vel_indx ', vel_indx

    ndata = 0
    do i = 1, map%n
        ndata = ndata + max(vel_indx(1,i),0)
    enddo

	allocate(vbuf(ndata),stat=ierr)
    start_address(1) = 1 
    do i = 1, map%n
        source =  map%rank_list(i) 
        ndata = vel_indx(1,i)
        start_address(i+1) = start_address(i) + ndata
        if (ndata > 0) then 
            ! Attention ncall could go over max tag value for long runs!!
			itag = mod(ncalls, MPI_TAG_UB)
            call MPI_irecv(vbuf(start_address(i)),ndata,MPI_DOUBLE_PRECISION,source,itag,&
                				COUPLER_ICOMM,req(i),ierr)
        else 
            req(i) = MPI_REQUEST_NULL
        endif
    enddo

	!print*, 'BEFORE WAIT STATEMENT'
    call MPI_waitall(map%n, req, status, ierr)
	!print'(2a,2i8,4f25.16)', 'ICP recv data',code_name(COUPLER_REALM),myid, & 
	!							 size(vbuf), maxval(vbuf),minval(vbuf),sum(vbuf),vbuf(10)
    !	write(0,*) 'MD getCFD vel wait',  myid, id, i, source, ierr

    ! Allocate atmp corresponding to reunion of mapped grid
    ! This must be thought as landing area for data coming from the other realm
    ig(1,ix) = minval(vel_indx(3,:))
    ig(2,ix) = maxval(vel_indx(4,:))
    ig(1,iy) = minval(vel_indx(5,:))
    ig(2,iy) = maxval(vel_indx(6,:))
    ig(1,iz) = minval(vel_indx(7,:))
    ig(2,iz) = maxval(vel_indx(8,:))

    !write(0,*) 'coupler recv atmp sizes ', vel_indx(2,1),ig

    allocate(atmp(vel_indx(2,1), ig(1,1):ig(2,1), ig(1,2):ig(2,2), ig(1,3):ig(2,3)))
    atmp(:,:,:,:) = 0.d0

    do i=1,map%n
        
        if(vel_indx(1,i) <= 0) cycle
        
        ! Get array indices from the recieved grid data
        ig(1,ix) = vel_indx(3,i)
        ig(2,ix) = vel_indx(4,i)
        ig(1,iy) = vel_indx(5,i)
        ig(2,iy) = vel_indx(6,i)
        ig(1,iz) = vel_indx(7,i)
        ig(2,iz) = vel_indx(8,i)
        
        if (present(accumulate)) then
            !write(0,*) 'coupler recv ig', ig
            atmp(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)) = &
                atmp(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)) &
                + reshape(vbuf(start_address(i):start_address(i+1)-1), &
                (/ vel_indx(2,i), ig(2,1)-ig(1,1)+1, ig(2,2)-ig(1,2)+1, ig(2,3)-ig(1,3)+1 /))
        else
            atmp(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)) = &
                reshape(vbuf(start_address(i):start_address(i+1)-1), &
                (/ vel_indx(2,i), ig(2,1)-ig(1,1)+1, ig(2,2)-ig(1,2)+1, ig(2,3)-ig(1,3)+1 /))
        endif

        ! call MPI_Error_string(status(MPI_ERROR,i), err_string, len(err_string), ierr);
	    !  write(0,*) 'MD getCFD vel err, myid, i ', myid, i, trim(err_string) 
        !call MPI_get_count(status(1,i),MPI_double_precision,ib,ierr)
		! write(0,*) 'MD recv ', myid, id, i, ib, ' DP'
    enddo

    ! Transfer data from the landing area to arecv 
    ! if the flag accumulate is present normalise
    ! data with the value stored on the last line of the first column
    ! and apply periodic boundary conditions

    ! Indices for overlapping window
    ! usually the landing area and arecv coincide
    ! some test and warning should be done here
    p1s = max(lbound(atmp,2) - bgt(1,1) + at(1,1),at(1,1))
    p1e = min(ubound(atmp,2) - bgt(1,1) + at(1,1),at(2,1))
    p2s = max(lbound(atmp,3) - bgt(1,2) + at(1,2),at(1,2))
    p2e = min(ubound(atmp,3) - bgt(1,2) + at(1,2),at(2,2))
    p3s = max(lbound(atmp,4) - bgt(1,3) + at(1,3),at(1,3))
    p3e = min(ubound(atmp,4) - bgt(1,3) + at(1,3),at(2,3))

    ! needs revision, I'm not sure that it works in general !!! ARE YOU FUCKING KIDDING ME!
    pt1s = p1s + bgt(1,1) - at(1,1)
    pt1e = p1e + bgt(1,1) - at(1,1)
    pt2s = p2s + bgt(1,2) - at(1,2)
    pt2e = p2e + bgt(1,2) - at(1,2)
    pt3s = p3s + bgt(1,3) - at(1,3)
    pt3e = p3e + bgt(1,3) - at(1,3)

    !write(0,*)' p1s ...', myid, p1s,p1e,p2s,p2e,p3s,p3e, pt1s, pt1e, pt2s,pt2e,pt3s,pt3e

    if(present(accumulate)) then
        if ( present(pbc)) then
			call halos(pbc) 
            !call set_pbc(pbc)
        endif
        
        ndata = ubound(atmp,1)
        do k = p3s , p3e
 		do j = p2s , p2e
 		do i = p1s , p1e
			ii = i + bgt(1,1) - at(1,1)
			jj = j + bgt(1,2) - at(1,2)
			kk = k + bgt(1,3) - at(1,3)
 			if (atmp(ndata, ii,jj,kk) > 0.d0) then
				arecv(:,i,j,k) = atmp(1:ndata-1,ii,jj,kk)/atmp(ndata,ii,jj,kk)
			else
				arecv(:,i,j,k) = 0.d0
			endif
		end do
		end do
        end do
    else 
        arecv(:,p1s:p1e,p2s:p2e,p3s:p3e) = atmp(:,pt1s:pt1e,pt2s:pt2e,pt3s:pt3e) 
    endif


	print*,'Size of recieved data',ig(1,1),ig(2,1),ig(1,2),ig(2,2),ig(1,3),ig(2,3),map%n, & 
									start_address(i),start_address(i+1)-1, &
					                vel_indx(2,i),ig(2,1)-ig(1,1)+1,ig(2,2)-ig(1,2)+1,ig(2,3)-ig(1,3)+1, &
									p1s,p1e,p2s,p2e,p3s,p3e,pt1s,pt1e,pt2s,pt2e,pt3s,pt3e

	!print'(2a,2i8,4f25.16)', 'ICP trfr data',code_name(COUPLER_REALM),myid, & 
	!						size(arecv), maxval(arecv),minval(arecv),sum(arecv),arecv(1,p1s,p2s,p3s)
            
contains

!-----------------------------------------------------------------------------
! Sets periodic boundary conditions
! adds together the sums of the first and last bins on the global grid (x and z directions)
! and puts sum in both boxes
! Note:
! attention coord array are shifted with 1 in CFD and MD
! to get correct result in MPI function such as mpi_cart_rank, subtract 1 
!-----------------------------------------------------------------------------
subroutine set_pbc(pbc)
    use coupler_internal_cfd, only : npx
    implicit none
    integer, intent(in) :: pbc

    real(kind(0.d0)), allocatable, dimension(:,:,:) :: x1, x2
    integer dest, ip(3), ierr, status(MPI_STATUS_SIZE)

    ! PBC works only for 3,2,1 coordinate transposition

    if ( .not. (ix == 2 .and. iy == 3 .and. iz == 1)) then
        write(0,*) " Periodic boundary condition not implemented for this layout of data"
        return
    endif

    ! Sanity check atmp must touch the boundary of the global grid in zin x directions
    ! I guess that the PBC are apllied at the ends of global grid, check with CFD people
    select case(pbc)
    case(3)
        ! first array dimension which corresponds to z direction for velocity arrays 
        allocate(x1(size(atmp,1),size(atmp,3),size(atmp,4)),x2(size(atmp,1),size(atmp,3),size(atmp,4)))
        x1 =  atmp(:,aks, :,:)
        x2 =  atmp(:,ake, :,:)
        atmp(:, aks, :, :) =   atmp(:, aks, :, :) + x2
        atmp(:, ake, :, :) =   atmp(:, ake, :, :) + x1
    case(1)
        ! second array dimension which correponds to x direction
        allocate(x1(size(atmp,1),size(atmp,2),size(atmp,4)),x2(size(atmp,1),size(atmp,2),size(atmp,4)))
        if ( npx == 1 )then 
            ! no MPI communication needed  
            !write(0,*) 'coupler internal cfd pbc', npx
            x1 =  atmp(:, :,lbound(atmp,3), :)                 
            x2 =  atmp(:, :,ubound(atmp,3), :)
            atmp(:, :, lbound(atmp,3), :) =  atmp(:, :, lbound(atmp,3), :) + x2
            atmp(:, :, ubound(atmp,3), :) =  atmp(:, :, ubound(atmp,3), :) + x1  
        else 
            call MPI_comm_rank(coupler_grid_comm,myid,ierr)
            ip = icoord(: ,myid+1)
            if (ip(1) == 1 .and. ip(2) == 1) then 
                x1 =  atmp(:, :,lbound(atmp,3),:)
                call MPI_cart_rank(coupler_grid_comm, (/ npx-1, 0, ip(3)-1 /), dest, ierr)
                call MPI_sendrecv(x1,size(x1),MPI_DOUBLE_PRECISION,dest,1,x2,size(x2),MPI_DOUBLE_PRECISION,&
                    dest,1,coupler_grid_comm,status,ierr)
                atmp(:, :, lbound(atmp,3), :) =  atmp(:, :, lbound(atmp,3), :) + x2
            else if (ip(1) == npx .and. ip(2) == 1) then 
                x2 =  atmp(:, :,ubound(atmp,3), :)
                call MPI_cart_rank(coupler_grid_comm, (/ 0, 0, ip(3) - 1 /), dest, ierr)
                call MPI_sendrecv(x2,size(x2),MPI_DOUBLE_PRECISION,dest,1,x1,size(x1),MPI_DOUBLE_PRECISION,&
                    dest,1,coupler_grid_comm,status,ierr)
                atmp(:, :, ubound(atmp,3), :) =  atmp(:, :, ubound(atmp,3), :) + x1
            endif
        endif
    end select
end subroutine set_pbc

!-----------------------------------------------------------------------------
! Halo boundary condition swapped or combined
!-----------------------------------------------------------------------------
subroutine halos(pbc)
    use coupler_internal_cfd, only : npx
    implicit none

    integer, intent(in) :: pbc
    real(kind(0.d0)), allocatable, dimension(:,:,:) :: x1, x2
    integer dest, ip(3), ierr, status(MPI_STATUS_SIZE)

    ! PBC works only for 3,2,1 coordinate transposition
    if ( .not. (ix == 2 .and. iy == 3 .and. iz == 1)) then
        write(0,*) " Periodic boundary condition not implemented for this layout of data"
        return
    endif

    ! Sanity check atmp must touch the boundary of the global grid in zin x directions
    ! I guess that the PBC are applied at the ends of global grid, check with CFD people
    select case(pbc)
    case(3)
        ! first array dimension which corresponds to z direction for velocity arrays 
        allocate(x1(size(atmp,1),size(atmp,3),size(atmp,4)),x2(size(atmp,1),size(atmp,3),size(atmp,4)))
        x1 =  atmp(:, aks, :, :)
        x2 =  atmp(:, ake, :, :)
        atmp(:, aks, :, :) =  x2 ! atmp(:, aks, :, :) + x2
        atmp(:, ake, :, :) =  x1 ! atmp(:, ake, :, :) + x1
    case(1)
        ! second array dimension which correponds to x direction
        allocate(x1(size(atmp,1),size(atmp,2),size(atmp,4)),x2(size(atmp,1),size(atmp,2),size(atmp,4)))
        if ( npx == 1 )then 
            ! no MPI communication needed  
            !write(0,*) 'coupler internal cfd pbc', npx
            x1 =  atmp(:, :,lbound(atmp,3), :)                 
            x2 =  atmp(:, :,ubound(atmp,3), :)
            atmp(:, :, lbound(atmp,3), :) = x2  !atmp(:, :, lbound(atmp,3), :) + x2
            atmp(:, :, ubound(atmp,3), :) = x1  !atmp(:, :, ubound(atmp,3), :) + x1  
        else 
            call MPI_comm_rank(coupler_grid_comm,myid,ierr)
            ip = icoord(: ,myid+1)
            if (ip(1) == 1 .and. ip(2) == 1) then 
                x1 =  atmp(:, :,lbound(atmp,3),:)
                call MPI_cart_rank(coupler_grid_comm, (/ npx-1, 0, ip(3)-1 /), dest, ierr)
                call MPI_sendrecv(x1,size(x1),MPI_DOUBLE_PRECISION,dest,1,x2,size(x2),MPI_DOUBLE_PRECISION,&
                    				dest,1,coupler_grid_comm,status,ierr)
                atmp(:, :, lbound(atmp,3), :) = x2		!atmp(:, :, lbound(atmp,3), :) + x2
            else if (ip(1) == npx .and. ip(2) == 1) then 
                x2 =  atmp(:, :,ubound(atmp,3), :)
                call MPI_cart_rank(coupler_grid_comm, (/ 0, 0, ip(3) - 1 /), dest, ierr)
                call MPI_sendrecv(x2,size(x2),MPI_DOUBLE_PRECISION,dest,1,x1,size(x1),MPI_DOUBLE_PRECISION,&
                    				dest,1,coupler_grid_comm,status,ierr)
                atmp(:, :, ubound(atmp,3), :) =  x1		!atmp(:, :, ubound(atmp,3), :) + x1
            endif
        endif
    end select
end subroutine halos

end subroutine coupler_recv_data_xd

!============================================================================
!
! utility functions and subroutines that extract parameters from internal modules 
!
!-----------------------------------------------------------------------------    

!-------------------------------------------------------------------------------
! return to the caller coupler parameters from cfd realm
!-------------------------------------------------------------------------------
subroutine coupler_cfd_get(jmax_overlap,jmin,jmino)
    use coupler_internal_cfd, jmax_overlap_ => jmax_overlap,jmin_=>jmin,jmino_=>jmino
    implicit none

    integer,optional,intent(out) :: jmax_overlap,jmin,jmino

    if(present(jmax_overlap)) then
        jmax_overlap = jmax_overlap_
    endif

    if(present(jmin))then 
        jmin = jmin_
    end if

    if(present(jmino))then
        jmino = jmino_
    end if

end subroutine coupler_cfd_get


!===============================================================================
! Return to the caller coupler parameters from MD realm
!-------------------------------------------------------------------------------
subroutine coupler_md_get(	xL_md,yL_md,zL_md, MD_initial_cellsize, top_dy, &
							overlap_with_continuum_force, overlap_with_top_cfd, ymin_continuum_force, &
    						ymax_continuum_force, xmin_cfd_grid, xmax_cfd_grid, zmin_cfd_grid, zmax_cfd_grid, &
    						dx_cfd, dz_cfd, cfd_box)
    use mpi
	use coupler_internal_md, only : xL_md_ =>  xL_md, yL_md_ =>  yL_md, zL_md_ => zL_md, &
		b => MD_initial_cellsize, x, y, z,j => jmax_overlap_cfd, dx, dz, bbox, half_domain_lengths, &
        cfd_box_ => cfd_box
	implicit none

	real(kind(0.d0)), optional, intent(out) :: xL_md, yL_md, zL_md, MD_initial_cellsize, top_dy,&
        ymin_continuum_force,ymax_continuum_force, xmin_cfd_grid, xmax_cfd_grid, &
		 zmin_cfd_grid, zmax_cfd_grid, dx_cfd, dz_cfd 
    logical, optional, intent(out) :: overlap_with_continuum_force, overlap_with_top_cfd
    type(cfd_grid_info), optional, intent(out) :: cfd_box

    integer ierr

    if (present(xL_md)) then
		xL_md = xL_md_
	endif

	if (present(yL_md)) then 
		yL_md = yL_md_
	endif

	if (present(zL_md)) then 
		zL_md = zL_md_
	endif

	if (present(MD_initial_cellsize)) then 
		MD_initial_cellsize = b
	endif

	if (present(top_dy)) then
    	top_dy = y(j) - y(j-1)
	end if

    if(present(overlap_with_continuum_force)) then
        ! check if the countinuum force domain is contained in one 
        ! domain along y direction
        if ( (y(j - 2) < bbox%bb(1,2) .and. &
              y(j - 1) > bbox%bb(1,2)) .or. &
             (y(j - 2) < bbox%bb(2,2) .and. &
              y(j - 1) > bbox%bb(2,2))) then
            write(0,*) " the region in which the continuum constraint force is applied "
            write(0,*) " spans over two domains. This case is not programmed, please investigate"
            call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_CONTINUUM_FORCE,ierr)
        endif
        
        if ( y(j - 1) <  bbox%bb(2,2) .and. y(j - 2) >= bbox%bb(1,2) ) then
            overlap_with_continuum_force = .true.
        else   
            overlap_with_continuum_force = .false.
        endif

    endif

     if(present(overlap_with_top_cfd)) then
        ! check if the MD domain overlaps with the top of CFD grid (along y)
        ! the MD constrain force is applyied top layer of cfd cells
        if ( y(j - 1) < bbox%bb(2,2) .and. y(j - 1) >= bbox%bb(1,2) ) then
            overlap_with_top_cfd = .true.
        else   
            overlap_with_top_cfd = .false.
        endif

    endif
 
     if(present(ymin_continuum_force)) then
         ymin_continuum_force = y(j - 2) - bbox%bb(1,2) - half_domain_lengths(2)
     endif
         
     if(present(ymax_continuum_force)) then
         ymax_continuum_force = y(j - 1) - bbox%bb(1,2) - half_domain_lengths(2)
     endif 

     if(present(xmin_cfd_grid)) then
         xmin_cfd_grid = x(bbox%is) - bbox%bb(1,1) - half_domain_lengths(1)
     endif

     if(present(xmax_cfd_grid)) then
         xmax_cfd_grid = x(bbox%ie) - bbox%bb(1,1) - half_domain_lengths(1)
     endif

     if(present(zmin_cfd_grid)) then
         zmin_cfd_grid = z(bbox%ks) - bbox%bb(1,3) - half_domain_lengths(3)
     endif

     if(present(zmax_cfd_grid)) then
         zmax_cfd_grid = z(bbox%ke) - bbox%bb(1,3) - half_domain_lengths(3)
     endif
 
     if (present(dx_cfd)) then 
         dx_cfd = dx
     endif

     if (present(dz_cfd)) then 
         dz_cfd = dz
     endif

     if (present(cfd_box))then
         cfd_box%gimin = cfd_box_%gimin
         cfd_box%gimax = cfd_box_%gimax
         cfd_box%gjmin = cfd_box_%gjmin
         cfd_box%gjmax = cfd_box_%gjmax
         cfd_box%gkmin = cfd_box_%gkmin
         cfd_box%gkmax = cfd_box_%gkmax
         cfd_box%imin  = cfd_box_%imin
         cfd_box%imax  = cfd_box_%imax
         cfd_box%jmin  = cfd_box_%jmin
         cfd_box%jmax  = cfd_box_%jmax
         cfd_box%kmin  = cfd_box_%kmin
         cfd_box%kmax  = cfd_box_%kmax
         cfd_box%imino = cfd_box_%imino
         cfd_box%imaxo = cfd_box_%imaxo
         cfd_box%jmino = cfd_box_%jmino
         cfd_box%jmaxo = cfd_box_%jmaxo
         cfd_box%kmino = cfd_box_%kmino
         cfd_box%kmaxo = cfd_box_%kmaxo
         cfd_box%xmin  = cfd_box_%xmin
         cfd_box%xmax  = cfd_box_%xmax
         cfd_box%dx    = cfd_box_%dx
         cfd_box%ymin  = cfd_box_%ymin
         cfd_box%ymax  = cfd_box_%ymax
         allocate( cfd_box%y(cfd_box%jmin:cfd_box%jmax))
         cfd_box%y(cfd_box%jmin:cfd_box%jmax) = cfd_box_%y(cfd_box%jmin:cfd_box%jmax)
         cfd_box%zmin = cfd_box_%zmin
         cfd_box%zmax = cfd_box_%zmax
         cfd_box%dz = cfd_box_%dz
     end if

end subroutine coupler_md_get

!-----------------------------------------------------------------------------
 
function coupler_md_get_save_period() result(p)
	use coupler_internal_md, only : save_period
	implicit none

	integer p
	p = save_period

end function coupler_md_get_save_period

!----------------------------------------------------------------------------- 

function coupler_md_get_average_period() result(p)
	use coupler_internal_md, only : average_period
	implicit none

	integer p
	p = average_period

end function coupler_md_get_average_period

!----------------------------------------------------------------------------- 

function coupler_md_get_md_steps_per_cfd_dt() result(p)
	use coupler_internal_md, only : md_steps_per_dt_cfd
	implicit none

	integer p
	p = md_steps_per_dt_cfd

end function coupler_md_get_md_steps_per_cfd_dt

!-----------------------------------------------------------------------------

function coupler_md_get_nsteps() result(n)
	 use coupler_internal_md, only : nsteps
	 implicit none 

	 integer n
	 n = nsteps

end function coupler_md_get_nsteps

!-----------------------------------------------------------------------------

function coupler_md_get_dt_cfd() result(t)
	 use coupler_internal_md, only : dt_CFD  
	 implicit none

	 real(kind=kind(0.d0)) t
	 t = dt_CFD

end function coupler_md_get_dt_cfd

!-----------------------------------------------------------------------------  
subroutine coupler_md_set(zL_md)
	use coupler_internal_md, only : zL => zL_md, z, dz, kmino, kmin_cfd, kmax_cfd, kmaxo
	implicit none
	
	real(kind(0.d0)), optional, intent(in) :: zL_md

	if ( present(zL_md) )then
		zL = zL_md
		if (allocated(z)) then 
			deallocate(z)
		endif
		! values need for the overlap map
		allocate(z(2))
		z(:)       = (/ 0.d0, zL_md /)
		dz	 = zL_md
		kmino      = 1
		kmin_cfd   = 1
		kmax_cfd   = 2
		kmaxo      = 2
	endif

end subroutine coupler_md_set


function coupler_md_get_density() result(r)
	use coupler_internal_md, only : density
	implicit none

	real(kind(0.d0)) r
	r = density

end function coupler_md_get_density

!------------------------------------------------------------------------------
function coupler_md_get_cfd_id() result(r)
    use coupler_internal_md, only : cfd_code_id
    implicit none

    integer r
    r = cfd_code_id

end function coupler_md_get_cfd_id

end module coupler
