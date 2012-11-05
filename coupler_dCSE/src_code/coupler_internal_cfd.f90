!=============================================================================
!!				   Coupler internal CFD   				   
!! Internal data and subroutines used by the coupler when working in CFD realm
!! It must not be used by subroutimes working in MD realm
! Subroutines include:
!
! coupler_cfd_init          (cfd)    initialises coupler with CFD data
! create_map_cfd 			Establish for all CFD processors the mapping (if any) 
! 							to coupled MD processors
! find_overlaps				Establish location and size of the overlap region 
! 							between the continuum and molecular simulations
! make_bbox					Make bbox which contains all domain information 
!							required in coupling process such as domain extents
! recv_vel_MD(vel,p1s,p1e, 	Get velocity fields from MD for 
!	   p2s,p2e,p3s,p3e,pbc) the boundary condiitons needed in CFD
! set_pbc(pbc)
!
!  Lucian Anton, November 2011
!
!=============================================================================


module coupler_internal_cfd
    implicit none

contains

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


subroutine coupler_cfd_init(nsteps,dt_cfd,icomm_grid,icoord,npxyz_cfd,xyzL,ncxyz, & 
							   density,ijkcmax,ijkcmin,iTmin,iTmax,jTmin, & 
							   jTmax,kTmin,kTmax,xg,yg,zg)
	use coupler, only : CPL_rank_map
    use mpi
    use coupler_module,	dt_cfd_=>dt_cfd,						&	!Simulation lengths	
						xg_=>xg,yg_=>yg,zg_=>zg				!CFD grid arrays
 	implicit none			

    integer,					    intent(in)	:: nsteps,icomm_grid 
    integer,dimension(3),		    intent(in)	:: ijkcmin,ijkcmax,npxyz_cfd,ncxyz
	integer,dimension(:),		    intent(in)	:: iTmin,iTmax,jTmin,jTmax,kTmin,kTmax
    integer,dimension(:,:),		    intent(in)	:: icoord
    real(kind(0.d0)),			    intent(in)	:: dt_cfd,density
    real(kind(0.d0)),dimension(3),  intent(in)	:: xyzL
    real(kind(0.d0)),dimension(:  ),intent(in)	:: zg
    real(kind(0.d0)),dimension(:,:),intent(in)	:: xg,yg


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

	write(99+rank_realm,*), 'CFD side',rank_realm,'CFD procs',npx_cfd,npy_cfd,npz_cfd,nproc_cfd

	! Receive & Store MD number of processors
	allocate(buf(3))
    call MPI_bcast(   buf   ,3,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)	!Receive
    npx_md = buf(1)
    npy_md = buf(2)
    npz_md = buf(3)
    nproc_md = npx_md * npy_md * npz_md
	deallocate(buf)

	write(99+rank_realm,*), 'CFD side',rank_realm,'MD procs',npx_md ,npy_md ,npz_md ,nproc_md

	! Store & Send CFD processor rank to coord
    allocate(rank2coord_cfd(3,nproc_cfd),stat=ierr); rank2coord_cfd = icoord
	iblock_realm=icoord(1,rank_realm); jblock_realm=icoord(2,rank_realm); kblock_realm=icoord(3,rank_realm)
	allocate(buf(3*nproc_cfd)); buf = reshape(icoord, (/ 3*nproc_cfd /) )
    call MPI_bcast(buf,3*nproc_cfd,MPI_INTEGER,source,CPL_INTER_COMM,ierr)	!Send
	deallocate(buf)

	call write_matrix_int(rank2coord_cfd,'cfd side, rank2coord_cfd=',99+rank_realm)

	! Receive & Store MD processor rank to coord
	allocate(buf(3*nproc_md))
    call MPI_bcast(buf,3*nproc_md ,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)	!Receive
    allocate(rank2coord_md (3,nproc_md),stat=ierr); rank2coord_md = reshape(buf,(/ 3,nproc_md /))
	deallocate(buf)

	call write_matrix_int(rank2coord_md,'cfd side, rank2coord_md=',99+rank_realm)

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

	write(99+rank_realm,*), 'CFD side',rank_realm,'coord2rank_cfd=',coord2rank_cfd

	! Receive & Store MD coordinate to rank mapping
	allocate(buf(nproc_md))
    call MPI_bcast(buf,nproc_md ,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)	!Receive
    allocate(coord2rank_md (npx_md,npy_md,npz_md)) 
	coord2rank_md = reshape(buf,(/ npx_md,npy_md,npz_md /))
	deallocate(buf)

	write(99+rank_realm,*), 'CFD side',rank_realm,'coord2rank_md=',coord2rank_md

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

	write(99+rank_realm,*), 'CFD side',rank_realm, 'rank_cfd2rank_world',      rank_cfdrealm2rank_world
	write(99+rank_realm,*), 'CFD side',rank_realm, 'rank_world2rank_cfdrealm', rank_world2rank_cfdrealm

	! Receive & Store MD mapping from realm to local rank from MD
	allocate(rank_mdrealm2rank_world(nproc_md))
	call MPI_bcast(rank_mdrealm2rank_world,nproc_md,MPI_integer,0,CPL_INTER_COMM,ierr)	!Receive

	write(99+rank_realm,*), 'CFD side',rank_realm, 'rank_md2rank_world',      rank_mdrealm2rank_world


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

	write(99+rank_realm,*), 'CFD side',rank_cart, 'rank_cfd2rank_world',      rank_cfdcart2rank_world
	write(99+rank_realm,*), 'CFD side',rank_cart, 'rank_world2rank_cfdcart', rank_world2rank_cfdcart

	! Receive & Store MD mapping from cart to local rank from MD
	allocate(rank_mdcart2rank_world(nproc_md))
	call MPI_bcast(rank_mdcart2rank_world,nproc_md,MPI_integer,0,CPL_INTER_COMM,ierr)	!Receive

	write(99+rank_realm,*), 'CFD side',rank_cart, 'rank_md2rank_world',      rank_mdcart2rank_world


	! ------------------ Timesteps and iterations ------------------------------
	! Store & send CFD nsteps and dt_cfd
	nsteps_cfd = nsteps
    call MPI_bcast(nsteps,1,MPI_integer,source,CPL_INTER_COMM,ierr)			!Send
	dt_cfd_ = dt_cfd
    call MPI_bcast(dt_cfd,1,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	write(99+rank_realm,*), 'CFD side',rank_realm,'CFD times', nsteps_cfd,dt_cfd

	! Receive & store MD timestep dt_md
    call MPI_bcast(dt_md,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive
    call MPI_bcast(nsteps_md,1,MPI_integer,     0,CPL_INTER_COMM,ierr)		!Receive

	write(99+rank_realm,*), 'CFD side',rank_realm,'MD times', nsteps_md,dt_md

	! ------------------ Send CFD grid extents ------------------------------

	! Store & send CFD density
	density_cfd = density
	call MPI_bcast(density_cfd,1,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	write(99+rank_realm,*), 'CFD side',rank_realm,'CFD density',density_cfd

	! Receive & store MD density
	call MPI_bcast(density_md,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive

	write(99+rank_realm,*), 'CFD side',rank_realm,'CFD density',density_md

	! Store & send CFD domain size
	xL_cfd = xyzL(1); yL_cfd = xyzL(2); zL_cfd = xyzL(3)
	call MPI_bcast(xyzL,3,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	write(99+rank_realm,*), 'CFD side',rank_realm,'CFD domain',xL_cfd,yL_cfd,zL_cfd	
	! Receive & store MD domain size
	allocate(rbuf(3))
	call MPI_bcast(rbuf,3,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive
	xL_md = rbuf(1); yL_md = rbuf(2); zL_md = rbuf(3);
	deallocate(rbuf)

	write(99+rank_realm,*), 'CFD side',rank_realm,'MD domain', xL_md,yL_md,zL_md

	! Store & send CFD grid extents
	icmin = ijkcmin(1); jcmin = ijkcmin(2); kcmin = ijkcmin(3)
	icmax = ijkcmax(1); jcmax = ijkcmax(2); kcmax = ijkcmax(3)
    call MPI_bcast((/ icmin,icmax,jcmin,jcmax,kcmin,kcmax /),6,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send

	write(99+rank_realm,*), 'CFD side',rank_realm,'CFD global extents',icmin,icmax,jcmin,jcmax,kcmin,kcmax

	! Store & send global number of cells in CFD
	ncx = ncxyz(1); ncy = ncxyz(2); ncz = ncxyz(3)
    call MPI_bcast(ncxyz,3,MPI_INTEGER,source,CPL_INTER_COMM,ierr)				!Send

	write(99+rank_realm,*), 'CFD side',rank_realm,'CFD global cells',ncx,ncy,ncz

	! Store & send array of global grid points
    allocate(xg_(size(xg+1,1)+1,size(xg,2)+1),stat=ierr); xg_ = xg
    allocate(yg_(size(yg+1,1)+1,size(yg,2)+1),stat=ierr); yg_ = yg
    allocate(zg_(size(zg+1,1)+1			   ),stat=ierr);  zg_ = zg
    call MPI_bcast(xg,size(xg),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(yg,size(yg),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(zg,size(zg),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send

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

	write(99+rank_realm,*), 'CFD side',rank_realm, & 
		'CFD local cells',icPmin_cfd,icPmax_cfd,jcPmin_cfd,jcPmax_cfd,kcPmin_cfd,kcPmax_cfd

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

	write(99+rank_realm,*), 'CFD side - y overlap',ncy_olap

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

	    ! test if MD_init_cell size is larger than CFD cell size
!		if (myid_realm .eq. rootid_realm) then
!		    if( MD_initial_cellsize .ge. dx .or. & 
!				MD_initial_cellsize .ge. dy .or. & 
!				MD_initial_cellsize .ge. dz .and. rank_realm == 0 ) then
!		        write(*,*)
!		        write(*,*) "********************************************************************"
!		        write(*,*) " WARNING ...WARNING ...WARNING ...WARNING ...WARNING ...WARNING ... "
!		        write(*,*) " MD initialisation cell size larger than CFD x,y cell sizes         "
!		        write(*,*) " MD_init_cellsize=",MD_initial_cellsize
!		        write(*,'(3(a,f10.5))') " dx=",xg(2,1)-xg(1,1),  & 
!										" dy=",yg(1,2)-yg(1,1),  & 
!										" dz=",zg(2  )-zg(1  )
!		        write(*,*) "********************************************************************"
!		        write(*,*)
!		    endif
!		endif

	end subroutine check_mesh

end subroutine coupler_cfd_init


end module coupler_internal_cfd
