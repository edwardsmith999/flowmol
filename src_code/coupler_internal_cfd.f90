!=============================================================================
!				   Coupler internal CFD   				   
! Internal data and subroutines used by the coupler when working in CFD realm
! It must not be used by subroutimes working in MD realm
! Subroutines include:
!
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
    save

    ! Internal data 
    ! bounding box type that holds on each processor the boundaries of the
    ! CFD grid sectors allocated to every processor
    type bbox
        integer, allocatable :: xbb(:,:), ybb(:,:), zbb(:,:)
    end type bbox

    type(bbox),target :: bbox_cfd, bbox_md

contains


! ----------------------------------------------------------------------------
! Initialisation routine for coupler - Every variable is sent and stored
! to ensure both md and cfd region have an identical list of parameters

subroutine coupler_cfd_init(nsteps,dt_cfd,icomm_grid,icoord,npxyz_cfd,xyzL,ncxyz, & 
							   density,ijkcmax,ijkcmin,iTmin,iTmax,jTmin, & 
							   jTmax,kTmin,kTmax,xg,yg,zg)
    use mpi
    use coupler_module,	dt_cfd_=>dt_cfd,						&	!Simulation lengths	
						xg_=>xg,yg_=>yg,zg_=>zg				!CFD grid arrays
	use coupler_input_data, only : cfd_coupler_input
 	implicit none			

    integer,					    intent(in)	:: nsteps,icomm_grid 
    integer,dimension(3),		    intent(in)	:: ijkcmin,ijkcmax,npxyz_cfd,ncxyz
	integer,dimension(:),		    intent(in)	:: iTmin,iTmax,jTmin,jTmax,kTmin,kTmax
    integer,dimension(:,:),		    intent(in)	:: icoord
    real(kind(0.d0)),			    intent(in)	:: dt_cfd,density
    real(kind(0.d0)),dimension(3),  intent(in)	:: xyzL
    real(kind(0.d0)),dimension(:  ),intent(in)	:: zg
    real(kind(0.d0)),dimension(:,:),intent(in)	:: xg,yg


    integer											:: i,ib,jb,kb,pcoords(3),root,source
    integer,dimension(:),allocatable				:: buf
	real(kind=kind(0.d0))							:: dxmin,dxmax,dzmin,dzmax
    real(kind=kind(0.d0)),dimension(:),allocatable 	:: rbuf

    ! Duplicate grid communicator for coupler use
    call MPI_comm_dup(icomm_grid,CPL_CART_COMM,ierr)
    call MPI_comm_rank(CPL_CART_COMM,myid_cart,ierr) 
    rank_cart = myid_cart + 1; rootid_cart = 0
	!Send only from root processor
    if ( rank_realm .eq. root ) then
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

	!Setup CFD mapping from coordinate to rank, store and send
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
    allocate(zg_(size(zg+1,1)+1			   ),stat=ierr); zg_ = zg
    call MPI_bcast(xg,size(xg),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(yg,size(yg),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send
    call MPI_bcast(zg,size(zg),MPI_double_precision,source,CPL_INTER_COMM,ierr) !Send

	call write_matrix(xg,'cfd side, xg=',50+rank_realm)
	call write_matrix(yg,'cfd side, yg=',50+rank_realm)
	write(50+rank_realm,*), 'CFD side',rank_realm,'zg',zg

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

	write(99+rank_realm,*), 'CFD side',rank_realm,'CFD local cells',icPmin_cfd,icPmax_cfd,jcPmin_cfd,jcPmax_cfd,kcPmin_cfd,kcPmax_cfd

	!Calculate the cell sizes dx,dy & dz
	dx = xL_cfd/ncx  !xg(2,1)-xg(1,1)
	dy = yL_cfd/ncy	 !yg(1,2)-yg(1,1)
	dz = zL_cfd/ncz  !zg(2  )-zg(1  )

    ! send number of overlap cells, read by input, to all processors
    ! Note : jcmax_overlap default is provided in coupler_internal_cfd
	! but for some reason it is only broadcast to CFD processors while 
	! the tag used in the if statement below is not broadcast at all...
	if (rank_realm .eq. root) then
	    if (cfd_coupler_input%overlap%tag == CPL) then
			ncx_olap = 0
    	    ncy_olap = cfd_coupler_input%overlap%y_overlap
			!ncz_olap = 0
		else
			call error_abort("j overlap not specified in COUPLER.in")
	    endif
	endif

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

	    ! test if MD_init_cell size is larger than CFD cell size
		if (rank_realm .eq. root) then
		    if( MD_initial_cellsize .ge. dx .or. & 
				MD_initial_cellsize .ge. dy .or. & 
				MD_initial_cellsize .ge. dz .and. rank_realm == 0 ) then
		        write(*,*)
		        write(*,*) "********************************************************************"
		        write(*,*) " WARNING ...WARNING ...WARNING ...WARNING ...WARNING ...WARNING ... "
		        write(*,*) " MD initialisation cell size larger than CFD x,y cell sizes         "
		        write(*,*) " MD_init_cellsize=",MD_initial_cellsize
		        write(*,'(3(a,f10.5))') " dx=",xg(2,1)-xg(1,1),  & 
										" dy=",yg(1,2)-yg(1,1),  & 
										" dz=",zg(2  )-zg(1  )
		        write(*,*) "********************************************************************"
		        write(*,*)
		    endif
		endif

	end subroutine check_mesh

end subroutine coupler_cfd_init



!=============================================================================
! Establish for all CFD processors the mapping (if any) 
! to coupled MD processors
!-----------------------------------------------------------------------------

subroutine create_map_cfd
    use mpi 
    use coupler_module !, only : CPL_REALM_COMM, CPL_CART_COMM, CPL_INTER_COMM, cfd_is_2d, map
    implicit none

    integer 				:: i, id_coord, color, noverlaps, ir, iaux(4)
    integer, allocatable 	:: md_grid_boxes(:,:), overlap_mask(:), ireq(:), overlap_box(:,:)
    real(kind(0.d0))		:: raux(2)

	call MPI_barrier(CPL_REALM_COMM,ierr)

    ! CPL_CART_COMM has the processor coordinates in cartesian topology
    call mpi_comm_rank(CPL_CART_COMM,id_coord,ierr)
    id_coord = id_coord+1 ! get in sync with icoord convention

    call make_bbox

	! Receive the overlapping box indices from all MD processors
    allocate(md_grid_boxes(6,0:nproc_md-1), overlap_mask(0:nproc_md-1), &
        	   overlap_box(6,0:nproc_md-1),         ireq(0:nproc_md-1))
	call mpi_allgather(MPI_BOTTOM,0,MPI_INTEGER,md_grid_boxes,6,MPI_INTEGER,CPL_INTER_COMM, ierr)
 
    ! Find overlaps and send domain overlap mask to all MD processors
    call find_overlaps
    call mpi_allgather(overlap_mask,nproc_md,MPI_INTEGER,MPI_BOTTOM,0,MPI_INTEGER,CPL_INTER_COMM,ierr)

    noverlaps = 0
    do i = 0, nproc_md - 1
        if ( overlap_mask(i) == 1) then 
            noverlaps = noverlaps + 1
        endif
    enddo

    map%n = noverlaps
    allocate ( map%rank_list(noverlaps), map%domains(6,noverlaps))

    !  Overlapping communicator
    if ( map%n > 0) then
        color = 1
    else 
        color = 0
    endif

    call mpi_comm_split(CPL_REALM_COMM, color, rank_realm, CPL_REALM_INTERSECTION_COMM, ierr)

    if (color == 0) then
        CPL_REALM_INTERSECTION_COMM = MPI_COMM_NULL
    endif

    ! Send the range of the overlaping domains
    ir = 0
    do i=0, nproc_md-1
        if (overlap_mask(i) == 1) then
            call mpi_isend(overlap_box(1,i),6,mpi_integer,i,2,CPL_INTER_COMM,ireq(i),ierr)
            ir = ir + 1
            map%rank_list(ir) = i
            map%domains(:,ir) = overlap_box(1:6,i)
        else
            ireq(i) = MPI_REQUEST_NULL
        endif
    enddo
    call mpi_waitall(nproc_md,ireq,MPI_STATUSES_IGNORE,ierr) 

contains

    !-----------------------------------------------------------------------------
    ! Make bbox which contains all domain information required in coupling process
    ! such as domain extents
    !-----------------------------------------------------------------------------

    subroutine make_bbox
            implicit none

        integer, parameter :: is = 1, ie = 2

        integer ixyz, i, nixyz(3), minxyz(3), npxyz(3)
        integer, pointer :: bb_ptr(:,:) => null()

		print*, icmax,icmin,npx_cfd,jcmax,jcmin,npy_cfd,kcmax,kcmin,npz_cfd
        ! number of grid per MPI task, remainder must be added !!!
        nixyz  = (/ (icmax - icmin) / npx_cfd + 1, (jcmax-jcmin) / npy_cfd + 1, (kcmax - kcmin) / npz_cfd + 1/)
        minxyz = (/ icmin,  jcmin,  kcmin /)
        npxyz  = (/ npx_cfd, npy_cfd, npz_cfd  /)

        allocate(bbox_cfd%xbb(2,npx_cfd),bbox_cfd%ybb(2,npy_cfd), bbox_cfd%zbb(2,npz_cfd))

        do ixyz = 1,3

            select case(ixyz)
            case(1)
                bb_ptr => bbox_cfd%xbb
            case(2)
                bb_ptr => bbox_cfd%ybb
            case(3) 
                bb_ptr => bbox_cfd%zbb
            end select

            bb_ptr(is,1) = minxyz(ixyz)
            bb_ptr(ie,1) = bb_ptr(is,1) + nixyz(ixyz) - 1

            do i = 2, npxyz(ixyz)
                bb_ptr(is, i) = bb_ptr(ie, i-1)
                bb_ptr(ie, i) = bb_ptr(is, i) + nixyz(ixyz) - 1
            enddo

        enddo

        ! set sizes of local grids
        nlgx_cfd = icPmax_cfd(rank2coord_cfd(1,rank_realm)) - icPmin_cfd(rank2coord_cfd(1,rank_realm))
        nlgy_cfd = jcPmax_cfd(rank2coord_cfd(2,rank_realm)) - jcPmin_cfd(rank2coord_cfd(2,rank_realm))
        nlgz_cfd = kcPmax_cfd(rank2coord_cfd(3,rank_realm)) - kcPmin_cfd(rank2coord_cfd(3,rank_realm))

		! THIS IS TERRIBLE NOTATION AS IT REDEFINES A CONTINUUM VALUE
        !nlx = bbox_cfd%xbb(2,rank2coord_cfd(1,id_coord)) - bbox_cfd%xbb(1,rank2coord_cfd(1,id_coord)) + 1
        !nly = min(bbox_cfd%ybb(2,rank2coord_cfd(2,id_coord)),jcmax_overlap) - bbox_cfd%ybb(1,rank2coord_cfd(2,id_coord)) + 1
        !nlz = bbox_cfd%zbb(2,rank2coord_cfd(3,id_coord)) - bbox_cfd%zbb(1,rank2coord_cfd(3,id_coord)) + 1

		!write(0,*)' CFD: bbox ', rank_realm, bbox_cfd%xbb,bbox_cfd%ybb,bbox_cfd%zbb, nlgx_cfd, nlgy_cfd, nlgz_cfd

    end subroutine make_bbox

    !-----------------------------------------------------------------------------
    ! Establish location and size of the overlap region between the continuum
    ! and molecular simulations
    !-----------------------------------------------------------------------------

    subroutine find_overlaps
        implicit none

        integer i, 	ibmin_cfd,ibmax_cfd,jbmin_cfd,jbmax_cfd,kbmin_cfd,kbmax_cfd, &
            		ibmin_md, ibmax_md, jbmin_md, jbmax_md, kbmin_md, kbmax_md

        ibmin_cfd = icPmin_cfd(iblock_realm) !bbox_cfd%xbb(1,rank2coord_cfd(1,id_coord))
        ibmax_cfd = icPmax_cfd(iblock_realm) !bbox_cfd%xbb(2,rank2coord_cfd(1,id_coord))
        jbmin_cfd = jcPmin_cfd(jblock_realm) !bbox_cfd%ybb(1,rank2coord_cfd(2,id_coord))
        jbmax_cfd = jcPmax_cfd(jblock_realm) !min(jcmax_overlap,bbox_cfd%ybb(2,rank2coord_cfd(2,id_coord))) 
        kbmin_cfd = kcPmin_cfd(kblock_realm) !bbox_cfd%zbb(1,rank2coord_cfd(3,id_coord))
        kbmax_cfd = kcPmax_cfd(kblock_realm) !bbox_cfd%zbb(2,rank2coord_cfd(3,id_coord))

        write(0, *) 'CFD: find overlap ibmin etc', rank_realm, ibmin_cfd, ibmax_cfd, jbmin_cfd, jbmax_cfd, kbmin_cfd, kbmax_cfd

        do i=0,nproc_md - 1

            ibmin_md = icPmin_md(iblock_realm)
            ibmax_md = icPmax_md(iblock_realm)
            jbmin_md = jcPmin_md(jblock_realm)
            jbmax_md = jcPmax_md(jblock_realm)
            kbmin_md = kcPmin_md(kblock_realm)
            kbmax_md = kcPmax_md(kblock_realm)

            if  ((( ibmin_md <  ibmin_cfd .and. ibmax_md > ibmin_cfd )	.or.  &
                (   ibmin_md >= ibmin_cfd .and. ibmin_md < ibmax_cfd )) .and. &
                ((  jbmin_md <  jbmin_cfd .and. jbmax_md > jbmin_cfd )	.or.  &
                (   jbmin_md >= jbmin_cfd .and. jbmin_md < jbmax_cfd))  .and. &   
                ((  kbmin_md <  kbmin_cfd .and. kbmax_md > kbmin_cfd )	.or.  &
                (   kbmin_md >= kbmin_cfd .and. kbmin_md < kbmax_cfd)))  then

				!This processor overlaps the MD domain
                overlap_mask(i) = 1

				!I think this ensures local min/max are replaced by global min/max values if appropriate
                overlap_box(1,i) = max(ibmin_cfd,ibmin_md)
                overlap_box(2,i) = min(ibmax_cfd,ibmax_md)
                overlap_box(3,i) = max(jbmin_cfd,jbmin_md)
                overlap_box(4,i) = min(ncy_olap,jbmax_md)
                overlap_box(5,i) = max(kbmin_cfd,kbmin_md)
                overlap_box(6,i) = min(kbmax_cfd,kbmax_md)

            else
                overlap_mask(i) = 0
                overlap_box(:,i) = -666

            endif

        enddo

    end subroutine find_overlaps

end subroutine create_map_cfd

end module coupler_internal_cfd
