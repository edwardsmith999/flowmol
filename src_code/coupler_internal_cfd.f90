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

subroutine coupler_cfd_init(nsteps,dt_cfd,icomm_grid,icoord,npxyz_cfd,xyzL,ngxyz,density, & 
							   ijkmax,ijkmin,iTmin,iTmax,jTmin,jTmax,kTmin,kTmax,xpg,ypg,zpg)
    use mpi
	use coupler, only : write_matrix,write_matrix_int
    use coupler_module,	dt_cfd_=>dt_cfd,						&	!Simulation lengths	
						iTmin_=>iTmin,iTmax_=>iTmax, 			&   !CFD local grid limits
						jTmin_=>jTmin,jTmax_=>jTmax, 			&   !CFD local grid limits
						kTmin_=>kTmin,kTmax_=>kTmax,			&   !CFD local grid limits
						xpg_=>xpg,ypg_=>ypg,zpg_=>zpg				!CFD grid arrays
	use coupler_input_data, only : cfd_coupler_input
 	implicit none			

    integer,					    intent(in)	:: nsteps,icomm_grid 
    integer,dimension(3),		    intent(in)	:: ijkmin,ijkmax,npxyz_cfd,ngxyz
	integer,dimension(:),		    intent(in)	:: iTmin,iTmax,jTmin,jTmax,kTmin,kTmax
    integer,dimension(:,:),		    intent(in)	:: icoord
    real(kind(0.d0)),			    intent(in)	:: dt_cfd,density
    real(kind(0.d0)),dimension(3),  intent(in)	:: xyzL
    real(kind(0.d0)),dimension(:  ),intent(in)	:: zpg
    real(kind(0.d0)),dimension(:,:),intent(in)	:: xpg,ypg


    integer											:: i, myid, source, ierr
    integer,dimension(:),allocatable				:: buf
	real(kind=kind(0.d0))							:: dxmin,dxmax,dzmin,dzmax
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

	write(99+myid,*), 'CFD side',myid,'CFD procs',npx_cfd,npy_cfd,npz_cfd,nproc_cfd

	! Receive & Store MD number of processors
	allocate(buf(3))
    call MPI_bcast(   buf   ,3,MPI_INTEGER,  0   ,COUPLER_ICOMM,ierr)	!Receive
    npx_md = buf(1)
    npy_md = buf(2)
    npz_md = buf(3)
    nproc_md = npx_md * npy_md * npz_md
	deallocate(buf)

	write(99+myid,*), 'CFD side',myid,'MD procs',npx_md ,npy_md ,npz_md ,nproc_md

	! Store & Send CFD processor topology
    allocate(icoord_cfd(3,nproc_cfd),stat=ierr); icoord_cfd = icoord
	allocate(buf(3*nproc_cfd)); buf = reshape(icoord, (/ 3*nproc_cfd /) )
    call MPI_bcast(buf,3*nproc_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr)	!Send
	deallocate(buf)

	call write_matrix_int(icoord_cfd,'cfd side, icoord_cfd=',99+myid)

	! Receive & Store MD processor topology
	allocate(buf(3*nproc_md))
    call MPI_bcast(buf,3*nproc_md ,MPI_INTEGER,  0   ,COUPLER_ICOMM,ierr)	!Receive
    allocate(icoord_md (3,nproc_md),stat=ierr); icoord_md = reshape(buf,(/ 3,nproc_md /))
	deallocate(buf)

	call write_matrix_int(icoord_md,'cfd side, icoord_md=',99+myid)

	! ------------------ Timesteps and iterations ------------------------------
	! Store & send CFD nsteps and dt_cfd
	nsteps_cfd = nsteps
    call MPI_bcast(nsteps,1,MPI_integer,source,COUPLER_ICOMM,ierr)			!Send
	dt_cfd_ = dt_cfd
    call MPI_bcast(dt_cfd,1,MPI_double_precision,source,COUPLER_ICOMM,ierr)	!Send

	write(99+myid,*), 'CFD side',myid,'CFD times', nsteps_cfd,dt_cfd

	! Receive & store MD timestep dt_md
    call MPI_bcast(dt_md,1,MPI_double_precision,0,COUPLER_ICOMM,ierr)		!Receive
    call MPI_bcast(nsteps_md,1,MPI_integer,     0,COUPLER_ICOMM,ierr)		!Receive

	write(99+myid,*), 'CFD side',myid,'MD times', nsteps_md,dt_md

	! ------------------ Send CFD grid extents ------------------------------

	! Store & send CFD density
	density_cfd = density
	call MPI_bcast(density_cfd,1,MPI_double_precision,source,COUPLER_ICOMM,ierr)	!Send

	write(99+myid,*), 'CFD side',myid,'CFD density',density_cfd

	! Receive & store MD density
	call MPI_bcast(density_md,1,MPI_double_precision,0,COUPLER_ICOMM,ierr)		!Receive

	write(99+myid,*), 'CFD side',myid,'CFD density',density_md

	! Store & send CFD domain size
	xL_cfd = xyzL(1); yL_cfd = xyzL(2); zL_cfd = xyzL(3)
	call MPI_bcast(xyzL,3,MPI_double_precision,source,COUPLER_ICOMM,ierr)	!Send

	write(99+myid,*), 'CFD side',myid,'CFD domain',xL_cfd,yL_cfd,zL_cfd	
	! Receive & store MD domain size
	allocate(rbuf(3))
	call MPI_bcast(rbuf,3,MPI_double_precision,0,COUPLER_ICOMM,ierr)		!Receive
	xL_md = rbuf(1); yL_md = rbuf(2); zL_md = rbuf(3);
	deallocate(rbuf)

	write(99+myid,*), 'CFD side',myid,'MD domain', xL_md,yL_md,zL_md

	! Store & send CFD grid extents
	imin = ijkmin(1); jmin = ijkmin(2); kmin = ijkmin(3)
	imax = ijkmax(1); jmax = ijkmax(2); kmax = ijkmax(3)
    call MPI_bcast((/ imin,imax,jmin,jmax,kmin,kmax /),6,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send

	write(99+myid,*), 'CFD side',myid,'CFD global extents',imin,imax,jmin,jmax,kmin,kmax

	! Store & send global number of cells in CFD
	ngx = ngxyz(1); ngy = ngxyz(2); ngz = ngxyz(3)
    call MPI_bcast(ngxyz,3,MPI_INTEGER,source,COUPLER_ICOMM,ierr)				!Send

	write(99+myid,*), 'CFD side',myid,'CFD global cells',ngx,ngy,ngz

	! Store & send array of global grid points
    allocate(xpg_(size(xpg,1),size(xpg,2)),stat=ierr); xpg_ = xpg
    allocate(ypg_(size(ypg,1),size(ypg,2)),stat=ierr); ypg_ = ypg
    allocate(zpg_(size(zpg,1)			 ),stat=ierr); zpg_ = zpg
    call MPI_bcast(xpg,size(xpg),MPI_double_precision,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(ypg,size(ypg),MPI_double_precision,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(zpg,size(zpg),MPI_double_precision,source,COUPLER_ICOMM,ierr) !Send

	call write_matrix(xpg,'cfd side, xpg=',50+myid)
	call write_matrix(ypg,'cfd side, ypg=',50+myid)
	write(50+myid,*), 'CFD side',myid,'zpg',zpg

    ! Store & Send local (processor) CFD grid extents
    allocate(iTmin_(npx_cfd),stat=ierr); iTmin_(:) = iTmin(:)
    allocate(iTmax_(npx_cfd),stat=ierr); iTmax_(:) = iTmax(:)
    allocate(jTmin_(npy_cfd),stat=ierr); jTmin_(:) = jTmin(:)
    allocate(jTmax_(npy_cfd),stat=ierr); jTmax_(:) = jTmax(:)
    allocate(kTmin_(npz_cfd),stat=ierr); kTmin_(:) = kTmin(:)
    allocate(kTmax_(npz_cfd),stat=ierr); kTmax_(:) = kTmax(:)
    call MPI_bcast(iTmin,npx_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(iTmax,npx_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(jTmin,npy_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(jTmax,npy_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(kTmin,npz_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send
    call MPI_bcast(kTmax,npz_cfd,MPI_INTEGER,source,COUPLER_ICOMM,ierr) !Send

	write(99+myid,*), 'CFD side',myid,'CFD local cells',iTmin,iTmax,jTmin,jTmax,kTmin,kTmax

    ! send CFD processor grid and overlap parameter
    ! Note: jmax_overlap default is provided in coupler_internal_cfd
    if (cfd_coupler_input%overlap%tag == CPL) then
        jmax_overlap =  jmin + cfd_coupler_input%overlap%y_overlap
    endif

	call MPI_barrier(MPI_COMM_WORLD,ierr)

	!Define cell sizes dx,dy & dz and check for grid stretching
	! - - x - -
	dx = xpg(2,1)-xpg(1,1)
	dxmax = maxval(xpg(2:ngx,2:ngy)-xpg(1:ngx-1,1:ngy-1))
	dxmin = minval(xpg(2:ngx,2:ngy)-xpg(1:ngx-1,1:ngy-1))
	if (dxmax-dx.gt.0.00001d0) call error_abort("ERROR - Grid stretching in x not supported")
	if (dx-dxmin.gt.0.00001d0) call error_abort("ERROR - Grid stretching in x not supported")
	! - - y - -
	dy = ypg(1,2)-ypg(1,1)
	dymax = maxval(ypg(2:ngx,2:ngy)-ypg(1:ngx-1,1:ngy-1))
	dymin = minval(ypg(2:ngx,2:ngy)-ypg(1:ngx-1,1:ngy-1))
	if (dymax-dy.gt.0.0001 .or. dy-dymin.gt.0.0001) then
        write(*,*) "********************************************************************"
        write(*,*) " Grid stretching employed in CFD domain - range of dy sizes:        "
		write(*,*) "dymin = ", dymin, " dy = ",dy, " dymax = ", dymax
        write(*,*) "********************************************************************"
        write(*,*)
	endif
	! - - z - -
	dz = zpg(2  )-zpg(1  )
	dzmax = maxval(zpg(2:ngz)-zpg(1:ngz-1))
	dzmin = minval(zpg(2:ngz)-zpg(1:ngz-1))
	if (dzmax-dz.gt.0.00001d0) call error_abort("ERROR - Grid stretching in z not supported")
	if (dz-dzmin.gt.0.00001d0) call error_abort("ERROR - Grid stretching in z not supported")

    ! test if MD_init_cell size is larger than CFD cell size
    if( MD_initial_cellsize .ge. dx .or. & 
		MD_initial_cellsize .ge. dy .or. & 
		MD_initial_cellsize .ge. dz .and. myid == 0 ) then
        write(*,*)
        write(*,*) "********************************************************************"
        write(*,*) " WARNING ...WARNING ...WARNING ...WARNING ...WARNING ...WARNING ... "
        write(*,*) " MD initialisation cell size larger than CFD x,y cell sizes         "
        write(*,*) " MD_init_cellsize=",MD_initial_cellsize
        write(*,'(3(a,f10.5))') " dx=",xpg(2,1)-xpg(1,1),  & 
								" dy=",ypg(1,2)-ypg(1,1),  & 
								" dz=",zpg(2  )-zpg(1  )
        write(*,*) "********************************************************************"
        write(*,*)
    endif

end subroutine coupler_cfd_init

!=============================================================================
! Establish for all CFD processors the mapping (if any) 
! to coupled MD processors
!-----------------------------------------------------------------------------

subroutine create_map_cfd
    use mpi 
    use coupler_module!, only : COUPLER_REALM_COMM, COUPLER_GRID_COMM, COUPLER_ICOMM, cfd_is_2d, map
    implicit none

    integer 				:: i, myid, id_coord, color, noverlaps, ir, iaux(4), ierr
    integer, allocatable 	:: md_grid_boxes(:,:), overlap_mask(:), ireq(:), overlap_box(:,:)
    real(kind(0.d0))		:: raux(2)

	call MPI_barrier(COUPLER_REALM_COMM,ierr)
    call mpi_comm_rank(COUPLER_REALM_COMM,myid,ierr)
	myid = myid + 1

    ! Coupler_grid_comm is the correct communicator to pick the processor coordinates in cartesian topology
    call mpi_comm_rank(coupler_grid_comm,id_coord,ierr)
    id_coord = id_coord+1 ! get in sync with icoord convention

    call make_bbox

	! Receive the overlapping box indices from all MD processors
    allocate(md_grid_boxes(6,0:nproc_md-1), overlap_mask(0:nproc_md-1), &
        	   overlap_box(6,0:nproc_md-1),         ireq(0:nproc_md-1))
	call mpi_allgather(MPI_BOTTOM,0,MPI_INTEGER,md_grid_boxes,6,MPI_INTEGER,COUPLER_ICOMM, ierr)
 
    ! Find overlaps and send domain overlap mask to all MD processors
    call find_overlaps
    call mpi_allgather(overlap_mask,nproc_md,MPI_INTEGER,MPI_BOTTOM,0,MPI_INTEGER,COUPLER_ICOMM,ierr)

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

    call mpi_comm_split(COUPLER_REALM_COMM, color, myid, CFD_COMM_OVERLAP, ierr)

    if (color == 0) then
        CFD_COMM_OVERLAP = MPI_COMM_NULL
    endif

    ! Send the range of the overlaping domains
    ir = 0
    do i=0, nproc_md-1
        if (overlap_mask(i) == 1) then
            call mpi_isend(overlap_box(1,i),6,mpi_integer,i,2,COUPLER_ICOMM,ireq(i),ierr)
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

		print*, imax,imin,npx_cfd,jmax,jmin,npy_cfd,kmax,kmin,npz_cfd
        ! number of grid per MPI task, remainder must be added !!!
        nixyz  = (/ (imax - imin) / npx_cfd + 1, (jmax-jmin) / npy_cfd + 1, (kmax - kmin) / npz_cfd + 1/)
        minxyz = (/ imin,  jmin,  kmin /)
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
        nlgx_cfd = iTmax(icoord_cfd(1,myid)) - iTmin(icoord_cfd(1,myid))
        nlgy_cfd = jTmax(icoord_cfd(2,myid)) - jTmin(icoord_cfd(2,myid))
        nlgz_cfd = kTmax(icoord_cfd(3,myid)) - kTmin(icoord_cfd(3,myid))

		! THIS IS TERRIBLE NOTATION AS IT REDEFINES A CONTINUUM VALUE
        !nlx = bbox_cfd%xbb(2,icoord_cfd(1,id_coord)) - bbox_cfd%xbb(1,icoord_cfd(1,id_coord)) + 1
        !nly = min(bbox_cfd%ybb(2,icoord_cfd(2,id_coord)),jmax_overlap) - bbox_cfd%ybb(1,icoord_cfd(2,id_coord)) + 1
        !nlz = bbox_cfd%zbb(2,icoord_cfd(3,id_coord)) - bbox_cfd%zbb(1,icoord_cfd(3,id_coord)) + 1

		write(0,*)' CFD: bbox ', myid, bbox_cfd%xbb,bbox_cfd%ybb,bbox_cfd%zbb, nlgx_cfd, nlgy_cfd, nlgz_cfd

    end subroutine make_bbox

    !-----------------------------------------------------------------------------
    ! Establish location and size of the overlap region between the continuum
    ! and molecular simulations
    !-----------------------------------------------------------------------------

    subroutine find_overlaps
        implicit none

        integer i, ibmin,ibmax,jbmin,jbmax,kbmin,kbmax, &
            	ibs, ibe, jbs, jbe, kbs, kbe

        ibmin = bbox_cfd%xbb(1,icoord_cfd(1,id_coord))
        ibmax = bbox_cfd%xbb(2,icoord_cfd(1,id_coord))
        jbmin = bbox_cfd%ybb(1,icoord_cfd(2,id_coord))
        jbmax = min(jmax_overlap,bbox_cfd%ybb(2,icoord_cfd(2,id_coord))) 
        kbmin = bbox_cfd%zbb(1,icoord_cfd(3,id_coord))
        kbmax = bbox_cfd%zbb(2,icoord_cfd(3,id_coord))

        write(0, *) 'CFD: find overlap ibmin etc', myid, ibmin, ibmax, jbmin, jbmax, kbmin, kbmax

        do i=0,nproc_md - 1

            ibs = md_grid_boxes(1,i)
            ibe = md_grid_boxes(2,i)
            jbs = md_grid_boxes(3,i)
            jbe = md_grid_boxes(4,i)
            kbs = md_grid_boxes(5,i)
            kbe = md_grid_boxes(6,i)

            if  ((( ibs <  ibmin .and. ibe > ibmin )	.or.  &
                (   ibs >= ibmin .and. ibs < ibmax ))  .and. &
                ((  jbs <  jbmin .and. jbe > jbmin )	.or.  &
                (   jbs >= jbmin .and. jbs < jbmax))   .and. &   
                ((  kbs <  kbmin .and. kbe > kbmin )	.or.  &
                (   kbs >= kbmin .and. kbs < kbmax)))  then

				!This processor overlaps the MD domain
                overlap_mask(i) = 1

				!I think this ensures local min/max are replaced by global min/max values if appropriate
                overlap_box(1,i) = max(ibmin,ibs)
                overlap_box(2,i) = min(ibmax,ibe)
                overlap_box(3,i) = max(jbmin,jbs)
                overlap_box(4,i) = min(jmax_overlap,jbe)
                overlap_box(5,i) = max(kbmin,kbs)
                overlap_box(6,i) = min(kbmax,kbe)

            else
                overlap_mask(i) = 0
                overlap_box(:,i) = -666

            endif

        enddo

        write(0,*)' CFD: overlap ', myid, overlap_mask, overlap_box

    end subroutine find_overlaps

end subroutine create_map_cfd

end module coupler_internal_cfd
