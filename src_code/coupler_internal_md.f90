!=============================================================================
!				   Coupler internal MD
! Internal data and subroutines used by the coupler when working in MD realm
! It must not be used by subroutimes working in CFD realm ONLY
! Subroutines include:
!
! create_map_md 			Establish for all MD processors the mapping (if any) 
! 							to coupled CFD processors
! make_bbox					Make bbox which contains all domain information 
!							required in coupling process such as domain extents
! get_CFDvel				Get velocity fields from CFD for the constraint 
!							force needed in MD
! send_vel(a,n1,n2,n3)   	Send the average MD bin velocities to the 
!							corresponding rank in the CFD realm
! write_overlap_map			Write the map of coupled processors for debugging
! function map_md2cfd(r)	Get global molecular position from local position
! 
!  Lucian Anton, November 2011
!
!=============================================================================

module coupler_internal_md
    use coupler_parameters
	implicit none
    save 

	! data structures to hold CFD - MD mapping 
	! bounding box of MD domain within CFD global domain
	type bbox_domain
		integer is,ie,js,je,ks,ke
        integer iso,ieo,jso,jeo,kso,keo ! overlap indices
		real(kind=kind(0.d0)) :: bb(2,3)
	end type bbox_domain

	type(bbox_domain), target :: bbox

	! thickness of the MD region below the CFD grid 
	real(kind=kind(0.d0)) DY_PURE_MD

	! local domain lenghts, and halves
	real(kind=kind(0.d0)) domain_lengths(3), half_domain_lengths(3)

    ! data made available to MD callers
    type(cfd_grid_info) cfd_box
 
	! write or not the overlap map
	logical :: dump_overlap_map = .true.

	type cfd_box_sum
		integer np
		real(kind=kind(0.d0))  v(3)
		real(kind=kind(0.d0))  a(3)
	end type cfd_box_sum

	real(kind=kind(0.d0)) :: FoP_time_ratio = 1.0   ! time ratio dt_CFD/dt_MD; to be fixed later
	real(kind=kind(0.d0)) :: fsig=1.0  		!Ratio betwen macroscopic unit lenght and molecular unit 
    integer               :: md_steps_per_dt_cfd = -1 ! number of steps per CFD step
													   ! if not defined in COUPLER.in 
                                                       ! nsteps is set dt_cfd/dt_md (in MD socket)

	! array for velocities from CFD, last dimension holds time indices 
	real(kind=kind(0.d0)),dimension(:,:,:,:,:), allocatable :: vel_fromCFD
	integer itm1,itm2

    ! CFD code id
    integer :: cfd_code_id = couette_parallel

contains 

! ----------------------------------------------------------------------------
! Initialisation routine for coupler - Every variable is sent and stored
! to ensure both md and cfd region have an identical list of parameters

subroutine coupler_md_init(nsteps,dt_md,icomm_grid,icoord,npxyz_md,globaldomain,density)
	use mpi
	use coupler_input_data, cfd_code_id_in => cfd_code_id, density_coupled => density
    use coupler_module, dt_md_=>dt_md		
	implicit none

	integer, intent(in)	  							:: nsteps, icomm_grid
	integer,dimension(3), intent(in)	  			:: npxyz_md	
	integer,dimension(:,:),allocatable,intent(in)	:: icoord
	real(kind(0.d0)),intent(in) 					:: dt_md,density
    real(kind=kind(0.d0)),dimension(3),intent(in) 	:: globaldomain

    integer											:: i,ib,jb,kb,pcoords(3),source
    integer,dimension(:),allocatable  				:: buf
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

	write(999+rank_realm,*), 'MD side',rank_realm,'CFD procs',npx_cfd,npy_cfd,npz_cfd,nproc_cfd

	! Store & Send MD number of processors
	npx_md = npxyz_md(1);	npy_md = npxyz_md(2);	npz_md = npxyz_md(3)	
	nproc_md = npx_md * npy_md * npz_md
	call MPI_bcast(npxyz_md,3,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send

	write(999+rank_realm,*), 'MD side',rank_realm,'MD procs',npx_md ,npy_md ,npz_md ,nproc_md

	! Receive & Store CFD processor rank to coord
	allocate(buf(3*nproc_cfd))
    call MPI_bcast(buf,3*nproc_cfd,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr) !Receive
    allocate(rank2coord_cfd(3,nproc_cfd),stat=ierr); rank2coord_cfd = reshape(buf,(/ 3,nproc_cfd /))
	deallocate(buf)

	call write_matrix_int(rank2coord_cfd,'MD side, rank2coord_cfd=',999+rank_realm)

	! Store & Send MD processor rank to coord
    allocate(rank2coord_md(3,nproc_md),stat=ierr); rank2coord_md = icoord
	iblock_realm=icoord(1,rank_realm); jblock_realm=icoord(2,rank_realm); kblock_realm=icoord(3,rank_realm)
	allocate(buf(3*nproc_md)); buf = reshape(icoord,(/ 3*nproc_md /))
    call MPI_bcast(buf ,3*nproc_md,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
	deallocate(buf)

	call write_matrix_int(rank2coord_md,'MD side, rank2coord_md=',999+rank_realm)

	! Receive & Store CFD coordinate to rank mapping
	allocate(buf(nproc_cfd))
    call MPI_bcast(buf,nproc_cfd ,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)	!Receive
    allocate(coord2rank_cfd (npx_cfd,npy_cfd,npz_cfd))
	coord2rank_cfd = reshape(buf,(/ npx_cfd,npy_cfd,npz_cfd /))
	deallocate(buf)

	write(999+rank_realm,*), 'MD side',rank_realm, 'coord2rank_cfd=', coord2rank_cfd

	!Setup MD mapping from coordinate to rank, store and send
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

	write(999+rank_realm,*), 'MD side',rank_realm, 'coord2rank_md=', coord2rank_md

	! ------------------ Timesteps and iterations ------------------------------
	! Receive & store CFD nsteps and dt_cfd
	call MPI_bcast(nsteps,1,MPI_integer,0,CPL_INTER_COMM,ierr)				!Receive
	call MPI_bcast(dt_cfd,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive

	! Store & send MD timestep to dt_md
	dt_MD_ = dt_MD
    call MPI_bcast(dt_md,1,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send
	nsteps_MD = nsteps
    call MPI_bcast(nsteps,1,MPI_integer,source,CPL_INTER_COMM,ierr)	!Send

	write(999+rank_realm,*), 'MD side',rank_realm,'CFD times',nsteps_cfd,dt_cfd
	write(999+rank_realm,*), 'MD side',rank_realm,'MD times', nsteps_MD,dt_md

	! ------------------ Receive CFD grid extents ------------------------------
	! Receive & store CFD density
	call MPI_bcast(density_cfd,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive

	write(999+rank_realm,*), 'MD side',rank_realm,'CFD density',density_cfd

	! Store & send MD density
	density_md = density 
	call MPI_bcast(density,1,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	write(999+rank_realm,*), 'MD side',rank_realm,'MD density',density_md

	! Receive & store CFD domain size
	allocate(rbuf(3))
	call MPI_bcast(rbuf,3,MPI_double_precision,0,CPL_INTER_COMM,ierr)				!Receive
	xL_cfd = rbuf(1); yL_cfd = rbuf(2); zL_cfd = rbuf(3)
	deallocate(rbuf)

	! Store & send MD domain size
	xL_md = globaldomain(1); yL_md = globaldomain(2); zL_md = globaldomain(3) 
	call MPI_bcast(globaldomain,3,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	write(999+rank_realm,*), 'MD side',rank_realm,'CFD domain',xL_cfd,yL_cfd,zL_cfd
	write(999+rank_realm,*), 'MD side',rank_realm,'MD domain', xL_md,yL_md,zL_md


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

	write(999+rank_realm,*), 'MD side',rank_realm,'CFD global extents',icmin,icmax,jcmin,jcmax,kcmin,kcmax 
	write(999+rank_realm,*), 'MD side',rank_realm,'CFD global cells',ncx,ncy,ncz

	! Receive & Store array of global grid points
	allocate(xg(ncx+1,ncy+1),yg(ncx+1,ncy+1),zg(ncz+1))
	call MPI_bcast(xg,size(xg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive
	call MPI_bcast(yg,size(yg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive
	call MPI_bcast(zg,size(zg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive

	call write_matrix(xg,'MD side, xg=',500+rank_realm)
	call write_matrix(yg,'MD side, yg=',500+rank_realm)
	write(500+rank_realm,*), 'MD side',rank_realm,'zg',zg

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

	write(999+rank_realm,*), 'md side',rank_realm,'CFD local cells', & 
							icPmin_cfd,icPmax_cfd,jcPmin_cfd,jcPmax_cfd,kcPmin_cfd,kcPmax_cfd

	!Calculate the cell sizes dx,dy & dz
	dx = xL_cfd/ncx	  !xg(2,1)-xg(1,1)
	dy = yL_cfd/ncy	  !yg(1,2)-yg(1,1)
	dz = zL_cfd/ncz	  !zg(2  )-zg(1  )

    ! Initialise other md module variables if data is provided in coupler.in
	if (md_average_period_tag == CPL) then 
		average_period = md_average_period
	endif
	if (md_save_period_tag == CPL) then
		save_period    = md_save_period
	endif

	!if (cfd_code_id_tag == CPL) then
	!	cfd_code_id  = cfd_code_id_in
	!endif
    
	if ( nsteps_md <= 0 ) then 
		write(0,*) "Number of MD steps per dt interval <= 0"
		write(0,*) "Coupler will not work, quitting ..."
		call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_INIT,ierr)
	endif

	!Receive overlap from the CFD
    call MPI_bcast(ncy_olap,1,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive

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

!=============================================================================
! Establish for all MD processors the mapping (if any) 
! to coupled CFD processors
!-----------------------------------------------------------------------------

subroutine create_map_md
	use mpi
	use coupler_module, jcmax_overlap=> ncy_olap !, only : CPL_INTER_COMM, CPL_REALM_COMM, cfd_is_2d, map
	implicit none

	integer  i, ir, ireq(nproc_cfd), noverlaps, source, jcmino
	integer, allocatable :: overlap_mask(:,:)

	call MPI_barrier(CPL_REALM_COMM,ierr)

	jcmino = jcmin-1

	! compute the boundaries of this MD domain in the CFD global domain.
	! assume that all domains have identical sides 
	domain_lengths(:) = (/ xL_md/npx_md, yL_md/npy_md, zL_md/npz_md /)
	half_domain_lengths(:) = 0.5d0 * domain_lengths(:)

	! bounding boxes coordinates start from x(icmin), z(kcmin) and y(jcmino)-DY_PURE_MD
	bbox%bb(1,:) = (rank2coord_md(:,rank_realm)-1) * domain_lengths(:) &
					+ (/ xg(1,icmin), yg(jcmin,1)-DY_PURE_MD, zg(kcmin) /)
	bbox%bb(2,:) =  bbox%bb(1,:) + domain_lengths(:)

	call make_bbox

	write(0,*) 'MD: bbox%is ', rank_realm, jcmax_overlap, bbox%is, bbox%ie, bbox%js, bbox%je, bbox%ks, bbox%ke 

	! Send the overlapping box indices to CFD processors
	call mpi_allgather((/ bbox%iso, bbox%ieo, bbox%jso, bbox%jeo, bbox%kso, bbox%keo /), & 
							 6,MPI_INTEGER,MPI_BOTTOM,0,MPI_INTEGER,CPL_INTER_COMM,ierr)

	! Receive domain overlap mask from CFD processors
	allocate(overlap_mask(0:nproc_md-1,0:nproc_cfd-1))
	call mpi_allgather(MPI_BOTTOM,0,MPI_INTEGER,overlap_mask,nproc_md,MPI_INTEGER,CPL_INTER_COMM,ierr)
	noverlaps = 0
	do i = 0, nproc_cfd - 1
		if ( overlap_mask(rank_realm,i) == 1) then 
			noverlaps = noverlaps + 1
		endif
	enddo

	! sort out which CFD ranks hold non-void domains for this MD rank
	map%n = noverlaps
	allocate ( map%rank_list(noverlaps), map%domains(6,noverlaps))
	ir=0
	do i=0, nproc_cfd - 1
		if (overlap_mask(rank_realm,i) == 1) then
			ir = ir + 1
			map%rank_list(ir) = i
		endif
	enddo

    ! Receive the range of the overlapping domains
	do i =1, ir
		call mpi_irecv(map%domains(1,i), 6, MPI_INTEGER,map%rank_list(i),2, CPL_INTER_COMM,ireq(i),ierr)
	enddo
	call mpi_waitall(ir,ireq,MPI_STATUSES_IGNORE,ierr)

	if ( dump_overlap_map) then
		call write_overlap_map
	endif

    ! set the CFD box information for MD application
    cfd_box%gicmin = icmin
    cfd_box%gicmax = icmax
    cfd_box%gjcmin = jcmin
    cfd_box%gjcmax = jcmax
    cfd_box%gkcmin = kcmin
    cfd_box%gkcmax = kcmax
    cfd_box%icmin = bbox%is
    cfd_box%icmax = bbox%ie
    cfd_box%jcmin = bbox%js
    cfd_box%jcmax = bbox%je
    cfd_box%kcmin = bbox%ks
    cfd_box%kcmax = bbox%ke
    cfd_box%icmino = bbox%iso
    cfd_box%icmaxo = bbox%ieo
    cfd_box%jcmino = bbox%jso
    cfd_box%jcmaxo = bbox%jeo
    cfd_box%kcmino = bbox%kso
    cfd_box%kcmaxo = bbox%keo
    cfd_box%xmin = xg(1,bbox%is) - bbox%bb(1,1) - half_domain_lengths(1)
    cfd_box%xmax = xg(1,bbox%ie) - bbox%bb(1,1) - half_domain_lengths(1)
    cfd_box%dx   = xg(1,icmin+1) - xg(1,icmin)
    cfd_box%ymin = bbox%bb(1,2)
    cfd_box%ymax = bbox%bb(2,2)
    allocate(cfd_box%y(cfd_box%jcmin:cfd_box%jcmax))
    cfd_box%y(cfd_box%jcmin:cfd_box%jcmax) = yg(cfd_box%jcmin:cfd_box%jcmax,1)-bbox%bb(1,2) - half_domain_lengths(2)
    cfd_box%zmin = zg(bbox%ks) - bbox%bb(1,3) - half_domain_lengths(3)
    cfd_box%zmax = zg(bbox%ke) - bbox%bb(1,3) - half_domain_lengths(3)
    cfd_box%dz   = zg(kcmin+1) - zg(kcmin)

contains 

!-----------------------------------------------------------------------------
! Make bbox which contains all domain information required in coupling process
! such as domain extents
!-----------------------------------------------------------------------------

	subroutine make_bbox
		implicit none

		integer, parameter :: is = 1, ie = 2

		type grid_pointer
			real(kind=kind(0.d0)) , pointer :: p(:) => null()
		end type grid_pointer

		type bbox_pointer
			integer, pointer :: starto => null(), start => null(), endo => null(), end => null()
		end type bbox_pointer

		type(grid_pointer) grid_ptr(3)
		type(bbox_pointer) bbox_ptr(3)

		integer ixyz, ngp, grid_sizes(2,3), halo_size(2,3), idmin(3)
		real(kind=kind(0.d0)) pl,pr,eps  ! left right grid points
		logical found_start

		! indices covering the CFD physical domain
		grid_sizes(:, :) = reshape((/ icmin, icmax, jcmin, jcmax_overlap, kcmin, kcmax /),(/2,3/))

		! starting indices (first in physical domain) in all directions
		idmin = (/ icmin, jcmin, kcmin /)

		! how large is the halo in each direction, depending also on the MD domain position
        ! this is to ensure that all MD domain is covered, also useful to collect quantites
        ! from particles that have left the MD domain
		halo_size(:,:) = 0

		! pointer to grid coordinates. Helpful to loop through dimensions
		grid_ptr(1)%p => xg(1,icmin:icmax)
		grid_ptr(2)%p => yg(jcmin:jcmax_overlap,1)
		grid_ptr(3)%p => zg(kcmin:kcmax)

		bbox_ptr(1)%starto => bbox%iso
        bbox_ptr(1)%start  => bbox%is
		bbox_ptr(1)%endo   => bbox%ieo
        bbox_ptr(1)%end    => bbox%ie
		bbox_ptr(2)%starto => bbox%jso
        bbox_ptr(2)%start  => bbox%js
		bbox_ptr(2)%endo   => bbox%jeo
        bbox_ptr(2)%end    => bbox%je
		bbox_ptr(3)%starto => bbox%kso
        bbox_ptr(3)%start  => bbox%ks
		bbox_ptr(3)%endo   => bbox%keo
        bbox_ptr(3)%end    => bbox%ke


		!Loop through all 3 Dimensions
        do ixyz=1,3

			!Define small number
			eps = 1.d-2 * (grid_ptr(ixyz)%p(2) - grid_ptr(ixyz)%p(1))

			!Check for special case of first processor containing first cell in global domain
			found_start = .false.
			if ( grid_ptr(ixyz)%p(1) >= bbox%bb(1,ixyz) - eps &
				.and. grid_ptr(ixyz)%p(1) < bbox%bb(2,ixyz) ) then 
				found_start = .true.
				bbox_ptr(ixyz)%starto = idmin(ixyz) - halo_size(1,ixyz)
                bbox_ptr(ixyz)%start  = idmin(ixyz)
				!write(0,*) "MD make box l", rank_realm,ixyz,bbox_ptr(ixyz)%start
			endif

			! On each processor, loop through all (global) cells in the domain (ngp) and store,
			! when found, the start and end cells covered by current processes' local domain extents
			ngp = grid_sizes(2,ixyz) - grid_sizes(1,ixyz) + 1
			do i=2, ngp
				pl = grid_ptr(ixyz)%p(i-1)
				pr = grid_ptr(ixyz)%p(i) 
				!write(0,*), 'MD make bbox ', rank_realm, ixyz,ngp,i,pl,pr, bbox%bb(:,ixyz) 
				!If processor starting cell is not yet found, continue to check
				if (.not. found_start )then
					!Check if current cell is processor starting cell
					if ( pl < bbox%bb(1,ixyz)  .and. pr >=  bbox%bb(1,ixyz) - eps ) then 
						found_start = .true.
						bbox_ptr(ixyz)%starto = idmin(ixyz) + i - 1 - halo_size(1,ixyz)
                        bbox_ptr(ixyz)%start  = idmin(ixyz) + i - 1
						!write(0,*), 'MD make bbox l', rank_realm, ixyz,i, pl, pr, bbox_ptr(ixyz)%start
					endif
				!Else, check if final cell in domain is found yet
				else
					if ( (i < ngp  .and. pl <= bbox%bb(2,ixyz) + eps  .and. pr > bbox%bb(2,ixyz))) then
						bbox_ptr(ixyz)%endo = idmin(ixyz) + i - 1 - 1 + halo_size(2,ixyz)
                        bbox_ptr(ixyz)%end  = idmin(ixyz) + i - 1 - 1
						!write(0,*), 'MD make bbox r', rank_realm, ixyz, i, pl, pr ,  bbox_ptr(ixyz)%end		     
						exit
					!Special case of final processor in domain containing last cell
					else if (i == ngp  .and. abs( pr - bbox%bb(2,ixyz)) < eps ) then 
						bbox_ptr(ixyz)%endo = idmin(ixyz) + ngp - 1 + halo_size(2,ixyz)
                        bbox_ptr(ixyz)%end  = idmin(ixyz) + ngp - 1
						!write(0,*), 'MD make bbox r', rank_realm, ixyz,i, pl, pr, bbox_ptr(ixyz)%end
					endif
				endif
			enddo

		enddo

        ! Sanity check
        ! Test if the local grid contains local MD domain
        ! that is: x(bbox%is) < bboxbb(1,1) ; bbox%bb(2,1) < x(bbox%ie) ...


		! Set maximum local grid sizes
		nlgx_md = bbox%ieo - bbox%iso + 1
		nlgy_md = bbox%jeo - bbox%jso + 1
		nlgz_md = bbox%keo - bbox%kso + 1

		write(0,*)' MD: bbox ', rank_realm, bbox, nlgx_md, nlgy_md, nlgz_md

	end subroutine make_bbox

end subroutine create_map_md

!=============================================================================
! Get MD position in CFD coordinate frame
!-----------------------------------------------------------------------------
function  map_md2cfd(r) result(rg)
	implicit none

	real(kind(0.d0)),intent(in) :: r(3)
	real(kind(0.d0)) rg(3)

	rg(:) = r(:) + half_domain_lengths(:) + bbox%bb(1,:)

end function map_md2cfd


!=============================================================================
! Get CFD position in MD coordinate frame
!-----------------------------------------------------------------------------
function map_cfd2md(r) result(rg)
	implicit none

	real(kind(0.d0)),intent(in) :: r(3)
	real(kind(0.d0)) rg(3)

	rg(:) = r(:) - half_domain_lengths(:) - bbox%bb(1,:)
	!print'(a,12f10.5)', 'MAP', r(:), half_domain_lengths(:), bbox%bb(1,:), rg

end function map_cfd2md


!=============================================================================
! Write to file the map of coupled processors for debugging
!-----------------------------------------------------------------------------

subroutine write_overlap_map
	use mpi
	use coupler_module
	implicit none

	integer, parameter :: max_msg_length = 1024
	character(len=*),parameter :: filename = "coupler_map.log"
	character(len=*),parameter :: nl = achar(10) ! new line; portable?
	
	character(len=max_msg_length) msg, rec
	character(len=100) fmtstr
	
	integer fh, i, n

	call mpi_file_open(CPL_REALM_COMM,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
	call mpi_file_set_size(fh,0_MPI_OFFSET_KIND,ierr) ! discard previous data

	n = map%n
	write(rec,'(a,I5,a,I5)')  'rank_world', rank_realm, ' overlaps', map%n
	msg = rec(1:len_trim(rec))//nl
	write(fmtstr,*) "( a,",1,"I5,","a,", 6,"I5)"
	
	do i = 1, n
		write(rec,fmtstr) 'rank', map%rank_list(i), &
			  ' indicies',map%domains(1:6,i)
	 
		msg = msg(1:len_trim(msg))//rec(1:len_trim(rec))//nl
	enddo
	
	msg = msg(1:len_trim(msg))//nl
	n = len_trim(msg)

	call mpi_file_write_shared(fh,msg(1:n),n,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr) 
	call mpi_barrier(CPL_REALM_COMM,ierr)
	call mpi_file_close(fh,ierr)

end subroutine write_overlap_map


end module coupler_internal_md
