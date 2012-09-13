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

	! local mpi ranks in realm and in topology communicators
	integer myid, myid_grid

	! MD data	
	integer npx, npy, npz, nproc
	integer, allocatable :: icoord(:,:)
	real(kind=kind(0.d0)) :: dt_MD

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
	logical :: dump_ovlerlap_map = .true.

	! local grid sizes
	integer nlx, nly, nlz

	type cfd_box_sum
		integer np
		real(kind=kind(0.d0))  v(3)
		real(kind=kind(0.d0))  a(3)
	end type cfd_box_sum

	integer			, dimension(:,:,:,:), allocatable :: mflux
	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: uc_bin, vc_bin, wc_bin
	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: uvwbin 

	real(kind=kind(0.d0)) :: FoP_time_ratio = 1.0   ! time ratio dt_CFD/dt_MD; to be fixed later
	real(kind=kind(0.d0)) :: xL_md, yL_md, zL_md 	! macroscopic sizes of MD domain. needed?
	real(kind=kind(0.d0)) :: fsig=1.0  		!Ratio betwen macroscopic unit lenght and molecular unit 
    integer               :: md_steps_per_dt_cfd = -1 ! number of steps per CFD step
													   ! if not defined in COUPLER.in 
                                                       ! nsteps is set dt_cfd/dt_md (in MD socket)

	! array for velocities from CFD, last dimension holds time indices 
	real(kind=kind(0.d0)),dimension(:,:,:,:,:), allocatable :: vel_fromCFD
	integer itm1,itm2

	! data from CFD
	integer :: npx_cfd, npy_cfd, npz_cfd, jmax_overlap_cfd, nproc_cfd
	! CFD mesh data
	integer :: imino, imin_cfd, imax_cfd,imaxo, jmino, jmin_cfd, jmax_cfd, jmaxo,&
		kmino, kmin_cfd, kmax_cfd, kmaxo
	real(kind=kind(0.d0)), allocatable, target :: x(:), y(:), z(:)
	real(kind=kind(0.d0)) :: dx, dz, dt_CFD, density

    ! CFD code id
    integer :: cfd_code_id = couette_parallel

    !inital MD cell size computed in CFD and used to adjust CFD domain sizes
    ! to integer multiples of cell size  
    real(kind(0.d0)) MD_initial_cellsize

	! nsteps from CFD
	integer :: nsteps
	! average period for averages ( it must come from CFD !!!)      
	integer :: average_period = 1
    ! save period ( corresponts to tplot in CFD, revise please !!!)
    integer :: save_period = 10


contains 

!=============================================================================
! Establish for all MD processors the mapping (if any) 
! to coupled CFD processors
!-----------------------------------------------------------------------------

subroutine create_map_md
	use mpi
	use coupler_internal_common, only : COUPLER_ICOMM, cfd_is_2d, map
	implicit none

	integer  i, ir, ireq(nproc_cfd), noverlaps, source, ierr
	integer, allocatable :: overlap_mask(:,:)

	! compute the boundaries of this MD domain in the CFD global domain.
	! assume that all domains have identical sides 
	domain_lengths(:) = (/ xL_md/npx, yL_md/npy, zL_md/npz /)
	half_domain_lengths(:) = 0.5d0 * domain_lengths(:)

	! bounding boxes coordinates start from x(imin), z(kmin) and y(jmino)-DY_PURE_MD
	bbox%bb(1,:) = (icoord(:,myid_grid)-1) * domain_lengths(:) &
					+ (/ x(imin_cfd), y(jmino)-DY_PURE_MD, z(kmin_cfd) /)
	bbox%bb(2,:) =  bbox%bb(1,:) + domain_lengths(:)

	!print*, 'Create Maps', xL_md/npx, yL_md/npy, zL_md/npz,myid_grid, (icoord(:,myid_grid)-1),  & 
	!					   (icoord(:,myid_grid)-1)*domain_lengths(:), x(imin_cfd), y(jmino)-DY_PURE_MD, z(kmin_cfd), domain_lengths, &
	!						bbox%bb(1,:),bbox%bb(1,:)

	! for 2D CFD problem one must broadcast  back the zL_md,z,dz
    ! because MD sets it
    if (cfd_is_2d) then
        if ( myid == 0 ) then
            source=MPI_ROOT
        else
            source=MPI_PROC_NULL
        endif
        call mpi_bcast(z, size(z), MPI_DOUBLE_PRECISION, source, COUPLER_ICOMM,ierr)
        call mpi_bcast((/kmino,kmin_cfd,kmax_cfd,kmaxo/), 4, MPI_INTEGER, source, COUPLER_ICOMM,ierr)
    endif
	!write(0,*) 'MD: bbox%bb ', myid, bbox%bb !, domain, npx, npy, npz, xL_md, yL_md, zL_md, icoord_md

	call make_bbox

	write(0,*) 'MD: bbox%is ', myid, jmax_overlap_cfd, bbox%is, bbox%ie, bbox%js, bbox%je, bbox%ks, bbox%ke 

	! Send the overlapping box indices to CFD processors
	call mpi_allgather((/ bbox%iso, bbox%ieo, bbox%jso, bbox%jeo, bbox%kso, bbox%keo /), & 
							 6,MPI_INTEGER,MPI_BOTTOM,0,MPI_INTEGER,COUPLER_ICOMM,ierr)

	! Receive domain overlap mask from CFD processors
	allocate(overlap_mask(0:nproc-1,0:nproc_cfd-1))
	call mpi_allgather(MPI_BOTTOM,0,MPI_INTEGER,overlap_mask,nproc,MPI_INTEGER,COUPLER_ICOMM,ierr)

	!		write(0,'(a,32I3)') 'MD, overlap mask: ', overlap_mask

	noverlaps = 0
	do i = 0, nproc_cfd - 1
		if ( overlap_mask(myid,i) == 1) then 
			noverlaps = noverlaps + 1
		endif
	enddo

	!		write(0,'(a,32I3)') 'MD, noverlaps: ', myid, noverlaps

	! sort out which CFD ranks hold non-void domains for this MD rank
	map%n = noverlaps
	allocate ( map%rank_list(noverlaps), map%domains(6,noverlaps))
	ir=0
	do i=0, nproc_cfd - 1
		if (overlap_mask(myid,i) == 1) then
			ir = ir + 1
			map%rank_list(ir) = i
		endif
	enddo

	do i =1, ir
		call mpi_irecv(map%domains(1,i), 6, MPI_INTEGER,map%rank_list(i),2, COUPLER_ICOMM,ireq(i),ierr)
	enddo
	call mpi_waitall(ir,ireq,MPI_STATUSES_IGNORE,ierr)

	if ( dump_ovlerlap_map) then
		call write_overlap_map
	endif

	! allocate array for CFD velocities and initialize the time indices; 
	!allocate(vel_fromCFD(3,bbox%ie - bbox%is + 1, bbox%je - bbox%js + 1 , bbox%ke - bbox%ks + 1 ,2))
	!vel_fromCFD = 0.d0
	!itm1 = 2; itm2 = 1

    ! set the CFD box information for MD application
    cfd_box%gimin = imin_cfd
    cfd_box%gimax = imax_cfd
    cfd_box%gjmin = jmin_cfd
    cfd_box%gjmax = jmax_cfd
    cfd_box%gkmin = kmin_cfd
    cfd_box%gkmax = kmax_cfd
    cfd_box%imin = bbox%is
    cfd_box%imax = bbox%ie
    cfd_box%jmin = bbox%js
    cfd_box%jmax = bbox%je
    cfd_box%kmin = bbox%ks
    cfd_box%kmax = bbox%ke
    cfd_box%imino = bbox%iso
    cfd_box%imaxo = bbox%ieo
    cfd_box%jmino = bbox%jso
    cfd_box%jmaxo = bbox%jeo
    cfd_box%kmino = bbox%kso
    cfd_box%kmaxo = bbox%keo
    cfd_box%xmin = x(bbox%is) - bbox%bb(1,1) - half_domain_lengths(1)
    cfd_box%xmax = x(bbox%ie) - bbox%bb(1,1) - half_domain_lengths(1)
    cfd_box%dx   = x(imin_cfd+1) - x(imin_cfd)
    cfd_box%ymin = bbox%bb(1,2)
    cfd_box%ymax = bbox%bb(2,2)
    allocate(cfd_box%y(cfd_box%jmin:cfd_box%jmax))
    cfd_box%y(cfd_box%jmin:cfd_box%jmax) = y(cfd_box%jmin:cfd_box%jmax)-bbox%bb(1,2) - half_domain_lengths(2)
    cfd_box%zmin = z(bbox%ks) - bbox%bb(1,3) - half_domain_lengths(3)
    cfd_box%zmax = z(bbox%ke) - bbox%bb(1,3) - half_domain_lengths(3)
    cfd_box%dz   = z(kmin_cfd+1) - z(kmin_cfd)

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

		integer id, ngp, grid_sizes(2,3), halo_size(2,3), idmin(3),ierr
		real(kind=kind(0.d0)) pl,pr,eps  ! left right grid points
		logical found_start

		! indices covering the CFD physical domain
		grid_sizes(:, :) = reshape((/ imin_cfd, imax_cfd, jmin_cfd, jmax_overlap_cfd, kmin_cfd, kmax_cfd /),(/2,3/))

		! starting indices (first in physical domain) in all directions
		idmin = (/ imin_cfd, jmin_cfd, kmin_cfd /)

		! write(0,*) 'MD: make box grid_sizes, idmin ', grid_sizes, idmin

		! how large is the halo in each direction, depending also on the MD domain position
        ! this is to ensure that all MD domain is covered, also useful to collect quantites
        ! from particles that have left the MD domain
		halo_size(:,:) = 1

		! special values for the boundaries
		if ( icoord(2,myid_grid) == 1 ) then 
			halo_size(1,2) = 0
		endif
		if ( icoord(2,myid_grid) == npy ) then
			halo_size(2,2) = 0
		endif
		if (icoord(1, myid_grid) == 1) then
			halo_size(1,1) = 0
		endif
		if (icoord(1, myid_grid) == npx) then 
			halo_size(2,1) = 0
		endif
		if (icoord(3, myid_grid) == 1) then 
			halo_size(1,3) = 0
		endif
		if (icoord(3, myid_grid) == npz) then 
			halo_size(2,3) = 0
		endif

		! pointer to grid coordinates. Helpful to loop through dimensions
		grid_ptr(1)%p => x(imin_cfd:imax_cfd)
		grid_ptr(2)%p => y(jmin_cfd:jmax_overlap_cfd)
		grid_ptr(3)%p => z(kmin_cfd:kmax_cfd)

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

        do id=1,3

        	!If there is only one CFD cell in a direction then things can be made simpler.
        	!One CFD cell means that direction is used to improve MD statistics. 
            if ( grid_sizes(2,id) - grid_sizes(1,id) == 1) then
				bbox_ptr(id)%starto = grid_sizes(1,id)
                bbox_ptr(id)%start  = grid_sizes(1,id)
				bbox_ptr(id)%endo   = grid_sizes(2,id)
                bbox_ptr(id)%end    = grid_sizes(2,id)
    
				cycle
			endif

			eps = 1.d-2 * (grid_ptr(id)%p(2) - grid_ptr(id)%p(1))
			!write(0,*) "MD: make box, grid step", myid, id, eps, grid_ptr(id)%p(1),bbox%bb(:,id)!   , grid_step !,x,y,z
			found_start = .false.

			if ( grid_ptr(id)%p(1) >= bbox%bb(1,id) - eps &
				.and. grid_ptr(id)%p(1) < bbox%bb(2,id) ) then 
				found_start = .true.
				bbox_ptr(id)%starto = idmin(id) - halo_size(1,id)
                bbox_ptr(id)%start  = idmin(id)
				!write(0,*) "MD make box l", myid,id,bbox_ptr(id)%start
			endif

			ngp = grid_sizes(2,id) - grid_sizes(1,id) + 1

			do i=2, ngp
				pl = grid_ptr(id)%p(i-1)
				pr = grid_ptr(id)%p(i) 
				!write(0,*), 'MD make bbox ', myid, id,ngp,i,pl,pr, bbox%bb(:,id) 
				if (.not. found_start )then
					if ( pl < bbox%bb(1,id)  .and. pr >=  bbox%bb(1,id) - eps ) then 
						found_start = .true.

						! here comes the decision of how much one want to cover
						bbox_ptr(id)%starto = idmin(id) + i - 1 - halo_size(1,id)
                        bbox_ptr(id)%start  = idmin(id) + i - 1
						!write(0,*), 'MD make bbox l', myid, id,i, pl, pr, bbox_ptr(id)%start
					endif
				else
					if ( (i < ngp  .and. pl <= bbox%bb(2,id) + eps  .and. pr > bbox%bb(2,id))) then
						bbox_ptr(id)%endo = idmin(id) + i - 1 - 1 + halo_size(2,id)
                        bbox_ptr(id)%end  = idmin(id) + i - 1 - 1
						!write(0,*), 'MD make bbox r', myid, id, i, pl, pr ,  bbox_ptr(id)%end		     
						exit
					else if (i == ngp  .and. abs( pr - bbox%bb(2,id)) < eps ) then 
						bbox_ptr(id)%endo = idmin(id) + ngp - 1 + halo_size(2,id)
                        bbox_ptr(id)%end  = idmin(id) + ngp - 1
						!write(0,*), 'MD make bbox r', myid, id,i, pl, pr, bbox_ptr(id)%end
					endif
				endif
			enddo

		enddo

        ! Sanity check
        ! Test if the local grid contains local MD domain
        ! that is: x(bbox%is) < bboxbb(1,1) ; bbox%bb(2,1) < x(bbox%ie) ...


		! Set maximum local grid sizes
		nlx = bbox%ieo - bbox%iso + 1
		nly = bbox%jeo - bbox%jso + 1
		nlz = bbox%keo - bbox%kso + 1

		write(0,*)' MD: bbox ', myid, bbox, nlx, nly, nlz

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
	use coupler_internal_common
	implicit none

	integer, parameter :: max_msg_length = 1024
	character(len=*),parameter :: filename = "coupler_map.log"
	character(len=*),parameter :: nl = achar(10) ! new line; portable?
	
	character(len=max_msg_length) msg, rec
	character(len=100) fmtstr
	
	integer fh, i, n, ierr

	call mpi_file_open(COUPLER_REALM_COMM,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
	call mpi_file_set_size(fh,0_MPI_OFFSET_KIND,ierr) ! discard previous data

	n = map%n
	write(rec,'(a,I5,a,I5)')  'myid', myid, ' overlaps', map%n
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
	call mpi_barrier(COUPLER_REALM_COMM,ierr)
	call mpi_file_close(fh,ierr)

end subroutine write_overlap_map


end module coupler_internal_md
