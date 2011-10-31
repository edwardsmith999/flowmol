module coupler_md_setup
	! 
	! contains data needed by CFD and MD
	! 
	use coupler_md_global_data
	implicit none
	save

	integer npx_cfd, npy_cfd, npz_cfd, nproc_cfd
	integer, allocatable :: icoord_xz_cfd(:,:) ! coordinates in XZ CFD communicator
	! FD boundary position; needed in several communication subroutines
	real(kind=kind(0.d0)) Y_boundary

	integer map_shift(3)     ! shift for the MD map to position MD slab 
	! with respect to CFD domain
 
	real(kind=kind(0.d0)) dt_CFD	      ! CFD time step

	! number of CFD boxes hold by local domain in each direction
	! the value in y direction is for the CFD grid section imersed in MD domain ( no boundaries)
	integer :: jmax_overlap_cfd  ! largest j in CFD grid overlaped by MS

	! data structures to hold CFD - MD mapping 

	type bbox_domain
	! bounding box of MD domain winthin FD global domain
	     integer is,ie,js,je,ks,ke
	     real(kind=kind(0.d0)) :: bb(2,3)
	end type bbox_domain
	type(bbox_domain), target :: bbox

	! thicknes of the MD region between the wall and CFD grid
	real(kind=kind(0.d0)) DY_PURE_MD

	! local domain lenghts, and halves
	real(kind=kind(0.d0)) domain_lengths(3), half_domain_lengths(3)
	
	type cfd_domain_map
	 integer n ! number of CFD ranks that overlap with this MD domain
	 integer, allocatable :: rank_list(:) ! rank list of overlpping CFD blocks
	 integer, allocatable :: domains(:,:) ! indices range of overlapping with CFD grid 
	end type cfd_domain_map

	type(cfd_domain_map) cfd_map

	! write or not the overlap map
	logical :: dump_ovlerlap_map = .true.

	! local grid sizes
	integer nlx, nly, nlz

contains

subroutine exchange_grid_data
	use messenger, only : npx, npy, npz
	use mpi
	implicit none

	integer i, ierr, myid, source, ngrid(3), ny_md, nblock(3), &
		iaux(12)
	  
	real(kind=kind(0.d0)) ra(2)

	! write(0,*) 'MD exchange grid data'

	! get CFD processor grid and the number of block in j direction
	call mpi_bcast( iaux, 4, MPI_INTEGER,&
		0, CFD_MD_ICOMM,ierr)
	npx_cfd = iaux(1)
	npy_cfd = iaux(2)
	npz_cfd = iaux(3)
	jmax_overlap_cfd = iaux(4)
	nproc_cfd = npx_cfd * npy_cfd * npz_cfd

	! send MD processor grid 
	call mpi_comm_rank(MD_COMM,myid,ierr)
	if ( myid == 0 ) then
		source=MPI_ROOT
	else
	 source=MPI_PROC_NULL
	endif

	call mpi_bcast((/ npx_md, npy_md, npz_md /), 3, MPI_INTEGER,&
	 source, CFD_MD_ICOMM,ierr)

	!! Test
	!	call mpi_comm_rank(MD_COMM,myid,ierr)!
	!
	!	write(0,*) 'exchange_grid_data: MD side', myid, bbox_cfd%xbb(1:2,1:npx_cfd),bbox_cfd%zbb(1:2,1:npz_cfd),&
	!	 bbox_md%xbb(1:2,1:npx_md),bbox_md%zbb(1:2,1:npz_md)    

	! CFD mesh data 
	call mpi_bcast(iaux, 12, MPI_INTEGER, 0, CFD_MD_ICOMM,ierr) 
	imino = iaux(1); imin_cfd = iaux(2);  imax_cfd = iaux(3);  imaxo = iaux(4)
	jmino = iaux(5); jmin_cfd = iaux(6);  jmax_cfd = iaux(7);  jmaxo = iaux(8)
	kmino = iaux(9); kmin_cfd = iaux(10); kmax_cfd = iaux(11); kmaxo = iaux(12)
	allocate (x(imino:imaxo),y(jmino:jmaxo),z(kmino:kmaxo))


	call mpi_bcast(x,size(x),mpi_double_precision,0,CFD_MD_ICOMM,ierr)
	call mpi_bcast(y,size(y),mpi_double_precision,0,CFD_MD_ICOMM,ierr)
	call mpi_bcast(z,size(z),mpi_double_precision,0,CFD_MD_ICOMM,ierr)
	call mpi_bcast(ra,2,mpi_double_precision,0,CFD_MD_ICOMM,ierr)

	! rescale all lengths to MD units to avoid disasters
	x = fsig * x; y = fsig * y; z = fsig * z 

	dx = fsig * ra(1); dz = fsig * ra(2)

!		write(0,*) 'MD exchange grid data: imin0, imin ...', imino, imin_cfd, imax_cfd,imaxo, &
!		 x(imino),x(imin_cfd), x(imax_cfd), x(imaxo)
!		write(0,*) 'MD exchage grid data: recv dx, dz ', dx, dz, jmax_overlap_cfd

	! get CFD nsteps
	call mpi_bcast(nsteps,1,mpi_integer,0,CFD_MD_ICOMM&
	 		&,ierr)

	! get CFD dt
	call mpi_bcast(dt_CFD,1,mpi_double_precision,0&
	 		&,CFD_MD_ICOMM,ierr)
	! should dt_CFD be scaled ?  see to it later
	dt_CFD = dt_CFD * FoP_time_ratio

	! set the sizes of MD box

	DY_PURE_MD = (y(jmin_cfd) - y(jmino)) ! prototype

	xL_md = (x(imax_cfd) - x(imin_cfd))
	yL_md = (y(jmax_overlap_cfd) - y(jmino) +&
	 & DY_PURE_MD)		   
	zL_md = (z(kmax_cfd) - z(kmin_cfd))

	!	write(0,*) 'MD: exchange_grid... xL_md, yL_md, zL_md',&
	!	 & myid, xL_md, yL_md, zL_md
	!	write(0,*) 'MD: nsteps (CFD) ', nsteps

end subroutine exchange_grid_data


subroutine create_map_cfd_md
		use mpi
		use messenger, only : icoord_md => icoord, myid, icoord, npx, npy, npz, nproc
		implicit none
		integer  i, ir, ireq(nproc), noverlaps, ierr
		integer, allocatable :: overlap_mask(:,:)


! compute the boundaries of this MD domain in the CFD global domain.
! assume that all domains have identical sides 
		
		domain_lengths(:) = (/ xL_md/npx, yL_md/npy, zL_md/npz /)
		half_domain_lengths(:) = 0.5d0 * domain_lengths(:)

! bounding boxes coordinates start from x(imin), z(kmin) and y(jmino)-DY_PURE_MD
		bbox%bb(1,:) = (icoord_md(:,myid+1)-1) * domain_lengths(:) &
		 + (/ x(imin_cfd), y(jmino)-DY_PURE_MD, z(kmin_cfd) /)
		bbox%bb(2,:) =  bbox%bb(1,:) + domain_lengths(:)

!		write(0,*) 'MD: bbox%bb ', myid, bbox%bb !, domain, npx, npy, npz, xL_md, yL_md, zL_md, icoord_md
		
		call make_bbox

!		write(0,*) 'MD: bbox%is ', myid, bbox%is, bbox%ie, bbox%js, bbox%je, bbox%ks, bbox%ke 
!
!  send box coordinate to CFD
!

		call mpi_allgather((/ bbox%is, bbox%ie, bbox%js, bbox%je, bbox%ks, bbox%ke /), 6, MPI_INTEGER,&
		 			MPI_BOTTOM,0,MPI_INTEGER,CFD_MD_ICOMM,ierr)

!!$		do i = 0, nproc_md-1
!!$
!!$		 if ( myid == i ) then
!!$		  source=MPI_ROOT
!!$		 else
!!$		  source=MPI_PROC_NULL
!!$		 endif
!!$		 
!!$		 call mpi_bcast((/ bbox%is, bbox%ie, bbox%js, bbox%je, bbox%ks, bbox%ke /), 6, MPI_INTEGER,&
!!$		  source, CFD_MD_ICOMM,ierr)
!!$
!!$		enddo


! get the domain overlap mask from cfd

		allocate(overlap_mask(0:nproc-1,0:nproc_cfd-1))


		call mpi_allgather(MPI_BOTTOM,0, MPI_INTEGER,overlap_mask, &
		 nproc,MPI_INTEGER,CFD_MD_ICOMM,ierr)

!		write(0,'(a,32I3)') 'MD, overlap mask: ', overlap_mask

		noverlaps = 0
		do i = 0, nproc_cfd - 1
		 if ( overlap_mask(myid,i) == 1) then 
		  noverlaps = noverlaps + 1
		 endif
		enddo

!		write(0,'(a,32I3)') 'MD, noverlaps: ', myid, noverlaps
		
! sort out which CFD ranks hold non-void domains for this MD rank

		cfd_map%n = noverlaps
		allocate ( cfd_map%rank_list(noverlaps), cfd_map%domains(6,noverlaps))

		ir=0
		do i=0, nproc_cfd - 1
		 
		 if (overlap_mask(myid,i) == 1) then
		  ir = ir + 1
		  cfd_map%rank_list(ir) = i
		 endif
		enddo

		do i =1, ir
		 call mpi_irecv(cfd_map%domains(1,i), 6, MPI_INTEGER,cfd_map%rank_list(i),2, CFD_MD_ICOMM,ireq(i),ierr)
		enddo
		call mpi_waitall( ir,ireq,MPI_STATUSES_IGNORE,ierr)

		 if ( dump_ovlerlap_map) then
		  call write_overlap_map
		 endif

! allocate array for CFD velocities and initialize the time indices; 
		allocate(vel_fromCFD(3,bbox%ie - bbox%is, bbox%je - bbox%js , bbox%ke - bbox%ks ,2))
		itm1 = 1; itm2 = 2
		
!		write(0,*) 'MD: end of create_map_cfd_md', myid

	contains 

		subroutine make_bbox
		 use messenger, only : icoord, myid
		 implicit none
		 
		 integer, parameter :: is = 1, ie = 2

		 type grid_pointer
		  real(kind=kind(0.d0)) , pointer :: p(:) => null()
		 end type grid_pointer

		 type bbox_pointer
		  integer, pointer :: start => null(), end => null()
		 end type bbox_pointer

		 type(grid_pointer) grid_ptr(3)
		 type(bbox_pointer) bbox_ptr(3)

		 integer id, ngp, grid_sizes(2,3), halo_size(2,3), idmin(3)
		 real(kind=kind(0.d0)) pl,pr,eps  ! left right grid points
		 logical found_start

! indices covering the CFD physical domain

		 grid_sizes(:, :) = reshape((/ imin_cfd, imax_cfd, jmin_cfd, jmax_cfd, kmin_cfd, kmax_cfd /),(/2,3/))

! starting indices (first in physical domain) in all directions

		 idmin = (/ imin_cfd, jmin_cfd, kmin_cfd /)

!		 write(0,*) 'MD: make box grid_sizes, idmin ', grid_sizes, idmin

! how large is the halo in each direction, depending also on the MD domain position

		 halo_size(:,:) = 1

! specical values for the boundaries

		if ( icoord(2,myid + 1) == 1 ) then 
		 
		 halo_size(1,2) = 0

		endif

		if ( icoord(2,myid + 1) == npy ) then
		 
		 halo_size(2,2) = 0

		endif

		if (icoord(1, myid + 1) == 1) then

		 halo_size(1,1) = 0

		endif

		if (icoord(1, myid + 1) == npx) then 

		 halo_size(2,1) = 0

		endif

		if (icoord(3, myid + 1) == 1) then 

		 halo_size(1,3) = 0

		endif

		if (icoord(3, myid + 1) == npz) then 

		 halo_size(2,3) = 0

		endif
		


! pointer to grid coordinates. Helpful to loop through dimensions
		 grid_ptr(1)%p => x(imin_cfd:imax_cfd)
		 grid_ptr(2)%p => y(jmin_cfd:jmax_cfd)
		 grid_ptr(3)%p => z(kmin_cfd:kmax_cfd)

		 bbox_ptr(1)%start => bbox%is
		 bbox_ptr(1)%end   => bbox%ie
		 bbox_ptr(2)%start => bbox%js
		 bbox_ptr(2)%end   => bbox%je
		 bbox_ptr(3)%start => bbox%ks
		 bbox_ptr(3)%end   => bbox%ke


		 do id=1,3
		 
		  eps = 1.d-2 * (grid_ptr(id)%p(2) - grid_ptr(id)%p(1))

!		  write(0,*) "MD: make box, grid step", myid, id, eps, grid_ptr(id)%p(1),bbox%bb(:,id)!   , grid_step !,x,y,z

		  found_start = .false.

		  if ( grid_ptr(id)%p(1) >= bbox%bb(1,id) - eps &
				  .and. grid_ptr(id)%p(1) < bbox%bb(2,id) ) then 
		     found_start = .true.
		     bbox_ptr(id)%start = idmin(id) - halo_size(1,id)

!		     write(0,*) "MD make box l", myid,id,bbox_ptr(id)%start
		  endif

		  ngp = grid_sizes(2,id) - grid_sizes(1,id) + 1

		  do i=2, ngp
		   
		   pl = grid_ptr(id)%p(i-1)
		   pr = grid_ptr(id)%p(i) 

!		    write(0,*), 'MD make bbox ', myid, id,ngp,i,pl,pr, bbox%bb(:,id) 
		   
		   if (.not. found_start )then
		    if ( pl < bbox%bb(1,id)  .and. pr >=  bbox%bb(1,id) - eps ) then 
		     found_start = .true.
		    
		     
! here comes the decision of how much one want to cover
		     
		      bbox_ptr(id)%start = idmin(id) + i - 1 - halo_size(1,id)
		      
 !		     write(0,*), 'MD make bbox l', myid, id,i, pl, pr, bbox_ptr(id)%start
		     endif
		      
		   else
		   
		   
		    
		    if ( (i < ngp  .and. pl <= bbox%bb(2,id) + eps  .and. pr > bbox%bb(2,id))) then

		     bbox_ptr(id)%end = idmin(id) + i - 1 -1 + halo_size(2,id)

!		    write(0,*), 'MD make bbox r', myid, id, i, pl, pr ,  bbox_ptr(id)%end		     
		     exit
		    
		    else if (i == ngp  .and. abs( pr - bbox%bb(2,id)) < eps ) then 
		     
		     bbox_ptr(id)%end = idmin(id) + ngp - 1 + halo_size(2,id)
		     
 !		    write(0,*), 'MD make bbox r', myid, id,i, pl, pr, bbox_ptr(id)%end

		   endif
		   
		  endif

		  enddo
		 enddo

! Set local grid sizes
		 nlx = bbox%ie - bbox%is + 1
		 nly = bbox%je - bbox%js + 1
		 nlz = bbox%ke - bbox%ks + 1

		end subroutine make_bbox



!!$		use mpi
!!$		use messenger, only : icomm_grid_md => icomm_grid, icoord_md => icoord
!!$		use computational_constants_MD,    only : nproc_md => nproc, domain
!!$		implicit none
!!$
!!$		integer i, j, k, ip, ierr, myid_cfd, myid_md, myid_xz_cfd, iroot,&
!!$		 myid_world, source, istat
!!$		integer ibmin_loc, ibmax_loc, kbmin_loc, kbmax_loc, ibs, ibe, &
!!$		 kbs, kbe, four_dp_type
!!$		integer a_sizes(3), a_subsizes(3), a_starts(3)
!!$		integer, allocatable :: overlap_aux(:,:)
!!$! ATTENTION Assumes at the moment that npy_md=1
!!$
!!$! get the rank of jblock 1
!!$
!!$		call mpi_bcast(iroot, 1, mpi_integer, 0, CFD_MD_ICOMM, ierr)
!!$
!!$! Send to CFD the smalest MPI_COMM_WORLD rank to avoid any implicit assumption that
!!$! might be broken e.g. the MD ranks start after CFD ranks in MPI_COMM_WORLD
!!$
!!$		call mpi_comm_rank(MD_COMM, myid_md, ierr)
!!$		call mpi_comm_rank(MPI_COMM_WORLD, myid_world, ierr)
!!$		if (myid_md == 0) then
!!$		 source = MPI_ROOT
!!$		else
!!$		 source = MPI_PROC_NULL
!!$		endif
!!$
!!$		call mpi_bcast(myid_world, 1, MPI_INTEGER, source,CFD_MD_ICOMM, ierr)
!!$
!!$!  The MD part should be restricted to the processors that 
!!$!  contain the CFD boundary; I'll see later to it
!!$
!!$		call mpi_intercomm_create(icomm_grid_md, 0, MPI_COMM_WORLD, iroot, 2,&
!!$		 CFD_BDY_MD_ICOMM, ierr)
!!$
!!$		allocate(icoord_xz_cfd(2,npx_cfd*npz_cfd))
!!$		call mpi_bcast(icoord_xz_cfd,2*npx_cfd*npz_cfd,mpi_integer,0,&
!!$		 CFD_BDY_MD_ICOMM,ierr)
!!$
!!$
!!$		call mpi_comm_rank(MD_COMM,myid_md,ierr)
!!$
!!$		if( myid_md == 0 ) then
!!$		 source = MPI_ROOT
!!$		else
!!$		 source = MPI_PROC_NULL
!!$		endif
!!$
!!$		write(0,*) 'MD side ', myid_md, ' icoord_xz_cfd ', icoord_xz_cfd
!!$
!!$		call mpi_bcast(icoord_md, 3*nproc_md, mpi_integer, source, &
!!$		 CFD_BDY_MD_ICOMM , ierr)
!!$		allocate(overlap_aux(9,npx_cfd*npz_cfd))
!!$		overlap_aux(:,:) = -1
!!$
!!$!   three doubles for the speed values
!!$!		call mpi_type_contiguous(4, MPI_DOUBLE_PRECISION, four_dp_type,ierr)
!!$!		call mpi_type_commit(four_dp_type, ierr)
!!$		four_dp_type = MPI_DOUBLE_PRECISION
!!$
!!$		ibmin_loc = bbox_md%xbb(1,icoord_md(1,myid_md+1))
!!$		ibmax_loc = bbox_md%xbb(2,icoord_md(1,myid_md+1))
!!$		kbmin_loc = bbox_md%zbb(1,icoord_md(3,myid_md+1))
!!$		kbmax_loc = bbox_md%zbb(2,icoord_md(3,myid_md+1))
!!$
!!$		novr = 0
!!$		do j= 1, npx_cfd*npz_cfd
!!$		 i = icoord_xz_cfd(1,j)
!!$		 k = icoord_xz_cfd(2,j)
!!$		 
!!$		 ibs = bbox_cfd%xbb(1,i)
!!$		 ibe = bbox_cfd%xbb(2,i)
!!$		 kbs = bbox_cfd%zbb(1,k)
!!$		 kbe = bbox_cfd%zbb(2,k)
!!$		 
!!$
!!$		 if( ((    ibs <  ibmin_loc  .and. ibe > ibmin_loc ) .or. &
!!$				 ( ibs >= ibmin_loc  .and. ibs < ibmax_loc )) .and. &
!!$				 ((kbs <  kbmin_loc  .and. kbe > kbmin_loc ) .or. &
!!$				 ( kbs >= kbmin_loc  .and. kbs < kbmax_loc )) ) then
!!$
!!$				 novr = novr+1
!!$
!!$				 overlap_aux(1,novr) = j
!!$				 overlap_aux(2,novr) = i
!!$				 overlap_aux(3,novr) = k
!!$				 overlap_aux(4,novr) = max(ibmin_loc,ibs)
!!$				 overlap_aux(5,novr) = min(ibmax_loc,ibe)
!!$				 overlap_aux(6,novr) = max(kbmin_loc,kbs)
!!$				 overlap_aux(7,novr) = min(kbmax_loc,kbe)
!!$! Creat subarray types for data transfers
!!$! MD needs  subsizes correspondig to the overlaping FD domains
!!$! at the moment this is for four dimensional vectors on a three dimensional grid
!!$
!!$				 a_sizes(:) = (/ kbmax_loc - kbmin_loc, ibmax_loc - ibmin_loc, 1 /) ! follows index order of uc, vc, wc
!!$				 a_subsizes(:) = (/ overlap_aux(7,novr) - overlap_aux(6,novr), &
!!$				  overlap_aux(5,novr) - overlap_aux(4,novr), 1/)
!!$				 a_starts(:) = (/ overlap_aux(6,novr) - kbmin_loc, &
!!$				  overlap_aux(4,novr) - ibmin_loc, 0 /)
!!$
!!$				  write(0,*) 'MD side setup subarray :', myid_md, a_sizes, a_subsizes, a_starts
!!$
!!$				 call mpi_type_create_subarray(3, a_sizes, a_subsizes, &
!!$				  a_starts, MPI_ORDER_FORTRAN, four_dp_type, overlap_aux(8,novr), ierr)
!!$				 call mpi_type_commit(overlap_aux(8,novr), ierr)
!!$
!!$! now for array that hold grid ( plabes ) points. They are longer with one and the plane between teo processors
!!$! must go to both of them (typical situation is CFD -> MD as MD domain are smaller)				  
!!$				 a_sizes(:) = (/ kbmax_loc - kbmin_loc+1, ibmax_loc - ibmin_loc+1, 1 /) ! follows index order of uc, vc, wc
!!$				 a_subsizes(:) = (/ overlap_aux(7,novr) - overlap_aux(6,novr) + 1, &
!!$				  overlap_aux(5,novr) - overlap_aux(4,novr) + 1, 1/)
!!$				 a_starts(:) = (/ overlap_aux(6,novr) - kbmin_loc, &
!!$				  overlap_aux(4,novr) - ibmin_loc, 0 /)
!!$				 call mpi_type_create_subarray(3, a_sizes, a_subsizes, &
!!$				  a_starts, MPI_ORDER_FORTRAN, four_dp_type, overlap_aux(9,novr), ierr)
!!$				 call mpi_type_commit(overlap_aux(9,novr), ierr)
!!$
!!$		 endif
!!$		enddo
!!$
!!$		write(0,*) 'MD side ', myid_md, ' novr ', novr
!!$! copy everthing in the map_overlap
!!$		allocate(map_overlap(novr),stat = istat)
!!$
!!$		do i=1,novr
!!$		 map_overlap(i)%rank = overlap_aux(1,i)
!!$		 map_overlap(i)%coord(:) = overlap_aux(2:3,i)
!!$		 map_overlap(i)%ib_range(:) = overlap_aux(4:5,i)
!!$		 map_overlap(i)%kb_range(:) = overlap_aux(6:7,i)
!!$		 map_overlap(i)%dp2d_type = overlap_aux(8,i)
!!$		 map_overlap(i)%plane_type = overlap_aux(9,i)
!!$		enddo
!!$
!!$		deallocate(overlap_aux)
!!$
!!$! Some other needed initialisation before the computation starts
!!$! number of boxes for the common CFD, MD region
!!$		nib = bbox_md%xbb(2,icoord_md(1,myid_md+1)) - bbox_md%xbb(1,icoord_md(1,myid_md+1))
!!$		nkb = bbox_md%zbb(2,icoord_md(3,myid_md+1)) - bbox_md%zbb(1,icoord_md(3,myid_md+1))
!!$! allocate array for CFD velocities and initialize the time indices; includes halos in x,z directions
!!$		allocate(vel_fromCFD(3,0:nib+1,njb,0:nkb+1,2))
!!$		itm1 = 1; itm2 = 2
!!$
!!$! where CFD bottom is positioned in MD box 
!!$! this corresponds to y(1) in CFD
!!$		Y_boundary = 2.0 * fsig * (y(1) - y(0))
!!$
!!$		write(200+myid_md,*) 'id ', myid_cfd, ' novr ', novr
!!$		write(200+myid_md,*) 'ibmin/max, kbmin/max ',ibmin_loc, ibmax_loc, kbmin_loc, kbmax_loc 
!!$		write(200+myid_md,*) 'icoord_md ' , icoord_md
!!$		write(200+myid_md,'(8I5)') map_overlap

	end subroutine create_map_cfd_md


	subroutine write_overlap_map
		use mpi
		use messenger, only : myid
		implicit none

		integer, parameter :: max_msg_length = 1024
		character(len=*),parameter :: filename = "coupler_map.log"
		character(len=*),parameter :: nl = achar(10) ! new line; portable?
		
		character(len=max_msg_length) msg, rec
		character(len=100) fmtstr
		
		integer fh, i, n, ierr

		call mpi_file_open(MD_COMM,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)

		n = cfd_map%n

		write(rec,'(a,I5,a,I5)')  'myid', myid, ' overlaps', cfd_map%n

		msg = rec(1:len_trim(rec))//nl

		write(fmtstr,*) "( a,",1,"I5,","a,", 6,"I5)"
		
		do i = 1, n
		 write(rec,fmtstr) 'rank', cfd_map%rank_list(i), &
		  ' indicies',cfd_map%domains(1:6,i)
		 
		 msg = msg(1:len_trim(msg))//rec(1:len_trim(rec))//nl
		enddo
		
		msg = msg(1:len_trim(msg))//nl
		n = len_trim(msg)

		call mpi_file_write_shared(fh,msg(1:n),n,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr) 
		
		call mpi_barrier(MD_COMM,ierr)

		call mpi_file_close(fh,ierr)

	end subroutine write_overlap_map


end module coupler_md_setup
