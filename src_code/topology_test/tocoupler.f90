

!==============================================================================

module coupler

	! Overloading
    interface CPL_send
        module procedure CPL_send_3d, CPL_send_4d
    end interface

    interface CPL_recv
        module procedure CPL_recv_3d, CPL_recv_4d
    end interface

    private CPL_send_3d, CPL_send_4d, &
        CPL_send_xd, CPL_recv_3d, CPL_recv_4d,&
        CPL_recv_xd

contains


!================================
!
!	COUPLER SETUP
!
!================================


subroutine get_cell_ranges_md
	use coupler_module
	implicit none

	integer :: n
	integer :: olap_jmin_mdcoord
	integer :: ncxl, ncyl, nczl
	integer :: ncy_mdonly, ncy_md, ncyP_md
	integer, parameter :: cnull = -666 

	icPmin_md = cnull
	jcPmin_md = cnull
	kcPmin_md = cnull
	icPmax_md = cnull
	jcPmax_md = cnull
	kcPmax_md = cnull
	
	ncxl = ncx / npx_md
	do n=1,npx_md
		icPmax_md(n) = n * ncxl
		icPmin_md(n) = icPmax_md(n) - ncxl + 1
	end do	

	nczl = ncz / npz_md
	do n=1,npz_md
		kcPmax_md(n) = n * nczl
		kcPmin_md(n) = kcPmax_md(n) - nczl + 1
	end do

	ncy_md   = nint(yL_md/dy)
	ncy_mdonly = ncy_md - ncy_olap
	ncyP_md = ncy_md / npy_md
	olap_jmin_mdcoord = npy_md - floor(dble(ncy_olap)/dble(ncyP_md))	 
	do n = olap_jmin_mdcoord,npy_md
		jcPmax_md(n) = n * ncyP_md - ncy_mdonly
		jcPmin_md(n) = jcPmax_md(n) - ncyP_md + 1
		if (jcPmin_md(n).le.0) jcPmin_md(n) = 1
	end do 

	if (myid_world.eq.0) call print_map_md

contains
	subroutine print_map_md
		use coupler_module
		implicit none

		integer :: n

		write(6000+myid_world,*), ''
		write(6000+myid_world,*), '==========================================='
		write(6000+myid_world,*), '------------ M D   M A P ------------------'
		write(6000+myid_world,*), '==========================================='
		write(6000+myid_world,*), 'npx_md = ', npx_md
		write(6000+myid_world,*), 'ncx    = ', ncx
		write(6000+myid_world,*), 'ncxl   = ', ncxl
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), '  icoord_md   icPmin_md(n)    icPmax_md(n) '
		write(6000+myid_world,*), '-------------------------------------------'
		do n=1,npx_md
			write(6000+myid_world,'(1x,3i11)'), n, icPmin_md(n), icPmax_md(n)
		end do	
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), 'npy_md     = ', npy_md
		write(6000+myid_world,*), 'ncy_md     = ', ncy_md
		write(6000+myid_world,*), 'ncyP_md    = ', ncyP_md 
		write(6000+myid_world,*), 'ncy_olap   = ', ncy_olap
		write(6000+myid_world,*), 'ncy_mdonly = ', ncy_mdonly
		write(6000+myid_world,*), 'olap_jmin_mdcoord = ', olap_jmin_mdcoord
		write(6000+myid_world,*), 'dy         = ', dy
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), '  jcoord_md   jcPmin_md(n)    jcPmax_md(n) '
		write(6000+myid_world,*), '-------------------------------------------'
		do n = 1,npy_md	
			write(6000+myid_world,'(1x,3i11)'), n, jcPmin_md(n), jcPmax_md(n)
		end do
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), 'npz_md = ', npz_md
		write(6000+myid_world,*), 'ncz    = ', ncz
		write(6000+myid_world,*), 'nczl   = ', nczl
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), '  kcoord_md   kcPmin_md(n)    kcPmax_md(n) '
		write(6000+myid_world,*), '-------------------------------------------'
		do n=1,npz_md
			write(6000+myid_world,'(1x,3i11)'), n, kcPmin_md(n), kcPmax_md(n)
		end do
		write(6000+myid_world,*), '-------------------------------------------'

	end subroutine print_map_md
	
end subroutine get_cell_ranges_md



subroutine get_overlap_blocks_cfd
	use coupler_module
	implicit none

	integer :: n,i,endproc
	integer :: nolapsx, nolapsy, nolapsz, nolaps
	integer, parameter :: olap_null = -666

	cfd_icoord2olap_md_icoords = olap_null
	cfd_jcoord2olap_md_jcoords = olap_null
	cfd_kcoord2olap_md_kcoords = olap_null
	
	nolapsx = npx_md/npx_cfd	
	do n = 1, npx_cfd
		do i = 1,nolapsx	
			cfd_icoord2olap_md_icoords(n,i) = (n-1)*nolapsx + i
		end do
	end do
	
	nolapsz = npz_md/npz_cfd	
	do n = 1, npz_cfd
		do i = 1,nolapsz	
			cfd_kcoord2olap_md_kcoords(n,i) = (n-1)*nolapsz + i
		end do
	end do

	nolapsy = ceiling(yL_olap/yLl_md)
	endproc = ceiling(yL_olap/yLl_cfd)
	do n = 1,endproc
		do i = 1,nolapsy
			cfd_jcoord2olap_md_jcoords(n,i) = (n-1)*nolapsy + i &
			                                  + (npy_md - nolapsy)
		end do
	end do
	
	nolaps = nolapsx*nolapsy*nolapsz	
	if(myid_world.eq.0) call print_map_cfd

contains

	subroutine print_map_cfd
		use coupler_module
		implicit none

		integer :: n,i
		
		write(6000+myid_world,*), ''
		write(6000+myid_world,*), '==========================================='
		write(6000+myid_world,*), '------------ C F D   M A P ----------------'
		write(6000+myid_world,*), '==========================================='
		write(6000+myid_world,*), 'npx_cfd = ', npx_cfd
		write(6000+myid_world,*), 'nolapsx = ', nolapsx
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), '  icoord_cfd       olapmin     olapmax     ' 
		write(6000+myid_world,*), '-------------------------------------------'
		do n=1,npx_cfd
			write(6000+myid_world,'(1x,3i11)'), n,               &
				  cfd_icoord2olap_md_icoords(n,1),         &
				  cfd_icoord2olap_md_icoords(n,nolapsx)
		end do	
		write(6000+myid_world,*), '-------------------------------------------'

		write(6000+myid_world,*), 'npy_cfd = ', npy_cfd
		write(6000+myid_world,*), 'nolapsy = ', nolapsy
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), '  jcoord_cfd       olapmin     olapmax     ' 
		write(6000+myid_world,*), '-------------------------------------------'
		do n=1,npy_cfd
			write(6000+myid_world,'(1x,3i11)'), n,               &
				  cfd_jcoord2olap_md_jcoords(n,1),         &
				  cfd_jcoord2olap_md_jcoords(n,nolapsy)
		end do	
		write(6000+myid_world,*), '-------------------------------------------'

		write(6000+myid_world,*), 'npz_cfd = ', npz_cfd
		write(6000+myid_world,*), 'nolapsz = ', nolapsz
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), '  kcoord_cfd       olapmin     olapmax     ' 
		write(6000+myid_world,*), '-------------------------------------------'
		do n=1,npz_cfd
			write(6000+myid_world,'(1x,3i11)'), n,               &
				  cfd_kcoord2olap_md_kcoords(n,1),         &
				  cfd_kcoord2olap_md_kcoords(n,nolapsz)
		end do	
		write(6000+myid_world,*), '-------------------------------------------'

	end subroutine print_map_cfd
	
end subroutine get_overlap_blocks_cfd





!=========================================================================

subroutine prepare_overlap_comms
	use coupler_module
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

	allocate(mdicoords(npx_md/npx_cfd))
	allocate(mdjcoords(npy_md/npy_cfd))
	allocate(mdkcoords(npz_md/npz_cfd))

	!Set default values, must be done because coord2rank_md cannot
	!take "null" coordinates.
	group(:) = olap_null
	olap_mask(:) = 0
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
			olap_mask(trank_world) = 1
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

			olap_mask(trank_world) = 1
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
	if (myid_olap .eq. CFDid_olap) testval = group(rank_world)
	call MPI_bcast(testval,1,MPI_INTEGER,CFDid_olap,CPL_OLAP_COMM,ierr)

	! Set all non-overlapping processors to MPI_COMM_NULL
	if (olap_mask(rank_world).eq.0) then
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
	
	if (realm.eq.md_realm) call write_overlap_comms_md

end subroutine prepare_overlap_comms


!=========================================================================
!Setup topology graph of overlaps between CFD & MD processors

subroutine CPL_overlap_topology
	use coupler_module
	use mpi
	implicit none

	integer								:: i, n, nneighbors, nconnections
	integer, dimension(:),allocatable	:: index, edges, id_neighbors
	logical								:: reorder

	!Allow optimisations of ordering
	reorder = .true.

	!Get number of processors in communicating overlap region 
	if (olap_mask(rank_world).eq.1) then

		!CFD processor is root and has mapping to all MD processors
		allocate(index(nproc_olap))			!Index for each processor
		allocate(edges(2*(nproc_olap)-1))	!nproc_olap-1 for CFD and one for each of nproc_olap MD processors
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

        ! <><><><><><>  TEST <><><><><><>  TEST <><><><><><>  TEST <><><><><><> 
        !Get number of neighbours
        call MPI_comm_rank( CPL_GRAPH_COMM, myid_graph, ierr)
        call MPI_Graph_neighbors_count( CPL_GRAPH_COMM, myid_graph, nneighbors, ierr)
        allocate(id_neighbors(nneighbors))
        !Get neighbours
        call MPI_Graph_neighbors( CPL_GRAPH_COMM, myid_graph, nneighbors, id_neighbors,ierr )
        select case(realm)
        case(cfd_realm)
                write(3000+myid_world,*), realm_name(realm),' My graph', & 
								myid_world,myid_graph,myid_olap, & 
								rank2coord_cfd(:,rank_realm), nneighbors, id_neighbors
        case(md_realm)
                write(3000+myid_world,*), realm_name(realm),' My graph', & 
								myid_world,myid_graph,myid_olap, & 
								rank2coord_md(:,rank_realm), nneighbors, id_neighbors
        end select
		! <><><><><><>  TEST <><><><><><>  TEST <><><><><><>  TEST <><><><><><> 

	else
		CPL_GRAPH_COMM = MPI_COMM_NULL
	endif

	! Setup graph map
	call CPL_rank_map(CPL_GRAPH_COMM,rank_graph,nproc_olap, & 
					 rank_graph2rank_world,rank_world2rank_graph,ierr)
	myid_graph = rank_graph - 1

end subroutine CPL_overlap_topology





!================================
!
!	COUPLER MAIN
!
!================================




!==============================================================================

!------------------------------------------------------------------------------
!                              CPL_gather                                     -
!------------------------------------------------------------------------------

! Perform gather operation on CPL_OLAP_COMM communicator. The CFD processor
! is the root process.  

! - - - Synopsis - - -

! CPL_gather(gatherarray,npercell)

! - - - Input Parameters - - -

!gatherarray
!	assumed shape array of data to be gathered (double precision)
!npercell
!	number of data points per cell to be gathered (integer)
!	note - should be the same as size(gatherarray(1))

! - - - Output Parameters - - -
! - NONE -

subroutine CPL_gather(gatherarray,npercell)
	use mpi
	use coupler_module
	implicit none

	integer, intent(in) :: npercell
	real(kind(0.d0)), dimension(:,:,:,:), intent(in) :: gatherarray
	real(kind(0.d0)), dimension(:), allocatable :: sendbuf 

	integer :: sendcount
	integer, dimension(:), allocatable :: recvcounts, displs
	real(kind(0.d0)), dimension(:), allocatable :: recvbuf 
		
	call prepare_gatherv_parameters	
	
	if (realm.eq.md_realm) call pack_sendbuf
	!if (realm.eq.md_realm) call CPL_pack(gatherarray,sendbuf)

	call MPI_gatherv(sendbuf,sendcount,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,&
	                 displs,MPI_DOUBLE_PRECISION,CFDid_olap,CPL_OLAP_COMM,ierr)

	if (realm.eq.cfd_realm) call unpack_recvbuf 
	!if (realm.eq.cfd_realm) call CPL_unpack(gatherbuf,recvarray)

	call deallocate_gather_u
	
contains

	subroutine prepare_gatherv_parameters
		implicit none

		integer :: coord(3), extents(6)
		integer :: ncells
		integer :: trank_olap,trank_world,trank_cart,tid_olap
		character(len=128) :: errmessage

		! Check if CFD processor has tried to "send" anything
		if (myid_olap.eq.CFDid_olap .and. any(shape(gatherarray).ne.0)) then
			call error_abort('CFD proc input to CPL_gather: '          // &
							 'gatherarray has nonzero size. Aborting ' // &
							 'from prepare_gatherv_parameters' )
		end if

		! Allocate send buffer
		allocate(sendbuf(size(gatherarray)))

		! Allocate array of sendbuffer sizes and populate it
		allocate(recvcounts(nproc_olap))
		do trank_olap = 1,nproc_olap
			tid_olap = trank_olap - 1
			if (tid_olap .eq. CFDid_olap) then
				recvcounts(trank_olap) = 0 
			else
				call CPL_Cart_coords(CPL_OLAP_COMM,trank_olap,md_realm,3, &
				                     coord,ierr) 
				call CPL_proc_extents(coord,md_realm,extents,ncells)
				recvcounts(trank_olap) = npercell*ncells 
			end if
		end do
	
		! Grab own sendbuffer size
		sendcount = size(gatherarray)
		! Sanity check
		if (sendcount .ne. recvcounts(rank_olap)) then
			call error_abort('Send buffer sizes calculated incorrectly '//&
			                 'in prepare_gatherv_parameters. Aborting.')
		end if	

		! Allocate recvbuffer on CFD proc
		allocate(recvbuf(sum(recvcounts)))

		! Calculate displacements for each proc in array recvbuffer
		allocate(displs(nproc_olap))
		displs(1) = 0
		do trank_olap=2,nproc_olap
			displs(trank_olap) = sum(recvcounts(1:trank_olap-1))	
		end do

	end subroutine prepare_gatherv_parameters

	subroutine pack_sendbuf
		implicit none

		integer :: coord(3), extents(6)
		integer :: pos, ixyz, icell, jcell, kcell

		call CPL_cart_coords(CPL_OLAP_COMM,rank_olap,md_realm,3,coord,ierr)
		call CPL_proc_extents(coord,md_realm,extents)

		pos = 1
		do ixyz  = 1,size(gatherarray,1)
		do icell = 1,size(gatherarray,2)
		do jcell = 1,size(gatherarray,3)
		do kcell = 1,size(gatherarray,4)
			sendbuf(pos) = gatherarray(ixyz,icell,jcell,kcell)
			pos = pos + 1
		end do
		end do
		end do
		end do
	
	end subroutine pack_sendbuf
	
	subroutine unpack_recvbuf
		implicit none

		integer :: coord(3), extents(6)
		integer :: trank_olap, trank_world, trank_cart, tid_olap
		integer :: pos,ixyz,icell,jcell,kcell
		real(kind(0.d0)), dimension(:,:,:,:), allocatable :: recvarray 

		!coord(:) = rank2coord_cfd(:,rank_cart)
		! Get CFD proc coords and extents, allocate suitable array
		call CPL_cart_coords(CPL_OLAP_COMM,rank_olap,cfd_realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,extents)
		allocate(recvarray(3,extents(1):extents(2), &
		                     extents(3):extents(4), &
		                     extents(5):extents(6)))

		! Loop over all processors in overlap comm
		do trank_olap = 1,nproc_olap

			tid_olap = trank_olap - 1
			if (tid_olap .eq. CFDid_olap) cycle

			call CPL_Cart_coords(CPL_OLAP_COMM,trank_olap,md_realm,3,coord,ierr)
			call CPL_proc_extents(coord,md_realm,extents)

			! Set position and unpack MD proc's part of recvbuf to
			! correct region of recvarray	
			pos = displs(trank_olap) + 1	
			write(8000+myid_olap,'(a)'), 'recvarray(ixyz,icell,jcell,kcell)'
			do ixyz = 1,3
			do icell = extents(1),extents(2)
			do jcell = extents(3),extents(4)
			do kcell = extents(5),extents(6)

				recvarray(ixyz,icell,jcell,kcell) = recvbuf(pos)
				pos = pos + 1

				write(8000+myid_olap,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
				      'recvarray(',ixyz,',',icell,',',jcell,',',kcell,') =', &
				       recvarray(ixyz,icell,jcell,kcell)

			end do	
			end do	
			write(8000+myid_olap,'(a)'), 'recvarray(ixyz,icell,jcell,kcell)'
			end do
			end do
					
		end do

		deallocate(recvarray)

	end subroutine unpack_recvbuf
	
	subroutine deallocate_gather_u
		implicit none

		if(allocated(recvcounts)) deallocate(recvcounts)
		if(allocated(displs))     deallocate(displs)
		if(allocated(sendbuf))    deallocate(sendbuf)
		if(allocated(recvbuf))    deallocate(recvbuf)

	end subroutine deallocate_gather_u

end subroutine CPL_gather

!------------------------------------------------------------------------------
!                              CPL_gather                                     -
!------------------------------------------------------------------------------

! Scatter cell-wise data from CFD processor to corresponding MD processors
! on the overlap communicator CPL_OLAP_COMM.

! - - - Synopsis - - -

! CPL_scatter(scatterarray,npercell)

! - - - Input Parameters - - -

!scatterarray
!	assumed shape array of data to be scattered (double precision)
!npercell
!	number of data points per cell to be gathered (integer)
!	note - should be the same as size(gatherarray(1))

! - - - Output Parameters - - -
! - NONE -

subroutine CPL_scatter(scatterarray,npercell)
	use coupler_module
	use mpi
	implicit none

	integer, intent(in) :: npercell
	real(kind(0.d0)), dimension(:,:,:,:), intent(in) :: scatterarray

	integer :: recvcount
	integer, dimension(:), allocatable :: displs,sendcounts
	real(kind(0.d0)), dimension(:), allocatable :: recvbuf 
	real(kind(0.d0)), dimension(:), allocatable :: scatterbuf

	call prepare_scatterv_parameters

	if (realm.eq.cfd_realm) call pack_scatterbuf

	call MPI_scatterv(scatterbuf,sendcounts,displs,MPI_DOUBLE_PRECISION, &
	                  recvbuf,recvcount,MPI_DOUBLE_PRECISION,CFDid_olap, &
	                  CPL_OLAP_COMM,ierr) 

	if (realm.eq.md_realm)  call unpack_scatterbuf

	call deallocate_scatter_s

contains

	subroutine prepare_scatterv_parameters
		implicit none

		integer :: ncxl,ncyl,nczl
		integer :: ncells
		integer :: icell,jcell,kcell
		integer :: coord(3),extents(6)
		integer :: bufsize
		integer :: trank_olap, tid_olap, trank_world, trank_cart

		! THIS ASSUMES ONE OVERLAP CFD PROCESSOR IN Y DIRECTION
		ncxl = icmax_olap - icmin_olap + 1
		ncyl = jcmax_olap - jcmin_olap + 1
		nczl = kcmax_olap - kcmin_olap + 1

		if (realm.eq.cfd_realm) bufsize = npercell*ncxl*ncyl*nczl
		if (realm.eq.md_realm)  bufsize = 0

		allocate(scatterbuf(bufsize))
		allocate(sendcounts(nproc_olap))
		allocate(displs(nproc_olap))

		! Loop over all procs in overlap comm
		do trank_olap = 1,nproc_olap
			tid_olap = trank_olap - 1
			! Calculate number of data points to scatter to each proc
			if (tid_olap .eq. CFDid_olap) then
				sendcounts(trank_olap) = 0 
			else
				call CPL_Cart_coords(CPL_OLAP_COMM,trank_olap,md_realm,3, &
				                     coord, ierr) 
				call CPL_proc_extents(coord,md_realm,extents,ncells)
				sendcounts(trank_olap) = npercell*ncells 
			end if
		end do

		! Get number of data points this MD proc will receive
		recvcount = sendcounts(rank_olap)

		! Calculate starting positions of each MD proc region in
		! scatterbuf array
		displs(1) = 0
		do trank_olap=2,nproc_olap
			displs(trank_olap) = sum(sendcounts(1:trank_olap-1))	
		end do

		! Allocate space to receive data
		allocate(recvbuf(sum(sendcounts)))

	end subroutine prepare_scatterv_parameters

	subroutine pack_scatterbuf
		implicit none
		
		integer :: pos, n
		integer :: trank_world, trank_cart
		integer :: coord(3), extents(6)
		integer :: ixyz, icell, jcell, kcell

		! CFD proc is rank 1, loop over MD procs in olap comm and
		! pack scatter buffer in separate regions for each MD proc
		pos = 1
		do n = 2,nproc_olap

			call CPL_Cart_coords(CPL_OLAP_COMM,n,md_realm,3,coord,ierr)
			call CPL_proc_extents(coord,md_realm,extents)

			do ixyz = 1,npercell
			do icell= extents(1),extents(2)
			do jcell= extents(3),extents(4)
			do kcell= extents(5),extents(6)
				scatterbuf(pos) = 0.1d0*ixyz + 1*icell + &
				                            1000*jcell + &
				                         1000000*kcell
				pos = pos + 1
			end do
			end do
			end do
			end do

		end do
	
	end subroutine pack_scatterbuf	

	subroutine unpack_scatterbuf
		implicit none

		integer :: pos, n, ierr
		integer :: trank_world, trank_cart
		integer :: coord(3), extents(6)
		integer :: ixyz, icell, jcell, kcell
		real(kind(0.d0)), dimension(:,:,:,:), allocatable :: recvarray

		call CPL_cart_coords(CPL_OLAP_COMM,rank_olap,md_realm,3,coord,ierr)
		call CPL_proc_extents(coord,realm,extents)

		allocate(recvarray(npercell,extents(1):extents(2), &
		                            extents(3):extents(4), &
		                            extents(5):extents(6)))

		write(7000+myid_olap,'(a)'), 'recvarray(ixyz,icell,jcell,kcell)'
		pos = 1
		do ixyz  = 1,npercell
		do icell= extents(1),extents(2)
		do jcell= extents(3),extents(4)
		do kcell= extents(5),extents(6)
			recvarray(ixyz,icell,jcell,kcell) = recvbuf(pos)
			write(7000+myid_realm,'(i4,a,i4,a,i4,a,i4,a,i4,a,f20.1)'),        &
				  rank_cart,' recvarray(',ixyz,',',icell,',',jcell,',',kcell, &
				  ') =',recvarray(ixyz,icell,jcell,kcell)
			pos = pos + 1
		end do	
		end do	
		write(7000+myid_olap,'(a)'), 'recvarray(ixyz,icell,jcell,kcell)'
		end do
		end do

		if(allocated(recvarray))  deallocate(recvarray)

	end subroutine unpack_scatterbuf
	
	subroutine deallocate_scatter_s
		implicit none

		if(allocated(displs))     deallocate(displs)
		if(allocated(scatterbuf)) deallocate(scatterbuf)
		if(allocated(recvbuf))    deallocate(recvbuf)
		if(allocated(sendcounts)) deallocate(sendcounts)

	end subroutine deallocate_scatter_s

end subroutine CPL_scatter


!=============================================================================
! CPL_send_data wrapper for 3d arrays
! see CPL_send_xd for input description
!-----------------------------------------------------------------------------
subroutine CPL_send_3d(temp,jcmax_send,jcmin_send,index_transpose)
    implicit none
 
 	integer, intent(in), optional						:: jcmax_send,jcmin_send

    integer, optional, intent(in) :: index_transpose(3)    
    real(kind=kind(0.d0)),dimension(:,:,:), intent(in) :: temp
    
    integer n1,n2,n3,n4
    real(kind=kind(0.d0)),dimension(:,:,:,:),allocatable :: asend
   
    n1 = 1
    n2 = size(temp,1)
    n3 = size(temp,2)
    n4 = size(temp,3)

	!Add padding column to 3D array to make it 4D
	allocate(asend(n1,n2,n3,n4))
	asend(1,:,:,:) = temp(:,:,:)
  
    call CPL_send_xd(asend,jcmax_send,jcmin_send,index_transpose)

end subroutine CPL_send_3d

!=============================================================================
! CPL_send_data wrapper for 4d arrays
! see CPL_send_xd for input description
!-----------------------------------------------------------------------------
subroutine CPL_send_4d(asend,jcmax_send,jcmin_send,index_transpose)
    implicit none
 
 	integer, intent(in), optional						:: jcmax_send,jcmin_send
	integer, optional, intent(in)  :: index_transpose(3)   
	real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in) :: asend
    
    integer n1,n2,n3,n4
 
    call CPL_send_xd(asend,jcmax_send,jcmin_send,index_transpose)

end subroutine CPL_send_4d

!=============================================================================
! Send data from the local grid to the associated ranks from the other 
! realm
!-----------------------------------------------------------------------------
subroutine CPL_send_xd(asend,jcmax_send,jcmin_send,index_transpose)
	use mpi
	use coupler_module
	implicit none

    ! Minimum and maximum values of j to send
 	integer, intent(in), optional						:: jcmax_send,jcmin_send
    ! Specify the order of dimensions in asend default is ( x, y, z) but some CFD solvers use (z,x,y)   
	integer, intent(in), dimension(3), optional			:: index_transpose
   ! Array containing data distributed on the grid
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in):: asend
   

	!Number of halos
	integer	:: nh = 1    

	!Neighbours
	integer								:: nneighbors   
	integer,dimension(:),allocatable	:: id_neighbors

    ! local indices 
	integer	:: jcmin_lim,jcmax_lim
    integer :: ix,iy,iz,bicmin,bicmax,bjcmin,bjcmax,bkcmin,bkcmax
	integer	:: npercell

    ! auxiliaries 
    integer								:: nbr,ndata,itag,destid
    integer,dimension(3)				:: pcoords
	real(kind=kind(0.d0)), allocatable 	:: vbuf(:)

	! This local CFD domain is outside MD overlap zone 
	!if (CPL_OLAP_COMM .eq. MPI_COMM_NULL) return
	if (olap_mask(rank_world) .eq. 0) return

	!Revert to default j domain sending - top of overlap for CFD and bottom of overlap for MD
	if ((present(jcmax_send)) .and. & 
		(present(jcmin_send))) then
			jcmax_lim = jcmax_send
			jcmin_lim = jcmin_send
	elseif ((.not. present(jcmax_send)) .and. & 
		    (.not. present(jcmin_send))) then
		if (realm .eq. cfd_realm) then
			jcmax_lim = jcmax_olap
			jcmin_lim = jcmax_olap
		elseif (realm .eq. md_realm) then
			jcmax_lim = jcmin_olap!+1
			jcmin_lim = jcmin_olap
		endif
	else
		call error_abort("CPL_send error - both maximum and minimum j limits required and only one supplied")
	endif

    ! Re-order indices of x,y,z coordinates (handles transpose arrays)
    ix = 1 ; iy = 2; iz = 3
    if( present(index_transpose)) then
        ix = index_transpose(1)
        iy = index_transpose(2)
        iz = index_transpose(3)
    endif

    ! Number of components at each grid point
    npercell = size(asend,1)

	!Get neighbours
	call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
	allocate(id_neighbors(nneighbors))
	call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr )

	!print*, 'send data',rank_world,rank_realm,rank_olap,olap_mask(rank_world),bicmax,bjcmax,bjcmin,bkcmax,bkcmin
   
    ! loop through the maps and send the corresponding sections of asend
    do nbr = 1, nneighbors

		!Get taget processor from mapping
        destid = id_neighbors(nbr)

	    ! Get size of data to Send
		if (realm .eq. cfd_realm) then
			!Data to send is based on destination processor
			call CPL_Cart_coords(CPL_GRAPH_COMM, destid+1,  md_realm, 3, pcoords, ierr) 
	        bicmin = icPmin_md(pcoords(1))
	        bicmax = icPmax_md(pcoords(1)) 
	        bjcmin = jcmin_lim
	        bjcmax = jcmax_lim 
	        bkcmin = kcPmin_md(pcoords(3)) 
	        bkcmax = kcPmax_md(pcoords(3)) 

	        ! Amount of data to be sent
	        ndata = npercell * (bicmax-bicmin+1) * (bjcmax-bjcmin+1) * (bkcmax-bkcmin+1)
			if (allocated(vbuf)) deallocate(vbuf); allocate(vbuf(ndata))

			vbuf(1:ndata) = reshape(asend(:,bicmin:bicmax,bjcmin:bjcmax,bkcmin:bkcmax), (/ ndata /))

		elseif (realm .eq. md_realm) then
			!Data to send is based on current processor
	        bicmin = icPmin_md(iblock_realm)
	        bicmax = icPmax_md(iblock_realm) 
	        bjcmin = jcmin_lim
	        bjcmax = jcmax_lim
	        bkcmin = kcPmin_md(kblock_realm) 
	        bkcmax = kcPmax_md(kblock_realm) 

        	! Amount of data to be sent
	        ndata = npercell * (bicmax-bicmin+1) * (bjcmax-bjcmin+1) * (bkcmax-bkcmin+1)
			if (allocated(vbuf)) deallocate(vbuf); allocate(vbuf(ndata))

			vbuf(1:ndata) = reshape(asend, (/ ndata /))
		endif

		!print'(a,17i4,f20.5)', 'send data',rank_world,rank_realm,rank_olap,ndata, & 
								!size(asend),iblock_realm,jblock_realm,kblock_realm,pcoords, & 
								!bicmax,bicmin,bjcmax,bjcmin,bkcmax,bkcmin,vbuf(10)

        ! Send data 
        itag = 0 !mod( ncalls, MPI_TAG_UB) !Attention ncall could go over max tag value for long runs!!
		call MPI_send(vbuf, ndata, MPI_DOUBLE_PRECISION, destid, itag, CPL_GRAPH_COMM, ierr)
    enddo

end subroutine CPL_send_xd

!=============================================================================
! CPL_recv_data wrapper for 3d arrays
! see CPL_recv_xd for input description
!-----------------------------------------------------------------------------
subroutine CPL_recv_3d(temp,jcmax_recv,jcmin_recv,index_transpose)
    implicit none

 	integer, optional, intent(in)				      :: jcmax_recv,jcmin_recv

    integer, optional, intent(in) :: index_transpose(3)  
    real(kind(0.d0)),dimension(:,:,:),intent(inout) :: temp 
                                                        
	integer n1, n2, n3, n4
	real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: arecv    

    n1 = 1 
    n2 = size(temp,1)
    n3 = size(temp,2)
    n4 = size(temp,3)

	!Add padding column to 3D array to make it 4D
	allocate(arecv(n1,n2,n3,n4))
    call CPL_recv_xd(arecv,jcmax_recv,jcmin_recv,index_transpose)
	temp(:,:,:) = 	arecv(1,:,:,:) 

end subroutine CPL_recv_3d

!=============================================================================
! CPL_recv_data wrapper for 4d arrays
! see CPL_recv_xd for input description
!-----------------------------------------------------------------------------
subroutine CPL_recv_4d(arecv,jcmax_recv,jcmin_recv,index_transpose)

    implicit none

 	integer, optional, intent(in)				      :: jcmax_recv,jcmin_recv
    integer, optional, intent(in) :: index_transpose(3)  
    real(kind(0.d0)),dimension(:,:,:,:),intent(inout) :: arecv

    call CPL_recv_xd(arecv,jcmax_recv,jcmin_recv,index_transpose)

end subroutine CPL_recv_4d

!=============================================================================
! Receive data from to local grid from the associated ranks from the other 
! realm
!-----------------------------------------------------------------------------
subroutine CPL_recv_xd(arecv,jcmax_recv,jcmin_recv,index_transpose)
    use mpi
	use coupler_module
    implicit none

    ! Minimum and maximum values of j to receive
 	integer, optional, intent(in)				      :: jcmax_recv,jcmin_recv
    ! specify the order of dimensions in asend default is ( x, y, z) but some CFD solvers use (z,x,y)   
    integer, optional, dimension(3), intent(in) 	  :: index_transpose(3) 
    ! Array that recieves grid distributed data 
    real(kind(0.d0)), dimension(:,:,:,:),intent(inout):: arecv     

	!Neighbours
	integer								:: nneighbors   
	integer,dimension(:),allocatable	:: id_neighbors
                                                         
    ! local indices 
	integer	:: jcmax_lim, jcmin_lim
    integer :: n,nbr,i,j,k,ix,iy,iz,bicmin,bicmax,bjcmin,bjcmax,bkcmin,bkcmax
	integer	:: ncl(2,3),pcoords(3),recvsize,startbuf,endbuf,npercell

    ! auxiliaries 
    integer	:: ndata, itag, sourceid,source_realm,start_address
	integer,dimension(:),allocatable   :: req
	integer,dimension(:,:),allocatable :: status
	integer,dimension(:,:,:),allocatable :: ncl_recv
    real(kind(0.d0)),dimension(:), allocatable ::  vbuf,vbuf2
 
	! This local CFD domain is outside MD overlap zone 
	if (olap_mask(rank_world).eq.0) return
	!if (CPL_OLAP_COMM .eq. MPI_COMM_NULL) return

	!Revert to default j domain receive - top of overlap for CFD and bottom of overlap for MD
	if ((present(jcmax_recv)) .and. & 
		(present(jcmin_recv))) then
			jcmax_lim = jcmax_recv
			jcmin_lim = jcmin_recv
	elseif ((.not. present(jcmax_recv)) .and. & 
		    (.not. present(jcmin_recv))) then
		if (realm .eq. cfd_realm) then
			jcmax_lim = jcmax_olap
			jcmin_lim = jcmax_olap
		elseif (realm .eq. md_realm) then
			jcmax_lim = jcmin_olap
			jcmin_lim = jcmin_olap
		endif
	else
		call error_abort("CPL_recv error - both maximum and minimum j limits required and only one supplied")
	endif

    ! Local grid box
    ! Get local grid box ranges seen by this rank for either CFD or MD
    if (realm .eq. cfd_realm) then 
		!Load CFD cells per processor
        bicmin = icPmin_cfd(iblock_realm)
        bicmax = icPmax_cfd(iblock_realm) 
		bjcmin = jcmin_lim 
        bjcmax = jcmax_lim 
        bkcmin = kcPmin_cfd(kblock_realm) 
        bkcmax = kcPmax_cfd(kblock_realm)

		! Amount of data to receive from all MD processors
		npercell = size(arecv,1)
		ndata = npercell * (bicmax-bicmin+1) * (bjcmax-bjcmin+1) * (bkcmax-bkcmin+1)
		allocate(vbuf(ndata)); vbuf = 0.d0

	endif

    ! Get the indices in x,y,z direction from transpose array
    ix = 1 ; iy = 2; iz = 3
    if( present(index_transpose)) then
        ix = index_transpose(1)
        iy = index_transpose(2)
        iz = index_transpose(3)
    endif
	!ncl(1,ix) = bicmin;	ncl(2,ix) = bicmax
	!ncl(1,iy) = bjcmin;	ncl(2,iy) = bjcmax
	!ncl(1,iz) = bkcmin;	ncl(2,iz) = bkcmax

	!Get neighbours
	call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
	allocate(id_neighbors(nneighbors))
	call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr )

    ! Receive from all attached processors
	allocate(req(nneighbors))
	allocate(status(MPI_STATUS_SIZE,nneighbors))
	allocate(ncl_recv(2,3,nneighbors))
    start_address = 1 
    do nbr = 1, nneighbors

		!Get source processor from topology graph
        sourceid =  id_neighbors(nbr)

	    ! Get size of data to receive from source processors
		if (realm .eq. cfd_realm) then

			!CFD realm receives data based on size of MD processor domain
			call CPL_Cart_coords(CPL_GRAPH_COMM, sourceid+1, md_realm, 3, pcoords, ierr) 

	        bicmin = icPmin_md(pcoords(1))
	        bicmax = icPmax_md(pcoords(1)) 
	        bjcmin = jcmin_lim 
	        bjcmax = jcmax_lim 
	        bkcmin = kcPmin_md(pcoords(3)) 
	        bkcmax = kcPmax_md(pcoords(3))

			! Amount of data to receive
			npercell = size(arecv,1)
			ndata = npercell * (bicmax-bicmin+1) * (bjcmax-bjcmin+1) * (bkcmax-bkcmin+1)

			!increment pointer ready to receive next piece of data		
			start_address = 1+((pcoords(1)-1) + (pcoords(3)-1)*npx_md)*ndata

		elseif (realm .eq. md_realm) then

			!MD realm receives data as big as own processor domain
	        bicmin = icPmin_md(iblock_realm)
	        bicmax = icPmax_md(iblock_realm) 
	        bjcmin = jcmin_lim
	        bjcmax = jcmax_lim 
	        bkcmin = kcPmin_md(kblock_realm) 
	        bkcmax = kcPmax_md(kblock_realm)

			! Amount of data to receive
			npercell = size(arecv,1)
			ndata = npercell * (bicmax-bicmin+1) * (bjcmax-bjcmin+1) * (bkcmax-bkcmin+1) 
			start_address = 1

			if (size(arecv) .ne. ndata) then
				!print*, 'size of expected recv data = ',  size(arecv), 'actual size of recv data = ', ndata
				call error_abort("domain size mismatch in recv data")
			endif
			! Amount of data to receive
			allocate(vbuf(ndata)); vbuf = 0.d0

		endif
		! Receive section of data
		itag = 0
		!print'(a,8i8)', 'recv data',realm,sourceid+1,ndata,size(arecv),start_address,pcoords
		call MPI_irecv(vbuf(start_address),ndata,MPI_DOUBLE_PRECISION,sourceid,itag,&
            						CPL_GRAPH_COMM,req(nbr),ierr)

    enddo
    call MPI_waitall(nneighbors, req, status, ierr)


	if (realm .eq. cfd_realm) then
		do nbr=1, nneighbors

			!CFD realm receives data based on size of MD processor domain
			call CPL_Cart_coords(CPL_GRAPH_COMM, sourceid+1, md_realm, 3, pcoords, ierr) 

	        bicmin = icPmin_md(pcoords(1))
	        bicmax = icPmax_md(pcoords(1)) 
	        bjcmin = jcmin_lim 
	        bjcmax = jcmax_lim 
	        bkcmin = kcPmin_md(pcoords(3)) 
	        bkcmax = kcPmax_md(pcoords(3))

			! Amount of data to receive
			npercell = size(arecv,1)
			ndata = npercell * (bicmax-bicmin+1) * (bjcmax-bjcmin+1) * (bkcmax-bkcmin+1)
			!increment pointer ready to receive next piece of data		
			start_address = 1+((pcoords(1)-1) + (pcoords(3)-1)*npx_md)*ndata

			arecv(:,bicmin:bicmax,bjcmin:bjcmax,bkcmin:bkcmax) =  & 
						reshape(vbuf(start_address:start_address+ndata-1), & 
								(/npercell,bicmax-bicmin+1,bjcmax-bjcmin+1,bkcmax-bkcmin+1 /))	
		enddo
	endif


	!do nbr=1, nneighbors
	!	arecv(:,ncl_recv(1,1,nbr):ncl_recv(2,1,nbr), & 
	!			1:1, & 
	!			ncl_recv(1,3,nbr):ncl_recv(2,3,nbr)) = reshape(vbuf,(/ 3,ncl_recv(1,1,nbr)-ncl_recv(2,1,nbr), & 
	!															 	 1, & 
	!															 	 ncl_recv(1,3,nbr)-ncl_recv(2,3,nbr)/))
	!enddo

	!arecv = reshape(vbuf,(/ 3, ncl(2,1)-ncl(1,1), ncl(2,2)-ncl(1,2), ncl(2,3)-ncl(1,3) /))

	do n = 1,size(vbuf)
		print*, 'vbuf', n, vbuf(n)
	enddo

	!do n = 1,size(arecv,1)
	!do i = 1,size(arecv,2)
	!do j = 1,size(arecv,3)
	!do k = 1,size(arecv,4)
	!	print*, 'arecv', n,i,j,k, arecv(n,i,j,k)
	!enddo
	!enddo
	!enddo
	!enddo
           
end subroutine CPL_recv_xd

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
	use coupler_module, only : rank_world, nproc_world, CPL_WORLD_COMM, VOID
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

!-------------------------------------------------------------------

subroutine CPL_pack(unpacked,packed,realm)
	use coupler_module, only : CPL_CART_COMM,rank_cart,md_realm,cfd_realm, & 
	                           error_abort,CPL_GRAPH_COMM,myid_graph,realm_name
	implicit none

	integer, intent(in)											:: realm
	real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in)		:: unpacked
	real(kind=kind(0.d0)),dimension(:),allocatable, intent(out)	:: packed

	integer                          :: pos,n,nbr,id_nbr,icell,jcell,kcell,ierr
	integer                          :: npercell,ncells,nneighbors
	integer,dimension(3)			 :: coord
	integer,dimension(6)			 :: extents,gextents
	integer,dimension(:),allocatable :: id_neighbors

	!Amount of data per cell
	npercell = size(unpacked,1)

	!Allocate packing buffer
	if (allocated(packed)) deallocate(packed)
	allocate(packed(size(unpacked)))

	! Get neighbour topology to determine ordering of packed data
	call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
    allocate(id_neighbors(nneighbors))
    call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr)

	! Loop through all neighbours which will be sent to and order the data 
	! appropriately to send each correctly
	do nbr = 1,nneighbors

		if (realm .eq. cfd_realm) then
			! Get MD neighbour
			id_nbr = id_neighbors(nbr)
			call CPL_Cart_coords(CPL_GRAPH_COMM,id_nbr+1,md_realm,3,coord,ierr) 
			call CPL_proc_extents(coord,md_realm,extents,ncells)
			! Get offset of neighbouring processor
			pos = id_nbr * npercell * ncells


			!print*,'Pack',rank_cart,realm_name(realm),coord,nbr,id_nbr,extents,coord

		elseif (realm .eq. md_realm) then
			!Get own processor coordinates 
			call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
			call CPL_proc_extents(coord,realm,gextents,ncells)
			! Get local extents
			extents(1) = 1;	extents(2) = gextents(2)-gextents(1)
			extents(3) = 1;	extents(4) = gextents(4)-gextents(3)
			extents(5) = 1;	extents(6) = gextents(6)-gextents(5)
			pos = 1
		endif



		! Pack array into buffer
		do kcell=extents(5),extents(6)
		do jcell=extents(3),extents(4)
		do icell=extents(1),extents(2)
		do n = 1,npercell

			packed(pos) = unpacked(n,icell,jcell,kcell)
			pos = pos + 1

		end do
		end do
		end do
		end do

	end do

	!Sanity check
	if (size(packed) .ne. npercell*ncells) then
		!print*, 'data size', size(packed), 'expected size', npercell*ncells
		call error_abort("CPL_pack error - cell array does not match expected extents")
	endif

end subroutine CPL_pack


!-------------------------------------------------------------------

subroutine CPL_unpack(packed,unpacked,realm)
	use coupler_module, only : CPL_CART_COMM,rank_cart,md_realm,cfd_realm, & 
	                           error_abort,CPL_GRAPH_COMM,myid_graph
	implicit none

	integer, intent(in)											      :: realm
	real(kind=kind(0.d0)),dimension(:,:,:,:),allocatable, intent(out) :: unpacked
	real(kind=kind(0.d0)),dimension(:),allocatable, intent(inout)     :: packed

	integer                          :: pos,n,nbr,id_nbr,icell,jcell,kcell,ierr
	integer                          :: npercell,ncells,nneighbors
	integer,dimension(3)			 :: coord
	integer,dimension(6)			 :: extents,gextents
	integer,dimension(:),allocatable :: id_neighbors

	call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
	call CPL_proc_extents(coord,realm,extents,ncells)

	!Amount of data per cell
	npercell = size(packed)/ncells

	!Allocate packing buffer
	if (allocated(unpacked)) deallocate(unpacked)
	allocate(unpacked(npercell,1:extents(2)-extents(1), &
	                 		   1:extents(4)-extents(3), &
	                 		   1:extents(6)-extents(5)))

	! Get neighbour topology to determine ordering of packed data
	call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
    allocate(id_neighbors(nneighbors))
    call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr)

	! Loop through all neighbours which will be sent to and order the data 
	! appropriately to send each correctly
	do nbr = 1,nneighbors

		if (realm .eq. cfd_realm) then
			! Get MD neighbour
			id_nbr = id_neighbors(nbr)
			call CPL_Cart_coords(CPL_GRAPH_COMM,id_nbr,md_realm,3,coord,ierr) 
			call CPL_proc_extents(coord,md_realm,extents,ncells)
			! Get offset of neighbouring processor
			pos = id_nbr * npercell * ncells	!ASSUMES all ncell the same!!
		elseif (realm .eq. md_realm) then
			!Get own processor coordinates 
			call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
			call CPL_proc_extents(coord,realm,gextents,ncells)
			! Get local extents
			extents(1) = 1;	extents(2) = gextents(2)-gextents(1)
			extents(3) = 1;	extents(4) = gextents(4)-gextents(3)
			extents(5) = 1;	extents(6) = gextents(6)-gextents(5)
			pos = 1
		endif

		! Unpack buffer into array
		do kcell=extents(5),extents(6)
		do jcell=extents(3),extents(4)
		do icell=extents(1),extents(2)
		do n = 1,npercell

			unpacked(n,icell,jcell,kcell) = packed(pos)
			pos = pos + 1

		end do
		end do
		end do
		end do

	end do

	!Deallocate packed buffer
	deallocate(packed)

end subroutine CPL_unpack


!-------------------------------------------------------------------
! 					CPL_proc_extents  			      -
!-------------------------------------------------------------------

! Get maximum and minimum cells for current communicator

! - - - Synopsis - - -

! CPL_proc_extents(coord,realm,extents,ncells)

! - - - Input Parameters - - -

!coord
!    processor cartesian coordinate (3 x integer) 
!realm
!    cfd_realm (1) or md_realm (2) (integer) 
!
! - - - Output Parameter - - -

!extents
!	 Six components array which defines processor extents
!	 xmin,xmax,ymin,ymax,zmin,zmax (6 x integer) 
!ncells (optional)
!    number of cells on processor (integer) 


subroutine CPL_proc_extents(coord,realm,extents,ncells)
	use mpi
	use coupler_module, only: md_realm,      cfd_realm,      &
	                          icPmin_md,     icPmax_md,      &
	                          jcPmin_md,     jcPmax_md,      &
	                          kcPmin_md,     kcPmax_md,      &
	                          icPmin_cfd,    icPmax_cfd,     &
	                          jcPmin_cfd,    jcPmax_cfd,     &
	                          kcPmin_cfd,    kcPmax_cfd,     &
	                          error_abort
	implicit none

	integer, intent(in)  :: coord(3), realm
	integer, intent(out) :: extents(6)
	integer, optional, intent(out) :: ncells

	select case(realm)
	case(md_realm)
		extents = (/icPmin_md(coord(1)),icPmax_md(coord(1)), & 
		            jcPmin_md(coord(2)),jcPmax_md(coord(2)), & 
		            kcPmin_md(coord(3)),kcPmax_md(coord(3))/)
	case(cfd_realm)
		extents = (/icPmin_cfd(coord(1)),icPmax_cfd(coord(1)), & 
		            jcPmin_cfd(coord(2)),jcPmax_cfd(coord(2)), & 
		            kcPmin_cfd(coord(3)),kcPmax_cfd(coord(3))/)

	case default
		call error_abort('Wrong realm in rank_cart_to_cell_extents')
	end select

	if (present(ncells)) then
		ncells = (extents(2) - extents(1) + 1) * &
				 (extents(4) - extents(3) + 1) * &
				 (extents(6) - extents(5) + 1)
	end if

end subroutine CPL_proc_extents


!-------------------------------------------------------------------
! 					CPL_Cart_coords								   -
!-------------------------------------------------------------------

! Determines process coords in appropriate realm's cartesian topology 
! given a rank in any communicator

! - - - Synopsis - - -

! CPL_Cart_coords(COMM, rank, realm, maxdims, coords, ierr)

! - - - Input Parameters - - -

!comm
!    communicator with cartesian structure (handle) 
!realm
!    cfd_realm (1) or md_realm (2) (integer) 
!rank
!    rank of a process within group of comm (integer) 
!    NOTE fortran convention rank=1 to nproc
!maxdims
!    length of vector coords in the calling program (integer) 

! - - - Output Parameter - - -

!coords
!    integer array (of size ndims) containing the Cartesian coordinates 
!    of specified process (integer) 
!ierr
!    error flag

subroutine CPL_Cart_coords(COMM, rank, realm, maxdims, coords, ierr)
	use coupler_module, only :  CPL_WORLD_COMM, CPL_REALM_COMM, CPL_INTER_COMM, & 
								CPL_CART_COMM, CPL_OLAP_COMM, CPL_GRAPH_COMM,   &
								CPL_REALM_INTERSECTION_COMM, md_realm,cfd_realm, &
								rank_world2rank_mdrealm,rank_world2rank_mdcart,     &
								rank_world2rank_cfdrealm,rank_world2rank_cfdcart,    &
								rank_world2rank_olap,rank_world2rank_graph,      &
								rank_world2rank_inter,rank_mdrealm2rank_world,    &
								rank_mdcart2rank_world,rank_cfdrealm2rank_world,    &
								rank_cfdcart2rank_world,rank_olap2rank_world,    &
								rank_graph2rank_world,rank_inter2rank_world,	&
								rank2coord_cfd,	rank2coord_md, &
								COUPLER_ERROR_CART_COMM, VOID, nproc_cfd,nproc_md, & 
								error_abort
	implicit none

	integer, intent(in)		:: COMM, realm, rank, maxdims
	integer, intent(out)	:: coords(maxdims), ierr

	integer					:: worldrank, cartrank

	!Get rank in world COMM from current COMM
	if (COMM .eq. CPL_WORLD_COMM) then
		! -  -  -  -  -  -  -  -  -  -  -  -  -
		worldrank = rank
		! -  -  -  -  -  -  -  -  -  -  -  -  -
	elseif(COMM .eq. CPL_REALM_COMM) then
		! -  -  -  -  -  -  -  -  -  -  -  -  -
		if (realm .eq. cfd_realm) then
			if (allocated(rank_cfdrealm2rank_world) .eqv. .false.) then
				call error_abort("CPL_Cart_coords Error - " // &
								 "Setup not complete for CFD CPL_REALM_COMM")
			elseif (rank .gt. size(rank_cfdrealm2rank_world)) then
				print*, 'rank = ', rank, 'comm size = ', size(rank_cfdrealm2rank_world)
				call error_abort("CPL_Cart_coords Error -Specified rank is not in CFD realm")
			endif
			worldrank = rank_cfdrealm2rank_world(rank)
		elseif (realm .eq. md_realm) then
			if (allocated(rank_mdrealm2rank_world) .eqv. .false.) then
				call error_abort("CPL_Cart_coords Error - Setup not " // &
								 "complete for MD CPL_REALM_COMM")
			elseif (rank .gt. size(rank_mdrealm2rank_world)) then
				print*, 'rank = ', rank, 'comm size = ', size(rank_cfdrealm2rank_world)
				call error_abort("CPL_Cart_coords Error -Specified rank " // &
								 "is not in MD realm")
			endif
			worldrank = rank_mdrealm2rank_world(rank)
		endif
		! -  -  -  -  -  -  -  -  -  -  -  -  -
	elseif(COMM .eq. CPL_CART_COMM) then
		! -  -  -  -  -  -  -  -  -  -  -  -  -
		if (realm .eq. cfd_realm) then
			if (allocated(rank2coord_cfd) .eqv. .false.) &
				call error_abort("CPL_Cart_coords Error - Setup not complete" // &
								 " for CFD CPL_CART_COMM")
			coords = rank2coord_cfd(:,rank)
		elseif (realm .eq. md_realm) then
			if (allocated(rank2coord_md) .eqv. .false.) &
				call error_abort("CPL_Cart_coords Error - Setup not complete " // &
								 "for MD CPL_CART_COMM")
			coords = rank2coord_md(:,rank)
		endif
		return
		! -  -  -  -  -  -  -  -  -  -  -  -  -
	elseif(COMM .eq. CPL_OLAP_COMM) then
		! -  -  -  -  -  -  -  -  -  -  -  -  -
		if (allocated(rank_olap2rank_world) .eqv. .false.) then
			call error_abort("CPL_Cart_coords Error - Setup not complete for CPL_OLAP_COMM")
		elseif (rank .gt. size(rank_olap2rank_world)) then
			print*, 'rank = ', rank, 'comm size = ', size(rank_olap2rank_world)
			call error_abort("CPL_Cart_coords Error - Specified rank is not in overlap")
		endif
		worldrank = rank_olap2rank_world(rank)
		! -  -  -  -  -  -  -  -  -  -  -  -  -
	elseif(COMM .eq. CPL_GRAPH_COMM) then
		! -  -  -  -  -  -  -  -  -  -  -  -  -
		if (allocated(rank_graph2rank_world) .eqv. .false.) then
			call error_abort("CPL_Cart_coords Error - Setup not complete for CPL_GRAPH_COMM")
		elseif (rank .gt. size(rank_graph2rank_world)) then
			call error_abort("CPL_Cart_coords Error - Specified rank is not in graph")
		endif
		worldrank = rank_graph2rank_world(rank)
		! -  -  -  -  -  -  -  -  -  -  -  -  -
	elseif(COMM .eq. CPL_REALM_INTERSECTION_COMM) then
		! -  -  -  -  -  -  -  -  -  -  -  -  -
		call error_abort("CPL_Cart_coords Error - Intersection COMM not programmed")
		! -  -  -  -  -  -  -  -  -  -  -  -  -

	elseif(COMM .eq. CPL_INTER_COMM) then
		! -  -  -  -  -  -  -  -  -  -  -  -  -
		call error_abort("CPL_Cart_coords Error - Intercomm between realms id" // &
								 " - use realm comms instead")
		ierr = COUPLER_ERROR_CART_COMM 
		return
		! -  -  -  -  -  -  -  -  -  -  -  -  -
	else 
		! -  -  -  -  -  -  -  -  -  -  -  -  -
		call error_abort("CPL_Cart_coords Error - Unknown COMM")
		ierr = COUPLER_ERROR_CART_COMM 
		return
		! -  -  -  -  -  -  -  -  -  -  -  -  -
	endif
	
	!Get rank in realm cartesian communicator
	if (realm .eq. cfd_realm) then
		if (allocated(rank_world2rank_cfdcart) .eqv. .false.) then
			call error_abort("CPL_Cart_coords Error - world to cart mapping " // &
								 "not initialised correctly")
		endif
		cartrank = rank_world2rank_cfdcart(worldrank)
		if (cartrank .eq. VOID) call error_abort("CPL_Cart_coords Error - void element in mapping")
		if (cartrank .gt. nproc_cfd) call error_abort("CPL_Cart_coords Error - rank not in cfd realm")
	elseif (realm .eq. md_realm) then
		if (allocated(rank_world2rank_mdcart) .eqv. .false.) then
			call error_abort("CPL_Cart_coords Error - world to cart " // &
								 "mapping not initialised correctly")
		endif
		cartrank = rank_world2rank_mdcart(worldrank)
		if (cartrank .eq. VOID) call error_abort("CPL_Cart_coords Error - void element in mapping")
		if (cartrank .gt. nproc_md) call error_abort("CPL_Cart_coords Error - rank not in md realm")
	endif

	!Get cartesian coordinate in appropriate realm
	if (realm .eq. cfd_realm) then
		coords = rank2coord_cfd(:,cartrank)
	elseif (realm .eq. md_realm) then
		coords = rank2coord_md(:,cartrank)
	endif

	!Success
	ierr = 0
 
end subroutine CPL_Cart_coords


end module coupler
