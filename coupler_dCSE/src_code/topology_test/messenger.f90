subroutine create_realms
	use coupler_module
	use mpi
	implicit none

	integer :: &
		gridsize(3), &
		coord(3)
	logical, dimension(3), parameter :: &
		periodicity = (/.true.,.false.,.true./)

	if (rank_world.le.nproc_md) then
		realm = md_realm
		gridsize = (/npx_md,npy_md,npz_md/)	
	else
		realm = cfd_realm
		gridsize = (/npx_cfd,npy_cfd,npz_cfd/)	
	end if

	call MPI_comm_split(CPL_WORLD_COMM,realm,myid_world,CPL_REALM_COMM,ierr)
	call MPI_cart_create(CPL_REALM_COMM,3,gridsize,periodicity,.true., &
						 CPL_CART_COMM,ierr)

	call MPI_comm_rank(CPL_REALM_COMM,myid_realm,ierr)
	rank_realm = myid_realm + 1
	call MPI_comm_rank(CPL_CART_COMM,myid_cart,ierr)
	rank_cart = myid_cart + 1

	call MPI_cart_coords(CPL_CART_COMM,myid_cart,3,coord,ierr)
	coord(:) = coord(:) + 1

	if (realm .eq. md_realm) then
		coord2rank_md(coord(1),coord(2),coord(3)) = rank_realm
		rank2coord_md(:,rank_realm)    = coord(:)
		rank_md2rank_world(rank_realm) = rank_world
	else if (realm .eq. cfd_realm) then
		coord2rank_cfd(coord(1),coord(2),coord(3)) = rank_realm
		rank2coord_cfd(:,rank_realm)    = coord(:)
		rank_cfd2rank_world(rank_realm) = rank_world
	end if

	call collect_coord2ranks
	call collect_rank2coords
	call collect_rank2ranks
	call write_realm_info

end subroutine create_realms

subroutine collect_coord2ranks
	use mpi
	use coupler_module
	implicit none

	integer, dimension(:), allocatable :: mbuf, cbuf

	allocate(mbuf(nproc_md))
	allocate(cbuf(nproc_cfd))

	call MPI_allreduce(coord2rank_md,mbuf,nproc_md,MPI_INTEGER,MPI_SUM,   &
	                   CPL_WORLD_COMM,ierr)	
	call MPI_allreduce(coord2rank_cfd,cbuf,nproc_cfd,MPI_INTEGER,MPI_SUM, &
	                   CPL_WORLD_COMM,ierr)	

	coord2rank_md  = reshape(mbuf,(/npx_md,npy_md,npz_md/))                  
	coord2rank_cfd = reshape(cbuf,(/npx_cfd,npy_cfd,npz_cfd/))

	deallocate(mbuf)
	deallocate(cbuf)

end subroutine collect_coord2ranks

subroutine collect_rank2coords
	use mpi
	use coupler_module
	implicit none

	!integer :: n
	integer, dimension(:), allocatable :: mbuf, cbuf

	allocate(mbuf(3*nproc_md))
	allocate(cbuf(3*nproc_cfd))

	call MPI_allreduce(rank2coord_md,mbuf,3*nproc_md,MPI_INTEGER,MPI_SUM,   &
	                   CPL_WORLD_COMM,ierr)	
	call MPI_allreduce(rank2coord_cfd,cbuf,3*nproc_cfd,MPI_INTEGER,MPI_SUM, &
	                   CPL_WORLD_COMM,ierr)	

	rank2coord_md  = reshape(mbuf,(/3,nproc_md/))                  
	rank2coord_cfd = reshape(cbuf,(/3,nproc_cfd/))
	
!	if (myid_world.eq.0) then
!		do n = 1,nproc_cfd
!			print('(4i5)'), n, rank2coord_cfd(1,n),rank2coord_cfd(2,n),rank2coord_cfd(3,n)
!		end do
!	end if

	deallocate(mbuf)
	deallocate(cbuf)

end subroutine collect_rank2coords

subroutine collect_rank2ranks
	use mpi
	use coupler_module
	implicit none
	
	integer, dimension(:), allocatable :: mbuf, cbuf

	allocate(mbuf(nproc_md))
	allocate(cbuf(nproc_cfd))

	call MPI_allreduce(rank_md2rank_world,mbuf,nproc_md,MPI_INTEGER,  &
	                   MPI_SUM, CPL_WORLD_COMM,ierr)	
	call MPI_allreduce(rank_cfd2rank_world,cbuf,nproc_cfd,MPI_INTEGER, &
	                   MPI_SUM, CPL_WORLD_COMM,ierr)	

	rank_md2rank_world  = mbuf
	rank_cfd2rank_world = cbuf

	deallocate(mbuf)
	deallocate(cbuf)

end subroutine collect_rank2ranks

subroutine prepare_overlap_comms
	use coupler_module
	use mpi
	implicit none

	!loop over cfd cart ranks
	! find cfd cart coords from cfd cart rank
	!  find md cart coords (from cfd_icoord2olap_md_jcoords)
	!   find md cart rank from md cart coords (coord2rank_md) 
	!    find md world rank from md cart rank (rank_md2rank_world)
	!     set group(md_world_rank) to cfd cart rank
	!      split comm to groups
	!       if group(world_rank) == 0, set olap_comm to null 

	integer :: i,j,k,ic,jc,kc
	integer :: trank_md, trank_cfd, trank_world
	integer, dimension(:), allocatable :: mdicoords, mdjcoords, mdkcoords
	integer, parameter :: olap_null = -666
	integer :: group(nproc_world)
	integer :: cfdcoord(3)

	allocate(mdicoords(npx_md/npx_cfd))
	allocate(mdjcoords(npy_md/npy_cfd))
	allocate(mdkcoords(npz_md/npz_cfd))

	group = olap_null

	do trank_cfd = 1,nproc_cfd

		cfdcoord(:) = rank2coord_cfd(:,trank_cfd)
		mdicoords(:) = cfd_icoord2olap_md_icoords(cfdcoord(1),:)
		mdjcoords(:) = cfd_jcoord2olap_md_jcoords(cfdcoord(2),:)
		mdkcoords(:) = cfd_kcoord2olap_md_kcoords(cfdcoord(3),:)

		do i = 1,size(mdicoords)
		do j = 1,size(mdjcoords)
		do k = 1,size(mdkcoords)
			ic = mdicoords(i)
			jc = mdjcoords(j)
			kc = mdkcoords(k)
			if (any((/ic,jc,kc/).eq.olap_null)) cycle
			trank_md = coord2rank_md(ic,jc,kc)
			trank_world = rank_md2rank_world(trank_md)
			group(trank_world) = trank_cfd
		end do
		end do	
		end do
			
		trank_world = rank_cfd2rank_world(trank_cfd)
		group(trank_world) = trank_cfd

	end do

	! Split world Comm into a set of comms for overlapping processors
	call MPI_comm_split(CPL_WORLD_COMM,group(rank_world),realm, &
	                    CPL_OLAP_COMM,ierr)

	call MPI_comm_rank(CPL_OLAP_COMM,myid_olap,ierr)
	rank_olap = myid_olap + 1	
	if (myid_olap .eq. 0) testval = group(rank_world)
	call MPI_bcast(testval,1,MPI_INTEGER,0,CPL_OLAP_COMM,ierr)

	! Set all non-overlapping processors to MPI_COMM_NULL
	call MPI_comm_size(CPL_OLAP_COMM, nproc_olap, ierr)
	if (group(rank_world).eq.olap_null) then
		CPL_OLAP_COMM = MPI_COMM_NULL
	elseif (nproc_olap .eq. 1) then
		CPL_OLAP_COMM = MPI_COMM_NULL
	endif

	deallocate(mdicoords)
	deallocate(mdjcoords)
	deallocate(mdkcoords)

	if (realm.eq.md_realm) call write_overlap_comms_md

	!TODO OLAP COMM NULL

end subroutine prepare_overlap_comms



!Setup topology graph of overlaps between CFD & MD processors

subroutine CPL_overlap_topology
	use coupler_module
	use mpi
	implicit none

	integer								:: i, n, nneighbors, myid_graph, nconnections
	integer, dimension(:),allocatable	:: index, edges, neighbors
	logical								:: reorder

	!Allow optimisations of ordering
	reorder = .true.

	!Get number of processors in communicating overlap region 
	!call MPI_comm_rank(CPL_OLAP_COMM, myid_graph, ierr)
	if (CPL_OLAP_COMM .ne. MPI_COMM_NULL) then

		!call MPI_comm_size(CPL_OLAP_COMM, nproc_olap, ierr)

		!CFD processor is root and has mapping to all MD processors
		allocate(index(nproc_olap))			!Index for each processor
		allocate(edges(2*(nproc_olap)-1))	!nproc_olap-1 for CFD and one for each of nproc_olap MD processors
		index = 0; 	edges = 0


		!select case(realm)
		!case(cfd_realm)
		!	print'(2a,i5,a,i5,a,4i5)', code_name(realm),' World procs ',myid_world+1, & 
		!								' Is olap procs ', myid_olap+1, ' of ',nproc_olap, rank2coord_cfd(:,rank_realm)
		!case(md_realm)
		!	print'(2a,i5,a,i5,a,4i5)', code_name(realm),' World procs ',myid_world+1, & 
		!								' Is olap procs ', myid_olap+1, ' of ',nproc_olap,rank2coord_md(:,rank_realm)
		!end select

		!CFD processor has connections to nproc_olap MD processors
		nconnections = nproc_olap-1
		index(1) = nconnections
		do n = 1,nconnections
			edges(n) = n !CFD connected to all MD processors 1 to nconnections
			!if (myid_olap .eq. 0) &
			!print*, 'CFD graph info',code_name(realm),myid_realm,myid_olap, index(1), edges(n)
		enddo

		!MD processor has a single connection to CFD
		nconnections = 1; i = 2
		do n = nproc_olap+1,2*(nproc_olap)-1
			index(i) = index(i-1) + nconnections !Each successive index incremented by one
			edges(n) = CFDid_olap	!Connected to CFD processor
			!if (myid_olap .eq. 0) &
			!print*, 'MD graph info',code_name(realm),myid_realm,myid_olap,n,i,index(i), edges(n)
			i = i + 1
		enddo

		!if (myid_olap .eq. 0) &
		!print*, 'final values  ', ' index = ',index, 'Edges = ',edges

		!Create graph topology for overlap region
		call MPI_Graph_create(CPL_OLAP_COMM, nproc_olap,index,edges,reorder,CPL_GRAPH_COMM,ierr)

		! - - TEST - - TEST - -
		!Get number of neighbours
		call MPI_comm_rank( CPL_GRAPH_COMM, myid_graph, ierr)
		call MPI_Graph_neighbors_count( CPL_GRAPH_COMM, myid_graph, nneighbors, ierr)
		allocate(neighbors(nneighbors))
		!Get neighbours
		call MPI_Graph_neighbors( CPL_GRAPH_COMM, myid_graph, nneighbors, neighbors,ierr )

		print*, 'My graph',myid_world,myid_graph,myid_olap, nneighbors, neighbors

	endif



end subroutine CPL_overlap_topology

subroutine gatherscatter
	use mpi
	use coupler_module
	implicit none
	
	integer :: icell, jcell, kcell, ixyz
	integer :: coord(3)
	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: u
	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: s

	integer :: sendcount
	
	if (realm .eq. md_realm) then
		
		call MPI_cart_coords(CPL_CART_COMM,myid_cart,3,coord,ierr)
		coord(:) = coord(:) + 1
	
		allocate(u(icPmin_md(coord(1)):icPmax_md(coord(1)), &
		           jcPmin_md(coord(2)):jcPmax_md(coord(2)), &
		           kcPmin_md(coord(3)):kcPmax_md(coord(3)), &
		           3))

		! Populate dummy u
		do icell=icPmin_md(coord(1)),icPmax_md(coord(1))
		do jcell=jcPmin_md(coord(2)),jcPmax_md(coord(2))
		do kcell=kcPmin_md(coord(3)),kcPmax_md(coord(3))
			do ixyz = 1,3
				u(icell,jcell,kcell,ixyz) = 1000*icell + &
				                            100 *jcell + &
				                            10  *kcell + &
				                            1   *ixyz
			end do
		end do
		end do
		end do

		!Get sendcount
		sendcount = size(u)
		!Get length of avgu on CFD side
		!Get recvcounts
		!Get displacements
		!Get CFD root ID

	else
						
		!Get sendcount
		sendcount = 0
		!Get length of avgu on CFD side
		!Get recvcounts
		!Get displacements
		!Get CFD root ID

	end if

!	print*, 'sendcount = ', sendcount	
!	call MPI_gatherv(u,-SENDCNT-,MPI_DOUBLE_PRECISION,-AVGU-,-RECVCNTS-, &
!	                 -DISPLS-,MPI_DOUBLE_PRECISION,-CFDROOTID-,CPL_OLAP_COMM)

	!CFD gets u(:) from md rank
	!convert md rank to coords
	!grab icPmin/max, jcPmin/max, kcPmin/max from coords
	!slot into correct part of array

	if (realm.eq.md_realm) deallocate(u)

end subroutine gatherscatter

subroutine write_realm_info
	use coupler_module
	use mpi
	implicit none

	integer :: coord(3)

	if (myid_world.eq.0) then
		write(1000+rank_world,*), '---------- REALM INFORMATION --------------'
		write(1000+rank_world,*), ' wrank  realm  realmrank       cart coords '
		write(1000+rank_world,*), '-------------------------------------------'
	end if

	call MPI_cart_coords(CPL_CART_COMM,myid_cart,3,coord,ierr)
	coord(:) = coord(:) + 1

	write(1000+rank_world,'(3i6,a10,3i5)'),rank_world, realm, rank_realm,'',&
	                                       coord(1), coord(2), coord(3)
	
	if (myid_world.eq.nproc_world) then
		write(1000+rank_world,*), '------------ END REALM INFO ---------------'
		write(1000+rank_world,*), '==========================================='
	end if
	
end subroutine write_realm_info

subroutine write_overlap_comms_md
	use coupler_module
	use mpi
	implicit none

	integer :: coord(3)

	if (myid_realm.eq.0) then
		write(2000+rank_realm,*),'rank_realm,rank_olap,   mdcoord,'  &
		                        ,'   overlapgroup' 
	end if
	
	call MPI_cart_coords(CPL_CART_COMM,myid_cart,3,coord,ierr)
	coord(:) = coord(:) + 1

	write(2000+rank_realm,'(2i7,a5,3i5,a5,i10)'), &
		rank_realm,rank_olap,'',coord,'',testval 	

end subroutine write_overlap_comms_md

! ++ UNINTERESTING ++ ========================================================
subroutine initialise
	use mpi
	use coupler_module
	implicit none
	
	call MPI_init(ierr)
	CPL_WORLD_COMM = MPI_COMM_WORLD

	call MPI_comm_size(CPL_WORLD_COMM, nproc_world, ierr)

	call MPI_comm_rank(CPL_WORLD_COMM, myid_world, ierr)	
	rank_world = myid_world + 1

end subroutine initialise

subroutine finalise
	use mpi
	implicit none
	
	integer :: ierr
	call MPI_finalize(ierr)

end subroutine finalise

subroutine barrier
	use mpi
	implicit none
	
	integer :: ierr	
	call MPI_barrier(MPI_COMM_WORLD,ierr)

end subroutine barrier

subroutine lasterrorcheck
	use mpi 
	implicit none
	
	integer :: ierr
	integer :: resultlen
	character*12 err_buffer

	call MPI_Error_string(ierr,err_buffer,resultlen,ierr)
	print*, err_buffer

end subroutine lasterrorcheck
