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
		iblock_realm=rank2coord_md(1,rank_realm)
		jblock_realm=rank2coord_md(2,rank_realm)
		kblock_realm=rank2coord_md(3,rank_realm)	
	else if (realm .eq. cfd_realm) then
		coord2rank_cfd(coord(1),coord(2),coord(3)) = rank_realm
		rank2coord_cfd(:,rank_realm)    = coord(:)
		rank_cfd2rank_world(rank_realm) = rank_world
		iblock_realm=rank2coord_cfd(1,rank_realm)
		jblock_realm=rank2coord_cfd(2,rank_realm)
		kblock_realm=rank2coord_cfd(3,rank_realm)	
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

	integer, dimension(:), allocatable :: mbuf, cbuf

	allocate(mbuf(3*nproc_md))
	allocate(cbuf(3*nproc_cfd))

	call MPI_allreduce(rank2coord_md,mbuf,3*nproc_md,MPI_INTEGER,MPI_SUM,   &
	                   CPL_WORLD_COMM,ierr)	
	call MPI_allreduce(rank2coord_cfd,cbuf,3*nproc_cfd,MPI_INTEGER,MPI_SUM, &
	                   CPL_WORLD_COMM,ierr)	

	rank2coord_md  = reshape(mbuf,(/3,nproc_md/))                  
	rank2coord_cfd = reshape(cbuf,(/3,nproc_cfd/))

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
	!split comm to groups
	!if group(world_rank) == 0, set olap_comm to null 

	integer :: i,j,k,ic,jc,kc
	integer :: trank_md, trank_cfd, trank_world
	integer, dimension(:), allocatable :: mdicoords, mdjcoords, mdkcoords
	integer, parameter :: olap_null = -666
	integer :: group(nproc_world)
	integer :: cfdcoord(3)

	allocate(mdicoords(npx_md/npx_cfd))
	allocate(mdjcoords(npy_md/npy_cfd))
	allocate(mdkcoords(npz_md/npz_cfd))

	!Set default values, must be done because coord2rank_md cannot
	!take "null" coordinates.
	group(:) = olap_null
	olap_mask(:) = 0

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

			trank_world = rank_cfd2rank_world(trank_cfd)
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
			trank_world = rank_md2rank_world(trank_md)

			olap_mask(trank_world) = 1
			group    (trank_world) = trank_cfd

		end do
		end do	
		end do

	end do

	! Split world Comm into a set of comms for overlapping processors
	call MPI_comm_split(CPL_WORLD_COMM,group(rank_world),realm, &
	                    CPL_OLAP_COMM,ierr)

	call MPI_comm_size(CPL_OLAP_COMM,nproc_olap,ierr)
	call MPI_comm_rank(CPL_OLAP_COMM,myid_olap,ierr)
	rank_olap = myid_olap + 1	
	if (myid_olap .eq. 0) testval = group(rank_world)
	call MPI_bcast(testval,1,MPI_INTEGER,0,CPL_OLAP_COMM,ierr)

	! Set all non-overlapping processors to MPI_COMM_NULL
	if (olap_mask(rank_world).eq.0) then
		myid_olap = olap_null
		rank_olap = olap_null
		CPL_OLAP_COMM = MPI_COMM_NULL
	end if

	deallocate(mdicoords)
	deallocate(mdjcoords)
	deallocate(mdkcoords)

	if (realm.eq.md_realm) call write_overlap_comms_md

end subroutine prepare_overlap_comms



!Setup topology graph of overlaps between CFD & MD processors

subroutine CPL_overlap_topology
	use coupler_module
	use mpi
	implicit none

	integer								:: i, n, nneighbors, nconnections
	integer, dimension(:),allocatable	:: index, edges, neighbors
	logical								:: reorder

	!Allow optimisations of ordering
	reorder = .true.

	!Get number of processors in communicating overlap region 
	!call MPI_comm_rank(CPL_OLAP_COMM, myid_graph, ierr)
	!if (CPL_OLAP_COMM .ne. MPI_COMM_NULL) then
	if (olap_mask(rank_world).eq.1) then

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
		call MPI_comm_rank(CPL_GRAPH_COMM,myid_graph,ierr)

		! - - TEST - - TEST - -
		!Get number of neighbours
		call MPI_comm_rank( CPL_GRAPH_COMM, myid_graph, ierr)
		call MPI_Graph_neighbors_count( CPL_GRAPH_COMM, myid_graph, nneighbors, ierr)
		allocate(neighbors(nneighbors))
		!Get neighbours
		call MPI_Graph_neighbors( CPL_GRAPH_COMM, myid_graph, nneighbors, neighbors,ierr )
		select case(realm)
		case(cfd_realm)
			write(3000+myid_world,*), code_name(realm),' My graph',myid_world,myid_graph,myid_olap,rank2coord_cfd(:,rank_realm), nneighbors, neighbors
		case(md_realm)
			write(3000+myid_world,*), code_name(realm),' My graph',myid_world,myid_graph,myid_olap,rank2coord_md(:,rank_realm), nneighbors, neighbors
		end select
		! - - TEST - - TEST - -


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


module coupler


    interface coupler_send
        module procedure coupler_send_3d, coupler_send_4d
    end interface

    interface coupler_recv
        module procedure coupler_recv_3d, coupler_recv_4d
    end interface

    private coupler_send_3d, coupler_send_4d, &
        coupler_send_xd, coupler_recv_3d, coupler_recv_4d,&
        coupler_recv_xd

contains


!=============================================================================
! coupler_send_data wrapper for 3d arrays
! see coupler_send_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_send_3d(temp,index_transpose)
    implicit none
 
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
  
    call coupler_send_xd(asend,index_transpose)

end subroutine coupler_send_3d

!=============================================================================
! coupler_send_data wrapper for 4d arrays
! see coupler_send_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_send_4d(asend,index_transpose)
    implicit none
 
	integer, optional, intent(in)  :: index_transpose(3)   
	real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in) :: asend
    
    integer n1,n2,n3,n4
 
    call coupler_send_xd(asend,index_transpose)

end subroutine coupler_send_4d


!=============================================================================
! Send data from the local grid to the associated ranks from the other 
! realm
!-----------------------------------------------------------------------------

subroutine coupler_send_xd(asend,index_transpose)
	use mpi
	use coupler_module
	implicit none
 
    ! array containing data distributed on the grid
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in) :: asend
    ! specify the order of dimensions in asend default is ( x, y, z) 
	! but some CFD solvers use (z,x,y)     
    integer, optional, intent(in) :: index_transpose(3)    ! rule of transposition for the coordinates 

	!Number of halos
	integer	:: nh = 1    

	!Neighbours
	integer								:: nneighbors   
	integer,dimension(:),allocatable	:: neighbors

    ! local indices 
    integer :: ix,iy,iz,bicmin,bicmax,bjcmin,bjcmax,bkcmin,bkcmax
	integer	:: ig(2,3),a_components

    ! auxiliaries 
    integer	:: i,ndata,itag,dest
	real(kind=kind(0.d0)), allocatable :: vbuf(:)
	integer, save :: ncalls = 0

	! This local CFD domain is outside MD overlap zone 
	!if (CPL_OLAP_COMM .eq. MPI_COMM_NULL) return
	if (olap_mask(rank_world).eq.0) return

	ncalls = ncalls + 1
    ! Get local grid box ranges seen by this rank for either CFD or MD
    if (realm .eq. cfd_realm) then 
		!Load CFD cells per processor
        bicmin = icPmin_cfd(iblock_realm)
        bicmax = icPmax_cfd(iblock_realm) 
        bjcmin = 1 !jcPmin_cfd(jblock_realm) 
        bjcmax = 2 !jcPmax_cfd(jblock_realm) 
        bkcmin = kcPmin_cfd(kblock_realm) 
        bkcmax = kcPmax_cfd(kblock_realm) 
    elseif (realm .eq. md_realm) then 
        ! Load MD cells per processor
        bicmin = icPmin_md(iblock_realm)
        bicmax = icPmax_md(iblock_realm) 
        bjcmin = 1 !jcPmin_md(jblock_realm) 
        bjcmax = 2 !jcPmax_md(jblock_realm) 
        bkcmin = kcPmin_md(kblock_realm) 
        bkcmax = kcPmax_md(kblock_realm) 
	endif

    ! Re-order indices of x,y,z coordinates (handles transpose arrays)
    ix = 1 ; iy = 2; iz = 3
    if( present(index_transpose)) then
        ix = index_transpose(1)
        iy = index_transpose(2)
        iz = index_transpose(3)
    endif

    ! Number of components at each grid point
    a_components = size(asend,1)

	!Get neighbours
	call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
	allocate(neighbors(nneighbors))
	call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,neighbors,ierr )
   
    ! loop through the maps and send the corresponding sections of asend
    do i = 1, nneighbors

		!Get taget processor from mapping
        dest = neighbors(i)

        ! Amount of data to be sent
        ndata = a_components * (bicmax-bicmin) * (bjcmax-bjcmin) * (bkcmax-bkcmin)
		if (allocated(vbuf)) deallocate(vbuf); allocate(vbuf(ndata))
		!print*, 'seg fault', ndata,size(asend,1),size(asend,2),size(asend,3),size(asend,4),size(asend), &
 		!		a_components,bicmax,bicmin,bjcmax,bjcmin,bkcmax,bkcmin,(bicmax-bicmin),(bjcmax-bjcmin),&
		!		(bkcmax-bkcmin)
		vbuf(1:ndata) = reshape(asend, (/ ndata /))

        ! Send data 
        itag = 0 !mod( ncalls, MPI_TAG_UB) !Attention ncall could go over max tag value for long runs!!
		call MPI_send(vbuf, ndata, MPI_DOUBLE_PRECISION, dest, itag, CPL_GRAPH_COMM, ierr)
    enddo

end subroutine coupler_send_xd

!=============================================================================
! coupler_recv_data wrapper for 3d arrays
! see coupler_recv_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_recv_3d(temp,index_transpose)
    implicit none

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
    call coupler_recv_xd(arecv,index_transpose)
	temp(:,:,:) = 	arecv(1,:,:,:) 

end subroutine coupler_recv_3d

!=============================================================================
! coupler_recv_data wrapper for 4d arrays
! see coupler_recv_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_recv_4d(arecv,index_transpose)

    implicit none

    integer, optional, intent(in) :: index_transpose(3)  
    real(kind(0.d0)),dimension(:,:,:,:),intent(inout) :: arecv

    call coupler_recv_xd(arecv,index_transpose)

end subroutine coupler_recv_4d

!=============================================================================
! Receive data from to local grid from the associated ranks from the other 
! realm
!-----------------------------------------------------------------------------
subroutine coupler_recv_xd(arecv,index_transpose)
    use mpi
	use coupler_module
    implicit none

    ! specify the order of dimensions in asend default is ( x, y, z) 
	! but some CFD solvers use (z,x,y)   
    integer, optional, intent(in) :: index_transpose(3) 
     
    ! Array that recieves grid distributed data 
    real(kind(0.d0)), dimension(:,:,:,:),intent(inout) :: arecv     

	!Neighbours
	integer								:: nneighbors   
	integer,dimension(:),allocatable	:: neighbors
                                                         
    ! local indices 
    integer :: ix,iy,iz,bicmin,bicmax,bjcmin,bjcmax,bkcmin,bkcmax
	integer	:: i,ig(2,3),pcoords(3),recvdata,startbuf,endbuf,a_components

    ! auxiliaries 
    integer	:: ndata, itag, source,start_address
	integer,dimension(:),allocatable   :: req
	integer,dimension(:,:),allocatable :: vel_indx, status
    real(kind(0.d0)),dimension(:), allocatable ::  vbuf,vbuf2
 
	! This local CFD domain is outside MD overlap zone 
	if (CPL_OLAP_COMM .eq. MPI_COMM_NULL) return

    ! Local grid box
    ! Get local grid box ranges seen by this rank for either CFD or MD
    if (realm .eq. cfd_realm) then 
		!Load CFD cells per processor
        bicmin = icPmin_cfd(iblock_realm)
        bicmax = icPmax_cfd(iblock_realm) 
        bjcmin = 1 !jcPmin_cfd(jblock_realm) 
        bjcmax = 2 !jcPmax_cfd(jblock_realm) 
        bkcmin = kcPmin_cfd(kblock_realm) 
        bkcmax = kcPmax_cfd(kblock_realm)
    elseif (realm .eq. md_realm) then 
        ! Load MD cells per processor
        bicmin = icPmin_md(iblock_realm)
        bicmax = icPmax_md(iblock_realm) 
        bjcmin = 1 !jcPmin_md(jblock_realm) 
        bjcmax = 2 !jcPmax_md(jblock_realm) 
        bkcmin = kcPmin_md(kblock_realm) 
        bkcmax = kcPmax_md(kblock_realm) 
	endif

    ! Get the indices in x,y,z direction from transpose array
    ix = 1 ; iy = 2; iz = 3
    if( present(index_transpose)) then
        ix = index_transpose(1)
        iy = index_transpose(2)
        iz = index_transpose(3)
    endif
	ig(1,ix) = bicmin;	ig(2,ix) = bicmax
	ig(1,iy) = bjcmin;	ig(2,iy) = bjcmax
	ig(1,iz) = bkcmin;	ig(2,iz) = bkcmax

    ! Amount of data to receive
	a_components = size(arecv,1)
    ndata = a_components * (bicmax-bicmin) * (bjcmax-bjcmin) * (bkcmax-bkcmin)
	!print*, 'vbuf size', ndata , a_components , (bicmax-bicmin + 1) , (bjcmax-bjcmin + 1) , (bkcmax-bkcmin + 1)
	allocate(vbuf(ndata),stat=ierr) 

	!Get neighbours
	call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
	allocate(neighbors(nneighbors))
	call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,neighbors,ierr )

    ! Receive from all attached processors
	allocate(req(nneighbors))
	allocate(vel_indx(8,nneighbors), status(MPI_STATUS_SIZE,nneighbors))
    start_address = 1 

    do i = 1, nneighbors

		!Get source processor from mapping
        source =  neighbors(i)

	    ! Get size of data to receive from source processors
		if (realm .eq. cfd_realm) then
			pcoords(1)=rank2coord_md(1,source)
			pcoords(2)=rank2coord_md(2,source)
			pcoords(3)=rank2coord_md(3,source)
			recvdata = a_components * (icPmax_md(pcoords(1))-icPmin_md(pcoords(1))) & 
									* (jcPmax_md(pcoords(2))-jcPmin_md(pcoords(2))) & 
									* (kcPmax_md(pcoords(3))-kcPmin_md(pcoords(3)))
		elseif (realm .eq. md_realm) then
			pcoords(1)=rank2coord_cfd(1,source)
			pcoords(2)=rank2coord_cfd(2,source)
			pcoords(3)=rank2coord_cfd(3,source)
			recvdata = a_components * (icPmax_cfd(pcoords(1))-icPmin_cfd(pcoords(1))) & 
									* (jcPmax_cfd(pcoords(2))-jcPmin_cfd(pcoords(2))) & 
									* (kcPmax_cfd(pcoords(3))-kcPmin_cfd(pcoords(3)))
		endif

		if (recvdata .gt. ndata) then
			! If data received is greater than required, discard excess
			allocate(vbuf2(recvdata))
			itag = 0 !mod(ncalls, MPI_TAG_UB) ! Attention ncall could go over max tag value for long runs!!
			call MPI_irecv(vbuf2(:),recvdata,MPI_DOUBLE_PRECISION,source,itag,&
                				CPL_GRAPH_COMM,req(i),ierr)
			startbuf = neighbors(rank_realm) * ndata
			endbuf   = neighbors(rank_realm) * (ndata + 1 )
			vbuf(1:ndata) = vbuf2(startbuf:endbuf)
			deallocate(vbuf2)
		else
			! Otherwise Receive section of data and increment pointer 
			! ready to receive next piece of data
	        start_address = start_address + recvdata
			itag = 0 !mod(ncalls, MPI_TAG_UB) ! Attention ncall could go over max tag value for long runs!!
			call MPI_irecv(vbuf(start_address),recvdata,MPI_DOUBLE_PRECISION,source,itag,&
                				CPL_GRAPH_COMM,req(i),ierr)
		endif
    enddo
    call MPI_waitall(nneighbors, req, status, ierr)

	print*, 'RECV reshape', ig(2,1),ig(1,1), ig(2,2),ig(1,2), ig(2,3),ig(1,3) , & 
			ig(2,1)-ig(1,1)+1, ig(2,2)-ig(1,2)+1, ig(2,3)-ig(1,3)+1 
	arecv = reshape(vbuf,(/ ndata, ig(2,1)-ig(1,1), ig(2,2)-ig(1,2), ig(2,3)-ig(1,3) /))
           
end subroutine coupler_recv_xd

end module coupler


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
		write(2000+rank_realm,*),'rank_realm,rank_olap,  mdcoord,'  &
		                        ,'   overlapgroup,   olap_mask,  CPL_OLAP_COMM'
	end if
	
	call MPI_cart_coords(CPL_CART_COMM,myid_cart,3,coord,ierr)
	coord(:) = coord(:) + 1

	write(2000+rank_realm,'(2i7,a5,3i5,a5,2i10,a5,i)'), &
		rank_realm,rank_olap,'',coord,'',testval,olap_mask(rank_world), &
		'',CPL_OLAP_COMM

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
