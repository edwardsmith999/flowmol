
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
	integer	:: a_components

    ! auxiliaries 
    integer	:: i,ndata,itag,destid
	real(kind=kind(0.d0)), allocatable :: vbuf(:)
	integer, save :: ncalls = 0

	! This local CFD domain is outside MD overlap zone 
	!if (CPL_OLAP_COMM .eq. MPI_COMM_NULL) return
	if (olap_mask(rank_world) .eq. 0) return

    ! Get local grid box ranges seen by this rank for either CFD or MD
    if (realm .eq. cfd_realm) then 
		!Load CFD cells per processor
        bicmin = icPmin_cfd(iblock_realm)
        bicmax = icPmax_cfd(iblock_realm) 
        bjcmin = 1 !jcPmin_cfd(jblock_realm) 
        bjcmax = 1 !jcPmax_cfd(jblock_realm) 
        bkcmin = kcPmin_cfd(kblock_realm) 
        bkcmax = kcPmax_cfd(kblock_realm) 
    elseif (realm .eq. md_realm) then 
        ! Load MD cells per processor
        bicmin = icPmin_md(iblock_realm)
        bicmax = icPmax_md(iblock_realm) 
        bjcmin = 1 !jcPmin_md(jblock_realm) 
        bjcmax = 1 !jcPmax_md(jblock_realm) 
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

	!print*, 'send data',rank_world,rank_realm,rank_olap,olap_mask(rank_world),bicmax,bjcmax,bjcmin,bkcmax,bkcmin
   
    ! loop through the maps and send the corresponding sections of asend
    do i = 1, nneighbors

		!Get taget processor from mapping
        destid = neighbors(i)

        ! Amount of data to be sent
        ndata = a_components * (bicmax-bicmin+1) * (bjcmax-bjcmin+1) * (bkcmax-bkcmin+1)
		if (allocated(vbuf)) deallocate(vbuf); allocate(vbuf(ndata))
		vbuf(1:ndata) = reshape(asend, (/ ndata /))
		!print'(a,8i8,f20.5)', 'send data',rank_world,rank_realm,rank_olap,ndata, & 
		!						size(asend),iblock_realm,jblock_realm,kblock_realm,vbuf(10)

        ! Send data 
        itag = 0 !mod( ncalls, MPI_TAG_UB) !Attention ncall could go over max tag value for long runs!!
		call MPI_send(vbuf, ndata, MPI_DOUBLE_PRECISION, destid, itag, CPL_GRAPH_COMM, ierr)
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
    integer :: n,i,j,k,ix,iy,iz,bicmin,bicmax,bjcmin,bjcmax,bkcmin,bkcmax
	integer	:: ncl(2,3),pcoords(3),recvsize,startbuf,endbuf,a_components

    ! auxiliaries 
    integer	:: ndata, itag, sourceid,source_realm,start_address
	integer,dimension(:),allocatable   :: req
	integer,dimension(:,:),allocatable :: status
	integer,dimension(:,:,:),allocatable :: ncl_recv
    real(kind(0.d0)),dimension(:), allocatable ::  vbuf,vbuf2
 
	! This local CFD domain is outside MD overlap zone 
	if (olap_mask(rank_world).eq.0) return
	!if (CPL_OLAP_COMM .eq. MPI_COMM_NULL) return

    ! Local grid box
    ! Get local grid box ranges seen by this rank for either CFD or MD
    if (realm .eq. cfd_realm) then 
		!Load CFD cells per processor
        bicmin = icPmin_cfd(iblock_realm)
        bicmax = icPmax_cfd(iblock_realm) 
        bjcmin = 1 !jcPmin_cfd(jblock_realm) 
        bjcmax = 1 !jcPmax_cfd(jblock_realm) 
        bkcmin = kcPmin_cfd(kblock_realm) 
        bkcmax = kcPmax_cfd(kblock_realm)
    elseif (realm .eq. md_realm) then 
        ! Load MD cells per processor
        bicmin = icPmin_md(iblock_realm)
        bicmax = icPmax_md(iblock_realm) 
        bjcmin = 1 !jcPmin_md(jblock_realm) 
        bjcmax = 1 !jcPmax_md(jblock_realm) 
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
	ncl(1,ix) = bicmin;	ncl(2,ix) = bicmax
	ncl(1,iy) = bjcmin;	ncl(2,iy) = bjcmax
	ncl(1,iz) = bkcmin;	ncl(2,iz) = bkcmax

    ! Amount of data to receive
	a_components = size(arecv,1)
    ndata = a_components * (bicmax-bicmin+1) * (bjcmax-bjcmin+1) * (bkcmax-bkcmin+1)
	allocate(vbuf(ndata),stat=ierr); vbuf = 0.d0

	if (size(arecv) .ne. ndata) stop "domain size mismatch in recv data"

	!Get neighbours
	call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
	allocate(neighbors(nneighbors))
	call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,neighbors,ierr )

    ! Receive from all attached processors
	allocate(req(nneighbors))
	allocate(status(MPI_STATUS_SIZE,nneighbors))
	allocate(ncl_recv(2,3,nneighbors))
    start_address = 1 
    do i = 1, nneighbors

		!Get source processor from topology graph
        sourceid =  neighbors(i)

	    ! Get size of data to receive from source processors
		if (realm .eq. cfd_realm) then
			source_realm = rank_olap2rank_realm(sourceid+1)
			pcoords(1)=rank2coord_md(1,source_realm)
			pcoords(2)=rank2coord_md(2,source_realm)
			pcoords(3)=rank2coord_md(3,source_realm)
			ncl_recv(1,ix,sourceid) = icPmin_md(pcoords(1)); ncl_recv(2,ix,sourceid) = icPmax_md(pcoords(1))
			ncl_recv(1,iy,sourceid) = jcPmin_md(pcoords(2)); ncl_recv(2,iy,sourceid) = jcPmax_md(pcoords(2))
			ncl_recv(1,iz,sourceid) = kcPmin_md(pcoords(3)); ncl_recv(2,iz,sourceid) = kcPmax_md(pcoords(3))
			recvsize = a_components * (icPmax_md(pcoords(1))-icPmin_md(pcoords(1))+1) & 
									* 1 & !(jcPmax_md(pcoords(2))-jcPmin_md(pcoords(2))) & 
									* (kcPmax_md(pcoords(3))-kcPmin_md(pcoords(3))+1)
			!print'(a,12i8)', 'cfd source data',realm,sourceid,source_realm, & 
			!				olap_mask(rank_mdcart2rank_world(source_realm)),pcoords,recvsize,a_components, & 
			!				icPmax_md(pcoords(1))-icPmin_md(pcoords(1))+1,jcPmax_md(pcoords(2))- & 
			!				jcPmin_md(pcoords(2))+1,kcPmax_md(pcoords(3))-kcPmin_md(pcoords(3))+1
		elseif (realm .eq. md_realm) then
			source_realm = rank_olap2rank_realm(sourceid)
			pcoords(1)=rank2coord_cfd(1,source_realm)
			pcoords(2)=rank2coord_cfd(2,source_realm)
			pcoords(3)=rank2coord_cfd(3,source_realm)
			recvsize = a_components * (icPmax_cfd(pcoords(1))-icPmin_cfd(pcoords(1))) & 
									* 1 & !(jcPmax_cfd(pcoords(2))-jcPmin_cfd(pcoords(2))) & 
									* (kcPmax_cfd(pcoords(3))-kcPmin_cfd(pcoords(3)))
			!print'(a,14i8)', 'md  source data',realm,sourceid,source_realm, & 
			!				olap_mask(rank_cfdcart2rank_world(source_realm)),pcoords,recvsize, & 
			!				icPmax_md(pcoords(1)),icPmin_md(pcoords(1)),jcPmax_md(pcoords(2)), & 
			!				jcPmin_md(pcoords(2)),kcPmax_md(pcoords(3)),kcPmin_md(pcoords(3))
		endif

		if (recvsize .gt. ndata) then
			! If data received is greater than required, discard excess
			allocate(vbuf2(recvsize))
			itag = 0 !mod(ncalls, MPI_TAG_UB) ! Attention ncall could go over max tag value for long runs!!
			call MPI_irecv(vbuf2(:),recvsize,MPI_DOUBLE_PRECISION,sourceid,itag,&
                				CPL_GRAPH_COMM,req(i),ierr)
			startbuf = neighbors(rank_realm) * ndata
			endbuf   = neighbors(rank_realm) *(ndata + 1 )
			vbuf(1:ndata) = vbuf2(startbuf:endbuf)
			deallocate(vbuf2)
		else
			! Otherwise Receive section of data and increment pointer 
			! ready to receive next piece of data
			itag = 0 !mod(ncalls, MPI_TAG_UB) ! Attention ncall could go over max tag value for long runs!!
			start_address = 1+((pcoords(1)-1) + (pcoords(3)-1)*npx_md)*recvsize
			!print'(a,10i8)', 'recv data',realm,sourceid+1,source_realm,ndata,recvsize,size(arecv),start_address,pcoords
			!start_address=(sourceid-1)*recvsize+1
			call MPI_irecv(vbuf(start_address),recvsize,MPI_DOUBLE_PRECISION,sourceid,itag,&
                				CPL_GRAPH_COMM,req(i),ierr)

	        !start_address = start_address + recvsize

		endif
    enddo
    call MPI_waitall(nneighbors, req, status, ierr)

	!do i=1, nneighbors
	!	arecv(:,ncl_recv(1,1,i):ncl_recv(2,1,i), & 
	!			1:1, & 
	!			ncl_recv(1,3,i):ncl_recv(2,3,i)) = reshape(vbuf,(/ 3,ncl_recv(1,1,i)-ncl_recv(2,1,i), & 
	!															 	 1, & 
	!															 	 ncl_recv(1,3,i)-ncl_recv(2,3,i)/))
	!enddo

	!arecv = reshape(vbuf,(/ 3, ncl(2,1)-ncl(1,1), ncl(2,2)-ncl(1,2), ncl(2,3)-ncl(1,3) /))

	!do n = 1,size(vbuf),recvsize
	!	print*, 'vbuf', n, vbuf(n)
	!enddo

	!do n = 1,size(arecv,1)
	!do i = 1,size(arecv,2)
	!do j = 1,size(arecv,3)
	!do k = 1,size(arecv,4)
	!	print*, 'arecv', n,i,j,k, arecv(n,i,j,k)
	!enddo
	!enddo
	!enddo
	!enddo
           
end subroutine coupler_recv_xd


!-------------------------------------------------------------------
! 					CPL_comm_map								   -
!-------------------------------------------------------------------

! Get COMM map for current communicator and relationship to 
! world rank used to link to others in the coupler hierachy

! - - - Synopsis - - -

! CPL_comm_rank(COMM, rank, comm2world, world2comm, ierr)

! - - - Input Parameters - - -

!comm
!    communicator with cartesian structure (handle) 
!rank
!    rank of a process within group of comm (integer) 

! - - - Output Parameter - - -

!comm2world
!	Array of size nproc_world which for element at 
!	world_rank has local rank in COMM
!world2comm
!	Array of size nproc_COMM which for element at 
!	for local ranks in COMM has world rank 
!ierr
!    error flag


subroutine CPL_comm_map(COMM,rank,nproc,comm2world,world2comm,ierr)
	use coupler_module, only : rank_world, nproc_world, CPL_WORLD_COMM, VOID
	use mpi
	implicit none

	integer, intent(in)								:: COMM
	integer, intent(out)							:: rank,nproc,ierr
	integer, dimension(:),allocatable,intent(out)	:: comm2world,world2comm

	!Mapping from world rank to comm rank
	allocate(world2comm( nproc_world))
	world2comm( nproc_world) = VOID

	if (COMM .ne. MPI_COMM_NULL) then

		!Mapping from comm rank to world rank
		call MPI_comm_rank(COMM,rank,ierr)
		call MPI_comm_size(COMM,nproc,ierr)
		allocate(comm2world(nproc))

		!print*, 'Inside CPL_comm_map', rank, nproc

		call MPI_allgather(rank_world,1,MPI_INTEGER, & 
						   comm2world,1,MPI_INTEGER,COMM,ierr)
		!print'(a,11i8)', 'after comm2world in CPL_comm_map', rank_world,rank+1,comm2world

	else
		rank = VOID
	endif

	call MPI_allgather(rank      ,1,MPI_INTEGER, & 
					   world2comm,1,MPI_INTEGER,CPL_WORLD_COMM,ierr)
	!if (rank .ne. VOID) print*, 'after world2comm in CPL_comm_map', rank_world,rank+1,world2comm(rank+1)

end subroutine CPL_comm_map


end module coupler



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
		coord2rank_md(coord(1),coord(2),coord(3)) = rank_cart
		rank2coord_md(:,rank_cart)    = coord(:)
		rank_mdcart2rank_world(rank_cart) = rank_world
		iblock_realm=rank2coord_md(1,rank_cart)
		jblock_realm=rank2coord_md(2,rank_cart)
		kblock_realm=rank2coord_md(3,rank_cart)	
	else if (realm .eq. cfd_realm) then
		coord2rank_cfd(coord(1),coord(2),coord(3)) = rank_cart
		rank2coord_cfd(:,rank_cart)    = coord(:)
		rank_cfdcart2rank_world(rank_cart) = rank_world
		iblock_realm=rank2coord_cfd(1,rank_cart)
		jblock_realm=rank2coord_cfd(2,rank_cart)
		kblock_realm=rank2coord_cfd(3,rank_cart)	
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
	
	integer							   :: buf
	integer, dimension(:), allocatable :: mbuf, cbuf

	! Realm to world
	allocate(mbuf(nproc_md))
	allocate(cbuf(nproc_cfd))

	call MPI_allreduce(rank_mdcart2rank_world,mbuf,nproc_md,MPI_INTEGER,  &
	                   MPI_SUM, CPL_WORLD_COMM,ierr)	
	call MPI_allreduce(rank_cfdcart2rank_world,cbuf,nproc_cfd,MPI_INTEGER, &
	                   MPI_SUM, CPL_WORLD_COMM,ierr)	

	rank_mdcart2rank_world  = mbuf
	rank_cfdcart2rank_world = cbuf

	deallocate(mbuf)
	deallocate(cbuf)

	allocate(rank_world2rank_cfdcart(nproc_world))
	allocate(rank_world2rank_mdcart(nproc_world))
	call MPI_allgather(rank_cart		      ,1,MPI_INTEGER, & 
					   rank_world2rank_cfdcart,1,MPI_INTEGER,CPL_WORLD_COMM,ierr)
	call MPI_allgather(rank_cart		      ,1,MPI_INTEGER, & 
					   rank_world2rank_mdcart, 1,MPI_INTEGER,CPL_WORLD_COMM,ierr)
	!Collect on own realm
	!if (realm .eq. md_realm) then
	!	call CPL_comm_map(CPL_CART_COMM,rank_cart,nproc_cart, & 
	!					 rank_mdcart2rank_world,rank_world2rank_mdcart,ierr)
	!else if (realm .eq. cfd_realm) then
	!	call CPL_comm_map(CPL_CART_COMM,rank_cart,nproc_cart, & 
	!					 rank_cfdcart2rank_world,rank_world2rank_cfdcart,ierr)
	!endif

	!Exchange across intercomm
	!call bcast(icomm

end subroutine collect_rank2ranks

subroutine collect_rank2ranks_olap
	use mpi
	use coupler_module
	implicit none
	
	integer :: buf

	buf = rank_world
	call MPI_allgather(        buf          ,1,MPI_INTEGER, & 
						rank_olap2rank_world,1,MPI_INTEGER,CPL_OLAP_COMM,ierr)
	buf = rank_realm
	call MPI_allgather(        buf          ,1,MPI_INTEGER, & 
						rank_olap2rank_realm,1,MPI_INTEGER,CPL_OLAP_COMM,ierr)

end subroutine collect_rank2ranks_olap

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
	call MPI_comm_size(CPL_OLAP_COMM,nproc_olap,ierr)
	call MPI_comm_rank(CPL_OLAP_COMM,myid_olap,ierr)
	rank_olap = myid_olap + 1	

	!Setup rank_olap_2_rank_realm array
	allocate(rank_olap2rank_realm(nproc_olap))
	allocate(rank_olap2rank_world(nproc_olap))
	call collect_rank2ranks_olap

	!Store rank_olap_2_rank_realm from temp counting array
	!allocate(rank_olap2rank_realm(nolap))
	!rank_olap2rank_realm = rank_olap2rank_realm_temp(1:nolap)

	! USED ONLY FOR OUTPUT/TESTING??
	if (myid_olap .eq. 0) testval = group(rank_world)
	call MPI_bcast(testval,1,MPI_INTEGER,0,CPL_OLAP_COMM,ierr)

	! Set all non-overlapping processors to MPI_COMM_NULL
	if (olap_mask(rank_world).eq.0) then
		myid_olap = olap_null
		rank_olap = olap_null
		CPL_OLAP_COMM = MPI_COMM_NULL
	end if

	!call CPL_comm_map(CPL_OLAP_COMM,rank_olap,nproc_olap, & 
	!				 rank_olap2rank_world,rank_world2rank_olap,ierr)
	!myid_olap = rank_olap - 1

	deallocate(mdicoords)
	deallocate(mdjcoords)
	deallocate(mdkcoords)
	
	if (realm.eq.md_realm) call write_overlap_comms_md

end subroutine prepare_overlap_comms



!Setup topology graph of overlaps between CFD & MD processors

subroutine CPL_overlap_topology
	use coupler_module
	use coupler, only : CPL_comm_map
	use mpi
	implicit none

	integer								:: i, n, nneighbors, nconnections
	integer, dimension(:),allocatable	:: index, edges, neighbors
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


		!call MPI_comm_rank(CPL_GRAPH_COMM,myid_graph,ierr)

		!print'(a,11i8)', 'graph',realm,rank_world, iblock_realm,jblock_realm,kblock_realm,& 
		!							icPmin_md(iblock_realm),icPmax_md(iblock_realm), &
		!           				jcPmin_md(jblock_realm),jcPmax_md(jblock_realm), &
		!           				kcPmin_md(kblock_realm),kcPmax_md(kblock_realm)
   

        !   TEST   TEST  
        !Get number of neighbours
        !call MPI_comm_rank( CPL_GRAPH_COMM, myid_graph, ierr)
        call MPI_Graph_neighbors_count( CPL_GRAPH_COMM, myid_graph, nneighbors, ierr)
        allocate(neighbors(nneighbors))
        !Get neighbours
        call MPI_Graph_neighbors( CPL_GRAPH_COMM, myid_graph, nneighbors, neighbors,ierr )
        select case(realm)
        case(cfd_realm)
                write(3000+myid_world,*), realm_name(realm),' My graph', & 
								myid_world,myid_graph,myid_olap, & 
								rank2coord_cfd(:,rank_realm), nneighbors, neighbors
        case(md_realm)
                write(3000+myid_world,*), realm_name(realm),' My graph', & 
								myid_world,myid_graph,myid_olap, & 
								rank2coord_md(:,rank_realm), nneighbors, neighbors
        end select
        !   TEST   TEST  

	else
		CPL_GRAPH_COMM = MPI_COMM_NULL
	endif

	call CPL_comm_map(CPL_GRAPH_COMM,rank_graph,nproc_olap, & 
					 rank_graph2rank_world,rank_world2rank_graph,ierr)
	myid_graph = rank_graph - 1

	!print*, 'graph comm', rank_world,rank_graph,nproc_olap,rank_graph2rank_world,rank_world2rank_graph

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
