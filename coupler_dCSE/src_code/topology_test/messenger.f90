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
	call print_realm_info

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

	group = -666

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

	call MPI_comm_split(CPL_WORLD_COMM,group(rank_world),realm, &
	                    CPL_OLAP_COMM,ierr)

	call MPI_comm_rank(CPL_OLAP_COMM,myid_olap,ierr)
	rank_olap = myid_olap + 1	

	if (myid_olap .eq. 0) testval = group(rank_world)
	call MPI_bcast(testval,1,MPI_INTEGER,0,CPL_OLAP_COMM,ierr)

	deallocate(mdicoords)
	deallocate(mdjcoords)
	deallocate(mdkcoords)

!	integer :: trank_md, trank_cfd, trank_world, tid_md, tid_world
!	integer :: ncxl_cfd
!	integer :: ncyl_cfd
!	integer :: nczl_cfd
!	integer :: ncyP_md
!	integer :: ncy_md
!	integer :: ncy_olap
!	integer :: olap_jmin_mdcoord
!	integer, parameter :: coordnull = -666
!	integer :: mdcoord(3)
!	integer :: cfdcoord(3)
!	integer :: group(nproc_world)
!
!	ncxl_cfd = ncx / npx_cfd
!	ncyl_cfd = ncy / npy_cfd
!	nczl_cfd = ncz / npz_cfd

!	do trank_md = 1,nproc_md
!		
!		trank_world = rank_md2rank_world(trank_md)
!
!		tid_md = trank_md - 1
!
!		if (realm.eq.md_realm) then
!			call MPI_cart_coords(CPL_CART_COMM,tid_md,3,mdcoord,ierr)
!			call MPI_comm_rank(CPL_WORLD_COMM,tid_world,ierr)
!		end if
!
!		call MPI_bcast(mdcoord,3,MPI_INTEGER,tid_world,CPL_WORLD_COMM,ierr)
!		mdcoord = mdcoord + 1
!
!		!CALCULATE OLAP_JMIN_MDCOORD
!		mdcoord(2) = mdcoord(2) - (olap_jmin_mdcoord - 1) 
!
!		cfdcoord(1) = ceiling(mdcoord(1)*dble(npx_cfd)/dble(npx_md))	
!		cfdcoord(2) = ceiling(mdcoord(2)*dble(npy_cfd)/dble(npy_md))
!		cfdcoord(3) = ceiling(mdcoord(3)*dble(npz_cfd)/dble(npz_md))
!		trank_cfd   = coord2rank_cfd(cfdcoord(1),cfdcoord(2),cfdcoord(3))
!
!		if (olap_mask(trank_world) .eq. 1) then
!			group(trank_world) = trank_cfd
!		else
!			group(trank_world) = 0
!		end if
!
!	end do
!	
!	do trank_cfd = 1,nproc_cfd 
!
!		trank_world = rank_cfd2rank_world(trank_cfd)
!		if (olap_mask(trank_world) .eq. 1) then
!			group(trank_world) = trank_cfd
!		else
!			group(trank_world) = 0
!		end if
!
!	end do
!	
!	call MPI_comm_split(CPL_WORLD_COMM,group(rank_world),realm, &
!	                    CPL_OLAP_COMM,ierr)
!
!	if (olap_mask(rank_world).eq.1) then
!
!		call MPI_comm_rank(CPL_OLAP_COMM,myid_olap,ierr)
!		rank_olap = myid_olap + 1	
!		if (myid_olap .eq. 0) testval = group(rank_world)
!		call MPI_bcast(testval,1,MPI_INTEGER,0,CPL_OLAP_COMM,ierr)
!
!	else
!
!		CPL_OLAP_COMM = MPI_COMM_NULL
!
!	end if


	call print_overlap_comms



end subroutine prepare_overlap_comms

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


subroutine print_realm_info
	use coupler_module
	use mpi
	implicit none

	integer :: trank
	integer :: coord(3)

	call barrier

	if (myid_world.eq.0) then
		print*, '---------- REALM INFORMATION --------------'
		print*, ' wrank    realm    realmrank   cart coords '
		print*, '-------------------------------------------'
	end if

	call MPI_cart_coords(CPL_CART_COMM,myid_cart,3,coord,ierr)
	coord(:) = coord(:) + 1

	do trank = 1,nproc_world
		if (rank_world.eq.trank) then
			if (realm .eq. md_realm) then
				print('(i5,2x,a10,i8,4x,3i5)'), rank_world, 'md_realm',  &
				      rank_realm, coord(1), coord(2), coord(3)
			else if (realm .eq. cfd_realm) then
				print('(i5,2x,a10,i8,4x,3i5)'), rank_world, 'cfd_realm', &
				      rank_realm, coord(1), coord(2), coord(3)
			end if
		end if
		call barrier
	end do
	
	if (myid_world.eq.0) then
		print*, '------------ END REALM INFO ---------------'
		print*, '==========================================='
	end if
	
end subroutine print_realm_info

subroutine print_overlap_comms
	use coupler_module
	use mpi
	implicit none

	integer :: trank

	call barrier

	if (myid_world.eq.0) then
		print*, ''
		print*, '----------- OVERLAP COMMS INFO ------------'
		print*, '-------------------------------------------'
		print*, '        RANKS              BROADCAST TEST  '
		print*, '  world  realm  olap      testval( = group)'
		print*, '-------------------------------------------'
	end if
	
	do trank = 1,nproc_world
		if (rank_world.eq.trank) then
			print('(3i7,i16)'), rank_world,rank_realm, &
			                    rank_olap, testval 	
		end if
		call barrier
	end do

	if (myid_world.eq.0) then
		print*, '-------- END OVERLAP COMMS INFO  ----------'
		print*, '==========================================='
	end if
	
end subroutine print_overlap_comms

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
