program create_map
	use coupler_module, only : olap_mask, rank_world
	implicit none

	call initialise
	call setup_input_and_arrays
	call get_cell_ranges_md 
	call get_overlap_blocks_cfd 
	call create_realms
	call prepare_overlap_comms
	call CPL_overlap_topology
	!call test_COMMS
	!call test_packing

!	call test_send_recv_MD2CFD
!	call test_send_recv_CFD2MD
	if (olap_mask(rank_world).eq.1) call test_gather_scatter

	call finalise

end program create_map

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

subroutine setup_input_and_arrays
	use coupler_module
	implicit none

	open(unit=1,file='procs',form='formatted')
		read(1,*) npx_md 
		read(1,*) npy_md 
		read(1,*) npz_md 
		read(1,*) npx_cfd 
		read(1,*) npy_cfd 
		read(1,*) npz_cfd 
	close(1)

	nproc_md  = npx_md*npy_md*npz_md
	nproc_cfd = npx_cfd*npy_cfd*npz_cfd
	nproc_world = nproc_md + nproc_cfd

	allocate(icPmin_md(npx_md))
	allocate(icPmax_md(npx_md))
	allocate(jcPmin_md(npy_md))
	allocate(jcPmax_md(npy_md))
	allocate(kcPmin_md(npz_md))
	allocate(kcPmax_md(npz_md))
	allocate(icPmin_cfd(npx_cfd))
	allocate(icPmax_cfd(npx_cfd))
	allocate(jcPmin_cfd(npy_cfd))
	allocate(jcPmax_cfd(npy_cfd))
	allocate(kcPmin_cfd(npz_cfd))
	allocate(kcPmax_cfd(npz_cfd))
	allocate(cfd_icoord2olap_md_icoords(npx_cfd,npx_md/npx_cfd))
	allocate(cfd_jcoord2olap_md_jcoords(npy_cfd,npy_md/npy_cfd))
	allocate(cfd_kcoord2olap_md_kcoords(npz_cfd,npz_md/npz_cfd))

	allocate(olap_mask(nproc_world))

	allocate(coord2rank_md(npx_md,npy_md,npz_md))
	allocate(coord2rank_cfd(npx_cfd,npy_cfd,npz_cfd))
	allocate(rank2coord_cfd(3,nproc_cfd))
	allocate(rank2coord_md(3,nproc_md))
	allocate(rank_mdcart2rank_world(nproc_md))
	allocate(rank_cfdcart2rank_world(nproc_cfd))

	coord2rank_md       = 0
	coord2rank_cfd      = 0
	rank2coord_md       = 0
	rank2coord_cfd      = 0
	rank_mdcart2rank_world  = 0
	rank_cfdcart2rank_world = 0

	! Inputs
	ncx = 128
	ncy = 128
	ncz = 8	
	xL_cfd = 2048.d0
	yL_cfd = 2048.d0
	zL_cfd = 128.d0
	icmin_olap = 1
	icmax_olap = ncx
	jcmin_olap = 1
	jcmax_olap = 23 
	kcmin_olap = 1
	kcmax_olap = ncz

	ncx_olap = icmax_olap - icmin_olap + 1
	ncy_olap = jcmax_olap - jcmin_olap + 1
	ncz_olap = kcmax_olap - kcmin_olap + 1

	! Calcs	
	dy = yL_cfd / ncy
	yL_md = yL_cfd / 2.d0
	yL_olap = (jcmax_olap - jcmin_olap + 1) * dy
	yLl_md = yL_md/npy_md
	yLl_cfd = yL_cfd/npy_cfd

	call setup_CFD_procs

end subroutine setup_input_and_arrays

subroutine setup_CFD_procs
	use coupler_module
	implicit none

	integer	:: n
	integer :: ncxl, ncyl, nczl

	ncxl = ncx / npx_cfd
	do n=1,npx_cfd
		icPmax_cfd(n) = n * ncxl
		icPmin_cfd(n) = icPmax_cfd(n) - ncxl + 1
	end do	

	ncyl = ncy / npy_cfd
	do n=1,npy_cfd
		jcPmax_cfd(n) = n * ncyl
		jcPmin_cfd(n) = jcPmax_cfd(n) - ncyl + 1
	end do

	nczl = ncz / npz_cfd
	do n=1,npz_cfd
		kcPmax_cfd(n) = n * nczl
		kcPmin_cfd(n) = kcPmax_cfd(n) - nczl + 1
	end do

end subroutine setup_CFD_procs


! ----------------------------------------------
! Test the packing routines from coupler

subroutine test_packing
	use coupler_module
	use coupler
	implicit none


	integer 										:: ncxl, ncyl, nczl
	integer 										:: coord(3), extents(6)
	integer											:: ixyz, icell, jcell, kcell
	double precision,dimension(:),allocatable		:: outbuf
	double precision,dimension(:,:,:,:),allocatable	:: packbuf,testbuf

	! Test Packing in the overlap region only
	if (olap_mask(rank_world).ne.1) return
				   
	if (realm .eq. md_realm) then		   

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_proc_extents(coord,md_realm,extents)

		ncxl = extents(2)-extents(1)+1
		ncyl = extents(4)-extents(3)+1
		nczl = extents(6)-extents(5)+1

		allocate(packbuf(3,ncxl,ncyl,nczl),testbuf(3,ncxl,ncyl,nczl))

		! Populate dummy packbuf
		do ixyz = 1,3
		do icell=1,extents(2)-extents(1)
		do jcell=1,extents(4)-extents(3)
		do kcell=1,extents(6)-extents(5)

			packbuf(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                       			  1000*jcell + &
			                    			  1000000*kcell

		end do
		end do
		end do
		end do


	else if (realm .eq. cfd_realm) then	

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,extents)

		ncxl = extents(2)-extents(1)+1
		ncyl = extents(4)-extents(3)+1
		nczl = extents(6)-extents(5)+1

		allocate(packbuf(3,ncxl,ncyl,nczl),testbuf(3,ncxl,ncyl,nczl))

		! Populate dummy packbuf
		do ixyz = 1,3
		do icell=1,extents(2)-extents(1)
		do jcell=1,extents(4)-extents(3)
		do kcell=1,extents(6)-extents(5)

			packbuf(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                       			  1000*jcell + &
			                    			  1000000*kcell

		end do
		end do
		end do
		end do

	end if

	!print*,'Test pack',rank_world,myid_cart,realm,olap_mask(rank_world), extents, coord

	call CPL_pack(packbuf,outbuf,realm)
	call CPL_unpack(outbuf,testbuf,realm)

	print'(a,3f10.5)', 'Error in pack/unpack = ', maxval(testbuf-packbuf), & 
												  minval(testbuf-packbuf), & 
												     sum(testbuf-packbuf)

end subroutine test_packing

! ----------------------------------------------
! Test the send and recv routines from coupler

subroutine test_send_recv_MD2CFD
	use coupler_module
	use coupler
	implicit none

	integer :: ncxl, ncyl, nczl
	double precision,dimension(:,:,:,:),allocatable	:: sendbuf,recvbuf

	! Test Sending from CFD to MD							   
	if (realm .eq. md_realm) then		   
		ncxl = icPmax_md(iblock_realm) - icPmin_md(iblock_realm) + 1
		ncyl = 1 !jcPmax_md(jblock_realm) - jcPmin_md(jblock_realm) 
		nczl = kcPmax_md(kblock_realm) - kcPmin_md(kblock_realm) + 1
		!print*, 'sent size',realm_name(realm),3*ncxl*ncyl*nczl 
		allocate(sendbuf(3,ncxl,ncyl,nczl))
		sendbuf = 0.d0
		sendbuf = rank_realm*1000 + iblock_realm*100 + jblock_realm*10 + kblock_realm*1
		call CPL_send(sendbuf,1,1)		   
	else if (realm .eq. cfd_realm) then	   
		ncxl = ncx/npx_cfd 	!icPmax_cfd(iblock_realm) - icPmin_cfd(iblock_realm)
		ncyl = 1	   		!jcPmax_cfd(jblock_realm) - jcPmin_cfd(jblock_realm) 
		nczl = ncz/npz_cfd 	!kcPmax_cfd(kblock_realm) - kcPmin_cfd(kblock_realm) 
		!print*, 'recv size', realm_name(realm),3*ncxl*ncyl*nczl 
		allocate(recvbuf(3,ncxl,ncyl,nczl))
		recvbuf = 0.d0
		call CPL_recv(recvbuf,1,1)
	end if								   

	 if (realm .eq.  md_realm) write(4000+myid_world,*),myid_world, 'BUF=', sendbuf
	 if (realm .eq. cfd_realm) write(5000+myid_world,*),myid_world, 'BUF=', recvbuf
	
end subroutine test_send_recv_MD2CFD

! Test Sending from MD to CFD

subroutine test_send_recv_CFD2MD
	use coupler_module
	use coupler
	implicit none

	integer :: ncxl, ncyl, nczl
	double precision,dimension(:,:,:,:),allocatable	:: sendbuf,recvbuf

	! Test Sending from CFD to MD							   
	if (realm .eq. md_realm) then		   
		ncxl = icPmax_md(iblock_realm) - icPmin_md(iblock_realm) + 1
		ncyl = 1 !jcPmax_md(jblock_realm) - jcPmin_md(jblock_realm) 
		nczl = kcPmax_md(kblock_realm) - kcPmin_md(kblock_realm) + 1
		!print*, 'recv size', realm_name(realm),3*ncxl*ncyl*nczl 
		allocate(recvbuf(3,ncxl,ncyl,nczl)); recvbuf = 0.d0
		call CPL_recv(recvbuf,1,1)		   
	else if (realm .eq. cfd_realm) then	   
		ncxl = ncx/npx_cfd 	!icPmax_cfd(iblock_realm) - icPmin_cfd(iblock_realm)
		ncyl = 1	   		!jcPmax_cfd(jblock_realm) - jcPmin_cfd(jblock_realm) 
		nczl = ncz/npz_cfd 	!kcPmax_cfd(kblock_realm) - kcPmin_cfd(kblock_realm) 
		!print*, 'sent size',realm_name(realm),3*ncxl*ncyl*nczl 
		allocate(sendbuf(3,ncxl,ncyl,nczl))
		sendbuf = rank_realm*1000 + iblock_realm*100 + jblock_realm*10 + kblock_realm*1
		call CPL_send(sendbuf,1,1)
	end if								   

	 if (realm .eq.  md_realm) write(9000+myid_world,*),myid_world, 'BUF=', recvbuf
	 if (realm .eq. cfd_realm) write(11000+myid_world,*),myid_world, 'BUF=', sendbuf
	
end subroutine test_send_recv_CFD2MD



subroutine test_gather_scatter
	use coupler_module
	use coupler
	implicit none

	double precision,dimension(:,:,:,:),allocatable	:: u,stress
	integer :: coord(3), extents(6), npercell
	integer :: pos, ixyz, icell, jcell, kcell

	if (realm .eq. md_realm) then	

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,md_realm,3,coord,ierr)
		call CPL_proc_extents(coord,md_realm,extents)
		npercell = 3
		allocate(u(npercell,extents(1):extents(2), &
		                    extents(3):extents(4), &
		                    extents(5):extents(6)))
		allocate(stress(0,0,0,0))

!		print*, rank_world, 'main extents', extents

		! Populate dummy gatherbuf
		pos = 1
		do ixyz = 1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)

			u(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                                      1000*jcell + &
			                                   1000000*kcell
			pos = pos + 1

		end do
		end do
		end do
		end do

	else if (realm .eq. cfd_realm) then	  
		
		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,extents)
		npercell = 9
		allocate(u(0,0,0,0))
		allocate(stress(npercell,extents(1):extents(2), &
		                         extents(3):extents(4), &
		                         extents(5):extents(6)))

		! Populate dummy gatherbuf
		pos = 1
		do ixyz = 1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)

			stress(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                                           1000*jcell + &
			                                        1000000*kcell
			pos = pos + 1

		end do
		end do
		end do
		end do

	endif

	if (olap_mask(rank_world).eq.1) call CPL_gather(u,3)
	if (olap_mask(rank_world).eq.1) call CPL_scatter(stress,9)

	
end subroutine test_gather_scatter

				 
subroutine test_comms
	use coupler_module
	use mpi
	implicit none
	integer					:: i
	integer, dimension(3)	:: coords
	double precision		:: rand


	!if (realm .eq. md_realm) then
	!	call CPL_Cart_coords(CPL_WORLD_COMM, myid_world+1, md_realm, 3, coords, ierr)
	!	print'(a,5i8)', 'CPL_CART_COORDS WORLD_COMM MD ', realm, myid_world+1, coords
	!elseif (realm .eq. cfd_realm) then
	!	call CPL_Cart_coords(CPL_WORLD_COMM, myid_world+1, cfd_realm, 3, coords, ierr)
	!	print'(a,5i8)', 'CPL_CART_COORDS WORLD_COMM CFD', realm, myid_world+1, coords
	!endif

	!Test loop on a single random processor
	
	! Get a random proccessor
	if (myid_world .eq. 0) then
		CALL random_seed()
		call random_number(rand)
	endif
	call MPI_bcast(rand,1,MPI_DOUBLE_PRECISION,0,CPL_WORLD_COMM,ierr)

	!Only loop on random processor within world processors
	if (myid_world .eq. floor(rand*nproc_world)) then

		!World
		do i=1,nproc_md
			call CPL_Cart_coords(CPL_WORLD_COMM, i, md_realm, 3, coords, ierr)
			print'(a,6i8)', 'CPL_WORLD_MD ',myid_world, md_realm, i, coords
		enddo
		do i=nproc_md+1,nproc_world
			call CPL_Cart_coords(CPL_WORLD_COMM, i, cfd_realm, 3, coords, ierr)
			print'(a,6i8)', 'CPL_WORLD_CFD', myid_world,cfd_realm, i, coords
		enddo

		!Realm
		do i=1,nproc_cfd
			call CPL_Cart_coords(CPL_REALM_COMM, i, cfd_realm, 3, coords, ierr)
			print'(a,6i8)', 'CPL_REALM_CFD', myid_world,cfd_realm, i, coords
		enddo
		do i=1,nproc_md
			call CPL_Cart_coords(CPL_REALM_COMM, i, md_realm, 3, coords, ierr)
			print'(a,6i8)', 'CPL_REALM_MD ',myid_world, md_realm, i, coords
		enddo

		!Cart
		do i=1,nproc_cfd
			call CPL_Cart_coords(CPL_CART_COMM, i, cfd_realm, 3, coords, ierr)
			print'(a,6i8)', 'CPL_CART_CFD', myid_world,cfd_realm, i, coords
		enddo
		do i=1,nproc_md
			call CPL_Cart_coords(CPL_CART_COMM, i, md_realm, 3, coords, ierr)
			print'(a,6i8)', 'CPL_CART_MD ',myid_world, md_realm, i, coords
		enddo

	endif

	call barrier

	!Only within random overlap processors
	if (myid_olap .eq. floor(rand*nproc_olap)) then

		!Olap
		call CPL_Cart_coords(CPL_OLAP_COMM, 1, cfd_realm, 3, coords, ierr)
		print'(a,6i8)', 'CPL_OLAP', myid_world,cfd_realm, 1, coords
		do i=2,nproc_olap
			call CPL_Cart_coords(CPL_OLAP_COMM, i, md_realm, 3, coords, ierr)
			print'(a,6i8)', 'CPL_OLAP', myid_world,cfd_realm, i, coords
		enddo

		!Graph
		call CPL_Cart_coords(CPL_GRAPH_COMM, 1, cfd_realm, 3, coords, ierr)
		print'(a,6i8)', 'CPL_GRAPH', myid_world,cfd_realm, 1, coords
		do i=2,nproc_olap
			call CPL_Cart_coords(CPL_GRAPH_COMM, i, md_realm, 3, coords, ierr)
			print'(a,6i8)', 'CPL_GRAPH', myid_world,cfd_realm, i, coords
		enddo

	endif

end subroutine test_comms


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
				call error_abort("CPL_Cart_coords Error - Setup not complete for CFD CPL_REALM_COMM")
			elseif (rank .gt. size(rank_cfdrealm2rank_world)) then
				print*, 'rank = ', rank, 'comm size = ', size(rank_cfdrealm2rank_world)
				call error_abort("CPL_Cart_coords Error -Specified rank is not in CFD realm")
			endif
			worldrank = rank_cfdrealm2rank_world(rank)
		elseif (realm .eq. md_realm) then
			if (allocated(rank_mdrealm2rank_world) .eqv. .false.) then
				call error_abort("CPL_Cart_coords Error - Setup not complete for MD CPL_REALM_COMM")
			elseif (rank .gt. size(rank_mdrealm2rank_world)) then
				print*, 'rank = ', rank, 'comm size = ', size(rank_cfdrealm2rank_world)
				call error_abort("CPL_Cart_coords Error -Specified rank is not in MD realm")
			endif
			worldrank = rank_mdrealm2rank_world(rank)
		endif
		! -  -  -  -  -  -  -  -  -  -  -  -  -
	elseif(COMM .eq. CPL_CART_COMM) then
		! -  -  -  -  -  -  -  -  -  -  -  -  -
		if (realm .eq. cfd_realm) then
			if (allocated(rank2coord_cfd) .eqv. .false.) &
				call error_abort("CPL_Cart_coords Error - Setup not complete for CFD CPL_CART_COMM")
			coords = rank2coord_cfd(:,rank)
		elseif (realm .eq. md_realm) then
			if (allocated(rank2coord_md) .eqv. .false.) &
				call error_abort("CPL_Cart_coords Error - Setup not complete for MD CPL_CART_COMM")
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
		call error_abort("CPL_Cart_coords Error - Intercomm between realms id - use realm comms instead")
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
			call error_abort("CPL_Cart_coords Error - world to cart mapping not initialised correctly")
		endif
		cartrank = rank_world2rank_cfdcart(worldrank)
		if (cartrank .eq. VOID) call error_abort("CPL_Cart_coords Error - void element in mapping")
		if (cartrank .gt. nproc_cfd) call error_abort("CPL_Cart_coords Error - rank not in cfd realm")
	elseif (realm .eq. md_realm) then
		if (allocated(rank_world2rank_mdcart) .eqv. .false.) then
			call error_abort("CPL_Cart_coords Error - world to cart mapping not initialised correctly")
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


!subroutine get_overlap_gridstretch
!implicit none
!	do j=1,ngy
!		if (in_box(x(j),y(j),z(j))) then
!			olap_jmin = j
!			do j=olap_minj+1,ngy
!				if (.not. in_box(x(j),y(j),z(j))) then
!					olap_jmax = j - 1
!					exit
!				end if
!			end do
!			exit
!		end if
!		olap_jmin = olap_null
!		olap_jmax = olap_null
!	end do
!
!end subroutine get_overlap_gridstretch

!=================================================================
!module md_proc_box
!use cpl_common
!
!	!6 planes define MD processor's domain "position" and size
!	double precision :: &
!		max_plane(3),   &
!		min_plane(3)
!	
!	max_plane(1) = picoord*xL_md
!	min_plane(1) = (picoord-1)*xL_md	
!
!	max_plane(2) = pjcoord*yL_md
!	min_plane(2) = (pjcoord-1)*yL_md	
!
!	max_plane(3) = pkcoord*zL_md
!	min_plane(3) = (pkcoord-1)*zL_md	
!
!contains
!
!	logical function in_box(x,y,z)
!		real(kind(0.d0)) :: x,y,z
!		logical in_box
!			
!		if (x.le.max_plane(1) .and. &
!		    x.gt.min_plane(1) .and. &
!		    y.le.max_plane(2) .and. &
!		    y.gt.min_plane(2) .and. &
!		    z.le.max_plane(3) .and. &
!		    z.gt.min_plane(3)) then
!			in_box = .true.
!		else
!			in_box = .false.
!		end if
!		return
!	end
!	
!end module md_proc_box
!
