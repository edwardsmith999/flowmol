program test_coupler
	use coupler_module, only : olap_mask, rank_world
	use coupler
	implicit none


	!SETUP PROCESS
	call initialise                        ! FROM TEST
	call test_setup_input_and_arrays       ! FROM TEST
	call get_cfd_cell_ranges               ! FROM TEST
	call create_realms                     ! FROM TEST
	call CPL_create_map	                   ! FROM COUPLER
	!call get_md_cell_ranges                ! FROM COUPLER        
	!call get_overlap_blocks                ! FROM COUPLER
	!call prepare_overlap_comms             ! FROM COUPLER
	!call CPL_overlap_topology              ! FROM COUPLER

	call barrier

	!EXCHANGE/TESTING ROUTINES
	!call test_COMMS
	!call test_packing
	call test_gather_scatter           ! FROM TEST
	call test_send_recv_MD2CFD         ! FROM TEST
	call test_send_recv_CFD2MD         ! FROM TEST


	call finalise                      ! FROM TEST

end program test_coupler

!===========================================
!
!		SETUPS USED FOR DUMMY COUPLER
!
!===========================================

subroutine test_setup_input_and_arrays
	use coupler_module
	implicit none

	open(unit=1,file='input',form='formatted')
		read(1,*) npx_md 
		read(1,*) npy_md 
		read(1,*) npz_md 
		read(1,*) npx_cfd 
		read(1,*) npy_cfd 
		read(1,*) npz_cfd 
		read(1,*) ncx 
		read(1,*) ncy 
		read(1,*) ncz 
		read(1,*) xL_cfd 
		read(1,*) yL_cfd 
		read(1,*) zL_cfd 
		read(1,*) icmin_olap 
		read(1,*) icmax_olap 
		read(1,*) jcmin_olap 
		read(1,*) jcmax_olap 
		read(1,*) kcmin_olap 
		read(1,*) kcmax_olap 
	close(1)


	! Calcs	
	ncx_olap = icmax_olap - icmin_olap + 1
	ncy_olap = jcmax_olap - jcmin_olap + 1
	ncz_olap = kcmax_olap - kcmin_olap + 1

	nproc_md  = npx_md*npy_md*npz_md
	nproc_cfd = npx_cfd*npy_cfd*npz_cfd
	nproc_world = nproc_md + nproc_cfd

	dx = xL_cfd / ncx
	dy = yL_cfd / ncy
	dz = zL_cfd / ncz

	! TODO TODO TODO TODO TODO TODO
	xL_md = xL_cfd
	yL_md = yL_cfd / 2.d0
	zL_md = zL_cfd

	yL_olap = (jcmax_olap - jcmin_olap + 1) * dy

end subroutine test_setup_input_and_arrays

subroutine get_cfd_cell_ranges
	use coupler_module
	implicit none

	integer	:: n
	integer :: ncxl, ncyl, nczl

	allocate(icPmin_cfd(npx_cfd))
	allocate(icPmax_cfd(npx_cfd))
	allocate(jcPmin_cfd(npy_cfd))
	allocate(jcPmax_cfd(npy_cfd))
	allocate(kcPmin_cfd(npz_cfd))
	allocate(kcPmax_cfd(npz_cfd))

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

end subroutine get_cfd_cell_ranges 

!=========================================================================
subroutine create_realms
	use coupler_module
	use mpi
	implicit none

	integer :: gridsize(3),coord(3)
	integer	::  callingrealm,ibuf(2),jbuf(2),remote_leader,comm,comm_size
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
	call MPI_comm_rank(CPL_REALM_COMM,myid_realm,ierr)
	rank_realm = myid_realm + 1; rootid_realm = 0

	! Get the MPI_comm_world ranks that hold the largest ranks in cfd_comm and md_comm
	call MPI_comm_size(CPL_REALM_COMM,comm_size,ierr)
	ibuf(:) = -1
	jbuf(:) = -1
	if ( myid_realm .eq. comm_size - 1) then
		ibuf(realm) = myid_world
	endif

	call MPI_allreduce( ibuf ,jbuf, 2, MPI_INTEGER, MPI_MAX, &
						CPL_WORLD_COMM, ierr)

	!Set this largest rank on each process to be the inter-communicators (WHY NOT 0??)
	select case (realm)
	case (cfd_realm)
		remote_leader = jbuf(md_realm)
	case (md_realm)
		remote_leader = jbuf(cfd_realm)
	end select

	call MPI_intercomm_create(CPL_REALM_COMM, comm_size - 1, CPL_WORLD_COMM,&
									remote_leader, 1, CPL_INTER_COMM, ierr)

!	write(0,*) 'did (inter)communicators ', realm_name(realm), myid_world

	!Setup cartesian topology
	call MPI_cart_create(CPL_REALM_COMM,3,gridsize,periodicity,.true., &
						 CPL_CART_COMM,ierr)
	call MPI_comm_rank(CPL_CART_COMM,myid_cart,ierr)
	rank_cart = myid_cart + 1

	call collect_rank_mappings

end subroutine create_realms

subroutine collect_rank_mappings
	use coupler_module
	implicit none

	integer :: coord(3)

	allocate(coord2rank_md(npx_md,npy_md,npz_md))
	allocate(coord2rank_cfd(npx_cfd,npy_cfd,npz_cfd))
	allocate(rank2coord_cfd(3,nproc_cfd))
	allocate(rank2coord_md(3,nproc_md))
	allocate(rank_mdcart2rank_world(nproc_md))
	allocate(rank_cfdcart2rank_world(nproc_cfd))

	coord2rank_md           = 0
	coord2rank_cfd          = 0
	rank2coord_md           = 0
	rank2coord_cfd          = 0
	rank_mdcart2rank_world  = 0
	rank_cfdcart2rank_world = 0

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

end subroutine collect_rank_mappings

subroutine collect_coord2ranks
	use mpi
	use coupler_module
	implicit none

	integer :: coord(3)
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
	use coupler, only : CPL_rank_map
	implicit none
	
	integer							   :: buf, source, nproc
	integer, dimension(:), allocatable :: mbuf, cbuf
	integer, dimension(:), allocatable :: rank_cart2rank_world,rank_world2rank_cart
	integer, dimension(:), allocatable :: rank_realm2rank_world,rank_world2rank_realm

	!------------------------ Cart------------------
	call CPL_rank_map(CPL_CART_COMM,rank_cart,nproc, & 
					 rank_cart2rank_world,rank_world2rank_cart,ierr)

	!World to rank
	allocate(rank_world2rank_cfdcart(nproc_world))
	allocate(rank_world2rank_mdcart(nproc_world))
	rank_world2rank_cfdcart = rank_world2rank_cart
	rank_world2rank_mdcart  = rank_world2rank_cart

	!print*, 'world to cart', nproc, nproc_cfd, rank_world2rank_mdcart

	! Rank to world
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

	! - - Collect on own realm intracomm - -
	call CPL_rank_map(CPL_REALM_COMM,rank_realm,nproc, & 
					 rank_realm2rank_world,rank_world2rank_realm,ierr)

	!World to rank
	allocate(rank_world2rank_cfdrealm(nproc_world))
	allocate(rank_world2rank_mdrealm(nproc_world))
	rank_world2rank_cfdrealm = rank_world2rank_realm
	rank_world2rank_mdrealm  = rank_world2rank_realm

	!print*, 'world to realm', nproc, nproc_cfd, rank_world2rank_mdrealm

	!Rank to world
	if (realm .eq. cfd_realm) then	
		allocate(rank_cfdrealm2rank_world(nproc))
		rank_cfdrealm2rank_world = rank_realm2rank_world
		!print*, 'CFD realm to world', nproc, nproc_cfd, rank_cfdrealm2rank_world
	elseif (realm .eq. md_realm) then
		allocate(rank_mdrealm2rank_world(nproc))
		rank_mdrealm2rank_world = rank_realm2rank_world
		!print*, 'MD realm to world', nproc, nproc_md, rank_mdrealm2rank_world
	endif

	!  - - Exchange across intercomm sending only from root processor  - - 
    if (myid_realm .eq. rootid_realm ) then
        source = MPI_ROOT
    else
        source = MPI_PROC_NULL
    endif

	if (realm .eq. cfd_realm) then
		call MPI_bcast(rank_cfdrealm2rank_world,nproc_cfd, & 
								MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
		allocate(rank_mdrealm2rank_world(nproc_md))
		call MPI_bcast(rank_mdrealm2rank_world,nproc_md, & 
								MPI_INTEGER,0,CPL_INTER_COMM,ierr) 		!Receive
		!print*, 'CFD realm to world MD', rank_mdrealm2rank_world
	elseif (realm .eq. md_realm) then
		allocate(rank_cfdrealm2rank_world(nproc_cfd))
		call MPI_bcast(rank_cfdrealm2rank_world,nproc_cfd, & 
								MPI_INTEGER,0,CPL_INTER_COMM,ierr) 		!Receive
		call MPI_bcast(rank_mdrealm2rank_world,nproc_md, & 
								MPI_INTEGER,source,CPL_INTER_COMM,ierr)	!Send
		!print*, 'MD realm to world CFD', rank_cfdrealm2rank_world
	endif

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

	write(2000+rank_realm,'(2i7,a5,3i5,a5,2i10,a5,i20)'), &
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




!===========================================
!
!		TESTING USED IN DUMMY COUPLER
!
!===========================================




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

! ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬ஜ۩۞۩ஜ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
! Test the send and recv routines from coupler

subroutine test_send_recv_MD2CFD
	use coupler_module
	use coupler
	implicit none

	logical	:: send_flag,recv_flag
	integer :: ncxl,ncyl,nczl,ixyz,icell,jcell,kcell
	integer	:: jcmin_send,jcmax_send,jcmin_recv,jcmax_recv,npercell,coord(3),extents(6)
	double precision,dimension(:,:,:,:),allocatable	:: sendbuf,recvbuf

	npercell = 3
	jcmax_send=1; jcmin_send=1; 
	jcmax_recv = jcmax_send
	jcmin_recv = jcmin_send
	if (olap_mask(rank_world) .eq. 0) return

	! Test Sending from MD to CFD							   
	if (realm .eq. md_realm) then	

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_olap_extents(coord,realm,extents)

		allocate(sendbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))

		!print'(2a,11i7)', 'sent size',realm_name(realm),extents,size(sendbuf),shape(sendbuf)

		! Populate dummy gatherbuf
		sendbuf = -333.d0 ! 0.d0
		do ixyz = 1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)
			sendbuf(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                                            1000*jcell + &
			                                         1000000*kcell
		end do
		end do
		end do
		end do

		call CPL_send(sendbuf,jcmax_send=jcmax_send,jcmin_send=jcmin_send,send_flag=send_flag)	

		if (send_flag .eqv. .true.) then
			do kcell=extents(5),extents(6)
			do jcell=jcmin_send,jcmax_send
			do icell=extents(1),extents(2)
			do ixyz =1,npercell
				write(4000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
				      'send MD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
				       sendbuf(ixyz,icell,jcell,kcell)
			end do
			end do
			end do
			end do
		endif

	else if (realm .eq. cfd_realm) then	 

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_proc_extents(coord,realm,extents)
		!print'(2a,8i7)', 'proc extents', realm_name(realm),rank_world,rank_cart,extents
		call CPL_olap_extents(coord,realm,extents)
		!print'(2a,8i7)', 'olap extents', realm_name(realm),rank_world,rank_cart,extents

		allocate(recvbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))

		!print'(2a,11i7)', 'recv size', realm_name(realm),extents,size(recvbuf),shape(recvbuf)
		recvbuf = -444.d0
		call CPL_recv(recvbuf,jcmax_recv=jcmax_recv,jcmin_recv=jcmin_recv,recv_flag=recv_flag)

		if (recv_flag .eqv. .true.) then
			do kcell=extents(5),extents(6)
			do jcell=jcmin_recv,jcmax_recv  !extents(3),extents(4)
			do icell=extents(1),extents(2)
			do ixyz =1,npercell
				!if ( recvbuf(ixyz,icell,jcell,kcell) .ne. -444.d0) then
				!print'(a,i4,a,i4,a,i4,a,i4,a,f20.1)',   &
				!      'recv CFD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
				!       recvbuf(ixyz,icell,jcell,kcell)
					write(5000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
					      'recv CFD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
					       recvbuf(ixyz,icell,jcell,kcell)
				!endif

			end do
			end do
			end do
			end do
		endif
	end if								   

	! if (realm .eq.  md_realm) write(4000+myid_world,*),myid_world, 'BUF=', sendbuf
	! if (realm .eq. cfd_realm) write(5000+myid_world,*),myid_world, 'BUF=', recvbuf
	
end subroutine test_send_recv_MD2CFD


! ۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩
! Test Sending from MD to CFD

subroutine test_send_recv_CFD2MD
	use coupler_module
	use coupler
	implicit none

	logical	:: send_flag,recv_flag
	integer	:: jcmin_send,jcmax_send,jcmin_recv,jcmax_recv
	integer :: ncxl,ncyl,nczl,ixyz,icell,jcell,kcell,npercell,coord(3),extents(6)
	double precision,dimension(:,:,:,:),allocatable	:: sendbuf,recvbuf

	npercell = 3
	jcmax_send=1; jcmin_send=1; 
	jcmax_recv = jcmax_send
	jcmin_recv = jcmin_send
	if (olap_mask(rank_world) .eq. 0) return

	! Test Sending from CFD to MD							   
	if (realm .eq. md_realm) then		   

		coord = (/iblock_realm,jblock_realm,kblock_realm /)
		call CPL_olap_extents(coord,realm,extents)

		allocate(recvbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))
		recvbuf = -444

		!print*, 'recv size', realm_name(realm),extents, size(recvbuf),shape(recvbuf)
		call CPL_recv(recvbuf,jcmax_recv=1,jcmin_recv=1,recv_flag=recv_flag)   

		if (recv_flag .eqv. .true.) then
			do kcell=extents(5),extents(6)
			do jcell=jcmin_send,jcmax_send
			do icell=extents(1),extents(2)
			do ixyz = 1,npercell
				write(2000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
			      	'recv MD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
			      	 recvbuf(ixyz,icell,jcell,kcell)
			end do
			end do
			end do
			end do
		endif

	else if (realm .eq. cfd_realm) then	   

		coord = (/iblock_realm,jblock_realm,kblock_realm /)
		call CPL_olap_extents(coord,realm,extents)
		allocate(sendbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))

		do ixyz =1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)
			sendbuf(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*(icell) + &
			                       			  1000*(jcell) + &
			                    			  1000000*(kcell)

		end do
		end do
		end do
		end do

		!print*, 'sent size',realm_name(realm),3*ncxl*ncyl*nczl,size(sendbuf)
		call CPL_send(sendbuf,jcmax_send=1,jcmin_send=1,send_flag=send_flag)

		do kcell=extents(5),extents(6)
		do jcell=jcmin_send,jcmax_send
		do icell=extents(1),extents(2)
		do ixyz = 1,npercell
			write(9000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
			      'send CFD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
			       sendbuf(ixyz,icell,jcell,kcell)
		end do
		end do
		end do
		end do
	end if								   
	
end subroutine test_send_recv_CFD2MD

subroutine test_gather_scatter
	use coupler_module
	use coupler
	implicit none

	double precision,dimension(:,:,:,:),allocatable	:: u,stress,gatheru,scatterstress
	integer :: coord(3), extents(6), gatherlims(6), scatterlims(6), npercell
	integer :: pos, ixyz, icell, jcell, kcell
	integer :: ncxl,ncyl,nczl
	integer :: i,j,k

 	if (olap_mask(rank_world).ne.1) return

	if (realm .eq. md_realm) then	

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,md_realm,3,coord,ierr)
		call CPL_proc_extents(coord,md_realm,extents)
		npercell = 3
		allocate(u(npercell,extents(1):extents(2), &
		                    extents(3):extents(4), &
		                    extents(5):extents(6)))
		allocate(stress(0,0,0,0))

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


	! Allocate test arrays over local domain
	if (realm.eq.cfd_realm) then
		call CPL_cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,extents)
		ncxl = extents(2) - extents(1) + 1
		ncyl = extents(4) - extents(3) + 1
		nczl = extents(6) - extents(5) + 1
		allocate(gatheru(3,ncxl,ncyl,nczl))
		gatheru = 0.d0
	else if (realm.eq.md_realm) then
		call CPL_cart_coords(CPL_CART_COMM,rank_cart,md_realm,3,coord,ierr)
		call CPL_proc_extents(coord,md_realm,extents)
		ncxl = extents(2) - extents(1) + 1
		ncyl = extents(4) - extents(3) + 1
		nczl = extents(6) - extents(5) + 1
		allocate(scatterstress(9,ncxl,ncyl,nczl))
		scatterstress = 0.d0
	end if




	!gatherlims  = (/1,1,1,1,1,1/)
	!scatterlims = (/1,1,1,1,1,1/)
	!================== PERFORM GATHER/SCATTER =============================!	
	gatherlims  = (/1,ncx, jcmax_olap, jcmax_olap , 1, ncz/)
	scatterlims = (/1,ncx, 1, 1, 1,ncz/)
	if (olap_mask(rank_world).eq.1) call CPL_gather(u,3,gatherlims,gatheru)
	if (olap_mask(rank_world).eq.1) call CPL_scatter(stress,9,scatterlims, &
	                                                 scatterstress)

	! Print results to file
	if (realm.eq.cfd_realm) then

		do ixyz  = 1,size(gatheru,1)
		do icell = 1,size(gatheru,2)
		do jcell = 1,size(gatheru,3)
		do kcell = 1,size(gatheru,4)

			i = icell + extents(1) - 1
			j = jcell + extents(3) - 1
			k = kcell + extents(5) - 1

			if (gatheru(ixyz,icell,jcell,kcell).lt.0.0001) then
				!write(8000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
				!	  'gatheru(',0,',',0,',',0,',',0,') =', 0.d0
			else
				write(8000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
					  'gatheru(',ixyz,',',i,',',j,',',k,') =', &
					   gatheru(ixyz,icell,jcell,kcell)
			end if

		end do	
		end do	
		end do
		end do

	else if (realm.eq.md_realm) then

		do ixyz  = 1,size(scatterstress,1)
		do icell = 1,size(scatterstress,2)
		do jcell = 1,size(scatterstress,3)
		do kcell = 1,size(scatterstress,4)

			i = icell + extents(1) - 1
			j = jcell + extents(3) - 1
			k = kcell + extents(5) - 1

			if (scatterstress(ixyz,icell,jcell,kcell).lt.0.0001) then
				!write(7000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
				!	  'scatterstress(',0,',',0,',',0,',',0,') =', 0.d0
			else
				write(7000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
					  'scatterstress(',ixyz,',',i,',',j,',',k,') =', &
					   scatterstress(ixyz,icell,jcell,kcell)
			end if

		end do	
		end do	
		end do
		end do
	
	end if
	
end subroutine test_gather_scatter

				 
subroutine test_comms
	use coupler_module
	use mpi
	use coupler, only : CPL_Cart_coords
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

