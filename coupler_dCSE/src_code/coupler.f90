!=============================================================================
!
!				  Coupler 
!
! Routines accessible from application ( molecular or continuum ) after 
! the name, in parenthesis, is the realm in which each routine must be called
!
! SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP
!
! coupler_create_comm	      (cfd+md)   splits MPI_COMM_WORLD, create inter - 
!										 communicator between CFD and MD
! coupler_create_map	      (cfd+md)   creates correspondence maps between 
!										 the CFD grid and MD domains
! coupler_cfd_init              (cfd)    initialises coupler with CFD data
! coupler_md_init               (cfd)    initialises coupler and set MD 
!										 parameters with using data from CFD 
!										 or COUPLER.in
! coupler_cfd_adjust_domain     (cfd)    adjust CFD tomain to an integer number 
!										 FCC or similar MD initial layout
!
! SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION
!
! coupler_send_data        	  (cfd+md)   sends grid data exchanged between 
!										 realms ( generic interface)
! coupler_recv_data        	  (cfd+md)   receives data exchanged between realms 
!										 ( generic interface)
! coupler_cfd_get               (cfd)    returns coupler internal parameters 
!										 for CFD realm
! coupler_md_get                 (md)    returns coupler internal parameters 
!										 for MD realm
! coupler_md_get_save_period     (md)    auxiliary used for testing
! coupler_md_get_average_period  (md)    returns average period of BC
! coupler_md_get_md_per_cfd_dt   (md) 	 returns the number of step MD does for 
										 !each CFD step
! coupler_md_get_nsteps          (md)    returm CFD nsteps  
! coupler_md_get_dt_cfd          (md)    returns MD dt
! coupler_md_set                 (md)    sets zL if CFD is 2D
! coupler_md_get_density         (md)    gets CFD density
! coupler_md_get_cfd_id          (md)    id for CFD code, possible values set 
!										 in coupler_parameters
!
!  Lucian Anton, November 2011
!  Revised April 2012
!
!=============================================================================

module coupler
	use coupler_parameters
    implicit none
    save

    interface CPL_send
        module procedure CPL_send_3d, CPL_send_4d
    end interface

    interface CPL_recv
        module procedure CPL_recv_3d, CPL_recv_4d
    end interface

    private CPL_send_3d, CPL_send_4d, &
        CPL_send_xd, CPL_recv_3d, CPL_recv_4d,&
        CPL_recv_xd

    interface coupler_send_data
        module procedure coupler_send_data_3d, coupler_send_data_4d
    end interface

    interface coupler_recv_data
        module procedure coupler_recv_data_3d, coupler_recv_data_4d
    end interface

    private coupler_send_data_3d, coupler_send_data_4d, &
        coupler_send_data_xd, coupler_recv_data_3d, coupler_recv_data_4d,&
        coupler_recv_data_xd

contains

!=============================================================================
! Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup
!
!								SETUP
!
! Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup
!=============================================================================

!=============================================================================
! 							coupler_create_comm	      	
! (cfd+md) Splits MPI_COMM_WORLD in both the CFD and MD code respectively
! 		   and create intercommunicator between CFD and MD
!-----------------------------------------------------------------------------

subroutine coupler_create_comm(callingrealm, RETURNED_REALM_COMM, ierror)
	use mpi
	use coupler_module, only : rank_world,myid_world,rootid_world,nproc_world, & 
							   realm, rank_realm,myid_realm,rootid_realm, ierr, & 
								CPL_WORLD_COMM, CPL_REALM_COMM, CPL_INTER_COMM
	use coupler_input_data
	implicit none

	integer, intent(in) :: callingrealm ! CFD or MD
	integer, intent(out):: RETURNED_REALM_COMM, ierror

	!Get processor id in world across both realms
	call MPI_comm_rank(MPI_COMM_WORLD,myid_world,ierr)
	rank_world = myid_world + 1; rootid_world = 0
	call MPI_comm_size(MPI_COMM_WORLD,nproc_world,ierr)

	! test if we have a CFD and a MD realm
	ierror=0
	! Test realms are assigned correctly
	call test_realms	

	! Create intercommunicator between realms		
	realm = callingrealm
	call create_comm			

contains

!-----------------------------------------------------------------------------
!	Test if CFD and MD realms are assigned correctly
!-----------------------------------------------------------------------------

subroutine test_realms
	implicit none

	integer 			 :: i, root, nproc, ncfd, nmd
	integer, allocatable :: realm_list(:)

	! Allocate and gather array with realm (MD or CFD) of each 
	! processor in the coupled domain on the root processor
	root = 1
	if (rank_world .eq. root) then
		call MPI_comm_size(MPI_comm_world, nproc, ierr)
		allocate(realm_list(nproc))
	endif
	call MPI_gather(callingrealm,1,MPI_INTEGER,realm_list,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	!Check through array of processors on both realms
	!and return error if wrong values or either is missing
	if (rank_world .eq. root) then
		ncfd = 0; nmd = 0
		do i =1, nproc
			if ( realm_list(i) .eq. cfd_realm ) then 
				ncfd = ncfd + 1
			else if ( realm_list(i) .eq. md_realm ) then
				nmd = nmd +1
			else
				ierror = COUPLER_ERROR_REALM
				write(*,*) "wrong realm value in coupler_create_comm"
				call MPI_abort(MPI_COMM_WORLD,ierror,ierr)
			endif
		enddo

		if ( ncfd .eq. 0 .or. nmd .eq. 0) then 
			ierror = COUPLER_ERROR_ONE_REALM
			write(*,*) "CFD or MD realm is missing in MPI_COMM_WORLD"
			call MPI_abort(MPI_COMM_WORLD,ierror,ierr)
		endif

	endif

end subroutine test_realms

!-----------------------------------------------------------------------------
! Create communicators for each realm and inter-communicator
!-----------------------------------------------------------------------------

subroutine create_comm
	implicit none

	integer	::  callingrealm,ibuf(2),jbuf(2),remote_leader,comm,comm_size

	callingrealm = realm
	! Split MPI COMM WORLD ready to establish two communicators
	! 1) A global intra-communicator in each realm for communication
	! internally between CFD processes or between MD processes
	! 2) An inter-communicator which allows communication between  
	! the 'groups' of processors in MD and the group in the CFD 
	call MPI_comm_dup(MPI_COMM_WORLD,CPL_WORLD_COMM,ierr)
	RETURNED_REALM_COMM	= MPI_COMM_NULL
	CPL_REALM_COMM 		= MPI_COMM_NULL

	!------------ create realm intra-communicators -----------------------
	! Split MPI_COMM_WORLD into an intra-communicator for each realm 
	! (used for any communication within each realm - e.g. broadcast from 
	!  an md process to all other md processes) 
	call MPI_comm_split(CPL_WORLD_COMM,callingrealm,myid_world,RETURNED_REALM_COMM,ierr)

	!------------ create realm inter-communicators -----------------------
	! Create intercommunicator between the group of processor on each realm
	! (used for any communication between realms - e.g. md group rank 2 sends
	! to cfd group rank 5). inter-communication is by a single processor on each group
	! Split duplicate of MPI_COMM_WORLD
	call MPI_comm_split(CPL_WORLD_COMM,callingrealm,myid_world,CPL_REALM_COMM,ierr)
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

	!print*,color, jbuf, remote_leader

	call MPI_intercomm_create(CPL_REALM_COMM, comm_size - 1, CPL_WORLD_COMM,&
									remote_leader, 1, CPL_INTER_COMM, ierr)

	write(0,*) 'did (inter)communicators ', realm_name(realm), myid_world

end subroutine create_comm

end subroutine coupler_create_comm

!=============================================================================
! Establish for all MD processors the mapping (if any) 
! to coupled CFD processors
!-----------------------------------------------------------------------------

subroutine coupler_create_map
	use mpi
	use coupler_module
	implicit none

	integer		:: n

	!Get ranges of cells on each MD processor
	call get_md_cell_ranges

    ! Find overlaps between processors
    !call find_overlaps

	!Get overlapping mapping for MD to CFD
	call get_overlap_blocks

	!Setup overlap communicators
	call prepare_overlap_comms

	!Setup graph topology
	call CPL_overlap_topology


	!Write Debug information
	!if (realm .eq. cfd_realm) then
	!	call write_map_cfd
	!else if (realm .eq. md_realm) then
	!	call write_map_md
	!else
	!	write(*,*) "Wrong realm in coupler_create_map"
	!	call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_REALM,ierr)
	!end if
	!call write_map_olap

	call MPI_barrier(CPL_WORLD_COMM,ierr)

end subroutine coupler_create_map

!------------------------------------------------------------
!Calculate processor cell ranges of MD code on all processors
	
subroutine get_md_cell_ranges
	use coupler_input_data, only : cfd_coupler_input
	use coupler_module
	implicit none

	integer :: n
	integer :: olap_jmin_mdcoord
	integer :: ncxl, ncyl, nczl
	integer :: ncy_mdonly, ncy_md, ncyP_md

	allocate(icPmin_md(npx_md)); icPmin_md = VOID
	allocate(jcPmin_md(npy_md)); jcPmin_md = VOID
	allocate(kcPmin_md(npz_md)); kcPmin_md = VOID
	allocate(icPmax_md(npx_md)); icPmax_md = VOID
	allocate(jcPmax_md(npy_md)); jcPmax_md = VOID
	allocate(kcPmax_md(npz_md)); kcPmax_md = VOID

	! - - x - -
	ncxl = ceiling(dble(ncx)/dble(npx_md))
	do n=1,npx_md
		icPmax_md(n) = n * ncxl
		icPmin_md(n) = icPmax_md(n) - ncxl + 1
	end do	

	! - - y - -
	ncy_md   = nint(yL_md/dy)
	ncy_mdonly = ncy_md - ncy_olap
	ncyP_md = ncy_md / npy_md
	olap_jmin_mdcoord = npy_md - floor(dble(ncy_olap)/dble(ncyP_md))	 
	do n = olap_jmin_mdcoord,npy_md
		jcPmax_md(n) = n * ncyP_md - ncy_mdonly
		jcPmin_md(n) = jcPmax_md(n) - ncyP_md + 1
		if (jcPmin_md(n).le.0) jcPmin_md(n) = 1
	end do  

	! - - z - -
	nczl = ceiling(dble(ncz)/dble(npz_md))
	do n=1,npz_md
		kcPmax_md(n) = n * nczl
		kcPmin_md(n) = kcPmax_md(n) - nczl + 1
	end do

	if (myid_world.eq.0) then

		write(6000+myid_world,*), ''
		write(6000+myid_world,*), '==========================================='
		write(6000+myid_world,*), '------------ M D   M A P ------------------'
		write(6000+myid_world,*), '==========================================='
		write(6000+myid_world,*), 'npx_md = ', npx_md
		write(6000+myid_world,*), 'ncx    = ', ncx
		write(6000+myid_world,*), 'ncxl   = ', ncxl
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), '  icoord_md     icPmin_md     icPmax_md    '
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
		write(6000+myid_world,*), '  jcoord_md     jcPmin_md       jcPmax_md  '
		write(6000+myid_world,*), '-------------------------------------------'
		do n = 1,npy_md	
			write(6000+myid_world,'(1x,3i11)'), n, jcPmin_md(n), jcPmax_md(n)
		end do
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), 'npz_md = ', npz_md
		write(6000+myid_world,*), 'ncz    = ', ncz
		write(6000+myid_world,*), 'nczl   = ', nczl
		write(6000+myid_world,*), '-------------------------------------------'
		write(6000+myid_world,*), '  kcoord_md     kcPmin_md       kcPmax_md  '
		write(6000+myid_world,*), '-------------------------------------------'
		do n=1,npz_md
			write(6000+myid_world,'(1x,3i11)'), n, kcPmin_md(n), kcPmax_md(n)
		end do
		write(6000+myid_world,*), '-------------------------------------------'

	endif

end subroutine get_md_cell_ranges

!------------------------------------------------------------
!Calculate processor overlap between CFD/MD on all processors

subroutine get_overlap_blocks
	use coupler_module
	implicit none

	integer 			:: i,n,endproc,nolapsx,nolapsy,nolapsz
	integer				:: yLl_md,yLl_cfd
	integer,dimension(3):: pcoords

	!Get cartesian coordinate of overlapping md cells & cfd cells
	allocate(cfd_icoord2olap_md_icoords(npx_cfd,npx_md/npx_cfd)) 
	allocate(cfd_jcoord2olap_md_jcoords(npy_cfd,npy_md/npy_cfd)) 
	allocate(cfd_kcoord2olap_md_kcoords(npz_cfd,npz_md/npz_cfd)) 
	cfd_icoord2olap_md_icoords = VOID
	cfd_jcoord2olap_md_jcoords = VOID
	cfd_kcoord2olap_md_kcoords = VOID

	! - - x - -
	nolapsx = nint(dble(npx_md)/dble(npx_cfd))
	do n = 1,npx_cfd
	do i = 1,nolapsx	
		cfd_icoord2olap_md_icoords(n,i) = (n-1)*nolapsx + i
	end do
	end do

	! - - y - -
	yL_olap = (jcmax_olap - jcmin_olap + 1) * dy
	yLl_md  = yL_md/npy_md
	yLl_cfd = yL_cfd/npy_cfd
	nolapsy = ceiling(yL_olap/yLl_md)
	endproc = ceiling(yL_olap/yLl_cfd)
	do n = 1,endproc
	do i = 1,nolapsy
			cfd_jcoord2olap_md_jcoords(n,i) = (n-1)*nolapsy + i &
											  + (npy_md - nolapsy)
	end do
	end do

	! - - z - -
	nolapsz = nint(dble(npz_md)/dble(npz_cfd))
	do n = 1,npz_cfd
	do i = 1,nolapsz	
		cfd_kcoord2olap_md_kcoords(n,i) = (n-1)*nolapsz + i
	end do
	end do

	!nolaps = nolapsx*nolapsy*nolapsz

	if(myid_world.eq.0) then 
		
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

	endif
			
end subroutine get_overlap_blocks

subroutine intersect_comm
	use coupler_module
	use mpi
	implicit none

	integer :: n
	integer,dimension(3)   :: pcoords

	!Create communicator for all intersecting processors
	if (realm .eq. cfd_realm) then
		map%n = npx_md/npx_cfd
		allocate(map%rank_list(map%n))
		!Get rank(s) of overlapping MD processor(s)
		do n = 1,map%n
			pcoords(1)=cfd_icoord2olap_md_icoords(rank2coord_cfd(1,rank_realm),n)
			pcoords(2)=cfd_jcoord2olap_md_jcoords(rank2coord_cfd(2,rank_realm),n)
			pcoords(3)=cfd_kcoord2olap_md_kcoords(rank2coord_cfd(3,rank_realm),n)
			if (any(pcoords(:).eq.VOID)) then
				map%n = 0; map%rank_list(:) = VOID
			else
				map%rank_list(n) = coord2rank_md(pcoords(1),pcoords(2),pcoords(3))
			endif
			write(250+rank_realm,'(2a,6i5)'), 'overlap',realm_name(realm),rank_realm,map%n,map%rank_list(n),pcoords
		enddo
	else if (realm .eq. md_realm) then
		map%n = 1
		allocate(map%rank_list(map%n))	
		!Get rank of overlapping CFD processor
		pcoords(1) = rank2coord_md(1,rank_realm)*(dble(npx_cfd)/dble(npx_md))
		pcoords(2) = npy_cfd !rank2coord_md(2,rank_realm)*(dble(npy_cfd)/dble(npy_md))
		pcoords(3) = rank2coord_md(3,rank_realm)*(dble(npz_cfd)/dble(npz_md))
		map%rank_list(1) = coord2rank_cfd(pcoords(1),pcoords(2),pcoords(3))
		write(300+rank_realm,'(2a,6i5)'), 'overlap',realm_name(realm),rank_realm,map%n,map%rank_list(1),pcoords
	endif

end subroutine intersect_comm

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
	
	allocate(olap_mask(nproc_world))

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

subroutine print_overlap_comms
	use coupler_module
	use mpi
	implicit none

	integer :: trank

	if (myid_world.eq.0) then
		write(7500+rank_realm,*), ''
		write(7500+rank_realm,*), '----------- OVERLAP COMMS INFO ------------'
		write(7500+rank_realm,*), '-------------------------------------------'
		write(7500+rank_realm,*), '        RANKS              BROADCAST TEST  '
		write(7500+rank_realm,*), '  world  realm  olap      testval( = group)'
		write(7500+rank_realm,*), '-------------------------------------------'
	end if
	
	do trank = 1,nproc_world
		if (rank_world.eq.trank) then
			write(7500+rank_realm,'(3i7,i16)'), rank_world,rank_realm, &
								rank_olap, testval 	
		end if
	end do

	if (myid_world.eq.0) then
		write(7500+rank_realm,*), '-------- END OVERLAP COMMS INFO  ----------'
		write(7500+rank_realm,*), '==========================================='
	end if
	
end subroutine print_overlap_comms

!=============================================================================
!	Adjust CFD domain size to an integer number of lattice units used by  
!	MD if sizes are given in sigma units
!-----------------------------------------------------------------------------

subroutine coupler_cfd_adjust_domain(xL, yL, zL, nx, ny, nz, density_cfd)
    use mpi
    use coupler_input_data
    use coupler_module, only : b => MD_initial_cellsize, CPL_REALM_COMM, rank_realm, ierr
    implicit none

    integer, optional, intent(inout) 			:: nx, ny, nz
    real(kind(0.d0)), optional, intent(inout) 	:: xL,yL,zL
    real(kind(0.d0)), optional, intent(inout) 	:: density_cfd

    ! Internal variables
    integer										:: ierror, root
    real(kind=kind(0.d0)), pointer 				:: xyz_ptr => null()
	character(1)				   				:: direction
	logical										:: changed

	!Define root processes
	root = 1

    ! Local rank, useful for messeges
    if (density_tag == CPL) then 
        density_cfd = density
    else
        density     = density_cfd
        density_tag = CFD
    endif

	! Coupler input parameters are set in CFD 
    if ( cfd_coupler_input%domain%tag == VOID) then
        cfd_coupler_input%domain%tag = CFD
	end if

	print*, 'cell type =', cfd_coupler_input%domain%cell_type

    select case (cfd_coupler_input%domain%cell_type)
    case ("FCC","Fcc","fcc")
		b = (4.d0/density)**(1.d0/3.d0)
    case default
		write(*,*) "Wrong unit cell type in coupler_cfd_adjust_domain. Stopping ... "
		ierror = COUPLER_ERROR_INIT
 		call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
    end select

	! Check CFD domain and MD domain are compatible sizes to allow a
	! stable initial MD lattice structure - resize if possible or
	! stop code and demand a regeneration of grid if vary by more than 0.01
	changed = .false.
    if (present(xL)) then
 		xyz_ptr => cfd_coupler_input%domain%x
        call init_length(xyz_ptr,xL,resize=.true.,direction='x',print_warning=changed)
    endif

    if (present(yL)) then
		! No need to adjust y because we can adjust DY in MD to
		! have an integer number of FCC units.
		xyz_ptr => cfd_coupler_input%domain%y
		call init_length(xyz_ptr,yL,resize=.false.,direction='y',print_warning=changed)
    endif

    if (present(zL)) then
        xyz_ptr => cfd_coupler_input%domain%z
        call init_length(xyz_ptr,zL,resize=.true.,direction='z',print_warning=changed)
    endif

	if (changed .and. cfd_code_id .ne. couette_serial) then
		print*, "Regenerate Grid with corrected sizes as above"
		call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
	endif

    ! set CFD number of cells
    if (present(nx)) then
        if (cfd_coupler_input%ncells%tag == CPL ) then
            nx = cfd_coupler_input%ncells%x
        else
            if(rank_realm .eq. root) then
                write(0,*)"WARNING: nx is present in coupler_cfd_adjust_domain argument list"
                write(0,*)"         but coupler_input%ncells tag is void."
                write(0,*)"         Using CFD input value!"
            endif
        endif

    endif

    if (present(ny)) then
        if (cfd_coupler_input%ncells%tag == CPL ) then
            ny = cfd_coupler_input%ncells%y
        else
            if(rank_realm .eq. root) then
                write(0,*)"WARNING: ny is present in coupler_cfd_adjust_domain argument list"
                write(0,*)"         but coupler_input%ncells tag is void."
                write(0,*)"         Using CFD input value!"
            endif
        endif
    endif
        

    if (present(nz)) then
        if (cfd_coupler_input%ncells%tag == CPL ) then
            nz = cfd_coupler_input%ncells%z
        else
            if(rank_realm .eq. root) then
                write(0,*)"WARNING: nz is present in coupler_cfd_adjust_domain argument list"
                write(0,*)"         but coupler_input%ncells tag is void."
                write(0,*)"         Using CFD input value!"
            endif
        endif
    endif

    ! check id CFD cell sizes are larger than 2*sigma 
    call test_cfd_cell_sizes

contains

!-----------------------------------------------------------------------------

subroutine init_length(rin,rout,resize,direction,print_warning)
	implicit none
            
    real(kind=kind(0.d0)), intent(in)    :: rin
    real(kind=kind(0.d0)), intent(inout) :: rout
    logical, intent(in)                  :: resize
	character(*),intent(in) 			 :: direction
    logical,intent(out)					 :: print_warning

    real(kind=kind(0.d0)) :: rinit  ! store the initial value of rout or rin needed for print


    print_warning=.false.

    select case (cfd_coupler_input%domain%tag)
    case (CPL)
        select case (cfd_coupler_input%domain%units ) 
        case ("CELLSIDE","CellSide","Cellside","cellside")
            rout = b * rin
        case("SIGMA", "Sigma", "sigma")
            select case (cfd_coupler_input%domain%cell_type)
            case("FCC","Fcc","fcc")
                if (resize) then 
                rinit = rin 
                rout = real(floor(rin/b),kind(0.d0))*b
                print_warning = .true.
                endif
            case default
                write(*,*) "wrong unit cell type in coupler_cfd_adjust_domain. Stopping ... "
                ierror = COUPLER_ERROR_INIT
                call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
            end select
        end select
    case (CFD) 
        if(resize) then
            rinit = rout
            rout = real(nint(rout/b),kind(0.d0))*b
            print_warning = .true.
        endif
    case default
        write(*,*) "Wrong domain tag in coupler_cfd_adjust_domain. Stopping ... "
        ierror = COUPLER_ERROR_INIT
        call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
    end select

    if (print_warning) then 
        if (rank_realm .eq. root) then 
            write(*,'(3(a,/),3a,/,2(a,f20.10),/a,/,a)') &
                    "*********************************************************************",	&
                    "WARNING - this is a coupled run which resets CFD domain size         ",	&
                    " to an integer number of MD initial cells:		                      ", 	&
                    "	Domain resized in the ", direction, " direction			          ",	&
                    " inital size =", rinit, " resized ", rout,									&
                    "								                                      ",	& 
                    "*********************************************************************"   
        endif
		!If resize is insignificant then return flag print_warning as false
		if (abs(rinit-rout) .lt. 0.01) print_warning = .false.
    end if
    
end subroutine init_length

!-----------------------------------------------------------------------------


subroutine test_cfd_cell_sizes
    implicit none

    integer, pointer :: ndim => null()

    ndim => cfd_coupler_input%domain%ndim

    if (rank_realm .eq. root) then
        if ((present(xL) .and. present(nx)) .or. &
            (cfd_coupler_input%domain%tag == CPL .and. &
            cfd_coupler_input%ncells%tag == CPL)) then
            if (xL/nx < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in x direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          xL=",xL,"nx=",nx
            endif
        endif

        if ((present(yL) .and. present(ny)) .or. &
            (cfd_coupler_input%domain%tag == CPL .and. & 
            cfd_coupler_input%ncells%tag == CPL .and. ndim > 1)) then
            if (yL/ny < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in y direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          yL=",yL,"nx=",ny
            endif
        endif

        if ((present(zL) .and. present(nz)) .or. &
            (cfd_coupler_input%domain%tag == CPL .and. &
            cfd_coupler_input%ncells%tag == CPL .and. ndim > 2 )) then
            if (zL/nz < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in z direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          zL=",zL,"nx=",nz
            endif
        endif
    end if
        
end subroutine test_cfd_cell_sizes
            
end subroutine coupler_cfd_adjust_domain

!=============================================================================
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!
!							SIMULATION
!
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!=============================================================================




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

		integer :: pos, ixyz, icell, jcell, kcell

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

		integer :: maxi,mini,maxj,minj,maxk,mink
		integer :: ncxl,ncyl,nczl
		integer :: ncells
		integer :: icell,jcell,kcell
		integer :: coord(3),extents(6)
		integer :: bufsize
		integer :: trank_olap, tid_olap, trank_world, trank_cart

		if (realm.eq.cfd_realm) then
			call CPL_Cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
			call CPL_olap_extents(coord,cfd_realm,extents,ncells)
			bufsize = npercell*ncells
		else
			bufsize = 0
		end if

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
		integer :: tid_olap
		integer :: ncxl, ncyl, nczl
		integer :: trank_world, trank_cart
		integer :: coord(3),cfdextents(6),mdextents(6)
		integer :: ixyz, icell, jcell, kcell
		integer :: i,j,k

		! Grab CFD proc extents to be used for local cell mapping (when an 
		! array is an input to subroutine it has lower bound 1 by default)
		call CPL_cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,cfdextents)

		! Loop over procs in olap comm and pack scatter buffer 
		! in separate regions for each MD proc
		pos = 1
		do n = 1,nproc_olap

			tid_olap = n - 1
			if (tid_olap.eq.CFDid_olap) cycle

			call CPL_Cart_coords(CPL_OLAP_COMM,n,md_realm,3,coord,ierr)
			call CPL_proc_extents(coord,md_realm,mdextents)
			ncxl = mdextents(2) - mdextents(1) + 1
			ncyl = mdextents(4) - mdextents(3) + 1
			nczl = mdextents(6) - mdextents(5) + 1

			do ixyz = 1,npercell
			do icell= mdextents(1),mdextents(2)
			do jcell= mdextents(3),mdextents(4)
			do kcell= mdextents(5),mdextents(6)
				i = icell - cfdextents(1) + 1
				j = jcell - cfdextents(3) + 1
				k = kcell - cfdextents(5) + 1
				scatterbuf(pos) = scatterarray(ixyz,i,j,k)
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

		pos = 1
		do ixyz = 1,npercell
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
		call error_abort('Wrong realm in CPL_proc_extents')
	end select

	if (present(ncells)) then
		ncells = (extents(2) - extents(1) + 1) * &
				 (extents(4) - extents(3) + 1) * &
				 (extents(6) - extents(5) + 1)
	end if

end subroutine CPL_proc_extents

subroutine CPL_olap_extents(coord,realm,extents,ncells)
	use mpi
	use coupler_module, only: md_realm,      cfd_realm,      &
	                          icPmin_md,     icPmax_md,      &
	                          jcPmin_md,     jcPmax_md,      &
	                          kcPmin_md,     kcPmax_md,      &
	                          icPmin_cfd,    icPmax_cfd,     &
	                          jcPmin_cfd,    jcPmax_cfd,     &
	                          kcPmin_cfd,    kcPmax_cfd,     &
	                          cfd_icoord2olap_md_icoords,    &
	                          cfd_jcoord2olap_md_jcoords,    &
	                          cfd_kcoord2olap_md_kcoords,    &
	                          icmin_olap,                    & 
	                          icmax_olap,                    & 
	                          jcmin_olap,                    & 
	                          jcmax_olap,                    & 
	                          kcmin_olap,                    & 
	                          kcmax_olap,                    & 
	                          error_abort
	implicit none

	integer, intent(in)  :: coord(3), realm
	integer, intent(out) :: extents(6)
	integer, optional, intent(out) :: ncells
!	integer :: mini, maxi, minj, maxj, mink, maxk

	select case(realm)
	case(md_realm)
		!call CPL_proc_extents(coord,md_realm,extents)
		extents(1) = max(icPmin_md(coord(1)),icmin_olap)
		extents(2) = min(icPmax_md(coord(1)),icmax_olap) 
		extents(3) = max(jcPmin_md(coord(2)),jcmin_olap) 
		extents(4) = min(jcPmax_md(coord(2)),jcmax_olap) 
		extents(5) = max(kcPmin_md(coord(3)),kcmin_olap) 
		extents(6) = min(kcPmax_md(coord(3)),kcmax_olap) 
	case(cfd_realm)
		!maxi = maxval(cfd_icoord2olap_md_icoords(coord(1),:))
		!mini = minval(cfd_icoord2olap_md_icoords(coord(1),:))
		!maxj = maxval(cfd_jcoord2olap_md_jcoords(coord(2),:))
		!minj = minval(cfd_jcoord2olap_md_jcoords(coord(2),:))
		!maxk = maxval(cfd_kcoord2olap_md_kcoords(coord(3),:))
		!mink = minval(cfd_kcoord2olap_md_kcoords(coord(3),:))
		!extents(1) = icPmin_md(mini)
		!extents(2) = icPmax_md(maxi)
		!extents(3) = jcPmin_md(minj)
		!extents(4) = jcPmax_md(maxj)
		!extents(5) = kcPmin_md(mink)
		!extents(6) = kcPmax_md(maxk)
		extents(1) = max(icPmin_cfd(coord(1)),icmin_olap)
		extents(2) = min(icPmax_cfd(coord(1)),icmax_olap) 
		extents(3) = max(jcPmin_cfd(coord(2)),jcmin_olap) 
		extents(4) = min(jcPmax_cfd(coord(2)),jcmax_olap) 
		extents(5) = max(kcPmin_cfd(coord(3)),kcmin_olap) 
		extents(6) = min(kcPmax_cfd(coord(3)),kcmax_olap) 
	case default
		call error_abort('Wrong realm in CPL_olap_extents')
	end select

	if (present(ncells)) then
		ncells = (extents(2) - extents(1) + 1) * &
				 (extents(4) - extents(3) + 1) * &
				 (extents(6) - extents(5) + 1)
	end if

end subroutine CPL_olap_extents

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

!=============================================================================
! coupler_send_data wrapper for 3d arrays
! see coupler_send_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_send_data_3d(temp,index_transpose,asend_lbound,&
    asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)
    implicit none

    integer , optional, intent(in) :: glower(3), gupper(3) 
    integer, optional, intent(in) :: asend_lbound(3)       
    integer, optional, intent(in) :: asend_grid_start(3)   
    integer, optional, intent(in) :: asend_grid_end(3)     
    integer, optional, intent(in) :: index_transpose(3)    
    logical, optional, intent(in) :: use_overlap_box(3)
    real(kind=kind(0.d0)),dimension(:,:,:), intent(in) :: temp
    
    integer n1,n2,n3,n4
    real(kind=kind(0.d0)),dimension(:,:,:,:),allocatable :: asend

	!print*,'3DTemp size', size(temp,1),size(temp,2),size(temp,3)
    
    n1 = 1
    n2 = size(temp,1)
    n3 = size(temp,2)
    n4 = size(temp,3)

	!Add padding column to 3D array to make it 4D
	allocate(asend(n1,n2,n3,n4))
	asend(1,:,:,:) = temp(:,:,:)
  
    call coupler_send_data_xd(asend,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)

end subroutine coupler_send_data_3d

!=============================================================================
! coupler_send_data wrapper for 4d arrays
! see coupler_send_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_send_data_4d(asend,index_transpose,asend_lbound,asend_grid_start,&
    asend_grid_end,glower,gupper,use_overlap_box)
    implicit none

    integer , optional, intent(in) :: glower(3), gupper(3)
    integer, optional, intent(in)  :: asend_lbound(3)      
    integer, optional, intent(in)  :: asend_grid_start(3)  
    integer, optional, intent(in)  :: asend_grid_end(3)    
    integer, optional, intent(in)  :: index_transpose(3)   
    logical, optional, intent(in)  :: use_overlap_box(3)
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in) :: asend
    
    integer n1,n2,n3,n4
 
    call coupler_send_data_xd(asend,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)

end subroutine coupler_send_data_4d


!=============================================================================
! Send data from the local grid to the associated ranks from the other 
! realm
! The entities of this operation are as follow:
!  1) the local box of the global grid held by this rank
!  2) the data associate with a subset of local box which are held in some region of asend
!  3) the associated ranks from the other realm which overlap the same region of
!     the global grid (described in map)
!
!  In index range we have as follow:
!  the input array has the limits
!  a(n,ias:iae,jas:jae,kas:kae)
!  contains local grid data in the sector
!  a(n,iags:iage,jags:jage,kags:kage)
! 
!  which corresponds to the global index region
!  igs:ige jgs:jge, kgs:kge
!
!  which must be contain in the local grid box
!  ibs:ibs, jbs:jbe, kb:kbe 
!
!  so far the transfers done correspond to ias=iags=igs=ibs, etc 
!
! Main steps of the algorithm:
!  1) Using the arguments one must find the range of the grid data in asend
!     that is exclude halos, etc
!     This can be done by passing array sector of asend 
!     e.g. asend(igrid_start:igrid_end,jgrid_start_j:grid_end,kgrid_start:kgrid_end)
!     So far only this method was tested
!  2) Use the coupler map to find which sector of the grid go to which rank from the other
!     realm
!  3) send a quick description of the data ( an array of 8 integers )
!  4) send the data if there is anything to send
!-----------------------------------------------------------------------------

subroutine coupler_send_data_xd(asend,index_transpose,asend_lbound,&
    asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)
	use mpi
	use coupler_module
	use coupler_internal_cfd, only: bbox_cfd
    use coupler_internal_md, only : bbox_md => bbox
	implicit none
 
    ! array containing data distributed on the grid
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in) :: asend
    ! lower and upper corners of the grid to which the data of asend
    ! must be associated
    ! if not present the local grid box limits are used
    ! any negative value is discared for the default
    integer , optional, intent(in) :: glower(3), gupper(3)
    ! starting index of asend if is not (1,1,1) 
    integer, optional, intent(in) :: asend_lbound(3)
    ! starting index coordinated of grid data if is not (1,1,1)      
    integer, optional, intent(in) :: asend_grid_start(3)   
    ! upper corner of grid data stored in asend
    ! if not present the the local grid box limits are  used
    integer, optional, intent(in) :: asend_grid_end(3)
    ! specify the order of dimensions in asend
    ! default is ( x, y, z) but e.g. cartesian velocity
    ! in CFD solver uses (z,x,y)     
    integer, optional, intent(in) :: index_transpose(3)    ! rule of transposition for the coordinates 
    ! whether to use on MD side the extended grid of the map ( which
    ! has a halo of CFD cells around MD local box)
    ! this is useful to trap molecules that left 
    ! the local domain or if the CFD cells sixe is not
    ! conmensurate with sizes of MD local box
    logical, optional, intent(in) :: use_overlap_box(3)

	!Number of halos
	integer	:: nh = 1    

    ! local indices 
    integer :: ix,iy,iz,bicmin,bicmax,bjcmin,bjcmax,bkcmin,bkcmax
	integer	:: micmin,micmax,mjcmin,mjcmax,mkcmin,mkcmax
	integer	:: aicmin,ajcmin,akcmin,at(2,3),ig(2,3)
    ! auxiliaries 
    integer	:: i,ndata,itag, dest, req(map%n), vel_indx(8,map%n), &
		        start_address(map%n+1), a_components
	real(kind=kind(0.d0)), allocatable :: vbuf(:)
	integer, save :: ncalls = 0

	! This local CFD domain is outside MD overlap zone 
	if ( map%n .eq. 0 ) return 

	ncalls = ncalls + 1
    ! Get local grid box ranges seen by this rank for either CFD or MD
    if (realm .eq. cfd_realm) then 
		!Load CFD cells per processor
        bicmin = icPmin_cfd(iblock_realm)
        bicmax = icPmax_cfd(iblock_realm) 
        bjcmin = jcPmin_cfd(jblock_realm) 
        bjcmax = jcPmax_cfd(jblock_realm) 
        bkcmin = kcPmin_cfd(kblock_realm) 
        bkcmax = kcPmax_cfd(kblock_realm) 
    elseif (realm .eq. md_realm) then 
        ! Load MD cells per processor
        bicmin = icPmin_md(iblock_realm)
        bicmax = icPmax_md(iblock_realm) 
        bjcmin = jcPmin_md(jblock_realm) 
        bjcmax = jcPmax_md(jblock_realm) 
        bkcmin = kcPmin_md(kblock_realm) 
        bkcmax = kcPmax_md(kblock_realm) 
		! use_overlap_box includes halos on the MD side
        if(present(use_overlap_box)) then 
            if (use_overlap_box(1)) then
                bicmin = icPmin_md(iblock_realm)-nh
                bicmax = icPmax_md(iblock_realm)-nh
            endif
            if(use_overlap_box(2)) then
                bjcmin = jcPmin_md(jblock_realm)-nh
                bjcmax = jcPmax_md(jblock_realm)-nh
            endif
            if(use_overlap_box(3)) then
                bkcmin = kcPmin_md(kblock_realm)-nh
                bkcmax = kcPmax_md(kblock_realm)-nh
            endif
       endif
	endif

    ! Re-order indices of x,y,z coordinates (handles transpose arrays)
    ix = 1 ; iy = 2; iz = 3
    if( present(index_transpose)) then
        ix = index_transpose(1)
        iy = index_transpose(2)
        iz = index_transpose(3)
    endif

	! If asend is less than number of points in domain
	! change amount of data to send
    bicmin = bicmin
    bicmax = min(bicmin + size(asend,ix+1) - 1,bicmax)
    bjcmin = bjcmin
    bjcmax = min(bjcmin + size(asend,iy+1) - 1,bjcmax)
    bkcmin = bkcmin
    bkcmax = min(bkcmin + size(asend,iz+1) - 1,bkcmax)

    ! Warning if asend goes over local grid box
	if (bicmax-bicmin .ne. size(asend,ix+1)) & 
		write(0,'(3(a,i8),a)') "Proc=",rank_realm, " Sent datasize ",size(asend,ix+1), & 
								" not equal to gridsize ",bicmax-bicmin," in x" 
	if (bjcmax-bjcmin .ne. size(asend,iy+1)) & 
		write(0,'(3(a,i8),a)') "Proc=",rank_realm, " Sent datasize ", size(asend,iy+1), & 
								" not equal to gridsize ",bjcmax-bjcmin," in y" 
	if (bkcmax-bkcmin .ne. size(asend,iz+1)) & 
		write(0,'(3(a,i8),a)') "Proc=",rank_realm, " Sent datasize ",size(asend,iz+1), & 
								" not equal to gridsize ",bkcmax-bkcmin," in z"  

    ! grid data in asend data can be mapped to a specific region
    ! of the local grid -- negative values are discarded in favor of defaults
    if (present(glower))then 
        if (glower(1) > 0) bicmin = glower(1)
        if (glower(2) > 0) bjcmin = glower(2)
        if (glower(3) > 0) bkcmin = glower(3)
    endif

    !  upper limit exceed for grid data point stored in asend
    if (present(gupper))then 
        if (gupper(1) > 0) bicmax = gupper(1)
        if (gupper(2) > 0) bjcmax = gupper(2)
        if (gupper(3) > 0) bkcmax = gupper(3)
    endif

    ! sanity check
    if (present(gupper) .and. present(glower))then 
		if (glower(1) .gt. gupper(1)) print*, 'Proc=',rank_realm, 'Lower bounds of send', & 
								   	glower(1),'greater than upper',gupper(1)
		if (glower(2) .gt. gupper(2)) print*,'Proc=',rank_realm, 'Lower bounds of send', & 
								  	glower(2),'greater than upper',gupper(2)
		if (glower(3) .gt. gupper(3)) print*,'Proc=',rank_realm, 'Lower bounds of send', & 
									glower(3),'greater than upper',gupper(3)
	endif

    ! Array indices in asend for the data mapped on grid
    ! +1 shift is because asend grid indices start from 2
    aicmin = 1
    ajcmin = 1
    akcmin = 1
    if (present(asend_lbound)) then 
        if (.not. present(asend_grid_start))then 
            write(0,*) "because the asend lower bound is not default asend_grid_start argument must be provided"
            call MPI_Abort(MPI_COMM_WORLD,COUPLER_ABORT_SEND_CFD,ierr)
        endif
        aicmin = asend_grid_start(ix) - asend_lbound(ix) + 1
        ajcmin = asend_grid_start(iy) - asend_lbound(iy) + 1
        akcmin = asend_grid_start(iz) - asend_lbound(iz) + 1
    endif

    ! insanity checks 

    ! Number of components at each grid point
    a_components = size(asend,dim=1)
   
    ! loop through the maps and send the corresponding sections of asend
    do i = 1, map%n
        dest = map%rank_list(i)

		! map%domains contains the limits of the overlap - 
		! if all cells do not overlap in the x-z plane then the
		! limits defined in map domains are used instead
        micmin = max(bicmin, map%domains(1,i)) 
        micmax = min(bicmax, map%domains(2,i))
        mjcmin = max(bjcmin, map%domains(3,i))
        mjcmax = min(bjcmax, map%domains(4,i))
        mkcmin = max(bkcmin, map%domains(5,i)) 
        mkcmax = min(bkcmax, map%domains(6,i))	

        ! Amount of data to be sent
        ndata = a_components * (micmax-micmin + 1) * (mjcmax-mjcmin + 1) * (mkcmax-mkcmin + 1)

        if ( ndata > 0) then 
            if(allocated(vbuf)) deallocate(vbuf)
            allocate(vbuf(ndata))
            
            ! Location in asend of the domain to be sent
            ig(1,ix) = micmin - bicmin + aicmin 
            ig(2,ix) = micmax - bicmax + aicmin
            ig(1,iy) = mjcmin - bjcmin + ajcmin
            ig(2,iy) = mjcmax - bjcmax + ajcmin
            ig(1,iz) = mkcmin - bkcmin + akcmin
            ig(2,iz) = mkcmax - bkcmax + akcmin
 
            write(0,'(a,14i5)') ' coupler send a*, ig ...', ncalls, rank_realm,aicmin,ajcmin,akcmin,ig 

            vbuf(1:ndata) = reshape(asend(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)), (/ ndata /) )
        endif
        ! Attention ncall could go over max tag value for long runs!!
        itag = mod( ncalls, MPI_TAG_UB)
        ! send the info about the data to come
        call MPI_send((/ndata,a_components,micmin,micmax,mjcmin,mjcmax,mkcmin,mkcmax/),8,MPI_INTEGER,&
            			dest, itag, CPL_INTER_COMM, ierr)
        ! send data only if there is anything to send
        if (ndata > 0) then 
			!vbuf = 1.234567891011121314151617d0
			!print'(2a,2i8,4f25.16)', 'ICP send data',realm_name(realm), rank_realm, & 
			!							 size(vbuf), maxval(vbuf),minval(vbuf),sum(vbuf),vbuf(10)
            call MPI_send(vbuf, ndata, MPI_DOUBLE_PRECISION, dest, itag, CPL_INTER_COMM, ierr)
        endif
    enddo

end subroutine coupler_send_data_xd

!=============================================================================
! coupler_recv_data wrapper for 3d arrays
! see coupler_recv_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_recv_data_3d(temp,index_transpose,asend_lbound,asend_grid_start,&
    asend_grid_end,glower,gupper,accumulate,pbc)
    implicit none

    integer, optional, intent(in) :: asend_lbound(3)
    integer, optional, intent(in) :: asend_grid_start(3) 
    integer, optional, intent(in) :: asend_grid_end(3)   
    integer, optional, intent(in) :: index_transpose(3)  
    integer, optional, intent(in) :: glower(3), gupper(3)
    logical, optional, intent(in) :: accumulate          
    integer, optional, intent(in) :: pbc  
    real(kind(0.d0)),dimension(:,:,:),intent(inout) :: temp 
                                                        
	integer n1, n2, n3, n4
	real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: arecv    

    n1 = 1 
    n2 = size(temp,1)
    n3 = size(temp,2)
    n4 = size(temp,3)

	!Add padding column to 3D array to make it 4D
	allocate(arecv(n1,n2,n3,n4))
    call coupler_recv_data_xd(arecv,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,accumulate,pbc)
	temp(:,:,:) = 	arecv(1,:,:,:) 

end subroutine coupler_recv_data_3d

!=============================================================================
! coupler_recv_data wrapper for 4d arrays
! see coupler_recv_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_recv_data_4d(arecv,index_transpose,asend_lbound,&
    asend_grid_start,asend_grid_end,glower,gupper,accumulate,pbc)

    implicit none

    integer, optional, intent(in) :: asend_lbound(3)     
    integer, optional, intent(in) :: asend_grid_start(3) 
    integer, optional, intent(in) :: asend_grid_end(3)   
    integer, optional, intent(in) :: index_transpose(3)  
    integer, optional, intent(in) :: glower(3), gupper(3)
    logical, optional, intent(in) :: accumulate          
    integer, optional, intent(in) :: pbc  
    real(kind(0.d0)),dimension(:,:,:,:),intent(inout) :: arecv

    call coupler_recv_data_xd(arecv,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,accumulate,pbc)

end subroutine coupler_recv_data_4d

!=============================================================================
! Receive data from to local grid from the associated ranks from the other 
! realm
! The entities of this operation are as follow:
!  1) the local box of the global grid held by this rank
!  2) the grid data region to be held in some region of arecv
!  3) the associated ranks from the other realm which overlap the same region of
!     the global grid
!
! Main steps of the algorithm:
!  1) Using the arguments find the range of the grid data in arecv
!     that is: exclude halos, etc
!     This can be done by passing array sector of arecv that is meant to hold grid data 
!     e.g. arecv(igrid_start:igrid_end,jgrid_start_j:grid_end,kgrid_start:kgrid_end)
!     So far only this method was tested
!  2) Use the coupler map to find where the arriving data should land in a temporary array
!  3) copy the whole lot, after the receives completes, to the corresponding location
!     of arecv
!
!  Notes: 
!  1) This subroutine can perform accumulation of data, if one chooses to 
!     send sum of cell quantities and respective numbers. Could be useful to 
!     collect contribution from molecule leaving local MD box.
!  2) Periodic boundary condition can be applied, useful for averages with
!     staggered boxes in which the end boxes collect only form half box
!     at the moment this works only with (z,x,y) array in z and x directions
!   
!    
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine coupler_recv_data_xd(arecv,index_transpose,a_lbound,&
    									a_grid_start,a_grid_end,glower,gupper,accumulate,pbc)
    use mpi
	use coupler_internal_cfd, only : bbox_cfd
    use coupler_internal_md, only :  bbox_md => bbox
	use coupler_module
    implicit none

    ! array that recieves grid distributed data 
    real(kind(0.d0)), dimension(:,:,:,:),intent(inout) :: arecv
    ! lower and upper corners of the grid to which the data of asend
    ! must be associated
    ! if not present the local grid box limits are used
    ! any negative value is discared for the default
    integer, optional, intent(in) :: glower(3), gupper(3)
    ! starting index of asend if is not (1,1,1)
    integer, optional, intent(in) :: a_lbound(3) 
    ! starting index coordinated of grid data if is not (1,1,1)    
    integer, optional, intent(in) :: a_grid_start(3)
    ! upper corner of grid data stored in asend
    ! if not present the the local grid box limits are used
    integer, optional, intent(in) :: a_grid_end(3)  
    ! specifyle the order of dimensions in arecv
    ! default is ( x, y, z) but e.g. cartesian velocity
    ! in CFD solver uses (z,x,y)
    integer, optional, intent(in) :: index_transpose(3) 
    ! if present add the receiving value if they land in the same
    ! grid position
    logical, optional, intent(in) :: accumulate   
    ! if present applies periodic boundary condition
    ! posible values 1 for x direction, 3 for z direction
    integer, optional, intent(in) :: pbc                 
                                                         
    ! local indices 
    integer is,ie,js,je,ks,ke,bicmin,bicmax,bjcmin,bjcmax,bkcmin,bkcmax,ix,iy,iz,aicmin,aicmax,ajcmin,ajcmax,akcmin,akcmax,&
        at(2,3),ig(2,3)
    ! auxiliaries 
    integer i,j,k,ii,jj,kk,ndata,itag, source,req(map%n), vel_indx(8,map%n), &
        p1s,p1e,p2s,p2e,p3s,p3e,pt1s,pt1e,pt2s,pt2e,pt3s,pt3e, bgt(2,3),& 
        start_address(map%n+1), status(MPI_STATUS_SIZE,map%n)
    real(kind(0.d0)), allocatable ::  vbuf(:), atmp(:,:,:,:)
    integer, save :: ncalls = 0

	! This local CFD domain is outside MD overlap zone 
	if ( map%n .eq. 0 ) return 
 
	ncalls = ncalls + 1

    ! Local grid box
    ! Get local grid box ranges seen by this rank for either CFD or MD
    if (realm .eq. cfd_realm) then 
		!Load CFD cells per processor
        bicmin = icPmin_cfd(iblock_realm)
        bicmax = icPmax_cfd(iblock_realm) 
        bjcmin = jcPmin_cfd(jblock_realm) 
        bjcmax = jcPmax_cfd(jblock_realm)
        bkcmin = kcPmin_cfd(kblock_realm) 
        bkcmax = kcPmax_cfd(kblock_realm) 
    elseif (realm .eq. md_realm) then 
        ! Load MD cells per processor
        bicmin = icPmin_md(iblock_realm)
        bicmax = icPmax_md(iblock_realm) 
        bjcmin = jcPmin_md(jblock_realm) 
        bjcmax = jcPmax_md(jblock_realm) 
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

    ! grid data in asend data can be mapped to a specific region
    ! of the local grid -- negative values are discarded in favor of defaults
    if (present(glower))then 
        if (glower(1) > 0) bicmin = glower(1)
        if (glower(2) > 0) bjcmin = glower(2)
        if (glower(3) > 0) bkcmin = glower(3)
    endif

    ! sanity check is needed here

    !  upper limit exteed for grid data point stored in asend
    if (present(gupper))then 
        if (gupper(1) > 0) bicmax = gupper(1)
        if (gupper(2) > 0) bjcmax = gupper(2)
        if (gupper(3) > 0) bkcmax = gupper(3)
    endif

    ! sanity check is needed here

    ! put the mapped block limits in a transposed array
    bgt(1,ix) = bicmin
    bgt(2,ix) = bicmax
    bgt(1,iy) = bjcmin
    bgt(2,iy) = bjcmax
    bgt(1,iz) = bkcmin
    bgt(2,iz) = bkcmax

    ! Array indices in arecv for the data mapped on grid
    ! +1 shift is because asend grid indices start from 2
    aicmin = 1
    ajcmin = 1
    akcmin = 1
    if (present(a_lbound)) then 
        if (.not. present(a_grid_start))then 
            write(0,*) "because the arecv lower bound is not default as_grid_start argument must be provided"
            call MPI_Abort(MPI_COMM_WORLD,COUPLER_ABORT_SEND_CFD,ierr)
        endif
        aicmin = a_grid_start(1) - a_lbound(1) + 1
        ajcmin = a_grid_start(2) - a_lbound(2) + 1
        akcmin = a_grid_start(3) - a_lbound(3) + 1
    endif

	!Use smallest of expected grid size or receive array size
    aicmax = min(aicmin + (bicmax-bicmin),size(arecv,ix+1)) 
    ajcmax = min(ajcmin + (bjcmax-bjcmin),size(arecv,iy+1))
    akcmax = min(akcmin + (bkcmax-bkcmin),size(arecv,iz+1))

    ! sanity checks are needed 

    
    ! Store the transposition for grid boundaries limits
    at(1,ix) = aicmin
    at(2,ix) = aicmax
    at(1,iy) = ajcmin
    at(2,iy) = ajcmax
    at(1,iz) = akcmin
    at(2,iz) = akcmax

    ! First get info about what's comming
    do i = 1, map%n
        source =  map%rank_list(i)
        itag = mod( ncalls, MPI_TAG_UB)
        call MPI_irecv(vel_indx(1,i),8,MPI_Integer,source,itag,&
            			CPL_INTER_COMM,req(i),ierr)
    enddo
    call MPI_waitall(map%n, req, status, ierr)

    write(0,'(a,16i5)') 'coupler recv vel_indx ', vel_indx

    ndata = 0
    do i = 1, map%n
        ndata = ndata + max(vel_indx(1,i),0)
    enddo

	allocate(vbuf(ndata),stat=ierr)
    start_address(1) = 1 
    do i = 1, map%n
        source =  map%rank_list(i) 
        ndata = vel_indx(1,i)
        start_address(i+1) = start_address(i) + ndata
        if (ndata > 0) then 
            ! Attention ncall could go over max tag value for long runs!!
			itag = mod(ncalls, MPI_TAG_UB)
            call MPI_irecv(vbuf(start_address(i)),ndata,MPI_DOUBLE_PRECISION,source,itag,&
                				CPL_INTER_COMM,req(i),ierr)
        else 
            req(i) = MPI_REQUEST_NULL
        endif
    enddo

	!print*, 'BEFORE WAIT STATEMENT'
    call MPI_waitall(map%n, req, status, ierr)
	!print'(2a,2i8,4f25.16)', 'ICP recv data',realm_name(realm),rank_realm, & 
	!							 size(vbuf), maxval(vbuf),minval(vbuf),sum(vbuf),vbuf(10)
    !	write(0,*) 'MD getCFD vel wait',  rank_realm, id, i, source, ierr

    ! Allocate atmp corresponding to reunion of mapped grid
    ! This must be thought as landing area for data coming from the other realm
    ig(1,ix) = minval(vel_indx(3,:))
    ig(2,ix) = maxval(vel_indx(4,:))
    ig(1,iy) = minval(vel_indx(5,:))
    ig(2,iy) = maxval(vel_indx(6,:))
    ig(1,iz) = minval(vel_indx(7,:))
    ig(2,iz) = maxval(vel_indx(8,:))

    !write(0,*) 'coupler recv atmp sizes ', vel_indx(2,1),ig

    allocate(atmp(vel_indx(2,1), ig(1,1):ig(2,1), ig(1,2):ig(2,2), ig(1,3):ig(2,3)))
    atmp(:,:,:,:) = 0.d0

    do i=1,map%n
        
        if(vel_indx(1,i) <= 0) cycle
        
        ! Get array indices from the recieved grid data
        ig(1,ix) = vel_indx(3,i)
        ig(2,ix) = vel_indx(4,i)
        ig(1,iy) = vel_indx(5,i)
        ig(2,iy) = vel_indx(6,i)
        ig(1,iz) = vel_indx(7,i)
        ig(2,iz) = vel_indx(8,i)
        
        if (present(accumulate)) then
            !write(0,*) 'coupler recv ig', ig
            atmp(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)) = &
                atmp(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)) &
                + reshape(vbuf(start_address(i):start_address(i+1)-1), &
                (/ vel_indx(2,i), ig(2,1)-ig(1,1)+1, ig(2,2)-ig(1,2)+1, ig(2,3)-ig(1,3)+1 /))
        else
            atmp(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)) = &
                reshape(vbuf(start_address(i):start_address(i+1)-1), &
                (/ vel_indx(2,i), ig(2,1)-ig(1,1)+1, ig(2,2)-ig(1,2)+1, ig(2,3)-ig(1,3)+1 /))
        endif

        ! call MPI_Error_string(status(MPI_ERROR,i), err_string, len(err_string), ierr);
	    !  write(0,*) 'MD getCFD vel err, rank_realm, i ', rank_realm, i, trim(err_string) 
        !call MPI_get_count(status(1,i),MPI_double_precision,ib,ierr)
		! write(0,*) 'MD recv ', rank_realm, id, i, ib, ' DP'
    enddo

    ! Transfer data from the landing area to arecv 
    ! if the flag accumulate is present normalise
    ! data with the value stored on the last line of the first column
    ! and apply periodic boundary conditions

    ! Indices for overlapping window
    ! usually the landing area and arecv coincide
    ! some test and warning should be done here
    p1s = max(lbound(atmp,2) - bgt(1,1) + at(1,1),at(1,1))
    p1e = min(ubound(atmp,2) - bgt(1,1) + at(1,1),at(2,1))
    p2s = max(lbound(atmp,3) - bgt(1,2) + at(1,2),at(1,2))
    p2e = min(ubound(atmp,3) - bgt(1,2) + at(1,2),at(2,2))
    p3s = max(lbound(atmp,4) - bgt(1,3) + at(1,3),at(1,3))
    p3e = min(ubound(atmp,4) - bgt(1,3) + at(1,3),at(2,3))

    ! needs revision, I'm not sure that it works in general
    pt1s = p1s + bgt(1,1) - at(1,1)
    pt1e = p1e + bgt(1,1) - at(1,1)
    pt2s = p2s + bgt(1,2) - at(1,2)
    pt2e = p2e + bgt(1,2) - at(1,2)
    pt3s = p3s + bgt(1,3) - at(1,3)
    pt3e = p3e + bgt(1,3) - at(1,3)

    !write(0,*)' p1s ...', rank_realm, p1s,p1e,p2s,p2e,p3s,p3e, pt1s, pt1e, pt2s,pt2e,pt3s,pt3e

    if(present(accumulate)) then
        if ( present(pbc)) then
			call halos(pbc) 
            !call set_pbc(pbc)
        endif
        
        ndata = ubound(atmp,1)
        do k = p3s , p3e
 		do j = p2s , p2e
 		do i = p1s , p1e
			ii = i + bgt(1,1) - at(1,1)
			jj = j + bgt(1,2) - at(1,2)
			kk = k + bgt(1,3) - at(1,3)
 			if (atmp(ndata, ii,jj,kk) > 0.d0) then
				arecv(:,i,j,k) = atmp(1:ndata-1,ii,jj,kk)/atmp(ndata,ii,jj,kk)
			else
				arecv(:,i,j,k) = 0.d0
			endif
		end do
		end do
        end do
    else 
        arecv(:,p1s:p1e,p2s:p2e,p3s:p3e) = atmp(:,pt1s:pt1e,pt2s:pt2e,pt3s:pt3e) 
    endif


	print*,'Size of recieved data',ig(1,1),ig(2,1),ig(1,2),ig(2,2),ig(1,3),ig(2,3),map%n, & 
									start_address(i),start_address(i+1)-1, &
					                vel_indx(2,i),ig(2,1)-ig(1,1)+1,ig(2,2)-ig(1,2)+1,ig(2,3)-ig(1,3)+1, &
									p1s,p1e,p2s,p2e,p3s,p3e,pt1s,pt1e,pt2s,pt2e,pt3s,pt3e

	!print'(2a,2i8,4f25.16)', 'ICP trfr data',realm_name(realm),rank_realm, & 
	!						size(arecv), maxval(arecv),minval(arecv),sum(arecv),arecv(1,p1s,p2s,p3s)
            
contains

!-----------------------------------------------------------------------------
! Sets periodic boundary conditions
! adds together the sums of the first and last bins on the global grid (x and z directions)
! and puts sum in both boxes
! Note:
! attention coord array are shifted with 1 in CFD and MD
! to get correct result in MPI function such as mpi_cart_rank, subtract 1 
!-----------------------------------------------------------------------------

subroutine halos(pbc)
    use coupler_module, only : npx_md
    implicit none

    integer, intent(in) :: pbc
    real(kind(0.d0)), allocatable, dimension(:,:,:) :: x1, x2
    integer dest, ip(3), status(MPI_STATUS_SIZE)

    ! PBC works only for 3,2,1 coordinate transposition
    if ( .not. (ix == 2 .and. iy == 3 .and. iz == 1)) then
        write(0,*) " Periodic boundary condition not implemented for this layout of data"
        return
    endif

    ! Sanity check atmp must touch the boundary of the global grid in zin x directions
    ! I guess that the PBC are applied at the ends of global grid, check with CFD people
    select case(pbc)
    case(3)
        ! first array dimension which corresponds to z direction for velocity arrays 
        allocate(x1(size(atmp,1),size(atmp,3),size(atmp,4)),x2(size(atmp,1),size(atmp,3),size(atmp,4)))
        x1 =  atmp(:, akcmin, :, :)
        x2 =  atmp(:, akcmax, :, :)
        atmp(:, akcmin, :, :) =  x2 ! atmp(:, aks, :, :) + x2
        atmp(:, akcmax, :, :) =  x1 ! atmp(:, ake, :, :) + x1
    case(1)
        ! second array dimension which correponds to x direction
        allocate(x1(size(atmp,1),size(atmp,2),size(atmp,4)),x2(size(atmp,1),size(atmp,2),size(atmp,4)))
        if ( npx_md == 1 )then 
            ! no MPI communication needed  
            !write(0,*) 'coupler internal cfd pbc', npx_md
            x1 =  atmp(:, :,lbound(atmp,3), :)                 
            x2 =  atmp(:, :,ubound(atmp,3), :)
            atmp(:, :, lbound(atmp,3), :) = x2  !atmp(:, :, lbound(atmp,3), :) + x2
            atmp(:, :, ubound(atmp,3), :) = x1  !atmp(:, :, ubound(atmp,3), :) + x1  
        else 
            ip = rank2coord_md(: ,rank_realm)
            if (ip(1) == 1 .and. ip(2) == 1) then 
                x1 =  atmp(:, :,lbound(atmp,3),:)
                call MPI_cart_rank(CPL_CART_COMM, (/ npx_md-1, 0, ip(3)-1 /), dest, ierr)
                call MPI_sendrecv(x1,size(x1),MPI_DOUBLE_PRECISION,dest,1,x2,size(x2),MPI_DOUBLE_PRECISION,&
                    				dest,1,CPL_CART_COMM,status,ierr)
                atmp(:, :, lbound(atmp,3), :) = x2		!atmp(:, :, lbound(atmp,3), :) + x2
            else if (ip(1) == npx_md .and. ip(2) == 1) then 
                x2 =  atmp(:, :,ubound(atmp,3), :)
                call MPI_cart_rank(CPL_CART_COMM, (/ 0, 0, ip(3) - 1 /), dest, ierr)
                call MPI_sendrecv(x2,size(x2),MPI_DOUBLE_PRECISION,dest,1,x1,size(x1),MPI_DOUBLE_PRECISION,&
                    				dest,1,CPL_CART_COMM,status,ierr)
                atmp(:, :, ubound(atmp,3), :) =  x1		!atmp(:, :, ubound(atmp,3), :) + x1
            endif
        endif
    end select
end subroutine halos

end subroutine coupler_recv_data_xd

!============================================================================
!
! Utility functions and subroutines that extract parameters from internal modules 
!
!-----------------------------------------------------------------------------    

!-------------------------------------------------------------------------------
! return to the caller coupler parameters from cfd realm
!-------------------------------------------------------------------------------
subroutine coupler_cfd_get(jcmax_overlap,jcmin)
    use coupler_module, jcmax_overlap_ => ncy_olap,jcmin_=>jcmin
    implicit none

    integer,optional,intent(out) :: jcmax_overlap,jcmin

    if(present(jcmax_overlap)) then
        jcmax_overlap = jcmax_overlap_
    endif

    if(present(jcmin))then 
        jcmin = jcmin_
    end if

end subroutine coupler_cfd_get


!===============================================================================
! Return to the caller coupler parameters from MD realm
!-------------------------------------------------------------------------------
subroutine coupler_md_get(	xL_md,yL_md,zL_md, MD_initial_cellsize, top_dy, &
							overlap_with_continuum_force, overlap_with_top_cfd, ymin_continuum_force, &
    						ymax_continuum_force, xmin_cfd_grid, xmax_cfd_grid, zmin_cfd_grid, zmax_cfd_grid, &
    						dx_cfd, dz_cfd, cfd_box)
    use mpi
	use coupler_internal_md, only : bbox,cfd_box_ => cfd_box
	!use coupler_internal_md, only : xL_md_ =>  xL_md, yL_md_ =>  yL_md, zL_md_ => zL_md, &
	!	b => MD_initial_cellsize, x, y, z,j => jcmax_overlap_cfd, dx, dz, bbox, half_domain_lengths, &
    !    cfd_box_ => cfd_box
	use coupler_module, xL_md_=>xL_md,yL_md_=>yL_md,zL_md_=>zL_md, b => MD_initial_cellsize, j=> ncy_olap
	implicit none

	real(kind(0.d0)), optional, intent(out) :: xL_md, yL_md, zL_md, MD_initial_cellsize, top_dy,&
        ymin_continuum_force,ymax_continuum_force, xmin_cfd_grid, xmax_cfd_grid, &
		 zmin_cfd_grid, zmax_cfd_grid, dx_cfd, dz_cfd 
    logical, optional, intent(out) :: overlap_with_continuum_force, overlap_with_top_cfd
    type(cfd_grid_info), optional, intent(out) :: cfd_box

	real(kind(0.d0)),dimension(3)	:: half_domain_lengths

	half_domain_lengths = 0.5d0 * (/ xL_md_, yL_md_, zL_md_ /)

    if (present(xL_md)) then
		xL_md = xL_md_
	endif

	if (present(yL_md)) then 
		yL_md = yL_md_
	endif

	if (present(zL_md)) then 
		zL_md = zL_md_
	endif

	if (present(MD_initial_cellsize)) then 
		MD_initial_cellsize = b
	endif

	if (present(top_dy)) then
    	top_dy = yg(1,j) - yg(1,j-1)
	end if

    if(present(overlap_with_continuum_force)) then
        ! check if the countinuum force domain is contained in one 
        ! domain along y direction
        if ( (yg(1,j - 2) < bbox%bb(1,2) .and. &
              yg(1,j - 1) > bbox%bb(1,2)) .or. &
             (yg(1,j - 2) < bbox%bb(2,2) .and. &
              yg(1,j - 1) > bbox%bb(2,2))) then
            write(0,*) " the region in which the continuum constraint force is applied "
            write(0,*) " spans over two domains. This case is not programmed, please investigate"
            call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_CONTINUUM_FORCE,ierr)
        endif
        
        if ( yg(1,j - 1) <  bbox%bb(2,2) .and. yg(1,j - 2) >= bbox%bb(1,2) ) then
            overlap_with_continuum_force = .true.
        else   
            overlap_with_continuum_force = .false.
        endif

    endif

     if(present(overlap_with_top_cfd)) then
        ! check if the MD domain overlaps with the top of CFD grid (along y)
        ! the MD constrain force is applyied top layer of cfd cells
        if ( yg(1,j - 1) < bbox%bb(2,2) .and. yg(1,j - 1) >= bbox%bb(1,2) ) then
            overlap_with_top_cfd = .true.
        else   
            overlap_with_top_cfd = .false.
        endif

    endif
 
     if(present(ymin_continuum_force)) then
         ymin_continuum_force = yg(1,j - 2) - bbox%bb(1,2) - half_domain_lengths(2)
     endif
         
     if(present(ymax_continuum_force)) then
         ymax_continuum_force = yg(1,j - 1) - bbox%bb(1,2) - half_domain_lengths(2)
     endif 

     if(present(xmin_cfd_grid)) then
         xmin_cfd_grid = xg(bbox%is,1) - bbox%bb(1,1) - half_domain_lengths(1)
     endif

     if(present(xmax_cfd_grid)) then
         xmax_cfd_grid = xg(bbox%ie,1) - bbox%bb(1,1) - half_domain_lengths(1)
     endif

     if(present(zmin_cfd_grid)) then
         zmin_cfd_grid = zg(bbox%ks) - bbox%bb(1,3) - half_domain_lengths(3)
     endif

     if(present(zmax_cfd_grid)) then
         zmax_cfd_grid = zg(bbox%ke) - bbox%bb(1,3) - half_domain_lengths(3)
     endif
 
     if (present(dx_cfd)) then 
         dx_cfd = dx
     endif

     if (present(dz_cfd)) then 
         dz_cfd = dz
     endif

     if (present(cfd_box))then
         cfd_box%gicmin = cfd_box_%gicmin
         cfd_box%gicmax = cfd_box_%gicmax
         cfd_box%gjcmin = cfd_box_%gjcmin
         cfd_box%gjcmax = cfd_box_%gjcmax
         cfd_box%gkcmin = cfd_box_%gkcmin
         cfd_box%gkcmax = cfd_box_%gkcmax
         cfd_box%icmin  = cfd_box_%icmin
         cfd_box%icmax  = cfd_box_%icmax
         cfd_box%jcmin  = cfd_box_%jcmin
         cfd_box%jcmax  = cfd_box_%jcmax
         cfd_box%kcmin  = cfd_box_%kcmin
         cfd_box%kcmax  = cfd_box_%kcmax
         cfd_box%icmino = cfd_box_%icmino
         cfd_box%icmaxo = cfd_box_%icmaxo
         cfd_box%jcmino = cfd_box_%jcmino
         cfd_box%jcmaxo = cfd_box_%jcmaxo
         cfd_box%kcmino = cfd_box_%kcmino
         cfd_box%kcmaxo = cfd_box_%kcmaxo
         cfd_box%xmin  = cfd_box_%xmin
         cfd_box%xmax  = cfd_box_%xmax
         cfd_box%dx    = cfd_box_%dx
         cfd_box%ymin  = cfd_box_%ymin
         cfd_box%ymax  = cfd_box_%ymax
         allocate( cfd_box%y(cfd_box%jcmin:cfd_box%jcmax))
         cfd_box%y(cfd_box%jcmin:cfd_box%jcmax) = cfd_box_%y(cfd_box%jcmin:cfd_box%jcmax)
         cfd_box%zmin = cfd_box_%zmin
         cfd_box%zmax = cfd_box_%zmax
         cfd_box%dz = cfd_box_%dz
     end if

end subroutine coupler_md_get

!-----------------------------------------------------------------------------
 
function coupler_md_get_save_period() result(p)
	use coupler_module, only : save_period
	implicit none

	integer p
	p = save_period

end function coupler_md_get_save_period

!----------------------------------------------------------------------------- 

function coupler_md_get_average_period() result(p)
	use coupler_module, only : average_period
	implicit none

	integer p
	p = average_period

end function coupler_md_get_average_period

!----------------------------------------------------------------------------- 

function coupler_md_get_md_steps_per_cfd_dt() result(p)
	use coupler_input_data, only : md_steps_per_dt_cfd
	implicit none

	integer p
	p = md_steps_per_dt_cfd

end function coupler_md_get_md_steps_per_cfd_dt

!-----------------------------------------------------------------------------

function coupler_md_get_nsteps() result(n)
	use coupler_module, only :  nsteps_md
	implicit none 

	 integer n
	 n = nsteps_md

end function coupler_md_get_nsteps

!-----------------------------------------------------------------------------

function coupler_md_get_dt_cfd() result(t)
	use coupler_module, only : dt_CFD  
	implicit none

	real(kind=kind(0.d0)) t
	t = dt_CFD

end function coupler_md_get_dt_cfd

!------------------------------------------------------------------------------

function coupler_md_get_density() result(r)
	use coupler_module, only : density_md
	implicit none

	real(kind(0.d0)) r
	r = density_md

end function coupler_md_get_density

!------------------------------------------------------------------------------
function coupler_md_get_cfd_id() result(r)
	use coupler_module, only : cfd_code_id
    implicit none

    integer r
    r = cfd_code_id

end function coupler_md_get_cfd_id

end module coupler
