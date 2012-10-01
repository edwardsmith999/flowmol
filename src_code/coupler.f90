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

    interface coupler_send
        module procedure coupler_send_3d, coupler_send_4d
    end interface

    interface coupler_recv
        module procedure coupler_recv_3d, coupler_recv_4d
    end interface

    private coupler_send_3d, coupler_send_4d, &
        coupler_send_xd, coupler_recv_3d, coupler_recv_4d,&
        coupler_recv_xd

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

subroutine coupler_create_comm(callingrealm, REALM_COMM, ierror)
	use mpi
	use coupler_module, only : rank_world,myid_world,rootid_world, & 
							   realm, rank_realm,myid_realm,rootid_realm, ierr, & 
								CPL_WORLD_COMM, CPL_REALM_COMM, CPL_INTER_COMM
	use coupler_input_data
	implicit none

	integer, intent(in) :: callingrealm ! CFD or MD
	integer, intent(out):: REALM_COMM, ierror

	!Get processor id in world across both realms
	call MPI_comm_rank(MPI_COMM_WORLD,myid_world,ierr)
	rank_world = myid_world + 1; rootid_world = 0

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
	REALM_COMM		= MPI_COMM_NULL
	CPL_REALM_COMM 	= MPI_COMM_NULL

	!------------ create realm intra-communicators -----------------------
	! Split MPI_COMM_WORLD into an intra-communicator for each realm 
	! (used for any communication within each realm - e.g. broadcast from 
	!  an md process to all other md processes) 
	! THIS COMMUNICATOR IS REDEFINED TO BE THE CARTESIAN COMMUNICATOR OF 
	! BOTH CODES DURING THE INITIALISATION
	call MPI_comm_split(CPL_WORLD_COMM,callingrealm,myid_world,REALM_COMM,ierr)

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

	write(0,*) 'did (inter)communicators ', code_name(realm), myid_world

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

	!Setup graph topology
	!call CPL_overlap_topology


	!Write Debug information
	if (realm .eq. cfd_realm) then
		call write_map_cfd
	else if (realm .eq. md_realm) then
		call write_map_md
	else
		write(*,*) "Wrong realm in coupler_create_map"
		call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_REALM,ierr)
	end if
	call write_map_olap

	call MPI_barrier(CPL_WORLD_COMM,ierr)

contains

	!------------------------------------------------------------
	!Calculate processor cell ranges of MD code on all processors
		
	subroutine get_md_cell_ranges
	use coupler_input_data, only : cfd_coupler_input
	implicit none

		integer :: startproc
		integer, parameter :: Tnull = -666

		allocate(icPmin_md(npx_md)); icPmin_md = Tnull
		allocate(jcPmin_md(npy_md)); jcPmin_md = Tnull
		allocate(kcPmin_md(npz_md)); kcPmin_md = Tnull
		allocate(icPmax_md(npx_md)); icPmax_md = Tnull
		allocate(jcPmax_md(npy_md)); jcPmax_md = Tnull
		allocate(kcPmax_md(npz_md)); kcPmax_md = Tnull

		! - - x - -
		nlgx_md = ceiling(dble(ncx)/dble(npx_md))
		do n=1,npx_md
			icPmax_md(n) = n * nlgx_md
			icPmin_md(n) = icPmax_md(n) - nlgx_md + 1
		end do	

		! - - y - -
		nlgy_md = ceiling(dble(ncy)/dble(npy_md))
		yL_olap = ncy_olap * dy
		yL_puremd = yL_md - yL_olap
		ncy_puremd = yL_puremd / dy
		yLl_md = yL_md / npy_md
		startproc = ceiling(yL_puremd/yLl_md)
		do n = startproc,npy_md
			jcPmax_md(n) = n * nlgy_md - ncy_puremd
			jcPmin_md(n) = jcPmax_md(n) - nlgy_md + 1
			if (jcPmin_md(n).le.0) jcPmin_md(n) = 1
		end do 

		! - - z - -
		nlgz_md = ceiling(dble(ncz)/dble(npz_md))
		do n=1,npz_md
			kcPmax_md(n) = n * nlgz_md
			kcPmin_md(n) = kcPmax_md(n) - nlgz_md + 1
		end do

	end subroutine get_md_cell_ranges

	!------------------------------------------------------------
	!Calculate processor overlap between CFD/MD on all processors
	
	subroutine get_overlap_blocks
	implicit none

		integer 			:: iblock_realm,jblock_realm,kblock_realm
		integer 			:: i,endproc,nolapsx,nolapsy,nolapsz
		integer,dimension(3):: pcoords
		integer, parameter 	:: olap_null = -666

		!Get cartesian coordinate of overlapping md cells & cfd cells
		allocate(cfd_icoord2olap_md_icoords(npx_cfd,npx_md/npx_cfd)); cfd_icoord2olap_md_icoords = olap_null
		allocate(cfd_jcoord2olap_md_jcoords(npy_cfd,npy_md/npy_cfd)); cfd_jcoord2olap_md_jcoords = olap_null
		allocate(cfd_kcoord2olap_md_kcoords(npz_cfd,npz_md/npz_cfd)); cfd_kcoord2olap_md_kcoords = olap_null

		! - - x - -
		nolapsx = ceiling(dble(npx_md)/dble(npx_cfd))
		do n = 1,npx_cfd
		do i = 1,nolapsx	
			cfd_icoord2olap_md_icoords(n,i) = (n-1)*nolapsx + i
		end do
		end do

		! - - y - -
		yLl_md = yL_md/dble(npy_md)
		nolapsy = ceiling(yL_olap/yLl_md)
		yLl_cfd = yL_cfd / npy_cfd
		endproc = ceiling(yL_olap/yLl_cfd)
		do n = 1,endproc
		do i = 1,nolapsy
			cfd_jcoord2olap_md_jcoords(n,i) = (n-1)*nolapsy + i + (npy_md - nolapsy)
		end do
		end do

		! - - z - -
		nolapsz = ceiling(dble(npz_md)/dble(npz_cfd))
		do n = 1,npz_cfd
		do i = 1,nolapsz	
			cfd_kcoord2olap_md_kcoords(n,i) = (n-1)*nolapsz + i
		end do
		end do

		!Create communicator for overlapping processors
		if (realm .eq. cfd_realm) then
			map%n = npx_md/npx_cfd
			allocate(map%rank_list(map%n))
			!Get rank(s) of overlapping MD processor(s)
			do n = 1,map%n
				pcoords(1)=cfd_icoord2olap_md_icoords(rank2coord_cfd(1,rank_realm),n)
				pcoords(2)=cfd_jcoord2olap_md_jcoords(rank2coord_cfd(2,rank_realm),n)
				pcoords(3)=cfd_kcoord2olap_md_kcoords(rank2coord_cfd(3,rank_realm),n)
				map%rank_list(n) = coord2rank_md(pcoords(1),pcoords(2),pcoords(3))
				map%rank_list(n) = map%rank_list(n) + 1
				write(250+rank_realm,'(2a,6i5)'), 'overlap',code_name(realm),rank_realm,map%n,map%rank_list(n),pcoords
			enddo
		else if (realm .eq. md_realm) then
			map%n = 1
			allocate(map%rank_list(map%n))	
			!Get rank of overlapping CFD processor
			pcoords(1) = rank2coord_md(1,rank_realm)*(dble(npx_cfd)/dble(npx_md))
			pcoords(2) = 1 !rank2coord_md(2,rank_realm)*(dble(npy_cfd)/dble(npy_md))
			pcoords(3) = rank2coord_md(3,rank_realm)*(dble(npz_cfd)/dble(npz_md))
			map%rank_list(1) = coord2rank_cfd(pcoords(1),pcoords(2),pcoords(3))
			map%rank_list(1) = map%rank_list(1) + 1
			write(300+rank_realm,'(2a,6i5)'), 'overlap',code_name(realm),rank_realm,map%n,map%rank_list(1),pcoords
		endif

		call MPI_Barrier(MPI_COMM_WORLD,ierr)
			
	end subroutine get_overlap_blocks

	!Setup topology graph of overlaps between CFD & MD processors

!	subroutine CPL_overlap_topology
!	use coupler_module, only : nproc_olap
!	implicit none

!		integer	:: n, nneighbors
!		logical :: reorder

		!Allow optimisations of ordering
!		reorder = .true.

		!Get number of processors in communicating overlap region 
!		call MPI_comm_size(MPI_OLAP_COMM, nproc_olap, ierr)

		!CFD processor is root and has mapping to all MD processors
!		allocate(index(nproc_olap))			!Index for each processor
!		allocate(edges(2*(nproc_olap-1)))	!nproc_olap-1 for CFD and one for each of nproc_olap MD processors

		!CFD processor has connections to nproc_olap MD processors
!		index(rootid+1) = nproc_olap
!		do n = 1,nproc_olap
!			edges(n) = n
!		enddo

		!MD processor has a single connection to CFD
!		i = 2
!		do n = nproc_olap,2*(nproc_olap-1)
!			index(i) = index(i-1) + 1
!			edges(n) = rootid
!			i = i + 1
!		enddo

		!Create graph topology for overlap region
!		call MPI_Graph_create(MPI_OLAP_COMM, nproc_olap, index, edges,reorder, MPI_OLAP_GRAPH_COMM, ierr)

		! - - TEST - - TEST - -
		!Get number of neighbours
!		call MPI_Graph_neighbors_count( MPI_OLAP_COMM, rank_world, nneighbors, ierr)
		!Get neighbours
!		call MPI_Graph_neighbors( MPI_OLAP_COMM, rank_world, nneighbors, neighbors,ierr )

!	end subroutine CPL_overlap_topology

	!------------------------------------------------------------
	! Write MD cell mapping to file

	subroutine write_map_md
	implicit none

		write(5000+rank_realm,*), '==========================================='
		write(5000+rank_realm,*), '------------  M D   M A P  ----------------'
		write(5000+rank_realm,*), '==========================================='
		write(5000+rank_realm,*), 'npx_md = ', npx_md
		write(5000+rank_realm,*), 'ncx   = ', ncx
		write(5000+rank_realm,*), 'nlgx_md   = ', nlgx_md
		write(5000+rank_realm,*), '-------------------------------------------'
		write(5000+rank_realm,*), '     n        icPmin_md(n)    icPmax_md(n)   '
		write(5000+rank_realm,*), '-------------------------------------------'
		do n=1,npx_md
			write(5000+rank_realm,'(1x,3i11)'), n, icPmin_md(n), icPmax_md(n)
		end do	
		write(5000+rank_realm,*), '-------------------------------------------'
		write(5000+rank_realm,*), 'npy_md = ', npy_md
		write(5000+rank_realm,*), 'nlgy_md = ', nlgy_md 
		write(5000+rank_realm,*), 'yL_pmd = ', yL_puremd
		write(5000+rank_realm,*), 'ncy_puremd = ', ncy_puremd
		write(5000+rank_realm,*), '-------------------------------------------'
		write(5000+rank_realm,*), '     n        jcPmin_md(n)    jcPmax_md(n)   '
		write(5000+rank_realm,*), '-------------------------------------------'
		do n = 1,npy_md	
			write(5000+rank_realm,'(1x,3i11)'), n, jcPmin_md(n), jcPmax_md(n)
		end do
		write(5000+rank_realm,*), '-------------------------------------------'
		write(5000+rank_realm,*), 'npz_md = ', npz_md
		write(5000+rank_realm,*), 'ncz-1   = ', ncz
		write(5000+rank_realm,*), 'nlgz_md   = ', nlgz_md
		write(5000+rank_realm,*), '-------------------------------------------'
		write(5000+rank_realm,*), '     n        kcPmin_md(n)    kcPmax_md(n)   '
		write(5000+rank_realm,*), '-------------------------------------------'
		do n=1,npz_md
			write(5000+rank_realm,'(1x,3i11)'), n, kcPmin_md(n), kcPmax_md(n)
		end do
		write(5000+rank_realm,*), '-------------------------------------------'

	end subroutine write_map_md

	!------------------------------------------------------------
	! Write CFD cell mapping to file

	subroutine write_map_cfd
	implicit none

		write(10000+rank_realm,*), '==========================================='
		write(10000+rank_realm,*), '------------ C F D   M A P ----------------'
		write(10000+rank_realm,*), '==========================================='
		write(10000+rank_realm,*), 'npx_cfd = ', npx_cfd
		write(10000+rank_realm,*), 'ncx   = ', ncx
		write(10000+rank_realm,*), 'nlgx_cfd   = ', nlgx_cfd
		write(10000+rank_realm,*), '-------------------------------------------'
		write(10000+rank_realm,*), '     n        icPmin_cfd(n)    icPmax_cfd(n)   '
		write(10000+rank_realm,*), '-------------------------------------------'
		do n=1,npx_cfd
			write(10000+rank_realm,'(1x,3i11)'), n, icPmin_cfd(n), icPmax_cfd(n)
		end do	
		write(10000+rank_realm,*), '-------------------------------------------'
		write(10000+rank_realm,*), 'npy_cfd = ', npy_cfd
		write(10000+rank_realm,*), 'nlgy_cfd = ', nlgy_cfd 
		write(10000+rank_realm,*), '-------------------------------------------'
		write(10000+rank_realm,*), '     n        jcPmin_cfd(n)    jcPmax_cfd(n)   '
		write(10000+rank_realm,*), '-------------------------------------------'
		do n = 1,npy_cfd	
			write(10000+rank_realm,'(1x,3i11)'), n, jcPmin_cfd(n), jcPmax_cfd(n)
		end do
		write(10000+rank_realm,*), '-------------------------------------------'
		write(10000+rank_realm,*), 'npz_cfd = ', npz_cfd
		write(10000+rank_realm,*), 'ncz   = ', ncz
		write(10000+rank_realm,*), 'nlgz_cfd   = ', nlgz_cfd
		write(10000+rank_realm,*), '-------------------------------------------'
		write(10000+rank_realm,*), '     n        kcPmin_cfd(n)    kcPmax_cfd(n)   '
		write(10000+rank_realm,*), '-------------------------------------------'
		do n=1,npz_cfd
			write(10000+rank_realm,'(1x,3i11)'), n, kcPmin_cfd(n), kcPmax_cfd(n)
		end do
		write(10000+rank_realm,*), '-------------------------------------------'

	end subroutine write_map_cfd

	!------------------------------------------------------------
	! Write overlap block mapping to file

	subroutine write_map_olap
		implicit none

		integer	:: nolapsx,nolapsy,nolapsz

		nolapsx = npx_md/npx_cfd
		nolapsy = ceiling(yL_olap/yLl_md)
		print*, nolapsy,yL_olap,yLl_md
		nolapsz = npz_md/npz_cfd

		write(7500+rank_realm,*), ''
		write(7500+rank_realm,*), '==========================================='
		write(7500+rank_realm,*), '-------- O V E R L A P   M A P ------------'
		write(7500+rank_realm,*), '==========================================='
		write(7500+rank_realm,*), 'nolapsx = ', nolapsx
		write(7500+rank_realm,*), '-------------------------------------------'
		write(7500+rank_realm,*), '  rank2coord_cfd       olapmin     olapmax     ' 
		write(7500+rank_realm,*), '-------------------------------------------'
		do n=1,npx_cfd
			write(7500+rank_realm,'(1x,3i11)'), n,               &
	              cfd_icoord2olap_md_icoords(n,1),         &
			      cfd_icoord2olap_md_icoords(n,nolapsx)
		end do	
		write(7500+rank_realm,*), '-------------------------------------------'

		write(7500+rank_realm,*), 'nolapsy = ', nolapsy
		write(7500+rank_realm,*), '-------------------------------------------'
		write(7500+rank_realm,*), '  jcoord_cfd       olapmin     olapmax     ' 
		write(7500+rank_realm,*), '-------------------------------------------'
		do n=1,npy_cfd
			write(7500+rank_realm,'(1x,3i11)'), n,               &
	              cfd_jcoord2olap_md_jcoords(n,1),         &
			      cfd_jcoord2olap_md_jcoords(n,nolapsy)
		end do	
		write(7500+rank_realm,*), '-------------------------------------------'

		write(7500+rank_realm,*), 'nolapsx = ', nolapsz
		write(7500+rank_realm,*), '-------------------------------------------'
		write(7500+rank_realm,*), '  kcoord_cfd       olapmin     olapmax     ' 
		write(7500+rank_realm,*), '-------------------------------------------'
		do n=1,npz_cfd
			write(7500+rank_realm,'(1x,3i11)'), n,               &
	              cfd_kcoord2olap_md_kcoords(n,1),         &
			      cfd_kcoord2olap_md_kcoords(n,nolapsz)
		end do	
		write(7500+rank_realm,*), '-------------------------------------------'

	end subroutine write_map_olap

end subroutine coupler_create_map

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
	use coupler_internal_cfd, only: bbox_cfd
    use coupler_internal_md, only : bbox_md => bbox
	implicit none
 
    ! array containing data distributed on the grid
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in) :: asend
    ! specify the order of dimensions in asend default is ( x, y, z) 
	! but some CFD solvers use (z,x,y)     
    integer, optional, intent(in) :: index_transpose(3)    ! rule of transposition for the coordinates 

	!Number of halos
	integer	:: nh = 1    

    ! local indices 
    integer :: ix,iy,iz,bicmin,bicmax,bjcmin,bjcmax,bkcmin,bkcmax
	integer	:: ig(2,3),a_components

    ! auxiliaries 
    integer	:: i,ndata,itag,dest
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
	endif

    ! Re-order indices of x,y,z coordinates (handles transpose arrays)
    ix = 1 ; iy = 2; iz = 3
    if( present(index_transpose)) then
        ix = index_transpose(1)
        iy = index_transpose(2)
        iz = index_transpose(3)
    endif

    ! Number of components at each grid point
    a_components = size(asend,dim=1)
   
    ! loop through the maps and send the corresponding sections of asend
    do i = 1, map%n

		!Get taget processor from mapping
        dest = map%rank_list(i)

        ! Amount of data to be sent
        ndata = a_components * (bicmax-bicmin + 1) * (bjcmax-bjcmin + 1) * (bkcmax-bkcmin + 1)
		if(allocated(vbuf)) deallocate(vbuf); allocate(vbuf(ndata))
		vbuf(1:ndata) = reshape(asend, (/ ndata /) )

        ! Send data 
        itag = mod( ncalls, MPI_TAG_UB) !Attention ncall could go over max tag value for long runs!!
		call MPI_send(vbuf, ndata, MPI_DOUBLE_PRECISION, dest, itag, CPL_INTER_COMM, ierr)
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
                                                         
    ! local indices 
    integer ::ix,iy,iz,bicmin,bicmax,bjcmin,bjcmax,bkcmin,bkcmax
	integer	::i,ig(2,3),pcoords(3),recvdata,startbuf,endbuf,a_components

    ! auxiliaries 
    integer ndata,itag,source,req(map%n),vel_indx(8,map%n), &
         start_address(map%n+1), status(MPI_STATUS_SIZE,map%n)
    real(kind(0.d0)),dimension(:), allocatable ::  vbuf,vbuf2
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
	ig(1,ix) = bicmin;	ig(2,ix) = bicmax
	ig(1,iy) = bjcmin;	ig(2,iy) = bjcmax
	ig(1,iz) = bkcmin;	ig(2,iz) = bkcmax

    ! Amount of data to receive
	a_components = size(arecv,dim=1)
    ndata = a_components * (bicmax-bicmin + 1) * (bjcmax-bjcmin + 1) * (bkcmax-bkcmin + 1)
	print*, 'vbuf size', ndata , a_components , (bicmax-bicmin + 1) , (bjcmax-bjcmin + 1) , (bkcmax-bkcmin + 1)
	allocate(vbuf(ndata),stat=ierr) 

    ! Receive from all attached processors
    start_address(1) = 1 
    do i = 1, map%n

		!Get source processor from mapping
        source =  map%rank_list(i)

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
			itag = mod(ncalls, MPI_TAG_UB) ! Attention ncall could go over max tag value for long runs!!
			call MPI_irecv(vbuf2(:),recvdata,MPI_DOUBLE_PRECISION,source,itag,&
                				CPL_INTER_COMM,req(i),ierr)
			startbuf = map%rank_list(rank_realm) * ndata
			endbuf   = map%rank_list(rank_realm) * (ndata + 1 )
			vbuf(1:ndata) = vbuf2(startbuf:endbuf)
			deallocate(vbuf2)
		else
			! Otherwise Receive section of data and increment pointer 
			! ready to receive next piece of data
	        start_address(i+1) = start_address(i) + recvdata
			print*, start_address(i+1),start_address(i),recvdata
			itag = mod(ncalls, MPI_TAG_UB) ! Attention ncall could go over max tag value for long runs!!
			call MPI_irecv(vbuf(start_address(i)),recvdata,MPI_DOUBLE_PRECISION,source,itag,&
                				CPL_INTER_COMM,req(i),ierr)
		endif
    enddo
    call MPI_waitall(map%n, req, status, ierr)

	arecv = reshape(vbuf,(/ ndata, ig(2,1)-ig(1,1)+1, ig(2,2)-ig(1,2)+1, ig(2,3)-ig(1,3)+1 /))
           
end subroutine coupler_recv_xd

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
			!print'(2a,2i8,4f25.16)', 'ICP send data',code_name(realm), rank_realm, & 
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
	!print'(2a,2i8,4f25.16)', 'ICP recv data',code_name(realm),rank_realm, & 
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

	!print'(2a,2i8,4f25.16)', 'ICP trfr data',code_name(realm),rank_realm, & 
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
