!=============================================================================
!
!				  Coupler 
!
!! Routines accessible from application ( molecular or continuum ) after 
!! the name, in parenthesis, is the realm in which each routine must be called
!!
! SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP
!
!! - CPL_create_comm	      (cfd+md)   splits MPI_COMM_WORLD, create inter - 
!!                                   communicator between CFD and MD
!!
!! - CPL_create_map	      (cfd+md)   creates correspondence maps between 
!!                                      the CFD grid and MD domains
!!
!! - CPL_cfd_adjust_domain     (cfd)    adjust CFD tomain to an integer number 
!!                                      FCC or similar MD initial layout
!!
! SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION SIMULATION
!!
!! - CPL_send_data        	  (cfd+md)   sends grid data exchanged between 
!!                                      realms ( generic interface)
!!
!! - CPL_recv_data        	  (cfd+md)   receives data exchanged between realms 
!!                                      ( generic interface)
!!
!! - CPL_cfd_get               (cfd)    returns coupler internal parameters 
!!                                      for CFD realm
!!
!! - CPL_md_get                 (md)    returns coupler internal parameters 
!!                                      for MD realm
!!
!! - CPL_md_get_save_period     (md)    auxiliary used for testing
!!
!! - CPL_md_get_average_period  (md)    returns average period of BC
!!
!! - CPL_md_get_md_per_cfd_dt   (md) 	returns the number of step MD does for 
!!                                      each CFD step
!!
!! - CPL_md_get_nsteps          (md)    returm CFD nsteps  
!!
!! - CPL_md_get_dt_cfd          (md)    returns MD dt
!!
!! - CPL_md_set                 (md)    sets zL if CFD is 2D
!!
!! - CPL_md_get_density         (md)    gets CFD density
!!
!! - CPL_md_get_cfd_id          (md)    id for CFD code, possible values set 
!!                                      in coupler_parameters
!!
!! @author  Lucian Anton, November 2011  
!! @author Edward Smith, Dave Trevelyan September 2012
!! @see coupler_module
!! @see coupler_internal_cfd
!! @see coupler_internal_md
!! @see coupler_parameters
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

contains

!=============================================================================
!					 _____      _               
!					/  ___|    | |              
!					\ `--.  ___| |_ _   _ _ __  
!					 `--. \/ _ \ __| | | | '_ \ 
!					/\__/ /  __/ |_| |_| | |_) |
!					\____/ \___|\__|\__,_| .__/ 
!					                     | |    
!					                     |_|    
!=============================================================================

!=============================================================================
! 							coupler_create_comm	      	
!! (cfd+md) Splits MPI_COMM_WORLD in both the CFD and MD code respectively
!! 		   and create intercommunicator between CFD and MD
!-----------------------------------------------------------------------------

subroutine CPL_create_comm(callingrealm, RETURNED_REALM_COMM, ierror)
	use mpi
	use coupler_module, only : rank_world,myid_world,rootid_world,nproc_world,&
	                           realm, rank_realm,myid_realm,rootid_realm,ierr,& 
	                           CPL_WORLD_COMM, CPL_REALM_COMM, CPL_INTER_COMM
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
	print*, 'did (inter)communicators ', realm_name(realm), myid_world

end subroutine create_comm

end subroutine CPL_create_comm

!=============================================================================
!! Establish for all MD processors the mapping (if any) 
!! to coupled CFD processors
!-----------------------------------------------------------------------------

subroutine CPL_create_map
	use mpi
	use coupler_module
	implicit none

	! Check (as best as one can) that the inputs will work
	call check_config_feasibility

	! Get ranges of cells on each MD processor
	call get_md_cell_ranges

	! Get overlapping mapping for MD to CFD
	call get_overlap_blocks

	! Setup overlap communicators
	call prepare_overlap_comms

	! Setup graph topology
	call CPL_overlap_topology

contains

!------------------------------------------------------------
!Calculate processor cell ranges of MD code on all processors
	
subroutine get_md_cell_ranges
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
	integer             :: xLl_md, yLl_md, zLl_md
	integer				:: yLl_cfd
	integer,dimension(3):: pcoords

	xL_olap = ncx_olap * dx 
	yL_olap = ncy_olap * dy 
	zL_olap = ncz_olap * dz 

	xLl_md  = xL_md / npx_md
	yLl_md  = yL_md / npy_md
	zLl_md  = zL_md / npz_md

	if (realm .eq. md_realm) then
		xLl = xLl_md; yLl = yLl_md ; zLl = zLl_md 
	endif

	nolapsx = nint( dble( npx_md ) / dble( npx_cfd ) )
	nolapsy = ceiling ( yL_olap / yLl_md ) 
	nolapsz = nint( dble( npz_md ) / dble( npz_cfd ) )

	!Get cartesian coordinate of overlapping md cells & cfd cells
	allocate(cfd_icoord2olap_md_icoords(npx_cfd,nolapsx)) 
	allocate(cfd_jcoord2olap_md_jcoords(npy_cfd,nolapsy)) 
	allocate(cfd_kcoord2olap_md_kcoords(npz_cfd,nolapsz)) 
	cfd_icoord2olap_md_icoords = VOID
	cfd_jcoord2olap_md_jcoords = VOID
	cfd_kcoord2olap_md_kcoords = VOID

	! - - x - -
	do n = 1,npx_cfd
	do i = 1,nolapsx	
		cfd_icoord2olap_md_icoords(n,i) = (n-1)*nolapsx + i
	end do
	end do

	! - - y - -
	yLl_cfd = yL_cfd/npy_cfd
	endproc = ceiling(yL_olap/yLl_cfd)
	do n = 1,endproc
	do i = 1,nolapsy
		cfd_jcoord2olap_md_jcoords(n,i) =   (n-1)*nolapsy + i &
										  + (npy_md - nolapsy)
	end do
	end do

	! - - z - -
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

subroutine check_config_feasibility
	implicit none

	integer :: ival
	character(len=256) :: string

	! Check there is only one overlap CFD proc in y
	ival = nint( dble(ncy) / dble(npy_cfd) )
	if (ncy_olap .gt. ival) then

	    string = "This coupler will not work if there is more than one "// &
		         "CFD (y-coordinate) in the overlapping region. "       // &
		         "Aborting simulation."

		call error_abort(string)

	end if


	! Check whether ncx,ncy,ncz are an integer multiple of npx_md, etc.
	! - N.B. no need to check ncy/npy_md.
	ival = 0
	ival = ival + mod(ncx,npx_cfd)
	ival = ival + mod(ncy,npy_cfd)
	ival = ival + mod(ncz,npz_cfd)
	ival = ival + mod(ncx,npx_md)
	ival = ival + mod(ncz,npz_md)
	if (ival.ne.0) then	
		string = "The number of cells in the cfd domain is not an "    // &
		         "integer multiple of the number of processors in "    // &
		         "the x and z directions. Aborting simulation."	
		call error_abort(string)
	end if

	! Check that the MD region is large enough to cover overlap
	ival = 0
	ival = ival + floor(xL_olap/xL_md)
	ival = ival + floor(yL_olap/yL_md)
	ival = ival + floor(zL_olap/zL_md)
	if (ival.ne.0) then	
		string = "Overlap region is larger than the MD region. "       // &
		         "Aborting simulation."
		call error_abort(string)
	end if


	! Check overlap cells are within CFD extents
	ival = 0
	if (icmin_olap.lt.icmin) ival = ival + 1		
	if (icmax_olap.gt.icmax) ival = ival + 1		
	if (jcmin_olap.lt.jcmin) ival = ival + 1		
	if (jcmax_olap.gt.jcmax) ival = ival + 1		
	if (kcmin_olap.lt.kcmin) ival = ival + 1		
	if (kcmax_olap.gt.kcmax) ival = ival + 1		
	if (ival.ne.0) then
		string = "Overlap region has been specified outisde of the "  // &
		         "CFD region. Aborting simulation."
		call error_abort(string)
	end if
	
		
	! Check MD/CFD ratios are integers in x and z
	if (mod(npx_md,npx_cfd) .ne. 0) then

		print'(a,i8,a,i8)', ' number of MD processors in x ', npx_md,     & 
							' number of CFD processors in x ', npx_cfd

		call error_abort("get_overlap_blocks error - number of MD "    // & 
						 "processors in x must be an integer multiple "// &
						 "of number of CFD processors in x")

	elseif (mod(npz_md,npz_cfd) .ne. 0) then

		print'(a,i8,a,i8)', ' number of MD processors in z ', npz_md,     & 
							' number of CFD processors in z ', npz_cfd

		call error_abort("get_overlap_blocks error - number of MD "    // &
						 "processors in z must be an integer multiple "// &
						 "of number of CFD processors in z")

	endif

end subroutine check_config_feasibility

!subroutine intersect_comm
!	use coupler_module
!	use mpi
!	implicit none

!	integer :: n
!	integer,dimension(3)   :: pcoords

	!Create communicator for all intersecting processors
!	if (realm .eq. cfd_realm) then
!		map%n = npx_md/npx_cfd
!		allocate(map%rank_list(map%n))
		!Get rank(s) of overlapping MD processor(s)
!		do n = 1,map%n
!			pcoords(1)=cfd_icoord2olap_md_icoords(rank2coord_cfd(1,rank_realm),n)
!			pcoords(2)=cfd_jcoord2olap_md_jcoords(rank2coord_cfd(2,rank_realm),n)
!			pcoords(3)=cfd_kcoord2olap_md_kcoords(rank2coord_cfd(3,rank_realm),n)
!			if (any(pcoords(:).eq.VOID)) then
!				map%n = 0; map%rank_list(:) = VOID
!			else
!				map%rank_list(n) = coord2rank_md(pcoords(1),pcoords(2),pcoords(3))
!			endif
!			write(250+rank_realm,'(2a,6i5)'), 'overlap',realm_name(realm),rank_realm,map%n,map%rank_list(n),pcoords
!		enddo
!	else if (realm .eq. md_realm) then
!		map%n = 1
!		allocate(map%rank_list(map%n))	
		!Get rank of overlapping CFD processor
!		pcoords(1) = rank2coord_md(1,rank_realm)*(dble(npx_cfd)/dble(npx_md))
!		pcoords(2) = npy_cfd !rank2coord_md(2,rank_realm)*(dble(npy_cfd)/dble(npy_md))
!		pcoords(3) = rank2coord_md(3,rank_realm)*(dble(npz_cfd)/dble(npz_md))
!!		map%rank_list(1) = coord2rank_cfd(pcoords(1),pcoords(2),pcoords(3))
!		write(300+rank_realm,'(2a,6i5)'), 'overlap',realm_name(realm),rank_realm,map%n,map%rank_list(1),pcoords
!	endif

!end subroutine intersect_comm

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
	integer :: tempsize

	tempsize = size(cfd_icoord2olap_md_icoords,2)
	allocate(mdicoords(tempsize))
	tempsize = size(cfd_jcoord2olap_md_jcoords,2)
	allocate(mdjcoords(tempsize))
	tempsize = size(cfd_kcoord2olap_md_kcoords,2)
	allocate(mdkcoords(tempsize))
	
	allocate(olap_mask(nproc_world))

	!Set default values, must be done because coord2rank_md cannot
	!take "null" coordinates.
	group(:) = olap_null
	olap_mask(:) = .false.
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
			olap_mask(trank_world) = .true.
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

			olap_mask(trank_world) = .true.
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
	!if (myid_olap .eq. CFDid_olap) testval = group(rank_world)
	!call MPI_bcast(testval,1,MPI_INTEGER,CFDid_olap,CPL_OLAP_COMM,ierr)

	! Set all non-overlapping processors to MPI_COMM_NULL
	if (olap_mask(rank_world).eq..false.) then
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
	
	!if (realm.eq.md_realm) call write_overlap_comms_md

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
	if (olap_mask(rank_world).eq..true.) then

		!CFD processor is root and has mapping to all MD processors
		allocate(index(nproc_olap))			! Index for each processor
		allocate(edges(2*(nproc_olap)-1))	! nproc_olap-1 for CFD and one for
		                                    ! each of nproc_olap MD processors
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


end subroutine CPL_create_map

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

!=============================================================================
!	Adjust CFD domain size to an integer number of lattice units used by  
!	MD if sizes are given in sigma units
!-----------------------------------------------------------------------------

subroutine CPL_cfd_adjust_domain(xL, yL, zL, nx, ny, nz, density_output)
    use mpi
    use coupler_module, only : density_cfd,CPL_REALM_COMM, rank_realm, ierr
    implicit none

    integer, optional, intent(inout) 			:: nx, ny, nz
    real(kind(0.d0)), optional, intent(inout) 	:: xL,yL,zL
    real(kind(0.d0)), optional, intent(inout)  	:: density_output

    ! Internal variables
    integer										:: ierror, root
	character(1)				   				:: direction
	logical										:: changed

	!Define root processes
	root = 1

    density_output = density_cfd

	! Check CFD domain and MD domain are compatible sizes to allow a
	! stable initial MD lattice structure - resize if possible or
	! stop code and demand a regeneration of grid if vary by more than 0.01
	changed = .false.
    if (present(xL)) then
        call init_length(xL,resize=.true.,direction='x', &
		                 print_warning=changed)
    endif

	! No need to adjust y because we can adjust DY in MD to
	! have an integer number of FCC units. ??????? What

    if (present(zL)) then
        call init_length(zL,resize=.true.,direction='z', &
		                 print_warning=changed)
    endif

	if ( changed ) then
		print*, "Regenerate Grid with corrected sizes as above"
		call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
	endif

    ! check id CFD cell sizes are larger than 2*sigma 
    call test_cfd_cell_sizes

contains

!-----------------------------------------------------------------------------

subroutine init_length(rout,resize,direction,print_warning)
	use coupler_module, only: dx,dy,dz,error_abort
	implicit none
            
    real(kind=kind(0.d0)), intent(inout) :: rout
    logical, intent(in)                  :: resize
	character(*),intent(in)              :: direction
    logical,intent(out)                  :: print_warning

	real(kind(0.d0)) :: dxyz  ! dx, dy or dz
    real(kind(0.d0)) :: rinit ! initial val of rout or rin for print

    print_warning=.false.

	select case (direction)
	case('x','X')
		dxyz = dx
	case('y','Y')
		dxyz = dy
	case('z','Z')
		dxyz = dz
	case default
		call error_abort('Wrong direction specified in init_length')
	end select

	if ( resize ) then

		rinit = rout
		rout = real(nint(rout/dxyz),kind(0.d0))*dxyz
		print_warning = .true.
		print*, direction, 'dxyz = ', dxyz 

	endif

    if (print_warning) then 

        !if (rank_realm .eq. root) then 

            write(*,'(3(a,/),3a,/,2(a,f20.10),/a,/,a)') &
                    "*********************************************************************",	&
                    "WARNING - this is a coupled run which resets CFD domain size         ",	&
                    " to an integer number of MD initial cells:		                      ", 	&
                    "	Domain resized in the ", direction, " direction			          ",	&
                    " inital size =", rinit, " resized ", rout,									&
                    "								                                      ",	& 
                    "*********************************************************************"   
        !endif

		!If resize is insignificant then return flag print_warning as false
		if (abs(rinit-rout) .lt. 0.01) print_warning = .false.

    end if
    
end subroutine init_length

!-----------------------------------------------------------------------------

subroutine test_cfd_cell_sizes
    implicit none

    if (rank_realm .eq. root) then
        if (present(xL) .and. present(nx)) then
            if (xL/nx < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in x direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          xL=",xL,"nx=",nx
            endif
        endif

        if (present(yL) .and. present(ny)) then
            if (yL/ny < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in y direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          yL=",yL,"nx=",ny
            endif
        endif

        if (present(zL) .and. present(nz)) then
            if (zL/nz < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in z direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          zL=",zL,"nx=",nz
            endif
        endif
    end if
        
end subroutine test_cfd_cell_sizes
            
end subroutine CPL_cfd_adjust_domain

!=============================================================================
!				 _____ _                 _       _   _             
!				/  ___(_)               | |     | | (_)            
!				\ `--. _ _ __ ___  _   _| | __ _| |_ _  ___  _ __  
!				 `--. \ | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_ \ 
!				/\__/ / | | | | | | |_| | | (_| | |_| | (_) | | | |
!				\____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
!
!------------------------------------------------------------------------------
!                              CPL_gather                                     -
!------------------------------------------------------------------------------
!! @author David Trevelyan
!
!! Perform gather operation on CPL_OLAP_COMM communicator. The CFD processor
!! is the root process.  
!!
!! - Synopsis
!!
!!  - CPL_gather(gatherarray,npercell,limits,recvarray)
!!
!! - Input
!!
!!  - gatherarray
!!   - Assumed shape array of data to be gathered (double precision)
!!
!!  - limits
!!   - Integer array of length 6, specifying the global cell extents of the
!!     region to be gathered, is the same on all processors.
!!
!!  - npercell
!!   - number of data points per cell to be gathered (integer).
!!     Note: should be the same as size(gatherarray(1)) for MD proc
!!
!! - Input/Output
!!  - recvarray
!!   - The array in which the gathered values are stored
!!
!! - Output Parameters
!!  - NONE
!!
subroutine CPL_gather(gatherarray,npercell,limits,recvarray)!todo better name than recvarray
	use mpi
	use coupler_module
	implicit none

	integer, intent(in) :: npercell
	integer, intent(in) :: limits(6)
	real(kind(0.d0)), dimension(:,:,:,:), intent(in)    :: gatherarray
	real(kind(0.d0)), dimension(:,:,:,:), intent(inout) :: recvarray

	integer :: sendcount
	integer, dimension(:), allocatable :: recvcounts, displs
	real(kind(0.d0)), dimension(:), allocatable :: sendbuf 
	real(kind(0.d0)), dimension(:), allocatable :: recvbuf 
	
	call prepare_gatherv_parameters	
	
	if (realm.eq.md_realm) call pack_sendbuf

	call MPI_gatherv(sendbuf,sendcount,MPI_DOUBLE_PRECISION,recvbuf,    &
	                 recvcounts,displs,MPI_DOUBLE_PRECISION,CFDid_olap, &
	                 CPL_OLAP_COMM,ierr)

	if (realm.eq.cfd_realm) call unpack_recvbuf 

	call deallocate_gather_u
	
contains

	subroutine prepare_gatherv_parameters
		implicit none

		integer :: coord(3),portion(6)
		integer :: ncells,bufsize
		integer :: trank_olap,tid_olap
		
		! Check send limits are inside overlap region
		if (limits(1) .lt. icmin .or. &
		    limits(2) .gt. icmax .or. &
		    limits(3) .lt. jcmin .or. &
		    limits(4) .gt. jcmax .or. &
		    limits(5) .lt. kcmin .or. &
		    limits(6) .lt. kcmax) then
			
			call error_abort("Gather limits are outside global domain. " // &
			                 "Aborting simulation.")
			
		end if 
		

		! Check if CFD processor has tried to "send" anything
		if (myid_olap.eq.CFDid_olap .and. any(shape(gatherarray).ne.0)) then
			call error_abort('CFD proc input to CPL_gather: '          // &
							 'gatherarray has nonzero size. Aborting ' // &
							 'from prepare_gatherv_parameters' )
		end if

		! Allocate send buffer
		if (realm.eq.md_realm) then
			call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
			call CPL_proc_portion(coord,md_realm,limits,portion,ncells)
			bufsize = npercell*ncells
		else
			bufsize = 0
		end if
		allocate(sendbuf(bufsize))

		! Allocate array of sendbuffer sizes and populate it
		allocate(recvcounts(nproc_olap))
		do trank_olap = 1,nproc_olap
			tid_olap = trank_olap - 1
			if (tid_olap .eq. CFDid_olap) then
				recvcounts(trank_olap) = 0 
			else
				call CPL_Cart_coords(CPL_OLAP_COMM,trank_olap,md_realm,3, &
				                     coord,ierr) 
				call CPL_proc_portion(coord,md_realm,limits,portion,ncells)
				recvcounts(trank_olap) = npercell*ncells 
			end if
		end do
	
		! Own sendbuffer size
		sendcount = size(sendbuf)
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
		integer :: i,j,k
		integer :: coord(3),portion(6),mdextents(6)

		! Get MD processor extents and cells portion of send region
		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
		call CPL_proc_extents(coord,realm,mdextents)
		call CPL_proc_portion(coord,realm,limits,portion)

		! If MD proc has nothing to send, exit
		if (any(portion.eq.VOID)) return

		pos = 1
		do ixyz  = 1,npercell 
		do icell = portion(1),portion(2) 
		do jcell = portion(3),portion(4) 
		do kcell = portion(5),portion(6) 
			! Map to local coords (assumed shape array has
			! lower bound 1 by default when input to subroutine) 
			i = icell - mdextents(1) + 1
			j = jcell - mdextents(3) + 1
			k = kcell - mdextents(5) + 1
			sendbuf(pos) = gatherarray(ixyz,i,j,k)
			pos = pos + 1
		end do
		end do
		end do
		end do
	
	end subroutine pack_sendbuf
	
	subroutine unpack_recvbuf
		implicit none

		integer :: coord(3),portion(6),cfdextents(6)
		integer :: trank_olap, tid_olap
		integer :: pos,ixyz,icell,jcell,kcell
		integer :: i,j,k

		integer :: tempextents(6)

		! Get CFD proc coords and extents, allocate suitable array
		call CPL_cart_coords(CPL_OLAP_COMM,rank_olap,cfd_realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,cfdextents)
		call CPL_proc_portion(coord,cfd_realm,limits,portion)

		! Loop over all processors in overlap comm
		do trank_olap = 1,nproc_olap

			tid_olap = trank_olap - 1
			if (tid_olap .eq. CFDid_olap) cycle

			call CPL_Cart_coords(CPL_OLAP_COMM,trank_olap,md_realm,3,coord,ierr)
			call CPL_proc_portion(coord,md_realm,limits,portion)

			call CPL_proc_extents(coord,md_realm,tempextents)
			
			if (any(portion.eq.VOID)) cycle

			! Set position and unpack MD proc's part of recvbuf to
			! correct region of recvarray	
			pos = displs(trank_olap) + 1	
			do ixyz = 1,npercell
			do icell = portion(1),portion(2)
			do jcell = portion(3),portion(4)
			do kcell = portion(5),portion(6)

				i = icell - cfdextents(1) + 1 
				j = jcell - cfdextents(3) + 1 
				k = kcell - cfdextents(5) + 1 

				recvarray(ixyz,i,j,k) = recvbuf(pos)
				pos = pos + 1

!				write(8000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
!				      'recvarray(',ixyz,',',icell,',',jcell,',',kcell,') =', &
!				       recvarray(ixyz,i,j,k)

			end do	
			end do	
			end do
			end do
					
		end do

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
!                              CPL_scatter                                    -
!------------------------------------------------------------------------------
!! @author David Trevelyan
!!
!! Scatter cell-wise data from CFD processor to corresponding MD processors
!! on the overlap communicator CPL_OLAP_COMM.
!!
!! - Synopsis
!!
!!  - CPL_scatter(scatterarray,npercell,limits,recvarray)
!!
!! - Input
!!
!!  - scatterarray
!!   - assumed shape array of data to be scattered (double precision)
!!
!!  - limits
!!   - integer array of length 6, specifying the global cell extents of the
!!     region to be scattered, is the same on all processors.
!!
!!  - npercell
!!   - number of data points per cell to be scattered (integer).
!!     Note: should be the same as size(scatterarray(1)) for CFD proc
!!
!! - Input/Output
!!  - recvarray
!!   - the array in which the scattered values are stored
!!
!! - Output
!!  - NONE
!! 
subroutine CPL_scatter(scatterarray,npercell,limits,recvarray)
	use coupler_module
	use mpi
	implicit none

	integer,intent(in) :: npercell
	integer,intent(in) :: limits(6)
	real(kind(0.d0)),dimension(:,:,:,:),intent(in)    :: scatterarray
	real(kind(0.d0)),dimension(:,:,:,:),intent(inout) :: recvarray 

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

		integer :: ncells
		integer :: coord(3),portion(6)
		integer :: bufsize
		integer :: trank_olap, tid_olap

		! Check send limits are inside overlap region
		if (limits(1) .lt. icmin .or. &
		    limits(2) .gt. icmax .or. &
		    limits(3) .lt. jcmin .or. &
		    limits(4) .gt. jcmax .or. &
		    limits(5) .lt. kcmin .or. &
		    limits(6) .lt. kcmax) then
			
			call error_abort("Scatter limits are outside global domain. " // &
			                 "Aborting simulation.")
			
		end if 

		if (realm.eq.cfd_realm) then
			call CPL_Cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
			call CPL_proc_portion(coord,cfd_realm,limits,portion,ncells)
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
				call CPL_proc_portion(coord,md_realm,limits,portion,ncells)
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
		integer :: coord(3),cfdextents(6),portion(6)
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
			call CPL_proc_portion(coord,md_realm,limits,portion)
			if (any(portion.eq.VOID)) cycle

			do ixyz = 1,npercell
			do icell= portion(1),portion(2)
			do jcell= portion(3),portion(4)
			do kcell= portion(5),portion(6)
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
		integer :: coord(3),portion(6),extents(6)
		integer :: ixyz,icell,jcell,kcell
		integer :: i,j,k

		call CPL_cart_coords(CPL_OLAP_COMM,rank_olap,md_realm,3,coord,ierr)
		call CPL_proc_extents(coord,realm,extents)
		call CPL_proc_portion(coord,realm,limits,portion)
		if (any(portion.eq.VOID)) return

		pos = 1
		do ixyz = 1,npercell
		do icell= portion(1),portion(2)
		do jcell= portion(3),portion(4)
		do kcell= portion(5),portion(6)
			i = icell - extents(1) + 1
			j = jcell - extents(3) + 1
			k = kcell - extents(5) + 1
			recvarray(ixyz,i,j,k) = recvbuf(pos)
!			write(7000+myid_realm,'(i4,a,i4,a,i4,a,i4,a,i4,a,f20.1)'),        &
!				  rank_cart,' recvarray(',ixyz,',',icell,',',jcell,',',kcell, &
!				  ') =',recvarray(ixyz,i,j,k)
			pos = pos + 1
		end do	
		end do	
		end do
		end do

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
!! CPL_send_data wrapper for 3d arrays
!! see CPL_send_xd for input description
!! @see coupler#subroutine_CPL_send_xd
!-----------------------------------------------------------------------------
subroutine CPL_send_3d(temp,icmin_send,icmax_send,jcmin_send, & 
							jcmax_send,kcmin_send,kcmax_send,send_flag)
	use coupler_module, only : icmin_olap,icmax_olap, & 
							   jcmin_olap,jcmax_olap, &
							   kcmin_olap,kcmax_olap,error_abort
    implicit none
 
	logical, intent(out), optional						:: send_flag
 	integer, intent(in), optional						:: icmax_send,icmin_send
 	integer, intent(in), optional						:: jcmax_send,jcmin_send
 	integer, intent(in), optional						:: kcmax_send,kcmin_send
    real(kind=kind(0.d0)),dimension(:,:,:), intent(in)  :: temp
    
    integer :: n1,n2,n3,n4
	integer	:: icmin,icmax,jcmin,jcmax,kcmin,kcmax
    real(kind=kind(0.d0)),dimension(:,:,:,:),allocatable :: asend


	!Revert to default i domain sending - top of overlap to bottom of overlap
	if ((present(icmax_send)) .and. (present(icmin_send))) then
			icmax = icmax_send;	icmin = icmin_send
	elseif ((.not. present(icmax_send)).and.(.not. present(icmin_send))) then
			icmax = icmax_olap;	icmin = icmin_olap
	else
		call error_abort("CPL_send error - both max and min i limits " // &
						 "required and only one supplied")
	endif

	!Revert to default j domain sending - top of overlap to bottom of overlap
	if ((present(jcmax_send)) .and. (present(jcmin_send))) then
			jcmax = jcmax_send;	jcmin = jcmin_send
	elseif ((.not. present(jcmax_send)).and.(.not. present(jcmin_send))) then
			jcmax = jcmax_olap;	jcmin = jcmin_olap
	else
		call error_abort("CPL_send error - both max and min j limits " // &
						 "required and only one supplied")
	endif

	!Revert to default k domain sending - top of overlap to bottom of overlap
	if ((present(kcmax_send)) .and. (present(kcmin_send))) then
			kcmax = kcmax_send; kcmin = kcmin_send
	elseif ((.not. present(kcmax_send)).and.(.not. present(kcmin_send))) then
			kcmax = kcmax_olap;	kcmin = kcmin_olap
	else
		call error_abort("CPL_send error - both max and min k limits " // &
						 "required and only one supplied")
	endif

   
    n1 = 1
    n2 = size(temp,1)
    n3 = size(temp,2)
    n4 = size(temp,3)

	!Add padding column to 3D array to make it 4D
	allocate(asend(n1,n2,n3,n4))
	asend(1,:,:,:) = temp(:,:,:)

    call CPL_send_xd(asend,icmax,icmin,jcmax, & 
						   jcmin,kcmax,kcmin,send_flag )

end subroutine CPL_send_3d

!=============================================================================
!! CPL_send_data wrapper for 4d arrays
!! see CPL_send_xd for input description
!! @see coupler#subroutine_CPL_send_xd
!-----------------------------------------------------------------------------
subroutine CPL_send_4d(asend,icmin_send,icmax_send,jcmin_send, & 
							 jcmax_send,kcmin_send,kcmax_send,send_flag)
	use coupler_module, only : icmin_olap,icmax_olap, & 
							   jcmin_olap,jcmax_olap, &
							   kcmin_olap,kcmax_olap,error_abort
    implicit none
 
	logical, intent(out), optional						:: send_flag
 	integer, intent(in), optional						:: icmax_send,icmin_send
 	integer, intent(in), optional						:: jcmax_send,jcmin_send
 	integer, intent(in), optional						:: kcmax_send,kcmin_send
	real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in) :: asend
    
    integer	:: npercell,send_extents(4)
	integer	:: icmin,icmax,jcmin,jcmax,kcmin,kcmax

	npercell = size(asend,1)
	!if ((present(icmax_send))) print*, 'icmax_send', icmax_send
	!if ((present(icmin_send))) print*, 'icmin_send', icmin_send
	!if ((present(jcmax_send))) print*, 'jcmax_send', jcmax_send
	!if ((present(jcmin_send))) print*, 'jcmin_send', jcmin_send
	!if ((present(kcmax_send))) print*, 'kcmax_send', kcmax_send
	!if ((present(kcmin_send))) print*, 'kcmin_send', kcmin_send

	!Revert to default i domain sending - top of overlap to bottom of overlap
	if ((present(icmax_send)) .and. (present(icmin_send))) then
			icmax = icmax_send;	icmin = icmin_send
	elseif ((.not. present(icmax_send)).and.(.not. present(icmin_send))) then
			icmax = icmax_olap;	icmin = icmin_olap
	else
		call error_abort("CPL_send error - both max and min i limits " // &
						 "required and only one supplied")
	endif

	!Revert to default j domain sending - top of overlap to bottom of overlap
	if ((present(jcmax_send)) .and. (present(jcmin_send))) then
			jcmax = jcmax_send;	jcmin = jcmin_send
	elseif ((.not. present(jcmax_send)).and.(.not. present(jcmin_send))) then
			jcmax = jcmax_olap;	jcmin = jcmin_olap
	else
		call error_abort("CPL_send error - both max and min j limits " // &
						 "required and only one supplied")
	endif

	!Revert to default k domain sending - top of overlap to bottom of overlap
	if ((present(kcmax_send)) .and. (present(kcmin_send))) then
			kcmax = kcmax_send; kcmin = kcmin_send
	elseif ((.not. present(kcmax_send)).and.(.not. present(kcmin_send))) then
			kcmax = kcmax_olap;	kcmin = kcmin_olap
	else
		call error_abort("CPL_send error - both max and min k limits " // &
						 "required and only one supplied")
	endif

	!send_extents = (/npercell,icmax-icmin+1,jcmax-jcmin+1,kcmax-kcmin+1 /) 
	!if (any(shape(asend) .lt. send_extents)) then
	!	print'(2(a,4i5))', '  Shape of input array = ', shape(asend), & 
	!					  '  Passed range = ',send_extents 
	!	call error_abort("CPL_send error - Specified send range is greater" // &
	!					 "than the number of cells on processor")
	!endif
 
    call CPL_send_xd(asend,icmin,icmax,jcmin, & 
						   jcmax,kcmin,kcmax,send_flag )

end subroutine CPL_send_4d

!=============================================================================
!						CPL_send_xd
!
!! Send data from the local grid to the associated ranks from the other 
!! realm
!!
!! - Synopsis
!!
!!  - CPL_send_xd(asend,icmin_send,icmax_send,jcmin_send,  
!!						     jcmax_send,kcmin_send,kcmax_send,send_flag)
!!
!! - Input Parameters
!!
!!   - asend
!!
!!   - jcmax_recv
!!
!!   - jcmin_recv
!!
!! - Output Parameter
!!
!!   - send_flag
!!
!! @author Edward Smith
! ----------------------------------------------------------------------------
subroutine CPL_send_xd(asend,icmin_send,icmax_send,jcmin_send, & 
						     jcmax_send,kcmin_send,kcmax_send,send_flag)
	use mpi
	use coupler_module, only : CPL_CART_COMM,rank_cart,md_realm,cfd_realm, & 
	                           error_abort,CPL_GRAPH_COMM,myid_graph,olap_mask, &
							   rank_world,rank_realm,realm,rank_olap, & 
							   iblock_realm,jblock_realm,kblock_realm,ierr
	implicit none

	!Flag set if processor has passed data
	logical, intent(out), optional	:: send_flag

    ! Minimum and maximum values to send
 	integer, intent(in)	:: icmax_send,icmin_send, & 
						   jcmax_send,jcmin_send, & 
						   kcmax_send,kcmin_send

   ! Array containing data distributed on the grid
    real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in):: asend
   
	!Number of halos
	integer	:: nh = 1    

	!Neighbours
	integer								:: nneighbors   
	integer,dimension(:),allocatable	:: id_neighbors

    ! local indices 
	integer	:: icell,jcell,kcell
    integer :: ix,iy,iz,icmin,icmax,jcmin,jcmax,kcmin,kcmax
	integer :: n,pos,iclmin,iclmax,jclmin,jclmax,kclmin,kclmax
	integer	:: npercell

    ! auxiliaries 
    integer								:: nbr,ndata,itag,destid,ncells
    integer                             :: pcoords(3)
	integer,dimension(6)				:: extents,limits,portion
	real(kind=kind(0.d0)), allocatable 	:: vbuf(:)

	! This local CFD domain is outside MD overlap zone 
	if (olap_mask(rank_world) .eq. .false.) return

	! Save limits array of Minimum and maximum values to send
	limits = (/ icmin_send,icmax_send,jcmin_send,jcmax_send,kcmin_send,kcmax_send /)

    ! Get local grid box ranges seen by this rank for CFD
    if (realm .eq. cfd_realm) then 
		!Load extents of CFD processor
		pcoords = (/iblock_realm,jblock_realm,kblock_realm /)
		call CPL_proc_extents(pcoords,cfd_realm,extents)
	endif

    ! Number of components at each grid point
    npercell = size(asend,1)

	!Get neighbours
	call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
	allocate(id_neighbors(nneighbors))
	call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr )

	!Set sendflag to false and only change if anything is sent
	if (present(send_flag)) send_flag = .false.
  
    ! loop through the maps and send the corresponding sections of asend
    do nbr = 1, nneighbors

		!Get taget processor from mapping
        destid = id_neighbors(nbr)

		! ----------------- pack data for destid-----------------------------
		if (realm .eq. cfd_realm) then
			!Get extents of nbr MD processor to send to
			call CPL_Cart_coords(CPL_GRAPH_COMM, destid+1,  md_realm, 3, pcoords, ierr) 
		elseif (realm .eq. md_realm) then
			!Data to send is based on current processor
			pcoords = (/iblock_realm,jblock_realm,kblock_realm /)
			!Get extents of current processor
			call CPL_proc_extents(pcoords,md_realm,extents)
		endif

		! If limits passed to send routine, use these instead
		! of overlap/processor limits
		call CPL_proc_portion(pcoords,md_realm,limits,portion,ncells)

        ! Amount of data to be sent
		if (any(portion.eq.VOID)) then
			!print*, 'VOID send qqqq',realm_name(realm),rank_world,rank_realm
			ndata = 0
		else

			! Get data range on processor's local extents
			iclmin = portion(1)-extents(1)+1;	iclmax = portion(2)-extents(1)+1
			jclmin = portion(3)-extents(3)+1;	jclmax = portion(4)-extents(3)+1
			kclmin = portion(5)-extents(5)+1;	kclmax = portion(6)-extents(3)+1

	        ndata = npercell * ncells
			if (allocated(vbuf)) deallocate(vbuf); allocate(vbuf(ndata))
			if (present(send_flag)) send_flag = .true.
			! Pack array into buffer
			pos = 1
			!print'(a,5i4,2i6,i4,24i4)', 'send qqqq',rank_world,rank_realm,rank_olap,ndata,nbr,destid, & 
			!						size(asend),pos,&
			!						iclmin,   iclmax,   jclmin,   jclmax,   kclmin,   kclmax,     &
			!						icmin_send,icmax_send,jcmin_send,jcmax_send,kcmin_send,kcmax_send, & 
			!						portion, extents
			do kcell=kclmin,kclmax
			do jcell=jclmin,jclmax
			do icell=iclmin,iclmax
			do n = 1,npercell
				vbuf(pos) = asend(n,icell,jcell,kcell)
				!write(98000+destid+1+10*rank_world,'(3i8,f20.5)') rank_world,destid+1,n, vbuf(pos)
				pos = pos + 1
			end do
			end do
			end do
			end do
			! ----------------- pack data for destid -----------------------------

 			! Send data 
	        itag = 0 !mod( ncalls, MPI_TAG_UB) !Attention ncall could go over max tag value for long runs!!
			call MPI_send(vbuf, ndata, MPI_DOUBLE_PRECISION, destid, itag, CPL_GRAPH_COMM, ierr)

		endif

    enddo

end subroutine CPL_send_xd

!=============================================================================
!! CPL_recv_xd wrapper for 3d arrays
!! see CPL_recv_xd for input description
!! @see coupler#subroutine_CPL_recv_xd
!-----------------------------------------------------------------------------
subroutine CPL_recv_3d(temp,icmin_recv,icmax_recv,jcmin_recv, & 
							jcmax_recv,kcmin_recv,kcmax_recv,recv_flag)
	use coupler_module, only : icmin_olap,icmax_olap, & 
							   jcmin_olap,jcmax_olap, &
							   kcmin_olap,kcmax_olap,error_abort
    implicit none

	logical, intent(out), optional					:: recv_flag
  	integer, intent(in), optional					:: icmax_recv,icmin_recv
 	integer, intent(in), optional					:: jcmax_recv,jcmin_recv
 	integer, intent(in), optional					:: kcmax_recv,kcmin_recv
    real(kind(0.d0)),dimension(:,:,:),intent(inout) :: temp 
                                                          
    integer	:: n1,n2,n3,n4
	integer	:: icmin,icmax,jcmin,jcmax,kcmin,kcmax
	real(kind(0.d0)),dimension(:,:,:,:),allocatable	 :: arecv

	!if ((present(icmax_recv))) print*, 'icmax_recv', icmax_recv
	!if ((present(icmin_recv))) print*, 'icmin_recv', icmin_recv
	!if ((present(jcmax_recv))) print*, 'jcmax_recv', jcmax_recv
	!if ((present(jcmin_recv))) print*, 'jcmin_recv', jcmin_recv
	!if ((present(kcmax_recv))) print*, 'kcmax_recv', kcmax_recv
	!if ((present(kcmin_recv))) print*, 'kcmin_recv', kcmin_recv

	!Revert to default i domain sending - top of overlap to bottom of overlap
	if ((present(icmax_recv)) .and. (present(icmin_recv))) then
			icmax = icmax_recv;	icmin = icmin_recv
	elseif ((.not. present(icmax_recv)).and.(.not. present(icmin_recv))) then
			icmax = icmax_olap;	icmin = icmin_olap
	else
		call error_abort("CPL_recv error - both max and min i limits " // &
						 "required and only one supplied")
	endif

	!Revert to default j domain sending - top of overlap to bottom of overlap
	if ((present(jcmax_recv)) .and. (present(jcmin_recv))) then
			jcmax = jcmax_recv;	jcmin = jcmin_recv
	elseif ((.not. present(jcmax_recv)).and.(.not. present(jcmin_recv))) then
			jcmax = jcmax_olap;	jcmin = jcmin_olap
	else
		call error_abort("CPL_recv error - both max and min j limits " // &
						 "required and only one supplied")
	endif

	!Revert to default k domain sending - top of overlap to bottom of overlap
	if ((present(kcmax_recv)) .and. (present(kcmin_recv))) then
			kcmax = kcmax_recv; kcmin = kcmin_recv
	elseif ((.not. present(kcmax_recv)).and.(.not. present(kcmin_recv))) then
			kcmax = kcmax_olap;	kcmin = kcmin_olap
	else
		call error_abort("CPL_recv error - both max and min k limits " // &
						 "required and only one supplied")
	endif
 
    n1 = 1 
    n2 = size(temp,1)
    n3 = size(temp,2)
    n4 = size(temp,3)

	!Add padding column to 3D array to make it 4D
	allocate(arecv(n1,n2,n3,n4))
    call CPL_recv_xd(arecv,icmin,icmax,jcmin, & 
						   jcmax,kcmin,kcmax,recv_flag)
	temp(:,:,:) = 	arecv(1,:,:,:) 

end subroutine CPL_recv_3d

!=============================================================================
!! CPL_recv_xd  wrapper for 4d arrays
!! See CPL_recv_xd for input description
!! @see coupler#subroutine_CPL_recv_xd
!-----------------------------------------------------------------------------
subroutine CPL_recv_4d(arecv,icmin_recv,icmax_recv,jcmin_recv, & 
							 jcmax_recv,kcmin_recv,kcmax_recv,recv_flag)
	use coupler_module, only : icmin_olap,icmax_olap, & 
							   jcmin_olap,jcmax_olap, &
							   kcmin_olap,kcmax_olap,error_abort
    implicit none
 
	logical, intent(out), optional						  :: recv_flag
 	integer, intent(in), optional						  :: icmax_recv,icmin_recv
 	integer, intent(in), optional						  :: jcmax_recv,jcmin_recv
 	integer, intent(in), optional						  :: kcmax_recv,kcmin_recv
	real(kind=kind(0.d0)),dimension(:,:,:,:), intent(out) :: arecv
    
    integer	:: n1,n2,n3,n4
	integer	:: icmin,icmax,jcmin,jcmax,kcmin,kcmax

	!if ((present(icmax_recv))) print*, 'icmax_recv', icmax_recv
	!if ((present(icmin_recv))) print*, 'icmin_recv', icmin_recv
	!if ((present(jcmax_recv))) print*, 'jcmax_recv', jcmax_recv
	!if ((present(jcmin_recv))) print*, 'jcmin_recv', jcmin_recv
	!if ((present(kcmax_recv))) print*, 'kcmax_recv', kcmax_recv
	!if ((present(kcmin_recv))) print*, 'kcmin_recv', kcmin_recv

	!Revert to default i domain sending - top of overlap to bottom of overlap
	if ((present(icmax_recv)) .and. (present(icmin_recv))) then
			icmax = icmax_recv;	icmin = icmin_recv
	elseif ((.not. present(icmax_recv)).and.(.not. present(icmin_recv))) then
			icmax = icmax_olap;	icmin = icmin_olap
	else
		call error_abort("CPL_recv error - both max and min i limits " // &
						 "required and only one supplied")
	endif

	!Revert to default j domain sending - top of overlap to bottom of overlap
	if ((present(jcmax_recv)) .and. (present(jcmin_recv))) then
			jcmax = jcmax_recv;	jcmin = jcmin_recv
	elseif ((.not. present(jcmax_recv)).and.(.not. present(jcmin_recv))) then
			jcmax = jcmax_olap;	jcmin = jcmin_olap
	else
		call error_abort("CPL_recv error - both max and min j limits " // &
						 "required and only one supplied")
	endif

	!Revert to default k domain sending - top of overlap to bottom of overlap
	if ((present(kcmax_recv)) .and. (present(kcmin_recv))) then
			kcmax = kcmax_recv; kcmin = kcmin_recv
	elseif ((.not. present(kcmax_recv)).and.(.not. present(kcmin_recv))) then
			kcmax = kcmax_olap;	kcmin = kcmin_olap
	else
		call error_abort("CPL_recv error - both max and min k limits " // &
						 "required and only one supplied")
	endif
 
    call CPL_recv_xd(arecv,icmin,icmax,jcmin, & 
						   jcmax,kcmin,kcmax,recv_flag)

end subroutine CPL_recv_4d


!=============================================================================
!						CPL_recv_xd
!
!! Receive data from to local grid from the associated ranks from the other 
!! realm
!!
!! - Synopsis
!!
!!  - CPL_recv_xdarecv,icmin_recv,icmax_recv,jcmin_recv,  
!!						     jcmax_recv,kcmin_recv,kcmax_recv,recv_flag)
!!
!! - Input Parameters
!!
!!   - arecv
!!
!!   - jcmax_recv
!!
!!   - jcmin_recv
!!
!!   - index_transpose
!!
!! - Output Parameter
!!
!!   - recv_flag
!!
!! @author Edward Smith
! ----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine CPL_recv_xd(arecv,icmin_recv,icmax_recv,jcmin_recv, & 
						     jcmax_recv,kcmin_recv,kcmax_recv,recv_flag)
	use mpi
	use coupler_module, only : CPL_CART_COMM,rank_cart,md_realm,cfd_realm, & 
							   rank_graph, rank_graph2rank_world, &
	                           error_abort,CPL_GRAPH_COMM,myid_graph,olap_mask, &
							   rank_world,rank_realm,realm,rank_olap, & 
							   iblock_realm,jblock_realm,kblock_realm,ierr
	implicit none

	!Flag set if processor has received data
	logical, intent(out), optional					:: recv_flag

    ! Minimum and maximum values of j to send
 	integer, intent(in)	:: icmax_recv,icmin_recv, & 
						   jcmax_recv,jcmin_recv, & 
						   kcmax_recv,kcmin_recv

    ! Array that recieves grid distributed data 
    real(kind(0.d0)), dimension(:,:,:,:),intent(inout):: arecv     

	!Neighbours
	integer								:: nneighbors   
	integer,dimension(:),allocatable	:: id_neighbors
                                                         
    ! local indices 
    integer :: n,nbr,i,j,k,icell,jcell,kcell,icmin,icmax,jcmin,jcmax,kcmin,kcmax
	integer :: pos,iclmin,iclmax,jclmin,jclmax,kclmin,kclmax
	integer	:: pcoords(3),npercell,ndata,ncells
	integer,dimension(6) :: extents,mdportion,portion,limits

    ! auxiliaries 
    integer	:: itag, sourceid,source_realm,start_address
	integer,dimension(:),allocatable   :: req
	integer,dimension(:,:),allocatable :: status
    real(kind(0.d0)),dimension(:), allocatable ::  vbuf
 
	! This local CFD domain is outside MD overlap zone 
	if (olap_mask(rank_world).eq. .false.) return

	! Save limits array of Minimum and maximum values to recv
	limits = (/ icmin_recv,icmax_recv,jcmin_recv,jcmax_recv,kcmin_recv,kcmax_recv /)

    ! Number of components at each grid point
	npercell = size(arecv,1)

    ! Get local grid box ranges seen by this rank for CFD and allocate buffer
    if (realm .eq. cfd_realm) then 

		!Load CFD cells per processor
		call CPL_Cart_coords(CPL_GRAPH_COMM, rank_graph, cfd_realm, 3, pcoords, ierr) 
		call CPL_proc_extents(pcoords,cfd_realm,extents)

		! If limits passed to recv routine, use these instead
		! of overlap/processor limits
		call CPL_proc_portion(pcoords,cfd_realm,limits,portion,ncells)

		! Amount of data to receive from all MD processors
		ndata = npercell * ncells
		allocate(vbuf(ndata)); vbuf = 0.d0

	endif

	!Get neighbours
	call MPI_Graph_neighbors_count(CPL_GRAPH_COMM,myid_graph,nneighbors,ierr)
	allocate(id_neighbors(nneighbors))
	call MPI_Graph_neighbors(CPL_GRAPH_COMM,myid_graph,nneighbors,id_neighbors,ierr )

    ! Receive from all attached processors
	allocate(req(nneighbors))
	allocate(status(MPI_STATUS_SIZE,nneighbors))
    start_address = 1 
    do nbr = 1, nneighbors

		!Get source processor from topology graph
        sourceid =  id_neighbors(nbr)

	    ! Get size of data to receive from source processors
		if (realm .eq. cfd_realm) then
			!CFD realm receives data based on size of MD processor domain
			call CPL_Cart_coords(CPL_GRAPH_COMM, sourceid+1, md_realm, 3, pcoords, ierr) 
		elseif (realm .eq. md_realm) then
			!MD realm receives data as big as own processor domain
			pcoords = (/iblock_realm,jblock_realm,kblock_realm /)
		endif

		! If limits passed to recv routine, use these instead
		! of overlap/processor limits
		call CPL_proc_portion(pcoords,md_realm,limits,portion,ncells)

		!Only receive if overlapping
		if (any(portion.eq.VOID)) then
			ndata = 0
			req(nbr) = MPI_REQUEST_NULL
			if (present(recv_flag)) recv_flag = .false.
		else
			! Amount of data to receive
	        ndata = npercell * ncells
			if (present(recv_flag)) recv_flag = .true.
			
			if (realm .eq. md_realm) then
				!Allocate receive buffer for MD processors
				if (allocated(vbuf)) deallocate(vbuf); allocate(vbuf(ndata)); vbuf=0.d0
			endif

			! Receive section of data
			itag = 0
			call MPI_irecv(vbuf(start_address),ndata,MPI_DOUBLE_PRECISION,sourceid,itag,&
            						CPL_GRAPH_COMM,req(nbr),ierr)

		endif

		!Increment pointer ready to receive next piece of data		
		start_address = start_address + ndata

    enddo
    call MPI_waitall(nneighbors, req, status, ierr)

	!if (rank_world .eq. 33) then
	!	do n = 1,size(vbuf)
	!		write(98000+rank_world,*) rank_world,n, vbuf(n)
	!	enddo
	!endif

	! ----------------- Unpack data -----------------------------
	start_address = 1
	do nbr = 1, nneighbors

		!Get source processor from topology graph
        sourceid =  id_neighbors(nbr)

	    ! Get size of data to receive from source processors
		if (realm .eq. cfd_realm) then
			!CFD always receives data
			if (present(recv_flag)) recv_flag = .true.
			!CFD realm receives data based on size of MD processor domain
			call CPL_Cart_coords(CPL_GRAPH_COMM, sourceid+1, md_realm, 3, pcoords, ierr) 
		elseif (realm .eq. md_realm) then
			!MD realm receives data as big as own processor domain
			pcoords = (/iblock_realm,jblock_realm,kblock_realm /)
			!Get extents of current processor/overlap region
			call CPL_proc_extents(pcoords,md_realm,extents)
		endif

		! If limits passed to recv routine, use these instead
		! of overlap/processor limits
		call CPL_proc_portion(pcoords,md_realm,limits,portion,ncells)
				
		! Unpack array into buffer
		if (any(portion.eq.VOID)) then
			!print*, 'VOID recv qqqq',realm_name(realm),rank_world,rank_realm,rank_graph2rank_world(sourceid+1),recv_flag
			ndata = 0
		else
			! Get local extents in received region
			pos = start_address; ndata = npercell * ncells
			iclmin = portion(1)-extents(1)+1;	iclmax = portion(2)-extents(1)+1
			jclmin = portion(3)-extents(3)+1;	jclmax = portion(4)-extents(3)+1
			kclmin = portion(5)-extents(5)+1;	kclmax = portion(6)-extents(5)+1
			!print'(a,5i4,2i6,i4,18i4,l)', 'recv qqqq',rank_world,rank_realm,rank_olap,ndata,nbr, & 
			!							rank_graph2rank_world(sourceid+1),size(arecv),start_address,&
			!							iclmin,iclmax,jclmin,jclmax,kclmin,kclmax, & 
			!							portion,extents,recv_flag
			do kcell=kclmin,kclmax
			do jcell=jclmin,jclmax
			do icell=iclmin,iclmax
			do n = 1,npercell
				arecv(n,icell,jcell,kcell) = vbuf(pos)
				pos = pos + 1
			end do
			end do
			end do
			end do


		endif

		!Increment reading position of packed data
		start_address = start_address + ndata

	enddo
	! ----------------- Unpack data -----------------------------
           
end subroutine CPL_recv_xd

!-------------------------------------------------------------------

subroutine CPL_pack(unpacked,packed,realm,icmax_pack,icmin_pack,jcmax_pack, & 
	                                      jcmin_pack,kcmax_pack,kcmin_pack    )
	use coupler_module, only : CPL_CART_COMM,rank_cart,md_realm,cfd_realm, & 
	                           error_abort,CPL_GRAPH_COMM,myid_graph,realm_name
	implicit none

	integer, intent(in)											:: realm
	real(kind=kind(0.d0)),dimension(:,:,:,:), intent(in)		:: unpacked
	real(kind=kind(0.d0)),dimension(:),allocatable, intent(out)	:: packed

    ! Optional minimum and maximum values to pack
 	integer, intent(in), optional    :: icmax_pack,icmin_pack
 	integer, intent(in), optional    :: jcmax_pack,jcmin_pack
 	integer, intent(in), optional    :: kcmax_pack,kcmin_pack

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
			call CPL_olap_extents(coord,md_realm,extents,ncells)

			! Get offset of neighbouring processor
			pos = id_nbr * npercell * ncells

			!print*,'Pack',rank_cart,realm_name(realm),coord,nbr,id_nbr,extents,coord

		elseif (realm .eq. md_realm) then
			!Get own processor coordinates 
			call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr) 
			call CPL_olap_extents(coord,realm,gextents,ncells)
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
!! @author David Trevelyan
!!
!! Gets maximum and minimum cells for processor coordinates
!!
!! - Synopsis
!!
!!  - CPL_proc_extents(coord,realm,extents,ncells)
!!
!! - Input
!!
!!  - coord
!!   - processor cartesian coordinate (3 x integer) 
!!
!!  - realm
!!   - cfd_realm (1) or md_realm (2) (integer) 
!!
!! - Input/Output
!!  - NONE
!!
!! - Output
!!
!!  - extents
!!   - Six components array which defines processor extents
!!     xmin,xmax,ymin,ymax,zmin,zmax (6 x integer) 
!!
!!  - ncells (optional)
!!   - number of cells on processor (integer) 
!!
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

!-------------------------------------------------------------------
! 					CPL_olap_extents  			      -
!-------------------------------------------------------------------
!! @author David Trevelyan
!!
!! Get maximum and minimum cells for current communicator within
!! the overlapping region only
!!
!! - Synopsis 
!!
!!  - CPL_olap_extents(coord,realm,extents,ncells)
!!
!! - Input 
!!
!!  - coord
!!   - processor cartesian coordinate (3 x integer) 
!!
!!  - realm
!!   - cfd_realm (1) or md_realm (2) (integer) 
!!
!! - Input/Output
!!  - NONE
!!
!! - Output 
!!
!!  - extents
!!   - Six components array which defines processor extents within
!!     the overlap region only: xmin,xmax,ymin,ymax,zmin,zmax (6 x integer) 
!!
!!  - ncells (optional)
!!   - number of cells on processor (integer) 
!!
subroutine CPL_olap_extents(coord,realm,extents,ncells)
	use mpi
	use coupler_module, only: md_realm,      cfd_realm,      &
	                          icPmin_md,     icPmax_md,      &
	                          jcPmin_md,     jcPmax_md,      &
	                          kcPmin_md,     kcPmax_md,      &
	                          icPmin_cfd,    icPmax_cfd,     &
	                          jcPmin_cfd,    jcPmax_cfd,     &
	                          kcPmin_cfd,    kcPmax_cfd,     &
	                          icmin_olap,    icmax_olap,     & 
	                          jcmin_olap,    jcmax_olap,     & 
	                          kcmin_olap,    kcmax_olap,     & 
	                          error_abort
	implicit none

	integer, intent(in)  :: coord(3), realm
	integer, intent(out) :: extents(6)
	integer, optional, intent(out) :: ncells

	select case(realm)
	case(md_realm)
		extents(1) = max(icPmin_md(coord(1)),icmin_olap)
		extents(2) = min(icPmax_md(coord(1)),icmax_olap) 
		extents(3) = max(jcPmin_md(coord(2)),jcmin_olap) 
		extents(4) = min(jcPmax_md(coord(2)),jcmax_olap) 
		extents(5) = max(kcPmin_md(coord(3)),kcmin_olap) 
		extents(6) = min(kcPmax_md(coord(3)),kcmax_olap) 
	case(cfd_realm)
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
! 					CPL_proc_portion  			      -
!-------------------------------------------------------------------
!! @author David Trevelyan
! Get maximum and minimum cells (a 'portion') based on the INPUT 
! limits within the current processor and overlapping region

! - - - Synopsis - - -

! CPL_proc_portion(coord,realm,limits,portion,ncells)

! - - - Input Parameters - - -

!coord
!    processor cartesian coordinate (3 x integer) 
!realm
!    cfd_realm (1) or md_realm (2) (integer) 
!limits
!	 Six components array which, together with processor and
!    overlap extents, determines the returned extents (6 x integer) 
!
! - - - Output Parameter - - -

!extents
!	 Six components array which defines processor extents within
!    the overlap region only 
!	 xmin,xmax,ymin,ymax,zmin,zmax (6 x integer) 
!ncells (optional)
!    number of cells in portion (integer) 


subroutine CPL_proc_portion(coord,realm,limits,portion,ncells)
	use mpi
	use coupler_module, only: md_realm,      cfd_realm,      &
	                          error_abort
	implicit none

	integer, intent(in)  :: coord(3), limits(6),realm
	integer, intent(out) :: portion(6)
	integer, optional, intent(out) :: ncells
	integer :: extents(6)

	call CPL_proc_extents(coord,realm,extents)

	if (extents(1).gt.limits(2) .or. &
		extents(2).lt.limits(1) .or. &
		extents(3).gt.limits(4) .or. &
	    extents(4).lt.limits(3) .or. &
	    extents(5).gt.limits(6) .or. &
	    extents(6).lt.limits(5)) then

		portion(:) = VOID
		if(present(ncells)) ncells = 0
		return
		
	end if

	portion(1) = max(extents(1),limits(1))							
	portion(2) = min(extents(2),limits(2))							
	portion(3) = max(extents(3),limits(3))							
	portion(4) = min(extents(4),limits(4))							
	portion(5) = max(extents(5),limits(5))							
	portion(6) = min(extents(6),limits(6))							

	if (present(ncells)) then
		ncells = (portion(2) - portion(1) + 1) * &
				 (portion(4) - portion(3) + 1) * &
				 (portion(6) - portion(5) + 1)
	end if

end subroutine CPL_proc_portion 

!-------------------------------------------------------------------
! 					CPL_Cart_coords								   -
!-------------------------------------------------------------------

!! Determines process coords in appropriate realm's cartesian topology 
!! given a rank in any communicator
!!
!! - Synopsis
!!
!!  - CPL_Cart_coords(COMM, rank, realm, maxdims, coords, ierr)
!!
!! - Input Parameters
!!
!!  - comm
!!   - communicator with cartesian structure (handle) 
!!
!!  - realm
!!   - cfd_realm (1) or md_realm (2) (integer) 
!!
!!  - rank
!!   - rank of a process within group of comm (integer) 
!!      NOTE fortran convention rank=1 to nproc
!!
!!  - maxdims
!!   - length of vector coords in the calling program (integer) 
!!
!! - Output Parameter
!!
!!  - coords
!!   - integer array (of size ndims) containing the Cartesian coordinates 
!!     of specified process (integer) 
!!
!!  - ierr
!!   - error flag
!! @author Edward Smith

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

!-------------------------------------------------------------------
! 					CPL_get_rank					   -
!-------------------------------------------------------------------
!! Return rank of current processor in specified COMM 
!!
!! - Synopsis
!!
!!  - CPL_get_rank(COMM, rank)
!!
!! - Input Parameters
!!
!!  - comm
!!   - communicator with cartesian structure (handle) 
!!
!! - Output Parameter
!!
!!  - rank
!!   - rank of a process within group of comm (integer) 
!!      NOTE fortran convention rank=1 to nproc
!!
!! @author Edward Smith

subroutine CPL_get_rank(COMM,rank)
	use coupler_module, only :  CPL_WORLD_COMM, CPL_REALM_COMM, CPL_INTER_COMM, & 
								CPL_CART_COMM,  CPL_OLAP_COMM,  CPL_GRAPH_COMM, &
								CPL_REALM_INTERSECTION_COMM,rank_world,			&
								rank_realm,rank_cart,rank_olap,rank_graph,error_abort

	integer, intent(in)		:: COMM
	integer, intent(out) 	:: rank

	!Get rank in world COMM from current COMM
	if (COMM .eq. CPL_WORLD_COMM) then
		rank = rank_world
		return
	elseif(COMM .eq. CPL_REALM_COMM) then
		rank = rank_realm
		return
	elseif(COMM .eq. CPL_CART_COMM) then
		rank = rank_cart
		return
	elseif(COMM .eq. CPL_OLAP_COMM) then
		rank = rank_olap
		return
	elseif(COMM .eq. CPL_GRAPH_COMM) then
		rank = rank_graph
		return
	elseif(COMM .eq. CPL_REALM_INTERSECTION_COMM) then
		call error_abort("CPL_Cart_coords Error - Intersection COMM not programmed")
	elseif(COMM .eq. CPL_INTER_COMM) then
		call error_abort("CPL_Cart_coords Error - No rank in Intercomm" // &
								 " - use realm comms instead")
	else 
		call error_abort("CPL_Cart_coords Error - Unknown COMM")
	endif
	

end subroutine CPL_get_rank


!-------------------------------------------------------------------
! 					CPL_olap_check										   -
!-------------------------------------------------------------------
!! Check if current processor is in the overlap region
!!
!! - Synopsis
!!
!!  - CPL_olap_check()
!!
!! - Input Parameters
!!
!!  - NONE
!!
!! - Returns
!!
!!  - CPL_olap_check
!!	 - True if calling processor is in the overlap region
!!     and false otherwise
!!
!! @author Edward Smith


function CPL_overlap() result(p)
	use coupler_module, only : olap_mask, rank_world
	implicit none

	logical	:: p

	p = olap_mask(rank_world)

end function CPL_overlap


!============================================================================
!
! Utility functions and subroutines that extract parameters from internal modules 
!
!-----------------------------------------------------------------------------    

!-------------------------------------------------------------------
! 					CPL_get										   -
!-------------------------------------------------------------------
!! Wrapper to retrieve (read only) parameters from the coupler_module 
!! Note - this ensures all variable in the coupler are protected
!! from corruption by either CFD or MD codes
!!
!! - Synopsis
!!
!!  - CPL_get([see coupler_module])
!!
!! - Input Parameters
!!
!!  - NONE
!!
!! - Output Parameter
!!
!!  - @see coupler_module
!!
!! @author Edward Smith


subroutine CPL_get(icmax_olap,icmin_olap,jcmax_olap,jcmin_olap,  & 
				   kcmax_olap,kcmin_olap,density_cfd,density_md, &
				   dt_cfd,dt_MD,dx,dy,dz,ncx,ncy,ncz,xg,yg,zg,	 &
				   xL_md,xL_cfd,yL_md,yL_cfd,zL_md,zL_cfd       )
	use coupler_module, only : 	icmax_olap_=>icmax_olap, &
							    icmin_olap_=>icmin_olap, &
								jcmax_olap_=>jcmax_olap, &		
								kcmax_olap_=>kcmax_olap, &		
								density_cfd_=>density_cfd, &		
								jcmin_olap_=>jcmin_olap, &		
								kcmin_olap_=>kcmin_olap, &		
								density_md_=>density_md, &		
								dt_cfd_=>dt_cfd,dt_MD_=>dt_MD, &	
								dx_=>dx,dy_=>dy,dz_=>dz, &
								ncx_=>ncx,ncy_=> ncy,ncz_=> ncz, &
								xL_md_ =>xL_md, &		
								yL_md_ =>yL_md, &		
								zL_md_ =>zL_md, &		
								xL_cfd_=> xL_cfd, &		
								yL_cfd_=> yL_cfd, &		
								zL_cfd_=> zL_cfd, &		
								xg_=>xg, yg_=> yg, zg_=> zg
	implicit none

	integer, optional, intent(out)			:: icmax_olap ,icmin_olap
	integer, optional, intent(out)			:: jcmax_olap ,jcmin_olap
	integer, optional, intent(out)			:: kcmax_olap ,kcmin_olap
	integer, optional, intent(out)			:: ncx,ncy,ncz
	real(kind(0.d0)), optional, intent(out)	:: density_cfd,density_md
	real(kind(0.d0)), optional, intent(out) :: dt_cfd,dt_MD
	real(kind(0.d0)), optional, intent(out) :: dx,dy,dz
	real(kind(0.d0)), optional, intent(out) :: xL_md,xL_cfd
	real(kind(0.d0)), optional, intent(out) :: yL_md,yL_cfd
	real(kind(0.d0)), optional, intent(out) :: zL_md,zL_cfd
	real(kind(0.d0)),dimension(:,:),allocatable,optional,intent(out) :: xg,yg
	real(kind(0.d0)),dimension(:)  ,allocatable,optional,intent(out) :: zg

	!Overlap extents
	if (present(icmax_olap)) icmax_olap = icmax_olap_
	if (present(icmin_olap)) icmin_olap = icmin_olap_
	if (present(jcmax_olap)) jcmax_olap = jcmax_olap_
	if (present(jcmin_olap)) jcmin_olap = jcmin_olap_
	if (present(kcmax_olap)) kcmax_olap = kcmax_olap_
	if (present(kcmin_olap)) kcmin_olap = kcmin_olap_

	!Density
	if (present(density_cfd)) density_cfd = density_cfd_
	if (present(density_md )) density_md  = density_md_

	!Cells
	if (present(ncx)) ncx = ncx_
	if (present(ncy)) ncy = ncy_
	if (present(ncz)) ncz = ncz_
			
	!Temporal and spatial steps
	if (present(dt_cfd)) dt_cfd= dt_cfd_
	if (present(dt_MD )) dt_MD = dt_MD_
	if (present(dx)) dx = dx_
	if (present(dy)) dy = dy_	
	if (present(dz)) dz = dz_

	!Domain sizes
	if (present(xL_md )) xL_md = xL_md_
	if (present(xL_cfd)) xL_cfd= xL_cfd_
	if (present(yL_md )) yL_md = yL_md_
	if (present(yL_cfd)) yL_cfd= yL_cfd_
	if (present(zL_md )) zL_md = zL_md_
	if (present(zL_cfd)) zL_cfd= zL_cfd_
			
	!The mesh
	if (present(xg)) then
		allocate(xg(size(xg_,1),size(xg_,2)))
		xg = xg_
	endif
	if (present(yg)) then
		allocate(yg(size(yg_,1),size(yg_,2)))
		yg = yg_
	endif
	if (present(zg)) then
		allocate(zg(size(zg_,1)))
		zg = zg_
	endif

end subroutine CPL_get


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
	use coupler_module, only : timestep_ratio 
	implicit none

	integer p
	p = timestep_ratio 

end function coupler_md_get_md_steps_per_cfd_dt

!-----------------------------------------------------------------------------

function coupler_md_get_nsteps() result(p)
	use coupler_module, only :  nsteps_md
	implicit none 

	 integer p
	 p = nsteps_md

end function coupler_md_get_nsteps

!-----------------------------------------------------------------------------

function coupler_md_get_dt_cfd() result(p)
	use coupler_module, only : dt_CFD  
	implicit none

	real(kind=kind(0.d0)) p
	p = dt_CFD

end function coupler_md_get_dt_cfd

!------------------------------------------------------------------------------

end module coupler
