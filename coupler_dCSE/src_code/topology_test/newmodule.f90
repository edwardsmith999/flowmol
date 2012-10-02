!=============================================================================
! COUPLER MODULE: 
!
!	- Error handling
!	- MPI communicators
!	- Simulation realms
!	- MPI processor IDs
!	- Processor topologies
!	- Processor cartesian coords
!	- Global cell grid parameters
!	- Processor cell ranges
!	- Domain and cell dimensions
!	- Positions of CFD grid lines
!	- CFD to MD processor mapping
!	- Simulation parameters
!   - PROBABLY OBSOLETE STUFF
!
!=============================================================================

module coupler_module
!	use coupler_parameters

	!REMOVE THIS, ONLY REQUIRED FOR TESTING!TODO
	integer :: testval

	! Error handling
	integer :: ierr
	
	! MPI Communicators
	integer :: &
		CPL_WORLD_COMM,   &
		CPL_INTER_COMM,   &
		CPL_REALM_COMM,   &
		CPL_CART_COMM,    &
		CPL_OLAP_COMM,    &
		CPL_REALM_INTERSECTION_COMM

	! Simulation realms
	integer :: &
		realm
	integer, parameter :: &
		cfd_realm = 1,    &
		md_realm  = 2

	! MPI processor IDs
	integer :: &
		myid_world,       &
		rank_world,       &
		rootid_world,     &
		myid_realm,       &
		rank_realm,       &
		rootid_realm,     &
		myid_cart,        &
		rank_cart,        &
		rootid_cart,      &
		myid_olap,        &
		rank_olap,        &
		rootid_olap
	integer, dimension(:), allocatable :: &
		rank_cfd2rank_world, &
		rank_md2rank_world

	! Processor topologies
	integer :: &
		nproc_world,      &
		nproc_md,         &
		nproc_cfd,        &
		npx_md,           &
		npy_md,           &
		npz_md,           &
		npx_cfd,          &
		npy_cfd,          &
		npz_cfd
	integer, dimension(:,:), allocatable :: &
		rank2coord_cfd,   &
		rank2coord_md
	integer, dimension(:,:,:), allocatable :: &
		coord2rank_cfd,   &
		coord2rank_md

	! Processor cartesian coords	
	integer :: &
		iblock_realm,     &
		jblock_realm,     &
		kblock_realm

	! Global cell grid parameters
	integer :: &
		ncx,              &
		ncy,              &
		ncz,              &
		icmin,            &
		icmax,            &
		jcmin,            &
		jcmax,            &
		kcmin,            &
		kcmax,            &
		icmin_olap,       &
		icmax_olap,       &
		jcmin_olap,       &
		jcmax_olap,       &
		kcmin_olap,       &
		kcmax_olap
	
	! Processor cell ranges 
	integer, dimension(:), allocatable :: &
		icPmin_md,        &
		icPmax_md,        &
		jcPmin_md,        &
		jcPmax_md,        &
		kcPmin_md,        &
		kcPmax_md,        &
		icPmin_cfd,       &
		icPmax_cfd,       &
		jcPmin_cfd,       &
		jcPmax_cfd,       &
		kcPmin_cfd,       &
		kcPmax_cfd
	
	! Domain and cell dimensions
	real(kind(0.d0)) :: &
		xL_md,            &
		yL_md,            &
		zL_md,            &
		xL_cfd,           &
		yL_cfd,           &
		zL_cfd,           &
		xL_olap,          &
		yL_olap,          &
		zL_olap,          &
		dx,               &
		dy,               &
		dz,               &
		dymin,            &
		dymax

	! Positions of CFD grid lines
	real(kind(0.d0)), dimension(:,:), allocatable :: &
		xg,               &
		yg
	real(kind(0.d0)), dimension(:)  , allocatable :: &
		zg

	! CFD to MD processor mapping
	integer, dimension(:,:), allocatable :: &
		cfd_icoord2olap_md_icoords, &
		cfd_jcoord2olap_md_jcoords, &
		cfd_kcoord2olap_md_kcoords

	! Simulation parameters
	integer :: &
		nsteps_md,        &
		nsteps_cfd,       &
		nsteps_coupled,   &
		average_period,   &
		save_period
	real(kind(0.d0)) :: &
		dt_md,            &
		dt_cfd,           &
		density_md,       &
		density_cfd

	! PROBABLY OBSOLETE STUFF ------------------------------------------------!	
	real(kind(0.d0)) :: MD_initial_cellsize                                   !
    type overlap_map                                                          !
        integer                             :: n                              ! 
        integer,dimension(:), allocatable   :: rank_list                      !
        integer,dimension(:,:),allocatable  :: domains                        ! 
    end type overlap_map                                                      !
    type(overlap_map)   :: map                                                !
	!-------------------------------------------------------------------------!

end module coupler_module
