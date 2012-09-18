! A single coupler module for both codes - this contains the same information 
! on both md and cfd side 
module coupler_module

	! local mpi ranks in realm and in topology communicators
	!integer myid, myid_grid
    ! communicator for tasks that overlap  MD region 
    integer CFD_COMM_OVERLAP

	!Total number of steps for coupled simulation
	integer	:: nsteps

	! CFD number of cells
    integer	:: ngx,ngy,ngz
	! CFD grid indices
    integer	:: imin,imax,jmin,jmax,kmin,kmax 
	! CFD local grid indices (start and end of grid per CFD process)
    integer,dimension(:),allocatable :: iTmin,iTmax,jTmin,jTmax,kTmin,kTmax

	! Number of processor in CFD grid
	integer :: npx_cfd, npy_cfd, npz_cfd, nproc_cfd	
    ! Number of processor in MD grid
    integer :: npx_md,  npy_md,  npz_md,  nproc_md


	! COUPLER VARIABLES
    integer :: jmax_overlap = 5 ! maximum j index ( in y direction) which MD 
	real(kind(0.d0)) :: MD_initial_cellsize

    ! Coordinates of MD/CFD topologies
    ! ATTENTION the values are shifted with +1, FORTRAN style
    ! if this array is passed by an MPI function, remove the shift
    integer,dimension(:,:),allocatable :: icoord_cfd, icoord_md

	!Timesteps
	real(kind(0.d0)) :: dt_md , dt_cfd

	!Domain sizes
	real(kind(0.d0)) ::	md_xL,md_yL,md_zL,cfd_xL,cfd_yL,cfd_zL

	!CFD grid 
    real(kind(0.d0)),dimension(:),  allocatable :: zpg
    real(kind(0.d0)),dimension(:,:),allocatable :: xpg,ypg

end module coupler_module
