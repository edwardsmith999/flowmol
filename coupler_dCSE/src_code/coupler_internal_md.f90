!=============================================================================
!				   Coupler internal MD
! Internal data and subroutines used by the coupler when working in MD realm
! It must not be used by subroutimes working in CFD realm ONLY
! Subroutines include:
! coupler_md_init           initialises coupler and set MD 
!							parameters with using data from CFD 
!							or COUPLER.in
! create_map_md 			Establish for all MD processors the mapping (if any) 
! 							to coupled CFD processors
! make_bbox					Make bbox which contains all domain information 
!							required in coupling process such as domain extents
! get_CFDvel				Get velocity fields from CFD for the constraint 
!							force needed in MD
! send_vel(a,n1,n2,n3)   	Send the average MD bin velocities to the 
!							corresponding rank in the CFD realm
! write_overlap_map			Write the map of coupled processors for debugging
! function map_md2cfd(r)	Get global molecular position from local position
! 
!  Lucian Anton, November 2011
!
!=============================================================================

module coupler_internal_md
	implicit none

contains 

!------------------------------------------------------------------------------
!                              coupler_md_init                               -
!------------------------------------------------------------------------------
!
!! Initialisation routine for coupler module - Every variable is sent and stored
!! to ensure both md and cfd region have an identical list of parameters
!!
!! - Synopsis
!!
!!  - coupler_mf_init(nsteps,dt_md,icomm_grid,icoord,npxyz_md,globaldomain,density)
!!
!! - Input
!!
!!  - nsteps
!!   - Number of time steps the MD code is expected to run for (integer)
!!  - dt_md
!!   - MD timestep (dp real)
!!  - icomm_grid
!!   - The MPI communicator setup by the MPI_CART_CREATE command in the 
!!     CFD region (integer)
!!  - icoord
!!   - The three coordinate for each rank in the domain (integer array nproc by 3)
!!  - npxyz_md
!!   - Number of processors in each cartesian dimension (integer array 3)
!!  - globaldomain
!!   - Size of domain in each cartesian dimension (dp real array 3)
!!  - density
!!   - Density of the CFD simulation (dp_real)
!!
!! - Input/Output
!!  - NONE
!!
!! - Output
!!  - NONE
!! 
!! @author Edward Smith
!
! ----------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! Initialisation routine for coupler - Every variable is sent and stored
! to ensure both md and cfd region have an identical list of parameters

subroutine coupler_md_init(nsteps,dt_md,icomm_grid,icoord,npxyz_md,globaldomain,density)
	use mpi
	use coupler, only : CPL_rank_map
    use coupler_module, dt_md_=>dt_md		
	implicit none

	integer, intent(in)	  							:: nsteps, icomm_grid
	integer,dimension(3), intent(in)	  			:: npxyz_md	
	integer,dimension(:,:),allocatable,intent(in)	:: icoord
	real(kind(0.d0)),intent(in) 					:: dt_md,density
    real(kind=kind(0.d0)),dimension(3),intent(in) 	:: globaldomain

    integer											:: i,ib,jb,kb,pcoords(3),source,nproc
    integer,dimension(:),allocatable  				:: buf,rank_world2rank_realm,rank_world2rank_cart
    real(kind=kind(0.d0)),dimension(:),allocatable 	:: rbuf

    ! Duplicate grid communicator for coupler use
    call MPI_comm_dup(icomm_grid,CPL_CART_COMM,ierr)
    call MPI_comm_rank(CPL_CART_COMM,myid_cart,ierr) 
    rank_cart = myid_cart + 1; rootid_cart = 0	
	!Send only from root processor
	if ( myid_realm .eq. rootid_realm ) then
		source=MPI_ROOT
	else
		source=MPI_PROC_NULL
	endif

	! ================ Exchange and store Data ==============================
	! Data is stored to the coupler module with the same name in both realms
	! Note - MPI Broadcast between intercommunicators is only supported by MPI-2

	! ------------------------ Processor Topology ---------------------------
	! Receive & Store CFD number of processors
    allocate(buf(3))
	call MPI_bcast(  buf   ,3,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr) !Receive
	npx_cfd = buf(1); npy_cfd = buf(2); npz_cfd = buf(3)
	nproc_cfd = npx_cfd * npy_cfd * npz_cfd
	deallocate(buf)

	write(999+rank_realm,*), 'MD side',rank_realm,'CFD procs',npx_cfd,npy_cfd,npz_cfd,nproc_cfd

	! Store & Send MD number of processors
	npx_md = npxyz_md(1);	npy_md = npxyz_md(2);	npz_md = npxyz_md(3)	
	nproc_md = npx_md * npy_md * npz_md
	call MPI_bcast(npxyz_md,3,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send

	write(999+rank_realm,*), 'MD side',rank_realm,'MD procs',npx_md ,npy_md ,npz_md ,nproc_md

	! Receive & Store CFD processor rank to coord
	allocate(buf(3*nproc_cfd))
    call MPI_bcast(buf,3*nproc_cfd,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr) !Receive
    allocate(rank2coord_cfd(3,nproc_cfd),stat=ierr); rank2coord_cfd = reshape(buf,(/ 3,nproc_cfd /))
	deallocate(buf)

	call write_matrix_int(rank2coord_cfd,'MD side, rank2coord_cfd=',999+rank_realm)

	! Store & Send MD processor rank to coord
    allocate(rank2coord_md(3,nproc_md),stat=ierr); rank2coord_md = icoord
	iblock_realm=icoord(1,rank_realm); jblock_realm=icoord(2,rank_realm); kblock_realm=icoord(3,rank_realm)
	allocate(buf(3*nproc_md)); buf = reshape(icoord,(/ 3*nproc_md /))
    call MPI_bcast(buf ,3*nproc_md,MPI_INTEGER,source,CPL_INTER_COMM,ierr) !Send
	deallocate(buf)

	call write_matrix_int(rank2coord_md,'MD side, rank2coord_md=',999+rank_realm)

	! Receive & Store CFD coordinate to rank mapping
	allocate(buf(nproc_cfd))
    call MPI_bcast(buf,nproc_cfd ,MPI_INTEGER,  0   ,CPL_INTER_COMM,ierr)	!Receive
    allocate(coord2rank_cfd (npx_cfd,npy_cfd,npz_cfd))
	coord2rank_cfd = reshape(buf,(/ npx_cfd,npy_cfd,npz_cfd /))
	deallocate(buf)

	write(999+rank_realm,*), 'MD side',rank_realm, 'coord2rank_cfd=', coord2rank_cfd

	! Setup MD mapping from coordinate to rank, 
	! Store & Send MD mapping from coordinate to rank to CFD
	allocate(coord2rank_md(npx_md,npy_md,npz_md))
	do ib = 1,npx_md
	do jb = 1,npy_md
	do kb = 1,npz_md
		pcoords = (/ ib, jb, kb /)-1
		call MPI_Cart_rank(CPL_CART_COMM,pcoords,i,ierr)
		coord2rank_md(ib,jb,kb) = i + 1
	enddo
	enddo
	enddo
	allocate(buf(nproc_md)); buf = reshape(coord2rank_md, (/ nproc_md /) )
    call MPI_bcast(coord2rank_md,nproc_md,MPI_INTEGER,source,CPL_INTER_COMM,ierr)	!Send
	deallocate(buf)

	write(999+rank_realm,*), 'MD side',rank_realm, 'coord2rank_md=', coord2rank_md

	! Setup MD mapping between realm & world rank 
	allocate(rank_mdrealm2rank_world(nproc_md))
	allocate(rank_world2rank_realm(nproc_world))
	call CPL_rank_map(CPL_REALM_COMM,rank_realm,nproc, & 
						rank_mdrealm2rank_world,rank_world2rank_realm,ierr)

	!World to rank is the same on both realms
	allocate(rank_world2rank_cfdrealm(nproc_world))
	allocate(rank_world2rank_mdrealm(nproc_world))
	rank_world2rank_cfdrealm = rank_world2rank_realm
	rank_world2rank_mdrealm  = rank_world2rank_realm

	! Receive & Store CFD_mapping from realm to local rank
	allocate(rank_cfdrealm2rank_world(nproc_cfd))
	call MPI_bcast(rank_cfdrealm2rank_world,nproc_cfd,MPI_integer,0,CPL_INTER_COMM,ierr)	!Receive

	write(999+rank_realm,*), 'MD side',rank_realm, 'rank_cfdrealm2rank_world', rank_cfdrealm2rank_world

	! Send MD mapping from realm to local rank to CFD
	call MPI_bcast(rank_mdrealm2rank_world,nproc_md,MPI_integer,source,CPL_INTER_COMM,ierr)	 !send

	write(999+rank_realm,*), 'MD side',rank_realm, 'rank_mdrealm2rank_world', rank_mdrealm2rank_world
	write(999+rank_realm,*), 'MD side',rank_realm, 'rank_world2rank_mdrealm', rank_world2rank_mdrealm

	! Setup MD mapping between cartesian topology & world rank 
	allocate(rank_mdcart2rank_world(nproc_md))
	allocate(rank_world2rank_cart(nproc_world))
	call CPL_rank_map(CPL_CART_COMM,rank_cart,nproc, & 
						rank_mdcart2rank_world,rank_world2rank_cart,ierr)

	!World to rank is the same on both realms cart
	allocate(rank_world2rank_cfdcart(nproc_world))
	allocate(rank_world2rank_mdcart(nproc_world))
	rank_world2rank_cfdcart = rank_world2rank_cart
	rank_world2rank_mdcart  = rank_world2rank_cart

	! Receive & Store CFD_mapping from cart to local rank
	allocate(rank_cfdcart2rank_world(nproc_cfd))
	call MPI_bcast(rank_cfdcart2rank_world,nproc_cfd,MPI_integer,0,CPL_INTER_COMM,ierr)	!Receive

	write(999+rank_realm,*), 'MD side',rank_cart, 'rank_cfdcart2rank_world', rank_cfdcart2rank_world

	! Send MD mapping from cart to local rank to CFD
	call MPI_bcast(rank_mdcart2rank_world,nproc_md,MPI_integer,source,CPL_INTER_COMM,ierr)	 !send

	write(999+rank_realm,*), 'MD side',rank_cart, 'rank_mdcart2rank_world', rank_mdcart2rank_world
	write(999+rank_realm,*), 'MD side',rank_cart, 'rank_world2rank_mdcart', rank_world2rank_mdcart

	! ------------------ Timesteps and iterations ------------------------------
	! Receive & store CFD nsteps and dt_cfd
	call MPI_bcast(nsteps_cfd,1,MPI_integer,0,CPL_INTER_COMM,ierr)				!Receive
	call MPI_bcast(dt_cfd,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive

	! Store & send MD timestep to dt_md
	dt_MD_ = dt_MD
    call MPI_bcast(dt_md,1,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send
	nsteps_MD = nsteps
    call MPI_bcast(nsteps,1,MPI_integer,source,CPL_INTER_COMM,ierr)	!Send

	write(999+rank_realm,*), 'MD side',rank_realm,'CFD times',nsteps_cfd,dt_cfd
	write(999+rank_realm,*), 'MD side',rank_realm,'MD times', nsteps_MD,dt_md

	! ------------------ Receive CFD grid extents ------------------------------
	! Receive & store CFD density
	call MPI_bcast(density_cfd,1,MPI_double_precision,0,CPL_INTER_COMM,ierr)		!Receive

	write(999+rank_realm,*), 'MD side',rank_realm,'CFD density',density_cfd

	! Store & send MD density
	density_md = density 
	call MPI_bcast(density,1,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	write(999+rank_realm,*), 'MD side',rank_realm,'MD density',density_md

	! Receive & store CFD domain size
	allocate(rbuf(3))
	call MPI_bcast(rbuf,3,MPI_double_precision,0,CPL_INTER_COMM,ierr)				!Receive
	xL_cfd = rbuf(1); yL_cfd = rbuf(2); zL_cfd = rbuf(3)
	deallocate(rbuf)

	! Store & send MD domain size
	xL_md = globaldomain(1); yL_md = globaldomain(2); zL_md = globaldomain(3) 
	call MPI_bcast(globaldomain,3,MPI_double_precision,source,CPL_INTER_COMM,ierr)	!Send

	write(999+rank_realm,*), 'MD side',rank_realm,'CFD domain',xL_cfd,yL_cfd,zL_cfd
	write(999+rank_realm,*), 'MD side',rank_realm,'MD domain', xL_md,yL_md,zL_md


	! Receive & Store global CFD grid extents
	allocate(buf(6))
	call MPI_bcast(buf, 6, MPI_INTEGER, 0, CPL_INTER_COMM,ierr) !Send
	icmin = buf(1); icmax = buf(2)
	jcmin = buf(3); jcmax = buf(4)
	kcmin = buf(5); kcmax = buf(6)
	deallocate(buf)

	! Receive & Store array of global number of cells in CFD
	allocate(buf(3))
	call MPI_bcast(buf, 3, MPI_INTEGER, 0, CPL_INTER_COMM,ierr) !Receive
	ncx = buf(1); ncy = buf(2); ncz = buf(3)
	deallocate(buf)

	write(999+rank_realm,*), 'MD side',rank_realm,'CFD global extents',icmin,icmax,jcmin,jcmax,kcmin,kcmax 
	write(999+rank_realm,*), 'MD side',rank_realm,'CFD global cells',ncx,ncy,ncz

	! Receive & Store array of global grid points
	allocate(xg(ncx+1,ncy+1),yg(ncx+1,ncy+1),zg(ncz+1))
	call MPI_bcast(xg,size(xg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive
	call MPI_bcast(yg,size(yg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive
	call MPI_bcast(zg,size(zg),MPI_double_precision,0,CPL_INTER_COMM,ierr) !Receive

	!call write_matrix(xg,'MD side, xg=',500+rank_realm)
	!call write_matrix(yg,'MD side, yg=',500+rank_realm)
	!write(500+rank_realm,*), 'MD side',rank_realm,'zg',zg

	! Receive & Store local (processor) CFD grid extents
    allocate(icPmin_cfd(npx_cfd)); allocate(icPmax_cfd(npx_cfd));  
	allocate(jcPmin_cfd(npy_cfd)); allocate(jcPmax_cfd(npy_cfd));
    allocate(kcPmin_cfd(npz_cfd)); allocate(kcPmax_cfd(npz_cfd));
    call MPI_bcast(icPmin_cfd,npx_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(icPmax_cfd,npx_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(jcPmin_cfd,npy_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(jcPmax_cfd,npy_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(kcPmin_cfd,npz_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive
    call MPI_bcast(kcPmax_cfd,npz_cfd,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive

	write(999+rank_realm,*), 'md side',rank_realm,'CFD local cells', & 
							icPmin_cfd,icPmax_cfd,jcPmin_cfd,jcPmax_cfd,kcPmin_cfd,kcPmax_cfd

	!Calculate the cell sizes dx,dy & dz
	dx = xL_cfd/ncx	  !xg(2,1)-xg(1,1)
	dy = yg(1,2)-yg(1,1) ! yL_cfd/ncy
	dz = zL_cfd/ncz	  !zg(2  )-zg(1  )

	ncx_olap = icmax_olap - icmin_olap + 1
	ncy_olap = jcmax_olap - jcmin_olap + 1
	ncz_olap = kcmax_olap - kcmin_olap + 1

    ! Initialise other md module variables if data is provided in coupler.in
!	if (md_average_period_tag == CPL) then 
!		average_period = md_average_period
!	endif
!	if (md_save_period_tag == CPL) then
!		save_period    = md_save_period
!	endif

	!if (cfd_code_id_tag == CPL) then
	!	cfd_code_id  = cfd_code_id_in
	!endif
    
	if ( nsteps_md <= 0 ) then 
		write(0,*) "Number of MD steps per dt interval <= 0"
		write(0,*) "Coupler will not work, quitting ..."
		call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_INIT,ierr)
	endif

	!Receive overlap from the CFD
    !call MPI_bcast(ncy_olap,1,MPI_INTEGER,0,CPL_INTER_COMM,ierr) !Receive

	write(999+rank_realm,*), 'MD side - y overlap',ncy_olap

	! ================ Apply domain setup  ==============================
	! --- set the sizes of the MD domain ---
	!Fix xL_md domain size to continuum
	!xL_md = x(icmax_cfd) - x(icmin_cfd)

    ! yL_md is adjusted to an integer number of initialisation cells in the following steps
   ! if (md_ly_extension_tag == CPL) then 
    !    DY_PURE_MD = md_ly_extension
    !else 
    !    DY_PURE_MD = y(jcmin) - y(jcmino)
    !end if
    !yL_md = y(jcmax_overlap) - y(jcmino) + DY_PURE_MD
    !yL_md = real(floor(yL_md/b),kind(0.d0))*b
    !DY_PURE_MD = yL_md - (y(jcmax_overlap) - y(jcmino))

	!Fix zL_md domain size to continuum
	!zL_md = z(kcmax_cfd) - z(kcmin)

end subroutine coupler_md_init

!=============================================================================
! Get molecule's global position from position local to processor.
!-----------------------------------------------------------------------------
function globalise(r) result(rg)
	use coupler_module, only : 	xLl,iblock_realm,npx_md, & 
								yLl,jblock_realm,npy_md, & 
								zLl,kblock_realm,npz_md
	implicit none

	real(kind(0.d0)),intent(in) :: r(3)
	real(kind(0.d0)) rg(3)

	rg(1) = r(1) - xLl*(iblock_realm-1)+0.5d0*xLl*(npx_md-1)
	rg(2) = r(2) - yLl*(jblock_realm-1)+0.5d0*yLl*(npy_md-1)
	rg(3) = r(3) - zLl*(kblock_realm-1)+0.5d0*zLl*(npz_md-1)

end function globalise

!=============================================================================
! Get local position on processor from molecule's global position.
!-----------------------------------------------------------------------------
function localise(r) result(rg)
	use coupler_module, only : 	xLl,iblock_realm,npx_md, & 
								yLl,jblock_realm,npy_md, & 
								zLl,kblock_realm,npz_md
	implicit none

	real(kind(0.d0)),intent(in) :: r(3)
	real(kind(0.d0)) rg(3)

	!Global domain has origin at centre
	rg(1) = r(1) - xLl*(iblock_realm-1)+0.5d0*xLl*(npx_md-1)
	rg(2) = r(2) - yLl*(jblock_realm-1)+0.5d0*yLl*(npy_md-1)
	rg(3) = r(3) - zLl*(kblock_realm-1)+0.5d0*zLl*(npz_md-1)

end function localise

!=============================================================================
! Map global MD position to global CFD coordinate frame
!-----------------------------------------------------------------------------
function map_md2cfd_global(r) result(rg)
	use coupler_module, only : 	xL_md,xg,icmin_olap,icmax_olap, & 
								yL_md,yg,jcmin_olap,jcmax_olap, & 
								zL_md,zg,kcmin_olap,kcmax_olap
	implicit none

	real(kind(0.d0)),intent(in) :: r(3)
	real(kind(0.d0)):: md_only(3), rg(3)

	!Get size of MD domain which has no CFD cells overlapping
	!This should be general enough to include grid stretching
	!and total overlap in any directions 
	md_only(1) = xL_md-(xg(icmax_olap+1,1) - xg(icmin_olap,1))
	md_only(2) = yL_md-(yg(1,jcmax_olap+1) - yg(1,jcmin_olap))
	md_only(3) = zL_md-(zg( kcmax_olap+1 ) - zg( kcmin_olap ))

	! CFD has origin at bottom left while MD origin at centre
	rg(1) = r(1) + 0.5d0*xL_md - md_only(1)
	rg(2) = r(2) + 0.5d0*yL_md - md_only(2)
	rg(3) = r(3) + 0.5d0*zL_md - md_only(3)

end function map_md2cfd_global


!=============================================================================
! Map global CFD position in global MD coordinate frame
!-----------------------------------------------------------------------------
function map_cfd2md_global(r) result(rg)
	use coupler_module, only : 	xL_md,xg,icmin_olap,icmax_olap, & 
								yL_md,yg,jcmin_olap,jcmax_olap, & 
								zL_md,zg,kcmin_olap,kcmax_olap
	implicit none

	real(kind(0.d0)),intent(in) :: r(3)
	real(kind(0.d0)) :: md_only(3), rg(3)

	!Get size of MD domain which has no CFD cells overlapping
	!This should be general enough to include grid stretching
	!and total overlap in any directions 
	md_only(1) = xL_md-(xg(icmax_olap+1,1) - xg(icmin_olap,1))
	md_only(2) = yL_md-(yg(1,jcmax_olap+1) - yg(1,jcmin_olap))
	md_only(3) = zL_md-(zg( kcmax_olap+1 ) - zg( kcmin_olap ))

	! CFD has origin at bottom left while MD origin at centre
	rg(1) = r(1) - 0.5d0*xL_md + md_only(1)
	rg(2) = r(2) - 0.5d0*yL_md + md_only(2)
	rg(3) = r(3) - 0.5d0*zL_md + md_only(3)

	!print'(a,13f8.3)', 'md only', r,md_only,rg,yL_md,(yg(1,jcmax_olap+1),yg(1,jcmin_olap)),(yg(1,jcmax_olap+1) - yg(1,jcmin_olap))

end function map_cfd2md_global


end module coupler_internal_md
