!=============================================================================
!				   	CFD coupler Socket
! Routine which interface with the coupler to the MD code
!
! socket_coupler_init				Passes CFD initialisation variables to 
!									coupler_cfd_init
! socket_coupler_send_velocity		Send velocity from CFD for MD constraint
! socket_coupler_get_md_BC			Receive BC from MD
!
!=============================================================================


module continuum_coupler_socket

#if USE_COUPLER
	implicit none

contains

!=============================================================================
! Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup
!
!								SETUP
!
! Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup
!=============================================================================

!=============================================================================
! Invoke messenger and split inter-communicators
!-----------------------------------------------------------------------------
subroutine socket_coupler_invoke
	use messenger
	use CPL, only : CPL_create_comm, cfd_realm
	use computation_parameters, only : prefix_dir
	implicit none

	call CPL_create_comm(cfd_realm,CFD_COMM,ierr)
	prefix_dir ="./couette_data/"

end subroutine socket_coupler_invoke

!=============================================================================
! Call coupler initialise to swap all setup data
!-----------------------------------------------------------------------------
subroutine socket_coupler_init
	use CPL, only : coupler_cfd_init,CPL_create_map,error_abort, density_cfd,CPL_write_header
    use messenger, only : icomm_grid, icoord
    use data_export, only : imin,imax,jmin,jmax,kmin,kmax, &
							iTmin_1,iTmax_1,jTmin_1,jTmax_1, &
                            ngx,ngy,ngz,xL,yL,zL, &
                            npx,npy,npz,dt,xpg,ypg,zpg
	use mesh_export, only : xL, yL, zL
    implicit none

    integer							:: nsteps
	integer,dimension(1)			:: kTmin_1,kTmax_1
    integer,dimension(3)		   	:: ijkmin,ijkmax,npxyz,ngxyz
    real(kind(0.d0)),dimension(3)	:: xyzL

    call readInt("nsteps", nsteps)
	kTmin_1(1) = kmin; kTmax_1(1) = kmax-1

	!Define compound arrays to make passing more concise
	ijkmin = (/ imin, jmin, kmin /)	
	ijkmax = (/ imax, jmax, kmax /)-1 !Minus 1 as coupler requires max cells NOT max cell vertices
	npxyz  = (/ npx , npy , npz  /)
	ngxyz  = (/ ngx , ngy , ngz  /)-1 !Minus 1 as coupler requires no. of cells NOT no. cell vertices
	xyzL   = (/  xL ,  yL ,  zL  /)

    call coupler_cfd_init(nsteps,dt,icomm_grid,icoord,npxyz,xyzL,ngxyz, & 
	                      density_cfd,ijkmax,ijkmin,iTmin_1,iTmax_1,jTmin_1,&
	                      jTmax_1,kTmin_1,kTmax_1,xpg,ypg,zpg)

	!Write coupler information to header file
	call CPL_write_header('./results/coupler_header')

end subroutine socket_coupler_init


!=============================================================================
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!
!							SIMULATION
!
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!=============================================================================

!=============================================================================
!  __  __  ____     ___      ___  ____  ____     ____   ___ 
! (  \/  )(  _ \   (__ \    / __)( ___)(  _ \   (  _ \ / __)
!  )    (  )(_) )   / _/   ( (__  )__)  )(_) )   ) _ <( (__ 
! (_/\/\_)(____/   (____)   \___)(__)  (____/   (____/ \___)
!
! Get Boundary condition for continuum from average of MD 

subroutine  socket_coupler_get_md_BC(uc,vc,wc)
    use CPL, only : CPL_get,CPL_recv,CPL_realm,error_abort,cpl_proc_extents,VOID, & 
					printf,CPL_OLAP_COMM,xg,xL_cfd,yg,yL_cfd,zg,zL_cfd !DEBUG DEBUG
	use data_export, only : nixb,niyb,nlx,ngz,nizb,iblock,jblock,kblock, & 
							i1_u,i2_u,j1_u,j2_u, & 
							i1_v,i2_v,j1_v,j2_v, & 
							i1_w,i2_w,j1_w,j2_w, & 
							i1_T,i2_T,j1_T,j2_T,k1_T,k2_T, &
							ibmin,jbmin,ibmax,jbmax, &
							imap_1,jmap_1,ibmap_1,jbmap_1, &
							npx, npy, ngx, ngy, ngzm
	use messenger, only : icomm_grid
	use mpi
    implicit none

    real(kind(0.d0)),dimension(0:,0:,0:),intent(out)  :: uc,vc,wc 

	logical		  								      :: recv_flag, MD_BC_SLICE_average, exchangehalo
	integer											  :: icount, isource, idest, ierr
    logical, save 								      :: firsttime = .true.
	integer											  :: i,j,k,ii,ib,jb
	integer											  :: nclx,ncly,nclz,pcoords(3),extents(6)
	integer											  :: i1,i2,j1,j2,k1,k2
	integer											  :: jcmin_recv,jcmax_recv
    real(kind(0.d0)), allocatable, dimension(:,:,:,:) :: uvw_md, ucvcwc_md
    real(kind(0.d0)), allocatable, dimension(:,:,:,:) :: buf1, buf2
	real											  :: uvw_BC(4)

	integer		:: bufsize, jcmin_olap
	character	:: str_bufsize

	!Setup extents
	call CPL_get(jcmin_olap=jcmin_olap)
	jcmin_recv = jcmin_olap; jcmax_recv = jcmin_olap

	!Get number of CFD number of cells ready to receive data
	pcoords = (/ iblock,jblock,kblock /)
	call CPL_proc_extents(pcoords,CPL_realm(),extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1

	!Copy recieved data into appropriate cells
	MD_BC_SLICE_average = .false. 	!Average or cell by cell values
	if (MD_BC_SLICE_average) then

		! Use a single averaged value for the BC

		!Allocate array to CFD number of cells ready to receive data
		allocate(uvw_md(4,nclx,ncly,nclz)); uvw_md = VOID

		!Receive data
		call CPL_recv(uvw_md,jcmax_recv=jcmax_recv,jcmin_recv=jcmin_recv,recv_flag=recv_flag)
		if (any(uvw_md(:,:,jcmin_recv:jcmax_recv,:) .eq. VOID)) & 
			call error_abort("socket_coupler_get_md_BC error - VOID value copied to uc,vc or wc")

		!Average all cells on a CFD processor to give a single BC
		uvw_BC(1) = sum(uvw_md(1,:,jcmin_recv,:)) 
		uvw_BC(2) = sum(uvw_md(2,:,jcmin_recv,:))
		uvw_BC(3) = sum(uvw_md(3,:,jcmin_recv,:))
		uvw_BC(4) = sum(uvw_md(4,:,jcmin_recv,:)) 

		!Average in x so all processor have same global average BC
		call globalDirSum(uvw_BC,4,1)

		!Set full extent of halos to zero and set domain portion to MD values
		uc(:,:,0) = 0.d0; vc(:,:,1) = 0.d0; wc(:,:,0) = 0.d0
		uc(:,:,0) = uvw_BC(1)/uvw_BC(4)
		vc(:,:,1) = uvw_BC(2)/uvw_BC(4)
		wc(:,:,0) = uvw_BC(3)/uvw_BC(4)
	else

		! Use cell by cell value as BC

		! u interval [i1_u, i2_u], or [2,  ngx ] ??? 4 Procs = [2 32][3 34][3 34][3 35] ??? 30 31 32
		! v interval [i1_v, i2_v], or [1, ngx-1] ??? 4 Procs = [1 32][3 34][3 34][3 34] ???
		! w interval [i1_w, i2_w], or [1, ngx-1] ??? 4 Procs = [1 32][3 34][3 34][3 34] ???
		!uc(ngz  ,nlx+1,nly  )    vc(ngz  ,nlx  ,nly+1)   wc(ngz+1,nlx  ,nly  )

		!print'(a,28i3)', 'array extents',iblock,jblock,kblock,shape(uvw_md),nixb, niyb, nizb, & 
		!														   i1_u,i2_u,j1_u,j2_u, & 
		!														   i1_v,i2_v,j1_v,j2_v, & 
		!														   i1_w,i2_w,j1_w,j2_w, & 
		!														   i1_T,i2_T,j1_T,j2_T,k1_T,k2_T

		!Interpolate internal points and extrapolate halos OR use MPI to exchange halos?
		exchangehalo = .true. 
		if (exchangehalo) then

			!Allocate array to CFD number of cells ready to receive data including halos
			allocate(uvw_md(4,nclx+2,ncly,nclz+2)); uvw_md = VOID

			!Receive data from MD
			call CPL_recv(uvw_md(:,2:nclx+1,:,2:nclz+1),jcmax_recv=jcmax_recv,jcmin_recv=jcmin_recv,recv_flag=recv_flag)
			if (any(uvw_md(:,2:nclx+1,jcmin_recv:jcmax_recv,2:nclz+1) .eq. VOID)) & 
				call error_abort("socket_coupler_get_md_BC error - VOID value copied to uc,vc or wc")

			! - - - - Exchange cell halos in x direction - - - 
			!Determine size of send buffer and copy to buffer data to pass to lower neighbour
			allocate(buf1(4,1,ncly,nclz), buf2(4,1,ncly,nclz))
			icount = 4*ncly*nclz
			buf2 = 0.d0
			buf1(:,1,:,:) = uvw_md(:,2,:,2:nclz+1)

			! Send to lower neighbor
			call MPI_Cart_shift(icomm_grid, 0, -1, isource, idest, ierr)
			call MPI_sendrecv(buf1, icount, MPI_double_precision, idest, 0, &
			                  buf2, icount, MPI_double_precision, isource, 0, &
			                  icomm_grid, MPI_STATUS_IGNORE, ierr)

			!Copy to buffer data to pass to upper neighbour & save recieved data from upper neighbour 
			buf1(:,1,:,:)        = uvw_md(:,nclx+1,:,2:nclz+1)
			uvw_md(:,nclx+2,:,2:nclz+1) = buf2(:,1,:,:)

			! Send to upper neighbor
			call MPI_Cart_shift(icomm_grid, 0, +1, isource, idest, ierr)
			call MPI_sendrecv(buf1, icount, MPI_double_precision, idest, 0, &
			                  buf2, icount, MPI_double_precision, isource, 0, &
			                  icomm_grid, MPI_STATUS_IGNORE, ierr)

			!Save recieved data from lower neighbour
			uvw_md(:,1,:,2:nclz+1) = buf2(:,1,:,:)
			deallocate(buf1, buf2)
			! - - - END Exchange cell halos in x direction - - - 

			! - - - Exchange cell halos in z direction - - - 
			uvw_md(:,2:nclx+1,:,  1   ) = uvw_md(:,2:nclx+1,:,nclz+1)
			uvw_md(:,2:nclx+1,:,nclz+2) = uvw_md(:,2:nclx+1,:,  2   ) 

			!Allocate cell surface velocities
			allocate(ucvcwc_md(3,nclx+1,ncly,nclz+1))
			ucvcwc_md = 0.d0

			! - - - Interpolate cell surfaces including halos- - -
			! X direction
			do i=1, nclx+1
			do j=jcmin_recv,jcmax_recv
			do k=2, nclz+1
				ucvcwc_md(1,i,j,k-1) = 0.5d0*(uvw_md(1,i,j,k)/uvw_md(4,i,j,k) + uvw_md(1,i+1,j,k)/uvw_md(4,i+1,j,k))
			enddo
			enddo
			enddo
			! Y direction
			do i=2, nclx+1
			do j=jcmin_recv,jcmax_recv
			do k=2, nclz+1
				ucvcwc_md(2,i-1,j,k-1) =      uvw_md(2,i,j,k)/uvw_md(4,i,j,k)
			enddo
			enddo
			enddo
			! Z direction
			do i=2, nclx+1
			do j=jcmin_recv,jcmax_recv
			do k=1, nclz+1
				ucvcwc_md(3,i-1,j,k) = 0.5d0*(uvw_md(3,i,j,k)/uvw_md(4,i,j,k) + uvw_md(3,i,j,k+1)/uvw_md(4,i,j,k+1))
			enddo
			enddo
			enddo

		else

			!Allocate array to CFD number of cells ready to receive data including halos
			allocate(uvw_md(4,nclx,ncly,nclz)); uvw_md = VOID

			!Receive data from MD
			call CPL_recv(uvw_md,jcmax_recv=jcmax_recv,jcmin_recv=jcmin_recv,recv_flag=recv_flag)
			if (any(uvw_md(:,:,jcmin_recv:jcmax_recv,:) .eq. VOID)) & 
				call error_abort("socket_coupler_get_md_BC error - VOID value copied to uc,vc or wc")

			!Allocate cell surface velocities
			allocate(ucvcwc_md(3,nclx+1,ncly,nclz+1))
			ucvcwc_md = 0.d0

			! - - - Interpolate cell surfaces - - -
			! X direction
			do i=1, nclx-1
			do j=jcmin_recv,jcmax_recv
			do k=1, nclz
				ucvcwc_md(1,i+1,j,k) = 0.5d0*(uvw_md(1,i,j,k)/uvw_md(4,i,j,k) + uvw_md(1,i+1,j,k)/uvw_md(4,i+1,j,k))
			enddo
			enddo
			enddo
			! Y direction
			do i=1, nclx
			do j=jcmin_recv,jcmax_recv
			do k=1, nclz
				ucvcwc_md(2,i,j,k) =          uvw_md(2,i,j,k)/uvw_md(4,i,j,k)
			enddo
			enddo
			enddo
			! Z direction
			do i=1, nclx
			do j=jcmin_recv,jcmax_recv
			do k=1, nclz-1
				ucvcwc_md(3,i,j,k+1) = 0.5d0*(uvw_md(3,i,j,k)/uvw_md(4,i,j,k) + uvw_md(3,i,j,k+1)/uvw_md(4,i,j,k+1))
			enddo
			enddo
			enddo

			! - - - Extrapolate bottom surface - - -
			! X direction
			ucvcwc_md(1,  1   ,jcmin_recv,1:nclz)= 1.5d0*uvw_md(1,1,jcmin_recv,1:nclz)/uvw_md(4,1,jcmin_recv,1:nclz) & 
												  -0.5d0*uvw_md(1,2,jcmin_recv,1:nclz)/uvw_md(4,2,jcmin_recv,1:nclz)
			! Z direction
			ucvcwc_md(3,1:nclx,jcmin_recv,  1   )= 1.5d0*uvw_md(3,1:nclx,jcmin_recv,1)/uvw_md(4,1:nclx,jcmin_recv,1) & 
												  -0.5d0*uvw_md(3,1:nclx,jcmin_recv,2)/uvw_md(4,1:nclx,jcmin_recv,2)
			!Extrapolate top surface
			! X direction
			ucvcwc_md(1,nclx+1,jcmin_recv,1:nclz)= 1.5d0*uvw_md(1,nclx  ,jcmin_recv,1:nclz)/uvw_md(4,nclx  ,jcmin_recv,1:nclz) & 
												  -0.5d0*uvw_md(1,nclx-1,jcmin_recv,1:nclz)/uvw_md(4,nclx-1,jcmin_recv,1:nclz)
			! Z direction
			ucvcwc_md(3,1:nclx,jcmin_recv,nclz+1)= 1.5d0*uvw_md(3,1:nclx,jcmin_recv,nclz  )/uvw_md(4,1:nclx,jcmin_recv,nclz-1) & 
												  -0.5d0*uvw_md(3,1:nclx,jcmin_recv,nclz-1)/uvw_md(4,1:nclx,jcmin_recv,nclz-1)

		endif

		! ************** DEBUG ***************
		!do i=1, nclx+1
		!	if (i .le. nclx ) then
		!		print'(a,4i6,3f14.5)','qqqq',iblock,jblock,kblock,i, ucvcwc_md(1,i,jcmin_recv,4), & 
		!																uvw_md(1,i,jcmin_recv,4)/uvw_md(4,i,jcmin_recv,4), &
		!																uvw_md(1,i,jcmin_recv,4)
		!	else
		!		print'(a,4i6,3f14.5)','qqqq',iblock,jblock,kblock,i, ucvcwc_md(1,i,jcmin_recv,4), 0.d0, 0.d0
		!	endif
		!enddo

		!do k=1, nclz+2
		!	if (k .le. nclz+1 ) then
		!		print'(a,4i6,3f14.5)','qqqq',iblock,jblock,kblock,k, ucvcwc_md(3,4,jcmin_recv,k), & 
		!																uvw_md(3,4,jcmin_recv,k)/uvw_md(4,4,jcmin_recv,k), &
		!																uvw_md(3,4,jcmin_recv,k)
		!	else
		!		print'(a,4i6,3f14.5)','qqqq',iblock,jblock,kblock,k, 0.d0, & 
		!																uvw_md(3,4,jcmin_recv,k)/uvw_md(4,4,jcmin_recv,k), &
		!																uvw_md(3,4,jcmin_recv,k)
		!	endif
		!enddo
		!print'(a,13i7)', 'cells per proc', iblock,jblock,kblock,nixb, nizb,nlx,ngz, i2_u,i1_u, i2_w,i1_w, i1_T,i2_T
		! ************** DEBUG ***************


		!Set CFD halo - internal cell surfaces from 2 to nclx
		! X direction
		uc(:,:,0) = 0.d0; 
		do i=1, nclx +1
			if (iblock .eq. 1) then
				ii = i + i1_u - 2
			else
				ii = i + i1_u - 1
			endif
			if (ii .gt. i2_u) cycle
		do j=jcmin_recv,jcmax_recv
		do k=1, nclz
			uc(k,ii,j-1) = ucvcwc_md(1,i,j,k)
			!print'(a,11i5,3f14.6)', 'uc', iblock,jblock,kblock,i1_u,i2_u,i,j,k,ii,imap_1(ii),ibmap_1(ii), & 
			!				      ucvcwc_md(1,i,j,k),uc(k,ii,j-1),sin(2.d0*3.141592654d0*(xg(ibmap_1(ii),1)-0.5d0*xL_cfd)/xL_cfd)!!,uvw_md(1,i,j,k)!,vc(k,ii,j ),wc(k,ii,j-1)
		enddo
		enddo
		enddo
		! Y direction
		vc(:,:,1) = 0.d0; vc(:,:,0) = 0.d0
		do i=1, nclx
			ii = i + i1_v - 1
			if (ii .gt. i2_v) cycle
		do j=jcmin_recv,jcmax_recv
		do k=1, nclz
			vc(k,ii,j  ) = ucvcwc_md(2,i,j,k)
			!print'(a,11i5,3f14.6)', 'vc', iblock,jblock,kblock,i1_v,i2_v,i,j,k,ii,imap_1(ii),ibmap_1(ii), & 
			!				      ucvcwc_md(2,i,j,k),vc(k,ii,j  ),yg(1,j)-0.5d0*yL_cfd
		enddo
		enddo
		enddo
		! Z direction
		wc(:,:,0) = 0.d0
		do i=1, nclx
			ii = i + i1_w - 1
			if (ii .gt. i2_w) cycle
		do j=jcmin_recv,jcmax_recv
		do k=1, nclz + 1
			wc(k,ii,j-1) = ucvcwc_md(3,i,j,k)
			!print'(a,11i5,3f14.6)', 'wc', iblock,jblock,kblock,i1_w,i2_w,i,j,k,ii,imap_1(ii),ibmap_1(ii), & 
			!				      ucvcwc_md(3,i,j,k),wc(k,ii,j-1),sin(2.d0*3.141592654d0*(zg(k)-0.5d0*zL_cfd)/zL_cfd)
		enddo
		enddo
		enddo

	endif
	
end subroutine socket_coupler_get_md_BC

!=============================================================================
!   ___  _____  _  _  ___  ____  ____    __    ____  _  _  ____ 
!  / __)(  _  )( \( )/ __)(_  _)(  _ \  /__\  (_  _)( \( )(_  _)
! ( (__  )(_)(  )  ( \__ \  )(   )   / /(__)\  _)(_  )  (   )(  
!  \___)(_____)(_)\_)(___/ (__) (_)\_)(__)(__)(____)(_)\_) (__) 
!

subroutine  socket_coupler_send
	use CPL, only : error_abort, CPL_get
	implicit none

	integer :: constraint_algorithm
	integer :: OT, NCER, Flekkoy, off

	call CPL_get(	constraint_algo	      = constraint_algorithm, & 
					constraint_OT         = OT,        & 
					constraint_NCER       = NCER,      &
					constraint_Flekkoy    = Flekkoy,   &
					constraint_off        = off          )

	if ( constraint_algorithm .eq. off ) then
		return
	else if ( constraint_algorithm .eq. OT ) then
		call error_abort("OT constraint force not yet implemented")
	else if ( constraint_algorithm .eq. NCER ) then
		call socket_coupler_send_mass
		call socket_coupler_send_velocity
	else if ( constraint_algorithm .eq. Flekkoy ) then
		call socket_coupler_send_stress
	else
		call error_abort("Unrecognised constraint algorithm flag")
	end if	

end subroutine socket_coupler_send

!---------------------------------------------------------------------
! Send continuum velocity in x direction to be used by MD

subroutine socket_coupler_send_mass
    use CPL, only : CPL_send,CPL_olap_extents,CPL_overlap,CPL_get,CPL_realm
 	use data_export, only : uc,i1_v,ngz,vc,nlx,nlxb,ibmin_1,ibmax_1,iblock,jblock,kblock
    implicit none

	logical	:: send_flag
    integer	:: i,j,ii,ixyz,icell,jcell,kcell,npercell,nclx,ncly,nclz
    integer	:: coord(3),extents(6),cnstd(6)
    real(kind(0.d0)),dimension(:,:,:,:), allocatable :: sendbuf

	! Check processor is inside MD/CFD overlap zone 
	if (.not.(CPL_overlap())) return

	!Number of cells to package and send
	call CPL_get( icmin_cnst=cnstd(1),icmax_cnst=cnstd(2), &
	              jcmin_cnst=cnstd(3),jcmax_cnst=cnstd(4), &
	              kcmin_cnst=cnstd(5),kcmax_cnst=cnstd(6)  )

	!One mass component
	npercell = 1

	!Allocate array for size of data on local processor
	coord = (/iblock,jblock,kblock /)
	call CPL_olap_extents(coord,CPL_realm(),extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1
	allocate(sendbuf(npercell,nclx,ncly,nclz))

	!Interpolate cell centres using surfaces
	sendbuf(:,:,:,:) = 0.d0
	do j=cnstd(4),cnstd(4)
	do i=1,nclx
		ii = i + i1_v - 1
		sendbuf(1,i,j,:) = vc(:,ii,j)
	enddo
	!call printf(sendbuf(1,40,j,:))
	enddo

	call CPL_send( sendbuf,                                 &
	               icmin_send=cnstd(1),icmax_send=cnstd(2), &
	               jcmin_send=cnstd(4),jcmax_send=cnstd(4), &
	               kcmin_send=cnstd(5),kcmax_send=cnstd(6), &
	               send_flag=send_flag                        )

end subroutine socket_coupler_send_mass


!---------------------------------------------------------------------
! Send continuum velocity in x direction to be used by MD

subroutine socket_coupler_send_velocity
    use CPL, only : CPL_send,CPL_olap_extents,CPL_overlap,CPL_get,CPL_realm
 	use data_export, only : uc,i1_u,i2_u,ngz,nlx,nlxb,ibmin_1,ibmax_1,iblock,jblock,kblock
    implicit none

	logical	:: send_flag
    integer	:: i,j,ii,ixyz,icell,jcell,kcell,npercell,nclx,ncly,nclz
    integer	:: coord(3),extents(6),cnstd(6)
    real(kind(0.d0)),dimension(:,:,:,:), allocatable :: sendbuf

	! Check processor is inside MD/CFD overlap zone 
	if (.not.(CPL_overlap())) return

	!Number of cells to package and send
	call CPL_get( icmin_cnst=cnstd(1),icmax_cnst=cnstd(2), &
	              jcmin_cnst=cnstd(3),jcmax_cnst=cnstd(4), &
	              kcmin_cnst=cnstd(5),kcmax_cnst=cnstd(6)  )

	!Three velocity components
	npercell = 3

	!Allocate array for size of data on local processor
	coord = (/iblock,jblock,kblock /)
	call CPL_olap_extents(coord,CPL_realm(),extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1
	allocate(sendbuf(npercell,nclx,ncly,nclz))

	!Interpolate cell centres using surfaces
	sendbuf(:,:,:,:) = 0.d0
	do j=cnstd(3),cnstd(4)
	do i=1,nclx
		ii = i + i1_u - 1
		sendbuf(1,i,j,:) = 0.5d0*(uc(:,ii,j) + uc(:,ii+1,j))
	enddo
	!call printf(sendbuf(1,40,j,:))
	enddo

	call CPL_send( sendbuf,                                 &
	               icmin_send=cnstd(1),icmax_send=cnstd(2), &
	               jcmin_send=cnstd(3),jcmax_send=cnstd(4), &
	               kcmin_send=cnstd(5),kcmax_send=cnstd(6), &
	               send_flag=send_flag                        )

end subroutine socket_coupler_send_velocity


!---------------------------------------------------------------------
! Send continuum Stress in all directions to be used by MD

subroutine socket_coupler_send_stress
    use CPL, only : CPL_send,CPL_olap_extents,CPL_overlap,CPL_get,CPL_realm
 	use data_export, only : uc,vc,wc,P,i1_T,i2_T,ngz,nlx,nlxb, &
							ibmin_1,ibmax_1,iblock,jblock,kblock,xpg,ibmap_1
    implicit none

	logical	:: send_flag
    integer	:: i,j,k,ii,ixyz,icell,jcell,kcell,npercell,nclx,ncly,nclz
    integer	:: coord(3),extents(6),cnstd(6)
    real(kind(0.d0)),dimension(:,:,:,:), allocatable 	:: sendbuf
	real(kind(0.d0)),dimension(:,:,:,:,:), allocatable 	:: stress

	! Check processor is inside MD/CFD overlap zone 
	if (.not.(CPL_overlap())) return
	!Number of cells to package and send
	call CPL_get( icmin_cnst=cnstd(1),icmax_cnst=cnstd(2), &
	              jcmin_cnst=cnstd(3),jcmax_cnst=cnstd(4), &
	              kcmin_cnst=cnstd(5),kcmax_cnst=cnstd(6)  )

	!Three by three stress tensor components
	npercell = 9

	!Allocate array for size of data on local processor
	coord = (/iblock,jblock,kblock /)
	call CPL_olap_extents(coord,CPL_realm(),extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1
	allocate(sendbuf(npercell,nclx,ncly,nclz))

	!Get stress tensor at cell centers and store in buffer
	call Evaluate_stress(uc,vc,wc,P,stress)

	!DEBUG DEBUG - set passed values to location in space
	!do i=1,nclx
	!	if (iblock .eq. 1) then
	!		ii = i + i1_u - 2
	!	else
	!		ii = i + i1_u - 1
	!	endif
	!	if (ii .gt. i2_u) cycle
	!	ii = i + i1_T - 1
	!	if (ii .gt. i2_T) cycle
	!	stress(1,1,:,ii,:) = xpg(ibmap_1(ii),1)
	!enddo
	!DEBUG DEBUG - set passed values to location in space

	sendbuf = 0.d0
	do i=1,nclx
	do j=cnstd(3),cnstd(4)
	do k=1,nclz
		!if (iblock .eq. 1) then
		!	ii = i + i1_u - 2
		!else
		!	ii = i + i1_u - 1
		!endif
		ii = i + i1_T - 1
		if (ii .gt. i2_T) cycle
		sendbuf(:,i,j,k) = reshape(stress(:,:,k,ii,j),(/ npercell /))
		!print'(a,14i5,9f16.8)','sending', iblock,jblock,kblock,i,j,k,ii,cnstd,ibmap_1(ii),sendbuf(1,i,j,k)
	enddo
	enddo
	enddo

	!stress(1:3,ibin(1),ibin(2)-minbin(2)+2,ibin(3)) = sin(2*pi*(globalise(ibin(:)*dxyz(:)-halfdomain(:)-0.5*dxyz(:)))/globaldomain(:))

	!print*, 'sent stress cfd',iblock,jblock,kblock,nclz,size(sendbuf(:,i-1,j-1,k-1))
	!do k=1,nclz
	!    call printf(stress(2,2,k,3,:))
	!enddo

	!Send stress tensor to MD code
	call CPL_send( sendbuf,                                 &
	               icmin_send=cnstd(1),icmax_send=cnstd(2), &
	               jcmin_send=cnstd(3),jcmax_send=cnstd(4), &
	               kcmin_send=cnstd(5),kcmax_send=cnstd(6), &
	               send_flag=send_flag                        )

end subroutine socket_coupler_send_stress

!---------------------------------------------------------------------
! 	Gets cell centered stress from velocity and pressure field


subroutine Evaluate_stress(uc,vc,wc,P,stress)
	use data_export, only : ibmin,ibmax,jbmin,jbmax,kbmin,kbmax, &
						   	imap_1,jmap_1,ibmap_1,jbmap_1,npx,npy, & 
							ngx,ngy,ngzm,iblock,jblock,kblock,vp, &
							suxix,suxiy,svetax,svetay,spz,visc,nlx,nly,ngz
	implicit none

	real(kind(0.d0)),intent(in),dimension(:,:,:)					:: uc,vc,wc,P
	real(kind(0.d0)),intent(out),dimension(:,:,:,:,:),allocatable	:: stress

	integer												:: i,j,k,ib,jb,i1,i2,j1,j2
	real(kind(0.d0))									:: P_guage
	real(kind(0.d0)),dimension(:,:,:,:,:),allocatable	:: dUidxj

	call Evaluate_strain(uc,vc,wc,dUidxj)

	i1 = ibmin ; if (iblock ==  1 ) i1 =  1   ; i1 = imap_1(i1)
	i2 = ibmax ; if (iblock == npx) i2 = ngx-1; i2 = imap_1(i2)
	j1 = jbmin ; if (jblock ==  1 ) j1 =  1   ; j1 = jmap_1(j1)
	j2 = jbmax ; if (jblock == npy) j2 = ngy-1; j2 = jmap_1(j2)

	!Allocate array to store stress tensor
	P_guage = 4.d0
	if (allocated(stress)) deallocate(stress)
	allocate(stress(3,3,0:ngz,0:nlx,0:nly)); stress = 0.d0

	do j=j1,j2
	do i=i1,i2
	do k=1,ngzm
		stress(:,:,k,i,j) = -(P(k,i,j) + P_guage) * IDM(3) & 
							-(2.d0/3.d0)*visc*trace(dUidxj(:,:,k,i,j)) * IDM(3) &
						    +visc*(dUidxj(:,:,k,i,j)+transpose(dUidxj(:,:,k,i,j)))

	enddo
	enddo
	enddo

end subroutine Evaluate_stress


!---------------------------------------------------------------------
! 			Gets cell centered strain from velocity field

subroutine Evaluate_strain(uc,vc,wc,dUidxj)
	use data_export, only : ibmin,ibmax,jbmin,jbmax,kbmin,kbmax, &
						   	imap_1,jmap_1,ibmap_1,jbmap_1,npx,npy, & 
							ngx,ngy,ngzm,iblock,jblock,kblock,vp, &
							suxix,suxiy,svetax,svetay,spz,nlx,nly,ngz
	implicit none

	real(kind(0.d0)),intent(in),dimension(0:,0:,0:)					:: uc,vc,wc
	real(kind(0.d0)),intent(out),dimension(:,:,:,:,:),allocatable	:: dUidxj

	integer			:: i,j,k,ib,jb,i1,i2,j1,j2
	real(kind(0.d0)):: voli, ds
	real(kind(0.d0)):: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

	! Average in the homogeneous z direction
	ds = 1./ngzm

	!-----------------------------------------------------------
	! Derivatives: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
	!-----------------------------------------------------------
	! Dissipation:  epsilon=\nu * average ( duidxj * duidxj )
	!-----------------------------------------------------------

	i1 = ibmin ; if (iblock ==  1 ) i1 =  1   ; i1 = imap_1(i1)
	i2 = ibmax ; if (iblock == npx) i2 = ngx-1; i2 = imap_1(i2)
	j1 = jbmin ; if (jblock ==  1 ) j1 =  1   ; j1 = jmap_1(j1)
	j2 = jbmax ; if (jblock == npy) j2 = ngy-1; j2 = jmap_1(j2)

	!Allocate array to store strain tensor
	if (allocated(dUidxj)) deallocate(dUidxj)
	allocate(dUidxj(3,3,0:ngz,0:nlx,0:nly)); dUidxj = 0.d0

	do j=j1,j2
		jb = jbmap_1(j)
	do i=i1,i2
		ib = ibmap_1(i)
		voli = 1./vp(ib,jb)
	do k=1,ngzm
		!------------------------ du/dx -------------------------
		dudx= voli  * (  uc(k,i+1,j)*suxix(ib+1,jb) &
						-uc(k, i ,j)*suxix(ib  ,jb) &
				+0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j+1)+uc(k,i+1,j+1))*svetax(ib,jb+1) &
				-0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j-1)+uc(k,i+1,j-1))*svetax(ib,jb  )   )

		dudy= voli  * (  uc(k,i+1,j)*suxiy(ib+1,jb) &
						-uc(k, i ,j)*suxiy( ib ,jb) &
				+0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j+1)+uc(k,i+1,j+1))*svetay(ib,jb+1) &
				-0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j-1)+uc(k,i+1,j-1))*svetay(ib, jb )   )

		dudz= voli  *	spz(ib,jb)*( 0.25*(uc( k ,i,j)+uc( k ,i+1,j)+uc(k+1,i,j)+uc(k+1,i+1,j)) &
					    -0.25*(uc(k-1,i,j)+uc(k-1,i+1,j)+uc( k ,i,j)+uc( k ,i+1,j))  )

		!------------------------ dv/dx -------------------------
		dvdx= voli  * (  vc(k,i,j+1)*svetax(ib,jb+1) &
						-vc(k,i, j )*svetax(ib, jb ) &
				+0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i+1,j)+vc(k,i+1,j+1))*suxix(ib+1,jb) &
				-0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i-1,j)+vc(k,i-1,j+1))*suxix( ib ,jb)   )

		dvdy= voli  * (  vc(k,i,j+1)*svetay(ib,jb+1) &
						-vc(k,i, j )*svetay(ib, jb ) &
				+0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i+1,j)+vc(k,i+1,j+1))*suxiy(ib+1,jb) &
				-0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i-1,j)+vc(k,i-1,j+1))*suxiy( ib ,jb)   )

		dvdz= voli  *	spz(ib,jb)*( 0.25*(vc( k ,i,j)+vc( k ,i,j+1)+vc(k+1,i,j)+vc(k+1,i,j+1)) &
					    			-0.25*(vc(k-1,i,j)+vc(k-1,i,j+1)+vc( k ,i,j)+vc( k ,i,j+1))   )

		!------------------------ dw/dx -------------------------
		dwdx= voli  *  ( 0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i+1,j  )+wc(k+1,i+1,j  ))*suxix(ib+1,jb) &
						-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i-1,j  )+wc(k+1,i-1,j  ))*suxix( ib ,jb) &
						+0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j+1)+wc(k+1,i  ,j+1))*svetax(ib,jb+1) &
						-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j-1)+wc(k+1,i  ,j-1))*svetax(ib, jb )   )

		dwdy= voli  *  ( 0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i+1,j  )+wc(k+1,i+1,j  ))*suxiy(ib+1,jb) &
						-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i-1,j  )+wc(k+1,i-1,j  ))*suxiy( ib ,jb) &
						+0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j+1)+wc(k+1,i  ,j+1))*svetay(ib,jb+1) &
						-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j-1)+wc(k+1,i  ,j-1))*svetay(ib, jb )   )

		dwdz= voli  *	spz(ib,jb)*(wc(k+1,i,j)-wc(k,i,j))

		!--------------------- Store strain rate tensor ----------------------
		dUidxj(1,:,k,i,j) = (/dudx,dudy,dudz/)
		dUidxj(2,:,k,i,j) = (/dvdx,dvdy,dvdz/)
		dUidxj(3,:,k,i,j) = (/dwdx,dwdy,dwdz/)

	end do
	end do
	end do

end subroutine Evaluate_strain


!  A function that returns the trace of a square matrix

function trace(A)
   implicit none

	real(kind(0.d0)) 					:: trace
	real(kind(0.d0)) , dimension(:,:) 	:: a

	integer :: j            ! local loop index
   
	trace = 0.d0   ! Set trace to zero
	if (size(A,1) .ne. size(A,2)) stop "Trace error - must be square matrix"
	do j = 1,size(A,1)                      
	    trace = trace + A(j,j)  ! Add along diagonal 
	enddo

end function trace

!  A function that produces an Identity Matrix of dimension (n,n)

function IDM(n)
	implicit none
	integer, intent(in) 				:: n
	real(kind(0.d0)) , dimension(n,n)	:: IDM

	integer :: j            		! local loop index

	IDM = 0.d0   			! Set all elements to zero
	do j = 1,n
	    IDM(j,j) = 1.d0  	! Change value of diagonal elements to one
	enddo

end function IDM



!▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
!  ۩ ۩ ۩  ۩ ۩ ۩  ۩ ۩ ۩  ۩ ۩ ۩  ۩ ۩ ۩  ۩ ۩ ۩  ۩ ۩ ۩  ۩ ۩ ۩  ۩ ۩ ۩  ۩ ۩ ۩  ۩ ۩ 
!▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

! ▬▬▬▬▬▬▬▬(●̪•) ▬▬▬▬▬▬(●̪•)▬▬▬▬▬▬(-_-)▬▬▬▬▬▬▬ஜ۩۞۩ஜ▬▬▬▬▬(●̪•)▬▬▬▬▬▬(-_-)▬▬▬▬▬(-_-) ▬▬▬▬▬
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

	call CPL_Cart_coords(CPL_WORLD_COMM,rank_world,realm,3,coord,ierr)
	!print'(2a,5i8)', 'CFD SIDE',realm_name(realm), rank_world, CPL_overlap,coord

	if (.not.(CPL_overlap())) return

	! Test Sending from MD to CFD							   
	if (realm .eq. md_realm) then	



		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_olap_extents(coord,realm,extents)

		allocate(sendbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))

		print'(2a,11i7)', 'sent size',realm_name(realm),extents,size(sendbuf),shape(sendbuf)

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

		print'(2a,11i7)', 'recv size', realm_name(realm),extents,size(recvbuf),shape(recvbuf)
		recvbuf = -444.d0
		call CPL_recv(recvbuf,jcmax_recv=jcmax_recv,jcmin_recv=jcmin_recv,recv_flag=recv_flag)

		if (recv_flag .eqv. .true.) then
			do kcell=extents(5),extents(6)
			do jcell=jcmin_recv,jcmax_recv  !extents(3),extents(4)
			do icell=extents(1),extents(2)
			do ixyz =1,npercell
					write(5000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
					      'recv CFD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
					       recvbuf(ixyz,icell,jcell,kcell)
			end do
			end do
			end do
			end do
		endif
	end if								   
	
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
	if (.not.(CPL_overlap())) return

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
				write(11000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
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

 	if (.not.(CPL_overlap())) return

	!print*, 'test_gather_scatter called on CFD proc ID:', rank_realm, rank_world

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
	gatherlims  = (/1,85,15,21, 3, 4/)
	scatterlims = (/1,85, 2, 9, 1, 8/)
	call CPL_gather(u,3,gatherlims,gatheru)
	call CPL_scatter(stress,9,scatterlims,scatterstress)

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

	!print*, 'test_gather_scatter finished on CFD proc ID:', rank_realm, rank_world
	
end subroutine test_gather_scatter

! ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬

#if COUPLER_DEBUG_LA
!---------------------------------------------------------------------
! debug subroutine that writes uc for gnuplot
    subroutine socket_coupler_write_uc(ntime)
        use data, only : uc,ngz,i1_u,i2_u, j1_u,j2_u
        use data_export, only : kmin, kmax
        use messenger, only : icomm_grid, icoord
        use computation_parameters, only : prefix_dir
		use coupler_module, only : rank_realm
        implicit none

        integer, intent(in) :: ntime

        integer, parameter :: psave = 10 
        integer, save      :: icount=0

        integer i,j,k,is, ie, ks, ke, ierr
        logical, save :: firsttime=.true.
        character(len=32) file_position

        icount = icount + 1
        if ( mod(icount,psave) /= 0) return

        ! pick the z index range of the slab to write
        ! for beginig we look at the midle layer
        ! this should be fixed with coupler input 
        ! parameter
        ks =kmin + (kmax-kmin)/2-1; ke = ks

        if (firsttime) then
            firsttime = .false.
            file_position = "rewind"
        else
            file_position = "append"
        endif

        open (unit=1003, file=trim(prefix_dir)//"results/continuum_uc.txt",position=file_position)
        
        is = i1_u
        ie = i2_u


        if (icoord(1,rank_realm)==1) then
            is = i1_u-1
        endif

        write(1003,'(a,5I6)')'# step', ntime, is, ie, j1_u, j2_u
        do k= ks,ke
            do i=is,ie
                do j=j1_u-1,j2_u+1
                    write(1003,'(1000E12.4)') uc(k,i,j)
                enddo
                 write(1003,'(1x)')
             enddo
             write(1003,'(1x/1x)')
         enddo
        close(1003)
        
  end subroutine socket_coupler_write_uc
#endif

#endif
end module continuum_coupler_socket
