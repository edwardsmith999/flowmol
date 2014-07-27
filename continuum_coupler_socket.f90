module continuum_coupler_socket
        implicit none
#if USE_COUPLER
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
	implicit none

	call CPL_create_comm(cfd_realm,CFD_COMM,ierr)
	prefix_dir ="./couette_serial/"

end subroutine socket_coupler_invoke

    
    
subroutine socket_coupler_init
    use computational_constants, only : nx, ny, & 
					nsteps => continuum_Nsteps,&
        			continuum_delta_t, lx, ly, lz
    use physical_constants, only : rho
    use grid_arrays, only : mx, my
    use continuum_data_export, only : npx,npy,npz,icoord
    use CPL, only : coupler_cfd_init, CPL_write_header, VOID
	use messenger, only : icomm_grid
    implicit none

	integer							:: i,j
	integer,dimension(1)			:: iTmin,iTmax,jTmin,jTmax,kTmin,kTmax
    integer,dimension(3)		   	:: ijkmin,ijkmax,npxyz,ngxyz
    real(kind(0.d0)),dimension(3)	:: xyzL
    ! 2D problem, z direction parameters are set to trivial values
    real(kind(0.d0)),dimension(:,:),allocatable :: x0,y0
    real(kind(0.d0)) :: z0(1)
          
	!Define compound arrays to make passing more concise
	ijkmin = (/ 1, 1, 1 /)
	!Minus 1 as coupler requires max cells NOT max cell vertices
	ijkmax = (/ nx+2, ny+2, 1 /)
	iTmin = ijkmin(1); iTmax = ijkmax(1)
	jTmin = ijkmin(2); jTmax = ijkmax(2)
	kTmin = ijkmin(3); kTmax = ijkmax(3)
	npxyz  = (/ npx , npy , npz  /)
	!Minus 1 as coupler requires no. of cells NOT no. cell vertices
	ngxyz  = (/ nx , ny , 1   /) 
	xyzL   = (/ lx , ly , lz  /)
    z0 = (/0.d0/)

	allocate(x0(nx+1,ny+1))
	allocate(y0(nx+1,ny+1))
	do i =1,nx+1
	do j =1,ny+1
		x0(i,j) = mx(i)
		y0(i,j) = my(j)
       ! print'(a,2i5,2f10.5)', 'continuum socket_coupler_init grid', i,j,x0(i,j),y0(i,j)
	enddo
	enddo

    ! nsteps = nsteps+1 for the intialisation step in setup_continuum
    call coupler_cfd_init(nsteps=nsteps+1,& 
						  dt=continuum_delta_t, & 
						  icomm_grid=icomm_grid, & 
						  icoord=icoord, & 
						  npxyz_cfd=npxyz, & 
						  xyzL=xyzL, & 
						  ncxyz=ngxyz, & 
						  density=rho, & 
						  ijkcmax=ijkmax, & 
						  ijkcmin=ijkmin, & 
						  iTmin=iTmin, & 
						  iTmax=iTmax, & 
						  jTmin=jTmin, & 
						  jTmax=jTmax, & 
						  kTmin=kTmin, & 
						  kTmax=kTmax, & 
						  xgrid=x0, & 
						  ygrid=y0, & 
						  zgrid=z0)

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

subroutine socket_coupler_send_CFD_to_MD
    use grid_arrays, only : uc,vc,tau_xx,tau_xy,tau_yx,tau_yy
    use CPL, only : CPL_get, error_abort
    implicit none

	integer :: constraint_algorithm
	integer :: OT, NCER, Flekkoy, CV, off

	call CPL_get(	constraint_algo	      = constraint_algorithm, & 
					constraint_OT         = OT,        & 
					constraint_NCER       = NCER,      &
					constraint_Flekkoy    = Flekkoy,   &
                    constraint_CV         = CV,        &
					constraint_off        = off          )
	
	if ( constraint_algorithm .eq. off ) then
		return
	else if ( constraint_algorithm .eq. OT ) then
		call error_abort("OT constraint force not yet implemented")
	else if ( constraint_algorithm .eq. NCER ) then
		call socket_coupler_send_velocity(uc,vc)
	else if ( constraint_algorithm .eq. Flekkoy ) then
		call socket_coupler_send_stress(tau_xx,tau_xy,tau_yx,tau_yy)
	else if ( constraint_algorithm .eq. CV ) then
		call socket_coupler_send_velocity(uc,vc)
		call socket_coupler_send_stress(tau_xx,tau_xy,tau_yx,tau_yy)
	else
		call error_abort("Unrecognised constraint algorithm flag")
	end if	

end subroutine socket_coupler_send_CFD_to_MD

!---------------------------------------------------------------------
! Send continuum velocity in x direction to be used by MD

subroutine socket_coupler_send_velocity(u,v)
    use CPL, only : CPL_send,CPL_olap_extents,CPL_overlap,CPL_get,CPL_realm
    implicit none

    real(kind(0.d0)), intent(in) :: u(:,:), v(:,:)

	logical	:: send_flag
    integer	:: i,j,k,jj,ixyz,icell,jcell,kcell,npercell,nclx,ncly,nclz
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
	coord = (/1,1,1/)
	call CPL_olap_extents(coord,CPL_realm(),extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1
	allocate(sendbuf(npercell,nclx,cnstd(3):cnstd(4),nclz))

	!Interpolate cell centres using surfaces
    !print*, 'output', shape(sendbuf),cnstd
	sendbuf(:,:,:,:) = 0.d0
	do i=1,nclx
	do j=cnstd(3),cnstd(4)
	do k=1,nclz
		sendbuf(1,i,j,k) = u(i,j)
		sendbuf(2,i,j,k) = v(i,j)
		sendbuf(3,i,j,k) = 0.d0
	enddo
	enddo
	enddo

	call CPL_send( sendbuf,                                 &
	               icmin_send=cnstd(1),icmax_send=cnstd(2), &
	               jcmin_send=cnstd(3)-cnstd(3)+1,jcmax_send=cnstd(4)-cnstd(3)+1, &
	               kcmin_send=cnstd(5),kcmax_send=cnstd(6), &
	               send_flag=send_flag                        )

end subroutine socket_coupler_send_velocity


!---------------------------------------------------------------------
! Send continuum Stress in all directions to be used by MD

subroutine socket_coupler_send_stress(tau_xx,tau_xy,tau_yx,tau_yy)
    use CPL, only : CPL_send,CPL_olap_extents,CPL_overlap,CPL_get,CPL_realm
    implicit none

    real(kind(0.d0)), intent(in) :: tau_xx(:,:,:), tau_xy(:,:,:), &
									tau_yx(:,:,:), tau_yy(:,:,:)

	logical	:: send_flag
    integer	:: i,j,k,jj,ixyz,icell,jcell,kcell,npercell,nclx,ncly,nclz
    integer	:: coord(3),extents(6),cnstd(6)
    real(kind(0.d0)),dimension(:,:,:,:), allocatable 	:: sendbuf
	real(kind(0.d0)),dimension(:,:,:,:,:), allocatable 	:: stress

	! Check processor is inside MD/CFD overlap zone 
	if (.not.(CPL_overlap())) return
	!Number of cells to package and send
	call CPL_get( icmin_cnst=cnstd(1),icmax_cnst=cnstd(2), &
	              jcmin_cnst=cnstd(3),jcmax_cnst=cnstd(4), &
	              kcmin_cnst=cnstd(5),kcmax_cnst=cnstd(6)  )

	!Allocate array for size of data on local processor
	npercell = 18
	coord = (/1,1,1 /)
	call CPL_olap_extents(coord,CPL_realm(),extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1
	allocate(sendbuf(npercell,nclx,cnstd(3):cnstd(4),nclz))

	! Pack stresses into a dummy 3D cube with 6 surfaces and 3 stresses per surface
	allocate(stress(3,6,nclx,ncly,nclz)); stress = 0.d0
	stress(1,1,:,:,1) = tau_xx(:,:,1)
	stress(1,3,:,:,1) = tau_xx(:,:,3)
	stress(2,1,:,:,1) = tau_yx(:,:,1)
	stress(2,3,:,:,1) = tau_yx(:,:,3)
	stress(1,2,:,:,1) = tau_xy(:,:,2)
	stress(1,4,:,:,1) = tau_xy(:,:,4)
	stress(2,2,:,:,1) = tau_yy(:,:,2)
	stress(2,4,:,:,1) = tau_yy(:,:,4)

	sendbuf = 0.d0
	do i=1,nclx
	do j=cnstd(3),cnstd(4)
	do k=1,nclz
		jj = j - cnstd(3)+1
		sendbuf(:,i,j,k) = reshape(stress(:,:,i,jj,k),(/ npercell /))
	enddo
	enddo
	enddo

	!Send stress tensor to MD code
	call CPL_send( sendbuf,                                 &
	               icmin_send=cnstd(1),icmax_send=cnstd(2), &
	               jcmin_send=cnstd(3)-cnstd(3)+1,jcmax_send=cnstd(4)-cnstd(3)+1, &
	               kcmin_send=cnstd(5),kcmax_send=cnstd(6), &
	               send_flag=send_flag                        )

end subroutine socket_coupler_send_stress


!==================================================================================
! Get Boundary condition for continuum from average of MD 
! receives average velocities between jmin0:jmin and copies the in boundary 
! condition arrays  
!----------------------------------------------------------------------------------

subroutine socket_coupler_get_md_BC(u,v)
    use CPL, only : CPL_get,CPL_recv,CPL_realm,error_abort,cpl_proc_extents,VOID, & 
					printf,CPL_OLAP_COMM,xg,xL_cfd,yg,yL_cfd,zg,zL_cfd !DEBUG DEBUG
    implicit none

    real(kind(0.d0)),dimension(:),intent(out)  	      :: u,v

	logical		  								      :: recv_flag
	integer											  :: i,j,k,ii,ib,jb
	integer											  :: nclx,ncly,nclz,pcoords(3),extents(6)
	integer											  :: i1,i2,j1,j2,k1,k2
	integer											  :: jcmin_recv,jcmax_recv
    real(kind(0.d0))								  :: dy, dvdy
    real(kind(0.d0)), allocatable, dimension(:,:,:,:) :: uvw_md, ucvcwc_md
    real(kind(0.d0)), allocatable, dimension(:,:,:,:) :: buf1, buf2
	real											  :: uvw_BC(4)

	integer		:: bufsize, jcmin_olap
	character	:: str_bufsize

	!Setup extents -- halo is one below minimum cell in domain
	call CPL_get(jcmin_olap=jcmin_olap)
	jcmin_recv = jcmin_olap; jcmax_recv = jcmin_olap

	!Get number of CFD number of cells ready to receive data
	pcoords = (/ 1,1,1 /)
	call CPL_proc_extents(pcoords,CPL_realm(),extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1

	!Allocate array to CFD number of cells ready to receive data including halos
	allocate(uvw_md(4,nclx,ncly,nclz)); uvw_md = VOID

	!Receive data from MD
	call CPL_recv(uvw_md,jcmax_recv=jcmax_recv,jcmin_recv=jcmin_recv,recv_flag=recv_flag)
	if (any(uvw_md(:,2:nclx+1,jcmin_recv:jcmax_recv,2:nclz+1) .eq. VOID)) & 
		call error_abort("socket_coupler_get_md_BC error - VOID value copied to uc,vc or wc")

	u(:) = uvw_md(1,:,1,1)/uvw_md(4,:,1,1)
	v(:) = uvw_md(2,:,1,1)/uvw_md(4,:,1,1)

end subroutine socket_coupler_get_md_BC
#endif
end module continuum_coupler_socket

