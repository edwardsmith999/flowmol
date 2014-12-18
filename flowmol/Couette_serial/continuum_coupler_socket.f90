module continuum_coupler_socket
        implicit none

contains
    !LINSPACE Linearly spaced vector.
    !
    !   LINSPACE(X1, X2, N) generates N points between X1 and X2.
    !   For N = 1, LINSPACE returns X2.
    !   Based on the MATLAB function

    function linspace(d1, d2, n)
        implicit none

        integer,intent(in)              :: n
        double precision,intent(in)     :: d1,d2
        double precision,dimension(:),allocatable    :: linspace

        integer                                     :: i, n1

        n1 = n-1
        if (n1 .eq. 0) stop "Error in linspace -- Ensure n > 1"
        
        allocate(linspace(n))
        do i = 0,n1
            linspace(i+1) = d1 + i * (d2-d1)/(n1)
        enddo

       ! print'(20f7.2)',linspace

    end function linspace


#if USE_COUPLER

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
    use computational_constants, only : nx, ny, nz, & 
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
    real(kind(0.d0)) :: z0(nz+1)
          
	!Define compound arrays to make passing more concise
	ijkmin = (/ 1, 1, 1 /)
	!Minus 1 as coupler requires max cells NOT max cell vertices
	ijkmax = (/ nx, ny, nz /)
	iTmin = ijkmin(1); iTmax = ijkmax(1)
	jTmin = ijkmin(2); jTmax = ijkmax(2)
	kTmin = ijkmin(3); kTmax = ijkmax(3)
	npxyz  = (/ npx , npy , npz  /)
	!Minus 1 as coupler requires no. of cells NOT no. cell vertices
	ngxyz  = (/ nx , ny , nz   /) 
	xyzL   = (/ lx , ly , lz  /)
    z0 = linspace(0.d0,lz,nz+1)

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
	use computational_constants, only : continuum_iter
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
	allocate(sendbuf(npercell,extents(1):extents(2), &
		                      extents(3):extents(4), &
		                      extents(5):extents(6)))

	!do j = 1,10
    !    print'(a,4i6,3e27.10)','all  vel  ',continuum_iter, 0,j,0, u(2,j),v(2,j),0.d0
	!enddo

	!Copy cell centered velocity to buffer
	sendbuf(:,:,:,:) = 0.d0
	do i=cnstd(1),cnstd(2)
	do j=cnstd(3),cnstd(4)
	do k=cnstd(5),cnstd(6)
		sendbuf(1,i,j,k) = u(i+1,j+1)
		sendbuf(2,i,j,k) = v(i+1,j+1)
		sendbuf(3,i,j,k) = 0.d0
		!print'(a,3i8,3f10.5)', 'Pre send vel  ', i+1,j+1,k+1,sendbuf(:,i,j,k)
	enddo
	enddo
	enddo

!	call CPL_send( sendbuf,jcmin_send=cnstd(3), & 
!                           jcmax_send=cnstd(4),send_flag=send_flag)
	call CPL_send( sendbuf,                                 &
	               icmin_send=cnstd(1),icmax_send=cnstd(2), &
	               jcmin_send=cnstd(3),jcmax_send=cnstd(4), &
	               kcmin_send=cnstd(5),kcmax_send=cnstd(6), &
	               send_flag=send_flag                        )

end subroutine socket_coupler_send_velocity


!---------------------------------------------------------------------
! Send continuum Stress in all directions to be used by MD

subroutine socket_coupler_send_stress(tau_xx,tau_xy,tau_yx,tau_yy)
    use CPL, only : CPL_send,CPL_olap_extents,CPL_overlap,CPL_get,CPL_realm
    use physical_constants, only : Re
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
	allocate(sendbuf(npercell,extents(1):extents(2), &
		                      extents(3):extents(4), &
		                      extents(5):extents(6)))

	! Pack stresses into a dummy 3D cube with 6 surfaces and 3 stresses per surface
	allocate(stress(extents(1):extents(2), &
		            extents(3):extents(4), &
		            extents(5):extents(6),3,6))
    stress = 0.d0
	sendbuf = 0.d0
	do i=cnstd(1),cnstd(2)
	do j=cnstd(3),cnstd(4)
	do k=cnstd(5),cnstd(6)
	    stress(i,j,k,1,1) = tau_xx(i+1,j+1,1)
	    stress(i,j,k,1,4) = tau_xx(i+1,j+1,3)
	    stress(i,j,k,2,1) = tau_yx(i+1,j+1,1)
	    stress(i,j,k,2,4) = tau_yx(i+1,j+1,3)
	    stress(i,j,k,1,2) = tau_xy(i+1,j+1,2)
	    stress(i,j,k,1,5) = tau_xy(i+1,j+1,4)
	    stress(i,j,k,2,2) = tau_yy(i+1,j+1,2)
	    stress(i,j,k,2,5) = tau_yy(i+1,j+1,4)
		sendbuf(:,i,j,k) = (1.d0/Re)*reshape(stress(i,j,k,:,:),(/ npercell /))
	enddo
	enddo
	enddo

	!Send stress tensor to MD code
	call CPL_send( sendbuf,                                 &
	               icmin_send=cnstd(1),icmax_send=cnstd(2), &
	               jcmin_send=cnstd(3),jcmax_send=cnstd(4), &
	               kcmin_send=cnstd(5),kcmax_send=cnstd(6), &
	               send_flag=send_flag                        )

end subroutine socket_coupler_send_stress


!==================================================================================
! Get Boundary condition for continuum from average of MD 
! receives average velocities between jmin0:jmin and copies the in boundary 
! condition arrays  
!----------------------------------------------------------------------------------

subroutine socket_coupler_get_md_BC(u,v)
    use CPL, only : CPL_get,CPL_recv,CPL_realm,error_abort,cpl_proc_extents,VOID
    implicit none

    real(kind(0.d0)),dimension(:),intent(out)  	      :: u,v

	logical		  								      :: recv_flag
	integer											  :: i
	integer											  :: nclx,ncly,nclz,pcoords(3),extents(6)
	integer											  :: jcmin_recv,jcmax_recv
    real(kind(0.d0)), allocatable, dimension(:,:,:,:) :: uvw_md
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
	if (any(uvw_md(:,:,jcmin_recv:jcmax_recv,:) .eq. VOID)) & 
		call error_abort("socket_coupler_get_md_BC error - VOID value copied to uc,vc or wc")

	do i = 1,nclx
		if (uvw_md(4,i,1,1) .ne. 0) then
			u(i) = uvw_md(1,i,1,1)/uvw_md(4,i,1,1)
			v(i) = uvw_md(2,i,1,1)/uvw_md(4,i,1,1)
			!w(i) = uvw_md(3,i,1,1)/uvw_md(4,i,1,1)
		else
			u(i) = 0.d0
			v(i) = 0.d0
			!w(i) = 0.d0
		endif
	enddo


end subroutine socket_coupler_get_md_BC
#endif

end module continuum_coupler_socket

