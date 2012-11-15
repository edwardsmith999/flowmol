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
!  Read coupler input files
!-----------------------------------------------------------------------------
subroutine socket_read_coupler_input
	use messenger
	use CPL, only : read_coupler_input
	implicit none

	call read_coupler_input		! Read COUPLER.in input file

end subroutine socket_read_coupler_input

!=============================================================================
! Call coupler initialise to swap all setup data
!-----------------------------------------------------------------------------
subroutine socket_coupler_init
	use CPL, only : coupler_cfd_init,CPL_create_map,error_abort, density_cfd
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
   ! write(0,*) 'CFD socket nsteps, dt ', nsteps, dt
	kTmin_1(1) = kmin; kTmax_1(1) = kmax-1

	!Define compound arrays to make passing more concise
	ijkmin = (/ imin, jmin, kmin /)	
	ijkmax = (/ imax, jmax, kmax /)-1 !Minus 1 as coupler requires max cells NOT max surfaces
	npxyz  = (/ npx , npy , npz  /)
	ngxyz  = (/ ngx , ngy , ngz  /)-1 !Minus 1 as coupler requires no. of cells NOT no. surfaces
	xyzL   = (/  xL ,  yL ,  zL  /)

    call coupler_cfd_init(nsteps,dt,icomm_grid,icoord,npxyz,xyzL,ngxyz, & 
	                      density_cfd,ijkmax,ijkmin,iTmin_1,iTmax_1,jTmin_1,&
	                      jTmax_1,kTmin_1,kTmax_1,xpg,ypg,zpg)

	! Establish mapping between MD an CFD
	call CPL_create_map

end subroutine socket_coupler_init


!=============================================================================
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!
!							SIMULATION
!
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!=============================================================================

!---------------------------------------------------------------------
! Get Boundary condition for continuum from average of MD 

subroutine  socket_coupler_get_md_BC(uc,vc,wc)
    use CPL, only : CPL_recv,CPL_realm,error_abort,cpl_proc_extents,VOID
	use data_export, only : nixb, niyb, nizb,iblock,jblock,kblock, & 
							i1_u,i2_u,j1_u,j2_u, & 
							i1_v,i2_v,j1_v,j2_v, & 
							i1_w,i2_w,j1_w,j2_w, & 
							i1_T,i2_T,j1_T,j2_T,k1_T,k2_T
    implicit none

    real(kind(0.d0)),dimension(0:,0:,0:),intent(out)  :: uc,vc,wc 

	logical		  								      :: recv_flag
    logical, save 								      :: firsttime = .true.
	integer											  :: i,j,k,nclx,ncly,nclz,pcoords(3),extents(6)
	integer											  :: jcmin_recv,jcmax_recv
    real(kind(0.d0)), allocatable, dimension(:,:,:,:) :: uvw_md
	real											  :: uvw_BC(4)

	integer		:: bufsize
	character	:: str_bufsize

	jcmin_recv = 1; jcmax_recv = 1

	!Allocate array to CFD number of cells ready to receive data
	pcoords = (/ iblock,jblock,kblock /)
	call CPL_proc_extents(pcoords,CPL_realm(),extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1
	allocate(uvw_md(4,nclx,ncly,nclz)); uvw_md = VOID

	call CPL_recv(uvw_md,jcmax_recv=jcmax_recv,jcmin_recv=jcmin_recv,recv_flag=recv_flag)
	!call CPL_gather(uvw_md,3)

	!Set full extent of halos to zero and set domain portion to MD values
	uc(:,:,0) = 0.d0; vc(:,:,1) = 0.d0; wc(:,:,0) = 0.d0

	!Average all cells on a CFD processor to give a single BC
	uvw_BC(1) = sum(uvw_md(1,:,1,:)) 
	uvw_BC(2) = sum(uvw_md(2,:,1,:))
	uvw_BC(3) = sum(uvw_md(3,:,1,:))
	uvw_BC(4) = sum(uvw_md(4,:,1,:)) 

	!Average in x so all processor have same global average BC
	call globalDirSum(uvw_BC,4,1)

	uc(:,:,0) = uvw_BC(1)/uvw_BC(4)
	vc(:,:,1) = uvw_BC(2)/uvw_BC(4)
	wc(:,:,0) = uvw_BC(3)/uvw_BC(4)

	!if (rank_realm .eq. 1) then
	!	print'(i4,a,7f10.3)', rank_world,'global average BC',uc(5,10,0),vc(5,10,1),wc(5,10,0),uvw_BC
	!endif

	if (any(uc .eq. VOID) .or. &
		any(uc .eq. VOID) .or. &
		any(uc .eq. VOID)) call error_abort("socket_coupler_get_md_BC error - VOID value copied to uc,vc or wc")

	!print'(a,26i5)', 'array extents',rank_world,shape(uvw_md),nixb, niyb, nizb, & 
	!														   i1_u,i2_u,j1_u,j2_u, & 
	!														   i1_v,i2_v,j1_v,j2_v, & 
	!														   i1_w,i2_w,j1_w,j2_w, & 
	!														   i1_T,i2_T,j1_T,j2_T,k1_T,k2_T

	! u interval [i1_u, i2_u], or [2,  ngx ] ??? 4 Procs = [2 32][3 34][3 34][3 35] ???
	! v interval [i1_v, i2_v], or [1, ngx-1] ??? 4 Procs = [1 32][3 34][3 34][3 34] ???
	! w interval [i1_w, i2_w], or [1, ngx-1] ??? 4 Procs = [1 32][3 34][3 34][3 34] ???
	!uc(ngz  ,nlx+1,nly  )  
	!vc(ngz  ,nlx  ,nly+1)
	!wc(ngz+1,nlx  ,nly  )

	!call printf(uc(:,jcmax_recv,4))
	!call printf(vc(:,jcmax_recv,4))
	!call printf(wc(:,jcmax_recv,4))

	!Transposed indices
	!ix = 2; iy = 3; iz = 1

	!Coupler passes cell centered values so set CFD halo directly
	!do i=0,size(uc,2)-1
	!do j=jcmin_recv-1,jcmax_recv-1
	!do k=0,size(uc,1)-1
	!	uc(k,i,j) = uvw_md(1,i+1,j+1,k+1)
	!	print'(3i8,2f20.8)', i,j,k,uc(k,i,j),uvw_md(1,i+1,j+1,k+1)
	!enddo
	!enddo
	!enddo


end subroutine socket_coupler_get_md_BC


!---------------------------------------------------------------------
! Send continuum velocity in x direction to be used by MD

subroutine socket_coupler_send_velocity
    use CPL, only : CPL_send,CPL_olap_extents,CPL_overlap,CPL_get,CPL_realm
 	use data_export, only : uc,i1_u,i2_u,ngz,nlx,nlxb,ibmin_1,ibmax_1,iblock,jblock,kblock
    implicit none

	logical	:: send_flag
    integer	:: i,j,n,ixyz,icell,jcell,kcell,npercell,nclx,ncly,nclz
    integer	:: coord(3),extents(6),cnstnd_cells,jcmin_send,jcmax_send,jcmax_olap
    real(kind(0.d0)),dimension(:,:,:,:), allocatable 	:: sendbuf
	real(kind(0.d0)) :: gradient

	! Check processor is inside MD/CFD overlap zone 
	if (.not.(CPL_overlap())) return
	call CPL_get(jcmax_olap=jcmax_olap)

	!Three velocity components
	npercell = 3

	!Allocate array for size of data on local processor
	coord = (/iblock,jblock,kblock /)
	call CPL_olap_extents(coord,CPL_realm(),extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1
	allocate(sendbuf(npercell,nclx,ncly,nclz))
 
	!Number of cells to package and send
	cnstnd_cells = 1
	jcmin_send = jcmax_olap-1-cnstnd_cells
	jcmax_send = jcmax_olap-1

	!Interpolate cell centres using all surface
	sendbuf = 0.d0
	do j=jcmin_send,jcmax_send
	do i=1,nclx
		n = i + i1_u - 1
		sendbuf(1,i,j,:) = 0.5d0*(uc(:,n,j) + uc(:,n+1,j))
	enddo
	!call printf(sendbuf(1,40,j,:))
	enddo

	!Send data to MD
	call CPL_send(sendbuf,jcmax_send=jcmax_send, & 
					      jcmin_send=jcmin_send,send_flag=send_flag)

end subroutine socket_coupler_send_velocity


!---------------------------------------------------------------------
! Send continuum Stress in all directions to be used by MD

subroutine socket_coupler_send_stress
    use CPL, only : CPL_send,CPL_olap_extents,CPL_overlap,CPL_get,CPL_realm
 	use data_export, only : uc,vc,wc,visc,iblock,jblock,kblock
    implicit none

	logical			 	  :: send_flag
    integer			 	  :: i,n,ixyz,icell,jcell,kcell,npercell,nclx,ncly,nclz
    integer			 	  :: coord(3),extents(6),jcmax_olap
	real(kind(0.d0))	  :: density_cfd
    real(kind(0.d0)),dimension(:,:,:,:), allocatable 	:: sendbuf
	real(kind(0.d0)),dimension(:,:,:,:,:), allocatable 	:: dUidxj

	! Check processor is inside MD/CFD overlap zone 
	if (.not.(CPL_overlap())) return
	call CPL_get(jcmax_olap=jcmax_olap,density_cfd=density_cfd)

	!Three by three stress tensor components
	npercell = 9

	!Allocate array for size of data on local processor
	coord = (/iblock,jblock,kblock /)
	call CPL_olap_extents(coord,CPL_realm(),extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1
	allocate(sendbuf(npercell,nclx,ncly,nclz))

	!Get strain tensor at cell centers
	call Evaluate_strain(uc,vc,wc,dUidxj)

	!Get stress and store in buffer
	! QUESTION -- if this visc is kinematic (rho/meu) then 
	! is this the stress? Does DNS assume unit density?
	sendbuf = visc*reshape(dUidxj,(/ 9,nclx,ncly,nclz /))

	!Send stress tensor to MD code
	call CPL_send(sendbuf,jcmax_send=jcmax_olap-1, & 
					      jcmin_send=jcmax_olap-1,send_flag=send_flag)

end subroutine socket_coupler_send_stress

!---------------------------------------------------------------------
! 			Gets cell centered strain from velocity field

subroutine Evaluate_strain(uc,vc,wc,dUidxj)
	use data_export, only : ibmin,ibmax,jbmin,jbmax,kbmin,kbmax, &
						   	imap_1,jmap_1,ibmap_1,jbmap_1,npx,npy, & 
							ngx,ngy,ngzm,iblock,jblock,kblock,vp, &
							suxix,suxiy,svetax,svetay,spz
	implicit none

	real(kind(0.d0)),intent(in),dimension(:,:,:)					:: uc,vc,wc
	real(kind(0.d0)),intent(out),dimension(:,:,:,:,:),allocatable	:: dUidxj

	integer			:: i,j,k,ib,jb,i1,i2,j1,j2
	real(kind(0.d0)):: voli, ds, Eval_eps
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
	allocate(dUidxj(3,3,ibmin:ibmax,jbmin:jbmax,ngzm)); dUidxj = 0.d0

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
		dUidxj(1,:,ib,jb,k) = (/dudx,dudy,dudz/)
		dUidxj(2,:,ib,jb,k) = (/dvdx,dvdy,dvdz/)
		dUidxj(3,:,ib,jb,k) = (/dwdx,dwdy,dwdz/)

	end do
	end do
	end do

end subroutine Evaluate_strain

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
