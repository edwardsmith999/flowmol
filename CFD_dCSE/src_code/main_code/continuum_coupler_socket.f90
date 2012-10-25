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
	use coupler
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
	use coupler_input_data, only : read_coupler_input
    use coupler_module, only : request_stop
	implicit none

	call read_coupler_input		! Read COUPLER.in input file

    ! stop if requested ( useful for development )
	call request_stop("create_comm") ! stops here if in COUPLER.in stop requestis set to "create_comm"

end subroutine socket_read_coupler_input

!=============================================================================
! Call coupler initialise to swap all setup data
!-----------------------------------------------------------------------------
subroutine socket_coupler_init
	use coupler_internal_cfd, only : coupler_cfd_init
    use messenger, only : icomm_grid, icoord
    use data_export, only : imin,imax,jmin,jmax,kmin,kmax, &
							iTmin_1,iTmax_1,jTmin_1,jTmax_1, &
                            ngx,ngy,ngz,xL,yL,zL, &
                            npx,npy,npz,dt,xpg,ypg,zpg
	use mesh_export, only : xL, yL, zL
	use coupler_input_data, only : density_tag, density
	use coupler_module, only : error_abort
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

	!Density is NOT defined for the DNS code
	if (density_tag == CPL) then 
		!Do nothing
		density = density
	else
		call error_abort("Density not specified in coupler")
	endif

    call coupler_cfd_init(nsteps,dt,icomm_grid,icoord,npxyz,xyzL,ngxyz,density, & 
							   ijkmax,ijkmin,iTmin_1,iTmax_1,jTmin_1,jTmax_1,kTmin_1,kTmax_1,xpg,ypg,zpg)

end subroutine socket_coupler_init


!=============================================================================
! Establish mapping between CFD and MD
!-----------------------------------------------------------------------------
subroutine socket_create_map
	implicit none

	!Note coupler cannot be called directly so this socket is needed
	call CPL_create_map

end subroutine socket_create_map


!=============================================================================
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!
!							SIMULATION
!
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!=============================================================================


! ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬(-_-)▬▬ஜ۩۞۩ஜ▬▬(●̪•)▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬ 
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
	!print'(2a,5i8)', 'CFD SIDE',realm_name(realm), rank_world, olap_mask(rank_world),coord

	if (olap_mask(rank_world) .eq. 0) return

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
	if (olap_mask(rank_world) .eq. 0) return

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

 	if (olap_mask(rank_world).ne.1) return

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

subroutine socket_coupler_send_velocity
    use coupler
    use data, only : uc, i1_u, i2_u, j1_u, j2_u, ibmap_1, jbmap_1, ngz
    use messenger, only : icomm_grid, icoord
	use coupler_module, only : rank_realm
    implicit none

    integer i1b_u, j1b_u, i, i1, jmax_ovr,jo,js,je, ierr
    real(kind(0.d0)),allocatable :: buff(:,:,:,:)

    i1b_u=ibmap_1(i1_u)
    j1b_u=jbmap_1(j1_u)
    i1   = i1_u
    
    call coupler_cfd_get(jcmax_overlap=jmax_ovr)

    !uc(1:ngz-1,i1_u:i2_u,j1_u+2) = rank_world+13.d0

    ! this is to catch the boundary condtion at x=0 (uc start from i=2 in global grid)
    ! needs some further discusion
    if (icoord(1,rank_realm)==1) then
        i1 = i1_u-1
        !uc(1:ngz-1,i1:i1,j1_u+2) = rank_realm+66.d0
    endif
    
    je = j1_u + jmax_ovr-jbmap_1(j1_u)-2
    js = je

	!print*, 'CFD constraint', js,j1_u,jmax_ovr,jbmap_1(j1_u)
    !write(0,*)'cfd socket, jo:',jo
    !call CPL_send(uc(1:ngz-1,i1:i2_u,js:je),index_transpose=(/2,3,1/))
    !call CPL_send(uc,index_transpose=(/2,3,1/))
    !do i=i1, i2_u
    !    write(600+rank_realm,'(1000(E11.3,1x))') uc(1:ngz-1,i,j1_u+2)
    !enddo
    ! write(600+rank_realm,'(1000(E11.3,1x))')
    ! call flush(600+rank_realm)

end subroutine socket_coupler_send_velocity



!---------------------------------------------------------------------
! Get Boundary condition for continuum from average of MD 

subroutine  socket_coupler_get_md_BC(uc,vc,wc)
    use coupler, only : CPL_recv
	use coupler_module, only : rank_realm,olap_mask,rank_world,printf, & 
							   iblock_realm,jblock_realm,kblock_realm,error_abort
	use data_export, only : nixb, niyb, nizb, & 
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
    real(kind(0.d0)), allocatable, dimension(:,:,:,:) :: uvwbuff
	real											  :: uvw_BC(4)

	integer		:: bufsize
	character	:: str_bufsize

	jcmin_recv = 1; jcmax_recv = 1

	!Allocate array to CFD number of cells ready to receive data
	pcoords = (/ iblock_realm,jblock_realm,kblock_realm /)
	call CPL_proc_extents(pcoords,cfd_realm,extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1
	allocate(uvwbuff(4,nclx,ncly,nclz)); uvwbuff = VOID

	call CPL_recv(uvwbuff,jcmax_recv=jcmax_recv,jcmin_recv=jcmin_recv,recv_flag=recv_flag)
	!call CPL_gather(uvwbuff,3)
	!call printf(uvwbuff(1,:,jcmax_recv,4))
	!call printf(uvwbuff(4,:,jcmax_recv,4))

	!Set full extent of halos to zero and set domain portion to MD values
	uc(:,:,0) = 0.d0; vc(:,:,1) = 0.d0; wc(:,:,0) = 0.d0

	!Average all cells on a CFD processor to give a single BC
	uvw_BC(1) = sum(uvwbuff(1,:,1,:)) 
	uvw_BC(2) = sum(uvwbuff(2,:,1,:))
	uvw_BC(3) = sum(uvwbuff(3,:,1,:))
	uvw_BC(4) = sum(uvwbuff(4,:,1,:)) 
	print'(i4,a,4f10.3)', rank_world,'per proc BC',uvw_BC

	!Average in x so all processor have same global average BC
	call globalDirSum(uvw_BC,4,1)

	uc(:,:,0) = uvw_BC(1)/uvw_BC(4)
	vc(:,:,1) = uvw_BC(2)/uvw_BC(4)
	wc(:,:,0) = uvw_BC(3)/uvw_BC(4)

	if (rank_realm .eq. 1) then
		print'(i4,a,7f10.3)', rank_world,'global average BC',uc(5,10,0),vc(5,10,1),wc(5,10,0),uvw_BC
	endif

	if (any(uc .eq. VOID) .or. &
		any(uc .eq. VOID) .or. &
		any(uc .eq. VOID)) call error_abort("socket_coupler_get_md_BC error - VOID value copied to uc,vc or wc")

	!print'(a,26i5)', 'array extents',rank_world,shape(uvwbuff),nixb, niyb, nizb, & 
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
	!	uc(k,i,j) = uvwbuff(1,i+1,j+1,k+1)
	!	print'(3i8,2f20.8)', i,j,k,uc(k,i,j),uvwbuff(1,i+1,j+1,k+1)
	!enddo
	!enddo
	!enddo


end subroutine socket_coupler_get_md_BC

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
