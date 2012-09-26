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

	call coupler_create_comm(COUPLER_CFD,CFD_COMM,ierr)
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
    write(0,*) 'CFD socket nsteps, dt ', nsteps, dt
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
	call coupler_create_map

end subroutine socket_create_map


!=============================================================================
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!
!							SIMULATION
!
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!=============================================================================

subroutine socket_coupler_send_velocity
    use coupler
    use data, only : uc, i1_u, i2_u, j1_u, j2_u, ibmap_1, jbmap_1, ngz
    use messenger, only : icomm_grid, icoord
	use coupler_module, only : myid_grid
    implicit none

    integer i1b_u, j1b_u, i, i1, jmax_ovr,jo,js,je, ierr
    real(kind(0.d0)),allocatable :: buff(:,:,:,:)

    i1b_u=ibmap_1(i1_u)
    j1b_u=jbmap_1(j1_u)
    i1   = i1_u
    
    call coupler_cfd_get(jmax_overlap=jmax_ovr)

    !uc(1:ngz-1,i1_u:i2_u,j1_u+2) = myid+13.d0

    ! this is to catch the boundary condtion at x=0 (uc start from i=2 in global grid)
    ! needs some further discusion
    if (icoord(1,myid_grid)==1) then
        i1 = i1_u-1
        !uc(1:ngz-1,i1:i1,j1_u+2) = myid_grid+66.d0
    endif
    
    je = j1_u + jmax_ovr-jbmap_1(j1_u)-2
    js = je

	!print*, 'CFD constraint', js,j1_u,jmax_ovr,jbmap_1(j1_u)
    !write(0,*)'cfd socket, jo:',jo
    call coupler_send(uc(1:ngz-1,i1:i2_u,js:je),index_transpose=(/2,3,1/))

    !do i=i1, i2_u
    !    write(600+myid_grid,'(1000(E11.3,1x))') uc(1:ngz-1,i,j1_u+2)
    !enddo
    ! write(600+myid_grid,'(1000(E11.3,1x))')
    ! call flush(600+myid_grid)

end subroutine socket_coupler_send_velocity


subroutine  socket_coupler_get_md_BC(uc,vc,wc)
    use data, only : ngz,i1_u,i2_u,i1_v,i2_v,i1_w,i2_w, ibmap_1
    use messenger, only : icomm_grid, icoord
    use coupler
	use coupler_module, only : myid_grid
    implicit none

    real(kind(0.d0)), intent(inout) :: uc(0:,0:,0:),vc(0:,0:,0:),wc(0:,0:,0:) 
    integer j, icell
    integer, save :: i1g_u,i1g_v,i1g_w,i1_ul,ierr
    logical, save :: firsttime = .true.

	integer, dimension(3)	:: indices
    real(kind(0.d0)), allocatable, dimension(:,:,:,:) :: uvwbuff

	allocate(uvwbuff(4,size(wc,1),size(uc,2),size(vc,3)))

    if ( firsttime) then
        firsttime = .false.
        i1g_u = ibmap_1(i1_u)
        i1g_v = ibmap_1(i1_v)
        i1g_w = ibmap_1(i1_w)
        !write(0,*) 'i indices ngz, i1_u,i2_u,i1_v,i2_v,i1_w,i2_w ', ngz, i1_u,i2_u,i1_v,i2_v,i1_w,i2_w 
        !write(0,*) 'size uc ...', size(uc,2), size(vc,2), size(wc,2)
        !write(0,*) 'global indices i1g_u, i1g_v, i1g_w ',  i1g_u, i1g_v, i1g_w
        
        ! this is to catch the boundary condtion at x=0 (uc start from i=2 in global grid)
        ! needs some further discusion !!!!!!!!!! WHEN?
        i1_ul = i1_u
        if (icoord(1,myid_grid)==1) then
            i1_ul = i1_u-1
        endif

        endif
	!print'(2a,2i8,4f25.16)', 'CFD befr data',code_name(COUPLER_REALM), myid_grid, & 
	!							size(uc(:,:,0)), maxval(uc(:,:,0)),minval(uc(:,:,0)),sum(uc(:,:,0)),uc(3,3,0)
    call coupler_recv(uc,index_transpose=(/2,3,1/))
 	!print*, 'Extents of array', ngz-1, i1_ul,i2_u
	!uc(:,:,0:0) = 0.d0
	indices = (/2,3,1/)

	!========UNDER DEVELOPMENT=======
	!call coupler_recv_data(uvwbuff,index_transpose=indices)

	!print*, 'extents of CFD BC loop', i1_ul,i2_u
	!do icell=i1_ul,i2_u
	!	write(110+myid_grid,'(a,2i5,11f10.4)') '1 CFD',myid_grid,icell,uvwbuff(1,:,icell,1)
	!	write(110+myid_grid,'(a,2i5,11f10.4)') '2 CFD',myid_grid,icell,uvwbuff(2,:,icell,1)
	!	write(110+myid_grid,'(a,2i5,11f10.4)') '3 CFD',myid_grid,icell,uvwbuff(3,:,icell,1)
	!	write(110+myid_grid,'(a,2i5,11f10.4)') '4 CFD',myid_grid,icell,uvwbuff(4,:,icell,1)
		!print'(a,2i8,6f10.5)', 'CFD grid location',irank,icell,xpg(icell,5),xpu(icell,5),x(icell),ypg(icell,5),ypu(icell,5),y(icell)
	!enddo

	!uc(1:ngz-1,i1_ul:i2_u,0:0) 	= uvwbuff(1,1:ngz-1,i1_ul:i2_u,:)/uvwbuff(4,1:ngz-1,i1_ul:i2_u,:)
	!vc(1:ngz-1,i1_v :i2_v,0:1) 	= uvwbuff(2,1:ngz-1,i1_v :i2_v,:)
	!wc(1:ngz  ,i1_w :i2_w,0:0) 	= uvwbuff(3,1:ngz-1,i1_w :i2_w,:)
	!========UNDER DEVELOPMENT=======

    !call coupler_recv_data(uc(1:ngz-1,i1_ul:i2_u,0:0),index_transpose=indices,accumulate=.true.,pbc=1)

	!print'(a,2i8,4f25.16)', 'CFD recv MD     ', myid_grid,size(uc(:,:,0)), & 
	!					maxval(uc(1:ngz-1,i1_ul:i2_u,0:0)),minval(uc(1:ngz-1,i1_ul:i2_u,0:0)),sum(uc(1:ngz-1,i1_ul:i2_u,0:0)),uc(3,3,0)
    !call coupler_recv_data(vc(1:ngz-1,i1_v:i2_v,0:1),index_transpose=indices,accumulate=.true.)
    !call coupler_recv_data(wc(1:ngz,i1_w:i2_w,0:0),index_transpose=indices,accumulate=.true.,pbc=3)

	do icell=i1_ul,i2_u
		write(110+myid_grid,'(a,2i5,11f10.4)') '1 CFDuc',myid_grid,icell,uc(1:ngz-1,icell,0)
	enddo

	do icell=i1_v,i2_v
		write(110+myid_grid,'(a,2i5,11f10.4)') '1 CFDvc',myid_grid,icell,vc(1:ngz-1,icell,0)
		write(110+myid_grid,'(a,2i5,11f10.4)') '2 CFDvc',myid_grid,icell,vc(1:ngz-1,icell,1)
	enddo

	do icell=i1_w,i2_w
		write(110+myid_grid,'(a,2i5,12f10.4)') '1 CFDwc',myid_grid,icell,wc(1:ngz,icell,0)
	enddo

    !debug writes
    !do j=i1_ul, i2_u
    !    write(200+myid_grid,'(1000(E11.3,1x))') uc(1:ngz-1,j,0)
    !enddo
    !write(200+myid_grid,*) 
    !call flush(200+myid_grid)
    !do j=i1_v,i2_v
    !    write(300+myid_grid,'(1000(E11.3,1x))') vc(1:ngz-1,j,0)
    !enddo
    ! write(300+myid_grid,*)
    !do j=i1_v,i2_v
    !    write(400+myid_grid,'((1000(E11.3,1x)))') vc(1:ngz-1,j,1)
    !enddo
    ! write(400+myid_grid,*)
    !do j=i1_w,i2_w
    !    write(500+myid_grid,'((1000(E11.3,1x)))') wc(1:ngz,j,0)
    !enddo
    !write(500+myid_grid,*)

end subroutine socket_coupler_get_md_BC

#if COUPLER_DEBUG_LA
!---------------------------------------------------------------------
! debug subroutine that writes uc for gnuplot
    subroutine socket_coupler_write_uc(ntime)
        use data, only : uc,ngz,i1_u,i2_u, j1_u,j2_u
        use data_export, only : kmin, kmax
        use messenger, only : icomm_grid, icoord
        use computation_parameters, only : prefix_dir
		use coupler_module, only : myid_grid
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


        if (icoord(1,myid_grid)==1) then
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
