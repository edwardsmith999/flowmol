module continuum_coupler_socket
    implicit none
#if USE_COUPLER

contains


    subroutine socket_coupler_init 
        use messenger, only : icomm_grid, icoord
        use data_export, only : imino, imin, imax, imaxo,&
                                jmino, jmin, jmax, jmaxo,& 
                                kmino, kmin, kmax, kmaxo,&
                                npx, npy, npz, dt
        use mesh_export, only : x, y, z, dx, dz
!        use simulation, only  : nsteps
        use coupler, only : coupler_cfd_init, coupler_create_map 
        implicit none
 
        integer nsteps

        call readInt("nsteps", nsteps)
       
        write(0,*) 'CFD socket nsteps, dt ', nsteps, dt
        
        call coupler_cfd_init(icomm_grid,imino=imino,imin=imin,imax=imax,&
            imaxo=imaxo,jmino=jmino,jmin=jmin,jmax=jmax,jmaxo=jmaxo,&
            kmino=kmino,kmin=kmin,kmax=kmax,kmaxo=kmaxo,nsteps=nsteps,&
            x=x,y=y,z=z,dx=dx,dz=dz,npx=npx,npy=npy,npz=npz,&
            icoord=icoord,dt=dt)
        
        write(0,*) 'CFD socket, after cfd init'

        call coupler_create_map

        write(0,*) 'CFD socket, after create map'
                
    end subroutine socket_coupler_init


    subroutine socket_coupler_send_velocity
        use coupler
        use data, only : uc, i1_u, i2_u, j1_u, j2_u, ibmap_1, jbmap_1, ngz
        use messenger, only : icomm_grid, icoord
        implicit none

        integer i1b_u, j1b_u, i, i1, jmax_ovr,jo, myid,ierr
        real(kind(0.d0)),allocatable :: buff(:,:,:,:)

        i1b_u=ibmap_1(i1_u)
        j1b_u=jbmap_1(j1_u)
        i1   = i1_u
        
        call coupler_cfd_get(jmax_overlap=jmax_ovr)
        call mpi_comm_rank(icomm_grid,myid,ierr)

        !uc(1:ngz-1,i1_u:i2_u,j1_u+2) = myid+13.d0

        ! this is to catch the boundary condtion at x=0 (uc start from i=2 in global grid)
        ! needs some further discusion
        if (icoord(1,myid+1)==1) then
            i1 = i1_u-1
            !uc(1:ngz-1,i1:i1,j1_u+2) = myid+66.d0
        endif
        
        jo = j1_u + jmax_ovr-jbmap_1(j1_u)-2 
        !write(0,*)'cfd socket, jo:',jo
        call coupler_send_grid_data(uc(1:ngz-1,i1:i2_u,jo:jo),index_transpose=(/2,3,1/))


        !do i=i1, i2_u
        !    write(600+myid,'(1000(E11.3,1x))') uc(1:ngz-1,i,j1_u+2)
        !enddo
        ! write(600+myid,'(1000(E11.3,1x))')
        ! call flush(600+myid)

    end subroutine socket_coupler_send_velocity


    subroutine  socket_coupler_get_md_BC(uc,vc,wc)
        use data, only : ngz,i1_u,i2_u,i1_v,i2_v,i1_w,i2_w, ibmap_1
        use messenger, only : icomm_grid, icoord
        use coupler
        implicit none

        real(kind(0.d0)), intent(inout) :: uc(0:,0:,0:),vc(0:,0:,0:),wc(0:,0:,0:) 
        integer j
        integer, save :: i1g_u,i1g_v,i1g_w,i1_ul,myid,ierr
        logical, save :: firsttime = .true.

        real(kind(0.d0)), allocatable, dimension(:,:,:,:) :: ucbuff,vcbuff,wcbuff

        call mpi_comm_rank(icomm_grid,myid,ierr)

        if ( firsttime) then
            firsttime = .false.
            i1g_u = ibmap_1(i1_u)
            i1g_v = ibmap_1(i1_v)
            i1g_w = ibmap_1(i1_w)
            !write(0,*) 'i indices ngz, i1_u,i2_u,i1_v,i2_v,i1_w,i2_w ', ngz, i1_u,i2_u,i1_v,i2_v,i1_w,i2_w 
            
            !write(0,*) 'size uc ...', size(uc,2), size(vc,2), size(wc,2)

            !write(0,*) 'global indices i1g_u, i1g_v, i1g_w ',  i1g_u, i1g_v, i1g_w
            
            ! this is to catch the boundary condtion at x=0 (uc start from i=2 in global grid)
            ! needs some further discusion
            i1_ul = i1_u
            if (icoord(1,myid+1)==1) then
                i1_ul = i1_u-1
            endif

         endif
		!print'(2a,2i8,4f25.16)', 'CFD befr data',code_name(COUPLER_REALM), myid, & 
		!							size(uc(:,:,0)), maxval(uc(:,:,0)),minval(uc(:,:,0)),sum(uc(:,:,0)),uc(3,3,0)
		!print*, 'cfd b4', ngz-1,i1_ul,i2_u
        call coupler_recv_grid_data(uc(1:ngz-1,i1_ul:i2_u,0:0),index_transpose=(/2,3,1/),accumulate=.true.,pbc=1)
		!print'(2a,2i8,4f25.16)', 'CFD recv data',code_name(COUPLER_REALM), myid, & 
		!							size(uc(:,:,0)), maxval(uc(:,:,0)),minval(uc(:,:,0)),sum(uc(:,:,0)),uc(3,3,0)
        call coupler_recv_grid_data(vc(1:ngz-1,i1_v:i2_v,0:1),index_transpose=(/2,3,1/),accumulate=.true.)
        call coupler_recv_grid_data(wc(1:ngz,i1_w:i2_w,0:0),index_transpose=(/2,3,1/),accumulate=.true.,pbc=3)

        !debug writes
        !do j=i1_ul, i2_u
        !    write(200+myid,'(1000(E11.3,1x))') uc(1:ngz-1,j,0)
        !enddo
        !write(200+myid,*) 
        !call flush(200+myid)
        !do j=i1_v,i2_v
        !    write(300+myid,'(1000(E11.3,1x))') vc(1:ngz-1,j,0)
        !enddo
        ! write(300+myid,*)
        !do j=i1_v,i2_v
        !    write(400+myid,'((1000(E11.3,1x)))') vc(1:ngz-1,j,1)
        !enddo
        ! write(400+myid,*)
        !do j=i1_w,i2_w
        !    write(500+myid,'((1000(E11.3,1x)))') wc(1:ngz,j,0)
        !enddo
        !write(500+myid,*)

    end subroutine socket_coupler_get_md_BC

#endif
end module continuum_coupler_socket
