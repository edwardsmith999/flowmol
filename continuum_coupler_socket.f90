module continuum_coupler_socket
        implicit none
if USE_COUPLER
contains

    subroutine continuum_coupler_adjust_domain
        use coupler
        use computational_constants, only : lx,ly,nx,ny
        use physical_constants,      only : rho
        implicit none
        
        call coupler_cfd_adjust_domain(xL=lx,yL=ly,nx=nx,ny=ny,density_cfd=rho)
        
    end subroutine continuum_coupler_adjust_domain
    
    
    subroutine continuum_coupler_init 
        use computational_constants, only : nx, ny, nsteps => continuum_Nsteps,&
            continuum_delta_t
        use physical_constants, only : rho
        use  grid_arrays, only : mx, my
        use continuum_data_export, only : npx,npy,npz,icoord
        use coupler, only : coupler_cfd_init, coupler_create_map 
        use continuum_messenger
        implicit none
        
        real(kind(0.d0)) :: z0(1)=(/0.d0/)
            
        ! nsteps = nsteps+1 for the intialisation step in setup_continuum
        
        ! 2D problem, z direction parameters are set o trivial values
        
        call coupler_cfd_init(icomm_grid,imino=1,imin=2,imax=nx+2,&
            imaxo=nx+3,jmino=1,jmin=2,jmax=ny+2,jmaxo=ny+3,&
            kmino=1,kmin=1,kmax=1,kmaxo=1,nsteps=nsteps+1,&
            x=mx,y=my,z=z0,dx=mx(2)-mx(1),dz=0.d0,npx=npx,npy=npy,npz=npz,&
            icoord=icoord,dt=continuum_delta_t)
        
        call coupler_create_map
        
    end subroutine continuum_coupler_init


!====================================================================================
! sends uc in the layer jmax_overlap-2:jmax_overlap-1 for continuum force calculation
!------------------------------------------------------------------------------------
    subroutine send_CFDvel
        use computational_constants, only : nx
        use grid_arrays, only : uc, vc
        use  coupler
        implicit none
        
        integer i, jmax_ovr
        real(kind(0.d0)) uc3d(1,2:nx+1,1)!, & !hack for 2d parallelism, z dimension independent of nx
                                              !vc3d(1,size(vc,dim=1),size(vc,dim=2))
 
        call coupler_cfd_get(jmax_overlap=jmax_ovr)
  
        do i=2,nx+1
            uc3d(1,i,1) = uc(i,jmax_ovr-2) 
        enddo
        !do i=1,size(vc3d,dim=1)
        !    vc3d(i,:,:) = vc(:,:)
        !enddo

        call coupler_send_grid_data(uc3d,index_transpose=(/2,3,1/))
        
        !do i = 2,nx+1
        !    write(300,*) uc3d(1,i,1)
        !enddo
        !write(300,*)
        
    end subroutine send_CFDvel

!==================================================================================
! receives average velocities between jmin0:jmin and copies the in boundary 
! condition arrays  
!----------------------------------------------------------------------------------
subroutine MD_continuum_BC(u,v)
    use coupler
    implicit none

    real(kind(0.d0)), intent(out) :: u(:), v(:)

    integer i
    !integer, save :: ncall = 0
    real(kind(0.d0))  u3(1,size(u),1),v3(1,size(v),1)


    call coupler_recv_grid_data(u3,index_transpose=(/2,3,1/),accumulate=.true.)
    call coupler_recv_grid_data(v3,index_transpose=(/2,3,1/),accumulate=.true.)

    u(:) =  u3(1,:,1)
    v(:) =  v3(1,:,1)
	!debug 
	!                 u(:) = 0.d0
	!                 v(:) = 0.d0

	!                ncall = ncall + 1
    !do i = 1,size(u)
    !        write(200,*) u(i),v(i)
    !enddo
    !write(200,*)

end subroutine MD_continuum_BC
#endif
end module continuum_coupler_socket

