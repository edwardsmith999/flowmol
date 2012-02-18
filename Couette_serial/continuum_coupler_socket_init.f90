module continuum_coupler_socket_init
        implicit none
#if USE_COUPLER
        logical, parameter :: use_coupling = .true.
	! attention nx, ny are the number of cell in the fluid domain !!!

	! we need a nz and z array for MD even CFD is 2D. 
	! the simple solution is to use a small cell in z direction
	! but a larger domain can be used to get better averages and
	! 2d MD topology. 

contains

         subroutine init_coupler(comm)
                use coupler
                implicit none
                
                integer, intent(out) :: comm
                integer ierr

                call coupler_create_comm(COUPLER_CFD,comm,ierr)

        end subroutine init_coupler


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
                implicit none

                real(kind(0.d0)) :: z0(1)=(/0.d0/)
                
                                
                ! nsteps = nsteps+1 for the intialisation step in setup_continuum

                ! 2D problem, z direction parameters are set o trivial values

                call coupler_cfd_init(imino=1,imin=2,imax=nx+2,&
                        imaxo=nx+3,jmino=1,jmin=2,jmax=ny+2,jmaxo=ny+3,&
                        kmino=1,kmin=1,kmax=1,kmaxo=1,nsteps=nsteps+1,&
                        x=mx,y=my,z=z0,dx=mx(2)-mx(1),dz=0.d0,npx=npx,npy=npy,npz=npz,&
                        icoord=icoord,dt=continuum_delta_t)

                call coupler_create_map
                
        end subroutine continuum_coupler_init

#endif
end module continuum_coupler_socket_init


