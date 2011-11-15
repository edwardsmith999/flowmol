module continuum_coupler_socket_init
        implicit none

        logical, parameter :: use_coupling = .true.
! attention nx, ny are the number of cell in the fluid domain !!!

contains

         subroutine create_communicators(comm)
                use coupler
                implicit none
                
                integer, intent(out) :: comm
                integer ierr

                call coupler_create_comm(COUPLER_CFD,ierr)
                comm = COUPLER_COMM

       end subroutine create_communicators

        subroutine continuum_coupler_init 
                use computational_constants, only : nx, ny, nsteps => continuum_Nsteps,&
                        continuum_delta_t
                use  grid_arrays, only : mx, my
                use continuum_data_export, only : npx,npy,npz,icoord
                use coupler, only : coupler_get_cfd_info, coupler_create_map 
                implicit none

!                integer kmino, kmin, kmax, kmaxo
                real(kind(0.d0)) z(2), dz
                
                z(1) =  0.d0   ! thincknes of MD simulation
                z(2) = 34.199518933533936d0

                call coupler_get_cfd_info(imino=1,imin=2,imax=nx+2,&
                        imaxo=nx+3,jmino=1,jmin=2,jmax=ny+2,jmaxo=nx+3,&
                        kmino=1,kmin=1,kmax=2,kmaxo=2,nsteps=nsteps,&
                        x=mx,y=my,z=z,dx=mx(2)-mx(1),dz=z(2)-z(1),npx=npx,npy=npy,npz=npz,&
                        icoord=icoord,dt=continuum_delta_t)

                call coupler_create_map
                
        end subroutine continuum_coupler_init

end module continuum_coupler_socket_init


