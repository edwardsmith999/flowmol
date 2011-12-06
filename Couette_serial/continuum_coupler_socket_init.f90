module continuum_coupler_socket_init
        implicit none

        logical, parameter :: use_coupling = .true.
	! attention nx, ny are the number of cell in the fluid domain !!!

	! we need a nz and z array for MD even CFD is 2D. 
	! the simple solution is to use a small cell in z direction
	! but a larger domain can be used to get better averages and
	! 2d MD topology. 

        integer :: nz = 4 			!test, go as x direction
        real(kind(0.d0)), allocatable :: z(:)

contains

         subroutine init_coupler(comm)
                use coupler
                implicit none
                
                integer, intent(out) :: comm
                integer ierr

                call coupler_create_comm(COUPLER_CFD,comm,ierr)

                allocate(z(nz+1)) ! nz is the number of cells, to keep with the other dimensions

        end subroutine init_coupler

        subroutine continuum_coupler_init 
                use computational_constants, only : nx, ny, nsteps => continuum_Nsteps,&
                        continuum_delta_t
                use  grid_arrays, only : mx, my
                use continuum_data_export, only : npx,npy,npz,icoord
                use coupler, only : coupler_get_cfd_info, coupler_create_map 
                implicit none

!                integer kmino, kmin, kmax, kmaxo
                real(kind(0.d0))  dz
                
                !z(1)    =  0.d0    ! thincknes of MD simulation
                !z(2) =  54.7189 ! 34.199518933533936d0
                
                z(:) = mx(1:nz+1) ! assumes nz < nx !!!
                
                ! nsteps = nsteps+1 for the intialisation step in setup_continuum

                call coupler_get_cfd_info(imino=1,imin=2,imax=nx+2,&
                        imaxo=nx+3,jmino=1,jmin=2,jmax=ny+2,jmaxo=ny+3,&
                        kmino=1,kmin=1,kmax=nz+1,kmaxo=nz+1,nsteps=nsteps+1,&
                        x=mx,y=my,z=z,dx=mx(2)-mx(1),dz=z(2)-z(1),npx=npx,npy=npy,npz=npz,&
                        icoord=icoord,dt=continuum_delta_t)

                call coupler_create_map
                
        end subroutine continuum_coupler_init

end module continuum_coupler_socket_init


