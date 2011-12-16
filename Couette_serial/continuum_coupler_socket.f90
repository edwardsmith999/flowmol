!=============================================================================
!				   
! Routine which interface with the coupler to the CFD code
! Corresponding empty shell dummy routine used for uncoupled calculation
!
!  Lucian Anton, November 2011
!
!-----------------------------------------------------------------------------

module continuum_coupler_socket
	implicit none

contains

!=============================================================================
!  call coupler routine to build CFD communicator                                  
!-----------------------------------------------------------------------------

subroutine init_coupler(comm)
	use coupler
	implicit none
	
	integer, intent(out) :: comm
	integer ierr

	call coupler_create_comm(COUPLER_CFD,comm,ierr)

end subroutine init_coupler

!=============================================================================
! Call coupler routine to create map from CFD to MD                          
!-----------------------------------------------------------------------------

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
		icoord=icoord,dt=continuum_delta_t,density=rho)
	call coupler_create_map
	
end subroutine continuum_coupler_init

!=============================================================================
!  Package velocities and pass to coupler to send to MD                                
!-----------------------------------------------------------------------------

subroutine send_CFDvel
	use grid_arrays, only : uc, vc
	use  coupler, only : coupler_cfd_send_velocity
	implicit none

		integer i
	real(kind(0.d0)) uc3d(1,size(uc,dim=1),size(uc,dim=2)), & !hack for 2d parallelism, z dimension independent of nx
			 vc3d(1,size(vc,dim=1),size(vc,dim=2))

	do i=1,size(uc3d,dim=1)
		uc3d(i,:,:) = uc (:,:)
	enddo
	do i=1,size(vc3d,dim=1)
		vc3d(i,:,:) = vc(:,:)
	enddo

	call coupler_cfd_send_velocity(uc3d,vc3d)

end subroutine send_CFDvel

!=============================================================================
!  Get MD velocities from coupler and unpack for CFD boundary conditions                                   
!-----------------------------------------------------------------------------

subroutine MD_continuum_BC(u,v)
	use coupler, only : coupler_cfd_get_velocity
	implicit none

	real(kind(0.d0)), intent(out) :: u(:), v(:)
	real(kind(0.d0)) w(1,1,1), u3(1,size(u),1),v3(1,size(v),1)

	!A polymorphic interface is needed here
	call coupler_cfd_get_velocity(u3,v3,w)

	u(:) =  u3(1,:,1)
	v(:) =  v3(1,:,1)

end subroutine MD_continuum_BC
end module continuum_coupler_socket

