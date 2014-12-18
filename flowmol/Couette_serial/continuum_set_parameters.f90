!=============================================================================
!                                   Set parameters
! Setup values used for the simulation 
!
!-----------------------------------------------------------------------------

module module_continuum_set_parameters
	
	use physical_constants
	use computational_constants
	use grid_arrays

end module module_continuum_set_parameters
!----------------------------------------------------------------------------------

subroutine continuum_set_parameters
	use module_continuum_set_parameters
	!use calculated_properties_MD
	implicit none

	integer			:: i, j
	double precision	:: U

	!Add one to each to include halo cells
	allocate(mx(nx+3))
	allocate(my(ny+3))
	allocate(px(nx+2))
	allocate(py(ny+2))
	allocate(delta_x(nx+2))
	allocate(delta_y(ny+2))
	allocate(sx(nx+2,4))
	allocate(sy(ny+2,4))

	allocate(solid(nx,ny))
	allocate(uc(nx+2,ny+2))
	allocate(vc(nx+2,ny+2))
	allocate(uc_t_minus_1(nx+2,ny+2))
	allocate(vc_t_minus_1(nx+2,ny+2))
	allocate(xresidual(nx+2,ny+2))
	allocate(yresidual(nx+2,ny+2))
	allocate(xresidual_t_minus_1(nx+2,ny+2))
	allocate(yresidual_t_minus_1(nx+2,ny+2))
	allocate(flux_xx(nx+2,ny+2))
	allocate(flux_xy(nx+2,ny+2))
	allocate(flux_yx(nx+2,ny+2))
	allocate(flux_yy(nx+2,ny+2))
	allocate(vcell(nx+2,ny+2))

	allocate(sflux_xx(nx+2,ny+2,4))
	allocate(sflux_xy(nx+2,ny+2,4))
	allocate(sflux_yx(nx+2,ny+2,4))
	allocate(sflux_yy(nx+2,ny+2,4))

	allocate(tau_xx(nx+2,ny+2,4))
	allocate(tau_xy(nx+2,ny+2,4))
	allocate(tau_yx(nx+2,ny+2,4))
	allocate(tau_yy(nx+2,ny+2,4))

	allocate(xfluxdots(nx+2,ny+2,4))
	allocate(yfluxdots(nx+2,ny+2,4))

	!Set initial step of simulation
	continuum_initialstep = 0

	!Establish Reynolds Number
	!Length scale of problem is the molecular length scale and so L = 1!
	U = maxval(uwall_BC)
	Re = (rho * 1.d0 * U) / meu
	
	!Routine to inport Mesh from file
	!open(4,file='mesh_gen/mesh',form='binary')
	!do i=1,nx
	!do j=1,ny
	!	read(4) solid(i,j)
	!enddo
	!enddo
	!close(4,status='keep')

end subroutine continuum_set_parameters
