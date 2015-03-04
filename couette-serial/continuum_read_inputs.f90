!=============================================================================
!                                   inputs
! Read Input values used to set up the simulation 
!
!-----------------------------------------------------------------------------

module module_continuum_inputs

	use physical_constants
	use computational_constants
	use grid_arrays
	use continuum_data_export, only : prefix_dir

end module module_continuum_inputs
!----------------------------------------------------------------------------------

subroutine continuum_read_inputs
	use module_continuum_inputs
	implicit none

	!double precision	:: i	!dummy variable as lx and ly are define by coupling

	open(1001,file=trim(prefix_dir)//'input_continuum')

	!Input computational co-efficients
	read(1001,* ) continuum_Nsteps		!Number of simulation steps
	read(1001,* ) continuum_tplot		!Output ever number of steps
	read(1001,* ) nx					!Number of cells in Domain in x - nx
	read(1001,* ) ny          			!Number of cells in Domain in y - ny
	read(1001,* ) nz          			!Number of cells in Domain in z - nz
	read(1001,* ) lx		       		!Domain size in x - lx
	read(1001,* ) ly		       		!Domain size in y - ly
	read(1001,* ) lz		       		!DUMMY Domain size in z - lz
	read(1001,* ) continuum_delta_t  	!Timestep
	read(1001,* ) meu					!Viscosity
	read(1001,* ) rho					!Density
	read(1001,* ) solver				!Choice of solver
	read(1001,* ) continuum_vflag		!OUTPUT FLAG - Velocity record 0=off 1,2=x,y slice 3=3D bins

	!Read boundary conditions
	read(1001,* ) BC_flag(1)			!Right wall flag  0=periodic, 1=Dirichlet, 2=Neumann
	read(1001,* ) uwall_BC(1)			!Right wall Boundary value (not used for periodic)
	read(1001,* ) BC_flag(2)			!Top wall flag    0=periodic, 1=Dirichlet, 2=Neumann
	read(1001,* ) uwall_BC(2)			!Top wall Boundary value (not used for periodic)
	read(1001,* ) BC_flag(3)			!Left wall flag   0=periodic, 1=Dirichlet, 2=Neumann 
	read(1001,* ) uwall_BC(3)			!Left wall Boundary value (not used for periodic) 
	read(1001,* ) BC_flag(4)			!Bottom wall flag 0=periodic, 1=Dirichlet, 2=Neumann
	read(1001,* ) uwall_BC(4)			!Bottom wall Boundary value (not used for periodic)

	close(1001,status='keep')      !Close input file

end subroutine continuum_read_inputs
