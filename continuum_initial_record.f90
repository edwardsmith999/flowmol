!=============================================================================
!                               Setup initial record
! Print initial output for continuum simulation 
!
!-----------------------------------------------------------------------------

module module_continuum_initial_record

	use physical_constants
	use computational_constants
	use grid_arrays
	!use calculated_properties_MD
	!use computational_constants_MD

end module module_continuum_initial_record
!----------------------------------------------------------------------------------

subroutine continuum_initial_record
	use module_continuum_initial_record
	implicit none	

	print*, '=================Continuum Simulation Parameters======================'
	print'(3(a,f10.5))', 'Specified Size in x ', lx , ' & y ', ly, ' Volume of domain ', domain_volume
	print'(2(a,f10.5))', 'Resulting computaional domain size in x ', mx(nx+2)-mx(2),' & y ',my(ny+2)-my(2)
	print'(2(a,i8))', 'Number of points in x ', nx, ' Number of points in y ', ny
	print'(2(a,f10.5))','Cell size in x ', delta_x(1), ' Cell size in y ', delta_y(1)
	print*, 'number of timesteps', continuum_Nsteps
	print*, 'frequency of plots', continuum_tplot
	print*, 'meu', meu
	print*, 'rho', rho

	print*, 'continuum Reynolds number =', Re
	print*, '======================================================================='

end subroutine continuum_initial_record
