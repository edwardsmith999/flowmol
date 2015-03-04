subroutine continuum_CFL
use physical_constants
use computational_constants
implicit none

	double precision	::	CFL

	CFL = continuum_delta_t / (Re*(minval(delta_y)*minval(delta_y)))
	print*, 'CFL number for this case', CFL
	if (CFL .lt. 0.d0) stop "Error in inputs - negative value for delta_t, delta_y or Reynolds number"
	if (CFL .gt. 0.5d0) stop "Input is unstable for diffusion equation with forward Euler scheme - CFL must be less than 0.5"
	if (CFL .gt. 0.695d0) stop "Input is unstable for diffusion equation with 4th order Runge Kutta scheme"

end subroutine continuum_CFL
