!=============================================================================
!                               Advance_time
! Calculate next timestep using fluxes from previous
!
!-----------------------------------------------------------------------------

module module_advance_time

	use physical_constants
	use computational_constants
	use grid_arrays

end module module_advance_time
!----------------------------------------------------------------------------------

!================
!=Explicit Euler=
!================

subroutine continuum_advance_time
	use module_advance_time
	implicit none

	integer	:: i, j

	!Apply a force to all cells
	!call continuum_apply_force(xresidual)

	!Store previous timesteps
	uc_t_minus_1 = uc

	!Advance Time
	do i = 2, nx+1
	do j = 2, ny+1

		!Finite difference
		!uc(i,j) = uc(i,j) + continuum_delta_t * xresidual(i,j)

		!Finite volume
		uc(i,j) = uc(i,j) + continuum_delta_t * xresidual(i,j) / vcell(i,j)
		vc(i,j) = vc(i,j) + continuum_delta_t * yresidual(i,j) / vcell(i,j)
	
	enddo
	enddo

end subroutine continuum_advance_time


!=======================
!=Four Step Runge Kutta=
!=======================
subroutine continuum_advance_time_RK
	use module_advance_time
	implicit none

	integer	:: i, j
	double precision :: alpha_1, alpha_2, alpha_3, alpha_4
	double precision :: beta_1, beta_2, beta_3, beta_4
	double precision, dimension(nx+2,ny+2):: uc_2, vc_2, xres_2, yres_2
	double precision, dimension(nx+2,ny+2):: uc_3, vc_3, xres_3, yres_3
	double precision, dimension(nx+2,ny+2):: uc_4, vc_4, xres_4, yres_4

	!Define Co-efficients
	alpha_1 = 0.5d0
	alpha_2 = 0.5d0
	alpha_3 = 1.d0
	beta_1 = 1.d0/6.d0
	beta_2 = 1.d0/3.d0
	beta_3 = 1.d0/3.d0
	beta_4 = 1.d0/6.d0

	!call continuum_apply_force(xresidual)

	!Calculate second step

	do i = 2, nx+1
	do j = 2, ny+1
		uc_2(i,j) = uc(i,j) + (continuum_delta_t * alpha_1 * xresidual(i,j))
		vc_2(i,j) = vc(i,j) + (continuum_delta_t * alpha_1 * yresidual(i,j))
	enddo
	enddo
	call continuum_calculate_residual(uc_2, vc_2, xres_2, yres_2)
	!call continuum_apply_force(xres_2)

	!Calculate third step

	do i = 2, nx+1
	do j = 2, ny+1
		uc_3(i,j) = uc(i,j) + (continuum_delta_t * alpha_2 *  xres_2(i,j))
		vc_3(i,j) = vc(i,j) + (continuum_delta_t * alpha_2 *  yres_2(i,j))
	enddo
	enddo		
	call continuum_calculate_residual(uc_3, vc_3, xres_3, yres_3)
	!call continuum_apply_force(xres_3)

	!Calculate fourth step

	do i = 2, nx+1
	do j = 2, ny+1
		uc_4(i,j) = uc(i,j) + (continuum_delta_t * alpha_3 *  xres_3(i,j))
		vc_4(i,j) = vc(i,j) + (continuum_delta_t * alpha_3 *  yres_3(i,j))
	enddo
	enddo
	call continuum_calculate_residual(uc_4, vc_4, xres_4, yres_4)
	!call continuum_apply_force(xres_4)

	!Advance velocity one timestep
	do i = 2, nx+1
	do j = 2, ny+1
		uc(i,j) = uc(i,j) + (continuum_delta_t *(beta_1 * xresidual(i,j) + & 
							 beta_2 * xres_2(i,j)    + &
					       		 beta_3 * xres_3(i,j)    + &
							 beta_4 * xres_4(i,j)))
		vc(i,j) = vc(i,j) + (continuum_delta_t *(beta_1 * yresidual(i,j) + &
							 beta_2 * yres_2(i,j)    + &
							 beta_3 * yres_3(i,j)    + &
							 beta_4 * yres_4(i,j)))
	enddo
	enddo

end subroutine continuum_advance_time_RK

!---------------------------------------------------------------------------------

subroutine continuum_apply_force(xres_t)	
	use module_advance_time
	implicit none

	integer	:: i, j
	double precision, dimension(nx+2,ny+2)	:: xres_t

	!Apply a force to all cells
	do i = 2, nx+1
	do j = 2, ny+1
		xres_t(i,j) = xres_t(i,j) + 0.1d0
				
	enddo
	enddo

end subroutine continuum_apply_force
