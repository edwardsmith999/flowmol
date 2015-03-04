!=============================================================================
!                                   Finish
! Write final state and close file 
!
!-----------------------------------------------------------------------------

module module_finish

	use computational_constants
	use grid_arrays

end module module_finish
!----------------------------------------------------------------------------------

subroutine continuum_finish
	use module_advance_time
	implicit none

	deallocate(mx)
	deallocate(my)
	deallocate(px)
	deallocate(py)
	deallocate(delta_x)
	deallocate(delta_y)
	deallocate(sx)
	deallocate(sy)

	deallocate(solid)
	deallocate(uc)
	deallocate(vc)
	deallocate(uc_t_minus_1)
	deallocate(vc_t_minus_1)
	deallocate(xresidual)
	deallocate(yresidual)
	deallocate(xresidual_t_minus_1)
	deallocate(yresidual_t_minus_1)
	deallocate(flux_xx)
	deallocate(flux_xy)
	deallocate(flux_yx)
	deallocate(flux_yy)
	deallocate(vcell)

	deallocate(sflux_xx)
	deallocate(sflux_xy)
	deallocate(sflux_yx)
	deallocate(sflux_yy)

	deallocate(tau_xx)
	deallocate(tau_xy)
	deallocate(tau_yx)
	deallocate(tau_yy)

	deallocate(xfluxdots)
	deallocate(yfluxdots)


end subroutine continuum_finish
