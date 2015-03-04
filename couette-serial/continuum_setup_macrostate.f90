!=============================================================================
!                               Setup Macrostate
! Setup initial state of continuum solver
!
!-----------------------------------------------------------------------------

module module_setup_macrostate

	use computational_constants
	use grid_arrays

end module module_setup_macrostate
!----------------------------------------------------------------------------------

subroutine continuum_setup_macrostate
	use module_setup_macrostate
	implicit none

	integer	:: i, j
	
	do i = 1, nx+2
	do j = 1, ny+2

		uc(i,j) = 0.d0 
		vc(i,j) = 0.d0

	enddo
	enddo

end subroutine continuum_setup_macrostate
