!=============================================================================
!                                   mesh_gen
!
! Generate numerical Mesh used in simulation
!
!-----------------------------------------------------------------------------

module module_mesh_gen

	use computational_constants
	use grid_arrays

end module module_mesh_gen
!----------------------------------------------------------------------------------

subroutine continuum_mesh_gen
	use module_mesh_gen
	implicit none

	integer 		::  i , j

	!domain_volume = (nx*ny) / grid_density
	!lx = domain_volume**(nx/dble((nx+ny)))
	!ly = domain_volume**(ny/dble((nx+ny)))	

	do i = 1, nx+3
		mx(i) = ((i-1)*(lx + 2 * lx/nx))/(nx+2)
	enddo

	do j = 1, ny+3	
		my(j) = ((j-1)*(ly + 2 * ly/ny))/(ny+2)
	enddo

	do i = 1, nx+2
		px(i) = (mx(i) + mx(i+1))/ 2.d0
	enddo

	do j = 1, ny+2
		py(j) = (my(j) + my(j+1))/ 2.d0
	enddo

	do i = 1, nx+2
		!Grid spacing
		delta_x(i) = mx(i+1) - mx(i)

		!Surface vector dotted with length to give surface area
		sx(i,1) = 1 * delta_x(i)
		sx(i,2) = 0 * delta_x(i)
		sx(i,3) =-1 * delta_x(i)
		sx(i,4) = 0 * delta_x(i)

	enddo

	do j = 1, ny+2
		!Grid spacing
		delta_y(j) = my(j+1) - my(j)

		!Surface vector dotted with length to give surface area
		sy(j,1) = 0 * delta_y(j)
		sy(j,2) = 1 * delta_y(j)
		sy(j,3) = 0 * delta_y(j)
		sy(j,4) =-1 * delta_y(j)

	enddo

	do i = 1, nx+2
	do j = 1, ny+2

		vcell(i,j) = (delta_x(i) * delta_y(j))
	enddo
	enddo

	domain_volume = lx*ly

end subroutine continuum_mesh_gen
