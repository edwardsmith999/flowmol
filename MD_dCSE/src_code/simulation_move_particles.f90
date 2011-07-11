!-----------------------------------------------------------------------------
!
!                                Move Particles
! Move particles as a result of forces using the verlet time integration
!
!-----------------------------------------------------------------------------

module module_move_particles

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD

end module module_move_particles
!----------------------------------------------------------------------------------
!Move molecules using leapfrog algorithm with no thermostats. Arrays of integers
!are used to fix or add sliding velocities to molecules based on initial tagging process

subroutine simulation_move_particles
use module_move_particles
implicit none

	integer :: n, ixyz

	do ixyz = 1,nd        !Step through each dimension ixyz
	do n = 1,np        	!Step through each particle n

		!Check for tethering force and correct applied force accordingly
		if (tag(n).eq. 3) then
			call tether_force(n)
		endif

		!Leapfrog mean velocity calculated here at v(t+0.5delta_t) = v(t-0.5delta_t) + a*delta_t 
		!Leapfrog mean position calculated here at r(t+delta_t) = r(t) + v(t+0.5delta_t)*delta_t
		v(n,ixyz) =(v(n,ixyz) + delta_t*a(n,ixyz))*fix(n,ixyz) 	!Velocity calculated from acceleration
		v(n,ixyz) = v(n,ixyz) + slidev(n,ixyz)			!Add velocity of sliding wall
		r(n,ixyz) = r(n,ixyz) + delta_t*v(n,ixyz)		!Position calculated from velocity

	enddo
	enddo

end subroutine simulation_move_particles
