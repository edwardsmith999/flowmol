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

!----------------------------------------------------------------------------------
!Move molecules and apply constraints or thermostatting as determined by tagging system

subroutine simulation_move_particles_tag
use module_move_particles
use calculated_properties_MD
implicit none

	integer 		:: n, ixyz, thermostatnp
	double precision	:: vel, slice_momentum2, dzeta_dt, massheatbath
	double precision	:: vreduce, relaxfactor, ascale, bscale



	!Check if any molecules are thermostatted and calculate appropriate coefficients
	if (maxval(tag) .ge. 4) then
		v2sum = 0.d0      ! Reset all sums
		thermostatnp = 0
		do n = 1, np    ! Loop over all particles
			!CHECK IF YOU WANT THIS LINE IN!!
			if (tag(n) .lt. 4) cycle	!Only include thermostatted molecules
			!WALL THERMOSTATS MAY NEED TO KNOW DOMAIN TEMPERATURE
			do ixyz = 1, nd    ! Loop over all dimensions
				vel = v(n,ixyz) - 0.5d0*a(n,ixyz)*delta_t
				v2sum = v2sum + vel**2 !Add up all molecules' velocity squared components  
			enddo
			thermostatnp = thermostatnp + 1
		enddo

		massheatbath = thermostatnp * delta_t

		call globalSum(v2sum)	!Obtain global velocity sum

		dzeta_dt = (v2sum - (nd*thermostatnp + 1)*inputtemperature) / massheatbath
		zeta = zeta + delta_t*dzeta_dt

		bscale=1.0/(1.0+0.5*delta_t*zeta)
		ascale=(1-0.5*delta_t*zeta)*bscale
	endif


	!Step through each particle n
	do n = 1,np        

		select case (tag(n))
		case (0)
			!Leapfrog mean velocity calculated here at v(t+0.5delta_t) = v(t-0.5delta_t) + a*delta_t 
			!Leapfrog mean position calculated here at r(t+delta_t) = r(t) + v(t+0.5delta_t)*delta_t
			v(n,:) = v(n,:) + delta_t * a(n,:) 	!Velocity calculated from acceleration
			r(n,:) = r(n,:) + delta_t * v(n,:)	!Position calculated from velocity
		case (1)
			!Fixed Molecules - no movement r(n+1) = r(n)
		case (2)
			!Fixed with constant sliding speed
			r(n,:) = r(n,:) + delta_t*slidev(n,:)	!Position calculated from velocity
		case (3)
			!Tethered molecules
			call tether_force(n)
			v(n,:) = v(n,:) + delta_t * a(n,:) 	!Velocity calculated from acceleration
			r(n,:) = r(n,:) + delta_t * v(n,:)	!Position calculated from velocity
		case (4)
			!Nose Hoover Thermostatted Molecule
	        	!v(n,1) = v(n,1)*ascale + a(n,1)*delta_t*bscale
			!r(n,1) = r(n,1)    +     v(n,1)*delta_t			
	        	!v(n,2) = v(n,2)*ascale + a(n,2)*delta_t*bscale
			!r(n,2) = r(n,2)    + 	 v(n,2)*delta_t				
	        	!v(n,3) = v(n,3)*ascale + a(n,3)*delta_t*bscale
			!r(n,3) = r(n,3)    +     v(n,3)*delta_t	

			!Nose Hoover Thermostatted Molecule z direction only
			v(n,1) = v(n,1) + delta_t * a(n,1) 	!Velocity calculated from acceleration
			r(n,1) = r(n,1) + delta_t * v(n,1)	!Position calculated from velocity
			v(n,2) = v(n,2) + delta_t * a(n,2) 	!Velocity calculated from acceleration
			r(n,2) = r(n,2) + delta_t * v(n,2)	!Position calculated from velocity
	        	v(n,3) = v(n,3)*ascale + a(n,3)*delta_t*bscale
			r(n,3) = r(n,3)    +     v(n,3)*delta_t	

		case (5)
			!Thermostatted Tethered molecules unfixed with no sliding velocity
			call tether_force(n)
	        	v(n,1) = v(n,1)*ascale + a(n,1)*delta_t*bscale
			r(n,1) = r(n,1)    +     v(n,1)*delta_t			
	        	v(n,2) = v(n,2)*ascale + a(n,2)*delta_t*bscale
			r(n,2) = r(n,2)    + 	 v(n,2)*delta_t				
	        	v(n,3) = v(n,3)*ascale + a(n,3)*delta_t*bscale
			r(n,3) = r(n,3)    +     v(n,3)*delta_t	
		case (6)
			!Tethered molecules with sliding velocity
			call tether_force(n)
			v(n,:) = v(n,:) + delta_t * a(n,:) 	!Velocity calculated from acceleration
			r(n,:) = r(n,:) + delta_t * v(n,:) + delta_t*slidev(n,:)	!Position calculated from velocity+slidevel
		case (7)
			!Thermostatted Tethered molecules unfixed with sliding velocity
			call tether_force(n)
	        	v(n,1) = v(n,1)*ascale + a(n,1)*delta_t*bscale
			r(n,1) = r(n,1)    +     v(n,1)*delta_t	+ slidev(n,1)*delta_t		
	        	v(n,2) = v(n,2)*ascale + a(n,2)*delta_t*bscale
			r(n,2) = r(n,2)    + 	 v(n,2)*delta_t	+ slidev(n,2)*delta_t			
	        	v(n,3) = v(n,3)*ascale + a(n,3)*delta_t*bscale
			r(n,3) = r(n,3)    +     v(n,3)*delta_t	+ slidev(n,3)*delta_t

		case default
			stop "Invalid molecular Tag"
		end select
	enddo


end subroutine simulation_move_particles_tag
