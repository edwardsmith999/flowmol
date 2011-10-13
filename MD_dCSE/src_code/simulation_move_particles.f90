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


!----------------------------------------------------------------------------------
! SLLOD_move based on the extensive literature from Hoover, Evans etc
! The shear rate specifies the velocity applied at the top point of the domain.
! SLLOD force is also applied to fixed atoms ~ they move at a slower speed
! than free liquid atoms as there is no effect of bulk flow. 

!A Nose Hoover thermostat is also applied


subroutine SLLOD_move
use module_move_particles
use calculated_properties_MD
implicit none

	integer			:: n, ixyz
	double precision	:: vel, slice_momentum2, dzeta_dt, massheatbath
	double precision	:: vreduce, relaxfactor
	double precision	:: shear_rate

	v2sum = 0.d0      ! Reset all sums
	massheatbath = globalnp * delta_t

	!Temperature used is of the wall atoms only
	do n = 1, np    ! Loop over all particles
	do ixyz = 1, nd    ! Loop over all dimensions
		vel = (v(n,ixyz) - 0.5d0*a(n,ixyz)*delta_t )*thermostat(n,ixyz)
		v2sum = v2sum + vel**2 !Add up all molecules' velocity squared components  
	enddo
	enddo

	!print*, 'wall temp', v2sum/ (real(nd,kind(0.d0))*real(globalnp,kind(0.d0)))

	call globalSum(v2sum)	!Obtain global velocity sum

	dzeta_dt = (v2sum - (nd*globalnp + 1)*inputtemperature) / massheatbath
	zeta = zeta + delta_t*dzeta_dt

	shear_rate= 0.0005d0/domain(2)
	
	do n = 1,np        !Step through each particle n

		!Check for tethering force and correct applied force accordingly
		if (tag(n).eq. 3) then
			call tether_force(n)
		endif

		!Velocity calculated from acceleration
		v(n,1) =(v(n,1) + delta_t*a(n,1))*fix(n,1)		&	!Fixed Molecules ~> a=0
				+ shear_rate*(r(n,2)+halfdomain(2))*(1-thermostat(n,1))	&	!SLLOD force on x direction
				+ slidev(n,1)				&	!Add sliding velocity
				- zeta*v(n,1)*delta_t*thermostat(n,1)! *         &	!Thermostat

		!Position calculated from velocity
		r(n,1) = r(n,1) + delta_t*v(n,1)

		!print*, shear_rate*(r(n,2)+halfdomain(2))

		!Velocity calculated from acceleration
		v(n,2) =(v(n,2) + delta_t*a(n,2))*fix(n,2)	  &	!Fixed Molecules ~> a=0
				+ slidev(n,2)			  &	!Add sliding velocity
				- zeta*v(n,2)*delta_t*thermostat(n,2)! *   &	!Thermostat

		!Position calculated from velocity
		r(n,2) = r(n,2) + delta_t*v(n,2)				

		!Velocity calculated from acceleration
		v(n,3) =(v(n,3) + delta_t*a(n,3))*fix(n,3) 	  &	!Fixed Molecules ~> a=0
				+ slidev(n,3)			  &	!Add sliding velocity
				- zeta*v(n,3)*delta_t*thermostat(n,3)! *   &	!Thermostat

		!Position calculated from velocity
		r(n,3) = r(n,3) + delta_t*v(n,3)				
		
	enddo
	
end subroutine SLLOD_move


!----------------------------------------------------------------------------------
! SLLOD_move based on the extensive literature from Hoover, Evans etc
! The shear rate specifies the velocity applied at the top point of the domain.
! SLLOD force is also applied to fixed atoms ~ they move at a slower speed
! than free liquid atoms as there is no effect of bulk flow. 

!A Profile Unbias [nose hoover] Thermostat (PUT) has also been used


subroutine SLLOD_move_PUT
use module_move_particles
use calculated_properties_MD
implicit none

	integer				:: n, ixyz
	integer				:: ibin, jbin, kbin
	double precision		:: vel, slice_momentum2, dzeta_dt, massheatbath
	double precision		:: vreduce, relaxfactor
	double precision		:: shear_rate
	double precision,dimension(3)	:: vmean, slicebinsize
	
	v2sum = 0.d0      ! Reset all sums
	massheatbath = globalnp * delta_t

	slicebinsize(:) = domain(:) / nbins(:)

	!Velocity sum for fluctuations only so remove mean flow
	do n = 1, np    ! Loop over all particles

		!Determine bins using integer division
		ibin = ceiling((r(n,1)+halfdomain(1))/slicebinsize(1)) !Establish current bin
		if (ibin > nbins(1)) ibin = nbins(1) 		!Prevents out of range values
		if (ibin < 1 ) ibin = 1        		!Prevents out of range values
		jbin = ceiling((r(n,2)+halfdomain(2))/slicebinsize(2)) !Establish current bin
		if (jbin > nbins(2)) jbin = nbins(2) 		!Prevents out of range values
		if (jbin < 1 ) jbin = 1        		!Prevents out of range values
		kbin = ceiling((r(n,3)+halfdomain(3))/slicebinsize(3)) !Establish current bin
		if (kbin > nbins(3)) kbin = nbins(3) 		!Prevents out of range values
		if (kbin < 1 ) kbin = 1        		!Prevents out of range values

		!Establish average velocity of current cell
		vmean(:) = slice_momentumbin(ibin,jbin,kbin,:)/(slice_massbin(ibin,jbin,kbin)+1)

		do ixyz = 1, nd    ! Loop over all dimensions
			vel = v(n,ixyz) - 0.5d0*a(n,ixyz)*delta_t - vmean(ixyz)
			v2sum = v2sum + vel**2 !Add up all molecules' velocity squared components  
		enddo

	enddo

	call globalSum(v2sum)	!Obtain global velocity sum

	dzeta_dt = (v2sum - (nd*globalnp + 1)*inputtemperature) / massheatbath
	zeta = zeta + delta_t*dzeta_dt

	shear_rate= 0.005d0/domain(2)
	
	do n = 1,np        !Step through each particle n

		!Check for tethering force and correct applied force accordingly
		if (tag(n).eq. 3) then
			call tether_force(n)
		endif

		!Velocity calculated from acceleration
		v(n,1) =(v(n,1) + delta_t*a(n,1))*fix(n,1)		&	!Fixed Molecules ~> a=0
				+ shear_rate*(r(n,2)+halfdomain(2))	&	!SLLOD force on x direction
				+ slidev(n,1)				&	!Add sliding velocity
				- zeta*v(n,1)*delta_t				!Thermostat

		!Position calculated from velocity
		r(n,1) = r(n,1) + delta_t*v(n,1)

		!Velocity calculated from acceleration
		v(n,2) =(v(n,2) + delta_t*a(n,2))*fix(n,2)	  &	!Fixed Molecules ~> a=0
				+ slidev(n,2)			  &	!Add sliding velocity
				- zeta*v(n,2)*delta_t			!Thermostat

		!Position calculated from velocity
		r(n,2) = r(n,2) + delta_t*v(n,2)				

		!Velocity calculated from acceleration
		v(n,3) =(v(n,3) + delta_t*a(n,3))*fix(n,3) 	  &	!Fixed Molecules ~> a=0
				+ slidev(n,3)			  &	!Add sliding velocity
				- zeta*v(n,3)*delta_t			!Thermostat

		!Position calculated from velocity
		r(n,3) = r(n,3) + delta_t*v(n,3)				
		
	enddo
	
end subroutine SLLOD_move_PUT

