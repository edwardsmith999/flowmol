!----------------------------------------------------------------------------------
!                                Initialise Microstate
! Set up position and velocity of each molecule (the microstate)
!
!---------------------------------------------------------------------------------

module module_initialise_microstate

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use calculated_properties_MD

	double precision                  :: angle, rand  !Define variables
	double precision		  :: v13 !Magnitude of v1 and v3 vectors
	
end module module_initialise_microstate
!----------------------------------------------------------------------------------

subroutine setup_initialise_microstate
use module_initialise_microstate
implicit none
	
	if (potential_flag.eq.0) call setup_initialise_microstate_LJ
	if (potential_flag.eq.1) call setup_initialise_microstate_FENE
	if (potential_flag.gt.1) stop 'Potential flag not recognised!'

end subroutine setup_initialise_microstate

!----------------------------------------------------------------------------------
!==================================================================================
subroutine setup_initialise_microstate_LJ
	use module_initialise_microstate
	implicit none

	integer		::	n

	!call setup_initialise_position       	!Setup initial positions
	call setup_initialise_parallel_position !Setup initial positions in //el
	call setup_tag				!Setup location of fixed molecules
	do n = 1,np
		call read_tag(n)		!Read tag and assign properties
	enddo
	call setup_initialise_velocities     	!Setup initial velocities
	!call setup_initialise_velocities_test


end subroutine setup_initialise_microstate_LJ

subroutine setup_initialise_microstate_FENE
	use module_initialise_microstate
	use polymer_info_MD
	implicit none

	integer		::	n

	!call setup_initialise_position       		 !Setup initial positions
	call setup_initialise_position_FENE 		 !Setup initial positions in //el
	call setup_tag								 !Setup location of fixed molecules
	do n = 1,np
		call read_tag(n)						 !Read tag and assign properties
	enddo
	call setup_initialise_velocities     		 !Setup initial velocities
	!call setup_initialise_velocities_test


end subroutine setup_initialise_microstate_FENE
!==================================================================================
!----------------------------------------------------------------------------------
!Initialise Positions
!Set up the intial position of the particles

subroutine setup_initialise_position
	use module_initialise_microstate
	implicit none

	integer :: j, ixyz, n, nx, ny, nz
	double precision, dimension (nd) :: c !Temporary variable

	do ixyz=1,nd
		initialunitsize(ixyz) = domain(ixyz) / initialnunits(ixyz)
	enddo

	!Molecules per unit FCC structure (3D)
	n=1  		!Reset n for start of loop
	do nz=1,initialnunits(3)	!Loop over z column
	c(3) = (nz - 0.75d0)*initialunitsize(3) - halfdomain(3) 
		do ny=1,initialnunits(2)	!Loop over y column
		c(2) = (ny - 0.75d0)*initialunitsize(2) - halfdomain(2) 
			do nx=1,initialnunits(1)	!Loop over all x elements of y column
			c(1) = (nx - 0.75d0)*initialunitsize(1) - halfdomain(1)
				do j=1,4	!4 Molecules per cell
				do ixyz=1,nd	!For all dimensions
					if(j.eq.ixyz .or. j.eq.4) then
						r(n,ixyz)=c(ixyz)
					else
						r(n,ixyz)=c(ixyz)+0.5d0*initialunitsize(ixyz)
					endif
				enddo
				n = n + 1   !Move to next molecule
				enddo
			enddo
		enddo
	enddo

	rinitial = r !Record initial position of all molecules

end subroutine setup_initialise_position

subroutine setup_initialise_position_FENE
	use module_initialise_microstate
	use polymer_info_MD
	implicit none

	integer :: j, ixyz, n, nx, ny, nz
	integer	:: modcheck
	integer	:: chainID, subchainID
	double precision, dimension (nd) :: c !Temporary variable

	do ixyz=1,nd
		initialunitsize(ixyz) = domain(ixyz) / initialnunits(ixyz)
	enddo

	modcheck = 0 + mod(np,chain_length) + mod(4*initialnunits(1),chain_length)
	if (modcheck.ne.0) stop 'Number of molecules must be exactly divisible by &
	& the polymer chain length. Please change the chain length in the input file. &
	& A chain length of 4 should (hopefully) always work.'
	
	!Molecules per unit FCC structure (3D)
	n=1  !Reset n for start of loop
	do nz=1,initialnunits(3) !Loop over z column
	c(3) = (nz - 0.75d0)*initialunitsize(3) - halfdomain(3) 
		do ny=1,initialnunits(2) !Loop over y column
		c(2) = (ny - 0.75d0)*initialunitsize(2) - halfdomain(2) 
			do nx=1,initialnunits(1) !Loop over all x elements of y column
			c(1) = (nx - 0.75d0)*initialunitsize(1) - halfdomain(1)
				do j=1,4 !4 Molecules per cell
					r(n,:) = c(:)
					if (j.eq.2) then
						r(n,1) = c(1) + 0.5d0*initialunitsize(1)
						r(n,3) = c(3) + 0.5d0*initialunitsize(3)
					else if (j.eq.3) then
						r(n,2) = c(2) + 0.5d0*initialunitsize(2)
						r(n,3) = c(3) + 0.5d0*initialunitsize(3)
					else if (j.eq.4) then
						r(n,1) = c(1) + 0.5d0*initialunitsize(1)
						r(n,2) = c(2) + 0.5d0*initialunitsize(2)
					end if
					
					chainID = ceiling(dble(n)/chain_length)	 !Set chain ID of mol n
					subchainID = mod(n,chain_length) !Beads are numbered 1 to chain_length
					if (subchainID.eq.0) subchainID = chain_length !Correct for mod returning 0

					polyinfo_mol(n)%chainID = chainID
 					polyinfo_mol(n)%subchainID = subchainID

					polyinfo_mol(n)%left = n-1
					polyinfo_mol(n)%right= n+1
					if (subchainID.eq.1) polyinfo_mol(n)%left = 0 !Flag for beginning of chain
					if (subchainID.eq.chain_length) polyinfo_mol(n)%right = 0	!Flag for end of chain

					n = n + 1  !Move to next molecule
				enddo
			enddo
		enddo
	enddo

	rinitial = r !Record initial position of all molecules

end subroutine setup_initialise_position_FENE

!----------------------------------------------------------------------------------
!Initialise Positions
!Set up the intial position of the particles

subroutine setup_initialise_parallel_position
	use module_initialise_microstate
        use messenger
	use coupler

	implicit none

	integer 			:: j, ixyz, n, nl, nx, ny, nz
	integer,dimension(nd) 		:: p_units_lb, p_units_ub, nfcc_max
	double precision 		:: CFD_region, removed_height
	double precision, dimension (nd):: rc, c !Temporary variable

        p_units_lb(1) = (iblock-1)*floor(initialnunits(1)/real((npx),kind(0.d0)))
        p_units_ub(1) =  iblock *ceiling(initialnunits(1)/real((npx),kind(0.d0)))
        p_units_lb(2) = (jblock-1)*floor(initialnunits(2)/real((npy),kind(0.d0)))
        p_units_ub(2) =  jblock *ceiling(initialnunits(2)/real((npy),kind(0.d0)))
        p_units_lb(3) = (kblock-1)*floor(initialnunits(3)/real((npz),kind(0.d0)))
        p_units_ub(3) =  kblock *ceiling(initialnunits(3)/real((npz),kind(0.d0)))

	!Set CFD region to top of domain initially
	CFD_region = domain(2)/2.d0
	if (jblock .eq. npy) then
	if ( coupler_is_active ) then
		removed_height = 2*cellsidelength(2)
		CFD_region = domain(2)/2.d0 - removed_height
	endif
	endif

       ! if ( use_coupling ) then
                ! Let the last two cells free at top of the domain in y direction
                ! as this is the region the constrain force acts on
        !        if (jblock == npy) CFD_region = domain(2)/2.d0 - 2*cellsidelength(2)
	!else
	!	CFD_region = domain(2)/2.d0	!Top of processor domain
	!endif

               ! p_units_lb(:) = floor( (/ iblock-1, jblock -1, kblock - 1 /) &
               !         * domain(:) / initialunitsize(:) ) 
               ! p_units_ub(:) = ceiling(  (/ iblock,   jblock,    kblock /) &
               !         * domain(:) / initialunitsize(:) )
               ! nfcc_max = floor( globaldomain(:) / initialunitsize(:))

               ! do j = 1, nd
               !         p_units_ub(j) = min( p_units_ub(j), nfcc_max(j))
               ! enddo

               ! if (jblock == npy) then 
               !        p_units_ub(2) = p_units_ub(2) - 1
               ! endif

        !else
        !        p_units_lb(1) = (iblock-1)*floor(initialnunits(1)/real((npx),kind(0.d0)))
        !        p_units_ub(1) =  iblock *ceiling(initialnunits(1)/real((npx),kind(0.d0)))
        !        p_units_lb(2) = (jblock-1)*floor(initialnunits(2)/real((npy),kind(0.d0)))
        !        p_units_ub(2) =  jblock *ceiling(initialnunits(2)/real((npy),kind(0.d0)))
        !        p_units_lb(3) = (kblock-1)*floor(initialnunits(3)/real((npz),kind(0.d0)))
        !        p_units_ub(3) =  kblock *ceiling(initialnunits(3)/real((npz),kind(0.d0)))
       ! endif
	!print*, p_units_lb, p_units_ub

	!Molecules per unit FCC structure (3D)
	n  = 0  	!Reset n
	nl = 0		!Reset nl
	do nz=p_units_lb(3),p_units_ub(3)	!Loop over z column
	c(3) = (nz - 0.75d0)*initialunitsize(3) - halfdomain(3) 
	do ny=p_units_lb(2),p_units_ub(2)	!Loop over y column
	c(2) = (ny - 0.75d0)*initialunitsize(2) - halfdomain(2) 
	do nx=p_units_lb(1),p_units_ub(1)	!Loop over all x elements of y column
	c(1) = (nx - 0.75d0)*initialunitsize(1) - halfdomain(1)
		do j=1,4	!4 Molecules per cell
			do ixyz=1,nd	!For all dimensions
				if(j.eq.ixyz .or. j.eq.4) then
					rc(ixyz)=c(ixyz)
				else
					rc(ixyz)=c(ixyz)+0.5d0*initialunitsize(ixyz)
				endif
			enddo

			n = n + 1	!Move to next molecule

			!Remove molecules from top of domain if constraint applied
			if (rc(2)-domain(2)*(jblock-1) .gt.  CFD_region) cycle 

			!Check if molecule is in domain of processor
			if(rc(1).lt.-halfdomain(1)+domain(1)*(iblock-1)) cycle
			if(rc(1).ge. halfdomain(1)+domain(1)*(iblock-1)) cycle
			if(rc(2).lt.-halfdomain(2)+domain(2)*(jblock-1)) cycle
			if(rc(2).ge. halfdomain(2)+domain(2)*(jblock-1)) cycle
			if(rc(3).lt.-halfdomain(3)+domain(3)*(kblock-1)) cycle
			if(rc(3).ge. halfdomain(3)+domain(3)*(kblock-1)) cycle

			!If molecules is in the domain then add to total
			nl = nl + 1 !Local molecule count

			!Correct to local coordinates
			r(nl,1) = rc(1)-domain(1)*(iblock-1)
			r(nl,2) = rc(2)-domain(2)*(jblock-1)
			r(nl,3) = rc(3)-domain(3)*(kblock-1)
		enddo
	enddo
	enddo
	enddo

	np = nl			!Correct local number of particles on processor
	rinitial = r !Record initial position of all molecules

	!Establish global number of particles on current process
	globalnp = np
	call globalSumInt(globalnp)

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

	if (coupler_is_active .and. myid .eq. 0) then
		print*, '*********************************************************************'
		print*, '*WARNING - TOP LAYER OF DOMAIN REMOVED IN LINE WITH CONSTRAINT FORCE*'
		print*, 'Removed from', CFD_region, 'to Domain top', globaldomain(2)/2.d0
		print*, 'Number of molecules reduced from',  & 
			 4*initialnunits(1)*initialnunits(2)*initialnunits(3), 'to', np
		print*, '*********************************************************************'

                !print*, 'microstate ', minval(r(:,1)), maxval(r(:,1)),minval(r(:,2)), maxval(r(:,2)),minval(r(:,3)), maxval(r(:,3))
	endif


end subroutine setup_initialise_parallel_position

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!Initialise Velocities
! Set up the intial velocities of the particles using velocity magnitude calcualted
! from intitial temperature and giving random vectorial components

subroutine setup_initialise_velocities
	use module_initialise_microstate
	implicit none

	integer                            :: n,i 
	double precision, dimension (nd)   :: netv   !Overall momentum of system

	!Use definition of temperature and re-arrange to define an average velocity
	initialvel = sqrt(nd * (1.d0 - 1.d0/globalnp)*inputtemperature)

	v = 0.d0	!Set velocity initially to zero
	i = 0		!Zero number of molecules with velocity assigned
	netv=0.d0	!Set net velocity of system to zero initially

	do n=1,np			      		!Step through each molecule
		if (fix(n,1) .eq. 1) then     		!For x component as fix occurs in all 3 dimensions
			call random_number(rand)	!Generate a random number for each dimension/particle
			angle  = 2.d0*pi*rand          	!Random angle between 0 and 2pi
			v(n,2) = initialvel*sin(angle)	!Y component of velocity magnitude for random angle
			v13    = initialvel*cos(angle)	!Magnitude of x and z component
			call random_number(rand)    	!Generate a new random number
			angle  = 2.d0*pi*rand          	!Random angle between 0 and 2pi		
			v(n,1) = v13*sin(angle)        	!X component of velocity magnitude for random angle
			v(n,3) = v13*cos(angle)        	!Z component of velocity magnitude for random angle
			i = i + 1			!Count number of molecules with velocity assigned
		else
			v(n,:) = 0.d0			!Don't assign velocity if molecule is fixed
		endif
		netv(:)= netv(:) + v(n,:)      !Sum up overall momentum of system due to random movement
	enddo

	call globalSumVect(netv, nd)	       !Sum net velocity on all processors
	call globalSumInt(i)	       	       !Sum number of molecules assigned velocity on all processors

	if(i .ne. 0) netv(:) = netv(:)/i       !Divide overall momentum by number of particles

	do n=1,np
		!reducing all non-fixed particles by same amount
		if (fix(n,1) .eq. 1) v(n,:)= v(n,:)-netv(:)
			       
	enddo

end subroutine setup_initialise_velocities


!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!Initialise Velocities
! Set up the intial velocities of the particles using velocity magnitude calcualted
! from intitial temperature and giving random vectorial components

subroutine setup_initialise_velocities_test
	use module_initialise_microstate
	implicit none

	integer                            :: i, n
	double precision, dimension (nd)   :: netv   !Overall momentum of system

	!Use definition of temperature and re-arrange to define an average velocity
	!initialvel = sqrt(nd * (1.d0 - 1.d0/np)*inputtemperature)
	v = 0.d0	!Set velocity initially to zero
	!i = 0		!Zero number of molecules with velocity assigned
	!netv=0		!Set net velocity of system to zero initially
	!zeta=0.d0	!Set Nose Hoover thermostat scaling property to zero
	
	!v(1,1) = 2.d0

	do n=1,np			      		!Step through each molecule
		!r(1,:) = halfdomain(:)
		!v(n,1) = 1.0d0
		v(n,2) = 0.1d0
		!v(n,3) = 0.0d0

		!r(1,:) = -halfdomain(:)
		!v(n,1) = -0.0d0 
		!v(n,2) = -0.0d0
		!v(n,3) = -0.0d0
	
	enddo

end subroutine setup_initialise_velocities_test

