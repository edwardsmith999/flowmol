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

	integer		::	n

	!call setup_initialise_position       	!Setup initial positions
	call setup_initialise_parallel_position !Setup initial positions in //el
	call setup_tag				!Setup location of fixed molecules
	do n = 1,np
		call read_tag(n)		!Read tag and assign properties
	enddo
	call setup_initialise_velocities     	!Setup initial velocities
	!call setup_initialise_velocities_test


end subroutine setup_initialise_microstate

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

!----------------------------------------------------------------------------------
!Initialise Positions
!Set up the intial position of the particles

subroutine setup_initialise_parallel_position
	use module_initialise_microstate
        use messenger
	use coupler_md_global_data, only : use_coupling
	implicit none

	integer :: j, ixyz, n, nl, nx, ny, nz
	integer,dimension(nd) :: p_units_lb, p_units_ub, nfcc_max
	double precision, dimension (nd) :: rc, c !Temporary variable

        if ( use_coupling ) then
                p_units_lb(:) = floor( (/ iblock-1, jblock -1, kblock - 1 /) &
                        * domain(:) / initialunitsize(:) ) 
                p_units_ub(:) = ceiling(  (/ iblock,   jblock,    kblock /) &
                        * domain(:) / initialunitsize(:) )
                nfcc_max = floor( globaldomain(:) / initialunitsize(:))
                do j = 1, nd
                        p_units_ub(j) = min( p_units_ub(j), nfcc_max(j))
                enddo
        else
                p_units_lb(1) = (iblock-1)*floor(initialnunits(1)/real((npx),kind(0.d0)))
                p_units_ub(1) =  iblock *ceiling(initialnunits(1)/real((npx),kind(0.d0)))
                p_units_lb(2) = (jblock-1)*floor(initialnunits(2)/real((npy),kind(0.d0)))
                p_units_ub(2) =  jblock *ceiling(initialnunits(2)/real((npy),kind(0.d0)))
                p_units_lb(3) = (kblock-1)*floor(initialnunits(3)/real((npz),kind(0.d0)))
                p_units_ub(3) =  kblock *ceiling(initialnunits(3)/real((npz),kind(0.d0)))
        endif
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

! Should one shift r initial because initialunitsize is not conmensurate with domain size ?


	!Establish global number of particles on current process
	globalnp = np
	call globalSumInt(globalnp)

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp


end subroutine setup_initialise_parallel_position

!----------------------------------------------------------------------------------
!Create an array of fixed molecules

subroutine setup_fix
	use module_initialise_microstate
	implicit none

	integer :: n, mol_layers
	double precision, dimension(3) :: fixdist

	!For initialunitsize "a"
	!		 [  o     o ]
	!a (1 cell size) [     o    ]  a/2 (distance between molcules)	
	!		 [  o     o
	!		  __________]  a/4 (distance from bottom of domain)
	!
	!So use (0.20+0.5d0*mol_layers)*initialunitsize(ixyz)

	
	mol_layers = 2

	!+0.2 to include gap between wall and 1st molecule
	fixdist(1) = 0.d0 !initialunitsize(1)
	fixdist(2) = 2.d0 !1.7029d0 !3.4058197d0  !(0.20+0.5d0*mol_layers)*initialunitsize(2) 
	fixdist(3) = 0.d0 !initialunitsize(3)

	!Set all molecules fix equal to one (unfixed)
	fix = 1

	!Fix x bottom
	if (iblock .eq. 1) then
		do n = 1 , np
			if(r(n,1).lt.-halfdomain(1)+fixdist(1)) fix(n,:) = 0
		enddo
	endif
		
	!Fix x top	
	if (iblock .eq. npx) then
		do n = 1 , np
			if(r(n,1).ge. halfdomain(1)-fixdist(1)) fix(n,:) = 0
		enddo
	endif

	!Fix y bottom
	if (jblock .eq. 1) then
		do n = 1 , np
			if(r(n,2).lt.-halfdomain(2)+fixdist(2)) fix(n,:) = 0
			if(r(n,2).lt.-halfdomain(2)+fixdist(2)) slidev(n,1) = 3.0d0
		enddo
	endif

	!Fix y top
	if (jblock .eq. npy) then
		do n = 1 , np
			if(r(n,2).ge. halfdomain(2)-fixdist(2)) fix(n,:) = 0 
		enddo
	endif

	!Fix z bottom
	if (kblock .eq. 1) then
		do n = 1 , np
			if(r(n,3).lt.-halfdomain(3)+fixdist(3)) fix(n,:) = 0
		enddo
	endif

	!Fix z top
	if (kblock .eq. npz) then
		do n = 1 , np
			if(r(n,3).ge. halfdomain(3)-fixdist(3)) fix(n,:) = 0
		enddo
	endif

end subroutine setup_fix

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
	initialvel = sqrt(nd * (1.d0 - 1.d0/np)*inputtemperature)
	v = 0.d0	!Set velocity initially to zero
	i = 0		!Zero number of molecules with velocity assigned
	netv=0		!Set net velocity of system to zero initially
	zeta=0.d0	!Set Nose Hoover thermostat scaling property to zero
	
	v(1,1) = 2.d0

	!do n=1,np			      		!Step through each molecule
		!r(1,:) = halfdomain(:)
		!v(n,1) = 1.0d0
		!v(n,2) = 0.0d0
		!v(n,3) = 0.0d0

		!r(1,:) = -halfdomain(:)
		!v(n,1) = -0.0d0 
		!v(n,2) = -0.0d0
		!v(n,3) = -0.0d0
	
	!enddo

end subroutine setup_initialise_velocities_test

