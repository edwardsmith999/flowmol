!------------------------------------------------------------------------------
!                           Initialise Microstate
!        Set up position and velocity of each molecule (the microstate)
!
!------------------------------------------------------------------------------

module module_initialise_microstate

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use calculated_properties_MD

	double precision			:: angle, rand  !Define variables
	double precision			:: v13          !Mag of v1 and v3 vectors
	
end module module_initialise_microstate
!------------------------------------------------------------------------------

subroutine setup_initialise_microstate
use interfaces
use module_initialise_microstate
implicit none

	integer		::	n

	select case(initial_config_flag)
	case(0)
		call setup_initialise_lattice        !Setup FCC lattice
	case(1)
		select case (config_special_case)
		case('fene_melt')
			call setup_initialise_lattice    !Numbering for FENE bonds
			call setup_initialise_polyinfo   !Chain IDs, etc
		case default
			call error_abort('Unidentified configuration special case')
		end select
	case(2)
		call error_abort('Initial configuration file input not yet developed')
	case default
		call error_abort('Unidentified initial configuration flag')	
	end select

!	select case(potential_flag)
!	case(0)	
!		call setup_initialise_parallel_position       !Setup initial pos
!	case(1) 
!		call setup_initialise_parallel_position_FENE  !Numbering for FENE bonds
!		call setup_initialise_polyinfo                !Chain IDs, etc
!	case default
!		call error_abort('Potential flag not recognised!')
!	end select

	do n=1,np    !Initialise global true positions
		rtrue(1,n) = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
		rtrue(2,n) = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
		rtrue(3,n) = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
	end do
	!rinitial = rtrue                                 !Store initial true pos

	call setup_tag                                    !Setup locn of fixed mols
	rtether = r                                       !Init tether pos
	do n = 1,np
		call read_tag(n)                              !Read tag, assign props
	enddo
	call setup_initialise_velocities				  !Setup initial velocities

end subroutine setup_initialise_microstate

!==================================================================================
!----------------------------------------------------------------------------------
!Initialise Positions
!Set up the intial position of the particles in an FCC lattice

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
						r(ixyz,n)=c(ixyz)
					else
						r(ixyz,n)=c(ixyz)+0.5d0*initialunitsize(ixyz)
					endif
				enddo
				n = n + 1   !Move to next molecule
				enddo
			enddo
		enddo
	enddo

	!rinitial = r !Record initial position of all molecules

end subroutine setup_initialise_position

subroutine setup_initialise_position_FENE
    use interfaces
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

	modcheck = 0 + mod(np,nmonomers) + mod(4*initialnunits(1),nmonomers)
	if (modcheck.ne.0) call error_abort('Number of molecules must be exactly divisible by &
	& the polymer chain length. Please change the chain length in the input file. &
	& A chain length of 4 should (hopefully) always work.')
	
	!Molecules per unit FCC structure (3D)
	n=1  !Reset n for start of loop
	do nz=1,initialnunits(3) !Loop over z column
	c(3) = (nz - 0.75d0)*initialunitsize(3) - halfdomain(3) 
		do ny=1,initialnunits(2) !Loop over y column
		c(2) = (ny - 0.75d0)*initialunitsize(2) - halfdomain(2) 
			do nx=1,initialnunits(1) !Loop over all x elements of y column
			c(1) = (nx - 0.75d0)*initialunitsize(1) - halfdomain(1)
				do j=1,4 !4 Molecules per cell
					r(:,n) = c(:)

					if (j.eq.2) then
						r(1,n) = c(1) + 0.5d0*initialunitsize(1)
						r(3,n) = c(3) + 0.5d0*initialunitsize(3)
					else if (j.eq.3) then
						r(2,n) = c(2) + 0.5d0*initialunitsize(2)
						r(3,n) = c(3) + 0.5d0*initialunitsize(3)
					else if (j.eq.4) then
						r(1,n) = c(1) + 0.5d0*initialunitsize(1)
						r(2,n) = c(2) + 0.5d0*initialunitsize(2)
					end if

					n = n + 1  !Move to next molecule

				enddo
			enddo
		enddo
	enddo

	!rinitial = r !Record initial position of all molecules

end subroutine setup_initialise_position_FENE

!----------------------------------------------------------------------------------
!Initialise Positions
!Set up the intial position of the particles

subroutine setup_initialise_parallel_position
	use module_initialise_microstate
	use messenger
#if USE_COUPLER
	use md_coupler_socket, only: socket_get_domain_top
#endif

	implicit none

	integer 						:: j, ixyz, n, nl, nx, ny, nz
	integer,dimension(nd) 			:: p_units_lb, p_units_ub, nfcc_max
	double precision 				:: domain_top, removed_height
	double precision, dimension (nd):: rc, c !Temporary variable

    p_units_lb(1) = (iblock-1)*floor(initialnunits(1)/real((npx),kind(0.d0)))
    p_units_ub(1) =  iblock *ceiling(initialnunits(1)/real((npx),kind(0.d0)))
    p_units_lb(2) = (jblock-1)*floor(initialnunits(2)/real((npy),kind(0.d0)))
    p_units_ub(2) =  jblock *ceiling(initialnunits(2)/real((npy),kind(0.d0)))
    p_units_lb(3) = (kblock-1)*floor(initialnunits(3)/real((npz),kind(0.d0)))
    p_units_ub(3) =  kblock *ceiling(initialnunits(3)/real((npz),kind(0.d0)))

	!Set CFD region to top of domain initially
	domain_top = domain(2)/2.d0

#if USE_COUPLER
	if (jblock .eq. npy) then
		domain_top = socket_get_domain_top()
	end if
#endif

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
			if (rc(2)-domain(2)*(jblock-1) .gt. domain_top) cycle 

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
			r(1,nl) = rc(1)-domain(1)*(iblock-1)
			r(2,nl) = rc(2)-domain(2)*(jblock-1)
			r(3,nl) = rc(3)-domain(3)*(kblock-1)
		enddo
	enddo
	enddo
	enddo

	np = nl			 !Correct local number of particles on processor
	!rinitial = rtrue !Record initial position of all molecules

	!Establish global number of particles on current process
	globalnp = np
	call globalSumInt(globalnp)

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

#if USE_COUPLER
	if (myid .eq. 0) then
		print*, '*********************************************************************'
		print*, '*WARNING - TOP LAYER OF DOMAIN REMOVED IN LINE WITH CONSTRAINT FORCE*'
		print*, 'Removed from', domain_top, 'to Domain top', globaldomain(2)/2.d0
		print*, 'Number of molecules reduced from this',  & 
			 4*initialnunits(1)*initialnunits(2)*initialnunits(3), 'to', globalnp
		print*, '*********************************************************************'

                !print*, 'microstate ', minval(r(1,:)), maxval(r(1,:)),minval(r(2,:)), maxval(r(2,:)),minval(r(3,:)), maxval(r(3,:))
	endif
#endif


end subroutine setup_initialise_parallel_position

!--------------------------------------------------------------------------------
!FENE equivalent
subroutine setup_initialise_lattice
	use module_initialise_microstate
	use messenger
#if USE_COUPLER
	use coupler
	use md_coupler_socket, only: socket_get_domain_top
#endif
	implicit none

	integer	:: j, n, nl, nx, ny, nz
	integer, dimension(nd) :: p_units_lb, p_units_ub 
	double precision :: domain_top
	double precision, dimension (nd):: rc, c !Temporary variable

	p_units_lb(1) = (iblock-1)*floor(initialnunits(1)/real((npx),kind(0.d0)))
	p_units_ub(1) =  iblock *ceiling(initialnunits(1)/real((npx),kind(0.d0)))
	p_units_lb(2) = (jblock-1)*floor(initialnunits(2)/real((npy),kind(0.d0)))
	p_units_ub(2) =  jblock *ceiling(initialnunits(2)/real((npy),kind(0.d0)))
	p_units_lb(3) = (kblock-1)*floor(initialnunits(3)/real((npz),kind(0.d0)))
	p_units_ub(3) =  kblock *ceiling(initialnunits(3)/real((npz),kind(0.d0)))

	!Set top of domain initially
	domain_top = domain(2)/2.d0

#if USE_COUPLER

	if (jblock .eq. npy) then
		domain_top = socket_get_domain_top()
	endif

#endif

	!Molecules per unit FCC structure (3D)
	n  = 0  	!Initialise global np counter n
	nl = 0		!Initialise local np counter nl

	!Inner loop in y (useful for setting connectivity)
	do nz=p_units_lb(3),p_units_ub(3)
	c(3) = (nz - 0.75d0)*initialunitsize(3) !- halfdomain(3) 
	do nx=p_units_lb(1),p_units_ub(1)
	c(1) = (nx - 0.75d0)*initialunitsize(1) !- halfdomain(1)
	do ny=p_units_lb(2),p_units_ub(2)
	c(2) = (ny - 0.75d0)*initialunitsize(2) !- halfdomain(2) 

		do j=1,4	!4 Molecules per cell

			rc(:) = c(:)
			select case(j)
			case(2)
				rc(1) = c(1) + 0.5d0*initialunitsize(1)
				rc(3) = c(3) + 0.5d0*initialunitsize(3)
			case(3)
				rc(2) = c(2) + 0.5d0*initialunitsize(1)
				rc(3) = c(3) + 0.5d0*initialunitsize(3)
			case(4)
				rc(1) = c(1) + 0.5d0*initialunitsize(1)
				rc(2) = c(2) + 0.5d0*initialunitsize(3)
			case default
			end select

			n = n + 1	!Move to next particle
			
			!Remove molecules from top of domain if constraint applied
			if (jblock .eq. npy) then
				if (rc(2)-domain(2)*(jblock-1)-halfdomain(2) .gt.  domain_top) cycle 
			endif

			!Check if molecule is in domain of processor
			if(rc(1).lt. domain(1)*(iblock-1)) cycle
			if(rc(1).ge. domain(1)*(iblock  )) cycle
			if(rc(2).lt. domain(2)*(jblock-1)) cycle
			if(rc(2).ge. domain(2)*(jblock  )) cycle
			if(rc(3).lt. domain(3)*(kblock-1)) cycle
			if(rc(3).ge. domain(3)*(kblock  )) cycle

			!If molecules is in the domain then add to total
			nl = nl + 1 !Local molecule count

			!Correct to local coordinates
			r(1,nl) = rc(1)-domain(1)*(iblock-1)-halfdomain(1)
			r(2,nl) = rc(2)-domain(2)*(jblock-1)-halfdomain(2)
			r(3,nl) = rc(3)-domain(3)*(kblock-1)-halfdomain(3)

		enddo

	enddo
	enddo
	enddo

	!Correct local number of particles on processor
	np = nl

	!Establish global number of particles on current process
	globalnp = np
	call globalSumInt(globalnp)

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

#if USE_COUPLER

	if (jblock .eq. npy .and. iblock .eq. 1 .and. kblock .eq. 1) then
		print*, '*********************************************************************'
		print*, '*WARNING - TOP LAYER OF DOMAIN REMOVED IN LINE WITH CONSTRAINT FORCE*'
		print*, 'Removed from', domain_top, 'to Domain top', globaldomain(2)/2.d0
		print*, 'Number of molecules reduced from',  & 
		         4*initialnunits(1)*initialnunits(2)*initialnunits(3), 'to', globalnp
		print*, '*********************************************************************'
		!print*, 'microstate ', minval(r(1,:)), maxval(r(1,:)),minval(r(2,:)),  &
		!         maxval(r(2,:)),minval(r(3,:)), maxval(r(3,:))
	endif

#endif

end subroutine setup_initialise_lattice
!-------------------------------------------------------------------------------
!Assign chainIDs, subchainIDs, global molecule numbers, etc...

subroutine setup_initialise_polyinfo
	use interfaces
	use polymer_info_MD
	use messenger
	use physical_constants_MD, only: np
	implicit none

	integer :: i,n
	integer :: chainID
	integer :: subchainID
	integer :: modcheck
	integer :: solvent_selector
	integer, dimension(nproc) :: proc_chains, proc_nps
	
	proc_chains(:)         = 0
	proc_nps(:)            = 0

	intbits = bit_size(monomer(1)%bin_bflag(1))
	
	modcheck = 0 + mod(np,nmonomers) + mod(4*initialnunits(1),nmonomers)
	if (modcheck.ne.0) call error_abort('Number of molecules must be exactly divisible by &
	& the polymer chain length. Please change the chain length in the input file. &
	& A chain length of 4 should (hopefully) always work.')

	do n=1,np+extralloc
		monomer(n)%chainID      = 0
		monomer(n)%subchainID   = 0
		monomer(n)%funcy        = 0
		monomer(n)%glob_no      = 0
		monomer(n)%bin_bflag(:) = 0
	end do	

	chainID    = 1
	subchainID = 0
	do n=1,np

		select case (solvent_flag)
		case (0)
			solvent_selector = 0
		case (1,2)
			solvent_selector = mod((n-1)/nmonomers,solvent_ratio)
		case default
		end select
	
		select case (solvent_selector)
		case(0) !POLYMER
			
			subchainID = subchainID + 1
			if (subchainID .gt. nmonomers) then
				subchainID = 1
				chainID    = chainID + 1
			end if
	
			monomer(n)%chainID    = chainID
			monomer(n)%subchainID = subchainID
			monomer(n)%glob_no    = n	
			
			if (subchainID.eq.1) then
				call connect_to_monomer(subchainID+1,n)
			else if (subchainID.eq.nmonomers) then
				call connect_to_monomer(subchainID-1,n)				
			else
				call connect_to_monomer(subchainID+1,n)
				call connect_to_monomer(subchainID-1,n)
			end if
			
		case(1:) !SOLVENT

			!SET MOLECULES TO BE ATHERMAL SOLVENT MOLECULES
			monomer(n)%chainID     = 0
			monomer(n)%subchainID  = 1
			monomer(n)%glob_no     = n
			monomer(n)%funcy       = 0
			bond(:,n)              = 0
		
		case default
		end select
			
	!	print*, '-n,c,sc,f,g------------------------------------------------'
	!	print*, n, monomer(n)%chainID,monomer(n)%subchainID,monomer(n)%funcy,monomer(n)%glob_no
	!	print*, 'bin_bflag -', monomer(n)%bin_bflag
	!	print*, '==========================================================='

	end do

	proc_chains(irank) = chainID
	proc_nps(irank)    = np
	call globalSumIntVect(proc_chains,nproc)
	call globalSumIntVect(proc_nps,nproc)
	
	do n=1,np
		if (monomer(n)%chainID.ne.0) then
			monomer(n)%chainID     = monomer(n)%chainID + sum(proc_chains(1:irank)) - proc_chains(irank)
		end if
		monomer(n)%glob_no     = monomer(n)%glob_no + sum(proc_nps(1:irank))    - proc_nps(irank)
	end do

	nchains = sum(proc_chains)

end subroutine setup_initialise_polyinfo

!=============================================================================
!Initialise branched polymer simulation
subroutine setup_initialise_polyinfo_singlebranched
	use interfaces
	use polymer_info_MD
	use messenger
	use physical_constants_MD, only: np
	use arrays_MD, only:r
	implicit none

	integer :: i,n
	integer :: chainID
	integer :: subchainID
	integer :: modcheck
	integer :: solvent_selector
	integer :: scIDbranch, glob_n_branch, nbranch, branchmonomers
	integer, dimension(nproc) :: proc_chains, proc_nps
	double precision, dimension(3) :: rij
	double precision :: rij2
	
	proc_chains(:)         = 0
	proc_nps(:)            = 0
	branchmonomers         = nmonomers    !Branch length specified in input
	nmonomers              = 2*nmonomers  !Two branches per polymer

	modcheck = 0 + mod(np,nmonomers) + mod(4*initialnunits(1),nmonomers)
	if (modcheck.ne.0) call error_abort('Number of molecules must be exactly &
	                                    &divisible by the polymer chain      &
	                                    &length. Please change the chain     &
	                                    &length in the input file. A chain   &
	                                    &length of 4 should (hopefully)      &
	                                    &work.')

	!Set all to solvents
	do n=1,np+extralloc
		monomer(n)%chainID      = 0
		monomer(n)%subchainID   = 1
		monomer(n)%funcy        = 0
		monomer(n)%glob_no      = n
		monomer(n)%bin_bflag(:) = 0
	end do	

	scIDbranch = branchmonomers/2            ! Choose subchainID from which
	                                         ! to create a branch
	chainID    = 1                           ! Initialise
	subchainID = 0                           ! Initialise
	!First branch
	do n=1,branchmonomers                    ! Loop over first branch monomers

		chainID = 1                          ! Set all polymers to same chain
		subchainID = subchainID + 1          
		monomer(n)%chainID = chainID
		monomer(n)%subchainID = subchainID
		if (subchainID.eq.1) then            ! The usual chain connections...
			call connect_to_monomer(subchainID+1,n)
		else if (subchainID.eq.branchmonomers) then
			call connect_to_monomer(subchainID-1,n)				
		else
			call connect_to_monomer(subchainID+1,n)
			call connect_to_monomer(subchainID-1,n)
		end if

		if (subchainID .eq. scIDbranch) then ! Store global mol number of
			glob_n_branch = n                ! branching monomer
		end if

	end do

	!Locate new monomer close to branching monomer
	do n=branchmonomers+1,np
		rij(:) = r(glob_n_branch,:) - r(:,n) ! Check separation
		rij2 = dot_product(rij,rij)
		if (rij2.lt.R_0**2.d0) then          ! If within spring max elongation
			nbranch = n                      ! Store global mol no
			exit                             ! Nothing else needed
		end if
	end do

	!Second branch
	do n=nbranch,nbranch+branchmonomers-1    ! Start from monomer close to
		chainID = 1                          ! branching monomer. Same cID.
		subchainID = subchainID + 1
		monomer(n)%chainID = chainID
		monomer(n)%subchainID = subchainID
		if (subchainID.eq.branchmonomers+1) then
			call connect_to_monomer(scIDbranch,n)   ! Connect monomer to
		                                            ! new branch 
			call connect_to_monomer(subchainID,glob_n_branch) ! And vice-versa
			call connect_to_monomer(subchainID+1,n) ! Continue chain
		else if (subchainID .gt. branchmonomers+1 .and. &
		         subchainID .lt. nmonomers) then
			call connect_to_monomer(subchainID+1,n)
			call connect_to_monomer(subchainID-1,n)
		else	
			call connect_to_monomer(subchainID-1,n) ! End chain
		end if
	end do

	proc_chains(irank) = chainID
	proc_nps(irank)    = np
	call globalSumIntVect(proc_chains,nproc)
	call globalSumIntVect(proc_nps,nproc)
	
	do n=1,np
		if (monomer(n)%chainID.ne.0) then
			monomer(n)%chainID = monomer(n)%chainID + &
			                     sum(proc_chains(1:irank)) - proc_chains(irank)
		end if
		monomer(n)%glob_no     = monomer(n)%glob_no + &
		                         sum(proc_nps(1:irank))    - proc_nps(irank)
	end do

	nchains = sum(proc_chains)

end subroutine setup_initialise_polyinfo_singlebranched

!-----------------------------------------------------------------------------
!Connect monomer to subchainID bscID
subroutine connect_to_monomer(bscID,n)
use polymer_info_MD
implicit none

	integer, intent(in) :: bscID,n
	integer :: group, expo
	integer :: otherglob
	
	intbits = bit_size(monomer(1)%bin_bflag(1))
	group = ceiling(real(bscID)/real(intbits))
	expo  = mod(bscID,intbits) - 1
	if(expo.eq.-1) expo = intbits - 1
	monomer(n)%funcy            = monomer(n)%funcy + 1
	monomer(n)%bin_bflag(group) = monomer(n)%bin_bflag(group) + 2**(expo)
	
end subroutine connect_to_monomer
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

	do n=1,np			      				!Step through each molecule
		if (fix(1,n) .eq. 1) then     		!For x component as fix occurs in all 3 dimensions
			call random_number(rand)		!Generate a random number for each dimension/particle
			angle  = 2.d0*pi*rand          	!Random angle between 0 and 2pi
			v(2,n) = initialvel*sin(angle)	!Y component of velocity magnitude for random angle
			v13    = initialvel*cos(angle)	!Magnitude of x and z component
			call random_number(rand)    	!Generate a new random number
			angle  = 2.d0*pi*rand          	!Random angle between 0 and 2pi		
			v(1,n) = v13*sin(angle)       	!X component of velocity magnitude for random angle
			v(3,n) = v13*cos(angle)        	!Z component of velocity magnitude for random angle
			i = i + 1						!Count number of molecules with velocity assigned
		else
			v(:,n) = 0.d0					!Don't assign velocity if molecule is fixed
		endif
		netv(:)= netv(:) + v(:,n)      		!Sum up overall momentum of system due to random movement
	enddo

	call globalSumVect(netv, nd)			!Sum net velocity on all processors
	call globalSumInt(i)					!Sum number of molecules assigned velocity on all processors

	if(i .ne. 0) netv(:) = netv(:)/i		!Divide overall momentum by number of particles

	do n=1,np
		!reducing all non-fixed particles by same amount
		if (fix(1,n) .eq. 1) v(:,n)= v(:,n)-netv(:) 
			       
	enddo

end subroutine setup_initialise_velocities

!Initialise Velocities
! Set up the intial velocities of the particles using Taylor Green velocity
! contours

subroutine setup_initialise_velocities_TG
	use module_initialise_microstate
	implicit none

	integer                            :: n,i 
	double precision				   :: x,y,z,Lx,Ly,Lz
	double precision, dimension (nd)   :: netv   !Overall momentum of system

	!Use definition of temperature and re-arrange to define an average velocity
	initialvel = sqrt(nd * (1.d0 - 1.d0/globalnp)*inputtemperature)

	v = 0.d0	!Set velocity initially to zero
	i = 0		!Zero number of molecules with velocity assigned
	netv=0.d0	!Set net velocity of system to zero initially

	do n=1,np			      				!Step through each molecule
		x  = r(1,n);    y  = r(2,n);    z  = r(3,n);
		Lx = halfdomain(1); Ly = halfdomain(2); Lz = halfdomain(3);	!Domain should be cubic...
		v(1,n) =  initialvel*sin(pi*x/Lx)*cos(pi*y/Ly)*cos(pi*z/Lz)
		v(2,n) = -initialvel*cos(pi*x/Lx)*sin(pi*y/Ly)*cos(pi*z/Lz)
		v(3,n) =  initialvel*cos(pi*x/Lx)*cos(pi*y/Ly)*sin(pi*z/Lz)
		netv(:)= netv(:) + v(:,n)      		!Sum up overall momentum of system due to random movement
	enddo

	call globalSumVect(netv, nd)			!Sum net velocity on all processors
	netv(:) = netv(:)/np		!Divide overall momentum by number of particles

	do n=1,np
		!reducing all particles by same amount
		v(:,n)= v(:,n) - netv(:) 
			       
	enddo

end subroutine setup_initialise_velocities_TG


subroutine setup_initialise_velocities_TG_parallel
	use module_initialise_microstate
	implicit none

	integer                            :: n,i 
	double precision				   :: x,y,z,Lx,Ly,Lz
	double precision, dimension (nd)   :: netv   !Overall momentum of system

	!Use definition of temperature and re-arrange to define an average velocity
	initialvel = sqrt(nd * (1.d0 - 1.d0/globalnp)*inputtemperature)

	v = 0.d0	!Set velocity initially to zero
	i = 0		!Zero number of molecules with velocity assigned
	netv=0.d0	!Set net velocity of system to zero initially

	do n=1,np			      				!Step through each molecule
		x  = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
	    y  = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
		z  = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
		Lx = 0.5d0*globaldomain(1); Ly = 0.5d0*globaldomain(2); Lz = 0.5d0*globaldomain(3);	!Domain should be cubic...
		v(1,n) =  initialvel*sin(pi*x/Lx)*cos(pi*y/Ly)*cos(pi*z/Lz)
		v(2,n) = -initialvel*cos(pi*x/Lx)*sin(pi*y/Ly)*cos(pi*z/Lz)
		v(3,n) =  initialvel*cos(pi*x/Lx)*cos(pi*y/Ly)*sin(pi*z/Lz)
		netv(:)= netv(:) + v(:,n)      		!Sum up overall momentum of system due to random movement
	enddo

	print*, 'before sum', irank, netv, nd

	call globalSumVect(netv, nd)			!Sum net velocity on all processors
	netv(:) = netv(:)/np		!Divide overall momentum by number of particles

	print*, 'after sum', irank, netv, nd

	do n=1,np
		!reducing all particles by same amount
		v(:,n)= v(:,n) - netv(:) 
			       
	enddo

end subroutine setup_initialise_velocities_TG_parallel

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

!	do n=1,np			      		!Step through each molecule
		!r(1,:) = halfdomain(:)
!		v(1,n) = 1.0d0
!		v(2,n) = 1.0d0
!		v(3,n) = 0.0d0

		!r(1,:) = -halfdomain(:)
		!v(1,n) = -0.0d0 
		!v(2,n) = -0.0d0
		!v(3,n) = -0.0d0
	
!	enddo
	
!	v(1,:) = 0.5d0
	
!	v(7,3) = -0.5d0
!	v(4,3) = 0.5d0

end subroutine setup_initialise_velocities_test

