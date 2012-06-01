!-----------------------------------------------------------------------------
!
!                            Compute_forces                               
! Establish the forces applied to the particles as a result of interactions
! using the lennard jones potential in dimensionless units
! For a cell, molecule by molecule, the interactions with ALL of the 
! neighbouring cells' molecules are determined.
!
!-----------------------------------------------------------------------------

module module_compute_forces

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use calculated_properties_MD
	use linked_list

	double precision                :: rij2        	!rij dotted with itself (one dimensional)
	double precision                :: invrij2     	!inverse of rij2 
	double precision                :: accijmag    	!Non directional component of acceleration
	double precision,dimension(3)   :: ri, rj 		!Position of molecule i and j
	double precision,dimension(3)   :: rij    		!vector between particles i and j
	double precision,dimension(3)   :: fij    		!force between particles i and j

end module module_compute_forces

!========================================================================
! Force calculation using the force list method and potential type
! specified in the input file

subroutine simulation_compute_forces
	use interfaces
	use module_compute_forces
	use polymer_info_MD, only: solvent_flag
	implicit none
	
	a					= 0.d0	!Reset acceleration matrix before force calculations
	potenergymol		= 0.d0	!Reset potential energy per molecule before calculation
	potenergymol_LJ		= 0.d0	!Reset LJ energy per molecule before calculation
	potenergysum		= 0.d0  !Reset total potential energy sum before calculation
	potenergysum_LJ		= 0.d0  !Reset LJ potential energy sum before calculation
	virial				= 0.d0	!Reset virial sum before calculation
	virialmol			= 0.d0	!Reset virial sum per molecule before calculation

	!Choose between all pairs, cell list or neighbour list
	select case(force_list)
	case(0)
		!Forces calculated using all pairs
		select case(potential_flag)
		case(0)
			call simulation_compute_forces_LJ_AP
		case default
			call error_abort("Potential flag/force_list incompatible - only LJ available with all pairs")
		end select

	case(1)
		!Forces calculated using cell lists
		select case(potential_flag)
		case(0)					!If simple LJ fluid
			call simulation_compute_forces_LJ_cells
		case default								
			call error_abort("Potential flag/force_list incompatible - only LJ available with cell lists")
		end select

	case(2)
		!Forces calculated using neighbour lists with all interactions
		select case(potential_flag)
		case(0)					!If simple LJ fluid
			call simulation_compute_forces_LJ_neigbr
		case default								
			call error_abort("Potential flag/force_list incompatible - &
                        &only LJ available with all interactions neighbour lists")
		end select
	case(3)
		!Forces calculated using neighbour lists optimised using 
		!Newton's 3rd law to count only half of the interactions
		select case(potential_flag)
		case(0)					!If simple LJ fluid
			call simulation_compute_forces_LJ_neigbr_halfint
		case(1)					!If FENE polymer
			potenergymol_FENE	= 0.d0
			potenergysum_FENE	= 0.d0
			select case(solvent_flag)
			case(0:1)
				call simulation_compute_forces_LJ_neigbr_halfint	!Compute LJ bead interactions
				call simulation_compute_forces_FENE					!Add on FENE spring interactions
			case(2)
				!call simulation_compute_forces_Soddemann_AP
				call simulation_compute_forces_Soddemann_neigbr_halfint
				call simulation_compute_forces_FENE
			case default
				call error_abort('Solvent flag not recognised!')
			end select
		case default								
			call error_abort('Potential flag not recognised!')
		end select
	case default
	end select

end subroutine simulation_compute_forces

!========================================================================
!Compute forces using cells instead of neighbourlist

subroutine simulation_compute_forces_LJ_AP
	use module_compute_forces
	implicit none

	integer                         :: i, j,ixyz   !Define dummy index

	do i = 1,np						!Step through each particle in list 
		ri = r(i,:)         		!Retrieve ri
		do j = i+1,np				!Step through all j for each i
			rj = r(j,:)				!Retrieve rj

			!Calculate rij using nearest image convention, i.e. if more than half
			! a domain betwen molecules, must be closer over periodic boundaries  
			rij2 = 0.d0
			do ixyz=1,nd
				rij(ixyz) = r (i,ixyz) - r(j,ixyz)          !Evaluate distance between particle i and j
    				if (abs(rij(ixyz)) > halfdomain(ixyz)) then	
					rij(ixyz) = rij(ixyz) - sign(domain(ixyz),rij(ixyz)) 
				endif
				rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
			enddo

			if (rij2 < rcutoff2) then

				!Linear magnitude of acceleration for each molecule
				invrij2 = 1.d0/rij2                 !Invert value
				accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)

				!Sum of forces on particle i added for each j
				a(i,1)= a(i,1) + accijmag*rij(1)
				a(i,2)= a(i,2) + accijmag*rij(2)
				a(i,3)= a(i,3) + accijmag*rij(3)

				!Sum of forces on particle j added for each i
				a(j,1)= a(j,1) - accijmag*rij(1)
				a(j,2)= a(j,2) - accijmag*rij(2)
				a(j,3)= a(j,3) - accijmag*rij(3)

				!Only calculate properties when required for output
				if (mod(iter,tplot) .eq. 0) then
					!Record potential energy total to use for output later (potshift=-1 for WCA)
					potenergymol_LJ(i)=potenergymol_LJ(i) & 
						     +4.d0*(invrij2**6-invrij2**3)-potshift
					potenergymol_LJ(j)=potenergymol_LJ(j) & 
						     +4.d0*(invrij2**6-invrij2**3)-potshift
					!Virial expression used to obtain pressure
					virialmol(i) = virialmol(i) + accijmag*rij2
					virialmol(j) = virialmol(j) + accijmag*rij2
				endif
			endif
		enddo
	enddo

	!Total used with other potentials (e.g. FENE)
	potenergymol = potenergymol + potenergymol_LJ

end subroutine simulation_compute_forces_LJ_AP

subroutine simulation_compute_forces_Soddemann_AP
use interfaces
use module_compute_forces
use polymer_info_MD
implicit none

	integer                         :: n,i,j,ixyz
	integer                         :: p_i, p_j, ptot
	double precision, parameter     :: sod_a    = 3.1730728678
	double precision, parameter     :: sod_b    = -0.85622864544
	double precision, parameter     :: wca_cut  = 1.12246204830937
!	double precision, parameter     :: wca_cut2 = 1.25992104989486
	double precision, parameter     :: wca_cut2 = 1.25992104989487
	double precision                :: eps

	do i = 1,np
		ri = r(i,:)
		do j = i+1,np                !Step through all pairs

			rj     = r(j,:)
			rij    = ri - rj
			rij(:) = rij(:) - domain(:)*anint(rij(:)/domain(:))
			rij2   = dot_product(rij,rij)

			if (rij2 .lt. sod_cut2) then

				p_i = 0
				p_j = 0
				if (monomer(i)%chainID .ne. 0) p_i = 1
				if (monomer(j)%chainID .ne. 0) p_j = 1
				ptot = p_i + p_j
			
				!Linear magnitude of acceleration for each bead---------------
				invrij2  = 1.d0/rij2             !Invert value
				select case (ptot)
				case(0)
					eps = eps_ss                 !Solvent-solvent interaction
				case(1)
					eps = eps_ps                 !Polymer-solvent interaction
				case(2)
					eps = eps_pp                 !Polymer-polymer interaction
				                                 !(no FENE)
				case default
					call error_abort("Undetermined interaction in &
				                      &compute_forces_Soddemann")
				end select
			
				if (rij2 .lt. wca_cut2) then
					accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)
				else
					accijmag = (eps*sod_a*invrij2**2.d0)*sin(invrij2*sod_a &
				                                          + sod_b)
				end if          
				!-------------------------------------------------------------
	
				!Sum of forces on particle i added for each j
				a(i,1)= a(i,1) + accijmag*rij(1)
				a(i,2)= a(i,2) + accijmag*rij(2)
				a(i,3)= a(i,3) + accijmag*rij(3) 

				!Sum of forces on particle j added for each i
				a(j,1)= a(j,1) - accijmag*rij(1)
				a(j,2)= a(j,2) - accijmag*rij(2)
				a(j,3)= a(j,3) - accijmag*rij(3) 

				!Only calculate properties when required for output
				if (mod(iter,tplot) .eq. 0) then
					!Record potential energy total to use for output later
					if (rij2 .lt. wca_cut2) then
						potenergymol_LJ(i) = potenergymol_LJ(i)   &
					    + 4.d0*(invrij2**6.d0 - invrij2**3.d0 + 0.25d0) - eps 
						potenergymol_LJ(j) = potenergymol_LJ(j)   &
					    + 4.d0*(invrij2**6.d0 - invrij2**3.d0 + 0.25d0) - eps 
					else
						potenergymol_LJ(i) = potenergymol_LJ(i)   &
						+ 0.5d0*eps*(cos(sod_a*rij2 + sod_b) - 1.d0)
						potenergymol_LJ(j) = potenergymol_LJ(j)   &
						+ 0.5d0*eps*(cos(sod_a*rij2 + sod_b) - 1.d0)
					end if
					!Virial expression used to obtain pressure
					virialmol(i) = virialmol(i) + accijmag*rij2
					virialmol(j) = virialmol(j) + accijmag*rij2
				endif
			endif
		enddo
	enddo

	!Total used with other potentials (e.g. FENE)
	potenergymol = potenergymol + potenergymol_LJ

end subroutine simulation_compute_forces_Soddemann_AP

!========================================================================
!Compute forces using cells instead of neighbourlist

subroutine simulation_compute_forces_LJ_cells
	use module_compute_forces
	implicit none

	integer                         :: i,j,n, ixyz   	!Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp
	integer							:: molnoi, molnoj, memloc
	type(node), pointer 	        :: oldi, currenti, oldj, currentj


!TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
!	do n = 1,2
!	 	a = 0.d0
!		potenergymol		= 0.d0	!Reset potential energy per molecule before calculation
!		potenergymol_LJ		= 0.d0	!Reset LJ energy per molecule before calculation
!		potenergysum		= 0.d0  !Reset total potential energy sum before calculation
!		potenergysum_LJ		= 0.d0  !Reset LJ potential energy sum before calculation
!		virial				= 0.d0	!Reset virial sum before calculation
!		virialmol			= 0.d0	!Reset virial sum per molecule before calculation
!TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP

	do kcell=2, ncells(3)+1
	do jcell=2, ncells(2)+1
	do icell=2, ncells(1)+1

		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp					!Step through each particle in list 
			molnoi = oldi%molno 	 	!Number of molecule
			ri = r(molnoi,:)         	!Retrieve ri

			do kcellshift = -1,1
			do jcellshift = -1,1
			do icellshift = -1,1
				oldj => cell%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
				adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

				do j = 1,adjacentcellnp			!Step through all j for each i

					molnoj = oldj%molno			!Number of molecule
					rj = r(molnoj,:)			!Retrieve rj

					currentj => oldj
					oldj => currentj%next		!Use pointer in datatype to obtain next item in list

					if(molnoi==molnoj) cycle	!Check to prevent interaction with self

					rij2=0						!Set rij^2 to zero
					rij(:) = ri(:) - rj(:)		!Evaluate distance between particle i and j
					rij2 = dot_product(rij,rij)	!Square of vector calculated

					if (rij2 < rcutoff2) then

						!Linear magnitude of acceleration for each molecule
						invrij2 = 1.d0/rij2                 !Invert value
						accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)
	
						!Sum of forces on particle i added for each j
						a(molnoi,1)= a(molnoi,1) + accijmag*rij(1)
						a(molnoi,2)= a(molnoi,2) + accijmag*rij(2)
						a(molnoi,3)= a(molnoi,3) + accijmag*rij(3)

						!CV stress an force calculations
!TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP#
					!if (n .eq. 2) then
!TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP#
						if (vflux_outflag .eq. 4) then
							if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
								!if (molnoj .gt. np .or. molnoi .gt. np) then
								!	fij = 2.d0*accijmag*rij(:)
									!call Control_Volume_Forces(fij,ri,rj,molnoi,molnoj)
								!	call Control_Volume_stresses(fij,ri,rj,molnoi,molnoj)
								!else
									fij = accijmag*rij(:)
									!call Control_Volume_Forces(fij,ri,rj,molnoi,molnoj)
									call Control_Volume_stresses(fij,ri,rj,molnoi,molnoj)
								!endif
							endif
						endif
!TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP#
					!endif
!TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP#

						!Only calculate properties when required for output
						if (mod(iter,tplot) .eq. 0) then
							!Record potential energy total to use for output later (potshift=-1 for WCA)
							potenergymol_LJ(molnoi)=potenergymol_LJ(molnoi) & 
								     +4.d0*(invrij2**6-invrij2**3)-potshift
							!Virial expression used to obtain pressure
							virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
							if (pressure_outflag .eq. 1) call pressure_tensor_forces(molnoi,rij,accijmag)
							!if (pressure_outflag .eq. 2) call pressure_tensor_forces_VA(ri,rj,rij,accijmag)
							if (vflux_outflag.eq.1)	call pressure_tensor_forces_MOP(1,ri(:),rj(:),rij(:),accijmag)
							if (vflux_outflag.eq.2)	call pressure_tensor_forces_MOP(2,ri(:),rj(:),rij(:),accijmag)
							if (vflux_outflag.eq.3)	call pressure_tensor_forces_MOP(3,ri(:),rj(:),rij(:),accijmag)
						endif
					endif
				enddo
			enddo
			enddo
			enddo
			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo
	enddo
	enddo

	nullify(oldi)      		!Nullify as no longer required
	nullify(oldj)      		!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

	!Total used with other potentials (e.g. FENE)
	potenergymol = potenergymol + potenergymol_LJ

!TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP#
!		aold = a
!	enddo
!TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP

end subroutine simulation_compute_forces_LJ_cells

!-----------------------------------------------------------------

subroutine simulation_compute_forces_LJ_neigbr
	use module_compute_forces
	implicit none

	integer                         :: i, j, ixyz   !Define dummy index
	integer							:: molnoi, molnoj
	integer							:: noneighbrs
	type(neighbrnode), pointer		:: old, current

	do molnoi = 1, np

        noneighbrs = neighbour%noneighbrs(molnoi)	!Determine number of elements in neighbourlist
		old => neighbour%head(molnoi)%point			!Set old to head of neighbour list
		ri(:) = r(molnoi,:)							!Retrieve ri

		do j = 1,noneighbrs							!Step through all pairs of neighbours molnoi and j

			molnoj = old%molnoj						!Number of molecule j
			rj(:) = r(molnoj,:)						!Retrieve rj
			rij(:) = ri(:) - rj(:)   				!Evaluate distance between particle i and j
			rij2 = dot_product(rij,rij)				!Square of vector calculated

			if (rij2 < rcutoff2) then
				invrij2  = 1.d0/rij2                !Invert value
				accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4) ! (-dU/dr)*(1/|r|)
!				if (molnoi.eq.1) print*, accijmag
				!Sum of forces on particle i added for each j
				a(molnoi,1)= a(molnoi,1) + accijmag*rij(1)
				a(molnoi,2)= a(molnoi,2) + accijmag*rij(2)
				a(molnoi,3)= a(molnoi,3) + accijmag*rij(3)

				!CV stress an force calculations
				if (vflux_outflag .eq. 4) then
					if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
					!if (molnoj .gt. np .or. molnoi .gt. np) then
					!	fij = 2.d0*accijmag*rij(:)
						!call Control_Volume_Forces(fij,ri,rj,molnoi,molnoj)	
					!	call control_volume_stresses(fij,ri,rj,molnoi,molnoj)
					!else
						fij = accijmag*rij(:)
						!call Control_Volume_Forces(fij,ri,rj,molnoi,molnoj)
						call control_volume_stresses(fij,ri,rj,molnoi,molnoj)
					endif
				endif

				!Only calculate properties when required for output
				if (mod(iter,tplot) .eq. 0) then

					!Record potential energy total to use for output later (potshift=-1 for WCA)
					potenergymol_LJ(molnoi)=potenergymol_LJ(molnoi) & 
						     +4.d0*(invrij2**6-invrij2**3)-potshift
					!Virial expression used to obtain pressure
					virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2

					if (pressure_outflag .eq. 1) call pressure_tensor_forces(molnoi,rij,accijmag)
					!if (pressure_outflag .eq. 2) call pressure_tensor_forces_VA(ri,rj,rij,accijmag)
					if (vflux_outflag.eq.1)	call pressure_tensor_forces_MOP(1,ri(:),rj(:),rij(:),accijmag)
					if (vflux_outflag.eq.2)	call pressure_tensor_forces_MOP(2,ri(:),rj(:),rij(:),accijmag)
					if (vflux_outflag.eq.3)	call pressure_tensor_forces_MOP(3,ri(:),rj(:),rij(:),accijmag)
				endif
			endif
			current => old
			old => current%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo

	!Total used with other potentials (e.g. FENE)
	potenergymol = potenergymol + potenergymol_LJ

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine simulation_compute_forces_LJ_neigbr

!========================================================================
!Compute forces using neighbourlist with only half the interactions

subroutine simulation_compute_forces_LJ_neigbr_halfint
	use module_compute_forces
	implicit none

	integer                         :: n, j, ixyz   !Define dummy index
	integer							:: molnoi, molnoj
	integer							:: noneighbrs
	type(neighbrnode), pointer		:: old, current

	do molnoi = 1, np

	    noneighbrs = neighbour%noneighbrs(molnoi)	!Determine number of elements in neighbourlist
		old => neighbour%head(molnoi)%point			!Set old to head of neighbour list
		ri(:) = r(molnoi,:)							!Retrieve ri

		do j = 1,noneighbrs							!Step through all pairs of neighbours i and j

			molnoj = old%molnoj			!Number of molecule j
			rj(:) = r(molnoj,:)			!Retrieve rj
			rij(:)= ri(:) - rj(:)   	!Evaluate distance between particle i and j
			rij2  = dot_product(rij,rij)!Square of vector calculated
	
			if (rij2 .lt. rcutoff2) then

				!Linear magnitude of acceleration for each molecule
				invrij2  = 1.d0/rij2                 !Invert value
				accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)
	
				!Sum of forces on particle i added for each j
				a(molnoi,1)= a(molnoi,1) + accijmag*rij(1)
				a(molnoi,2)= a(molnoi,2) + accijmag*rij(2)
				a(molnoi,3)= a(molnoi,3) + accijmag*rij(3) 

				!Sum of forces on particle j added for each i
				a(molnoj,1)= a(molnoj,1) - accijmag*rij(1)
				a(molnoj,2)= a(molnoj,2) - accijmag*rij(2)
				a(molnoj,3)= a(molnoj,3) - accijmag*rij(3) 
				if (vflux_outflag.eq.4) then
					if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
						fij = 2.d0*accijmag*rij(:)
						!call Control_Volume_Forces(fij,ri,rj,molnoi,molnoj)
						call control_volume_stresses(fij,ri,rj,molnoi,molnoj)
					endif
				endif

				!Only calculate properties when required for output
				if (mod(iter,tplot) .eq. 0) then

					!Record potential energy total to use for output later (potshift=-1 for WCA)
					potenergymol_LJ(molnoi)=potenergymol_LJ(molnoi)+4.d0*(invrij2**6-invrij2**3)-potshift
					potenergymol_LJ(molnoj)=potenergymol_LJ(molnoj)+4.d0*(invrij2**6-invrij2**3)-potshift
			
					!Virial expression used to obtain pressure
					virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
					virialmol(molnoj) = virialmol(molnoj) + accijmag*rij2

					!Factor of two for half interaction virial calculation
					if (pressure_outflag .eq. 1) then
						call pressure_tensor_forces(molnoi,rij,accijmag)
						if (molnoj .le. np) call pressure_tensor_forces(molnoj,rij,accijmag)
					endif
					!if (pressure_outflag .eq. 2) call pressure_tensor_forces_VA(ri,rj,rij,accijmag)
					if (vflux_outflag.eq.1)	call pressure_tensor_forces_MOP(1,ri(:),rj(:),rij(:),accijmag)
					if (vflux_outflag.eq.2)	call pressure_tensor_forces_MOP(2,ri(:),rj(:),rij(:),accijmag)
					if (vflux_outflag.eq.3)	call pressure_tensor_forces_MOP(3,ri(:),rj(:),rij(:),accijmag)
				endif
			endif
			current => old
			old => current%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo

	!Total used with other potentials (e.g. FENE)
	potenergymol = potenergymol + potenergymol_LJ

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine simulation_compute_forces_LJ_neigbr_halfint

!========================================================================
!Compute polymer FENE potential forces using monomer bond lists

subroutine simulation_compute_forces_FENE
	use module_compute_forces
	use polymer_info_MD
	implicit none

	integer	:: molnoi,molnoj						!Current LJ bead
	integer :: b
	
	do molnoi=1,np

		ri(:) = r(molnoi,:)							!Retrieve ri(:)

		do b=1,monomer(molnoi)%funcy

			molnoj = bond(molnoi,b)
			if (molnoj.eq.0) cycle

			rj(:)  = r(molnoj,:)
			rij(:) = ri(:) - rj(:)
			rij2   = dot_product(rij,rij)

			if(rij2.ge.R_0**2)	call polymer_bond_error(molnoj)

			accijmag = -k_c/(1-(rij2/(R_0**2)))			!(-dU/dr)*(1/|r|)

			a(molnoi,1)= a(molnoi,1) + accijmag*rij(1)	!Add components of acceleration
			a(molnoi,2)= a(molnoi,2) + accijmag*rij(2)
			a(molnoi,3)= a(molnoi,3) + accijmag*rij(3)

			if (mod(iter,tplot) .eq. 0) then
				if (pressure_outflag .eq. 1) then
					call pressure_tensor_forces(molnoi,rij,accijmag)
				endif
				potenergymol_FENE(molnoi) = potenergymol_FENE(molnoi) - 0.5d0*k_c*R_0*R_0*dlog(1.d0-(rij2/(R_0**2)))
				potenergymol(molnoi)      = potenergymol(molnoi)      + potenergymol_FENE(molnoi)
				virialmol(molnoi)         = virialmol(molnoi)         + accijmag*rij2
			endif

		end do	

	end do

contains

subroutine polymer_bond_error(molnoX)
	use interfaces
	implicit none

	integer :: molnoX

	print*, 'irank: ', irank
	print '(a,i6,a,i4,a,i4,a,f8.5,a,f8.5,a,f12.5)', & 
			 	'Bond broken at iter ',iter,': atoms ',molnoi,' and ',molnoX,' are separated by ', &
				rij2**0.5,', which is greater than the allowed limit of ', R_0, &
				'. Stopping simulation, total time elapsed = ', iter*delta_t
	print '(a)', 'Atomic positions:'
	print '(a,i4,a,f10.5,a,f10.5,a,f10.5)', 'Atom ',molnoi,' is located at ',r(molnoi,1),' ',r(molnoi,2),' ',r(molnoi,3) 
	print '(a,i4,a,f10.5,a,f10.5,a,f10.5)', 'Atom ',molnoX,' is located at ',r(molnoX,1),' ',r(molnoX,2),' ',r(molnoX,3) 

	print '(a,i4,a)', 'Monomer information for atom ', molnoi,':'
	print '(a,i8)', 'ChainID: '   , monomer(molnoi)%chainID
	print '(a,i8)', 'SubchainID: ', monomer(molnoi)%subchainID
	print '(a,i8)', 'Funcy: '     , monomer(molnoi)%funcy
	print '(a,i8)', 'Glob_no: '   , monomer(molnoi)%glob_no
	print '(a,i8)', 'Bin_bflag: ' , monomer(molnoi)%bin_bflag

	print '(a,i4,a)', 'Monomer information for atom ', molnoX,':'
	print '(a,i8)', 'ChainID: '   , monomer(molnoX)%chainID
	print '(a,i8)', 'SubchainID: ', monomer(molnoX)%subchainID
	print '(a,i8)', 'Funcy: '     , monomer(molnoX)%funcy
	print '(a,i8)', 'Glob_no: '   , monomer(molnoX)%glob_no
	print '(a,i8)', 'Bin_bflag: ' , monomer(molnoX)%bin_bflag

	if (molnoX.gt.np) then
		call error_abort('Halo!')
	else
		call error_abort('')
	end if

end subroutine polymer_bond_error

end subroutine simulation_compute_forces_FENE
!==============================================================================
!Compute solvent-solvent, solvent-monomer, monomer-monomer (-FENE) forces
subroutine simulation_compute_forces_Soddemann_neigbr_halfint
use interfaces
use module_compute_forces
use polymer_info_MD
implicit none

	integer                         :: p_i, p_j, ptot
	integer							:: molnoi, molnoj, j
	integer							:: noneighbrs
	type(neighbrnode), pointer		:: old, current
	double precision, parameter     :: sod_a    = 3.1730728678
	double precision, parameter     :: sod_b    = -0.85622864544
	double precision, parameter     :: wca_cut  = 1.12246204830937
	double precision, parameter     :: wca_cut2 = 1.25992104989487
	double precision                :: eps

	do molnoi = 1, np

	    noneighbrs = neighbour%noneighbrs(molnoi)	!elements in neighbour list
		old => neighbour%head(molnoi)%point			!old>head of neighbour list
		ri(:) = r(molnoi,:)							!Retrieve ri

		do j = 1,noneighbrs                         !Step through all pairs
		                                            !of neighbours i and j
			molnoj = old%molnoj                     !Number of molecule j
			rj     = r(molnoj,:)                    !Position
			rij    = ri - rj                        !Difference
			rij(:) = rij(:) - domain(:)*anint(rij(:)/domain(:))!Min image
			rij2   = dot_product(rij,rij)           !Square

			if (rij2 .lt. sod_cut2) then            !If within potential range

				p_i = 0                             !Init as solvent
				p_j = 0                             !Init as solvent
				if (monomer(molnoi)%chainID .ne. 0) p_i = 1 !Flag polymer
				if (monomer(molnoj)%chainID .ne. 0) p_j = 1 !Flag polymer
				ptot = p_i + p_j                    !Find flag total
			
				!Linear magnitude of acceleration for each bead---------------
				select case (ptot)
				case(0)
					eps = eps_ss                 !Solvent-solvent interaction
				case(1)
					eps = eps_ps                 !Polymer-solvent interaction
				case(2)
					eps = eps_pp                 !Polymer-polymer interaction
				                                 !(no FENE)
				case default
					call error_abort("Undetermined interaction in &
				                      &compute_forces_Soddemann")
				end select
			
				invrij2  = 1.d0/rij2             !Useful value
				if (rij2 .lt. wca_cut2) then
					accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)
				else
					accijmag = (eps*sod_a*invrij2**2.d0)*sin(invrij2*sod_a &
				                                             + sod_b)
				end if          
				!-------------------------------------------------------------
	
				!Sum of forces on particle i added for each j
				a(molnoi,1)= a(molnoi,1) + accijmag*rij(1)
				a(molnoi,2)= a(molnoi,2) + accijmag*rij(2)
				a(molnoi,3)= a(molnoi,3) + accijmag*rij(3) 

				!Sum of forces on particle j added for each i
				a(molnoj,1)= a(molnoj,1) - accijmag*rij(1)
				a(molnoj,2)= a(molnoj,2) - accijmag*rij(2)
				a(molnoj,3)= a(molnoj,3) - accijmag*rij(3) 

				!Only calculate properties when required for output
				if (mod(iter,tplot) .eq. 0) then
					!Record potential energy total to use for output later
					if (rij2 .lt. wca_cut2) then
						potenergymol_LJ(molnoi) = potenergymol_LJ(molnoi)   &
					    + 4.d0*(invrij2**6.d0 - invrij2**3.d0 + 0.25d0) - eps 
						potenergymol_LJ(molnoj) = potenergymol_LJ(molnoj)   &
					    + 4.d0*(invrij2**6.d0 - invrij2**3.d0 + 0.25d0) - eps 
					else
						potenergymol_LJ(molnoi) = potenergymol_LJ(molnoi)   &
						+ 0.5d0*eps*(cos(sod_a*rij2 + sod_b) - 1.d0)
						potenergymol_LJ(molnoj) = potenergymol_LJ(molnoj)   &
						+ 0.5d0*eps*(cos(sod_a*rij2 + sod_b) - 1.d0)
					end if
					!Virial expression used to obtain pressure
					virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
					virialmol(molnoj) = virialmol(molnoj) + accijmag*rij2
				endif
			endif
			current => old
			old => current%next                  !obtain next item in list
		enddo
	enddo

	!Total used with other potentials (e.g. FENE)
	potenergymol = potenergymol + potenergymol_LJ
	
	nullify(current)       	                     !no longer required
	nullify(old)                                 !no longer required


end subroutine simulation_compute_forces_Soddemann_neigbr_halfint

!========================================================================
!Compute Volume Averaged stress using all cells including halos

subroutine simulation_compute_rfbins(imin, imax, jmin, jmax, kmin, kmax)
use module_compute_forces
implicit none

	integer                         :: i, j, ixyz ,jxyz  !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp, cellsperbin
	integer							:: molnoi, molnoj, memloc
	integer							:: minbin, maxbin, imin, jmin, kmin, imax, jmax, kmax
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	!rfbin = 0.d0
	rijsum = 0.d0

	!Calculate bin to cell ratio
	cellsperbin = ceiling(ncells(1)/dble(nbins(1)))

	do kcell=(kmin-1)*cellsperbin+1, kmax*cellsperbin
	do jcell=(jmin-1)*cellsperbin+1, jmax*cellsperbin
	do icell=(imin-1)*cellsperbin+1, imax*cellsperbin
	
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp					!Step through each particle in list 
			molnoi = oldi%molno 	 	!Number of molecule
			ri = r(molnoi,:)         	!Retrieve ri

			do kcellshift = -1,1
			do jcellshift = -1,1
			do icellshift = -1,1

				!Prevents out of range values in i
				if (icell+icellshift .lt. imin) cycle
				if (icell+icellshift .gt. imax) cycle
				!Prevents out of range values in j
				if (jcell+jcellshift .lt. jmin) cycle
				if (jcell+jcellshift .gt. jmax) cycle
				!Prevents out of range values in k
				if (kcell+kcellshift .lt. kmin) cycle
				if (kcell+kcellshift .gt. kmax) cycle

				oldj => cell%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
				adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

				do j = 1,adjacentcellnp          !Step through all j for each i

					molnoj = oldj%molno 	 !Number of molecule
					rj = r(molnoj,:)         !Retrieve rj

					currentj => oldj
					oldj => currentj%next    !Use pointer in datatype to obtain next item in list

					if(molnoi==molnoj) cycle !Check to prevent interaction with self

					rij2=0                   !Set rij^2 to zero
					rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j

					do ixyz=1,nd
						rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
					enddo

					if (rij2 < rcutoff2) then
						!Add current distance to rijsum for molecules i and j
						rijsum(molnoi,:) = rijsum(molnoi,:) + 0.5d0*rij(:)
						rijsum(molnoj,:) = rijsum(molnoj,:) + 0.5d0*rij(:)
						!Linear magnitude of acceleration for each molecule
						invrij2 = 1.d0/rij2                 !Invert value
						accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)
						call pressure_tensor_forces_VA(ri,rj,rij,accijmag)
						!call pressure_tensor_forces_H(ri,rj,rij,accijmag)
						!if (vflux_outflag.eq.4) then
						!	fij = accijmag*rij(:)
							!call Control_Volume_Forces(fij,ri,rj,molnoi,molnoj)
						!	call control_volume_stresses(fij,ri,rj,molnoi,molnoj)
						!endif
					endif
				enddo
			enddo
			enddo
			enddo
			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo
	enddo
	enddo

	nullify(oldi)      	!Nullify as no longer required
	nullify(oldj)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

end subroutine simulation_compute_rfbins
