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

contains

	!========================================================================
	!Cell list computations of potential and force on "would-be" molecules
	subroutine compute_force_and_potential_at(input_pos,Usum,f,extra_pos) 
		use linked_list, only : node, cell
		use physical_constants_MD, only : rcutoff2, potshift
		use computational_constants_MD, only: halfdomain, cellsidelength, nh, ncells
		use arrays_MD, only : r
		implicit none

		real(kind(0.d0)),dimension(3), intent(in)	:: input_pos
		real(kind(0.d0)), intent(out) 				:: Usum
		real(kind(0.d0)),dimension(3), intent(out)	:: f
		!Optional array of extra molecular positions to check against
		double precision,dimension(:,:),allocatable,optional,intent(in)  :: extra_pos

		integer :: i,j 
		integer :: icell, jcell, kcell
		integer :: icellshift, jcellshift, kcellshift
		integer :: cellnp
		integer :: molno
		type(node), pointer :: current, temp
		real(kind(0.d0)) :: fmol(3), Umol, rij(3), rij2,invrij2

		!print*, 'compute_force_and_potential_at', present(extra_pos)

		! Init	
		Usum = 0.d0
		f = 0.d0 

		!Find cell, adding nh for halo(s)
		icell = ceiling((input_pos(1)+halfdomain(1))/cellsidelength(1)) + nh
		jcell = ceiling((input_pos(2)+halfdomain(2))/cellsidelength(2)) + nh
		kcell = ceiling((input_pos(3)+halfdomain(3))/cellsidelength(3)) + nh

		!Return Usum and f zero if position is outside the domain
		if ( icell .lt. 2 .or. icell .gt. ncells(1)+1 .or. &
			 jcell .lt. 2 .or. jcell .gt. ncells(2)+1 .or. &
			 kcell .lt. 2 .or. kcell .gt. ncells(3)+1      ) then
			!print*, 'Warning - attempted to calculated force and potential'
			!print*, 'outside of the domain. Returning Usum=f=0.'
			return
		end if

		do kcellshift = -1,1
		do jcellshift = -1,1
		do icellshift = -1,1

			current => cell%head  (icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
			cellnp  =  cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)
		
			do j = 1,cellnp	

				rij(:) = input_pos(:) - r(:,current%molno)
				rij2   = dot_product(rij,rij)

		        if (rij2 .eq. 0.d0) then
		            !print*, 'Warning, computing potential with zero separation...', current%molno
					current => current%next
		            cycle
		        end if

		        !Linear magnitude of acceleration for each molecule
		        invrij2 = 1.d0/rij2

				if (rij2 < rcutoff2) then

					!Linear magnitude of acceleration for each molecule
					invrij2 = 1.d0/rij2

					!Find molecule's contribution to f and Usum
					fmol = 48.d0*( invrij2**7 - 0.5d0*invrij2**4 )*rij
					Umol = 4.d0*( invrij2**6 - invrij2**3 )-potshift

					!Add to totals
					f = f + fmol
					Usum = Usum + Umol

				endif

				current => current%next

			enddo

		enddo
		enddo
		enddo

		nullify(current)      	!Nullify as no longer required


		! If array of extra_positions is present, loop through and compare
		! current molecule to them

		if (present(extra_pos)) then
			do j = 1,size(extra_pos,2)

				rij(:) = input_pos(:) - extra_pos(:,j)
				rij2   = dot_product(rij,rij)

		        if (rij2 .eq. 0.d0) then
		            !print*, 'Warning, computing potential with zero separation...', current%molno
					current => current%next
		            cycle
		        end if

		        !Linear magnitude of acceleration for each molecule
		        invrij2 = 1.d0/rij2

				if (rij2 < rcutoff2) then

					!Linear magnitude of acceleration for each molecule
					invrij2 = 1.d0/rij2

					!Find molecule's contribution to f and Usum
					fmol = 48.d0*( invrij2**7 - 0.5d0*invrij2**4 )*rij
					Umol = 4.d0*( invrij2**6 - invrij2**3 )-potshift

					!Add to totals
					f = f + fmol
					Usum = Usum + Umol

				endif

			enddo
		endif


	end subroutine compute_force_and_potential_at 

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
			call error_abort("Potential flag/force_list incompatible - only LJ available with full interaction neighbour lists")
			!call simulation_compute_forces_LJ_neigbr
            !call simulation_compute_forces_FENE					!Add on FENE spring interactions
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
	case(-1)
		!Do Nothing -- don't calculate the force
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
		ri = r(:,i)         		!Retrieve ri
		do j = i+1,np				!Step through all j for each i
			rj = r(:,j)				!Retrieve rj

			!Calculate rij using nearest image convention, i.e. if more than half
			! a domain betwen molecules, must be closer over periodic boundaries  
			rij2 = 0.d0
			do ixyz=1,nd
				rij(ixyz) = r (ixyz,i) - r(ixyz,j)          !Evaluate distance between particle i and j
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
				a(1,i)= a(1,i) + accijmag*rij(1)
				a(2,i)= a(2,i) + accijmag*rij(2)
				a(3,i)= a(3,i) + accijmag*rij(3)

				!Sum of forces on particle j added for each i
				a(1,j)= a(1,j) - accijmag*rij(1)
				a(2,j)= a(2,j) - accijmag*rij(2)
				a(3,j)= a(3,j) - accijmag*rij(3)

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
				else if (mod(iter,teval) .eq. 0) then
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

	integer                         :: i,j
	integer                         :: p_i, p_j, ptot
	double precision, parameter     :: sod_a    = 3.1730728678
	double precision, parameter     :: sod_b    = -0.85622864544
	double precision, parameter     :: wca_cut  = 1.12246204830937
!	double precision, parameter     :: wca_cut2 = 1.25992104989486
	double precision, parameter     :: wca_cut2 = 1.25992104989487
	double precision                :: eps

	do i = 1,np
		ri = r(:,i)
		do j = i+1,np                !Step through all pairs

			rj     = r(:,j)
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
				a(1,i)= a(1,i) + accijmag*rij(1)
				a(2,i)= a(2,i) + accijmag*rij(2)
				a(3,i)= a(3,i) + accijmag*rij(3) 

				!Sum of forces on particle j added for each i
				a(1,j)= a(1,j) - accijmag*rij(1)
				a(2,j)= a(2,j) - accijmag*rij(2)
				a(3,j)= a(3,j) - accijmag*rij(3) 

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
				else if (mod(iter,teval) .eq. 0) then
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

	integer                         :: i,j,repeatno !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp
	integer							:: molnoi, molnoj
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	do kcell=2, ncells(3)+1
	do jcell=2, ncells(2)+1
	do icell=2, ncells(1)+1

		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp					!Step through each particle in list 

			molnoi = oldi%molno 	 	!Number of molecule
			ri = r(:,molnoi)         	!Retrieve ri
            
			do kcellshift = -1,1
			do jcellshift = -1,1
			do icellshift = -1,1

				oldj => cell%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
				adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

				do j = 1,adjacentcellnp			!Step through all j for each i

					molnoj = oldj%molno			!Number of molecule
					rj = r(:,molnoj)			!Retrieve rj

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
						a(1,molnoi)= a(1,molnoi) + accijmag*rij(1)
						a(2,molnoi)= a(2,molnoi) + accijmag*rij(2)
						a(3,molnoi)= a(3,molnoi) + accijmag*rij(3)

						!CV stress and force calculations
						if (vflux_outflag .eq. 4 .or. eflux_outflag.eq.4) then
							if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
								fij = accijmag*rij(:)
								call Control_Volume_stresses(fij,ri,rj)
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
						else if (mod(iter,teval) .eq. 0) then
							!Virial expression used to obtain pressure
							virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
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

end subroutine simulation_compute_forces_LJ_cells

!-----------------------------------------------------------------

subroutine simulation_compute_forces_LJ_neigbr
	use module_compute_forces
	implicit none

	integer                         :: j  !Define dummy index
	integer							:: molnoi, molnoj
	integer							:: noneighbrs
	type(neighbrnode), pointer		:: old, current

	do molnoi = 1, np

        noneighbrs = neighbour%noneighbrs(molnoi)	!Determine number of elements in neighbourlist
		old => neighbour%head(molnoi)%point			!Set old to head of neighbour list
		ri(:) = r(:,molnoi)							!Retrieve ri

		do j = 1,noneighbrs							!Step through all pairs of neighbours molnoi and j

			molnoj = old%molnoj						!Number of molecule j
			rj(:) = r(:,molnoj)						!Retrieve rj
			rij(:) = ri(:) - rj(:)   				!Evaluate distance between particle i and j
			rij2 = dot_product(rij,rij)				!Square of vector calculated

			if (rij2 < rcutoff2) then
				invrij2  = 1.d0/rij2                !Invert value
				accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4) ! (-dU/dr)*(1/|r|)

				!Sum of forces on particle i added for each j
				a(1,molnoi)= a(1,molnoi) + accijmag*rij(1)
				a(2,molnoi)= a(2,molnoi) + accijmag*rij(2)
				a(3,molnoi)= a(3,molnoi) + accijmag*rij(3)

				!CV stress and force calculations
				if (vflux_outflag .eq. 4 .or. eflux_outflag.eq.4) then
					if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
						fij(:) = accijmag*rij(:)
						call control_volume_stresses(fij,ri,rj)
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
				else if (mod(iter,teval) .eq. 0) then
					!Virial expression used to obtain pressure
					virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
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
	!use CV_objects, only : CV_sphere_momentum
	implicit none

	integer                         :: j, ixyz 
	integer							:: molnoi, molnoj
	integer							:: noneighbrs
	type(neighbrnode), pointer		:: old, current

	do molnoi = 1, np

	    noneighbrs = neighbour%noneighbrs(molnoi)	!Determine number of elements in neighbourlist
		old => neighbour%head(molnoi)%point			!Set old to head of neighbour list
		ri(:) = r(:,molnoi)							!Retrieve ri

		do j = 1,noneighbrs							!Step through all pairs of neighbours i and j

			molnoj = old%molnoj			!Number of molecule j
			rj(:) = r(:,molnoj)			!Retrieve rj
			rij(:)= ri(:) - rj(:)   	!Evaluate distance between particle i and j
			rij2  = dot_product(rij,rij)!Square of vector calculated

			if (rij2 .lt. rcutoff2) then

				!Linear magnitude of acceleration for each molecule
				invrij2  = 1.d0/rij2                 !Invert value
				accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)
	
				!Sum of forces on particle i added for each j
				a(1,molnoi)= a(1,molnoi) + accijmag*rij(1)
				a(2,molnoi)= a(2,molnoi) + accijmag*rij(2)
				a(3,molnoi)= a(3,molnoi) + accijmag*rij(3) 

				!Sum of forces on particle j added for each i
				a(1,molnoj)= a(1,molnoj) - accijmag*rij(1)
				a(2,molnoj)= a(2,molnoj) - accijmag*rij(2)
				a(3,molnoj)= a(3,molnoj) - accijmag*rij(3) 

				if (vflux_outflag.eq.4 .or. eflux_outflag.eq.4) then
					if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
						if (molnoj .gt. np .or. molnoi .gt. np) then
							fij = accijmag*rij(:)
							call control_volume_stresses(fij,ri,rj)
						    !call CV_sphere_momentum%Add_spherical_CV_forces(fij,ri,rj)
						else
							fij = 2.d0*accijmag*rij(:)
							call control_volume_stresses(fij,ri,rj)
							!call CV_sphere_momentum%Add_spherical_CV_forces(fij,ri,rj)
						endif
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
					if (vflux_outflag.ne.0 .and. vflux_outflag.ne.4)	then
						call pressure_tensor_forces_MOP(vflux_outflag,ri(:),rj(:),rij(:),accijmag)
					endif

				else if (mod(iter,teval) .eq. 0) then

					!Virial expression used to obtain pressure
					virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
					virialmol(molnoj) = virialmol(molnoj) + accijmag*rij2

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

		ri(:) = r(:,molnoi)							!Retrieve ri(:)

		do b=1,monomer(molnoi)%funcy

			molnoj = bond(b,molnoi)
			if (molnoj.eq.0) cycle

			rj(:)  = r(:,molnoj)
			rij(:) = ri(:) - rj(:)
			rij2   = dot_product(rij,rij)

			if(rij2.ge.R_0**2)	call polymer_bond_error(molnoj)

			accijmag = -k_c/(1-(rij2/(R_0**2)))			!(-dU/dr)*(1/|r|)

			a(1,molnoi)= a(1,molnoi) + accijmag*rij(1)	!Add components of acceleration
			a(2,molnoi)= a(2,molnoi) + accijmag*rij(2)
			a(3,molnoi)= a(3,molnoi) + accijmag*rij(3)

			if (mod(iter,tplot) .eq. 0) then
				if (pressure_outflag .eq. 1) then
					call pressure_tensor_forces(molnoi,rij,accijmag)
				endif
				potenergymol_FENE(molnoi) = potenergymol_FENE(molnoi) - 0.5d0*k_c*R_0*R_0*dlog(1.d0-(rij2/(R_0**2)))
				potenergymol(molnoi)      = potenergymol(molnoi)      + potenergymol_FENE(molnoi)
				virialmol(molnoi)         = virialmol(molnoi)         + accijmag*rij2

			else if (mod(iter,teval) .eq. 0) then
			
				virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2

			endif

		end do	

	end do

contains

	subroutine polymer_bond_error(molnoX)
		use interfaces
		implicit none

		integer, intent(in) :: molnoX
		real(kind(0.d0)) :: rglobi(3),rglobX(3)

		rglobi(1) = r(1,molnoi) - halfdomain(1)*(npx-1) + domain(1)*(iblock - 1)   
		rglobi(2) = r(2,molnoi) - halfdomain(2)*(npy-1) + domain(2)*(jblock - 1)   
		rglobi(3) = r(3,molnoi) - halfdomain(3)*(npz-1) + domain(3)*(kblock - 1)   

		rglobX(1) = r(1,molnoX) - halfdomain(1)*(npx-1) + domain(1)*(iblock - 1)   
		rglobX(2) = r(2,molnoX) - halfdomain(2)*(npy-1) + domain(2)*(jblock - 1)   
		rglobX(3) = r(3,molnoX) - halfdomain(3)*(npz-1) + domain(3)*(kblock - 1)   

		print*, 'irank: ', irank
		print '(a,i6,a,i8,a,i8,a,f8.5,a,f8.5,a,f12.5)', & 
					'Bond broken at iter ',iter,': atoms ',molnoi,' and ',molnoX,' are separated by ', &
					rij2**0.5,', which is greater than the allowed limit of ', R_0, &
					'. Stopping simulation, total time elapsed = ', iter*delta_t
		print '(a)', 'Atomic positions:'
		print '(a,i8,a,f10.5,a,f10.5,a,f10.5)', 'Atom ',molnoi,' is located at global position', &
		                                        rglobi(1),' ',rglobi(2),' ',rglobi(3) 
		print '(a,i8,a,f10.5,a,f10.5,a,f10.5)', 'Atom ',molnoX,' is located at global position', &
		                                        rglobX(1),' ',rglobX(2),' ',rglobX(3) 

		print '(a,i8,a)', 'Monomer information for atom ', molnoi,':'
		print '(a,i8)', 'ChainID: '   , monomer(molnoi)%chainID
		print '(a,i8)', 'SubchainID: ', monomer(molnoi)%subchainID
		print '(a,i8)', 'Funcy: '     , monomer(molnoi)%funcy
		print '(a,i8)', 'Glob_no: '   , monomer(molnoi)%glob_no
		print '(a,i8)', 'Bin_bflag: ' , monomer(molnoi)%bin_bflag
		print '(a,i8)', 'Tag: '       , tag(molnoi) 

		print '(a,i8,a)', 'Monomer information for atom ', molnoX,':'
		print '(a,i8)', 'ChainID: '   , monomer(molnoX)%chainID
		print '(a,i8)', 'SubchainID: ', monomer(molnoX)%subchainID
		print '(a,i8)', 'Funcy: '     , monomer(molnoX)%funcy
		print '(a,i8)', 'Glob_no: '   , monomer(molnoX)%glob_no
		print '(a,i8)', 'Bin_bflag: ' , monomer(molnoX)%bin_bflag
		print '(a,i8)', 'Tag: '       , tag(molnoX) 

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
		ri(:) = r(:,molnoi)							!Retrieve ri

		do j = 1,noneighbrs                         !Step through all pairs
		                                            !of neighbours i and j
			molnoj = old%molnoj                     !Number of molecule j
			rj     = r(:,molnoj)                    !Position
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
				a(1,molnoi)= a(1,molnoi) + accijmag*rij(1)
				a(2,molnoi)= a(2,molnoi) + accijmag*rij(2)
				a(3,molnoi)= a(3,molnoi) + accijmag*rij(3) 

				!Sum of forces on particle j added for each i
				a(1,molnoj)= a(1,molnoj) - accijmag*rij(1)
				a(2,molnoj)= a(2,molnoj) - accijmag*rij(2)
				a(3,molnoj)= a(3,molnoj) - accijmag*rij(3) 

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

				else if (mod(iter,teval) .eq. 0 ) then
				
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

	integer,intent(in)  			:: imin, jmin, kmin, imax, jmax, kmax

	integer                         :: i, j, ixyz !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp 
	integer							:: molnoi, molnoj
	integer							:: ibinmin,jbinmin,kbinmin,ibinmax,jbinmax,kbinmax
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	double precision,dimension(3)	:: cellsperbin

	!rfbin = 0.d0
	!allocate(rijsum(nd,np+extralloc)) !Sum of rij for each i, used for SLLOD algorithm
	!rijsum = 0.d0

	!Calculate bin to cell ratio
	cellsperbin = 1.d0/binspercell !ceiling(ncells(1)/dble(nbins(1)))
    ! Still need to loop over every cell (i.e. get all interactions) if
    ! bins are bigger than cells
	where (cellsperbin .lt. 1.d0) cellsperbin = 1.d0

!	do kcell=(kmin-1)*cellsperbin(3)+1, kmax*cellsperbin(3)
!	do jcell=(jmin-1)*cellsperbin(2)+1, jmax*cellsperbin(2)
!	do icell=(imin-1)*cellsperbin(1)+1, imax*cellsperbin(1)

	ibinmin = (imin-1)*cellsperbin(1)+2+(1-cellsperbin(1))
	ibinmax =  imax   *cellsperbin(1)-cellsperbin(1)
	jbinmin = (jmin-1)*cellsperbin(2)+2+(1-cellsperbin(2))
	jbinmax =  jmax   *cellsperbin(2)-cellsperbin(2)
	kbinmin = (kmin-1)*cellsperbin(3)+2+(1-cellsperbin(3))
	kbinmax =  kmax   *cellsperbin(3)-cellsperbin(3)

	do kcell=kbinmin, kbinmax
	do jcell=jbinmin, jbinmax 
	do icell=ibinmin, ibinmax 
	
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp					!Step through each particle in list 
			molnoi = oldi%molno 	 	!Number of molecule
			ri = r(:,molnoi)         	!Retrieve ri

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
					rj = r(:,molnoj)         !Retrieve rj

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
						!rijsum(:,molnoi) = rijsum(:,molnoi) + 0.5d0*rij(:)
						!rijsum(:,molnoj) = rijsum(:,molnoj) + 0.5d0*rij(:)
						!Linear magnitude of acceleration for each molecule
						invrij2 = 1.d0/rij2                 !Invert value
						accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)

						!Select requested configurational line partition methodology
						select case(VA_calcmethod)
						case(0)
							call pressure_tensor_forces_H(ri,rj,rij,accijmag)
						case(1)
							call pressure_tensor_forces_VA_trap(ri,rj,accijmag)
						case(2)
							call pressure_tensor_forces_VA(ri,rj,rij,accijmag)
						case default
							stop "Error - VA_calcmethod incorrect"
						end select 

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

    ! Add FENE contribution if it's there
    if (potential_flag .eq. 1) then
        call add_FENE_contribution
    end if

	nullify(oldi)      	!Nullify as no longer required
	nullify(oldj)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

contains

    subroutine add_FENE_contribution
        use polymer_info_MD
        implicit none

        integer :: b

        do molnoi=1,np

            ri(:) = r(:,molnoi) !Retrieve ri(:)
            do b=1,monomer(molnoi)%funcy

                molnoj = bond(b,molnoi)
                if (molnoj.eq.0) cycle

                rj(:)  = r(:,molnoj)
                rij(:) = ri(:) - rj(:)
                rij2   = dot_product(rij,rij)
                accijmag = -k_c/(1-(rij2/(R_0**2)))	!(-dU/dr)*(1/|r|)

                call pressure_tensor_forces_VA_trap(ri,rj,accijmag)

            end do	

        end do

    end subroutine

end subroutine simulation_compute_rfbins

!Cylindrical polar version
subroutine simulation_compute_rfbins_cpol(imin, imax, jmin, jmax, kmin, kmax)
	use module_compute_forces
	use librarymod, only: cpolariser
	use messenger, only: globalise
	implicit none

	integer, intent(in) :: imin, jmin, kmin, imax, jmax, kmax

	integer :: i, j
	integer	:: icell, jcell, kcell
	integer :: icellshift, jcellshift, kcellshift
	integer :: cellnp, adjacentcellnp
	integer	:: molnoi, molnoj
	type(node), pointer :: oldi, currenti, oldj, currentj

	double precision,dimension(3)	:: cellsperbin

	!Calculate bin to cell ratio
	cellsperbin = 1.d0/binspercell !ceiling(ncells(1)/dble(nbins(1)))
	where (cellsperbin .gt. 1.d0) cellsperbin = 1.d0

	do kcell=(kmin-1)*cellsperbin(3)+1, kmax*cellsperbin(3)
	do jcell=(jmin-1)*cellsperbin(2)+1, jmax*cellsperbin(2)
	do icell=(imin-1)*cellsperbin(1)+1, imax*cellsperbin(1)
	
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point

		do i = 1,cellnp

			molnoi = oldi%molno
			ri = r(:,molnoi)

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

				oldj => cell%head(icell+icellshift, &
				                  jcell+jcellshift, &
				                  kcell+kcellshift) % point

				adjacentcellnp = cell%cellnp(icell+icellshift, &
				                             jcell+jcellshift, &
				                             kcell+kcellshift)

				do j = 1,adjacentcellnp

					molnoj = oldj%molno
					rj = r(:,molnoj)

					currentj => oldj
					oldj => currentj%next 

					! Prevent interaction with self
					if ( molnoi == molnoj ) cycle 

					rij(:) = ri(:) - rj(:)
					rij2 = dot_product(rij,rij)

					if (rij2 < rcutoff2) then

						invrij2 = 1.d0/rij2
						accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)

						call pressure_tensor_forces_VA_trap_cpol(ri,rj,accijmag)

					endif

				enddo

			enddo
			enddo
			enddo

			currenti => oldi
			oldi => currenti%next

		enddo

	enddo
	enddo
	enddo

    ! Add FENE contribution if it's there
    if (potential_flag .eq. 1) then
        call add_FENE_contribution
    end if

	nullify(oldi)
	nullify(oldj)
	nullify(currenti)
	nullify(currentj)

contains

    subroutine add_FENE_contribution
        use polymer_info_MD
        implicit none

        integer :: b

        do molnoi=1,np

            ri(:) = r(:,molnoi) !Retrieve ri(:)
            do b=1,monomer(molnoi)%funcy

                molnoj = bond(b,molnoi)
                if (molnoj.eq.0) cycle

                rj(:)  = r(:,molnoj)
                rij(:) = ri(:) - rj(:)
                rij2   = dot_product(rij,rij)
                accijmag = -k_c/(1-(rij2/(R_0**2)))	!(-dU/dr)*(1/|r|)

                call pressure_tensor_forces_VA_trap_cpol(ri,rj,accijmag)

            end do	

        end do

    end subroutine

end subroutine simulation_compute_rfbins_cpol



!========================================================================
! Compute Energy for specified range of CV, requires current acceleration
! to already be calculated for all molecules so must be called after
! simulation computer forces

subroutine simulation_compute_power(imin, imax, jmin, jmax, kmin, kmax)
	use module_compute_forces
	implicit none

	integer,intent(in)				:: imin, jmin, kmin, imax, jmax, kmax


	integer                         :: i, j, ixyz !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp 
	integer							:: molnoi, molnoj
	integer							:: ibinmin,jbinmin,kbinmin,ibinmax,jbinmax,kbinmax
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	double precision,dimension(3)	:: cellsperbin

	!Calculate bin to cell ratio
	cellsperbin = 1.d0/binspercell !ceiling(ncells(1)/dble(nbins(1)))

    ! Still need to loop over every cell (i.e. get all interactions) if
    ! bins are bigger than cells
	where (cellsperbin .lt. 1.d0) cellsperbin = 1.d0

!	do kcell=(kmin-1)*cellsperbin(3)+1, kmax*cellsperbin(3)
!	do jcell=(jmin-1)*cellsperbin(2)+1, jmax*cellsperbin(2)
!	do icell=(imin-1)*cellsperbin(1)+1, imax*cellsperbin(1)

	ibinmin = (imin-1)*cellsperbin(1)+2+(1-cellsperbin(1))
	ibinmax =  imax   *cellsperbin(1)-cellsperbin(1)
	jbinmin = (jmin-1)*cellsperbin(2)+2+(1-cellsperbin(2))
	jbinmax =  jmax   *cellsperbin(2)-cellsperbin(2)
	kbinmin = (kmin-1)*cellsperbin(3)+2+(1-cellsperbin(3))
	kbinmax =  kmax   *cellsperbin(3)-cellsperbin(3)

	do kcell=kbinmin, kbinmax
	do jcell=jbinmin, jbinmax 
	do icell=ibinmin, ibinmax 

		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp					!Step through each particle in list 
			molnoi = oldi%molno 	 	!Number of molecule
			ri = r(:,molnoi)         	!Retrieve ri

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
					rj = r(:,molnoj)         !Retrieve rj

					currentj => oldj
					oldj => currentj%next    !Use pointer in datatype to obtain next item in list

					if(molnoi==molnoj) cycle !Check to prevent interaction with self

					rij2=0                   !Set rij^2 to zero
					rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j

					do ixyz=1,nd
						rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
					enddo

					if (rij2 < rcutoff2) then

						!Linear magnitude of acceleration for each molecule
						invrij2 = 1.d0/rij2                 !Invert value
						accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)

						!CV stress and force calculations
						if (eflux_outflag.eq. 4) then
							if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
								fij = accijmag*rij(:)
								call control_volume_power(fij,ri,rj,molnoi,a(:,molnoi))
							endif
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

	nullify(oldi)      	!Nullify as no longer required
	nullify(oldj)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

end subroutine simulation_compute_power

subroutine collect_bforce_pdf_data
    use boundary_MD
    use physical_constants_MD, only: np, rcutoff2
    use computational_constants_MD, only: cellsidelength, nh, halfdomain, ncells
    use linked_list, only: neighbrnode, neighbour, cell, node
    use librarymod, only: heaviside, normal_dist
    use arrays_MD, only: r
	implicit none

	integer :: i, j, ixyz 
	integer :: molnoi, molnoj
	integer :: noneighbrs
    integer :: ycell_i, ycell_j, ysubcell

    integer :: icell, jcell, kcell, icellshift, jcellshift, kcellshift
    integer :: cellnp, adjacentcellnp
    type(node), pointer :: oldi, oldj, currenti, currentj

    real(kind(0.d0)) :: ri(3), rj(3), rij(3), rij2, invrij2, accijmag
    real(kind(0.d0)) :: bforce(3)
    logical :: bflag
	type(neighbrnode), pointer :: old, current

	do kcell=2, ncells(3)+1
	do jcell=2, ncells(2)+1
	do icell=2, ncells(1)+1

		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point

		do i = 1,cellnp					!Step through each particle in list 

            bforce = 0.d0
            bflag = .true.

			molnoi = oldi%molno 	 	!Number of molecule
			ri = r(:,molnoi)         	!Retrieve ri
            
			do kcellshift = -1,1
			do jcellshift = -1,1
			do icellshift = -1,1

				oldj => cell%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
				adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

				do j = 1,adjacentcellnp			!Step through all j for each i

					molnoj = oldj%molno			!Number of molecule
					rj = r(:,molnoj)			!Retrieve rj

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

                        ycell_i = ceiling((r(2,molnoi)+halfdomain(2)) &
                        /cellsidelength(2))!+nh !Add nh due to halo(s)
                        ycell_j = ceiling((r(2,molnoj)+halfdomain(2)) &
                        /cellsidelength(2))!+nh !Add nh due to halo(s)
                       
                        if (ycell_i .lt. ycell_j) then
                            bflag = .true.
                            bforce(1) = bforce(1) + accijmag*rij(1)
                            bforce(2) = bforce(2) + accijmag*rij(2)
                            bforce(3) = bforce(3) + accijmag*rij(3)
                        end if 

					endif

				enddo

			enddo
			enddo
			enddo

            if (bflag) then
                ysubcell = ceiling((real(bforce_pdf_nsubcells,kind(0.d0))*( &
                           r(2,molnoi) + halfdomain(2)))/cellsidelength(2))
                ysubcell = mod(ysubcell, bforce_pdf_nsubcells)
                if (ysubcell .eq. 0) then
                    ysubcell = bforce_pdf_nsubcells 
                end if

                do ixyz = 1,3
                    call bforce_pdf(ixyz, ysubcell)%update( &
                                    (/bforce(ixyz)/))
                end do

            end if

			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list

		enddo

	enddo
	enddo
	enddo

!	do molnoi = 1, np
!
!        virtual_plane_force = 0.d0
!
!	    noneighbrs = neighbour%noneighbrs(molnoi)	!Determine number of elements in neighbourlist
!		old => neighbour%head(molnoi)%point			!Set old to head of neighbour list
!		ri(:) = r(:,molnoi)							!Retrieve ri
!
!		do j = 1,noneighbrs							!Step through all pairs of neighbours i and j
!
!			molnoj = old%molnoj			!Number of molecule j
!			rj(:) = r(:,molnoj)			!Retrieve rj
!			rij(:)= ri(:) - rj(:)   	!Evaluate distance between particle i and j
!			rij2  = dot_product(rij,rij)!Square of vector calculated
!
!			if (rij2 .lt. rcutoff2) then
!
!				!Linear magnitude of acceleration for each molecule
!				invrij2  = 1.d0/rij2                 !Invert value
!				accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)
!
!                ycell_i = ceiling((r(2,molnoi)+halfdomain(2)) &
!                /cellsidelength(2))!+nh !Add nh due to halo(s)
!                ycell_j = ceiling((r(2,molnoj)+halfdomain(2)) &
!                /cellsidelength(2))!+nh !Add nh due to halo(s)
!               
!                if (ycell_i .lt. ycell_j) then
!                !if (ycell_i .ne. ycell_j) then
!                    virtual_plane_force(:) = virtual_plane_force(:) &
!                      + accijmag*rij(:)*heaviside(ycell_j-ycell_i)!  & 
!                      !- accijmag*rij(:)*heaviside(ycell_i-ycell_j)
!                end if 
!
!			endif
!
!			current => old
!			old => current%next !Use pointer in datatype to obtain next item in list
!
!		enddo
!
!        ysubcell = ceiling((real(bforce_pdf_nsubcells,kind(0.d0))*( &
!                   r(2,molnoi) + halfdomain(2)))/cellsidelength(2))
!        ysubcell = mod(ysubcell, bforce_pdf_nsubcells)
!        if (ysubcell .eq. 0) then
!            ysubcell = bforce_pdf_nsubcells 
!        end if
!    
!        !print'(a,f12.5,a,f12.5,a,f12.5, a, f12.5, a, i8, a, i8)', &
!        !'yL/2: ', halfdomain(2), &
!        !'  r2 + yL/2: ', r(2,molnoi) + halfdomain(2), &
!        !'  subcellratio1: ',(real(bforce_pdf_nsubcells,kind(0.d0))*(r(2,molnoi) + halfdomain(2)))/cellsidelength(2), &
!        !'  cellratio1: ', (r(2,molnoi) + halfdomain(2))/cellsidelength(2), &
!        !'  ycell_i ', ycell_i, &
!        !'  ysubcell: ', ysubcell
!
!        do ixyz = 1,3
!            !print*, virtual_plane_force
!            call bforce_pdf(ixyz, ysubcell)%update( &
!                            (/virtual_plane_force(ixyz)/))
!            !call bforce_pdf(ixyz,ysubcell)%update((/normal_dist()/))
!        end do
!
!	enddo

	nullify(current)
	nullify(old)

end subroutine collect_bforce_pdf_data 

!========================================================================
!Cell list computations of potential and force on "would-be" molecules
subroutine compute_force_and_potential_at(input_pos,Usum,f) 
	use module_compute_forces
	implicit none

	real(kind(0.d0)), intent(in)  :: input_pos(3)
	real(kind(0.d0)), intent(out) :: Usum, f(3)

	integer :: i,j 
	integer :: icell, jcell, kcell
	integer :: icellshift, jcellshift, kcellshift
	integer :: cellnp
	integer :: molno
	type(node), pointer :: current, temp
	real(kind(0.d0)) :: fmol(3), Umol

	! Init	
	Usum = 0.d0
	f = 0.d0 

	!Find cell, adding nh for halo(s)
    icell = ceiling((input_pos(1)+halfdomain(1))/cellsidelength(1)) + nh
    jcell = ceiling((input_pos(2)+halfdomain(2))/cellsidelength(2)) + nh
    kcell = ceiling((input_pos(3)+halfdomain(3))/cellsidelength(3)) + nh

	!Return Usum and f zero if position is outside the domain
	if ( icell .lt. 2 .or. icell .gt. ncells(1)+1 .or. &
	     jcell .lt. 2 .or. jcell .gt. ncells(2)+1 .or. &
	     kcell .lt. 2 .or. kcell .gt. ncells(3)+1      ) then
		print*, 'Warning - attempted to calculated force and potential'
		print*, 'outside of the domain. Returning Usum=f=0.'
		return
	end if

	do kcellshift = -1,1
	do jcellshift = -1,1
	do icellshift = -1,1

		current => cell%head  (icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
		cellnp  =  cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)
    
		do j = 1,cellnp	

			rij(:) = input_pos(:) - r(:,current%molno)
			rij2   = dot_product(rij,rij)

            if (rij2 .eq. 0.d0) then
                !print*, 'Warning, computing potential with zero separation...', current%molno
			    current => current%next
                cycle
            end if

            !Linear magnitude of acceleration for each molecule
            invrij2 = 1.d0/rij2

			if (rij2 < rcutoff2) then

				!Linear magnitude of acceleration for each molecule
				invrij2 = 1.d0/rij2

				!Find molecule's contribution to f and Usum
				fmol = 48.d0*( invrij2**7 - 0.5d0*invrij2**4 )*rij
				Umol = 4.d0*( invrij2**6 - invrij2**3 )-potshift

				!Add to totals
				f = f + fmol
				Usum = Usum + Umol

			endif

			current => current%next

		enddo

	enddo
	enddo
	enddo

	nullify(current)      	!Nullify as no longer required

end subroutine compute_force_and_potential_at 
