!-----------------------------------------------------------------------------
!
!                            Compute_forces                               
! Establish the forces applied to the particles as a result of interactions
! using the lennard jones potential in diemsnionless units
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

	double precision                :: rij2        !rij dotted with itself (one dimensional)
	double precision                :: invrij2     !inverse of rij2 
	double precision                :: accijmag    !Non directional component of acceleration
	double precision,dimension(3)   :: ri, rj !Position of molecule i and j
	double precision,dimension(3)   :: rij    !vector between particles i and j

end module module_compute_forces
!========================================================================
!Force calculation using re-ordered neighbour list for all interactions

subroutine simulation_compute_forces
use module_compute_forces
implicit none

	integer                         :: j, ixyz   !Define dummy index
	integer				:: molnoi, molnoj
	integer				:: noneighbrs
	type(neighbrnode), pointer	:: old, current

	a	    	= 0.d0	!Reset acceleration matrix before force calculations
	potenergymol	= 0.d0	!Reset potential energy per molecule before calculation
	potenergysum	= 0.d0	!Reset potential energy sum before calculation
	virial 	    	= 0.d0	!Reset virial sum before calculation
	virialmol   	= 0.d0	!Reset virial sum per molecule before calculation

	do molnoi = 1, np

	        noneighbrs = neighbour%noneighbrs(molnoi)!Determine number of elements in neighbourlist
		old => neighbour%head(molnoi)%point	 !Set old to head of neighbour list
		ri(:) = r(molnoi,:)		!Retrieve ri

		do j = 1,noneighbrs			!Step through all pairs of neighbours i and j

			molnoj = old%molnoj		!Number of molecule j
			rj(:) = r(molnoj,:)		!Retrieve rj
			rij(:) = ri(:) - rj(:)   	!Evaluate distance between particle i and j

			rij2=0                   !Set rij^2 to zero
			do ixyz=1,nd
				rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
			enddo

			if (rij2 < rcutoff2) then

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

				!Only calculate properties when required for output
				if (mod(iter,tplot) .eq. 0) then

					!Record potential energy total to use for output later (potshift=-1 for WCA)
					potenergymol(molnoi)=potenergymol(molnoi) & 
						     +4.d0*(invrij2**6-invrij2**3)-potshift
					potenergymol(molnoj)=potenergymol(molnoj) & 
						     +4.d0*(invrij2**6-invrij2**3)-potshift

					!Virial expression used to obtain pressure
					virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
					virialmol(molnoj) = virialmol(molnoj) + accijmag*rij2

					!Factor of two for half interaction virial calculation
					if (pressure_outflag .eq. 1) then
						call pressure_tensor_forces(molnoi,rij,accijmag)
						if (molnoj .le. np) call pressure_tensor_forces(molnoj,rij,accijmag)
					endif
					!if (pressure_outflag .eq. 2) call pressure_tensor_forces_VA(ri,rj,rij,accijmag)
					if (pressure_outflag .eq. 3) call pressure_tensor_forces_MOP(2,ri(:),rj(:),rij(:),2*accijmag)

				endif
			endif
			current => old
			old => current%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine simulation_compute_forces
