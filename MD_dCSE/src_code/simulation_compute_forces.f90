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
!----------------------------------------------------------------------------------
!========================================================================
!Force calculation using re-ordered neighbour list for all interactions

subroutine simulation_compute_forces
use module_compute_forces
implicit none

	integer                         :: i, j, ixyz   !Define dummy index
	integer				:: molnoi, molnoj
	integer				:: noneighbrs
	type(neighbrnode), pointer	:: old, current

	a=0.d0            !Reset acceleration matrix before force calculations
	potenergymol=0.d0 !Reset potential energy per molecule before calculation
	potenergysum=0.d0 !Reset potential energy sum before calculation
	virial = 0.d0     !Reset virial sum before calculation
	virialmol = 0.d0  !Reset virial sum per molecule before calculation

	do molnoi = 1, np

	        noneighbrs = neighbour%noneighbrs(molnoi)	!Determine number of elements in neighbourlist
		old => neighbour%head(molnoi)%point		!Set old to head of neighbour list

		do j = 1,noneighbrs			!Step through all pairs of neighbours molnoi and j

			ri(:) = r(molnoi,:)			!Retrieve ri
			molnoj = old%molnoj		!Number of molecule j
			rj(:) = r(molnoj,:)		!Retrieve rj
			rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j

			rij2=0                   !Set rij^2 to zero
			do ixyz=1,nd
				rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
			enddo

			if (rij2 < rcutoff2) then
				invrij2 = 1.d0/rij2                 !Invert value
				!Linear magnitude of acceleration for each molecule
				accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)

				!Sum of forces on particle i added for each j
				a(molnoi,1)= a(molnoi,1) + accijmag*rij(1)
				a(molnoi,2)= a(molnoi,2) + accijmag*rij(2)
				a(molnoi,3)= a(molnoi,3) + accijmag*rij(3)

				!Only calculate properties when required for output
				if (mod(iter,tplot) .eq. 0) then
					!Record potential energy total to use for output later (potshift=-1 for WCA)
					potenergymol(molnoi)=potenergymol(molnoi) & 
						     +4.d0*(invrij2**6-invrij2**3)-potshift
					!Virial expression used to obtain pressure
					virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
					if (pressure_outflag .eq. 1) then
						if (pressure_outflag .eq. 1) call pressure_tensor_forces(molnoi,rij,accijmag)
						if (pressure_outflag .eq. 2) then 
							write(2,'(a)') 'VOLUME AVERAGE NOT AVAILABLE WITH ALL INTERACTIONS -using MOP'
							pressure_outflag = 3
						endif
						if (pressure_outflag .eq. 3) call pressure_tensor_forces_MOP(2,ri(:),rj(:),rij(:),accijmag)
					endif
				endif
			endif
			current => old
			old => current%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine simulation_compute_forces

!========================================================================
!Compute forces using only half the interactions

subroutine simulation_compute_forces_halfint
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

				!call Control_Volume_Forces(accijmag*rij(:),ri,rj,molnoi,molnoj)
				if (vflux_outflag .ne. 0) call Control_Volume_stresses(2.d0*accijmag*rij(:),ri,rj,molnoi,molnoj)

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

	!Calculation of stress tensor Pxy

	!if (mod(iter,tplot) .eq. 0) then
	!	if (pressure_outflag .eq. 2) call  simulation_compute_rfbins(1,nbins(1)+2,1,nbins(1)+2,1,nbins(1)+2)
	!endif

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine simulation_compute_forces_halfint


!========================================================================
!Compute forces using cells instead of neighbourlist

subroutine simulation_compute_forces_cells
use module_compute_forces
implicit none

	integer                         :: i, j, ixyz   !Define dummy index
	integer				:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp
	integer				:: molnoi, molnoj, memloc
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	a=0.d0            !Reset acceleration matrix before force calculations
	potenergymol=0.d0 !Reset potential energy per molecule before calculation
	potenergysum=0.d0 !Reset potential energy sum before calculation
	virial = 0.d0     !Reset virial sum before calculation
	virialmol = 0.d0  !Reset virial sum per molecule before calculation

	do kcell=2, ncells(1)+1
	do jcell=2, ncells(2)+1
	do icell=2, ncells(3)+1

		!print*, icell, jcell, kcell

		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp                  !Step through each particle in list 
			molnoi = oldi%molno 	 !Number of molecule
			ri = r(molnoi,:)         !Retrieve ri

			do kcellshift = -1,1
			do jcellshift = -1,1
			do icellshift = -1,1
				oldj => cell%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
				adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

				!print*, icell+icellshift,jcell+jcellshift,kcell+kcellshift

				do j = 1,adjacentcellnp          !Step through all j for each i

					molnoj = oldj%molno 	 !Number of molecule
					rj = r(molnoj,:)         !Retrieve rj

					!do ixyz = 1,3
						!print'(a,7i,f10.5)', 'CPU', &
						!icell+icellshift,jcell+jcellshift,kcell+kcellshift,i,j,ixyz,molnoj,rj(ixyz)
					!enddo

					currentj => oldj
					oldj => currentj%next    !Use pointer in datatype to obtain next item in list

					if(molnoi==molnoj) cycle !Check to prevent interaction with self

					rij2=0                   !Set rij^2 to zero
					rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j

					do ixyz=1,nd
						rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
					enddo

					if (rij2 < rcutoff2) then

						!print'(2(i,3f10.5))',i,ri,j,rj

						invrij2 = 1.d0/rij2                 !Invert value
						!Linear magnitude of acceleration for each molecule
						accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)
	
						!Sum of forces on particle i added for each j
						a(molnoi,1)= a(molnoi,1) + accijmag*rij(1)
						a(molnoi,2)= a(molnoi,2) + accijmag*rij(2)
						a(molnoi,3)= a(molnoi,3) + accijmag*rij(3)

						!Only calculate properties when required for output
						if (mod(iter,tplot) .eq. 0) then
							!Record potential energy total to use for output later (potshift=-1 for WCA)
							potenergymol(molnoi)=potenergymol(molnoi) & 
								     +4.d0*(invrij2**6-invrij2**3)-potshift
							!Virial expression used to obtain pressure
							virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
							if (pressure_outflag .eq. 1) call pressure_tensor_forces(molnoi,rij,accijmag)
							if (pressure_outflag .eq. 2) & 
								call pressure_tensor_forces_VA(ri,rj,rij,accijmag)
							if (pressure_outflag .eq. 3) & 
								call pressure_tensor_forces_MOP(2,ri(:),rj(:),rij(:),accijmag)
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

end subroutine simulation_compute_forces_cells

!========================================================================
!Compute Volume Averaged stress using all cells including halos

subroutine simulation_compute_rfbins(imin, imax, jmin, jmax, kmin, kmax)
use module_compute_forces
implicit none

	integer                         :: i, j, ixyz ,jxyz  !Define dummy index
	integer				:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp, cellsperbin
	integer				:: molnoi, molnoj, memloc
	integer				:: minbin, maxbin, imin, jmin, kmin, imax, jmax, kmax
	type(node), pointer 	        :: oldi, currenti, oldj, currentj


	rfbin = 0.d0
	rijsum = 0.d0

	!Calculate bin to cell ratio
	cellsperbin = ceiling(ncells(1)/dble(nbins(1)))

	do kcell=(kmin-1)*cellsperbin+1, kmax*cellsperbin
	do jcell=(jmin-1)*cellsperbin+1, jmax*cellsperbin
	do icell=(imin-1)*cellsperbin+1, imax*cellsperbin
	
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp                  !Step through each particle in list 
			molnoi = oldi%molno 	 !Number of molecule
			ri = r(molnoi,:)         !Retrieve ri

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
