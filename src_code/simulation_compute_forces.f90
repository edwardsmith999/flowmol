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
	
	a					= 0.d0							!Initialise forces to zero
	potenergymol		= 0.d0
	potenergymol_LJ		= 0.d0
	potenergysum		= 0.d0
	potenergysum_LJ		= 0.d0
	virial				= 0.d0
	virialmol			= 0.d0

	if (potential_flag.eq.1) then
		potenergymol_FENE	= 0.d0
		potenergysum_FENE	= 0.d0
	end if

	if (potential_flag.eq.0) then					!If simple LJ fluid
		call simulation_compute_forces_LJ_halfint	!Compute LJ forces
	else if (potential_flag.eq.1) then				!If FENE polymer
		call simulation_compute_forces_LJ_halfint	!Compute LJ bead interactions
		call simulation_compute_forces_FENE			!Add on FENE spring interactions
	else									
		stop 'Potential flag not recognised!'
	end if

end subroutine simulation_compute_forces

subroutine simulation_compute_forces_LJ
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

				!Only calculate properties when required for output
				if (mod(iter,tplot) .eq. 0) then
					!Record potential energy total to use for output later (potshift=-1 for WCA)
					potenergymol_LJ(molnoi)=potenergymol_LJ(molnoi) & 
						     +4.d0*(invrij2**6-invrij2**3)-potshift
					potenergymol(molnoi) = potenergymol_LJ(molnoi)
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

end subroutine simulation_compute_forces_LJ

subroutine simulation_compute_forces_FENE
use module_compute_forces
use polymer_info_MD
implicit none

	integer				:: molnoi								!Current LJ bead
	integer				:: molnoL								!Bead on left
	integer				:: molnoR								!Bead on right
	
	do molnoi=1,np
		
		ri(:) = r(molnoi,:)							!Retrieve ri(:)
		molnoL = polyinfo_mol(molnoi)%left
		molnoR = polyinfo_mol(molnoi)%right
		
		if (molnoL.ne.0) then						!If there is a bead connected to the left of i

			rj(:)  = r(molnoL,:)					
			rij(:) = ri(:) - rj(:)
			rij2 = dot_product(rij,rij)
			if(rij2.ge.R_0**2) then
				print '(a,i6,a,i4,a,i4,a,f8.5,a,f8.5,a,f12.5)', & 
                                              'Bond broken at iter ',iter,': atoms ',molnoi,' and ',molnoL,' are separated by ', &
					      rij2**0.5,', which is greater than the allowed limit of ', R_0, &
                                              '. Stopping simulation, total time elapsed = ', iter*delta_t
				print '(a)', 'Atomic positions:'
				print '(a,i4,a,f10.5,a,f10.5,a,f10.5)', 'Atom ',molnoi,' is located at ',r(molnoi,1),' ',r(molnoi,2),' ',r(molnoi,3) 
				print '(a,i4,a,f10.5,a,f10.5,a,f10.5)', 'Atom ',molnoL,' is located at ',r(molnoL,1),' ',r(molnoL,2),' ',r(molnoL,3) 
				if (molnoL.gt.np) print*, 'Halo!'
				stop
			end if
	
			accijmag = -k_c/(1-(rij2/(R_0**2)))			!(-dU/dr)*(1/|r|)
			a(molnoi,1)= a(molnoi,1) + accijmag*rij(1)	!Add components of acceleration
			a(molnoi,2)= a(molnoi,2) + accijmag*rij(2)
			a(molnoi,3)= a(molnoi,3) + accijmag*rij(3)
			
			if (mod(iter,tplot) .eq. 0) then
				potenergymol_FENE(molnoi)=potenergymol_FENE(molnoi)-0.5d0*k_c*R_0*R_0*dlog(1.d0-(rij2/(R_0**2)))
				potenergymol(molnoi) = potenergymol(molnoi) + potenergymol_FENE(molnoi)
				virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
			endif
	
		end if
		
		if (molnoR.ne.0) then

			rj(:)  = r(molnoR,:)	
			rij(:) = ri(:) - rj(:)
			rij2 = dot_product(rij,rij)
			if(rij2.ge.R_0**2) then
				print '(a,i6,a,i4,a,i4,a,f8.5,a,f8.5,a,f12.5)',&
                                      'Bond broken at iter ',iter,': atoms ',molnoi,' and ',molnoR,' are separated by ', &
				      rij2**0.5,', which is greater than the allowed limit of ', R_0,&
                                      '. Stopping simulation, total time elapsed = ', iter*delta_t
				print '(a)', 'Atomic positions:'
				print '(a,i4,a,f10.5,a,f10.5,a,f10.5)', 'Atom ',molnoi,' is located at ',r(molnoi,1),' ',r(molnoi,2),' ',r(molnoi,3) 
				print '(a,i4,a,f10.5,a,f10.5,a,f10.5)', 'Atom ',molnoR,' is located at ',r(molnoR,1),' ',r(molnoR,2),' ',r(molnoR,3) 
				if (molnoR.gt.np) print*, 'Halo!'
				stop
			end if
			accijmag = -k_c/(1-(rij2/(R_0**2)))
			a(molnoi,1)= a(molnoi,1) + accijmag*rij(1)
			a(molnoi,2)= a(molnoi,2) + accijmag*rij(2)
			a(molnoi,3)= a(molnoi,3) + accijmag*rij(3)
	
			if (mod(iter,tplot) .eq. 0) then
				potenergymol_FENE(molnoi)=potenergymol_FENE(molnoi)-0.5d0*k_c*R_0*R_0*dlog(1.d0-(rij2/(R_0**2)))
				potenergymol(molnoi) = potenergymol(molnoi) + potenergymol_FENE(molnoi)
				virialmol(molnoi) = virialmol(molnoi) + accijmag*rij2
			endif
		end if	
		
	end do

end subroutine simulation_compute_forces_FENE


!========================================================================
!Compute forces using only half the interactions

subroutine simulation_compute_forces_LJ_halfint
use module_compute_forces
implicit none

	integer                         :: j, ixyz   !Define dummy index
	integer							:: molnoi, molnoj
	integer							:: noneighbrs
	type(neighbrnode), pointer		:: old, current

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
					potenergymol_LJ(molnoi)=potenergymol_LJ(molnoi)+4.d0*(invrij2**6-invrij2**3)-potshift
					potenergymol_LJ(molnoj)=potenergymol_LJ(molnoj)+4.d0*(invrij2**6-invrij2**3)-potshift

					potenergymol(molnoi) = potenergymol(molnoi) + 4.d0*(invrij2**6-invrij2**3) - potshift
					potenergymol(molnoj) = potenergymol(molnoj) + 4.d0*(invrij2**6-invrij2**3) - potshift
					
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

end subroutine simulation_compute_forces_LJ_halfint


!========================================================================
!Compute forces using cells instead of neighbourlist

subroutine simulation_compute_forces_cells
use module_compute_forces
implicit none

	integer                         :: i, j, ixyz   !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp
	integer							:: molnoi, molnoj, memloc
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	a=0.d0            !Reset acceleration matrix before force calculations
	potenergymol=0.d0 !Reset potential energy per molecule before calculation
	potenergysum=0.d0 !Reset potential energy sum before calculation
	virial = 0.d0     !Reset virial sum before calculation
	virialmol = 0.d0  !Reset virial sum per molecule before calculation

	do kcell=2, ncells(3)+1
	do jcell=2, ncells(2)+1
	do icell=2, ncells(1)+1

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
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp, cellsperbin
	integer							:: molnoi, molnoj, memloc
	integer							:: minbin, maxbin, imin, jmin, kmin, imax, jmax, kmax
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
