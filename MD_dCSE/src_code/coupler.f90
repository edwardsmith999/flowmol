subroutine coupler

!--------------------------------------------------------------------------------------
!Apply force to match velocity to continuum

subroutine simulation_apply_continuum_forces
use module_external_forces
use grid_arrays
use computational_constants
use computational_constants_MD
implicit none

	integer         			:: n, molno, ixyz
	integer         			:: cbin, binnp
	integer					:: ibin, jbin, kbin
	double precision			:: t_fract
	double precision			:: isumvel, isumacc
	double precision, dimension(overlap)	:: continuum_u
	type(node), pointer			:: old, current

	call cumulative_bin_velocity(2)

	!Fraction of continuum timestep which has passed
	t_fract = dble((iter - initialstep)) / dble((Nsteps-initialstep))

	!Linear extrapolation between velocity at t and t+1 and save in continuum_u
	!Taken for cell at nx/2 and top domain cell (ny + 1) to (ny + 1) - overlap
	!N.B. domain runs from 2 to ny + 1 due to halos
	do cbin = 1, overlap
		continuum_u(cbin) = uc_t_minus_1(nx/2.,ny+2-cbin)*(1.d0-t_fract)  &
			 + uc(nx/2.,ny+2-cbin) * t_fract
	enddo

	!Apply force to top three bins in y
	do ibin=1, 1
	do jbin=nbins(1)-overlap+1, nbins(1)	!Top bins 
	do kbin=1, 1
		binnp = bin%cellnp(ibin,jbin,kbin)
		old => bin%head(ibin,jbin,kbin)%point
		cbin = jbin - nbins(1) + overlap !cbin runs from 1 to overlap, jbin from 1 to nbins(1)

		isumvel = 0.d0
		isumacc = 0.d0

		!Calculate averages for bin
		do n = 1, binnp    ! Loop over all particles
			molno = old%molno !Number of molecule
			!Assign to bins using integer division
			isumvel = isumvel + v(molno,1) 	!Add streamwise velocity to current bin
			isumacc = isumacc + a(molno,1) 	!Add acceleration to current bin
			current => old
			old => current%next 
		enddo

		!Reset pointer to head of bin
		old => bin%head(ibin,jbin,kbin)%point
			
		!Apply coupling force as Nie, Chen and Robbins (2004), using
		!Linear extrapolation of velocity
		do n = 1, binnp    ! Loop over all particles
			molno = old%molno !Number of molecule
			a(molno,1)= a(molno,1)-isumacc/binnp - & 
				     (isumvel/binnp &
				      -continuum_u(cbin))/delta_t
			!average = -isumacc/binnp - & 
			!	     (isumvel/binnp &
			!	      -continuum_u(cbin))/delta_t
			!averagecount = averagecount + 1
			!print'(4f10.5)',a(molno,1), isumacc/binnp, (isumvel/binnp)/delta_t, continuum_u(cbin)/delta_t
			current => old
			old => current%next 
		enddo
	enddo
	enddo
	enddo

	!print'(a,f10.5,a,f10.5)', 'average applied force', average/averagecount, 'average acceleration', sum(a)/np

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required
	
end subroutine simulation_apply_continuum_forces


end subroutine coupler
