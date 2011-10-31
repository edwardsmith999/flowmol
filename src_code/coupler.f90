
module module_coupler

	use module_external_forces
	use grid_arrays
	use computational_constants
	use computational_constants_MD
	use calculated_properties_MD

	implicit none 

end module module_coupler
!==========================================================================

!--------------------------------------------------------------------------------------
! Apply force to match velocity to continuum based on Nie et al 2004 paper

subroutine simulation_apply_continuum_forces
	use module_coupler
	implicit none

	integer         			:: n, molno, ixyz, cellnp
	integer         			:: cbin, averagecount
	integer					:: icell, jcell, kcell
	integer					:: isummol
	double precision			:: isumvel, isumacc
	double precision			:: t_fract, average
	double precision, dimension(4)		:: continuum_u
	type(node), pointer 	        	:: old, current

	!Fraction of continuum timestep which has passed
	t_fract = dble((iter - initialstep)) / dble((Nsteps-initialstep))

	!Linear extrapolation between velocity at t and t+1 and save in continuum_u
	!Taken for cell at nx/2 and top domain cell (ny + 1) to (ny + 1) - overlap
	!N.B. domain runs from 2 to ny + 1 due to halos
	continuum_u(1) = uc_t_minus_1(nint(nx/2.d0),3)*(1.d0-t_fract)  &
		 	         + uc(nint(nx/2.d0),3)*      t_fract
	continuum_u(2) = uc_t_minus_1(nint(nx/2.d0),3)*(1.d0-t_fract)  &
		 	         + uc(nint(nx/2.d0),3)*      t_fract
	continuum_u(3) = uc_t_minus_1(nint(nx/2.d0),4)*(1.d0-t_fract)  &
		 	         + uc(nint(nx/2.d0),4)*      t_fract
	continuum_u(4) = uc_t_minus_1(nint(nx/2.d0),4)*(1.d0-t_fract)  &
		 	         + uc(nint(nx/2.d0),4)*      t_fract


	average = 0.d0
	averagecount = 0

	!Apply force to top three bins in y
	!ASSUME Cell same size as bins and one continuum cell is two MD cells
	do jcell= (ncells(2)+1)-3,(ncells(2)+1)			!Loop through 4 y cells in controlled region
 		cbin = jcell - (ncells(2)+1)+4			!Local no. of overlap cell from 1 to overlap

		!print'(3i8,4f20.7)', iter, jcell, cbin, continuum_u

		!Reset acceleration and velocity sums
		isummol = 0
		isumvel = 0.d0
		isumacc = 0.d0

		do icell=2 , ncells(1)+1	!Loop through all x cells
		do kcell=2 , ncells(3)+1	!Loop through all z cells

			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list

			!Calculate averages for bin
			do n = 1, cellnp    ! Loop over all particles
				molno = old%molno 	 !Number of molecule
				!Assign to bins using integer division
				isumvel = isumvel + v(molno,1) 	!Add streamwise velocity to current bin
				isumacc = isumacc + a(molno,1) 	!Add acceleration to current bin
				isummol = isummol + 1
				current => old
				old => current%next 
			enddo

		enddo
		enddo

		!Get average velocity and acceleration in bin
		if (isummol .ne. 0) then
			isumacc = isumacc/real(isummol,kind(0.d0))
		 	isumvel = isumvel/real(isummol,kind(0.d0))
		endif

		!print*, isumacc, isumvel, isummol

		do icell=2 , ncells(1)+1
		do kcell=2 , ncells(3)+1

			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list
			
			!Apply coupling force as Nie, Chen and Robbins (2004), using
			!Linear extrapolation of velocity
			do n = 1, cellnp    ! Loop over all particles
				molno = old%molno !Number of molecule

				a(molno,1)= a(molno,1) - isumacc   &
					    -(isumvel-continuum_u(cbin))/delta_t

				current => old
				old => current%next 

				if (molno .gt. np) stop "Force applied to halo molecules"

				!average = average - isumacc   &
				!	    -(isumvel-continuum_u(cbin))/delta_t
				!averagecount = averagecount + 1

			enddo

		enddo
		enddo
	enddo

	!print'(a,f10.5,a,f18.9)', 'MD_velocity ',sum(continuum_u(:))/3 , & 
	!			' average force applied to MD molecules ', average/(averagecount)

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine simulation_apply_continuum_forces

!-------------------------------------------------------------------------------------

subroutine MD_continuum_BC
	use module_coupler
	use module_setup_BC
	implicit none

	double precision	 	:: binvolume 
	double precision,dimension(3) 	:: slicebinsize

	slicebinsize(:) = domain(:) / nbins(:)
	binvolume = slicebinsize(1)*slicebinsize(2)*slicebinsize(3)

	select case(velocity_outflag)
	case(2)
		u_MD = sum(slice_momentum((ncells(2)+1)-2*overlap-2:(ncells(2)+1)-2*overlap-1,1)) & 
		      /sum(slice_mass(    (ncells(2)+1)-2*overlap-2:(ncells(2)+1)-2*overlap-1  ))
	case(4)
		!Calculate molecular velocity by averaging overlap region
		stop "BINS NOT CORRECT FOR CONTINUUM BC - check this first"
		u_MD = 	 sum(volume_momentum(2:nbins(1),(nbins(2)+1-overlap),2:nbins(3),1)) & 
		       /(sum(volume_mass(2:nbins(1),(nbins(2)+1-overlap),2:nbins(3)))*nbins(1)*nbins(3))
	case default
		stop "Error input for velocity MD_continuum BC coupling incorrect - MD velocity_outflag must be non-zero"
	end select

end subroutine MD_continuum_BC
