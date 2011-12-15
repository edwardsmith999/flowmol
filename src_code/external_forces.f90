!				EXTERNAL FORCES
! APPLY EXTERNAL FORCES FOR THERMOSTATTING AND CONSTRAINING EQUATIONS
!
! --Applied Forces--
! subroutine simulation_apply_constraint_forces()
! subroutine simulation_apply_linear_forces()
! subroutine simulation_apply_continuum_forces()
!
! --Thermostats--
! subroutine vthermostat_move()
! subroutine NHthermostat()
! 
!--------------------------------------------------------------------------------------

module module_external_forces

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use calculated_properties_MD
	use linked_list

end module module_external_forces

!--------------------------------------------------------------------------------------
!Apply force to prevent molecules leaving domain using form suggested by Nie, Chen and
!Robbins (2004)

subroutine simulation_apply_constraint_forces
use module_external_forces
use coupler, only : coupler_md_constrain_forces
implicit none

! call the coupler version.
! This is needed because the coupler has the information about the continuum grid

	pressure = 2.5d0

	call coupler_md_constrain_forces(np,pressure,r,a)

!!$	integer	:: n
!!$	double precision :: delta_y
!!$	double precision :: y, y0, y1, y2, y3
!!$
!!$	if (jblock .eq. npy) then
!!$
!!$		delta_y = globaldomain(2)/6.d0
!!$	
!!$		y0 = 0.d0
!!$		y1 = y0 + delta_y
!!$		y2 = y1 + delta_y
!!$		y3 = globaldomain(2)/2.d0
!!$	
!!$		do n = 1, np
!!$			if (r(n,2) .gt. y2) then
!!$		
!!$				y = r(n,2)
!!$				!Initial pressure is -ve - this line prevents problems
!!$				if (pressure .lt. 0.d0) then
!!$					a(n,2)= a(n,2) - (y-y2)/(1-(y-y2)/(y3-y2))
!!$				else
!!$					a(n,2)= a(n,2) - (y-y2)/(1-(y-y2)/(y3-y2))*pressure
!!$				endif
!!$			endif
!!$		enddo
!!$
!!$	endif
	
end subroutine simulation_apply_constraint_forces

!--------------------------------------------------------------------------------------
!Apply force to give linear profile

subroutine simulation_apply_linear_forces
use module_external_forces
implicit none

	!double precision :: y, y0, y1, y2, y3
	integer         			:: n, molno, ixyz
	integer         			:: cbin, binnp
	integer					:: ibin, jbin, kbin
	integer					:: averagecount
	double precision 			:: F, fixdist, slicebinsize
	double precision			:: isumvel, isumacc
	double precision, dimension(overlap)	:: continuum_u
	type(node), pointer:: old, current

	if (jblock .eq. npy) then

		!fixdist = 2.d0

		F = 0.1d0 !uwall/(domain(2)-fixdist)

		!slicebinsize = domain(2) / nbins(1)
		!y0 = halfdomain(2)-3*slicebinsize
		!y1 = halfdomain(2)-2*slicebinsize
		!y2 = halfdomain(2)-1*slicebinsize
		!y3 = halfdomain(2)

		!Apply force to top three bins in y
		do ibin=1, 1
		do jbin=nbins(1)-overlap+1, nbins(1)	!Top bins 
		do kbin=1, 1
			binnp = bin%cellnp(ibin,jbin,kbin)
			old => bin%head(ibin,jbin,kbin)%point
			cbin = jbin - nbins(1) + overlap !cbin runs from 1 to overlap, jbin from 1 to nbins(1)

			continuum_u(cbin) = 0.1d0*F*cbin + 2.d0*F

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
				current => old
				old => current%next 
			enddo
		enddo
		enddo
		enddo
	endif

end subroutine simulation_apply_linear_forces


!--------------------------------------------------------------------------------------
!Apply a force with the same value to all molecules

subroutine simulation_apply_constant_force(ixyz,F_const)
use module_external_forces
implicit none

	integer         			:: ixyz
	double precision 			:: F_const

	a(:,ixyz) = a(:,ixyz) + F_const

end subroutine simulation_apply_constant_force
