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

#if USE_COUPLER

!--------------------------------------------------------------------------------------
!Apply force to prevent molecules leaving domain using form suggested by Nie, Chen and
!Robbins (2004)

subroutine simulation_apply_boundary_forces
!	use module_external_forces
!	use coupler, only : coupler_md_boundary_forces

	implicit none

!	call coupler_md_boundary_forces(np,pressure,r,a)
    call top_boundary_constraint_force

contains

    subroutine top_boundary_constraint_force
        use coupler
        use interfaces
        use calculated_properties_MD, only : pressure
        use computational_constants_MD, only : nh, ncells, cellsidelength, domain, halfdomain, delta_rneighbr
        use linked_list, only : cell, node
        use arrays_MD, only : r, a
        implicit none

        integer i, j, k, js, je, ip, m
        real(kind(0.d0)) dy, y2, y3, yc, p
        type(node), pointer :: current => null()
        logical :: overlap, firsttime = .true.
        save :: y2, y3, dy, js, je, firsttime
        
        call coupler_md_get(overlap_with_top_cfd=overlap)

        if(.not. overlap) return

        if (firsttime) then
            firsttime = .false.
            call coupler_md_get(top_dy=dy)
            
            y2 = halfdomain(2) - dy
            y3 = halfdomain(2)
            
            ! get the range of j cell index in y direction
            js = ceiling(y2/cellsidelength(2)) + nh
            je = min(ncells(2),ceiling(domain(2)/cellsidelength(2))) + nh

            write(0,*) 'boundary constraint force ', dy, y2, y3, js, je, ncells, (js-nh)*cellsidelength(2)

            ! check if molecules from below cells can get in constraint region ( heuristic condition )
            if ( y2 - (js - nh - 1) * cellsidelength(2) < delta_rneighbr ) then  
               js = js - 1
            endif

            ! sanity check
			if ( y2 > (js - nh) * cellsidelength(2)  ) then  
				write(0,*) "wrong value of js in top_boundary_constraint_force", js
			!    call error_abort("wrong value of js in top_boundary_constraint_force", js)
			endif

        endif

        p = max(pressure,1.d0)

        do k = nh + 1, nh + 1 + ncells(3) 
            do j = js, je
                do i = nh + 1, nh + 1 + ncells(1) 
                    current => cell%head(i,j,k)%point
                    do ip = 1, cell%cellnp(i,j,k) 
                        m = current%molno
                        yc = r(2,m)
                        if ( y2 <= yc .and. yc < y3 ) then
                            a(2,m)= a(2,m) - p*(yc-y2)/(1.d0-(yc-y2)/(y3-y2))
                        end if
                        current => current%next
                    end do
                end do
            end do
        end do
                        
    end subroutine top_boundary_constraint_force

end subroutine simulation_apply_boundary_forces

#endif

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
				isumvel = isumvel + v(1,molno) 	!Add streamwise velocity to current bin
				isumacc = isumacc + a(1,molno) 	!Add acceleration to current bin
				current => old
				old => current%next 
			enddo

			!Reset pointer to head of bin
			old => bin%head(ibin,jbin,kbin)%point
			
			!Apply coupling force as Nie, Chen and Robbins (2004), using
			!Linear extrapolation of velocity
			do n = 1, binnp    ! Loop over all particles
				molno = old%molno !Number of molecule
				a(1,molno)= a(1,molno)-isumacc/binnp - & 
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

	a(ixyz,:) = a(ixyz,:) + F_const

end subroutine simulation_apply_constant_force
