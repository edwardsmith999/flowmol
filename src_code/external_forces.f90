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
!Apply a force with the same value to all molecules
! ixyz - direction of force
! F_const - force magnitude

subroutine simulation_apply_global_force(ixyz,F_const)
	use module_external_forces
	implicit none

	integer,intent(in)     			:: ixyz
	double precision,intent(in)		:: F_const

	integer				 		 :: n
	double precision,dimension(3):: F_vector
	
	!Put directional results into a vector 
	F_vector = 0.d0
	F_vector(ixyz) = F_const

	do n=1,np

		if (any(tag(n) .eq. tether_tags)) cycle
		a(ixyz,n) = a(ixyz,n) + F_vector(ixyz)

		if (vflux_outflag .eq. 4) then
			if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
					!Add constant force to cells based on number of molecules in each
					call record_external_forces(F_vector(:),r(:,n))
			endif
		endif

	enddo

end subroutine simulation_apply_global_force

!--------------------------------------------------------------------------------------
! Apply a localised force at a specified location
! ixyz - direction of force
! F_const - force magnitude
! xmin,xmax,ymin,ymax,zmin,zmax - extents of forced region in global coordinates

subroutine simulation_apply_local_force(ixyz,F_const,xmin,xmax,ymin,ymax,zmin,zmax)
	use module_external_forces
	use messenger, only : localise
	implicit none

	integer, intent(in) 		 :: ixyz
	double precision, intent(in) :: F_const
	double precision, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax

	integer				 		 :: n
	double precision,dimension(3):: lmin,lmax, F_vector

	!Get local coordinates
	lmin = (/ xmin, ymin, zmin /)
	lmax = (/ xmax, ymax, zmax /)
	lmin = localise(lmin); lmax = localise(lmax)

	!Put directional results into a vector 
	F_vector = 0.d0
	F_vector(ixyz) = F_const

	do n=1,np
		if (r(1,n) .lt. lmin(1)) cycle
		if (r(1,n) .gt. lmax(1)) cycle
		if (r(2,n) .lt. lmin(2)) cycle
		if (r(2,n) .gt. lmax(2)) cycle
		if (r(3,n) .lt. lmin(3)) cycle
		if (r(3,n) .gt. lmax(3)) cycle
		if (any(tag(n) .eq. tether_tags)) cycle

		!print'(3i4,6f10.5)',iblock,jblock,kblock, xmin,lmin(1),r(1,n),lmax(1),xmax, F_const

		a(ixyz,n) = a(ixyz,n) + F_vector(ixyz)

		if (vflux_outflag .eq. 4) then
			if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
				call record_external_forces(F_vector(:),r(:,n))
			endif
		endif

	enddo

	!print'(15f7.3)', xmin,ymin,zmin,lmin,xmax,ymax,zmax,lmax,maxval(a(1,:)),maxval(a(2,:)),maxval(a(3,:))

end subroutine simulation_apply_local_force

!--------------------------------------------------------------------------------------
!Apply a force to prevent molecules from escaping the domain

subroutine apply_boundary_force
	use computational_constants_MD, only: bforce_flag, bforce_dxyz, &
	                                      bforce_off, bforce_NCER, bforce_OT, &
	                                      bforce_Flekkoy, bforce_rdf
	use interfaces, only: error_abort
	use computational_constants_MD, only: periodic
#if USE_COUPLER
	use md_coupler_socket, only: socket_get_constraint_info
#endif
	implicit none

#if USE_COUPLER

		integer :: constraint_algorithm,OT,NCER,Flekkoy,off

		call socket_get_constraint_info(constraint_algorithm,OT=OT, &
                                        NCER=NCER,Flekkoy=Flekkoy,off=off)
		if ( constraint_algorithm .eq. off ) then
			return
		else if ( constraint_algorithm .eq. OT ) then
			call error_abort("OT boundary force not yet implemented")
		else if ( constraint_algorithm .eq. NCER ) then
			call coupled_apply_boundary_force_NCER
			!call coupled_apply_boundary_force(bforce_flag,bforce_dxyz)
		else if ( constraint_algorithm .eq. Flekkoy ) then
			!call simulation_apply_boundary_force(bforce_flag,(/ 0 0 0 2.d0 0 0 /)) 
			!Flekkoy boundary force applied by constraint
			return
		else
			call error_abort("Unrecognised constraint algorithm flag")
		end if	
	
#else

		if (all(periodic.ne.0)) then
			return
		else
			call simulation_apply_boundary_force(bforce_flag,bforce_dxyz) 
		end if

#endif	

end subroutine apply_boundary_force

subroutine simulation_apply_boundary_force(flags,dists)
	use arrays_MD,  only: r,a
	use computational_constants_MD, only: iblock,jblock,kblock,npx,npy,npz, &
	                                      domain,irank
	use physical_constants_MD, only: np, procnp
	use calculated_properties_MD, only: pressure
	use interfaces, only: error_abort
	implicit none

	integer,          dimension(6), intent(in) :: flags
	real(kind(0.d0)), dimension(6), intent(in) :: dists 

	integer :: n,ixyz,flag,block(3),npxyz(3)
	real(kind(0.d0)) :: xyz,thresh,hdom
	real(kind(0.d0)), dimension(3) :: tops,bottoms

	tops    = (/domain(1)/2.d0 - dists(2), &
				domain(2)/2.d0 - dists(4), &
				domain(3)/2.d0 - dists(6)  /)

	bottoms = (/dists(1) - domain(1)/2.d0, &
				dists(3) - domain(2)/2.d0, &
				dists(5) - domain(3)/2.d0  /)

	block  = (/iblock,jblock,kblock/)
	npxyz  = (/npx,npy,npz/)

	do n = 1, np
	do ixyz = 1, 3
			
		if ( r(ixyz,n)   .gt. tops(ixyz)  .and. &
		     block(ixyz) .eq. npxyz(ixyz)       ) then
			
			xyz    = r(ixyz,n)
			thresh = tops(ixyz)
			hdom   = domain(ixyz)/2.d0
			flag   = flags(2*ixyz)

		else if ( r(ixyz,n)   .lt. bottoms(ixyz) .and. &
		          block(ixyz) .eq. 1                   ) then
				
			xyz    = r(ixyz,n)
			thresh = bottoms(ixyz)
			hdom   = -domain(ixyz)/2.d0
			flag   = flags(2*ixyz - 1)
		
		else

			cycle

		end if

		call apply_bforce(a(ixyz,n),xyz,thresh,hdom,flag)

	end do
	end do

end subroutine simulation_apply_boundary_force

subroutine apply_bforce(a_in,xyz,thresh,hdom,flag)
	use computational_constants_MD, only: bforce_NCER,bforce_off
	use calculated_properties_MD,   only: pressure
	use interfaces,                 only: error_abort
	implicit none

	real(kind(0.d0)), intent(inout) :: a_in	  !Accel of which to add bforce
	real(kind(0.d0)), intent(in) :: xyz       !Position, 1D
	real(kind(0.d0)), intent(in) :: thresh    !Threshold for bforce, 1D
	real(kind(0.d0)), intent(in) :: hdom      !Domain edge
	integer, intent(in) :: flag               !Type of bforce 

	real(kind(0.d0)) :: numer,denom,ratio,P
	character(128)   :: string

	select case ( flag )
	case ( bforce_off  )

		return

	case ( bforce_NCER )

		P     = max(pressure,1.d0)                !Account for negative pressure
		numer = xyz - thresh                      !+ve or -ve dependent on pos in domain
		denom = 1.d0 - (xyz-thresh)/(hdom-thresh) !denom always +ve	
		ratio = numer / denom	

		a_in  = a_in - ratio*P	

	case default

		string="MD uncoupled boundary force only developed for NCER case"
		call error_abort(string)

	end select

end subroutine apply_bforce


#if USE_COUPLER

!--------------------------------------------------------------------------------------
!Apply force to prevent molecules leaving domain using form suggested by Nie, Chen and
!Robbins (2004)

subroutine simulation_apply_boundary_forces
!	use module_external_forces
!	use coupler, only : coupler_md_boundary_forces
	implicit none

	!call coupler_md_boundary_forces(np,pressure,r,a)
	!call top_boundary_constraint_force
	!call simulation_apply_boundary_forces_NCER

contains

    subroutine top_boundary_constraint_force
        use interfaces
        use calculated_properties_MD, only : pressure
        use computational_constants_MD, only : nh, ncells, cellsidelength, domain, halfdomain, delta_rneighbr
        use linked_list, only : cell, node
        use arrays_MD, only : r, a
		use md_coupler_socket, only: socket_get_bottom_of_top_boundary, &
		                             socket_get_overlap_status
        implicit none

        integer i, j, k, js, je, ip, m
        real(kind(0.d0)) y2, y3, yc, p
        type(node), pointer :: current => null()
        logical :: overlap, firsttime = .true.
        save :: y2, y3, js, je, firsttime

		overlap = socket_get_overlap_status()
		if (.not. overlap) return
        
        if (firsttime) then
            firsttime = .false.

			y2 = socket_get_bottom_of_top_boundary() !todo better name
            y3 = halfdomain(2)
            
            ! get the range of j cell index in y direction
            js = ceiling(y2/cellsidelength(2)) + nh
            je = min(ncells(2),ceiling(domain(2)/cellsidelength(2))) + nh

            ! check if molecules from below cells can get in constraint region ( heuristic condition )
            if ( y2 - (js - nh - 1) * cellsidelength(2) < delta_rneighbr ) then  
               js = js - 1
            endif

            ! sanity check
			if ( y2 > (js - nh) * cellsidelength(2)  ) then  
           		print'(2(a,f10.5),5i6,f10.5)', "WARNING - cnst region from", & 
									 y2, 'to', y3, js, je, ncells, (js-nh)*cellsidelength(2)
			    call error_abort("THIS ROUTINE IS OBSOLETE")
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

!--------------------------------------------------------------------------------------
!Apply force to prevent molecules leaving domain using form suggested by O'Connell
!and Thompson (1995)
subroutine simulation_apply_boundary_forces_OT
	use module_external_forces
	use md_coupler_socket, only: socket_get_dy
	implicit none

	integer				:: n
	double precision 	:: alpha
	double precision 	:: y, y2, delta_y

	delta_y = socket_get_dy()

	if (jblock .eq. npy) then

		y2 = halfdomain(2) - delta_y
	
		!Set alpha and pressure to fixed value in O'Connell and Thompson
		alpha = 2
		pressure = 3.d0

		do n = 1, np
			if (r(n,2) .gt. y2) then
				y = r(n,2)
				a(n,2)= a(n,2) - alpha * pressure * density**(-2.d0/3.d0)*(y-y2)
			endif
		enddo

	endif
	
end subroutine simulation_apply_boundary_forces_OT


!--------------------------------------------------------------------------------------
!Apply force to prevent molecules leaving domain using form suggested by Nie, Chen and
!Robbins (2004)

subroutine coupled_apply_boundary_force_NCER
	use module_external_forces
	use md_coupler_socket, only: socket_get_dy
	implicit none

	integer	:: n
	double precision 	:: delta_y
	double precision 	:: y, y2, y3

	if (jblock .eq. npy) then

		delta_y = socket_get_dy()

		y3 = halfdomain(2)
		y2 = y3 - delta_y
	
		!Prevents negative pressure (during initialisation) from releasing molecules
        pressure = max(pressure,1.d0)

		do n = 1, np
			if (r(2,n) .gt. y2) then
				y = r(2,n)
				a(2,n)= a(2,n) - (y-y2)/(1-(y-y2)/(y3-y2))*pressure
			endif
		enddo

	endif

end subroutine coupled_apply_boundary_force_NCER 

!-------------------------------------------------------------------
! Apply boundary forces based on RDF from T. Werder et al. J Comp. 
! Phys 205 (2005) 373–390. NOTE molecules will still escape

subroutine simulation_apply_boundary_forces_Werder
	use module_external_forces
	use md_coupler_socket, only: socket_get_dy
	implicit none

	integer	:: n
	double precision 	:: delta_y
	double precision 	:: y, y2, y3

	if (jblock .eq. npy) then

		delta_y = socket_get_dy()

		y3 = halfdomain(2)
		y2 = y3 - delta_y
	
		do n = 1, np
			if (r(2,n) .gt. y2) then
				y = r(2,n)
				a(2,n)= a(2,n) + approx_RDF(y,y2)
			endif
		enddo

	endif

contains
	
	! An approximation for the Radial Distribution Function given by curve 
	! fitting to simulation data at density = 0.6 and tempterature=1.8 given in
	! Appendix A of T. Werder et al. J Comp. Phys 205 (2005) 373–390
	function approx_RDF(r,rmin) result(F)
		use librarymod, only : heaviside 

		real(kind=kind(0.d0)), intent(in)	:: r,rmin
		real(kind=kind(0.d0))				:: F, rhat

		!Define local coordinate 
		rhat = r - rmin

		F = ( 10.8007 + 0.860717*rhat    - 172.468* rhat**2 					& 
					  + 86.9134* rhat**3 - 140.214* rhat**4	 ) 					& 
					 *(heaviside(rhat-0.2975)-heaviside(rhat)		) 			&
		   +(-3621.30 + 44657.4* rhat    - 204844.0*rhat**2 					&
					  + 414123.0*rhat**3 - 311674.0*rhat**4  ) 					& 
					 *(heaviside(rhat-0.3475)-heaviside(rhat-0.2975)) 			&
		   +( 4331.63 - 45188.5* rhat    + 176236.0*rhat**2 					&
					  - 305157.0*rhat**3 + 198111.0*rhat**4	) 					&
					 *(heaviside(rhat-0.3975)-heaviside(rhat-0.3475)) 			&
		   +(-94.4796 + 576.282* rhat    - 1436.11* rhat**2 					&
					  + 1804.53* rhat**3 - 1133.47* rhat**4 + 283.244*rhat**5) 	& 
					 *(heaviside(rhat-1.0000)-heaviside(rhat-0.3975))
	end function

end subroutine simulation_apply_boundary_forces_Werder

#endif




!=============================================================================
! Testing routine to apply a generic Flekkøy (2004) force 
!-----------------------------------------------------------------------------
subroutine apply_flekkoy_test
	use physical_constants_MD, only : np
	use computational_constants_MD, only : globaldomain,halfdomain, iter, &
											irank, jblock, npy
	use linked_list
	use librarymod, only : get_new_fileunit
	implicit none

	integer					:: box_np,length, fileunit
	integer,allocatable 	:: list(:)
	real(kind(0.d0))		:: dx,dy,dz, box_average
	real(kind(0.d0))		:: shear_flekkoy,pressure_flekkoy

	!Only apply force on top processor
	if (jblock .ne. npy) return

	dx = globaldomain(1);	dy = 4.d0; dz = globaldomain(3)

	!Read Shear pressure from file...
	fileunit = get_new_fileunit()
	inquire(iolength=length) shear_flekkoy
	open(unit=fileunit,file='./couette_stress_analy',form='unformatted', access='direct',recl=length)
	read(fileunit,rec=iter) shear_flekkoy !Divided by the Reynolds number
	close(fileunit,status='keep')
	pressure_flekkoy = 4.d0
	
	!Get average over current cell and apply constraint forces
	call average_over_bin
	call apply_force

contains

!=============================================================================
! Average molecules in overlap region to obtain values for 
! constrained dynamics algorithms
!-----------------------------------------------------------------------------

subroutine average_over_bin
	use computational_constants_MD, only : nhb
	use arrays_MD, only : r, v, a
	implicit none

	integer				:: n
	double precision	:: g


	!Zero box averages
	box_np = 0
	box_average = 0.d0

	!find the maximum number of molecules and allocate a list array	   
    allocate(list(np))

	do n = 1,np

		if ( r(2,n) .gt. halfdomain(2)-dy .and. r(2,n) .lt. halfdomain(2)) then

			!Add molecule to overlap list
			box_np   =  box_np   + 1
			list(box_np) =  n
			g = flekkoy_gweight(r(2,n),halfdomain(2)-dy,halfdomain(2))
			box_average =  box_average  + g

			!write(8888,'(4i8,5f15.9)'), irank,iter,n,box_np,g,box_average,halfdomain(2)-dy,r(2,n),halfdomain(2)
			
		endif

	enddo

	call SubcommSum(box_average,1)
	call SubcommSum(box_average,3)

end subroutine average_over_bin

!=============================================================================
! Apply force to molecules in overlap region
!-----------------------------------------------------------------------------

subroutine apply_force
	use arrays_MD, only : r,v,a
	use physical_constants_MD, only : density
	implicit none

	integer					:: i, n
	real(kind=kind(0.d0)) 	:: g, gsum, dA, dV

	dA = dx*dz
	dV = dx*dy*dz

	!Loop over all molecules and apply constraint
	do i = 1, box_np

		n = list(i)
		g = flekkoy_gweight(r(2,n),halfdomain(2)-dy,halfdomain(2))

		!Gsum is replaced with the fixed value based on density and volume
		gsum = box_average
		if (gsum .eq. 0.d0) cycle

		a(1,n) = a(1,n) + (g/gsum) * dA * shear_flekkoy
		a(2,n) = a(2,n) + (g/gsum) * dA * pressure_flekkoy

		!write(7777,'(4i8,2f10.5,f18.12,5f10.5)'), irank, iter,n,box_np, shear_flekkoy, dA, g, gsum, halfdomain(2)-dy,r(2,n),halfdomain(2), a(1,n)

        if (g .ne. 0.d0) then
			if (iter .lt. 1000) then
				write(1234,'(i3,2i7,5f12.6)'),irank,iter,n, &
						 					  r(2,n),v(2,n),a(2,n),g, & 
											 (g/gsum) * dA * pressure_flekkoy
			endif
        endif
	enddo

    !write(99999,'(i2,i7,i7,2f10.2,f6.1,3f9.3,6f12.4)'), rank_world,iter,np_overlap,sum(box_average(:,:,:)%a(2)),  &
    !                gsumcheck, gratio, stress_cfd(:,2,ib,jb+jcmin_recv-extents(3),kb), ave_a,ave_a_consrnt

end subroutine apply_force

! -----------------------------------------------------------
! Function returns Flekkoy weighting for given y and max/min

function flekkoy_gweight(y,ymin,ymax) result (g)

	real(kind=kind(0.d0)), intent(in)	:: y, ymin, ymax
	real(kind=kind(0.d0))				:: g, L, yhat

	!Define local coordinate as const runs from 0 < y < L/2
	L = ymax - ymin
	yhat = y - ymin - 0.5*L

    !Sanity Check and exceptions
    if (yhat .lt. 0.d0) then
        g = 0
        return
    elseif (yhat .gt. 0.5*L) then
		stop " flekkoy_gweight error - input y cannot be greater than ymax"
    endif

	!Calculate weighting function
	g = 2*( 1/(L-2*yhat) - 1/L - 2*yhat/(L**2))

end function

end subroutine apply_flekkoy_test




!--------------------------------------------------------------------------------------
!Apply force to give linear profile

subroutine simulation_apply_linear_forces
	use module_external_forces
	implicit none

	!double precision :: y, y0, y1, y2, y3
	integer         			:: n, molno
	integer         			:: cbin, binnp
	integer						:: ibin, jbin, kbin
	double precision 			:: F
	double precision			:: isumvel, isumacc
	double precision, dimension(overlap)	:: continuum_u
	type(node), pointer			:: old, current

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

!=========================================================================
! Routine to test flux calculations agree by various methods

subroutine apply_CV_force(iter)
	use control_volume
	use computational_constants_MD, only : irank, jblock, npy, globaldomain, CVforce_flag
	use calculated_properties_MD, only : pressure
	use librarymod, only : get_new_fileunit
	use module_external_forces, only : np, momentum_flux, irank,nbins, nbinso,domain,delta_t,Nvflux_ave
	implicit none
	
	integer,intent(in)				:: iter

	logical							:: apply_CVforce
	integer							:: M, box_np, m_bin1, m_bin2
	integer,allocatable 			:: list(:)
	!real(kind(0.d0))				:: dx,dy,dz
	real(kind(0.d0)),dimension(3)	:: binsize,MD_Pi_dS,MD_rhouu_dS,F_constraint
	real(kind(0.d0)),dimension(3)	:: u_bin,F_bin,u_bin1, u_bin2, F_bin1, F_bin2
	real(kind(0.d0)),dimension(3)	:: CFD_Pi_dS,CFD_rhouu_dS

	if (CVforce_flag .eq. 0) return

	binsize = domain/nbins

	!dx = globaldomain(1);	dy = 4.d0; dz = globaldomain(3)

	!Get average over current cell and apply constraint forces
	!call get_continuum_values
	call get_test_values(CVforce_flag)

	!Exchange CV data ready to apply force
	call update_CV_halos

	!Only apply force on top processor
	if (jblock .ne. npy) return

	!Retreive CV data and calculate force to apply
	call average_over_bin

	!Apply the force
	!call apply_force
	call apply_force_tests(apply_CVforce)

contains

! Apply arbitary forces for testing purposes
subroutine get_test_values(flag)
	use physical_constants_MD, only : pi
	implicit none

	integer,intent(in)	:: flag

	integer				:: i,length,fileunit

	!When velocity is near zero, apply force from then on...
	!if (abs(CV%X(i,j,k,1)) .lt. 0.001) then
	if (iter .gt. 100) then !111 is the time force is close to zero
		apply_CVforce = .true.
		F_constraint = 0.d0
	endif

	select case(flag)
	case(0)
		!No Force applied
		apply_CVforce = .false.
	case(1)
		!Zero Force applied
		CFD_Pi_dS = 0.d0
		CFD_rhouu_dS = 0.d0	
	case(2)
		!Constant function 
		CFD_Pi_dS = 1.d0
		CFD_rhouu_dS = 0.d0
	case(3)
		!Sin function 
		CFD_Pi_dS = sin(2*pi*iter/100)*10.0
		CFD_rhouu_dS = 0.d0
	case(4)
		! Get continuum values of surface stresses, etc
		!Read Shear pressure from file...
		CFD_Pi_dS(:) = 0.d0
		fileunit = get_new_fileunit()
		inquire(iolength=length) CFD_Pi_dS(1)
		open(unit=fileunit,file='./F_hist',form='unformatted', access='direct',recl=length)
		read(fileunit,rec=iter) CFD_Pi_dS(1)
		close(fileunit,status='keep')
		CFD_rhouu_dS = 0.d0
	end select
	
end subroutine get_test_values

! Exchange all halo values for CV values ready to apply forces
subroutine update_CV_halos
	use CV_objects, only :  CV2 => CVcheck_momentum2
	use mpi
	implicit none

	call CV2%swap_halos(nbinso)

end subroutine update_CV_halos

!=============================================================================
! Average molecules in overlap region to obtain values for 
! constrained dynamics algorithms
!-----------------------------------------------------------------------------

subroutine average_over_bin
	use computational_constants_MD, only : nhb, iblock
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	use CV_objects, only : CVcheck_mass, CV  => CVcheck_momentum, & 
										 CV2 => CVcheck_momentum2
	implicit none

	integer					:: n,i,j,k,molno,cellnp
	integer,dimension(3)	:: ibin
	type(node), pointer		:: old, current
	real(kind(0.d0))		:: dx,dy,dz

	logical,save			:: apply_force = .false.

	!Only apply force on processor 1 in x
	!if (iblock .ne. 1) return

	!Test case focuses on a single CV
	i = 3; j = 3; k = 3

	!Find the molecules in constraint region and add to list array	   
	!Zero box averages
	box_np = 0
	allocate(list(np))
	do n = 1,np
		ibin(:) = ceiling((r(:,n)+0.5d0*domain(:))/binsize(:))+1
		if (ibin(1) .eq. i .and. ibin(2) .eq. j .and. ibin(3) .eq. k) then
			!Add molecule to overlap list
			box_np   =  box_np   + 1
			list(box_np) =  n
		endif
	enddo

	! - - - - - - - - - - Mass  - - - - - - - - - - - - 
	M = box_np

	! - - - - - - - - - Momentum  - - - - - - - - - - - 
    !Total CV flux
	MD_rhouu_dS  =		((CV2%flux(i,j,k,:,1)+CV2%flux(i,j,k,:,4)) &
	          	 		+(CV2%flux(i,j,k,:,2)+CV2%flux(i,j,k,:,5)) &
	          	 		+(CV2%flux(i,j,k,:,3)+CV2%flux(i,j,k,:,6)))/delta_t
	!Total surface stresses
	MD_Pi_dS  =	0.25d0 *((CV2%Pxy(i,j,k,:,1)-CV2%Pxy(i,j,k,:,4)) &
	              		+(CV2%Pxy(i,j,k,:,2)-CV2%Pxy(i,j,k,:,5)) &
	              		+(CV2%Pxy(i,j,k,:,3)-CV2%Pxy(i,j,k,:,6)))
	CV2%flux = 0.d0; CV2%Pxy = 0.d0


	if (M .ne. 0 .and. apply_CVforce) then
		F_constraint = MD_Pi_dS-MD_rhouu_dS + CFD_rhouu_dS-CFD_Pi_dS
	else
		F_constraint = 0.d0
	endif

	!Debugging plots of applied and remaining force
	!u_bin=0.d0; F_bin = 0.d0; F_bin2 = 0.d0
	!do i = 1, box_np
	!	n = list(i)
	!	u_bin  = u_bin  + v(:,n)
	!	F_bin  = F_bin  + a(:,n)
	!	F_bin2 = F_bin2 + a(:,n) - F_constraint(:)/(dble(M))
	!enddo

	!u_bin1 = CV%dXdt(5,5,5,:)
	!u_bin2 = CV%X(3,3,3,:)
	!if (M .ne. 0) then
	!	print'(2i4,12f10.4)', iter,M, u_bin, F_bin, F_bin2, F_constraint(:)/(dble(M))
	!else
	!	print'(2i4,12f10.5)', iter,M, u_bin, F_bin, F_bin2, (/ 0.d0, 0.d0, 0.d0 /)
	!endif

end subroutine average_over_bin

!=============================================================================
! Apply force to molecules in overlap region
!-----------------------------------------------------------------------------

subroutine apply_force
	use arrays_MD, only : r,v,a
	use physical_constants_MD, only : density
	implicit none

	integer								:: i, n
	double precision,dimension(3)		:: F_vector

	!Loop over all molecules and apply constraint
	if (M .ne. 0) F_vector = F_constraint/dble(M)
	do i = 1, box_np
		n = list(i)

		a(:,n) = a(:,n) - F_vector

		!Add external force to CV total
		call record_external_forces(F_vector,r(:,n))

        !if (any(F_constraint .ne. 0.d0)) then
			if (iter .lt. 1000) then
				write(1200+irank,'(i3,3i7,11f12.6)'),irank,iter,m_bin1,m_bin2, &
						 					  delta_t*F_vector,u_bin1,u_bin2
			endif
        !endif
	enddo

end subroutine apply_force


subroutine apply_force_tests(apply_the_force)
	use arrays_MD, only : r,v,a
	use physical_constants_MD, only : density
	implicit none

	logical, intent(in) 						:: apply_the_force

	integer										:: i, n
	double precision,dimension(3)				:: F_vector
	double precision,dimension(:,:),allocatable	:: v_temp,a_temp

	!Check evolution without constraint
	allocate(v_temp(3,np),a_temp(3,np))
	v_temp = 0.d0
	do n = 1,np
		v_temp(:,n) = v(:,n) + delta_t * a(:,n) 	!Velocity calculated from acceleration
	enddo
	m_bin2 = box_np
	u_bin2(:) = 0.d0
	F_bin2(:) = 0.d0
	do i = 1, box_np
		n = list(i)
		u_bin2(:) = u_bin2(:) + v_temp(:,n)
		F_bin2(:) = F_bin2(:) + a(:,n)
	enddo

	!Loop over all molecules and apply constraint
	a_temp(:,1:np) = a(:,1:np)
	if (M .ne. 0) F_vector = F_constraint/dble(M)
	do i = 1, box_np
		n = list(i)
		!print'(2(a,i4),2i6,9f10.5)', 'acceleration mol',i,'of',box_np, iter,n, a(:,n),a(:,n) - F_constraint/dble(M),F_constraint/dble(M)
		a_temp(:,n) = a(:,n) - F_vector

		!if (apply_the_force) then
			a(:,n) = a(:,n) - F_vector
			!Add external force to CV total
			call record_external_forces(F_vector,r(:,n))
		!endif
	enddo
	v_temp = 0.d0
	do n = 1,np
		v_temp(:,n) = v(:,n) + delta_t * a_temp(:,n) 	!Velocity calculated from acceleration
	enddo
	m_bin1 = box_np
	u_bin1(:) = 0.d0
	F_bin1(:) = 0.d0
	do i = 1, box_np
		n = list(i)
		u_bin1(:) = u_bin1(:) + v_temp(:,n)
		F_bin1(:) = F_bin1(:) + a_temp(:,n)
	enddo
	deallocate(v_temp)
	deallocate(a_temp)

	if (M .ne. 0) then
		!print'(2i4,15f8.3)', iter,m_bin1,u_bin1,u_bin2,F_bin1,F_bin2,F_constraint/dble(M)
		!if (iter .lt. 1000) then
			write(1200+irank,'(i3,3i7,9f12.6)'),irank,iter,m_bin1,m_bin2, &
					 					  delta_t*F_vector,u_bin1,u_bin2
										
		!endif
	else
		!print'(2i4,12f10.5)', iter,m_bin1,u_bin1,u_bin2,F_bin1,F_bin2, (/ 0.d0, 0.d0, 0.d0 /)
		!if (iter .lt. 1000) then
			write(1200+irank,'(i3,3i7,9f12.6)'),irank,iter,m_bin1,m_bin2, &
					 					 (/ 0.d0, 0.d0, 0.d0 /),u_bin1,u_bin2
		!endif
	endif



end subroutine apply_force_tests



end subroutine apply_CV_force



#if USE_COUPLER

!=============================================================================
! Apply force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
! Adapted serial version written by ES including cells and (originally) fully verified
!-----------------------------------------------------------------------------

subroutine apply_continuum_forces_ES(iter)
	use computational_constants_MD, only : delta_t,nh,halfdomain,ncells, & 
											cellsidelength,initialstep,Nsteps, & 
											npx,npy,npz
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	use coupler_module, only : icPmin_md,icPmax_md,jcPmin_md,jcPmax_md,kcPmin_md,kcPmax_md, & 
								iblock_realm,jblock_realm,kblock_realm
	implicit none

	integer, intent(in) 				:: iter ! iteration step, it assumes that each MD average starts from iter = 1

	integer         					:: n, molno, ixyz, cellnp
	integer								:: ii,jj,kk,icell,jcell,kcell
	integer								:: ibmin_md,ibmax_md,jbmin_md,jbmax_md,kbmin_md,kbmax_md
	integer								:: isummol
	double precision					:: inv_dtCFD
	double precision					:: isumvel, isumacc
	double precision, dimension(:,:,:),allocatable	:: u_continuum
	type(node), pointer 	        	:: old, current

	integer         					:: averagecount
	double precision					:: average

	!allocate(u_continuum(icPmin_md(iblock_realm):icPmax_md(iblock_realm), & 
	!					 jcPmin_md(jblock_realm):jcPmax_md(jblock_realm), & 
	!					 kcPmin_md(kblock_realm):kcPmax_md(kblock_realm)))

	!	print'(a,6i8)', 'limits', icPmin_md(iblock_realm),icPmax_md(iblock_realm),jcPmin_md(jblock_realm),jcPmax_md(jblock_realm),kcPmin_md(kblock_realm),kcPmax_md(kblock_realm)


	!do ii=icPmin_md(iblock_realm),icPmax_md(iblock_realm)
	!do jj=jcPmin_md(jblock_realm),jcPmax_md(jblock_realm)
	!do kk=kcPmin_md(kblock_realm),kcPmax_md(kblock_realm)

		allocate(u_continuum(1,1,1))
		u_continuum = 1.d0
		ii = 1; jj=1; kk=1

		! For each continuum cell get MD cells to average over
		!call CFD_cells_to_MD_compute_cells(ii,jj,kk,ibmin_md,ibmax_md,jbmin_md,jbmax_md,kbmin_md,kbmax_md)
		!Choose a cube in the centre of the domain
		ibmin_md=ceiling(ncells(1)/2.d0)-1; ibmax_md=ceiling(ncells(1)/2.d0)+1
		jbmin_md=ceiling(ncells(2)/2.d0)-1; jbmax_md=ceiling(ncells(2)/2.d0)+1
		kbmin_md=ceiling(ncells(3)/2.d0)-1; kbmax_md=ceiling(ncells(3)/2.d0)+1

		!Reset acceleration and velocity sums
		isummol = 0
		isumvel = 0.d0
		isumacc = 0.d0

		do icell = ibmin_md+nh, ibmax_md+nh
		do jcell = jbmin_md+nh, jbmax_md+nh
		do kcell = kbmin_md+nh, kbmax_md+nh
	 
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list

			!Calculate averages for bin
			do n = 1, cellnp    ! Loop over all particles
				molno = old%molno 	 !Number of molecule

				isumvel = isumvel + v(1,molno) 	!Add streamwise velocity to current bin
				isumacc = isumacc + a(1,molno) 	!Add acceleration to current bin
				isummol = isummol + 1
				current => old
				old => current%next 
			enddo

		enddo
		enddo
		enddo

		!Get average velocity and acceleration in bin
		if (isummol .ne. 0) then
			isumacc = isumacc/real(isummol,kind(0.d0))
		 	isumvel = isumvel/real(isummol,kind(0.d0))
		endif

		inv_dtCFD = 1/delta_t

		!Reset force averages
		average = 0.d0
		averagecount = 0

		do icell = ibmin_md+nh, ibmax_md+nh
		do jcell = jbmin_md+nh, jbmax_md+nh
		do kcell = kbmin_md+nh, kbmax_md+nh

			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list
			
			!Apply coupling force as Nie, Chen and Robbins (2004), using
			!Linear extrapolation of velocity
			do n = 1, cellnp    ! Loop over all particles
				molno = old%molno !Number of molecule

				a(1,molno)= a(1,molno) - isumacc   &
					    -(isumvel-u_continuum(ii,jj,kk))*inv_dtCFD

				current => old
				old => current%next 

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

	!enddo
	!enddo
	!enddo

end subroutine apply_continuum_forces_ES

!=============================================================================
! Apply force to match velocity to continuum using CV formulation of Nie et al 2004
! which results in matching of force and flux on a CV.
!-----------------------------------------------------------------------------

subroutine apply_continuum_forces_CV(iter)
	use computational_constants_MD, only : delta_t,nh,halfdomain,ncells, & 
											cellsidelength,initialstep,Nsteps, & 
											npx,npy,npz
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	use coupler_module, only : icPmin_md,icPmax_md,jcPmin_md,jcPmax_md,kcPmin_md,kcPmax_md, & 
							   iblock_realm,jblock_realm,kblock_realm
	use md_coupler_socket
	implicit none

	integer, intent(in) 				:: iter ! iteration step, it assumes that each MD average starts from iter = 1

	integer         					:: n, molno, ixyz, cellnp
	integer								:: ii,jj,kk,icell,jcell,kcell
	integer								:: ibmin_md,ibmax_md,jbmin_md,jbmax_md,kbmin_md,kbmax_md
	integer								:: isummol
	double precision					:: inv_dtCFD
	double precision					:: isumvel, isumacc
	double precision					:: isumflux, isumforce
	double precision, dimension(:,:,:),allocatable :: continuum_res, continuum_Fs
	type(node), pointer 	        	:: old, current

	integer         					:: averagecount
	double precision					:: average

	!allocate(u_continuum(icPmin_md(iblock_realm):icPmax_md(iblock_realm), & 
	!					 jcPmin_md(jblock_realm):jcPmax_md(jblock_realm), & 
	!					 kcPmin_md(kblock_realm):kcPmax_md(kblock_realm)))

	!	print'(a,6i8)', 'limits', icPmin_md(iblock_realm),icPmax_md(iblock_realm),jcPmin_md(jblock_realm),jcPmax_md(jblock_realm),kcPmin_md(kblock_realm),kcPmax_md(kblock_realm)


	!do ii=icPmin_md(iblock_realm),icPmax_md(iblock_realm)
	!do jj=jcPmin_md(jblock_realm),jcPmax_md(jblock_realm)
	!do kk=kcPmin_md(kblock_realm),kcPmax_md(kblock_realm)
		allocate(continuum_Fs(1,1,1))
		continuum_Fs = 1.d0
		ii = 1; jj=1; kk=1

		! For each continuum cell get MD cells to average over
		!call CFD_cells_to_MD_compute_cells(ii,jj,kk,ibmin_md,ibmax_md,jbmin_md,jbmax_md,kbmin_md,kbmax_md)
		!Choose a cube in the centre of the domain
		ibmin_md=ceiling(ncells(1)/2.d0)-1; ibmax_md=ceiling(ncells(1)/2.d0)+1
		jbmin_md=ceiling(ncells(2)/2.d0)-1; jbmax_md=ceiling(ncells(2)/2.d0)+1
		kbmin_md=ceiling(ncells(3)/2.d0)-1; kbmax_md=ceiling(ncells(3)/2.d0)+1

		!Reset acceleration and velocity sums
		isummol = 0
		isumvel = 0.d0
		isumacc = 0.d0

		do icell = ibmin_md+nh, ibmax_md+nh
		do jcell = jbmin_md+nh, jbmax_md+nh
		do kcell = kbmin_md+nh, kbmax_md+nh

			!Reset flux and force sums
			isummol   = 0

			!Calculate flux and force averages for bin
			!call compute_bin_surface_flux(icell,jcell,kcell,isumflux)
			!call compute_force_surrounding_bins(icell,jcell,kcell,isumforce)

			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list

			isumvel = 0.d0

			!Calculate velocity and molecular averages for bin
			do n = 1, cellnp    ! Loop over all particles

				molno = old%molno	!Number of molecule
				isumvel = isumvel + v(1,molno) 	!Add streamwise velocity to current bin
				isummol = isummol + 1

				current => old
				old => current%next !Use pointer in datatype to obtain next item in list
			enddo

		enddo
		enddo
		enddo

		!All constraint terms are multiplied by { (m_j \theta_j)/(\sum_i^N m_i^2 \theta_i^2) }
		!Which for unit mass reduces to { 1/(sum inside volume) }
		if (isummol .ne. 0) then
			isumflux     = isumflux     / real(isummol,kind(0.d0))
		 	isumforce    = isumforce    / real(isummol,kind(0.d0))
			continuum_Fs = continuum_res/ real(isummol,kind(0.d0))
			!if(icell .eq. 5 .and. kcell .eq. 5 .and. jcell .eq. 11 ) then 
				!print'(a,i8,4f10.5)','bin and forces',cbin,continuum_Fs
				!write(99,'(3i8,4f18.9)')iter,jcell,isummol, isumflux, isumforce, isumvel,continuum_Fs(cbin)
			!endif
		endif

		do icell = ibmin_md+nh, ibmax_md+nh
		do jcell = jbmin_md+nh, jbmax_md+nh
		do kcell = kbmin_md+nh, kbmax_md+nh

			old => cell%head(icell,jcell,kcell)%point 	!Reset old to first molecule in list
			cellnp = cell%cellnp(icell,jcell,kcell)

			!Apply coupling force for CV form of Nie, Chen and Robbins (2004)
			do n = 1, cellnp    ! Loop over all particles
				molno = old%molno !Number of molecule

				a(1,molno) = a(1,molno) + ((isumflux - isumforce) + continuum_Fs(icell,jcell,kcell))

				current => old
				old => current%next 

			enddo
	 	
		enddo
		enddo
		enddo

		nullify(current)        !Nullify current as no longer required
		nullify(old)            !Nullify old as no longer required

end subroutine apply_continuum_forces_CV

!------------------------------------------
!subroutine CFD_cells_to_MD_compute_cells(ibmin_cfd,icmin_cfd,jbmin_cfd,jcmin_cfd,kbmin_cfd,kcmin_cfd, & 
!										  ibmin_md, icmin_md, jbmin_md, jcmin_md, kbmin_md, kcmin_md)
!	implicit none

!	integer,intent(in)		:: ibmin_cfd,icmin_cfd,jbmin_cfd,jcmin_cfd,kbmin_cfd,kcmin_cfd
!	integer,intent(out)		:: ibmin_md,icmin_md,jbmin_md,jcmin_md,kbmin_md,kcmin_md

subroutine CFD_cells_to_MD_compute_cells(ii_cfd,jj_cfd,kk_cfd, & 
										  ibmin_md, ibmax_md, jbmin_md, jbmax_md, kbmin_md, kbmax_md)
	use coupler_module, only : xg, yg, zg, xL_md,yL_md,zL_md, iblock_realm,jblock_realm,kblock_realm
	use computational_constants_MD, only : cellsidelength
	implicit none

	integer,intent(in)		:: ii_cfd,jj_cfd,kk_cfd
	integer,intent(out)		:: ibmin_md,ibmax_md,jbmin_md,jbmax_md,kbmin_md,kbmax_md
	
	double precision		:: xL_min,xL_max,yL_min,yL_max,zL_min,zL_max

	! Get minimum point in processors domain
	xL_min = xL_md*(iblock_realm-1); xL_max = xL_md*(iblock_realm)
	yL_min = yL_md*(jblock_realm-1); yL_max = yL_md*(jblock_realm)
	zL_min = zL_md*(kblock_realm-1); zL_max = zL_md*(kblock_realm)

	! Get range of cells to check so that top and bottom of current CFD cell are covered
	ibmin_md = (xg(ii_cfd  ,jj_cfd  )-xL_min)/cellsidelength(1)+1
	ibmax_md = (xg(ii_cfd+1,jj_cfd  )-xL_min)/cellsidelength(1)+1
	jbmin_md = (yg(ii_cfd  ,jj_cfd  )-yL_min)/cellsidelength(2)+1
	jbmax_md = (yg(ii_cfd  ,jj_cfd+1)-yL_min)/cellsidelength(2)+1
	kbmin_md = (zg(     kk_cfd      )-zL_min)/cellsidelength(3)+1
	kbmax_md = (zg(     kk_cfd+1    )-zL_min)/cellsidelength(3)+1

	print'(a,9i8)','indices', ii_cfd,ibmin_md,ibmax_md,jj_cfd,jbmin_md,jbmax_md,kk_cfd,kbmin_md,kbmax_md
	print*,'xcells', xg(ii_cfd  ,jj_cfd  ),(xg(ii_cfd  ,jj_cfd  )-xL_min)/cellsidelength(1)+1, (xg(ii_cfd+1,jj_cfd  )-xL_min)/cellsidelength(1)+1
	print*,'ycells', yg(ii_cfd  ,jj_cfd  ),(yg(ii_cfd  ,jj_cfd  )-yL_min)/cellsidelength(2)+1, (yg(ii_cfd+1,jj_cfd  )-yL_min)/cellsidelength(2)+1
	print*,'zcells', zg(kk_cfd  ),(zg(kk_cfd)-zL_min)/cellsidelength(3)+1, (zg(kk_cfd+1)-zL_min)/cellsidelength(3)+1

end subroutine CFD_cells_to_MD_compute_cells

#endif



subroutine pointsphere(centre,targetradius,start_iter)
	use arrays_MD, only : r,a
	use physical_constants_MD, only : np, pi
	use computational_constants_MD, only : iter
	implicit none

	integer,optional,intent(in)					:: start_iter
	double precision, intent(in)				:: targetradius
	double precision, dimension(3), intent(in)	:: centre

	integer							:: n, start
	double precision 				:: radius, radius2, Fapplied, magnitude
	double precision 				:: rspherical,rspherical2
	double precision, dimension(3)	:: rmapped

	!Grow pore slowly
	if (present(start_iter)) then
		start = start_iter
	else
		start = 1
	endif
	radius = min((iter-start)/1000.d0,targetradius)
	magnitude = 100.d0

	!Define square of radius
	radius2 = radius**2

	! Loop through all molecules
	do n=1,np
		rmapped = r(:,n)-centre						!Map to origin of sphere
		rspherical2 = dot_product(rmapped,rmapped)	!Get position in spherical coordinates
		rspherical = sqrt(rspherical2)
    	!theta 	   = acos(rmapped(3)/rspherical)
    	!phi 	   = atan(rmapped(2)/rmapped(1))
		if (rspherical .lt. radius2) then			!Check if inside sphere
			!if (rspherical .lt. 1.d0) rspherical = 1.d0			!Bounded force
			Fapplied = magnitude *(1.d0 + 1.d0/exp(rspherical)**2.d0)			!Apply force which diverges at centre
			!print'(i8,13f10.5)', n,rmapped,rmapped(3)/rspherical,rspherical,theta,phi, Fapplied*sin(theta)*cos(phi), Fapplied*sin(theta)*sin(phi),Fapplied*cos(theta),a(:,n)
    		a(1,n)= a(1,n) + Fapplied * rmapped(1)!*sin(theta)*cos(phi)
			a(2,n)= a(2,n) + Fapplied * rmapped(2)!*sin(theta)*sin(phi)
    		a(3,n)= a(3,n) + Fapplied * rmapped(3)!*cos(theta)
		endif
	enddo

end subroutine pointsphere
