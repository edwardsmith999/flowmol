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

	use module_record_external_forces
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
	use boundary_MD, only: bforce_flag, bforce_dxyz, &
	                                      bforce_off, bforce_NCER, bforce_OT, &
	                                      bforce_Flekkoy
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

	integer :: n,ixyz,flag,block(3),nPxy(3)
	real(kind(0.d0)) :: xyz,thresh,hdom
	real(kind(0.d0)), dimension(3) :: tops,bottoms

	tops    = (/domain(1)/2.d0 - dists(2), &
				domain(2)/2.d0 - dists(4), &
				domain(3)/2.d0 - dists(6)  /)

	bottoms = (/dists(1) - domain(1)/2.d0, &
				dists(3) - domain(2)/2.d0, &
				dists(5) - domain(3)/2.d0  /)

	block  = (/iblock,jblock,kblock/)
	nPxy  = (/npx,npy,npz/)

	do n = 1, np
	do ixyz = 1, 3
			
		if ( r(ixyz,n)   .gt. tops(ixyz)  .and. &
		     block(ixyz) .eq. nPxy(ixyz)       ) then
			
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
    
		call apply_bforce(a(ixyz,n),ixyz,xyz,thresh,hdom,flag)

	end do
	end do

end subroutine simulation_apply_boundary_force

subroutine apply_bforce(a_in,ixyz,xyz,thresh,hdom,flag)
	use boundary_MD, only: bforce_NCER, bforce_off, bforce_pdf_input, &
                           bforce_pdf_nsubcells
	use calculated_properties_MD,   only: pressure
    use computational_constants_MD, only: cellsidelength 
	use interfaces,                 only: error_abort
	implicit none

	real(kind(0.d0)), intent(inout) :: a_in	  !Accel to which add bforce
	real(kind(0.d0)), intent(in) :: xyz       !Position, 1D
	real(kind(0.d0)), intent(in) :: thresh    !Threshold for bforce, 1D
	real(kind(0.d0)), intent(in) :: hdom      !Domain edge
	integer, intent(in) :: ixyz               !Direction 
	integer, intent(in) :: flag               !Type of bforce 

	real(kind(0.d0)) :: numer,denom,ratio,P,f,dxyz
	character(128)   :: string
    integer :: subcell

	select case ( flag )
	case ( bforce_off  )

		return

	case ( bforce_NCER )

		P     = max(pressure,1.d0)                !Account for negative pressure
		numer = xyz - thresh                      !+ve or -ve dependent on pos in domain
		denom = 1.d0 - (xyz-thresh)/(hdom-thresh) !denom always +ve	
		ratio = numer / denom	

		a_in  = a_in - ratio*P	

    case ( bforce_pdf_input )

        dxyz = abs(hdom - xyz)
        subcell = ceiling(dxyz/cellsidelength(ixyz))*bforce_pdf_nsubcells
        f = pull_from_bforce_pdf(subcell,ixyz)
        !print*, 'subcell, ixyz, f', subcell, ixyz, f
        a_in = a_in + f

	case default

		string="MD uncoupled boundary force only developed for NCER case"
		call error_abort(string)

	end select

contains

    function pull_from_bforce_pdf(subcell, component) result(F)
        use boundary_MD, only: bforce_pdf_min, bforce_pdf_max, &
                               bforce_pdf_nbins, bforce_pdf_input_data, &
                               bforce_pdf_binsize
        implicit none
      
        integer, intent(in) :: subcell, component

        real(kind(0.d0)) :: randu, F, PF, maxP
        integer :: bin , n 
        logical :: success
        
        real(kind(0.d0)) :: s, p, nmax
        integer :: maxattempts

        ! Choose maxattempts so that 99% of pull calls will at
        ! least hit the 0 bin if every other bin has probability 0
        ! This is done assuming a geometric distribution:
        !
        !   P(Y=k) = p(1-p)^k, 
        ! 
        ! where Y is the number of trials (random numbers generated)
        ! before success (is in the 0 bin, where prob of force there
        ! is guaranteed to be non-zero). p here is 1/histbins (uniform
        ! sampling). We seek nmax from
        !
        !   s = 0.99 = \sum_{k=0}^{n-1} p(1-p)^k 
        !            = p(1-(1-p)^n)/(1-(1-p))
        !            = 1 - (1-p)^n

        s = 0.99d0
        p = 1.d0/real(bforce_pdf_nbins,kind(0.d0))
        nmax = log(1.d0 - s)/log(1.d0 - p)
        maxattempts = ceiling(nmax)

        ! Maximum value of the PDF
        maxP = maxval(bforce_pdf_input_data(subcell, :, component))

        success = .false.
        do n=1,maxattempts

            call random_number(randu)
            F = bforce_pdf_min + randu*(bforce_pdf_max-bforce_pdf_min)
            bin = ceiling((F-bforce_pdf_min)/bforce_pdf_binsize)

            call random_number(randu) 
            PF = randu * maxP

            if ( PF .le. bforce_pdf_input_data(subcell, bin, component)) then
                success = .true.
                exit
            endif

        end do
       
        if (success .eqv. .false.) then
            print*, 'Failed to pull bforce, applying 0.d0!'
            F = 0.d0
        end if 

    end function pull_from_bforce_pdf 

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
			g = approx_RDF(r(2,n),halfdomain(2)-dy,halfdomain(2))
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
		g = approx_RDF(r(2,n),halfdomain(2)-dy,halfdomain(2))

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

!Obtain an approximation for the Radial Distribution Function using the
!expression from Appendix A of T. Werder et al. J Comp. Phys 205 (2005) 373–390
! NOTE - this assumes a density of 0.6 and temp of 1.7
! also, the units are in nano-meters so a conversion is required
function approx_RDF(r,rmin,rmax) result(g)
	use librarymod, only : heaviside 

	real(kind=kind(0.d0)), intent(in)	:: r,rmin,rmax
	real(kind=kind(0.d0))				:: g, rhat, L
	real(kind=kind(0.d0))				:: LJunits2nm

	!Define local coordinate as const runs from 0 < y < L/2
	L = rmax - rmin
	rhat = (r - rmin - 0.5*L)

    !Sanity Check and exceptions
    if (rhat .lt. 0.d0) then
        g = 0
        return
    elseif (rhat .gt. 0.5*L) then
		stop " Werder et al RDF error - input y cannot be greater than rmax"
    endif

	!Conversion factors
	LJunits2nm = 2.94d0
	rhat = rhat/LJunits2nm

	g = ( 10.8007 + 0.860717*rhat    - 172.468* rhat**2 					& 
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

	g = -g

	if (g .ne. 0.d0) print'(5f10.5,5i5)',  rhat, g, rmin ,r, rmax, & 
							(heaviside(rhat-0.2975)-heaviside(rhat)		)   , & 
							(heaviside(rhat-0.3475)-heaviside(rhat-0.2975)) , &
							(heaviside(rhat-0.3475)-heaviside(rhat-0.2975)) , &
							(heaviside(rhat-0.3975)-heaviside(rhat-0.3475)) , &
							(heaviside(rhat-1.0000)-heaviside(rhat-0.3975))

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



!==================================================================
!  			A P P L Y      C V      C O N S T R A I N T
!==================================================================

module apply_CV_force_mod
	use calculated_properties_MD, only : nbins
	use computational_constants_MD, only: domain, & 
                                          CVforce_starttime, &
                                          CVforce_flag
	use physical_constants_MD, only: np
	implicit none

	real(kind(0.d0))				::	volume
	real(kind(0.d0)),dimension(3)	::	Fbinsize

    logical                         :: couette_timeevolve=.false.

contains 

! ---------------------------------------------------------
! Take the difference between top and bottom stresses for
! a range of CV in 3 dimensions to obtain force
!   Note bin limits here should be local to processor

subroutine flux2force(flux,bl,force)
	implicit none

	integer,intent(in)             		 :: bl(6)
	real(kind(0.d0)),intent(in), & 
		allocatable,dimension(:,:,:,:,:) :: flux
	real(kind(0.d0)),intent(inout), & 
		allocatable,dimension(:,:,:,:) 	 :: force

	force(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:) & 
		  =  ( flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,1)  & 
			  -flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,4)) &
		    +( flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,2)  & 
			  -flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,5)) &
		    +( flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,3)  & 
			  -flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,6))

end subroutine

! ---------------------------------------------------------

subroutine get_MD_stresses(MD_stress)
	use CV_objects, only : 	CV_constraint
	implicit none

	real(kind(0.d0)),intent(inout), & 
		allocatable,dimension(:,:,:,:,:) 	:: MD_stress

    integer :: ibin,jbin,kbin,ixyz

	MD_stress(:,:,:,:,:) = 0.25d0 * CV_constraint%Pxy(:,:,:,:,:)
    CV_constraint%Pxy = 0.d0

end subroutine get_MD_stresses

! ---------------------------------------------------------

subroutine get_MD_fluxes(MD_flux, r_in, v_in)
	use CV_objects, only : 	CV_constraint
	use arrays_MD, only : r, v, a
	use computational_constants_MD, only : delta_t
	use calculated_properties_MD, only : nbinso
	use cumulative_momentum_flux_mod, only : cumulative_momentum_flux
	implicit none

	real(kind(0.d0)),intent(inout), & 
		allocatable,dimension(:,:,:,:,:) 	:: MD_flux
	real(kind(0.d0)),intent(in),optional, &
		allocatable, dimension(:,:)			:: r_in, v_in

	integer	:: n
	real(kind(0.d0)),dimension(:,:),allocatable	:: r_temp,v_temp

	allocate(r_temp(3,np),v_temp(3,np))
	if (present(r_in) .and. present(v_in)) then
		if (.not. allocated(r_in) .or. .not. allocated(v_in)) then
			stop "Error in get_MD_fluxes -- r_in or v_in not allocated"
		endif
		r_temp = r_in
		v_temp = v_in
	else
		!Get velocity/positions at next timestep without constraint
		do n = 1,np
			v_temp(:,n) = v(:,n) + delta_t * a(:,n) 		!Velocity calculated from acceleration
			r_temp(:,n) = r(:,n) + delta_t * v_temp(:,n)	!Position calculated from velocity
		enddo
	endif
	!Get momentum flux using velocity and project backwards
	CV_constraint%flux = 0.d0
	call cumulative_momentum_flux(r_temp,v_temp,CV_constraint%flux)

	! Exchange all halo values for CV values ready to apply forces
	call CV_constraint%swap_halos(nbinso)
	MD_flux(:,:,:,:,:) = CV_constraint%flux(:,:,:,:,:)/delta_t

end subroutine get_MD_fluxes

! ---------------------------------------------------------

subroutine get_boxnp(boxnp)
	use arrays_MD, only : r
	implicit none

	integer,intent(out), & 
		allocatable,dimension(:,:,:):: boxnp

	integer							:: n
	integer,dimension(3)			:: bin

	allocate(boxnp(nbins(1)+2,nbins(2)+2,nbins(3)+2)); boxnp = 0
	do n = 1,np
		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1
		!Add molecule to overlap list
		boxnp( bin(1),bin(2),bin(3)) = boxnp( bin(1),bin(2),bin(3)) + 1
	enddo


end subroutine get_boxnp

! ---------------------------------------------------------


subroutine get_CFD_velocity(binlimits, & 
							u_CFD)
    use messenger, only : localise_bin
    use computational_constants_MD, only : irank, iter, globaldomain, delta_t
    use physical_constants_MD, only : pi, tethereddisttop,tethereddistbottom, density
    use calculated_properties_MD, only : gnbins
    use librarymod, only : couette_analytical_fn, linspace
	implicit none

	integer,dimension(6),intent(in)		:: binlimits(6)
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:,:,:) 	:: u_CFD

	integer,dimension(6)				 	:: lbl
	integer									:: appliedbins, jb, j
	integer									:: wallbintop,wallbinbottom
	real(kind(0.d0))             			:: sin_mag, sin_period
	real(kind(0.d0))             			:: t,Re,Uwall,H
	real(kind(0.d0)), & 
		allocatable,dimension(:) 			:: Utemp
	allocate(u_CFD(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); u_CFD=0.d0

    lbl(1:5:2) = max(localise_bin((/binlimits(1), & 
									binlimits(3), & 
									binlimits(5)/)),(/2,2,2 /))
    lbl(2:6:2) = min(localise_bin((/binlimits(2), & 
									binlimits(4), & 
									binlimits(6)/)),nbins(:)+1)

#if USE_COUPLER
	! v v v v v   This should be replaced by coupler call   v v v v v 
    !call socket_get_velocity(u_CFD)
	! ^ ^ ^ ^ ^   This should be replaced by coupler call   ^ ^ ^ ^ ^
#else
    select case (CVforce_flag)
    case(0)
        !Do nothing, zero CFD velocity
    case(1)
        ! - - - - -  Sinusoid of velocity corresponding to CFD stresses - - - - - - - - - 
        !The correction parameter is obtained from curve fitting based on the 
        !resulting sinusoid after application of the surface flux based force. I have
        !no idea why there is a difference of 1.07205...
        sin_mag = 2.0d0;  sin_period = 100.d0
        u_CFD(lbl(1):lbl(2),lbl(3):lbl(4),lbl(5):lbl(6),1) = & 
                    1.07205d0*sin_mag*(sin_period/(2.d0*pi)) & 
                    *sin(2.d0*pi*((iter-CVforce_starttime-0.5)/sin_period))
    case(4)
        sin_mag = 2.0d0;  sin_period = 100.d0
        u_CFD(lbl(1):lbl(2),lbl(3):lbl(4),lbl(5):lbl(6),1) = & 
                    1.07205d0*sin_mag*(sin_period/(2.d0*pi)) & 
                    *(sin(2.d0*pi*((iter-CVforce_starttime-0.5)/sin_period))-0.5d0)
        

    case(5)
        !Number of bins -- use global liquid domain
        wallbintop    = ceiling(tethereddisttop(2)   /Fbinsize(2))
        wallbinbottom = ceiling(tethereddistbottom(2)/Fbinsize(2))
	    appliedbins = gnbins(2) - wallbintop - wallbinbottom !lbl(4)-lbl(3)+1
	    allocate(Utemp(appliedbins))	

	    if (couette_timeevolve) then
	        !Time evolving solution
            Re = 1.4d0
            Uwall = 1.d0
            H = globaldomain(2) - tethereddistbottom(2) - tethereddisttop(2)
            t = (iter-CVforce_starttime)*delta_t
            Utemp = couette_analytical_fn(t,Re,Uwall,H,appliedbins,0)
	    else
		    ! Steady state solution
		    Utemp = linspace(0.d0, Uwall, appliedbins)
	    endif



	    !Copy u to CFD velocity
	    do jb = lbl(3),lbl(4)
            j = jb+wallbinbottom
		    u_CFD(lbl(1):lbl(2),jb,lbl(5):lbl(6),1) = Utemp(j)
            if (mod(iter,100) .eq. 0) then
    		    write(8000+iter/100,'(a,2i8,2f18.12)'), 'uvel', jb,j, Utemp(j),u_CFD(3,jb,3,1)
            endif
	    enddo
    end select
#endif

end subroutine get_CFD_velocity

subroutine get_CFD_stresses_fluxes(binlimits, & 
								   CFD_stress, & 
								   CFD_flux)
    use messenger, only : localise_bin
    use computational_constants_MD, only : irank, iter, globaldomain,delta_t
    use physical_constants_MD, only : pi, tethereddisttop,tethereddistbottom, density
    use calculated_properties_MD, only : gnbins
    use librarymod, only : couette_analytical_stress_fn
	implicit none

	integer,dimension(6),intent(in)			:: binlimits(6)
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:,:,:,:) 	:: CFD_stress,CFD_flux


	integer									:: appliedbins, taumin,taumax,ib,jb,kb,j
	integer									:: wallbintop,wallbinbottom
	integer,dimension(6)				 	:: lbl
	real(kind(0.d0))             			:: sin_mag, sin_period
	real(kind(0.d0))             			:: Re,Uwall,H,t
	real(kind(0.d0)), & 
		allocatable,dimension(:) 			:: tautemp

	allocate(CFD_stress(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6)); CFD_stress=0.d0
	allocate(  CFD_flux(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6)); CFD_flux=0.d0

    lbl(1:5:2) = max(localise_bin((/binlimits(1), & 
									binlimits(3), & 
									binlimits(5)/)),(/2,2,2 /))
    lbl(2:6:2) = min(localise_bin((/binlimits(2), & 
									binlimits(4), & 
									binlimits(6)/)),nbins(:)+1)

#if USE_COUPLER
	! v v v v v   This should be replaced by coupler call   v v v v v 
    !call socket_get_fluxes_and_stresses(CFD_stress,CFD_flux)
	! ^ ^ ^ ^ ^   This should be replaced by coupler call   ^ ^ ^ ^ ^
#else
    select case (CVforce_flag)
    case(0)
        !Do nothing, zero CFD constraint
    case(1)
        !Sinusoidally fluctuating stress
        CFD_stress(lbl(1):lbl(2),lbl(3):lbl(4),lbl(5):lbl(6),:,:) = 0.d0
        !Sin function shifted by pi/2 so velocity which is given
        !by the integral (cos) alternates around zero
        sin_mag = 2.0d0;  sin_period = 100.d0
        CFD_stress(lbl(1):lbl(2), & 
                   lbl(3):lbl(4), & 
                   lbl(5):lbl(6),1,2) =  0.5d0*sin_mag * cos(2.d0*pi*((iter-CVforce_starttime)/sin_period))
        CFD_stress(lbl(1):lbl(2), & 
                   lbl(3):lbl(4), & 
                   lbl(5):lbl(6),1,5) = -0.5d0*sin_mag * cos(2.d0*pi*((iter-CVforce_starttime)/sin_period))

    case(2)

        !SINGLE CELL VORTEX
        if (lbl(1) .ne. lbl(2) .or. lbl(3) .ne. lbl(4) .or. lbl(5) .ne. lbl(6)) then
            stop "Error in get_CFD_stresses_fluxes -- for CVforce_flag 3, only 1 CV allowed in x and y"
        endif
	    !Top
        CFD_stress(lbl(1),lbl(3),lbl(5),:,1) = (/  0.d0,  1.d0,  0.d0 /)
        CFD_stress(lbl(1),lbl(3),lbl(5),:,2) = (/ -1.d0,  0.d0,  0.d0 /)
        CFD_stress(lbl(1),lbl(3),lbl(5),:,3) = (/  0.d0,  0.d0,  0.d0 /)
    	!Bottom
        CFD_stress(lbl(1),lbl(3),lbl(5),:,4) = (/  0.d0,  1.d0,  0.d0 /)
        CFD_stress(lbl(1),lbl(3),lbl(5),:,5) = (/ -1.d0,  0.d0,  0.d0 /)
        CFD_stress(lbl(1),lbl(3),lbl(5),:,6) = (/  0.d0,  0.d0,  0.d0 /)
        !Adjacent surfaces
        CFD_stress(lbl(1)-1,lbl(3),lbl(5),:,1) = -CFD_stress(lbl(1),lbl(3),lbl(5),:,4)
        CFD_stress(lbl(1)+1,lbl(3),lbl(5),:,4) = -CFD_stress(lbl(1),lbl(3),lbl(5),:,1)
        CFD_stress(lbl(1),lbl(3)-1,lbl(5),:,2) = -CFD_stress(lbl(1),lbl(3),lbl(5),:,5)
        CFD_stress(lbl(1),lbl(3)+1,lbl(5),:,5) = -CFD_stress(lbl(1),lbl(3),lbl(5),:,2)
        CFD_stress(lbl(1),lbl(3),lbl(5)-1,:,3) = -CFD_stress(lbl(1),lbl(3),lbl(5),:,6)
        CFD_stress(lbl(1),lbl(3),lbl(5)+1,:,6) = -CFD_stress(lbl(1),lbl(3),lbl(5),:,3)

    case(3)

        !CELL VORTEX AND STRETCHED ELONGATION in all z
        if (lbl(1) .ne. lbl(2) .or. lbl(3) .ne. lbl(4)) then
            stop "Error in get_CFD_stresses_fluxes -- for CVforce_flag 3, only 1 CV allowed in x and y"
        endif
        do kb = lbl(5),lbl(6)
	        !Top
            CFD_stress(lbl(1),lbl(3),kb,:,1) = (/  0.d0,  1.d0,  0.d0 /)
            CFD_stress(lbl(1),lbl(3),kb,:,2) = (/  1.d0,  0.d0,  0.d0 /)
            CFD_stress(lbl(1),lbl(3),kb,:,3) = (/  0.d0,  0.d0,  0.d0 /)
	        !Bottom
            CFD_stress(lbl(1),lbl(3),kb,:,4) = (/  0.d0,  1.d0,  0.d0 /)
            CFD_stress(lbl(1),lbl(3),kb,:,5) = (/  1.0d0, 0.d0,  0.d0 /)
            CFD_stress(lbl(1),lbl(3),kb,:,6) = (/  0.d0,  0.d0,  0.d0 /)

	        !Top
            CFD_stress(lbl(1)+1,lbl(3),kb,:,1) = (/  0.d0, -1.0d0, 0.d0 /)
            CFD_stress(lbl(1)+1,lbl(3),kb,:,2) = (/  1.d0,  0.d0,  0.d0 /)
            CFD_stress(lbl(1)+1,lbl(3),kb,:,3) = (/  0.d0,  0.d0,  0.d0 /)
	        !Bottom
            CFD_stress(lbl(1)+1,lbl(3),kb,:,4) = (/  0.d0, -1.0d0, 0.d0 /)
            CFD_stress(lbl(1)+1,lbl(3),kb,:,5) = (/  1.0d0, 0.d0,  0.d0 /)
            CFD_stress(lbl(1)+1,lbl(3),kb,:,6) = (/  0.d0,  0.d0,  0.d0 /)

            CFD_stress(lbl(1),lbl(3),kb-1,:,3) = -CFD_stress(lbl(1),lbl(3),kb,:,6)
            CFD_stress(lbl(1),lbl(3),kb+1,:,6) = -CFD_stress(lbl(1),lbl(3),kb,:,3)

            CFD_stress(lbl(1)+1,lbl(3),kb-1,:,3) = -CFD_stress(lbl(1)+1,lbl(3),kb,:,6)
            CFD_stress(lbl(1)+1,lbl(3),kb+1,:,6) = -CFD_stress(lbl(1)+1,lbl(3),kb,:,3)
        enddo

    case(4)

        ! VORTEX GENERATOR 
        do kb = lbl(5),lbl(6)
	        !Top
            CFD_stress(lbl(1)+1,lbl(3),kb,:,1) = (/  0.d0, -1.0d0, 0.d0 /)
            CFD_stress(lbl(1)+1,lbl(3),kb,:,2) = (/  1.d0,  0.d0,  0.d0 /)
            CFD_stress(lbl(1)+1,lbl(3),kb,:,3) = (/  0.d0,  0.d0,  0.d0 /)
	        !Bottom
            CFD_stress(lbl(1)+1,lbl(3),kb,:,4) = (/  0.d0, -1.0d0, 0.d0 /)
            CFD_stress(lbl(1)+1,lbl(3),kb,:,5) = (/  1.0d0, 0.d0,  0.d0 /)
            CFD_stress(lbl(1)+1,lbl(3),kb,:,6) = (/  0.d0,  0.d0,  0.d0 /)

            CFD_stress(lbl(1),lbl(3),kb-1,:,3) = -CFD_stress(lbl(1),lbl(3),kb,:,6)
            CFD_stress(lbl(1),lbl(3),kb+1,:,6) = -CFD_stress(lbl(1),lbl(3),kb,:,3)

            CFD_stress(lbl(1)+1,lbl(3),kb-1,:,3) = -CFD_stress(lbl(1)+1,lbl(3),kb,:,6)
            CFD_stress(lbl(1)+1,lbl(3),kb+1,:,6) = -CFD_stress(lbl(1)+1,lbl(3),kb,:,3)

            !Sin function shifted by pi/2 so velocity which is given
            !by the integral (cos) alternates around zero
            sin_mag = 2.0d0;  sin_period = 100.d0
            CFD_stress(lbl(1):lbl(2), & 
                       lbl(3):lbl(4), & 
                       kb,1,4) = sin_mag * cos(2.d0*pi*((iter-CVforce_starttime)/sin_period))
        enddo

    
    case(5)
        !Number of bins -- use global liquid domain
        wallbintop    = ceiling(tethereddisttop(2)   /Fbinsize(2))
        wallbinbottom = ceiling(tethereddistbottom(2)/Fbinsize(2))
	    appliedbins = gnbins(2) - wallbintop - wallbinbottom !lbl(4)-lbl(3)+1
        !Add one extra record to get values at surfaces
        appliedbins = appliedbins + 1
	    allocate(tautemp(appliedbins))	

	    ! Couette Analytical solution
	    if (couette_timeevolve) then
		    !Time evolving solution
            Re = 1.4d0
            Uwall = 1.d0
            H = globaldomain(2) - tethereddistbottom(2) - tethereddisttop(2)
            t = (iter-CVforce_starttime)*delta_t
		    tautemp = couette_analytical_stress_fn(t,Re,Uwall,H,appliedbins,0)
            tautemp = tautemp/density
	    else
		    ! Steady state solution
		    tautemp = Uwall/domain(2)
	    endif

	    !Copy tau to CFD surfaces
	    do jb = lbl(3),lbl(4)
            j = jb-1
		    CFD_stress(lbl(1):lbl(2),jb,lbl(5):lbl(6),1,2) = tautemp(j)	    !Top
		    CFD_stress(lbl(1):lbl(2),jb,lbl(5):lbl(6),1,5) = tautemp(j-1)	!Bottom
            if (mod(iter,100) .eq. 0) then
    		    write(6000+iter/100,'(a,2i8,4f18.12)'), 'Stress', jb,j, tautemp(j),tautemp(j-1),CFD_stress(3,jb,3,1,2),CFD_stress(3,jb,3,1,5)
            endif
	    enddo

    case default
        stop "Error in get_CFD_stresses_fluxes -- CVforce_flag is unknown case"
    end select

	CFD_stress = CFD_stress*volume
	CFD_flux   = CFD_flux*volume

#endif

end subroutine get_CFD_stresses_fluxes

! ---------------------------------------------------------

subroutine get_MD_stresses_fluxes(binlimits, & 
								  MD_stress, & 
								  MD_flux)

	!use , only : nbins
	implicit none

	integer,intent(in)             			:: binlimits(6)
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:,:,:,:) 	:: MD_stress,MD_flux

	allocate( MD_stress(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6)); MD_stress=0.d0
	allocate(   MD_flux(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6)); MD_flux=0.d0

    !Call get fluxes first which exchanges halos
	call get_MD_fluxes(MD_flux)
	call get_MD_stresses(MD_stress)

end subroutine get_MD_stresses_fluxes


! ---------------------------------------------------------

subroutine get_Fcfdflux_CV(CFD_flux,  & 
						   binlimits, &
						   Fcfdflux_CV)
    use messenger, only : localise_bin
	implicit none

	integer,intent(in)             		 :: binlimits(6)
	real(kind(0.d0)),intent(in), & 
		allocatable,dimension(:,:,:,:,:) :: CFD_flux
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:,:,:) 	 :: Fcfdflux_CV

	integer,dimension(6)				 :: localbinlimits

    localbinlimits(1:5:2) = max(localise_bin((/binlimits(1), & 
											   binlimits(3), & 
											   binlimits(5)/)),(/2,2,2 /))
    localbinlimits(2:6:2) = min(localise_bin((/binlimits(2), & 
											   binlimits(4), & 
											   binlimits(6)/)),nbins(:)+1)

	allocate( Fcfdflux_CV(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); Fcfdflux_CV=0.d0
	call flux2force(CFD_flux, localbinlimits, Fcfdflux_CV)

end subroutine get_Fcfdflux_CV

! ---------------------------------------------------------

subroutine get_Fstresses_CV(CFD_stress, & 
							MD_stress,  &
							binlimits,  &
							Fstresses_CV)
    use messenger, only : localise_bin
	implicit none

	integer,intent(in)             		 :: binlimits(6)
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:,:,:)   :: Fstresses_CV
	real(kind(0.d0)),intent(in), & 
		allocatable,dimension(:,:,:,:,:) :: CFD_stress,MD_stress

	integer,dimension(6)				:: localbinlimits
	real(kind(0.d0)),allocatable,dimension(:,:,:,:,:) :: Pxy

    localbinlimits(1:5:2) = max(localise_bin((/binlimits(1), & 
											   binlimits(3), & 
											   binlimits(5)/)),(/2,2,2 /))
    localbinlimits(2:6:2) = min(localise_bin((/binlimits(2), & 
											   binlimits(4), & 
											   binlimits(6)/)),nbins(:)+1)

	allocate(Fstresses_CV(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); Fstresses_CV=0.d0
	allocate(Pxy( nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6)); Pxy = 0.d0

	Pxy = MD_stress-CFD_stress
	call flux2force(Pxy, localbinlimits, Fstresses_CV)

end subroutine get_Fstresses_CV

! ---------------------------------------------------------

subroutine get_Fstresses_mol(CFD_stress,   & 
							 MD_stress,    &
						  	 boxnp,	 	   &
						  	 flag,	 	   &
							 binlimits,    & 
							 Fstresses_mol)
	use arrays_MD, only: r
	use librarymod, only : linearsurface_weight, lagrange_poly_weight
	use computational_constants_MD, only : iter
    use calculated_properties_MD, only : nbins
	implicit none

	integer,intent(in)             		 :: binlimits(6), flag
	integer,intent(in), & 
		allocatable,dimension(:,:,:)	 :: boxnp
	real(kind(0.d0)),intent(in), & 
		allocatable,dimension(:,:,:,:,:) :: CFD_stress,MD_stress
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:)   	 :: Fstresses_mol

	integer								 :: n, bin(3), ixyz, s, ibin,jbin,kbin
	real(kind(0.d0)),allocatable,& 
		dimension(:,:,:,:,:) 			 :: Pxy
	real(kind(0.d0)), & 
		allocatable,dimension(:,:,:,:)   :: Fstresses_CV

	allocate(Fstresses_mol(3,np));  Fstresses_mol=0.d0
	allocate(Fstresses_CV( nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); Fstresses_CV=0.d0
	allocate(Pxy( nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6)); Pxy = 0.d0

	!Get F_stresses on the CV
	call get_Fstresses_CV(CFD_stress, MD_stress, binlimits, Fstresses_CV)

    !Get stress to use in weighting, either 
    if (flag .eq. 1) then
    	Pxy = -CFD_stress
    elseif (flag .eq. 2) then
    	Pxy = MD_stress-CFD_stress
    endif

    !Add sign convention to take into account surface normals and divide by boxnp
    do ibin = 2,nbins(1)+1
    do jbin = 2,nbins(2)+1
    do kbin = 2,nbins(3)+1
        do ixyz = 1,3
	        Pxy(ibin,jbin,kbin,ixyz,1) =  Pxy(ibin,jbin,kbin,ixyz,1)/dble(boxnp(ibin,jbin,kbin))
	        Pxy(ibin,jbin,kbin,ixyz,2) =  Pxy(ibin,jbin,kbin,ixyz,2)/dble(boxnp(ibin,jbin,kbin))
	        Pxy(ibin,jbin,kbin,ixyz,3) =  Pxy(ibin,jbin,kbin,ixyz,3)/dble(boxnp(ibin,jbin,kbin))
	        Pxy(ibin,jbin,kbin,ixyz,4) = -Pxy(ibin,jbin,kbin,ixyz,4)/dble(boxnp(ibin,jbin,kbin))
	        Pxy(ibin,jbin,kbin,ixyz,5) = -Pxy(ibin,jbin,kbin,ixyz,5)/dble(boxnp(ibin,jbin,kbin))
	        Pxy(ibin,jbin,kbin,ixyz,6) = -Pxy(ibin,jbin,kbin,ixyz,6)/dble(boxnp(ibin,jbin,kbin))
        enddo
    enddo
    enddo
    enddo

!	Fstresses_mol(:,:) = linearsurface_weight(Pxy,r(:,1:np),Fbinsize,domain, &
!											  shiftmean=flag,meanvalue=Fstresses_CV)
	Fstresses_mol(:,:) = lagrange_poly_weight(Pxy,r(:,1:np),Fbinsize,domain,order=2, &
                                              node_averages=0,shiftmean=2, & 
                                              meanvalue=Fstresses_CV)

end subroutine  get_Fstresses_mol


! ---------------------------------------------------------
! Apply a correction to velocity of the form 
!  u_correct = (u_CV - CFD_u_cnst)*tanh(t/t_correct)
! as a differential constraint in the form of a CV force

subroutine get_F_correction_CV(u_CFD, &
							   binlimits, & 
							   Fcorrect_CV)
	use calculated_properties_MD, only : nbins
	use arrays_MD, only : r, v
	use computational_constants_MD, only : iter, initialstep, delta_t
    use messenger, only : localise_bin
	use physical_constants_MD, only: pi
	use librarymod, only : linspace
	implicit none

	!Perhaps they need a good talking to, if you don't mind my saying so. 
	!Perhaps a bit more. My girls, sir, they didn't care for the Overlook 
	!at first. One of them actually stole a pack of matches, and tried to 
	!burn it down. But I "corrected" them sir. 
	!And when my wife tried to prevent me from doing my duty,
	!I "corrected" her.

	integer,intent(in)             		:: binlimits(6)
	real(kind(0.d0)),intent(in), & 
		allocatable,dimension(:,:,:,:)   :: u_CFD
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:,:,:)   :: Fcorrect_CV

	!While currently being "corrected", don't apply further constraints
	integer												 :: n,i,j,k, rel_iter
    integer                                              :: CVforce_correctime
	integer,parameter									 :: iter_correct=50
	integer,dimension(3)								 :: bmin,bmax,bin
	integer,allocatable,dimension(:,:,:),save			 :: correctstart
	logical,allocatable,dimension(:,:,:),save			 :: correction_lock
	!Machine Precision -- ε such that 1 = 1+ε on a computer
	double precision,parameter							 :: tol = 2.2204460E-16 
	double precision									 :: t,t_correct,that
	double precision,save					 			 :: normalise
	double precision,dimension(3),save					 :: F_sum
	double precision,allocatable,dimension(:,:,:,:) 	 :: u_CV
	double precision,allocatable,dimension(:,:,:,:),save :: u_error

    bmin = max(localise_bin((/binlimits(1),binlimits(3),binlimits(5)/)),(/2,2,2 /))
    bmax = min(localise_bin((/binlimits(2),binlimits(4),binlimits(6)/)),nbins(:)+1)

	!Allocate only on first call
	if (iter-initialstep+1 .eq. CVforce_starttime) then
		allocate(correctstart(nbins(1)+2,nbins(2)+2,nbins(3)+2))
		allocate(correction_lock(nbins(1)+2,nbins(2)+2,nbins(3)+2))
		allocate(u_error(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
		!As the integral is based on discrete points. We need to 
		!normalise by sum to ensure it will be zero for the points we use
		normalise = sum(1.d0/(cosh(linspace(-pi,pi,iter_correct+1))**2.d0))
		correction_lock=.false.
	endif

	allocate(Fcorrect_CV(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); Fcorrect_CV=0.d0
	allocate(u_CV(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); u_CV = 0.d0

	!Start correcting only after constraint is applied
	if (iter .le. CVforce_starttime) then
		return
	endif

	!Get velocity of current CV cells
	u_CV = 0.d0
	do n = 1,np
		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1
		u_CV( bin(1),bin(2),bin(3),:) = u_CV( bin(1),bin(2),bin(3),:) + v(:,n)
	enddo
    !print'(a,4i6,6f10.5)', 'Cell ', iter,3,3,3,u_CFD(3,3,3,:),u_CV(3,3,3,:)

    !Maximum correcting period before stop
!    CVforce_correctime = CVforce_starttime+iter_correct
!	if (iter-initialstep+1 .ge. CVforce_correctime) then
!        Fcorrect_CV = 0.d0
!        return
!    endif

	!Result force du_correct/dt
	Fcorrect_CV = 0.d0
	do i = bmin(1),bmax(2)
	do j = bmin(2),bmax(2)
	do k = bmin(3),bmax(3)
		!Check if CV is currently being corrected
		if (.not. correction_lock(i,j,k)) then
			u_error(i,j,k,:) = u_CV(i,j,k,:) - u_CFD(i,j,k,:)
			if (any(u_error(i,j,k,:) .gt. tol*1000.d0)) then
				correction_lock(i,j,k) = .true.
				correctstart(i,j,k) = iter
                !print'(a,3i6,6f10.5)', 'Correcting ', i,j,k,u_CFD(i,j,k,:),u_CV(i,j,k,:)
			endif
		endif

		!If correction lock is on, calculate correction
		if (correction_lock(i,j,k)) then
			!Get time passed relative to start of correction
			rel_iter = (iter - correctstart(i,j,k))
			t = rel_iter*delta_t
			t_correct = iter_correct*delta_t
			that = t/t_correct-0.5d0
			Fcorrect_CV(i,j,k,:) = (u_error(i,j,k,:))/(delta_t*normalise*cosh(2.d0*pi*that)**2.d0)
			!Switch lock off if end of correction period
			if (rel_iter .eq. iter_correct) correction_lock(i,j,k) = .false.
		endif
	enddo
	enddo
	enddo

end subroutine get_F_correction_CV

! ---------------------------------------------------------

subroutine get_Fmdflux_CV(F_CV, 	 & 
						  F_mol, 	 &
						  MD_flux,   &
						  boxnp,	 &
						  binlimits, &
						  Fmdflux_CV)

	use arrays_MD, only : r, v, a
	use computational_constants_MD, only : delta_t, iter
    use messenger, only : localise_bin
	implicit none

	integer,intent(in)             		:: binlimits(6)
	integer,intent(in), & 
		allocatable,dimension(:,:,:)	:: boxnp
	real(kind(0.d0)),intent(in), & 
		allocatable,dimension(:,:,:,:) 	:: F_CV
	real(kind(0.d0)),intent(in), & 
		allocatable,dimension(:,:)     	:: F_mol
	real(kind(0.d0)),intent(inout), & 
		allocatable,dimension(:,:,:,:,:):: MD_flux
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:,:,:) 	:: Fmdflux_CV

	logical								:: converged
	integer								:: n
	integer								:: maxiter, attempt
	integer,dimension(3)			    :: bin
	integer,dimension(6)				:: localbinlimits
	real(kind(0.d0)),dimension(3)		:: F_iext
	real(kind(0.d0)),allocatable, & 
		dimension(:,:)  				:: r_temp, v_temp
	real(kind(0.d0)),allocatable, & 
		dimension(:,:,:,:) 				:: Fmdflux_CV_prev

    localbinlimits(1:5:2) = max(localise_bin((/binlimits(1), & 
											   binlimits(3), & 
											   binlimits(5)/)),(/2,2,2 /))
    localbinlimits(2:6:2) = min(localise_bin((/binlimits(2), & 
											   binlimits(4), & 
											   binlimits(6)/)),nbins(:)+1)

	allocate( Fmdflux_CV(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); Fmdflux_CV=0.d0
	allocate( Fmdflux_CV_prev(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); Fmdflux_CV_prev=0.d0

	!Get initial value of flux based force
	call flux2force(MD_flux, localbinlimits, Fmdflux_CV)

	!Iterate until momentum flux and applied CV force are consistent
	converged = .false.; maxiter = 100; 
	allocate(r_temp(3,np),v_temp(3,np))
	do attempt = 1,maxiter

    	!Get velocity at next timestep with constraint
    	do n = 1,np
    		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1

			!Get total force on a molecule from sum of CV forces and molecular forces
			F_iext = F_mol(:,n) + ( Fmdflux_CV(bin(1),bin(2),bin(3),:)  &
						              +   F_CV(bin(1),bin(2),bin(3),:) ) & 
								   /dble(boxnp(bin(1),bin(2),bin(3)))

			!Velocity and positions calculated  		
    		v_temp(:,n) = v(:,n) + delta_t * ( a(:,n) - F_iext ) 
    		r_temp(:,n) = r(:,n) + delta_t * v_temp(:,n) 

    	enddo

		!Check convergence
 		if (converged) then
			exit
		elseif (attempt .eq. maxiter-1) then
			print*, "Warning -- CV Force could not converge"
		endif

		!Update momentum flux and force
		call get_MD_fluxes(MD_flux, r_in=r_temp, v_in=v_temp)
		Fmdflux_CV_prev = Fmdflux_CV
        Fmdflux_CV = 0.d0
		call flux2force(MD_flux,localbinlimits,Fmdflux_CV)
		call check_convergence

	enddo
	deallocate(r_temp,v_temp)

contains


	subroutine check_convergence

		use computational_constants_MD, only : iter
		implicit none

		integer					:: global_converged

		!Machine Precision -- ε such that 1 = 1+ε on a computer
		double precision,parameter	:: tol = 2.2204460E-16 

		!If flux force have converged (to tolerence) then exit after calculating final force
		if (all(abs(Fmdflux_CV_prev - Fmdflux_CV) .lt. tol)) then
			converged = .true.
			!Wait for all processors to converge
			global_converged = 1
			call globalMinInt(global_converged)
			if (global_converged .ne. 1) then
				converged = .false.
			else
                !call print_convergence()
			endif

		else
			!Tell other processors "I'm not converged"
			global_converged = 0
			call globalMinInt(global_converged)

			!Ouput convergence if it looks like it's going to be a problem
			if (attempt .gt. 25) then
                call print_convergence()
			endif
		endif

	end subroutine check_convergence

    subroutine print_convergence()
    use messenger, only : globalise_bin, localise_bin
        implicit none

		integer, dimension(3)	:: bmin, bmax
		integer					:: converge_cells
		integer					:: ib, jb, kb
		integer,save			:: convergence_count=0
		double precision		:: convergence

        bmin = max(localise_bin((/binlimits(1),binlimits(3),binlimits(5)/)),(/2,2,2 /))
        bmax = min(localise_bin((/binlimits(2),binlimits(4),binlimits(6)/)),nbins(:)+1)

		convergence = 0.d0; converge_cells = 0
        do ib = bmin(1),bmax(1)
        do jb = bmin(2),bmax(2)
        do kb = bmin(3),bmax(3)
			if (sum(Fmdflux_CV_prev(ib,jb,kb,:)) .ne. sum(Fmdflux_CV(ib,jb,kb,:))) then
				convergence = convergence + sum(abs(Fmdflux_CV_prev(ib,jb,kb,:) - Fmdflux_CV(ib,jb,kb,:)))
                converge_cells = converge_cells + 1
			endif
		enddo
		enddo
		enddo
		convergence_count = convergence_count + 1
		print'(a,2i5,2i8,f28.17)', 'convergence ', iter, attempt, convergence_count, converge_cells, convergence
    end subroutine print_convergence

end subroutine get_Fmdflux_CV

! ---------------------------------------------------------

subroutine apply_force(F_CV,  &
					   F_mol, &
					   boxnp, &
					   binlimits)
	use arrays_MD, only: r, v, a
	use computational_constants_MD, only : eflux_outflag, delta_t, & 
                                           iter,CVweighting_flag
	use module_record_external_forces, only : record_external_forces
	implicit none

	integer,intent(in)             		:: binlimits(6)
	integer,intent(in), & 
		allocatable,dimension(:,:,:)	:: boxnp
	real(kind(0.d0)),intent(in), & 
		allocatable,dimension(:,:,:,:) 	:: F_CV
	real(kind(0.d0)),intent(in), & 
		allocatable,dimension(:,:)     	:: F_mol

	integer								:: n
	integer,dimension(3)				:: bin
	real(kind(0.d0)),dimension(3)		:: F_iext, velvect

	!Loop over all molecules and apply constraint
	do n = 1, np
		!Get bin
		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1

		!Get total force on a molecule from sum of CV forces and molecular forces
		F_iext(:) = F_mol(:,n) + (F_CV(bin(1),bin(2),bin(3),:))/dble(boxnp(bin(1),bin(2),bin(3)))

        !Debug -- write out force applied to molecules every 100 timesteps
		if (mod(iter,100) .eq. 0) then
            if (CVweighting_flag .eq. 0) then
                write(5000+iter,'(2i8,9f18.8)'), iter, n, r(:,n), F_iext, a(:,n)
            elseif (CVweighting_flag .ge. 1) then
                write(5000+iter,'(2i8,9f18.8)'), iter, n, r(:,n), F_mol(:,n), a(:,n)
            endif
		endif
						
		!Apply force by adding to total
		a(:,n) = a(:,n) - F_iext(:)

		!Add external force to CV total
		if (eflux_outflag .eq. 4) then
			velvect(:) = v(:,n) + 0.5d0*delta_t*a(:,n)
			call record_external_forces(F_iext,r(:,n),velvect)
		else
			call record_external_forces(F_iext,r(:,n))
		endif

	enddo

end subroutine apply_force

! ---------------------------------------------------------

subroutine check_limits(CV_limits)
    use computational_constants_MD, only : VOID
    use calculated_properties_MD, only : gnbins
    implicit none

    integer,dimension(6),intent(inout)	:: CV_limits

    !Check for undefined or incorrect limits
    if (any(CV_limits .eq. VOID)) then
	    CV_limits(1) = 1
	    CV_limits(2) = gnbins(1)
	    CV_limits(3) = 1
	    CV_limits(4) = gnbins(2)
	    CV_limits(5) = 1
	    CV_limits(6) = gnbins(3)
    endif
    if (CV_limits(1) .lt. 1        ) CV_limits(1) = 1
    if (CV_limits(2) .gt. gnbins(1)) CV_limits(2) = gnbins(1)
    if (CV_limits(3) .lt. 1        ) CV_limits(3) = 1
    if (CV_limits(4) .gt. gnbins(2)) CV_limits(4) = gnbins(2)
    if (CV_limits(5) .lt. 1        ) CV_limits(5) = 1
    if (CV_limits(6) .gt. gnbins(3)) CV_limits(6) = gnbins(3)

    !Negative maximum are taken from top as with python array syntax
    if (CV_limits(2) .le. 0) CV_limits(2) = gnbins(1)+CV_limits(2)
    if (CV_limits(4) .le. 0) CV_limits(4) = gnbins(2)+CV_limits(4)
    if (CV_limits(6) .le. 0) CV_limits(6) = gnbins(3)+CV_limits(6)

end subroutine check_limits


end module apply_CV_force_mod


subroutine apply_CV_force
	use apply_CV_force_mod
	use computational_constants_MD, only : iter,initialstep,  F_CV_limits, & 
										   CVforce_correct, CVweighting_flag, & 
										   VOID, CVforce_flag
	implicit none

	integer,allocatable,dimension(:,:,:) 				:: boxnp
	real(kind(0.d0)),allocatable,dimension(:,:)     	:: Fstresses_mol
	real(kind(0.d0)),allocatable,dimension(:,:,:,:) 	:: F_CV, u_CFD
	real(kind(0.d0)),allocatable,dimension(:,:,:,:) 	:: Fmdflux_CV, Fstresses_CV
	real(kind(0.d0)),allocatable,dimension(:,:,:,:) 	:: Fcfdflux_CV, Fcorrect_CV
	real(kind(0.d0)),allocatable,dimension(:,:,:,:,:) 	:: CFD_stress,CFD_flux
	real(kind(0.d0)),allocatable,dimension(:,:,:,:,:) 	:: MD_stress,MD_flux

	!Check if CV force is turned on
	if (CVforce_flag .eq. VOID) return

	!Check for start time 
	if (iter-initialstep+1 .lt. CVforce_starttime) then
        return
    elseif (iter-initialstep+1 .eq. CVforce_starttime) then
	    Fbinsize = domain/nbins
	    volume = product(Fbinsize)
        call check_limits(F_CV_limits)
    endif

	!Get CFD stresses and fluxes
	call get_CFD_velocity(F_CV_limits, u_CFD)
	call get_CFD_stresses_fluxes(F_CV_limits, CFD_stress, CFD_flux)

	!Get molecules per bin, MD stresses and MD fluxes
	call get_boxnp(boxnp)
	call get_MD_stresses_fluxes(F_CV_limits, MD_stress, MD_flux)

	! Get CV force based on CFD flux
	call get_Fcfdflux_CV(CFD_flux, F_CV_limits, Fcfdflux_CV)

	! Get force required to correct velocity setpoint if drift occurs
	if (CVforce_correct .eq. 0) then
		allocate(Fcorrect_CV(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); Fcorrect_CV = 0.d0
	elseif (CVforce_correct .eq. 1) then
		call get_F_correction_CV(u_CFD, F_CV_limits, Fcorrect_CV)
	endif

	!Apply CV stresses as constant force to whole volume or distribute
	allocate(F_CV(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); F_CV = 0.d0

	! Force per molecule is zero and Fstresses added to CV force
	if (CVweighting_flag .eq. 0) then
		! Get CV force based on difference in CFD and MD stresses
		call get_Fstresses_CV(CFD_stress, MD_stress, F_CV_limits, Fstresses_CV)
		allocate(Fstresses_mol(3,np)); Fstresses_mol = 0.d0
		F_CV = Fcfdflux_CV + Fstresses_CV + Fcorrect_CV

	! Get Fstresses per molecule based on distribution of stress 
	! and don't add Fstresses_CV to total
	elseif (CVweighting_flag .eq. 1 .or. & 
			CVweighting_flag .eq. 2) then
		call get_Fstresses_mol(CFD_stress, MD_stress, boxnp, CVweighting_flag, F_CV_limits, Fstresses_mol)
		F_CV = Fcfdflux_CV  + Fcorrect_CV
    else
        stop "Error in apply_CV_force - incorrect CVweighting_flag"
	endif

	! Get CV force based on MD fluxes (Using all other forces and iterating until
	! 								   MD_flux and force are consistent)
	call get_Fmdflux_CV(F_CV, Fstresses_mol, MD_flux, boxnp, F_CV_limits, Fmdflux_CV)

	!Add up all forces and apply 
	F_CV = F_CV + Fmdflux_CV
	call apply_force(F_CV, Fstresses_mol, boxnp, F_CV_limits)

end subroutine apply_CV_force



!----------------------------------------------------------------------------------
! Set velocity of a reange of bins to prescribed value

subroutine set_mol_velocity(mol_list, listnp, velocity)
	use linked_list, only : cell, node
	use arrays_MD, only : v, r  !TEMP remove r here
	use computational_constants_MD, only : binspercell,nhb,iter  !TEMP remove nhb,iter here
	implicit none

	integer,intent(in)							:: listnp
	integer,dimension(listnp), intent(in)       :: mol_list
	double precision,dimension(3),intent(in)	:: velocity   !Overall momentum of system

	integer										:: i,molno,binNsum
	double precision,dimension(3)				:: binvsum, vcorrection

	if (listnp .eq. 0) return

	binNsum = 0; binvsum = 0.d0
	do i = 1,listnp	!Step through each particle in list 
		molno = mol_list(i) 	!Number of molecule
		binNsum = binNsum + 1    
		binvsum(:) = binvsum(:) + v(:,molno)

		!print'(i5,a,5i6,6f10.5)',iter,' velocities ',i,size(mol_list,1),molno,binNsum,listnp, r(:,molno), v(:,molno)
	enddo
	
	!Calculate velocity correction per molecule
	vcorrection(:) = binvsum(:)/binNsum - velocity(:)

	!print'(3(a,3f10.5))', 'applied v ', velocity, ' bin v ', binvsum(:)/binNsum, ' v correct ', vcorrection

	!Apply velocity correction per molecule to bin
	binNsum = 0; binvsum = 0.d0
	do i = 1,listnp	!Step through each particle in list 
		molno = mol_list(i) 	!Number of molecule
		v(:,molno) =  v(:,molno) - vcorrection
		binNsum = binNsum + 1    
		binvsum(:) = binvsum(:) + v(:,molno)
	enddo

	!print'(a,3f10.5)', 'Corrected velocity is then ',  binvsum(:)/binNsum

end subroutine set_mol_velocity



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

	integer,optional,intent(in)			:: start_iter
	double precision, intent(in)			:: targetradius
	double precision, dimension(3), intent(in)	:: centre

	integer						:: n, start
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


