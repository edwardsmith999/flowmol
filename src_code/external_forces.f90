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
	use computational_constants_MD, only: domain,CVforce_starttime
	use physical_constants_MD, only: np
	implicit none

	real(kind(0.d0))				::	volume
	real(kind(0.d0)),dimension(3)	::	Fbinsize

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

	!Total CFD flux
	!print*, 'flux2force', bl, shape(flux), shape(force), & 
	!		shape(flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,1)), shape(force(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:))
	force(:,:,:,:) =  ( flux(:,:,:,:,1)-flux(:,:,:,:,4)) &
					 +( flux(:,:,:,:,2)-flux(:,:,:,:,5)) &
					 +( flux(:,:,:,:,3)-flux(:,:,:,:,6))
		 

!	force(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:) & 
!		  =  ( flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,1)  & 
!			  -flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,4)) &
!		    +( flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,2)  & 
!			  -flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,5)) &
!		    +( flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,3)  & 
!			  -flux(bl(1):bl(2),bl(3):bl(4),bl(5):bl(6),:,5))

end subroutine

! ---------------------------------------------------------

subroutine get_MD_stresses(MD_stress)
	use CV_objects, only : 	CV_constraint
	implicit none

	real(kind(0.d0)),intent(inout), & 
		allocatable,dimension(:,:,:,:,:) 	:: MD_stress

	MD_stress = 0.25d0 * CV_constraint%Pxy
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
	MD_flux = CV_constraint%flux/delta_t

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
	!use , only : nbins
	implicit none

	integer,dimension(6),intent(in)		:: binlimits(6)
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:,:,:) 	:: u_CFD

	allocate(u_CFD(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); u_CFD=0.d0

	!This should be replaced by coupler call
	u_CFD = u_CFD

end subroutine get_CFD_velocity

subroutine get_CFD_stresses_fluxes(binlimits, & 
								   CFD_stress, & 
								   CFD_flux)
    use messenger, only : localise_bin
	implicit none

	integer,dimension(6),intent(in)			:: binlimits(6)
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:,:,:,:) 	:: CFD_stress,CFD_flux


	logical 								:: timeevolve = .false.
	integer									:: appliedbins, taumin,taumax,jb
	integer,dimension(6)				 	:: localbinlimits
	real(kind(0.d0)), & 
		allocatable,dimension(:) 			:: tautemp

	allocate(CFD_stress(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6)); CFD_stress=0.d0
	allocate(  CFD_flux(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6)); CFD_flux=0.d0

	! v v v v v   This should be replaced by coupler call   v v v v v 
!    localbinlimits(1:5:2) = max(localise_bin((/binlimits(1), & 
!											   binlimits(3), & 
!											   binlimits(5)/)),(/2,2,2 /))
!    localbinlimits(2:6:2) = min(localise_bin((/binlimits(2), & 
!											   binlimits(4), & 
!											   binlimits(6)/)),nbins(:)+1)

!	taumin = 2.d0 * (binlimits(3)-1) - 1
!	taumax = 2.d0 * (binlimits(4)  ) - 1
!	appliedbins = taumax-taumin+1
!	allocate(tautemp(appliedbins))	!Double resolution for surfaces

!	! Couette Analytical solution
!	if (timeevolve) then
!		!Time evolving solution
!		!tautemp = couette_analytical_stress_fn((iter-CVforce_starttime)*delta_t,Re,1.d0,H,appliedbins,2)
!	else
!		! Steady state solution
!		tautemp = 2.d0/domain(2)
!	endif

!	!Copy tau to CFD surfaces
!	CFD_flux = 0.d0
!	CFD_stress  = 0.d0
!	do jb = binlimits(3),binlimits(4)
!		!CFD_stress(:,jb,:,1,2) = tautemp(2.d0*(jb-1)-1)	!Top
!		!CFD_stress(:,jb,:,1,5) = tautemp(2.d0*(jb  )-1)	!Bottom
!	    if (mod(jb,2) .eq. 0) then
!			CFD_stress(:,jb,:,1,2) =  tautemp(2.d0*(jb-1)-1)	!Top
!			CFD_stress(:,jb,:,1,5) =  0.d0						!Bottom
!		else
!			CFD_stress(:,jb,:,1,2) =  0.d0						!Top
!			CFD_stress(:,jb,:,1,5) =  tautemp(2.d0*(jb  )-1)	!Bottom
!		endif
!		!print*, jb, tautemp(2.d0*(jb-1)-1), CFD_stress(3,jb,3,1,2)
!	enddo

	! ^ ^ ^ ^ ^   This should be replaced by coupler call   ^ ^ ^ ^ ^

	CFD_stress = CFD_stress*volume
	CFD_flux   = CFD_flux*volume


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
	allocate(Pxy( nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6))

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
	use librarymod, only : linearsurface_weight
	use computational_constants_MD, only : iter
	implicit none

	integer,intent(in)             		 :: binlimits(6), flag
	integer,intent(in), & 
		allocatable,dimension(:,:,:)	:: boxnp
	real(kind(0.d0)),intent(in), & 
		allocatable,dimension(:,:,:,:,:) :: CFD_stress,MD_stress
	real(kind(0.d0)),intent(out), & 
		allocatable,dimension(:,:)   	 :: Fstresses_mol

	integer								 :: n, bin(3), ixyz, s
	real(kind(0.d0)),allocatable,& 
		dimension(:,:,:,:,:) 			 :: Pxy
	real(kind(0.d0)), & 
		allocatable,dimension(:,:,:,:)   :: Fstresses_CV

	allocate(Fstresses_mol(3,np));  Fstresses_mol=0.d0

	allocate(Fstresses_CV( nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); Fstresses_CV=0.d0
	allocate(Pxy( nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6))

	!Get F_stresses on the CV
	call get_Fstresses_CV(CFD_stress, MD_stress, binlimits, Fstresses_CV)
	Fstresses_CV(:,:,:,1)=Fstresses_CV(:,:,:,1)*dble(boxnp)
	Fstresses_CV(:,:,:,2)=Fstresses_CV(:,:,:,2)*dble(boxnp)
	Fstresses_CV(:,:,:,3)=Fstresses_CV(:,:,:,3)*dble(boxnp)

	Pxy = MD_stress-CFD_stress
	print*, 'WARNING -- STRSS IN  get_Fstresses_mol HAS SIGN REVERSED'
	Pxy(:,:,:,:,1:3) = -Pxy(:,:,:,:,1:3)
	Pxy(:,:,:,:,4:6) = Pxy(:,:,:,:,4:6)
	print*, 'WARNING -- STRSS IN  get_Fstresses_mol HAS SIGN REVERSED'
	Fstresses_mol(:,:) = linearsurface_weight(Pxy,r(:,1:np),Fbinsize,domain, &
											  shiftmean=flag,meanvalue=Fstresses_CV)

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

	!Result force du_correct/dt
	Fcorrect_CV = 0.d0
	do i = binlimits(1),binlimits(2)
	do j = binlimits(3),binlimits(4)
	do k = binlimits(5),binlimits(6)
		!Check if CV is currently being corrected
		if (.not. correction_lock(i,j,k)) then
			u_error(i,j,k,:) = u_CV(i,j,k,:) - u_CFD(i,j,k,:)
			if (any(u_error(i,j,k,:) .gt. tol*1000.d0)) then
				correction_lock(i,j,k) = .true.
				correctstart(i,j,k) = iter
!				if (i .eq. 3 .and. j .eq. 3 .and. k .eq. 3) then
!					F_sum = 0.d0
!					print'(a,4i5,7f10.5)', 'New error detected',iter,i,j,k,u_CV(i,j,k,:),u_CFD(i,j,k,:),tol*1000.d0
!				endif
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
!			if (i .eq. 3 .and. j .eq. 3 .and. k .eq. 3) then
!				F_sum = F_sum + u_error(i,j,k,:)/(normalise*cosh(2.d0*pi*that)**2.d0)
!				print'(a,2i5,2f10.5,4i3,6f10.4,4f12.5)', 'Corr',iter,rel_iter,that,1.d0/cosh(2.d0*pi*that)**2.d0, & 
!														  correctstart(i,j,k),i,j,k,u_error(i,j,k,:),Fcorrect_CV(i,j,k,:),F_sum(:),normalise
!			endif
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
	use computational_constants_MD, only : delta_t
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
			F_iext = (F_mol(:,n) + F_CV(bin(1),bin(2),bin(3),:)  &
						   + Fmdflux_CV(bin(1),bin(2),bin(3),:))/dble(boxnp(bin(1),bin(2),bin(3)))

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
		!print*, 'converge check', attempt, Fmdflux_CV(3,3,3,:), Fmdflux_CV_prev(3,3,3,:)
		call flux2force(MD_flux,localbinlimits,Fmdflux_CV)
		call check_convergence

	enddo
	deallocate(r_temp,v_temp)

contains


	subroutine check_convergence
	    use messenger, only : globalise_bin, localise_bin
		use computational_constants_MD, only : iter
		implicit none

		integer					:: ib, jb, kb
		integer					:: global_converged
		integer					:: converge_cells
		integer,save			:: convergence_count=0
		integer, dimension(3)	:: bmin, bmax
		!Machine Precision -- ε such that 1 = 1+ε on a computer
		double precision			:: convergence
		double precision,parameter	:: tol = 2.2204460E-16 

		bmin = max(localise_bin((/binlimits(1),binlimits(3),binlimits(5)/)),(/2,2,2 /))
		bmax = min(localise_bin((/binlimits(2),binlimits(4),binlimits(6)/)),nbins(:)+1)

		!If flux force have converged (to tolerence) then exit after calculating final force
		if (all(abs(Fmdflux_CV_prev - Fmdflux_CV) .lt. tol)) then
			converged = .true.
			!Wait for all processors to converge
			global_converged = 1
			call globalMinInt(global_converged)
			if (global_converged .ne. 1) then
				converged = .false.
			else
!	            convergence = 0.d0; converge_cells = 0
!		        do ib = bmin(1),bmax(1)
!		        do jb = bmin(2),bmax(2)
!		        do kb = bmin(3),bmax(3)
!				    convergence = convergence + sum(abs(Fmdflux_CV_prev(ib,jb,kb,:) - Fmdflux_CV(ib,jb,kb,:)))
!    	            converge_cells = converge_cells + 1
!			    enddo
!			    enddo
!			    enddo
!				print'(a,2i5,2i8,f28.17)', 'converged ', iter, attempt, convergence_count, converge_cells, convergence
			endif

		else
			!Tell other processors "I'm not converged"
			global_converged = 0
			call globalMinInt(global_converged)

			!Ouput convergence if it looks like it's going to be a problem
			if (attempt .gt. 25) then
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
			endif
		endif

!		if (attempt .eq. maxiter-1 .and. converged .eqv. .false.) then
!			print*, "Warning -- CV Force could not converge"
!			continue
!		endif

	end subroutine check_convergence


end subroutine get_Fmdflux_CV

! ---------------------------------------------------------

subroutine apply_force(F_CV,  &
					   F_mol, &
					   boxnp, &
					   binlimits)
	use arrays_MD, only: r, v, a
	use computational_constants_MD, only : eflux_outflag, delta_t,iter
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
		F_iext(:) = (F_mol(:,n) + F_CV(bin(1),bin(2),bin(3),:))/dble(boxnp(bin(1),bin(2),bin(3)))

		if (mod(iter,100) .eq. 0) then
			bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1
			write(5000+iter,'(2i8,9f18.8)'), iter, n, r(:,n), F_iext, a(:,n)
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


end module apply_CV_force_mod


subroutine apply_CV_force
	use apply_CV_force_mod
	use computational_constants_MD, only : iter,initialstep,  F_CV_limits, & 
										   CVforce_correct, CVweighting_flag, & 
										   VOID, CVforce_flag
	use calculated_properties_MD, only : gnbins
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
	if (iter-initialstep+1 .lt. CVforce_starttime) return

	!Test case focuses on a single CV
	Fbinsize = domain/nbins
	volume = product(Fbinsize)

	if (any(F_CV_limits .eq. VOID)) then
		F_CV_limits(1) = 1
		F_CV_limits(2) = gnbins(1)
		F_CV_limits(3) = 1
		F_CV_limits(4) = gnbins(2)
		F_CV_limits(5) = 1
		F_CV_limits(6) = gnbins(3)
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
	endif

	! Get CV force based on MD fluxes (Using all other forces and iterating until
	! 								   MD_flux and force are consistent)
	call get_Fmdflux_CV(F_CV, Fstresses_mol, MD_flux, boxnp, F_CV_limits, Fmdflux_CV)

	!Add up all forces and apply 
	F_CV = F_CV + Fmdflux_CV
	call apply_force(F_CV, Fstresses_mol, boxnp, F_CV_limits)

end subroutine apply_CV_force

























!=========================================================================
! Routine to test flux calculations by various methods
! with iteration over all bins to convergence

subroutine apply_CV_force_multibin(iter)
    use physical_constants_MD, only : density
	use computational_constants_MD, only : irank, npy, globaldomain, CVforce_flag,CVforce_starttime, & 
										   CVweighting_flag, VOID,  Nsteps, iblock, jblock, kblock, F_CV_limits
	use calculated_properties_MD, only : pressure, nbins, gnbins
	use librarymod, only : get_new_fileunit, couette_analytical_fn, lagrange_poly_weight_Nmol, linearsurface_weight
	use module_external_forces, only : np, momentum_flux, irank,nbins, nbinso,domain,delta_t,Nvflux_ave
	use CV_objects, only : 	CV_constraint
	use set_bin_velocity_mod
	implicit none
	
	integer,intent(in)				:: iter

	logical							:: apply_CVforce = .false.
	integer							:: ibin,jbin,kbin,starttime
	integer                     	:: igmin,jgmin,kgmin,igmax,jgmax,kgmax
	real(kind(0.d0))								:: volume, y_loc
	real(kind(0.d0)),dimension(3)					:: Fbinsize, u_bin
	real(kind(0.d0)),allocatable,dimension(:,:)     :: F_dist
	real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: F_constraint, F_correct
	real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: CFD_Pi_dS,CFD_rhouu_dS,CFD_u_cnst

	!Used for force applied only to non-crossing molecules
	integer,dimension(:),allocatable				:: notcrossing

	if (CVforce_flag .eq. VOID) return

	!Test case focuses on a single CV
	Fbinsize = domain/nbins
	volume = product(Fbinsize)

	starttime = CVforce_starttime  !Nsteps/2.d0

	if (any(F_CV_limits .eq. VOID)) then
		F_CV_limits(1) = 1
		F_CV_limits(2) = gnbins(1)
		F_CV_limits(3) = 1
		F_CV_limits(4) = gnbins(2)
		F_CV_limits(5) = 1
		F_CV_limits(6) = gnbins(3)
	endif

    igmin = 1; jgmin = 1; kgmin = 1;
    igmax = gnbins(1); jgmax = nint(gnbins(2)/2.d0)+1; kgmax = gnbins(3)

	!if (iter .lt. starttime) return

!    igmin = nint(gnbins(1)/2.d0); jgmin = nint(gnbins(2)/2.d0); kgmin = nint(gnbins(3)/2.d0);
!    igmax = nint(gnbins(1)/2.d0); jgmax = nint(gnbins(2)/2.d0); kgmax = nint(gnbins(3)/2.d0)

	!Exchange CV data ready to apply force
	call CV_constraint%swap_halos(nbinso)

	!Get average over current cell and apply constraint forces
	!call get_continuum_values
	call get_test_forces(CVforce_flag, F_CV_limits(1),F_CV_limits(3),F_CV_limits(5), & 
									   F_CV_limits(2),F_CV_limits(4),F_CV_limits(6))
	call get_distributed_force(CVweighting_flag,  F_CV_limits(1),F_CV_limits(3),F_CV_limits(5), & 
									   			  F_CV_limits(2),F_CV_limits(4),F_CV_limits(6),F_dist)

	! Apply correction between current velocity and required velocity by adding
	! the derivative of a sigmoid type function to the CV constraint force
	call add_velocity_correction_forces(CVforce_flag, igmin,jgmin,kgmin,igmax,jgmax,kgmax)

	!Set velocity to required value and then apply constraint 
	!force to adjust its evolution from then on...
	if (iter .eq. CVforce_starttime .and. apply_CVforce .ne. .false.) then
		!This routine takes the global bin number which is one less than the local 
		!used in other routines as these include halos
!  		call set_bin_velocity(ibin-1, ibin-1, & 
!  							  jbin-1, jbin-1, & 
!  							  kbin-1, kbin-1, & 
!   							  (/ 0.0d0,0.d0,0.d0 /),0)

		do jbin = F_CV_limits(3), F_CV_limits(4)
    	do ibin = F_CV_limits(1), F_CV_limits(2)
    	do kbin = F_CV_limits(5), F_CV_limits(6)

!			y_loc = (jbin-1)*Fbinsize(2)
            if (CVforce_flag .eq. 5) then
                if (jbin .eq. jgmin) then
                    u_bin = (/ -1.d0*density,0.d0, 0.d0 /)
                elseif(jbin .eq. jgmax) then
                    u_bin = (/  1.d0*density,0.d0, 0.d0 /)
                else
	    		    u_bin = (/  0.d0,0.d0, 0.d0 /)
                endif
            else
       		    u_bin = (/  0.d0,0.d0, 0.d0 /)
            endif

!			u_bin = (/ (jbin-1)*2.d0/(nbins(2)-1)-1.0 ,0.d0, 0.d0 /)
!			!u_bin = (/ 0.25*sin(2.d0*3.14159*y_loc/globaldomain(2)),0.d0, 0.d0 /)
!			call set_bin_velocity(ibin, ibin, & 
!								  jbin, jbin, & 
!								  kbin, kbin, & 
! 								  u_bin,1)
		enddo
		enddo
		enddo

	endif

	if (iter .ge. CVforce_starttime) then

		!Retreive CV data and calculate force to apply
		call average_over_allbins_iter(CVforce_flag, igmin,jgmin,kgmin,igmax,jgmax,kgmax)
		!call average_over_allbins_noFcrossing(CVforce_flag)

		!Apply the force
		!call apply_force
		call apply_force_weighted

	endif

	!Reset CV force values
	CV_constraint%flux = 0.d0; CV_constraint%Pxy = 0.d0

contains

! Return a distributed force based on random or stress based distributions

subroutine get_distributed_force(flag,igmin,jgmin,kgmin,igmax,jgmax,kgmax,F_dist_)
	use arrays_MD, only : r
    use messenger, only : globalise_bin, localise_bin
	implicit none

	integer,intent(in)	:: flag,igmin,jgmin,kgmin,igmax,jgmax,kgmax
	real(kind(0.d0)),allocatable,dimension(:,:),intent(out)     :: F_dist_

	integer         	:: bmin(3),bmax(3)
	double precision	:: Re, H

    bmin = max(localise_bin((/igmin,jgmin,kgmin/)),(/2,2,2 /))
    bmax = min(localise_bin((/igmax,jgmax,kgmax/)),nbins(:)+1)

	if (allocated(F_dist_)) deallocate(F_dist_)
	select case(CVweighting_flag)
	case(0)
		!No weighting -- apply evenly. 
		allocate(F_dist_,source=get_uniform_F_dist())
	case(1)
		!Random weighting
		allocate(F_dist_,source=get_random_F_dist())
	case(2)
		!MD stress weighting
		allocate(F_dist_,source=get_MDstress_F_dist(r,Fbinsize,domain))
	case(3)
		!CFD Couette flow weighting only
        Re = 0.33333d0
        H = globaldomain(2)
		allocate(F_dist_,source=get_Couette_stress_F_dist(r,Re,H,bmin,bmax,Fbinsize,domain,timeevolve=.false.))
	case(4)
		!MD stress and Couette flow weighting
        Re = 0.33333d0
        H = globaldomain(2)
		allocate(F_dist_,source=get_MD_and_Couette_stress_F_dist(r,Re,H,bmin,bmax,Fbinsize,domain))
	case(5)
		stop "DEBUG CASE 5 -- Doesn't work at present"
		!allocate(F_dist_,source=get_parabolic_F_dist_(r,Fbinsize,domain))
	case default 
		stop "CVweighting_flag Error -- must be 0 to 3"
	end select

end subroutine get_distributed_force


	!No distribution added to the uniform F_constraint
	function get_uniform_F_dist() result(w)
		implicit none

		double precision,dimension(:,:),allocatable	:: w

		allocate(w(3,np))
		w = 0.d0
	
	end function get_uniform_F_dist

	!Random distribution added to the uniform F_constraint
	function get_random_F_dist() result(w)
		implicit none

		double precision,dimension(:,:),allocatable	:: w

		allocate(w(3,np))
		call random_number(w)

		!Normalise total contribution to zero
		w(1,:) = w(1,:) - w(1,:)/sum(w(1,:))
		w(2,:) = w(2,:) - w(2,:)/sum(w(2,:))
		w(3,:) = w(3,:) - w(3,:)/sum(w(3,:))
	
	end function get_random_F_dist

	!MD stress based distribution added to the uniform F_constraint
	function get_MDstress_F_dist(r_in,binsize,domain) result(w)
		use arrays_MD, only : r,v,a
		use cumulative_momentum_flux_mod, only : cumulative_momentum_flux
		use computational_constants_MD, only : halfdomain
		implicit none

		double precision,dimension(3),intent(in)	:: binsize,domain
		double precision,dimension(:,:),intent(in)	:: r_in

		integer	:: n,i,j,k, bin(3)
		double precision,dimension(3)               ::    rand, w_ave
		double precision,dimension(:,:),allocatable	:: w,r_temp,v_temp
		double precision,allocatable,dimension(:,:,:,:) :: wsum
		double precision,allocatable,dimension(:,:,:,:,:) :: MD_stress

		! N.B. Only the configurational stress is used to determine the distribution of forces as
		! application of forces to compensate for instantanous molecules crossings 
		! doesn't make physical sense (I believe)

		allocate(w(3,np))
		allocate(MD_stress(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6))
		MD_stress = 0.25d0*CV_constraint%Pxy
		w(:,:) = linearsurface_weight(MD_stress,r_in(:,1:np),binsize,domain,shiftmean=2)

		!Call linear lagrange polynomial to distribute forces
		!w(:,:) = lagrange_poly_weight(MD_stress,r_in(:,1:np),binsize,domain,2)

	end function get_MDstress_F_dist

	function get_Couette_stress_F_dist(r_in,Re,H,bmin,bmax,binsize,domain,timeevolve) result(w)
		use librarymod, only : couette_analytical_stress_fn
		implicit none

		integer,dimension(3),intent(in)	:: bmin,bmax
		double precision,dimension(3),intent(in)	:: binsize,domain
		double precision,intent(in)		:: Re, H
		double precision,dimension(:,:),intent(in)	:: r_in
		logical,optional							:: timeevolve

		integer						:: i,j,k,appliedbins, ib,jb,kb,taumin, taumax

		double precision,dimension(:),allocatable	:: tautemp
		double precision,dimension(:,:),allocatable	:: w
		double precision,allocatable,dimension(:,:,:,:,:) :: CFDstress
		real(kind(0.d0)),allocatable,dimension(:,:,:,:,:) :: CFD_Pxy,CFD_flux

		if (.not. present(timeevolve)) then
			timeevolve = .false.
		endif

		allocate(CFD_Pxy( nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6))
		allocate(CFD_flux(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6))
		allocate(CFDstress(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6))

		!appliedbins = 2.d0*gnbins(2)+1
		taumin = 2.d0 * (bmin(2)-1) - 1
		taumax = 2.d0 * (bmax(2)  ) - 1
		appliedbins = taumax-taumin+1
		allocate(tautemp(appliedbins))	!Double resolution for surfaces

		! Couette Analytical solution
		if (timeevolve) then
			!Time evolving solution
			tautemp = couette_analytical_stress_fn((iter-CVforce_starttime)*delta_t,Re,1.d0,H,appliedbins,2)
		else
			! Steady state solution
			tautemp = 2.d0/domain(2)
		endif

		!Copy tau to CFD surfaces
		CFD_flux = 0.d0
		CFD_Pxy  = 0.d0
		do jb = bmin(2),bmax(2)
			!CFD_Pxy(:,jb,:,1,2) =  tautemp(2.d0*(jb-1)-1)	!Top
			!CFD_Pxy(:,jb,:,1,5) =  tautemp(2.d0*(jb  )-1)	!Bottom
		    if (mod(jb,2) .eq. 0) then
				CFD_Pxy(:,jb,:,1,2) =  tautemp(2.d0*(jb-1)-1)	!Top
				CFD_Pxy(:,jb,:,1,5) =  0.d0						!Bottom
			else
				CFD_Pxy(:,jb,:,1,2) =  0.d0						!Top
				CFD_Pxy(:,jb,:,1,5) =  tautemp(2.d0*(jb  )-1)	!Bottom
			endif
		enddo
	!    do ib = bmin(1),bmax(1)
	!    do kb = bmin(3),bmax(3)
	!		CFD_Pxy(ib,bmin(2):bmax(2),kb,1,2) =  tautemp(taumin:taumax-2:2)	!Top
	!		CFD_Pxy(ib,bmin(2):bmax(2),kb,1,5) =  tautemp(taumin+2:taumax:2)	!Bottom
	!    enddo
	!    enddo

		!Take CFD stress field and use to determine weighting function
		CFDstress(:,:,:,:,:) = -(CFD_Pxy(:,:,:,:,:)+CFD_flux(:,:,:,:,:))*volume

		!Call linear lagrange polynomial to distribute forces
		allocate(w(3,np))
		w(:,:) = linearsurface_weight(CFDstress,r_in(:,1:np),binsize,domain,shiftmean=2)
		!w(:,:) = lagrange_poly_weight_Nmol(CFDstress,r_in(:,1:np),binsize,domain,2,zeromean=.true.)
		w(2:3,:) = 0.d0
	
	end function get_Couette_stress_F_dist


	function get_MD_and_Couette_stress_F_dist(r_in,Re,H,bmin,bmax,binsize,domain) result(w)
		use librarymod, only : couette_analytical_stress_fn
		implicit none

		integer,dimension(3),intent(in)	:: bmin,bmax
		double precision,dimension(3),intent(in)	:: binsize,domain
		double precision,intent(in)		:: Re, H
		double precision,dimension(:,:),intent(in)	:: r_in
		double precision,dimension(:,:),allocatable	:: w
		stop "get_MD_and_Couette_stress_F_dist IS NOT DEVELOPED YET"
	
	end function get_MD_and_Couette_stress_F_dist


! Apply arbitary forces for testing purposes
! min and max values are specified in the global bin numbering
subroutine get_test_forces(flag,igmin,jgmin,kgmin,igmax,jgmax,kgmax)
    use physical_constants_MD, only : pi, tethereddisttop, tethereddistbottom, density
    use computational_constants_MD, only : npx, npy, npz, Nvel_ave, tplot,CVweighting_flag
	use arrays_MD, only : r
    use messenger, only : globalise_bin, localise_bin
	implicit none

	integer,intent(in)	:: flag,igmin,jgmin,kgmin,igmax,jgmax,kgmax

	integer         	:: bmin(3),bmax(3),appliedbins
	integer				:: ib,jb,kb,length,fileunit
	double precision	:: sin_mag, sin_period, Re, H

	double precision,dimension(:),allocatable	:: u_t,u_mdt,utemp

    bmin = max(localise_bin((/igmin,jgmin,kgmin/)),(/2,2,2 /))
    bmax = min(localise_bin((/igmax,jgmax,kgmax/)),nbins(:)+1)

	allocate(CFD_Pi_dS(   nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); CFD_Pi_dS=0.d0
	allocate(CFD_rhouu_dS(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); CFD_rhouu_dS=0.d0
	allocate(CFD_u_cnst(  nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); CFD_u_cnst=0.d0

	!Set velocity to required value and then apply constraint 
	!force to adjust its evolution from then on...
	if (iter .eq. CVforce_starttime .and. flag .ne. 0) then
		apply_CVforce = .true.
	elseif (iter .le. CVforce_starttime) then
		apply_CVforce = .false.
	else
		!Do nothing
	endif

    select case(flag)
    case(0)
	    !No Force applied
	    apply_CVforce = .false.
    case(1)
	    !Zero Force applied
	    CFD_Pi_dS(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),:) = 0.d0
	    CFD_rhouu_dS(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),:) = 0.d0	
	    !Velocity is a constant (assumed zero)
	    CFD_u_cnst(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),:) = 0.d0
    case(2)
	    !Constant function 
	    CFD_Pi_dS(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),:) = 1.d0
	    CFD_rhouu_dS(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),:) = 0.d0
	    !Velocity is a linear fn of time
	    CFD_u_cnst(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),1) = iter
    case(3)
	    CFD_Pi_dS(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),:) = 0.d0
	    !Sin function shifted by pi/2 so velocity which is given
	    !by the integral (cos) alternates around zero
	    sin_mag = 0.1d0;  sin_period = 100.d0
	    CFD_Pi_dS(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),1) = sin_mag * cos(2.d0*pi*((iter-CVforce_starttime)/sin_period))! - 0.505*pi)
	    CFD_rhouu_dS = 0.d0

	    ! - - - - -  Sinusoid of velocity with NCER special discretisation - - - - - - - - - 
	    CFD_u_cnst(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),:) = 0.d0
	    CFD_u_cnst(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),1) = & 
			1.25d0*delta_t*(sin_mag*sin_period/(2.d0*pi))* sin(2.d0*pi*((iter-CVforce_starttime)/sin_period))

    case(4)
	    ! Get continuum values of surface stresses, etc
	    !Read Shear pressure from file...
	    CFD_Pi_dS(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),:) = 0.d0
	    fileunit = get_new_fileunit()
	    inquire(iolength=length) CFD_Pi_dS(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),1)
	    open(unit=fileunit,file='./F_hist',form='unformatted', access='direct',recl=length)
	    read(fileunit,rec=iter) CFD_Pi_dS(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),1)
	    close(fileunit,status='keep')
	    CFD_rhouu_dS = 0.d0
	    ! Read velocity too for NCER style proportional constraint
	    fileunit = get_new_fileunit()
	    inquire(iolength=length) CFD_u_cnst(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),1)
	    open(unit=fileunit,file='./u_hist',form='unformatted', access='direct',recl=length)
	    read(fileunit,rec=iter) CFD_u_cnst(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3),1)
	    close(fileunit,status='keep')

    case(5)

	    ! Couette Analytical solution applied as 2 seperate regions to top and bottom of domain...
        Re = 0.33333d0

        !Range of bins applied force is spread over
        appliedbins = nint(gnbins(2)/2.d0)+1  !Half domain
        H = globaldomain(2)/2.d0 !Half domain

        allocate(utemp(gnbins(2)),u_t(bmin(2):bmax(2)),u_mdt(bmin(2):bmax(2)))
        CFD_rhouu_dS = 0.d0
	    CFD_Pi_dS    = 0.d0
        CFD_u_cnst   = 0.d0

        ! Get velocity from analytical solution first
        if (iter .gt. CVforce_starttime+1) then
            utemp(:) = couette_analytical_fn( (iter-CVforce_starttime) *delta_t,Re,1.d0,H,appliedbins,2)
            u_t(:)   = utemp(jgmin:jgmax)*density
            utemp(:) = couette_analytical_fn((iter-CVforce_starttime-1)*delta_t,Re,1.d0,H,appliedbins,2)
            u_mdt(:) = utemp(jgmin:jgmax)*density

            !Differentiate w.r.t time
            do ib = bmin(1),bmax(1)
            do kb = bmin(3),bmax(3)
        	    CFD_u_cnst(ib,bmin(2):bmax(2),kb,1) =  u_t(:)
	            CFD_Pi_dS( ib,bmin(2):bmax(2),kb,1) = (u_t(:) - u_mdt(:))/delta_t
            enddo
            enddo

!            if (mod(iter,tplot*Nvel_ave) .eq. 0) then
!                do jb = bmin(2),bmax(2)
!                    write(irank*10000+iter,'(5i6,4f18.12)'),iter, iblock, jblock, kblock, jb, u_t(jb), & 
!                                             minval(CFD_Pi_dS(bmin(1):bmax(1),jb,bmin(3):bmax(3),1)), & 
!                            CFD_Pi_dS(ceiling(0.5d0*(bmin(1)+bmax(1))),jb,ceiling(0.5d0*(bmin(3)+bmax(3))),1), & 
!                                             maxval(CFD_Pi_dS(:,jb,:,1))
!                enddo
!            endif
!            close(irank*10000+iter,status='keep')

!        else

!            if (mod(iter,tplot*Nvel_ave) .eq. 0) then
!                do jb = bmin(2),bmax(2)
!                    write(irank*10000+iter,'(5i6,4f18.12)'),iter, iblock, jblock, kblock, jb, 0.d0,0.d0,0.d0,0.d0
!            enddo
!        endif
!        close(irank*10000+iter,status='keep')
        endif

	case(6)

		stop "NO case 6 for CV test forces"


    end select

end subroutine get_test_forces


!sigmoid type function to drive velocity to required value
subroutine add_velocity_correction_forces(flag,igmin,jgmin,kgmin,igmax,jgmax,kgmax)
	use arrays_MD, only : r, v
	implicit none

	!Perhaps they need a good talking to, if you don't mind my saying so. 
	!Perhaps a bit more. My girls, sir, they didn't care for the Overlook 
	!at first. One of them actually stole a pack of matches, and tried to 
	!burn it down. 
	!But I "corrected" them sir. 
	!And when my wife tried to prevent me from doing my duty,
	!I "corrected" her.

	integer,intent(in)	:: flag,igmin,jgmin,kgmin,igmax,jgmax,kgmax

	!While currently being "corrected", don't apply further constraints
	logical,allocatable,dimension(:,:,:),save		:: correction_lock

	integer											:: n,i,j,k,t, bin(3)
	integer,parameter								:: t_correct=10
	integer,allocatable,dimension(:,:,:),save		:: correctstart
	!Machine Precision -- ε such that 1 = 1+ε on a computer
	double precision,parameter						:: tol = 2.2204460E-16 
	double precision,allocatable,dimension(:,:,:,:) :: u_CV
	double precision,allocatable,dimension(:,:,:,:),save :: u_error

	if (iter .lt. CVforce_starttime .or. apply_CVforce .eq. .false.) return

	if (iter .eq. CVforce_starttime .and. apply_CVforce .eq. .true.) then
		allocate(F_correct(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
		allocate(correctstart(nbins(1)+2,nbins(2)+2,nbins(3)+2))
		allocate(correction_lock(nbins(1)+2,nbins(2)+2,nbins(3)+2))
		allocate(u_error(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
		correction_lock=.false.
	endif

	allocate(u_CV(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))

	!Get velocity of current CV cells
	u_CV = 0.d0
	do n = 1,np
		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1
		!Add molecule to overlap list
		u_CV( bin(1),bin(2),bin(3),:) = u_CV( bin(1),bin(2),bin(3),:) + v(:,n)
	enddo

	!Difference in velocity
	!u_correct = (u_CV - CFD_u_cnst)*tanh(t/t_correct)

	!Result force du_correct/dt
	F_correct = 0.d0
	do i = 2,nbins(1)+1
	do j = 2,nbins(2)+1
	do k = 2,nbins(3)+1
		!Check if CV is currently being corrected
		if (.not. correction_lock(i,j,k)) then
			u_error(i,j,k,:) = u_CV(i,j,k,:) - CFD_u_cnst(i,j,k,:)
			if (any(u_error(i,j,k,:) .gt. tol*1000.d0)) then
				correction_lock(i,j,k) = .true.
				correctstart(i,j,k) = iter
				IF (i .eq. 3 .and. j .eq. 3 .and. k .eq. 3) then
					print'(a,4i5,7f10.5)', 'New error detected',iter,i,j,k,u_CV(i,j,k,:),CFD_u_cnst(i,j,k,:),tol*1000.d0
				endif
			endif
		endif

		!If correction lock is on, calculate correction
		if (correction_lock(i,j,k)) then
			!Get time passed relative to start of correction
			t = (iter - correctstart(i,j,k))
			F_correct(i,j,k,:) = u_error(i,j,k,:)/(t_correct*cosh(dble(t)/dble(t_correct))**2.d0)
			!Switch lock off if end of correction period
			if (t .eq. t_correct) correction_lock(i,j,k) = .false.
			!print'(a,6i5,6f10.5)', 'Correction force',iter,t,correctstart(i,j,k),i,j,k,u_error(i,j,k,:),F_correct(i,j,k,:)
		endif
	enddo
	enddo
	enddo

end subroutine add_velocity_correction_forces


!=============================================================================
! Average molecules in overlap region to obtain values for 
! constrained dynamics algorithms
! AND iterate until constraint is consistent
!-----------------------------------------------------------------------------

subroutine average_over_allbins_iter(flag,igmin,jgmin,kgmin,igmax,jgmax,kgmax)
	use computational_constants_MD, only : nhb, iblock,CVforce_testcaseflag
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	use physical_constants_MD, only : pi
	use CV_objects, only :  CV_constraint, CVE => CVcheck_energy
	use cumulative_momentum_flux_mod, only : cumulative_momentum_flux
	use cumulative_energy_flux_mod, only : cumulative_energy_flux
    use messenger, only : globalise_bin, localise_bin
	implicit none

	integer,intent(in)	                            :: flag,igmin,jgmin,kgmin,igmax,jgmax,kgmax

	integer         	                            :: bmin(3),bmax(3)
	integer											:: ib,jb,kb, maxiter,global_converged, converge_cells
	integer											:: n,i,molno,cellnp, attempt
	integer,dimension(3)							:: bin
	integer,allocatable,dimension(:,:,:) 			:: box_np_bfr,box_np_afr

	double precision								:: convergence,w(3),wmol(3)
	double precision,parameter						:: tol = 2.2204460E-16 !Machine Precision -- ε such that 1 = 1+ε on a computer
	double precision,dimension(:,:),allocatable 	:: r_temp,v_temp,converge_history
	double precision,allocatable,dimension(:,:,:,:) :: u_bin_afr,u_bin_bfr

	real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: MD_Pi_dS,MD_rhouu_dS,MD_rhouu_dS_prev

	type(node), pointer								:: old, current

	integer,save									:: convergence_count = 0
	logical											:: converged
	logical,save									:: apply_force = .false., first_time = .true.
	double precision, allocatable, dimension(:,:,:,:)   :: wsum

    bmin = max(localise_bin((/igmin,jgmin,kgmin/)),(/2,2,2 /))
    bmax = min(localise_bin((/igmax,jgmax,kgmax/)),nbins(:)+1)

	!delta_t = 0.005

	!Allocate and zero arrays
	allocate(r_temp(3,np),v_temp(3,np))
	r_temp = 0.d0; v_temp = 0.d0
	allocate(box_np_bfr(nbins(1)+2,nbins(2)+2,nbins(3)+2))
	allocate(box_np_afr(nbins(1)+2,nbins(2)+2,nbins(3)+2))
	allocate(u_bin_bfr(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(u_bin_afr(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(MD_Pi_dS(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(MD_rhouu_dS(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(MD_rhouu_dS_prev(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(F_constraint(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(wsum(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); 	wsum = 0.d0

	!Get velocity/positions before force applied
	box_np_bfr = 0; u_bin_bfr = 0.d0
	do n = 1,np
		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1
		!Add molecule to overlap list
		box_np_bfr(bin(1),bin(2),bin(3))  = box_np_bfr(bin(1),bin(2),bin(3))   + 1
		u_bin_bfr( bin(1),bin(2),bin(3),:)= u_bin_bfr( bin(1),bin(2),bin(3),:) + v(:,n)

        !Get weighting tally
        wsum(bin(1),bin(2),bin(3),:) = wsum(bin(1),bin(2),bin(3),:) + F_dist(:,n)

		if (iter .eq. 1500 .and. irank .eq. 1) then
		    write(unit=200,fmt='(3i6,6f12.7)'),bin,r(:,n), F_dist(:,n)
		endif

	enddo

	!TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP
	!Write out for debugging purposes (this must be removed at some point)
	!do n=1,np
	!	bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1
	!	if (bin(1) .eq. 2 .and. bin(2) .eq. 2 .and. bin(3) .eq. 2 ) then
	!		write(5000,'(a,2i6,3f10.5,6f14.5)'), 'Weighting functions', iter,n, r(:,n), F_dist(:,n), wsum(bin(1),bin(2),bin(3),:)
	!	endif
!		!Set to zero if sum of values are zero
!		if (abs(wsum(bin(1),bin(2),bin(3),1) .lt. 0.0001d0)) then
!			F_dist(1,n) = 1.d0; wsum(bin(1),bin(2),bin(3),1) = wsum(bin(1),bin(2),bin(3),1) + F_dist(1,n)
!			stop "Wsum x is too small"
!		endif
!		if (abs(wsum(bin(1),bin(2),bin(3),2) .lt. 0.0001d0)) then
!			F_dist(2,n) = 1.d0; wsum(bin(1),bin(2),bin(3),2) = wsum(bin(1),bin(2),bin(3),2) + F_dist(2,n)
!			stop "Wsum y is too small"
!		endif
!		if (abs(wsum(bin(1),bin(2),bin(3),3) .lt. 0.0001d0)) then
!			F_dist(3,n) = 1.d0; wsum(bin(1),bin(2),bin(3),3) = wsum(bin(1),bin(2),bin(3),3) + F_dist(3,n)
!			stop "Wsum z is too small"
!		endif
	!enddo
	!TEMP TEMP TEMP TEMP TEMP TEMP TEMP TEMP

	!if (any(abs(wsum(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)) .lt. 0.0001d0)) then	
	!	print*, wsum(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:), F_dist
	!	stop "Wsum is too small"
	!endif
	
	!Get velocity/positions at next timestep without constraint
	box_np_afr = 0; u_bin_afr = 0.d0
	do n = 1,np

		v_temp(:,n) = v(:,n) + delta_t * a(:,n) 		!Velocity calculated from acceleration
		r_temp(:,n) = r(:,n) + delta_t * v_temp(:,n)	!Position calculated from velocity
		bin(:) = ceiling((r_temp(:,n)+0.5d0*domain(:))/Fbinsize(:))+1

		!Add molecule to overlap list
		box_np_afr(bin(1),bin(2),bin(3)) =  box_np_afr(bin(1),bin(2),bin(3))  + 1
		u_bin_afr(bin(1),bin(2),bin(3),:) = u_bin_afr(bin(1),bin(2),bin(3),:) + v_temp(:,n)
	enddo

	! - - - - -  Get CV Momentum Totals - - - - - - - - - 
	!Total surface stresses
	MD_Pi_dS(:,:,:,:)=	0.25d0 *((CV_constraint%Pxy(:,:,:,:,1)-CV_constraint%Pxy(:,:,:,:,4)) &
	              				+(CV_constraint%Pxy(:,:,:,:,2)-CV_constraint%Pxy(:,:,:,:,5)) &
		    	              	+(CV_constraint%Pxy(:,:,:,:,3)-CV_constraint%Pxy(:,:,:,:,6)))

	!Update momentum flux using velocity at next timestep
	CV_constraint%flux = 0.d0
	call cumulative_momentum_flux(r_temp,v_temp,CV_constraint%flux)

	!Total CV flux
	MD_rhouu_dS(:,:,:,:)=	((CV_constraint%flux(:,:,:,:,1)-CV_constraint%flux(:,:,:,:,4)) &
	          	 			+(CV_constraint%flux(:,:,:,:,2)-CV_constraint%flux(:,:,:,:,5)) &
	          	 			+(CV_constraint%flux(:,:,:,:,3)-CV_constraint%flux(:,:,:,:,6)))/delta_t

	F_constraint = 0.d0; converged = .false.; maxiter = 100; 
	!allocate(converge_history(maxiter,3))
	!Iterate until momentum flux and applied CV force are consistent
	do attempt = 1,maxiter

    	do ib = bmin(1),bmax(1)
    	do jb = bmin(2),bmax(2)
    	do kb = bmin(3),bmax(3)

			if (box_np_afr(ib,jb,kb) .eq. 0) then
				print'(a,3i6)', "Error -- no molecules left in box after timestep", ib,jb,kb
				stop
			endif

        	!if (box_np_bfr(ib,jb,kb) .ne. 0 .and. apply_CVforce .and. &
			!	((ib .eq. 3 .and. jb .eq. 3 .and. kb .eq. 3) .or. & 
			!	 (ib .eq. 3 .and. jb .eq. 4 .and. kb .eq. 3))) then

        	if (box_np_bfr(ib,jb,kb) .ne. 0 .and. apply_CVforce) then

        		select case (CVforce_testcaseflag)
        		case(1)
        			!Apply differential CV flux constraint
        			F_constraint(ib,jb,kb,:)  =	       MD_Pi_dS(ib,jb,kb,:) + MD_rhouu_dS(ib,jb,kb,:) &
    										  - ( -CFD_rhouu_dS(ib,jb,kb,:) + CFD_Pi_dS(ib,jb,kb,:) )*volume


                    !if (ib .eq. 3 .and. jb .eq. 3 .and. kb .eq. 3) then
                    !    print'(a,i5,2(a,3f18.10))','iter = ', iter, ' L(u)= ', F_constraint(ib,jb,kb,:), ' L(û)= ', wsum(bin(1),bin(2),bin(3),:)
                   ! endif

        		case default
    				stop "Error - Only full constraint possbile with average_over_allbins_iter"
        		end select

        		!Divide constraint per molecules
        		!F_constraint(ib,jb,kb,:) = F_constraint(ib,jb,kb,:)/dble(box_np_bfr(ib,jb,kb))

!				if (ib .eq. 2 .and. jb .eq. 7 .and. kb .eq. 2 .or. &
!					ib .eq. 2 .and. jb .eq. 7 .and. kb .eq. 3 .or. &
!					ib .eq. 3 .and. jb .eq. 7 .and. kb .eq. 3 ) then
! 				if (any(abs(F_constraint(ib,jb,kb,:)) .ne. 0.000000)) then
!					print'(a,9f16.8,i4,3f16.8)', 'Forces ', MD_Pi_dS(ib,jb,kb,:),MD_rhouu_dS(ib,jb,kb,:), F_constraint(ib,jb,kb,:),box_np_afr(ib,jb,kb),u_bin_afr(ib,jb,kb,:)
! 					print'(3i8,3f10.6)', ib,jb,kb,F_constraint(ib,jb,kb,:)
! 				endif

        	else
        		F_constraint(ib,jb,kb,:) = 0.d0
        	endif

    	enddo
    	enddo
    	enddo

    	!Get velocity at next timestep with constraint
    	box_np_afr = 0; u_bin_afr = 0.d0
    	do n = 1,np
    		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1

    		!Velocity calculated from acceleration
            !wmol(:) = F_dist(:,n)/wsum(bin(1),bin(2),bin(3),:)
			!write(5000,'(a,2i6,3f10.5,6f14.5)'), 'Weighting functions', iter,n, r(:,n), F_dist(:,n), F_constraint(bin(1),bin(2),bin(3),:)

    		!v_temp(:,n) = v(:,n) + delta_t * (a(:,n) - F_constraint(bin(1),bin(2),bin(3),:)*wmol(:))
		
    		v_temp(:,n) = v(:,n) + delta_t * (a(:,n) - (F_constraint(bin(1),bin(2),bin(3),:)+F_dist(:,n))/dble(box_np_bfr(bin(1),bin(2),bin(3))))

    		!Position calculated from velocity
    		r_temp(:,n) = r(:,n) + delta_t * v_temp(:,n)

    		bin(:) = ceiling((r_temp(:,n)+0.5d0*domain(:))/Fbinsize(:))+1
    		box_np_afr(bin(1),bin(2),bin(3))  = box_np_afr(bin(1),bin(2),bin(3))   + 1
    		u_bin_afr( bin(1),bin(2),bin(3),:)= u_bin_afr( bin(1),bin(2),bin(3),:) + v_temp(:,n)
    	enddo

! 		print'(i6,3i3,a,3i6,6f12.6)', iter,3,3,3,' Attempt=', attempt,box_np_bfr(ib,jb,kb), & 
! 			box_np_afr(ib,jb,kb),sum(MD_rhouu_dS(3,3,3,:)/(product(Fbinsize))), & 
! 				sum(MD_rhouu_dS_prev(3,3,3,:)/(product(Fbinsize))) , sum(MD_Pi_dS), &
! 				sum(F_constraint(3,3,3,:))/product(Fbinsize), & 
! 				sum(CFD_rhouu_dS(3,3,3,:)),sum(CFD_Pi_dS(3,3,3,:))

		!print'(i4,l,12f10.5)', attempt, converged,MD_rhouu_dS_prev(3,3,3,:), MD_rhouu_dS(3,3,3,:),MD_rhouu_dS_prev(3,4,3,:), MD_rhouu_dS(3,4,3,:)

		!Check convergence
 		if (converged) then
			exit
		elseif (attempt .eq. maxiter-1) then
			print*, "Warning -- CV Force could not converge"
			continue
			!stop "Error -- CV Force could not converge"
		else
			!Keep going
		endif

		!Update momentum flux using velocity at next timestep
		CV_constraint%flux = 0.d0
		call cumulative_momentum_flux(r_temp,v_temp,CV_constraint%flux)
		! Exchange all halo values for CV values ready to apply forces
		call CV_constraint%swap_halos(nbinso)

    	! - - - - -  Get CV Momentum Totals - - - - - - - - - 
    	!Total CV flux
		MD_rhouu_dS_prev = MD_rhouu_dS
    	MD_rhouu_dS(:,:,:,:)=((CV_constraint%flux(:,:,:,:,1)-CV_constraint%flux(:,:,:,:,4)) &
    	          	 		 +(CV_constraint%flux(:,:,:,:,2)-CV_constraint%flux(:,:,:,:,5)) &
    	          	 		 +(CV_constraint%flux(:,:,:,:,3)-CV_constraint%flux(:,:,:,:,6)))/delta_t

!		converge_history(attempt,:) = sum(sum(sum(abs(MD_rhouu_dS_prev(2:nbins(1)+1,   &
!                                                                       2:nbins(2)+1,   &
! 																	   2:nbins(3)+1,:) &
!													  	-  MD_rhouu_dS(2:nbins(1)+1,   &
!                                                                       2:nbins(2)+1,   &
! 																       2:nbins(3)+1,:)),1),1),1)

		!If fluxes have converged (exactly) then exit after calculating final force
		!if (all(MD_rhouu_dS_prev .eq. MD_rhouu_dS)) then
		if (all(abs(MD_rhouu_dS_prev-MD_rhouu_dS) .lt. tol)) then
			converged = .true.
			!Wait for all processors to converge
			global_converged = 1
			call globalMinInt(global_converged)
			if (global_converged .ne. 1) converged = .false.
            convergence = 0.d0; converge_cells = 0
	        do ib = bmin(1),bmax(1)
	        do jb = bmin(2),bmax(2)
	        do kb = bmin(3),bmax(3)
			    convergence = convergence + sum(abs(MD_rhouu_dS_prev(ib,jb,kb,:) - MD_rhouu_dS(ib,jb,kb,:)))
                converge_cells = converge_cells + 1
		    enddo
		    enddo
		    enddo
			!print'(a,2i5,2i8,f28.17)', 'convergence finished', iter, attempt, convergence_count, converge_cells, convergence

		else
			!Tell other processors "I'm not converged"
			global_converged = 0
			call globalMinInt(global_converged)

			!Ouput convergence if it looks like it's going to be a problem
			if (attempt .gt. 25) then
				convergence = 0.d0; converge_cells = 0
				!if (mod(attempt,20)) delta_t = delta_t/2.d0

    	        do ib = bmin(1),bmax(1)
    	        do jb = bmin(2),bmax(2)
    	        do kb = bmin(3),bmax(3)
					if (sum(MD_rhouu_dS_prev(ib,jb,kb,:)) .ne. sum(MD_rhouu_dS(ib,jb,kb,:))) then
						convergence = convergence + sum(abs(MD_rhouu_dS_prev(ib,jb,kb,:) - MD_rhouu_dS(ib,jb,kb,:)))
                        converge_cells = converge_cells + 1
						!print'(5i6,7f14.6)', iter, attempt, ib, jb, kb, sum(abs(MD_rhouu_dS_prev(ib,jb,kb,:) - MD_rhouu_dS(ib,jb,kb,:))), MD_rhouu_dS_prev(ib,jb,kb,:),MD_rhouu_dS(ib,jb,kb,:)
					endif
				enddo
				enddo
				enddo
				convergence_count = convergence_count + 1
				print'(a,2i5,2i8,f28.17)', 'convergence ', iter, attempt, convergence_count, converge_cells, convergence
                !print*, converge_history 
			endif
		endif

	enddo

	deallocate(r_temp,v_temp,box_np_afr,box_np_bfr,u_bin_afr,u_bin_bfr)

end subroutine average_over_allbins_iter




subroutine average_over_allbins_iter_optimised(flag,igmin,jgmin,kgmin,igmax,jgmax,kgmax)
	use computational_constants_MD, only : nhb, iblock,CVforce_testcaseflag
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	use physical_constants_MD, only : pi
	use CV_objects, only :  CV_constraint, CVE => CVcheck_energy
	use cumulative_momentum_flux_mod, only : cumulative_momentum_flux
	use cumulative_energy_flux_mod, only : cumulative_energy_flux
    use messenger, only : globalise_bin, localise_bin
	implicit none

	integer,intent(in)	                            :: flag,igmin,jgmin,kgmin,igmax,jgmax,kgmax

	integer         	                            :: bmin(3),bmax(3)
	integer											:: ib,jb,kb, maxiter,global_converged, converge_cells
	integer											:: n,i,molno,cellnp, attempt
	integer,dimension(3)							:: bin
	integer,allocatable,dimension(:,:,:) 			:: box_np_bfr,box_np_afr

	double precision								:: convergence,w(3),wmol(3)
	double precision,parameter						:: tol = 2.2204460E-16 !Machine Precision -- ε such that 1 = 1+ε on a computer
	double precision,dimension(:,:),allocatable 	:: r_temp,v_temp,converge_history
	double precision,allocatable,dimension(:,:,:,:) :: u_bin_afr,u_bin_bfr

	real(kind(0.d0)),allocatable,dimension(:,:,:,:) :: MD_Pi_dS,MD_rhouu_dS,MD_rhouu_dS_prev

	type(node), pointer								:: old, current

	integer,save									:: convergence_count = 0
	logical											:: converged
	logical,save									:: apply_force = .false., first_time = .true.
	double precision, allocatable, dimension(:,:,:,:)   :: wsum

	if (CVforce_testcaseflag .ne. 1) then
		stop "Error - Only full constraint possbile with average_over_allbins_iter"
	endif

    bmin = max(localise_bin((/igmin,jgmin,kgmin/)),(/2,2,2 /))
    bmax = min(localise_bin((/igmax,jgmax,kgmax/)),nbins(:)+1)

	!Allocate and zero arrays
	allocate(r_temp(3,np),v_temp(3,np))
	r_temp = 0.d0; v_temp = 0.d0
	allocate(box_np_bfr(nbins(1)+2,nbins(2)+2,nbins(3)+2))
	allocate(box_np_afr(nbins(1)+2,nbins(2)+2,nbins(3)+2))
	allocate(u_bin_afr(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(MD_Pi_dS(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(MD_rhouu_dS(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(MD_rhouu_dS_prev(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(F_constraint(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))

	!Get velocity/positions before force applied
	box_np_bfr = 0; u_bin_bfr = 0.d0
	do n = 1,np
		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1
		!Add molecule to overlap list
		box_np_bfr(bin(1),bin(2),bin(3))  = box_np_bfr(bin(1),bin(2),bin(3))   + 1
	enddo

	!Get velocity/positions at next timestep without constraint
	box_np_afr = 0; u_bin_afr = 0.d0
	do n = 1,np
		v_temp(:,n) = v(:,n) + delta_t * a(:,n) 		!Velocity calculated from acceleration
		r_temp(:,n) = r(:,n) + delta_t * v_temp(:,n)	!Position calculated from velocity
	enddo
	!Update momentum flux using velocity at next timestep
	CV_constraint%flux = 0.d0
	call cumulative_momentum_flux(r_temp,v_temp,CV_constraint%flux)

	! - - - - -  Get CV Momentum Totals - - - - - - - - - 
	!Total surface stresses
	MD_Pi_dS(:,:,:,:)=	0.25d0 *((CV_constraint%Pxy(:,:,:,:,1)-CV_constraint%Pxy(:,:,:,:,4)) &
	              				+(CV_constraint%Pxy(:,:,:,:,2)-CV_constraint%Pxy(:,:,:,:,5)) &
		    	              	+(CV_constraint%Pxy(:,:,:,:,3)-CV_constraint%Pxy(:,:,:,:,6)))
	!Total CV flux
	MD_rhouu_dS(:,:,:,:)=  -((CV_constraint%flux(:,:,:,:,1)-CV_constraint%flux(:,:,:,:,4)) &
	          	 			+(CV_constraint%flux(:,:,:,:,2)-CV_constraint%flux(:,:,:,:,5)) &
	          	 			+(CV_constraint%flux(:,:,:,:,3)-CV_constraint%flux(:,:,:,:,6)))/delta_t

	! Stress difference between MD and CFD
	! F_stress_CV(ib,jb,kb,:)  = MD_Pi_dS(ib,jb,kb,:) - CFD_Pi_dS(ib,jb,kb,:)*volume
	!	F_stress_CV(ib,jb,kb,:) => F_stress_mol(n,:) 

	! CFD fluxes are constant
	! F_CFDflux_CV(ib,jb,kb,:) = CFD_rhouu_dS(ib,jb,kb,:)*volume
		
	! MD fluxes are iterated
	! F_MDflux_CV(ib,jb,kb,:)  = MD_rhouu_dS(ib,jb,kb,:)

	!Iterate until momentum flux and applied CV force are consistent
	F_constraint = 0.d0; converged = .false.; maxiter = 100; 
	do attempt = 1,maxiter

    	do ib = bmin(1),bmax(1)
    	do jb = bmin(2),bmax(2)
    	do kb = bmin(3),bmax(3)

        	if (box_np_bfr(ib,jb,kb) .ne. 0 .and. apply_CVforce) then
       			!Apply differential CV flux constraint
       			F_constraint(ib,jb,kb,:)  =	   -MD_rhouu_dS(ib,jb,kb,:) +  MD_Pi_dS(ib,jb,kb,:)  &
   										  - ( -CFD_rhouu_dS(ib,jb,kb,:) + CFD_Pi_dS(ib,jb,kb,:) )*volume
        	else
        		F_constraint(ib,jb,kb,:) = 0.d0
        	endif

    	enddo
    	enddo
    	enddo

    	!Get velocity at next timestep with constraint
    	box_np_afr = 0; u_bin_afr = 0.d0
    	do n = 1,np
    		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1

    		!Velocity calculated from acceleration
    		v_temp(:,n) = v(:,n) + delta_t * (a(:,n) - (F_constraint(bin(1),bin(2),bin(3),:)+F_dist(:,n))/dble(box_np_bfr(bin(1),bin(2),bin(3))))

    		!Position calculated from velocity
    		r_temp(:,n) = r(:,n) + delta_t * v_temp(:,n)

    		bin(:) = ceiling((r_temp(:,n)+0.5d0*domain(:))/Fbinsize(:))+1
    		box_np_afr(bin(1),bin(2),bin(3))  = box_np_afr(bin(1),bin(2),bin(3))   + 1
    		u_bin_afr( bin(1),bin(2),bin(3),:)= u_bin_afr( bin(1),bin(2),bin(3),:) + v_temp(:,n)
    	enddo

		if (any(box_np_afr(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3)) .eq. 0)) then
			print'(a,3i6)', "Error -- no molecules left in box after timestep", ib,jb,kb
			stop
		endif

		!Check convergence
 		if (converged) then
			exit
		elseif (attempt .eq. maxiter-1) then
			print*, "Warning -- CV Force could not converge"
			continue
			!stop "Error -- CV Force could not converge"
		else
			!Keep going
		endif

		!Update momentum flux using velocity at next timestep
		CV_constraint%flux = 0.d0
		call cumulative_momentum_flux(r_temp,v_temp,CV_constraint%flux)
		! Exchange all halo values for CV values ready to apply forces
		call CV_constraint%swap_halos(nbinso)

    	! - - - - -  Get CV Momentum Totals - - - - - - - - - 
    	!Total CV flux
		MD_rhouu_dS_prev = MD_rhouu_dS
    	MD_rhouu_dS(:,:,:,:)=-((CV_constraint%flux(:,:,:,:,1)-CV_constraint%flux(:,:,:,:,4)) &
    	          	 		  +(CV_constraint%flux(:,:,:,:,2)-CV_constraint%flux(:,:,:,:,5)) &
    	          	 		  +(CV_constraint%flux(:,:,:,:,3)-CV_constraint%flux(:,:,:,:,6)))/delta_t

		!If fluxes have converged (to tolerence) then exit after calculating final force
		if (all(abs(MD_rhouu_dS_prev-MD_rhouu_dS) .lt. tol)) then
			converged = .true.
			!Wait for all processors to converge
			global_converged = 1
			call globalMinInt(global_converged)
			if (global_converged .ne. 1) then
				converged = .false.
			else
!	            convergence = 0.d0; converge_cells = 0
!		        do ib = bmin(1),bmax(1)
!		        do jb = bmin(2),bmax(2)
!		        do kb = bmin(3),bmax(3)
!				    convergence = convergence + sum(abs(MD_rhouu_dS_prev(ib,jb,kb,:) - MD_rhouu_dS(ib,jb,kb,:)))
!    	            converge_cells = converge_cells + 1
!			    enddo
!			    enddo
!			    enddo
!				print'(a,2i5,2i8,f28.17)', 'converged ', iter, attempt, convergence_count, converge_cells, convergence
			endif

		else
			!Tell other processors "I'm not converged"
			global_converged = 0
			call globalMinInt(global_converged)

			!Ouput convergence if it looks like it's going to be a problem
			if (attempt .gt. 25) then
				convergence = 0.d0; converge_cells = 0
				!if (mod(attempt,20)) delta_t = delta_t/2.d0

    	        do ib = bmin(1),bmax(1)
    	        do jb = bmin(2),bmax(2)
    	        do kb = bmin(3),bmax(3)
					if (sum(MD_rhouu_dS_prev(ib,jb,kb,:)) .ne. sum(MD_rhouu_dS(ib,jb,kb,:))) then
						convergence = convergence + sum(abs(MD_rhouu_dS_prev(ib,jb,kb,:) - MD_rhouu_dS(ib,jb,kb,:)))
                        converge_cells = converge_cells + 1
					endif
				enddo
				enddo
				enddo
				convergence_count = convergence_count + 1
				print'(a,2i5,2i8,f28.17)', 'convergence ', iter, attempt, convergence_count, converge_cells, convergence
			endif
		endif

	enddo

	deallocate(r_temp,v_temp,box_np_afr,box_np_bfr,u_bin_afr,u_bin_bfr)

end subroutine average_over_allbins_iter_optimised




!=============================================================================
! Apply force to molecules in overlap region
!-----------------------------------------------------------------------------

subroutine apply_force
	use arrays_MD, only : r,v,a
	use physical_constants_MD, only : density
	use computational_constants_MD, only : eflux_outflag
	use module_record_external_forces
	implicit none

	integer								:: n, molno
	integer,dimension(3)				:: bin
	double precision,dimension(3)		:: velvect

	!Loop over all molecules and apply constraint
	do n = 1, np
		!Get bin
		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1
		!Apply force
		a(:,n) = a(:,n) - F_constraint(bin(1),bin(2),bin(3),:)
		!Add external force to CV total
		if (eflux_outflag .eq. 4) then
			velvect(:) = v(:,n) + 0.5d0*delta_t*a(:,n)
			call record_external_forces(F_constraint(bin(1),bin(2),bin(3),:),r(:,n),velvect)
		else
			call record_external_forces(F_constraint(bin(1),bin(2),bin(3),:),r(:,n))
		endif

	enddo

end subroutine apply_force



!=============================================================================
! Apply force to molecules in overlap region
!-----------------------------------------------------------------------------

subroutine apply_force_weighted
	use arrays_MD, only : r,v,a
	use physical_constants_MD, only : density
	use computational_constants_MD, only : eflux_outflag
	use module_record_external_forces
	implicit none

	integer									:: n, molno
	integer,dimension(3)					:: bin
	integer, allocatable, dimension(:,:,:)	:: boxnp
	double precision	:: rand
	double precision,dimension(3)			:: w,wmol,velvect,Ftest,Fwtest

	double precision, allocatable, dimension(:,:,:,:)   :: wsum

	allocate(boxnp(nbins(1)+2,nbins(2)+2,nbins(3)+2)); 	boxnp = 0
	allocate(wsum(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)); wsum = 0.d0

    !Get weighting function -- linear across domain from -0.5 to 0.5
	do n = 1,np

!		!Get bin
		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1

!		!Add molecule weighting to sum
		boxnp(bin(1),bin(2),bin(3))  = boxnp(bin(1),bin(2),bin(3)) + 1
!		!wsum( bin(1),bin(2),bin(3),:)= wsum( bin(1),bin(2),bin(3),:) + F_dist(:,n) 

	enddo

	!Loop over all molecules and apply constraint
	Ftest = 0.d0; Fwtest = 0.d0
	do n = 1, np

		!Get bin
		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))+1

		!Apply force
        !wmol(:) = F_dist(:,n)/wsum(bin(1),bin(2),bin(3),:)
		!a(:,n) = a(:,n) - F_constraint(bin(1),bin(2),bin(3),:) * wmol(:)
		a(:,n) = a(:,n) - (F_constraint(bin(1),bin(2),bin(3),:)+F_dist(:,n))/dble(boxnp(bin(1),bin(2),bin(3)))


		call random_number(rand)
		if (iter .eq. 150 ) then
			if ( bin(1) .eq. 3 .and.  bin(3) .eq. 3) then
			if (mod(bin(2),2) .eq. 0) then
				write(5000,'(a,2i6,3f10.5,6f14.5)'), 'Weight', iter,n, r(:,n), F_dist(:,n), (F_constraint(bin(1),bin(2),bin(3),:)/dble(boxnp(bin(1),bin(2),bin(3))))
			else
				write(5001,'(a,2i6,3f10.5,6f14.5)'), 'Weight', iter,n, r(:,n), F_dist(:,n), (F_constraint(bin(1),bin(2),bin(3),:)/dble(boxnp(bin(1),bin(2),bin(3))))
			endif
			endif
		endif

		!Add external force to CV total
		if (eflux_outflag .eq. 4) then
			velvect(:) = v(:,n) + 0.5d0*delta_t*a(:,n)
			call record_external_forces((F_constraint(bin(1),bin(2),bin(3),:)+F_dist(:,n))/dble(boxnp(bin(1),bin(2),bin(3))),r(:,n),velvect)
		else
			call record_external_forces((F_constraint(bin(1),bin(2),bin(3),:)+F_dist(:,n))/dble(boxnp(bin(1),bin(2),bin(3))),r(:,n))
		endif

	enddo

end subroutine apply_force_weighted


end subroutine apply_CV_force_multibin



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




! ===================================================================
! I think multi-bin replaces all the code below this point using 
! the special case where only one bin specified...
! ===================================================================



!!=========================================================================
!! Routine to test flux calculations by various methods

!subroutine apply_CV_force(iter)
!	use computational_constants_MD, only : irank, jblock, npy, globaldomain, CVforce_flag, & 
!										   CVforce_starttime, VOID,  Nsteps
!	use calculated_properties_MD, only : pressure, nbins
!	use librarymod, only : get_new_fileunit
!	use module_external_forces, only : np, momentum_flux, irank,nbins, nbinso,domain,delta_t,Nvflux_ave
!	use CV_objects, only : 	CV_constraint
!	implicit none
!	
!	integer,intent(in)				:: iter

!	logical							:: apply_CVforce = .false.
!	integer							:: ibin,jbin,kbin
!	integer							:: box_np_bfr,box_np_afr
!	integer,allocatable 			:: list_bfr(:),list_afr(:)

!	real(kind(0.d0))				:: volume
!	real(kind(0.d0)),dimension(3)	:: binsize,MD_Pi_dS,MD_rhouu_dS,F_constraint
!	real(kind(0.d0)),dimension(3)	:: MD_rhouu_dS_afr, MD_rhouu_dS_prev
!	real(kind(0.d0)),dimension(3)	:: u_bin,F_bin,u_bin_afr, u_bin_bfr, F_bin1, F_bin2
!	real(kind(0.d0)),dimension(3)	:: CFD_Pi_dS,CFD_rhouu_dS,CFD_u_cnst

!	if (CVforce_flag .eq. VOID) return

!	!Test case focuses on a single CV
!	binsize = domain/nbins
!	volume = product(binsize)
!	!ibin = 3; jbin = 3; kbin = 3

!	!Exchange CV data ready to apply force
!	call update_CV_halos

!	!Only apply force on top processor
!	if (jblock .ne. npy) return

!	do ibin = 3,3
!	do jbin = 3,3
!	do kbin = 3,3

!		!Get average over current cell and apply constraint forces
!		!call get_continuum_values
!		call get_test_forces(CVforce_flag,ibin,jbin,kbin)

!		!Retreive CV data and calculate force to apply
!		call average_over_bin_iter(CVforce_flag,ibin,jbin,kbin)

!		!Apply the force
!		!call apply_force
!		call apply_force_tests(apply_CVforce)

!		!Set velocity to required value and then apply constraint 
!		!force to adjust its evolution from then on...
!		if (iter .eq. CVforce_starttime .and. apply_CVforce .ne. 0) then
!			!This routine takes the global bin number which is one less than the local 
!			!used in other routines as these include halos
! 			call set_bin_velocity(ibin-1, ibin-1, & 
! 								  jbin-1, jbin-1, & 
! 								  kbin-1, kbin-1, & 
! 								  (/ 0.d0,0.d0,0.d0 /),0)
!		endif

!	enddo
!	enddo
!	enddo

!	!Reset CV force values
!	CV_constraint%flux = 0.d0; CV_constraint%Pxy = 0.d0

!	!Recalculate flux values in line with applied force constraint
!	

!contains

!! Apply arbitary forces for testing purposes
!subroutine get_test_forces(flag,ibin_,jbin_,kbin_)
!	use physical_constants_MD, only : pi
!	implicit none

!	integer,intent(in)	:: flag,ibin_,jbin_,kbin_

!	integer				:: length,fileunit
!	double precision	:: sin_mag, sin_period

!	!Set velocity to required value and then apply constraint 
!	!force to adjust its evolution from then on...
!	if (iter .eq. CVforce_starttime .and. flag .ne. 0) then
!		apply_CVforce = .true.
!	elseif (iter .le. CVforce_starttime) then
!		apply_CVforce = .false.
!	else
!		!Do nothing
!	endif

!	select case(flag)
!	case(0)
!		!No Force applied
!		apply_CVforce = .false.
!	case(1)
!		!Zero Force applied
!		CFD_Pi_dS = 0.d0
!		CFD_rhouu_dS = 0.d0	
!		!Velocity is a constant (assumed zero)
!		CFD_u_cnst = 0.d0
!	case(2)
!		!Constant function 
!		CFD_Pi_dS = 1.d0
!		CFD_rhouu_dS = 0.d0
!		!Velocity is a linear fn of time
!		CFD_u_cnst = iter
!	case(3)
!		CFD_Pi_dS(:) = 0.d0
!		!Sin function shifted by pi/2 so velocity which is given
!		!by the integral (cos) alternates around zero
!		sin_mag = 0.1d0;  sin_period = 100.d0
!		CFD_Pi_dS(1) = sin_mag * cos(2.d0*pi*((iter-CVforce_starttime)/sin_period))! - 0.505*pi)
!		CFD_rhouu_dS = 0.d0

!		! - - - - -  Sinusoid of velocity with NCER special discretisation - - - - - - - - - 
!		CFD_u_cnst = 0.d0
!		CFD_u_cnst(1) = 1.25d0*delta_t*(sin_mag*sin_period/(2.d0*pi))* sin(2.d0*pi*((iter-CVforce_starttime)/sin_period))

!	case(4)
!		! Get continuum values of surface stresses, etc
!		!Read Shear pressure from file...
!		CFD_Pi_dS(:) = 0.d0
!		fileunit = get_new_fileunit()
!		inquire(iolength=length) CFD_Pi_dS(1)
!		open(unit=fileunit,file='./F_hist',form='unformatted', access='direct',recl=length)
!		read(fileunit,rec=iter) CFD_Pi_dS(1)
!		close(fileunit,status='keep')
!		CFD_rhouu_dS = 0.d0
!		! Read velocity too for NCER style proportional constraint
!		fileunit = get_new_fileunit()
!		inquire(iolength=length) CFD_u_cnst(1)
!		open(unit=fileunit,file='./u_hist',form='unformatted', access='direct',recl=length)
!		read(fileunit,rec=iter) CFD_u_cnst(1)
!		close(fileunit,status='keep')
!	end select
!	
!end subroutine get_test_forces

!! Exchange all halo values for CV values ready to apply forces
!subroutine update_CV_halos
!	use CV_objects, only :  CV_constraint
!	use mpi
!	implicit none

!	call CV_constraint%swap_halos(nbinso)

!end subroutine update_CV_halos

!!=============================================================================
!! Average molecules in overlap region to obtain values for 
!! constrained dynamics algorithms
!!-----------------------------------------------------------------------------

!subroutine average_over_bin(flag,ib,jb,kb)
!	use computational_constants_MD, only : nhb, iblock,CVforce_testcaseflag
!	use arrays_MD, only : r, v, a
!	use linked_list, only : node, cell
!	use physical_constants_MD, only : pi
!	use CV_objects, only :  CV_constraint
!	use cumulative_momentum_flux_mod, only : cumulative_momentum_flux
!	implicit none

!	integer,intent(in)		:: flag,ib,jb,kb

!	integer					:: n,i,molno,cellnp
!	integer,dimension(3)	:: bin
!	type(node), pointer		:: old, current
!	!double precision,dimension(:,:),allocatable	:: v_temp
!	double precision,dimension(:,:),allocatable :: r_temp,v_temp
!	logical,save			:: apply_force = .false., first_time = .true.

!	!Get velocity at next timestep without constraint
!	!Zero box averages and add to list array
!	box_np_bfr = 0
!	allocate(list_bfr(np))
!	allocate(r_temp(3,np),v_temp(3,np))
!	r_temp = 0.d0; v_temp = 0.d0
!	do n = 1,np
!		v_temp(:,n) = v(:,n) + delta_t * a(:,n) 		!Velocity calculated from acceleration
!		r_temp(:,n) = r(:,n) + delta_t * v_temp(:,n)	!Position calculated from velocity
!		bin(:) = ceiling((r_temp(:,n)+0.5d0*domain(:))/binsize(:))+1
!		if (bin(1) .eq. ib .and. bin(2) .eq. jb .and. bin(3) .eq. kb) then
!			!Add molecule to overlap list
!			box_np_bfr   =  box_np_bfr + 1
!			list_bfr(box_np_bfr) =  n
!		endif
!	enddo

!	!Find the molecules in constraint region
!	u_bin_bfr(:) = 0.d0
!	F_bin2(:) = 0.d0
!	do i = 1,box_np_bfr
!		n = list_bfr(i)
!		!Save sum of velocity and forces without constraint
!		u_bin_bfr(:) = u_bin_bfr(:) + v_temp(:,n)
!		F_bin2(:)    = F_bin2(:) + a(:,n)
!		!if (iter .eq. 1001) print'(i5,a,3i8,a,i8,3f10.5)', iter, ' bin ',ib,jb,kb,' molno ', n, r(:,n)
!	enddo

!	!Update momentum flux using velocity at next timestep
!	CV_constraint%flux = 0.d0
!	call cumulative_momentum_flux(r_temp,v_temp,CV_constraint%flux)
!	!deallocate(r_temp,v_temp)

!	! - - - - -  Get CV Momentum Totals - - - - - - - - - 
!	!Total CV flux
!	MD_rhouu_dS  =		((CV_constraint%flux(ib,jb,kb,:,1)-CV_constraint%flux(ib,jb,kb,:,4)) &
!	          	 		+(CV_constraint%flux(ib,jb,kb,:,2)-CV_constraint%flux(ib,jb,kb,:,5)) &
!	          	 		+(CV_constraint%flux(ib,jb,kb,:,3)-CV_constraint%flux(ib,jb,kb,:,6)))/delta_t
!	!Total surface stresses
!	MD_Pi_dS  =	0.25d0 *((CV_constraint%Pxy(ib,jb,kb,:,1)-CV_constraint%Pxy(ib,jb,kb,:,4)) &
!	              		+(CV_constraint%Pxy(ib,jb,kb,:,2)-CV_constraint%Pxy(ib,jb,kb,:,5)) &
!	              		+(CV_constraint%Pxy(ib,jb,kb,:,3)-CV_constraint%Pxy(ib,jb,kb,:,6)))

!	if (box_np_bfr .ne. 0 .and. apply_CVforce) then

!		select case (CVforce_testcaseflag)
!		case(1)
!			!Apply differential CV flux constraint
!			F_constraint = MD_Pi_dS+MD_rhouu_dS + (CFD_rhouu_dS-CFD_Pi_dS)*volume
!		!DEBUG CASES WITH TERMS MISSING 
!		case(2)
!			!No molecular convective fluxes
!			F_constraint = MD_Pi_dS + (CFD_rhouu_dS-CFD_Pi_dS)*volume
!		case(3)
!			!CFD fluxes only (no MD convective or stresses)
!			F_constraint = (CFD_rhouu_dS-CFD_Pi_dS)*volume 
!		case(4)
!			!No force at all (debug case)
!			F_constraint = 0.d0
!		!NCER SPECIAL DISCRETISATION
!		case(5)
!			!Apply proportional velocity constraint (Borg et al)
!			F_constraint = - (CFD_u_cnst*dble(box_np_bfr) - u_bin_bfr)/delta_t
!			if (first_time) then 
!				first_time = .false.
!				F_constraint = 0.d0	
!			endif
!		case(6)
!			!NCER special discretisation -- proportional 
!			!velocity constraint including forces
!			F_constraint = MD_Pi_dS - (CFD_u_cnst*dble(box_np_bfr) - u_bin_bfr)/delta_t
!		case(7)
!			!NCER special discretisation with convective fluxes too
!			F_constraint = MD_Pi_dS-MD_rhouu_dS - (CFD_u_cnst*dble(box_np_bfr) - u_bin_bfr)/delta_t
!		end select

!		!Divide constraint per molecules
!		F_constraint = F_constraint/dble(box_np_bfr)

!	else
!		F_constraint = 0.d0
!	endif

!end subroutine average_over_bin


!!=============================================================================
!! Average molecules in overlap region to obtain values for 
!! constrained dynamics algorithms
!! AND iterate until constraint is consistent
!!-----------------------------------------------------------------------------

!subroutine average_over_bin_iter(flag,ib,jb,kb)
!	use computational_constants_MD, only : nhb, iblock,CVforce_testcaseflag
!	use arrays_MD, only : r, v, a
!	use linked_list, only : node, cell
!	use physical_constants_MD, only : pi
!	use CV_objects, only :  CV_constraint, CVE => CVcheck_energy
!	use cumulative_momentum_flux_mod, only : cumulative_momentum_flux
!	use cumulative_energy_flux_mod, only : cumulative_energy_flux
!	implicit none

!	integer,intent(in)		:: flag,ib,jb,kb

!	integer					:: n,i,molno,cellnp, attempt
!	integer,dimension(3)	:: bin
!	type(node), pointer		:: old, current
!	double precision,dimension(:,:),allocatable :: r_temp,v_temp
!	logical					:: converged
!	logical,save			:: apply_force = .false., first_time = .true.

!	!Allocate and zero arrays
!	allocate(list_bfr(np),list_afr(np)); 
!	allocate(r_temp(3,np),v_temp(3,np))
!	r_temp = 0.d0; v_temp = 0.d0

!	!Get velocity/positions before force applied
!	box_np_bfr = 0; u_bin_bfr = 0.d0
!	do n = 1,np
!		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/binsize(:))+1
!		if (bin(1) .eq. ib .and. bin(2) .eq. jb .and. bin(3) .eq. kb) then
!			!Add molecule to overlap list
!			box_np_bfr   =  box_np_bfr + 1
!			list_bfr(box_np_bfr) =  n
!			u_bin_bfr(:) = u_bin_bfr(:) + v_temp(:,n)
!		endif
!	enddo
!	
!	!Get velocity/positions at next timestep without constraint
!	box_np_afr = 0; u_bin_afr = 0.d0
!	do n = 1,np
!		v_temp(:,n) = v(:,n) + delta_t * a(:,n) 		!Velocity calculated from acceleration
!		r_temp(:,n) = r(:,n) + delta_t * v_temp(:,n)	!Position calculated from velocity
!		bin(:) = ceiling((r_temp(:,n)+0.5d0*domain(:))/binsize(:))+1
!		if (bin(1) .eq. ib .and. bin(2) .eq. jb .and. bin(3) .eq. kb) then
!			!Add molecule to overlap list
!			box_np_afr   =  box_np_afr + 1
!			list_afr(box_np_afr) =  n
!			u_bin_afr(:) = u_bin_afr(:) + v_temp(:,n)
!		endif
!	enddo

!	! - - - - -  Get CV Momentum Totals - - - - - - - - - 
!	!Total surface stresses
!	MD_Pi_dS  =	0.25d0 *((CV_constraint%Pxy(ib,jb,kb,:,1)-CV_constraint%Pxy(ib,jb,kb,:,4)) &
!	              		+(CV_constraint%Pxy(ib,jb,kb,:,2)-CV_constraint%Pxy(ib,jb,kb,:,5)) &
!	              		+(CV_constraint%Pxy(ib,jb,kb,:,3)-CV_constraint%Pxy(ib,jb,kb,:,6)))

!	!Update momentum flux using velocity at next timestep
!	CV_constraint%flux = 0.d0
!	call cumulative_momentum_flux(r_temp,v_temp,CV_constraint%flux)

!	!Total CV flux
!	MD_rhouu_dS  =	((CV_constraint%flux(ib,jb,kb,:,1)-CV_constraint%flux(ib,jb,kb,:,4)) &
!	          	 	+(CV_constraint%flux(ib,jb,kb,:,2)-CV_constraint%flux(ib,jb,kb,:,5)) &
!	          	 	+(CV_constraint%flux(ib,jb,kb,:,3)-CV_constraint%flux(ib,jb,kb,:,6)))/delta_t


!	!Iterate until momentum flux and applied CV force are consistent
!	converged = .false.
!	do attempt = 1,100

!    	if (box_np_afr .ne. 0 .and. apply_CVforce) then

!    		select case (CVforce_testcaseflag)
!    		case(1)
!    			!Apply differential CV flux constraint
!    			F_constraint = MD_Pi_dS+MD_rhouu_dS + (CFD_rhouu_dS-CFD_Pi_dS)*volume
!    		!DEBUG CASES WITH TERMS MISSING 
!    		case(2)
!    			!No molecular convective fluxes
!    			F_constraint = MD_Pi_dS + (CFD_rhouu_dS-CFD_Pi_dS)*volume
!				converged = .true.
!    		case(3)
!    			!CFD fluxes only (no MD convective or stresses)
!    			F_constraint = (CFD_rhouu_dS-CFD_Pi_dS)*volume 
!    		case(4)
!    			!No force at all (debug case)
!    			F_constraint = 0.d0
!				converged = .true.
!    		!NCER SPECIAL DISCRETISE
!    		case(5)
!    			!Apply proportional velocity constraint (Borg et al)
!    			F_constraint = - (CFD_u_cnst*dble(box_np_afr) - u_bin_afr)/delta_t
!    			if (first_time) then 
!    				first_time = .false.
!    				F_constraint = 0.d0	
!    			endif
!				converged = .true.
!    		case(6)
!    			!NCER special discretisation -- proportional 
!    			!velocity constraint including forces
!    			F_constraint = MD_Pi_dS - (CFD_u_cnst*dble(box_np_afr) - u_bin_afr)/delta_t
!				converged = .true.
!    		case(7)
!    			!NCER special discretisation with convective fluxes too
!    			F_constraint = MD_Pi_dS+MD_rhouu_dS - (CFD_u_cnst*dble(box_np_afr) - u_bin_afr)/delta_t
!    		end select

!    		!Divide constraint per molecules
!    		F_constraint = F_constraint/dble(box_np_bfr)

!        	!Get velocity at next timestep with constraint
!			box_np_afr = 0; u_bin_afr = 0.d0
!        	do n = 1,np
!        		bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/binsize(:))+1
!        		if (bin(1) .eq. ib .and. bin(2) .eq. jb .and. bin(3) .eq. kb) then
!    				!Velocity calculated from acceleration
!        			v_temp(:,n) = v(:,n) + delta_t * (a(:,n) - F_constraint)
!        		else
!    				!Velocity calculated from acceleration
!        			v_temp(:,n) = v(:,n) + delta_t * a(:,n) 
!        		endif
!    			!Position calculated from velocity
!        		r_temp(:,n) = r(:,n) + delta_t * v_temp(:,n)

!				bin(:) = ceiling((r_temp(:,n)+0.5d0*domain(:))/binsize(:))+1
!        		if (bin(1) .eq. ib .and. bin(2) .eq. jb .and. bin(3) .eq. kb) then
!					box_np_afr   =  box_np_afr + 1
!        			list_afr(box_np_afr) =  n
!        			u_bin_afr(:) = u_bin_afr(:) + v_temp(:,n)
!				endif
!        	enddo

! 			if (attempt .gt. 3) then
!	 			print'(i6,a,3i6,6f12.6)', iter,	' Attempt=', attempt,box_np_bfr,box_np_afr, sum(MD_rhouu_dS/(product(binsize))), & 
!  									sum(MD_rhouu_dS_prev/(product(binsize))) , sum(MD_Pi_dS), &
!  									sum(F_constraint)/product(binsize),sum(CFD_rhouu_dS),sum(CFD_Pi_dS)
! 			endif


!			if (attempt == 99) stop "Error -- CV Force could not converge"


!			if (converged) exit

!        	!Update momentum flux using velocity at next timestep
!        	CV_constraint%flux = 0.d0
!        	call cumulative_momentum_flux(r_temp,v_temp,CV_constraint%flux)

!			!CVE%flux = 0.d0
!			!call cumulative_energy_flux(r_temp,v_temp,CVE%flux)

!        	! - - - - -  Get CV Momentum Totals - - - - - - - - - 
!        	!Total CV flux
!    		MD_rhouu_dS_prev = MD_rhouu_dS
!        	MD_rhouu_dS  =		((CV_constraint%flux(ib,jb,kb,:,1)-CV_constraint%flux(ib,jb,kb,:,4)) &
!        	          	 		+(CV_constraint%flux(ib,jb,kb,:,2)-CV_constraint%flux(ib,jb,kb,:,5)) &
!        	          	 		+(CV_constraint%flux(ib,jb,kb,:,3)-CV_constraint%flux(ib,jb,kb,:,6)))/delta_t

!        	!if (abs(sum(MD_rhouu_dS_prev)) .ne. abs(sum(MD_rhouu_dS))) then	
!    		!	print'(a,i8,2f18.7)', 'Updated', iter, sum(MD_rhouu_dS), sum(MD_rhouu_dS_afr)
!    		!endif
!			!if (all(MD_rhouu_dS_prev .ne. MD_rhouu_dS)) then	

!			!endif

!    		!If fluxes have converged (exactly) then exit
!    		if (all(MD_rhouu_dS_prev .eq. MD_rhouu_dS)) converged = .true.

!    	else
!    		F_constraint = 0.d0
!			!list_afr = list_bfr
!			!box_np_afr = box_np_bfr
!    	endif

!	enddo
!	deallocate(r_temp,v_temp)

!end subroutine average_over_bin_iter

!!=============================================================================
!! Apply force to molecules in overlap region
!!-----------------------------------------------------------------------------

!subroutine apply_force
!	use arrays_MD, only : r,v,a
!	use physical_constants_MD, only : density
!	use computational_constants_MD, only : eflux_outflag
!	use module_record_external_forces
!	implicit none

!	integer								:: n, molno
!	double precision,dimension(3)		:: velvect

!	!Loop over all molecules and apply constraint
!	do n = 1, box_np_bfr
!		!Get molecule number
!		molno = list_bfr(n)
!		!Apply force
!		a(:,molno) = a(:,molno) - F_constraint
!		!Add external force to CV total
!		if (eflux_outflag .eq. 4) then
!			velvect(:) = v(:,molno) + 0.5d0*delta_t*a(:,molno)
!			call record_external_forces(F_constraint(:),r(:,molno),velvect)
!		else
!			call record_external_forces(F_constraint(:),r(:,molno))
!		endif

!	enddo
!	deallocate(list_bfr)

!end subroutine apply_force

!! Testing routine to apply the force and compare evolution with the case
!! which doesn't evolve in time
!subroutine apply_force_tests(apply_the_force)
!	use arrays_MD, only : r,v,a
!	use physical_constants_MD, only : density, pi
!	use calculated_properties_MD, only : temperature
!	use module_set_parameters, only : velPDF
!	use computational_constants_MD, only : eflux_outflag
!	use module_record_external_forces
!	implicit none

!	logical, intent(in) 						:: apply_the_force
!	logical										:: skip
!	integer										:: i,j, n,m, unitno
!	integer,dimension(3)						:: bin
!	double precision,dimension(3)				:: velvect
!	double precision,dimension(:,:),allocatable	:: r_temp,v_temp,a_temp
!	double precision,dimension(:),allocatable 	:: vmagnitude,normalisedvfd_bin,binloc

!	!Calculate acceleration with constraint
!	allocate(r_temp(3,np),v_temp(3,np),a_temp(3,np))
!	a_temp(:,1:np) = a(:,1:np)
!	do i = 1, box_np_bfr
!		n = list_bfr(i)
!		a_temp(:,n) = a(:,n) - F_constraint
!		if (apply_the_force) then
!			!print'(a,i8,a,6f10.4)', 'Force2mol ', n, ' mag ', F_constraint, & 
!			!		dble(box_np_bfr)*F_constraint/product(binsize)
!			a(:,n) = a(:,n) - F_constraint
!			!Add external force to CV total
!			if (eflux_outflag .eq. 4) then
!				velvect(:) = v(:,n) + 0.5d0*delta_t*a(:,n)
!				call record_external_forces(F_constraint(:),r(:,n),velvect)
!			else
!				call record_external_forces(F_constraint(:),r(:,n))
!			endif
!		endif
!	enddo

!	box_np_afr = 0
!	if (.not. allocated(list_afr)) allocate(list_afr(np))

!	!Update velocity, positions and molecules in box
!	v_temp = 0.d0; r_temp = 0.d0
!	do n = 1,np
!		v_temp(:,n) = v(:,n) + delta_t * a_temp(:,n) 	!Velocity calculated from acceleration
!		r_temp(:,n) = r(:,n) + delta_t * v_temp(:,n)    !Position calculated from velocity

!		bin(:) = ceiling((r_temp(:,n)+0.5d0*domain(:))/binsize(:))+1
!		if (bin(1) .eq. 3 .and. bin(2) .eq. 3 .and. bin(3) .eq. 3) then
!			!Add new molecule to overlap list
!			box_np_afr   =  box_np_afr   + 1
!			list_afr(box_np_afr) =  n
!		endif
!	enddo

!	!Calculate acceleration on any new molecules with constraint
!! 	do i = 1, box_np_afr
!! 		n = list_afr(i)
!! 		skip = .false.
!! 		do j = 1,box_np_bfr
!! 			m = list_bfr(j)
!! 			if (m .eq. n) skip = .true.
!! 		enddo
!! 		if (skip) cycle
!! 		print'(a,4i8)', 'qqqq', iter, n, box_np_bfr, box_np_afr
!!    		!a_temp(:,n) = a(:,n) - F_constraint
!!    		!if (apply_the_force) then
!!    		!	a(:,n) = a(:,n) - F_constraint
!!    			!Add external force to CV total
!!    		!	if (eflux_outflag .eq. 4) then
!!    		!		call record_external_forces(F_constraint(:),r(:,n),v(:,n))
!!    		!	else
!!    		!		call record_external_forces(F_constraint(:),r(:,n))
!!    		!	endif
!!    		!endif
!! 	enddo


!	!Next save evolution with constraint applied (and PDF function)
!	u_bin_afr(:) = 0.d0
!	F_bin1(:) = 0.d0
!! 	allocate(vmagnitude(box_np_afr))
!! 	do i = 1, box_np_afr
!! 		n = list_afr(i)
!! 		u_bin_afr(:) = u_bin_afr(:) + v_temp(:,n)
!! 		F_bin1(:) = F_bin1(:) + a_temp(:,n)
!! 		vmagnitude(i) = v_temp(1,n)
!! 	enddo
!	allocate(vmagnitude(box_np_bfr))
!	do i = 1, box_np_bfr
!		n = list_bfr(i)
!		u_bin_afr(:) = u_bin_afr(:) + v_temp(:,n)
!		F_bin1(:) = F_bin1(:) + a_temp(:,n)
!		vmagnitude(i) = v_temp(1,n)
!	enddo
!	deallocate(v_temp)
!	deallocate(a_temp)
!	call velPDF%update(vmagnitude(:))
!	deallocate(vmagnitude)

!	!Normalise PDF bins and write out
!	if (mod(iter,nint(0.25*Nsteps)) .eq. 0) then
!		allocate(normalisedvfd_bin,source=velPDF%normalise())
!		allocate(binloc,source=velPDF%binvalues())
!		do n=1,size(normalisedvfd_bin,1) 
!			write(12,'(4(f10.5))') binloc(n), normalisedvfd_bin(n), &
!									sqrt(2/pi)*((binloc(n)**2)*exp((-binloc(n)**2) & 
!									/(2*temperature))/(temperature**(3.d0/2.d0))), & 
!								    (1.d0/(sqrt(temperature)*sqrt(2.d0*pi))) & 
!									*exp( -((binloc(n))**2.d0-u_bin_afr(1)/volume) & 
!									/(2.d0*temperature) ) 
!		enddo
!		velPDF%hist = 0
!	endif

!	!Write values of constrained and unconstrained velocity
!	if (box_np_afr .ne. 0) then
!		unitno = irank*1000000+ibin*10000+jbin*100+kbin
!		write(unitno,'(i3,6i5,9f11.6)'),irank,iter,ibin,jbin,kbin,box_np_bfr,box_np_afr, &
!					 					  delta_t*F_constraint,u_bin_afr/volume,u_bin_bfr/volume
!										
!	else
!		unitno = irank*1000000+ibin*10000+jbin*100+kbin
!		write(unitno,'(i3,6i5,9f11.6)'),irank,iter,ibin,jbin,kbin,box_np_bfr,box_np_afr, &
!					 					 (/ 0.d0, 0.d0, 0.d0 /),u_bin_afr/volume,u_bin_bfr/volume
!	endif
!	deallocate(list_bfr, list_afr)

!end subroutine apply_force_tests


!end subroutine apply_CV_force



