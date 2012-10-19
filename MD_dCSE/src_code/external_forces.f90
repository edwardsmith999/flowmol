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
	integer						:: ibin, jbin, kbin
	integer						:: averagecount
	double precision 			:: F, fixdist, slicebinsize
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
			call compute_bin_surface_flux(icell,jcell,kcell,isumflux)
			call compute_force_surrounding_bins(icell,jcell,kcell,isumforce)

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


	
subroutine get_cell_ranges
	use coupler_module
	use computational_constants_MD, only : npx,npy,npz, domain
	implicit none

	integer				:: nlgx_md,nlgy_md, nlgz_md, yLl_md,yLl_cfd
	integer 			:: i,n,startproc
	integer, parameter 	:: Tnull = -666

	allocate(icPmin_md(npx)); icPmin_md = Tnull
	allocate(jcPmin_md(npy)); jcPmin_md = Tnull
	allocate(kcPmin_md(npz)); kcPmin_md = Tnull
	allocate(icPmax_md(npx)); icPmax_md = Tnull
	allocate(jcPmax_md(npy)); jcPmax_md = Tnull
	allocate(kcPmax_md(npz)); kcPmax_md = Tnull

	xL_cfd = 181.1; yL_cfd = 181.1; zL_cfd = 10.0
	ncx = 64;	ncy = 64;	ncz = 8
	ncy_olap = 5
	yL_md = domain(2)

  	allocate( xg(ncx+1,ncy+1) , yg(ncx+1,ncy+1), zg(ncz+1) )
	!----- X grid ------
	dx = xL_cfd  / ncx
	do i=1,ncx
		xg(i,:) = (i-1.) * dx
	end do

	!----- Y grid ------
	dy = yL_cfd  / ncy
	do i=1,ncy
		yg(i,:) = (i-1.) * dy
	end do

	!----- Z grid ------
	dz = zL_cfd  / ncz
	do i=1,ncz
		zg(i) = (i-1.) * dz
	end do
		
	! - - x - -
	nlgx_md = nint(dble(ncx)/dble(npx))
	do n=1,npx
		icPmax_md(n) = n * nlgx_md
		icPmin_md(n) = icPmax_md(n) - nlgx_md + 1
	end do	

	! - - y - -
	nlgy_md = ceiling(dble(ncy)/dble(npy))
	yL_olap = ncy_olap * dy
	yL_puremd = yL_md - yL_olap
	ncy_puremd = yL_puremd / dy
	yLl_md = yL_md / npy
	startproc = ceiling(yL_puremd/yLl_md)
	do n = startproc,npy
		jcPmax_md(n) = n * nlgy_md - ncy_puremd
		jcPmin_md(n) = jcPmax_md(n) - nlgy_md + 1
		if (jcPmin_md(n).le.0) jcPmin_md(n) = 1
	end do 

	! - - z - -
	nlgz_md = nint(dble(ncz)/dble(npz))
	do n=1,npz
		kcPmax_md(n) = n * nlgz_md
		kcPmin_md(n) = kcPmax_md(n) - nlgz_md + 1
	end do

end subroutine get_cell_ranges

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
