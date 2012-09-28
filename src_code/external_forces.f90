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


module coupler_module

	! <=><=><=><=> Grid and domain data <=><=><=><=> 
	! Density
	real(kind(0.d0))	:: density_cfd, density_md
	! CFD/MD number of cells
    integer	:: ngx,ngy,ngz						!Global
	integer	:: nlgx_cfd,nlgy_cfd,nlgz_cfd, & 	!Local CFD
			   nlgx_md ,nlgy_md ,nlgz_md 		!Local MD
	! Overlap cells 
    integer :: i_olap,j_olap,k_olap	
	! MD grid indices
	integer	:: ngy_puremd
	! CFD grid indices
    integer	:: imin,imax,jmin,jmax,kmin,kmax 
	! CFD/MD local grid indices (start and end of grid per CFD/MD process)
    integer,dimension(:),allocatable :: iTmin_cfd,iTmax_cfd,jTmin_cfd,jTmax_cfd,kTmin_cfd,kTmax_cfd, & 
										iTmin_md ,iTmax_md ,jTmin_md ,jTmax_md ,kTmin_md ,kTmax_md 
	! Domain sizes
	real(kind(0.d0)) ::	xL_md  ,yL_md  ,zL_md , & 
						xL_cfd ,yL_cfd ,zL_cfd, &
						xL_olap,yL_olap,zL_olap,&
								yL_puremd
	! Local Domain
	real(kind(0.d0)) :: yLl_md, yLl_cfd

	!CFD cells sizes 
	real(kind(0.d0)) 								   :: dx,dymin,dy,dymax,dz
    real(kind(0.d0)),dimension(:),  allocatable,target :: zpg
    real(kind(0.d0)),dimension(:,:),allocatable,target :: xpg,ypg

	! <=><=><=><=> Processor Topology & MPI <=><=><=><=> 
	!Global processor number across both realms
	integer	:: myid
	!Processor id in grid
	integer	:: myid_grid,iblock,jblock,kblock
	! Number of processor in CFD grid
	integer :: npx_cfd, npy_cfd, npz_cfd, nproc_cfd	
    ! Number of processor in MD grid
    integer :: npx_md,  npy_md,  npz_md,  nproc_md
    ! Coordinates of MD/CFD topologies
	integer,dimension(:,:),allocatable 	 :: rank2coord_cfd, rank2coord_md
    integer,dimension(:,:,:),allocatable :: coord2rank_cfd, coord2rank_md
	!Mapping between CFD processors and MD processors
    integer,dimension(:,:),allocatable :: imap_olap,jmap_olap,kmap_olap

end module coupler_module

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

!=============================================================================
! Apply force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
! Adapted serial version written by ES including cells and (originally) fully verified
!-----------------------------------------------------------------------------

subroutine apply_continuum_forces_ES(iter)
	use computational_constants_MD, only : delta_t,nh,halfdomain,ncells, & 
											cellsidelength,initialstep,Nsteps,iblock,jblock,kblock, & 
											npx,npy,npz
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	use coupler_module, only : iTmin_md,iTmax_md,jTmin_md,jTmax_md,kTmin_md,kTmax_md
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

	!allocate(u_continuum(iTmin_md(iblock):iTmax_md(iblock), & 
	!					 jTmin_md(jblock):jTmax_md(jblock), & 
	!					 kTmin_md(kblock):kTmax_md(kblock)))

	!	print'(a,6i8)', 'limits', iTmin_md(iblock),iTmax_md(iblock),jTmin_md(jblock),jTmax_md(jblock),kTmin_md(kblock),kTmax_md(kblock)


	!do ii=iTmin_md(iblock),iTmax_md(iblock)
	!do jj=jTmin_md(jblock),jTmax_md(jblock)
	!do kk=kTmin_md(kblock),kTmax_md(kblock)

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
											cellsidelength,initialstep,Nsteps,iblock,jblock,kblock, & 
											npx,npy,npz
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	use coupler_module, only : iTmin_md,iTmax_md,jTmin_md,jTmax_md,kTmin_md,kTmax_md
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

	!allocate(u_continuum(iTmin_md(iblock):iTmax_md(iblock), & 
	!					 jTmin_md(jblock):jTmax_md(jblock), & 
	!					 kTmin_md(kblock):kTmax_md(kblock)))

	!	print'(a,6i8)', 'limits', iTmin_md(iblock),iTmax_md(iblock),jTmin_md(jblock),jTmax_md(jblock),kTmin_md(kblock),kTmax_md(kblock)


	!do ii=iTmin_md(iblock),iTmax_md(iblock)
	!do jj=jTmin_md(jblock),jTmax_md(jblock)
	!do kk=kTmin_md(kblock),kTmax_md(kblock)
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

	integer 			:: i,n,startproc
	integer, parameter 	:: Tnull = -666

	allocate(iTmin_md(npx)); iTmin_md = Tnull
	allocate(jTmin_md(npy)); jTmin_md = Tnull
	allocate(kTmin_md(npz)); kTmin_md = Tnull
	allocate(iTmax_md(npx)); iTmax_md = Tnull
	allocate(jTmax_md(npy)); jTmax_md = Tnull
	allocate(kTmax_md(npz)); kTmax_md = Tnull

	xL_cfd = 181.1; yL_cfd = 181.1; zL_cfd = 10.0
	ngx = 64;	ngy = 64;	ngz = 8
	j_olap = 5
	yL_md = domain(2)

  	allocate( xpg(ngx+1,ngy+1) , ypg(ngx+1,ngy+1), zpg(ngz+1) )
	!----- X grid ------
	dx = xL_cfd  / ngx
	do i=1,ngx
		xpg(i,:) = (i-1.) * dx
	end do

	!----- Y grid ------
	dy = yL_cfd  / ngy
	do i=1,ngy
		ypg(i,:) = (i-1.) * dy
	end do

	!----- Z grid ------
	dz = zL_cfd  / ngz
	do i=1,ngz
		zpg(i) = (i-1.) * dz
	end do
		
	! - - x - -
	nlgx_md = ceiling(dble(ngx)/dble(npx))
	do n=1,npx
		iTmax_md(n) = n * nlgx_md
		iTmin_md(n) = iTmax_md(n) - nlgx_md + 1
	end do	

	! - - y - -
	nlgy_md = ceiling(dble(ngy)/dble(npy))
	yL_olap = j_olap * dy
	yL_puremd = yL_md - yL_olap
	ngy_puremd = yL_puremd / dy
	yLl_md = yL_md / npy
	startproc = ceiling(yL_puremd/yLl_md)
	do n = startproc,npy
		jTmax_md(n) = n * nlgy_md - ngy_puremd
		jTmin_md(n) = jTmax_md(n) - nlgy_md + 1
		if (jTmin_md(n).le.0) jTmin_md(n) = 1
	end do 

	! - - z - -
	nlgz_md = ceiling(dble(ngz)/dble(npz))
	do n=1,npz
		kTmax_md(n) = n * nlgz_md
		kTmin_md(n) = kTmax_md(n) - nlgz_md + 1
	end do

end subroutine get_cell_ranges

!------------------------------------------
!subroutine CFD_cells_to_MD_compute_cells(ibmin_cfd,imin_cfd,jbmin_cfd,jmin_cfd,kbmin_cfd,kmin_cfd, & 
!										  ibmin_md, imin_md, jbmin_md, jmin_md, kbmin_md, kmin_md)
!	implicit none

!	integer,intent(in)		:: ibmin_cfd,imin_cfd,jbmin_cfd,jmin_cfd,kbmin_cfd,kmin_cfd
!	integer,intent(out)		:: ibmin_md,imin_md,jbmin_md,jmin_md,kbmin_md,kmin_md

subroutine CFD_cells_to_MD_compute_cells(ii_cfd,jj_cfd,kk_cfd, & 
										  ibmin_md, ibmax_md, jbmin_md, jbmax_md, kbmin_md, kbmax_md)
	use coupler_module, only : xpg, ypg, zpg, xL_md,yL_md,zL_md, iblock,jblock,kblock
	use computational_constants_MD, only : cellsidelength
	implicit none

	integer,intent(in)		:: ii_cfd,jj_cfd,kk_cfd
	integer,intent(out)		:: ibmin_md,ibmax_md,jbmin_md,jbmax_md,kbmin_md,kbmax_md
	
	double precision		:: xL_min,xL_max,yL_min,yL_max,zL_min,zL_max

	! Get minimum point in processors domain
	xL_min = xL_md*(iblock-1); xL_max = xL_md*(iblock)
	yL_min = yL_md*(jblock-1); yL_max = yL_md*(jblock)
	zL_min = zL_md*(kblock-1); zL_max = zL_md*(kblock)

	! Get range of cells to check so that top and bottom of current CFD cell are covered
	ibmin_md = (xpg(ii_cfd  ,jj_cfd  )-xL_min)/cellsidelength(1)+1
	ibmax_md = (xpg(ii_cfd+1,jj_cfd  )-xL_min)/cellsidelength(1)+1
	jbmin_md = (ypg(ii_cfd  ,jj_cfd  )-yL_min)/cellsidelength(2)+1
	jbmax_md = (ypg(ii_cfd  ,jj_cfd+1)-yL_min)/cellsidelength(2)+1
	kbmin_md = (zpg(     kk_cfd      )-zL_min)/cellsidelength(3)+1
	kbmax_md = (zpg(     kk_cfd+1    )-zL_min)/cellsidelength(3)+1

	print'(a,9i8)','indices', ii_cfd,ibmin_md,ibmax_md,jj_cfd,jbmin_md,jbmax_md,kk_cfd,kbmin_md,kbmax_md
	print*,'xcells', xpg(ii_cfd  ,jj_cfd  ),(xpg(ii_cfd  ,jj_cfd  )-xL_min)/cellsidelength(1)+1, (xpg(ii_cfd+1,jj_cfd  )-xL_min)/cellsidelength(1)+1
	print*,'ycells', ypg(ii_cfd  ,jj_cfd  ),(ypg(ii_cfd  ,jj_cfd  )-yL_min)/cellsidelength(2)+1, (ypg(ii_cfd+1,jj_cfd  )-yL_min)/cellsidelength(2)+1
	print*,'zcells', zpg(kk_cfd  ),(zpg(kk_cfd)-zL_min)/cellsidelength(3)+1, (zpg(kk_cfd+1)-zL_min)/cellsidelength(3)+1

end subroutine CFD_cells_to_MD_compute_cells
