!=============================================================================
!				   	MD coupler Socket
! Routine which interface with the coupler to the CFD code
!
! socket_coupler_init						Passes MD initialisation variables 
!											to coupler_md_init
! average_and_send_MD_to_CFD				Average MD and pass to CFD as BC
! 	coupler_md_boundary_cell_average		Calls routines for uc,vc and wc
! 	  compute_uc_average					Accumulate averages and send
!	  compute_vc_average					Accumulate averages and send
!	  compute_wc_average					Accumulate averages and send
! socket_coupler_apply_continuum_forces		Apply CFD constraint on MD 
!	apply_continuum_forces					Calls routines for setup, average 
!											& application of constraint
!	  call setup_CFD_box					Setup type storing average for cell
!	  call average_over_bin					Accumulate averages
!	  call apply_force						Apply force
!
!=============================================================================


module md_coupler_socket

#if USE_COUPLER 

	use coupler
	implicit none

    ! CFD id
    integer cfd_code_id 

	! type used in applying continuum force
	type cfd_box_sum
		integer np
		real(kind=kind(0.d0))  v(3)
		real(kind=kind(0.d0))  a(3)
	end type cfd_box_sum


	type(cfd_box_sum),allocatable :: box_average(:,:,:)

	integer			, dimension(:,:,:,:), allocatable :: mflux
	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: uvw_md, uvw_cfd

contains

!=============================================================================
! Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup
!
!								SETUP
!
! Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup
!=============================================================================

!=============================================================================
! Invoke messenger and split inter-communicators 
!-----------------------------------------------------------------------------
subroutine socket_coupler_invoke
	use messenger
	implicit none

	call CPL_create_comm(md_realm,MD_COMM,ierr)
    prefix_dir = "./md_data/"

end subroutine socket_coupler_invoke

!=============================================================================
!  Read coupler input files
!-----------------------------------------------------------------------------
subroutine socket_read_coupler_input
	use messenger
	use coupler_input_data, only : read_coupler_input
    use coupler_module, only : request_stop
	implicit none

	call read_coupler_input		! Read COUPLER.in input file

    ! stop if requested ( useful for development )
	call request_stop("create_comm") ! stops here if in COUPLER.in stop requestis set to "create_comm"

end subroutine socket_read_coupler_input

!=============================================================================
! Setup initial times based on coupled calculation
!-----------------------------------------------------------------------------
subroutine socket_coupler_init
    use interfaces
	use coupler_input_data, only : md_cfd_match_cellsize
	use computational_constants_MD, only : npx,npy,npz,delta_t,elapsedtime, & 
										   Nsteps,initialstep,delta_t, & 
										   globaldomain,initialnunits
	use physical_constants_MD, only : nd,density
	use messenger, only	 :  myid, icoord, icomm_grid
	use coupler_internal_md, only : coupler_md_init
	implicit none

 	integer			 :: naverage, nsteps_cfd, ixyz
	real(kind(0.d0)) :: delta_t_CFD

	!Establish Domain size from MD inputs
	do ixyz=1,nd
		globaldomain(ixyz) = initialnunits(ixyz) & 	
								/((density/4.d0)**(1.d0/nd))
	enddo

	! If coupled calculation prepare exchange layout
	call coupler_md_init(nsteps,delta_t,icomm_grid,icoord,(/ npx,npy,npz /),globaldomain,density)

	! Setup the domain and cells in the MD based on coupled data
	call set_parameters_global_domain_coupled
	if (md_cfd_match_cellsize .eq. 0) then
		call set_parameters_cells
 	else
		call set_parameters_cells_coupled
	endif

	! Setup timesteps and simulation timings based on CFD/coupler
	call set_coupled_timing(Nsteps)

end subroutine socket_coupler_init

!=============================================================================
! get the global domain lenghts from x, y, z array of CFD realm

subroutine set_parameters_global_domain_coupled
	use computational_constants_MD
	use physical_constants_MD
	use coupler 
	use coupler_module, b0 => MD_initial_cellsize
	implicit none

	integer          ixyz, n0(3)

    ! fix the numner of FCC cells starting from CFD density
    !density = coupler_md_get_density()

    ! size of cubic FCC cell
    b0=(4.d0/density)**(1.0d0/3.0d0)
    !call coupler_md_get(xL_md=xL_md,yL_md=yL_md,zL_md=zL_md,MD_initial_cellsize=b0)
    n0(:) = nint( (/ xL_md, yL_md, zL_md/) / b0)
    initialunitsize(1:3) = b0
    initialnunits(1:3) 	 = n0(:)
   
	!Set MD domain values
	globaldomain(1) = xL_md
	globaldomain(2) = yL_md
	globaldomain(3) = zL_md
	volume   = xL_md*yL_md*zL_md

    ! no need to fix globalnp if we have it already
    if(.not. restart) then
		globalnp=1      !Set number of particles to unity for loop below
		do ixyz=1,nd
			globalnp = globalnp*initialnunits(ixyz)		!One particle per unit cell
		enddo
		globalnp=4*globalnp   !FCC structure in 3D had 4 molecules per unit cell
    endif

    np = globalnp / nproc	

	domain(1) = globaldomain(1) / real(npx, kind(0.d0))			!determine domain size per processor
	domain(2) = globaldomain(2) / real(npy, kind(0.d0))			!determine domain size per processor
	domain(3) = globaldomain(3) / real(npz, kind(0.d0))			!determine domain size per processor

	do ixyz=1,nd
		halfdomain(ixyz) = 0.5d0*domain(ixyz)			!Useful definition
	enddo

	!write(0,*) 'set_parameter_global_domain_hybrid ', globalnp, np, domain, initialunitsize

    if(myid_world .eq. rootid_world) then
        write(*,'(a/a/a,f5.2,a,f5.2,/,a,3(f5.2),a,/,a,3(I6),/,a)') &
                "**********************************************************************", &
                "WARNING - this is a coupled run which resets the following parameters:", &
                " density         =", density ,                                           & 
                " hence the cubic FCC side  is b=", b0 ,                                  &
                " initialunitsize =", initialunitsize(:)/b0," in b units ",               &     
                " initialnunits   =", initialnunits(:),                                   &
                "**********************************************************************"
    endif

end subroutine set_parameters_global_domain_coupled

!-----------------------------------------------------------------------------
! Adjust MD cells to match to continuum

subroutine set_parameters_cells_coupled
	use interfaces
	use computational_constants_MD
	use physical_constants_MD
	use polymer_info_MD
	use coupler
	use coupler_module, only : icmax,icmin,jcmax,jcmin,kcmax,kcmin,dx,dy,dz,rank_realm
	implicit none

	integer 						:: ixyz
	integer,dimension(3) 			:: max_ncells,cfd_ncells,cfd_md_cell_ratio
	double precision 				:: rneighbr
	double precision ,dimension(3) 	:: cfd_cellsidelength, maxdelta_rneighbr
    type(cfd_grid_info)				:: cfd_box

	!In coupled simulation, passed properties are calculated from cell lists
	!for efficiency. The size of the cells should therefore be a multiple
	!of the continuum cellsizes
    cfd_ncells(1) = icmax - icmin
    cfd_ncells(2) = jcmax - jcmin
    cfd_ncells(3) = kcmax - kcmin

	!Check number of cells based on rcutoff and neighbourlist size
	max_ncells= floor(domain/rcutoff)

	!Calculate ratio of CFD to maximum numbers of cells
	cfd_cellsidelength(1) = dx
	cfd_cellsidelength(2) = dy
	cfd_cellsidelength(3) = dz
	cellsidelength(:) 	  = domain(:)/max_ncells(:)
	cfd_md_cell_ratio(:)  = floor(cfd_cellsidelength(:)/cellsidelength(:))

	!Determine side length of cells after rounding and MD ncells
	cellsidelength = cfd_cellsidelength/cfd_md_cell_ratio
	ncells = ceiling(domain/cellsidelength)
	!print*, 'IS this ok?',cfd_cellsidelength, ncells, cellsidelength,  domain,size(x), size(y)

	!Ensure domain allows MD cells to be a multiple of CFD cell sizes...
	if (any(abs(domain(:)/cellsidelength(:)-nint(domain(:)/cellsidelength(:))) .gt. 0.01)) & 
		call error_abort("ERROR - CFD cellsize and MD cellsize not compatible - Adjust domain size to      &
						&  correct this or remove MD_CFD_MATCH_CELLSIZE from COUPLER input")
	cellsidelength = domain/ncells

	!Recalculate required delta_rneighbr to ensure integer numbers of cells for both domains 
	delta_rneighbr = minval(domain/ncells-rcutoff)

	!Calculate size of neighbour list region
	rneighbr  = rcutoff + delta_rneighbr
	rneighbr2 = (rcutoff + delta_rneighbr)**2

	!print'(a,6f10.5)', 'domains', domain, x(icmax)-x(icmin),y(jcmax)-y(jcmin),z(kcmax)-z(kcmin)
	!print'(a,12i8)',      'cell', cfd_md_cell_ratio,cfd_ncells,ncells,max_ncells
	!print'(a,6f10.5)', 'cellsize', cfd_cellsidelength(:),cellsidelength(:)

	if (potential_flag .eq. 1) then
		select case(solvent_flag)
		case(0:1)
			if (rneighbr < R_0) then
				rneighbr = R_0 
				rneighbr2 = R_0**2
				print*, 'Neighbour list distance rneighbr set to &
						& maximum elongation of polymer spring, ',R_0
			end if
		case(2)
			if (rneighbr < sod_cut) then
				rcutoff   = sod_cut
				rcutoff2  = sod_cut2
				rneighbr  = rcutoff + delta_rneighbr
				rneighbr2 = rneighbr**2.d0
			end if
		case default
			call error_abort('Unrecognised solvent_flag in set_parameters_cells')
		end select
	endif

	if (ncells(1)<3 .or. ncells(2)<3 .or. ncells(3)<3) then
		print*, ncells(1),'    in x and ', ncells(2), '    in y' , ncells(3), '    in z' 
		call  error_abort( "ERROR - DOMAIN SHOULD HAVE AT LEAST 3 CELLS, &
		 					& IN X, Y AND Z - INCREASE NUMBER OF UNITS IN INPUT")
	endif

    if(rank_realm .eq. 1) then
        write(*,'(a/a/a,f8.6,/,a,3i8,/,a,3(f8.5),/,a,3(i8),/,a)') &
                    "**********************************************************************", &
                    "WARNING - this is a coupled run which resets the following parameters:", &
	    	        " Extra cell size for neighbourlist =", delta_rneighbr  ,                 & 
                    " MD computational cells per CFD cell = ",cfd_md_cell_ratio,  			  &
					" cellsize =", cellsidelength(:),     									  &     
                    " no of cell  =", ncells(:),                                   &
                    "**********************************************************************"
    endif

	!Ensure CFD cells are not smaller than minimum possible MD cells
	if (any(max_ncells .lt. cfd_ncells)) & 
		call error_abort("ERROR - CFD cellsize smaller than minimum MD computational/averaging cell")

end subroutine set_parameters_cells_coupled

!-----------------------------------------------------------------------------
!	THIS ROUTINE SHOULD BE PART OF THE INITIALISATION OF THE COUPLER
!	WHERE THE INITIAL TIMES (STIME FOR CFD & ELAPSEDTIME FOR MD) ARE
!	COMPARED TO SEE IF RESTART IS CONSISTENT AND NUMBER OF STEPS ON
!	BOTH SIDES OF THE COUPLER ARE CALCULATED AND STORED
! Setup ratio of CFD to MD timing and total number of timesteps

subroutine set_coupled_timing(Nsteps_md)
	use computational_constants_MD, only : initialstep,elapsedtime
	use coupler_input_data, only : md_steps_per_dt_cfd,md_steps_per_dt_cfd_tag
	use coupler_module, only : dt_cfd,dt_MD,Nsteps_cfd,Nsteps_md_old=>Nsteps_md, & 
							   rank_realm,Nsteps_coupled
	implicit none

	integer,intent(out)		:: Nsteps_md
	integer					:: Nsteps_MDperCFD

	!Set number of MD timesteps per CFD using ratio of timestep or coupler value
	if(md_steps_per_dt_cfd_tag == CPL) then
		Nsteps_MDperCFD = md_steps_per_dt_cfd
	else 
		Nsteps_MDperCFD = int(dt_cfd/dt_MD)
	endif
	Nsteps_coupled = Nsteps_cfd

 	!Set number of steps in MD simulation, Nsteps
	Nsteps_md   = initialstep + Nsteps_cfd * Nsteps_MDperCFD
	elapsedtime = elapsedtime + Nsteps_cfd * Nsteps_MD * dt_MD

	if (rank_realm .eq. 1) then 
		write(*,'(2(a,/),a,i7,a,i7,/a,i7,a,i7,/a,f12.4,a,/a)') &
			"*********************************************************************", 		&
			" WARNING - WARNING - WARNING - WARNING - WARNING - WARNING - WARNING  ", 		&
			" Input number of timesteps from MD: ",Nsteps_md_old," & CFD: ", Nsteps_cfd,	&
			" is set in this coupled run to MD: ", Nsteps_md, ",CFD/Coupled: ", Nsteps_cfd,	&
			" At the end of the run, elapsed time will be: ", elapsedtime, " LJ time units ", 		&
			"*********************************************************************"   
	endif 

end subroutine set_coupled_timing
!=============================================================================
! Establish mapping between MD an CFD
!-----------------------------------------------------------------------------
subroutine socket_create_map
	implicit none

	!Note coupler cannot be called directly so this socket is needed
	call CPL_create_map

end subroutine socket_create_map

!=============================================================================
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!
!							SIMULATION
!
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!=============================================================================

!=============================================================================
! Take average of x,y and z components of MD velocity to 
! calculate all components of velocity to pass to contiunuum region 
!
!-----------------------------------------------------------------------------
subroutine average_and_send_MD_to_CFD(iter)
	use computational_constants_MD, only : initialstep, delta_t, nhb
	use calculated_properties_MD, only : nbins
	use physical_constants_MD, only : np
	use arrays_MD, only :r,v
   	use coupler_module, only : staggered_averages, ncx,ncy,ncz, & 
							   dx,dz,yg,cfd_code_id,md_realm,rank_realm
	use coupler_input_data, only : md_steps_per_dt_cfd
	use messenger, only : MD_COMM
	implicit none

	integer, intent(in) :: iter
	
	integer :: ixyz,icell,jcell,kcell,pcoords(3),extents(6)
	integer :: iter_cfd, iter_average, save_period, average_period
	integer	:: ierr
	logical, save :: first_time=.true.
	save  average_period

	!Setup arrays on first call
    if (first_time) then 
	    first_time	= .false.
		call setup_velocity_average
    endif
    iter_average = mod(iter-1, md_steps_per_dt_cfd)+1			! current step
    iter_cfd     = (iter-initialstep)/md_steps_per_dt_cfd +1 	! CFD corresponding step

	!Collect uc data every save_period cfd iteration but discard the first one which cfd uses for initialisation
    if ( mod(iter_average,average_period) .eq. 0 ) then
	    call cumulative_velocity_average
	endif

	!Send accumulated results to CFD at the end of average cycle 
    if  (iter_average .eq. md_steps_per_dt_cfd) then
	    call send_velocity_average
	endif

contains

	!THIS SHOULD BE DONE IN THE SETUP!!!!
	subroutine setup_velocity_average
	   	use coupler_module, only : iblock_realm,jblock_realm,kblock_realm 
		implicit none

		integer		:: nclx,ncly,nclz
		integer		:: pcoords(3),extents(6)

	    average_period = coupler_md_get_average_period() 	! collection interval in the average cycle

		!Allocate array to size of cells in current processor
		pcoords = (/ iblock_realm,jblock_realm,kblock_realm /)
		call CPL_proc_extents(pcoords,md_realm,extents)
		nclx = extents(2)-extents(1)+1
		ncly = extents(4)-extents(3)+1
		nclz = extents(6)-extents(5)+1

		! Setup averaging array on first call
		select case(staggered_averages(1))
		case(.true.)
			allocate( mflux(6,nclx,ncly,nclz))
			mflux = 0
		case(.false.)
			allocate(uvw_md(4,nclx,ncly,nclz))
			uvw_md = 0.d0
		end select
	end subroutine setup_velocity_average

!------------------------------------------
	subroutine cumulative_velocity_average
		use coupler, only : CPL_Cart_coords, CPL_proc_extents
		use coupler_module, only : xL_md,zL_md,xg,yg,zg,jcmin_olap, & 
								   dx,dy,dz,ncx,ncy,ncz,cfd_code_id, CPL_REALM_COMM,rank_world
		use coupler_internal_md, only : map_md2cfd_global,map_cfd2md_global,globalise,localise	
		use computational_constants_MD, only : iter, ncells,domain,halfdomain
		use librarymod, only : heaviside, imaxloc
		implicit none

		!Limits of cells to average

		integer							:: n, ixyz,extents(6),pcoords(3)
		integer,dimension(3)			:: ibin,ibin1,ibin2,minbin,maxbin,crossplane,cfdbins
		double precision	 			:: xbcmin,xbcmax, ybcmin,ybcmax, zbcmin,zbcmax
		double precision,dimension(3) 	:: cfd_cellsidelength
		double precision,dimension(3)	:: dxyz,ri1,ri2,avrg_top,avrg_bot,rd,rd2

		!Specify BC region to average molecules
		ybcmin = yg(1,jcmin_olap)-dy;    ybcmax = yg(1,jcmin_olap)

		!Velocity measurement for 3D bins throughout the domain
		!Determine bin size
		dxyz = (/ dx, dy, dz /)

		!Eliminate processors outside of passing region
		call CPL_Cart_coords(CPL_REALM_COMM,rank_realm,md_realm,3,pcoords,ierr)
		call CPL_proc_extents(pcoords,md_realm,extents)
		if (any(extents .eq. VOID)) return
		if ((yg(1,extents(3)) .gt. ybcmax) .or. (yg(1,extents(4)+1) .lt. ybcmin)) return

		!Get local extents on processor(s) of interest
		rd(:) = (/ 0.d0 , ybcmin, 0.d0 /)	!Bottom of cell below domain
		avrg_bot = localise(map_cfd2md_global(rd))
		rd2(:) = (/ xL_md , ybcmax , zL_md   /)   !Top of cell below domain (= bottom of domain)
		avrg_top = localise(map_cfd2md_global(rd2))

		minbin = ceiling((avrg_bot+halfdomain(:))/dxyz(:)) + nhb
		maxbin = ceiling((avrg_top+halfdomain(:))/dxyz(:)) + nhb

		!print'(a,16i5)','extents',rank_world, ceiling(map_cfd2md_global(rd)/dy), minbin,extents(1),extents(3),extents(5),maxbin,extents(2),extents(4),extents(6)

		select case(staggered_averages(1))	
		!- - - - - - - - - - - - - - - - - - - -
		!Record velocity flux over surface
		!- - - - - - - - - - - - - - - - - - - -
		case(.true.)
			!Determine bin size
			do n = 1,np
				ri1(:) = r(:,n) 							!Molecule i at time t
				ri2(:) = r(:,n)	-delta_t*v(:,n)				!Molecule i at time t-dt
				!Assign to bins before and after using integer division
				ibin1(:) = ceiling((ri1+halfdomain(:))/dxyz(:)) + nhb
				ibin2(:) = ceiling((ri2+halfdomain(:))/dxyz(:)) + nhb
					
				!Exclude molecules outside of processor domain
				if (	  ibin1(2) .lt. minbin(2) .or. ibin1(2) .ge. maxbin(2) &
					.and. ibin2(2) .lt. minbin(2) .or. ibin2(2) .ge. maxbin(2) ) cycle

				!Replace Signum function with this functions which gives a
				!check for plane crossing and the correct sign 
				crossplane(:) =  ibin1(:) - ibin2(:)
				if (sum(abs(crossplane(:))) .ne. 0) then
					!Find which direction the surface is crossed
					!For simplicity, if more than one surface has been crossed surface fluxes of intermediate cells
					!are not included. This assumption => more reasonable as Delta_t => 0
					ixyz = imaxloc(abs(crossplane))
					!Add mass flux to the new bin surface count and take from the old
					mflux(ixyz+3*heaviside(-dble(crossplane(ixyz))),ibin1(3),ibin1(1),ibin1(2)) = & 
						mflux(ixyz+3*heaviside(-dble(crossplane(ixyz))),ibin1(3),ibin1(1),ibin1(2)) & 
							+ abs(crossplane(ixyz))
					mflux(ixyz+3*heaviside(dble(crossplane(ixyz))),ibin2(3),ibin2(1),ibin2(2)) = & 
						mflux(ixyz+3*heaviside(dble(crossplane(ixyz))),ibin2(3),ibin2(1),ibin2(2)) &
							- abs(crossplane(ixyz))
				endif
			enddo
		!- - - - - - - - - - - - - - - - - - - -
		!Record velocity in cell centre
		!- - - - - - - - - - - - - - - - - - - -
		!Add up current volume mass and momentum densities
		case(.false.)
			do n = 1,np
				!Get bin containing molecule
				ibin(:) = ceiling((r(:,n)+halfdomain(:))/dxyz(:)) + nhb

				!Exclude molecules outside of averaging region
				!if (any(ibin(:).lt.minbin(:)) .or. any(ibin(:).ge.maxbin(:))) cycle
				if (r(2,n).lt.avrg_bot(2) .or. r(2,n).gt.avrg_top(2) ) cycle

				!Exclude molecules which leave processor between rebuilds
				if (ibin(1).lt.1 .or. ibin(1).gt.nbins(1)) cycle
				if (ibin(3).lt.1 .or. ibin(3).gt.nbins(3)) cycle

				!Add velocity and molecular count to bin
				uvw_md(1:3,ibin(1),ibin(2)-minbin(2)+1,ibin(3)) = &
					uvw_md(1:3,ibin(1),ibin(2)-minbin(2)+1,ibin(3)) + v(:,n)
				uvw_md(4,  ibin(1),ibin(2)-minbin(2)+1,ibin(3)) = & 	
					uvw_md(4,  ibin(1),ibin(2)-minbin(2)+1,ibin(3)) + 1.d0

				!print'(a,6i4,3f9.5)', 'binning',ibin(1),ibin(2),ibin(3),minbin(2), & 
				!						rank_realm,n, r(2,n),avrg_bot(2),avrg_top(2)
			enddo
		case default
			call error_abort('Unknown case in staggered_averages')
		end select

	end subroutine cumulative_velocity_average

!------------------------------------------
	subroutine send_velocity_average
		use coupler_module, only : olap_mask,rank_world, & 
								   iblock_realm,jblock_realm,kblock_realm, &
								   icmin_olap,icmax_olap,kcmin_olap,kcmax_olap,printf
		implicit none

		logical :: send_flag
        logical :: ovr_box_x
		integer :: jcmin_send,jcmax_send,limits(6)

		jcmin_send = 1; jcmax_send = 1

		!Send data to CFD if send_data flag is set
		select case(staggered_averages(1))	
		! Send velocity flux over surface
		case(.true.)
            call CPL_send(dble(mflux))!,index_transpose=(/2,3,1/))
			mflux = 0
		! Send velocity in cell centre
		case(.false.)
			!limits = (/ icmin_olap,icmax_olap,jcmin_send,jcmax_send,kcmin_olap,kcmax_olap /)
			!call CPL_gather(uvw_md,limits)

            call CPL_send(uvw_md,jcmax_send=jcmax_send,jcmin_send=jcmin_send,send_flag=send_flag)
			uvw_md = 0.d0

		end select

	end subroutine send_velocity_average

end subroutine average_and_send_MD_to_CFD






!=============================================================================
! Apply coupling forces so MD => CFD
! Force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
!-----------------------------------------------------------------------------
subroutine socket_apply_continuum_forces(iter)
	use physical_constants_MD, only : np
	use computational_constants_MD, only : delta_t, nh, ncells, & 
										cellsidelength, halfdomain, delta_rneighbr
	use coupler_input_data, only : md_steps_per_dt_cfd
	use coupler_module, only : rank_world,olap_mask, icmin_olap,icmax_olap, & 
								jcmin_olap,jcmax_olap,kcmin_olap,kcmax_olap, printf
	use linked_list
	implicit none

	integer, intent(in) 	:: iter ! iteration step, it assumes that each MD average starts from iter = 1


	integer 				:: iter_average, limits(6)
	integer					:: i,j,k,n,np_overlap
	integer,allocatable 	:: list(:,:)
	real(kind=kind(0.d0))	:: inv_dtCFD,t_fract,CFD_box(6)

	integer,save			:: cnstnd_cells,jcmin_recv,jcmax_recv
	logical,save			:: recv_flag
	logical, save 	 		:: first_time=.true.
	save CFD_box

	! Check processor is inside MD/CFD overlap zone 
	if (olap_mask(rank_world) .eq. 0) return

	if (first_time) then
		first_time = .false.
		!Number of cells to receive
		cnstnd_cells = 17	!~10% of the total domain
		jcmin_recv = jcmax_olap-1-cnstnd_cells
		jcmax_recv = jcmax_olap-1
		limits = (/ icmin_olap,icmax_olap, jcmin_recv,jcmax_recv, kcmin_olap,kcmax_olap  /)
		call setup_CFD_box(limits,CFD_box,recv_flag)
		!At first CFD step we don't have two values to extrapolate CFD velocities, set inv_dtCFD=0
		inv_dtCFD = 0.0
	else
		inv_dtCFD = 1.0/coupler_md_get_dt_cfd()
	endif
	iter_average = mod(iter-1, md_steps_per_dt_cfd)+1

	if (rank_world .eq. 1) print'(i8,a,i4,a,i4)', iter,'step', iter_average, 'of', md_steps_per_dt_cfd

	! Receive value of CFD velocities at first timestep of md_steps_per_dt_cfd
	if (iter_average .eq. 1) then
			call CPL_recv(uvw_cfd,jcmax_recv=jcmax_recv, & 
						          jcmin_recv=jcmin_recv,recv_flag=recv_flag)
	else
		!Linear extrapolation between velocity at t and t+1
	endif

	!Get average over current cell and apply constraint forces
	if (recv_flag .eq. .true.) then
		call average_over_bin
		call apply_force
	endif

contains

!===================================================================================
! Run through the particle, check if they are in the overlap region and
! find the CFD box to which the particle belongs		 

subroutine setup_CFD_box(limits,CFD_box,recv_flag)
	use coupler_module, only : iblock_realm,jblock_realm,kblock_realm, &
							   xg,yg,zg, CPL_REALM_COMM
	use coupler, only : CPL_recv
	use coupler_internal_md, only : localise,map_cfd2md_global
	implicit none

	!Limits of CFD box to receive data in
	integer,dimension(6) ,intent(in)	:: limits
	!Flag to check if limits cover current processor
	logical				 ,intent(out)	:: recv_flag
	!Returned spacial limits of CFD box to receive data
	real(kind=kind(0.d0)),dimension(6) :: CFD_box

	integer 	  		  :: pcoords(3),extents(6),portion(6)
	integer		  		  :: nclx,ncly,nclz,ncbax,ncbay,ncbaz,ierr
	integer, save 		  :: ncalls = 0
    logical, save 		  :: firsttime=.true.
	real(kind=kind(0.d0)),dimension(3) :: xyzmin,xyzmax

	! Get total number of CFD cells on each processor
	pcoords= (/ iblock_realm,jblock_realm,kblock_realm   /)
	call CPL_proc_extents(pcoords,md_realm,extents)
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1

	!Allocate CFD received box
	allocate(uvw_cfd(3,nclx,ncly,nclz))
	uvw_cfd = 0.d0

	!Get limits of constraint region in which to receive data
	call CPL_proc_portion(pcoords,md_realm,limits,portion)

	! Get physical extents of received region on MD processor
	!print'(a,19i4)', 'portion',rank_world, portion,limits,extents
	if (all(portion .ne. VOID)) then
		recv_flag = .true.
		xyzmin(1) = xg(portion(1)  ,1); xyzmax(1) = xg(portion(2)+1,1)
		xyzmin(2) = yg(1,portion(3)  ); xyzmax(2) = yg(1,portion(4)+1)
		xyzmin(3) = zg(  portion(5)  );	xyzmax(3) = zg(  portion(6)+1)

		!Map to local MD processor
		xyzmin = localise(map_cfd2md_global(xyzmin))
		xyzmax = localise(map_cfd2md_global(xyzmax))

		!Store in return array
		CFD_box(1) = xyzmin(1); CFD_box(2) = xyzmax(1)
		CFD_box(3) = xyzmin(2); CFD_box(4) = xyzmax(2)
		CFD_box(5) = xyzmin(3); CFD_box(6) = xyzmax(3)

		!Allocate MD averaging box to size of received region on MD processor
		ncbax = portion(2)-portion(1)+1
		ncbay = portion(4)-portion(3)+1
		ncbaz = portion(6)-portion(5)+1

		!print'(a,4i4,6f9.3)', 'Box average',rank_world, ncbax,ncbay,ncbaz,xmin,xmax,ymin,ymax,zmin,zmax
        allocate(box_average(ncbax,ncbay,ncbaz))
	else
		recv_flag = .false.
	endif

end subroutine setup_CFD_box

!=============================================================================
! Average molecules in overlap region to obtain values for 
! constrained dynamics algorithms
!-----------------------------------------------------------------------------

subroutine average_over_bin
	use computational_constants_MD, only : nhb
	use arrays_MD, only : r, v, a
	use coupler_module, only : dx,dy,dz,CPL_OLAP_COMM
	implicit none

	integer	:: ib,jb,kb,n

	!Zero box averages
	do kb = 1, ubound(box_average,dim=3)
	do jb = 1, ubound(box_average,dim=2)
	do ib = 1, ubound(box_average,dim=1)
		box_average(ib,jb,kb)%np   = 0
		box_average(ib,jb,kb)%v(:) = 0.0d0
		box_average(ib,jb,kb)%a(:) = 0.0d0
	enddo
	enddo
	enddo

	!find the maximum number of molecules and allocate a list array	   
	np_overlap = 0 ! number of particles in overlapping reg
    allocate(list(4,np))

	do n = 1,np

		if ( r(2,n) .gt. CFD_box(3) .and. r(2,n) .lt. CFD_box(4)) then

			ib = ceiling((r(1,n)+halfdomain(1))/dx)
			jb = ceiling((r(2,n)-CFD_box(3)   )/dy)
			kb = ceiling((r(3,n)+halfdomain(3))/dz)

			!print'(8i5,5f10.5)', rank_world,ceiling((r(2,n)+halfdomain(2))/dy),jb,jcmin_recv, & 
			!						jcmax_recv,jcmin_olap,jcmax_olap,cnstnd_cells,r(2,n),      & 
			!						CFD_box(3),CFD_box(4),dy,halfdomain(2)

			!Exlude out of domain molecules
			if (ib.lt.1 .or. ib.gt.size(box_average,1)) cycle
			if (kb.lt.1 .or. kb.gt.size(box_average,3)) cycle

			!print'(a,i4,5f8.3,5i4,l)', 'COUPLED AVERAGE',rank_world, CFD_box(3), CFD_box(4), halfdomain, &
			! 		np_overlap,n, ib,jb,kb,  (r(1,n) >= -halfdomain(1)  .and. r(1,n) < halfdomain(1) .and. &
            ! 					   r(3,n) >= -halfdomain(3)  .and. r(3,n) < halfdomain(3))
           
			np_overlap = np_overlap + 1
			list(1:4, np_overlap) = (/ n, ib, jb, kb /)

			box_average(ib,jb,kb)%np   =  box_average(ib,jb,kb)%np   + 1
			box_average(ib,jb,kb)%v(:) =  box_average(ib,jb,kb)%v(:) + v(:,n)
			box_average(ib,jb,kb)%a(:) =  box_average(ib,jb,kb)%a(:) + a(:,n)
		endif
	enddo

end subroutine average_over_bin

!=============================================================================
! Apply force to molecules in overlap region
!-----------------------------------------------------------------------------

subroutine apply_force
	use arrays_MD, only : a
	implicit none

	integer ib, jb, kb, i, ip, n
	real(kind=kind(0.d0)) alpha(3), u_cfd_t_plus_dt(3), inv_dtMD, acfd

	! set the continnum constraints for the particle in the bin
	! speed extrapolation add all up
	inv_dtMD =1.d0/delta_t

	do i = 1, np_overlap
		ip = list(1,i)
		ib = list(2,i)
		jb = list(3,i)
		kb = list(4,i)

		!print'(a,5i8,3f10.5)','Molecule in constraint', i,ib,jb,kb,ip,r(ip,:)

		n = box_average(ib,jb,kb)%np
		if ( n .eq. 0 ) stop "This cycle statement makes NO SENSE!! (in md_coupler_socket_apply_force)"

		! using the following extrapolation formula for continuum velocity
		! y = (y2-y1)/(x2-x1) * (x-x2) +y2
        !alpha(1) = inv_dtCFD*(uvw_cfd(1,ib,1,kb) - &
        !        			  uvw_cfd(1,ib,1,kb))

		!u_cfd_t_plus_dt(1) = alpha(1) * (iter_average + 1)*delta_t + uvw_cfd(1,ib,1,kb) 
		if (uvw_cfd(1,ib,jb+jcmin_recv-1,kb) .eq. 0.00) then
			print*,rank_world,ib,jb+jcmin_recv-1,kb, uvw_cfd(1,ib,jb+jcmin_recv-1,kb)
		endif

		acfd =	- box_average(ib,jb,kb)%a(1) / n - inv_dtMD * & 
				( box_average(ib,jb,kb)%v(1) / n - uvw_cfd(1,ib,jb+jcmin_recv-1,kb) )
		!if (ib .eq. 8 .and. kb .eq. 5) then
		!	print'(a,2i5,i10,4i4,5f9.4)', 'FORCE OUT', rank_world, i, np_overlap, box_average(ib,jb,kb)%np, & 
		!							   			   ib,jb,kb,box_average(ib,jb,kb)%a(1),box_average(ib,jb,kb)%v(1), & 
		!									       uvw_cfd(1,ib,size(uvw_cfd,3)-1,kb), a(1,ip), acfd
		!endif
		a(1,ip) = a(1,ip) + acfd

	enddo

end subroutine apply_force






!=============================================================================
! Average molecules using cell lists in overlap region to obtain values for 
! constrained dynamics algorithms
!**************************************
!The cell based version does not work 
!*************************************
!-----------------------------------------------------------------------------

subroutine average_over_bin_cells
	use arrays_MD, only : r, v, a
	implicit none

	integer		:: n,icell,jcell,kcell,js,je
	integer		:: cellnp, molno
	type(node), pointer :: current => null(),old => null()

	!Zero box averages
	do kcell = 1, ubound(box_average,dim=3)
	do jcell = 1, ubound(box_average,dim=2)
	do icell = 1, ubound(box_average,dim=1)
		box_average(icell,jcell,kcell)%np   = 0
		box_average(icell,jcell,kcell)%v(:) = 0.0d0
		box_average(icell,jcell,kcell)%a(:) = 0.0d0
	enddo
	enddo
	enddo

	! get the range of j cell index in y direction
	js = min(ncells(2),ceiling((CFD_box(3)+halfdomain(2))/cellsidelength(2))) + nh
	je = min(ncells(2),ceiling((CFD_box(4)+halfdomain(2))/cellsidelength(2))) + nh

	!find the maximum number of molecules and allocate a list array	   
	np_overlap = 0 ! number of particles in overlapping reg
    allocate(list(4,np))

	do kcell = nh+1,ncells(3)+1+nh
	do jcell = js  ,   je
	do icell = nh+1,ncells(1)+1+nh

		cellnp = cell%cellnp(icell,jcell,kcell)
		old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list

		!Calculate averages for bin
		do n = 1, cellnp    ! Loop over all particles
			molno = old%molno 	 !Number of molecule
			!Add streamwise velocity & acceleration to current bin
			box_average(icell,jcell,kcell)%np   = box_average(icell,jcell,kcell)%np   + 1
			box_average(icell,jcell,kcell)%v(:) = box_average(icell,jcell,kcell)%v(:) + v(molno,1)
			box_average(icell,jcell,kcell)%a(:) = box_average(icell,jcell,kcell)%a(:) + a(molno,1) 

			box_average(icell,jcell,kcell)%np   =  & 
							box_average(icell,jcell,kcell)%np   + 1
			box_average(icell,jcell,kcell)%v(:) =  & 
							box_average(icell,jcell,kcell)%v(:) + v(1,molno) !Add streamwise velocity to current bin
			box_average(icell,jcell,kcell)%a(:) =  & 
							box_average(icell,jcell,kcell)%a(:) + a(1,molno) !Add acceleration to current bin

			current => old
			old => current%next 
		enddo

	enddo
	enddo
	enddo

end subroutine average_over_bin_cells

end subroutine socket_apply_continuum_forces

!=============================================================================
! Apply force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
! Adapted serial version written by ES including cells and (originally) fully verified
!-----------------------------------------------------------------------------

subroutine socket_CPL_apply_continuum_forces(iter)
	use coupler_input_data, only : md_steps_per_dt_cfd
	implicit none
	
	integer, intent(in) :: iter

	integer :: iter_average,  average_period
	real(kind(0.d0)) :: delta_t_CFD
	logical, save :: first_time=.true.
	save  average_period

	!Setup arrays on first call
    if (first_time) then 
	    first_time	= .false.
	    average_period = coupler_md_get_average_period() 	! collection interval in the average cycle
    endif
    iter_average = mod(iter-1, md_steps_per_dt_cfd)+1			! current step

	!Receive results from CFD at exchange times
    if  (iter_average .eq. md_steps_per_dt_cfd) then
	    call receive_CFD_velocity
	else
		call interpolate_CFD_velocity
	endif

	!Apply the force to the molecular dynamics region
    if ( mod(iter_average,average_period) .eq. 0 ) then
		call average_over_bin
	    call apply_force
	endif


contains
!------------------------------------------
subroutine CFD_cells_to_MD_compute_cells(cfdis,cfdie,cfdjs,cfdje,cfdks,cfdke, & 
											  mdis, mdie, mdjs, mdje, mdks, mdke)
	implicit none
	integer,intent(in)		:: cfdis,cfdie,cfdjs,cfdje,cfdks,cfdke
	integer,intent(out)		:: mdis,mdie,mdjs,mdje,mdks,mdke

	mdis = cfdis
	mdie = cfdie
	mdjs = cfdjs 
	mdje = cfdje
	mdks = cfdks 
	mdke = cfdke 

end subroutine CFD_cells_to_MD_compute_cells
!-----------------------------------------------
subroutine receive_CFD_velocity
	implicit none
end subroutine receive_CFD_velocity
!-----------------------------------------------
subroutine interpolate_CFD_velocity
	implicit none
end subroutine interpolate_CFD_velocity
!-----------------------------------------------
subroutine average_over_bin
	implicit none
end subroutine average_over_bin
!-----------------------------------------------
subroutine apply_force
	implicit none
end subroutine apply_force
!-----------------------------------------------
end subroutine socket_CPL_apply_continuum_forces




subroutine socket_apply_continuum_forces_ES(iter)
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

	allocate(u_continuum(icPmin_md(iblock_realm):icPmax_md(iblock_realm), & 
						 jcPmin_md(jblock_realm):jcPmax_md(jblock_realm), & 
						 kcPmin_md(kblock_realm):kcPmax_md(kblock_realm)))
	u_continuum = 1.d0

	do ii=icPmin_md(iblock_realm),icPmax_md(iblock_realm)
	do jj=jcPmin_md(jblock_realm),jcPmax_md(jblock_realm)
	do kk=kcPmin_md(kblock_realm),kcPmax_md(kblock_realm)

		! For each continuum cell get MD cells to average over
		call CFD_cells_to_MD_compute_cells(ii,jj,kk,ibmin_md,ibmax_md,jbmin_md,jbmax_md,kbmin_md,kbmax_md)

		!call CFD_cells_to_MD_compute_cells(ii,jj,kk,ibmin_md,ibmax_md,jbmin_md,jbmax_md,kbmin_md,kbmax_md)
		!Choose a cube in the centre of the domain
		!ibmin_md=ceiling(ncells(1)/2.d0)-1; ibmax_md=ceiling(ncells(1)/2.d0)+1
		!jbmin_md=ceiling(ncells(2)/2.d0)-1; jbmax_md=ceiling(ncells(2)/2.d0)+1
		!kbmin_md=ceiling(ncells(3)/2.d0)-1; kbmax_md=ceiling(ncells(3)/2.d0)+1

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

	enddo
	enddo
	enddo

end subroutine socket_apply_continuum_forces_ES

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





!----------------------------------------------------------------
! Attempt to apply continuum forces based on Flekkoy 

subroutine apply_continuum_forces_flekkoy(iter)
	use computational_constants_MD, only : delta_t,nh,halfdomain,ncells,cellsidelength,initialstep,Nsteps
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	use coupler_module, only : md_steps_per_dt_cfd
	implicit none

	integer, intent(in) 				:: iter ! iteration step, it assumes that each MD average starts from iter = 1

	integer         					:: js,je,n, molno, cellnp
	integer         					:: cbin, iter_average
	integer								:: icell, jcell, kcell
	double precision					:: weight, cfd_stress

	type(node), pointer 	        	:: old, current

	iter_average = mod(iter-1, md_steps_per_dt_cfd)+1

	! Receive value of CFD velocities at first timestep of md_steps_per_dt_cfd
	if (iter_average .eq. 1) then
			call CPL_recv(uvw_cfd,jcmax_recv=jcmax_recv, & 
						          jcmin_recv=jcmin_recv,recv_flag=recv_flag)
	else

	!Apply force to top three bins in y
	!ASSUME Cell same size as bins and one continuum cell is two MD cells
	do jcell= js,je	!Loop through all overlap y cells
 		cbin = jcell - js+1					!Local no. of overlap cell from 1 to overlap
    	do kcell = nh+1,ncells(3)+1+nh !Loop through all x cells
		do icell = nh+1,ncells(1)+1+nh !Loop through all z cells

			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list
			
			!Apply coupling force as Nie, Chen and Robbins (2004), using
			!Linear extrapolation of velocity
			do n = 1, cellnp    ! Loop over all particles
				molno = old%molno !Number of molecule

				!weight = weight(r(:,n))
				a(1,molno)= a(1,molno) - weight * cfd_stress

				current => old
				old => current%next 

			enddo

		enddo
		enddo
	enddo

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine apply_continuum_forces_flekkoy






subroutine apply_continuum_forces_ES(iter)
	use computational_constants_MD, only : delta_t,nh,halfdomain,ncells,cellsidelength,initialstep,Nsteps
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	implicit none

	integer, intent(in) 				:: iter ! iteration step, it assumes that each MD average starts from iter = 1

    type(cfd_grid_info)                 :: cfd_box
	integer								:: i,j,k,js,je,ib,jb,kb,nib,njb,nkb,ip,m,np_overlap
	integer         					:: n, molno, ixyz, cellnp
	integer         					:: cbin, averagecount
	integer								:: icell, jcell, kcell
	integer								:: isummol
	double precision					:: inv_dtCFD,ymin,ymax,xmin,xmax,zmin,zmax,dx_cfd,dy_cfd,dz_cfd
	double precision					:: isumvel, isumacc
	double precision					:: t_fract, average
	double precision, dimension(4)		:: continuum_u
    integer, save                       :: ncalls = 0, itm1=1, itm2=2
	logical, save						:: firsttime, overlap
	type(node), pointer 	        	:: old, current

	! run through the particle, check if they are in the overlap region
	! find the CFD box to which the particle belongs		  
	! attention to the particle that have left the domain boundaries 
	!call setup_CFD_box(iter,xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD,itm1,itm2)

	! get the range of j cell index in y direction
	js = min(ncells(2),ceiling(ymin/cellsidelength(2))) + nh
	je = min(ncells(2),ceiling(ymax/cellsidelength(2))) + nh

	!Fraction of continuum timestep which has passed
	t_fract = dble((iter - initialstep)) / dble((Nsteps-initialstep))

	!Linear extrapolation between velocity at t and t+1 and save in continuum_u
	!Taken for cell at nx/2 and top domain cell (ny + 1) to (ny + 1) - overlap
	!N.B. domain runs from 2 to ny + 1 due to halos
	!continuum_u(1) = uc_t_minus_1(nint(nx/2.d0),3)*(1.d0-t_fract)  &
	!	 	         + uc(nint(nx/2.d0),3)*      t_fract
	!continuum_u(2) = uc_t_minus_1(nint(nx/2.d0),3)*(1.d0-t_fract)  &
	!	 	         + uc(nint(nx/2.d0),3)*      t_fract
	!continuum_u(3) = uc_t_minus_1(nint(nx/2.d0),4)*(1.d0-t_fract)  &
	!	 	         + uc(nint(nx/2.d0),4)*      t_fract
	!continuum_u(4) = uc_t_minus_1(nint(nx/2.d0),4)*(1.d0-t_fract)  &
	!	 	         + uc(nint(nx/2.d0),4)*      t_fract

	average = 0.d0
	averagecount = 0

	!Apply force to top three bins in y
	!ASSUME Cell same size as bins and one continuum cell is two MD cells
	!do jcell= (ncells(2)+1)-3,(ncells(2)+1)			!Loop through 4 y cells in controlled region
	do jcell= js,je
 		cbin = jcell - js+1					!Local no. of overlap cell from 1 to overlap
 		!cbin = jcell - (ncells(2)+1)+4		!Local no. of overlap cell from 1 to overlap

		!print'(3i8,4f20.7)', iter, jcell, cbin, continuum_u

		!Reset acceleration and velocity sums
		isummol = 0
		isumvel = 0.d0
		isumacc = 0.d0

		!do icell=2 , ncells(1)+1	!Loop through all x cells
		!do kcell=2 , ncells(3)+1	!Loop through all z cells
    	do kcell = nh+1,ncells(3)+1+nh
		do icell = nh+1,ncells(1)+1+nh

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

		!Get average velocity and acceleration in bin
		if (isummol .ne. 0) then
			isumacc = isumacc/real(isummol,kind(0.d0))
		 	isumvel = isumvel/real(isummol,kind(0.d0))
		endif

		!print*, isumacc, isumvel, isummol

		!do icell=2 , ncells(1)+1	!Loop through all x cells
		!do kcell=2 , ncells(3)+1	!Loop through all z cells
    	do kcell = nh+1,ncells(3)+1+nh
		do icell = nh+1,ncells(1)+1+nh

			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list
			
			!Apply coupling force as Nie, Chen and Robbins (2004), using
			!Linear extrapolation of velocity
			do n = 1, cellnp    ! Loop over all particles
				molno = old%molno !Number of molecule

				a(1,molno)= a(1,molno) - isumacc   &
					    -(isumvel-uvw_cfd(1,icell,1,kcell))/delta_t

				current => old
				old => current%next 

				!if (molno .gt. np) stop "Force applied to halo molecules"

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

end subroutine apply_continuum_forces_ES

!--------------------------------------------------------------------------------------
! Apply force to match velocity to continuum using CV formulation of Nie et al 2004
! which results in matching of force and flux on a CV.

subroutine simulation_apply_continuum_forces_CV(iter)
	use computational_constants_MD, only : delta_t,nh,domain,halfdomain,ncells,cellsidelength,initialstep,Nsteps
	use calculated_properties_MD, only : nbins
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	implicit none

	integer, intent(in) 				:: iter ! iteration step, it assumes that each MD average starts from iter = 1

	integer								:: i,j,k,js,je,ib,jb,kb,nib,njb,nkb,ip,m,np_overlap
	integer         					:: n, molno, ixyz, cellnp
	integer         					:: cbin, averagecount
	integer								:: icell, jcell, kcell
	integer								:: isummol
	integer, save 						:: ncalls = 0, itm1=1, itm2=2
	logical, save						:: firsttime, overlap
	double precision					:: inv_dtCFD,ymin,ymax,xmin,xmax,zmin,zmax,dx_cfd,dy_cfd,dz_cfd
	double precision					:: isumflux, isumforce, isumvel
	double precision					:: t_fract, average
	double precision, dimension(3)		:: ri
	double precision, dimension(4)		:: continuum_res, continuum_Fs
	double precision, dimension(4)		:: continuum_u
	double precision, dimension(:,:,:,:,:),allocatable 	:: uvw_cfd

	!logical								:: overlap
	type(node), pointer 	        	:: old, current
    save nib,njb,nkb,xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD

	! run through the particle, check if they are in the overlap region
	! find the CFD box to which the particle belongs		  
	! attention to the particle that have left the domain boundaries 
	!call setup_CFD_box(iter,xmin,xmax,ymin,ymax,zmin,zmax,inv_dtCFD,itm1,itm2)

	! get the range of j cell index in y direction
	js = min(ncells(2),ceiling(ymin/cellsidelength(2))) + nh
	je = min(ncells(2),ceiling(ymax/cellsidelength(2))) + nh

	!Linear extrapolation between force at t and t+1 and save in continuum_F
	!Taken for cell at nx/2 and top domain cell (ny + 1) to (ny + 1)-overlap
	!N.B. domain runs from 2 to ny + 1 due to halos
	!continuum_res(1) = xresidual_t_minus_1(nint(nx/2.d0),3)*(1.d0-t_fract)  &
	!	 	           + xresidual(nint(nx/2.d0),3)*      t_fract
	!continuum_res(2) = xresidual_t_minus_1(nint(nx/2.d0),3)*(1.d0-t_fract)  &
	!	 	           + xresidual(nint(nx/2.d0),3)*      t_fract
	!continuum_res(3) = xresidual_t_minus_1(nint(nx/2.d0),4)*(1.d0-t_fract)  &
	!	 	           + xresidual(nint(nx/2.d0),4)*      t_fract
	!continuum_res(4) = xresidual_t_minus_1(nint(nx/2.d0),4)*(1.d0-t_fract)  &
	!	 	           + xresidual(nint(nx/2.d0),4)*      t_fract

	!Fraction of continuum timestep which has passed
	t_fract = dble((iter - initialstep)) / dble((Nsteps-initialstep))

	!Apply force to top three bins in y
	!ASSUME Cell same size as bins and one continuum cell is two MD cells
	do jcell= js , je				!Loop through y cells in controlled region
 	cbin = jcell - je					!Local no. of overlap cell from 1 to overlap
	!cbin = jcell - (ncells(2)+1)+4		!Local no. of overlap cell from 1 to overlap
	do icell=2 , ncells(1)+1		!Loop through all x cells
	do kcell=2 , ncells(3)+1		!Loop through all z cells

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

		!All constraint terms are multiplied by { (m_j \theta_j)/(\sum_i^N m_i^2 \theta_i^2) }
		!Which for unit mass reduces to { 1/(sum inside volume) }
		if (isummol .ne. 0) then
			isumflux     = isumflux     / real(isummol,kind(0.d0))
		 	isumforce    = isumforce    / real(isummol,kind(0.d0))
			continuum_Fs = continuum_res/ real(isummol,kind(0.d0))
			if(icell .eq. 5 .and. kcell .eq. 5 .and. jcell .eq. 11 ) then 
				!print'(a,i8,4f10.5)','bin and forces',cbin,continuum_Fs
				!write(99,'(3i8,4f18.9)')iter,jcell,isummol, isumflux, isumforce, isumvel,continuum_Fs(cbin)
			endif
		endif

		old => cell%head(icell,jcell,kcell)%point 	!Reset old to first molecule in list
		
		!Apply coupling force for CV form of Nie, Chen and Robbins (2004)
		do n = 1, cellnp    ! Loop over all particles
			molno = old%molno !Number of molecule

			a(1,molno) = a(1,molno) + ((isumflux - isumforce) + continuum_Fs(cbin))

			current => old
			old => current%next 

			!if (molno .gt. np) stop "Force applied to halo molecules"
		enddo

	enddo
	enddo
	enddo

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine simulation_apply_continuum_forces_CV

!===================================================================================
! Flux of molecules over bin surface requires all molecules in current cell and 
! surrounding cells to be checked

subroutine compute_bin_surface_flux(icell,jcell,kcell,isumflux)
	use computational_constants_MD, only : delta_t,nh,domain,halfdomain,ncells,cellsidelength,initialstep,Nsteps
	use calculated_properties_MD, only : nbins
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	implicit none

	integer,intent(in)					:: icell,jcell,kcell
	double precision,intent(out)		:: isumflux

	integer								:: n,molno
	integer								:: icellshift,jcellshift,kcellshift,adjacentcellnp
	integer,dimension(3)				:: ibin1, ibin2
	double precision,dimension(3)		:: Fbinsize,ri1,ri2,ri12,Fsurface,velvec
	type(node), pointer		 	        :: old, current

	isumflux = 0.d0

	!Determine bin size
	Fbinsize(:) = domain(:) / nbins(:)

	!Calculate bin surface Forces
	do kcellshift = -1,1
	do jcellshift = -1,1
	do icellshift = -1,1
		old =>  cell%head(icell+icellshift, & 
				  jcell+jcellshift, & 
				  kcell+kcellshift)%point
		adjacentcellnp = cell%cellnp(icell+icellshift, & 
					     jcell+jcellshift, & 
					     kcell+kcellshift)

		do n = 1,adjacentcellnp          !Step through all adjacent cells' molecules

			molno = old%molno			!Number of molecule

			velvec(:) = v(:,molno) + delta_t *a(:,molno) 	!Velocity at t calculated from acceleration
			ri1(:)    = r(:,molno) + delta_t * velvec(:)	!Position at t calculated from velocity
			ri2(:)    = r(:,molno)				!Molecule i at time t-dt

			current => old
			old => current%next    !Use pointer in datatype to obtain next item in list

			!Check if molecule crosses surface, if it has, 
			!add to surface flux bin total

			Fsurface = 0.d0
			! *********************************************************************************
			!Calculate flux over surface only if molecule is entering/leaving bin of interest
			!ibin1(:) = ceiling((ri1(:)+halfdomain(:))/Fbinsize(:))+1	!Establish current bin
			!ibin2(:) = ceiling((ri2(:)+halfdomain(:))/Fbinsize(:))+1 !Establish current bin

			!if (icell.eq.ibin1(1).and.jcell.eq.ibin1(2).and.kcell.eq.ibin1(3).or.&
			!    icell.eq.ibin2(1).and.jcell.eq.ibin2(2).and.kcell.eq.ibin2(3)) then
			!	call get_Fsurface(icell,jcell,kcell,molno,molno,velvec,ri1,ri2,Fsurface)
			!	if(icell .eq. 5 .and. kcell .eq. 5 .and. jcell .eq. 11 ) then 
			!!		if (any(Fsurface .ne. 0.d0)) print'(2i8,6f10.5)', iter,molno, & 
			!!			dble(ibin1(:)-ibin2(:)),velvec
			!!	endif
			!endif
			! *********************************************************************************
			call get_Flux(icell,jcell,kcell,molno,molno,velvec,ri1,ri2,Fsurface)

			!if(icell .eq. 5 .and. kcell .eq. 5 ) then 
			!	if (any(Fsurface .ne. 0.d0)) print'(5i8,6f10.5)', iter,icell,jcell,kcell,molno,velvec,Fsurface
			!endif
			isumflux = isumflux + Fsurface(1)

		enddo
	enddo
	enddo
	enddo

contains

	!===================================================================================
	!Flux over the surface of a bin

	subroutine get_Flux(icell,jcell,kcell,molnoi,molnoj,velvect,ri1,ri2,Fsurface)
		use librarymod, only : heaviside
		implicit none

		integer,intent(in)							:: icell,jcell,kcell,molnoi,molnoj
		double precision,dimension(3),intent(in)	:: ri1, ri2, velvect
		double precision,dimension(3),intent(out)	:: Fsurface

		integer										:: ixyz,jxyz,kxyz,i,j,k,n
		integer										:: planeno
		integer										:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
		integer		,dimension(1)					:: imaxloc
		integer		,dimension(3)					:: ibin1,ibin2,cbin
		double precision							:: crosstime,crossplane,rplane,shift
		double precision,dimension(3)				:: mbinsize,crossface
		double precision,dimension(3)				:: ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
		double precision,dimension(3,3,3,3)			:: Fsurfacebins

		!CV momentum flux
		!Determine bin size
		mbinsize(:) = domain(:) / nbins(:)

		ri12   = ri1 - ri2		!Molecule i trajectory between t-dt and t
		where (ri12 .eq. 0.d0) ri12 = 0.000001d0

		!Assign to bins before and after using integer division
		ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:))+1
		ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:))+1

		!Replace Signum function with this functions which gives a
		!check for plane crossing and the correct sign 
		crossface(:) =  ibin1(:) - ibin2(:)

		!if(icell .eq. 5 .and. kcell .eq. 5 ) then 
		!	if (any(crossface .eq. 1)) print*, ibin1(:),ibin2(:)
		!endif

		if (sum(abs(crossface(:))) .ne. 0) then

			Fsurfacebins = 0.d0

			do i = ibin1(1),ibin2(1),sign(1,ibin2(1)-ibin1(1))
			do j = ibin1(2),ibin2(2),sign(1,ibin2(2)-ibin1(2))
			do k = ibin1(3),ibin2(3),sign(1,ibin2(3)-ibin1(3))

				cbin(1) = i; cbin(2) = j; cbin(3) = k

				bintop(:) = (cbin(:)-1)*mbinsize(:)-halfdomain(:)
				binbot(:) = (cbin(:)-2)*mbinsize(:)-halfdomain(:)

				!Calculate the plane intersect of trajectory with surfaces of the cube
				Pxt=(/ 							  bintop(1),		 & 
						ri1(2)+(ri12(2)/ri12(1))*(bintop(1)-ri1(1)), & 
						ri1(3)+(ri12(3)/ri12(1))*(bintop(1)-ri1(1))  	/)
				Pxb=(/ 							  binbot(1),		 & 
						ri1(2)+(ri12(2)/ri12(1))*(binbot(1)-ri1(1)), & 
						ri1(3)+(ri12(3)/ri12(1))*(binbot(1)-ri1(1))  	/)
				Pyt=(/	ri1(1)+(ri12(1)/ri12(2))*(bintop(2)-ri1(2)), & 
												  bintop(2),		 & 
						ri1(3)+(ri12(3)/ri12(2))*(bintop(2)-ri1(2))  	/)
				Pyb=(/	ri1(1)+(ri12(1)/ri12(2))*(binbot(2)-ri1(2)), &
												  binbot(2),		 & 
						ri1(3)+(ri12(3)/ri12(2))*(binbot(2)-ri1(2))  	/)
				Pzt=(/	ri1(1)+(ri12(1)/ri12(3))*(bintop(3)-ri1(3)), & 
						ri1(2)+(ri12(2)/ri12(3))*(bintop(3)-ri1(3)), &
											 	  bintop(3) 			/)
				Pzb=(/	ri1(1)+(ri12(1)/ri12(3))*(binbot(3)-ri1(3)), &
						ri1(2)+(ri12(2)/ri12(3))*(binbot(3)-ri1(3)), & 
												  binbot(3) 			/)

				onfacexb =0.5d0*(sign(1.d0,binbot(1) - ri2(1)) 	 & 
					       - sign(1.d0,binbot(1) - ri1(1)))* &
							(heaviside(bintop(2) - Pxb(2)) 	 &
					       - heaviside(binbot(2) - Pxb(2)))* &
							(heaviside(bintop(3) - Pxb(3)) 	 &
					       - heaviside(binbot(3) - Pxb(3)))
				onfaceyb =0.5d0*(sign(1.d0,binbot(2) - ri2(2))   &
					       - sign(1.d0,binbot(2) - ri1(2)))* &
							(heaviside(bintop(1) - Pyb(1))   &
					       - heaviside(binbot(1) - Pyb(1)))* &
							(heaviside(bintop(3) - Pyb(3))   &
					       - heaviside(binbot(3) - Pyb(3)))
				onfacezb =0.5d0*(sign(1.d0,binbot(3) - ri2(3))   &
					       - sign(1.d0,binbot(3) - ri1(3)))* &
							(heaviside(bintop(1) - Pzb(1))   &
					       - heaviside(binbot(1) - Pzb(1)))* &
							(heaviside(bintop(2) - Pzb(2))   &
					       - heaviside(binbot(2) - Pzb(2)))

				onfacext =0.5d0*(sign(1.d0,bintop(1) - ri2(1))   &
					       - sign(1.d0,bintop(1) - ri1(1)))* &
							(heaviside(bintop(2) - Pxt(2))   &
					       - heaviside(binbot(2) - Pxt(2)))* &
			            	(heaviside(bintop(3) - Pxt(3))   &
					       - heaviside(binbot(3) - Pxt(3)))
				onfaceyt =0.5d0*(sign(1.d0,bintop(2) - ri2(2))   &
					       - sign(1.d0,bintop(2) - ri1(2)))* &
							(heaviside(bintop(1) - Pyt(1))   &
					       - heaviside(binbot(1) - Pyt(1)))* &
							(heaviside(bintop(3) - Pyt(3))   &
					       - heaviside(binbot(3) - Pyt(3)))
				onfacezt =0.5d0*(sign(1.d0,bintop(3) - ri2(3))   &
					       - sign(1.d0,bintop(3) - ri1(3)))* &
							(heaviside(bintop(1) - Pzt(1))   &
					       - heaviside(binbot(1) - Pzt(1)))* &
							(heaviside(bintop(2) - Pzt(2))   &
					       - heaviside(binbot(2) - Pzt(2)))


				!Add Momentum flux over face
				Fsurfacebins(modulo(cbin(1),3)+1, 	& 
					     modulo(cbin(2),3)+1, 	& 
					     modulo(cbin(3),3)+1,:) = 	& 
						 Fsurfacebins(modulo(cbin(1),3)+1, 	& 
							      modulo(cbin(2),3)+1, 	& 
							      modulo(cbin(3),3)+1,:) 	& 
							 		- velvect(:)*dble(onfacexb - onfacext) & 
							 		- velvect(:)*dble(onfaceyb - onfaceyt) & 
							 		- velvect(:)*dble(onfacezb - onfacezt)

				!if(icell .eq. 5 .and. kcell .eq. 5 ) then 
				!	print'(12i8)', cbin, ibin1,ibin2,icell,jcell,kcell
				!	print'(12i8)', modulo(cbin(:),3)+1,modulo(ibin1(:),3)+1,modulo(ibin2(:),3)+1 & 
				!	      	      ,modulo(icell,3)+1,modulo(jcell,3)+1,modulo(kcell,3)+1
				!	print'(6f10.5,6i8)', velvect(:),Fsurface,onfacexb,onfacext,onfaceyb,onfaceyt,onfacezb,onfacezt
				!endif

			enddo
			enddo
			enddo

			!if(icell .eq. 5 .and. kcell .eq. 5 ) then 
			!	print'(6f10.5,6i8)', velvect(:),Fsurface,onfacexb,onfacext,onfaceyb,onfaceyt,onfacezb,onfacezt
			!endif
			

			Fsurface(:) = Fsurfacebins( modulo(icell,3)+1, 	& 
										modulo(jcell,3)+1, 	& 
										modulo(kcell,3)+1,:)

		endif

	end subroutine get_Flux

end subroutine compute_bin_surface_flux

!===================================================================================
! Forces between current bin and surrounding bins

subroutine compute_force_surrounding_bins(icell,jcell,kcell,isumforce,Traction)
	use physical_constants_MD, only : nd, np, rcutoff2
	use computational_constants_MD, only : domain,halfdomain
	use calculated_properties_MD, only : nbins
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	use librarymod, only : heaviside
	implicit none

	integer,intent(in)					:: icell,jcell,kcell
	double precision,intent(out)		:: isumforce
	double precision,dimension(3,6),optional,intent(inout)	:: Traction

	integer								:: n,j,ixyz,molnoi,molnoj
	integer								:: icellshift,jcellshift,kcellshift,cellnp,adjacentcellnp
	integer,dimension(3)				:: ibin, jbin
	double precision					:: rij2, invrij2, accijmag
	double precision,dimension(3)		:: ri,rj,rij,fij,Fsurface
	type(node), pointer		 	        :: oldi, currenti, oldj, currentj

	!print'(a,4i8,4f10.5)', 'Before input', icell,jcell,kcell,molnoi,ri(:),isumforce

	isumforce = 0.d0

	cellnp = cell%cellnp(icell,jcell,kcell)
	oldi => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list

	!Calculate averages for bin
	do n = 1, cellnp    ! Loop over all particles

		molnoi = oldi%molno	!Number of molecule
		ri = r(:,molnoi)	!Retrieve ri

		!Calculate bin surface Forces
		do kcellshift = -1,1
		do jcellshift = -1,1
		do icellshift = -1,1
			oldj => cell%head(icell+icellshift, & 
					  jcell+jcellshift, & 
					  kcell+kcellshift)%point
			adjacentcellnp = cell%cellnp(icell+icellshift, & 
						     jcell+jcellshift, & 
						     kcell+kcellshift)

			do j = 1,adjacentcellnp          !Step through all j for each i

				molnoj = oldj%molno 	 !Number of molecule
				rj = r(:,molnoj)         !Retrieve rj

				currentj => oldj
				oldj => currentj%next    !Use pointer in datatype to obtain next item in list

				if(molnoi==molnoj) cycle !Check to prevent interaction with self

				rij2=0                   !Set rij^2 to zero
				rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j

				!rij2 = dot_product(rij)
				do ixyz=1,nd
					rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
				enddo

				if (rij2 < rcutoff2) then

					invrij2 = 1.d0/rij2                 !Invert value

					!Linear magnitude of acceleration for each molecule
					accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)

					!Get force and add to bin total
					fij = accijmag*rij(:)
					Fsurface = 0.d0
					!print*, 'FORCE', fij
					if (present(Traction)) then
						call get_Traction(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface,Traction)
					else
						call get_Fsurface(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface)
					endif
					isumforce = isumforce +  Fsurface(1)
					call get_Fsurface(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface)
					isumforce = isumforce +  Fsurface(1)

				endif
			enddo
		enddo
		enddo
		enddo

		currenti => oldi
		oldi => currenti%next !Use pointer in datatype to obtain next item in list
	enddo

contains

	!===================================================================================
	!Forces over the surface of a bin

	subroutine get_Fsurface(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface)
		implicit none

		integer,intent(in)							:: icell,jcell,kcell,molnoi,molnoj
		double precision,dimension(3),intent(in)	:: ri, rj, fij
		double precision,dimension(3),intent(out)	:: Fsurface

		integer										:: ixyz
		integer,dimension(3)						:: ibin, jbin
		double precision,dimension(3)				:: Fbinsize, bintopi, binboti, bintopj, binbotj, crossplane
	
		!Determine bin size
		Fbinsize(:) = domain(:) / nbins(:)

		!Assign to bins using integer division
		ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+1	!Establish current bin
		jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+1 	!Establish current bin

		bintopi(:) = (ibin(:)-1)*Fbinsize(:)-halfdomain(:)
		binboti(:) = (ibin(:)-2)*Fbinsize(:)-halfdomain(:)
		bintopj(:) = (jbin(:)-1)*Fbinsize(:)-halfdomain(:)
		binbotj(:) = (jbin(:)-2)*Fbinsize(:)-halfdomain(:)

		!Add for molecule i
		if(molnoi .le. np) then
			Fsurface = fij(:)* dble((heaviside(bintopi(1)-ri(1))-heaviside(binboti(1)-ri(1)))* & 
							  		(heaviside(bintopi(2)-ri(2))-heaviside(binboti(2)-ri(2)))* & 
							  		(heaviside(bintopi(3)-ri(3))-heaviside(binboti(3)-ri(3)))- & 
							  		(heaviside(bintopi(1)-rj(1))-heaviside(binboti(1)-rj(1)))* & 
							  		(heaviside(bintopi(2)-rj(2))-heaviside(binboti(2)-rj(2)))* & 
							  		(heaviside(bintopi(3)-rj(3))-heaviside(binboti(3)-rj(3))))
		endif

		!Add for molecule j
		if(molnoj .le. np) then
			Fsurface = fij(:)* dble((heaviside(bintopj(1)-ri(1))-heaviside(binbotj(1)-ri(1)))* & 
							  		(heaviside(bintopj(2)-ri(2))-heaviside(binbotj(2)-ri(2)))* & 
							  		(heaviside(bintopj(3)-ri(3))-heaviside(binbotj(3)-ri(3)))- & 
							  		(heaviside(bintopj(1)-rj(1))-heaviside(binbotj(1)-rj(1)))* & 
							  		(heaviside(bintopj(2)-rj(2))-heaviside(binbotj(2)-rj(2)))* & 
							  		(heaviside(bintopj(3)-rj(3))-heaviside(binbotj(3)-rj(3))))
		endif

	end subroutine get_Fsurface



!-----------------------------------------------------------------------------------
! Tractions on one surface of a bin

	subroutine get_Traction(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface,Traction)
		implicit none

		integer,intent(in)										:: molnoi,molnoj
		double precision,dimension(3),intent(in)				:: ri,rj,fij
		double precision,dimension(3),intent(out)				:: Fsurface
		double precision,dimension(3,6),intent(inout),optional	:: Traction

		integer									:: i,j,k,ixyz,n,tempi
		integer									:: icell,jcell,kcell
		integer									:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
		integer,dimension(3)					:: cbin, ibin, jbin
		double precision						:: binforce
		double precision,dimension(3)			:: rij,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
		double precision,dimension(3)			:: Fbinsize, bintop, binbot
		double precision,dimension(3,3,3,3,6)	:: Tractionbins
		!Calculate rij
		rij = ri - rj
		!Prevent Division by zero
		do ixyz = 1,3
			if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
		enddo

		!Determine bin size
		Fbinsize(:) = domain(:) / nbins(:)

		!Assign to bins using integer division
		ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+1	!Establish current bin
		jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+1 	!Establish current bin

		do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
		do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
		do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

			cbin(1) = i; cbin(2) = j; cbin(3) = k

			bintop(:) = (ibin(:)-1)*Fbinsize(:)-halfdomain(:)
			binbot(:) = (ibin(:)-2)*Fbinsize(:)-halfdomain(:)

			if(present(Traction)) then

				!Calculate the plane intersect of line with surfaces of the cube
				Pxt=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
				Pxb=(/ binbot(1),ri(2)+(rij(2)/rij(1))*(binbot(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(binbot(1)-ri(1))  /)
				Pyt=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
				Pyb=(/ri(1)+(rij(1)/rij(2))*(binbot(2)-ri(2)), binbot(2),ri(3)+(rij(3)/rij(2))*(binbot(2)-ri(2))  /)
				Pzt=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)
				Pzb=(/ri(1)+(rij(1)/rij(3))*(binbot(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(binbot(3)-ri(3)), binbot(3) /)

				onfacexb =   (sign(1.d0,binbot(1)- rj(1))   &
							- sign(1.d0,binbot(1)- ri(1)))* &
							 (heaviside(bintop(2)-Pxb(2))   &
							- heaviside(binbot(2)-Pxb(2)))* &
							 (heaviside(bintop(3)-Pxb(3))   & 
							- heaviside(binbot(3)-Pxb(3)))
				onfaceyb =	 (sign(1.d0,binbot(2)- rj(2))	&
							- sign(1.d0,binbot(2)- ri(2)))* &
							 (heaviside(bintop(1)-Pyb(1))	& 
							- heaviside(binbot(1)-Pyb(1)))* &
							 (heaviside(bintop(3)-Pyb(3))	&
							- heaviside(binbot(3)-Pyb(3)))
				onfacezb =	 (sign(1.d0,binbot(3)- rj(3)) 	&
							- sign(1.d0,binbot(3)- ri(3)))* &
							 (heaviside(bintop(1)-Pzb(1))	&
							- heaviside(binbot(1)-Pzb(1)))* &
							 (heaviside(bintop(2)-Pzb(2))	&
							- heaviside(binbot(2)-Pzb(2)))

				onfacext =	 (sign(1.d0,bintop(1)- rj(1))	&
							- sign(1.d0,bintop(1)- ri(1)))* &
							 (heaviside(bintop(2)-Pxt(2))	&
							- heaviside(binbot(2)-Pxt(2)))* &
			            	 (heaviside(bintop(3)-Pxt(3))	&
							- heaviside(binbot(3)-Pxt(3)))
				onfaceyt =	 (sign(1.d0,bintop(2)- rj(2))	&
							- sign(1.d0,bintop(2)- ri(2)))* &
							 (heaviside(bintop(1)-Pyt(1))	&
							- heaviside(binbot(1)-Pyt(1)))* &
							 (heaviside(bintop(3)-Pyt(3))	&
							- heaviside(binbot(3)-Pyt(3)))
				onfacezt =	 (sign(1.d0,bintop(3)- rj(3))	&
							- sign(1.d0,bintop(3)- ri(3)))* &
							 (heaviside(bintop(1)-Pzt(1))	&
							- heaviside(binbot(1)-Pzt(1)))* &
							 (heaviside(bintop(2)-Pzt(2))	&
							- heaviside(binbot(2)-Pzt(2)))

				!Prevent halo molecules from being included but include molecule which have left domain 
				!before rebuild has been triggered.
				if (molnoi .gt. np .or. molnoj .gt. np) then
					if (cbin(1) .lt. 2 .or. cbin(1) .gt. nbins(1)+1) cycle
					if (cbin(2) .lt. 2 .or. cbin(2) .gt. nbins(2)+1) cycle
					if (cbin(3) .lt. 2 .or. cbin(3) .gt. nbins(3)+1) cycle
				endif

				!Stress acting on face over volume
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,1) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,1) + fij(:)*dble(onfacexb)
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,2) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,2) + fij(:)*dble(onfaceyb)
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,3) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,3) + fij(:)*dble(onfacezb)
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,4) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,4) + fij(:)*dble(onfacext)
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,5) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,5) + fij(:)*dble(onfaceyt)
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,6) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,6) + fij(:)*dble(onfacezt)	


				!Force applied to volume
				!fsurface(:) = 0.25d0*fij(:)*dble(onfacexb - onfacext)
				!fsurface(:) = 0.25d0*fij(:)*dble(onfaceyb - onfaceyt)
				!fsurface(:) = 0.25d0*fij(:)*dble(onfacezb - onfacezt)

				!Add surface force to current bin
				!Traction(:,1) = Traction(:,1) + 0.25d0*fij(:)*dble(onfacexb)
				!Traction(:,2) = Traction(:,2) + 0.25d0*fij(:)*dble(onfaceyb)
				!Traction(:,3) = Traction(:,3) + 0.25d0*fij(:)*dble(onfacezb)
				!Traction(:,4) = Traction(:,4) + 0.25d0*fij(:)*dble(onfacext)
				!Traction(:,5) = Traction(:,5) + 0.25d0*fij(:)*dble(onfaceyt)
				!Traction(:,6) = Traction(:,6) + 0.25d0*fij(:)*dble(onfacezt)

				!Add for molecule i
				Fsurface = fij(:)* dble((heaviside(bintop(1)-ri(1))-heaviside(binbot(1)-ri(1)))* & 
					  					(heaviside(bintop(2)-ri(2))-heaviside(binbot(2)-ri(2)))* & 
					  					(heaviside(bintop(3)-ri(3))-heaviside(binbot(3)-ri(3)))- & 
					  					(heaviside(bintop(1)-rj(1))-heaviside(binbot(1)-rj(1)))* & 
					  					(heaviside(bintop(2)-rj(2))-heaviside(binbot(2)-rj(2)))* & 
					  					(heaviside(bintop(3)-rj(3))-heaviside(binbot(3)-rj(3))))

				!if (onfaceyb.ne.0.or.onfaceyt.ne.0) print'(9i8)', iter, molnoi, molnoj,np, ibin,onfaceyb,onfaceyt

			else
				!Add for molecule i
				Fsurface = fij(:)* dble((heaviside(bintop(1)-ri(1))-heaviside(binbot(1)-ri(1)))* & 
					  					(heaviside(bintop(2)-ri(2))-heaviside(binbot(2)-ri(2)))* & 
					  					(heaviside(bintop(3)-ri(3))-heaviside(binbot(3)-ri(3)))- & 
					  					(heaviside(bintop(1)-rj(1))-heaviside(binbot(1)-rj(1)))* & 
					  					(heaviside(bintop(2)-rj(2))-heaviside(binbot(2)-rj(2)))* & 
					  					(heaviside(bintop(3)-rj(3))-heaviside(binbot(3)-rj(3))))
			endif

		enddo
		enddo
		enddo

		!Take flux from central bin only
		if (present(Traction))then
			Traction(:,:) = Traction(:,:) +  Tractionbins(modulo(icell,3)+1, 	& 
							     		     			  modulo(jcell,3)+1, 	& 
			    						     			  modulo(kcell,3)+1,:,:)


			if (icell .eq. 5 .and. kcell .eq. 3) then
				print'(3i8,2f10.5)',icell,jcell,kcell, Tractionbins(modulo(icell,3)+1, 	& 
							     		     			  modulo(jcell,3)+1, 	& 
			    						     			  modulo(kcell,3)+1,1,2),& 
											 Tractionbins(modulo(icell,3)+1, 	& 
							     		     			  modulo(jcell,3)+1, 	& 
			    						     			  modulo(kcell,3)+1,1,5)
			endif
		endif
			
	end subroutine get_Traction


end subroutine compute_force_surrounding_bins



! ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬ஜ۩۞۩ஜ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬
! Test the send and recv routines from coupler

subroutine test_send_recv_MD2CFD
	use coupler_module
	use coupler
	use computational_constants_MD, only : Nsteps
	implicit none

	logical	:: send_flag,recv_flag
	integer :: ncxl,ncyl,nczl,ixyz,icell,jcell,kcell
	integer	:: jcmin_send,jcmax_send,jcmin_recv,jcmax_recv,npercell,coord(3),extents(6)
	double precision,dimension(:,:,:,:),allocatable	:: sendbuf,recvbuf

	npercell = 3
	jcmax_send=1; jcmin_send=1; 
	jcmax_recv = jcmax_send
	jcmin_recv = jcmin_send

	call CPL_Cart_coords(CPL_WORLD_COMM,rank_world,realm,3,coord,ierr)
	!print'(2a,5i8)', 'MD SIDE',realm_name(realm), rank_world, olap_mask(rank_world),coord

	if (olap_mask(rank_world) .eq. 0) return

	! Test Sending from MD to CFD							   
	if (realm .eq. md_realm) then

		if (Nsteps_cfd .ne. Nsteps) then
			call error_abort("test_send_recv_MD2CFD error - MD_STEPS_PER_DT_CFD must be 1 in COUPLER.in for this testcase")
		endif

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_olap_extents(coord,realm,extents)

		allocate(sendbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))

		!print'(2a,11i7)', 'sent size',realm_name(realm),extents,size(sendbuf),shape(sendbuf)

		! Populate dummy gatherbuf
		sendbuf = -333.d0 ! 0.d0
		do ixyz = 1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)
			sendbuf(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                                            1000*jcell + &
			                                         1000000*kcell
		end do
		end do
		end do
		end do

		call CPL_send(sendbuf,jcmax_send=jcmax_send,jcmin_send=jcmin_send,send_flag=send_flag)	

		if (send_flag .eqv. .true.) then
			do kcell=extents(5),extents(6)
			do jcell=jcmin_send,jcmax_send
			do icell=extents(1),extents(2)
			do ixyz =1,npercell
				write(4000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
				      'send MD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
				       sendbuf(ixyz,icell,jcell,kcell)
			end do
			end do
			end do
			end do
		endif

	else if (realm .eq. cfd_realm) then	 

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,realm,3,coord,ierr)
		call CPL_proc_extents(coord,realm,extents)
		!print'(2a,8i7)', 'proc extents', realm_name(realm),rank_world,rank_cart,extents
		call CPL_olap_extents(coord,realm,extents)
		!print'(2a,8i7)', 'olap extents', realm_name(realm),rank_world,rank_cart,extents

		allocate(recvbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))

		!print'(2a,11i7)', 'recv size', realm_name(realm),extents,size(recvbuf),shape(recvbuf)
		recvbuf = -444.d0
		call CPL_recv(recvbuf,jcmax_recv=jcmax_recv,jcmin_recv=jcmin_recv,recv_flag=recv_flag)

		if (recv_flag .eqv. .true.) then
			do kcell=extents(5),extents(6)
			do jcell=jcmin_recv,jcmax_recv  !extents(3),extents(4)
			do icell=extents(1),extents(2)
			do ixyz =1,npercell
					write(5000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
					      'recv CFD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
					       recvbuf(ixyz,icell,jcell,kcell)
			end do
			end do
			end do
			end do
		endif
	end if								   
	
end subroutine test_send_recv_MD2CFD


! ۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩۩
! Test Sending from MD to CFD

subroutine test_send_recv_CFD2MD
	use coupler_module
	use coupler
	use computational_constants_MD, only : Nsteps
	implicit none

	logical	:: send_flag,recv_flag
	integer	:: jcmin_send,jcmax_send,jcmin_recv,jcmax_recv
	integer :: ncxl,ncyl,nczl,ixyz,icell,jcell,kcell,npercell,coord(3),extents(6)
	double precision,dimension(:,:,:,:),allocatable	:: sendbuf,recvbuf

	npercell = 3
	jcmax_send=1; jcmin_send=1; 
	jcmax_recv = jcmax_send
	jcmin_recv = jcmin_send
	if (olap_mask(rank_world) .eq. 0) return

	! Test Sending from CFD to MD							   
	if (realm .eq. md_realm) then	

		if (Nsteps_cfd .ne. Nsteps) then
			call error_abort("test_send_recv_MD2CFD error - MD_STEPS_PER_DT_CFD must be 1 in COUPLER.in for this testcase")
		endif	   

		coord = (/iblock_realm,jblock_realm,kblock_realm /)
		call CPL_olap_extents(coord,realm,extents)

		allocate(recvbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))
		recvbuf = -444

		!print*, 'recv size', realm_name(realm),extents, size(recvbuf),shape(recvbuf)
		call CPL_recv(recvbuf,jcmax_recv=1,jcmin_recv=1,recv_flag=recv_flag)   

		if (recv_flag .eqv. .true.) then
			do kcell=extents(5),extents(6)
			do jcell=jcmin_send,jcmax_send
			do icell=extents(1),extents(2)
			do ixyz = 1,npercell
				write(11000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
			      	'recv MD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
			      	 recvbuf(ixyz,icell,jcell,kcell)
			end do
			end do
			end do
			end do
		endif

	else if (realm .eq. cfd_realm) then	   

		coord = (/iblock_realm,jblock_realm,kblock_realm /)
		call CPL_olap_extents(coord,realm,extents)
		allocate(sendbuf(npercell,extents(1):extents(2), &
		                          extents(3):extents(4), &
		                          extents(5):extents(6)))

		do ixyz =1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)
			sendbuf(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*(icell) + &
			                       			  1000*(jcell) + &
			                    			  1000000*(kcell)

		end do
		end do
		end do
		end do

		!print*, 'sent size',realm_name(realm),3*ncxl*ncyl*nczl,size(sendbuf)
		call CPL_send(sendbuf,jcmax_send=1,jcmin_send=1,send_flag=send_flag)

		do kcell=extents(5),extents(6)
		do jcell=jcmin_send,jcmax_send
		do icell=extents(1),extents(2)
		do ixyz = 1,npercell
			write(9000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
			      'send CFD(',ixyz,',',icell,',',jcell,',',kcell,') =', &
			       sendbuf(ixyz,icell,jcell,kcell)
		end do
		end do
		end do
		end do
	end if								   
	
end subroutine test_send_recv_CFD2MD

subroutine test_gather_scatter
	use coupler_module
	use coupler
	use computational_constants_MD, only : Nsteps
	implicit none

	double precision,dimension(:,:,:,:),allocatable	:: u,stress,gatheru,scatterstress
	integer :: coord(3), extents(6), gatherlims(6), scatterlims(6), npercell
	integer :: pos, ixyz, icell, jcell, kcell
	integer :: ncxl,ncyl,nczl
	integer :: i,j,k

 	if (olap_mask(rank_world).ne.1) return

	!print*, 'test_gather_scatter called on MD proc ID:', rank_realm, rank_world

	if (realm .eq. md_realm) then	

		if (Nsteps_cfd .ne. Nsteps) then
			call error_abort("test_send_recv_MD2CFD error - MD_STEPS_PER_DT_CFD must be 1 in COUPLER.in for this testcase")
		endif

		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,md_realm,3,coord,ierr)
		call CPL_proc_extents(coord,md_realm,extents)
		npercell = 3
		allocate(u(npercell,extents(1):extents(2), &
		                    extents(3):extents(4), &
		                    extents(5):extents(6)))
		allocate(stress(0,0,0,0))

		! Populate dummy gatherbuf
		pos = 1
		do ixyz = 1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)

			u(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                                      1000*jcell + &
			                                   1000000*kcell
			pos = pos + 1

		end do
		end do
		end do
		end do

	else if (realm .eq. cfd_realm) then	  
		
		call CPL_Cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,extents)
		npercell = 9
		allocate(u(0,0,0,0))
		allocate(stress(npercell,extents(1):extents(2), &
		                         extents(3):extents(4), &
		                         extents(5):extents(6)))

		! Populate dummy gatherbuf
		pos = 1
		do ixyz = 1,npercell
		do icell=extents(1),extents(2)
		do jcell=extents(3),extents(4)
		do kcell=extents(5),extents(6)

			stress(ixyz,icell,jcell,kcell) = 0.1d0*ixyz + 1*icell + &
			                                           1000*jcell + &
			                                        1000000*kcell
			pos = pos + 1

		end do
		end do
		end do
		end do

	endif


	! Allocate test arrays over local domain
	if (realm.eq.cfd_realm) then
		call CPL_cart_coords(CPL_CART_COMM,rank_cart,cfd_realm,3,coord,ierr)
		call CPL_proc_extents(coord,cfd_realm,extents)
		ncxl = extents(2) - extents(1) + 1
		ncyl = extents(4) - extents(3) + 1
		nczl = extents(6) - extents(5) + 1
		allocate(gatheru(3,ncxl,ncyl,nczl))
		gatheru = 0.d0
	else if (realm.eq.md_realm) then
		call CPL_cart_coords(CPL_CART_COMM,rank_cart,md_realm,3,coord,ierr)
		call CPL_proc_extents(coord,md_realm,extents)
		ncxl = extents(2) - extents(1) + 1
		ncyl = extents(4) - extents(3) + 1
		nczl = extents(6) - extents(5) + 1
		allocate(scatterstress(9,ncxl,ncyl,nczl))
		scatterstress = 0.d0
	end if




	!gatherlims  = (/1,1,1,1,1,1/)
	!scatterlims = (/1,1,1,1,1,1/)
	!================== PERFORM GATHER/SCATTER =============================!	
	gatherlims  = (/1,85,15,21, 3, 4/)
	scatterlims = (/1,85, 2, 9, 1, 8/)
	if (olap_mask(rank_world).eq.1) call CPL_gather(u,3,gatherlims,gatheru)
	if (olap_mask(rank_world).eq.1) call CPL_scatter(stress,9,scatterlims, &
	                                                 scatterstress)

	! Print results to file
	if (realm.eq.cfd_realm) then

		do ixyz  = 1,size(gatheru,1)
		do icell = 1,size(gatheru,2)
		do jcell = 1,size(gatheru,3)
		do kcell = 1,size(gatheru,4)

			i = icell + extents(1) - 1
			j = jcell + extents(3) - 1
			k = kcell + extents(5) - 1

			if (gatheru(ixyz,icell,jcell,kcell).lt.0.0001) then
				!write(8000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
				!	  'gatheru(',0,',',0,',',0,',',0,') =', 0.d0
			else
				write(8000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
					  'gatheru(',ixyz,',',i,',',j,',',k,') =', &
					   gatheru(ixyz,icell,jcell,kcell)
			end if

		end do	
		end do	
		end do
		end do

	else if (realm.eq.md_realm) then

		do ixyz  = 1,size(scatterstress,1)
		do icell = 1,size(scatterstress,2)
		do jcell = 1,size(scatterstress,3)
		do kcell = 1,size(scatterstress,4)

			i = icell + extents(1) - 1
			j = jcell + extents(3) - 1
			k = kcell + extents(5) - 1

			if (scatterstress(ixyz,icell,jcell,kcell).lt.0.0001) then
				!write(7000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
				!	  'scatterstress(',0,',',0,',',0,',',0,') =', 0.d0
			else
				write(7000+myid_world,'(a,i4,a,i4,a,i4,a,i4,a,f20.1)'),   &
					  'scatterstress(',ixyz,',',i,',',j,',',k,') =', &
					   scatterstress(ixyz,icell,jcell,kcell)
			end if

		end do	
		end do	
		end do
		end do
	
	end if

	!print*, 'test_gather_scatter finished on MD proc ID:', rank_realm, rank_world
	
end subroutine test_gather_scatter

! ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬


#if COUPLER_DEBUG_LA
	! dump debug data 
	subroutine write_uc(iter)
	    use coupler
	    use computational_constants_MD, only : initialstep
	    use physical_constants_MD, only : np
		use arrays_MD, only :r,v
	    implicit none
	    integer,intent(in) :: iter

	    integer :: iter_cfd, iter_average, Naverage, save_period, average_period
		logical, save :: first_time=.true.
		save  save_period, average_period, Naverage

	    if (first_time) then 
		    first_time     = .false.
	        save_period    = coupler_md_get_save_period()   
		    average_period = coupler_md_get_average_period() 	! collection interval in the average cycle
		    Naverage = coupler_md_get_md_steps_per_cfd_dt()  	! number of steps in MD average cycle
	    endif

	    iter_average = mod(iter-1, Naverage)+1			! current step
	    iter_cfd     = (iter-initialstep)/Naverage +1 	! CFD corresponding step

	    if ( mod(iter_cfd,save_period) == 0) then 
	        if (mod(iter_average,average_period) == 0 ) then
	            call coupler_uc_average_test(np,r,v,lwrite=.false.)
	        endif
	        if (iter_average == Naverage) then
	            call coupler_uc_average_test(np,r,v,lwrite=.true.)
	        endif
	    endif
	end subroutine write_uc
#endif

#endif

end module md_coupler_socket
