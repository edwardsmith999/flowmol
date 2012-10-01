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

	real(kind(0.d0)), allocatable :: vel_cfd(:,:,:,:,:), vbuff(:,:,:,:)
	type(cfd_box_sum),allocatable :: box_average(:,:,:)

	integer			, dimension(:,:,:,:), allocatable :: mflux
	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: uc_bin, vc_bin, wc_bin
	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: uvwbin 

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

	call coupler_create_comm(md_realm,MD_COMM,ierr)
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

	! Setup timesteps and simulation timings
	call set_coupled_timing

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

	write(0,*) 'set_parameter_global_domain_hybrid ', globalnp, np, domain, initialunitsize

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

    if(rank_realm == 0) then
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
!Setup ratio of CFD to MD timing and total number of timesteps

subroutine set_coupled_timing
	use computational_constants_MD
	use coupler_input_data
	use coupler_module
	use coupler
	implicit none

	integer		:: naverage

	!Set number of MD timesteps per CFD using ratio of timestep or coupler value
	if(md_steps_per_dt_cfd_tag == CPL) then
		nsteps_MD = md_steps_per_dt_cfd
		nsteps_coupled = nsteps_cfd
	else 
		nsteps_MD = int(dt_cfd/dt_MD)
		nsteps_coupled = nsteps_cfd
	endif

	! fix NSTEPS for the coupled case
    nsteps_cfd = coupler_md_get_nsteps()
    naverage   = coupler_md_get_md_steps_per_cfd_dt()
	nsteps_coupled = nsteps_cfd
    
	Nsteps = initialstep + nsteps_cfd * naverage
	elapsedtime = elapsedtime + nsteps_cfd * naverage * delta_t

	if (rank_realm .eq. 0) then 
		write(*,'(2(a,/),a,i7,a,i7,/a,/a,i8,a/,a,f10.5,/a)') &
				"*********************************************************************", 	&
 				"WARNING - WARNING - WARNING - WARNING - WARNING - WARNING - WARNING  ", 	&
				" Current input timesteps from MD", nsteps_md, "and CFD", nsteps_cfd   ,	&
				" this is a coupled run which resets the number of extrasteps to:     ", 	&
				"								     ", 									&
											nsteps_coupled*naverage,  							&
				"								     ", 									& 
				" The elapsed time was changed accordingly to: ", elapsedtime, " s    ", 	&
				" The value of NSTEPS parameter form input file was discarded.	", 			&
				"*********************************************************************"   
	endif 

end subroutine set_coupled_timing
!=============================================================================
! Establish mapping between MD an CFD
!-----------------------------------------------------------------------------
subroutine socket_create_map
	implicit none

	!Note coupler cannot be called directly so this socket is needed
	call coupler_create_map

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
	use coupler_internal_md, only : map_md2cfd,map_cfd2md
	use calculated_properties_MD, only : nbins
	use physical_constants_MD, only : np
	use arrays_MD, only :r,v
   	use coupler_module, only : staggered_averages, ncx,ncy,ncz,dx,dz,yg,cfd_code_id
	use messenger, only : MD_COMM
	implicit none

	integer, intent(in) :: iter
	
	integer :: iter_cfd, iter_average, Naverage, save_period, average_period
	integer	:: myid, ierr
	logical, save :: first_time=.true.
	save  average_period, Naverage

	call MPI_COMM_RANK(MD_COMM,myid,ierr)

    if (first_time) then 
	    first_time	= .false.
		call setup_velocity_average
    endif

    iter_average = mod(iter-1, Naverage)+1			! current step
    iter_cfd     = (iter-initialstep)/Naverage +1 	! CFD corresponding step

	!print*, 'iter counts', iter, iter_average, Naverage

    if ( mod(iter_average,average_period) .eq. 0 ) then
		!Collect uc data every save_period cfd iteration but discard the first one which cfd uses for initialisation
	    call cumulative_velocity_average
	endif
    if  (iter_average .eq. Naverage) then
		!Send accumulated results to CFD at the end of average cycle 
	    call send_velocity_average
	endif

contains

!------------------------------------------
	!subroutine CFD_cells_to_MD_compute_cells(cfdis,cfdie,cfdjs,cfdje,cfdks,cfdke, & 
	!										  mdis, mdie, mdjs, mdje, mdks, mdke)
	!	implicit none

	!	integer,intent(in)		:: cfdis,cfdie,cfdjs,cfdje,cfdks,cfdke
	!	integer,intent(out)		:: mdis,mdie,mdjs,mdje,mdks,mdke

	!	ncells(:)

	!end subroutine CFD_cells_to_MD_compute_cells
	
!------------------------------------------
	!THIS SHOULD BE DONE IN THE SETUP!!!!
	subroutine setup_velocity_average
		implicit none

	    average_period = coupler_md_get_average_period() 	! collection interval in the average cycle
	    Naverage = coupler_md_get_md_steps_per_cfd_dt()  	! number of steps in MD average cycle

		! Setup averaging array on first call
		select case(staggered_averages(1))
		case(.true.)
			allocate( mflux(6,ncz,ncx,1))
			mflux = 0
		case(.false.)
			allocate(uvwbin(4,ncz,ncx,1))
			uvwbin = 0.d0
		end select
	end subroutine setup_velocity_average

!------------------------------------------
	subroutine cumulative_velocity_average
		use coupler_internal_md, only : bbox
		use coupler_module, only : icmax,icmin,jcmax,jcmin,kcmax,kcmin,dx,dy,dz,ncx,ncy,ncz,cfd_code_id
		use computational_constants_MD, only : iter, ncells,domain,halfdomain
		use librarymod, only : heaviside, imaxloc
		implicit none

		!Limits of cells to average

		integer							:: n, ixyz
		integer,dimension(3) 			:: cfd_ncells
		integer,dimension(3)			:: ibin,ibin1,ibin2,minbin,maxbin,crossplane,cfdbins
		double precision,dimension(3) 	:: cfd_cellsidelength
		double precision,dimension(3)	:: Vbinsize, ri1, ri2, cnst_top, cnst_bot, rd

		!Velocity measurement for 3D bins throughout the domain
		!Determine bin size
    	cfdbins(1) = ncx
    	cfdbins(2) = ncy
    	cfdbins(3) = ncz

		Vbinsize(1) = dx
		Vbinsize(2) = dy
		Vbinsize(3) = dz

		!print*,'cellsizes',  domain, cfdbins, domain(:) / cfdbins(:),Vbinsize

		rd(:) = (/ 0.d0 , yg(jcmin,1) , 0.d0 /)
		cnst_bot = map_cfd2md(rd)
		rd(:) = (/ 0.d0 , yg(jcmin+1,1) , 0.d0 /)
		cnst_top = map_cfd2md(rd)

		minbin = ceiling((cnst_bot+halfdomain(:))/Vbinsize(:)) + nhb
		maxbin = ceiling((cnst_top+halfdomain(:))/Vbinsize(:)) + nhb

		!print*,'extents', minbin,icmin,jcmin,kcmin,maxbin,icmax,jcmax,kcmax

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
				ibin1(:) = ceiling((ri1+halfdomain(:))/Vbinsize(:)) + nhb
				ibin2(:) = ceiling((ri2+halfdomain(:))/Vbinsize(:)) + nhb
					
				!Exclude molecules outside of domain
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
		case(.false.)
			do n = 1,np
				!Add up current volume mass and momentum densities
				ibin(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
				!Exclude molecules outside of domain
				if (	 ibin(2) .lt. minbin(2) .or. ibin(2) .ge. maxbin(2)) cycle
				!print'(a,12i8,3f10.5)', 'bins & ting',iter,myid,n, minbin, maxbin, ibin, r(:,n)
				uvwbin(1:3,ibin(3),ibin(1),ibin(2)-minbin(2)+1) = uvwbin(1:3,ibin(3),ibin(1),ibin(2)-minbin(2)+1) + v(:,n)
				uvwbin(4,  ibin(3),ibin(1),ibin(2)-minbin(2)+1) = uvwbin(4,  ibin(3),ibin(1),ibin(2)-minbin(2)+1) + 1.d0
				!print*, 'binning part', ibin(3),ibin(1),ibin(2)-minbin(2)+1, uvwbin(4,  ibin(3),ibin(1),ibin(2)-minbin(2)+1), v(:,n)
			enddo
		case default
			call error_abort('Unknown case in staggered_averages')
		end select

	end subroutine cumulative_velocity_average

!------------------------------------------
	subroutine send_velocity_average
		implicit none

		integer		:: icell
        logical 	:: ovr_box_x

	    if (cfd_code_id == couette_parallel) then 
	        ovr_box_x = .true.
	    else
	        ovr_box_x = .false.
	    endif

		!Send data to CFD if send_data flag is set
		select case(staggered_averages(1))	
		! Send velocity flux over surface
		case(.true.)
            call coupler_send(dble(mflux),index_transpose=(/2,3,1/))
			mflux = 0
		! Send velocity in cell centre
		case(.false.)
            call coupler_send(uvwbin,index_transpose=(/2,3,1/))
			uvwbin = 0.d0
		end select

	end subroutine send_velocity_average

end subroutine average_and_send_MD_to_CFD


!===================================================================================
! Run through the particle, check if they are in the overlap region and
! find the CFD box to which the particle belongs		 

subroutine setup_CFD_box(iter,xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD,itm1,itm2)
	implicit none

	!iteration step, it assumes that each MD average starts from iter = 1
	integer, intent(in) 				:: iter 
	integer, intent(inout) 				:: itm1,itm2 
	real(kind=kind(0.d0)),intent(out)	:: xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD

	integer i, j, k, js, je, ib, jb, kb, nib, njb, nkb, ip,m, np_overlap
	integer, save :: ncalls = 0
    logical, save :: firsttime=.true., overlap
	type(cfd_grid_info) cfd_box

    save nib,njb,nkb

	! check if this MD domain overlap with the region in which continuum
	! force is applied
    if (firsttime)then
        firsttime=.false.
        call coupler_md_get(overlap_with_continuum_force=overlap)

        if (overlap) then 
			! get cfd cell sizes
			call coupler_md_get(cfd_box=cfd_box)
			cfd_code_id = coupler_md_get_cfd_id()

			! number of CFD cells in each direction
			nib = cfd_box%icmax - cfd_box%icmin
			njb = 1
			nkb = cfd_box%kcmax - cfd_box%kcmin
		
			! layer extend in local coordinates, i.e. centered on MD box
			xmin = cfd_box%xmin
			xmax = cfd_box%xmax
			ymin = cfd_box%y(cfd_box%jcmax-2)
			ymax = cfd_box%y(cfd_box%jcmax-1)
			zmin = cfd_box%zmin
			zmax = cfd_box%zmax
        
			dx_cfd = cfd_box%dx
			dz_cfd = cfd_box%dz

			allocate(vel_cfd(3,nib,njb,nkb,2),vbuff(1,nkb,nib+1,njb))
			vel_cfd = 0.d0
            allocate(box_average(nib,njb,nkb))
			! vel_fromCFD cell index from which continum velocity is collected
			!jb_constrain =	  njb - 1 ! the second row of cells from the top
		
        endif
    endif
	
	if (iter .eq. 1) then
		! get the previous value of CFD velocities
		! this call must be global at the moment
		! because the whole CFD grid is transferred
		! it will be optimised later
        call coupler_recv(vbuff,index_transpose=(/2,3,1/))
	
        if ( .not. overlap) return

        !swap the time indices and put the cfd velocities work array
        itm1 = mod(itm1,2)+1
        itm2 = mod(itm2,2)+1

        select case (cfd_code_id)
        case (couette_parallel)
			do k=1,nkb
			do j=1,njb
			do i=1,nib
				vel_cfd(1,i,j,k,itm1)= 0.5d0*(vbuff(1,k,i,j)+vbuff(1,k,i+1,j))
			enddo
 			enddo
			enddo
        case (couette_serial)
			do k=1,nkb
			do j=1,njb
			do i=1,nib
				vel_cfd(1,i,j,k,itm1)= vbuff(1,k,i,j)
 			enddo
			enddo
			enddo
        end select

        ! at first CFD step we don't have two values to extrapolate CFD velocities, set inv_dtCFD=0
        if (ncalls .eq. 0) then
			inv_dtCFD = 0.0
        else
			inv_dtCFD = 1.0/coupler_md_get_dt_cfd()
        endif
        ncalls = ncalls + 1

    endif

end subroutine setup_CFD_box

!=============================================================================
!	Apply coupling forces so MD => CFD
!-----------------------------------------------------------------------------
subroutine socket_coupler_apply_continuum_forces(iter)
	use computational_constants_MD, only : delta_t
	use physical_constants_MD, only : np
	use arrays_MD, only :r,v,a
	implicit none
	
	integer, intent(in) :: iter

	integer :: iter_average, Naverage
	real(kind(0.d0)) :: delta_t_CFD
	logical, save :: first_time=.true.
	save Naverage

	if (first_time) then
		first_time = .false.
		Naverage = coupler_md_get_md_steps_per_cfd_dt()
	endif
	iter_average = mod(iter-1, Naverage)+1

	!  use this module version
	call apply_continuum_forces(iter_average)

	
end subroutine socket_coupler_apply_continuum_forces

!=============================================================================
! Apply force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
!-----------------------------------------------------------------------------
subroutine apply_continuum_forces(iter)
    use physical_constants_MD, only : np
	use computational_constants_MD, only : delta_t, nh, ncells, cellsidelength, halfdomain, delta_rneighbr
	use linked_list
	use arrays_MD, only : r, v, a
    use messenger, only : myid
	implicit none

	integer, intent(in) 	:: iter ! iteration step, it assumes that each MD average
								! starts from iter = 1

	integer					:: i, j, k, js, je, ib, jb, kb, nib, njb, nkb, ip,m, np_overlap
	integer, save 			:: ncalls = 0, itm1=1, itm2=2

	integer,allocatable 	:: list(:,:)
	real(kind=kind(0.d0))	:: inv_dtCFD,ymin,ymax,xmin,xmax,zmin,zmax,dx_cfd,dz_cfd

    save xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD

	! here we should have the cell coordinates for the particle ip which is 
	! in the overlap region
	! one has to treat the particle that have left the domain separatley
	! compute the average force for each bin
	call setup_CFD_box(iter,xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD,itm1,itm2)
	call average_over_bin
	call apply_force

contains

!=============================================================================
! Average molecules in overlap region to obtain values for 
! constrained dynamics algorithms
!-----------------------------------------------------------------------------

subroutine average_over_bin
	implicit none

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

	! get the range of j cell index in y direction
	js = min(ncells(2),ceiling((ymin+halfdomain(2))/cellsidelength(2))) + nh
	je = min(ncells(2),ceiling((ymax+halfdomain(2))/cellsidelength(2))) + nh

	!find the maximum number of molecules and allocate a list array	   
	np_overlap = 0 ! number of particles in overlapping reg
    allocate(list(4,np))
	
	do m = 1,np
		if ( r(2,m) >  ymin            .and. r(2,m) < ymax          .and. &
             r(1,m) >= -halfdomain(1)  .and. r(1,m) < halfdomain(1) .and. &
             r(3,m) >= -halfdomain(3)  .and. r(3,m) < halfdomain(3) ) then
			ib = ceiling( (r(1,m) -	xmin) / dx_cfd)
			jb = 1
			kb = ceiling( (r(3,m) - zmin) / dz_cfd)

			!print*, 'COUPLED AVERAGE',myid, xmin, halfdomain(1), zmin, halfdomain(3), &
			! 		m, ib,jb,kb, dx_cfd
           
			np_overlap = np_overlap + 1
			list(1:4, np_overlap) = (/ m, ib, jb, kb /)

			box_average(ib,jb,kb)%np   =  box_average(ib,jb,kb)%np   + 1
			box_average(ib,jb,kb)%v(:) =  box_average(ib,jb,kb)%v(:) + v(:,m)
			box_average(ib,jb,kb)%a(:) =  box_average(ib,jb,kb)%a(:) + a(:,m)
		endif
	enddo

end subroutine average_over_bin

!=============================================================================
! Average molecules using cell lists in overlap region to obtain values for 
! constrained dynamics algorithms
!**************************************
!The cell based version does not work 
!*************************************
!-----------------------------------------------------------------------------

subroutine average_over_bin_cells
	implicit none

	integer		:: n, icell, jcell, kcell
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
	js = min(ncells(2),ceiling((ymin+halfdomain(2))/cellsidelength(2))) + nh
	je = min(ncells(2),ceiling((ymax+halfdomain(2))/cellsidelength(2))) + nh

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

!=============================================================================
! Apply force to molecules in overlap region
!-----------------------------------------------------------------------------

subroutine apply_force
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

		! using the following exptrapolation formula for continuum velocity
		! y = (y2-y1)/(x2-x1) * (x-x2) +y2
        alpha(1) = inv_dtCFD*(vel_cfd(1,ib,1,kb,itm1) - &
                			  vel_cfd(1,ib,1,kb,itm2))

		u_cfd_t_plus_dt(1) = alpha(1) * (iter + 1)*delta_t + vel_cfd(1,ib,1,kb,itm1) 

		acfd =	- box_average(ib,jb,kb)%a(1) / n - inv_dtMD * & 
				( box_average(ib,jb,kb)%v(1) / n - u_cfd_t_plus_dt(1) )
		a(1,ip) = a(1,ip) + acfd

	enddo

end subroutine apply_force

end subroutine apply_continuum_forces

!=============================================================================
! Apply force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
! Adapted serial version written by ES including cells and (originally) fully verified
!-----------------------------------------------------------------------------

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
    save nib,njb,nkb,xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD

	! run through the particle, check if they are in the overlap region
	! find the CFD box to which the particle belongs		  
	! attention to the particle that have left the domain boundaries 
	call setup_CFD_box(iter,xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD,itm1,itm2)

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
 		cbin = jcell - je					!Local no. of overlap cell from 1 to overlap
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
					    -(isumvel-vel_cfd(1,icell,1,kcell,1))/delta_t

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
	double precision, dimension(:,:,:,:,:),allocatable 	:: vel_cfd

	!logical								:: overlap
	type(node), pointer 	        	:: old, current
    save nib,njb,nkb,xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD

	! run through the particle, check if they are in the overlap region
	! find the CFD box to which the particle belongs		  
	! attention to the particle that have left the domain boundaries 
	call setup_CFD_box(iter,xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD,itm1,itm2)

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
				write(99,'(3i8,4f18.9)')iter,jcell,isummol, isumflux, isumforce, isumvel,continuum_Fs(cbin)
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
