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
	implicit none

    ! CFD id
    integer cfd_code_id 

	! type used in applying continuum force
	type cfd_box_sum
		integer np
		real(kind=kind(0.d0))  v(3)
		real(kind=kind(0.d0))  a(3)
	end type cfd_box_sum


	type(cfd_box_sum),dimension(:,:,:),allocatable 		:: box_average

	integer			, dimension(:,:,:,:), allocatable 	:: mflux
	real(kind(0.d0)), dimension(:,:,:,:), allocatable 	:: uvw_md, uvw_cfd
	real(kind(0.d0)), dimension(:,:,:,:,:), allocatable :: stress_cfd

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
	use CPL, only : CPL_create_comm, md_realm
	implicit none

	call CPL_create_comm(md_realm,MD_COMM,ierr)
    prefix_dir = "./md_data/"

end subroutine socket_coupler_invoke

!=============================================================================
! Setup initial times based on coupled calculation
!-----------------------------------------------------------------------------
subroutine socket_coupler_init
    use interfaces
	use CPL, only : coupler_md_init, CPL_create_map, set_coupled_timing, CPL_get,CPL_write_header
	use computational_constants_MD, only : npx,npy,npz,delta_t,elapsedtime, & 
										   Nsteps,initialstep,delta_t, & 
										   globaldomain,initialnunits
	use physical_constants_MD, only : nd,density
	use messenger, only	 :  myid, icoord, icomm_grid
	implicit none

 	integer			 :: naverage, nsteps_cfd, ixyz
	integer          :: match_cellsize
	real(kind(0.d0)) :: delta_t_CFD

	!Establish Domain size from MD inputs
	do ixyz=1,nd
		globaldomain(ixyz) = initialnunits(ixyz) & 	
								/((density/4.d0)**(1.d0/nd))
	enddo

	! If coupled calculation prepare exchange layout
	call coupler_md_init(Nsteps,initialstep,delta_t,icomm_grid,icoord, &
	                     (/ npx,npy,npz /),globaldomain,density)

	!update elapsedtime using new 
	elapsedtime = Nsteps * delta_t

	! Setup the domain and cells in the MD based on coupled data
	call set_params_globdomain_cpld
	call CPL_get(md_cfd_match_cellsize=match_cellsize)
	if (match_cellsize .eq. 0) then
		call set_parameters_cells
 	else
		call set_parameters_cells_coupled
	endif

	!Write coupler information to header file
	call CPL_write_header('./results/coupler_header')

end subroutine socket_coupler_init

!=============================================================================
! get the global domain lenghts from x, y, z array of CFD realm

subroutine set_params_globdomain_cpld
	use computational_constants_MD
	use physical_constants_MD, only: globalnp,volume,nd,np,density
	use CPL, only : CPL_get
	!use coupler_module, only: xL_md, yL_md, zL_md
	implicit none

	integer :: ixyz, n0(3)
	real(kind(0.d0)) :: b0, xL_md, yL_md, zL_md,density_cfd

	!Get domain size and density from coupler
	call CPL_get(xL_md=xL_md, yL_md=yL_md, zL_md=zL_md, density_cfd=density_cfd)

    ! size of cubic FCC cell
    b0 = (4.d0/density)**(1.0d0/3.0d0)
    n0(:) = nint( (/ xL_md, yL_md, zL_md/) / b0)
    initialunitsize(1:3) = b0
    initialnunits(1:3) 	 = n0(:)
   
	!Set MD domain values
	globaldomain(1) = xL_md
	globaldomain(2) = yL_md
	globaldomain(3) = zL_md
	volume = xL_md*yL_md*zL_md

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

    if(irank .eq. iroot) then
        write(*,'(a/a/a,f5.2,a,f5.2,a,f5.2,/,a,3(f5.2),a,/,a,3(I6),/,a)'), &
                "**********************************************************************", &
                "WARNING - this is a coupled run which resets the following parameters:", &
                " density from MD =", density , ' changed to ', density_cfd,              & 
                " hence the cubic FCC side =", b0 ,                  	                  &
                " initialunitsize =", initialunitsize(:)/b0," in b units ",               &     
                " initialnunits   =", initialnunits(:),                                   &
                "**********************************************************************"
    endif
	density = density_cfd

end subroutine set_params_globdomain_cpld

!-----------------------------------------------------------------------------
! Adjust MD cells to match to continuum

subroutine set_parameters_cells_coupled
	use interfaces
	use computational_constants_MD
	use physical_constants_MD
	use polymer_info_MD
	use CPL, only : error_abort, CPL_get
	implicit none

	integer 						:: ixyz,icmax,icmin,jcmax,jcmin,kcmax,kcmin
	integer,dimension(3) 			:: max_ncells,cfd_ncells,cfd_md_cell_ratio
	double precision 				:: rneighbr,dx,dy,dz
	double precision ,dimension(3) 	:: cfd_cellsidelength, maxdelta_rneighbr


	call CPL_get(icmax_olap=icmax,icmin_olap=icmin,jcmax_olap=jcmax,jcmin_olap=jcmin, & 
				 kcmax_olap=kcmax,kcmin_olap=kcmin,dx=dx,dy=dy,dz=dz)

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

    if(irank .eq. iroot) then
        write(*,'(a/a/a,f9.6,/,a,3i8,/,a,3(f8.5),/,a,3(i8),/,a)') &
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


subroutine socket_check_cell_sizes
	use computational_constants_MD, only: cellsidelength,irank,iroot
	use CPL, only: CPL_get
	implicit none

	double precision							:: dx,dy,dz
	double precision,dimension(:),allocatable	:: zg
	double precision,dimension(:,:),allocatable	:: xg,yg

	call CPL_get(dx=dx,dy=dy,dz=dz,xg=xg,yg=yg,zg=zg)
	if (irank .eq. iroot) then
		if( cellsidelength(1) .ge. dx .or. & 
			cellsidelength(2) .ge. dy .or. & 
			cellsidelength(3) .ge. dz 		 ) then
			write(*,*), ""
			write(*,*), "********************************************************************"
			write(*,*), " WARNING ...WARNING ...WARNING ...WARNING ...WARNING ...WARNING ... "
			write(*,*), " MD cell size larger than CFD x,y cell sizes         "
			write(*,*), " cellsidelength = ", cellsidelength
			write(*,'(3(a,f10.5))'),   " dx=",  xg(2,1) - xg(1,1),  & 
									   " dy=",  yg(1,2) - yg(1,1),  & 
									   " dz=",  zg(2  ) - zg(1  )
			write(*,*), "********************************************************************"
			write(*,*), ""
       endif
   endif

end subroutine socket_check_cell_sizes

!=============================================================================
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!
!							SIMULATION
!
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!=============================================================================


!=============================================================================
!  __  __  ____     ___      ___  ____  ____     ____   ___ 
! (  \/  )(  _ \   (__ \    / __)( ___)(  _ \   (  _ \ / __)
!  )    (  )(_) )   / _/   ( (__  )__)  )(_) )   ) _ <( (__ 
! (_/\/\_)(____/   (____)   \___)(__)  (____/   (____/ \___)
!
! Take average of x,y and z components of MD velocity to 
! calculate all components of velocity to pass to contiunuum region 
!-----------------------------------------------------------------------------

subroutine average_and_send_MD_to_CFD(iter)
	use computational_constants_MD, only : initialstep, delta_t, nhb,iblock,jblock,kblock
	use calculated_properties_MD, only : nbins
	use physical_constants_MD, only : np
	use arrays_MD, only :r,v
   	use CPL, only : CPL_get, CPL_realm,coupler_md_get_average_period,CPL_proc_extents
	implicit none

	integer, intent(in) :: iter
	
	integer :: ixyz,icell,jcell,kcell,pcoords(3)
	integer :: iter_cfd, iter_average, save_period 
	integer	:: ierr

	logical, save :: first_time=.true.,staggered_averages(3)
	integer, save :: ncx, ncy, ncz, average_period, jcmin_olap,timestep_ratio, extents(6)
	real(kind(0.d0)),save :: dx, dy, dz
	real(kind(0.d0)),dimension(:),allocatable,save 		:: zg
	real(kind(0.d0)),dimension(:,:),allocatable,save 	:: xg, yg

	!Setup arrays on first call
    if (first_time) then 
	    first_time	= .false.

		!Get processor extents
		pcoords=(/ iblock,jblock,kblock /)
		call CPL_proc_extents(pcoords,CPL_realm(),extents)
		call setup_velocity_average
		call CPL_get(ncx=ncx,ncy=ncy,ncz=ncz,dx=dx,dy=dy,dz=dz,xg=xg,yg=yg,zg=zg, & 
						staggered_averages=staggered_averages,timestep_ratio=timestep_ratio, &
						jcmin_olap=jcmin_olap)
	    average_period = coupler_md_get_average_period() 	! collection interval in the average cycle
    endif
    iter_average = mod(iter-1, timestep_ratio)+1			! current step
    iter_cfd     = (iter-initialstep)/timestep_ratio +1 	! CFD corresponding step

	!Collect uc data every save_period cfd iteration but discard the first one which cfd uses for initialisation
    if ( mod(iter_average,average_period) .eq. 0 ) then
	    call cumulative_velocity_average
	endif

	!Send accumulated results to CFD at the end of average cycle 
    if  (iter_average .eq. timestep_ratio) then
	    call send_velocity_average
	endif

contains

!=============================================================================
! Setup arrays and constants used in average
!-----------------------------------------------------------------------------

	subroutine setup_velocity_average
		use CPL, only : CPL_proc_extents
		use computational_constants_MD, only : iblock,jblock,kblock 
		implicit none

		integer		:: nclx,ncly,nclz

		!Allocate array to size of cells in current processor
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

!=============================================================================
! Cumulativly average velocity 
!-----------------------------------------------------------------------------

	subroutine cumulative_velocity_average
		use CPL, only : CPL_Cart_coords, CPL_realm, VOID, &
						map_md2cfd_global,map_cfd2md_global,globalise,localise
		use computational_constants_MD, only : iter,ncells,domain,halfdomain, & 
												globaldomain,iblock,jblock,kblock,irank
		use librarymod, only : heaviside, imaxloc
		use physical_constants_MD, only : pi
		implicit none

		!Limits of cells to average

		integer							:: i,j,k,n, ixyz
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

		if (any(extents .eq. VOID)) return
		if ((yg(1,extents(3)) .gt. ybcmax) .or. (yg(1,extents(4)+1) .lt. ybcmin)) return

		!Get local extents on processor(s) of interest
		rd(:) = (/ 0.d0 , ybcmin, 0.d0 /)	!Bottom of cell below domain
		avrg_bot = localise(map_cfd2md_global(rd))
		rd2(:) = (/ globaldomain(1) , ybcmax , globaldomain(3)   /)   !Top of cell below domain (= bottom of domain)
		avrg_top = localise(map_cfd2md_global(rd2))

		minbin = ceiling((avrg_bot+halfdomain(:))/dxyz(:)) + nhb
		maxbin = ceiling((avrg_top+halfdomain(:))/dxyz(:)) + nhb

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
				ibin(:) = ceiling((r(:,n)+halfdomain(:))/dxyz(:)) !+ nhb

				!Exclude molecules outside of averaging region
				!if (any(ibin(:).lt.minbin(:)) .or. any(ibin(:).ge.maxbin(:))) cycle
				if (r(2,n).lt.avrg_bot(2) .or. r(2,n).gt.avrg_top(2) ) cycle

				!Exclude molecules which leave processor between rebuilds
				if (ibin(1).lt.1) ibin(1) = 1; if (ibin(1).gt.extents(2)-extents(1)+1) ibin(1) = extents(2)-extents(1)+1
				if (ibin(3).lt.1) ibin(3) = 1; if (ibin(3).gt.extents(6)-extents(5)+1) ibin(3) = extents(6)-extents(5)+1

				!Add velocity and molecular count to bin
				uvw_md(1:3,ibin(1),ibin(2)-minbin(2)+2,ibin(3)) = &
					uvw_md(1:3,ibin(1),ibin(2)-minbin(2)+2,ibin(3)) + v(:,n)
				uvw_md(4,  ibin(1),ibin(2)-minbin(2)+2,ibin(3)) = & 	
					uvw_md(4,  ibin(1),ibin(2)-minbin(2)+2,ibin(3)) + 1.d0

				!DEBUG - each cell is loaded with its physical location in space
				!uvw_md(1:3,ibin(1),ibin(2)-minbin(2)+2,ibin(3)) = globalise(ibin(:)*dxyz(:)-halfdomain(:)-0.5*dxyz(:))
				!uvw_md(1:3,ibin(1),ibin(2)-minbin(2)+2,ibin(3)) = sin(2*pi*(globalise(ibin(:)*dxyz(:)-halfdomain(:)-0.5*dxyz(:)))/globaldomain(:))
				!uvw_md(4,  ibin(1),ibin(2)-minbin(2)+2,ibin(3)) = 1.d0
			enddo

		case default
			call error_abort('Unknown case in staggered_averages')
		end select

	end subroutine cumulative_velocity_average

!=============================================================================
! Send cumulativly average velocity to CFD code
!-----------------------------------------------------------------------------

	subroutine send_velocity_average
		use CPL, only : CPL_send
		use computational_constants_MD, only : iblock,jblock,kblock
		implicit none

		logical :: send_flag,ovr_box_x
		integer :: jcmin_send,jcmax_send, i,k

		!Define arbitary range to send -- TODO move to input file --
		jcmin_send = jcmin_olap; jcmax_send = jcmin_olap

		!Send data to CFD if send_data flag is set
		select case(staggered_averages(1))	
		! Send velocity flux over surface
		case(.true.)
            call CPL_send(dble(mflux),jcmax_send=jcmax_send,jcmin_send=jcmin_send,send_flag=send_flag)
			mflux = 0
		! Send velocity in cell centre
		case(.false.)
            call CPL_send(uvw_md,jcmax_send=jcmax_send,jcmin_send=jcmin_send,send_flag=send_flag)
			uvw_md = 0.d0
		end select

	end subroutine send_velocity_average

end subroutine average_and_send_MD_to_CFD

!==============================================================================
!  ____  _  _  ___  ____  ____  ____    __  __  _____  __    ___
! (_  _)( \( )/ __)( ___)(  _ \(_  _)  (  \/  )(  _  )(  )  / __)
!  _)(_  )  ( \__ \ )__)  )   /  )(     )    (  )(_)(  )(__ \__ \
! (____)(_)\_)(___/(____)(_)\_) (__)   (_/\/\_)(_____)(____)(___/
! 
! Receive mass flux from CFD and create or remove molecules as required
!-----------------------------------------------------------------------------

subroutine insert_remove_molecules
	use physical_constants_MD, only : np
	use computational_constants_MD, only : delta_t, nh, ncells, iter, & 
										   cellsidelength, halfdomain, &
	                                       delta_rneighbr,iblock,jblock,kblock
	use CPL, only : CPL_overlap, CPL_recv, CPL_proc_extents, & 
					CPL_realm, CPL_get, coupler_md_get_dt_cfd
	use linked_list
	use particle_insertion
	implicit none

	logical,save			:: recv_flag, first_time=.true.
	integer 				:: iter_average
	integer					:: icell,jcell,kcell,n,molno, dir
	integer,save			:: cnstd(6),pcoords(3),extents(6),timestep_ratio,nclx,ncly,nclz
	integer,dimension(:,:,:),allocatable				:: mols_change
	integer,dimension(:,:,:),allocatable,save 			:: total_mols_change
	double precision,dimension(3)						:: rin, vin, u
	double precision,dimension(:,:,:),allocatable		:: mass_cfd
	double precision,dimension(:,:,:),allocatable,save	:: mass_change

	! Check processor is inside MD/CFD overlap zone 
	if (.not.(CPL_overlap())) return

	if (first_time) then
		first_time = .false.
		!Save extents of current processor
		pcoords= (/ iblock,jblock,kblock   /)
		call CPL_proc_extents(pcoords,CPL_realm(),extents)

		! Get total number of CFD cells on each processor
		nclx = extents(2)-extents(1)+1
		ncly = extents(4)-extents(3)+1
		nclz = extents(6)-extents(5)+1

		!Allocate array to keep track of running total of mass change
		allocate(mass_change(nclx,ncly,nclz))
		uvw_cfd = 0.d0

		!Get local copies of required simulation parameters
		call CPL_get(icmin_cnst=cnstd(1),icmax_cnst=cnstd(2), & 
	                 jcmin_cnst=cnstd(3),jcmax_cnst=cnstd(4), & 
					 kcmin_cnst=cnstd(5),kcmax_cnst=cnstd(6), &
	                 timestep_ratio=timestep_ratio 				)

	endif

	iter_average = mod(iter-1, timestep_ratio)+1

	! Receive value of CFD mass fluxes at first timestep of timestep_ratio
	if (iter_average .eq. 1) then
		allocate(mass_cfd(nclx,ncly,nclz))
		call CPL_recv(mass_cfd,                                 & 
		              icmin_recv=cnstd(1),icmax_recv=cnstd(2), &
		              jcmin_recv=cnstd(3),jcmax_recv=cnstd(4), &
		              kcmin_recv=cnstd(5),kcmax_recv=cnstd(6), &
		              recv_flag=recv_flag                       )

		!Get total number of molecules to insert before next exchange
		total_mols_change = floor(mass_change)
	endif

	!Get number of molecules per cell to insert/remove
	allocate(mols_change(nclx,ncly,nclz))
	!call get_no_mols_to_append(mass_cfd,mass_change,mols_change)

	!Loop through all cells inserting/removing molecules
	dir = 4; u=(/0.d0,0.d0,0.d0 /)
	jcell = cnstd(4)
	do icell=cnstd(1),cnstd(2)
	do kcell=cnstd(5),cnstd(6)
		do n=1,abs(mass_cfd(icell,jcell,kcell))
			!Check if molecule is to insert or remove
			if (mass_cfd(icell,jcell,kcell) .gt. 0) then
				!Create molecule and insert
				call create_position(icell,jcell,kcell,rin,dir,u)
				call create_velocity(u,vin)
				call insert_molecule(rin,vin)
			else
				!Choose molecule to remove and then remove it
				!call choose_molecule(icell,jcell,kcell,molno)
				call remove_molecule(molno)
			endif
		enddo
	enddo
	enddo

!contains

	!subroutine get_no_mols_to_append(mass_change,mols_change,new_mass_cfd)
	!implicit none

	!Default is to not insert/remove anything
	!mols_change = 0

	!Add any new mass to cumulative mass change
	!if (present(new_mass_cfd)) then
	!	mass_change = mass_change + new_mass_cfd
	!endif

	!Check if multiple molecules are to be inserted each timestep
	!if(any(total_mols_change .gt. timestep_ratio)) then
		
	!endif

	!Loop over all cells
	!jcell = cnstd(4)
	!do icell=cnstd(1),cnstd(2)
	!do kcell=cnstd(5),cnstd(6)

		!Space out insertion of molecules over the time interval
	!	insert_time = mod(iter-1, (timestep_ratio/total_mols_change(icell,jcell,kcell)))+1

	!	if (insert_time .eq. 1) then
	!		mols_change(icell,jcell,kcell) = ceiling(total_mols_change(icell,jcell,kcell)/timestep_ratio)
	!		mass_change(icell,jcell,kcell) = mass_change(icell,jcell,kcell) - mols_change(icell,jcell,kcell)
	!	endif

	!enddo
	!enddo
	
	!end subroutine get_no_mols_to_append

end subroutine insert_remove_molecules


!=============================================================================
!   ___  _____  _  _  ___  ____  ____    __    ____  _  _  ____ 
!  / __)(  _  )( \( )/ __)(_  _)(  _ \  /__\  (_  _)( \( )(_  _)
! ( (__  )(_)(  )  ( \__ \  )(   )   / /(__)\  _)(_  )  (   )(   S
!  \___)(_____)(_)\_)(___/ (__) (_)\_)(__)(__)(____)(_)\_) (__) 
!
! Apply coupling forces so MD => CFD
! Uses value from input flag to choose appropriate routine
!-----------------------------------------------------------------------------

subroutine socket_apply_continuum_forces
	use interfaces, only: error_abort
	use CPL, only: CPL_get
	implicit none

	integer :: constraint_algorithm
	integer :: OT, NCER, Flekkoy, off

	call CPL_get(	constraint_algo	      = constraint_algorithm, & 
					constraint_OT         = OT,        & 
					constraint_NCER       = NCER,      &
					constraint_Flekkoy    = Flekkoy,   &
					constraint_off        = off          )
	
	if ( constraint_algorithm .eq. off ) then
		return
	else if ( constraint_algorithm .eq. OT ) then
		call error_abort("OT constraint force not yet implemented")
	else if ( constraint_algorithm .eq. NCER ) then
		call apply_continuum_forces_NCER
	else if ( constraint_algorithm .eq. Flekkoy ) then
		call apply_continuum_forces_flekkoy
	else
		call error_abort("Unrecognised constraint algorithm flag")
	end if	

end subroutine socket_apply_continuum_forces


!=============================================================================
! Apply coupling forces so MD => CFD
! Force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
!-----------------------------------------------------------------------------

subroutine apply_continuum_forces_NCER
	use physical_constants_MD, only : np
	use computational_constants_MD, only : delta_t, nh, ncells, iter, & 
										   cellsidelength, halfdomain, &
	                                       delta_rneighbr,iblock,jblock,kblock
	use CPL, only : CPL_overlap, CPL_recv, CPL_proc_extents, & 
					CPL_realm, CPL_get, coupler_md_get_dt_cfd
	use linked_list
	implicit none

	integer 				:: iter_average
	integer					:: i,j,k,n,np_overlap
	integer,allocatable 	:: list(:,:)
	real(kind=kind(0.d0))	:: inv_dtCFD,t_fract,CFD_box(6)
	integer,save			:: cnstd(6),pcoords(3),extents(6),timestep_ratio
	logical,save			:: recv_flag, first_time=.true.
	save CFD_box

	! Check processor is inside MD/CFD overlap zone 
	if (.not.(CPL_overlap())) return

	if (first_time) then
		first_time = .false.
		!Save extents of current processor
		pcoords= (/ iblock,jblock,kblock   /)
		call CPL_proc_extents(pcoords,CPL_realm(),extents)

		!Get local copies of required simulation parameters
		call CPL_get(icmin_cnst=cnstd(1),icmax_cnst=cnstd(2), & 
	                 jcmin_cnst=cnstd(3),jcmax_cnst=cnstd(4), & 
					 kcmin_cnst=cnstd(5),kcmax_cnst=cnstd(6), &
	                 timestep_ratio=timestep_ratio              )

		call setup_CFD_box(cnstd,CFD_box,recv_flag)
		!At first CFD step we don't have two values to extrapolate CFD velocities, set inv_dtCFD=0
		inv_dtCFD = 0.0

	else
		inv_dtCFD = 1.0/coupler_md_get_dt_cfd()
	endif
	iter_average = mod(iter-1, timestep_ratio)+1

	! Receive value of CFD velocities at first timestep of timestep_ratio
	if (iter_average .eq. 1) then
		call CPL_recv(uvw_cfd,                                 & 
		              icmin_recv=cnstd(1),icmax_recv=cnstd(2), &
		              jcmin_recv=cnstd(3),jcmax_recv=cnstd(4), &
		              kcmin_recv=cnstd(5),kcmax_recv=cnstd(6), &
		              recv_flag=recv_flag                       )
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
!-----------------------------------------------------------------------------------
subroutine setup_CFD_box(limits,CFD_box,recv_flag)
	use CPL, only : CPL_recv,CPL_proc_portion,localise,map_cfd2md_global,CPL_get,VOID
	implicit none

	!Limits of CFD box to receive data in
	integer,dimension(6) ,intent(in)	:: limits
	!Flag to check if limits cover current processor
	logical				 ,intent(out)	:: recv_flag
	!Returned spacial limits of CFD box to receive data
	real(kind=kind(0.d0)),dimension(6)  :: CFD_box

	integer 	  		  :: portion(6)
	integer		  		  :: nclx,ncly,nclz,ncbax,ncbay,ncbaz,ierr
    logical, save 		  :: firsttime=.true.
	real(kind=kind(0.d0)),dimension(3) 			:: xyzmin,xyzmax
	real(kind(0.d0)),dimension(:),allocatable 	:: zg
	real(kind(0.d0)),dimension(:,:),allocatable :: xg, yg

	! Get total number of CFD cells on each processor
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1

	!Allocate CFD received box
	allocate(uvw_cfd(3,nclx,ncly,nclz))
	uvw_cfd = 0.d0

	!Get limits of constraint region in which to receive data
	call CPL_proc_portion(pcoords,CPL_realm(),limits,portion)

	if (all(portion .ne. VOID)) then
		!Get CFD overlapping grid arrays
		call CPL_get(xg=xg,yg=yg,zg=zg)
		recv_flag = .true.

		! Get physical extents of received region on MD processor
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
	use computational_constants_MD, only : nhb, irank
	use arrays_MD, only : r, v, a
	use CPL, only : CPL_get
	implicit none

	integer				:: ib,jb,kb,n,ixyz
	real(kind(0.d0))	:: dx,dy,dz

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

	!Find the maximum number of molecules and allocate a list array	   
	np_overlap = 0 ! number of particles in overlapping reg
    allocate(list(4,np))
	call CPL_get(dx=dx,dy=dy,dz=dz)

	do n = 1,np

		if ( r(2,n) .gt. CFD_box(3) .and. r(2,n) .lt. CFD_box(4)) then

			ib = ceiling((r(1,n)-CFD_box(1))/dx)
			jb = ceiling((r(2,n)-CFD_box(3))/dy) !CFD_box already in -halfdom -> +halfdom system 
			kb = ceiling((r(3,n)-CFD_box(5))/dz)

			!Add out of domain molecules to nearest cell on domain
			if (ib.lt.1) ib = 1; if (ib.gt.size(box_average,1)) ib = size(box_average,1)
			if (kb.lt.1) kb = 1; if (kb.gt.size(box_average,3)) kb = size(box_average,3)

			np_overlap = np_overlap + 1
			list(1:4, np_overlap) = (/ n, ib, jb, kb /)

			box_average(ib,jb,kb)%np   = box_average(ib,jb,kb)%np   + 1
			box_average(ib,jb,kb)%v(:) = box_average(ib,jb,kb)%v(:) + v(:,n)
			box_average(ib,jb,kb)%a(:) = box_average(ib,jb,kb)%a(:) + a(:,n)

		endif

	enddo

    !Get single average value for slice and store in slice
    !do jb = 1,size(box_average,2)
	!	box_average(:,jb,:)%np  =  sum(box_average(:,jb,:)%np)
	!	do ixyz =1,3
	!		box_average(:,jb,:)%v(ixyz) = sum(box_average(:,jb,:)%v(ixyz))
	!		box_average(:,jb,:)%a(ixyz) = sum(box_average(:,jb,:)%a(ixyz))
	!	enddo
	!enddo

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

	!Loop over all molecules and apply constraint
	do i = 1, np_overlap
		ip = list(1,i)
		ib = list(2,i)
		jb = list(3,i)
		kb = list(4,i)

		n = box_average(ib,jb,kb)%np

		acfd =	- box_average(ib,jb,kb)%a(1) / n - inv_dtMD * & 
				( box_average(ib,jb,kb)%v(1) / n - uvw_cfd(1,ib,jb+cnstd(3)-extents(3),kb) )
		! ib,jb,kb are indicators of which CFD cell in !!!constrained!!! region
		! uvw_cfd is allocated by number of CFD cells on !!!MD!!! processor
		! box_avg is allocated by number of CFD cells in !!!constrained!!! region

		a(1,ip) = a(1,ip) + acfd

	enddo

end subroutine apply_force

end subroutine apply_continuum_forces_NCER

!=============================================================================
! Apply coupling forces so MD => CFD
! Force from Flekkøy (2004) paper to apply continuum stress to molecular region
! inside the overlap region. 
!-----------------------------------------------------------------------------
subroutine apply_continuum_forces_flekkoy
	use physical_constants_MD, only : np
	use computational_constants_MD, only : delta_t, nh, ncells,iter, & 
										   cellsidelength, halfdomain, &
	                                       delta_rneighbr,iblock,jblock,kblock,irank
	use CPL, only : CPL_overlap, CPL_recv, CPL_proc_extents,globalise, & 
					CPL_realm, CPL_get, coupler_md_get_dt_cfd
	use linked_list
	implicit none

	integer 				:: iter_average
	integer					:: i,j,k,n,np_overlap
	integer					:: icmin_olap,icmax_olap,jcmin_olap,jcmax_olap,kcmin_olap,kcmax_olap
	integer,allocatable 	:: list(:,:)
	real(kind=kind(0.d0))	:: inv_dtCFD,t_fract,CFD_box(6)
	real(kind=kind(0.d0)),allocatable,dimension(:,:,:,:)	:: recv_buf
	integer,save			:: cnstd(6),pcoords(3),extents(6),timestep_ratio
	logical,save			:: recv_flag, first_time=.true.
	save CFD_box

	! Check processor is inside MD/CFD overlap zone 
	if (.not.(CPL_overlap())) return

	if (first_time) then
		first_time = .false.
		!Save extents of current processor
		pcoords= (/ iblock,jblock,kblock   /)
		call CPL_proc_extents(pcoords,CPL_realm(),extents)
		!Get local copies of required simulation parameters
		call CPL_get(icmin_cnst=cnstd(1),icmax_cnst=cnstd(2), & 
	                 jcmin_cnst=cnstd(3),jcmax_cnst=cnstd(4), & 
					 kcmin_cnst=cnstd(5),kcmax_cnst=cnstd(6), &
	                 timestep_ratio=timestep_ratio              )
		call setup_CFD_box(cnstd,CFD_box,recv_flag)

	endif

	iter_average = mod(iter-1, timestep_ratio)+1

	! Receive value of CFD velocities at first timestep of timestep_ratio
	if (iter_average .eq. 1 .or. first_time) then
		allocate(recv_buf(9,size(stress_cfd,3),size(stress_cfd,4),size(stress_cfd,5)))
		recv_buf = -666.d0
		call CPL_recv(recv_buf,                                 & 
		              icmin_recv=cnstd(1),icmax_recv=cnstd(2), &
		              jcmin_recv=cnstd(3),jcmax_recv=cnstd(4), &
		              kcmin_recv=cnstd(5),kcmax_recv=cnstd(6), &
		              recv_flag=recv_flag                       )
		stress_cfd = reshape(recv_buf,(/ 3,3,size(recv_buf,2),size(recv_buf,3),size(recv_buf,4) /) )
		deallocate(recv_buf)

		!do i=1,size(stress_cfd,3)
		!do j=1,size(stress_cfd,4)
		!do k=1,size(stress_cfd,5)
		!	write(100+irank,'(a,7i5,4f16.8)'),'recving', iter,iblock,jblock,kblock,i,j,k,stress_cfd(1,1,i,j,k),globalise(i*dx-halfdomain(1)-0.5d0*dx)
		!enddo	
		!enddo	
		!enddo		
	else
		!Linear extrapolation between velocity at t and t+1
	endif

	!Get average over current cell and apply constraint forces
	if (recv_flag .eqv. .true.) then
		call average_over_bin
		call apply_force
	endif

contains

!===================================================================================
! Run through the particle, check if they are in the overlap region and
! find the CFD box to which the particle belongs		 

subroutine setup_CFD_box(limits,CFD_box,recv_flag)
	use CPL, only : CPL_recv,CPL_proc_portion,localise,map_cfd2md_global, & 
					CPL_get,VOID
	implicit none

	!Limits of CFD box to receive data in
	integer,dimension(6) ,intent(in)	:: limits
	!Flag to check if limits cover current processor
	logical				 ,intent(out)	:: recv_flag
	!Returned spacial limits of CFD box to receive data
	real(kind=kind(0.d0)),dimension(6) :: CFD_box

	integer 	  		  :: portion(6)
	integer		  		  :: nclx,ncly,nclz,ncbax,ncbay,ncbaz,ierr
    logical, save 		  :: firsttime=.true.
	real(kind=kind(0.d0)),dimension(3) 			:: xyzmin,xyzmax
	real(kind(0.d0)),dimension(:),allocatable 	:: zg
	real(kind(0.d0)),dimension(:,:),allocatable :: xg, yg

	! Get total number of CFD cells on each processor
	nclx = extents(2)-extents(1)+1
	ncly = extents(4)-extents(3)+1
	nclz = extents(6)-extents(5)+1

	!Allocate CFD received box
	allocate(stress_cfd(3,3,nclx,ncly,nclz))
	stress_cfd = 0.d0

	!Get limits of constraint region in which to receive data
	call CPL_proc_portion(pcoords,CPL_realm(),limits,portion)

	if (all(portion .ne. VOID)) then
		!Get CFD overlapping grid arrays
		call CPL_get(xg=xg,yg=yg,zg=zg)
		recv_flag = .true.

		! Get physical extents of received region on MD processor
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
	use computational_constants_MD, only : nhb, irank
	use arrays_MD, only : r, v, a
	use CPL, only : CPL_get
	implicit none

	integer				:: ib,jb,kb,n
	real(kind(0.d0))	:: dx,dy,dz

	!Zero box averages
	do kb = 1, ubound(box_average,dim=3)
	do jb = 1, ubound(box_average,dim=2)
	do ib = 1, ubound(box_average,dim=1)
		box_average(ib,jb,kb)%np = 0
		box_average(ib,jb,kb)%a	 = 0.d0
	enddo
	enddo
	enddo

	!find the maximum number of molecules and allocate a list array	   
	np_overlap = 0 ! number of particles in overlapping reg
    allocate(list(4,np))
	call CPL_get(dx=dx,dy=dy,dz=dz)

	do n = 1,np

		if ( r(2,n) .gt. CFD_box(3) .and. r(2,n) .lt. CFD_box(4)) then

			ib = ceiling((r(1,n)-CFD_box(1))/dx)
			jb = ceiling((r(2,n)-CFD_box(3))/dy)
			kb = ceiling((r(3,n)-CFD_box(5))/dz)

			!Add out of domain molecules to nearest cell on domain
			if (ib.lt.1) ib = 1; if (ib.gt.size(box_average,1)) ib = size(box_average,1)
			if (kb.lt.1) kb = 1; if (kb.gt.size(box_average,3)) kb = size(box_average,3)

			!Add molecule to overlap list
			np_overlap = np_overlap + 1
			list(1:4, np_overlap) = (/ n, ib, jb, kb /)
			box_average(ib,jb,kb)%np   =  box_average(ib,jb,kb)%np   + 1
			box_average(ib,jb,kb)%a(2) =  box_average(ib,jb,kb)%a(2)  + flekkoy_gweight(r(2,n),CFD_box(3),CFD_box(4))

		endif

	enddo

    !Get single average value for slice and store in slice
    !do jb = 1,size(box_average,2)
	!	box_average(:,jb,:)%np   =  sum(box_average(:,jb,:)%np)
	!	box_average(:,jb,:)%a(2) =  sum(box_average(:,jb,:)%a(2))
    !enddo

end subroutine average_over_bin

!=============================================================================
! Apply force to molecules in overlap region
!-----------------------------------------------------------------------------

subroutine apply_force
	use arrays_MD, only : r,v,a
	use physical_constants_MD, only : density
	use computational_constants_MD, only : irank
	use CPL, only :  rank_world
	implicit none

	integer					:: ib, jb, kb, i, ip, n
	real(kind=kind(0.d0)) 	:: alpha(3), u_cfd_t_plus_dt(3), g, gsum, dx, dy, dz, dA, dV

	real(kind=kind(0.d0)) 	:: 	gsumcheck,gratio, ave_a(3), ave_a_consrnt(3)


	call CPL_get(dx=dx,dy=dy,dz=dz)
	dA = dx*dz
	dV = dx*dy*dz

	!Loop over all molecules and apply constraint
	do i = 1, np_overlap
		ip = list(1,i)
		ib = list(2,i)
		jb = list(3,i)
		kb = list(4,i)

		n = box_average(ib,jb,kb)%np
		g = flekkoy_gweight(r(2,ip),CFD_box(3),CFD_box(4))

		!Gsum is replaced with the fixed value based on density and volume
		gsum = density*dV
		!gsum = box_average(ib,jb,kb)%a(2)

		if (gsum .eq. 0.d0) cycle

		a(:,ip) = a(:,ip) + (g/gsum) * dA * stress_cfd(:,2,ib,jb+cnstd(3)-extents(3),kb) 

        !if (g .ne. 0.d0) then
		!	if (iter .lt. 1000) then
		!		write(1234,'(i3,2i7,3i4,5f12.6)'),rank_world,iter,ip,ib,jb,kb, &
		!				 					  r(2,ip),v(2,ip),a(2,ip),g, & 
		!									 (g/gsum)*dA*stress_cfd(2,2,ib,jb+cnstd(3)-extents(3),kb)
		!	endif
        !endif
	enddo

    !write(99999,'(i2,i7,i7,2f10.2,f6.1,3f9.3,6f12.4)'), rank_world,iter,np_overlap,sum(box_average(:,:,:)%a(2)),  &
    !                gsumcheck, gratio, stress_cfd(:,2,ib,jb+jcmin_recv-extents(3),kb), ave_a,ave_a_consrnt

end subroutine apply_force

! -----------------------------------------------------------
! Function returns Flekkoy weighting for given y and max/min

function flekkoy_gweight(y,ymin,ymax) result (g)
	use CPL, only : error_abort

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
		call error_abort(" flekkoy_gweight error - input y cannot be greater than ymax")
    endif

	!Calculate weighting function
	g = 2*( 1/(L-2*yhat) - 1/L - 2*yhat/(L**2))

end function

end subroutine apply_continuum_forces_flekkoy


!=============================================================================
! 					INQUIRY ROUTINES
!
!=============================================================================

! Get constraint info from CPL module
subroutine socket_get_constraint_info(algorithm,OT,NCER,Flekkoy,off)
	use CPL, only: CPL_get
	implicit none

	integer, intent(out) :: algorithm
	integer, intent(out), optional :: OT,NCER,Flekkoy,off
	
	call CPL_get(	constraint_algo    = algorithm, & 
					constraint_OT      = OT,        & 
					constraint_NCER    = NCER,      &
					constraint_Flekkoy = Flekkoy,   &
					constraint_off     = off          )

end subroutine socket_get_constraint_info


! Get overlap status
function socket_get_overlap_status result(olap)
	use coupler, only: CPL_overlap
	implicit none

	logical :: olap

	olap = CPL_overlap()
	
end function socket_get_overlap_status

! Get domain top minus removed molecules (if appropriate for choice of coupling scheme) 
function socket_get_domain_top() result(top)
	use CPL, only: CPL_get, error_abort
	implicit none

	real(kind(0.d0)) :: yL_md, dy, top, removed_dist
	integer :: algorithm
	integer :: OT,NCER,Flekkoy,off

	call CPL_get(dy=dy,yL_md=yL_md, &
				 constraint_algo    = algorithm, & 
				 constraint_OT      = OT,        & 
				 constraint_NCER    = NCER,      &
				 constraint_Flekkoy = Flekkoy,   &
				 constraint_off     = off          )

	!Specifiy size of removed distance as half a cell
	removed_dist = dy/2.d0

	if ( algorithm .eq. off ) then
		top = yL_md/2.d0
	else if ( algorithm .eq. OT ) then
		top = yL_md/2.d0 - dy/2.d0
	else if ( algorithm .eq. NCER ) then
		top = yL_md/2.d0 - dy/2.d0
	else if ( algorithm .eq. Flekkoy ) then
		top = yL_md/2.d0
	else
		call error_abort("Error in socket_get_domain_top - Unrecognised constraint algorithm flag")
	end if	

end function socket_get_domain_top

! Get domain top minus dy
function socket_get_bottom_of_top_boundary() result(bottom)
	use computational_constants_MD, only: halfdomain
	use coupler_module, only: dy
	implicit none

	real(kind(0.d0)) :: bottom

	bottom = halfdomain(2) - dy

end function socket_get_bottom_of_top_boundary

! Get dy
function socket_get_dy() result(dy_)
	use coupler_module, only: dy
	implicit none

	real(kind(0.d0)) :: dy_

	dy_ = dy

end function socket_get_dy

!======================================================================
!======================================================================
!================================モデル=================================
!======================================================================
!======================================================================
!
!
! TESTING ROUTINES TESTING ROUTINES TESTING ROUTINES TESTING ROUTINES
!
!
!=========✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘✘==========ƦƦƦƦƦƦƦƦ=======
!=============================モデル======================================= ϟƘƦƖןןΣ✘
!======================================================================
!======================================================================

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
	!print'(2a,5i8)', 'MD SIDE',realm_name(realm), rank_world, CPL_overlap,coord

	if (.not.(CPL_overlap())) return

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
	if (.not.(CPL_overlap())) return

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

 	if (.not.(CPL_overlap())) return

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
	if ((CPL_overlap())) call CPL_gather(u,3,gatherlims,gatheru)
	if ((CPL_overlap())) call CPL_scatter(stress,9,scatterlims, &
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

! ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ ♥ 


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
