!=============================================================================
!                   MD coupler Socket
! Routine which interface with the coupler to the CFD code
!
! socket_coupler_init                       Passes MD initialisation variables 
!                                           to coupler_md_init
! average_and_send_MD_to_CFD                Average MD and pass to CFD as BC
!   coupler_md_boundary_cell_average        Calls routines for uc,vc and wc
!     compute_uc_average                    Accumulate averages and send
!     compute_vc_average                    Accumulate averages and send
!     compute_wc_average                    Accumulate averages and send
! socket_coupler_apply_continuum_forces     Apply CFD constraint on MD 
!   apply_continuum_forces                  Calls routines for setup, average 
!                                           & application of constraint
!     call setup_CFD_box                    Setup type storing average for cell
!     call average_over_bin                 Accumulate averages
!     call apply_force                      Apply force
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


    type(cfd_box_sum),dimension(:,:,:),allocatable      :: box_average

    integer         , dimension(:,:,:,:), allocatable   :: mflux
    real(kind(0.d0)), dimension(:,:,:,:), allocatable   :: uvw_md, uvw_cfd,uvw_cfdt_m_dt
    real(kind(0.d0)), dimension(:,:,:,:,:), allocatable :: stress_cfd

contains

!==============================================================================
! Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup
!
!                               SETUP
!
! Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup  Setup
!==============================================================================

!==============================================================================
! Invoke messenger and split inter-communicators 
!------------------------------------------------------------------------------
subroutine socket_coupler_invoke
    use messenger
    use CPL, only : CPL_init, md_realm
    implicit none

    call CPL_init(md_realm, MD_COMM, ierr)
    prefix_dir = "./flowmol/"

end subroutine socket_coupler_invoke

!==============================================================================
! Setup initial times based on coupled calculation
!------------------------------------------------------------------------------
subroutine socket_coupler_init
    !use interfaces
    use CPL, only : CPL_setup_md, CPL_get
    use computational_constants_MD, only : npx, npy, npz, delta_t, &
                                           elapsedtime, Nsteps, initialstep, &
                                           delta_t, globaldomain, initialnunits
    use physical_constants_MD, only : nd,density
    use messenger, only  :  myid, icoord, icomm_grid
    implicit none

    integer          :: naverage, nsteps_cfd, ixyz
    integer          :: match_cellsize
    real(kind(0.d0)) :: delta_t_CFD

    !Establish Domain size from MD inputs
    do ixyz=1,nd
        globaldomain(ixyz) = initialnunits(ixyz) &  
                                /((density/4.d0)**(1.d0/nd))
    enddo

    call CPL_setup_md(icomm_grid, globaldomain, -0.5*globaldomain)


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
    initialnunits(1:3)   = n0(:)
   
    !Set MD domain values
    globaldomain(1) = xL_md
    globaldomain(2) = yL_md
    globaldomain(3) = zL_md
    volume = xL_md*yL_md*zL_md

    ! no need to fix globalnp if we have it already
    if(.not. restart) then
        globalnp=1      !Set number of particles to unity for loop below
        do ixyz=1,nd
            globalnp = globalnp*initialnunits(ixyz)     !One particle per unit cell
        enddo
        globalnp=4*globalnp   !FCC structure in 3D had 4 molecules per unit cell
    endif

    np = globalnp / nproc   
    domain(1) = globaldomain(1) / real(npx, kind(0.d0))         !determine domain size per processor
    domain(2) = globaldomain(2) / real(npy, kind(0.d0))         !determine domain size per processor
    domain(3) = globaldomain(3) / real(npz, kind(0.d0))         !determine domain size per processor

    do ixyz=1,nd
        halfdomain(ixyz) = 0.5d0*domain(ixyz)           !Useful definition
    enddo

    if (config_special_case .eq. 'solid_liquid') then
        liquid_density = density_cfd
        density = density
        if(irank .eq. iroot) then
            write(*,'(a/a/a,f5.2,a,f5.2,a,f5.2,/,a,f5.2,/,a,3(f5.2),a,/,a,3(I6),/,a)') &
                    "**********************************************************************", &
                    "WARNING - this is a coupled run which resets the following parameters:", &
                    " liquid density MD =", liquid_density , ' changed to ', density_cfd,     & 
                    " solid density MD =", density,    & 
                    " hence the cubic FCC side =", b0 ,                                       &
                    " initialunitsize =", initialunitsize(:)/b0," in b units ",               &     
                    " initialnunits   =", initialnunits(:),                                   &
                    "**********************************************************************"
        endif      
    else
        density = density_cfd
        if(irank .eq. iroot) then
            write(*,'(a/a/a,f5.2,a,f5.2,a,f5.2,/,a,3(f5.2),a,/,a,3(I6),/,a)') &
                    "**********************************************************************", &
                    "WARNING - this is a coupled run which resets the following parameters:", &
                    " density from MD =", density , ' changed to ', density_cfd,              & 
                    " hence the cubic FCC side =", b0 ,                                       &
                    " initialunitsize =", initialunitsize(:)/b0," in b units ",               &     
                    " initialnunits   =", initialnunits(:),                                   &
                    "**********************************************************************"
        endif
    endif

end subroutine set_params_globdomain_cpld

!-----------------------------------------------------------------------------
! Adjust MD cells to match to continuum

subroutine set_parameters_cells_coupled
    !use interfaces
    use computational_constants_MD
    use physical_constants_MD
    use polymer_info_MD
    use CPL, only : error_abort, CPL_get
    implicit none

    integer                         :: ixyz,icmax,icmin,jcmax,jcmin,kcmax,kcmin
    integer,dimension(3)            :: max_ncells,cfd_ncells,cfd_md_cell_ratio
    real(kind(0.d0))                :: rneighbr,dx,dy,dz
    real(kind(0.d0)) ,dimension(3)  :: cfd_cellsidelength, maxdelta_rneighbr


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
    cellsidelength(:)     = domain(:)/max_ncells(:)
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
    rneighbr  = rcutoff + (delta_rneighbr(1)+delta_rneighbr(2)+delta_rneighbr(3))/3.d0
    rneighbr2 = rneighbr**2

    !print'(a,6f10.5)', 'domains', domain, x(icmax)-x(icmin),y(jcmax)-y(jcmin),z(kcmax)-z(kcmin)
    !print'(a,12i8)',      'cell', cfd_md_cell_ratio,cfd_ncells,ncells,max_ncells
    !print'(a,6f10.5)', 'cellsize', cfd_cellsidelength(:),cellsidelength(:)

    if (potential_flag .eq. 1) then
        select case(solvent_flag)
        case(0)
            if (rneighbr < R_0) then
                rneighbr = R_0 
                rneighbr2 = R_0**2
                print*, 'Neighbour list distance rneighbr set to &
                        & maximum elongation of polymer spring, ',R_0
            end if
        case(1)
            if (rneighbr < sod_cut) then
                rcutoff   = sod_cut
                rcutoff2  = sod_cut2
                !rneighbr  = rcutoff + delta_rneighbr
                !rneighbr2 = rneighbr**2.d0
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
        write(*,'(a/a/a,3f9.6,/,a,3i8,/,a,3(f8.5),/,a,3(i8),/,a)') &
                    "**********************************************************************", &
                    "WARNING - this is a coupled run which resets the following parameters:", &
                    " Extra cell size for neighbourlist =", delta_rneighbr  ,                 & 
                    " MD computational cells per CFD cell = ",cfd_md_cell_ratio,              &
                    " cellsize =", cellsidelength(:),                                         &     
                    " no of cell  =", ncells(:),                                   &
                    "**********************************************************************"
    endif

    !Ensure CFD cells are not smaller than minimum possible MD cells
    if (any(max_ncells .lt. cfd_ncells)) then
        print*, "max_ncells", max_ncells, "cfd_ncells", cfd_ncells
        call error_abort("ERROR - CFD cellsize smaller than minimum MD computational/averaging cell")
    endif

end subroutine set_parameters_cells_coupled


subroutine socket_check_cell_sizes
    use computational_constants_MD, only: cellsidelength,irank,iroot
    use CPL, only: CPL_get
    implicit none

    real(kind(0.d0))                            :: dx,dy,dz
    real(kind(0.d0)),dimension(:,:,:),allocatable :: xg, yg, zg

    call CPL_get(dx=dx,dy=dy,dz=dz,xg=xg,yg=yg,zg=zg)
    if (irank .eq. iroot) then
        if( cellsidelength(1) .ge. dx .or. & 
            cellsidelength(2) .ge. dy .or. & 
            cellsidelength(3) .ge. dz        ) then
            write(*,*) ""
            write(*,*) "********************************************************************"
            write(*,*) " WARNING ...WARNING ...WARNING ...WARNING ...WARNING ...WARNING ... "
            write(*,*) " MD cell size larger than CFD x,y cell sizes         "
            write(*,*) " cellsidelength = ", cellsidelength
            write(*,'(3(a,f10.5))')   " dx=",  xg(2,1,1) - xg(1,1,1),  & 
                                      " dy=",  yg(1,2,1) - yg(1,1,1),  & 
                                      " dz=",  zg(1,1,2) - zg(1,1,1)
            write(*,*) "********************************************************************"
            write(*,*) ""
       endif
   endif

end subroutine socket_check_cell_sizes

!=============================================================================
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!
!                           SIMULATION
!
! Simulation  Simulation  Simulation  Simulation  Simulation  Simulation  
!=============================================================================


!=============================================================================
!  __  __  ____     ___      ___  ____  ____     ____   ___ 
! (  \/  )(  _ \   (__ \    / __)( ___)(  _ \   (  _ \ / __)
!  )    (  )(_) )   / _/   ( (__  )__)  )(_) )   ) _ <( (__ 
! (_/\/\_)(____/   (____)   \___)(__)  (____/   (____/ \___)
!-----------------------------------------------------------------------------

subroutine average_and_send_MD_to_CFD(iter)
    implicit none

    integer, intent(in) :: iter
    logical :: debug = .false.

    integer :: sendtype=1
    integer, parameter :: velocity=1, stress=2

    if (debug) then
        call debug_check_send(iter)
    else
        select case (sendtype)
        case (velocity)
            call average_and_send_velocity_MD_to_CFD(iter)
        case (stress)
            call send_stress_to_CFD(iter)
        end select
    endif

end subroutine average_and_send_MD_to_CFD

!=============================================================================
! Debug routine package dummy data and send
!-----------------------------------------------------------------------------

subroutine debug_check_send(iter)
    use cpl, only : CPL_get_olap_limits, CPL_my_proc_portion, &
                    CPL_get_no_cells, CPL_send
    implicit none

    integer, intent(in) :: iter

    logical :: send_flag
    integer :: i,j,k,ii,jj,kk,ierr
    integer, dimension(3) :: Ncells
    integer, dimension(6) :: portion, limits
    double precision, dimension(:,:,:,:), & 
         allocatable  :: send_array

    !Get detail for grid
    call CPL_get_olap_limits(limits)
    call CPL_my_proc_portion(limits, portion)
    call CPL_get_no_cells(portion, Ncells)

    ! Pack send_array with cell coordinates. Each cell in the array carries
    ! its global cell number within the overlap region.
    allocate(send_array(4, Ncells(1), Ncells(2), Ncells(3)))
    do i = 1,Ncells(1)
    do j = 1,Ncells(2)
    do k = 1,Ncells(3)
        ! -2 indices to match c++ and python indexing in portion and i,j,k
        ii = i + portion(1) - 2
        jj = j + portion(3) - 2
        kk = k + portion(5) - 2

        send_array(1,i,j,k) = ii
        send_array(2,i,j,k) = jj
        send_array(3,i,j,k) = kk
        send_array(4,i,j,k) = iter
    enddo
    enddo
    enddo

    call CPL_send(send_array, limits, send_flag)

end subroutine debug_check_send

!=============================================================================
! Take average of x,y and z components of MD velocity to 
! calculate all components of velocity to pass to contiunuum region 
!-----------------------------------------------------------------------------

subroutine average_and_send_velocity_MD_to_CFD(iter)
    use computational_constants_MD, only : initialstep, delta_t, nhb,iblock,jblock,kblock
    use calculated_properties_MD, only : nbins
    use physical_constants_MD, only : np
    use arrays_MD, only :r,v
    use CPL, only : CPL_get, CPL_realm, CPL_proc_extents
    implicit none

    integer, intent(in) :: iter
    
    integer :: ixyz,icell,jcell,kcell,pcoords(3)
    integer :: iter_cfd, iter_average, save_period 
    integer :: ierr

    logical, save :: first_time=.true.,staggered_averages(3)
    integer, save :: ncx, ncy, ncz, jcmin_olap,timestep_ratio, extents(6)
    real(kind(0.d0)),save :: dx, dy, dz
    real(kind(0.d0)),dimension(:),allocatable,save      :: zg
    real(kind(0.d0)),dimension(:,:),allocatable,save    :: xg, yg

    !Setup arrays on first call
    if (first_time) then 
        first_time  = .false.

        !Get processor extents
        pcoords=(/ iblock,jblock,kblock /)
        call setup_velocity_average()
        call CPL_get(timestep_ratio=timestep_ratio)

    endif
    iter_average = mod(iter-1, timestep_ratio)+1            ! current step
    iter_cfd     = (iter-initialstep)/timestep_ratio +1     ! CFD corresponding step

    !Collect uc data every save_period cfd iteration but discard the first one which cfd uses for initialisation
    call cumulative_velocity_average()

    !Send accumulated results to CFD at the end of average cycle 
    !if  (iter_average .eq. timestep_ratio) then
        call send_velocity_average
    !endif

contains

!=============================================================================
! Setup arrays and constants used in average
!-----------------------------------------------------------------------------

    subroutine setup_velocity_average()
        use CPL, only : CPL_get_bnry_limits, CPL_my_proc_portion, CPL_get_no_cells
        implicit none

        integer, dimension(3) :: cells
        integer, dimension(6) :: limits, portion

        !Get detail for grid
        call CPL_get_bnry_limits(limits)
        call CPL_my_proc_portion(limits, portion)
        call CPL_get_no_cells(portion, cells)

        ! Setup averaging array on first call
        select case(staggered_averages(1))
        case(.true.)
            allocate( mflux(6,cells(1),cells(2),cells(3)))
            mflux = 0
        case(.false.)
            allocate(uvw_md(4,cells(1),cells(2),cells(3)))
            uvw_md = 0.d0
        end select
    end subroutine setup_velocity_average

!=============================================================================
! Cumulativly average velocity 
!-----------------------------------------------------------------------------

    subroutine cumulative_velocity_average()
        use CPL, only : CPL_map_coord2cell, CPL_map_glob2loc_cell, &
                        CPL_get_bnry_limits, CPL_my_proc_portion, CPL_get
        use computational_constants_MD, only : iter, iblock, jblock, kblock, domain
        use calculated_properties_MD, only : nbins
        use computational_constants_MD, only :  iter,irank
        use messenger, only : globalise                                               
        implicit none

        !Limits of cells to average
        logical :: inbc, inproc
        integer :: n
        integer, dimension(3) :: ibin, bin
        integer, dimension(6) :: limits, portion
        real(kind(0.d0)) :: dy
        real(kind(0.d0)),dimension(3) :: Fbinsize, rglob

        !Get coupler cell size
        call CPL_get(dy=dy)
        
        Fbinsize = domain(:)/nbins(:)
        call CPL_get_bnry_limits(limits)
        call CPL_my_proc_portion(limits, portion)
        !- - - - - - - - - - - - - - - - - - - -
        !Record velocity in cell centre
        !- - - - - - - - - - - - - - - - - - - -
        !Add up current volume mass and momentum densities
        do n = 1,np

            rglob = globalise(r(:,n))
            !We need molecular from half a cell lower to be consistent with OpenFOAM convention here
            rglob(2) = rglob(2) + dy*0.5d0
            inbc = CPL_map_coord2cell(rglob(1), rglob(2), rglob(3), bin)
            inproc = CPL_map_glob2loc_cell(portion, bin, ibin)
            !if (inbc) print'(3i8, 3f10.5,6i8,2l8)', iblock, jblock, kblock, rglob, bin, ibin, inbc, inproc

            if (inproc) then
                !bin(:) = ceiling((r(:,n)+0.5d0*domain(:))/Fbinsize(:))
                !if (bin(1) .ne. ibin(1) .or. bin(3) .ne. ibin(3)) then
                !    print'(a,6i8,6f10.5)', 'CPL_map_coord2cell error:', bin, ibin, Fbinsize, r(:,n)
                !    ibin(1) = bin(1); ibin(3) = bin(3)
                !endif
                !print'(2i8, l5, 7i8,4f10.5)', iter, n, inbc, iblock, jblock, kblock, ibin, limits(4), r(:,n), & 
                !                            uvw_md(1,ibin(1),ibin(2),ibin(3))/uvw_md(4,ibin(1),ibin(2),ibin(3))
                ! Add velocity and molecular count to bin
                uvw_md(1:3,ibin(1),ibin(2),ibin(3)) = uvw_md(1:3,ibin(1),ibin(2),ibin(3)) + v(:,n)
                uvw_md(4,  ibin(1),ibin(2),ibin(3)) = uvw_md(4,  ibin(1),ibin(2),ibin(3)) + 1.d0

            endif

        enddo

        !print*, iter, 'uvw', uvw_md(1,2,:,2)/uvw_md(4,2,:,2)

    end subroutine cumulative_velocity_average

!=============================================================================
! Send cumulativly average velocity to CFD code
!-----------------------------------------------------------------------------

    subroutine send_velocity_average
        use CPL, only : CPL_send, CPL_get_bnry_limits, CPL_overlap, &
                        CPL_my_proc_portion, error_abort    
        use librarymod, only : couette_analytical_fn
        use computational_constants_MD, only :  iter,irank
        use messenger_data_exchange, only : PlaneSum
        implicit none

        logical :: send_flag
        integer :: limits(6), portion(6), cells, send_type
        real(kind(0.d0))  :: v_sum
        real(kind(0.d0)), dimension(:), allocatable  :: u

        if (.not. CPL_overlap()) return

        !Get detail for grid
        call CPL_get_bnry_limits(limits)
        call CPL_my_proc_portion(limits, portion)

        !TEMPTEMPTEMPTEMPTEMPTEMPTEMP
!        allocate(u(26))
!        u = couette_analytical_fn(iter*delta_t,0.65d0,1.d0,44.45d0,26,0)
!        write(100,'(26f16.8)'),u
!        uvw_md(1,:,1,:) = u(11)
!        uvw_md(2,:,1,:) = 0.d0
!        uvw_md(3,:,1,:) = 0.d0
!        uvw_md(4,:,1,:) = 1.d0
        !print*, "WARNING FROM send_velocity_average, setting v component to zero"
    
        send_type = 3

        select case(send_type)
        case(1)
            !Local processor sum (THIS WILL FAIL IF MD/CFD PROCESS BOUNDARIES DON'T LINE UP)
            v_sum = sum(uvw_md(2,:,1,:)/uvw_md(4,:,1,:)) !N.B. velocity usum/msum must be zero
            cells = (portion(2)-portion(1)+1)*(portion(4)-portion(3)+1)*(portion(6)-portion(5)+1)
            uvw_md(2,:,1,:) = uvw_md(2,:,1,:) - uvw_md(4,:,1,:)*v_sum/dble(cells)

        case(2)
            !Global sum (As check is local, this doesn't work)
            v_sum = sum(uvw_md(2,:,1,:)/uvw_md(4,:,1,:))
            cells = (limits(2)-limits(1)+1)*(limits(4)-limits(3)+1)*(limits(6)-limits(5)+1)
            call PlaneSum(v_sum, 2)
            uvw_md(2,:,1,:) = uvw_md(2,:,1,:) - uvw_md(4,:,1,:)*v_sum/dble(cells)

        case(3)
            !Global sum first 
            v_sum = sum(uvw_md(2,:,1,:)/uvw_md(4,:,1,:))
            cells = (limits(2)-limits(1)+1)*(limits(4)-limits(3)+1)*(limits(6)-limits(5)+1)
            call PlaneSum(v_sum, 2)
            uvw_md(2,:,1,:) = uvw_md(2,:,1,:) - uvw_md(4,:,1,:)*v_sum/dble(cells)

            !then local to satisfy OpenFOAM local only checks
            v_sum = sum(uvw_md(2,:,1,:)/uvw_md(4,:,1,:)) !N.B. velocity usum/msum must be zero
            cells = (portion(2)-portion(1)+1)*(portion(4)-portion(3)+1)*(portion(6)-portion(5)+1)
            uvw_md(2,:,1,:) = uvw_md(2,:,1,:) - uvw_md(4,:,1,:)*v_sum/dble(cells)

        case(4)
            !Set v components to zero
            uvw_md(2,:,:,:) = 0.d0
            cells = 1
        end select

        !Check
        !uvw_sum = sum(uvw_md(2,:,1,:))
        !call PlaneSum(uvw_sum, 2)
        !print*, "v_sum", iter,irank, uvw_sum/dble(cells)

        !Send data to CFD if send_data flag is set
        call CPL_send(uvw_md, limits=portion, send_flag=send_flag)
        uvw_md = 0.d0

    end subroutine send_velocity_average

end subroutine average_and_send_velocity_MD_to_CFD

!=============================================================================
! Take calculated and send to CFD
!-----------------------------------------------------------------------------

subroutine send_stress_to_CFD(iter)
        use CPL, only : CPL_send, CPL_get_bnry_limits, &
                        CPL_my_proc_portion, error_abort    
        use messenger, only : localise_bin
        use calculated_properties_MD, only : nbins, Pxybin, volume_momentum
        implicit none

        integer, intent(in) :: iter

        logical :: send_flag
        integer :: limits(6), portion(6), nd
        double precision, allocatable, dimension(:,:,:,:) :: send_buf

        !Get detail for grid
        call CPL_get_bnry_limits(limits)
        call CPL_my_proc_portion(limits, portion)

        !Copy three components of stress
        nd = 3
        allocate(send_buf(nd, nbins(1), nbins(2), nbins(3)))
        send_buf(1,:,:,:) = Pxybin(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,1,2)
        send_buf(2,:,:,:) = Pxybin(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,1,3)
        send_buf(3,:,:,:) = Pxybin(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,2,3)

        call CPL_send(send_buf, limits=portion, send_flag=send_flag)

end subroutine send_stress_to_CFD

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
                    CPL_realm, CPL_get, CPL_cfd_dt
    use particle_insertion, only : insert_molecule, create_position, &
                                    create_velocity, remove_molecule
    implicit none

    logical,save            :: recv_flag, first_time=.true.
    integer                 :: iter_average, mass
    integer                 :: icell,jcell,kcell,n,molno, dir
    integer,save            :: cnstd(6),pcoords(3),extents(6),timestep_ratio,nclx,ncly,nclz
    integer,dimension(:,:,:),allocatable                :: mols_change
    integer,dimension(:,:,:),allocatable,save           :: total_mols_change
    real(kind(0.d0)),dimension(3)                       :: rin, vin, u
    real(kind(0.d0)),dimension(:,:,:,:),allocatable       :: mass_cfd
    real(kind(0.d0)),dimension(:,:,:),allocatable,save  :: mass_change

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
                     timestep_ratio=timestep_ratio              )

    endif

    iter_average = mod(iter-1, timestep_ratio)+1

    ! Receive value of CFD mass fluxes at first timestep of timestep_ratio
    if (iter_average .eq. 1) then
        allocate(mass_cfd(1, nclx,ncly,nclz))
        call CPL_recv(mass_cfd, limits = cnstd, recv_flag=recv_flag )

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
        mass = mass_cfd(1,icell,jcell,kcell)
        do n=1,abs(mass)
            !Check if molecule is to insert or remove
            if (mass .gt. 0) then
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

subroutine socket_apply_continuum_forces()
    use CPL, only: CPL_get, error_abort
    implicit none

    integer :: constraint_algorithm
    integer :: OT, NCER, Flekkoy, CV, off, debug=-666

    call CPL_get(   constraint_algo       = constraint_algorithm, & 
                    constraint_OT         = OT,        & 
                    constraint_NCER       = NCER,      &
                    constraint_Flekkoy    = Flekkoy,   &
                    constraint_CV         = CV,        &
                    constraint_off        = off          )
    
    if ( constraint_algorithm .eq. off ) then
        return
    else if ( constraint_algorithm .eq. OT ) then
        call error_abort("OT constraint force not yet implemented")
    else if ( constraint_algorithm .eq. NCER ) then
        call apply_continuum_forces_NCER
    else if ( constraint_algorithm .eq. Flekkoy ) then
        call apply_continuum_forces_flekkoy
    else if ( constraint_algorithm .eq. CV ) then
        call apply_CV_force()
        !call DEBUG_apply_continuum_forces_CV()
    else if ( constraint_algorithm .eq. debug ) then
        call debug_check_recv()
    else
        call error_abort("Unrecognised constraint algorithm flag")
    end if  

end subroutine socket_apply_continuum_forces


subroutine debug_check_recv()
    use cpl, only : CPL_get_olap_limits, CPL_my_proc_portion, &
                    CPL_get_no_cells, CPL_recv, CPL_overlap
    use computational_constants_MD, only :  iter,irank
    implicit none

    logical :: recv_flag,no_error
    integer :: i,j,k,ii,jj,kk,ierr
    integer, dimension(3) :: Ncells
    integer, dimension(6) :: portion, limits
    double precision, dimension(:,:,:,:), allocatable  :: recv_array

    !Get detail for grid
    call CPL_get_olap_limits(limits)
    call CPL_my_proc_portion(limits, portion)
    call CPL_get_no_cells(portion, Ncells)

    !Coupled Recieve and print
    allocate(recv_array(4, Ncells(1), Ncells(2), Ncells(3)))
    recv_array = 0.d0
    call CPL_recv(recv_array, limits, recv_flag)

    ! Check that every processor inside the overlap region receives the cell correctly
    ! number.  
    if (CPL_overlap()) then
        no_error = .true.
        do i = 1, Ncells(1)
        do j = 1, Ncells(2)
        do k = 1, Ncells(3)
            ! -2 indices to match c++ and python indexing in portion and i,j,k
            ii = i + portion(1) - 2
            jj = j + portion(3) - 2
            kk = k + portion(5) - 2

            if ((dble(ii) - recv_array(1,i,j,k)) .gt. 1e-8) then 
                print'(a,2i5,a,i5,a,i6,a,f10.5)', "ERROR -- portion in x: ", portion(1:2), & 
                       " MD rank: ", irank, " cell i: ",ii, & 
                       " recv_array: ", recv_array(1,i,j,k)
                no_error = .false.
            endif
            if ((dble(jj) - recv_array(2,i,j,k)) .gt. 1e-8) then 
                print'(a,2i5,a,i5,a,i6,a,f10.5)', "ERROR -- portion in y: ", portion(3:4), & 
                       " MD rank: ", irank, " cell j: ", jj , & 
                       " recv_array: ", recv_array(2,i,j,k)
                no_error = .false.  
            endif
            if ((dble(kk) - recv_array(3,i,j,k)) .gt. 1e-8) then 
                print'(a,2i5,a,i5,a,i6,a,f10.5)', "ERROR -- portion in z: ", portion(5:6), & 
                       " MD rank: ", irank, " cell k: ", kk , & 
                       " recv_array: ", recv_array(3,i,j,k)
                no_error = .false.
            endif
            if ((dble(iter)-1 - recv_array(4,i,j,k)) .gt. 1e-8) then 
                print'(a,f10.1,a,i10)', "ERROR -- recieved iter time: ", recv_array(4,i,j,k), &
                                      " does not match current iter: ",iter 
                no_error = .false.
            endif
        enddo
        enddo
        enddo
        if (no_error) then
            print*, iter-1, "Recv test passed"
        endif

    endif

end subroutine debug_check_recv

!=============================================================================
! Apply coupling forces so MD => CFD
! Force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
!-----------------------------------------------------------------------------
subroutine apply_continuum_forces_NCER()
    use physical_constants_MD, only : np
    use computational_constants_MD, only : delta_t, nh, ncells, iter, & 
                                           cellsidelength, halfdomain, &
                                           delta_rneighbr,iblock,jblock,kblock
    use CPL, only : CPL_overlap, CPL_recv, CPL_proc_extents, & 
                    CPL_realm, CPL_get, CPL_cfd_dt
    use linked_list
    implicit none

    integer                 :: iter_average
    integer                 :: i,j,k,n,np_overlap
    integer,allocatable     :: list(:,:)
    real(kind(0.d0))        :: inv_dtCFD,t_fract,CFD_box(6)
    integer,save            :: cnstd(6),pcoords(3),extents(6),timestep_ratio
    logical,save            :: recv_flag, first_time=.true.
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
        inv_dtCFD = 1.0/CPL_cfd_dt()
    endif
    iter_average = mod(iter-1, timestep_ratio)+1

    ! Receive value of CFD velocities at first timestep of timestep_ratio
    if (iter_average .eq. 1 .or. first_time) then
        uvw_cfdt_m_dt = uvw_cfd !Store previous timestep velocity
        call CPL_recv(uvw_cfd, limits=cnstd, recv_flag=recv_flag )

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
!-----------------------------------------------------------------------------------
subroutine setup_CFD_box(limits,CFD_box,recv_flag)
    use CPL, only : CPL_proc_portion,CPL_cfd2md,CPL_get,VOID
    use messenger, only : localise
    implicit none

    !Limits of CFD box to receive data in
    integer,dimension(6) ,intent(in)    :: limits
    !Flag to check if limits cover current processor
    logical              ,intent(out)   :: recv_flag
    !Returned spacial limits of CFD box to receive data
    real(kind(0.d0)),dimension(6)  :: CFD_box

    integer               :: portion(6)
    integer               :: nclx,ncly,nclz,ncbax,ncbay,ncbaz,ierr
    logical, save         :: firsttime=.true.
    real(kind(0.d0)),dimension(3)           :: xyzmin,xyzmax
    real(kind(0.d0)),dimension(:,:,:),allocatable :: xg, yg, zg

    ! Get total number of CFD cells on each processor
    nclx = extents(2)-extents(1)+1
    ncly = extents(4)-extents(3)+1
    nclz = extents(6)-extents(5)+1

    !Allocate CFD received box
    allocate(uvw_cfd(3,nclx,ncly,nclz))
    uvw_cfd = 0.d0
    allocate(uvw_cfdt_m_dt(3,nclx,ncly,nclz))
    uvw_cfdt_m_dt = 0.d0

    !Get limits of constraint region in which to receive data
    call CPL_proc_portion(pcoords,CPL_realm(),limits,portion)

    if (all(portion .ne. VOID)) then
        !Get CFD overlapping grid arrays
        call CPL_get(xg=xg, yg=yg, zg=zg)
        recv_flag = .true.

        ! Get physical extents of received region on MD processor
        xyzmin(1) = xg(portion(1),1,1); xyzmax(1) = xg(portion(2)+1,1,1)
        xyzmin(2) = yg(1,portion(3),1); xyzmax(2) = yg(1,portion(4)+1,1)
        xyzmin(3) = zg(1,1,portion(5)); xyzmax(3) = zg(1,1,portion(6)+1)

        !Map to local MD processor
        xyzmax = localise(CPL_cfd2md(xyzmin))
        xyzmax = localise(CPL_cfd2md(xyzmax))

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

subroutine average_over_bin()
    use computational_constants_MD, only : nhb, irank
    use arrays_MD, only : r, v, a
    use CPL, only : CPL_get, cpl_md_bc_slice
    implicit none

    integer             :: ib,jb,kb,n,ixyz
    real(kind(0.d0))    :: dx,dy,dz

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

            !Get cell -- N.B. CFD_box already in -halfdom -> +halfdom system 
            ib = ceiling((r(1,n)-CFD_box(1))/dx)
            jb = ceiling((r(2,n)-CFD_box(3))/dy) 
            kb = ceiling((r(3,n)-CFD_box(5))/dz)

            !Add out of domain molecules to nearest cell on domain
            if (ib.lt.1) ib = 1; if (ib.gt.size(box_average,1)) ib = size(box_average,1)
            if (kb.lt.1) kb = 1; if (kb.gt.size(box_average,3)) kb = size(box_average,3)

            !Add molecule to overlap list
            np_overlap = np_overlap + 1
            list(1:4, np_overlap) = (/ n, ib, jb, kb /)

            box_average(ib,jb,kb)%np   = box_average(ib,jb,kb)%np   + 1
            box_average(ib,jb,kb)%v(:) = box_average(ib,jb,kb)%v(:) + v(:,n)
            box_average(ib,jb,kb)%a(:) = box_average(ib,jb,kb)%a(:) + a(:,n)

        endif

    enddo

    if (cpl_md_bc_slice .eq. 1) then
        !Get single average value for slice and store in slice
        do jb = 1,size(box_average,2)
            box_average(:,jb,:)%np  =  sum(box_average(:,jb,:)%np)
            do ixyz =1,3
                box_average(:,jb,:)%v(ixyz) = sum(box_average(:,jb,:)%v(ixyz))
                uvw_cfd(ixyz,:,jb+cnstd(3)-extents(3),:) = &
                    sum(uvw_cfd(ixyz,:,jb+cnstd(3)-extents(3),:))&
                    /real((size(uvw_cfd,2)*size(uvw_cfd,4)))
                box_average(:,jb,:)%a(ixyz) = sum(box_average(:,jb,:)%a(ixyz))
            enddo
        enddo
    end if


end subroutine average_over_bin

!=============================================================================
! Apply force to molecules in overlap region
!-----------------------------------------------------------------------------

subroutine apply_force
    use arrays_MD, only : r, a
    use computational_constants_MD, only : irank, vflux_outflag, CV_conserve, tplot, initialstep
    use CPL, only : error_abort
    use module_record_external_forces, only : record_external_forces
    implicit none

    integer :: NCER_type
    integer ib, jb, kb, i, j, k, molno, n
    real(kind(0.d0)) alpha(3), u_cfd_t_plus_dt(3), inv_dtMD, acfd(3)

    ! set the continnum constraints for the particle in the bin
    ! speed extrapolation add all up
    inv_dtMD =1.d0/delta_t

    !Loop over all molecules and apply constraint
    do i = 1, np_overlap
        molno = list(1,i)
        ib = list(2,i)
        jb = list(3,i)
        kb = list(4,i)

        n = box_average(ib,jb,kb)%np

        ! ib,jb,kb are indicators of which CFD cell in !!!constrained!!! region
        ! uvw_cfd is allocated by number of CFD cells on !!!MD!!! processor
        ! box_avg is allocated by number of CFD cells in !!!constrained!!! region
        NCER_type = 2
        select case(NCER_type)
        case(0)
            ! Difference in velocity in continuum (with no average force term from constrained dynamics)
            acfd(:) = - box_average(ib,jb,kb)%a(:) / n  &
                      - inv_dtMD * (   uvw_cfdt_m_dt(:,ib,jb+cnstd(3)-extents(3),kb) & 
                                     - uvw_cfd      (:,ib,jb+cnstd(3)-extents(3),kb) )
        case(1)
            ! NCER with no force term but including correct special "discretisation" to apply propertional
            ! constraint to equations of motion (This is same as control scheme use by 
            ! Borg, M. K. and Macpherson, G. and Reese, J. (2010) Mol. Sim., 36 (10). pp. 745-757.) 
            acfd(:) =   - inv_dtMD * ( box_average(ib,jb,kb)%v(:) / n - uvw_cfd(:,ib,jb+cnstd(3)-extents(3),kb) )
        case(2)
            ! Full NCER including correct special "discretisation" to apply proportional
            ! constraint to equations of motion
            !print('(a, 3f10.3, a, i10)'), 'a: ', box_average(ib,jb,kb)%a(:), 'n: ', n 
            acfd(:) =   - box_average(ib,jb,kb)%a(:) / real(n,kind(0.d0)) & 
                        - inv_dtMD * (    box_average(ib,jb,kb)%v(:) / real(n,kind(0.d0)) & 
                                        - uvw_cfd(:,ib,jb+cnstd(3)-extents(3),kb) )
        case default 
            call error_abort("Incorrect case in apply_continuum_forces_NCER")
        end select

        if (vflux_outflag .eq. 4) then
            if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
                call record_external_forces( acfd(:) , r(:,molno))
            endif
        endif

        a(:,molno) = a(:,molno) + acfd(:)

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
    use CPL, only : CPL_overlap, CPL_recv, CPL_proc_extents, & 
                    CPL_realm, CPL_get, CPL_cfd_dt, CPL_scatter, &
                    comm_style, comm_style_gath_scat, comm_style_send_recv, &
                    error_abort, CPL_get_olap_limits, CPL_my_proc_portion, CPL_get_no_cells
    use linked_list
    implicit none

    integer :: iter_average
    integer :: i, j, k, n, np_overlap
    integer :: icmin_olap, icmax_olap
    integer :: jcmin_olap, jcmax_olap
    integer :: kcmin_olap, kcmax_olap
    integer, allocatable :: list(:,:)
    real(kind(0.d0)) :: inv_dtCFD, t_fract
    real(kind(0.d0)) :: emptybuf(0, 0, 0, 0)
    real(kind(0.d0)), allocatable, dimension(:,:,:,:) :: recv_buf

    integer, save :: timestep_ratio
    integer, dimension(3), save :: pcoords, recv_cells
    integer, dimension(6), save :: cnstd, limits, portion, extents

    logical, save :: recv_flag, first_time=.true.
    real(kind(0.d0)), save :: CFD_box(6)

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

        !Get detail for grid
        call CPL_get_olap_limits(limits)
        call CPL_my_proc_portion(limits, portion)
        call CPL_get_no_cells(portion, recv_cells)
        allocate(recv_buf(9, recv_cells(1), recv_cells(2), recv_cells(3)))
        recv_buf = -666.d0
        call CPL_recv(recv_buf, limits=limits, recv_flag=recv_flag)

        stress_cfd = reshape(recv_buf, &
            (/3, 3, size(recv_buf,2), size(recv_buf,3), size(recv_buf,4)/))

        deallocate(recv_buf)

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
!-----------------------------------------------------------------------------------

subroutine setup_CFD_box(limits,CFD_box,recv_flag)
    use CPL, only : CPL_proc_portion,CPL_cfd2md, & 
                    CPL_get,VOID
    use messenger, only : localise
    implicit none

    !Limits of CFD box to receive data in
    integer,dimension(6) ,intent(in)    :: limits
    !Flag to check if limits cover current processor
    logical              ,intent(out)   :: recv_flag
    !Returned spacial limits of CFD box to receive data
    real(kind(0.d0)),dimension(6)       :: CFD_box

    integer               :: portion(6)
    integer               :: nclx,ncly,nclz,ncbax,ncbay,ncbaz,ierr
    logical, save         :: firsttime=.true.
    real(kind(0.d0)),dimension(3)               :: xyzmin,xyzmax
    real(kind(0.d0)),dimension(:,:,:),allocatable :: xg, yg, zg

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
        xyzmin(1) = xg(portion(1),1,1); xyzmax(1) = xg(portion(2)+1,1,1)
        xyzmin(2) = yg(1,portion(3),1); xyzmax(2) = yg(1,portion(4)+1,1)
        xyzmin(3) = zg(1,1,portion(5)); xyzmax(3) = zg(1,1,portion(6)+1)

        !Map to local MD processor
        xyzmin = localise(CPL_cfd2md(xyzmin))
        xyzmax = localise(CPL_cfd2md(xyzmax))

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

    integer             :: ib,jb,kb,n
    real(kind(0.d0))    :: dx,dy,dz

    !Zero box averages
    do kb = 1, ubound(box_average,dim=3)
    do jb = 1, ubound(box_average,dim=2)
    do ib = 1, ubound(box_average,dim=1)
        box_average(ib,jb,kb)%np = 0
        box_average(ib,jb,kb)%a  = 0.d0
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
            box_average(ib,jb,kb)%a(2) =  box_average(ib,jb,kb)%a(2) + &
                    flekkoy_gweight(r(2,n),CFD_box(3),CFD_box(4))

        endif

    enddo

end subroutine average_over_bin

!=============================================================================
! Apply force to molecules in overlap region
!-----------------------------------------------------------------------------

subroutine apply_force
    use arrays_MD, only : r,v,a
    use physical_constants_MD, only : density
    use computational_constants_MD, only : irank
    use calculated_properties_MD, only: pressure
    use CPL, only :  rank_world
    implicit none

    integer                 :: ib, jb, kb, i, molno, n
    real(kind=kind(0.d0))   :: alpha(3), u_cfd_t_plus_dt(3), g, gsum, dx, dy, dz, dA, dV

    real(kind=kind(0.d0))   ::  gsumcheck,gratio, ave_a(3), ave_a_consrnt(3)


    call CPL_get(dx=dx,dy=dy,dz=dz)
    dA = dx*dz
    dV = dx*dy*dz

    !Loop over all molecules and apply constraint
    do i = 1, np_overlap
        molno = list(1,i)
        ib = list(2,i)
        jb = list(3,i)
        kb = list(4,i)

        n = box_average(ib,jb,kb)%np
        g = flekkoy_gweight(r(2,molno),CFD_box(3),CFD_box(4))

        !Sum of all gsum for all molecules
        gsum = box_average(ib,jb,kb)%a(2)

        if (gsum .eq. 0.d0) cycle
        if (n .eq. 0) cycle

        a(:,molno) = a(:,molno) + (g/gsum) * dA * stress_cfd(:,2,ib,jb+cnstd(3)-extents(3),kb) 
        a(2,molno) = a(2,molno) - (g/gsum) * dA * pressure

    enddo

end subroutine apply_force

! -----------------------------------------------------------
! Function returns Flekkoy weighting for given y and max/min

function flekkoy_gweight(y,ymin,ymax) result (g)
    use CPL, only : error_abort

    real(kind=kind(0.d0)), intent(in)   :: y, ymin, ymax
    real(kind=kind(0.d0))               :: g, L, yhat

    !Define local coordinate as const runs from 0 < y < L/2
    L = ymax - ymin
    yhat = y - ymin - 0.5*L

    !Sanity Check and exceptions
    if (yhat .lt. 0.d0) then
        g = 0.d0
        return
    elseif (yhat .gt. 0.5*L) then
        call error_abort(" flekkoy_gweight error - input y cannot be greater than ymax")
    endif

    !Calculate weighting function
    g = 2.d0*( 1.d0/(L-2.d0*yhat) - 1.d0/L - 2.d0*yhat/(L**2.d0))

end function

end subroutine apply_continuum_forces_flekkoy



subroutine DEBUG_apply_continuum_forces_CV()
    use computational_constants_MD, only :iblock,jblock,kblock, iter
    use calculated_properties_MD, only : nbins
    use CPL, only : CPL_proc_extents,CPL_realm, & 
                    CPL_get
    implicit none

    integer :: i,j,k
    integer :: nclx,ncly,nclz
    integer :: pcoords(3),extents(6),cnstd(6)
    real(kind(0.d0)),allocatable,dimension(:,:,:,:)    :: u_CFD 
    real(kind(0.d0)),allocatable,dimension(:,:,:,:,:)  :: CFD_stress,CFD_flux

    !Save extents of current processor
    pcoords = (/ iblock,jblock,kblock /)
    call CPL_proc_extents(pcoords,CPL_realm(),extents)

    ! Get total number of CFD cells on each processor
    nclx = extents(2)-extents(1)+1
    ncly = extents(4)-extents(3)+1
    nclz = extents(6)-extents(5)+1

    call CPL_get(icmin_cnst=cnstd(1),icmax_cnst=cnstd(2), & 
                 jcmin_cnst=cnstd(3),jcmax_cnst=cnstd(4), & 
                 kcmin_cnst=cnstd(5),kcmax_cnst=cnstd(6)    )

    !Allocate CFD received boxes
    allocate(u_CFD(nclx,ncly,nclz,3)); u_CFD = 0.d0
    allocate(CFD_stress(nclx,ncly,nclz,3,6)); CFD_stress = 0.d0
    allocate(CFD_flux(nclx,ncly,nclz,3,6)); CFD_flux =0.d0

    call socket_get_velocity_debug(u_CFD)
    call socket_get_fluxes_and_stresses_debug(CFD_stress,CFD_flux)

    do i=1,nclx
    do j=1,ncly
    do k=1,nclz
        if (any(abs(u_CFD(i,j,k,:)) .gt. 0.0000001)) then
            print'(a,4i6,3f27.10)','recv vel  ',iter, i,j,k,u_CFD(i,j,k,:)
            !print*,'recv vel  ',iter, i,j,k,u_CFD(i,j,k,:)
        endif
        if (any(abs(CFD_stress(i,j,k,:,:)) .gt. 0.0000001)) then
            print'(a,4i6,2f27.10)','recv stress',iter, i,j,k,CFD_stress(i,j,k,1,2),CFD_stress(i,j,k,1,4)
            !print*,'recv stress',iter, i,j,k,CFD_stress(i,j,k,1,2),CFD_stress(i,j,k,1,4)
        endif
    enddo
    enddo
    enddo

contains

    subroutine socket_get_velocity_debug(u_CFD)
        use CPL, only : CPL_recv, CPL_scatter, comm_style, &
                        comm_style_send_recv, comm_style_gath_scat
        implicit none

        real(kind(0.d0)),intent(inout), & 
            allocatable,dimension(:,:,:,:)      :: u_CFD

        logical :: recv_flag
        integer :: npercell
        real(kind(0.d0)) :: emptybuf(0, 0, 0, 0)
        real(kind(0.d0)), & 
            allocatable,dimension(:,:,:,:)      :: recv_buf

        npercell = 3
        allocate(recv_buf(npercell,nclx,ncly,nclz))
        recv_buf = -666.d0
        print*, "socket_get_velocity_debug"
        call CPL_recv(recv_buf, limits=cnstd, recv_flag=recv_flag)
        do i=cnstd(1),cnstd(2)
        do j=cnstd(3),cnstd(4)
        do k=cnstd(5),cnstd(6)
            u_CFD(i,j,k,:) = recv_buf(:,i,j,k)
            !print'(a,4i3,3e27.10)','recv vel  ',iter, i,j,k,recv_buf(:,i,j,k)
        enddo
        enddo
        enddo
        deallocate(recv_buf)

    end subroutine socket_get_velocity_debug


    subroutine socket_get_fluxes_and_stresses_debug(CFD_stress,CFD_flux)
        use CPL, only : CPL_recv, CPL_scatter, comm_style, &
                        comm_style_send_recv, comm_style_gath_scat
        implicit none

        real(kind(0.d0)),intent(inout), & 
            allocatable,dimension(:,:,:,:,:)    :: CFD_stress,CFD_flux

        logical :: recv_flag
        integer :: npercell
        real(kind(0.d0)) :: emptybuf(0, 0, 0, 0)
        real(kind(0.d0)), & 
            allocatable,dimension(:,:,:,:)      :: recv_buf

        npercell = 18
        allocate(recv_buf(npercell,nclx,ncly,nclz))

        ! vvvvvvv FLUX EXCHANGE COMING SOON vvvvvvv
        !Get Fluxes 
        CFD_flux = 0.d0
        recv_buf = -666.d0
        ! ^^^^^^^ FLUX EXCHANGE COMING SOON ^^^^^^^

        !Get Stresses
        CFD_stress = 0.d0
        recv_buf = -666.d0
        print*, "socket_get_fluxes_and_stresses_debug"
        call CPL_recv(recv_buf, limits=cnstd, recv_flag=recv_flag )
        do i=cnstd(1),cnstd(2)
        do j=cnstd(3),cnstd(4)
        do k=cnstd(5),cnstd(6)
            CFD_stress(i,j,k,:,:) = reshape(recv_buf(:,i,j,k), (/ 3,6 /))
            !print*,'recv stress',iter,i,j,k,CFD_stress(i,j,k,1,2),CFD_stress(i,j,k,1,4)
        enddo
        enddo
        enddo
        deallocate(recv_buf)

    end subroutine socket_get_fluxes_and_stresses_debug

end subroutine DEBUG_apply_continuum_forces_CV



subroutine socket_get_velocity(u_CFD, lbl)
    use CPL, only : CPL_recv, CPL_scatter, comm_style, &
                    comm_style_send_recv, comm_style_gath_scat, &
                    CPL_proc_portion, CPL_get_cnst_limits, CPL_my_proc_portion, CPL_get_no_cells
    use computational_constants_MD, only :iblock,jblock,kblock,iter,nhb
    use calculated_properties_MD, only : nbins
    use CPL, only : CPL_proc_extents,CPL_realm, & 
                    CPL_get,error_abort
    implicit none

	integer,dimension(6),intent(in)		    :: lbl(6)
    real(kind(0.d0)),intent(inout),allocatable,dimension(:,:,:,:)    :: u_CFD 

    logical :: recv_flag
    integer :: i,j,k,ii,jj,kk,ncbax,ncbay,ncbaz
    ! Extra plus one here for y as CPL_library cannot simulate case
    ! where MD goes beyond top of domain
    integer :: y_MDcells_above_CFD = 1
    integer :: Ncells(3),cnstd(6),portion(6)
    integer :: npercell,nclx,ncly,nclz
    real(kind(0.d0)) :: emptybuf(0, 0, 0, 0)
    real(kind(0.d0)), & 
        allocatable,dimension(:,:,:,:)      :: recv_buf

    ! Get limits of constrained region
    call CPL_get_cnst_limits(cnstd)
    call CPL_my_proc_portion(cnstd, portion)
    call CPL_get_no_cells(portion, Ncells)

    !Check CFD receiving array
    if (size(u_CFD,1)-2*nhb(1) .ne. Ncells(1)) then
        print*, 'nbins(1) = ', nbins(1), 'nclx = ', Ncells(1)
        call error_abort("u_CFD -- nbins(1) not equal to x CFD cells")
    endif
    !if (size(u_CFD,2)-2*nhb(2) .ne. Ncells(2))then
    !    print*, 'nbins(2) = ', nbins(2), 'ncly = ', Ncells(2)
    !     call error_abort("u_CFD -- nbins(2) not equal to y CFD cells")
    !endif
    if (size(u_CFD,3)-2*nhb(3) .ne. Ncells(3)) then
        print*, 'nbins(3) = ', nbins(3), 'nclz = ', Ncells(3)
        call error_abort("u_CFD -- nbins(3) not equal to z CFD cells")
    endif
    if (size(u_CFD,4) .ne. 3) call error_abort("Fourth index of velocity should be size three")
    u_CFD = 0.d0

    npercell = 3
    allocate(recv_buf(npercell,Ncells(1),Ncells(2),Ncells(3)))
    recv_buf = -666.d0
    call CPL_recv(recv_buf, limits=portion, recv_flag=recv_flag )

    if (lbl(2)-lbl(1)+1 .ne. Ncells(1)) then
        print'(2(a,2i8))', "F CV constraint cells in x = ",lbl(1)-nhb(1), lbl(2)-nhb(1), & 
                           " portion of domain = ", portion(1), portion(2)
        call error_abort("socket_get_velocity Error -- constraint region not applied to recieved region of domain")
    endif
    if (lbl(4)-lbl(3)+1 .ne. Ncells(2)) then
        print'(2(a,2i8))', "F CV constraint cells in y = ",lbl(3)-nhb(2), lbl(4)-nhb(2), & 
                           " portion of domain = ", portion(3), portion(4)
        call error_abort("socket_get_velocity Error -- constraint region not applied to recieved region of domain")
    endif
    if (lbl(6)-lbl(5)+1 .ne. Ncells(3)) then
        print'(2(a,2i8))', "F CV constraint cells in z = ",lbl(5)-nhb(3), lbl(6)-nhb(3), & 
                            " portion of domain = ", portion(5), portion(6)
        call error_abort("socket_get_velocity Error -- constraint region not applied to recieved region of domain")
    endif


    !Map cells to location in array
    do i = 1, Ncells(1)
    do j = 1, Ncells(2)
    do k = 1, Ncells(3)
        ! -2 indices to match c++ and python indexing in portion and i,j,k
        ! NOTE, constraint region information is not used at all here!!!
        ! this is deliberate so MD code specifies where to apply constraint
        ii = lbl(1) + i - 1; jj = lbl(3) + j - 1; kk = lbl(5) + k - 1
!        ii = i + portion(1) - 2*nhb(1) + 1
!        jj = j + portion(3) - 2*nhb(2) + 1
!        kk = k + portion(5) - 2*nhb(3) + 1

        u_CFD(ii,jj,kk,:) = recv_buf(:,i,j,k)
        !print'(6i6,3f10.5)', i,j,k,ii,jj,kk,u_CFD(ii,jj,kk,:)

    enddo
    enddo
    enddo
    deallocate(recv_buf)


!    !print*, cnstd(3),cnstd(4),ncbax,ncbaz, lbl
!    do i=1,ncbax
!    do j=cnstd(3),cnstd(4)
!    do k=1,ncbaz

!        !Map F_CV_limits to correct cells in recv array
!        ii = lbl(1) + i - 1; jj = lbl(3) + j - cnstd(3); kk = lbl(5) + k - 1
!        u_CFD(ii,jj,kk,:) = recv_buf(:,i,j,k)

!    enddo
!    enddo
!    enddo
!    deallocate(recv_buf)

end subroutine socket_get_velocity


subroutine socket_get_fluxes_and_stresses(CFD_stress,CFD_flux)
    use CPL, only : CPL_recv, CPL_scatter, comm_style, &
                    comm_style_send_recv, comm_style_gath_scat, &
                    CPL_proc_portion, CPL_proc_extents,CPL_realm, & 
                    CPL_get,error_abort
    use computational_constants_MD, only :iblock,jblock,kblock,iter,nhb
    use calculated_properties_MD, only : nbins
    implicit none

    real(kind(0.d0)),intent(inout), & 
        allocatable,dimension(:,:,:,:,:)    :: CFD_stress,CFD_flux

    logical :: recv_flag
    integer :: npercell,nclx,ncly,nclz
    integer :: i,j,k,ii,jj,kk,ncbax,ncbay,ncbaz
    ! Extra plus one here for y as CPL_library cannot simulate case
    ! where MD goes beyond top of domain
    integer :: y_MDcells_above_CFD = 1
    integer :: pcoords(3),extents(6),cnstd(6), portion(6)
    real(kind(0.d0)) :: emptybuf(0, 0, 0, 0)
    real(kind(0.d0)), & 
        allocatable,dimension(:,:,:,:)      :: recv_buf

    !Save extents of current processor
    pcoords = (/ iblock,jblock,kblock /)
    call CPL_proc_extents(pcoords,CPL_realm(),extents)

    ! Get total number of CFD cells on each processor
    nclx = extents(2)-extents(1)+1
    ncly = extents(4)-extents(3)+1
    nclz = extents(6)-extents(5)+1

    call CPL_get(icmin_cnst=cnstd(1),icmax_cnst=cnstd(2), & 
                 jcmin_cnst=cnstd(3),jcmax_cnst=cnstd(4), & 
                 kcmin_cnst=cnstd(5),kcmax_cnst=cnstd(6)    )

    !Allocate CFD received boxes
    if (size(CFD_stress,1)-2*nhb(1) .ne. nclx) then
        print*, 'nbins(1) = ', nbins(1), 'nclx = ', nclx
        call error_abort("CFD_stress -- nbins(1) not equal to x CFD cells")
    endif
    !if (size(CFD_stress,2)-2*nhb(2) .ne. ncly)then
    !    print*, 'nbins(2) = ', nbins(2), 'ncly = ', ncly
    !     call error_abort("CFD_stress -- nbins(2) not equal to y CFD cells")
    !endif
    if (size(CFD_stress,3)-2*nhb(3) .ne. nclz) then
        print*, 'nbins(3) = ', nbins(3), 'nclz = ', nclz
        call error_abort("CFD_stress -- nbins(3) not equal to z CFD cells")
    endif
    if (size(CFD_stress,4) .ne. 3) call error_abort("Fourth index of stress should be size three")
    if (size(CFD_stress,5) .ne. 6) call error_abort("Fifth index of stress should be size six")
    if (any(shape(CFD_stress) .ne. shape(CFD_flux))) call error_abort("CFD_Flux should match CFD_stress")

    npercell = 18
    allocate(recv_buf(npercell,nclx,ncly,nclz))

    ! vvvvvvv COMING SOON vvvvvvv
    !Get Fluxes 
    CFD_flux = 0.d0
    recv_buf = -666.d0
    ! ^^^^^^^ COMING SOON ^^^^^^^

    !Get Stresses
    CFD_stress = 0.d0
    recv_buf = -666.d0
    print*, "socket_get_fluxes_and_stresses"
    call CPL_recv(recv_buf, limits=cnstd, recv_flag=recv_flag )
    call CPL_proc_portion(pcoords,CPL_realm(),cnstd,portion)
    ncbax = portion(2)-portion(1)+1
    ncbay = portion(4)-portion(3)+1
    ncbaz = portion(6)-portion(5)+1
    do i=1,ncbax
    do j=cnstd(3),cnstd(4)
    do k=1,ncbaz
        ! Extra plus one here for y as CPL_library cannot simulate case
        ! where MD goes beyond top of domain
        ii = i + nhb(1); jj = j + nhb(2)+y_MDcells_above_CFD; kk = k + nhb(3)
        CFD_stress(ii,jj,kk,:,:) = reshape(recv_buf(:,i,j,k), (/ 3,6 /))
        !if (abs(CFD_stress(ii,jj,kk,1,2)) .gt. 0.00001) then
        !    print'(a,7i5,2f27.10)','recv stress',iter,i,j,k,ii,jj,kk,CFD_stress(ii,jj,kk,1,2),CFD_stress(ii,jj,kk,1,4)
        !endif
    enddo
    enddo
    enddo
    deallocate(recv_buf)

end subroutine socket_get_fluxes_and_stresses


!=============================================================================
!                   INQUIRY ROUTINES
!
!=============================================================================

! Get constraint info from CPL module
subroutine socket_get_constraint_info(algorithm,OT,NCER,CV,Flekkoy,off)
    use CPL, only: CPL_get
    implicit none

    integer, intent(out) :: algorithm
    integer, intent(out), optional :: OT,NCER,Flekkoy,CV,off
    
    call CPL_get(   constraint_algo    = algorithm, & 
                    constraint_OT      = OT,        & 
                    constraint_NCER    = NCER,      &
                    constraint_Flekkoy = Flekkoy,   &
                    constraint_CV      = CV,        &
                    constraint_off     = off          )

end subroutine socket_get_constraint_info


! Get overlap status
function socket_get_overlap_status() result(olap)
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
    integer :: OT,NCER,Flekkoy,CV,off,debug=-666

    call CPL_get(dy=dy,yL_md=yL_md, &
                 constraint_algo    = algorithm, & 
                 constraint_OT      = OT,        & 
                 constraint_NCER    = NCER,      &
                 constraint_Flekkoy = Flekkoy,   &
                 constraint_CV      = CV,   &
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
    else if ( algorithm .eq. CV ) then
        top = yL_md/2.d0
    else if ( algorithm .eq. debug ) then
        top = yL_md/2.d0
    else
        call error_abort("Error in socket_get_domain_top - Unrecognised constraint algorithm flag")
    end if  

end function socket_get_domain_top


! Get domain top minus removed molecules (if appropriate for choice of coupling scheme) 
function socket_get_domain_bottom() result(bottom)
    use CPL, only: CPL_get, error_abort
    implicit none

    real(kind(0.d0)) :: yL_md, dy, bottom, removed_dist
    integer :: algorithm
    integer :: OT,NCER,Flekkoy,CV,off,debug=-666

    call CPL_get(dy=dy,yL_md=yL_md, &
                 constraint_algo    = algorithm, & 
                 constraint_OT      = OT,        & 
                 constraint_NCER    = NCER,      &
                 constraint_Flekkoy = Flekkoy,   &
                 constraint_CV      = CV,   &
                 constraint_off     = off          )

    !Specifiy size of removed distance as half a cell
    removed_dist = dy/2.d0

    if ( algorithm .eq. off ) then
        bottom = -yL_md/2.d0
    else if ( algorithm .eq. OT ) then
        bottom = -yL_md/2.d0
    else if ( algorithm .eq. NCER ) then
        bottom = -yL_md/2.d0
    else if ( algorithm .eq. Flekkoy ) then
        bottom = -yL_md/2.d0
    else if ( algorithm .eq. CV ) then
        bottom = -yL_md/2.d0
    else if ( algorithm .eq. debug ) then
        bottom = -yL_md/2.d0
    else
        call error_abort("Error in socket_get_domain_bottom - Unrecognised constraint algorithm flag")
    end if  

end function socket_get_domain_bottom

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

#endif

end module md_coupler_socket
