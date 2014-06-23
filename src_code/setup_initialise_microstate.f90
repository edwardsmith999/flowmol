!------------------------------------------------------------------------------
!                           Initialise Microstate
!        Set up position and velocity of each molecule (the microstate)
!
!------------------------------------------------------------------------------

module module_initialise_microstate

    use physical_constants_MD
    use computational_constants_MD
    use arrays_MD
    use calculated_properties_MD

    double precision            :: angle, rand  !Define variables
    double precision            :: v13          !Mag of v1 and v3 vectors
    
end module module_initialise_microstate
!------------------------------------------------------------------------------

subroutine setup_initialise_microstate
    use interfaces
    use module_initialise_microstate
	use module_read_input, only : COUETTE_t,COUETTE_Re,COUETTE_Uwall, & 
								  COUETTE_H,COUETTE_slidewall,COUETTE_ixyz
    implicit none

    integer     ::  n

    !Choose initial molecular positions using configurational flag
    select case(initial_config_flag)
    case(0)
        call setup_initialise_lattice          !Setup FCC lattice
        call setup_location_tags               !Setup locn of fixed mols
    case(1)
        select case (config_special_case)
        case('sparse_fene')
            call setup_initialise_sparse_FENE
            call setup_location_tags           !Setup locn of fixed mols
        case('dense_fene')
            call setup_initialise_lattice      !Numbering for FENE bonds
            call setup_lattice_dense_FENE_info !Chain IDs, etc
            call setup_location_tags           !Setup locn of fixed mols
        case('solid_liquid')
            call setup_initialise_solid_liquid
            call setup_location_tags               !Setup locn of fixed mols
        case('polymer_brush')
            call setup_initialise_polymer_brush
            call setup_location_tags
        case('rubber_liquid')
            call setup_initialise_solid_liquid
            call setup_lattice_dense_FENE_info !Chain IDs, etc
        case('concentric_cylinders')
            call setup_initialise_concentric_cylinders
            tag = free 
        case('fill_cylinders')
            call parallel_io_import_cylinders
            call setup_initialise_fill_cylinders
            call setup_cylinder_tags
        case('fill_cylinders_fene_solution')
            call parallel_io_import_cylinders
            call setup_initialise_fill_cylinders
            call setup_cylinder_tags!_equilibrate
            call setup_cylinder_FENE_solution
        case('rotate_cylinders')
            call setup_cylinder_tags
        case default
            call error_abort('Unidentified configuration special case')
        end select
    case(2)
        call error_abort('Initial configuration file input not yet developed')
    case default
        call error_abort('Unidentified initial configuration flag') 
    end select

    if (rtrue_flag.eq.1) then
        do n=1,np    !Initialise global true positions
            rtrue(1,n) = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
            rtrue(2,n) = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
            rtrue(3,n) = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
        end do
    endif

    if (ensemble .eq. tag_move) then
        rtether = r                                   !Init tether pos
        do n = 1,np
            call read_tag(n)                          !Read tag, assign props
        enddo
    end if

    !Choose initial molecular velocities using velocity flag
    select case(initial_velocity_flag)
    case(0)
        call setup_initialise_velocities                 !Setup initial velocities
    case(1)
        select case (trim(velocity_special_case))
        case('debug')
            call setup_initialise_velocities_test
        case('taylor_green')
            call setup_initialise_velocities_TG_parallel
        case('couette_analytical')
            call setup_initialise_velocities                 !Setup initial velocities
            call set_velocity_field_from_couette_analytical(COUETTE_t,COUETTE_Re, & 
															COUETTE_Uwall,COUETTE_H, &
															COUETTE_slidewall,COUETTE_ixyz)
        case('dns')
            call setup_initialise_velocities                 !Setup initial velocities
            call set_velocity_field_from_DNS_restart(trim(DNS_filename),DNS_ngx,DNS_ngy,DNS_ngz)
        case default
            call error_abort('Unidentified initial velocities_special_case')    
        end select
    case default
        call error_abort('Unidentified initial velocity flag')  
    end select

end subroutine setup_initialise_microstate

!==================================================================================
!----------------------------------------------------------------------------------
!Initialise Positions
!Set up the intial position of the particles in an FCC lattice
subroutine setup_initialise_lattice
    use module_initialise_microstate
    use messenger
    use messenger_data_exchange, only : globalSum
#if USE_COUPLER
    use coupler
    use md_coupler_socket, only: socket_get_domain_top
#endif
    implicit none

    integer :: j, n, nl, nx, ny, nz
    integer, dimension(nd) :: p_units_lb, p_units_ub 
    double precision :: domain_top
    double precision, dimension (nd):: rc, c !Temporary variable

    p_units_lb(1) = (iblock-1)*floor(initialnunits(1)/real((npx),kind(0.d0)))
    p_units_ub(1) =  iblock *ceiling(initialnunits(1)/real((npx),kind(0.d0)))
    p_units_lb(2) = (jblock-1)*floor(initialnunits(2)/real((npy),kind(0.d0)))
    p_units_ub(2) =  jblock *ceiling(initialnunits(2)/real((npy),kind(0.d0)))
    p_units_lb(3) = (kblock-1)*floor(initialnunits(3)/real((npz),kind(0.d0)))
    p_units_ub(3) =  kblock *ceiling(initialnunits(3)/real((npz),kind(0.d0)))

    !Set top of domain initially
    domain_top = globaldomain(2)/2.d0

#if USE_COUPLER

    if (jblock .eq. npy) then
        domain_top = socket_get_domain_top()
    endif

#endif

    !Molecules per unit FCC structure (3D)
    n  = 0      !Initialise global np counter n
    nl = 0      !Initialise local np counter nl

    !Inner loop in y (useful for setting connectivity)
    do nz=p_units_lb(3),p_units_ub(3)
    c(3) = (nz - 0.75d0)*initialunitsize(3) !- halfdomain(3) 
    do nx=p_units_lb(1),p_units_ub(1)
    c(1) = (nx - 0.75d0)*initialunitsize(1) !- halfdomain(1)
    do ny=p_units_lb(2),p_units_ub(2)
    c(2) = (ny - 0.75d0)*initialunitsize(2) !- halfdomain(2) 

        do j=1,4    !4 Molecules per cell

            rc(:) = c(:)
            select case(j)
            case(2)
                rc(1) = c(1) + 0.5d0*initialunitsize(1)
                rc(3) = c(3) + 0.5d0*initialunitsize(3)
            case(3)
                rc(2) = c(2) + 0.5d0*initialunitsize(2)
                rc(3) = c(3) + 0.5d0*initialunitsize(3)
            case(4)
                rc(1) = c(1) + 0.5d0*initialunitsize(1)
                rc(2) = c(2) + 0.5d0*initialunitsize(2)
            case default
            end select

            n = n + 1   !Move to next particle
            
            !Remove molecules from top of domain if constraint applied
            if (jblock .eq. npy) then
                ! Note rc is in global coordinates from 0 to globaldomain while Domaintop
                ! is given in global coordinates from -halfglobaldomain to halfglobaldomain
                !print*, "MOLECULAR REMOVAL TURNED OFF IN setup_initialise_lattice"
                if (rc(2)-globaldomain(2)/2.d0 .gt. domain_top) cycle 
            endif

            !Check if molecule is in domain of processor
            if(rc(1).lt. domain(1)*(iblock-1)) cycle
            if(rc(1).ge. domain(1)*(iblock  )) cycle
            if(rc(2).lt. domain(2)*(jblock-1)) cycle
            if(rc(2).ge. domain(2)*(jblock  )) cycle
            if(rc(3).lt. domain(3)*(kblock-1)) cycle
            if(rc(3).ge. domain(3)*(kblock  )) cycle

            !If molecules is in the domain then add to total
            nl = nl + 1 !Local molecule count

            !Correct to local coordinates
            r(1,nl) = rc(1)-domain(1)*(iblock-1)-halfdomain(1)
            r(2,nl) = rc(2)-domain(2)*(jblock-1)-halfdomain(2)
            r(3,nl) = rc(3)-domain(3)*(kblock-1)-halfdomain(3)

        enddo

    enddo
    enddo
    enddo


    !Correct local number of particles on processor
    np = nl

    !Establish global number of particles on current process
    globalnp = np
    call globalSum(globalnp)

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp

#if USE_COUPLER

    if (jblock .eq. npy .and. iblock .eq. 1 .and. kblock .eq. 1) then
        print*, '*********************************************************************'
        print*, '*WARNING - TOP LAYER OF DOMAIN REMOVED IN LINE WITH CONSTRAINT FORCE*'
        print*, 'Removed from', domain_top, 'to Domain top', globaldomain(2)/2.d0
        print*, 'Number of molecules reduced from',  & 
                 4*initialnunits(1)*initialnunits(2)*initialnunits(3), 'to', globalnp
        print*, '*********************************************************************'
    endif

#endif

end subroutine setup_initialise_lattice
!--------------------------------------------------------------------------------
!FENE info
!-------------------------------------------------------------------------------
!Assign chainIDs, subchainIDs, global molecule numbers, etc...
subroutine setup_lattice_dense_FENE_info
    use interfaces
    use polymer_info_MD
    use messenger
    use messenger_data_exchange, only : globalSum
    use physical_constants_MD, only: np 
    implicit none

    integer :: n
    integer :: chainID
    integer :: subchainID
    integer :: modcheck
    integer :: solvent_selector
    integer, dimension(nproc) :: proc_chains, proc_nps
    character(256) :: string
    
    proc_chains(:)         = 0
    proc_nps(:)            = 0

    !intbits = bit_size(monomer(1)%bin_bflag(1))
    
    modcheck = 0 + mod(np,nmonomers) + mod(4*initialnunits(2)/npy,nmonomers)
    if (modcheck.ne.0) then
        string = 'Number of molecules in the y direction on each processor      &
                  must be exactly divisible by the polymer chain length. Please & 
                  change the chain length in the input file. A chain length of  &
                  4 should (hopefully) always work.'
        call error_abort(string)
    end if

    do n=1,np+extralloc
        monomer(n)%chainID      = 0
        monomer(n)%subchainID   = 0
        monomer(n)%funcy        = 0
        monomer(n)%glob_no      = 0
        monomer(n)%bin_bflag(:) = 0
    end do  

    chainID    = 1
    subchainID = 0
    do n=1,np

        select case (solvent_flag)
        case (0)
            solvent_selector = 0
        case (1,2)
            solvent_selector = mod((n-1)/nmonomers,solvent_ratio)
        case default
        end select

        !if (foam_tag(n).eq.foam) then
        !   solvent_selector = 0
        !else
        !   solvent_selector = 1
        !end if
        
        select case (solvent_selector)
        case(0) !POLYMER
            
            subchainID = subchainID + 1
            if (subchainID .gt. nmonomers) then
                subchainID = 1
                chainID    = chainID + 1
            end if
    
            monomer(n)%chainID    = chainID
            monomer(n)%subchainID = subchainID
            monomer(n)%glob_no    = n   
            
            
            if (subchainID.eq.1) then
                call connect_to_monomer(subchainID+1,n)
            else if (subchainID.eq.nmonomers) then
                call connect_to_monomer(subchainID-1,n)             
            else
                call connect_to_monomer(subchainID+1,n)
                call connect_to_monomer(subchainID-1,n)
            end if

            
        case(1:) !SOLVENT

            !SET MOLECULES TO BE ATHERMAL SOLVENT MOLECULES
            monomer(n)%chainID     = 0
            monomer(n)%subchainID  = 1
            monomer(n)%glob_no     = n
            monomer(n)%funcy       = 0
            bond(:,n)              = 0
        
        case default
        end select
            
        !print*, '-n,c,sc,f,g------------------------------------------------'
        !print*, n, monomer(n)%chainID,monomer(n)%subchainID,monomer(n)%funcy,monomer(n)%glob_no
        !print*, 'bin_bflag -', monomer(n)%bin_bflag
        !print*, '==========================================================='

    end do

    proc_chains(irank) = chainID
    proc_nps(irank)    = np
    call globalSum(proc_chains,nproc)
    call globalSum(proc_nps,nproc)
    
    do n=1,np
        if (monomer(n)%chainID.ne.0) then
            monomer(n)%chainID     = monomer(n)%chainID + sum(proc_chains(1:irank)) - proc_chains(irank)
        end if
        monomer(n)%glob_no     = monomer(n)%glob_no + sum(proc_nps(1:irank))    - proc_nps(irank)
    end do

    nchains = sum(proc_chains)

end subroutine setup_lattice_dense_FENE_info

!--------------------------------------------------------------------------------
!FENE info
!-------------------------------------------------------------------------------
!Assign chainIDs, subchainIDs, global molecule numbers, etc...
subroutine setup_cylinder_FENE_solution
    use interfaces
    use polymer_info_MD
    use messenger
    use messenger_data_exchange, only : globalSum
    use arrays_MD, only: r
    use physical_constants_MD, only: np, density
    use concentric_cylinders, only: r_oi, r_io
    use librarymod, only: cpolariser
    implicit none

    logical :: connectable 
    integer :: n
    integer :: maxchainID
    integer :: nchainsremove, fluid_np, cyl_np
    integer, dimension(nproc) :: proc_chains, proc_nps
    real(kind(0.d0)) :: check, concentration

    call initialise_info
    
    ! Connect all chains, removed later
    call connect_all_possible_chains(maxchainID)

    proc_chains(irank) = maxchainID
    proc_nps(irank) = np

    if (irank.eq.1) then
        print*, 'Target concentration: ', targetconc
    end if

    ! Remove chains to get target concentration (as close as possible) 
    concentration = real(nmonomers*proc_chains(irank))/real(fluid_np)
    nchainsremove = nint((concentration - targetconc)*real(fluid_np) &
                    /real(nmonomers))

    ! If nchains is too many, remove, otherwise we can't do anything now
    if (nchainsremove .gt. 0) then
        call remove_chains(nchainsremove)
        proc_chains(irank) = proc_chains(irank) - nchainsremove 
    else
        print*, 'Concentration before chain removal:', concentration
        print*, 'Density, nchains connectable: ', density, proc_chains(irank) 
        print*, 'The target concentration is not obtainable by removing &
                 existing chains.'
        call error_abort('Aborting') 
    end if

    ! Print actual concentration once chains are removed
    concentration = real(nmonomers*proc_chains(irank)) / real(fluid_np)
    if (irank .eq. 1) then
        print*, 'Actual concentration: ', concentration
    end if

    ! Shift the z-positions of chain monomers so that they are all separated
    ! by their equilibrium distance.
    call contract_chains_to_equil_sep

    ! Relabel chainIDs globally
    call globalSum(proc_chains,nproc)
    call globalSum(proc_nps,nproc)
    do n=1,np
        if (monomer(n)%chainID.ne.0) then
            monomer(n)%chainID = monomer(n)%chainID + sum(proc_chains(1:irank)) - proc_chains(irank)
        end if
        monomer(n)%glob_no = monomer(n)%glob_no + sum(proc_nps(1:irank)) - proc_nps(irank)
    end do

    nchains = sum(proc_chains)

    !if (irank .eq. 1) then 
    !do n = 1, np
    !    print*, '-n,c,sc,f,g------------------------------------------------'
    !    print*, n, monomer(n)%chainID,monomer(n)%subchainID,monomer(n)%funcy,monomer(n)%glob_no
    !    print*, 'bin_bflag -', monomer(n)%bin_bflag
    !    print*, '==========================================================='
    !end do 
    !end if

contains

    subroutine initialise_info
        implicit none

        !integer :: intbits
        real(kind(0.d0)) :: check, rpol(3)
        character(256) :: string

        proc_chains(:)         = 0
        proc_nps(:)            = 0

        ! Check chains will fit
        check = (4.0*real(initialnunits(3))/real(npz))/real(nmonomers)
        if (check.lt.1.d0) then
            print*, 'np, nmonomers, nmols in z direction =  ', np, &
                     nmonomers, 4*initialnunits(3)/npz 
            string = 'Number of molecules in the z direction on each processor &
                      must not be more than 4*ncells(3). Please & 
                      change the chain length in the input file.'
            call error_abort(string)
        end if
    
        ! Check that density is not too low to connect molecules
        !check = 0.5*sqrt(initialunitsize(1)**2.0 + initialunitsize(2)**2.0 + &
        !                 initialunitsize(3)**2.0)
        !if (check .gt. R_0) then
        !    print*, 'Density is too low for polymer connections to be made:  &
        !            the FCC lattice atoms are separated by', check,' which is&
        !            more than the polymer bond maximum extension:', R_0
        !    call error_abort('Aborting.')
        !end if

        ! Count cylinder mols and initialise monomer info
        fluid_np = 0
        cyl_np = 0 
        do n=1,np+extralloc

            ! Count cylinder and fluid nps
            if (n .le. np) then
                rpol(:) = cpolariser(globalise(r(:,n)))
                if (rpol(1) .lt. r_oi .or. &
                    rpol(1) .gt. r_io) then 
                     cyl_np = cyl_np + 1 
                else
                    fluid_np = fluid_np + 1
                end if
            end if 

            ! Initialise monomerinfo
            monomer(n)%chainID      = 0
            monomer(n)%subchainID   = 0
            monomer(n)%funcy        = 0
            monomer(n)%glob_no      = 0
            monomer(n)%bin_bflag(:) = 0

        end do  

    end subroutine 

    subroutine connect_all_possible_chains(maxchainID)
        implicit none

        integer, intent(out) :: maxchainID
        integer :: subchainID, chainID, n, molno

        chainID = 1
        n = 1
        do
            ! Exit if past np 
            if ( n .gt. np ) exit

            ! Test to see if we can build a whole chain
            call test_for_connectability(n, nmonomers, connectable)

            ! If possible, build whole chain, otherwise mark as solvent and move on 
            if (connectable) then

                ! Connect 1 to 2
                monomer(n)%chainID    = chainID
                monomer(n)%subchainID = 1 
                monomer(n)%glob_no    = n   
                call connect_to_monomer(2,n)

                ! Connect middles
                do subchainID = 2, nmonomers-1
                    molno = n+subchainID-1
                    monomer(molno)%chainID    = chainID
                    monomer(molno)%subchainID = subchainID 
                    monomer(molno)%glob_no    = molno !corrected at bottom
                    call connect_to_monomer(subchainID-1,molno) 
                    call connect_to_monomer(subchainID+1,molno)             
                end do

                ! Connect end-1 to end
                molno = n+nmonomers-1
                monomer(molno)%chainID    = chainID
                monomer(molno)%subchainID = subchainID 
                monomer(molno)%glob_no    = molno !corrected at bottom
                call connect_to_monomer(nmonomers-1,molno)
                chainID = chainID + 1

                n = n + nmonomers

            else

                monomer(n)%chainID     = 0
                monomer(n)%subchainID  = 1
                monomer(n)%glob_no     = n
                monomer(n)%funcy       = 0
                bond(:,n)              = 0

                n = n + 1

            end if

        end do

        maxchainID = chainID

    end subroutine

    subroutine remove_chains(nremove)
        implicit none

        integer, intent(in) :: nremove

        integer :: m, cnt, chainID
        integer, dimension(:), allocatable :: removeIDs

        ! Spread removed chains throughout proc domain
        allocate(removeIDs(nremove))
        removeIDs(:) = 0
        do n = 1, nremove
            chainID = floor(real(n*proc_chains(irank))/real(nremove))
            removeIDs(n) = chainID
        end do

        ! Remove chains labelled with removeIDs      
        do n = 1,nremove
            call mark_chain_as_solvent(removeIDs(n)) 
        end do

        ! Relabel remaining chains
        cnt = 0
        n = 1
        do 

            if (monomer(n)%chainID .gt. cnt) then

                do m = 0, nmonomers-1
                    monomer(n+m)%chainID = cnt + 1
                end do

                n = n + nmonomers
                cnt = cnt + 1

            else
                
                n = n + 1

            end if
            
            if (n .gt. np) exit

        end do

        deallocate(removeIDs)

    end subroutine

    subroutine contract_chains_to_equil_sep 
        implicit none

        integer :: m
        real(kind(0.d0)) :: rscale
        real(kind(0.d0)) :: rmiddle(3)
        real(kind(0.d0)) :: oldsep, equil_sep

        equil_sep = 0.9608971929802091
        if (irank.eq.iroot) then
            print*, "Warning: equilibrium separation distance of FENE chain "&
            "set to ", equil_sep, ", based on R_0 = 1.5 and k = 30, with LJ "&
            "cutoff 2^1/6"
        end if

        n = 1
        do 

            if (n .gt. np) exit

            if (monomer(n)%chainID .ne. 0) then

                oldsep = sqrt(dot_product(r(:,n+1) - r(:,n),r(:,n+1) - r(:,n)))
                rmiddle(:) = (r(:,n) + r(:,n+nmonomers-1))/2.d0
                rscale = equil_sep / oldsep 
                do m = 0, nmonomers-1
                    r(:,n+m) = rmiddle(:) + rscale*(r(:,n+m) - rmiddle(:)) 
                end do

                n = n + nmonomers

            else
                
                n = n + 1

            end if

        end do


    end subroutine

    subroutine mark_chain_as_solvent(ID)
        implicit none

        integer, intent(in) :: ID

        integer :: molno

        do molno = 1, np

            if (monomer(molno)%chainID .eq. ID) then
                monomer(molno)%chainID = 0
                monomer(molno)%subchainID = 1
                monomer(molno)%funcy = 0
                monomer(molno)%bin_bflag(:) = 0
                bond(:,molno) = 0
            end if

        end do

    end subroutine mark_chain_as_solvent

    subroutine test_for_connectability(molno, seekahead, success)
        implicit none
        
        integer, intent(in) :: molno, seekahead
        logical, intent(out) :: success

        integer :: m
        real(kind(0.d0)) :: rn(3), rm(3), rnm(3), rnmmag
        real(kind(0.d0)) :: rpol(3)

        success = .true.

        if (molno + seekahead - 1 .gt. np) then
            success = .false. 
        end if

        do m = 0,seekahead-1
            
            rn = r(:,molno+m) 
            rm = r(:,molno+m+1)
            rnm = rm - rn 
            rnmmag = sqrt(dot_product(rnm,rnm))

            ! Check that the molecule isn't a cylinder one
            rpol(:) = cpolariser(globalise(rm))
            if (rpol(1) .lt. r_oi) success = .false.! ; return
            if (rpol(1) .gt. r_io) success = .false.! ; return 

            ! Check if molecule is in domain of processor
            if (abs(rm(1)) .ge.  halfdomain(1) .or. &
                abs(rm(2)) .ge.  halfdomain(2) .or. &
                abs(rm(3)) .ge.  halfdomain(3)) then
                success = .false. 
            end if

            ! Check if molecule pair is in bond range
            if (rnmmag .ge. R_0) success = .false. !; return

        end do


    end subroutine test_for_connectability

end subroutine setup_cylinder_FENE_solution


subroutine setup_initialise_sparse_FENE
    use computational_constants_MD, only: domain,globaldomain,halfdomain,  &
                                          irank,iroot,potential_flag,      &
                                          iblock,jblock,kblock, extralloc, &
                                          nproc, npy
    use polymer_info_MD, only: nchains, nmonomers,monomer,intbits
    use physical_constants_MD, only: np,globalnp,rcutoff
    use interfaces, only: error_abort
    use arrays_MD, only: r
    use messenger_data_exchange, only : globalSum
#if USE_COUPLER
    use coupler
    use md_coupler_socket, only: socket_get_domain_top
#endif
    implicit none

    integer :: maxmonx, maxmony, maxmonz, maxmonxyz
    integer :: n, chainID, subchainID
    integer, dimension(:), allocatable :: proc_chains, proc_nps
    real(kind(0.d0)) :: rc(3)
    real(kind(0.d0)) :: equil_sep
    real(kind(0.d0)) :: dir(3)
    real(kind(0.d0)) :: xL,yL,zL,domain_top
    character(len=200) :: string

    !intbits = bit_size(monomer(1)%bin_bflag(1))

    !Set top of domain initially
    domain_top = domain(2)/2.d0

#if USE_COUPLER

    if (jblock .eq. npy) then
        domain_top = socket_get_domain_top()
    endif

#endif

    ! Allocate arrays to store np and nchains per proc
    allocate(proc_chains(nproc))
    allocate(proc_nps(nproc))
    proc_chains(:) = 0
    proc_nps(:) = 0

    ! Initialise monomer data
    do n=1,np+extralloc
        monomer(n)%chainID      = 0
        monomer(n)%subchainID   = 0
        monomer(n)%glob_no      = 0
        monomer(n)%funcy        = 0
        monomer(n)%bin_bflag(:) = 0
    end do  

    ! Use equil sep as cubic lattice parameter
    equil_sep = 0.9608971929802091
    if (irank.eq.iroot) then
        print*, "Warning: equilibrium separation distance of FENE chain set to", &
        equil_sep, ", based on R_0 = 1.5 and k = 30, with LJ cutoff 2^1/6"
    end if

    ! Max num of monomers separated by equil_sep in x,y,z   
    xL = domain(1)
    yL = domain(2) - (halfdomain(2) - domain_top) 
    zL = domain(3)
    maxmonx = xL / equil_sep        
    maxmony = yL / equil_sep        
    maxmonz = zL / equil_sep 

    ! Max monomers in proc volume
    maxmonxyz = maxmonx * maxmony * maxmonz 
    ! Set number of particles per processor
    np = floor(real(nchains)/real(nproc))*nmonomers

    if (maxmonxyz .lt. np) then
        string = "Domain is incorrectly sized for the specified number " &
                  // "of chains and monomers per chain."
        call error_abort(string)
    end if

    ! Start in bottom corner and walk on cubic lattice, turning
    ! at the walls like a snake!
    rc(:) = 0.5*equil_sep 
    dir(:) = 1.d0
    chainID = 1
    subchainID = 0
    do n=1,np
    
        subchainID = subchainID + 1
        !Correct to -halfdomain to halfdomain coords
        r(:,n) = rc(:) - halfdomain(:)
        if (subchainID .gt. nmonomers) then
            subchainID = 1
            chainID    = chainID + 1
        end if
        monomer(n)%chainID      = chainID 
        monomer(n)%subchainID   = subchainID
        monomer(n)%glob_no      = n
        if (subchainID.eq.1) then
            call connect_to_monomer(subchainID+1,n)
        else if (subchainID.eq.nmonomers) then
            call connect_to_monomer(subchainID-1,n)             
        else
            call connect_to_monomer(subchainID+1,n)
            call connect_to_monomer(subchainID-1,n)
        end if

        !Walk rc
        rc(2)  = rc(2) + dir(2)*equil_sep
        ! Check if turn needed in y direction
        if ( dir(2) .gt. 0 ) then
            if ( rc(2) .ge. domain_top - 0.5*equil_sep) then 
                rc(2)  = rc(2) - equil_sep
                rc(1)  = rc(1) + dir(1)*equil_sep
                dir(2) = dir(2) * -1.d0
            end if
        else 
            if ( rc(2) .lt. 0.5*equil_sep) then 
                rc(2) = rc(2) + equil_sep
                rc(1) = rc(1) + dir(1)*equil_sep
                dir(2) = dir(2) * -1.d0
            end if
        end if
        ! Check if turn needed in x direction
        if ( dir(1) .gt. 0 ) then
            if ( rc(1) .ge. domain(1)-0.5*equil_sep) then
                rc(1) = rc(1) - equil_sep
                rc(3) = rc(3) + equil_sep
                dir(1) = dir(1) * -1.d0
            end if
        else
            if ( rc(1) .lt. 0.5*equil_sep) then
                rc(1) = rc(1) + equil_sep
                rc(3) = rc(3) + equil_sep
                dir(1) = dir(1) * -1.d0
            end if
        end if              

        !!if (irank.eq.2) then
        !print*, '-n,c,sc,f,g------------------------------------------------'
        !print*, n, monomer(n)%chainID,monomer(n)%subchainID,monomer(n)%funcy,monomer(n)%glob_no
        !print*, 'bin_bflag -', monomer(n)%bin_bflag
        !print*, '==========================================================='
        !!end if

    end do

    !Establish global number of particles on current process
    globalnp = np
    call globalSum(globalnp)

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp

    proc_chains(irank) = chainID
    proc_nps(irank)    = np

    call globalSum(proc_chains,nproc)
    call globalSum(proc_nps,nproc)
    
    do n=1,np
        if (monomer(n)%chainID.ne.0) then
            monomer(n)%chainID = monomer(n)%chainID + sum(proc_chains(1:irank)) - proc_chains(irank)
        end if
        monomer(n)%glob_no = monomer(n)%glob_no + sum(proc_nps(1:irank))    - proc_nps(irank)
    end do

#if USE_COUPLER

    if (jblock .eq. npy .and. iblock .eq. 1 .and. kblock .eq. 1) then
        print*, '*********************************************************************'
        print*, '*WARNING - TOP LAYER OF DOMAIN REMOVED IN LINE WITH CONSTRAINT FORCE*'
        print*, 'Removed from', domain_top, 'to Domain top', globaldomain(2)/2.d0
        print*, '*********************************************************************'
    endif

#endif

end subroutine setup_initialise_sparse_FENE

!-----------------------------------------------------------------------------
!Connect monomer to subchainID bscID
subroutine connect_to_monomer(bscID,n)
use polymer_info_MD
implicit none


    integer, intent(in) :: bscID,n
    integer :: group, expo

    !intbits = bit_size(monomer(1)%bin_bflag(1))
    group = ceiling(real(bscID)/real(intbits))
    expo  = mod(bscID,intbits) - 1
    if(expo.eq.-1) expo = intbits - 1

    monomer(n)%funcy            = monomer(n)%funcy + 1
    monomer(n)%bin_bflag(group) = monomer(n)%bin_bflag(group) + 2**(expo)
    
end subroutine connect_to_monomer

subroutine setup_initialise_polymer_brush
    use module_initialise_microstate
    use polymer_info_MD
    use messenger
    use interfaces, only: error_abort
    use messenger_data_exchange, only : globalSum

    double precision :: solid_bottom(3), solid_top(3)
    double precision :: grafting_density_real, density_ratio
    integer :: nchainsremove, proc_units_xz
    integer :: wall_np, fluid_np, maxchainID
    integer, dimension(nproc) :: proc_chains, proc_nps

    ! Set positions on a lattice, connectable increasing pos in y
    ! Checks and intialisations
    ! Connect all the chains we can, bonds to be removed later 
    ! Store maximum number of chains
    solid_density = density
    call initialise_lattice_positions
    call initialise_info

    ! Remove chains to get target grafting density (as close as possible) 
    if (jblock .eq. 1) then
        call connect_all_possible_chains(maxchainID)
        proc_chains(irank) = maxchainID
        if (irank.eq.1) then
            print('(a,f10.5)'), 'Target grafting density: ', grafting_density
        end if
        grafting_density_real = 1.d0 
        nchainsremove = nint((grafting_density_real - grafting_density) &
                              *real(proc_chains(irank)))
        if (nchainsremove .gt. 0) then
            call remove_chains_random(nchainsremove)
            proc_chains(irank) = proc_chains(irank) - nchainsremove 
        else if (nchainsremove .eq. 0) then
            print*, "Don't need to remove any chains on irank:", irank
        else
            print*, 'Grafting density before chain removal:', grafting_density_real 
            print*, 'Target grafting density, nchains connectable: ', &
                     grafting_density, proc_chains(irank) 
            print*, 'The target grafting density is not obtainable by removing &
                     existing chains.'
            call error_abort('Aborting') 
        end if
        proc_units_xz = initialnunits(1)*initialnunits(3)/(npx*npz) 
        grafting_density = real(proc_chains(irank))/real(proc_units_xz)
        if (irank .eq. 1) then
            print('(a,f10.5)'), 'Actual grafting density: ', grafting_density 
        end if
        ! Shift the y-positions of chain monomers so that they are all separated
        ! by their equilibrium distance.
        call contract_chains_to_equil_sep
    else
        maxchainID = 0
        proc_chains(irank) = 0
    end if

    ! Remove solvent molecules so that target density is acquired, store 
    ! new number of particles
    density_ratio = liquid_density/solid_density
    nmolsremove = nint(real(np)*(1.d0 - density_ratio))
    !print*, irank, 'before'
    call remove_solvent_mols(nmolsremove)
    !print*, irank, 'after'
    proc_nps(irank) = np

    ! Relabel chainIDs globally
    call globalSum(proc_chains,nproc)
    call globalSum(proc_nps,nproc)
    do n=1,np
        if (monomer(n)%chainID.ne.0) then
            monomer(n)%chainID = monomer(n)%chainID + sum(proc_chains(1:irank)) - proc_chains(irank)
        end if
        monomer(n)%glob_no = monomer(n)%glob_no + sum(proc_nps(1:irank)) - proc_nps(irank)
    end do
    nchains = sum(proc_chains)

    !Establish global number of particles on current process
    globalnp = np
    call globalSum(globalnp)

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp

#if USE_COUPLER

    if (jblock .eq. npy .and. iblock .eq. 1 .and. kblock .eq. 1) then
        print*, '*********************************************************************'
        print*, '*WARNING - TOP LAYER OF DOMAIN REMOVED IN LINE WITH CONSTRAINT FORCE*'
        print*, 'Removed from', domain_top, 'to Domain top', globaldomain(2)/2.d0
        print*, 'Number of molecules reduced from',  & 
                 4*initialnunits(1)*initialnunits(2)*initialnunits(3), 'to', globalnp
        print*, '*********************************************************************'
        !print*, 'microstate ', minval(r(1,:)), maxval(r(1,:)),minval(r(2,:)),  &
        !         maxval(r(2,:)),minval(r(3,:)), maxval(r(3,:))
    endif

#endif

contains

    subroutine initialise_lattice_positions
#if USE_COUPLER
        use coupler
        use md_coupler_socket, only: socket_get_domain_top
#endif
        implicit none

        integer :: j, n, nl, nx, ny, nz
        integer, dimension(nd) :: p_units_lb, p_units_ub 
        double precision :: domain_top, solid_density, density_ratio
        double precision, dimension (nd):: rc, c

        p_units_lb(1) = (iblock-1)*floor(initialnunits(1)/real((npx),kind(0.d0)))
        p_units_ub(1) =  iblock *ceiling(initialnunits(1)/real((npx),kind(0.d0)))
        p_units_lb(2) = (jblock-1)*floor(initialnunits(2)/real((npy),kind(0.d0)))
        p_units_ub(2) =  jblock *ceiling(initialnunits(2)/real((npy),kind(0.d0)))
        p_units_lb(3) = (kblock-1)*floor(initialnunits(3)/real((npz),kind(0.d0)))
        p_units_ub(3) =  kblock *ceiling(initialnunits(3)/real((npz),kind(0.d0)))

        !Set top of domain initially
        domain_top = globaldomain(2)/2.d0

        !Setup solid/liquid properties
        solid_density = density
        density_ratio = liquid_density/solid_density

#if USE_COUPLER

        if (jblock .eq. npy) then
            domain_top = socket_get_domain_top()
        endif

#endif

        !Molecules per unit FCC structure (3D)
        n  = 0      !Initialise global np counter n
        nl = 0      !Initialise local np counter nl

        !Inner loop in y (useful for setting connectivity)
        do nz=p_units_lb(3),p_units_ub(3)
        c(3) = (nz - 0.75d0)*initialunitsize(3) !- halfdomain(3) 
        do nx=p_units_lb(1),p_units_ub(1)
        c(1) = (nx - 0.75d0)*initialunitsize(1) !- halfdomain(1)
        do ny=p_units_lb(2),p_units_ub(2)
        c(2) = (ny - 0.75d0)*initialunitsize(2) !- halfdomain(2) 

            do j=1,4    !4 Molecules per cell

                rc(:) = c(:)
                select case(j)
                case(2)
                    rc(1) = c(1) + 0.5d0*initialunitsize(1)
                    rc(3) = c(3) + 0.5d0*initialunitsize(3)
                case(3)
                    rc(2) = c(2) + 0.5d0*initialunitsize(2)
                    rc(3) = c(3) + 0.5d0*initialunitsize(3)
                case(4)
                    rc(1) = c(1) + 0.5d0*initialunitsize(1)
                    rc(2) = c(2) + 0.5d0*initialunitsize(2)
                case default
                end select

                n = n + 1   !Move to next particle
                
                !Remove molecules from top of domain if constraint applied
                if (jblock .eq. npy) then
                    ! Note rc is in global coordinates from 0 to globaldomain while Domaintop
                    ! is given in global coordinates from -halfglobaldomain to halfglobaldomain
                    if (rc(2)-globaldomain(2)/2.d0 .gt. domain_top) cycle 
                endif

                !Check if molecule is in domain of processor
                if(rc(1).lt. domain(1)*(iblock-1)) cycle
                if(rc(1).ge. domain(1)*(iblock  )) cycle
                if(rc(2).lt. domain(2)*(jblock-1)) cycle
                if(rc(2).ge. domain(2)*(jblock  )) cycle
                if(rc(3).lt. domain(3)*(kblock-1)) cycle
                if(rc(3).ge. domain(3)*(kblock  )) cycle

                !Solid region given by wall textures or tethered region
                call wall_textures(texture_type,(rc(:)-globaldomain(:)/2.d0),solid_bottom,solid_top)

                !If molecules is in the domain then add to total
                nl = nl + 1 !Local molecule count

                !Correct to local coordinates
                r(1,nl) = rc(1)-domain(1)*(iblock-1)-halfdomain(1)
                r(2,nl) = rc(2)-domain(2)*(jblock-1)-halfdomain(2)
                r(3,nl) = rc(3)-domain(3)*(kblock-1)-halfdomain(3)

            enddo

        enddo
        enddo
        enddo

        !Correct local number of particles on processor
        np = nl

    end subroutine initialise_lattice_positions

    subroutine initialise_info
        implicit none

        integer :: n
        real(kind(0.d0)) :: check 
        character(256) :: string

        proc_chains(:)         = 0
        proc_nps(:)            = 0

        ! Check chains will fit
        check = (4.0*real(initialnunits(2))/real(npy))/real(nmonomers)
        if (check.lt.1.d0) then
            print*, 'np, nmonomers, nmols in y direction =  ', np, &
                     nmonomers, 4*initialnunits(2)/npy
            string = 'Number of molecules in the y direction on each processor &
                      must not be more than 4*ncells(2). Please & 
                      change the chain length in the input file.'
            call error_abort(string)
        end if

        ! Count cylinder mols and initialise monomer info
        fluid_np = 0
        wall_np = 0 
        do n=1,np+extralloc

            ! Count cylinder and fluid nps
            if (n .le. np) then

                if (any(r(:,n) .lt. solid_bottom(:)) &
                    .or. any(r(:,n) .gt. solid_top(:))) then

                    wall_np = wall_np + 1 

                else

                    fluid_np = fluid_np + 1

                end if 

            end if 

            ! Initialise monomerinfo
            monomer(n)%chainID      = 0
            monomer(n)%subchainID   = 1
            monomer(n)%funcy        = 0
            monomer(n)%glob_no      = n
            monomer(n)%bin_bflag(:) = 0

        end do  

    end subroutine 

    subroutine connect_all_possible_chains(maxchainID)
        implicit none

        integer, intent(out) :: maxchainID
        integer :: subchainID, chainID, n, molno
        logical :: connectable

        chainID = 1
        n = 1
        do
           ! Exit if past np 
            if ( n .gt. np ) exit

            ! Test to see if we can build a whole chain
            call test_for_connectability(n, nmonomers, connectable)

            ! If possible, build whole chain, otherwise mark as solvent and move on 
            if (connectable) then

                ! Connect 1 to 2
                monomer(n)%chainID    = chainID
                monomer(n)%subchainID = 1 
                monomer(n)%glob_no    = n   
                call connect_to_monomer(2,n)

                ! Connect middles
                do subchainID = 2, nmonomers-1
                    molno = n+subchainID-1
                    monomer(molno)%chainID    = chainID
                    monomer(molno)%subchainID = subchainID 
                    monomer(molno)%glob_no    = molno !corrected at bottom
                    call connect_to_monomer(subchainID-1,molno) 
                    call connect_to_monomer(subchainID+1,molno)             
                end do

                ! Connect end-1 to end
                molno = n+nmonomers-1
                monomer(molno)%chainID    = chainID
                monomer(molno)%subchainID = subchainID 
                monomer(molno)%glob_no    = molno !corrected at bottom
                call connect_to_monomer(nmonomers-1,molno)

                chainID = chainID + 1
                n = n + nmonomers

            else

                monomer(n)%chainID     = 0
                monomer(n)%subchainID  = 1
                monomer(n)%glob_no     = n
                monomer(n)%funcy       = 0
                bond(:,n)              = 0

                n = n + 1

            end if

        end do

        maxchainID = chainID - 1

    end subroutine

    subroutine remove_chains_random(nremove)
        implicit none

        integer, intent(in) :: nremove

        integer :: m, cnt, chainID
        integer, dimension(:), allocatable :: removeIDs, removeflags
        double precision :: random

        allocate(removeIDs(nremove))
        allocate(removeflags(proc_chains(irank)))
        removeflags(:) = 0
        cnt = 0
        do 
            call random_number(random)
            chainID = ceiling(random*proc_chains(irank))
            if (removeflags(chainID) .eq. 0) then
                cnt = cnt + 1
                removeflags(chainID) = 1
                removeIDs(cnt) = chainID 
            end if

            if (cnt .eq. nremove) then
                exit
            end if  

        end do 

        ! Remove chains labelled with removeIDs      
        do n = 1,nremove
            call mark_chain_as_solvent(removeIDs(n)) 
        end do

        ! Relabel remaining chains
        cnt = 0
        n = 1
        do 

            if (monomer(n)%chainID .gt. cnt) then

                do m = 0, nmonomers-1
                    monomer(n+m)%chainID = cnt + 1
                end do

                n = n + nmonomers
                cnt = cnt + 1

            else
                
                n = n + 1

            end if
            
            if (n .gt. np) exit

        end do

        deallocate(removeIDs)

    end subroutine

    subroutine contract_chains_to_equil_sep 
       implicit none

        integer :: m
        real(kind(0.d0)) :: rscale
        real(kind(0.d0)) :: rfirst(3)
        real(kind(0.d0)) :: oldsep, equil_sep

        equil_sep = 0.9608971929802091
        if (irank.eq.iroot) then
            print*, "Warning: equilibrium separation distance of FENE chain "&
            "set to ", equil_sep, ", based on R_0 = 1.5 and k = 30, with LJ "&
            "cutoff 2^1/6"
        end if

        n = 1
        do 

            if (n .gt. np) exit

            if (monomer(n)%chainID .ne. 0) then

                oldsep = sqrt(dot_product(r(:,n+1) - r(:,n),r(:,n+1) - r(:,n)))
                rfirst(:) = r(:,n)
                rscale = equil_sep / oldsep 
                do m = 0, nmonomers-1
                    r(:,n+m) = rfirst(:) + rscale*(r(:,n+m) - rfirst(:)) 
                end do

                n = n + nmonomers

            else
                
                n = n + 1

            end if

        end do


    end subroutine

    subroutine mark_chain_as_solvent(ID)
        implicit none

        integer, intent(in) :: ID

        integer :: molno

        do molno = 1, np

            if (monomer(molno)%chainID .eq. ID) then
                monomer(molno)%chainID = 0
                monomer(molno)%subchainID = 1
                monomer(molno)%funcy = 0
                monomer(molno)%bin_bflag(:) = 0
                bond(:,molno) = 0
            end if

        end do

    end subroutine mark_chain_as_solvent

    subroutine test_for_connectability(molno, seekahead, success)
        implicit none
        
        integer, intent(in) :: molno, seekahead
        logical, intent(out) :: success

        integer :: m
        real(kind(0.d0)) :: rn(3), rm(3), rnm(3), rnmmag

        success = .true.

        ! If first molecule is in fluid region, don't connect
        rn = r(:,molno)
        call wall_textures(texture_type,globalise(rn),solid_bottom,solid_top)
        if (all(globalise(rn) .ge. solid_bottom-globaldomain/2.d0) .and. & 
            all(globalise(rn) .le. (globaldomain/2.d0)-solid_top)) then
            success = .false.
        end if

        if (molno + seekahead - 1 .gt. np) then
            success = .false. 
        end if

        do m = 0,seekahead-1
            
            rn = r(:,molno+m) 
            rm = r(:,molno+m+1)
            rnm = rm - rn 
            rnmmag = sqrt(dot_product(rnm,rnm))

            ! Check that the molecule isn't a cylinder one
            call wall_textures(texture_type,globalise(rm),solid_bottom,solid_top)
            if (any(globalise(rm) .lt. solid_bottom-globaldomain/2.d0)) success = .false.! ; return
            if (any(globalise(rm) .gt. (globaldomain/2.d0)-solid_top)) success = .false.! ; return 

            ! Check if molecule is in domain of processor
            if (abs(rm(1)) .ge.  halfdomain(1) .or. &
                abs(rm(2)) .ge.  halfdomain(2) .or. &
                abs(rm(3)) .ge.  halfdomain(3)) then
                success = .false. 
            end if

            ! Check if molecule pair is in bond range
            if (rnmmag .ge. R_0) success = .false. !; return

        end do
        
    end subroutine test_for_connectability


    subroutine remove_solvent_mols(nremove)
        implicit none
   
        integer, intent(in) :: nremove
        
        integer :: n, molno, failcount
        double precision :: random

        n = 0
        failcount = 0
        do 

            if (failcount .gt. 100*np) then
                print*, "Can't remove enough molecules from irank ", irank
                call error_abort()
            end if

            failcount = failcount + 1
            if (n .eq. nremove) exit

            call random_number(random)
            molno = ceiling(random*real(np))

            call wall_textures(texture_type,globalise(r(:,molno)),solid_bottom,solid_top)
            if (any(globalise(r(:,molno)) .lt. solid_bottom-globaldomain/2.d0)) cycle
            if (any(globalise(r(:,molno)) .gt. (globaldomain/2.d0)-solid_top)) cycle

            if (monomer(molno)%chainID .eq. 0) then 

                n = n + 1
                r(:,molno) = r(:,np)
                monomer(molno) = monomer(np)
                monomer(molno)%glob_no = molno
                tag(molno) = tag(np)
                np = np - 1 
                failcount = 0
            
            end if

        
        end do 

    end subroutine remove_solvent_mols

end subroutine setup_initialise_polymer_brush

!-----------------------------------------------------------------------------
subroutine setup_initialise_solid_liquid
    use physical_constants_MD, only : fixdistbottom
    use module_initialise_microstate
    use messenger
    use messenger_data_exchange, only : globalSum
#if USE_COUPLER
    use coupler
    use md_coupler_socket, only: socket_get_domain_top
#endif
    implicit none

    integer :: j, n, nl, nx, ny, nz
    integer, dimension(nd) :: p_units_lb, p_units_ub 
    double precision :: domain_top, solid_density, density_ratio
    double precision, dimension (nd):: solid_bottom,solid_top, rc, c

    p_units_lb(1) = (iblock-1)*floor(initialnunits(1)/real((npx),kind(0.d0)))
    p_units_ub(1) =  iblock *ceiling(initialnunits(1)/real((npx),kind(0.d0)))
    p_units_lb(2) = (jblock-1)*floor(initialnunits(2)/real((npy),kind(0.d0)))
    p_units_ub(2) =  jblock *ceiling(initialnunits(2)/real((npy),kind(0.d0)))
    p_units_lb(3) = (kblock-1)*floor(initialnunits(3)/real((npz),kind(0.d0)))
    p_units_ub(3) =  kblock *ceiling(initialnunits(3)/real((npz),kind(0.d0)))

    !Set top of domain initially
    domain_top = globaldomain(2)/2.d0

    !Setup solid/liquid properties
    solid_density = density
    density_ratio = liquid_density/solid_density

#if USE_COUPLER

    if (jblock .eq. npy) then
        domain_top = socket_get_domain_top()
    endif

#endif

    !Molecules per unit FCC structure (3D)
    n  = 0      !Initialise global np counter n
    nl = 0      !Initialise local np counter nl

    !Inner loop in y (useful for setting connectivity)
    do nz=p_units_lb(3),p_units_ub(3)
    c(3) = (nz - 0.75d0)*initialunitsize(3) !- halfdomain(3) 
    do nx=p_units_lb(1),p_units_ub(1)
    c(1) = (nx - 0.75d0)*initialunitsize(1) !- halfdomain(1)
    do ny=p_units_lb(2),p_units_ub(2)
    c(2) = (ny - 0.75d0)*initialunitsize(2) !- halfdomain(2) 

        do j=1,4    !4 Molecules per cell

            rc(:) = c(:)
            select case(j)
            case(2)
                rc(1) = c(1) + 0.5d0*initialunitsize(1)
                rc(3) = c(3) + 0.5d0*initialunitsize(3)
            case(3)
                rc(2) = c(2) + 0.5d0*initialunitsize(2)
                rc(3) = c(3) + 0.5d0*initialunitsize(3)
            case(4)
                rc(1) = c(1) + 0.5d0*initialunitsize(1)
                rc(2) = c(2) + 0.5d0*initialunitsize(2)
            case default
            end select

            n = n + 1   !Move to next particle
            
            !Remove molecules from top of domain if constraint applied
            if (jblock .eq. npy) then
                ! Note rc is in global coordinates from 0 to globaldomain while Domaintop
                ! is given in global coordinates from -halfglobaldomain to halfglobaldomain
                if (rc(2)-globaldomain(2)/2.d0 .gt. domain_top) cycle 
            endif

            !Check if molecule is in domain of processor
            if(rc(1).lt. domain(1)*(iblock-1)) cycle
            if(rc(1).ge. domain(1)*(iblock  )) cycle
            if(rc(2).lt. domain(2)*(jblock-1)) cycle
            if(rc(2).ge. domain(2)*(jblock  )) cycle
            if(rc(3).lt. domain(3)*(kblock-1)) cycle
            if(rc(3).ge. domain(3)*(kblock  )) cycle

            !Solid region given by wall textures or tethered region
            call wall_textures(texture_type,(rc(:)-globaldomain(:)/2.d0),solid_bottom,solid_top)

            !If outside solid region, randomly remove molecules from lattice
            if (rc(2)-globaldomain(2)/2.d0 .gt. (solid_bottom(2) - globaldomain(2)/2.d0) .and. & 
                rc(2)-globaldomain(2)/2.d0 .lt. (globaldomain(2)/2.d0  -  solid_top(2))) then
                !cycle
                call random_number(rand)
                if (rand .gt. density_ratio) cycle
            endif

            !If molecules is in the domain then add to total
            nl = nl + 1 !Local molecule count

            !Correct to local coordinates
            r(1,nl) = rc(1)-domain(1)*(iblock-1)-halfdomain(1)
            r(2,nl) = rc(2)-domain(2)*(jblock-1)-halfdomain(2)
            r(3,nl) = rc(3)-domain(3)*(kblock-1)-halfdomain(3)

        enddo

    enddo
    enddo
    enddo

    !Correct local number of particles on processor
    np = nl

    !Establish global number of particles on current process
    globalnp = np
    call globalSum(globalnp)

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp

#if USE_COUPLER

    if (jblock .eq. npy .and. iblock .eq. 1 .and. kblock .eq. 1) then
        print*, '*********************************************************************'
        print*, '*WARNING - TOP LAYER OF DOMAIN REMOVED IN LINE WITH CONSTRAINT FORCE*'
        print*, 'Removed from', domain_top, 'to Domain top', globaldomain(2)/2.d0
        print*, 'Number of molecules reduced from',  & 
                 4*initialnunits(1)*initialnunits(2)*initialnunits(3), 'to', globalnp
        print*, '*********************************************************************'
        !print*, 'microstate ', minval(r(1,:)), maxval(r(1,:)),minval(r(2,:)),  &
        !         maxval(r(2,:)),minval(r(3,:)), maxval(r(3,:))
    endif

#endif

end subroutine setup_initialise_solid_liquid

!=============================================================================
!Initialise "solid" concentric cylinders
subroutine setup_initialise_concentric_cylinders
    use module_initialise_microstate
    use concentric_cylinders
    use messenger
    use messenger_data_exchange, only : globalSum
    use interfaces, only: error_abort
    implicit none

    integer :: j, n, nl, nx, ny, nz
    integer :: p_units_lb(nd), p_units_ub(nd)
    real(kind(0.d0)) :: rr,rx,ry             !Radial pos (cylindrical polar)
    double precision, dimension (nd):: rc, c !Temporary variable

    p_units_lb(1) = (iblock-1)*floor(initialnunits(1)/real((npx),kind(0.d0)))
    p_units_ub(1) =  iblock *ceiling(initialnunits(1)/real((npx),kind(0.d0)))
    p_units_lb(2) = (jblock-1)*floor(initialnunits(2)/real((npy),kind(0.d0)))
    p_units_ub(2) =  jblock *ceiling(initialnunits(2)/real((npy),kind(0.d0)))
    p_units_lb(3) = (kblock-1)*floor(initialnunits(3)/real((npz),kind(0.d0)))
    p_units_ub(3) =  kblock *ceiling(initialnunits(3)/real((npz),kind(0.d0)))

    !Molecules per unit FCC structure (3D)
    n  = 0      !Initialise global np counter n
    nl = 0      !Initialise local np counter nl

    !Inner loop in y (useful for setting connectivity)
    do nz=p_units_lb(3),p_units_ub(3)
    c(3) = (nz - 0.75d0)*initialunitsize(3) !- halfdomain(3) 
    do nx=p_units_lb(1),p_units_ub(1)
    c(1) = (nx - 0.75d0)*initialunitsize(1) !- halfdomain(1)
    do ny=p_units_lb(2),p_units_ub(2)
    c(2) = (ny - 0.75d0)*initialunitsize(2) !- halfdomain(2) 

        do j=1,4    !4 Molecules per cell

            rc(:) = c(:)
            select case(j)
            case(2)
                rc(1) = c(1) + 0.5d0*initialunitsize(1)
                rc(3) = c(3) + 0.5d0*initialunitsize(3)
            case(3)
                rc(2) = c(2) + 0.5d0*initialunitsize(2)
                rc(3) = c(3) + 0.5d0*initialunitsize(3)
            case(4)
                rc(1) = c(1) + 0.5d0*initialunitsize(1)
                rc(2) = c(2) + 0.5d0*initialunitsize(2)
            case default
            end select

            n = n + 1   !Move to next particle

            !Check if molecule is in domain of processor
            if(rc(1).lt. domain(1)*(iblock-1)) cycle
            if(rc(1).ge. domain(1)*(iblock  )) cycle
            if(rc(2).lt. domain(2)*(jblock-1)) cycle
            if(rc(2).ge. domain(2)*(jblock  )) cycle
            if(rc(3).lt. domain(3)*(kblock-1)) cycle
            if(rc(3).ge. domain(3)*(kblock  )) cycle
             
            rx = rc(1) - globaldomain(1)/2.d0 !rc between 0->domain, 
            ry = rc(2) - globaldomain(2)/2.d0 !r between -half->half

            ! Cylindrical polar radial position
            rr = sqrt(rx*rx + ry*ry)

            ! Cycle if not in cylinder walls
            if (rr .lt. r_ii) cycle
            if (rr .gt. r_oo) cycle
            if (rr .gt. r_oi .and. rr .lt. r_io) cycle

            !If molecules is in the domain then add to total
            nl = nl + 1 !Local molecule count
            !Correct to local coordinates
            r(1,nl) = rc(1)-domain(1)*(iblock-1)-halfdomain(1)
            r(2,nl) = rc(2)-domain(2)*(jblock-1)-halfdomain(2)
            r(3,nl) = rc(3)-domain(3)*(kblock-1)-halfdomain(3)

        enddo

    enddo
    enddo
    enddo

    !Correct local number of particles on processor
    np = nl

    !Establish global number of particles on current process
    globalnp = np
    call globalSum(globalnp)

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp

end subroutine setup_initialise_concentric_cylinders

subroutine setup_initialise_fill_cylinders
    use module_initialise_microstate
    use concentric_cylinders
    use messenger
    use messenger_data_exchange, only : globalSum
    use boundary_MD, only: specular_flag, specular_flat, specular_wall
    use interfaces, only: error_abort
    implicit none

    integer :: j, nl, nx, ny, nz
    integer :: p_units_lb(nd), p_units_ub(nd)
    real(kind(0.d0)) :: rr,rx,ry,dr          !Radial pos (cylindrical polar)
    double precision, dimension (nd):: rc, c !Temporary variable

    dr = rcutoff ! no-overlap tolerance

    p_units_lb(1) = (iblock-1)*floor(initialnunits(1)/real((npx),kind(0.d0)))
    p_units_ub(1) =  iblock *ceiling(initialnunits(1)/real((npx),kind(0.d0)))
    p_units_lb(2) = (jblock-1)*floor(initialnunits(2)/real((npy),kind(0.d0)))
    p_units_ub(2) =  jblock *ceiling(initialnunits(2)/real((npy),kind(0.d0)))
    p_units_lb(3) = (kblock-1)*floor(initialnunits(3)/real((npz),kind(0.d0)))
    p_units_ub(3) =  kblock *ceiling(initialnunits(3)/real((npz),kind(0.d0)))

    !Molecules per unit FCC structure (3D)
    nl = 0      !Initialise local np counter nl

    !Inner loop in z (useful for setting connectivity)
    do nx=p_units_lb(1),p_units_ub(1)
    c(1) = (nx - 0.75d0)*initialunitsize(1) !- halfdomain(1)
    do ny=p_units_lb(2),p_units_ub(2)
    c(2) = (ny - 0.75d0)*initialunitsize(2) !- halfdomain(2) 
    do nz=p_units_lb(3),p_units_ub(3)
    c(3) = (nz - 0.75d0)*initialunitsize(3) !- halfdomain(3) 

        do j=1,4    !4 Molecules per cell

            rc(:) = c(:)
            select case(j)
            case(2)
                rc(1) = c(1) + 0.5d0*initialunitsize(1)
                rc(2) = c(2) + 0.5d0*initialunitsize(2)
            case(3)
                rc(1) = c(1) + 0.5d0*initialunitsize(1)
                rc(3) = c(3) + 0.5d0*initialunitsize(3)
            case(4)
                rc(2) = c(2) + 0.5d0*initialunitsize(2)
                rc(3) = c(3) + 0.5d0*initialunitsize(3)
            case default
            end select

            !Check if molecule is in domain of processor
            if(rc(1).lt. domain(1)*(iblock-1)) cycle
            if(rc(1).ge. domain(1)*(iblock  )) cycle
            if(rc(2).lt. domain(2)*(jblock-1)) cycle
            if(rc(2).ge. domain(2)*(jblock  )) cycle
            if(rc(3).lt. domain(3)*(kblock-1)) cycle
            if(rc(3).ge. domain(3)*(kblock  )) cycle
    
            rx = rc(1) - globaldomain(1)/2.d0 !rc between 0->domain, 
            ry = rc(2) - globaldomain(2)/2.d0 !r between -half->half

            ! Cylindrical polar radial position
            rr = sqrt(rx*rx + ry*ry)

            ! Cycle if not between cylinders
            if (rr .lt. r_oi + dr) cycle
            if (rr .gt. r_io - dr) cycle

            if (specular_flag .eq. specular_flat) then
                if (rc(3) .lt. specular_wall(3)) cycle
                if (rc(3) .ge. globaldomain(3) - specular_wall(3)) cycle
            end if

            !If molecules is in the domain then add to total
            nl = nl + 1 !Local molecule count
            !Correct to local coordinates
            r(1,np+nl) = rc(1)-domain(1)*(iblock-1)-halfdomain(1)
            r(2,np+nl) = rc(2)-domain(2)*(jblock-1)-halfdomain(2)
            r(3,np+nl) = rc(3)-domain(3)*(kblock-1)-halfdomain(3)

        enddo

    enddo
    enddo
    enddo

    ! Correct local number of particles on processor to now include
    ! inner fluid molecules
    np = np + nl

    !Establish global number of particles on current process
    globalnp = np
    call globalSum(globalnp)

    !Build array of number of particles on neighbouring
    !processe's subdomain on current proccess
    call globalGathernp

end subroutine setup_initialise_fill_cylinders

!=============================================================================
!Initialise branched polymer simulation
subroutine setup_initialise_polyinfo_singlebranched
    use interfaces
    use polymer_info_MD
    use messenger
    use messenger_data_exchange, only : globalSum
    use physical_constants_MD, only: np
    use arrays_MD, only:r
    implicit none

    integer :: n
    integer :: chainID
    integer :: subchainID
    integer :: modcheck
    integer :: scIDbranch, glob_n_branch, nbranch, branchmonomers
    integer, dimension(nproc) :: proc_chains, proc_nps
    double precision, dimension(3) :: rij
    double precision :: rij2
    
    proc_chains(:)         = 0
    proc_nps(:)            = 0
    branchmonomers         = nmonomers    !Branch length specified in input
    nmonomers              = 2*nmonomers  !Two branches per polymer

    modcheck = 0 + mod(np,nmonomers) + mod(4*initialnunits(1),nmonomers)
    if (modcheck.ne.0) call error_abort('Number of molecules must be exactly &
                                        &divisible by the polymer chain      &
                                        &length. Please change the chain     &
                                        &length in the input file. A chain   &
                                        &length of 4 should (hopefully)      &
                                        &work.')

    !Set all to solvents
    do n=1,np+extralloc
        monomer(n)%chainID      = 0
        monomer(n)%subchainID   = 1
        monomer(n)%funcy        = 0
        monomer(n)%glob_no      = n
        monomer(n)%bin_bflag(:) = 0
    end do  

    scIDbranch = branchmonomers/2            ! Choose subchainID from which
                                             ! to create a branch
    chainID    = 1                           ! Initialise
    subchainID = 0                           ! Initialise
    !First branch
    do n=1,branchmonomers                    ! Loop over first branch monomers

        chainID = 1                          ! Set all polymers to same chain
        subchainID = subchainID + 1          
        monomer(n)%chainID = chainID
        monomer(n)%subchainID = subchainID
        if (subchainID.eq.1) then            ! The usual chain connections...
            call connect_to_monomer(subchainID+1,n)
        else if (subchainID.eq.branchmonomers) then
            call connect_to_monomer(subchainID-1,n)             
        else
            call connect_to_monomer(subchainID+1,n)
            call connect_to_monomer(subchainID-1,n)
        end if

        if (subchainID .eq. scIDbranch) then ! Store global mol number of
            glob_n_branch = n                ! branching monomer
        end if

    end do

    !Locate new monomer close to branching monomer
    do n=branchmonomers+1,np
        rij(:) = r(glob_n_branch,:) - r(:,n) ! Check separation
        rij2 = dot_product(rij,rij)
        if (rij2.lt.R_0**2.d0) then          ! If within spring max elongation
            nbranch = n                      ! Store global mol no
            exit                             ! Nothing else needed
        end if
    end do

    !Second branch
    do n=nbranch,nbranch+branchmonomers-1    ! Start from monomer close to
        chainID = 1                          ! branching monomer. Same cID.
        subchainID = subchainID + 1
        monomer(n)%chainID = chainID
        monomer(n)%subchainID = subchainID
        if (subchainID.eq.branchmonomers+1) then
            call connect_to_monomer(scIDbranch,n)   ! Connect monomer to
                                                    ! new branch 
            call connect_to_monomer(subchainID,glob_n_branch) ! And vice-versa
            call connect_to_monomer(subchainID+1,n) ! Continue chain
        else if (subchainID .gt. branchmonomers+1 .and. &
                 subchainID .lt. nmonomers) then
            call connect_to_monomer(subchainID+1,n)
            call connect_to_monomer(subchainID-1,n)
        else    
            call connect_to_monomer(subchainID-1,n) ! End chain
        end if
    end do

    proc_chains(irank) = chainID
    proc_nps(irank)    = np
    call globalSum(proc_chains,nproc)
    call globalSum(proc_nps,nproc)
    
    do n=1,np
        if (monomer(n)%chainID.ne.0) then
            monomer(n)%chainID = monomer(n)%chainID + &
                                 sum(proc_chains(1:irank)) - proc_chains(irank)
        end if
        monomer(n)%glob_no     = monomer(n)%glob_no + &
                                 sum(proc_nps(1:irank))    - proc_nps(irank)
    end do

    nchains = sum(proc_chains)

end subroutine setup_initialise_polyinfo_singlebranched
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!Initialise Velocities




!----------------------------------------------------------------------------------
! Set velocity of a range of bins to prescribed value
! note that i,j and k bin values are in the global bin coordinate system
! which runs from 1 to gnbins. The local nbins include halos!

module set_bin_velocity_mod

contains


! 15/04/14 -- THERE APPEARS TO BE A BUG IN THIS ROUTINE WHICH OCCURS WHEN
!             CERTAIN CONFIGURATIONS OF CELLS PER BIN ARE USED. IT MAY BE 
!             WHEN NON-EVEN NUMBERS OF CELLS PER BIN ARE APPLIED... ALWAYS
!             CHECK THE CORRECT VELOCITY IS OBTAINED

subroutine set_bin_velocity(imin, imax, jmin, jmax, kmin, kmax, velocity,veltype)
	use linked_list, only : cell, node
	use arrays_MD, only : r,v,a
	use computational_constants_MD, only : binspercell, iblock, jblock, kblock, globaldomain, halfdomain, & 
										   npx, npy, npz, iter, irank, ncells, delta_t
	use calculated_properties_MD, only : gnbins, nbins
	implicit none

	integer, intent(in)                			:: veltype ! 0 = vbin/Nbin, 1 = vbin/binvolume
	integer, intent(in)                			:: imin, imax, jmin, jmax, kmin, kmax
	double precision,dimension(3),intent(in)	:: velocity   !Overall momentum of system

	integer			                			:: iminl, imaxl, jminl, jmaxl, kminl, kmaxl
	integer										:: ibinmin,jbinmin,kbinmin,ibinmax,jbinmax,kbinmax
	integer										:: i,icell,jcell,kcell,molno,binNsum,cellnp
	integer	,dimension(3)						:: p_lb, p_ub
	double precision,dimension(3)				:: binvsum, vcorrection,cellsperbin
	double precision,dimension(3)				:: r_temp,v_temp,binsize,binmin,binmax
	type(node), pointer 	        			:: old, current

	if (imin .ne. imax) stop "Error set_bin_velocity -- bin indices imin and imax currently must be the same"
	if (jmin .ne. jmax) stop "Error set_bin_velocity -- bin indices jmin and jmax currently must be the same"
	if (kmin .ne. kmax) stop "Error set_bin_velocity -- bin indices kmin and kmax currently must be the same"

	p_lb(1) = (iblock-1)*floor(gnbins(1)/real((npx),kind(0.d0)))
	p_ub(1) =  iblock *ceiling(gnbins(1)/real((npx),kind(0.d0)))
	p_lb(2) = (jblock-1)*floor(gnbins(2)/real((npy),kind(0.d0)))
	p_ub(2) =  jblock *ceiling(gnbins(2)/real((npy),kind(0.d0)))
	p_lb(3) = (kblock-1)*floor(gnbins(3)/real((npz),kind(0.d0)))
	p_ub(3) =  kblock *ceiling(gnbins(3)/real((npz),kind(0.d0)))

	!Convert to local bin number from input which is global bin number
	if (imin .gt. p_lb(1) .and. imax .le. p_ub(1)) then 
		iminl = imin - p_lb(1)+1
		imaxl = imax - p_lb(1)+1
	else
		return
	endif
	if (jmin .gt. p_lb(2) .and. jmax .le. p_ub(2)) then 
		jminl = jmin - p_lb(2)+1
		jmaxl = jmax - p_lb(2)+1
	else
		return
	endif
	if (kmin .gt. p_lb(3) .and. kmax .le. p_ub(3)) then 
		kminl = kmin - p_lb(3)+1
		kmaxl = kmax - p_lb(3)+1
	else
		return
	endif

	!Calculate bin to cell ratio
	cellsperbin = 1.d0/binspercell !ceiling(ncells(1)/dble(nbins(1)))
	where (cellsperbin .lt. 1.d0) cellsperbin = 1.d0
	!Safety check for non-integer cell ratios
	if (any(abs(ncells/nbins - dble(ncells)/dble(nbins)) .gt. 0.000000000001d0)) then
		stop "ERROR in set_bin_velocity -- Specified bin/cell ratio results in non-integer number of cells!"
	endif
	binsize = globaldomain/gnbins
	!Get bin extents -- minus one due to halo bins
	binmin(1) = (iminl-2) * binsize(1) - halfdomain(1)
	binmax(1) = (imaxl-1) * binsize(1) - halfdomain(1)
 	binmin(2) = (jminl-2) * binsize(2) - halfdomain(2)
	binmax(2) = (jmaxl-1) * binsize(2) - halfdomain(2)	
	binmin(3) = (kminl-2) * binsize(3) - halfdomain(3)
	binmax(3) = (kmaxl-1) * binsize(3) - halfdomain(3)

	!Get cell number from bin numbers
	ibinmin = (iminl-1)*cellsperbin(1)+1+(1-cellsperbin(1))
	ibinmax =  imaxl   *cellsperbin(1)  +(1-cellsperbin(1))
	jbinmin = (jminl-1)*cellsperbin(2)+1+(1-cellsperbin(2))
	jbinmax = jmaxl    *cellsperbin(2)  +(1-cellsperbin(2))
	kbinmin = (kminl-1)*cellsperbin(3)+1+(1-cellsperbin(3))
	kbinmax = kmaxl    *cellsperbin(3)  +(1-cellsperbin(3))

    !print'(18i4,6f9.4)', imin, imax, jmin, jmax, kmin, kmax,iminl, imaxl, jminl, jmaxl, kminl, kmaxl, ibinmin,jbinmin,kbinmin,ibinmax,jbinmax,kbinmax, binmin, binmax

	!Calculate velocity in bin
	binNsum = 0; binvsum = 0.d0
	do kcell=kbinmin, kbinmax
	do jcell=jbinmin, jbinmax 
	do icell=ibinmin, ibinmax 

		!print'(2i5,a,3i4,2(a,3i4),a,3f10.6,2(a,3i4),i5)', iter, iblock,' Cells =', icell,jcell,kcell,' Bins= ',imin,jmin,kmin,' Binsl= ',iminl,jminl,kminl,' cellperbin= ',cellsperbin, 'nbins =', nbins , ' gnbins =', gnbins
	
		cellnp = cell%cellnp(icell,jcell,kcell)
		old => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp					!Step through each particle in list 
			molno = old%molno 	 	!Number of molecule
			binNsum = binNsum + 1    
			binvsum(:) = binvsum(:) + v(:,molno)

			!GET VELOCITY AT NEXT TIMESTEP
			v_temp(:) = v(:,molno) + delta_t * a(:,molno) 	
			r_temp(:) = r(:,molno) + delta_t * v_temp(:) 

			!BIN VELOCITY
!			if (r_temp(1)  .lt. binmin(1) .or. & 
!			    r_temp(1)  .gt. binmax(1) .or. & 
!			    r(1,molno) .lt. binmin(1) .or. & 
!			    r(1,molno) .gt. binmax(1)) print'(a,i4,6(a,f9.4))', "set_bin_vel -- Mol Outside x bin ", iminl, & 
!															 " min ", binmin(1), &
!															 " r before = ", r(1,molno), " r after = ", r_temp(1), & 
!															 " max ", binmax(1), & 
!															 " v before = ", v(1,molno), " v after = ", v_temp(1)
!			if (r_temp(2)  .lt. binmin(2) .or. & 
!			    r_temp(2)  .gt. binmax(2) .or. & 
!			    r(2,molno) .lt. binmin(2) .or. & 
!			    r(2,molno) .gt. binmax(2)) print'(a,i4,6(a,f9.4))', "set_bin_vel -- Mol Outside y bin ", jminl, & 
!															 " min ", binmin(2), &
!															 " r before = ", r(2,molno), " r after = ", r_temp(2), & 
!															 " max ", binmax(2), & 
!															 " v before = ", v(2,molno), " v after = ", v_temp(2)
!			if (r_temp(3)  .lt. binmin(3) .or. & 
!			    r_temp(3)  .gt. binmax(3) .or. & 
!			    r(3,molno) .lt. binmin(3) .or. & 
!			    r(3,molno) .gt. binmax(3)) print'(a,i4,6(a,f9.4))', "set_bin_vel -- Mol Outside z bin ", kminl, & 
!															 " min ", binmin(3), &
!															 " r before = ", r(3,molno), " r after = ", r_temp(3), & 
!															 " max ", binmax(3), & 
!															 " v before = ", v(3,molno), " v after = ", v_temp(3)

			!print'(i5,a,7i6,6f10.5)',iter,' velocities ',i,cellnp,molno,binNsum,icell,jcell,kcell, r(:,molno), v(:,molno)

			current => old
			old => current%next !Use pointer in datatype to obtain next item in list
		enddo

	enddo
	enddo
	enddo

	!Calculate velocity correction per molecule
	if (binNsum .eq. 0) then
		print*, "No molecules in bin ", imin, jmin, kmin
		return
	endif

	if (veltype .eq. 0) then
		vcorrection(:) = binvsum(:)/binNsum - velocity(:)
		!print'(3(a,3f10.5))', 'applied vel ', velocity, ' bin v ', binvsum(:)/binNsum, ' v correct ', vcorrection
	elseif (veltype .eq. 1) then
		vcorrection(:) = binvsum(:) - velocity(:)*product(binsize)
		vcorrection = vcorrection/binNsum
		!print'(3(a,3f10.5))', 'applied mom ', velocity, ' bin v ', binvsum(:)/product(binsize), ' v correct ', vcorrection
	else
		stop "Error in set_bin_velocity -- velocity type not correctly specifiy (must be 0 or 1)"
	endif



	!Apply velocity correction per molecule to bin
	binNsum = 0; binvsum = 0.d0
	do kcell=kbinmin, kbinmax
	do jcell=jbinmin, jbinmax 
	do icell=ibinmin, ibinmax 
	
		cellnp = cell%cellnp(icell,jcell,kcell)
		old => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list
!        print*, icell,jcell,kcell,cell%cellnp(icell,jcell,kcell)
!        if (jmax .eq. nint(gnbins(2)/2.d0) + 1) then
!            print'(i5,a,7i6)',iter, ' corrected_vel ',i,cellnp,molno,binNsum,icell,jcell,kcell
!        endif

		do i = 1,cellnp					!Step through each particle in list 
			molno = old%molno 	 	!Number of molecule
			v(:,molno) =  v(:,molno) - vcorrection

            !print'(i5,a,7i6,6f10.5)',iter, ' corrected_vel ',i,cellnp,molno,binNsum,icell,jcell,kcell,r(:,molno),v(:,molno)

			binNsum = binNsum + 1    
			binvsum(:) = binvsum(:) + v(:,molno)

			current => old
			old => current%next !Use pointer in datatype to obtain next item in list
		enddo

	enddo
	enddo
	enddo

 	if (veltype .eq. 0) then
        if (any(abs(binvsum(:)/binNsum-velocity) .gt. 0.000000001d0)) then
         	print'(i8,a,3f20.16,a,3i4,a,3f10.5)', jblock, ' Corrected velocity is then ',  binvsum(:)/binNsum, & 
 	    										 ' in Bin= ',imin,jmin,kmin, ' should be ', velocity
        endif
 	elseif (veltype .eq. 1) then
        if (any(abs(binvsum(:)/product(binsize)-velocity) .gt. 0.000000001d0)) then
     	    print'(i8,a,3f20.16,a,3i4,a,3f10.5)', jblock, ' Corrected momentum : ',  binvsum(:)/product(binsize), & 
 											 ' in Bin= ',imin,jmin,kmin, ' should be ', velocity
        endif
 	else
 		stop "Error in set_bin_velocity -- velocity type not correctly specifiy (must be 0 or 1)"
 	endif

end subroutine set_bin_velocity


end module set_bin_velocity_mod






! Set up the intial velocities of the particles using velocity magnitude calculated
! from intitial temperature and giving random vectorial components

subroutine setup_initialise_velocities
    use module_initialise_microstate
    use messenger_data_exchange, only : globalSum
    implicit none

    integer                            :: n,i 
    double precision, dimension (nd)   :: netv   !Overall momentum of system

    !If ensemble is not tag based, set all molecules to unfixed
    if (ensemble .ne. tag_move) then
        allocate(fix(3,np)); fix = 1
    endif

    !Use definition of temperature and re-arrange to define an average velocity
    initialvel = sqrt(nd * (1.d0 - 1.d0/globalnp)*inputtemperature)

    v = 0.d0    !Set velocity initially to zero
    i = 0       !Zero number of molecules with velocity assigned
    netv=0.d0   !Set net velocity of system to zero initially

    do n=1,np                               !Step through each molecule
        if (fix(1,n) .eq. 1) then           !For x component as fix occurs in all 3 dimensions
            call random_number(rand)        !Generate a random number for each dimension/particle
            angle  = 2.d0*pi*rand           !Random angle between 0 and 2pi
            v(2,n) = initialvel*sin(angle)  !Y component of velocity magnitude for random angle
            v13    = initialvel*cos(angle)  !Magnitude of x and z component
            call random_number(rand)        !Generate a new random number
            angle  = 2.d0*pi*rand           !Random angle between 0 and 2pi     
            v(1,n) = v13*sin(angle)         !X component of velocity magnitude for random angle
            v(3,n) = v13*cos(angle)         !Z component of velocity magnitude for random angle
            i = i + 1                       !Count number of molecules with velocity assigned
        else
            v(:,n) = 0.d0                   !Don't assign velocity if molecule is fixed
        endif
        netv(:)= netv(:) + v(:,n)           !Sum up overall momentum of system due to random movement
    enddo

    call globalSum(netv, nd)            !Sum net velocity on all processors
    call globalSum(i)                    !Sum number of molecules assigned velocity on all processors

    if(i .ne. 0) netv(:) = netv(:)/i        !Divide overall momentum by number of particles

    do n=1,np
        !reducing all non-fixed particles by same amount
        if (fix(1,n) .eq. 1) v(:,n)= v(:,n)-netv(:) 
                   
    enddo

    !If ensemble is not tag based, set all molecules to unfixed
    if (ensemble .ne. tag_move) then
        deallocate(fix)
    endif

end subroutine setup_initialise_velocities

!Initialise Velocities
! Set up the intial velocities of the particles using Taylor Green velocity
! contours

subroutine setup_initialise_velocities_TG
    use module_initialise_microstate
    use messenger_data_exchange, only : globalSum
    implicit none

    integer                            :: n,i 
    double precision                   :: x,y,z,Lx,Ly,Lz
    double precision, dimension (nd)   :: netv   !Overall momentum of system

    !Use definition of temperature and re-arrange to define an average velocity
    initialvel = sqrt(nd * (1.d0 - 1.d0/globalnp)*inputtemperature)

    v = 0.d0    !Set velocity initially to zero
    i = 0       !Zero number of molecules with velocity assigned
    netv=0.d0   !Set net velocity of system to zero initially

    do n=1,np                               !Step through each molecule
        x  = r(1,n);    y  = r(2,n);    z  = r(3,n);
        Lx = halfdomain(1); Ly = halfdomain(2); Lz = halfdomain(3); !Domain should be cubic...
        v(1,n) =  initialvel*sin(pi*x/Lx)*cos(pi*y/Ly)*cos(pi*z/Lz)
        v(2,n) = -initialvel*cos(pi*x/Lx)*sin(pi*y/Ly)*cos(pi*z/Lz)
        v(3,n) =  initialvel*cos(pi*x/Lx)*cos(pi*y/Ly)*sin(pi*z/Lz)
        netv(:)= netv(:) + v(:,n)           !Sum up overall momentum of system due to random movement
    enddo

    call globalSum(netv, nd)            !Sum net velocity on all processors
    netv(:) = netv(:)/np        !Divide overall momentum by number of particles

    do n=1,np
        !reducing all particles by same amount
        v(:,n)= v(:,n) - netv(:) 
                   
    enddo

end subroutine setup_initialise_velocities_TG


subroutine setup_initialise_velocities_TG_parallel
    use module_initialise_microstate
    use messenger_data_exchange, only : globalSum
    implicit none

    integer                            :: n,i 
    double precision                   :: x,y,z,Lx,Ly,Lz
    double precision, dimension (nd)   :: netv   !Overall momentum of system

    !Use definition of temperature and re-arrange to define an average velocity
    initialvel = sqrt(nd * (1.d0 - 1.d0/globalnp)*inputtemperature)

    v = 0.d0    !Set velocity initially to zero
    i = 0       !Zero number of molecules with velocity assigned
    netv=0.d0   !Set net velocity of system to zero initially

    do n=1,np                               !Step through each molecule
        x  = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
        y  = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
        z  = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
        Lx = 0.5d0*globaldomain(1); Ly = 0.5d0*globaldomain(2); Lz = 0.5d0*globaldomain(3); !Domain should be cubic...
        v(1,n) =  initialvel*sin(pi*x/Lx)*cos(pi*y/Ly)*cos(pi*z/Lz)
        v(2,n) = -initialvel*cos(pi*x/Lx)*sin(pi*y/Ly)*cos(pi*z/Lz)
        v(3,n) =  initialvel*cos(pi*x/Lx)*cos(pi*y/Ly)*sin(pi*z/Lz)
        netv(:)= netv(:) + v(:,n)           !Sum up overall momentum of system due to random movement
    enddo

    print*, 'before sum', irank, netv, nd

    call globalSum(netv, nd)            !Sum net velocity on all processors
    netv(:) = netv(:)/np        !Divide overall momentum by number of particles

    print*, 'after sum', irank, netv, nd

    do n=1,np
        !reducing all particles by same amount
        v(:,n)= v(:,n) - netv(:) 
                   
    enddo

end subroutine setup_initialise_velocities_TG_parallel

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!Initialise Velocities TEST
! Set up the intial velocities to test values for debugging

subroutine setup_initialise_velocities_test
    use module_initialise_microstate
    implicit none

    !integer                            :: i, n
    !double precision, dimension (nd)   :: netv   !Overall momentum of system

    !Use definition of temperature and re-arrange to define an average velocity
    !initialvel = sqrt(nd * (1.d0 - 1.d0/np)*inputtemperature)
    v = 0.d0    !Set velocity initially to zero
    !i = 0      !Zero number of molecules with velocity assigned
    !netv=0     !Set net velocity of system to zero initially
    !zeta=0.d0  !Set Nose Hoover thermostat scaling property to zero
    
    !v(1,1) = 2.d0
	!Fextbinsize = domain/nbins

    !do n=1,np                      !Step through each molecule
        !r(1,:) = halfdomain(:)
        !v(1,100) = -1.0d0
        !v(2,100) = -0.0d0
        !v(3,100) = -0.0d0

        v(1,99) = +1.0d0
        v(2,99) = +0.0d0
        v(3,99) = +0.0d0
        
        !v(:,99:100) = 2.d0


        !r(1,:) = -halfdomain(:)
        !v(1,n) = -0.0d0 
        !v(2,n) = -0.0d0
        !v(3,n) = -0.0d0
    
    !enddo
    
!   v(1,:) = 0.5d0
    
!   v(7,3) = -0.5d0
!   v(4,3) = 0.5d0

end subroutine setup_initialise_velocities_test




subroutine set_velocity_field_from_couette_analytical(t,Re,Uwall,H,slidewall,ixyz)
    use interfaces, only : error_abort
    use calculated_properties_MD, only : gnbins
	use computational_constants_MD, only : irank
	use librarymod, only : couette_analytical_fn
	use set_bin_velocity_mod
    implicit none

	integer, intent(in)				:: ixyz, slidewall
	double precision, intent(in)	:: t, Re, Uwall, H


	integer										:: ibin, jbin, kbin,appliedbins
	double precision,dimension(3)				:: binvel
	double precision,dimension(:),allocatable	:: utemp

	appliedbins = gnbins(ixyz)
	allocate(utemp(gnbins(ixyz))); utemp = 0.d0
	utemp = couette_analytical_fn(t,Re,Uwall,H,appliedbins,slidewall)

	!Create cell lists to be used in specifying velocity!
    call assign_to_cell

    !Set MD velocity
    do jbin=2,gnbins(2)+1
    do ibin=2,gnbins(1)+1
    do kbin=2,gnbins(3)+1
        binvel(1) =  utemp(jbin-1)
        binvel(2) =  0.d0
        binvel(3) =  0.d0
        call set_bin_velocity(ibin-1, ibin-1, & 
							  jbin-1, jbin-1, & 
							  kbin-1, kbin-1, binvel,0)
    enddo
    enddo
    enddo


end subroutine set_velocity_field_from_couette_analytical




subroutine set_velocity_field_from_DNS_restart(filename,ngx,ngy,ngz)
    use interfaces, only : error_abort
    use calculated_properties_MD, only : gnbins
    use librarymod, only : read_DNS_velocity_files
	use set_bin_velocity_mod
    implicit none

    integer,intent(in)      :: ngx, ngy, ngz
    character(*),intent(in):: filename

    integer     :: i,j,k,ibin,jbin,kbin
    logical, dimension(3)       :: bin_error
    double precision, dimension(3)  :: binvel
    double precision, dimension(:,:,:),allocatable  :: uc, vc, wc

    !Read DNS data into arrays
    call read_DNS_velocity_files(trim(filename),ngx,ngy,ngz,uc,vc,wc)

    !Check DNS size vs number of bins
    bin_error = .false.
    if (gnbins(1) .ne. ngx-1) then
        bin_error(1) = .true.
        print'(2(a,i8))', ' gnbinsx = ',gnbins(1),' DNS restart file bins = ', ngx-1
    endif
    if (gnbins(2) .ne. ngy-1) then
        bin_error(2) = .true.
        print'(2(a,i8))', ' gnbinsy = ',gnbins(2),' DNS restart file bins = ', ngy-1
    endif
    if (gnbins(3) .ne. ngz-1) then
        bin_error(3) = .true.
        print'(2(a,i8))', ' gnbinsz = ',gnbins(3),' DNS restart file bins = ', ngz-1
    endif
    if (any(bin_error)) then
        print'(3(a,l))', ' Error in x ', bin_error(1),' Error in y ', bin_error(2),' Error in z ', bin_error(3)
        call error_abort("Error -- number of bins disagrees with DNS initial velocity file")
    endif

	!Create cell lists to be used in specifying velocity!
    call assign_to_cell

    !Set MD velocity
    do k=2,size(wc,1)-2
    do i=2,size(uc,2)-2
    do j=2,size(vc,3)-2
        binvel(1) =  0.5d0*(uc(k,i,j)+ uc(k,i+1,j))
        binvel(2) =  0.5d0*(vc(k,i,j)+ vc(k,i,j+1))
        binvel(3) =  0.5d0*(wc(k,i,j)+ wc(k+1,i,j))
        ibin = i - 1; jbin = j - 1; kbin = k - 1
        !print'(a,3i8,4f10.5)', 'BIN NUMBER & VEL = ', i,j,k,binvel,uc(k,i,j)
        call set_bin_velocity(ibin, ibin, jbin, jbin, kbin, kbin, binvel,0)
    enddo
    enddo
    enddo

    call linklist_deallocate_cells            !Deallocate all linklist components
    deallocate(uc,vc,wc)

end subroutine set_velocity_field_from_DNS_restart

