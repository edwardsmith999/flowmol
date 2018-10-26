!=====================================================================================
!------------------------------DEFINE MODULES-----------------------------------------
!=====================================================================================

!-----------------------------Physical Constants---------------------------------
!Constants which represent physical values and consequently change the physical
!conditions of the simulation
module physical_constants_MD

	integer, parameter 				:: nd = 3		   		!Number of dimensions
	integer            				:: np                   !Number of particles
	integer		   					:: globalnp             !Global number of particles
	integer                         :: tethernp             !Number of tethered particles
	integer		   					:: halo_np              !Number of molecules in halo
    integer                         :: reinsertnp           !Number of molecules to reinsert (leaving over surfaces)
	integer,dimension(:),allocatable:: procnp 				!Array of all processors np
	integer,dimension(:),allocatable:: proctethernp 		!Array of all processors np
	double precision   				:: volume, density      !Define constant volume and density
	double precision   				:: rcutoff, halfrcutoff !Cut off distance for particle interactions
	double precision   				:: rcutoff2             !Cut off distance for particle interactions squared
	double precision   				:: potential_sLRC 		!Long range potential correction 
	double precision   				:: pressure_sLRC 		!Long range pressure correction 
	double precision   				:: inputtemperature     !Define initial temperature
	double precision   				:: initialunitcell      !Initial size of unit cell
	double precision   				:: initialvel           !Initial velocity of particles
	double precision,parameter 		:: pi=4.d0*atan(1.d0)
	double precision, dimension(:),allocatable :: thermostattemperature!Define thermostat setpoint temperature

	!Simulation setup conditions
	real(kind(0.d0)), dimension(3)	:: fixdisttop, slidedisttop, fixdistbottom, slidedistbottom, wallslidev
	real(kind(0.d0)), dimension(3)	:: tethereddisttop, tethereddistbottom, thermstattop,thermstatbottom
	real(kind(0.d0)), dimension(3)	:: emptydistbottom, emptydisttop

end module physical_constants_MD

!------------------------------------------------------------------------------
!-----------------------------Computational Constants--------------------------
! Constants used to determine the computational parameters of the simulation
module computational_constants_MD

	!VOID value for data initialisation
 	integer, parameter :: VOID=-666			

    integer :: finalstate_version_no = huge(0)

	!Command-line arguments
	logical					:: restart
	character(len=200) 		:: input_file, initial_microstate_file

	!Force and Potential flags
	integer						  :: force_list	   		 !flag for neighbr/cell list
	integer						  :: potential_flag  	 !Choose LJ or Polymer potential
    integer                       :: Mie_potential       !Choose if LJ/Polymer replaced by more general Mie potential
    integer                       :: default_moltype=1   !Default moltype for mie potential -- taken to be LJ (=1) here
	integer                 	  :: tether_flag     	 !True if there exists 
	integer                 	  :: external_force_flag !Apply external forces?
	integer                 	  :: F_ext_ixyz			 !Direction of external forces
	real(kind(0.d0))        	  :: F_ext				 !Magnitude of external forces
    real(kind(0.d0)),dimension(2) :: eij_wall            !Interaction potential strength for substrate wall
	real(kind(0.d0)),dimension(6) :: F_ext_limits		 !Limits of region external forces applied to


    !Rebuild check flags
    integer             :: rebuild_criteria    !Choice of rebuild criteria
    integer             :: fixed_rebuild      !Fixed rebuild frequency
	real(kind(0.d0)) 	:: rneighbr2         !Square of rcuttoff+delta_rneighbr
	real(kind(0.d0)),dimension(3) 	:: delta_rneighbr    !Radius used for neighbour list construction

	! Move particle tags
	logical 		   :: tag_thermostat_active
	logical 		   :: dynamically_update_tags = .false. !Tags are always reset based on position each timestep
	integer, parameter :: free = 0
	integer, parameter :: fixed = 1
	integer, parameter :: fixed_slide = 2
	integer, parameter :: teth = 3
	integer, parameter :: thermo = 4
	integer, parameter :: teth_thermo = 5
	integer, parameter :: teth_slide = 6
	integer, parameter :: teth_thermo_slide = 7
	integer, parameter :: PUT_thermo = 8
	integer, parameter :: z_thermo = 9
	integer, parameter :: cyl_teth_thermo_rotate = 10
	integer, dimension(5), parameter :: tether_tags=(/teth,teth_thermo,teth_slide,teth_thermo_slide,cyl_teth_thermo_rotate/)
	integer, dimension(6), parameter :: thermo_tags=(/thermo,teth_thermo,teth_thermo_slide,PUT_thermo,z_thermo,cyl_teth_thermo_rotate/)
	integer, dimension(2), parameter :: fixed_tags=(/fixed, fixed_slide/)
	integer, dimension(3), parameter :: slide_tags=(/teth_slide,teth_thermo_slide, fixed_slide/)

	! Wall texture flags
	integer			   :: texture_type, texture_therm
    integer            :: nthermo
	integer, parameter :: posts = 1
	integer, parameter :: roughness = 2
	integer, parameter :: converge_diverge = 3
	integer, parameter :: fractal = 4
	real(kind(0.d0))   :: texture_intensity
	real(kind(0.d0))   :: local_heat_region(6)
	

	!Integration algorithm
	integer                 :: integration_algorithm
	integer                 :: ensemble
	integer, parameter      :: & 
		leap_frog_verlet = 0, &
		velocity_verlet  = 1, &
		nve         = 0, &
		nvt_NH      = 1, &
		nvt_GIK     = 2, &
		nvt_PUT_NH  = 3, &
		nvt_pwa_NH  = 4, &
		nvt_DPD     = 5, &
		tag_move    = 6, &
		SLLOD       = 7

    !Misc flags
    integer :: global_numbering         !Enforce global molecular numbering
    integer :: reset_tags_on_restart !Reset tags on restart

	!Initial configuration selection
	integer           	:: initial_config_flag
	integer         	:: lg_direction     !Direction in which the domain is split into liquid and gas
	character(len=128)	:: config_special_case
	real(kind(0.d0))	:: liquid_density	!Density of liquid if solid/liquid case used
	real(kind(0.d0))	:: gas_density	    !Density of liquid if gas/liquid case used
	real(kind(0.d0))	:: lg_fract	        !Fraction of the domain which is liquid (0 = all gas, 1 = all liquid)
	real(kind(0.d0))	:: dropletH =0.d0,dropletHLratio=0.d0   !Droplet height and H to length ratio
	real(kind(0.d0))	:: rbubble, rcentre(3)   !Radius of bubble
    logical             :: Twophase_from_file = .false.
	character(len=128)	:: FEA_filename

	!Initial velocity selection
	integer           	:: initial_velocity_flag
	character(len=128)	:: velocity_special_case
	character(len=128)	:: DNS_filename
	integer           	:: DNS_ngx,DNS_ngy,DNS_ngz


	!Write a separate file for each timestep
	logical ::	separate_outfiles = .false.
    !Start file numbering at beginning again?
    logical :: restart_numbering = .true.

	!Input (on or off) flags
	integer	:: & 
		vmd_outflag, vmd_skip, vmd_start, vmd_finish, vmd_count=1, &
		macro_outflag, &
		sLRC_flag,	&
		mass_outflag, &
		momentum_outflag, &
		temperature_outflag, &
		energy_outflag, &
		pressure_outflag, &
		heatflux_outflag, &
		viscosity_outflag, &
		centre_of_mass_outflag, &
		rdf_outflag, &
		rtrue_flag, &
		prev_rtrue_flag, &
		ssf_outflag, &
		vPDF_flag, & 
		cv_conserve,	&
		mflux_outflag, &
		vflux_outflag, &
		eflux_outflag, &
        msurf_outflag, &
		proc_reorder, &				    !Reorder processors at restart
		pass_vhalo = 0, &
		peculiar_flag, &	 			!Take streaming velocity away from temperature 	
		CVforce_flag = VOID, & 			!Type of CV force to apply
		CVweighting_flag = 0, &			!Distribution of CV forces
		CVforce_correct = 0, &			!Apply force to correct CV velocity to setpoint
        CVforce_correct_nsteps, &       !Number of steps to apply velocity correction over
		CVforce_starttime, &			!Start time of applied force
        cluster_analysis_outflag        !Cluster analysis with interface tracking

    logical,dimension(3)    :: CVforce_applied_dir  !Which directions to apply force

	!Add debugging CV flags
	integer,dimension(6) :: F_CV_limits		 !Limits of CV forces
	integer	             :: CV_debug=0
    integer,dimension(3) :: debug_CV

	integer, dimension(3)	:: periodic

	!Tethered molecules spring potential coefficients
	real(kind(0.d0))    :: teth_k2,teth_k4,teth_k6

	!Parameters
	integer	:: & 
		nproc, npx, npy, npz, 	&	!Number of MPI ranks and cartesian topology sizes
		iter, 					&	!Global simulation iteration count
		tplot, 					&	!Frequency at which to record results
		teval, 					&	!Frequency at which to evaluate results
		Nmass_ave, 				&	!Number of averages for each mass average
		Nvel_ave, 				&	!Number of averages for each velocity average
		NTemp_ave, 				&	!Number of averages for each temperature measurement
		Nenergy_ave,			&	!Number of averages for each energy measurement
		Ncom_ave,			    &	!Number of averages for each centre of mass measurement
		Nstress_ave, 			&	!Number of averages for VA or virial stress calculation
		Nheatflux_ave, 			&	!Number of averages for VA heat flux calculation
		split_kin_config, 		&	!Flag to determine if kinetic and configurational stress separated
		split_hfkin_config, 	&	!Flag to determine if kinetic and configurational heatflux are separated
		split_pol_sol_stats,    &	!Flag to determine if polymer and solvent stats separated
		Nvisc_ave, 				&	!Number of samples for viscosity measurement
		Nmflux_ave, 			&	!Number of averages for each mass flux
		Nvflux_ave, 			&	!Number of averages for each velocity flux
		Neflux_ave, 			&	!Number of averages for each energy flux
		Nsurfm_ave, 			&	!Number of averages for each surface mass
		initialstep, 			&	!Initial step of simulation
		finalstep,              &   !Final step of simulation
		Nsteps, 				&	!Total number of computational steps
		initialise_steps, 		&	!Initialisation steps to run before simulation start
		extralloc, 				&	!Extra allocation space to include copied halos
		overlap, 				&	!Size of overlap region used to apply force to molecular region
		Nvmd_intervals,         &
		rdf_nbins,              &   !Number of discrete "bins" used to calculate g(r)
		ssf_ax1,                &   !1st projection axis for static structure factor
		ssf_ax2,                &
		ssf_nmax                    !Maximum wavenumber for S(k) calculation

	integer,dimension(:,:),allocatable	:: vmd_intervals			!Multiple intervals for vmd record

	real(kind(0.d0)) 	:: delta_t           !Size of timestep for each computational step
	real(kind(0.d0)) 	:: elapsedtime       !Total time elapsed including restarts
	real(kind(0.d0)) 	:: simtime=0.d0      !Total incremented simulation time
	real(kind(0.d0))    :: rdf_rmax          !Maximum radius for radial distribution function
	real(kind(0.d0))	:: rescue_snapshot_freq	!Rescue snapshot output frequency in seconds

    !Constants for probability density function
	integer             :: NvPDF_ave   !Number of averages for each velocity PDF 
	integer             :: NPDFbins    !Number of histogram bins for velocity PDF 
	real(kind(0.d0))	:: PDFvlims    !Velocity Probability density functions min/max value


	real(kind(0.d0)), dimension(3)			:: binspercell     !Number of avergaing bins per computational cell
	
	!Store surface bins of processes subdomain for outputs over periodic boundaries
	integer									:: nsurfacebins     !Number of surface bins
	integer,allocatable,dimension(:,:)	    :: surfacebins		!Surface Bins

	integer									:: nsurfacecells     !Number of surface bins
	integer,allocatable,dimension(:,:)	    :: surfacecells		!Surface Bins

	integer		 								:: nhalobins 	!Number of halo bins
	integer,allocatable,dimension(:,:),target	:: halobins		!halo Bins

	integer		 								:: nhalocells 	!Number of halo bins
	integer,allocatable,dimension(:,:),target	:: halocells		!halo Bins

	integer		 								:: nhalocellbins 	!Minimum of halo bins or cells
	integer,allocatable,dimension(:,:),target	:: halocellbins		!Minimum values of halo bins or cells

	!Number and size of unit used for initial setup of molecules (i.e. FCC unit)
	integer,          dimension(3)		:: initialnunits
	real(kind(0.d0)), dimension(3) 		:: initialunitsize

	!Size of global computational domain and domain per processor
	real(kind(0.d0)), dimension(3)	:: globaldomain
	real(kind(0.d0)), dimension(3)	:: domain, halfdomain

	!Number and size of cells used for domain subdivision into cells of size ~rcutoff
	integer,          dimension(3)	:: ncells
	real(kind(0.d0)), dimension(3)	:: cellsidelength, halfcellsidelength

	!Array Storing mapping from cell number to point on 3D
	!Hilbert curve
	integer									:: sort_flag,sort_freq,sortblocksize
	integer,allocatable,dimension(:,:,:)    :: Hcurve		

	!Setup seed for random number generation
	integer,          dimension(:), allocatable :: seed 

	! Block/Process ID
	integer irank, iroot, ierr
	integer irankx, iranky, irankz
	integer iblock, jblock, kblock

 	!Set number of halos
	integer, parameter :: nh = 1
	integer, dimension(3) :: nhb


	!Number of cells per processor
	integer ncellxl, ncellyl, ncellzl    !Inc Halos
	integer nicellxl, nicellyl, nicellzl !Inner only

	!Directory that holds input/output files, useful in coupling mode
 	character(len=128) :: prefix_dir = "./"	

	!Calcultion, 0 -- Harasima contour (half per bin), 1 -- Line length per bin trapizium rule
	! and 2 -- Line length per bin explicit calculation 
	integer	:: VA_calcmethod, VA_heatflux_calcmethod

	!Number of samples used to calculate l_ij in VA stress calculation
	integer :: VA_line_samples, VA_heatflux_line_samples

end module computational_constants_MD

!-------------------------------------------------------------------------------------
!-------------------------------Shearing BCs------------------------------------------
module shear_info_MD

	integer 			:: define_shear_as
	integer				:: wrap_integer
	integer				:: le_i0                    !Lees-Edwards iter0
	integer				:: le_sp                    !Shear plane
	integer				:: le_sd                    !Shear direction
	integer				:: le_rp                    !Shear remaining plane
	real(kind(0.d0)) 	:: le_sv                    !Shear velocity
	real(kind(0.d0)) 	:: le_sr                    !Shear rate
	real(kind(0.d0))	:: le_sx                    !Shear distance
	real(kind(0.d0))	:: le_st                    !Shear time

	integer, dimension(:), allocatable	:: mol_wrap_integer		

end module shear_info_MD

!------------------------------------------------------------------------------
!----------------------------------Arrays--------------------------------------

module arrays_MD

	integer,          dimension(:),   allocatable, target	:: & 
        tag,     &   !Mol tags
        moltype, &   !Type used for interactions
        glob_no      !Global molecular number

	integer, 	  	  dimension(:,:), allocatable, target	:: fix  !Fixed molecules
	real(kind(0.d0)), dimension(:),   allocatable, target 	:: &
		potenergymol, 		&		!Potential energy of each molecule
		potenergymol_LJ, 	&		!LJ Potential energy of each molecule
		potenergymol_POLY,	&		!FENE Potential energy of each molecule
		potenergymol_mdt,	&		!Potential energy of each molecule at previous timestep
		virialmol,			&		!Virial of each molecule
		recvbuffer
	real(kind(0.d0)), dimension(:,:),   allocatable			:: &
		rtrue, 		&      			!Positions with no period BC
		vtrue,      &               !Corresponding velocities
		rtether, 	&
		rijsum, 	&				!Sum of all molecular rij values
		theta, 		&
		aD,aR
	real(kind(0.d0)), dimension(:,:),   allocatable, target 	:: &
		r, 		&        		  	!Positions
		v, 		&        		  	!Velocity
		a, 		&					!Accelerations
		a_old,	&					!Accelerations at last timestep
		slidev, &					!Speed for sliding molecules
		U							!Streaming velocity of molecules

end module arrays_MD


!------------------------------------------------------------------------------
!-------------------------------Taylor-Couette---------------------------------
module concentric_cylinders 

	real(kind(0.d0)) :: r_oo,r_oi,r_io,r_ii  !Outer and inner radii of cylinders
	                                         !r_ab - a=face, b=cyl

	integer :: cyl_units_oo                  !FCC units outer face of outer cyl
	integer :: cyl_units_io
	integer :: cyl_units_oi
	integer :: cyl_units_ii
	integer :: cyl_units_z

	integer :: cyl_np

	integer :: cyltag
	integer, parameter :: cyl_outer = 2, cyl_inner = 1, cyl_off = 0
	
	real(kind(0.d0)) :: cyl_density          !Density of cylinder walls

	real(kind(0.d0)) :: omega                !Actual angular vel
	real(kind(0.d0)) :: omega_i, omega_f     !Initial and final angular vel
	real(kind(0.d0)) :: omega_ramplength     !Time over which to ramp omega
	integer :: omega_rampiters               !Iterations over which to ramp omega

	character(len=200) :: cyl_file           !Output/input file
	
	integer :: cpol_bins(3)                  !Number of averaging bins in cpol
	integer :: cpol_binso(3)                 !Including halos
	integer :: gcpol_bins(3)                 !Glob averaging bins in cpol
	integer :: cpol_nhbz                     !Number of halo bins in z direction

	real(kind(0.d0)), dimension(:,:,:),   allocatable :: cyl_mass
	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: cyl_mom
	real(kind(0.d0)), dimension(:,:,:), allocatable :: cyl_KE

end module concentric_cylinders 

!-------------------------------------------------------------------------------------
!----------------------------------Polymer--------------------------------------------
module polymer_info_MD

	integer				:: etevtcf_outflag			!End-to-end time correlation function output flag
	integer				:: etevtcf_iter0			!Iteration from which to begin averaging
	integer             :: nchains                  !Number of FENE chains in domain
	integer				:: nmonomers				!Number of LJ beads per FENE chain
    integer             :: angular_potential        !Calculate and apply angular potential
	integer             :: r_gyration_outflag       !Radius of gyration outflag
	integer             :: r_gyration_iter0         !Iteration at which to start recording R_g

	integer             :: solvent_flag             !Solvent on/off flag
    real(kind(0.d0))    :: targetconc               !Solvent target concentration
    real(kind(0.d0))    :: surface_surfactant_layer !Solvent thickness at liquid surface
	real(kind(0.d0))    :: eps_pp, eps_ps, eps_ss   !Solvent parameters
	real(kind(0.d0)), parameter :: sod_cut =1.5d0   !Soddemann potential cutoff
	real(kind(0.d0)), parameter :: sod_cut2=2.25d0  ! 
	real(kind(0.d0)), parameter :: sod_a    = 3.1730728678
	real(kind(0.d0)), parameter :: sod_b    = -0.85622864544
	real(kind(0.d0)), parameter :: wca_cut  = 1.12246204830937
!	real(kind(0.d0)), parameter :: wca_cut2 = 1.25992104989486
	real(kind(0.d0)), parameter :: wca_cut2 = 1.25992104989487

	real(kind(0.d0)) 	:: k_c, R_0					!Spring constant and max elongation of bonds	
	real(kind(0.d0)) 	:: etevtcf                  !End-to-end vector time correlation function
	real(kind(0.d0))    :: R_g                      !Radius of gyration
	
	integer, parameter   :: max_funcy=4             !Maximum functionality of monomers in system
	integer, allocatable :: bond(:,:)               !Chain neighbour molnos
	integer, allocatable :: bondcount(:)            !Transient bond counter for linklist rebuild
	
	type monomer_info
		SEQUENCE									!For MPI convenience
		integer :: chainID                          !Integer value: chain number 
		integer :: subchainID                       !Integer value: bead number
		integer :: funcy                            !Functionality of the monomer
		integer :: glob_no                          !Global molecule number
		integer :: bin_bflag(4)                     !Integer for bit manipulation to find bond flags...
                                                    !...for more info see function get_bondflag
		! THE TOTAL NUMBER OF ITEMS IN THIS DATA TYPE MUST ALSO BE STORED IN THE VARIABLE nsdmi
	end type monomer_info

	type(monomer_info), dimension(:), allocatable,target :: monomer
	!eg. to access chainID of mol 23, call monomer(23)%chainID

	integer, parameter :: nsdmi=8                   !Number of sent doubles for monomer_info 

	type chain_info
		SEQUENCE									!For MPI convenience
		integer :: chainID                          !Integer value: chain number 
        integer :: nmonomers                        !Number of monomers in chain
        integer,allocatable,dimension(:) :: beads   !Molno of all beads in chain
	end type chain_info

	type(chain_info), dimension(:), allocatable :: chain

	integer, parameter :: intbits=bit_size(1)-1     !Size in bits of integer on this machine
    !-1, because 2**31 + 2**(something<31) not possible to represent with 32 bit integer
    !See get_bondflag explanation for a better idea.
    
    real(kind(0.d0)) :: grafting_density                                                            



	real(kind(0.d0)), dimension(:,:), allocatable :: etev
    !End-to-end vectors for polymer chains at iter=etevtcf_iter0 (size of np,only stored for leftmost chain mols) 
	real(kind(0.d0)), dimension(:,:), allocatable :: etev_0 

contains

	!----------------------------------------------------------------------
	! function get_bondflag(molno,scID)
	! author: David Trevelyan
	! 
	! get_bondflag returns:
	! 		1: if molecule MOLNO is connected to the monomer with
	!		   subchainID SCID
	!		0: otherwise
	!
	! get_bondflag is intended for use AFTER establishing that, for
	! example, two monomers are on the same chain:
	!
	!       if (monomer(i)%chainID .eq. monomer(j)%chainID) then
	!			bflag = get_bondflag(i,monomer(j)%subchainID))
	!			select case(bflag)
	!			case(0)
	!				THEY ARE NOT CONNECTED
	!			case(1)
	!				THEY ARE CONNECTED
	!			end select
	!		end if
	!
	! monomer(n)%bin_bflag is an array of 4 integers that, in binary
	! representation, represent a longer array of 0s and 1s. For example,
	! working with 4-bit integers to keep things simple, a bin_bflag
	!
	!		5 0 0 0 
	!
	! would represent an array that looks like:
	!
	!	0 0 0 0   0 0 0 0   0 0 0 0   0 1 0 1
	!      4         3         2         1
	!
	! ...where we count from right to left (from LSB to MSB).
	! 
	! Similarly, a bin_bflag
	!
	!		0 0 10 0
	!
	! represents
	!
	!	0 0 0 0   1 0 1 0   0 0 0 0   0 0 0 0
	!      4         3         2         1
	!
	! ... where a 1 in position scID signifies a connection to
	! the monomer in the same chain as monomer(n) with subchainID scID. 
	! 
	! Again, by example, to explain the meanings of "group" and "expo" we
	! assume 4-bit integer size. Given a subchainID of 9, we would seek
	! the bond flag at position 9 in the array of 16:
	!
	!  _  _  _  _   _  _  _ ?   _ _ _ _   _ _ _ _
	! 16 15 14 13  12 11 10 9   8 7 6 5   4 3 2 1
	!
	! So, we find the integer in bin_bflag(1:4) that stores the information
	! about subchainID 9 and call it "group". In this case:
	! 
	!	group = ceiling(real(9)/real(4)) = 3
	!
	! We now know that subchainID 9 is in bin_bflag(3), and then need to
	! find the position within "group" 3 that corresponds to subchainID 9,
	! and the corresponding binary exponent with which to bit-test our 3rd
	! binary number in bin_bflag.
	!
	!                               ___    ___    ___    _?_
	! position:                      4      3      2      1
	! binary exponent:               3      2      1      0
	!
	! In our case, the position is 1, and the corresponding binary exponent
	! is 0. That is, adding 2**0 (and only 2**0) to bin_bflag(3) would 
	! result in a 1 at position 1 in binary representation. We can
	! calculate this binary exponent using the MOD intrinsic and if
	! statement in the function below.
	!
	! If a 1 is found under bit test, then monomer(molno) is assumed to
	! be connected to the monomer in the same chain with subchainID scID.
    !
    ! Because the maximum possible integer is 2**31 with a 32 bit 
    ! representation (default Fortran integer), we set "intbits" to be
    ! 31, so the maximum exponent becomes 30. This allows something to be
    ! bonded to subchainID corresponding with exponent=intbits AND
    ! something else. Otherwise, the integer bflag would attempt to
    ! represent 2**31 + 2**(something less than 31), which is greater than
    ! 2**31 and therefore not representable. 
    !
	!----------------------------------------------------------------------
	integer function get_bondflag(molno,scID)
	implicit none
	
		integer, intent(in)  :: molno, scID
		integer :: group, expo
	
		group = ceiling(real(scID)/real(intbits))
		expo  = mod(scID,intbits) - 1
		if (expo .eq. -1) expo = intbits - 1
		
		get_bondflag  = 0	
		if (btest(monomer(molno)%bin_bflag(group),expo)) get_bondflag = 1

	end function get_bondflag

    ! Errors
    subroutine missing_bond_error(molno, bondno)
    use interfaces
    use computational_constants_MD, only: irank
    implicit none

        integer, intent (in) :: molno, bondno
        print('(a,i2,a,i10,a,i8.8,a)'), 'Bond', bondno, &
            ' missing from the bond list of monomer with glob_no', &
            monomer(molno)%glob_no, &
            '. See monomer_info_', irank,' for monomer information. Aborting.'
        call write_monomer_info(molno)
        call error_abort()

    end subroutine

	subroutine write_monomer_info(molno)
        use computational_constants_MD, only: halfdomain, npx, npy, npz, &
        domain, iblock, jblock, kblock, irank, prefix_dir
        use physical_constants_MD, only: np
        use arrays_MD, only: tag, r
        use librarymod, only: get_new_fileunit
		implicit none

		integer, intent(in) :: molno

        integer :: f
		real(kind(0.d0)) :: rglob(3)
        character(len=8) :: rankstr
        logical	:: op

        ! Get new unopened file unit
        f = get_new_fileunit() 

        write(rankstr,'(i8.8)') irank 
        open(unit=f,file=trim(prefix_dir)//"monomer_info_"//rankstr, &
        status='replace')

		rglob(1) = r(1,molno) - halfdomain(1)*(npx-1) + domain(1)*(iblock - 1)   
		rglob(2) = r(2,molno) - halfdomain(2)*(npy-1) + domain(2)*(jblock - 1)   
		rglob(3) = r(3,molno) - halfdomain(3)*(npz-1) + domain(3)*(kblock - 1)   

        write(f,'(a)') '-----------------------------------------------------'
		write(f, '(a,i8,a)') 'Monomer information for atom ', molno,':'
        write(f,'(a)') '-----------------------------------------------------'
		write(f, '(a,i8)') 'irank: ', irank
		write(f, '(a,i8)') 'molno:', molno
		write(f, '(a,f10.5,a,f10.5,a,f10.5)') 'Local position: ', &
		                                   r(1,molno),' ',r(2,molno),' ',r(3,molno) 
		write(f, '(a,f10.5,a,f10.5,a,f10.5)') 'Global position: ', &
		                                   rglob(1),' ',rglob(2),' ',rglob(3) 
		write(f, '(a,i8)') 'ChainID: '   , monomer(molno)%chainID
		write(f, '(a,i8)') 'SubchainID: ', monomer(molno)%subchainID
		write(f, '(a,i8)') 'Funcy: '     , monomer(molno)%funcy
		write(f, '(a,i8)') 'Glob_no: '   , monomer(molno)%glob_no
		write(f, '(a,4i8)') 'Bin_bflag: ' , monomer(molno)%bin_bflag
		write(f, '(a,i8)') 'Tag: '       , tag(molno) 

		if (molno.gt.np) then
		    write(f,'(a,i8,a)') 'Atom ', molno, ' is a halo particle.'
		end if

        write(f, '(a)') '____________________________________________________'
        write(f, '(a)') '____________________________________________________'
        
        close(f, status='keep')

	end subroutine write_monomer_info

    subroutine mark_chain_as_solvent(ID)
        use physical_constants_MD, only: np
        use computational_constants_MD, only : Mie_potential
        use arrays_MD, only : moltype
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
                if (Mie_potential .eq. 1) then
                    moltype(molno) = 3 !WATER
                endif

            end if

        end do

    end subroutine mark_chain_as_solvent

end module polymer_info_MD

!-------------------------------------------------------------------------------------
!----------------------------------CUDA---------------------------------------------

module CUDA_MD


	integer									:: ngpusurfacecells  !Number of suface cells on surface of GPU subdomain
	integer									:: pcellnp, size_ncells
	integer,allocatable,dimension(:),target	:: r_cellnp, r_cellnpH, r_cellnpS
	integer,allocatable,dimension(:),target	:: molnorecord, memlocrecord
	integer,allocatable,dimension(:),target	:: molnorecordS, memlocrecordS
	integer,allocatable,dimension(:),target	:: molnorecordH, memlocrecordH
	integer,allocatable,dimension(:,:)		:: gpusurfacecell	!Surface cells
	real*4 ,allocatable,dimension(:),target	:: r_cell   !Padded Domain position of molecule i
	real*4 ,allocatable,dimension(:),target	:: v_cell   !Padded Domain velocity of molecule i
	real*4 ,allocatable,dimension(:),target	:: r_cellH  !Padded Halo position of molecule i
	real*4 ,allocatable,dimension(:),target	:: v_cellH  !Padded Halo velocity of molecule i
	real*4 ,allocatable,dimension(:),target	:: r_cellS  !Padded GPU surface position of molecule i
	real*4 ,allocatable,dimension(:),target	:: v_cellS  !Padded GPU surface velocity of molecule i


end module CUDA_MD

!-------------------------------------------------------------------------------------
!---------------------------------Calculated Properties-------------------------------

module calculated_properties_MD
    !use librarymod, only: PDF

	integer	:: nplanes,gnplanes     !Number of planes used for MOP

	integer,dimension(3) 					:: nbins,nbinso,gnbins  !Number of bins to store molecular properties
	integer,dimension(:), allocatable       :: rdf_hist             !Array to keep tally of radial distribution
	integer,dimension(:,:), allocatable	    :: rdf3d_hist           !Array to keep tally of radial distribution
	real(kind(0.d0)),dimension(:), allocatable 		 :: slice_mass	    	!Array to keep tally of molecules in slice
	real(kind(0.d0)),dimension(:,:,:), allocatable 	 :: slice_massbin 		!Recorded molecules in a bin
	real(kind(0.d0)),dimension(:,:,:), allocatable	 :: volume_mass			!Mass in a control volume at time t
	real(kind(0.d0)),dimension(:,:,:), allocatable	 :: volume_mass_s		!Solvnt mass in a control volume at time t
	real(kind(0.d0)),dimension(:,:,:), allocatable	 :: volume_mass_p		!Polymer mass in a control volume at time t
	real(kind(0.d0)),dimension(:,:,:), allocatable	 :: volume_mass_pdt		!Mass in a control volume at time t - dt
	real(kind(0.d0)),dimension(:,:,:), allocatable	 :: dmdt				!Mass change in control volume from t-dt to t
	real(kind(0.d0)),dimension(:,:,:,:), allocatable :: mass_flux  			!Flow of mass over a control volume surface
	real(kind(0.d0)),dimension(:,:,:,:), allocatable :: surface_density  	!Mass located on a control volume surface

	real(kind(0.d0)) :: 	&
		binsize(3),			&		!Size of each bin
		planespacing,		&		!Spacing between planes for MOP
		vsum, mv2sum,		&		!velocity/mv2 sum
		potenergysum,		&		!potential energy sum for system
		potenergysum_LJ,	&		!potential energy sum from LJ
		potenergysum_POLY,	&		!potential energy sum from FENE
		potenergy_LJ,		&		!potential energy from LJ
		potenergy_POLY,		&		!potential energy FENE
		kinenergy, 			&		!Energies
		potenergy,			&		!Energies
		totenergy,			&  		!Energies
		virial,				&		!System properties
		temperature,		&		!System properties
		pressure,			&   	!System properties
		initialenergy,		&		!Intial energy of system
		gamma						!Parameter used in Nose Hoover shearostat

	real(kind(0.d0)), dimension(3,3) 			:: gamma_xy	 !Parameter used in Nose Hoover tensor stressostat 
	real(kind(0.d0)), dimension(:), allocatable :: &
		planes,				&  		!Location of planes used in MOP
		rdf,				&		!Radial distribution function
		diffusion,			&		!Diffusion of molecules
		meandiffusion,		&		!Time averaged diffusion of molecules
		Pxycorrel,			&     	!Sum of correlations of Pxy and Pxyzero
		slice_temperature,	&		!Temperature in a domain slice
		Pxyv_plane 	 				!Energy on plane for MOP

	real(kind(0.d0)), dimension(:,:), allocatable 	:: & 
		ssf,                &       !Static structure factor
		rdf3d,				&		!Radial distribution function
		ssf_hist,           &       !Running total for structure factor calculation
		slice_momentum,		&		!Mean velocity used in velocity slice
		Pxy_plane,			&		!Stress on plane for MOP
		Pxy,				&  		!Stress tensor for whole domain
		Pxyzero,            & 	    !Stress tensor at start of sample
        rough_array                 !Array used to setup fractal surface

	real(kind(0.d0)), dimension(:,:,:), allocatable :: & 
		rfmol,				&  		!Position(x)Force tensor per molecule
		Pxymol,				&  		!Stress tensor per molecule
		zeta_array,			&		!Local Nose Hoover Thermostat strength
		volume_temperature, &		!Temperature in a control volume at time t
		volume_energy,      &       !Energy in a control volume at time t
		Fv_ext_bin					!Power due to external forces in bins

	real(kind(0.d0)), dimension(:,:,:,:), allocatable	:: &
		volume_momentum,	& 		!Momentum in a control volume at time t
		volume_momentum_s,	& 		!Solvent momentum in a control volume at time t
		volume_momentum_p,	& 		!Polymer momentum in a control volume at time t
		centre_of_mass,	    & 		!Centre of mass of control volume
		energy_flux,		&		!Flow of energy over a control volume surface
		Pxyvface,			&		!Power tensor on bin face
		Pxyvface_mdt,		&		!Power tensor on bin face at previous timestep
		Pxyvface_integrated,&		!Integrated form of Power
		F_ext_bin,          &		!External Force per bin
		evbin,              &  		!velocity energy per bin
        heatfluxbin

	real(kind(0.d0)), dimension(:,:,:,:,:), allocatable :: & 
		volume_force,  		& 		!Force acting over control volume surface 
		momentum_flux, 		&		!Flow of momentum over a control volume surface
		rfbin, 				& 		!Position(x)Force tensor per bin
		rfvbin,				& 		!Position(x)Force dot v per bin
		vvbin, 				& 		!velocity(x)velocity tensor per bin
		Pxybin, 			&		!Stress tensor per bin
		Pxyface, 			&		!Stress tensor on bin face
		Gxybins	    				!Parameter used in Nose Hoover stressostat

	!real(kind(0.d0)),dimension(2,3,44)	:: shiftVAstress
contains
 
	function get_mass_slices(ixyz)
	use computational_constants_MD, only: domain,halfdomain
	use physical_constants_MD, only: np
	use arrays_MD, only: r
	implicit none

		integer, intent(in) 				:: ixyz
		integer								:: bin,n
		real(kind(0.d0)) 					:: binsize
		integer, dimension(nbins(ixyz)) 	:: get_mass_slices
		
		binsize = domain(ixyz)/nbins(ixyz)
		get_mass_slices = 0
		do n=1,np
			bin = ceiling((r(ixyz,n)+halfdomain(ixyz))/binsize)
			if (bin.lt.1) 				bin = 1
			if (bin.gt.nbins(ixyz)) 	bin = nbins(ixyz)
			get_mass_slices(bin) = get_mass_slices(bin) + 1
		end do

	end function get_mass_slices
 
	function get_velo_slices(ixyz)
	use arrays_MD, only: r,v,a
	use computational_constants_MD
	use physical_constants_MD, only: np,nd
	implicit none

		integer, intent(in) 				:: ixyz
		integer								:: bin,n
		real(kind(0.d0)) 					:: binsize
		real(kind(0.d0)), dimension(nbins(ixyz),nd) :: get_velo_slices
		
		binsize = domain(ixyz)/nbins(ixyz)
		get_velo_slices = 0
		do n=1,np
			bin = ceiling((r(ixyz,n)+halfdomain(ixyz))/binsize)
			if (bin.lt.1) 				bin = 1
			if (bin.gt.nbins(ixyz)) 	bin = nbins(ixyz)
			select case (integration_algorithm)
			case (leap_frog_verlet) 
				get_velo_slices(bin,:) = get_velo_slices(bin,:) + v(:,n) + 0.5d0*a(:,n)*delta_t
			case (velocity_verlet)
				get_velo_slices(bin,:) = get_velo_slices(bin,:) + v(:,n)
			end select
		end do

	end function get_velo_slices

	real(kind(0.d0)) function get_temperature_PUT()
	use computational_constants_MD, only: domain, halfdomain,delta_t
	use physical_constants_MD, only: np,nd,globalnp
	use arrays_MD, only: r,v,a
	use shear_info_MD, only: le_sp
	implicit none
	
		integer	:: slicebin,n	
		integer, dimension(:), allocatable 	:: m_slice
		real(kind(0.d0)) :: pec_mv2sum
		real(kind(0.d0)), dimension(nd) 	:: slicebinsize, vel
		real(kind(0.d0)), dimension(:,:), allocatable :: v_slice,v_avg
		
		slicebinsize(:) = domain(:)/nbins(:)				! Get bin size
		pec_mv2sum 		= 0.d0								! Initialise

        print*, "WARNING -- this get_temperature_PUT routine is old -- mass not included"
		
		!todo reevaluate
		allocate(m_slice(nbins(le_sp)))	
		allocate(v_slice(nbins(le_sp),nd))
		allocate(v_avg(nbins(le_sp),nd))
		m_slice = get_mass_slices(le_sp)				! Get total mass in all slices
		v_slice = get_velo_slices(le_sp)				! Get total velocity in all slices
		do slicebin=1,nbins(le_sp)
			v_avg(slicebin,:) = v_slice(slicebin,:)/m_slice(slicebin)
		end do

		do n=1,np
			slicebin 			= ceiling((r(le_sp,n)+halfdomain(le_sp))/slicebinsize(le_sp))
			if (slicebin > nbins(le_sp)) slicebin = nbins(le_sp)	! Prevent out-of-range values
			if (slicebin < 1) slicebin = 1										! Prevent out-of-range values
			vel(:) 			 	= v(:,n) - v_avg(slicebin,:) - 0.5d0*a(:,n)*delta_t
			pec_mv2sum 			= pec_mv2sum + dot_product(vel,vel)
		end do

		call globalSum_(pec_mv2sum)
	
		get_temperature_PUT = pec_mv2sum / real(nd*globalnp,kind(0.d0))

		deallocate(m_slice)
		deallocate(v_slice)
		deallocate(v_avg)

	end function get_temperature_PUT

end module calculated_properties_MD

module boundary_MD
#if __INTEL_COMPILER > 1200
    use librarymod, only: PDF
#endif
	!Boundary force flag and parameters
	integer,          dimension(6) :: bforce_flag
	real(kind(0.d0)), dimension(6) :: bforce_dxyz
    integer                        :: specular_wall_flag
	real(kind(0.d0)), dimension(3) :: specular_wall

	! Specular wall easy-read parameters
	integer :: specular_flag
	integer, parameter ::       &
		specular_off = 0,       &
		specular_flat = 1,      &
		specular_radial = 2

    ! Open Boundary flags
    integer, dimension(6) :: open_boundary 

	! Boundary force easy-read parameters
	integer,          parameter    :: &
		bforce_off = 0,         &
		bforce_OT = 1,          &
		bforce_NCER = 2,        &
		bforce_Flekkoy = 3,     &
        bforce_pdf_input = 4,   &
        substrate_force = 5

    ! Measure boundary force with pdf
    integer :: bforce_pdf_measure
    integer :: bforce_pdf_nsubcells
    integer :: bforce_pdf_nbins
    integer :: bforce_pdf_Nave
    real(kind(0.d0)) :: bforce_pdf_min
    real(kind(0.d0)) :: bforce_pdf_max
    real(kind(0.d0)) :: bforce_pdf_binsize
#if __INTEL_COMPILER > 1200
    type(PDF), dimension(:,:), allocatable :: bforce_pdf
#endif
    real(kind(0.d0)), allocatable :: bforce_pdf_input_data(:,:,:)

end module boundary_MD

