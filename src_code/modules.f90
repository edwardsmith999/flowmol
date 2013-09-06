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
	integer,dimension(:),allocatable:: procnp 				!Array of all processors np
	integer,dimension(:),allocatable:: proctethernp 		!Array of all processors np
	double precision   				:: volume, density      !Define constant volume and density
	double precision   				:: rcutoff, halfrcutoff !Cut off distance for particle interactions
	double precision   				:: rcutoff2             !Cut off distance for particle interactions squared
	double precision   				:: potshift		   		!Shift in Lennard Jones potential due to cutoff
	double precision   				:: potential_sLRC 		!Long range potential correction 
	double precision   				:: pressure_sLRC 		!Long range pressure correction 
	double precision   				:: inputtemperature     !Define initial temperature
	double precision   				:: initialunitcell      !Initial size of unit cell
	double precision   				:: initialvel           !Initial velocity of particles
	double precision,parameter 		:: pi=4.d0*atan(1.d0)

	!Simulation setup conditions
	double precision, dimension(3)	:: fixdisttop, slidedisttop, fixdistbottom, slidedistbottom, wallslidev
	double precision, dimension(3)	:: tethereddisttop, tethereddistbottom, thermstattop,thermstatbottom

end module physical_constants_MD

!------------------------------------------------------------------------------
!-----------------------------Computational Constants--------------------------
! Constants used to determine the computational parameters of the simulation
module computational_constants_MD

	!Command-line arguments
	logical					:: restart
	character(len=200) 		:: input_file, initial_microstate_file

	!Force and Potential flags
	integer						  :: force_list	   		 !flag for neighbr/cell list
	integer						  :: potential_flag  	 !Choose LJ or Polymer potential
	integer                 	  :: tether_flag     	 !True if there exists 
	integer                 	  :: external_force_flag !Apply external forces?
	integer                 	  :: F_ext_ixyz			 !Direction of external forces
	double precision        	  :: F_ext				 !Magnitude of external forces
	double precision,dimension(6) :: F_ext_limits		 !Limits of region external forces applied to


	! Move particle tags
	logical 		   :: tag_thermostat_active
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

	! Wall texture flags
	integer			   :: texture_type, texture_therm
	integer, parameter :: posts = 1
	integer, parameter :: roughness = 2
	integer, parameter :: converge_diverge = 3
	double precision   :: texture_intensity
	

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
		tag_move    = 6

	!Initial configuration selection
	integer            :: initial_config_flag
	character(len=128) :: config_special_case

	double precision	:: liquid_density	!Density of liquid if solid/liquid case used

	!Input (on or off) flags
	integer	:: & 
		vmd_outflag, vmd_start, vmd_finish, vmd_count=1, &
		macro_outflag, &
		sLRC_flag,	&
		mass_outflag, &
		velocity_outflag, &
		temperature_outflag, &
		pressure_outflag, &
		viscosity_outflag, &
		rdf_outflag, &
		rtrue_flag, &
		prev_rtrue_flag, &
		ssf_outflag, &
		vdist_flag, & 
		cv_conserve,	&
		mflux_outflag, &
		vflux_outflag, &
		eflux_outflag, &
		proc_reorder, &				!Reorder processors at restart
		pass_vhalo = 0, &
		fixed_rebuild_flag, &		!Fixed rebuild flag
		peculiar_flag	 			!Take streaming velocity away from temperature 	

	integer, dimension(3)	:: periodic

	!Boundary force flag and parameters
	integer,          dimension(6) :: bforce_flag
	real(kind(0.d0)), dimension(6) :: bforce_dxyz
	real(kind(0.d0)), dimension(3) :: specular_wall

	! Specular wall easy-read parameters
	integer :: specular_flag
	integer, parameter ::       &
		specular_off = 0,       &
		specular_flat = 1,      &
		specular_radial = 2

	! Boundary force easy-read parameters
	integer,          parameter    :: &
		bforce_off = 0,         &
		bforce_OT = 1,          &
		bforce_NCER = 2,        &
		bforce_Flekkoy = 3,     &
		bforce_rdf = 4

	!Tethered molecules spring potential coefficients
	double precision    :: teth_k2,teth_k4,teth_k6

	!Parameters
	integer	:: & 
		nproc, npx, npy, npz, 	&	!Number of MPI ranks and cartesian topology sizes
		iter, 					&	!Global simulation iteration count
		tplot, 					&	!Frequency at which to record results
		teval, 					&	!Frequency at which to evaluate results
		Nmass_ave, 				&	!Number of averages for each mass average
		Nvel_ave, 				&	!Number of averages for each velocity average
		NTemp_ave, 				&	!Number of averages for each temperature measurement
		Nstress_ave, 			&	!Number of bins for stress calculation
		split_kin_config, 		&	!Flag to determine if kinetic and configurational stress seperated
		Nvisc_ave, 				&	!Number of samples for viscosity measurement
		Nmflux_ave, 			&	!Number of averages for each mass flux
		Nvflux_ave, 			&	!Number of averages for each velocity flux
		Neflux_ave, 			&	!Number of averages for each energy flux
		initialstep, 			&	!Initial step of simulation
		finalstep,              &   !Final step of simulation
		Nsteps, 				&	!Total number of computational steps
		initialise_steps, 		&	!Initialisation steps to run before simulation start
		fixed_rebuild,          &   !Fixed rebuild frequency
		extralloc, 				&	!Extra allocation space to include copied halos
		overlap, 				&	!Size of overlap region used to apply force to molecular region
		Nvmd_intervals,         &
		rdf_nbins,              &   !Number of discrete "bins" used to calculate g(r)
		ssf_ax1,                &   !1st projection axis for static structure factor
		ssf_ax2,                &
		ssf_nmax                    !Maximum wavenumber for S(k) calculation

	integer,dimension(:,:),allocatable	:: vmd_intervals			!Multiple intervals for vmd record

	double precision 	:: delta_t           !Size of timestep for each computational step
	double precision 	:: elapsedtime       !Total time elapsed including restarts
	double precision 	:: simtime=0.d0      !Total incremented simulation time
	double precision 	:: rneighbr2         !Square of rcuttoff+delta_rneighbr
	double precision 	:: delta_rneighbr    !Radius used for neighbour list construction
	double precision    :: rdf_rmax          !Maximum radius for radial distribution function
	double precision	:: rescue_snapshot_freq	!Rescue snapshot output frequency in seconds
	
	!Store surface bins of processes subdomain for outputs over periodic boundaries
	integer									:: nsurfacebins     !Number of surface bins
	integer,allocatable,dimension(:,:)	    :: surfacebins		!Surface Bins

	integer									:: nsurfacecells     !Number of surface bins
	integer,allocatable,dimension(:,:)	    :: surfacecells		!Surface Bins

	integer		 							:: nhalobins 	!Number of halo bins
	integer,allocatable,dimension(:,:)	    :: halobins		!halo Bins

	integer		 							:: nhalocells 	!Number of halo bins
	integer,allocatable,dimension(:,:)	    :: halocells		!halo Bins

	!Number and size of unit used for initial setup of molecules (i.e. FCC unit)
	integer,          dimension(3)		:: initialnunits
	double precision, dimension(3) 		:: initialunitsize

	!Size of global computational domain and domain per processor
	double precision, dimension(3)	:: globaldomain
	double precision, dimension(3)	:: domain, halfdomain

	!Number and size of cells used for domain subdivision into cells of size ~rcutoff
	integer,          dimension(3)	:: ncells
	double precision, dimension(3)	:: cellsidelength, halfcellsidelength

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

	!Number of samples used to calculate l_ij in VA stress calculation
	integer :: VA_line_samples

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
	double precision 	:: le_sv                    !Shear velocity
	double precision 	:: le_sr                    !Shear rate
	double precision	:: le_sx                    !Shear distance
	double precision	:: le_st                    !Shear time

	integer, dimension(:), allocatable	:: mol_wrap_integer		

end module shear_info_MD

!------------------------------------------------------------------------------
!----------------------------------Arrays--------------------------------------

module arrays_MD

	integer,          dimension(:),   allocatable, target	:: tag !Mol tags
	integer, 	  	  dimension(:,:), allocatable, target	:: &
		fix                         !Fixed molecules
	double precision, dimension(:),   allocatable, target 	:: &
		potenergymol, 		&		!Potential energy of each molecule
		potenergymol_LJ, 	&		!LJ Potential energy of each molecule
		potenergymol_FENE,	&		!FENE Potential energy of each molecule
		virialmol,			&		!Virial of each molecule
		recvbuffer
	double precision, dimension(:,:),   allocatable			:: &
		rtrue, 		&      			!Positions with no period BC
		vtrue,      &               !Corresponding velocities
		rtether, 	&
		rijsum, 	&				!Sum of all molecular rij values
		theta, 		&
		aD,aR
	double precision, dimension(:,:),   allocatable, target 	:: &
		r, 		&        		  	!Positions
		v, 		&        		  	!Velocity
		a, 		&					!Accelerations
		slidev, &					!Speed for sliding molecules
		U							!Streaming velocity of molecules

end module arrays_MD


!------------------------------------------------------------------------------
!----------------------------------Linked List---------------------------------

module linked_list

	!3D Information for each cell and pointer to cell/bin linked list
	type cellinfo
		integer  , dimension(:,:,:), allocatable :: cellnp
		type(ptr), dimension(:,:,:), pointer     :: head
	end type cellinfo
	
	!Cannot create array of pointers in Fortran so create an array of
	!type ptr which contains a pointer
	type ptr
		type(node), pointer :: point
	end type ptr

	!A Node of the linked list which is built for each cell/bin
	type node
		double precision,dimension(:), pointer 	:: rp
		double precision,dimension(:), pointer 	:: vp
		double precision,dimension(:), pointer	:: ap
		integer             :: molno
		type(node), pointer :: next, previous
	end type node

	! Neighbourlist with cell head 
	!Information for pointer to neighbour list
	type neighbrinfo
		integer,dimension(:),allocatable		:: noneighbrs
		type(ptr_nbr), dimension(:), pointer	:: head
	end type neighbrinfo
	
	!Cannot create array of pointers in Fortran so create an array of
	!type ptr which contains a pointer
	type ptr_nbr
		type(neighbrnode), pointer 	:: point
	end type ptr_nbr

	!A Node of the linked list which is built for the neighbourlist
	type neighbrnode
		integer			   			:: molnoj
		type(neighbrnode), pointer 	:: next, previous
	end type neighbrnode


	!Information for pass list
	type passinfo
		integer					:: sendnp
		type(passnode), pointer :: head
	end type passinfo

	!A Node of the passed molecule list
	type passnode
		integer                 :: molno
		integer					:: ipass
		integer					:: jpass
		integer					:: kpass
		type(passnode), pointer :: next, previous
	end type passnode

	!Global type cell used to record location of all cells linked lists
	!N.B. Type cell is not a pointer & nothing points to it
	type(cellinfo)    	:: cell
	type(cellinfo)    	:: bin	!Used to bin molecules for statistics
	type(neighbrinfo) 	:: neighbour
	type(passinfo)    	:: pass

end module linked_list

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

	integer, dimension(:,:,:),   allocatable :: cyl_mass
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
	integer             :: r_gyration_outflag       !Radius of gyration outflag
	integer             :: r_gyration_iter0         !Iteration at which to start recording R_g

	integer             :: solvent_flag             !Solvent on/off flag
	integer             :: solvent_ratio            !Solvent concentration
	double precision    :: eps_pp, eps_ps, eps_ss   !Solvent parameters
	double precision, parameter :: sod_cut =1.5d0   !Soddemann potential cutoff
	double precision, parameter :: sod_cut2=2.25d0  ! 
	
	double precision 	:: k_c, R_0					!Spring constant and max elongation of bonds	
	double precision 	:: etevtcf                  !End-to-end vector time correlation function
	double precision    :: R_g                      !Radius of gyration
	
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

	integer, parameter :: nsdmi=8                   !Number of sent doubles for monomer_info 

	integer :: intbits                              !Size in bits of integer on this machine

	type(monomer_info), dimension(:), allocatable :: monomer
	!eg. to access chainID of mol 23, call monomer(23)%chainID

	double precision, dimension(:,:), allocatable :: etev
	double precision, dimension(:,:), allocatable :: etev_0 !End-to-end vectors for polymer chains at iter=etevtcf_iter0 (size of np,only stored for leftmost chain mols) 

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

	integer	:: nplanes,gnplanes     !Number of planes used for MOP

	integer,dimension(3) :: nbins,nbinso,gnbins	                    !Number of bins to store molecular properties
	integer,dimension(:), allocatable 		:: vfd_bin        		!Array to keep tally of molecular velocity distribution
	integer,dimension(:), allocatable       :: rdf_hist             !Array to keep tally of radial distribution
	integer,dimension(:,:), allocatable	    :: rdf3d_hist           !Array to keep tally of radial distribution
	integer,dimension(:), allocatable 		:: slice_mass	    	!Array to keep tally of molecules in slice
	integer,dimension(:,:,:), allocatable 	:: slice_massbin 		!Recorded molecules in a bin
	integer,dimension(:,:,:), allocatable	:: volume_mass			!Mass in a control volume at time t
	integer,dimension(:,:,:), allocatable	:: volume_mass_pdt		!Mass in a control volume at time t - dt
	integer,dimension(:,:,:), allocatable	:: dmdt					!Mass change in control volume from t-dt to t
	integer,dimension(:,:,:,:), allocatable	:: mass_flux  			!Flow of mass over a control volume surface

	double precision :: 	&
		binsize,			&		!Size of each bin
		planespacing,		&		!Spacing between planes for MOP
		vsum, v2sum,		&		!velocity sum
		potenergysum,		&		!potential energy sum for system
		potenergysum_LJ,	&		!potential energy sum from LJ
		potenergysum_FENE,	&		!potential energy sum from FENE
		potenergy_LJ,		&		!potential energy from LJ
		potenergy_FENE,		&		!potential energy FENE
		kinenergy, 			&		!Energies
		potenergy,			&		!Energies
		totenergy,			&  		!Energies
		virial,				&		!System properties
		temperature,		&		!System properties
		pressure,			&   	!System properties
		initialenergy,		&		!Intial energy of system
		zeta,				&		!Parameter used in Nose Hoover thermostat
		gamma						!Parameter used in Nose Hoover shearostat

	double precision, dimension(3,3) 			:: gamma_xy	 !Parameter used in Nose Hoover tensor stressostat 
	double precision, dimension(:), allocatable :: &
		planes,				&  		!Location of planes used in MOP
		rdf,				&		!Radial distribution function
		normalisedvfd_bin,	&		!Bin normalised so sum of all bins is one
		diffusion,			&		!Diffusion of molecules
		meandiffusion,		&		!Time averaged diffusion of molecules
		Pxycorrel,			&     	!Sum of correlations of Pxy and Pxyzero
		slice_temperature,	&		!Temperature in a domain slice
		Pxyv_plane 	 				!Energy on plane for MOP

	double precision, dimension(:,:), allocatable 	:: & 
		ssf,                &       !Static structure factor
		rdf3d,				&		!Radial distribution function
		ssf_hist,           &       !Running total for structure factor calculation
		slice_momentum,		&		!Mean velocity used in velocity slice
		Pxy_plane,			&		!Stress on plane for MOP
		Pxy,				&  		!Stress tensor for whole domain
		Pxyzero   					!Stress tensor at start of sample

	double precision, dimension(:,:,:), allocatable :: & 
		rfmol,				&  		!Position(x)Force tensor per molecule
		Pxymol,				&  		!Stress tensor per molecule
		zeta_array,			&		!Local Nose Hoover Thermostat strength
		volume_temperature			!Temperature in a control volume at time t

	double precision, dimension(:,:,:,:), allocatable	:: &
		volume_momentum,	& 		!Momentum in a control volume at time t
		volume_momentum_pdt,& 		!Momentum in a control volume at time t-dt
		dmvdt,				&		!Momentum change in control volume from t-dt to t
		energy_flux,		&		!Flow of energy over a control volume surface
		Pxyvface,			&		!Stress tensor on bin face
		F_ext_bin					!External Force per bin

	double precision, dimension(:,:,:,:,:), allocatable :: & 
		volume_force,  		& 		!Force acting over control volume surface 
		momentum_flux, 		&		!Flow of momentum over a control volume surface
		rfbin, 				& 		!Position(x)Force tensor per bin
		vvbin, 				& 		!velocity(x)velocity tensor per bin
		Pxybin, 			&		!Stress tensor per bin
		Pxyface, 			&		!Stress tensor on bin face
		Gxybins	    				!Parameter used in Nose Hoover stressostat
	


	!double precision,dimension(2,3,44)	:: shiftVAstress
contains
 
	function get_mass_slices(ixyz)
	use computational_constants_MD, only: domain,halfdomain
	use physical_constants_MD, only: np
	use arrays_MD, only: r
	implicit none

		integer, intent(in) 				:: ixyz
		integer								:: bin,n
		double precision 					:: binsize
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
		double precision 					:: binsize
		double precision, dimension(nbins(ixyz),nd) :: get_velo_slices
		
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

	double precision function get_temperature_PUT()
	use computational_constants_MD, only: domain, halfdomain,delta_t
	use physical_constants_MD, only: np,nd,globalnp
	use arrays_MD, only: r,v,a
	use shear_info_MD, only: le_sp
	implicit none
	
		integer	:: slicebin,n	
		integer, dimension(:), allocatable 	:: m_slice
		double precision :: pec_v2sum
		double precision, dimension(nd) 	:: slicebinsize, vel
		double precision, dimension(:,:), allocatable :: v_slice,v_avg
		
		slicebinsize(:) = domain(:)/nbins(:)				! Get bin size
		pec_v2sum 		= 0.d0								! Initialise
		
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
			pec_v2sum 			= pec_v2sum + dot_product(vel,vel)
		end do

		call globalSum(pec_v2sum)
	
		get_temperature_PUT = pec_v2sum / real(nd*globalnp,kind(0.d0))

		deallocate(m_slice)
		deallocate(v_slice)
		deallocate(v_avg)

	end function get_temperature_PUT

end module calculated_properties_MD



!Fortran object oriented solution to CV conservation checking
module CV_objects
	implicit none

type :: check_CV_mass
	integer,dimension(:,:,:,:),allocatable 	:: flux
	integer,dimension(:,:,:),allocatable	:: dXdt, X, X_minus_t, X_minus_2t
	contains
		procedure :: initialise  => initialise_mass
		procedure :: update_dXdt => update_dXdt_mass
		procedure :: check_error => check_error_mass
end type check_CV_mass

type :: check_CV_momentum
	double precision,dimension(:,:,:,:,:),allocatable 	:: flux, Pxy,  Pxy_minus_t
	double precision,dimension(:,:,:,:),allocatable		:: dXdt, X, X_minus_t, X_minus_2t, F_ext
	contains
		procedure :: initialise  => initialise_momentum
		procedure :: update_dXdt => update_dXdt_momentum
		procedure :: update_flux => update_flux_momentum
		procedure :: update_Pxy  => update_Pxy
		procedure :: update_F_ext=> update_F_ext
		procedure :: check_error => check_error_momentum
end type check_CV_momentum

	!Check CV conservation
	logical	:: CV_debug=.false.
	type(check_CV_mass)		:: CVcheck_mass		! declare an instance of CV checker
	type(check_CV_momentum)	:: CVcheck_momentum		! declare an instance of CV checker

contains

!===================================================
!	   M a s s   C o n t r o l   V o l u m e
!===================================================

	!Constructor for object
	subroutine initialise_mass(self, nb)
		implicit none

		! initialize shape objects
		class(check_CV_mass) 		:: self

		integer, dimension(3),intent(in) :: nb

		allocate(self%flux(nb(1),nb(2),nb(3),6))
		allocate(self%dXdt(nb(1),nb(2),nb(3)))
		allocate(self%X(nb(1),nb(2),nb(3)))
		allocate(self%X_minus_t(nb(1),nb(2),nb(3)))
		allocate(self%X_minus_2t(nb(1),nb(2),nb(3)))

		self%flux 		= 0
		self%dXdt 		= 0
		self%X 			= 0
		self%X_minus_t 	= 0
		self%X_minus_2t = 0

	end subroutine initialise_mass

	!Update time evolution and store previous two values
	subroutine update_dXdt_mass(self, X)
		implicit none

		! initialize shape objects
		class(check_CV_mass) :: self

		integer,dimension(:,:,:),intent(in) :: X

		!self%X_minus_2t = self%X_minus_t
		self%X_minus_t  = self%X
		self%X 		  = X

		self%dXdt = self%X - self%X_minus_t

	end subroutine update_dXdt_mass


	!Check error for specified range of bins
	subroutine check_error_mass(self,imin,imax,jmin,jmax,kmin,kmax,iter,irank)
		implicit none

		! initialize shape objects
		class(check_CV_mass) :: self

		integer,intent(in) :: iter,irank,imin,imax,jmin,jmax,kmin,kmax
		integer :: i,j,k

		logical,save :: first_time = .true., check_ok = .false.

		!First call doesn't have difference in time yet so skip
		if (first_time) then
			first_time = .false.
			return
		endif

		check_ok = .true.

		do i = imin,imax
		do j = jmin,jmax
		do k = kmin,kmax

			if(sum(self%flux(i,j,k,:))-self%dXdt(i,j,k) .ne. 0) then
				print'(a,i8,4i4,4i8)','Error in mass flux', iter,irank,i,j,k, & 
					sum(self%flux(i,j,k,:)),self%dXdt(i,j,k),self%X_minus_t(i,j,k),self%X(i,j,k)
				check_ok = .false.
			endif

		enddo
		enddo
		enddo

		if (check_ok .eq. .false.) then
			stop "Error in mass flux"
		endif

	end subroutine check_error_mass

	!Returns a CV which is the size of the specified range of CV
	!subroutine coarse_grain_CV_mass(self,imin,imax,jmin,jmax,kmin,kmax,iter,irank)
	!	implicit none

	!	class(check_CV_mass) :: self
	!	integer,intent(in) :: iter,irank,imin,imax,jmin,jmax,kmin,kmax



	!end subroutine coarse_grain_CV_mass


!===================================================
!	M o m e n t u m   C o n t r o l   V o l u m e
!===================================================

	!Constructor for object
	subroutine initialise_momentum(self, nb)
		implicit none

		! initialize shape objects
		class(check_CV_momentum) 		:: self

		integer, dimension(3),intent(in) :: nb

		allocate(self%flux(nb(1),nb(2),nb(3),3,6))
		allocate(self%Pxy(nb(1),nb(2),nb(3),3,6))
		allocate(self%Pxy_minus_t(nb(1),nb(2),nb(3),3,6))
		allocate(self%dXdt(nb(1),nb(2),nb(3),3))
		allocate(self%X(nb(1),nb(2),nb(3),3))
		allocate(self%X_minus_t(nb(1),nb(2),nb(3),3))
		allocate(self%F_ext(nb(1),nb(2),nb(3),3))
		!allocate(self%X_minus_2t(nb(1),nb(2),nb(3),3))

		self%flux 		= 0.d0
		self%Pxy  		= 0.d0
		self%Pxy_minus_t= 0.d0
		self%dXdt 		= 0.d0
		self%X 			= 0.d0
		self%X_minus_t 	= 0.d0
		self%F_ext		= 0.d0
		!self%X_minus_2t = 0.d0

	end subroutine initialise_momentum

	!Update time evolution and store previous two values
	subroutine update_dXdt_momentum(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_momentum) :: self

		double precision,dimension(:,:,:,:),intent(in) :: X

		!self%X_minus_2t = self%X_minus_t
		self%X_minus_t  = self%X
		self%X 		  = X

		self%dXdt = self%X - self%X_minus_t

	end subroutine update_dXdt_momentum

	!Update time evolution and store previous two values
	subroutine update_F_ext(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_momentum) :: self

		double precision,dimension(:,:,:,:),intent(in) :: X

		self%F_ext = X

	end subroutine update_F_ext

	!Update time evolution and store previous two values
	subroutine update_flux_momentum(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_momentum) :: self

		double precision,dimension(:,:,:,:,:),intent(in) :: X

		self%flux = X

	end subroutine update_flux_momentum

	!Update time evolution and store previous two values
	subroutine update_Pxy(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_momentum) :: self

		double precision,dimension(:,:,:,:,:),intent(in) :: X

		self%Pxy_minus_t = self%Pxy
		self%Pxy = X

	end subroutine update_Pxy


	!Check error for specified range of bins
	subroutine check_error_momentum(self,imin,imax,jmin,jmax,kmin,kmax,iter,irank)
		! initialize shape objects
		use computational_constants_MD, only : domain,delta_t,Nvflux_ave
		use calculated_properties_MD, only : nbins
		implicit none

		class(check_CV_momentum) :: self

		integer,intent(in) :: iter,irank,imin,imax,jmin,jmax,kmin,kmax

		integer :: i,j,k
		logical,save 					:: first_time = 0
		double precision				:: conserved
		double precision,dimension(3)	:: binsize,totalpressure,totalflux,F_ext,dvelocitydt

		!First call doesn't have difference in time yet so skip
		if (first_time .lt. 2) then
			first_time = first_time + 1
			return
		endif

		binsize = domain/nbins

		!print'(2i4,f13.5,3i8)', iter, irank, maxval(self%F_ext), maxloc(self%F_ext)
		!print'(a,4i8,4f13.8)', 'Inside  object', irank,6,6,6, self%F_ext(6,6,6,:),sum(self%F_ext(6,6,6,:))

		do i = imin,imax
		do j = jmin,jmax
		do k = kmin,kmax


		    !Calculate total CV flux and change in mass
		    totalflux =(self%flux(i,j,k,:,1)+self%flux(i,j,k,:,4))/binsize(1) &
		              +(self%flux(i,j,k,:,2)+self%flux(i,j,k,:,5))/binsize(2) &
		              +(self%flux(i,j,k,:,3)+self%flux(i,j,k,:,6))/binsize(3)

		    !Totalpressure = totalpressure*delta_t
		    totalpressure =(self%Pxy_minus_t(i,j,k,:,1)-self%Pxy_minus_t(i,j,k,:,4))/binsize(1) &
		                  +(self%Pxy_minus_t(i,j,k,:,2)-self%Pxy_minus_t(i,j,k,:,5))/binsize(2) &
		                  +(self%Pxy_minus_t(i,j,k,:,3)-self%Pxy_minus_t(i,j,k,:,6))/binsize(3)
			F_ext = self%F_ext(i,j,k,:)/product(binsize)

			!drhou/dt
		    dvelocitydt =  self%dxdt(i,j,k,:)/(delta_t*Nvflux_ave)

		    !Verify that CV momentum is exactly conservative
		    conserved = sum(totalpressure-totalflux-dvelocitydt-F_ext)
			if(conserved .gt. 0.000000001d0) then
				!print'(a,i8,3f13.8)', 'Inside  object', irank, self%F_ext(i+1,j+1,k+1,:)
				!if (maxval(abs(self%F_ext(i-1:i+1,j-1:j+1,k-1:k+1,:))) .ne. 0.0000000001) &
				!print'(a,i8,4f13.8,3i8)', 'Inside  object', irank, maxval(abs(self%F_ext(i-1:i+1,j-1:j+1,k-1:k+1,1))), & 
				!							    maxval(abs(self%F_ext(i-1:i+1,j-1:j+1,k-1:k+1,2))), & 
				!							    maxval(abs(self%F_ext(i-1:i+1,j-1:j+1,k-1:k+1,3))), &
				!								   sum(self%F_ext(i-1:i+1,j-1:j+1,k-1:k+1,:)), &
				!								maxloc(abs(self%F_ext(i-1:i+1,j-1:j+1,k-1:k+1,1)))
				print'(a,i8,4i4,7f10.5)','Error in momentum flux', iter,irank,i,j,k, & 
					 conserved, sum(totalpressure),-sum(totalflux),sum(dvelocitydt), & 
					+sum(F_ext), sum(self%X(i,j,k,:)),   & 
					 sum(self%X_minus_t(i,j,k,:))
			endif

		enddo
		enddo
		enddo

	end subroutine check_error_momentum


	!Returns a CV which is the size of the specified range of CV
    ! Calculate total CV flux and change in momentum - NOTE binsize used
    ! as during sum all dx/dy/dz are identical and cancel leaving
    ! a single dx,dy or dz
	!subroutine coarse_grain_CV_momentum(self,imin,imax,jmin,jmax,kmin,kmax,iter,irank)
	!	implicit none

	!	class(check_CV_momentum) :: self
	!	integer,intent(in) :: iter,irank,imin,imax,jmin,jmax,kmin,kmax


		

	!end subroutine coarse_grain_CV_momentum


end module CV_objects
