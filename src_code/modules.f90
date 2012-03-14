!=====================================================================================
!------------------------------DEFINE MODULES-----------------------------------------
!=====================================================================================

!-------------------------------------------------------------------------------------
!-----------------------------Physical Constants---------------------------------
!Constants which represent physical values and consequently change the physical
!conditions of the simulation
module physical_constants_MD

	integer, parameter 				:: nd = 3		   		!Number of dimensions
	integer            				:: np                   !Number of particles
	integer		   					:: globalnp             !Global number of particles
	integer		   					:: halo_np              !Number of molecules in halo
	integer,dimension(:),allocatable:: procnp 				!Array of all processors np
	double precision   				:: volume, density      !Define constant volume and density
	double precision   				:: rcutoff, halfrcutoff !Cut off distance for particle interactions
	double precision   				:: rcutoff2             !Cut off distance for particle interactions squared
	double precision   				:: potshift		   		!Shift in Lennard Jones potential due to cutoff
	double precision   				:: inputtemperature     !Define initial temperature
	double precision   				:: initialunitcell      !Initial size of unit cell
	double precision   				:: initialvel           !Initial velocity of particles
	double precision,parameter 		:: pi=4.d0*atan(1.d0)

	!Simulation setup conditions
	double precision, dimension(3)	:: fixdisttop, slidedisttop, fixdistbottom, slidedistbottom, wallslidev
	double precision, dimension(3)	:: tethereddisttop, tethereddistbottom, thermstattop,thermstatbottom

end module physical_constants_MD

!-------------------------------------------------------------------------------------
!-----------------------------Computational Constants---------------------------------
!Constants used to determine the computational parameters of the simulation

module computational_constants_MD

	!Command-line arguments
	logical					:: restart
	character(len=200) 		:: input_file, initial_microstate_file

	!Force and Potential flags
	integer					:: force_list		!flag to switch between neighbour and cell list
	integer					:: potential_flag	!Choose LJ or Polymer Potential 

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

	!Thermostat flag
	integer					:: thermstat_flag

	!Input (on or off) flags
	integer		:: & 
		vmd_outflag, vmd_start, vmd_finish, vmd_count=1, &
		macro_outflag, &
		mass_outflag, &
		velocity_outflag, &
		temperature_outflag, &
		pressure_outflag, &
		viscosity_outflag, &
		mflux_outflag, &
		vflux_outflag, &
		eflux_outflag, &
		proc_reorder					!Reorder processors at restart
	integer, dimension(3)	:: periodic


	!Parameters
	integer	:: & 
		nproc, npx, npy, npz, 	&	!Number of MPI ranks and cartesian topology sizes
		iter, 					&	!Global simulation iteration count
		tplot, 					&	!Frequency at which to record results
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
		Nsteps, 				&	!Total number of computational steps
		initialise_steps, 		&	!Initialisation steps to run before simulation start
		extralloc, 				&	!Extra allocation space to include copied halos
		overlap, 				&	!Size of overlap region used to apply force to molecular region
		Nvmd_intervals
	integer,dimension(:,:),allocatable	:: vmd_intervals			!Multiple intervals for vmd record

	double precision 	:: delta_t           !Size of timestep for each computational step
	double precision 	:: elapsedtime       !Total time elapsed including restarts
	double precision 	:: rneighbr2         !Square of rcuttoff+delta_rneighbr
	double precision 	:: delta_rneighbr    !Radius used for neighbour list construction
	double precision 	:: rd                !Radius used for radial distribution function

	!Store surface bins of processes subdomain for outputs over periodic boundaries
	integer									:: nsurfacebins     !Number of surface bins
	integer,allocatable,dimension(:,:)	    :: surfacebins		!Surface Bins

	integer		 							:: nhalobins !Number of halo bins
	integer,allocatable,dimension(:,:)	    :: halobins		!halo Bins

	!Number and size of unit used for initial setup of molecules (i.e. FCC unit)
	integer,          dimension(3)		    	:: initialnunits
	double precision, dimension(:), allocatable :: initialunitsize

	!Size of global computational domain and domain per processor
	double precision, dimension(3)		    	:: globaldomain
	double precision, dimension(:), allocatable :: domain, halfdomain

	!Number and size of cells used for domain subdivision into cells of size ~rcutoff
	integer,          dimension(:), allocatable :: ncells
	double precision, dimension(:), allocatable :: cellsidelength, halfcellsidelength

	!Setup seed for random number generation
	integer,          dimension(:), allocatable :: seed 

	! Block/Process ID
	integer irank, iroot, ierr
	integer irankx, iranky, irankz
	integer iblock, jblock, kblock
	integer, allocatable :: ibmin(:), ibmax(:), ibmino(:), ibmaxo(:), &
							jbmin(:), jbmax(:), jbmino(:), jbmaxo(:), &
							kbmin(:), kbmax(:), kbmino(:), kbmaxo(:)

 	!Set number of halos
	integer, parameter :: nh = 1

	!Number of cells per processor
	integer ncellxl, ncellyl, ncellzl    !Inc Halos
	integer nicellxl, nicellyl, nicellzl !Inner only

	!Directory that holds input/output files, useful in coupling mode
 	character(len=128) :: prefix_dir = "./"	

end module computational_constants_MD

!-------------------------------------------------------------------------------------
!----------------------------------Polymer--------------------------------------------
module polymer_info_MD

	integer				:: etevtcf_outflag			!End-to-end time correlation function output flag
	integer				:: etevtcf_iter0			!Iteration from which to begin averaging
	integer             :: nchains                  !Number of FENE chains in domain
	integer				:: nmonomers				!Number of LJ beads per FENE chain
	integer             :: r_gyration_outflag       !Radius of gyration outflag
	integer             :: r_gyration_iter0         !Iteration at which to start recording R_g
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
		integer, allocatable :: bondflag(:)         !Integer array of length nmonomers describing bonds
	end type monomer_info

	type(monomer_info), dimension(:), allocatable :: monomer
	!eg. to access chainID of mol 23, call monomer(23)%chainID

	double precision, dimension(:,:), allocatable :: etev
	double precision, dimension(:,:), allocatable :: etev_0 !End-to-end vectors for polymer chains at iter=etevtcf_iter0 (size of np,only stored for leftmost chain mols) 

end module polymer_info_MD

!-------------------------------------------------------------------------------------
!-------------------------------Shearing BCs------------------------------------------
module shear_info_MD

	integer 			:: define_shear_as
	integer				:: shear_iter0
	integer				:: wrap_integer
	integer				:: shear_plane
	integer				:: shear_direction
	integer				:: shear_remainingplane
	double precision 	:: shear_velocity
	double precision 	:: shear_rate
	double precision	:: shear_distance
	double precision	:: shear_time

	integer, dimension(:), allocatable	:: mol_wrap_integer		

end module shear_info_MD

!-------------------------------------------------------------------------------------
!----------------------------------Arrays---------------------------------------------

module arrays_MD

	integer,          dimension(:),   allocatable, target	:: tag       	!Molecular Tags
	integer, 	  	  dimension(:,:), allocatable, target	:: &
		fix, &      				!Fixed molecules
		thermostat					!Thermostatted molecules
	double precision, dimension(:),   allocatable, target 	:: &
		potenergymol, 		&		!Potential energy of each molecule
		potenergymol_LJ, 	&		!LJ Potential energy of each molecule
		potenergymol_FENE,	&		!FENE Potential energy of each molecule
		virialmol,			&		!Virial of each molecule
		vmagnitude, 		&
		recvbuffer
	double precision, dimension(:,:),   allocatable			:: &
		rtrue, 		&      			!Positions with no period BC
		rinitial, 	&
		rijsum, 	&				!Sum of all molecular rij values
		theta, 		&
		aD,aR
	double precision, dimension(:,:),   allocatable, target 	:: &
		r, 		&        		  	!Positions
		v, 		&        		  	!Velocity
		a, 		&					!Accelerations
		aold,	&					!TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
		slidev						!Speed for sliding molecules

end module arrays_MD

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
!----------------------------------Linked List----------------------------------------

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

!-------------------------------------------------------------------------------------
!---------------------------------Calculated Properties-------------------------------

module calculated_properties_MD

	integer		 							:: nshells			    !Number of shells used to store radial distribution
	integer		 							:: nplanes			    !Number of planes used for MOP
	integer,dimension(3) 					:: nbins, globalnbins	!Number of groups or bins to store frequency of molecular velocity
	integer,dimension(:), allocatable 		:: vfd_bin        		!Array to keep tally of molecular velocity distribution
	integer,dimension(:), allocatable 		:: shell          		!Array to keep tally of radial distribution
	integer,dimension(:), allocatable 		:: slice_mass	    	!Array to keep tally of molecules in slice
	integer,dimension(:,:,:), allocatable 	:: slice_massbin 		!Recorded molecules in a bin
	integer,dimension(:,:,:), allocatable	:: volume_mass			!Mass in a control volume at time t
	integer,dimension(:,:,:,:), allocatable	:: mass_flux  			!Flow of mass over a control volume surface

	double precision :: 	&
		binsize,			&		!Size of each bin
		delta_r	,			&		!Size of each shell
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
		RDF,				&		!Array to keep tally of radial distribution
		normalisedvfd_bin,	&		!Bin normalised so sum of all bins is one
		diffusion,			&		!Diffusion of molecules
		meandiffusion,		&		!Time averaged diffusion of molecules
		Pxycorrel,			&     	!Sum of correlations of Pxy and Pxyzero
		slice_temperature,	&		!Temperature in a domain slice
		Pxyv_plane 	 				!Energy on plane for MOP

	double precision, dimension(:,:), allocatable 	:: & 
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
		slice_momentumbin,	&		!Mean velocity in a bin
		volume_momentum,	& 		!Momentum in a control volume at time t
		energy_flux,		&		!Flow of energy over a control volume surface
		Pxyvface					!Stress tensor on bin face

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
			bin = ceiling((r(n,ixyz)+halfdomain(ixyz))/binsize)
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
			bin = ceiling((r(n,ixyz)+halfdomain(ixyz))/binsize)
			if (bin.lt.1) 				bin = 1
			if (bin.gt.nbins(ixyz)) 	bin = nbins(ixyz)
			select case (integration_algorithm)
			case (leap_frog_verlet) 
				get_velo_slices(bin,:) = get_velo_slices(bin,:) + v(n,:) + 0.5d0*a(n,:)*delta_t
			case (velocity_verlet)
				get_velo_slices(bin,:) = get_velo_slices(bin,:) + v(n,:)
			end select
		end do

	end function get_velo_slices

	double precision function get_temperature_PUT()
	use computational_constants_MD, only: domain, halfdomain,delta_t
	use physical_constants_MD, only: np,nd,globalnp
	use arrays_MD, only: r,v,a
	use shear_info_MD, only: shear_plane
	implicit none
	
		integer	:: slicebin,n	
		integer, dimension(:), allocatable 	:: m_slice
		double precision :: pec_v2sum
		double precision, dimension(nd) 	:: slicebinsize, vel
		double precision, dimension(:,:), allocatable :: v_slice,v_avg
		
		slicebinsize(:) = domain(:)/nbins(:)				! Get bin size
		pec_v2sum 		= 0.d0								! Initialise
		
		!todo reevaluate
		allocate(m_slice(nbins(shear_plane)))	
		allocate(v_slice(nbins(shear_plane),nd))
		allocate(v_avg(nbins(shear_plane),nd))
		m_slice = get_mass_slices(shear_plane)				! Get total mass in all slices
		v_slice = get_velo_slices(shear_plane)				! Get total velocity in all slices
		do slicebin=1,nbins(shear_plane)
			v_avg(slicebin,:) = v_slice(slicebin,:)/m_slice(slicebin)
		end do

		do n=1,np
			slicebin 			= ceiling((r(n,shear_plane)+halfdomain(shear_plane))/slicebinsize(shear_plane))
			if (slicebin > nbins(shear_plane)) slicebin = nbins(shear_plane)	! Prevent out-of-range values
			if (slicebin < 1) slicebin = 1										! Prevent out-of-range values
			vel(:) 			 	= v(n,:) - v_avg(slicebin,:) - 0.5d0*a(n,:)*delta_t
			pec_v2sum 			= pec_v2sum + dot_product(vel,vel)
		end do

		call globalSum(pec_v2sum)
	
		get_temperature_PUT = pec_v2sum / real(nd*globalnp,kind(0.d0))

		deallocate(m_slice)
		deallocate(v_slice)
		deallocate(v_avg)

	end function get_temperature_PUT

end module calculated_properties_MD
