!--------------------------------------------------------------------------------
!
!                              Read Input files
! Reads input file
!
!--------------------------------------------------------------------------------

module module_read_input

	use interfaces
	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use polymer_info_MD
	use concentric_cylinders
	use shear_info_MD
    use calculated_properties_MD
    use boundary_MD

	! These variables are put here as a test -- The plan is to put EVERYTHING here
	! that could possibily be read in from the input file and proctect it! 
	! It is then called throughout the code by using this module...
	integer			 :: COUETTE_slidewall,  COUETTE_ixyz
	real(kind(0.d0)) :: COUETTE_t, COUETTE_Re, COUETTE_Uwall, COUETTE_H
	real(kind(0.d0)) :: initial_u, initial_v, initial_w

end module module_read_input
!----------------------------------------------------------------------------------

subroutine setup_read_input
	use module_read_input
	use librarymod, only :locate, linspace, find3factors
	implicit none

	logical					:: found_in_input, error, empty
	integer 				:: ios, ixyz, n, Nvmd_interval_size
	double precision,dimension(1000) :: temp
    character(256)          :: str

	! Open input file
	open(1,file=input_file)

	! ########################################################################
	! ## NEWPAGE - SYSTEM SETUP
	! ########################################################################
	call locate(1,'NEWPAGE_SYSTEM_SETUP',.false.,found_in_input)
	if (found_in_input) then
		!read(1,*) newpage
		print*, "The keyword OUTPUT does nothing, "
		print*, "it is used to denote start of output section flowmol_inputs" 
	endif


	! #########################################################################
	! # Ensemble selector
	! # [1] 0 - NVE
	! # [1] 1 - NVT (Nosé-Hoover thermostat)
	! # [1] 2 - NVT (Gaussian iso-kinetic thermostat) - only availble with VV
	! # [1] 3 - NVT (Profile unbiased Nosé-Hoover thermostat)
	! # [1] 4 - NVT (Pairwise additive Nosé-Hoover thermostat by Allen & Schmid)
	! # [1] 5 - NVT (DPD thermostat by Soddemann)
	! # [1] 6 - Tagged Move System
	! # Dynamically update tags based on location of molecules
	! # [2] 0 - Off
	! # [2] 1 - On
	! # -----------------------------------------------------------------------
	call locate(1,'ENSEMBLE',.true.)
	read(1,*) ensemble
	if (ensemble .eq. 6) then
		read(1,*,iostat=ios) dynamically_update_tags
		if (ios .ne. 0) dynamically_update_tags  = .false.
	endif

	! #########################################################################
	! # Number density of the particles in the system:
	! # [1] float - Density (float)
	! # -----------------------------------------------------------------------
	call locate(1,'DENSITY',.true.)
	read(1,*) density

	! #########################################################################
	! # Cut-off distance for particle interaction.
	! # Longer gives a slower calculation, notable values inclue:
	! # WCA cutofff is 2^(1/6) = 1.12246204830937
	! # Typical value is 2.5
	! # For surface tension which matches experiements ~4.5
	! # [1] float - Cutoff distance (float)
	! # -----------------------------------------------------------------------
	call locate(1,'RCUTOFF',.true.)
	read(1,*) rcutoff

	! #########################################################################
	! # Number of unit cells in the x,y,z directions.
	! # For initialunitsize "a" (1 cell size)
	! #  __________
	! #	|  o     o |
	! #	|     o    |  a/2 (distance between molcules)	
	! #	|  o     o |
	! #	|__________|
	! #
	! #   [1] int - number of FCC units in x
	! #   [2] int - number of FCC units in y 
	! #   [3] int - number of FCC units in z   
	! # -----------------------------------------------------------------------
	call locate(1,'INITIALNUNITS',.true.)
	read(1,*) initialnunits(1)		!x dimension split into number of cells
	read(1,*) initialnunits(2)		!y dimension split into number of cells
	read(1,*) initialnunits(3)		!z dimension split into number of cells

	! #########################################################################
	! # Initial temperature to use for the random molecule velocities picked
	! # from a Maxwell Boltzmann distribution:
	! # [1] float - Temperature (float)
	! # -----------------------------------------------------------------------
	call locate(1,'INPUTTEMPERATURE',.true.)
	read(1,*) inputtemperature

	! #########################################################################
	! # Thermostat target temperature (if not specified, uses input temperature)
	! # Can be specified as multiple arguments if using top/bottom walls
	! # where first argument is then bottom wall and second top
	! # [1] float - Temperature (float)
	! # [2] float - Temperature of top wall (float) 
	! # -----------------------------------------------------------------------	
	call locate(1,'THERMOSTATTEMPERATURE',.false.,found_in_input)
	if (found_in_input) then
        nthermo = 1
		read(1,*,iostat=ios) temp(1)
		if (ios .ne. 0) then
            call error_abort( "Error -- specify at least one value for THERMOSTATTEMPERATURE")
		endif
        do n = 2,1000
            read(1,*,iostat=ios) temp(n)
		    if (ios .ne. 0) then
                exit
            else
                nthermo = nthermo + 1
            endif
        enddo
        if (nthermo .lt. 1) then
            call error_abort( "Error -- specify at least one value for THERMOSTATTEMPERATURE")
        endif
        allocate(thermostattemperature(nthermo))
        thermostattemperature = temp(1:nthermo)
    else
        nthermo = 1
        allocate(thermostattemperature(nthermo))
        thermostattemperature = inputtemperature
    endif


	! #########################################################################
	! # Initial configuration
	! # [1] 0 - FCC Lattice
	! # [1] 1 - Special case, must be followed by string on next line. 
	! # [1] 2 - Configuration file, e.g. pdb or lammps input. Under developement.
	! #
	! #        This string must be lower case, so the same string in capitals is
	! #        reserved as an input flag for information about the special
	! #        case. Options:
	! #
	! # [2] solid_liquid - tethered walls of different density to fluid
	! # [2] dense_fene - connect monomers on an FCC lattice with a specified 
	! #                  chain length and FENE potential parameters. 
	! # [2] sparse_fene - connect monomers on a cubic lattice
	! #                   separated by FENE equilibrium distance
	! # [2] single_fene - A single fene chain
	! # [2] fene_solution - fene chains in an explicit solvent
	! # [2] droplet2D - A cylindrical droplet on a surface
	! # [2] droplet3D - A sherical droplet on a surface
	! # [2] 2phase - A two phase coexisitence
	! # [2] 2phase_LJ - A two phase coexisitence using Lennard Jones
	! # [2] bubble - A bubble in a liquid
	! # [2] film - A solid with a liquid film and gas region above
	! # [2] polymer_brush - A solid with polymer chains attached
	! # [2] 2phase_surfactant_solution - A fluid with SAFT gamma Mie surfactants
	! # [2] 2phase_surfactant_atsurface - A fluid with SAFT gamma Mie surfactants at surface
	! # [2] concentric_cylinders - build concentric cylinders from an FCC lattice
	! #                 and melt them between specular walls
	! # [2] fill_cylinders - tether concentric cylinders from previous simulation
	! #                  to their initial sites, and fill them with fluid
	! # [2] fill_cylinders_fene_solution - concentric cylinders with FENE 
	! # [2] rotate_cylinders - restart from filled cylinder and rotate 
	! # 						inner cylinder with specified angular velocity 
	! # -----------------------------------------------------------------------
	call locate(1,'INITIAL_CONFIG_FLAG',.true.)
	read(1,*) initial_config_flag
	if (initial_config_flag .eq. 0) then
		potential_flag = 0
	endif
	if (initial_config_flag .eq. 1)	then

		read(1,*) config_special_case
        print*, "initial config special case ", trim(config_special_case)

		select case (trim(config_special_case))
		case('sparse_fene')	

			potential_flag = 1	
			rcutoff = 2.d0**(1.d0/6.d0)

			! #########################################################################
			! # Sparse FENE special case info: designed for cases where the density
			! # is lower than the density at which the average separation between
			! # monomers on an FCC lattice is greater than the FENE maximum bond 
			! # elongation.
			! #
			! #	[1] int - number of LJ beads for each polymer chain (nmonomers)
			! #	[2] float - spring constant, k_c
			! #	[3] float - maximum spring elongation
			! # [4] float - domain length in x
			! # [5] float - domain length in y
			! # [6] float - domain length in z
			! # [7] int  - nchains
			! # -----------------------------------------------------------------------
			call locate(1,'SPARSE_FENE',.true.)
			read(1,*) nmonomers
			read(1,*) k_c
			read(1,*) r_0
			read(1,*) globaldomain(1)
			read(1,*) globaldomain(2)
			read(1,*) globaldomain(3)
			read(1,*) nchains

			density = nmonomers * nchains / product(globaldomain(1:3))	

		case('dense_fene')	

			potential_flag = 1	
			rcutoff = 2.d0**(1.d0/6.d0)

			! #########################################################################
			! # Dense FENE special case information
			! #
			! #	[1] int - number of LJ beads for each polymer chain (nmonomers)
			! #	[2] float - spring constant, k_c
			! #	[3] float - maximum spring elongation
			! # [4] float - monomer density
			! # [5] int - FCC (4 monomer) units in x
			! # [6] int - FCC (4 monomer) units in y
			! # [7] int - FCC (4 monomer) units in z
			! # -----------------------------------------------------------------------
			call locate(1,'DENSE_FENE',.true.)
			read(1,*) nmonomers
			read(1,*) k_c
			read(1,*) r_0
			read(1,*) density
			read(1,*) initialnunits(1)
			read(1,*) initialnunits(2)
			read(1,*) initialnunits(3)	

        case('fene_solution','single_fene')

            potential_flag = 1
            rcutoff = 2.d0**(1.d0/6.d0)

			! #########################################################################
			! # FENE_fill special case information
			! #
			! #	[1] int - number of LJ beads for each polymer chain (nmonomers)
			! #	[2] float - spring constant, k_c
			! #	[3] float - maximum spring elongation
			! # [4] float - monomer density
			! # [5] float - target concentration
			! # -----------------------------------------------------------------------
            call locate(1,'FENE_SOLUTION',.true.)
			read(1,*) nmonomers
			read(1,*) k_c
			read(1,*) r_0
			read(1,*) density
            read(1,*) targetconc

			! #########################################################################
			! # Solvent information:
			! #	[1] 0 - All bead interactions equal (plus FENE)
			! #	[1]	1 - Solvent of variable quality (Soddemann)
			! #	[2] float - solvent energy parameter eps_pp (polymer-polymer)
			! #	[3] float - solvent energy parameter eps_ps (polymer-solvent)
			! #	[4] float - solvent energy parameter eps_ss (solvent-solvent)
			! # -----------------------------------------------------------------------
			call locate(1,'SOLVENT_INFO',.true.)
			read(1,*) solvent_flag
			select case(solvent_flag)
			case(0)
			case(1)
				read(1,*) eps_pp
				read(1,*) eps_ps
				read(1,*) eps_ss
			case default
				call error_abort('Unrecognised solvent flag!')
			end select
            
		case('solid_liquid')

			potential_flag = 0

			! #########################################################################
			! # Number density of the untethered or liquid particles in the system
			! # Use with solid_liquid config flag
			! # [1] float - Density (float)
			! # -----------------------------------------------------------------------
			call locate(1,'LIQUIDDENSITY',.true.)
			read(1,*) liquid_density

		case('droplet2D','droplet3D','2phase','2phase_LJ', 'bubble', 'film')

            potential_flag = 0

			! #########################################################################
			! # Number density of the liquid part of the domain obtained by	
			! # randomly removing particles until target density is obtained:
			! # [1] float - Density (float)
			! # -----------------------------------------------------------------------
			call locate(1,'LIQUIDDENSITY',.true.) 
    	    read(1,*) liquid_density

			! #########################################################################
			! # Number density of the gas part of the domain obtained by	
			! # randomly removing particles until target density is obtained:
			! # [1] float - Density (float)
			! # -----------------------------------------------------------------------
			call locate(1,'GASDENSITY',.true.) 
    	    read(1,*) gas_density

            if (config_special_case .eq. 'droplet2D') then
			    call locate(1,'DROPHEIGHT',.false.,found_in_input)
	            if (found_in_input) then
                    read(1,*) dropletH
                endif
			    call locate(1,'DROPHLRATIO',.false.,found_in_input)
	            if (found_in_input) then
                    read(1,*) dropletHLratio
                endif
            endif

            if (config_special_case .eq. '2phase' .or. &
                config_special_case .eq. '2phase_LJ' .or. &
                config_special_case .eq. 'film') then
				! #########################################################################
				! # Name of output from FEA code 
				! # [1] str - Name of input file in correct format (see flowmol sourcecode)
				! # -----------------------------------------------------------------------
			    call locate(1,'FEA_FILENAME',.false.,found_in_input) 
	            if (found_in_input) then
                    Twophase_from_file = .true.
                    read(1,*) FEA_filename
                endif
            endif

			! #########################################################################
			! # Fraction of the domain filled with liquid (from 0 to 1)
			! # [1] float - Fraction (float)
			! # [2] 1 - x direction
			! # [2] 2 - y direction
			! # [2] 3 - z direction
			! # -----------------------------------------------------------------------
			call locate(1,'LIQUID_FRACTION',.false.,found_in_input) 
			if (found_in_input) then
				read(1,*) lg_fract
				read(1,*,iostat=ios) lg_direction
				if (ios .ne. 0) then
					print*, "Default direction not given for LIQUID_FRACTION, assuming x"
					lg_direction = 1
				endif
			endif

            if (config_special_case .eq. "bubble") then
				! #########################################################################
				! # Location and radius of bubble
				! # [1] float - Radius of bubble
				! # [2] float - centre in x direction
				! # [3] float - centre in y direction
				! # [4] float - centre in z direction
				! # -----------------------------------------------------------------------
                call locate(1,'BUBBLERADIUS',.true.)
                read(1,*) rbubble
                read(1,*) rcentre(1)
                read(1,*) rcentre(2)
                read(1,*) rcentre(3)
            endif

		case('concentric_cylinders')
			
			potential_flag = 0
			rcutoff = 2.d0**(1.d0/6.d0)

			call locate(1,'CONCENTRIC_CYLINDERS',.true.)
			read(1,*) cyl_density
			read(1,*) cyl_units_oo 
			read(1,*) cyl_units_io 
			read(1,*) cyl_units_oi 
			read(1,*) cyl_units_ii 
			read(1,*) cyl_units_z 
		
			initialnunits(1) = cyl_units_oo	+ 1 !To leave a gap
			initialnunits(2) = cyl_units_oo	+ 1 !To leave a gap
			initialnunits(3) = cyl_units_z

		case('fill_cylinders')
			
			potential_flag = 0
            if (ensemble .ne. tag_move) then
                call error_abort("Error in setup_read_input -- special case "&
                     //"fill_cylinders needs ENSEMBLE=6 for tag_move_system")
            endif

		case('fill_cylinders_fene_solution')

			potential_flag = 1	    
            if (ensemble .ne. tag_move) then
                call error_abort("Error in setup_read_input -- special case "&
                //"fill_cylinders_fene_solution needs ENSEMBLE=6 for tag_move_system")
            endif
			rcutoff = 2.d0**(1.d0/6.d0)

			!#########################################################################
			!# FENE_fill special case information
			! #
			! #	[1] int - number of LJ beads for each polymer chain (nmonomers)
			! #	[2] float - spring constant, k_c
			! #	[3] float - maximum spring elongation
			! # [4] float - monomer density
			! # [5] float - target concentration
			!# -----------------------------------------------------------------------
            call locate(1,'FENE_SOLUTION',.true.)
			read(1,*) nmonomers
			read(1,*) k_c
			read(1,*) r_0
			read(1,*) density
            read(1,*) targetconc

		case('rotate_cylinders')
			
            if (ensemble .ne. tag_move) then
                call error_abort("Error in setup_read_input -- special case "&
                //"rotate_cylinders needs ENSEMBLE=6 for tag_move_system")
            endif

			call locate(1,'POTENTIAL_FLAG',.true.)
            read(1,*) potential_flag

			call locate(1,'ROTATE_CYLINDERS',.true.)
			read(1,*) omega_i
			read(1,*) omega_f 
			read(1,*) omega_ramplength

		case('polymer_brush')

            potential_flag = 1
            if (ensemble .ne. tag_move) then
                call error_abort("Error in setup_read_input -- special case "&
                //"polymer_brush needs ENSEMBLE=6 for tag_move_system")
            endif
            rcutoff = 2.d0**(1.d0/6.d0)

			! #########################################################################
			! # Polymer brush special case
			! # 
			! #  [1] int - N Number of LJ beads for each chain
			! #  [2] float - k Spring constant
			! #  [3] float - R_0 Maximum bond elongation
			! #  [4] float - Grafting density
			! # -----------------------------------------------------------------------
            call locate(1,'POLYMER_BRUSH',.true.)
            read(1,*) nmonomers
            read(1,*) k_c
            read(1,*) r_0
            read(1,*) grafting_density


			call locate(1,'SOLVENT_INFO',.true.)
			read(1,*) solvent_flag
			select case(solvent_flag)
			case(0)
			case(1)
				read(1,*) eps_pp
				read(1,*) eps_ps
				read(1,*) eps_ss
			case default
				call error_abort('Unrecognised solvent flag!')
			end select

            call locate(1,'LIQUIDDENSITY',.true.)
            read(1,*) liquid_density

		case('2phase_surfactant_solution', '2phase_surfactant_atsurface')

            !Specifiy more general potential than LJ
	        call locate(1,'MIE_POTENTIAL',.false.,found_in_input) 
	        if (found_in_input) then
                read(1,*) Mie_potential
!                print*, "Specifying MIE_POTENTIAL in the input is depreciated"
!                print*, "Must be specified as a parameter in the module file"
!                print*, "as this results in a 30% efficiency difference"
                read(1,*,iostat=ios) default_moltype
				if (ios .ne. 0) then
                    print*, "Default moltype not given -- assuming Argon (=1)"
                    default_moltype = 1
                endif
            else
                Mie_potential = 0
            endif

			call locate(1,'POTENTIAL_FLAG',.true.)
            read(1,*) potential_flag

            if (mie_potential .ne. 1) then
                call error_abort("2phase_surfactant_solution initial case "&
                //"used but Mie flag is off -- aborting")
            endif
            if (potential_flag .ne. 1) then
                print*, "2phase_surfactant_solution initial case used but "&
                //"potential flag is off -- switching on"
                potential_flag = 1
            endif
            if (ensemble .ne. tag_move) then
                call error_abort("Error in setup_read_input -- special case "&
                //"2phase_surfactant_solution needs ENSEMBLE=6 for tag_move_system")
            endif

            call locate(1,'2PHASE_SURFACTANT',.true.)
			read(1,*) nmonomers
            read(1,*) targetconc
            read(1,*) angular_potential
            !if (config_special_case .eq. '2phase_surfactant_atsurface') then
			read(1,*,iostat=ios) surface_surfactant_layer
			if (ios .ne. 0) surface_surfactant_layer = 4.d0
            !endif

			!print*, "2phase conc =", targetconc

            if (targetconc .gt. 1.d0) then
                call error_abort("ERROR in 2PHASE_SURFACTANT input -- "&
                //"targetconc must be between 0.0 and 1.0")
            end if

            call locate(1,'LIQUIDDENSITY',.true.)
            read(1,*) liquid_density
            if (density .lt. liquid_density) then
                call error_abort("ERROR in 2PHASE_SURFACTANT input -- "&
                //"DENSITY must be greater than LIQUIDDENSITY")
            end if
			call locate(1,'GASDENSITY',.true.) 
    	    read(1,*) gas_density
            if (liquid_density .lt. gas_density) then
                call error_abort("ERROR in 2PHASE_SURFACTANT input -- "&
                //"LIQUIDDENSITY must be greater than GASDENSITY")
            end if

		    call locate(1,'LIQUID_FRACTION',.false.,found_in_input) 
            if (found_in_input) then
                read(1,*) lg_fract
                read(1,*,iostat=ios) lg_direction
                if (ios .ne. 0) then
                    print*, "Default direction not given for LIQUID_FRACTION, assuming x"
                    lg_direction = 1
                endif
            endif

		case default
            print*, "Config special case string", trim(config_special_case)
			call error_abort("ERROR in setup_read_input -- Unrecognised special case string")
		end select

	endif
	if (initial_config_flag .eq. 2) then
		call error_abort("ERROR in setup_read_input -- Generic input file read in is not developed yet")
	endif 


	call locate(1,'INITIAL_VELOCITY_FLAG',.false.,found_in_input) 
	if (found_in_input) then
		read(1,*) initial_velocity_flag 
		if (initial_velocity_flag .ne. 0) then
	   		read(1,*) velocity_special_case	
	   		select case (trim(velocity_special_case	))
	   		case('debug')
				!Nothing must be read in
	   		case('taylor_green')
				!Nothing must be read in
	   		case('couette_analytical')
		   		read(1,*) COUETTE_t
		   		read(1,*) COUETTE_Re
		   		read(1,*) COUETTE_Uwall
		   		read(1,*) COUETTE_h
		   		read(1,*) COUETTE_slidewall
		   		read(1,*) COUETTE_ixyz
	   		case('constant')
		   		read(1,*) initial_u
		   		read(1,*) initial_v
		   		read(1,*) initial_w
	   		case('dns')
		   		read(1,*) DNS_filename
		   		read(1,*) DNS_ngx
		   		read(1,*) DNS_ngy
		   		read(1,*) DNS_ngz
	   		case default
	   			call error_abort('Unidentified initial velocities_special_case')	
	   		end select
		endif
	else
		initial_velocity_flag = 0
	endif

	! ########################################################################
	! ## NEWPAGE - COMPUTATIONAL PARAMETERS
	! ########################################################################
	call locate(1,'NEWPAGE_COMPUTATIONAL_PARAMETERS',.false.,found_in_input)
	if (found_in_input) then
		!read(1,*) newpage
		print*, "The keyword OUTPUT does nothing, "
		print*, "it is used to denote start of output section flowmol_inputs" 
	endif


	! #########################################################################
	! # Integration algorithm
	! # [1] 0 - Leap-frog Verlet
	! # [1] 1 - Velocity Verlet
	! # [1] 2 - Other
	! # -----------------------------------------------------------------------
	call locate(1,'INTEGRATION_ALGORITHM',.true.)
	read(1,*) integration_algorithm

	! #########################################################################
	! # Force list used
	! # [1] 0 - All Pairs
	! # [1] 1 - Cell List
	! # [1] 2 - Neighbour list with all interactions
	! # [1] 3 - Neighbour list using 3rd law optimisation (half interactions)
	! # -----------------------------------------------------------------------
	call locate(1,'FORCE_LIST',.true.)	!LJ or FENE potential
	read(1,*) force_list
	if (force_list .ne. 3 .and. ensemble .eq. 4) then
		call error_abort("Half int neighbour list only is compatible with pwa_terms_pwaNH thermostat")
	endif
	if (force_list .ne. 3 .and. ensemble .eq. 5) then 
		call error_abort("Half int neighbour list only is compatible with DPD thermostat")
	endif


	! #########################################################################
	! # Total number of timesteps:
	! # [1] int - Number of steps
	! # -----------------------------------------------------------------------
	call locate(1,'NSTEPS',.true.)
	read(1,*) Nsteps 		!Number of computational steps

	! #########################################################################
	! # Timestep:
	! # [1] float - timestep delta t (usually 0.005)
	! # -----------------------------------------------------------------------
	call locate(1,'DELTA_T',.true.)
	read(1,*) delta_t 		!Size of time step

	! #########################################################################
	! # Frequency to collect statistics in number of timesteps:
	! # [1] int - Frequency at which to record results
	! # -----------------------------------------------------------------------
	call locate(1,'TPLOT',.true.)
	read(1,*) tplot 

	! #########################################################################
	! # Number of steps before collecting statistics (doesn't appear to be used)
	! # 
	! # [1] int - Number of initialisation steps for simulation
	! # -----------------------------------------------------------------------
	call locate(1,'INITIALISE_STEPS',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) initialise_steps 
	else
		initialise_steps = 0
	endif

	! #########################################################################
	! # Extra distance to add to rcutoff to use in neighbourlist
	! # building, bigger gets more interactions but requires rebuild
	! # to occur less often. Can also be used to tweak cellsizes which
	! # in term define binsizes for averaging. A single value
	! # is used for the radius, while 3 values are used in each x,y and z-direction:
	! # [1] float - extra distance (r or x)
	! # [2] float - extra distance (y)
	! # [3] float - extra distance (z)
	! # -----------------------------------------------------------------------	
	call locate(1,'DELTA_RNEIGHBR',.true.) 
	read(1,*) delta_rneighbr(1)
	read(1,*,iostat=ios) delta_rneighbr(2)
    if (ios .ne. 0) then
        delta_rneighbr(2) = delta_rneighbr(1)
        delta_rneighbr(3) = delta_rneighbr(1)
    else
    	read(1,*,iostat=ios) delta_rneighbr(3)
        if (ios .ne. 0) then
            call error_abort("Error -- DELTA_RNEIGHBR should be one value or three for x, y and z")
        endif
    endif

	! #########################################################################
	! # rebuild CRITERIA flag
	! # [1] 0 - Use Rapaport delta_rneighbr to calculate rebuild with vmax each step
	! # [1] 1 - Use fixed frequency specifed on next line
	! # [1] 2 - Use exact displacement of all molecules to decide rebuild
	! # [1] 3 - mixed - don't rebuild before frequency specifed on next line then
	! #			 use delta_rneighbr to calculate rebuild (NOT AVAILABLE YET!!)
	! # [2] int - fixed_rebuild requency
	! # -----------------------------------------------------------------------
	call locate(1,'REBUILD_CRITERIA',.false.,found_in_input) 
	if (found_in_input) then
		read(1,*) rebuild_criteria   	!Choice of rebuild criteria
        if (rebuild_criteria .eq. 1) then
    		read(1,*) fixed_rebuild 		!Fixed rebuild frequency
        endif
	else
		rebuild_criteria = 0			!Rebuild uses neighbourcell
	endif

	! #########################################################################
	! # Frequency (in seconds) to write rescue snapshot 'final_state' file, 
	! # e.g. 5min=3600 or 6 hours=21600 (which is the defaul)
	! # [1] int - Seconds 
	! # [2] .false. - Save individual "interim_state" snapshots
	! # [2] .true. - Overwrite "final_state" (default)
	! # -----------------------------------------------------------------------
	call locate(1,'RESCUE_SNAPSHOT_FREQ',.false.,found_in_input) 
	if (found_in_input) then
		read(1,*) rescue_snapshot_freq 	!Rescue snapshot frequency in seconds
        read(1,*,iostat=ios) overwrite_rescue_snapshot
		if (ios .ne. 0) overwrite_rescue_snapshot = .true.
	else
		rescue_snapshot_freq = 21600	!Every 6 hours
        overwrite_rescue_snapshot = .true.
	endif

	! #########################################################################
	! # Sort molecules using a Riemann spacing filling curve which 
	! # take reorders molecular position, velocities, etc to maximise
	! # memory locality and increase cache hits. Shows some improvement
	! # 
	! # [1] 0 - Off
	! # [1] 1 - On
	! # [2] int - Frequency to sort in number of timesteps
	! # [3] int - Size of block used for molecules
	! # -----------------------------------------------------------------------
	call locate(1,'SORT_FLAG',.false.,found_in_input) 
	if (found_in_input) then
		read(1,*) sort_flag
		read(1,*) sort_freq
		read(1,*) sortblocksize
	else !Default to switched off at the moment for back compatibility
		sort_flag = 0
		sort_freq = 0
		sortblocksize = 0
	endif

	! #########################################################################
	! # Global numbering system, creates another array to keep track of the
	! # numbering of the molecules
	! # [1] 0 - Off
	! # [1] 1 - On
	! # -----------------------------------------------------------------------
	call locate(1,'GLOBAL_NUMBERING',.false.,found_in_input) 
	if (found_in_input) then
		read(1,*) global_numbering  	!Include global molecular numbering
        if (potential_flag .eq. 1 .and. global_numbering .eq. 1) then
            print*, "POTENTIAL_FLAG is on so GLOBAL_NUMBERING not needed, use mononmer(n)%glob_no instead"
            global_numbering = 0
        endif
        if (global_numbering .eq. 1 .and. sort_flag .ne. 0 ) then
            call error_abort("Global number does not work with sort flag on")
        endif
    else
        global_numbering = 0
    endif


	! #########################################################################
	! # Random number generator seed system, needs to be two numbers 
	! # and defaults to 1 and 2 if not specified. Set both to the same number
	! # to trigger a random seed based on the current date and time.
	! # [1] int - Seed number 1
	! # [2] int - Seed number 2
	! # -----------------------------------------------------------------------
	call locate(1,'SEED',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) seed(1) 	!Random number seed value 1
		read(1,*) seed(2) 	!Random number seed value 2
	else
		seed(1) = 1		!Fixed default seed for repeatability
		seed(2) = 2		!Fixed default seed for repeatability
	endif

	! #########################################################################
	! # Specifiy more general potential than LJ, including powers 
	! # 12 and 6 with lambda_a and lambda_r to allow tuning to more 
	! # general molecular models with mol_type keeping track of the 
	! # parameter which chooses a pre-defined molecule type and
	! # allows varying wetting potential/wall interaction eij
	! # [1] 0 - Lennard Jones Potential
	! # [1] 1 - Mie Potential
	! # [2] int - Default value of Mie molecule type
	! #-----------------------------------------------------------------------
	call locate(1,'MIE_POTENTIAL',.false.,found_in_input) 
	if (found_in_input) then
        read(1,*) Mie_potential
        read(1,*,iostat=ios) default_moltype
		if (ios .ne. 0) then
            print*, "Default moltype not given -- assuming Argon (=1)"
            default_moltype = 1
        endif
        !IF Mie potential, check for EIJ wall
        if (Mie_potential .ne. 0) then
            call locate(1,'EIJ_WALL',.false.,found_in_input)
            if (found_in_input) then
                read(1,*,iostat=ios) eij_wall(1)
        		if (ios .ne. 0) then
                    call error_abort('Input Error -- EIJ_WALL no specified value') 
                endif
                read(1,*,iostat=ios) eij_wall(2)
        		if (ios .ne. 0) eij_wall(2) = eij_wall(1)
            else
                eij_wall = 1.d0
            endif
        endif
    else
        Mie_potential = 0
    endif

    call locate(1,'PROCESSORS',.false.,found_in_input)
	if (found_in_input) then
	    read(1,*) npx
   		read(1,*) npy
	    read(1,*) npz
	    !check if npx*npy*npz=nproc
	    if (npx * npy * npz .ne. nproc ) then
			print*, npx, npy, npz , nproc
			call error_abort(' Wrong specification for processor topology, nproc not equal to npx*npy*npz')
	    endif
	else
		!Assign (using prime factors) to each dimension if not specified
		call find3factors(nproc,npx,npy,npz)
		print*, 'WARNING - Number of processors not specified - Arbitrarily assigned as follows:'
		print*, 'npx = ', npx, 'npy = ', npy, 'npz = ', npz

	endif


	! ########################################################################
	! ## NEWPAGE - BOUNDARY CONDITIONS
	! ########################################################################
	call locate(1,'NEWPAGE_BOUNDARY_CONDITIONS',.false.,found_in_input)
	if (found_in_input) then
		!read(1,*) newpage
		print*, "The keyword OUTPUT does nothing, "
		print*, "it is used to denote start of output section flowmol_inputs" 
	endif


	!Flags to determine if periodic boundaries are on or shearing Lees Edwards
	call locate(1,'PERIODIC',.true.)
	read(1,*) periodic(1)
	read(1,*) periodic(2)
	read(1,*) periodic(3)

	if (any(periodic.eq.0)) then

		bforce_flag(:) = 0
		bforce_dxyz(:) = 0.0

		! #########################################################################
		! # Boundary force flag and distances, formatted:
		! #   bforce_flag(1) (x direction, bottom)           (integer)
		! #   bforce_flag(2) (x direction, top)              (integer)
		! #   bforce_flag(3) (y direction, bottom)           (integer)
		! #   bforce_flag(4) (y direction, top)              (integer)
		! #   bforce_flag(5) (z direction, bottom)           (integer)
		! #   bforce_flag(6) (z direction, top)              (integer)
		! #   bforce_dxyz(1) (x distance, bottom)            (real)
		! #   bforce_dxyz(2) (x distance, top)               (real)
		! #   bforce_dxyz(3) (y distance, bottom)            (real)
		! #   bforce_dxyz(4) (y distance, top)               (real)
		! #   bforce_dxyz(5) (z distance, bottom)            (real)
		! #   bforce_dxyz(6) (z distance, top)               (real)
		! #
		! # If periodic boundaries in x, y or z directions are turned off, the
		! # option to apply a boundary force to keep molecules within the domain
		! # is provided via the bforce flag. Choose from the following in each
		! # direction:
		! #
		! #    0 = off 
		! #    1 = O'Connell & Thompson     1995  (Phys. Rev. E 52, 5792)
		! #    2 = Nie, Chen, E & Robbins   2003  (JFM 500 pp. 55-64)
		! #    3 = Flekkøy, Wagner & Feder  2000  (Europhys. Lett. 52 3, p271)
		! #    4 = boundary force pdf input file ./bforce.input 
		! #
		! # bforce_dxyz defines the distance from each edge of the domain beyond 
		! # which any particle will no longer feel the effect of the boundary 
		! # force. Low values correspond to hard walls, while larger values will
		! # make the boundaries softer.
		! # -----------------------------------------------------------------------
		call locate(1,'BFORCE',.false.,found_in_input)
		if (found_in_input) then

			read(1,*) bforce_flag(1)			
			read(1,*) bforce_flag(2)			
			read(1,*) bforce_flag(3)
			read(1,*) bforce_flag(4)			
			read(1,*) bforce_flag(5)			
			read(1,*) bforce_flag(6)
			read(1,*) bforce_dxyz(1)
			read(1,*) bforce_dxyz(2)
			read(1,*) bforce_dxyz(3)
			read(1,*) bforce_dxyz(4)
			read(1,*) bforce_dxyz(5)
			read(1,*) bforce_dxyz(6)

			! Correct any bforce_flags if periodic boundaries on
			if (periodic(1).ne.0) then
!				if (irank.eq.iroot) print*, 'Warning, resetting bforce_flag(1:2)' , &
!				' because periodic boundaries are on in the x-direction'
				bforce_flag(1:2) = 0
			end if
			if (periodic(2).ne.0) then 
!				if (irank.eq.iroot) print*, 'Warning, resetting bforce_flag(3:4)' , &
!				' because periodic boundaries are on in the y-direction'
				bforce_flag(3:4) = 0
			end if
			if (periodic(3).ne.0) then 
!				if (irank.eq.iroot) print*, 'Warning, resetting bforce_flag(5:6)' , &
!				' because periodic boundaries are on in the z-direction'
				bforce_flag(5:6) = 0
			end if
#if __INTEL_COMPILER > 1200
            if (any(bforce_flag.eq.bforce_pdf_input)) then
                call load_bforce_pdf
            end if
#endif            
            do n=1,6
                if (bforce_flag(n) .eq. bforce_pdf_input) then
                    bforce_dxyz(n) = rcutoff
                end if
            end do

            if (any(bforce_flag .eq. substrate_force)) then
                call locate(1,'EIJ_WALL',.false.,found_in_input)
                if (found_in_input) then
                    read(1,*,iostat=ios) eij_wall(1)
            		if (ios .ne. 0) then
                        call error_abort('Input Error -- EIJ_WALL no specified value') 
                    endif
                    read(1,*,iostat=ios) eij_wall(2)
            		if (ios .ne. 0) eij_wall(2) = eij_wall(1)
                else
                    eij_wall = 1.d0
                endif
            endif

		end if

		! #########################################################################
		! # Open boundary flags, formatted:
		! #   open_boundary(1) (x direction, bottom)           (integer)
		! #   open_boundary(2) (x direction, top)              (integer)
		! #   open_boundary(3) (y direction, bottom)           (integer)
		! #   open_boundary(4) (y direction, top)              (integer)
		! #   open_boundary(5) (z direction, bottom)           (integer)
		! #   open_boundary(6) (z direction, top)              (integer)
		! #
		! # If periodic boundaries in x, y or z directions are turned off, the
		! # option to catch molecules leaving the domain is provided via the 
		! # open_boundary flag. Choose from the following in each
		! # direction:
		! #
		! #    0 = off 
		! #    1 = on 
		! # 
		! # Molecules due to be sent to MPI_PROC_NULL across the open boundary are 
		! # caught and counted in the messenger. This number can then be used to
		! # reinsert the same number of molecules after rebuild.
		! # -----------------------------------------------------------------------
		open_boundary(:) = 0
		call locate(1,'OPEN_BOUNDARY',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) open_boundary(1) 
			read(1,*) open_boundary(2) 
			read(1,*) open_boundary(3) 
			read(1,*) open_boundary(4) 
			read(1,*) open_boundary(5) 
			read(1,*) open_boundary(6)
		end if


		error = .false.
		do ixyz=1,3
			if (periodic(ixyz).eq.1) then
				if (open_boundary(2*ixyz - 1).ne.0) then
					print'(a,i6,a,i6,a)',  'Open boundary ',  2*ixyz - 1, ' is on but periodic Boundary ', ixyz, ' is also stil on '
					error = .true.
				endif
				if (open_boundary(2*ixyz).ne.0) then
					print'(a,i6,a,i6,a)',  'Open boundary ',  2*ixyz , ' is on but periodic Boundary ', ixyz, ' is also stil on '
					error = .true.
				endif
			endif
		enddo
		if (error) call error_abort("ERROR - Periodic Boundaries must be turned off for OPEN_BOUNDARY to be on ")

		call locate(1,'MEASURE_BFORCE_PDF',.false.,found_in_input) 
		if (found_in_input) then
			read(1,*) bforce_pdf_measure
			if (bforce_pdf_measure .ne. 0 .and. any(open_boundary.ne.0)) then
				call error_abort('Cannot measure the PDF of the boundary '//&
				'force if any boundaries are open or forced. Aborting.') 
			end if
			read(1,*) bforce_pdf_nsubcells
			read(1,*) bforce_pdf_nbins
			read(1,*) bforce_pdf_min
			read(1,*) bforce_pdf_max
			read(1,*) bforce_pdf_Nave
		else
			bforce_pdf_measure = 0 
		endif

	end if



	! #########################################################################
	! # Apply external force to region of space
	! # 0 - off
	! # 1 - apply to all molecules, requires
	! # 		- F_ext direction x=1,y=2,z=3
	! #		- F_ext magnitude
	! # 2 - apply only to specified region
	! # 		- F_ext direction x=1,y=2,z=3
	! #		- F_ext magnitude
	! #		- xmin - minimum x coordinate in global system 
	! #		- xmax - maximum x coordinate in global system 
	! #		- ymin - minimum y coordinate in global system 
	! #		- ymax - maximum y coordinate in global system 
	! #		- zmin - minimum z coordinate in global system 
	! #		- zmax - maximum z coordinate in global system 
	! #		NOTE : min/max of globaldomain and specified extent is used
	! # 			   so simply specifiy large numbers if you want
	! #			   region to extend to edge of globaldomain
	! # -----------------------------------------------------------------------
	call locate(1,'EXTERNAL_FORCE',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) external_force_flag
		if (external_force_flag .ge. 1) then
			read(1,*) F_ext_ixyz
			read(1,*) F_ext

			if (external_force_flag .eq. 2) then
				read(1,*) F_ext_limits(1)
				read(1,*) F_ext_limits(2)
				read(1,*) F_ext_limits(3)
				read(1,*) F_ext_limits(4)
				read(1,*) F_ext_limits(5)
				read(1,*) F_ext_limits(6)
			endif
        else
            external_force_flag = 0
            F_ext_ixyz = 1
            F_ext = 0.d0
		endif
	else
		external_force_flag = 0
	endif


	! #########################################################################
	! #
	! # Apply a CV based force with a number of options
	! # Apply force, either:
	! #   NON-COUPLED 
	! #      -1 - Debug, apply nothing
	! #	    0 - Zero force (MD values taken off and no force applied)
	! #	    1 - Sinusoidal force (MD values taken off)
	! #	    2 - vortex
	! #	    3 - vortex and elongation
	! #       4 - vortex generator
	! #	    5 - Couette Analytical solution from function_lib
	! #   COUPLED
	! #      -1 - Debug, apply nothing
	! #       0 - Zero force
	! #       1 - Coupled CFD
	! # Weighting Function
	! #	0 - No weighting -- apply evenly. 
	! #	1 - Weighting based on continuum stress only
	! #	2 - Weighting based on MD and continuum stress
	! # Start time of CV constraint
	! #	200 - default value
	! # Range of application
	! #	- xmin - minimum x coordinate in global system 
	! #	- xmax - maximum x coordinate in global system 
	! #	- ymin - minimum y coordinate in global system 
	! #	- ymax - maximum y coordinate in global system 
	! #	- zmin - minimum z coordinate in global system 
	! #	- zmax - maximum z coordinate in global system 
	! # Correction to velocity value
	! #	0 - Off. 
	! #	1 - On
	! # Number of steps to apply correction for
	! #   Nsteps - default value
	! # Direction to apply correction in
	! #   .true. or .false. - x
	! #   .true. or .false. - y
	! #   .true. or .false. - z
	! # -----------------------------------------------------------------------
	call locate(1,'CV_FORCES',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) CVforce_flag
		if (CVforce_flag .ne. VOID) then
			if (CV_debug .eq. 0) call error_abort("Input ERROR -- CV_FORCES true so CV_CONSERVE should be set to 1 and debugging set to > 1")
			if (vflux_outflag .ne. 4) call error_abort("Input ERROR -- CV_FORCES .true. but VFLUX_OUTFLAG not set to 4 (CV averages)")
		endif
		read(1,*,iostat=ios) CVweighting_flag
		if (ios .ne. 0) CVweighting_flag = 0
		read(1,*,iostat=ios) CVforce_starttime
		if (ios .ne. 0) CVforce_starttime = 200
        if (CVforce_starttime .le. 0) call error_abort("Error in read inputs -- CVforce starttime must be greater than 0")
		read(1,*,iostat=ios) F_CV_limits(1)
		if (ios .ne. 0) F_CV_limits(1) = VOID
		read(1,*,iostat=ios) F_CV_limits(2)
		if (ios .ne. 0) F_CV_limits(2) = VOID
		read(1,*,iostat=ios) F_CV_limits(3)
		if (ios .ne. 0) F_CV_limits(3) = VOID
		read(1,*,iostat=ios) F_CV_limits(4)
		if (ios .ne. 0) F_CV_limits(4) = VOID
		read(1,*,iostat=ios) F_CV_limits(5)
		if (ios .ne. 0) F_CV_limits(5) = VOID
		read(1,*,iostat=ios) F_CV_limits(6)
		if (ios .ne. 0) F_CV_limits(6) = VOID
		read(1,*,iostat=ios) CVforce_correct
		if (ios .ne. 0) CVforce_correct = 0
        if (CVforce_correct .eq. 1) then
    		read(1,*,iostat=ios) CVforce_correct_nsteps
	    	if (ios .ne. 0) CVforce_correct_nsteps = NSTEPS
        endif

        do ixyz = 1,3
    		read(1,*,iostat=ios) CVforce_applied_dir(ixyz)
            if (ios .ne. 0) CVforce_applied_dir(ixyz) = .true.
        enddo

        !print*, CVforce_flag, CVweighting_flag, CVforce_correct,  CVforce_starttime
	endif

	!Define specular wall location (if any)
	specular_wall = 0.d0
	call locate(1,'SPECULAR_WALL',.false.,found_in_input)
	if (found_in_input) then
		specular_flag = specular_flat
		read(1,*) specular_wall(1)			
		read(1,*) specular_wall(2)			
		read(1,*) specular_wall(3)
    	read(1,*,iostat=ios) specular_wall_flag
		if (ios .ne. 0) specular_wall_flag = 0
	endif

	call locate(1,'DEFINE_SHEAR',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) le_sd
		read(1,*) le_i0
		read(1,*) define_shear_as
		if (define_shear_as.eq.0) then
			read(1,*) le_sv
		endif
		if (define_shear_as.eq.1) then
			read(1,*) le_sr
		endif
		if (define_shear_as.gt.1) call error_abort( 'Poorly defined shear in input file')
	endif


	! #########################################################################
	! # Ignore tags in restart file and use location to set tags again
	! # [1] 0 - Off
	! # [1] 1 - On
	! # Debug mode write all tags to files thermo_tags, teth_tags__, etc
	! # [2] 0 - Off
	! # [2] 1 - On
	! # -----------------------------------------------------------------------
	call locate(1,'RESET_TAGS_ON_RESTART',.false.,found_in_input)
	if (found_in_input) then
        read(1,*) reset_tags_on_restart
        !Optional argument to check tags are setup correctly
        read(1,*,iostat=ios) debug_tags
		if (ios .ne. 0) debug_tags = .false.
    else
        reset_tags_on_restart = 0
    endif

	!-------------------------------------
	!Flag to determine molecular tags
	!-------------------------------------
	!Note: For initialunitsize "a"
	!                [  o     o ]
	!a (1 cell size) [     o    ]  a/2 (distance between molcules)	
	!                [  o     o ]
	!                [__________]  a/4 (distance from bottom of domain)
	!
	!So use (0.20+0.5d0*mol_layers)*initialunitsize(ixyz)

	!Set all to zero if no specifiers
	!Setup wall speeds
	wallslidev = 0.d0
	!Setup fixed molecules
	fixdistbottom = 0.d0;	fixdisttop = 0.d0
	!Setup sliding molecules
	slidedistbottom = 0.d0; slidedisttop = 0.d0
	!Setup molecules with tethered potentials
	tethereddistbottom = 0.d0; tethereddisttop = 0.d0
	!Setup thermostatted molecules
	thermstatbottom = 0.d0; thermstattop = 0.d0 
	!Setup regions to remove molecules (used for some boundary forces)
    emptydistbottom = 0.d0; emptydisttop = 0.d0

	if (ensemble .eq. 6) then
		! ########################################################################
		! # Velocity of sliding molecules in wall
		! # [1] float - x component of velocity of wall molecules
		! # [2] float - y component of velocity of wall molecules
		! # [3] float - z component of velocity of wall molecules
		! # ----------------------------------------------------------------------
		call locate(1,'WALLSLIDEV',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) wallslidev(1)
			read(1,*) wallslidev(2)
			read(1,*) wallslidev(3)
		endif

		! ########################################################################
		! # Distance from domain bottom to Fix molecules, i.e. force v=0 
		! # unless sliding
		! # [1] float - distance from x bottom
		! # [2] float - distance from y bottom
		! # [3] float - distance from z bottom
		! # ----------------------------------------------------------------------
		call locate(1,'FIXDISTBOTTOM',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) fixdistbottom(1)
			read(1,*) fixdistbottom(2)
			read(1,*) fixdistbottom(3)
		endif
		! ########################################################################
		! # Distance from domain top to Fix molecules, i.e. force v=0 
		! # unless sliding
		! # [1] float - distance from x top
		! # [2] float - distance from y top
		! # [3] float - distance from z top
		call locate(1,'FIXDISTTOP',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) fixdisttop(1)
			read(1,*) fixdisttop(2)
			read(1,*) fixdisttop(3)
		endif
		! ########################################################################
		! # Distance from domain bottom to apply sliding velocity to molecules 
		! # where applied velocity v=WALLSLIDEV	
		! # [1] float - distance from x bottom
		! # [2] float - distance from y bottom
		! # [3] float - distance from z bottom
		! # ----------------------------------------------------------------------
		call locate(1,'SLIDEDISTBOTTOM',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) slidedistbottom(1)
			read(1,*) slidedistbottom(2)
			read(1,*) slidedistbottom(3)
		endif
		! ########################################################################
		! # Distance from domain top to apply sliding velocity to molecules 
		! # where applied velocity v=WALLSLIDEV	
		! # [1] float - distance from x top
		! # [2] float - distance from y top
		! # [3] float - distance from z top
		! # ----------------------------------------------------------------------
		call locate(1,'SLIDEDISTTOP',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) slidedisttop(1)
			read(1,*) slidedisttop(2)
			read(1,*) slidedisttop(3)
		endif
		! ########################################################################
		! # Distance from domain bottom to tethered molecules using spring like
		! # restoring forces
		! # [1] float - distance from x top
		! # [2] float - distance from y top
		! # [3] float - distance from z top
		! # ----------------------------------------------------------------------
		call locate(1,'TETHEREDDISTBOTTOM',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) tethereddistbottom(1)
			read(1,*) tethereddistbottom(2)
			read(1,*) tethereddistbottom(3)
		endif
		! ########################################################################
		! # Distance from domain top to tethered molecules using spring like
		! # restoring forces
		! # [1] float - distance from x top
		! # [2] float - distance from y top
		! # [3] float - distance from z top
		! # ----------------------------------------------------------------------
		call locate(1,'TETHEREDDISTTOP',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) tethereddisttop(1)
			read(1,*) tethereddisttop(2)
			read(1,*) tethereddisttop(3)
		endif
			
		if (max(maxval(tethereddistbottom),maxval(tethereddisttop)).gt.0.d0) then
			tether_flag = 1
		else
			tether_flag = 0
		end if

		! #######################################################################
		! # Specifiy cooefficients of potential equation 
		! # phi= - k2*rio^2 - k4*rio^4 - k6*rio^6 in format
		! #  [1] float - k2
		! #  [2] float - k4
		! #  [3] float - k6
		! #
		! # Possible known combinations include:
		! #
		! # a)	Petravich and Harrowell (2006) J. Chem. Phys.124, 014103. 
		! # 		with constants  ( k2 = 0, k4 = 5,000, k6 = 5,000,000)
		! # b)	B. D. Todd, Peter J. Daivis, and Denis J. Evans (1995) 
		! #		Phys. Rev. E. 52, 5 with constants  (k2 = 28.575, k4 = 0, k6 = 0)
		! # c)  S. Y. Liem, D. Brown, and J. H. R. Clarke (1992) 
		! #		Phys. Rev. A. 45, 6 with constants  (k2 = 36.0,   k4 = 0, k6 = 0)
		! # Default Force constants (k2 = 0, k4 = 5,000, k6 = 5,000,000)  
		! # from Petravich and Harrowell (2006) J. Chem. Phys.124, 014103.
		! # ---------------------------------------------------------------------
		call locate(1,'TETHERCOEFFICIENTS',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) teth_k2
			read(1,*) teth_k4
			read(1,*) teth_k6
		else
			!Define default strength of tethering potential
			! phi = - k2*rio^2 - k4*rio^4 - k6*rio^6
			!Default Force constants (k2 = 0, k4 = 5,000, k6 = 5,000,000)  
			!from Petravich and Harrowell (2006) J. Chem. Phys.124, 014103.
			teth_k2=0.d0
			teth_k4=5000.d0
			teth_k6=5000000.d0
		endif

		! ########################################################################
		! # Apply wall texture - note must be used with tag move system
		! # and replaces all tethered wall specifications with a texture
		! # Different options have different meaning for parameters, for
		! # example 
		! # post, notch and roughness/fractal use texture_intensity
		! # to set depth of features
		! # Converging diverging nozzle:
		! # texture_intensity - size of throat
		! # Opt1 - fraction of domain for outlet region
		! # Opt2 - fraction of domain for diverging nozzle
		! # Opt3 - Set outlet region to still have wall
		! # [1] 0 - Off
		! # [1] 1 - posts
		! # [1] 2 - random spikes using sines/cosines
		! # [1] 3 - converging - diverging channel
		! # [1] 4 - fractal roughness
		! # [1] 5 - triangular notch
		! # [2] float - texture intensity
		! # [3] 0 - Thermostat applied to distance THERMSTATTOP and THERMSTATBOTTOM
		! # [3] 1 - Thermostat whole wall following texture
		! # [3] 2 - Slide whole wall including texture with WALLSLIDEV
		! # [4] float - wall texture parameter 1
		! # [5] float - wall texture parameter 2
		! # [6] float - wall texture parameter 3
		! # ----------------------------------------------------------------------
		call locate(1,'WALL_TEXTURE',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) texture_type
			read(1,*,iostat=ios) texture_intensity
			if (ios .ne. 0) then
				texture_intensity = 0.5d0
			endif
			read(1,*,iostat=ios) texture_therm
			if (ios .ne. 0) then
				texture_therm = 0
			endif
			read(1,*,iostat=ios) tex_opt1
			if (ios .ne. 0) then
				tex_opt1 = -666.d0
			endif
			read(1,*,iostat=ios) tex_opt2
			if (ios .ne. 0) then
				tex_opt2 = -666.d0
			endif
			read(1,*,iostat=ios) tex_opt3
			if (ios .ne. 0) then
				tex_opt3 = -666.d0
			endif
		else
			texture_type = 0
		endif
		! ########################################################################
		! # Distance from domain bottom to apply thermostat 
		! # [1] float - distance from x top
		! # [2] float - distance from y top
		! # [3] float - distance from z top
		! # ----------------------------------------------------------------------
		call locate(1,'THERMSTATBOTTOM',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) thermstatbottom(1)
			read(1,*) thermstatbottom(2)
			read(1,*) thermstatbottom(3)
		endif
		! ########################################################################
		! # Distance from domain top to apply thermostat 
		! # [1] float - distance from x top
		! # [2] float - distance from y top
		! # [3] float - distance from z top
		! # ----------------------------------------------------------------------
		call locate(1,'THERMSTATTOP',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) thermstattop(1)
			read(1,*) thermstattop(2)
			read(1,*) thermstattop(3)
		endif
		! ########################################################################
		! # Distance from domain top to remove particles
		! # leaving a buffer region which prevents molecules escaping if tethered
		! # or to strip out outer wall molecules to save computer time
		! # [1] float - distance from x top
		! # [2] float - distance from y top
		! # [3] float - distance from z top
		! # ----------------------------------------------------------------------
		call locate(1,'EMPTYDISTTOP',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) emptydisttop(1)
			read(1,*) emptydisttop(2)
			read(1,*) emptydisttop(3)
		endif
		! ########################################################################
		! # Distance from domain bottom to remove particles
		! # leaving a buffer region which prevents molecules escaping if tethered
		! # or to strip out outer wall molecules to save computer time
		! # [1] float - distance from x top
		! # [2] float - distance from y top
		! # [3] float - distance from z top
		! # ----------------------------------------------------------------------
		call locate(1,'EMPTYDISTBOTTOM',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) emptydistbottom(1)
			read(1,*) emptydistbottom(2)
			read(1,*) emptydistbottom(3)
		endif

		! ########################################################################
		! # Apply heating in a local region
		! # [1] float - start of heated region in x
		! # [2] float - end of heated region in x
		! # [3] float - start of heated region in y
		! # [4] float - end of heated region in y
		! # [5] float - start of heated region in z
		! # [6] float - end of heated region in z
		! # ----------------------------------------------------------------------
		call locate(1,'LOCAL_HEAT',.false.,found_in_input)
		if (found_in_input) then
			read(1,*) local_heat_region(1)
			read(1,*) local_heat_region(2)
			read(1,*) local_heat_region(3)
			read(1,*) local_heat_region(4)
			read(1,*) local_heat_region(5)
			read(1,*) local_heat_region(6)
		else
			local_heat_region = -666.d0
		endif

	endif

	! ########################################################################
	! ## NEWPAGE - Outputs
	! ########################################################################
	call locate(1,'NEWPAGE_OUTPUTS',.false.,found_in_input)
	if (found_in_input) then
		!read(1,*) newpage
		print*, "The keyword OUTPUT does nothing, "
		print*, "it is used to denote start of output section flowmol_inputs" 
	endif


    !print*, "INTRINSIC_INTERFACE inputs", intrinsic_interface_outflag, II_normal, II_alpha, II_tau, II_eps, II_ns      


	!Flag to determine if output is switched on
	call locate(1,'VMD_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) vmd_outflag
		if (vmd_outflag .ne. 0) then
			read(1,*,iostat=ios) Nvmd_intervals	!Number of vmd intervals
            !If zero intervals or not specified, switch on for entire time
			if (Nvmd_intervals .eq. 0 .or. ios .ne. 0) then
				allocate(vmd_intervals(2,1))
				vmd_intervals(1,1) = 1; vmd_intervals(2,1) = Nsteps
				Nvmd_intervals = 1
            !Otherwise, try to read intervals from next line
			else
				allocate(vmd_intervals(2,Nvmd_intervals))
                !Check if interval in form of comma seperated list of inputs
                read(1,'(a)',iostat=ios) str
                backspace(1)
                if (scan(str, ",").gt.0) then
                    read(1,*,iostat=ios) vmd_intervals
                !Otherwise, use Nvmd_interval_size to specify linearly space records
                else
                    !See if a single interval size is specified, otherwise use 1000
                    read(1,*,iostat=ios) Nvmd_interval_size
                    if (ios .ne. 0) Nvmd_interval_size = 1000
                    vmd_intervals(1,:) = linspace(0.d0, & 
                                                  real(Nsteps,kind(0.d0)),&
                                                  Nvmd_intervals)
                    vmd_intervals(2,:) = vmd_intervals(1,:) + Nvmd_interval_size
                    !Note, convention to have last interval from Nsteps 
                    !interval_size to Nsteps
                    vmd_intervals(1,size(vmd_intervals,2)) = & 
                        vmd_intervals(1,size(vmd_intervals,2)) - Nvmd_interval_size
                    vmd_intervals(2,size(vmd_intervals,2)) = & 
                        vmd_intervals(2,size(vmd_intervals,2)) - Nvmd_interval_size
                endif

#if USE_COUPLER
				!Coupler total simulation time is setup later so defer this check
				!until later
#else
				if (maxval(vmd_intervals) .gt. Nsteps) then
					print'(2(a,i8))', 'Value specified for end of final vmd_interval = ' & 
									, maxval(vmd_intervals), 'but Nsteps = ', Nsteps 
					call error_abort("Specified VMD interval greater than Nsteps")
				endif
#endif
			endif
		endif
	else
		!If not switched on in input then VMD set to off
		vmd_outflag = 0
	endif

    call locate(1,'VMD_SKIP',.false.,found_in_input)
    if (found_in_input) then
        read(1,*) vmd_skip  
        if (vmd_skip .lt. 1) then
            call error_abort('VMD_SKIP cannot be less than 1')
        end if
    else
        vmd_skip = 1
    end if

	call locate(1,'SEPERATE_OUTFILES',.false.,found_in_input)
	if (found_in_input) then
        call error_abort('SEPERATE_OUTFILES error, sepArate is misspelt')
    endif
    
	! #########################################################################
	! # Save output files as a single file or seperate file at each timestep
	! # [1] .false. - Output a single file per output
	! # [1] .true. - Output a value per timestep and per output type 
	! # [2] .false. - On restart, use numbering starting from current record (initialstep)
	! # [2] .true. - Delete all previous files and create file numbes
	! # -----------------------------------------------------------------------
	call locate(1,'SEPARATE_OUTFILES',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) separate_outfiles
        if (separate_outfiles) then
    		read(1,*,iostat=ios) restart_numbering
			if (ios .ne. 0) restart_numbering = .false.
        endif
	endif

	! #########################################################################
	! # Define the number of output bins in terms of the compuational cells
	! # This constraint is useful for efficiency (cell lists) and consistency
	! # while not being particularly restrictive (bins must be integer no in
	! # each process and so must cells). The bin cell ratio in the x,y,z 
	! # directions is specified in real format with e.g. 2.0 is 2 bins per
	! # cell or 0.5 is 2 cells per bin. Care must be taken that cells is an
	! # even number if splitting.
	! # [1] float - x ratio of bins per cell
	! # [2] float - y ratio of bins per cell
	! # [3] float - z ratio of bins per cell
	! # -----------------------------------------------------------------------
	call locate(1,'BIN2CELLRATIO',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) binspercell(1)
		read(1,*) binspercell(2)
		read(1,*) binspercell(3)
!		if (any(binspercell .gt. 1.d0) .and. & 
!			any(binspercell .lt. 1.d0)) then
!			print*, "WARNING -- BIN2CELLRATIO input - cannot specify multiple ", &
!					"bins per cell in one direction and multiple cells per bin ", &
!					"in the other -- setting minimum binspercell to 1"
!			where(binspercell .gt. 1.d0) binspercell = 1.d0
!		endif
	else
		binspercell = 1.d0
	endif

	! #########################################################################
	! # Output flag for macroscopic properties:
	! # [1] 0 - off
	! # [1] 1 - high precision > stdout
	! # [1] 2 - high precision > stdout + results/macroscopic_properties
	! # [1] 3 - concise        > stdout
	! # [1] 4 - concise        > stdout + results/macroscopic_properties
	! # -----------------------------------------------------------------------
	call locate(1,'MACRO_OUTFLAG',.false.,found_in_input)
	if (found_in_input) read(1,*) macro_outflag

	! #########################################################################
	! # Add correction to energy and virial pressure in macro output based
	! # standard long range correction (see an MD textbook e.g. Rapaport or 
	! # Allen and Tildesley)
	! # [1] 0 - Off
	! # [1] 1 - On
	! # -----------------------------------------------------------------------
	call locate(1,'SLRC_FLAG',.false.,found_in_input)
	if (found_in_input) read(1,*) sLRC_flag
	!Test for WCA potential and switch LRC off
	if (abs(rcutoff-2.d0**(1.d0/6.d0)) .lt. 0.0001) then
		if (sLRC_flag .eq. 1) then
			print*, "WCA potential used - switching sLRC off"
			sLRC_flag = 0
		endif
	endif

	! #########################################################################
	! # Output flag for mass record:
	! # [1] 0 - off 
	! # [1] 4 - 3D grid of bins
	! # [2] int - No. of samples for mass average
	! # -----------------------------------------------------------------------
	call locate(1,'MASS_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) mass_outflag
		if (mass_outflag .ne. 0) then
			read(1,*) Nmass_ave
		endif
		if (mass_outflag .eq. 5) then
			call locate(1,'CPOL_BINS',.true.)
			read(1,*) gcpol_bins(1)	
			read(1,*) gcpol_bins(2)	
			read(1,*) gcpol_bins(3)	
		end if
	endif


	! #########################################################################
	! # Output flag for momentum record:
	! # [1] 0 - off 
	! # [1] 4 - 3D grid of bins
	! # [2] int - No. of samples for momentum average
	! # -----------------------------------------------------------------------
	call locate(1,'MOMENTUM_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) momentum_outflag
		if (momentum_outflag .ne. 0) then
			read(1,*) Nvel_ave
		endif
		if (momentum_outflag .eq. 5) then
			call locate(1,'CPOL_BINS',.true.)
			read(1,*) gcpol_bins(1)	
			read(1,*) gcpol_bins(2)	
			read(1,*) gcpol_bins(3)	
		end if
	endif
	! #########################################################################
	! # Output flag for temperature record:
	! # [1] 0 - off 
	! # [1] 4 - 3D grid of bins
	! # [2] int - No. of samples for temperature average
	! # [3] 0 - Peculiar velocity off
	! # [3] 1 - Peculiar velocity on
	! # -----------------------------------------------------------------------
	call locate(1,'TEMPERATURE_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) temperature_outflag
		if (temperature_outflag .ne. 0) then
			read(1,*) NTemp_ave
			read(1,*,iostat=ios) peculiar_flag
			if (ios .ne. 0) peculiar_flag = 0 !default to zero if value not found
		endif
	endif
	! #########################################################################
	! # Output flag for energy record:
	! # [1] 0 - off 
	! # [1] 4 - 3D grid of bins
	! # [2] int - No. of samples for energy average
	! # -----------------------------------------------------------------------
	call locate(1,'ENERGY_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) energy_outflag
		if (energy_outflag .ne. 0) then
			read(1,*) Nenergy_ave
		endif
	endif
	! #########################################################################
	! # Output flag for bin centre of mass record:
	! # [1] 0 - off 
	! # [1] 4 - 3D grid of bins
	! # [2] int - No. of samples for centre of mass average
	! # -----------------------------------------------------------------------
	call locate(1,'CENTRE_OF_MASS_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) centre_of_mass_outflag
		if (centre_of_mass_outflag .ne. 0) then
			read(1,*) Ncom_ave
		endif
	endif
	call locate(1,'PRESSURE_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) pressure_outflag
		if (pressure_outflag .ne. 0) then
			!Stress averaging
			read(1,*) Nstress_ave
			!Split kinetic/config
			read(1,*,iostat=ios) split_kin_config
			if (ios .ne. 0) split_kin_config = 0 !default to zero if value not found
			!Check other options such as calculation method and 
			if (pressure_outflag .eq. 2 .or. pressure_outflag .eq. 3) then
    			read(1,*,iostat=ios) VA_calcmethod
				if (ios .ne. 0) VA_calcmethod = 1
    			if (any(binspercell .gt. 1.d0) .and. VA_calcmethod .eq. 2) then
    				print*, "WARNING -- Cannot specify multiple bins per cell with exact VA "
    				print*, "   calculation (method 2), switching to trapizium (method 1)   "
					VA_calcmethod = 1
    			endif
				if (VA_calcmethod .eq. 1) then
					read(1,*,iostat=ios) VA_line_samples
    				if (ios .ne. 0)	VA_line_samples = 20
				endif
			endif
		endif
		if (pressure_outflag .eq. 3) then
			call locate(1,'CPOL_BINS',.true.)
			read(1,*) gcpol_bins(1)	
			read(1,*) gcpol_bins(2)	
			read(1,*) gcpol_bins(3)	
		end if
	endif

	call locate(1,'HEATFLUX_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) heatflux_outflag
		if (heatflux_outflag .ne. 0) then
			!Stress averaging
			read(1,*) Nheatflux_ave
			!Split kinetic/config
			read(1,*,iostat=ios) split_hfkin_config
			if (ios .ne. 0) split_hfkin_config = 0 !default to zero if value not found
			!Check other options such as calculation method and 
			if (heatflux_outflag .eq. 2) then
				if (pressure_outflag .ne. 2) then
					call error_abort("Error -- Volume average stress must be recorded for heat flux")
				endif
			endif
			if (heatflux_outflag .eq. 2 .or. heatflux_outflag .eq. 3) then
    			read(1,*,iostat=ios) VA_heatflux_calcmethod
				if (ios .ne. 0) VA_heatflux_calcmethod = 1
    			if (any(binspercell .gt. 1.d0) .and. VA_heatflux_line_samples .eq. 2) then
    				print*, "WARNING -- Cannot specify multiple bins per cell with exact VA "
    				print*, "   calculation (method 2), switching to trapizium (method 1)   "
					VA_heatflux_calcmethod = 1
    			endif
				if (VA_heatflux_calcmethod .eq. 1) then
					read(1,*,iostat=ios) VA_heatflux_line_samples
    				if (ios .ne. 0)	VA_heatflux_line_samples = 20
				endif
			endif
		endif
		if (heatflux_outflag .eq. 3) then
			call error_abort("Error -- heat flux not developed for cpol bins")
		end if
	endif

	call locate(1,'VISCOSITY_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) viscosity_outflag
		if ( viscosity_outflag .ne. 0) then
			read(1,*) Nvisc_ave
		endif
	endif
	! #########################################################################
	! # Control Volume Conservation Averaging
	! #	CV Conservation averaging (0=off 1=on) - take mass, momentum or 
	! #	energy flux measure every step to ensureflux will be equal to change 
	! #	in snapshots. Note, this increases computational cost somewhat
	! #   Debug mode .true. or .false. can be used to enforce in code CV 
	! #	conservation checking
	! #
	! # -----------------------------------------------------------------------
	call locate(1,'CV_CONSERVE',.false.,found_in_input)
	cv_conserve = 0
	if (found_in_input) then
		read(1,*) cv_conserve
        if (cv_conserve .ne. 0) then
		    read(1,*,iostat=ios) CV_debug
		    if (ios .ne. 0) CV_debug = 0
            if (CV_debug .eq. 2) then
		        read(1,*,iostat=ios) debug_CV(1)
                read(1,*,iostat=ios) debug_CV(2)
                read(1,*,iostat=ios) debug_CV(3)
        		if (ios .ne. 0) then
                    print*, "Warning - CV_debug = 2 so CV number should be specified, setting to 3,3,3"
                    debug_CV = (/ 3, 3, 3 /)
                endif
				!Try to keep reading to see if range specified
		        read(1,*,iostat=ios) debug_CV_range(1)				
				read(1,*,iostat=ios) debug_CV_range(2)
				read(1,*,iostat=ios) debug_CV_range(3)
        		if (ios .ne. 0) then
                    debug_CV_range = (/ 1, 1, 1 /)
                endif
            endif
        endif
	endif

	call locate(1,'PASS_VHALO',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) pass_vhalo
    else
        pass_vhalo = 0
	endif

	! #########################################################################
	! # Output flag for mass flux control volume analysis
	! # Will save mass snapshots msnap and surface fluxes mflux.
	! # [1] 0 - Off
	! # [1] 1 - On (3D grid using local control volumes)
	! # [2] int - No. of samples for mass flux & interval for CV mass snapshot
	! #
	! # -----------------------------------------------------------------------
	call locate(1,'MFLUX_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) mflux_outflag
		if (mflux_outflag .ne. 0) then
			read(1,*) Nmflux_ave
		endif
	endif
	! #########################################################################
	! # Output flag for momentum flux over a surface for control volume analysis
	! # Will save momentum snapshots vsnap and surface fluxes vflux as well as	
	! # force interactions crossing the surface of the control volume  
	! # [1] 0 - Off
	! # [1] 1 - Method of planes in x (depreciated)
	! # [1] 2 - Method of planes in y (depreciated)
	! # [1] 3 - Method of planes in z (depreciated)
	! # [1] 4 - 3D grid using local control volumes
	! # [2] int - No. of samples for momentum flux & interval for snapshot
	! # -----------------------------------------------------------------------
	call locate(1,'VFLUX_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) vflux_outflag
		if (vflux_outflag .ne. 0) then
			read(1,*) Nvflux_ave
		endif
		if (mflux_outflag .eq. 0) Nmflux_ave = Nvflux_ave !Mass set to same as velocity
	endif
	! #########################################################################
	! # Output flag for energy flux/power & CV energy snapshots over 
	! # a surface for control volume analysis
	! # Will save energy snapshots esnap and surface fluxes eflux as well as	
	! # force time velocity interactions crossing the surface of the control volume  
	! # [1] 0 - Off
	! # [1] 1 - Method of planes in x (depreciated)
	! # [1] 2 - Method of planes in y (depreciated)
	! # [1] 3 - Method of planes in z (depreciated)
	! # [1] 4 - 3D grid using local control volumes
	! # [2] int - No. of samples for energy flux & interval for snapshot
	! # -----------------------------------------------------------------------
	call locate(1,'EFLUX_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) eflux_outflag
		if (eflux_outflag .ne. 0) then
			read(1,*) Neflux_ave
			pass_vhalo = 1		!Turn on passing of velocities for halo images
		endif
	endif
	if (CV_debug .gt. 0) then
		if (mflux_outflag .eq. 0 .and. & 
			vflux_outflag .eq. 0 .and. & 
			eflux_outflag .eq. 0) then
			call error_abort("If CV_debug is true, mass/momentum/energy flux must be turned on")
		endif
	endif
	call locate(1,'MSURF_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
        read(1,*) msurf_outflag
        read(1,*) Nsurfm_ave
    else
        msurf_outflag = 0
        Nsurfm_ave = 0
    endif

	call locate(1,'SURF_EVO_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
        read(1,*) Nsurfevo_outflag
        read(1,*) Nsurfevo_ave
    else
        Nsurfevo_outflag = 0
        Nsurfevo_ave = 0
    endif

	call locate(1,'ETEVTCF_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) etevtcf_outflag
		if (etevtcf_outflag.ne.0) then
			read(1,*) etevtcf_iter0
	
			if (mod(etevtcf_iter0,tplot).ne.0) then
				etevtcf_iter0 = etevtcf_iter0 + (tplot - mod(etevtcf_iter0,tplot))
				print*, 'Etevtcf must be a multiple of tplot, resetting etevtcf to ', etevtcf_iter0
			end if
		end if
	endif

	call locate(1,'RTRUE_FLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) rtrue_flag
	else
		rtrue_flag = 0
	endif

	call locate(1,'R_GYRATION_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) r_gyration_outflag
		read(1,*) r_gyration_iter0
        if (r_gyration_outflag .ne. 0) rtrue_flag = 1
	endif

	call locate(1,'RDF_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) rdf_outflag
		read(1,*) rdf_rmax
		read(1,*) rdf_nbins
	endif

	call locate(1,'VPDF',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) vPDF_flag
		read(1,*) NvPDF_ave
		read(1,*) NPDFbins 
		read(1,*) PDFvlims
    else
        vPDF_flag = 0
	endif
	
	call locate(1,'STRUCT_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) ssf_outflag
		read(1,*) ssf_ax1
		read(1,*) ssf_ax2 
		read(1,*) ssf_nmax 
	endif

    split_pol_sol_stats = 0
    call locate(1,'SPLIT_POL_SOL_STATS',.false.,found_in_input)
    if (found_in_input) then
        read(1,*) split_pol_sol_stats
    else
        split_pol_sol_stats = 0
    end if


	! ########################################################################
	! # Turn on cluster analysis to build clusters from molecules
	! # a linked list of all molecules within a cutoff distance of each other
	! # [1] 0 - Cluster Analysis Off 
	! # [1] 1 - Cluster Analysis On - allows intrinsic interface
	! # [1] 2 - Cluster Analysis On - for average mass bin interface (OBSOLETE)
	! # [2] float - Cutoff length for cluster search
	! # [3] int - Minimum number of neighbours for inclusion in cluster
	! # [4] 1 - Write interface as an xyz file (for VMD) with all molecuels in cluster.xyz 
	! # [4] 2 - Write interface as an obj file (for Blender, etc)
	! # [4] 3 - Write interface as surface.grid, a simple ascii grid of center locations 
	! # [4] 4 - Write modes to ascii file 
	! # [5] int - Resolution to write interface (Number of points in each direction)
	! # ----------------------------------------------------------------------
	call locate(1,'CLUSTER_ANALYSIS',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) cluster_analysis_outflag
        if (cluster_analysis_outflag .ne. 0) then
            if (nproc .ne. 1) then
                call error_abort("Cluster Analysis only works with one processor")
            end if
		    read(1,*,iostat=ios) CA_rd   ! Cutoff length for cluster search
            if (ios .ne. 0) CA_rd = 1.5
		    read(1,*,iostat=ios) CA_min_nghbr   ! Minimum number of neighbours
            if (ios .ne. 0) CA_min_nghbr = 0  ! Set to zero (i.e. default no minimum)
		    read(1,*,iostat=ios) CA_generate_xyz   ! Output xyz files for vmd
            if (ios .ne. 0) CA_generate_xyz = 0  ! Set to zero (i.e. default no output)
            if (CA_generate_xyz .ne. 0) then
    		    read(1,*,iostat=ios) CA_generate_xyz_res   ! Resolution for output xyz files 
                if (ios .ne. 0) CA_generate_xyz_res = 0  
            endif
            ! If interface cutoff is less that interaction rcutoff
            ! then we can use the neighbourlist to get molecules in 
            ! interface region (N.B. need to use all interations)
            if (CA_rd .gt. (rcutoff + minval(delta_rneighbr))) then
                call error_abort("Error in build cluster -- rd must be less than neighbourlist cutoff")
            endif

            if ((force_list .lt. 1) .or. (force_list .gt. 2)) then
                call error_abort("Error in build_from_cellandneighbour_lists -- full "//&
                                 "neightbour list should be used with interface tracking."//&
                                 " Set FORCE_LIST to 2 in input file.")
            end if
        endif
    else
        cluster_analysis_outflag = 0
	endif

    !print*, "CLUSTER_ANALYSIS inputs", CA_rd, CA_min_nghbr

    !#########################################################################
    !# Fit an intrinsic surface to the outside of the cluster
    !#   flag   1 - Intrinsic sine/cosine  2 - sine/cosine with bilinear approx  
    !#       - normal     Surface normal direction
    !#       - alpha      Smallest wavelength
    !#       - tau        Search radius around surface
    !#       - omega      Weight for surface energy minimising constraint
    !#       - ns         Target density of surface Npivots/Area
    !#   flag       3  - Linear and cubic surfaces
    !#       - No flags yet
    !# -----------------------------------------------------------------------
	call locate(1,'INTRINSIC_INTERFACE',.false.,found_in_input)
	if (found_in_input) then
        if (cluster_analysis_outflag .eq. 0) then
            call error_abort("Cluster Analysis must be on to use intrinsic interface")
        endif
		read(1,*) intrinsic_interface_outflag
        if (intrinsic_interface_outflag .ne. 0) then
    		read(1,*,iostat=ios) II_normal   ! Surface normal direction
            if (ios .ne. 0) II_normal = 3
            read(1,*,iostat=ios) II_alpha    ! Smallest wavelength
            if (ios .ne. 0) II_alpha = 0.5d0
            read(1,*,iostat=ios) II_tau      ! Search radius around surface
            if (ios .ne. 0) II_tau = 1.d0
            read(1,*,iostat=ios) II_eps    ! Weight for surface energy minimising constraint
            if (ios .ne. 0) II_eps = 0.00000001d0
            read(1,*,iostat=ios) II_ns       ! Target density of surface Npivots/Area
            if (ios .ne. 0) II_ns = 0.8d0
            read(1,*,iostat=ios) II_topbot       !Top =1 or bottom=2
            if (ios .ne. 0) II_topbot = 1
        endif
    else
        intrinsic_interface_outflag = 0
	endif


	close(1,status='keep')      !Close input file

end subroutine setup_read_input
