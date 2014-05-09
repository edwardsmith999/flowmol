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

end module module_read_input
!----------------------------------------------------------------------------------

subroutine setup_read_input
	use module_read_input
	use librarymod, only : locate
	implicit none

	logical					:: found_in_input, error
	integer 				:: ios, ixyz


	! Open input file
	open(1,file=input_file)

	!Simulation ensemble choice
	call locate(1,'ENSEMBLE',.true.)
	read(1,*) ensemble

	!Setup initial configuration
	call locate(1,'INITIAL_CONFIG_FLAG',.true.)
	read(1,*) initial_config_flag 
	select case (initial_config_flag)
	case(0)
	
		potential_flag = 0

		call locate(1,'DENSITY',.true.)
		read(1,*) density
		call locate(1,'RCUTOFF',.true.)
		read(1,*) rcutoff
		call locate(1,'INITIALNUNITS',.true.)
		read(1,*) initialnunits(1)		!x dimension split into number of cells
		read(1,*) initialnunits(2)		!y dimension split into number of cells
		read(1,*) initialnunits(3)		!z dimension split into number of cells

	case(1)	

		read(1,*) config_special_case	
		select case (trim(config_special_case))
		case('sparse_fene')	

			potential_flag = 1	
			rcutoff = 2.d0**(1.d0/6.d0)

			call locate(1,'SPARSE_FENE',.true.)
			read(1,*) nmonomers
			read(1,*) k_c
			read(1,*) R_0
			read(1,*) globaldomain(1)
			read(1,*) globaldomain(2)
			read(1,*) globaldomain(3)
			read(1,*) nchains

			density = nmonomers * nchains / product(globaldomain(1:3))	

			call locate(1,'SOLVENT_INFO',.true.)
			read(1,*) solvent_flag
			select case(solvent_flag)
			case(0)
			case(1)
				read(1,*) solvent_ratio
			case(2)
				read(1,*) solvent_ratio
				read(1,*) eps_pp
				read(1,*) eps_ps
				read(1,*) eps_ss
			case default
				call error_abort('Unrecognised solvent flag!')
			end select

		case('dense_fene')	

			potential_flag = 1	
			rcutoff = 2.d0**(1.d0/6.d0)

			call locate(1,'DENSE_FENE',.true.)
			read(1,*) nmonomers
			read(1,*) k_c
			read(1,*) R_0
			read(1,*) density
			read(1,*) initialnunits(1)
			read(1,*) initialnunits(2)
			read(1,*) initialnunits(3)	

			call locate(1,'SOLVENT_INFO',.true.)
			read(1,*) solvent_flag
			select case(solvent_flag)
			case(0)
			case(1)
				read(1,*) solvent_ratio
			case(2)
				read(1,*) solvent_ratio
				read(1,*) eps_pp
				read(1,*) eps_ps
				read(1,*) eps_ss
			case default
				call error_abort('Unrecognised solvent flag!')
			end select

		case('solid_liquid')

			potential_flag = 0

			call locate(1,'DENSITY',.true.)
			read(1,*) density
			call locate(1,'LIQUIDDENSITY',.true.)
			read(1,*) liquid_density
			call locate(1,'RCUTOFF',.true.)
			read(1,*) rcutoff
			call locate(1,'INITIALNUNITS',.true.)
			read(1,*) initialnunits(1)		!x dimension split into number of cells
			read(1,*) initialnunits(2)		!y dimension split into number of cells
			read(1,*) initialnunits(3)		!z dimension split into number of cells


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
			ensemble = tag_move

			call locate(1,'DENSITY',.true.)
			read(1,*) density
			call locate(1,'RCUTOFF',.true.)
			read(1,*) rcutoff

		case('fill_cylinders_fene_solution')

			potential_flag = 1	    
            ensemble = tag_move
			rcutoff = 2.d0**(1.d0/6.d0)

            call locate(1,'FENE_SOLUTION',.true.)
			read(1,*) nmonomers
			read(1,*) k_c
			read(1,*) R_0
			read(1,*) density
            read(1,*) targetconc

		case('rotate_cylinders')
			
			ensemble = tag_move

			!call locate(1,'RCUTOFF',.true.)
			!read(1,*) rcutoff
			call locate(1,'POTENTIAL_FLAG',.true.)
            read(1,*) potential_flag

			call locate(1,'ROTATE_CYLINDERS',.true.)
			read(1,*) omega_i
			read(1,*) omega_f 
			read(1,*) omega_ramplength

		case default

			call error_abort("ERROR -- Unrecognised special case string")

		end select

	case(2)
		call error_abort("ERROR -- Generic input file read in is not developed yet")
	case default
		call error_abort("ERROR -- Unrecognised initial_config_flag")
	end select 

	!Read in initial temperature
	call locate(1,'INPUTTEMPERATURE',.true.)
	read(1,*) inputtemperature

	!Setup velocity initial condition
	call locate(1,'INITIAL_VELOCITY_FLAG',.false.,found_in_input) 
	if (found_in_input) then
		read(1,*) initial_velocity_flag 
   		read(1,*) velocity_special_case	
   		select case (trim(velocity_special_case	))
   		case('debug')
			!Nothing must be read in
   		case('taylor_green')
			!Nothing must be read in
   		case('dns')
	   		read(1,*) DNS_filename
	   		read(1,*) DNS_ngx
	   		read(1,*) DNS_ngy
	   		read(1,*) DNS_ngz
   		case default
   			call error_abort('Unidentified initial velocities_special_case')	
   		end select
	else
		initial_velocity_flag = 0
	endif

	call locate(1,'INTEGRATION_ALGORITHM',.true.)
	read(1,*) integration_algorithm
	call locate(1,'FORCE_LIST',.true.)	!LJ or FENE potential
	read(1,*) force_list
	if (force_list .ne. 3 .and. ensemble .eq. 4) &
		call error_abort("Half int neighbour list only is compatible with pwa_terms_pwaNH thermostat")
	if (force_list .ne. 3 .and. ensemble .eq. 5) & 
		call error_abort("Half int neighbour list only is compatible with DPD thermostat")
	!Input computational co-efficients
	call locate(1,'NSTEPS',.true.)
	read(1,*) Nsteps 		!Number of computational steps
	call locate(1,'DELTA_T',.true.)
	read(1,*) delta_t 		!Size of time step
	call locate(1,'TPLOT',.true.)
	read(1,*) tplot 		!Frequency at which to record results
	call locate(1,'INITIALISE_STEPS',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) initialise_steps 	!Number of initialisation steps for simulation
	else
		initialise_steps = 0
	endif
	call locate(1,'DELTA_RNEIGHBR',.true.) 
	read(1,*) delta_rneighbr 	!Extra distance used for neighbour cell

	call locate(1,'FIXED_REBUILD_FLAG',.false.,found_in_input) 
	if (found_in_input) then
		read(1,*) fixed_rebuild_flag 	!Fixed rebuild flag
		read(1,*) fixed_rebuild 		!Fixed rebuild frequency
	else
		fixed_rebuild_flag = 0			!Rebuild uses neighbourcell
	endif

	call locate(1,'RESCUE_SNAPSHOT_FREQ',.false.,found_in_input) 
	if (found_in_input) then
		read(1,*) rescue_snapshot_freq 	!Rescue snapshot frequency in seconds
	else
		rescue_snapshot_freq = 21600	!Every 6 hours
	endif

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

	call locate(1,'SEED',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) seed(1) 	!Random number seed value 1
		read(1,*) seed(2) 	!Random number seed value 2
	else
		seed(1) = 1		!Fixed default seed for repeatability
		seed(2) = 2		!Fixed default seed for repeatability
	endif

	!Flags to determine if periodic boundaries are on or shearing Lees Edwards
	call locate(1,'PERIODIC',.true.)
	read(1,*) periodic(1)
	read(1,*) periodic(2)
	read(1,*) periodic(3)

	if (any(periodic.eq.0)) then

		bforce_flag(:) = 0
		bforce_dxyz(:) = 0.0
	
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
				if (irank.eq.iroot) print*, 'Warning, resetting bforce_flag(1:2)' , &
				' because periodic boundaries are on in the x-direction'
				bforce_flag(1:2) = 0
			end if
			if (periodic(2).ne.0) then 
				if (irank.eq.iroot) print*, 'Warning, resetting bforce_flag(3:4)' , &
				' because periodic boundaries are on in the y-direction'
				bforce_flag(3:4) = 0
			end if
			if (periodic(3).ne.0) then 
				if (irank.eq.iroot) print*, 'Warning, resetting bforce_flag(5:6)' , &
				' because periodic boundaries are on in the z-direction'
				bforce_flag(5:6) = 0
			end if

		end if

	end if

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

	!Apply force to region in space
	call locate(1,'EXTERNAL_FORCE',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) external_force_flag
		if (external_force_flag .eq. 1) then
			read(1,*) F_ext_ixyz
			read(1,*) F_ext
		elseif (external_force_flag .eq. 2) then
			read(1,*) F_ext_ixyz
			read(1,*) F_ext
			read(1,*) F_ext_limits(1)
			read(1,*) F_ext_limits(2)
			read(1,*) F_ext_limits(3)
			read(1,*) F_ext_limits(4)
			read(1,*) F_ext_limits(5)
			read(1,*) F_ext_limits(6)
		endif
	else
		external_force_flag = 0
	endif

	!Define specular wall location (if any)
	specular_wall = 0.d0
	call locate(1,'SPECULAR_WALL',.false.,found_in_input)
	if (found_in_input) then
		specular_flag = specular_flat
		read(1,*) specular_wall(1)			
		read(1,*) specular_wall(2)			
		read(1,*) specular_wall(3)
	endif

	call locate(1,'DEFINE_SHEAR',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) le_sd
		read(1,*) le_i0
		read(1,*) define_shear_as
		if (define_shear_as.eq.0) read(1,*) le_sv
		if (define_shear_as.eq.1) read(1,*) le_sr
		if (define_shear_as.gt.1) call error_abort( 'Poorly defined shear in input file')
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
	
	call locate(1,'WALLSLIDEV',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) wallslidev(1)
		read(1,*) wallslidev(2)
		read(1,*) wallslidev(3)
	endif
	call locate(1,'FIXDISTBOTTOM',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) fixdistbottom(1)
		read(1,*) fixdistbottom(2)
		read(1,*) fixdistbottom(3)
	endif
	call locate(1,'FIXDISTTOP',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) fixdisttop(1)
		read(1,*) fixdisttop(2)
		read(1,*) fixdisttop(3)
	endif
	call locate(1,'SLIDEDISTBOTTOM',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) slidedistbottom(1)
		read(1,*) slidedistbottom(2)
		read(1,*) slidedistbottom(3)
	endif
	call locate(1,'SLIDEDISTTOP',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) slidedisttop(1)
		read(1,*) slidedisttop(2)
		read(1,*) slidedisttop(3)
	endif
	call locate(1,'TETHEREDDISTBOTTOM',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) tethereddistbottom(1)
		read(1,*) tethereddistbottom(2)
		read(1,*) tethereddistbottom(3)
	endif
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
	else
		texture_type = 0
	endif
	call locate(1,'THERMSTATBOTTOM',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) thermstatbottom(1)
		read(1,*) thermstatbottom(2)
		read(1,*) thermstatbottom(3)
	endif
	call locate(1,'THERMSTATTOP',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) thermstattop(1)
		read(1,*) thermstattop(2)
		read(1,*) thermstattop(3)
	endif

	!Flag to determine if output is switched on
	call locate(1,'VMD_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) vmd_outflag
		if (vmd_outflag .ne. 0) then
			read(1,*,iostat=ios) Nvmd_intervals	!Number of vmd intervals
			if (Nvmd_intervals .gt. 20) then
				print*, "Number of VMD intervals greater than 20 or not specified, setting on for all simulation"
				Nvmd_intervals = 0
			endif
			if (Nvmd_intervals .eq. 0 .or. ios .ne. 0) then
				allocate(vmd_intervals(2,1))
				vmd_intervals(1,1) = 1; vmd_intervals(2,1) = Nsteps
				Nvmd_intervals = 1
			else
				allocate(vmd_intervals(2,Nvmd_intervals))
				!write(readin_format,'(a,i5,a)') '(',2*Nvmd_intervals,'i)'
				!read(1,trim(readin_format)) vmd_intervals
				read(1,*) vmd_intervals
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

	call locate(1,'SEPARATE_OUTFILES',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) separate_outfiles
	endif


	call locate(1,'BIN2CELLRATIO',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) binspercell(1)
		read(1,*) binspercell(2)
		read(1,*) binspercell(3)
		if (any(binspercell .gt. 1.d0) .and. & 
			any(binspercell .lt. 1.d0)) then
			print*, "WARNING -- BIN2CELLRATIO input - cannot specify multiple ", &
					"bins per cell in one direction and multiple cells per bin ", &
					"in the other -- setting minimum binspercell to 1"
			where(binspercell .gt. 1.d0) binspercell = 1.d0
		endif
	else
		binspercell = 1.d0
	endif
	call locate(1,'MACRO_OUTFLAG',.false.,found_in_input)
	if (found_in_input) read(1,*) macro_outflag
	call locate(1,'SLRC_FLAG',.false.,found_in_input)
	if (found_in_input) read(1,*) sLRC_flag
	!Test for WCA potential and switch LRC off
	if (abs(rcutoff-2.d0**(1.d0/6.d0)) .lt. 0.0001) then
		if (sLRC_flag .eq. 1) then
			print*, "WCA potential used - switching sLRC off"
			sLRC_flag = 0
		endif
	endif
	call locate(1,'MASS_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) mass_outflag
		if (mass_outflag .ne. 0) read(1,*) Nmass_ave
		if (mass_outflag .eq. 5) then
			call locate(1,'CPOL_BINS',.true.)
			read(1,*) gcpol_bins(1)	
			read(1,*) gcpol_bins(2)	
			read(1,*) gcpol_bins(3)	
		end if
	endif
	call locate(1,'VELOCITY_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) velocity_outflag
		if (velocity_outflag .ne. 0) read(1,* ) Nvel_ave
		if (velocity_outflag .eq. 5) then
			call locate(1,'CPOL_BINS',.true.)
			read(1,*) gcpol_bins(1)	
			read(1,*) gcpol_bins(2)	
			read(1,*) gcpol_bins(3)	
		end if
	endif
	call locate(1,'TEMPERATURE_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) temperature_outflag
		if (temperature_outflag .ne. 0)	then
			read(1,* ) NTemp_ave
			read(1,*,iostat=ios) peculiar_flag
			if (ios .ne. 0) peculiar_flag = 0 !default to zero if value not found
		endif
	endif
	call locate(1,'PRESSURE_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) pressure_outflag
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
		!if (pressure_outflag .eq. 2 .or. pressure_outflag .eq. 3) then
		!	VA_line_samples = 20
		!	call locate(1,'VA_LINE_SAMPLES',.false.,found_in_input)
		!	if (found_in_input) then
		!		read(1,*) VA_line_samples	
		!	end if
		!end if
	endif
	call locate(1,'VISCOSITY_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) viscosity_outflag
		if ( viscosity_outflag .ne. 0)	read(1,* ) Nvisc_ave
	endif
	call locate(1,'CV_CONSERVE',.false.,found_in_input)
	cv_conserve = 0
	if (found_in_input) then
		read(1,* ) cv_conserve
		read(1,*,iostat=ios) CV_debug
		if (ios .ne. 0) CV_debug = .false.
	endif
	call locate(1,'MFLUX_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) mflux_outflag
		if (mflux_outflag .ne. 0)	read(1,* ) Nmflux_ave
	endif
	call locate(1,'VFLUX_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) vflux_outflag
		if (vflux_outflag .ne. 0)	read(1,* ) Nvflux_ave
		if (mflux_outflag .eq. 0) Nmflux_ave = Nvflux_ave !Mass set to same as velocity
	endif
	call locate(1,'EFLUX_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) eflux_outflag
		if (eflux_outflag .ne. 0) then
			read(1,* ) Neflux_ave
			pass_vhalo = 1		!Turn on passing of velocities for halo images
		endif
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
	endif
	
	call locate(1,'STRUCT_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) ssf_outflag
		read(1,*) ssf_ax1
		read(1,*) ssf_ax2 
		read(1,*) ssf_nmax 
	endif

	call locate(1,'CV_FORCES',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) CVforce_flag
		if (CVforce_flag .ne. VOID) then
			if (.not. CV_debug) call error_abort("Input ERROR -- CV_FORCES true so CV_CONSERVE should be set to 1 and debugging set to .true.")
			if (vflux_outflag .ne. 4) call error_abort("Input ERROR -- CV_FORCES .true. but VFLUX_OUTFLAG not set to 4 (CV averages)")
		endif
		read(1,*,iostat=ios) CVforce_testcaseflag
		if (ios .ne. 0) CVforce_testcaseflag = 1
		read(1,*,iostat=ios) CVweighting_flag
		if (ios .ne. 0) CVweighting_flag = 0
	endif

	close(1,status='keep')      !Close input file

end subroutine setup_read_input
