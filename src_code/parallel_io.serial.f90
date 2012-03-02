!======================================================================
!	        	Serial Read Write Subroutines               

! Serial emulation of parallel i/o

! --- INPUT ROUTINES ---
! setup_command_arguments		Determine restart file (if any) and input file
! setup_inputs_locate			Read input file
! setup_restart_inputs			Read restart files and input file
! setup_restart_microstate		Read Initial configuration from restart file

! --- OUTPUT ROUTINES ---
! simulation_header				Write all variables in ; seperated variable form "Description; name; variable"
! update_simulation_progress    Record current iteration
! parallel_io_final_state 		Write final configuration and simulation details for restart
! parallel_io_vmd				Output for VMD
! parallel_io_vmd_sl			Output for VMD with seperate solid/lquid regions
! parallel_io_vmd_halo			Output Halos

! CV AVERAGING
! mass_slice_io					Write out mass bins for a single dimension slice through domain
! mass_bin_io					Write out mass bins in all 3 dimensions of domain
! velocity_slice_io				Write out velocity bins for a single dimension slice through domain
! velocity_bin_io				Write out velocity bins in all 3 dimensions of domain
!
! FLUX AVERAGING
! virial_stress_io				Write out virial stress
! VA_stress_io					Write out Volume Averaged stress throughout domain
! viscosity_io					Write out viscosity
! mass_flux_io					Write out flux of mass through bin surfaces
! momentum_flux_io				Write out flux of momnetum through bin surfaces
! MOP_stress_io					Write out stress on single plane through domain
! surface_stress_io				Write out stress on all surfaces of bins in domain
! 
!======================================================================

module module_parallel_io
	use interfaces
	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use polymer_info_MD
	use shear_info_MD
	use interfaces
end module

!======================================================================
!	        		INPUTS			              =
!======================================================================

!=============================================================================
!							command_arguments
! Checks for command-line arguments passed to the program and assigns the 
! relevant values.
!-----------------------------------------------------------------------------
subroutine setup_command_arguments
use module_parallel_io
implicit none
	
	integer				:: i,argcount
	logical 			:: restart_file_exists, input_file_exists
	character(len=200) 	:: arg,nextarg

	!Set default values in case they aren't specified by user
	restart = .false.			
	input_file_exists = .false.
	restart_file_exists = .false.
	input_file = 'MD.in'
	initial_microstate_file = 'results/final_state'

	argcount = command_argument_count()
	
	if (argcount.gt.0) then									!If more than 0 arguments	
			
		do i=1,argcount									!Loop through all arguments...
			
			call get_command_argument(i,arg)				!Reading two at once
            if (i < argcount) then  
				call get_command_argument(i+1,nextarg)
			else
				nextarg=''
			endif

			if (trim(arg).eq.'-r' .and. nextarg(1:1).ne.'-') then
				initial_microstate_file = trim(nextarg)
				inquire(file=initial_microstate_file, exist=restart_file_exists) 	!Check file exists
				if (restart_file_exists) restart = .true.                         
			end if			

			if (trim(arg).eq.'-i' .and. nextarg(1:1).ne.'-') then
				input_file = trim(nextarg)
			end if
		end do

	end if

	inquire(file=input_file, exist=input_file_exists)					!Check file exists
	if (input_file.eq.'MD.in'.and..not. input_file_exists) input_file = 'default.in'
	inquire(file=input_file, exist=input_file_exists)
	if(.not. input_file_exists) then
		print*, 'Input file ', trim(input_file), ' not found. Stopping simulation.'
		call error_abort
	end if

end subroutine setup_command_arguments

!-----------------------------------------------------------------------------
! Subroutine:	setup_inputs
! Author(s):	David Trevelyan & Edward Smith
! Description:
!		The input file contains capitalised keywords followed by
!		numerical values. The "locate" subroutine rewinds to the beginning
!		of the input file and scans each line until the keyword is matched.
!
!		Consequently, the file position is set for the next statement to read
!		the line underneath the previously "located" keyword. 
!-----------------------------------------------------------------------------
subroutine setup_inputs
	use module_parallel_io
	use librarymod, only : locate
	implicit none

	logical					:: found_in_input
	integer 				:: i, n , ios
	integer,dimension(8)	:: tvalue
	character(20)			:: readin_format

	!call random_seed
	call random_seed(size=n)
	allocate(seed(n))

	open(1,file=input_file)

	!Input physical co-efficients
	call locate(1,'DENSITY',.true.)
	read(1,*) density
	call locate(1,'RCUTOFF',.true.)
	read(1,*) rcutoff
	call locate(1,'INPUTTEMPERATURE',.true.)
	read(1,*) inputtemperature
	call locate(1,'INITIALNUNITS',.true.)
	read(1,*) initialnunits(1)		!x dimension split into number of cells
	read(1,*) initialnunits(2)		!y dimension split into number of cells
	read(1,*) initialnunits(3)		!z dimension split into number of cells
	call locate(1,'INTEGRATION_ALGORITHM',.true.)
	read(1,*) integration_algorithm
	call locate(1,'ENSEMBLE',.true.)
	read(1,*) ensemble
	call locate(1,'FORCE_LIST',.true.)	!AP, Cells, neighbrs..
	read(1,*) force_list
	if (force_list .ne. 3 .and. ensemble .eq. 4) &
		stop "Half int neighbour list only is compatible with pwa_terms_pwaNH thermostat"
	if (force_list .ne. 3 .and. ensemble .eq. 5) & 
		stop "Half int neighbour list only is compatible with DPD thermostat"
	call locate(1,'POTENTIAL_FLAG',.true.)	!LJ or FENE potential
	read(1,*) potential_flag
	if (potential_flag.eq.1) then
		call locate(1,'FENE_INFO',.true.)
		read(1,*) nmonomers
		read(1,*) k_c
		read(1,*) R_0
	end if	
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
	call locate(1,'SEED',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) seed(1) 	!Random number seed value 1
		read(1,*) seed(2) 	!Random number seed value 2
	else
		seed(1) = 1		!Fixed default seed for repeatability
		seed(2) = 2		!Fixed default seed for repeatability
	endif

	!Flags to determine if periodic boundaries are on	
	call locate(1,'PERIODIC',.true.)
	read(1,*) periodic(1)
	read(1,*) periodic(2)
	read(1,*) periodic(3)

	call locate(1,'DEFINE_SHEAR',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) shear_direction
		read(1,*) shear_iter0
		read(1,*) define_shear_as
		if (define_shear_as.eq.0) read(1,*) shear_velocity
		if (define_shear_as.eq.1) read(1,*) shear_rate
		if (define_shear_as.gt.1) then 
                	call error_abort( 'Poorly defined shear in input file')
        endif
	endif

	!-------------------------------------
	!Flag to determine molecular tags
	!-------------------------------------
	!Note: For initialunitsize "a"
	!		 [  o     o ]
	!a (1 cell size) [     o    ]  a/2 (distance between molcules)	
	!		 [  o     o
	!		  __________]  a/4 (distance from bottom of domain)
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

	call locate(1,'THERMSTAT_FLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) thermstat_flag
		select case(thermstat_flag)
		case(0)
			if (abs(maxval(thermstattop   )).ne.0.0 & 
		        .or.abs(maxval(thermstatbottom)).ne.0.0) call error_abort( & 
			 'THERMSTATTOP or THERMSTATBOTTOM non zero but THERMSTAT_FLAG_INFO &
                         & set to off (THERMSTAT_FLAG=0)')
			thermstatbottom = 0.d0; thermstattop = 0.d0 
		case(1)	!N-H thermostat all molecules
			thermstattop 	= initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
			thermstatbottom = initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
		case(2) !N-H PUT all molecules
			thermstattop 	= initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
			thermstatbottom = initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
		case(3)
			if (abs(maxval(thermstattop   )).eq.0.0 & 
		       .and.abs(maxval(thermstatbottom)).eq.0.0) & 
			call error_abort('THERMSTATTOP or THERMSTATBOTTOM must also be specified')
		end select
	endif

	!Flag to determine if output is switched on
	call locate(1,'VMD_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) vmd_outflag
		if (vmd_outflag .ne. 0) then
			read(1,*) Nvmd_intervals	!Number of vmd intervals
			if (Nvmd_intervals .gt. 20) then
				print*, 'Number of VMD intervals greater than 20 or not specified, setting on for all simualtion'
				Nvmd_intervals = 0
			endif
			if (Nvmd_intervals .eq. 0) then
				allocate(vmd_intervals(2,1))
				vmd_intervals(1,1) = 1; vmd_intervals(2,1) = huge(1)
			else
				allocate(vmd_intervals(2,Nvmd_intervals))
				!write(readin_format,'(a,i,a)') '(',2*Nvmd_intervals,'i)'
				!read(1,trim(readin_format)) vmd_intervals
                                read(1,*) vmd_intervals
#if USE_COUPLER
				!NEED SOME SORT OF coupler total simulation time retrival here!!
				print*, 'WARNING - CHECK VMD INTERVALS is not greater than coupled number of steps'
#else
				if (maxval(vmd_intervals) .gt. Nsteps) &
                                    call error_abort('Specified VMD interval greater than Nsteps')
#endif
			endif
		endif
	endif

	call locate(1,'MACRO_OUTFLAG',.false.,found_in_input)
	if (found_in_input) read(1,*) macro_outflag
	call locate(1,'MASS_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) mass_outflag
		if (mass_outflag .ne. 0) 	read(1,*) Nmass_ave
	endif
	call locate(1,'VELOCITY_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) velocity_outflag
		if (velocity_outflag .ne. 0)	read(1,* ) Nvel_ave
	endif
	call locate(1,'TEMPERATURE_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) temperature_outflag
		if (temperature_outflag .ne. 0)	read(1,* ) NTemp_ave
	endif
	call locate(1,'PRESSURE_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) pressure_outflag
		if (pressure_outflag .ne. 0) then
			read(1,* ) Nstress_ave
			read(1,*,iostat=ios) 	split_kin_config
			if (ios .ne. 0) split_kin_config = 0 !default to zero if value not found
		endif
	endif
	call locate(1,'VISCOSITY_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) viscosity_outflag
		if ( viscosity_outflag .ne. 0)	read(1,* ) Nvisc_ave
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
		if (eflux_outflag .ne. 0)	read(1,* ) Neflux_ave
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
	
	call locate(1,'R_GYRATION_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) r_gyration_outflag
		read(1,*) r_gyration_iter0
	end if

	close(1,status='keep')      !Close input file

	rcutoff2= rcutoff**2         !Useful definition to save computational time
	initialstep = 0   	     !Set initial step to one to start

	if (seed(1)==seed(2)) then
		call random_seed(get=seed(1:n))
                call date_and_time(values=tvalue)
                seed = IEOR(tvalue(8)+irank,seed)
        else 
                seed = irank
	endif
	
	!Assign seed to random number generator
	call random_seed(put=seed(1:n))

	elapsedtime = 1.d0*delta_t*Nsteps !Set elapsed time to end of simualtion

end subroutine setup_inputs

!----------------------------------------------------------------------
!
!                    Restart Simulation inputs
! Set up inputs based on the final state of a previous simulation
!
!----------------------------------------------------------------------

subroutine setup_restart_inputs
	use module_parallel_io
	use librarymod, only : locate
	implicit none

	logical					:: found_in_input
	integer					:: n, k, ios
	integer 				:: extrasteps
	integer 				:: checkint
	double precision 		:: checkdp
	integer,dimension(8)	:: tvalue
	character(20)			:: readin_format

 	integer(kind=selected_int_kind(18)) header_pos, end_pos ! 8 byte integer for header address

	!Allocate random number seed
	call random_seed(size=n)
	allocate(seed(n))

	!Open file to read integers
	open(2,file=initial_microstate_file,form='unformatted', access='stream',position='append')

    inquire(2,POS=end_pos) ! go the end of file
    read(2,pos=end_pos-8) header_pos ! header start is in the last 8 bytes
    header_pos = header_pos +1 ! for compatibility with MPI IO we store header_pos - 1 in final state 

	read(2,pos=header_pos) np
	globalnp = np			!Global np and local np same in serial
	read(2) initialnunits(1)
	read(2) initialnunits(2)
	read(2) initialnunits(3)
	read(2) Nsteps  	    
	read(2) tplot
	read(2) seed(1)
	read(2) seed(2)
	read(2) periodic(1)
	read(2) periodic(2)
	read(2) periodic(3)
	read(2) potential_flag
	read(2) nmonomers
	read(2) npx
	read(2) npy
	read(2) npz
	npx = 1
	npy = 1
	npz = 1

	read(2) density			!Density of system
	read(2) rcutoff			!Cut off distance for particle interaction
	rcutoff2= rcutoff**2         				!Useful definition to save computational time
	!read(2) inputtemperature	!Define initial temperature
	read(2) delta_t			!Timestep
	read(2) elapsedtime   		!Elapsed simulation time to date
	read(2)  k_c			!Polymer spring constant
	read(2)  R_0			!Polymer max bond elongation
	read(2)  delta_rneighbr		!Extra distance used for neighbour cell size

	close(2,status='keep')

	!Check if values from input file are different and alert user - all processors have
	!read the same file so only need to check on one processor

	open(1,file=input_file)

	call locate(1,'DENSITY',.true.)
	read(1,* ) checkdp                    !Density of system
	if (checkdp .ne. density)             print*, 'Discrepancy between system density', &
										  'in input & restart file - restart file will be used'
	call locate(1,'RCUTOFF',.true.)
	read(1,* ) checkdp                    !Cut off distance for particle interaction
	if (checkdp .ne. rcutoff)             print*, 'Discrepancy between cut off radius', &
										  'in input & restart file - restart file will be used'
	!call locate(1,'INPUTTEMPERATURE',.true.)
	!read(1,* ) checkdp                    !Define initial temperature
	!if (checkdp .ne. inputtemperature)    print*, 'Discrepancy between initial temperature', &
	!                                      'in input & restart file - restart file will be used'
	call locate(1,'INITIALNUNITS',.true.)
	read(1,* ) checkint                   !x dimension split into number of cells
	if (checkint .ne. initialnunits(1))   print*, 'Discrepancy between x domain size', &
										  'in input & restart file - restart file will be used'
	read(1,* ) checkint                   !y dimension box split into number of cells
	if (checkint .ne. initialnunits(2))   print*, 'Discrepancy between y domain size', &
										  'in input & restart file - restart file will be used'
	if (nd == 3) then	
	  read(1,* ) checkint                 !z dimension box split into number of cells
	  if (checkint .ne. initialnunits(3)) print*, 'Discrepancy between z domain size', &
										  'in input & restart file - restart file will be used'
	endif
	
	call locate(1,'POTENTIAL_FLAG',.true.)
	read(1,*) checkint                    !LJ or FENE potential
	if (checkint .ne. potential_flag)     print*, 'Discrepancy between potential_flag', &
										  'in input & restart file - restart file will be used'
	if (potential_flag.eq.1) then
	  call locate(1,'FENE_INFO',.true.) 
	  read(1,*) checkint                  !nmonomers - number of beads per chain
	  if (checkint.ne.nmonomers)          print*, 'Discrepancy between nmonomers', &
										  'in input & restart file - restart file will be used'
	  read(1,*) checkdp                   !k_c - FENE spring constant
	  if (checkint.ne.k_c)                print*, 'Discrepancy between k_c', &
										  'in input & restart file - restart file will be used'
	  read(1,*) checkdp                   !k_c - FENE spring constant
	  if (checkint.ne.R_0)                print*, 'Discrepancy between R_0', &
										  'in input & restart file - restart file will be used'
	end if	

	!Check periodic BC and shear
	call locate(1,'PERIODIC',.true.)
	read(1,*) periodic(1)
	read(1,*) periodic(2)
	read(1,*) periodic(3)

	!Get number of extra steps, timestep and plot frequency from input file	
	call locate(1,'NSTEPS',.true.)
	read(1,* ) extrasteps       !Number of computational steps
	call locate(1,'DELTA_T',.true.)
	read(1,* ) delta_t          !Size of time step
	call locate(1,'TPLOT',.true.)
	read(1,* ) tplot            !Frequency at which to record results
	call locate(1,'DELTA_RNEIGHBR',.true.)
	read(1,* ) delta_rneighbr   !Extra distance used for neighbour cell
	call locate(1,'SEED',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) seed(1) 	!Random number seed value 1
		read(1,*) seed(2) 	!Random number seed value 2
	else
		seed(1) = 1		!Fixed default seed for repeatability
		seed(2) = 2		!Fixed default seed for repeatability
	endif


	!Choose integration algorithm
	call locate(1,'INTEGRATION_ALGORITHM',.true.)
	read(1,*) integration_algorithm
	
	call locate(1,'ENSEMBLE',.true.)
	read(1,*) ensemble

	call locate(1,'FORCE_LIST',.true.)	!LJ or FENE potential
	read(1,*) force_list
	
	call locate(1,'INPUTTEMPERATURE',.true.)
	read(1,* ) inputtemperature         !Define initial temperature
	

	call locate(1,'DEFINE_SHEAR',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) shear_direction
		read(1,*) shear_iter0
		read(1,*) define_shear_as
		if (define_shear_as.eq.0) read(1,*) shear_velocity
		if (define_shear_as.eq.1) read(1,*) shear_rate
		if (define_shear_as.gt.1) then 
                call error_abort( 'Poorly defined shear in input file')
        endif
	endif

	!-------------------------------------
	!Flag to determine molecular tags
	!-------------------------------------
	!Note: For initialunitsize 'a'
	!		 		 [  o     o ]
	!a (1 cell size) [     o    ]  a/2 (distance between molcules)	
	!		 		 [  o     o ]
	!		  		 [__________]  a/4 (distance from bottom of domain)
	!
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

	call locate(1,'THERMSTAT_FLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) thermstat_flag
		select case(thermstat_flag)
		case(0)
			if (abs(maxval(thermstattop   )).ne.0.0 & 
		        .or.abs(maxval(thermstatbottom)).ne.0.0) call error_abort( & 
			'THERMSTATTOP or THERMSTATBOTTOM non zero but THERMSTAT_FLAG_INFO &
                        & set to off (THERMSTAT_FLAG=0)')
			thermstatbottom = 0.d0; thermstattop = 0.d0 
		case(1)
			thermstattop 	= initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
			thermstatbottom = initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
		case(2)
			if (abs(maxval(thermstattop   )).eq.0.0 & 
		       .and.abs(maxval(thermstatbottom)).eq.0.0)  call error_abort( &
			 'THERMSTATTOP or THERMSTATBOTTOM must also be specified')
		end select
	endif


	!Flag to determine if output is switched on
	call locate(1,'VMD_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) vmd_outflag
		if (vmd_outflag .ne. 0) then
			read(1,*) Nvmd_intervals	!Number of vmd intervals
			if (Nvmd_intervals .gt. 20) then
				print*, 'Number of VMD intervals greater than 20 or not specified, setting on for all simualtion'
				Nvmd_intervals = 0
			endif
			if (Nvmd_intervals .eq. 0) then
				allocate(vmd_intervals(2,1))
				vmd_intervals(1,1) = 1; vmd_intervals(2,1) = huge(1)
			else
				allocate(vmd_intervals(2,Nvmd_intervals))
				!write(readin_format,'(a,i,a)') '(',2*Nvmd_intervals,'i)'
				!read(1,trim(readin_format)) vmd_intervals
                                read(1,*) vmd_intervals
#if USE_COUPLER
				!NEED SOME SORT OF coupler total simulation time retrival here!!
				print*, 'WARNING - CHECK VMD INTERVALS is not greater than coupled number of steps'
#else
				if (maxval(vmd_intervals) .gt. Nsteps) &
                                    call error_abort('Specified VMD interval greater than Nsteps')
#endif
			endif
		endif
	endif

	call locate(1,'MACRO_OUTFLAG',.false.,found_in_input)
	if (found_in_input) read(1,*) macro_outflag
	call locate(1,'MASS_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) mass_outflag
		if (mass_outflag .ne. 0) 	read(1,*) Nmass_ave
	endif
	call locate(1,'VELOCITY_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) velocity_outflag
		if (velocity_outflag .ne. 0)	read(1,* ) Nvel_ave
	endif
	call locate(1,'PRESSURE_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) pressure_outflag
		if (pressure_outflag .ne. 0) then
			read(1,* ) Nstress_ave
			read(1,*,iostat=ios) 	split_kin_config
			if (ios .ne. 0) split_kin_config = 0 !default to zero if value not found
		endif
	endif
	call locate(1,'VISCOSITY_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,* ) viscosity_outflag
		if ( viscosity_outflag .ne. 0)	read(1,* ) Nvisc_ave
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
	endif

	elapsedtime = elapsedtime + delta_t*extrasteps !Set elapsed time to end of simualtion
	initialstep = Nsteps         !Set plot count to final plot of last
	iter = initialstep			 !Set iter to initialstep so that initial record is performed correctly at restart
	Nsteps = Nsteps + extrasteps !Establish final iteration step based on previous

	call locate(1,'ETEVTCF_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) etevtcf_outflag
		if (etevtcf_outflag.ne.0) then
			read(1,*) etevtcf_iter0
			if (etevtcf_iter0.lt.iter) then
				print*, 'Restart functionality has not yet been implemented for time correlation functions.', &
					'iter_0 for calculation of the end-to-end vector time correlation function has been reset to ', iter
				etevtcf_iter0 = iter	
			end if
			if (mod(etevtcf_iter0,tplot).ne.0) then
				etevtcf_iter0 = etevtcf_iter0 + (tplot - mod(etevtcf_iter0,tplot))
				print*, 'Etevtcf must be a multiple of tplot, resetting etevtcf to ', etevtcf_iter0
			end if
		end if
	endif
	
	call locate(1,'R_GYRATION_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) r_gyration_outflag
		read(1,*) r_gyration_iter0
	end if

	close(1,status='keep')      !Close input file

end subroutine setup_restart_inputs


!----------------------------------------------------------------------------------
!
!                                Restart Microstate
! Set up position and velocity of each molecule (the microstate) based on the final
! state of a previous simulation
!
!---------------------------------------------------------------------------------

subroutine setup_restart_microstate
	use module_parallel_io
	implicit none

	integer 				:: ixyz, n

	!Open file at first recorded value
	open(2,file=initial_microstate_file, form='unformatted', access='stream',position='rewind')
	!Read positions
	select case (potential_flag)
	case(0)
		do n=1,globalnp
			!Read position from file
			do ixyz=1,nd
				read(2) r(n,ixyz)
			enddo
			do ixyz=1,nd
				read(2) v(n,ixyz)
			enddo
		enddo
	case(1)
		do n=1,globalnp
			read(2) r(n,:)
			read(2) v(n,:)
			read(2) monomer(n)%chainID
			read(2) monomer(n)%subchainID
			read(2) monomer(n)%funcy
			read(2) monomer(n)%glob_no
			read(2) monomer(n)%bondflag(1:nmonomers)
		end do
	case default
	end select

!	select case (integration_algorithm)
!		case(leap_frog_verlet)
!			call simulation_compute_forces
!			do n=1,globalnp
!				v(n,:) = v(n,:) - 0.5d0*delta_t*a(n,:)		
!			end do
!		case(velocity_verlet)
!			!Nothing
!	end select

	call setup_tag				!Setup location of fixed molecules
	do n = 1,np
		call read_tag(n)		!Read tag and assign properties
	enddo

	close(2,status='keep') 		!Close final state file

end subroutine setup_restart_microstate

!======================================================================
!	        		OUTPUTS			              =
!======================================================================

!---------------------------------------------------------------------------------
! Write Simulation Parameter File contain all data required to completely recreate
! simulation and to be used for post processing
!
subroutine simulation_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	Character(8)  		:: the_date
	Character(10)  		:: the_time

	call date_and_time(the_date, the_time)

	open(3,file=trim(prefix_dir)//'results/simulation_header')

	write(3,*) 'Simulation run on Date;  sim_date ;', the_date
	write(3,*) 'Simulation start time ;  sim_start_time ;', the_time
	write(3,*) 'Number of Dimensions ;  nd ;', nd
	write(3,*) 'Number of Particles ;  globalnp ;', globalnp
	write(3,*) 'Time Step - delta t ;   delta_t ;',  delta_t
	write(3,*) 'Total number of steps ; Nsteps;',  Nsteps - initialstep
	write(3,*) 'Integration algorithm ; integration_algorithm;', integration_algorithm
	select case(potential_flag)
	case(0)
		write(3,*) 'Potential flag ; potential_flag;', potential_flag
	case(1)
		write(3,*) 'Potential flag ; potential_flag;', potential_flag
		write(3,*) 'Number of LJ beads per FENE chain ; nmonomers;', nmonomers
		write(3,*) 'Number of FENE chains in domain ; nchains;', nchains
		write(3,*) 'FENE bond maximum elongation ; R_0;', R_0
		write(3,*) 'FENE spring stiffness ; k_c;', k_c
	end select	
	write(3,*) 'Starting step of simulation ;  initialstep ;', initialstep
	write(3,*) 'Generate output file every steps ;   tplot ;',  tplot
	write(3,*) 'Density ; density ;',density
	write(3,*) 'Initial Temperature ;   inputtemperature ;',  inputtemperature
	write(3,*) 'Cut off distance ;  rcutoff ;', rcutoff
	write(3,*) 'Neighbour List Delta r ;  delta_rneighbr ;', delta_rneighbr
	write(3,*) 'Initial FCC unit size in x ;  initialunitsize(1) ;', initialunitsize(1)
	write(3,*) 'Initial FCC unit size in y ;  initialunitsize(2) ;', initialunitsize(2)
	write(3,*) 'Initial FCC unit size in z ;  initialunitsize(3) ;', initialunitsize(3)
	write(3,*) 'Domain in x ;  globaldomain(1)  ;', globaldomain(1) 
	write(3,*) 'Domain in y ;  globaldomain(2)  ;', globaldomain(2) 
	write(3,*) 'Domain in z ;  globaldomain(3)  ;', globaldomain(3) 
	write(3,*) 'Domain volume ;  volume ;', volume
	write(3,*) 'Periodic Boundary Conditions in x ;  periodic(1) ;', periodic(1)
	write(3,*) 'Periodic Boundary Conditions in y ;  periodic(2) ;', periodic(2)
	write(3,*) 'Periodic Boundary Conditions in z ;  periodic(3) ;', periodic(3)
	write(3,*) 'Dist frm bot Fixed Mol in x; fixdistbot(1);', fixdistbottom(1)
	write(3,*) 'Dist frm bot Fixed Mol in y; fixdistbot(2);', fixdistbottom(2)
	write(3,*) 'Dist frm bot Fixed Mol in z; fixdistbot(3);', fixdistbottom(3)
	write(3,*) 'Dist frm top Fixed Mol in x; fixdisttop(1);', fixdisttop(1)
	write(3,*) 'Dist frm top Fixed Mol in y; fixdisttop(2);', fixdisttop(2)
	write(3,*) 'Dist frm top Fixed Mol in z; fixdisttop(3);', fixdisttop(3)
	write(3,*) 'Dist frm bot Tethered Mol in x; tethdistbot(1);', tethereddistbottom(1)
	write(3,*) 'Dist frm bot Tethered Mol in y; tethdistbot(2);', tethereddistbottom(2)
	write(3,*) 'Dist frm bot Tethered Mol in z; tethdistbot(3);', tethereddistbottom(3)
	write(3,*) 'Dist frm top Tethered Mol in x; tethdisttop(1);', tethereddisttop(1)
	write(3,*) 'Dist frm top Tethered Mol in y; tethdisttop(2);', tethereddisttop(2)
	write(3,*) 'Dist frm top Tethered Mol in z; tethdisttop(3);', tethereddisttop(3)
	write(3,*) 'Dist frm bot Sliding Mol in x; slidedistbot(1);', slidedistbottom(1)
	write(3,*) 'Dist frm bot Sliding Mol in y; slidedistbot(2);', slidedistbottom(2)
	write(3,*) 'Dist frm bot Sliding Mol in z; slidedistbot(3);', slidedistbottom(3)
	write(3,*) 'Dist frm top Sliding Mol in x; slidedisttop(1);', slidedisttop(1)
	write(3,*) 'Dist frm top Sliding Mol in y; slidedisttop(2);', slidedisttop(2)
	write(3,*) 'Dist frm top Sliding Mol in z; slidedisttop(3);', slidedisttop(3)
	write(3,*) 'Sliding velocity of wall in x; wallslidev(1);', wallslidev(1)
	write(3,*) 'Sliding velocity of wall in y; wallslidev(2);', wallslidev(2)
	write(3,*) 'Sliding velocity of wall in z; wallslidev(3);', wallslidev(3)
	write(3,*) 'Dist frm bot NH Thermstat Mol in x; thermstatbot(1);', thermstatbottom(1)
	write(3,*) 'Dist frm bot NH Thermstat Mol in y; thermstatbot(2);', thermstatbottom(2)
	write(3,*) 'Dist frm bot NH Thermstat Mol in z; thermstatbot(3);', thermstatbottom(3)
	write(3,*) 'Dist frm top NH Thermstat Mol in x; thermstattop(1);', thermstattop(1)
	write(3,*) 'Dist frm top NH Thermstat Mol in y; thermstattop(2);', thermstattop(2)
	write(3,*) 'Dist frm top NH Thermstat Mol in z; thermstattop(3);', thermstattop(3)
	write(3,*) 'Computational cells in x ;  globalncells(1) ;',  ncells(1)*npx
	write(3,*) 'Computational cells in y ;  globalncells(2)  ;', ncells(2)*npy 
	write(3,*) 'Computational cells in z ;  globalncells(3)  ;', ncells(3)*npz 
	write(3,*) 'Of size in x ;  cellsidelength(1) ;', cellsidelength(1)
	write(3,*) 'Of size in y ;  cellsidelength(2) ;', cellsidelength(2)
	write(3,*) 'Of size in z ;  cellsidelength(3) ;', cellsidelength(3)
	write(3,*) 'Number of processors in x ;  npx ;', npx
	write(3,*) 'Number of processors in y ;  npy ;', npy
	write(3,*) 'Number of processors in z ;  npz ;', npz
	write(3,*) 'Cells per Processor in x ;  nicellxl ;', nicellxl
	write(3,*) 'Cells per Processor in y ;  nicellyl ;', nicellyl
	write(3,*) 'Cells per Processor in z ;  nicellzl ;', nicellzl
	write(3,*) 'Cells per Processor including Halos in x ;  ncellxl ;', ncellxl
	write(3,*) 'Cells per Processor including Halos in y ;  ncellyl ;', ncellyl
	write(3,*) 'Cells per Processor including Halos in z ;  ncellzl ;', ncellzl
	write(3,*) '1st Random seed ;  seed_1 ;', seed(1)
	write(3,*) '2nd Random seed ;  seed_2 ;', seed(2)
	write(3,*)  'VMD flag ;  vmd_outflag ;', vmd_outflag
	write(3,*)  'macro flag ;  macro_outflag	 ;', macro_outflag
	write(3,*)  'mass flag ;  mass_outflag ;', mass_outflag	
	write(3,*)  'velocity flag ;  velocity_outflag ;', velocity_outflag
	write(3,*)  'temperature flag ;  temperature_outflag ;', temperature_outflag
	write(3,*)  'Pressure flag ;  pressure_outflag ;', pressure_outflag
	write(3,*)  'viscosity flag ;  viscosity_outflag ;', viscosity_outflag
	write(3,*)  'mass flux flag ;  mflux_outflag ;', mflux_outflag
	write(3,*)  'velocity flux flag ;  vflux_outflag ;', vflux_outflag
	write(3,*)  'mass average steps ;  Nmass_ave ;', Nmass_ave
	write(3,*)  'velocity average steps ;  Nvel_ave ;', Nvel_ave
	write(3,*)  'Temperature average steps ;  NTemp_ave ;', NTemp_ave
	write(3,*)  'pressure average steps ;  Nstress_ave ;', Nstress_ave
	write(3,*)  'viscosity average samples ;  Nvisc_ave ;', Nvisc_ave
	write(3,*)  'mass flux average steps ;  Nmflux_ave ;', Nmflux_ave
	write(3,*)  'velocity flux average steps ;  Nvflux_ave ;', Nvflux_ave
	write(3,*)  'Velocity/stress Averaging Bins in x ;  globalnbins(1) ;', globalnbins(1)
	write(3,*)  'Velocity/stress Averaging Bins in y ;  globalnbins(2) ;', globalnbins(2)
	write(3,*)  'Velocity/stress Averaging Bins in z ;  globalnbins(3) ;', globalnbins(3)
	write(3,*)  'Of size in x ;  binsize(1)  ;', globaldomain(1)/globalnbins(1) 
	write(3,*)  'Of size in y ;  binsize(2)  ;', globaldomain(2)/globalnbins(2) 
	write(3,*)  'Of size in z ;  binsize(3)  ;', globaldomain(3)/globalnbins(3) 
	write(3,*)  'Bins per Processor in x ;  nbins(1) ;', nbins(1)
	write(3,*)  'Bins per Processor in y ;  nbins(2) ;', nbins(2)
	write(3,*)  'Bins per Processor in z ;  nbins(3) ;', nbins(3)
	write(3,*)  'Number of Bins on outer Surface of each processor ;  nsurfacebins ;', nsurfacebins
	write(3,*)  'Number of Bins in halo of each processor ;  nhalobins ;', nhalobins
	write(3,*)  'Domain split into Planes for Pressure Averaging ; nplanes  ;',nplanes 
	write(3,*)  'Separated by distance ;  planespacing  ;', planespacing 
	write(3,*)  'with first plane at ;  planes ;', planes(1)
	write(3,*)	'Shear direction ; shear_direction;', shear_direction
	write(3,*)  'Leapfrog or Velocity-Verlet ; integration_algorithm ;', integration_algorithm
	write(3,*)  'Force calculation list methodd ; force_list ;', force_list
	!write(3,*)  'Ensemble; ensemble; ', ensemble		!MATLAB input functions can't deal with words...
	write(3,*)	'Shear direction ; shear_direction;', shear_direction

	close(3,status='keep')

end subroutine simulation_header

!------------------------------------------------------------------------------
!Serial version of parallel code to print final_state for restart

subroutine parallel_io_final_state
	use module_parallel_io
	implicit none

	integer 								:: ixyz,n
	integer, parameter                      :: zero=0
	integer, dimension(np) 					:: chainID, subchainID,right,left
	integer 								:: int_filesize,dp_filesize
	integer(kind=selected_int_kind(18))		:: header_pos ! 8 byte integer for header address

	!Rebuild simulation before recording final state
	call sendmols			   			!Exchange particles between processors
	call linklist_deallocateall	   		!Deallocate all linklist components
	call assign_to_cell	  	   			!Re-build linklist every timestep
	call messenger_updateborders(1)	   	!Update borders between processors
	call assign_to_neighbourlist	   	!Setup neighbourlist

	!Remove previous final state file
	open(2,file=trim(prefix_dir)//'results/final_state')
	close(2,status='delete')
	!open(3,file=trim(prefix_dir)//'results/finalvelocities',status='replace')

!	select case (integration_algorithm)
!		case(leap_frog_verlet)
!			call simulation_compute_forces
!			do n=1,globalnp
!				v(n,:) = v(n,:) + 0.5d0*delta_t*a(n,:)		
!			end do
!		case(velocity_verlet)
!			!Nothing
!	end select
	
	!Written in this form so each molecule's information is together to allow 
	!re-allocation to seperate processors
	open(2,file=trim(prefix_dir)//'results/final_state', form='unformatted',access='stream')
	select case (potential_flag)
	case(0)
		do n=1,np
			!Write particle n's positions and speed
			do ixyz=1,nd
				write(2) r(n,ixyz) 
			enddo
			!Write particle n's velocities
			do ixyz=1,nd
				write(2) v(n,ixyz)   
			enddo
		enddo
	case(1)
		do n=1,np
			write(2) r(n,:)
			write(2) v(n,:)
			write(2) monomer(n)%chainID
			write(2) monomer(n)%subchainID
			write(2) monomer(n)%funcy
			write(2) monomer(n)%glob_no
			write(2) monomer(n)%bondflag(1:nmonomers)
		end do
	case default
	end select	
 
	close(2,status='keep')

	!Write integer data at end of file	
	open(2,file=trim(prefix_dir)//'results/final_state', form='unformatted',access='stream',position='append')

	! header address
	inquire(2,POS=header_pos)
        
	write(2) np               	!Number of particles
	write(2) initialnunits(1) 	!x dimension split into number of cells
	write(2) initialnunits(2) 	!y dimension box split into number of cells
	write(2) initialnunits(3) 	!z dimension box split into number of cells
	write(2) Nsteps           	!Number of computational steps
	write(2) tplot            	!Frequency at which to record results
	write(2) seed(1)          	!Random number seed value 1
	write(2) seed(2)          	!Random number seed value 2
	write(2) periodic(1)  	 	!Boundary condition flags
	write(2) periodic(2)   		!Boundary condition flags
	write(2) periodic(3)	   	!Boundary condition flags
	write(2) potential_flag   	!Polymer/LJ potential flag
	write(2) nmonomers		   	!Polymer chain length
	write(2) 0                  !Processors (npx) flag for serial
	write(2) 0                  !Processors (npy) flag for serial
	write(2) 0                  !Processors (npz) flag for serial

	write(2) density          !Density of system
	write(2) rcutoff          !Cut off distance for particle interaction
	write(2) delta_t          !Size of time step
	write(2) elapsedtime      !Total elapsed time of all restarted simulations
	write(2) k_c		   	  !FENE spring constant
	write(2) R_0		      !FENE spring max elongation
	write(2) delta_rneighbr	  !Extra distance used for neighbour list cell size
    write(2) header_pos-1     ! -1 for MPI IO compatibility
	close(2,status='keep') 	  !Close final_state file

end subroutine parallel_io_final_state

subroutine parallel_io_final_state_old
	use module_parallel_io
	implicit none

	integer :: ixyz,n

	!Rebuild simulation before recording final state
	call sendmols			   !Exchange particles between processors
	call linklist_deallocateall	   !Deallocate all linklist components
	call assign_to_cell	  	   !Re-build linklist every timestep
	call messenger_updateborders(1)	   !Update borders between processor
	call assign_to_neighbourlist	   !Setup neighbourlist

	!Remove previous final state file
	open(2,file=trim(prefix_dir)//'results/final_state')
	close(2,status='delete')
	!open(3,file=trim(prefix_dir)//'results/finalvelocities',status='replace')

	open(2,file=trim(prefix_dir)//'results/final_state', form='unformatted',access='direct',recl=2)
	!Written in this form so each molecule's information is together to allow 
	!re-allocation to seperate processors
	do n=1,np
		!Write particle n's positions
		do ixyz=1,nd
			write(2,rec=6*(n-1)+(ixyz-1)+1) r(n,ixyz)	
		enddo
		!Write particle n's velocities
		do ixyz=1,nd
			write(2,rec=6*(n-1)+(ixyz-1)+4) v(n,ixyz)   
	!		write(3,*) v(n,ixyz)
		enddo

	enddo

	!close(3)
	
	write(2,rec=6*np+1) density           !Density of system
	write(2,rec=6*np+2) rcutoff           !Cut off distance for particle interaction
	write(2,rec=6*np+3) inputtemperature  !Define initial temperature
	write(2,rec=6*np+4) delta_t           !Size of time step
	write(2,rec=6*np+5) elapsedtime       !Total elapsed time of all restarted simulations
	write(2,rec=6*np+6) k_c			   	  !FENE spring constant
	write(2,rec=6*np+7) R_0			      !FENE spring max elongation
	write(2,rec=6*np+8) delta_rneighbr	  !Extra distance used for neighbour list cell size

	close(2,status='keep') !Close final_state file

	!Re-open file with different record length to write Integer Data
	open(2,file=trim(prefix_dir)//'results/final_state', form='unformatted',access='direct',recl=1)

	write(2,rec=12*np+17) np               !Number of particles
	write(2,rec=12*np+18) initialnunits(1) !x dimension split into number of cells
	write(2,rec=12*np+19) initialnunits(2) !y dimension box split into number of cells
	write(2,rec=12*np+20) initialnunits(3) !z dimension box split into number of cells
	write(2,rec=12*np+21) Nsteps           !Number of computational steps
	write(2,rec=12*np+22) tplot            !Frequency at which to record results
	write(2,rec=12*np+23) seed(1)          !Random number seed value 1
	write(2,rec=12*np+24) seed(2)          !Random number seed value 2
	write(2,rec=12*np+25) periodic(1)	   !Boundary condition flags
	write(2,rec=12*np+26) periodic(2)	   !Boundary condition flags
	write(2,rec=12*np+27) periodic(3)	   !Boundary condition flags
	write(2,rec=12*np+28) potential_flag   !Polymer/LJ potential flag
	write(2,rec=12*np+29) nmonomers	   !Polymer chain length
	write(2,rec=12*np+30) 0				   !Dummy to make even filesize

	close(2,status='keep') !Close final_state file

end subroutine parallel_io_final_state_old

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd(start, finish,interval_no)
	use module_parallel_io
	implicit none

	integer,intent(in)			:: start, finish,interval_no
	integer						:: i, n, length
	integer                     :: molno
	real,dimension(3*np)		:: buf
	real,dimension(np)          :: Xbuf, Ybuf, Zbuf

	select case (potential_flag)
	case(0)
		buf(1     :  np) = r(:,1)
		buf(np+1  :2*np) = r(:,2)
		buf(2*np+1:3*np) = r(:,3)
	case(1)
		do n=1,np
			molno       = monomer(n)%glob_no
			Xbuf(molno) = r(n,1)
			Ybuf(molno) = r(n,2)
			Zbuf(molno) = r(n,3)
		end do
		buf(1     :np  ) = Xbuf
		buf(np+1  :2*np) = Ybuf
		buf(2*np+1:3*np) = Zbuf
	case default
	end select		

	!If intervals set to zero then full simulation recorded
	if (Nvmd_intervals.eq.0) then
		i = (iter-initialstep)/tplot+1
	else
		!Otherwise, calculate number of previous intervals
		i = vmd_count
		!do n=1,interval_no-1
		!	i = i + (vmd_intervals(2,n)-vmd_intervals(1,n))/tplot
		!enddo
		!i = i + ((iter-start)/tplot)+1
	endif

	inquire(iolength=length) buf
	open (unit=4, file=trim(prefix_dir)//'results/vmd_temp.dcd',access='direct',recl=length)
	write(4,rec=i) buf(:)
	close(4,status='keep')

end subroutine parallel_io_vmd

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd_sl(start, finish,interval_no)
	use module_parallel_io
	implicit none

	integer,intent(in)		:: start, finish,interval_no
	integer					:: i, n
	real,dimension(np)		:: Xbuf, Ybuf, Zbuf
	real,dimension(3)		:: rhalfdomain

	Xbuf(:) = r(:,1)
	Ybuf(:) = r(:,2)
	Zbuf(:) = r(:,3)

	rhalfdomain(1) = -halfdomain(1)
	rhalfdomain(2) = -halfdomain(2)
	rhalfdomain(3) = -halfdomain(3)

	!If number of intervals is equal to zero then run a full simulation recorded
	if (Nvmd_intervals.eq.0) then
		i = (iter-initialstep)/tplot+1
	else
		!Calculate number of previous intervals and start writing from here
		i = vmd_count
		!i = 0
		!do n=1,interval_no-1
		!	i = i + (vmd_intervals(2,n)-vmd_intervals(1,n))/tplot
		!enddo
		!i = i + ((iter-start)/tplot)+1
	endif

	!---Write liquid molecules---
	open (unit=4, file=trim(prefix_dir)//'results/vmd_liquid_temp.dcd',access='direct',recl=1)
	
	do n = 1,np
		select case (tag(n))
		case(0, 4) 	!Liquid Molecules
			write(4,rec=(i-1)*nd*np+n) Xbuf(n)
			write(4,rec=(i-1)*nd*np+n+np) Ybuf(n)
			write(4,rec=(i-1)*nd*np+n+2*np) Zbuf(n)
		case default	!Solid molecules
			write(4,rec=(i-1)*nd*np+n) rhalfdomain(1)
			write(4,rec=(i-1)*nd*np+n+np) rhalfdomain(2)
			write(4,rec=(i-1)*nd*np+n+2*np) rhalfdomain(3)
		end select
	enddo

	close(4,status='keep')

	!---Write solid molecules---
	open (unit=4, file=trim(prefix_dir)//'results/vmd_solid_temp.dcd',access='direct',recl=1)

	do n = 1,np
		select case (tag(n))
		case(0, 4) 	!Liquid Molecules
			write(4,rec=(i-1)*nd*np+n) rhalfdomain(1)
			write(4,rec=(i-1)*nd*np+n+np) rhalfdomain(2)
			write(4,rec=(i-1)*nd*np+n+2*np) rhalfdomain(3)
		case default	!Solid molecules
			write(4,rec=(i-1)*nd*np+n) Xbuf(n)
			write(4,rec=(i-1)*nd*np+n+np) Ybuf(n)
			write(4,rec=(i-1)*nd*np+n+2*np) Zbuf(n)
		end select
	enddo

	close(4,status='keep')

end subroutine parallel_io_vmd_sl

!------------------------------------------------------------------------
!Write positions of molecules in halo to a file

subroutine parallel_io_vmd_halo(start, finish,interval_no)
	use module_parallel_io
	implicit none

	integer,intent(in)			:: start, finish,interval_no
	integer						:: i, n
	real,dimension(halo_np)		:: Xbuf, Ybuf, Zbuf
	real,dimension(3)			:: rhalfdomain

	Xbuf(:) = r(np+1:np+halo_np,1)
	Ybuf(:) = r(np+1:np+halo_np,2)
	Zbuf(:) = r(np+1:np+halo_np,3)

	rhalfdomain(1) = -halfdomain(1)
	rhalfdomain(2) = -halfdomain(2)
	rhalfdomain(3) = -halfdomain(3)

	!If finish eq zero then full simulation recorded
	if (Nvmd_intervals.eq.0) then
		i = (iter-initialstep)/tplot+1
	else
		!Calculate number of previous intervals
		i = vmd_count
		!do n=1,interval_no-1
		!	i = i + (vmd_intervals(2,n)-vmd_intervals(1,n))/tplot
		!enddo
		!i = i + ((iter-start)/tplot)+1
	endif

	open (unit=4, file=trim(prefix_dir)//'results/vmd_halo_temp.dcd',access='direct',recl=1)
	
	do n = 1,halo_np
		write(4,rec=(i-1)*nd*extralloc+n) Xbuf(n)
		write(4,rec=(i-1)*nd*extralloc+n+extralloc) Ybuf(n)
		write(4,rec=(i-1)*nd*extralloc+n+2*extralloc) Zbuf(n)
	enddo

	!Write up to extralloc to give constant number of molecule in vmd file
	do n = halo_np,halo_np+extralloc
		write(4,rec=(i-1)*nd*extralloc+n+3*extralloc) rhalfdomain(1)
		write(4,rec=(i-1)*nd*extralloc+n+4*extralloc) rhalfdomain(2)
		write(4,rec=(i-1)*nd*extralloc+n+5*extralloc) rhalfdomain(3)
	enddo

	close(4,status='keep')

end subroutine parallel_io_vmd_halo

!---------------------------------------------------------------------------------
! Write value of last output iteration

subroutine update_simulation_progress_file
	use module_parallel_io
	implicit none

	open (unit=99999, file=trim(prefix_dir)//'results/simulation_progress')
	write(99999,*) iter
	close(99999,status='keep')

end subroutine update_simulation_progress_file

!---------------------------------------------------------------------------------
! Record mass in a slice through the domain

subroutine mass_slice_io(ixyz)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer		:: ixyz, m, length

	!Write mass slice to file
	m = iter/(tplot*Nmass_ave)
	inquire(iolength=length) slice_mass(1:nbins(ixyz))
	open (unit=5, file=trim(prefix_dir)//'results/mslice',form='unformatted',access='direct',recl=length)
	write(5,rec=m) slice_mass(1:nbins(ixyz))
	close(5,status='keep')

end subroutine mass_slice_io

!---------------------------------------------------------------------------------
! Record mass in 3D bins throughout domain

subroutine mass_bin_io(CV_mass_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: i,j,k,n
	integer				:: m,length
	integer				:: CV_mass_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	integer				:: buf(nbins(1),nbins(2),nbins(3))
	character(4)		:: io_type
	character(13)		:: filename

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) 'results/m', io_type

	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nhalobins
		i = halobins(n,1); j = halobins(n,2); k = halobins(n,3)  
		!Change in number of Molecules in halo cells
		CV_mass_out(modulo((i-2),nbins(1))+2, & 
			    	modulo((j-2),nbins(2))+2, & 
			    	modulo((k-2),nbins(3))+2) = & 
			    CV_mass_out(modulo((i-2),nbins(1))+2, & 
							modulo((j-2),nbins(2))+2, & 
							modulo((k-2),nbins(3))+2) &
									+ CV_mass_out(i,j,k)
	enddo

	if (io_type .eq. 'snap') then
		m = iter/(Nmflux_ave) + 1 !Initial snapshot taken
	else
		m = iter/(tplot*Nmass_ave)
	endif
	!Write mass to file
	buf = CV_mass_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1)
	inquire(iolength=length) buf
	open (unit=5,file=filename,form='unformatted',access='direct',recl=length)
	write(5,rec=m) buf
	close(5,status='keep')

end subroutine mass_bin_io

!---------------------------------------------------------------------------------
! Record velocity in a slice through the domain

subroutine velocity_slice_io(ixyz)
use module_parallel_io
use calculated_properties_MD
implicit none

	integer		:: ixyz, m, length

	!Write mass
	call mass_slice_io(ixyz)

	!Write velocity to file
	m = iter/(tplot*Nvel_ave)
	inquire(iolength=length) slice_momentum(1:nbins(ixyz),:)
	open (unit=6, file=trim(prefix_dir)//'results/vslice',form='unformatted',access='direct',recl=length)
	write(6,rec=m) slice_momentum(1:nbins(ixyz),:)
	close(6,status='keep')

end subroutine velocity_slice_io

!---------------------------------------------------------------------------------
! Record velocity in 3D bins throughout domain

subroutine velocity_bin_io(CV_mass_out,CV_momentum_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer					:: n,m,i,j,k
	integer					:: length
	integer					:: CV_mass_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	double precision		:: CV_momentum_out(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)
	double precision		:: buf(nbins(1),nbins(2),nbins(3),3)
	character(4)			:: io_type
	character(13)			:: filename

	!Write mass bins
	call mass_bin_io(CV_mass_out,io_type)

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) 'results/v', io_type

	!---------------Correct for surface fluxes on halo cells---------------
	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nhalobins
		i = halobins(n,1); j = halobins(n,2); k = halobins(n,3)  

		!Change in Momentum in halo cells
		CV_momentum_out(modulo((i-2),nbins(1))+2, & 
			      	modulo((j-2),nbins(2))+2, & 
			      	modulo((k-2),nbins(3))+2,:) = & 
				CV_momentum_out(modulo((i-2),nbins(1))+2,& 
						modulo((j-2),nbins(2))+2,&
						modulo((k-2),nbins(3))+2,:) & 
							+ CV_momentum_out(i,j,k,:)
	enddo

	if (io_type .eq. 'snap') then
		!CV_momentum_out = CV_momentum_out / (tplot*Nvflux_ave)
		m = iter/(Nvflux_ave) + 1 !Initial snapshot taken
	else
		!CV_momentum_out = CV_momentum_out / (tplot*Nvel_ave)
		m = iter/(tplot*Nvel_ave)
	endif
	!Write velocity to file
	buf = CV_momentum_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	inquire(iolength=length) buf
	open (unit=6, file=filename,form='unformatted',access='direct',recl=length)
	write(6,rec=m) buf
	close(6,status='keep')

end subroutine velocity_bin_io


!---------------------------------------------------------------------------------
! Record temperature in a slice through the domain

subroutine temperature_slice_io(ixyz)
use module_parallel_io
use calculated_properties_MD
implicit none

	integer		:: ixyz, m, length

	!Write mass
	call mass_slice_io(ixyz)

	!Write temperature to file
	m = iter/(tplot*NTemp_ave)
	inquire(iolength=length) slice_temperature(1:nbins(ixyz))
	open (unit=6, file=trim(prefix_dir)//'results/Tslice',form='unformatted',access='direct',recl=length)
	write(6,rec=m) slice_temperature(1:nbins(ixyz))
	close(6,status='keep')

end subroutine temperature_slice_io

!---------------------------------------------------------------------------------
! Record temperature in 3D bins throughout domain

subroutine temperature_bin_io(CV_mass_out,CV_temperature_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer					:: n,m,i,j,k
	integer					:: length
	integer					:: CV_mass_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	double precision		:: CV_temperature_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	double precision		:: buf(nbins(1),nbins(2),nbins(3))
	character(4)			:: io_type
	character(13)			:: filename

	!Write mass bins
	call mass_bin_io(CV_mass_out,io_type)

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) 'results/T', io_type

	!---------------Correct for surface fluxes on halo cells---------------
	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nhalobins
		i = halobins(n,1); j = halobins(n,2); k = halobins(n,3)  

		!Change in temperature in halo cells
		CV_temperature_out(modulo((i-2),nbins(1))+2, & 
			      	modulo((j-2),nbins(2))+2, & 
			      	modulo((k-2),nbins(3))+2) = & 
				CV_temperature_out(modulo((i-2),nbins(1))+2,& 
						modulo((j-2),nbins(2))+2,&
						modulo((k-2),nbins(3))+2) & 
							+ CV_temperature_out(i,j,k)
	enddo

	if (io_type .eq. 'snap') then
		!CV_temperature_out = CV_temperature_out / (tplot*Nvflux_ave)
		stop "Error in temp_io - Temperature SNAP called and NTflux_ave not defined"
		!m = iter/(NTflux_ave) + 1 !Initial snapshot taken
	else
		!CV_temperature_out = CV_temperature_out / (tplot*Nvel_ave)
		m = iter/(tplot*NTemp_ave)
	endif
	!Write temperature to file
	buf = CV_temperature_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1)
	inquire(iolength=length) buf
	open (unit=6, file=filename,form='unformatted',access='direct',recl=length)
	write(6,rec=m) buf
	close(6,status='keep')

end subroutine temperature_bin_io

!---------------------------------------------------------------------------------
! Record energy in 3D bins throughout domain

subroutine energy_bin_io(CV_energy_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer					:: n,m,i,j,k
	integer					:: length
	double precision		:: CV_energy_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	double precision		:: buf(nbins(1),nbins(2),nbins(3))
	character(4)			:: io_type
	character(13)			:: filename

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) 'results/e', io_type

	!---------------Correct for surface fluxes on halo cells---------------
	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nhalobins
		i = halobins(n,1); j = halobins(n,2); k = halobins(n,3)  

		!Change in energy in halo cells
		CV_energy_out(modulo((i-2),nbins(1))+2, & 
			      	modulo((j-2),nbins(2))+2, & 
			      	modulo((k-2),nbins(3))+2) = & 
				CV_energy_out(modulo((i-2),nbins(1))+2,& 
						modulo((j-2),nbins(2))+2,&
						modulo((k-2),nbins(3))+2) & 
							+ CV_energy_out(i,j,k)
	enddo

	if (io_type .eq. 'snap') then
		!CV_energy_out = CV_energy_out / (tplot*Nvflux_ave)
		m = iter/(Neflux_ave) + 1 !Initial snapshot taken
	else
		!CV_energy_out = CV_energy_out / (tplot*Nvel_ave)
		m = iter/(tplot*Neflux_ave)
	endif
	!Write velocity to file
	buf = CV_energy_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1)
	inquire(iolength=length) buf
	open (unit=10, file=filename,form='unformatted',access='direct',recl=length)
	write(10,rec=m) buf
	close(10,status='keep')

end subroutine energy_bin_io

!---------------------------------------------------------------------------------
!Calculate Virial Stress in volume

subroutine virial_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer		:: m, length

	!Write virial pressure to file
	m = iter/(tplot*Nstress_ave)
	inquire(iolength=length) Pxy
	open (unit=7, file=trim(prefix_dir)//'results/pvirial',form='unformatted',access='direct',recl=length)
	write(7,rec=m) Pxy
	close(7,status='keep')

end subroutine virial_stress_io

!---------------------------------------------------------------------------------

subroutine VA_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: ixyz, jxyz, m, length
	double precision	:: binvolume
	double precision	:: buf(nbins(1),nbins(2),nbins(3),3,3)

	!Add kinetic and configurational to Pxybin total
	Pxybin(:,:,:,:,:) = 	vvbin(:,:,:,:,:) 		& 
					      + rfbin(  2:nbins(1)+1, 	& 
						      		2:nbins(2)+1, 	& 
						      		2:nbins(3)+1,:,:)/2.d0

	!Average over samples
	Pxybin = Pxybin / Nstress_ave
	vvbin  = vvbin  / Nstress_ave 
	rfbin  = rfbin  / Nstress_ave 

	!Calculate Virial from sum of Volume Averaged Stresses
	do ixyz = 1,3
	do jxyz = 1,3
		Pxy(ixyz,jxyz) = sum(Pxybin(:,:,:,ixyz,jxyz))
	enddo
	enddo
	Pxy = Pxy / volume

	!Write out Virial stress
	call virial_stress_io

	!VA pressure per bin
	binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
	Pxybin = Pxybin / binvolume
	vvbin  = vvbin  / binvolume
	rfbin  = rfbin  / (2.d0*binvolume)

	!Write VA pressure to file
	m = iter/(tplot*Nstress_ave)
	select case (split_kin_config)
	case(0)
		!Write sum of kinetic and configurational
		inquire(iolength=length) Pxybin
		open (unit=7, file=trim(prefix_dir)//'results/pVA',form='unformatted',access='direct',recl=length)
		write(7,rec=m) Pxybin
		close(7,status='keep')
	case(1)
		!Kinetic
		inquire(iolength=length) vvbin
		open (unit=7, file=trim(prefix_dir)//'results/pVA_k',form='unformatted',access='direct',recl=length)
		write(7,rec=m) vvbin
		close(7,status='keep')
		!Configurational
		buf = rfbin(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
		inquire(iolength=length) buf
		open (unit=7, file=trim(prefix_dir)//'results/pVA_c',form='unformatted',access='direct',recl=length)
		write(7,rec=m) 	buf
		close(7,status='keep')
	case default
		stop 'Error in VA/virial extra flag to split_kinetic_& configuartional parts'
	end select

end subroutine VA_stress_io

!===================================================================================
!Integrate virial pressure to get autocorrelations (Green Kubo) viscosity

subroutine viscosity_io
	use module_parallel_io
	use physical_constants_MD
	use calculated_properties_MD
	use librarymod, only : intergrate_trap
	implicit none

	integer				:: m, length
	double precision	:: viscosity

	print*, Pxycorrel

	call intergrate_trap(Pxycorrel,tplot*delta_t,Nstress_ave,viscosity)

	viscosity = (viscosity*volume)/(3.0*Nstress_ave*Nvisc_ave*inputtemperature)

	!Write viscosity to file
	m = iter/(tplot*Nstress_ave*Nvisc_ave)
	inquire(iolength=length) viscosity
	open (unit=7, file=trim(prefix_dir)//'results/visc',form='unformatted',access='direct',recl=length)
	write(7,rec=m) viscosity
	close(7,status='keep')

	Pxycorrel = 0.d0	!Reset Pxycorrel to zero

end subroutine viscosity_io

!=================================================================================
! Record Fluxes accross surfaces of Control Volumes
!=================================================================================
!---------------------------------------------------------------------------------
! Record mass fluxes accross surfaces of Control Volumes

subroutine mass_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: i,j,k,n,m,length
	integer				:: buf(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,1:6)

	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nhalobins
		i = halobins(n,1); j = halobins(n,2); k = halobins(n,3)  
		!Flux over halo cells
		mass_flux(modulo((i-2),nbins(1))+2, & 
				  modulo((j-2),nbins(2))+2, & 
				  modulo((k-2),nbins(3))+2,:) = & 
			mass_flux(modulo((i-2),nbins(1))+2, & 
					  modulo((j-2),nbins(2))+2, &
					  modulo((k-2),nbins(3))+2,:) + mass_flux(i,j,k,:)
	enddo

	!Write six CV surface mass fluxes to file
	m = iter/(Nmflux_ave)
	buf = mass_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	inquire(iolength=length) buf
	open (unit=8, file=trim(prefix_dir)//'results/mflux',form='unformatted',access='direct',recl=length)
	write(8,rec=m) buf
	close(8,status='keep')

end subroutine mass_flux_io

!---------------------------------------------------------------------------------
! Record momentum fluxes accross surfaces of Control Volumes

subroutine momentum_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer					:: ixyz,i,j,k,n,m,length
	double precision		:: buf(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,1:3,1:6)
	double precision		:: binface

	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nhalobins
		i = halobins(n,1); j = halobins(n,2); k = halobins(n,3)  
		!Flux over halo cells
		momentum_flux(	modulo((i-2),nbins(1))+2, & 
			      		modulo((j-2),nbins(2))+2, & 
			      		modulo((k-2),nbins(3))+2,:,:) = & 
				momentum_flux(	modulo((i-2),nbins(1))+2,& 
								modulo((j-2),nbins(2))+2,&
								modulo((k-2),nbins(3))+2,:,:) & 
									+ momentum_flux(i,j,k,:,:)
	enddo

	do ixyz = 1,3
		binface	      = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
			     		(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
		momentum_flux(:,:,:,:,ixyz  )=momentum_flux(:,:,:,:,ixyz  )/(binface) !Bottom
		momentum_flux(:,:,:,:,ixyz+3)=momentum_flux(:,:,:,:,ixyz+3)/(binface) !Top
	enddo

	!Divide momentum flux by averaing period tau=delta_t*Nvflux_ave
	momentum_flux = momentum_flux/(delta_t*Nvflux_ave)

	!Write momentum flux to file
	m = iter/(Nvflux_ave)
	buf = momentum_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	inquire(iolength=length) buf
	open (unit=9, file=trim(prefix_dir)//'results/vflux',form='unformatted',access='direct',recl=length)
	write(9,rec=m) buf
	close(9,status='keep')

end subroutine momentum_flux_io

!---------------------------------------------------------------------------------
! Record stress accross plane

subroutine MOP_stress_io(ixyz)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer		:: ixyz, m, length

	!Divide by number of samples taken
	Pxy_plane = Pxy_plane/(Nstress_ave)

	!Divide by area of domain and factor of 4 for interactions
	Pxy_plane = Pxy_plane/(4*domain(1)*domain(3))

	!Write plane pressures to file
	m = iter/(tplot*Nvflux_ave)
	inquire(iolength=length) Pxy_plane
	open (unit=9, file=trim(prefix_dir)//'results/pplane',form='unformatted',access='direct',recl=length)
	write(9,rec=m) Pxy_plane
	close(9,status='keep')

end subroutine MOP_stress_io

!---------------------------------------------------------------------------------
! Record stress accross surfaces of Control Volumes

subroutine surface_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer							:: ixyz,i,j,k,n,m,length
	double precision				:: buf(nbins(1),nbins(2),nbins(3),1:3,1:6)
	double precision,dimension(3)	:: binface

	!Include halo surface stresses to get correct values for all cells
	do n = 1, nhalobins
		i = halobins(n,1); j = halobins(n,2); k = halobins(n,3)  
		!Set Stresses to value of halo cells
		Pxyface(	modulo((i-2),nbins(1))+2, & 
			      	modulo((j-2),nbins(2))+2, & 
			      	modulo((k-2),nbins(3))+2,:,:) = & 
				Pxyface(	modulo((i-2),nbins(1))+2,& 
						modulo((j-2),nbins(2))+2,&
						modulo((k-2),nbins(3))+2,:,:) & 
							     + Pxyface(i,j,k,:,:)

	enddo

	do ixyz = 1,3
		binface(ixyz) = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
			     		(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
		Pxyface(:,:,:,:,ixyz  ) = 0.25d0 * Pxyface(:,:,:,:,ixyz  )/binface(ixyz) !Bottom
		Pxyface(:,:,:,:,ixyz+3) = 0.25d0 * Pxyface(:,:,:,:,ixyz+3)/binface(ixyz) !Top
	enddo

	!Integration of stress using trapizium rule requires multiplication by timestep
	!so delta_t cancels upon division by tau=elta_t*Nvflux_ave resulting in division by Nvflux_ave
	Pxyface = Pxyface/Nvflux_ave

	!Write surface pressures to file
	m = iter/(Nvflux_ave)
	buf = Pxyface(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	inquire(iolength=length) buf
	open (unit=9, file=trim(prefix_dir)//'results/psurface',form='unformatted',access='direct',recl=length)
	write(9,rec=m) buf
	close(9,status='keep')

end subroutine surface_stress_io

!---------------------------------------------------------------------------------
! Record  energy accross plane

subroutine MOP_energy_io(ixyz)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer		:: ixyz, m, length

	!Divide by number of samples taken
	Pxyv_plane = Pxyv_plane/(Nstress_ave)

	!Divide by area of domain and factor of 4 for interactions
	Pxyv_plane = Pxyv_plane/(4*domain(1)*domain(3))

	!Write plane pressures to file
	m = iter/(tplot*Nvflux_ave)
	inquire(iolength=length) Pxy_plane
	open (unit=10, file=trim(prefix_dir)//'results/eplane',form='unformatted',access='direct',recl=length)
	write(10,rec=m) Pxy_plane
	close(10,status='keep')

end subroutine MOP_energy_io

!---------------------------------------------------------------------------------
! Record energy fluxes accross surfaces of Control Volumes

subroutine energy_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer					:: ixyz,i,j,k,n,m,length
	double precision		:: binface
	double precision		:: buf(nbins(1),nbins(2),nbins(3),6)

	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nhalobins
		i = halobins(n,1); j = halobins(n,2); k = halobins(n,3)  
		!Flux over halo cells
		energy_flux(	modulo((i-2),nbins(1))+2, & 
			      	modulo((j-2),nbins(2))+2, & 
			      	modulo((k-2),nbins(3))+2,:) = & 
				energy_flux(	modulo((i-2),nbins(1))+2,& 
						modulo((j-2),nbins(2))+2,&
						modulo((k-2),nbins(3))+2,:) & 
							+ energy_flux(i,j,k,:)
	enddo

	do ixyz = 1,3
		binface	      = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
			     		(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
		energy_flux(:,:,:,ixyz  )=energy_flux(:,:,:,ixyz  )/(binface) !Bottom
		energy_flux(:,:,:,ixyz+3)=energy_flux(:,:,:,ixyz+3)/(binface) !Top
	enddo

	energy_flux = energy_flux/(delta_t*Neflux_ave)
	!Write energy flux to file
	m = iter/(Neflux_ave)
	buf = energy_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	inquire(iolength=length) buf 
	open (unit=10, file=trim(prefix_dir)//'results/eflux',form='unformatted',access='direct',recl=length)
	write(10,rec=m) buf 
	close(10,status='keep')

end subroutine energy_flux_io

!---------------------------------------------------------------------------------
! Record stress times velocity (power) accross surfaces of Control Volumes

subroutine surface_power_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer							:: ixyz,i,j,k,n,m,length
	double precision				:: buf(nbins(1),nbins(2),nbins(3),6)
	double precision,dimension(3)	:: binface


	!Include halo surface stresses to get correct values for all cells
	do n = 1, nhalobins
		i = halobins(n,1); j = halobins(n,2); k = halobins(n,3)  
		!Set Stresses to value of halo cells
		Pxyvface(	modulo((i-2),nbins(1))+2, & 
			      	modulo((j-2),nbins(2))+2, & 
			      	modulo((k-2),nbins(3))+2,:) = & 
				Pxyvface(	modulo((i-2),nbins(1))+2,& 
						modulo((j-2),nbins(2))+2,&
						modulo((k-2),nbins(3))+2,:) & 
							     + Pxyvface(i,j,k,:)

	enddo

	do ixyz = 1,3
		binface(ixyz) = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
			     		(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
		Pxyvface(:,:,:,ixyz  ) = 0.25d0 * Pxyvface(:,:,:,ixyz  )/binface(ixyz) !Bottom
		Pxyvface(:,:,:,ixyz+3) = 0.25d0 * Pxyvface(:,:,:,ixyz+3)/binface(ixyz) !Top
	enddo

	!Integration of stress using trapizium rule requires multiplication by timestep
	Pxyvface = Pxyvface/Neflux_ave

	!Write surface pressures * velocity to file
	m = iter/(Neflux_ave)
	buf = Pxyvface(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	inquire(iolength=length) buf
	open (unit=10, file=trim(prefix_dir)//'results/esurface',form='unformatted',access='direct',recl=length)
	write(10,rec=m) buf
	close(10,status='keep')

end subroutine surface_power_io

!-----------------------------------------------------------------------------
! Write macroscopic properties to file
!-----------------------------------------------------------------------------
subroutine macroscopic_properties_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none
	
	open(unit=10,file=trim(prefix_dir)//'results/macroscopic_properties',status='replace')
	
	if (potential_flag.eq.0) then
		write(10,'(2a)') &
		'Iteration; 	   VSum;        V^2Sum;        Temp;', &
		'         KE;        PE;         TE;        Pressure;'
		!Print initial conditions for simulations at iteration 0
		write(10,'(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.4)') &
		initialstep,';',vsum,';', v2sum,';', temperature,';', &
		kinenergy,';',potenergy,';',totenergy,';',pressure
	else if (potential_flag.eq.1) then
		write(10,'(2a)') &
		'Iteration; 	   VSum;        V^2Sum;        Temp;', &
		'       KE;     PE (LJ);  PE (FENE); PE (Tot);    TE;       Pressure;'
		!Print initial conditions for simulations at iteration 0
		write(10, '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.4)') &
		initialstep,';',vsum,';', v2sum,';', temperature,';', &
		kinenergy,';',potenergy_LJ,';',potenergy_FENE,';',potenergy,';',totenergy,';',pressure
	end if

end subroutine macroscopic_properties_header

subroutine macroscopic_properties_record
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	if (potential_flag.eq.0) then	
		!write(10,'(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.4)') &
		!iter,';',vsum,';', v2sum,';', temperature,';', &
		!kinenergy,';',potenergy,';',totenergy,';',pressure
		write(10,'(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f10.4,a,f10.4)'), &
		iter,';',vsum,';', v2sum,';', temperature,';', &
		totenergy,';',((density/(globalnp*nd))*virial/2)**2,';',pressure**2,';',(density/(globalnp*nd))*virial/2,';',pressure
	else if (potential_flag.eq.1) then
		write(10,'(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.4)') &
		iter,';',vsum,';', v2sum,';', temperature,';', &
		kinenergy,';',potenergy_LJ,';',potenergy_FENE,';',potenergy,';',totenergy,';',pressure
	end if

end subroutine macroscopic_properties_record

!-----------------------------------------------------------------------------
! Write end-to-end vector time correlation function
!-----------------------------------------------------------------------------
subroutine etev_io
	use module_parallel_io
	implicit none
		
	integer :: m
	integer :: length

	m = (iter-etevtcf_iter0)/tplot + 1
	inquire(iolength=length) etev
	
	if (iter.eq.etevtcf_iter0) then
		open  (15,file=trim(prefix_dir)//'results/etev_xyz',status='replace',form='unformatted',access='direct',recl=length)
		write (15,rec=m) etev 
	else if (iter.gt.etevtcf_iter0) then
		open  (15,file=trim(prefix_dir)//'results/etev_xyz',form='unformatted',access='direct',recl=length)
		write (15,rec=m) etev
	end if
	
	close(15,status='keep')

end subroutine etev_io

subroutine etevtcf_io
use module_parallel_io
implicit none
	
	integer :: m
	integer :: length

	m = (iter-etevtcf_iter0)/tplot + 1
	inquire(iolength=length) etevtcf
	
	if (iter.eq.etevtcf_iter0) then
		open(14,file=trim(prefix_dir)//'results/etevtcf',status='replace',form='unformatted',access='direct',recl=length)
		write(14,rec=m) etevtcf
	else if (iter.gt.etevtcf_iter0) then
		open(14,file=trim(prefix_dir)//'results/etevtcf',form='unformatted',access='direct',recl=length)
		write(14,rec=m) etevtcf
	end if
	
	close(14,status='keep')

end subroutine etevtcf_io

!Write radius of gyration to output file
subroutine r_gyration_io
use module_parallel_io
implicit none

	if (iter.eq.r_gyration_iter0) then
		open(15,file=trim(prefix_dir)//'results/r_gyration',status='replace')
		write(15,'(i8,f15.8)') iter, R_g
	else if (iter.gt.r_gyration_iter0) then
		open(15,file=trim(prefix_dir)//'results/r_gyration',position='append')
		write(15,'(i8,f15.8)') iter, R_g
	end if
	
	close(15,status='keep')

end subroutine r_gyration_io
