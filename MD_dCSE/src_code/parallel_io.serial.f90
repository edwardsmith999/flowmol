!======================================================================
!	        	Serial Read Write Subroutines               

! Serial emulation of parallel i/o

! --- INPUT ROUTINES ---
! setup_command_arguments	Determine restart file (if any) and input file
! setup_inputs_locate		Read input file
! setup_restart_inputs		Read restart files and input file
! setup_restart_microstate	Read Initial configuration from restart file

! --- OUTPUT ROUTINES ---
! simulation_header				Write all variables in ; seperated variable form "Description; name; variable"
! update_simulation_progress    Record current iteration
! parallel_io_final_state 		Write final configuration and simulation details for restart
! parallel_io_vmd				Output for VMD
! parallel_io_vmd_sl			Output for VMD with seperate solid/lquid regions
! parallel_io_vmd_halo			Output Halos

! CV AVERAGING
! mass_slice_io			Write out mass bins for a single dimension slice through domain
! mass_bin_io			Write out mass bins in all 3 dimensions of domain
! velocity_slice_io		Write out velocity bins for a single dimension slice through domain
! velocity_bin_io		Write out velocity bins in all 3 dimensions of domain
!
! FLUX AVERAGING
! virial_stress_io		Write out virial stress
! VA_stress_io			Write out Volume Averaged stress throughout domain
! viscosity_io			Write out viscosity
! mass_flux_io			Write out flux of mass through bin surfaces
! momentum_flux_io		Write out flux of momnetum through bin surfaces
! MOP_stress_io			Write out stress on single plane through domain
! surface_stress_io		Write out stress on all surfaces of bins in domain
! 
!======================================================================

module module_parallel_io
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
	
	integer	:: i,argcount
	logical :: restart_file_exists, input_file_exists
	character(len=200) :: arg,nextarg

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

	logical			:: from_input
	integer 		:: i, n 
	integer,dimension(8)	:: tvalue

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

	call locate(1,'POTENTIAL_FLAG',.true.)	!LJ or FENE potential
	read(1,*) potential_flag
	if (potential_flag.eq.1) then
		call locate(1,'FENE_INFO',.true.)
		read(1,*) chain_length
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
	call locate(1,'INITISE_STEPS',.false.,from_input)
	if (from_input) then
		read(1,*) initise_steps 	!Number of initialisation steps for simulation
	else
		initise_steps = 0
	endif
	call locate(1,'DELTA_RNEIGHBR',.true.) 
	read(1,*) delta_rneighbr 	!Extra distance used for neighbour cell
	call locate(1,'SEED',.false.,from_input)
	if (from_input) then
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

	call locate(1,'DEFINE_SHEAR',.false.,from_input)
	if (from_input) then
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
	
	call locate(1,'WALLSLIDEV',.false.,from_input)
	if (from_input) then
		read(1,*) wallslidev(1)
		read(1,*) wallslidev(2)
		read(1,*) wallslidev(3)
	endif
	call locate(1,'FIXDISTBOTTOM',.false.,from_input)
	if (from_input) then
		read(1,*) fixdistbottom(1)
		read(1,*) fixdistbottom(2)
		read(1,*) fixdistbottom(3)
	endif
	call locate(1,'FIXDISTTOP',.false.,from_input)
	if (from_input) then
		read(1,*) fixdisttop(1)
		read(1,*) fixdisttop(2)
		read(1,*) fixdisttop(3)
	endif
	call locate(1,'SLIDEDISTBOTTOM',.false.,from_input)
	if (from_input) then
		read(1,*) slidedistbottom(1)
		read(1,*) slidedistbottom(2)
		read(1,*) slidedistbottom(3)
	endif
	call locate(1,'SLIDEDISTTOP',.false.,from_input)
	if (from_input) then
		read(1,*) slidedisttop(1)
		read(1,*) slidedisttop(2)
		read(1,*) slidedisttop(3)
	endif
	call locate(1,'TETHEREDDISTBOTTOM',.false.,from_input)
	if (from_input) then
		read(1,*) tethereddistbottom(1)
		read(1,*) tethereddistbottom(2)
		read(1,*) tethereddistbottom(3)
	endif
	call locate(1,'TETHEREDDISTTOP',.false.,from_input)
	if (from_input) then
		read(1,*) tethereddisttop(1)
		read(1,*) tethereddisttop(2)
		read(1,*) tethereddisttop(3)
	endif

	call locate(1,'THERMSTATBOTTOM',.false.,from_input)
	if (from_input) then
		read(1,*) thermstatbottom(1)
		read(1,*) thermstatbottom(2)
		read(1,*) thermstatbottom(3)
	endif
	call locate(1,'THERMSTATTOP',.false.,from_input)
	if (from_input) then
		read(1,*) thermstattop(1)
		read(1,*) thermstattop(2)
		read(1,*) thermstattop(3)
	endif

	call locate(1,'THERMSTAT_FLAG',.false.,from_input)
	if (from_input) then
		read(1,*) thermstat_flag
		select case(thermstat_flag)
		case(0)
			if (abs(maxval(thermstattop   )).ne.0.0 & 
		        .or.abs(maxval(thermstatbottom)).ne.0.0) stop & 
			 "THERMSTATTOP or THERMSTATBOTTOM non zero but THERMSTAT_FLAG_INFO set to off (THERMSTAT_FLAG=0)"
			thermstatbottom = 0.d0; thermstattop = 0.d0 
		case(1)	!N-H thermostat all molecules
			thermstattop 	= initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
			thermstatbottom = initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
		case(2) !N-H PUT all molecules
			thermstattop 	= initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
			thermstatbottom = initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
		case(3)
			if (abs(maxval(thermstattop   )).eq.0.0 & 
		       .and.abs(maxval(thermstatbottom)).eq.0.0) stop & 
			"THERMSTATTOP or THERMSTATBOTTOM must also be specified"
		end select
	endif

	!Flag to determine which outputs are switched on
	call locate(1,'VMD_OUTFLAG',.false.,from_input)
	if (from_input) read(1,*) vmd_outflag
	call locate(1,'MACRO_OUTFLAG',.false.,from_input)
	if (from_input) read(1,*) macro_outflag
	call locate(1,'MASS_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,*) mass_outflag
		if (mass_outflag .ne. 0) 	read(1,*) Nmass_ave
	endif
	call locate(1,'VELOCITY_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,* ) velocity_outflag
		if (velocity_outflag .ne. 0)	read(1,* ) Nvel_ave
	endif
	call locate(1,'PRESSURE_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,* ) pressure_outflag
		if (pressure_outflag .ne. 0)	read(1,* ) Nstress_ave
	endif
	call locate(1,'VISCOSITY_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,* ) viscosity_outflag
		if ( viscosity_outflag .ne. 0)	read(1,* ) Nvisc_ave
	endif
	call locate(1,'MFLUX_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,* ) mflux_outflag
		if (mflux_outflag .ne. 0)	read(1,* ) Nmflux_ave
	endif
	call locate(1,'VFLUX_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,* ) vflux_outflag
		if (vflux_outflag .ne. 0)	read(1,* ) Nvflux_ave
	endif

	call locate(1,'ETEVTCF_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,*) etevtcf_outflag
		if (etevtcf_outflag.ne.0) then
			read(1,*) etevtcf_iter0
	
			if (mod(etevtcf_iter0,tplot).ne.0) then
				etevtcf_iter0 = etevtcf_iter0 + (tplot - mod(etevtcf_iter0,tplot))
				print*, 'Etevtcf must be a multiple of tplot, resetting etevtcf to ', etevtcf_iter0
			end if
		end if
	endif

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

	logical				:: from_input
	integer				:: n, k
	integer 			:: extrasteps
	integer 			:: checkint
	double precision 		:: checkdp

         integer(kind=selected_int_kind(18)) header_pos, end_pos ! 8 byte integer for header address

	!Allocate random number seed
!	call random_seed
	call random_seed(size=n)
	allocate(seed(n))


	!Open file to read integers
	open(2,file=initial_microstate_file,form="unformatted", access="stream",position="append")
        inquire(2,POS=end_pos) ! go the end of file
 
        read(2,pos=end_pos-8) header_pos ! header start is in the last 8 bytes

        header_pos = header_pos +1 ! for compatibility with MPI IO 
                                   ! we store header_pos - 1 in final state 

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
	read(2) chain_length

	read(2) density			!Density of system
	read(2) rcutoff			!Cut off distance for particle interaction
	rcutoff2= rcutoff**2         				!Useful definition to save computational time
	read(2) inputtemperature	!Define initial temperature
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
	read(1,* ) checkdp          !Density of system
	if (checkdp .ne. density) print*, 'Discrepancy between system density', &
	'in input & restart file - restart file will be used'

	call locate(1,'RCUTOFF',.true.)
	read(1,* ) checkdp          !Cut off distance for particle interaction
	if (checkdp .ne. rcutoff) print*, 'Discrepancy between cut off radius', &
	'in input & restart file - restart file will be used'

	call locate(1,'INPUTTEMPERATURE',.true.)
	read(1,* ) checkdp	     !Define initial temperature
	if (checkdp .ne. inputtemperature) then
		inputtemperature = checkdp	
		print*, 'Discrepancy between initial temperature', &
				'in input & restart file - input file will be used'
	end if

	call locate(1,'INITIALNUNITS',.true.)
	read(1,* ) checkint	     	!x dimension split into number of cells
	if (checkint .ne. initialnunits(1)) print*, 'Discrepancy between x domain size', &
				'in input & restart file - restart file will be used'
	read(1,* ) checkint         !y dimension box split into number of cells
	if (checkint .ne. initialnunits(2)) print*, 'Discrepancy between y domain size', &
				'in input & restart file - restart file will be used'
	read(1,* ) checkint	     	!z dimension box split into number of cells
	if (checkint .ne. initialnunits(3)) print*, 'Discrepancy between z domain size', &
				'in input & restart file - restart file will be used'

	!Get number of extra steps, timestep and plot frequency from input file	
	call locate(1,'NSTEPS',.true.)
	read(1,* ) extrasteps       !Number of computational steps
	call locate(1,'DELTA_T',.true.)
	read(1,* ) delta_t          !Size of time step
	call locate(1,'TPLOT',.true.)
	read(1,* ) tplot            !Frequency at which to record results
	call locate(1,'DELTA_RNEIGHBR',.true.)
	read(1,* ) delta_rneighbr   !Extra distance used for neighbour cell
	call locate(1,'SEED',.false.,from_input)
	if (from_input) then
		read(1,*) seed(1) 	!Random number seed value 1
		read(1,*) seed(2) 	!Random number seed value 2
	else
		seed(1) = 1		!Fixed default seed for repeatability
		seed(2) = 2		!Fixed default seed for repeatability
	endif

	!Check periodic BC and shear
	call locate(1,'PERIODIC',.true.)
	read(1,*) periodic(1)
	read(1,*) periodic(2)
	read(1,*) periodic(3)

	call locate(1,'DEFINE_SHEAR',.false.,from_input)
	if (from_input) then
		read(1,*) shear_direction
		read(1,*) shear_iter0
		read(1,*) define_shear_as
		if (define_shear_as.eq.0) read(1,*) shear_velocity
		if (define_shear_as.eq.1) read(1,*) shear_rate
		if (define_shear_as.gt.1) then 
                	call error_abort( 'Poorly defined shear in input file')
        	endif
	endif

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

	call locate(1,'WALLSLIDEV',.false.,from_input)
	if (from_input) then
		read(1,*) wallslidev(1)
		read(1,*) wallslidev(2)
		read(1,*) wallslidev(3)
	endif
	call locate(1,'FIXDISTBOTTOM',.false.,from_input)
	if (from_input) then
		read(1,*) fixdistbottom(1)
		read(1,*) fixdistbottom(2)
		read(1,*) fixdistbottom(3)
	endif
	call locate(1,'FIXDISTTOP',.false.,from_input)
	if (from_input) then
		read(1,*) fixdisttop(1)
		read(1,*) fixdisttop(2)
		read(1,*) fixdisttop(3)
	endif
	call locate(1,'SLIDEDISTBOTTOM',.false.,from_input)
	if (from_input) then
		read(1,*) slidedistbottom(1)
		read(1,*) slidedistbottom(2)
		read(1,*) slidedistbottom(3)
	endif
	call locate(1,'SLIDEDISTTOP',.false.,from_input)
	if (from_input) then
		read(1,*) slidedisttop(1)
		read(1,*) slidedisttop(2)
		read(1,*) slidedisttop(3)
	endif
	call locate(1,'TETHEREDDISTBOTTOM',.false.,from_input)
	if (from_input) then
		read(1,*) tethereddistbottom(1)
		read(1,*) tethereddistbottom(2)
		read(1,*) tethereddistbottom(3)
	endif
	call locate(1,'TETHEREDDISTTOP',.false.,from_input)
	if (from_input) then
		read(1,*) tethereddisttop(1)
		read(1,*) tethereddisttop(2)
		read(1,*) tethereddisttop(3)
	endif
	call locate(1,'THERMSTATBOTTOM',.false.,from_input)
	if (from_input) then
		read(1,*) thermstatbottom(1)
		read(1,*) thermstatbottom(2)
		read(1,*) thermstatbottom(3)
	endif
	call locate(1,'THERMSTATTOP',.false.,from_input)
	if (from_input) then
		read(1,*) thermstattop(1)
		read(1,*) thermstattop(2)
		read(1,*) thermstattop(3)
	endif
	
	call locate(1,'THERMSTAT_FLAG',.false.,from_input)
	if (from_input) then
		read(1,*) thermstat_flag
		select case(thermstat_flag)
		case(0)
			if (abs(maxval(thermstattop   )).ne.0.0 & 
		        .or.abs(maxval(thermstatbottom)).ne.0.0) stop & 
			"THERMSTATTOP or THERMSTATBOTTOM non zero but THERMSTAT_FLAG_INFO set to off (THERMSTAT_FLAG=0)"
			thermstatbottom = 0.d0; thermstattop = 0.d0 
		case(1)
			thermstattop 	= initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
			thermstatbottom = initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
		case(2)
			if (abs(maxval(thermstattop   )).eq.0.0 & 
		       .and.abs(maxval(thermstatbottom)).eq.0.0) stop &
			 "THERMSTATTOP or THERMSTATBOTTOM must also be specified"
		end select
	endif


	!Flag to determine which outputs are switched on
	call locate(1,'VMD_OUTFLAG',.false.,from_input)
	if (from_input) read(1,*) vmd_outflag
	call locate(1,'MACRO_OUTFLAG',.false.,from_input)
	if (from_input) read(1,*) macro_outflag
	call locate(1,'MASS_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,*) mass_outflag
		if (mass_outflag .ne. 0) 	read(1,*) Nmass_ave
	endif
	call locate(1,'VELOCITY_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,* ) velocity_outflag
		if (velocity_outflag .ne. 0)	read(1,* ) Nvel_ave
	endif
	call locate(1,'PRESSURE_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,* ) pressure_outflag
		if (pressure_outflag .ne. 0)	read(1,* ) Nstress_ave
	endif
	call locate(1,'VISCOSITY_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,* ) viscosity_outflag
		if ( viscosity_outflag .ne. 0)	read(1,* ) Nvisc_ave
	endif
	call locate(1,'MFLUX_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,* ) mflux_outflag
		if (mflux_outflag .ne. 0)	read(1,* ) Nmflux_ave
	endif
	call locate(1,'VFLUX_OUTFLAG',.false.,from_input)
	if (from_input) then
		read(1,* ) vflux_outflag
		if (vflux_outflag .ne. 0)	read(1,* ) Nvflux_ave
	endif

	elapsedtime = elapsedtime + delta_t*extrasteps !Set elapsed time to end of simualtion
	initialstep = Nsteps         !Set plot count to final plot of last
	iter = initialstep			 !Set iter to initialstep so that initial record is performed correctly at restart
	Nsteps = Nsteps + extrasteps !Establish final iteration step based on previous

	call locate(1,'ETEVTCF_OUTFLAG',.false.,from_input)
	if (from_input) then
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

	integer :: ixyz, n
	integer, dimension(np) :: chainID, subchainID,left,right

	!Open file at first recorded value
	open(2,file=initial_microstate_file, form='unformatted', access='stream',position="rewind")
	!Read positions
	do n=1,globalnp
		!Read position from file
		do ixyz=1,nd
			read(2) r(n,ixyz)
		enddo
                do ixyz=1,nd
			read(2) v(n,ixyz)
		enddo
	enddo
	call setup_tag				!Setup location of fixed molecules
	do n = 1,np
		call read_tag(n)		!Read tag and assign properties
	enddo

	if (potential_flag.eq.1) then
		read(2) chainID(:)				
		read(2) subchainID(:)				
		read(2) left(:)
		read(2) right(:)
		do n=1,np
			polyinfo_mol(n)%chainID = chainID(n)
			polyinfo_mol(n)%subchainID = subchainID(n)
			polyinfo_mol(n)%left = left(n)
			polyinfo_mol(n)%right= right(n)
		end do

	end if
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
	write(3,*)  'Pressure flag ;  pressure_outflag ;', pressure_outflag
	write(3,*)  'viscosity flag ;  viscosity_outflag ;', viscosity_outflag
	write(3,*)  'mass flux flag ;  mflux_outflag ;', mflux_outflag
	write(3,*)  'velocity flux flag ;  vflux_outflag ;', vflux_outflag
	write(3,*)  'mass average steps ;  Nmass_ave ;', Nmass_ave
	write(3,*)  'velocity average steps ;  Nvel_ave ;', Nvel_ave
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
	write(3,*)  'Domain split into Planes for Pressure Averaging ; nplanes  ;',nplanes 
	write(3,*)  'Separated by distance ;  planespacing  ;', planespacing 
	write(3,*)  'with first plane at ;  planes ;', planes(1)
	write(3,*)	'Shear direction ; shear_direction;', shear_direction

	close(3,status='keep')

end subroutine simulation_header

!------------------------------------------------------------------------------
!Serial version of parallel code to print final_state for restart

subroutine parallel_io_final_state
	use module_parallel_io
	implicit none

	integer :: ixyz,n
	integer, dimension(np) :: chainID, subchainID,right,left
	integer :: int_filesize,dp_filesize
        integer(kind=selected_int_kind(18)) header_pos ! 8 byte integer for header address

	!Rebuild simulation before recording final state
	call sendmols			   			!Exchange particles between processors
	call linklist_deallocateall	   		!Deallocate all linklist components
	call assign_to_cell	  	   			!Re-build linklist every timestep
	call messenger_updateborders(1)	   	!Update borders between processors
	call assign_to_neighbourlist	   	!Setup neighbourlist

	!Remove previous final state file
	open(2,file='results/final_state')
	close(2,status='delete')
	!open(3,file='results/finalvelocities',status='replace')

	
	!Written in this form so each molecule's information is together to allow 
	!re-allocation to seperate processors
	open(2,file='results/final_state', form='unformatted',access='stream')
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
!	close(2,status='keep')

	if (potential_flag.eq.1) then
		do n=1,np
			chainID(n)    = polyinfo_mol(n)%chainID
			subchainID(n) = polyinfo_mol(n)%subchainID
			left(n)       = polyinfo_mol(n)%left
			right(n)	  = polyinfo_mol(n)%right
		end do
!		open(2,file='results/final_state', form='unformatted',access='direct',recl=np)
		write(2) chainID(:)
		write(2) subchainID(:)
		write(2) left(:)
		write(2) right(:)

	end if
 
        close(2,status='keep')

	!Write integer data at end of file	
	open(2,file='results/final_state', form='unformatted',access='stream',position="append")

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
	write(2) periodic(1)	   	 	!Boundary condition flags
	write(2) periodic(2)	   		!Boundary condition flags
	write(2) periodic(3)	   	!Boundary condition flags
	write(2) potential_flag   	!Polymer/LJ potential flag
	write(2) chain_length	   	!Polymer chain length

	write(2) density           !Density of system
	write(2) rcutoff           !Cut off distance for particle interaction
	write(2) inputtemperature  !Define initial temperature
	write(2) delta_t           !Size of time step
	write(2) elapsedtime       !Total elapsed time of all restarted simulations
	write(2) k_c			   	  !FENE spring constant
	write(2) R_0			      !FENE spring max elongation
	write(2) delta_rneighbr	  !Extra distance used for neighbour list cell size
        write(2) header_pos-1     ! -1 for MPI IO compatibility
	close(2,status='keep') !Close final_state file

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
	open(2,file='results/final_state')
	close(2,status='delete')
	!open(3,file='results/finalvelocities',status='replace')

	open(2,file='results/final_state', form='unformatted',access='direct',recl=2)
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
	open(2,file='results/final_state', form='unformatted',access='direct',recl=1)

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
	write(2,rec=12*np+29) chain_length	   !Polymer chain length
	write(2,rec=12*np+30) 0				   !Dummy to make even filesize

	close(2,status='keep') !Close final_state file

end subroutine parallel_io_final_state_old

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd
	use module_parallel_io
	implicit none

	integer				:: i, n, length
	real,dimension(3*np)		:: buf

	buf(1     :  np) = r(:,1)
	buf(np+1  :2*np) = r(:,2)
	buf(2*np+1:3*np) = r(:,3)

	i = iter / tplot

	inquire(iolength=length) buf
	open (unit=4, file="results/vmd_temp.dcd",access='direct',recl=length)
	write(4,rec=i) buf(:)
	close(4,status='keep')

end subroutine parallel_io_vmd

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd_sl
	use module_parallel_io
	implicit none

	integer				:: i, n
	real,dimension(np)		:: Xbuf, Ybuf, Zbuf
	real,dimension(3)		:: rhalfdomain

	Xbuf(:) = r(:,1)
	Ybuf(:) = r(:,2)
	Zbuf(:) = r(:,3)

	rhalfdomain(1) = -halfdomain(1)
	rhalfdomain(2) = -halfdomain(2)
	rhalfdomain(3) = -halfdomain(3)

	i = iter / tplot

	!---Write liquid molecules---
	open (unit=4, file="results/vmd_liquid_temp.dcd",access='direct',recl=1)
	
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
	open (unit=4, file="results/vmd_solid_temp.dcd",access='direct',recl=1)

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

subroutine parallel_io_vmd_halo
	use module_parallel_io
	implicit none

	integer				:: i, n
	real,dimension(halo_np)		:: Xbuf, Ybuf, Zbuf
	real,dimension(3)		:: rhalfdomain

	Xbuf(:) = r(np+1:np+halo_np,1)
	Ybuf(:) = r(np+1:np+halo_np,2)
	Zbuf(:) = r(np+1:np+halo_np,3)

	rhalfdomain(1) = -halfdomain(1)
	rhalfdomain(2) = -halfdomain(2)
	rhalfdomain(3) = -halfdomain(3)

	i = iter / tplot

	open (unit=4, file="results/vmd_halo_temp.dcd",access='direct',recl=1)
	
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

	open (unit=99999, file="results/simulation_progress")
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
	open (unit=5, file="results/mslice",form="unformatted",access='direct',recl=length)
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
	character(4)			:: io_type
	character(13)			:: filename

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) 'results/m', io_type

	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nsurfacebins
		i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
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
	!print*, m,iter,tplot,Nmass_ave, io_type, size(CV_mass_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1))
	!Write mass to file
	inquire(iolength=length) CV_mass_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1)
	open (unit=5,file=filename,form="unformatted",access='direct',recl=length)
	write(5,rec=m) CV_mass_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1)
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
	open (unit=6, file="results/vslice",form="unformatted",access='direct',recl=length)
	write(6,rec=m) slice_momentum(1:nbins(ixyz),:)
	close(6,status='keep')


end subroutine velocity_slice_io

!---------------------------------------------------------------------------------
! Record velocity in 3D bins throughout domain

subroutine velocity_bin_io(CV_mass_out,CV_momentum_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: n,m,i,j,k
	integer				:: length
	integer				:: CV_mass_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	double precision		:: CV_momentum_out(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)
	character(4)			:: io_type
	character(13)			:: filename

	!Write mass bins
	call mass_bin_io(CV_mass_out,io_type)

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) 'results/v', io_type

	!---------------Correct for surface fluxes on halo cells---------------
	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nsurfacebins
		i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  

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
	inquire(iolength=length) CV_momentum_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	open (unit=6, file=filename,form="unformatted",access='direct',recl=length)
	write(6,rec=m) CV_momentum_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	close(6,status='keep')

end subroutine velocity_bin_io

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
	open (unit=7, file="results/pvirial",form="unformatted",access='direct',recl=length)
	write(7,rec=m) Pxy
	close(7,status='keep')

end subroutine virial_stress_io

!---------------------------------------------------------------------------------

subroutine VA_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer			:: ixyz, jxyz, m, length
	double precision	:: binvolume

	!Average over samples
	Pxybin = Pxybin / Nstress_ave 

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

	!Write VA pressure to file
	m = iter/(tplot*Nstress_ave)
	inquire(iolength=length) Pxybin
	open (unit=7, file="results/pVA",form="unformatted",access='direct',recl=length)
	write(7,rec=m) Pxybin
	close(7,status='keep')

end subroutine VA_stress_io

!===================================================================================
!Integrate virial pressure to get autocorrelations (Green Kubo) viscosity

subroutine viscosity_io
use module_parallel_io
use physical_constants_MD
use calculated_properties_MD
implicit none

	integer			:: m, length
	double precision	:: viscosity

	print*, Pxycorrel

	call intergrate_trap(Pxycorrel,tplot*delta_t,Nstress_ave,viscosity)

	viscosity = (viscosity*volume)/(3.0*Nstress_ave*Nvisc_ave*inputtemperature)

	!Write viscosity to file
	m = iter/(tplot*Nstress_ave*Nvisc_ave)
	inquire(iolength=length) viscosity
	open (unit=7, file="results/visc",form="unformatted",access='direct',recl=length)
	write(7,rec=m) viscosity
	close(7,status='keep')

	Pxycorrel = 0.d0	!Reset Pxycorrel to zero

end subroutine viscosity_io

!=================================================================================
! Record Fluxes accross surfaces of Control Volumes


!---------------------------------------------------------------------------------
! Record mass fluxes accross surfaces of Control Volumes

subroutine mass_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer		:: i,j,k,n,m,length

	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nsurfacebins
		i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
		!Flux over halo cells
		mass_flux(modulo((i-2),nbins(1))+2,modulo((j-2),nbins(2))+2,modulo((k-2),nbins(3))+2,:) = & 
		mass_flux(modulo((i-2),nbins(1))+2,modulo((j-2),nbins(2))+2,modulo((k-2),nbins(3))+2,:) + mass_flux(i,j,k,:)
	enddo

	!Write six CV surface mass fluxes to file
	m = iter/(Nmflux_ave)
	inquire(iolength=length) mass_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	open (unit=8, file="results/mflux",form="unformatted",access='direct',recl=length)
	write(8,rec=m) mass_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	close(8,status='keep')

end subroutine mass_flux_io

!---------------------------------------------------------------------------------
! Record momentum fluxes accross surfaces of Control Volumes

subroutine momentum_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: ixyz,i,j,k,n,m,length
	double precision		:: binface

	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nsurfacebins
		i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
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

	!momentum_flux = momentum_flux/(delta_t*Nvflux_ave)
	!Write momnetum flux to file
	m = iter/(Nvflux_ave)
	inquire(iolength=length) momentum_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	open (unit=9, file="results/vflux",form="unformatted",access='direct',recl=length)
	write(9,rec=m) momentum_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
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
	open (unit=9, file="results/pplane",form="unformatted",access='direct',recl=length)
	write(9,rec=m) Pxy_plane
	close(9,status='keep')

end subroutine MOP_stress_io


!---------------------------------------------------------------------------------
! Record stress accross surfaces of Control Volumes

subroutine surface_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: ixyz,i,j,k,n,m,length
	double precision,dimension(3)	:: binface

	!Include halo surface stresses to get correct values for all cells
	do n = 1, nsurfacebins
		i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
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
	Pxyface = delta_t*Pxyface!/Nvflux_ave

	!Write surface pressures to file
	m = iter/(Nvflux_ave)
	inquire(iolength=length) Pxyface(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	open (unit=9, file="results/psurface",form="unformatted",access='direct',recl=length)
	write(9,rec=m) Pxyface(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	close(9,status='keep')

end subroutine surface_stress_io

!-----------------------------------------------------------------------------
! Write macroscopic properties to file
!-----------------------------------------------------------------------------
subroutine macroscopic_properties_header
use module_parallel_io
use calculated_properties_MD
implicit none
	
	open(unit=10,file='results/macroscopic_properties',status='replace')
	
	if (potential_flag.eq.0) then
		write(10,'(2a)'), &
		'Iteration; 	   VSum;        V^2Sum;        Temp;', &
		'         KE;        PE;         TE;        Pressure;'
		!Print initial conditions for simulations at iteration 0
		write(10,'(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.4)'), &
		initialstep,';',vsum,';', v2sum,';', temperature,';', &
		kinenergy,';',potenergy,';',totenergy,';',pressure
	else if (potential_flag.eq.1) then
		write(10,'(2a)'), &
		'Iteration; 	   VSum;        V^2Sum;        Temp;', &
		'       KE;     PE (LJ);  PE (FENE); PE (Tot);    TE;       Pressure;'
		!Print initial conditions for simulations at iteration 0
		write(10, '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.4)'), &
		initialstep,';',vsum,';', v2sum,';', temperature,';', &
		kinenergy,';',potenergy_LJ,';',potenergy_FENE,';',potenergy,';',totenergy,';',pressure
	end if
	

end subroutine macroscopic_properties_header

subroutine macroscopic_properties_record
use module_parallel_io
use calculated_properties_MD
implicit none

	if (potential_flag.eq.0) then	
		write(10,'(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.4)'), &
		iter,';',vsum,';', v2sum,';', temperature,';', &
		kinenergy,';',potenergy,';',totenergy,';',pressure
	else if (potential_flag.eq.1) then
		write(10,'(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.4)'), &
		iter,';',vsum,';', v2sum,';', temperature,';', &
		kinenergy,';',potenergy_LJ,';',potenergy_FENE,';',potenergy,';',totenergy,';',pressure
	end if

end subroutine macroscopic_properties_record

!-----------------------------------------------------------------------------
! Write end-to-end vector time correlation function
!-----------------------------------------------------------------------------
subroutine etevtcf_io
use module_parallel_io
implicit none
	
	integer :: m
	integer :: length

	m = (iter-etevtcf_iter0)/tplot + 1
	inquire(iolength=length) etevtcf
	
	if (iter.eq.etevtcf_iter0) then
		open(14,file='results/etevtcf',status='replace',form='unformatted',access='direct',recl=length)
		write(14,rec=m) etevtcf
	else if (iter.gt.etevtcf_iter0) then
		open(14,file='results/etevtcf',form='unformatted',access='direct',recl=length)
		write(14,rec=m) etevtcf
	end if
	
	close(14,status='keep')

end subroutine etevtcf_io
