!======================================================================
!	        	Parallel Read Write Subroutines               

! Parallel i/o

! --- INPUT ROUTINES ---
! setup_inputs				Read input file
! setup_restart_inputs		Read restart files and input file
! setup_restart_microstate	Read Initial configuration from restart file

! --- OUTPUT ROUTINES ---
! parallel_io_final_state 	Write final configuration and simulation details for restart
! parallel_io_vmd			Output for VMD
! parallel_io_vmd_sl		Output for VMD with seperate solid/lquid regions
! parallel_io_vmd_halo		Output Halos
! CV AVERAGING
! mass_averaging			slice/CV
! cumulative_mass			slice/CV
! mass_snapshot				CV
! momentum_averaging		slice/CV
! cumulative_momentum		slice/CV
! momentum_snapshot			CV
! pressure_averaging		domain/CV
! cumulative_pressure		domain/CV
! FLUX AVERAGING
! mass_flux_averaging		CV_surface
! cumulative_mass_flux		CV_surface
! momentum_flux_averaging	plane(MOP)/CV_surface
! cumulative_momentum_flux	plane(MOP)/CV_surface
! surface_pressure			plane(MOP)/CV_surface
! cumulative_pressure		plane(MOP)/CV_surface
! 
!======================================================================

module module_parallel_io
	use mpi
	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use polymer_info_MD
	use shear_info_MD
	use calculated_properties_MD
	use messenger, only : MD_COMM
	use interfaces

	integer		:: restartfileid, fileid !File name used for parallel i/o

end module

!======================================================================
!	        		INPUTS			              =
!======================================================================


!=============================================================================
! setup_command_arguments
! Checks for command-line arguments passed to the program and assigns the 
! relevant values.
!-----------------------------------------------------------------------------
subroutine setup_command_arguments
use module_parallel_io
implicit none
	
	integer	:: i,argcount
	logical :: restart_file_exists, input_file_exists
	character(len=200) :: arg,nextarg

	if (irank .eq. iroot) then

		!Set default values in case they aren't specified by user
		restart = .false.			
		input_file_exists = .false.
		restart_file_exists = .false.
!		input_file = 'MD.in'
!		initial_microstate_file = 'results/final_state'
		input_file = trim(prefix_dir)//'MD.in'
		initial_microstate_file = trim(prefix_dir)//'final_state'

		argcount = command_argument_count()

		if (argcount.gt.0) then			!If more than 0 arguments	
				
			do i=1,argcount 		!Loop through all arguments...
				
				call get_command_argument(i,arg)	!Reading two at once
                                if ( i < argcount) then 
                                        call get_command_argument(i+1,nextarg)
                                else
                                        nextarg=''
                                endif
				
				if (trim(arg).eq.'-r' .and. nextarg(1:1).ne.'-') then
					initial_microstate_file = trim(nextarg)
					inquire(file=initial_microstate_file, exist=restart) !Check file exists
				end if			

				if (trim(arg).eq.'-i' .and. nextarg(1:1).ne.'-') then
					input_file = trim(nextarg)
				end if

			end do

		end if

		inquire(file=input_file, exist=input_file_exists)	!Check file exists
		if (input_file.eq.'./MD.in'.and..not. input_file_exists) input_file = './default.in'
		inquire(file=input_file, exist=input_file_exists)

		if(.not. input_file_exists) then
			print*, 'Input file ', trim(input_file), ' not found. Stopping simulation.'
			call error_abort
		end if

	endif
	

end subroutine setup_command_arguments

!-----------------------------------------------------------------------------
! Subroutine:	setup_inputs
! Author(s):	David Trevelyan & Edward Smith
! Description:
!		The input file MD.in contains capitalised keywords followed by
!		numerical values. The "locate" subroutine rewinds to the beginning
!		of MD.in and scans each line until the keyword is matched.
!
!		Consequently, the file position is set for the next statement to read
!		the line underneath the previously "located" keyword. 
!-----------------------------------------------------------------------------

subroutine setup_inputs
	use module_parallel_io
	use librarymod, only : locate
	implicit none
	
	integer 			:: k, n, tvalue(8), ios
	logical 			:: found_in_input
	character(20)		:: readin_format

!	call random_seed
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
	call locate(1,'FORCE_LIST',.true.)	!LJ or FENE potential
	read(1,*) force_list
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
		        .or.abs(maxval(thermstatbottom)).ne.0.0) call error_abort(&
                             "THERMSTATTOP or THERMSTATBOTTOM non zero but THERMSTAT_FLAG_INFO set&
                             & to off (THERMSTAT_FLAG=0)")
			thermstatbottom = 0.d0; thermstattop = 0.d0 
		case(1) !N-H thermostat all molecules
			thermstattop 	= initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
			thermstatbottom = initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
		case(2) !N-H PUT all molecules
			thermstattop 	= initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
			thermstatbottom = initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
		case(3)
			if (abs(maxval(thermstattop   )).eq.0.0 & 
		       .and.abs(maxval(thermstatbottom)).eq.0.0) call error_abort(& 
			"THERMSTATTOP or THERMSTATBOTTOM must also be specified")
		end select
        else 
                ! default initialisation
                thermstatbottom = 0.d0; thermstattop = 0.d0
	endif
	!Flag to determine if output is switched on
	call locate(1,'VMD_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) vmd_outflag
		if (vmd_outflag .ne. 0) then
			read(1,*) Nvmd_intervals	!Number of vmd intervals
			if (Nvmd_intervals .gt. 20) then
				print*, "Number of VMD intervals greater than 20 or not specified, setting on for all simualtion"
				Nvmd_intervals = 0
			endif
			if (Nvmd_intervals .eq. 0) then
				allocate(vmd_intervals(2,1))
				vmd_intervals(1,1) = 1; vmd_intervals(2,1) = huge(1)
			else
				allocate(vmd_intervals(2,Nvmd_intervals))
				write(readin_format,'(a,i5,a)') '(',2*Nvmd_intervals,'i)'
				read(1,trim(readin_format)) vmd_intervals
#if USE_COUPLER
				!NEED SOME SORT OF coupler total simulation time retrival here!!
				print*, "WARNING - CHECK VMD INTERVALS is not greater than coupled number of steps"
#else
				if (maxval(vmd_intervals) .gt. Nsteps) &
                                     call error_abort("Specified VMD interval greater than Nsteps")
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
!		call random_seed
!               randomisations 
		call random_seed(get=seed(1:n))
                call date_and_time(values=tvalue)
                seed=IEOR(tvalue(8)+irank,seed)
        else 
                !Assign different random number seed to each processor
	        seed =  irank
	endif

	!Assign seed to random number generator
	call random_seed(put=seed(1:n))

	elapsedtime = 1.d0*delta_t*Nsteps !Set elapsed time to end of simualtion

end subroutine setup_inputs

!-------------------------------------------------------------------------------------
!                    Restart Simulation inputs
! Set up inputs on every processor, based on the final state of a previous simulation
!-------------------------------------------------------------------------------------

subroutine setup_restart_inputs
	use module_parallel_io
	use librarymod, only : locate
	implicit none

	logical							:: found_in_input
	integer							:: n, k
	integer 						:: extrasteps
	integer 						:: checkint
    integer(MPI_OFFSET_KIND)        :: ofs, header_ofs
    integer(selected_int_kind(18))  :: header_pos
	double precision 				:: checkdp
	character(20)					:: readin_format

	!Allocate random number seed
!	call random_seed
	call random_seed(size=n)
	allocate(seed(n))

	!=====================================================================================================!
	!========================   R E A D    R E S T A R T    H E A D E R   ================================!

	!Open on a single process and broadcast
	if (irank .eq. iroot) then
            
	    call MPI_File_open(MPI_COMM_SELF, initial_microstate_file, & 
		MPI_MODE_RDONLY , MPI_INFO_NULL, restartfileid, ierr)

	    ! read the size of offset, just in case we hit a system that still uses 32 bits addresses
	    ! read in a 8 bit integer
	    ofs = -8
	    call mpi_file_seek(restartfileid,ofs,mpi_seek_end,ierr)
	    call mpi_file_read(restartfileid,header_pos	   ,1,mpi_integer8,MPI_STATUS_IGNORE,ierr)

	    header_ofs = header_pos
	    
	    call MPI_FILE_SEEK(restartfileid,header_ofs,MPI_SEEK_SET,ierr)

	    call MPI_File_read(restartfileid,globalnp        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,initialnunits   ,3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,Nsteps          ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,tplot           ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,seed            ,2,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,periodic        ,3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,potential_flag  ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,nmonomers       ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,npx             ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,npy             ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,npz             ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		
	    call MPI_File_read(restartfileid,density         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,rcutoff         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,delta_t         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,elapsedtime     ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,k_c             ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,R_0             ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,delta_rneighbr  ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

		call MPI_File_close(restartfileid,ierr)

		!Check if values from input file are different and alert user - all processors have
		!read the same file so only need to check on one processor
		!But aparently more stuff must be read from input file, this is not very nice at the moment

		open(1,file=input_file)

		call locate(1,'DENSITY',.true.)
		read(1,* ) checkdp                    !Density of system
		if (checkdp .ne. density)             print*, 'Discrepancy between system density', &
                                              'in input & restart file - restart file will be used'
		call locate(1,'RCUTOFF',.true.)
		read(1,* ) checkdp                    !Cut off distance for particle interaction
		if (checkdp .ne. rcutoff)             print*, 'Discrepancy between cut off radius', &
		                                      'in input & restart file - restart file will be used'
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
		  if (checkdp.ne.k_c)                 print*, 'Discrepancy between k_c', &
		                                      'in input & restart file - restart file will be used'
		  read(1,*) checkdp                   !k_c - FENE spring constant
		  if (checkdp.ne.R_0)                 print*, 'Discrepancy between R_0', &
		                                      'in input & restart file - restart file will be used'
		end if	
		
		call locate(1,'PROCESSORS',.true.)
		if (npx .eq. 0 .and. npy .eq. 0 .and. npz .eq. 0) then
			read(1,*) npx
			read(1,*) npy
			read(1,*) npz
			!call setup_restart_microstate_p_to_s  !todo
		else
			read(1,*) checkint
			if (checkint .ne. npx) call error_abort('Number of processors in &
									      			 input file does not match the restart file.')
			read(1,*) checkint
			if (checkint .ne. npy) call error_abort('Number of processors in &
								           			 input file does not match the restart file.')
			read(1,*) checkint
			if (checkint .ne. npz) call error_abort('Number of processors in &
		                                             input file does not match the restart file.')
			!TODO: CHANGE SETUP_RESTART_MICROSTATE TO ONLY READ DATA LOCAL TO PROCESSOR
		end if

		close(1,status='keep')

	endif

	!=============  E N D    R E A D    R E S T A R T    H E A D E R   ==============================!
	!================================================================================================!
        
	! read some more
	open(1,file=input_file)

	!Check periodic BC and shear
	call locate(1,'PERIODIC',.true.)
	read(1,*) periodic(1)
	read(1,*) periodic(2)
	read(1,*) periodic(3)
	
	call locate(1,'NSTEPS',.true.)
	read(1,* ) extrasteps                                     !Number of computational steps
	call locate(1,'DELTA_T',.true.)
	read(1,* ) delta_t                                        !Size of time step
	call locate(1,'TPLOT',.true.)
	read(1,* ) tplot                                          !Frequency at which to record results
	call locate(1,'INITIALISE_STEPS',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) initialise_steps                            !Number of initialisation steps
	else
		initialise_steps = 0
	endif
	call locate(1,'DELTA_RNEIGHBR',.true.)
	read(1,* ) delta_rneighbr                                 !Extra distance used for neighbour cell
	call locate(1,'SEED',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) seed(1)                                     !Random number seed value 1
		read(1,*) seed(2)                                     !Random number seed value 2
	else
		seed(1) = 1                                           !Fixed default seed for repeatability
		seed(2) = 2                                           !Fixed default seed for repeatability
	endif
	
	call locate(1,'INTEGRATION_ALGORITHM',.true.)
	read(1,*) integration_algorithm
	call locate(1,'ENSEMBLE',.true.)
	read(1,*) ensemble
	call locate(1,'FORCE_LIST',.true.)                        !LJ or FENE potential
	read(1,*) force_list

	call locate(1,'INPUTTEMPERATURE',.true.)
	read(1,* ) inputtemperature                               !Define initial temperature
	
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
	!		 		 [  o     o ]
	!a (1 cell size) [     o    ]  a/2 (distance between molcules)	
	!		 		 [  o     o ]
	!		  		 [__________]  a/4 (distance from bottom of domain)
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
			 "THERMSTATTOP or THERMSTATBOTTOM non zero but THERMSTAT_FLAG_INFO&
			 &  set to off (THERMSTAT_FLAG=0)")
			thermstatbottom = 0.d0; thermstattop = 0.d0 
		case(1)
			thermstattop 	= initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
			thermstatbottom = initialnunits(:)/((density/4)**(1.d0/nd))	!Whole domain size
		case(2)
			if (abs(maxval(thermstattop   )).eq.0.0 & 
			 .and.abs(maxval(thermstatbottom)).eq.0.0) & 
				call error_abort("THERMSTATTOP or THERMSTATBOTTOM must also be specified")
		end select
	endif
	
	!Flag to determine if output is switched on
	call locate(1,'VMD_OUTFLAG',.false.,found_in_input)
	if (found_in_input) then
		read(1,*) vmd_outflag
		if (vmd_outflag .ne. 0) then
			read(1,*) Nvmd_intervals	!Number of vmd intervals
			if (Nvmd_intervals .gt. 20) then
				print*, "Number of VMD intervals greater than 20 or not specified, setting on for all simualtion"
				Nvmd_intervals = 0
			endif
			if (Nvmd_intervals .eq. 0) then
				allocate(vmd_intervals(2,1))
				vmd_intervals(1,1) = 1; vmd_intervals(2,1) = huge(1)
			else
				   allocate(vmd_intervals(2,Nvmd_intervals))
				   write(readin_format,'(a,i5,a)') '(',2*Nvmd_intervals,'i)'
				   read(1,trim(readin_format)) vmd_intervals
#if USE_COUPLER
				   !NEED SOME SORT OF coupler total simulation time retrival here!!
				   print*, "WARNING - CHECK VMD INTERVALS is not greater than coupled number of steps"
#else
				   if (maxval(vmd_intervals) .gt. Nsteps) &
				   call error_abort("Specified VMD interval greater than Nsteps")
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
		  if (pressure_outflag .ne. 0)	read(1,* ) Nstress_ave
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
           
           
		!---------------Broadcast data read by root to all other processors-------------------------!
		! temporary np, exact value will be fixed after reading r and v arrays
		! np is needed in set_parameters_outputs that is called before microstate initialisation
		! when restart is true
		call MPI_BCAST(globalnp,          1,MPI_integer,iroot-1,MD_COMM,ierr)
		np = globalnp/nproc
		call MPI_BCAST(rcutoff,           1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
		rcutoff2 = rcutoff**2             !Useful definition to save computational time

		call MPI_BCAST(initialnunits(1),  3,MPI_integer,iroot-1,MD_COMM,ierr)
		call MPI_BCAST(Nsteps,            1,MPI_integer,iroot-1,MD_COMM,ierr)
		call MPI_BCAST(potential_flag,    1,MPI_integer,iroot-1,MD_COMM,ierr)
		call MPI_BCAST(nmonomers,         1,MPI_integer,iroot-1,MD_COMM,ierr)
		call MPI_BCAST(density,           1,MPI_double_precision,iroot-1,MD_COMM,ierr)
		call MPI_BCAST(inputtemperature,  1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
		call MPI_BCAST(elapsedtime,       1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
		call MPI_BCAST(k_c,               1,MPI_double_precision,iroot-1,MD_COMM,ierr)
		call MPI_BCAST(R_0,               1,MPI_double_precision,iroot-1,MD_COMM,ierr)
	   !call MPI_BCAST(delta_rneighbr,1,MPI_DOUBLE_PRECISION,iroot-1,MD_COMM,ierr)
	   !call MPI_BCAST(extrasteps,1,MPI_integer,iroot-1,MD_COMM,ierr)
	   !call MPI_BCAST(nmonomers,1,MPI_integer,iroot-1,MD_COMM,ierr)
	   !call MPI_BCAST(delta_t,1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
	   
		! t plot, seed, periodic, potential flag nmonomers are read from input as well. Pleese fix ! 
	   !call MPI_BCAST(tplot,1,MPI_integer,iroot-1,MD_COMM,ierr)
	   ! seed ?
	   ! periodic ?
	   
	   ! data from bellow is read by all ranks from input file, no need to broadcase
	   !call MPI_BCAST(delta_rneighbr,1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
	   !call MPI_BCAST(vmd_outflag,1,MPI_integer,iroot-1,MD_COMM,ierr)
	   !call MPI_BCAST(macro_outflag,1,MPI_integer,iroot-1,MD_COMM,ierr) 
	   !call MPI_BCAST(velocity_outflag,1,MPI_integer,iroot-1,MD_COMM,ierr) 
	   !call MPI_BCAST(Nvel_ave,1,MPI_integer,iroot-1,MD_COMM,ierr) 
	   !call MPI_BCAST(pressure_outflag,1,MPI_integer,iroot-1,MD_COMM,ierr)
	   !call MPI_BCAST(Nstress_ave,1,MPI_integer,iroot-1,MD_COMM,ierr)  
	   !call MPI_BCAST(Nvisc_ave,1,MPI_integer,iroot-1,MD_COMM,ierr)
	   
	   elapsedtime = elapsedtime + delta_t*extrasteps !Set elapsed time to end of simualtion
	   initialstep = Nsteps         !Set plot count to final plot of last
	   Nsteps = Nsteps + extrasteps !Establish final iteration step based on previous
           
end subroutine setup_restart_inputs

!------------------------------------------------------------------------------

subroutine setup_restart_microstate
	use module_parallel_io
	implicit none
	!include 'mpif.h'

	integer 							:: n, nl, ixyz, procassign
	integer								:: dp_datasize
	integer(kind=MPI_OFFSET_KIND)   	:: disp
	double precision, dimension (2*nd)	:: rvc !Temporary variable
	integer, dimension(:), allocatable  :: monomerc
	
	allocate(monomerc(4+nmonomers))

	!Determine size of datatypes
  	call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

	!Open restart file on all processor
	call MPI_FILE_OPEN(MD_COMM, initial_microstate_file, & 
		MPI_MODE_RDONLY , MPI_INFO_NULL, restartfileid, ierr)

	nl = 0		!Reset local molecules count nl

	!read positions
	do n=1,globalnp

		!---------- For all molecule positions ------------

		!Move through location of position co-ordinates
!		disp =  2 * (n-1) * nd * dp_datasize	

		!Set each processor to that location and write particlewise
!		call MPI_FILE_SET_VIEW(restartfileid, disp, MPI_double_precision, & 
! 					MPI_double_precision, 'native', MPI_INFO_NULL, ierr)
		select case(potential_flag)
		case(0)
			call MPI_FILE_READ_ALL(restartfileid, rvc(:), 2*nd, MPI_double_precision, & 
			                       MPI_STATUS_IGNORE, ierr) !Read position from file
			!Use integer division to determine which processor to assign molecule to
			procassign = ceiling((rvc(1)+globaldomain(1)/2.d0)/domain(1))
			if (procassign .ne. iblock) cycle
			procassign = ceiling((rvc(2)+globaldomain(2)/2.d0)/domain(2))
			if (procassign .ne. jblock) cycle
			procassign = ceiling((rvc(3)+globaldomain(3)/2.d0)/domain(3))
			if (procassign .ne. kblock) cycle

			!If molecules is in the domain then add to processor's total
			nl = nl + 1 !Local molecule count

			!Correct to local coordinates
			r(nl,1) = rvc(1)-domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
			r(nl,2) = rvc(2)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
			r(nl,3) = rvc(3)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)

			!print'(2(i8,a),3(f10.5,a))', iblock, ';', n, ';', r(nl,1), ';', r(nl,2), ';', r(nl,3), ';'

			!Read corresponding velocities
	!		call MPI_FILE_READ(restartfileid, vc(:), nd, MPI_double_precision, & 
	!					MPI_STATUS_IGNORE, ierr) !Read velocity from file
			v(nl,:) = rvc(nd+1:)

		case(1)
			call MPI_FILE_READ_ALL(restartfileid, rvc(:), 2*nd, MPI_double_precision, & 
			                       MPI_STATUS_IGNORE, ierr)
			call MPI_FILE_READ_ALL(restartfileid, monomerc, 4+nmonomers, MPI_integer, &
			                       MPI_STATUS_IGNORE, ierr)
			
			!Use integer division to determine which processor to assign molecule to
			procassign = ceiling((rvc(1)+globaldomain(1)/2.d0)/domain(1))
			if (procassign .ne. iblock) cycle
			procassign = ceiling((rvc(2)+globaldomain(2)/2.d0)/domain(2))
			if (procassign .ne. jblock) cycle
			procassign = ceiling((rvc(3)+globaldomain(3)/2.d0)/domain(3))
			if (procassign .ne. kblock) cycle

			!If molecule is in the domain then add to processor's total
			nl = nl + 1 !Local molecule count

			!Correct to local coordinates
			r(nl,1) = rvc(1)-domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
			r(nl,2) = rvc(2)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
			r(nl,3) = rvc(3)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)
			
			!Assign corresponding velocity
			v(nl,:)     = rvc(nd+1:)
			
			!Assign corresponding monomer info
			monomer(nl)%chainID     = monomerc(1)
			monomer(nl)%subchainID  = monomerc(2)
			monomer(nl)%funcy       = monomerc(3)
			monomer(nl)%glob_no     = monomerc(4) 
			monomer(nl)%bondflag(:) = monomerc(5:)
		
		case default
		end select
		
	enddo

	np = nl	!Correct local number of particles on processor

	!Close file used to load initial state and remove if call "final_state" 
	!to prevent confusion with final state of current run
	call MPI_FILE_CLOSE(restartfileid, ierr)
	if (initial_microstate_file .eq. './results/final_state') then
		call MPI_FILE_DELETE(initial_microstate_file, MPI_INFO_NULL, ierr)
	endif

	call setup_tag				!Setup location of fixed molecules
	do n = 1,np
		call read_tag(n)		!Read tag and assign properties
	enddo
		
	deallocate(monomerc)

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

	close(3,status='keep')

end subroutine simulation_header


!------------------------------------------------------------------------
!Write positions and velocities of molecules to a file for restart

subroutine parallel_io_final_state
	use module_parallel_io
	use polymer_info_MD
	implicit none
	!include 'mpif.h'

	integer				   					:: n, i
	integer 			   					:: dp_datasize,int_datasize
	integer(kind=MPI_OFFSET_KIND)      		:: disp, procdisp, filesize
    integer(kind=selected_int_kind(18))     :: header_pos
	!type(monomer_info)                      :: monomerwrite
	integer, dimension(:), allocatable      :: monomerwrite
	double precision, dimension(nd)			:: Xwrite	!Temporary variable used in write
	double precision, dimension(nd,2*np)	:: buf		!Temporary variable used in write

	!allocate(monomerwrite%bondflag(nmonomers))
	allocate(monomerwrite(4+nmonomers))

	!Rebuild simulation before recording final state
	call linklist_deallocateall	   !Deallocate all linklist components
	call sendmols			   !Exchange particles between processors
	call assign_to_cell	  	   !Re-build linklist every timestep
	call messenger_updateborders(1)	   !Update borders between processors
	call assign_to_neighbourlist	   !Setup neighbourlist

	!Build array of number of particles on neighbouring
	!process' subdomains on current proccess
	call globalGathernp

	!Determine size of datatypes
  	call MPI_type_size(MPI_double_precision,dp_datasize,ierr)
	call MPI_type_size(MPI_integer,int_datasize,ierr)

	!Adjust r according to actual location for storage according
	!to processor topology with r = 0 at centre
	r(:,1) = r(:,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
	r(:,2) = r(:,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
	r(:,3) = r(:,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

	!Obtain displacement of each processor using all other procs' np
	!with 3 position and 3 velocity components for each molecule
	procdisp = 0
	select case (potential_flag)
	case(0)
		do i=1,irank-1
			procdisp = procdisp + 2*nd*procnp(i)*dp_datasize
		enddo
	case(1)
		do i=1,irank-1
			procdisp = procdisp + 2*nd*procnp(i)*dp_datasize + (4+nmonomers)*procnp(i)*int_datasize
		enddo
	case default
	end select
		
	!Remove previous final state file
	call MPI_FILE_DELETE(trim(prefix_dir)//'results/final_state', MPI_INFO_NULL, ierr)

	!Open file on all processors
	call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/final_state', & 
			MPI_MODE_RDWR + MPI_MODE_CREATE, & 
			MPI_INFO_NULL, restartfileid, ierr)

	!-------------Write coordinates--------------------

	!Obtain location to write in file
	disp =  procdisp

	!Set each processor to that location and write particlewise
	call MPI_FILE_SET_VIEW(restartfileid, disp, MPI_double_precision, & 
 		MPI_double_precision, 'native', MPI_INFO_NULL, ierr)

	select case (potential_flag)
	case(0)
		!do n=1,np
		!	Xwrite = r(n,:) !Load into temp in case r dimensions are non contiguous
		!	call MPI_FILE_WRITE(restartfileid, Xwrite, nd, MPI_double_precision, & 
		!				MPI_STATUS_IGNORE, ierr) 
		!	Xwrite = v(n,:) !Load into temp in case v dimensions are non contiguous
		!	call MPI_FILE_WRITE(restartfileid, Xwrite, nd, MPI_double_precision, & 
		!				MPI_STATUS_IGNORE, ierr) 
		!enddo

	 	do n = 1, np
	 		buf(:,2*n-1) = r(n,:)
	 		buf(:,2*n  ) = v(n,:)
	 	enddo
	 	call MPI_FILE_WRITE(restartfileid, buf,2*np*nd, & 
	 	 						MPI_double_precision, MPI_STATUS_IGNORE, ierr)

	case(1)
		do n=1,np
			Xwrite = r(n,:) !Load into temp in case r dimensions are non contiguous
			call MPI_FILE_WRITE(restartfileid, Xwrite, nd, MPI_double_precision, & 
						MPI_STATUS_IGNORE, ierr) 
			Xwrite = v(n,:) !Load into temp in case v dimensions are non contiguous
			call MPI_FILE_WRITE(restartfileid, Xwrite, nd, MPI_double_precision, & 
						MPI_STATUS_IGNORE, ierr) 
			
			monomerwrite(1)  = monomer(n)%chainID
			monomerwrite(2)  = monomer(n)%subchainID
			monomerwrite(3)  = monomer(n)%funcy
			monomerwrite(4)  = monomer(n)%glob_no
			monomerwrite(5:) = monomer(n)%bondflag(:)
			
			call MPI_FILE_WRITE(restartfileid, monomerwrite, 4+nmonomers, MPI_integer, &
			                    MPI_STATUS_IGNORE, ierr)
		end do
	case default
	end select

	!Close file on all processors
	call MPI_FILE_CLOSE(restartfileid, ierr)

	!This barrier is needed inn order to get the correct file size in the next write
	call MPI_Barrier(MD_COMM, ierr)	
 
    !call MPI_FILE_CLOSE(restartfileid, ierr)


	!----------------Write header------------------------
	!Written at the end for performance and simplicity reasons 
	!(See Gropp, lusk & Thakur Using MPI-2)

	!Write the header with one processor only
	if (irank .eq. iroot) then

        call MPI_FILE_OPEN(MPI_COMM_SELF,trim(prefix_dir)//'results/final_state', & 
			MPI_MODE_WRONLY, MPI_INFO_NULL, restartfileid, ierr)
                
        if (ierr /= 0) then 
                write(0,*) "MD parallel_io: error in MPI_File open"
        endif

        call MPI_File_get_size(restartfileid,filesize,ierr)
        
        disp = filesize

        call MPI_FILE_SET_VIEW(restartfileid, disp, MPI_BYTE, & 
 			MPI_BYTE, 'native', MPI_INFO_NULL, ierr)

        call MPI_File_write(restartfileid,sum(procnp)   ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,initialnunits ,3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,Nsteps        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,tplot         ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,seed          ,2,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,periodic      ,3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,potential_flag,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,nmonomers     ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		call MPI_File_write(restartfileid,npx           ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		call MPI_File_write(restartfileid,npy           ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		call MPI_File_write(restartfileid,npz           ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

		call MPI_File_write(restartfileid,density       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,rcutoff       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,delta_t       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,elapsedtime   ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,k_c           ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,R_0           ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,delta_rneighbr,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

        header_pos = filesize ! just in case offset kind is 32 bit, rather improbable these days  !!!
        call MPI_File_write(restartfileid,header_pos,1,MPI_INTEGER8,MPI_STATUS_IGNORE,ierr)
        call MPI_FILE_CLOSE(restartfileid, ierr)
   
		call MPI_FILE_OPEN(MPI_COMM_SELF,trim(prefix_dir)//'results/final_state', & 
		                   MPI_MODE_RDONLY, MPI_INFO_NULL, restartfileid, ierr)
        call MPI_File_get_size(restartfileid,filesize,ierr)
        call MPI_File_close(restartfileid, ierr)

	endif

end subroutine parallel_io_final_state

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer							:: procdisp
	integer							:: i, datasize
	real,dimension(np)				:: Xbuf, Ybuf, Zbuf
	real,dimension(globalnp)        :: Xbufglob
	real,dimension(globalnp)        :: Ybufglob
	real,dimension(globalnp)        :: Zbufglob
	integer                         :: n,globmolno
	integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

	Xbufglob = 0.0           !Initialise to zero so that global array...
	Ybufglob = 0.0           !...may be found by summation in parallel
	Zbufglob = 0.0           !------------------------------------------

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

	!Determine size of real datatype
 	call MPI_type_size(MPI_real,datasize,ierr)

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre
	select case(potential_flag)
	case(0)
		Xbuf(:) = r(1:np,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
		Ybuf(:) = r(1:np,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
		Zbuf(:) = r(1:np,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

		procdisp = 0
		!Obtain displacement of each processor using all other procs' np
		do i = 1, irank -1
			procdisp = procdisp + procnp(i)*datasize
		enddo

		!Open file on all processors
		call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/vmd_temp.dcd', & 
			MPI_MODE_RDWR + MPI_MODE_CREATE, & 
			MPI_INFO_NULL, fileid, ierr)

		!-------------Write X coordinates--------------------

		!Obtain location to write in file
		disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
			+ procdisp				  	!Processor location

		call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
			MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
		!Write information to file
		call MPI_FILE_WRITE_ALL(fileid, Xbuf, np, MPI_REAL, & 
			MPI_STATUS_IGNORE, ierr) 

		!-------------Write Y coordinates--------------------

		!Obtain location to write in file
		disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
			+ procdisp &				  	!Processor location
			+ globalnp * datasize			  	!Y Coordinate location

		call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
			MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
		!Write information to file
		call MPI_FILE_WRITE_ALL(fileid, Ybuf, np, MPI_REAL, & 
			MPI_STATUS_IGNORE, ierr) 

		!-------------Write Z coordinates--------------------

		!Obtain location to write in file
		disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
			+ procdisp &				  	!Processor location
			+ 2 * globalnp * datasize		  	!Z Coordinate location

		call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
			MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
		!Write information to file
		call MPI_FILE_WRITE_ALL(fileid, Zbuf, np, MPI_REAL, & 
			MPI_STATUS_IGNORE, ierr) 
	
		call MPI_FILE_CLOSE(fileid, ierr) 
	
	case(1)

		!Build sparse individual "global" buffers according to global molecular ID of each monomer
		do n=1,np
			globmolno           = monomer(n)%glob_no
			Xbufglob(globmolno) = r(n,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
			Ybufglob(globmolno) = r(n,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
			Zbufglob(globmolno) = r(n,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
		end do

		call globalSumVectReal(Xbufglob,globalnp)  !Global summation to complete global buffer
		call globalSumVectReal(Ybufglob,globalnp)
		call globalSumVectReal(Zbufglob,globalnp)
	
		call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/vmd_temp.dcd', & 
		                   MPI_MODE_RDWR + MPI_MODE_CREATE,       & 
		                   MPI_INFO_NULL, fileid, ierr)

		disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize   !Current iteration

		!Write X positions---------------------------------------------
		call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, MPI_REAL,      &      !Find position
		                       'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(fileid, Xbufglob, globalnp, MPI_REAL, &      !Write buffer
		                    MPI_STATUS_IGNORE, ierr)
		disp = disp + globalnp*datasize                                      !Update file disp

		!Write Y positions---------------------------------------------
		call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, MPI_REAL,      &        
		                       'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(fileid, Ybufglob, globalnp, MPI_REAL, &
		                    MPI_STATUS_IGNORE, ierr)
		disp = disp + globalnp*datasize

		!Write Z positions---------------------------------------------
		call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, MPI_REAL,      &
		                       'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(fileid, Zbufglob, globalnp, MPI_REAL, &
		                    MPI_STATUS_IGNORE, ierr)
			
		call MPI_FILE_CLOSE(fileid, ierr) 
	
	case default
	end select


end subroutine parallel_io_vmd

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd_sl
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer							:: procdisp
	integer							:: i,n, datasize
	real,dimension(np)				:: Xbuf, Ybuf, Zbuf
	integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

	!Determine size of real datatype
 	call MPI_type_size(MPI_real,datasize,ierr)

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre
	do n = 1, np
		select case(tag(n))
		case(0, 4) 	!Liquid Molecules
			Xbuf(:) = r(:,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
			Ybuf(:) = r(:,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
			Zbuf(:) = r(:,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
		case default	!Solid molecules
			Xbuf(n) = -halfdomain(1)
			Ybuf(n) = -halfdomain(2)
			Zbuf(n) = -halfdomain(3)
		end select
	enddo

	procdisp = 0
	!Obtain displacement of each processor using all other procs' np
	do i = 1, irank -1
		procdisp = procdisp + procnp(i)*datasize
	enddo

	!================================
	!    Write liquid Molecules	=
	!================================
	!Open file on all processors
	call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/vmd_liquid_temp.dcd', & 
		MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		MPI_INFO_NULL, fileid, ierr)

	!-------------Write X coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp				  	!Processor location

	!print*, irank, 'x disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Xbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Y coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp &				  	!Processor location
		+ globalnp * datasize			  	!Y Coordinate location

	!print*, irank, 'y disp', disp
	
	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Ybuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Z coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp &				  	!Processor location
		+ 2 * globalnp * datasize		  	!Z Coordinate location

	!print*, irank, 'z disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Zbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!Close file on all processors
	call MPI_FILE_CLOSE(fileid, ierr) 

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre
	do n = 1, np
		select case(tag(n))
		case(0, 4) 	!Liquid Molecules
			Xbuf(n) = -halfdomain(1)
			Ybuf(n) = -halfdomain(2)
			Zbuf(n) = -halfdomain(3)
		case default	!Solid molecules
			Xbuf(n) = r(n,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
			Ybuf(n) = r(n,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
			Zbuf(n) = r(n,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
		end select
	enddo

	!================================
	!    Write Solid Molecules	=
	!================================

	!Open file on all processors
	call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/vmd_solid_temp.dcd', & 
		MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		MPI_INFO_NULL, fileid, ierr)

	!-------------Write X coordinates--------------------
		
	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp				  	!Processor location

	!print*, irank, 'x disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Xbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Y coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp &				  	!Processor location
		+ globalnp * datasize			  	!Y Coordinate location

	!print*, irank, 'y disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Ybuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Z coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp &				  	!Processor location
		+ 2 * globalnp * datasize		  	!Z Coordinate location

	!print*, irank, 'z disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Zbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!Close file on all processors
	call MPI_FILE_CLOSE(fileid, ierr) 

end subroutine parallel_io_vmd_sl

!------------------------------------------------------------------------
!Write positions of molecules to a file using one single write instead of 3

subroutine parallel_io_vmd_optimised
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer							:: procdisp
	integer							:: i, datasize
	real,dimension(np*3)			:: buf
	integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

	!Determine size of real datatype
 	call MPI_type_size(MPI_real,datasize,ierr)

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre

	!buf(1:np) = r(:,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
	!buf((np+1):(2*np)) = r(:,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
	!buf((2*np+1):(3*np)) = r(:,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

	do i=1,np
		buf(3*(i-1)+1) = r(i,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
		buf(3*(i-1)+2) = r(i,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
		buf(3*(i-1)+3) = r(i,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
	enddo

	procdisp = 0
	!Obtain displacement of each processor using all other procs' np
	do i = 1, irank -1
		procdisp = procdisp + procnp(i)*datasize
	enddo

	!Open file on all processors
	call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/vmd_temp.dcd', & 
		MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		MPI_INFO_NULL, fileid, ierr)

	!-------------Write X coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp				  	!Processor location

	!print*, irank, 'x disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, buf, np*3, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!Close file on all processors
	call MPI_FILE_CLOSE(fileid, ierr) 

end subroutine parallel_io_vmd_optimised

!------------------------------------------------------------------------
!Write positions of molecules in halo to a file

subroutine parallel_io_vmd_halo
	use module_parallel_io
	implicit none

	call error_abort("Cannot print vmd halos in parallel simulation - run in serial")

end subroutine parallel_io_vmd_halo

!---------------------------------------------------------------------------------
! Write value of last output iteration

subroutine update_simulation_progress_file
	use module_parallel_io
	implicit none

	if (irank .eq. iroot) then
		open (unit=99999, file=trim(prefix_dir)//"results/simulation_progress")
		write(99999,*) iter
		close(99999,status='keep')
	endif

end subroutine update_simulation_progress_file


!---------------------------------------------------------------------------------
! Record mass in a slice through the domain

subroutine mass_slice_io(ixyz)
	use module_parallel_io
	use calculated_properties_MD
	use messenger
	use interfaces
	implicit none

	integer							:: ixyz,jxyz,kxyz
	integer							:: slicefileid, int_datasize
	integer(kind=MPI_OFFSET_KIND)   :: disp

	!Get two directions orthogonal to slice direction
	kxyz = mod(ixyz,3)+1
	jxyz = mod(ixyz+1,3)+1

	!Sum over all bins using directional sub communicators and gather on {ijk}block=1
	call SubcommSum(slice_mass, nbins(ixyz), jxyz)
	call SubcommSum(slice_mass, nbins(ixyz), kxyz)

	!Only root processor in each directional subcomm writes data
	if (icoord(jxyz,irank) .eq. 1 .and. icoord(kxyz,irank) .eq. 1) then

		!Determine size of datatypes
		call MPI_type_size(MPI_Integer,int_datasize,ierr)

		!Only processors on directional subcomm write
		call MPI_FILE_OPEN(icomm_xyz(ixyz), trim(prefix_dir)//'./results/mslice', & 
				   			MPI_MODE_WRONLY+ MPI_MODE_CREATE , & 
				   			MPI_INFO_NULL, slicefileid, ierr)

		!Obtain displacement of current record
		disp =   (iter/(tplot*Nmass_ave) - 1) 	  	&		!Current iteration
		       * globalnbins(ixyz)*int_datasize 	&		!Record size
		       + nbins(ixyz)*int_datasize*(jblock-1)		!Processor location

		call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_INTEGER, & 
	 				MPI_INTEGER, 'native', MPI_INFO_NULL, ierr)

		call MPI_FILE_WRITE_ALL(slicefileid,slice_mass,nbins(ixyz), MPI_INTEGER, & 
						MPI_STATUS_IGNORE, ierr)

		call MPI_FILE_CLOSE(slicefileid, ierr)

	endif 

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
	write(filename, '(a9,a4)' ) trim(prefix_dir)//'results/m', io_type

	!Include halo surface fluxes to get correct values for all cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!do n = 1, nsurfacebins
	!	i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
		!Change in number of Molecules in halo cells
	!	CV_mass_out(modulo((i-2),nbins(1))+2, & 
	!		    modulo((j-2),nbins(2))+2, & 
	!		    modulo((k-2),nbins(3))+2) = & 
	!		    CV_mass_out(modulo((i-2),nbins(1))+2, & 
	!				modulo((j-2),nbins(2))+2, & 
	!				modulo((k-2),nbins(3))+2) &
	!					+ CV_mass_out(i,j,k)
	!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (io_type .eq. 'snap') then
		m = iter/(Nmflux_ave) + 1 !Initial snapshot taken
	else
		m = iter/(tplot*Nmass_ave)
	endif

	!Write mass to file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!inquire(iolength=length) CV_mass_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1)
	!open (unit=5,file=filename,form="unformatted",access='direct',recl=length)
	!write(5,rec=m) CV_mass_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1)
	!close(5,status='keep')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine mass_bin_io

!---------------------------------------------------------------------------------
! Record velocity in a slice through the domain

subroutine velocity_slice_io(ixyz)
	use module_parallel_io
	use calculated_properties_MD
	use messenger
	implicit none

	integer							:: ixyz,jxyz,kxyz
	integer							:: slicefileid, dp_datasize
	integer,dimension(3)			:: idims
	integer(kind=MPI_OFFSET_KIND)   :: disp

	!Write mass
	call mass_slice_io(ixyz)

	!Get two directions orthogonal to slice direction
	kxyz = mod(ixyz,3)+1
	jxyz = mod(ixyz+1,3)+1
	idims(1) = npx; idims(2) = npy; idims(3) = npz

	!Sum over all bins using directional sub communicators and gather on root
	call SubcommSum(slice_momentum(:,1), nbins(ixyz), jxyz)
	call SubcommSum(slice_momentum(:,1), nbins(ixyz), kxyz)
	call SubcommSum(slice_momentum(:,2), nbins(ixyz), jxyz)
	call SubcommSum(slice_momentum(:,2), nbins(ixyz), kxyz)
	call SubcommSum(slice_momentum(:,3), nbins(ixyz), jxyz)
	call SubcommSum(slice_momentum(:,3), nbins(ixyz), kxyz)

	!Only root processor in each directional subcomm writes data
	if (icoord(jxyz,irank) .eq. 1 .and. icoord(kxyz,irank) .eq. 1) then

		call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

		!Only processors on directional subcomm write
		call MPI_FILE_OPEN(icomm_xyz(ixyz), trim(prefix_dir)//'./results/vslice', & 
				   MPI_MODE_WRONLY + MPI_MODE_CREATE , & 
				   MPI_INFO_NULL, slicefileid, ierr)

		!Obtain displacement of x record
		disp =   (iter/(tplot*Nmass_ave) - 1) 	  	&	!Current iteration
		       * nd*globalnbins(ixyz)*dp_datasize 	&	!times record size
		       + nbins(ixyz)*dp_datasize*(jblock-1)		!Processor location

		call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
	 				MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,1),nbins(ixyz), & 
					MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

		!Obtain displacement of y record
		disp =   (iter/(tplot*Nmass_ave) - 1) 	  	&	!Current iteration
		       * nd*globalnbins(ixyz)*dp_datasize 	&	!Record size
		       + nbins(ixyz)*dp_datasize*(jblock-1)	&	!Processor location
		       + nbins(ixyz)*dp_datasize*idims(ixyz)		!after x data 

		call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
	 				MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,2),nbins(ixyz), & 
					MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

		!Obtain displacement of z record
		disp =   (iter/(tplot*Nmass_ave) - 1) 	  	&	!Current iteration
		       * nd*globalnbins(ixyz)*dp_datasize 	&	!Record size
		       + nbins(ixyz)*dp_datasize*(jblock-1)	&	!Processor location
		       + 2*nbins(ixyz)*dp_datasize*idims(ixyz)		!after x & y data 

		call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
	 				MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,3),nbins(ixyz), & 
					MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

		call MPI_FILE_CLOSE(slicefileid, ierr)

	endif 

end subroutine velocity_slice_io

!------------------------------------------------------------------------
!A large scale routine with each proc writing its own bins in binary
!Write velocity slice information to a file


subroutine parallel_slice_io_large_scale
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer								:: slicefileid, int_datasize, dp_datasize
	integer(kind=MPI_OFFSET_KIND)   	:: disp!, resultsize

	!Determine size of datatypes
	call MPI_type_size(MPI_Integer,int_datasize,ierr)
	call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

	call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/vslice', & 
			MPI_MODE_WRONLY , & 
			MPI_INFO_NULL, slicefileid, ierr)

	!WRITE BINNED VELOCITY SUMMATION

	!Obtain displacement of current record
	disp =(iter/real((tplot*Nvel_ave),kind(0.d0))-1) *       &	!Current iteration
	      globalnbins(1) *(int_datasize+3*dp_datasize) & 		!Current iteration	
	      + nbins(1)*(3*dp_datasize)*(jblock-1)			!Processor location
	!Set current processes view
	call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_double_precision, & 
 			MPI_double_precision, 'native', MPI_INFO_NULL, ierr)
	!Write data
	call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,:),nd*nbins(1), MPI_double_precision, & 
			MPI_STATUS_IGNORE, ierr) !Velocity bins

	!WRITE BINNED MOLECULAR COUNT

	!Obtain displacement of current record
	disp =(iter/real((tplot*Nvel_ave),kind(0.d0))-1) *		& !Current iteration
	      globalnbins(1) * (int_datasize+dp_datasize)  & !Current iteration	
		+ nbins(1)*(int_datasize)*(jblock-1) 	& !Processor location
		+ globalnbins(1)* dp_datasize 		  !After binned velocity

	!Set current processes view
	call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_Integer, & 
				MPI_Integer, 'native', MPI_INFO_NULL, ierr)

	!Write data
	call MPI_FILE_WRITE_ALL(slicefileid,slice_mass(:),nbins(1), MPI_Integer, & 
			MPI_STATUS_IGNORE, ierr) !molecular tally bins

	call MPI_FILE_CLOSE(slicefileid, ierr) 

end subroutine parallel_slice_io_large_scale

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Work out correct filename for i/o type
	!write(filename, '(a9,a4)' ) 'results/v', io_type

	!---------------Correct for surface fluxes on halo cells---------------
	!Include halo surface fluxes to get correct values for all cells
	!do n = 1, nsurfacebins
	!	i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  

	!	!Change in Momentum in halo cells
	!	CV_momentum_out(modulo((i-2),nbins(1))+2, & 
	!		      	modulo((j-2),nbins(2))+2, & 
	!		      	modulo((k-2),nbins(3))+2,:) = & 
	!			CV_momentum_out(modulo((i-2),nbins(1))+2,& 
	!					modulo((j-2),nbins(2))+2,&
	!					modulo((k-2),nbins(3))+2,:) & 
	!						+ CV_momentum_out(i,j,k,:)
	!enddo

	if (io_type .eq. 'snap') then
		!CV_momentum_out = CV_momentum_out / (tplot*Nvflux_ave)
		m = iter/(Nvflux_ave) + 1 !Initial snapshot taken
	else
		!CV_momentum_out = CV_momentum_out / (tplot*Nvel_ave)
		m = iter/(tplot*Nvel_ave)
	endif
	!Write velocity to file
	!inquire(iolength=length) CV_momentum_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	!open (unit=6, file=filename,form="unformatted",access='direct',recl=length)
	!write(6,rec=m) CV_momentum_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	!close(6,status='keep')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine velocity_bin_io

!---------------------------------------------------------------------------------
! Record temperature in a slice through the domain

subroutine temperature_slice_io(ixyz)
use module_parallel_io
use calculated_properties_MD
implicit none

	integer		:: ixyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine temperature_slice_io

!---------------------------------------------------------------------------------
! Record temperature in 3D bins throughout domain

subroutine temperature_bin_io(CV_mass_out,CV_temperature_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer					:: CV_mass_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	double precision		:: CV_temperature_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	character(4)			:: io_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine temperature_bin_io



!---------------------------------------------------------------------------------
! Record velocity in 3D bins throughout domain

subroutine energy_bin_io(CV_energy_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	double precision		:: CV_energy_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	character(4)			:: io_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine energy_bin_io


!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine virial_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

end subroutine virial_stress_io

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine VA_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!call intergrate_trap(Pxycorrel,tplot*delta_t,Nstress_ave,viscosity)

	!viscosity = (viscosity*volume)/(3.0*Nstress_ave*Nvisc_ave*inputtemperature)

	!Write viscosity to file
	!m = iter/(tplot*Nstress_ave*Nvisc_ave)
	!inquire(iolength=length) viscosity
	!open (unit=7, file="results/visc",form="unformatted",access='direct',recl=length)
	!write(7,rec=m) viscosity
	!close(7,status='keep')

	!Pxycorrel = 0.d0	!Reset Pxycorrel to zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Include halo surface fluxes to get correct values for all cells
	!do n = 1, nsurfacebins
	!	i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
		!Flux over halo cells
	!	mass_flux(modulo((i-2),nbins(1))+2,modulo((j-2),nbins(2))+2,modulo((k-2),nbins(3))+2,:) = & 
	!	mass_flux(modulo((i-2),nbins(1))+2,modulo((j-2),nbins(2))+2,modulo((k-2),nbins(3))+2,:) + mass_flux(i,j,k,:)
	!enddo

	!Write six CV surface mass fluxes to file
	!m = iter/(Nmflux_ave)
	!inquire(iolength=length) mass_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	!open (unit=8, file="results/mflux",form="unformatted",access='direct',recl=length)
	!write(8,rec=m) mass_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	!close(8,status='keep')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine mass_flux_io

!---------------------------------------------------------------------------------
! Record momentum fluxes accross surfaces of Control Volumes

subroutine momentum_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: ixyz,i,j,k,n,m,length
	double precision		:: binface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Include halo surface fluxes to get correct values for all cells
	!do n = 1, nsurfacebins
	!	i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
		!Flux over halo cells
	!	momentum_flux(	modulo((i-2),nbins(1))+2, & 
	!		      	modulo((j-2),nbins(2))+2, & 
	!		      	modulo((k-2),nbins(3))+2,:,:) = & 
	!			momentum_flux(	modulo((i-2),nbins(1))+2,& 
	!					modulo((j-2),nbins(2))+2,&
	!					modulo((k-2),nbins(3))+2,:,:) & 
	!						+ momentum_flux(i,j,k,:,:)
	!enddo

	!do ixyz = 1,3
	!	binface	      = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
	!		     	(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
	!	momentum_flux(:,:,:,:,ixyz  )=momentum_flux(:,:,:,:,ixyz  )/(binface) !Bottom
	!	momentum_flux(:,:,:,:,ixyz+3)=momentum_flux(:,:,:,:,ixyz+3)/(binface) !Top
	!enddo

	!momentum_flux = momentum_flux/(delta_t*Nvflux_ave)
	!Write momnetum flux to file
	!m = iter/(Nvflux_ave)
	!inquire(iolength=length) momentum_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	!open (unit=9, file="results/vflux",form="unformatted",access='direct',recl=length)
	!write(9,rec=m) momentum_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	!close(9,status='keep')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine momentum_flux_io

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MOP_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none


end subroutine MOP_stress_io

!---------------------------------------------------------------------------------
! Record stress accross surfaces of Control Volumes

subroutine surface_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: ixyz,i,j,k,n,m,length
	double precision,dimension(3)	:: binface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Include halo surface stresses to get correct values for all cells
	!do n = 1, nsurfacebins
	!	i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
		!Set Stresses to value of halo cells
	!	Pxyface(	modulo((i-2),nbins(1))+2, & 
	!		      	modulo((j-2),nbins(2))+2, & 
	!		      	modulo((k-2),nbins(3))+2,:,:) = & 
	!			Pxyface(	modulo((i-2),nbins(1))+2,& 
	!					modulo((j-2),nbins(2))+2,&
	!					modulo((k-2),nbins(3))+2,:,:) & 
	!						     + Pxyface(i,j,k,:,:)

	!enddo


	!do ixyz = 1,3
	!	binface(ixyz) = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
	!		     	(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
	!	Pxyface(:,:,:,:,ixyz  ) = 0.25d0 * Pxyface(:,:,:,:,ixyz  )/binface(ixyz) !Bottom
	!	Pxyface(:,:,:,:,ixyz+3) = 0.25d0 * Pxyface(:,:,:,:,ixyz+3)/binface(ixyz) !Top
	!enddo

	!Integration of stress using trapizium rule requires multiplication by timestep
	!Pxyface = delta_t*Pxyface!/Nvflux_ave

	!Write surface pressures to file
	!m = iter/(Nvflux_ave)
	!inquire(iolength=length) Pxyface(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	!open (unit=9, file="results/psurface",form="unformatted",access='direct',recl=length)
	!write(9,rec=m) Pxyface(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	!close(9,status='keep')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine surface_stress_io

!---------------------------------------------------------------------------------
! Record energy fluxes accross surfaces of Control Volumes

subroutine energy_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine energy_flux_io

!---------------------------------------------------------------------------------
! Record  energy accross plane

subroutine MOP_energy_io(ixyz)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer		:: ixyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine MOP_energy_io


!---------------------------------------------------------------------------------
! Record stress times velocity (power) accross surfaces of Control Volumes

subroutine surface_power_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine surface_power_io


!-----------------------------------------------------------------------------
! Write macroscopic properties to file
!-----------------------------------------------------------------------------
subroutine macroscopic_properties_header
use module_parallel_io
use calculated_properties_MD
implicit none

	if (irank .eq. iroot) then	
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
	endif

end subroutine macroscopic_properties_header


subroutine macroscopic_properties_record
use module_parallel_io
use calculated_properties_MD
implicit none

	if (irank .eq. iroot) then
		if (potential_flag.eq.0) then	
			write(10,'(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.4)') &
			iter,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy,';',totenergy,';',pressure
		else if (potential_flag.eq.1) then
			write(10,'(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.4)') &
			iter,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy_LJ,';',potenergy_FENE,';',potenergy,';',totenergy,';',pressure
		end if
	endif

end subroutine macroscopic_properties_record

!-----------------------------------------------------------------------------
! Write end-to-end vector time correlation function
!-----------------------------------------------------------------------------
subroutine etev_io
	use module_parallel_io
	implicit none

	integer :: etev_fid

	if (irank .eq. iroot) then
	
		if (iter.eq.etevtcf_iter0) then
			call MPI_FILE_DELETE(trim(prefix_dir)//'results/etev', MPI_INFO_NULL, ierr)
			call MPI_FILE_OPEN(MPI_COMM_SELF, trim(prefix_dir)//'results/etev', &
			                   MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, etev_fid, ierr)
			call MPI_FILE_WRITE(etev_fid,etev,nchains*nd,MPI_DOUBLE_PRECISION, &
			                    MPI_STATUS_IGNORE, ierr)
		else if (iter.gt.etevtcf_iter0) then
			call MPI_FILE_OPEN(MPI_COMM_SELF, trim(prefix_dir)//'results/etev', &
			                   MPI_MODE_WRONLY + MPI_MODE_APPEND, MPI_INFO_NULL, etev_fid, ierr)
			call MPI_FILE_WRITE(etev_fid,etev,nchains*nd,MPI_DOUBLE_PRECISION, &
			                    MPI_STATUS_IGNORE, ierr)
		end if
		
		call MPI_FILE_CLOSE(etev_fid,ierr)

	endif

end subroutine etev_io

!Write radius of gyration to output file
subroutine r_gyration_io
use module_parallel_io
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine r_gyration_io
