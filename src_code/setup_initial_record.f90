!--------------------------------------------------------------------------------
!
!                                Initial Record
! Prints the headers for the output table
!
!--------------------------------------------------------------------------------

module module_initial_record

	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use calculated_properties_MD

end module module_initial_record
!----------------------------------------------------------------------------------

subroutine setup_initial_record
	use interfaces
	use module_initial_record
	use polymer_info_MD
	use shear_info_MD
	use concentric_cylinders, only: gcpol_bins
	use librarymod, only : get_version_number
	implicit none

	integer					:: i
	logical 				:: file_exist
	character				:: ixyz_char
	character(8)			:: the_date
	character(10)			:: the_time
	character(23),parameter :: file_names(25) = &
								(/ "mslice      ", "mbins       ", "msnap   ",&
								   "vslice      ", "vbins       ", "vsnap   ",&
								   "pvirial     ", "pVA         ", "pVA_k   ",& 
								   "pVA_c       ", "visc        ", "mflux   ",& 
								   "vflux       ", "pplane      ", "psurface",&
								   "esnap       ", "eflux       ", "eplane  ",&
								   "esurface    ", "viscometrics", "rdf     ",&
	                               "rdf3d       ", "ssf         ", "Fext    ",&
								   "Tbins       " /) 


	!Delete all files from previous run
	if (irank.eq.iroot) then
		do i=1,size(file_names)
			inquire(file=trim(prefix_dir)//'results/'//file_names(i),exist=file_exist)
			if(file_exist) then
				open (unit=23, file=trim(prefix_dir)//'results/'//file_names(i))
				close(23,status='delete')
			endif
		enddo
	endif
	call messenger_syncall()

	!Evaluate system properties on all processes
	call simulation_compute_forces
	call evaluate_macroscopic_properties

	if (potential_flag.eq.1) then
		select case(etevtcf_outflag)
		case (1)
			call etevtcf_calculate_parallel
		case (2)
			call etevtcf_calculate_parallel
			call etev_io
		case default
		end select
		
		select case(r_gyration_outflag)
		case (1)
			call r_gyration_calculate_parallel
		case (2)
			call r_gyration_calculate_parallel
			call r_gyration_io
		case default
		end select
	end if

	if (vmd_outflag.ne.0 .and. potential_flag.eq.1) call build_psf
	
	!Calculate Control Volume starting state
	call initial_control_volume

	if (irank .eq. iroot) then

		call date_and_time(the_date, the_time)
		print*, "              ==================================="
		print*, "              |    _____      _                 |"
		print*, "              |   /  ___|    | |                |"
		print*, "              |   \ `--.  ___| |_ _   _ _ __    |"
		print*, "              |    `--. \/ _ \ __| | | | '_ \   |"
		print*, "              |   /\__/ /  __/ |_| |_| | |_) |  |"
		print*, "              |   \____/ \___|\__|\__,_| .__/   |"
		print*, "              |                        | |      |"
		print*, "              |                        |_|      |"
		print*, "              ==================================="

		!Display all parameters required to describe simulation
		print*, 'Simulation run on Date: ', the_date
		print*, 'Simulation start time: ', the_time
		print*, 'Subversion revision number: ', get_version_number()
		print*, '================= Molecular Simulation Parameters ===================='
		print*, 'Number of Dimensions: ', nd
		print*, 'Number of Particles: ', globalnp
		print*, 'Time Step - delta t: ',  delta_t
		print*, 'Total number of steps: ',  Nsteps - initialstep
		select case(integration_algorithm)
		case(0)
			print*, 'Integration algorithm: Leapfrog-Verlet'
		case(1)
			print*, 'Integration algorithm: Velocity-Verlet'
		end select
		select case(potential_flag)
		case(0)
			print*, 'Interatomic potential: LJ'
		case(1)
			print*, 'Interatomic potential: LJ + FENE'
		end select
		select case(ensemble)
		case(0)
			print*, 'NVE ensemble'
		case(1)
			print*, 'NVT (Nosé-Hoover thermostat)'
		case(2)
			print*, 'NVT (Gaussian iso-kinetic thermostat) - only availble with VV'
		case(3)
			print*, 'NVT (Profile unbiased Nosé-Hoover thermostat)'
		case(4)
			print*, 'NVT (Pairwise additive Nosé-Hoover thermostat by Allen & Schmid)'
		case(5)
			print*, 'NVT (DPD thermostat by Soddemann)' 
		case(6)
			print*, 'Tagged molecule system:' 
			print'(a,3f10.5)', ' Distance from bottom of Fixed Molecules in x,y and z:', 	fixdistbottom
			print'(a,3f10.5)', ' Distance from top of Fixed Molecules in x,y and z:', 	fixdisttop
			print'(a,3f10.5)', ' Distance from bottom of Tethered Molecules in x,y and z:', tethereddistbottom
			print'(a,3f10.5)', ' Distance from top of Tethered Molecules in x,y and z:', 	tethereddisttop
			print'(a,3f10.5)', ' Distance from bottom of Sliding Molecules in x,y and z:', 	slidedistbottom
			print'(a,3f10.5)', ' Distance from top of Sliding Molecules in x,y and z:', 	slidedisttop
			print'(a,3f10.5)', ' Velocity of Sliding Molecules in x,y and z:', 	wallslidev
			print'(a,3f10.5)', ' Distance from bottom of NH Themostatted Molecules in x,y and z:', 	thermstatbottom
			print'(a,3f10.5)', ' Distance from top of NH Themostatted Molecules in x,y and z:', 	thermstattop
		end select

		if (external_force_flag .ne. 0) then
			if ( F_ext_ixyz .eq. 1) then
				ixyz_char = 'x'
			elseif (F_ext_ixyz .eq. 2) then 
				ixyz_char = 'y'
			else 
				ixyz_char = 'z'
			endif
			if (external_force_flag .eq. 1) then
				print'(a,f10.5,3a)', ' External force of magnitude ', F_ext, ' applied in the ', ixyz_char, ' direction to all molecules'
			elseif (external_force_flag .eq. 2) then
				print'(a,f10.5,3a)', ' External force of magnitude ', F_ext, ' applied in the ', ixyz_char, ' direction to '
				print'(a,6f10.5)', ' molecules between ', F_ext_limits
			endif

		endif

 		select case(force_list) 
		case(0)
			print*,'All pairs force calculation (using Newton`s 3rd law)'
		case(1)
			print*,'Cell list force calculation'
		case(2)
			print*,'Neighbour list force calculation'
		case(3)
			print*,'Neighbour list force calculation (using Newton`s 3rd law)'
		end select
		print*, 'Starting step of simulation:', initialstep
		print*, 'Generate output file every: ',  tplot, 'steps'
		print*, 'Density: ',              density
		print*, 'Initial Temperature: ',  inputtemperature
		if (fixed_rebuild_flag .eq. 0) then
			print'(a,f19.15,a,f10.5)', ' Cut off distance:  ', rcutoff, &
				    '  Neighbour List Delta r:  ', delta_rneighbr
		else
			print'(a,f19.15a,f8.4,a,i5,a)', ' Cut off ', rcutoff, ' nbr: ', delta_rneighbr, &
				    ' & fixed rebuild every ', fixed_rebuild, ' steps'
		endif
		print*, 'Initial unit size (FCC unit) in x,y and z:'
		print*,	initialunitsize(1), initialunitsize(2), initialunitsize(3)
		print*, 'Domain in x,y and z: '
		print*,  globaldomain(1), globaldomain(2), globaldomain(3)
		print*, 'Domain volume: ', volume
		print'(a,3i8)', ' Periodic Boundary Conditions in x,y and z:', periodic
		if (any(periodic.gt.1)) then
			print*, 'Lees-Edwards sliding boundary conditions ON'
			print*, 'Shear velocity: ', le_sv
			print*, 'Shear rate: ', le_sr
			print*, 'Shear plane: ', le_sp
			print*, 'Shear direction: ', le_sd
			print*, 'Shear iter0: ', le_i0
		end if

		print*, '==================== Computational Parameters ========================='
		print'(a,3i8)', ' Domain split into computational cells in x,y and z:', & 
					 ncells(1)*npx, ncells(2)*npy, ncells(3)*npz
		print'(a,3f10.5)', ' Each of size:', cellsidelength
		print*, 'Number of processors in x,y and z: ', npx, npy, npz
		print*, 'Cells per Processor:', ncells(1), ncells(2), ncells(3)
		print*, 'Cells per Processor including Halos: ', ncells(1)+2*nh, ncells(2)+2*nh, ncells(3)+2*nh
		print*, 'Random seed used for initial velocity generation:', seed
		print'(a,3(i4,a))', ' Rescue backup file saved every ', floor(rescue_snapshot_freq/3600.d0), ' hours ', &
				 floor(mod(rescue_snapshot_freq,3600.d0)/60.d0), ' minutes and ', floor(mod(rescue_snapshot_freq,60.d0)), ' seconds'
		print*, '======================== Output Parameters ============================'
	
		select case(vmd_outflag)
		case(0)
			print*, 'VMD output off'
		case(1)
			print*, 'VMD output every:', tplot, 'iterations'
		case(2)
			print*, 'VMD output with seperate solid/liquid every:', tplot, 'iterations'
		case(3)
			print*, 'VMD output Domain+halos every:', tplot, 'iterations'
		case(4)
			print*, 'VMD output true unwrapped positions every:', tplot, 'iterations'
		case default
			call error_abort("Invalid VMD output flag in input file")
		end select

		select case(macro_outflag)
		case(0)
			print*, 'No Macroscopic Properties printed to screen'
		case(1,3)
			print*, 'Macroscopic properties printed to screen every:', tplot, 'iterations'
		case(2,4)
			call macroscopic_properties_header
			print*, 'Macroscopic properties printed to results/macroscopic_properties every:', tplot, 'iterations.'
		case default
			call error_abort("Invalid Macroscopic properties output flag in input file")
		end select

		select case(mass_outflag)
		case(0)
			print*, 'Mass record off'
			print*, ''
			print*, ''
		case(1:3)
			if (mass_outflag .eq. 1) then
				ixyz_char = 'x'
			elseif (mass_outflag .eq. 2) then 
				ixyz_char = 'y'
			else 
				ixyz_char = 'z'
			endif

			print'(3(a,i8),a)', ' Mass slice recorded every:', & 
					tplot,' x ',Nmass_ave,' = ',tplot*Nmass_ave,' iterations'
			print'(a,i8,2a)', ' Domain split into',gnbins(mass_outflag) ,' mass Averaging Slices in  ', ixyz_char  
			print'(a,f10.5)', ' With each averaging slice of width:', & 
						globaldomain(mass_outflag)/gnbins(mass_outflag)
		case(4)
			print'(3(a,i8),a)', ' mass 3D bins recorded every:', &
					tplot,' x ',Nvel_ave,' = ',tplot*Nvel_ave,' iterations'
			print'(a,3i8)', ' Domain split into mass Averaging Bins in x,y and z:', gnbins
			print'(a,3f10.5)', ' Each of size:', & 
			globaldomain(1)/gnbins(1), globaldomain(2)/gnbins(2),globaldomain(3)/gnbins(3)
		case default
			call error_abort("Invalid Mass output flag in input file")
		end select

		select case(velocity_outflag)
		case(0)
			print*, 'Velocity record off'
			print*, ''
			print*, ''
		case(1:3)
			if (velocity_outflag .eq. 1) then
				ixyz_char = 'x'
			elseif (velocity_outflag .eq. 2) then 
				ixyz_char = 'y'
			else 
				ixyz_char = 'z'
			endif

			print'(3(a,i8),a)', ' Velocity slice recorded every:', & 
					tplot,' x ',Nvel_ave,' = ',tplot*Nvel_ave,' iterations'
			print'(a,i8,2a)', ' Domain split into',gnbins(velocity_outflag) , & 
					  ' Velocity Averaging Slices in  ', ixyz_char  
			print'(a,f10.5)', ' With each averaging slice of width:', & 
						globaldomain(velocity_outflag)/gnbins(velocity_outflag)
		case(4)
			print'(3(a,i8),a)', ' Velocity 3D bins recorded every:', &
					tplot,' x ',Nvel_ave,' = ',tplot*Nvel_ave,' iterations'
			print'(a,3i8)', ' Domain split into Velocity Averaging Bins in x,y and z:', gnbins
			print'(a,3f10.5)', ' Each of size:', & 
			globaldomain(1)/gnbins(1), globaldomain(2)/gnbins(2),globaldomain(3)/gnbins(3)
		case(5)
			print'(3(a,i8),a)', ' Velocity 3D Cylindrical Polar bins recorded every:', &
					tplot,' x ',Nvel_ave,' = ',tplot*Nvel_ave,' iterations'
			print'(a,3i8)', ' Domain split into Velocity Averaging Bins in r,theta and z:', gcpol_bins
		case default
			call error_abort("Invalid Velocity output flag in input file")
		end select

		select case(Temperature_outflag)
		case(0)
			print*, 'Temperature record off'
			print*, ''
			print*, ''
		case(1:3)
			if (Temperature_outflag .eq. 1) then
				ixyz_char = 'x'
			elseif (Temperature_outflag .eq. 2) then 
				ixyz_char = 'y'
			else 
				ixyz_char = 'z'
			endif

			print'(3(a,i8),a)', ' Temperature slice recorded every:', & 
					tplot,' x ',Nvel_ave,' = ',tplot*Nvel_ave,' iterations'
			print'(a,i8,2a)', ' Domain split into',gnbins(Temperature_outflag) , & 
					  ' Temperature Averaging Slices in  ', ixyz_char  
			print'(a,f10.5)', ' With each averaging slice of width:', & 
						globaldomain(Temperature_outflag)/gnbins(Temperature_outflag)
		case(4)
			print'(3(a,i8),a)', ' Temperature 3D bins recorded every:', &
					tplot,' x ',Nvel_ave,' = ',tplot*NTemp_ave,' iterations'
			print'(a,3i8)', ' Domain split into Temperature Averaging Bins in x,y and z:', gnbins
			if (peculiar_flag .eq. 1) then
				print'(a,3f10.5,a)', ' Each of size:', & 
				globaldomain(1)/gnbins(1), globaldomain(2)/gnbins(2),globaldomain(3)/gnbins(3), &
				' and peculiar mometum is used '
			else
				print'(a,3f10.5)', ' Each of size:', & 
				globaldomain(1)/gnbins(1), globaldomain(2)/gnbins(2),globaldomain(3)/gnbins(3)
			endif
		case default
			call error_abort("Invalid Velocity output flag in input file")
		end select

		select case(pressure_outflag)
		case(0)
			print*, 'Pressure tensor off'
			print*, ''
			print*, ''
		case(1)
			print'(3(a,i8),a)', ' Pressure tensor Virial recorded every', & 
					tplot,' x ',Nstress_ave,' = ',tplot*Nstress_ave,' iterations'
			print*, 'Single Value for Whole Domain'
			print*, ''
		case(2)
			if (split_kin_config .eq. 0) then
				print'(3(a,i8),a)', ' Pressure tensor Volume Averaged recorded every', & 
						tplot,' x ',Nstress_ave,' = ',tplot*Nstress_ave,' iterations'
				print'(a,3i8)', ' Domain split into Pressure Volume Averaging Bins in x,y and z:', gnbins
				print'(a,3f10.5)', ' Each of size:', & 
				globaldomain(1)/gnbins(1), globaldomain(2)/gnbins(2),globaldomain(3)/gnbins(3)
			else
				print'(3(a,i8),a)', ' Seperate Kinetic/Configurational Pressure tensor Volume Averaged recorded every', & 
						tplot,' x ',Nstress_ave,' = ',tplot*Nstress_ave,' iterations'
				print'(a,3i8)', ' Domain split into Pressure Volume Averaging Bins in x,y and z:', gnbins
				print'(a,3f10.5)', ' Each of size:', & 
				globaldomain(1)/gnbins(1), globaldomain(2)/gnbins(2),globaldomain(3)/gnbins(3)
			endif
		case default
			call error_abort("Invalid Pressure tensor output flag in input file")
		end select

		select case(viscosity_outflag)
		case(0)
			print*, 'Viscosity calculation off'
			print*, ''
			print*, ''
		case(1)
			print'(a)', ' Viscosity calculated from Virial pressure using Green Kubo autocorrelation'
			print'(2(a,i8),a)', 'Taking', Nstress_ave, 'bins in time with a new bin written every', & 
				tplot*Nstress_ave,'iterations'
			print*, ' with each bin averaged', Nvisc_ave,' times'
		case default
			call error_abort("Invalid viscosity output flag in input file")
		end select

		select case(mflux_outflag)
		case(0)
			print*, 'Mass flux record off'
			print*, ''
			print*, ''
		case(1)
			if (cv_conserve .eq. 0) then
				print'(3(a,i8),a)', ' Mass flux over surface of 3D bins and snapshots recorded every:', &
						tplot,' x ',Nmflux_ave,' = ',tplot*Nmflux_ave,' iterations'
			else
				print'(a,i8,a)', ' Mass flux over surface of 3D bins and snapshots recorded every:', &
						Nmflux_ave,' iterations'
			endif
			print'(a,3i8)', ' Domain split into bins in x,y and z:', gnbins
			print'(a,3f10.5)', ' Each of size:', & 
			globaldomain(1)/gnbins(1), globaldomain(2)/gnbins(2),globaldomain(3)/gnbins(3)
		case default
			call error_abort("Invalid mass flux output flag in input file")
		end select

		select case(vflux_outflag)
		case(0)
			print*, 'Momentum Flux record off'
			print*, ''
			print*, ''
		case(1:3)
			if (vflux_outflag .eq. 1) then
				ixyz_char = 'x'
			elseif (vflux_outflag .eq. 2) then 
				ixyz_char = 'y'
			else 
				ixyz_char = 'z'
			endif
			print'(2a,3(a,i8),a)', ' Pressure tensor Method Of Planes in', ixyz_char, 'recorded every', &
					tplot,' x ',Nvflux_ave,' = ',tplot*Nvflux_ave,' iterations'
			print'(a,i8,a)', ' Domain split into ',nplanes,' Planes for Pressure Averaging' 
			print'(2(a,f10.5))', ' Seperated by distance:', planespacing, ' with first plane at ', & 
						 planes(1)
		case(4)
			if (cv_conserve .eq. 0) then
				print'(3(a,i8),a)', ' Momentum flux over surface of 3D bins and snapshots recorded every:', &
						tplot,' x ',Nvflux_ave,' = ',tplot*Nvflux_ave,' iterations'
			else
				print'(a,i8,a)', ' Momentum flux over surface of 3D bins and snapshots recorded every:', &
						Nvflux_ave,' iterations'
			endif
			print'(a,3i8)', ' Domain split into bins in x,y and z:', gnbins
			print'(a,3f10.5)', ' Each of size:', & 
			globaldomain(1)/gnbins(1), globaldomain(2)/gnbins(2),globaldomain(3)/gnbins(3)
		case default
			call error_abort("Invalid velocity flux output flag in input file")
		end select


		select case(eflux_outflag)
		case(0)
			print*, 'Energy Flux record off'
			print*, ''
			print*, ''
		case(1:3)
			call error_abort("Energy MOP not developed - change eflux output flag in input file")
			!if (eflux_outflag .eq. 1) then
			!	ixyz_char = 'x'
			!elseif (eflux_outflag .eq. 2) then 
			!	ixyz_char = 'y'
			!else 
			!	ixyz_char = 'z'
			!endif
			!print'(2a,3(a,i8),a)', ' Energy Method Of Planes in', ixyz_char, 'recorded every', &
			!		tplot,' x ',Neflux_ave,' = ',tplot*Neflux_ave,' iterations'
			!print'(a,i8,a)', ' Domain split into ',nplanes,' Planes for Energy Averaging' 
			!print'(2(a,f10.5))', ' Seperated by distance:', planespacing, ' with first plane at ', & 
			!			 planes(1)
		case(4)
			print'(3(a,i8),a)', ' Energy flux over surface of 3D bins and snapshots recorded every:', &
					tplot,' x ',Neflux_ave,' = ',tplot*Neflux_ave,' iterations'
			print'(a,3i8)', ' Domain split into bins in x,y and z:', gnbins
			print'(a,3f10.5)', ' Each of size:', & 
			globaldomain(1)/gnbins(1), globaldomain(2)/gnbins(2),globaldomain(3)/gnbins(3)
		case default
			call error_abort("Invalid energy flux output flag in input file")
		end select

		!print*, 'Bins per Processor:', nbins
		!print*, 'Number of Bins on outer Surface of each processor', nsurfacebins
		print*, 'Initial condition:'
		select case(potential_flag)
		case(0)
			select case(macro_outflag)
			case(1:2)
				print '(2a)', &
				'     iter;   simtime;      VSum;    V^2Sum;   Temp;', &
				'          KE;                 PE;                  TE;               P'
			case(3:4)
				print '(2a)', &
				'   iter;   simtime;    VSum;      T;', & 
				'      KE;      PE;      TE;       P'
			case default
			end select
		case(1)
			select case(macro_outflag)
			case(1:2)
				print '(2a)', &
				'     iter;   simtime;      VSum;     V^2Sum;   Temp;', & 
				'        KE;            PE;             TE;               P;   Rtcf;     R_g'
			case(3:4)
				print '(2a)', &
				'    iter; simtime;    VSum;   Temp;', & 
				'     KE;     PE;     TE;      P;  Rtcf;  R_g'
			case default
			end select
		case default
			call error_abort("Invalid potential flag in input file")
		end select

       select case(integration_algorithm)
       case(leap_frog_verlet)
           call print_macroscopic_properties(initialstep)
       case(velocity_verlet)
           call print_macroscopic_properties(initialstep)
       end select

		print*, "Simulation:"
!		print*, "   =========================================================="
!		print*,	"   |    _____ _                 _       _   _               |"
!		print*,	"   |   /  ___(_)               | |     | | (_)              |"
!		print*,	"   |   \ `--. _ _ __ ___  _   _| | __ _| |_ _  ___  _ __    |"
!		print*,	"   |    `--. \ | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_ \   |"
!		print*,	"   |   /\__/ / | | | | | | |_| | | (_| | |_| | (_) | | | |  |"
!		print*,	"   |   \____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|  |"
!		print*, "   |                                                        |"
!		print*, "   =========================================================="
!		select case(potential_flag)
!		case(0)
!			select case(macro_outflag)
!			case(1:2)
!				print '(2a)', &
!				'     iter;   simtime;      VSum;    V^2Sum;   Temp;', &
!				'          KE;                 PE;                  TE;               P'
!			case(3:4)
!				print '(2a)', &
!				'   iter;   simtime;    VSum;      T;', & 
!				'      KE;      PE;      TE;       P'
!			case default
!			end select
!		case(1)
!			select case(macro_outflag)
!			case(1:2)
!				print '(2a)', &
!				'     iter;   simtime;      VSum;    V^2Sum;   Temp;', & 
!				'        KE;            PE;             TE;               P;   Rtcf;      R_g'
!			case(3:4)
!				print '(2a)', &
!				'    iter; simtime;   VSum;   Temp;', & 
!				'     KE;     PE;     TE;      P;  Rtcf;   R_g'
!			case default
!			end select
!		case default
!			call error_abort("Invalid potential flag in input file")
!		end select

	
		!Obtain and record radial distributions
		select case (rdf_outflag)
		case(1)
			call evaluate_properties_rdf
			call rdf_io
		case(2)
			call evaluate_properties_rdf3d
			call rdf3d_io
		case default
		end select

		select case (ssf_outflag)
		case(1)
			call evaluate_properties_ssf
			call ssf_io
		case default
		end select

		call simulation_header

	endif

end subroutine setup_initial_record

!---------------------------------------------------------------------------------
! Write Simulation Parameter File contain all data required to completely recreate
! simulation and to be used for post processing
!
subroutine simulation_header
	use module_parallel_io
	use calculated_properties_MD
	use concentric_cylinders
	use librarymod, only : get_new_fileunit,get_version_number
	implicit none

	integer				:: fileunit, i
	Character(8)  		:: the_date
	Character(10)  		:: the_time, version_no

	call date_and_time(the_date, the_time)
	version_no = get_version_number()
	fileunit = get_new_fileunit()

	open(fileunit,file=trim(prefix_dir)//'results/simulation_header')

	write(fileunit,*) 'Simulation run on Date;  sim_date ;', the_date
	write(fileunit,*) 'Simulation start time ;  sim_start_time ;', the_time
	write(fileunit,*) 'Subversion revision number ;  svn_version_number ;', version_no
	write(fileunit,*) 'Number of Dimensions ;  nd ;', nd
	write(fileunit,*) 'Number of Particles ;  globalnp ;', globalnp
	write(fileunit,*) 'Time Step - delta t ;   delta_t ;',  delta_t
	write(fileunit,*) 'Total number of steps ; Nsteps;',  Nsteps - initialstep
	write(fileunit,*) 'Integration algorithm ; integration_algorithm;', integration_algorithm
	select case(potential_flag)
	case(0)
		write(fileunit,*) 'Potential flag ; potential_flag;', potential_flag
	case(1)
		write(fileunit,*) 'Potential flag ; potential_flag;', potential_flag
		write(fileunit,*) 'Number of LJ beads per FENE chain ; nmonomers;', nmonomers
		write(fileunit,*) 'Number of FENE chains in domain ; nchains;', nchains
		write(fileunit,*) 'FENE bond maximum elongation ; R_0;', R_0
		write(fileunit,*) 'FENE spring stiffness ; k_c;', k_c
	end select	
	write(fileunit,*) 'Starting step of simulation ;  initialstep ;', initialstep
	write(fileunit,*) 'Generate output file every steps ;   tplot ;',  tplot
	write(fileunit,*) 'Density ; density ;',density
	write(fileunit,*) 'Initial Temperature ;   inputtemperature ;',  inputtemperature
	write(fileunit,*) 'Cut off distance ;  rcutoff ;', rcutoff
	write(fileunit,*) 'Neighbour List Delta r ;  delta_rneighbr ;', delta_rneighbr
	write(fileunit,*) 'Initial FCC unit size in x ;  initialunitsize(1) ;', initialunitsize(1)
	write(fileunit,*) 'Initial FCC unit size in y ;  initialunitsize(2) ;', initialunitsize(2)
	write(fileunit,*) 'Initial FCC unit size in z ;  initialunitsize(3) ;', initialunitsize(3)
	write(fileunit,*) 'Domain in x ;  globaldomain(1)  ;', globaldomain(1) 
	write(fileunit,*) 'Domain in y ;  globaldomain(2)  ;', globaldomain(2) 
	write(fileunit,*) 'Domain in z ;  globaldomain(3)  ;', globaldomain(3) 
	write(fileunit,*) 'Domain volume ;  volume ;', volume
	write(fileunit,*) 'Periodic Boundary Conditions in x ;  periodic(1) ;', periodic(1)
	write(fileunit,*) 'Periodic Boundary Conditions in y ;  periodic(2) ;', periodic(2)
	write(fileunit,*) 'Periodic Boundary Conditions in z ;  periodic(3) ;', periodic(3)
	write(fileunit,*) 'Dist frm bot Fixed Mol in x; fixdistbot(1);', fixdistbottom(1)
	write(fileunit,*) 'Dist frm bot Fixed Mol in y; fixdistbot(2);', fixdistbottom(2)
	write(fileunit,*) 'Dist frm bot Fixed Mol in z; fixdistbot(3);', fixdistbottom(3)
	write(fileunit,*) 'Dist frm top Fixed Mol in x; fixdisttop(1);', fixdisttop(1)
	write(fileunit,*) 'Dist frm top Fixed Mol in y; fixdisttop(2);', fixdisttop(2)
	write(fileunit,*) 'Dist frm top Fixed Mol in z; fixdisttop(3);', fixdisttop(3)
	write(fileunit,*) 'Dist frm bot Tethered Mol in x; tethdistbot(1);', tethereddistbottom(1)
	write(fileunit,*) 'Dist frm bot Tethered Mol in y; tethdistbot(2);', tethereddistbottom(2)
	write(fileunit,*) 'Dist frm bot Tethered Mol in z; tethdistbot(3);', tethereddistbottom(3)
	write(fileunit,*) 'Dist frm top Tethered Mol in x; tethdisttop(1);', tethereddisttop(1)
	write(fileunit,*) 'Dist frm top Tethered Mol in y; tethdisttop(2);', tethereddisttop(2)
	write(fileunit,*) 'Dist frm top Tethered Mol in z; tethdisttop(3);', tethereddisttop(3)
	write(fileunit,*) 'Dist frm bot Sliding Mol in x; slidedistbot(1);', slidedistbottom(1)
	write(fileunit,*) 'Dist frm bot Sliding Mol in y; slidedistbot(2);', slidedistbottom(2)
	write(fileunit,*) 'Dist frm bot Sliding Mol in z; slidedistbot(3);', slidedistbottom(3)
	write(fileunit,*) 'Dist frm top Sliding Mol in x; slidedisttop(1);', slidedisttop(1)
	write(fileunit,*) 'Dist frm top Sliding Mol in y; slidedisttop(2);', slidedisttop(2)
	write(fileunit,*) 'Dist frm top Sliding Mol in z; slidedisttop(3);', slidedisttop(3)
	write(fileunit,*) 'Sliding velocity of wall in x; wallslidev(1);', wallslidev(1)
	write(fileunit,*) 'Sliding velocity of wall in y; wallslidev(2);', wallslidev(2)
	write(fileunit,*) 'Sliding velocity of wall in z; wallslidev(3);', wallslidev(3)
	write(fileunit,*) 'Dist frm bot NH Thermstat Mol in x; thermstatbot(1);', thermstatbottom(1)
	write(fileunit,*) 'Dist frm bot NH Thermstat Mol in y; thermstatbot(2);', thermstatbottom(2)
	write(fileunit,*) 'Dist frm bot NH Thermstat Mol in z; thermstatbot(3);', thermstatbottom(3)
	write(fileunit,*) 'Dist frm top NH Thermstat Mol in x; thermstattop(1);', thermstattop(1)
	write(fileunit,*) 'Dist frm top NH Thermstat Mol in y; thermstattop(2);', thermstattop(2)
	write(fileunit,*) 'Dist frm top NH Thermstat Mol in z; thermstattop(3);', thermstattop(3)
	write(fileunit,*) 'Computational cells in x ;  globalncells(1) ;',  ncells(1)*npx
	write(fileunit,*) 'Computational cells in y ;  globalncells(2)  ;', ncells(2)*npy 
	write(fileunit,*) 'Computational cells in z ;  globalncells(3)  ;', ncells(3)*npz 
	write(fileunit,*) 'Of size in x ;  cellsidelength(1) ;', cellsidelength(1)
	write(fileunit,*) 'Of size in y ;  cellsidelength(2) ;', cellsidelength(2)
	write(fileunit,*) 'Of size in z ;  cellsidelength(3) ;', cellsidelength(3)
	write(fileunit,*) 'Number of processors in x ;  npx ;', npx
	write(fileunit,*) 'Number of processors in y ;  npy ;', npy
	write(fileunit,*) 'Number of processors in z ;  npz ;', npz
	write(fileunit,*) 'Cells per Processor in x ;  nicellxl ;', nicellxl
	write(fileunit,*) 'Cells per Processor in y ;  nicellyl ;', nicellyl
	write(fileunit,*) 'Cells per Processor in z ;  nicellzl ;', nicellzl
	write(fileunit,*) 'Cells per Processor including Halos in x ;  ncellxl ;', ncellxl
	write(fileunit,*) 'Cells per Processor including Halos in y ;  ncellyl ;', ncellyl
	write(fileunit,*) 'Cells per Processor including Halos in z ;  ncellzl ;', ncellzl
	write(fileunit,*) '1st Random seed ;  seed_1 ;', seed(1)
	write(fileunit,*) '2nd Random seed ;  seed_2 ;', seed(2)
	write(fileunit,*)  'VMD flag ;  vmd_outflag ;', vmd_outflag
	do i=1,Nvmd_intervals
		if (i .lt. 10) then
			write(fileunit,'(a,i1,a,i1,a,i8)')  'VMD interval ',i,' start ;  vmd_start_0',i,' ;', vmd_intervals(1,i)
			write(fileunit,'(a,i1,a,i1,a,i8)')  'VMD interval ',i,' end   ;    vmd_end_0',i,' ;', vmd_intervals(2,i)
		else
			write(fileunit,'(a,i2,a,i2,a,i8)')  'VMD interval ',i,' start ;  vmd_start_',i,' ;', vmd_intervals(1,i)
			write(fileunit,'(a,i2,a,i2,a,i8)')  'VMD interval ',i,' end   ;    vmd_end_',i,' ;', vmd_intervals(2,i)
		endif
	enddo
	write(fileunit,*)  'macro flag ;  macro_outflag	 ;', macro_outflag
	write(fileunit,*)  'mass flag ;  mass_outflag ;', mass_outflag	
	write(fileunit,*)  'velocity flag ;  velocity_outflag ;', velocity_outflag
	write(fileunit,*)  'temperature flag ;  temperature_outflag ;', temperature_outflag
	write(fileunit,*)  'Peculiar flag ;  peculiar_flag ;', peculiar_flag
	write(fileunit,*)  'Pressure flag ;  pressure_outflag ;', pressure_outflag
	write(fileunit,*)  'viscosity flag ;  viscosity_outflag ;', viscosity_outflag
	write(fileunit,*)  'cv conservation flag ;  cv_conserve ;', cv_conserve
	write(fileunit,*)  'mass flux flag ;  mflux_outflag ;', mflux_outflag
	write(fileunit,*)  'velocity flux flag ;  vflux_outflag ;', vflux_outflag
	write(fileunit,*)  'energy flux flag ;  eflux_outflag ;', eflux_outflag
	write(fileunit,*)  'mass average steps ;  Nmass_ave ;', Nmass_ave
	write(fileunit,*)  'velocity average steps ;  Nvel_ave ;', Nvel_ave
	write(fileunit,*)  'Temperature average steps ;  NTemp_ave ;', NTemp_ave
	write(fileunit,*)  'pressure average steps ;  Nstress_ave ;', Nstress_ave
	write(fileunit,*)  'viscosity average samples ;  Nvisc_ave ;', Nvisc_ave
	write(fileunit,*)  'mass flux average steps ;  Nmflux_ave ;', Nmflux_ave
	write(fileunit,*)  'velocity flux average steps ;  Nvflux_ave ;', Nvflux_ave
	write(fileunit,*)  'energy flux average steps ;  Neflux_ave ;', Neflux_ave
	if (velocity_outflag .eq. 5 .or. mass_outflag .eq. 5 ) then
		write(fileunit,*)  'Velocity/stress Averaging Bins in r ;      gnbins(1) ;', gcpol_bins(1)
		write(fileunit,*)  'Velocity/stress Averaging Bins in theta ;  gnbins(2) ;', gcpol_bins(2)
		write(fileunit,*)  'Velocity/stress Averaging Bins in z ;      gnbins(3) ;', gcpol_bins(3)
		write(fileunit,*)  'Outer radius of outer cylinder ;   r_oo   ;', r_oo
		write(fileunit,*)  'Inner radius of outer cylinder ;   r_io   ;', r_io
		write(fileunit,*)  'Outer radius of inner cylinder ;   r_oi   ;', r_oi
		write(fileunit,*)  'Inner radius of inner cylinder ;   r_ii   ;', r_ii
	else
		write(fileunit,*)  'Velocity/stress Averaging Bins in x ;  gnbins(1) ;', gnbins(1)
		write(fileunit,*)  'Velocity/stress Averaging Bins in y ;  gnbins(2) ;', gnbins(2)
		write(fileunit,*)  'Velocity/stress Averaging Bins in z ;  gnbins(3) ;', gnbins(3)
		write(fileunit,*)  'Of size in x ;  binsize(1)  ;', globaldomain(1)/gnbins(1) 
		write(fileunit,*)  'Of size in y ;  binsize(2)  ;', globaldomain(2)/gnbins(2) 
		write(fileunit,*)  'Of size in z ;  binsize(3)  ;', globaldomain(3)/gnbins(3) 
		write(fileunit,*)  'Bins per Processor in x ;  nbins(1) ;', nbins(1)
		write(fileunit,*)  'Bins per Processor in y ;  nbins(2) ;', nbins(2)
		write(fileunit,*)  'Bins per Processor in z ;  nbins(3) ;', nbins(3)
		write(fileunit,*)  'Number of Bins on outer Surface of each processor ;  nsurfacebins ;', nsurfacebins
		write(fileunit,*)  'Number of Bins in halo of each processor ;  nhalobins ;', nhalobins
	end if
	if (vflux_outflag .gt.0 .and. vflux_outflag  .lt. 3) then
		write(fileunit,*)  'Domain split into Planes for Pressure Averaging ; nplanes  ;',nplanes 
		write(fileunit,*)  'Separated by distance ;  planespacing  ;', planespacing 
		write(fileunit,*)  'with first plane at ;  planes ;', planes(1)
	endif
	write(fileunit,*)  'Leapfrog or Velocity-Verlet ; integration_algorithm ;', integration_algorithm
	write(fileunit,*)  'Force calculation list methodd ; force_list ;', force_list
	!write(fileunit,*)  'Ensemble; ensemble; ', ensemble		!MATLAB input functions can't deal with words...
	write(fileunit,*)	'Shear direction ; le_sd;', le_sd
	write(fileunit,*)	'Shear plane ; le_sp;', le_sp
	write(fileunit,*)	'Shear remaining plane ; le_rp;', le_rp
	write(fileunit,*)	'Shear rate; le_sr;', le_sr
	write(fileunit,*)	'Shear velocity; le_sv;', le_sv
	write(fileunit,*)	'Shear iter0 ; le_i0;', le_i0
	write(fileunit,*)  'g(r) nbins ; rdf_nbins;', rdf_nbins
	write(fileunit,*)  'g(r) rmax  ; rdf_rmax;', rdf_rmax
	write(fileunit,*)  'S(k) nmax  ; ssf_nmax;', ssf_nmax
	write(fileunit,*)  'S(k) axis1 ; ssf_ax1;', ssf_ax1
	write(fileunit,*)  'S(k) axis2 ; ssf_ax2;', ssf_ax2

	close(fileunit,status='keep')

end subroutine simulation_header

!----------------------------------------------------------------------------------
!Calculate Initial kinetic and potential energy as well as temperature and pressure

subroutine initial_macroscopic_properties
use module_initial_record
use interfaces
implicit none

	integer          :: n, ixyz
	double precision :: vel

	vsum  = 0.d0      ! Reset all sums
	v2sum = 0.d0      ! Reset all sums

	!Calculate forces to obtain initial potential energies and virial
	call simulation_compute_forces
	
	select case(potential_flag)
	case(0)
		potenergysum = sum(potenergymol(1:np))
		call globalSum(potenergysum)
	case(1)
		potenergysum_LJ = sum(potenergymol_LJ(1:np))
		potenergysum_FENE = sum(potenergymol_FENE(1:np))
		potenergysum = sum(potenergymol_LJ(1:np) + potenergymol_FENE(1:np))
		call globalSum(potenergysum_LJ)
		call globalSum(potenergysum_FENE)
 		call globalSum(potenergysum)
	case default
		call error_abort("Unrecognised potential flag in initial_macroscopic_properties")
	end select

	virial = sum(virialmol(1:np))
    call globalSum(virial)  

	select case (integration_algorithm)
	case(leap_frog_verlet)
		do ixyz = 1, nd   ! Loop over all dimensions
		do n = 1, np    ! Loop over all particles
			!Velocity component must be shifted back half a timestep to determine 
			!velocity of interest - required due to use of the leapfrog method
			vel = v(ixyz,n) + 0.5d0*a(ixyz,n)*delta_t
			vsum = vsum + vel         !Add up all molecules' velocity components
			v2sum = v2sum + vel**2.d0 !Add up all molecules' velocity squared components  
		enddo
		enddo
		call globalSum(vsum)
		call globalSum(v2sum)
	case(velocity_verlet)
  		do ixyz = 1, nd   ! Loop over all dimensions
 		do n = 1, np    ! Loop over all particles  
			vsum  = vsum  + v(ixyz,n)
			v2sum = v2sum +  v(ixyz,n) * v(ixyz,n)
 		enddo
		enddo
		call globalSum(vsum)
		call globalSum(v2sum)
	end select

	!Obtain global sums for all parameters
	kinenergy   = (0.5d0 * v2sum) / real(globalnp,kind(0.d0))
	potenergy   = potenergysum /(2.d0*real(globalnp,kind(0.d0))) !N.B. extra 1/2 as all interactions calculated
	if (potential_flag.eq.1) then
		potenergy_LJ= potenergysum_LJ/(2.d0*real(globalnp,kind(0.d0)))
		potenergy_FENE= potenergysum_FENE/(2.d0*real(globalnp,kind(0.d0)))
	end if
	totenergy   = kinenergy + potenergy
	temperature = v2sum / real(nd*globalnp,kind(0.d0))
	if (any(periodic.gt.1)) temperature = get_temperature_PUT()
	pressure    = (density/(globalnp*nd))*(v2sum+virial/2) !N.B. virial/2 as all interactions calculated
	
	!WARNING ABOUT THERMOSTATTED RESTARTS! WILL CONSERVE TEMPERATURE OF LAST ITERATION
	!For velocity rescaling thermostat
	initialenergy = (potenergysum+v2sum)/(np)

end subroutine initial_macroscopic_properties

!----------------------------------------------------------------------------------
!Add up initial control volume mass densities
subroutine initial_control_volume
use module_initial_record
implicit none

	!Obtain and record velocity and mass
	if (velocity_outflag .ne. 0) then

		!Set mass record frequency to same as velocity
		Nmass_ave = Nvel_ave
		!Call first record		
		call velocity_averaging(velocity_outflag)
	endif

	!Obtain and record mass only
	if (velocity_outflag .eq. 0 .and. mass_outflag .ne. 0) then
		!Call first record
		call mass_averaging(mass_outflag)
	endif

	if (temperature_outflag .ne. 0) then
		!Call first record
		call temperature_averaging(temperature_outflag)
	endif

	if (mflux_outflag .ne. 0) call mass_snapshot
	if (vflux_outflag .eq. 4) call momentum_snapshot

	!Open pressure tensor and viscosity record file 
	!if (pressure_outflag .ne. 0) call stress_header
	!call control_volume_header

end subroutine initial_control_volume

!----------------------------------------------------------------------------------
! TODO -- COMMENT
subroutine build_psf
	use module_initial_record
	use polymer_info_MD
	implicit none

	integer :: i,j,n,item,molno,sc
	integer :: NTITLE, NATOM, NBONDS
	integer	:: write_items
	integer :: group, bin_expo
	integer :: plot_mod
	integer, allocatable, dimension(:,:) :: bonds	
	integer, allocatable, dimension(:)	 ::	res_ID
	integer, allocatable, dimension(:)	 ::	glob_sc
	integer, allocatable, dimension(:,:) ::	glob_bf
	character(len=4), allocatable, dimension(:) :: seg_name, res_name, atom_name, atom_type
	real, allocatable, dimension(:)	:: charge, mass

	if (irank.eq.iroot) then
		print*, 'Generating polymer topology file polymer_topol.psf...'
	end if
	
	NTITLE = 1												! How many 'REMARKS' lines you want
	NATOM  = globalnp										! Determine total number of atoms
	NBONDS = 0
	do n=1,np

		do sc=1,nmonomers
            j=0
			group = ceiling(real(sc)/real(intbits))
			bin_expo = mod(sc,intbits)-1
			if (bin_expo.eq.-1) bin_expo = intbits - 1
            if ( btest(monomer(n)%bin_bflag(group),bin_expo) ) j = 1
            NBONDS = NBONDS + j
		end do

	end do
	call globalSumInt(NBONDS)
	NBONDS = int(NBONDS/2)

	allocate(seg_name(NATOM))								! Determine segment names for each atom
	allocate(res_ID(NATOM))									! Determine molecule ID for each atom
	allocate(glob_sc(NATOM))								! Determine molecule ID for each atom
	allocate(glob_bf(4,NATOM))                              ! Determine molecule ID for each atom
	allocate(res_name(NATOM))								! Determine name for each molecule
	allocate(atom_name(NATOM))								! Determine name for each atom
	allocate(atom_type(NATOM))								! Determine type for each atom
	allocate(charge(NATOM))									! Determine charge for each atom
	allocate(mass(NATOM))									! Determine mass of each atom

	res_ID(:) = 0
	glob_sc(:) = 0
	glob_bf(:,:) = 0
	do n=1,np
		molno            = monomer(n)%glob_no
		res_ID(molno)    = monomer(n)%chainID 
		glob_sc(molno)   = monomer(n)%subchainID
		glob_bf(:,molno) = monomer(n)%bin_bflag(:)
	end do

	call globalSumIntVect(res_ID,globalnp)
	call globalSumIntVect(glob_sc,globalnp)
	call globalSumIntTwoDim(glob_bf,4,globalnp)

	do n=1,globalnp
		select case (res_ID(n))
		case(0)
			seg_name(n)  = 'N'
			res_name(n)  = 'SOL'
			atom_name(n) = 'N'
			atom_type(n) = 'N'
			mass(n)      = 1.00794
		case(1:)
			seg_name(n)  = 'C'
			res_name(n)  = 'POL'
			atom_name(n) = 'C'
			atom_type(n) = 'C'
			mass(n)      = 1.00794
		case default
		end select
		charge(n)    = 0.00000
	end do

	if (irank.eq.iroot) then
		
		open(unit=1, file=trim(prefix_dir)//"results/polymer_topol.psf", status='replace', form='formatted')

		! Header
		write(1,'(a3)') 'PSF'
		write(1,'(a1)')
		write(1,'(i8,a)') NTITLE, ' !NTITLE'
		write(1,'(a9,a)') "REMARKS ","FENE polymer protein structure file, written by MD_dCSE"
			
		! Atoms
		write(1,'(a1)')
		write(1,'(i8,a)') NATOM, ' !NATOM'

		do i=1,NATOM

			write(1,'(i8,3a,i4,6a,2f10.5,i1)') i,' ',seg_name(i),' ',res_ID(i),'&
					 & ',res_name(i),' ',atom_name(i),' ',atom_type(i),charge(i),mass(i),0 

		end do

		! Bonds
		write(1,'(a1)')
		write(1,'(i8,a)') NBONDS, ' !NBONDS'
	
		write_items = 4*ceiling(NBONDS/4.0)!todo change for branched polymers
		allocate(bonds(write_items,2))
		bonds(:,:)=0

		! Get ready to print progress bar
		!Open unit 6 (stdout) with fortran carriage control 
		open (unit=6, carriagecontrol='fortran')  
		plot_mod = max(1,globalnp/100)
		write(*,'(a)') ' Gathering PSF bond data: ' 

		item=1                                                   !Initialise bond item number
		do i=1,globalnp                                          !Loop through all molecules

			do j=1,i-1                                           !Avoid double counting
				if (res_ID(i).eq.res_ID(j)) then                 !If same global chainID
					
					!If j is bonded to i then add pair to items	
					do n=1,nmonomers
						group    = ceiling(real(n)/real(intbits))
						bin_expo = mod(n,intbits)-1
						if (bin_expo.eq.-1) bin_expo = intbits - 1

                        if(btest(glob_bf(group,i),bin_expo) .and. glob_sc(j) .eq. n) then
								bonds(item,1) = i
								bonds(item,2) = j
								item=item+1
						end if

					end do
	
				end if
			end do

			if (mod(i,plot_mod).eq.0) call progress(100*i/globalnp)

		end do

		!Write all bonds to .psf file
		do i=1,write_items,4
			write(1,'(8i8)')bonds(i,  1), bonds(i,  2),&
							bonds(i+1,1), bonds(i+1,2),&
							bonds(i+2,1), bonds(i+2,2),&
							bonds(i+3,1), bonds(i+3,2)
		end do

		! Angles
		! Dihedrals
		! Impropers
		! Donors
		! Acceptors
		! NNB
		! NGRP

		close(1, status='keep')	

		deallocate(bonds)

	end if
	
	deallocate(seg_name)								! Determine segment names for each atom
	deallocate(res_ID)									! Determine molecule ID for each atom
	deallocate(res_name)								! Determine name for each molecule
	deallocate(atom_name)								! Determine name for each atom
	deallocate(atom_type)								! Determine type for each atom
	deallocate(charge)									! Determine charge for each atom
	deallocate(mass)									! Determine mass of each atom
	if(allocated(bonds)) deallocate(bonds)
	deallocate(glob_sc)
	deallocate(glob_bf)

	if (irank.eq.iroot) print*, 'Finished generating polymer topology file.'

end subroutine build_psf
