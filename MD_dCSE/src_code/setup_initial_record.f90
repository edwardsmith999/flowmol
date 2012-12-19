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
	implicit none

	integer					:: i
	logical 				:: file_exist
	character				:: ixyz_char
	character(8)			:: the_date
	character(10)			:: the_time
	character(23),parameter :: file_names(23) = &
								(/ "mslice      ", "mbins       ", "msnap   ",&
								   "vslice      ", "vbins       ", "vsnap   ",&
								   "pvirial     ", "pVA         ", "pVA_k   ",& 
								   "pVA_c       ", "visc        ", "mflux   ",& 
								   "vflux       ", "pplane      ", "psurface",&
								   "esnap       ", "eflux       ", "eplane  ",&
								   "esurface    ", "viscometrics", "rdf     ",&
	                               "rdf3d       ", "ssf         "            /) 
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
		end select
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
		print'(a,3f10.5)', ' Distance from bottom of Fixed Molecules in x,y and z:', 	fixdistbottom
		print'(a,3f10.5)', ' Distance from top of Fixed Molecules in x,y and z:', 	fixdisttop
		print'(a,3f10.5)', ' Distance from bottom of Tethered Molecules in x,y and z:', tethereddistbottom
		print'(a,3f10.5)', ' Distance from top of Tethered Molecules in x,y and z:', 	tethereddisttop
		print'(a,3f10.5)', ' Distance from bottom of Sliding Molecules in x,y and z:', 	slidedistbottom
		print'(a,3f10.5)', ' Distance from top of Sliding Molecules in x,y and z:', 	slidedisttop
		print'(a,3f10.5)', ' Velocity of Sliding Molecules in x,y and z:', 	wallslidev
		print'(a,3f10.5)', ' Distance from bottom of NH Themostatted Molecules in x,y and z:', 	thermstatbottom
		print'(a,3f10.5)', ' Distance from top of NH Themostatted Molecules in x,y and z:', 	thermstattop
		print*, 'Molecular Reynolds number = ', (density * maxval(slidev(1,1:np)) * domain(2))/1.d0
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
			print'(3(a,i8),a)', ' Mass flux over surface of 3D bins and snapshots recorded every:', &
					tplot,' x ',Nmflux_ave,' = ',tplot*Nmflux_ave,' iterations'
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
			print'(3(a,i8),a)', ' Momentum flux over surface of 3D bins and snapshots recorded every:', &
					tplot,' x ',Nvflux_ave,' = ',tplot*Nvflux_ave,' iterations'
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
	write(3,*)  'cv conservation flag ;  cv_conserve ;', cv_conserve
	write(3,*)  'mass flux flag ;  mflux_outflag ;', mflux_outflag
	write(3,*)  'velocity flux flag ;  vflux_outflag ;', vflux_outflag
	write(3,*)  'energy flux flag ;  eflux_outflag ;', eflux_outflag
	write(3,*)  'mass average steps ;  Nmass_ave ;', Nmass_ave
	write(3,*)  'velocity average steps ;  Nvel_ave ;', Nvel_ave
	write(3,*)  'Temperature average steps ;  NTemp_ave ;', NTemp_ave
	write(3,*)  'pressure average steps ;  Nstress_ave ;', Nstress_ave
	write(3,*)  'viscosity average samples ;  Nvisc_ave ;', Nvisc_ave
	write(3,*)  'mass flux average steps ;  Nmflux_ave ;', Nmflux_ave
	write(3,*)  'velocity flux average steps ;  Nvflux_ave ;', Nvflux_ave
	write(3,*)  'energy flux average steps ;  Neflux_ave ;', Neflux_ave
	write(3,*)  'Velocity/stress Averaging Bins in x ;  gnbins(1) ;', gnbins(1)
	write(3,*)  'Velocity/stress Averaging Bins in y ;  gnbins(2) ;', gnbins(2)
	write(3,*)  'Velocity/stress Averaging Bins in z ;  gnbins(3) ;', gnbins(3)
	write(3,*)  'Of size in x ;  binsize(1)  ;', globaldomain(1)/gnbins(1) 
	write(3,*)  'Of size in y ;  binsize(2)  ;', globaldomain(2)/gnbins(2) 
	write(3,*)  'Of size in z ;  binsize(3)  ;', globaldomain(3)/gnbins(3) 
	write(3,*)  'Bins per Processor in x ;  nbins(1) ;', nbins(1)
	write(3,*)  'Bins per Processor in y ;  nbins(2) ;', nbins(2)
	write(3,*)  'Bins per Processor in z ;  nbins(3) ;', nbins(3)
	write(3,*)  'Number of Bins on outer Surface of each processor ;  nsurfacebins ;', nsurfacebins
	write(3,*)  'Number of Bins in halo of each processor ;  nhalobins ;', nhalobins
	if (vflux_outflag .gt.0 .and. vflux_outflag  .lt. 3) then
		write(3,*)  'Domain split into Planes for Pressure Averaging ; nplanes  ;',nplanes 
		write(3,*)  'Separated by distance ;  planespacing  ;', planespacing 
		write(3,*)  'with first plane at ;  planes ;', planes(1)
	endif
	write(3,*)  'Leapfrog or Velocity-Verlet ; integration_algorithm ;', integration_algorithm
	write(3,*)  'Force calculation list methodd ; force_list ;', force_list
	!write(3,*)  'Ensemble; ensemble; ', ensemble		!MATLAB input functions can't deal with words...
	write(3,*)	'Shear direction ; le_sd;', le_sd
	write(3,*)	'Shear plane ; le_sp;', le_sp
	write(3,*)	'Shear remaining plane ; le_rp;', le_rp
	write(3,*)	'Shear rate; le_sr;', le_sr
	write(3,*)	'Shear velocity; le_sv;', le_sv
	write(3,*)	'Shear iter0 ; le_i0;', le_i0
	write(3,*)  'g(r) nbins ; rdf_nbins;', rdf_nbins
	write(3,*)  'g(r) rmax  ; rdf_rmax;', rdf_rmax
	write(3,*)  'S(k) nmax  ; ssf_nmax;', ssf_nmax
	write(3,*)  'S(k) axis1 ; ssf_ax1;', ssf_ax1
	write(3,*)  'S(k) axis2 ; ssf_ax2;', ssf_ax2

	close(3,status='keep')

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

	if (mflux_outflag .ne. 0) call mass_snapshot
	if (vflux_outflag .eq. 4) call momentum_snapshot

	!Open pressure tensor and viscosity record file 
	!if (pressure_outflag .ne. 0) call stress_header
	!call control_volume_header

end subroutine initial_control_volume
