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
implicit none

        integer i
        logical file_exist
	character		:: ixyz_char
	Character(8)  		:: the_date
	Character(10)  		:: the_time
        Character(10),parameter :: file_names(17) = (/"mslice    ", "mbins     ", "msnap     ",&
                                                      "vslice    ", "vbins     ", "vsnap     ",&
                                                      "pvirial   ", "pVA       ", "visc      ",&
                                                      "mflux     ", "vflux     ", "pplane    ",&
                                                      "psurface  ", "esnap     ", "eflux     ",&
                                                      "eplane    ", "esurface  "/) 
       

	if (irank .eq. iroot) then
		!Delete existing files
		open (unit=5, file=trim(prefix_dir)//'results/mslice')
		close(5,status='delete')
		open (unit=5, file=trim(prefix_dir)//'results/mbins')
		close(5,status='delete')
		open (unit=5, file=trim(prefix_dir)//'results/msnap')
		close(5,status='delete')
		open (unit=6, file=trim(prefix_dir)//'results/vslice')
		close(6,status='delete')
		open (unit=6, file=trim(prefix_dir)//'results/vbins')
		close(6,status='delete')
		open (unit=6, file=trim(prefix_dir)//'results/vsnap')
		close(6,status='delete')
		open (unit=6, file=trim(prefix_dir)//'results/vslice')
		close(6,status='delete')
		open (unit=6, file=trim(prefix_dir)//'results/vbins')
		close(6,status='delete')
		open (unit=7, file=trim(prefix_dir)//'results/pvirial')
		close(7,status='delete')
		open (unit=7, file=trim(prefix_dir)//'results/pVA')
		close(7,status='delete')
		open (unit=7, file=trim(prefix_dir)//'results/pVA_k')
		close(7,status='delete')
		open (unit=7, file=trim(prefix_dir)//'results/pVA_c')
		close(7,status='delete')
		open (unit=7, file=trim(prefix_dir)//'results/visc')
		close(7,status='delete')
		open (unit=8, file=trim(prefix_dir)//'results/mflux')
		close(8,status='delete')
		open (unit=9, file=trim(prefix_dir)//'results/vflux')
		close(9,status='delete')
		open (unit=9, file=trim(prefix_dir)//'results/pplane')
		close(9,status='delete')
		open (unit=9, file=trim(prefix_dir)//'results/psurface')
		close(9,status='delete')
		open (unit=10, file=trim(prefix_dir)//'results/esnap')
		close(10,status='delete')
		open (unit=10, file=trim(prefix_dir)//'results/eflux')
		close(10,status='delete')
		open (unit=10, file=trim(prefix_dir)//'results/eplane')
		close(10,status='delete')
		open (unit=10, file=trim(prefix_dir)//'results/esurface')
		close(10,status='delete')
	endif

        do i=1,size(file_names)
            inquire(file=trim(prefix_dir)//'results/'//file_names(i),exist=file_exist)
            if(file_exist) then
               open (unit=23, file=trim(prefix_dir)//'results/'//file_names(i))
               close(23,status='delete')
            endif
        enddo

!!$	if (irank .eq. iroot) then
!!$		!Delete existing files
!!$		open (unit=5, file=trim(prefix_dir)//'results/mslice')
!!$		close(5,status='delete')
!!$		open (unit=5, file=trim(prefix_dir)//'results/mbins')
!!$		close(5,status='delete')
!!$		open (unit=5, file=trim(prefix_dir)//'results/msnap')
!!$		close(5,status='delete')
!!$		open (unit=6, file=trim(prefix_dir)//'results/vslice')
!!$		close(6,status='delete')
!!$		open (unit=6, file=trim(prefix_dir)//'results/vbins')
!!$		close(6,status='delete')
!!$		open (unit=6, file=trim(prefix_dir)//'results/vsnap')
!!$		close(6,status='delete')
!!$		open (unit=6, file=trim(prefix_dir)//'results/vslice')
!!$		close(6,status='delete')
!!$		open (unit=6, file=trim(prefix_dir)//'results/vbins')
!!$		close(6,status='delete')
!!$		open (unit=7, file=trim(prefix_dir)//'results/pvirial')
!!$		close(7,status='delete')
!!$		open (unit=7, file=trim(prefix_dir)//'results/pVA')
!!$		close(7,status='delete')
!!$		open (unit=7, file=trim(prefix_dir)//'results/visc')
!!$		close(7,status='delete')
!!$		open (unit=8, file=trim(prefix_dir)//'results/mflux')
!!$		close(8,status='delete')
!!$		open (unit=9, file=trim(prefix_dir)//'results/vflux')
!!$		close(9,status='delete')
!!$		open (unit=9, file=trim(prefix_dir)//'results/pplane')
!!$		close(9,status='delete')
!!$		open (unit=9, file=trim(prefix_dir)//'results/psurface')
!!$		close(9,status='delete')
!!$		open (unit=10, file=trim(prefix_dir)//'results/esnap')
!!$		close(10,status='delete')
!!$		open (unit=10, file=trim(prefix_dir)//'results/eflux')
!!$		close(10,status='delete')
!!$		open (unit=10, file=trim(prefix_dir)//'results/eplane')
!!$		close(10,status='delete')
!!$		open (unit=10, file=trim(prefix_dir)//'results/esurface')
!!$		close(10,status='delete')
!!$	endif

	!Evaluate system properties on all processes
	call initial_macroscopic_properties
	
	!Calculate Control Volume starting state
	call initial_control_volume
	if (irank .eq. iroot) then

		call date_and_time(the_date, the_time)

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
			print*, 'Interatomic potential: LJ only'
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
		print'(a,f19.15,a,f10.5)', ' Cut off distance:  ', rcutoff, &
				    '  Neighbour List Delta r:  ', delta_rneighbr
		print*, 'Initial unit size (FCC unit) in x,y and z:'
		print*,	initialunitsize(1), initialunitsize(2), initialunitsize(3)
		print*, 'Domain in x,y and z: '
		print*,  globaldomain(1), globaldomain(2), globaldomain(3)
		print*, 'Domain volume: ', volume
		print'(a,3i8)', ' Periodic Boundary Conditions in x,y and z:', periodic
		print'(a,3f10.5)', ' Distance from bottom of Fixed Molecules in x,y and z:', 	fixdistbottom
		print'(a,3f10.5)', ' Distance from top of Fixed Molecules in x,y and z:', 	fixdisttop
		print'(a,3f10.5)', ' Distance from bottom of Tethered Molecules in x,y and z:', tethereddistbottom
		print'(a,3f10.5)', ' Distance from top of Tethered Molecules in x,y and z:', 	tethereddisttop
		print'(a,3f10.5)', ' Distance from bottom of Sliding Molecules in x,y and z:', 	slidedistbottom
		print'(a,3f10.5)', ' Distance from top of Sliding Molecules in x,y and z:', 	slidedisttop
		print'(a,3f10.5)', ' Velocity of Sliding Molecules in x,y and z:', 	wallslidev
		print'(a,3f10.5)', ' Distance from bottom of NH Themostatted Molecules in x,y and z:', 	thermstatbottom
		print'(a,3f10.5)', ' Distance from top of NH Themostatted Molecules in x,y and z:', 	thermstattop
		print*, 'Molecular Reynolds number = ', (density * maxval(slidev(1:np,1)) * domain(2))/1.d0
		print*, '==================== Computational Parameters ========================='
		print'(a,3i8)', ' Domain split into computational cells in x,y and z:', & 
					 ncells(1)*npx, ncells(2)*npy, ncells(3)*npz
		print'(a,3f10.5)', ' Each of size:', cellsidelength
		print*, 'Number of processors in x,y and z: ', npx, npy, npz
		print*, 'Cells per Processor:', nicellxl, nicellyl, nicellzl
		print*, 'Cells per Processor including Halos: ', ncellxl, ncellyl, ncellzl
		print*, 'Random seed used for initial velocity generation:', seed
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
		case default
			call error_abort("Invalid VMD output flag in input file")
		end select

		select case(macro_outflag)
		case(0)
			print*, 'No Macroscopic Properties printed to screen'
		case(1)
			print*, 'Macroscopic properties printed to screen every:', tplot, 'iterations'
		case(2)
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
			print'(a,i8,2a)', ' Domain split into',globalnbins(mass_outflag) ,' mass Averaging Slices in  ', ixyz_char  
			print'(a,f10.5)', ' With each averaging slice of width:', & 
						globaldomain(mass_outflag)/globalnbins(mass_outflag)
		case(4)
			print'(3(a,i8),a)', ' mass 3D bins recorded every:', &
					tplot,' x ',Nvel_ave,' = ',tplot*Nvel_ave,' iterations'
			print'(a,3i8)', ' Domain split into mass Averaging Bins in x,y and z:', globalnbins
			print'(a,3f10.5)', ' Each of size:', & 
			globaldomain(1)/globalnbins(1), globaldomain(2)/globalnbins(2),globaldomain(3)/globalnbins(3)
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
			print'(a,i8,2a)', ' Domain split into',globalnbins(velocity_outflag) , & 
					  ' Velocity Averaging Slices in  ', ixyz_char  
			print'(a,f10.5)', ' With each averaging slice of width:', & 
						globaldomain(velocity_outflag)/globalnbins(velocity_outflag)
		case(4)
			print'(3(a,i8),a)', ' Velocity 3D bins recorded every:', &
					tplot,' x ',Nvel_ave,' = ',tplot*Nvel_ave,' iterations'
			print'(a,3i8)', ' Domain split into Velocity Averaging Bins in x,y and z:', globalnbins
			print'(a,3f10.5)', ' Each of size:', & 
			globaldomain(1)/globalnbins(1), globaldomain(2)/globalnbins(2),globaldomain(3)/globalnbins(3)
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
				print'(a,3i8)', ' Domain split into Pressure Volume Averaging Bins in x,y and z:', globalnbins
				print'(a,3f10.5)', ' Each of size:', & 
				globaldomain(1)/globalnbins(1), globaldomain(2)/globalnbins(2),globaldomain(3)/globalnbins(3)
			else
				print'(3(a,i8),a)', ' Seperate Kinetic/Configurational Pressure tensor Volume Averaged recorded every', & 
						tplot,' x ',Nstress_ave,' = ',tplot*Nstress_ave,' iterations'
				print'(a,3i8)', ' Domain split into Pressure Volume Averaging Bins in x,y and z:', globalnbins
				print'(a,3f10.5)', ' Each of size:', & 
				globaldomain(1)/globalnbins(1), globaldomain(2)/globalnbins(2),globaldomain(3)/globalnbins(3)
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
			print'(a,3i8)', ' Domain split into bins in x,y and z:', globalnbins
			print'(a,3f10.5)', ' Each of size:', & 
			globaldomain(1)/globalnbins(1), globaldomain(2)/globalnbins(2),globaldomain(3)/globalnbins(3)
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
			print'(a,3i8)', ' Domain split into bins in x,y and z:', globalnbins
			print'(a,3f10.5)', ' Each of size:', & 
			globaldomain(1)/globalnbins(1), globaldomain(2)/globalnbins(2),globaldomain(3)/globalnbins(3)
		case default
			call error_abort("Invalid velocity flux output flag in input file")
		end select

		!print*, 'Bins per Processor:', nbins
		!print*, 'Number of Bins on outer Surface of each processor', nsurfacebins
		print*, '======================================================================='

		select case(potential_flag)
		case(0)
			print '(2a)', &
			'Iteration; 	   VSum;        V^2Sum;        Temp;', & 
			'          KE;                 PE;                  TE;          Pressure;'
			!Print initial conditions for simulations at iteration 0
			print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f10.4)', &
			initialstep,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy,';',totenergy,';',pressure
		case(1)
			print '(2a)', &
			'Iteration; 	   VSum;        V^2Sum;        Temp;', & 
			'       KE;     PE (LJ);  PE (FENE); PE (Tot);    TE;       Pressure;'
			!Print initial conditions for simulations at iteration 0
			print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.4)', &
			initialstep,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy_LJ,';',potenergy_FENE,';',potenergy,';',totenergy,';',pressure
		case default
			call error_abort("Invalid potential flag in input file")
		end select

		call simulation_header

	endif

	!Initialise etevtcf calculation if etevtcf_iter0 = 0
	if (potential_flag.eq.1) then
		if (etevtcf_outflag.ne.0) call etevtcf_calculate_parallel
		if (etevtcf_outflag.eq.2) call etev_io
		
		if (r_gyration_outflag.ne.0) call r_gyration_calculate
		if (r_gyration_outflag.eq.2) call r_gyration_io
	end if

end subroutine setup_initial_record

!----------------------------------------------------------------------------------
!Calculate Initial kinetic and potential energy as well as temperature and pressure

subroutine initial_macroscopic_properties
use module_initial_record
implicit none

	integer          :: n, ixyz
	double precision :: vel

	vsum  = 0.d0      ! Reset all sums
	v2sum = 0.d0      ! Reset all sums

	!Calculate forces to obtain initial potential energies and virial
	call simulation_compute_forces
	
	do n = 1, np    ! Loop over all particles

		select case(potential_flag)
		case(0)
			potenergysum	= potenergysum + potenergymol(n)
		case(1)
			potenergysum_LJ = potenergysum_LJ + potenergymol_LJ(n)
			potenergysum_FENE = potenergysum_FENE + potenergymol_FENE(n)
			potenergysum = potenergysum + potenergymol_LJ(n) + potenergymol_FENE(n)
		end select

		virial = virial + virialmol(n)

		select case (integration_algorithm)
		case(leap_frog_verlet)
			do ixyz = 1, nd   ! Loop over all dimensions
				!Velocity component must be shifted back half a timestep to determine 
				!velocity of interest - required due to use of the leapfrog method
				vel = v(n,ixyz) + 0.5d0*a(n,ixyz)*delta_t
				vsum = vsum + vel      !Add up all molecules' velocity components
				v2sum = v2sum + vel**2 !Add up all molecules' velocity squared components  
			enddo
		case(velocity_verlet)
			vsum = vsum + sum(v(n,:))
			v2sum = v2sum + dot_product(v(n,:),v(n,:))
		end select

	enddo

	!Obtain global sums for all parameters
	call globalSum(vsum)
	call globalSum(v2sum)
	call globalSum(virial)
	call globalSum(potenergysum)
	call globalSum(potenergysum_LJ)
	call globalSum(potenergysum_FENE)

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
	!kinenergy   = (0.5d0 * v2sum) / globalnp 
	!potenergy   = potenergysum /(2*globalnp) !N.B. extra 1/2 as all interactions calculated
	!if (potential_flag.eq.1) then
	!	potenergy_LJ = potenergysum_LJ/(2.d0*real(globalnp,kind(0.d0)))
	!	potenergy_FENE = potenergysum_FENE/(2.d0*real(globalnp,kind(0.d0)))
	!end if
	!totenergy   = kinenergy + potenergy
	!temperature = v2sum / (nd * globalnp)
	!if (thermstat_flag.eq.2) temperature = get_temperature_PUT()
	!pressure    = (density/(globalnp*nd))*(v2sum+virial/2) !N.B. virial/2 as all interactions calculated
	
	!WARNING ABOUT THERMOSTATTED RESTARTS! WILL CONSERVE TEMPERATURE OF LAST ITERATION
	!For velocity rescaling thermostat
	initialenergy = (potenergysum+v2sum)/(np)

end subroutine initial_macroscopic_properties

!----------------------------------------------------------------------------------
!Add up initial control volume mass densities

subroutine initial_control_volume
use module_initial_record
implicit none

	integer							:: n
	integer, dimension(3)			:: ibin
	double precision, dimension(3)	:: mbinsize

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


