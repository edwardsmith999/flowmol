!==========================================================================
!                     CALCULATE AND RECORD PROPERTIES
!==========================================================================
! Evaluate properties of interest from the simulation. 
! This includes macroscopic properties including sum of velocity/velocity^2, 
! kinetic energy, potential energy, temperature and pressure.
! Also calculated is the diffusion; velocity distributions and Boltzmanns
! H function; the radial distribution function; the velocity magnitude by cells 
! to be used for contours; plane based uni-directional velocity slices; 
! the stress tensor calculated from a virial (bulk or local for a homogenous fluid),
! Volume averaged (local assignment based on proportion of force/momentum is in a 
! given cell) or Method Of Planes (force/momentum flux over a plane).
! Many of these outputs are interpreted using MATLAB scripts included
!
! simulation_record  			   	Top level function choosing logging based on inputs
! evaluate_macroscopic_properties_parallel 	Macroscopic properties
! evaluate_properties_vdistribution		Maxwell Boltzmann distribution
! evaluate_properties_radialdist		Radial distribution
! evaluate_properties_diffusion			Diffusion

! mass_averaging					slice/CV
! cumulative_mass					slice/CV
! velocity_averaging				slice/CV
! cumulative_velocity				slice/CV
! pressure_averaging				virial/CV
! cumulative_pressure				virial/CV
! simulation_compute_kinetic_VA		
! simulation_compute_kinetic_VA_cells
! mass_flux_averaging				CV_surface
! cumulative_mass_flux				CV_surface
! mass_snapshot						CV
! momentum_flux_averaging			MOP_plane/CV_surface
! cumulative_momentum_flux			MOP_plane/CV_surface
! momentum_snapshot					CV
! pressure_tensor_forces			virial
! pressure_tensor_forces_VA			CV
! control_volume_forces				CV_surface
! control_volume_stresses			CV_surface
! pressure_tensor_forces_MOP		MOP_plane
! evaluate_properties_cellradialdist

!===================================================================================
! 			FULL DOMAIN VIRIAL STRESS CALCULATION
! Calculate pressure_tensor for entire domain- average over all xy, yz and xz 
! components performed when calculating the interactions of particles
!===================================================================================
!===================================================================================
! 			VOLUME AVERAGE STRESS TENSOR
! Stress calculated for each cell by assinging to cells based on location of 
! molecules (kinetic) and fraction of pairwise interaction in cell (potential)
! First used in J. F. Lutsko, J. Appl. Phys 64(3), 1152 (1988) and corrected in
! J. Cormier, J.M. Rickman and T. J. Delph (2001) Journal of Applied Physics,
! Volume 89, No. 1 "Stress calculations in atomistic simulations of perfect 
! and imperfect solids" 
!===================================================================================
!===================================================================================
!	PRESSURE TENSOR CALCULATED USING METHOD OF PLANES (MOP)
! See B.D.Todd, D.J.Evans and P.J.Daivis (1995) Phys Review E, Vol 52,2
! Method first presented by Lutzko (1988) J. AppL Phys, Vol. 64, No.3,1 1154	!!!NO LUTSKO FIRST VA USING FOURIER TRANSFORMS!!!
!===================================================================================

module module_record

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use calculated_properties_MD
	use polymer_info_MD
        use librarymod

	double precision :: vel

end module module_record
!==========================================================================


subroutine simulation_record
	use module_record
	implicit none

	!Only record every tplot iterations
	if (mod(iter,tplot) .ne. 0) return

	!Parallel output for molecular positions
	if (vmd_outflag .eq. 1) call parallel_io_vmd
	if (vmd_outflag .eq. 2) call parallel_io_vmd_sl
	if (vmd_outflag .eq. 3) then
		call parallel_io_vmd
		call parallel_io_vmd_halo
	end if

	!call evaluate_macroscopic_properties
	!Obtain each processe's subdomain's macroscopic 
	!properties; gather on root process and record
	if (macro_outflag .eq. 1) call evaluate_macroscopic_properties_parallel
	if (macro_outflag .eq. 2) then
		call evaluate_macroscopic_properties_parallel
		call macroscopic_properties_record
	end if

	!Obtain and record velocity distributions
	!call evaluate_properties_vdistribution

	!Obtain and record radial distributions
	!call evaluate_properties_radialdist

	!Obtain and record velocity and mass
	if (velocity_outflag .ne. 0) call velocity_averaging(velocity_outflag)

	!Obtain and record mass only
	if (velocity_outflag .eq. 0 .and. &
		mass_outflag .ne. 0) call mass_averaging(mass_outflag)


	!Obtain and record molecular diffusion
	!call evaluate_properties_diffusion

	!Calculate pressure tensor
	if (pressure_outflag .ne. 0) call pressure_averaging(pressure_outflag)

	!Check if backup file should be saved based on elapsed iterations
	!and average size of current system per processor
	if (mod((iter*globalnp/nproc),1000000000).eq.0) then
		call messenger_syncall
		call parallel_io_final_state
	endif
	
	if (potential_flag.eq.1) then
		if (etevtcf_outflag.ne.0) call etevtcf_calculate
		if (etevtcf_outflag.eq.2) call etevtcf_io
	end if

	call update_simulation_progress_file

end subroutine simulation_record

!==========================================================================
!Calculate kinetic and potential energy as well as temperature and pressure

subroutine evaluate_macroscopic_properties_parallel
	use module_record
	implicit none

	integer			        :: n, ixyz

	vsum  = 0.d0      ! Reset all sums
	v2sum = 0.d0      ! Reset all sums

	do n = 1, np    ! Loop over all particles

		if (potential_flag.eq.0) then
			potenergysum	= potenergysum + potenergymol(n)
		else if (potential_flag.eq.1) then
			potenergysum_LJ = potenergysum_LJ + potenergymol_LJ(n)
			potenergysum_FENE = potenergysum_FENE + potenergymol_FENE(n)
			potenergysum = potenergysum + potenergymol_LJ(n) + potenergymol_FENE(n)
		end if

		virial = virial + virialmol(n)

		do ixyz = 1, nd   ! Loop over all dimensions
			!Velocity component must be shifted back half a timestep to determine 
			!velocity of interest - required due to use of the leapfrog method
			vel = v(n,ixyz) + 0.5d0*a(n,ixyz)*delta_t
			vsum = vsum + vel      !Add up all molecules' velocity components
			v2sum = v2sum + vel**2 !Add up all molecules' velocity squared components  
		enddo

	enddo

	!Obtain global sums for all parameters
	call globalSum(vsum)
	call globalSum(v2sum)
	call globalSum(virial)
	call globalSum(potenergysum)

	!Root processes prints results
	if (irank .eq. iroot) then
		

		kinenergy   = (0.5d0 * v2sum) / real(globalnp,kind(0.d0))
		potenergy   = potenergysum /(2.d0*real(globalnp,kind(0.d0))) !N.B. extra 1/2 as all interactions calculated
		if (potential_flag.eq.1) then
			potenergy_LJ= potenergysum_LJ/(2.d0*real(globalnp,kind(0.d0)))
			potenergy_FENE= potenergysum_FENE/(2.d0*real(globalnp,kind(0.d0)))
		end if
		totenergy   = kinenergy + potenergy
		temperature = v2sum / real(nd*globalnp,kind(0.d0))
		if (thermstat_flag.eq.2) temperature = get_temperature_PUT()
		pressure    = (density/(globalnp*nd))*(v2sum+virial/2) !N.B. virial/2 as all interactions calculated
	
		if (potential_flag.eq.0) then	
			print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f10.4)', &
			iter,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy,';',totenergy,';',pressure
		else if (potential_flag.eq.1) then
			print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.4)', &
			iter,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy_LJ,';',potenergy_FENE,';',potenergy,';',totenergy,';',pressure
		end if
		!print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.4,a,f10.4)', &
		!iter,';',vsum,';', v2sum,';', temperature,';', &
		!kinenergy,';',potenergy,';',totenergy,';',(density/(globalnp*nd))*(v2sum),';',(density/(globalnp*nd))*virial/2

	endif

	!Broacast pressure to all processes for constraint force
	call globalbroadcast(pressure,1,iroot)

end subroutine evaluate_macroscopic_properties_parallel

!===================================================================================
!Molecules grouped into velocity ranges to give vel frequency distribution graph
!and calculate Boltzmann H function

subroutine evaluate_properties_vdistribution
	use module_record
	implicit none

	integer          :: n
	integer          :: cbin
	double precision :: Hfunction

	!Calculate matrix of velocity magnitudes
	do n = 1, np    ! Loop over all particles
		!Assign to bins using integer division
 		cbin = ceiling(vmagnitude(n)/binsize)   !Establish current bin
		if (cbin > nbins(1)) cbin = nbins(1) 		!Prevents out of range values
		if (cbin < 1 ) cbin = 1        		!Prevents out of range values
		vfd_bin(cbin) = vfd_bin(cbin)+1		!Add one to current bin
	enddo

	!Normalise bins to use for output and H function - so all bins add up to one
	!Total stored molecules must be equal to number of molecules (np) times the
	!number of times bins have been assigned (iter/tplot) with the 1.d0 to avoid integer
	!division
	normalisedvfd_bin=0
	do n=1,nbins(1)
		normalisedvfd_bin(n) = vfd_bin(n)/((np*iter)*1.d0/tplot) 
	enddo

	!Calculate Boltzmann H function using discrete defintion as in
	!Rapaport p37. N.B. velocity at middle or range is used for vn
	Hfunction = 0.d0 !Reset H function before re-calculating
	do n=1,nbins(1)
		if (normalisedvfd_bin(n) .ne. 0) then
			Hfunction=Hfunction+normalisedvfd_bin(n)*log(normalisedvfd_bin(n)/(((n-0.5d0)*binsize)**(nd-1)))
		endif
	enddo
	write(12,'(a,i5, a, f20.10)') 'Boltzmann H function at iteration ', iter , ' is ', Hfunction

	!Write values of bin to file to follow evolution of distribution function
	write(12,'(a)') 'Velocity frequency distribution'
	do n=1,nbins(1) 
		write(12,'(2(f10.5))') (n-0.5d0)*binsize, normalisedvfd_bin(n) 
	enddo

end subroutine evaluate_properties_vdistribution

!===================================================================================
!Calculate Radial distribution function (RDF) using all pairs

subroutine evaluate_properties_radialdist
	use module_record
	implicit none

	integer                         :: n, i, j, ixyz   !Define dummy index
	integer				:: currentshell
	double precision		:: rmagnitude, rd2
	double precision, dimension(nd) :: rij       !'nd' directional vector between particles i and j
	
	rd2 = rd**2    !Calculate rd2 once to reduce computation

	do i = 1,np-1			!Step through each particle i 
	do j = i+1,np        		!Step through half of particles j for each i
		rmagnitude=0		!Set rmagnitude to zero

		do ixyz=1,nd
			rij(ixyz) = r (i,ixyz) - r(j,ixyz)          !Evaluate distance between particle i and j
			!If more than half a domain betwen molecules, must be closer over periodic boundaries      
			if (abs(rij(ixyz)) > halfdomain(ixyz)) then  !N.B array operator-corresponding element compared 
				rij(ixyz) = rij(ixyz) - sign(domain(ixyz),rij(ixyz))    !Sign of rij applied to domain
			endif
			rmagnitude = rmagnitude+rij(ixyz)*rij(ixyz) !Square of vector calculated
		enddo

		!Assign to shell using inter division
		if (rmagnitude < rd2) then
			rmagnitude = sqrt(rmagnitude)
			currentshell = ceiling(rmagnitude/delta_r)  !Assign to a current shell
			shell(currentshell) = shell(currentshell)+1 !Add one to that current shell
		endif
	enddo
	enddo

	!Discrete RDF using binned particles
	!In 3D g(rn)=(V*hn)/(2*pi*np^2*rn^2*delta_r*measurementno.) where
	!rn=(n-0.5)*delta_r and measurementno. = iter/tplot
	RDF = 0		!Set RFD ro zero ready for next loop
	if (nd==2) then
		do n = 1, nshells
			RDF(n) = shell(n) * volume / (pi*(np**2.d0)*((n-0.5d0)*delta_r)*delta_r*(iter/tplot))
		enddo
	else
		do n = 1, nshells
			RDF(n) = shell(n) * volume / (2.d0*pi*(np**2.d0)*(((n-0.5d0)*delta_r)**2.d0)*delta_r*(iter/tplot))
		enddo
	endif

	write(12,'(a)') 'Radial distribution function'
	do n = 1, nshells
		write(12,'(2(f10.5))') (n-0.5d0)*delta_r, RDF(n)
	enddo

end subroutine evaluate_properties_radialdist

!===================================================================================
!Diffusion function calculated

subroutine evaluate_properties_diffusion
	use module_record
	implicit none

	integer          :: n

	diffusion = 0
	do n=1,np
		diffusion(:)=diffusion(:)+(rtrue(n,:) - rinitial(n,:))**2
	enddo
	meandiffusion(:) = diffusion(:)/(2.d0*nd*np*(delta_t*iter/tplot))
	!print*, 'Instantanous diffusion', diffusion
	print*, 'time average of diffusion', meandiffusion

	!Using the auto-correlation function
	do n=1,np
		diffusion(:) = v(n,:) - 0.5 * a(n,:) * delta_t
	enddo

end subroutine evaluate_properties_diffusion

!===================================================================================
!Calculate end-to-end time correlation function of FENE chain
subroutine etevtcf_calculate 
use module_record
implicit none
	
	integer 			:: i,molL,molR
	double precision	:: etev_prod, etev_prod_sum
	double precision	:: etev2, etev2_sum
	double precision, dimension(nd) :: etev

	if (iter.eq.etevtcf_iter0) then						!Initialise end-to-end vectors at t_0
		do i=1,np														
			if (polyinfo_mol(i)%left.eq.0) then
				molL = i
				molR = i+(chain_length-1)										!Right end of chain
				etev_0(i,:) = r(molR,:) - r(molL,:)								!End-to-end vector
				etev_0(i,:) = etev_0(i,:) - domain(:)*anint(etev_0(i,:)/domain(:)) !Minimum image	
			end if
		end do
	end if

	etev_prod_sum	= 0.d0
	etev2_sum 		= 0.d0

	do i=1,np
		if (polyinfo_mol(i)%left.eq.0) then
			
			molL = i
			molR = i+(chain_length-1)

			etev(:) 		= r(molR,:) - r(molL,:)								!End-to-end vector
			etev(:) 		= etev(:) - domain(:)*anint(etev(:)/domain(:))		!Minimum image

			etev_prod		= dot_product(etev(:),etev_0(i,:))
			etev2			= dot_product(etev,etev)			
		
			etev_prod_sum	= etev_prod_sum + etev_prod
			etev2_sum		= etev2_sum + etev2		
		
		end if
	end do
	
	!etev_prod_mean = etev_prod_sum/dble(samplecount)
	!etev2_mean = etev2_sum/dble(samplecount)
	etevtcf = etev_prod_sum/etev2_sum											!Sample counts cancel
	if (etevtcf_outflag.eq.1) print*, 'ETEVTCF = ', etevtcf

end subroutine etevtcf_calculate

!===================================================================================
!		RECORD MASS AT LOCATION IN SPACE
! Either by binning or taking slices on a plane, the mass in the molecular
! system is recorded and output
!===================================================================================

!Calculate averaged velocity components of each bin or slice with 2D slice in 
!ixyz = 1,2,3 or in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine mass_averaging(ixyz)
	use module_record
	implicit none

	integer						:: ixyz, n,i,j,k
	integer, save					:: average_count=-1
	integer		,dimension(3)			:: ibin
	double precision,dimension(3)			:: ri, mbinsize

	average_count = average_count + 1
	call cumulative_mass(ixyz)
	if (average_count .eq. Nmass_ave) then
		average_count = 0
		!Determine bin size

		select case(ixyz)
		case(1:3)
			call mass_slice_io(ixyz)
			!Reset mass slice
			slice_mass = 0
		case(4)
			call mass_bin_io(volume_mass,'bins')
			!Reset mass slice
			volume_mass = 0
		case default
			stop "Error input for velocity averaging incorrect"
		end select

		!Collect mass for next step
		call cumulative_mass(ixyz)

	endif


end subroutine mass_averaging

!-----------------------------------------------------------------------------------
!Add mass to running total, with 2D slice in ixyz = 1,2,3 or
!in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine cumulative_mass(ixyz)
	use module_record
	use linked_list
	implicit none

	integer         		:: n, ixyz
	integer         		:: cbin
	integer,dimension(3)		:: ibin
	double precision		:: slicebinsize
	double precision,dimension(3) 	:: mbinsize 

	select case(ixyz)
	!Mass measurement is a number of 2D slices through the domain
	case(1:3)
		slicebinsize = domain(ixyz) / nbins(ixyz)
		do n = 1, np    ! Loop over all particles
			!Assign to bins using integer division
			cbin = ceiling((r(n,ixyz)+halfdomain(ixyz))/slicebinsize)!Establish current bin
			if (cbin > nbins(ixyz)) cbin = nbins(ixyz) 		 !Prevents out of range values
			if (cbin < 1 ) cbin = 1        				 !Prevents out of range values
			slice_mass(cbin)= slice_mass(cbin)+1      			 !Add one to current bin
		enddo
	!Mass measurement for 3D bins throughout the domain
	case(4)
		mbinsize(:) = domain(:) / nbins(:) 
		do n = 1,np
			!Add up current volume mass densities
			ibin(:) = ceiling((r(n,:)+halfdomain(:))/mbinsize(:)) + 1
			volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + 1
		enddo
	case default 
		stop "Mass Binning Error"
	end select
 
end subroutine cumulative_mass

!===================================================================================
!		RECORD VELOCITY AT LOCATION IN SPACE
! Either by binning or taking slices on a plane, the velocity field in the molecular
! system is recorded and output
!===================================================================================

!Calculate averaged velocity components of each bin or slice with 2D slice in 
!ixyz = 1,2,3 or in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine velocity_averaging(ixyz)
	use module_record
	use linked_list
	implicit none

	integer			:: ixyz
	integer, save		:: average_count=-1
	
	average_count = average_count + 1
	call cumulative_velocity(ixyz)
	if (average_count .eq. Nvel_ave) then
		average_count = 0

		select case(ixyz)
		case(1:3)
			call velocity_slice_io(ixyz)
			!Reset velocity slice
			slice_mass = 0
			slice_momentum  = 0.d0
		case(4)
			call velocity_bin_io(volume_mass,volume_momentum,'bins')
			!Reset velocity bins
			volume_mass = 0
			volume_momentum = 0.d0
		case default
			stop "Error input for velocity averaging incorrect"
		end select

		!Collect velocities for next step
		call cumulative_velocity(ixyz)

	endif

end subroutine velocity_averaging

!-----------------------------------------------------------------------------------
!Add velocities to running total, with 2D slice in ixyz = 1,2,3 or
!in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine cumulative_velocity(ixyz)
	use module_record
	use linked_list
	implicit none

	integer				:: n,m,i,j,k,ixyz
	integer         		:: cbin
	integer		,dimension(3)	:: ibin
	double precision		:: slicebinsize
	double precision,dimension(3) 	:: Vbinsize 

	select case(ixyz)
	!Velocity measurement is a number of 2D slices through the domain
	case(1:3)

		slicebinsize = domain(ixyz) / nbins(ixyz)

		do n = 1, np    ! Loop over all particles
			!Assign to bins using integer division
			cbin = ceiling((r(n,ixyz)+halfdomain(ixyz))/slicebinsize)!Establish current bin
			if (cbin > nbins(ixyz)) cbin = nbins(ixyz) 		 !Prevents out of range values
			if (cbin < 1 ) cbin = 1        				 !Prevents out of range values
			slice_mass(cbin)= slice_mass(cbin)+1      			 !Add one to current bin
			slice_momentum(cbin,:) = slice_momentum(cbin,:)+v(n,:) + slidev(n,:) 	 !Add streamwise velocity to current bin
		enddo

	!Velocity measurement for 3D bins throughout the domain
	case(4)
		!Determine bin size
		Vbinsize(:) = domain(:) / nbins(:)

		!Reset Control Volume momentum 
		do n = 1,np
			!Add up current volume mass and momentum densities
			ibin(:) = ceiling((r(n,:)+halfdomain(:))/Vbinsize(:)) + 1
			volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + 1
			volume_momentum(ibin(1),ibin(2),ibin(3),:) = volume_momentum(ibin(1),ibin(2),ibin(3),:) & 
										+ v(n,:) + slidev(n,:)
		enddo

	case default 
		stop "Velocity Binning Error"
	end select
	 
end subroutine cumulative_velocity

!===================================================================================
!		RECORD PRESSURE AT LOCATION IN SPACE
! Either by Volume Average style binning or virial expresion for the whole plane
!===================================================================================

subroutine pressure_averaging(ixyz)
	use module_record
	implicit none

	integer		:: ixyz
	integer, save	:: sample_count, average_count

	sample_count = sample_count + 1
	call cumulative_pressure(ixyz,sample_count)
	if (sample_count .eq. Nstress_ave) then
		sample_count = 0
		Pxyzero = Pxy		!Update Pxy(0) value

		select case(ixyz)
		case(1)
		!FULL DOMAIN VIRIAL STRESS CALCULATION
			call virial_stress_io
		case(2)
		!VA STRESS CALCULATION
			call VA_stress_io
		case default 
			stop "Average Pressure Binning Error"
		end select
		if(viscosity_outflag .eq. 1) then
			average_count = average_count+1
			if (average_count .eq. Nvisc_ave) then
				call viscosity_io
				average_count = 0
			endif
		endif
	endif

end subroutine pressure_averaging

!===================================================================================
!Add pressure_tensor to running total

subroutine cumulative_pressure(ixyz,sample_count)
	use module_record
	implicit none

	integer					:: sample_count,i,n,ixyz,jxyz,kxyz
	integer					:: imin, jmin, kmin, imax, jmax, kmax
	double precision, dimension(3)		:: velvect

	Pxy = 0.d0
	Pxymol = 0.d0

	!Factor of 2 as every interaction calculated
	rfmol = rfmol / 2.d0

	select case(ixyz)
	case(1)
		!FULL DOMAIN VIRIAL STRESS CALCULATION
		do n = 1, np    ! Calculate pressure tensor for all molecules

			!Calculate velocity at time t (v is at t+0.5delta_t due to use of verlet algorithm)
			velvect(:) = v(n,:) - 0.5d0 * a(n,:) * delta_t

			do jxyz = 1,3
			do kxyz = 1,3
				Pxymol(n,kxyz,jxyz) = Pxymol(n,kxyz,jxyz)    &
						      + velvect(kxyz)*velvect(jxyz) &
						      + rfmol(n,kxyz,jxyz)
			enddo
			enddo

			!Sum of molecules to obtain pressure tensor for entire domain
			Pxy(:,:) = Pxy(:,:) + Pxymol(n,:,:)

		enddo

		!Sum pressure tensor over all processors
		call globalSumVect(Pxy(:,1), nd)
		call globalSumVect(Pxy(:,2), nd)
		call globalSumVect(Pxy(:,3), nd)

		!Divide sum of stress by volume
		Pxy = Pxy / volume

		!Reset position force tensor before calculation
		rfmol = 0.d0
	case(2)
		!VOLUME AVERAGE STRESS CALCULATION
		!Calculate Position (x) force for configurational part of stress tensor
		call  simulation_compute_rfbins(1,nbins(1)+2,1,nbins(2)+2,1,nbins(3)+2)
		!Calculate velocity (x) velocity for kinetic part of stress tensor
		if (nbins(1) .eq. ncells(1) .and. & 
		    nbins(2) .eq. ncells(2) .and. & 
		    nbins(3) .eq. ncells(3)) then
			call  simulation_compute_kinetic_VA_cells(2,nbins(1)+1,2,nbins(2)+1,2,nbins(3)+1)
		else
			call  simulation_compute_kinetic_VA(2,nbins(1)+1,2,nbins(2)+1,2,nbins(3)+1)
		endif

		!Add results to cumulative total
		Pxybin(:,:,:,:,:) = Pxybin(:,:,:,:,:) + vvbin(:,:,:,:,:) 	& 
						      + rfbin(2:nbins(1)+1, 	& 
							      2:nbins(2)+1, 	& 
							      2:nbins(3)+1,:,:)/2.d0

		!NOTE: REBUILD AT (mod(iter+1,tplot) .eq. 0) WHEN RECORD AFTER FORCE CALCULATION
		!Reset bin force tensor before next calculation
	  	rfbin = 0.d0; vvbin = 0.d0
		Pxy(:,:) = 0.d0

		!Calculate Virial from sum of Volume Averaged Stresses
		do kxyz = 1,3
		do jxyz = 1,3
			Pxy(kxyz,jxyz) = sum(Pxybin(:,:,:,kxyz,jxyz))
		enddo
		enddo
		Pxy = Pxy / volume

	case default 
		stop "Cumulative Pressure Averaging Error"
	end select

	!Calculate correrlation between current value of Pxy and intial value of sample
	!Average 3 components of Pxycorrel to improve statistics
	if(viscosity_outflag .eq. 1) then
		if(sample_count .ne. 0) then
			Pxycorrel(sample_count) = Pxycorrel(sample_count) + Pxy(2,3)*Pxyzero(2,3)
			Pxycorrel(sample_count) = Pxycorrel(sample_count) + Pxy(1,3)*Pxyzero(1,3)
			Pxycorrel(sample_count) = Pxycorrel(sample_count) + Pxy(1,2)*Pxyzero(1,2)
		endif
	endif


end subroutine cumulative_pressure


!----------------------------------------------------------------------------------
!Compute kinetic part of stress tensor

subroutine simulation_compute_kinetic_VA(imin,imax,jmin,jmax,kmin,kmax)
	use module_record
	use physical_constants_MD
	implicit none

	integer					:: imin, jmin, kmin, imax, jmax, kmax
	integer         			:: n, ixyz,jxyz
	integer 				:: ibin, jbin, kbin
	double precision, dimension(3)		:: VAbinsize, velvect

	vvbin = 0.d0

	!Determine bin size
	VAbinsize(:) = domain(:) / nbins(:)

	! Add kinetic part of pressure tensor for all molecules
	do n = 1, np

		!Assign to bins using integer division
		ibin = ceiling((r(n,1)+halfdomain(1))/VAbinsize(1))	!Establish current bin
		if (ibin .gt. imax) cycle ! ibin = maxbin		!Prevents out of range values
		if (ibin .lt. imin) cycle ! ibin = minbin		!Prevents out of range values
		jbin = ceiling((r(n,2)+halfdomain(2))/VAbinsize(2)) 	!Establish current bin
		if (jbin .gt. jmax) cycle ! jbin = maxbin 		!Prevents out of range values
		if (jbin .lt. jmin) cycle ! jbin = minbin		!Prevents out of range values
		kbin = ceiling((r(n,3)+halfdomain(3))/VAbinsize(3)) 	!Establish current bin
		if (kbin .gt. kmax) cycle ! kbin = maxbin		!Prevents out of range values
		if (kbin .lt. kmin) cycle ! kbin = minbin		!Prevents out of range values

		!Calculate velocity at time t (v is at t+0.5delta_t due to verlet algorithm)
		velvect(:) = v(n,:) - 0.5 * a(n,:) * delta_t

		do ixyz = 1,3
		do jxyz = 1,3
			vvbin(ibin,jbin,kbin,ixyz,jxyz) = vvbin(ibin,jbin,kbin,ixyz,jxyz)	&
					       		  + velvect(ixyz) * velvect(jxyz)
		enddo
		enddo

	enddo

end subroutine simulation_compute_kinetic_VA

!----------------------------------------------------------------------------------
!Compute kinetic part of stress tensor ONLY IF BINSIZE = CELLSIZE

subroutine simulation_compute_kinetic_VA_cells(imin,imax,jmin,jmax,kmin,kmax)
	use module_record
	use physical_constants_MD
	use linked_list
	implicit none

	integer                         :: i, j, ixyz, jxyz   !Define dummy index
	integer 			:: ibin, jbin, kbin,icell, jcell, kcell
	integer                         :: cellnp, molnoi
	integer				:: imin, jmin, kmin, imax, jmax, kmax
	double precision, dimension(3)	:: velvect
	type(node), pointer 	        :: oldi, currenti

	vvbin = 0.d0

	! Add kinetic part of pressure tensor for all cells
	do kcell=kmin, kmax
	do jcell=jmin, jmax
	do icell=imin, imax

		!ASSUME Cell same size as bins
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list
		ibin = icell-1; jbin = jcell-1; kbin = kcell-1	

		do i = 1,cellnp                  !Step through each particle in list 
			molnoi = oldi%molno 	 !Number of molecule

			!Calculate velocity at time t (v is at t+0.5delta_t due to verlet algorithm)
			velvect(:) = v(molnoi,:) - 0.5 * a(molnoi,:) * delta_t

			do ixyz = 1,3
			do jxyz = 1,3
				vvbin(ibin,jbin,kbin,ixyz,jxyz) = vvbin(ibin,jbin,kbin,ixyz,jxyz) &
								  + velvect(ixyz) * velvect(jxyz)
			enddo
			enddo

			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo
	enddo
	enddo

	nullify(oldi)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required

end subroutine simulation_compute_kinetic_VA_cells

!===================================================================================
!		CONTROL VOLUME FORMULATION OF MASS AND STRESS
! Based on the control volume formulation of mass and momentum, it is possible to
! define change in properties in terms of surface fluxes over a volume in space
! Fluxes need to be called every timestep!!!
!===================================================================================

!!-----------------------------------------------------------------------------------
! Control Volume mass continuity
!-----------------------------------------------------------------------------------

subroutine mass_flux_averaging
	use module_record
	implicit none

	integer, save	:: sample_count

	!Only average if mass averaging turned on
	if (mflux_outflag .eq. 0) return

	call cumulative_mass_flux
	sample_count = sample_count + 1
	if (sample_count .eq. Nmflux_ave) then
		call mass_flux_io
		!call control_volume_mass_evolution
		sample_count = 0
		mass_flux = 0
		call mass_snapshot
	endif

end subroutine mass_flux_averaging

!===================================================================================
! Mass Flux over a surface of a bin

subroutine cumulative_mass_flux
	use module_record
	implicit none

	integer				:: ixyz, n
	integer		,dimension(3)	:: ibin1,ibin2,crossplane
	double precision,dimension(3)	:: ri, mbinsize

	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	do n = 1,np

		!Assign to bins before and after using integer division
		ibin1(:) = ceiling((r(n,:)+halfdomain(:))/mbinsize(:))+1
		ibin2(:) = ceiling((r(n,:)-delta_t*v(n,:)+halfdomain(:))/mbinsize(:))+1

		!Replace Signum function with this functions which gives a
		!check for plane crossing and the correct sign 
		crossplane(:) =  ibin1(:) - ibin2(:)

		if (sum(abs(crossplane(:))) .ne. 0) then

			!Find which direction the surface is crossed
			!For simplicity, if more than one surface has been crossed surface fluxes of intermediate cells
			!are not included. This assumption => more reasonable as Delta_t => 0
			!imaxloc = maxloc(abs(crossplane))
			ixyz = imaxloc(abs(crossplane)) 

			!Add mass flux to the new bin surface count and take from the old
			mass_flux(ibin1(1),ibin1(2),ibin1(3),ixyz+3*heaviside(dble(crossplane(ixyz)))) = & 
				mass_flux(ibin1(1),ibin1(2),ibin1(3),ixyz+3*heaviside(dble(crossplane(ixyz)))) & 
				+ abs(crossplane(ixyz))
			mass_flux(ibin2(1),ibin2(2),ibin2(3),ixyz+3*heaviside(-dble(crossplane(ixyz)))) = & 
				mass_flux(ibin2(1),ibin2(2),ibin2(3),ixyz+3*heaviside(-dble(crossplane(ixyz)))) &
				- abs(crossplane(ixyz))
		endif

	enddo

end subroutine cumulative_mass_flux

!===================================================================================
! Control Volume snapshot of the mass in a given bin

subroutine mass_snapshot
	use module_record
	implicit none

	integer						:: n
	integer		,dimension(3)			:: ibin
	integer		,dimension(:,:,:)  ,allocatable	:: volume_mass_temp
	double precision,dimension(3)			:: mbinsize

	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	!Allocate temporary array for mass and momentum in volume
	allocate(volume_mass_temp(nbins(1)+2,nbins(2)+2,nbins(3)+2))

	!Reset Control Volume momentum 
	volume_mass_temp = 0
	do n = 1,np
		!Add up current volume momentum densities
		ibin(:) = ceiling((r(n,:)+halfdomain(:))/mbinsize(:)) + 1
		volume_mass_temp(ibin(1),ibin(2),ibin(3)) = volume_mass_temp(ibin(1),ibin(2),ibin(3)) + 1
		!print*, n,r(n,:), ibin, sum(volume_mass_temp)
	enddo

	!Output Control Volume momentum change and fluxes
	call mass_bin_io(volume_mass_temp,'snap')

	deallocate(volume_mass_temp)

end subroutine mass_snapshot


!===================================================================================
! Control Volume Momentum continuity
!===================================================================================

subroutine momentum_flux_averaging(ixyz)
	use module_record
	implicit none

	integer			:: ixyz, i
	integer, save		:: sample_count
	double precision	:: binface

	if (vflux_outflag .eq. 0) return

	call cumulative_momentum_flux(ixyz)
	sample_count = sample_count + 1
	if (sample_count .eq. Nvflux_ave) then

		select case(ixyz)
		case(1:3)
		!MOP momentum flux and stresses
			call MOP_stress_io(ixyz)
			Pxy_plane = 0.d0
		case(4)
		!CV momentum flux and stress
			call momentum_flux_io
			call surface_stress_io
			momentum_flux = 0.d0
			Pxyface = 0.d0
			call momentum_snapshot
		case default 
			stop "Momentum flux and pressure averaging Error"
		end select

		sample_count = 0

	endif

end subroutine momentum_flux_averaging

!===================================================================================
! Momentum Flux over a surface of a bin including all intermediate bins

subroutine cumulative_momentum_flux(ixyz)
	use module_record
	implicit none

	integer				:: ixyz,jxyz,kxyz,i,j,k,n
	integer				:: planeno
	integer				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	!integer		,dimension(1)	:: imaxloc
	integer		,dimension(3)	:: ibin1,ibin2,cbin
	double precision		:: crosstime,crossplane,rplane,shift
	double precision,dimension(3)	:: mbinsize,velvect,crossface
	double precision,dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

	select case(ixyz)
	case(1:3)
		!MOP momentum flux
		!Shift by half difference between value rounded down and actual value
		!to ensure same distance to top and bottom plane from domain edge
		shift = 0.5d0*(domain(ixyz) - planespacing * (nplanes-1))

		!Add molecular velocities to the configuration stresses
		do n=1,np

			!Replace Signum function with this functions which gives a
			!check for plane crossing and the correct sign 
			crossplane = ceiling((r(n,ixyz)+halfdomain(ixyz)-shift)/planespacing) & 
				    -ceiling((r(n,ixyz)-delta_t*v(n,ixyz)+halfdomain(ixyz)-shift)/planespacing)

			if (crossplane .ne. 0) then

				!Obtain nearest plane number by integer division 
				!and retrieve location of plane from array
				planeno = ceiling((r(n,ixyz)+halfdomain(ixyz)-shift) 	& 
					  /planespacing)-heaviside(dble(crossplane))+1
				if (planeno .lt. 1) planeno = 1
				if (planeno .gt. nplanes) planeno = nplanes
				rplane = planes(planeno)

				!Calculate velocity at time of intersection
				!crosstime = (r(n,ixyz) - rplane)/v(n,ixyz)
				velvect(:) = v(n,:) !- a(n,:) * crosstime

				if (crosstime/delta_t .gt. 1.d0) stop "error in kinetic MOP"

				!Obtain stress for three components on y plane
				Pxy_plane(:,planeno) = Pxy_plane(:,planeno) + velvect(:)!*crossplane

			endif
		enddo

	case(4)
		!CV momentum flux
		!Determine bin size
		mbinsize(:) = domain(:) / nbins(:)

		do n = 1,np

			ri1(:) = r(n,:)			!Molecule i at time  t
			ri2(:) = r(n,:)-delta_t*v(n,:)	!Molecule i at time t-dt
			ri12   = ri1 - ri2		!Molecule i trajectory between t-dt and t
			where (ri12 .eq. 0.d0) ri12 = 0.000001d0

			!Assign to bins before and after using integer division
			ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:))+1
			ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:))+1

			!Replace Signum function with this functions which gives a
			!check for plane crossing and the correct sign 
			crossface(:) =  ibin1(:) - ibin2(:)

			if (sum(abs(crossface(:))) .ne. 0) then

				do i = ibin1(1),ibin2(1),sign(1,ibin2(1)-ibin1(1))
				do j = ibin1(2),ibin2(2),sign(1,ibin2(2)-ibin1(2))
				do k = ibin1(3),ibin2(3),sign(1,ibin2(3)-ibin1(3))

					cbin(1) = i; cbin(2) = j; cbin(3) = k

					bintop(:) = (cbin(:)-1)*mbinsize(:)-halfdomain(:)
					binbot(:) = (cbin(:)-2)*mbinsize(:)-halfdomain(:)

					!Calculate the plane intersect of trajectory with surfaces of the cube
					Pxt=(/ 			bintop(1), 		     & 
						ri1(2)+(ri12(2)/ri12(1))*(bintop(1)-ri1(1)), & 
						ri1(3)+(ri12(3)/ri12(1))*(bintop(1)-ri1(1))  	/)
					Pxb=(/ 			binbot(1), 		     & 
						ri1(2)+(ri12(2)/ri12(1))*(binbot(1)-ri1(1)), & 
						ri1(3)+(ri12(3)/ri12(1))*(binbot(1)-ri1(1))  	/)
					Pyt=(/	ri1(1)+(ri12(1)/ri12(2))*(bintop(2)-ri1(2)), & 
								bintop(2), 		     & 
						ri1(3)+(ri12(3)/ri12(2))*(bintop(2)-ri1(2))  	/)
					Pyb=(/	ri1(1)+(ri12(1)/ri12(2))*(binbot(2)-ri1(2)), &
								binbot(2), 		     & 
						ri1(3)+(ri12(3)/ri12(2))*(binbot(2)-ri1(2))  	/)
					Pzt=(/	ri1(1)+(ri12(1)/ri12(3))*(bintop(3)-ri1(3)), & 
						ri1(2)+(ri12(2)/ri12(3))*(bintop(3)-ri1(3)), &
								bintop(3) 			/)
					Pzb=(/	ri1(1)+(ri12(1)/ri12(3))*(binbot(3)-ri1(3)), &
						ri1(2)+(ri12(2)/ri12(3))*(binbot(3)-ri1(3)), & 
								binbot(3) 			/)

					onfacexb =0.5d0*(sign(1.d0,binbot(1) - ri2(1)) 	 & 
						       - sign(1.d0,binbot(1) - ri1(1)))* &
							(heaviside(bintop(2) - Pxb(2)) 	 &
						       - heaviside(binbot(2) - Pxb(2)))* &
							(heaviside(bintop(3) - Pxb(3)) 	 &
						       - heaviside(binbot(3) - Pxb(3)))
					onfaceyb =0.5d0*(sign(1.d0,binbot(2) - ri2(2))   &
						       - sign(1.d0,binbot(2) - ri1(2)))* &
							(heaviside(bintop(1) - Pyb(1))   &
						       - heaviside(binbot(1) - Pyb(1)))* &
							(heaviside(bintop(3) - Pyb(3))   &
						       - heaviside(binbot(3) - Pyb(3)))
					onfacezb =0.5d0*(sign(1.d0,binbot(3) - ri2(3))   &
						       - sign(1.d0,binbot(3) - ri1(3)))* &
							(heaviside(bintop(1) - Pzb(1))   &
						       - heaviside(binbot(1) - Pzb(1)))* &
							(heaviside(bintop(2) - Pzb(2))   &
						       - heaviside(binbot(2) - Pzb(2)))

					onfacext =0.5d0*(sign(1.d0,bintop(1) - ri2(1))   &
						       - sign(1.d0,bintop(1) - ri1(1)))* &
							(heaviside(bintop(2) - Pxt(2))   &
						       - heaviside(binbot(2) - Pxt(2)))* &
				            		(heaviside(bintop(3) - Pxt(3))   &
						       - heaviside(binbot(3) - Pxt(3)))
					onfaceyt =0.5d0*(sign(1.d0,bintop(2) - ri2(2))   &
						       - sign(1.d0,bintop(2) - ri1(2)))* &
							(heaviside(bintop(1) - Pyt(1))   &
						       - heaviside(binbot(1) - Pyt(1)))* &
							(heaviside(bintop(3) - Pyt(3))   &
						       - heaviside(binbot(3) - Pyt(3)))
					onfacezt =0.5d0*(sign(1.d0,bintop(3) - ri2(3))   &
						       - sign(1.d0,bintop(3) - ri1(3)))* &
							(heaviside(bintop(1) - Pzt(1))   &
    						       - heaviside(binbot(1) - Pzt(1)))* &
							(heaviside(bintop(2) - Pzt(2))   &
						       - heaviside(binbot(2) - Pzt(2)))

					!imaxloc = maxloc(abs(crossface))
					jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer

					!Calculate velocity at time of intersection
					!crosstime = (r(n,jxyz) - rplane)/v(n,jxyz)
					velvect(:) = v(n,:) !- a(n,:) * crosstime
					!Change in velocity at time of crossing is not needed as velocity assumed constant 
					!for timestep and changes when forces are applied.

					!Add Momentum flux over face
					momentum_flux(cbin(1),cbin(2),cbin(3),:,1) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,1) & 
					      - velvect(:)*dble(onfacexb)*abs(crossface(jxyz))
					momentum_flux(cbin(1),cbin(2),cbin(3),:,2) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,2) & 
					      - velvect(:)*dble(onfaceyb)*abs(crossface(jxyz))
					momentum_flux(cbin(1),cbin(2),cbin(3),:,3) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,3) &
					      - velvect(:)*dble(onfacezb)*abs(crossface(jxyz))
					momentum_flux(cbin(1),cbin(2),cbin(3),:,4) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,4) &
					      + velvect(:)*dble(onfacext)*abs(crossface(jxyz))
					momentum_flux(cbin(1),cbin(2),cbin(3),:,5) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,5) &
					      + velvect(:)*dble(onfaceyt)*abs(crossface(jxyz))
					momentum_flux(cbin(1),cbin(2),cbin(3),:,6) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,6) &
					      + velvect(:)*dble(onfacezt)*abs(crossface(jxyz))

				enddo
				enddo
				enddo

			endif

		enddo
	case default 
		stop "Cumulative Momentum flux Error"
	end select

end subroutine cumulative_momentum_flux

!===================================================================================
! Control Volume snapshot of momentum in a given bin

subroutine momentum_snapshot
	use module_record
	implicit none

	integer						:: n,i,j,k,ixyz,m
	integer		,dimension(3)			:: ibin
	integer		,dimension(:,:,:)  ,allocatable	:: volume_mass_temp
	double precision				:: binvolume
	double precision,dimension(3)			:: ri, mbinsize, bincentre
	double precision,dimension(:,:,:,:),allocatable :: volume_momentum_temp

	mbinsize(:) = domain(:) / nbins(:)

	!Allocate temporary array for mass and momentum in volume
	allocate(volume_mass_temp(nbins(1)+2,nbins(2)+2,nbins(3)+2))
	allocate(volume_momentum_temp(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))

	!Reset Control Volume momentum 
	volume_mass_temp = 0
	volume_momentum_temp = 0.d0
	do n = 1,np
		!Add up current volume momentum densities
		ibin(:) = ceiling((r(n,:)+halfdomain(:))/mbinsize(:)) + 1
		volume_mass_temp(ibin(1),ibin(2),ibin(3)) = volume_mass_temp(ibin(1),ibin(2),ibin(3)) + 1
		volume_momentum_temp(ibin(1),ibin(2),ibin(3),:) = volume_momentum_temp(ibin(1),ibin(2),ibin(3),:) + v(n,:)
	enddo

	binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
	volume_momentum_temp = volume_momentum_temp/binvolume

	!Output Control Volume momentum change and fluxes
	call velocity_bin_io(volume_mass_temp,volume_momentum_temp,'snap')

	deallocate(volume_mass_temp)
	deallocate(volume_momentum_temp)

end subroutine momentum_snapshot

!====================================================================================
!
!	Force based statistics called from computation of forces routine
!
!====================================================================================
! Called from force routine - adds up all molecular interactions

subroutine pressure_tensor_forces(molno, rij, accijmag)
	use module_record
	implicit none

	integer						:: ixyz, jxyz, i
	integer,intent(in)				:: molno
	double precision,intent(in)     		:: accijmag    !Non directional component of acceleration
	double precision,dimension(3),intent(in)   	:: rij         !vector between particles i and j

	do ixyz = 1,3
	do jxyz = 1,3
		rfmol(molno,ixyz,jxyz) = rfmol(molno,ixyz,jxyz) + accijmag*rij(ixyz)*rij(jxyz)
	enddo
	enddo

end subroutine pressure_tensor_forces


!====================================================================================
! Use a configurational expression, similar to the virial with interaction partitioned 
! between bins based on the share of the interaction between two molecules in a given bin
! N.B. Assume no more than 1 bin between bins containing the molecules

subroutine pressure_tensor_forces_VA(ri,rj,rij,accijmag)
	use module_record
	implicit none

	integer						:: ixyz, jxyz,i,j,k,l,n
	integer						:: diff
	integer,dimension(3)				:: ibin, jbin, bindiff 
	integer,dimension(:,:)		   ,allocatable	:: interbin, interbindiff
	double precision,intent(in)            		:: accijmag    !Non directional component of acceleration
	double precision				:: shift
	double precision,dimension(3), intent(in)	:: rij, ri, rj
	double precision,dimension(3)			:: VAbinsize, normal, p
	double precision,dimension(:,:)	   ,allocatable	:: intersection
	double precision,dimension(:,:,:)  ,allocatable	:: MLfrac !Magnitude of fraction of stress
	double precision,dimension(:,:,:,:),allocatable	:: Lfrac  !Fraction of stress in a given cell


	!================================================================
	!= Establish bins for molecules & check number of required bins	=
	!================================================================

	do ixyz = 1,nd

		VAbinsize(ixyz) = domain(ixyz) / nbins(ixyz)
		if (VAbinsize(ixyz) .lt. cellsidelength(ixyz)) stop "Binsize bigger than cellsize ~ Not ok for volume averaging"

		!Determine current bins using integer division
		ibin(ixyz) = ceiling((ri(ixyz)+halfdomain(ixyz))/VAbinsize(ixyz)) + 1 !Establish current i bin
		jbin(ixyz) = ceiling((rj(ixyz)+halfdomain(ixyz))/VAbinsize(ixyz)) + 1 !Establish current j bin

		!Check number of bins between molecules
		bindiff(ixyz) = abs(ibin(ixyz) - jbin(ixyz)) + 1

	enddo

	!================================================================
	!=			Assign to bins				=
	!================================================================

	!Check difference between bins i and j
	diff = bindiff(1)+bindiff(2)+bindiff(3)


	!Ignore values outside of domain resulting from shifted bin 
	if (minval(ibin) .lt. 1) diff = 0
	if (minval(jbin) .lt. 1) diff = 0
	if (maxval(ibin) .gt. maxval(nbins)+2) diff = 0
	if (maxval(jbin) .gt. maxval(nbins)+2) diff = 0

	select case(diff)
	!================Skip Force addition===========================
	case(0)
	!Do Nothing
	!print*, maxval(ri(:)+halfdomain(:)), 'is outside of domain', domain(1)
	!================Molecules in same bin===========================
	case(3)

		!------Add molecules to bin-----
		do ixyz = 1,nd
		do jxyz = 1,nd
			!Factor of two as molecule i and molecule j are both counted in bin i
			rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) = rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) &
				  + accijmag*rij(ixyz)*rij(jxyz)
		enddo
		enddo

	!===================Interactions split over 2 cells only==========
	case(4)

		!Allocate array for vectors from molecules to bin walls
		allocate(Lfrac(bindiff(1),bindiff(2),bindiff(3),nd))
		!Allocate arrays for magnitude of vector Lfrac
		allocate(MLfrac(bindiff(1),bindiff(2),bindiff(3)))
		!Allocate array for one bin intersection point
		allocate(intersection(3,1))
		Lfrac = 0.d0
		MLfrac = 0.d0

		do ixyz = 1,nd
			!Test to see over which coordinate molecules are in different bins
			if (bindiff(ixyz) .ne. 1) then

				!Set normal to plane for ixyz
				normal = 0.d0
				normal(ixyz) = 1.d0

				!Establish location of plane between ri and rj
				if (ri(ixyz) .lt. rj(ixyz)) then
					p(ixyz) = (ibin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				elseif (ri(ixyz) .gt. rj(ixyz)) then
					p(ixyz) = (jbin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				endif

				!Calculate location of intersection of rij and plane
				call plane_line_intersect(intersection(:,1),normal,p,ri,rj)

				!Calculate vectors from ri & rj to intersect
				Lfrac(1,1,1,:) = ri(:)-intersection(:,1)
				Lfrac(bindiff(1),bindiff(2),bindiff(3),:) = rj(:)-intersection(:,1)
				
			endif
		enddo

		!Calculate magnitude of 3 vector components
		do ixyz = 1,3
			MLfrac(:,:,:) = MLfrac(:,:,:) + Lfrac(:,:,:,ixyz)**2
		enddo
		MLfrac(:,:,:) = MLfrac(:,:,:)**0.5d0
		!Normalise to one
		MLfrac(:,:,:) = MLfrac(:,:,:)/magnitude(rij(:))

		!------Add molecules to bin-----
		!Molecule i and j contribution split between bins
		do ixyz = 1,nd
		do jxyz = 1,nd
			!-----------Molecule i bin-----------
			rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) = rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) &
				  + accijmag*rij(ixyz)*rij(jxyz)*MLfrac(1,1,1)

			!-----------Molecule j bin-----------
			rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) = rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) &
				  + accijmag*rij(ixyz)*rij(jxyz)*MLfrac(bindiff(1),bindiff(2),bindiff(3))

		enddo
		enddo

		deallocate(intersection)
		deallocate(Lfrac)
		deallocate(MLfrac)

		
	!==============Interactions split over intermediate cell===========
	case(5)

		!Allocate array for vectors from molecules to bin walls
		allocate(Lfrac(bindiff(1),bindiff(2),bindiff(3),nd))
		!Allocate arrays for magnitude of vector Lfrac
		allocate(MLfrac(bindiff(1),bindiff(2),bindiff(3)))
		!Allocate array for bin intersection points
		allocate(intersection(3,2))
		!Allocate array for intersection points
		allocate(interbin(3,1))
		allocate(interbindiff(3,1))

		!Set intersection location array to zero
		n = 1
		intersection = 0.d0
		Lfrac = 0.d0
		MLfrac = 0.d0

		do ixyz = 1,nd

			!Test to see over which coordinate molecules are in different bins
			if (bindiff(ixyz) .ne. 1) then

				!Set normal to plane for ixyz
				normal = 0.d0
				normal(ixyz) = 1.d0

				!Establish location of plane between ri and rj
				if (ri(ixyz) .lt. rj(ixyz)) then
					p(ixyz) = (ibin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				elseif (ri(ixyz) .ge. rj(ixyz)) then
					p(ixyz) = (jbin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				endif

				!Calculate location of intersection of rij and plane
				call plane_line_intersect(intersection(:,n),normal,p,ri,rj)

				n = n + 1

			endif
		enddo

		!Take average of 2 cell side intersections to determine intermediate cell
		do ixyz=1,3
			interbin(ixyz,1) = ceiling((0.5d0*(intersection(ixyz,1) &
					+intersection(ixyz,2))+halfdomain(ixyz))/VAbinsize(ixyz))+1
			interbindiff(ixyz,1) = abs(interbin(ixyz,1)-ibin(ixyz)) + 1
		enddo



		!Check which plane is closest to i and which corresponds to j
		if (magnitude(ri(:)-intersection(:,1)) .le.  & 
  	       	    magnitude(ri(:)-intersection(:,2))) then
			i = 1
			j = 2
		else
			i = 2
			j = 1
		endif

		!Fix for vectors going directly through vertex of bin
		if (all(interbindiff(:,1) .eq. 1)) then
			!Ensure not in same bin as 1
			if (bindiff(1)+bindiff(2) .eq. 3) then
				interbindiff(3,1) = 2 
			else 
				interbindiff(1,1) = 2
			endif
		endif
		if (all(interbindiff(:,1) .eq. bindiff(:))) then
			!Ensure not in same bin as bindiff
			if (bindiff(1)+bindiff(2) .eq. 3) then
				interbindiff(3,1) = 1 
			else 
				interbindiff(1,1) = 1
			endif
		endif

		!Calculate vectors from ri to intersect and rj to intersect
		Lfrac(1,1,1,:) = ri(:)-intersection(:,i)
		Lfrac(bindiff(1),bindiff(2),bindiff(3),:) = rj(:)-intersection(:,j)
		Lfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1),:) = intersection(:,i)-intersection(:,j)

		!Calculate magnitude of 3 vector components
		do ixyz = 1,3
			MLfrac(:,:,:) = MLfrac(:,:,:) + Lfrac(:,:,:,ixyz)**2
		enddo
		MLfrac(:,:,:) = MLfrac(:,:,:)**0.5d0

		!Normalise to one
		MLfrac(:,:,:) = MLfrac(:,:,:)/magnitude(rij(:))

		!------Add stress component to bin weighted by line length-----
		!Intermediate bin is either in domain or in halo. 
		!For halo bins the stress is added to the cell the halo represents. 
		!For domain cells it is added directly to that cell
		do ixyz = 1,nd
		do jxyz = 1,nd

			!-----------Molecule i bin-----------
			rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) = rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(1,1,1)

			!-----------Intermediate Bin-----------
			!If both bins are in the domain then contribution is added for both molecules
			rfbin(interbin(1,1),interbin(2,1),interbin(3,1),ixyz,jxyz) = &
				rfbin(interbin(1,1),interbin(2,1),interbin(3,1),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1))

			!-----------Molecule j bin-----------
			rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) = rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(bindiff(1),bindiff(2),bindiff(3))

		enddo
		enddo

		deallocate(Lfrac)
		deallocate(MLfrac)
		deallocate(interbin)
		deallocate(intersection)
		deallocate(interbindiff)

	!===========Interactions split over 2 intermediate cell===============
	case (6)

		!Allocate array for vectors from molecules to bin walls
		allocate(Lfrac(bindiff(1),bindiff(2),bindiff(3),nd))
		!Allocate arrays for magnitude of vector Lfrac
		allocate(MLfrac(bindiff(1),bindiff(2),bindiff(3)))
		!Allocate array for bin intersection points
		allocate(intersection(3,3))
		!Allocate array for intersection points
		allocate(interbin(3,2))
		allocate(interbindiff(3,2))

		!Set intersection location array to zero
		n = 1
		intersection = 0.d0

		Lfrac = 0.d0
		MLfrac = 0.d0

		do ixyz = 1,nd

			!Test to see over which coordinate molecules are in different bins
			if (bindiff(ixyz) .ne. 1) then

				!Set normal to plane for ixyz
				normal = 0.d0
				normal(ixyz) = 1.d0

				!Establish location of plane between ri and rj
				if (ri(ixyz) .lt. rj(ixyz)) then
					p(ixyz) = (ibin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				elseif (ri(ixyz) .gt. rj(ixyz)) then
					p(ixyz) = (jbin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				endif

				!Calculate location of intersection of rij and plane
				call plane_line_intersect(intersection(:,n),normal,p,ri,rj)

				n = n+1

			endif
		enddo

		!Determine which intersection covers both intermediate cells
		! |(1)-(2)| > |(1)-(3)| > |(2)-(3)| then 1,2 must cover
		! both cells while 1,3 and 2,3 are the cell intercepts

		if (magnitude(intersection(:,1)-intersection(:,2)) .gt. & 
		    magnitude(intersection(:,3)-intersection(:,2))) then
			if (magnitude(intersection(:,1)-intersection(:,2)) .gt. & 
			    magnitude(intersection(:,1)-intersection(:,3))) then
				k = 1
				l = 2
			else
				k = 1
				l = 3
			endif
		else
			if (magnitude(intersection(:,3)-intersection(:,2)) .gt. & 
			    magnitude(intersection(:,1)-intersection(:,3))) then
				k = 2
				l = 3
			else
				k = 1
				l = 3
			endif
		endif

		!Take average of cell side intersections to determine intermediate bins
		!k and l are used to determine intermediate bins
		do ixyz=1,3
		
			interbin(ixyz,1) = ceiling((0.5d0*(intersection(ixyz,(6-k-l)) &
					    +intersection(ixyz,l))+halfdomain(ixyz))/VAbinsize(ixyz))+1
			interbindiff(ixyz,1) = abs(interbin(ixyz,1)-ibin(ixyz)) + 1

			interbin(ixyz,2) = ceiling((0.5d0*(intersection(ixyz,(6-k-l)) &
					    +intersection(ixyz,k))+halfdomain(ixyz))/VAbinsize(ixyz))+1
			interbindiff(ixyz,2) = abs(interbin(ixyz,2)-ibin(ixyz)) + 1

		enddo

		!Check which plane is closest to i 
		if (magnitude(ri(:)-intersection(:,1)) .lt.  & 
  	       	    magnitude(ri(:)-intersection(:,2))) then
			if (magnitude(ri(:)-intersection(:,1)) .lt.  & 
  		       	    magnitude(ri(:)-intersection(:,3))) then
				i = 1
			else
				i = 3
			endif
		else
			if (magnitude(ri(:)-intersection(:,2)) .lt.  & 
  	       	  	    magnitude(ri(:)-intersection(:,3))) then
				i = 2
			else
				i = 3
			endif
		endif

		!Check which plane is closest to j 
		if (magnitude(rj(:)-intersection(:,1)) .lt.  & 
  	       	    magnitude(rj(:)-intersection(:,2))) then
			if (magnitude(rj(:)-intersection(:,1)) .lt.  & 
  		       	    magnitude(rj(:)-intersection(:,3))) then
				j = 1
			else
				j = 3
			endif
		else
			if (magnitude(rj(:)-intersection(:,2)) .lt.  & 
  	       	  	    magnitude(rj(:)-intersection(:,3))) then
				j = 2
			else
				j = 3
			endif
		endif

		!Calculate vectors from ri to intersect & rj to intersect
		Lfrac(1,1,1,:) = ri(:)-intersection(:,i)
		Lfrac(bindiff(1),bindiff(2),bindiff(3),:) = rj(:)-intersection(:,j)

		!Calculate vectors in two intermediate cells
		Lfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1),:) = intersection(:,l)-intersection(:,(6-k-l))
		Lfrac(interbindiff(1,2),interbindiff(2,2),interbindiff(3,2),:) = intersection(:,k)-intersection(:,(6-k-l))

		!Calculate magnitude of 3 vector components
		do ixyz = 1,3
			MLfrac(:,:,:) = MLfrac(:,:,:) + Lfrac(:,:,:,ixyz)**2
		enddo
		MLfrac(:,:,:) = MLfrac(:,:,:)**0.5d0

		!Normalise to one
		MLfrac(:,:,:) = MLfrac(:,:,:)/magnitude(rij(:))

		!------Add stress component to bin weighted by line length-----
		!Intermediate bin is either in domain or in halo. 
		!For halo bins the stress is added to the cell the halo represents. 
		!For domain cells it is added directly to that cell
		do ixyz = 1,nd
		do jxyz = 1,nd
			!-----------Molecule i bin-----------
			rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) = rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(1,1,1)

			!-----------1st Intermediate Bin-----------
			!Intermediate and i bin in domain, j in halo - add intermediate for molecule i
			rfbin(interbin(1,1),interbin(2,1),interbin(3,1),ixyz,jxyz) = &
				rfbin(interbin(1,1),interbin(2,1),interbin(3,1),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1))

			!-----------2nd Intermediate Bin-----------
			!Intermediate and i bin in domain, j in halo - add intermediate for molecule i
			rfbin(interbin(1,2),interbin(2,2),interbin(3,2),ixyz,jxyz) = &
				rfbin(interbin(1,2),interbin(2,2),interbin(3,2),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(interbindiff(1,2),interbindiff(2,2),interbindiff(3,2))

			!-----------Molecule j bin-----------
			rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) = rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(bindiff(1),bindiff(2),bindiff(3))
		enddo
		enddo

		deallocate(intersection)
		deallocate(Lfrac)
		deallocate(MLfrac)
		deallocate(interbin)
		deallocate(interbindiff)
		
	case default

		stop "VOLUME AVERAGING ERROR"

	end select

end subroutine pressure_tensor_forces_VA

!===================================================================================
!Forces over the surface of a Volume

subroutine control_volume_forces(fij,ri,rj,molnoi,molnoj)
use module_record
implicit none

	integer				:: ixyz, n, molnoi, molnoj
	integer,dimension(3)		:: ibin, jbin
	double precision		:: binforce
	double precision,dimension(3)	:: ri, rj, fij,crossplane,fsurface
	double precision,dimension(3)	:: Fbinsize, bintopi, binboti, bintopj, binbotj

	!Determine bin size
	Fbinsize(:) = domain(:) / nbins(:)

	!Assign to bins using integer division
	ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+1	!Establish current bin
	jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+1 	!Establish current bin

	crossplane(:) =  dble(ibin(:)-jbin(:))

	bintopi(:) = (ibin(:)-1)*Fbinsize(:)-halfdomain(:)
	binboti(:) = (ibin(:)-2)*Fbinsize(:)-halfdomain(:)
	bintopj(:) = (jbin(:)-1)*Fbinsize(:)-halfdomain(:)
	binbotj(:) = (jbin(:)-2)*Fbinsize(:)-halfdomain(:)

	!Add for molecule i
	if(molnoi .le. np) then
		fsurface = fij(:)* dble((heaviside(bintopi(1)-ri(1))-heaviside(binboti(1)-ri(1)))* & 
			  		(heaviside(bintopi(2)-ri(2))-heaviside(binboti(2)-ri(2)))* & 
			  		(heaviside(bintopi(3)-ri(3))-heaviside(binboti(3)-ri(3)))- & 
			  		(heaviside(bintopi(1)-rj(1))-heaviside(binboti(1)-rj(1)))* & 
			  		(heaviside(bintopi(2)-rj(2))-heaviside(binboti(2)-rj(2)))* & 
			  		(heaviside(bintopi(3)-rj(3))-heaviside(binboti(3)-rj(3))))
		volume_force(ibin(1),ibin(2),ibin(3),:,1) = volume_force(ibin(1),ibin(2),ibin(3),:,1) + fsurface*delta_t
	endif

	!Add for molecule j
	if(molnoj .le. np) then
		fsurface = fij(:)* dble((heaviside(bintopj(1)-ri(1))-heaviside(binbotj(1)-ri(1)))* & 
			  		(heaviside(bintopj(2)-ri(2))-heaviside(binbotj(2)-ri(2)))* & 
			  		(heaviside(bintopj(3)-ri(3))-heaviside(binbotj(3)-ri(3)))- & 
			  		(heaviside(bintopj(1)-rj(1))-heaviside(binbotj(1)-rj(1)))* & 
			  		(heaviside(bintopj(2)-rj(2))-heaviside(binbotj(2)-rj(2)))* & 
			  		(heaviside(bintopj(3)-rj(3))-heaviside(binbotj(3)-rj(3))))
		volume_force(jbin(1),jbin(2),jbin(3),:,1) = volume_force(jbin(1),jbin(2),jbin(3),:,1) + fsurface*delta_t
	endif

end subroutine control_volume_forces

!===================================================================================
!Forces over the surface of a Volume

subroutine control_volume_stresses(fij,ri,rj,molnoi,molnoj)
use module_record
implicit none

	integer				:: i,j,k,ixyz,n,molnoi,molnoj,tempi
	integer				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	integer,dimension(3)		:: cbin, ibin, jbin
	double precision		:: binforce
	double precision,dimension(3)	:: ri,rj,rij,fij,fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
	double precision,dimension(3)	:: Fbinsize, bintop, binbot, temp

	!Calculate rij
	rij = ri - rj
	!Prevent Division by zero
	do ixyz = 1,3
		if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
	enddo

	!Determine bin size
	Fbinsize(:) = domain(:) / nbins(:)

	!Assign to bins using integer division
	ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+1	!Establish current bin
	jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+1 	!Establish current bin

	do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
	do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
	do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

		cbin(1) = i; cbin(2) = j; cbin(3) = k

		bintop(:) = (cbin(:)-1)*Fbinsize(:)-halfdomain(:)
		binbot(:) = (cbin(:)-2)*Fbinsize(:)-halfdomain(:)

		!Calculate the plane intersect of line with surfaces of the cube
		Pxt=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
		Pxb=(/ binbot(1),ri(2)+(rij(2)/rij(1))*(binbot(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(binbot(1)-ri(1))  /)
		Pyt=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
		Pyb=(/ri(1)+(rij(1)/rij(2))*(binbot(2)-ri(2)), binbot(2),ri(3)+(rij(3)/rij(2))*(binbot(2)-ri(2))  /)
		Pzt=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)
		Pzb=(/ri(1)+(rij(1)/rij(3))*(binbot(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(binbot(3)-ri(3)), binbot(3) /)

		onfacexb =  	(sign(1.d0,binbot(1)- rj(1)) - sign(1.d0,binbot(1)- ri(1)))* &
				(heaviside(bintop(2)-Pxb(2)) - heaviside(binbot(2)-Pxb(2)))* &
				(heaviside(bintop(3)-Pxb(3)) - heaviside(binbot(3)-Pxb(3)))
		onfaceyb =  	(sign(1.d0,binbot(2)- rj(2)) - sign(1.d0,binbot(2)- ri(2)))* &
				(heaviside(bintop(1)-Pyb(1)) - heaviside(binbot(1)-Pyb(1)))* &
				(heaviside(bintop(3)-Pyb(3)) - heaviside(binbot(3)-Pyb(3)))
		onfacezb =  	(sign(1.d0,binbot(3)- rj(3)) - sign(1.d0,binbot(3)- ri(3)))* &
				(heaviside(bintop(1)-Pzb(1)) - heaviside(binbot(1)-Pzb(1)))* &
				(heaviside(bintop(2)-Pzb(2)) - heaviside(binbot(2)-Pzb(2)))

		onfacext =  	(sign(1.d0,bintop(1)- rj(1)) - sign(1.d0,bintop(1)- ri(1)))* &
				(heaviside(bintop(2)-Pxt(2)) - heaviside(binbot(2)-Pxt(2)))* &
	            		(heaviside(bintop(3)-Pxt(3)) - heaviside(binbot(3)-Pxt(3)))
		onfaceyt = 	(sign(1.d0,bintop(2)- rj(2)) - sign(1.d0,bintop(2)- ri(2)))* &
				(heaviside(bintop(1)-Pyt(1)) - heaviside(binbot(1)-Pyt(1)))* &
				(heaviside(bintop(3)-Pyt(3)) - heaviside(binbot(3)-Pyt(3)))
		onfacezt =  	(sign(1.d0,bintop(3)- rj(3)) - sign(1.d0,bintop(3)- ri(3)))* &
				(heaviside(bintop(1)-Pzt(1)) - heaviside(binbot(1)-Pzt(1)))* &
				(heaviside(bintop(2)-Pzt(2)) - heaviside(binbot(2)-Pzt(2)))

		!Prevent halo molecules from being included but include molecule which have left domain 
		!before rebuild has been triggered.
		if (molnoi .gt. np .or. molnoj .gt. np) then
			if (cbin(1) .lt. 2 .or. cbin(1) .gt. nbins(1)+1) cycle
			if (cbin(2) .lt. 2 .or. cbin(2) .gt. nbins(2)+1) cycle
			if (cbin(3) .lt. 2 .or. cbin(3) .gt. nbins(3)+1) cycle
		endif

		!Stress acting on face over volume
		Pxyface(cbin(1),cbin(2),cbin(3),:,1) = Pxyface(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*dble(onfacexb)
		Pxyface(cbin(1),cbin(2),cbin(3),:,2) = Pxyface(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*dble(onfaceyb)
		Pxyface(cbin(1),cbin(2),cbin(3),:,3) = Pxyface(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*dble(onfacezb)
		Pxyface(cbin(1),cbin(2),cbin(3),:,4) = Pxyface(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*dble(onfacext)
		Pxyface(cbin(1),cbin(2),cbin(3),:,5) = Pxyface(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*dble(onfaceyt)
		Pxyface(cbin(1),cbin(2),cbin(3),:,6) = Pxyface(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*dble(onfacezt)

		!Force applied to volume
		fsurface(:) = 0.d0
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfacexb - onfacext)
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfaceyb - onfaceyt)
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfacezb - onfacezt)
		volume_force(cbin(1),cbin(2),cbin(3),:,1) = volume_force(cbin(1),cbin(2),cbin(3),:,1) + fsurface*delta_t

	enddo
	enddo
	enddo

end subroutine control_volume_stresses

!====================================================================================

subroutine pressure_tensor_forces_MOP(pnxyz,ri,rj,rij,accijmag)
	use module_record
	implicit none

	integer,save			:: i,j,k
	integer				:: n
	integer				:: pnxyz	 !Plane normal direction
	integer				:: molno, planenoi,planenoj
	double precision                :: shift, plane !Plane normal components i and j
	double precision                :: accijmag      !Non directional component of acceleration
	double precision,dimension(3)   :: ri, rj, rij   !Vector between particles i and j
	double precision,dimension(3)   :: Pyb           !Location of intercept with plane

	!Shift by half difference between value rounded down and actual value
	!to ensure same distance to top and bottom plane from domain edge
	shift = 0.5d0*(domain(pnxyz)-planespacing*(nplanes-1))

	!Obtain nearest plane number by integer division
	planenoi = nint((ri(2)+halfdomain(pnxyz)-shift)/planespacing)+1
	if (planenoi .lt. 1)	   planenoi = 1
	if (planenoi .gt. nplanes) planenoi = nplanes
	planenoj = nint((rj(2)+halfdomain(pnxyz)-shift)/planespacing)+1
	if (planenoj .lt. 1)	   planenoj = 1
	if (planenoj .gt. nplanes) planenoj = nplanes

	!Use calculated plane numbers and check all planes between
	do n = planenoi,planenoj,sign(1,planenoj-planenoi)
		plane = planes(n)

		!Caluclate intersection with plane to check if interaction is within domain
		Pyb=(/ri(1)+(rij(1)/rij(2))*(plane-ri(2)), plane,ri(3)+(rij(3)/rij(2))*(plane-ri(2)) /)

		!Using sign function (Heaviside function is less efficient as no inbuilt Fortran Heaviside)
		Pxy_plane(:,n) =  Pxy_plane(:,n) + accijmag * rij(:) * & 
				( sign(1.d0,ri(2)-plane)-sign(1.d0,rj(2)-plane) )* &
				(heaviside(halfdomain(1)-Pyb(1)) - heaviside(-halfdomain(1)-Pyb(1)))* &
				(heaviside(halfdomain(3)-Pyb(3)) - heaviside(-halfdomain(3)-Pyb(3)))

	enddo

end subroutine pressure_tensor_forces_MOP


!===================================================================================
!Calculate Radial distribution function (RDF) using cell method
!subroutine evaluate_properties_cellradialdist(molnoi)
!use module_record
!implicit none
!
!	integer                         :: j, ixyz   !Define dummy index
!	integer                         :: icell, jcell
!	integer                         :: icellshift, jcellshift
!	integer                         :: adjacentcellnp
!	integer				:: molnoi, molnoj
!	type(node), pointer 	        :: oldj, currentj

	!Find cell location of specified molecules
!	icell = ceiling((r(molnoi,1)+halfdomain(1))/cellsidelength(1))+1 !Add 1 due to halo
!	jcell = ceiling((r(molnoi,2)+halfdomain(2))/cellsidelength(2))+1 !Add 1 due to halo
!	ri = r(molnoi,:)

!	do jcellshift = -1,1
!	do icellshift = -1,1
!		oldj => cell%head(icell+icellshift,jcell+jcellshift)%point
!		adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift)
!
!		do j = 1,adjacentcellnp          !Step through all j for each i
!			rij2=0                   !Set rij^2 to zero
!			molnoj = oldj%molno !Number of molecule
!			rj = r(molnoj,:)         !Retrieve rj
!					
!			currentj => oldj
!			oldj => currentj%next    !Use pointer in datatype to obtain next item in list

!			if(molnoi==molnoj) cycle !Check to prevent interaction with self

!			do ixyz=1,nd
!				rij(ixyz) = ri(ixyz) - rj(ixyz)                 !Evaluate distance between particle i and j
!				if (abs(rij(ixyz)) > halfdomain(ixyz)) then  !N.B array operator-corresponding element compared 
!					rij(ixyz) = rij(ixyz) - sign(domain(ixyz),rij(ixyz))  !Sign of rij applied to domain
!				endif
!				rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
!			enddo
!		enddo
!	enddo
!	enddo

!	nullify(currentj)        !Nullify current as no longer required
!	nullify(oldj)            !Nullify old as no longer required

!end subroutine evaluate_properties_cellradialdist
