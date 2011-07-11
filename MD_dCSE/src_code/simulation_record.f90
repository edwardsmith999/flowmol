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
! simulation_record  !Top level function choosing logging based on inputs
! evaluate_macroscopic_properties_parallel !Macroscopic properties
! evaluate_properties_vdistribution
! evaluate_properties_radialdist
! evaluate_properties_diffusion
! velocity_contour
! velocity_slice
! 	cumulative_bin_velocity
! 	record_cumulative_velocity
! pressure_tensor
! 	pressure_tensor_forces
! 	cumulative_pressure_tensor
! 	record_pressure_tensor
! pressure_tensor_MOP
! 	pressure_tensor_forces_MOP
! 	cumulative_pressure_tensor_MOP
!	record_pressure_tensor_MOP
!
!==========================================================================

module module_record

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use calculated_properties_MD

	double precision :: vel

end module module_record
!==========================================================================


subroutine simulation_record
use module_record
implicit none

	!Parallel output for molecular positions
	if (vmd_outflag .eq. 1) call parallel_io_vmd
	if (vmd_outflag .eq. 2) call parallel_io_vmd_sl
	if (vmd_outflag .eq. 3) call parallel_io_vmd_halo

	!call evaluate_macroscopic_properties
	!Obtain each processe's subdomain's macroscopic 
	!properties; gather on root process and record
	if (macro_outflag .eq. 1) call evaluate_macroscopic_properties_parallel

	!Obtain and record velocity distributions
	!call evaluate_properties_vdistribution

	!Obtain and record radial distributions
	!call evaluate_properties_radialdist

	!Obtain and record veolcity contour
	!call velocity_contour

	!Obtain and record velocity across input direction
	if (slice_outflag .eq. 1) call velocity_slice(2)
	if (slice_outflag .eq. 2) call velocity_averaging
	if (slice_outflag .eq. 3) then
		call velocity_averaging
		call velocity_slice(2)
	endif
	!Obtain and record molecular diffusion
	!call evaluate_properties_diffusion

	!Calculate pressure tensor
	if (pressure_outflag .eq. 1) call pressure_tensor
	if (pressure_outflag .eq. 2) call pressure_tensor_VA
	if (pressure_outflag .eq. 3) call pressure_tensor_MOP

	!Check if backup file should be saved based on elapsed iterations
	!and average size of current system per processor
	if (mod((iter*globalnp/nproc),1000000000).eq.0) then
		call messenger_syncall
		call parallel_io_final_state
	endif

end subroutine simulation_record

!==========================================================================
!Calculate kinetic and potential energy as well as temperature and pressure

subroutine evaluate_macroscopic_properties_parallel
use module_record
implicit none

	integer          :: n, ixyz

	vsum  = 0.d0      ! Reset all sums
	v2sum = 0.d0      ! Reset all sums

	do n = 1, np    ! Loop over all particles
		potenergysum = potenergysum + potenergymol(n)
		virial = virial + virialmol(n)
		do ixyz = 1, nd    ! Loop over all dimensions
			!Velocity component must be shifted forward half a timestep to determine 
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
		totenergy   = kinenergy + potenergy
		temperature = v2sum / real(nd*globalnp,kind(0.d0))
		pressure    = (density/(globalnp*nd))*(v2sum+virial/2) !N.B. virial/2 as all interactions calculated

		print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f10.4)', &
		iter,';',vsum,';', v2sum,';', temperature,';', &
		kinenergy,';',potenergy,';',totenergy,';',pressure

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
	!integer,dimension(

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
!		RECORD VELOCITY AT LOCATION IN SPACE
! Either by binning, taking slices or on a plane, the velocity field in the molecular
! system is recorded and output
!===================================================================================
!===================================================================================
!Calculate averaged velocity components of each bin

subroutine velocity_averaging
use module_record
use linked_list
implicit none

	integer			:: ixyz
	integer			:: ibin, jbin, kbin
	integer, save		:: average_count

	average_count = average_count + 1
	call cumulative_bin_velocity_3D
	if (average_count .eq. Nslice_ave) then
		average_count = 0

		!Scale Mean velocity based on number of molecules
		do ixyz = 1,3
			!Add 1 to prevent division by zero
			meanvelbin(:,:,:,ixyz) = meanvelbin(:,:,:,ixyz)/(countvelbin(:,:,:)+1)
		enddo

		if (irank .eq. iroot) print*, 'Recording Velocity Field @ iter', iter
		call velocity_average_io

		!Reset velocity records
		meanvelbin  = 0.d0
		countvelbin = 0

		!Collect velocities for next step
		call cumulative_bin_velocity_3D

	endif

end subroutine velocity_averaging

!===================================================================================
!Add velocities to running total

subroutine cumulative_bin_velocity_3D
use module_record
use linked_list
implicit none

	integer         		:: n, ixyz
	integer 			:: ibin, jbin, kbin
	double precision,dimension(3) 	:: slicebinsize

	slicebinsize(1) = domain(1) / nbins(1)
	slicebinsize(2) = domain(2) / nbins(2)
	slicebinsize(3) = domain(3) / nbins(3)

	do   n 	= 1, np    ! Loop over all particles
		!Assign to bins using integer division
		ibin = ceiling((r(n,1)+halfdomain(1))/slicebinsize(1)) 	!Establish current bin
		if (ibin > nbins(1)) ibin = nbins(1) 			!Prevents out of range values
		if (ibin < 1 ) ibin = 1        				!Prevents out of range values
		jbin = ceiling((r(n,2)+halfdomain(2))/slicebinsize(2)) 	!Establish current bin
		if (jbin > nbins(2)) jbin = nbins(2) 			!Prevents out of range values
		if (jbin < 1 ) jbin = 1        				!Prevents out of range values
		kbin = ceiling((r(n,3)+halfdomain(3))/slicebinsize(3)) 	!Establish current bin
		if (kbin > nbins(3)) kbin = nbins(3) 			!Prevents out of range values
		if (kbin < 1 ) kbin = 1        				!Prevents out of range values
		!call linklist_checkpush_bin(ibin,jbin,kbin,n) !Adds molecule to bin
		countvelbin(ibin,jbin,kbin)  = countvelbin(ibin,jbin,kbin) + 1     	!Add one to current bin
		meanvelbin(ibin,jbin,kbin,:) = meanvelbin(ibin,jbin,kbin,:)+ v(n,:) 	!Add streamwise velocity to current bin
	enddo

end subroutine cumulative_bin_velocity_3D


!===================================================================================
!Calculate Velocity magnitude of each cell to allow contour plot

subroutine velocity_contour
use module_record
use linked_list
implicit none

	integer             :: j
	integer             :: molno, cellnp
	integer             :: icell, jcell, kcell
	double precision, dimension(ncells(1)+2,ncells(2)+2,ncells(3)+2) :: vcontour
	type(node), pointer :: old, current

	write(14,'(i10)'), ncells(1)
	write(14,'(i10)'), ncells(2)
	write(14,'(i10)'), ncells(3)

	vcontour = 0.d0

	do icell = 2, ncells(1)+1
	do jcell = 2, ncells(2)+1
	do kcell = 2, ncells(3)+1
		!Obtain molecular number and top item of link list from cell head
		old => cell%head(icell,jcell,kcell)%point
		cellnp = cell%cellnp(icell,jcell, kcell)

		if(cellnp == 0) vcontour(icell, jcell, kcell) = 0.d0

		current => old ! make current point to head of list
		do j=1,cellnp
			molno = current%molno
			vcontour(icell, jcell, kcell) = vcontour(icell, jcell, kcell) + vmagnitude(molno)
			if (associated(old%next) .eqv. .true. ) then !Exit if null
				old => current%next ! Use pointer in datatype to obtain next item in list
				current => old      ! make current point to old - move alone one
			endif
		enddo
		vcontour(icell, jcell, kcell) = vcontour(icell, jcell, kcell)/ cellnp

		write(14,'(f10.5)'), vcontour(icell, jcell, kcell)
	enddo
	enddo
	enddo

	nullify(current)                    !Nullify current as no longer required
	nullify(old)                        !Nullify old as no longer required

end subroutine velocity_contour

!===================================================================================
!Calculate Velocities of each cell to allow velocity profile - assumes uniform in
!Streamwise and spanwise directions; change in wall normal only

subroutine velocity_slice(ixyz)
use module_record
implicit none

	integer		:: ixyz
	integer, save	:: sample_count

	call cumulative_bin_velocity(ixyz)
	sample_count = sample_count + 1
	if (sample_count .eq. Nslice_ave) then
		if (irank .eq. iroot) print*, 'recording velocity @ iter', iter
		call record_cumulative_velocity(ixyz)
		sample_count = 0
	endif

end subroutine velocity_slice

!===================================================================================
!Add velocities to running total

subroutine cumulative_bin_velocity(ixyz)
use module_record
use linked_list
implicit none

	integer         	:: n, ixyz
	integer         	:: cbin
	integer 		:: ibin, jbin, kbin
	double precision 	:: slicebinsize

	!Deallocate all previous bin's linklists
	call linklist_deallocate_bins

	slicebinsize = domain(ixyz) / nbins(ixyz)

	ibin = 1	!Set to 1 as only binning in y direction
	kbin = 1	!Set to 1 as only binning in y direction

	do n = 1, np    ! Loop over all particles
		!Assign to bins using integer division
		cbin = ceiling((r(n,ixyz)+halfdomain(ixyz))/slicebinsize) !Establish current bin
		if (cbin > nbins(ixyz)) cbin = nbins(ixyz) 		!Prevents out of range values
		if (cbin < 1 ) cbin = 1        		!Prevents out of range values
		call linklist_checkpush_bin(ibin,cbin,kbin,n) !Adds molecule to bin
		meanvel(cbin) = meanvel(cbin)+v(n,1) 	!Add streamwise velocity to current bin
		countvel(cbin)= countvel(cbin)+1      	!Add one to current bin
	enddo

end subroutine cumulative_bin_velocity

!===================================================================================
!Write out running total and reset binned velocity

subroutine record_cumulative_velocity(ixyz)
use module_record
use linked_list
implicit none

	integer		:: ixyz, jxyz, kxyz

	!Get two directions orthogonal to slice direction
	kxyz = mod(ixyz,3)+1
	jxyz = mod(ixyz+1,3)+1

	!Sum over all bins using directional sub communicators
	call SubcommSumVect(meanvel, nbins(jxyz), jxyz)
	call SubcommSumVect(meanvel, nbins(kxyz), kxyz)
	call SubcommSumIntVect(countvel, nbins(jxyz), jxyz)
	call SubcommSumIntVect(countvel, nbins(kxyz), kxyz)

	!Write for processes
	call parallel_slice_io

	!Reset bins
	meanvel = 0.d0
	countvel = 0

end subroutine record_cumulative_velocity

!===================================================================================
!		CONTROL VOLUME FORMULATION OF MASS AND STRESS
! Based on the control volume formulation of mass and momentum, it is possible to
! define change in properties in terms of surface fluxes over a volume in space
! Fluxes need to be called every timestep!!!
!===================================================================================

!===================================================================================
! Control Volume mass continuity
!-----------------------------------------------------------------------------------

subroutine control_volume_mass
use module_record
implicit none

	integer, save	:: sample_count

	call control_volume_mass_flux
	sample_count = sample_count + 1
	if (sample_count .eq. Nslice_ave) then
		call control_volume_mass_evolution
		sample_count = 0
		mass_flux = 0
	endif

end subroutine control_volume_mass

!===================================================================================
! Control Volume evolution of mass in a given bin

subroutine control_volume_mass_evolution
use module_record
implicit none

	integer						:: n,i,j,k
	integer,dimension(:,:,:),allocatable 		:: volume_mass_temp
	integer		,dimension(3)			:: ibin
	double precision,dimension(3)			:: ri, mbinsize

	!Allocate temporary array for mass in volume
	allocate(volume_mass_temp(nbins(1)+2,nbins(2)+2,nbins(3)+2))

	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	!Create copy of old Control Volume mass
	volume_mass_temp = volume_mass

	!Reset Control Volume mass 
	volume_mass = 0
	do n = 1,np
		!Add up current volume mass densities
		ibin(:) = ceiling((r(n,:)+halfdomain(:))/mbinsize(:)) + 1
		volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + 1
	enddo

	!Calculate change in volume mass
	volume_mass_temp = volume_mass - volume_mass_temp

	!Output Control Volume mass change and fluxes
	!call record_volume_mass(volume_mass_temp)
	
	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nsurfacecells
		i = surfacecell(n,1); j = surfacecell(n,2); k = surfacecell(n,3)  

		!Flux over halo cells
		mass_flux(modulo((i-2),nbins(1))+2,modulo((j-2),nbins(2))+2,modulo((k-2),nbins(3))+2,:) = & 
		mass_flux(modulo((i-2),nbins(1))+2,modulo((j-2),nbins(2))+2,modulo((k-2),nbins(3))+2,:) + mass_flux(i,j,k,:)
		!Change in number of Molecules in halo cells
		volume_mass_temp(modulo((i-2),nbins(1))+2,modulo((j-2),nbins(2))+2,modulo((k-2),nbins(3))+2) = & 
		volume_mass_temp(modulo((i-2),nbins(1))+2,modulo((j-2),nbins(2))+2,modulo((k-2),nbins(3))+2) + volume_mass_temp(i,j,k)
	enddo

	!print*, sum(mass_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:))

	do i = 2,nbins(3)+1
	do j = 2,nbins(1)+1
	do k = 2,nbins(2)+1
		if(sum(mass_flux(i,j,k,:))-volume_mass_temp(i,j,k) .ne. 0) print'(a,6i8)', & 
				'Error in mass flux', iter,i,j,k,sum(mass_flux(i,j,k,:)),volume_mass_temp(i,j,k)
	enddo
	enddo
	enddo

	deallocate(volume_mass_temp)

end subroutine control_volume_mass_evolution

!===================================================================================
! Mass Flux over a surface of a bin

subroutine control_volume_mass_flux
use module_record
implicit none

	integer				:: ixyz, heaviside, n
	integer		,dimension(1)	:: imaxloc
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
			imaxloc = maxloc(abs(crossplane))
			ixyz = imaxloc(1)

			!Add mass flux to the new bin surface count and take from the old
			mass_flux(ibin1(1),ibin1(2),ibin1(3),ixyz+3*heaviside(dble(crossplane(ixyz)))) = & 
				mass_flux(ibin1(1),ibin1(2),ibin1(3),ixyz+3*heaviside(dble(crossplane(ixyz)))) & 
				+ abs(crossplane(ixyz))
			mass_flux(ibin2(1),ibin2(2),ibin2(3),ixyz+3*heaviside(-dble(crossplane(ixyz)))) = & 
				mass_flux(ibin2(1),ibin2(2),ibin2(3),ixyz+3*heaviside(-dble(crossplane(ixyz)))) &
				- abs(crossplane(ixyz))
		endif

	enddo

end subroutine control_volume_mass_flux

!===================================================================================
! Record change in mass flux

subroutine record_volume_mass(dvolume_massdt)
use module_record
implicit none

	double precision	:: dvolume_massdt

	print*, 'This is where the change in mass will be recorded'

end subroutine record_volume_mass

!===================================================================================
! Control Volume Momentum continuity
!===================================================================================

subroutine control_volume_momentum
use module_record
implicit none

	integer			:: ixyz
	integer, save		:: sample_count
	double precision	:: binface

	call control_volume_momentum_flux
	sample_count = sample_count + 1
	!Copy volume force to storage array ready to collect next measurement
	volume_force(:,:,:,:,2) = volume_force(:,:,:,:,1)
	if (sample_count .eq. Nslice_ave) then
		call control_volume_momentum_evolution
		!Calculate bin face area and corresponding stress
		Pxyface = Pxyface / sample_count
		do ixyz = 1,3
			binface = (domain(modulo(ixyz,3)+1)/nbins(modulo(ixyz,3)+1))* & 
				  (domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
			!print*, binface, domain(modulo(ixyz,3)+1)*domain(modulo(ixyz+1,3)+1)
			Pxyface(:,:,:,:,ixyz)   = Pxyface(:,:,:,:,ixyz  )/(2.d0*binface) !Bottom
			Pxyface(:,:,:,:,ixyz+3) = Pxyface(:,:,:,:,ixyz+3)/(2.d0*binface) !Top
		enddo
		!Record Control Volume Results
		call Control_Volume_io

		sample_count = 0
		momentum_flux = 0.d0
		volume_force(:,:,:,:,1) = 0.d0
		Pxyface = 0.d0
	endif

end subroutine control_volume_momentum

!===================================================================================
! Control Volume evolution of momentum in a given bin

subroutine control_volume_momentum_evolution
use module_record
implicit none

	integer						:: n,i,j,k, ixyz,m
	integer		,dimension(3)			:: ibin
	double precision,dimension(3)			:: ri, mbinsize
	double precision,dimension(:,:,:,:),allocatable :: volume_momentum_temp

	!Allocate temporary array for momentum in volume
	allocate(volume_momentum_temp(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))

	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	!Create copy of old Control Volume momentum
	volume_momentum_temp = volume_momentum

	!Reset Control Volume momentum 
	volume_momentum = 0.d0
	do n = 1,np
		!Add up current volume momentum densities
		ibin(:) = ceiling((r(n,:)+halfdomain(:))/mbinsize(:)) + 1
		volume_momentum(ibin(1),ibin(2),ibin(3),:) = volume_momentum(ibin(1),ibin(2),ibin(3),:) + v(n,:)
		!print'(4i8,3f10.5)', n, ibin(:), v(n,:)
	enddo

	!Calculate change in volume momentum
	volume_momentum_temp = volume_momentum - volume_momentum_temp

	!Output Control Volume momentum change and fluxes
	!call record_volume_momentum(volume_momentum_temp)
	
	!Include halo surface fluxes to get correct values for all cells
	do n = 1, nsurfacecells
		i = surfacecell(n,1); j = surfacecell(n,2); k = surfacecell(n,3)  
		!Flux over halo cells
		momentum_flux(	modulo((i-2),nbins(1))+2, & 
			      	modulo((j-2),nbins(2))+2, & 
			      	modulo((k-2),nbins(3))+2,:,:) = & 
				momentum_flux(	modulo((i-2),nbins(1))+2,& 
						modulo((j-2),nbins(2))+2,&
						modulo((k-2),nbins(3))+2,:,:) & 
							+ momentum_flux(i,j,k,:,:)

		!Set forces to value of halo cells
		volume_force(	modulo((i-2),nbins(1))+2, & 
			      	modulo((j-2),nbins(2))+2, & 
			      	modulo((k-2),nbins(3))+2,:,2) = & 
				volume_force(	modulo((i-2),nbins(1))+2,& 
						modulo((j-2),nbins(2))+2,&
						modulo((k-2),nbins(3))+2,:,2) & 
							+ volume_force(i,j,k,:,2)

		!Change in Momentum in halo cells
		volume_momentum_temp(	modulo((i-2),nbins(1))+2, & 
			      	modulo((j-2),nbins(2))+2, & 
			      	modulo((k-2),nbins(3))+2,:) = & 
				volume_momentum_temp(	modulo((i-2),nbins(1))+2,& 
						modulo((j-2),nbins(2))+2,&
						modulo((k-2),nbins(3))+2,:) & 
							+ volume_momentum_temp(i,j,k,:)
	enddo

	do i = 2,nbins(3)+1
	do j = 2,nbins(1)+1
	do k = 2,nbins(2)+1

		if (abs(sum(momentum_flux(i,j,k,:,1))+volume_force(i,j,k,1,2)-volume_momentum_temp(i,j,k,1)) .gt. 0.0001d0) then
		!if (i .eq. 4 .and. j .eq. 4 .and. k .eq. 4) then
		print'(a,4i5,9f10.5)','Error in momentum flux ', &
				    iter,i,j,k, sum(momentum_flux(i,j,k,:,1)),volume_force(i,j,k,1,2),volume_momentum_temp(i,j,k,1) & 
					       ,sum(momentum_flux(i,j,k,:,2)),volume_force(i,j,k,2,2),volume_momentum_temp(i,j,k,2) & 
					       ,sum(momentum_flux(i,j,k,:,3)),volume_force(i,j,k,3,2),volume_momentum_temp(i,j,k,3)
		endif

	enddo
	enddo
	enddo

	deallocate(volume_momentum_temp)

end subroutine control_volume_momentum_evolution

!===================================================================================
! Momentum Flux over a surface of a bin

subroutine control_volume_momentum_flux
use module_record
implicit none

	integer				:: ixyz,jxyz,kxyz,heaviside,n,planeno
	integer		,dimension(1)	:: imaxloc
	integer		,dimension(3)	:: ibin1,ibin2,crossplane
	double precision		:: crosstime, rplane
	double precision,dimension(3)	:: ri, mbinsize,velvect,bintop1,binbot1,bintop2,binbot2

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
			imaxloc = maxloc(abs(crossplane))
			ixyz = imaxloc(1)	!Integer array of size 1 copied to integer

			!Calculate velocity at time of intersection
			!crosstime = (r(n,ixyz) - rplane)/v(n,ixyz)
			velvect(:) = v(n,:) !- a(n,:) * crosstime

			!Add momentum flux to the new bin surface count and take from the old
			momentum_flux(ibin1(1),ibin1(2),ibin1(3),ixyz+3*heaviside(dble(crossplane(ixyz))),:) = & 
				momentum_flux(ibin1(1),ibin1(2),ibin1(3),ixyz+3*heaviside(dble(crossplane(ixyz))),:) & 
				+ velvect(:)*abs(crossplane(ixyz))
			momentum_flux(ibin2(1),ibin2(2),ibin2(3),ixyz+3*heaviside(-dble(crossplane(ixyz))),:) = & 
				momentum_flux(ibin2(1),ibin2(2),ibin2(3),ixyz+3*heaviside(-dble(crossplane(ixyz))),:) &
				- velvect(:)*abs(crossplane(ixyz))

		endif

	enddo

end subroutine control_volume_momentum_flux

!===================================================================================
!Forces over the surface of a Volume

subroutine control_volume_forces(fij,ri,rj,molnoi,molnoj)
use module_record
implicit none

	integer				:: ixyz, heaviside, n, molnoi, molnoj
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

	integer				:: i,j,k,ixyz,n,molnoi,molnoj,tempi,heaviside
	integer				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	integer,dimension(3)		:: cbin, ibin, jbin
	double precision		:: binforce
	double precision,dimension(3)	:: ri,rj,rij,fij,fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
	double precision,dimension(3)	:: Fbinsize, bintop, binbot, temp

	!Calculate rij
	rij = ri - rj
	!Prevent Division vy zero
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
		fsurface(:) = fsurface(:) + 0.5d0*fij(:)*dble(onfacexb - onfacext)
		fsurface(:) = fsurface(:) + 0.5d0*fij(:)*dble(onfaceyb - onfaceyt)
		fsurface(:) = fsurface(:) + 0.5d0*fij(:)*dble(onfacezb - onfacezt)
		volume_force(cbin(1),cbin(2),cbin(3),:,1) = volume_force(cbin(1),cbin(2),cbin(3),:,1) + fsurface*delta_t

	enddo
	enddo
	enddo

end subroutine control_volume_stresses

!===================================================================================
! Record change in momentum flux

subroutine record_volume_momentum
use module_record
implicit none

	if (irank .eq. iroot) print*, 'recording control volume momentum @ iter', iter

end subroutine record_volume_momentum


!===================================================================================
! 			FULL DOMAIN VIRIAL STRESS CALCULATION
! Calculate pressure_tensor for entire domain- average over all xy, yz and xz 
! components performed when calculating the interactions of particles
!===================================================================================

subroutine pressure_tensor
use module_record
implicit none

	integer, save	:: sample_count, average_count

	sample_count = sample_count + 1
	call cumulative_pressure_tensor(sample_count)
	if (sample_count .eq. viscsample) then
		sample_count = 0
		Pxyzero = Pxy		!Update Pxy(0) value
		if (irank .eq. iroot) print*, 'recording whole domain virial pressure tensor @ iter', iter
		call virial_stress_io
		average_count = average_count+1
		if (average_count .eq. Nvisc_ave) then
			print*, 'Calculating viscosity from domain virial stress @ iter', iter
			call record_pressure_tensor
			average_count = 0
		endif
	endif

end subroutine pressure_tensor

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


!===================================================================================
!Add pressure_tensor to running total

subroutine cumulative_pressure_tensor(sample_count)
use module_record
implicit none

	integer					:: sample_count
	integer         			:: n, ixyz,jxyz
	double precision, dimension(3)		:: velvect

	Pxy = 0.d0
	Pxymol = 0.d0

	!Factor of 2 as every interaction calculated
	rfmol = rfmol / 2.d0

	do n = 1, np    ! Calculate pressure tensor for all molecules

		!Calculate velocity at time t (v is at t+0.5delta_t due to use of verlet algorithm)
		velvect(:) = v(n,:) - 0.5d0 * a(n,:) * delta_t

		do ixyz = 1,3
		do jxyz = 1,3
			Pxymol(n,ixyz,jxyz) = Pxymol(n,ixyz,jxyz)    &
					      + velvect(ixyz)*velvect(jxyz) &
					      + rfmol(n,ixyz,jxyz)
		enddo
		enddo

		!Sum of molecules to obtain pressure tensor for entire domain
		Pxy(:,:) = Pxy(:,:) + Pxymol(n,:,:)	
	enddo

	!Sum pressure tensor over all processors
	call globalSumVect(Pxy(:,1), nd)
	call globalSumVect(Pxy(:,2), nd)
	call globalSumVect(Pxy(:,3), nd)

	!Devide sum of stress by volume
	Pxy = Pxy / volume

	!Calculate correrlation between current value of Pxy and intial value of sample
	!Average 3 components of Pxycorrel to improve statistics
	if(sample_count .ne. 0) then
		Pxycorrel(sample_count) = Pxycorrel(sample_count) + Pxy(2,3)*Pxyzero(2,3)
		Pxycorrel(sample_count) = Pxycorrel(sample_count) + Pxy(1,3)*Pxyzero(1,3)
		Pxycorrel(sample_count) = Pxycorrel(sample_count) + Pxy(1,2)*Pxyzero(1,2)
	endif

	!Reset position force tensor before calculation
	rfmol = 0.d0	  

end subroutine cumulative_pressure_tensor

!===================================================================================
!Write out pressure_tensor total and reset binned pressure_tensor

subroutine record_pressure_tensor
use module_record
use physical_constants_MD
implicit none

	double precision	::	viscosity

	call intergrate_trap(Pxycorrel,tplot*delta_t,viscsample,viscosity)

	viscosity = (viscosity*volume)/(3.0*Nvisc_ave*inputtemperature)
	
	if (irank .eq. iroot) then
		write(2,'(a)') 'Correlation'
		write(2,'(200f10.5)') Pxycorrel
		write(2,'(a)') 'viscosity record'
		write(2,'(i8,f18.9)') iter, viscosity
	endif

	Pxycorrel = 0.d0	!Reset Pxycorrel to zero

end subroutine record_pressure_tensor

!===================================================================================
! 			VOLUME AVERAGE STRESS TENSOR
! Stress calculated for each cell by assinging to cells based on location of 
! molecules (kinetic) and fraction of pairwise interaction in cell (potential)
! First used in Hardy (1982) J. Chem. Phys. 76(1),1 and simplified in paper by   ?!?!?!LUKSKO (1987) NOT HARDY?!?!?
! J. Cormier, J.M. Rickman and T. J. Delph (2001) Journal of Applied Physics,
! Volume 89, No. 1 "Stress calculations in atomistic simulations of perfect 
! and imperfect solids" 
!===================================================================================

subroutine pressure_tensor_VA
use module_record
implicit none

	integer			:: ixyz, jxyz
	integer			:: ibin,jbin,kbin
	integer, save		:: sample_count, average_count
	double precision	:: binvolume

	sample_count = sample_count + 1
	call cumulative_pressure_tensor_VA(1,nbins(1),1,nbins(2),1,nbins(3),sample_count)
	!call cumulative_pressure_tensor_VA(5,5,5,5,5,5,sample_count)
	if (sample_count .eq. viscsample) then
		Pxybin = Pxybin / sample_count !Average over samples
		sample_count = 0
		!Pxyzero = Pxy		!Update Pxy(0) value

		Pxy(:,:) = 0.d0

		do ibin = 1, nbins(1)
		do jbin = 1, nbins(2)
		do kbin = 1, nbins(3)
			do ixyz = 1,3
			do jxyz = 1,3
				Pxy(ixyz,jxyz) = Pxy(ixyz,jxyz) + Pxybin(ibin,jbin,kbin,ixyz,jxyz)
			enddo
			enddo
		enddo
		enddo
		enddo

		!Domain pressure calculated using division by volume of cell
		Pxy = Pxy/volume
		print*, 'Average Direct Pressure', (Pxy(1,1) + Pxy(2,2) + Pxy(3,3))/3.d0
		binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
		Pxybin = Pxybin / binvolume

		if(irank .eq. iroot) print*, 'Recording Volume Averaged pressure tensor @ iter', iter
		call VA_stress_io

		average_count = average_count+1
		if (average_count .eq. Nvisc_ave) then
			if(irank .eq. iroot) print*, 'Recording Volume Averaged viscosity @ iter', iter
			call record_pressure_tensor
			average_count = 0
		endif

		!Reset Stress per Bin ready for next averaging period
		Pxybin = 0.d0

	endif

end subroutine pressure_tensor_VA


!===================================================================================
!Add pressure_tensor to running total

subroutine cumulative_pressure_tensor_VA(imin,imax,jmin,jmax,kmin,kmax,sample_count)
use module_record
use calculated_properties_MD
implicit none

	integer					:: imin, jmin, kmin, imax, jmax, kmax
	integer					:: sample_count

 	!rfbin = 0.d0; vvbin = 0.d0

	!Calculate Position (x) force for configurational part of stress tensor
	call  simulation_compute_rfbins(imin,imax+2,jmin,jmax+2,kmin,kmax+2)
	!Calculate velocity (x) velocity for kinetic part of stress tensor
	call  simulation_compute_kinetic(imin,imax,jmin,jmax,kmin,kmax)

	!Add results to cumulative total
	Pxybin(:,:,:,:,:) = Pxybin(:,:,:,:,:) + vvbin(:,:,:,:,:) + rfbin(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)/2.d0

	!NOTE: REBUILD AT (mod(iter+1,tplot) .eq. 0) WHEN RECORD AFTER FORCE CALCULATION
	!Pxybin(:,:,:,:,:) = Pxybin(:,:,:,:,:) + rfbin(1:nbins(1),1:nbins(1),1:nbins(1),:,:)/2.d0
	!Reset bin force tensor before next calculation
  	rfbin = 0.d0; vvbin = 0.d0

end subroutine cumulative_pressure_tensor_VA

!----------------------------------------------------------------------------------
!Compute kinetic part of stress tensor

subroutine simulation_compute_kinetic(imin,imax,jmin,jmax,kmin,kmax)
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
	do kcell=kmin+1, kmax+1
	do jcell=jmin+1, jmax+1
	do icell=imin+1, imax+1

		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list
		ibin = icell-1; jbin = jcell-1; kbin = kcell-1	!ASSUME Cell same size as bins

		do i = 1,cellnp                  !Step through each particle in list 
			molnoi = oldi%molno 	 !Number of molecule

			!Calculate velocity at time t (v is at t+0.5delta_t due to verlet algorithm)
			velvect(:) = v(molnoi,:) - 0.5 * a(molnoi,:) * delta_t

			do ixyz = 1,3
			do jxyz = 1,3
				vvbin(ibin,jbin,kbin,ixyz,jxyz) = vvbin(ibin,jbin,kbin,ixyz,jxyz)	&
								  + velvect(ixyz) * velvect(jxyz)
			enddo
			enddo

			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list

			!print'(i8,3f10.5,3i8)', molnoi, r(molnoi,:), ibin, jbin, kbin
		enddo
	enddo
	enddo
	enddo

	nullify(oldi)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required

end subroutine simulation_compute_kinetic
!===================================================================================
!Write out pressure_tensor total and reset binned pressure_tensor

subroutine record_pressure_tensor_VA
use module_record
use physical_constants_MD
implicit none

	double precision	::	viscosity

	call intergrate_trap(Pxycorrel,tplot*delta_t,viscsample,viscosity)

	viscosity = (viscosity*volume)/(3.0*Nvisc_ave*inputtemperature)
	
	if (irank .eq. iroot) then
		write(2,'(a)') 'Correlation'
		write(2,'(200f10.5)') Pxycorrel
		write(2,'(a)') 'viscosity record'
		write(2,'(i8,f18.9)') iter, viscosity
	endif

	Pxycorrel = 0.d0	!Reset Pxycorrel to zero

end subroutine record_pressure_tensor_VA

!====================================================================================
! Use an expression, similar to the virial per bin, to partition the interaction
! between two molecules to seperate bins. Partitioning is based on the share of the
! distance between two molecules in a given bin
! N.B. Assume no more than 1 bin between bins containing the molecules

subroutine pressure_tensor_forces_VA(ri,rj,rij,accijmag)
use module_record
implicit none

	integer						:: ixyz, jxyz,i,j,k,l,n
	integer						:: diff
	integer,dimension(3)				:: ibin, jbin, bindiff 
	integer,dimension(:,:)		   ,allocatable	:: interbin, interbindiff
	double precision,intent(in)            		:: accijmag    !Non directional component of acceleration
	double precision				:: magnitude
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

	select case(diff)
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

	!print'(a,3(f18.10,a), 7i8)', '; ; ; ; ', sum(rfbin(ibin(1),ibin(2),ibin(3),:,:)),';', sum(rfbin(jbin(1),jbin(2),jbin(3),:,:)), ';', accijmag, ';', diff, ibin, jbin

end subroutine pressure_tensor_forces_VA

!========================================================================
!Compute Volume Averaged stress using all cells including halos

subroutine simulation_compute_rfbins(imin, imax, jmin, jmax, kmin, kmax)
use module_compute_forces
implicit none

	integer                         :: i, j, ixyz   !Define dummy index
	integer				:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp, cellsperbin
	integer				:: molnoi, molnoj, memloc
	integer				:: minbin, maxbin, imin, jmin, kmin, imax, jmax, kmax
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	rfbin = 0.d0
	rijsum = 0.d0

	!Calculate bin to cell ratio
	cellsperbin = ceiling(ncells(1)/dble(nbins(1)))

	do kcell=(kmin-1)*cellsperbin+1, kmax*cellsperbin
	do jcell=(jmin-1)*cellsperbin+1, jmax*cellsperbin
	do icell=(imin-1)*cellsperbin+1, imax*cellsperbin
	
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp                  !Step through each particle in list 
			molnoi = oldi%molno 	 !Number of molecule
			ri = r(molnoi,:)         !Retrieve ri

			do kcellshift = -1,1
			do jcellshift = -1,1
			do icellshift = -1,1

				!Prevents out of range values in i
				if (icell+icellshift .lt. imin) cycle
				if (icell+icellshift .gt. imax) cycle
				!Prevents out of range values in j
				if (jcell+jcellshift .lt. jmin) cycle
				if (jcell+jcellshift .gt. jmax) cycle
				!Prevents out of range values in k
				if (kcell+kcellshift .lt. kmin) cycle
				if (kcell+kcellshift .gt. kmax) cycle

				oldj => cell%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
				adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

				do j = 1,adjacentcellnp          !Step through all j for each i

					molnoj = oldj%molno 	 !Number of molecule
					rj = r(molnoj,:)         !Retrieve rj

					currentj => oldj
					oldj => currentj%next    !Use pointer in datatype to obtain next item in list

					if(molnoi==molnoj) cycle !Check to prevent interaction with self

					rij2=0                   !Set rij^2 to zero
					rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j

					do ixyz=1,nd
						rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
					enddo

					if (rij2 < rcutoff2) then
						!Add current distance to rijsum for molecules i and j
						rijsum(molnoi,:) = rijsum(molnoi,:) + 0.5d0*rij(:)
						rijsum(molnoj,:) = rijsum(molnoj,:) + 0.5d0*rij(:)
						!Linear magnitude of acceleration for each molecule
						invrij2 = 1.d0/rij2                 !Invert value
						accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)
						call pressure_tensor_forces_VA(ri,rj,rij,accijmag)

					endif
				enddo
			enddo
			enddo
			enddo
			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo
	enddo
	enddo

	nullify(oldi)      	!Nullify as no longer required
	nullify(oldj)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

end subroutine simulation_compute_rfbins

!===================================================================================
!	PRESSURE TENSOR CALCULATED USING METHOD OF PLANES (MOP)
! See B.D.Todd, D.J.Evans and P.J.Daivis (1995) Phys Review E, Vol 52,2
! Method first presented by Lutzko (1988) J. AppL Phys, Vol. 64, No.3,1 1154	!!!NO LUTSKO FIRST VA USING FOURIER TRANSFORMS!!!
!===================================================================================
!Calculate pressure_tensor across each plane

subroutine pressure_tensor_MOP
use module_record
implicit none

	integer, save	:: sample_count, average_count

	sample_count = sample_count + 1
	call cumulative_pressure_tensor_MOP(sample_count)
	if (sample_count .eq. viscsample) then
		sample_count = 0
		average_count = average_count+1
		if (average_count .eq. Nvisc_ave) then
			print*, 'Recording Method of Planes pressure tensor @ iter', iter
			call record_pressure_tensor_MOP
			average_count = 0
		endif
	endif

end subroutine pressure_tensor_MOP

!====================================================================================

subroutine pressure_tensor_forces_MOP(pnxyz,ri,rj,rij,accijmag)
use module_record
implicit none

	integer,save			:: i,j,k
	integer				:: n, heaviside
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

		!if (n .eq. 40 .and. sign(1.d0,ri(2)-plane)-sign(1.d0,rj(2)-plane) .ne. 0) print'(i8,a,3f10.5)',iter, &
		!									 ' MOP ',ri(2),rj(2),plane 

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
!Kinetic Part of the Method of Planes Pressure Tensor

subroutine pressure_tensor_kinetic_MOP
use module_record
implicit none

	integer         			:: planeno, n, ixyz
	integer					:: crossplane, heaviside
	double precision			:: crosstime, rplane
	double precision			:: shift
	double precision, dimension(3)		:: velvect

	!Shift by half difference between value rounded down and actual value
	!to ensure same distance to top and bottom plane from domain edge
	shift = 0.5d0*(domain(2) - planespacing * (nplanes-1))

	!Add molecular velocities to the configuration stresses
	do n=1,np

		!Replace Signum function with this functions which gives a
		!check for plane crossing and the correct sign 
		crossplane = ceiling((r(n,2)+halfdomain(2)-shift)/planespacing) & 
			    -ceiling((r(n,2)-delta_t*v(n,2)+halfdomain(2)-shift)/planespacing)

		if (crossplane .ne. 0) then

			!Obtain nearest plane number by integer division 
			!and retrieve location of plane from array
			planeno = ceiling((r(n,2)+halfdomain(2)-shift)/planespacing)-heaviside(dble(crossplane))+1
			if (planeno .lt. 1) planeno = 1
			if (planeno .gt. nplanes) planeno = nplanes
			rplane = planes(planeno)

			!print'(3i8,3f10.5)',planeno,crossplane,heaviside(dble(crossplane)),r(n,2)-delta_t*v(n,2), rplane, r(n,2)

			!Calculate velocity at time of intersection
			crosstime = (r(n,2) - rplane)/v(n,2)
			velvect(:) = v(n,:) - a(n,:) * crosstime

			!print'(2i8,4f10.5)',n,crossplane,r(n,2)-rplane,v(n,2),crosstime/delta_t
			if (crosstime/delta_t .gt. 1.d0) stop "error in kinetic MOP"

			!Obtain stress for three components on y plane
			Pxy_plane(:,planeno) = Pxy_plane(:,planeno) + velvect(:)*crossplane

		endif

	enddo


end subroutine pressure_tensor_kinetic_MOP

!===================================================================================
!Add pressure_tensor to running total

subroutine cumulative_pressure_tensor_MOP(sample_count)
use module_record
implicit none

	integer					:: sample_count

	!Configuration part of stress tensor already included

	!Add Kinetic Component to stress tensor
	call pressure_tensor_kinetic_MOP


end subroutine cumulative_pressure_tensor_MOP

!===================================================================================
!Write out pressure_tensor total and reset binned pressure_tensor

subroutine record_pressure_tensor_MOP
use module_record
use physical_constants_MD
implicit none

	!Devide by number of samples taken
	Pxy_plane = Pxy_plane/(Nvisc_ave*viscsample)

	!Divide by area of domain and factor of 4 for interactions
	Pxy_plane = Pxy_plane/(4*domain(1)*domain(3))

	call MOP_stress_io

	Pxy_plane = 0.d0


end subroutine record_pressure_tensor_MOP
