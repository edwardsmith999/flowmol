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
use module_initial_record
implicit none

	!Evaluate system properties on all processes
	call initial_macroscopic_properties

	!Calculate Control Volume starting state
	call initial_control_volume

	if (irank .eq. iroot) then

		!Display all parameters required to describe simulation
		print*, '=================Molecular Simulation Parameters======================'
		print*, 'Number of Dimensions: ', nd
		print*, 'Number of Particles: ', globalnp
		print*, 'Time Step - delta t: ',  delta_t
		print*, 'Total number of steps: ',  Nsteps - initialstep
		print*, 'Starting step of simulation:', initialstep
		print*, 'Generate output file every: ',  tplot, 'steps'
		print*, 'Density: ',              density
		print*, 'Initial Temperature: ',  inputtemperature
		print'(a,f19.15,a,f10.5)', ' Cut off distance:  ', rcutoff, &
				    '  Neighbour List Delta r:  ', delta_rneighbr
		print*, 'Initial unit size (FCC unit) in x,y and z'
		print*,	initialunitsize(1), initialunitsize(2), initialunitsize(3)
		print*, 'Domain in x,y and z: '
		print*, globaldomain(1), globaldomain(2), globaldomain(3)
		print'(a,i8,a,i8,a,i8,a)', 'Split into', ncells(1)*npx, & 
				       '  cells in x, ', ncells(2)*npy, &
				       '  cells in y, ', ncells(3)*npz, &
				       '  cells in z'
		print*, cellsidelength
		print*, 'Domain volume: ', volume
		print*, 'Number of processors in x,y and z: ', npx, npy, npz
		print*, 'So cells per Processor', nicellxl, nicellyl, nicellzl
		print*, 'Cells per Processor including Halos: ', ncellxl, ncellyl, ncellzl
		print*, 'Random seed used for initial velocity generation', seed
		!Viscosity is set to 1, wall velocity is maximum in domain
		print*, 'Molecular Reynolds number = ', (density * maxval(slidev(:,1)) * domain(2))/1.d0

		print*, '======================================================================='

		!open(12,file='results/statistics')

		!Set up print out table
		print '(8a)', &
		'Iteration; 	   VSum;        V^2Sum;        Temp;         KE;                 PE;                  TE;          Pressure;'

		!Print initial conditions for simulations at iteration 0
		print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f10.4)', &
		initialstep,';',vsum,';', v2sum,';', temperature,';', &
		kinenergy,';',potenergy,';',totenergy,';',pressure

	endif

end subroutine setup_initial_record

!----------------------------------------------------------------------------------
!Calculate Initial kinetic and potential energy as well as temperature and pressure

subroutine initial_macroscopic_properties
use module_initial_record
implicit none

	integer          :: n, k
	double precision :: vel

	vsum  = 0.d0      ! Reset all sums
	v2sum = 0.d0      ! Reset all sums

	do n = 1, np    ! Loop over all particles
	do k = 1, nd    ! Loop over all dimensions
		!Velocity component must be shifted back half a timestep to determine 
		!velocity of interest - required due to use of the leapfrog method
		vel = v(n,k)
		vsum = vsum + vel      !Add up all molecules' velocity components
		v2sum = v2sum + vel**2 !Add up all molecules' velocity squared components  
	enddo
	enddo

	!Calculate forces to obtain initial potential energies and virial
	call simulation_compute_forces

	!Obtain global sums for all parameters
	call globalSum(vsum)
	call globalSum(v2sum)
	call globalSum(virial)
	call globalSum(potenergysum)

	kinenergy   = (0.5d0 * v2sum) / globalnp 
	potenergy   = potenergysum /(2*globalnp) !N.B. extra 1/2 as all interactions calculated
	totenergy   = kinenergy + potenergy
	temperature = v2sum / (nd * globalnp)
	pressure    = (density/(globalnp*nd))*(v2sum+virial/2) !N.B. virial/2 as all interactions calculated
	
	!WARNING ABOUT THERMOSTATTED RESTARTS! WILL CONSERVE TEMPERATURE OF LAST ITERATION
	!For velocity rescaling thermostat
	initialenergy = (potenergysum+v2sum)/(np)

end subroutine initial_macroscopic_properties

!Add up initial control volume mass densities

subroutine initial_control_volume
use module_initial_record
implicit none

	integer				:: n
	integer, dimension(3)		:: ibin
	double precision, dimension(3)	:: mbinsize

	!Setup Initial Mass in Control Volume
	mbinsize(:) = domain(:) / nbins(:)
	do n = 1,np
		ibin(:) = ceiling((r(n,:)+halfdomain(:))/mbinsize(:)) + 1
		volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + 1
	enddo

	!Setup Initial Momentum in Control Volume
	volume_momentum = 0.d0
	do n = 1,np
		!Add up current volume momentum densities
		ibin(:) = ceiling((r(n,:)+halfdomain(:))/mbinsize(:)) + 1
		volume_momentum(ibin(1),ibin(2),ibin(3),:) = volume_momentum(ibin(1),ibin(2),ibin(3),:) + v(n,:)
		!print'(4i8,3f10.5)', n, ibin(:), v(n,:)
	enddo

end subroutine initial_control_volume


