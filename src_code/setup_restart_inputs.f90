!----------------------------------------------------------------------------------
!
!                                Restart Simulation inputs
! Set up inputs based on the final state of a previous simulation
!
!---------------------------------------------------------------------------------

module module_restart_inputs

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	
end module module_restart_inputs
!----------------------------------------------------------------------------------

subroutine setup_restart_inputs
use module_restart_inputs
implicit none

	integer :: k, n
	integer :: extrasteps
	integer :: checkint
	double precision :: checkdp

	!Allocate random number seed
	call random_seed
	call random_seed(size=n)
	allocate(seed(n))

	!Open file at first recorded value
	open(13,file=trim(file_dir)//'final_state', position='rewind', form='unformatted')

	read(13) np		     !Number of particles
	read(13) density             !Density of system
	read(13) rcutoff             !Cut off distance for particle interaction
	rcutoff2= rcutoff**2         !Useful definition to save computational time
	read(13) inputtemperature    !Define initial temperature
	read(13) initialnunits(1)    !x dimension split into number of cells
	read(13) initialnunits(2)    !y dimension box split into number of cells
	if (nd == 3) then
	read(13) initialnunits(3)    !z dimension box split into number of cells
	else
	read(13) k		     !read into dummy variable as value not used
	endif

	!read computational co-efficients
	read(13) Nsteps              !Number of computational steps
	read(13) delta_t             !Size of time step
	read(13) tplot               !Frequency at which to record results
	initialstep = Nsteps         !Set plot count to final plot of last +1
	read(13) seed(1)	     !Random number seed value 1
	read(13) seed(2)	     !Random number seed value 2
	read(13) elapsedtime	     !Total elapsed time of all restarted simulations

	!Check if values from input file are different and alert user 
	open(11,file=trim(file_dir)//'input')
	
	read(11,* ) checkint         !Number of dimensions
	if (checkint .ne. nd) print*, 'Discrepancy between no. dimension', &
	'in input & restart file - restart file will be used'
	read(11,* ) checkdp          !Density of system
	if (checkdp .ne. density) print*, 'Discrepancy between system density', &
	'in input & restart file - restart file will be used'
	read(11,* ) checkdp          !Cut off distance for particle interaction
	if (checkdp .ne. rcutoff) print*, 'Discrepancy between cut off radius', &
	'in input & restart file - restart file will be used'
	read(11,* ) checkdp	     !Define initial temperature
	if (checkdp .ne. inputtemperature) print*, 'Discrepancy between initial temperature', &
	'in input & restart file - restart file will be used'
	read(11,* ) checkint	     !x dimension split into number of cells
	if (checkint .ne. initialnunits(1)) print*, 'Discrepancy between x domain size', &
	'in input & restart file - restart file will be used'
	read(11,* ) checkint         !y dimension box split into number of cells
	if (checkint .ne. initialnunits(2)) print*, 'Discrepancy between y domain size', &
	'in input & restart file - restart file will be used'
	read(11,* ) checkint	     !z dimension box split into number of cells
	if (nd == 3) then
	if (checkint .ne. initialnunits(3)) print*, 'Discrepancy between z domain size', &
	'in input & restart file - restart file will be used'
	endif

	!Get number of extra steps, timestep and plot frequency from input file
	
	read(11,* ) extrasteps       !Number of computational steps
	read(11,* ) delta_t          !Size of time step
	read(11,* ) tplot            !Frequency at which to record results
	read(11,* ) delta_rneighbr   !Extra distance used for neighbour cell
	read(11,* ) seed(1)	     !Random number seed value 1
	read(11,* ) seed(2)	     !Random number seed value 2

	!Flag to determine if output is switched on
	read(11,* ) vmd_outflag
	read(11,* ) macro_outflag	
	read(11,* ) slice_outflag
	read(11,* ) Nslice_ave
	read(11,* ) pressure_outflag
	read(11,* ) viscsample
	read(11,* ) Nvisc_ave

	close(11,status='keep')      !Close input file

	elapsedtime = elapsedtime + delta_t*extrasteps !Set elapsed time to end of simualtion

	Nsteps = Nsteps + extrasteps !Establish final iteration step based on previous

end subroutine setup_restart_inputs
