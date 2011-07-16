!=============================================================================
!                                   inputs
! Input values used to set up the simulation such as number of dimensions and
! number of particles
!
!-----------------------------------------------------------------------------

module module_inputs

	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
        save
       
end module module_inputs
!----------------------------------------------------------------------------------

subroutine setup_inputs
use module_inputs
implicit none
  
	integer :: k, n

	call random_seed
	call random_seed(size=n)
	allocate(seed(n))

	open(11,file=trim(file_dir)//'input')

	!Input physical co-efficients
	read(11,* ) density          !Density of system
	read(11,* ) rcutoff          !Cut off distance for particle interaction
	rcutoff2= rcutoff**2         !Useful definition to save computational time
	read(11,* ) inputtemperature !Define initial temperature
	read(11,* ) initialnunits(1) !x dimension split into number of cells
	read(11,* ) initialnunits(2) !y dimension box split into number of cells

	if (nd .eq. 3) then
	read(11,* ) initialnunits(3) !z dimension box split into number of cells
	else
	read(11,* ) k		     !Read into dummy variable as value not used
	endif
	!Input computational co-efficients
	read(11,* ) Nsteps           !Number of computational steps
	read(11,* ) delta_t          !Size of time step
	read(11,* ) tplot            !Frequency at which to record results
	initialstep = 0   	     !Set initial step to one to start
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

	!Flags to determine if periodic boundaries required
	read(11,* ) periodic(1)
	read(11,* ) periodic(2)
	read(11,* ) periodic(3)

	close(11,status='keep')      !Close input file

	if (seed(1)==seed(2)) then
		call random_seed
		call random_seed(get=seed(1:n))
	endif
	
	!Assign different random number seed to each processor
	seed = seed * irank
	!Assign seed to random number generator
	call random_seed(put=seed(1:n))

	elapsedtime = 1.d0*delta_t*Nsteps !Set elapsed time to end of simualtion

end subroutine setup_inputs
