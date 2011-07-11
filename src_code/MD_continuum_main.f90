!=============================================================================
!
!                                  MAIN PROGRAM
!
! MD solver coupled to a Viscous Laminar Flow solver
! Written by Edward Smith (unless otherwise commented), 21/01/10
! Updated 16/01/10
!
! program: CFD-MD_main()
! subroutine: setup()
! subroutine: simulation()
! subroutine: finish()
!
!=============================================================================

program MD_continuum_main
implicit none

	double precision :: t1, t2 !Values used to establish simulation time

	call cpu_time(t1)  !Get intial cpu time
	call setup         !Setup simulation 
	call simulation    !Calculate & write    
	call finish        !Write final values to file and clear all arrays
	call cpu_time(t2)  !Get final cpu time
		
	print*, 'Time taken by code is; ', t2 - t1, '; seconds' 
 
end program MD_continuum_main

!=============================================================================
!
!                                    SETUP  
!  
!-----------------------------------------------------------------------------
	
subroutine setup
implicit none

	call setup_MD		!Setup molecular simulation to obtain BC for continuum
	call setup_coupling	!Setup coupling parameters
	call setup_continuum	!Setup continuum simulation
	
end subroutine setup

!=============================================================================
!
!                                SIMULATION
!
!-----------------------------------------------------------------------------

subroutine simulation
use computational_constants
use computational_constants_MD
implicit none

	integer :: t	      !Time t used to determine frequency of plots
	integer	:: extrasteps !Extra steps to add on to MD simulation
	
	t=0						   !Set plot counter to one
	continuum_initialstep = continuum_initialstep + 1  !Increment initial step by one (restart)
	extrasteps = Nsteps				   !Set extra steps for MD simulation 

	do continuum_iter=continuum_initialstep,continuum_Nsteps   !Loop over all continuum steps
		
		call simulation_continuum	!Run single step of continuum simulation
		call simulation_MD		!Run Nsteps (specified in input) of molecular simulation
		
		initialstep = initialstep + extrasteps - 1 !Increment initialstep of molecular simulation
		Nsteps = Nsteps + extrasteps		   !Increment final step of molecular simulation
		
		t=t+1					   !Plot counter increased by one	

		if (t.eq.continuum_tplot) then
			call continuum_record		   !Evaluate & write properties to file
			t = 0				   !Reset plot counter
		endif

	enddo

	!Nsteps adjusted to final number of Nsteps for VMD output formatting
	initialstep = initialstep - extrasteps + 1
	Nsteps = Nsteps - extrasteps	

end subroutine simulation

!=============================================================================
!
!                                FINISH
!
!-----------------------------------------------------------------------------

subroutine finish
implicit none

	call continuum_finish
	call messenger_syncall			!Synchronizes all processors using a barrier
	call finish_final_record		!Write summary of simulation and close output files
	call finish_clear_all                   !Clear all arrays ready for next simulation
	call messenger_free			!Terminates MPI
	
end subroutine finish
