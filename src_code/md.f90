!=============================================================================
!
!                                  MAIN PROGRAM
!
! Molecular Dynamics Simulation Program
! Written by Edward Smith (unless otherwise commented), 10/08/09
! Updated 16/12/09
!
! program: md()
! subroutine: setup()
! subroutine: simulation()
! subroutine: finish()
!
!=============================================================================
!=============================================================================
!
!                                    SETUP  
!  
! Reads in inputs and initialises the simulation by calculating required parameters;
! establishing a starting molecular arrangement and spliting the domain into a
! series of cells each with a seperate linklist of local molecules.
! An initial print includes simulation information and output table headers.
!
!-----------------------------------------------------------------------------

subroutine setup_MD
use computational_constants_MD
implicit none

	logical :: restart
	
	call messenger_invoke	 		     !Initialises MPI

	!Check to see if simulation is a restart of a previous simualtion
	inquire(file='final_state', exist=restart)

	if (restart .eqv. .true.) then
		print*, 'Simulation restarted from "final_state" file'
		call messenger_init			!Establish processor topology
		call setup_restart_inputs		!Recover simulation inputs from file
		call setup_set_parameters		!Calculate parameters using input
		call setup_restart_microstate		!Recover position and velocities
	else
		call messenger_init			!Establish processor topology
		call setup_inputs			!Input simulation parameters
		call setup_set_parameters		!Calculate parameters using input
		call setup_initialise_microstate	!Setup position and velocities
	endif

	call assign_to_cell				!Assign molecules to cells
	call messenger_proc_topology			!Obtain topology of processors
	call messenger_updateborders			!Update borders between processors
	call assign_to_halocell				!Assign halo molecules to cells
	call assign_to_neighbourlist_halfint		!Build neighbourlist using cell list
	call setup_initial_record			!Setup print headers and output inital

end subroutine setup_MD

!=============================================================================
!
!                                SIMULATION
! Main part of the program, calculating forces and advancing the simulation
! Quantities of interest are calculated with all values recorded to file
! Each output step the program writes to file and then calculates for a number
! of iterations specified by freqency of plot: tplot
!
!-----------------------------------------------------------------------------

subroutine simulation_MD
use computational_constants_MD
implicit none
  
	integer :: rebuild    !Flag set by check rebuild to determine if linklist rebuild required
	
	initialstep = initialstep + 1			   	!Increment initial step by one

	do iter=initialstep,Nsteps			   	!Loop over specified output steps

		call simulation_compute_forces		 	!Calculate forces on particles

		if (mod(iter,tplot) .eq. 0) then
			call simulation_record		   	!Evaluate & write properties to file
		endif

		!call simulation_apply_constraint_forces   	!Apply force to prevent molecules leaving domain
		!call simulation_apply_continuum_forces	   	!Apply force based on Nie,Chen an Robbins coupling
		call simulation_move_particles			!Move particles as a result of forces
		call messenger_updateborders		   	!Update borders between processors
		call simulation_checkrebuild(rebuild)	   	!Determine if neighbourlist rebuild required

		if(rebuild .eq. 1) then
			call linklist_deallocateall	   	!Deallocate all linklist components
			call sendmols			   	!Exchange particles between processors
			call assign_to_cell	  	   	!Re-build linklist for domain cells
			call messenger_updateborders	   	!Update borders between processors
			call assign_to_halocell		   	!Re-build linklist for halo cells
			call assign_to_neighbourlist_halfint	!Setup neighbourlist
		endif

	enddo

end subroutine simulation_MD

!=============================================================================
!
!                                FINISH
! Record final outputs and deallocates all arrays, linklists and pointers
!
!-----------------------------------------------------------------------------

subroutine finish_MD
implicit none

	call messenger_syncall			!Synchronizes all processors using a barrier
	call finish_final_record		!Write summary of simulation and close output files
	call finish_clear_all                   !Clear all arrays ready for next simulation
	call messenger_free			!Terminates MPI

end subroutine finish_MD


