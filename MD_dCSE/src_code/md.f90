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
#if USE_COUPLER
	use md_coupler_socket, only : socket_coupler_invoke,  &
								  socket_coupler_init
#endif
	implicit none

	call messenger_invoke                   !Initialises MPI

#if USE_COUPLER
	call socket_coupler_invoke				!INITIALISES COUPLER INTERCOMMUNICATOR
#endif

	call setup_command_arguments            !Process command line arguments specifying restart and input files
	call messenger_init						!Establish processor topology

	!Check to see if simulation is a restart of a previous simualtion
	if (restart) then
		call setup_restart_inputs           !Recover simulation inputs from file
#if USE_COUPLER
		call socket_coupler_init			!SETUP COUPLER AND EXCHANGE WITH CFD
#endif
		call setup_set_parameters           !Calculate parameters using input
		call setup_restart_microstate       !Recover position and velocities
	else
		call setup_inputs                   !Input simulation parameters
#if USE_COUPLER
 		call socket_coupler_init			!SETUP COUPLER AND EXCHANGE WITH CFD
#endif
		call setup_set_parameters           !Calculate parameters using input
		call setup_initialise_microstate    !Setup position and velocities
	endif

	call assign_to_cell                     !Assign molecules to cells
	call messenger_proc_topology            !Obtain topology of processors
	call messenger_updateborders(1)         !Update borders between processors
	call assign_to_neighbourlist		    !Build neighbourlist using cell list
	call setup_initial_record               !Setup print headers and output inital

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
	use interfaces
	use computational_constants_MD
	use physical_constants_MD
#if USE_COUPLER
	use md_coupler_socket, only : socket_apply_continuum_forces, &
								  average_and_send_MD_to_CFD,apply_continuum_forces_flekkoy
#endif
	implicit none
  
	integer :: rebuild                              !Flag set by check rebuild to determine if linklist rebuild require

	initialstep = initialstep + 1                   !Increment initial step by one 

	do iter = initialstep, Nsteps                   !Loop over specified output steps 

		select case(integration_algorithm)
		case(leap_frog_verlet)
			call md_advance_lfv                     !Advance simulation (leap-frog Verlet algorithm)
		case(velocity_verlet)
			call md_advance_vv                      !Advance simulation (velocity Verlet algorithm)
		case default
			call error_abort('Incorrect integration algorithm specification')
		end select

 	end do

contains

	!---------------------------------------------
	! 		Leapfrom integration routines
	subroutine md_advance_lfv
		implicit none
		
		call simulation_compute_forces         !Calculate forces on particles	
		call simulation_record                 !Evaluate & write properties

#if USE_COUPLER
		call apply_boundary_force              ! Apply boundary force to prevent molecules leaving domain     
		call socket_apply_continuum_forces     ! CFD=> MD Apply CFD based coupling forces on MD
		call average_and_send_MD_to_CFD(iter)  ! MD=>CFD Calculate averages of MD to pass to CFD
#else
		call apply_boundary_force
#endif

		call simulation_move_particles_lfv                  !Move particles as a result of forces
		call messenger_updateborders(0)                     !Update borders between processors
		call simulation_checkrebuild(rebuild)
		if (rebuild .eq. 1) then
			call linklist_deallocateall             		!Deallocate all linklist components
			call sendmols                           		!Exchange particles between processors
			call sort_mols									!Reorder molecules to improve cache efficency
			call assign_to_cell                     		!Re-build linklist for domain cells
			call messenger_updateborders(1)         		!Update borders between processors
			call assign_to_neighbourlist		    		!Setup neighbourlist
		endif

		end subroutine md_advance_lfv

	!---------------------------------------------
	!Velocity Verlet integration routines
	subroutine md_advance_vv
		implicit none

		call simulation_move_particles_vv(1)        !Find r(t+dt) and v(t+dt/2)
		call messenger_updateborders(0)             !Update borders between processors

		call simulation_checkrebuild(rebuild)
		if (rebuild .eq. 1) then
			call linklist_deallocateall             !Deallocate all linklist components
			call sendmols                           !Exchange particles between processors
			call sort_mols							!Reorder molecules to improve cache efficency
			call assign_to_cell                     !Re-build linklist for domain cells
			call messenger_updateborders(1)         !Update borders between processors
			call assign_to_neighbourlist		    !Setup neighbourlist
		endif

		call simulation_compute_forces              !Calculate forces on all particles
		call simulation_move_particles_vv(2)        !Find v(t+dt)
		call simulation_record                      !Evaluate and write properties 

	end subroutine md_advance_vv

end subroutine simulation_MD

!=============================================================================
!
!                                FINISH
! Record final outputs and deallocates all arrays, linklists and pointers
!
!-----------------------------------------------------------------------------

subroutine finish_MD
	implicit none

	call messenger_syncall          !Synchronizes all processors using a barrier
	call finish_final_record        !Write summary of simulation and close output files
	call finish_clear_all           !Clear all arrays ready for next simulation
	call messenger_free             !Terminates MPI

end subroutine finish_MD

