
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
use coupler_md_setup, only : exchange_grid_data, create_map_cfd_md
implicit none

	call messenger_invoke	 		     	!Initialises MPI

	! if coupled calculation prepare exchange layout 
	call exchange_grid_data

	!Check to see if simulation is a restart of a previous simualtion
	call setup_command_arguments			!Establish command line arguments specifying restart and input files
	
	if (restart .eqv. .true.) then
		print*, 'Simulation restarted from file: ', initial_microstate_file
		call messenger_init			!Establish processor topology
		call setup_restart_inputs		!Recover simulation inputs from file
		call setup_set_parameters		!Calculate parameters using input
		call setup_restart_microstate		!Recover position and velocities
	else
		call messenger_init   			!Establish processor topology
		call setup_inputs_locate		!Input simulation parameters
		call setup_set_parameters		!Calculate parameters using input
		call setup_initialise_microstate	!Setup position and velocities
	endif

	call assign_to_cell				!Assign molecules to cells
	call messenger_proc_topology			!Obtain topology of processors
	call messenger_updateborders(1)			!Update borders between processors
	call assign_to_neighbourlist_halfint		!Build neighbourlist using cell list
	call setup_initial_record			!Setup print headers and output inital

	! if coupled
	call create_map_cfd_md

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
	use mpi
	use computational_constants_MD
	use physical_constants_MD, only : np
	use arrays_MD, only :r,v
	use coupler_md_global_data, only : use_coupling, nsteps_cfd => nsteps, average_period
	use coupler_md_communication, only : boundary_box_average, simulation_apply_continuum_forces, &
	                                     coupler_uc_average_test
	use messenger, only : myid
	implicit none
  
	integer :: rebuild    				!Flag set by check rebuild to determine if linklist rebuild required
	integer :: icfd, save_period

	if (.not. use_coupling) nsteps_cfd = 1
	save_period = 10 				!This value should come from CFD
	initialstep = initialstep + 1			!Increment initial step by one 

	do icfd = 1, nsteps_cfd + initise_steps   	!Initise_steps - number of Nsteps times to run code in initialisation
		do iter = initialstep, Nsteps		!Loop over specified output steps
				
			call simulation_compute_forces	!Calculate forces on particles
															
			if (mod(iter,tplot) .eq. 0) then
					call simulation_record		   	!Evaluate & write properties to file
			endif
			if (mflux_outflag .ne. 0) then
					call mass_flux_averaging
			endif
			!call simulation_apply_constraint_forces  	!Apply force to prevent molecules leaving domain
			!call simulation_apply_continuum_forces(iter)	!Apply force based on Nie,Chen an Robbins coupling
			call simulation_move_particles_tag				!Move particles as a result of forces
			
			if (vflux_outflag .ne. 0) then
					call momentum_flux_averaging(vflux_outflag)
			endif
			
			if ( mod(iter,average_period) == 0 ) then
					 call boundary_box_average(send_data=.false.) ! accumlate velocities
					 if ( mod(icfd-1,save_period) == 0 .and. icfd > 1) then
							 call coupler_uc_average_test(np,r,v,lwrite=.false.)
					endif
			endif
		
			call messenger_updateborders(0)			!Update borders between processors
			call simulation_checkrebuild(rebuild)		!Determine if neighbourlist rebuild required
			
			if(rebuild .eq. 1 .or. &
							(use_coupling .and. iter .eq. Nsteps)) then
				call linklist_deallocateall	   	!Deallocate all linklist components
				call sendmols			   	!Exchange particles between processors
				call assign_to_cell	  	   	!Re-build linklist for domain cells
				call messenger_updateborders(rebuild)	!Update borders between processors
				call assign_to_neighbourlist_halfint	!Setup neighbourlist
			endif
		
		enddo

		! The speed for CFD boundary must be computed somewhere here
 		if ( use_coupling ) then
			write(0,*) ' md: myid, icfd: ', myid,  icfd
 			call  boundary_box_average(send_data=.false.)
		endif

		! Average the boundary velocity and send the results to CFD
		call boundary_box_average(send_data=.true.)
		if ( mod(icfd-1,save_period) == 0 .and. icfd > 1) then
			call coupler_uc_average_test(np,r,v,lwrite=.true.)
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
	call finish_clear_all			!Clear all arrays ready for next simulation
	call messenger_free   			!Terminates MPI

end subroutine finish_MD

