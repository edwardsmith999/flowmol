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
	use coupler
	implicit none

	call messenger_invoke	 		     	!Initialises MPI

	!Check to see if simulation is a restart of a previous simualtion
	call setup_command_arguments			!Establish command line arguments specifying restart and input files
	
	if (restart .eqv. .true.) then
		print*, 'Simulation restarted from file: ', initial_microstate_file
		call messenger_init			!Establish processor topology
		call setup_restart_inputs		!Recover simulation inputs from file
                call socket_coupler_init
		call setup_set_parameters		!Calculate parameters using input
		call setup_restart_microstate		!Recover position and velocities
	else
		call messenger_init   			!Establish processor topology
		call setup_inputs			!Input simulation parameters
                call socket_coupler_init
		call setup_set_parameters		!Calculate parameters using input
		call setup_initialise_microstate	!Setup position and velocities
	endif

	call assign_to_cell				!Assign molecules to cells
	call messenger_proc_topology			!Obtain topology of processors
	call messenger_updateborders(1)			!Update borders between processors
	call assign_to_neighbourlist_halfint		!Build neighbourlist using cell list
	call setup_initial_record			!Setup print headers and output inital

	! if coupled
	call coupler_create_map

contains 

	subroutine socket_coupler_init
		use coupler
		use computational_constants_MD, only : npx,npy,npz,delta_t,elapsedtime 
		use messenger, only                  :  myid, icoord
		implicit none

                integer nsteps_cfd
                real(kind(0.d0)) :: delta_t_CFD

                if (.not. coupler_is_active) return

		! if coupled calculation prepare exchange layout
		call coupler_md_init(npx,npy,npz,icoord,delta_t)
 

                ! fix NSTEPS for the coupled case
                delta_t_CFD = coupler_md_get_dt_cfd()
                nsteps_cfd  = coupler_get_nsteps()

                Nsteps = initialstep + nsteps_cfd * int(delta_t_cfd/delta_t)
                elapsedtime = elapsedtime + nsteps_cfd * delta_t_CFD

                if (myid .eq. 0) then 
			write(*,'(4(a,/),a,I7,/a,/a,E10.4,a,/,a,/a)') &
                                       "*********************************************************************", &
                                       "WARNING - this is a coupled run which resets the number              ", &
                                       " of extrasteps to:                                                   ", &
                                       "                                                                     ", &
                                       " nstep_cfd*(delta_t_cdf/delta_t_md) = ", nsteps_cfd * int(delta_t_cfd/delta_t), &
                                       "                                                                     ", & 
                                       " The elapsed time was changed accordingly to: ", elapsedtime, " s    ", &
 				       " The value of NSTEPS parameter form input file was discarded.        ", &
                                       "*********************************************************************"
                                              
               endif 

       end subroutine socket_coupler_init

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
  
	integer :: rebuild    				!Flag set by check rebuild to determine if linklist rebuild require

	initialstep = initialstep + 1		!Increment initial step by one 

	do iter = initialstep, Nsteps		!Loop over specified output steps 
			
		call simulation_compute_forces	!Calculate forces on particles	
		call simulation_record		!Evaluate & write properties to file
		call mass_flux_averaging	!Average mass flux before movement of particles

		call simulation_apply_constraint_forces				!Apply force to prevent molecules leaving domain
		call socket_coupler_apply_continuum_forces(iter)

		call simulation_move_particles					!Move particles as a result of forces
		call momentum_flux_averaging(vflux_outflag)			!!Average momnetum flux after movement of particles
		
                call socket_coupler_average(iter)

		call messenger_updateborders(0)			!Update borders between processors
		call simulation_checkrebuild(rebuild)		!Determine if neighbourlist rebuild required
		
		if(rebuild .eq. 1) then
			call linklist_deallocateall	   	!Deallocate all linklist components
			call sendmols			   	!Exchange particles between processors
			call assign_to_cell	  	   	!Re-build linklist for domain cells
			call messenger_updateborders(rebuild)	!Update borders between processors
			call assign_to_neighbourlist_halfint	!Setup neighbourlist
		endif
  
 	enddo

contains 

        subroutine socket_coupler_apply_continuum_forces(iter)
                use physical_constants_MD, only : np
                use arrays_MD, only :r,v,a
                use coupler
                implicit none
                
                integer, intent(in) :: iter

                integer :: iter_average, Naverage
                real(kind(0.d0)) :: delta_t_CFD
                logical, save :: first_time=.true.
                save Naverage

                if (.not. coupler_is_active) return 

                if (first_time) then
                        first_time = .false.
                        Naverage = int(coupler_md_get_dt_cfd()/delta_t)
                endif
                
                iter_average = mod(iter-1, Naverage)+1
                
                call coupler_md_apply_continuum_forces(np,r,v,a,iter_average)        !Apply force based on Nie,Chen an Robbins coupling
                
        end subroutine socket_coupler_apply_continuum_forces
       
        
        subroutine socket_coupler_average(iter)
                use physical_constants_MD, only : np
                use arrays_MD, only :r,v
                use coupler
                implicit none

                integer, intent(in) :: iter
             
                integer :: iter_cfd, iter_average, Naverage, save_period, average_period
                logical, save :: first_time=.true.
                save save_period, average_period, Naverage

                 if (.not. coupler_is_active) return

                 if (first_time) then 
                         first_time     = .false.
                         save_period 	= coupler_get_save_period()    ! period to save uc data (for debugging, testing)
                         average_period = coupler_get_average_period() ! collection interval in the average cycle
                         Naverage = int(coupler_md_get_dt_cfd() /delta_t)           ! number of steps in MD average cycle
                 endif

                 iter_average = mod(iter-1, Naverage)+1        ! current step
                 
                 iter_cfd     = (iter-initialstep)/Naverage +1 ! CFD corresponding step
                
                 if ( mod(iter_average,average_period) .eq. 0 ) then
                         call coupler_md_boundary_cell_average(np,r,v,send_data=.false.) ! accumlate velocities
                         if ( mod(iter_cfd,save_period) .eq. 0) then
                                 call coupler_uc_average_test(np,r,v,lwrite=.false.)
                         endif
                 endif
                 
                 ! Send accumulated results to CFD at the end of average cycle 
                 if (iter_average .eq. Naverage) then 
                         call coupler_md_boundary_cell_average(np,r,v,send_data=.true.)
                         if (mod(iter_cfd,save_period) .eq. 0 ) then
                                 call coupler_uc_average_test(np,r,v,lwrite=.true.)
                         endif
                 endif
                 
                
         end subroutine socket_coupler_average

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

