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

	logical :: restart
	
	call messenger_invoke	 		     !Initialises MPI

! if coupled calculation prepare exchange layout 
        call exchange_grid_data

	!Check to see if simulation is a restart of a previous simualtion
	inquire(file=trim(file_dir)//'final_state', exist=restart)

	if (restart .eqv. .true.) then
		print*, 'Simulation restarted from "final_state" file'
		call messenger_init			!Establish processor topology
		call setup_restart_inputs_locate	!Recover simulation inputs from file
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
	call messenger_updateborders			!Update borders between processors
	call assign_to_halocell				!Assign halo molecules to cells
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
use coupler_md_communication, only : boundary_box_average, simulation_apply_continuum_forces
use messenger, only : myid
implicit none

	integer :: rebuild    !Flag set by check rebuild to determine if linklist rebuild required
        integer icfd, save_period

        save_period = 10 ! this value should come from CFD

        initialstep = initialstep + 1			   	!Increment initial step by one 

        do icfd = 1, nsteps_cfd + initise_steps   ! + initise_steps for the initialisation step of CFD
                print*, icfd, nsteps_cfd + initise_steps
                do iter = initialstep, Nsteps		 	        	!Loop over specified output steps
                        
                        call simulation_compute_forces_halfint	 	!Calculate forces on particles
                        
                        if (mod(iter,tplot) .eq. 0) then
                                call simulation_record		   	!Evaluate & write properties to file
                        endif
                        if (mflux_outflag .ne. 0) then
                                call mass_flux_averaging
                        endif
                        
!                        call simulation_apply_constraint_forces  	!Apply force to prevent molecules leaving domain
!                        call simulation_apply_continuum_forces(iter)	!Apply force based on Nie,Chen an Robbins coupling
                        call simulation_move_particles_tag		!Move particles as a result of forces
                        
                        if (vflux_outflag .ne. 0) then
                                call momentum_flux_averaging(vflux_outflag)
                        endif
                        
                        if ( mod(iter,average_period) == 0 ) then
                                 call boundary_box_average(send_data=.false.) ! accumlate velocities
                                 if ( mod(icfd-1,save_period) == 0 .and. icfd > 1) then
                                         call coupler_uc_average_test(.false.)
                                endif
                        endif

                        call messenger_updateborders		   	!Update borders between processors
                        call simulation_checkrebuild(rebuild)	   	!Determine if neighbourlist rebuild required
                        
                        if(rebuild .eq. 1 )  then
                                call linklist_deallocateall	   	!Deallocate all linklist components
                                call sendmols			   	!Exchange particles between processors
                                call assign_to_cell	  	   	!Re-build linklist for domain cells
                                call messenger_updateborders	   	!Update borders between processors
                                call assign_to_halocell		   	!Re-build linklist for halo cells
                                call assign_to_neighbourlist_halfint	!Setup neighbourlist
                        endif
                        
                enddo
                
! Average the boundary velocity and send the results to CFD
                call boundary_box_average(send_data=.true.)
                if ( mod(icfd-1,save_period) == 0 .and. icfd > 1) then
                        call coupler_uc_average_test(.true.)
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
	call messenger_free   			!Terminates MPI

end subroutine finish_MD


subroutine coupler_uc_average_test(lwrite)
        use physical_constants_MD, only : np
        use computational_constants_MD, only : halfdomain, domain, npx
        use messenger, only : myid
        use arrays_MD, only : r,v
        implicit none

        logical, intent(in) :: lwrite

        integer ib, kb, jb, ip, nlx,nlz, ierr
        real(kind=kind(0.d0)) rd(3), y0,ymin, ymax, dx, dy,dz
        real(kind=kind(0.d0)),allocatable, save :: uc_bin(:,:,:,:)
        logical,save :: firsttime=.true.
        character(len=64), save :: fn

         
        nlx = 11
        nlz = 2
        
        dx = 5.d0
        dy = 5.d0
        dz= 35.d0

        y0 = 15.d0
        
        ymin = y0 - 2.d0 * dy
        ymax = y0 + 2.d0 * dy 

        if(firsttime)then
                firsttime = .false.
                write(fn,'(a,i0,a,i0,a)')'md_vel_np',npx,'_r',myid,'.txt'
                write(0,*) 'file name ', fn
                write(0,*) 'domain ', domain

                allocate(uc_bin(2,nlz-1,nlx-1,4))
                uc_bin = 0.d0
                
                open(45, file=fn,position="rewind")
                write(45,*)'# dx,dy,dz ', dx,dy,dz
                close(45)

        endif

        if (lwrite) then
                call write_data
                return
        endif

        
 

!        write(0,*)'MD uc test', np, dy, ymin,ymax
        
        do ip = 1, np
! using global particle coordinates
                rd(:) = r(ip,:) + halfdomain(:) + (/ myid, 0, 0 /)*domain(:)
                
                if ( rd(2) < ymin .or. rd(2) > ymax ) then
! molecule outside the average layer
                        cycle
                endif
                
                ib = ceiling((rd(1) - myid*domain(1)) / dx) 
                kb = ceiling( rd(3)                   / dz)      
                jb = ceiling((rd(2) - ymin    )       / dy)
                
                if ( ib > 0 .and. ib    <  nlx  .and. &
                        kb > 0 .and. kb <   nlz  ) then 
!  this particle are in this ranks domain
                        uc_bin(1,kb,ib,jb) = uc_bin(1,kb,ib,jb) + v(ip,1)
                        uc_bin(2,kb,ib,jb) = uc_bin(2,kb,ib,jb) + 1.d0 
                else 
!                                       write(0,*) 'MD uc_average, outside domain rd', rd, ' bbox%bb ', bbo
                endif
        end do
        
! debug   
!                         do i = 1, size(uc_bin,dim=2)
!                          write(0, '(a,I4,64F7.1)') 'MD myid uc_bin(2,..',myid,uc_bin(2,1,:)
!                         enddo
        
        
        
!                        write(0,*) 'MD uc sum in boxes', myid
!                        do i = 1, size(uc_bin,dim=2)
!                                write(0, '(a,I4,64E11.4)') 'MD myid uc_bin(1,..',myid, uc_bin(1,1,:)
!                        enddo
! send it to CFD        

contains 

      subroutine write_data
               use mpi
               implicit none

               integer i

                open(45,file=fn,position="append")
                do i = 1,4
                 write(45, '(100(E12.4,1x))') uc_bin(:,:,:,i)
                enddo
                write(45, '(1x/1x)')
                close(45)
               
               uc_bin = 0.d0

       end subroutine write_data
end subroutine coupler_uc_average_test
