!=============================================================================
!				   
! Routine which interface with the coupler to the CFD code
! Corresponding empty shell dummy routine used for uncoupled calculation
!
!  Lucian Anton, November 2011
!
! Gave up dummy subroutines for preprocessor flags
!
! LA, February 2012
!
!-----------------------------------------------------------------------------


module md_coupler_socket
	implicit none

#if USE_COUPLER 

contains

!=============================================================================
! Setup initial times based on coupled calculation
!-----------------------------------------------------------------------------
subroutine socket_coupler_init
	use coupler
	use computational_constants_MD, only : npx,npy,npz,delta_t,elapsedtime,Nsteps,initialstep, delta_t
	use messenger, only	      :  myid, icoord
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
				"*********************************************************************", 		&
 				"WARNING - this is a coupled run which resets the number	      ", 		&
				" of extrasteps to:						   ", 		&
				"								     ", 		&
				" nstep_cfd*(delta_t_cdf/delta_t_md) = ",nsteps_cfd*int(delta_t_cfd/delta_t), 	&
				"								     ", 		& 
				" The elapsed time was changed accordingly to: ", elapsedtime, " s    ", 		&
				" The value of NSTEPS parameter form input file was discarded.	", 		&
				"*********************************************************************"   
	endif 

end subroutine socket_coupler_init

!=============================================================================
!Apply coupling forces so MD => CFD
!-----------------------------------------------------------------------------
subroutine socket_coupler_apply_continuum_forces(iter)
	use computational_constants_MD, only : delta_t
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

	call coupler_md_apply_continuum_forces(np,r,v,a,iter_average)	!Apply force based on Nie,Chen an Robbins coupling
	
end subroutine socket_coupler_apply_continuum_forces

!=============================================================================
! Calculate averages of MD to pass to CFD
!-----------------------------------------------------------------------------
subroutine socket_coupler_average(iter)
	use computational_constants_MD, only : initialstep, delta_t
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
	    Naverage = int(coupler_md_get_dt_cfd() /delta_t)	   ! number of steps in MD average cycle
    endif

    iter_average = mod(iter-1, Naverage)+1	! current step
    
    iter_cfd     = (iter-initialstep)/Naverage +1 ! CFD corresponding step

    if ( mod(iter_average,average_period) .eq. 0 ) then
	    call coupler_md_boundary_cell_average(np,r,v,send_data=.false.) ! accumlate velocities

            ! collect uc data every save_perion cfd iteration but discard the first one which cfd uses for initialisation
	    if ( mod(iter_cfd-1,save_period) .eq. 0 .and. iter_cfd > 1) then
		    call coupler_uc_average_test(np,r,v,lwrite=.false.)
	    endif
    endif

	! Send accumulated results to CFD at the end of average cycle 
	if (iter_average .eq. Naverage) then 
	    call coupler_md_boundary_cell_average(np,r,v,send_data=.true.)
	    if (mod(iter_cfd-1,save_period) .eq. 0 .and. iter_cfd > 1) then
		    call coupler_uc_average_test(np,r,v,lwrite=.true.)
	    endif
	endif
	    	
end subroutine socket_coupler_average

#endif

end module md_coupler_socket
