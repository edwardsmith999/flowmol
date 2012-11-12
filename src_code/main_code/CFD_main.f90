!=======================================================================
! Couette flow simulation
! 
! Tom Jelly (after Tamer Zaki (after Robert Jacobs))
!
program BLAYERCODE

	call main_init()          ! Initialize
	call main_restore()       ! Read restart data
	call simulation_run()     ! Run the simulation
	call main_save()          ! Save the results
	call main_free()          ! Clean up

	stop "Exited normally"
end
