program continuum
use computational_constants
use messenger
implicit none

	double precision :: t1, t2 !Values used to establish simulation time

        call messenger_invoke()

	call cpu_time(t1)  	!Get intial cpu time

	call setup_continuum    !Setup simulation
	continuum_initialstep = continuum_initialstep + 1  !Increment initial step by one (restart)

	do continuum_iter=continuum_initialstep,continuum_Nsteps   !Loop over all continuum steps
		!Run single step of continuum simulation
		call simulation_continuum	

		if (mod(continuum_iter,continuum_tplot).eq.0) then
			call continuum_record		   !Evaluate & write properties to file
			!continuum_tplot = continuum_tplot*2
		endif
	enddo

	call finish_continuum        !Write final values to file and clear all arrays

        call messenger_free

	call cpu_time(t2)  !Get final cpu time
	print*, 'Time taken by continuum code is; ', t2 - t1, '; seconds' 
 
end program continuum
