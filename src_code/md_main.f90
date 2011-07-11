
program md_main
implicit none

	double precision :: t1, t2 !Values used to establish simulation time

	call setup_MD   !Setup simulation 
	call cpu_time(t1)  !Get intial cpu time
	call simulation_MD !Calculate & write forces & movement of particles    
	call cpu_time(t2)  !Get final cpu time
	call finish_MD     !Write final values to file and clear all arrays

	print*, 'Time taken by code is; ', t2 - t1, '; seconds'

end program md_main

