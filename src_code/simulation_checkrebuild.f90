!-----------------------------------------------------------------------------
!
!                            Check_Rebuild                               
! Check to see if neighbourlist needs to be rebuit and returns flag true
!
!-----------------------------------------------------------------------------
module module_checkrebuild

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD

end module module_checkrebuild
!----------------------------------------------------------------------------------

subroutine simulation_checkrebuild(rebuild)
use module_checkrebuild
implicit none
	
	integer                :: n, ixyz
	integer, save		   :: rb_count, total_rb=0
	integer, intent(out)   :: rebuild
	double precision       :: vmax
	double precision, save :: rmax = 0.d0

	double precision, save :: average_rb_count

	rebuild = 0

	!Trigger rebuild if record to be taken on next timestep
	if (mod(iter+1,tplot) .eq. 0) rebuild = 1

	!Trigger rebuild if FIXED REBUILD specified in inpuy
	!if (fixed_rebuild_flag .eq. 1) then
	!	if (mod(iter,fixed_rebuild) .eq. 0) rebuild = 1
	!endif

	!Evaluate velocity magnitude
	do n = 1, np    ! Loop over all molecules
		vmagnitude(n) = dot_product(v(:,n),v(:,n)) 
	enddo
	!Obtain maximum velocity
	vmax = maxval(vmagnitude)
	vmax = sqrt(vmax)

	!Calculate maximum possible displacement based on maximum velocity
	rmax = rmax + vmax * delta_t

	!Check if maximum displacment has exceeded extra distance
	if (rmax .gt. 0.5d0*delta_rneighbr) rebuild = 1

	!Ensure all processor flags are the same - if any rebuild
	!then all should rebuild
	call globalMaxInt(rebuild)

	!If rebuild is happening ...
	if(rebuild .eq. 1) then
		!Set local processor rmax to zero as rebuild is global
		rmax = 0.d0
		!Reset count and store average rebuild frequency
		average_rb_count = average_rb_count + dble(rb_count)
		rb_count = 0; total_rb = total_rb + 1
		if (irank .eq. iroot) then
			if (mod(total_rb,1000) .eq. 0) then	!Print rebuild stats every 1000 rebuilds
		   		!open(30,file=trim(prefix_dir)//'/results/rebuild_stats',access='append')
				!write(30,'(a,i7,a,f10.5,a)'), 'Up to iteration ', iter, & 
				!	' On average, rebuild happened every ', average_rb_count/total_rb, ' steps'
				!close(30,status='keep')
				print('(a,i7,a,f10.5,a)'), 'Up to iteration ', iter, & 
					' On average, rebuild happened every ', average_rb_count/total_rb, ' steps'
			endif
		endif
	else
		rb_count = rb_count + 1
	endif

end subroutine simulation_checkrebuild
