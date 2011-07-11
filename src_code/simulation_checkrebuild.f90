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
	integer, intent(out)   :: rebuild
	double precision       :: vmax
	double precision, save :: rmax

	rebuild = 0

	!Trigger rebuild if record to be taken on next timestep
	if (mod(iter+1,tplot) .eq. 0) then
		rebuild = 1
	endif

	!Evaluate velocity magnitude
	vmagnitude = 0.d0
	do n = 1, np    ! Loop over all molecules
	do ixyz = 1, nd    ! Loop over all dimensions
		vmagnitude(n) = vmagnitude(n) + v(n,ixyz)*v(n,ixyz)
	enddo
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

	!Set local processor rmax to zero as rebuild is global
	if(rebuild .eq. 1) then
		rmax = 0.d0
	endif

end subroutine simulation_checkrebuild
