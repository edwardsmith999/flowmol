!----------------------------------------------------------------------------------
!
!                                Restart Microstate
! Set up position and velocity of each molecule (the microstate) based on the final
! state of a previous simulation
!
!---------------------------------------------------------------------------------

module module_restart_microstate

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD

end module module_restart_microstate
!----------------------------------------------------------------------------------

subroutine setup_restart_microstate
use module_restart_microstate
implicit none

	integer :: ixyz, n, nl
	double precision, dimension (nd) :: rc !Temporary variable

	nl = 0		!Reset local molecules count nl

	!read positions
	do n=1,globalnp
		do ixyz=1,nd
			read(13) rc(ixyz)	!Read position from file
		enddo

		!Checixyz if molecule is in domain of processor
		if(rc(1).lt.-halfdomain(1)+domain(1)*(iblock-1)) cycle
		if(rc(1).ge. halfdomain(1)+domain(1)*(iblock-1)) cycle
		if(rc(2).lt.-halfdomain(2)+domain(2)*(jblock-1)) cycle
		if(rc(2).ge. halfdomain(2)+domain(2)*(jblock-1)) cycle
		if(rc(3).lt.-halfdomain(3)+domain(3)*(kblock-1)) cycle
		if(rc(3).ge. halfdomain(3)+domain(3)*(kblock-1)) cycle

		!If molecules is in the domain then add to processor's total
		nl = nl + 1 !Local molecule count

		!Correct to local coordinates
		r(nl,1) = rc(1)-domain(1)*(iblock-1)
		r(nl,2) = rc(2)-domain(2)*(jblock-1)
		r(nl,3) = rc(3)-domain(3)*(kblock-1)

		!Read corresponding velocities
		do ixyz=1,nd
			read(13) v(nl,ixyz)
		enddo

		!print'(i8,6f10.5)', nl, r(nl,:), v(nl,:)
	enddo

	np = nl	!Correct local number of particles on processor

	close(13,status='delete') !Close final_state file

end subroutine setup_restart_microstate
