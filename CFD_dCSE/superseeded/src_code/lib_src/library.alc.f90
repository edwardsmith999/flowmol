!=======================================================================
! Commonly-used, general-purpose routines for SGI platforms
!
! Charles Pierce, August 1996
!
! funct random()
! randomSeed(iseed)
!
! funct realClock()
! funct timeRemaining()
!

!=======================================================================
function random()
	! Fortran 90 intrinsic
	call random_number(r)
	random = r
	return
end

subroutine randomSeed(iseed)
	integer, allocatable :: iseeds(:)
	! Fortran 90 intrinsic
	call random_seed(size=n)
	allocate (iseeds(n))
	iseeds = 0
	iseeds(1) = iseed
	call random_seed(put=iseeds)
	deallocate (iseeds)
	return
end

!=======================================================================
function realClock()
	include "mpif.h"
	real realClock
	realClock = MPI_Wtime()
	return
end

function timeRemaining()
	include "mpif.h"
	real, save :: wallTime0 = -9.
	real, save :: whours, wallTime

	! MPI_Wtime doesn't start at 0.
	if (wallTime0 == -9.) then
		wallTime0 = MPI_Wtime()
		call readFloat("whours",whours)
		wallTime = whours*3600.
	endif

	! Give 600 second safety factor (to read/write data files)
	timeRemaining = wallTime - (MPI_Wtime()-wallTime0) - 600.

	return
end

