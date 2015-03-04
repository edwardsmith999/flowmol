!=======================================================================
! Commonly-used, general-purpose routines
!
! Charles Pierce, March 1996
!
! funct iopen()
! funct iclose(iu)
!
! funct rmaxabs(f, n)
! funct rmaxval(f, n)
! funct rminval(f, n)
! funct rsum(f, n)
!
! sectionTitle(iunit, name)
! signalShell(sig)
!
! funct lowerIndex(x, xvalue, nx)
!

!=======================================================================
function iopen()
	common /iopenclose_/ index, iunit(128)
	integer, save :: icall = 1

	if (icall == 1) then
		index = 1
		icall = 0
	end if

	iunit(index) = 0
	do i=1,index
		if (iunit(i) == 0) exit
	end do
	if (i == index) then
		index = index + 1
		if (index >= 128) stop "iopen: maximum unit number exceeded"
	end if
	iunit(i) = 1
	iopen = i + 10

	return
end

function iclose(iu)
	common /iopenclose_/ index, iunit(128)

	iu = iu - 10
	if (iu > 0 .and. iu < index) then
		iunit(iu) = 0
		iclose = iu + 10
	else
		iclose = -1
	end if

	return
end

!=======================================================================
function rmaxabs(f, n)
	real f(n)
	fmax = 0.
	do i=1,n
		fmax = max(fmax,abs(f(i)))
	end do
	rmaxabs = fmax
	return
end

function rmaxval(f, n)
	real f(n)
	fmax = f(1)
	do i=2,n
		fmax = max(fmax,f(i))
	end do
	rmaxval = fmax
	return
end

function rminval(f, n)
	real f(n)
	fmin = f(1)
	do i=2,n
		fmin = min(fmin,f(i))
	end do
	rminval = fmin
	return
end

function rsum(f, n)
	real f(n)
	x = 0.
	do i=1,n
		x = x + f(i)
	end do
	rsum = x
	return
end

!=======================================================================
subroutine sectionTitle(iunit, name)
	character*(*) name
	character*60 titleBar
	titleBar = "------------------------------&
	           &------------------------------"

	titleBar(6:7+len(name)) = " "//name//" "
	write (iunit, 10) titleBar

10  format (/ 2x, a /)

	return
end

!=======================================================================
! Typical signals: done, inflow, checkpoint, zdata, resubmit
!
subroutine signalShell(sig)
	character*(*) sig

	iunit = iopen()
	open(iunit, file="signal")
	write (iunit, *) sig
	close(iclose(iunit))

	return
end

!=======================================================================
function lowerIndex(x, xvalue, nx)
	real x(nx)
	do i=2,nx-1
		if (xvalue < x(i)) exit
	end do
	lowerIndex = i-1
	return
end

