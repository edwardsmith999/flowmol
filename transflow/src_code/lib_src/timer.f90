!=======================================================================
! Provides a general facility for timing parts of a program
!
! Charles Pierce, April 1996
!
! timer_init()
! realTimer(ctag, ipart, istart)
! cpuTimer(ctag, ipart, istart)
! timer_getMean(ctag, ipart, t)
! timer_globalize()
! timer_indexForTag(ctag, ipart, index, iwrite)*
! timer_info(iunit)
!

module timer
	parameter (maxTag=16, maxPart=256)

	character*16 tag(maxTag)         ! name tags for timers
	integer ntag                     ! number of tags

	real    time_(maxTag*maxPart)
	real    tmean(maxTag*maxPart)
	real    t2mean(maxTag*maxPart)
	real    tmin(maxTag*maxPart)
	real    tmax(maxTag*maxPart)
	integer nsample(maxTag*maxPart)

end module

!=======================================================================
subroutine timer_init()
	use timer

	tag = ""
	ntag = 0

	time_ = 0.
	tmean = 0.
	t2mean = 0.
	tmin = 0.
	tmax = 0.
	nsample = 0

	return
end

subroutine realTimer(ctag, ipart, istart)
	use timer
	character*(*) ctag

	! Wall clock time measured using realClock()

	call timer_indexForTag(ctag, ipart, i, istart)
	if (i .ne. 0) then

		if (istart == 1) then
			time_(i) = realClock()

		else if (istart == 0) then
			time = realClock() - time_(i)

			tmean(i) = (nsample(i)*tmean(i) + time) / (nsample(i) + 1)
			t2mean(i) = (nsample(i)*t2mean(i) + time**2) / (nsample(i) + 1)
			tmin(i) = min(tmin(i), time)
			tmax(i) = max(tmax(i), time)
			if (nsample(i) == 0) tmin(i) = time
			nsample(i) = nsample(i) + 1

		end if

	else
		print *, "timer: name not found: ", ctag
		stop

	end if

	return
end

subroutine cpuTimer(ctag, ipart, istart)
	use timer
	character*(*) ctag

	! CPU usage measured using cpuClock()

	return
end

subroutine timer_getMean(ctag, ipart, t)
	use timer
	character*(*) ctag

	call timer_indexForTag(ctag, ipart, index, istart)
	t = tmean(index)

	return
end

subroutine timer_globalize()
	use timer
	real r(maxTag*maxPart)

	call globalMax(tmean, maxTag*maxPart)
	call globalMax(t2mean, maxTag*maxPart)
	call globalMax(tmin, maxTag*maxPart)
	call globalMax(tmax, maxTag*maxPart)

	r = real(nsample)
	call globalMax(r, maxTag*maxPart)
	nsample = nint(r)

	return
end

subroutine timer_indexForTag(ctag, ipart, index, iwrite)
	use timer
	character*(*) ctag

	if (ipart <= 0 .or. ipart > maxPart) &
		stop "timer: invalid ipart reference"

	do i=1,ntag
		if (tag(i) == ctag) then
			index = ((i-1)*maxPart) + ipart
			return
		end if
	end do

	if (iwrite == 1) then
		! Add new tag
		if (ntag < maxTag) then
			ntag = ntag + 1
			tag(ntag) = ctag
			index = ((ntag-1)*maxPart) + ipart
		else
			stop "timer: maxTag insufficient"
		end if
	else
		index = 0
	end if

	return
end

!=======================================================================
subroutine timer_info(iunit)
	use timer

	if (ntag > 0) then

		call sectionTitle(iunit, "Timer")

		write (iunit, 11)
		do n=1,ntag
			write (iunit, 12) tag(n)
			total = 0.
			do j=1,maxPart
				i = ((n-1)*maxPart) + j
				total = total + tmean(i)*nsample(i)
			end do

			do j=1,maxPart
				i = ((n-1)*maxPart) + j
				if (nsample(i) == 0) cycle
				ttime = tmean(i)*nsample(i)
				write (iunit, 13) j, nsample(i), &
				                  ttime, tmean(i), &
				                  sqrt(abs(t2mean(i) - tmean(i)**2)), &
				                  tmin(i), tmax(i), &
				                  100.*ttime/total
			end do
		end do

	end if

11  format (2x, "Name", "       N", &
	        "      Total", "       Avg.", "       Std.", &
	        "        Min", "        Max", "   %time" /)
12  format (2x, a16)
13  format (3x, i3, 1x, i7, 5(2x, es9.3), 1x, f7.2)

	return
end

