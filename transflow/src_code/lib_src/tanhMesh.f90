!=======================================================================
! Hyperbolic-tangent mesh generation
!
! Charles Pierce, June 1995
!
! Uses: triDiagonal, atanh, acosh
!
! tanhMesh(ymatch, ratio, imatch, iyc, nseg, y, ny)
! tanhMeshSeg(yp, ratio, jp, iyc, np, y, ny)
! tanhMatch(y, ratio, x, n)
! tanhStretch(x, ratio, yc, iend, slope)
!

subroutine tanhMesh(ymatch, ratio, imatch, iyc, nseg, y, ny)
	real    ymatch(0:nseg), ratio(nseg)
	integer imatch(0:nseg), iyc(2)
	real    y(ny)
	real    xmatch(0:nseg)

	! Calculate xmatch
	call tanhMatch(ymatch, ratio, xmatch, nseg)

	! Determine imatch from xmatch
	imatch(0) = 1
	do i=1,nseg-1
		percent = xmatch(i) - xmatch(i-1)
		imatch(i) = imatch(i-1) + percent*(ny-1)
	end do
	imatch(nseg) = ny

	do i=1,nseg
		if (i == 1 .and. iyc(1) == 1) then
			yc = 0.
			iend = 2
		else if (i == nseg .and. iyc(2) == 1) then
			yc = 1.
			iend = 1
		else
			yc = 0.5
			iend = 1
		end if
		y1 = ymatch(i-1)
		y2 = ymatch(i)
		j1 = imatch(i-1)
		j2 = imatch(i)
		do j=j1,j2
			if (j1 == j2) cycle
			yy = float(j-j1)/(j2-j1)
			call tanhStretch (yy, ratio(i), yc, iend, slope)
			y(j) = y1 + (y2-y1)*yy
		end do
	end do

	return
end

subroutine tanhMeshSeg(yp, ratio, jp, iyc, np, y, ny)
	real    yp(0:np), ratio(np)
	integer jp(0:np), iyc(2)
	real y(ny)

	do n=1,np
		if (n == 1 .and. iyc(1) == 1) then
			yc = 0.
			iend = 2
		else if (n == np .and. iyc(2) == 1) then
			yc = 1.
			iend = 1
		else
			yc = 0.5
			iend = 1
		end if

		y1 = yp(n-1)
		y2 = yp(n)
		j1 = jp(n-1)
		j2 = jp(n)
		do j=j1,j2
			yy = float(j-j1)/(j2-j1)
			call tanhStretch (yy, ratio(n), yc, iend, slope)
			y(j) = y1 + (y2-y1)*yy
		end do

	end do

	return
end

subroutine tanhMatch(y, ratio, x, n)
	real y(0:n), ratio(n), x(0:n)
	! y and ratio must be ordered
	real A(0:n,3)

	! Enforce continuity of slopes
	! 1./slope(i) = ratio(i)*(x(i)-x(i-1))/(x(n)-x(0))
	!                       *(y(n)-y(0))/(y(i)-y(i-1))

	! Set up equations
	A(0,1) = 0.
	A(0,2) = 1.
	A(0,3) = 0.
	A(n,1) = 0.
	A(n,2) = 1.
	A(n,3) = 0.
	x(0) = 0.
	x(n) = 1.

	! Each equation has the form:
	! (x(i)-x(i-1))*ratio(i)/(x(n)-x(0))*(y(n)-y(0))/(y(i)-y(i-1)) = 
	! (x(i+1)-x(i))*ratio(i+1)/(x(n)-x(0))*(y(n)-y(0))/(y(i+1)-y(i))

	do i=1,n-1
		st1 = ratio(i)/(x(n)-x(0))*(y(n)-y(0))/(y(i)-y(i-1))
		st2 = ratio(i+1)/(x(n)-x(0))*(y(n)-y(0))/(y(i+1)-y(i))
		A(i,1) = -st1
		A(i,2) = st1 + st2
		A(i,3) = -st2
		x(i) = 0.
	end do

	call triDiagonal(A(0,1), A(0,2), A(0,3), x, n+1)

	return
end

subroutine tanhStretch(x, ratio, yc, iend, slope)
	! yc gives the y-location of the point-of-inflection
	! which is the point where the stretching is the most coarse
	! It will usually be set to (0., 0.5, or 1.0)
	! iend determines to which end of the domain ratio is to be applied
	if (iend == 1) then
		x1 = -acosh(sqrt(ratio))
		y1 = tanh(x1)
		yc1 = y1/(y1 - 1.)
		if (yc <= yc1) then
			print *, "yc <= yc1: ", yc, " <= ", yc1
			stop "yc out of range"
		end if
		y2 = y1*(1. - 1./yc)
		x2 = atanh(y2)
	else if (iend == 2) then
		x2 = acosh(sqrt(ratio))
		y2 = tanh(x2)
		yc2 = 1./(1. + y2)
		if (yc >= yc2) then
			print *, "yc >= yc2: ", yc, " >= ", yc2
			stop "yc out of range"
		end if
		y1 = y2*(yc/(yc - 1.))
		x1 = atanh(y1)
	else
		stop "invalid value for iend"
	end if
	x = (x2 - x1)*x + x1
	slope = (y2 - y1)/(x2 - x1)/cosh(x)**2
	x = tanh(x)
	x = (x - y1)/(y2 - y1)
	return
end

!=======================================================================
! triDiagonal
! Solves a tri-diagonal system of algebraic equations using gauss 
! elimination.  The banded matrix given by (a,b,c) is destroyed, 
! and the rhs vector r is replaced by the solution vector.
!
! Single system, variable a,b,c.
!
subroutine triDiagonal(a, b, c, r, n)
	real a(n), b(n), c(n), r(n)

	! Forward elimination
	do i=2,n
		const = a(i)/b(i-1)
		b(i) = b(i) - c(i-1)*const
		r(i) = r(i) - r(i-1)*const
	end do

	! Back-substitution
	r(n) = r(n)/b(n)
	do i=n-1,1,-1
		r(i) = (r(i) - c(i)*r(i+1))/b(i)
	end do

	return
end

function atanh(y)
	if (y < -1. .or. y > 1.) stop "atanh: argument out of range"
	x = 0.
	dx = 1.
	do while (abs(dx) > 2.*spacing(x))
		dx = (y - tanh(x))*(cosh(x)**2)
		x = x + dx
	end do
	atanh = x
	return
end

function acosh(y)
	if (y < 1.) stop "acosh: argument out of range"
	if (y == 1.) then
		acosh = 0.
		return
	end if
	x = 10.
	dx = 1.
	do while (abs(dx) > 2.*spacing(x))
		dx = (y - cosh(x))/sinh(x)
		x = x + dx
	end do
	acosh = abs(x)
	return
end

