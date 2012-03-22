!=======================================================================
! Piecewise-linear mesh generation
!
! Charles Pierce, June 1995
!
! linearMesh(yp, r, n, y, ny)
!
! yp  locations of segment matching points (returns xp)
! r   relative stretching to be applied at matching points
! n   number of matching points (including end points)
! y   mesh coordinates
! ny  number of mesh points
!
subroutine linearMesh(yp, r, n, y, ny)
	real yp(n), r(n)
	real y(ny)
	real xp(n)

	xp(1) = 0.
	do i=2,n
		a = r(i-1)
		b = (r(i)-r(i-1))/(yp(i)-yp(i-1))
		xp(i) = xp(i-1) + (1./b)*log((a+b*(yp(i)-yp(i-1)))/a)
	end do

	do j=1,ny
		x = (float(j-1)/(ny-1))*xp(n)
		i = lowerIndex(xp, x, n)
		a = r(i)
		b = (r(i+1)-r(i))/(yp(i+1)-yp(i))
		y(j) = (a/b)*(exp(b*(x-xp(i))) - 1.) + yp(i)
	end do

	xp = xp/xp(n)

	! Return xp-coordinates in yp
	yp = xp

	return
end

