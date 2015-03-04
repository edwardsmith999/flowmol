!=======================================================================
! Gauss quadrature integration rule
!
! Charles Pierce, March 1996
!
! gaussQuadrature(n, x, w)
! funct rLegendre(n, x)
!

!=======================================================================
! gaussQuadrature
! Integrates exactly a polynomial of order 2n-1 using n points
! on the standard interval [-1,1]
!
subroutine gaussQuadrature(n, x, w)
	parameter (itmax=200)
	real x(n), w(n)
	real xe(n/2 + 1)

	! Quadrature points are roots of Legendre polynomials
	! Search for roots, using Chebychev points as initial guess
	PI = acos(-1.)
	nroot = 0
	np = 2*n
	do while (nroot < n/2)
		nroot = 0
		np = 2*np
		x1 = -1.
		p1 = rLegendre(n, x1)
		do i=1,np
			x2 = -cos(0.5*(i-0.5)*PI/np)
			p2 = rLegendre(n, x2)
			if (p1*p2 < 0.) then
				! Sign changed; save interval
				nroot = nroot + 1
				x(nroot) = x1
				xe(nroot) = x2
			end if
			x1 = x2
			p1 = p2
		end do
	end do

	! Find roots using iterative technique
	! Use modified secant method to avoid getting stuck
	do i=1,nroot
		x1 = x(i)
		x2 = xe(i)
		p1 = rLegendre(n, x1)
		p2 = rLegendre(n, x2)
		do k=1,itmax
			if (abs(p1) .gt. abs(p2)) then
				x3 = x1 - p1/(p2-p1)*(x2-x1)
			else
				x3 = x2 - p2/(p2-p1)*(x2-x1)
			end if
			if (mod(k,4) .eq. 0) x3 = 0.5*(x1+x2)
			p3 = rLegendre(n, x3)
			if (p3*p1 <= 0.) then
				x2 = x3
				p2 = p3
			else
				x1 = x3
				p1 = p3
			end if
			if (x2-x1 .le. spacing(x3)) exit
		end do
		if (k >= itmax) stop "quadrature: convergence failed"
		x(i) = x3
		w(i) = 2.*(1.-x3**2)/(n*(rLegendre(n-1, x3) - x3*p3))**2
		x(n-i+1) = -x(i)
		w(n-i+1) = w(i)
	end do
	if (mod(n,2) .eq. 1) then
		x(n/2+1) = 0.
		w(n/2+1) = 2./(n*rLegendre(n-1, 0.))**2
	end if

	return
end

function rLegendre(n, x)
	y1 = 0.
	y = 1.
	do i=0,n-1
		y0 = y1
		y1 = y
		y = ((2*i+1)*x*y1 - i*y0)/(i+1)
	end do
	rLegendre = y
	return
end

