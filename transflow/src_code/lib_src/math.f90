!=======================================================================
! Math routines
!
! Charles Pierce, August 1997
!
! funct betai(a, b, x)
! funct betacf(a, b, x)*
! funct gammaln(x)
!

!=======================================================================
! Incomplete beta function
! See Numerical Recipes, s.6.4, p.219
!
function betai(a, b, x)

	if (x < 0. .or. x > 1.) stop "Bad argument x in betai"
	if (x == 0. .or. x == 1.) then
		bt = 0.
	else
		bt = exp(gammaln(a+b) - gammaln(a) - gammaln(b) &
		         + a*log(x) + b*log(1.-x))
	end if
	if (x < (a+1.)/(a+b+2.)) then
		betai = bt*betacf(a, b, x)/a
		return
	else
		betai = 1. - bt*betacf(b, a, 1.-x)/b
		return
	end if

end

!=======================================================================
! Auxiliary continued fraction computation for beta function
! See Numerical Recipes, s.6.4, p.219
!
function betacf(a, b, x)
	parameter (MAXIT=100, EPS=3.e-7, FPMIN=1.e-30)

	qab = a + b
	qap = a + 1.
	qam = a - 1.
	c = 1.
	d = 1. - qab * x / qap

	if (abs(d) < FPMIN) d = FPMIN
	d = 1. / d
	h = d
	do m=1,MAXIT
		m2 = 2 * m
		aa = m*(b - m)*x / ((qam + m2)*(a + m2))
		d = 1. + aa*d
		if (abs(d) < FPMIN) d = FPMIN
		c = 1. + aa/c
		if (abs(c) < FPMIN) c = FPMIN
		d = 1./d
		h = h*d*c
		aa = -(a + m)*(qab + m) * x / ((a + m2)*(qap + m2))
		d = 1. + aa*d
		if (abs(d) < FPMIN) d = FPMIN
		c = 1. + aa/c
		if (abs(c) < FPMIN) c = FPMIN
		d = 1./d
		del = d*c
		h = h*del
		if (abs(del-1.) < EPS) exit
	end do

	if (m >= MAXIT) stop "betacf: a or b too big or MAXIT too small"

	betacf = h

	return
end

!=======================================================================
! Returns the value ln[Gamma(x)] for x > 0.
! See Numerical Recipes, s.6.1, p.206
!
function gammaln(x)
	real x
	real coef(6)
	data coef / 76.18009172947146, -86.50532032941677, &
	            24.01409824083091, -1.231739572450155, &
	            0.1208650973866179e-2, -0.5395239384953e-5 /
	data sqrt2pi / 2.5066282746310005 /

	tmp = x + 5.5
	tmp = (x + 0.5)*log(tmp) - tmp

	y = x
	series = 1.000000000190015
	do j=1,6
		y = y + 1.
		series = series + coef(j)/y
	end do

	gammaln = tmp + log(sqrt2pi*series / x)

	return
end

