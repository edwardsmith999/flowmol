!=======================================================================
! Sine and Cosine transforms using Real FFT
!
! Charles Pierce, July 1996
! See _Numerical Recipes_ for explanation
!
! cosFFTM_staggered
!

!=======================================================================
! cosFFTM_staggered
! Multiple 1D staggered cosine transforms using FFT
!
subroutine cosFFTM_staggered(x, n, lot, isign)
	real x(0:n+1,lot)
	real wr(0:n/2-1), wi(0:n/2-1), wi1(0:n/2-1)

	PI = acos(-1.)
	do i=0,n/2-1
		wr(i) = cos(PI*i/n)
		wi(i) = sin(PI*i/n)
	end do
	do i=0,n/2-1
		wi1(i) = sin(PI*(i+0.5)/n)
	end do

	if (isign == 1) then

		! Preprocess data values
		do i=0,n/2-1
			do k=1,lot
				x1 = 0.5*(x(i,k) + x(n-i-1,k))
				x2 = wi1(i)*(x(i,k) - x(n-i-1,k))
				x(i,k) = x1 + x2
				x(n-i-1,k) = x1 - x2
			end do
		end do

		! Fourier transform
		call realFFTM(x, n, lot, +1)
		x = x * n

		! Postprocess transform values
		do i=0,n/2-1
			do k=1,lot
				x1 = wr(i)*x(2*i,k) + wi(i)*x(2*i+1,k)
				x2 = wi(i)*x(2*i,k) - wr(i)*x(2*i+1,k)
				x(2*i,k) = x1
				x(2*i+1,k) = x2
			end do
		end do

		do k=1,lot
			sum = 0.5*x(n,k)
			do i=n-1,1,-2
				sum1 = sum
				sum = sum + x(i,k)
				x(i,k) = sum1
			end do
		end do

		x = x * (2./n)
		x(0,:) = 0.5*x(0,:)

	else if (isign == -1) then

		x(0,:) = 2.*x(0,:)
		x = x * (n/2.)

		do k=1,lot
			tmp = x(n-1,k)
			do i=n-1,3,-2
				x(i,k) = x(i-2,k) - x(i,k)
			end do
			x(n,k) = 2.*tmp
		end do

		do i=0,n/2-1
			do k=1,lot
				x1 = wr(i)*x(2*i,k) + wi(i)*x(2*i+1,k)
				x2 = - wi(i)*x(2*i,k) + wr(i)*x(2*i+1,k)
				x(2*i,k) = x1
				x(2*i+1,k) = -x2
			end do
		end do

		! Inverse transform
		x = x * (1./n)
		call realFFTM(x, n, lot, -1)

		do i=0,n/2-1
			do k=1,lot
				x1 = 0.5*(x(i,k) + x(n-i-1,k))
				x2 = (0.25/wi1(i))*(x(i,k) - x(n-i-1,k))
				x(i,k) = x1 + x2
				x(n-i-1,k) = x1 - x2
			end do
		end do

	end if

	return
end

