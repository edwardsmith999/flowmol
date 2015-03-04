!=======================================================================
! Fast Fourier transform routines for COMPAQ platforms
!=======================================================================
! All these routines have extra padding to resemble SGI data structure
! It is not needed, but makes it easier to use the same data in the main code.
!=======================================================================
! "TIGHT" version of the code with the Fourier Transforms in pure COMPAQ format.
!	lib/transforms.compaq.tight.f90	:	real	:: x(0:n-1,:)
!						complex	:: x(0:n-1,:)	! real	:: x(0:2n-1,:)
!	lib/transforms.sgi.f90		:	real	:: x(0:n+1,:)
!						complex :: x(0:n  ,:)	! real	:: x(0:2n+1,:)
! Following versions will have all the routines running:
!						real	:: x(0:n+1,:)
!						complex	:: x(0:n  ,:)	! real	:: x(0:2n+1,:)
!=======================================================================
!
!
! realFFTM(x, n, lot, isign)
! realFFTM_A(x, n, lot, isign)
! SGI :: realFFT2DM_A(x, n1, n2, lot, isign)
! SGI :: realFFT3D(x, n1, n2, n3, isign)
! SGI :: realFFT3D_A(x, n1, n2, n3, isign)
! complexFFTM(x, n, lot, isign)
! SGI :: complexFFTM_A(x, n, lot, isign)
!

!=======================================================================
! realFFTM
! Multiple 1D Fast Fourier transforms
!
! Uses COMPAQ subroutines DFFT
!
! Structure of Fourier coefficient array is as follows:
! i = 0    :n/2:1   -> k = 0    :n/2 (real part)
! i = n/2+1:n-1:1   -> k = n/2-1:1   (imaginary part)
! Stores real first, then imag. parts
! from {0 ... n/2} forllowed by {n/2-1 .... 1}  wavenumbers (in that order)
!
!
!slow subroutine realFFTM(x, n, lot, isign)
!slow         real    :: x(0:n+1,lot)
!slow         integer :: STATUS
!slow 
!slow 	if (isign == +1) then
!slow 		! Forward-Fourier transform
!slow 		STATUS =  DFFT ('R','R','F',x(0,1),x(0,1),n,1)
!slow 		x = x * (1./n)
!slow 
!slow 	else if (isign == -1) then
!slow 		! Inverse-Fourier transform
!slow                 x = x * (1.*n)
!slow 		STATUS =  DFFT ('R','R','B',x(0,1),x(0,1),n,1)
!slow 
!slow 	end if
!slow 
!slow 	return
!slow end


!=======================================================================
! realFFTM
! Multiple 1D Fast Fourier transforms
!

subroutine realFFTM(x, n, lot, isign)
	integer FFTW_FORWARD,FFTW_BACKWARD
	parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

	integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
	parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

	integer FFTW_ESTIMATE,FFTW_MEASURE
	parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

	integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
	parameter (FFTW_OUT_OF_PLACE=0)
	parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

	integer FFTW_THREADSAFE
	parameter (FFTW_THREADSAFE=128)

	INTEGER          ::  n, lot, isign, i
	REAL             ::  x(0:n+1,lot), y(0:n+1,lot)
	INTEGER*8,SAVE   ::  plan
	INTEGER  ,SAVE   ::  nn = 0, iisign = 0

	!LAMBDA_LAHEY		ml_external rfftw_f77_one		!name is case-sensitive.
	!LAMBDA_LAHEY		ml_external rfftw_f77_create_plan	!name is case-sensitive.
	!LAMBDA_LAHEY		ml_external rfftw_f77_destroy_plan	!name is case-sensitive.

	if (n == 1) return

	! Initialize coefficients
	if ( nn.ne.n  .or.  iisign.ne.isign ) then
		if (iisign .ne.  0)  call rfftw_f77_destroy_plan(plan)
		if ( isign .eq. +1)  call rfftw_f77_create_plan(plan, n, FFTW_FORWARD , FFTW_MEASURE)
		if ( isign .eq. -1)  call rfftw_f77_create_plan(plan, n, FFTW_BACKWARD, FFTW_MEASURE)
		iisign = isign
		nn = n
	end if

	y = x
	do i=1,lot
		call rfftw_f77_one(plan, y(0,i), x(0,i))
	end do

	if (isign .eq. +1)   x = x * (1./n)

	return
end

!=======================================================================
! realFFTM_A
! Multiple 1D Fast Fourier transforms
! Lot on inner index
!
!FFTW??	subroutine realFFTM_A(x, n, lot, isign)
!FFTW??		REAL         ::   x(lot,0:n+1)
!FFTW??		RECORD /DXML_D_FFT_STRUCTURE/ FFT_STRUCT
!FFTW??		INTEGER,SAVE ::   nn = 0
!FFTW??		INTEGER      ::  STATUS, STRIDE
!FFTW??		CHARACTER*1 DIRECTION
!FFTW??	
!FFTW??		if (n == 1) return
!FFTW??	
!FFTW??	        STRIDE = lot
!FFTW??		! Initialize coefficients
!FFTW??		if (nn .ne. n) then
!FFTW??			STATUS = DFFT_INIT(n,FFT_STRUCT,.FALSE.)
!FFTW??			nn = n
!FFTW??		end if
!FFTW??	
!FFTW??		if (isign == +1) then
!FFTW??			! Forward-Fourier transform
!FFTW??			DIRECTION = 'F'
!FFTW??			do i=1,lot
!FFTW??	                   STATUS = DFFT_APPLY('R','R',DIRECTION,x(i,0),x(i,0),FFT_STRUCT,STRIDE)
!FFTW??			end do
!FFTW??	                x = x * (1./n)
!FFTW??	
!FFTW??		else if (isign == -1) then
!FFTW??			! Inverse-Fourier transform
!FFTW??			DIRECTION = 'B'
!FFTW??	                x = x * (1.*n)
!FFTW??			do i=1,lot
!FFTW??	                   STATUS = DFFT_APPLY('R','R',DIRECTION,x(i,0),x(i,0),FFT_STRUCT,STRIDE)
!FFTW??			end do
!FFTW??		end if
!FFTW??	
!FFTW??		return
!FFTW??	end
!FFTW??	

!=======================================================================
! realFFT2DM_A
! Multiple 2D Fast Fourier transforms
!
! Uses SGI subroutines DZFFTM1DU and ZDFFTM1DU
!
! Double application of 1-D real transform
! Structure of 2D Fourier coefficient array is as follows:
! i = 0:n1:2, 1:n1+1:2 -> k1 = 0:n1/2
! j = 0:n2:2, 1:n2+1:2 -> k2 = 0:n2/2
!
! NOTE: The values are not actual Fourier coefficients.
! This routine is suitable for Fourier decomposition in a poisson solver.
!
!SGI subroutine realFFT2DM_A(x, n1, n2, lot, isign)
!SGI 	real x(0:n1+1,lot,0:n2+1)
!SGI 	real, save, allocatable :: coeff1(:), coeff2(:)
!SGI 	integer, save :: nn1 = 0, nn2 = 0
!SGI 
!SGI 	! Initialize coefficients
!SGI 	if (nn1 .ne. n1 .or. nn2 .ne. n2) then
!SGI 		if (allocated(coeff1)) deallocate (coeff1, coeff2)
!SGI 		allocate (coeff1(n1+15), coeff2(n2+15))
!SGI 		call DZFFTM1DUI (n1, coeff1)
!SGI 		call DZFFTM1DUI (n2, coeff2)
!SGI 		nn1 = n1
!SGI 		nn2 = n2
!SGI 	end if
!SGI 
!SGI 	if (isign == 1) then
!SGI 		! Forward-Fourier transform
!SGI 		call DZFFTM1DU (-1, n1, lot*(n2+2), x, 1, n1+2, coeff1)
!SGI 		call DZFFTM1DU (-1, n2, (n1+2)*lot, x, (n1+2)*lot, 1, coeff2)
!SGI 		x = x * (1./n1/n2)
!SGI 
!SGI 	else if (isign == -1) then
!SGI 		! Inverse-Fourier transform
!SGI 		call ZDFFTM1DU (+1, n1, lot*(n2+2), x, 1, n1+2, coeff1)
!SGI 		call ZDFFTM1DU (+1, n2, (n1+2)*lot, x, (n1+2)*lot, 1, coeff2)
!SGI 
!SGI 	end if
!SGI 
!SGI 	return
!SGI end

!=======================================================================
! realFFT3D
! Single 3D Fast Fourier transform
!
! Uses SGI subroutines DZFFTM1DU, ZDFFTM1DU, and ZFFTM1D
!
! Triple application of 1-D real and complex transforms
! Structure of 3D Fourier coefficient array is as follows:
! i = 0:n1:2, 1:n1+1:2 -> k1 = 0:n1/2
! j = 0:n2-1           -> k2 = 0:n2/2, -n2/2:-1
! k = 0:n3-1           -> k3 = 0:n3/2, -n3/2:-1
! Array may be redimensioned as complex(0:n1/2,0:n2+1,0:n3+1)
!
!SGI subroutine realFFT3D(x, n1, n2, n3, isign)
!SGI 	real x(0:n1+1,0:n2+1,0:n3+1)
!SGI 	real, save, allocatable :: coeff1(:), coeff2(:), coeff3(:)
!SGI 	integer, save :: nn1 = 0, nn2 = 0, nn3 = 0
!SGI 
!SGI 	! Initialize coefficients
!SGI 	if (nn1 .ne. n1 .or. nn2 .ne. n2 .or. nn3 .ne. n3) then
!SGI 		if (allocated(coeff1)) deallocate (coeff1, coeff2, coeff3)
!SGI 		allocate (coeff1(n1+15), coeff2(2*(n2+15)), coeff3(2*(n3+15)))
!SGI 		call DZFFTM1DUI (n1, coeff1)
!SGI 		call ZFFTM1DI (n2, coeff2)
!SGI 		call ZFFTM1DI (n3, coeff3)
!SGI 		nn1 = n1
!SGI 		nn2 = n2
!SGI 		nn3 = n3
!SGI 	end if
!SGI 
!SGI 	if (isign == +1) then
!SGI 		! Forward-Fourier transform
!SGI 		call DZFFTM1DU (-1, n1, (n2+2)*(n3+2), x, 1, n1+2, coeff1)
!SGI 		call ZFFTM1D (-1, n3, (n1+2)*(n2+2)/2, x, (n1+2)*(n2+2)/2, 1, coeff3)
!SGI 		do k=0,n3+1
!SGI 			call ZFFTM1D (-1, n2, (n1+2)/2, x(0,0,k), (n1+2)/2, 1, coeff2)
!SGI 		end do
!SGI 		x = x * (1./n1/n2/n3)
!SGI 
!SGI 	else if (isign == -1) then
!SGI 		! Inverse-Fourier transform
!SGI 		do k=0,n3+1
!SGI 			call ZFFTM1D (+1, n2, (n1+2)/2, x(0,0,k), (n1+2)/2, 1, coeff2)
!SGI 		end do
!SGI 		call ZFFTM1D (+1, n3, (n1+2)*(n2+2)/2, x, (n1+2)*(n2+2)/2, 1, coeff3)
!SGI 		call ZDFFTM1DU (+1, n1, (n2+2)*(n3+2), x, 1, n1+2, coeff1)
!SGI 
!SGI 	end if
!SGI 
!SGI 	return
!SGI end

!=======================================================================
! realFFT3D_A
! Single 3D Fast Fourier transform
!
! Uses SGI subroutines DZFFTM1DU, ZDFFTM1DU
!
! Triple application of 1-D real transform
! Structure of 3D Fourier coefficient array is as follows:
! i = 0:n1:2, 1:n1+1:2 -> k1 = 0:n1/2
! j = 0:n2:2, 1:n2+1:2 -> k2 = 0:n2/2
! k = 0:n3:2, 1:n3+1:2 -> k3 = 0:n3/2
!
! NOTE: The values are not actual Fourier coefficients.
! This routine is suitable for Fourier decomposition in a poisson solver.
!
!SGI subroutine realFFT3D_A(x, n1, n2, n3, isign)
!SGI 	real x(0:n1+1,0:n2+1,0:n3+1)
!SGI 	real, save, allocatable :: coeff1(:), coeff2(:), coeff3(:)
!SGI 	integer, save :: nn1 = 0, nn2 = 0, nn3 = 0
!SGI 
!SGI 	! Initialize coefficients
!SGI 	if (nn1 .ne. n1 .or. nn2 .ne. n2 .or. nn3 .ne. n3) then
!SGI 		if (allocated(coeff1)) deallocate (coeff1, coeff2, coeff3)
!SGI 		allocate (coeff1(n1+15), coeff2(n2+15), coeff3(n3+15))
!SGI 		call DZFFTM1DUI (n1, coeff1)
!SGI 		call DZFFTM1DUI (n2, coeff2)
!SGI 		call DZFFTM1DUI (n3, coeff3)
!SGI 		nn1 = n1
!SGI 		nn2 = n2
!SGI 		nn3 = n3
!SGI 	end if
!SGI 
!SGI 	if (isign == +1) then
!SGI 		! Forward-Fourier transform
!SGI 		call DZFFTM1DU (-1, n1, (n2+2)*(n3+2), x, 1, n1+2, coeff1)
!SGI 		call DZFFTM1DU (-1, n3, (n1+2)*(n2+2), x, (n1+2)*(n2+2), 1, coeff3)
!SGI 		do k=0,n3+1
!SGI 			call DZFFTM1DU (-1, n2, n1+2, x(0,0,k), n1+2, 1, coeff2)
!SGI 		end do
!SGI 		x = x * (1./n1/n2/n3)
!SGI 
!SGI 	else if (isign == -1) then
!SGI 		! Inverse-Fourier transform
!SGI 		do k=0,n3+1
!SGI 			call ZDFFTM1DU (+1, n2, n1+2, x(0,0,k), n1+2, 1, coeff2)
!SGI 		end do
!SGI 		call ZDFFTM1DU (+1, n3, (n1+2)*(n2+2), x, (n1+2)*(n2+2), 1, coeff3)
!SGI 		call ZDFFTM1DU (+1, n1, (n2+2)*(n3+2), x, 1, n1+2, coeff1)
!SGI 
!SGI 	end if
!SGI 
!SGI 	return
!SGI end

!=======================================================================
! complexFFTM
! Multiple 1D Fast Fourier transforms
!
! Structure of Fourier coefficient array is as follows:
! i =   0   : n/2 -> k = 0:n/2
! i = n/2+1 : n-1 -> k = -(n/2-1):-1
!
!=======================================================================
!FFTW??	subroutine complexFFTM_RR(x, n, lot, isign)
!FFTW??		REAL	     ::   x(0:2*n+1,lot)
!FFTW??		RECORD /DXML_Z_FFT_STRUCTURE/ FFT_STRUCT
!FFTW??		INTEGER,SAVE ::   nn = 0
!FFTW??		INTEGER      ::  STATUS, STRIDE
!FFTW??		CHARACTER*1 DIRECTION
!FFTW??	
!FFTW??		if (n == 1) return
!FFTW??	
!FFTW??	        STRIDE = 1
!FFTW??		! Initialize coefficients
!FFTW??		if (nn .ne. n) then
!FFTW??			STATUS = ZFFT_INIT(n,FFT_STRUCT,.TRUE.)
!FFTW??			nn = n
!FFTW??		end if
!FFTW??	
!FFTW??		if (isign == +1) then
!FFTW??			! Forward-Fourier transform
!FFTW??			DIRECTION = 'F'
!FFTW??			do i=1,lot
!FFTW??	                   STATUS = ZFFT_APPLY('R','R',DIRECTION,x(0,i),x(n,i),x(0,i),x(n,i),FFT_STRUCT,STRIDE)
!FFTW??			end do
!FFTW??	                x = x * (1./n)
!FFTW??	
!FFTW??		else if (isign == -1) then
!FFTW??			! Inverse-Fourier transform
!FFTW??			DIRECTION = 'B'
!FFTW??	                x = x * (1.*n)
!FFTW??			do i=1,lot
!FFTW??	                   STATUS = ZFFT_APPLY('R','R',DIRECTION,x(0,i),x(n,i),x(0,i),x(n,i),FFT_STRUCT,STRIDE)
!FFTW??			end do
!FFTW??		end if
!FFTW??	
!FFTW??		return
!FFTW??	end
!FFTW??	
!=======================================================================
! complexFFTM
! Multiple 1D Fast Fourier transforms
!
! Structure of Fourier coefficient array is as follows:
! i =   0   : n/2 -> k = 0:n/2
! i = n/2+1 : n-1 -> k = -(n/2-1):-1
!
!=======================================================================
subroutine complexFFTM(x, n, lot, isign)
	integer FFTW_FORWARD,FFTW_BACKWARD
	parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

	integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
	parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

	integer FFTW_ESTIMATE,FFTW_MEASURE
	parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

	integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
	parameter (FFTW_OUT_OF_PLACE=0)
	parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

	integer FFTW_THREADSAFE
	parameter (FFTW_THREADSAFE=128)


	INTEGER      ::   n, lot, isign, i
	REAL	     ::   x(0:2*n+1,lot), y(0:2*n+1,lot)	! complex :: x(0:n,lot)
	INTEGER*8    ::   plan
	INTEGER,SAVE ::   nn = 0, iisign = 0

	!LAMBDA_LAHEY		ml_external fftw_f77_one		!name is case-sensitive.
	!LAMBDA_LAHEY		ml_external fftw_f77_create_plan	!name is case-sensitive.
	!LAMBDA_LAHEY		ml_external fftw_f77_destroy_plan	!name is case-sensitive.

	if (n == 1) return

	! Initialize coefficients
	if ( nn.ne.n  .or.  iisign.ne.isign ) then
		if (iisign .ne.  0)  call fftw_f77_destroy_plan(plan)
		if ( isign .eq. +1)  call fftw_f77_create_plan(plan, n, FFTW_FORWARD , FFTW_MEASURE)
		if ( isign .eq. -1)  call fftw_f77_create_plan(plan, n, FFTW_BACKWARD, FFTW_MEASURE)
		iisign = isign
		nn = n
	end if

	y = x
	do i=1,lot
		call fftw_f77_one(plan, y(0,i), x(0,i))
	end do

	if (isign .eq. +1)   x = x * (1./n)

	return
end


!=======================================================================
! complexFFTM
! Multiple 1D Fast Fourier transforms
!
! Structure of Fourier coefficient array is as follows:
! i =   0   : n/2 -> k = 0:n/2
! i = n/2+1 : n-1 -> k = -(n/2-1):-1
!
!=======================================================================
subroutine complexFFTM_CC(x, n, lot, isign)
	integer FFTW_FORWARD,FFTW_BACKWARD
	parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

	integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
	parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

	integer FFTW_ESTIMATE,FFTW_MEASURE
	parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

	integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
	parameter (FFTW_OUT_OF_PLACE=0)
	parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

	integer FFTW_THREADSAFE
	parameter (FFTW_THREADSAFE=128)


	INTEGER          ::  n, lot, isign, i
	COMPLEX	         ::  x(0:n,lot), y(0:n, lot)
	INTEGER*8,SAVE   ::  plan
	INTEGER  ,SAVE   ::  nn = 0, iisign = 0

	!LAMBDA_LAHEY		ml_external fftw_f77_one		!name is case-sensitive.
	!LAMBDA_LAHEY		ml_external fftw_f77_create_plan	!name is case-sensitive.
	!LAMBDA_LAHEY		ml_external fftw_f77_destroy_plan	!name is case-sensitive.

	if (n == 1) return

	! Initialize coefficients
	if ( nn.ne.n  .or.  iisign.ne.isign ) then
		if (iisign .ne.  0)  call fftw_f77_destroy_plan(plan)
		if ( isign .eq. +1)  call fftw_f77_create_plan(plan, n, FFTW_FORWARD , FFTW_MEASURE)
		if ( isign .eq. -1)  call fftw_f77_create_plan(plan, n, FFTW_BACKWARD, FFTW_MEASURE)
		iisign = isign
		nn = n
	end if

	y = x
	do i=1,lot
		call fftw_f77_one(plan, y(0,i), x(0,i))
	end do

	if (isign .eq. +1)   x = x * (1./n)

	return
end


!=======================================================================
! complexFFTM_A
! Multiple 1D Fast Fourier transforms
! Lot on inner index
!
! Uses SGI subroutine ZFFTM1D
!
! Structure of Fourier coefficient array is as follows:
! i = 0:n/2 -> k = 0:n/2
! i = n/2:n -> k = -n/2:-1
!
!SGI  subroutine complexFFTM_A(x, n, lot, isign)
!SGI  	real x(2*lot,0:n) ! complex x(lot,0:n)
!SGI  	real, save, allocatable :: coeff(:)
!SGI  	integer, save :: nn = 0
!SGI  
!SGI  	if (n == 1) return
!SGI  
!SGI  	! Initialize coefficients
!SGI  	if (nn .ne. n) then
!SGI  		if (allocated(coeff)) deallocate (coeff)
!SGI  		allocate (coeff(2*(n+15)))
!SGI  		call ZFFTM1DI (n, coeff)
!SGI  		nn = n
!SGI  	end if
!SGI  
!SGI  	if (isign == +1) then
!SGI  		! Forward-Fourier transform
!SGI  		call ZFFTM1D (-1, n, lot, x, lot, 1, coeff)
!SGI  		x = x * (1./n)
!SGI  
!SGI  	else if (isign == -1) then
!SGI  		! Inverse-Fourier transform
!SGI  		call ZFFTM1D (+1, n, lot, x, lot, 1, coeff)
!SGI  
!SGI  	end if
!SGI  
!SGI  	return
!SGI  end

