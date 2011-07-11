!=======================================================================
! Fast Fourier transform routines for INTEL_MKL platforms
!=======================================================================
! All these routines resemble SGI data structure
!
!   -  The extra padding in array sizes is need
!       for [Real ----to---> Complex] transform
!
!   -  The extra padding in array sizes in NOT needed
!       for [Complex --to--> Complex] transform, but is 
!       included for consistency/convenience.
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
! Uses MKL subroutines FFT
!
! Structure of Fourier coefficient array is as follows:
! i = 0:n:2   -> k = 0:n/2 (real part)
! i = 1:n+1:2 -> k = 0:n/2 (imaginary part)
! Both real and imag. parts are stored
! from 0 ... n/2 wavenumbers (in that order)
! Array may be redimensioned as complex(0:n/2)
!
!FFT  subroutine realFFTM(x, n, lot, direc)
!FFT          real x(0:n+1,lot)
!FFT          real, save, allocatable :: wsave(:)
!FFT          integer, save :: nn = 0
!FFT          integer isign
!FFT  
!FFT          if (n == 1) return
!FFT  
!FFT          ! Initialize coefficients
!FFT          if (nn .ne. n) then
!FFT                  if (allocated(wsave)) deallocate (wsave)
!FFT                  allocate (wsave(2*n+4))
!FFT                  isign = 0
!FFT                  call DZFFT1D(x(0,1),n, 0, wsave)
!FFT                  nn = n
!FFT          end if
!FFT  
!FFT          if (direc == +1) then
!FFT                  ! Forward-Fourier transform
!FFT                  isign = -1
!FFT                  do i = 1,lot
!FFT                          call DZFFT1D(x(0,i),n, isign, wsave)
!FFT                  end do
!FFT                  x = x * (1./n)
!FFT  
!FFT          else if (direc == -1) then
!FFT                  ! Inverse-Fourier transform
!FFT                  isign = +1
!FFT                  x = x * (1.*n)
!FFT                  do i = 1,lot
!FFT                          call ZDFFT1D(x(0,i),n, isign, wsave)
!FFT                  end do
!FFT  
!FFT          end if
!FFT  
!FFT          return
!FFT  end

!=======================================================================
! realFFTM
! Multiple 1D Fast Fourier transforms
!
! Uses MKL subroutines DFT
!
! Structure of Fourier coefficient array is as follows:
! i = 0:n:2   -> k = 0:n/2 (real part)
! i = 1:n+1:2 -> k = 0:n/2 (imaginary part)
! Both real and imag. parts are stored
! from 0 ... n/2 wavenumbers (in that order)
! Array may be redimensioned as complex(0:n/2)
!
!
subroutine realFFTM(x, n, lot, direc)
	use MKL_DFTI
	INTEGER      ::   n, lot, direc
	REAL         ::   x((n+2)*lot)		! real :: x(0:n+1,lot)

	type(DFTI_DESCRIPTOR), SAVE, POINTER :: MyDescHandle
	INTEGER,SAVE          :: nn = 0
	INTEGER               :: DIMS, STATUS
	INTEGER, DIMENSION(2) :: STRIDE

	if (n == 1) return

        DIMS      = 1
        STRIDE   = (/0,1/)
	! Initialize coefficients
	if (nn .ne. n) then
		STATUS = DftiFreeDescriptor   (MyDescHandle)
		STATUS = DftiCreateDescriptor (MyDescHandle, DFTI_DOUBLE, DFTI_REAL, DIMS, n)
		STATUS = DftiSetValue         (MyDescHandle, DFTI_NUMBER_OF_TRANSFORMS, lot )
		STATUS = DftiSetValue         (MyDescHandle, DFTI_INPUT_DISTANCE , n+2 )
		STATUS = DftiSetValue         (MyDescHandle, DFTI_OUTPUT_DISTANCE, n+2 )
		STATUS = DftiSetValue         (MyDescHandle, DFTI_INPUT_STRIDES  , STRIDE)
		STATUS = DftiSetValue         (MyDescHandle, DFTI_OUTPUT_STRIDES , STRIDE)
		STATUS = DftiCommitDescriptor (MyDescHandle)
		nn = n
	end if

	if (direc == +1) then
		! Forward-Fourier transform
		STATUS = DftiComputeForward   (MyDescHandle , x)
                x = x * (1./n)

	else if (direc == -1) then
		! Inverse-Fourier transform
		STATUS = DftiComputeBackward  (MyDescHandle , x)
	end if

	return
end

!=======================================================================
! realFFTM_A
! Multiple 1D Fast Fourier transforms
! Lot on inner index
!
subroutine realFFTM_A(x, n, lot, direc)
	use MKL_DFTI
	INTEGER      ::   n, lot, direc
	REAL         ::   x(lot*(n+2))		! real :: x(lot,0:n+1)

	type(DFTI_DESCRIPTOR), SAVE, POINTER :: MyDescHandle
	INTEGER,SAVE          :: nn = 0
	INTEGER               :: DIMS, STATUS
	INTEGER, DIMENSION(2) :: STRIDE

	if (n == 1) return

        DIMS      = 1
        STRIDE   = (/0,lot/)
	! Initialize coefficients
	if (nn .ne. n) then
		STATUS = DftiFreeDescriptor   (MyDescHandle)
		STATUS = DftiCreateDescriptor (MyDescHandle, DFTI_DOUBLE, DFTI_REAL, DIMS, n)
		STATUS = DftiSetValue         (MyDescHandle, DFTI_NUMBER_OF_TRANSFORMS, lot )
		STATUS = DftiSetValue         (MyDescHandle, DFTI_INPUT_DISTANCE , 1 )
		STATUS = DftiSetValue         (MyDescHandle, DFTI_OUTPUT_DISTANCE, 1 )
		STATUS = DftiSetValue         (MyDescHandle, DFTI_INPUT_STRIDES  , STRIDE)
		STATUS = DftiSetValue         (MyDescHandle, DFTI_OUTPUT_STRIDES , STRIDE)
		STATUS = DftiCommitDescriptor (MyDescHandle)
		nn = n
	end if

	if (direc == +1) then
		! Forward-Fourier transform
		STATUS = DftiComputeForward   (MyDescHandle , x)
                x = x * (1./n)

	else if (direc == -1) then
		! Inverse-Fourier transform
		STATUS = DftiComputeBackward  (MyDescHandle , x)
	end if

	return
end


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
!           * FUNCTIONALITY NOT IMPLEMENTED IN MKL YET
!=======================================================================
!MKL_NOT_AVAILABLE	subroutine complexFFTM_RR(x, n, lot, direc)
!MKL_NOT_AVAILABLE	use MKL_DFTI
!MKL_NOT_AVAILABLE	INTEGER      ::   n, lot, direc
!MKL_NOT_AVAILABLE	REAL         ::   x(0:2*n+1,lot)
!MKL_NOT_AVAILABLE
!MKL_NOT_AVAILABLE	type(DFTI_DESCRIPTOR), SAVE, POINTER :: MyDescHandle
!MKL_NOT_AVAILABLE	INTEGER,SAVE          :: nn = 0
!MKL_NOT_AVAILABLE	INTEGER               :: DIMS, STATUS
!MKL_NOT_AVAILABLE	INTEGER, DIMENSION(2) :: STRIDE
!MKL_NOT_AVAILABLE
!MKL_NOT_AVAILABLE	if (n == 1) return
!MKL_NOT_AVAILABLE
!MKL_NOT_AVAILABLE        DIMS      = 1
!MKL_NOT_AVAILABLE        STRIDE   = (/0,1/)
!MKL_NOT_AVAILABLE	! Initialize coefficients
!MKL_NOT_AVAILABLE	if (nn .ne. n) then
!MKL_NOT_AVAILABLE		STATUS = DftiFreeDescriptor   (MyDescHandle)
!MKL_NOT_AVAILABLE		STATUS = DftiCreateDescriptor (MyDescHandle, DFTI_DOUBLE, DFTI_COMPLEX, DIMS, n)
!MKL_NOT_AVAILABLE		STATUS = DftiSetValue         (MyDescHandle, DFTI_NUMBER_OF_TRANSFORMS, lot )
!MKL_NOT_AVAILABLE		STATUS = DftiSetValue         (MyDescHandle, DFTI_COMPLEX_STORAGE, DFTI_REAL_REAL)
!MKL_NOT_AVAILABLE		STATUS = DftiSetValue         (MyDescHandle, DFTI_INPUT_DISTANCE , 2*n+2 )
!MKL_NOT_AVAILABLE		STATUS = DftiSetValue         (MyDescHandle, DFTI_OUTPUT_DISTANCE, 2*n+2 )
!MKL_NOT_AVAILABLE		STATUS = DftiSetValue         (MyDescHandle, DFTI_INPUT_STRIDES  , STRIDE)
!MKL_NOT_AVAILABLE		STATUS = DftiSetValue         (MyDescHandle, DFTI_OUTPUT_STRIDES , STRIDE)
!MKL_NOT_AVAILABLE		STATUS = DftiCommitDescriptor (MyDescHandle)
!MKL_NOT_AVAILABLE		nn = n
!MKL_NOT_AVAILABLE	end if
!MKL_NOT_AVAILABLE
!MKL_NOT_AVAILABLE	if (direc == +1) then
!MKL_NOT_AVAILABLE		! Forward-Fourier transform
!MKL_NOT_AVAILABLE		STATUS = DftiComputeForward   (MyDescHandle , x(0,1),x(n,1))
!MKL_NOT_AVAILABLE                x = x * (1./n)
!MKL_NOT_AVAILABLE
!MKL_NOT_AVAILABLE	else if (direc == -1) then
!MKL_NOT_AVAILABLE		! Inverse-Fourier transform
!MKL_NOT_AVAILABLE		STATUS = DftiComputeBackward  (MyDescHandle , x(0,1),x(n,1))
!MKL_NOT_AVAILABLE	end if
!MKL_NOT_AVAILABLE
!MKL_NOT_AVAILABLE	return
!MKL_NOT_AVAILABLE      end
!MKL_NOT_AVAILABLE

!=======================================================================
! complexFFTM
! Multiple 1D Fast Fourier transforms
!
! Structure of Fourier coefficient array is as follows:
! i =   0   : n/2 -> k = 0:n/2
! i = n/2+1 : n-1 -> k = -(n/2-1):-1
!
!=======================================================================
!    * FUNCTIONALITY IMPLEMENTED IN MKL but I am not sure it works
!    * It is best to use the routine with Complex numbers declarations,
!    instead of complex numbers saved in real data strings
!=======================================================================
!    *  I will rename (complexFFTM_CC) to (complexFFTM) till this
!    routine is safe to use and then we can separate them.
!=======================================================================
!NotBest	subroutine complexFFTM(x, n, lot, direc)
!NotBest		use MKL_DFTI
!NotBest		INTEGER      ::   n, lot, direc
!NotBest		REAL         ::   x((2*n+2)*lot)	! real    :: x(0:2*n+1,lot)
!NotBest							! complex :: x(0:n    ,lot)
!NotBest	
!NotBest		type(DFTI_DESCRIPTOR), SAVE, POINTER :: MyDescHandle
!NotBest		INTEGER,SAVE          :: nn = 0
!NotBest		INTEGER               :: DIMS, STATUS
!NotBest		INTEGER, DIMENSION(2) :: STRIDE
!NotBest	
!NotBest		if (n == 1) return
!NotBest	
!NotBest	        DIMS      = 1
!NotBest	        STRIDE   = (/0,1/)
!NotBest		! Initialize coefficients
!NotBest		if (nn .ne. n) then
!NotBest			STATUS = DftiFreeDescriptor   (MyDescHandle)
!NotBest			STATUS = DftiCreateDescriptor (MyDescHandle, DFTI_DOUBLE, DFTI_COMPLEX, DIMS, n)
!NotBest			STATUS = DftiSetValue         (MyDescHandle, DFTI_NUMBER_OF_TRANSFORMS, lot )
!NotBest			STATUS = DftiSetValue         (MyDescHandle, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_REAL)
!NotBest			STATUS = DftiSetValue         (MyDescHandle, DFTI_INPUT_DISTANCE , 2*n+2 )
!NotBest			STATUS = DftiSetValue         (MyDescHandle, DFTI_OUTPUT_DISTANCE, 2*n+2 )
!NotBest			STATUS = DftiSetValue         (MyDescHandle, DFTI_INPUT_STRIDES  , STRIDE)
!NotBest			STATUS = DftiSetValue         (MyDescHandle, DFTI_OUTPUT_STRIDES , STRIDE)
!NotBest			STATUS = DftiCommitDescriptor (MyDescHandle)
!NotBest			nn = n
!NotBest		end if
!NotBest	
!NotBest		if (direc == +1) then
!NotBest			! Forward-Fourier transform
!NotBest			STATUS = DftiComputeForward   (MyDescHandle , x)
!NotBest	                x = x * (1./n)
!NotBest	
!NotBest		else if (direc == -1) then
!NotBest			! Inverse-Fourier transform
!NotBest			STATUS = DftiComputeBackward  (MyDescHandle , x)
!NotBest		end if
!NotBest	
!NotBest		return
!NotBest	end
!NotBest	

!=======================================================================
! complexFFTM
! Multiple 1D Fast Fourier transforms
!
! Structure of Fourier coefficient array is as follows:
! i =   0   : n/2 -> k = 0:n/2
! i = n/2+1 : n-1 -> k = -(n/2-1):-1
!
!=======================================================================
!    * STANDARD / DEFAULT IMPLEMENTION IN MKL
!=======================================================================
!    *  I renamed (complexFFTM_CC) to (complexFFTM)
!    to use till as default since it is the best / default
!    MKL implementation
!=======================================================================
!RENAME	subroutine complexFFTM_CC(x, n, lot, direc)
	subroutine complexFFTM   (x, n, lot, direc)

	use MKL_DFTI
	INTEGER      ::   n, lot, direc
	COMPLEX      ::   x(n*lot)		! complex :: x(1:n, lot)

	type(DFTI_DESCRIPTOR), SAVE, POINTER :: MyDescHandle
	INTEGER,SAVE          :: nn = 0
	INTEGER               :: DIMS, STATUS
	INTEGER, DIMENSION(2) :: STRIDE

	if (n == 1) return

        DIMS      = 1
        STRIDE   = (/0,1/)
	! Initialize coefficients
	if (nn .ne. n) then
		STATUS = DftiFreeDescriptor   (MyDescHandle)
		STATUS = DftiCreateDescriptor (MyDescHandle, DFTI_DOUBLE, DFTI_COMPLEX, DIMS, n)
		STATUS = DftiSetValue         (MyDescHandle, DFTI_NUMBER_OF_TRANSFORMS, lot )
		STATUS = DftiSetValue         (MyDescHandle, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX) ! detault
		STATUS = DftiSetValue         (MyDescHandle, DFTI_INPUT_DISTANCE , n )
		STATUS = DftiSetValue         (MyDescHandle, DFTI_OUTPUT_DISTANCE, n )
		STATUS = DftiSetValue         (MyDescHandle, DFTI_INPUT_STRIDES  , STRIDE)
		STATUS = DftiSetValue         (MyDescHandle, DFTI_OUTPUT_STRIDES , STRIDE)
		STATUS = DftiCommitDescriptor (MyDescHandle)
		nn = n
	end if

	if (direc == +1) then
		! Forward-Fourier transform
		STATUS = DftiComputeForward   (MyDescHandle , x)
                x = x * (1./n)

	else if (direc == -1) then
		! Inverse-Fourier transform
		STATUS = DftiComputeBackward  (MyDescHandle , x)
	end if

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

