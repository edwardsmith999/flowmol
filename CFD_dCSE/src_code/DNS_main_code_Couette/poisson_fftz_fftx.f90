!=======================================================================
! Controls the poisson solver
!
! poisson_init_dct()
! poisson_fftz_dctx()
!

module poisson_dct
	use data_export
	include "mpif.h"

	real, dimension(:), allocatable :: z_dct
	real*8                          :: array_fft(0:ngz),coeff_fft(ngzm+15)
	real*8                          :: array_dct(0:ngx)
        INTERFACE
                real*8 function realClock()
                end function realClock
        END INTERFACE

	double precision                :: p_left_halo (0:ngz,0:nly)
        double precision                :: p_right_halo(0:ngz,0:nly)

end module

!=======================================================================
subroutine poisson_init_dct()
	use poisson_dct

	real    dx2, dz2, pi, twopi, arg
	integer i, k, m

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Compute and store the modified wavenumber (vkz)
	!  for central differencing (used in FFT of Spanwise
	!  direction --> Page 40)
	!  (This replaces "subroutine waveno" in xiahua's code)
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	twopi=2.*2.*acos(0.)
	nzh=ngz/2+1
	dz2=(alz/float(ngz-1))**2
	c=twopi/float(ngzm)
	do k=1,ngzm
		if (k .lt. nzh) then
			m = k-1
		else
			m = k-ngzm-1
		end if
		arg = c*float(m)
		vkz(k)=2.*(1.-cos(arg))/dz2
	end do

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Compute and store the modified wavenumber (vkx)
	!  for central differencing (used in discrete cosine
	!  transform(DCT) in the streamwise direction)
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	twopi=2.*2.*acos(0.)
	nxh=ngx/2+1
	write(6,*) 'xpg(ngx,1)-xpg(1,1) = ',xpg(ngx,1)-xpg(1,1)
	dx2=((xpg(ngx,1)-xpg(1,1))/float(ngxm))**2
	c=twopi/float(ngxm)
	do i=1,ngxm
		if (i .lt. nxh) then
			m = i-1
		else
			m = i-ngxm-1
		end if
		arg = c*float(m)
		vkx(i)=2.*(1.-cos(arg))/dx2
	end do

	!cccccccccccccccccccccccccccccccccccccccc
	!  Initialize the multi-grid parameters
	!cccccccccccccccccccccccccccccccccccccccc
	call readInt("ngrid" , ngrid )

	!ccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Calculate the coefficients of poisson equation
	!ccccccccccccccccccccccccccccccccccccccccccccccccc
	call preprocess_dct()

	return
end


subroutine poisson_fftz_dctx(cpu_poisson_total)
	!cccccccccccccccccccccccccc
	! solve poisson equation
	!cccccccccccccccccccccccccc
	use poisson_dct
        real    :: cpu_poisson_total

        cpu_poisson_total=0.
	array_fft=0.
	coeff_fft=0.

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!    evaluate the source term using velocity gradient
	!    Source term is nothing but the Divergence of Ust  (D u*)
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	before = realClock() ! cpu_source
        ka = 1
        kb = ngzm
        ia = i1_T
        ib = i2_T
        ja = j1_T
        jb = j2_T

	qr(ka:kb, ia:ib, ja:jb)=ust(ka:kb, ia+1:ib+1, ja:jb)-ust(ka:kb, ia:ib, ja:jb)
	qr(ka:kb, ia:ib, ja:jb)= qr(ka:kb, ia:ib    , ja:jb)  &
				+vst(ka:kb, ia:ib, ja+1:jb+1)-vst(ka:kb, ia:ib, ja:jb)
	qr(ka:kb, ia:ib, ja:jb)= qr(ka:kb, ia:ib, ja:jb)  &
				+wst(ka+1:kb+1, ia:ib, ja:jb)-wst(ka:kb, ia:ib, ja:jb)
	
	qr(ka:kb, ia:ib, ja:jb)=qr(ka:kb, ia:ib, ja:jb)/dt

	after = realClock() ! cpu_source
	cpu_source= (after - before) - t_over

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Fourier transform the source term q(i,j,k) in z-direction 
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	before = realClock() ! cpu_fft1
	denom1 = 1./ float(ngzm)
	do j = j1_T,j2_T
	do i = i1_T,i2_T
		array_fft(0:ngzm-1) = qr(1:ngzm,i,j) 

		call realFFTM(array_fft,ngzm,1,+1)
		ka = ngzm/2 + 1
		kb = 2*(ngzm/2)
		qr(1:ka   , i, j) = array_fft(0:ngzm:2)
		qr(ka+1:kb, i, j) = array_fft(3:kb-1:2)

		!J------------------------------------------------------
		!J	call dzfft1dui(ngzm,coeff_fft)
		!J	call dzfft1du(-1,ngzm,array_fft,1,coeff_fft)
		!J	!-----------------------------------------------
		!J	!   First 'half' of (qr) has the real values.
		!J	!   Second 'half' has the imaginary component.   
		!J	!-----------------------------------------------
		!J	ka = ngzm/2 + 1
		!J	kb = 2*(ngzm/2)
		!J	qr(1:ka   , i, j) = array_fft(0:ngzm:2)*denom1
		!J	qr(ka+1:kb, i, j) = array_fft(3:kb-1:2)*denom1
		!J------------------------------------------------------
	end do
	end do

	after  = realClock() ! cpu_fft1
	cpu_fft1= (after - before) - t_over
	before = realClock() ! cpu_tranpose

	!------------------------------------
	!   Time for Transposing qr in qT
	!------------------------------------
	allocate(qT(0:ngx, 0:ngy, nlzm))
	call transpose_qr_qT(qr,qT)

	if (irank.eq.1) then
		tmp = sum(qT(1:ngxm,1:ngym,1));
		print*,'Total sum of RHS for (k=0) = ', iblock,jblock,tmp
	end if
	after  = realClock() ! cpu_transpose
	cpu_transpose= (after - before) - t_over
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Discrete cosine transform of qT in the x direction 
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	before = realClock() ! cpu_dct1
        ka = kTmin_2(irank)
        kb = kTmax_2(irank)
	do k = ka,kb
	kk = k-ka+1
	do j = 1,ngym
	        array_dct=0.0
		array_dct(0:ngxm-1) = qT(1:ngxm,j,kk)	
	
		call realFFTM(array_dct,ngxm,1,+1) ! See SinCosFFT.f90 in lib directory
		ia = ngxm/2 + 1
		ib = 2*(ngxm/2)
		qT(  1:ia, j,kk) = array_dct(0:ngxm:2)
		qT(ia+1:ib,j,kk) = array_dct(3:ib-1:2)
		
	end do
	end do
	after  = realClock() ! cpu_dct1
	cpu_dct1= (after - before) - t_over
	before = realClock() ! cpu_helmholtz

!------------------------------------------------------------------------------------------
!			START SOLUTION LOOP OVER WAVE NUMBERS (Kz)
!------------------------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled system in wave space 
!       real part first (k < nzm/2+1)
!cccccccccccccccccccccccccccccccccccccccccccccc
ka = kTmin_2(irank)
kb = kTmax_2(irank)

do 20 k = ka,kb
	kk = k-ka+1
	allocate(z_dct(memlenap_z))
       
	if(k.le.ngzm/2+1) then
		goto 8000
	else
		ki=k-ngzm/2
		goto 8001
	end if
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled system in wave space - real part 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 8000   nc=0
        !--------------------------
        !  Set up right-hand side
        !--------------------------
        do j=0,ngy
        do i=0,ngx
                ij=i+(ngxm+2)*j
                z_dct(irhs_mt(ngrid)+ij)=qT(i,j,kk)
                z_dct(ip_mt  (ngrid)+ij)=0.0
        end do
        end do
	nnx=ngxm
	nny=ngym

	!--------------------------------
	!  Solve the equation using tdma
	!--------------------------------
	mm=ngrid
	if(k==1) then
	call triDiagonal_y_k0(z_dct (ip_mt  (mm)),z_dct (irhs_mt(mm)),nnx,nny,             &
			      zap_mt(iaps_mt(mm)),zap_mt(iapn_mt(mm)),zap_mt(iapp_mt(mm)), &
			      zap_mt(ivp_mt (mm)),vkz(k))
	else
	call triDiagonal_y(   z_dct (ip_mt  (mm)),z_dct (irhs_mt(mm)),nnx,nny,k,           &
			      zap_mt(iaps_mt(mm)),zap_mt(iapn_mt(mm)),zap_mt(iapp_mt(mm)), &
			      zap_mt(ivp_mt (mm)),vkz(k))
	end if

	!--------------------------------
	!  Recover "phatr" from "z_dct"
	!--------------------------------
	do j=0,ngy
	do i=0,ngx
		ij=i+(ngxm+2)*j
		phatr(i,j,kk)=z_dct(ip_mt(ngrid)+ij)
	end do
	end do
	goto 21

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled system in wave space - imaginary part 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 8001   nc=0
	!--------------------------------
	!  Set up right-hand side
	!--------------------------------
	do j=0,ngy 
	do i=0,ngx
		ij=i+(ngxm+2)*j 
                z_dct(irhs_mt(ngrid)+ij)=qT(i,j,kk)
                z_dct(ip_mt  (ngrid)+ij)=0.0
	end do
	end do 
	nnx=ngxm
	nny=ngym
	!--------------------------------
	!  Solve the equation using tdma
	!--------------------------------
	mm=ngrid
	call triDiagonal_y(z_dct (ip_mt  (mm)),z_dct (irhs_mt(mm)),nnx,nny,ki,          &
			   zap_mt(iaps_mt(mm)),zap_mt(iapn_mt(mm)),zap_mt(iapp_mt(mm)), &
			   zap_mt(ivp_mt (mm)),vkz(ki))

	!-------------------------------
	!  Recover "phatr" from "z_dct"
	!-------------------------------
	do j=0,ngy
	do i=0,ngx
		ij=i+(ngxm+2)*j
		phatr(i,j,kk)=z_dct(ip_mt(ngrid)+ij)
	end do
	end do

 21   continue
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	deallocate(z_dct)
 20   continue
!------------------------------------------------------------------------------------------
!				FINISHED ALL WAVENUMBERS OF KZ
!------------------------------------------------------------------------------------------
	after  = realClock() ! cpu_helmholtz
	cpu_helmholtz= (after - before) - t_over
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Inverse discrete cosine transform of phatr in the x direction 
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	before = realClock() ! cpu_dct2
        ka = kTmin_2(irank)
        kb = kTmax_2(irank)
	do k = ka,kb
	kk = k-ka+1
	do j = 1,ngym
	        array_dct=0.0

		ia = ngxm/2 + 1
		ib = 2*(ngxm/2)
		array_dct(0:ngxm:2) = phatr( 1:ia,  j,kk)
		array_dct(3:ib-1:2) = phatr(ia+1:ib,j,kk)
		call realFFTM(array_dct,ngxm,1,-1) ! See SinCosFFT.f90 in lib directory
	        phatr(1:ngxm,j,kk)=array_dct(0:ngxm-1)

	end do
	end do
	after  = realClock() ! cpu_dct2
	cpu_dct2= (after - before) - t_over
	before = realClock() ! cpu_transpose

	!------------------------------------
	!  Time for Transposing phatr in p 
	!------------------------------------
	deallocate(qT)

	call transpose_phat_p(phatr,p)
	after  = realClock() ! cpu_transpose
	cpu_tranpose= cpu_transpose + (after - before) - t_over
	if (irank.eq.1) write(6,*) 'Transpose time =', cpu_transpose
	before = realClock() ! cpu_fft2

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! Fourier transform back to physical space, z-dir first
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do j = j1_T-1, j2_T+1
	do i = i1_T-1, i2_T+1
 		ka = ngzm/2 + 1
 		kb = 2*(ngzm/2)

		array_fft = 0.
		array_fft(0:ngzm:2) = p(1:ka   , i,j)
		array_fft(3:kb-1:2) = p(ka+1:kb, i,j)
		!for2D  array_fft(1)  = 0.
		!for2D  array_fft(kb+1) = 0.
	
		call realFFTM(array_fft,ngzm,1,-1)
		!J---------------------------------------------------
		!J	call dzfft1dui(ngzm,coeff_fft)
		!J	call zdfft1du(+1,ngzm,array_fft,1,coeff_fft)
		!J---------------------------------------------------
		p(1:ngzm,i,j) = array_fft(0:ngzm-1)*dt
	end do
	end do

	after = realClock() ! cpu_fft2
	cpu_fft2= (after - before) - t_over
	before= realClock() ! cpu_period

        !cccccccccccccccccccccccccccccccccccccccccccc
        ! TJ: make pressure periodic in x directions
        !cccccccccccccccccccccccccccccccccccccccccccc
        p_left_halo  = 0.
        p_right_halo = 0.

        if (iblock.eq.1 .or. iblock.eq.npx) then
	        if (npx.eq.1) then
       	        	p_right_halo(:,:) = p(:,  1,   :)
       	        	p_left_halo(:,:)  = p(:,nlxb-1,:)
        	else
                	if (iblock.eq.1) then
                        	p_right_halo(:,:) = p(:,  1,  :)
                        	call poisson_xperiodic(p_right_halo(0,0),p_left_halo(0,0),((ngzm+2)*(nly+1)))
                	else if (iblock.eq.npx) then
                        	p_left_halo(:,:) = p(:,nlxb-1,:)
                        	call poisson_xperiodic(p_left_halo(0,0),p_right_halo(0,0),((ngzm+2)*(nly+1)))
                	end if
        	end if
	end if

	if (iblock.eq.1) then
		p(:,  0, :) = p_left_halo(:,:)
	end if
	if (iblock.eq.npx) then
		p(:,nlxb,:) = p_right_halo(:,:)
	end if

	!ccccccccccccccccccccccccccccccccccccccccc
	! make pressure periodic in z directions
	!ccccccccccccccccccccccccccccccccccccccccc
	p(ngz,:,:)=p(1   ,:,:)
	p(0  ,:,:)=p(ngzm,:,:)
	after = realClock() ! cpu_period
	cpu_period= (after - before) - t_over
	cpu_poisson_total=cpu_source+cpu_fft1+cpu_transpose+cpu_dct1+cpu_helmholtz+cpu_dct2+cpu_fft2+cpu_period

	return
end
!--------------------------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled equation by TDMA in the y direction
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine triDiagonal_y(p_ta,rhs_ta,nnx,nny,k_ta, &
                       aps_ta,apn_ta,app_ta,vp_ta,vkz_ta)
	use data_export
        real, dimension(0:nnx+1,0:nny+1) :: rhs_ta,  p_ta, vp_ta
        real, dimension(  nnx  ,  nny  ) :: aps_ta,apn_ta,app_ta
        real, dimension(  nnx  ,  nny  ) ::   a_ta,  b_ta,  c_ta, d_ta
        integer :: nnx,nny,k_ta

        p_ta=0.0
        !-------------------------------------
        ! calculate the coefficients and rhs
        !-------------------------------------
        do j=1,nny
        do i=1,nnx

		!ccccccccccccccccccccccccccccccccccccccccccccc
		! SY & TJ: distinguish real and imag modified
		!          wavnumbers in x-direction
		!ccccccccccccccccccccccccccccccccccccccccccccc
	
		if(i.le.ngxm/2+1) then
			ii = i
		else
			ii = i - ngxm/2 
		end if	

		!ccccccccccccccccccccccccccccccccccccccccccccc
		! SY & TJ: solve using TDMA
		!ccccccccccccccccccccccccccccccccccccccccccccc
		
                d_ta(i,j)=-rhs_ta(i,j)
                c_ta(i,j)=-apn_ta(i,j)
                b_ta(i,j)=-app_ta(i,j)+vp_ta(i,j)*(vkz_ta+vkx(ii))
                a_ta(i,j)=-aps_ta(i,j)
	
        end do
        end do
        call TDMA(a_ta,b_ta,c_ta,d_ta,1,nny,1,nnx)
        !-------------------------------------
        ! update p_ta
        !-------------------------------------
        do j=1,nny
        do i=1,nnx
                p_ta(i,j)=d_ta(i,j)
        end do
        end do

	return
end
!--------------------------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled equation by TDMA in the y direction
! for the zeroth wavenumber of Kz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine triDiagonal_y_k0(p_ta,rhs_ta,nnx,nny, &
                          aps_ta,apn_ta,app_ta,vp_ta,vkz_ta)
	use data_export
        real, dimension(0:nnx+1,0:nny+1) :: rhs_ta,  p_ta, vp_ta
        real, dimension(  nnx  ,  nny  ) :: aps_ta,apn_ta,app_ta
        real, dimension(  nnx  ,  nny  ) ::   a_ta,  b_ta,  c_ta, d_ta
        integer :: nnx,nny

        p_ta=0.0
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! special treatment to remove singularity when i=k=1
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !-------------------------------------
        ! calculate the coefficients and rhs
        !-------------------------------------
        i=1
        do j=1,nny-1
                d_ta(i,j)=-rhs_ta(i,j)
                c_ta(i,j)=-apn_ta(i,j)
                b_ta(i,j)=-app_ta(i,j)+vp_ta(i,j)*(vkz_ta+vkx(i))
                a_ta(i,j)=-aps_ta(i,j)
        end do
        call TDMA(a_ta,b_ta,c_ta,d_ta,1,nny-1,1,1)
        !-------------------------------------
        ! update p_ta
        !-------------------------------------
        do j=1,nny-1
                p_ta(i,j)=d_ta(i,j)
        end do
        !ccccccccccccccccccccccc
        ! solve for i > 1
        !ccccccccccccccccccccccc
        !-------------------------------------
        ! calculate the coefficients and rhs
        !-------------------------------------
        do j=1,nny
        do i=2,nnx

		!ccccccccccccccccccccccccccccccccccccccccccccc
		! SY & TJ: distinguish real and imag modified
		!          wavnumbers in x-direction
		!ccccccccccccccccccccccccccccccccccccccccccccc
	
		if(i.le.ngxm/2+1) then
			ii = i
		else
			ii = i - ngxm/2 
		end if	

		!ccccccccccccccccccccccccccccccccccccccccccccc
		! SY & TJ: solve using TDMA
		!ccccccccccccccccccccccccccccccccccccccccccccc

                d_ta(i,j)=-rhs_ta(i,j)
                c_ta(i,j)=-apn_ta(i,j)
                b_ta(i,j)=-app_ta(i,j)+vp_ta(i,j)*(vkz_ta+vkx(ii))
                a_ta(i,j)=-aps_ta(i,j)

        end do
        end do
        call TDMA(a_ta,b_ta,c_ta,d_ta,1,nny,2,nnx)
        !-------------------------------------
        ! update p_ta
        !-------------------------------------
        do j=1,nny
        do i=2,nnx
                p_ta(i,j)=d_ta(i,j)
        end do
        end do

	return
end
!--------------------------------------------------------

!ccccccccccccccccccccccccccccccccccccccccc
! solve the equation by Thomas algorithm
!ccccccccccccccccccccccccccccccccccccccccc
subroutine TDMA(a,b,c,d,ns,nf,nis,nif)
	use data_export
        real, dimension(ngxm,ngym) :: a,b,c,d,f,g
        integer :: ns,nf,nis,nif

        !cccccccccccccccccccccccccccccccccccc
        !      Thomas algorithm begin
        !cccccccccccccccccccccccccccccccccccc
        do i = nis, nif
                a(i,ns) = 0.
                c(i,nf) = 0.
                f(i,ns) = c(i,ns) / b(i,ns)
                g(i,ns) = d(i,ns) / b(i,ns)
        end do
        !cccccccccccccccccccccccccccccccccccc
        !       forward elimination
        !cccccccccccccccccccccccccccccccccccc
        do j = ns+1, nf
        do i = nis , nif
                f(i,j) = c(i,j) /(b(i,j) - a(i,j)*f(i,j-1))
                g(i,j) =(d(i,j)-a(i,j)*g(i,j-1))/(b(i,j)-a(i,j)*f(i,j-1))
        end do
        end do
        !ccccccccccccccccccccccccccccccccccccccccccccccccc
        !   back substitution (solution put into d(i,j))
        !ccccccccccccccccccccccccccccccccccccccccccccccccc
        do i = nis , nif
                d(i,nf) = g(i,nf)
        end do
        do j = nf-1, ns, -1
        do i = nis , nif
                d(i,j) = g(i,j) - f(i,j)*d(i,j+1)
        end do
        end do

	return
end
