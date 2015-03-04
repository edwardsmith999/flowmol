!=======================================================================
! Controls the poisson solver
!
! poisson_init()
! waveno_init()
! poisson_fftz_mtgxy()
!

module poisson
	use data_export
	include "mpif.h"

	real :: start_time, end_time

	real, dimension(:),    allocatable :: z_mt
	real, dimension(:,:),  allocatable :: p_temp
	real, dimension(:),    allocatable :: a_mt,b_mt,c_mt,d_mt,f_mt,g_mt
	dimension                          :: myidp(nzm),cpu_wave(nzm)
	real*8                             :: array_fft(0:ngzm+2-1),coeff_fft(ngzm+15)
        INTERFACE
                real*8 function realClock()
                end function realClock
        END INTERFACE

end module

!=======================================================================
subroutine poisson_init()
	use poisson

	real    dz2, twopi, arg
	real    omegak1, omegakn
	integer k, m

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Compute and store the modified wavenumber (vkz)
	!  for central differencing (used in FFT of Spanwise
	!  direction --> Page 40)
	!  (This replaces "subroutine waveno" in xiahua's code
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

	!cccccccccccccccccccccccccccccccccccccccc
	!  Initialize the multi-grid parameters
	!cccccccccccccccccccccccccccccccccccccccc

	call readInt("ngrid" , ngrid )
	call readInt("ncycle", ncycle)
	call readInt("npre"  , npre  )

	call readInt("nitmax" , nitmax )
	call readInt("nitmaxb", nitmaxb)

	call readFloat("resmin" , resmin )
	call readFloat("resminb", resminb)

	call readFloat("omegak1", omegak1)
	call readFloat("omegakn", omegakn)
	omegak(1)       = omegak1
	omegak(2:ngz-1) = omegakn

	call preprocess_mt()

	return
end


subroutine poisson_fftz_mtgxy()
	!cccccccccccccccccccccccccc
	! solve poisson equation
	!cccccccccccccccccccccccccc
	use poisson

	array_fft=0.
	coeff_fft=0.

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!    evaluate the source term using velocity gradient
	!    Source term is nothing but the Divergence of Ust  (D u*)
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	before = realClock()
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

	!T tmp = sum(qr(ka:kb,ia:ib,ja:jb));
	!T print*,'Total sum of RHS  = ', tmp
	!T call globalSum(tmp,1);
	!T print*,'Global sum = ', tmp

	after = realClock()
	cpu_source= (after - before) - t_over

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!  Fourier transform the source term q(i,j,k) in z-direction 
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	before = realClock()
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

	after = realClock()
	cpu_fft1= (after - before) - t_over
	if(ntime-ntime_.eq.0) cpu_wavet=0.0
	before = realClock()

	!------------------------------------
	!   Time for Transposing qr in qT
	!------------------------------------
	allocate(qT(0:ngx, 0:ngy, nlzm))
	start_time = MPI_Wtime()
	call transpose_qr_qT(qr,qT)
	end_time = MPI_Wtime()
	if (irank.eq.1) write(6,*) 'Transpose time =', 2*(end_time-start_time)


	if (irank.eq.1) then
		tmp = sum(qT(1:ngxm,1:ngym,1));
		print*,'Total sum of RHS for (k=0) = ', iblock,jblock,tmp
	end if

!------------------------------------------------------------------------------------------
!			START SOLUTION LOOP OVER WAVE NUMBERS
!------------------------------------------------------------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled system in wave space 
!       real part first (k < nzm/2+1)
!cccccccccccccccccccccccccccccccccccccccccccccc
ka = kTmin_2(irank)
kb = kTmax_2(irank)
do 20 k = ka,kb
	kk = k-ka+1
	allocate(z_mt(memlenap_z))
	allocate(p_temp(0:ngx,0:ngy))
	allocate(a_mt(ngxm),b_mt(ngxm),c_mt(ngxm),d_mt(ngxm),f_mt(ngxm),g_mt(ngxm))
       
	beforeK=realClock()
!	if(k.eq.1) call system("ps -l -u tzaki | grep 'LPTrun'")

	!cccccccccccccccccccccccccccccccccccccccccccccc
	! allocate memory for solution, source term,
	!  and residual at all levels
	!cccccccccccccccccccccccccccccccccccccccccccccc
	if(k.le.ngzm/2+1) then
		goto 8000
	else
		ki=k-ngzm/2
		goto 8001
	end if
 8000   nc=0
	iprob_mt=iprob
	resminb_mt=resminb
	nitmaxb_mt=nitmaxb
	npre_mt=npre
	!cccccccccccccccccccccccccccccccccccccccccccccc
	! set up right-hand side at the finest level
	!cccccccccccccccccccccccccccccccccccccccccccccc
	do j=0,ngy 
	do i=0,ngx
		ij=i+(ngxm+2)*j 
		z_mt(ires_mt(ngrid)+ij)=0.0
		z_mt(irhs_mt(ngrid)+ij)=qT   (i,j,kk)
		z_mt(ip_mt  (ngrid)+ij)=phatr(i,j,kk)
	end do
	end do 
	nnx=ngxm
	nny=ngym

	rhs_mtmax=-1.e7
	do j=1,ngy-1
	do i=1,ngx-1
		rhs_mtmax=max(abs(qT(i,j,kk)),rhs_mtmax)
	end do
	end do
	if(k.eq.1.and.mod(ntime-ntime_,ipout).eq.0) print*, 'max rhs (1:ngx-1, 1:ngy-1) = ' , rhs_mtmax

	rhs_mtmax=-1.e7
	do j=0,ngy
	do i=0,ngx
		rhs_mtmax=max(abs(qT(i,j,kk)),rhs_mtmax)
	end do
	end do
	if(k.eq.1.and.mod(ntime-ntime_,ipout).eq.0) print*, 'max rhs (0:ngx, 0:ngy) = ' , rhs_mtmax

	!cccccccccccccccccccccccccc
	! start V-cycle loop here
	!cccccccccccccccccccccccccc
	do 15 mcycle=1,ncycle
		nfx=nnx
		nfy=nny 
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! going downhill from the current level to the coarset level
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          	do mm=ngrid,2,-1
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! pre-smoothing, compute solution for npre iterations
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			call relax(	z_mt  (ip_mt   (mm)),z_mt  (irhs_mt (mm)),z_mt  (ires_mt (mm)),nfx,nfy,k, &
					zap_mt(iapw_mt (mm)),zap_mt(iape_mt (mm)),zap_mt(iaps_mt (mm)), &
					zap_mt(iapn_mt (mm)),zap_mt(iapws_mt(mm)),zap_mt(iapwn_mt(mm)), &
					zap_mt(iapes_mt(mm)),zap_mt(iapen_mt(mm)),zap_mt(iapp_mt (mm)), &
					zap_mt(ivp_mt  (mm)),zap_mt(iedge_mt(mm)),omegak(k),vkz(k), &
					p_temp(0,0),a_mt(1),b_mt(1),c_mt(1),d_mt(1),f_mt(1),g_mt(1), &
					iprob_mt,npre_mt)

			call resid_mt(	z_mt  (ires_mt (mm)),z_mt  (ip_mt   (mm)),z_mt  (irhs_mt (mm)),nfx,nfy, &
					zap_mt(iapw_mt (mm)),zap_mt(iape_mt (mm)),zap_mt(iaps_mt (mm)), &
					zap_mt(iapn_mt (mm)),zap_mt(iapws_mt(mm)),zap_mt(iapwn_mt(mm)), &
					zap_mt(iapes_mt(mm)),zap_mt(iapen_mt(mm)),zap_mt(iapp_mt (mm)), &
					zap_mt(ivp_mt(mm)),vkz(k))
			nfx=nfx/2
			nfy=nfy/2
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! restriction of residual is the rhs_mt in next coarse level
			! i.e. determine RHS of coarser grid using residual of fine grid.
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            		call restrict(z_mt(irhs_mt(mm-1)),z_mt(ires_mt(mm)),nfx,nfy) 
			
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! zero for initial guess in (correction) next relaxation
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			do jfill=0,(nfx+2)*(nfy+2)-1 
				z_mt(ip_mt(mm-1)+jfill)=0.0
			end do
		end do
		!ccccccccccccccccccccccccccccccccccccccccc
		! V-cycle reaches the bottom level
		!ccccccccccccccccccccccccccccccccccccccccc
		mm = 1 
		call relax_b(	z_mt  (ip_mt   (mm)),z_mt  (irhs_mt (mm)),z_mt  (ires_mt (mm)),nfx,nfy,k, &
				zap_mt(iapw_mt (mm)),zap_mt(iape_mt (mm)),zap_mt(iaps_mt (mm)), &
				zap_mt(iapn_mt (mm)),zap_mt(iapws_mt(mm)),zap_mt(iapwn_mt(mm)), &
				zap_mt(iapes_mt(mm)),zap_mt(iapen_mt(mm)),zap_mt(iapp_mt (mm)), &
				zap_mt(ivp_mt  (mm)),zap_mt(iedge_mt(mm)),omegak(k),vkz(k), &
				p_temp(0,0),a_mt(1),b_mt(1),c_mt(1),d_mt(1),f_mt(1),g_mt(1), &
				iprob_mt,resminb_mt,nitmaxb_mt)
		!ccccccccccccccccccccccccccccccccccccccccc
		! V-cycle going upward from the bottom
		!ccccccccccccccccccccccccccccccccccccccccc
		do mm=2,ngrid
			nfx=2*nfx
			nfy=2*nfy
			!ccccccccccccccccccccccccccccccccccccccccc
			! give coarse grid solution to fine grid
			!ccccccccccccccccccccccccccccccccccccccccc
			call addint(z_mt(ip_mt(mm)),z_mt(ip_mt(mm-1)), &
			z_mt(ires_mt(mm)),nfx,nfy,zap_mt(iedge_mt(mm)),p_temp(0,0),iprob_mt)
			!ccccccccccccccccccccccccccccccccccccccccc
			! do post-smoothing at the fine level
			!ccccccccccccccccccccccccccccccccccccccccc
			call relax(	z_mt  (ip_mt   (mm)),z_mt  (irhs_mt (mm)),z_mt  (ires_mt (mm)),nfx,nfy,k, &
					zap_mt(iapw_mt (mm)),zap_mt(iape_mt (mm)),zap_mt(iaps_mt (mm)), &
					zap_mt(iapn_mt (mm)),zap_mt(iapws_mt(mm)),zap_mt(iapwn_mt(mm)), &
					zap_mt(iapes_mt(mm)),zap_mt(iapen_mt(mm)),zap_mt(iapp_mt (mm)), &
					zap_mt(ivp_mt  (mm)),zap_mt(iedge_mt(mm)),omegak(k),vkz(k), &
					p_temp(0,0),a_mt(1),b_mt(1),c_mt(1),d_mt(1),f_mt(1),g_mt(1), &
					iprob_mt,npre_mt)
		end do
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! end do 15 below is number of V-cycles at j level, (E-S)-(E-S) is gamma=2
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		mm=ngrid
		call resid_mt(	z_mt  (ires_mt (mm)),z_mt  (ip_mt   (mm)),z_mt  (irhs_mt (mm)),nfx,nfy, &
				zap_mt(iapw_mt (mm)),zap_mt(iape_mt (mm)),zap_mt(iaps_mt (mm)), &
				zap_mt(iapn_mt (mm)),zap_mt(iapws_mt(mm)),zap_mt(iapwn_mt(mm)), &
				zap_mt(iapes_mt(mm)),zap_mt(iapen_mt(mm)),zap_mt(iapp_mt (mm)), &
				zap_mt(ivp_mt(mm)),vkz(k))

		res_mtmax=-1.e7
		do j=0,(ngym+2)*(ngxm+2)-1
				if (abs(z_mt(ires_mt(ngrid)+j)).gt.res_mtmax) then
					res_mtmax=abs(z_mt(ires_mt(ngrid)+j))
					iit = mod(j,ngxm+2)
					jit = j / (ngxm+2)
					rhs_mtmax=qT(iit,jit,kk)
				end if
		end do

	if(k.eq.1.and.mod(ntime-ntime_,ipout).eq.0) write(6,991) mcycle,res_mtmax, rhs_mtmax, iit, jit
 991    format(1x, 'cycle counter= No.   ', 5x, I5, 5x,       & 
              'maximum residual at original level=', 5x, 1pe12.5, 5x, 1pe12.5, I10, I10)
	if(res_mtmax.le.resmin) goto 17 
  15    continue

  17    continue

	!ccccccccccccccccccccccccccccccccc
	!  Recover "phatr" from "z_mt"
	!ccccccccccccccccccccccccccccccccc

	do j=0,ngy
	do i=0,ngx
		ij=i+(ngxm+2)*j
		phatr(i,j,kk)=z_mt(ip_mt(ngrid)+ij)
	end do
	end do
	goto 21

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! solve the decoupled system in wave space       imaginary part 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 8001   nc=0
	iprob_mt=iprob
	resminb_mt=resminb
	nitmaxb_mt=nitmaxb
	npre_mt=npre
	!cccccccccccccccccccccccccccccccccccccccccccccc
	! set up right-hand side at the finest level
	!cccccccccccccccccccccccccccccccccccccccccccccc
	do j=0,ngy 
	do i=0,ngx
		ij=i+(ngxm+2)*j 
		z_mt(ires_mt(ngrid)+ij)=0.0
		z_mt(irhs_mt(ngrid)+ij)=qT(i,j,kk)
		z_mt(ip_mt(ngrid)+ij)=phatr(i,j,kk)
	end do
	end do 
	nnx=ngxm
	nny=ngym
	!ccccccccccccccccccccccccccccccc
	! start V-cycle loop here
	!ccccccccccccccccccccccccccccccc
	do 215 mcycle=1,ncycle
		nfx=nnx
		nfy=nny 
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! going downhill from the current level to the coarset level
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		do mm=ngrid,2,-1
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! pre-smoothing, compute solution for npre iterations
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			call relax(	z_mt  (ip_mt   (mm)),z_mt  (irhs_mt (mm)),z_mt  (ires_mt (mm)),nfx,nfy,ki, &
					zap_mt(iapw_mt (mm)),zap_mt(iape_mt (mm)),zap_mt(iaps_mt (mm)), &
					zap_mt(iapn_mt (mm)),zap_mt(iapws_mt(mm)),zap_mt(iapwn_mt(mm)), &
					zap_mt(iapes_mt(mm)),zap_mt(iapen_mt(mm)),zap_mt(iapp_mt (mm)), &
					zap_mt(ivp_mt  (mm)),zap_mt(iedge_mt(mm)),omegak(ki),vkz(ki), &
					p_temp(0,0),a_mt(1),b_mt(1),c_mt(1),d_mt(1),f_mt(1),g_mt(1), &
					iprob_mt,npre_mt)
			
			call resid_mt(	z_mt  (ires_mt (mm)),z_mt  (ip_mt   (mm)),z_mt  (irhs_mt (mm)),nfx,nfy, &
					zap_mt(iapw_mt (mm)),zap_mt(iape_mt (mm)),zap_mt(iaps_mt (mm)), &
					zap_mt(iapn_mt (mm)),zap_mt(iapws_mt(mm)),zap_mt(iapwn_mt(mm)), &
					zap_mt(iapes_mt(mm)),zap_mt(iapen_mt(mm)),zap_mt(iapp_mt (mm)), &
					zap_mt(ivp_mt  (mm)),vkz(ki))
			nfx=nfx/2
			nfy=nfy/2
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! restriction of residual is the rhs_mt in next coarse level
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			call restrict(z_mt(irhs_mt(mm-1)),z_mt(ires_mt(mm)),nfx,nfy) 
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! zero for initial guess in (correction) next relaxation
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			do jfill=0,(nfx+2)*(nfy+2)-1 
				z_mt(ip_mt(mm-1)+jfill)=0.0
			end do
		end do
		!cccccccccccccccccccccccccccccccccccc
		! V-cycle reaches the bottom level
		!cccccccccccccccccccccccccccccccccccc
		mm = 1 
		call relax_b(	z_mt  (ip_mt   (mm)),z_mt  (irhs_mt (mm)),z_mt  (ires_mt (mm)),nfx,nfy,ki, &
				zap_mt(iapw_mt (mm)),zap_mt(iape_mt (mm)),zap_mt(iaps_mt (mm)), &
				zap_mt(iapn_mt (mm)),zap_mt(iapws_mt(mm)),zap_mt(iapwn_mt(mm)), &
				zap_mt(iapes_mt(mm)),zap_mt(iapen_mt(mm)),zap_mt(iapp_mt (mm)), &
				zap_mt(ivp_mt  (mm)),zap_mt(iedge_mt(mm)),omegak(ki),vkz(ki), &
				p_temp(0,0),a_mt(1),b_mt(1),c_mt(1),d_mt(1),f_mt(1),g_mt(1), &
				iprob_mt,resminb_mt,nitmaxb_mt)
		!ccccccccccccccccccccccccccccccccccccccccc
		! V-cycle going upward from the bottom
		!ccccccccccccccccccccccccccccccccccccccccc
		do mm=2,ngrid
			nfx=2*nfx
			nfy=2*nfy
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! give coarse grid solution to fine grid
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			call addint(	z_mt(ip_mt(mm)),z_mt(ip_mt(mm-1)),               &  
					z_mt(ires_mt(mm)),nfx,nfy,zap_mt(iedge_mt(mm)),p_temp(0,0),iprob_mt)
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! do post-smoothing at the fine level
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			call relax(	z_mt  (ip_mt   (mm)),z_mt  (irhs_mt (mm)),z_mt  (ires_mt (mm)),nfx,nfy,ki, &
					zap_mt(iapw_mt (mm)),zap_mt(iape_mt (mm)),zap_mt(iaps_mt (mm)), &
					zap_mt(iapn_mt (mm)),zap_mt(iapws_mt(mm)),zap_mt(iapwn_mt(mm)), &
					zap_mt(iapes_mt(mm)),zap_mt(iapen_mt(mm)),zap_mt(iapp_mt (mm)), &
					zap_mt(ivp_mt  (mm)),zap_mt(iedge_mt(mm)),omegak(ki),vkz(ki),p_temp(0,0), &
					a_mt(1),b_mt(1),c_mt(1),d_mt(1),f_mt(1),g_mt(1), &
					iprob_mt,npre_mt)
		end do
          	mm=ngrid
		call resid_mt(	z_mt  (ires_mt (mm)),z_mt  (ip_mt   (mm)),z_mt  (irhs_mt (mm)),nfx,nfy, &
				zap_mt(iapw_mt (mm)),zap_mt(iape_mt (mm)),zap_mt(iaps_mt (mm)), &
				zap_mt(iapn_mt (mm)),zap_mt(iapws_mt(mm)),zap_mt(iapwn_mt(mm)), &
				zap_mt(iapes_mt(mm)),zap_mt(iapen_mt(mm)),zap_mt(iapp_mt (mm)), &
				zap_mt(ivp_mt  (mm)),vkz(ki))

		res_mtmax=-1.e7
		do j=0,(ngym+2)*(ngxm+2)-1
				res_mtmax=max(abs(z_mt(ires_mt(ngrid)+j)),res_mtmax)
		end do

	if(ki.eq.2.and.mod(ntime-ntime_,ipout).eq.0) write(6,992) mcycle,res_mtmax
  992  	format(1x, 'cycle counter= No_i   ', 5x, I5, 5x,       &
                 'maximum residual at original level=', 5x, 1pe12.5)
	if(res_mtmax.le.resmin) goto 217 
  215   continue
  217   continue

	!ccccccccccccccccccccccccccccccc
	!  Recover "phatr" from "z_mt"
	!ccccccccccccccccccccccccccccccc

	do j=0,ngy
	do i=0,ngx
		ij=i+(ngxm+2)*j
		phatr(i,j,kk)=z_mt(ip_mt(ngrid)+ij)
	end do
	end do

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! (21 continue) is finishing a wave number (either real or imaginary) 
	! and need to deallocate memory and calculate computational time.
	! (It does NOT an end for a (do-loop) start).
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 21     continue

	after=realClock()
	cpu_wave(k)=after-beforeK-t_over
	cpu_wavet(k)=cpu_wavet(k)+cpu_wave(k)
	deallocate(a_mt,b_mt,c_mt,d_mt,f_mt,g_mt)
	deallocate(p_temp,z_mt)
	
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (20 continue) is finishing all wave numbers, real and imaginary
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 20   continue


!------------------------------------------------------------------------------------------
!				FINISHED ALL WAVENUMBERS
!------------------------------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! TJ: enforce x-periodicity on phatr here
!	- physical pressures are defined on the x-interval [1,ngx-1]
!	-   halo   pressures are defined at (i.eq.0) and (i.eq.ngx)  
!	- phatr(0:ngx, 0:ngy, nlzm  )
!	- the global dimensions of phatr means we can enforce
!	  x-periodicity without calling routines from the messenger
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	phatr( 0,  :, :) = phatr(ngx-1,:, :)	!left  halo
	phatr(ngx, :, :) = phatr(  1,  :, :)	!right halo

	after = realClock()
	cpu_helmholtz= (after - before) - t_over
	before = realClock()

	!------------------------------------
	!  Time for Transposing phatr in p 
	!------------------------------------
	deallocate(qT)

	call transpose_phat_p(phatr,p)

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

	after = realClock()
	cpu_fft2= (after - before) - t_over
	!ccccccccccccccccccccccccccccccccccccccccc
	! make pressure periodic in z directions
	!ccccccccccccccccccccccccccccccccccccccccc
	p(ngz,:,:)=p(1   ,:,:)
	p(0  ,:,:)=p(ngzm,:,:)

	return
end


