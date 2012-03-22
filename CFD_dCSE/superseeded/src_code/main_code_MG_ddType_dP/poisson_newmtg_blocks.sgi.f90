!=======================================================================
! Controls the poisson solver
!
! poisson_init()
! waveno_init()
! poisson_fftz_mtgxy()
! poisson_free()
!

module poisson
	use data_export
	include "mpif.h"

	real, dimension(:),    allocatable :: a_mt,b_mt,c_mt,d_mt,f_mt,g_mt
	dimension                          :: myidp(nzm),cpu_wave(nzm)

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!SY 1D array for fft 
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8                             :: array_fft(0:ngzm+2-1),coeff_fft(ngzm+15)

        INTERFACE
                real*8 function realClock()
                end function realClock
        END INTERFACE

end module

!=======================================================================
subroutine poisson_init()
	use poisson

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !  parameters for the modified wavenumber
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real    dz2, twopi, arg
        integer k, m

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !  parameters for relaxation factor
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real    omegak1, omegakn

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

        !cccccccccccccccccccccccccccccccccccccccc
        !SY  relaxation factor
        !cccccccccccccccccccccccccccccccccccccccc
	omegak(1)       = omegak1
	omegak(2:ngz-1) = omegakn

        !ccccccccccccccccccccccccccccccccccccccccccccccccccc
        !SY  allocate the coefficients of Poisson equation
        !ccccccccccccccccccccccccccccccccccccccccccccccccccc
        call allocate_mt()

        !ccccccccccccccccccccccccccccccccccccccccccccccccccc
        !SY  calculate the coefficients of Poisson equation
        !ccccccccccccccccccccccccccccccccccccccccccccccccccc
        call preprocess_mt()

	return
end
!----------------------------------------

subroutine allocate_mt()
        use poisson

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! allocate the coefficients of poisson equation
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do m=ngrid,1,-1
                mgl(m)%nx = nixp*2**float(m-ngrid)
                mgl(m)%ny = niyp*2**float(m-ngrid)

                nxl = mgl(m)%nx
                nyl = mgl(m)%ny
                allocate(mgl(m)%zp   (0:nxl+1,0:nyl+1))
                allocate(mgl(m)%zrhs (0:nxl+1,0:nyl+1))
                allocate(mgl(m)%zres (0:nxl+1,0:nyl+1))
                allocate(mgl(m)%zapw (nxl,nyl), mgl(m)%zape (nxl,nyl))
                allocate(mgl(m)%zaps (nxl,nyl), mgl(m)%zapn (nxl,nyl))
                allocate(mgl(m)%zapws(nxl,nyl), mgl(m)%zapwn(nxl,nyl))
                allocate(mgl(m)%zapes(nxl,nyl), mgl(m)%zapen(nxl,nyl))
                allocate(mgl(m)%zapp (nxl,nyl))
                allocate(mgl(m)%zvp  (0:nxl+1,0:nyl+1))
                allocate(mgl(m)%zsux (0:nxl+2,0:nyl+1))
                allocate(mgl(m)%zsuy (0:nxl+2,0:nyl+1))
                allocate(mgl(m)%zsvx (0:nxl+1,0:nyl+2))
                allocate(mgl(m)%zsvy (0:nxl+1,0:nyl+2))
        end do

	return
end
!----------------------------------------

subroutine poisson_fftz_mtgxy()

	!cccccccccccccccccccccccccc
	! solve poisson equation
	!cccccccccccccccccccccccccc
	use poisson
	integer :: nlxymax

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

	after = realClock()
	cpu_source= (after - before) - t_over

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!   Fourier transform the source term q(i,j,k) in z-direction
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

	!------------------------------------------------------------------------------------------
	!			START SOLUTION LOOP OVER WAVE NUMBERS
	!------------------------------------------------------------------------------------------

	!cccccccccccccccccccccccccccccccccccccccccccccc
	! solve the decoupled system in wave space
	!       real part first (k < nzm/2+1)
	!cccccccccccccccccccccccccccccccccccccccccccccc
	ka = 1
	kb = ngzm/2+1
	do 20 k = ka,kb

	nlxymax = max(nlxm,nlym)
	allocate(a_mt(nlxymax),b_mt(nlxymax),c_mt(nlxymax),d_mt(nlxymax),f_mt(nlxymax),g_mt(nlxymax))

	beforeK=realClock()
	!if(k.eq.1) call system("ps -l -u tzaki | grep 'LPTrun'")

	iprob_mt=iprob
	resminb_mt=resminb
	nitmaxb_mt=nitmaxb
	npre_mt=npre
	resmin_mt=resmin
	nitmax_mt=nitmax

	!ccccccccccccccccccccccccccccccccccccccccccccccccc
	! set up right-hand side at the finest level
	!ccccccccccccccccccccccccccccccccccccccccccccccccc
	do j=0,niyp+1
	do i=0,nixp+1
		ja = j-0+(j1_T-1)
		ia = i-0+(i1_T-1)
		mgl(ngrid)%zres(i,j)=0.0
		mgl(ngrid)%zrhs(i,j)=qr(k,ia,ja)
		mgl(ngrid)%zp  (i,j)=phatr(k,i,j)
	end do
	end do

	!cccccccccccccccccccccccccc
	! start V-cycle loop here
	!cccccccccccccccccccccccccc
	do 15 mcycle=1,ncycle

		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! going downhill from the current level to the coarset level
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          	do mm=ngrid,2,-1

          	nfx=mgl(mm)%nx
          	nfy=mgl(mm)%ny
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! pre-smoothing, compute solution for npre iterations
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			call relax(	mgl(mm)%zp   , mgl(mm)%zrhs , mgl(mm)%zres , nfx, nfy, k, &
					mgl(mm)%zapw , mgl(mm)%zape , mgl(mm)%zaps ,              &
			                mgl(mm)%zapn , mgl(mm)%zapws, mgl(mm)%zapwn,              &
			                mgl(mm)%zapes, mgl(mm)%zapen, mgl(mm)%zapp ,              &
		 	                mgl(mm)%zvp  , mgl(mm)%zedge, omegak(k)    , vkz(k),      &
					a_mt(1), b_mt(1), c_mt(1), d_mt(1), f_mt(1), g_mt(1),     &
					iprob_mt, npre_mt, resmin_mt, nitmax_mt, mm )

			call resid_mt(	mgl(mm)%zres , mgl(mm)%zp   , mgl(mm)%zrhs , nfx, nfy, &
					mgl(mm)%zapw , mgl(mm)%zape , mgl(mm)%zaps ,           &
					mgl(mm)%zapn , mgl(mm)%zapws, mgl(mm)%zapwn,           &
					mgl(mm)%zapes, mgl(mm)%zapen, mgl(mm)%zapp ,           &
					mgl(mm)%zvp  , vkz(k) )

			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! restriction of residual is the rhs_mt in next coarse level
			! i.e. determine RHS of coarser grid using residual of fine grid.
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          	        nfx=mgl(mm-1)%nx
          	        nfy=mgl(mm-1)%ny

			call restrict(mgl(mm-1)%zrhs, mgl(mm)%zres, nfx, nfy)

			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! zero for initial guess in (correction) next relaxation
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			mgl(mm-1)%zp(0:nfx+1, 0:nfy+1)=0.0

		end do
		!ccccccccccccccccccccccccccccccccccccccccc
		! V-cycle reaches the bottom level
		!ccccccccccccccccccccccccccccccccccccccccc
		mm = 1
			call relax_b(	mgl(mm)%zp   , mgl(mm)%zrhs , mgl(mm)%zres , nfx, nfy, k, &
					mgl(mm)%zapw , mgl(mm)%zape , mgl(mm)%zaps ,              &
			                mgl(mm)%zapn , mgl(mm)%zapws, mgl(mm)%zapwn,              &
			                mgl(mm)%zapes, mgl(mm)%zapen, mgl(mm)%zapp ,              &
			                mgl(mm)%zvp  , mgl(mm)%zedge, omegak(k)    , vkz(k),      &
					a_mt(1), b_mt(1), c_mt(1), d_mt(1), f_mt(1), g_mt(1),     &
					iprob_mt     , resminb_mt, nitmaxb_mt )

		!ccccccccccccccccccccccccccccccccccccccccc
		! V-cycle going upward from the bottom
		!ccccccccccccccccccccccccccccccccccccccccc
		do mm=2,ngrid

         	nfx=mgl(mm)%nx
         	nfy=mgl(mm)%ny
			!ccccccccccccccccccccccccccccccccccccccccc
			! give coarse grid solution to fine grid
			!ccccccccccccccccccccccccccccccccccccccccc
			call addint(	mgl(mm)%zp  , mgl(mm-1)%zp, &
					mgl(mm)%zres, nfx, nfy, mgl(mm)%zedge, iprob_mt)

			!ccccccccccccccccccccccccccccccccccccccccc
			! do post-smoothing at the fine level
			!ccccccccccccccccccccccccccccccccccccccccc
			call relax(	mgl(mm)%zp   , mgl(mm)%zrhs , mgl(mm)%zres , nfx, nfy, k, &
					mgl(mm)%zapw , mgl(mm)%zape , mgl(mm)%zaps ,              &
			                mgl(mm)%zapn , mgl(mm)%zapws, mgl(mm)%zapwn,              &
			                mgl(mm)%zapes, mgl(mm)%zapen, mgl(mm)%zapp ,              &
			                mgl(mm)%zvp  , mgl(mm)%zedge, omegak(k)    , vkz(k),      &
					a_mt(1), b_mt(1), c_mt(1), d_mt(1), f_mt(1), g_mt(1),     &
					iprob_mt, npre_mt, resmin_mt, nitmax_mt, mm )

		end do
		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! end do 15 below is number of V-cycles at j level, (E-S)-(E-S) is gamma=2
		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		mm=ngrid

		call resid_mt(	mgl(mm)%zres , mgl(mm)%zp   , mgl(mm)%zrhs , nfx, nfy, &
				mgl(mm)%zapw , mgl(mm)%zape , mgl(mm)%zaps ,           &
				mgl(mm)%zapn , mgl(mm)%zapws, mgl(mm)%zapwn,           &
				mgl(mm)%zapes, mgl(mm)%zapen, mgl(mm)%zapp ,           &
				mgl(mm)%zvp  , vkz(k) )

		!cccccccccccccccccccccccccccccccccccccccccccccc
		! find the maximum of residual
		!cccccccccccccccccccccccccccccccccccccccccccccc
		res_mtmax=-1.e7

		call find_max_res(mgl(mm)%zres, res_mtmax)

	if(irank.eq.1 .and. k.eq.1 .and. mod(ntime-ntime_,ipout).eq.0) write(6,991) mcycle,res_mtmax

 991    format(1x, 'cycle counter= No.   ', 5x, I5, 5x,       &
              'maximum residual at original level=', 5x, 1pe12.5)
	if(res_mtmax.le.resmin) exit
  15    continue

	!ccccccccccccccccccccccccccccccccc
	!  Recover "phatr" from "mgl%zp"
	!ccccccccccccccccccccccccccccccccc
	do j=0,niyp+1
	do i=0,nixp+1
	   phatr(k,i,j)=mgl(ngrid)%zp(i,j)
	end do
	end do

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! (21 continue) is finishing a wave number (real)
        ! and need to deallocate memory and calculate computational time.
        ! (It does NOT an end for a (do-loop) start).
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   21   continue

        after=realClock()
        cpu_wave(k)=after-beforeK-t_over
        cpu_wavet(k)=cpu_wavet(k)+cpu_wave(k)
        deallocate(a_mt,b_mt,c_mt,d_mt,f_mt,g_mt)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (20 continue) is finishing all wave numbers, real
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   20   continue

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! solve the decoupled system in wave space
	!       imaginary part
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ka = ngzm/2+2
        kb = ngzm
        do 220 k = ka,kb
        ki=k-ngzm/2

        nlxymax = max(nlxm,nlym)
        allocate(a_mt(nlxymax),b_mt(nlxymax),c_mt(nlxymax),d_mt(nlxymax),f_mt(nlxymax),g_mt(nlxymax))

        beforeK=realClock()

	iprob_mt=iprob
	resminb_mt=resminb
	nitmaxb_mt=nitmaxb
	npre_mt=npre
	resmin_mt=resmin
	nitmax_mt=nitmax

	!cccccccccccccccccccccccccccccccccccccccccccccc
	! set up right-hand side at the finest level
	!cccccccccccccccccccccccccccccccccccccccccccccc
	do j=0,niyp+1
	do i=0,nixp+1
		ja = j-0+(j1_T-1)
		ia = i-0+(i1_T-1)
		mgl(ngrid)%zres(i,j)=0.0
		mgl(ngrid)%zrhs(i,j)=qr(k,ia,ja)
		mgl(ngrid)%zp  (i,j)=phatr(k,i,j)
	end do
	end do

	!ccccccccccccccccccccccccccccccc
	! start V-cycle loop here
	!ccccccccccccccccccccccccccccccc
	do 215 mcycle=1,ncycle

		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! going downhill from the current level to the coarset level
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		do mm=ngrid,2,-1
		    	nfx=mgl(mm)%nx
		        nfy=mgl(mm)%ny

			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! pre-smoothing, compute solution for npre iterations
			!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			call relax(	mgl(mm)%zp   , mgl(mm)%zrhs , mgl(mm)%zres , nfx, nfy, ki, &
					mgl(mm)%zapw , mgl(mm)%zape , mgl(mm)%zaps ,               &
			                mgl(mm)%zapn , mgl(mm)%zapws, mgl(mm)%zapwn,               &
			                mgl(mm)%zapes, mgl(mm)%zapen, mgl(mm)%zapp ,               &
			                mgl(mm)%zvp  , mgl(mm)%zedge, omegak(ki)   , vkz(ki),      &
					a_mt(1), b_mt(1), c_mt(1), d_mt(1), f_mt(1), g_mt(1),      &
					iprob_mt, npre_mt, resmint_mt, nitmax_mt, mm )

			call resid_mt(	mgl(mm)%zres , mgl(mm)%zp   , mgl(mm)%zrhs , nfx, nfy, &
					mgl(mm)%zapw , mgl(mm)%zape , mgl(mm)%zaps ,           &
					mgl(mm)%zapn , mgl(mm)%zapws, mgl(mm)%zapwn,           &
					mgl(mm)%zapes, mgl(mm)%zapen, mgl(mm)%zapp ,           &
					mgl(mm)%zvp  , vkz(ki) )

			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! restriction of residual is the rhs_mt in next coarse level
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			nfx=mgl(mm-1)%nx
			nfy=mgl(mm-1)%ny

			call restrict(mgl(mm-1)%zrhs, mgl(mm)%zres, nfx, nfy)

			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! zero for initial guess in (correction) next relaxation
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			mgl(mm-1)%zp(0:nfx+1, 0:nfy+1)=0.0

		end do

		!cccccccccccccccccccccccccccccccccccc
		! V-cycle reaches the bottom level
		!cccccccccccccccccccccccccccccccccccc
			mm = 1
			call relax_b(	mgl(mm)%zp   , mgl(mm)%zrhs , mgl(mm)%zres , nfx, nfy, ki, &
					mgl(mm)%zapw , mgl(mm)%zape , mgl(mm)%zaps ,               &
			                mgl(mm)%zapn , mgl(mm)%zapws, mgl(mm)%zapwn,               &
			                mgl(mm)%zapes, mgl(mm)%zapen, mgl(mm)%zapp ,               &
			                mgl(mm)%zvp  , mgl(mm)%zedge, omegak(ki)   , vkz(ki),      &
					a_mt(1), b_mt(1), c_mt(1), d_mt(1), f_mt(1), g_mt(1),      &
					iprob_mt     , resminb_mt, nitmaxb_mt )

		!ccccccccccccccccccccccccccccccccccccccccc
		! V-cycle going upward from the bottom
		!ccccccccccccccccccccccccccccccccccccccccc
		do mm=2,ngrid
			nfx=mgl(mm)%nx
			nfy=mgl(mm)%ny

			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! give coarse grid solution to fine grid
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			call addint(	mgl(mm)%zp  , mgl(mm-1)%zp, &
					mgl(mm)%zres, nfx, nfy, mgl(mm)%zedge, iprob_mt)

			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			! do post-smoothing at the fine level
			!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
			call relax(	mgl(mm)%zp   , mgl(mm)%zrhs , mgl(mm)%zres , nfx, nfy, ki, &
					mgl(mm)%zapw , mgl(mm)%zape , mgl(mm)%zaps ,               &
			                mgl(mm)%zapn , mgl(mm)%zapws, mgl(mm)%zapwn,               &
			                mgl(mm)%zapes, mgl(mm)%zapen, mgl(mm)%zapp ,               &
			                mgl(mm)%zvp  , mgl(mm)%zedge, omegak(ki)   , vkz(ki),      &
					a_mt(1), b_mt(1), c_mt(1), d_mt(1), f_mt(1), g_mt(1),      &
					iprob_mt, npre_mt, resmin_mt, nitmax_mt, mm )

		end do
		mm=ngrid
			call resid_mt(	mgl(mm)%zres , mgl(mm)%zp   , mgl(mm)%zrhs , nfx, nfy, &
					mgl(mm)%zapw , mgl(mm)%zape , mgl(mm)%zaps ,           &
					mgl(mm)%zapn , mgl(mm)%zapws, mgl(mm)%zapwn,           &
					mgl(mm)%zapes, mgl(mm)%zapen, mgl(mm)%zapp ,           &
					mgl(mm)%zvp  , vkz(ki) )

		!cccccccccccccccccccccccccccccccccccccccccccccc
		! find the maximum of residual
		!cccccccccccccccccccccccccccccccccccccccccccccc
		call find_max_res(mgl(mm)%zres, res_mtmax)

	if(irank.eq.1 .and. ki.eq.2 .and. mod(ntime-ntime_,ipout).eq.0) write(6,992) mcycle,res_mtmax
  992   format(1x, 'cycle counter= No_i   ', 5x, I5, 5x,       &
              'maximum residual at original level=', 5x, 1pe12.5)
	if(res_mtmax.le.resmin) exit
  215   continue

	!ccccccccccccccccccccccccccccccc
	!  Recover "phatr" from "mgl%zp"
	!ccccccccccccccccccccccccccccccc
	do j=0,niyp+1
	do i=0,nixp+1
	   phatr(k,i,j)=mgl(ngrid)%zp(i,j)
	end do
	end do

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! (221 continue) is finishing a wave number (imaginary)
        ! and need to deallocate memory and calculate computational time.
        ! (It does NOT an end for a (do-loop) start).
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  221   continue

        after=realClock()
        cpu_wave(k)=after-beforeK-t_over
        cpu_wavet(k)=cpu_wavet(k)+cpu_wave(k)
        deallocate(a_mt,b_mt,c_mt,d_mt,f_mt,g_mt)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (220 continue) is finishing all wave numbers, imaginary
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  220   continue

!------------------------------------------------------------------------------------------
!				FINISHED ALL WAVENUMBERS
!------------------------------------------------------------------------------------------

	after = realClock()
	cpu_helmholtz= (after - before) - t_over
	before = realClock()

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! Fourier transform back to physical space, z-dir first
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do j = 0, niyp+1
	do i = 0, nixp+1
		ia = i+i1_T-1
		ja = j+j1_T-1

		ka = ngzm/2 + 1
		kb = 2*(ngzm/2)

		array_fft = 0.
		array_fft(0:ngzm:2) = phatr(1:ka   , i,j)
		array_fft(3:kb-1:2) = phatr(ka+1:kb, i,j)
		!for2D  array_fft(1)  = 0.
		!for2D  array_fft(kb+1) = 0.

		call realFFTM(array_fft,ngzm,1,-1)
		!J---------------------------------------------------
		!J	call dzfft1dui(ngzm,coeff_fft)
		!J	call zdfft1du(+1,ngzm,array_fft,1,coeff_fft)
		!J---------------------------------------------------
		p(1:ngzm,ia,ja) = array_fft(0:ngzm-1)*dt
	end do
	end do

	after = realClock()
	cpu_fft2= (after - before) - t_over

	!P-----------------------------------------------------------------
	!P   Update boundary cells for pressure (Don't think it is needed)
	!P-----------------------------------------------------------------
	!P    !T-- (x,y) directions
	!P       call updateBorder_lim(p, ngz  ,nlx  ,nly  , id_x, 2, i1_T, i2_T, 1)
	!P       call updateBorder_lim(p, ngz  ,nlx  ,nly  , id_y, 3, j1_T, j2_T, 1)


	!ccccccccccccccccccccccccccccccccccccccccc
	! make pressure periodic in z directions
	!ccccccccccccccccccccccccccccccccccccccccc
	p(ngz,:,:)=p(1   ,:,:)
	p(0  ,:,:)=p(ngzm,:,:)

	return
end
!----------------------------------------

subroutine poisson_free()
        use poisson

        call deallocate_mt()

        return
end
!----------------------------------------

subroutine deallocate_mt()
        use poisson

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !SY Deallocate the coefficients of the Poisson equation
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do m=ngrid,1,-1
                deallocate(mgl(m)%zp   )
                deallocate(mgl(m)%zrhs )
                deallocate(mgl(m)%zres )
                deallocate(mgl(m)%zapw )
                deallocate(mgl(m)%zape )
                deallocate(mgl(m)%zapn )
                deallocate(mgl(m)%zaps )
                deallocate(mgl(m)%zapwn)
                deallocate(mgl(m)%zapws)
                deallocate(mgl(m)%zapen)
                deallocate(mgl(m)%zapes)
                deallocate(mgl(m)%zapp )
                deallocate(mgl(m)%zvp  )
                deallocate(mgl(m)%zsux )
                deallocate(mgl(m)%zsuy )
                deallocate(mgl(m)%zsvx )
                deallocate(mgl(m)%zsvy )
        end do

        return
end
