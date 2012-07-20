!=======================================================================
! Calculating turbulence statistics
!
! statistics_init()
! statistics_read()
! statistics_write()
! statistics_sample()
! spectra(dt, S, zspectra)*
! statistics_info(iunit)
!
! This code combines statistics.f90.Ez and statistics.f90.uhatwxyz (including scalar)
! This code computes:
!	uhat(omega,x,y,z) at specified x,y and all z
!	Ez(x,y,kz) = uhat uhat*
! Post-processing can then compute:
!	Using ~/research/codes/utils/uhatwxyz.f90):
!		uhat(omega,x,ky,kz) if computed at all y -> can decompose into OS modes
!		uhat(omega,x,y) after averaging over z
!		E(omega,x,y) averaged over z
!		E(omega,x,y,kz)
!	Using ~/research/codes/utils/Ezstat.f90
!		Ez(x,y,kz) -- 1D energy spectrum vs kz
!		Ruu(x,y,dz), Rvv, or Rww depending on value of iuvw
!=======================================================================
module statistics
	use data_export
	use mesh_export

	integer     istat            ! interval to calculate statistics
	integer     nstat            ! number of calls to statistics
	integer     nsample          ! number of statistics samples
	integer     nsample_         ! number of samples from previous run
	real        sampleTime       ! total integration time
	real        sampleTime_      ! previous integration time
	real	    VT		     ! molecular viscosity

	!TAZ real Ui(nx_,ny_,6)           ! velocity statistics
	!TAZ real Rij(nx_,ny_,6)          ! Reynolds stresses
	!TAZ real Tij(nx_,ny_,6)          ! viscous stresses
	real Ui (0:nlx+1, 0:nly+1, 6)          ! velocity statistics
	real Rij(0:nlx+1, 0:nly+1, 6)          ! Reynolds stresses
	!LES ONLY	real Tij(0:nlx+1, 0:nly+1, 6)          ! viscous stresses

	!------- Terms in TKE equation not calculated by Jacobs ------
	real eps(0:nlx+1, 0:nly+1, 1)		! Dissipation (\epsilon)
	real UiP(0:nlx+1, 0:nly+1, 3)		! Pressure Corrolations
	real UiUj2(0:nlx+1, 0:nly+1, 3)		! TKE transport (u*TKE, v*TKE, w*TKE)

	!------- Write direct access file, or to archive -----
	logical :: Direct_IO = .true.

	! Scalar statistics
	! mean (T), variance(t^2), turbulent scalar fluxes (-u_i t),
	! fluctuating-pressure/temperature gradient correlation (p dt/dx_i)
	!TAZ parameter (nscalar=9)
	!TAZ real SCij(nx_,ny_,nscalar)         ! scalar statistics
	!T_SC ONLY	parameter (nscalar=9)
	!T_SC ONLY	real SCij(0:nlx+1, 0:nly+1, nscalar)         ! scalar statistics

	! 1D energy spectrum E(kz) -- post-process with Ezstat.f90
	logical :: sampleEz=.true.	! Compute one-dimensional spectrum vs kz
	allocatable Ez(:,:,:)		! z energy spectra

	! Spectra variables
	parameter (inx = 4, iny = 4)	! # locations per processor at which to compute uhat,
					! should have nixb divisible by inx, etc.
	integer indx(inx), indy(iny)	! index for x,y locations where uhst is computed
	complex, allocatable :: uhat(:,:,:,:)	! uhat
	parameter (nfreq=52360)		! # frequencies at which uhat could be computed
	parameter (nfreq_=  78)		! # frequencies at which to actually compute uhat
	integer ifreq_(nfreq_)		! integer frequencies (see below for definition)
	integer :: iuvw=3		! Compute Ez, uhat with U(1),V(2) or W(3)
	integer nterms			! # terms computed in Fourier transform
	integer nft			! # complete Fourier transforms computed
	integer ilength			! direct access file record length
	real Tlow			! Period of lowest frequency to resolve
	real pi
	complex cfac, ii		! -i*2pi/T, i
	complex ctmp			! Temporary complex number
	character*20 fname		! File name for uhat

	! real TKE(nx_,ny_,6)          ! TKE equation statistics
	! real RSTij(nx_,ny_,6)        ! RST equation statistics

end module

!=======================================================================

subroutine statistics_init()
	use statistics

	real	:: Re, theta0, x0, LengthInv, LengthScale
	real	:: ReOS, omegaOS
	integer	:: iunit

	call readFloat("Re",Re)
	call readFloat("x0",x0)
	call readFloat("theta0", theta0)
	call readFloat("LengthInv",LengthInv)

        ! apparently this is not needed
	!LengthScale = 1.0/LengthInv

	VT = 1./Re
	!-------------------------------------
	!T  Used to be this way,
	!T  but theta0 is not really needed
	!-------------------------------------
	!Blasius	x0 = Re*(theta0/.664)**2		! For Blasius profile
	!Blasius	LengthScale = sqrt(2. / Re * x0)	! If change, remember yOS below too
	!Blasius		! print*,' x0 = Re (theta0/.664)^2  = ', x0
	!Blasius		! print*,'LengthScale = ',LengthScale
	!Blasius	
	!--------------------------------------------
	!T  New way of calculating (x0, LengthScale)
	!--------------------------------------------
	!Blasius	LengthScale = 1.0/LengthInv
	!Blasius	x0 = 0.5 * Re * LengthScale * LengthScale
	!Blasius		! print*,'LengthScale = ',LengthScale
	!Blasius		! print*,' x0 = 0.5 R LengthScale^2 = ', x0
	!Blasius	
	!--------------------------------------------

	!iunit = iopen()
	!open (iunit, file="os2D.dat", form="formatted", &
	!	status="old", iostat=ierr)
	!if (ierr .ne. 0) stop "An Orr-Sommerfeld data file is required"
	!read(iunit,*) ReOS, omegaOS
	!close (iclose(iunit))
	!omegaOS = omegaOS / LengthScale

	call readInt("istat", istat)

	nstat = 0
	nsample = 0
	nsample_ = nsample
	sampleTime = 0.
	sampleTime_ = sampleTime

	Ui	= 0.
	UiP	= 0.
	Rij	= 0.
	eps	= 0.
	UiUj2	= 0.
	!LES ONLY	Tij = 0.
	!T_SC ONLY	SCij = 0.

	! Compute Tlow to resolve down to Fmin = omega nu 1e6/Uinf^2 = 3.0
	! Set nfreq to resolve desired Fmax = Fmin*(nfreq/2)
	! For Fmax = Fmax_DNS -> Use nfreq = nw*2pi/T0*nu*1e6/Fmin (use Fmin from here)
	! Using nfreq small than this will result in aliasing (because the highest
	! frequencies won't be resolved and so will be aliased back in)
	! 	DNS: Fmin = 2pi/T0 * nu * 1e6 : Fmax = Fmin * nw/2
	!	Use dt = Tlow / nfreq
	! Often dt > dt_stability so use dt that is stable
	! and use istat = dt_statistics/dt_dns
	! BOTTOM LINE:
	!	nfreq sets dt and Fmax (need nfreq timesteps for complete transform)
	!	so set Fmin, nfreq = # timesteps for transform --> compute Fmax
	! 	or set Fmin -> find dt_stability -> set nfreq = Tlow/dt from this
	!	(then make nfreq an integer -> dt=Tlow/nfreq)
	! Memory use for uhat(niz,nfreq_,inx,iny)) = nproc*16*inx*iny*nfreq_*niz bytes

	Fmin = 3.	! User defined (also, nfreq above);  In DNS units {F = w 10^6 / Re}
	pi = 4.*atan(1.)
	call readFloat("Re",Re)
	Tlow = 2.*pi *1.e6 / (Re * Fmin)
	ii = cmplx(0.,1.)
	cfac = -ii*2.*pi/Tlow

	! Set frequencies at which to compute uhat
	! Could compute up to nfreq/2+1 frequencies
	! To compute all, set nfreq_ = nfreq/2 and (ifreq_(i)=i for i=1,nfreq/2)
	! No ifreq_ value should be above nfreq/2
	! omega = 2pi/Tlow * ifreq_
	! Fine on 0...ifreq1, coarser on ifreq1...ifreq2, even coarser on ifreq2..ifreq3
	ifreq_(1) = 0
	ifreq1 = 60	! integer frequencies to interpolate between
	ifreq2 = 210
	ifreq3 = 4050
	!  ifreq3 = nfreq/2

	if ((istat > 0) .and. (ifreq2 > ifreq3) .or. (ifreq3 > nfreq/2)) stop "Check frequencies"
	nfreq1 = 21	! # freq. between 0 and ifreq1 (inclusive)
	nfreq2 = 25	! # freq. between ifreq1 and ifreq2, increment=(ifreq2-ifreq1)/nfreq2
	nfreq3 = 32	! # freq. between ifreq2 and ifreq3, increment=(ifreq3-ifreq2)/nfreq3
	if ((istat>0) .and. nfreq_ .ne. (nfreq1 + nfreq2 + nfreq3)) stop "Wrong number of frequencies"
	do i=1,nfreq1
		ifreq_(i) = (i-1)*ifreq1/(nfreq1-1)
	enddo
	do i=1,nfreq2
		ifreq_(nfreq1+i) = ifreq1 + (i*(ifreq2-ifreq1))/nfreq2
	enddo
	do i=1,nfreq3
		ifreq_(nfreq1+nfreq2+i) = ifreq2 + (i*(ifreq3-ifreq2))/nfreq3
	enddo

	!  if (irank==iroot) print*,'i    ifreq_'
	do i=1,nfreq_
		!  if (irank==iroot) print*,i,ifreq_(i)
		if (irank==iroot .and. ifreq_(i) > nfreq/2) print*,'ifreq_(i) = ',i,ifreq_(i)
		if (irank==iroot .and. ifreq_(i) == 0 .and. i > 1) print*,'multiple zero frequencies'
	enddo

	! Set x,y locations (in local coordinates) at which to compute uhat
	do i=1,inx
		indx(i) = i1_T+nixb/inx*(i-1)
	enddo
	if (iblock .eq. 1) indx(1) = indx(1) ! - 1	! Compute at inlet

	! If npy=1, then easier to control sample locations
	! Interpolate j sample locations from jmin to jtop
	!TAZ  jtop = 132		! max. j location where energy statistics are taken
	jtop = j1_T+niyb	! max. j location where energy statistics are taken
	if (jtop > j1_T+niyb) stop "jtop too large"
	do j=1,iny
		indy(j) = j1_T+(jtop-j1_T)/(iny-1)*(j-1)
	enddo

	! Use record length of niz complex numbers
	inquire(iolength=ilength) ctmp
	ilength = ilength*niz

	if ((istat > 0) .and. (inx > nix) .or. (iny > niy)) stop "Check inx and iny"

	return
end

subroutine statistics_read()
	use statistics
	real, allocatable :: S(:,:,:)
	real rtmp
	integer npnts

	if (istat .ne. 0) then
		allocate (uhat(ngzm,nfreq_,inx,iny))
		uhat = 0.

		! 1D energy spectrum
	        if (sampleEz) then
                        allocate (Ez(0:ngzm/2, 0:nlx+1, 0:nly+1))
	                Ez = 0.
	        endif

		call archive_isDefined("nstat", iflag)
		if (iflag .ne. 0) then
			call readInt("nstat", nstat)
			call readInt("nsample", nsample)
			call readFloat("sampleTime", sampleTime)
			nsample_ = nsample
			sampleTime_ = sampleTime

			!----------------------------------------------------
			!	Read in statistics from Direct Access File,
			!		or from Arvchive
			!---------------------------------------------------- 
			if (Direct_IO) then
				call Read_Stats()
			else
				!---- Local limits projected to Background ---
				i1 = ibmap_1(0)	; 	i2 = ibmap_1(nlxb+1)
				j1 = jbmap_1(0)	;	j2 = jbmap_1(nlyb+1)

				!---- Read global, & copy out local arrays ----
				allocate (S(0:ngx+1, 0:ngy+1, 6))		! b/c Rij has 6-elements

				npnts = (ngx+2)*(ngy+2)*6
				call readArray("Ui", S, npnts)
				!NO_MPI  call scatterXY(S, Ui, 6)
				Ui(0:nlxb+1, 0:nlyb+1, :) = S(i1:i2 , j1:j2 , 1:6)

				npnts = (ngx+2)*(ngy+2)*3
				call readArray("UiP", S, npnts)
				!NO_MPI  call scatterXY(S, UiP, 3)
				UiP(0:nlxb+1, 0:nlyb+1, :) = S(i1:i2 , j1:j2 , 1:3)

				npnts = (ngx+2)*(ngy+2)*6
				call readArray("Rij", S, npnts)
				!NO_MPI  call scatterXY(S, Rij, 6)
				Rij(0:nlxb+1, 0:nlyb+1, :) = S(i1:i2 , j1:j2 , 1:6)

				npnts = (ngx+2)*(ngy+2)*1
				call readArray("eps", S, npnts)
				!NO_MPI  call scatterXY(S, eps, 1)
				eps(0:nlxb+1, 0:nlyb+1, 1) = S(i1:i2 , j1:j2 , 1)

				npnts = (ngx+2)*(ngy+2)*3
				call readArray("UiUj2", S, npnts)
				!NO_MPI  call scatterXY(S, UiUj2, 3)
				UiUj2(0:nlxb+1, 0:nlyb+1, 1:3) = S(i1:i2 , j1:j2 , 1:3)

				!LES ONLY	call readArray("Tij", S, npnts)
				!LES ONLY	!NO_MPI  call scatterXY(S, Tij, 6)
				!LES ONLY	Tij(0:nlxb+1, 0:nlyb+1, :) = S(i1:i2 , j1:j2 , :)

				deallocate (S)

				!T_SC ONLY	allocate (S(0:ngx+1, 0:ngy+1, nscalar))
				!T_SC ONLY	npnts = (ngx+2)*(ngy+2)*nscalar
				!T_SC ONLY	
				!T_SC ONLY	call readArray("SCij", S, npnts)
				!T_SC ONLY	!NO_MPI  call scatterXY(S, SCij, nscalar)
				!T_SC ONLY	SCij(0:nlxb+1, 0:nlyb+1, :) = S(i1:i2 , j1:j2 , :)
				!T_SC ONLY	
				!T_SC ONLY	deallocate (S)
			end if

			! Frequency spectra
			call archive_isDefined("nterms", iflag)
			if (iflag .ne. 0) then
				call readInt("nterms", nterms)
			else
				nterms = 0
			endif
			call archive_isDefined("nft", iflag)
			if (iflag .ne. 0) then
				call readInt("nft", nft)
			else
				nft = 0
			endif

			if (nterms > 0) then	! Continuing a Fourier transform in time
				if (irank==iroot) print*,'Continuing transform'
				! Direct access file of "uhat" is read using MPI2 IO
				fname = "uhat"
				CALL Read_uhat()
			else	! Beginning a new Fourier transform
				if (irank==iroot) print*,'Starting new transform'
			endif

			call archive_isDefined("iEz", iflag)
			if (sampleEz .and. (iflag .ne. 0)) then
				fname = "Ez"
				CALL Read_Ez()
			endif

		end if
	end if

	return
end

subroutine statistics_write()
	use statistics
        use messenger
        include "mpif.h"
	real, allocatable :: S(:,:,:)
	integer npnts

	call writeInt("istat", istat)

	if (istat .ne. 0) then
		call writeInt("nstat", nstat)
		call writeInt("nsample", nsample)
		call writeFloat("sampleTime", sampleTime)

		!----------------------------------------------------
		!	Write statistics out to Direct Access File,
		!		or to Arvchive
		!---------------------------------------------------- 
		if (Direct_IO) then
			call Write_Stats()
		else
			! Globalize statistics and save to archive

			allocate (S(0:ngx+1, 0:ngy+1, 6))	! b/c Rij has 6-elements

			npnts = (ngx+2) * (ngy+2) * 6
			call stat_allGatherXY(Ui, S, 6)
			call writeArray("Ui", S, npnts)

			npnts = (ngx+2) * (ngy+2) * 4
			call stat_allGatherXY(UiP, S, 4)
			call writeArray("Ui", S, npnts)

			npnts = (ngx+2) * (ngy+2) * 6
			call stat_allGatherXY(Rij, S, 6)
			! Zero invalid values
			S(imino,:,1) = 0.
			S(imaxo,:,1) = 0.
			S(:,jmino,2) = 0.
			S(:,jmaxo,2) = 0.
			S(imaxo,:,4) = 0.
			S(:,jmaxo,4) = 0.
			S(:,jmaxo,5) = 0.
			S(imaxo,:,6) = 0.
			call writeArray("Rij", S, npnts)

			npnts = (ngx+2) * (ngy+2) * 1
			call stat_allGatherXY(eps, S, 1)
			call writeArray("eps", S, npnts)

			npnts = (ngx+2) * (ngy+2) * 3
			call stat_allGatherXY(UiUj2, S, 3)
			call writeArray("UiUj2", S, npnts)

			!LES ONLY	call stat_allGatherXY(Tij, S, 6)
			!LES ONLY	! Zero invalid values
			!LES ONLY	S(imino,:,1) = 0.
			!LES ONLY	S(imaxo,:,1) = 0.
			!LES ONLY	S(:,jmino,2) = 0.
			!LES ONLY	S(:,jmaxo,2) = 0.
			!LES ONLY	S(:,jmino,3) = 0.
			!LES ONLY	S(:,jmaxo,3) = 0.
			!LES ONLY	S(imaxo,:,4) = 0.
			!LES ONLY	S(:,jmaxo,4) = 0.
			!LES ONLY	S(:,jmaxo,5) = 0.
			!LES ONLY	S(imaxo,:,6) = 0.
			!LES ONLY	call writeArray("Tij", S, npnts)

			deallocate (S)

			!T_SC ONLY	allocate (S(0:ngx+1, 0:ngy+1 ,nscalar))
			!T_SC ONLY	npnts = (ngx+2) * (ngy+2) * nscalar
			!T_SC ONLY	
			!T_SC ONLY	call stat_allGatherXY(SCij, S, nscalar)
			!T_SC ONLY	! Zero invalid values
			!T_SC ONLY	S(imaxo,:,3) = 0.
			!T_SC ONLY	S(:,jmaxo,4) = 0.
			!T_SC ONLY	S(imaxo,:,6) = 0.
			!T_SC ONLY	S(:,jmaxo,7) = 0.
			!T_SC ONLY	S(imaxo,:,8) = 0.
			!T_SC ONLY	S(:,jmaxo,9) = 0.
			!T_SC ONLY	call writeArray("SCij", S, npnts)
			!T_SC ONLY	
			!T_SC ONLY	deallocate (S)
		end if

		! Frequency spectra
		call writeInt("nterms", nterms)
		call writeInt("nft", nft)
		call writeInt("nfreq", nfreq)
		call writeInt("nfreq_", nfreq_)
		call writeFloat("Tlow", Tlow)
		call writeIntArray("ifreq_",ifreq_,nfreq_)
		call writeInt("inx", inx)
		call writeInt("iny", iny)
!		call writeIntArray("indx", indx, inx)
!		call writeIntArray("indy", indy, iny)


		! Write uhat using direct access
		! Same direct access usage as in Parallel I/O

		fname = "uhat"
		CALL Write_uhat()
		call MPI_Barrier(icomm_grid,ierr)

		fname = "Ez"
		if (sampleEZ) call Write_Ez()
		call MPI_Barrier(icomm_grid,ierr)

		deallocate(uhat)
		deallocate(Ez)

		! Write flag that says Ez file was written
		if (sampleEz) call writeInt("iEz", 1)

		!--- Output for Ez (use for small grid only) ----
		!print*,"start write Ez",irank
		!		if (sampleEz) then
		!			allocate (S(nx,ny,0:niz/2))
		!			call stat_allGatherXY(Ez, S, niz/2+1)
		!			deallocate (Ez)
		!			call writeArray("Ez", S, nx*ny*(niz/2+1))
		!			deallocate (S)
		!		endif
		!print*,"end write Ez",irank

	end if

	return
end

!=======================================================================
subroutine statistics_sample()
	use statistics
        use messenger
        include "mpif.h"
	real, allocatable :: Uavg(:,:)
	real :: umid, vmid, wmid	! velocities at pressure node (center of cell C.V.)
	real :: TKEmid			! Total Kinetic Energy at pressure node (center of cell C.V.)

	if (istat == 0) return

	sampleTime = sampleTime + dt
	nstat = nstat + 1

	if (mod(nstat,istat) .ne. 0) return

	nsample = nsample + 1
	nterms = nterms + 1

	! Average in the homogeneous z direction
	ds = 1./ngzm

	!------------------------------------------------------------------
	!	i1 = ibmin ; if (iblock ==  1 ) i1 =  0  ; i1 = imap_1(i1)
	!	i2 = ibmax ; if (iblock == npx) i2 = ngx ; i2 = imap_1(i2)
	!	j1 = jbmin ; if (jblock ==  1 ) j1 =  0  ; j1 = jmap_1(j1)
	!	j2 = jbmax ; if (jblock == npy) j2 = ngy ; j2 = jmap_1(j2)
	!------------------------------------------------------------------
	i1 = ibmin ; if (iblock ==  1 ) i1 =  1  ; i1 = imap_1(i1)
	i2 = ibmax ; if (iblock == npx) i2 = ngx ; i2 = imap_1(i2)
	j1 = jbmin ; if (jblock ==  1 ) j1 =  1  ; j1 = jmap_1(j1)
	j2 = jbmax ; if (jblock == npy) j2 = ngy ; j2 = jmap_1(j2)

	do j=j1,j2
	do i=i1,i2
	do k=1,ngzm
		! Velocities at center of P cell C.V.
		umid = 0.5 * (uc(k,i,j) + uc(k,i+1,j))
		vmid = 0.5 * (vc(k,i,j) + vc(k,i,j+1))
		wmid = 0.5 * (wc(k,i,j) + wc(k+1,i,j))
		TKEmid = umid**2 + vmid**2 + wmid**2

		Ui(i,j,1) = Ui(i,j,1) + ds*uc(k,i,j)
		Ui(i,j,2) = Ui(i,j,2) + ds*vc(k,i,j)
		Ui(i,j,3) = Ui(i,j,3) + ds*wc(k,i,j)
		Ui(i,j,4) = Ui(i,j,4) + ds*p(k,i,j)
		Ui(i,j,5) = Ui(i,j,5) + ds*VT
		Ui(i,j,6) = Ui(i,j,6) + ds*p(k,i,j)**2

		!------------------------------------------------------------------
		!	! Pressure Corelations = UiP
		!	UiP(i,j,1) = UiP(i,j,1) + ds*( umid * p(k,i,j) )
		!	UiP(i,j,2) = UiP(i,j,2) + ds*( vmid * p(k,i,j) )
		!	UiP(i,j,3) = UiP(i,j,3) + ds*( wmid * p(k,i,j) )
		!------------------------------------------------------------------

		! Pressure Corelations = UiP
		UiP(i,j,1) = UiP(i,j,1) + ds*( uc(k,i,j) * 0.5*(p(k  ,i-1,j  )+p(k,i,j)) )
		UiP(i,j,2) = UiP(i,j,2) + ds*( vc(k,i,j) * 0.5*(p(k  ,i  ,j-1)+p(k,i,j)) )
		UiP(i,j,3) = UiP(i,j,3) + ds*( wc(k,i,j) * 0.5*(p(k-1,i  ,j  )+p(k,i,j)) )

		!------------------------------------------------------------------
		!	! Reynolds stresses = UiUj
		!	Rij(i,j,1) = Rij(i,j,1) + ds*( -1. *  umid**2 )
		!	Rij(i,j,2) = Rij(i,j,2) + ds*( -1. *  vmid**2 )
		!	Rij(i,j,3) = Rij(i,j,3) + ds*( -1. *  wmid**2 )
		!	Rij(i,j,4) = Rij(i,j,4) + ds*( -1. * umid*vmid)
		!	Rij(i,j,5) = Rij(i,j,5) + ds*( -1. * vmid*wmid)
		!	Rij(i,j,6) = Rij(i,j,6) + ds*( -1. * umid*wmid)
		!------------------------------------------------------------------

		! Reynolds stresses = UiUj
		Rij(i,j,1) = Rij(i,j,1) + ds*( -1. *  uc(k,i,j)**2 )
		Rij(i,j,2) = Rij(i,j,2) + ds*( -1. *  vc(k,i,j)**2 )
		Rij(i,j,3) = Rij(i,j,3) + ds*( -1. *  wc(k,i,j)**2 )
		Rij(i,j,4) = Rij(i,j,4) + ds*( -1. * 0.5*(uc(k,i,j-1)+uc(k,i,j)) * 0.5*(vc(k,i-1,j)+vc(k,i,j)) )
		Rij(i,j,5) = Rij(i,j,5) + ds*( -1. * 0.5*(vc(k-1,i,j)+vc(k,i,j)) * 0.5*(wc(k,i,j-1)+wc(k,i,j)) )
		Rij(i,j,6) = Rij(i,j,6) + ds*( -1. * 0.5*(uc(k-1,i,j)+uc(k,i,j)) * 0.5*(wc(k,i-1,j)+wc(k,i,j)) )


		! TKE transport (u*TKE, v*TKE, w*TKE) = UiUj2
		UiUj2(i,j,1) = UiUj2(i,j,1) + ds*( umid * TKEmid)
		UiUj2(i,j,2) = UiUj2(i,j,2) + ds*( vmid * TKEmid) 
		UiUj2(i,j,3) = UiUj2(i,j,3) + ds*( wmid * TKEmid)

		! Viscous stresses = 2*vt*Sij
		!LES ONLY	Tij(i,j,1) = Tij(i,j,1) + &
		!LES ONLY			ds*(+2.*VT *    dxi*(uc(k,i+1,j) - uc(k,i,j)) )
		!LES ONLY	Tij(i,j,2) = Tij(i,j,2) + &
		!LES ONLY			ds*(+2.*VT * dyi(j)*(vc(k,i,j+1) - vc(k,i,j)) )
		!LES ONLY	Tij(i,j,3) = Tij(i,j,3) + &
		!LES ONLY			ds*(+2.*VT *    dzi*(wc(k+1,i,j) - wc(k,i,j)) )
		!LES ONLY	Tij(i,j,4) = Tij(i,j,4) + &
		!LES ONLY			ds*(+0.5*( cym1(j)*(VT + VT) + cym2(j)*(VT + VT) ) &
		!LES ONLY	    			*( dymi(j)*(uc(k,i,j) - uc(k,i,j-1)) &
		!LES ONLY	      			      +dxi*(vc(k,i,j) - vc(k,i-1,j))  )      )
		!LES ONLY	Tij(i,j,5) = Tij(i,j,5) + &
		!LES ONLY			ds*(+0.5*( cym1(j)*(VT + VT) + cym2(j)*(VT + VT) ) &
		!LES ONLY	    			*( dymi(j)*(wc(k,i,j) - wc(k,i,j-1)) &
		!LES ONLY	      			      +dzi*(vc(k,i,j) - vc(k-1,i,j))  )      )
		!LES ONLY	Tij(i,j,6) = Tij(i,j,6) + &
		!LES ONLY			ds*(+0.25*( VT + VT  +VT + VT ) &
		!LES ONLY	   			*( dxi*(wc(k,i,j) - wc(k,i-1,j)) &
		!LES ONLY	      			  +dzi*(uc(k,i,j) - uc(k-1,i,j))  ) )

		!T_SC ONLY	! Scalar
		!T_SC ONLY	! BC used is such that SC = (T-Tinf)/(T0-Tinf)
		!T_SC ONLY	! I want (T-T0)/(T0-Tinf) so use 1-SC
		!T_SC ONLY	! 1: Mean (T)
		!T_SC ONLY	! 2: Variance (t^2)
		!T_SC ONLY	! 3-5: Turbulent scalar fluxes (-u_i t)
		!T_SC ONLY	! 6-7: fluctuating-pressure/temp. gradient (p dt/dx_i)
		!T_SC ONLY	! 8-9: fluctuating-temperature gradient RMS (dt/dx_i)^2
		!T_SC ONLY	SCij(i,j,1) = SCij(i,j,1) + ds*(1.-SC(i,j,k))
		!T_SC ONLY	SCij(i,j,2) = SCij(i,j,2) + ds*(1.-SC(i,j,k))**2
		!T_SC ONLY	SCij(i,j,3) = SCij(i,j,3) + ds*(-U(i,j,k)) &
		!T_SC ONLY			*0.5*((1.-SC(i+1,j,k)) + (1.-SC(i,j,k)))
		!T_SC ONLY	SCij(i,j,4) = SCij(i,j,4) + ds*(-V(i,j,k)) &
		!T_SC ONLY			*( cym1(j)*(1.-SC(i,j+1,k)) &
		!T_SC ONLY			+cym2(j)*(1.-SC(i,j,k)) )
		!T_SC ONLY	SCij(i,j,5) = SCij(i,j,5) + ds*(-W(i,j,k)) &
		!T_SC ONLY			*0.5*((1.-SC(i,j,k+1)) + (1.-SC(i,j,k)))
		!T_SC ONLY	SCij(i,j,6) = SCij(i,j,6) + ds*0.5*(P(i+1,j,k)+P(i,j,k))* &
		!T_SC ONLY			dxi*((1.-SC(i+1,j,k)) - (1.-SC(i,j,k)))
		!T_SC ONLY	SCij(i,j,7) = SCij(i,j,7) +  &
		!T_SC ONLY			ds*( cym1(j)*P(i,j+1,k) + &
		!T_SC ONLY			     cym2(j)*P(i,j,k) ) * &
		!T_SC ONLY			dymi(j)*((1.-SC(i,j+1,k)) - (1.-SC(i,j,k)))
		!T_SC ONLY	SCij(i,j,8) = SCij(i,j,8) + &
		!T_SC ONLY		ds*( dxi*((1.-SC(i+1,j,k)) - (1.-SC(i,j,k))) )**2
		!T_SC ONLY	SCij(i,j,9) = SCij(i,j,9) + &
		!T_SC ONLY		ds*( dymi(j)*((1.-SC(i,j+1,k)) - (1.-SC(i,j,k))) )**2

	end do
	end do
	end do

	! Compute Dissipation
	call Evaluate_eps()

	! Compute mean using average over span as estimate
	allocate(Uavg(0:nlx+1, 0:nly+1))
	Uavg = 0.

	i1 = ibmin ; if (iblock ==  1 ) i1 =  0  ; i1 = imap_1(i1)
	i2 = ibmax ; if (iblock == npx) i2 = ngx ; i2 = imap_1(i2)
	j1 = jbmin ; if (jblock ==  1 ) j1 =  0  ; j1 = jmap_1(j1)
	j2 = jbmax ; if (jblock == npy) j2 = ngy ; j2 = jmap_1(j2)

	do j=j1,j2
	do i=i1,i2
		do k=1,ngzm
			if (iuvw==1) Uavg(i,j) = Uavg(i,j) + uc(k,i,j)
			if (iuvw==2) Uavg(i,j) = Uavg(i,j) + vc(k,i,j)
			if (iuvw==3) Uavg(i,j) = Uavg(i,j) + wc(k,i,j)
		enddo
	end do
	end do
	Uavg = Uavg * ds

	! One dimensional energy spectra versus frequency E(omega)
	if (iuvw==1) call spectra(uc,  ngz , nlx+1,  nly , Uavg)
	if (iuvw==2) call spectra(vc,  ngz ,  nlx , nly+1, Uavg)
	if (iuvw==3) call spectra(wc, ngz+1,  nlx ,  nly , Uavg)
	deallocate(Uavg)

	! Write out uhat(omega,x,y,z) after each complete Fourier transform
	if (nterms==nfreq) then
		nft = nft + 1
		nterms = 0
		if (irank == iroot) then
			! Check dt
			if (abs ((dt*istat - Tlow/nfreq)/(Tlow/nfreq)) > .01) then
				print*,"Time step should be Tlow/nfreq = ",Tlow/nfreq
				print*,"The time step is ",dt*istat
			endif
		endif
		! Recall that uhat(-n) = uhat*(n)

		! Same direct access usage as in "Output_mod.f90"
		fname = "uhatwxyz.xxx"
		write(fname(10:12),"(i3.3)") nft
		CALL Write_uhat()
		call MPI_Barrier(icomm_grid,ierr)

		uhat = 0.	! reset uhat for next Fourier transform
	endif

	return
end

subroutine spectra(S,nk,ni,nj,Uavg)
	use statistics
	integer nk,ni,nj, npnts
	real S(0:nk, 0:ni, 0:nj)
	real fw(  ngzm  ,  inx  ,  iny  )	! u(z,x,y)
	!COMPAQ_TIGHT	real fz(0:ngzm-1,0:nlx+1,0:nly+1)	! 1D energy spectrum uhat(kz,x,y)
	real fz(0:ngzm+1,0:nlx+1,0:nly+1)	! (Extra padding needed for SGI)
	real Uavg(0:nlx+1, 0:nly+1)

	i1 = ibmin ; if (iblock ==  1 ) i1 =  0  ; i1 = imap_1(i1)
	i2 = ibmax ; if (iblock == npx) i2 = ngx ; i2 = imap_1(i2)
	j1 = jbmin ; if (jblock ==  1 ) j1 =  0  ; j1 = jmap_1(j1)
	j2 = jbmax ; if (jblock == npy) j2 = ngy ; j2 = jmap_1(j2)

	! Shift and transpose indices (only use a subset of x,y)
	do j=1,iny
	do i=1,inx
	do k=1,ngzm
		fw(k,i,j) = S(k,indx(i),indy(j)) - Uavg(indx(i),indy(j))
	end do
	end do
	end do

	! Shift and transpose indices (use entire x,y range)
	do j=j1,j2
	do i=i1,i2
	do k=1,ngzm
		fz(k-1,i,j) = S(k,i,j) - Uavg(i,j)
!		fz(k-1,i,j) = S(k,i,j)
	end do
	end do
	end do

	! Fourier transform data in t
	! Each call to spectra adds one more term to the brute force transform
	! Timestep should be set to dt = Tlow / nfreq
	! To get one complete period you need nfreq = Tlow/dt timesteps
	! cfac = -ii*2pi/Tlow
	hanning = 0.5*(1.-cos(2.*pi*(nterms-1)/nfreq))
	do j=1,iny
	do i=1,inx
	do iw=1,nfreq_
	do k=1,ngzm
		uhat(k,iw,i,j) = uhat(k,iw,i,j) + &
			fw(k,i,j) * exp(cfac*ifreq_(iw)*sampleTime) * hanning
	end do
	end do
	end do
	end do

	if (sampleEz) then
		
		!TJ: FFT each spanwise line for each ii, jj
		do i=0,nlx+1
		do j=0,nly+1
			call realFFTM(fz(0,i,j),ngzm,1,+1)
		end do
		end do
		
		!cccccccccccccccccccccccccccccccccccccccccccc	
		! Fourier transform data in z
		!npnts = (nlx+2)*(nly+2)
		!call realFFTM(fz, ngzm, npnts, +1)
		!cccccccccccccccccccccccccccccccccccccccccccc
		! 1-dimensional energy spectra in z
		! Energy = (real**2 + imag**2)
		! Ez stores from 0 wavenumber to niz/2 wavenumber
		! To plot Ez should combine energy from negative
		! wavenumbers into positive ones
		!       e.g. Ez(1:ngzm/2-1,:,:) = 2.* Ez(1:ngzm/2-1,:,:)

		!----- Intel fft ordering ------
		ds = 1.
		do j=j1,j2
		do i=i1,i2
			Ez(0,i,j) = Ez(0,i,j) &
					+ ds*(fz(0,i,j)**2)
			do k=1,ngzm/2-1
				Ez(k,i,j) = Ez(k,i,j) &
					+ds*(fz(2*k,i,j)**2 + fz(2*k+1,i,j)**2)
			end do
			Ez(ngzm/2,i,j) = Ez(ngzm/2,i,j) &
					+ ds*(fz(ngzm,i,j)**2)
		end do
		end do

		!----- Compaq fft ordering ------
		!	ds = 1.
		!	do j=j1,j2
		!	do i=i1,i2
		!		Ez(0,i,j) = Ez(0,i,j) &
		!				+ ds*(fz(0,i,j)**2)
		!		do k=1,ngzm/2-1
		!			Ez(k,i,j) = Ez(k,i,j) &
		!				+ds*(fz(k,i,j)**2 + fz(ngzm-k,i,j)**2)
		!		end do
		!		Ez(ngzm/2,i,j) = Ez(ngzm/2,i,j) &
		!				+ ds*(fz(ngzm/2,i,j)**2)
		!	end do
		!	end do
		!------------------------------
	endif

	return
end

!=======================================================================
subroutine statistics_info(iunit)
	use statistics

	call sectionTitle(iunit, "Statistics")
	if (istat == 0) then
		write (iunit, 11)
	else
		write (iunit, 12) istat, &
		                  nsample_, sampleTime_, &
		                  (nsample - nsample_), &
		                  (sampleTime - sampleTime_), &
		                  nsample, sampleTime
	end if

11  format (4x, "Statistics not computed")
12  format (4x, "Statistics level: ", i2 / &
	        / &
	        4x, "Previous number of samples      : ", i13 / &
	        4x, "Previous statistics period      : ", f13.6 / &
	        / &
	        4x, "Additional samples this run     : ", i13 / &
	        4x, "Additional time period this run : ", f13.6 / &
	        / &
	        4x, "Cumulative number of samples   : ", i13 / &
	        4x, "Cumulative statistics period   : ", f13.6)

	return
end


subroutine Evaluate_eps()
	use statistics
	real	:: voli, ds, Eval_eps

	! Average in the homogeneous z direction
	ds = 1./ngzm

	!-----------------------------------------------------------
	! Derivatives: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
	!-----------------------------------------------------------
	! Dissipation:  epsilon=\nu * average ( duidxj * duidxj )
	!-----------------------------------------------------------

	i1 = ibmin ; if (iblock ==  1 ) i1 =  1   ; i1 = imap_1(i1)
	i2 = ibmax ; if (iblock == npx) i2 = ngx-1; i2 = imap_1(i2)
	j1 = jbmin ; if (jblock ==  1 ) j1 =  1   ; j1 = jmap_1(j1)
	j2 = jbmax ; if (jblock == npy) j2 = ngy-1; j2 = jmap_1(j2)

	do j=j1,j2
		jb = jbmap_1(j)
	do i=i1,i2
		ib = ibmap_1(i)
		voli = 1./vp(ib,jb)
	do k=1,ngzm
		!------------------------ du/dx -------------------------
		dudx= voli  * (  uc(k,i+1,j)*suxix(ib+1,jb) &
				-uc(k, i ,j)*suxix(ib  ,jb) &
				+0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j+1)+uc(k,i+1,j+1))*svetax(ib,jb+1) &
				-0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j-1)+uc(k,i+1,j-1))*svetax(ib,jb  )   )

		dudy= voli  * (  uc(k,i+1,j)*suxiy(ib+1,jb) &
				-uc(k, i ,j)*suxiy( ib ,jb) &
				+0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j+1)+uc(k,i+1,j+1))*svetay(ib,jb+1) &
				-0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j-1)+uc(k,i+1,j-1))*svetay(ib, jb )   )

		dudz= voli  *	spz(ib,jb)*( 0.25*(uc( k ,i,j)+uc( k ,i+1,j)+uc(k+1,i,j)+uc(k+1,i+1,j)) &
					    -0.25*(uc(k-1,i,j)+uc(k-1,i+1,j)+uc( k ,i,j)+uc( k ,i+1,j))  )

		!------------------------ dv/dx -------------------------
		dvdx= voli  * (  vc(k,i,j+1)*svetax(ib,jb+1) &
				-vc(k,i, j )*svetax(ib, jb ) &
				+0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i+1,j)+vc(k,i+1,j+1))*suxix(ib+1,jb) &
				-0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i-1,j)+vc(k,i-1,j+1))*suxix( ib ,jb)   )

		dvdy= voli  * (  vc(k,i,j+1)*svetay(ib,jb+1) &
				-vc(k,i, j )*svetay(ib, jb ) &
				+0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i+1,j)+vc(k,i+1,j+1))*suxiy(ib+1,jb) &
				-0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i-1,j)+vc(k,i-1,j+1))*suxiy( ib ,jb)   )

		dvdz= voli  *	spz(ib,jb)*( 0.25*(vc( k ,i,j)+vc( k ,i,j+1)+vc(k+1,i,j)+vc(k+1,i,j+1)) &
					    -0.25*(vc(k-1,i,j)+vc(k-1,i,j+1)+vc( k ,i,j)+vc( k ,i,j+1))   )

		!------------------------ dw/dx -------------------------
		dwdx= voli  *  ( 0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i+1,j  )+wc(k+1,i+1,j  ))*suxix(ib+1,jb) &
				-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i-1,j  )+wc(k+1,i-1,j  ))*suxix( ib ,jb) &
				+0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j+1)+wc(k+1,i  ,j+1))*svetax(ib,jb+1) &
				-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j-1)+wc(k+1,i  ,j-1))*svetax(ib, jb )   )

		dwdy= voli  *  ( 0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i+1,j  )+wc(k+1,i+1,j  ))*suxiy(ib+1,jb) &
				-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i-1,j  )+wc(k+1,i-1,j  ))*suxiy( ib ,jb) &
				+0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j+1)+wc(k+1,i  ,j+1))*svetay(ib,jb+1) &
				-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j-1)+wc(k+1,i  ,j-1))*svetay(ib, jb )   )

		dwdz= voli  *	spz(ib,jb)*(wc(k+1,i,j)-wc(k,i,j))

		!------------------------ Calculate \epsilon -------------------------
		! dUidxj(1,:) = (/dudx,dudy,dudz/)
		! dUidxj(2,:) = (/dvdx,dvdy,dvdz/)
		! dUidxj(3,:) = (/dwdx,dwdy,dwdz/)

		Eval_eps  =  VT*( dudx*dudx+dudy*dudy+dudz*dudz &
			 	 +dvdx*dvdx+dvdy*dvdy+dvdz*dvdz &
			 	 +dwdx*dwdx+dwdy*dwdy+dwdz*dwdz  )

		! Dissipation (\epsilon)
		eps(i,j,1) = eps(i,j,1) + ds*Eval_eps

	end do
	end do
	end do

	return
end


!--------------------------------- BELOW IS BACKUP FOR COMPUTING DISSIPATION (\epsilon) ------------------
!
!	!-----------------------------------------------------------
!	! Derivatives: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
!	!-----------------------------------------------------------
!	! Dissipation:  epsilon=\nu * average ( duidxj * duidxj )
!	!-----------------------------------------------------------
!
!	! Average in the homogeneous z direction
!	ds = 1./ngzm
!
!	i1 = ibmin ; if (iblock ==  1 ) i1 =  0  ; i1 = imap_1(i1)
!	i2 = ibmax ; if (iblock == npx) i2 = ngx ; i2 = imap_1(i2)
!	j1 = jbmin ; if (jblock ==  1 ) j1 =  0  ; j1 = jmap_1(j1)
!	j2 = jbmax ; if (jblock == npy) j2 = ngy ; j2 = jmap_1(j2)
!
!	do j=j1,j2	!KALITZIN	do j=1,ny-1 ==> do j=1,ngym
!	do i=i1,i2	!KALITZIN	do i=1,nx-1 ==> do i=1,ngxm
!	do k=1,ngzm	!KALITZIN	do k=1,nz-1 ==> do k=1,ngzm
!
!		!------------------------ 1/V(i,j) -------------------------
!		voli = 1./vp(i,j)
!		!------------------------ du/dx -------------------------
!		dudx= voli  * (  uc(k,i+1,j)*suxix(i+1,j) &
!				-uc(k, i ,j)*suxix( i ,j) &
!				+0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j+1)+uc(k,i+1,j+1))*svetax(i,j+1) &
!				-0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j-1)+uc(k,i+1,j-1))*svetax(i, j )   )
!
!		dudy= voli  * (  uc(k,i+1,j)*suxiy(i+1,j) &
!				-uc(k, i ,j)*suxiy( i ,j) &
!				+0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j+1)+uc(k,i+1,j+1))*svetay(i,j+1) &
!				-0.25*(uc(k,i,j)+uc(k,i+1,j)+uc(k,i,j-1)+uc(k,i+1,j-1))*svetay(i, j )   )
!
!		dudz= voli  *	spz(i,j)*( 0.25*(uc( k ,i,j)+uc( k ,i+1,j)+uc(k+1,i,j)+uc(k+1,i+1,j)) &
!					  -0.25*(uc(k-1,i,j)+uc(k-1,i+1,j)+uc( k ,i,j)+uc( k ,i+1,j))  )
!
!		!------------------------ dv/dx -------------------------
!		dvdx= voli  * (  vc(k,i,j+1)*svetax(i,j+1) &
!				-vc(k,i, j )*svetax(i, j ) &
!				+0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i+1,j)+vc(k,i+1,j+1))*suxix(i+1,j) &
!				-0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i-1,j)+vc(k,i-1,j+1))*suxix( i ,j)   )
!
!		dvdy= voli  * (  vc(k,i,j+1)*svetay(i,j+1) &
!				-vc(k,i, j )*svetay(i, j ) &
!				+0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i+1,j)+vc(k,i+1,j+1))*suxiy(i+1,j) &
!				-0.25*(vc(k,i,j)+vc(k,i,j+1)+vc(k,i-1,j)+vc(k,i-1,j+1))*suxiy( i ,j)   )
!
!		dvdz= voli  *	spz(i,j)*( 0.25*(vc( k ,i,j)+vc( k ,i,j+1)+vc(k+1,i,j)+vc(k+1,i,j+1)) &
!					  -0.25*(vc(k-1,i,j)+vc(k-1,i,j+1)+vc( k ,i,j)+vc( k ,i,j+1))   )
!
!		!------------------------ dw/dx -------------------------
!		dwdx= voli  *  ( 0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i+1,j  )+wc(k+1,i+1,j  ))*suxix(i+1,j) &
!				-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i-1,j  )+wc(k+1,i-1,j  ))*suxix( i ,j) &
!				+0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j+1)+wc(k+1,i  ,j+1))*svetax(i,j+1) &
!				-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j-1)+wc(k+1,i  ,j-1))*svetax(i, j )   )
!
!		dwdy= voli  *  ( 0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i+1,j  )+wc(k+1,i+1,j  ))*suxiy(i+1,j) &
!				-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i-1,j  )+wc(k+1,i-1,j  ))*suxiy( i ,j) &
!				+0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j+1)+wc(k+1,i  ,j+1))*svetay(i,j+1) &
!				-0.25*(wc(k,i,j)+wc(k+1,i,j)+wc(k,i  ,j-1)+wc(k+1,i  ,j-1))*svetay(i, j )   )
!
!		dwdz= voli  *	spz(i,j)*(wc(k+1,i,j)-wc(k,i,j))
!
!		!------------------------ epsilon -------------------------
!		epsa  (i,j) = epsa  (i,j) + VT*	(dudx*dudx+dudy*dudy+dudz*dudz &
!						+dvdx*dvdx+dvdy*dvdy+dvdz*dvdz &
!						+dwdx*dwdx+dwdy*dwdy+dwdz*dwdz)
!
!	end do
!	end do
!	end do
!
!	epsat  =epsat  +epsa
!
!-------------------------------------------------------------------------------------------------------

