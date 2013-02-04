!=========================================================================
! TJ: Initial fields for x-z periodic channel flows
!	- Orr-Sommerfeld fluctuations can be imposed on baseflow from
!	  the "osdata" file that can be created in MATLAB
!	- Simple spanwise flucations can be imposed on the baseflow
!	- Zero-mean white noise can be imposed on the baseflow for
!	  turbulent channel flows
!=========================================================================

module initialField
        use mesh_parameters
!	include "mesh.inc"
	allocatable U(:,:,:)
	allocatable V(:,:,:)
	allocatable W(:,:,:)

	real			:: Re
	integer			:: BC_bottom, BC_top
	integer			:: ifield, iOSS, iPeriodic, ierr
end module

!=========================================================================
subroutine initialField_define()
	use initialField

	allocate(U(0:niz  ,0:nix+1,0:niy  ))
	allocate(V(0:niz  ,0:nix  ,0:niy+1))
	allocate(W(0:niz+1,0:nix  ,0:niy  ))

	! Initialize variables to zero
	U = 0.
	V = 0.
	W = 0.

	! Initialize some flow parameters
	Re	= 0.
	
	! Read Reynolds here
	call readFloat("Re",Re)

	! Get type of initial field from archive
	call readInt("ifield", ifield)

	select case (ifield)
	case (1)
		! zero initial field: do nothing
	case (2)
		call initialField_uniform()
	case (3)
		! STOP "Option for PARABOLIC initial field is not available"
		call initialField_parabolic()
	case(4) 
		call initialField_dummy()
	case(5)
		call initialField_read_uy()
	case(6)
		call initialField_linear()
 	end select

	!---- Fluctuations ----
	call readInt("iturb", iturb)
	if (iturb .ne. 0) then
		! STOP "Option for initial fluctuations is not available"
		  !call initialField_fluct_v2()
		  call initialField_fluct_v3()
	end if
	
	call readInt("iOSS", iOSS)
	if (iOSS .ne. 0) then
		call initialField_OSSfluct()
	end if

	call readInt("iPeriodic", iPeriodic)
	if (iPeriodic .ne. 0) then
		call initialField_periodic_fluct()
	end if

	!---- bottom & top BCs ----
	call readInt("BC_bottom", BC_bottom)
	call readInt("BC_top",    BC_top)

	!---- perform a check of the BCs ----
	if (BC_bottom.ne.BC_top) then
		stop 'BC_top must equal BC_bottom'
	end if

	!---- stagger velocities ----
	call Stagger_initialField()

	!---- impose boundary conditions ----
	call initialField_boundaries()

	!---- FFT checks for periodic fluctuations ----
	if (iPeriodic.ne.0) then
		call FFT_check
	end if

	! Create streamwise baseflow on collocated grid    
	do j=jmin,jmax
		print'(a,i8,2(a,f10.5))', 'Init Vel field ux(y) for cell ', j, ' at ', ypg(1,j), ' is' , U(5,5,j)
	enddo

	return
end

!=======================================================================
subroutine initialField_uniform()
	use initialField

	call readFloat("ubulk", ubulk)
	call readFloat("wbulk", wbulk)

	! Streamwise flow
	do j=jmin,jmax
		U(:,:,j) = ubulk
	end do

	! Wall-normal flow
	V = 0.

	! Spanwise flow
	do j=jmin,jmax
		W(:,:,j) = wbulk
	end do

	return
end

!=======================================================================
subroutine initialField_parabolic()
	use initialField

	double precision :: deltay

	call readFloat("ubulk", ubulk)
	call readFloat("wbulk", wbulk)

	! Streamwise flow
	y1 = y(jmin)		!lower wall
	y2 = y(jmax)		!upper wall
	yc = 0.5*(y2+y1)+y1	!centre-line position
	yd = 0.5*(y2-y1)	!domain half-height

	! Create streamwise baseflow on collocated grid    
	do j=jmin,jmax
		U(:,:,j) = 1.5*ubulk*(1.0 - (ypg(1,j)-yc)**2)
		!print'(a,i8,a,f10.5)', 'Initial Velcity field ux(y) for cell ', j, ' is' , U(5,5,j)
	enddo

	! Wall-normal flow
	V = 0.

	! Spanwise flow
	y1 = y(jmin)		!lower wall
	y2 = y(jmax)		!upper wall
	yc = 0.5*(y2+y1)+y1	!centre-line position
	yd = 0.5*(y2-y1)	!domain half-height

	! Create spanwise baseflow on collocated grid
	do j=jmin,jmax
			W(:,:,j) = 1.5*wbulk*(1.0 - (ypg(1,j)-yc)**2)
	enddo

	return
end

!==========================================================================
subroutine initialField_dummy()
	use initialField

	!--- define pi ---
	twopi = 2*2.*acos(0.) 

	!--- compute dz ---
	dz = zL/dble(niz-1)

	!--- create a dummy baseflow ---
	!do k=1,niz-1
	do k=1,niz
		U(k,:,:) =    SIN(16.*twopi*(k-1)*dz/zL)
		V(k,:,:) = 2.*SIN(16.*twopi*(k-1)*dz/zL)
		W(k,:,:) =    COS(16.*twopi*(k-1)*dz/zL)
	end do

	print*,W(1:niz-1,1,1)

	!print*,MAXVAL(U),MAXVAL(V),MAXVAL(W)
	
	print*,W(:,64,64)

	return
end


!=======================================================================
subroutine initialField_read_uy()
	use initialField

	integer				:: j, ilength
	double precision	:: read_Uanaly

	iunit = iopen()
	inquire(iolength=ilength) read_Uanaly
	open(iunit, file='uy_input', status='old',form='unformatted',access="direct",recl=ilength, iostat=ierr)

	! Read Couette profile from input file    
	
	do j=jmin,jmax
		read(iunit,rec=j) read_Uanaly
		U(:,:,j) = read_Uanaly
		!print'(a,i8,a,f10.5)', 'Initial Velcity field ux(y) for cell ', j, ' is' , U(5,5,j)
	enddo

	close(iunit)

	! Wall-normal flow
	V = 0.

	! Spanwise flow
	W = 0.

	return
end



!=======================================================================
subroutine initialField_linear()
	use initialField

	double precision :: deltay

	call readFloat("uwall_top",    uwall_top)
	call readFloat("uwall_bottom", uwall_bottom)

	! Streamwise flow
	y1 = y(jmin)		!lower wall
	y2 = y(jmax)		!upper wall
	L  = y2-y1			!Channel height
	yd = 0.5*(y2-y1)	!domain half-height
	uwall_diff = uwall_top - uwall_bottom

	! Create streamwise baseflow on collocated grid    
	do j=jmin,jmax
		U(:,:,j) = uwall_diff*ypg(1,j)/L - uwall_top
		!print'(a,i8,a,f10.5)', 'Initial Velcity field ux(y) for cell ', j, ' is' , U(5,5,j)
	enddo

	! Wall-normal flow
	V = 0.d0

	! Spanwise baseflow on collocated grid
	W = 0.d0

	return
end



!=======================================================================
subroutine initialField_OSSfluct()
	use initialField

	double precision, allocatable	:: uOS_real(:), uOS_imag(:)
	double precision, allocatable	:: vOS_real(:), vOS_imag(:)
	double precision, allocatable	:: wOS_real(:), wOS_imag(:)
	double precision, allocatable	:: yOS(:)

	complex, allocatable		:: uOS(:)	
	complex, allocatable		:: vOS(:)
	complex, allocatable		:: wOS(:)	

	complex				:: ii

	double precision		:: eps_oss, alp_oss, beta_oss    !input file parameters
	double precision		:: eps_oss2, alp_oss2, beta_oss2 !input file parameters
	double precision		:: ReOS, alphaOS,  betaOS         !MATLAB file parameters
	double precision		::       alphaOS2, betaOS2         !MATLAB file parameters
	integer				:: nyOS, nmodes

	!define unit imaginary number
	ii = cmplx(0.,1.)

	!read perturbation variables from input
	call readFloat("alp_oss", alp_oss)
	call readFloat("beta_oss", beta_oss)
	call readFloat("eps_oss", eps_oss)

	!read in Reynolds number, wavenumbers and number of points in y
	iunit = iopen()
		open(iunit, file='osdata', status='old', iostat=ierr)
		if (ierr .ne. 0) Stop "The Orr-Sommerfeld data file is missing!"

		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		!	file structure:
		!
		!	line #	parameter						     type
		!	  1	Reynolds number 					    (real)
		!	  2	Alpha							    (real)
		!	  3	Beta							    (real)
		!	  4	ngym+1						 	    (integer)
		!	  5 	u_real(:) u_imag(:) v_real(:) v_imag(:) w_real(:) w_imag(:) (double)
		!
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		!
		!	NOTE: u_real( 1 ) corresponds to BOTTOM wall
		!	      u_real(niy) corresponds to   TOP  wall
		!
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		!read Reynolds number, alpha, beta and number of points in y	
		read(iunit, *) ReOS, alphaOS, betaOS, nyOS
		
		!Perform some checks on input parameters
		if (ReOS.ne.Re)         Stop "ReOS and Re must be equal"
		if (alphaOS.ne.alp_oss) Stop "alphaOS and alp_oss must be equal"	
		if (betaOS.ne.beta_oss) Stop "betaOS and beta_oss must be equal"
		if (nyOS.ne.niy)	Stop "nyOS and niy must be equal"
	
		!Parameters now match: allocate arrays for perturbation velocities
		allocate(uOS_real(nyOS),uOS_imag(nyOS))
		allocate(vOS_real(nyOS),vOS_imag(nyOS))
		allocate(wOS_real(nyOS),wOS_imag(nyOS))
	
		!initialize arrays
		uOS_real = 0. ; uOS_imag = 0.
		vOS_real = 0. ; vOS_imag = 0.
		wOS_real = 0. ; wOS_imag = 0.

		!read in normalized perturbation velocities into their arrays
		do j=1,nyOS
			read(iunit,*) uOS_real(j), uOS_imag(j), vOS_real(j), vOS_imag(j), wOS_real(j), wOS_imag(j)
		end do	
	close(iunit)

	!allocate complex arraus
	allocate(uOS(nyOS),vOS(nyOS),wOS(nyOS))

	!initialize complex arrays
	uOS = 0. ; vOS = 0. ; wOS = 0.
	
	!populate complex arrays
	uOS(:) = uOS_real(:) + ii*uOS_imag(:)
	vOS(:) = vOS_real(:) + ii*vOS_imag(:)
	wOS(:) = wOS_real(:) + ii*wOS_imag(:)

	!print perturbation parameters to screen
	print*, ' '
	print*, 'Orr-Sommerfeld mode parameters:'
	print*,'ReOS    = ',ReOS
	print*,'alphaOS = ',alphaOS
	print*,'betaOS  = ',betaOS
	print*,'nyOS    = ',nyOS
	print*, ' '

	!superimpose perturbation on to the baseflow
	do k=kmin,kmax
	do i=imin,imax
	do j=jmin,jmax
	
		U(k,i,j) = U(k,i,j) + eps_oss*REAL( uOS(j)*EXP( ii*(alp_oss*x(i)+beta_oss*z(k)) ) )
		V(k,i,j) = V(k,i,j) + eps_oss*REAL( vOS(j)*EXP( ii*(alp_oss*x(i)+beta_oss*z(k)) ) )
		W(k,i,j) = W(k,i,j) + eps_oss*REAL( wOS(j)*EXP( ii*(alp_oss*x(i)+beta_oss*z(k)) ) )
	
	end do
	end do
	end do
	
	!deallocate arrays here
	deallocate(uOS_real,uOS_imag,uOS)
	deallocate(vOS_real,vOS_imag,vOS)
	deallocate(wOS_real,wOS_imag,wOS)

	!TJ: SUPERIMPOSE 2ND DISTURBANCE IF NECCESARY
	call readInt("nmodes", nmodes)
	if (nmodes.gt.1) then

	        !read perturbation variables from input
	        call readFloat("alp_oss2", alp_oss2)
       		call readFloat("beta_oss2", beta_oss2)
        	call readFloat("eps_oss2", eps_oss2)

		print*,'Superimposing 2nd mode on initial field.'

	        !read in Reynolds number, wavenumbers and number of points in y
       		iunit = iopen()
                	open(iunit, file='osdata2', status='old', iostat=ierr)
                	if (ierr .ne. 0) Stop "The 2nd Orr-Sommerfeld data file is missing!"

                	!read Reynolds number, alpha, beta and number of points in y
                	read(iunit, *) ReOS, alphaOS2, betaOS2, nyOS
		
                	!Perform some checks on input parameters
                	if (ReOS.ne.Re)         Stop "ReOS and Re must be equal"
                	if (alphaOS2.ne.alp_oss2) Stop "alphaOS and alp_oss must be equal"
                	if (betaOS2.ne.beta_oss2) Stop "betaOS and beta_oss must be equal"
                	if (nyOS.ne.niy)        Stop "nyOS and niy must be equal"
		
                	!Parameters now match: allocate arrays for perturbation velocities
                	allocate(uOS_real(nyOS),uOS_imag(nyOS))
                	allocate(vOS_real(nyOS),vOS_imag(nyOS))
                	allocate(wOS_real(nyOS),wOS_imag(nyOS))
			
                	!initialize arrays
                	uOS_real = 0. ; uOS_imag = 0.
                	vOS_real = 0. ; vOS_imag = 0.
                	wOS_real = 0. ; wOS_imag = 0.
			
                	!read in normalized perturbation velocities into their arrays
                	do j=1,nyOS
                        	read(iunit,*) uOS_real(j), uOS_imag(j), vOS_real(j), vOS_imag(j), wOS_real(j), wOS_imag(j)
                	end do
        	close(iunit)

	        !allocate complex arraus
	        allocate(uOS(nyOS),vOS(nyOS),wOS(nyOS))

       		!initialize complex arrays
        	uOS = 0. ; vOS = 0. ; wOS = 0.
		
        	!populate complex arrays
        	uOS(:) = uOS_real(:) + ii*uOS_imag(:)
        	vOS(:) = vOS_real(:) + ii*vOS_imag(:)
        	wOS(:) = wOS_real(:) + ii*wOS_imag(:)
		
        	!print perturbation parameters to screen
        	print*, ' '
        	print*, 'Orr-Sommerfeld mode parameters:'
        	print*,'ReOS     = ',ReOS
        	print*,'alphaOS2 = ',alphaOS2
        	print*,'betaOS2  = ',betaOS2
        	print*,'nyOS     = ',nyOS
        	print*, ' '
	
        	!superimpose perturbation on to the baseflow
        	do k=kmin,kmax
        	do i=imin,imax
        	do j=jmin,jmax

                	U(k,i,j) = U(k,i,j) + eps_oss2*REAL( uOS(j)*EXP( ii*(alp_oss2*x(i)+beta_oss2*z(k)) ) )
                	V(k,i,j) = V(k,i,j) + eps_oss2*REAL( vOS(j)*EXP( ii*(alp_oss2*x(i)+beta_oss2*z(k)) ) )
                	W(k,i,j) = W(k,i,j) + eps_oss2*REAL( wOS(j)*EXP( ii*(alp_oss2*x(i)+beta_oss2*z(k)) ) )

	        end do
        	end do
        	end do
	
		!deallocate arrays here
		deallocate(uOS_real,uOS_imag,uOS)
		deallocate(vOS_real,vOS_imag,vOS)
		deallocate(wOS_real,wOS_imag,wOS)

	end if

        !TJ: SUPERIMPOSE 3RD DISTURBANCE IF NECCESARY
        call readInt("nmodes", nmodes)
        if (nmodes.eq.3) then

                !read perturbation variables from input
                call readFloat("alp_oss3", alp_oss3)
                call readFloat("beta_oss3", beta_oss3)
                call readFloat("eps_oss3", eps_oss3)

                print*,'Superimposing 3nd mode on initial field.'

                !read in Reynolds number, wavenumbers and number of points in y
                iunit = iopen()
                        open(iunit, file='osdata3', status='old', iostat=ierr)
                        if (ierr .ne. 0) Stop "The 2nd Orr-Sommerfeld data file is missing!"

                        !read Reynolds number, alpha, beta and number of points in y
                        read(iunit, *) ReOS, alphaOS3, betaOS3, nyOS

                        !Perform some checks on input parameters
                        if (ReOS.ne.Re)         Stop "ReOS and Re must be equal"
                        if (alphaOS3.ne.alp_oss3) Stop "alphaOS and alp_oss must be equal"
                        if (betaOS3.ne.beta_oss3) Stop "betaOS and beta_oss must be equal"
                        if (nyOS.ne.niy)        Stop "nyOS and niy must be equal"

                        !Parameters now match: allocate arrays for perturbation velocities
                        allocate(uOS_real(nyOS),uOS_imag(nyOS))
                        allocate(vOS_real(nyOS),vOS_imag(nyOS))
                        allocate(wOS_real(nyOS),wOS_imag(nyOS))

                        !initialize arrays
                        uOS_real = 0. ; uOS_imag = 0.
                        vOS_real = 0. ; vOS_imag = 0.
                        wOS_real = 0. ; wOS_imag = 0.

                        !read in normalized perturbation velocities into their arrays
                        do j=1,nyOS
                                read(iunit,*) uOS_real(j), uOS_imag(j), vOS_real(j), vOS_imag(j), wOS_real(j), wOS_imag(j)
                        end do
                close(iunit)

                !allocate complex arraus
                allocate(uOS(nyOS),vOS(nyOS),wOS(nyOS))

                !initialize complex arrays
                uOS = 0. ; vOS = 0. ; wOS = 0.

                !populate complex arrays
                uOS(:) = uOS_real(:) + ii*uOS_imag(:)
                vOS(:) = vOS_real(:) + ii*vOS_imag(:)
                wOS(:) = wOS_real(:) + ii*wOS_imag(:)

                !print perturbation parameters to screen
                print*, ' '
                print*, 'Orr-Sommerfeld mode parameters:'
                print*,'ReOS     = ',ReOS
                print*,'alphaOS3 = ',alphaOS3
                print*,'betaOS3  = ',betaOS3
                print*,'nyOS     = ',nyOS
                print*, ' '

                !superimpose perturbation on to the baseflow
                do k=kmin,kmax
                do i=imin,imax
                do j=jmin,jmax

                        U(k,i,j) = U(k,i,j) + eps_oss3*REAL( uOS(j)*EXP( ii*(alp_oss3*x(i)+beta_oss3*z(k)) ) )
                        V(k,i,j) = V(k,i,j) + eps_oss3*REAL( vOS(j)*EXP( ii*(alp_oss3*x(i)+beta_oss3*z(k)) ) )
                        W(k,i,j) = W(k,i,j) + eps_oss3*REAL( wOS(j)*EXP( ii*(alp_oss3*x(i)+beta_oss3*z(k)) ) )

                end do
                end do
                end do

                !deallocate arrays here
                deallocate(uOS_real,uOS_imag,uOS)
                deallocate(vOS_real,vOS_imag,vOS)
                deallocate(wOS_real,wOS_imag,wOS)

        end if

	U(:,:,jmax+1) = -U(:,:,jmax-1)
	U(:,:,  0   ) = -U(:,:,   2  )

	return
end

!=======================================================================
subroutine initialField_periodic_fluct()
	use initialField

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! TJ: baseflow consists of a poiseuille flow plus a periodic
	!     perturbation
	! 
	!	 U(y) = U_pois(y) + eps*U_pois(y)*cos(twopi*alp *x/xL)
	!
	!				OR
	!
	!	 U(y) = U_pois(y) + eps*U_pois(y)*cos(twopi*beta*z/zL)
	!
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	double precision	:: xfluc(imin:imax), zfluc(kmin:kmax)
	double precision	:: yxfluc(imin:imax,jmin:jmax)
	double precision	:: eps_x, eps_z, twopi
	integer			:: alp, beta

	!define twopi
	twopi = 2.*2.*acos(0.) 

	!read perturbation variables from input
	call readInt("alp", alp)
	call readInt("beta", beta)
	call readFloat("eps_x", eps_x)
	call readFloat("eps_z", eps_z)

	!initialize perturbation arrays
	xfluc  = 0.
	yfluc  = 0.
	yxfluc = 0.

	!add a periodic perturbation to the baseflow	
	do j=jmin,jmax
	do k=kmin,kmax
	do i=imin,imax
		
		xfluc(i)    = cos(twopi*alp* ((x(i)-x(1))/xL))	
		zfluc(k)    = cos(twopi*beta*((z(k)-z(1))/zL))

	        !yxfluc(i,j) = cos(twopi*((y(j)-y(1))/yL))*cos(twopi*alp*((x(i)-x(1))/xL))
	        yxfluc(i,j) = cos(twopi*alp*((x(i)-x(1))/xL))

		SELECT CASE(BC_BOTTOM)
		CASE(0)
			U(k,i,j) = U(k,i,j) + eps_x*U(k,i,j)*xfluc(i)	
			U(k,i,j) = U(k,i,j) + eps_z*U(k,i,j)*zfluc(k)
			W(k,i,j) = W(k,i,j) + eps_x*W(k,i,j)*xfluc(i)
			W(k,i,j) = W(k,i,j) + eps_z*W(k,i,j)*zfluc(k)
		CASE(1)
			W(k,i,j) = W(k,i,j) + eps_x*yxfluc(i,j)
		END SELECT

	end do
	end do
	end do

	return
end

!=======================================================================
subroutine initialField_fluct_v3()
	use initialField

	        ! Add random fluctuations

        ! Initialize random generator
        iseed = 25
        call randomSeed(iseed)

        ! Fluctuation amplitude is typically 30%
        amplitude = 0.8

        do j=jmin,jmax
        do k=kmin,kmax
        do i=imin,imax
                umean = max(abs(U(k,i,j)), abs(V(k,i,j)), abs(W(k,i,j)))
                fluct1 = umean*amplitude*2.*(random()-0.5)
                fluct2 = umean*amplitude*2.*(random()-0.5)
                fluct3 = umean*amplitude*2.*(random()-0.5)
                U(k,i,j) = U(k,i,j) + fluct1
                V(k,i,j) = V(k,i,j) + fluct2
                W(k,i,j) = W(k,i,j) + fluct3
        end do
        end do
        end do

	return
end
!=======================================================================
subroutine initialField_fluct_v2()
	use initialField

	double precision :: fluct1(kmin:kmax,imin:imax,jmin:jmax)
	double precision :: fluct2(kmin:kmax,imin:imax,jmin:jmax)
	double precision :: fluct3(kmin:kmax,imin:imax,jmin:jmax)

	double precision :: fluct1_mean23(imax) !the i'th y-z mean 
	double precision :: fluct2_mean13(jmax) !the j'th x-z mean
	double precision :: fluct3_mean12(kmax) !the k-th x-y mean

	double precision :: deltax, deltay, deltaz, deltaS

	fluct1 = 0.
	fluct2 = 0.
	fluct3 = 0.

	! Initialize random generator
	iseed = 7
	call randomSeed(iseed)

	! Fluctuation amplitude is typically 30%
	amplitude = 0.3

	! TJ: generate the random fluctuations here that satisfy the BCs
	do k=kmin,kmax
	do i=imin,imax
	do j=jmin,jmax

		umean = max(U(k,i,j), V(k,i,j), W(k,i,j))

		fluct1(k,i,j) = umean*amplitude*2.*(random()-0.5)
		fluct2(k,i,j) = umean*amplitude*2.*(random()-0.5)
		fluct3(k,i,j) = umean*amplitude*2.*(random()-0.5)

	end do
	end do
	end do

	! TJ: remove net mass flux in x-direction
	fluct1_mean23 = 0.
	do i=imin,imax
	do j=jmin,jmax
	do k=kmin,kmax
		deltay = y(j+1)-y(j)
		deltaz = z(k+1)-z(k)
		deltaS = deltay*deltaz
		fluct1_mean23(i) = fluct1_mean23(i) + fluct1(k,i,j)*deltaS		
	end do
	end do
	end do

	fluct1_mean23 = (1./DBLE(jmax*kmax))*(1./(yL*zL))*fluct1_mean23

	do i=imin,imax
	do k=kmin,kmax
	do j=jmin,jmax	
		fluct1(k,i,j) = fluct1(k,i,j) - fluct1_mean23(i)	
	end do
	end do
	end do

	! TJ: remove net mass flux in y-direction
	fluct2_mean13 = 0.
	do j=jmin,jmax
	do i=imin,imax
	do k=kmin,kmax
		deltax = x(i+1)-x(i)
		deltaz = z(k+1)-z(k)
		deltaS = deltax*deltaz
		fluct2_mean13(j) = fluct2_mean13(j) + fluct2(k,i,j)*deltaS
	end do
	end do
	end do

	fluct2_mean13 = (1./DBLE(imax*kmax))*(1./(xL*zL))*fluct2_mean13

	do j=jmin,jmax
	do i=imin,imax
	do k=kmin,kmax
		fluct2(k,i,j) = fluct2(k,i,j) - fluct2_mean13(j)
	end do
	end do
	end do

	! TJ: remove net mass flux in z-direction
	fluct3_mean12 = 0.
	do k=kmin,kmax
	do i=imin,imax
	do j=jmin,jmax
		deltax = x(i+1)-x(i)
		deltay = y(j+1)-y(j)
		deltaS = deltax*deltay
		fluct3_mean12(k) = fluct3_mean12(k) + fluct3(k,i,j)*deltaS
	end do
	end do
	end do

	fluct3_mean12 = (1./DBLE(imax*jmax))*(1./(xL*yL))*fluct3_mean12

	do k=kmin,kmax
	do i=imin,imax
	do j=jmin,jmax
		fluct3(k,i,j) = fluct3(k,i,j) - fluct3_mean12(k)
	end do
	end do
	end do

	! TJ: add my random noise to the base-flow
	do k=kmin,kmax
	do i=imin,imax
	do j=jmin,jmax

               U(k,i,j) = U(k,i,j) + fluct1(k,i,j)
               V(k,i,j) = V(k,i,j) + fluct2(k,i,j)
               W(k,i,j) = W(k,i,j) + fluct3(k,i,j)

	end do
	end do
	end do

	return
end
!=======================================================================
subroutine initialField_fluct()
	use initialField

	! TJ: declare 3D fluctuation array
	double precision	:: fluct_123(kmin:kmax,imin:imax,jmin:jmax)

	! TJ: declare mean-related variables
	double precision	:: fluct_mean1, fluct_mean3
	double precision	:: fluct_mean12,fluct_mean13,fluct_mean23

	! TJ: initialize the fluctuation ararys here
	fluct_123 = 0.

	! Initialize random generator
	iseed = 7
	call randomSeed(iseed)

	! Fluctuation amplitude is typically 30%
	amplitude = 0.3

	! TJ: generate the random fluctuations here that satisfy the BCs
	do k=kmin,kmax
	do i=imin,imax
	do j=jmin,jmax

		umean = max(U(k,i,j), V(k,i,j), W(k,i,j))

		fluct_123(k,i,j) = umean*amplitude*2.*(random()-0.5)

	end do
	end do
	end do

	! TJ: subtract mean value in x-direction
	do j=jmin,jmax	
	do k=kmin,kmax
		fluct_mean1 = 0.
	do i=imin,imax 
		fluct_mean1 = fluct_mean1 + fluct_123(k,i,j) !summation in i-direction
		!multiply fluct_123 by dy(j)*dz(k)
	end do
		fluct_mean1      = (1./DBLE(imax))*fluct_mean1
		!divive by (Ly*Lz) ---> for fluct_1 only
		fluct_123(k,:,j) = fluct_123(k,:,j)-fluct_mean1
	end do
	end do

	! TJ: subtract mean value in z-direction
	do j=jmin,jmax
	do i=imin,imax
		fluct_mean3 = 0.
	do k=kmin,kmax
		fluct_mean3 = fluct_mean3 + fluct_123(k,i,j)	!summation in k-direction
		!multiply fluct_123 by dx(i)*dy(j)
	end do
		fluct_mean3      = (1./DBLE(kmax))*fluct_mean3
		fluct_123(:,i,j) = fluct_123(:,i,j)-fluct_mean3
	end do
	end do

	! TJ: is the mean of fluct_123 across the z-y plane zero for every i?
	do i=imin,imax
		
		fluct_mean23 = SUM(fluct_123(i,:,:))
		fluct_mean23 = (1./DBLE(jmax*kmax))*fluct_mean23
		if (ABS(fluct_mean23).gt.1e-16) Stop 'the mean of fluct_123 across the y-z plane is non-zero'
		!print*,'(1./DBLE(jmax*kmax))*SUM(fluct_123(i,:,:) = ', fluct_mean23

	end do
		
	! TJ: is the mean of fluct_123 across the x-z plane zero for every j?
	do j=jmin,jmax

		fluct_mean13 = SUM(fluct_123(:,:,j))
		fluct_mean13 = (1./DBLE(imax*kmax))*fluct_mean13
		if (ABS(fluct_mean13).gt.1e-16) Stop 'the mean of fluct_123 across the x-z plane is non-zero'
		!print*,'(1./DBLE(imax*kmax))*SUM(fluct_123(:,:,j) = ', fluct_mean13

	end do

	! TJ: is the mean of fluct_123 across the x-y plane zero for every k?
	do k=kmin,kmax
		fluct_mean12 = SUM(fluct_123(k,:,:))
		fluct_mean12 = (1./DBLE(imax*jmax))*fluct_mean12
		if (ABS(fluct_mean12).gt.1e-16) Stop 'the mean of fluct_123 across the x-y plane is non-zero'
		!print*,'(1./DBLE(imax*jmax))*SUM(fluct_123(k,:,:) = ', fluct_mean12
	end do

	! TJ: add fluct1, fluct2, fluct3 to base-flow
	do j=jmin,jmax
		U(:,:,j) = U(:,:,j) + fluct_123(:,:,j)
		V(:,:,j) = V(:,:,j) + fluct_123(:,:,j)
		W(:,:,j) = W(:,:,j) + fluct_123(:,:,j)
	end do

	

	return
end

!=======================================================================
subroutine initialField_write(iunit)
	use initialField

	!T  Use a record length of all (u+w+w)
	inquire(iolength=ilength) U(1,1,1)

	ilength = ilength * (	  (nix+2)*(niy+1)*(niz+1) &
				+ (nix+1)*(niy+2)*(niz+1) &
				+ (nix+1)*(niy+1)*(niz+2)  )

	! Reopen file as direct access
	close(iunit)
	open (iunit, file="ucvcwc.dble.000000", form="unformatted", &
		access="direct",recl=ilength)

	ni = nix
	nj = niy
	nk = niz
	nvar = 3

	write (iunit,rec = 1) U,V,W

	return
end

!=======================================================================
subroutine initialMeanField_write(iunit)
	use initialField
	
	inquire(iolength=ilength) U(1,1,1)
	
        ilength = ilength * ( (nix+2)*(niy+1)*(niz+1) )

	! Reopen as direct access
	close(iunit)
	open (iunit, file="ucmean.dble.000000", form="unformatted", &
		access="direct",recl=ilength)
	
	write(iunit,rec=1) U(:,:,:)

	return
end

!=======================================================================
subroutine initialMeanUy_write(iunit)
        use initialField

        inquire(iolength=ilength) U(1,1,1)

        ilength = ilength * (niy+1)

        ! Reopen as direct access
        close(iunit)
        open (iunit, file="uymean.dble.000000", form="unformatted", &
                access="direct",recl=ilength)

        write(iunit,rec=1) U(1,1,:)

        return
end

!=======================================================================
subroutine initial_Flux_write(iunit)
	use initialField

	!T  Use a record length of (u+w+w)
	inquire(iolength=ilength) U(1,1,1)
	ilength = ilength * (	  (nix+2)*(niy+1)*(niz+1) &
				+ (nix+1)*(niy+2)*(niz+1) &
				+ (nix+1)*(niy+1)*(niz+2)  )

	! Reopen file as direct access
	close(iunit)
	open (iunit, file="uuvvww.dble.000000", form="unformatted", &
		access="direct",recl=ilength)


	U = 0.0
	V = 0.0
	W = 0.0

	write (iunit,rec = 1) U,V,W

	return
end

!=======================================================================
subroutine initial_Con_write(iunit)
	use initialField
	!real, dimension(niz-1,nix-1,niy-1) :: con
	
	real, dimension(niz-1+1,nix-1+1,niy-1+1) :: con

	con = 0.0
	!TJ:  Use a record length of (niz-1+1)
	inquire(iolength=ilength) con(1,1,1)
	ilength = ilength * 3 * ( (nix-1+1)*(niy-1+1)*(niz-1+1) )

	! Reopen file as direct access
	close(iunit)
	open (iunit, file="conold.dble.000000", form="unformatted", &
		access="direct",recl=ilength)

	write (iunit,rec = 1) con, con, con

	return
end


!=======================================================================
subroutine initial_phatr_write(iunit)
        use initialField
        real, dimension(0:nix, 0:niy , niz-1) :: phatr

        phatr = 0.0
        !T  Use a record length of (nix+1)*(niy+1)*(niz-1)
        inquire(iolength=ilength) phatr(:,:,:)

        ! Reopen file as direct access
        close(iunit)
        open (iunit, file="pressure_ph.000000", form="unformatted", &
                access="direct",recl=ilength)

        write (iunit,rec = 1) phatr

        return
end

!=======================================================================

subroutine initialField_info(iunit)
	use initialField

	return
end


!=======================================================================

subroutine Stagger_initialField()

        use initialField

        !--- Stagger U in y-direction ---
        do i=1,nix
                U(:,i,1:niy-1) = 0.5 * (U(:,i,1:niy-1) + U(:,i,2:niy))
        end do

        !--- Stagger V in x-direction ---
        do j=1,niy
                V(:,1:nix-1,j) = 0.5 * (V(:,1:nix-1,j) + V(:,2:nix,j))
        end do
        
	!--- Stagger W in x,y-directions ---
	do k=1,niz
		W(k,1:nix-1,1:niy-1) = 0.25*(W(k,1:nix-1,1:niy-1) + W(k,2:nix,1:niy-1)  &
 				           + W(k,1:nix-1, 2:niy ) + W(k,2:nix, 2:niy ))
	end do

        return

end

!=======================================================================

subroutine Inflow_BC_write(iunit)
        use initialField

        real UatLeftBC  (0:niy+1,3)     ! Left boundary Mean velocity (U,V,W)
        real UppatLeftBC(0:niy+1)       ! Inlet Mean velocity (U'')

        !--- copy inflow velocity ---
        UatLeftBC(0:niy  ,1) = U(1,1,0:niy  )
        UatLeftBC(0:niy+1,2) = V(1,1,0:niy+1)
        UatLeftBC(0:niy  ,3) = 0.0
       
	!--- extend ---
        UatLeftBC(  niy+1,1) = UatLeftBC(  niy  ,1)
        UatLeftBC(  niy+1,3) = UatLeftBC(  niy  ,3)

        !------------------------------------------------
        !--- compute U" from UatLeftBC ------------------
        !This formula won't work for a non-uniform grid--
        !UppatLeftBC(1:niy-1) = (UatLeftBC(2:niy,1)-2.0*UatLeftBC(1:niy-1,1)+UatLeftBC(0:niy-2,1)) / ((ypg(1,2:niy)-ypg(1,1:niy-1))**2)
        !------------------------------------------------

        !--- compute UppatLeftBC ---
        call compute_UppatLeftBC(UatLeftBC(:,1),UppatLeftBC)

        UppatLeftBC(0      ) = 2.0*UppatLeftBC(1) - 1.0*UppatLeftBC(2)
        UppatLeftBC(  niy  ) = UppatLeftBC( niy-1)
        UppatLeftBC(  niy+1) = UppatLeftBC( niy-1)

        !-- Compute record length ---
        inquire(iolength=ilength) UatLeftBC(:,1)
        ilength = ilength * 4

        ! Reopen file as direct access
        close(iunit)
        open(iunit,file='Inflow.mean.000000',form='unformatted',access='direct',recl=ilength)
        write(iunit,rec=1) UatLeftBC(:,1), UatLeftBC(:,2), UatLeftBC(:,3), UppatLeftBC(:)

        return
end

!=======================================================================

subroutine compute_UppatLeftBC(UatLeftBC,UppatLeftBC)

        use initialField

        double precision, dimension(1:niy-1)    :: delta_y
        double precision                        :: C1, C2, C3, C4, C5
        double precision, dimension(0:niy+1)    :: UppatLeftBC
        double precision, dimension(0:niy+1)    :: UatLeftBC

        !--- read dns mesh ---
        iunit = iopen()
        call mesh_read(iunit,2)

        !--- compute non-uniform delta_y from the dns grid ---
        do j=1,niy-1
                delta_y(j) = ypg(imin,j+1)-ypg(imin,j)
        end do

        !-----------------------------------------------------------------------------------
        ! 2nd Order Central differencing scheme (CDS) on a non-uniform grid:
        !-----------------------------------------------------------------------------------
        !
        !  d2U_{j}             2                          dy(j-1)        dy(j-1)
        !  -------  = --------------------- *[U(j-1)-(1 + -------)*U(j)+(-------)*U(j+1)]
        !    dy2      (dy(j)**2+dy(j-1)**2)                dy(j)          dy(j)
        !
        ! OR:
        !
        !  d2U_{j}
        !  -------  = C1*[U(j-1)-C2**U(j)+C3*U(j+1)]
        !    dy2
        !
        ! where:
        !                    2
        !       C1 = ---------------------
        !             dy(j)**2+dy(j-1)**2
        !
        !                  dy(j-1)
        !       C2 = 1 + -----------
        !                   dy(j)
        !
        !             dy(j-1)
        !       C3 = ---------
        !              dy(j)
        !
        !-----------------------------------------------------------------------------------
        ! For a constant dy, the CDS reduces to the standard form:
        !-----------------------------------------------------------------------------------
        !
        !  d2U_{j}       1
        !  -------  = --------*[U(j-1)-2*U(j)+U(j+1)]
        !    dy2      (dy**2)
        !
        !-----------------------------------------------------------------------------------

        !--- central-differencing scheme for interior points ---
        do j=2,jmax-1

                !--- coefficients account for non-uniform wall-normal grid ---
                C1 = 2.0/(delta_y(j)**2+delta_y(j-1)**2)
                C2 = 1 + (delta_y(j-1)/delta_y(j))
                C3 =      delta_y(j-1)/delta_y(j)

                !--- compute UppatLeftBC across interior points ---
                UppatLeftBC(j) = C1*(UatLeftBC(j-1) - C2*UatLeftBC(j) + C3*UatLeftBC(j+1))

        enddo

        !-----------------------------------------------------------------------------------
        ! 2nd Order, One-sided, backward difference scheme (BDS) on a non-uniform grid:
        !-----------------------------------------------------------------------------------
        !
        !  d2U_{wall}             2                           dy(2)        dy(2)
        !  ----------  = ----------------------- *[U(3)-(1 + ------)*U(2)+(-----)*U(1)]
        !     dy2        dy(1)*dy(2) + dy(2)**2               dy(1)        dy(1)
        !
        !
        ! OR:
        !
        !   d2U_{wall}
        !  ----------  = C4*[U(3)-C5*U(2)+(C5-1)*U(1)]
        !     dy2
        !
        ! where:
        !                     2
        !       C4 = ----------------------
        !            dy(1)*dy(2) + dy(2)**2
        !
        !                 dy(2)
        !       C5 = 1 + -------
        !                 dy(1)
        !
        !-----------------------------------------------------------------------------------
        ! For a constant dy, the BDS reduces to the standard form:
        !-----------------------------------------------------------------------------------
        !
        !  d2U_{wall}       1
        !  ----------  = -------*[U(3)-2*U(2)+U(1)]
        !     dy2        (dy**2)
        !
        !-----------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------
        !TJ: This isn't needed for a ZPG BL because Upp_wall = 0.0
        !    But Upp_wall can be non-zero for PG / slipping profiles!
        !
        !       C4 =     delta_y(1)*delta_y(2) + delta_y(2)**2
        !       C5 = 1 +(delta_y(2)/delta_y(1))
        !
        !       UppatLeftBC(jmin) = (2.0/C4)*(UatLeftBC(3)-C5*UatLeftBC(2)+(C5-1)*UatLeftBC(1))
        !--------------------------------------------------------------------------------------

        !--- upp is zero at the wall for ZPG, Lslip=0.0 ---
        UppatLeftBC(jmin) = 0.0

        !--- upp is zero in the free-stream ---
        UppatLeftBC(jmax) = 0.0

        return

end subroutine

!=======================================================================

subroutine initialField_boundaries()
        use initialField

        !---- Wall BCs ----
	V(:,:, jmin ) =  0.0
	V(:,:, jmax ) =  0.0

        !---- Wall BCs ----
	SELECT CASE(BC_bottom)
	CASE(0)
	        U(:,:,jminm1) = -U(:,:, jmin )
		V(:,:,jminm1) =  V(:,:,jminp1)
	        W(:,:,jminm1) = -W(:,:, jmin ) 
	CASE(1)
        	U(:,:,jminm1) =  U(:,:, jmin )
		V(:,:,jminm1) =  V(:,:,jminp1)
        	W(:,:,jminm1) =  W(:,:, jmin )
	END SELECT

        !---- Wall BCs ----
	SELECT CASE(BC_top)
	CASE(0)
	        U(:,:, jmax ) = -U(:,:,jmaxm1)
		V(:,:,jmaxp1) =  V(:,:,jmaxm1)
	        W(:,:, jmax ) = -W(:,:,jmaxm1)
	CASE(1)
	        U(:,:, jmax ) =  U(:,:,jmaxm1)
		V(:,:,jmaxp1) =  V(:,:,jmaxm1)
	        W(:,:, jmax ) =  W(:,:,jmaxm1)
	END SELECT	

        !V(:,:,jminm1) = -V(:,:,jminp1)
	!V(:,:,jmaxp1) = -V(:,:,jmaxm1)

	!---- Periodicity in x ----
	U(:,iminm1,:) = U(:,imaxm1,:) 	
	U(:,imaxp1,:) = U(:,iminp1,:)    
	V(:,iminm1,:) = V(:,imaxm1,:)   
	V(:, imax, :) = V(:, imin, :)   
	W(:,iminm1,:) = W(:,imaxm1,:)   
	W(:, imax, :) = W(:, imin, :)  

        !---- Watch out: The following is important ----
	!U(:, imax, :) = 0.5*(U(:, imax, :) + U(:, imin, :))
	U(:, imin, :) =      U(:, imax, :)

        !---- Periodicity in z ----
        U(kminm1,:,:) = U(kmaxm1,:,:) 
        U( kmax ,:,:) = U( kmin ,:,:)
        V(kminm1,:,:) = V(kmaxm1,:,:) 
        V( kmax ,:,:) = V( kmin ,:,:)
        W(kminm1,:,:) = W(kmaxm1,:,:) 
        W(kmaxp1,:,:) = W(kminp1,:,:)

        !---- Watch out: The following is important ----
        !W( kmax ,:,:) = 0.5 * (W(kmax,:,:) + W(kmin,:,:))
        W( kmin ,:,:) =        W(kmax,:,:)

        return
end

!=======================================================================
subroutine FFT_check
	use initialField

	double precision	:: array_fftx(0:imax-1+2-1) !CCS storage
	double precision	:: array_fftz(0:kmax-1+2-1) !CCS storage
	integer			:: ipick, jpick, kpick

	!---- define ipick, jpick, kpick ----
	ipick = nix/2
	jpick = niy/2
	kpick = niz/2

	!---- initialize arrays here ----
	array_fftx = 0.
	array_fftz = 0.

	!---- populate arrays here ----
	array_fftx(0:imaxm1-1) = U(kpick,1:nix-1,jpick)
	array_fftz(0:kmaxm1-1) = U(1:niz-1,ipick,jpick)

		open(11,file='U_xfluc.dat')
			do i=0,imax-2
				write(11,*)array_fftx(i)
			end do
		close(11)
		
		open(12,file='U_zfluc.dat')
			do k=0,kmax-2
				write(12,*)array_fftz(k)
			end do
		close(12)

	!---- Forward FFT ----
	call realFFTM(array_fftx,nix-1,1,+1)
	call realFFTM(array_fftz,niz-1,1,+1)
		
		open(13,file='U_xfluc_fft.dat')
			do i=0,imax-2
				write(13,*)array_fftx(i)
			end do
		close(13)
		
		open(14,file='U_zfluc_fft.dat')
			do k=0,kmax-2
				write(14,*)array_fftz(k)
			end do
		close(14)

	!---- Inverse FFT -----
	call realFFTM(array_fftx,nix-1,1,-1)
	call realFFTM(array_fftz,niz-1,1,-1)

		open(15,file='U_xfluc_ifft.dat')
			do i=0,imax-2
				write(15,*)array_fftx(i)
			end do
		close(15)
		
		open(16,file='U_zfluc_ifft.dat')
			do k=0,kmax-2
				write(16,*)array_fftz(k)
			end do
		close(16)
	
	return
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

