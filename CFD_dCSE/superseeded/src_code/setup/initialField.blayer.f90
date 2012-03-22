!=======================================================================
! Creating initial fields for boundary layer flows
!

module initialField
	include "mesh.inc"
	allocatable U(:,:,:)
	allocatable V(:,:,:)
	allocatable W(:,:,:)

	real			:: Re
end module

!=======================================================================
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
	case (4)
		STOP "Option for BLASIUS initial field is not available"
		! call initialField_blasius()
	case (5)
		! STOP "Option for TAYLOR_GREEN initial field is not available"
		  call initialField_TaylorGreen()
	end select

	call readInt("iturb", iturb)
	if (iturb .ne. 0) then
		! STOP "Option for initial fluctuations is not available"
		  call initialField_fluct()
	end if

	call initialField_boundaries()

	return
end



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



subroutine initialField_parabolic()
	use initialField

	call readFloat("ubulk", ubulk)
	call readFloat("wbulk", wbulk)

	! Streamwise flow
	y1 = y(jminm1)
	y2 = y(jmax)
	yc = y2
	yd = y2 - y1
	do j=jmin,jmax
		U(:,:,j) = 1.5*ubulk*(1. - ((ym(j) - yc)/yd)**2)
	end do

	! Wall-normal flow
	V = 0.

	! Spanwise flow
	y1 = y(jminm1)
	y2 = y(jmax)
	yc = y2
	yd = y2 - y1
	do j=jmin,jmax
		W(:,:,j) = 1.5*wbulk*(1. - ((ym(j) - yc)/yd)**2)
	end do

	return
end



!=======================================================================
subroutine initialField_TaylorGreen()
	use initialField

	do k=kminm1,kmax
	do j=jminm1,jmax
	do i=iminm1,imax
		U(k,i,j) = -cos(x(i)) * sin(ym(j))
		V(k,i,j) = sin(xm(i)) * cos(y(j))
		W(k,i,j) = 0.
	enddo
	enddo
	enddo

	return
end

!=======================================================================
subroutine initialField_fluct()
	use initialField

	! Add random fluctuations

	! Initialize random generator
	iseed = 7
	call randomSeed(iseed)

	! Fluctuation amplitude is typically 30%
	amplitude = 0.3

	do j=jmin,jmax
	do k=kmin,kmax
	do i=imin,imax
		umean = max(U(k,i,j), V(k,i,j), W(k,i,j))
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

subroutine initialField_boundaries()
	use initialField

	! Apply boundary conditions

	! Periodicity in z
	U(kminm1,:,:) = U(kmax-1,:,:)
	U( kmax ,:,:) = U( kmin ,:,:)
	V(kminm1,:,:) = V(kmax-1,:,:)
	V( kmax ,:,:) = V( kmin ,:,:)
	W(kminm1,:,:) = W(kmax-1,:,:)
	W( kmax ,:,:) = W( kmin ,:,:)
	W(kmaxp1,:,:) = W(kmin+1,:,:)

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
	real, dimension(niz-1,nix-1,niy-1) :: con

	con = 0.0
	!T  Use a record length of (niz-1)
	inquire(iolength=ilength) con(1,1,1)
	ilength = ilength * 3 * ( (nix-1)*(niy-1)*(niz-1) )

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

