!=======================================================================
! Mesh definitions for cartesian coordinates on staggered grid
!
! mesh_init()
! mesh_write()
!
! Subclass must implement:
! subMesh_init()
! subMesh_write()
!

module mesh
	use data_export
	use mesh_export
#if USE_COUPLER
    use coupler
#endif
end module

!=======================================================================
! Mesh Definition
!
! Domain size
!     xL, yL, zL
! Grid size
!     nx, ny, nz
! Nodal coordiates
!     x(i), xm(i) 
!     y(j), ym(j) 
!     z(k), zm(k)
! Shorthand notation
!     dx, dxm, dxi, dxmi
!     dy, dym, yi, ymi, dyi, dymi
!     dz, dxm, dzi, dzmi
!
subroutine mesh_init()
	use mesh

    real(kind(0.d0)) rho_loc

	! Domain size
	call readFloat("xL", xL)
	call readFloat("yL", yL)
	call readFloat("zL", zL)

#if USE_COUPLER
    call CPL_cfd_adjust_domain(xL=xL,zL=zL,density_cfd=rho_loc)
#endif

	! Grid size specified in data module

	!-------------------------------------------------------------------
	! X-mesh
    allocate (x(imino:imaxo), xm(imino:imax))
	allocate (y (jmino:jmaxo), ym (jmino:jmax), dy (jmino:jmax), dym (jmin:jmax))
	allocate (yi(jmino:jmaxo), ymi(jmino:jmax), dyi(jmino:jmax), dymi(jmin:jmax))
	allocate (z(kmino:kmaxo), zm(kmino:kmax))
	allocate (cym1(jmino:jmaxo), cym2(jmino:jmaxo))
	allocate (ayu(jmino:jmax), ayv(jmin:jmax), ayw(jmino:jmax), ayp(jmino:jmax))

	! Uniform X-mesh !T includes halos
	! Non-periodic boundaries
	dx = xL/(imax - imin)
	dxi = 1./dx
	dxm = dx
	dxmi = dxi
	do i=imino,imaxo
		x(i) = (i-imin)*dx		!T  [x(0) = -dx]
	end do
        xm = 0.5*(x(imino:imax)+x(imin:imaxo))

	!-------------------------------------------------------------------
	! Y-mesh
	!T
	!T allocate (y (jmino:jmaxo), ym (jmino:jmax), dy (jmino:jmax), dym (jmin:jmax))
	!T allocate (yi(jmino:jmaxo), ymi(jmino:jmax), dyi(jmino:jmax), dymi(jmin:jmax))

	! Read the Y-mesh from archive
	call readArray("y", y, ny)
	call readArray("ym", ym, nym)

	! Shorthand definitions
	dy = 0.
	dyi = 0.
	dym = 0.
	dymi = 0.
	yi = 0.
	ymi = 0.

	do j=jmino,jmax
		dy(j) = y(j+1) - y(j)
		if (dy(j) .ne. 0.) dyi(j) = 1./dy(j)
	end do

	do j=jmin,jmax
		dym(j) = ym(j) - ym(j-1)
		if (dym(j) .ne. 0.) dymi(j) = 1./dym(j)
	end do

	do j=jmino,jmax
		if (y(j) .ne. 0.) yi(j) = 1./y(j)
		if (ym(j) .ne. 0.) ymi(j) = 1./ym(j)
	end do
	if (y(jmaxo) .ne. 0.) yi(jmaxo) = 1./y(jmaxo)
	! For cylindrical coordinates note that 1/y at origin is set to zero

	!-------------------------------------------------------------------
	! Z-mesh
	!T
	!T  allocate (z(kmino:kmaxo), zm(kmino:kmax))

	! Uniform Z-mesh
	! Periodic boundaries
	dz = zL/(kmax - kmin)
	dzi = 1./dz
	dzm = dz
	dzmi = dzi
	do k=kmino,kmaxo
		z(k) = (k-kmin)*dz
	end do
	zm = 0.5*(z(kmino:kmax)+z(kmin:kmaxo))

	!-------------------------------------------------------------------
	! Weights and metrics for finite-volume calculations

	! Interpolation weights
	!T
	!T  allocate (cym1(jmino:jmaxo), cym2(jmino:jmaxo))

	cym1 = 0.
	cym2 = 0.

	! Interpolation from ym to y
	do j=jmin,jmax
		if (ym(j)-ym(j-1) == 0.) cycle
		cym1(j) = (y(j)  - ym(j-1))/(ym(j) - ym(j-1))
		cym2(j) = (ym(j) - y(j)   )/(ym(j) - ym(j-1))
	end do

	! Interpolation from y to ym
	! ym grid is defined so as to be centered between y grid
	cy1 = 0.5
	cy2 = 0.5

	! Finite-volume metrics
	!T
	!T  allocate (ayu(jmino:jmax), ayv(jmin:jmax), ayw(jmino:jmax), ayp(jmino:jmax))
	
	axu = dxmi; ayu = dyi ; azu = dzi
	axv = dxi ; ayv = dymi; azv = dzi
	axw = dxi ; ayw = dyi ; azw = dzmi
	axp = dxi ; ayp = dyi ; azp = dzi

	!-------------------------------------------------------------------
	! Remap metric coefficient arrays in local indices
	!TAZ  DANGER ZONE (done)
	do j=1,nlyb+1
		cym1(j) = cym1(jbmap_1(j))
		cym2(j) = cym2(jbmap_1(j))
	end do
	do j=1,nlyb
		ayu(j) = ayu(jbmap_1(j))
		ayv(j) = ayv(jbmap_1(j))
		ayw(j) = ayw(jbmap_1(j))
		ayp(j) = ayp(jbmap_1(j))
		dy(j) = dy(jbmap_1(j))
		dym(j) = dym(jbmap_1(j))
		dyi(j) = dyi(jbmap_1(j))
		dymi(j) = dymi(jbmap_1(j))
	end do

	!-------------------------------------------------------------------
	! Subclass responsibility
	call subMesh_init()

	return
end

subroutine mesh_write()
	use mesh

	call writeFloat("xL", xL)
	call writeFloat("yL", yL)
	call writeFloat("zL", zL)
	call writeArray("y", y, ny)
	call writeArray("ym", ym, nym)

	! Subclass responsibility
	call subMesh_write()

	return
end

