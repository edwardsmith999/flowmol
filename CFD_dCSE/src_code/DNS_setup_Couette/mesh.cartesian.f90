!=======================================================================
! Mesh generation for cartesian coordinates on staggered grid
!
! mesh_define()
! mesh_staggered(xc, yc, zc)
! mesh_centered(xc, yc, zc)
! mesh_read (iunit, iformat)
! mesh_write(iunit, iformat)
! mesh_info(iunit)
!
! Subclass must implement:
! subMesh_define()
! subMesh_write(iunit, iformat)
! subMesh_info(iunit)
!

module mesh
        use mesh_parameters
	!include "mesh.inc"
       
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
!
subroutine mesh_define()
	use mesh

	! Domain size specified by subclass

	! Number of interior data points
	! (nix=ngx) , (niy=ngy), (niz=ngz)
	call readInt("nix", nix)
	call readInt("niy", niy)
	call readInt("niz", niz)
	nixyz = nix*niy*niz

	! Grid size
	nx = nix + 2
	ny = niy + 2
	nz = niz + 2
	nxyz = nx*ny*nz

        nxm = nx-1
        nym = ny-1
        nzm = nz-1

	! Define data index limits
        !T I made everything 0-->ngx+1
	imin = 1
	jmin = 1
	kmin = 1
	imax = nix
	jmax = niy
	kmax = niz

	! Some convenient definitions
	iminm1 = imin-1
	jminm1 = jmin-1
	kminm1 = kmin-1
	imaxp1 = imax+1
	jmaxp1 = jmax+1
	kmaxp1 = kmax+1

	! TJ: some additional covenient definitions (added to mesh.inc)
	iminp1 = imin+1
	jminp1 = jmin+1
	kminp1 = kmin+1
	imaxm1 = imax-1
	jmaxm1 = jmax-1
	kmaxm1 = kmax-1

	! Allocate arrays (!T (xm,ym,zm) are pressure locations)
	allocate (x(0:nix+1), xm(0:nix), y(0:niy+1), ym(0:niy), z(0:niz+1), zm(0:niz))
	allocate (xpg(nix,niy), ypg(nix,niy))

	!-------------------------------------------------------------------
	! Read in "grid.data"
				iunit = iopen()
	!(Open in mesh_read)	open (iunit, file="grid.data", form="unformatted")
				call mesh_read (iunit, 2)
	!(Open in mesh_read)	close (iclose(iunit))

	! Subclass must generate the mesh
	call subMesh_define()

	return
end

subroutine mesh_staggered(xc, yc, zc)
	use mesh
	real xc(nx,ny,nz), yc(nx,ny,nz), zc(nx,ny,nz)

	do i=1,nx
		xc(i,:,:) = x(i-1)
	end do

	do j=1,ny
		yc(:,j,:) = y(j-1)
	end do

	do k=1,nz
		zc(:,:,k) = z(k-1)
	end do

	return
end

subroutine mesh_centered(xc, yc, zc)
!T      I subtracted the one from (nx,ny,nz) for my mesh....!!!!
	use mesh
	real xc(nxm,nym,nzm), yc(nxm,nym,nzm), zc(nxm,nym,nzm)

	do i=1,nxm
		xc(i,:,:) = xm(i-1)
	end do

	do j=1,nym
		yc(:,j,:) = ym(j-1)
	end do

	do k=1,nzm
		zc(:,:,k) = zm(k-1)
	end do

	return
end

!=======================================================================
subroutine mesh_read (iunit, iformat)
	use mesh

        inquire(iolength=ilength) xpg
        iunit = iopen()
        open(iunit, file="grid.data", form="unformatted", access="direct", &
                recl=ilength, status="old", iostat=ierr)
        if (ierr .ne. 0) stop "grid.data file is required"
        read(iunit, rec=1) xpg
        read(iunit, rec=2) ypg
        close(iclose(iunit))

	return
end

!=======================================================================
subroutine mesh_write(iunit, iformat)
	use mesh
	allocatable xc(:,:,:), yc(:,:,:), zc(:,:,:)

	select case (iformat)
	case (0) ! Archive
		call writeInt("nx", nx)
		call writeInt("ny", ny)
		call writeInt("nz", nz)
		call writeFloat("xL", xL)
		call writeFloat("yL", yL)
		call writeFloat("zL", zL)

	case (1) ! Formatted
		write (iunit, *) "X-Mesh"
		write (iunit, *) iminm1, x(iminm1), x(imin)-x(iminm1), 1.
		write (iunit, *) imin, x(imin), x(imin)-x(iminm1), 1.
		do i=imin+1,imaxp1
			write (iunit, *) i, x(i), x(i)-x(i-1), &
			                 (x(i)-x(i-1))/(x(i-1)-x(i-2))
		end do

		write (iunit, *)
		write (iunit, *) "Y-Mesh"
		write (iunit, *) jminm1, y(jminm1), y(jmin)-y(jminm1), 1.
		write (iunit, *) jmin, y(jmin), y(jmin)-y(jminm1), 1.
		do j=jmin+1,jmaxp1
			write (iunit, *) j, y(j), y(j)-y(j-1), &
			                 (y(j)-y(j-1))/(y(j-1)-y(j-2))
		end do

	case (2) ! Unformatted
		!T---- Replace enclosed with direct access ----
		!T  write (iunit) nix
		!T  write (iunit) ((x(i), i=1,nix), j=1,niy)
		!T  
		!T  write (iunit) niy
		!T  write (iunit) ((y(j), i=1,nix), j=1,niy)
		!T---- Replace enclosed with direct access ----

		inquire(iolength=ilength) x(1)
		!Reopen file as direct access
		close(iunit)
		open (iunit, file="grid.data", form="unformatted", &
			access="direct", recl=ilength)

		nrec = 0
		!T  write (iunit) nix
		do j=1,niy
		do i=1,nix
			nrec = nrec+1
			write (iunit, rec=nrec) x(i)
		end do
		end do

		!T  write (iunit) niy
		do j=1,niy
		do i=1,nix
			nrec = nrec + 1
			write (iunit, rec=nrec) y(j)
		end do
		end do

	case (3) ! Unformatted PLOT3D grid file (staggered grid)
		allocate(xc(nx,ny,nz), yc(nx,ny,nz), zc(nx,ny,nz))
		call mesh_staggered(xc, yc, zc)
		write (iunit) nx, ny, nz
		write (iunit) xc
		write (iunit) yc
		write (iunit) zc
		deallocate(xc, yc, zc)

	case (4) ! Unformatted PLOT3D grid file (centered grid)
                 !T  Subtracted one from (nx,ny,nz) for my mesh....!!!!
		allocate(xc(nxm,nym,nzm), yc(nxm,nym,nzm), zc(nxm,nym,nzm))
		call mesh_centered(xc, yc, zc)
		write (iunit) nxm, nym, nzm
		write (iunit) xc
		write (iunit) yc
		write (iunit) zc
		deallocate(xc, yc, zc)

	case (5) ! Binary

        case (6) ! Formatted PLOT2D grid file
		allocate(xc(nx,ny,nz), yc(nx,ny,nz), zc(nx,ny,nz))
                call mesh_staggered(xc, yc, zc)
                write (iunit,'(2(i5))') nx, ny
                write (iunit,'(f12.6)') xc(:,:,nz/2)
                write (iunit,'(f12.6)') yc(:,1:ny-1,nz/2),yc(:,ny-1,nz/2)*1.01
		deallocate(xc, yc, zc)

        case (7) ! Formatted PLOT2D grid file (y-z plane for inlet)
		allocate(xc(nx,ny,nz), yc(nx,ny,nz), zc(nx,ny,nz))
                call mesh_staggered(xc, yc, zc)
                write (iunit,'(2(i5))') nz, ny ! so y direction is up on contour plot
		do j=jminm1,jmaxp1
		do k=kminm1,kmaxp1
                	write (iunit,'(f12.6)') zc(iminm1,j,k)
		enddo
		enddo
		do j=jminm1,jmaxp1
		do k=kminm1,kmaxp1
                	write (iunit,'(f12.6)') yc(iminm1,j,k)
		enddo
		enddo
		deallocate(xc, yc, zc)

	end select

	! Subclass responsibility
	call subMesh_write(iunit, iformat)

	return
end

subroutine mesh_info(iunit)
	use mesh

	call sectionTitle(iunit, "Mesh")
	write (iunit, 11) nx, ny, nz, &
	                  xL, yL, zL

	! Subclass responsibility
	call subMesh_info(iunit)

11  format (4x, "Mesh size   : ", 3i11 / &
	        4x, "Domain size : ", 3f11.6 /)

	return
end

