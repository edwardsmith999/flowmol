!=======================================================================
! Mesh generation for boundary-layer
!
! Uniform grid in X and Z
! Nonuniform in Y
! Staggered mesh
!
! Mesh refinement at walls
!
! subMesh_define()
! mesh_x()
! mesh_z()
! mesh_y()
! mesh_y_linear()
! mesh_y_tanh()
! mesh_y_tanh_seg()
! subMesh_write(iunit, iformat)
! subMesh_info(iunit)
!

module subMesh
	include "mesh.inc"
end module

!=======================================================================
subroutine subMesh_define()
	use subMesh
	real	:: Diff_y, MaxDiff

	! Domain size
	call readFloat("xL", xL)
	call readFloat("yL", yL)
	call readFloat("zL", zL)

	! Mesh generators
	call mesh_x()
	call mesh_y()
	call mesh_z()

	! Make sure y-mesh (at inlet) is same as ypg(1,:)
	MaxDiff = 0.
	do j=1,niy
		Diff_y = abs(y(j)-ypg(1,j))
		if (Diff_y.gt.MaxDiff)	MaxDiff = Diff_y
	end do
	print*,'Max. | (ypg(1,:) - (y(:) mesh.blayer.f90) | = ', MaxDiff

	! Define ym coordinates as centered between y points
        !T Changed by Tamer to make "ym" in inner and "y" outer
	do j=jminm1,jmax
		ym(j) = 0.5*(y(j+1) + y(j))
	end do

	!T   ! ym points coincide with y points on wall boundaries
	!T   ym(jminm1) = y(jminm1)
	!T 
	!T   ! ym is staggered on upper boundary
	!T   ym(jmaxp1) = y(jmax) + 0.5*(y(jmax) - y(jmax-1))
	!T 
	!T   ! y(jmaxp1) is never used, but set it equal to ym(jmaxp1)
	!T   y(jmaxp1) = ym(jmaxp1)

	return
end

subroutine mesh_x()
	use subMesh

	! Uniform X-mesh
	! Non-periodic boundaries  (!T changed the denominator)
	dx = xL/(imax - imin)
	do i=iminm1,imax
		x(i) = (i-imin)*dx
		xm(i) = (i-imin+0.5)*dx
	end do
	x(imaxp1) = (imaxp1-imin)*dx

	return
end

subroutine mesh_z()
	use subMesh

	! Uniform Z-mesh
	! Periodic boundaries  (!T changed the denominator)
	dz = zL/(kmax - kmin)
	do k=kminm1,kmax
		z(k) = (k-kmin)*dz
		zm(k) = (k-kmin+0.5)*dz
	end do
	z(kmaxp1) = (kmaxp1-kmin)*dz

	return
end

subroutine mesh_y()
	use subMesh

	call readInt("imeshY", imeshY)

	select case (imeshY)

	case (0) ! Uniform

		! Non-periodic boundaries
		dy = yL/(jmax - jmin)
		do j=jminm1,jmax
			y(j) = (j-jmin)*dy
		end do
                y(jmaxp1) = (jmaxp1-jmin)*dy

	case (1) ! Piecewise-linear
		y(1:niy) = ypg(1,:)
		y(  0  ) = y( 1 ) - (y( 2 )-y(  1  ))
		y(niy+1) = y(niy) + (y(niy)-y(niy-1))
		! call mesh_y_linear()

	case (2) ! Hyperbolic-tangent
		y(1:niy) = ypg(1,:)
		y(  0  ) = y( 1 ) - (y( 2 )-y(  1  ))
		y(niy+1) = y(niy) + (y(niy)-y(niy-1))
		! call mesh_y_tanh()

	case default
		stop "invalid value for imeshY"

	end select

	return
end

subroutine mesh_y_linear()
	use subMesh
	parameter(np = 2)
	real yp(np), ratio(np)

	call readFloat("ratio", ratio1)

	yp(1) = 0.
	yp(2) = yL

	ratio(1) = 1.
	ratio(2) = ratio1

       !T  I changed the number of points to be (niy) and "y -> y(jmin)"
	call linearMesh(yp, ratio, np, y(jmin), jmax-jmin+1)
        y(jminm1) = y(jmin) - (y(jmin+1)-y(jmin  ))
        y(jmaxp1) = y(jmax) + (y(jmax  )-y(jmax-1))

	return
end

subroutine mesh_y_tanh()
	use subMesh
	parameter (np = 1)
	real    yp(0:np), ratio(np)
	integer jp(0:np), iyc(2)

	call readFloat("ratio", ratio1)

	ratio(1) = ratio1

	yp(0) = 0.
	yp(1) = yL

	iyc(1) = 0
	iyc(2) = 1

       !T  I changed the number of points to be (niy) and "y -> y(jmin)"
	call tanhMesh (yp, ratio, jp, iyc, np, y(jmin), jmax-jmin+1)
        y(jminm1) = y(jmin) - (y(jmin+1)-y(jmin  ))
        y(jmaxp1) = y(jmax) + (y(jmax  )-y(jmax-1))

	return
end

subroutine mesh_y_tanh_seg()
	use subMesh
	parameter (np = 1)
	real    yp(0:np), ratio(np)
	integer jp(0:np), iyc(2)

	call readFloat("ratio", ratio1)

	ratio(1) = ratio1

	yp(0) = 0.
	yp(1) = yL

	jp(0) = jminm1
	jp(1) = jmax

	iyc(1) = 0
	iyc(2) = 1

       !T  I changed the number of points to be (niy) and "y -> y(jmin)"
	call tanhMeshSeg(yp, ratio, jp, iyc, np, y(jmin), jmax-jmin+1)
        y(jminm1) = y(jmin) - (y(jmin+1)-y(jmin  ))
        y(jmaxp1) = y(jmax) + (y(jmax  )-y(jmax-1))

	return
end

!=======================================================================
subroutine subMesh_write(iunit, iformat)
	use subMesh
	allocatable iblank(:,:,:)

	select case (iformat)
	case (0) ! Archive
                call writeArray("y", y, ny)
                call writeArray("ym", ym, nym)

	case (3,4) ! Unformatted PLOT3D grid file
		! Append iBLANK array
		allocate(iblank(nx,ny,nz))
		iblank = 1

		! Periodicity in z
		iblank(:,:,1) = -1
		iblank(:,:,nz-1) = -1
		iblank(:,:,nz) = 0

		write (iunit) iblank
		deallocate(iblank)

	end select

	return
end

subroutine subMesh_info(iunit)
	use subMesh

	! Minimum and maximum mesh spacing
	dymin = minval(y(jmin:jmax)-y(jmin-1:jmax-1))
	dymax = maxval(y(jmin:jmax)-y(jmin-1:jmax-1))

	write (iunit, 11) dx, &
	                  dymin, dymax, &
	                  dz

11  format (4x, "Mesh spacing" / &
	        8x, "    dx = ", f11.6 / &
	        8x, "min dy = ", f11.6, 4x, "max dy = ", f11.6 / &
	        8x, "    dz = ", f11.6 /)

	return
end

