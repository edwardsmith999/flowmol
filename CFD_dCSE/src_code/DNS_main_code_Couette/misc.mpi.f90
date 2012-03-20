! These are miscellaneous routines primarily for use with
! freestream.f90 and fluctuations.f90
! They are written to implement MPI on a parallel machine
!
!T	RECALL SIZES of B:
!T		allocate(uPrime(ny_1,ngz,nw/2))
!T		allocate(vPrime(ny_1,ngz,nw/2))
!T		allocate(wPrime(ny_1,ngz,nw/2))

!=======================================================================
subroutine scatterY(A, B, ni)
	use messenger
	include "mpif.h"

	complex A(ngy,ngz,ni), B(ny_1,ngz,ni)
	complex sendbuf(npy*ny_1*ngz*ni)
	integer j1(npy), j2(npy)
!	integer idims(3), id

	! Scatter a complex XY-array among the processors on one x-plane
	! including overlapping borders
	! Scattered index is first
	
	!TAZ DANGER ZONE:---------------------------------
	! This routine assumes the processes are ordered
	! in the y-direction according to (0,1,2,...,npy) 
	! on the left most boundary of the domain
	!-------------------------------------------------
	if (myid == idroot) then
		ipy = 0
		do ip=1,nproc
			if (icoord(1,ip)==1) then
				ipy = ipy + 1
				j1(ipy) = jbmin_1(icoord(2,ip))
			endif
		end do
		j2 = j1 + ny_1 - 1

		!T ----- My addition (because he starts at wall, and I have a halo)
                do ip=1,ipy
		   if (j1(ip) .eq.  1 ) j1(ip) = 1
                   if (j2(ip) .eq. ngy) j2(ip) = ngy
		end do

		L = 0
		do ip=1,npy
			do i=1,ni
			do k=1,ngz
			do j=j1(ip),j2(ip)
				L = L + 1
				sendbuf(L) = A(j,k,i)
			end do
			end do
			end do
		end do
	end if

	!TAZ DANGER ZONE:  I think the multiplicatin by two is for double precision
	call MPI_scatter (sendbuf, 2*ny_1*ngz*ni, MPI_COMPLEX, &
	                  B, 2*ny_1*ngz*ni, MPI_COMPLEX, &
	                  idroot, icomm_xyz(id_y), ierr)
	return
end

subroutine computeOmega(omega,nw)
	use messenger
	include "mpif.h"

	real omega(nw/2)

	call MPI_Bcast(omega,nw/2,MPI_REAL8,idroot,icomm_xyz(id_y),ierr)

	return
end

subroutine BCAST_TOP_EDGE(A,npnts,ixyz)
	use messenger
	include "mpif.h"

	integer	:: npnts, ixyz, root_id
	real	:: A(npnts)

	root_id = 0
	call MPI_Bcast(A,npnts,MPI_REAL8,root_id,icomm_xyz(ixyz),ierr)

	return
end

