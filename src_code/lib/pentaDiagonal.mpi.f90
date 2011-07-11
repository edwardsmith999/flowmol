!=======================================================================
! Parallel pentaDiagonal solvers using MPI
!
! Charles Pierce, January 1997
!
! MPI_pentaDiagonal
! MPI_pentaDiagonalM
! MPI_pentaPeriodicM
!

!=======================================================================
subroutine MPI_pentaDiagonal(a, b, c, d, e, r, n, icomm)
	include "mpif.h"
	real a(n), b(n), c(n), d(n), e(n), r(n)
	!T  include "mpif.inc"
	real s1(n,2), s2(n,2)
	real r1(2), r2(2), swap(2)
	allocatable sendbuf(:), recvbuf(:)
	allocatable AA(:,:)

	if (n < 1) stop "MPI_pentaDiagonal: n < 1"

	! Get communicator info
	call MPI_comm_size (icomm, nproc, ierr)
	call MPI_comm_rank (icomm, irank, ierr)

	if (n > 1) then ! do normal stuff

	! Initialize boundary values
	if (irank == 0) then
		s1(1:2,1:2) = 0.
	else
		s1(1,1) = a(1)
		s1(1,2) = b(1)
		s1(2,1) = 0.
		s1(2,2) = a(2)
	end if
	if (irank == nproc-1) then
		s2(n-1:n,1:2) = 0.
	else
		s2(n-1,1) = e(n-1)
		s2(n-1,2) = 0.
		s2(n,1) = d(n)
		s2(n,2) = e(n)
	end if

	! Forward elimination
	! Upper boundary in s1(i,1:2)
	do i=1,n-2
		! Eliminate a(i+2)
		const = a(i+2)/c(i)
		b(i+2) = b(i+2) - d(i)*const
		c(i+2) = c(i+2) - e(i)*const
		r(i+2) = r(i+2) - r(i)*const
		s1(i+2,1) = -s1(i,1)*const
		s1(i+2,2) = -s1(i,2)*const

		! Eliminate b(i+1)
		const = b(i+1)/c(i)
		c(i+1) = c(i+1) - d(i)*const
		d(i+1) = d(i+1) - e(i)*const
		r(i+1) = r(i+1) - r(i)*const
		s1(i+1,1) = s1(i+1,1) - s1(i,1)*const
		s1(i+1,2) = s1(i+1,2) - s1(i,2)*const
	end do
	! Eliminate b(n)
	const = b(n)/c(n-1)
	c(n) = c(n) - d(n-1)*const
	r(n) = r(n) - r(n-1)*const
	s1(n,1) = s1(n,1) - s1(n-1,1)*const
	s1(n,2) = s1(n,2) - s1(n-1,2)*const
	s2(n,1) = s2(n,1) - s2(n-1,1)*const

	! Backward elimination
	! Lower boundary in s2(i,1:2)
	do i=n,3,-1
		! Eliminate e(i-2)
		const = e(i-2)/c(i)
		r(i-2) = r(i-2) - r(i)*const
		s1(i-2,1) = s1(i-2,1) - s1(i,1)*const
		s1(i-2,2) = s1(i-2,2) - s1(i,2)*const
		s2(i-2,1) = -s2(i,1)*const
		s2(i-2,2) = -s2(i,2)*const

		! Eliminate d(i-1)
		const = d(i-1)/c(i)
		r(i-1) = r(i-1) - r(i)*const
		s1(i-1,1) = s1(i-1,1) - s1(i,1)*const
		s1(i-1,2) = s1(i-1,2) - s1(i,2)*const
		s2(i-1,1) = s2(i-1,1) - s2(i,1)*const
		s2(i-1,2) = s2(i-1,2) - s2(i,2)*const
	end do
	! Eliminate d(1)
	const = d(1)/c(2)
	r(1) = r(1) - r(2)*const
	s1(1,1) = s1(1,1) - s1(2,1)*const
	s1(1,2) = s1(1,2) - s1(2,2)*const
	s2(1,1) = s2(1,1) - s2(2,1)*const
	s2(1,2) = s2(1,2) - s2(2,2)*const

	end if ! n > 1

	! All dependence has been shifted to the boundary elements
	! Communicate boundary values to root process
	! and solve reduced 11-diagonal system
	! Use of 11-diagonal system is more robust than the
	! reordered (removes zeros) 7-diagonal system
	iroot = 0

	! Send rows of 11-diagonal system
	! (0, 0, 0, a, b, c, 0, 0, 0, d, e; r)
	!    (0, 0, a, b, 0, c, 0, 0, d, e, 0; r)
	!       (0, a, b, 0, 0, c, 0, d, e, 0, 0; r)
	!          (a, b, 0, 0, 0, c, d, e, 0, 0, 0; r)
	! For efficiency, only send non-zero elements
	ibufsize = 4*2*(2+1) ! 24
	allocate(sendbuf(ibufsize))

	if (n > 1) then ! do normal stuff

		L = 1
		do i=1,2
			sendbuf(L+0) = c(i)
			sendbuf(L+1) = r(i)
			sendbuf(L+2) = c(n-2+i)
			sendbuf(L+3) = r(n-2+i)
			L = L + 4
		end do
		do i=1,2
		do j=1,2
			sendbuf(L+0) = s1(i,j)
			sendbuf(L+1) = s2(i,j)
			sendbuf(L+2) = s1(n-2+i,j)
			sendbuf(L+3) = s2(n-2+i,j)
			L = L + 4
		end do
		end do

	else ! n == 1 special case

		sendbuf(1) = c(1)
		sendbuf(2) = r(1)
		sendbuf(3) = 1.
		sendbuf(4) = 0.
		sendbuf(5) = 1.
		sendbuf(6) = 0.
		sendbuf(7) = c(1)
		sendbuf(8) = r(1)
		sendbuf(9) = a(1)
		sendbuf(10) = d(1)
		sendbuf(11) = 0.
		sendbuf(12) = 0.
		sendbuf(13) = b(1)
		sendbuf(14) = e(1)
		sendbuf(15) = -1.
		sendbuf(16) = 0.
		sendbuf(17) = 0.
		sendbuf(18) = -1.
		sendbuf(19) = a(1)
		sendbuf(20) = d(1)
		sendbuf(21) = 0.
		sendbuf(22) = 0.
		sendbuf(23) = b(1)
		sendbuf(24) = e(1)

	end if

	if (irank == iroot) allocate(recvbuf(nproc*ibufsize))
	call MPI_gather (sendbuf, ibufsize, MPI_REAL8, &
	                 recvbuf, ibufsize, MPI_REAL8, &
	                 iroot, icomm, ierr)

	if (irank == iroot) then
		! Build reduced matrix
		allocate(AA(-5:5 + 1, 4*nproc))
		AA = 0.
		L = 1
		do k=1,nproc
			m = 4*(k-1)
			do i=1,2
				AA(0, m+i)   = recvbuf(L+0)         ! c(i)
				AA(6, m+i)   = recvbuf(L+1)         ! r(i)
				AA(0, m+i+2) = recvbuf(L+2)         ! c(n-2+i)
				AA(6, m+i+2) = recvbuf(L+3)         ! r(n-2+i)
				L = L + 4
			end do
			do i=1,2
			do j=1,2
				AA(-2+j-i, m+i)   = recvbuf(L+0)    ! s1(i,j)
				AA( 4+j-i, m+i)   = recvbuf(L+1)    ! s2(i,j)
				AA(-4+j-i, m+i+2) = recvbuf(L+2)    ! s1(n-2+i,j)
				AA( 2+j-i, m+i+2) = recvbuf(L+3)    ! s2(n-2+i,j)
				L = L + 4
			end do
			end do
		end do

		! Solve reduced system
		call polyDiagonal_1(5, AA(-5,3), (2*nproc-2)*2)

		! Move solution to beginning of recvbuf
		recvbuf(1:4*nproc) = AA(6,:)
		deallocate(AA)

		! Permute the order
		do i=1,nproc-1
			swap(:) = recvbuf(4*i-1:4*i)
			recvbuf(4*i-1:4*i) = recvbuf(4*i+1:4*i+2)
			recvbuf(4*i+1:4*i+2) = swap(:)
		end do
	end if

	! Scatter back the solution
	call MPI_scatter (recvbuf, 4, MPI_REAL8, &
	                  sendbuf, 4, MPI_REAL8, &
	                  iroot, icomm, ierr)
	if (irank == iroot) deallocate(recvbuf)

	r1 = sendbuf(1:2)
	r2 = sendbuf(3:4)

	deallocate(sendbuf)

	if (irank == 0) r1 = 0.
	if (irank == nproc-1) r2 = 0.

	if (n > 1) then ! do normal stuff

		do j=1,2
			r(:) = r(:) - s1(:,j)*r1(j) - s2(:,j)*r2(j)
		end do
		r = r / c

	else ! n == 1 special case

		r(1) = ( r(1) - a(1)*r1(1) - b(1)*r1(2) &
		        -d(1)*r2(1) - e(1)*r2(2) )/c(1)

	end if

	return
end

!=======================================================================
subroutine MPI_pentaDiagonalM(a, b, c, d, e, r, n, lot, icomm)
	include "mpif.h"
	real a(lot,n), b(lot,n), c(lot,n), d(lot,n), e(lot,n), r(lot,n)
	!T  include "mpif.inc"
	real const(lot)
	real s1(lot,n,2), s2(lot,n,2)
	real r1(lot,2), r2(lot,2), swap(lot,2)
	allocatable sendbuf(:,:), recvbuf(:,:)
	allocatable AA(:,:,:)
	allocatable ngroup(:)

	if (n < 1) stop "MPI_pentaDiagonalM: n < 1"

	! Get communicator info
	call MPI_comm_size (icomm, nproc, ierr)
	call MPI_comm_rank (icomm, irank, ierr)

	! Partition the lot
	if (lot < nproc) &
		stop "MPI_pentaDiagonalM: lot must be >= nproc"
	allocate(ngroup(nproc))
	ngroup(:) = lot/nproc
	nremain = mod(lot,nproc)
	ngroup(1:nremain) = ngroup(1:nremain) + 1
	nlot = ngroup(1)

	allocate(sendbuf(nlot,24*nproc))
	allocate(recvbuf(nlot,24*nproc))

	if (n > 1) then ! do normal stuff

	! Initialize boundary values
	if (irank == 0) then
		s1(:,1:2,1:2) = 0.
	else
		s1(:,1,1) = a(:,1)
		s1(:,1,2) = b(:,1)
		s1(:,2,1) = 0.
		s1(:,2,2) = a(:,2)
	end if
	if (irank == nproc-1) then
		s2(:,n-1:n,1:2) = 0.
	else
		s2(:,n-1,1) = e(:,n-1)
		s2(:,n-1,2) = 0.
		s2(:,n,1) = d(:,n)
		s2(:,n,2) = e(:,n)
	end if

	! Forward elimination
	! Upper boundary in s1(:,i,1:2)
	do i=1,n-2
		! Eliminate a(i+2)
		const(:) = a(:,i+2)/c(:,i)
		b(:,i+2) = b(:,i+2) - d(:,i)*const(:)
		c(:,i+2) = c(:,i+2) - e(:,i)*const(:)
		r(:,i+2) = r(:,i+2) - r(:,i)*const(:)
		s1(:,i+2,1) = -s1(:,i,1)*const(:)
		s1(:,i+2,2) = -s1(:,i,2)*const(:)

		! Eliminate b(i+1)
		const(:) = b(:,i+1)/c(:,i)
		c(:,i+1) = c(:,i+1) - d(:,i)*const(:)
		d(:,i+1) = d(:,i+1) - e(:,i)*const(:)
		r(:,i+1) = r(:,i+1) - r(:,i)*const(:)
		s1(:,i+1,1) = s1(:,i+1,1) - s1(:,i,1)*const(:)
		s1(:,i+1,2) = s1(:,i+1,2) - s1(:,i,2)*const(:)
	end do
	! Eliminate b(n)
	const(:) = b(:,n)/c(:,n-1)
	c(:,n) = c(:,n) - d(:,n-1)*const(:)
	r(:,n) = r(:,n) - r(:,n-1)*const(:)
	s1(:,n,1) = s1(:,n,1) - s1(:,n-1,1)*const(:)
	s1(:,n,2) = s1(:,n,2) - s1(:,n-1,2)*const(:)
	s2(:,n,1) = s2(:,n,1) - s2(:,n-1,1)*const(:)

	! Backward elimination
	! Lower boundary in s2(:,i,1:2)
	do i=n,3,-1
		! Eliminate e(i-2)
		const(:) = e(:,i-2)/c(:,i)
		r(:,i-2) = r(:,i-2) - r(:,i)*const(:)
		s1(:,i-2,1) = s1(:,i-2,1) - s1(:,i,1)*const(:)
		s1(:,i-2,2) = s1(:,i-2,2) - s1(:,i,2)*const(:)
		s2(:,i-2,1) = -s2(:,i,1)*const(:)
		s2(:,i-2,2) = -s2(:,i,2)*const(:)

		! Eliminate d(i-1)
		const(:) = d(:,i-1)/c(:,i)
		r(:,i-1) = r(:,i-1) - r(:,i)*const(:)
		s1(:,i-1,1) = s1(:,i-1,1) - s1(:,i,1)*const(:)
		s1(:,i-1,2) = s1(:,i-1,2) - s1(:,i,2)*const(:)
		s2(:,i-1,1) = s2(:,i-1,1) - s2(:,i,1)*const(:)
		s2(:,i-1,2) = s2(:,i-1,2) - s2(:,i,2)*const(:)
	end do
	! Eliminate d(1)
	const(:) = d(:,1)/c(:,2)
	r(:,1) = r(:,1) - r(:,2)*const(:)
	s1(:,1,1) = s1(:,1,1) - s1(:,2,1)*const(:)
	s1(:,1,2) = s1(:,1,2) - s1(:,2,2)*const(:)
	s2(:,1,1) = s2(:,1,1) - s2(:,2,1)*const(:)
	s2(:,1,2) = s2(:,1,2) - s2(:,2,2)*const(:)

	end if ! n > 1

	! All dependence has been shifted to the boundary elements
	! Communicate boundary values to root process
	! and solve reduced 11-diagonal system
	! Use of 11-diagonal system is more robust than the
	! reordered (removes zeros) 7-diagonal system

	! Send rows of 11-diagonal system
	! (0, 0, 0, a, b, c, 0, 0, 0, d, e; r)
	!    (0, 0, a, b, 0, c, 0, 0, d, e, 0; r)
	!       (0, a, b, 0, 0, c, 0, d, e, 0, 0; r)
	!          (a, b, 0, 0, 0, c, d, e, 0, 0, 0; r)
	! For efficiency, only send non-zero elements

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		if (n > 1) then ! do normal stuff

			do i=1,2
				sendbuf(1:nk, L+0) = c(k1:k2, i)
				sendbuf(1:nk, L+1) = r(k1:k2, i)
				sendbuf(1:nk, L+2) = c(k1:k2, n-2+i)
				sendbuf(1:nk, L+3) = r(k1:k2, n-2+i)
				L = L + 4
			end do
			do i=1,2
			do j=1,2
				sendbuf(1:nk, L+0) = s1(k1:k2, i,j)
				sendbuf(1:nk, L+1) = s2(k1:k2, i,j)
				sendbuf(1:nk, L+2) = s1(k1:k2, n-2+i,j)
				sendbuf(1:nk, L+3) = s2(k1:k2, n-2+i,j)
				L = L + 4
			end do
			end do

		else ! n == 1 special case

			sendbuf(1:nk, L+0) = c(k1:k2, 1)
			sendbuf(1:nk, L+1) = r(k1:k2, 1)
			sendbuf(1:nk, L+2) = 1.
			sendbuf(1:nk, L+3) = 0.
			sendbuf(1:nk, L+4) = 1.
			sendbuf(1:nk, L+5) = 0.
			sendbuf(1:nk, L+6) = c(k1:k2, 1)
			sendbuf(1:nk, L+7) = r(k1:k2, 1)
			sendbuf(1:nk, L+8) = a(k1:k2, 1)
			sendbuf(1:nk, L+9) = d(k1:k2, 1)
			sendbuf(1:nk, L+10) = 0.
			sendbuf(1:nk, L+11) = 0.
			sendbuf(1:nk, L+12) = b(k1:k2, 1)
			sendbuf(1:nk, L+13) = e(k1:k2, 1)
			sendbuf(1:nk, L+14) = -1.
			sendbuf(1:nk, L+15) = 0.
			sendbuf(1:nk, L+16) = 0.
			sendbuf(1:nk, L+17) = -1.
			sendbuf(1:nk, L+18) = a(k1:k2, 1)
			sendbuf(1:nk, L+19) = d(k1:k2, 1)
			sendbuf(1:nk, L+20) = 0.
			sendbuf(1:nk, L+21) = 0.
			sendbuf(1:nk, L+22) = b(k1:k2, 1)
			sendbuf(1:nk, L+23) = e(k1:k2, 1)
			L = L + 24

		end if

		k1 = k2 + 1
	end do

	! Gather the boundary data
	call MPI_AllToAll (sendbuf, nlot*24, MPI_REAL8, &
	                   recvbuf, nlot*24, MPI_REAL8, &
	                   icomm, ierr)
	deallocate(sendbuf)

	! Build reduced matrix
	allocate(AA(nlot, -5:5 + 1, 4*nproc))
	AA = 0.
	L = 1
	do k=1,nproc
		m = 4*(k-1)

		do i=1,2
			AA(:, 0, m+i)   = recvbuf(:, L+0)    ! c(i)
			AA(:, 6, m+i)   = recvbuf(:, L+1)    ! r(i)
			AA(:, 0, m+i+2) = recvbuf(:, L+2)    ! c(n-2+i)
			AA(:, 6, m+i+2) = recvbuf(:, L+3)    ! r(n-2+i)
			L = L + 4
		end do
		do i=1,2
		do j=1,2
			AA(:, -2+j-i, m+i)   = recvbuf(:, L+0)    ! s1(i,j)
			AA(:,  4+j-i, m+i)   = recvbuf(:, L+1)    ! s2(i,j)
			AA(:, -4+j-i, m+i+2) = recvbuf(:, L+2)    ! s1(n-2+i,j)
			AA(:,  2+j-i, m+i+2) = recvbuf(:, L+3)    ! s2(n-2+i,j)
			L = L + 4
		end do
		end do

	end do

	! Clear unused values
	nk = ngroup(irank+1)
	AA(nk+1:nlot, :, :) = 0.
	AA(nk+1:nlot, 0, :) = 1.

	! Solve reduced systems
	call polyDiagonalM_1(5, AA(1,-5,3), (2*nproc-2)*2, nlot)

	! Move solution to beginning of recvbuf
	recvbuf(:, 1:4*nproc) = AA(:, 6, :)
	deallocate(AA)

	! Permute the order
	do i=1,nproc-1
		swap(:,:) = recvbuf(:, 4*i-1:4*i)
		recvbuf(:, 4*i-1:4*i) = recvbuf(:, 4*i+1:4*i+2)
		recvbuf(:, 4*i+1:4*i+2) = swap(:,:)
	end do

	! Scatter back the solution
	allocate(sendbuf(nlot,4*nproc))
	call MPI_AllToAll (recvbuf, nlot*4, MPI_REAL8, &
	                   sendbuf, nlot*4, MPI_REAL8, &
	                   icomm, ierr)

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		r1(k1:k2, :) = sendbuf(1:nk, L+0:L+1)
		r2(k1:k2, :) = sendbuf(1:nk, L+2:L+3)
		L = L + 4

		k1 = k2 + 1
	end do

	if (irank == 0) r1 = 0.
	if (irank == nproc-1) r2 = 0.

	if (n > 1) then ! do normal stuff

		do j=1,2
		do i=1,n
				r(:,i) = r(:,i) - s1(:,i,j)*r1(:,j) - s2(:,i,j)*r2(:,j)
		end do
		end do
		r = r / c

	else ! n == 1 special case

		r(:,1) = ( r(:,1) - a(:,1)*r1(:,1) - b(:,1)*r1(:,2) &
		        -d(:,1)*r2(:,1) - e(:,1)*r2(:,2) )/c(:,1)

	end if

	deallocate(sendbuf)
	deallocate(recvbuf)
	deallocate(ngroup)

	return
end

!=======================================================================
subroutine MPI_pentaPeriodicM(a, b, c, d, e, r, n, lot, icomm)
	include "mpif.h"
	real a(lot,n), b(lot,n), c(lot,n), d(lot,n), e(lot,n), r(lot,n)
	!T  include "mpif.inc"
	real const(lot)
	real s1(lot,n,2), s2(lot,n,2)
	real r1(lot,2), r2(lot,2), swap(lot,2)
	allocatable sendbuf(:,:), recvbuf(:,:)
	allocatable AA(:,:,:)
	allocatable ngroup(:)

	if (n < 1) stop "MPI_pentaPeriodicM: n < 1"

	! Get communicator info
	call MPI_comm_size (icomm, nproc, ierr)
	call MPI_comm_rank (icomm, irank, ierr)

	! Partition the lot
	if (lot < nproc) &
		stop "MPI_pentaDiagonalM: lot must be >= nproc"
	allocate(ngroup(nproc))
	ngroup(:) = lot/nproc
	nremain = mod(lot,nproc)
	ngroup(1:nremain) = ngroup(1:nremain) + 1
	nlot = ngroup(1)

	allocate(sendbuf(nlot,24*nproc))
	allocate(recvbuf(nlot,24*nproc))

	if (n > 1) then ! do normal stuff

	! Initialize boundary values
	s1(:,1,1) = a(:,1)
	s1(:,1,2) = b(:,1)
	s1(:,2,1) = 0.
	s1(:,2,2) = a(:,2)
	s2(:,n-1,1) = e(:,n-1)
	s2(:,n-1,2) = 0.
	s2(:,n,1) = d(:,n)
	s2(:,n,2) = e(:,n)

	! Forward elimination
	! Upper boundary in s1(:,i,1:2)
	do i=1,n-2
		! Eliminate a(i+2)
		const(:) = a(:,i+2)/c(:,i)
		b(:,i+2) = b(:,i+2) - d(:,i)*const(:)
		c(:,i+2) = c(:,i+2) - e(:,i)*const(:)
		r(:,i+2) = r(:,i+2) - r(:,i)*const(:)
		s1(:,i+2,1) = -s1(:,i,1)*const(:)
		s1(:,i+2,2) = -s1(:,i,2)*const(:)

		! Eliminate b(i+1)
		const(:) = b(:,i+1)/c(:,i)
		c(:,i+1) = c(:,i+1) - d(:,i)*const(:)
		d(:,i+1) = d(:,i+1) - e(:,i)*const(:)
		r(:,i+1) = r(:,i+1) - r(:,i)*const(:)
		s1(:,i+1,1) = s1(:,i+1,1) - s1(:,i,1)*const(:)
		s1(:,i+1,2) = s1(:,i+1,2) - s1(:,i,2)*const(:)
	end do
	! Eliminate b(n)
	const(:) = b(:,n)/c(:,n-1)
	c(:,n) = c(:,n) - d(:,n-1)*const(:)
	r(:,n) = r(:,n) - r(:,n-1)*const(:)
	s1(:,n,1) = s1(:,n,1) - s1(:,n-1,1)*const(:)
	s1(:,n,2) = s1(:,n,2) - s1(:,n-1,2)*const(:)
	s2(:,n,1) = s2(:,n,1) - s2(:,n-1,1)*const(:)

	! Backward elimination
	! Lower boundary in s2(:,i,1:2)
	do i=n,3,-1
		! Eliminate e(i-2)
		const(:) = e(:,i-2)/c(:,i)
		r(:,i-2) = r(:,i-2) - r(:,i)*const(:)
		s1(:,i-2,1) = s1(:,i-2,1) - s1(:,i,1)*const(:)
		s1(:,i-2,2) = s1(:,i-2,2) - s1(:,i,2)*const(:)
		s2(:,i-2,1) = -s2(:,i,1)*const(:)
		s2(:,i-2,2) = -s2(:,i,2)*const(:)

		! Eliminate d(i-1)
		const(:) = d(:,i-1)/c(:,i)
		r(:,i-1) = r(:,i-1) - r(:,i)*const(:)
		s1(:,i-1,1) = s1(:,i-1,1) - s1(:,i,1)*const(:)
		s1(:,i-1,2) = s1(:,i-1,2) - s1(:,i,2)*const(:)
		s2(:,i-1,1) = s2(:,i-1,1) - s2(:,i,1)*const(:)
		s2(:,i-1,2) = s2(:,i-1,2) - s2(:,i,2)*const(:)
	end do
	! Eliminate d(1)
	const(:) = d(:,1)/c(:,2)
	r(:,1) = r(:,1) - r(:,2)*const(:)
	s1(:,1,1) = s1(:,1,1) - s1(:,2,1)*const(:)
	s1(:,1,2) = s1(:,1,2) - s1(:,2,2)*const(:)
	s2(:,1,1) = s2(:,1,1) - s2(:,2,1)*const(:)
	s2(:,1,2) = s2(:,1,2) - s2(:,2,2)*const(:)

	end if ! n > 1

	! All dependence has been shifted to the boundary elements
	! Communicate boundary values to root process
	! and solve reduced 11-diagonal system
	! Use of 11-diagonal system is more robust than the
	! reordered (removes zeros) 7-diagonal system

	! Send rows of 11-diagonal system
	! (0, 0, 0, a, b, c, 0, 0, 0, d, e; r)
	!    (0, 0, a, b, 0, c, 0, 0, d, e, 0; r)
	!       (0, a, b, 0, 0, c, 0, d, e, 0, 0; r)
	!          (a, b, 0, 0, 0, c, d, e, 0, 0, 0; r)
	! For efficiency, only send non-zero elements

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		if (n > 1) then ! do normal stuff

			do i=1,2
				sendbuf(1:nk, L+0) = c(k1:k2, i)
				sendbuf(1:nk, L+1) = r(k1:k2, i)
				sendbuf(1:nk, L+2) = c(k1:k2, n-2+i)
				sendbuf(1:nk, L+3) = r(k1:k2, n-2+i)
				L = L + 4
			end do
			do i=1,2
			do j=1,2
				sendbuf(1:nk, L+0) = s1(k1:k2, i,j)
				sendbuf(1:nk, L+1) = s2(k1:k2, i,j)
				sendbuf(1:nk, L+2) = s1(k1:k2, n-2+i,j)
				sendbuf(1:nk, L+3) = s2(k1:k2, n-2+i,j)
				L = L + 4
			end do
			end do

		else ! n == 1 special case

			sendbuf(1:nk, L+0) = c(k1:k2, 1)
			sendbuf(1:nk, L+1) = r(k1:k2, 1)
			sendbuf(1:nk, L+2) = 1.
			sendbuf(1:nk, L+3) = 0.
			sendbuf(1:nk, L+4) = 1.
			sendbuf(1:nk, L+5) = 0.
			sendbuf(1:nk, L+6) = c(k1:k2, 1)
			sendbuf(1:nk, L+7) = r(k1:k2, 1)
			sendbuf(1:nk, L+8) = a(k1:k2, 1)
			sendbuf(1:nk, L+9) = d(k1:k2, 1)
			sendbuf(1:nk, L+10) = 0.
			sendbuf(1:nk, L+11) = 0.
			sendbuf(1:nk, L+12) = b(k1:k2, 1)
			sendbuf(1:nk, L+13) = e(k1:k2, 1)
			sendbuf(1:nk, L+14) = -1.
			sendbuf(1:nk, L+15) = 0.
			sendbuf(1:nk, L+16) = 0.
			sendbuf(1:nk, L+17) = -1.
			sendbuf(1:nk, L+18) = a(k1:k2, 1)
			sendbuf(1:nk, L+19) = d(k1:k2, 1)
			sendbuf(1:nk, L+20) = 0.
			sendbuf(1:nk, L+21) = 0.
			sendbuf(1:nk, L+22) = b(k1:k2, 1)
			sendbuf(1:nk, L+23) = e(k1:k2, 1)
			L = L + 24

		end if

		k1 = k2 + 1
	end do

	! Gather the boundary data
	call MPI_AllToAll (sendbuf, nlot*24, MPI_REAL8, &
	                   recvbuf, nlot*24, MPI_REAL8, &
	                   icomm, ierr)
	deallocate(sendbuf)

	! Build reduced matrix
	allocate(AA(nlot, -5:5 + 1, 4*nproc))
	AA = 0.
	L = 1
	do k=1,nproc
		m = 4*(k-1)

		do i=1,2
			AA(:, 0, m+i)   = recvbuf(:, L+0)    ! c(i)
			AA(:, 6, m+i)   = recvbuf(:, L+1)    ! r(i)
			AA(:, 0, m+i+2) = recvbuf(:, L+2)    ! c(n-2+i)
			AA(:, 6, m+i+2) = recvbuf(:, L+3)    ! r(n-2+i)
			L = L + 4
		end do
		do i=1,2
		do j=1,2
			AA(:, -2+j-i, m+i)   = recvbuf(:, L+0)    ! s1(i,j)
			AA(:,  4+j-i, m+i)   = recvbuf(:, L+1)    ! s2(i,j)
			AA(:, -4+j-i, m+i+2) = recvbuf(:, L+2)    ! s1(n-2+i,j)
			AA(:,  2+j-i, m+i+2) = recvbuf(:, L+3)    ! s2(n-2+i,j)
			L = L + 4
		end do
		end do

	end do

	! Clear unused values
	nk = ngroup(irank+1)
	AA(nk+1:nlot, :, :) = 0.
	AA(nk+1:nlot, 0, :) = 1.

	! Solve reduced systems
	call polyPeriodicM_1(5, AA(1,-5,1), (2*nproc)*2, nlot)

	! Move solution to beginning of recvbuf
	recvbuf(:, 1:4*nproc) = AA(:, 6, :)
	deallocate(AA)

	! Permute the order
	do i=1,nproc-1
		swap(:,:) = recvbuf(:, 4*i-1:4*i)
		recvbuf(:, 4*i-1:4*i) = recvbuf(:, 4*i+1:4*i+2)
		recvbuf(:, 4*i+1:4*i+2) = swap(:,:)
	end do
	! Don't forget end points
	swap(:,:) = recvbuf(:, 4*nproc-1:4*nproc)
	recvbuf(:, 4*nproc-1:4*nproc) = recvbuf(:, 1:2)
	recvbuf(:, 1:2) = swap(:,:)

	! Scatter back the solution
	allocate(sendbuf(nlot,4*nproc))
	call MPI_AllToAll (recvbuf, nlot*4, MPI_REAL8, &
	                   sendbuf, nlot*4, MPI_REAL8, &
	                   icomm, ierr)

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		r1(k1:k2, :) = sendbuf(1:nk, L+0:L+1)
		r2(k1:k2, :) = sendbuf(1:nk, L+2:L+3)
		L = L + 4

		k1 = k2 + 1
	end do

	if (n > 1) then ! do normal stuff

		do j=1,2
		do i=1,n
				r(:,i) = r(:,i) - s1(:,i,j)*r1(:,j) - s2(:,i,j)*r2(:,j)
		end do
		end do
		r = r / c

	else ! n == 1 special case

		r(:,1) = ( r(:,1) - a(:,1)*r1(:,1) - b(:,1)*r1(:,2) &
		        -d(:,1)*r2(:,1) - e(:,1)*r2(:,2) )/c(:,1)

	end if

	deallocate(sendbuf)
	deallocate(recvbuf)
	deallocate(ngroup)

	return
end

