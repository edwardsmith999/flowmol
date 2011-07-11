!=======================================================================
! Parallel triDiagonal solvers using MPI
!
! MPI_triDiagonal
! MPI_triDiagonalM
! MPI_triDiagonalLUM
! MPI_triLUSolveM
! MPI_triPeriodicM
!

!=======================================================================
subroutine MPI_triDiagonal(a, b, c, r, n, icomm)
	include "mpif.h"
	real a(n), b(n), c(n), r(n)
	!T  include "mpif.inc"
	real s1(n), s2(n)
	allocatable sendbuf(:), recvbuf(:)

	! Get communicator info
	call MPI_comm_size (icomm, nproc, ierr)
	call MPI_comm_rank (icomm, irank, ierr)

	! Initialize boundary values
	s1(1) = a(1)
	s2(n) = c(n)
	if (irank == 0) s1(1) = 0.
	if (irank == nproc-1) s2(n) = 0.

	! Forward elimination
	! Upper boundary in s1(i)
	do i=2,n
		const = a(i)/b(i-1)
		b(i) = b(i) - c(i-1)*const
		r(i) = r(i) - r(i-1)*const
		s1(i) = -s1(i-1)*const
	end do

	! Backward elimination
	! Lower boundary in s2(i)
	do i=n-1,1,-1
		const = c(i)/b(i+1)
		r(i) = r(i) - r(i+1)*const
		s1(i) = s1(i) - s1(i+1)*const
		s2(i) = -s2(i+1)*const
	end do

	! All dependence has been shifted to the boundary elements
	! Communicate boundary values to root process
	! and solve reduced pentadiagonal system
	! Use of pentadiagonal system is more robust than the
	! reordered (removes zeros) tridiagonal system
	iroot = 0

	! Send rows of pentadiagonal system
	! (0, s1, b, 0, s2; r)
	!    (s1, 0, b, s2, 0; r)
	allocate(sendbuf(12))
	sendbuf(1) = 0.
	sendbuf(2) = s1(1)
	sendbuf(3) = b(1)
	sendbuf(4) = 0.
	sendbuf(5) = s2(1)
	sendbuf(6) = r(1)
	sendbuf(7) = s1(n)
	sendbuf(8) = 0.
	sendbuf(9) = b(n)
	sendbuf(10) = s2(n)
	sendbuf(11) = 0.
	sendbuf(12) = r(n)

	if (irank == iroot) allocate(recvbuf(12*nproc))
	call MPI_gather (sendbuf, 12, MPI_REAL8, &
	                 recvbuf, 12, MPI_REAL8, &
	                 iroot, icomm, ierr)

	if (irank == iroot) then
		! Solve reduced system
		call pentaDiagonal_1(recvbuf(1+6), 2*nproc-2)

		! Move solution to first slot
		do i=1,2*nproc
			recvbuf(i) = recvbuf(6*i)
		end do

		! Permute the order
		do i=1,nproc-1
			swap = recvbuf(2*i)
			recvbuf(2*i) = recvbuf(2*i+1)
			recvbuf(2*i+1) = swap
		end do
	end if

	! Scatter back the solution
	call MPI_scatter (recvbuf, 2, MPI_REAL8, &
	                  sendbuf, 2, MPI_REAL8, &
	                  iroot, icomm, ierr)
	if (irank == iroot) deallocate(recvbuf)

	r1 = sendbuf(1)
	r2 = sendbuf(2)

	deallocate(sendbuf)

	if (irank == 0) r1 = 0.
	if (irank == nproc-1) r2 = 0.

	do i=1,n
		r(i) = (r(i) - s1(i)*r1 - s2(i)*r2)/b(i)
	end do

	return
end

!=======================================================================
subroutine MPI_triDiagonalM(a, b, c, r, n, lot, icomm, ibc)
	include "mpif.h"
	real a(lot,n), b(lot,n), c(lot,n), r(lot,n)
	!T  include "mpif.inc"
	real const(lot)
	real s1(lot,n), s2(lot,n)
	real r1(lot), r2(lot)
	allocatable sendbuf(:), recvbuf(:,:,:)
	allocatable ngroup(:)

	! Get communicator info
	call MPI_comm_size (icomm, nproc, ierr)
	call MPI_comm_rank (icomm, irank, ierr)

	! Partition the lot
	allocate(ngroup(nproc))
	ngroup(:) = lot/nproc
	nremain = mod(lot,nproc)
	ngroup(1:nremain) = ngroup(1:nremain) + 1
	nlot = ngroup(1)

	allocate(sendbuf(nlot*12*nproc))
	allocate(recvbuf(nlot,6,2*nproc))

	! Initialize boundary values
	s1(:,1) = a(:,1)
	s2(:,n) = c(:,n)
	if (irank == 0) s1(:,1) = 0.
	if (irank == nproc-1) s2(:,n) = 0.

	! Forward elimination
	! Upper boundary in s1(i)
	do i=2,n
		const(:) = a(:,i)/b(:,i-1)
		b(:,i) = b(:,i) - c(:,i-1)*const(:)
		r(:,i) = r(:,i) - r(:,i-1)*const(:)
		s1(:,i) = -s1(:,i-1)*const(:)
	end do

	! Backward elimination
	! Lower boundary in s2(i)
	do i=n-1,1,-1
		const(:) = c(:,i)/b(:,i+1)
		r(:,i) = r(:,i) - r(:,i+1)*const(:)
		s1(:,i) = s1(:,i) - s1(:,i+1)*const(:)
		s2(:,i) = -s2(:,i+1)*const(:)
	end do

	! All dependence has been shifted to the boundary elements
	! Communicate boundary values to root process
	! and solve reduced pentadiagonal system
	! Use of pentadiagonal system is more robust than the
	! reordered (removes zeros) tridiagonal system

	! Send rows of pentadiagonal system
	! (0, s1, b, 0, s2; r)
	!    (s1, 0, b, s2, 0; r)

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = s1(k1:k2,1) ; L = L + nlot
		sendbuf(L:L+nk-1) = b(k1:k2,1)  ; L = L + nlot
		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = s2(k1:k2,1) ; L = L + nlot
		sendbuf(L:L+nk-1) = r(k1:k2,1)  ; L = L + nlot
		sendbuf(L:L+nk-1) = s1(k1:k2,n) ; L = L + nlot
		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = b(k1:k2,n)  ; L = L + nlot
		sendbuf(L:L+nk-1) = s2(k1:k2,n) ; L = L + nlot
		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = r(k1:k2,n)  ; L = L + nlot

		k1 = k2 + 1
	end do

	! Gather the boundary data
	call MPI_AllToAll (sendbuf, nlot*12, MPI_REAL8, &
	                   recvbuf, nlot*12, MPI_REAL8, &
	                   icomm, ierr)

	! Clear unused values
	nk = ngroup(irank+1)
	recvbuf(nk+1:nlot, :, :) = 0.
	recvbuf(nk+1:nlot, 3, :) = 1.

	! Solve reduced systems
	call pentaDiagonalM_1(recvbuf(1,1,2), 2*nproc-2, nlot)

	! Move solution to first slot
	do i=1,2*nproc
		recvbuf(:,i,1) = recvbuf(:,6,i)
	end do

	! Permute the order
	do i=1,nproc-1
		const(1:nlot) = recvbuf(:,2*i,1)
		recvbuf(:,2*i,1) = recvbuf(:,2*i+1,1)
		recvbuf(:,2*i+1,1) = const(1:nlot)
	end do

	! Scatter back the solution
	call MPI_AllToAll (recvbuf, nlot*2, MPI_REAL8, &
	                   sendbuf, nlot*2, MPI_REAL8, &
	                   icomm, ierr)

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		r1(k1:k2) = sendbuf(L:L+nk-1) ; L = L + nlot
		r2(k1:k2) = sendbuf(L:L+nk-1) ; L = L + nlot

		k1 = k2 + 1
	end do

	if (irank == 0) r1 = 0.
	if (irank == nproc-1) r2 = 0.

	do i=1,n
		r(:,i) = (r(:,i) - s1(:,i)*r1(:) - s2(:,i)*r2(:))/b(:,i)
	end do

!	if (ibc == 1) then
!		if (irank .ne. 0) r(:,0) = r1(:)
!		if (irank .ne. nproc-1) r(:,n+1) = r2(:)
!	end if

	deallocate(sendbuf)
	deallocate(recvbuf)
	deallocate(ngroup)

	return
end


!=======================================================================
subroutine MPI_triDiagonalLUM(a, b, c, s1, s2, rbuf, n, lot, icomm)
	include "mpif.h"
	real a(lot,n), b(lot,n), c(lot,n)
	real s1(lot,n), s2(lot,n), rbuf(lot*20)
	!T  include "mpif.inc"
	real const(lot)
	real r1(lot), r2(lot)
	allocatable sendbuf(:), recvbuf(:,:,:)
	allocatable ngroup(:)

	! Get communicator info
	call MPI_comm_size (icomm, nproc, ierr)
	call MPI_comm_rank (icomm, irank, ierr)

	! Partition the lot
	if (lot < nproc) &
		stop "MPI_triDiagonalLUM: lot must be >= nproc"
	allocate(ngroup(nproc))
	ngroup = lot/nproc
	nremain = mod(lot,nproc)
	ngroup(1:nremain) = ngroup(1:nremain) + 1
	nlot = ngroup(1)

	allocate(sendbuf(nlot*10*nproc))
	allocate(recvbuf(nlot,5,2*nproc))

	! Initialize boundary values
	s1(:,1) = a(:,1)
	s2(:,n) = c(:,n)
	if (irank == 0) s1(:,1) = 0.
	if (irank == nproc-1) s2(:,n) = 0.

	! Forward elimination
	! Upper boundary in s1(i)
	do i=2,n
		a(:,i) = a(:,i)/b(:,i-1)
		b(:,i) = b(:,i) - c(:,i-1)*a(:,i)
		s1(:,i) = -s1(:,i-1)*a(:,i)
	end do

	! Backward elimination
	! Lower boundary in s2(i)
	do i=n-1,1,-1
		c(:,i) = c(:,i)/b(:,i+1)
		s1(:,i) = s1(:,i) - s1(:,i+1)*c(:,i)
		s2(:,i) = -s2(:,i+1)*c(:,i)
	end do

	! Send rows of pentadiagonal system
	! (0, s1, b, 0, s2; r)
	!    (s1, 0, b, s2, 0; r)

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = s1(k1:k2,1) ; L = L + nlot
		sendbuf(L:L+nk-1) = b(k1:k2,1)  ; L = L + nlot
		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = s2(k1:k2,1) ; L = L + nlot
		sendbuf(L:L+nk-1) = s1(k1:k2,n) ; L = L + nlot
		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = b(k1:k2,n)  ; L = L + nlot
		sendbuf(L:L+nk-1) = s2(k1:k2,n) ; L = L + nlot
		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot

		k1 = k2 + 1
	end do

	! Gather the boundary data
	call MPI_AllToAll (sendbuf, nlot*10, MPI_REAL8, &
	                   recvbuf, nlot*10, MPI_REAL8, &
	                   icomm, ierr)

	! Clear unused values
	nk = ngroup(irank+1)
	recvbuf(nk+1:nlot, :, :) = 0.
	recvbuf(nk+1:nlot, 3, :) = 1.

	! LU decompose reduced systems
	call pentaDiagonalLUM_1(recvbuf(1,1,2), 2*nproc-2, nlot)

	! Save the result
	do L=1,nlot*10*(nproc-1)
		rbuf(L) = recvbuf(L,1,2)
	end do

	! Pre-divide diagonal
	do i=1,n
		b(:,i) = 1./b(:,i)
	end do

	deallocate(sendbuf)
	deallocate(recvbuf)
	deallocate(ngroup)

	return
end


!=======================================================================
subroutine MPI_triLUSolveM(a, b, c, s1, s2, rbuf, r, n, lot, icomm, ibc)
	include "mpif.h"
	real a(lot,n), b(lot,n), c(lot,n), r(lot,n)
	real s1(lot,n), s2(lot,n), rbuf(lot*20)
	!T  include "mpif.inc"
	real r1(lot), r2(lot)
	allocatable sendbuf(:), recvbuf(:,:)
	allocatable ngroup(:)

	! Get communicator info
	call MPI_comm_size (icomm, nproc, ierr)
	call MPI_comm_rank (icomm, irank, ierr)

	! Partition the lot
	allocate(ngroup(nproc))
	ngroup = lot/nproc
	nremain = mod(lot,nproc)
	ngroup(1:nremain) = ngroup(1:nremain) + 1
	nlot = ngroup(1)

	allocate(sendbuf(nlot*2*nproc))
	allocate(recvbuf(nlot,2*nproc))

	! Forward substitution
	! Upper boundary in s1(i)
	do i=2,n
		r(:,i) = r(:,i) - r(:,i-1)*a(:,i)
	end do

	! Backward elimination
	! Lower boundary in s2(i)
	do i=n-1,1,-1
		r(:,i) = r(:,i) - r(:,i+1)*c(:,i)
	end do

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		sendbuf(L:L+nk-1) = r(k1:k2,1)  ; L = L + nlot
		sendbuf(L:L+nk-1) = r(k1:k2,n)  ; L = L + nlot

		k1 = k2 + 1
	end do

	! Gather the boundary data
	call MPI_AllToAll (sendbuf, nlot*2, MPI_REAL8, &
	                   recvbuf, nlot*2, MPI_REAL8, &
	                   icomm, ierr)

	! Clear unused values
	nk = ngroup(irank+1)
	recvbuf(nk+1:nlot,:) = 0.

	! Solve reduced systems
	call pentaLUSolveM_1(rbuf, recvbuf(1,2), 2*nproc-2, nlot)

	! Permute the order
	do i=1,nproc-1
		r1(1:nlot) = recvbuf(:,2*i)
		recvbuf(:,2*i) = recvbuf(:,2*i+1)
		recvbuf(:,2*i+1) = r1(1:nlot)
	end do

	! Scatter back the solution
	call MPI_AllToAll (recvbuf, nlot*2, MPI_REAL8, &
	                   sendbuf, nlot*2, MPI_REAL8, &
	                   icomm, ierr)

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		r1(k1:k2) = sendbuf(L:L+nk-1) ; L = L + nlot
		r2(k1:k2) = sendbuf(L:L+nk-1) ; L = L + nlot

		k1 = k2 + 1
	end do

	if (irank == 0) r1 = 0.
	if (irank == nproc-1) r2 = 0.

	! Diagonal has been pre-divided
	do i=1,n
		r(:,i) = (r(:,i) - s1(:,i)*r1(:) - s2(:,i)*r2(:))*b(:,i)
	end do

!	if (ibc == 1) then
!		if (irank .ne. 0) r(:,0) = r1(:)
!		if (irank .ne. nproc-1) r(:,n+1) = r2(:)
!	end if

	deallocate(sendbuf)
	deallocate(recvbuf)
	deallocate(ngroup)

	return
end


!=======================================================================
subroutine MPI_triPeriodicM(a, b, c, r, n, lot, icomm)
	include "mpif.h"
	real a(lot,n), b(lot,n), c(lot,n), r(lot,n)
	!T  include "mpif.inc"
	real const(lot)
	real s1(lot,n), s2(lot,n)
	real r1(lot), r2(lot)
	allocatable sendbuf(:), recvbuf(:,:,:)
	allocatable ngroup(:)

	! Get communicator info
	call MPI_comm_size (icomm, nproc, ierr)
	call MPI_comm_rank (icomm, irank, ierr)

	! Partition the lot
	if (lot < nproc) &
		stop "MPI_triDiagonalM: lot must be >= nproc"
	allocate(ngroup(nproc))
	ngroup(:) = lot/nproc
	nremain = mod(lot,nproc)
	ngroup(1:nremain) = ngroup(1:nremain) + 1
	nlot = ngroup(1)

	allocate(sendbuf(nlot*12*nproc))
	allocate(recvbuf(nlot,6,2*nproc))

	! Initialize boundary values
	s1(:,1) = a(:,1)
	s2(:,n) = c(:,n)

	! Forward elimination
	! Upper boundary in s1(i)
	do i=2,n
		const(:) = a(:,i)/b(:,i-1)
		b(:,i) = b(:,i) - c(:,i-1)*const(:)
		r(:,i) = r(:,i) - r(:,i-1)*const(:)
		s1(:,i) = -s1(:,i-1)*const(:)
	end do

	! Backward elimination
	! Lower boundary in s2(i)
	do i=n-1,1,-1
		const(:) = c(:,i)/b(:,i+1)
		r(:,i) = r(:,i) - r(:,i+1)*const(:)
		s1(:,i) = s1(:,i) - s1(:,i+1)*const(:)
		s2(:,i) = -s2(:,i+1)*const(:)
	end do

	! All dependence has been shifted to the boundary elements
	! Communicate boundary values to root process
	! and solve reduced pentadiagonal system
	! Use of pentadiagonal system is more robust than the
	! reordered (removes zeros) tridiagonal system

	! Send rows of pentadiagonal system
	! (0, s1, b, 0, s2; r)
	!    (s1, 0, b, s2, 0; r)

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = s1(k1:k2,1) ; L = L + nlot
		sendbuf(L:L+nk-1) = b(k1:k2,1)  ; L = L + nlot
		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = s2(k1:k2,1) ; L = L + nlot
		sendbuf(L:L+nk-1) = r(k1:k2,1)  ; L = L + nlot
		sendbuf(L:L+nk-1) = s1(k1:k2,n) ; L = L + nlot
		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = b(k1:k2,n)  ; L = L + nlot
		sendbuf(L:L+nk-1) = s2(k1:k2,n) ; L = L + nlot
		sendbuf(L:L+nk-1) = 0.          ; L = L + nlot
		sendbuf(L:L+nk-1) = r(k1:k2,n)  ; L = L + nlot

		k1 = k2 + 1
	end do

	! Gather the boundary data
	call MPI_AllToAll (sendbuf, nlot*12, MPI_REAL8, &
	                   recvbuf, nlot*12, MPI_REAL8, &
	                   icomm, ierr)

	! Clear unused values
	nk = ngroup(irank+1)
	recvbuf(nk+1:nlot, :, :) = 0.
	recvbuf(nk+1:nlot, 3, :) = 1.

	! Solve reduced systems
	call pentaPeriodicM_1(recvbuf, 2*nproc, nlot)

	! Move solution to first slot
	do i=1,2*nproc
		recvbuf(:,i,1) = recvbuf(:,6,i)
	end do

	! Permute the order
	do i=1,nproc-1
		const(1:nlot) = recvbuf(:,2*i,1)
		recvbuf(:,2*i,1) = recvbuf(:,2*i+1,1)
		recvbuf(:,2*i+1,1) = const(1:nlot)
	end do
	! Don't forget end points
	const(1:nlot) = recvbuf(:,1,1)
	recvbuf(:,1,1) = recvbuf(:,2*nproc,1)
	recvbuf(:,2*nproc,1) = const(1:nlot)

	! Scatter back the solution
	call MPI_AllToAll (recvbuf, nlot*2, MPI_REAL8, &
	                   sendbuf, nlot*2, MPI_REAL8, &
	                   icomm, ierr)

	L = 1
	k1 = 1
	do igroup=1,nproc
		k2 = k1+ngroup(igroup)-1
		nk = k2-k1+1

		r1(k1:k2) = sendbuf(L:L+nk-1) ; L = L + nlot
		r2(k1:k2) = sendbuf(L:L+nk-1) ; L = L + nlot

		k1 = k2 + 1
	end do

	do i=1,n
		r(:,i) = (r(:,i) - s1(:,i)*r1(:) - s2(:,i)*r2(:))/b(:,i)
	end do

	deallocate(sendbuf)
	deallocate(recvbuf)
	deallocate(ngroup)

	return
end

