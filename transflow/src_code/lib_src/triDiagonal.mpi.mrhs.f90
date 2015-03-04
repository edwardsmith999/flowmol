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
subroutine MPI_triDiag_mrhs(a, b, c, r, nrhs, n, icomm)
	include "mpif.h"
	
        real*8 a(n), b(n), c(n), r(nrhs,n)
        real*8 s1(n), s2(n),const
        real*8 r1(nrhs), r2(nrhs)
        real*8, allocatable, dimension(:) :: sendbuf,sendbuf2, recvbuf
        real*8, allocatable :: recvbuf2(:,:)
        real*8, allocatable :: swap(:)
        integer, allocatable:: ngroup(:)
        integer :: nproc,icomm,irank,n,ierr,ibuftot
        integer :: nrhs, L, nlot, igroup, nk, nremain
        integer :: i, k1, k2

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
                b(i)   = b(i)   - c(i-1)  *const
                r(:,i) = r(:,i) - r(:,i-1)*const
                s1(i)  = -s1(i-1)*const
        end do

        ! Backward elimination
        ! Lower boundary in s2(i)
        do i=n-1,1,-1
                const = c(i)/b(i+1)
                r(:,i) = r(:,i) - r(:,i+1)*const
                s1(i)  = s1(i)  - s1(i+1) *const
                s2(i)  = -s2(i+1)*const
        end do

        ! All dependence has been shifted to the boundary elements
        ! Communicate boundary values to root processes
        ! and solve reduced pentadiagonal system
        ! Use of pentadiagonal system is more robust than the
        ! reordered (removes zeros) tridiagonal system
        iroot = 0

        ! Scatter the Coefficients of pentadiagonal system
        ! (0, s1, b, 0, s2)
        !    (s1, 0, b, s2, 0)
        allocate(sendbuf(10))
        sendbuf(1) = 0.
        sendbuf(2) = s1(1)
        sendbuf(3) = b(1)
        sendbuf(4) = 0.
        sendbuf(5) = s2(1)

        sendbuf(6) = s1(n)
        sendbuf(7) = 0.
        sendbuf(8) = b(n)
        sendbuf(9) = s2(n)
        sendbuf(10) = 0.

        ibuftot=10*nproc
        allocate(recvbuf(ibuftot))
        !OR  call MPI_Allgather (sendbuf, 10, MPI_DOUBLE_PRECISION, &
        !OR                      recvbuf, 10, MPI_DOUBLE_PRECISION, &
        !OR                      icomm, ierr)
        call MPI_gather (sendbuf, 10, MPI_DOUBLE_PRECISION, &
                            recvbuf, 10, MPI_DOUBLE_PRECISION, &
                            0,icomm, ierr)
        call MPI_BCAST(recvbuf,ibuftot,MPI_DOUBLE_PRECISION,0,icomm,ierr)

        ! Scatter the RHSides
        !--------------------
        ! Partition the nrhs
        ! Each processor will be responsible for (nrhs/nproc)
        ! The first couple of processors might have more because mod(nrhs/nproc)
        ! 'ngroup' has the number of lines in a group
        ! 'nlot'   has the maximum number of lines in a group

        allocate(ngroup(nproc))
        ngroup(:) = nrhs/nproc
        nremain = mod(nrhs,nproc)
        ngroup(1:nremain) = ngroup(1:nremain) + 1
        nlot = ngroup(1)

        ! Each process sends (2 rhs /equation)*(nlot equations)*(nproc in the comm)
        allocate(sendbuf2(nlot*2*nproc))
        allocate(recvbuf2(nlot,2*nproc))

        ! (K1,K2) are the limits of the equations index
        ! 'K1' lower equation number in the nrhs
        ! 'K2' Upper equation number in the nrhs
        ! 'sendbuf' has (nlot rhs1, nlot rhs2 ) for process 1
        !           followed by (nlot rhs1, nlot rhs1, ...) for process 2

        L = 1
        k1 = 1
        do igroup=1,nproc
                k2 = k1+ngroup(igroup)-1
                nk = k2-k1+1

                sendbuf2(L:L+nk-1) = r(k1:k2,1)  ; L = L + nlot
                sendbuf2(L:L+nk-1) = r(k1:k2,n)  ; L = L + nlot

                k1 = k2 + 1
        end do

        ! Gather the boundary data
        call MPI_AllToAll (sendbuf2, nlot*2, MPI_REAL8, &
                           recvbuf2, nlot*2, MPI_REAL8, &
                           icomm, ierr)

        ! Clear unused values
        nk = ngroup(irank+1)
        recvbuf2(nk+1:nlot,:) = 0.

        call pentaDiagonal_1_Mrhs(recvbuf(5+1),recvbuf2(1,2),2*nproc-2,nlot)
        recvbuf2(nk+1:nlot,:) = 0.

        allocate(swap(1:nlot))
        ! Permute the order
        do i=1,nproc-1
                swap(1:nlot)    = recvbuf2(:,2*i)
                recvbuf2(:,2*i) = recvbuf2(:,2*i+1)
                recvbuf2(:,2*i+1) = swap(1:nlot)
        end do
        deallocate(swap)

        ! Scatter back the solution
        call MPI_AllToAll(recvbuf2, nlot*2, MPI_REAL8, &
                          sendbuf2, nlot*2, MPI_REAL8, &
                          icomm, ierr)

        L = 1
        k1 = 1
        do igroup=1,nproc
                k2 = k1+ngroup(igroup)-1
                nk = k2-k1+1

                r1(k1:k2) = sendbuf2(L:L+nk-1) ; L = L + nlot
                r2(k1:k2) = sendbuf2(L:L+nk-1) ; L = L + nlot

                k1 = k2 + 1
        end do

        if (irank == 0) r1 = 0.
        if (irank == nproc-1) r2 = 0.

        do i=1,n
                r(:,i) = (r(:,i) - s1(i)*r1(:) - s2(i)*r2(:))/b(i)
        end do

!        if (ibc == 1) then
!                if (irank .ne. 0) r(:,0) = r1(:)
!                if (irank .ne. nproc-1) r(:,n+1) = r2(:)
!        end if

        deallocate(sendbuf)
        deallocate(recvbuf)
        deallocate(sendbuf2)
        deallocate(recvbuf2)
        deallocate(ngroup)

        return
end


!=======================================================================


!=======================================================================
subroutine MPI_triDiagM_mrhs(a, b, c, r, n, lot, nrhs, icomm, ibc)
	include "mpif.h"

        integer :: n, lot, nrhs, icomm, ibc, ierr, nlot
        integer :: nproc, irank, nremain, igroup
        integer :: L, k1,k2, i, nkr, M
        real*8 a(lot,n), b(lot,n), c(lot,n), r(nrhs,lot,n)
        real*8 const(lot)
        real*8 s1(lot,n), s2(lot,n)
        real*8 r1(nrhs,lot), r2(nrhs,lot)
        real*8, allocatable :: sendbuf(:) , recvbuf(:,:,:)
        real*8, allocatable :: sendbuf2(:), recvbuf2(:,:,:)
        real*8, allocatable :: swap(:,:)
        real*8, allocatable :: ngroup(:)

        ! Get communicator info
        call MPI_comm_size (icomm, nproc, ierr)
        call MPI_comm_rank (icomm, irank, ierr)

        ! Partition the lot
        ! Each processor will be responsible for (lot/nproc)
        ! The first couple of processors might have more because mod(lot/nproc)
        ! 'ngroup' has the number of lines in a group
        ! 'nlot'   has the maximum number of lines in a group
        allocate(ngroup(nproc))
        ngroup(:) = lot/nproc
        nremain = mod(lot,nproc)
        ngroup(1:nremain) = ngroup(1:nremain) + 1
        nlot = ngroup(1)

        ! Each process sends (10/equation)*(nlot equations)*(nproc in the comm)
        allocate(sendbuf(nlot*10*nproc))
        allocate(recvbuf(nlot,5,2*nproc))

        ! Allocation for rhs (recall there are multiple rhs per equations)
        allocate(sendbuf2(nlot*2*nrhs*nproc))
        allocate(recvbuf2(nrhs,nlot,2*nproc))


        ! Initialize boundary values
        s1(:,1) = a(:,1)
        s2(:,n) = c(:,n)
        if (irank == 0) s1(:,1) = 0.
        if (irank == nproc-1) s2(:,n) = 0.

        ! Forward elimination
        ! Upper boundary in s1(i)
        do i=2,n
                const(:) = a(:,i)/b(:,i-1)
                b(:,i)   = b(:,i) - c(:,i-1)*const(:)
                do j=1,lot
                   r(:,j,i) = r(:,j,i) - r(:,j,i-1)*const(j)
                end do
                s1(:,i)  = -s1(:,i-1)*const(:)
        end do

        ! Backward elimination
        ! Lower boundary in s2(i)
        do i=n-1,1,-1
                const(:) = c(:,i)/b(:,i+1)
                do j=1,lot
                   r(:,j,i) = r(:,j,i) - r(:,j,i+1)*const(j)
                end do
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

        ! (K1,K2) are the limits of the equations index
        ! 'K1' lower equation number in the lot
        ! 'K2' Upper equation number in the lot
        ! 'sendbuf' has (nlot S1, nlot b, ....) for process 1
        !           followed by (nlot S1, nlot S2, ...) for process 2

        L = 1
        k1 = 1
        do igroup=1,nproc
                k2 = k1+ngroup(igroup)-1
                nk = k2-k1+1

                sendbuf(L:L+nk-1) = 0.            ; L = L + nlot
                sendbuf(L:L+nk-1) = s1(k1:k2,1)   ; L = L + nlot
                sendbuf(L:L+nk-1) = b(k1:k2,1)    ; L = L + nlot
                sendbuf(L:L+nk-1) = 0.            ; L = L + nlot
                sendbuf(L:L+nk-1) = s2(k1:k2,1)   ; L = L + nlot
                sendbuf(L:L+nk-1) = s1(k1:k2,n)   ; L = L + nlot
                sendbuf(L:L+nk-1) = 0.            ; L = L + nlot
                sendbuf(L:L+nk-1) = b(k1:k2,n)    ; L = L + nlot
                sendbuf(L:L+nk-1) = s2(k1:k2,n)   ; L = L + nlot
                sendbuf(L:L+nk-1) = 0.            ; L = L + nlot

                k1 = k2 + 1
        end do

        ! Gather the boundary data
        call MPI_AllToAll (sendbuf, nlot*10, MPI_REAL8, &
                           recvbuf, nlot*10, MPI_REAL8, &
                           icomm, ierr)
        M = 1
        L = 1
        k1 = 1
        do igroup=1,nproc
                k2 = k1+ngroup(igroup)-1
                nk = k2-k1+1
                nkr = nk * nrhs

                do i=1,nk
                   sendbuf2(M:M+nrhs-1) = r(:,k1+i-1,1)
                   M = M+nrhs
                end do
                L = L + nlot*nrhs
                M = L
                do i=1,nk
                   sendbuf2(M:M+nrhs-1) = r(:,k1+i-1,n)
                   M = M+nrhs
                end do
                L = L + nlot*nrhs
                M = L

                k1 = k2 + 1
        end do

        ! Gather the boundary data
        call MPI_AllToAll (sendbuf2, nrhs*nlot*2, MPI_REAL8, &
                           recvbuf2, nrhs*nlot*2, MPI_REAL8, &
                           icomm, ierr)

        ! Clear unused values
        nk = ngroup(irank+1)
        recvbuf(nk+1:nlot, :,  :) = 0.
        recvbuf(nk+1:nlot, 3,  :) = 1.
        recvbuf2(:, nk+1:nlot, :) = 0.


        ! Solve reduced systems
        call pentaDiagonalM_1_Mrhs(recvbuf(1,1,2),recvbuf2(1,1,2),2*nproc-2,nlot,nrhs)
        recvbuf2(:,nk+1:nlot,:)=0.


        allocate(swap(nrhs,nlot))
        ! Permute the order
        do i=1,nproc-1
                swap(1:nrhs,1:nlot) = recvbuf2(:,:,2*i)
                recvbuf2(:,:,2*i   ) = recvbuf2(:,:,2*i+1)
                recvbuf2(:,:,2*i+1 ) = swap(1:nrhs,1:nlot)
        end do
        deallocate(swap)

        ! Scatter back the solution
        call MPI_AllToAll (recvbuf2, nrhs*nlot*2, MPI_REAL8, &
                           sendbuf2, nrhs*nlot*2, MPI_REAL8, &
                           icomm, ierr)

        M = 1
        L = 1
        k1 = 1
        do igroup=1,nproc
                k2 = k1+ngroup(igroup)-1
                nk = k2-k1+1
                nkr = nk*nrhs

                do i=1,nk
                   r1(:, k1+i-1) = sendbuf2(M:M+nrhs-1)
                   M = M+nrhs
                end do
                L = L + nlot*nrhs
                M = L

                do i=1,nk
                   r2(:, k1+i-1) = sendbuf2(M:M+nrhs-1)
                   M = M+nrhs
                end do
                L = L + nlot*nrhs
                M = L

                k1 = k2 + 1
        end do

        if (irank == 0) r1 = 0.
        if (irank == nproc-1) r2 = 0.

        do j=1,n
            do i=1,lot
                r(:,i,j) = (r(:,i,j) - s1(i,j)*r1(:,i) - s2(i,j)*r2(:,i))/b(i,j)
            end do
        end do

!        if (ibc == 1) then
!                if (irank .ne. 0) r(:,0) = r1(:)
!                if (irank .ne. nproc-1) r(:,n+1) = r2(:)
!        end if

        deallocate(sendbuf)
        deallocate(recvbuf)
        deallocate(sendbuf2)
        deallocate(recvbuf2)
        deallocate(ngroup)

        return
end


!=======================================================================


