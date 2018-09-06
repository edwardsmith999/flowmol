module lapack_fns 

contains 

subroutine pinverse(A, PINV)
    implicit none

    double complex, allocatable, dimension(:,:), intent(in) :: A
    double complex, allocatable, dimension(:,:), intent(out) :: PINV

    integer :: i, j, M, N, K, L, LWORK, INFO

    double precision, allocatable, dimension(:) :: S
    double complex, allocatable, dimension(:) :: WORK, RWORK
    double complex, allocatable, dimension(:,:) :: U, VT, BUFF

    M = size(A,1)
    N = size(A,2)
    K = min(M,N)
    L = max(M,N)
    LWORK = MAX(1,2*K+L)

    allocate(PINV(N,M), U(M,K), BUFF(N,N))
    allocate(VT(K,N), WORK(LWORK), RWORK(5*K), S(K))

    !Compute the SVD of A1
    call ZGESVD( 'S', 'S', M, N, A, M, S, U, & 
                M, VT, K, WORK, LWORK, RWORK, INFO)

    !Compute PINV = VT**T * SIGMA * U**T in two steps
    do j = 1, K
        call ZSCAL( M, dcmplx(1 / S( j )), U( 1, j ), 1 )
    end do
    call ZGEMM( 'C', 'C', N, M, K, dcmplx(1.0), & 
                VT, K, U, M, dcmplx(0.0), PINV, N)

end subroutine pinverse 

end module lapack_fns

!program PseudoInverse
!    use lapack_fns
!    Implicit none
!    
!    external ZLANGE
!    double precision ZLANGE
!    
!    integer i, j, M, N, K, L, LWORK, INFO
!    parameter (M=15)
!    parameter (N=10)

!    parameter (K = MIN(M,N))
!    parameter (L = MAX(M,N))

!    parameter (LWORK = MAX(1,2*K+L))

!    double complex, allocatable, dimension(:,:) :: A3, pinvr   
!    double complex, dimension(M,N) :: A1, A2, SIGMA
!    double complex, dimension(N,M) :: PINV
!    double complex, dimension(M,K) :: U
!    double complex, dimension(K,N) :: VT
!    double complex, dimension(N,N) :: BUFF
!    double complex, dimension(LWORK) :: WORK
!    double precision, dimension(5*K) :: RWORK
!    double precision, dimension(K) :: S
!    integer, dimension(4) :: ISEED

!    double precision :: normA, normAPA, normPAP

!    data ISEED/0,0,0,1/

!    !Fill A1 with random values and copy into A2
!    allocate(A3(M,N))
!    call ZLARNV( 1, ISEED, M*N, A1 )
!    do i=1,M
!        do j=1,N
!            A3(i,j) = A1(i,j)
!            print*, A3(i,j)
!        end do
!    end do

!    call pinverse(A3, pinvr)

!    do i=1,M
!        do j=1,N
!            print*, pinvr(i,j)
!        end do
!    end do

!    ! check the result
!    normA = ZLANGE( 'F', M, N, A1, M)
!    call ZGEMM( 'N', 'N', N, N, M, dcmplx(1.0), & 
!                 PINVr, N, A1, M, dcmplx(0.0), BUFF, N )
!    call ZGEMM( 'N', 'N', M, N, N, dcmplx(-1.0), & 
!                 A1, M, BUFF, N, dcmplx(1.0), A1, M )
!    normAPA = ZLANGE( 'F', M, N, A1, M )

!    call ZGEMM( 'N', 'N', N, M, N, dcmplx(-1.0), &
!                  BUFF, N, PINVr, N, dcmplx(1.0), PINVr, N )
!    normPAP = ZLANGE( 'F', N, M, PINVr, N )

!    write(*,"(A, e10.4)") '|| A - A*P*A || = ', normAPA/normA
!    write(*,"(A, e10.4)") '|| P - P*A*P || = ', normPAP/normA

!end program PseudoInverse

