!=======================================================================
! Tri-diagonal solvers
!
! Charles Pierce, June 1995
!
! triDiagonal_mrhs
! triDiagonalC_mrhs
! triDiagonalM_mrhs
!

!=======================================================================
! triDiagonal_mrhs
! Solves a tri-diagonal system of algebraic equations using gauss 
! elimination.  The banded matrix given by (a,b,c) is destroyed, 
! and the rhs vector r is replaced by the solution vector.
!
! Single system, variable a,b,c.
!
subroutine triDiag_mrhs(a, b, c, r, nrhs, n)
	real a(n), b(n), c(n), r(nrhs,n)

	! Forward elimination
	do i=2,n
		const = a(i)/b(i-1)
		b(i) = b(i) - c(i-1)*const
		r(:,i) = r(:,i) - r(:,i-1)*const
	end do

	! Back substitution
	r(:,n) = r(:,n)/b(n)
	do i=n-1,1,-1
		r(:,i) = (r(:,i) - c(i)*r(:,i+1))/b(i)
	end do

	return
end

!=======================================================================
! triDiagonalM
! Solves an array of tri-diagonal systems of algebraic equations using 
! gauss elimination.  The banded matrices given by (a,b,c) are 
! destroyed, and the rhs vectors r are replaced by the solution vectors.
!
! Multiple systems, different (a,b,c) for each system.
!
subroutine triDiagM_mrhs(a, b, c, r, n, lot, nrhs)
	real a(lot,n), b(lot,n), c(lot,n), r(nrhs,lot,n)
	real const(lot)

	! Forward elimination
	do i=2,n
		const(:) = a(:,i)/b(:,i-1)
		b(:,i) = b(:,i) - c(:,i-1)*const(:)
            do j=1,lot
		r(:,j,i) = r(:,j,i) - r(:,j,i-1)*const(j)
            end do
	end do

	! Back-substitution
        do j=1,lot
	   r(:,j,n) = r(:,j,n)/b(j,n)
        end do

	do i=n-1,1,-1
           do j=1,lot
		r(:,j,i) = (r(:,j,i) - c(j,i)*r(:,j,i+1))/b(j,i)
           end do
	end do

	return
end


!=======================================================================
! triDiagonalC_mrhs
! Solves a tri-diagonal system of algebraic equations using gauss
! elimination.  The banded matrix given by (a,b,c) is destroyed,
! and the rhs vector r is replaced by the solution vector.
!
! Single system, variable a,b,c.
! All Variables are Complex
!
subroutine triDiagC_mrhs(a, b, c, r, nrhs, n)
	complex a(n), b(n), c(n), r(nrhs,n)
	complex const

	! Forward elimination
	do i=2,n
		const = a(i)/b(i-1)
		b(i) = b(i) - c(i-1)*const
		r(:,i) = r(:,i) - r(:,i-1)*const
	end do

	! Back substitution
	r(:,n) = r(:,n)/b(n)
	do i=n-1,1,-1
		r(:,i) = (r(:,i) - c(i)*r(:,i+1))/b(i)
	end do

return
end


