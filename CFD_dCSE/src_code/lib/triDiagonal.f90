!=======================================================================
! Tri-diagonal solvers
!
! Charles Pierce, June 1995
!
! triDiagonal
! triDiagonalM
! triDiagonalMC
! triDiagonalLUM
! triLUSolveM
! triPeriodic
! triPeriodicM
! triDiagonal_1
! triDiagonalM_1
! triDiagonalLUM_1
! triLUSolveM_1
! triDiagMultM
! triPeriMultM
!

!=======================================================================
! triDiagonal
! Solves a tri-diagonal system of algebraic equations using gauss 
! elimination.  The banded matrix given by (a,b,c) is destroyed, 
! and the rhs vector r is replaced by the solution vector.
!
! Single system, variable a,b,c.
!
subroutine triDiagonal(a, b, c, r, n)
	real a(n), b(n), c(n), r(n)
        real const
        integer n, i

	! Forward elimination
	do i=2,n
		const = a(i)/b(i-1)
		b(i) = b(i) - c(i-1)*const
		r(i) = r(i) - r(i-1)*const
	end do

	! Back substitution
	r(n) = r(n)/b(n)
	do i=n-1,1,-1
		r(i) = (r(i) - c(i)*r(i+1))/b(i)
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
subroutine triDiagonalM(a, b, c, r, n, lot)
	real a(lot,n), b(lot,n), c(lot,n), r(lot,n)
	real const(lot)

	! Forward elimination
	do i=2,n
		const(:) = a(:,i)/b(:,i-1)
		b(:,i) = b(:,i) - c(:,i-1)*const(:)
		r(:,i) = r(:,i) - r(:,i-1)*const(:)
	end do

	! Back-substitution
	r(:,n) = r(:,n)/b(:,n)
	do i=n-1,1,-1
		r(:,i) = (r(:,i) - c(:,i)*r(:,i+1))/b(:,i)
	end do

	return
end

!=======================================================================
! triDiagonalMC
! Solves a tri-diagonal system of algebraic equations.  The banded
! matrix given by (a,b,c) is destroyed, and the rhs vector r is
! replaced by the solution vector.
! All variables are complex.
!
! Multiple systems, different (a,b,c) for each system.
!
subroutine triDiagonalMC(a, b, c, r, n, lot)
	complex a(lot,n), b(lot,n), c(lot,n), r(lot,n)
	complex const(lot)

	! Forward elimination
	do i=2,n
		const(:) = a(:,i)/b(:,i-1)
		b(:,i) = b(:,i) - c(:,i-1)*const(:)
		r(:,i) = r(:,i) - r(:,i-1)*const(:)
	end do

	! Back-substitution
	r(:,n) = r(:,n)/b(:,n)
	do i=n-1,1,-1
		r(:,i) = (r(:,i) - c(:,i)*r(:,i+1))/b(:,i)
	end do

	return
end

!=======================================================================
! triDiagonalLUM
! Computes the LU decomposition if a tri-diagonal matrix using 
! Guass elimination.
!
! Multiple systems
!
subroutine triDiagonalLUM(a, b, c, n, lot)
	real a(lot,n), b(lot,n), c(lot,n)
	real const(lot)

	! Forward elimination
	do i=2,n
		const(:) = a(:,i)/b(:,i-1)
		a(:,i) = const(:)
		b(:,i) = b(:,i) - c(:,i-1)*const(:)
	end do

	! Pre-divide diagonal
	do i=1,n
		b(:,i) = 1./b(:,i)
	end do

	return
end

!=======================================================================
! triLUSolveM
! Solves a tri-diagonal system, given the LU decomposition.
!
! Multiple systems
!
subroutine triLUSolveM(a, b, c, r, n, lot)
	real a(lot,n), b(lot,n), c(lot,n), r(lot,n)

	! Forward-substitution
	do i=2,n
		r(:,i) = r(:,i) - a(:,i)*r(:,i-1)
	end do

	! Back-substitution
	r(:,n) = r(:,n)*b(:,n)
	do i=n-1,1,-1
		r(:,i) = (r(:,i) - c(:,i)*r(:,i+1))*b(:,i)
	end do

	return
end

!=======================================================================
! triPeriodic
! Solves a periodic, tri-diagonal system of algebraic equations 
! using gauss elimination.  The periodic banded matrix given by 
! (a,b,c) is destroyed, and the rhs vector r is replaced by the 
! solution vector.
!
! Single system, periodic diagonals, variable a,b,c.
!
subroutine triPeriodic(a, b, c, r, n)
	real a(n), b(n), c(n), r(n)

	if (n == 1) then
		r = r/(a + b + c)
		return
	else if (n == 2) then
		! Solve 2x2 system
		c(1) = c(1) + a(1)
		a(2) = a(2) + c(2)
		const = a(2)/b(1)
		b(2) = b(2) - c(1)*const
		r(2) = r(2) - r(1)*const
		r(2) = r(2)/b(2)
		r(1) = (r(1) - c(1)*r(2))/b(1)
		return
	end if

	! Forward elimination
	do i=2,n-2
		const = a(i)/b(i-1)
		b(i) = b(i) - c(i-1)*const
		r(i) = r(i) - r(i-1)*const
		! Boundary is stored in a(i)
		a(i) = -a(i-1)*const
	end do
	i=n-1
		const = a(i)/b(i-1)
		b(i) = b(i) - c(i-1)*const
		r(i) = r(i) - r(i-1)*const
		a(i) = c(i) - a(i-1)*const
	i=n
		const = a(i)/b(i-1)
		r(i) = r(i) - r(i-1)*const
		a(i) = b(i) - a(i-1)*const

	! Backward elimination
	do i=n-2,1,-1
		const = c(i)/b(i+1)
		r(i) = r(i) - r(i+1)*const
		a(i) = a(i) - a(i+1)*const
	end do

	! Eliminate oddball
	const = c(n)/b(1)
	r(n) = r(n) - r(1)*const
	a(n) = a(n) - a(1)*const

	! Backward substitution
	r(n) = r(n)/a(n)
	do i=n-1,1,-1
		r(i) = (r(i) - a(i)*r(n))/b(i)
	end do

	return
end

!=======================================================================
! triPeriodicM
! Solves an array of periodic, tri-diagonal systems of algebraic equations 
! using gauss elimination.  The periodic banded matrix given by 
! (a,b,c) is destroyed, and the rhs vector r is replaced by the 
! solution vector.
!
! Multiple systems, periodic diagonals, variable a,b,c.
! Different a,b,c for each system
!
subroutine triPeriodicM(a, b, c, r, n, lot)
	real a(lot,n), b(lot,n), c(lot,n), r(lot,n)
	real const(lot)

	if (n == 1) then
		r = r/(a + b + c)
		return
	else if (n == 2) then
		! Solve 2x2 system
		c(:,1) = c(:,1) + a(:,1)
		a(:,2) = a(:,2) + c(:,2)
		const(:) = a(:,2)/b(:,1)
		b(:,2) = b(:,2) - c(:,1)*const(:)
		r(:,2) = r(:,2) - r(:,1)*const(:)
		r(:,2) = r(:,2)/b(:,2)
		r(:,1) = (r(:,1) - c(:,1)*r(:,2))/b(:,1)
		return
	end if

	! Forward elimination
	do i=2,n-2
		const(:) = a(:,i)/b(:,i-1)
		b(:,i) = b(:,i) - c(:,i-1)*const(:)
		r(:,i) = r(:,i) - r(:,i-1)*const(:)
		! Boundary is stored in a(i)
		a(:,i) = -a(:,i-1)*const(:)
	end do
	i=n-1
		const(:) = a(:,i)/b(:,i-1)
		b(:,i) = b(:,i) - c(:,i-1)*const(:)
		r(:,i) = r(:,i) - r(:,i-1)*const(:)
		a(:,i) = c(:,i) - a(:,i-1)*const(:)
	i=n
		const(:) = a(:,i)/b(:,i-1)
		r(:,i) = r(:,i) - r(:,i-1)*const(:)
		a(:,i) = b(:,i) - a(:,i-1)*const(:)

	! Backward elimination
	do i=n-2,1,-1
		const(:) = c(:,i)/b(:,i+1)
		r(:,i) = r(:,i) - r(:,i+1)*const(:)
		a(:,i) = a(:,i) - a(:,i+1)*const(:)
	end do

	! Eliminate oddball
	const(:) = c(:,n)/b(:,1)
	r(:,n) = r(:,n) - r(:,1)*const(:)
	a(:,n) = a(:,n) - a(:,1)*const(:)

	! Backward substitution
	r(:,n) = r(:,n)/a(:,n)
	do i=n-1,1,-1
		r(:,i) = (r(:,i) - a(:,i)*r(:,n))/b(:,i)
	end do

	return
end

!=======================================================================
! triDiagonal_1
!
! Combined matrix-rhs into single array
!
subroutine triDiagonal_1(A, n)
	real A(4,n)

	if (n == 1) then
		! Solve 1x1 system
		A(4,1) = A(4,1)/A(2,1)
		return
	end if

	! Forward elimination
	do i=2,n
		const = A(1,i)/A(2,i-1)
		A(2,i) = A(2,i) - A(3,i-1)*const
		A(4,i) = A(4,i) - A(4,i-1)*const
	end do

	! Back-substitution
	A(4,n) = A(4,n)/A(2,n)
	do i=n-1,1,-1
		A(4,i) = (A(4,i) - A(3,i)*A(4,i+1))/A(2,i)
	end do

	return
end

subroutine triDiagonalM_1(A, n, lot)
	real A(lot,4,n)
	real const(lot)

	if (n == 1) then
		! Solve 1x1 system
		A(:,4,1) = A(:,4,1)/A(:,2,1)
		return
	end if

	! Forward elimination
	do i=2,n
		const(:) = A(:,1,i)/A(:,2,i-1)
		A(:,2,i) = A(:,2,i) - A(:,3,i-1)*const(:)
		A(:,4,i) = A(:,4,i) - A(:,4,i-1)*const(:)
	end do

	! Back-substitution
	A(:,4,n) = A(:,4,n)/A(:,2,n)
	do i=n-1,1,-1
		A(:,4,i) = (A(:,4,i) - A(:,3,i)*A(:,4,i+1))/A(:,2,i)
	end do

	return
end

subroutine triDiagonalLUM_1(A, n, lot)
	real A(lot,3,n)
	real const(lot)

	! Forward elimination
	do i=2,n
		const(:) = A(:,1,i)/A(:,2,i-1)
		A(:,1,i) = const(:)
		A(:,2,i) = A(:,2,i) - A(:,3,i-1)*const(:)
	end do

	! Pre-divide diagonal
	do i=1,n
		A(:,2,i) = 1./A(:,2,i)
	end do

	return
end

subroutine triLUSolveM_1(A, r, n, lot)
	real A(lot,3,n), r(lot,n)

	! Forward-substitution
	do i=2,n
		r(:,i) = r(:,i) - A(:,1,i)*r(:,i-1)
	end do

	! Back-substitution
	r(:,n) = r(:,n)*A(:,2,n)
	do i=n-1,1,-1
		r(:,i) = (r(:,i) - A(:,3,i)*r(:,i+1))*A(:,2,i)
	end do

	return
end

!=======================================================================
! triDiagMultM
! Multiplies a vector by a tri-diagonal matrix.
!
subroutine triDiagMultM(a, b, c, r, n, lot)
	real a(lot,n), b(lot,n), c(lot,n), r(lot,n)
	real rr(lot,n)

	if (n == 1) then
		rr = b*r
	else
		rr(:,1) = b(:,1)*r(:,1) + c(:,1)*r(:,2)
		do i=2,n-1
			rr(:,i) = a(:,i)*r(:,i-1) + b(:,i)*r(:,i) + c(:,i)*r(:,i+1)
		end do
		rr(:,n) = a(:,n)*r(:,n-1) + b(:,n)*r(:,n)
	end if

	r = rr

	return
end

!=======================================================================
! triPeriMultM
! Multiplies a vector by periodic tri-diagonal matrix.
!
subroutine triPeriMultM(a, b, c, r, n, lot)
	real a(lot,n), b(lot,n), c(lot,n), r(lot,n)
	real rr(lot,n)

	if (n==1) then
		rr = r * (a + b + c)
	else
		rr(:,1) = a(:,1)*r(:,n) + b(:,1)*r(:,1) + c(:,1)*r(:,2)
		do i=2,n-1
			rr(:,i) = a(:,i)*r(:,i-1) + b(:,i)*r(:,i) + c(:,i)*r(:,i+1)
		end do
		rr(:,n) = a(:,n)*r(:,n-1) + b(:,n)*r(:,n) + c(:,n)*r(:,1)
	end if

	r = rr

	return
end

