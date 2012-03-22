!=======================================================================
! Penta-diagonal solvers
!
! Charles Pierce, June 1995
!
! pentaDiagonalM
! pentaDiagonalLUM
! pentaLUsolveM
! pentaPeriodicM
! pentaDiagonal_1
! pentaDiagonalM_1
! pentaDiagonalLUM_1
! pentaLUsolveM_1
! pentaPeriodicM_1
! pentaDiagonalMC (complex variable version of pentaDiagonalM)
!

!=======================================================================
! pentaDiagonalM
! Solves a penta-diagonal system of algebraic equations.  The banded 
! matrix given by (a,b,c,d,e) is destroyed, and the rhs vector r is 
! replaced by the solution vector.
!
! Multiple systems, different (a,b,c,d,e) for each system.
!
subroutine pentaDiagonalM(a, b, c, d, e, r, n, lot)
	real a(lot,n), b(lot,n), c(lot,n), &
         d(lot,n), e(lot,n), r(lot,n)
	real const(lot)

	if (n == 1) then
		! Solve 1x1 system
		r(:,1) = r(:,1)/c(:,1)
		return
	else if (n == 2) then
		! Solve 2x2 system
		const(:) = b(:,2)/c(:,1)
		c(:,2) = c(:,2) - d(:,1)*const(:)
		r(:,2) = r(:,2) - r(:,1)*const(:)
		r(:,2) = r(:,2)/c(:,2)
		r(:,1) = (r(:,1) - d(:,1)*r(:,2))/c(:,1)
		return
	end if

	! Forward elimination
	do i=1,n-2
		! Eliminate b(i+1)
		const(:) = b(:,i+1)/c(:,i)
		c(:,i+1) = c(:,i+1) - d(:,i)*const(:)
		d(:,i+1) = d(:,i+1) - e(:,i)*const(:)
		r(:,i+1) = r(:,i+1) - r(:,i)*const(:)

		! Eliminate a(i+2)
		const(:) = a(:,i+2)/c(:,i)
		b(:,i+2) = b(:,i+2) - d(:,i)*const(:)
		c(:,i+2) = c(:,i+2) - e(:,i)*const(:)
		r(:,i+2) = r(:,i+2) - r(:,i)*const(:)
	end do
	! Eliminate b(n)
	const(:) = b(:,n)/c(:,n-1)
	c(:,n) = c(:,n) - d(:,n-1)*const(:)
	r(:,n) = r(:,n) - r(:,n-1)*const(:)

	! Back-substitution
	r(:,n) = r(:,n)/c(:,n)
	r(:,n-1) = (r(:,n-1) - d(:,n-1)*r(:,n))/c(:,n-1)
	do i=n-2,1,-1
		r(:,i) = (r(:,i) - d(:,i)*r(:,i+1) - e(:,i)*r(:,i+2))/c(:,i)
	end do

	return
end

!=======================================================================
! pentaDiagonalLUM
! Computes the LU decomposition if a penta-diagonal matrix using 
! Guass elimination.
!
! Multiple systems
!
subroutine pentaDiagonalLUM(a, b, c, d, e, n, lot)
	real a(lot,n), b(lot,n), c(lot,n), d(lot,n), e(lot,n)

	if (n == 1) then
		! 1x1 system
		c = 1./c
		return
	else if (n == 2) then
		! 2x2 system
		b(:,2) = b(:,2)/c(:,1)
		c(:,2) = c(:,2) - d(:,1)*b(:,2)
		c = 1./c
		return
	end if

	! Forward elimination
	do i=1,n-2
		! Eliminate b(i+1)
		b(:,i+1) = b(:,i+1)/c(:,i)
		c(:,i+1) = c(:,i+1) - d(:,i)*b(:,i+1)
		d(:,i+1) = d(:,i+1) - e(:,i)*b(:,i+1)

		! Eliminate a(i+2)
		a(:,i+2) = a(:,i+2)/c(:,i)
		b(:,i+2) = b(:,i+2) - d(:,i)*a(:,i+2)
		c(:,i+2) = c(:,i+2) - e(:,i)*a(:,i+2)
	end do
	! Eliminate b(n)
	b(:,n) = b(:,n)/c(:,n-1)
	c(:,n) = c(:,n) - d(:,n-1)*b(:,n)

	! Pre-divide diagonal
	c = 1./c

	return
end

!=======================================================================
! pentaLUSolveM
! Solves a penta-diagonal system, given the LU decomposition.
!
! Multiple systems
!
subroutine pentaLUSolveM(a, b, c, d, e, r, n, lot)
	real a(lot,n), b(lot,n), c(lot,n), &
         d(lot,n), e(lot,n), r(lot,n)

	if (n == 1) then
		! Solve 1x1 system
		r(:,1) = r(:,1)*c(:,1)
		return
	else if (n == 2) then
		! Solve 2x2 system
		r(:,2) = (r(:,2) - r(:,1)*b(:,2))*c(:,2)
		r(:,1) = (r(:,1) - d(:,1)*r(:,2))*c(:,1)
		return
	end if

	! Forward substitution
	do i=1,n-2
		r(:,i+1) = r(:,i+1) - r(:,i)*b(:,i+1)
		r(:,i+2) = r(:,i+2) - r(:,i)*a(:,i+2)
	end do
	r(:,n) = r(:,n) - r(:,n-1)*b(:,n)

	! Backward substitution
	! Diagonal has been pre-divided
	r(:,n) = r(:,n)*c(:,n)
	r(:,n-1) = (r(:,n-1) - d(:,n-1)*r(:,n))*c(:,n-1)
	do i=n-2,1,-1
		r(:,i) = (r(:,i) - d(:,i)*r(:,i+1) - e(:,i)*r(:,i+2))*c(:,i)
	end do

	return
end

!=======================================================================
! pentaPeriodicM
! Solves multiple periodic penta-diagonal systems
!
subroutine pentaPeriodicM(a, b, c, d, e, r, n, lot)
	real a(lot,n), b(lot,n), c(lot,n), d(lot,n), e(lot,n), r(lot,n)
	real s1(lot,n), s2(lot,n)
	real const(lot)

	if (n == 1) then
		! Solve 1x1 system
		r = r/(a + b + c + d + e)
		return
	else if (n == 2) then
		! Solve 2x2 system
		c = c + a + e
		d(:,1) = d(:,1) + b(:,1)
		b(:,2) = b(:,2) + d(:,2)
		const(:) = b(:,2)/c(:,1)
		c(:,2) = c(:,2) - d(:,1)*const(:)
		r(:,2) = r(:,2) - r(:,1)*const(:)
		r(:,2) = r(:,2)/c(:,2)
		r(:,1) = (r(:,1) - d(:,1)*r(:,2))/c(:,1)
		return
	else if (n == 3) then
		b = b + e
		d = d + a
		call triPeriodicM(b, c, d, r, n, lot)
		return
	else if (n == 4) then
		a = a + e
		e = 0.
	end if

	! Initialize boundary data
	s1 = 0.
	s1(:,1) = a(:,1)
	s1(:,n-3) = s1(:,n-3) + e(:,n-3)
	s1(:,n-2) = d(:,n-2)
	s1(:,n-1) = c(:,n-1)
	s1(:,n) = b(:,n)
	s2 = 0.
	s2(:,1) = b(:,1)
	s2(:,2) = a(:,2)
	s2(:,n-2) = s2(:,n-2) + e(:,n-2)
	s2(:,n-1) = d(:,n-1)
	s2(:,n) = c(:,n)

	! Forward elimination
	do i=1,n-2
		! Eliminate b(i+1)
		const(:) = b(:,i+1)/c(:,i)
		c(:,i+1) = c(:,i+1) - d(:,i)*const(:)
		d(:,i+1) = d(:,i+1) - e(:,i)*const(:)
		r(:,i+1) = r(:,i+1) - r(:,i)*const(:)
		s1(:,i+1) = s1(:,i+1) - s1(:,i)*const(:)
		s2(:,i+1) = s2(:,i+1) - s2(:,i)*const(:)

		! Eliminate a(i+2)
		const(:) = a(:,i+2)/c(:,i)
		b(:,i+2) = b(:,i+2) - d(:,i)*const(:)
		c(:,i+2) = c(:,i+2) - e(:,i)*const(:)
		r(:,i+2) = r(:,i+2) - r(:,i)*const(:)
		s1(:,i+2) = s1(:,i+2) - s1(:,i)*const(:)
		s2(:,i+2) = s2(:,i+2) - s2(:,i)*const(:)
	end do

	! Backward elimination
	do i=n-2,3,-1
		! Eliminate d(i-1)
		const(:) = d(:,i-1)/c(:,i)
		r(:,i-1) = r(:,i-1) - r(:,i)*const(:)
		s1(:,i-1) = s1(:,i-1) - s1(:,i)*const(:)
		s2(:,i-1) = s2(:,i-1) - s2(:,i)*const(:)

		! Eliminate e(i-2)
		const(:) = e(:,i-2)/c(:,i)
		r(:,i-2) = r(:,i-2) - r(:,i)*const(:)
		s1(:,i-2) = s1(:,i-2) - s1(:,i)*const(:)
		s2(:,i-2) = s2(:,i-2) - s2(:,i)*const(:)
	end do
	i=2
		! Eliminate d(i-1)
		const(:) = d(:,i-1)/c(:,i)
		r(:,i-1) = r(:,i-1) - r(:,i)*const(:)
		s1(:,i-1) = s1(:,i-1) - s1(:,i)*const(:)
		s2(:,i-1) = s2(:,i-1) - s2(:,i)*const(:)

	! Eliminate oddball region
	const(:) = e(:,n-1)/c(:,1)
	r(:,n-1) = r(:,n-1) - r(:,1)*const(:)
	s1(:,n-1) = s1(:,n-1) - s1(:,1)*const(:)
	s2(:,n-1) = s2(:,n-1) - s2(:,1)*const(:)

	const(:) = d(:,n)/c(:,1)
	r(:,n) = r(:,n) - r(:,1)*const(:)
	s1(:,n) = s1(:,n) - s1(:,1)*const(:)
	s2(:,n) = s2(:,n) - s2(:,1)*const(:)

	const(:) = e(:,n)/c(:,2)
	r(:,n) = r(:,n) - r(:,2)*const(:)
	s1(:,n) = s1(:,n) - s1(:,2)*const(:)
	s2(:,n) = s2(:,n) - s2(:,2)*const(:)

	! Eliminate corner region
	const(:) = s1(:,n)/s1(:,n-1)
	r(:,n) = r(:,n) - r(:,n-1)*const(:)
	s2(:,n) = s2(:,n) - s2(:,n-1)*const(:)

	r(:,n) = r(:,n)/s2(:,n)
	r(:,n-1) = (r(:,n-1) - s2(:,n-1)*r(:,n))/s1(:,n-1)
	do i=n-2,1,-1
		r(:,i) = (r(:,i) - s1(:,i)*r(:,n-1) - s2(:,i)*r(:,n))/c(:,i)
	end do

	return
end

!=======================================================================
! pentaDiagonal_1
!
! Combined matrix-rhs into single array
!
subroutine pentaDiagonal_1(A, n)
	real A(6,n)

	if (n == 1) then
		! Solve 1x1 system
		A(6,1) = A(6,1)/A(3,1)
		return
	else if (n == 2) then
		! Solve 2x2 system
		const = A(2,2)/A(3,1)
		A(3,2) = A(3,2) - A(4,1)*const
		A(6,2) = A(6,2) - A(6,1)*const
		A(6,2) = A(6,2)/A(3,2)
		A(6,1) = (A(6,1) - A(4,1)*A(6,2))/A(3,1)
		return
	end if

	! Forward elimination
	do i=1,n-2
		! Eliminate A(2,i+1)
		const = A(2,i+1)/A(3,i)
		A(3,i+1) = A(3,i+1) - A(4,i)*const
		A(4,i+1) = A(4,i+1) - A(5,i)*const
		A(6,i+1) = A(6,i+1) - A(6,i)*const

		! Eliminate A(1,i+2)
		const = A(1,i+2)/A(3,i)
		A(2,i+2) = A(2,i+2) - A(4,i)*const
		A(3,i+2) = A(3,i+2) - A(5,i)*const
		A(6,i+2) = A(6,i+2) - A(6,i)*const
	end do
	! Eliminate A(2,n)
	const = A(2,n)/A(3,n-1)
	A(3,n) = A(3,n) - A(4,n-1)*const
	A(6,n) = A(6,n) - A(6,n-1)*const

	! Back-substitution
	A(6,n) = A(6,n)/A(3,n)
	A(6,n-1) = (A(6,n-1) - A(4,n-1)*A(6,n))/A(3,n-1)
	do i=n-2,1,-1
		A(6,i) = (A(6,i) - A(4,i)*A(6,i+1) - A(5,i)*A(6,i+2))/A(3,i)
	end do

	return
end

!=======================================================================
! pentaDiagonalM_1
!
! Combined matrix-rhs into single array
!
subroutine pentaDiagonalM_1(A, n, lot)
	real A(lot,6,n)
	real const(lot)

	if (n == 1) then
		! Solve 1x1 system
		A(:,6,1) = A(:,6,1)/A(:,3,1)
		return
	else if (n == 2) then
		! Solve 2x2 system
		const(:) = A(:,2,2)/A(:,3,1)
		A(:,3,2) = A(:,3,2) - A(:,4,1)*const(:)
		A(:,6,2) = A(:,6,2) - A(:,6,1)*const(:)
		A(:,6,2) = A(:,6,2)/A(:,3,2)
		A(:,6,1) = (A(:,6,1) - A(:,4,1)*A(:,6,2))/A(:,3,1)
		return
	end if

	! Forward elimination
	do i=1,n-2
		! Eliminate A(2,i+1)
		const(:) = A(:,2,i+1)/A(:,3,i)
		A(:,3,i+1) = A(:,3,i+1) - A(:,4,i)*const(:)
		A(:,4,i+1) = A(:,4,i+1) - A(:,5,i)*const(:)
		A(:,6,i+1) = A(:,6,i+1) - A(:,6,i)*const(:)

		! Eliminate A(1,i+2)
		const(:) = A(:,1,i+2)/A(:,3,i)
		A(:,2,i+2) = A(:,2,i+2) - A(:,4,i)*const(:)
		A(:,3,i+2) = A(:,3,i+2) - A(:,5,i)*const(:)
		A(:,6,i+2) = A(:,6,i+2) - A(:,6,i)*const(:)
	end do
	! Eliminate A(2,n)
	const(:) = A(:,2,n)/A(:,3,n-1)
	A(:,3,n) = A(:,3,n) - A(:,4,n-1)*const(:)
	A(:,6,n) = A(:,6,n) - A(:,6,n-1)*const(:)

	! Back-substitution
	A(:,6,n) = A(:,6,n)/A(:,3,n)
	A(:,6,n-1) = (A(:,6,n-1) - A(:,4,n-1)*A(:,6,n))/A(:,3,n-1)
	do i=n-2,1,-1
		A(:,6,i) = (A(:,6,i) - A(:,4,i)*A(:,6,i+1) - A(:,5,i)*A(:,6,i+2))/A(:,3,i)
	end do

	return
end

!=======================================================================
! pentaDiagonalLUM_1
!
! Combined matrix-rhs into single array
!
subroutine pentaDiagonalLUM_1(A, n, lot)
	real A(lot,5,n)

	if (n == 1) then
		! 1x1 system
		A(:,3,:) = 1./A(:,3,:)
		return
	else if (n == 2) then
		! 2x2 system
		A(:,2,2) = A(:,2,2)/A(:,3,1)
		A(:,3,2) = A(:,3,2) - A(:,4,1)*A(:,2,2)
		A(:,3,:) = 1./A(:,3,:)
		return
	end if

	! Forward elimination
	do i=1,n-2
		! Eliminate b(i+1)
		A(:,2,i+1) = A(:,2,i+1)/A(:,3,i)
		A(:,3,i+1) = A(:,3,i+1) - A(:,4,i)*A(:,2,i+1)
		A(:,4,i+1) = A(:,4,i+1) - A(:,5,i)*A(:,2,i+1)

		! Eliminate a(i+2)
		A(:,1,i+2) = A(:,1,i+2)/A(:,3,i)
		A(:,2,i+2) = A(:,2,i+2) - A(:,4,i)*A(:,1,i+2)
		A(:,3,i+2) = A(:,3,i+2) - A(:,5,i)*A(:,1,i+2)
	end do
	! Eliminate b(n)
	A(:,2,n) = A(:,2,n)/A(:,3,n-1)
	A(:,3,n) = A(:,3,n) - A(:,4,n-1)*A(:,2,n)

	! Pre-divide diagonal
	A(:,3,:) = 1./A(:,3,:)

	return
end

!=======================================================================
! pentaLUSolveM_1
!
! Combined matrix-rhs into single array
!
subroutine pentaLUSolveM_1(A, r, n, lot)
	real A(lot,5,n), r(lot,n)

	if (n == 1) then
		! Solve 1x1 system
		r(:,1) = r(:,1)*A(:,3,1)
		return
	else if (n == 2) then
		! Solve 2x2 system
		r(:,2) = (r(:,2) - r(:,1)*A(:,2,2))*A(:,3,2)
		r(:,1) = (r(:,1) - A(:,4,1)*r(:,2))*A(:,3,1)
		return
	end if

	! Forward substitution
	do i=1,n-2
		r(:,i+1) = r(:,i+1) - r(:,i)*A(:,2,i+1)
		r(:,i+2) = r(:,i+2) - r(:,i)*A(:,1,i+2)
	end do
	r(:,n) = r(:,n) - r(:,n-1)*A(:,2,n)

	! Backward substitution
	! Diagonal has been pre-divided
	r(:,n) = r(:,n)*A(:,3,n)
	r(:,n-1) = (r(:,n-1) - A(:,4,n-1)*r(:,n))*A(:,3,n-1)
	do i=n-2,1,-1
		r(:,i) = (r(:,i) - A(:,4,i)*r(:,i+1) - A(:,5,i)*r(:,i+2))*A(:,3,i)
	end do

	return
end

!=======================================================================
! pentaPeriodicM_1
! Solves multiple periodic penta-diagonal systems
!
! Combined matrix-rhs into single array
!
subroutine pentaPeriodicM_1(A, n, lot)
	real A(lot,6,n)
	real s1(lot,n), s2(lot,n)
	real const(lot)

	if (n == 1) then
		! Solve 1x1 system
		A(:,6,:) = A(:,6,:)/( A(:,1,:) + A(:,2,:) + A(:,3,:) &
		                     +A(:,4,:) + A(:,5,:) )
		return
	else if (n == 2) then
		! Solve 2x2 system
		A(:,3,:) = A(:,3,:) + A(:,1,:) + A(:,5,:)
		A(:,4,1) = A(:,4,1) + A(:,2,1)
		A(:,2,2) = A(:,2,2) + A(:,4,2)
		const(:) = A(:,2,2)/A(:,3,1)
		A(:,3,2) = A(:,3,2) - A(:,4,1)*const(:)
		A(:,6,2) = A(:,6,2) - A(:,6,1)*const(:)
		A(:,6,2) = A(:,6,2)/A(:,3,2)
		A(:,6,1) = (A(:,6,1) - A(:,4,1)*A(:,6,2))/A(:,3,1)
		return
	else if (n == 3) then
		A(:,2,:) = A(:,2,:) + A(:,5,:)
		A(:,4,:) = A(:,4,:) + A(:,1,:)
		call triPeriodicM(A(:,2,:), A(:,3,:), A(:,4,:), A(:,6,:), n, lot)
		return
	else if (n == 4) then
		A(:,1,:) = A(:,1,:) + A(:,5,:)
		A(:,5,:) = 0.
	end if

	! Initialize boundary data
	s1 = 0.
	s1(:,1) = A(:,1,1)
	s1(:,n-3) = s1(:,n-3) + A(:,5,n-3)
	s1(:,n-2) = A(:,4,n-2)
	s1(:,n-1) = A(:,3,n-1)
	s1(:,n) = A(:,2,n)
	s2 = 0.
	s2(:,1) = A(:,2,1)
	s2(:,2) = A(:,1,2)
	s2(:,n-2) = s2(:,n-2) + A(:,5,n-2)
	s2(:,n-1) = A(:,4,n-1)
	s2(:,n) = A(:,3,n)

	! Forward elimination
	do i=1,n-2
		! Eliminate b(i+1)
		const(:) = A(:,2,i+1)/A(:,3,i)
		A(:,3,i+1) = A(:,3,i+1) - A(:,4,i)*const(:)
		A(:,4,i+1) = A(:,4,i+1) - A(:,5,i)*const(:)
		A(:,6,i+1) = A(:,6,i+1) - A(:,6,i)*const(:)
		s1(:,i+1) = s1(:,i+1) - s1(:,i)*const(:)
		s2(:,i+1) = s2(:,i+1) - s2(:,i)*const(:)

		! Eliminate a(i+2)
		const(:) = A(:,1,i+2)/A(:,3,i)
		A(:,2,i+2) = A(:,2,i+2) - A(:,4,i)*const(:)
		A(:,3,i+2) = A(:,3,i+2) - A(:,5,i)*const(:)
		A(:,6,i+2) = A(:,6,i+2) - A(:,6,i)*const(:)
		s1(:,i+2) = s1(:,i+2) - s1(:,i)*const(:)
		s2(:,i+2) = s2(:,i+2) - s2(:,i)*const(:)
	end do

	! Backward elimination
	do i=n-2,3,-1
		! Eliminate d(i-1)
		const(:) = A(:,4,i-1)/A(:,3,i)
		A(:,6,i-1) = A(:,6,i-1) - A(:,6,i)*const(:)
		s1(:,i-1) = s1(:,i-1) - s1(:,i)*const(:)
		s2(:,i-1) = s2(:,i-1) - s2(:,i)*const(:)

		! Eliminate e(i-2)
		const(:) = A(:,5,i-2)/A(:,3,i)
		A(:,6,i-2) = A(:,6,i-2) - A(:,6,i)*const(:)
		s1(:,i-2) = s1(:,i-2) - s1(:,i)*const(:)
		s2(:,i-2) = s2(:,i-2) - s2(:,i)*const(:)
	end do
	i=2
		! Eliminate d(i-1)
		const(:) = A(:,4,i-1)/A(:,3,i)
		A(:,6,i-1) = A(:,6,i-1) - A(:,6,i)*const(:)
		s1(:,i-1) = s1(:,i-1) - s1(:,i)*const(:)
		s2(:,i-1) = s2(:,i-1) - s2(:,i)*const(:)

	! Eliminate oddball region
	const(:) = A(:,5,n-1)/A(:,3,1)
	A(:,6,n-1) = A(:,6,n-1) - A(:,6,1)*const(:)
	s1(:,n-1) = s1(:,n-1) - s1(:,1)*const(:)
	s2(:,n-1) = s2(:,n-1) - s2(:,1)*const(:)

	const(:) = A(:,4,n)/A(:,3,1)
	A(:,6,n) = A(:,6,n) - A(:,6,1)*const(:)
	s1(:,n) = s1(:,n) - s1(:,1)*const(:)
	s2(:,n) = s2(:,n) - s2(:,1)*const(:)

	const(:) = A(:,5,n)/A(:,3,2)
	A(:,6,n) = A(:,6,n) - A(:,6,2)*const(:)
	s1(:,n) = s1(:,n) - s1(:,2)*const(:)
	s2(:,n) = s2(:,n) - s2(:,2)*const(:)

	! Eliminate corner region
	const(:) = s1(:,n)/s1(:,n-1)
	A(:,6,n) = A(:,6,n) - A(:,6,n-1)*const(:)
	s2(:,n) = s2(:,n) - s2(:,n-1)*const(:)

	A(:,6,n) = A(:,6,n)/s2(:,n)
	A(:,6,n-1) = (A(:,6,n-1) - s2(:,n-1)*A(:,6,n))/s1(:,n-1)
	do i=n-2,1,-1
		A(:,6,i) = (A(:,6,i) - s1(:,i)*A(:,6,n-1) - s2(:,i)*A(:,6,n))/A(:,3,i)
	end do

	return
end

!=======================================================================
! pentaDiagonalMC
! Solves a penta-diagonal system of algebraic equations.  The banded 
! matrix given by (a,b,c,d,e) is destroyed, and the rhs vector r is 
! replaced by the solution vector.
! All variables are complex.
!
! Multiple systems, different (a,b,c,d,e) for each system.
!
subroutine pentaDiagonalMC(a, b, c, d, e, r, n, lot)
	complex a(lot,n), b(lot,n), c(lot,n), &
         d(lot,n), e(lot,n), r(lot,n)
	complex const(lot)

	if (n == 1) then
		! Solve 1x1 system
		r(:,1) = r(:,1)/c(:,1)
		return
	else if (n == 2) then
		! Solve 2x2 system
		const(:) = b(:,2)/c(:,1)
		c(:,2) = c(:,2) - d(:,1)*const(:)
		r(:,2) = r(:,2) - r(:,1)*const(:)
		r(:,2) = r(:,2)/c(:,2)
		r(:,1) = (r(:,1) - d(:,1)*r(:,2))/c(:,1)
		return
	end if

	! Forward elimination
	do i=1,n-2
		! Eliminate b(i+1)
		const(:) = b(:,i+1)/c(:,i)
		c(:,i+1) = c(:,i+1) - d(:,i)*const(:)
		d(:,i+1) = d(:,i+1) - e(:,i)*const(:)
		r(:,i+1) = r(:,i+1) - r(:,i)*const(:)

		! Eliminate a(i+2)
		const(:) = a(:,i+2)/c(:,i)
		b(:,i+2) = b(:,i+2) - d(:,i)*const(:)
		c(:,i+2) = c(:,i+2) - e(:,i)*const(:)
		r(:,i+2) = r(:,i+2) - r(:,i)*const(:)
	end do
	! Eliminate b(n)
	const(:) = b(:,n)/c(:,n-1)
	c(:,n) = c(:,n) - d(:,n-1)*const(:)
	r(:,n) = r(:,n) - r(:,n-1)*const(:)

	! Back-substitution
	r(:,n) = r(:,n)/c(:,n)
	r(:,n-1) = (r(:,n-1) - d(:,n-1)*r(:,n))/c(:,n-1)
	do i=n-2,1,-1
		r(:,i) = (r(:,i) - d(:,i)*r(:,i+1) - e(:,i)*r(:,i+2))/c(:,i)
	end do

	return
end
