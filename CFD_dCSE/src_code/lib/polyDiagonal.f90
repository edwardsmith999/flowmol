!=======================================================================
! Poly-diagonal solvers
!
! Charles Pierce, June 1995
!
! polyDiagonal(nd, d, r, n)
! polyDiagonalM(nd, d, r, n, lot)
! polyDiagonalLU
! polyLUSolve
! polyPeriodicM
! polyDiagonal_1
! polyDiagonalM_1
! polyPeriodicM_1
!

!=======================================================================
! polyDiagonal
! Solves an N-diagonal banded system of algebraic equations using 
! Gauss elimation.
!
subroutine polyDiagonal(nd, d, r, n)
	real d(n,-nd:nd), r(n)

	! Forward elimination
	do i=1,n-1
		pivot = 1./d(i,0)
		do j=1,min(nd,n-i)
			const = d(i+j,-j)*pivot
			do k=1,min(nd,n-i)
				d(i+j,-j+k) = d(i+j,-j+k) - d(i,k)*const
			end do
			r(i+j) = r(i+j) - r(i)*const
		end do
	end do

	! Back-substitution
	do i=n,1,-1
		do j=1,min(nd,n-i)
			r(i) = r(i) - d(i,j)*r(i+j)
		end do
		r(i) = r(i)/d(i,0)
	end do

	return
end

!=======================================================================
! polyDiagonalM
! Solves multiple N-diagonal banded systems of algebraic equations 
! using Gauss elimation.
!
subroutine polyDiagonalM(nd, d, r, n, lot)
	real d(lot,n,-nd:nd), r(lot,n)
	real const(lot), pivot(lot)

	! Forward elimination
	do i=1,n-1
		pivot(:) = 1./d(:,i,0)
		do j=1,min(nd,n-i)
			const(:) = d(:,i+j,-j)*pivot(:)
			do k=1,min(nd,n-i)
				d(:,i+j,-j+k) = d(:,i+j,-j+k) - const(:)*d(:,i,k)
			end do
			r(:,i+j) = r(:,i+j) - const(:)*r(:,i)
		end do
	end do

	! Back-substitution
	do i=n,1,-1
		do j=1,min(nd,n-i)
			r(:,i) = r(:,i) - d(:,i,j)*r(:,i+j)
		end do
		r(:,i) = r(:,i)/d(:,i,0)
	end do

	return
end

!=======================================================================
! polyDiagonalLU
! Computes the LU decomposition of an N-diagonal banded matrix using 
! Gauss elimation.
!
subroutine polyDiagonalLU(nd, d, n)
	real d(-nd:nd,n)

	! Forward elimination
	do i=1,n-1
		pivot = 1./d(0,i)
		do j=1,min(nd,n-i)
			const = d(-j,i+j)*pivot
			d(-j,i+j) = const
			do k=1,min(nd,n-i)
				d(-j+k,i+j) = d(-j+k,i+j) - const*d(k,i)
			end do
		end do
	end do

	return
end

!=======================================================================
! polyLUSolve
! Solves an N-diagonal banded system given the LU decomposition.
!
subroutine polyLUSolve(nd, d, r, n)
	real d(-nd:nd,n), r(n)

	! Forward-substitution
	do i=1,n-1
		do j=1,min(nd,n-i)
			r(i+j) = r(i+j) - d(-j,i+j)*r(i)
		end do
	end do

	! Back-substitution
	do i=n,1,-1
		r(i) = r(i)/d(0,i)
		do j=1,min(nd,i-1)
			r(i-j) = r(i-j) - d(j,i-j)*r(i)
		end do
	end do

	return
end

!=======================================================================
! polyPeriodicM
! Solves multiple N-diagonal periodic banded systems 
! of algebraic equations using Gauss elimation.
!
subroutine polyPeriodicM(nd, d, r, n, lot)
	real d(lot,n,-nd:nd), r(lot,n)
	! Storage for the "bridge" or "boundary"
	real b(lot,n,nd)
	real const(lot), pivot(lot)

	if (n==1) then
		const(:) = 0.
		do j=-nd,nd
			const(:) = const(:) + d(:,1,j)
		end do
		r(:,1) = r(:,1) / const(:)
		return
	end if

	! Setup bridge array
	b = 0.
	do i=1,nd
		do j=i,nd
			b(:,i,j) = d(:,i,-nd+j-i)
			b(:,n-nd-i+1,j-i+1) = d(:,n-nd-i+1,j)
		end do
	end do
	do i=n-nd+1,n
		do j=1,nd
			b(:,i,j) = d(:,i,-nd+j-(i-n))
		end do
	end do

	! Forward elimination
	do i=1,n-nd
		pivot(:) = 1./d(:,i,0)
		do j=1,nd
			const(:) = d(:,i+j,-j)*pivot(:)
			do k=1,nd ! min(nd,n-nd-i)
				d(:,i+j,-j+k) = d(:,i+j,-j+k) - const(:)*d(:,i,k)
			! end do
			! do k=1,nd
				b(:,i+j,k) = b(:,i+j,k) - const(:)*b(:,i,k)
			end do
			r(:,i+j) = r(:,i+j) - const(:)*r(:,i)
		end do
	end do

	! Backward elimination
	do i=n-nd,1,-1
		pivot(:) = 1./d(:,i,0)
		do j=1,min(nd,i-1)
			const(:) = d(:,i-j,j)*pivot(:)
			do k=1,nd
				b(:,i-j,k) = b(:,i-j,k) - const(:)*b(:,i,k)
			end do
			r(:,i-j) = r(:,i-j) - const(:)*r(:,i)
		end do
	end do

	! Eliminate oddball region
	do i=1,nd
		pivot(:) = 1./d(:,i,0)
		do j=i,nd
			const(:) = d(:,n-j+i,j)*pivot(:)
			do k=1,nd
				b(:,n-j+i,k) = b(:,n-j+i,k) - const(:)*b(:,i,k)
			end do
			r(:,n-j+i) = r(:,n-j+i) - const(:)*r(:,i)
		end do
	end do

	! Elimination for corner matrix
	i0 = n-nd
	do i=1,nd-1
		pivot(:) = 1./b(:,i+i0,i)
		do j=i+1,nd
			const(:) = b(:,j+i0,i)*pivot(:)
			do k=i+1,nd
				b(:,j+i0,k) = b(:,j+i0,k) - const(:)*b(:,i+i0,k)
			end do
			r(:,j+i0) = r(:,j+i0) - const(:)*r(:,i+i0)
		end do
	end do

	! Back-substitution for corner matrix
	i0 = n-nd
	do i=nd,1,-1
		do j=i+1,nd
			r(:,i+i0) = r(:,i+i0) - b(:,i+i0,j)*r(:,j+i0)
		end do
		r(:,i+i0) = r(:,i+i0)/b(:,i+i0,i)
	end do

	! Back-substitution for bridge
	do i=n-nd,1,-1
		do j=1,nd
			r(:,i) = r(:,i) - b(:,i,j)*r(:,n-nd+j)
		end do
		r(:,i) = r(:,i)/d(:,i,0)
	end do

	return
end

!=======================================================================
subroutine polyDiagonal_1(nd, A, n)
	real A(-nd:nd + 1,n)

	! Forward elimination
	do i=1,n-1
		pivot = 1./A(0,i)
		do j=1,min(nd,n-i)
			const = A(-j,i+j)*pivot
			do k=1,min(nd,n-i)
				A(-j+k,i+j) = A(-j+k,i+j) - A(k,i)*const
			end do
			A(nd+1,i+j) = A(nd+1,i+j) - A(nd+1,i)*const
		end do
	end do

	! Back-substitution
	do i=n,1,-1
		do j=1,min(nd,n-i)
			A(nd+1,i) = A(nd+1,i) - A(j,i)*A(nd+1,i+j)
		end do
		A(nd+1,i) = A(nd+1,i)/A(0,i)
	end do

	return
end

!=======================================================================
subroutine polyDiagonalM_1(nd, A, n, lot)
	real A(lot, -nd:nd + 1, n)
	real const(lot), pivot(lot)

	! Forward elimination
	do i=1,n-1
		pivot(:) = 1./A(:,0,i)
		do j=1,min(nd,n-i)
			const(:) = A(:,-j,i+j)*pivot(:)
			do k=1,min(nd,n-i)
				A(:,-j+k,i+j) = A(:,-j+k,i+j) - A(:,k,i)*const(:)
			end do
			A(:,nd+1,i+j) = A(:,nd+1,i+j) - A(:,nd+1,i)*const(:)
		end do
	end do

	! Back-substitution
	do i=n,1,-1
		do j=1,min(nd,n-i)
			A(:,nd+1,i) = A(:,nd+1,i) - A(:,j,i)*A(:,nd+1,i+j)
		end do
		A(:,nd+1,i) = A(:,nd+1,i)/A(:,0,i)
	end do

	return
end

!=======================================================================
subroutine polyPeriodicM_1(nd, A, n, lot)
	real d(lot,-nd:nd,n), r(lot,n)
	real A(lot, -nd:nd + 1, n)
	! Storage for the "bridge" or "boundary"
	real b(lot,nd,n)
	real const(lot), pivot(lot)

	if (n==1) then
		const(:) = 0.
		do j=-nd,nd
			const(:) = const(:) + A(:,1,j)
		end do
		A(:,nd+1,1) = A(:,nd+1,1) / const(:)
		return
	end if

	! Setup bridge array
	b = 0.
	do i=1,nd
		do j=i,nd
			b(:,j,i) = A(:,-nd+j-i,i)
			b(:,j-i+1,n-nd-i+1) = A(:,j,n-nd-i+1)
		end do
	end do
	do i=n-nd+1,n
		do j=1,nd
			b(:,j,i) = A(:,-nd+j-(i-n),i)
		end do
	end do

	! Forward elimination
	do i=1,n-nd
		pivot(:) = 1./A(:,0,i)
		do j=1,nd
			const(:) = A(:,-j,i+j)*pivot(:)
			do k=1,nd ! min(nd,n-nd-i)
				A(:,-j+k,i+j) = A(:,-j+k,i+j) - A(:,k,i)*const(:)
			! end do
			! do k=1,nd
				b(:,k,i+j) = b(:,k,i+j) - b(:,k,i)*const(:)
			end do
			A(:,nd+1,i+j) = A(:,nd+1,i+j) - A(:,nd+1,i)*const(:)
		end do
	end do

	! Backward elimination
	do i=n-nd,1,-1
		pivot(:) = 1./A(:,0,i)
		do j=1,min(nd,i-1)
			const(:) = A(:,j,i-j)*pivot(:)
			do k=1,nd
				b(:,k,i-j) = b(:,k,i-j) - b(:,k,i)*const(:)
			end do
			A(:,nd+1,i-j) = A(:,nd+1,i-j) - A(:,nd+1,i)*const(:)
		end do
	end do

	! Eliminate oddball region
	do i=1,nd
		pivot(:) = 1./A(:,0,i)
		do j=i,nd
			const(:) = A(:,j,n-j+i)*pivot(:)
			do k=1,nd
				b(:,k,n-j+i) = b(:,k,n-j+i) - b(:,k,i)*const(:)
			end do
			A(:,nd+1,n-j+i) = A(:,nd+1,n-j+i) - A(:,nd+1,i)*const(:)
		end do
	end do

	! Elimination for corner matrix
	i0 = n-nd
	do i=1,nd-1
		pivot(:) = 1./b(:,i,i+i0)
		do j=i+1,nd
			const(:) = b(:,i,j+i0)*pivot(:)
			do k=i+1,nd
				b(:,k,j+i0) = b(:,k,j+i0) - b(:,k,i+i0)*const(:)
			end do
			A(:,nd+1,j+i0) = A(:,nd+1,j+i0) - A(:,nd+1,i+i0)*const(:)
		end do
	end do

	! Back-substitution for corner matrix
	i0 = n-nd
	do i=nd,1,-1
		do j=i+1,nd
			A(:,nd+1,i+i0) = A(:,nd+1,i+i0) - b(:,j,i+i0)*A(:,nd+1,j+i0)
		end do
		A(:,nd+1,i+i0) = A(:,nd+1,i+i0)/b(:,i,i+i0)
	end do

	! Back-substitution for bridge
	do i=n-nd,1,-1
		do j=1,nd
			A(:,nd+1,i) = A(:,nd+1,i) - b(:,j,i)*A(:,nd+1,n-nd+j)
		end do
		A(:,nd+1,i) = A(:,nd+1,i)/A(:,0,i)
	end do

	return
end

