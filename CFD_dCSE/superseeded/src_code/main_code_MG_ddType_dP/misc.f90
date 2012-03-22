! These are miscellaneous routines primarily for use with
! freestream.f90 and fluctuations.f90
! They are written to emulate MPI in a serial code

!=======================================================================
subroutine scatterY(A, B, ni)
	use messenger
	complex A(ny,nz,ni), B(ny_1,nz,ni)

	B(1:ny,1:nz,:) = A

	return
end

subroutine globalFSTSum(A, na)
	use messenger
	real A(na)

	return
end

subroutine computeOmega(omega,nw)
	use messenger
	real omega(nw/2)

	return
end

subroutine BCAST_TOP_EDGE(A,npnts,ixyz)
        use messenger
        integer :: npnts, ixyz, root_id
        real    :: A(npnts)

        root_id = 0

        return
end

