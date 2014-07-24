!=======================================================================
! Message-passing routines for parallel processing using MPI
!
! messenger_invoke()
! messenger_init()
! messenger_free()
!
! globalSum(A, na)
! globalMax(A, na)
! globalMin(A, na)
! globalDirSum(A, na, ixyz)
!
! scatter(A, B)
! gather(A, B)
! allGather(A, B)
! scatterXY(A, B, nk)
! gatherXY(A, B, nk)
! allGatherXY(A, B, nk)
! stat_gatherXY(A, B, nk)
! stat_allGatherXY(A, B, nk)
! transpose(R1, R2, isign)
! transposeA(R1, R2, isign)
! transpose_Mine(R1, R2, isign)
! transpose_qr_qT(R1,R2)
! transpose_phat_p(R1,R2)
! updateBorder(A, n1,n2,n3, ixyz, ijk)
!
! triDiagonal_ (ixyz, a, b, c, r, n)
! triDiagonalM_(ixyz, a, b, c, r, n, lot)
! triDiagonalLUM_(ixyz, a, b, c, work, n, lot)
! triLUSolveM_(ixyz, a, b, c, r, work, n, lot)
!

module messenger
	use continuum_data_export
#if USE_COUPLER
    use coupler
#endif
    save
       
	integer myid                      ! my process rank
	integer idroot                    ! rank of root process
    integer CFD_COMM                  ! CFD communicator
    integer ierr                  	  ! error flag

	! Grid topology
	integer icomm_grid                ! comm for grid topology
	integer icomm_xyz(3)              ! directional subcomms

	real(kind(0.d0)) wallTime

end module

!=======================================================================
subroutine messenger_invoke()
    use mpi
#if USE_COUPLER
    use coupler
#endif
	use messenger
!    use continuum_coupler_socket_init

    call MPI_init (ierr)

#if USE_COUPLER
            !call init_coupler(CFD_COMM)
            call coupler_create_comm(COUPLER_CFD, CFD_COMM, ierror)
            prefix_dir ="./couette_data/"
#else
            CFD_COMM = MPI_COMM_WORLD
            prefix_dir = "./"
#endif

    wallTime = mpi_wtime()

	! If coupling is used MPI initialisation is done at the top level
	! of coupled program

	return

end subroutine messenger_invoke


subroutine messenger_init()
    use mpi
    use continuum_data_export, only : icoord
	use messenger

	integer idims(3)
	logical Lperiodic(3)
	logical Lremain_dims(3)
        integer np, ndims, ip, ixyz

	! Initialize MPI
	call MPI_comm_size (CFD_COMM, np, ierr)
	call MPI_comm_rank (CFD_COMM, myid, ierr)
    if (np .ne. nproc) then 
            write(0,'(a,I0,a)') "rank ", myid, "Wrong number of processors in CFD_COMM"
            call mpi_abort(CFD_COMM,1,ierr)
    endif

	! Grid topology
	ndims = 3
	idims(1) = npx
	idims(2) = npy
	idims(3) = 1
	Lperiodic = .true.
	call MPI_Cart_create(CFD_COMM, ndims, idims, Lperiodic, .true., &
	                     icomm_grid, ierr)
	do ip=1,nproc
		call MPI_Cart_coords(icomm_grid, ip-1, ndims, icoord(1,ip), ierr)
	end do
	icoord = icoord + 1
	call MPI_comm_rank (icomm_grid, irank, ierr)
	irank = irank + 1
	iblock = icoord(1, irank)
	jblock = icoord(2, irank)
	kblock = icoord(3, irank)

	! Directional subcomms
	do ixyz=1,3
		Lremain_dims(:) = .false.
		Lremain_dims(ixyz) = .true.
		call MPI_Cart_sub (icomm_grid, Lremain_dims, icomm_xyz(ixyz), ierr)
	enddo
	call MPI_comm_rank (icomm_xyz(1), irankx, ierr)
	call MPI_comm_rank (icomm_xyz(2), iranky, ierr)
	call MPI_comm_rank (icomm_xyz(3), irankz, ierr)
        ! print *, " Old comm=",irank-1,"new x,y,z:",irankx,iranky

	! Root process at coordinates (0,0,0)
	idims = 0
	call MPI_Cart_rank(icomm_grid, idims, idroot, ierr)
	iroot = idroot + 1

end

subroutine messenger_free()
	use messenger
	use mpi

	! Report time used
	print "(a,f8.2)", "time: ", MPI_wtime() - wallTime

	! Finalize MPI
	call MPI_finalize (ierr)

	return
end

!=======================================================================
subroutine globalSum(A, na)
	use messenger
	use mpi

        integer na
	real*8 A(na)
	real*8 buf(na)

	call MPI_AllReduce (A, buf, na, MPI_REAL8, &
	                    MPI_SUM, CFD_COMM, ierr)
	A = buf

	return
end

subroutine globalMax(A, na)
	use messenger
	use mpi

        integer na
	real*8 A(na)
	real*8 buf(na)

	call MPI_AllReduce (A, buf, na, MPI_REAL8, &
	                    MPI_MAX, CFD_COMM, ierr)
	A = buf

	return
end

subroutine globalMin(A, na)
	use messenger
	use mpi

        integer na
	real*8 A(na)
	real*8 buf(na)

	call MPI_AllReduce (A, buf, na, MPI_REAL8, &
	                    MPI_MIN, CFD_COMM, ierr)
	A = buf

	return
end

subroutine globalDirSum(A, na, ixyz)
	use messenger
	use mpi

        integer na
	real*8 A(na)
	real*8 buf(na)

	call MPI_AllReduce (A, buf, na, MPI_REAL8, &
	                    MPI_SUM, icomm_xyz(ixyz), ierr)
	A = buf

	return
end

