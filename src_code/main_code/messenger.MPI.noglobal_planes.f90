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
! triPeriodic_(ixyz, a, b, c, r, nl)
! triPeriodicM_(ixyz, a, b, c, r, nl, lot)
!

module messenger
	use data_export

    integer CFD_COMM                  ! CFD communicator, replaces MPI_COMM_WORLD
	integer myid                      ! my process rank
	integer idroot                    ! rank of root process

	! Grid topology
	integer icomm_grid                ! comm for grid topology
	integer icomm_xyz(3)              ! directional subcomms
	integer icoord(3,nproc)           ! proc grid coordinates

	real wallTime

end module

!=======================================================================
subroutine messenger_invoke()
    use computation_parameters
	use messenger
	use mpi

	!Initialise MPI
    call MPI_init (ierr)

#if (USE_COUPLER == 0)
    CFD_COMM = MPI_COMM_WORLD 
	prefix_dir = "./"
#endif

	return
end


subroutine messenger_init()
	use messenger
	use mpi

	integer	:: idims(3)
	logical	:: Lperiodic(3)
	logical :: Lremain_dims(3)
	integer	:: np, ndims, ip, ixyz

	! Initialize MPI
	call MPI_comm_size (CFD_COMM, np, ierr)
	call MPI_comm_rank (CFD_COMM, myid, ierr)
	if (np .ne. nproc) stop "Wrong number of processors"

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
	end do
	call MPI_comm_rank (icomm_xyz(1), irankx, ierr)
	call MPI_comm_rank (icomm_xyz(2), iranky, ierr)
	call MPI_comm_rank (icomm_xyz(3), irankz, ierr)
	! print *, " Old comm=",irank-1,"new x,y,z:",irankx,iranky

	! Root process at coordinates (0,0,0)
	idims = 0
	call MPI_Cart_rank(icomm_grid, idims, idroot, ierr)
	iroot = idroot + 1

	! Save current time
	wallTime = MPI_wtime()

	return
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

!=======================================================================
subroutine scatter(A, B)
	use messenger
	use mpi

	real*8 A(nx,ny,nz), B(nx_1,ny_1,nz_1)
 	real*8 sendbuf(nproc*nxyz_1)
	!MEM  real*8, allocatable :: sendbuf(:)
	integer i1(nproc), i2(nproc), &
	        j1(nproc), j2(nproc), &
	        k1(nproc), k2(nproc)
        integer L, ip

	!MEM  allocate (sendbuf(nproc*nxyz_1))

	! Scatter an array among the processors
	! including overlapping borders

	if (myid == idroot) then
		do ip=1,nproc
			i1(ip) = ibmino_1(icoord(1,ip))
			j1(ip) = jbmino_1(icoord(2,ip))
			k1(ip) = kbmino_1(icoord(3,ip))
		end do
		i2 = i1 + nx_1 - 1
		j2 = j1 + ny_1 - 1
		k2 = k1 + nz_1 - 1

		L = 0
		do ip=1,nproc
			do k=k1(ip),k2(ip)
			do j=j1(ip),j2(ip)
			do i=i1(ip),i2(ip)
				L = L + 1
				sendbuf(L) = A(i,j,k)
			end do
			end do
			end do
		end do
	end if

	call MPI_scatter (sendbuf, nxyz_1, MPI_REAL8, &
	                  B, nxyz_1, MPI_REAL8, &
	                  idroot, icomm_grid, ierr)
	!MEM  deallocate (sendbuf)

	return
end

subroutine gather(A, B)
	use messenger
	use mpi

	real*8 A(nx_1,ny_1,nz_1), B(nx,ny,nz)
 	real*8 recvbuf(nx_1,ny_1,nz_1,nproc)
	!MEM   real*8, allocatable :: recvbuf(:,:,:,:)
	integer i1(nproc), i2(nproc), &
	        j1(nproc), j2(nproc), &
	        k1(nproc), k2(nproc)
        integer ip, ii,jj,kk

	!MEM  allocate (recvbuf(nx_1,ny_1,nz_1,nproc))

	! Gather a distributed array to the root process
	! including overlapping borders

	call MPI_gather (A, nxyz_1, MPI_REAL8, &
	                 recvbuf, nxyz_1, MPI_REAL8, &
	                 idroot, CFD_COMM, ierr)

	if (myid == idroot) then
		do ip=1,nproc
			i1(ip) = ibmino_1(icoord(1,ip))
			j1(ip) = jbmino_1(icoord(2,ip))
			k1(ip) = kbmino_1(icoord(3,ip))
			i2(ip) = ibmaxo_1(icoord(1,ip))
			j2(ip) = jbmaxo_1(icoord(2,ip))
			k2(ip) = kbmaxo_1(icoord(3,ip))
		end do

		do ip=1,nproc
			do k=k1(ip),k2(ip)
			do j=j1(ip),j2(ip)
			do i=i1(ip),i2(ip)
				kk = k-k1(ip)+1
				jj = j-j1(ip)+1
				ii = i-i1(ip)+1
				B(i,j,k) = recvbuf(ii,jj,kk,ip)
			end do
			end do
			end do
		end do
	end if
	!MEM   deallocate (recvbuf)

	return
end

subroutine allGather(A, B)
	use messenger
	use mpi

	real*8 A(nx_1,ny_1,nz_1), B(nx,ny,nz)
	real*8 recvbuf(nx_1,ny_1,nz_1,nproc)
	!MEM  real*8, allocatable :: recvbuf(:,:,:,:)
	integer i1(nproc), i2(nproc), &
	        j1(nproc), j2(nproc), &
	        k1(nproc), k2(nproc)
	integer ip,ii,jj,kk

	! Gather a distributed array to the root process
	! including overlapping borders

	!MEM  allocate (recvbuf(nx_1,ny_1,nz_1,nproc))

	call MPI_allGather (A, nxyz_1, MPI_REAL8, &
	                    recvbuf, nxyz_1, MPI_REAL8, &
	                    CFD_COMM, ierr)

	do ip=1,nproc
		i1(ip) = ibmino_1(icoord(1,ip))
		j1(ip) = jbmino_1(icoord(2,ip))
		k1(ip) = kbmino_1(icoord(3,ip))
		i2(ip) = ibmaxo_1(icoord(1,ip))
		j2(ip) = jbmaxo_1(icoord(2,ip))
		k2(ip) = kbmaxo_1(icoord(3,ip))
	end do

	do ip=1,nproc
		do k=k1(ip),k2(ip)
		do j=j1(ip),j2(ip)
		do i=i1(ip),i2(ip)
			kk = k-k1(ip)+1
			jj = j-j1(ip)+1
			ii = i-i1(ip)+1
			B(i,j,k) = recvbuf(ii,jj,kk,ip)
		end do
		end do
		end do
	end do

	!MEM   deallocate (recvbuf)

	return
end

!=======================================================================
subroutine scatterXY(A, B, nk)
	use messenger
	use mpi

	real*8 A(nx,ny,nk), B(nx_1,ny_1,nk)
	real*8 sendbuf(nproc*nx_1*ny_1*nk)
	integer i1(nproc), i2(nproc), &
	        j1(nproc), j2(nproc)
 	integer L, ip

	! Scatter an XY-array among the processors
	! including overlapping borders

	if (myid == idroot) then
		do ip=1,nproc
			i1(ip) = ibmino_1(icoord(1,ip))
			j1(ip) = jbmino_1(icoord(2,ip))
		end do
		i2 = i1 + nx_1 - 1
		j2 = j1 + ny_1 - 1

		L = 0
		do ip=1,nproc
			do k=1,nk
			do j=j1(ip),j2(ip)
			do i=i1(ip),i2(ip)
				L = L + 1
				sendbuf(L) = A(i,j,k)
			end do
			end do
			end do
		end do
	end if

	call MPI_scatter (sendbuf, nx_1*ny_1*nk, MPI_REAL8, &
	                  B, nx_1*ny_1*nk, MPI_REAL8, &
	                  idroot, icomm_grid, ierr)

	return
end

subroutine gatherXY(A, B, nk)
	use messenger
	use mpi

	integer nk, ip, ii, jj
	real*8 A(nx_1,ny_1,nk), B(nx,ny,nk)
	real*8 recvbuf(nx_1,ny_1,nk,nproc)
	integer i1(nproc), i2(nproc), &
	        j1(nproc), j2(nproc)

	! Gather an distributed XY-array to the root process
	! including overlapping borders

	call MPI_gather (A, nx_1*ny_1*nk, MPI_REAL8, &
	                 recvbuf, nx_1*ny_1*nk, MPI_REAL8, &
	                 idroot, CFD_COMM, ierr)

	if (myid == idroot) then
		do ip=1,nproc
			i1(ip) = ibmino_1(icoord(1,ip))
			j1(ip) = jbmino_1(icoord(2,ip))
			i2(ip) = ibmaxo_1(icoord(1,ip))
			j2(ip) = jbmaxo_1(icoord(2,ip))
		end do

		do ip=1,nproc
			do k=1,nk
			do j=j1(ip),j2(ip)
			do i=i1(ip),i2(ip)
				jj = j-j1(ip)+1
				ii = i-i1(ip)+1
				B(i,j,k) = recvbuf(ii,jj,k,ip)
			end do
			end do
			end do
		end do
	end if

	return
end

subroutine allGatherXY(A, B, nk)
	use messenger
	use mpi

	integer nk, ip, ii,jj
	real*8 A(nx_1,ny_1,nk), B(nx,ny,nk)
	real*8 recvbuf(nx_1,ny_1,nk,nproc)
	!MEM   real*8, allocatable :: recvbuf(:,:,:,:)
	integer i1(nproc), i2(nproc), &
	        j1(nproc), j2(nproc)

	!MEM  allocate (recvbuf(nx_1,ny_1,nk,nproc))

	call updateBorder(A, nx_1,ny_1,nk, id_x, 1)
	call updateBorder(A, nx_1,ny_1,nk, id_y, 2)

	! Gather a distributed XY-array to the root process
	! including overlapping borders
	call MPI_allGather (A, nx_1*ny_1*nk, MPI_REAL8, &
	                    recvbuf, nx_1*ny_1*nk, MPI_REAL8, &
	                    CFD_COMM, ierr)

	do ip=1,nproc
		i1(ip) = ibmino_1(icoord(1,ip))
		j1(ip) = jbmino_1(icoord(2,ip))
		i2(ip) = ibmaxo_1(icoord(1,ip))
		j2(ip) = jbmaxo_1(icoord(2,ip))
	end do

	do ip=1,nproc
		do k=1,nk
		do j=j1(ip),j2(ip)
		do i=i1(ip),i2(ip)
			jj = j-j1(ip)+1
			ii = i-i1(ip)+1
			B(i,j,k) = recvbuf(ii,jj,k,ip)
		end do
		end do
		end do
	end do

	!MEM  deallocate (recvbuf)

	return
end

subroutine stat_gatherXY(A, B, nk)
	use messenger
	use mpi

	integer nk, ip, ii, jj, isendcount
	real*8 A(0:nlx+1, 0:nly+1, nk), B(0:ngx+1, 0:ngy+1, nk)
	real*8 recvbuf(0:nlx+1, 0:nly+1, nk, nproc)
	integer iblk, jblk
	integer i1_loc(nproc), j1_loc(nproc)
	integer i1(nproc), i2(nproc), &
	        j1(nproc), j2(nproc)

	! Gather an distributed XY-array to the root process
	! including overlapping borders

	isendcount = (nlx+2)*(nly+2)*nk
	call MPI_gather (A, isendcount, MPI_REAL8, &
	                 recvbuf, isendcount, MPI_REAL8, &
	                 idroot, CFD_COMM, ierr)

	if (myid == idroot) then
		do ip=1,nproc
			iblk = icoord(1,ip)
			jblk = icoord(2,ip)

			i1(ip) = iTmin_1(iblk)
			j1(ip) = jTmin_1(jblk)
			i2(ip) = iTmax_1(iblk)
			j2(ip) = jTmax_1(jblk)

			i1_loc(ip) = nox1
			j1_loc(ip) = noy1

			if (iblk.eq.1) then
				i1(ip)     = 0
				i1_loc(ip) = 0
			end if
			if (jblk.eq.1) then
				j1(ip)     = 0
				j1_loc(ip) = 0
			end if
			if (iblk.eq.npx) i2(ip) = ngx+1
			if (jblk.eq.npy) j2(ip) = ngy+1
		end do

		do ip=1,nproc
			do k=1,nk
			do j=j1(ip),j2(ip)
			do i=i1(ip),i2(ip)
				jj = j-j1(ip)+j1_loc(ip)
				ii = i-i1(ip)+i1_loc(ip)
				B(i,j,k) = recvbuf(ii,jj,k,ip)
			end do
			end do
			end do
		end do
	end if

	return
end

subroutine stat_allGatherXY(A, B, nk)
	use messenger
	use mpi

	integer nk, ip, ii,jj, isendcount
	real*8 A(0:nlx+1, 0:nly+1, nk), B(0:ngx+1, 0:ngy+1, nk)
	real*8 recvbuf(0:nlx+1, 0:nly+1, nk, nproc)
	!MEM   real*8, allocatable :: recvbuf(:,:,:,:)
	integer iblk, jblk
	integer i1_loc(nproc), j1_loc(nproc)
	integer i1(nproc), i2(nproc), &
	        j1(nproc), j2(nproc)

	!MEM  allocate (recvbuf(nx_1,ny_1,nk,nproc))

	!TAMER_NEEDED?  call updateBorder(A, nx_1,ny_1,nk, id_x, 1)
	!TAMER_NEEDED?  call updateBorder(A, nx_1,ny_1,nk, id_y, 2)

	! Gather a distributed XY-array to the root process
	! including overlapping borders
	isendcount = (nlx+2)*(nly+2)*nk
	call MPI_allGather (A, isendcount, MPI_REAL8, &
	                    recvbuf, isendcount, MPI_REAL8, &
	                    CFD_COMM, ierr)

	do ip=1,nproc
		iblk = icoord(1,ip)
		jblk = icoord(2,ip)

		i1(ip) = iTmin_1(iblk)
		j1(ip) = jTmin_1(jblk)
		i2(ip) = iTmax_1(iblk)
		j2(ip) = jTmax_1(jblk)

		i1_loc(ip) = nox1
		j1_loc(ip) = noy1

		if (iblk.eq.1) then
			i1(ip)     = 0
			i1_loc(ip) = 0
		end if
		if (jblk.eq.1) then
			j1(ip)     = 0
			j1_loc(ip) = 0
		end if
		if (iblk.eq.npx) i2(ip) = ngx+1
		if (jblk.eq.npy) j2(ip) = ngy+1
	end do

	do ip=1,nproc
		do k=1,nk
		do j=j1(ip),j2(ip)
		do i=i1(ip),i2(ip)
			jj = j-j1(ip)+j1_loc(ip)
			ii = i-i1(ip)+i1_loc(ip)
			B(i,j,k) = recvbuf(ii,jj,k,ip)
		end do
		end do
		end do
	end do

	!MEM  deallocate (recvbuf)

	return
end

!=======================================================================
!  subroutine transpose(R1, R2, isign)
!  	use messenger
!  	use mpi
!
!	real   R1(nz,nix_1,niy_1), R2(nx,niz_2,niy_2)
!	real   sendbuf(nix_1,niy_1,niz_2,nbz_2), &
!	       recvbuf(nix_1,niy_2,niz_2,nbx_1)
!
!       integer isign, ip, ii,jj,kk, ib,kb
!	integer i1(nbx_1), i2(nbx_1), &
!	        k1(nbz_2), k2(nbz_2)
!
!	logical Lremain_dims(3)
!	integer, save :: icomm = 0
!
!	! Transpose for the poisson equation
!	! A(block,block,*) ==> B(*,block,block)
!
!	if (icomm == 0) then
!		Lremain_dims(1) = .true.
!		Lremain_dims(2) = .false.
!		Lremain_dims(3) = .true.
!		call MPI_Cart_sub (icomm_grid, Lremain_dims, icomm, ierr)
!	end if
!
!	do ib=1,nbx_1
!		i1(ib) = ibmin_1(ib)
!	end do
!	do kb=1,nbz_2
!		k1(kb) = kbmin_2(kb)
!	end do
!	i2 = i1 + nix_1 - 1
!	k2 = k1 + niz_2 - 1
!
!	isendcount = nix_1*niy_1*niz_2
!	irecvcount = nix_1*niy_2*niz_2
!
!	if (isign == +1) then
!		do ip=1,npx
!			do k=k1(ip),k2(ip)
!			do j=1,niy_1
!			do i=1,nix_1
!				kk = k-k1(ip)+1
!				sendbuf(i,j,kk,ip) = R1(k,i,j)
!			end do
!			end do
!			end do
!		end do
!
!		call MPI_AllToAll (sendbuf, isendcount, MPI_REAL8, &
!		                   recvbuf, irecvcount, MPI_REAL8, &
!		                   icomm, ierr)
!
!		do ip=1,npx
!			do k=1,niz_2
!			do j=1,niy_2
!			do i=i1(ip),i2(ip)
!				ii = i-i1(ip)+1
!				R2(i-imin+1,k,j) = recvbuf(ii,j,k,ip)
!			end do
!			end do
!			end do
!		end do
!
!	else if (isign == -1) then
!		do ip=1,npx
!			do k=1,niz_2
!			do j=1,niy_2
!			do i=i1(ip),i2(ip)
!				ii = i-i1(ip)+1
!				recvbuf(ii,j,k,ip) = R2(i-imin+1,k,j)
!			end do
!			end do
!			end do
!		end do
!
!		call MPI_AllToAll (recvbuf, irecvcount, MPI_REAL8, &
!		                   sendbuf, isendcount, MPI_REAL8, &
!		                   icomm, ierr)
!
!		do ip=1,npx
!			do k=k1(ip),k2(ip)
!			do j=1,niy_1
!			do i=1,nix_1
!				kk = k-k1(ip)+1
!				R1(k,i,j) = sendbuf(i,j,kk,ip)
!			end do
!			end do
!			end do
!		end do
!
!	end if
!
!	return
!end
!
!=======================================================================
!subroutine transposeA(R1, R2, isign)
!	use messenger
!	use mpi
!
!	real   R1(nix_1,niy_1,nz), R2(nx,niz_2,niy_2)
!	real   sendbuf(nix_1,niy_1,niz_2,nbz_2), &
!	       recvbuf(nix_1,niy_2,niz_2,nbx_1)
!        integer ip, ib,kb, ii,jj,kk, isign
!	integer i1(nbx_1), i2(nbx_1), &
!	        k1(nbz_2), k2(nbz_2)
!	logical Lremain_dims(3)
!	integer, save :: icomm = 0
!
!	! Transpose for the poisson equation
!	! A(block,block,*) ==> B(*,block,block)
!
!	if (icomm == 0) then
!		Lremain_dims(1) = .true.
!		Lremain_dims(2) = .false.
!		Lremain_dims(3) = .true.
!		call MPI_Cart_sub (icomm_grid, Lremain_dims, icomm, ierr)
!	end if
!
!	do ib=1,nbx_1
!		i1(ib) = ibmin_1(ib)
!	end do
!	do kb=1,nbz_2
!		k1(kb) = kbmin_2(kb)
!	end do
!	i2 = i1 + nix_1 - 1
!	k2 = k1 + niz_2 - 1
!
!	isendcount = nix_1*niy_1*niz_2
!	irecvcount = nix_1*niy_2*niz_2
!
!	if (isign == +1) then
!		do ip=1,npx
!			do k=k1(ip),k2(ip)
!			do j=1,niy_1
!			do i=1,nix_1
!				kk = k-k1(ip)+1
!				sendbuf(i,j,kk,ip) = R1(i,j,k)
!			end do
!			end do
!			end do
!		end do
!
!		call MPI_AllToAll (sendbuf, isendcount, MPI_REAL8, &
!		                   recvbuf, irecvcount, MPI_REAL8, &
!		                   icomm, ierr)
!
!		do ip=1,npx
!			do k=1,niz_2
!			do j=1,niy_2
!			do i=i1(ip),i2(ip)
!				ii = i-i1(ip)+1
!				R2(i-imin+1,k,j) = recvbuf(ii,j,k,ip)
!			end do
!			end do
!			end do
!		end do
!
!	else if (isign == -1) then
!		do ip=1,npx
!			do k=1,niz_2
!			do j=1,niy_2
!			do i=i1(ip),i2(ip)
!				ii = i-i1(ip)+1
!				recvbuf(ii,j,k,ip) = R2(i-imin+1,k,j)
!			end do
!			end do
!			end do
!		end do
!
!		call MPI_AllToAll (recvbuf, irecvcount, MPI_REAL8, &
!		                   sendbuf, isendcount, MPI_REAL8, &
!		                   icomm, ierr)
!
!		do ip=1,npx
!			do k=k1(ip),k2(ip)
!			do j=1,niy_1
!			do i=1,nix_1
!				kk = k-k1(ip)+1
!				R1(i,j,k) = sendbuf(i,j,kk,ip)
!			end do
!			end do
!			end do
!		end do
!
!	end if
!
!	return
!end
!
!=======================================================================
!subroutine transpose_Mine(R1, R2, isign)
!    use messenger
!    include "mpif.h"
!
!    real   R1(nix_1, niy_1, nz), R2(nx,ny,niz_2)
!    real   sendbuf(nix_1,niy_1,niz_2,nbz_2), &
!           recvbuf(nix_1,niy_1,niz_2,nbz_2)
!    integer isign, ib, jb, kb, ii,jj,kk, ip
!    integer i1(nbx_1), i2(nbx_1), &
!            j1(nby_1), j2(nby_1), &
!            k1(nbz_2), k2(nbz_2)
!    logical Lremain_dims(3)
!    integer, save :: icomm = 0
!
!    ! Transpose for the poisson equation
!    ! A(*,block1,block2) <==> B(*,*,block1*block2)
!
!    ! Obtain icomm
!    if (icomm == 0) then
!        Lremain_dims(1) = .true.
!        Lremain_dims(2) = .true.
!        Lremain_dims(3) = .false.
!        call MPI_Cart_sub (icomm_grid, Lremain_dims, icomm, ierr)
!    end if
!
!    do ib=1,nbx_1
!        i1(ib) = ibmin_1(ib)
!        i2(ib) = ibmax_1(ib)
!    end do
!
!    do jb=1,nby_1
!        j1(jb) = jbmin_1(jb)
!        j2(jb) = jbmax_1(jb)
!    end do
!
!    do kb=1,nbz_2
!        k1(kb) = kbmin_2(kb)
!        k2(kb) = kbmax_2(kb)
!    end do
!
!    isendcount = nix_1*niy_1*niz_2
!    irecvcount = nix_1*niy_1*niz_2
!
!    ! Load R2
!    if (isign == +1) then
!        do ip=1,nbz_2
!            do k=k1(ip),k2(ip)
!                kk = k-k1(ip)+1
!                do j=1,niy_1
!                do i=1,nix_1
!                    sendbuf(i,j,kk,ip) = R1(i,j,k)
!                end do
!                end do
!            end do
!        end do
!
!        call MPI_AllToAll (sendbuf, isendcount, MPI_REAL8, &
!                           recvbuf, irecvcount, MPI_REAL8, &
!                           icomm, ierr)
!
!        R2 = 0.
!        do ip=1,nbz_2
!            ib = icoord(1, ip)
!            jb = icoord(2, ip)
!            do k=1,niz_2
!            do j=j1(jb),j2(jb)
!                jj = j-j1(jb)+1
!            do i=i1(ib),i2(ib)
!                ii = i-i1(ib)+1
!                R2(i,j,k) = recvbuf(ii,jj,k,ip)
!            end do
!            end do
!            end do
!        end do
!
!    ! Load R1
!    else if (isign == -1) then
!        do ip=1,nbz_2
!            ib = icoord(1, ip)
!            jb = icoord(2, ip)
!            do k=1,niz_2
!            do j=j1(jb), j2(jb)
!                jj = j-j1(jb)+1
!            do i=i1(ib), i2(ib)
!                ii = i-i1(ib)+1
!                recvbuf(ii,jj,k,ip) = R2(i,j,k)
!            end do
!            end do
!            end do
!        end do
!
!        call MPI_AllToAll (recvbuf, irecvcount, MPI_REAL8, &
!                           sendbuf, isendcount, MPI_REAL8, &
!                           icomm, ierr)
!
!        R1 = 0.0
!        do ip=1,nbz_2
!            do k=k1(ip),k2(ip)
!               kk = k-k1(ip)+1
!               do j=1,niy_1
!               do i=1,nix_1
!                   R1(i,j,k) = sendbuf(i,j,kk,ip)
!               end do
!               end do
!            end do
!        end do
!
!    end if
!
!    return
!end
!

!=======================================================================
! TJ: new transpose subroutine for 2DFFT
!	- (R1.eq.qr_gbl)
!	- (R2.eq.qT)
!	- replace nlx with ngx
!=======================================================================
subroutine transpose_qr_qT_gbl(R1, R2)
    use messenger
    include "mpif.h"

    real    R1(0:ngx,1:ngz+1,0:nly), R2(0:ngx,0:ngy,1:nlzm) 
    real   sendbuf(ngx+1,nly+1,nlzm,nbz_2), &
           recvbuf(ngx+1,nly+1,nlzm,nbz_2)
    integer ip, ib,jb,kb, ii,jj,kk
    integer i1(nbx_1), i2(nbx_1), &
            j1(nby_1), j2(nby_1), &
            k1(nbz_2), k2(nbz_2)
    integer isendcount, irecvcount
    logical Lremain_dims(3)
    integer, save :: icomm = 0

    ! Transpose for the poisson equation
    ! A(*,block1,block2) <==> B(*,*,block1*block2)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! Obtain icomm, which is a subcomm
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! TJ: MPI_Cart_sub partitions the cartesian topology generated
    !     by MPI_Cart_create into a "sub-grid" cartesian topology: 
    !
    !		- Lremain_dims(:)  = (.true.,.true.,.false.)
    !           - process topology = (npx,npy,npz)
    !
    !     Calling MPI_Cart_sub will create npz new communicators
    !     called icomm in a 2D (npx,npy) cartesian topology
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    if (icomm == 0) then
        Lremain_dims(1) = .true.
        Lremain_dims(2) = .true.
        Lremain_dims(3) = .false.
        call MPI_Cart_sub (icomm_grid, Lremain_dims, icomm, ierr)
    end if
    
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! TJ: note that (nbx_1.eq.npx), (nby_1.eq.npy), (nbz_2.eq.nproc)
    ! TJ: these loops determine the starting and ending indices
    !     in terms of global coordinates (see data.export_planes.f90) 
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    do ib=1,nbx_1
        i1(ib) = iTmin_1(ib)
        i2(ib) = iTmax_1(ib)
    end do

    do jb=1,nby_1
        j1(jb) = jTmin_1(jb)
        j2(jb) = jTmax_1(jb)
    end do

    do kb=1,nbz_2
        k1(kb) = kTmin_2(kb)
        k2(kb) = kTmax_2(kb)
    end do

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! TJ: define the size of message to be sent and received
    !	- these must be equal in size
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    isendcount = (ngx+1)*nlzm*(nly+1) !qr_gbl(0:ngx,1:ngz+1,0:nly)
    irecvcount = (ngx+1)*(nly+1)*nlzm !    qT(0:ngx:0:ngy,nlzm)

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! TJ: copy R1 (qr_gbl) into sendbuf here
    !	- note: (R1.eq.qr_bgbl)
    !	- note: (   R2.eq.qT  )
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! Load R2
        do ip=1,nbz_2		! ip=1,nproc
            do k=k1(ip),k2(ip)	! global indices
                kk = k-k1(ip)+1 ! local  indices
            do j=j1_T,j2_T
                jj = j-j1_T+1
            do i=i1_T,i2_T
                ii = i-i1_T+1
                sendbuf(ii,jj,kk,ip) = R1(i,k,j) !qr_gbl(i,k,j)
            end do
            end do
            end do
	end do

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! TJ: call MPI_AlltoAll here. 
    !	- this sends data from all to all processes
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        call MPI_AllToAll (sendbuf, isendcount, MPI_REAL8, &
                           recvbuf, irecvcount, MPI_REAL8, &
                           icomm, ierr)

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! TJ: receive data here
    !	- where does the tranpose actually take place? 
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        R2 = 0.
        do ip=1,nbz_2	       !ip = 1,nproc
            ib = icoord(1, ip) !iblock = icoord(1,irank)
            jb = icoord(2, ip) !jblock = icoord(2,irank)
            do k=1,nlzm
            do j=j1(jb),j2(jb)
                jj = j-j1(jb)+1
            do i=i1(ib),i2(ib)
                ii = i-i1(ib)+1
                R2(i,j,k) = recvbuf(ii,jj,k,ip) !qT(0:ngx,0:ngy,nlzm)
            end do
            end do
            end do
        end do

    return
end

!=======================================================================
subroutine transpose_qr_qT(R1, R2)
    use messenger
    include "mpif.h"

    real   R1(ngz+1,0:nlx,0:nly), R2(0:ngx,0:ngy,1:nlzm)
    real   sendbuf(nlx+1,nly+1,nlzm,nbz_2), &
           recvbuf(nlx+1,nly+1,nlzm,nbz_2)
    integer ip, ib,jb,kb, ii,jj,kk
    integer i1(nbx_1), i2(nbx_1), &
            j1(nby_1), j2(nby_1), &
            k1(nbz_2), k2(nbz_2)
    integer isendcount, irecvcount
    logical Lremain_dims(3)
    integer, save :: icomm = 0

    ! Transpose for the poisson equation
    ! A(*,block1,block2) <==> B(*,*,block1*block2)

    ! Obtain icomm
    if (icomm == 0) then
        Lremain_dims(1) = .true.
        Lremain_dims(2) = .true.
        Lremain_dims(3) = .false.
        call MPI_Cart_sub (icomm_grid, Lremain_dims, icomm, ierr)
    end if

    do ib=1,nbx_1
        i1(ib) = iTmin_1(ib)
        i2(ib) = iTmax_1(ib)
    end do

    do jb=1,nby_1
        j1(jb) = jTmin_1(jb)
        j2(jb) = jTmax_1(jb)
    end do

    do kb=1,nbz_2
        k1(kb) = kTmin_2(kb)
        k2(kb) = kTmax_2(kb)
    end do

    isendcount = (nlx+1)*(nly+1)*nlzm
    irecvcount = (nlx+1)*(nly+1)*nlzm

    ! Load R2
        do ip=1,nbz_2
            do k=k1(ip),k2(ip)
                kk = k-k1(ip)+1
            do j=j1_T,j2_T
                jj = j-j1_T+1
            do i=i1_T,i2_T
                ii = i-i1_T+1
                sendbuf(ii,jj,kk,ip) = R1(k,i,j)
            end do
            end do
            end do
        end do

        call MPI_AllToAll (sendbuf, isendcount, MPI_REAL8, &
                           recvbuf, irecvcount, MPI_REAL8, &
                           icomm, ierr)

        R2 = 0.
        do ip=1,nbz_2
            ib = icoord(1, ip)
            jb = icoord(2, ip)
            do k=1,nlzm
            do j=j1(jb),j2(jb)
                jj = j-j1(jb)+1
            do i=i1(ib),i2(ib)
                ii = i-i1(ib)+1
                R2(i,j,k) = recvbuf(ii,jj,k,ip)
            end do
            end do
            end do
        end do

    return
end

!=======================================================================
subroutine transpose_phat_p(R1, R2)
    use messenger
    include "mpif.h"

    real   R1(0:ngx,0:ngy,1:nlzm) , R2(0:ngz,0:nlx,0:nly)
    real   sendbuf(nlx+1,nly+1,nlzm,nbz_2), &
           recvbuf(nlx+1,nly+1,nlzm,nbz_2)
    integer ip, ib,jb,kb, ii,jj,kk
    integer i1(nbx_1), i2(nbx_1), &
            j1(nby_1), j2(nby_1), &
            k1(nbz_2), k2(nbz_2)
    integer isendcount, irecvcount
    logical Lremain_dims(3)
    integer, save :: icomm = 0

    ! Transpose for the poisson equation
    ! A(*,block1,block2) <==> B(*,*,block1*block2)

    ! Obtain icomm
    if (icomm == 0) then
        Lremain_dims(1) = .true.
        Lremain_dims(2) = .true.
        Lremain_dims(3) = .false.
        call MPI_Cart_sub (icomm_grid, Lremain_dims, icomm, ierr)
    end if

    do ib=1,nbx_1
        i1(ib) = iTmin_1(ib)
        i2(ib) = iTmax_1(ib)
    end do

    do jb=1,nby_1
        j1(jb) = jTmin_1(jb)
        j2(jb) = jTmax_1(jb)
    end do

    do kb=1,nbz_2
        k1(kb) = kTmin_2(kb)
        k2(kb) = kTmax_2(kb)
    end do

    isendcount = (nlx+1)*(nly+1)*nlzm
    irecvcount = (nlx+1)*(nly+1)*nlzm

    ! Load R1
    ! Use (nox,noy) b/c I need only ONE extra P cell for vel projection
        do ip=1,nbz_2
            ib = icoord(1, ip)
            jb = icoord(2, ip)
            do k=1,nlzm
            do j=j1(jb)-noy, j2(jb)+noy
                jj = j - (j1(jb)-noy) +1
            do i=i1(ib)-nox, i2(ib)+nox
                ii = i -(i1(ib)-nox) +1
                sendbuf(ii,jj,k,ip) = R1(i,j,k)
            end do
            end do
            end do
        end do

        call MPI_AllToAll (sendbuf, isendcount, MPI_REAL8, &
                           recvbuf, irecvcount, MPI_REAL8, &
                           icomm, ierr)
        R2 = 0.
        do ip=1,nbz_2
            do k=k1(ip),k2(ip)
                kk = k-k1(ip)+1
            do j=j1_T-noy, j2_T+noy
                jj = j - (j1_T-noy) +1
            do i=i1_T-nox, i2_T+nox
                ii = i - (i1_T-nox) +1
                R2(k,i,j) = recvbuf(ii,jj,kk,ip)
            end do
            end do
            end do
        end do

    return
end

!=======================================================================
! call transpose_uvwp(uc, ngz,nlx+1,nly, uc_t, ngx+1,ngy,nlz, 1, +/-)
!
! call transpose_uvwp(vc, ngz,nlx,nly+1, vc_t, ngx,ngy+1,nlz, 2, +/-)
!
! call transpose_uvwp(wc, ngz+1,nlx,nly, wc_t, ngx,ngy,nlz+1, 3, +/-)
!
! call transpose_uvwp(p , ngz,nlx,nly,   p_t , ngx,ngy,nlz,   4, +/-)
!=======================================================================
subroutine transpose_uvwp(R1, n1,n2,n3 , R2, m1,m2,m3, iuvwp, isign)
    use messenger
    include "mpif.h"

    integer :: n1,n2,n3, m1,m2,m3, iuvwp, isign
    real    :: R1(0:n1, 0:n2, 0:n3), R2(0:m1, 0:m2, 0:m3)
    real    :: sendbuf(n2+1, n3+1, m3+1, nbz_2), &
               recvbuf(n2+1, n3+1, m3+1, nbz_2)

    integer :: ip, ib,jb,kb, ii,jj,kk
    integer :: i1(nbx_1), i2(nbx_1), &
               j1(nby_1), j2(nby_1), &
               k1(nbz_2), k2(nbz_2)
    integer :: isendcount, irecvcount
    logical :: Lremain_dims(3)
    integer, save :: icomm = 0

    ! Transpose for I/O
    ! A(*,block1,block2) <==> B(*,*,block1*block2)

    ! Obtain icomm
    if (icomm == 0) then
        Lremain_dims(1) = .true.
        Lremain_dims(2) = .true.
        Lremain_dims(3) = .false.
        call MPI_Cart_sub (icomm_grid, Lremain_dims, icomm, ierr)
    end if

   !-- Local Limits --
    ilow = i1_T
    iup  = i2_T
    if (iblock_1.eq.1    ) ilow = 0
    if (iblock_1.eq.nbx_1) iup  = nlxb

    jlow = j1_T
    jup  = j2_T
    if (jblock_1.eq.1    ) jlow = 0
    if (jblock_1.eq.nby_1) jup  = nlyb

    klow = k1_T
    kup  = k2_T
    if (kblock_2.eq.1    ) klow = 0
    if (kblock_2.eq.nbz_2) kup  = nlzb

   !-- Global Limits --
    do ib=1,nbx_1
        i1(ib) = iTmin_1(ib)
        i2(ib) = iTmax_1(ib)
    end do
    i1(1)     = 0
    i2(nbx_1) = ngx

    do jb=1,nby_1
        j1(jb) = jTmin_1(jb)
        j2(jb) = jTmax_1(jb)
    end do
    j1(1)     = 0
    j2(nby_1) = ngy

    do kb=1,nbz_2
        k1(kb) = kTmin_2(kb)
        k2(kb) = kTmax_2(kb)
    end do
    k1(1)     = 0
    k2(nbz_2) = ngz

   !-- Edit Limits --
    select case (iuvwp)
    case (1)
            i2(nbx_1) = ngx + 1
            if (iblock_1 .eq. nbx_1)  iup = nlxb+1
    case (2)
            j2(nby_1) = ngy + 1
            if (jblock_1 .eq. nby_1)  jup = nlyb+1
    case (3)
            k2(nbz_2) = ngz + 1
            if (kblock_2 .eq. nbz_2)  kup = nlzb+1
    case (4)
    case default
            stop "Transpose_uvwp: invalid value for iuvwp"
    end select


    isendcount = (n2+1)*(n3+1)*(m3+1)
    irecvcount = isendcount

    if (isign.eq.+1) then
            !--Load R1
            do ip=1,nbz_2
                do k=k1(ip),k2(ip)
                    kk = k-k1(ip)+1
                do j=jlow,jup
                    jj = j-jlow+1
                do i=ilow,iup
                    ii = i-ilow+1
                    sendbuf(ii,jj,kk,ip) = R1(k,i,j)
                end do
                end do
                end do
            end do

            call MPI_AllToAll (sendbuf, isendcount, MPI_REAL8, &
                               recvbuf, irecvcount, MPI_REAL8, &
                               icomm, ierr)
            R2 = 0.
            do ip=1,nbz_2
                ib = icoord(1, ip)
                jb = icoord(2, ip)
                do k=klow,kup
                    kk = k-klow+1
                do j=j1(jb),j2(jb)
                    jj = j-j1(jb)+1
                do i=i1(ib),i2(ib)
                    ii = i-i1(ib)+1
                    R2(i,j,k) = recvbuf(ii,jj,kk,ip)
                end do
                end do
                end do
            end do
    else if (isign.eq.-1) then
            !--Load R2
            do ip=1,nbz_2
                ib = icoord(1, ip)
                jb = icoord(2, ip)
                do k=klow,kup
                    kk = k-klow+1
                do j=j1(jb),j2(jb)
                    jj = j-j1(jb)+1
                do i=i1(ib),i2(ib)
                    ii = i-i1(ib)+1
                    sendbuf(ii,jj,kk,ip) = R2(i,j,k)
                end do
                end do
                end do
            end do

            call MPI_AllToAll (sendbuf, isendcount, MPI_REAL8, &
                               recvbuf, irecvcount, MPI_REAL8, &
                               icomm, ierr)
            R1 = 0.
            do ip=1,nbz_2
                do k=k1(ip),k2(ip)
                    kk = k-k1(ip)+1
                do j=jlow,jup
                    jj = j-jlow+1
                do i=ilow,iup
                    ii = i-ilow+1
                    R1(k,i,j) = recvbuf(ii,jj,kk,ip)
                end do
                end do
                end do
            end do
    end if

    return
end

!=======================================================================
! call transpose_con(conx  , ngz-1, nlx-1, nly-1,
!                    conx_t,  ngx-1, ngy-1, nlz-1, 1, +/-)
!
! call transpose_con(cony  ,  ngz-1, nlx-1, nly-1,
!                    cony_t,  ngx-1, ngy-1, nlz-1, 2, +/-)
!
! call transpose_con(conz  ,  ngz-1, nlx-1, nly-1,
!                    conz_t,  ngx-1, ngy-1, nlz-1, 3, +/-)
!=======================================================================
subroutine transpose_con(R1, n1,n2,n3 , R2, m1,m2,m3, ixyz, isign)
    use messenger
    include "mpif.h"

    integer :: n1,n2,n3,  m1,m2,m3
    real    :: R1( n1 , n2 , n3 ), R2( m1 , m2 , m3 )
    real    :: sendbuf(n2, n3, m3, nbz_2), recvbuf(n2, n3, m3, nbz_2)

    integer :: ixyz, isign

    integer :: ip, ib,jb,kb, ii,jj,kk
    integer :: i1(nbx_1), i2(nbx_1), &
               j1(nby_1), j2(nby_1), &
               k1(nbz_2), k2(nbz_2)
    integer :: isendcount, irecvcount
    logical :: Lremain_dims(3)
    integer, save :: icomm = 0

    ! Transpose for the poisson equation
    ! A(*,block1,block2) <==> B(*,*,block1*block2)

    ! Obtain icomm
    if (icomm == 0) then
        Lremain_dims(1) = .true.
        Lremain_dims(2) = .true.
        Lremain_dims(3) = .false.
        call MPI_Cart_sub (icomm_grid, Lremain_dims, icomm, ierr)
    end if

   !-- Local Limits --
    ilow = i1_T
    iup  = i2_T

    jlow = j1_T
    jup  = j2_T

    klow = k1_T
    kup  = k2_T

   !-- Global Limits --
    do ib=1,nbx_1
        i1(ib) = iTmin_1(ib)
        i2(ib) = iTmax_1(ib)
    end do

    do jb=1,nby_1
        j1(jb) = jTmin_1(jb)
        j2(jb) = jTmax_1(jb)
    end do

    do kb=1,nbz_2
        k1(kb) = kTmin_2(kb)
        k2(kb) = kTmax_2(kb)
    end do


   !-- Edit Limits --
    select case (ixyz)
    case (1)
            i1(1) = 2
            if (iblock_1 .eq. 1)  ilow = 2
    case (2)
            j1(1) = 2
            if (jblock_1 .eq. 1)  jlow = 2
    case (3)

    case default
            stop "Transpose_uvwp: invalid value for iuvwp"
    end select


    isendcount = n2*n3*m3
    irecvcount = isendcount
  
    if (isign .eq. +1) then
          !--Load R1
            do ip=1,nbz_2
                do k=k1(ip),k2(ip)
                    kk = k-k1(ip)+1
                do j=jlow,jup
                    jj = j-jlow+1
                do i=ilow,iup
                    ii = i-ilow+1
                    sendbuf(ii,jj,kk,ip) = R1(k,i,j)
                end do
                end do
                end do
            end do

            call MPI_AllToAll (sendbuf, isendcount, MPI_REAL8, &
                               recvbuf, irecvcount, MPI_REAL8, &
                               icomm, ierr)
            R2 = 0.
            do ip=1,nbz_2
                ib = icoord(1, ip)
                jb = icoord(2, ip)
                do k=klow,kup
                do j=j1(jb),j2(jb)
                    jj = j-j1(jb)+1
                do i=i1(ib),i2(ib)
                    ii = i-i1(ib)+1
                    R2(i,j,k) = recvbuf(ii,jj,k,ip)
                end do
                end do
                end do
            end do
    else if (isign .eq. -1) then
           !--Load R2
            do ip=1,nbz_2
                ib = icoord(1, ip)
                jb = icoord(2, ip)
                do k=klow,kup
                do j=j1(jb),j2(jb)
                    jj = j-j1(jb)+1
                do i=i1(ib),i2(ib)
                    ii = i-i1(ib)+1
                    sendbuf(ii,jj,k,ip) = R2(i,j,k)
                end do
                end do
                end do
            end do

            call MPI_AllToAll (sendbuf, isendcount, MPI_REAL8, &
                               recvbuf, irecvcount, MPI_REAL8, &
                               icomm, ierr)
            R1 = 0.
            do ip=1,nbz_2
                do k=k1(ip),k2(ip)
                    kk = k-k1(ip)+1
                do j=jlow,jup
                    jj = j-jlow+1
                do i=ilow,iup
                    ii = i-ilow+1
                    R1(k,i,j) = recvbuf(ii,jj,kk,ip)
                end do
                end do
                end do
            end do
    end if
    return
end



subroutine updateBorder(A, n1,n2,n3, ixyz, ijk)
	use messenger
	use mpi

	real A(n1,n2,n3)
	real*8, allocatable :: buf1(:,:,:), buf2(:,:,:)
	integer istatus(MPI_STATUS_SIZE)
        integer n1,n2,n3, ixyz, ijk, no, nn, ia,ib
        integer i1,i2,j1,j2,k1,k2, icount
        integer ilower, iupper, isource, idest

	! Update overlapping borders

	i1 = imap_1(ibmin_1(iblock_1))
	i2 = imap_1(ibmax_1(iblock_1))
	j1 = jmap_1(jbmin_1(jblock_1))
	j2 = jmap_1(jbmax_1(jblock_1))
	k1 = kmap_1(kbmin_1(kblock_1))
	k2 = kmap_1(kbmax_1(kblock_1))

	select case (ixyz)
	case (1)
		no = nox
		nn = nx_1
		ia = i1
		ib = i2
		ilower = 0
		iupper = 0
		if (iblock_1 .ne. 1) ilower = 1
		if (iblock_1 .ne. nbx_1) iupper = 1

	case (2)
		no = noy
		nn = ny_1
		ia = j1
		ib = j2
		ilower = 0
		iupper = 0
		if (jblock_1 .ne. 1) ilower = 1
		if (jblock_1 .ne. nby_1) iupper = 1

	case (3)
		no = noz
		nn = nz_1
		ia = k1
		ib = k2
		! Periodic direction
		ilower = 1
		iupper = 1

	case default
		stop "updateBorder: invalid value for ixyz"
	end select

	select case (ijk)
	case (1)
		if (n1 .ne. nn) stop "updateBorder: bad value for n1"
		allocate(buf1(no,n2,n3), buf2(no,n2,n3))
		icount = no*n2*n3

		buf1 = A(ia:ia+no-1,:,:)
	case (2)
		if (n2 .ne. nn) stop "updateBorder: bad value for n2"
		allocate(buf1(n1,no,n3), buf2(n1,no,n3))
		icount = n1*no*n3

		buf1 = A(:,ia:ia+no-1,:)
	case (3)
		if (n3 .ne. nn) stop "updateBorder: bad value for n3"
		allocate(buf1(n1,n2,no), buf2(n1,n2,no))
		icount = n1*n2*no

		buf1 = A(:,:,ia:ia+no-1)
	case default
		stop "updateBorder: invalid value for ijk"
	end select

	! Send to lower neighbor
	call MPI_Cart_shift(icomm_grid, ixyz-1, -1, isource, idest, ierr)
	call MPI_sendrecv(buf1, icount, MPI_REAL8, idest, 0, &
	                  buf2, icount, MPI_REAL8, isource, 0, &
	                  icomm_grid, istatus, ierr)

	select case (ijk)
	case (1)
		if (iupper == 1) &
		A(ib+1:ib+no,:,:) = buf2
		A(ib+no+1:n1,:,:) = 0.
		buf1 = A(ib-no+1:ib,:,:)
	case (2)
		if (iupper == 1) &
		A(:,ib+1:ib+no,:) = buf2
		A(:,ib+no+1:n2,:) = 0.
		buf1 = A(:,ib-no+1:ib,:)
	case (3)
		if (iupper == 1) &
		A(:,:,ib+1:ib+no) = buf2
		A(:,:,ib+no+1:n3) = 0.
		buf1 = A(:,:,ib-no+1:ib)
	end select

	! Send to upper neighbor
	call MPI_Cart_shift(icomm_grid, ixyz-1, +1, isource, idest, ierr)
	call MPI_sendrecv(buf1, icount, MPI_REAL8, idest, 0, &
	                  buf2, icount, MPI_REAL8, isource, 0, &
	                  icomm_grid, istatus, ierr)

	select case (ijk)
	case (1)
		if (ilower == 1) &
		A(ia-no:ia-1,:,:) = buf2
	case (2)
		if (ilower == 1) &
		A(:,ia-no:ia-1,:) = buf2
	case (3)
		if (ilower == 1) &
		A(:,:,ia-no:ia-1) = buf2
	end select

	deallocate(buf1, buf2)

	return
end

!=======================================================================
! Border exchange between adjacent processor
! Prototype
! 	updateBorder_lim(A, n1,n2,n3, ixyz, ijk, ia, ib, no)
! Inputs
!	A 			- array to exchange
!   n1,n2,n3 	- size of array A  (starting from 0 to n1,n2 or n3)
!	ixyz 		- direction of exchange
!	ijk  		- array indice to exchange (Note required due to k,i,j ordering)
!	ia 			- minimum value to send (subsection of A - 0 to n1,n2 or n3)
!	ib 			- maximum value to send (subsection of A - 0 to n1,n2 or n3)
!	no 			- number of halos
!
subroutine updateBorder_lim(A, n1,n2,n3, ixyz, ijk, ia, ib, no)
	use messenger
	use mpi

	real*8 A(0:n1,0:n2,0:n3)
	real*8, allocatable :: buf1(:,:,:), buf2(:,:,:)
	integer istatus(MPI_STATUS_SIZE)
	integer n1,n2,n3, ixyz,ijk,ia,ib,no
	integer ilower, iupper, icount, idest, isource

	! Update overlapping borders

	select case (ixyz)
	case (1)
        ilower = 0
        iupper = 0
        if (iblock_1 .ne. 1) ilower = 1
        if (iblock_1 .ne. nbx_1) iupper = 1

	case (2)
        ilower = 0
        iupper = 0
        if (jblock_1 .ne. 1) ilower = 1
        if (jblock_1 .ne. nby_1) iupper = 1

	case (3)
        ! Periodic direction
        ilower = 1
        iupper = 1

	case default
	        stop "updateBorder: invalid value for ixyz"
	end select

	select case (ijk)
	case (1)
        allocate(buf1(no,0:n2,0:n3), buf2(no,0:n2,0:n3))
        icount = no*(n2+1)*(n3+1)

        buf1 = A(ia:ia+no-1,:,:)
	case (2)
        allocate(buf1(0:n1,no,0:n3), buf2(0:n1,no,0:n3))
        icount = (n1+1)*no*(n3+1)

        buf1 = A(:,ia:ia+no-1,:)
	case (3)
        allocate(buf1(0:n1,0:n2,no), buf2(0:n1,0:n2,no))
        icount = (n1+1)*(n2+1)*no

        buf1 = A(:,:,ia:ia+no-1)
	case default
	        stop "updateBorder: invalid value for ijk"
	end select

	! Send to lower neighbor
	call MPI_Cart_shift(icomm_grid, ixyz-1, -1, isource, idest, ierr)
	call MPI_sendrecv(buf1, icount, MPI_REAL8, idest, 0, &
	                    buf2, icount, MPI_REAL8, isource, 0, &
	                    icomm_grid, istatus, ierr)

	select case (ijk)
	case (1)
        if (iupper == 1) &
        A(ib+1:ib+no,:,:) = buf2
        A(ib+no+1:n1,:,:) = 0.
        buf1 = A(ib-no+1:ib,:,:)
	case (2)
        if (iupper == 1) &
        A(:,ib+1:ib+no,:) = buf2
        A(:,ib+no+1:n2,:) = 0.
        buf1 = A(:,ib-no+1:ib,:)
	case (3)
        if (iupper == 1) &
        A(:,:,ib+1:ib+no) = buf2
        A(:,:,ib+no+1:n3) = 0.
        buf1 = A(:,:,ib-no+1:ib)
	end select

	! Send to upper neighbor
	call MPI_Cart_shift(icomm_grid, ixyz-1, +1, isource, idest, ierr)
	call MPI_sendrecv(buf1, icount, MPI_REAL8, idest, 0, &
	                    buf2, icount, MPI_REAL8, isource, 0, &
	                    icomm_grid, istatus, ierr)

	select case (ijk)
	case (1)
        if (ilower == 1) &
        A(ia-no:ia-1,:,:) = buf2
	case (2)
        if (ilower == 1) &
        A(:,ia-no:ia-1,:) = buf2
	case (3)
        if (ilower == 1) &
        A(:,:,ia-no:ia-1) = buf2
	end select

	deallocate(buf1, buf2)

	return
end

!=======================================================================
subroutine updateBorder_xperiodic(flag)
	use messenger
	use mpi

	integer flag, src_dst, icount
	real*8, allocatable :: buf1(:,:,:), buf2(:,:,:)
	integer istatus(MPI_STATUS_SIZE)
	integer idims(3)

	idims = 0

	if (iblock.eq.1  ) idims(1) = npx-1
	if (iblock.eq.npx) idims(1) = 0
	idims(2) = jblock-1
	idims(3) = 0

	call MPI_Cart_rank(icomm_grid, idims, src_dst, ierr)

	if (flag.eq.1) then

	    allocate(buf1(1:ngzm,0:nly,4))
	    allocate(buf2(1:ngzm,0:nly,4))
	    icount = ngzm*(nly+1)*4

	    if (iblock.eq.1) then

		    buf1(:,:,1) = uc(1:ngzm,1,0:nly)
		    buf1(:,:,2) = uc(1:ngzm,2,0:nly)
		    buf1(:,:,3) = vc(1:ngzm,1,0:nly)
		    buf1(:,:,4) = wc(1:ngzm,1,0:nly)

			call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
					 	      buf2, icount, MPI_REAL8, src_dst, 0, &
					    	  icomm_grid, istatus, ierr)

		    urightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
		    urightbc_temp(1:ngzm,0:nly,2) = buf2(:,:,2)
		    vrightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,3)
		    wrightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,4)

	    endif

	    if (iblock.eq.npx) then

		    buf1(:,:,1) = uc(1:ngzm,nlxb-1,0:nly)
		    buf1(:,:,2) = uc(1:ngzm, nlxb, 0:nly)
		    buf1(:,:,3) = vc(1:ngzm,nlxb-1,0:nly)
		    buf1(:,:,4) = wc(1:ngzm,nlxb-1,0:nly)


			call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
				    		  buf2, icount, MPI_REAL8, src_dst, 0, &
				    		  icomm_grid, istatus, ierr)

		    uleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
		    uleftbc_temp(1:ngzm,0:nly,2) = buf2(:,:,2)
		    vleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,3)
		    wleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,4)

	    endif

	else if (flag.eq.2) then

		allocate(buf1(1:ngzm,0:nly,4))
		allocate(buf2(1:ngzm,0:nly,4))
		icount = ngzm*(nly+1)*4

		if (iblock.eq.1) then

			buf1(:,:,1) = u(1:ngzm,1,0:nly)	
			buf1(:,:,2) = u(1:ngzm,2,0:nly)	
			buf1(:,:,3) = v(1:ngzm,1,0:nly)	
			buf1(:,:,4) = w(1:ngzm,1,0:nly)	

		   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
				     buf2, icount, MPI_REAL8, src_dst, 0, &
				     icomm_grid, istatus, ierr)

			urightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
			urightbc_temp(1:ngzm,0:nly,2) = buf2(:,:,2)
			vrightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,3)
			wrightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,4)

		endif  

		if (iblock.eq.npx) then

			buf1(:,:,1) = u(1:ngzm,nlxb-1,0:nly)
			buf1(:,:,2) = u(1:ngzm, nlxb, 0:nly)
			buf1(:,:,3) = v(1:ngzm,nlxb-1,0:nly)
			buf1(:,:,4) = w(1:ngzm,nlxb-1,0:nly)

		   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
				     		 buf2, icount, MPI_REAL8, src_dst, 0, &
				     						icomm_grid, istatus, ierr)

			uleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
			uleftbc_temp(1:ngzm,0:nly,2) = buf2(:,:,2)
			vleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,3)
			wleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,4)

		endif

	else if (flag.eq.3) then

			allocate(buf1(1:ngzm,0:nly,1))
			allocate(buf2(1:ngzm,0:nly,1))
			icount = ngzm*(nly+1)*1

		if (iblock.eq.1) then

			buf1(:,:,1) = ust(1:ngzm,2,0:nly)	!this can arbitrary

		   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
				     		 buf2, icount, MPI_REAL8, src_dst, 0, &
				     						icomm_grid, istatus, ierr)

			urightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
		end if

		if (iblock.eq.npx) then
	
			buf1(:,:,1) = ust(1:ngzm,nlxb,0:nly)

		   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
				     		 buf2, icount, MPI_REAL8, src_dst, 0, &
				     						icomm_grid, istatus, ierr)
		
			uleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)	
	
		end if

	else if (flag.eq.4) then

			allocate(buf1(1:ngzm,0:nly,1))
			allocate(buf2(1:ngzm,0:nly,1))
			icount = ngzm*(nly+1)*1

		if (iblock.eq.1) then

			buf1(:,:,1) = p(1:ngzm,1,0:nly)	

		   call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
				     		 buf2, icount, MPI_REAL8, src_dst, 0, &
				     						icomm_grid, istatus, ierr)

			urightbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)
		end if

		if (iblock.eq.npx) then
			buf1(:,:,1) = p(1:ngzm,nlxb-1,0:nly)

			call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
				     		  buf2, icount, MPI_REAL8, src_dst, 0, &
				    		  icomm_grid, istatus, ierr)
		
			uleftbc_temp(1:ngzm,0:nly,1) = buf2(:,:,1)	

		end if

	endif

	deallocate(buf1)
	deallocate(buf2)

    return
end

!=======================================================================
subroutine poisson_xperiodic(buf1,buf2,nl)
        use messenger
        include "mpif.h"

        integer nl, src_dst
        real*8  buf1(nl), buf2(nl)
        integer istatus(MPI_STATUS_SIZE)
        integer idims(3)

        idims = 0

        idims(2) = jblock-1
        idims(3) = 0
        if (iblock.eq.1  ) idims(1) = npx-1
        if (iblock.eq.npx) idims(1) = 0

        call MPI_Cart_rank(icomm_grid, idims, src_dst, ierr)

        call MPI_sendrecv(buf1, nl, MPI_REAL8, src_dst, 0, &
                          buf2, nl, MPI_REAL8, src_dst, 0, &
                          icomm_grid, istatus, ierr)

     return
end

!==============================================================================
subroutine advancex_xperiodic(buf1,buf2,nl)
        use messenger
        include "mpif.h"

        integer nl, src_dst
        real*8  buf1(nl), buf2(nl)
        integer istatus(MPI_STATUS_SIZE)
        integer idims(3)

        idims = 0

        idims(2) = jblock-1
        idims(3) = 0
        if (iblock.eq.1  ) idims(1) = npx-1
        if (iblock.eq.npx) idims(1) = 0

        call MPI_Cart_rank(icomm_grid, idims, src_dst, ierr)

        call MPI_sendrecv(buf1, nl, MPI_REAL8, src_dst, 0, &
                          buf2, nl, MPI_REAL8, src_dst, 0, &
                          icomm_grid, istatus, ierr)

     return
 end

!==============================================================================
subroutine projection_xperiodic(buf1,buf2,nl)
        use messenger
        include "mpif.h"

        integer nl, src_dst
        real*8  buf1(nl), buf2(nl)
        integer istatus(MPI_STATUS_SIZE)
        integer idims(3)

        idims = 0

        idims(2) = jblock-1
        idims(3) = 0
        if (iblock.eq.1  ) idims(1) = npx-1
        if (iblock.eq.npx) idims(1) = 0

        call MPI_Cart_rank(icomm_grid, idims, src_dst, ierr)

        call MPI_sendrecv(buf1, nl, MPI_REAL8, src_dst, 0, &
                          buf2, nl, MPI_REAL8, src_dst, 0, &
                          icomm_grid, istatus, ierr)

     return
 end

!=======================================================================
subroutine updateBorder_yperiodic(flag)
        use messenger
        include "mpif.h"

        integer flag, src_dst, icount
        real*8, allocatable :: buf1(:,:,:), buf2(:,:,:)
        integer istatus(MPI_STATUS_SIZE)
        integer idims(3)

        idims = 0

        idims(1) = iblock-1
        idims(3) = 0
        if (jblock.eq.1  ) idims(2) = npy-1
        if (jblock.eq.npy) idims(2) = 0

        call MPI_Cart_rank(icomm_grid, idims, src_dst, ierr)

        if (flag.eq.1) then
            allocate(buf1(1:ngzm,0:nlx,4))
            allocate(buf2(1:ngzm,0:nlx,4))
            icount = ngzm*(nlx+1)*4
    
            if (jblock.eq.1) then
               buf1(:,:,1) = vc(1:ngzm,0:nlx,1)
               buf1(:,:,2) = vc(1:ngzm,0:nlx,2)
               buf1(:,:,3) = uc(1:ngzm,0:nlx,1)
               buf1(:,:,4) = wc(1:ngzm,0:nlx,1)
    
               call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                 buf2, icount, MPI_REAL8, src_dst, 0, &
                                 icomm_grid, istatus, ierr)
    
               vcbcn(1:ngzm,0:nlx,1) = buf2(:,:,1)
               vcbcn(1:ngzm,0:nlx,2) = buf2(:,:,2)
               ucbcn(1:ngzm,0:nlx,1) = buf2(:,:,3)
               wcbcn(1:ngzm,0:nlx,1) = buf2(:,:,4)
            end if
    
            if (jblock.eq.npy) then
               buf1(:,:,1) = vc(1:ngzm,0:nlx,nlyb-1)
               buf1(:,:,2) = vc(1:ngzm,0:nlx,nlyb  )
               buf1(:,:,3) = uc(1:ngzm,0:nlx,nlyb-1)
               buf1(:,:,4) = wc(1:ngzm,0:nlx,nlyb-1)
    
               call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                 buf2, icount, MPI_REAL8, src_dst, 0, &
                                 icomm_grid, istatus, ierr)
    
               vcbcs(1:ngzm,0:nlx,1) = buf2(:,:,1)
               vcbcs(1:ngzm,0:nlx,2) = buf2(:,:,2)
               ucbcs(1:ngzm,0:nlx,1) = buf2(:,:,3)
               wcbcs(1:ngzm,0:nlx,1) = buf2(:,:,4)
            end if

        else if (flag.eq.2) then

            allocate(buf1(1:ngzm,0:nlx,2))
            allocate(buf2(1:ngzm,0:nlx,2))
            icount = ngzm*(nlx+1)*2
     
            if (jblock.eq.1) then
               buf1(:,:,1) = v(1:ngzm,0:nlx,1)
               buf1(:,:,2) = v(1:ngzm,0:nlx,2)
     
               call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                 buf2, icount, MPI_REAL8, src_dst, 0, &
                                 icomm_grid, istatus, ierr)
     
               vcbcn(1:ngzm,0:nlx,1) = buf2(:,:,1)
               vcbcn(1:ngzm,0:nlx,2) = buf2(:,:,2)
            end if
     
            if (jblock.eq.npy) then
               buf1(:,:,1) = v(1:ngzm,0:nlx,nlyb-1)
               buf1(:,:,2) = v(1:ngzm,0:nlx,nlyb  )
     
               call MPI_sendrecv(buf1, icount, MPI_REAL8, src_dst, 0, &
                                 buf2, icount, MPI_REAL8, src_dst, 0, &
                                 icomm_grid, istatus, ierr)
     
               vcbcs(1:ngzm,0:nlx,1) = buf2(:,:,1)
               vcbcs(1:ngzm,0:nlx,2) = buf2(:,:,2)
            end if

        end if

        deallocate(buf1)
        deallocate(buf2)

     return
 end

!=======================================================================

subroutine poisson_yperiodic(buf1,buf2,nl)
        use messenger
        include "mpif.h"

        integer nl, src_dst
        real*8  buf1(nl), buf2(nl)
        integer istatus(MPI_STATUS_SIZE)
        integer idims(3)

        idims = 0

        idims(1) = iblock-1
        idims(3) = 0
        if (jblock.eq.1  ) idims(2) = npy-1
        if (jblock.eq.npy) idims(2) = 0

        call MPI_Cart_rank(icomm_grid, idims, src_dst, ierr)

        call MPI_sendrecv(buf1, nl, MPI_REAL8, src_dst, 0, &
                          buf2, nl, MPI_REAL8, src_dst, 0, &
                          icomm_grid, istatus, ierr)

     return
 end




!=======================================================================

subroutine poisson_zebra(sendbuf, recvbuf, nl, ixyz, flag)
        use messenger
        include "mpif.h"

        integer nl, ixyz, flag, src_dst
        real*8  sendbuf(nl), recvbuf(nl), buf(nl)
        integer istatus(MPI_STATUS_SIZE)
        integer ilower, iupper, icount, idest, isource

        ! Update Borders Of Zebra Scheme

	ilower = 0
	iupper = 0

	if (flag.eq.0) then
        	! Send to lower neighbor
        	call MPI_Cart_shift(icomm_grid, ixyz-1, -1, isource, idest, ierr)
		if (ixyz.eq.2   .and.  jblock_1 .ne. nby_1) iupper = 1
		if (ixyz.eq.1   .and.  iblock_1 .ne. nbx_1) iupper = 1
	else if (flag.eq.1) then
		! Send to upper neighbor
		call MPI_Cart_shift(icomm_grid, ixyz-1, +1, isource, idest, ierr)
		if (ixyz.eq.2   .and.  jblock_1 .ne.   1  ) ilower = 1
		if (ixyz.eq.1   .and.  iblock_1 .ne.   1  ) ilower = 1
	else
		print*,'Illegal value of flag = ', flag
	end if

        call MPI_sendrecv(sendbuf, nl, MPI_REAL8, idest, 0, &
                	  buf    , nl, MPI_REAL8, isource, 0, &
        	       	  icomm_grid, istatus, ierr)

	if (iupper.eq.1 .or. ilower.eq.1)    recvbuf = buf

	return
end



!=============================================================================
!SY     Update borders for zebra scheme in the poisson equation (3D multigrid)
!=============================================================================
subroutine poisson_zebra_3D(sendbuf, recvbuf, lot, nl, ixyz, flag)
        use messenger
        include "mpif.h"

        integer lot, nl, ixyz, flag, src_dst
        real*8  sendbuf(lot,nl), recvbuf(lot, nl)
        real*8, allocatable ::  buf(:,:)
        integer istatus(MPI_STATUS_SIZE)
        integer ilower, iupper, icount, idest, isource

        ! Update Borders Of Zebra Scheme

	ilower = 0
	iupper = 0

	if (flag.eq.0) then
        	! Send to lower neighbor
        	call MPI_Cart_shift(icomm_grid, ixyz-1, -1, isource, idest, ierr)
		if (ixyz.eq.2   .and.  jblock_1 .ne. nby_1) iupper = 1
		if (ixyz.eq.1   .and.  iblock_1 .ne. nbx_1) iupper = 1
	else if (flag.eq.1) then
		! Send to upper neighbor
		call MPI_Cart_shift(icomm_grid, ixyz-1, +1, isource, idest, ierr)
		if (ixyz.eq.2   .and.  jblock_1 .ne.   1  ) ilower = 1
		if (ixyz.eq.1   .and.  iblock_1 .ne.   1  ) ilower = 1
	else
		print*,'Illegal value of flag = ', flag
	end if

	allocate(buf(lot,nl))
	icount=(lot*nl)

        call MPI_sendrecv(sendbuf, icount, MPI_REAL8, idest, 0, &
                	  buf    , icount, MPI_REAL8, isource, 0, &
        	       	  icomm_grid, istatus, ierr)

	if (iupper.eq.1 .or. ilower.eq.1)    recvbuf = buf

	deallocate(buf)

	return
end

!=======================================================================

subroutine triDiagonal_mrhs(ixyz, a, b, c, r, nrhs, nl)
        use messenger
        include "mpif.h"

        integer nl, nrhs, ixyz, isize
        real*8 :: a(nl),b(nl),c(nl),r(nrhs,nl)

        call MPI_comm_size (icomm_xyz(ixyz), isize, ierr)

        if (isize == 1) then
                call triDiag_mrhs(a, b, c, r, nrhs, nl)
        else
! Changed to work on T3E (talked with Chuck 6/18/98)
!               call MPI_triDiagonal (a, b, c, r, nl, icomm_xyz(ixyz))
!               call MPI_triDiagonal (a, b, c, r, nl, icomm_xyz(ixyz),0)

                call MPI_Barrier(icomm_xyz(ixyz),ierr)
                call MPI_triDiag_mrhs(a, b, c, r, nrhs, nl, icomm_xyz(ixyz))
        end if

        return
end



subroutine triDiagonalM_mrhs(ixyz, a, b, c, r, nl, lot, nrhs)
        use messenger
        include "mpif.h"

        integer nl, lot, nrhs, ixyz, isize
        real    a(lot,nl), b(lot,nl), c(lot,nl), r(nrhs,lot,nl)

        call MPI_comm_size (icomm_xyz(ixyz), isize, ierr)

        if (isize == 1) then
                call triDiagM_mrhs(a, b, c, r, nl, lot, nrhs)
        else
! Changed to work on T3E (talked with Chuck 6/18/98)
!               call MPI_triDiagonalM (a, b, c, r, nl, lot, icomm_xyz(ixyz))
!               call MPI_triDiagonalM (a, b, c, r, nl, lot, icomm_xyz(ixyz),0)
                call MPI_triDiagM_mrhs(a, b, c, r, nl, lot, nrhs, icomm_xyz(ixyz), 0)
        end if

        return
end


!=======================================================================
subroutine triPeriodic_(ixyz, a, b, c, r, nl)
	use messenger
	use mpi

	integer nl, isize, ixyz
	real a(nl), b(nl), c(nl), r(nl)


	call MPI_comm_size (icomm_xyz(ixyz), isize, ierr)

	!-- isize should be 1 !
		call triPeriodic(a, b, c, r, nl)

	return
end


subroutine triPeriodicM_(ixyz, a, b, c, r, nl, lot)
	use messenger
	use mpi

	integer lot, nl, ixyz,isize
	real    a(lot,nl), b(lot,nl), c(lot,nl), r(lot,nl)

	call MPI_comm_size (icomm_xyz(ixyz), isize, ierr)

	if (isize == 1) then
		call triPeriodicM(a, b, c, r, nl, lot)
	else
	! Changed to work on T3E (talked with Chuck 6/18/98)
		call MPI_triPeriodicM(a, b, c, r, nl, lot, icomm_xyz(ixyz))
	end if

	return
end

subroutine messenger_lasterrorcheck
	!use messenger
	use mpi

	integer resultlen
	character*12 err_buffer

	call MPI_Error_string(ierr,err_buffer,resultlen,ierr)
	print*, err_buffer

end subroutine messenger_lasterrorcheck




!!=======================================================================
!subroutine triDiagonal_(ixyz, a, b, c, r, nl)
!       use messenger
!       include "mpif.h"
!
!        integer nl, isize, ixyz
!        real a(nl), b(nl), c(nl), r(nl)
!
!
!       call MPI_comm_size (icomm_xyz(ixyz), isize, ierr)
!
!       if (isize == 1) then
!               call triDiagonal(a, b, c, r, nl)
!       else
!  ! Changed to work on T3E (talked with Chuck 6/18/98)
!  !            call MPI_triDiagonal (a, b, c, r, nl, icomm_xyz(ixyz))
!               call MPI_triDiagonal (a, b, c, r, nl, icomm_xyz(ixyz),0)
!       end if
!
!       return
!end



!subroutine triDiagonalM_(ixyz, a, b, c, r, nl, lot)
!        use messenger
!        include "mpif.h"
!
!        integer lot, nl, ixyz,isize
!        real    a(lot,nl), b(lot,nl), c(lot,nl), r(lot,nl)
!
!        call MPI_comm_size (icomm_xyz(ixyz), isize, ierr)
!
!        if (isize == 1) then
!                call triDiagonalM (a, b, c, r, nl, lot)
!        else
!  ! Changed to work on T3E (talked with Chuck 6/18/98)
!  !             call MPI_triDiagonalM (a, b, c, r, nl, lot, icomm_xyz(ixyz))
!                call MPI_triDiagonalM (a, b, c, r, nl, lot, icomm_xyz(ixyz),0)
!        end if
!
!        return
!end


!subroutine triDiagonalLUM_(ixyz, a, b, c, work, nl, lot)
!	use messenger
!	use mpi
!
!       integer nl, ixyz, isize, lot
!	real    work(lot,2*nl+6)
!
!	call MPI_comm_size (icomm_xyz(ixyz), isize, ierr)
!
!	if (isize == 1) then
!		call triDiagonalLUM (a, b, c, nl, lot)
!	else
!		call MPI_triDiagonalLUM (a, b, c, work(1,1), work(1,1+nl), &
!		                         work(1,1+2*nl), nl, lot, icomm_xyz(ixyz))
!	end if
!
!	return
!end

!subroutine triLUSolveM_(ixyz, a, b, c, r, work, nl, lot)
!	use messenger
!	use mpi
!
!       integer ixyz, nl, lot
!	real    work(lot,2*nl+6)
!
!	call MPI_comm_size (icomm_xyz(ixyz), isize, ierr)
!
!	if (isize == 1) then
!		call triLUSolveM (a, b, c, r, nl, lot)
!	else
!		call MPI_triLUSolveM (a, b, c, work(1,1), work(1,1+nl), &
!		                      work(1,1+2*nl), r, nl, lot, icomm_xyz(ixyz))
!	end if
!
!	return
!end
!

