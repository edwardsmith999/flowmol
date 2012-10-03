program create_map
	use coupler_module
	implicit none

	call initialise
	call setup_input_and_arrays
	call get_cell_ranges_md 
	call get_overlap_blocks_cfd 
	call create_realms
	call prepare_overlap_comms
	call CPL_overlap_topology

	call gatherscatter
	
	call finalise

end program create_map

subroutine get_cell_ranges_md
	use coupler_module
	implicit none

	integer :: n
	integer :: olap_jmin_mdcoord
	integer :: ncxl, ncyl, nczl
	integer :: ncy_mdonly, ncy_md, ncyP_md
	integer, parameter :: cnull = -666 

	icPmin_md = cnull
	jcPmin_md = cnull
	kcPmin_md = cnull
	icPmax_md = cnull
	jcPmax_md = cnull
	kcPmax_md = cnull
	
	ncxl = ncx / npx_md
	do n=1,npx_md
		icPmax_md(n) = n * ncxl
		icPmin_md(n) = icPmax_md(n) - ncxl + 1
	end do	

	nczl = ncz / npz_md
	do n=1,npz_md
		kcPmax_md(n) = n * nczl
		kcPmin_md(n) = kcPmax_md(n) - nczl + 1
	end do

	ncy_md   = nint(yL_md/dy)
	ncy_mdonly = ncy_md - ncy_olap
	ncyP_md = ncy_md / npy_md
	olap_jmin_mdcoord = npy_md - floor(dble(ncy_olap)/dble(ncyP_md))	 
	do n = olap_jmin_mdcoord,npy_md
		jcPmax_md(n) = n * ncyP_md - ncy_mdonly
		jcPmin_md(n) = jcPmax_md(n) - ncyP_md + 1
		if (jcPmin_md(n).le.0) jcPmin_md(n) = 1
	end do 

	if (myid_world.eq.0) call print_map_md

contains
	subroutine print_map_md
		use coupler_module
		implicit none

		integer :: n

		print*, ''
		print*, '==========================================='
		print*, '------------ M D   M A P ------------------'
		print*, '==========================================='
		print*, 'npx_md = ', npx_md
		print*, 'ncx    = ', ncx
		print*, 'ncxl   = ', ncxl
		print*, '-------------------------------------------'
		print*, '  icoord_md   icPmin_md(n)    icPmax_md(n) '
		print*, '-------------------------------------------'
		do n=1,npx_md
			print('(1x,3i11)'), n, icPmin_md(n), icPmax_md(n)
		end do	
		print*, '-------------------------------------------'
		print*, 'npy_md     = ', npy_md
		print*, 'ncy_md     = ', ncy_md
		print*, 'ncyP_md    = ', ncyP_md 
		print*, 'ncy_olap   = ', ncy_olap
		print*, 'ncy_mdonly = ', ncy_mdonly
		print*, 'olap_jmin_mdcoord = ', olap_jmin_mdcoord
		print*, 'dy         = ', dy
		print*, '-------------------------------------------'
		print*, '  jcoord_md   jcPmin_md(n)    jcPmax_md(n) '
		print*, '-------------------------------------------'
		do n = 1,npy_md	
			print('(1x,3i11)'), n, jcPmin_md(n), jcPmax_md(n)
		end do
		print*, '-------------------------------------------'
		print*, 'npz_md = ', npz_md
		print*, 'ncz    = ', ncz
		print*, 'nczl   = ', nczl
		print*, '-------------------------------------------'
		print*, '  kcoord_md   kcPmin_md(n)    kcPmax_md(n) '
		print*, '-------------------------------------------'
		do n=1,npz_md
			print('(1x,3i11)'), n, kcPmin_md(n), kcPmax_md(n)
		end do
		print*, '-------------------------------------------'

	end subroutine print_map_md
	
end subroutine get_cell_ranges_md

subroutine get_overlap_blocks_cfd
	use coupler_module
	implicit none

	integer :: n,i,endproc
	integer :: nolapsx, nolapsy, nolapsz, nolaps
	integer, parameter :: olap_null = -666

	cfd_icoord2olap_md_icoords = olap_null
	cfd_jcoord2olap_md_jcoords = olap_null
	cfd_kcoord2olap_md_kcoords = olap_null
	
	nolapsx = npx_md/npx_cfd	
	do n = 1, npx_cfd
		do i = 1,nolapsx	
			cfd_icoord2olap_md_icoords(n,i) = (n-1)*nolapsx + i
		end do
	end do
	
	nolapsz = npz_md/npz_cfd	
	do n = 1, npz_cfd
		do i = 1,nolapsz	
			cfd_kcoord2olap_md_kcoords(n,i) = (n-1)*nolapsz + i
		end do
	end do

	nolapsy = ceiling(yL_olap/yLl_md)
	endproc = ceiling(yL_olap/yLl_cfd)
	do n = 1,endproc
		do i = 1,nolapsy
			cfd_jcoord2olap_md_jcoords(n,i) = (n-1)*nolapsy + i &
			                                  + (npy_md - nolapsy)
		end do
	end do
	
	nolaps = nolapsx*nolapsy*nolapsz	
	if(myid_world.eq.0) call print_map_cfd

contains

	subroutine print_map_cfd
		use coupler_module
		implicit none

		integer :: n,i
		
		print*, ''
		print*, '==========================================='
		print*, '------------ C F D   M A P ----------------'
		print*, '==========================================='
		print*, 'npx_cfd = ', npx_cfd
		print*, 'nolapsx = ', nolapsx
		print*, '-------------------------------------------'
		print*, '  icoord_cfd       olapmin     olapmax     ' 
		print*, '-------------------------------------------'
		do n=1,npx_cfd
			print('(1x,3i11)'), n,               &
				  cfd_icoord2olap_md_icoords(n,1),         &
				  cfd_icoord2olap_md_icoords(n,nolapsx)
		end do	
		print*, '-------------------------------------------'

		print*, 'npy_cfd = ', npy_cfd
		print*, 'nolapsy = ', nolapsy
		print*, '-------------------------------------------'
		print*, '  jcoord_cfd       olapmin     olapmax     ' 
		print*, '-------------------------------------------'
		do n=1,npy_cfd
			print('(1x,3i11)'), n,               &
				  cfd_jcoord2olap_md_jcoords(n,1),         &
				  cfd_jcoord2olap_md_jcoords(n,nolapsy)
		end do	
		print*, '-------------------------------------------'

		print*, 'npz_cfd = ', npz_cfd
		print*, 'nolapsz = ', nolapsz
		print*, '-------------------------------------------'
		print*, '  kcoord_cfd       olapmin     olapmax     ' 
		print*, '-------------------------------------------'
		do n=1,npz_cfd
			print('(1x,3i11)'), n,               &
				  cfd_kcoord2olap_md_kcoords(n,1),         &
				  cfd_kcoord2olap_md_kcoords(n,nolapsz)
		end do	
		print*, '-------------------------------------------'

	end subroutine print_map_cfd
	
end subroutine get_overlap_blocks_cfd

subroutine setup_input_and_arrays
	use coupler_module
	implicit none

	open(unit=1,file='procs',form='formatted')
		read(1,*) npx_md 
		read(1,*) npy_md 
		read(1,*) npz_md 
		read(1,*) npx_cfd 
		read(1,*) npy_cfd 
		read(1,*) npz_cfd 
	close(1)

	nproc_md  = npx_md*npy_md*npz_md
	nproc_cfd = npx_cfd*npy_cfd*npz_cfd
	nproc_world = nproc_md + nproc_cfd

	allocate(icPmin_md(npx_md))
	allocate(icPmax_md(npx_md))
	allocate(jcPmin_md(npy_md))
	allocate(jcPmax_md(npy_md))
	allocate(kcPmin_md(npz_md))
	allocate(kcPmax_md(npz_md))
	allocate(cfd_icoord2olap_md_icoords(npx_cfd,npx_md/npx_cfd))
	allocate(cfd_jcoord2olap_md_jcoords(npy_cfd,npy_md/npy_cfd))
	allocate(cfd_kcoord2olap_md_kcoords(npz_cfd,npz_md/npz_cfd))

	allocate(olap_mask(nproc_world))

	allocate(coord2rank_md(npx_md,npy_md,npz_md))
	allocate(coord2rank_cfd(npx_cfd,npy_cfd,npz_cfd))
	allocate(rank2coord_cfd(3,nproc_cfd))
	allocate(rank2coord_md(3,nproc_md))
	allocate(rank_md2rank_world(nproc_md))
	allocate(rank_cfd2rank_world(nproc_cfd))

	coord2rank_md       = 0
	coord2rank_cfd      = 0
	rank2coord_md       = 0
	rank2coord_cfd      = 0
	rank_md2rank_world  = 0
	rank_cfd2rank_world = 0

	! Inputs
	ncx = 128
	ncy = 128
	ncz = 8	
	xL_cfd = 2048.d0
	yL_cfd = 2048.d0
	zL_cfd = 128.d0
	icmin_olap = 1
	icmax_olap = ncx
	jcmin_olap = 1
	jcmax_olap = 23 
	kcmin_olap = 1
	kcmax_olap = ncz

	ncx_olap = icmax_olap - icmin_olap + 1
	ncy_olap = jcmax_olap - jcmin_olap + 1
	ncz_olap = kcmax_olap - kcmin_olap + 1

	! Calcs	
	dy = yL_cfd / ncy
	yL_md = yL_cfd / 2.d0
	yL_olap = (jcmax_olap - jcmin_olap + 1) * dy
	yLl_md = yL_md/npy_md
	yLl_cfd = yL_cfd/npy_cfd

end subroutine setup_input_and_arrays

!subroutine get_overlap_gridstretch
!implicit none
!	do j=1,ngy
!		if (in_box(x(j),y(j),z(j))) then
!			olap_jmin = j
!			do j=olap_minj+1,ngy
!				if (.not. in_box(x(j),y(j),z(j))) then
!					olap_jmax = j - 1
!					exit
!				end if
!			end do
!			exit
!		end if
!		olap_jmin = olap_null
!		olap_jmax = olap_null
!	end do
!
!end subroutine get_overlap_gridstretch

!=================================================================
!module md_proc_box
!use cpl_common
!
!	!6 planes define MD processor's domain "position" and size
!	double precision :: &
!		max_plane(3),   &
!		min_plane(3)
!	
!	max_plane(1) = picoord*xL_md
!	min_plane(1) = (picoord-1)*xL_md	
!
!	max_plane(2) = pjcoord*yL_md
!	min_plane(2) = (pjcoord-1)*yL_md	
!
!	max_plane(3) = pkcoord*zL_md
!	min_plane(3) = (pkcoord-1)*zL_md	
!
!contains
!
!	logical function in_box(x,y,z)
!		real(kind(0.d0)) :: x,y,z
!		logical in_box
!			
!		if (x.le.max_plane(1) .and. &
!		    x.gt.min_plane(1) .and. &
!		    y.le.max_plane(2) .and. &
!		    y.gt.min_plane(2) .and. &
!		    z.le.max_plane(3) .and. &
!		    z.gt.min_plane(3)) then
!			in_box = .true.
!		else
!			in_box = .false.
!		end if
!		return
!	end
!	
!end module md_proc_box
!
