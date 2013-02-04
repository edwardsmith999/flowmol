!
!  Tamer Zaki
!  Mechanical Engineering Department
!  Stanford University
!  tzaki@stanford.edu
!
!	This code reads in a grid, and computes a symmetric version about the top plane
!
!
!	setenv ABSOFT_RT_FLAGS -littleendian	(changes it to little_endian)
!	f90 -N113 -YEXT_NAMES=LCS -o Gen_grid.data.exe -s main.f90 mesh_tanh_stretch.f90
!
!	f95 -r8 -o Gen_sym_grid.exe sym_grid.f90
!
!	ifc -w95 -cm -r8 -o Gen_sym_grid.exe main_sym_grid.f90
!
!	ifort        -r8 -o Gen_sym_grid.exe main_sym_grid.f90
!


program main

	implicit none

	integer	:: i,j
	integer	:: iunit, ilength, ierr

	!--------INPUT PARAMETERS------
	real	:: Lx, Ly
	integer	:: ngx , ngy , imeshY
	integer	:: ngx2, ngy2
	real	:: ratio

	!------ GEOMETRY VARIABLES -------
	real :: y_lower, y_upper
	real :: dx, dy
	real, allocatable, dimension( : ) ::  x ,  y
	real, allocatable, dimension(:,:) :: xpg , ypg
	real, allocatable, dimension(:,:) :: xpg2, ypg2

	open(17, file="input.file")
	read(17, *) Lx
	read(17, *) Ly
	read(17, *) ngx
	read(17, *) ngy
	read(17, *) imeshY
	read(17, *) ratio
	close(17)

	ngx2 = ngx
	ngy2 = 2*(ngy-1)+1

  	allocate(  x (  ngx  ) ,  y (  ngy  ) )
  	allocate( xpg(ngx,ngy) , ypg(ngx,ngy) )

  	allocate( xpg2( ngx2 , ngy2 ) , ypg2( ngx2 , ngy2 ) )

	!----print out y-grid to screen----
	inquire(iolength=ilength) xpg
	iunit = 19
	open(iunit, file="grid.data", form="unformatted", access="direct", &
		recl=ilength, iostat=ierr)
	if (ierr .ne. 0) stop "grid.data file can not be created"
	read(iunit, rec=1) xpg
	read(iunit, rec=2) ypg
	close(iunit)

	do j=1,ngy
		ypg2(:,j) = ypg(:,j)
	end do
	do j=ngy,ngy2
		ypg2(:, j ) = Ly + ypg(:,ngy) - ypg(:,ngy2-j+1)
	end do
	do j=1,ngy2
		xpg2(:, j ) = xpg(:,1)
	end do


	inquire(iolength=ilength) xpg2
	iunit = 13
	open(iunit, file="grid.data.2", form="unformatted", access="direct", &
		recl=ilength, iostat=ierr)
	if (ierr .ne. 0) stop "grid.data.2 file can not be created"
	write(iunit, rec=1) xpg2
	write(iunit, rec=2) ypg2
	close(iunit)

	print*,'Mirrored y mesh y(j=1,ngy)'
	do j=1,ngy2
		print'(a,i4,a,f10.5)', 'Cell ', j, ' in y is ', ypg2(1,j)
	end do

end program
