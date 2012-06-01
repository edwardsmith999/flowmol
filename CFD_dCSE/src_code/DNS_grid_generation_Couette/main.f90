!
!  Tamer Zaki
!  Mechanical Engineering Department
!  Stanford University
!  tzaki@stanford.edu
!
!	MAC	:	setenv ABSOFT_RT_FLAGS -littleendian	(changes it to little_endian)
!			f90 -N113 -YEXT_NAMES=LCS -o Gen_grid.data.exe -s main.f90 mesh_tanh_stretch.f90
!
!	COMPAQ	:	f95    -r8       -o Gen_grid.data.exe  main.f90 mesh_tanh_stretch.f90	(Alpha)
!	LAMBDA	:	lf95 --dbl       -o Gen_grid.data.exe  main.f90 mesh_tanh_stretch.f90	(Lahey)
!	ALC	:	ifc -w95 -cm -r8 -o Gen_grid.data.exe  main.f90 mesh_tanh_stretch.f90	(Intel)
!	CX1	:	ifc -w95 -cm -r8 -o Gen_grid.data.exe  main.f90 mesh_tanh_stretch.f90	(Intel)
!	CX1	:	ifort        -r8 -o Gen_grid.data.exe  main.f90 mesh_tanh_stretch.f90	(Intel)
!


program main

	implicit none

	integer	:: i,j
	integer	:: iunit, ilength, ierr
	double precision :: PI, arg
	integer	:: isort

	!--------INPUT PARAMETERS------
	real	:: Lx, Ly
	integer	:: ngx, ngy, imeshY, direction
	real	:: ratio

	!------ GEOMETRY VARIABLES -------
	real :: y_lower, y_upper
	real :: dx, dy
	real, allocatable, dimension( : ) ::  x ,  y
	real, allocatable, dimension(:,:) :: xpg, ypg

	open(17, file="input.file")
	read(17, *) Lx
	read(17, *) Ly
	read(17, *) ngx
	read(17, *) ngy
	read(17, *) imeshY
	read(17, *) ratio
	read(17, *) direction
	close(17)

	!Check grid stretching is not greater than 3% change per cell
	if(imeshY .ne. 0 .and. ratio .gt. 1.d0/(0.97d0**(ngy-1))) then
		print'(a,f18.8)', "WARNING - Grid stretching greater than 3% - maximum allowable ratio with specified ngy =", 1.d0/(0.97d0**(ngy-1))
	endif

  	allocate(  x (  ngx  ) ,  y (  ngy  ) )
  	allocate( xpg(ngx,ngy) , ypg(ngx,ngy) )

	!----- X grid ------
	dx = Lx  / (ngx-1.)
	do i=1,ngx
		x(i) = (i-1.) * dx
		xpg(i,:) = x(i)
	end do

	!----- Y grid ------
	y_lower = 0.;
	y_upper = Ly;
	if (imeshY.eq.1) then
        	call tanh_stretch(y_lower,y_upper,ngy,ratio, y,direction)
		do i=1,ngx
			ypg(i,:) = y(:)
		end do

	else if (imeshY.eq.2) then
	
		!--- isort sorts the points such that y(jmin) = 0., y(ngy) = yL ---
		isort = -1
		
		!--- Define PI here ---
		PI = 2.*ACOS(0.)

		!--- Generate Y-grid ---
		do j=1,ngy
			arg = (j-1)*PI
			y(j) = isort*COS(arg/DBLE(ngy-1))+1.0
		end do

		!--- Populate 2D grid ---
		do i=1,ngx
			ypg(i,:) = y(:)
		end do

	else
		dy = Ly  / (ngy-1.)
		do j=1,ngy
			y(j) = (j-1.) * dy
			ypg(:,j) = y(j)
		end do
	end if

	!----print out y-grid to screen----
	print*,'Print out y(j=1,ngy)'
	do j=1,ngy
		print*,y(j)
	end do

	inquire(iolength=ilength) xpg
	iunit = 19
	open(iunit, file="grid.data", form="unformatted", access="direct", &
		recl=ilength, status='replace', iostat=ierr)
	if (ierr .ne. 0) stop "grid.data file can not be created"
	write(iunit, rec=1) xpg
	write(iunit, rec=2) ypg
	close(iunit)

end program
