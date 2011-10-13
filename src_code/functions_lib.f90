!-------------------------------------------------------------------------------------
!
!				Library Functions 
! Used to calculate least squares fits, intergration under a line, plane line intersections
! Heaviside function
!
! --Misc Mathematics--
! least_squares(x_interval,y,npoints,intercept, gradient)
! intergrate_trap(y,npoints,x_interval,integeral)
! heaviside(a)
! plane_line_intersect(intersection,normal,p,ri,rj)

! --Extended Intrinsics
! swap(a,b)   	  Swap a and b
! imaxloc(A)	  Return integer location of maximum array element
! outerprod(a,b)  Outer product of two vectors
! crossprod(a,b)  Cross product of two vectors
!
! --Linear Algebra--
! LUdcmp(a,indx,d) LU decomposition
! LUbksb(a,indx,b) LU backsubstituion
!
! --File handling--
! get_file_size(filename,file_size)
! locate(fileid,keyword)
!
! --Check Routines--
! check()

!-------------------------------------------------------------------------------------


module librarymod
	interface
		function magnitude(a)
			double precision			:: magnitude
			double precision,dimension(3),intent(in):: a
		end function magnitude
	end interface
	interface
		subroutine swap(a,b)
			double precision,dimension(:),intent(inout)	:: a, b
		end subroutine swap
	end interface
	interface
		function outerprod(a,b)
			double precision,dimension(:),intent(in)	:: a, b
			double precision,dimension(size(a),size(b))	:: outerprod
		end function outerprod
	end interface
	interface
		function crossprod(a,b)
			double precision,dimension(3),intent(in)	:: a, b
			double precision,dimension(3)			:: crossprod
		end function crossprod
	end interface
	interface
		function imaxloc(A)
			integer 					:: imaxloc
			integer,dimension(1)				:: imax
			double precision,dimension(:),intent(in)	:: A
		end function imaxloc
	end interface
	interface
		function nonzero(A)
			integer 					:: nonzero
			integer,dimension(1)				:: imax
			integer,dimension(:),intent(in)			:: A
		end function nonzero
	end interface
	interface
		subroutine LUdcmp(a,indx,d)
			integer, dimension(:), intent(out)		:: indx
			double precision				:: d
			double precision, dimension(:,:), intent(inout)	:: A
		end subroutine LUdcmp
	end interface
	interface
		subroutine lubksb(a,indx,b)
			integer						:: i, n, ii, ll
			integer, dimension(:), intent(out)		:: indx
			double precision				:: summ
			double precision, dimension(:,:), intent(in)	:: A
			double precision, dimension(:), intent(inout)	:: b
		end subroutine lubksb
	end interface
end module librarymod

!--------------------------------------------------------------------------------------
!Returns the magnitude of an 3 dimensional vector

function magnitude(a)

	
	double precision			:: magnitude
	double precision,dimension(3),intent(in):: a

	magnitude = (a(1)*a(1)+a(2)*a(2)+a(3)*a(3))**0.5d0

end function magnitude

!--------------------------------------------------------------------------------------
!Returns the magnitude of an N dimensional vector

function magnitudeN(a,n)

	
	integer,intent(in)		:: n
	double precision		:: magnitudeN
	double precision,intent(in)	:: a(n)

	magnitudeN = 0.d0

	do i=1,n
		magnitudeN = magnitudeN + a(i)*a(i)
	enddo

	magnitudeN = magnitudeN**0.5d0

end function  magnitudeN

!--------------------------------------------------------------------------------------
!Subroutine for calculating least squares straight line with a uniform interval between x values

subroutine least_squares(y,x_interval,npoints,lstsqrsinter, lstsqrsgrad)
implicit none

	integer						:: n
	integer		, intent(in)			:: npoints
	double precision				:: lstsqrsx,lstsqrsy,lstsqrsx2,lstsqrsxy
	double precision, intent(in)			:: x_interval
	double precision, intent(out)			:: lstsqrsgrad, lstsqrsinter
	double precision, dimension(npoints), intent(in):: y

	!Calculate molecular velocity using least squares to fit line
	!and extrapolate down to point below overlap corresponding to continuum halo
	lstsqrsy  = 0.d0
	lstsqrsx  = 0.d0
	lstsqrsx2 = 0.d0
	lstsqrsxy = 0.d0

	do n = 1, npoints
		lstsqrsx  = lstsqrsx  + n*x_interval
		lstsqrsx2 = lstsqrsx2 + (n*x_interval)*(n*x_interval)
		lstsqrsy  = lstsqrsy  + y(n)
		lstsqrsxy = lstsqrsxy + (n*x_interval)*y(n)
	enddo

	!print*, lstsqrsx, lstsqrsx2, lstsqrsy, lstsqrsxy

	!Calculate gradient
	lstsqrsgrad = (npoints*lstsqrsxy-lstsqrsx*lstsqrsy) / &
		      (npoints*lstsqrsx2-lstsqrsx*lstsqrsx)

	!Calculate intercept
	lstsqrsinter =(lstsqrsy*lstsqrsx2 - lstsqrsx*lstsqrsxy) / &
		      (npoints*lstsqrsx2 - lstsqrsx*lstsqrsx)


	print'(a,f10.5,a,f10.5)', 'y = ', lstsqrsgrad, ' x + ', lstsqrsinter

end subroutine

!-------------------------------------------------------------------------------------
!Subrountine used to intergrate a function over a uniform grid using the trapizium rule

subroutine intergrate_trap(y,x_interval,npoints,s)
implicit none

	integer						:: n
	integer		, intent(in)			:: npoints
	double precision, intent(in)			:: x_interval
	double precision, intent(out)			:: s
	double precision, dimension(npoints), intent(in):: y

	!Zero integral before
	s = 0.d0

	!Trapizium Rule to calculate area under line
	s = 0.5d0 * x_interval * (y(1) + y(npoints))
	do n = 2, npoints-1
		s = s + x_interval * y(n)
	enddo

end subroutine

!-------------------------------------------------------------------------------------
!Returns the heaviside function for input a

function heaviside(a)

	integer				:: heaviside
	double precision,intent(in)	:: a

	heaviside = ceiling(sign(0.5d0,a))

end function heaviside

!--------------------------------------------------------------------------------------
! Subroutine computes the intersection of a plane and a straight line
!Inputs: 
!       normal - normal vector of the Plane 
!       p - any point that belongs to the Plane 
!       ri- end point 1
!       rj- end point 2

subroutine plane_line_intersect(intersection,normal,p,ri,rj)
implicit none

	double precision, dimension(3)			:: rij
	double precision, dimension(3), intent(in)	:: ri, rj, normal, p
	double precision, dimension(3), intent(out)	:: intersection

	rij = rj-ri

	!Use vector identity to calculate point of intersection.
	intersection = ri + (-dot_product(normal,(ri-p)) / dot_product(normal,rij))*rij

end subroutine plane_line_intersect

!--------------------------------------------------------------------------------------
!Swap column a and column b
subroutine swap(a,b)
	implicit none

	double precision,dimension(:),intent(inout)	:: a, b
	double precision,dimension(size(a,1))		:: temp

	if (size(a) .ne. size(b)) stop "Array sizes different in swap"

	temp = a
	a = b
	b = temp
	
end subroutine swap

!--------------------------------------------------------------------------------------
!Calculate outer product of two vectors (uses Fortran intrinsic 'spread')
function outerprod(a,b)
	implicit none

	double precision,dimension(:),intent(in)	:: a, b
	double precision,dimension(size(a),size(b))	:: outerprod

	outerprod = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
	
end function outerprod

!--------------------------------------------------------------------------------------
!Calculate cross product of two 3D vectors 
function crossprod(a,b)
	implicit none

	double precision,dimension(3),intent(in)	:: a, b
	double precision,dimension(3)			:: crossprod

	if (size(a) .ne. 3 .or. size(b) .ne. 3) stop "Error - vectors must be 3 Dimensional for cross product"

	crossprod(1) = a(2) * b(3) - a(3) * b(2)
	crossprod(2) = a(3) * b(1) - a(1) * b(3)
	crossprod(3) = a(1) * b(2) - a(2) * b(1)

end function crossprod

!--------------------------------------------------------------------------------------
!Get an integer from fortran intrinsic 'maxloc' which returns a rank 1 array
function imaxloc(a)
	implicit none

	integer 					:: imaxloc
	integer,dimension(1)				:: imax
	double precision,dimension(:),intent(in)	:: a

	imax = maxloc(a(:)) 
	imaxloc = imax(1)

end function imaxloc

!--------------------------------------------------------------------------------------
! Find location of only non-zero element in integer array and return an integer of its
! location

function nonzero(a)
	implicit none

	integer 					:: nonzero
	integer,dimension(1)				:: imax
	integer,dimension(:),intent(in)			:: a

	imax = maxloc(a(:)) 
	nonzero = imax(1)

end function nonzero

!--------------------------------------------------------------------------------------
!LU decomposition taken from Fortran 90 Numerical Recipes
!a - input NxN matrix and output as LU matrix
!indx - output vector storing row permutation

subroutine LUdcmp(A,indx,d)
	use librarymod
	implicit none

	double precision, dimension(:,:), intent(inout)	:: A
	integer, dimension(:), intent(out)		:: indx
	double precision				:: d
	double precision, dimension(size(A,1))		:: vv
	integer						:: j, n, imax

	n = size(A,1)			!Store size of matrix
	d=1.d0 				!Row interchange flag - zero as no row interchanges yet.
	vv = maxval(abs(A),dim=2)	!Scaling information

	!Check matrix is square
	if (size(A,1) .ne. size(A,2)) stop "Matrix not square in LU dcmp"
	!Check matrix is singular
	if (any(vv .eq. 0.0)) stop "Singular matrix in LU dcmp"
	vv = 1.d0 / vv 			!Invert vv

	do j=1,n
		imax=(j-1)+imaxloc(vv(j:n)*abs(A(j:n,j)))
		if (j .ne. imax) then
			call swap(A(imax,:),A(j,:))
			d = -d
			vv(imax) = vv(j)
		endif
		indx(j)= imax
		if (A(j,j) .eq. 0.d0) A(j,j) = 1.0e-20 !Replace zero with a tiny number to prevent problems
		A(j+1:n,j)   =   A(j+1:n,j)/A(j,j)	!Divide by pivot element
		A(j+1:n,j+1:n) = A(j+1:n,j+1:n)-outerprod(A(j+1:n,j),A(j,j+1:n))
	enddo

end subroutine LUdcmp

!--------------------------------------------------------------------------------------
!Linear equation solver Ax=b taken from Fortran 90 Numerical Recipes
!Matrix A must be in the form of an LU decomposition

subroutine LUbksb(A,indx,b)
	use librarymod
	implicit none

	integer						:: i, n, ii, ll
	integer, dimension(:), intent(in)		:: indx
	double precision				:: summ
	double precision, dimension(:,:), intent(in)	:: A
	double precision, dimension(:), intent(inout)	:: b

	n = size(a,1)
	ii = 0
	do i = 1,n
		ll = indx(i)
		summ=b(ll)
		b(ll) = b(i)
		if (ii .ne. 0) then		!Skip zero elements
			summ = summ - dot_product(a(i,ii:i-1),b(ii:i-1))
		elseif(summ .ne. 0.d0) then
			ii = i			!Non zero element encountered so begin sum
		endif
		b(i) = summ
	enddo
	!Perform Back substitution
	do i=n,1,-1
		b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
	enddo

end subroutine LUbksb
!--------------------------------------------------------------------------------------
!Subroutine used to get size of file
subroutine get_file_size(filename,file_size)
implicit none

	integer				:: unit_no
	integer,intent(out)		:: file_size
	integer,dimension(13)		:: SArray(13)
	logical				:: op
	character(len=*), intent(in)	:: filename

	!Check if unit number is used and assign unique number
	do unit_no = 1,1000
		inquire(unit_no,opened=op)
		if (op .eqv. .false.) exit
	enddo

	!Use Fstat to obtain size of fule ## MAY BE COMPATABILITY ISSUE ##
	open (unit=unit_no, file=filename)
	call FStat(unit_no, SArray)
	close(unit_no,status="keep")

	file_size = SArray(8)

end subroutine


!-----------------------------------------------------------------------------
! Subroutine:	locate(keyword)
! Author(s):	David Trevelyan
! Description:
!		The file opened with 'fileid' is scanned sequentially by line until the 
!		character string 'keyword' matches the beginning of what is read.
!		The file position is thus set to the line after the matched keyword.
!
!		If a keyword is not found, it is assumed that the file is
!		incomplete and the program is stopped.
!-----------------------------------------------------------------------------

subroutine locate(fileid,keyword)
implicit none
	
	character*(*)	:: keyword		! Input keyword	
	character*(100)	:: linestring		! First 100 characters in a line
	integer		:: keyword_length	! Length of input keyword
	integer		:: io			! File status flag
	integer		:: fileid		! File unit number

	keyword_length = len(keyword)
	
	rewind(fileid)	! Go to beginning of input file
	do
		read (fileid,*,iostat=io) linestring			! Read first 100 characters
		if (linestring(1:keyword_length).eq.keyword) exit	! If the first characters match keyword then exit loop

		if (io.ne.0) then	! If end of file is reached
			print*, "Input parser didn't find ", keyword," in file. Stopping simulation."
			stop 	
		end if

	end do	

end subroutine locate

!======================================================================
!Ensures functions work correctly
subroutine check
	use librarymod
	implicit none

	integer						:: i,j,np
	integer						:: n,npoints
	integer,dimension(3)				:: indx
	double precision				:: d
	double precision				:: integral
	double precision				:: rand, x_interval, intercept, gradient
	double precision,dimension(3)			:: b
	double precision,dimension(3,3)			:: a
	double precision, dimension(:), allocatable	:: y

	!Check least squares line

	x_interval = 0.1d0
	npoints = 20
	allocate(y(npoints))

	do n = 1, npoints
		call random_number(rand)
		y(n) = 2*(n * x_interval)! + (rand-0.5d0))
	enddo

	call least_squares(y, x_interval , npoints ,intercept, gradient)

	print'(a,f10.5,a,f10.5)', 'y = ', gradient, 'x + ', intercept

	print*, 'points'
	do n = 1, npoints
		print*, n * x_interval, y(n)
	enddo

	!Check intergration under that line

	call intergrate_trap(y,x_interval,npoints,integral)

	print*, 'integral', integral

	print*, 'check area under triangle', &
	 0.5d0*gradient*(npoints*x_interval)*(npoints*x_interval) + (npoints*x_interval)*intercept

	n = 3

	!do i=1,n
	!do j=1,n
	!	a(i,j) = i + (j-1)*3
	!	b(i) = i
	!enddo
	!enddo

	a = reshape(source= (/3,1,5,2,-3,-1,-5,2,4/), shape = (/ 3,3 /))
	b(1)=12 ;b(2)=-13 ; b(3)=10

	do i=1,n
	do j=1,n
		print*, 'A', i,j, a(i,j)
	enddo
	enddo

	print*, 'b', b

	!call swap(a(1,:),a(2,:))

	!print*, 'A', a
	!print*, imaxloc(a(1,:))
	!print*, outerprod(a(1,:),b)

	call LUdcmp(a,indx,d)

	do i=1,n
	do j=1,n
		print*, 'LU(A)', i,j, a(i,j)
	enddo
	enddo

	call lubksb(a,indx,b)

	print*, 'soln', b

end
