!------------------------------------------------------------------------------
!
!				Library Functions 
!
! Used to calculate least squares fits, intergration under a line, plane line 
! intersections Heaviside function
!
! --Misc Mathematics--
! least_squares(x_interval,y,npoints,intercept, gradient)
! integrate_trap(y,npoints,x_interval,integeral)
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
!
! Written by E d W a R d S m i t h unless otherwise specified

!------------------------------------------------------------------------------

module librarymod

	! use same name for integer or double precision args versions of imaxloc
	interface imaxloc
		module procedure imaxloc_int, imaxloc_dp
	end interface

	! use same name for the same logical operation; consider BLAS
	interface magnitude
		module procedure magnitude3, magnitudeN
	end interface

	!Various Heavisides
	interface heaviside
		module procedure int_heaviside, int_array_heaviside, dp_heaviside, &
		                 dp_array_heaviside
	end interface

    private int_heaviside, int_array_heaviside, dp_heaviside, &
	        dp_array_heaviside
	
contains

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

subroutine locate(fileid,keyword,required,input_present)
use interfaces
implicit none
	
	character*(*),intent(in)	:: keyword		! Input keyword	
	integer,intent(in)		:: fileid		! File unit number
	logical,intent(in)		:: required		! Flag to check if input is required
	logical,intent(out),optional	:: input_present	! Optional flag passed back if present in intput

	character*(100)			:: linestring		! First 100 characters in a line
	integer				:: keyword_length	! Length of input keyword
	integer				:: io			! File status flag

	!Check if required, if not then check if input check variable is include
	!and terminate if missing
	if (.not. required) then
		if(present(input_present)) then
			input_present = .true.	!Assume true until not found
		else
			call error_abort( "ERROR IN LOCATE - If input not required, extra logical argument&
				& must be included to check if variable is specified in input file") 
			
		endif
	endif

	keyword_length = len(keyword)
	rewind(fileid)	! Go to beginning of input file
	do
		read (fileid,*,iostat=io) linestring			! Read first 100 characters
		if (linestring(1:keyword_length).eq.keyword) exit	! If the first characters match keyword then exit loop

		if (io.ne.0) then	! If end of file is reached
			if (.not. required) then
				!print*, keyword, ' - Not specified in input file - default value taken'
				input_present = .false.
				return
			else
				call error_abort("ERROR IN LOCATE -  Required input "//trim(keyword)//"& 
					& not found in input file. Stopping simulation.")
			endif
		end if

	end do	

end subroutine locate

!--------------------------------------------------------------------------------------
!Returns the magnitude of an 3 dimensional vector

function magnitude3(a)

	
	double precision			:: magnitude3
	double precision,dimension(3),intent(in):: a

    ! this should use BLAS library 
    ! magnitude = dnorm2(3,a,1)
 
	magnitude3 = sqrt((a(1)*a(1)+a(2)*a(2)+a(3)*a(3)))

end function magnitude3

!--------------------------------------------------------------------------------------
!Returns the magnitude of an N dimensional vector

function magnitudeN(a,n)

	
	integer,intent(in)			:: n
	double precision			:: magnitudeN
	double precision,intent(in)	:: a(n)

    ! simpler with a BLAS call
    ! magnituneN = dnorm2(n,a,1)

	magnitudeN = 0.d0

	do i=1,n
		magnitudeN = magnitudeN + a(i)*a(i)
	enddo

	magnitudeN = sqrt(magnitudeN)

end function  magnitudeN

!--------------------------------------------------------------------------------------
!Subroutine for calculating least squares straight line with a uniform interval between x values

subroutine least_squares(y,x_interval,npoints,lstsqrsinter, lstsqrsgrad)
implicit none

	integer											:: n
	integer		, intent(in)						:: npoints
	double precision								:: lstsqrsx,lstsqrsy,lstsqrsx2,lstsqrsxy
	double precision, intent(in)					:: x_interval
	double precision, intent(out)					:: lstsqrsgrad, lstsqrsinter
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

end subroutine least_squares

!-------------------------------------------------------------------------------------
!Subrountine used to intergrate a function over a uniform grid using the trapizium rule

subroutine integrate_trap(y,x_interval,npoints,s)
implicit none

	integer		, intent(in)						:: npoints
	double precision, intent(in)					:: x_interval
	double precision, intent(out)					:: s
	double precision, dimension(npoints), intent(in):: y

	integer											:: n

	!Zero integral before
	s = 0.d0

	!Trapizium Rule to calculate area under line
	s = 0.5d0 * x_interval * (y(1) + y(npoints))
	do n = 2, npoints-1
		s = s + x_interval * y(n)
	enddo

end subroutine integrate_trap

!-------------------------------------------------------------------------------------
!Returns the heaviside function for input x

function int_heaviside(x)
	implicit none

	integer							:: int_heaviside
	integer	,intent(in)				:: x

	int_heaviside = 0.5*sign(1,x)+1

end function

function int_array_heaviside(x)
	implicit none

	integer,dimension(:),intent(in)	:: x
	integer,dimension(size(x))		:: int_array_heaviside

	int_array_heaviside = 0.5*sign(1,x(:))+1

end function int_array_heaviside

function dp_heaviside(x)
	implicit none

	integer						:: dp_heaviside
	double precision,intent(in)	:: x

	dp_heaviside = ceiling(sign(0.5d0,x))
	!heaviside = 0.5*sign(1.d0,x)+1

end function dp_heaviside

function dp_array_heaviside(x)
	implicit none

	double precision,dimension(:),intent(in)	:: x
	integer,dimension(size(x))					:: dp_array_heaviside

	dp_array_heaviside = ceiling(sign(0.5d0,x(:)))

end function dp_array_heaviside

!------------------------------------------------------------------------------
! Subroutine computes the intersection of a plane and a straight line
!Inputs: 
!       normal - normal vector of the Plane 
!       p - any point that belongs to the Plane 
!       ri- end point 1
!       rj- end point 2

subroutine plane_line_intersect(intersection,normal,p,ri,rj)
implicit none

	double precision, dimension(3)				:: rij
	double precision, dimension(3), intent(in)	:: ri, rj, normal, p
	double precision, dimension(3), intent(out)	:: intersection

	rij = rj-ri

	!Use vector identity to calculate point of intersection.
	intersection = ri + (-dot_product(normal,(ri-p)) / dot_product(normal,rij))*rij

end subroutine plane_line_intersect

!--------------------------------------------------------------------------------------
!Swap column a and column b
subroutine swap(a,b)
        use interfaces
	implicit none

	double precision,dimension(:),intent(inout)	:: a, b
	double precision,dimension(size(a,1))		:: temp

	if (size(a) .ne. size(b)) call error_abort("Array sizes different in swap")

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
	use interfaces
	implicit none

	double precision,dimension(3),intent(in)	:: a, b
	double precision,dimension(3)				:: crossprod

	if (size(a) .ne. 3 .or. size(b) .ne. 3) call error_abort("Error - vectors must be 3 Dimensional for cross product")

	crossprod(1) = a(2) * b(3) - a(3) * b(2)
	crossprod(2) = a(3) * b(1) - a(1) * b(3)
	crossprod(3) = a(1) * b(2) - a(2) * b(1)

end function crossprod

!--------------------------------------------------------------------------------------
!Get an integer from fortran intrinsic 'maxloc' which returns a rank 1 array for real or integer arrrays
function imaxloc_dp(a)
	implicit none

	integer 									:: imaxloc_dp
	integer,dimension(1)						:: imax
	double precision,dimension(:),intent(in)	:: a

	imax = maxloc(a(:)) 
	imaxloc_dp = imax(1)

end function imaxloc_dp


function imaxloc_int(a)
	implicit none

	integer 						:: imaxloc_int
	integer,dimension(1)			:: imax
	integer,dimension(:),intent(in)	:: a

	imax = maxloc(a(:)) 
	imaxloc_int = imax(1)

end function imaxloc_int

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
	use interfaces
	implicit none

	double precision, dimension(:,:), intent(inout)	:: A
	integer, dimension(:), intent(out)				:: indx
	double precision								:: d
	double precision, dimension(size(A,1))			:: vv
	integer											:: j, n, imax

    ! Lapack version
    ! call DGETRF( M, N, A, LDA, IPIV, INFO )

	n = size(A,1)				!Store size of matrix
	d=1.d0 						!Row interchange flag - zero as no row interchanges yet.
	vv = maxval(abs(A),dim=2)	!Scaling information

	!Check matrix is square
	if (size(A,1) .ne. size(A,2)) call error_abort("Matrix not square in LU dcmp")
	!Check matrix is singular
	if (any(vv .eq. 0.0)) call error_abort("Singular matrix in LU dcmp")
	vv = 1.d0 / vv 			!Invert vvs

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
	implicit none

	integer						:: i, n, ii, ll
	integer, dimension(:), intent(in)		:: indx
	double precision				:: summ
	double precision, dimension(:,:), intent(in)	:: A
	double precision, dimension(:), intent(inout)	:: b

        ! Lapack version, factorisation included
        ! call DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

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
! Return unused fileunit by checking all exisiting
function get_new_fileunit() result (f)
	implicit none

	logical	:: op
	integer	:: f

	f = 1
	do 
		inquire(f,opened=op)
		if (op .eqv. .false.) exit
		f = f + 1
	enddo

end function

!--------------------------------------------------------------------------------------
!Subroutine used to get size of file
subroutine get_file_size(filename,file_size)
	implicit none

	integer							:: unit_no
	integer,intent(out)				:: file_size
	!integer,dimension(13)			:: SArray(13)
	logical							:: op
	character(len=*), intent(in)	:: filename

	!Check if unit number is used and assign unique number
	do unit_no = 1,1000
		inquire(unit_no,opened=op)
		if (op .eqv. .false.) exit
	enddo

	! simpler version using inquire
	inquire(file=filename,size=file_size)


	!Use Fstat to obtain size of file ## MAY BE COMPATABILITY ISSUE ##
	!open (unit=unit_no, file=filename)
	!call FStat(unit_no, SArray)
	!close(unit_no,status="keep")

	!file_size = SArray(8)

end subroutine get_file_size

 
!-----------------------------------------------------
! Build array of dimension ncells which maps to 3D Hilbert curve.
subroutine build_hilbert(ncells,Hcurve)
	use interfaces, only : error_abort
	implicit none

	integer,dimension(3),intent(in)					 :: ncells
	integer,dimension(:,:,:),allocatable,intent(out) :: Hcurve

	integer	:: icell,jcell,kcell,maxcells
	integer	:: hindex,i,n,m
	double precision	:: shift
	double precision,dimension(:),allocatable	:: x,y,z

	!Calculate max size of cell to determine required 3d Hilbert cube dimensions
	maxcells=maxval(ncells(:))
	!Log base 2 of max size is index used in hilbert 3 calc
	hindex = ceiling(log(dble(maxcells))/log(2.d0))
	call hilbert3(hindex,x,y,z)     !Get Hilbert Curve

	!Check if number of cells is a power of two - if not round max cell up to 
	!power of two, generate curves and discard rest of domain greater than ncells
	if (hindex-log(dble(maxcells))/log(2.d0) .gt. 0.0001d0 ) then
		 maxcells = ceiling(8.d0**(dble(hindex)/3.d0))
	endif

	!Get length of hilbert curve and loop through storing m at all cell indices 
	n=max(size(x),size(y),size(z)); m =0
	shift=0.00001d0-min(minval(x),minval(y),minval(z))      !Shift to just above zero then round up
	allocate(Hcurve(ncells(1),ncells(2),ncells(3)))
	Hcurve = -666			      !Debug values
	do i=1,n
		 !Get icell,jcell,kcell at index i
		 icell = ceiling((x(i)+shift)*maxcells)
		 jcell = ceiling((y(i)+shift)*maxcells)
		 kcell = ceiling((z(i)+shift)*maxcells)
		 !Trim off part of cubic domain outside limits of ncells
	    if (icell .le. ncells(1) .and. & 
			  jcell .le. ncells(2) .and. & 
			  kcell .le. ncells(3)) then
			  !Store value of index for icell,jcell & kcell
			  m = m + 1
		 Hcurve(icell,jcell,kcell) = m
		 endif

	enddo

	!Error catching
	if (minval(Hcurve) .eq. -666) then
		 call error_abort('Error in generating Hilbert curve - change sort flag')
	endif

contains

	!-----------------------------------------------------
	! Calculate and return the 3D Hilbert curve.
	recursive subroutine hilbert3(n,x,y,z)
		implicit none

		integer,intent(in)										:: n
		double precision,dimension(:),allocatable,intent(out)   :: x,y,z
		double precision,dimension(:),allocatable				:: xo,yo,zo


		if (n .le. 0) then
			!Bottom level of recursion - set x,y and z to zero
			allocate(x(1)); x = 0.d0
			allocate(y(1)); y = 0.d0
			allocate(z(1)); z = 0.d0
		else
			! Recursively call the fractal Hilbert curve code to concatenate 8 previous
			! smaller Hilbert curves connected by a higher level Hilbert curve
			call hilbert3(n-1,xo,yo,zo)
			allocate(x(8*size(xo)))
			allocate(y(8*size(yo)))
			allocate(z(8*size(zo)))
		    x = 0.5d0 * (/ 0.5d0+zo,  0.5d0+yo, -0.5d0+yo, -0.5d0-xo, & 
						  -0.5d0-xo, -0.5d0-yo,  0.5d0-yo,  0.5d0+zo /)
		    y = 0.5d0 * (/ 0.5d0+xo,  0.5d0+zo,  0.5d0+zo,  0.5d0+yo, & 
						  -0.5d0+yo, -0.5d0-zo, -0.5d0-zo, -0.5d0-xo /)
		    z = 0.5d0 * (/ 0.5d0+yo, -0.5d0+xo, -0.5d0+xo,  0.5d0-zo, & 
						   0.5d0-zo, -0.5d0+xo, -0.5d0+xo,  0.5d0-yo /)
		endif

	end subroutine hilbert3

end subroutine build_hilbert

!--------------------------------------------------------------------------------------
! Split an integer into three factors minimising value of each

subroutine find3factors(n,nx,ny,nz)
	implicit none

	integer, intent(in)		:: n
	integer, intent(out)	:: nx,ny,nz

	integer :: nfactors,  minval1, minval2
	integer :: i, f
	integer,dimension(1)	:: minlocation
	integer, allocatable, dimension(:) :: factors, nonzerof
	integer,parameter :: large=99999999

	!Check for sensible configuration
	nx = ceiling( dble(n)**(1.d0/3.d0))
	ny = ceiling((dble(n)/dble(nx))**(1.d0/2.d0))
	nz = ceiling( dble(n)/(dble(nx)*dble(ny)))
	if (nx * ny * nz .eq. n) return

	!Otherwise find it out the hard way
	if (n .ge. 4) then
		! Find all prime factors
		allocate(factors(n/2))
		factors = 0
		call primefactors(n,factors,nfactors)
	else
		! First 3 numbers are primes
		allocate(factors(1))
		factors = n
		nfactors = 1
	endif

	! Reduce/increase number of factors to three
	if (nfactors .eq. 1) then
		nx = factors(1); ny = 1         ; nz = 1
	elseif (nfactors .eq. 2) then
		nx = factors(1); ny = factors(2); nz = 1
	elseif (nfactors .eq. 3) then
		nx = factors(1); ny = factors(2); nz = factors(3)		
	elseif (nfactors .gt. 3) then
		allocate(nonzerof(nfactors))
		nonzerof = factors(1:nfactors)

		do while (nfactors .gt. 3)
			!Multiple two minimum values and store in 
			minlocation = minloc(nonzerof)
			minval1 = nonzerof(minlocation(1))
			nonzerof(minlocation(1)) = LARGE
			minlocation = minloc(nonzerof)
			minval2 = nonzerof(minlocation(1))
			nonzerof(minlocation(1)) = LARGE

			nonzerof(minlocation(1))=minval1*minval2
			nfactors = nfactors - 1
		enddo

		minlocation = minloc(nonzerof)
		nz = nonzerof(minlocation(1)); nonzerof(minlocation(1))=LARGE
		minlocation = minloc(nonzerof)
		ny = nonzerof(minlocation(1)); nonzerof(minlocation(1))=LARGE
		minlocation = minloc(nonzerof)
		nx = nonzerof(minlocation(1)); nonzerof(minlocation(1))=LARGE
 		
	endif

	if (n - nx*ny*nz .ne. 0) stop "ERROR in find3factors"

end subroutine find3factors

!--------------------------------------------------------------------------------------
!		SUBROUTINE TO FIND THE PRIME FACTORS OF A NUMBER
! Start with 2, check whether 2 is a factor by seeing if MOD(<input_number>,2)
! is zero. If it is zero, then 2 becomes a factor. If not, check with the next number.
! When a factor is found, divide the given number with the factor found. However,
! do not move to the next possible factor - a number can occur more than once as a factor

subroutine primefactors(num, factors, f)
	implicit none

	integer, intent(in) :: num  !input number
	integer,intent(out), dimension((num/2))::factors !array to store factors
	integer, intent(inout) :: f
	integer :: i, n

	i = 2  !eligible factor
	f = 1  !number of factors
	n = num !store input number into a temporary variable
	do
		if (mod(n,i) == 0) then !if i divides 2, it is a factor
			factors(f) = i
			f = f+1
			n = n/i
		else
			i = i+1     !not a factor. move to next number
		end if
		if (n == 1) then		
			f = f-1		!its value will be one more than the number of factors
			exit
		end if
	end do

end subroutine primefactors

!--------------------------------------------------------------------------------------
! Split an integer into three factors (works well for some but not as robust as
! prime based version above

subroutine find3factors_old(n,nx,ny,nz)
	implicit none

	integer, intent(in)		:: n
	integer, intent(out)	:: nx,ny,nz

	nx = ceiling( dble(n)**(1.d0/3.d0))
	ny = ceiling((dble(n)/dble(nx))**(1.d0/2.d0))
	nz = ceiling( dble(n)/(dble(nx)*dble(ny)))

    !check if npx*npy*npz=nproc and correct
	do while (nx * ny * nz .ne. n )
		if (nz .eq. 1 .and. ny .eq. 1) then
			nx = nx + 1
		endif
		if (nz .eq. 1 .and. ny .ne. 1) then
			nx = nx + 1
			ny = ny - 1
		endif
		if (nz .ne. 1) then
			ny = ny + 1
			nz = nz - 1
		endif
	end do

end subroutine find3factors_old

!--------------------------------------------------------------------------------------
!Write matrix in correct format
subroutine write_matrix(a)
	implicit none

	integer					:: i,j
	real, dimension(:,:)	:: a

	write(*,*)
	do i = lbound(a,1), ubound(a,1)
	    write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
	end do

end subroutine write_matrix

!--------------------------------------------------------------------------------------
! Prints formatted debug statements
subroutine printf(buf,dplaces_in)
	implicit none

	double precision,dimension(:),intent(in):: buf
	integer, intent(in), optional			:: dplaces_in

	integer				:: n,dplaces
	double precision	:: maxbuf,minbuf,order
	character*10	 	:: string
	character*42	 	:: buf_precision

	!Default number of decimal places if not supplied
	if (present(dplaces_in)) then
		if (dplaces_in .le. 9) then
			dplaces = dplaces_in
		else
			print*, 'Number of decimal places in printf if limited to 9'
			dplaces = 9 !Maximum
		endif
	else
		dplaces = 4
	endif

	!Find out required format to display maximum element in buffer
	maxbuf = maxval(buf); minbuf = minval(buf)
	maxbuf = max(maxbuf,10*abs(minbuf))	!10*Ensures extra space for minus sign
	order = 1.d0; n =1
	do while (max(maxbuf,order) .ne. order)
		order = order*10.d0
		n = n + 1
	enddo
	if (n+dplaces+2 .le. 9) then
		write(buf_precision,'(a,i1,a,i1)'), 'f',n+dplaces+2,'.', dplaces
	else
		write(buf_precision,'(a,i2,a,i1)'), 'f',n+dplaces+2,'.', dplaces
	endif

	! Build up format specifier string based on size of passed array
	string='(   ' // trim(buf_precision) // ')'
	write(string(2:4),'(i3)'), size(buf) 

	!Write formatted data 
	print(string), buf

end subroutine printf

!--------------------------------------------------------------------------------------
REAL FUNCTION random_normal()

	! Adapted from the following Fortran 77 code
	!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
	!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
	!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

	!  The function random_normal() returns a normally distributed pseudo-random
	!  number with zero mean and unit variance.

	!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
	!  and J.F. Monahan augmented with quadratic bounding curves.

	IMPLICIT NONE

	!     Local variables
	REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
	            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

	!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

	DO
	  CALL RANDOM_NUMBER(u)
	  CALL RANDOM_NUMBER(v)
	  v = 1.7156 * (v - 0.5d0)

	!     Evaluate the quadratic form
	  x = u - s
	  y = ABS(v) - t
	  q = x**2 + y*(a*y - b*x)

	!     Accept P if inside inner ellipse
	  IF (q < r1) EXIT
	!     Reject P if outside outer ellipse
	  IF (q > r2) CYCLE
	!     Reject P if outside acceptance region
	  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
	END DO

	!     Return ratio of P's coordinates as the normal deviate
	random_normal = v/u
	RETURN

END FUNCTION random_normal


end module librarymod

!======================================================================
!Ensures functions work correctly
subroutine check
	use librarymod
	implicit none

	integer						:: i,j
	integer						:: n,npoints
	integer,dimension(3)		:: indx
	double precision			:: d
	double precision			:: integral
	double precision			:: rand, x_interval, intercept, gradient
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

	!Check integration under that line

	call integrate_trap(y,x_interval,npoints,integral)

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

	deallocate(y)

end
