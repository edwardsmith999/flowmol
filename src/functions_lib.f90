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

module Weight_fn_mod

	! Vector shape function object of the form \vek{f}(x) = ax^2 + bx + c for 0 ≤ x ≤ 1.0  
	! Constructed subject to the following values :
	!	tBC   -- vector - top boundary condition, f(x=1) = tBC(1:3)
	!	bBC   -- vector - bottom boundary conditions, f(x=0) = bBC(1:3)
	!	sBC   -- vector - sum boundary condition, the sum of f(x) at a supplied selection of 
	!		     		  x values must satisfy \sum f(x_n) = sBC(1:3)
	!	sumx  -- scalar - sum of positions of points that must sum to sBC
	!	sumx2 -- scalar - sum of squares of positions of points that must sum to sBC	
	!	Np	  -- scalar interger - Number of points

	type :: Weight_1D_vector
		private
		real(kind(0.d0)),dimension(3) :: a,b,c

	contains

		procedure :: constructor => set_coeffs_1D_vector
		procedure :: forces  	 => get_force_1D_vector

	end type Weight_1D_vector

	! NOTE -- THIS DOESN'T WORK -- THE SUM CONSTRAINT SHOULD BE BASED ON THE POSITION IN 3D
	! SPACE AND CANNOT BE WORKED OUT FROM THE PRODUCT OF 3 x 1D SHAPE FUNCTIONS
	! 			INCLUDED FOR RECORD ONLY!!!!!!
	! Vector shape function object of the form \vek{Na}(x,y,z) =  \vek{f}(x)\vek{f}(y)\vek{f}(z)
	! where each \vek{f}(x)  = ax^2 + bx + c for 0 ≤ x ≤ 1.0  
	! Constructed subject to the following values :
	!	tBC   -- vector - top boundary condition, f(x=1) = tBC(1:3,1:3)
	!	bBC   -- vector - bottom boundary conditions, f(x=0) = bBC(1:3,1:3)
	!	sBC   -- vector - sum boundary condition, the sum of f(x) at a supplied selection of 
	!		     		  x values must satisfy \sum f(x_n) = sBC(1:3)
	!	sumx  -- scalar - sum of positions of points that must sum to sBC
	!	sumx2 -- scalar - sum of squares of positions of points that must sum to sBC	
	!	Np	  -- scalar interger - Number of points

	type :: Weight_3D_vector
		private
		type(Weight_1D_vector),dimension(3) 	:: Weight_1D

	contains

		procedure :: constructor => set_coeffs_3D_vector
		procedure :: forces  	 => get_force_3D_vector

	end type Weight_3D_vector

contains

subroutine set_coeffs_1D_vector(self, Np, sumx, sumx2, bBC, tBC, sBC)
	implicit none

	class(Weight_1D_vector) 		:: self

	integer,intent(in)						  :: Np
	real(kind(0.d0)),dimension(3),intent(in)  :: bBC, tBC, sBC
	real(kind(0.d0)),intent(in)				  :: sumx, sumx2

	!Set constant to bottom boundary condition
	self%c = bBC
	!Expression for b is given by:
	self%b = (sBC(:) - Np*self%c(:) - (tBC(:)-self%c(:))*sumx2)/(sumx - sumx2)
	!Expression for a is given by:
	self%a = tBC - self%b - self%c

end subroutine set_coeffs_1D_vector

function get_force_1D_vector(self,x) result(fx)
	implicit none

	class(Weight_1D_vector) 		:: self

	real(kind(0.d0)),dimension(:),allocatable,intent(in)	:: x

	integer													:: n
	real(kind(0.d0)),dimension(:,:),allocatable				:: fx

	!Get values of f(x_n)
	allocate(fx(size(x,1),3))
	do n = 1,size(x,1)
		fx(n,:) = self%a(:)*x(n)**2.d0 + self%b(:)*x(n) + self%c(:)
	enddo

end function get_force_1D_vector

! NOTE -- THIS DOESN'T WORK -- THE SUM CONSTRAINT SHOULD BE BASED ON THE POSITION IN 3D
! SPACE AND CANNOT BE WORKED OUT FROM THE PRODUCT OF 3 x 1D SHAPE FUNCTIONS
! 			INCLUDED FOR RECORD ONLY!!!!!!
subroutine set_coeffs_3D_vector(self, Np, sumx, sumx2, bBC, tBC, sBC)
	implicit none

	class(Weight_3D_vector) 				:: self

	integer,dimension(3),intent(in)			   :: Np
	real(kind(0.d0)),dimension(3),intent(in)   :: sumx, sumx2, sBC
	real(kind(0.d0)),dimension(3,3),intent(in) :: bBC, tBC

	integer	:: ixyz

	do ixyz=1,3
		print*, ixyz, bBC(:,ixyz)
		call self%Weight_1D(ixyz)%constructor(Np(ixyz), & 
											  sumx(ixyz), &
											  sumx2(ixyz), & 
						 			    	  bBC(:,ixyz), &
											  tBC(:,ixyz), &
											  sBC(:)**(1.d0/3.d0))
	enddo

end subroutine set_coeffs_3D_vector

function get_force_3D_vector(self,r) result(fx)
	implicit none

	class(Weight_3D_vector) 		:: self

	real(kind(0.d0)),dimension(:,:),allocatable,intent(in)	:: r

	real(kind(0.d0)),dimension(:),allocatable				:: x,y,z
	real(kind(0.d0)),dimension(:,:),allocatable				:: fx

	allocate(x(size(r,1)),y(size(r,1)),z(size(r,1)),fx(size(r,1),3))
	x = r(:,1); y = r(:,2); z = r(:,3)
	fx =  self%Weight_1D(1)%forces(x) &
		 *self%Weight_1D(2)%forces(y) &
		 *self%Weight_1D(3)%forces(z)

end function get_force_3D_vector

end module Weight_fn_mod



!*****************************************************************************
! Routine to call minpack in order to fit arbitary function

module minpack_fit_funcs_mod

    real(kind(0.d0)), dimension (:),allocatable :: xdat, ydat, zdat

    abstract interface
        subroutine mp_fn( m, n, x, fvec, iflag )
            integer ( kind = 4 ), intent(in) :: m
            integer ( kind = 4 ), intent(in) :: n
            integer ( kind = 4 ), intent(in) :: iflag
            real(kind(0.d0)), intent(in)    :: x(n)
            real(kind(0.d0)), intent(out)   :: fvec(m)
        end subroutine mp_fn
    end interface

   procedure (mp_fn),     pointer :: fn => null()

contains

! curve_fit fits an arbitary function, fn, to data in the form of xdata vs ydata
! using initial guess for coefficients p0.
! N.B. size of p0 must be consistent with unknowns in function fn.

subroutine curve_fit(fn, xdata, ydata, p0, fvec, zdata)
    implicit none

    real(kind(0.d0)), dimension(:), intent(in) :: xdata, ydata
    real(kind(0.d0)), dimension(:), intent(inout) :: p0
    real(kind(0.d0)), dimension(:), allocatable, intent(out) :: fvec

    real(kind(0.d0)), dimension(:), intent(in),optional :: zdata

    external fn
    integer :: m, n, iflag, info
    real(kind(0.d0)) :: tol

    tol = 0.00001D+00
    iflag = 1

    if (size(xdata) .ne. size(ydata)) then
        stop "Error in curve_fit -- xdata != ydata"
    endif
    m = size(xdata); n = size(p0);
    allocate(fvec(m))
    !For curve fitting, number of data points must be greater than
    ! or equal to number of unknowns
    if (m .lt. n) then
        !print*,   "Warning in curve_fit -- number of data points "// &
        !         +"must be greater than or equal to "// &
        !         +"number of unknowns. Returning zeros."
        print*,   "Warning in curve_fit -- number of data points must be greater than or equal to number of unknowns. Returning zeros."
        
        p0 = 0.d0
        fvec = 0.d0
        return
    endif

    !Set xdat, ydat and may be zdat arrays which are used in fn
    allocate(xdat(m),ydat(m))
    xdat = xdata; ydat = ydata
    if (present(zdata)) then
        if (size(xdata) .ne. size(zdata)) then
            stop "Error in curve_fit -- xdata != zdata"
        endif
        allocate(zdat(m))
        zdat = zdata
    endif

    call fn(m, n, p0, fvec, iflag )
    call lmdif1(fn, m, n, p0, fvec, tol, info)

    deallocate(xdat,ydat)
    if (present(zdata)) deallocate(zdat)

end subroutine curve_fit

! cubic_fn is a cubic function routine.
subroutine cubic_fn ( m, n, x, fvec, iflag )
    implicit none

    integer ( kind = 4 ), intent(in) :: m
    integer ( kind = 4 ), intent(in) :: n
    integer ( kind = 4 ), intent(in) :: iflag

    real(kind(0.d0)), intent(inout) :: x(n)

    real(kind(0.d0)), intent(out) :: fvec(m)

    !Cubic needs four coefficients
    if (n .ne. 4) then
        stop "Error in cubic_fn, number of coefficients 'n' should be 4"
    endif

    fvec(1:m) =   x(1) & 
                + x(2)*xdat(1:m) & 
                + x(3)*xdat(1:m)**2 & 
                + x(4)*xdat(1:m)**3 & 
                - ydat(1:m)

end subroutine cubic_fn

! cubic_fn2D is a cubic function routine.
subroutine cubic_fn2D ( m, n, x, fvec, iflag )
    implicit none

    integer ( kind = 4 ), intent(in) :: m
    integer ( kind = 4 ), intent(in) :: n
    integer ( kind = 4 ), intent(in) :: iflag

    real(kind(0.d0)), intent(inout) :: x(n)

    real(kind(0.d0)), intent(out) :: fvec(m)

    !Cubic needs four coefficients
    if (n .ne. 8) then
        stop "Error in cubic_fn, number of coefficients 'n' should be 4"
    endif

    fvec(1:m) = + x(1) & 
                + x(2)*xdat(1:m) & 
                + x(3)*xdat(1:m)**2 & 
                + x(4)*xdat(1:m)**3 & 
                + x(5) &
                + x(6)*ydat(1:m) &
                + x(7)*ydat(1:m)**2 &
                + x(8)*ydat(1:m)**3 &
                - zdat(1:m)

end subroutine cubic_fn2D

end module minpack_fit_funcs_mod





MODULE PolynomialRoots
    ! ---------------------------------------------------------------------------
    ! PURPOSE - Solve for the roots of a polynomial equation with real
    !   coefficients, up to quartic order. Retrun a code indicating the nature
    !   of the roots found.

    ! AUTHORS - Alfred H. Morris, Naval Surface Weapons Center, Dahlgren,VA
    !           William L. Davis, Naval Surface Weapons Center, Dahlgren,VA
    !           Alan Miller,  CSIRO Mathematical & Information Sciences
    !                         CLAYTON, VICTORIA, AUSTRALIA 3169
    !                         http://www.mel.dms.csiro.au/~alan
    !           Ralph L. Carmichael, Public Domain Aeronautical Software
    !                         http://www.pdas.com
    !     REVISION HISTORY
    !   DATE  VERS PERSON  STATEMENT OF CHANGES
    !    ??    1.0 AHM&WLD Original coding                     
    ! 27Feb97  1.1   AM    Converted to be compatible with ELF90
    ! 12Jul98  1.2   RLC   Module format; numerous style changes
    !  4Jan99  1.3   RLC   Made the tests for zero constant term exactly zero


      IMPLICIT NONE

      INTEGER,PARAMETER,PRIVATE:: SP=KIND(1.0), DP=KIND(1.0D0)
      REAL(DP),PARAMETER,PRIVATE:: ZERO=0.0D0, FOURTH=0.25D0, HALF=0.5D0
      REAL(DP),PARAMETER,PRIVATE:: ONE=1.0D0, TWO=2.0D0, THREE=3.0D0, FOUR=4.0D0
      COMPLEX(DP),PARAMETER,PRIVATE:: CZERO=(0.D0,0.D0)

      REAL(DP),PARAMETER,PRIVATE:: EPS=EPSILON(ONE)

      CHARACTER(LEN=*),PARAMETER,PUBLIC:: POLYROOTS_VERSION= "1.3 (4 Jan 1999)"
      INTEGER,PRIVATE:: outputCode
    !    =0 degenerate equation
    !    =1 one real root
    !    =21 two identical real roots
    !    =22 two distinct real roots
    !    =23 two complex roots
    !    =31 multiple real roots
    !    =32 one real and two complex roots
    !    =33 three distinct real roots
    !    =41
    !    =42 two real and two complex roots
    !    =43
    !    =44 four complex roots

      PRIVATE:: CubeRoot
      PUBLIC:: LinearRoot
      PRIVATE:: OneLargeTwoSmall
      PUBLIC:: QuadraticRoots
      PUBLIC:: CubicRoots
      PUBLIC:: QuarticRoots
      PUBLIC:: SolvePolynomial
    !----------------------------------------------------------------------------

      INTERFACE Swap
        MODULE PROCEDURE SwapDouble, SwapSingle
      END INTERFACE

    CONTAINS

    !+
    FUNCTION CubeRoot(x) RESULT(f)
    ! ---------------------------------------------------------------------------
    ! PURPOSE - Compute the Cube Root of a REAL(DP) number. If the argument is
    !   negative, then the cube root is also negative.

      REAL(DP),INTENT(IN) :: x
      REAL(DP):: f
    !----------------------------------------------------------------------------
      IF (x < ZERO) THEN
        f=-EXP(LOG(-x)/THREE)
      ELSE IF (x > ZERO) THEN
        f=EXP(LOG(x)/THREE)
      ELSE
        f=ZERO
      END IF
      RETURN
    END Function CubeRoot   ! ---------------------------------------------------

    !+
    SUBROUTINE LinearRoot(a, z)
    ! ---------------------------------------------------------------------------
    ! PURPOSE - COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
    !              A(1) + A(2)*Z 
    !     AND STORES THE RESULTS IN Z. It is assumed that a(2) is non-zero.
      REAL(DP),INTENT(IN),DIMENSION(:):: a
      REAL(DP),INTENT(OUT):: z
    !----------------------------------------------------------------------------
      IF (a(2)==0.0) THEN
        z=ZERO
      ELSE
        z=-a(1)/a(2)
      END IF
      RETURN
    END Subroutine LinearRoot   ! -----------------------------------------------

    !+
    SUBROUTINE OneLargeTwoSmall(a1,a2,a4,w, z)
    ! ---------------------------------------------------------------------------
    ! PURPOSE - Compute the roots of a cubic when one root, w, is known to be
    !   much larger in magnitude than the other two

      REAL(DP),INTENT(IN):: a1,a2,a4
      REAL(DP),INTENT(IN):: w
      COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z


      REAL(DP),DIMENSION(3):: aq
    !----------------------------------------------------------------------------
      aq(1)=a1
      aq(2)=a2+a1/w
      aq(3)=-a4*w
      CALL QuadraticRoots(aq, z)
      z(3)=CMPLX(w,ZERO,DP)
      
      IF (AIMAG(z(1)) == ZERO) RETURN
      z(3)=z(2)
      z(2)=z(1)
      z(1)=CMPLX(w,ZERO,DP)
      RETURN
    END Subroutine OneLargeTwoSmall   ! -----------------------------------------

    !+
    SUBROUTINE QuadraticRoots(a, z)
    ! ---------------------------------------------------------------------------
    ! PURPOSE - COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
    !              A(1) + A(2)*Z + A(3)*Z**2
    !     AND STORES THE RESULTS IN Z.  IT IS ASSUMED THAT A(3) IS NONZERO.

      REAL(DP),INTENT(IN),DIMENSION(:):: a
      COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z


      REAL(DP):: d, r, w, x, y
    !----------------------------------------------------------------------------
      IF(a(1)==0.0) THEN     ! EPS is a global module constant (private)
        z(1) = CZERO               ! one root is obviously zero
        z(2) = CMPLX(-a(2)/a(3), ZERO,DP)    ! remainder is a linear eq.
        outputCode=21   ! two identical real roots
        RETURN
      END IF

      d = a(2)*a(2) - FOUR*a(1)*a(3)             ! the discriminant
      IF (ABS(d) <= TWO*eps*a(2)*a(2)) THEN
        z(1) = CMPLX(-HALF*a(2)/a(3), ZERO, DP) ! discriminant is tiny
        z(2) = z(1)
        outputCode=22  ! two distinct real roots
        RETURN
      END IF

      r = SQRT(ABS(d))
      IF (d < ZERO) THEN
        x = -HALF*a(2)/a(3)        ! negative discriminant => roots are complex   
        y = ABS(HALF*r/a(3))
        z(1) = CMPLX(x, y, DP)
        z(2) = CMPLX(x,-y, DP)   ! its conjugate
        outputCode=23                        !  COMPLEX ROOTS
        RETURN
      END IF

      IF (a(2) /= ZERO) THEN              ! see Numerical Recipes, sec. 5.5
        w = -(a(2) + SIGN(r,a(2)))
        z(1) = CMPLX(TWO*a(1)/w,  ZERO, DP)
        z(2) = CMPLX(HALF*w/a(3), ZERO, DP)
        outputCode=22           ! two real roots
        RETURN
      END IF

      x = ABS(HALF*r/a(3))   ! a(2)=0 if you get here
      z(1) = CMPLX( x, ZERO, DP)
      z(2) = CMPLX(-x, ZERO, DP)
      outputCode=22
      RETURN
    END Subroutine QuadraticRoots   ! -------------------------------------------

    !+
    SUBROUTINE CubicRoots(a, z)
    !----------------------------------------------------------------------------
    ! PURPOSE - Compute the roots of the real polynomial
    !              A(1) + A(2)*Z + A(3)*Z**2 + A(4)*Z**3
      REAL(DP),INTENT(IN),DIMENSION(:):: a
      COMPLEX(DP),INTENT(OUT),DIMENSION(:):: z

      REAL(DP),PARAMETER:: RT3=1.7320508075689D0    ! (Sqrt(3)
      REAL (DP) :: aq(3), arg, c, cf, d, p, p1, q, q1
      REAL(DP):: r, ra, rb, rq, rt
      REAL(DP):: r1, s, sf, sq, sum, t, tol, t1, w
      REAL(DP):: w1, w2, x, x1, x2, x3, y, y1, y2, y3

    ! NOTE -   It is assumed that a(4) is non-zero. No test is made here.
    !----------------------------------------------------------------------------
      IF (a(1)==0.0) THEN
        z(1) = CZERO  ! one root is obviously zero
        CALL QuadraticRoots(a(2:4), z(2:3))   ! remaining 2 roots here
        RETURN
      END IF

      p = a(3)/(THREE*a(4))
      q = a(2)/a(4)
      r = a(1)/a(4)
      tol = FOUR*EPS

      c = ZERO
      t = a(2) - p*a(3)
      IF (ABS(t) > tol*ABS(a(2))) c = t/a(4)

      t = TWO*p*p - q
      IF (ABS(t) <= tol*ABS(q)) t = ZERO
      d = r + p*t
      IF (ABS(d) <= tol*ABS(r)) GO TO 110

    !           SET  SQ = (A(4)/S)**2 * (C**3/27 + D**2/4)

      s = MAX(ABS(a(1)), ABS(a(2)), ABS(a(3)))
      p1 = a(3)/(THREE*s)
      q1 = a(2)/s
      r1 = a(1)/s

      t1 = q - 2.25D0*p*p
      IF (ABS(t1) <= tol*ABS(q)) t1 = ZERO
      w = FOURTH*r1*r1
      w1 = HALF*p1*r1*t
      w2 = q1*q1*t1/27.0D0

      IF (w1 >= ZERO) THEN
        w = w + w1
        sq = w + w2
      ELSE IF (w2 < ZERO) THEN
        sq = w + (w1 + w2)
      ELSE
        w = w + w2
        sq = w + w1
      END IF

      IF (ABS(sq) <= tol*w) sq = ZERO
      rq = ABS(s/a(4))*SQRT(ABS(sq))
      IF (sq >= ZERO) GO TO 40

    !                   ALL ROOTS ARE REAL

      arg = ATAN2(rq, -HALF*d)
      cf = COS(arg/THREE)
      sf = SIN(arg/THREE)
      rt = SQRT(-c/THREE)
      y1 = TWO*rt*cf
      y2 = -rt*(cf + rt3*sf)
      y3 = -(d/y1)/y2

      x1 = y1 - p
      x2 = y2 - p
      x3 = y3 - p

      IF (ABS(x1) > ABS(x2)) CALL Swap(x1,x2)
      IF (ABS(x2) > ABS(x3)) CALL Swap(x2,x3)
      IF (ABS(x1) > ABS(x2)) CALL Swap(x1,x2)

      w = x3

      IF (ABS(x2) < 0.1D0*ABS(x3)) GO TO 70
      IF (ABS(x1) < 0.1D0*ABS(x2)) x1 = - (r/x3)/x2
      z(1) = CMPLX(x1, ZERO,DP)
      z(2) = CMPLX(x2, ZERO,DP)
      z(3) = CMPLX(x3, ZERO,DP)
      RETURN

    !                  REAL AND COMPLEX ROOTS

    40 ra =CubeRoot(-HALF*d - SIGN(rq,d))
      rb = -c/(THREE*ra)
      t = ra + rb
      w = -p
      x = -p
      IF (ABS(t) <= tol*ABS(ra)) GO TO 41
      w = t - p
      x = -HALF*t - p
      IF (ABS(x) <= tol*ABS(p)) x = ZERO
      41 t = ABS(ra - rb)
      y = HALF*rt3*t
      
      IF (t <= tol*ABS(ra)) GO TO 60
      IF (ABS(x) < ABS(y)) GO TO 50
      s = ABS(x)
      t = y/x
      GO TO 51
    50 s = ABS(y)
      t = x/y
    51 IF (s < 0.1D0*ABS(w)) GO TO 70
      w1 = w/s
      sum = ONE + t*t
      IF (w1*w1 < 0.01D0*sum) w = - ((r/sum)/s)/s
      z(1) = CMPLX(w, ZERO,DP)
      z(2) = CMPLX(x, y,DP)
      z(3) = CMPLX(x,-y,DP)
      RETURN

    !               AT LEAST TWO ROOTS ARE EQUAL

    60 IF (ABS(x) < ABS(w)) GO TO 61
      IF (ABS(w) < 0.1D0*ABS(x)) w = - (r/x)/x
      z(1) = CMPLX(w, ZERO,DP)
      z(2) = CMPLX(x, ZERO,DP)
      z(3) = z(2)
      RETURN
      61 IF (ABS(x) < 0.1D0*ABS(w)) GO TO 70
      z(1) = CMPLX(x, ZERO,DP)
      z(2) = z(1)
      z(3) = CMPLX(w, ZERO,DP)
      RETURN

    !     HERE W IS MUCH LARGER IN MAGNITUDE THAN THE OTHER ROOTS.
    !     AS A RESULT, THE OTHER ROOTS MAY BE EXCEEDINGLY INACCURATE
    !     BECAUSE OF ROUNDOFF ERROR.  TO DEAL WITH THIS, A QUADRATIC
    !     IS FORMED WHOSE ROOTS ARE THE SAME AS THE SMALLER ROOTS OF
    !     THE CUBIC.  THIS QUADRATIC IS THEN SOLVED.

    !     THIS CODE WAS WRITTEN BY WILLIAM L. DAVIS (NSWC).

    70 aq(1) = a(1)
      aq(2) = a(2) + a(1)/w
      aq(3) = -a(4)*w
      CALL QuadraticRoots(aq, z)
      z(3) = CMPLX(w, ZERO,DP)
      
      IF (AIMAG(z(1)) == ZERO) RETURN
      z(3) = z(2)
      z(2) = z(1)
      z(1) = CMPLX(w, ZERO,DP)
      RETURN
    !-----------------------------------------------------------------------


    !                   CASE WHEN D = 0

    110 z(1) = CMPLX(-p, ZERO,DP)
      w = SQRT(ABS(c))
      IF (c < ZERO) GO TO 120
      z(2) = CMPLX(-p, w,DP)
      z(3) = CMPLX(-p,-w,DP)
      RETURN

    120 IF (p /= ZERO) GO TO 130
      z(2) = CMPLX(w, ZERO,DP)
      z(3) = CMPLX(-w, ZERO,DP)
      RETURN

    130 x = -(p + SIGN(w,p))
      z(3) = CMPLX(x, ZERO,DP)
      t = THREE*a(1)/(a(3)*x)
      IF (ABS(p) > ABS(t)) GO TO 131
      z(2) = CMPLX(t, ZERO,DP)
      RETURN
    131 z(2) = z(1)
      z(1) = CMPLX(t, ZERO,DP)
      RETURN
    END Subroutine CubicRoots   ! -----------------------------------------------


    !+
    SUBROUTINE QuarticRoots(a,z)
    !----------------------------------------------------------------------------
    ! PURPOSE - Compute the roots of the real polynomial
    !               A(1) + A(2)*Z + ... + A(5)*Z**4

      REAL(DP), INTENT(IN)     :: a(:)
      COMPLEX(DP), INTENT(OUT) :: z(:)

      COMPLEX(DP) :: w
      REAL(DP):: b,b2, c, d, e, h, p, q, r, t
      REAL(DP),DIMENSION(4):: temp
      REAL(DP):: u, v, v1, v2, x, x1, x2, x3, y


    ! NOTE - It is assumed that a(5) is non-zero. No test is made here

    !----------------------------------------------------------------------------

      IF (a(1)==0.0) THEN
        z(1) = CZERO    !  one root is obviously zero
        CALL CubicRoots(a(2:), z(2:))
        RETURN
      END IF


      b = a(4)/(FOUR*a(5))
      c = a(3)/a(5)
      d = a(2)/a(5)
      e = a(1)/a(5)
      b2 = b*b

      p = HALF*(c - 6.0D0*b2)
      q = d - TWO*b*(c - FOUR*b2)
      r = b2*(c - THREE*b2) - b*d + e

    ! SOLVE THE RESOLVENT CUBIC EQUATION. THE CUBIC HAS AT LEAST ONE
    ! NONNEGATIVE REAL ROOT.  IF W1, W2, W3 ARE THE ROOTS OF THE CUBIC
    ! THEN THE ROOTS OF THE ORIGINIAL EQUATION ARE
    !     Z = -B + CSQRT(W1) + CSQRT(W2) + CSQRT(W3)
    ! WHERE THE SIGNS OF THE SQUARE ROOTS ARE CHOSEN SO
    ! THAT CSQRT(W1) * CSQRT(W2) * CSQRT(W3) = -Q/8.

      temp(1) = -q*q/64.0D0
      temp(2) = 0.25D0*(p*p - r)
      temp(3) =  p
      temp(4) = ONE
      CALL CubicRoots(temp,z)
      IF (AIMAG(z(2)) /= ZERO) GO TO 60

    !         THE RESOLVENT CUBIC HAS ONLY REAL ROOTS
    !         REORDER THE ROOTS IN INCREASING ORDER

      x1 = DBLE(z(1))
      x2 = DBLE(z(2))
      x3 = DBLE(z(3))
      IF (x1 > x2) CALL Swap(x1,x2)
      IF (x2 > x3) CALL Swap(x2,x3)
      IF (x1 > x2) CALL Swap(x1,x2)

      u = ZERO
      IF (x3 > ZERO) u = SQRT(x3)
      IF (x2 <= ZERO) GO TO 41
      IF (x1 >= ZERO) GO TO 30
      IF (ABS(x1) > x2) GO TO 40
      x1 = ZERO

    30 x1 = SQRT(x1)
      x2 = SQRT(x2)
      IF (q > ZERO) x1 = -x1
      temp(1) = (( x1 + x2) + u) - b
      temp(2) = ((-x1 - x2) + u) - b
      temp(3) = (( x1 - x2) - u) - b
      temp(4) = ((-x1 + x2) - u) - b
      CALL SelectSort(temp)
      IF (ABS(temp(1)) >= 0.1D0*ABS(temp(4))) GO TO 31
      t = temp(2)*temp(3)*temp(4)
      IF (t /= ZERO) temp(1) = e/t
    31 z(1) = CMPLX(temp(1), ZERO,DP)
      z(2) = CMPLX(temp(2), ZERO,DP)
      z(3) = CMPLX(temp(3), ZERO,DP)
      z(4) = CMPLX(temp(4), ZERO,DP)
      RETURN

    40 v1 = SQRT(ABS(x1))
    v2 = ZERO
    GO TO 50
    41 v1 = SQRT(ABS(x1))
    v2 = SQRT(ABS(x2))
    IF (q < ZERO) u = -u

    50 x = -u - b
    y = v1 - v2
    z(1) = CMPLX(x, y,DP)
    z(2) = CMPLX(x,-y,DP)
    x =  u - b
    y = v1 + v2
    z(3) = CMPLX(x, y,DP)
    z(4) = CMPLX(x,-y,DP)
    RETURN

    !                THE RESOLVENT CUBIC HAS COMPLEX ROOTS

    60 t = DBLE(z(1))
    x = ZERO
    IF (t < ZERO) THEN
      GO TO 61
    ELSE IF (t == ZERO) THEN
      GO TO 70
    ELSE
      GO TO 62
    END IF
    61 h = ABS(DBLE(z(2))) + ABS(AIMAG(z(2)))
    IF (ABS(t) <= h) GO TO 70
    GO TO 80
    62 x = SQRT(t)
    IF (q > ZERO) x = -x

    70 w = SQRT(z(2))
      u = TWO*DBLE(w)
      v = TWO*ABS(AIMAG(w))
      t =  x - b
      x1 = t + u
      x2 = t - u
      IF (ABS(x1) <= ABS(x2)) GO TO 71
      t = x1
      x1 = x2
      x2 = t
    71 u = -x - b
      h = u*u + v*v
      IF (x1*x1 < 0.01D0*MIN(x2*x2,h)) x1 = e/(x2*h)
      z(1) = CMPLX(x1, ZERO,DP)
      z(2) = CMPLX(x2, ZERO,DP)
      z(3) = CMPLX(u, v,DP)
      z(4) = CMPLX(u,-v,DP)
      RETURN

    80 v = SQRT(ABS(t))
      z(1) = CMPLX(-b, v,DP)
      z(2) = CMPLX(-b,-v,DP)
      z(3) = z(1)
      z(4) = z(2)
      RETURN

    END Subroutine QuarticRoots

    !+
    SUBROUTINE SelectSort(a)
    ! ---------------------------------------------------------------------------
    ! PURPOSE - Reorder the elements of in increasing order.
      REAL(DP),INTENT(IN OUT),DIMENSION(:):: a

      INTEGER:: j
      INTEGER,DIMENSION(1):: k
    ! NOTE - This is a n**2 method. It should only be used for small arrays. <25
    !----------------------------------------------------------------------------
      DO j=1,SIZE(a)-1
        k=MINLOC(a(j:))
        IF (j /= k(1)) CALL Swap(a(k(1)),a(j))
      END DO
      RETURN
    END Subroutine SelectSort   ! -----------------------------------------------

    !+
    SUBROUTINE SolvePolynomial(quarticCoeff, cubicCoeff, quadraticCoeff, &
      linearCoeff, constantCoeff, code, root1,root2,root3,root4)
    ! ---------------------------------------------------------------------------
      REAL(DP),INTENT(IN):: quarticCoeff
      REAL(DP),INTENT(IN):: cubicCoeff, quadraticCoeff
      REAL(DP),INTENT(IN):: linearCoeff, constantCoeff
      INTEGER,INTENT(OUT):: code
      COMPLEX(DP),INTENT(OUT):: root1,root2,root3,root4
      REAL(DP),DIMENSION(5):: a
      COMPLEX(DP),DIMENSION(5):: z
    !----------------------------------------------------------------------------
      a(1)=constantCoeff
      a(2)=linearCoeff
      a(3)=quadraticCoeff
      a(4)=cubicCoeff
      a(5)=quarticCoeff

      IF (quarticCoeff /= ZERO) THEN
        CALL QuarticRoots(a,z)  
      ELSE IF (cubicCoeff /= ZERO) THEN
        CALL CubicRoots(a,z)
      ELSE IF (quadraticCoeff /= ZERO) THEN
        CALL QuadraticRoots(a,z)
      ELSE IF (linearCoeff /= ZERO) THEN
        z(1)=CMPLX(-constantCoeff/linearCoeff, 0, DP)
        outputCode=1
      ELSE
        outputCode=0    !  { no roots }
      END IF

      code=outputCode
      IF (outputCode > 0) root1=z(1)
      IF (outputCode > 1) root2=z(2)
      IF (outputCode > 23) root3=z(3)
      IF (outputCode > 99) root4=z(4)
      RETURN
    END Subroutine SolvePolynomial   ! ------------------------------------------

    !+
    SUBROUTINE SwapDouble(a,b)
    ! ---------------------------------------------------------------------------
    ! PURPOSE - Interchange the contents of a and b
      REAL(DP),INTENT(IN OUT):: a,b
      REAL(DP):: t
    !----------------------------------------------------------------------------
      t=b
      b=a
      a=t
      RETURN
    END Subroutine SwapDouble   ! -----------------------------------------------

    !+
    SUBROUTINE SwapSingle(a,b)
    ! ---------------------------------------------------------------------------
    ! PURPOSE - Interchange the contents of a and b
      REAL(SP),INTENT(IN OUT):: a,b
      REAL(SP):: t
    !----------------------------------------------------------------------------
      t=b
      b=a
      a=t
      RETURN
    END Subroutine SwapSingle   ! -----------------------------------------------


END Module PolynomialRoots   ! ==============================================



module librarymod

	use Weight_fn_mod

	real(kind(0.d0)),parameter :: pi=4.d0*atan(1.d0)
	real(kind(0.d0)),parameter :: const0 = 0.d0, const1 = 1.d0

	! use same name for integer or real(kind(0.d0)) args versions of imaxloc
	interface imaxloc
		module procedure imaxloc_int, imaxloc_dp
	end interface

	! use same name for the same logical operation; consider BLAS
	interface magnitude
		module procedure magnitude3, magnitudeN
	end interface

	!Surfacewighting function
	interface linearsurface_weight
		module procedure linearsurface_weight_Nmol, linearsurface_weight_1mol
	end interface

    private linearsurface_weight_Nmol, linearsurface_weight_1mol

	interface lagrange_poly_weight
		module procedure lagrange_poly_weight_Nmol, lagrange_poly_weight_1mol
	end interface

    private lagrange_poly_weight_Nmol, lagrange_poly_weight_1mol


	interface surface_array_to_nodes
		module procedure surface_array_to_nodes_incell, surface_array_to_nodes_multicell
	end interface

    private surface_array_to_nodes_incell, surface_array_to_nodes_multicell
	

	!Various Heavisides
	interface heaviside
		module procedure int_heaviside, int_array_heaviside, dp_heaviside, &
		                 dp_array_heaviside
	end interface

    private int_heaviside, int_array_heaviside, dp_heaviside, &
	        dp_array_heaviside

	! use same name for 1D and 3D plotting
	interface PYplot
		module procedure PYplot_1D, PYplot_3D
	end interface

	private PYplot_1D, PYplot_3D

	! use same name for 1D and 3D functions to return value on a shape_fn
	interface get_weight_fn
		module procedure get_weight_fn_scalar, get_weight_fn_vector
	end interface

	private get_weight_fn_scalar, get_weight_fn_vector



	!Assembly language interface to Heaviside function
    interface
        real(c_double) function heaviside_a1(x)
    		!DEC$ ATTRIBUTES STDCALL :: heaviside_a1
            use iso_c_binding, only: c_double
            real(c_double), VALUE :: x
        end function heaviside_a1

        real(c_double) function heaviside_a2(arg,const0,const1)
    		!DEC$ ATTRIBUTES STDCALL :: heaviside_a2
            use iso_c_binding, only: c_double
            real(c_double), VALUE :: arg,const0,const1
        end function heaviside_a2

    end interface

#if __INTEL_COMPILER > 1200

	type :: PDF

		integer								:: nbins
		integer,dimension(:),allocatable 	:: hist
		real(kind(0.d0))					:: minvalue,maxvalue,binsize

	contains

		procedure :: initialise  => PDF_initialise
		procedure :: update      => PDF_update
		procedure :: binvalues   => PDF_binvalues
		procedure :: normalise   => PDF_normalise
		procedure :: Hfunction   => PDF_Hfunction
		procedure :: moments     => PDF_moments
		procedure :: GNUplot     => PDF_GNUplot
		procedure :: PYplot      => PDF_PYplot
		procedure :: destroy     => PDF_destroy

	end type PDF

	interface PDF
		module procedure PDF_constructor
	end interface PDF

#endif
    ! Define a local library version of error_abort so library
    ! exist in isolation without the main code -- the downside
    ! is that MPI_ABORT is not used 

    interface error_abort
          module procedure error_abort_s, error_abort_si
    end interface error_abort
	
contains

subroutine error_abort_s(msg)
    implicit none

    character(len=*), intent(in), optional :: msg
   
	integer errcode,ierr

    if (present(msg)) then 
        write(*,*) msg
    endif

    stop 

end subroutine error_abort_s

subroutine error_abort_si(msg,i)
    implicit none

    character(len=*), intent(in) :: msg
    integer, intent(in) :: i

    integer errcode,ierr

    write(*,*) msg,i

    stop

end subroutine error_abort_si

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
	implicit none
	
	character*(*),intent(in) :: keyword             ! Input keyword	
	integer, intent(in) :: fileid                   ! File unit number
	logical, intent(in) :: required                 ! Flag check if input reqd
	logical, intent(out), optional :: input_present ! Optl return flag 

	character*(100) :: linestring                   ! First 100 chars in a line
	integer :: keyword_length                       ! Length of input keyword
	integer :: io                                   ! File status flag

	! Terminate if input_present flag is missing for non-required input 
	if (.not. required) then

		if(present(input_present)) then
	
			!Assume true until not found
			input_present = .true.

		else

			call error_abort( "ERROR IN LOCATE - If input not required, &
			                  &extra logical argument must be included to &
			                  &check if variable is specified in input file") 
			
		endif

	endif

	! Bulk of the locate routine
	keyword_length = len(keyword)
	rewind(fileid)
	do

		read (fileid,*,iostat=io) linestring
		if (linestring(1:keyword_length).eq.keyword) exit

		! If EOF reached, return or abort depending on required flag
		if ( io .ne. 0 ) then
			if ( .not. required ) then
				input_present = .false.
				return
			else
				call error_abort("ERROR IN LOCATE -  &
				                 &Required input "//trim(keyword)//"not &
				                 &found in input file. Stopping simulation.")
			endif
		end if

	end do	

end subroutine locate

!------------------------------------------------------------------------------
!Returns the magnitude of an 3 dimensional vector

function magnitude3(a)
	implicit none
	
	real(kind(0.d0))						:: magnitude3
	real(kind(0.d0)),dimension(3),intent(in):: a

    ! this should use BLAS library 
    ! magnitude = dnorm2(3,a,1)
 
	magnitude3 = sqrt((a(1)*a(1)+a(2)*a(2)+a(3)*a(3)))

end function magnitude3

!--------------------------------------------------------------------------------------
!Returns the magnitude of an N dimensional vector

function magnitudeN(a,n)
	implicit none
	
	integer,intent(in)			:: n
	real(kind(0.d0)),intent(in)	:: a(:)

	integer						:: i
	real(kind(0.d0))			:: magnitudeN


    ! simpler with a BLAS call
    ! magnituneN = dnorm2(n,a,1)

	magnitudeN = 0.d0

	do i=1,n
		magnitudeN = magnitudeN + a(i)*a(i)
	enddo

	magnitudeN = sqrt(magnitudeN)

end function  magnitudeN

!------------------------------------------------------------------------------
! Functions to switch between cartesian and cylindrical polar coords

function cpolariser(rin) result(rcp)
implicit none

	real(kind(0.d0)), intent(in)  :: rin(3)
	real(kind(0.d0))              :: rcp(3)

	real(kind(0.d0)) :: x,y,z

	x = rin(1)
	y = rin(2)
	z = rin(3)

	rcp(1) = sqrt(x*x + y*y)
	rcp(2) = modulo(atan2(y,x),2.d0*pi)
	rcp(3) = z	

end function cpolariser

function cpolarisev(vcart,theta) result(vpol)
	implicit none

	real(kind(0.d0)), intent(in) :: vcart(3), theta
	real(kind(0.d0)) :: vpol(3)

	real(kind(0.d0)) :: vx,vy,vz

	vx = vcart(1)
	vy = vcart(2)
	vz = vcart(3)

	vpol(1) =  vx*cos(theta) + vy*sin(theta) 
	vpol(2) = -vx*sin(theta) + vy*cos(theta)
	vpol(3) =  vz	

end function cpolarisev

function cpolariseT(Tcart,theta) result(Tpol)
	implicit none

	real(kind(0.d0)), intent(in) :: Tcart(3,3), theta
	real(kind(0.d0)) :: R(3,3), Tpol(3,3)

	! Rotation matrix for cylindrical polar. Fortran has column-major
	! ordering so we need to transpose what we "see" below.
	R = transpose(reshape((/cos(theta),sin(theta),0.d0,&
	                       -sin(theta),cos(theta),0.d0,&
	                        0.d0,      0.d0,      1.d0/), (/3,3/)))

	! Tpol = R * Tcart * transpose(R)
	Tpol = matmul(R, Tcart)
	Tpol = matmul(Tpol, transpose(R))
	
end function cpolariseT

function cartesianiser(rin) result(rca)
	implicit none
	
	real(kind(0.d0)), intent(in) :: rin(3)
	real(kind(0.d0))             :: rca(3)

	real(kind(0.d0)) :: r,t,z

	r = rin(1)
	t = rin(2)
	z = rin(3)

	rca(1) = r*cos(t) 
	rca(2) = r*sin(t)
	rca(3) = z

end function cartesianiser

function cartesianisev(vpol,theta) result(vcart)
implicit none

	real(kind(0.d0)), intent(in) :: vpol(3),theta 
	real(kind(0.d0)) :: vcart(3)

	real(kind(0.d0)) :: vr,vt,vz

	vr = vpol(1)
	vt = vpol(2)
	vz = vpol(3)

	vcart(1) =  vr*cos(theta) - vt*sin(theta)
	vcart(2) =  vr*sin(theta) + vt*cos(theta)
	vcart(3) =  vz	

end function cartesianisev

!------------------------------------------------------------------------------
! Functions to switch between cartesian and spherical coords

function sphereiser(rin) result(rcp)
implicit none

    real(kind(0.d0)), intent(in)  :: rin(3)
    real(kind(0.d0))              :: rcp(3)

    real(kind(0.d0)) :: x,y,z

    x = rin(1); y = rin(2); z = rin(3)

    rcp(1) = sqrt(x*x + y*y + z*z)
	if (rcp(1) .lt. 0.00000001d0) then 
		rcp = 0.d0
		return
	endif
    rcp(2) = modulo(atan2(y,x),2.d0*pi)
	if (z/rcp(1) .eq. -1.d0) then
		rcp(3) = pi
	else
    	rcp(3) = modulo(acos(z/rcp(1)),pi)
	end if

end function sphereiser

function sphereisev(vcart,theta,phi) result(vpol)
	implicit none

	real(kind(0.d0)), intent(in) :: vcart(3), theta, phi
	real(kind(0.d0)) :: R(3,3), vpol(3)

	! Rotation matrix for cylindrical polar. Fortran has column-major
	! ordering so we need to transpose what we "see" below.
	R = transpose(reshape((/cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi), &
	                       -sin(theta)         , cos(theta)	        , 0.d0,		&
	                       cos(theta)*cos(phi) , sin(theta)*cos(phi), -sin(phi)/), (/3,3/)))


	! vpol = R * Tcart * transpose(R)
	vpol = matmul(R, vcart)
	
end function sphereisev

! function sphereisev(rin,vin) result(vs)
! 	implicit none

! 	real(kind(0.d0)), intent(in) :: rin(3),vin(3)
! 	real(kind(0.d0)) :: vs(3), x2y2, x2y2z2

! 	x2y2 = rin(1)*vin(1) + rin(2)*vin(2) 
! 	x2y2z2 = x2y2 + rin(3)*vin(3)

! 	vs(1) = (rin(1)*vin(1) + rin(2)*vin(2) + rin(3)*vin(3))/sqrt(x2y2z2)
! 	vs(2) = (rin(1)*vin(2) - rin(2)*vin(1))/x2y2
! 	vs(3) = (rin(3)*(rin(1)*vin(1) + rin(2)*vin(2))-x2y2*vin(3))/ &
! 			x2y2z2*sqrt(x2y2)

! end function sphereisev




! function jacobian(theta,phi) result(vs)
! 	implicit none

! 	real(kind(0.d0)), intent(in) :: theta, phi
! 	real(kind(0.d0)) 			 :: vs(3,3)

! 	vs(1,:) = (/ sin(theta)*cos(phi), + sin(theta)*sin(phi), + cos(theta) /)
! 	vs(2,:) = (/ cos(theta)*sin(phi), + cos(theta)*sin(phi), - sin(theta) /)
! 	vs(3,:) = (/      -sin(phi),            + cos(phi),          0.d0     /)

! end function jacobian



! function sphere2cart(rin) result(rcp)
! implicit none

! 	real(kind(0.d0)), intent(in)  :: rin(3)
! 	real(kind(0.d0))              :: rcp(3)

! 	real(kind(0.d0)) :: x,y,z

! 	x = rin(1); y = rin(2);	z = rin(3)

! 	rcp(1) = x*sin(y)*cos(z)
! 	rcp(2) = x*sin(y)*sin(z)
! 	rcp(3) = x*cos(y)

! end function sphere2cart


!LINSPACE Linearly spaced vector.
!
!   LINSPACE(X1, X2, N) generates N points between X1 and X2.
!   For N = 1, LINSPACE returns X2.
!   Based on the MATLAB function

function linspace(d1, d2, n)
    implicit none

    integer,intent(in)              :: n
    real(kind(0.d0)),intent(in)     :: d1,d2
    real(kind(0.d0)),dimension(:),allocatable    :: linspace

    integer                                     :: i, n1

    n1 = n-1
    if (n1 .eq. 0) stop "Error in linspace -- Ensure n > 1"
    
    allocate(linspace(n))
    do i = 0,n1
        linspace(i+1) = d1 + i * (d2-d1)/(n1)
    enddo

end function linspace


!Get surface crossings on a specified control volume surface
! r1 -- position of molecules 1
! r2 -- position of molecules 2
! CV -- Array of top and botton location of cubic CV
!       1) xbinbot, 2) ybinbot, 3) zbinbot
!       4) xbintop, 5) ybintop, 6) zbintop
! onface  -- onface of surfaces of CV
!       1) xbinbot, 2) ybinbot, 3) zbinbot
!       4) xbintop, 5) ybintop, 6) zbintop
subroutine CV_surface_crossing(r1, r2, CV, onface) bind (C)
    use, intrinsic :: ISO_C_BINDING, only : C_INT32_T, C_DOUBLE

    real(C_DOUBLE),dimension(3), intent(in)     :: r1, r2
    real(C_DOUBLE), dimension(6), intent(in)    :: CV
    real(C_DOUBLE), dimension(6), intent(out)   :: onface

    integer(C_INT32_T)  		:: jxyz
    real(C_DOUBLE),dimension(3)	:: r12,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

	r12   = r1 - r2							!Molecule i trajectory between t-dt and t
	where (r12 .eq. 0.d0) r12 = 0.000001d0

	!Calculate the plane intersect of trajectory with plane aligned with cube surfaces
	Pxt=(/ 			               CV(4), 	     & 
			r1(2)+(r12(2)/r12(1))*(CV(4)-r1(1)), & 
			r1(3)+(r12(3)/r12(1))*(CV(4)-r1(1))  	/)
	Pxb=(/ 			               CV(1), 	     & 
			r1(2)+(r12(2)/r12(1))*(CV(1)-r1(1)), & 
			r1(3)+(r12(3)/r12(1))*(CV(1)-r1(1))  	/)
	Pyt=(/	r1(1)+(r12(1)/r12(2))*(CV(5)-r1(2)), & 
				                   CV(5),	     & 
			r1(3)+(r12(3)/r12(2))*(CV(5)-r1(2))  	/)
	Pyb=(/	r1(1)+(r12(1)/r12(2))*(CV(2)-r1(2)), &
				                   CV(2),	     & 
			r1(3)+(r12(3)/r12(2))*(CV(2)-r1(2))  	/)
	Pzt=(/	r1(1)+(r12(1)/r12(3))*(CV(6)-r1(3)), & 
			r1(2)+(r12(2)/r12(3))*(CV(6)-r1(3)), &
                                   CV(6) 			/)
	Pzb=(/	r1(1)+(r12(1)/r12(3))*(CV(3)-r1(3)), &
			r1(2)+(r12(2)/r12(3))*(CV(3)-r1(3)), & 
				                   CV(3) 			/)

	!Use Heavisides to determine if interaction actually on the cube's surface 
	onface(1) =0.5d0*(sign(1.d0,CV(1) - r2(1)) 	 & 
			        - sign(1.d0,CV(1) - r1(1)))* &
					 (heaviside(CV(5) - Pxb(2)) 	 &
			        - heaviside(CV(2) - Pxb(2)))* &
					 (heaviside(CV(6) - Pxb(3)) 	 &
			        - heaviside(CV(3) - Pxb(3)))
	onface(2) =0.5d0*(sign(1.d0,CV(2) - r2(2))   &
			        - sign(1.d0,CV(2) - r1(2)))* &
					 (heaviside(CV(4) - Pyb(1))   &
			        - heaviside(CV(1) - Pyb(1)))* &
					 (heaviside(CV(6) - Pyb(3))   &
			        - heaviside(CV(3) - Pyb(3)))
	onface(3) =0.5d0*(sign(1.d0,CV(3) - r2(3))   &
			        - sign(1.d0,CV(3) - r1(3)))* &
					 (heaviside(CV(4) - Pzb(1))   &
			        - heaviside(CV(1) - Pzb(1)))* &
					 (heaviside(CV(5) - Pzb(2))   &
			        - heaviside(CV(2) - Pzb(2)))

	onface(4) =0.5d0*(sign(1.d0,CV(4) - r2(1))   &
			        - sign(1.d0,CV(4) - r1(1)))* &
					 (heaviside(CV(5) - Pxt(2))   &
			        - heaviside(CV(2) - Pxt(2)))* &
					 (heaviside(CV(6) - Pxt(3))   &
			        - heaviside(CV(3) - Pxt(3)))
	onface(5) =0.5d0*(sign(1.d0,CV(5) - r2(2))   &
			        - sign(1.d0,CV(5) - r1(2)))* &
					 (heaviside(CV(4) - Pyt(1))   &
			        - heaviside(CV(1) - Pyt(1)))* &
					 (heaviside(CV(6) - Pyt(3))   &
			        - heaviside(CV(3) - Pyt(3)))
	onface(6) =0.5d0*(sign(1.d0,CV(6) - r2(3))   &
			        - sign(1.d0,CV(6) - r1(3)))* &
			 		 (heaviside(CV(4) - Pzt(1))   &
				    - heaviside(CV(1) - Pzt(1)))* &
			 		 (heaviside(CV(5) - Pzt(2))   &
			        - heaviside(CV(2) - Pzt(2)))

end subroutine CV_surface_crossing


!Get flux/force over a specified control volume surface
! r1 -- position of molecules 1
! r2 -- position of molecules 2
! CV -- Array of top and botton location of cubic CV
!       1) xbinbot, 2) ybinbot, 3) zbinbot
!       4) xbintop, 5) ybintop, 6) zbintop
! Quantity  -- Quantity of interest (mass, momentum, energy, etc)
! N         -- Number of components in Quantity
! fluxes  -- Additional fluxes over surfaces of CV
!       1) xbinbot, 2) ybinbot, 3) zbinbot
!       4) xbintop, 5) ybintop, 6) zbintop
subroutine CV_surface_flux(r1, r2, CV, Quantity, N, fluxes) bind (C)
    use, intrinsic :: ISO_C_BINDING, only : C_INT32_T, C_DOUBLE

    integer(C_INT32_T), intent(in)				:: N
    real(C_DOUBLE),dimension(3), intent(in)     :: r1, r2
    real(C_DOUBLE), dimension(6), intent(in)    :: CV
    real(C_DOUBLE), dimension(N), intent(in)    :: Quantity
    real(C_DOUBLE), dimension(N,6), intent(out) :: fluxes

    integer(C_INT32_T)  		:: jxyz
    real(C_DOUBLE)				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
    real(C_DOUBLE),dimension(3)	:: r12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
    real(C_DOUBLE),dimension(6)	:: onface

	r12   = r1 - r2							!Molecule i trajectory between t-dt and t
	where (r12 .eq. 0.d0) r12 = 0.000001d0

    call CV_surface_crossing(r1, r2, CV, onface)

	!Add quantity over face
	fluxes(:,1) = Quantity(:)*nint(dble(onface(1)))
    fluxes(:,2) = Quantity(:)*nint(dble(onface(2)))
    fluxes(:,3) = Quantity(:)*nint(dble(onface(3)))
    fluxes(:,4) = Quantity(:)*nint(dble(onface(4)))
    fluxes(:,5) = Quantity(:)*nint(dble(onface(5)))
    fluxes(:,6) = Quantity(:)*nint(dble(onface(6)))

end subroutine CV_surface_flux



!*****************************************************************************80
!
! Subroutine for calculating least squares straight line
!
!  A formula for a line of the form y = m x + c is sought, which
!  will minimize the root-mean-square error to N data points ( x), y );
!
!  Parameters:
!
!    Input, real(kind(0.d0)) X(:), Y(:), data points.
!    Output, real(kind(0.d0)) A, B, the slope and y-intercept

subroutine least_squares( x, y, m, c )
    implicit none

    real (kind(0.d0)),dimension(:),allocatable,intent(in)  :: x
    real (kind(0.d0)),dimension(:),allocatable,intent(in)  :: y

    real (kind(0.d0)), intent(out) :: m
    real (kind(0.d0)), intent(out) :: c

    integer           :: n
    real (kind(0.d0)) :: xbar
    real (kind(0.d0)) :: ybar
    real (kind(0.d0)) :: bot
    real (kind(0.d0)) :: top

    !Sanity check
    if (size(x,1) .ne. size(y,1)) then
        print*, size(x,1), size(y,1)
        stop "Error in least_squares - x must be same size as y"
    endif

    !Get number of values
    n = size(x,1)

    !  Special case.
    if (n .eq. 0) then
        m = 0.d0; c=0.d0
        return
    elseif ( n .eq. 1 ) then
        m = 0.d0; c = y(1)
        return
    end if

    !  Average X and Y.
    xbar = sum ( x(1:n) ) / real(n, kind(0.d0))
    ybar = sum ( y(1:n) ) / real(n, kind(0.d0))

    !  Compute a and b
    top = dot_product ( x(1:n) - xbar, y(1:n) - ybar )
    bot = dot_product ( x(1:n) - xbar, x(1:n) - xbar )

    if (bot .lt. 1e-8) then
        m = 0.d0
    else
        m = top / bot
    endif
    c = ybar - m * xbar

end subroutine least_squares


!------------------------------------------------------------------------------
!Subroutine for calculating least squares straight line with a uniform interval between x values

!subroutine least_squares(y, x_interval, npoints, lstsqrsinter, lstsqrsgrad)
!implicit none

!	integer		, intent(in)						:: npoints
!	real(kind(0.d0)), dimension(npoints), intent(in):: y
!	real(kind(0.d0)), intent(in)					:: x_interval
!	real(kind(0.d0)), intent(out)					:: lstsqrsgrad, lstsqrsinter

!	integer											:: n
!	real(kind(0.d0))								:: ndx, lstsqrsx,lstsqrsy,lstsqrsx2,lstsqrsxy

!	!Calculate molecular velocity using least squares to fit line
!	!and extrapolate down to point below overlap corresponding to continuum halo
!	lstsqrsy  = 0.d0
!	lstsqrsx  = 0.d0
!	lstsqrsx2 = 0.d0
!	lstsqrsxy = 0.d0

!	do n = 1, npoints
!        ndx = n*x_interval
!		lstsqrsx  = lstsqrsx  + ndx
!		lstsqrsx2 = lstsqrsx2 + ndx*ndx
!		lstsqrsy  = lstsqrsy  + y(n)
!		lstsqrsxy = lstsqrsxy + ndx*y(n)
!	enddo

!	!print*, lstsqrsx, lstsqrsx2, lstsqrsy, lstsqrsxy

!	!Calculate gradient
!	lstsqrsgrad = (npoints*lstsqrsxy-lstsqrsx*lstsqrsy) / &
!		          (npoints*lstsqrsx2-lstsqrsx*lstsqrsx)

!	!Calculate intercept
!	lstsqrsinter =(lstsqrsy*lstsqrsx2 - lstsqrsx*lstsqrsxy) / &
!		           (npoints*lstsqrsx2 - lstsqrsx*lstsqrsx)


!	print'(a,f10.5,a,f10.5)', 'y = ', lstsqrsgrad, ' x + ', lstsqrsinter

!end subroutine least_squares

!-------------------------------------------------------------------------------------
!Subrountine used to intergrate a function over a uniform grid using the trapizium rule

subroutine integrate_trap(y,x_interval,npoints,s)
implicit none

	integer		, intent(in)						:: npoints
	real(kind(0.d0)), intent(in)					:: x_interval
	real(kind(0.d0)), intent(out)					:: s
	real(kind(0.d0)), dimension(npoints), intent(in):: y

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
! Very simple bubble sorting algorithm -- should only be used for small data sets
! Orders smallest to biggest

subroutine bubble_sort(vec)
	implicit none

	real(kind(0.d0)), dimension(:),allocatable,intent(inout)	:: vec 

	integer 			:: bubble, vsize, j
	real(kind(0.d0))	:: temp

	vsize = size(vec) 

	do while (vsize .gt. 1)
		bubble = 0
		do j = 1, (vsize-1)
			if (vec(j) .gt. vec(j+1)) then
				temp = vec(j)
				vec(j) = vec(j+1)
				vec(j+1) = temp
				bubble = j
			endif 
		enddo
		vsize = bubble   
	enddo

end subroutine bubble_sort

! Reverse of the above but ordering biggest to smallest

subroutine bubble_sort_r(vec)
	implicit none

	real(kind(0.d0)), dimension(:),allocatable,intent(inout)	:: vec 

	integer 			:: bubble, vsize, j
	real(kind(0.d0))	:: temp

	vsize = size(vec) 

	do while (vsize .gt. 1)
		bubble = 0
		do j = 2, vsize
			if (vec(j) .gt. vec(j-1)) then
				temp = vec(j)
				vec(j) = vec(j-1)
				vec(j-1) = temp
				bubble = j
			endif 
		enddo
		vsize = bubble   
	enddo

end subroutine bubble_sort_r

!-------------------------------------------------------------------------------------
!Returns the heaviside function for input x -- interface at top

function int_heaviside(x)
	implicit none

	real(kind(0.d0))				:: int_heaviside
	integer	,intent(in)				:: x

	int_heaviside = nint(0.5*sign(1,x)+1)

end function

function int_array_heaviside(x)
	implicit none

	integer,dimension(:),intent(in)	    :: x
	real(kind(0.d0)),dimension(size(x))	:: int_array_heaviside

	int_array_heaviside = nint(0.5*sign(1,x(:))+1)

end function int_array_heaviside

function dp_heaviside(x)
	implicit none

	real(kind(0.d0))			:: dp_heaviside
	real(kind(0.d0)),intent(in)	:: x

	dp_heaviside = ceiling(sign(0.5d0,x))
	!heaviside = 0.5*sign(1.d0,x)+1

end function dp_heaviside

function dp_array_heaviside(x)
	implicit none

	real(kind(0.d0)),dimension(:),intent(in)	:: x
	real(kind(0.d0)),dimension(size(x))			:: dp_array_heaviside

	dp_array_heaviside = ceiling(sign(0.5d0,x(:)))

end function dp_array_heaviside


function sphereCV(r,radius)

	real(kind(0.d0)),intent(in)					:: radius
	real(kind(0.d0)),dimension(3),intent(in)	:: r

	real(kind(0.d0))							:: sphereCV
	real(kind(0.d0)),dimension(3)				:: rs

	!Convert to spherical coordinates
	rs = sphereiser(r)

	!Check vs spherical CV function
	sphereCV = heaviside_a1(radius-rs(1))

end function sphereCV

! An approximation ot the Heaviside function using tanh 
! k adjusts how aggressive it is
function heaviside_dp_approx(x,k)

	real(kind(0.d0)),dimension(:),intent(in)	:: x, k
	real(kind(0.d0)),dimension(size(x))			:: heaviside_dp_approx

	heaviside_dp_approx = 0.5d0*(1.d0 + tanh(k*x))

end function heaviside_dp_approx



!------------------------------------------------------------------------------
! Subroutine computes the intersection of a plane and a straight line
!Inputs: 
!       normal - normal vector of the Plane 
!       p - any point that belongs to the Plane 
!       ri- end point 1
!       rj- end point 2

subroutine plane_line_intersect(intersection,normal,p,ri,rj)
implicit none

	real(kind(0.d0)), dimension(3)				:: rij
	real(kind(0.d0)), dimension(3), intent(in)	:: ri, rj, normal, p
	real(kind(0.d0)), dimension(3), intent(out)	:: intersection

	rij = rj-ri

	!Use vector identity to calculate point of intersection.
	intersection = ri + (-dot_product(normal,(ri-p)) / dot_product(normal,rij))*rij

end subroutine plane_line_intersect

!--------------------------------------------------------------------------------------
!Swap column a and column b
subroutine swap(a,b)
	implicit none

	real(kind(0.d0)),dimension(:),intent(inout)	:: a, b
	real(kind(0.d0)),dimension(size(a,1))		:: temp

	if (size(a) .ne. size(b)) call error_abort("Array sizes different in swap")

	temp = a
	a = b
	b = temp
	
end subroutine swap

!--------------------------------------------------------------------------------------
!Calculate outer product of two vectors (uses Fortran intrinsic 'spread')
function outerprod(a,b)
	implicit none

	real(kind(0.d0)),dimension(:),intent(in)	:: a, b
	real(kind(0.d0)),dimension(size(a),size(b))	:: outerprod

	outerprod = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
	
end function outerprod


subroutine get_outerprod(a,b,outerprod)
	implicit none

	real(kind(0.d0)),dimension(:),intent(in)	            :: a, b
	real(kind(0.d0)),dimension(:,:),allocatable,intent(out)	:: outerprod

    allocate(outerprod(size(a),size(b)))
	outerprod = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
	
end subroutine get_outerprod

!--------------------------------------------------------------------------------------
!Calculate cross product of two 3D vectors 
function crossprod(a,b)
	implicit none

	real(kind(0.d0)),dimension(3),intent(in)	:: a, b
	real(kind(0.d0)),dimension(3)				:: crossprod

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
	real(kind(0.d0)),dimension(:),intent(in)	:: a

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
	implicit none

	real(kind(0.d0)), dimension(:,:), intent(inout)	:: A
	integer, dimension(:), intent(out)				:: indx
	real(kind(0.d0))								:: d
	real(kind(0.d0)), dimension(size(A,1))			:: vv
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
	real(kind(0.d0))				:: summ
	real(kind(0.d0)), dimension(:,:), intent(in)	:: A
	real(kind(0.d0)), dimension(:), intent(inout)	:: b

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
	character(len=*), intent(in)	:: filename

	!Check if unit number is used and assign unique number
	!do unit_no = 1,1000
	!	inquire(unit_no,opened=op)
	!	if (op .eqv. .false.) exit
	!enddo
	unit_no = get_new_fileunit()

	! simpler version using inquire
	inquire(file=filename,size=file_size)

	!Use Fstat to obtain size of file ## MAY BE COMPATABILITY ISSUE ##
	!open (unit=unit_no, file=filename)
	!call FStat(unit_no, SArray)
	!close(unit_no,status="keep")

	!file_size = SArray(8)

end subroutine get_file_size

!------------------------------------------------------------------------------
!Pure fortran subroutine to return an updated filename by appending
!the current timestep to that file
subroutine get_Timestep_FileName(timestep,basename,filename)
		implicit none

		integer,intent(in) 			:: timestep
		character(*),intent(in) 	:: basename
		character(*),intent(out)	:: filename

        if(timestep.le.9                         		) &
        write(filename,'(a,a7,i1)') trim(basename),'.000000',timestep
        if(timestep.ge.10      .and. timestep.le.99     ) &
        write(filename,'(a,a6,i2)') trim(basename),'.00000' ,timestep
        if(timestep.ge.100     .and. timestep.le.999    ) &
        write(filename,'(a,a5,i3)') trim(basename),'.0000'  ,timestep
        if(timestep.ge.1000    .and. timestep.le.9999   ) &
        write(filename,'(a,a4,i4)') trim(basename),'.000'   ,timestep
        if(timestep.ge.10000   .and. timestep.le.99999  ) &
        write(filename,'(a,a3,i5)') trim(basename),'.00'    ,timestep
        if(timestep.ge.100000  .and. timestep.le.999999 ) &
        write(filename,'(a,a2,i6)') trim(basename),'.0'     ,timestep
        if(timestep.ge.1000000 .and. timestep.le.9999999) &
        write(filename,'(a,a1,i7)') trim(basename),'.'      ,timestep

		!Remove any surplus blanks
		filename = trim(filename)

end subroutine get_Timestep_FileName
!------------------------------------------------------------------------------
!Pure fortran subroutine to return an updated filename by appending
!the current no to that file
subroutine get_FileName(no,basename,filename)
		implicit none

		integer,intent(in) 			:: no
		character(*),intent(in) 	:: basename
		character(*),intent(out)	:: filename

        if(no.le.9                         		) &
        write(filename,'(a,a10,i1)') trim(basename),'.000000000',no
        if(no.ge.10      .and. no.le.99     ) &
        write(filename,'(a,a9,i2)') trim(basename),'.00000000' ,no
        if(no.ge.100     .and. no.le.999    ) &
        write(filename,'(a,a8,i3)') trim(basename),'.0000000'  ,no
        if(no.ge.1000    .and. no.le.9999   ) &
        write(filename,'(a,a7,i4)') trim(basename),'.000000'   ,no
        if(no.ge.10000   .and. no.le.99999  ) &
        write(filename,'(a,a6,i5)') trim(basename),'.00000'    ,no
        if(no.ge.100000  .and. no.le.999999 ) &
        write(filename,'(a,a5,i6)') trim(basename),'.0000'     ,no
        if(no.ge.1000000 .and. no.le.9999999) &
        write(filename,'(a,a4,i7)') trim(basename),'.000'      ,no
        if(no.ge.10000000 .and. no.le.99999999) &
        write(filename,'(a,a3,i8)') trim(basename),'.00'      ,no
        if(no.ge.100000000 .and. no.le.999999999) &
        write(filename,'(a,a2,i9)') trim(basename),'.0'      ,no
        if(no.ge.1000000000 .and. no.le.9999999999) &
        write(filename,'(a,a1,i10)') trim(basename),'.'      ,no

		!Remove any surplus blanks
		filename = trim(filename)

end subroutine get_FileName


!------------------------------------------------------------------------------
! Returns the version number of the current code from the version control
! system -- in this case subversion

function get_version_number()
        
        logical        :: file_exists
        integer        :: unit_no, statusno
        character(30)  :: get_version_number

        ! External system call -- this is almost certain not to
        ! work in general (e.g. not intel and not linux)

        get_version_number = "DISABLED"

!		call system("svnversion > ./subversion_no_temp")
!		!statusno = system("svnversion > ./subversion_no_temp")
!
!        !Check if system call has worked and file exists
!        inquire(file='./subversion_no_temp', exist=file_exists)
!
!		!Read file and store unit number if it exists, otherwise N/A
!        !if (file_exists .and. statusno .eq. 0) then
!        if (file_exists) then
!
!            !Check if unit number is used and assign unique number
!            unit_no = get_new_fileunit()
!            open(unit=unit_no, file='./subversion_no_temp')
!            read(unit_no,*,IOSTAT=statusno) get_version_number
!
!			!If nothing is written in file, set to N/A
!			if (statusno .ne. 0)  then
!				close(unit_no,status='delete')
!				get_version_number = 'N/A'
!			else
!				close(unit_no,status='delete')
!			endif
!		else
!			get_version_number = 'N/A'
!		endif

end function get_version_number
 
!-----------------------------------------------------
! Build array of dimension ncells which maps to 3D Hilbert curve.
subroutine build_hilbert(ncells,Hcurve)
	implicit none

	integer,dimension(3),intent(in)					 :: ncells
	integer,dimension(:,:,:),allocatable,intent(out) :: Hcurve

	integer	:: icell,jcell,kcell,maxcells
	integer	:: hindex,i,n,m
	real(kind(0.d0))	:: shift
	real(kind(0.d0)),dimension(:),allocatable	:: x,y,z

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
		real(kind(0.d0)),dimension(:),allocatable,intent(out)   :: x,y,z
		real(kind(0.d0)),dimension(:),allocatable				:: xo,yo,zo


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


!-----------------------------------------------------
! Return random numbers from distributions

! Normal distribtion
function normal_dist()
	implicit none

	real(kind(0.d0))			  :: normal_dist
	real(kind(0.d0)),dimension(2) :: rand, randn

	!Use box-muller to get normally distributed random numbers
	call random_number(rand)
	randn(1) = sqrt(-2*log(rand(1)))*cos(2*pi*rand(2))
	randn(2) = sqrt(-2*log(rand(1)))*sin(2*pi*rand(2))

	normal_dist = randn(1)

end function normal_dist

! Number picked from Rayleigh velocity distribution 
! m*vi/(kb*T) * exp(-m(vi - u)**2/(2*kb*T))
function Rayleigh_vel(T,u)
	implicit none

	real(kind(0.d0))			:: T, u, Rayleigh_vel
	real(kind(0.d0))			:: rand
	real(kind(0.d0)),parameter 	:: kB = 1.d0	!Boltzmann's constant

	! Rayleigh distributed number about 0 generated and added to
	! mean u to give required velocity
	call random_number(rand)
	Rayleigh_vel = u + sqrt(-2.d0 * kB*T * log(rand))

end function Rayleigh_vel

! Number picked from Maxwell Boltzmann velocity distribution 
! N.B. One component of a molecules velocity vector 
! so this is just the normal distribution
! m/(kb*T) * exp(-m(vi - u)**2/(2*kb*T))
function Maxwell_Boltzmann_vel(T,u)
	implicit none

	real(kind(0.d0))			 :: T, u, Maxwell_Boltzmann_vel
	real(kind(0.d0)),dimension(2):: rand

	!Use box-muller to get normally distributed random numbers for 1D Maxwell Boltzmann
	call random_number(rand)
	Maxwell_Boltzmann_vel = u + sqrt(-2*T*log(rand(1)))*cos(2*pi*rand(2))

end function Maxwell_Boltzmann_vel

! Number picked from Maxwell Boltzmann speed distribution 
! N.B. This is the magnitude of all three velocity vectors
! [m/(kb*T)]^(3/2) * 4 * pi * (vi-u)^2 * exp(-m(vi - u)**2/(2*kb*T))
function Maxwell_Boltzmann_speed(T,u)
	implicit none

	real(kind(0.d0))			 :: T, u, Maxwell_Boltzmann_speed
	real(kind(0.d0)),dimension(3):: randn
	real(kind(0.d0)),dimension(4):: rand
	real(kind(0.d0)),parameter 	 :: kB = 1.d0	!Boltzmann's constant

	!Use box-muller to get normally distributed random numbers
	call random_number(rand)
	randn(1) = sqrt(-2*log(rand(1)))*cos(2*pi*rand(2))
	randn(2) = sqrt(-2*log(rand(1)))*sin(2*pi*rand(2))
	randn(3) = sqrt(-2*log(rand(3)))*cos(2*pi*rand(4))

	!Convert to Maxwell Boltzmann
	Maxwell_Boltzmann_speed = u + sqrt(kB*T*dot_product(randn,randn))

end function Maxwell_Boltzmann_speed

!Use Maxwell Boltzmann distribution to pick a 3 component velocity vector
! m/(kb*T) * exp(-m(vi - u)**2/(2*kb*T))
function Maxwell_Boltzmann_vel3(T,u)
	implicit none

	real(kind(0.d0))			 :: T
	real(kind(0.d0)),dimension(3):: Maxwell_Boltzmann_vel3, u 
	real(kind(0.d0)),dimension(4):: rand
	real(kind(0.d0)),parameter 	 :: kB = 1.d0	!Boltzmann's constant

	!Use box-muller to get normally distributed random numbers
	call random_number(rand)
	Maxwell_Boltzmann_vel3(1) = u(1) + sqrt(-2*kB*T*log(rand(1)))*cos(2*pi*rand(2))
	Maxwell_Boltzmann_vel3(2) = u(2) + sqrt(-2*kB*T*log(rand(1)))*sin(2*pi*rand(2))
	Maxwell_Boltzmann_vel3(3) = u(3) + sqrt(-2*kB*T*log(rand(3)))*cos(2*pi*rand(4))

end function Maxwell_Boltzmann_vel3

!--------------------------------------------------------------------------------------
! Split an integer into three factors minimising value of each

subroutine find3factors(n,nx,ny,nz)
	implicit none

	integer, intent(in)		:: n
	integer, intent(out)	:: nx,ny,nz

	integer :: nfactors,  minval1, minval2
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

	real(kind(0.d0)),dimension(:),intent(in):: buf
	integer, intent(in), optional			:: dplaces_in

	integer				:: n,dplaces
	real(kind(0.d0))	:: maxbuf,minbuf,order
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
	if (maxbuf .lt. 0.d0 .and. maxbuf .gt. -1.d0) n = n + 1 !For the case of -0.something
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


subroutine GNUplot(x,y)
	implicit none

	integer										:: i,unitno1,unitno2
	real(kind(0.d0)),dimension(:),intent(in)	:: x,y

	if (size(x,1) .ne. size(y,1)) call error_abort("Error in GNUplot - array sizes differ")

	unitno1 = get_new_fileunit()
	unitno2 = get_new_fileunit()

	!Generate temp file of results
	open(unitno1,FILE='./tempout')
	do i =1,size(x,1)
		write(unitno1,*) x(i),y(i)
	enddo

	!Write file to output this
	open(unitno2,FILE='./plot_tempout')
	write(unitno2,*) "plot 'tempout' u 1:2 w l"

	!Close all previous gnuplot windows and create new
	call system("killall gnuplot_x11")
	call system('gnuplot --persist plot_tempout')

	!Clean up all temp files
	close(unitno1,status='delete')
	close(unitno2,status='delete')

end subroutine GNUplot


subroutine PYplot_1D(x,y,routine,savefig)
	implicit none

	integer										:: i,unitno1
	real(kind(0.d0)),dimension(:),intent(in)	:: x,y
	character(*)								:: routine
	character(100)								:: callpython
	logical										:: savefig_
	logical,intent(in),optional					:: savefig

	if (.not. present(savefig)) then
		savefig_ = .false.
	else
		savefig_ = savefig
	endif

	if (size(x,1) .ne. size(y,1)) call error_abort("Error in GNUplot - array sizes differ")
	unitno1 = get_new_fileunit()

	!Generate temp file of results
	open(unitno1,FILE='./tempout')
	do i =1,size(x,1)
		write(unitno1,*) x(i),y(i)
	enddo

	!Close all previous gnuplot windows and create new
	if (savefig_) then
		call system("killall eog")
		write(callpython,*) "python2.7 ", routine, " ; eog hist.png &"
	else
		write(callpython,*) "python2.7 ", routine
	endif
	call system(trim(callpython))

	!Clean up all temp files
	close(unitno1,status='delete')

end subroutine PYplot_1D


subroutine PYplot_3D(x,y,routine,savefig)
	implicit none

	integer										:: i,unitno1
	real(kind(0.d0)),dimension(:,:),intent(in)	:: x,y
	character(*)								:: routine
	character(100)								:: callpython
	logical										:: savefig_
	logical,intent(in),optional					:: savefig

	if (.not. present(savefig)) then
		savefig_ = .false.
	else
		savefig_ = savefig
	endif


	if (size(x,1) .ne. size(y,1)) stop "Error in GNUplot - array sizes differ"

	unitno1 = get_new_fileunit()

	!Generate temp file of results
	open(unitno1,FILE='./tempout')
	do i =1,size(x,1)
		print'(a,6f10.5)', 'writing', x(i,1),y(i,1),x(i,2),y(i,2),x(i,3),y(i,3)
		write(unitno1,'(6f18.12)') x(i,1),y(i,1),x(i,2),y(i,2),x(i,3),y(i,3)
	enddo

	!Close all previous gnuplot windows and create new
	if (savefig_) then
		call system("killall eog")
		write(callpython,*) "python2.7 ", routine, " ; eog hist.png &"
	else
		write(callpython,*) "python2.7 ", routine
	endif
	call system(trim(callpython))

	!Clean up all temp files
	close(unitno1,status='delete')

end subroutine PYplot_3D



#if __INTEL_COMPILER > 1200
!Constructor for PDF object
function PDF_constructor(nbins_in,minvalue_in, maxvalue_in)
	implicit none

	! initialize objects
	type(PDF) 		:: PDF_constructor

	integer, intent(in)			 :: nbins_in
	real(kind(0.d0)), intent(in) :: minvalue_in,maxvalue_in

	PDF_constructor%nbins     = nbins_in
	PDF_constructor%minvalue  = minvalue_in
	PDF_constructor%maxvalue  = maxvalue_in
	PDF_constructor%binsize   = (maxvalue_in-minvalue_in)/nbins_in

	allocate(PDF_constructor%hist(nbins_in))
	PDF_constructor%hist = 0

end function PDF_constructor


!Initialise for PDF object -- subroutine alternative to constructor 
subroutine PDF_initialise(self, nbins_in,minvalue_in, maxvalue_in)
	implicit none

	integer, intent(in)			 :: nbins_in
	real(kind(0.d0)), intent(in) :: minvalue_in,maxvalue_in

	! initialize objects
	class(PDF) 		:: self

	self%nbins = nbins_in
	self%minvalue  = minvalue_in
	self%maxvalue  = maxvalue_in
	self%binsize = (maxvalue_in-minvalue_in)/nbins_in

	allocate(self%hist(nbins_in))
	self%hist = 0

end subroutine PDF_initialise

!Update cumulative PDF count with contents of array
subroutine PDF_update(self, array, checkrange)
	implicit none

	real(kind(0.d0)), intent(in), dimension(:) 	:: array
	logical,intent(in),optional					:: checkrange

	logical										:: check
	integer										:: i, bin

	! initialize objects
	class(PDF) 		:: self

	!Assign array values to bins
	do i = 1,size(array,1)
		bin = ceiling((array(i)-self%minvalue)/self%binsize)
		if (bin .gt. self%nbins	) bin = self%nbins
		if (bin .lt. 1			) bin = 1
		self%hist(bin) = self%hist(bin) + 1

	enddo

	!Check if range is large enough and extend if not
	if (present(checkrange)) then
		check = checkrange
	else
		check = .false.
	endif
	if (check .and. self%hist(self%nbins) .gt. self%hist(self%nbins-1)) then
		print*, "Warning, maximum specified for PDF is too small - extending",self%maxvalue, 'by', self%binsize
		self%maxvalue = self%maxvalue + self%binsize
		print*, "Current PDF will be reset"
		self%hist = 1
	endif

end subroutine PDF_update

!Return location of bins in space
function PDF_binvalues(self)
	implicit none

	integer										:: i
	real(kind(0.d0)),dimension(:),allocatable	:: PDF_binvalues 

	! initialize objects
	class(PDF) 		:: self

	allocate(PDF_binvalues(self%nbins))
	do i =1, self%nbins
		PDF_binvalues(i) = self%minvalue + (i-0.5d0)*self%binsize
	enddo

end function PDF_binvalues

! !Return normaised PDF function
! function PDF_normalise(self)
! 	implicit none

! 	real(kind(0.d0)),dimension(:),allocatable	:: PDF_normalise 

! 	! initialize shape objects
! 	class(PDF) 		:: self

! 	allocate(PDF_normalise(self%nbins))
! 	PDF_normalise  = dble(self%hist)/dble(sum(self%hist))

! end function PDF_normalise


!Return normaised PDF function using area under pdf...
function PDF_normalise(self)
	implicit none

	integer										:: n

	real(kind(0.d0)),dimension(:),allocatable	:: PDF_normalise
	real(kind(0.d0))							:: normalise_factor
	!initialize objects
	class(PDF) 		:: self

	allocate(PDF_normalise(self%nbins))
	PDF_normalise = dble(self%hist)

	!Get total area under curve
	normalise_factor = 0.d0
	call integrate_trap(PDF_normalise,self%binsize,self%nbins,normalise_factor)

	!Normalise PDF by area
	if (normalise_factor .gt. 0.00000001d0) then
		PDF_normalise  = PDF_normalise/normalise_factor
	else
		print*, 'Warning - PDF_normalise: area under PDF is zero'
		PDF_normalise  = 0.d0
	endif

end function PDF_normalise



!Return normaised PDF function
function PDF_moments(self,momentno)
	implicit none

	integer,intent(in)							:: momentno

	integer										:: n
	real(kind(0.d0))							:: PDF_moments , zeromoment
	real(kind(0.d0)),dimension(:),allocatable	:: binslocs, normalisedPDF

	! initialize shape objects
	class(PDF) 		:: self

	!Get bin values and normalised PDF
	allocate(binslocs(self%nbins),normalisedPDF(self%nbins))
	binslocs = self%binvalues()
	normalisedPDF = self%normalise()

	!Get zeroth moment (mean) to use as centre
	zeromoment = 0.d0
	do n=1,self%nbins
		zeromoment = zeromoment + normalisedPDF(n)*binslocs(n)
	enddo
	zeromoment = zeromoment / sum(normalisedPDF)

	!Get required mean
	PDF_moments = 0.d0
	if (momentno .eq. 1) then
		PDF_moments = zeromoment
	else
		do n=1,self%nbins
			PDF_moments = PDF_moments + normalisedPDF(n)*(binslocs(n)-zeromoment)**dble(momentno)
		enddo
		PDF_moments = PDF_moments / sum(normalisedPDF)
	endif

end function PDF_moments

!Calculate Boltzmann H function using discrete defintion as in
!Rapaport p37. N.B. velocity at middle or range is used for vn
function PDF_Hfunction(self)
	implicit none

	integer										:: n
	real(kind(0.d0))							:: PDF_Hfunction
	real(kind(0.d0)),dimension(:),allocatable	:: normalisedPDF 

	! initialize shape objects
	class(PDF) 		:: self

	PDF_Hfunction = 0.d0
	normalisedPDF = self%normalise()
	allocate(normalisedPDF(self%nbins))

	do n=1,size(normalisedPDF,1)
		if (normalisedPDF(n) .ne. 0) then
			PDF_Hfunction = PDF_Hfunction + normalisedPDF(n) & 
							* log(normalisedPDF(n)/(((n-0.5d0)*self%binsize)**2))
		endif
	enddo

end function PDF_Hfunction

!Plot histogram using gnuplot
subroutine PDF_GNUplot(self)
	implicit none

	integer										:: i
	real(kind(0.d0)),allocatable,dimension(:)	:: binloc,out

	! initialize shape objects
	class(PDF) 		:: self

	allocate(binloc(self%nbins),out(self%nbins))
	binloc = self%binvalues()
	out = self%normalise()

	call GNUplot(binloc,out)

end subroutine PDF_GNUplot

!Plot histogram using python
subroutine PDF_PYplot(self)
	implicit none

	integer										:: i
	real(kind(0.d0)),allocatable,dimension(:)	:: binloc,out
	character(11)								:: routine

	! initialize shape objects
	class(PDF) 		:: self

	allocate(binloc(self%nbins),out(self%nbins))
	binloc = self%binvalues()
	out = self%normalise()

	routine = "histplot.py"
	call PYplot(binloc,out,routine)

end subroutine PDF_PYplot

!Object destructor
subroutine PDF_destroy(self)
	implicit none

	! initialize shape objects
	class(PDF) 		:: self

	deallocate(self%hist)

end subroutine PDF_destroy

#endif

!Read input file from the DNS codes

subroutine read_DNS_velocity_files(filename,ngx,ngy,ngz,uc,vc,wc)
	implicit none

	integer,intent(in)		:: ngx, ngy, ngz
	character(*),intent(in)	:: filename
	real(kind(0.d0)), dimension(:,:,:),allocatable, intent(out)	:: uc, vc, wc

	integer					::  ifieldlength, file_size, unitno

	allocate(uc(ngz+1,ngx+2,ngy+1))
	allocate(vc(ngz+1,ngx+1,ngy+2))
	allocate(wc(ngz+2,ngx+1,ngy+1))

	!Get size of data unit
	inquire(iolength=ifieldlength) uc(1,1,1)
	ifieldlength = 4*ifieldlength*(  (ngz+1)*(ngx+2)*(ngy+1) &
									+(ngz+1)*(ngx+1)*(ngy+2) &
									+(ngz+2)*(ngx+1)*(ngy+1)  )

	!Get size of file and check
	call get_file_size(filename,file_size)
	if (ifieldlength .ne. file_size) then
		print'(5(a,i8))', 'File size = ', file_size, & 
						  ' ngx = ', ngx, ' ngy = ', ngy, ' ngz = ', ngz,  & 
						' Size based on 4*dble*ngx*ngy*ngz ', ifieldlength
		call error_abort("Error in read_DNS_velocity_files --  Filesize not equal to specified data to read in")
	endif

	!---------------------------------------------------------
	!	   Read in velocity field
	!---------------------------------------------------------
	unitno = get_new_fileunit()
	open (unitno,file=filename,form="unformatted",access="direct",recl=ifieldlength)
	read (unitno,rec=1) uc, vc, wc
	close(unitno)

end subroutine read_DNS_velocity_files



!Read input file from the FEA codes


subroutine read_FEA_output_files(filename, HLratio, nNodes, X, Z, & 
                                 H, Hx, Hxx, Hxxx, U, S2, &
                                 S2x, S12_or_S1, S12x_or_S1x)
	implicit none

	character(*),intent(in)	:: filename

    integer,intent(out),optional          :: nNodes
    real(kind(0.d0)),intent(out),optional :: HLratio
    real(kind(0.d0)),intent(out),optional, &
        allocatable,dimension(:)       :: X, Z, H, Hx, Hxx, Hxxx, U, & 
                                          S2, S2x, S12_or_S1, S12x_or_S1x

    integer                               :: N
    real(kind(0.d0))                      :: Time_, HLratio_
    real(kind(0.d0)),allocatable,dimension(:) :: X_, Z_, H_, Hx_, Hxx_, Hxxx_, U_, & 
                                                 S2_, S2x_, S12_or_S1_, S12x_or_S1x_

	integer					:: unitno, i
	character(30)        	:: char_temp

	unitno = get_new_fileunit()

    !'======================'
    !  'Time = ', Time
    !'======================'
    open(unitno, file=filename,action='read')
    read(unitno,*) char_temp 
    read(unitno,*) char_temp
    read(unitno,*) char_temp

    !'X','Z','H','Hx','Hxx','Hxxx','U', &
	!'S2','S2x','S12_or_S1','S12x_or_S1x'
    read(unitno,'(11A14)') char_temp, char_temp, char_temp, char_temp,& 
                           char_temp, char_temp, char_temp, char_temp,& 
                           char_temp, char_temp, char_temp          

    !Should read number of elements and eps from SPARAMETERS
    N = 601
    HLratio_ = 0.42 !0.02d0 !0.005d0

    !Loop through and read all the data
    allocate(X_(N), Z_(N), H_(N), Hx_(N), Hxx_(N), Hxxx_(N), & 
             U_(N), S2_(N), S2x_(N), S12_or_S1_(N), S12x_or_S1x_(N))

    do i = 1, N
        read(unitno,'(1X,11ES14.6)') X_(i), Z_(i), H_(i), Hx_(i), Hxx_(i), Hxxx_(i), & 
                                     U_(i), S2_(i), S2x_(i), S12_or_S1_(i), S12x_or_S1x_(i)
    enddo
    close(unitno)

    if (present(HLratio)) HLratio = HLratio_
    if (present(Nnodes)) Nnodes = N
    !if (present(X)) allocate(X, source = X_)
    !if (present(Z)) allocate(Z, source = Z_)
    !if (present(H)) allocate(H, source = H_)
    !if (present(Hx)) allocate(Hx, source  = Hx_)
    !if (present(Hxx)) allocate(Hxx, source  = Hxx_)
    !if (present(Hxxx)) allocate(Hxxx, source  = Hxxx_)
    !if (present(U)) allocate(U, source  = U_)
    !if (present(S2)) allocate(S2, source  = S2_)
    !if (present(S2x)) allocate(S2x, source  = S2x_)
    !if (present(S12_or_S1)) allocate(S12_or_S1, source  = S12_or_S1_)
    !if (present(S12x_or_S1x)) allocate(S12x_or_S1x, source  = S12x_or_S1x_)
    if (present(X)) then
        allocate(X(size(X_)))
        X = X_
    endif
    if (present(Z)) then
        allocate(Z(size(Z_)))
 		Z = Z_
    endif
    if (present(H)) then
        allocate(H(size(H_)))
		H = H_
    endif
    if (present(Hx)) then
        allocate(Hx(size(Hx_)))
		Hx = Hx_
    endif
    if (present(Hxx)) then
        allocate(Hxx(size(Hxx_)))
		Hxx = Hxx_
    endif
    if (present(Hxxx)) then
        allocate(Hxxx(size(Hxxx_)))
		Hxxx = Hxxx_
    endif
    if (present(U)) then
        allocate(U(size(U_)))
		U = U_
    endif
    if (present(S2)) then
        allocate(S2(size(S2_)))
		S2 = S2_
    endif
    if (present(S2x)) then
        allocate(S2x(size(S2x_)))
		S2x = S2x_
    endif
    if (present(S12_or_S1)) then
        allocate(S12_or_S1(size(S12_or_S1_)))
		S12_or_S1 = S12_or_S1_
    endif
    if (present(S12x_or_S1x)) then
        allocate(S12x_or_S1x(size(S12x_or_S1x_)))
		S12x_or_S1x = S12x_or_S1x_
    endif

end subroutine read_FEA_output_files



! Ouput spectral solution for velocity in unsteady Couette flow u(y,t). Input in form
! "couette_analytical_fn(time,
!                        Reynolds_number,
!                        wall_velocity,
!                        Domain_height, 
!                        required_number_of_u_points_in_y
!                        Sliding_Wall 0=top, 1=bottom, 2=both  )"
! N.B should be used in the syntactically salty form 
!  allocate(ucouette, source=couette_analytical_fn(iter*delta_t,Re,wallslidev(1),globaldomain(2),gnbins(2),2))


function couette_analytical_fn(t,Re,U_wall,L,npoints,slidingwall) result (u)
    implicit none

    integer,intent(in)             :: npoints,slidingwall
    real(kind(0.d0)),intent(in)    :: t,Re,U_wall,L
    real(kind(0.d0)),dimension(:),allocatable  :: u
    
	integer,parameter				:: top=0, bottom=1,both=2
    integer                        :: nmodes, n
    real(kind(0.d0))               :: k, uinitial, lambda

    real(kind(0.d0)),dimension(:),allocatable :: y

    nmodes = 5000
    k = 1.d0/Re
    !allocate(y,source=linspace(0.d0,L,npoints))
    allocate(y(npoints)); y = linspace(0.d0,L,npoints)
    uinitial = 0.d0  !Initial condition

    allocate(u(npoints)); u = 0.d0
    select case (slidingwall)
    case(top)
        !Uwall at top
        do n = 1,nmodes
            lambda = (n*pi/L)**2
            u(:)=u(:)+(uinitial*exp(-lambda*k*t) - (-1)**n *(2.d0/(n*pi))*U_wall*(1.d0-exp(-lambda*k*t)))*sin(n*pi*y/L) 
        enddo
        u(npoints) =  U_wall
    case(bottom)
        !Uwall at bottom
        do n = 1,nmodes
            lambda = (n*pi/L)**2
            u(:)=u(:)-(uinitial*exp(-lambda*k*t)     +      (2.d0/(n*pi))*U_wall*(1.d0-exp(-lambda*k*t)))*sin(n*pi*y/L)
        enddo
        u(1)       = -U_wall
    case(both)
        !Uwall at bottom and top
        do n = 1,nmodes
            lambda = (n*pi/L)**2
            u(:)=u(:)+(uinitial*exp(-lambda*k*t) - (-1)**n *(2.d0/(n*pi))*U_wall*(1.d0-exp(-lambda*k*t)))*sin(n*pi*y(:)/L) &
                     -(uinitial*exp(-lambda*k*t)     +      (2.d0/(n*pi))*U_wall*(1.d0-exp(-lambda*k*t)))*sin(n*pi*y(:)/L)
        enddo
        u(npoints) = U_wall
        u(1)       = -U_wall
    end select

end function couette_analytical_fn

! Calculate the analytical expression for stress in couette flow

! Input of the form:
! tau=couette_analytical_stress_fn(t,Re,U_wall,L,npoints)
! t      - time
! Re     - Reynolds number (Distance units as fraction of domain height L)
! U_wall - wall velocity
! L      - domain height
! npoints- number of points required (size of returned array)


function couette_analytical_stress_fn(t,Re,U_wall,L,npoints,slidingwall) result (tau)
    implicit none

    integer,intent(in)             :: npoints,slidingwall
    real(kind(0.d0)),intent(in)    :: t,Re,U_wall,L
    real(kind(0.d0)),dimension(:),allocatable  :: tau
    
    integer                        :: nmodes, n
    real(kind(0.d0))               :: k, uinitial, lambda
    real(kind(0.d0)),dimension(:),allocatable :: y

    nmodes = 5000
    k = 1.d0/Re
    !allocate(y,source=linspace(0.d0,L,npoints))
    allocate(y(npoints)); y = linspace(0.d0,L,npoints)
    allocate(tau(npoints)); tau = 0.d0

    ! - - - Calculate strain - - -
    select case(slidingwall)
    case(0)
        !Add zero wavenumber
        tau(:) = tau(:) + U_wall/L
        !Add time specific modes
        do n = 1,nmodes
            lambda = (n*pi/L)**2
            tau(:) = tau(:) + (2.d0/L*((-1.d0)**n*U_wall)*(exp(-lambda*k*t))) * cos(n*pi*y/L)
        enddo
    case(1)
        !Add zero wavenumber
        tau(:) = tau(:) + U_wall/L
        !Add time specific modes
        do n = 1,nmodes
            lambda = (n*pi/L)**2
            tau(:) = tau(:) + (-1.d0)**n*(2.d0/L*((-1.d0)**n*U_wall)*(exp(-lambda*k*t))) * cos(n*pi*y/L)
        enddo
    case(2)
        !Add zero wavenumber
        tau(:) = tau(:) + (U_wall)*2.d0/L
        !Add time specific modes
        do n = 1,nmodes
            lambda = (n*pi/L)**2
            tau(:) = tau(:) +            (2.d0/L*((-1.d0)**dble(n)*U_wall)*(exp(-lambda*k*t))) * cos(dble(n)*pi*y/L) &
                            + (-1.d0)**n*(2.d0/L*((-1.d0)**dble(n)*U_wall)*(exp(-lambda*k*t))) * cos(dble(n)*pi*y/L)
        enddo
    case default
        stop "Couette Stress Analytical Solution Needs either top or bottom wall to be sliding"
    end select

    ! Shear stress is strain times viscosity (NOTE THIS ASSUMES UNIT DENSITY)
    tau = k * tau

end function couette_analytical_stress_fn




! shape function

function Na(r,a)
    implicit none

	integer,intent(in)							:: a
    real(kind(0.d0)),dimension(3),intent(in)	:: r
	
	real(kind(0.d0))							:: Na
	real(kind(0.d0)),dimension(3,8), & 
		parameter	::  ra = reshape((/ -1d0, -1.d0, -1.d0, & 
                                         1d0, -1.d0, -1.d0, & 
                                        -1d0,  1.d0, -1.d0, & 
                                         1d0,  1.d0, -1.d0, & 
                                        -1d0, -1.d0,  1.d0, & 
                                        -1d0,  1.d0, -1.d0, & 
                                        -1d0,  1.d0,  1.d0, & 
                                         1d0,  1.d0,  1.d0    /), &
                                      (/3,8/))

	Na = 0.125d0 * (1.d0 + ra(1,a)*r(1)) & 
				 * (1.d0 + ra(2,a)*r(2)) &
				 * (1.d0 + ra(3,a)*r(3))

end function Na


! Take surface forces on a cube's surfaces and return the 
! values at the node of intersection. 
! = = Surface data assumed to be stored as 
! 1 = top x;  2 = top y;  3 = top z
! 4 = bottom x;  5 = bottom y;  6 = bottom z
! = = Node data assumed to be stored as
! 1 = [0,0,0]; 2=[1,0,0]; 3 =[0,1,0]; 4 =[1,1,0]
! 5 = [0,0,1]; 6=[1,0,1]; 7 =[0,1,1]; 8 =[1,1,1]

! Note -- this only takes the 3 components in the cell
! 		  that the node is currently being considered for
!		  i.e. node(:,1) for cell 2 is not equal to
!			   node(:,2) for cell 1 

function surface_array_to_nodes_incell(surface) result(nodes)
    implicit none

    real(kind(0.d0)),dimension(3,6),intent(in):: surface
	real(kind(0.d0)),dimension(3,8)			  :: nodes

	nodes(:,1) = surface(:,4)+surface(:,5)+surface(:,6)
	nodes(:,2) = surface(:,1)+surface(:,5)+surface(:,6)
	nodes(:,3) = surface(:,2)+surface(:,4)+surface(:,6)
	nodes(:,4) = surface(:,1)+surface(:,2)+surface(:,6)
	nodes(:,5) = surface(:,3)+surface(:,4)+surface(:,5)
	nodes(:,6) = surface(:,1)+surface(:,3)+surface(:,5)
	nodes(:,7) = surface(:,2)+surface(:,3)+surface(:,4)
	nodes(:,8) = surface(:,1)+surface(:,2)+surface(:,3)

end function surface_array_to_nodes_incell


! Take surface forces on a cube's surfaces and return the 
! values at the node of intersection. 
! = = Surface data assumed to be stored as 
! 1 = top x;  2 = top y;  3 = top z
! 4 = bottom x;  5 = bottom y;  6 = bottom z
! = = Node data assumed to be stored as
! 1 = [0,0,0]; 2=[1,0,0]; 3 =[0,1,0]; 4 =[1,1,0]
! 5 = [0,0,1]; 6=[1,0,1]; 7 =[0,1,1]; 8 =[1,1,1]

! Note -- this only takes the 3 components in the cell
! 		  that the node is currently being considered for
!		  i.e. node(:,1) for cell 2 is not equal to
!			   node(:,2) for cell 1 


function surface_array_to_nodes_multicell(surface) result(nodes)
    implicit none

    real(kind(0.d0)),dimension(3,3,3,3,6),intent(in):: surface

    integer :: i,j,k
	real(kind(0.d0)),dimension(3,8)			  :: nodes


    !Each node is where 8 cubes join, each with three relevant faces. As
    !each face must touch another, the total no. of unique faces is 12. This
    !is expressed here as 3 on the central cube, 3 on it's diagonal and 2 on
    !each cube attached by a surface to the central cube

	nodes(:,1) = get_node(surface(1:2,1:2,1:2,:,:))
	nodes(:,2) = get_node(surface(2:3,1:2,1:2,:,:))
	nodes(:,3) = get_node(surface(1:2,2:3,1:2,:,:))
	nodes(:,4) = get_node(surface(2:3,2:3,1:2,:,:))
	nodes(:,5) = get_node(surface(1:2,1:2,2:3,:,:))
	nodes(:,6) = get_node(surface(2:3,1:2,2:3,:,:))
	nodes(:,7) = get_node(surface(1:2,2:3,2:3,:,:))
	nodes(:,8) = get_node(surface(2:3,2:3,2:3,:,:))

end function surface_array_to_nodes_multicell


function get_node(surface) result(node)
    implicit none

    real(kind(0.d0)),dimension(2,2,2,3,6),intent(in):: surface
	real(kind(0.d0)),dimension(3)			  :: node

    real(kind(0.d0)),parameter                :: divider = 1.d0
	node(:) =(+surface(1,1,1,:,1)+surface(1,1,1,:,2)+surface(1,1,1,:,3) &
			  +surface(2,2,2,:,4)+surface(2,2,2,:,5)+surface(2,2,2,:,6) &
              +surface(1,2,2,:,5)+surface(1,2,2,:,6) &
              +surface(2,2,1,:,4)+surface(2,2,1,:,5) & 
              +surface(1,1,2,:,1)+surface(2,1,1,:,3) )/divider

end function get_node

!===================================================
! Functions to return a polynomial weighting function
! at a given molecular position, based on the
! array of surface fluxes passed into the function

!Wrapper for single molecule
function linearsurface_weight_1mol(array,r_in,binsize,domain,shiftmean,meanvalue) result(weight)

    real(kind(0.d0)),dimension(3),intent(in)	:: domain,binsize
    real(kind(0.d0)),dimension(:),intent(in)	:: r_in
    real(kind(0.d0)),dimension(:,:,:,:,:),allocatable,intent(in)    :: array
	integer,intent(in),optional	:: shiftmean
	real(kind(0.d0)),dimension(:,:,:,:),allocatable,intent(in),optional	:: meanvalue

    real(kind(0.d0)),dimension(3,1)	    :: buf
    real(kind(0.d0)),dimension(:,:),allocatable	    :: weight
    
    buf(:,1) =  r_in(:)
	if (present(shiftmean)) then
		if (.not. present(meanvalue)) then
			weight = linearsurface_weight_Nmol(array,buf,binsize,domain,shiftmean)
			!stop "Error in linearsurface_weight_Nmol -- shiftmean requested and meanvalue not specified"
		else
			weight = linearsurface_weight_Nmol(array,buf,binsize,domain,shiftmean,meanvalue)
		endif
	else
		weight = linearsurface_weight_Nmol(array,buf,binsize,domain)
	endif


end function linearsurface_weight_1mol


! Use linear interpolation between the surfaces of the cube 
!(I think this is equivalent to nodes and lagrange interpolates below)

function linearsurface_weight_Nmol(array,r_in,binsize,domain,shiftmean,meanvalue) result(weight)
    implicit none

    real(kind(0.d0)),dimension(3),intent(in)	:: domain,binsize
    real(kind(0.d0)),dimension(:,:),intent(in)	:: r_in
    real(kind(0.d0)),dimension(:,:,:,:,:),allocatable,intent(in)    :: array
	integer,intent(in),optional	:: shiftmean
	real(kind(0.d0)),dimension(:,:,:,:),allocatable,intent(in),optional	:: meanvalue

    integer,dimension(:,:,:),allocatable :: nperbin
    real(kind(0.d0)),dimension(:,:),allocatable	    :: weight
    real(kind(0.d0)),dimension(:,:,:),allocatable 	:: sqr_term
    real(kind(0.d0)),dimension(:,:,:,:),allocatable :: wsum_bin, meanvalue_

	integer							:: npoints,n
    integer,dimension(3)            :: bin, order_nd, nbins
    real(kind(0.d0))				:: fxfyfz
    real(kind(0.d0)),dimension(3)   :: r_in_, Na
    real(kind(0.d0)),dimension(:,:),allocatable :: grid, rhat
    real(kind(0.d0)),dimension(3,6)	:: surfaces

	nbins = nint(domain/binsize)

	if (present(shiftmean)) then
		allocate(meanvalue_(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
		if (.not. present(meanvalue)) then
			meanvalue_ = 0.d0
		else
			meanvalue_ = meanvalue
		endif
		select case(shiftmean)
		case(0)
			!Do nothing, shiftmean not requested
		case(1)
			allocate(wsum_bin(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
			allocate(nperbin(nbins(1)+2,nbins(2)+2,nbins(3)+2))
			wsum_bin = 0.d0; nperbin = 0
		case(2)
			allocate(wsum_bin(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
			allocate(nperbin(nbins(1)+2,nbins(2)+2,nbins(3)+2))
			allocate(sqr_term(nbins(1)+2,nbins(2)+2,nbins(3)+2))
			sqr_term = 0.d0; wsum_bin = 0.d0; nperbin = 0
		case default
			stop "Error in linearsurface_weight_Nmol -- incorrect shiftmean value"
		end select
	endif

	!Setup polynomial and other functions
	allocate(  rhat(size(r_in,1),size(r_in,2)))
	allocate(weight(size(r_in,1),size(r_in,2)))

	do n =1,size(r_in,2)
    	! Shift to all positive, get bin
		! and map to bin local coordinate system (0 to +1)
		r_in_(:) = r_in(:,n)+0.5d0*domain(:)
		bin(:) = ceiling((r_in_)/binsize(:))+1
		rhat(:,n) = (r_in_(:)/binsize(:) - dble(bin(:)-2))

		if (any(rhat(:,n)-epsilon(rhat(:,n)) .gt. 1.d0) .or. & 
            any(rhat(:,n)+epsilon(rhat(:,n)) .lt. 0.d0)) then
		    stop "Error in lagrange_poly_weight_Nmol --rhat must satisfy 0 < rhat < 1"
		endif

		!Setup polynomial
        surfaces(:,:) = array(bin(1),bin(2),bin(3),:,:)
        Na(:) = rhat(:,n)
        weight(:,n) = ((surfaces(:,4)*(1.d0-Na(1)) + surfaces(:,1)*Na(1))+ &
                       (surfaces(:,5)*(1.d0-Na(2)) + surfaces(:,2)*Na(2))+ &
                       (surfaces(:,6)*(1.d0-Na(3)) + surfaces(:,3)*Na(3)))

		if (present(shiftmean)) then
			select case(shiftmean)
			case(0)
				!Do nothing, shiftmean not requested
			case(1)
				!Get sum of weightings to subtract so mean is zero by adjusting constant value
				wsum_bin(bin(1),bin(2),bin(3),:) = wsum_bin(bin(1),bin(2),bin(3),:) + weight(:,n)
				nperbin(bin(1),bin(2),bin(3)) 	 = nperbin(bin(1),bin(2),bin(3)) + 1
			case(2)
				! Zero mean by adding a 2nd order term preserving the B.C. and adjusting the 
				! 3D parabolicness until the sum is correct
				wsum_bin(bin(1),bin(2),bin(3),:) = wsum_bin(bin(1),bin(2),bin(3),:) + weight(:,n)
				nperbin(bin(1),bin(2),bin(3)) 	 = nperbin(bin(1),bin(2),bin(3)) + 1
				fxfyfz = ((Na(1)-0.5d0)**2-0.25d0) &
						*((Na(2)-0.5d0)**2-0.25d0) &
						*((Na(3)-0.5d0)**2-0.25d0)
				sqr_term(bin(1),bin(2),bin(3)) = sqr_term(bin(1),bin(2),bin(3)) + fxfyfz
			end select
		endif

	enddo

	if (present(shiftmean)) then
		select case(shiftmean)
		case(0)
			!Do nothing, shiftmean not requested
		case(1)
			!Zero mean by adjusting constant value
			do n =1,size(r_in,2)
				r_in_(:) = r_in(:,n)+0.5d0*domain(:)
				bin(:) = ceiling((r_in_)/binsize(:))+1
				weight(:,n) = weight(:,n) - (wsum_bin(bin(1),bin(2),bin(3),:) & 
							-meanvalue_(bin(1),bin(2),bin(3),:))/nperbin(bin(1),bin(2),bin(3))

			enddo
		case(2)
			! Zero mean by adding a 2nd order term preserving the B.C. and adjusting the 
			! parabolicness until the sum is correct
			do n =1,size(r_in,2)
				r_in_(:) = r_in(:,n)+0.5d0*domain(:)
				bin(:) = ceiling((r_in_)/binsize(:))+1
				Na(:) = rhat(:,n)
				fxfyfz = ((Na(1)-0.5d0)**2-0.25d0) &
						*((Na(2)-0.5d0)**2-0.25d0) &
						*((Na(3)-0.5d0)**2-0.25d0)
				weight(:,n) = weight(:,n) - fxfyfz * (wsum_bin(bin(1),bin(2),bin(3),:) & 
							-meanvalue_(bin(1),bin(2),bin(3),:))/sqr_term(bin(1),bin(2),bin(3))
			enddo

		end select
	endif

end function linearsurface_weight_Nmol


! Derive a weight function of the form f(x) = ax^2 + bx + c for 0 ≤ x ≤ 1.0  
! subject to the boundary conditions :
! INPUTS:
!	tBC -- top boundary condition, f(x=1) = tBC
!	bBC -- bottom boundary conditions, f(x=0) = bBC
!	sBC -- sum boundary condition, the sum of f(x) at a supplied selection of 
!		   x values must satisfy \sum f(x_n) = sBC
!	x_n -- points to evaluate that must sum to sBC
! OUTPUTS:
!	fx_n -- An array of values on the weight function at x_n, f(x_n) 


function get_weight_fn_scalar(tBC,bBC,sBC,x_n) result(fx_n)
	implicit none

	real(kind(0.d0)),intent(in)								:: tBC,bBC,sBC
	real(kind(0.d0)),dimension(:),allocatable,intent(in)	:: x_n
	real(kind(0.d0)),dimension(:),allocatable				:: fx_n

	integer													:: n, Np
	real(kind(0.d0))										:: sumx_n, sumx_n2
	real(kind(0.d0))										:: diffsumx
	real(kind(0.d0))										:: a,b,c

	!Get sum of x_n and sum of x_n^2
	sumx_n = 0.d0; sumx_n2 = 0.d0
	Np = size(x_n,1)	!Number of points
	do n = 1,Np
		sumx_n  = sumx_n  + x_n(n) 
		sumx_n2 = sumx_n2 + x_n(n)**2.d0
	enddo
	diffsumx = sumx_n - sumx_n2

	!Set constant to bottom boundary condition
	c = bBC

	!Expression for b is given by
	b = (sBC - Np*c - (tBC-c)*sumx_n2)/diffsumx

	!Expression for a is given by
	a = tBC - b - c

	!Get values of f(x_n)
	allocate(fx_n(Np))
	do n = 1,Np
		fx_n(n) = a*x_n(n)**2.d0 + b*x_n(n) + c
	enddo

end function get_weight_fn_scalar


! Derive a weight function of the form f(x) = ax^2 + bx + c for 0 ≤ x ≤ 1.0  
! subject to the boundary conditions :
!	tBC -- top boundary condition, f(x=1) = tBC
!	bBC -- bottom boundary conditions, f(x=0) = bBC
!	sBC -- sum boundary condition, the sum of f(x) at a supplied selection of 
!		   x values must satisfy \sum f(x_n) = sBC
!	x_n -- points that must sum to sBC

function get_weight_fn_vector(tBC,bBC,sBC,x_n) result(fx_n)
	implicit none

	real(kind(0.d0)),dimension(3),intent(in)				:: tBC,bBC,sBC
	real(kind(0.d0)),dimension(:),allocatable,intent(in)	:: x_n
	real(kind(0.d0)),dimension(:,:),allocatable				:: fx_n

	integer	:: n, Np
	real(kind(0.d0))	:: sumx_n, sumx_n2,diffsumx
	real(kind(0.d0)),dimension(3)	:: a,b,c

	!Get sum of x_n and sum of x_n^2
	sumx_n = 0.d0; sumx_n2 = 0.d0
	Np = size(x_n,1)	!Number of points
	do n = 1,Np
		sumx_n  = sumx_n  + x_n(n) 
		sumx_n2 = sumx_n2 + x_n(n)**2.d0
	enddo
	diffsumx = sumx_n - sumx_n2

	!Set constant to bottom boundary condition
	c = bBC

	!Expression for b is given by
	b = (sBC(:) - Np*c(:) - (tBC(:)-c(:))*sumx_n2)/diffsumx

	!Expression for a is given by
	a = tBC - b - c

	!Get values of f(x_n)
	allocate(fx_n(Np,3))
	do n = 1,Np
		fx_n(n,:) = a(:)*x_n(n)**2.d0 + b(:)*x_n(n) + c(:)
	enddo

end function get_weight_fn_vector


! Use lagrange interpolation between the nodes of the cube 

!Wrapper for single molecule
function lagrange_poly_weight_1mol(array, r_in, binsize, domain, & 
								   order, shiftmean, meanvalue) result(weight)

    real(kind(0.d0)),dimension(3),intent(in)	:: domain,binsize
    real(kind(0.d0)),dimension(:),intent(in)	:: r_in
    real(kind(0.d0)),dimension(:,:,:,:,:),allocatable,intent(in)    :: array
	integer,intent(in),optional	:: shiftmean, order
	real(kind(0.d0)),dimension(:,:,:,:),allocatable,intent(in),optional	:: meanvalue

	integer								:: order_
    real(kind(0.d0)),dimension(3,1)	    :: buf
    real(kind(0.d0)),dimension(:,:),allocatable	    :: weight
    
    buf(:,1) =  r_in(:)
	if (.not. present(order)) then
		order_ = 2	!Default linear (minimum order)
	endif
	if (present(shiftmean)) then
		if (.not. present(meanvalue)) then
			weight = lagrange_poly_weight_Nmol(array,buf,binsize,domain,order=order_,shiftmean=shiftmean)
			!stop "Error in lagrange_poly_weight_Nmol -- shiftmean requested and meanvalue not specified"
		else
			weight = lagrange_poly_weight_Nmol(array,buf,binsize,domain, & 
											   order=order_,shiftmean=shiftmean, & 
											   meanvalue=meanvalue)
		endif
	else
		weight = lagrange_poly_weight_Nmol(array,buf,binsize,domain,order=order_)
	endif


end function lagrange_poly_weight_1mol


function lagrange_poly_weight_Nmol(array, r_in, binsize, domain, & 
								   order, node_averages, shiftmean, meanvalue) result(weight)
    implicit none

    real(kind(0.d0)),dimension(3),intent(in)	:: domain,binsize
    real(kind(0.d0)),dimension(:,:),intent(in)	:: r_in
    real(kind(0.d0)),dimension(:,:,:,:,:),allocatable,intent(in)    :: array
	integer,intent(in),optional	:: shiftmean, order, node_averages
	real(kind(0.d0)),dimension(:,:,:,:),allocatable,intent(in),optional	:: meanvalue

	integer							:: npoints,n, node_ave
    integer,dimension(3)            :: bin, order_nd, nbins
    real(kind(0.d0))				:: fxfyfz
    real(kind(0.d0)),dimension(3)   :: r_in_, Na, lowerlim,upperlim
    real(kind(0.d0)),dimension(:,:),allocatable :: grid, rhat
    real(kind(0.d0)),dimension(3,6)	:: surfaces
    real(kind(0.d0)),dimension(3,8)	:: nodes,test
    real(kind(0.d0)),dimension(3,27):: nodes_order3

    integer,dimension(:,:,:),allocatable :: nperbin
    real(kind(0.d0)),dimension(:,:),allocatable	    :: weight
    real(kind(0.d0)),dimension(:,:,:),allocatable 	:: sqr_term
    real(kind(0.d0)),dimension(:,:,:,:),allocatable :: wsum_bin,meanvalue_

	!Setup polynomial and other functions
	allocate(  rhat(size(r_in,1),size(r_in,2)))
	allocate(weight(size(r_in,1),size(r_in,2)))
	order_nd(:) = order	!Assumed same order in all dimensions
	lowerlim(:) = (/ -1.d0,-1.d0,-1.d0 /)
	upperlim(:) = (/  1.d0, 1.d0, 1.d0 /)

    !Number of surrounding cells to average for node
    ! 0 -- in cell nodes only (nodes not equal between adjacent cells)
    ! 1 -- surrounding cells (nodes equal between cells but surfaces may not be)
    if (.not. present(node_averages)) then
        node_ave = 1
    else
        node_ave = node_averages
        if (node_ave .gt. 1) then
            stop "Error in lagrange_poly_weight -- node_averages should be zero or one"
        endif
    endif

	nbins = nint(domain/binsize)
	if (present(shiftmean)) then
		allocate(meanvalue_(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
		if (.not. present(meanvalue)) then
			meanvalue_ = 0.d0
		else
			meanvalue_ = meanvalue
		endif
		select case(shiftmean)
		case(0)
			allocate(wsum_bin(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)) !TEMP FOR DEBUG
			!Do nothing, shiftmean not requested
		case(1)
			allocate(wsum_bin(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
			allocate(nperbin(nbins(1)+2,nbins(2)+2,nbins(3)+2))
			wsum_bin = 0.d0; nperbin = 0
		case(2)
            if (node_ave .eq. 1) then
                print*, "Error in lagrange_poly_weight -- node_average should",&
                "be 0 (in cell, i.e. Finite Volume) with polynomial shiftmean"
                stop
            end if
			allocate(wsum_bin(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
			allocate(nperbin(nbins(1)+2,nbins(2)+2,nbins(3)+2))
			allocate(sqr_term(nbins(1)+2,nbins(2)+2,nbins(3)+2))
			sqr_term = 0.d0; wsum_bin = 0.d0; nperbin = 0
		case default
			print*, "Error in linearsurface_weight -- incorrect shiftmean value"
            stop
		end select
	endif

	do n =1,size(r_in,2)
    	! Shift to all positive, get bin
		! and map to bin local coordinate system (-1 to +1)
		r_in_(:) = r_in(:,n)+0.5d0*domain(:)
		bin(:) = ceiling((r_in_)/binsize(:))+1
		rhat(:,n) = 2.d0*(r_in_(:)/binsize(:) - dble(bin(:)-2))-1.d0

		!Setup polynomial
        if (order .eq. 0) then
            stop "Error in linearsurface_weight -- order must be greater than zero"
        elseif (order .eq. 1) then
            weight = 1.d0
            !stop "Error in lagrange_poly_weight -- First order weighting is a constant "
        elseif (order .eq. 2) then
            if (node_ave .eq. 0) then
                !I think a factor of 2 is needed here as each surface only counted once
           		nodes(:,:) = 2.d0*surface_array_to_nodes(array(bin(1),bin(2),bin(3),:,:))
            elseif (node_ave .eq. 1) then
                !I think a factor of 2 is needed here as each surface only counted once
           		nodes(:,:) = 2.d0*surface_array_to_nodes(array(bin(1)-node_ave:bin(1)+node_ave, & 
	            										       bin(2)-node_ave:bin(2)+node_ave, & 
	            										       bin(3)-node_ave:bin(3)+node_ave,:,:))
            endif
		    call lagrange_interp_nd_size ( 3, order_nd, npoints )
		    call lagrange_interp_nd_value( 3, order_nd, lowerlim, upperlim, npoints, nodes(1,:), 1, rhat(:,n), weight(1,n))
		    call lagrange_interp_nd_value( 3, order_nd, lowerlim, upperlim, npoints, nodes(2,:), 1, rhat(:,n), weight(2,n))
		    call lagrange_interp_nd_value( 3, order_nd, lowerlim, upperlim, npoints, nodes(3,:), 1, rhat(:,n), weight(3,n))
        elseif (order .ge. 3) then
            stop "Error in lagrange_poly_weight -- Higher orders than 2 are not currently possible"
        endif

		if (present(shiftmean)) then
			select case(shiftmean)
			case(0)
				wsum_bin(bin(1),bin(2),bin(3),:) = wsum_bin(bin(1),bin(2),bin(3),:) + weight(:,n) !TEMP FOR DEBUG
				!Do nothing, shiftmean not requested
			case(1)
				!Get sum of weightings to subtract so mean is zero by adjusting constant value
				wsum_bin(bin(1),bin(2),bin(3),:) = wsum_bin(bin(1),bin(2),bin(3),:) + weight(:,n)
				nperbin(bin(1),bin(2),bin(3)) 	 = nperbin(bin(1),bin(2),bin(3)) + 1
			case(2)
				! Zero mean by adding a 2nd order term preserving the B.C. and adjusting the 
				! 3D parabolicness until the sum is correct
				wsum_bin(bin(1),bin(2),bin(3),:) = wsum_bin(bin(1),bin(2),bin(3),:) + weight(:,n)
				nperbin(bin(1),bin(2),bin(3)) 	 = nperbin(bin(1),bin(2),bin(3)) + 1
        		Na(:) = (r_in_(:)/binsize(:) - dble(bin(:)-2))
				fxfyfz = ((Na(1)-0.5d0)**2-0.25d0) &
						*((Na(2)-0.5d0)**2-0.25d0) &
						*((Na(3)-0.5d0)**2-0.25d0)
				sqr_term(bin(1),bin(2),bin(3)) = sqr_term(bin(1),bin(2),bin(3)) + fxfyfz
			end select
		endif


		if (any(rhat(:,n)-epsilon(rhat(:,n)) .ge. 1.d0) .or. any(rhat(:,n)+epsilon(rhat(:,n)) .le. -1.d0)) then
		    stop "Error in lagrange_poly_weight_Nmol as -1 < rhat < 1 is not true"
		endif
	enddo

	if (present(shiftmean)) then
		select case(shiftmean)
		case(0)
			!Do nothing, shiftmean not requested
            !print'(a,6f18.4)', 'Error in integral cell 3,3,3  ',&
            !                     wsum_bin(3,3,3,1), meanvalue_(3,3,3,1), & 
            !                     wsum_bin(3,3,3,2), meanvalue_(3,3,3,2), & 
            !                     wsum_bin(3,3,3,3), meanvalue_(3,3,3,3) !TEMP FOR DEBUG
		case(1)
			!Zero mean by adjusting constant value
			do n =1,size(r_in,2)
				r_in_(:) = r_in(:,n)+0.5d0*domain(:)
				bin(:) = ceiling((r_in_)/binsize(:))+1
				weight(:,n) = weight(:,n) - (wsum_bin(bin(1),bin(2),bin(3),:) & 
							 -meanvalue_(bin(1),bin(2),bin(3),:))/nperbin(bin(1),bin(2),bin(3))

			enddo
            !print'(a,6f18.4)', 'Error in integral cell 3,3,3  ',&
            !                     wsum_bin(3,3,3,1), meanvalue_(3,3,3,1), & 
            !                     wsum_bin(3,3,3,2), meanvalue_(3,3,3,2), & 
            !                     wsum_bin(3,3,3,3), meanvalue_(3,3,3,3)
		case(2)
			! Zero mean by adding a 2nd order term preserving the B.C. and adjusting the 
			! parabolicness until the sum is correct
			do n =1,size(r_in,2)
				r_in_(:) = r_in(:,n)+0.5d0*domain(:)
				bin(:) = ceiling((r_in_)/binsize(:))+1
				Na(:) = (r_in_(:)/binsize(:) - dble(bin(:)-2))
				fxfyfz = ((Na(1)-0.5d0)**2-0.25d0) &
						*((Na(2)-0.5d0)**2-0.25d0) &
						*((Na(3)-0.5d0)**2-0.25d0)
				weight(:,n) = weight(:,n) - fxfyfz * (wsum_bin(bin(1),bin(2),bin(3),:) & 
							-meanvalue_(bin(1),bin(2),bin(3),:))/sqr_term(bin(1),bin(2),bin(3))
			enddo
           ! print'(a,6f18.4)', 'Error in integral cell 3,3,3  ', & 
           !                     wsum_bin(3,3,3,1), meanvalue_(3,3,3,1), & 
           !                     wsum_bin(3,3,3,2), meanvalue_(3,3,3,2), & 
           !                     wsum_bin(3,3,3,3), meanvalue_(3,3,3,3)

		end select

	endif

end function lagrange_poly_weight_Nmol

!!===================================================
!! Function to return a polynomial weighting function
!! at a given molecular position, based on the
!! array of surface fluxes passed into the function

!function lagrange_poly_weight(array,r_in,binsize,domain,order) result(weight)
!    implicit none

!    integer,dimension(nd),intent(in)    :: order
!    real(kind(0.d0)),dimension(nd),intent(in)    :: r_in,domain,binsize
!    real(kind(0.d0)),dimension(:,:,:,:,:),allocatable,intent(in)    :: array

!    real(kind(0.d0)),dimension(nd)		    :: weight

!	integer							:: npoints,i
!    integer,dimension(3)            :: bin, order_nd
!    real(kind(0.d0)),dimension(3)   :: r_in_, rhat, rhat_,lowerlim,upperlim
!    real(kind(0.d0)),dimension(:,:),allocatable :: grid
!    real(kind(0.d0)),dimension(3,8)	:: nodes,test

!    !Shift to all positive and get bin
!    r_in_ = r_in(:)+0.5d0*domain(:)
!    bin(:) = ceiling(r_in_/binsize(:))+1

!    !Map to bin local coordinate system (0 to +1)
!    rhat(:) = r_in_/binsize(:) - dble(bin(:)-2)
!    rhat_(:) = 1.d0 - rhat(:)

!	!Setup polynomial
!	order_nd(:) = order	!Assumed same order in all dimensions
!	lowerlim = (/ 0.d0,0.d0,0.d0 /)
!	upperlim = (/ 1.d0,1.d0,1.d0 /)

!	call lagrange_interp_nd_size ( nd, order_nd, npoints )
!	allocate (grid(nd,npoints))
!	call lagrange_interp_nd_grid ( nd, order_nd, lowerlim, upperlim, npoints, grid )
!	nodes(:,:) = surface_array_to_nodes(array(bin(1),bin(2),bin(3),:,:))
!	call lagrange_interp_nd_value( nd, order_nd, lowerlim, upperlim, npoints, nodes(1,:), 1, rhat(:), weight(1))
!	call lagrange_interp_nd_value( nd, order_nd, lowerlim, upperlim, npoints, nodes(2,:), 1, rhat(:), weight(2))
!	call lagrange_interp_nd_value( nd, order_nd, lowerlim, upperlim, npoints, nodes(3,:), 1, rhat(:), weight(3))

!!	!TEST POLYNOMIAL FUNCTIONS
!!	test = 0.d0
!!	call lagrange_interp (order_nd, lowerlim, upperlim, nodes(1,:), grid, test(1,:))
!!	call lagrange_interp (order_nd, lowerlim, upperlim, nodes(2,:), grid, test(2,:))
!!	call lagrange_interp (order_nd, lowerlim, upperlim, nodes(3,:), grid, test(3,:))
!!	do i =1,8
!!		print'(i6,9f10.5)',i, grid(:,i), nodes(:,i),test(:,i)
!!	enddo
!!    print'(3i6,12f10.5,2l)', bin, dble(bin(:)-2)*(binsize(:))-0.5d0*domain(:), r_in, rhat, weight(:),any(rhat-epsilon(rhat) .ge. 1.d0),any(rhat+epsilon(rhat) .le. 0.d0)
!!    if (any(rhat-epsilon(rhat) .ge. 1.d0) .or. any(rhat+epsilon(rhat) .le. 0.d0)) then
!!        print'(3f27.15)', rhat
!!        stop "Error in linearweight"
!!    endif

!end function lagrange_poly_weight

!==============================================================================================================
!  _                                       _               _____      _                             _       _ 
! | |                                     (_)             |  __ \    | |                           (_)     | |
! | |     __ _  __ _ _ __ __ _ _ __   __ _ _  __ _ _ __   | |__) |__ | |_   _ _ __   ___  _ __ ___  _  __ _| |
! | |    / _` |/ _` | '__/ _` | '_ \ / _` | |/ _` | '_ \  |  ___/ _ \| | | | | '_ \ / _ \| '_ ` _ \| |/ _` | |
! | |___| (_| | (_| | | | (_| | | | | (_| | | (_| | | | | | |  | (_) | | |_| | | | | (_) | | | | | | | (_| | |
! |______\__,_|\__, |_|  \__,_|_| |_|\__, |_|\__,_|_| |_| |_|   \___/|_|\__, |_| |_|\___/|_| |_| |_|_|\__,_|_|s
!               __/ |                 __/ |                              __/ |                                
!              |___/                 |___/                              |___/                                 
!
!==============================================================================================================

!In the simple case of a 1D quadratic, it's as simple as:

!Input
!       x -- three x positions of known values
!       y -- three values at the x positions y = f(x)
!       xeval -- positions to evaluate y at yeval = f(xeval)
!Ouput  yeval -- value at position x

subroutine quadratic_lagrange_interp(x,y,xeval,yeval)
    implicit none

    real(kind(0.d0)), intent(in) :: xeval
    real(kind(0.d0)), dimension(3),intent(in) :: x, y

    real(kind(0.d0)), intent(out) :: yeval

    yeval =   y(1) * (xeval-x(2))/(x(1)-x(2)) * (xeval-x(3))/(x(1)-x(3)) &
            + y(2) * (xeval-x(1))/(x(2)-x(1)) * (xeval-x(3))/(x(2)-x(3)) &
            + y(3) * (xeval-x(2))/(x(3)-x(2)) * (xeval-x(1))/(x(3)-x(1))

end subroutine quadratic_lagrange_interp


! Example usage of code
!

! containing the points, and the result is the vector V of dimension (N) containing the function values.
!
! Typical usage involves several steps: 
!
! 1) The size of the interpolant grid (nd) is determined by:
!
!        call lagrange_interp_nd_size ( M, ind, nd )
!
! where M is the spatial dimension and ind is an array of size M with the required order of polynomial in each dimension
! NOTE that 2 is a straight line
!
! 2) The interpolant grid (xd, size M by nd) is determined by:
!
!        call lagrange_interp_nd_grid ( M, ind, a, b, nd, xd )
!
! where M is the spatial dimension and ind is an array of size M with the required order of polynomial in each dimension
! A and B of size M with the top and bottom spatial limits and nd is the size of the interpolant grid 
! from lagrange_interp_nd_size
!
! The interpolation function needs data at the data points. It is assumed that this will be supplied 
! by a user specified function of the form
!
!        v = f ( m, n, x )
!
! where M is the spatial dimension, N is the number of points to be evaluated, X is a vector of dimension (M,N) 

!Once the interpolant has been defined, the user is free to evaluate it repeatedly, by specifying NI points XI, and requesting the interpolated values ZI by:
!        call lagrange_interp_nd_value ( m, ind, a, b, nd, zd, ni, xi, zi );
!


!*****************************************************************************80
!! LAGRANGE_INTERP evaluates an ND Lagrange interpolant.
!  Modified:
!    8 May 2014
!  Author:
!    Edward Smith using functions by John Burkardt
!  Parameters:
!    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!    in each dimension.
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
!    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant is 
!    to be evaluated.
!    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the 
!    points XI.
!	 Output, real ( kind = 8 ) XD(M,ND), the points at which data was sampled.
!*****************************************************************************80

subroutine lagrange_interp (n_1d, a, b, zd, xi, zi, xd )
	implicit none

	integer,intent(in) :: n_1d(:)
	real(kind(0.d0)),intent(in) :: a(:), b(:), zd(:), xi(:,:)

	real(kind(0.d0)),intent(out) :: zi(:)
	real(kind(0.d0)),intent(out),optional :: xd(:,:)

	integer :: m, nd, ni
	real(kind(0.d0)), allocatable,dimension(:) :: value,x_1d

	m = size(n_1d); nd = product(n_1d); ni = size(zd)
	call lagrange_interp_nd_value( m, n_1d, a, b, nd, zd, ni, xi, zi)
	if (present(xd)) then
		call lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )
	endif

end subroutine lagrange_interp


!*****************************************************************************80
!! LAGRANGE_INTERP_ND_SIZE sizes an M-dimensional Lagrange interpolant.
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    28 September 2012
!  Author:
!    John Burkardt
!  Parameters:
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!    in each dimension.
!    Output, integer ( kind = 4 ) ND, the number of points in the product grid.
!*****************************************************************************80

subroutine lagrange_interp_nd_size ( m, n_1d, nd )
	implicit none

	integer,intent(in)	::	m

	integer,intent(in)	::	n_1d(m)
	integer,intent(out) ::	nd

	!Determine the number of data points.
	nd = product ( n_1d(1:m) )

end subroutine lagrange_interp_nd_size

!*****************************************************************************80
!! LAGRANGE_INTERP_ND_GRID sets an M-dimensional Lagrange interpolant grid.
!	Licensing:
!		This code is distributed under the GNU LGPL license.
!	Modified:
!		29 September 2012
!	Author:
!		John Burkardt
!	Parameters:
!		Input, integer ( kind = 4 ) M, the spatial dimension.
!		Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!		in each dimension.
!		Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!		Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!		Output, real ( kind = 8 ) XD(M,ND), the points at which data was sampled.
!*****************************************************************************80

subroutine lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )
	implicit none

	integer, intent(in) ::	m, n_1d(m),nd
	real(kind(0.d0)),intent(in) :: a(m), b(m)

	real(kind(0.d0)),intent(out) :: xd(m,nd)

	integer ::	i, n
	real(kind(0.d0)), allocatable :: x_1d(:)

	!Compute the data points.
	xd(1:m,1:nd) = 0.0d0
	do i = 1, m
		n = n_1d(i)
		allocate ( x_1d(1:n) )
		call cc_compute_points ( n, x_1d )
		x_1d(1:n) = 0.5d0 * ( ( 1.0d0 - x_1d(1:n) ) * a(i) &
		                    + ( 1.0d0 + x_1d(1:n) ) * b(i) )
		call r8vec_direct_product ( i, n, x_1d, m, nd, xd )
		deallocate ( x_1d )
	end do

end subroutine lagrange_interp_nd_grid

!*****************************************************************************80
!! LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    28 September 2012
!  Author:
!    John Burkardt
!  Parameters:
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!    in each dimension.
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
!    Input, integer ( kind = 4 ) NI, the number of points at which the 
!    interpolant is to be evaluated.
!    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant is 
!    to be evaluated.
!    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the 
!    points XI.
!*****************************************************************************80

subroutine lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi, zi )
	implicit none

	integer,intent(in) :: m, n_1d(m), nd, ni
	real(kind(0.d0)),intent(in) :: a(m), b(m), zd(nd),xi(m,ni)
	real(kind(0.d0)),intent(out) :: zi(ni)

	integer :: i,j,n
	real(kind(0.d0)) :: w(nd)
	real(kind(0.d0)), allocatable,dimension(:) :: value,x_1d

	do j = 1, ni

		w(1:nd) = 1.0d0

		do i = 1, m
			n = n_1d(i)
			allocate ( x_1d(1:n) )
			allocate ( value(1:n) )
			call cc_compute_points ( n, x_1d )
			x_1d(1:n) = 0.5d0 * ( ( 1.0d0 - x_1d(1:n) ) * a(i) &
			                    + ( 1.0d0 + x_1d(1:n) ) * b(i) )
			call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
			call r8vec_direct_product2 ( i, n, value, m, nd, w )
			deallocate ( value )
			deallocate ( x_1d )
		end do

		zi(j) = dot_product ( w, zd )

	end do

end subroutine lagrange_interp_nd_value

!*****************************************************************************80
!! LAGRANGE_BASIS_1D evaluates a 1D Lagrange basis.
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    09 October 2012
!  Author:
!    John Burkardt
!  Parameters:
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    Input, real ( kind = 8 ) XD(ND), the interpolation nodes.
!    Input, integer ( kind = 4 ) NI, the number of evaluation points.
!    Input, real ( kind = 8 ) XI(NI), the evaluation points.
!    Output, real ( kind = 8 ) LB(NI,ND), the value, at the I-th point XI, 
!    of the Jth basis function.
!*****************************************************************************80

subroutine lagrange_basis_1d ( nd, xd, ni, xi, lb ) 
	implicit none

	integer,intent(in) ::	nd, ni
	real(kind(0.d0)),intent(in) :: xd(nd),xi(ni)

	real(kind(0.d0)),intent(out) :: lb(ni,nd)

	integer ::	i, j
	
	do i = 1, ni
		do j = 1, nd
			lb(i,j) = product ( ( xi(i) - xd(1:j-1)	)  / ( xd(j) - xd(1:j-1)	) ) &
			        * product ( ( xi(i) - xd(j+1:nd) ) / ( xd(j) - xd(j+1:nd) ) )
		end do
	end do

end subroutine lagrange_basis_1d


!*****************************************************************************80
!! CC_COMPUTE_POINTS: abscissas of a Clenshaw Curtis rule.
!	Discussion:
!    Our convention is that the abscissas are numbered from left to right.
!    The rule is defined on [-1,1].
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    08 October 2008
!  Author:
!    John Burkardt
!  Parameters:
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!    Output, real ( kind = 8 ) POINTS(N), the abscissas.
!*****************************************************************************

subroutine cc_compute_points ( n, points )
	implicit none

	integer,intent(in) ::	n
	real(kind(0.d0)),intent(out) :: points(n)

	integer ::	i

	if ( n < 1 ) then

		write ( *, '(a)' ) ' '
		write ( *, '(a)' ) 'CC_COMPUTE_POINTS - Fatal error!'
		write ( *, '(a,i8)' ) '	Illegal value of N = ', n
		stop

	else if ( n == 1 ) then

		points(1) = 0.0d0

	else

		do i = 1, n
			points(i) = cos ( real ( n - i, kind = 8 ) * pi &
			                / real ( n - 1, kind = 8 ) )
		end do

		points(1) = -1.0d0
		if ( mod ( n, 2 ) == 1 ) then
			points((n+1)/2) = 0.0d0
		end if
		points(n) = +1.0d0

	end if

	return
end subroutine cc_compute_points


!*****************************************************************************80
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!	Discussion:
!    An R8VEC is a vector of R8's.
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!    The product rule will be represented as a list of points and weights.
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!    This routine carries out the task involving the weights W.
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!  Example:
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    18 April 2009
!  Author:
!    John Burkardt
!  Parameters:
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!    Input/output, real ( kind = 8 ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!  Local Parameters:
!    Local, integer ( kind = 4 ) START, the first location of a block of values
!    to set.
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values
!    to set.
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!*****************************************************************************80

subroutine r8vec_direct_product ( factor_index, factor_order, factor_value,factor_num, point_num, x )
	implicit none

	integer ::	factor_num
	integer ::	factor_order
	integer ::	point_num

	integer , save :: contig
	integer ::	factor_index
	real ( kind = 8 ) factor_value(factor_order)
	integer ::	j
	integer ::	k
	integer , save :: rep
	integer , save :: skip
	integer ::	start
	real ( kind = 8 ) x(factor_num,point_num)

	if ( factor_index == 1 ) then
		contig = 1
		skip = 1
		rep = point_num
		x(1:factor_num,1:point_num) = 0.0d0
	end if

	rep = rep / factor_order
	skip = skip * factor_order

	do j = 1, factor_order

		start = 1 + ( j - 1 ) * contig

		do k = 1, rep
			x(factor_index,start:start+contig-1) = factor_value(j)
			start = start + skip
		end do

	end do

	contig = contig * factor_order

	return
end subroutine r8vec_direct_product


subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value,factor_num, point_num, w )
	implicit none

	integer ::	factor_num
	integer ::	factor_order
	integer ::	point_num

	integer , save :: contig
	integer ::	factor_index
	real ( kind = 8 ) factor_value(factor_order)
	integer ::	j
	integer ::	k
	integer , save :: rep
	integer , save :: skip
	integer ::	start
	real ( kind = 8 ) w(point_num)

	if ( factor_index == 1 ) then
		contig = 1
		skip = 1
		rep = point_num
		w(1:point_num) = 1.0d0
	end if

	rep = rep / factor_order
	skip = skip * factor_order

	do j = 1, factor_order

		start = 1 + ( j - 1 ) * contig

		do k = 1, rep
			w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
			start = start + skip
		end do

	end do

	contig = contig * factor_order

end subroutine r8vec_direct_product2




! ----------------------------------------------------------------------------
! Numerical diagonalization of 3x3 matrcies
! Copyright (C) 2006  Joachim Kopp
!
! If you use this code, please cite
! Joachim Kopp
! Efficient numerical diagonalization of hermitian 3x3 matrices
! Int. J. Mod. Phys. C 19 (2008) 523-548
! arXiv.org: physics/0610206
!
!
!
! ----------------------------------------------------------------------------
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
! ----------------------------------------------------------------------------


! ----------------------------------------------------------------------------
SUBROUTINE get_eigenval3x3(A, W)
    implicit none
    ! ----------------------------------------------------------------------------
    ! Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
    ! analytical algorithm.
    ! Only the diagonal and upper triangular parts of A are accessed. The access
    ! is read-only.
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   A: The symmetric input matrix
    !   W: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------
    !     .. Arguments ..
	real(kind(0.d0)), intent(in)   :: A(3,3)
	real(kind(0.d0)), intent(out)  :: W(3)

    !     .. Parameters ..
	real(kind(0.d0)) SQRT3
	PARAMETER	  ( SQRT3 = 1.73205080756887729352744634151D0 )

    !     .. Local Variables ..
	real(kind(0.d0)) M, C1, C0
	real(kind(0.d0)) DE, DD, EE, FF
	real(kind(0.d0)) P, SQRTP, Q, C, S, PHI
  
    !     Determine coefficients of characteristic poynomial. We write
    !	     | A   D   F  |
    !	A =  | D*  B   E  |
    !	     | F*  E*  C  |
	DE    = A(1,2) * A(2,3)
	DD    = A(1,2)**2
	EE    = A(2,3)**2
	FF    = A(1,3)**2
	M     = A(1,1) + A(2,2) + A(3,3)
	C1    = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) &
      	   - (DD + EE + FF)
	C0    = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) &
      	   - 2.0D0 * A(1,3)*DE

	P     = M**2 - 3.0D0 * C1
	Q     = M*(P - (3.0D0/2.0D0)*C1) - (27.0D0/2.0D0)*C0
	SQRTP = SQRT(ABS(P))

	PHI   = 27.0D0 * ( 0.25D0 * C1**2 * (P - C1) &
      	    + C0 * (Q + (27.0D0/4.0D0)*C0) )
	PHI   = (1.0D0/3.0D0) * ATAN2(SQRT(ABS(PHI)), Q)

	C     = SQRTP * COS(PHI)
	S     = (1.0D0/SQRT3) * SQRTP * SIN(PHI)

	W(2) = (1.0D0/3.0D0) * (M - C)
	W(3) = W(2) + S
	W(1) = W(2) + C
	W(2) = W(2) - S

END SUBROUTINE get_eigenval3x3


! ----------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
! matrix A using Cardano's method for the eigenvalues and an analytical
! method based on vector cross products for the eigenvectors.
! Only the diagonal and upper triangular parts of A need to contain
! meaningful values. However, all of A may be used as temporary storage
! and may hence be destroyed.
! ----------------------------------------------------------------------------
! Parameters:
!   A: The symmetric input matrix
!   Q: Storage buffer for eigenvectors
!   W: Storage buffer for eigenvalues
! ----------------------------------------------------------------------------
! Dependencies:
!   get_eigenval3x3()
! ----------------------------------------------------------------------------
! Version history:
!   v1.1 (12 Mar 2012): Removed access to lower triangualr part of A
!     (according to the documentation, only the upper triangular part needs
!     to be filled)
!   v1.0: First released version
! ----------------------------------------------------------------------------
SUBROUTINE get_eigenvec3x3(A, Q, W)
    implicit none
    !     .. Arguments ..
	real(kind(0.d0)), intent(inout)  :: A(3,3)
	real(kind(0.d0)), intent(out)    :: Q(3,3)
	real(kind(0.d0)), intent(out)    :: W(3)

    !     .. Parameters ..
	real(kind(0.d0)) EPS
	PARAMETER	  ( EPS = 2.2204460492503131D-16 )

    !     .. Local Variables ..
	real(kind(0.d0)) NORM, N1, N2, N1TMP, N2TMP
	real(kind(0.d0)) THRESH, ERROR, WMAX, F, T
	INTEGER	    I, J

    !     Calculate eigenvalues
	CALL get_eigenval3x3(A, W)

    !     --- The rest of this subroutine can be omitted if only the eigenvalues are desired ---
	WMAX   = MAX(ABS(W(1)), ABS(W(2)), ABS(W(3)))
	THRESH = (8.0D0 * EPS * WMAX)**2

    !     Prepare calculation of eigenvectors
	N1TMP   = A(1, 2)**2 + A(1, 3)**2
	N2TMP   = A(1, 2)**2 + A(2, 3)**2
	Q(1, 1) = A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)
	Q(1, 2) = Q(1, 1)
	Q(2, 1) = A(1, 3) * A(1, 2) - A(2, 3) * A(1, 1)
	Q(2, 2) = Q(2, 1)
	Q(3, 2) = A(1, 2)**2

    !     Calculate first eigenvector by the formula
    !	 v[0] = (A - lambda[0]).e1 x (A - lambda[0]).e2
	A(1, 1) = A(1, 1) - W(1)
	A(2, 2) = A(2, 2) - W(1)
	Q(1, 1) = Q(1, 2) + A(1, 3) * W(1)
	Q(2, 1) = Q(2, 2) + A(2, 3) * W(1)
	Q(3, 1) = A(1, 1) * A(2, 2) - Q(3, 2)
	NORM    = Q(1, 1)**2 + Q(2, 1)**2 + Q(3, 1)**2
	N1	= N1TMP + A(1, 1)**2
	N2	= N2TMP + A(2, 2)**2
	ERROR   = N1 * N2

    !     If the first column is zero, then (1, 0, 0) is an eigenvector
	IF (N1 .LE. THRESH) THEN
	  Q(1, 1) = 1.0D0
	  Q(2, 1) = 0.0D0
	  Q(3, 1) = 0.0D0
    !     If the second column is zero, then (0, 1, 0) is an eigenvector
	ELSE IF (N2 .LE. THRESH) THEN
	  Q(1, 1) = 0.0D0
	  Q(2, 1) = 1.0D0
	  Q(3, 1) = 0.0D0
    !     If angle between A(*,1) and A(*,2) is too small, don't use
    !     cross product, but calculate v ~ (1, -A0/A1, 0)
	ELSE IF (NORM .LT. (64.0D0 * EPS)**2 * ERROR) THEN
        T = ABS(A(1, 2))
        F = -A(1, 1) / A(1, 2)
        IF (ABS(A(2, 2)) .GT. T) THEN
            T = ABS(A(2, 2))
            F = -A(1, 2) / A(2, 2)
        END IF
        IF (ABS(A(2, 3)) .GT. T) THEN
            F = -A(1, 3) / A(2, 3)
        END IF
        NORM    = 1.0D0 / SQRT(1.0D0 + F**2)
        Q(1, 1) = NORM
        Q(2, 1) = F * NORM
        Q(3, 1) = 0.0D0
        !     This is the standard branch
	ELSE
        NORM = SQRT(1.0D0 / NORM)
        do J = 1, 3
        Q(J, 1) = Q(J, 1) * NORM
        enddo
	END IF
 
    !     Prepare calculation of second eigenvector     
	T = W(1) - W(2)

    !     Is this eigenvalue degenerate?
	IF (ABS(T) .GT. 8.0D0 * EPS * WMAX) THEN
    !	 For non-degenerate eigenvalue, calculate second eigenvector by
    !	 the formula
    !	   v[1] = (A - lambda[1]).e1 x (A - lambda[1]).e2
	  A(1, 1) = A(1, 1) + T
	  A(2, 2) = A(2, 2) + T
	  Q(1, 2) = Q(1, 2) + A(1, 3) * W(2)
	  Q(2, 2) = Q(2, 2) + A(2, 3) * W(2)
	  Q(3, 2) = A(1, 1) * A(2, 2) - Q(3, 2)
	  NORM    = Q(1, 2)**2 + Q(2, 2)**2 + Q(3, 2)**2
	  N1	= N1TMP + A(1, 1)**2
	  N2	= N2TMP + A(2, 2)**2
	  ERROR   = N1 * N2

	  IF (N1 .LE. THRESH) THEN
	    Q(1, 2) = 1.0D0
	    Q(2, 2) = 0.0D0
	    Q(3, 2) = 0.0D0
	  ELSE IF (N2 .LE. THRESH) THEN
	    Q(1, 2) = 0.0D0
	    Q(2, 2) = 1.0D0
	    Q(3, 2) = 0.0D0
	  ELSE IF (NORM .LT. (64.0D0 * EPS)**2 * ERROR) THEN
	    T = ABS(A(1, 2))
	    F = -A(1, 1) / A(1, 2)
	    IF (ABS(A(2, 2)) .GT. T) THEN
		T = ABS(A(2, 2))
		F = -A(1, 2) / A(2, 2)
	    END IF
	    IF (ABS(A(2, 3)) .GT. T) THEN
		F = -A(1, 3) / A(2, 3)
	    END IF
	    NORM    = 1.0D0 / SQRT(1.0D0 + F**2)
	    Q(1, 2) = NORM
	    Q(2, 2) = F * NORM
	    Q(3, 2) = 0.0D0
	  ELSE
	    NORM = SQRT(1.0D0 / NORM)
	    do J = 1, 3 
    		Q(J, 2) = Q(J, 2) * NORM
        enddo
	  END IF
	ELSE
    !	 For degenerate eigenvalue, calculate second eigenvector according to
    !	   v[1] = v[0] x (A - lambda[1]).e[i]
    !   
    !	 This would really get to complicated if we could not assume all of A to
    !	 contain meaningful values.
        A(2, 1) = A(1, 2)
        A(3, 1) = A(1, 3)
        A(3, 2) = A(2, 3)
        A(1, 1) = A(1, 1) + W(1)
        A(2, 2) = A(2, 2) + W(1)
        do I = 1, 3
            A(I, I) = A(I, I) - W(2)
            N1	= A(1, I)**2 + A(2, I)**2 + A(3, I)**2
            IF (N1 .GT. THRESH) THEN
                Q(1, 2) = Q(2, 1) * A(3, I) - Q(3, 1) * A(2, I)
                Q(2, 2) = Q(3, 1) * A(1, I) - Q(1, 1) * A(3, I)
                Q(3, 2) = Q(1, 1) * A(2, I) - Q(2, 1) * A(1, I)
                NORM    = Q(1,2)**2 + Q(2,2)**2 + Q(3,2)**2
                if (NORM .GT. (256.0D0 * EPS)**2 * N1) THEN
                    NORM = SQRT(1.0D0 / NORM)
                    do J = 1, 3
                        Q(J, 2) = Q(J, 2) * NORM
                    enddo
                    exit
                    !GO TO 60
                endif
            endif
        enddo
   
        ! This means that any vector orthogonal to v[0] is an EV.
        if (I .eq. 4) THEN
	        do J = 1, 3
            ! Find nonzero element of v[0] and swap it with the next one
		        if (Q(J, 1).NE.0.0D0) THEN
                    NORM = 1.0D0 / SQRT(Q(J, 1)**2 + Q(1 + MOD(J,3), 1)**2)
                    Q(J, 2)		  = Q(1 + MOD(J,3), 1) * NORM
                    Q(1 + MOD(J,3), 2)   = -Q(J, 1) * NORM
                    Q(1 + MOD(J+1,3), 2) = 0.0D0
                    exit
		            !GO TO 80
		        endif
            enddo
        endif
	endif

    !     Calculate third eigenvector according to
    !	 v[2] = v[0] x v[1]
    Q(1, 3) = Q(2, 1) * Q(3, 2) - Q(3, 1) * Q(2, 2)
    Q(2, 3) = Q(3, 1) * Q(1, 2) - Q(1, 1) * Q(3, 2)
    Q(3, 3) = Q(1, 1) * Q(2, 2) - Q(2, 1) * Q(1, 2)

END SUBROUTINE get_eigenvec3x3

end module librarymod



module circle_rectangle_intersection

    real(kind(0.d0)), parameter :: pi=4.d0*atan(1.d0)

contains

    function segment_area(ixs,r) result (area)

        real(kind(0.d0)), intent(in) :: ixs(2,2), r
        real(kind(0.d0)) :: area

        real(kind(0.d0)) :: theta

        theta = acos(dot_product(ixs(1,:),ixs(2,:))/(r**2.d0))
        area = 0.5d0*(r**2.d0)*(theta - sin(theta))
    
    end function segment_area

    function single_triangle_area(ixs,vx,vy,vertices_in) result(area)

        real(kind(0.d0)), intent(in) :: ixs(2,2)
        real(kind(0.d0)), intent(in) :: vx(4), vy(4)
        logical, intent(in) :: vertices_in(4)
        real(kind(0.d0)) :: area

        integer :: i
        real(kind(0.d0)) :: px(3), py(3)

        ! Store which vertex is in the radius        
        do i = 1,4

            if (vertices_in(i)) then

                px = (/vx(i), ixs(1,1), ixs(2,1)/)
                py = (/vy(i), ixs(1,2), ixs(2,2)/)
                area = triangle_area(px,py) 
                return

            end if

        end do

    end function single_triangle_area

    function single_triangle_outside_area(ixs,vx,vy,vertices_in) result(area)

        real(kind(0.d0)), intent(in) :: ixs(2,2)
        real(kind(0.d0)), intent(in) :: vx(4), vy(4)
        logical, intent(in) :: vertices_in(4)
        real(kind(0.d0)) :: area

        integer :: i
        real(kind(0.d0)) :: px(3), py(3)

        ! Store which vertex is in the radius        
        do i = 1,4

            if (vertices_in(i) .eqv. .false.) then

                px = (/vx(i), ixs(1,1), ixs(2,1)/)
                py = (/vy(i), ixs(1,2), ixs(2,2)/)
                area = triangle_area(px,py) 
                return

            end if

        end do

    end function single_triangle_outside_area

    function double_triangle_area(ixs,vx,vy,vertices_in) result(area)

        real(kind(0.d0)), intent(in) :: ixs(2,2)
        real(kind(0.d0)), intent(in) :: vx(4), vy(4)
        logical, intent(in) :: vertices_in(4)
        real(kind(0.d0)) :: area

        integer :: i, cnt, vertices(2)
        real(kind(0.d0)) :: px(3), py(3)
        real(kind(0.d0)) :: rvertex(2,2) 
        real(kind(0.d0)) :: first(2), second(2), third(2), fourth(2), temp(2)

        cnt = 0
        do i = 1,4
            if (vertices_in(i)) then
                cnt = cnt + 1
                vertices(cnt) = i 
            end if
        end do

        do cnt=1,2
            vertex = vertices(cnt)
            rvertex(cnt,:) = (/vx(vertex),vy(vertex)/)
        end do

        first = rvertex(1,:)
        second = ixs(1,:)
        third = rvertex(2,:)
        fourth = ixs(2,:)

        ! First and second must always be opposed 
        if ( (first(1) .eq. second(1)) .or. &
             (first(2) .eq. second(2)) ) then
            ! Swap second and fourth
            temp = second
            second = fourth 
            fourth = temp
        end if

        px = (/first(1), second(1), third(1)/)
        py = (/first(2), second(2), third(2)/)
        A1 = triangle_area(px,py)
        px = (/first(1), second(1), fourth(1)/)
        py = (/first(2), second(2), fourth(2)/)
        A2 = triangle_area(px,py)
        
        area = A1 + A2

    end function double_triangle_area

    function triangle_area(px, py) result(area)

        real(kind(0.d0)), intent(in) :: px(3), py(3) !Input points
        real(kind(0.d0)) :: area

        real(kind(0.d0)) :: u(2), v(2)

        u = (/px(2)-px(1), py(2)-py(1)/) 
        v = (/px(3)-px(1), py(3)-py(1)/) 
        area = 0.5*abs(u(1)*v(2) - u(2)*v(1))

    end function triangle_area

    subroutine circle_rectangle_intersection_area(vx, vy, r, npx, npy, A)
        implicit none

        real(kind(0.d0)), intent(in) :: vx(4), vy(4), r
        integer, intent(in) :: npx, npy
        real(kind(0.d0)), intent(out) :: A
        
        real(kind(0.d0)) :: ixs(8,2), A_seg, A_tri, A_dom
        logical :: vertices_in(4) 
        integer :: n_in, n_ixs

        call vertices_in_radius(vx,vy,r,vertices_in, n_in)
        call solve_block_intersections(vx,vy,r,ixs,n_ixs) 
        A = 0.d0

        if (n_in .eq. 0) then

            if (n_ixs .eq. 0) then

                if (npx .eq. 1 .and. npy .eq. 1) then
                    ! Whole circle inside rectangle
                    A = pi * r * r
                else
                    ! Circle not wholly contained
                    A = 0.0
                end if

            else

                stop 'Intersections appeared when n_in = 0. This is too complicated for me to calculate. '

            end if
    

        else if (n_in .eq. 1) then


            if (n_ixs .ne. 2) then
                stop 'One vertex inside circle but n intersections not 2!'
            else
                ! Corner only inside circle
                A_seg = segment_area(ixs(1:2,:),r)
                A_tri = single_triangle_area(ixs(1:2,:),vx,vy,vertices_in)
                A = A_seg + A_tri
            end if


        else if (n_in .eq. 2) then


            if (n_ixs .ne. 2) then
                stop 'Two vertices inside circle but n intersections not 2!'
            else
                ! Two vertices inside
                A_seg = segment_area(ixs(1:2,:),r)
                A_tri = double_triangle_area(ixs(1:2,:),vx,vy,vertices_in) 
                A = A_seg + A_tri
            end if


        else if (n_in .eq. 3) then

            if (n_ixs .ne. 2) then
                stop 'Three vertices inside circle but n intersections not 2!'
            else
                ! Two vertices inside
                A_seg = segment_area(ixs(1:2,:),r)
                A_tri = single_triangle_outside_area(ixs(1:2,:),vx,vy,vertices_in) 
                A = A_seg + A_tri
                !A_dom = domain(1)*domain(2)
                A_dom = (vx(3) - vx(1)) * (vy(2) - vy(1))
                A = A_dom - A_tri + A_seg
            end if
    
        else if (n_in .eq. 4) then
            
            A = (vx(3) - vx(1)) * (vy(2) - vy(1))

        else
            
            stop 'Unknown number of vertices within radius'
 
        end if

    end subroutine circle_rectangle_intersection_area
    
    subroutine solve_block_intersections(vx,vy,r,ixs,n_ixs)
        implicit none

        real(kind(0.d0)), intent(in) :: vx(4), vy(4), r 
        real(kind(0.d0)), intent(out) :: ixs(8,2)
        integer, intent(out) :: n_ixs

        real(kind(0.d0)) :: top, bottom, left, right
        real(kind(0.d0)) :: ixs_line(2,2)
        integer :: n_ixs_line, n

        n_ixs = 0
        ixs(:,:) = -666.d0

        top = vy(2)
        bottom = vy(1)
        left = vx(1)
        right = vx(3)

        ! Top face
        call solve_line_intersections(left, right, top, top, r, &
                                     ixs_line, n_ixs_line)
        call store_intersections

        ! Bottom face
        call solve_line_intersections(left, right, bottom, bottom, r, &
                                     ixs_line, n_ixs_line) 
        call store_intersections

        ! Left face
        call solve_line_intersections(left, left, bottom, top, r, &
                                     ixs_line, n_ixs_line) 
        call store_intersections

        ! Right face
        call solve_line_intersections(right, right, bottom, top, r, &
                                     ixs_line, n_ixs_line) 
        call store_intersections

        !ixs = ixs_top + ixs_bottom + ixs_left + ixs_right 
        
    contains
    
        subroutine store_intersections
            implicit none

            integer :: n

            if (n_ixs_line .gt. 0) then
                n_ixs = n_ixs + n_ixs_line 
                do n=1, n_ixs_line
                    ixs(n_ixs,:) = ixs_line(n,:)
                end do 
            end if
            
        end subroutine store_intersections

    end subroutine solve_block_intersections
       
    subroutine solve_line_intersections(x1, x2, y1, y2, r, line_ixs, n_ixs)
        implicit none

        real(kind(0.d0)), intent(in) :: x1, x2, y1, y2, r
        real(kind(0.d0)), intent(out) :: line_ixs(2,2) 
        integer, intent(out) :: n_ixs
        
        real(kind(0.d0)) :: dx, dy, dr, D, delta, eps
        real(kind(0.d0)) :: xplus, xminu, yplus, yminu
        real(kind(0.d0)) :: x, y, ix(2), xs(2), ys(2)
        integer :: loopx, loopy, i, j

        eps = 0.001 ! Small floating point error

        dx = x2 - x1 
        dy = y2 - y1 
        dr = sqrt(dx**2.0 + dy**2.0)
        D = x1*y2 - x2*y1 
        delta = (r*dr)**2.0 - D**2.0

        n_ixs = 0
        line_ixs(:,:) = -666.d0

        if (delta .lt. 0) then

            return 
        
        else if (delta == 0.0) then

            stop 'Incident boundary on cylinder, cannot compute.'

        else

            xplus = (D*dy + signum(dy)*dx*sqrt(delta)) / (dr**2.0)
            xminu = (D*dy - signum(dy)*dx*sqrt(delta)) / (dr**2.0)
            yplus = (-D*dx + abs(dy)*sqrt(delta)) / (dr**2.0)
            yminu = (-D*dx - abs(dy)*sqrt(delta)) / (dr**2.0)

            xs = (/xplus, xminu/)
            ys = (/yplus, yminu/)

            loopx = 2
            loopy = 2
            if (xs(1) .eq. xs(2)) loopx = 1
            if (ys(1) .eq. ys(2)) loopy = 1

            n_ixs = 0
            do i = 1,loopx 
            do j = 1,loopy

                x = xs(i)
                y = ys(j)
                ix = (/x,y/)

                if ( (x .ge. x1 - eps .and. x .le. x2 + eps) .and. &
                     (y .ge. y1 - eps .and. y .le. y2 + eps) ) then
                    n_ixs = n_ixs + 1
                    line_ixs(n_ixs,:) = ix
                end if

            end do
            end do

        end if

    end subroutine solve_line_intersections

    function signum(a) result(sig)
        implicit none
       
        real(kind(0.d0)), intent(in) :: a
        real(kind(0.d0)) :: sig

        if (a .lt. 0.d0) then
            sig = -1.d0
        else
            sig = 1.d0
        end if 

    end function signum

    subroutine vertices_in_radius(vx, vy, r, v_in, n_in)
        
        real(kind(0.d0)), intent(in) :: vx(4), vy(4), r
        logical, intent(out) :: v_in(4)
        integer, intent(out) :: n_in

        real(kind(0.d0)) :: x, y 
   
        v_in(:) = .false. 
        n_in = 0
        do vertex = 1,4
            x = vx(vertex)
            y = vy(vertex)
            if ((x**2.0 + y**2.0) .le. r**2.0) then
                v_in(vertex) = .true.   
                n_in = n_in + 1
            end if
        end do

    end subroutine vertices_in_radius

end module circle_rectangle_intersection

!======================================================================
!Ensures functions work correctly
subroutine check
	use librarymod
	implicit none

	integer						:: i,j
	integer						:: n,npoints
	integer,dimension(3)		:: indx
	real(kind(0.d0))			:: d
	real(kind(0.d0))			:: integral
	real(kind(0.d0))			:: rand, x_interval, intercept, gradient
	real(kind(0.d0)),dimension(3)			:: b
	real(kind(0.d0)),dimension(3,3)			:: a
	real(kind(0.d0)), dimension(:), allocatable	:: y

	!Check least squares line

	x_interval = 0.1d0
	npoints = 20
	allocate(y(npoints))

	do n = 1, npoints
		call random_number(rand)
		y(n) = 2*(n * x_interval)! + (rand-0.5d0))
	enddo

	!call least_squares(y, x_interval , npoints ,intercept, gradient)

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
