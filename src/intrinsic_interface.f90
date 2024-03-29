

!A module to provide object oriented version of intrinsic surface functions

module intrinsic_module

    double precision, parameter :: pi=3.141592653589793d0

    type :: intrinsic_surface_mock

        integer      :: normal, ixyz, jxyz, topbot
        integer, dimension(2)  :: modes_shape
        !integer, dimension(:), allocatable     :: pivots

        double precision  :: alpha, eps
        double precision, dimension(3) :: box

	    contains
		    procedure :: initialise  => initialise_mock
		    procedure :: update_surface => update_surface_mock
		    procedure :: get_surface => get_surface_mock

    end type intrinsic_surface_mock


    type :: intrinsic_surface_complex

        integer      :: normal, ixyz, jxyz, topbot
        integer, dimension(2)  :: modes_shape
        !integer, dimension(:), allocatable     :: pivots

        double precision  :: alpha, eps
        double precision, dimension(3) :: box
        double precision, dimension(:), allocatable     :: Q
        double precision, dimension(:,:), allocatable   :: Qxy
        double complex, dimension(:,:), allocatable     :: modes
        double precision, dimension(:,:,:), allocatable :: q_vectors

	    contains
		    procedure :: initialise  => compute_q_vectors
		    procedure :: update_surface => update_surface_modes
		    procedure :: get_surface => get_surface_complex
!		    procedure :: sample_surface => sample_intrinsic_surface

    end type intrinsic_surface_complex

    type :: intrinsic_surface_real

        integer      :: normal, ixyz, jxyz, topbot, n_waves2
		integer, dimension(2)  :: n_waves, qm
		integer, dimension(3)  :: nbins, nhb, periodic

        integer, dimension(:), allocatable     :: u, v
        integer, dimension(:), allocatable     :: pivots

        double precision  :: alpha, eps, area
        double precision, dimension(3) :: box, binsize
        double precision, dimension(:), allocatable     :: diag, coeff, coeff_mdt
        double precision, dimension(:,:), allocatable   :: diag_matrix, intrnsc_smple
        double precision, dimension(:,:,:,:), allocatable :: Abilinear


	    contains
		    procedure :: initialise  => initialise
			procedure :: get_A_b => get_A_b_fourier
		    procedure :: update_surface => update_real_surface
		    procedure :: get_surface => get_real_surface
		    procedure :: get_surface_derivative => get_real_surface_derivative
		    procedure :: get_metric_tensor => get_real_metric_tensor
		    procedure :: get_zero_mode => get_zero_mode
            procedure :: get_bin => get_bin_from_surface
			procedure :: get_tangent_bins => get_tangent_bins
		    procedure :: write_modes => write_modes
		    procedure :: sample_surface => sample_intrinsic_surface
		    procedure :: get_sampled_surface => get_sampled_surface
			procedure :: fit_intrinsic_surface => fit_intrinsic_surface_modes
			procedure :: index_to_vertex => index_to_vertex
			procedure :: get_crossings => get_crossings_real
			procedure :: get_flat_crossings => get_flat_crossings
		    procedure :: intrinsic_area => intrinsic_area
            procedure :: apply_metric_tensor_transform => apply_metric_tensor_transform

    end type intrinsic_surface_real
	
	type, extends(intrinsic_surface_real) :: intrinsic_surface_bilinear

		contains
			!Over ride existing routines
		    procedure :: get_surface => get_bilinear_surface
		    procedure :: get_surface_derivative => get_bilinear_surface_derivative
		    procedure :: get_metric_tensor => get_bilinear_metric_tensor
			procedure :: fit_intrinsic_surface => fit_intrinsic_surface_bilinear
			procedure :: get_crossings => get_crossings_bilinear

			!New specialist routines
    		procedure :: update_sampled_surface => update_sampled_surface!_opt
			procedure :: get_surface_bilinear => get_surface_bilinear
			procedure :: paramterise_bilinear_surface => paramterise_bilinear_surface
			procedure :: get_surface_derivative_from_bilinear => get_surface_derivative_from_bilinear
			procedure :: indices_to_points => indices_to_points

			procedure ::get_bilinear_patch_area => get_bilinear_patch_area
			procedure ::intrinsic_area_bilinear => intrinsic_area_bilinear

	end type intrinsic_surface_bilinear

	type, extends(intrinsic_surface_bilinear) :: intrinsic_surface_chebychev

		contains

			procedure :: get_A_b => get_A_b_chebychev
			procedure :: get_fuv => get_fuv
		    !procedure :: get_surface => get_chebychev_surface

	end type intrinsic_surface_chebychev

	!An object which allows localised surface calculations
	! type, extends(intrinsic_surface_real) :: intrinsic_surface_real_binwidth
		! contains
		    ! procedure :: get_surface => get_real_surface_binwidth
	! end type intrinsic_surface_real_binwidth

    !Wave function using array or single value
    interface wave_function
        module procedure wave_function_single, &
                         wave_function_multipoint, & 
                         wave_function_multipointandwave
    end interface

    private wave_function_single, &
            wave_function_multipoint, &
            wave_function_multipointandwave
contains 


subroutine initialise_mock(self, box, normal, alpha, eps, topbot)
    implicit none

	class(intrinsic_surface_mock) :: self

    integer, intent(in)         :: normal, topbot
    double precision, intent(in) :: alpha, eps
    double precision, intent(in), dimension(3) :: box

    integer :: i, j

    self%normal = normal
    self%box = box
    self%alpha = alpha
    self%eps = eps
    self%ixyz = mod(normal,3)+1
    self%jxyz = mod(normal+1,3)+1
    self%topbot = topbot

end subroutine initialise_mock


subroutine update_surface_mock(self, points)
    implicit none

	class(intrinsic_surface_mock) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points

end subroutine update_surface_mock

subroutine get_surface_mock(self, points, elevation)
    implicit none

	class(intrinsic_surface_mock) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(out), dimension(:), allocatable :: elevation

    allocate(elevation(size(points,1)))
    elevation = 0.d0

end subroutine get_surface_mock


!=========================================
!=========================================
!   Adapted from Alias surface code
!=========================================
!=========================================

!Function to allow domains where Lx != Ly
function check_uv(u, v)
    implicit none

    integer,dimension(:),intent(in)     :: u, v

    real(kind(0.d0)),dimension(size(u))    :: check_uv

    integer :: i

    do i =1,size(u,1)
        if ((abs(u(i))+abs(v(i))) .eq. 0) then
            check_uv(i) = 4.d0
        elseif (u(i)*v(i) .eq. 0) then
            check_uv(i) = 2.d0
        else
            check_uv(i) = 1.d0
        endif
    enddo

end function check_uv

!Functions for Fourier surface
function wave_function_single(x, u, Lx)
    implicit none

    integer,intent(in) :: u
    real(kind(0.d0)),intent(in) :: Lx
    real(kind(0.d0)),intent(in) :: x

    real(kind(0.d0))   :: wave_function_single

	if (u .ge. 0) then
		wave_function_single = cos(2.d0 * pi * u * x / Lx)
	else
		wave_function_single = sin(2.d0 * pi * abs(u) * x / Lx)
	endif

end function wave_function_single

function wave_function_multipoint(x, u, Lx)
    implicit none

    integer,intent(in) :: u
    real(kind(0.d0)),intent(in) :: Lx
    real(kind(0.d0)),intent(in),dimension(:) :: x

    real(kind(0.d0)),dimension(size(x))    :: wave_function_multipoint
    integer :: i

    do i =1,size(x,1)
        if (u .ge. 0) then
            wave_function_multipoint(i) = cos(2.d0 * pi * u * x(i) / Lx)
        else
            wave_function_multipoint(i) = sin(2.d0 * pi * abs(u) * x(i) / Lx)
        endif
    enddo

end function wave_function_multipoint


function wave_function_multipointandwave(x, u, Lx)
    implicit none

    integer,intent(in),dimension(:) :: u
    real(kind(0.d0)),intent(in) :: Lx
    real(kind(0.d0)),intent(in),dimension(:) :: x

    real(kind(0.d0)),dimension(size(x),size(u))  :: wave_function_multipointandwave

    integer :: i, j

    do i =1,size(x,1)
        do j =1,size(u,1)
            if (u(j) .ge. 0) then
                wave_function_multipointandwave(i,j) = cos(2.d0 * pi * u(j) * x(i) / Lx)
            else
                wave_function_multipointandwave(i,j) = sin(2.d0 * pi * abs(u(j)) * x(i) / Lx)
            endif
        enddo
    enddo

end function wave_function_multipointandwave


function derivative_wave_function(x, u, Lx)
    implicit none

    integer,intent(in) :: u
    real(kind(0.d0)),intent(in) :: Lx
    real(kind(0.d0)),intent(in),dimension(:) :: x

    real(kind(0.d0)),dimension(size(x))    :: derivative_wave_function
    integer :: i

    do i =1,size(x,1)
        if (u .ge. 0) then
            derivative_wave_function(i) = - (2.d0 * pi * u / Lx) & 
                                        *sin(2.d0 * pi * u * x(i) / Lx)
        else
            derivative_wave_function(i)=   (2.d0 * pi * abs(u) / Lx) & 
                                       *cos(2.d0 * pi * abs(u) * x(i) / Lx)
        endif
    enddo

end function derivative_wave_function


function chebyshev_function(x, u, Lx) result(T)
    implicit none

    integer,intent(in),dimension(:) :: u
    real(kind(0.d0)),intent(in) :: Lx
    real(kind(0.d0)),intent(in),dimension(:) :: x

    real(kind(0.d0)),dimension(size(x),size(u))  :: T

    integer :: n, M, un, unm1, unm2
	real(kind(0.d0)), dimension(:), allocatable :: xs
    real(kind(0.d0)),dimension(size(x))  :: Tm1, Tm2

	M = size(u)

	!Rescale to [-1, 1] from -Lx/2 to Lx/2
	allocate(xs(size(x)))
	xs = 2.d0 * x(:) / Lx

	do n=1,M
		un = u(n)
		!print*, n, un, unm1, unm2, un .eq. 0, un .eq. 1, un .eq. unm1
		select case(un)
		case (0)
			T(:,n) = 1
			Tm2 = T(:,n)
			unm2 = un
		case (1)
			T(:,n) = xs(:)
			Tm1 = T(:,n)
			unm1 = un
		case default
			if (un .eq. unm1) then
				T(:,n) = T(:,n-1)
			else
				T(:,n) = 2*xs(:)*Tm1 - Tm2
				Tm2 = Tm1	
				Tm1 = T(:,n)
				unm2 = unm1
			endif
			unm1 = un
		end select
        !Add scaling factor to make orthogonal
        !T(:,n) = T(:,n)*sqrt(1-xs(:)**2)
	enddo
	!Recurrence formula is used so advice is to always get all 
	!points together at same time for efficiency
	! T(1,:) = 1.d0
	! T(2,:) = xs(:)
	! do n=2,M-1
		! !print*, "chebyshev_function", n, u(n)
		! T(n+1,:) = 2.d0*xs(:)*T(n,:) - T(n-1,:)
	! enddo

end function chebyshev_function


!Start of intitialistion functions

subroutine initialise(self, box, normal, alpha, eps, & 
                      nbins, nhb, topbot, periodic)
    implicit none

	class(intrinsic_surface_real) :: self

    integer, intent(in)         :: normal, topbot
	
    integer, intent(in), dimension(3)  :: nbins, nhb, periodic
    double precision, intent(in) :: alpha, eps
    double precision, intent(in), dimension(3) :: box

    integer :: i, j

    self%normal = normal
    self%box = box
    self%alpha = alpha
    self%eps = eps
    self%ixyz = mod(normal,3)+1
    self%jxyz = mod(normal+1,3)+1
	self%nbins = nbins
	self%nhb = nhb
	self%binsize = box/float(nbins)
    self%topbot = topbot
    self%periodic = periodic

    self%area = box(self%ixyz)*box(self%jxyz)
    self%qm(1) = int(box(self%ixyz)/alpha)
    self%qm(2) = int(box(self%jxyz)/alpha)
    self%n_waves(1) = 2*self%qm(1)+1
    self%n_waves(2) = 2*self%qm(2)+1
    self%n_waves2 = self%n_waves(1)*self%n_waves(2)

    !print*, "initialise intrinsic function", box, normal, alpha, eps, & 
    !        nbins, nhb, topbot, self%qm, self%n_waves

    !Assume domain is square and nwaves same in both tangential directions
    if (abs(box(self%ixyz)-box(self%jxyz)) .gt. 1e-1) then
        print'(a,2(i6,a,f10.5))', "Intrinsic Inteface - Domain not square, using ", & 
            self%n_waves(1), " modes with sidelength ", box(self%ixyz), &
            self%n_waves(2), " modes with sidelength ", box(self%jxyz)
    endif

    allocate(self%u(self%n_waves2), self%v(self%n_waves2))
	allocate(self%coeff(self%n_waves2))
    allocate(self%diag(self%n_waves2))
	self%coeff = 0.d0

	!Note default behaviour in select type is 
	!to priorities type over class regardless of order
	!(i.e. if type b is derived from a, it will still
	! branch to type(b) even if class(a) is before)
    select type (self)
	type is (intrinsic_surface_chebychev)
		!Chebychev 
		do i=1,self%n_waves2
            if (self%periodic(self%ixyz) .eq. 1) then
			    self%u(i) = (i-1) / self%n_waves(2) - self%qm(1)
            else
    			self%u(i) = (i-1) / self%n_waves(2)
            endif
            if (self%periodic(self%jxyz) .eq. 1) then
			    self%v(i) = modulo((i-1), self%n_waves(2)) - self%qm(2)
            else
    			self%v(i) = modulo((i-1), self%n_waves(2))
            endif
			!print*, "setup u, v for chebychev", i, self%u(i), self%v(i)
		enddo
    class is (intrinsic_surface_real)
		!Numbers between -qm and +qm
		do i=1,self%n_waves2
			self%u(i) = (i-1) / self%n_waves(2) - self%qm(1)
			self%v(i) = modulo((i-1), self%n_waves(2)) - self%qm(2)
			!print*, "setup u, v for Fourier", i, self%u(i), self%v(i)
		enddo
	end select

	!Factor to account for domain not square from 
	!Longford et al (2018 J. Chem. Phys. 149, 234705 (2018); 
	!https://doi.org/10.1063/1.5055241
	self%diag = check_uv(self%u, self%v) & 
				 * (  self%u**2 * box(self%jxyz) / box(self%ixyz) & 
					+ self%v**2 * box(self%ixyz) / box(self%jxyz))
	
    allocate(self%diag_matrix(self%n_waves2,self%n_waves2))
    self%diag_matrix = 0.d0
    do i=1,self%n_waves2
    do j=1,self%n_waves2
        if (i .eq. j) then
            self%diag_matrix(i,j) = 4.d0 * pi**2 * eps * self%diag(i)
        endif
    enddo
    enddo

    select type (self)
    class is (intrinsic_surface_real)
	!	pass
    !type is (intrinsic_surface_complex)
	!	stop "Error - intrinsic_surface initialise -  Complex code is depricated"
    class is (intrinsic_surface_bilinear)
		allocate(self%intrnsc_smple(self%nbins(self%ixyz)+1+2*self%nhb(self%ixyz), & 
									self%nbins(self%jxyz)+1+2*self%nhb(self%jxyz)))
		allocate(self%Abilinear(2,2,self%nbins(self%ixyz)+2*self%nhb(self%ixyz), & 
									self%nbins(self%jxyz)+2*self%nhb(self%jxyz)))
    class default
      ! give error for unexpected/unsupported type
		stop "Error - intrinsic_surface initialise - unexpected type for intrinsic_surface"
    end select

end subroutine initialise

subroutine get_A_b_fourier(self, points, A, b)
#if USE_LAPACK
    use lapack_fns, only : multiply_by_tranpose
#endif
    implicit none

	class(intrinsic_surface_real) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points

    double precision, intent(out), dimension(:), allocatable :: b
    double precision, intent(out), dimension(:,:), allocatable :: A

    integer :: j
    double precision, dimension(:,:), allocatable :: fuv

    allocate(A(self%n_waves2, self%n_waves2))
    allocate(b(self%n_waves2))
    allocate(fuv(size(points,1), self%n_waves2))

    A = 0.d0
    b = 0.d0
    do j =1, self%n_waves2
        fuv(:,j) = wave_function(points(:,self%ixyz), self%u(j), self%box(self%ixyz)) & 
                  *wave_function(points(:,self%jxyz), self%v(j), self%box(self%jxyz))
        b(j) = b(j) + sum(points(:,self%normal) * fuv(:,j))
    enddo

    A = 0.d0
#if USE_LAPACK
    call multiply_by_tranpose(fuv, fuv, A)
#else
	A = matmul(transpose(fuv), fuv)
    !call error_abort("Error - must build with lapack (p_lapack or p_sys_lapack) to use intrinsic interface")
#endif

	!Add constraint on minimum size
    A = A + self%diag_matrix

end subroutine get_A_b_fourier


subroutine get_fuv(self, points, fuv)
    implicit none

	class(intrinsic_surface_chebychev) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(out), dimension(:,:), allocatable :: fuv

    allocate(fuv(size(points,1), self%n_waves2))

    !Tensor product of two surfaces
    !fuv = fu * fv
    if (self%periodic(self%ixyz) .eq. 1) then
    	fuv = wave_function(points(:,self%ixyz), self%u(:), self%box(self%ixyz))
    else
    	fuv = chebyshev_function(points(:,self%ixyz), self%u(:), self%box(self%ixyz))
    endif
    if (self%periodic(self%jxyz) .eq. 1) then
	    fuv = fuv*wave_function(points(:,self%jxyz), self%v(:), self%box(self%jxyz))
    else
	    fuv = fuv*chebyshev_function(points(:,self%jxyz), self%v(:), self%box(self%jxyz))
    endif

end subroutine get_fuv

subroutine get_A_b_chebychev(self, points, A, b)
#if USE_LAPACK
    use lapack_fns, only : multiply_by_tranpose
#endif
    implicit none

	class(intrinsic_surface_chebychev) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points

    double precision, intent(out), dimension(:), allocatable :: b
    double precision, intent(out), dimension(:,:), allocatable :: A

    integer :: j
    double precision, dimension(:,:), allocatable :: fuv

    allocate(A(self%n_waves2, self%n_waves2))
    allocate(b(self%n_waves2))
    !allocate(fuv(size(points,1), self%n_waves2))

    !print*, "get_A_b_chebychev, number of points", size(points,1), "number of coeffs",  self%n_waves2

    A = 0.d0
    b = 0.d0
    call self%get_fuv(points, fuv)
	!fuv = chebyshev_function(points(:,self%ixyz), self%u(:), self%box(self%ixyz)) & 
	!	 *chebyshev_function(points(:,self%jxyz), self%v(:), self%box(self%jxyz))
    do j =1, self%n_waves2
		b(j) = b(j) + sum(points(:,self%normal) * fuv(:,j))
	enddo

#if USE_LAPACK
    call multiply_by_tranpose(fuv, fuv, A)
#else
	A = matmul(transpose(fuv), fuv)
    !call error_abort("Error - must build with lapack (p_lapack or p_sys_lapack) to use intrinsic interface")
#endif

	!Add constraint on minimum size
    A = A + self%diag_matrix

end subroutine get_A_b_chebychev



subroutine get_A_b_mixed(self, points, A, b)
#if USE_LAPACK
    use lapack_fns, only : multiply_by_tranpose
#endif
    implicit none

	class(intrinsic_surface_chebychev) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points

    double precision, intent(out), dimension(:), allocatable :: b
    double precision, intent(out), dimension(:,:), allocatable :: A

    integer :: j
    double precision, dimension(:,:), allocatable :: fu, fv, fuv

    allocate(A(self%n_waves2, self%n_waves2))
    allocate(b(self%n_waves2))

    A = 0.d0
    b = 0.d0
    call self%get_fuv(points, fuv)
!    if (self%periodic(ixyz)) then
!    	fu = wave_function(points(:,self%ixyz), self%u(:), self%box(self%ixyz))
!    else
!    	fu = chebyshev_function(points(:,self%ixyz), self%u(:), self%box(self%ixyz)) & 
!    endif
!    if (self%periodic(jxyz)) then
!	    fv = wave_function(points(:,self%jxyz), self%v(:), self%box(self%jxyz))
!    else
!	    fv = chebyshev_function(points(:,self%jxyz), self%v(:), self%box(self%jxyz))
!    endif
!    fuv = fu * fv
    do j =1, self%n_waves2
		b(j) = b(j) + sum(points(:,self%normal) * fuv(:,j))
	enddo

#if USE_LAPACK
    call multiply_by_tranpose(fuv, fuv, A)
#else
	A = matmul(transpose(fuv), fuv)
#endif

	!Add constraint on minimum size
    A = A + self%diag_matrix

end subroutine get_A_b_mixed

subroutine update_real_surface(self, points)
#if USE_LAPACK
    use lapack_fns, only : solve, multiply_by_tranpose
#else
	use librarymod, only : LUdcmp, lubksb
    use interfaces, only : error_abort
#endif
    implicit none

	class(intrinsic_surface_real) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points

    integer :: i, j
	integer,dimension(:),allocatable		:: indx
    double precision, dimension(:), allocatable :: b
    double precision, dimension(:,:), allocatable :: A

    logical :: debug=.false.
    double precision, dimension(:), allocatable :: surf
    double precision, allocatable, dimension(:) :: x

	!Function to get the matrix A and RHS b
	call self%get_A_b(points, A, b)
	
    !This solves Ax = b
#if USE_LAPACK
    call solve(A, b, x)
#else
	call LUdcmp(A,indx)
	call lubksb(A,indx,b,x)
    ! call error_abort("Error - must build with lapack (p_lapack or p_sys_lapack) to use intrinsic interface")
#endif
	self%coeff(:) = x(:)

    !The Fejér summation (or averaging) to suppress oscillatory Gibbs phenomenon artifacts   
    !do j=1,size(x)
    !	self%coeff(j) = x(j)*(size(x)-j)/size(x)
    !enddo

    !Check surface we just fitted actually matches our points
    if (debug) then
        call self%get_surface(points, surf)
        do j = 1, size(points,1)
            print'(a, i5, 4f10.5, g20.10)', "update_real_surface DEBUG ON Elevation vs points = ", & 
					j, points(j,self%ixyz), points(j,self%jxyz), surf(j) , points(j,self%normal), &
                     surf(j)-points(j,self%normal)
        enddo
        !stop "DEBUG STOP - update_real_surface"
    endif

end subroutine update_real_surface


subroutine get_real_surface(self, points, elevation, include_zeromode, qu)
    implicit none

	class(intrinsic_surface_real) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(out), dimension(:), allocatable :: elevation

	logical, intent(in), optional :: include_zeromode
    integer, intent(in), dimension(2), optional ::  qu

    integer :: j, ui, vi
    integer, dimension(2) :: qu_
    double precision :: zeromode
	
	integer,save :: tcount
    !real(kind(0.d0)) :: t1, t2
    real(kind(0.d0)),save :: timing

	!call cpu_time(t1)

    if (present(qu)) then
        qu_ = qu
    else
        qu_ = self%qm
    endif

    !Get elevation at point from sum of modes
    allocate(elevation(size(points,1)))
    elevation = 0.d0
    do ui = -qu_(1), qu_(1)
    do vi = -qu_(2), qu_(2)
        !j = (2 * self%qm + 1) * (ui + self%qm) + (vi + self%qm) + 1
        j = (2*self%qm(2) + 1) * (ui + self%qm(1)) + (vi + self%qm(2)) + 1
        !print*, "get_real_surface", ui, vi, j, qu_, self%n_waves, self%n_waves2
        elevation = elevation + self%coeff(j) &
                     * wave_function(points(:,self%ixyz), ui, self%box(self%ixyz)) &
                     * wave_function(points(:,self%jxyz), vi, self%box(self%jxyz))
    enddo
    enddo
	if (present(include_zeromode)) then
		if (.not. include_zeromode) then
			call self%get_zero_mode(zeromode)
			elevation = elevation - zeromode
		endif
	endif
    !stop "STOP in get_real_surface"
	!call cpu_time(t2)
	!timing = timing + (t2 - t1)/size(points,1)
	!tcount = tcount + 1
	!if (mod(tcount,100000) .eq. 0) then
	!	print*, "time for 100000 iters of get_real_surface per point", timing
	!	timing = 0.d0
	!	tcount = 0
	!endif

end subroutine get_real_surface


subroutine get_real_surface_binwidth(self, points, elevation, include_zeromode, maprange, qu)
    implicit none

	class(intrinsic_surface_real) :: self

	logical, intent(in), optional :: include_zeromode
    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(in), optional ::  maprange
    integer, intent(in), dimension(2), optional ::  qu
    double precision, intent(out), dimension(:), allocatable :: elevation

    integer :: i, j, ui, vi
    integer, dimension(2) :: qu_
	double precision :: surface_location,  transition
	double precision :: intrinsic_shift_window_bot, intrinsic_shift_window_top
	double precision :: intrinsic_shift_total_bot, intrinsic_shift_total_top

    if (present(qu)) then
        qu_ = qu
    else
        qu_ = self%qm
    endif
	
	if (present(maprange) .or. present(include_zeromode)) then			
		call self%get_zero_mode(surface_location)
	endif

	
	if (present(maprange)) then
		!Loop and include map range
		transition = 3.d0
		call self%get_zero_mode(surface_location)
		intrinsic_shift_window_bot = surface_location-maprange
		intrinsic_shift_window_top = surface_location+maprange
		intrinsic_shift_total_bot = intrinsic_shift_window_bot-transition
		intrinsic_shift_total_top = intrinsic_shift_window_top+transition
	endif
	
	allocate(elevation(size(points,1)))
	elevation = 0.d0
	do i = 1, size(points,1)
	
		!We can exclude points outside of maprange from calculation of the intrinsic surface
		if (present(maprange)) then
			!print*, intrinsic_shift_total_bot, points(i,self%normal), intrinsic_shift_total_top
			if ((points(i,self%normal) .lt. intrinsic_shift_total_bot) .or. & 
			    (points(i,self%normal) .gt. intrinsic_shift_total_top)) then
				elevation(i) = surface_location
				cycle
			endif
		endif
		
		!Get elevation at point from sum of modes
        do ui = -qu_(1), qu_(1)
        do vi = -qu_(2), qu_(2)
            j = (self%qm(1) + self%qm(2) + 1) * (ui + self%qm(1)) + (vi + self%qm(2)) + 1
			elevation(i) = elevation(i) + self%coeff(j) & 
						* wave_function(points(i,self%ixyz), ui, self%box(self%ixyz)) & 
						* wave_function(points(i,self%jxyz), vi, self%box(self%jxyz)) 
		enddo
		enddo

		!Transition range to zero 
		if (present(maprange)) then
			if (points(i,self%normal) .lt. intrinsic_shift_window_bot .and. & 
				points(i,self%normal) .gt. intrinsic_shift_total_bot) then
				elevation(i) = (elevation(i)-surface_location) & 
								*(points(i,self%normal)-intrinsic_shift_total_bot)/transition & 
							   + surface_location
			else if (points(i,self%normal) .gt. intrinsic_shift_window_top .and. & 
					 points(i,self%normal) .le. intrinsic_shift_total_top) then
				elevation(i) = (elevation(i)-surface_location) & 
							   *(intrinsic_shift_window_top-points(i,self%normal))/transition & 
							   + surface_location
			endif
		endif

	enddo

	!Remove zero mode if needed
	if (present(include_zeromode)) then
		if (.not. include_zeromode) then
			elevation(:) = elevation(:) - surface_location
		endif
	endif

end subroutine get_real_surface_binwidth

subroutine get_real_surface_derivative(self, points, dSdr, qu)
    implicit none

	class(intrinsic_surface_real) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    integer, intent(in), dimension(2), optional ::  qu
    double precision, intent(out), dimension(:,:), allocatable :: dSdr

    integer :: j, ui, vi
    integer, dimension(2) :: qu_

    if (present(qu)) then
        qu_ = qu
    else
        qu_ = self%qm
    endif

    allocate(dSdr(size(points,1),2))
    dSdr = 0.d0
    do ui = -qu_(1), qu_(1)
    do vi = -qu_(2), qu_(2)
        j = (2*self%qm(2) + 1) * (ui + self%qm(1)) + (vi + self%qm(2)) + 1
        !j = (2 * self%qm + 1) * (ui + self%qm) + (vi + self%qm) + 1
        dSdr(:,1) = dSdr(:,1) + self%coeff(j) & 
                   * derivative_wave_function(points(:,self%ixyz), ui, & 
                                              self%box(self%ixyz)) & 
                   * wave_function(points(:,self%jxyz), vi, & 
                                   self%box(self%jxyz))
        dSdr(:,2) = dSdr(:,2) + self%coeff(j) & 
                  * wave_function(points(:,self%ixyz), ui, & 
                                  self%box(self%ixyz)) & 
                  * derivative_wave_function(points(:,self%jxyz), vi, & 
                                             self%box(self%jxyz))
    enddo
    enddo

end subroutine get_real_surface_derivative



subroutine get_chebychev_surface(self, points, elevation, include_zeromode, qu)
    implicit none

	class(intrinsic_surface_chebychev) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(out), dimension(:), allocatable :: elevation

	logical, intent(in), optional :: include_zeromode
    double precision, intent(in), optional ::  qu

    integer :: j, ui, vi, qu_
    double precision :: zeromode

    double precision, dimension(:,:), allocatable :: fuv

    if (present(qu)) then
        stop "Error in get_chebychev_surface, optional qu not supported"
    endif

    !Get array of chebychev terms per point
    call self%get_fuv(points, fuv)
	!fuv = chebyshev_function(points(:,self%ixyz), self%u(:), self%box(self%ixyz)) & 
	!	 *chebyshev_function(points(:,self%jxyz), self%v(:), self%box(self%jxyz))

    !Get elevation at point from sum of modes
    allocate(elevation(size(points,1)))
    elevation = 0.d0
	do j=1,self%n_waves2
		elevation(:) = elevation(:) + self%coeff(j)*fuv(:,j)
        !print*, "get_Cheby", j, self%u(j), self%v(j), self%coeff(j), fuv(1,j), elevation(1)
	enddo

	if (present(include_zeromode)) then
		if (.not. include_zeromode) then
            stop "Error in get_chebychev_surface, optional zeromode not supported"
			call self%get_zero_mode(zeromode)
			elevation = elevation - zeromode
		endif
	endif

end subroutine get_chebychev_surface


subroutine get_real_metric_tensor(self, points, g)
    implicit none

	class(intrinsic_surface_real) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(out), dimension(:,:,:), allocatable :: g

	integer :: i
    double precision, dimension(:,:), allocatable :: p, ds 

	allocate(g(size(points,1),2,2))
	allocate(p(1,3))
	do i =1,size(points,1)
		p(1,:) = points(i,:)
		call self%get_surface_derivative(p, ds)
        g(i,1,1) = 1 + ds(1,1)**2
        g(i,1,2) = ds(1,1)*ds(1,2)
        g(i,2,1) = ds(1,2)*ds(1,1)
        g(i,2,2) = 1 + ds(1,2)**2

	enddo

end subroutine get_real_metric_tensor





subroutine get_bilinear_surface(self, points, elevation, include_zeromode, qu)
    implicit none

	class(intrinsic_surface_bilinear) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(out), dimension(:), allocatable :: elevation

	logical, intent(in), optional :: include_zeromode
    integer, dimension(2), intent(in), optional ::  qu

	integer :: i, j, k, n, bins(2)
    double precision, dimension(:), allocatable :: e
    double precision, dimension(:,:), allocatable :: p

	allocate(elevation(size(points,1)))
	allocate(p(1,3))
	do i =1,size(points,1)
		bins = self%get_tangent_bins(points(i,:))
		j = bins(1); k = bins(2)
		p(1,:) = points(i,:)
		call self%get_surface_bilinear(p, self%Abilinear(:,:,j,k), e)
		elevation(i) = e(1)
		!if (elevation(i) .ne. self%intrnsc_smple(j+1,k+1)) then
		!	print*, "Error, intrinsic_smple != surface_bilinear", i,j,k,elevation(i), self%intrnsc_smple(j,k)
		!endif
		!elevation(i) = self%intrnsc_smple(j,k)
	enddo

end subroutine get_bilinear_surface


function intrinsic_area(self) result(A)
    implicit none

	class(intrinsic_surface_real) :: self

    real(kind(0.d0)) :: A

	A = 1.d0 + 0.5d0*pi**2*sum(check_uv(self%u, self%v)*self%coeff**2 & 
                               *(  self%u**2/self%box(self%ixyz)**2 & 
                                 + self%v**2/self%box(self%jxyz)**2))

    A = A * self%area

end function intrinsic_area



subroutine get_bilinear_surface_derivative(self, points, dSdr, qu)
    implicit none

	class(intrinsic_surface_bilinear) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    integer, dimension(2), intent(in), optional ::  qu
    double precision, intent(out), dimension(:,:), allocatable :: dSdr

	integer :: i, j, k, n, bins(2)
    double precision, dimension(:,:), allocatable :: ds 
    double precision, dimension(:,:), allocatable :: p

	allocate(dSdr(size(points,1),2))
	allocate(p(1,3))
	do i =1,size(points,1)
		bins = self%get_tangent_bins(points(i,:))
		j = bins(1); k = bins(2)
		p(1,:) = points(i,:)
		call self%get_surface_derivative_from_bilinear(P, self%Abilinear(:,:,j,k), ds)
		dSdr(i,:) = ds(1,:)
	enddo

end subroutine get_bilinear_surface_derivative


subroutine get_bilinear_metric_tensor(self, points, g)
    implicit none

	class(intrinsic_surface_bilinear) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(out), dimension(:,:,:), allocatable :: g

	integer :: i, j, k, n, bins(2)
    double precision, dimension(:,:), allocatable :: p, ds 

	allocate(g(size(points,1),2,2))
	allocate(p(1,3))
	do i =1,size(points,1)
		p(1,:) = points(i,:)
		call self%get_surface_derivative(p, ds)
        g(i,1,1) = 1 + ds(1,1)**2
        g(i,1,2) = ds(1,1)*ds(1,2)
        g(i,2,1) = ds(1,2)*ds(1,1)
        g(i,2,2) = 1 + ds(1,2)**2

        !g(i,1,1) = 1 !+ ds(1,1)**2
        !g(i,1,2) = 0.d0 !ds(1,1)*ds(1,2)
        !g(i,2,1) = 0.d0 !ds(1,2)*ds(1,1)
        !g(i,2,2) = 1 !+ ds(1,2)**2
	enddo

end subroutine get_bilinear_metric_tensor

!Apply mapping here based on Riemann metric
subroutine apply_metric_tensor_transform(self, A)
    use interfaces, only : error_abort
    implicit none

	class(intrinsic_surface_real) :: self

    real(kind(0.d0)),dimension(:,:,:,:,:),allocatable,intent(inout) :: A

	integer :: i, j, k, t1, t2
    real(kind(0.d0)) :: T11, T12, T21, T22
    real(kind(0.d0)),dimension(:,:), allocatable :: points, T
    real(kind(0.d0)),dimension(:,:,:), allocatable :: g

    if (self%normal .ne. 1) call error_abort("Error apply_metric_tensor_transform - only works for normal in x")
    t1 = self%ixyz; t2 = self%jxyz

!    if (size(A,3) .ne. 2 .or. size(A,4) .ne. 2) then
!        call error_abort("Error apply_metric_tensor_transform - last indices 2 by 2")
    if (size(A,2) .ne. self%nbins(t1) .or. size(A,3) .ne. self%nbins(t2)) then
        print*, "A2=",size(A,2), "nbins t1=", self%nbins(t1), & 
                "A3=",size(A,3), "nbins t2=", self%nbins(t2)
        call error_abort("Error apply_metric_tensor_transform - array A must be same size as bilinear surface")
    endif
    allocate(T(size(A,4),size(A,5)))
    allocate(points(1,2))

    do j = 1, self%nbins(t1)
    do k = 1, self%nbins(t2)

        !Get bin centrepoint and use metric at this point
        points(1,1) = (j-0.5d0)*self%binsize(t1)-0.5d0*self%box(t1)
        points(1,2) = (k-0.5d0)*self%binsize(t2)-0.5d0*self%box(t2)

        call self%get_metric_tensor(points, g)

        do i =1, size(A,1)
    
            !I think this is general with permutes of i,j,k
            if (self%normal .ne. 1) then
                T = A(i,j,k,:,:)
            else if (self%normal .ne. 2) then
                T = A(j,k,i,:,:)
            else 
                T = A(k,i,j,:,:)
            endif

            !Get two tangent components and multiply by metric tensir
            T11 = T(t1,t1)
            T12 = T(t2,t1)
            T21 = T(t1,t2)
            T22 = T(t2,t2)

            T(t1,t1) = T11*g(1,1,1) + T12*g(1,1,2)
            T(t2,t1) = T11*g(1,1,2) + T12*g(1,2,2)
            T(t1,t2) = T21*g(1,1,1) + T22*g(1,2,1)
            T(t2,t2) = T21*g(1,2,1) + T22*g(1,2,2)

            !Store back in array
            if (self%normal .ne. 1) then
                A(i,j,k,:,:) = T
            else if (self%normal .ne. 2) then
                A(j,k,i,:,:) = T
            else
                A(k,i,j,:,:) = T   
            endif

        enddo
!
!        T11 = A(j,k,1,1)
!        T12 = A(j,k,2,1)
!        T21 = A(j,k,1,2)
!        T22 = A(j,k,2,2)

!        A(j,k,1,1) = T11*g(1,1) + T12*g(1,2)
!        A(j,k,2,1) = T11*g(1,2) + T12*g(2,2)
!        A(j,k,1,2) = T21*g(1,1) + T22*g(2,1)
!        A(j,k,2,2) = T21*g(2,1) + T22*g(2,2)         

    enddo
    enddo

end subroutine apply_metric_tensor_transform

!subroutine apply_riemann_transform(A, ISR)
!    use module_record, only : riemann_transform, nhb, nbins
!    implicit none

!	real(kind(0.d0)),dimension(:,:,:,:,:), allocatable :: buf

!    !Apply mapping here based on Riemann metric
!    allocate(buf(nbins(1), nbins(2), nbins(3), size(A,4), size(A,5)))
!    buf = A(1+nhb(1):nbins(1)+nhb(1), & 
!            1+nhb(2):nbins(2)+nhb(2), &
!            1+nhb(3):nbins(3)+nhb(3),:,:)
!    call ISR%apply_metric_tensor_transform(buf)
!    A(1+nhb(1):nbins(1)+nhb(1), & 
!      1+nhb(2):nbins(2)+nhb(2), &
!      1+nhb(3):nbins(3)+nhb(3),:,:) = buf


!end subroutine apply_riemann_transform


    !Apply mapping here based on Riemann metric
    !if (riemann_transform .eq. 1) then
    !     call ISR_mdt%apply_metric_tensor_transform(momentum_flux_))!(1+nhb(1):nbins(1)+nhb(1), & 
                                                                   !1+nhb(2):nbins(2)+nhb(2), &
                                                                   !1+nhb(3):nbins(3)+nhb(3),:,:))
    !endif


!    subroutine metric_tensor_transform(self, A)

!	    real(kind(0.d0)),dimension(:,:,:,:),allocatable,intent(inout) :: A

!        !THIS ISR_mdt%get_metric FUNCTION WILL LOOK SOMETHING LIKE
!        t1 = self%ixyz; t2 = self%jxyz
!        do j = 1, nbins(t1)
!        do k = 1, nbins(t2)

!            points(1,1) = float(j-0.5)*binsize(t1)-0.5d0*self%box(t1)
!            points(1,2) = float(k-0.5)*binsize(t2)-0.5d0*self%box(t2)

!            call ISR_mdt%get_metric_tensor(points, g)

!            !These are stress in the surface, assuming x is normal
!            if (t1 .eq. 2 .and. t2 .eq. 3) then

!                !do i=1,nbins(1)
!                !momentum_flux_(i,j,k,2:3,2:3) = matmul(momentum_flux_(i,j,k,2:3,2:3), g(1,:,:))

!                T11 = A(j,k,2,2)
!                T12 = A(j,k,3,2)
!                T21 = A(j,k,2,3)
!                T22 = A(j,k,3,3)

!                A(j,k,2,2) = T11*g(1,1) + T12*g(1,2)
!                A(j,k,3,2) = T11*g(1,2) + T12*g(2,2)
!                A(j,k,2,3) = T21*g(1,1) + T22*g(2,1)
!                A(j,k,3,3) = T21*g(2,1) + T22*g(2,2)         

!            endif
!    end subroutine metric_tensor_transform


subroutine get_zero_mode(self, zeromode)
    implicit none

	class(intrinsic_surface_real) :: self

    integer :: j
    real(kind(0.d0)) :: zeromode

    j = (self%qm(1) + self%qm(2) + 1) * self%qm(1) + self%qm(2) + 1
    zeromode = self%coeff(j)

end subroutine get_zero_mode


!A version getting the explicit location from the intrinsic surface
function get_bin_from_surface(self, r, nbins, nhb) result(bin)
    implicit none

	class(intrinsic_surface_real) :: self

    real(kind(0.d0)),intent(in),dimension(3) :: r
    integer,dimension(3),intent(in),optional     :: nbins, nhb

    integer,dimension(3)	                 :: bin

    integer :: n, i, j
    integer,dimension(3) 	:: nbins_, nhb_
    real(kind(0.d0)) :: maprange, zeromode
    real(kind(0.d0)), dimension(3) :: halfdomain, binsize
    real(kind(0.d0)), allocatable, dimension(:) :: elevation, test
    real(kind(0.d0)), allocatable, dimension(:,:) :: points

	if (present(nbins) .and. present(nhb)) then
		binsize = self%box/float(nbins)
		nhb_ = nhb
		nbins_ = nbins
	else
		binsize = self%binsize
		nhb_ = self%nhb
		nbins_ = self%nbins
	endif
    halfdomain = 0.5*self%box
    n=self%normal
    i=self%ixyz
    j=self%jxyz

    allocate(points(1,3))
    points(1,:) = r(:)

    call self%get_surface(points, elevation, include_zeromode=.true.)

    !Added a shift by zero wavelength so surface is not at zero
    bin(n) = ceiling((r(n)+halfdomain(n)-elevation(1)+0.5d0*binsize(n))/binsize(n))+nhb_(n) !HALF SHIFT
    bin(i) = ceiling((r(i)+halfdomain(i))/binsize(i))+nhb_(i)
    bin(j) = ceiling((r(j)+halfdomain(j))/binsize(j))+nhb_(j)

	if (bin(n) > nbins_(n)+nhb_(n)) then
        bin(n) = nbins_(n)+nhb_(n)
    elseif (bin(n) < 1 ) then
        !select type (self)
        !type is (intrinsic_surface_chebychev)
        !    call get_chebychev_surface(self, points, test)
        !end select
        !"print*, "UNDERFLOW BIN", r, elevation, test, bin
        bin(n) = 1
    endif

end function get_bin_from_surface

!A version getting the explicit location from the intrinsic surface
function get_tangent_bins(self, r, nbins, nhb) result(bin)
    implicit none

	class(intrinsic_surface_real) :: self

    real(kind(0.d0)),intent(in),dimension(3) :: r
    integer,dimension(3),intent(in),optional    :: nbins, nhb

    integer,dimension(2)	                 :: bin

    integer :: n, i, j
    integer,dimension(3) 	:: nbins_, nhb_
    real(kind(0.d0)), dimension(3) :: halfdomain, binsize
	
	if (present(nbins) .and. present(nhb)) then
		binsize = self%box/float(nbins)
		nhb_ = nhb
		nbins_ = nbins
	else
		binsize = self%binsize
		nhb_ = self%nhb
		nbins_ = self%nbins
	endif
    
    halfdomain = 0.5*self%box
    i=self%ixyz
    j=self%jxyz
	
    bin(1) = ceiling((r(i)+halfdomain(i))/binsize(i))+nhb_(i)
    bin(2) = ceiling((r(j)+halfdomain(j))/binsize(j))+nhb_(j)

end function get_tangent_bins

function index_to_vertex(self, j, k) result(vertex)
    implicit none

	integer, intent(in) :: j, k
	double precision, dimension(4) :: vertex

	class(intrinsic_surface_real) :: self

	!vertex = (ysb, zsb, yst, zst)
	vertex(1) = float(j-1-self%nhb(self%ixyz))*self%binsize(self%ixyz)-0.5d0*self%box(self%ixyz)
	vertex(2) = float(k-1-self%nhb(self%jxyz))*self%binsize(self%jxyz)-0.5d0*self%box(self%jxyz)
	vertex(3) = float( j -self%nhb(self%ixyz))*self%binsize(self%ixyz)-0.5d0*self%box(self%ixyz)
	vertex(4) = float( k -self%nhb(self%jxyz))*self%binsize(self%jxyz)-0.5d0*self%box(self%jxyz)

end function index_to_vertex

function indices_to_points(self, i, j, k) result(points)
    implicit none

	integer, intent(in) :: i, j, k
	double precision, dimension(4,3) :: points

	class(intrinsic_surface_bilinear) :: self

	integer :: n, t1, t2
	double precision :: shift	
	double precision, dimension(4) :: vert	

	n = self%normal
	t1 = self%ixyz
	t2 = self%jxyz

	!Position in normal direction
	shift = (i-1*self%nhb(n)-0.5d0)*self%binsize(n)-0.5d0*self%box(n) 
	points(1,n) = self%intrnsc_smple( j ,  k ) + shift
	points(2,n) = self%intrnsc_smple(j+1,  k ) + shift
	points(3,n) = self%intrnsc_smple(j  , k+1) + shift
	points(4,n) = self%intrnsc_smple(j+1, k+1) + shift

	!position in tangential
	vert = self%index_to_vertex(j,k)
	points(1,t1) = vert(1); points(1,t2) = vert(2) !Bottom left
	points(2,t1) = vert(3); points(2,t2) = vert(2) !Bottom right
	points(3,t1) = vert(1); points(3,t2) = vert(4) !Top left
	points(4,t1) = vert(3); points(4,t2) = vert(4) !Top right	

end function indices_to_points


subroutine get_crossings_bilinear(self, r1, r2, bin1, bin2, n, rc, crossings, cbins)
	use bilnear_intersect, only : line_plane_intersect, line_patch_intersect
    use librarymod, only : Qsort, extend_array2d, extend_array1d_int
	implicit none

	class(intrinsic_surface_bilinear) :: self

	integer, intent(in)                      :: n   !normal
	integer, intent(in), dimension(3) 	     :: bin1, bin2
	real(kind(0.d0)),intent(in),dimension(3) :: r1, r2
	
	logical, intent(out)                    :: crossings
	integer,intent(out), optional, dimension(:), allocatable :: cbins
	real(kind(0.d0)),intent(out),dimension(:,:), allocatable :: rc

	integer                 :: t1, t2, t(2), i, j, k, m, c, mp1, jb, kb, ixyz, ind
	integer                 :: maxbin, minbin, minTbin, maxTbin
	real(kind(0.d0))        :: pt, r12(3)

	integer :: flag, ss, Ns, bc, tempsize, normbin, norm1, norm2
	integer, dimension(2) :: tbin
	integer, dimension(3) :: bin, bin_mdt, cbin, db, bt, bb
	integer, dimension(:), allocatable :: cbinstemp, indices
	integer, dimension(:,:), allocatable :: binstemp

	real(kind(0.d0)),save :: maxerror
	real(kind(0.d0)),parameter		:: eps = 1e-12
	real(kind(0.d0)) :: vert(4), s, dx, dy, dz, ds, shift
	real(kind(0.d0)) :: y(2,2), z(2,2), P(2,2,3), A(2,2)
 	real(kind(0.d0)), dimension(3) :: rb, rs, re, bs1, bs2, rc1, rc2
	real(kind(0.d0)), dimension(:), allocatable :: unordered_cross
	real(kind(0.d0)),dimension(3,2)  :: intersect, normal, intersect2
	real(kind(0.d0)),dimension(:,:), allocatable  :: points, temp,  cross

	logical, save :: first_bincount=.true.
	integer :: ibin, jbin, kbin, maxcross
	integer,dimension(:,:,:), allocatable, save  :: bincount

	if (all(bin1 .eq. bin2)) then
		crossings = .false.
		return
	endif

	t1=self%ixyz
	t2=self%jxyz
	r12   = r2 - r1		        !Molecule i trajectory between t-dt and t
	where (abs(r12) .lt. 0.000001d0) r12 = 0.000001d0

	!Allocate array for all normal crossings
	db = abs(bin1 - bin2)
	if (db(t1) .eq. 0) db(t1) = 1
	if (db(t2) .eq. 0) db(t2) = 1
	!Add two for two end points
	allocate(cross(db(t1)+db(t2)+2, 5))

	!Store start at ri
	cross(1, 1 ) = 0.d0
	cross(1,2:4) = r1
	cross(1, 5 ) = 0.d0 !normal crossing set to zero for end points

	!Get all crossings on flat surfaces
	m = 2
	t = (/ t1, t2 /)
	do i=1,2
		call self%get_flat_crossings(r1, r2, bin1, bin2, t(i), rc, crossings)
		if (crossings) then
			do j=1,size(rc,2)
				!Get s value of each crossing from rc(:) = r1(:) + s*r12(:)
				s = (rc(1,j)-r1(1))/r12(1)
				
				!Need to store which direction of surface it crosses...
				cross(m, 1 ) = s
				cross(m,2:4) = rc(:,j)
				cross(m, 5 ) = dble(t(i))
				m = m + 1
			enddo
		endif
	enddo

	!Store start at ri
	cross(m, 1 ) = 1.d0
	cross(m,2:4) = r2
	cross(m, 5 ) = 0.d0 !normal crossing set to zero for end points

	!Order crossings in increasng s
	allocate(unordered_cross(m))
	allocate(indices(m))
    do ind = 1, size(indices,1)
        indices(ind) = ind
    enddo
	unordered_cross = cross(:m,1)
	call Qsort(unordered_cross, indices)

	!Go through crossings in order
	bc = 0
	!Pre allocate arrays to large enough for possible crossings
	maxcross = 16
	allocate(temp(3,maxcross))
	if (present(cbins)) allocate(cbinstemp(maxcross))
	do i=1,size(indices,1)-1
		m = indices(i)
		mp1 = indices(i+1)

		!Retrieve and define values
		bs1 = (/0.,0.,0./)
		bs2 = (/0.,0.,0./)
		rc1 = cross(m  ,2:4)
		rc2 = cross(mp1,2:4)
		norm1 = cross(m  ,5) !normal from 1 to 3
		norm2 = cross(mp1,5) !normal from 1 to 3

		!Setup arrays to shift above/below plane
		!unless it's the particles (end points) with normal=0
		! Need sign r12 to specify direction of increasing s
		if (norm1 .ne. 0) bs1(norm1) = sign(1.d0,r12(norm1))
		if (norm2 .ne. 0) bs2(norm2) = sign(1.d0,r12(norm2))

		!Get bins either side of each crossing (\pm eps*normal of each surface)
		bt = self%get_bin(rc1+eps*bs1) 
		bb = self%get_bin(rc2-eps*bs2)

        !print*, "xbins=", bt(n)-bb(n), bt(n), bb(n), rc1(n), rc2(n), r1, r2

		!As we've got tangential crossings, bins here should only change 
		!in normal direction here so must be same in t1 and t2
		!if ((bt(t1) .ne. bb(t1)) .or. (bt(t2) .ne. bb(t2))) then
            !It is possible y and z change at the same time
            !if (abs(bt(t1)-bb(t1)) .ne. abs(bt(t2)-bb(t2))) then
			!    print*, "Warning in cumulative_flux_opt - bins not same in y and z", bt, bb, bs1, bs2
            !endif
			!stop "Error in cumulative_flux_opt bins, should be same in y and z"
		!endif

		!Loop over all x cells between  points
		!This needs to be increasing or decreasing as needed
		do normbin = bb(n), bt(n)-1, sign(1,bt(n)-bb(n))
			!Take lower surfaces of bins if moving in the 
			!negative direction
			if (sign(1,bt(n)-bb(n)) .eq. -1) then		
				j = normbin-1
			else
				j = normbin
			endif
			
			points = self%indices_to_points(j, bt(t1), bt(t2))

			!Get a patch of P values to use in bilinear patch & line calculation
			P(:,:,1) = reshape(points(:,n ), (/2,2/))
			P(:,:,2) = reshape(points(:,t1), (/2,2/))
			P(:,:,3) = reshape(points(:,t2), (/2,2/))

			flag = 0
			call line_patch_intersect(r1, r12, P, intersect, normal, flag)

			do c=1,flag

				if ((intersect(1,c) .gt. rc1(1) .and. &
					 intersect(1,c) .lt. rc2(1)) .or. &
					(intersect(1,c) .lt. rc1(1) .and. &
					 intersect(1,c) .gt. rc2(1))) then
					bc = bc + 1
					!Not trivial/cheap to get a priori size of temp 
					!so better to check and extend if needed
					if (bc .gt. maxcross) then
						call extend_array2d(temp)
						if (present(cbins)) call extend_array1d_int(cbinstemp)
						maxcross = size(temp,2)
					endif
					temp(:,bc) = intersect(:,c)
					if (present(cbins)) cbinstemp(bc) = j
				endif

			enddo

		enddo

	enddo

	!if (bb(n) .ne. bt(n) .and. bc .lt. 1) then
		!do i=1,size(indices,1)
		!	m = indices(i)
		!	print*, "Ordered crossings", i, m, cross(m,:)
		!enddo
		!print*, "Missed crossing", bt, bb, bin1, bin2, flag, intersect, r1, r2, points
		!write(345,*) 1, r1, r2, P
		!stop "Error in get_crossings_bilinear - Missed crossing"
	!endif

	!Copy array of intersects to rc
	if (bc .ge. 1) then
		if (allocated(rc)) deallocate(rc)
		allocate(rc(3, bc))
		rc = temp(:,1:bc)
		if (present(cbins)) then
			allocate(cbins(bc))
			cbins = cbinstemp(1:bc)
		endif
		crossings = .true.
	else
		if (allocated(rc)) deallocate(rc)
		allocate(rc(3, 1)); rc = 0.d0
		if (present(cbins)) then
			allocate(cbins(1))
			cbins = 1
		endif
		crossings = .false.
		!stop "Error in get_crossing_bilinear - No crossings found"
	endif

end subroutine get_crossings_bilinear



subroutine get_crossings_bilinear_old(self, r1, r2, bin1, bin2, n, rc, crossings, cbins)
	use bilnear_intersect, only : line_plane_intersect, line_patch_intersect
	implicit none

	class(intrinsic_surface_bilinear) :: self

	integer, intent(in)                      :: n   !normal
	integer, intent(in), dimension(3) 	     :: bin1, bin2
	real(kind(0.d0)),intent(in),dimension(3) :: r1, r2
	
	logical, intent(out)                    :: crossings
	integer,intent(out), optional, dimension(:), allocatable :: cbins
	real(kind(0.d0)),intent(out),dimension(:,:), allocatable :: rc

	integer                 :: t1, t2, i, j, k, jb, kb, ixyz
	integer                 :: maxbin, minbin, minTbin, maxTbin
	real(kind(0.d0))        :: pt, r12(3)

	integer :: flag, ss, Ns, cross, bc, tempsize
	integer, dimension(2) :: tbin
	integer, dimension(3) :: bin, bin_mdt, cbin
	integer, dimension(:), allocatable :: cbinstemp
	integer, dimension(:,:), allocatable :: binstemp

	real(kind(0.d0)) :: vert(4), s, dx, dy, dz, ds, shift
	real(kind(0.d0)),save :: maxerror
	real(kind(0.d0)) :: y(2,2), z(2,2), P(2,2,3), A(2,2)
 	real(kind(0.d0)), dimension(3) :: rb, rs, re
	real(kind(0.d0)), dimension(:), allocatable :: elevation
	real(kind(0.d0)),dimension(3,2)  :: intersect, normal, intersect2
	real(kind(0.d0)),dimension(:,:), allocatable  :: points, points2, temp

	!Tangents
	t1 = self%ixyz !mod(n,3)+1
	t2 = self%jxyz !mod(n+1,3)+1
	minbin = min(bin1(n), bin2(n))
	maxbin = max(bin1(n), bin2(n))
	r12   = r2 - r1  !Molecule i trajectory between t-dt and t
	where (abs(r12) .lt. 0.000001d0) r12 = 0.000001d0

	!Factor of two as can be maximum of two per surface
	tempsize = 1
	do i=1,3
		minTbin = min(bin1(i), bin2(i))
		maxTbin = max(bin1(i), bin2(i))
		tempsize = tempsize * (maxTbin-minTbin+1)
	enddo
	allocate(points(4,3))
	allocate(temp(3, 2*tempsize))
	if (present(cbins)) allocate(cbinstemp(2*tempsize))

	bc = 0

	do i=minbin-2, maxbin-1
	do j=min(bin1(t1), bin2(t1)),max(bin1(t1), bin2(t1))
	do k=min(bin1(t2), bin2(t2)),max(bin1(t2), bin2(t2))

		points = self%indices_to_points(i,j,k)

		!Get a patch of P values to use in bilinear patch & line calculation
		P(:,:,1) = reshape(points(:,n ), (/2,2/))
		P(:,:,2) = reshape(points(:,t1), (/2,2/))
		P(:,:,3) = reshape(points(:,t2), (/2,2/))

		call line_patch_intersect(r1, r12, P, intersect, normal, flag)


		do ixyz=1,flag
				bc = bc + 1
				temp(:,bc) = intersect(:,ixyz)
				if (present(cbins)) cbinstemp(bc) = i
		enddo
	enddo
	enddo
	enddo

	!Copy array of intersects to rc
	if (bc .ge. 1) then
		allocate(rc(3, bc))
		rc = temp(:,1:bc)
		if (present(cbins)) then
			allocate(cbins(bc))
			cbins = cbinstemp(1:bc)
		endif
		crossings = .true.
	else
		allocate(rc(3, 1)); rc = 0.d0
		if (present(cbins)) then
			allocate(cbins(1))
			cbins = 1
		endif
		crossings = .false.
		!stop "Error in get_crossing_bilinear - No crossings found"
	endif

end subroutine get_crossings_bilinear_old





subroutine get_crossings_real(self, r1, r2, bin1, bin2, n, rc, crossings, cbins)
	use bilnear_intersect, only : line_plane_intersect
	implicit none

	class(intrinsic_surface_real) :: self

	integer, intent(in)                      :: n   !normal
	integer, intent(in), dimension(3) 	     :: bin1, bin2
	real(kind(0.d0)),intent(in),dimension(3) :: r1, r2
	
	logical, intent(out)                    :: crossings
	integer,intent(out), optional, dimension(:), allocatable :: cbins
	real(kind(0.d0)),intent(out),dimension(:,:), allocatable :: rc

	integer                 :: t1, t2, i, j, k, jb, kb, ixyz
	integer                 :: maxbin, minbin, minTbin, maxTbin
	real(kind(0.d0))        :: pt, r12(3)

	integer :: flag, ss, Ns, cross, bc, tempsize
	integer, dimension(3) :: bin, bin_mdt
	integer, dimension(:), allocatable :: cbinstemp

	real(kind(0.d0)) :: yrb, yrt, zrb, zrt, s, ds
	real(kind(0.d0)),save :: maxerror
	real(kind(0.d0)) :: y(2,2), z(2,2), P(2,2,3), A(2,2), rb(3)
	real(kind(0.d0)), dimension(:), allocatable :: elevation
	real(kind(0.d0)),dimension(3,2)  :: intersect, normal, intersect2
	real(kind(0.d0)),dimension(:,:), allocatable  :: points, points2, temp

    select type (self)
    type is (intrinsic_surface_chebychev)
        stop "Error in get_crossings, using intrinsic_surface_chebychev in real routine"
    end select

	if (bin1(n) .eq. bin2(n)) then
		crossings = .false.
	else
		crossings = .true.

		!Tangents
		if (self%normal .ne. n) stop "Error in get_crossings_bilinear, normal not consistent with surface"
		if (self%ixyz .ne. mod(n,3)+1) stop "Error in get_crossings_bilinear, normal not consistent with surface"
		if (self%jxyz .ne. mod(n+1,3)+1) stop "Error in get_crossings_bilinear, normal not consistent with surface"
		t1 = self%ixyz !mod(n,3)+1
		t2 = self%jxyz !mod(n+1,3)+1

		minbin = min(bin1(n), bin2(n))
		maxbin = max(bin1(n), bin2(n))

		r12   = r2 - r1  !Molecule i trajectory between t-dt and t
		where (abs(r12) .lt. 0.000001d0) r12 = 0.000001d0

		!Get positions of molecules which
		!describe a bounding box
		yrb = min(r1(t1), r2(t1))
		yrt = max(r1(t1), r2(t1))
		zrb = min(r1(t2), r2(t2))
		zrt = max(r1(t2), r2(t2))

		allocate(points(4,3))
		allocate(points2(4,3))
		!Factor of two as can be maximum of two per surface
		tempsize = 1
		do i=1,3
			minTbin = min(bin1(i), bin2(i))
			maxTbin = max(bin1(i), bin2(i))
			tempsize = tempsize * (maxTbin-minTbin+1)
		enddo
		allocate(temp(3, 2*tempsize))
		if (present(cbins)) allocate(cbinstemp(2*tempsize))

		! First sample at r1 
		Ns = maxbin-minbin
		ds = 1.d0 / real(Ns, kind(0.d0))
		s = ds
		
		!loop over bin range
		bc = 1
		do i = minbin, maxbin-1

			! Loop over all samples, s varies from 0 to 1
			! Position of sample on line
			rb(:) = r1(:) + s*r12(:)
			yrb = min(r1(t1), rb(t1))
			yrt = max(r1(t1), rb(t1))
			zrb = min(r1(t2), rb(t2))
			zrt = max(r1(t2), rb(t2))
			!print'(a,3i5,10f10.5)', "patch", i, bc, Ns, s, r1(:), rb(:), r2(:)
			!print'(a,3i5,8f10.5)', "patch", i, bc, Ns, s, rb(:), yrb, yrt, zrb, zrt
			s = s + ds

			!Use distance between molecules
			!to create a bilinear patch
			points(:,n) = rb(n) !Not used
			points(1,t1) = yrb; points(1,t2) = zrb
			points(2,t1) = yrb; points(2,t2) = zrt 
			points(3,t1) = yrt; points(3,t2) = zrb 
			points(4,t1) = yrt; points(4,t2) = zrt

			!Get surface in x as a function of y and z
			call self%get_surface(points, elevation, include_zeromode=.true.)

			!Normal direction used in line plane intersect
			!Shift by current bin
			points(:,n) = elevation(:) & 
						  + (i-1*self%nhb(n)-0.5d0)*self%binsize(n) &
						  -0.5d0*self%box(n) 

			!Get a patch of P values to use in bilinear patch & line calculation
			P(:,:,1) = reshape(points(:,n ), (/2,2/))
			P(:,:,2) = reshape(points(:,t1), (/2,2/))
			P(:,:,3) = reshape(points(:,t2), (/2,2/))

			!print'(a,2i5,18f10.5)', "x_bilinear", i,j, p(1,1,:), P(1,2,:), P(2,1,:), p(2,2,:), r1(:), r12(:)

			!Special case of flat surface causes problems so need to handle separatly
			if (P(1,1,1) .eq. P(2,1,1) .and. &
				P(2,1,1) .eq. P(1,2,1) .and. &
				P(1,2,1) .eq. P(2,2,1)) then
				intersect = -666
				intersect(1,1) = P(1,1,1)
				intersect(2,1) = r1(t1)+(r12(t1)/r12(n))*(intersect(1,1)-r1(n))            
				intersect(3,1) = r1(t2)+(r12(t2)/r12(n))*(intersect(1,1)-r1(n))
				flag = 1
			else
				call line_plane_intersect(r1, r12, P, intersect, normal, flag)
			endif
			
			!Loop over intersects and add to temp
			do ixyz=1,size(intersect,2)
				if (all(abs(intersect(:,ixyz)+666) .gt. 1e-7)) then
					!if (size(intersect,2) .ne. 1) print'(a,3i5,9f10.5)', "intrsect", ixyz, i, bc, r1, intersect(:,ixyz), r2
					temp(:,bc) = intersect(:,ixyz)
					if (present(cbins)) cbinstemp(bc) = i
					bc = bc + 1
				endif
			enddo

		enddo

		!Copy array of intersects to rc
		if (bc .gt. 1) then
			allocate(rc(3, bc-1))
			rc = temp(:,1:bc-1)
			if (present(cbins)) then
				allocate(cbins(bc-1))
				cbins = cbinstemp(1:bc-1)
			endif
		else
			print*, "Warning in get_crossings_bilinear - no crossings found ", minbin, maxbin, r1, r2
			! Crossing of bilinear differs from Fourier surface
			! Walk line between two points and try to get location of crossing
			Ns = maxbin-minbin+1
			if (allocated(points)) deallocate(points)
			allocate(points(Ns+2,3))
			ds = 1.d0 / real(Ns, kind(0.d0))
			! First sample at r1 
			s = -ds
			! Loop over all samples, s varies from 0 to 1
			bin_mdt = bin1
			do ss = 1, Ns+2
				! Position of sample on line
				points(ss,:) = r1(:) + s*r12(:)
				bin = self%get_bin(points(ss,:), self%nbins, self%nhb)
				!print*, "Points on line", s, points(ss,:), bin, bin_mdt, &
				!bin(self%normal) .ne. bin_mdt(self%normal), self%normal
				!If bin changes then must be a crossing
				cross = bin(self%normal) - bin_mdt(self%normal) 
				if (cross .eq. 1) then
					allocate(rc(3, 1))
					rc(:, 1) = r1(:) + (s-0.5d0*ds)*r12(:)
					!print*, "Crossing found", rc, bin, s-0.5d0*ds
					if (present(cbins)) then
						allocate(cbins(1))
						if (cross .gt. 0) then
							cbins(1) = bin_mdt(self%normal)
						else 
							cbins(1) = bin(self%normal)
						endif
					endif
					exit
				endif
				bin_mdt = bin
				s = s + ds
			end do	
			!call self%get_surface(points, elevation)
			!print*, "Elevations", elevation+ (bin1(n)-1*nhb(n))*binsize(n)-0.5d0*self%box(self%normal)
			!stop "Error get_crossings_bilinear - interactions must be missed"
		endif
	endif

end subroutine get_crossings_real


subroutine get_flat_crossings(self, r1, r2, bin1, bin2, n, rc, crossings)
	implicit none

	class(intrinsic_surface_real) :: self

	integer, intent(in)                      :: n   !normal
	integer, intent(in), dimension(3) 	     :: bin1, bin2
	real(kind(0.d0)),intent(in),dimension(3) :: r1, r2
	
	logical, intent(out)                    :: crossings
	real(kind(0.d0)),intent(out),dimension(:,:), allocatable :: rc

	integer                 :: t1, t2, i, j, maxbin, minbin
	double precision        :: pt, r12(3), halfdomain(3)

    halfdomain = 0.5*self%box

	if (bin1(n) .eq. bin2(n)) then
		crossings = .false.
	else
		crossings = .true.

		!Tangents
		t1 = mod(n,3)+1
		t2 = mod(n+1,3)+1
		minbin = min(bin1(n), bin2(n))
		maxbin = max(bin1(n), bin2(n))

		r12   = r1 - r2							!Molecule i trajectory between t-dt and t
		where (abs(r12) .lt. 0.000001d0) r12 = 0.000001d0

		allocate(rc(3, maxbin-minbin))
		j = 1
		do i = minbin, maxbin-1
			pt = (i-1*self%nhb(n))*self%binsize(n)-halfdomain(n)
			rc(n,j) = pt
			rc(t1,j) = r1(t1)+(r12(t1)/r12(n))*(pt-r1(n))            
			rc(t2,j) = r1(t2)+(r12(t2)/r12(n))*(pt-r1(n))
			j = j + 1
		enddo

	endif

end subroutine get_flat_crossings

!subroutine get_real_surface_derivative(self, points, dir, dSdr, qu)
!    implicit none

!	class(intrinsic_surface_real) :: self

!    integer, intent(in) :: dir
!    double precision, intent(in), dimension(:,:), allocatable ::  points
!    double precision, intent(in), optional ::  qu
!    double precision, intent(out), dimension(:), allocatable :: dSdr

!    integer :: j, ui, vi, qu_

!    if (present(qu)) then
!        qu_ = qu
!    else
!        qu_ = self%qm
!    endif

!    allocate(dSdr(size(points,1)))
!    dSdr = 0.d0
!    do ui = -qu_, qu_
!    do vi = -qu_, qu_
!        j = (2 * self%qm + 1) * (ui + self%qm) + (vi + self%qm) + 1

!        if (dir .eq. self%ixyz) then
!            dSdr(:) = dSdr(:) + self%coeff(j) & 
!                      * derivative_wave_function(points(:,self%ixyz), ui, self%box(self%ixyz)) & 
!                      * wave_function(points(:,self%jxyz), vi, self%box(self%jxyz)) 
!        elseif (dir .eq. self%jxyz) then
!            dSdr(:) = dSdr(:) + self%coeff(j) & 
!                      * wave_function(points(:,self%ixyz), ui, self%box(self%ixyz)) & 
!                      * derivative_wave_function(points(:,self%jxyz), vi, self%box(self%jxyz))
!        else 
!            stop "Error in get_surface_derivative - normal direction requested"
!        endif
!    enddo
!    enddo

!end subroutine get_real_surface_derivative

!=========================================
!=========================================
!   Adapted from pytim surface code
!  Fits complex form of Fourier series
!=========================================
!=========================================

subroutine compute_q_vectors(self, box, normal, alpha, eps, topbot)
    use librarymod, only : meshgrid2d
    implicit none

	class(intrinsic_surface_complex) :: self

    integer, intent(in)         :: normal, topbot
    double precision, intent(in) :: alpha, eps
    double precision, intent(in), dimension(3) :: box

    integer :: i,j
    real(kind(0.d0)), dimension(:), allocatable :: qx, qy
    real(kind(0.d0)), dimension(:,:), allocatable :: qx_vectors,  qy_vectors

    !Compute the q-vectors compatible with the current box dimensions.

    !Inputs:
    !box       : List of domain dimensions [Lx, Ly, Lz]
    !normal    : Normal to surface x=1, y=2 and z=3
    !alpha     : Molecular scale cutoff, default 2.0
    !eps       : Constraint tolerance

    self%normal = normal
    self%box = box
    self%alpha = alpha
    self%eps = eps
    self%topbot = topbot
    !save in object:
    !q_vectors : two 2D arrays forming the grid of q-values, similar
    !            to a meshgrid
    !mode_shape: Number of modes
    !Qxy       : array of the different q-vectors
    !Q         : squared module of Qxy with the first element missing

    !Shift to allow normal that is any direction
    self%ixyz = mod(normal,3)+1
    self%jxyz = mod(normal+1,3)+1
    self%modes_shape(1) = ceiling(box(self%ixyz) / alpha)
    self%modes_shape(2) = ceiling(box(self%jxyz) / alpha)
    allocate(qx(self%modes_shape(1)), qy(self%modes_shape(2)))
    do i=1,self%modes_shape(1)
        qx(i) = dble((i-1)) * 2.d0 * pi / box(1)
    enddo
    do i=1,self%modes_shape(2)
        qy(i) = dble((i-1)) * 2.d0 * pi / box(2)
    enddo

    call meshgrid2D(qx, qy, qx_vectors, qy_vectors)
    allocate(self%q_vectors(2, self%modes_shape(1), self%modes_shape(2)))
    self%q_vectors(1,:,:) = qx_vectors
    self%q_vectors(2,:,:) = qy_vectors

    allocate(self%Qxy(self%modes_shape(1)*self%modes_shape(2),2))
    do i = 1, self%modes_shape(1)
    do j = 1, self%modes_shape(2)
        self%Qxy(j+(i-1)*self%modes_shape(2),1) = self%q_vectors(1,i,j)
        self%Qxy(j+(i-1)*self%modes_shape(2),2) = self%q_vectors(2,i,j)
    enddo
    enddo
    !Exclude zeroth mode
    allocate(self%Q(self%modes_shape(1)*self%modes_shape(2)-1))
    self%Q = sqrt(self%Qxy(2:,1)**2 + self%Qxy(2:,2)**2)


end subroutine compute_q_vectors

subroutine update_surface_modes(self, points)
#if USE_LAPACK
    use lapack_fns, only : pinverse
#else
    use interfaces, only : error_abort
#endif

    implicit none

	class(intrinsic_surface_complex) :: self

    double precision, intent(in), dimension(:,:), allocatable :: points

    logical :: debug=.false.
    integer :: i,j
    real(kind(0.d0)), parameter :: pi=3.141592653589793d0
    real(kind(0.d0)) :: az, Lx, Ly
    real(kind(0.d0)), dimension(:), allocatable :: z, surf
    real(kind(0.d0)), dimension(:,:), allocatable :: xy, QR
    double complex, dimension(:), allocatable :: s
    double complex, dimension(:,:), allocatable :: A, ph, pinv_ph, modesT
    double precision, dimension(:,:,:), allocatable :: q_vectors

    !Make some definitions
    allocate(xy(size(points,1),2))
    xy(:,1) = points(:,self%ixyz);  xy(:,2) = points(:,self%jxyz)
    allocate(z(size(points,1)))
    z = points(:,self%normal)
    az = sum(z)/size(z)
    z = z - az

    ! Qx is an unravelled array of modes, we dot each of the x-y parts with a point
    ! to get the contribution of each point to spectral basis set
    allocate(QR(size(xy,1), size(self%Qxy,1)))
    do i=1,size(QR,1)
    do j=1,size(QR,2)
        QR(i,j) = xy(i,1)*self%Qxy(j,1) + xy(i,2)*self%Qxy(j,2)
        !print'(2i6,5f10.5)', i,j,QR(i,j), xy(i,1), Qxy(j,1),xy(i,2), Qxy(j,2)
    enddo
    enddo

    ! We exclude the zero mode and add the mean back in as zero mode instead
    allocate(ph(size(QR,1), size(QR,2)-1))
    do j=1,size(ph,2)
        ph(:,j) = dcmplx(cos(QR(:,j+1)), sin(QR(:,j+1)))/ self%Q(j) 
    enddo

    ! Constraint here Least Squares solution
    Lx = 2.d0 * pi / self%Qxy(self%modes_shape(1)+1,2)    
    Ly = 2.d0 * pi / self%Qxy(2,1)  
    do j=1,size(ph,2)
        ph(:,j) = ph(:,j) + self%eps * dcmplx(self%Q(j), 0.d0)
    enddo

#if USE_LAPACK
    ! Least square solution solving ph*z = s
    call pinverse(ph, pinv_ph)
#else
    call error_abort("Error - must build with lapack (p_lapack or p_sys_lapack) to use intrinsic interface")
#endif

    !Multiply inverse with z values to get mode coefficients
    allocate(s(self%modes_shape(1)*self%modes_shape(2)))
    s(1) = dcmplx(az, 0.d0)
    do i=2,size(s,1)
        !s(i) = dot_product(conjg(pinv_ph(i-1,:)),z(:))/Q(i-1)
        s(i) = sum(pinv_ph(i-1,:)*z(:))/self%Q(i-1)
    enddo

    !Return modes in right shape (for some reason this needs to be done this way)
    modesT = reshape(s, (/self%modes_shape(2), self%modes_shape(1)/))
    if (.not. allocated(self%modes)) then
        allocate(self%modes(self%modes_shape(1),self%modes_shape(2)))
    elseif ((size(self%modes,1) .ne. self%modes_shape(1)) .or. &
            (size(self%modes,2) .ne. self%modes_shape(2))) then
        deallocate(self%modes)
        allocate(self%modes(self%modes_shape(1),self%modes_shape(2)))
    endif
    self%modes = transpose(modesT)

    !Check surface we just fitted actually matches our points
    if (debug) then
        call get_surface_complex(self, points, surf)
        do i = 1, size(points,1)
            print*, "Elevation vs points = ", i, surf(i) , points(i,self%normal), &
                     surf(i)-points(i,self%normal),   surf(i)/points(i,self%normal)
        enddo
        print'(a,10f10.5)', "Modes", self%modes(1,1), self%modes(1,2), & 
                    self%modes(2,1), self%modes(1,3), self%modes(3,1)
        !stop "Stop in debug of update_surface_modes" 
    endif

end subroutine update_surface_modes



subroutine get_surface_complex(self, points, elevation)
    implicit none

	class(intrinsic_surface_complex) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(out), dimension(:), allocatable :: elevation

    integer :: i, j, n

    double precision, dimension(:,:), allocatable :: dotp
    double complex, dimension(:,:), allocatable :: phase

    allocate(elevation(size(points,1)))
    allocate(dotp(size(self%q_vectors,2),size(self%q_vectors,3)))
    allocate(phase(size(self%q_vectors,2),size(self%q_vectors,3)))

    elevation = 0.d0
    do n = 1,size(points,1)
		do i = 1,size(self%q_vectors,2)
		do j = 1,size(self%q_vectors,3)
			!dotp(i,j) = dot_product(q_vectors(:,i,j), points(n, :))
			dotp(i,j) = self%q_vectors(1,i,j) * points(n,self%ixyz) & 
					  + self%q_vectors(2,i,j) * points(n,self%jxyz)
			phase(i,j) = dcmplx(cos(dotp(i,j)), sin(dotp(i,j)))
			!print'(3i7,8f10.5)', i,j,n,points(n, ixyz),points(n,jxyz),modes(i,j), & 
			!                     phase(i,j), real(phase(i,j) * modes(i,j)), elevation(n)
			elevation(n) = elevation(n) + real(phase(i,j) * self%modes(i,j))
		enddo
		enddo
    enddo

end subroutine get_surface_complex



subroutine sample_intrinsic_surface(self, vertices, nbins, writeiter)
    use librarymod, only : get_new_fileunit, get_Timestep_FileName
    implicit none

	class(intrinsic_surface_real) :: self

    integer,intent(in),dimension(3), optional :: nbins

    real(kind(0.d0)), dimension(:,:,:),intent(out), allocatable :: vertices
    integer, intent(in), optional :: writeiter

    integer,dimension(3) :: nbins_
    real(kind(0.d0)),dimension(3) :: binsize

    logical          :: debug=.false., writeobj=.false.
    logical, save    :: first_time = .true.
    integer          :: i, j, k, n, v, ixyz, jxyz, fileno
    real(kind(0.d0)) :: vert(4), area
    real(kind(0.d0)), dimension(:), allocatable :: elevation
    real(kind(0.d0)), dimension(:,:), allocatable :: points
    character(200) :: outfile_t, filename

	if (present(nbins)) then
		binsize = self%box/float(nbins)
        nbins_ = nbins
	else
		binsize = self%binsize
		nbins_ = self%nbins
	endif

    if (present(writeiter)) then      
        writeobj = .true.
    else
        writeobj = .false.
    endif

    allocate(points(4,3))
    allocate(vertices(nbins_(self%ixyz)*nbins_(self%jxyz),4,3))

    if (writeobj) then
        fileno = get_new_fileunit()
    	if (present(nbins)) then
            call get_Timestep_FileName(writeiter,"./results/surface",outfile_t)
        else
            call get_Timestep_FileName(writeiter,"./results/bilinear",outfile_t)
        endif
        write(filename,'(a,a4)') trim(outfile_t),'.obj'
        open(fileno, file=trim(filename))
    endif

    n = 1; v = 1
    do j = 1,nbins_(self%ixyz)
    do k = 1,nbins_(self%jxyz)

		!vert = self%index_to_vertex(j,k)
        vert(1) = float(j-1)*binsize(self%ixyz)-0.5d0*self%box(self%ixyz)
        vert(3) = float(j  )*binsize(self%ixyz)-0.5d0*self%box(self%ixyz)
        vert(2) = float(k-1)*binsize(self%jxyz)-0.5d0*self%box(self%jxyz)
        vert(4) = float(k  )*binsize(self%jxyz)-0.5d0*self%box(self%jxyz)

		points(1,self%ixyz) = vert(1); points(1,self%jxyz) = vert(2) !Bottom left
		points(2,self%ixyz) = vert(3); points(2,self%jxyz) = vert(2) !Bottom right
		points(3,self%ixyz) = vert(1); points(3,self%jxyz) = vert(4) !Top left
		points(4,self%ixyz) = vert(3); points(4,self%jxyz) = vert(4) !Top right
        !points(1,self%ixyz) = ysb; points(1,self%jxyz) = zsb !Bottom left
        !points(2,self%ixyz) = yst; points(2,self%jxyz) = zsb !Bottom right
        !points(3,self%ixyz) = ysb; points(3,self%jxyz) = zst !Top left
        !points(4,self%ixyz) = yst; points(4,self%jxyz) = zst !Top right

        select type (self)
        type is (intrinsic_surface_chebychev)
            call get_chebychev_surface(self, points, elevation)
        class is (intrinsic_surface_real)
            call get_real_surface(self, points, elevation)
        end select
        !call surface_from_modes(points, 3, q_vectors, modes, elevation)
        do i = 1, 4
            vertices(v,i,:) = (/elevation(i), points(i,self%ixyz), points(i,self%jxyz)/)
            !print*, i, j, k, v, vertices(v,i,:)
        enddo
       
        !Write to wavefunction obj format
        if (writeobj) then
            !Needs to be clockwise!
            do i =1, 4
                write(fileno,'(a1, 3f15.8)') "v", vertices(v,i,:)
            enddo
            write(fileno,'(a1, 4i8)') "f", n, n+1, n+3, n+2
            n = n + 4
        endif
        v = v + 1

    enddo
    enddo

    if (writeobj) close(fileno)

end subroutine sample_intrinsic_surface

!Update the saved array of samples of the surface to use for efficient 
!calculation of crossings and bin values


subroutine update_sampled_surface(self, nbins, nhb)
    implicit none

	class(intrinsic_surface_bilinear) :: self

    integer,intent(in),dimension(3), optional :: nbins, nhb

    logical          :: debug=.false., writeobj=.false.
    logical, save    :: first_time = .true.
    integer          :: j, k
    real(kind(0.d0)) :: v(4), binsize(3), A(2,2)
    real(kind(0.d0)), dimension(:), allocatable :: elevation, elevationtest
    real(kind(0.d0)), dimension(:,:), allocatable :: points

    integer,dimension(3) 	:: nbins_, nhb_
 	
	if (present(nbins) .and. present(nhb)) then
		binsize = self%box/float(nbins)
		nbins_ = nbins
		nhb_ = nhb
	else
		binsize = self%binsize
		nhb_ = self%nhb
		nbins_ = self%nbins
	endif
	
	allocate(points(4,3))

    do j = 1,self%nbins(self%ixyz)+2*self%nhb(self%ixyz)
    do k = 1,self%nbins(self%jxyz)+2*self%nhb(self%jxyz)

		v = self%index_to_vertex(j,k)
		points(1,self%ixyz) = v(1); points(1,self%jxyz) = v(2) !Bottom left
		points(2,self%ixyz) = v(3); points(2,self%jxyz) = v(2) !Bottom right
		points(3,self%ixyz) = v(1); points(3,self%jxyz) = v(4) !Top left
		points(4,self%ixyz) = v(3); points(4,self%jxyz) = v(4) !Top right

        select type (self)
        type is (intrinsic_surface_chebychev)
            call get_chebychev_surface(self, points, elevation, include_zeromode=.true.)
        class is (intrinsic_surface_bilinear)
            call get_real_surface(self, points, elevation, include_zeromode=.true.)
        end select

		points(:,self%normal) = elevation
        self%intrnsc_smple( j ,  k ) = elevation(1)
        self%intrnsc_smple(j+1,  k ) = elevation(2)
        self%intrnsc_smple(j  , k+1) = elevation(3)
        self%intrnsc_smple(j+1, k+1) = elevation(4)

		call self%paramterise_bilinear_surface(points, A)
		self%Abilinear(:,:,j,k) = A(:,:)

		!DEBUG
		!call self%get_surface_bilinear(points, self%Abilinear(:,:,j,k), elevationtest)
		!if (elevationtest(1) .ne. self%intrnsc_smple(j,k)) then
		!	print'(a,2i6,5f10.5)', "intrnsc_smple != surface_bilinear", j,k,elevationtest, self%intrnsc_smple(j,k)
		!endif
		!call self%get_surface_bilinear(points, A, elevationtest)
        !print*, elevation - elevationtest
		!print'(a,2i4,16f10.5)', "sampled_surf", j, k, points(:,self%normal)-elevationtest, & 
		!			points(:,self%ixyz), points(:,self%jxyz), elevation
		
    enddo
    enddo
	
end subroutine update_sampled_surface



subroutine get_sampled_surface(self, r, points)
    implicit none

	class(intrinsic_surface_real) :: self

    double precision,intent(in),dimension(3) :: r
    real(kind(0.d0)), dimension(:,:), intent(out), allocatable :: points

    integer          :: j, k, cells(2)
    real(kind(0.d0)) :: v(4)
    real(kind(0.d0)), dimension(:), allocatable :: elevation

	cells = self%get_tangent_bins(r)
	j = cells(1); k = cells(2)

    allocate(points(4,3))
	points(1, self%normal)=self%intrnsc_smple(j  ,k  )
	points(2, self%normal)=self%intrnsc_smple(j+1,k  )
	points(3, self%normal)=self%intrnsc_smple(j  ,k+1)
	points(4, self%normal)=self%intrnsc_smple(j+1,k+1)

	v = self%index_to_vertex(j,k)
	points(1,self%ixyz) = v(1); points(1,self%jxyz) = v(2) !Bottom left
	points(2,self%ixyz) = v(3); points(2,self%jxyz) = v(2) !Bottom right
	points(3,self%ixyz) = v(1); points(3,self%jxyz) = v(4) !Top left
	points(4,self%ixyz) = v(3); points(4,self%jxyz) = v(4) !Top right

	!ysb = float(j-1+self%nhb(self%ixyz))*self%binsize(self%ixyz)-0.5d0*self%box(self%ixyz)
	!yst = float(j  +self%nhb(self%ixyz))*self%binsize(self%ixyz)-0.5d0*self%box(self%ixyz)
	!zsb = float(k-1+self%nhb(self%jxyz))*self%binsize(self%jxyz)-0.5d0*self%box(self%jxyz)
	!zst = float(k  +self%nhb(self%jxyz))*self%binsize(self%jxyz)-0.5d0*self%box(self%jxyz)
	! points(1,self%ixyz) = ysb; points(1,self%jxyz) = zsb !Bottom left
	! points(2,self%ixyz) = yst; points(2,self%jxyz) = zsb !Bottom right
	! points(3,self%ixyz) = ysb; points(3,self%jxyz) = zst !Top left
	! points(4,self%ixyz) = yst; points(4,self%jxyz) = zst !Top right
	
	!print'(a,3f10.5,2i5,8f10.5)', "get_sampled_surface", r,j,k,ysb,yst,zsb,zst, points(:, self%normal)

end subroutine get_sampled_surface



subroutine write_modes(self, writeiter)
    use librarymod, only : get_new_fileunit, get_Timestep_FileName
    implicit none

	class(intrinsic_surface_real) :: self

    integer, intent(in) :: writeiter


    integer 		:: fileno
    integer 		:: j, ui, vi
    character(200)  :: outfile_t

	fileno = get_new_fileunit() 
	call get_Timestep_FileName(writeiter,"./results/surfacemodes",outfile_t)
	open(fileno, file=trim(outfile_t))

    do j=1,self%n_waves2
	    ui = self%u(j) 
        vi = self%v(j)
        write(fileno, '(3i12,f25.18)') ui, vi, j, self%coeff(j)
    enddo

!    do ui = -self%qm, self%qm
!    do vi = -self%qm, self%qm
!        j = (2 * self%qm + 1) * (ui + self%qm) + (vi + self%qm) + 1
!        write(fileno, '(3i12,f25.18)') ui, vi, j, self%coeff(j)
!    enddo
!    enddo

    close(fileno)

end subroutine write_modes


subroutine get_initial_pivots(points, ISR, pivots, ntarget)
    use librarymod, only : Qsort
    implicit none

	class(intrinsic_surface_real), intent(in) :: ISR	! declare an instance
    double precision, intent(in), dimension(:,:), allocatable ::  points
    integer, intent(in) ::  ntarget
    integer, intent(out), dimension(:), allocatable ::  pivots


    integer :: Npivots, maxpivots, i, ind, ratio
    integer, dimension(2) :: nxy, bins
    integer, dimension(:,:), allocatable :: sectors

    integer, dimension(:), allocatable :: indices
    double precision, dimension(3) :: binsize
    double precision, dimension(:), allocatable :: z

    ! Defines the initial pivots as a set of bins^2 particles, where
    ! each particle is in a distinct sector formed by dividing
    ! the macroscopic plane into bins(1)xbins(2) regions.
    ! In original work, value of 3 by 3 used (but small domains)
    ! For large domains or non flat domains this should be larger
    ! so initial guess is well fitted to general shape 
    ! Taking 50% of total target fittings
    bins = int(sqrt(ISR%box(ISR%jxyz)*0.5d0*ntarget/ISR%box(ISR%ixyz)))

    !Take into account non-equal domains
	ratio = int(bins(1)*ISR%box(ISR%ixyz)/ISR%box(ISR%jxyz))
	if (ratio .gt. 0) then
		bins(2) = ratio
	else
		bins(2) = 1
	endif

    !print*, "get_initial_pivots", bins, ntarget
    allocate(sectors(bins(1), bins(2)))
    sectors = 0
    Npivots = bins(1)*bins(2)
    maxpivots = Npivots
    allocate(pivots(Npivots))
    binsize(1) = ISR%box(ISR%ixyz)/dble(bins(1))
    binsize(2) = ISR%box(ISR%jxyz)/dble(bins(2))

    !Setup indices
    allocate(indices(size(points,1))) 
    do ind = 1, size(indices,1)
        indices(ind) = ind
    enddo

    !Sort ascending
    allocate(z(size(points,1)))
    z = points(:,ISR%normal)
    call Qsort(z, indices)

    if (ISR%topbot .eq. 1) then
        do ind = size(z,1), 1, -1
            nxy(1) = ceiling((points(indices(ind), ISR%ixyz) + 0.5d0*ISR%box(ISR%ixyz)) / binsize(1))
            nxy(2) = ceiling((points(indices(ind), ISR%jxyz) + 0.5d0*ISR%box(ISR%jxyz)) / binsize(2))

            if (nxy(1) .lt. 1 .or. nxy(2) .lt. 1 .or. & 
                nxy(1) .gt. bins(1) .or. nxy(2) .gt. bins(2)) then
                !print*, ind, nxy, points(indices(ind), :)
                cycle
            endif
    !        nxy(1) = floor(2.99999999999999999 * points(indices(ind), ixyz) / box(ixyz))+1
    !        nxy(2) = floor(2.99999999999999999 * points(indices(ind), jxyz) / box(jxyz))+1
            if (sectors(nxy(1), nxy(2)) .eq. 0) then
                pivots(Npivots) = indices(ind)
                sectors(nxy(1), nxy(2)) = 1
                Npivots = Npivots - 1
                !print*, ind, Npivots, pivots(Npivots+1), points(indices(ind), ixyz), points(indices(ind), jxyz)
            endif
            if (sum(sectors) .ge. maxpivots) then
                exit
            endif
        enddo

    elseif (ISR%topbot .eq. 2) then
        do ind = 1, size(z,1)
            nxy(1) = ceiling((points(indices(ind), ISR%ixyz) + 0.5d0*ISR%box(ISR%ixyz)) / binsize(1))
            nxy(2) = ceiling((points(indices(ind), ISR%jxyz) + 0.5d0*ISR%box(ISR%jxyz)) / binsize(2))

            if (nxy(1) .lt. 1 .or. nxy(2) .lt. 1 .or. & 
                nxy(1) .gt. bins(1) .or. nxy(2) .gt. bins(2)) then
                cycle
            endif
            if (sectors(nxy(1), nxy(2)) .eq. 0) then
                pivots(Npivots) = indices(ind)
                sectors(nxy(1), nxy(2)) = 1
                Npivots = Npivots - 1
                print*, "Bottom pivots", ind, points(indices(ind), :) 
            endif
            if (sum(sectors) .ge. maxpivots) then
                exit
            endif
        enddo
    endif

    if (Npivots .ne. 0) stop "Not all initial pivots found"

end subroutine get_initial_pivots



!subroutine update_pivots(points, ISR, pivots, tau, new_pivots)
!    use librarymod, only : Qsort
!    implicit none

!	class(intrinsic_surface_real), intent(in) :: ISR	! declare an instance
!    integer, intent(in), dimension(:), allocatable ::  pivots
!    double precision, intent(in) :: tau
!    double precision, intent(in), dimension(:,:), allocatable ::  points
!    integer, intent(out), dimension(:), allocatable ::  new_pivots

!    logical :: found_range
!    integer :: i, n, ixyz, jxyz, ind, nPivots, sp, qm, qu
!    double precision :: z_max, z_min, area
!    integer, dimension(:), allocatable :: candidates, updated_pivots, indices
!    double precision, dimension(:), allocatable :: z, surf
!    double precision, dimension(:,:), allocatable :: candidates_pos, pivot_pos

!    ! Searches for points within a distance tau from the
!    ! interface.
!    sp = size(pivots,1)
!    allocate(pivot_pos(sp,3)) 
!    do i =1, sp
!        pivot_pos(i,:) = points(pivots(i),:)
!    enddo
!    z_max = maxval(pivot_pos(:,ISR%normal)) + ISR%alpha * 2.d0
!    z_min = minval(pivot_pos(:,ISR%normal)) - ISR%alpha * 2.d0

!    z = points(:, ISR%normal)
!    allocate(indices(size(points,1))) 
!    do ind = 1, size(indices,1)
!        indices(ind) = ind
!    enddo
!    call Qsort(z, indices)
!    n = 0; found_range=.false.
!    allocate(candidates(size(points,1))) 
!    if (ISR%topbot .eq. 1) then
!        do i = size(z,1), 1, -1
!            if ((z(i) > z_min) .and. (z(i) < z_max)) then
!                n = n + 1
!                candidates(n) = indices(i)
!                found_range = .true.
!                !print*, "values", i, n, z_min, z(i), z_max, candidates(n)
!            else if (found_range) then
!                !If z is sorted and we've been inside the range,
!                !once we leave, no point looping anymore
!                exit
!            endif

!        enddo
!    elseif (ISR%topbot .eq. 2) then
!        do i = 1, size(z,1)
!            if ((z(i) > z_min) .and. (z(i) < z_max)) then
!                n = n + 1
!                candidates(n) = indices(i)
!                found_range = .true.
!                !print*, "values", i, n, z_min, z(i), z_max, candidates(n)
!            else if (found_range) then
!                !If z is sorted and we've been inside the range,
!                !once we leave, no point looping anymore
!                exit
!            endif

!        enddo
!    endif 

!    !Get positions from indices
!    allocate(candidates_pos(n,3))
!    do i =1, n
!        !print*, "values", i, n, candidates(i), points(candidates(i),:)
!        candidates_pos(i,:) = points(candidates(i),:)
!    enddo

!    !Recalculate surface at candidate pivot locations
!    select type (ISR)
!    type is (intrinsic_surface_chebychev)
!        call get_chebychev_surface(ISR, candidates_pos, surf)
!    class is (intrinsic_surface_real)
!        call get_real_surface(ISR, candidates_pos, surf)
!    end select

!    nPivots = 0
!    allocate(updated_pivots(n))
!    do i =1,n
!        if ((surf(i)-candidates_pos(i,ISR%normal))**2 .lt. tau**2) then
!            nPivots = nPivots + 1
!            updated_pivots(nPivots) = candidates(i)
!        endif
!    enddo
!    allocate(new_pivots(nPivots))
!    new_pivots = updated_pivots(1:nPivots)

!    !allocate(new_pivots(nPivots+sp))
!    !new_pivots(1:sp) = pivots(:)
!    !new_pivots(sp+1:) = updated_pivots(:nPivots)

!    !print*, new_pivots 

!end subroutine update_pivots


!An alternative way of updating pivots using a surface which 
!gets atoms a distance tau from existing surface
subroutine update_pivots_alt(points, ISR, pivots, tau, ns, new_pivots)
    use librarymod, only : Qsort
    use interfaces, only : error_abort
    implicit none

	class(intrinsic_surface_real), intent(in) :: ISR	! declare an instance
    integer, intent(in), dimension(:), allocatable ::  pivots
    double precision, intent(in) :: tau, ns
    double precision, intent(in), dimension(:,:), allocatable ::  points
    integer, intent(out), dimension(:), allocatable ::  new_pivots

    logical :: found_range
    integer :: i, n, ixyz, jxyz, ind, nmols, nPivots, sp, qm, qu
    double precision :: z_max, z_min, area
    integer, dimension(:), allocatable :: candidates, updated_pivots, indices
    double precision, dimension(:), allocatable :: z, surf
    double precision, dimension(:,:), allocatable :: candidates_pos, pivot_pos

    !Get surface for all molecular locations
    nmols = int(ns*ISR%area)

    select type (ISR)
    type is (intrinsic_surface_chebychev)
        call get_chebychev_surface(ISR, points, surf)
    class is (intrinsic_surface_real)
        call get_real_surface(ISR, points, surf)
    end select
    !call ISR%get_surface(points, surf)

    allocate(z(size(points,1)))
    if (ISR%topbot .eq. 1) then
        z = points(:,ISR%normal)-surf(:)
    else if (ISR%topbot .eq. 2) then
        z = surf(:)-points(:,ISR%normal)
        call error_abort("Only tested for top/right surface")
    endif 
    allocate(indices(size(points,1)))
    do ind = 1, size(indices,1)
        indices(ind) = ind
    enddo
    call Qsort(z, indices)
    n = 0; found_range=.false.
    allocate(candidates(size(points,1))) 
    do i = size(z,1), 1, -1
        !Test and add positive molecules first so preferentially
        !include all outer molecules before moving into cluster
        !Take five times fitting range to ensure we get all outer
        !molecules
        if (z(i) .lt. 5*tau .and. z(i) .gt. 0.d0) then
            n = n + 1
            candidates(n) = indices(i)
            !print*, "update_pivots_alt +ve range", i, n, z(i), tau, candidates(n)
        else if (abs(z(i)) .lt. tau) then
            n = n + 1
            candidates(n) = indices(i)
            found_range = .true.
            !print*, "update_pivots_alt -ve range", i, n, z(i), tau, candidates(n)
        else if (found_range) then
            !If z is sorted and we've been inside the range,
            !once we leave, no point looping anymore
            !print*, "=============================="
            !print*, "update_pivots_alt left range"
            !print*, "=============================="
            exit
        endif
        if (n .ge. nmols) exit
    enddo
    allocate(new_pivots(n))
    new_pivots = candidates(1:n)

!    do i = 1,size(pivots,1)
!        print*, i, pivots(i)
!    enddo
    !Get surface for all molecular locations
!    call ISR%get_surface(points, surf)
!    n = 0
!    allocate(candidates(size(points,1))) 
!    do i = 1, size(points,1)
!        if (abs(points(i,ISR%normal)-surf(i)) < tau) then
!            n = n + 1
!            candidates(n) = i !indices(i)
!        endif
!    enddo
!    new_pivots = candidates(1:n)


end subroutine update_pivots_alt



subroutine fit_intrinsic_surface_modes(self, points, tau, ns, pivots)
    !DEBUG DEBUGDEBUGDEBUG
    !use physical_constants_MD, only : np
    !use calculated_properties_MD, only : nbins
    !use computational_constants_MD, only : nhb
    !DEBUG DEBUGDEBUGDEBUG

    implicit none

	class(intrinsic_surface_real) :: self	! declare an instance

    double precision, intent(in) :: tau, ns
    double precision, intent(in), dimension(:,:), allocatable ::  points
    integer, dimension(:), allocatable, intent(inout) :: pivots

    integer :: i, j, ixyz, jxyz, sp, ntarget, try, maxtry=100
    integer, dimension(2) :: modes_shape
    integer, dimension(:), allocatable :: indices, new_pivots, initial_pivots
    !integer, dimension(:), allocatable, save :: pivots
    double precision, save :: savetau=0.d0
    double precision :: Error, tau_, rand, diff, maxpivot
    double precision, dimension(:), allocatable :: Q, z, d, surf
    double precision, dimension(:,:), allocatable :: Qxy, pivot_pos, initial_pivot_pos

    !Define things
    tau_ = max(tau, savetau)
    ntarget = int(ns*self%area)

    !Get initial pivots
    call get_initial_pivots(points, self, initial_pivots, ntarget)
    if (allocated(pivots)) deallocate(pivots)
    pivots = initial_pivots

    !Get positions and fits surface to initial pivots
    sp = size(pivots,1)
    allocate(pivot_pos(sp,3))
    do i =1, sp
        pivot_pos(i,:) = points(pivots(i),:)
        !print*, "initial pivot pos = ", i, pivot_pos(i,:)
    enddo

    call update_real_surface(self, pivot_pos)

!    do i=1,sp
!        print*, "Pivot error = ", pivot_pos(i,:), surf(i), pivot_pos(i,3)-surf(i)
!    enddo

    do try = 1, maxtry

        !Get new pivots
        !call update_pivots(points, self, pivots, tau_, new_pivots)
        call update_pivots_alt(points, self, pivots, tau_, ns, new_pivots)

        !Get new positions and new modes
        if (allocated(pivot_pos)) deallocate(pivot_pos)
        sp = size(new_pivots,1)
        allocate(pivot_pos(sp,3))
        do i =1, sp
            pivot_pos(i,:) = points(new_pivots(i),:)
            !print*, "new pivot pos = ", i, pivot_pos(i,:)
        enddo
        !call get_surface_modes(pivot_pos, Qxy, modes_shape, Q, modes, eps)
        call update_real_surface(self, pivot_pos)

        !DEBUG DEBUG DEBUG DEBUG
        !select type (self)
        !type is (intrinsic_surface_chebychev)
        !    call get_chebychev_surface(self, pivot_pos, surf)
        !class is (intrinsic_surface_real)
        !    call get_real_surface(self, pivot_pos, surf)
        !end select

        !do i =1, sp
        !    print*, "fi_surface_modes Difference = ", try, i, new_pivots(i), pivot_pos(i, self%normal),&
        !                 surf(i), pivot_pos(i, self%normal) - surf(i)
        !enddo
        !DEBUG DEBUG DEBUG DEBUG

        !Exit once we have converged to particles on surface / area = ns
        !print*, size(new_pivots)/area, size(new_pivots), area, ns 
        !if (size(new_pivots)/self%area .gt. ns)  then
        if (size(new_pivots) .ge. ntarget)  then
            deallocate(pivots)
            ! Truncate pivots to give ns as all positions in 
            ! order of increasing distance from max pivot mol
            allocate(pivots(ntarget))
            pivots = new_pivots(1:ntarget)

            !DEBUG DEBUG DEBUG DEBUG
            !Debugging print statement of surface molecules
!            allocate(initial_pivot_pos(9,3))
!            do i =1, 9
!                initial_pivot_pos(i,:) = points(initial_pivots(i),:)
!            enddo
!            maxpivot = maxval(initial_pivot_pos(:,self%normal))
!            do i =1, ntarget
!                diff = maxpivot - pivot_pos(i, self%normal)
!                print'(a,2i5,3f10.5,f18.8)', "new pivot pos = ", i, ntarget, pivot_pos(i,:), diff
!            enddo
            !DEBUG DEBUG DEBUG DEBUG

            exit
        !If stuck on same numbers of pivots, try increasing search range
        elseif (size(new_pivots,1) .eq. size(pivots,1)) then
            print'(a, f10.5, a, i5, a, i6, a, i6, a)', "Increasing Tau to ", tau_, " at try ",  try, & 
                   " with ", size(new_pivots), " pivots of target ", ntarget , " found."
            tau_ = tau_ + 0.1*tau
            savetau = tau_
        else
            deallocate(pivots)
            allocate(pivots(size(new_pivots,1)))
            pivots = new_pivots
            deallocate(new_pivots)
        endif

        !Some error handling here
        if (try .eq. 50) then

            !Plot updated surface error
            if (allocated(d)) deallocate(d)
            allocate(d(sp))

            select type (self)
            type is (intrinsic_surface_chebychev)
                call get_chebychev_surface(self, pivot_pos, surf)
            class is (intrinsic_surface_real)
                call get_real_surface(self, pivot_pos, surf)
            end select

            d = pivot_pos(:, self%normal) - surf(:)
            Error = sqrt(sum(d * d) / size(d))

            print*, "Try no. = ", try, "No. pivots = ", size(pivots), & 
                    "Error=", Error, "Area=", self%area, size(pivots)/self%area

            if (Error .gt. 1e2) then
                print*, "Solution appears to be diverging, reverting to initial pivots"
                deallocate(pivots)
                call get_initial_pivots(points, self, pivots, ntarget)
                deallocate(pivot_pos)
                sp = size(pivots,1)
                allocate(pivot_pos(sp,3))
                do i =1, sp
                    pivot_pos(i,:) = points(pivots(i),:)
                    !print*, "initial pivot pos = ", i, pivot_pos(i,:)
                enddo
                call update_real_surface(self, pivot_pos)

            endif
        endif

    enddo

	!This was handled by defining a constructor which called the parent
	!constructor but then type is parent
	!select type (self)
	!class is (intrinsic_surface_bilinear)
		!Update the saved surface at locations to use for quick surface 
		!interaction calculations 
	!	call self%update_sampled_surface(self%nbins, self%nhb)
	!end select

    !Plot updated surface error
    !DEBUG DEBUG DEBUG DEBUG
!    if (allocated(d)) deallocate(d)
!    allocate(d(sp))
!    call self%get_surface(pivot_pos, surf)

!    d = pivot_pos(:, self%normal) - surf(:)
!    Error = sqrt(sum(d * d) / size(d))
!    print*, "fit_intrinsic_surface_modes", try, Error, ntarget
    !DEBUG DEBUG DEBUG DEBUG

end subroutine fit_intrinsic_surface_modes

!Save coefficiients to matrix for bilinear 
! x - coordinates of corners x = [x1, x2]
! y - coordinate of corners y = [y1, y1]
! P - elevation at 4 corners P = [P11, P12, P21, P22]
!
!    2,1______2,2
!       |     |
!       |_____|
!    1,1      1,2
!    
! A - coefficient matrix for bilinear form 
! f = A0 + A1*x + A2*y + A3*x*y

subroutine get_bilinear_surface_coeff(x, y, P, A)
    implicit none

    real(kind(0.d0)),intent(in)  :: x(2), y(2)
    real(kind(0.d0)),intent(in)  :: P(2,2)
    real(kind(0.d0)),intent(out) :: A(2,2)
    real(kind(0.d0)) :: x1, x2, y1, y2, x12, y12, d
    x1 = x(1); x2 = x(2)
    y1 = y(1); y2 = y(2)
	!Determinate
    x12 = x2 - x1; y12 = y2 - y1
    d = x12*y12
    A(1,1) = ( P(1,1)*x2*y2 - P(1,2)*x2*y1 - P(2,1)*x1*y2 + P(2,2)*x1*y1)/d
    A(2,1) = (-P(1,1)*y2    + P(1,2)*y1    + P(2,1)*y2    - P(2,2)*y1   )/d
    A(1,2) = (-P(1,1)*x2    + P(1,2)*x2    + P(2,1)*x1    - P(2,2)*x1   )/d
    A(2,2) = ( P(1,1)       - P(1,2)       - P(2,1)       + P(2,2)      )/d

end subroutine get_bilinear_surface_coeff

subroutine paramterise_bilinear_surface(self, P, A)
    implicit none

	class(intrinsic_surface_bilinear), intent(in) :: self	! declare an instance

    real(kind(0.d0)),intent(in)  :: P(4,3)
    real(kind(0.d0)),intent(out) :: A(2,2)
    real(kind(0.d0)) :: x1, x2, y1, y2, x12, y12, d
    x1 = P(1,self%ixyz); x2 = P(4,self%ixyz)
    y1 = P(1,self%jxyz); y2 = P(4,self%jxyz)
	!Determinate
    x12 = x2 - x1; y12 = y2 - y1
    d = x12*y12
    A(1,1) = (  P(1,self%normal)*x2*y2 - P(3,self%normal)*x2*y1 &
			  - P(2,self%normal)*x1*y2 + P(4,self%normal)*x1*y1)/d
    A(2,1) = (- P(1,self%normal)*y2    + P(3,self%normal)*y1 &
			  + P(2,self%normal)*y2    - P(4,self%normal)*y1   )/d
    A(1,2) = (- P(1,self%normal)*x2    + P(3,self%normal)*x2 &
			  + P(2,self%normal)*x1    - P(4,self%normal)*x1   )/d
    A(2,2) = (  P(1,self%normal)       - P(3,self%normal) &
			  - P(2,self%normal)       + P(4,self%normal)      )/d

end subroutine paramterise_bilinear_surface


!Get surface position from positions for a known bin
subroutine get_surface_bilinear(self, points, A, elevation)
    implicit none

	class(intrinsic_surface_bilinear), intent(in) :: self	! declare an instance

    real(kind(0.d0)), intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(in), dimension(2,2) :: A
    double precision, intent(out), dimension(:), allocatable :: elevation

    integer :: n!, ixyz, jxyz

    allocate(elevation(size(points,1)))
    elevation = 0.d0
    do n=1,size(points,1)
        elevation(n) =  A(1,1) & 
					  + A(2,1)*points(n,self%ixyz) & 
					  + A(1,2)*points(n,self%jxyz) &
					  + A(2,2)*points(n,self%ixyz)*points(n,self%jxyz)
	enddo

    ! do n=1,size(points,1)
		! do ixyz=1,2
		! do jxyz=1,2
			! elevation(n) = elevation(n) + A(ixyz, jxyz) & 
						   ! *(points(n,self%ixyz)**(ixyz-1))*(points(n,self%jxyz)**(jxyz-1))
		! enddo
		! enddo
    ! enddo

end subroutine get_surface_bilinear


!Get surface position from positions for a known bin
subroutine get_surface_derivative_from_bilinear(self, points, A, dSdr)
    implicit none

	class(intrinsic_surface_bilinear), intent(in) :: self	! declare an instance

    real(kind(0.d0)), intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(in), dimension(2,2) :: A
    double precision, intent(out), dimension(:,:), allocatable :: dSdr

    integer :: n!, ixyz, jxyz

    allocate(dSdr(size(points,1),2))
    dSdr = 0.d0
    do n=1,size(points,1)
        dSdr(n,1) =  A(2,1)+ A(2,2)*points(n,self%jxyz)
		dSdr(n,2) =  A(1,2)+ A(2,2)*points(n,self%ixyz)
	enddo

end subroutine get_surface_derivative_from_bilinear



!Get surface area estimate by summing all bilinear patches
function intrinsic_area_bilinear(self) result(Area)
    implicit none

	class(intrinsic_surface_bilinear), intent(in) :: self

	integer :: i, j, k
    double precision	 :: Area, patchArea, dx, dy

    dx = self%box(self%ixyz)/self%nbins(self%ixyz)
    dy = self%box(self%jxyz)/self%nbins(self%jxyz)
    

	i = 1
	Area = 0.d0
    do j = 1,self%nbins(self%ixyz)
    do k = 1,self%nbins(self%jxyz)
        patchArea = self%get_bilinear_patch_area(i, j, k)
		Area = Area + patchArea
        !print*,"intrinsic_area_bilinear", j,k,patchArea,Area,dx*dy
	enddo	
	enddo

end function intrinsic_area_bilinear




! Bilinear patch we cannot get exactly so we use 
! A = int int sqrt[ (df/dx)^2 + (df/dy)^2 + 1] dx dy
! and use trapzium rule to integrate
function get_bilinear_patch_area(self, i, j, k) result(Area)
    use librarymod, only : linspace, integrate_trap2d
    implicit none

	class(intrinsic_surface_bilinear), intent(in) :: self

    double precision	 :: Area
	integer, intent(in) :: i, j, k

    integer :: n, m, xres, yres
    double precision                 :: xmin, ymin, xmax, ymax
    double precision, dimension(2,2) :: A
    double precision, dimension(:), allocatable :: x, y
    double precision, dimension(:,:), allocatable :: p, f

	!1=Bottom left, 2=Bottom right, 3=Top left, 4=Top right
	p = self%indices_to_points(i, j, k)
    xmin = p(1,self%ixyz); xmax = p(4,self%ixyz)
    ymin = p(1,self%jxyz); ymax = p(4,self%jxyz)

    xres = 3; yres = 3
    x = linspace(xmin, xmax, xres)
    y = linspace(ymin, ymax, yres)
    A = self%Abilinear(:,:,j,k)
    allocate(f(xres, yres))
    do n=1,xres
    do m=1,yres
        f(n,m) = sqrt( (A(1,2)+A(2,2)*x(n))**2 &
                      +(A(2,1)+A(2,2)*y(m))**2 + 1)
    enddo
    enddo
    call integrate_trap2d(f,x(2)-x(1),y(2)-y(1),xres,yres,Area)


    ! Bilinear patch we cannot get exactly so we use 
    ! A = int int sqrt[ (df/dx)^2 + (df/dy)^2 + 1] dx dy
    ! and use trapzium rule to integrate
    !allocate(f(size(p)/2, size(p)/2))
    !do n=1,4
    !    call self%get_surface_derivative_bilinear(p, self%Abilinear(:,:,j,k), dSdr)
    !    f(n/2,mod(n,2)) = sqrt( dSdr(n,1)**2 +dSdr(n,2)**2 + 1)

        !f(n/2,mod(n,2)) = sqrt( (self%Abilinear(1,2,j,k)+self%Abilinear(2,2,j,k)*x(n))**2 &
        !                       +(self%Abilinear(2,1,j,k)+self%Abilinear(2,2,j,k)*y(n))**2 + 1)
    !enddo
    !dx = self%box(self%ixyz)/self%nbins(self%ixyz)
    !dy = self%box(self%jxyz)/self%nbins(self%jxyz)
    !call integrate_trap2d(f,dx,dy,2,2,Area)

    !To a first approximation, can assume xy coefficient is zero
    !const = sqrt(self%Abilinear(1,2,j,k)**2 +  self%Abilinear(2,1,j,k)**2 + 1)
    !Area = dx*dy*const
    

	!p = self%indices_to_points(i, j, k)
	!x = p(:,self%ixyz); y = p(:,self%jxyz)
    !Area = ( x(1)*y(1) - x(2)*y(2) &
	!	    -x(3)*y(3) + x(4)*y(4) )*const

end function get_bilinear_patch_area

!Get volume under surface 
function get_bilinear_patch_volume(self, i, j, k) result(Area)
    implicit none

	class(intrinsic_surface_bilinear), intent(in) :: self

    double precision	 :: Area

	integer :: i, j, k
    double precision, dimension(4) :: x, y
    double precision, dimension(:), allocatable :: f
    double precision, dimension(:,:), allocatable :: p
    double precision, dimension(:,:), allocatable :: hp

	p = self%indices_to_points(i, j, k)
	!Integral of bilinear patch is [x*y*xi(x/2, y/2)] for all 4 corners
	!where xi(x,y) is the surface value
	allocate(hp(size(p,1),size(p,2)))
	hp = 0.5d0*p
	call self%get_surface_bilinear(hp, self%Abilinear(:,:,j,k), f)
	!1=Bottom left, 2=Bottom right, 3=Top left, 4=Top right
	x = p(:,self%ixyz); y = p(:,self%jxyz)
    print*, "get_bilinear_patch_volume", j,k, x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4), f

	Area = 	x(1)*y(1)*f(1) - x(2)*y(2)*f(2) &
		  - x(3)*y(3)*f(3) + x(4)*y(4)*f(4)

end function get_bilinear_patch_volume


! subroutine modes_surface_to_bilinear_surface(ISR, nbins, Abilinear)
    ! implicit none

	! type(intrinsic_surface_real), intent(in) :: ISR	! declare an instance

    ! integer,intent(in),dimension(3) :: nbins
    ! real(kind(0.d0)),intent(out),dimension(:,:,:,:),allocatable ::Abilinear

    ! logical          :: debug=.true.
    ! logical, save    :: first_time = .true.
    ! integer          :: i, j, k, n, ixyz, jxyz, gbin(3), fileno
    ! real(kind(0.d0)) :: ysb, yst, zsb, zst, binsize(2)
    ! real(kind(0.d0)) :: y(2,2), z(2,2), P(2,2), A(2,2)
    ! real(kind(0.d0)), dimension(:), allocatable :: elevation
    ! real(kind(0.d0)), dimension(:,:), allocatable :: points

    ! !Start from the largest (end) value
    ! ixyz = ISR%ixyz
    ! jxyz = ISR%jxyz

    ! allocate(Abilinear(2,2,nbins(ixyz), nbins(jxyz)))
    ! allocate(points(4,2))

    ! binsize(1) = ISR%box(ixyz)/dble(nbins(ixyz))
    ! binsize(2) = ISR%box(jxyz)/dble(nbins(jxyz))

    ! n = 1
    ! do j = 1,nbins(ixyz)
    ! do k = 1,nbins(jxyz)
        ! !gbin = globalise_bin((/ 1, j, k /))
        ! gbin(1) = j; gbin(2) = k
        ! ysb = float(gbin(1)-1)*binsize(1)-0.5d0*ISR%box(ixyz)
        ! yst = float(gbin(1)  )*binsize(1)-0.5d0*ISR%box(ixyz)
        ! zsb = float(gbin(2)-1)*binsize(2)-0.5d0*ISR%box(jxyz)
        ! zst = float(gbin(2)  )*binsize(2)-0.5d0*ISR%box(jxyz)

        ! points(1,1) = ysb; points(1,2) = zsb
        ! points(2,1) = yst; points(2,2) = zsb 
        ! points(3,1) = ysb; points(3,2) = zst 
        ! points(4,1) = yst; points(4,2) = zst 

        ! call ISR%get_surface(points, elevation, include_zeromode=.true.)
        ! !call surface_from_modes(points, 3, q_vectors, modes, elevation)
        ! print*, "elevation modes", elevation(:)
        ! P = reshape(elevation, (/2,2/))
        ! call get_bilinear_surface_coeff((/ysb, yst/), (/zsb, zst/), P, A)
        ! Abilinear(:,:,j,k) = A(:,:)

        ! !DEBUG code here
        ! call get_surface_bilinear(points, A, elevation)
        ! !fileno = get_new_fileunit() 
        ! !open(fileno, file="surface.obj")
        ! do i =1, 4
            ! if ((elevation(i)- P(mod(i+1,2)+1,int(ceiling(i/2.d0)))) .gt. 1e-10) then
                ! print'(a,3i5,f16.4, 5f10.5, e18.8)', "elevation bilinear",j,k,i,A(:,:),points(i,:), &
                                        ! elevation(i)-P(mod(i+1,2)+1,int(ceiling(i/2.d0)))
            ! endif

            ! !Write to wavefunction obj format
            ! !write(fileno,'(a1, f15.8)') "v", points(i,1), points(i,2), elevation(i)
        ! enddo
        ! !write(fileno,'(a1, i8)'), "f", n, n+1, n+2, n+3
       ! ! n = n + 4
        ! !close(fileno)

    ! enddo
    ! enddo

! end subroutine modes_surface_to_bilinear_surface


subroutine fit_intrinsic_surface_bilinear(self, points, tau, ns, pivots)
    implicit none

	class(intrinsic_surface_bilinear) :: self	! declare an instance

    integer, dimension(:), allocatable, intent(inout) :: pivots
    double precision, intent(in) :: tau, ns
    double precision, intent(in), dimension(:,:), allocatable ::  points

    call fit_intrinsic_surface_modes(self, points, tau, ns, pivots)
    !call self%fit_intrinsic_surface(points, tau, ns, pivots)
    !call self%intrinsic_surface_real%fit_intrinsic_surface_modes(points, tau, ns, pivots)

	!print*, "DEBUG in fit_intrinsic_surface_bilinear, setting coeff to zero"
	!self%coeff = 0.d0
	!self%coeff(314) = 2.5d0
	
	!Update the saved surface at locations to use for quick surface 
	!interaction calculations 
	call self%update_sampled_surface(self%nbins, self%nhb)

	!stop "end of fit_intrinsic_surface_bilinear"

end subroutine fit_intrinsic_surface_bilinear



end module intrinsic_module





!program test
!    use intrinsic_module
!    use fruit
!    implicit none

!	type(intrinsic_surface_real) :: ISR, ISR_x	! declare an instance

!    logical,save                    :: first_time=.true., first_time_coeff=.true.
!    character(32)                   :: filename, debug_outfile
!    integer                         :: n,i,j,ixyz, jxyz,resolution,fittype,normal, clustNo, bins(3)
!    integer, dimension(3)           :: nbins, bin, binx, temp

!    integer, dimension(:), allocatable :: pivots
!    double precision                :: alpha, tau, eps, ns, area
!    double precision, parameter :: tolerance = 1e-15
!    double precision, dimension(3)  :: bintop, binbot, box, globaldomain, binsize !Used in CV
!    double precision,dimension(6)   :: extents
!    double precision,dimension(4)   :: ptin, pbin
!    double precision,dimension(:),allocatable :: x,y,z,f,pt,pb, elevation, zelevation
!    double precision,dimension(:),allocatable,save :: coeffmdt
!    double precision,dimension(:,:),allocatable :: rnp, extents_grid, dSdx, dSdz
!    double precision,dimension(:,:),allocatable :: zpoints, points
!    double precision,dimension(:,:,:),allocatable :: vertices


!    !Intrinsic surface coefficients
!    alpha =  0.5d0
!    tau =  1.d0
!    eps =  0.00000001d0
!    ns =  0.8d0
!    globaldomain = (/47.622031559045986, 12.699208415745595, 12.699208415745595/)
!    box = (/globaldomain(2), globaldomain(3), globaldomain(1)/)
!    allocate(rnp(3,1795))
!    do i =1, size(rnp,2)
!        read(586410,*) n,rnp(:,i)
!    enddo

!    !Z normal case
!    normal = 3
!    call ISR%initialise(box, normal, alpha, eps)   ! initialise
!    allocate(zpoints(size(rnp,2), size(rnp,1)))
!    zpoints = 0.d0
!    do i =1, size(rnp,2)
!        zpoints(i,1) = rnp(2,i)
!        zpoints(i,2) = rnp(3,i)
!        zpoints(i,3) = rnp(1,i)
!    enddo
!    call fit_intrinsic_surface_modes(zpoints, ISR, tau, ns, pivots)

!    !X normal case
!    call ISR_x%initialise(globaldomain, 1, alpha, eps)   ! initialise
!    allocate(points(size(rnp,2), size(rnp,1)))
!    points = 0.d0
!    do i =1, size(rnp,2)
!        points(i,1) = rnp(1,i)
!        points(i,2) = rnp(2,i)
!        points(i,3) = rnp(3,i)
!    enddo
!    call fit_intrinsic_surface_modes(points, ISR_x, tau, ns, pivots)

!    !Unit testing in Fortran
!    call init_fruit

!    do i =1, size(ISR_x%coeff,1)
!        call assert_equals(ISR%coeff(i), ISR_x%coeff(i), tolerance)
!!        if (ISR%coeff(i) .ne. ISR_x%coeff(i)) then
!!            print*, "Error in coeff", i, ISR%coeff(i), ISR_x%coeff(i)
!!        endif
!    enddo

!    call ISR_x%get_surface(points, elevation)
!    call ISR%get_surface(zpoints, zelevation)
!    do i =1, size(rnp,2)
!        call assert_equals(zelevation(i), elevation(i), tolerance)
!!        if (zelevation(i) .ne. elevation(i)) then
!!            print*, "Error in get surface", i, points(i,:), zpoints(i,:), elevation(i), zelevation(i)
!!        endif
!    enddo
!                    
!    call ISR_x%get_surface_derivative(points, dSdx)
!    call ISR%get_surface_derivative(zpoints, dSdz)
!    do i =1, size(rnp,2)
!        call assert_equals(dSdx(i,1), dSdz(i,1), tolerance)
!        call assert_equals(dSdx(i,2), dSdz(i,2), tolerance)
!!        if (any(dSdx(i,:) .ne. dSdz(i,:))) then
!!            print*, "Error in get surface_derivative", i, points(i,:), zpoints(i,:), dSdx(i,:), dSdz(i,:)
!!        endif
!    enddo

!    nbins = (/10, 5, 5/)
!    binsize = globaldomain/nbins
!    do i =1, size(rnp,2)
!        binx = ISR_x%get_bin(rnp(:,i), binsize)
!        temp = ISR%get_bin((/rnp(2,i), rnp(3,i), rnp(1,i)/), & 
!                          (/binsize(2), binsize(3), binsize(1)/))
!        bin(:) = (/temp(3),temp(1),temp(2)/)
!        call assert_equals(binx(1), bin(1))
!        call assert_equals(binx(2), bin(2))
!        call assert_equals(binx(3), bin(3))
!    enddo

!    call fruit_summary
!    call fruit_finalize

!                   

!end program test



!subroutine update_sampled_surface_opt(self, nbins, nhb)
!    implicit none

!	class(intrinsic_surface_bilinear) :: self

!    integer,intent(in),dimension(3), optional :: nbins, nhb

!    logical          :: debug=.false., writeobj=.false.
!    logical, save    :: first_time = .true.
!    integer          :: n, j, k, Allbins
!    real(kind(0.d0)) :: v(4), binsize(3), A(2,2)
!    real(kind(0.d0)), dimension(:), allocatable :: elevation
!    real(kind(0.d0)), dimension(:,:), allocatable :: points, lpoints

!    integer,dimension(3) 	:: nbins_, nhb_
! 	
!    print*, "update_sampled_surface_opt"

!	if (present(nbins) .and. present(nhb)) then
!		binsize = self%box/float(nbins)
!		nbins_ = nbins
!		nhb_ = nhb
!	else
!		binsize = self%binsize
!		nhb_ = self%nhb
!		nbins_ = self%nbins
!	endif

!    Allbins =  (self%nbins(self%ixyz)+2*self%nhb(self%ixyz)) &
!	          *(self%nbins(self%jxyz)+2*self%nhb(self%jxyz))
!	allocate(points(4*Allbins,3))
!	allocate(lpoints(4,3))
!    n = 0
!    do j = 1,self%nbins(self%ixyz)+2*self%nhb(self%ixyz)
!    do k = 1,self%nbins(self%jxyz)+2*self%nhb(self%jxyz)

!		v = self%index_to_vertex(j,k)
!		points(n+1,self%ixyz) = v(1); points(n+1,self%jxyz) = v(2) !Bottom left
!		points(n+2,self%ixyz) = v(3); points(n+2,self%jxyz) = v(2) !Bottom right
!		points(n+3,self%ixyz) = v(1); points(n+3,self%jxyz) = v(4) !Top left
!		points(n+4,self%ixyz) = v(3); points(n+4,self%jxyz) = v(4) !Top right
!        n = n + 4

!    enddo
!    enddo

!    !Get all points in one go, might be more efficient
!    select type (self)
!    type is (intrinsic_surface_chebychev)
!        call get_chebychev_surface(self, points, elevation, include_zeromode=.true.)
!    class is (intrinsic_surface_bilinear)
!        call get_real_surface(self, points, elevation, include_zeromode=.true.)
!    end select

!    n = 0
!    do j = 1,self%nbins(self%ixyz)+2*self%nhb(self%ixyz)
!    do k = 1,self%nbins(self%jxyz)+2*self%nhb(self%jxyz)

!		lpoints(:,self%normal) = elevation(n+1:n+4)
!        self%intrnsc_smple( j ,  k ) = elevation(n+1)
!        self%intrnsc_smple(j+1,  k ) = elevation(n+2)
!        self%intrnsc_smple(j  , k+1) = elevation(n+3)
!        self%intrnsc_smple(j+1, k+1) = elevation(n+4)

!		call self%paramterise_bilinear_surface(lpoints, A)
!		self%Abilinear(:,:,j,k) = A(:,:)
!        n = n + 4

!    enddo
!    enddo

!end subroutine update_sampled_surface_opt

