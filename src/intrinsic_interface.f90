

!A module to provide object oriented version of intrinsic surface functions

module intrinsic_module

    double precision, parameter :: pi=3.141592653589793d0

    type :: intrinsic_surface_mock

        integer      :: normal, ixyz, jxyz
        integer, dimension(2)  :: modes_shape
        integer, dimension(:), allocatable     :: pivots

        double precision  :: alpha, eps
        double precision, dimension(3) :: box

	    contains
		    procedure :: initialise  => initialise_mock
		    procedure :: update_surface => update_surface_mock
		    procedure :: get_surface => get_surface_mock

    end type intrinsic_surface_mock


    type :: intrinsic_surface_complex

        integer      :: normal, ixyz, jxyz
        integer, dimension(2)  :: modes_shape
        integer, dimension(:), allocatable     :: pivots

        double precision  :: alpha, eps
        double precision, dimension(3) :: box
        double precision, dimension(:), allocatable     :: Q
        double precision, dimension(:,:), allocatable   :: Qxy
        double complex, dimension(:,:), allocatable     :: modes
        double precision, dimension(:,:,:), allocatable :: q_vectors

	    contains
		    procedure :: initialise  => compute_q_vectors
		    procedure :: update_surface => update_surface_modes
		    procedure :: get_surface => get_surface
!		    procedure :: sample_surface => sample_intrinsic_surface

    end type intrinsic_surface_complex

    type :: intrinsic_surface_real

        integer      :: normal, ixyz, jxyz, n_waves, n_waves2, qm
        integer, dimension(:), allocatable     :: u, v
        integer, dimension(:), allocatable     :: pivots

        double precision  :: alpha, eps, area
        double precision, dimension(3) :: box
        double precision, dimension(:), allocatable     :: diag, coeff, coeff_mdt
        double precision, dimension(:,:), allocatable   :: diag_matrix

	    contains
		    procedure :: initialise  => initialise
		    procedure :: update_surface => update_real_surface
		    procedure :: get_surface => get_real_surface
		    procedure :: get_surface_derivative => get_real_surface_derivative
		    procedure :: get_zero_mode => get_zero_mode
            procedure :: get_bin => get_bin_from_surface
		    procedure :: sample_surface => sample_intrinsic_surface

    end type intrinsic_surface_real


contains 


subroutine initialise_mock(self, box, normal, alpha, eps)
    implicit none

	class(intrinsic_surface_mock) :: self

    integer, intent(in)         :: normal
    double precision, intent(in) :: alpha, eps
    double precision, intent(in), dimension(3) :: box

    integer :: i, j

    self%normal = normal
    self%box = box
    self%alpha = alpha
    self%eps = eps
    self%ixyz = mod(normal,3)+1
    self%jxyz = mod(normal+1,3)+1

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

!Function 

function wave_function(x, u, Lx)
    implicit none

    integer,intent(in) :: u
    real(kind(0.d0)),intent(in) :: Lx
    real(kind(0.d0)),intent(in),dimension(:) :: x

    real(kind(0.d0)),dimension(size(x))    :: wave_function
    integer :: i

    do i =1,size(x,1)
        if (u .ge. 0) then
            wave_function(i) = cos(2.d0 * pi * u * x(i) / Lx)
        else
            wave_function(i) = sin(2.d0 * pi * abs(u) * x(i) / Lx)
        endif
    enddo

end function wave_function


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

subroutine initialise(self, box, normal, alpha, eps)
    implicit none

	class(intrinsic_surface_real) :: self

    integer, intent(in)         :: normal
    double precision, intent(in) :: alpha, eps
    double precision, intent(in), dimension(3) :: box

    integer :: i, j

    self%normal = normal
    self%box = box
    self%alpha = alpha
    self%eps = eps
    self%ixyz = mod(normal,3)+1
    self%jxyz = mod(normal+1,3)+1

    self%area = box(self%ixyz)*box(self%jxyz)
    self%qm = int(sqrt(self%area)/alpha)
    self%n_waves = 2*self%qm+1
    self%n_waves2 = self%n_waves**2
    allocate(self%u(self%n_waves2), self%v(self%n_waves2))
    do i=1,self%n_waves2
        self%u(i) = (i-1) / self%n_waves - self%qm
        self%v(i) = modulo((i-1), self%n_waves) - self%qm
    enddo
	
	allocate(self%coeff(self%n_waves2))
	self%coeff = 0.d0

    allocate(self%diag(self%n_waves2))
    self%diag = check_uv(self%u, self%v) & 
                 * (  self%u**2 * box(self%jxyz) / box(self%ixyz) & 
                    + self%v**2 * box(self%ixyz) / box(self%jxyz))
    allocate(self%diag_matrix(size(self%diag,1),size(self%diag,1)))
    self%diag_matrix = 0.d0
    do i=1,size(self%diag,1)
    do j=1,size(self%diag,1)
        if (i .eq. j) then
            self%diag_matrix(i,j) = 4.d0 * pi**2 * eps * self%diag(i)
        endif
    enddo
    enddo

end subroutine initialise

subroutine update_real_surface(self, points)
#if USE_LAPACK
    use lapack_fns, only : solve, multiply_by_tranpose
#endif
    implicit none

	class(intrinsic_surface_real) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points

    integer :: j
    double precision, dimension(:), allocatable :: b
    double precision, dimension(:,:), allocatable :: A, fuv

    logical :: debug=.false.
    double precision, dimension(:), allocatable :: surf
    double precision, allocatable, dimension(:) :: x
	
    allocate(A(self%n_waves2, self%n_waves2))
    allocate(b(self%n_waves2))
    allocate(fuv(size(points,1), self%n_waves2))

    !Save Previous solution
    !if (.not. allocated(self%coeff_mdt)) allocate(self%coeff_mdt)
    !self%coeff_mdt = self%coeff

    A = 0.d0
    b = 0.d0
    do j =1, self%n_waves2
        fuv(:,j) = wave_function(points(:,self%ixyz), self%u(j), self%box(self%ixyz)) & 
                  *wave_function(points(:,self%jxyz), self%v(j), self%box(self%jxyz))
        b(j) = b(j) + sum(points(:,self%normal) * fuv(:,j))
    enddo

    A = 0.d0
    call multiply_by_tranpose(fuv, fuv, A)

    !This solves Ax = b
    A = A + self%diag_matrix
    call solve(A, b, x)
	self%coeff(:) = x(:)

    !Check surface we just fitted actually matches our points
    if (debug) then
        call self%get_surface(points, surf)
        do j = 1, size(points,1)
            print*, "update_real_surface DEBUG ON Elevation vs points = ", & 
					j, surf(j) , points(j,self%normal), &
                     surf(j)-points(j,self%normal)
        enddo
!    do j=1, self%n_waves2
!        print*, "Coeffs = ", j, self%coeff(j)
!    enddo

    endif

end subroutine update_real_surface


subroutine get_real_surface(self, points, elevation, qu)
    implicit none

	class(intrinsic_surface_real) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(in), optional ::  qu
    double precision, intent(out), dimension(:), allocatable :: elevation

    integer :: j, ui, vi, qu_

    if (present(qu)) then
        qu_ = qu
    else
        qu_ = self%qm
    endif

    allocate(elevation(size(points,1)))
    elevation = 0.d0
    do ui = -qu_, qu_
    do vi = -qu_, qu_
        j = (2 * self%qm + 1) * (ui + self%qm) + (vi + self%qm) + 1
        elevation(:) = elevation(:) + self%coeff(j) & 
                     * wave_function(points(:,self%ixyz), ui, self%box(self%ixyz)) & 
                     * wave_function(points(:,self%jxyz), vi, self%box(self%jxyz)) 
    enddo
    enddo

end subroutine get_real_surface

subroutine get_real_surface_derivative(self, points, dSdr, qu)
    implicit none

	class(intrinsic_surface_real) :: self

    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(in), optional ::  qu
    double precision, intent(out), dimension(:,:), allocatable :: dSdr

    integer :: j, ui, vi, qu_

    if (present(qu)) then
        qu_ = qu
    else
        qu_ = self%qm
    endif

    allocate(dSdr(size(points,1),2))
    dSdr = 0.d0
    do ui = -qu_, qu_
    do vi = -qu_, qu_
        j = (2 * self%qm + 1) * (ui + self%qm) + (vi + self%qm) + 1
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


subroutine get_zero_mode(self, zeromode)
    implicit none

	class(intrinsic_surface_real) :: self

    integer :: j
    real(kind(0.d0)) :: zeromode

    j = (2 * self%qm + 1) * self%qm + self%qm + 1
    zeromode = self%coeff(j)

end subroutine get_zero_mode


!A version getting the explicit location from the intrinsic surface
function get_bin_from_surface(self, r, nbins, nhb) result(bin)
    implicit none

	class(intrinsic_surface_real) :: self

    real(kind(0.d0)),intent(in),dimension(3) :: r
    integer,dimension(3),intent(in)		     :: nbins, nhb

    integer,dimension(3)	                 :: bin

    integer :: n, i, j
    real(kind(0.d0)) :: maprange, zeromode
    real(kind(0.d0)), dimension(3) :: halfdomain, binsize
    real(kind(0.d0)), allocatable, dimension(:) :: elevation
    real(kind(0.d0)), allocatable, dimension(:,:) :: points

    binsize = self%box/float(nbins)
    halfdomain = 0.5*self%box
    n=self%normal
    i=self%ixyz
    j=self%jxyz
	
    !Add in a range over which the intrinsic deformation is applied
	!maprange = 5.d0
    call self%get_zero_mode(zeromode)
	!if ((r(n) .lt. zeromode-maprange) .or. & 
	!	(r(n) .gt. zeromode+maprange)) then
	!	bin(n) = ceiling((r(n)+halfdomain(n)-zeromode+0.5d0*binsize(n))/binsize(n))+nhb(n)
	!    bin(i) = ceiling((r(i)+halfdomain(i))/binsize(i))+nhb(i)
	!	bin(j) = ceiling((r(j)+halfdomain(j))/binsize(j))+nhb(j)
	!	return
	!endif
	
    allocate(points(1,3))
    points(1,:) = r(:)
    call self%get_surface(points, elevation)

    !Added a shift by zero wavelength so surface is not at zero
    bin(n) = ceiling((r(n)+halfdomain(n)-elevation(1)+zeromode+0.5d0*binsize(n))/binsize(n))+nhb(n) !HALF SHIFT
    bin(i) = ceiling((r(i)+halfdomain(i))/binsize(i))+nhb(i)
    bin(j) = ceiling((r(j)+halfdomain(j))/binsize(j))+nhb(j)

	if (bin(n) > nbins(n)+nhb(n)) then
        bin(n) = nbins(n)+nhb(n)
    elseif (bin(n) < 1 ) then
        bin(n) = 1
    endif

end function get_bin_from_surface

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
!=========================================
!=========================================

subroutine compute_q_vectors(self, box, normal, alpha, eps)
    use librarymod, only : meshgrid2d
    implicit none

	class(intrinsic_surface_complex) :: self

    integer, intent(in)         :: normal
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

    ! Least square solution solving ph*z = s
    call pinverse(ph, pinv_ph)

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
        call self%get_surface(points, surf)
        do i = 1, size(points,1)
            print*, "Elevation vs points = ", i, surf(i) , points(i,self%normal), &
                     surf(i)-points(i,self%normal),   surf(i)/points(i,self%normal)
        enddo
        print'(a,10f10.5)', "Modes", self%modes(1,1), self%modes(1,2), & 
                    self%modes(2,1), self%modes(1,3), self%modes(3,1)
        !stop "Stop in debug of update_surface_modes" 
    endif

end subroutine update_surface_modes



subroutine get_surface(self, points, elevation)
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

end subroutine get_surface



subroutine sample_intrinsic_surface(self, nbins, vertices, writeiter)
    use librarymod, only : get_new_fileunit, get_Timestep_FileName
    implicit none

	class(intrinsic_surface_real) :: self

    integer,intent(in),dimension(3) :: nbins

    real(kind(0.d0)), dimension(:,:,:),intent(out), allocatable :: vertices
    integer, intent(in), optional :: writeiter


    logical          :: debug=.true., writeobj=.false.
    logical, save    :: first_time = .true.
    integer          :: i, j, k, n, v, ixyz, jxyz, fileno, qm, qu
    real(kind(0.d0)) :: ysb, yst, zsb, zst, binsize(2), area
    real(kind(0.d0)), dimension(:), allocatable :: elevation
    real(kind(0.d0)), dimension(:,:), allocatable :: points
    character(200) :: outfile_t, filename

    if (present(writeiter)) then      
        writeobj = .true.
    else
        writeobj = .false.
    endif

    binsize(1) = self%box(self%ixyz)/dble(nbins(self%ixyz))
    binsize(2) = self%box(self%jxyz)/dble(nbins(self%jxyz))

    allocate(points(4,3))
    allocate(vertices(nbins(self%ixyz)*nbins(self%jxyz),4,3))

    if (writeobj) then
        fileno = get_new_fileunit() 
        call get_Timestep_FileName(writeiter,"./results/surface",outfile_t)
        write(filename,'(a,a4)') trim(outfile_t),'.obj'
        open(fileno, file=trim(filename))
    endif

    n = 1; v = 1
    do j = 1,nbins(self%ixyz)
    do k = 1,nbins(self%jxyz)
        ysb = float(j-1)*binsize(1)-0.5d0*self%box(self%ixyz)
        yst = float(j  )*binsize(1)-0.5d0*self%box(self%ixyz)
        zsb = float(k-1)*binsize(2)-0.5d0*self%box(self%jxyz)
        zst = float(k  )*binsize(2)-0.5d0*self%box(self%jxyz)

        points(1,self%ixyz) = ysb; points(1,self%jxyz) = zsb !Bottom left
        points(2,self%ixyz) = yst; points(2,self%jxyz) = zsb !Bottom right
        points(3,self%ixyz) = ysb; points(3,self%jxyz) = zst !Top left
        points(4,self%ixyz) = yst; points(4,self%jxyz) = zst !Top right

        call self%get_surface(points, elevation)
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


subroutine get_initial_pivots(points, ISR, pivots)
    use librarymod, only : Qsort
    implicit none

	type(intrinsic_surface_real), intent(in) :: ISR	! declare an instance
    double precision, intent(in), dimension(:,:), allocatable ::  points
    integer, intent(out), dimension(:), allocatable ::  pivots

    ! Defines the initial pivots as a set of 9 particles, where
    ! each particle is in a distinct sector formed by dividing
    ! the macroscopic plane into 3x3 regions.
    integer :: Npivots, maxpivots, i, ind
    integer, dimension(2) :: nxy, bins
    integer, dimension(:,:), allocatable :: sectors

    integer, dimension(:), allocatable :: indices
    double precision, dimension(3) :: binsize
    double precision, dimension(:), allocatable :: z

    !Define bins
    bins = 3
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

    !print*, ixyz, jxyz, box(:), box(ixyz)/dble(bins(1))
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

    if (Npivots .ne. 0) stop "Not all initial pivots found"

end subroutine get_initial_pivots



subroutine update_pivots(points, ISR, pivots, tau, new_pivots)
    use librarymod, only : Qsort
    implicit none

	type(intrinsic_surface_real), intent(in) :: ISR	! declare an instance
    integer, intent(in), dimension(:), allocatable ::  pivots
    double precision, intent(in) :: tau
    double precision, intent(in), dimension(:,:), allocatable ::  points
    integer, intent(out), dimension(:), allocatable ::  new_pivots

    logical :: found_range
    integer :: i, n, ixyz, jxyz, ind, nPivots, sp, qm, qu
    double precision :: z_max, z_min, area
    integer, dimension(:), allocatable :: candidates, updated_pivots, indices
    double precision, dimension(:), allocatable :: z, surf
    double precision, dimension(:,:), allocatable :: candidates_pos, pivot_pos

    ! Searches for points within a distance tau from the
    ! interface.
    sp = size(pivots,1)
    allocate(pivot_pos(sp,3)) 
    do i =1, sp
        pivot_pos(i,:) = points(pivots(i),:)
    enddo
    z_max = maxval(pivot_pos(:,ISR%normal)) + ISR%alpha * 2.d0
    z_min = minval(pivot_pos(:,ISR%normal)) - ISR%alpha * 2.d0

    z = points(:, ISR%normal)
    allocate(indices(size(points,1))) 
    do ind = 1, size(indices,1)
        indices(ind) = ind
    enddo
    call Qsort(z, indices)
    n = 0; found_range=.false.
    allocate(candidates(size(points,1))) 
    do i = size(z,1), 1, -1
        if ((z(i) > z_min) .and. (z(i) < z_max)) then
            n = n + 1
            candidates(n) = indices(i)
            found_range = .true.
            !print*, "values", i, n, z_min, z(i), z_max, candidates(n)
        else if (found_range) then
            !If z is sorted and we've been inside the range,
            !once we leave, no point looping anymore
            exit
        endif

    enddo

    !Get positions from indices
    allocate(candidates_pos(n,3))
    do i =1, n
        !print*, "values", i, n, candidates(i), points(candidates(i),:)
        candidates_pos(i,:) = points(candidates(i),:)
    enddo

    !Recalculate surface at candidate pivot locations
    call ISR%get_surface(candidates_pos, surf)

    nPivots = 0
    allocate(updated_pivots(n))
    do i =1,n
        if ((surf(i)-candidates_pos(i,ISR%normal))**2 .lt. tau**2) then
            nPivots = nPivots + 1
            updated_pivots(nPivots) = candidates(i)
        endif
    enddo
    allocate(new_pivots(nPivots))
    new_pivots = updated_pivots(1:nPivots)

    !allocate(new_pivots(nPivots+sp))
    !new_pivots(1:sp) = pivots(:)
    !new_pivots(sp+1:) = updated_pivots(:nPivots)

    !print*, new_pivots 

end subroutine update_pivots


!An alternative way of updating pivots using a surface which 
!gets atoms a distance tau from existing surface
subroutine update_pivots_alt(points, ISR, pivots, tau, ns, new_pivots)
    use librarymod, only : Qsort
    implicit none

	type(intrinsic_surface_real), intent(in) :: ISR	! declare an instance
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
    call ISR%get_surface(points, surf)
    allocate(z(size(points,1)))
    z = abs(points(:,ISR%normal)-surf(:))
    allocate(indices(size(points,1)))
    do ind = 1, size(indices,1)
        indices(ind) = ind
    enddo
    call Qsort(z, indices)
    n = 0; found_range=.false.
    allocate(candidates(size(points,1))) 
    do i = 1, size(z,1)
        if (z(i) .lt. tau) then
            n = n + 1
            candidates(n) = indices(i)
            found_range = .true.
            if (n .ge. nmols) exit
            !print*, "update_pivots_alt in range", i, n, z(i), tau, candidates(n)
        else if (found_range) then
            !If z is sorted and we've been inside the range,
            !once we leave, no point looping anymore
            !print*, "update_pivots_alt out range", i, n, z(i), tau, candidates(n)
            exit
        endif
    enddo
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



subroutine fit_intrinsic_surface_modes(points, ISR, tau, ns, pivots)
    !DEBUG DEBUGDEBUGDEBUG
    !use physical_constants_MD, only : np
    !use calculated_properties_MD, only : nbins
    !use computational_constants_MD, only : nhb
    !DEBUG DEBUGDEBUGDEBUG

    implicit none

	type(intrinsic_surface_real), intent(in) :: ISR	! declare an instance

    double precision, intent(in) :: tau, ns
    double precision, intent(in), dimension(:,:), allocatable ::  points
    integer, dimension(:), allocatable, intent(inout) :: pivots

    integer :: i, j, ixyz, jxyz, sp, ntarget, try, maxtry=100, qm, qu
    integer, dimension(2) :: modes_shape
    integer, dimension(:), allocatable :: indices, new_pivots, initial_pivots
    !integer, dimension(:), allocatable, save :: pivots
    double precision, save :: savetau=0.d0
    double precision :: Error, tau_, rand, diff, maxpivot
    double precision, dimension(:), allocatable :: Q, z, d, surf
    double precision, dimension(:,:), allocatable :: Qxy, pivot_pos, initial_pivot_pos

    !Define things
    tau_ = max(tau, savetau)
    ntarget = int(ns*ISR%area)

    !Get initial pivots
    call get_initial_pivots(points, ISR, initial_pivots)
    if (allocated(pivots)) deallocate(pivots)
    pivots = initial_pivots

    !Get positions and fits surface to initial pivots
    sp = size(pivots,1)
    allocate(pivot_pos(sp,3))
    do i =1, sp
        pivot_pos(i,:) = points(pivots(i),:)
        !print*, "initial pivot pos = ", i, pivot_pos(i,:)
    enddo

    call ISR%update_surface(pivot_pos)

!    do i=1,sp
!        print*, "Pivot error = ", pivot_pos(i,:), surf(i), pivot_pos(i,3)-surf(i)
!    enddo

    do try = 1, maxtry

        !Get new pivots
        !call update_pivots(points, ISR, pivots, tau_, new_pivots)
        call update_pivots_alt(points, ISR, pivots, tau_, ns, new_pivots)

        !Get new positions and new modes
        if (allocated(pivot_pos)) deallocate(pivot_pos)
        sp = size(new_pivots,1)
        allocate(pivot_pos(sp,3))
        do i =1, sp
            pivot_pos(i,:) = points(new_pivots(i),:)
            !print*, "new pivot pos = ", i, pivot_pos(i,:)
        enddo
        !call get_surface_modes(pivot_pos, Qxy, modes_shape, Q, modes, eps)
        call ISR%update_surface(pivot_pos)

        !Exit once we have converged to particles on surface / area = ns
        !print*, size(new_pivots)/area, size(new_pivots), area, ns 
        !if (size(new_pivots)/ISR%area .gt. ns)  then
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
!            maxpivot = maxval(initial_pivot_pos(:,ISR%normal))
!            do i =1, ntarget
!                diff = maxpivot - pivot_pos(i, ISR%normal)
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
            call ISR%get_surface(pivot_pos, surf)

            d = pivot_pos(:, ISR%normal) - surf(:)
            Error = sqrt(sum(d * d) / size(d))

            print*, "Try no. = ", try, "No. pivots = ", size(pivots), & 
                    "Error=", Error, "Area=", ISR%area, size(pivots)/ISR%area

            if (Error .gt. 1e2) then
                print*, "Solution appears to be diverging, reverting to initial pivots"
                deallocate(pivots)
                call get_initial_pivots(points, ISR, pivots)
                deallocate(pivot_pos)
                sp = size(pivots,1)
                allocate(pivot_pos(sp,3))
                do i =1, sp
                    pivot_pos(i,:) = points(pivots(i),:)
                    !print*, "initial pivot pos = ", i, pivot_pos(i,:)
                enddo
                call ISR%update_surface(pivot_pos)

            endif
        endif

    enddo

    !Plot updated surface error
    !DEBUG DEBUG DEBUG DEBUG
!    if (allocated(d)) deallocate(d)
!    allocate(d(sp))
!    call ISR%get_surface(pivot_pos, surf)

!    d = pivot_pos(:, ISR%normal) - surf(:)
!    Error = sqrt(sum(d * d) / size(d))
!    print*, "fit_intrinsic_surface_modes", try, Error, ntarget
    !DEBUG DEBUG DEBUG DEBUG

end subroutine fit_intrinsic_surface_modes

!Save coefficiients to matrix for bilinear 
subroutine get_bilinear_surface_coeff(x, y, P, A)
    implicit none

    real(kind(0.d0)),intent(in)  :: x(2), y(2)
    real(kind(0.d0)),intent(in)  :: P(2,2)
    real(kind(0.d0)),intent(out) :: A(2,2)
    real(kind(0.d0)) :: x1, x2, y1, y2, x12, y12, d
    x1 = x(1); x2 = x(2)
    y1 = y(1); y2 = y(2)
    x12 = x2 - x1; y12 = y2 - y1
    d = x12*y12
    A(1,1) = ( P(1,1)*x2*y2 - P(1,2)*x2*y1 - P(2,1)*x1*y2 + P(2,2)*x1*y1)/d
    A(2,1) = (-P(1,1)*y2    + P(1,2)*y1    + P(2,1)*y2    - P(2,2)*y1   )/d
    A(1,2) = (-P(1,1)*x2    + P(1,2)*x2    + P(2,1)*x1    - P(2,2)*x1   )/d
    A(2,2) = ( P(1,1)       - P(1,2)       - P(2,1)       + P(2,2)      )/d

end subroutine get_bilinear_surface_coeff


!Get surface position from positions for a known bin
subroutine get_surface_bilinear(points, A, elevation)
    implicit none

    real(kind(0.d0)), intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(in), dimension(2,2) :: A
    double precision, intent(out), dimension(:), allocatable :: elevation

    integer :: n, ixyz, jxyz

    allocate(elevation(size(points,1)))
    elevation = 0.d0
    do n=1,size(points,1)
    do ixyz=1,2
    do jxyz=1,2
        elevation(n) = elevation(n) + A(ixyz, jxyz) & 
                       *(points(n,1)**(ixyz-1))*(points(n,2)**(jxyz-1))
    enddo
    enddo
    enddo

end subroutine get_surface_bilinear


subroutine modes_surface_to_bilinear_surface(ISR, nbins, Abilinear)
    implicit none

	type(intrinsic_surface_real), intent(in) :: ISR	! declare an instance

    integer,intent(in),dimension(3) :: nbins
    real(kind(0.d0)),intent(out),dimension(:,:,:,:),allocatable ::Abilinear

    logical          :: debug=.true.
    logical, save    :: first_time = .true.
    integer          :: i, j, k, n, ixyz, jxyz, gbin(3), fileno
    real(kind(0.d0)) :: ysb, yst, zsb, zst, binsize(2)
    real(kind(0.d0)) :: y(2,2), z(2,2), P(2,2), A(2,2)
    real(kind(0.d0)), dimension(:), allocatable :: elevation
    real(kind(0.d0)), dimension(:,:), allocatable :: points

    !Start from the largest (end) value
    ixyz = mod(ISR%normal,3)+1
    jxyz = mod(ISR%normal+1,3)+1

    allocate(Abilinear(2,2,nbins(ixyz), nbins(jxyz)))
    allocate(points(4,2))

    binsize(1) = ISR%box(ixyz)/dble(nbins(ixyz))
    binsize(2) = ISR%box(jxyz)/dble(nbins(jxyz))

    n = 1
    do j = 1,nbins(ixyz)
    do k = 1,nbins(jxyz)
        !gbin = globalise_bin((/ 1, j, k /))
        gbin(1) = j; gbin(2) = k
        ysb = float(gbin(1)-1)*binsize(1)-0.5d0*ISR%box(ixyz)
        yst = float(gbin(1)  )*binsize(1)-0.5d0*ISR%box(ixyz)
        zsb = float(gbin(2)-1)*binsize(2)-0.5d0*ISR%box(jxyz)
        zst = float(gbin(2)  )*binsize(2)-0.5d0*ISR%box(jxyz)

        points(1,1) = ysb; points(1,2) = zsb
        points(2,1) = yst; points(2,2) = zsb 
        points(3,1) = ysb; points(3,2) = zst 
        points(4,1) = yst; points(4,2) = zst 

        call ISR%get_surface(points, elevation)
        !call surface_from_modes(points, 3, q_vectors, modes, elevation)
        !print*, "elevation modes", elevation(:)
        P = reshape(elevation, (/2,2/))
        call get_bilinear_surface_coeff((/ysb, yst/), (/zsb, zst/), P, A)
        Abilinear(:,:,j,k) = A(:,:)

        !DEBUG code here
        call get_surface_bilinear(points, A, elevation)
        !fileno = get_new_fileunit() 
        !open(fileno, file="surface.obj")
        do i =1, 4
            if ((elevation(i)- P(mod(i+1,2)+1,int(ceiling(i/2.d0)))) .gt. 1e-10) then
                print'(a,3i5,f16.4, 5f10.5, e18.8)', "elevation bilinear",j,k,i,A(:,:),points(i,:), &
                                        elevation(i)-P(mod(i+1,2)+1,int(ceiling(i/2.d0)))
            endif

            !Write to wavefunction obj format
            !write(fileno,'(a1, f15.8)') "v", points(i,1), points(i,2), elevation(i)
        enddo
        !write(fileno,'(a1, i8)'), "f", n, n+1, n+2, n+3
       ! n = n + 4
        !close(fileno)

    enddo
    enddo

end subroutine modes_surface_to_bilinear_surface


subroutine fit_intrinsic_surface(points, ISR, nbins, tau, ns, Abilinear)
    implicit none

	type(intrinsic_surface_real), intent(in) :: ISR	! declare an instance

    integer, dimension(3), intent(in) :: nbins
    double precision, intent(in) :: tau, ns
    double precision, intent(in), dimension(:,:), allocatable ::  points
    double precision, intent(out), dimension(:,:,:,:), allocatable :: Abilinear

    !Debug
    integer :: ixyz, jxyz, i,j,k
    double precision :: dx, dy
    double precision, allocatable, dimension(:,:) :: checkpoint
    integer, dimension(:), allocatable ::  pivots

    call fit_intrinsic_surface_modes(points, ISR, tau, ns, pivots)
    call modes_surface_to_bilinear_surface(ISR, nbins, Abilinear)

end subroutine fit_intrinsic_surface



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

