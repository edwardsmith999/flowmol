!=============================================================================
!
!				  Coupler 
!
! Routines accessible from application ( molecular or continuum )
! after the name, in parenthesis, is the realm in which each routine must be called
!
! coupler_create_comm	      	(cfd+md) splits MPI_COMM_WORLD, create intercommunicator between CFD and MD
!
! coupler_create_map	       	(cfd+md) creates correspondence maps between the CFD grid and MD domains
!
! coupler_cfd_init              (cfd)    initialises coupler with CFD data
!
! coupler_md_init               (cfd)    initialises coupler and set MD parameters with using dat from CFD or COUPLER.in
!
! coupler_cfd_adjust_domain     (cfd)    adjust CFD tomain to an integer number FCC or similar MD initial molecule layout
!
! coupler_send_grid_data        (cfd+md) sends grid data exchanged between realms ( generic interface)
!
! coupler_recv_grid_data        (cfd+md) receives data exchanged between realms ( generic interface)
!
! coupler_cfd_get               (cfd)    returns coupler internal parameters for CFD realm
!
! coupler_md_get                (md)     returns coupler internal parameters for MD realm
!
! coupler_md_get_save_period    (md)     auxiliary used for testing
!
! coupler_md_get_average_period (md)     gets average period of boundary condition
!
! coupler_md_get_md_steps_per_dt_cfd (md) returns the number of step MD does for ech CFD step
!
! coupler_md_get_nsteps         (md)     returm CFD nsteps
!      
! coupler_md_get_dt_cfd         (md)     returns MD dt
!
! coupler_md_set                (md)     sets zL if CFD is 2D
!
! coupler_md_get_density        (md)     gets CFD density
!
! coupler_md_get_cfd_id         (md)     id for CFD code, possible values set in coupler_parameters
!
!  Lucian Anton, November 2011
!  Revised April 2012
!
!=============================================================================
module coupler
	use coupler_parameters
    implicit none
    save

    interface coupler_send_grid_data
        module procedure coupler_send_grid_data_3d, coupler_send_grid_data_4d
    end interface

    interface coupler_recv_grid_data
        module procedure coupler_recv_grid_data_3d, coupler_recv_grid_data_4d
    end interface

    private coupler_send_grid_data_3d, coupler_send_grid_data_4d, &
        coupler_send_grid_data_xd, coupler_recv_grid_data_3d, coupler_recv_grid_data_4d,&
        coupler_recv_grid_data_xd

contains

!=============================================================================
! 							coupler_create_comm	      	
! (cfd+md) Splits MPI_COMM_WORLD in both the CFD and MD code respectively
! 		   and create intercommunicator between CFD and MD
!-----------------------------------------------------------------------------

subroutine coupler_create_comm(realm, realm_comm, ierror)
	use mpi
	use coupler_internal_common
	use coupler_input_data
	implicit none

	integer, intent(in) :: realm ! CFD or MD ?
	integer, intent(out):: realm_comm, ierror

	! test if we have a CFD and a MD realm
	ierror=0

	call test_realms
	COUPLER_REALM = realm
	call create_comm
	call  read_coupler_input

    ! stop if requested ( useful for development )
	call request_stop("create_comm") ! stops here if in COUPLER.in  stop requestis set to "create_comm"
contains

!-----------------------------------------------------------------------------
!	Test if CFD and MD realms are assigned correctly
!-----------------------------------------------------------------------------

subroutine test_realms
	use mpi
	implicit none

	integer i, myid, nproc, nf, np, ierr
	integer, allocatable :: ra(:)

	call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)

	if (myid .eq. 0) then
		call mpi_comm_size(mpi_comm_world, nproc, ierr)
		allocate(ra(nproc))
	else
					allocate(ra(1))	!Assign value arbitarily
	endif

	call mpi_gather(realm,1,MPI_INTEGER,ra,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	if (myid .eq. 0) then

		nf = 0
		np = 0

		do i =1, nproc
			if ( ra(i) .eq. COUPLER_CFD ) then 
				nf = nf + 1
			else if ( ra(i) .eq. COUPLER_MD ) then
				np = np +1
			else
				ierror = COUPLER_ERROR_REALM
				write(*,*) "wrong realm value in coupler_create_comm"
				call MPI_abort(MPI_COMM_WORLD,ierror,ierr)
			endif
		enddo

		if ( nf .eq. 0 .or. np .eq. 0) then 
			ierror = COUPLER_ERROR_ONE_REALM
			write(*,*) "CFD or MD  realm is missing in MPI_COMM_WORLD"
			call MPI_abort(MPI_COMM_WORLD,ierror,ierr)
		endif



	endif

end subroutine test_realms

!-----------------------------------------------------------------------------
! Create communicators for each realm and inter-communicator
!-----------------------------------------------------------------------------

subroutine create_comm
	use mpi
	use coupler_internal_common, only : COUPLER_REALM_COMM, COUPLER_GLOBAL_COMM, COUPLER_ICOMM
	implicit none

	integer ierr, myid, myid_comm, myid_comm_max, realm, &
		iaux(2), jaux(2), remote_leader, comm, comm_size

	realm = COUPLER_REALM

	! get a global internal communicator for coupler operations
	call mpi_comm_dup(MPI_COMM_WORLD,COUPLER_GLOBAL_COMM,ierr)
	call mpi_comm_rank(COUPLER_GLOBAL_COMM,myid,ierr)
	REALM_COMM         = MPI_COMM_NULL
	COUPLER_REALM_COMM = MPI_COMM_NULL

	! get a communicator for each realm
	call mpi_comm_split(MPI_COMM_WORLD,realm,myid,REALM_COMM,ierr)

	! get internal, realm specific communicator for coupler operations
	call mpi_comm_split(COUPLER_GLOBAL_COMM,realm,myid,COUPLER_REALM_COMM,ierr)

	comm = COUPLER_REALM_COMM ! shorthand

	! create inter-communicators

	! Get the mpi_comm_world ranks that hold the largest ranks in cfd_comm and md_comm
	call mpi_comm_rank(comm,myid_comm,ierr)
	call mpi_comm_size(comm,comm_size,ierr)

	iaux(:) = -1
	jaux(:) = -1
	if ( myid_comm .eq. comm_size - 1) then
		iaux(COUPLER_REALM) = myid
	endif
	call mpi_allreduce( iaux ,jaux, 2, MPI_INTEGER, MPI_MAX, &
		COUPLER_GLOBAL_COMM, ierr)  

	select case (COUPLER_REALM)
	case (COUPLER_CFD)
		remote_leader = jaux(COUPLER_MD)
	case (COUPLER_MD)
		remote_leader = jaux(COUPLER_CFD)
	end select

	!print*,color, jaux, remote_leader

	call mpi_intercomm_create(comm, comm_size - 1, COUPLER_GLOBAL_COMM,&
									remote_leader, 1, COUPLER_ICOMM, ierr)

	write(0,*) 'did (inter)communicators ', code_name(COUPLER_REALM), myid

end subroutine create_comm

end subroutine coupler_create_comm


!=============================================================================
!	Setup within each of the CFD/MD the mapping to the other realm
! --- CFD ---
! create_map_cfd is an internal CFD only routine to build the mapping
! only on CFD processors to the MD regions
! ---  MD  ---
! create_map_md is an internal MD only routine to build the mapping
! only on MD processors to the CFD regions
!-----------------------------------------------------------------------------

subroutine coupler_create_map
	use mpi
	use coupler_internal_cfd
	use coupler_internal_md
    use coupler_internal_common, only : request_stop
	implicit none

	integer ierr

	if (COUPLER_REALM .eq. COUPLER_CFD) then
		call create_map_cfd
	else if (COUPLER_REALM .eq. COUPLER_MD) then
		call create_map_md
	else
		write(*,*) "Wrong COUPLER_REALM in coupler_create_map"
		call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_REALM,ierr)
	end if

    call request_stop("create_map")

end subroutine coupler_create_map

!=============================================================================
! Get MD processor topology and timestep details on all CFD processors 
! and send mesh details
!-----------------------------------------------------------------------------

subroutine coupler_cfd_init(icomm_grid, imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,&
	kmax,kmaxo,nsteps,x,y,z,dx,dz,npx,npy,npz,icoord,dt)
	use mpi
	use coupler_internal_common
	use coupler_input_data
	use coupler_internal_cfd, only : imino_ => imino, imin_ => imin, imax_ => imax, &
        jmino_ => jmin, jmin_ => jmin, jmax_ => jmax, jmaxo_ => jmaxo, kmino_ => kmino, kmin_ => kmin, &
        kmax_ => kmax, kmaxo_ => kmaxo, nsteps_ => nsteps, x_ => x, y_ => y, z_ => z, dx_ => dx, dz_ => dz, &
		npx_ => npx, npy_ => npy, npz_ => npz, icoord_ => icoord, dt_ => dt, jmax_overlap, &
		npx_md, npy_md, npz_md, nproc_md, MD_initial_cellsize

	implicit none

	integer, intent(in) :: icomm_grid,imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,kmax,kmaxo,nsteps,&
		npx,npy,npz,icoord(:,:)
	real(kind(0.d0)), intent(in) :: x(:),y(:),z(:), dx, dz, dt

	integer i, myid, source, iaux(6),ierr

    ! duplicate grid communicator for coupler use
    call mpi_comm_dup(icomm_grid,coupler_grid_comm,ierr)

	imino_ = imino; imin_ = imin;  imax_ = imax;  jmino_ = jmin
	jmin_  = jmin;  jmax_ = jmax; jmaxo_ = jmaxo; kmino_ = kmino; kmin_ = kmin; kmax_ = kmax
	kmaxo_ = kmaxo; nsteps_ = nsteps;  dx_ = dx; dz_ = dz
	npx_ = npx; npy_ = npy; npz_ = npz 

    if(kmax == kmin) then
        cfd_is_2d = .true.
    endif

	allocate(x_(size(x)),stat=ierr); x_ = x
	allocate(y_(size(y)),stat=ierr); y_ = y
	allocate(z_(size(z)),stat=ierr); z_ = z
	allocate(icoord_(3,npx*npy*npz),stat=ierr)

	x_ = x; y_ = y; z_ = z;
	icoord_ = icoord

    call mpi_comm_rank(COUPLER_REALM_COMM,myid,ierr)

    ! test if MD_init_cell size is larger than CFD cell size
    if( (MD_initial_cellsize >= x(2)-x(1) .or. MD_initial_cellsize >= y(2)-y(1)) .and. myid == 0 ) then
        write(*,*) 
        write(*,*) "********************************************************************"
        write(*,*) " WARNING ...WARNING ...WARNING ...WARNING ...WARNING ...WARNING ... "
        write(*,*) " MD initialisation cell size larger than CFD x,y cell sizes         "
        write(*,*) " MD_init_cellsize=",MD_initial_cellsize
        write(*,*) " x(2)-x(1)=",x(2)-x(1), " y(2)-y(1)=",y(2)-y(1) 
        write(*,*) "********************************************************************"
        write(*,*)
    endif
  
	! send CFD processor grid and overlap parameter

	! Note: jmax_overlap default is provided in coupler_internal_cfd
    if (cfd_coupler_input%overlap%tag == CPL) then
        jmax_overlap =  jmin + cfd_coupler_input%overlap%y_overlap
    endif

	if ( myid .eq. 0 ) then
		source=MPI_ROOT
	else
		source=MPI_PROC_NULL
	endif

	call mpi_bcast((/ npx, npy, npz, jmax_overlap /), 4, MPI_INTEGER,&
						source, COUPLER_ICOMM,ierr)

	! receive MD processor grid 
	call mpi_bcast(iaux, 3, MPI_INTEGER,0, COUPLER_ICOMM,ierr)

	npx_md = iaux(1)
	npy_md = iaux(2)
	npz_md = iaux(3)
	nproc_md = npx_md * npy_md * npz_md

	! send CFD mesh data
	call mpi_bcast((/ imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,kmax,kmaxo /), 12,&
						MPI_INTEGER, source, COUPLER_ICOMM,ierr)
	call mpi_bcast(x,size(x),mpi_double_precision,source,COUPLER_ICOMM,ierr)
	call mpi_bcast(y,size(y),mpi_double_precision,source,COUPLER_ICOMM,ierr)
	call mpi_bcast(z,size(z),mpi_double_precision,source,COUPLER_ICOMM,ierr)
	call mpi_bcast((/ dx, dz /),2,mpi_double_precision,source,COUPLER_ICOMM,ierr)

	! send CFD nsteps and dt
	call mpi_bcast(nsteps,1,mpi_integer,source,COUPLER_ICOMM,ierr)
	call mpi_bcast( (/ dt, density, MD_initial_cellsize /),3,mpi_double_precision,source,COUPLER_ICOMM,ierr)

	! write(0,*)' CFD: did exchange grid data'

end subroutine coupler_cfd_init

!=============================================================================
! Get CFD mesh and timestep details on all MD processors and send MD processor
! topology to the CFD
!-----------------------------------------------------------------------------

subroutine coupler_md_init(npxin,npyin,npzin,icoordin,icomm_grid,dtin)
	use mpi
	use coupler_internal_common
	use coupler_input_data, cfd_code_id_in => cfd_code_id
	use coupler_internal_md, b => MD_initial_cellsize, md_density => density, &
        md_steps => md_steps_per_dt_cfd
	implicit none

	integer, intent(in)	  :: npxin, npyin, npzin, icoordin(:,:),icomm_grid
	real(kind(0.d0)), intent(in) :: dtin

	integer i, ierr, source, iaux(12)
	real(kind=kind(0.d0)) ra(3)

	call mpi_comm_rank(COUPLER_REALM_COMM,myid,ierr)
    call mpi_comm_dup(icomm_grid,coupler_grid_comm,ierr)
    call mpi_comm_rank(coupler_grid_comm,myid_grid,ierr) 
    myid_grid = myid_grid + 1

	npx = npxin; npy = npyin; npz = npzin
	nproc = npx * npy * npz

	dt_MD = dtin

	allocate(icoord(3,nproc),stat=ierr)
	icoord = icoordin

	! write(0,*) 'MD exchange grid data'

	! get CFD processor grid and the number of block in j direction
	call mpi_bcast( iaux, 4, MPI_INTEGER,0, COUPLER_ICOMM,ierr)
	npx_cfd = iaux(1)
	npy_cfd = iaux(2)
	npz_cfd = iaux(3)
	jmax_overlap_cfd = iaux(4)
	nproc_cfd = npx_cfd * npy_cfd * npz_cfd

	! send MD processor grid 
	if ( myid .eq. 0 ) then
		source=MPI_ROOT
	else
		source=MPI_PROC_NULL
	endif

	call mpi_bcast((/ npx, npy, npz /), 3, MPI_INTEGER,&
		source, COUPLER_ICOMM,ierr)

	!! Test
	!	call mpi_comm_rank(MD_COMM,myid,ierr)!
	!
	!	write(0,*) 'exchange_grid_data: MD side', myid, bbox_cfd%xbb(1:2,1:npx_cfd),bbox_cfd%zbb(1:2,1:npz_cfd),&
	!	 bbox_md%xbb(1:2,1:npx_md),bbox_md%zbb(1:2,1:npz_md)    

	! CFD mesh data 
	call mpi_bcast(iaux, 12, MPI_INTEGER, 0, COUPLER_ICOMM,ierr) 
	imino = iaux(1); imin_cfd = iaux(2);  imax_cfd = iaux(3);  imaxo = iaux(4)
	jmino = iaux(5); jmin_cfd = iaux(6);  jmax_cfd = iaux(7);  jmaxo = iaux(8)
	kmino = iaux(9); kmin_cfd = iaux(10); kmax_cfd = iaux(11); kmaxo = iaux(12)
	allocate (x(imino:imaxo),y(jmino:jmaxo),z(kmino:kmaxo))

	call mpi_bcast(x,size(x),mpi_double_precision,0,COUPLER_ICOMM,ierr)
	call mpi_bcast(y,size(y),mpi_double_precision,0,COUPLER_ICOMM,ierr)
	call mpi_bcast(z,size(z),mpi_double_precision,0,COUPLER_ICOMM,ierr)
	call mpi_bcast(ra,2,mpi_double_precision,0,COUPLER_ICOMM,ierr)

	! rescale all lengths to MD units
	x = fsig * x; y = fsig * y; z = fsig * z 
	dx = fsig * ra(1); dz = fsig * ra(2)

	!write(0,*) 'MD exchange grid data: imin0, imin ...', imino, imin_cfd, imax_cfd,imaxo, &
	! x(imino),x(imin_cfd), x(imax_cfd), x(imaxo)
	!write(0,*) 'MD exchage grid data: recv dx, dz ', dx, dz, jmax_overlap_cfd

	! get CFD nsteps
	call mpi_bcast(nsteps,1,mpi_integer,0,COUPLER_ICOMM,ierr)

	! get CFD dt
	call mpi_bcast(ra,3,mpi_double_precision,0&
 					&,COUPLER_ICOMM,ierr)
	! should dt_CFD be scaled ?  see to it later
	dt_CFD     = ra(1) * FoP_time_ratio
	md_density = ra(2)
	b          = ra(3)

	! --- set the sizes of the MD domain ---
	!Fix xL_md domain size to continuum
	xL_md = x(imax_cfd) - x(imin_cfd)

    ! yL_md is adjusted to an integer number of initialisation cells in the following steps
    if (md_ly_extension_tag == CPL) then 
        DY_PURE_MD = md_ly_extension
    else 
        DY_PURE_MD =  y(jmin_cfd) - y(jmino)
    end if
    yL_md = y(jmax_overlap_cfd) - y(jmino) + DY_PURE_MD
    yL_md = real(floor(yL_md/b),kind(0.d0))*b
    DY_PURE_MD = yL_md - (y(jmax_overlap_cfd) - y(jmino))

	!write(0,*) 'MD domain etc', DY_PURE_MD,y(jmin_cfd) , yL_md, y(jmax_overlap_cfd),y(jmino)

	!Fix zL_md domain size to continuum
	zL_md = z(kmax_cfd) - z(kmin_cfd)

    if( kmin_cfd == kmax_cfd) then
        cfd_is_2d = .true.
    endif

	!write(0,*) 'MD: exchange_grid... xL_md, yL_md, zL_md',&
	! 				& myid, xL_md, yL_md, zL_md
	!write(0,*) 'MD: nsteps (CFD) ', nsteps

    ! initialise other md module variables if data is provided in coupler.in
    if (md_average_period_tag == CPL) then 
        average_period = md_average_period
    endif
    if (md_save_period_tag == CPL) then
        save_period    = md_save_period
    endif

    if (cfd_code_id_tag == CPL) then
        cfd_code_id  = cfd_code_id_in
    endif
    
    if(md_steps_per_dt_cfd_tag == CPL) then
        md_steps = md_steps_per_dt_cfd
    else 
        md_steps = int(dt_cfd/dt_MD)
    endif
    if ( md_steps <= 0 ) then 
        write(0,*) "Number of MD steps per dt interval <= 0"
        write(0,*) "Coupler will not work, quitting ..."
        call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_INIT,ierr)
    endif

end subroutine coupler_md_init

!=============================================================================
!	Adjust CFD domain size to an integer number of lattice units used by  
!	MD if sizes are given in sigma units
!-----------------------------------------------------------------------------

subroutine coupler_cfd_adjust_domain(xL, yL, zL, nx, ny, nz, density_cfd)
    use mpi
    use coupler_internal_common
    use coupler_input_data
    use coupler_internal_cfd, only : b => MD_initial_cellsize
    implicit none

    integer, optional, intent(inout) 			:: nx, ny, nz
    real(kind(0.d0)),optional, intent(inout) 	:: xL,yL,zL
    real(kind(0.d0)), optional, intent(inout) 	:: density_cfd

    ! Internal variables
    integer										:: myid, ierror, ierr
    real(kind=kind(0.d0)), pointer 				:: xyz_ptr => null()
	character(1)				   				:: direction
	logical										:: changed

    ! Local rank, useful for messeges
    call mpi_comm_rank(COUPLER_REALM_COMM,myid,ierr)

    if (density_tag == CPL) then 
        density_cfd = density
    else
        density     = density_cfd
        density_tag = CFD
    endif

	! Coupler input parameters are set in CFD 
    if ( cfd_coupler_input%domain%tag == VOID) then
        cfd_coupler_input%domain%tag = CFD
	end if

    select case (cfd_coupler_input%domain%cell_type)
    case ("FCC","Fcc","fcc")
            b = (4.d0/density)**(1.d0/3.d0)
    case default
            write(*,*) "Wrong unit cell type in coupler_cfd_adjust_domain. Stopping ... "
            ierror = COUPLER_ERROR_INIT
            call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
    end select

	! Check CFD domain and MD domain are compatible sizes to allow a
	! stable initial MD lattice structure - resize if possible or
	! stop code and demand a regeneration of grid if vary by more than 0.01
	changed = .false.
    if (present(xL)) then
        xyz_ptr => cfd_coupler_input%domain%x
        call init_length(xyz_ptr,xL,resize=.true.,direction='x',print_warning=changed)
    endif

    if (present(yL)) then
		! No need to adjust y because we can adjust DY in MD to
		! have an integer number of FCC units.
		xyz_ptr => cfd_coupler_input%domain%y
		call init_length(xyz_ptr,yL,resize=.false.,direction='y',print_warning=changed)
    endif

    if (present(zL)) then
        xyz_ptr => cfd_coupler_input%domain%z
        call init_length(xyz_ptr,zL,resize=.true.,direction='z',print_warning=changed)
    endif

	if (changed .eq. .true. .and. cfd_code_id .ne. couette_serial) then
		print*, "Regenerate Grid with corrected sizes as above"
		call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
	endif

    ! set CFD number of cells
    if (present(nx)) then
        if (cfd_coupler_input%ncells%tag == CPL ) then
            nx = cfd_coupler_input%ncells%x
        else
            if(myid == 0) then
                write(0,*)"WARNING: nx is present in coupler_cfd_adjust_domain argument list"
                write(0,*)"         but coupler_input%ncells tag is void."
                write(0,*)"         Using CFD input value!"
            endif
        endif

    endif

    if (present(ny)) then
        if (cfd_coupler_input%ncells%tag == CPL ) then
            ny = cfd_coupler_input%ncells%y
        else
            if(myid == 0) then
                write(0,*)"WARNING: ny is present in coupler_cfd_adjust_domain argument list"
                write(0,*)"         but coupler_input%ncells tag is void."
                write(0,*)"         Using CFD input value!"
            endif
        endif
    endif
        

    if (present(nz)) then
        if (cfd_coupler_input%ncells%tag == CPL ) then
            nz = cfd_coupler_input%ncells%z
        else
            if(myid == 0) then
                write(0,*)"WARNING: nz is present in coupler_cfd_adjust_domain argument list"
                write(0,*)"         but coupler_input%ncells tag is void."
                write(0,*)"         Using CFD input value!"
            endif
        endif
    endif

    ! check id CFD cell sizes are larger than 2*sigma 
    call test_cfd_cell_sizes

contains

!-----------------------------------------------------------------------------

subroutine init_length(rin,rout,resize,direction,print_warning)
	implicit none
            
    real(kind=kind(0.d0)), intent(in)    :: rin
    real(kind=kind(0.d0)), intent(inout) :: rout
    logical, intent(in)                  :: resize
	character(*),intent(in) 			 :: direction
    logical,intent(out)					 :: print_warning

    real(kind=kind(0.d0)) :: rinit  ! store the initial value of rout or rin needed for print


    print_warning=.false.

    select case (cfd_coupler_input%domain%tag)
    case (CPL)
        select case (cfd_coupler_input%domain%units ) 
        case ("CELLSIDE","CellSide","Cellside","cellside")
            rout = b * rin
        case("SIGMA", "Sigma", "sigma")
            select case (cfd_coupler_input%domain%cell_type)
            case("FCC","Fcc","fcc")
                if (resize) then 
                rinit = rin 
                rout = real(floor(rin/b),kind(0.d0))*b
                print_warning = .true.
                endif
            case default
                write(*,*) "wrong unit cell type in coupler_cfd_adjust_domain. Stopping ... "
                ierror = COUPLER_ERROR_INIT
                call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
            end select
        end select
    case (CFD) 
        if(resize) then
            rinit = rout
            rout = real(nint(rout/b),kind(0.d0))*b
            print_warning = .true.
        endif
    case default
        write(*,*) "Wrong domain tag in coupler_cfd_adjust_domain. Stopping ... "
        ierror = COUPLER_ERROR_INIT
        call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
    end select

    if (print_warning) then 
        if (myid == 0) then 
            write(*,'(3(a,/),3a,/,2(a,f20.10),/a,/,a)') &
                    "*********************************************************************",	&
                    "WARNING - this is a coupled run which resets CFD domain size         ",	&
                    " to an integer number of MD initial cells:		                      ", 	&
                    "	Domain resized in the ", direction, " direction			          ",	&
                    " inital size =", rinit, " resized ", rout,									&
                    "								                                      ",	& 
                    "*********************************************************************"   
        endif
		!If resize is insignificant then return flag print_warning as false
		if (abs(rinit-rout) .lt. 0.01) print_warning = .false.
    end if
    
end subroutine init_length

!-----------------------------------------------------------------------------


subroutine test_cfd_cell_sizes
    implicit none

    integer, pointer :: ndim => null()

    ndim => cfd_coupler_input%domain%ndim

    if ( myid == 0) then
        if ((present(xL) .and. present(nx)) .or. &
            (cfd_coupler_input%domain%tag == CPL .and. &
            cfd_coupler_input%ncells%tag == CPL)) then
            if (xL/nx < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in x direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          xL=",xL,"nx=",nx
            endif
        endif

        if ((present(yL) .and. present(ny)) .or. &
            (cfd_coupler_input%domain%tag == CPL .and. & 
            cfd_coupler_input%ncells%tag == CPL .and. ndim > 1)) then
            if (yL/ny < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in y direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          yL=",yL,"nx=",ny
            endif
        endif

        if ((present(zL) .and. present(nz)) .or. &
            (cfd_coupler_input%domain%tag == CPL .and. &
            cfd_coupler_input%ncells%tag == CPL .and. ndim > 2 )) then
            if (zL/nz < 2.0d0) then
                write(0,*)" WARNING: CFD cell size in z direction is less that 2 * sigma. Does this make sense?" 
                write(0,*)"          zL=",zL,"nx=",nz
            endif
        endif
    end if
        
end subroutine test_cfd_cell_sizes
            
end subroutine coupler_cfd_adjust_domain


!=============================================================================
! coupler_send_grid_data wrapper for 3d arrays
! see coupler_send_grid_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_send_grid_data_3d(asend,index_transpose,asend_lbound,&
    asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)
    implicit none
    real(kind=kind(0.d0)), intent(in) :: asend(:,:,:)
    integer , optional, intent(in) :: glower(3), gupper(3) 
    integer, optional, intent(in) :: asend_lbound(3)       
    integer, optional, intent(in) :: asend_grid_start(3)   
    integer, optional, intent(in) :: asend_grid_end(3)     
    integer, optional, intent(in) :: index_transpose(3)    
    logical, optional, intent(in) :: use_overlap_box(3)
    
    integer n1,n2,n3,n4
    
    n1 = 1
    n2 = size(asend,1)
    n3 = size(asend,2)
    n4 = size(asend,3)
    
    call coupler_send_grid_data_xd(asend,n1,n2,n3,n4,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)

end subroutine coupler_send_grid_data_3d

!=============================================================================
! coupler_send_grid_data wrapper for 4d arrays
! see coupler_send_grid_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_send_grid_data_4d(asend,index_transpose,asend_lbound,asend_grid_start,&
    asend_grid_end,glower,gupper,use_overlap_box)
    implicit none
    real(kind=kind(0.d0)), intent(in) :: asend(:,:,:,:)
    integer , optional, intent(in) :: glower(3), gupper(3)
    integer, optional, intent(in) :: asend_lbound(3)      
    integer, optional, intent(in) :: asend_grid_start(3)  
    integer, optional, intent(in) :: asend_grid_end(3)    
    integer, optional, intent(in) :: index_transpose(3)   
    logical, optional, intent(in) :: use_overlap_box(3)
    
    integer n1,n2,n3,n4

    n1 = size(asend,1)
    n2 = size(asend,2)
    n3 = size(asend,3)
    n4 = size(asend,4)
    
    call coupler_send_grid_data_xd(asend,n1,n2,n3,n4,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)

end subroutine coupler_send_grid_data_4d


!=============================================================================
! Send data from the local grid to the associated ranks from the other 
! realm
! The entities of this operation are as follow:
!  1) the local box of the global grid held by this rank
!  2) the data associate with a subset of local box which are held in some region of asend
!  3) the associated ranks from the other realm which overlap the same region of
!     the global grid (described in map)
!
!  In index range we have as follow:
!  the input array has the limits
!  a(n,ias:iae,jas:jae,kas:kae)
!  contains local grid data in the sector
!  a(n,iags:iage,jags:jage,kags:kage)
! 
!  which correspods to the global index region
!  igs:ige jgs:jge, kgs:kge
!
!  which must be contain in the local grid box
!  ibs:ibs, jbs:jbe, kb:kbe 
!
!  so far the transfers done correspond to ias=iags=igs=ibs, etc 
!
! Main steps of the algorithm:
!  1) Using the arguments one must find the range of the grid data in asend
!     that is exclude halos, etc
!     This can be done by passing array sector of asend 
!     e.g. asend(igrid_start:igrid_end,jgrid_start_j:grid_end,kgrid_start:kgrid_end)
!     So far only this method was tested
!  2) Use the coupler map to find which sector of the grid go to which rank from the other
!     realm
!  3) send a quick description of the data ( an array of 8 integers )
!  4) send the data if there is anything to send
!-----------------------------------------------------------------------------
subroutine coupler_send_grid_data_xd(asend,n1,n2,n3,n4,index_transpose,asend_lbound,&
    asend_grid_start,asend_grid_end,glower,gupper,use_overlap_box)
	use mpi
	use coupler_internal_cfd, only :  CFD_COMM_OVERLAP, &
        bbox_cfd, jmax_overlap, dx, dz, icoord, nlz
    use coupler_internal_md, only : bbox_md => bbox
	use coupler_internal_common
	implicit none
 
    ! asend sizes
    integer, intent(in) :: n1,n2,n3,n4
    ! array containing data distributed on the grid
    real(kind=kind(0.d0)), intent(in) :: asend(n1,n2,n3,n4)
    ! lower and upper corners of the grid to which the data of asend
    ! must be associated
    ! if not present the local grid box limits are used
    ! any negative value is discared for the default
    integer , optional, intent(in) :: glower(3), gupper(3)
    ! starting index of asend if is not (1,1,1) 
    integer, optional, intent(in) :: asend_lbound(3)
    ! starting index coordinated of grid data if is not (1,1,1)      
    integer, optional, intent(in) :: asend_grid_start(3)   
    ! upper corner of grid data stored in asend
    ! if not present the the local grid box limits are  used
    integer, optional, intent(in) :: asend_grid_end(3)
    ! specifyle the order of dimensions in asend
    ! default is ( x, y, z) but e.g. cartesian velocity
    ! in CFD solver uses (z,x,y)     
    integer, optional, intent(in) :: index_transpose(3)    ! rule of transposition for the coordinates 
    ! whether to use on MD side the extended grid of the map ( which
    ! has a halo of CFD cells arounf MD local box)
    ! this is useful to trap molecules that left 
    ! the local domain or if the CFD cells sixe is not
    ! conmensurate with sizes of MD local box
    logical, optional, intent(in) :: use_overlap_box(3)
    

    ! local indices 
    integer is,ie,js,je,ks,ke,bis,bie,bjs,bje,bks,bke,ix,iy,iz,ais,aie,ajs,aje,aks,ake,&
        mis,mie,mjs,mje,mks,mke,at(2,3),ig(2,3)

    ! auxiliaries 
    integer i,np,itag,myid, dest, req(map%n), vel_indx(8,map%n), &
        start_address(map%n+1), a_components, ierr

	real(kind=kind(0.d0)), allocatable :: vbuf(:)

	integer, save :: ncalls = 0

	! This local CFD domain is outside MD overlap zone 
	if ( map%n .eq. 0 ) return 

	call mpi_comm_rank(coupler_grid_comm,myid,ierr)

	ncalls = ncalls + 1

    ! local grid box ranges seen by this rank
    if (COUPLER_REALM == COUPLER_CFD) then 
        bis = bbox_cfd%xbb(1,icoord(1,myid+1))
        bie = bbox_cfd%xbb(2,icoord(1,myid+1))
        bjs = bbox_cfd%ybb(1,icoord(2,myid+1))
        bje = bbox_cfd%ybb(2,icoord(2,myid+1))
        bks = bbox_cfd%zbb(1,icoord(3,myid+1))
        bke = bbox_cfd%zbb(2,icoord(3,myid+1))
    else
        ! using by default the indices of CFD grid  inside the MD domain
        bis = bbox_md%is
        bie = bbox_md%ie
        bjs = bbox_md%js
        bje = bbox_md%je
        bks = bbox_md%ks
        bke = bbox_md%ke
        if(present(use_overlap_box)) then 
            if (use_overlap_box(1)) then
                bis = bbox_md%iso
                bie = bbox_md%ieo
            endif
            if(use_overlap_box(2)) then
                bjs = bbox_md%jso
                bje = bbox_md%jeo
            endif
            if(use_overlap_box(3)) then
                bks = bbox_md%kso
                bke = bbox_md%keo
            endif
       endif
    endif

    ! get the index order of x,y,z coordinates (handles transpose arrays)
    ix = 1 ; iy = 2; iz = 3
    if( present(index_transpose)) then
        ix = index_transpose(1)
        iy = index_transpose(2)
        iz = index_transpose(3)
    endif


    ! default start/end points for grid data in asend
    ! +1 sift in size is because asend grid indices start from 2
    is = bis
    ie = min(is + size(asend,ix+1) - 1,bie)
    js = bjs
    je = min(js + size(asend,iy+1) - 1,bje)
    ks = bks
    ke = min(ks + size(asend,iz+1) - 1,bke)
    ! TODO: drop a warning if asend goes over local grid box
    

    ! grid data in asend data can be mapped to a specific region
    ! of the local grid 
    ! negative values are discarded in favor of defaults
    if (present(glower))then 
        if (glower(1) > 0) is = glower(1)
        if (glower(2) > 0) js = glower(2)
        if (glower(3) > 0) ks = glower(3)
    endif

    ! sanity check is needed here

    !  upper limit exteed for grid data point stored in asend
    if (present(gupper))then 
        if (gupper(1) > 0) ie = gupper(1)
        if (gupper(2) > 0) je = gupper(2)
        if (gupper(3) > 0) ke = gupper(3)
    endif

    ! sanity check is needed here

    ! array indices in asend for the data mapped on grid
    ! +1 sift is because asend grid indices start from 2
    ais = 1
    ajs = 1
    aks = 1
    if (present(asend_lbound)) then 
        if (.not. present(asend_grid_start))then 
            write(0,*) "because the asend lower bound is not default asend_grid_start argument must be provided"
            call MPI_Abort(MPI_COMM_WORLD,COUPLER_ABORT_SEND_CFD,ierr)
        endif
        ais = asend_grid_start(ix) - asend_lbound(ix) + 1
        ajs = asend_grid_start(iy) - asend_lbound(iy) + 1
        aks = asend_grid_start(iz) - asend_lbound(iz) + 1
    endif

    aie = ais + (ie-is) 
    aje = ajs + (je-js)
    ake = aks + (ke-ks)

    ! sanity checks are needed 
    
    ! number of components at each grid point
    a_components = size(asend,dim=1)

    ! loop trough the maps and send the corresponding sections of asend
    do i = 1, map%n
        dest = map%rank_list(i)

        ! indices for the grid sector to be sent to this dest
        mis = max(is, map%domains(1,i)) 
        mie = min(ie, map%domains(2,i))
        mjs = max(js, map%domains(3,i))
        mje = min(je, map%domains(4,i))
        mks = max(ks, map%domains(5,i)) 
        mke = min(ke, map%domains(6,i))	

        np = a_components * (mie - mis + 1) * (mje - mjs + 1) * (mke - mks + 1)

        !write(0,*) 'coupler send', myid, mis,mie,mjs,mje,mks,mke,a_components,np

        if ( np > 0) then 
            if(allocated(vbuf)) deallocate(vbuf)
            
            allocate(vbuf(np))
            
            ! location in asend of the domain to be send
            ig(1,ix) = mis - is + ais 
            ig(2,ix) = mie - is + ais
            ig(1,iy) = mjs - js + ajs
            ig(2,iy) = mje - js + ajs
            ig(1,iz) = mks - ks + aks
            ig(2,iz) = mke - ks + aks
 
            !write(0,*) ' coupler send a*, ig ...', ncalls, myid, ais, aie,ajs,aje,aks,ake,ig 

            vbuf(1:np) = reshape(asend(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)), (/ np /) )
        endif
        ! Attention ncall could go over max tag value for long runs!!
        itag = mod( ncalls, MPI_TAG_UB)
        ! send the info about the data to come
        call mpi_send((/np,a_components,mis,mie,mjs,mje,mks,mke/),8,MPI_INTEGER,&
            dest, itag, COUPLER_ICOMM, ierr)
        ! send data only if there is anything to send
        if (np > 0) then 
			!vbuf = 1.234567891011121314151617d0
			!print'(2a,2i8,4f25.16)', 'ICP send data',code_name(COUPLER_REALM), myid, & 
			!							 size(vbuf), maxval(vbuf),minval(vbuf),sum(vbuf),vbuf(10)
            call mpi_send(vbuf, np, MPI_DOUBLE_PRECISION, dest, itag, COUPLER_ICOMM, ierr)
        endif
    enddo

end subroutine coupler_send_grid_data_xd

!=============================================================================
! coupler_recv_grid_data wrapper for 3d arrays
! see coupler_recv_grid_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_recv_grid_data_3d(arecv,index_transpose,asend_lbound,asend_grid_start,&
    asend_grid_end,glower,gupper,accumulate,pbc)
    implicit none

    real(kind(0.d0)), intent(inout) :: arecv(:,:,:)

    integer, optional, intent(in) :: asend_lbound(3)
    integer, optional, intent(in) :: asend_grid_start(3) 
    integer, optional, intent(in) :: asend_grid_end(3)   
    integer, optional, intent(in) :: index_transpose(3)  
    integer, optional, intent(in) :: glower(3), gupper(3)
    logical, optional, intent(in) :: accumulate          
    integer, optional, intent(in) :: pbc                 
                                                         

    integer n1, n2, n3, n4

    n1 = 1 
    n2 = size(arecv,1)
    n3 = size(arecv,2)
    n4 = size(arecv,3)

    call coupler_recv_grid_data_xd(arecv,n1,n2,n3,n4,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,accumulate,pbc)

end subroutine coupler_recv_grid_data_3d

!=============================================================================
! coupler_recv_grid_data wrapper for 4d arrays
! see coupler_recv_grid_data_xd for input description
!-----------------------------------------------------------------------------
subroutine coupler_recv_grid_data_4d(arecv,index_transpose,asend_lbound,&
    asend_grid_start,asend_grid_end,glower,gupper,accumulate,pbc)

    implicit none

    real(kind(0.d0)), intent(inout) :: arecv(:,:,:,:)
    integer, optional, intent(in) :: asend_lbound(3)     
    integer, optional, intent(in) :: asend_grid_start(3) 
    integer, optional, intent(in) :: asend_grid_end(3)   
    integer, optional, intent(in) :: index_transpose(3)  
    integer, optional, intent(in) :: glower(3), gupper(3)
    logical, optional, intent(in) :: accumulate          
    integer, optional, intent(in) :: pbc                 
                                                         

    integer n1, n2, n3, n4

    n1 = size(arecv,1)
    n2 = size(arecv,2)
    n3 = size(arecv,3)
    n4 = size(arecv,4)

    call coupler_recv_grid_data_xd(arecv,n1,n2,n3,n4,index_transpose,asend_lbound,&
        asend_grid_start,asend_grid_end,glower,gupper,accumulate,pbc)

end subroutine coupler_recv_grid_data_4d

!=============================================================================
! Receive data from to local grid from the associated ranks from the other 
! realm
! The entities of this operation are as follow:
!  1) the local box of the global grid held by this rank
!  2) the grid data region to be held in some region of arecv
!  3) the associated ranks from the other realm which overlap the same region of
!     the global grid
!
! Main steps of the algorithm:
!  1) Using the arguments find the range of the grid data in arecv
!     that is: exclude halos, etc
!     This can be done by passing array sector of arecv that is meant to hold grid data 
!     e.g. arecv(igrid_start:igrid_end,jgrid_start_j:grid_end,kgrid_start:kgrid_end)
!     So far only this method was tested
!  2) Use the coupler map to find where the arriving data should land in a temporary array
!  3) copy the whole lot, after the receives completes, to the corresponding location
!     of arecv
!
!  Notes: 
!  1) This subroutine can perform accumulation of data, if one chooses to 
!     send sum of cell quantities and respective numbers. Could be useful to 
!     collect contribution from molecule leaving local MD box.
!  2) Periodic boundary condition can be applied, useful for averages with
!     staggered boxes in which the end boxes collect only form half box
!     at the moment this works only with (z,x,y) array in z and x directions
!   
!    
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine coupler_recv_grid_data_xd(arecv,n1,n2,n3,n4,index_transpose,a_lbound,&
    									a_grid_start,a_grid_end,glower,gupper,accumulate,pbc)
    use mpi
	use coupler_internal_cfd, only : bbox_cfd, jmax_overlap, dx, dz, icoord, nlz
    use coupler_internal_md, only :  bbox_md => bbox
	use coupler_internal_common
    implicit none

    ! arecv sizes
    integer, intent(in) :: n1, n2, n3, n4
    ! array that recieves grid distributed data 
    real(kind(0.d0)), intent(inout) :: arecv(n1,n2,n3,n4)
    ! lower and upper corners of the grid to which the data of asend
    ! must be associated
    ! if not present the local grid box limits are used
    ! any negative value is discared for the default
    integer, optional, intent(in) :: glower(3), gupper(3)
    ! starting index of asend if is not (1,1,1)
    integer, optional, intent(in) :: a_lbound(3) 
    ! starting index coordinated of grid data if is not (1,1,1)    
    integer, optional, intent(in) :: a_grid_start(3)
    ! upper corner of grid data stored in asend
    ! if not present the the local grid box limits are used
    integer, optional, intent(in) :: a_grid_end(3)  
    ! specifyle the order of dimensions in arecv
    ! default is ( x, y, z) but e.g. cartesian velocity
    ! in CFD solver uses (z,x,y)
    integer, optional, intent(in) :: index_transpose(3) 
    ! if present add the receiving value if they land in the same
    ! grid position
    logical, optional, intent(in) :: accumulate   
    ! if present applies periodic boundary condition
    ! posible values 1 for x direction, 3 for z direction
    integer, optional, intent(in) :: pbc                 
                                                         
    
    ! local indices 
    integer is,ie,js,je,ks,ke,bis,bie,bjs,bje,bks,bke,ix,iy,iz,ais,aie,ajs,aje,aks,ake,&
        at(2,3),ig(2,3)
    ! auxiliaries 
    integer i,j,k,ii,jj,kk,np,itag,myid, source,req(map%n), vel_indx(8,map%n), &
        p1s,p1e,p2s,p2e,p3s,p3e,pt1s,pt1e,pt2s,pt2e,pt3s,pt3e, bgt(2,3),& 
        start_address(map%n+1), status(MPI_STATUS_SIZE,map%n),ierr
    real(kind(0.d0)), allocatable ::  vbuf(:), atmp(:,:,:,:)
    integer, save :: ncalls = 0

	! This local CFD domain is outside MD overlap zone 
	if ( map%n .eq. 0 ) return 
 
	ncalls = ncalls + 1

	call mpi_comm_rank(coupler_grid_comm,myid,ierr)

    ! local grid box
    if (COUPLER_REALM == COUPLER_CFD) then 
         bis = bbox_cfd%xbb(1,icoord(1,myid+1))
         bie = bbox_cfd%xbb(2,icoord(1,myid+1))
         bjs = bbox_cfd%ybb(1,icoord(2,myid+1))
         bje = bbox_cfd%ybb(2,icoord(2,myid+1))
         bks = bbox_cfd%zbb(1,icoord(3,myid+1))
         bke = bbox_cfd%zbb(2,icoord(3,myid+1))
     else 
		! using by default the indices of CFD grid  inside the MD domain
        bis = bbox_md%is
        bie = bbox_md%ie
        bjs = bbox_md%js
        bje = bbox_md%je
        bks = bbox_md%ks
        bke = bbox_md%ke
    endif

    ! get the idices in x,y,z direction from transpose array
    ix = 1 ; iy = 2; iz = 3
    if( present(index_transpose)) then
        ix = index_transpose(1)
        iy = index_transpose(2)
        iz = index_transpose(3)
    endif
 
    
    ! default start/end points on grid for data in asend
	is = bis
    ie = bie
	js = bjs
    je = bje
	ks = bks
    ke = bke 

    ! grid data in asend data can be mapped to a specific region
    ! of the local grid 
    ! negative values are discarded in favor of defaults
    if (present(glower))then 
        if (glower(1) > 0) is = glower(1)
        if (glower(2) > 0) js = glower(2)
        if (glower(3) > 0) ks = glower(3)
    endif

    ! sanity check is needed here

    !  upper limit exteed for grid data point stored in asend
    if (present(gupper))then 
        if (gupper(1) > 0) ie = gupper(1)
        if (gupper(2) > 0) je = gupper(2)
        if (gupper(3) > 0) ke = gupper(3)
    endif


    ! sanity check is needed here

    ! put the mapped block limits in a transposed array
    bgt(1,ix) = is
    bgt(2,ix) = ie
    bgt(1,iy) = js
    bgt(2,iy) = je
    bgt(1,iz) = ks
    bgt(2,iz) = ke

    ! array indices in arecv for the data mapped on grid
    ! +1 sift is because asend grid indices start from 2
    ais = 1
    ajs = 1
    aks = 1
    if (present(a_lbound)) then 
        if (.not. present(a_grid_start))then 
            write(0,*) "because the arecv lower bound is not default as_grid_start argument must be provided"
            call MPI_Abort(MPI_COMM_WORLD,COUPLER_ABORT_SEND_CFD,ierr)
        endif
        ais = a_grid_start(1) - a_lbound(1) + 1
        ajs = a_grid_start(2) - a_lbound(2) + 1
        aks = a_grid_start(3) - a_lbound(3) + 1
    endif

    aie = min(ais + (ie-is),size(arecv,ix+1)) 
    aje = min(ajs + (je-js),size(arecv,iy+1))
    ake = min(aks + (ke-ks),size(arecv,iz+1))

    ! sanity checks are needed 

    
    ! store the transposition for grid boundaries limits
    at(1,ix) = ais
    at(2,ix) = aie
    at(1,iy) = ajs
    at(2,iy) = aje
    at(1,iz) = aks
    at(2,iz) = ake

    ! first get info about what's comming
    do i = 1, map%n
        source =  map%rank_list(i)
        itag = mod( ncalls, MPI_TAG_UB)
        call mpi_irecv(vel_indx(1,i),8,MPI_Integer,source,itag,&
            			COUPLER_ICOMM,req(i),ierr)
    enddo

    call mpi_waitall(map%n, req, status, ierr)

    !write(0,*) 'coupler recv vel_indx ', vel_indx

    np = 0
    do i = 1, map%n
        np = np + max(vel_indx(1,i),0)
    enddo

	allocate(vbuf(np),stat=ierr)
 
    start_address(1) = 1 
    do i = 1, map%n
        source =  map%rank_list(i) 
        
        np = vel_indx(1,i)
        
        start_address(i+1) = start_address(i) + np

        if (np > 0) then 
            
 			! Attention ncall could go over max tag value for long runs!!
			itag = mod(ncalls, MPI_TAG_UB)
            call mpi_irecv(vbuf(start_address(i)),np,MPI_DOUBLE_PRECISION,source,itag,&
                				COUPLER_ICOMM,req(i),ierr)
        else 
            req(i) = MPI_REQUEST_NULL
        endif

    enddo

    call mpi_waitall(map%n, req, status, ierr)
	!print'(2a,2i8,4f25.16)', 'ICP recv data',code_name(COUPLER_REALM),myid, & 
	!							 size(vbuf), maxval(vbuf),minval(vbuf),sum(vbuf),vbuf(10)
    !	write(0,*) 'MD getCFD vel wait',  myid, id, i, source, ierr

    ! allocate atmp corresponding to reunion of mapped grid
    ! This must be thought as landing area for data coming from the other realm
    ig(1,ix) = minval(vel_indx(3,:))
    ig(2,ix) = maxval(vel_indx(4,:))
    ig(1,iy) = minval(vel_indx(5,:))
    ig(2,iy) = maxval(vel_indx(6,:))
    ig(1,iz) = minval(vel_indx(7,:))
    ig(2,iz) = maxval(vel_indx(8,:))

    !write(0,*) 'coupler recv atmp sizes ', vel_indx(2,1),ig

    allocate(atmp(vel_indx(2,1), ig(1,1):ig(2,1), ig(1,2):ig(2,2), ig(1,3):ig(2,3)))
    atmp(:,:,:,:) = 0.d0

    do i=1,map%n
        
        if(vel_indx(1,i) <= 0) cycle
        
        ! get form grid the indices to the array ones
        ig(1,ix) = vel_indx(3,i)
        ig(2,ix) = vel_indx(4,i)
        ig(1,iy) = vel_indx(5,i)
        ig(2,iy) = vel_indx(6,i)
        ig(1,iz) = vel_indx(7,i)
        ig(2,iz) = vel_indx(8,i)
        
        if (present(accumulate)) then
            !write(0,*) 'coupler recv ig', ig
            atmp(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)) = &
                atmp(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)) &
                + reshape(vbuf(start_address(i):start_address(i+1)-1), &
                (/ vel_indx(2,i), ig(2,1)-ig(1,1)+1, ig(2,2)-ig(1,2)+1, ig(2,3)-ig(1,3)+1 /))
        else
            atmp(:,ig(1,1):ig(2,1),ig(1,2):ig(2,2),ig(1,3):ig(2,3)) = &
                reshape(vbuf(start_address(i):start_address(i+1)-1), &
                (/ vel_indx(2,i), ig(2,1)-ig(1,1)+1, ig(2,2)-ig(1,2)+1, ig(2,3)-ig(1,3)+1 /))
        endif

        ! call MPI_Error_string(status(MPI_ERROR,i), err_string, len(err_string), ierr);
	    !  write(0,*) 'MD getCFD vel err, myid, i ', myid, i, trim(err_string) 
        !call mpi_get_count(status(1,i),mpi_double_precision,ib,ierr)
                                ! write(0,*) 'MD recv ', myid, id, i, ib, ' DP'
    enddo

    ! transfer data from the landing area to arecv 
    ! if the flag accumulate is present normalise
    ! data with the value stored on the last line of the first column
    ! and apply periodic boundary conditions

    ! indices for overlapping window
    ! usually the landing arera and arecv coincide
    ! some test and warning should be done here
    p1s = max(lbound(atmp,2) - bgt(1,1) + at(1,1),at(1,1))
    p1e = min(ubound(atmp,2) - bgt(1,1) + at(1,1),at(2,1))
    p2s = max(lbound(atmp,3) - bgt(1,2) + at(1,2),at(1,2))
    p2e = min(ubound(atmp,3) - bgt(1,2) + at(1,2),at(2,2))
    p3s = max(lbound(atmp,4) - bgt(1,3) + at(1,3),at(1,3))
    p3e = min(ubound(atmp,4) - bgt(1,3) + at(1,3),at(2,3))

    ! needs revision, I'm not sure that it works in general !!!
    pt1s = p1s + bgt(1,1) - at(1,1)
    pt1e = p1e + bgt(1,1) - at(1,1)
    pt2s = p2s + bgt(1,2) - at(1,2)
    pt2e = p2e + bgt(1,2) - at(1,2)
    pt3s = p3s + bgt(1,3) - at(1,3)
    pt3e = p3e + bgt(1,3) - at(1,3)

    !write(0,*)' p1s ...', myid, p1s,p1e,p2s,p2e,p3s,p3e, pt1s, pt1e, pt2s,pt2e,pt3s,pt3e

    if(present(accumulate)) then
        if ( present(pbc)) then 
            call set_pbc(pbc)
        endif
        
        np = ubound(atmp,1)
        do k = p3s , p3e
 		do j = p2s , p2e
 		do i = p1s , p1e
			ii = i + bgt(1,1) - at(1,1)
			jj = j + bgt(1,2) - at(1,2)
			kk = k + bgt(1,3) - at(1,3)
 			if (atmp(np, ii,jj,kk) > 0.d0) then
				arecv(:,i,j,k) = atmp(1:np-1,ii,jj,kk)/atmp(np,ii,jj,kk)
			else
				arecv(:,i,j,k) = 0.d0
			endif
		end do
		end do
        end do
    else 
        arecv(:,p1s:p1e,p2s:p2e,p3s:p3e) = atmp(:,pt1s:pt1e,pt2s:pt2e,pt3s:pt3e) 
    endif

	!print'(2a,2i8,4f25.16)', 'ICP trfr data',code_name(COUPLER_REALM),myid, & 
	!						size(arecv), maxval(arecv),minval(arecv),sum(arecv),arecv(1,p1s,p2s,p3s)
            
contains

!-----------------------------------------------------------------------------
! Sets periodic boundary conditions
! adds together the sums of the first and last bins on the global grid (x and z directions)
! and puts sum in both boxes
! Note:
! attention coord array are shifted with 1 in CFD and MD
! to get correct result in MPI function such as mpi_cart_rank, subtract 1 
!-----------------------------------------------------------------------------
subroutine set_pbc(pbc)
    use coupler_internal_cfd, only : npx
    implicit none
    integer, intent(in) :: pbc

    real(kind(0.d0)), allocatable, dimension(:,:,:) :: x1, x2
    integer dest, ip(3), ierr, status(MPI_STATUS_SIZE)

    ! PBC works only for 3,2,1 coordinate transposition

    if ( .not. (ix == 2 .and. iy == 3 .and. iz == 1)) then
        write(0,*) " Periodic boundary condition not implements for this layout of data"
        return
    endif

    ! Sanity check atmp must touch the boundary of the global grid in zin x directions
    ! I guess that the PBC are apllied at the ends of global grid, check with CFD people
    select case(pbc)
    case(3)
        ! first array dimension which corresponds to z direction vor velocity arrays 
        allocate(x1(size(atmp,1),size(atmp,3),size(atmp,4)),x2(size(atmp,1),size(atmp,3),size(atmp,4)))
        x1 =  atmp(:,aks, :,:)
        x2 =  atmp(:,ake, :,:)
        atmp(:, aks, :, :) =   atmp(:, aks, :, :) + x2
        atmp(:, ake, :, :) =   atmp(:, ake, :, :) + x1
    case(1)
        ! second array dimension which correponds to x direction
        allocate(x1(size(atmp,1),size(atmp,2),size(atmp,4)),x2(size(atmp,1),size(atmp,2),size(atmp,4)))
        if ( npx == 1 )then 
            ! no MPI communication needed  
            !write(0,*) 'coupler internal cfd pbc', npx
            x1 =  atmp(:, :,lbound(atmp,3), :)                 
            x2 =  atmp(:, :,ubound(atmp,3), :)
            atmp(:, :, lbound(atmp,3), :) =  atmp(:, :, lbound(atmp,3), :) + x2
            atmp(:, :, ubound(atmp,3), :) =  atmp(:, :, ubound(atmp,3), :) + x1  
        else 
            call mpi_comm_rank(coupler_grid_comm,myid,ierr)
            ip = icoord(: ,myid+1)
            if (ip(1) == 1 .and. ip(2) == 1) then 
                x1 =  atmp(:, :,lbound(atmp,3),:)
                call mpi_cart_rank(coupler_grid_comm, (/ npx-1, 0, ip(3)-1 /), dest, ierr)
                call mpi_sendrecv(x1,size(x1),MPI_DOUBLE_PRECISION,dest,1,x2,size(x2),MPI_DOUBLE_PRECISION,&
                    dest,1,coupler_grid_comm,status,ierr)
                atmp(:, :, lbound(atmp,3), :) =  atmp(:, :, lbound(atmp,3), :) + x2
            else if (ip(1) == npx .and. ip(2) == 1) then 
                x2 =  atmp(:, :,ubound(atmp,3), :)
                call mpi_cart_rank(coupler_grid_comm, (/ 0, 0, ip(3) - 1 /), dest, ierr)
                call mpi_sendrecv(x2,size(x2),MPI_DOUBLE_PRECISION,dest,1,x1,size(x1),MPI_DOUBLE_PRECISION,&
                    dest,1,coupler_grid_comm,status,ierr)
                atmp(:, :, ubound(atmp,3), :) =  atmp(:, :, ubound(atmp,3), :) + x1
            endif
        endif
    end select
end subroutine set_pbc


end subroutine coupler_recv_grid_data_xd

!============================================================================
!
! utility functions and subroutines that extract parameters from internal modules 
!
!-----------------------------------------------------------------------------    

!===============================================================================
! return to the caller coupler parameters from cfd realm
!-------------------------------------------------------------------------------
subroutine coupler_cfd_get(jmax_overlap,jmin,jmino)
    use coupler_internal_cfd, jmax_overlap_ => jmax_overlap,jmin_=>jmin,jmino_=>jmino
    implicit none

    integer,optional,intent(out) :: jmax_overlap,jmin,jmino

    if(present(jmax_overlap)) then
        jmax_overlap = jmax_overlap_
    endif

    if(present(jmin))then 
        jmin = jmin_
    end if

    if(present(jmino))then
        jmino = jmino_
    end if

end subroutine coupler_cfd_get


!===============================================================================
! return to the caller coupler parameters from MD realm
!-------------------------------------------------------------------------------
subroutine coupler_md_get(xL_md,yL_md,zL_md, MD_initial_cellsize, top_dy, &
    overlap_with_continuum_force, overlap_with_top_cfd, ymin_continuum_force, &
    ymax_continuum_force, xmin_cfd_grid, xmax_cfd_grid, zmin_cfd_grid, zmax_cfd_grid, &
    dx_cfd, dz_cfd, cfd_box)
    use mpi
	use coupler_internal_md, only : xL_md_ =>  xL_md, yL_md_ =>  yL_md, zL_md_ => zL_md, &
		b => MD_initial_cellsize, x, y, z,j => jmax_overlap_cfd, dx, dz, bbox, half_domain_lengths, &
        cfd_box_ => cfd_box
	implicit none

	real(kind(0.d0)), optional, intent(out) :: xL_md, yL_md, zL_md, MD_initial_cellsize, top_dy,&
        ymin_continuum_force,ymax_continuum_force, xmin_cfd_grid, xmax_cfd_grid, &
		 zmin_cfd_grid, zmax_cfd_grid, dx_cfd, dz_cfd 
    logical, optional, intent(out) :: overlap_with_continuum_force, overlap_with_top_cfd
    type(cfd_grid_info), optional, intent(out) :: cfd_box

    integer ierr

    if (present(xL_md)) then
		xL_md = xL_md_
	endif

	if (present(yL_md)) then 
		yL_md = yL_md_
	endif

	if (present(zL_md)) then 
		zL_md = zL_md_
	endif

	if (present(MD_initial_cellsize)) then 
		MD_initial_cellsize = b
	endif

	if (present(top_dy)) then
    	top_dy = y(j) - y(j-1)
	end if

    if(present(overlap_with_continuum_force)) then
        ! check if the countinuum force domain is contained in one 
        ! domain along y direction
        if ( (y(j - 2) < bbox%bb(1,2) .and. &
              y(j - 1) > bbox%bb(1,2)) .or. &
             (y(j - 2) < bbox%bb(2,2) .and. &
              y(j - 1) > bbox%bb(2,2))) then
            write(0,*) " the region in which the continuum constraint force is applied "
            write(0,*) " spans over two domains. This case is not programmed, please investigate"
            call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_CONTINUUM_FORCE,ierr)
        endif
        
        if ( y(j - 1) <  bbox%bb(2,2) .and. y(j - 2) >= bbox%bb(1,2) ) then
            overlap_with_continuum_force = .true.
        else   
            overlap_with_continuum_force = .false.
        endif

    endif

     if(present(overlap_with_top_cfd)) then
        ! check if the MD domain overlaps with the top of CFD grid (along y)
        ! the MD constrain force is applyied top layer of cfd cells
        if ( y(j - 1) < bbox%bb(2,2) .and. y(j - 1) >= bbox%bb(1,2) ) then
            overlap_with_top_cfd = .true.
        else   
            overlap_with_top_cfd = .false.
        endif

    endif
 
     if(present(ymin_continuum_force)) then
         ymin_continuum_force = y(j - 2) - bbox%bb(1,2) - half_domain_lengths(2)
     endif
         
     if(present(ymax_continuum_force)) then
         ymax_continuum_force = y(j - 1) - bbox%bb(1,2) - half_domain_lengths(2)
     endif 

     if(present(xmin_cfd_grid)) then
         xmin_cfd_grid = x(bbox%is) - bbox%bb(1,1) - half_domain_lengths(1)
     endif

     if(present(xmax_cfd_grid)) then
         xmax_cfd_grid = x(bbox%ie) - bbox%bb(1,1) - half_domain_lengths(1)
     endif

     if(present(zmin_cfd_grid)) then
         zmin_cfd_grid = z(bbox%ks) - bbox%bb(1,3) - half_domain_lengths(3)
     endif

     if(present(zmax_cfd_grid)) then
         zmax_cfd_grid = z(bbox%ke) - bbox%bb(1,3) - half_domain_lengths(3)
     endif
 
     if (present(dx_cfd)) then 
         dx_cfd = dx
     endif

     if (present(dz_cfd)) then 
         dz_cfd = dz
     endif

     if (present(cfd_box))then
         cfd_box%gimin = cfd_box_%gimin
         cfd_box%gimax = cfd_box_%gimax
         cfd_box%gjmin = cfd_box_%gjmin
         cfd_box%gjmax = cfd_box_%gjmax
         cfd_box%gkmin = cfd_box_%gkmin
         cfd_box%gkmax = cfd_box_%gkmax
         cfd_box%imin  = cfd_box_%imin
         cfd_box%imax  = cfd_box_%imax
         cfd_box%jmin  = cfd_box_%jmin
         cfd_box%jmax  = cfd_box_%jmax
         cfd_box%kmin  = cfd_box_%kmin
         cfd_box%kmax  = cfd_box_%kmax
         cfd_box%imino  = cfd_box_%imino
         cfd_box%imaxo  = cfd_box_%imaxo
         cfd_box%jmino  = cfd_box_%jmino
         cfd_box%jmaxo  = cfd_box_%jmaxo
         cfd_box%kmino  = cfd_box_%kmino
         cfd_box%kmaxo  = cfd_box_%kmaxo
         cfd_box%xmin  = cfd_box_%xmin
         cfd_box%xmax  = cfd_box_%xmax
         cfd_box%dx    = cfd_box_%dx
         cfd_box%ymin  = cfd_box_%ymin
         cfd_box%ymax  = cfd_box_%ymax
         allocate( cfd_box%y(cfd_box%jmin:cfd_box%jmax))
         cfd_box%y(cfd_box%jmin:cfd_box%jmax) = cfd_box_%y(cfd_box%jmin:cfd_box%jmax)
         cfd_box%zmin = cfd_box_%zmin
         cfd_box%zmax = cfd_box_%zmax
         cfd_box%dz = cfd_box_%dz
     end if

end subroutine coupler_md_get

!-----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------- 
function coupler_md_get_save_period() result(p)
	use coupler_internal_md, only : save_period
	implicit none

	integer p

	p = save_period
end function coupler_md_get_save_period
!-----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------- 

function coupler_md_get_average_period() result(p)
	use coupler_internal_md, only : average_period
	implicit none

	integer p

	p = average_period
end function coupler_md_get_average_period
!-----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------- 

function coupler_md_get_md_steps_per_dt_cfd() result(p)
	use coupler_internal_md, only : md_steps_per_dt_cfd
	implicit none

	integer p

	p = md_steps_per_dt_cfd
end function coupler_md_get_md_steps_per_dt_cfd

!-----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------- 
function coupler_md_get_nsteps() result(n)
	 use coupler_internal_md, only : nsteps
	 implicit none 

	 integer n

	 n = nsteps
end function coupler_md_get_nsteps

!-----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------- 
function coupler_md_get_dt_cfd() result(t)
	 use coupler_internal_md, only : dt_CFD  
	 implicit none

	 real(kind=kind(0.d0)) t

	 t = dt_CFD
end function coupler_md_get_dt_cfd

!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------  
subroutine coupler_md_set(zL_md)
	use coupler_internal_md, only : zL => zL_md, z, dz, kmino, kmin_cfd, kmax_cfd, kmaxo
	implicit none
	
	real(kind(0.d0)), optional, intent(in) :: zL_md

	if ( present(zL_md) )then
		zL = zL_md
		if (allocated(z)) then 
			deallocate(z)
		endif
		! values need for the overlap map
		allocate(z(2))
		z(:)       = (/ 0.d0, zL_md /)
		dz	 = zL_md
		kmino      = 1
		kmin_cfd   = 1
		kmax_cfd   = 2
		kmaxo      = 2
		
	endif
end subroutine coupler_md_set

!-----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------- 

function coupler_md_get_density() result(r)
	use coupler_internal_md, only : density
	implicit none

	real(kind(0.d0)) r

	r = density

end function coupler_md_get_density

!==============================================================================
!
!------------------------------------------------------------------------------
function coupler_md_get_cfd_id() result(r)
    use coupler_internal_md, only : cfd_code_id
    implicit none

    integer r

    r = cfd_code_id

end function coupler_md_get_cfd_id




! Below are deprecated subroutine but still useful for
! testing or reference
!=============================================================================
! test subroutine
!-----------------------------------------------------------------------------

subroutine coupler_uc_average_test(np,r,v,lwrite)
	use coupler_internal_common
	use coupler_internal_md, only : nlx, nlz, nly, bbox, jmino, jmin => jmin_cfd,&
					global_r, x, dx, y, z, dz
	 
	implicit none

	integer, intent(in) :: np
	real(kind=kind(0.d0)), intent(in) :: r(:,:),v(:,:)
	logical, intent(in) :: lwrite

	integer ib, kb, jb, ip, myid, jbuff, ierr
	real(kind=kind(0.d0)) rd(3), ymin, ymax, dy
	real(kind=kind(0.d0)),allocatable, save :: uc_bin(:,:,:,:)
	logical,save :: firsttime=.true.
    save jbuff

	call mpi_comm_rank(COUPLER_REALM_COMM,myid,ierr)

	if(firsttime)then
		firsttime = .false.
       
        ! how many cells to collect along y
        jbuff = nly+1 !jmino - 1 to see the boundary effects


		allocate(uc_bin(2,nlz-1,nlx-1,jbuff))
		uc_bin = 0.d0

		if (myid .eq. 0) then 
			open(45, file="md_vel.txt",position="rewind")
			write(45,*)'# dx,dy,dz ', dx,y(jmin+1)-y(jmin),dz
            write(45,*)'# nlx, nly, nlz', nlx,nly,nlz
			close(45)
		endif
	endif

	if (lwrite) then
		call write_data
		return
	endif


	dy = y(jmin+1) - y(jmin)
	ymin = y(jmino) -  dy
	ymax = y(bbox%je - 1)

	!	write(0,*)'MD uc test', np, dy, ymin,ymax

	do ip = 1, np
		! using global particle coordinates
		rd(:)=r(ip,:)
		rd(:) = global_r(rd)

		if ( rd(2) < ymin .or. rd(2) > ymax ) then
			! molecule outside the boundary layer
			cycle
		endif

		ib = ceiling((rd(1) - x(bbox%is)) / dx) + 0      ! staggered !!!
		kb = ceiling((rd(3) - z(bbox%ks)) / dz)       ! the last z row unused
		jb = ceiling((rd(2) - ymin    )   / dy)

		if ( ib > 0 .and. ib < nlx  .and. &
			kb > 0 .and. kb < nlz  ) then 
			!  this particle are in this ranks domain
			uc_bin(1,kb,ib,jb) = uc_bin(1,kb,ib,jb) + v(ip,1)
			uc_bin(2,kb,ib,jb) = uc_bin(2,kb,ib,jb) + 1.d0 
		else 
			!				       write(0,*) 'MD uc_average, outside domain rd', rd, ' bbox%bb ', bbox 
		endif
	end do

	! debug   
	!			 do i = 1, size(uc_bin,dim=2)
	!			  write(0, '(a,I4,64F7.1)') 'MD myid uc_bin(2,..',myid,uc_bin(2,1,:)
	!			 enddo



	!			write(0,*) 'MD uc sum in boxes', myid
	!			do i = 1, size(uc_bin,dim=2)
	!				write(0, '(a,I4,64E11.4)') 'MD myid uc_bin(1,..',myid, uc_bin(1,1,:)
	!			enddo
	! send it to CFD	

contains 

!=============================================================================
! Write velocities from MD domain
!-----------------------------------------------------------------------------

	subroutine write_data
		use mpi
		use coupler_internal_md, only : nproc, imin_cfd, imax_cfd, kmin_cfd, kmax_cfd
		implicit none

		integer i, ibuff(2,2,0:nproc-1), ntot, nrecv, sa(nproc),req(nproc-1),  &
			ierr
		real(kind(0.d0)),allocatable :: buff(:,:,:,:),buff_recv(:)

		if(nproc > 1) then

			! works only for parallel decomposition in x and y direction
			call mpi_gather((/bbox%iso,bbox%ieo,bbox%kso,bbox%keo/),4,MPI_INTEGER,&
				ibuff,4,MPI_INTEGER,0,COUPLER_REALM_COMM,ierr)

			!	       write(0,*) "MD write test data", myid, ibuff

			if (myid .eq. 0) then

				! the local bit first
				allocate(buff(2,kmax_cfd-kmin_cfd,imin_cfd:imax_cfd-1,jbuff))

				buff = 0.d0
				buff(:,ibuff(1,2,0):ibuff(2,2,0)-1,ibuff(1,1,0):ibuff(2,1,0)-1,:) = &
					buff(:,ibuff(1,2,0):ibuff(2,2,0)-1,ibuff(1,1,0):ibuff(2,1,0)-1,:) + &
					uc_bin(:,1:nlz-1,1:nlx-1,:)


				ntot = 0
				do i=1,nproc-1
					ntot = ntot + 2*(ibuff(2,2,i)-ibuff(1,2,i))*(ibuff(2,1,i)-ibuff(1,1,i))*jbuff
				enddo

				allocate(buff_recv(ntot))
				buff_recv(ntot) = 0.d0

				sa(1)=1

				do i=1,nproc-1

					nrecv = 2*(ibuff(2,2,i)-ibuff(1,2,i))*(ibuff(2,1,i)-ibuff(1,1,i))*jbuff
					sa(i+1) = sa(i) + nrecv

					call mpi_irecv(buff_recv(sa(i)),nrecv,MPI_DOUBLE_PRECISION,&
						i,1,COUPLER_REALM_COMM,req(i),ierr)

				enddo

				call mpi_waitall(nproc-1,req,MPI_STATUSES_IGNORE,ierr)

				do i =1, nproc-1
					buff(:,ibuff(1,2,i):ibuff(2,2,i)-1,ibuff(1,1,i):ibuff(2,1,i)-1,:) = &
						buff(:,ibuff(1,2,i):ibuff(2,2,i)-1,ibuff(1,1,i):ibuff(2,1,i)-1,:) + &
						reshape(buff_recv(sa(i):sa(i+1)-1), (/ 2,ibuff(2,2,i)-ibuff(1,2,i),ibuff(2,1,i)-ibuff(1,1,i),jbuff /))
				enddo
			else
				call mpi_send(uc_bin,size(uc_bin),MPI_DOUBLE_PRECISION,0,1,COUPLER_REALM_COMM,ierr)
			endif
		endif

		if (nproc > 1) then
			if (myid .eq. 0 ) then
				open(45,file="md_vel.txt",position="append")
				do i = 1,jbuff
					write(45, '(100(E12.4,1x))') sum(buff(:,:,:,i),dim=2) 
				enddo
				write(45, '(1x/1x)')
				close(45)
			endif
		else
			open(45,file="md_vel.txt",position="append")
			do i = 1,jbuff
				write(45, '(100(E12.4,1x))') sum(uc_bin(:,:,:,i),dim=2)
			enddo
			write(45, '(1x/1x)')
			close(45)
		endif

		uc_bin = 0.d0

	end subroutine write_data

end subroutine coupler_uc_average_test


! Keep these subroutines for a while, useful references 

!=============================================================================
!	Apply force to prevent molecules leaving domain using form suggested by 
!	Nie, Chen, E and Robbins (2004)
!-----------------------------------------------------------------------------

!!$subroutine coupler_md_boundary_forces(np,pressure,r,a)
!!$	use coupler_internal_md, only : bbox, jmax_overlap_cfd, halfdomain => half_domain_lengths, y
!!$	implicit none 
!!$
!!$	integer, intent(in)		  :: np
!!$	real(kind=kind(0.d0)), intent(inout) :: a(:,:)
!!$	real(kind=kind(0.d0)), intent(in)    :: pressure, r(:,:)
!!$
!!$	! locals
!!$	real(kind=kind(0.d0)), parameter :: eps = 0.d-2 ! avoid singular forces for molecules that 
!!$													! happend to be very close to y(jmax_overlap)
!!$	integer n
!!$	real(kind=kind(0.d0)) p, yc, y2, y3
!!$
!!$
!!$	! Initial pressure is -ve - this line prevents problems	
!!$	p = pressure
!!$	if(pressure <= 0 ) then 
!!$		p= 1.d0
!!$	endif
!!$
!!$	y2 = y(jmax_overlap_cfd -1)
!!$	y3 = y(jmax_overlap_cfd) 
!!$
!!$	do n = 1, np
!!$		! get the global value of y coordinate	
!!$		yc  =  r(n,2) + halfdomain(2) + bbox%bb(1,2)
!!$		if (yc  < y3 .and. yc >= y2 ) then
!!$			a(n,2)= a(n,2) - p*(yc-y2)/(1.d0-(yc-y2)/(y3-y2)+eps)
!!$		endif
!!$	enddo
!!$
!!$
!!$end subroutine coupler_md_boundary_forces

!=============================================================================
! Apply force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
!-----------------------------------------------------------------------------

!!$subroutine coupler_md_apply_continuum_forces(np,r,v,a,iter)
!!$	use coupler_internal_common
!!$	use coupler_internal_md, only : cfd_box_sum, halfdomain => half_domain_lengths, &
!!$					x, y, z, dx, dz, global_r, jmin => jmin_cfd, &
!!$					nlx, nly, nlz, dt_CFD,bbox, get_CFDvel, &
!!$					jmax_overlap => jmax_overlap_cfd, myid
!!$	implicit none
!!$
!!$	real(kind=kind(0.d0)), dimension(:,:), intent(in) :: r,v
!!$	real(kind=kind(0.d0)), dimension(:,:), intent(inout) :: a 
!!$	integer, intent(in) :: np,iter  ! iteration step, it assumes that each MD average
!!$	! start from iter = 1
!!$
!!$	type(cfd_box_sum) :: box_average(bbox%ie - bbox%is,1, bbox%ke - bbox%ks)
!!$	integer j, ib, jb, kb, nib, njb, nkb, ip, np_overlap, jb_constrain
!!$	integer list(4,np)
!!$	real(kind=kind(0.d0)) inv_dtCFD, rd(3)
!!$	integer :: ncalls = 0
!!$
!!$	! run through the particle, check if they are in the overlap region
!!$	! find the CFD box to which the particle belongs	      
!!$	! attention to the particle that have left the domain boundaries 
!!$
!!$        ! This work is done only by the MD ranks that cover the constraint region
!!$        ! At the moment use only the second layer from the top of CFD cell 
!!$        if (  jmax_overlap - 2 < bbox%js .or. jmax_overlap - 2 >= bbox%je ) return
!!$
!!$	! number of CFD cells in each direction
!!$	nib = bbox%ie - bbox%is
!!$	njb = bbox%je - bbox%js
!!$	nkb = bbox%ke - bbox%ks
!!$
!!$	! vel_fromCFD cell index from which continum constrain is applied
!!$	jb_constrain =   njb - 1 ! the second row of cells from the top
!!$
!!$	if (iter .eq. 1) then
!!$		! get the previous value of CFD velocities
!!$		call  get_CFDvel
!!$	endif
!!$
!!$	! at first CFD step we don't have two values to extrapolate CFD velocities, set inv_dtCFD=0
!!$	if (ncalls .eq. 0) then
!!$		inv_dtCFD = 0.0
!!$	else
!!$		inv_dtCFD = 1.0/dt_CFD
!!$	endif
!!$	ncalls = ncalls + 1
!!$
!!$	np_overlap = 0 ! number of particles in overlapping reg
!!$
!!$	do kb = 1, ubound(box_average,dim=3)
!!$		do jb = 1, ubound(box_average,dim=2)
!!$			do ib = 1, ubound(box_average,dim=1)
!!$				box_average(ib,jb,kb)%np   = 0
!!$				box_average(ib,jb,kb)%v(:) = 0.0d0
!!$				box_average(ib,jb,kb)%a(:) = 0.0d0
!!$			enddo
!!$		enddo
!!$	enddo
!!$
!!$	do ip = 1, np
!!$
!!$		! we need global MD coordinates to check if the particle is in the extended box.
!!$		! bbox%bb(:,:) are ok for they were used to build the MD domains
!!$		rd(:) = r(ip,:)
!!$		rd(:) = global_r(rd)
!!$
!!$		! struggling with the bottom boundary, below average boxes but with particles
!!$		!  for the moment let's work with 1 layer of MD blocks in 1 D
!!$		if ( rd(2) <= y(jmax_overlap-2) .or.   rd(2) >= y(jmax_overlap-1) ) then
!!$			cycle 
!!$		else 
!!$                        jb = 1
!!$                        ! version to be analized later
!!$			! non uniform grid in j direction		
!!$			!			 do j =jmin+1, jmax_overlap
!!$			!				if( rd(2) <= y(j) ) then 
!!$			!					!this is my cell index, exit
!!$			!					jb = j - jmin
!!$			!					exit
!!$			!				endif
!!$			!			  enddo
!!$
!!$		endif
!!$
!!$		! get the CFD cell coordinates   
!!$		if (rd(1) < x(bbox%is) .or. rd(1) >= x(bbox%ie)) then
!!$			! this particle has left the domanin
!!$			!				write(0,*) 'particle lost in x direction'
!!$			cycle
!!$		else 
!!$			ib = ceiling((rd(1) -  x(bbox%is))/ dx)
!!$		endif
!!$
!!$		if (rd(3) < z(bbox%ks) .or. rd(3) >= z(bbox%ke) ) then
!!$			! this particle has left the domanin
!!$			!				write(0,*) 'particle lost in z direction'
!!$			cycle
!!$		else
!!$			kb = ceiling( (rd(3) - z(bbox%ks) ) / dz) 
!!$		endif
!!$
!!$		np_overlap = np_overlap + 1
!!$		list(1:4, np_overlap) = (/ ip, ib, jb, kb /)
!!$
!!$		box_average(ib,jb,kb)%np   =  box_average(ib,jb,kb)%np + 1
!!$		box_average(ib,jb,kb)%v(:) =  box_average(ib,jb,kb)%v(:) + v(ip,:)
!!$		box_average(ib,jb,kb)%a(:) =  box_average(ib,jb,kb)%a(:) + a(ip,:)
!!$
!!$	enddo
!!$
!!$	! here we should have the cell coordinates for the particle ip which is 
!!$	! in the overlap region
!!$	! one has to treat separatley the particle that have left the domain
!!$	! compute the average force for each bin
!!$
!!$	!write(0,*)'MD before average over bin. np_overlap', np_overlap
!!$
!!$	call average_over_bin
!!$
!!$	!write(0,*) 'MD: end simulation_apply_continuum_forces', myid
!!$
!!$contains
!!$
!!$!=============================================================================
!!$! Get velocity from CFD and apply force to molecule
!!$!-----------------------------------------------------------------------------
!!$
!!$	subroutine average_over_bin
!!$		use coupler_internal_md, only : dt_MD, myid, itm1, itm2, vel_fromCFD 
!!$		implicit none
!!$
!!$		integer ib, jb, kb, i, ip, n
!!$		real(kind=kind(0.d0)) alpha(3), u_cfd_t_plus_dt(3), inv_dtMD, acfd
!!$
!!$
!!$		! set the continnum constraints for the particle in the bin
!!$		! speed extrapolation 
!!$		! add all up
!!$		inv_dtMD =1.d0/dt_MD
!!$
!!$		!write(0,'(a,I7,2E12.4)') "MD continuum np, vel_fromCFD1 : ", np_overlap, &
!!$		!						  maxval(a(list(1,1:np_overlap),:)), &
!!$		!						  minval(a(list(1,1:np_overlap),:))
!!$
!!$		do i = 1, np_overlap  
!!$			ip = list(1,i)
!!$			ib = list(2,i)
!!$			jb = list(3,i)
!!$			kb = list(4,i)
!!$
!!$			n = box_average(ib,jb,kb)%np
!!$
!!$			!write(0,'(a,4I4,14E12.4)') "MD continuum force", ib,jb,kb,n,box_average(ib,jb,kb)%v(:), &
!!$			!	box_average(ib,jb,kb)%a(:),v(ip,:),a(ip,:),inv_dtMD,inv_dtCFD
!!$
!!$			if ( n .eq. 0 ) cycle
!!$
!!$			! using the following exptrapolation formula for continuum velocity
!!$			! y = (y2-y1)/(x2-x1) * (x-x2) +y2
!!$			alpha(1) = inv_dtCFD*(vel_fromCFD(1,ib,jb_constrain,kb,itm1) - &
!!$				vel_fromCFD(1,ib,jb_constrain,kb,itm2))
!!$
!!$			u_cfd_t_plus_dt(1) = alpha(1) * (iter + 1)*dt_MD + vel_fromCFD(1,ib,jb_constrain,kb,itm1) 
!!$
!!$			acfd =  - box_average(ib,jb,kb)%a(1) / n - inv_dtMD * & 
!!$				( box_average(ib,jb,kb)%v(1) / n - u_cfd_t_plus_dt(1) )
!!$			a(ip,1) = a(ip,1) + acfd
!!$
!!$			!	write(0,'(a,4I4,15E12.4)') "MD continuum force 2", ib,jb,kb,n, &
!!$			!	 alpha(1),u_cfd_t_plus_dt(1),vel_fromCFD(1,ib,jb+jb_offset,kb,itm1),&
!!$			!	 vel_fromCFD(1,ib,jb+jb_offset,kb,itm2),&
!!$			!	 a(ip,1),acfd, r(ip,2) 
!!$
!!$		enddo
!!$
!!$
!!$		!	write(400+10*ncalls+myid,'(a,I7,2E12.4)') "MD continuum np, vel_fromCFD 2: ", np_overlap, &
!!$		!						   maxval(a(list(1,1:np_overlap),:)), &
!!$		!						   minval(a(list(1,1:np_overlap),:))
!!$		!	write(400+10*ncalls+myid,'(a,2E12.4)')" inv_dtCFD, inv_dtMD ", inv_dtCFD, inv_dtMD
!!$		!	do kb=1,nkb
!!$		!	do jb=1,njb
!!$		!	do ib=1,nib
!!$		!		write(400+10*ncalls+myid,'(12E12.4,I7)') vel_fromCFD(:,ib,jb,kb,1), vel_fromCFD(:,ib,jb,kb,2),&
!!$		!					      box_average(ib,jb,kb)%v(:), box_average(ib,jb,kb)%a(:),&
!!$		!					      box_average(ib,jb,kb)%np
!!$		!	enddo
!!$		!	enddo
!!$		!	enddo
!!$
!!$	end subroutine average_over_bin
!!$
!!$end subroutine coupler_md_apply_continuum_forces

!============================================================================
! Send CFD velocities via the intercommunicator to the MD processors
! which are specified by the topological maps setup in create_map_cfd 
!
! Algoritm brief description
! 1) CFD caller passes the array of data to be sent, the index ranges
!    that contains data ( so far assumes 1:ngz in z ditection) and 
!    the start i,j start on the global grid
!
! 1') When a local quatity starts frm grif index > imin (as uc)
!     one has to exted the array index to the left to cover the 
!     whole global grid. ( Assumes that some halo are present in the
!     local quantity which hold boundary condition)
!     To be discussed further !!!
!     Later: this goes better to the socket where one can addapt
!     the local array at wish
!
! 2) Use map to figure with domain of local indices go to what MD core
!
! 3) Send the starting global indices (x,y) and the extend in all direction
!    followed by data 
!-----------------------------------------------------------------------------
!!$subroutine coupler_cfd_send_velocity(xc,i1_x,i2_x,i1b_x,j1_x,j2_x,j1b_x)
!!$	use mpi
!!$	use coupler_internal_cfd, only : md_map, nlx, nly, nlz, CFD_COMM_OVERLAP, &
!!$					 bbox_cfd, jmax_overlap, dx, dz, icoord, nlz, coupler_grid_comm
!!$	use coupler_internal_common
!!$	implicit none
!!$
!!$    integer , intent(in) :: i1_x,i2_x,i1b_x,j1_x,j2_x,j1b_x
!!$	real(kind=kind(0.d0)), intent(in) :: xc(0:,0:,0:)
!!$
!!$	real(kind=kind(0.d0)), allocatable :: vbuf(:)
!!$	integer i, j,k, is, ie, js, je, ks, ke,i1,i2,j1,j2,gsi,gsj, gsi_md, gsj_md, &
!!$        gsk_md, min_i, min_j, min_k, np, myid, &
!!$		itag, dest, type, req(md_map%n), ierr
!!$	integer status(MPI_STATUS_SIZE,md_map%n)
!!$	integer, save :: ncalls = 0
!!$	character(len=20) fname
!!$
!!$
!!$	! This local CFD domain is outside MD overlap zone 
!!$	if ( md_map%n .eq. 0 ) return 
!!$
!!$	call mpi_comm_rank(coupler_grid_comm,myid,ierr)
!!$
!!$	ncalls = ncalls + 1
!!$	!		write(0,*) "CFD, send_CFDvel: ", myid
!!$	is = bbox_cfd%xbb(1,icoord(1,myid+1))
!!$	!ie = bbox_cfd%xbb(2,icoord(1,myid+1))
!!$	js = bbox_cfd%ybb(1,icoord(2,myid+1))
!!$	!je = bbox_cfd%ybb(2,icoord(2,myid+1))
!!$	ks = bbox_cfd%zbb(1,icoord(3,myid+1))
!!$	!ke = bbox_cfd%zbb(2,icoord(3,myid+1))
!!$
!!$	min_i = minval(md_map%domains(1,:))
!!$	min_j = minval(md_map%domains(3,:))
!!$	min_k = minval(md_map%domains(5,:))
!!$ 
!!$
!!$    ! if i1b_x is not is we extend to the margin ( to catch the boundary condition)
!!$    ! however one has to check that we don't get out of xc index range
!!$    ! check it later 
!!$
!!$    i1  = i1_x
!!$    gsi = i1b_x
!!$    if (i1b_x > is ) then
!!$        if (i1u_x - (i1b_x-is) >= 0) then
!!$            i1 =i1_x - (i1b_x-is)
!!$            gsi = is
!!$            write(0,*) 'rank ', myid, ' xc exteded to the left ', is,i1b_x, i1_x,i2_x
!!$        else
!!$            write(0,*) 'ux index cannot be extended to the left'
!!$        endif
!!$    endif
!!$
!!$    i2 = i2_x
!!$
!!$     if (i1b_x + (i2_x-i1_x) < ie ) then
!!$        if ( i1_x + ie - i1b_x  <= ubound(ux,2)) then
!!$            i2  = i1_x + (ie - i1b_x)
!!$            write(0,*) 'rank ', myid, ' ux exteded to the right ', is,i1b_x, i1_x,i2_x
!!$        else
!!$            write(0,*) 'ux index cannot be extended to the right'
!!$        endif
!!$    endif
!!$
!!$    j1 = j1_x
!!$    gsj= j1b_x
!!$    j2 = min(j2_x, j1_x+jmax_overlap-j1b_x)
!!$
!!$    
!!$    !vaux(1:nzl, gsi:gsi+i2_x-i1x, gsj:gsj+j2_x-j1_x) = ux(1:nlz,i1_x:i2_x,j1_x:j2_x)
!!$
!!$    do i = 1, md_map%n
!!$        dest = md_map%rank_list(i)
!!$        is = max(i1,i1 + md_map%domains(1,i) - gsi) 
!!$        ie = min(i2,i1 + md_map%domains(2,i) - gsi)
!!$        js = max(j1,j1 + md_map%domains(3,i) - gsj)
!!$        je = min(j2,j1 + md_map%domains(4,i) - gsj)
!!$        ks = 1 + md_map%domains(5,i) - 1 
!!$        ke = 1 + md_map%domains(6,i) - 1				
!!$        np = (ie - is + 1) * (je - js + 1) * (ke - ks + 1)
!!$   
!!$        gsi_md = max(gsi,md_map%domains(1,i))
!!$        gsj_md = max(gsj,md_map%domains(3,i))
!!$        gsk_md = md_map%domains(5,i)
!!$   
!!$        if(allocated(vbuf)) deallocate(vbuf)
!!$
!!$        !np=size(xc(1:nlz,i1:i2,j1:j2))
!!$        allocate(vbuf(np))
!!$        
!!$        vbuf(1:np) = reshape(xc(ks:ke,is:ie,js:je), (/ np /) )
!!$
!!$        ! Attention ncall could go over max tag value for long runs!!
!!$        itag = mod( ncalls, MPI_TAG_UB)
!!$        call mpi_send((/gsk_md,ke-ks+1,gsi_md,ie-is+1,gsj_md,je-js+1/),6,MPI_INTEGER,&
!!$            dest, itag, COUPLER_ICOMM, ierr)
!!$        call mpi_send(vbuf, np, MPI_DOUBLE_PRECISION, dest, itag, COUPLER_ICOMM, ierr)
!!$        !				write(0,*) 'CFD sendCFD vel ', myid, j, i, itag,dest,np, is,ie,js,je,ks,ke,ierr	
!!$    enddo
!!$
!!$
!!$	! write(0,*) 'CFD send_CFDvel finish' , myid, ncalls
!!$	! This barrier does not work as this function si called inside an jblock .eq. 1
!!$	! condition. Is it really necessary?
!!$	!		call mpi_barrier(CFD_COMM, ierr)
!!$
!!$end subroutine coupler_cfd_send_velocity
!!$
!!$!============================================================================
!!$! Recieves MD velocities via the intercommunicator from the MD processors
!!$! which are specified by the topological maps setup in create_map_cfd 
!!$! 
!!$!  Algorithm description:
!!$!  1)MD data are received in a landing area which is indexed according to the 
!!$!    coupler map.
!!$!
!!$!  2)To this landing array periodic boundary condition are applied if needed.
!!$!  
!!$!  3)Using the local indexes and the global correspondent provided in the argumentlist 
!!$!    data from landing area are copied to the argument arrays
!!$!  
!!$!  Comment: this is somewhat wasteful in memory but very clear and easy to extend to 
!!$!  other quantities 
!!$!-----------------------------------------------------------------------------
!!$subroutine coupler_cfd_get_velocity(uc,iuls,iule,iugs,&
!!$    vc,ivls,ivle,ivgs,wc,iwls,iwle,iwgs)
!!$	use mpi
!!$	use coupler_internal_common
!!$	use coupler_internal_cfd, only : nlx, nlz, recv_vel_MD, bbox_cfd, icoord
!!$	implicit none
!!$    !
!!$    ! local start, end and global start for uc,vc,wc in x direction
!!$    integer, intent(in) :: iuls,iule,iugs,ivls,ivle,ivgs,iwls,iwle,iwgs
!!$    ! cartesian velocities
!!$	real(kind=kind(0.d0)),dimension(0:,0:,0:),intent(out) :: uc, vc, wc 
!!$
!!$	integer kuls,kule,kugs,kvls,kvle,kvgs,kwls,kwle,kwgs,pbc
!!$	integer  ierr
!!$
!!$    ! assume that in z direction the indexes are fixed
!!$    ! they should be moved in the argument list 
!!$    kuls = 1; kule = nlz-1; kugs = 1
!!$    kvls = 1; kvle = nlz-1; kvgs = 1
!!$    kwls = 1; kwle = nlz  ; kwgs = 1
!!$
!!$    !set periodic boundary condition
!!$    if (staggered_averages(1)) then 
!!$        pbc = 2
!!$    else
!!$        pbc = 0
!!$    endif
!!$	call  recv_vel_MD(uc, kuls,kule,kugs,iuls,iule,iugs, 0, 0, pbc)
!!$    pbc = 0
!!$	call  recv_vel_MD(vc, kvls,kvle,kvgs,ivls,ivle,ivgs, 0, 1, pbc) 
!!$    if (.not. cfd_is_2d) then
!!$        if (staggered_averages(3)) then 
!!$            pbc = 1
!!$        else
!!$            pbc = 0
!!$        endif
!!$        call  recv_vel_MD(wc,kwls,kwle,kwgs,iwls,iwle,iwgs, 0, 0, pbc)
!!$    endif
!!$
!!$end subroutine coupler_cfd_get_velocity
!!$
!!$!============================================================================
!!$!
!!$! Receives comtinuum velocity and copies into vel_cfd if overlap is true
!!$!
!!$!-----------------------------------------------------------------------------
!!$subroutine coupler_md_get_CFDvel(overlap, vel_cfd)
!!$    use coupler_internal_md, only : itm1, itm2, vel_fromCFD, bbox, nlx,nly,nlz
!!$    implicit none
!!$
!!$    logical, intent(in) :: overlap
!!$    real(kind(0.d0)), intent(out) :: vel_cfd(:,:,:,:,:)
!!$    ! locals    
!!$    integer jb_constrain
!!$    
!!$
!!$    !call get_CFDvel
!!$
!!$    call coupler_recv_grid_data(arecv,index_transpose,asend_lbound,asend_grid_start,asend_grid_end,&
!!$     glower,gupper,accumulate,pbc)
!!$
!!$    if (overlap) then 
!!$        
!!$        jb_constrain = bbox%je - bbox%js - 1
!!$        
!!$        call 
!!$
!!$        vel_cfd(:,:,1,:,1) = vel_fromCFD(:,1:nlx-1,jb_constrain,1:nlz-1,itm1)
!!$        vel_cfd(:,:,1,:,2) = vel_fromCFD(:,1:nlx-1,jb_constrain,1:nlz-1,itm2)
!!$
!!$    endif
!!$
!!$
!!$end subroutine coupler_md_get_CFDvel

end module coupler
