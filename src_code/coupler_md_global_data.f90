module coupler_md_global_data
        use computational_constants_MD, only : npx_md => npx, npy_md => npy, npz_md => npz
        implicit none
        save
! basic coupling data
        logical, parameter :: use_coupling = .true.

        integer MD_COMM                  ! MD communicator
        character(len=*), parameter :: code_name = "MD"
        integer, parameter :: CFD = 1, MD = 2

        integer CFD_MD_ICOMM, CFD_BDY_MD_ICOMM ! intercommunicates

! data needed by setup and communication modules

!
!  contains processor coordinates, and overlaping regions
!
        type map_data
                integer rank        ! the rank of the processor holding the overlaping domain
                                    ! these ranks start from 1 !!!
                integer coord(2)    ! coordinates ( which cart ?) of the above rank
                integer ib_range(2) ! min-max block coordinates in the overlaping range i dimension
                integer kb_range(2) ! min-max k dimension
                integer dp2d_type   ! mpi subarray type use to communicate data ( cell data) 
                integer plane_type  ! subarray to communicate plane values
        end type map_data
! integer, allocatable :: map_overlap(:,:)
        type(map_data),allocatable :: map_overlap(:)

! bounding box type that holds on each processor the boundaries of the
! CFD grid sectors allocate to every processor
!        type bbox
!                integer :: is, ie, js, je, ks, ke ! start and end CFD indices 
!                real(kind=kind(0.d0)) :: bb(2,3)  ! start and end of this domain in three directions
!        end type bbox

!        type(bbox), target :: bbox_cfd, bbox_md

! CFD mesh data
        integer imino, imin_cfd, imax_cfd,imaxo, jmino, jmin_cfd, jmax_cfd, jmaxo,&
                kmino, kmin_cfd, kmax_cfd, kmaxo
        real(kind=kind(0.d0)), allocatable, target :: x(:), y(:), z(:)
        real(kind=kind(0.d0)) dx, dz
        
! nteps from CFD
        integer nsteps
! average period for velocities        
        integer :: average_period = 5

! Y coordinate  of FD boundary in MD domains 
!        real :: Y_boundary_fraction = 0.5    ! prototype alternative
!        integer :: Y_domain_thickness = 4    ! MD global size along  
                                             ! y axis in (y(2)-y(1) units


        real(kind=kind(0.d0)) :: FoP_time_ratio = 1.0    ! time ratio dt_CFD/dt_MD; to be fixed later
        real(kind=kind(0.d0))    xL_md, yL_md, zL_md ! macroscopic sizes of MD domain. needed?
        real(kind=kind(0.d0)) :: fsig=1.0  !Ratio betwen macroscopic unit lenght and molecular unit 

! shift_x, shift_z are shifts for the average box from the center of FD cell
! default is that the average box overlaps the FD cell
        real :: shift_x=0.0, shift_z=0.0 ! shift of velocity average box

! array for velocities from CFD, last dimension holds time indices 
        real(kind=kind(0.d0)), allocatable :: vel_fromCFD(:,:,:,:,:)
        integer itm1,itm2


contains
        subroutine create_communicators
                use mpi
                use computational_constants_MD, only : file_dir
                implicit none

                integer :: color = MD

                integer ierr, myid, myid_comm, myid_comm_max,&
                        iaux(2), jaux(2), remote_leader, comm, comm_size

! set input output directory for this code
                 file_dir="md_data/"

                call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)

                comm = MPI_COMM_NULL

                call mpi_comm_split(MPI_COMM_WORLD,color,myid,comm,ierr)

                MD_COMM = comm

! create inter-communicators

! Get the mpi_comm_world ranks that hold the largest ranks in cfd_comm and md_comm
                call mpi_comm_rank(comm,myid_comm,ierr)
                call mpi_comm_size(comm,comm_size,ierr)

                iaux(:) = -1
                jaux(:) = -1
                if ( myid_comm == comm_size - 1) then
                        iaux(color) = myid
                endif
                call mpi_allreduce( iaux ,jaux, 2, MPI_INTEGER, MPI_MAX, &
                        MPI_COMM_WORLD, ierr)  

                select case (color)
                case (CFD)
                        remote_leader = jaux(MD)
                case (MD)
                        remote_leader = jaux(CFD)
                end select

                call mpi_intercomm_create(comm, comm_size - 1, MPI_COMM_WORLD,&
                        remote_leader, 1, CFD_MD_ICOMM, ierr)

                write(0,*) 'did(inter) communicators ', code_name, myid

        end subroutine create_communicators


end module coupler_md_global_data
