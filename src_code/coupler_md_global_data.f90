module coupler_md_global_data
        use computational_constants_MD, only : npx_md => npx, npy_md => npy, npz_md => npz
        implicit none
        save
! basic coupling data
        logical, parameter :: use_coupling = .true.

        integer MD_COMM                  ! MD communicator
        character(len=*), parameter :: code_name = "MD"
        integer, parameter :: CFD = 1, MD = 2

        integer CFD_MD_ICOMM, CFD_BDY_MD_ICOMM

! data needed by setup amd communication

!
!  contains processor coordinates, and overlaping regions
!
        type map_data
                integer rank        ! the rank of the processor holding the overlaping domain
                                    ! these ranks start from 1 !!!
                integer coord(2)    ! coordinates ( which cart ?) of the above rank
                integer ib_range(2) ! min-max block coordinates in the overlaping range i dimension
                integer kb_range(2) ! min-max k dimension
                integer dp2d_type   ! mpi subarray type use to communicate data 
        end type map_data
! integer, allocatable :: map_overlap(:,:)
        type(map_data),allocatable :: map_overlap(:)

! CFD grid 
        integer, allocatable :: ibmin_1(:), ibmax_1(:),&
                jbmin_1(:), jbmax_1(:), kbmin_1(:),&
                kbmax_1(:)

        integer novr ! number of overlaped blocks between CFD and MD 

! CFD grid boundary per MD MPI rank 
        integer ibmin_md(npx_md), ibmax_md(npx_md),jbmin_md(npx_md), jbmax_md(npx_md), kbmin_md(npz_md), kbmax_md(npz_md)

! CFD mesh data
        integer imino, imaxo, jmino, jmaxo, kmino, kmaxo
        real(kind=kind(0.d0)), allocatable :: x(:), y(:), z(:)
        real(kind=kind(0.d0)) dx, dz
        
! nteps from CFD
        integer nsteps

! Y coordinate  of FD boundary in MD domains 
        real :: Y_boundary_fraction = 0.5    ! prototype alternative

        real xL_md, yL_md, zL_md ! macroscopic sizes of MD domain. needed?
        real :: fsig=20.0  !Ratio betwen macroscopic unit lenght and molecular unit 

! shift_x, shift_z are shifts for the average box from the center of FD cell
! default is that the average box overlaps the FD cell
        real :: shift_x=0.0, shift_z=0.0 ! shift of velocity average box




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

                write(0,*) 'did communicators ', code_name, myid

        end subroutine create_communicators


end module coupler_md_global_data
