module coupler_cfd_global_data
        use data_export, only : npx_cfd => npx, npy_cfd => npy, npz_cfd => npz
        implicit none
        save 

        logical, parameter :: use_coupling = .true.

        integer CFD_COMM                  ! CFD communicator
        character(len=*), parameter :: code_name = "CFD"
        integer, parameter :: CFD = 1, MD = 2
        integer CFD_MD_ICOMM, CFD_BDY_MD_ICOMM


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

        integer icomm_xz, icoord_xz(2,npx_cfd*npz_cfd) ! obvious ?

        integer, allocatable :: ibmin_md(:), ibmax_md(:),&
                jbmin_md(:), jbmax_md(:), kbmin_md(:),&
                kbmax_md(:)

        integer novr ! number of overlaped blocks between CFD and MD 
! MD data

        integer npx_md, npy_md, npz_md, nproc_md
        integer, allocatable :: icoord_md(:,:)

! velocity average data
real, allocatable ::  vel_fromMD(:,:,:) 


contains

        subroutine create_communicators
                use mpi
                use data, only : file_dir
                implicit none
                
                integer :: color = CFD
                
                integer ierr, myid, myid_comm, myid_comm_max,&
                        iaux(2), jaux(2), remote_leader, comm, comm_size

! set input output directory for this code
                file_dir="cfd_data/"

                call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)

                comm = MPI_COMM_NULL

                call mpi_comm_split(MPI_COMM_WORLD,color,myid,comm,ierr)

                CFD_COMM = comm

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


end module coupler_cfd_global_data
