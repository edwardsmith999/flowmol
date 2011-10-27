module coupler_cfd_global_data
        implicit none
        save 

        logical, parameter :: use_coupling = .true.

        integer CFD_COMM                  ! CFD communicator
        character(len=*), parameter :: code_name = "CFD"
        integer, parameter :: CFD = 1, MD = 2
        integer CFD_MD_ICOMM, CFD_BDY_MD_ICOMM


! CFD data brough in by the constructor
        integer imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,&
                kmino,kmin,kmax,kmaxo,npx,npy,npz, nsteps
        integer, allocatable :: icoord(:,:)

        real(kind(0.d0)), allocatable :: x(:), y(:), z(:)
        real(kind(0.d0)) dx, dz, dt

! bounding box type that holds on each processor the boundaries of the
! CFD grid sectors allocate to every processor
        type bbox
                integer, allocatable :: xbb(:,:), ybb(:,:), zbb(:,:)
        end type bbox

        type(bbox),target :: bbox_cfd, bbox_md

        integer :: jmax_overlap = 4 ! maximum j index ( in y direction) which MD has to cover, on top of that it has to cover y(0):y(1) domain and a bit of room
                                    ! below

! MD data

        integer npx_md, npy_md, npz_md, nproc_md
        integer, allocatable :: icoord_md(:,:)

! velocity average data
real, allocatable ::  vel_fromMD(:,:,:,:) 


contains

        subroutine create_communicators(comm_out)
                use mpi
                implicit none

                integer, intent(out) :: comm_out
                
                integer :: color = CFD
                
                integer ierr, myid, myid_comm, myid_comm_max,&
                        iaux(2), jaux(2), remote_leader, comm, comm_size


                call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)

                comm = MPI_COMM_NULL

                call mpi_comm_split(MPI_COMM_WORLD,color,myid,comm,ierr)

                CFD_COMM = comm
                comm_out = comm

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

                write(0,*) 'did (inter)communicators ', code_name, myid

        end subroutine create_communicators


end module coupler_cfd_global_data
