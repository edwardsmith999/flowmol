module coupler_cfd_setup
        use coupler_cfd_global_data
        implicit none
        save

contains

        subroutine exchange_grid_data
                use mpi
                use coupler_cfd_global_data
                use data_export, only : ibmin_1, ibmax_1, kbmin_1, kbmax_1,&
                        imino, imaxo, jmino, jmaxo, kmino, kmaxo
                use mesh_export, only : x, y, z, dx, dz
                use simulation, only : nsteps
                implicit none

                 integer i, ierr, myid, source, ngrid(3), ny_md, nblock(3), &
                         iaux(6)

                 write(0,*) 'CFD exchange grid data'
! send CFD processor grid
                call mpi_comm_rank(CFD_COMM,myid,ierr)
                if ( myid == 0 ) then
                        source=MPI_ROOT
                else
                        source=MPI_PROC_NULL
                endif

                call mpi_bcast((/ npx_cfd, npy_cfd, npz_cfd /), 3, MPI_INTEGER,&
                        source, CFD_MD_ICOMM,ierr)

! send MD processor grid 
                call mpi_bcast(iaux, 3, MPI_INTEGER,&
                        0, CFD_MD_ICOMM,ierr)

                npx_md = iaux(1)
                npy_md = iaux(2)
                npz_md = iaux(3)
                nproc_md = npx_md * npy_md * npz_md

! send CFD grid info                
                call mpi_bcast(ibmin_1,npx_cfd,mpi_integer,source,&
                        CFD_MD_ICOMM,ierr)
                call mpi_bcast(ibmax_1,npx_cfd,mpi_integer,source,&
                        CFD_MD_ICOMM,ierr)
                call mpi_bcast(kbmin_1,npz_cfd,mpi_integer,source,&
                        CFD_MD_ICOMM,ierr)
                call mpi_bcast(kbmax_1,npz_cfd,mpi_integer,source,&
                        CFD_MD_ICOMM,ierr)
! send CFD mesh data

                call mpi_bcast((/ imino,imaxo,jmino,jmaxo,kmino,kmaxo /), 6,&
                        MPI_INTEGER, source, CFD_MD_ICOMM,ierr)
                call mpi_bcast(x,size(x),mpi_double_precision,source,CFD_MD_ICOMM,ierr)
                call mpi_bcast(y,size(y),mpi_double_precision,source,CFD_MD_ICOMM,ierr)
                call mpi_bcast(z,size(z),mpi_double_precision,source,CFD_MD_ICOMM,ierr)
                call mpi_bcast((/ dx, dz /),2,mpi_double_precision,source,CFD_MD_ICOMM,ierr)
! send CFD nsteps
                call mpi_bcast(nsteps,1,mpi_integer,source,CFD_MD_ICOMM,ierr)

! get MD grid decomposition
                allocate(ibmin_md(npx_md), ibmax_md(npx_md), &
                        kbmin_md(npz_md),kbmax_md(npz_md))
                call mpi_bcast(ibmin_md, npx_md,mpi_integer,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(ibmax_md, npx_md,mpi_integer,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(kbmin_md, npz_md,mpi_integer,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(kbmax_md, npz_md,mpi_integer,0,CFD_MD_ICOMM,ierr)

                write(0,*)' CFD: did exchange grid data' 

        end subroutine exchange_grid_data


        subroutine create_map_cfd_md
                use mpi
                use messenger, only   : icomm_grid
                use data_export, only : nproc_cfd => nproc, jblock, ibmin, ibmax, kbmin, kbmax, &
                        imin, imax, kmin, kmax,ibmap_1,kbmap_1,nix_1,niz_1,ibmin_1,ibmax_1, &
                        kbmin_1, kbmax_1
                use mesh_export, only : x, z

                implicit none
                integer i, j, k, ip, ierr, myid_cfd, myid_xz_cfd, myid_grid, &
                        myid_world, icomm2_cfd_leader, md_rank0, source, istat
                integer ibmin_loc, ibmax_loc, kbmin_loc, kbmax_loc, three_dp
                integer a_sizes(2), a_subsizes(2), a_starts(2), ic(3)
                integer, allocatable :: overlap_aux(:,:)


! build xz subcart
                call MPI_Cart_sub (icomm_grid, (/ .true., .false., .true. /),&
                        icomm_xz, ierr)

! set a local leader for the second itercommunicator, myid_cfd == 0 sends the value across
                call mpi_comm_rank(CFD_COMM, myid_cfd, ierr)
                call mpi_comm_rank(icomm_grid, myid_grid, ierr)
                call mpi_comm_rank(MPI_COMM_WORLD, myid_world, ierr)

                icomm2_cfd_leader = -1
                call MPI_Cart_coords(icomm_grid, myid_grid, 3, ic, ierr) 

                if ( ic(1) == 0 .and. ic(2) == 0 .and. ic(3) == 0) then
                        icomm2_cfd_leader = myid_world
                endif
  
                call mpi_allreduce(MPI_IN_PLACE, icomm2_cfd_leader, 1, MPI_INTEGER,&
                        MPI_MAX, CFD_COMM, ierr )

                if ( myid_cfd == 0 ) then
!  get the rank of first processor at the lower boundary
! a more general form is needed here
                        source=MPI_ROOT
                else
                        source=MPI_PROC_NULL
                endif

                call mpi_bcast(icomm2_cfd_leader, 1, mpi_integer, source, CFD_MD_ICOMM, ierr)

! Get the mpi_comm_world of MD 0 rank

                call mpi_bcast(md_rank0, 1, mpi_integer, 0, CFD_MD_ICOMM, ierr)

                if (jblock == 1) then
! intercommunicator between cfd boundary and md

                        call mpi_intercomm_create(icomm_xz, 0, MPI_COMM_WORLD,&
                                md_rank0, 2, CFD_BDY_MD_ICOMM, ierr)

! get coordinates in the xz plane

                        do ip=1,npx_cfd*npz_cfd
                                call MPI_Cart_coords(icomm_xz, ip-1, 2, &
                                        icoord_xz(1,ip), ierr)
                        end do
                        icoord_xz = icoord_xz + 1

!    broadcast icoord_xz(1)
                        call mpi_comm_rank(icomm_xz,myid_xz_cfd,ierr)
                        if( myid_xz_cfd == 0 ) then
                                source = MPI_ROOT
                        else
                                source = MPI_PROC_NULL
                        endif

                        call mpi_bcast(icoord_xz, 2*npx_cfd*npz_cfd, &
                                mpi_integer, source,CFD_BDY_MD_ICOMM , ierr)

! get MD icoord_md

                        allocate(icoord_md(3, nproc_md))
                        call mpi_bcast(icoord_md, 3*nproc_md, mpi_integer, &
                                0,CFD_BDY_MD_ICOMM , ierr)

                        write(0,*) 'CFD side ', myid_xz_cfd, ' icoord_md ',&
                                icoord_md



! sorting out the overlap betweeen CFD and MD boxes

                        allocate(overlap_aux(8,nproc_md))
                        overlap_aux(:,:) = -1


!  three doubles to keep the velocity value
                        call mpi_type_contiguous(3, MPI_DOUBLE_PRECISION, three_dp,ierr)
                        call mpi_type_commit(three_dp, ierr)

                        novr = 0

                        do j = 1, nproc_md

                                i = icoord_md(1,j)
                                k = icoord_md(3,j) 

                                if( ((  ibmin_md(i) <= ibmin .and. ibmax_md(i) > ibmin ) .or. &
                                        ( ibmin_md(i) > ibmin .and. ibmin_md(i) < ibmax )) .and. &
                                        ((  kbmin_md(k) <= kbmin .and. kbmax_md(k) > kbmin ) .or. &
                                        ( kbmin_md(k) > kbmin .and. kbmin_md(k) < kbmax )) ) then

                                        novr = novr + 1

                                        overlap_aux(1,novr) = j 
                                        overlap_aux(2,novr) = i
                                        overlap_aux(3,novr) = k
                                        overlap_aux(4,novr) = max(ibmin,ibmin_md(i))
                                        overlap_aux(5,novr) = min(ibmax,ibmax_md(i))
                                        overlap_aux(6,novr) = max(kbmin,kbmin_md(k))
                                        overlap_aux(7,novr) = min(kbmax,kbmax_md(k))

                                        a_sizes(:) = (/ ibmax-ibmin, kbmax-kbmin /)
                                        a_subsizes(:) = (/ overlap_aux(5,novr)-overlap_aux(4,novr), overlap_aux(7,novr)-overlap_aux(6,novr) /)
                                        a_starts(:) = (/ overlap_aux(4,novr) - ibmin, overlap_aux(6,novr) - kbmin/)
                                        call mpi_type_create_subarray(2, a_sizes, a_subsizes, a_starts, MPI_ORDER_FORTRAN, three_dp, overlap_aux(8,novr), ierr)
                                        call mpi_type_commit(overlap_aux(8,novr), ierr)

                                endif
                        enddo

! copy everthing in the map_overlap
                        allocate(map_overlap(novr),stat = istat)

                        do i=1,novr
                                map_overlap(i)%rank = overlap_aux(1,i)
                                map_overlap(i)%coord(:) = overlap_aux(2:3,i)
                                map_overlap(i)%ib_range(:) = overlap_aux(4:5,i)
                                map_overlap(i)%kb_range(:) = overlap_aux(6:7,i)
                                map_overlap(i)%dp2d_type = overlap_aux(8,i)
                        enddo

                        deallocate(overlap_aux)


                        write(100+myid_xz_cfd,*) 'id ',myid_xz_cfd, ' novr ', novr
                        do i = 1, novr
                                write(100+myid_xz_cfd,'(8I5)') map_overlap(i)
                        enddo
                end if
        end subroutine create_map_cfd_md





end module coupler_cfd_setup
