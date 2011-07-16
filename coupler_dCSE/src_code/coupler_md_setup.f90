module coupler_md_setup
! 
! contains data needed by CFD and MD
! 
        use coupler_md_global_data
        implicit none
        save

        integer npx_cfd, npy_cfd, npz_cfd
        integer, allocatable :: icoord_xz_cfd(:,:) ! coordinates in XZ CFD commuincator


        integer map_shift(3)     ! shift for the MD map to position MD slab 
! with respect to CFD domain
 




contains

        subroutine exchange_grid_data
                use mpi
                implicit none

                 integer i, ierr, myid, source, ngrid(3), ny_md, nblock(3), &
                         nproc_md(3), iaux(6)
                 real(kind=kind(0.d0)) ra(2)

                 write(0,*) 'MD exchange grid data'

! get CFD processor grid
                call mpi_bcast( iaux, 3, MPI_INTEGER,&
                        0, CFD_MD_ICOMM,ierr)
                npx_cfd = iaux(1)
                npy_cfd = iaux(2)
                npz_cfd = iaux(3)
! send MD processor grid 

                call mpi_comm_rank(MD_COMM,myid,ierr)
                if ( myid == 0 ) then
                        source=MPI_ROOT
                else
                        source=MPI_PROC_NULL
                endif

                call mpi_bcast((/ npx_md, npy_md, npz_md /), 3, MPI_INTEGER,&
                        source, CFD_MD_ICOMM,ierr)

! CFD grid info {i,k,j}b{min,max} arrays
                allocate(ibmin_1(npx_cfd), ibmax_1(npx_cfd), &
                        kbmin_1(npz_cfd),kbmax_1(npz_cfd))
                call mpi_bcast(ibmin_1,npx_cfd,mpi_integer,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(ibmax_1,npx_cfd,mpi_integer,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(kbmin_1,npz_cfd,mpi_integer,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(kbmax_1,npz_cfd,mpi_integer,0,CFD_MD_ICOMM,ierr)

! CFD mesh data 
                call mpi_bcast(iaux, 6, MPI_INTEGER, 0, CFD_MD_ICOMM,ierr) 
                imino = iaux(1) ; imaxo = iaux(2)
                jmino = iaux(3) ; jmaxo = iaux(4)
                kmino = iaux(5) ; kmaxo = iaux(6)
                allocate (x(imino:imaxo),y(jmino:jmaxo),z(kmino:kmaxo))


                call mpi_bcast(x,size(x),mpi_double_precision,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(y,size(y),mpi_double_precision,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(z,size(z),mpi_double_precision,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(ra,2,mpi_double_precision,0,CFD_MD_ICOMM,ierr)

                dx = ra(1); dz = ra(2)

                write(0,*) 'MD exchage grid data: recv dx, dz ', dx, dz

! get CFD nsteps
                call mpi_bcast(nsteps,1,mpi_integer,0,CFD_MD_ICOMM,ierr)



! grid points in CFD physical domain,
! no boundaries included
! This must be revised further, a first try 

                ny_md = 3  ! prototype 
                nproc_md = (/ npx_md, npy_md, npz_md /)
                ngrid = (/ ibmax_1(npx_cfd)-ibmin_1(1)+1,      ny_md,  &
                        kbmax_1(npz_cfd)-kbmin_1(1)+1  /)   
                nblock(:) = (ngrid(:)-1)/nproc_md(:) + 1 ! grid points per MD processor
! in each direction 


                call mpi_comm_rank(mpi_comm_world,myid,ierr)

                write(0,*) 'exchange_grid_data: MD side', myid, ngrid, nblock

!  x direction
!  bounding_box_md(1,1) = (iblock-1) * nblock(1) + 1
!  bounding_box_md(2,1) = iblock * nblock(1)

!  Use names similar to CFD. Helpful or confusing ?

                ibmin_md(1) = ibmin_1(1)
                ibmax_md(1) = ibmin_1(1) + nblock(1) - 1

                do i = 2, npx_md
                        ibmin_md(i) = ibmax_md(i-1)
                        ibmax_md(i) = ibmin_md(i) + nblock(1) - 1
                enddo

                kbmin_md(1) = kbmin_1(1)
                kbmax_md(1) = kbmin_1(1) + nblock(3) - 1

                do i = 2, npz_md
                        kbmin_md(i) = kbmax_md(i-1)
                        kbmax_md(i) = kbmin_md(i) + nblock(3)-1
                enddo

! prototype stage along y dimension
                jbmin_md(1) = 1
                jbmax_md(1) = ny_md


                write(0,*) 'exchange_grid... ibmin_md', myid, ibmin_md
                write(0,*) 'exchange_grid...  ibmax_md', myid, ibmax_md
                write(0,*) 'exchange_grid... kbmin_md', myid, kbmin_md
                write(0,*) 'exchange_grid... kbmax_md', myid, kbmax_md

! fixing the dimension of MD box

                xL_md = fsig*(x(ibmax_md(npx_md))-x(ibmin_md(1)))
                yL_md = fsig*2.0d0*(y(2)-y(1))                    ! prototype
                zL_md = fsig*(z(kbmax_md(npz_md))-z(kbmin_md(1)))

                write(0,*) 'exchange_grid... xL_md, yL_md, zL_md', myid, xL_md, yL_md, zL_md

! send grid data to CFD

                call mpi_bcast(ibmin_md,npx_md,mpi_integer,source,CFD_MD_ICOMM,ierr)
                call mpi_bcast(ibmax_md,npx_md,mpi_integer,source,CFD_MD_ICOMM,ierr)
                call mpi_bcast(kbmin_md,npz_md,mpi_integer,source,CFD_MD_ICOMM,ierr)
                call mpi_bcast(kbmax_md,npz_md,mpi_integer,source,CFD_MD_ICOMM,ierr)

        end subroutine exchange_grid_data


        subroutine create_map_cfd_md
                use mpi
                use messenger, only : icomm_grid_md => icomm_grid, icoord_md => icoord
                use computational_constants_MD,    only : nproc_md => nproc
                implicit none

                integer i, j, k, ip, ierr, myid_cfd, myid_md, myid_xz_cfd, iroot,&
                        myid_world, source, istat
                integer ibmin_loc, ibmax_loc, kbmin_loc, kbmax_loc, three_dp
                integer a_sizes(2), a_subsizes(2), a_starts(2)
                integer, allocatable :: overlap_aux(:,:)
! ATTENTION Assumes at the moment that npy_md=1

! get the rank of jblock 1

                call mpi_bcast(iroot, 1, mpi_integer, 0, CFD_MD_ICOMM, ierr)

! Send to CFD the smalest MPI_COMM_WORLD rank to avoid any implicit assumption that
! might be broken e.g. the MD ranks start after CFD ranks in MPI_COMM_WORLD

                call mpi_comm_rank(MD_COMM, myid_md, ierr)
                call mpi_comm_rank(MPI_COMM_WORLD, myid_world, ierr)
                if (myid_md == 0) then
                        source = MPI_ROOT
                else
                        source = MPI_PROC_NULL
                endif

                call mpi_bcast(myid_world, 1, MPI_INTEGER, source,CFD_MD_ICOMM, ierr)

!  The MD part should be restricted to the processors that 
!  contain the CFD boundary; I'll see later to it

                call mpi_intercomm_create(icomm_grid_md, 0, MPI_COMM_WORLD, iroot, 2,&
                        CFD_BDY_MD_ICOMM, ierr)

                allocate(icoord_xz_cfd(2,npx_cfd*npz_cfd))
                call mpi_bcast(icoord_xz_cfd,2*npx_cfd*npz_cfd,mpi_integer,0,&
                        CFD_BDY_MD_ICOMM,ierr)


                call mpi_comm_rank(MD_COMM,myid_md,ierr)

                if( myid_md == 0 ) then
                        source = MPI_ROOT
                else
                        source = MPI_PROC_NULL
                endif

                write(0,*) 'MD side ', myid_md, ' icoord_xz_cfd ', icoord_xz_cfd

                call mpi_bcast(icoord_md, 3*nproc_md, mpi_integer, source, &
                        CFD_BDY_MD_ICOMM , ierr)
                allocate(overlap_aux(8,npx_cfd*npz_cfd))
                overlap_aux(:,:) = -1

!   three doubles for the speed values
                call mpi_type_contiguous(3, MPI_DOUBLE_PRECISION, three_dp,ierr)
                call mpi_type_commit(three_dp, ierr)

                ibmin_loc = ibmin_md(icoord_md(1,myid_md+1))
                ibmax_loc = ibmax_md(icoord_md(1,myid_md+1))
                kbmin_loc = kbmin_md(icoord_md(3,myid_md+1))
                kbmax_loc = kbmax_md(icoord_md(3,myid_md+1))

                novr = 0
                do j= 1, npx_cfd*npz_cfd
                        i = icoord_xz_cfd(1,j)
                        k = icoord_xz_cfd(2,j)

                        if( (( ibmin_1(i) <= ibmin_loc .and. ibmax_1(i) > ibmin_loc ) .or. &
                                ( ibmin_1(i) > ibmin_loc .and. ibmin_1(i) < ibmax_loc )) .and. &
                                ((  kbmin_1(k) <= kbmin_loc .and. kbmax_1(k) > kbmin_loc ) .or. &
                                ( kbmin_1(k) > kbmin_loc .and. kbmin_1(k) < kbmax_loc )) ) then

                                novr = novr+1

                                overlap_aux(1,novr) = j
                                overlap_aux(2,novr) = i
                                overlap_aux(3,novr) = k
                                overlap_aux(4,novr) = max(ibmin_loc,ibmin_1(i))
                                overlap_aux(5,novr) = min(ibmax_loc,ibmax_1(i))
                                overlap_aux(6,novr) = max(kbmin_loc,kbmin_1(k))
                                overlap_aux(7,novr) = min(kbmax_loc,kbmax_1(k))
! Creat subarray types for data transfers
! MD needs  subsizes correspondig to the overlaping FD domains
! at the moment this is for three dimensional vectors on a two dimensional grid

                                a_sizes(:) = (/ ibmax_loc-ibmin_loc, kbmax_loc-kbmin_loc /)
                                a_subsizes(:) = (/ overlap_aux(5,novr)-overlap_aux(4,novr), &
                                        overlap_aux(7,novr)-overlap_aux(6,novr) /)
                                a_starts(:) = (/ overlap_aux(4,novr) - ibmin_loc, &
                                        overlap_aux(6,novr) - kbmin_loc/)
                                call mpi_type_create_subarray(2, a_sizes, a_subsizes, &
                                        a_starts, MPI_ORDER_FORTRAN, three_dp, overlap_aux(8,novr), ierr)
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

                write(200+myid_md,*) 'id ',myid_cfd, ' novr ', novr
                write(200+myid_md,*) 'ibmin/max, kbmin/max ',ibmin_loc, ibmax_loc, kbmin_loc, kbmax_loc 
                write(200+myid_md,*) 'icoord_md ' , icoord_md
                write(200+myid_md,'(8I5)') map_overlap

        end subroutine create_map_cfd_md



!!$ subroutine average_vel
!!$         implicit none
!!$
!!$! First stage in MD realm
!!$!  - because time average is needed this subroutine should be called
!!$!    somewhere in  MD simulation          
!!$
!!$! First step, prototype compute v_y average on the boundary
!!$
!!$! find the cells that overlap a dx*dy*dz bin around the average point
!!$
!!$! send the average velocity to FD realm
!!$
!!$!
!!$
!!$ end subroutine average_vel

end module coupler_md_setup
