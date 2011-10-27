module coupler_cfd_setup
        use coupler_cfd_global_data
        implicit none
        save

        type md_domain_map
         integer n ! number of MD ranks that overlap with this MD domain
         integer, allocatable :: rank_list(:) ! rank list of overlpping CFD blocks
         integer, allocatable :: domains(:,:) ! indices range of overlapping with CFD grid 
        end type md_domain_map
        type(md_domain_map) md_map

! local grid sizes ( no halos, taken from bbox, helps with simpler expression of local arrays sizes)
        integer nlx, nly, nlz

! communicator for tasks that overlap  MD region 
        integer CFD_COMM_OVERLAP

contains

        subroutine exchange_grid_data(imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,&
                kmax,kmaxo,nsteps,x,y,z,dx,dz,&
                npx,npy,npz,icoord,dt)
                use mpi
                use coupler_cfd_global_data, only : imino_ => imino, imin_ => imin, imax_ => imax, jmino_ => jmin, &
                        jmin_ => jmin, jmax_ => jmax, jmaxo_ => jmaxo, kmino_ => kmino, kmin_ => kmin, kmax_ => kmax, &
                        kmaxo_ => kmaxo, nsteps_ => nsteps, x_ => x, y_ => y, z_ => z, dx_ => dx, dz_ => dz, &
                        npx_ => npx, npy_ => npy, npz_ => npz, icoord_ => icoord, dt_ => dt
                
                implicit none


                integer, intent(in) :: imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,kmax,kmaxo,nsteps,&
                        npx,npy,npz,icoord(:,:)
                real(kind(0.d0)), intent(in)    :: x(:),y(:),z(:),dx,dz,dt

                 integer i, ierr, myid, source, ngrid(3), ny_md, nblock(3), &
                         iaux(6)

                 
                 imino_ = imino; imin_ = imin; imax_ = imax; jmino_ = jmin; &
                         jmin_ = jmin; jmax_ = jmax; jmaxo_ = jmaxo; kmino_ = kmino; kmin_ = kmin; kmax_ = kmax; &
                         kmaxo_ = kmaxo; nsteps_ = nsteps;  dx_ = dx; dz_ = dz; &
                         npx_ = npx; npy_ = npy; npz_ = npz 
                 
                 allocate(x_(size(x)),stat=ierr); x_ = x
                 allocate(y_(size(y)),stat=ierr); y_ = y
                 allocate(z_(size(z)),stat=ierr); z_ = z
                 allocate(icoord_(3,npx*npy*npz),stat=ierr); icoord_=icoord
                 

                 x_ = x; y_ = y; z_ = z;
                 icoord_ = icoord

!                 write(0,*) 'CFD exchange grid data'
! send CFD processor grid
                call mpi_comm_rank(CFD_COMM,myid,ierr)
                if ( myid == 0 ) then
                        source=MPI_ROOT
                else
                        source=MPI_PROC_NULL
                endif

                call mpi_bcast((/ npx, npy, npz, jmax_overlap /), 4, MPI_INTEGER,&
                        source, CFD_MD_ICOMM,ierr)

! receive MD processor grid 
                call mpi_bcast(iaux, 3, MPI_INTEGER,&
                        0, CFD_MD_ICOMM,ierr)

                npx_md = iaux(1)
                npy_md = iaux(2)
                npz_md = iaux(3)
                nproc_md = npx_md * npy_md * npz_md

! send CFD mesh data

                call mpi_bcast((/ imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,kmax,kmaxo /), 12,&
                        MPI_INTEGER, source, CFD_MD_ICOMM,ierr)
                call mpi_bcast(x,size(x),mpi_double_precision,source,CFD_MD_ICOMM,ierr)
                call mpi_bcast(y,size(y),mpi_double_precision,source,CFD_MD_ICOMM,ierr)
                call mpi_bcast(z,size(z),mpi_double_precision,source,CFD_MD_ICOMM,ierr)
                call mpi_bcast((/ dx, dz /),2,mpi_double_precision,source,CFD_MD_ICOMM,ierr)
! send CFD nsteps and dt
                call mpi_bcast(nsteps,1,mpi_integer,source,CFD_MD_ICOMM,ierr)
                call mpi_bcast(dt,1,mpi_double_precision,source,CFD_MD_ICOMM,ierr)

!                 write(0,*)' CFD: did exchange grid data'

        end subroutine exchange_grid_data
       
       
        subroutine create_map_cfd_md
               use mpi 
!               use coupler_cfd_global_data, only : CFD_COMM, CFD_MD_ICOMM, bbox_cfd
               implicit none
               
               integer i, myid, color, noverlaps, ir, ierr
               integer, allocatable :: md_grid_boxes(:,:), overlap_mask(:), ireq(:), overlap_box(:,:)

               call mpi_comm_rank(CFD_COMM,myid,ierr)

                 call make_bbox

!                 write(0,*)'CFD, bbox_cdf',  bbox_cfd%xbb, bbox_cfd%ybb, bbox_cfd%zbb

!  get the block boundaries cover by each MD domain
                
                allocate(md_grid_boxes(6,0:nproc_md - 1), overlap_mask(0:nproc_md - 1), &
                         overlap_box(6,0:nproc_md-1), ireq(0:nproc_md - 1))
                
                call mpi_allgather(MPI_BOTTOM, 0, MPI_INTEGER, md_grid_boxes, 6, MPI_INTEGER,CFD_MD_ICOMM, ierr)

!                write(0,*) ' CFD grid boxes ', myid, md_grid_boxes

                call find_overlaps

! send the overlap mask across to MD
                call mpi_allgather(overlap_mask,nproc_md,MPI_INTEGER,MPI_BOTTOM,0, MPI_INTEGER,CFD_MD_ICOMM,ierr)

                noverlaps = 0
                do i = 0, nproc_md - 1
                 if ( overlap_mask(i) == 1) then 
                  noverlaps = noverlaps + 1
                 endif
                enddo                

                 md_map%n = noverlaps
                allocate ( md_map%rank_list(noverlaps), md_map%domains(6,noverlaps))

!  overlaping communicator
                if ( md_map%n > 0) then
                 color = 1
                else 
                 color = 0
                endif

                call mpi_comm_split(CFD_COMM, color, myid, CFD_COMM_OVERLAP, ierr)

                if (color == 0) then
                 CFD_COMM_OVERLAP = MPI_COMM_NULL
                endif

! send the range of the overlaping domains

                ir = 0
                do i=0, nproc_md-1
                 if (overlap_mask(i) == 1) then
                  call mpi_isend(overlap_box(1,i),6,mpi_integer,i,2,CFD_MD_ICOMM,ireq(i),ierr)
                  ir = ir + 1
                  md_map%rank_list(ir) = i
                  md_map%domains(:,ir) = overlap_box(1:6,i)
                 else
                  ireq(i) = MPI_REQUEST_NULL
                 endif
                enddo
                call mpi_waitall(nproc_md,ireq,MPI_STATUSES_IGNORE,ierr) 

        contains
      
                subroutine find_overlaps
                       implicit none
                       
                       integer i, ibmin,ibmax,jbmin,jbmax,kbmin,kbmax, &
                              ibs, ibe, jbs, jbe, kbs, kbe
                       
                       ibmin = bbox_cfd%xbb(1,icoord(1,myid + 1))
                       ibmax = bbox_cfd%xbb(2,icoord(1,myid + 1))
                       jbmin = bbox_cfd%ybb(1,icoord(2,myid + 1))
                       jbmax = min(jmax_overlap,bbox_cfd%ybb(2,icoord(2,myid + 1))) 
                       kbmin = bbox_cfd%zbb(1,icoord(3,myid + 1))
                       kbmax = bbox_cfd%zbb(2,icoord(3,myid + 1))

!                       write(0, *) 'CFD: find overlap ibmin ...', myid, ibmin, ibmax, jbmin, jbmax, kbmin, kbmax

                       do i=0,nproc_md - 1
                       
                        ibs = md_grid_boxes(1,i)
                        ibe = md_grid_boxes(2,i)
                        jbs = md_grid_boxes(3,i)
                        jbe = md_grid_boxes(4,i)
                        kbs = md_grid_boxes(5,i)
                        kbe = md_grid_boxes(6,i)

                        if( ((   ibs <  ibmin .and. ibe > ibmin )        .or.  &
                              (  ibs >= ibmin .and. ibs < ibmax ))       .and. &
                              (( jbs <  jbmin .and. jbe > jbmin )        .or.  &
                               ( jbs >= jbmin .and. jbs < jbmax ))       .and. &   
                              (( kbs <  kbmin .and. kbe > kbmin )        .or.  &
                               ( kbs >= kbmin .and. kbs < kbmax )) ) then

                         overlap_mask(i) = 1

                         overlap_box(1,i) = max(ibmin,ibs)
                         overlap_box(2,i) = min(ibmax,ibe)
                         overlap_box(3,i) = max(jbmin,jbs)
                         overlap_box(4,i) = min(jmax_overlap,jbe)
                         overlap_box(5,i) = max(kbmin,kbs)
                         overlap_box(6,i) = min(kbmax,kbe)

                         else
                          overlap_mask(i) = 0
                          overlap_box(:,i) = -666

                        endif

                       enddo

!                       write(0,*)' CFD: find overlap ', myid, overlap_mask, overlap_box

                end subroutine find_overlaps
                

                subroutine make_bbox
!                        use coupler_cfd_global_data, only : npx, npy, npz, imin, jmin, kmin,icoord
                        implicit none

                        integer, parameter :: is = 1, ie = 2

                        integer ixyz, i, nixyz(3), minxyz(3), npxyz(3)
                        integer, pointer :: bb_ptr(:,:) => null()

! number of grid per MPI task, remainder must be added !!!
                        nixyz  = (/ (imax - imin) / npx + 1, (jmax-jmin) / npy + 1, (kmax - kmin) / npz + 1/)
                        minxyz = (/ imin,  jmin,  kmin     /)
                        npxyz  = (/ npx, npy, npz  /)

                        allocate(bbox_cfd%xbb(2,npx),bbox_cfd%ybb(2,npy), bbox_cfd%zbb(2,npz))

                        do ixyz = 1,3

                                select case(ixyz)
                                  case(1)
                                        bb_ptr => bbox_cfd%xbb
                                  case(2)
                                        bb_ptr => bbox_cfd%ybb
                                  case(3) 
                                        bb_ptr => bbox_cfd%zbb
                                end select

                                bb_ptr(is,1) = minxyz(ixyz)
                                bb_ptr(ie,1) = bb_ptr(is,1) + nixyz(ixyz) - 1

                                do i = 2, npxyz(ixyz)
                                        bb_ptr(is, i) = bb_ptr(ie, i-1)
                                        bb_ptr(ie, i) = bb_ptr(is, i) + nixyz(ixyz) - 1
                                enddo

                        enddo
! sizes of local grids

                        nlx = bbox_cfd%xbb(2,icoord(1,myid+1)) - bbox_cfd%xbb(1,icoord(1,myid + 1)) + 1
                        nly = min(bbox_cfd%ybb(2,icoord(2,myid+1)),jmax_overlap) - bbox_cfd%ybb(1,icoord(2,myid + 1)) + 1
                        nlz = bbox_cfd%zbb(2,icoord(3,myid+1)) - bbox_cfd%zbb(1,icoord(3,myid + 1)) + 1
                        

!!$! send box information accros to MD
!!$
!!$                        call mpi_bcast(bbox_cfd%xbb(1,1),2*npx_cfd,mpi_integer,source,&
!!$                                CFD_MD_ICOMM,ierr)
!!$                        call mpi_bcast(bbox_cfd%zbb(1,1),2*npz_cfd,mpi_integer,source,&
!!$                                CFD_MD_ICOMM,ierr)
!!$
!!$! receive MD grid decomposition
!!$                        allocate(bbox_md%xbb(2,npx_md), bbox_md%zbb(2,npz_md))
!!$                        call mpi_bcast(bbox_md%xbb(1,1), 2*npx_md,mpi_integer,0,CFD_MD_ICOMM,ierr)
!!$                        call mpi_bcast(bbox_md%zbb(1,1), 2*npz_md,mpi_integer,0,CFD_MD_ICOMM,ierr)


                end subroutine make_bbox



!!$                use data_export, only : nproc_cfd => nproc, jblock, &
!!$                        imin, imax, kmin, kmax,ibmap_1,kbmap_1
!!$                use mesh_export, only : x, z
!!$                
!!$                implicit none
!!$                integer i, j, k, ip, ierr, myid_cfd, myid_xz_cfd, myid_grid, &
!!$                        myid_world, icomm2_cfd_leader, md_rank0, source, istat, &
!!$                        ibs, ibe, kbs, kbe
!!$                integer ibmin, ibmax, kbmin, kbmax, four_dp_type
!!$                integer a_sizes(3), a_subsizes(3), a_starts(3), ic(3)
!!$                integer, allocatable :: overlap_aux(:,:)
!!$
!!$
!!$! build xz subcart
!!$                call MPI_Cart_sub (icomm_grid, (/ .true., .false., .true. /),&
!!$                        icomm_xz, ierr)
!!$
!!$! set a local leader for the second itercommunicator, myid_cfd == 0 sends the value across
!!$                call mpi_comm_rank(CFD_COMM, myid_cfd, ierr)
!!$                call mpi_comm_rank(icomm_grid, myid_grid, ierr)
!!$                call mpi_comm_rank(MPI_COMM_WORLD, myid_world, ierr)
!!$
!!$                icomm2_cfd_leader = -1
!!$                call MPI_Cart_coords(icomm_grid, myid_grid, 3, ic, ierr) 
!!$
!!$                if ( ic(1) == 0 .and. ic(2) == 0 .and. ic(3) == 0) then
!!$                        icomm2_cfd_leader = myid_world
!!$                endif
!!$  
!!$                call mpi_allreduce(MPI_IN_PLACE, icomm2_cfd_leader, 1, MPI_INTEGER,&
!!$                        MPI_MAX, CFD_COMM, ierr )
!!$
!!$                if ( myid_cfd == 0 ) then
!!$!  get the rank of first processor at the lower boundary
!!$! a more general form is needed here
!!$                        source=MPI_ROOT
!!$                else
!!$                        source=MPI_PROC_NULL
!!$                endif
!!$
!!$                call mpi_bcast(icomm2_cfd_leader, 1, mpi_integer, source, CFD_MD_ICOMM, ierr)
!!$
!!$! Get the mpi_comm_world of MD 0 rank
!!$
!!$                call mpi_bcast(md_rank0, 1, mpi_integer, 0, CFD_MD_ICOMM, ierr)
!!$
!!$                if (jblock == 1) then
!!$! intercommunicator between cfd boundary and md
!!$
!!$                        call mpi_intercomm_create(icomm_xz, 0, MPI_COMM_WORLD,&
!!$                                md_rank0, 2, CFD_BDY_MD_ICOMM, ierr)
!!$
!!$! get coordinates in the xz plane
!!$
!!$                        do ip=1,npx_cfd*npz_cfd
!!$                                call MPI_Cart_coords(icomm_xz, ip-1, 2, &
!!$                                        icoord_xz(1,ip), ierr)
!!$                        end do
!!$                        icoord_xz = icoord_xz + 1
!!$
!!$!    broadcast icoord_xz(1)
!!$                        call mpi_comm_rank(icomm_xz,myid_xz_cfd,ierr)
!!$                        if( myid_xz_cfd == 0 ) then
!!$                                source = MPI_ROOT
!!$                        else
!!$                                source = MPI_PROC_NULL
!!$                        endif
!!$
!!$                        call mpi_bcast(icoord_xz, 2*npx_cfd*npz_cfd, &
!!$                                mpi_integer, source,CFD_BDY_MD_ICOMM , ierr)
!!$
!!$! get MD icoord_md
!!$
!!$                        allocate(icoord_md(3, nproc_md))
!!$                        call mpi_bcast(icoord_md, 3*nproc_md, mpi_integer, &
!!$                                0,CFD_BDY_MD_ICOMM , ierr)
!!$
!!$                        write(0,*) 'CFD side ', myid_xz_cfd, ' icoord_md ',&
!!$                                icoord_md
!!$
!!$
!!$
!!$! sorting out the overlap betweeen CFD and MD boxes
!!$
!!$                        allocate(overlap_aux(9,nproc_md))
!!$                        overlap_aux(:,:) = -1
!!$
!!$
!!$!  three doubles to keep the velocity value
!!$!                        call mpi_type_contiguous(4, MPI_DOUBLE_PRECISION, four_dp_type,ierr)
!!$!                        call mpi_type_commit(four_dp_type, ierr)
!!$                        four_dp_type=MPI_DOUBLE_PRECISION
!!$
!!$                        ibmin = bbox_cfd%xbb(1,icoord_xz(1,myid_xz_cfd + 1))
!!$                        ibmax = bbox_cfd%xbb(2,icoord_xz(1,myid_xz_cfd + 1))
!!$                        kbmin = bbox_cfd%zbb(1,icoord_xz(2,myid_xz_cfd + 1))
!!$                        kbmax = bbox_cfd%zbb(2,icoord_xz(2,myid_xz_cfd + 1))
!!$
!!$                        novr = 0
!!$
!!$                        do j = 1, nproc_md
!!$
!!$                                i = icoord_md(1,j)
!!$                                k = icoord_md(3,j) 
!!$
!!$                                ibs = bbox_md%xbb(1,i)
!!$                                ibe = bbox_md%xbb(2,i)
!!$                                kbs = bbox_md%zbb(1,k)
!!$                                kbe = bbox_md%zbb(2,k)
!!$
!!$                                if( ((   ibs <  ibmin .and. ibe > ibmin ) .or. &
!!$                                      (  ibs >= ibmin .and. ibs < ibmax )) .and. &
!!$                                      (( kbs <  kbmin .and. kbe > kbmin ) .or. &
!!$                                       ( kbs >= kbmin .and. kbs < kbmax )) ) then
!!$
!!$                                        novr = novr + 1
!!$
!!$                                        overlap_aux(1,novr) = j 
!!$                                        overlap_aux(2,novr) = i
!!$                                        overlap_aux(3,novr) = k
!!$                                        overlap_aux(4,novr) = max(ibmin,ibs)
!!$                                        overlap_aux(5,novr) = min(ibmax,ibe)
!!$                                        overlap_aux(6,novr) = max(kbmin,kbs)
!!$                                        overlap_aux(7,novr) = min(kbmax,kbe)
!!$
!!$                                        a_sizes(:) = (/ kbmax-kbmin, ibmax-ibmin, 1  /)
!!$                                        a_subsizes(:) = (/ overlap_aux(7,novr)-overlap_aux(6,novr), overlap_aux(5,novr)-overlap_aux(4,novr), 1  /)
!!$                                        a_starts(:) = (/ overlap_aux(6,novr) - kbmin, overlap_aux(4,novr) - ibmin, 0 /)
!!$
!!$                                        write(0,*) 'CFD side setup:', myid_xz_cfd, a_sizes, a_subsizes, a_starts
!!$
!!$                                        call mpi_type_create_subarray(3, a_sizes, a_subsizes, a_starts, MPI_ORDER_FORTRAN, four_dp_type, overlap_aux(8,novr), ierr)
!!$                                        call mpi_type_commit(overlap_aux(8,novr), ierr)
!!$
!!$                                        a_sizes(:) = (/ kbmax-kbmin+1, ibmax-ibmin+1, 1  /)
!!$                                        a_subsizes(:) = (/ overlap_aux(7,novr)-overlap_aux(6,novr)+1, overlap_aux(5,novr)-overlap_aux(4,novr)+1, 1  /)
!!$                                        a_starts(:) = (/ overlap_aux(6,novr) - kbmin, overlap_aux(4,novr) - ibmin, 0 /)
!!$
!!$                                       
!!$
!!$                                        call mpi_type_create_subarray(3, a_sizes, a_subsizes, a_starts, MPI_ORDER_FORTRAN, four_dp_type, overlap_aux(9,novr), ierr)
!!$                                        call mpi_type_commit(overlap_aux(9,novr), ierr)
!!$                                        
!!$                                         
!!$
!!$                                endif
!!$                        enddo
!!$
!!$                          write(0,*) 'CFD side ', myid_xz_cfd, ' novr ', novr
!!$
!!$! copy everthing in the map_overlap
!!$                        allocate(map_overlap(novr),stat = istat)
!!$
!!$                        do i=1,novr
!!$                                map_overlap(i)%rank = overlap_aux(1,i)
!!$                                map_overlap(i)%coord(:) = overlap_aux(2:3,i)
!!$                                map_overlap(i)%ib_range(:) = overlap_aux(4:5,i)
!!$                                map_overlap(i)%kb_range(:) = overlap_aux(6:7,i)
!!$                                map_overlap(i)%dp2d_type = overlap_aux(8,i)
!!$                                map_overlap(i)%plane_type = overlap_aux(9,i)
!!$                        enddo
!!$
!!$                        deallocate(overlap_aux)
!!$
!!$
!!$! save the number of cfd cells in each local domain for future use
!!$                        
!!$                        nib = ibmax - ibmin 
!!$                        nkb = kbmax - kbmin 
!!$
!!$                        write(100+myid_xz_cfd,*) 'id ',myid_xz_cfd, ' novr ', novr
!!$                        do i = 1, novr
!!$                                write(100+myid_xz_cfd,'(8I5)') map_overlap(i)
!!$                        enddo
!!$                end if
        end subroutine create_map_cfd_md





end module coupler_cfd_setup
