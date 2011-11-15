!=============================================================================
!                                   
! Internal data and subroutines used by the coupler when working in CFD realm
! It must not be used by subroutimes working in MD realm
!
!-----------------------------------------------------------------------------
module coupler_internal_cfd
        implicit none
        save


! CFD data 
        integer imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,&
         kmino,kmin,kmax,kmaxo,npx,npy,npz, nsteps

        integer tplot               !Output every number of steps
        integer :: jmax_overlap = 5 ! maximum j index ( in y direction) which MD 
! has to cover, on top of that it has to cover 
! y(0):y(1) domain and a bit of room

        integer, allocatable :: icoord(:,:)

        real(kind(0.d0)), allocatable :: x(:), y(:), z(:)
        real(kind(0.d0)) dx, dz, dt

! Data recieved form MD

        integer npx_md, npy_md, npz_md, nproc_md

! Internal data 
! bounding box type that holds on each processor the boundaries of the
! CFD grid sectors allocate to every processor
        type bbox
         integer, allocatable :: xbb(:,:), ybb(:,:), zbb(:,:)
        end type bbox

        type(bbox),target :: bbox_cfd, bbox_md
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

        subroutine create_map_cfd
                use mpi 
               use coupler_parameters, only : COUPLER_COMM
               use coupler_internal_common, only : CFD_MD_ICOMM
                implicit none

                integer i, myid, color, noverlaps, ir, ierr
                integer, allocatable :: md_grid_boxes(:,:), overlap_mask(:), ireq(:), overlap_box(:,:)

                call mpi_comm_rank(COUPLER_COMM,myid,ierr)

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

                call mpi_comm_split(COUPLER_COMM, color, myid, CFD_COMM_OVERLAP, ierr)

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

        end subroutine create_map_cfd

        subroutine recv_vel_MD(a,p1s,p1e,p2s,p2e,p3s,p3e,pbc)
                use mpi
                use coupler_parameters, only : COUPLER_COMM
                use coupler_internal_common, only : CFD_MD_ICOMM
                implicit none
! the index ranges in z,x,y, periodic BC
                integer, intent(in) :: p1s,p1e,p2s,p2e,p3s,p3e,pbc
                real(kind=kind(0.d0)), intent(out) :: a(:,:,:)

                integer i,i1,j,k, is(md_map%n), ie(md_map%n), ks(md_map%n), ke(md_map%n), &
                 ny, ierr, source, myid, itag, type, req(md_map%n), &
                 start_address(md_map%n+1), min_i, min_j, min_k, np
                real(kind(0.d0)) x
                real(kind(0.d0)), allocatable :: vbuf(:), v1(:,:,:,:), vtot(:,:,:,:)
                integer, save :: ncalls = 0

! This local CFD domain is outside MD overlap zone 
                if ( md_map%n == 0 ) return 

                call mpi_comm_rank(COUPLER_COMM, myid, ierr)

                ncalls = ncalls + 1

                ny = p3e - p3s + 1 ! number of y planes

                min_i = minval(md_map%domains(1,:))
                min_j = minval(md_map%domains(3,:))
                min_k = minval(md_map%domains(5,:))

                np = 0
                do i = 1, md_map%n
                        np = np + 2 * (md_map%domains(2,i) - md_map%domains(1,i) + 0) &
                         * (md_map%domains(6,i) - md_map%domains(5,i) + 0) &
                         * ny

                        is(i) = md_map%domains(1,i) - min_i + 1
                        ie(i) = md_map%domains(2,i) - min_i + 1
                        ks(i) = md_map%domains(5,i) - min_k + 1
                        ke(i) = md_map%domains(6,i) - min_k + 1

                enddo

                allocate(vbuf(np),stat=ierr)
                vbuf=0.d0

!                write(0,*)'CFD recv_MDvel, vbuf size ', myid, np

                start_address(1) = 1
                do i = 1, md_map%n

                        source = md_map%rank_list(i)
                        source = md_map%rank_list(i)
                        np = 2 * (ke(i) - ks(i) + 0) * (ie(i) - is(i) + 0) * ny

                        start_address(i+1) = start_address(i) + np

! Attention ncall could go over max tag value for long runs!!
                        itag = mod(ncalls, MPI_TAG_UB)
                        call mpi_irecv(vbuf(start_address(i)),np, MPI_DOUBLE_PRECISION, source, itag, CFD_MD_ICOMM, &
                         req(i),ierr)
!                                write(0,*) 'CFD recv_MDvel  ', myid, i, itag,source,np, is,ie,ks,ke,ierr        
                enddo

                call mpi_waitall(md_map%n, req, MPI_STATUSES_IGNORE, ierr)

                allocate(vtot(2,nlz-1,nlx-1,ny),stat=ierr)
                vtot=0.d0

                do i = 1, md_map%n

                        if ( allocated(v1)) deallocate (v1)
                        allocate(v1(2,ks(i):ke(i)-1, is(i):ie(i)-1,ny))

                        v1(:,:,:,:) = reshape(vbuf(start_address(i):start_address(i+1)-1), &
                         (/ 2, ke(i)-ks(i)+0, ie(i)-is(i)+0, ny /))

                        vtot(:,ks(i):ke(i)-1, is(i):ie(i)-1,:) = vtot(:,ks(i):ke(i)-1, is(i):ie(i)-1,:) +  v1(:,:,:,:)

                enddo

! Periodic boundary condition?          

                call set_pbc(pbc)


                where (vtot(2,1:p1e-p1s+1,1:p2e-p2s+1,1:p3e-p3s+1) >= 0.d0) 
                        a(p1s:p1e,p2s:p2e,p3s:p3e) = vtot(1,1:p1e-p1s+1,1:p2e-p2s+1,1:p3e-p3s+1) &
                         /vtot(2,1:p1e-p1s+1,1:p2e-p2s+1,1:p3e-p3s+1)
                elsewhere
                        a(p1s:p1e,p2s:p2e,p3s:p3e) = 0.d0
                end where

        contains

                subroutine set_pbc(pbc)
                        implicit none
                        integer, intent(in) :: pbc

                        real(kind(0.d0)), dimension(2,p2s:p2e,size(vtot,dim=4)) :: x1, x2

                        select case(pbc)
                        case(1)
! In the case of PBC in z add the  sums of the boundary (half) cell
                                x1 =  vtot(:,1, 1:p2e-p2s+1,:)
                                x2 =  vtot(:,p1e-p1s+1, 1:p2e-p2s+1,:)
                                vtot(:,1, 1:p2e-p2s+1,:) =   vtot(:,1, 1:p2e-p2s+1,:) + x2
                                vtot(:,p1e-p1s+1, 1:p2e-p2s+1,:) =   vtot(:,p1e-p1s+1, 1:p2e-p2s+1,:) + x1

                        end select
                end subroutine set_pbc

        end subroutine recv_vel_MD



end module coupler_internal_cfd
