!=============================================================================
!				   
! Internal data and subroutines used by the coupler when working in CFD realm
! It must not be used by subroutimes working in MD realm
! Subroutines include:
!
! create_map_cfd 			Establish for all CFD processors the mapping (if any) 
! 							to coupled MD processors
! find_overlaps				Establish location and size of the overlap region 
! 							between the continuum and molecular simulations
! make_bbox					Make bbox which contains all domain information 
!							required in coupling process such as domain extents
! recv_vel_MD(vel,p1s,p1e, 	Get velocity fields from MD for 
!	   p2s,p2e,p3s,p3e,pbc) the boundary condiitons needed in CFD
! set_pbc(pbc)
!
!  Lucian Anton, November 2011
!
!-----------------------------------------------------------------------------

module coupler_internal_cfd
    implicit none
    save

    ! CFD data 

    ! CFD grid indices, number of processor in grid, number of steps
    integer imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,&
        kmino,kmin,kmax,kmaxo,npx,npy,npz, nsteps

    integer tplot	       !Output every number of steps
    integer :: jmax_overlap = 5 ! maximum j index ( in y direction) which MD 

    ! coordinates of CFD topologies
    ! ATTENYION the values are shifted with +1, FORTRAN style
    ! if this array is passed a MPI function, remove the shift !
    integer, allocatable :: icoord(:,:)

    real(kind(0.d0)), allocatable :: x(:), y(:), z(:)
    real(kind(0.d0)) dx, dz, dt

    ! size of initialisation MD cell  used to resize CFD domain, derived from density
    ! and cell tipe
    real(kind(0.d0)) MD_initial_cellsize

    ! Data recieved form MD
    integer npx_md, npy_md, npz_md, nproc_md

    ! Internal data 
    ! bounding box type that holds on each processor the boundaries of the
    ! CFD grid sectors allocate to every processor
    type bbox
        integer, allocatable :: xbb(:,:), ybb(:,:), zbb(:,:)
    end type bbox

    type(bbox),target :: bbox_cfd, bbox_md

    ! local grid sizes ( no halos, taken from bbox, helps with simpler expression of local arrays sizes)
    integer nlx, nly, nlz
    ! communicator for tasks that overlap  MD region 
    integer CFD_COMM_OVERLAP


contains

    !=============================================================================
    ! Establish for all CFD processors the mapping (if any) 
    ! to coupled MD processors
    !-----------------------------------------------------------------------------

    subroutine create_map_cfd
        use mpi 
        use coupler_internal_common, only : COUPLER_REALM_COMM, COUPLER_GRID_COMM, COUPLER_ICOMM, cfd_is_2d, map
        implicit none

        integer i, myid, id_coord, color, noverlaps, ir, iaux(4), ierr
        integer, allocatable :: md_grid_boxes(:,:), overlap_mask(:), ireq(:), overlap_box(:,:)
        real(kind(0.d0)) raux(2)

        call mpi_comm_rank(COUPLER_REALM_COMM,myid,ierr)

        ! coupler_grid_comm is the correct communicator to pick the processor coordinates in cartesian topology
        call mpi_comm_rank(coupler_grid_comm,id_coord,ierr)
        id_coord = id_coord+1 ! get in sync with icoord convention

        ! for 2D CFD problem one must broadcast  back the zL_md,z,dz
        if (cfd_is_2d) then
            call mpi_bcast( raux, 2, MPI_DOUBLE_PRECISION, 0, COUPLER_ICOMM,ierr)
            if (allocated(z)) then
                deallocate(z)
            endif
            allocate(z(2))
            z(:) = raux(:)
            dz   = z(2) - z(1)
            call mpi_bcast(iaux, 4, MPI_INTEGER, 0, COUPLER_ICOMM,ierr)
            kmino=iaux(1); kmin=iaux(2); kmax=iaux(3); kmaxo=iaux(4)
        endif

        call make_bbox

        ! Get the block boundaries covered by each MD domain
        allocate(md_grid_boxes(6,0:nproc_md - 1), overlap_mask(0:nproc_md - 1), &
            overlap_box(6,0:nproc_md-1), ireq(0:nproc_md - 1))

        call mpi_allgather(MPI_BOTTOM, 0, MPI_INTEGER, md_grid_boxes, 6, MPI_INTEGER,COUPLER_ICOMM, ierr)
        !write(0,*) ' CFD grid boxes ', myid, md_grid_boxes
        call find_overlaps
        ! send the overlap mask across to MD
        call mpi_allgather(overlap_mask,nproc_md,MPI_INTEGER,MPI_BOTTOM,0, MPI_INTEGER,COUPLER_ICOMM,ierr)

        noverlaps = 0
        do i = 0, nproc_md - 1
            if ( overlap_mask(i) == 1) then 
                noverlaps = noverlaps + 1
            endif
        enddo

        map%n = noverlaps
        allocate ( map%rank_list(noverlaps), map%domains(6,noverlaps))

        !  overlaping communicator
        if ( map%n > 0) then
            color = 1
        else 
            color = 0
        endif

        call mpi_comm_split(COUPLER_REALM_COMM, color, myid, CFD_COMM_OVERLAP, ierr)

        if (color == 0) then
            CFD_COMM_OVERLAP = MPI_COMM_NULL
        endif

        ! send the range of the overlaping domains
        ir = 0
        do i=0, nproc_md-1
            if (overlap_mask(i) == 1) then
                call mpi_isend(overlap_box(1,i),6,mpi_integer,i,2,COUPLER_ICOMM,ireq(i),ierr)
                ir = ir + 1
                map%rank_list(ir) = i
                map%domains(:,ir) = overlap_box(1:6,i)
            else
                ireq(i) = MPI_REQUEST_NULL
            endif
        enddo
        call mpi_waitall(nproc_md,ireq,MPI_STATUSES_IGNORE,ierr) 

    contains

        !-----------------------------------------------------------------------------
        ! Establish location and size of the overlap region between the continuum
        ! and molecular simulations
        !-----------------------------------------------------------------------------

        subroutine find_overlaps
            implicit none

            integer i, ibmin,ibmax,jbmin,jbmax,kbmin,kbmax, &
                ibs, ibe, jbs, jbe, kbs, kbe

            ibmin = bbox_cfd%xbb(1,icoord(1,id_coord))
            ibmax = bbox_cfd%xbb(2,icoord(1,id_coord))
            jbmin = bbox_cfd%ybb(1,icoord(2,id_coord))
            jbmax = min(jmax_overlap,bbox_cfd%ybb(2,icoord(2,id_coord))) 
            kbmin = bbox_cfd%zbb(1,icoord(3,id_coord))
            kbmax = bbox_cfd%zbb(2,icoord(3,id_coord))

            ! write(0, *) 'CFD: find overlap ibmin ...', myid, ibmin, ibmax, jbmin, jbmax, kbmin, kbmax

            do i=0,nproc_md - 1

                ibs = md_grid_boxes(1,i)
                ibe = md_grid_boxes(2,i)
                jbs = md_grid_boxes(3,i)
                jbe = md_grid_boxes(4,i)
                kbs = md_grid_boxes(5,i)
                kbe = md_grid_boxes(6,i)

                if	 ((( ibs <  ibmin .and. ibe > ibmin )	.or.  &
                    ( ibs >= ibmin .and. ibs < ibmax ))  .and. &
                    (( jbs <  jbmin .and. jbe > jbmin )	.or.  &
                    ( jbs >= jbmin .and. jbs < jbmax))   .and. &   
                    (( kbs <  kbmin .and. kbe > kbmin )	.or.  &
                    ( kbs >= kbmin .and. kbs < kbmax)))  then

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

            ! write(0,*)' CFD: find overlap ', myid, overlap_mask, overlap_box

        end subroutine find_overlaps

        !-----------------------------------------------------------------------------
        ! Make bbox which contains all domain information required in coupling process
        ! such as domain extents
        !-----------------------------------------------------------------------------

        subroutine make_bbox
            !use coupler_cfd_global_data, only : npx, npy, npz, imin, jmin, kmin,icoord
            implicit none

            integer, parameter :: is = 1, ie = 2

            integer ixyz, i, nixyz(3), minxyz(3), npxyz(3)
            integer, pointer :: bb_ptr(:,:) => null()

            ! number of grid per MPI task, remainder must be added !!!
            nixyz  = (/ (imax - imin) / npx + 1, (jmax-jmin) / npy + 1, (kmax - kmin) / npz + 1/)
            minxyz = (/ imin,  jmin,  kmin /)
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

            ! set sizes of local grids
            nlx = bbox_cfd%xbb(2,icoord(1,id_coord)) - bbox_cfd%xbb(1,icoord(1,id_coord)) + 1
            nly = min(bbox_cfd%ybb(2,icoord(2,id_coord)),jmax_overlap) - bbox_cfd%ybb(1,icoord(2,id_coord)) + 1
            nlz = bbox_cfd%zbb(2,icoord(3,id_coord)) - bbox_cfd%zbb(1,icoord(3,id_coord)) + 1

        end subroutine make_bbox


    end subroutine create_map_cfd

    !=============================================================================
    ! Get velocity sums from MD and write them to the specified addresses in vel arrays
    !-----------------------------------------------------------------------------
!!$    subroutine recv_vel_MD(vel,p1ls,p1le,p1gs,p2ls,p2le,p2gs,p3s,p3e,pbc)
!!$        use mpi
!!$        use coupler_internal_common, only : COUPLER_REALM_COMM, COUPLER_ICOMM
!!$        implicit none
!!$        ! the index ranges in z,x,y, periodic BC
!!$        integer, intent(in) :: p1ls,p1le,p1gs,p2ls,p2le,p2gs,p3s,p3e,pbc
!!$        real(kind=kind(0.d0)), intent(out) :: vel(0:,0:,0:)
!!$
!!$        integer i,i1,j,k, is(md_map%n), ie(md_map%n), ks(md_map%n), ke(md_map%n), &
!!$            ny, ierr, source, myid, itag, type, req(md_map%n), &
!!$            start_address(md_map%n+1), min_i, min_j, min_k, np, gs1,gs2
!!$        real(kind(0.d0)) x
!!$        real(kind(0.d0)), allocatable :: vbuf(:), v1(:,:,:,:), vtot(:,:,:,:)
!!$        integer, save :: ncalls = 0
!!$
!!$        ! This local CFD domain is outside MD overlap zone 
!!$        if ( md_map%n == 0 ) return 
!!$
!!$        call mpi_comm_rank(COUPLER_REALM_COMM, myid, ierr)
!!$
!!$        ncalls = ncalls + 1
!!$
!!$        ny = p3e - p3s + 1 ! number of y planes
!!$
!!$        min_i = minval(md_map%domains(1,:))
!!$        min_j = minval(md_map%domains(3,:))
!!$        min_k = minval(md_map%domains(5,:))
!!$
!!$        np = 0
!!$        do i = 1, md_map%n
!!$            np = np + 2 * (md_map%domains(2,i) - md_map%domains(1,i) + 1) &
!!$                * (md_map%domains(6,i) - md_map%domains(5,i) + 1) &
!!$                * ny
!!$
!!$            is(i) = md_map%domains(1,i) - min_i + 1
!!$            ie(i) = md_map%domains(2,i) - min_i + 1
!!$            ks(i) = md_map%domains(5,i) - min_k + 1
!!$            ke(i) = md_map%domains(6,i) - min_k + 1
!!$
!!$        enddo
!!$
!!$        allocate(vbuf(np),stat=ierr)
!!$        vbuf=0.d0
!!$
!!$        !write(0,*)'CFD recv_MDvel, vbuf size ', myid, np
!!$
!!$        start_address(1) = 1
!!$        do i = 1, md_map%n
!!$
!!$            source = md_map%rank_list(i)
!!$
!!$            np = 2 * (ke(i) - ks(i) + 1) * (ie(i) - is(i) + 1) * ny  ! +1 -> +0 for centered velocities
!!$
!!$            start_address(i+1) = start_address(i) + np
!!$
!!$            ! Attention ncall could go over max tag value for long runs!!
!!$            itag = mod(ncalls, MPI_TAG_UB)
!!$            call mpi_irecv(vbuf(start_address(i)),np, MPI_DOUBLE_PRECISION, source, itag, COUPLER_ICOMM, &
!!$                req(i),ierr)
!!$            !write(0,*) 'CFD recv_MDvel  ', myid, i, itag,source,np, is,ie,ks,ke,ierr	
!!$        enddo
!!$
!!$        call mpi_waitall(md_map%n, req, MPI_STATUSES_IGNORE, ierr)
!!$
!!$        allocate(vtot(2,nlz-0,nlx-0,ny),stat=ierr) ! for bulk uc 0->1
!!$        vtot=0.d0
!!$
!!$        do i = 1, md_map%n
!!$
!!$            if ( allocated(v1)) deallocate (v1)
!!$            allocate(v1(2,ks(i):ke(i)-0, is(i):ie(i)-0,ny))
!!$
!!$            v1(:,:,:,:) = reshape(vbuf(start_address(i):start_address(i+1)-1), &
!!$                (/ 2, ke(i)-ks(i)+1, ie(i)-is(i)+1, ny /))
!!$
!!$            vtot(:,ks(i):ke(i)-0, is(i):ie(i)-0,:) = vtot(:,ks(i):ke(i)-0, is(i):ie(i)-0,:) +  v1(:,:,:,:)
!!$
!!$        enddo
!!$
!!$
!!$        ! set k,i local indexes corresponding to the global ones
!!$        call mpi_comm_rank(coupler_grid_comm,myid,ierr)
!!$        gs1 = p1gs - bbox_cfd%zbb(1,icoord(3,myid+1)) + 1
!!$        gs2 = p2gs - bbox_cfd%xbb(1,icoord(1,myid+1)) + 1
!!$
!!$        ! Periodic boundary condition?	  
!!$        call set_pbc(pbc)
!!$
!!$
!!$        !OBSELETE check again and delete!!!
!!$        ! correction for z direction when used only for average ( 2d CFD)
!!$        !	if ( nlz > 2) then 
!!$        !		do j=1,ny
!!$        !		do i=1,nlx-1
!!$        !			x = sum(vtot(2,1:nlz-1,i,j))
!!$        !			do k=1,nlz-1
!!$        !				vtot(2,k,i,j) = x
!!$        !			enddo
!!$        !		enddo
!!$        !		enddo
!!$        !	endif
!!$
!!$
!!$
!!$        where (vtot(2,gs1:p1le-p1ls+gs1,gs2:p2le-p2ls+gs2,1:p3e-p3s+1) > 0.d0) 
!!$            vel(p1ls:p1le,p2ls:p2le,p3s:p3e) = vtot(1,gs1:p1le-p1ls+gs1,gs2:p2le-p2ls+gs2,1:p3e-p3s+1) &
!!$                /vtot(2,gs1:p1le-p1ls+gs1,gs2:p2le-p2ls+gs2,1:p3e-p3s+1)
!!$        elsewhere
!!$            vel(p1ls:p1le,p2ls:p2le,p3s:p3e) = 0.d0
!!$        end where
!!$
!!$    contains
!!$
!!$!-----------------------------------------------------------------------------
!!$! Sets periodic boundary conditions
!!$! adds together the sums of the first and last bins on the global grid (x and z directions)
!!$! and puts sum in both boxes
!!$! Note:
!!$! attention coord array are shifted with 1 in CFD and MD
!!$! to get correct result in MPI function such as mpi_cart_rank, subtract 1 
!!$!-----------------------------------------------------------------------------
!!$        subroutine set_pbc(pbc)
!!$            implicit none
!!$            integer, intent(in) :: pbc
!!$
!!$            real(kind(0.d0)), allocatable, dimension(:,:,:) :: x1, x2
!!$            integer dest, ip(3), ierr, status(MPI_STATUS_SIZE)
!!$
!!$            ! I guess that the PBC are apllied at the ends of global grid, check with CFD people
!!$            select case(pbc)
!!$            case(1)
!!$                ! first array dimension which corresponds to z direction vor velocity arrays 
!!$                allocate(x1(2,nlx,size(vtot,dim=4)),x2(2,nlx,size(vtot,dim=4)))
!!$                x1 =  vtot(:,1, 1:nlx,:)
!!$                x2 =  vtot(:,nlz, 1:nlx,:)
!!$                vtot(:, 1,   1:nlx, :) =   vtot(:, 1,   1:nlx, :) + x2
!!$                vtot(:, nlz, 1:nlx, :) =   vtot(:, nlz, 1:nlx, :) + x1
!!$            case(2)
!!$                ! second array dimension which correponds to x direction
!!$                allocate(x1(2,nlz,size(vtot,dim=4)),x2(2,nlz,size(vtot,dim=4)))
!!$                if ( npx == 1 )then 
!!$                    ! no MPI communication needed  
!!$                    write(0,*) 'coupler internal cfd pbc', npx
!!$                    x1 =  vtot(:, 1:nlz,1,  :)                 
!!$                    x2 =  vtot(:, 1:nlz,nlx,:)
!!$                    vtot(:, 1:nlz, 1,  :) =  vtot(:, 1:nlz, 1,  :) + x2
!!$                    vtot(:, 1:nlz, nlx,:) =  vtot(:, 1:nlz, nlx,:) + x1  
!!$                else 
!!$                    call mpi_comm_rank(coupler_grid_comm,myid,ierr)
!!$                    ip = icoord(: ,myid+1)
!!$                    if (ip(1) == 1 .and. ip(2) == 1) then 
!!$                        x1 =  vtot(:, 1:nlz, 1,:)
!!$                        call mpi_cart_rank(coupler_grid_comm, (/ npx-1, 0, ip(3)-1 /), dest, ierr)
!!$                        call mpi_sendrecv(x1,size(x1),MPI_DOUBLE_PRECISION,dest,1,x2,size(x2),MPI_DOUBLE_PRECISION,&
!!$                            dest,1,coupler_grid_comm,status,ierr)
!!$                        vtot(:, 1:nlz, 1, :) =  vtot(:, 1:nlz, 1, :) + x2
!!$                    else if (ip(1) == npx .and. ip(2) == 1) then 
!!$                        x2 =  vtot(:, 1:nlz, nlx,:)
!!$                        call mpi_cart_rank(coupler_grid_comm, (/ 0, 0, ip(3) - 1 /), dest, ierr)
!!$                        call mpi_sendrecv(x2,size(x2),MPI_DOUBLE_PRECISION,dest,1,x1,size(x1),MPI_DOUBLE_PRECISION,&
!!$                            dest,1,coupler_grid_comm,status,ierr)
!!$                        vtot(:, 1:nlz, nlx,:) =  vtot(:, 1:nlz, nlx,:) + x1
!!$                    endif
!!$                endif
!!$            end select
!!$        end subroutine set_pbc
!!$
!!$    end subroutine recv_vel_MD

end module coupler_internal_cfd
