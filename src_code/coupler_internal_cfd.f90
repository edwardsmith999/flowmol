!=============================================================================
!				   Coupler internal CFD   				   
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
!=============================================================================


module coupler_internal_cfd
    implicit none
    save

    ! CFD data 

    ! CFD grid indices, number of processor in grid, number of steps
    integer imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,&
        	kmino,kmin,kmax,kmaxo,npx,npy,npz,nsteps

    integer :: tplot	       	!Output every number of steps
    integer :: jmax_overlap = 5 ! maximum j index ( in y direction) which MD 

    ! coordinates of CFD topologies
    ! ATTENTION the values are shifted with +1, FORTRAN style
    ! if this array is passed by an MPI function, remove the shift !
    integer, allocatable :: icoord(:,:)

    real(kind(0.d0)), allocatable :: x(:), y(:), z(:)
    real(kind(0.d0)) dx, dz, dt

    ! size of initialisation MD cell  used to resize CFD domain, derived from density
    ! and cell type
    real(kind(0.d0)) MD_initial_cellsize

    ! Data recieved from MD
    integer npx_md, npy_md, npz_md, nproc_md

    ! Internal data 
    ! bounding box type that holds on each processor the boundaries of the
    ! CFD grid sectors allocated to every processor
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

        integer 				:: i, myid, id_coord, color, noverlaps, ir, iaux(4), ierr
        integer, allocatable 	:: md_grid_boxes(:,:), overlap_mask(:), ireq(:), overlap_box(:,:)
        real(kind(0.d0))		:: raux(2)

        call mpi_comm_rank(COUPLER_REALM_COMM,myid,ierr)

        ! Coupler_grid_comm is the correct communicator to pick the processor coordinates in cartesian topology
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

		! Receive the overlapping box indices from all MD processors
        allocate(md_grid_boxes(6,0:nproc_md-1), overlap_mask(0:nproc_md-1), &
            	   overlap_box(6,0:nproc_md-1),         ireq(0:nproc_md-1))
		call mpi_allgather(MPI_BOTTOM,0,MPI_INTEGER,md_grid_boxes,6,MPI_INTEGER,COUPLER_ICOMM, ierr)
        !write(0,*) ' CFD grid boxes ', myid, md_grid_boxes
        call find_overlaps
        ! Send domain overlap mask to all MD processors
        call mpi_allgather(overlap_mask,nproc_md,MPI_INTEGER,MPI_BOTTOM,0,MPI_INTEGER,COUPLER_ICOMM,ierr)

        noverlaps = 0
        do i = 0, nproc_md - 1
            if ( overlap_mask(i) == 1) then 
                noverlaps = noverlaps + 1
            endif
        enddo

        map%n = noverlaps
        allocate ( map%rank_list(noverlaps), map%domains(6,noverlaps))

        !  Overlapping communicator
        if ( map%n > 0) then
            color = 1
        else 
            color = 0
        endif

        call mpi_comm_split(COUPLER_REALM_COMM, color, myid, CFD_COMM_OVERLAP, ierr)

        if (color == 0) then
            CFD_COMM_OVERLAP = MPI_COMM_NULL
        endif

        ! Send the range of the overlaping domains
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
        ! Make bbox which contains all domain information required in coupling process
        ! such as domain extents
        !-----------------------------------------------------------------------------

        subroutine make_bbox
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

			write(0,*)' CFD: bbox ', myid, bbox_cfd%xbb,bbox_cfd%ybb,bbox_cfd%zbb, nlx, nly, nlz

        end subroutine make_bbox

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

            write(0, *) 'CFD: find overlap ibmin etc', myid, ibmin, ibmax, jbmin, jbmax, kbmin, kbmax

            do i=0,nproc_md - 1

                ibs = md_grid_boxes(1,i)
                ibe = md_grid_boxes(2,i)
                jbs = md_grid_boxes(3,i)
                jbe = md_grid_boxes(4,i)
                kbs = md_grid_boxes(5,i)
                kbe = md_grid_boxes(6,i)

                if  ((( ibs <  ibmin .and. ibe > ibmin )	.or.  &
                    (   ibs >= ibmin .and. ibs < ibmax ))  .and. &
                    ((  jbs <  jbmin .and. jbe > jbmin )	.or.  &
                    (   jbs >= jbmin .and. jbs < jbmax))   .and. &   
                    ((  kbs <  kbmin .and. kbe > kbmin )	.or.  &
                    (   kbs >= kbmin .and. kbs < kbmax)))  then

					!This processor overlaps the MD domain
                    overlap_mask(i) = 1

					!I think this ensures local min/max are replaced by global min/max values if appropriate
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

            write(0,*)' CFD: overlap ', myid, overlap_mask, overlap_box

        end subroutine find_overlaps

    end subroutine create_map_cfd

end module coupler_internal_cfd
