program TransFlow
 use mpi
 use hybrid_map, only : MD_maps
 use messenger, only    : use_coupling_cfd=>use_coupling, CFD_COMM
 use messenger_md, only : use_coupling_md =>use_coupling, MD_COMM
 use data_export, only  : nproc_cfd => nproc
 use cfd_control, only  : main_init, main_restore, main_save, main_free 
 implicit none

 integer ierr
 integer i

 call create_communicators

! Initialise CDF 

 if (CFD_COMM /= MPI_COMM_NULL) then
! CFD initialisation
         
         use_coupling_cfd = .true.   
         
         call main_init()    

         write(0,*) 'after main init cfd'
    
 endif


! build MD maps using CFD information
 call MD_maps
 write(0,*) 'after MD_maps MD'

! set up MD 
 if (MD_COMM /= MPI_COMM_NULL) then
! MD initialisation
        use_coupling_md = .true.

        call setup_MD
        write(0,*) 'after setup MD'

 endif

 call mpi_barrier(MPI_COMM_WORLD,ierr)

 

 call map_cfd_md           ! establish communication with the MD cells 
                           ! that overlap the CFD boundary 

 call mpi_barrier(MPI_COMM_WORLD,ierr)

 call test_boundary_communication

! work starts
 call mpi_comm_rank(mpi_comm_world, i, ierr)
 write(0, *) ' before start work ', i
 
 if ( .true. ) then
! debugging above, jump this
                        
  if (CFD_COMM /= MPI_COMM_NULL) then
   
   call main_restore()       ! Read restart data
   call simulation_run()     ! Run the simulation
   call main_save()          ! Save the results
   call main_free()          ! Clean up
   
  endif
  
  if (MD_COMM /= MPI_COMM_NULL) then
        
   call simulation_MD
   call finish_MD 
   
  endif
  
 endif
 
 call mpi_barrier(MPI_COMM_WORLD,ierr)
 
 call messenger_free()


contains

 subroutine create_communicators
  use hybrid_map, only : CFD_MD_ICOMM,  CFD_BDY_MD_ICOMM

  implicit none

  integer ierr, myid, color, new_comm

  call mpi_init(ierr)
  
  call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  
  if ( myid < nproc_cfd) then
   color = 1
  else 
   color = 2   
  endif
  
  CFD_COMM = MPI_COMM_NULL
  MD_COMM  = MPI_COMM_NULL
  
  call mpi_comm_split(MPI_COMM_WORLD,color,myid,new_comm,ierr)

  if ( color == 1) then
   CFD_COMM = new_comm
  else if ( color == 2) then
   MD_COMM  = new_comm
  else
!    you shoudn't be here
   write(0,*) 'wrong color in Transflow create_communicators'
   call mpi_abort(MPI_COMM_WORLD,1,ierr)

  endif

! create inter-communicators

  if ( color == 1) then
! CFD side
   call mpi_intercomm_create(CFD_COMM,0,MPI_COMM_WORLD,nproc_cfd,1,CFD_MD_ICOMM,ierr)
 
  else
! MD side

   call mpi_intercomm_create(MD_COMM, 0, MPI_COMM_WORLD,        0, 1, CFD_MD_ICOMM, ierr)

  end if
  
  write(0,*) 'did communicators ', myid

 end subroutine create_communicators


 subroutine  map_cfd_md
  use hybrid_map, only      : map_overlap, novr, CFD_MD_ICOMM, CFD_BDY_MD_ICOMM, ibmin_md, ibmax_md, kbmin_md, &
   kbmax_md, icomm_xz, icoord_xz, npx_cfd, npz_cfd
  use messenger, only   : icomm_grid
  use data_export, only : jblock, ibmin, ibmax, kbmin, kbmax, &
   imin, imax, kmin, kmax,ibmap_1,kbmap_1,nix_1,niz_1,ibmin_1,ibmax_1, &
   kbmin_1, kbmax_1
  use mesh_export, only : x, z
  use computational_constants_MD, only : npx_md => npx, npz_md => npz, nproc_md => nproc
  use messenger_md, only : icomm_grid_md => icomm_grid, icoord_md => icoord
  implicit none

  integer i, j, k, ip, ierr, myid_cfd, myid_md, myid_xz_cfd, iroot, source, istat
  integer ibmin_loc, ibmax_loc, kbmin_loc, kbmax_loc, three_dp
  integer a_sizes(2), a_subsizes(2), a_starts(2)
  integer, allocatable :: overlap_aux(:,:)


  if ( CFD_COMM /= MPI_COMM_NULL) then
! CFD side

! build xz subcart
   call MPI_Cart_sub (icomm_grid, (/ .true., .false., .true. /), icomm_xz, ierr)

   call mpi_comm_rank(CFD_COMM,myid_cfd,ierr)

   if ( myid_cfd == 0 ) then
!  get the rank of first processor at the lower boundary
! a more general form is needed here
    call mpi_cart_rank(icomm_grid, (/ 0, 0, 0 /), iroot, ierr) 
    source=MPI_ROOT
   else
    source=MPI_PROC_NULL
   endif
 
   call mpi_bcast(iroot,1,mpi_integer,source,CFD_MD_ICOMM,ierr)
 
   if (jblock == 1) then
! intercommunicator between cfd boundary and md

    call mpi_intercomm_create(icomm_xz, 0, MPI_COMM_WORLD, nproc_cfd, 2, CFD_BDY_MD_ICOMM, ierr)
    
    
! get coordinates in the xz plane

    do ip=1,npx_cfd*npz_cfd
     call MPI_Cart_coords(icomm_xz, ip-1, 2, icoord_xz(1,ip), ierr)
    end do
    icoord_xz = icoord_xz + 1

!    broadcast icoord_xz(1)
     call mpi_comm_rank(icomm_xz,myid_xz_cfd,ierr)
     if( myid_xz_cfd == 0 ) then
      source = MPI_ROOT
     else
      source = MPI_PROC_NULL
     endif

    call mpi_bcast(icoord_xz, 2*npx_cfd*npz_cfd, mpi_integer, source,CFD_BDY_MD_ICOMM , ierr)

    call mpi_bcast(icoord_md, 3*nproc_md, mpi_integer, 0,CFD_BDY_MD_ICOMM , ierr)

    write(0,*) 'CFD side ', myid_xz_cfd, ' icoord_md ', icoord_md


!    bbox(1,1)=ibmap_1(1)
!    bbox(2,1)=ibmap_1(nix_1)
!    bbox(1,2)=kbmap_1(1)
!    bbox(2,2)=kbmap_1(niz_1)
!    write(0,*) 'cfd rank ', myid_cfd, 'bbox_cfd ', bbox, ' imin/max, kmin/max ', imin,imax,kmin,kmax, 'x(bbox) z(bbox) ', x(bbox(1,1)),x(bbox(2,1)),z(bbox(1,2)),z(bbox(2,2))

!!$     write(0,*) 'cfd rank ', myid_cfd,' ibmin_1 ', ibmin_1
!!$     write(0,*) 'cfd rank ', myid_cfd,' ibmax_1 ', ibmax_1
!!$     write(0,*) 'cfd rank ', myid_cfd,' kbmin_1 ', kbmin_1
!!$     write(0,*) 'cfd rank ', myid_cfd,' kbmax_1 ', kbmax_1
!!$
!!$     write(0,*) 'cfd rank ', myid_cfd,' ibmin_md ', ibmin_md
!!$     write(0,*) 'cfd rank ', myid_cfd,' ibmax_md ', ibmax_md
!!$     write(0,*) 'cfd rank ', myid_cfd,' kbmin_md ', kbmin_md
!!$     write(0,*) 'cfd rank ', myid_cfd,' kbmax_md ', kbmax_md
!!$    
!!$!    call mpi_comm_rank(icomm_xyz(1), iroot, ierr)
!!$
!!$    allocate(bbox_all(2,2,npx_md*npz_md))

! non-symmetric allgather is the solution to this transfer

!!$    call mpi_allgather(bbox, 4, MPI_INTEGER, bbox_all, 4*npx_md*npz_md, MPI_INTEGER, CFD_BDY_MD_ICOMM, ierr)

   endif

  else
! MD side

! ATTENTION Assumes at the moment that npy_md=1

   ! get the rank of jblock 1

   call mpi_bcast(iroot, 1, mpi_integer, 0, CFD_MD_ICOMM, ierr)

!  The MD part should be restricted to the processors that 
!  contain the CFD boundary; I'll see later to it

   call mpi_intercomm_create(icomm_grid_md, 0, MPI_COMM_WORLD,     iroot, 2, CFD_BDY_MD_ICOMM, ierr)

   call mpi_bcast(icoord_xz,2*npx_cfd*npz_cfd,mpi_integer,0,CFD_BDY_MD_ICOMM,ierr)

   

   call mpi_comm_rank(MD_COMM,myid_md,ierr)

   if( myid_md == 0 ) then
    source = MPI_ROOT
   else
    source = MPI_PROC_NULL
   endif

   write(0,*) 'MD side ', myid_md, ' icoord_xz ', icoord_xz
   
   call mpi_bcast(icoord_md, 3*nproc_md, mpi_integer, source,CFD_BDY_MD_ICOMM , ierr)
   

!!$   bbox(1,1)=bounding_box_md(1,1)
!!$   bbox(2,1)=bounding_box_md(1,1)
!!$
!!$   bbox(1,2)=bounding_box_md(1,3)
!!$   bbox(2,2)=bounding_box_md(2,3)
!!$   
!!$   allocate(bbox_all(2,2,npx_cfd))
!!$
!!$   call mpi_allgather(bbox, 4, MPI_INTEGER, bbox_all, 4*npx_cfd, MPI_INTEGER, CFD_BDY_MD_ICOMM,ierr)  
!!$   
!!$
!!$   call mpi_comm_rank(MD_COMM,myid_md,ierr)   
!!$   if ( myid_md == 0 ) then
!!$
!!$    write(0,*) ' MD bbox_all ', bbox_all
!!$    
!!$   endif

  end if
   
! sorting out the overlap betweeen CFD and MD boxes

  if ( CFD_COMM /= MPI_COMM_NULL ) then

   if  ( jblock == 1 ) then
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

  else  !MD branch

   allocate(overlap_aux(8,npx_cfd*npz_cfd))
   overlap_aux(:,:) = -1

!   three doubles for the speed values
   call mpi_type_contiguous(3, MPI_DOUBLE_PRECISION, three_dp,ierr)
   call mpi_type_commit(three_dp, ierr)

   call mpi_comm_rank(MD_COMM, myid_md, ierr)

   ibmin_loc = ibmin_md(icoord_md(1,myid_md+1))
   ibmax_loc = ibmax_md(icoord_md(1,myid_md+1))
   kbmin_loc = kbmin_md(icoord_md(3,myid_md+1))
   kbmax_loc = kbmax_md(icoord_md(3,myid_md+1))

   novr = 0
    do j= 1, npx_cfd*npz_cfd
     i = icoord_xz(1,j)
     k = icoord_xz(2,j)
      
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
      a_subsizes(:) = (/ overlap_aux(5,novr)-overlap_aux(4,novr), overlap_aux(7,novr)-overlap_aux(6,novr) /)
      a_starts(:) = (/ overlap_aux(4,novr) - ibmin_loc, overlap_aux(6,novr) - kbmin_loc/)
      call mpi_type_create_subarray(2, a_sizes, a_subsizes, a_starts, MPI_ORDER_FORTRAN, three_dp, overlap_aux(8,novr), ierr)
      call mpi_type_commit(overlap_aux(8,novr) , ierr)
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

   endif
   
 end subroutine map_cfd_md


 subroutine test_boundary_communication
  use hybrid_map, only      : map_overlap, CFD_MD_ICOMM, CFD_BDY_MD_ICOMM, icomm_xz, novr
  use data_export, only : jblock
  use messenger_md, only : icomm_grid_md => icomm_grid
  implicit none

  integer i,src,myid, ierr, bout(novr), bin(novr), req(novr,2)

    if ( CFD_COMM /= MPI_COMM_NULL ) then

     if  ( jblock == 1 ) then

      call mpi_comm_rank(icomm_xz,myid,ierr)

      do i =1, novr

       bout(i) = 100+myid

       src = map_overlap(i)%rank-1

       call mpi_irecv(bin(i),1,MPI_INTEGER,src,1,CFD_BDY_MD_ICOMM,req(i,1),ierr)
       
       call mpi_isend(bout(i),1,MPI_INTEGER,src,2,CFD_BDY_MD_ICOMM,req(i,2),ierr)
      enddo

      call mpi_waitall(2*novr,req,MPI_STATUSES_IGNORE,ierr)

      write(0,*) 'CFD ', myid, novr, map_overlap(1:novr)%rank, bout,bin
      
     endif

    else ! MD side
     
     call mpi_comm_rank(MD_COMM,myid,ierr)

     do i =1, novr
      
      bout(i) = 200+myid
      
      src = map_overlap(i)%rank-1
      
      call mpi_irecv(bin(i),1,MPI_INTEGER,src,2,CFD_BDY_MD_ICOMM,req(i,1),ierr)
      
      call mpi_isend(bout(i),1,MPI_INTEGER,src,1,CFD_BDY_MD_ICOMM,req(i,2),ierr)
      
     enddo
     
     call mpi_waitall(2*novr,req,MPI_STATUSES_IGNORE,ierr)
     
     
     
     write(0,*) 'MD ', myid, novr, map_overlap(1:novr)%rank, bout,bin

    end if
     

 end subroutine test_boundary_communication


end program TransFlow



