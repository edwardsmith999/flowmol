module coupler_cfd_communication
        use coupler_cfd_global_data
        implicit none
        
contains
        subroutine recv_vel_MD
                use mpi
                use data_export, only : jblock, ibmax_1, ibmin_1, kbmax_1, kbmin_1
                implicit none

                integer i, is, ie, ks, ke, ierr, source, myid, type, req(novr)
! for debug
                integer, save :: ncalls = 0
                character(len=128) :: fname

                if ( jblock == 1 ) then 

                        ncalls = ncalls + 1

                        call mpi_comm_rank(icomm_xz, myid, ierr)


                        is = ibmin_1(icoord_xz(1,myid+1))
                        ie = ibmax_1(icoord_xz(1,myid+1)) 
                        ks = kbmin_1(icoord_xz(2,myid+1))
                        ke = kbmax_1(icoord_xz(2,myid+1))

                        write(0,*) ' in recv_del ', myid, is, ie, ks, ke,icoord_xz(:,myid+1)

                        if ( .not. allocated(vel_fromMD)) then
                                allocate(vel_fromMD(3, is+1:ie,ks+1:ke))
                                vel_fromMD(:, is+1:ie,ks+1:ke) = -101010
                                write(fname,'(a,i0)') 'cfd_vel',0
                                call write_vector(vel_fromMD,fname, icomm_xz)
                        endif


                        do i = 1, novr
                                source =  map_overlap(i)%rank - 1
                                type = map_overlap(i)%dp2d_type
! Attention ncall could go over max tag value for long runs!!
                                call mpi_irecv(vel_fromMD(1,is+1,ks+1),1,type,source,ncalls,CFD_BDY_MD_ICOMM,req(i),ierr)
                        enddo

                        call mpi_waitall(novr, req, MPI_STATUSES_IGNORE, ierr)

!  debug writes  

                        write(fname,'(a,i0)') 'cfd_vel',ncalls

                        call write_vector(vel_fromMD,fname, icomm_xz)

                endif

                call mpi_barrier(CFD_COMM, ierr)

        end subroutine recv_vel_MD


!
!     debugging subroutines below
!
 subroutine write_vector(q,fn, comm)
  use mpi
  use data_export, only : npx_cfd => npx,npz_cfd => npz,ibmax_1, ibmin_1, kbmax_1, kbmin_1
   implicit none

  real,dimension(:,:,:),intent(in) :: q
  character(*),intent(in) :: fn
  integer, intent(in) :: comm

  integer  gsizes(2),psizes(2),lsizes(2),lstart(2),nt, &
   ftype,fh,ierr, three_dp_type, ic, kc, myid
  integer(kind=MPI_OFFSET_KIND) disp


 
  call mpi_comm_rank(comm, myid, ierr)

   gsizes(:) = (/ ibmax_1(npx_cfd) - ibmin_1(1), kbmax_1(npz_cfd) -kbmin_1(1)   /)
   psizes(:) = (/ npx_cfd, npz_cfd /)
   ic        = icoord_xz(1,myid+1)
   kc        = icoord_xz(2,myid+1)
   lsizes(:) = (/ ibmax_1(ic) - ibmin_1(ic), kbmax_1(kc) - kbmin_1(kc)  /)
   lstart(:) = (/ ibmin_1(ic) - ibmin_1(1), kbmin_1(kc) - kbmin_1(1) /)

   disp = 0
  nt = lsizes(1) * lsizes(2)

  write(0,*) 'write q', myid, gsizes, psizes, lsizes, lstart

!  three doubles to keep the speed value
    call mpi_type_contiguous(3, MPI_DOUBLE_PRECISION, three_dp_type,ierr)
    call mpi_type_commit(three_dp_type, ierr)

  call mpi_type_create_subarray(2, gsizes, lsizes, lstart, &
   MPI_ORDER_FORTRAN, three_dp_type, ftype, ierr)
  call mpi_type_commit(ftype, ierr)
  call mpi_file_open(comm, fn, &
   MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, fh, ierr)
  call mpi_file_set_view(fh, disp, MPI_DOUBLE_PRECISION, ftype, &
   "native", MPI_INFO_NULL, IERR)
  call mpi_file_write_all(fh, q, nt, three_dp_type, &
   MPI_STATUS_IGNORE, ierr)
  call mpi_file_close(fh, ierr)
  call mpi_type_free(ftype, ierr)
  call mpi_type_free(three_dp_type,ierr)


 end subroutine write_vector



end module coupler_cfd_communication
