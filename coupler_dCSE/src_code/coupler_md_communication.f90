module coupler_md_communication
        use coupler_md_global_data, x_cfd => x, y_cfd => y, z_cfd => z, dx_cfd => dx, dz_cfd => dz
        implicit none
        save
type box_average
  real(kind=kind(0.d0)) sum(3)
  integer np
 end type box_average

 type(box_average), allocatable ::  vel_bin(:,:)
 real(kind=kind(0.d0)), allocatable :: vel_average(:,:,:)
        

contains

        subroutine boundary_box_average
!  computes MD average velocity in a box of size dx*dy*dz around a FD grid point 
                use mpi
! MD modules
                use physical_constants_MD,      only : np
                use computational_constants_MD, only : domain, halfdomain
                use arrays_MD,                  only : r, v
                use messenger,               only : myid, icomm_grid, icoord

                implicit none

! tolerance for boundary overlap, 100th of sigma should do
                real(kind=kind(0.d0)), parameter :: epsilon=1.0d-2

                integer ixyz, displ, nbuff, nbuff_in, ib, kb, ic, kc, i, ip, nib, nkb
                integer source, dest, ierr, istat
! pointers to particles tags 
                integer pbuff(np)

! array used to transfer particle coordinates
                real(kind=kind(0.d0)), allocatable :: rvbuff(:,:), rvbuff_in(:,:)

! box size in MD units
                real(kind=kind(0.d0)) dx, dy, dz
! FD boundary position
                real(kind=kind(0.d0)) Y_boundary

! shifted average boxes ?
                logical lshift



! Ignoring cells and halos at the moment, to keep things simpler.


                if ( abs(shift_x) > epsilon  .or. abs(shift_z) > epsilon ) then
                        lshift = .true.
                else if ( abs(shift_x) > epsilon  .and. abs(shift_z) > epsilon )then
                        write(*,*) ' in MD boundary_box_average '
                        write(*,*) 'shift_x * shift_z /=  0'
                        write(*,*) 'This case is not implemented yet. Quitting ...'
                        call MPI_Abort(MPI_COMM_WORLD,1,ierr)
                else
                        lshift = .false.
                endif


                if ( shift_x < -epsilon ) then
                        ixyz = 1
                        displ = -1
                else if ( shift_x > epsilon ) then
                        ixyz = 1
                        displ = -1
                else if ( shift_z < -epsilon ) then
                        ixyz = 3
                        displ   = -1
                else if ( shift_z > epsilon ) then
                        ixyz = 3
                        displ   = 1
                endif

! sort the particles in the bin, mark the one that
! need to be tranferred in case of shifted bins

                ic = icoord(1,myid+1)
                nib = ibmax_md(ic)-ibmin_md(ic) ! +1 is missing,  number of boxes !!!
                ic = icoord(3,myid+1)
                nkb = kbmax_md(ic)-kbmin_md(ic)

                Y_boundary = Y_boundary_fraction*domain(2)

                if( .not. allocated(vel_bin)) then
                        allocate(vel_bin(nib,nkb),stat=istat)
                endif
                do kb = 1, nkb
                        do ib = 1, nib
                                vel_bin(ib, kb)%sum(:) = 0.0
                                vel_bin(ib, kb)%np    = 0
                        enddo
                enddo

                nbuff = 0

!  size of FD cell along y dimension, to revise later
                dy = fsig * (y_cfd(2) - y_cfd(1))

                dx = fsig * dx_cfd
                dz = fsig * dz_cfd

                write(0,*) 'dx, dy, dz ', myid, dx,dy,dz,fsig

                do ip = 1, np
!  find the box to which this particle belongs

! first: is the molecule close to the FD boundary plane?
                        if ( r(ip,2) + halfdomain(2) > Y_boundary + 0.5d0 * dy  &
                                & .or. r(ip,2) + halfdomain(2) < Y_boundary - 0.5d0 * dy) then 
                                cycle
                        endif

! Attention at the sign of the shift

                        ib = ceiling((r(ip,1) + halfdomain(1) - shift_x)/dx)

                        if ( lshift .and. (ib > nib .or. ib < 1)) then
! this molecule goes to the neighbors
                                nbuff = nbuff + 1 
                                pbuff(nbuff)=ip
                                cycle
                        endif

                        kb = ceiling((r(ip,3) + halfdomain(3) - shift_z)/dx)


                        if ( lshift .and. (kb < 1 .or. kb > nkb)) then
! this molecule goes to the neighbors
                                nbuff = nbuff + 1 
                                pbuff(nbuff)=ip
                                cycle
                        endif

! sum the velocities of  the other molecules    
                        vel_bin(ib,kb)%sum(:) = vel_bin(ib,kb)%sum(:) + v(ip, :)
                        vel_bin(ib,kb)%np     = vel_bin(ib,kb)%np + 1

                end do

!  debug, how many particles are there ?
                ip = 0
                do kb = 1, nkb
                        do ib = 1, nib

                                ip = ip + vel_bin(ib,kb)%np 

                        enddo
                enddo

                write(0,*) 'total number of particles to average over ;', myid, ip


                if( lshift) then
! sort out to which rank we have to send coordinates and from which one has to receive   


                        call mpi_cart_shift(icomm_grid, ixyz-1, displ, source, dest, ierr)

                        call mpi_sendrecv(nbuff, 1, mpi_integer, dest, 0, nbuff_in, 1, mpi_integer,source, 0,icomm_grid,&
                                MPI_STATUS_IGNORE, ierr)

! allocate arrays for position and velocities

                        allocate(rvbuff(6,nbuff), rvbuff_in(6,nbuff_in),  stat=istat)

                        do ip = 1, nbuff
                                i = pbuff(ip)
                                rvbuff(1:3,i) = r(i,1:3)
                                rvbuff(4:6,i) = v(i,1:3)
                        enddo

! send-receive the coordinates and position of particles involed in average 

                        call mpi_sendrecv(rvbuff, 6*nbuff, mpi_double_precision, dest, 0, rvbuff_in, 6*nbuff_in, &
                                mpi_double_precision,source, 0, icomm_grid, MPI_STATUS_IGNORE, ierr)

! correct positions in received buffers to the current domain
! this block can be done in a smarter way  
                        if ( shift_x < -epsilon ) then
! receiving from left
                                rvbuff_in(1,1:nbuff_in)  =  rvbuff_in(1,1:nbuff_in)-domain(1) 
                        else if ( shift_x > epsilon ) then
! reciving from right
                                rvbuff_in(1,1:nbuff)  =  rvbuff_in(1,1:nbuff_in)+domain(1)
                        else if ( shift_z < -epsilon ) then
! z direction left
                                rvbuff_in(3,1:nbuff)  =  rvbuff_in(3,1:nbuff_in)-domain(3) 
                        else if ( shift_z > epsilon ) then
! receiving from right
                                rvbuff_in(3,1:nbuff)  =  rvbuff_in(3,1:nbuff_in)+domain(3)
                        endif


! complete the average with received particles

                        do ip = 1, nbuff_in

                                ib = ceiling((rvbuff_in(1,ip) + halfdomain(1) - shift_x)/dx)

                                kb = ceiling((rvbuff_in(3,ip) + halfdomain(3) - shift_z)/dx)

! sum in the box an store the numer of particles as well molecules    
                                vel_bin(ib,kb)%sum(:) = vel_bin(ib,kb)%sum(:) + rvbuff_in(4:6, ip)
                                vel_bin(ib,kb)%np     = vel_bin(ib,kb)%np + 1

                        enddo

                endif ! lshift

! compute the averages    
                if ( .not. allocated(vel_average))then
                        allocate(vel_average(3,nib,nkb), stat=istat)
                endif

                do kb = 1, nkb
                        do ib = 1, nib

                                if ( vel_bin(ib,kb)%np > 0 ) then
                                        vel_average(:,ib,kb) =  vel_bin(ib,kb)%sum(:)/real(vel_bin(ib,kb)%np)
                                else
                                        vel_average(:,ib,kb) = 0.d0
                                endif


! synthetic data to test communication

!    ic = icoord(1,myid+1)
!    kc = icoord(3,myid+1)

!    vel_average(:,ib,kb) = (ibmin_md(ic) - ibmin_md(1)) * (kbmin_md(kc) - kbmin_md(1))  + (kb - 1)*nib +ib

                        enddo
                enddo

                call send_vel

        end subroutine boundary_box_average


        subroutine send_vel
                use mpi
! send the average velocities to the corresponding rank in FD realm
                implicit none

                integer i, ierr, dest, type, req(novr)
! for debug
                integer, save :: ncalls = 0
                character(len=128) fname

                ncalls = ncalls + 1

                do i = 1, novr
                        dest =  map_overlap(i)%rank - 1
                        type = map_overlap(i)%dp2d_type

! Attention ncall could go over max tag value for long runs!!
                        call mpi_isend(vel_average(1,1,1),1,type,dest,ncalls,CFD_BDY_MD_ICOMM,req(i),ierr)

                enddo

                call mpi_waitall(novr, req, MPI_STATUSES_IGNORE, ierr)

                write(fname,'(a,i0)') 'md_vel',ncalls

                call write_vector(vel_average,fname, MD_COMM)


        end subroutine send_vel


        !
!     debugging subroutines below
        !
        subroutine write_vector(q,fn, comm)
                use messenger, only : icoord_md => icoord
                use mpi
                implicit none

                real(kind=kind(0.d0)),dimension(:,:,:),intent(in) :: q
                character(*),intent(in) :: fn
                integer, intent(in) :: comm

                integer  gsizes(2),psizes(2),lsizes(2),lstart(2),nt, &
                        ftype,fh,ierr, three_dp_type, ic, kc, myid
                integer(kind=MPI_OFFSET_KIND) disp



                call mpi_comm_rank(comm, myid, ierr)


                gsizes(:) = (/ ibmax_md(npx_md) - ibmin_md(1), kbmax_md(npz_md) - kbmin_md(1)  /)
                psizes(:) = (/ npx_md, npz_md /)
                ic        = icoord_md(1,myid+1)
                kc        = icoord_md(3,myid+1)
                lsizes(:) = (/ ibmax_md(ic) - ibmin_md(ic),  kbmax_md(kc) - kbmin_md(kc) /)
                lstart(:) = (/  ibmin_md(ic) - ibmin_md(1),  kbmin_md(kc) - kbmin_md(1) /)

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


end module coupler_md_communication
