module coupler_internal_md
        implicit none
        save 

	! a local mpi rank is useful
        integer myid

	! MD data        
        integer npx, npy, npz, nproc
        integer, allocatable :: icoord(:,:)
        real(kind=kind(0.d0)) :: dt_MD

	! data structures to hold CFD - MD mapping 

	! bounding box of MD domain within FD global domain
        type bbox_domain
                integer is,ie,js,je,ks,ke
                real(kind=kind(0.d0)) :: bb(2,3)
        end type bbox_domain

        type(bbox_domain), target :: bbox

	! thicknes of the MD region between the wall and CFD grid
        real(kind=kind(0.d0)) DY_PURE_MD

	! local domain lenghts, and halves
        real(kind=kind(0.d0)) domain_lengths(3), half_domain_lengths(3)

        type cfd_domain_map
                integer n ! number of CFD ranks that overlap with this MD domain
                integer, allocatable :: rank_list(:) ! rank list of overlpping CFD blocks
                integer, allocatable :: domains(:,:) ! indices range of overlapping 
                                                     ! with CFD grid 
        end type cfd_domain_map

        type(cfd_domain_map) cfd_map

	! write or not the overlap map
        logical :: dump_ovlerlap_map = .true.

	! local grid sizes
        integer nlx, nly, nlz

        type cfd_box_sum
                integer np
                real(kind=kind(0.d0))  v(3)
                real(kind=kind(0.d0))  a(3)
        end type cfd_box_sum

        real(kind(0.d0)), allocatable :: uc_bin(:,:,:), vc_bin(:,:,:,:), wc_bin(:,:,:) 

        real(kind=kind(0.d0)) :: FoP_time_ratio = 1.0   ! time ratio dt_CFD/dt_MD; to be fixed later
        real(kind=kind(0.d0)) :: xL_md, yL_md, zL_md 	! macroscopic sizes of MD domain. needed?
        real(kind=kind(0.d0)) :: fsig=1.0  		!Ratio betwen macroscopic unit lenght and molecular unit 

	! array for velocities from CFD, last dimension holds time indices 
        real(kind=kind(0.d0)), allocatable :: vel_fromCFD(:,:,:,:,:)
        integer itm1,itm2

	! data from CFD
        integer npx_cfd, npy_cfd, npz_cfd, jmax_overlap_cfd, nproc_cfd
	! CFD mesh data
        integer imino, imin_cfd, imax_cfd,imaxo, jmino, jmin_cfd, jmax_cfd, jmaxo,&
                kmino, kmin_cfd, kmax_cfd, kmaxo
        real(kind=kind(0.d0)), allocatable, target :: x(:), y(:), z(:)
        real(kind=kind(0.d0)) dx, dz, dt_CFD

	! nsteps from CFD
        integer nsteps
	! average period for averages ( it must come from CFD !!!)      
        integer :: average_period = 5
	! save period ( corresponts to tplot in CFD, revise please !!!)
        integer :: save_period = 10

contains 

        subroutine create_map_md
                use mpi
                use coupler_internal_common, only : CFD_MD_ICOMM
                implicit none
                integer  i, ir, ireq(nproc), noverlaps, ierr
                integer, allocatable :: overlap_mask(:,:)


		! compute the boundaries of this MD domain in the CFD global domain.
		! assume that all domains have identical sides 

                domain_lengths(:) = (/ xL_md/npx, yL_md/npy, zL_md/npz /)
                half_domain_lengths(:) = 0.5d0 * domain_lengths(:)

		! bounding boxes coordinates start from x(imin), z(kmin) and y(jmino)-DY_PURE_MD
                bbox%bb(1,:) = (icoord(:,myid+1)-1) * domain_lengths(:) &
                        + (/ x(imin_cfd), y(jmino)-DY_PURE_MD, z(kmin_cfd) /)
                bbox%bb(2,:) =  bbox%bb(1,:) + domain_lengths(:)

                !		write(0,*) 'MD: bbox%bb ', myid, bbox%bb !, domain, npx, npy, npz, xL_md, yL_md, zL_md, icoord_md

                call make_bbox

                !		write(0,*) 'MD: bbox%is ', myid, bbox%is, bbox%ie, bbox%js, bbox%je, bbox%ks, bbox%ke 
                !
		!  send box coordinate to CFD
                !

                call mpi_allgather((/ bbox%is, bbox%ie, bbox%js, bbox%je, bbox%ks, bbox%ke /), 6, MPI_INTEGER,&
                        MPI_BOTTOM,0,MPI_INTEGER,CFD_MD_ICOMM,ierr)

!!$		do i = 0, nproc_md-1
!!$
!!$		 if ( myid == i ) then
!!$		  source=MPI_ROOT
!!$		 else
!!$		  source=MPI_PROC_NULL
!!$		 endif
!!$		 
!!$		 call mpi_bcast((/ bbox%is, bbox%ie, bbox%js, bbox%je, bbox%ks, bbox%ke /), 6, MPI_INTEGER,&
!!$		  source, CFD_MD_ICOMM,ierr)
!!$
!!$		enddo


		! get the domain overlap mask from cfd

                allocate(overlap_mask(0:nproc-1,0:nproc_cfd-1))


                call mpi_allgather(MPI_BOTTOM,0, MPI_INTEGER,overlap_mask, &
                        nproc,MPI_INTEGER,CFD_MD_ICOMM,ierr)

                !		write(0,'(a,32I3)') 'MD, overlap mask: ', overlap_mask

                noverlaps = 0
                do i = 0, nproc_cfd - 1
                        if ( overlap_mask(myid,i) == 1) then 
                                noverlaps = noverlaps + 1
                        endif
                enddo

                !		write(0,'(a,32I3)') 'MD, noverlaps: ', myid, noverlaps

		! sort out which CFD ranks hold non-void domains for this MD rank

                cfd_map%n = noverlaps
                allocate ( cfd_map%rank_list(noverlaps), cfd_map%domains(6,noverlaps))

                ir=0
                do i=0, nproc_cfd - 1

                        if (overlap_mask(myid,i) == 1) then
                                ir = ir + 1
                                cfd_map%rank_list(ir) = i
                        endif
                enddo

                do i =1, ir
                        call mpi_irecv(cfd_map%domains(1,i), 6, MPI_INTEGER,cfd_map%rank_list(i),2, CFD_MD_ICOMM,ireq(i),ierr)
                enddo
                call mpi_waitall( ir,ireq,MPI_STATUSES_IGNORE,ierr)

                if ( dump_ovlerlap_map) then
                        call write_overlap_map
                endif

		! allocate array for CFD velocities and initialize the time indices; 
                allocate(vel_fromCFD(3,bbox%ie - bbox%is, bbox%je - bbox%js , bbox%ke - bbox%ks ,2))
                vel_fromCFD = 0.d0
                itm1 = 2; itm2 = 1

                !		write(0,*) 'MD: end of create_map_cfd_md', myid

        contains 

                subroutine make_bbox
                        implicit none

                        integer, parameter :: is = 1, ie = 2

                        type grid_pointer
                                real(kind=kind(0.d0)) , pointer :: p(:) => null()
                        end type grid_pointer

                        type bbox_pointer
                                integer, pointer :: start => null(), end => null()
                        end type bbox_pointer

                        type(grid_pointer) grid_ptr(3)
                        type(bbox_pointer) bbox_ptr(3)

                        integer id, ngp, grid_sizes(2,3), halo_size(2,3), idmin(3)
                        real(kind=kind(0.d0)) pl,pr,eps  ! left right grid points
                        logical found_start

			! indices covering the CFD physical domain

                        grid_sizes(:, :) = reshape((/ imin_cfd, imax_cfd, jmin_cfd, jmax_cfd, kmin_cfd, kmax_cfd /),(/2,3/))

			! starting indices (first in physical domain) in all directions

                        idmin = (/ imin_cfd, jmin_cfd, kmin_cfd /)

                        !		 write(0,*) 'MD: make box grid_sizes, idmin ', grid_sizes, idmin

			! how large is the halo in each direction, depending also on the MD domain position

                        halo_size(:,:) = 1

			! specical values for the boundaries

                        if ( icoord(2,myid + 1) == 1 ) then 

                                halo_size(1,2) = 0

                        endif

                        if ( icoord(2,myid + 1) == npy ) then

                                halo_size(2,2) = 0

                        endif

                        if (icoord(1, myid + 1) == 1) then

                                halo_size(1,1) = 0

                        endif

                        if (icoord(1, myid + 1) == npx) then 

                                halo_size(2,1) = 0

                        endif

                        if (icoord(3, myid + 1) == 1) then 

                                halo_size(1,3) = 0

                        endif

                        if (icoord(3, myid + 1) == npz) then 

                                halo_size(2,3) = 0

                        endif

			! pointer to grid coordinates. Helpful to loop through dimensions
                        grid_ptr(1)%p => x(imin_cfd:imax_cfd)
                        grid_ptr(2)%p => y(jmin_cfd:jmax_cfd)
                        grid_ptr(3)%p => z(kmin_cfd:kmax_cfd)

                        bbox_ptr(1)%start => bbox%is
                        bbox_ptr(1)%end   => bbox%ie
                        bbox_ptr(2)%start => bbox%js
                        bbox_ptr(2)%end   => bbox%je
                        bbox_ptr(3)%start => bbox%ks
                        bbox_ptr(3)%end   => bbox%ke


                        do id=1,3

                                eps = 1.d-2 * (grid_ptr(id)%p(2) - grid_ptr(id)%p(1))

                                !		  write(0,*) "MD: make box, grid step", myid, id, eps, grid_ptr(id)%p(1),bbox%bb(:,id)!   , grid_step !,x,y,z

                                found_start = .false.

                                if ( grid_ptr(id)%p(1) >= bbox%bb(1,id) - eps &
                                        .and. grid_ptr(id)%p(1) < bbox%bb(2,id) ) then 
                                        found_start = .true.
                                        bbox_ptr(id)%start = idmin(id) - halo_size(1,id)

                                        !		     write(0,*) "MD make box l", myid,id,bbox_ptr(id)%start
                                endif

                                ngp = grid_sizes(2,id) - grid_sizes(1,id) + 1

                                do i=2, ngp

                                        pl = grid_ptr(id)%p(i-1)
                                        pr = grid_ptr(id)%p(i) 

                                        !		    write(0,*), 'MD make bbox ', myid, id,ngp,i,pl,pr, bbox%bb(:,id) 

                                        if (.not. found_start )then
                                                if ( pl < bbox%bb(1,id)  .and. pr >=  bbox%bb(1,id) - eps ) then 
                                                        found_start = .true.


							! here comes the decision of how much one want to cover

                                                        bbox_ptr(id)%start = idmin(id) + i - 1 - halo_size(1,id)

                                                        !		     write(0,*), 'MD make bbox l', myid, id,i, pl, pr, bbox_ptr(id)%start
                                                endif

                                        else



                                                if ( (i < ngp  .and. pl <= bbox%bb(2,id) + eps  .and. pr > bbox%bb(2,id))) then

                                                        bbox_ptr(id)%end = idmin(id) + i - 1 -1 + halo_size(2,id)

                                                        !		    write(0,*), 'MD make bbox r', myid, id, i, pl, pr ,  bbox_ptr(id)%end		     
                                                        exit

                                                else if (i == ngp  .and. abs( pr - bbox%bb(2,id)) < eps ) then 

                                                        bbox_ptr(id)%end = idmin(id) + ngp - 1 + halo_size(2,id)

                                                        !		    write(0,*), 'MD make bbox r', myid, id,i, pl, pr, bbox_ptr(id)%end

                                                endif

                                        endif

                                enddo
                        enddo

! Set local grid sizes
                        nlx = bbox%ie - bbox%is + 1
                        nly = bbox%je - bbox%js + 1
                        nlz = bbox%ke - bbox%ks + 1

                end subroutine make_bbox

        end subroutine create_map_md


! get velocity fields from CFD for the force constrain needed in MD
        subroutine get_CFDvel
                use mpi
                use coupler_internal_common, only : CFD_MD_ICOMM
                implicit none

                real(kind=kind(0.d0)) vaux(nlz-1, nlx-1, nly-1)
                real(kind=kind(0.d0)), allocatable :: vbuf(:)
                integer i,j,k,id,ib,jb,kb,is,ie,js,je,ks,ke,np, start_address(cfd_map%n+1), &
                        di, dj, dk, itag, source, type, req(cfd_map%n), ierr
                integer status(MPI_STATUS_SIZE,cfd_map%n)
                character(len=MPI_MAX_ERROR_STRING) err_string, filename*20
                integer, save :: ncalls = 0

                ncalls = ncalls + 1

!                write(0,*) 'MD get_CFDvel' , myid, ncalls

!                call MPI_ERRHANDLER_SET(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr)

!               build derived types from cfd_map for data collections

!                write(0,*) ' MD_getCFDvel: nlz, nlx, nlz, cfd_map%n', nlz, nlx, nly, cfd_map%n

!              alocate the buffer, a bit too large, optimise later

                np = 0
                do i = 1, cfd_map%n
                        np = np + (cfd_map%domains(2,i) - cfd_map%domains(1,i)) &
                                *(cfd_map%domains(4,i) - cfd_map%domains(3,i)) &
                                *(cfd_map%domains(6,i) - cfd_map%domains(5,i))
                enddo


! two previous values of the CFD velocity field are needed for extrapolation
                i    = itm1  ! swap the time indices in vel_fromCFD
                itm1 = itm2
                itm2   = i

                allocate(vbuf(np),stat=ierr)

!                write(6000+10*ncalls+myid,'(3I4)') nlx,nly,nlz

!                write(0,*)'MD vbuf size ', myid, np

                do id = 1, 2

                        select case(id)
                        case(1)
                                di = 1; dj = 0; dk = 0
                        case(2)
                                di = 0; dj = 1; dk = 0
                        case(3)
                                di = 0; dj = 0; dk = 1
                        end select


                        start_address(1) = 1

                        do i = 1, cfd_map%n
                                source =  cfd_map%rank_list(i) 

                                is = cfd_map%domains(1,i) - bbox%is + 1
                                ie = cfd_map%domains(2,i) - bbox%is + 1
                                js = cfd_map%domains(3,i) - bbox%js + 1
                                je = cfd_map%domains(4,i) - bbox%js + 1
                                ks = cfd_map%domains(5,i) - bbox%ks + 1
                                ke = cfd_map%domains(6,i) - bbox%ks + 1
                                np = (ie - is + 0) * (je - js + 0) * (ke - ks +0)

                                start_address(i+1) = start_address(i) + np

! Attention ncall could go over max tag value for long runs!!
                                itag = mod(id * ncalls, MPI_TAG_UB)
                                call mpi_irecv(vbuf(start_address(i)),np,MPI_DOUBLE_PRECISION,source,itag,&
                                        CFD_MD_ICOMM,req(i),ierr)

!                                write(0,'(a,20(I4,x))') 'MD getCFD vel ',  myid, id, i, itag, source, np,is,ie,js,je,ks,ke, &
!                                        size(vbuf),start_address(i),ierr

                        enddo

                        call mpi_waitall(cfd_map%n, req, status, ierr)
!                        write(0,*) 'MD getCFD vel wait',  myid, id, i, source, ierr

                        do i=1,cfd_map%n

                                is = cfd_map%domains(1,i) - bbox%is + 1
                                ie = cfd_map%domains(2,i) - bbox%is + 1
                                js = cfd_map%domains(3,i) - bbox%js + 1
                                je = cfd_map%domains(4,i) - bbox%js + 1
                                ks = cfd_map%domains(5,i) - bbox%ks + 1
                                ke = cfd_map%domains(6,i) - bbox%ks + 1

                                vaux(ks:ke-1,is:ie-1,js:je-1) = reshape(vbuf(start_address(i):start_address(i+1)-1), &
                                        (/ ke-ks,ie-is,je-js /))

!                                 call MPI_Error_string(status(MPI_ERROR,i), err_string, len(err_string), ierr);
!                                  write(0,*) 'MD getCFD vel err, myid, i ', myid, i, trim(err_string) 
                                call mpi_get_count(status(1,i),mpi_double_precision,ib,ierr)
!                                 write(0,*) 'MD recv ', myid, id, i, ib, ' DP'
                        enddo

!                          write(2000+id,*) ks,ke,is,ie,js,je
!                          write(2000+id,*) vaux

!                         write(filename,'(a,i0,a,i0)') 'md_vel', 3*(ncalls-1),'_dim',id

!                         call write_vector(vaux(1,1,1),nlz,nlx,nly,filename, MD_COMM)
!                        call MPI_ERROR_STRING(status(MPI_ERROR,1),err_string,i,ierr)
!                        write(0,*) 'MD get_CFDvel error VC ', err_string(1:i)
!                         write(0,*) 'MD get_CFDvel loop' , myid, j                         


                        do kb = 1, nlz - 1 
                                do jb = 1, nly - 1 
                                        do ib = 1, nlx - 1 
                                                vel_fromCFD(id,ib,jb,kb,itm1) = vaux(kb,ib,jb) !&
!                                                          0.5d0 * (vaux(kb,ib,jb) + vaux(kb + dk,ib+ di,jb + dj))
!                                                   if ( id == 1) write(6000+10*ncalls+myid,'(3E12.4)')  vel_fromCFD(1,ib,jb,kb,itm1)
                                        enddo
                                enddo
                        enddo


                enddo ! id loop over dimensions

!                call flush(2001)
!                call flush(2002)
!                call flush(6000+10*ncalls+myid)

! average the velocity for the box 

! two previous values of the CFD velocity field are needed for extrapolation
!!$                i    = itm1  ! swap the time indices in vel_fromCFD
!!$                itm1 = itm2
!!$                itm2   = i
!!$
!!$                do kb = 1, nkb !bbox%ks+1, bbox%ke
!!$                        k = bbox%ks + kb
!!$                        do jb = 1, njb !bbox%js+1, bbox%je
!!$                                j = bbox%js + jb
!!$                                do ib = 1, nib !bbox%is+1, bbox%ie
!!$                                        i = bbox%is + ib
!!$                                        vel_fromCFD(1,ib,jb,kb,itm1) = 0.5d0 * (vaux(k,i-1,j,1) + vaux(k,i,j,1))
!!$                                        vel_fromCFD(2,ib,jb,kb,itm1) = 0.5d0 * (vaux(kb,ib,jb-1,2) + vaux(k,i,j,2))
!!$                                        vel_fromCFD(3,ib,jb,kb,itm1) = 0.5d0 * (vaux(k-1,i,j,3) + vaux(k,i,j,3))
!!$                                enddo
!!$                        enddo
!!$                enddo

! a halo of vel_fromCFD is needed in order to apply the continuum force constrains to the particles
! that have left the domain 
!                call exchange_vel_fromCFD_halo(itm1)x

!                write(0,*) 'MD get_CFDvel finish', myid, ncalls

        end subroutine get_CFDvel


        subroutine send_vel(a,n1,n2,n3)
                use mpi
                use coupler_internal_common
! send the average velocities to the corresponding rank in FD realm
                implicit none

                integer, intent(in) :: n1, n2, n3
                real(kind(0.d0)), intent(in) :: a(2,n1,n2,n3)

                integer i, is, ie, ks, ke, ierr, dest, type, itag, &
                        np, req(cfd_map%n), myid
                real(kind(0.d0)), allocatable :: buf(:)
! for debug
                integer, save :: ncalls = 0
                character(len=128) fname

                ncalls = ncalls + 1
                
                call mpi_comm_rank(COUPLER_COMM,myid,ierr)
!                write(0,*) 'MD: send_vel', myid, ncalls


                 do i = 1, cfd_map%n
                                dest =  cfd_map%rank_list(i) 
                                
                                is = cfd_map%domains(1,i) - bbox%is + 1
                                ie = cfd_map%domains(2,i) - bbox%is + 1
                                ks = cfd_map%domains(5,i) - bbox%ks + 1
                                ke = cfd_map%domains(6,i) - bbox%ks + 1
                                np = 2 * (ie - is + 0) * (ke - ks + 0) * n3

!                                write(0,*) 'MD: send_vel: dest, is, ie , ks, ke,..', myid, dest, is,ie,ks,ke

                                if( allocated(buf)) deallocate(buf)

                                allocate(buf(np), stat=ierr)

                                buf(:) = reshape( a(1:2,ks:ke-1,is:ie-1,1:n3), (/ np /)) 

                                 itag = mod( ncalls, MPI_TAG_UB)
                                call mpi_send(buf,np,MPI_DOUBLE_PRECISION,dest,itag,&
                                         CFD_MD_ICOMM,ierr)

                 end do

!                call mpi_waitall(cfd_map%n, req, MPI_STATUSES_IGNORE, ierr)

!                write(fname,'(a,i0)') 'md_vel',ncalls

!                call write_vector(a,n1,n2,n3,fname, MD_COMM)

        end subroutine send_vel

        
        function global_r(r) result(rg)
                implicit none

                real(kind(0.d0)),intent(in) :: r(3)
                real(kind(0.d0)) rg(3)

                rg(:) = r(:) + half_domain_lengths(:) + bbox%bb(1,:)

        end function global_r


	subroutine write_overlap_map
                use mpi
                use coupler_internal_common
		implicit none

		integer, parameter :: max_msg_length = 1024
		character(len=*),parameter :: filename = "coupler_map.log"
		character(len=*),parameter :: nl = achar(10) ! new line; portable?
		
		character(len=max_msg_length) msg, rec
		character(len=100) fmtstr
		
		integer fh, i, n, ierr

		call mpi_file_open(COUPLER_COMM,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)

		n = cfd_map%n

		write(rec,'(a,I5,a,I5)')  'myid', myid, ' overlaps', cfd_map%n

		msg = rec(1:len_trim(rec))//nl

		write(fmtstr,*) "( a,",1,"I5,","a,", 6,"I5)"
		
                do i = 1, n
                        write(rec,fmtstr) 'rank', cfd_map%rank_list(i), &
		                  ' indicies',cfd_map%domains(1:6,i)
		 
		        msg = msg(1:len_trim(msg))//rec(1:len_trim(rec))//nl
		enddo
		
		msg = msg(1:len_trim(msg))//nl
		n = len_trim(msg)

		call mpi_file_write_shared(fh,msg(1:n),n,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr) 
		
		call mpi_barrier(COUPLER_COMM,ierr)

		call mpi_file_close(fh,ierr)

	end subroutine write_overlap_map



end module coupler_internal_md
