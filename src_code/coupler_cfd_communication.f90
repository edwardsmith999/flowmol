module coupler_cfd_communication
        use coupler_cfd_global_data
        implicit none

contains

        subroutine md_vel(uc,vc,wc)
                use mpi
!                use data_export, only : ngz, i1_u, i2_u, j1_u, j2_u,&
!                        i1_v, i2_v, j1_v, j2_v, i1_w, i2_w, j1_w, j2_w
!                use messenger, only : myid, icoord
                use coupler_cfd_setup, only : nlx, nlz
                implicit none

                real(kind=kind(0.d0)),dimension(:,:,:),intent(out) :: uc, vc, wc 

                integer is, ie, iu_s, iu_e, iv_s, iv_e, iw_s, iw_e, myid
                integer  ierr

                call mpi_comm_rank(CFD_COMM, myid, ierr)
!                write(0,*) "CFD, md_vel: ", myid, i1_u, i2_u, i1_v, i2_v, i1_w, i2_w

                is = bbox_cfd%xbb(1,icoord(1,myid+1))
                ie = bbox_cfd%xbb(2,icoord(1,myid+1)) 
!                ks = bbox_cfd%zbb(1,icoord_xz(2,myid+1))
!                ke = bbox_cfd%zbb(2,icoord_xz(2,myid+1))

                iu_s = 1 !i1_u
                iv_s = 1 !i1_v
                iw_s = 1 !i1_w

                iu_e = iu_s + nlx - 1 - 1
                iv_e = iv_s + nlx - 1 - 1
                iw_e = iw_s + nlx - 1 - 1         

!                write(0,*) "CFD, md_vel 2: ", myid, iu_s, iu_e, iv_s, iv_e, iw_s, iw_e

                call  recv_vel_MD(uc, 1, 1, 1, nlx-1, 1, 1, 0)
                call  recv_vel_MD(vc, 1, 1, 1, nlx-1, 1, 1, 0)
!                call  recv_vel_MD(wc, 1, ngz-1, i1_w, i2_w, 0, 0, 1)

        end subroutine md_vel


        subroutine recv_vel_MD(a,p1s,p1e,p2s,p2e,p3s,p3e,pbc)
                use mpi
!                use data_export, only : jblock
                use coupler_cfd_setup, only : md_map, nlx, nlz, CFD_COMM_OVERLAP
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

                call mpi_comm_rank(CFD_COMM, myid, ierr)
                        
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
               
               v1(:,:,:,:) = reshape(vbuf(start_address(i):start_address(i+1)-1), (/ 2, ke(i)-ks(i)+0, ie(i)-is(i)+0, ny /))
               
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


        subroutine send_CFDvel(uc,vc)
                use mpi
                use coupler_cfd_setup, only : md_map, nlx, nly, nlz, CFD_COMM_OVERLAP
                implicit none
                
                real(kind=kind(0.d0)) uc(:,:),vc(:,:)

                real(kind=kind(0.d0)) vaux(nlz-1,nlx-1,nly-1), vbuf((nlz-1)*(nlx-1)*(nly-1))
                integer i, j,k, is, ie, js, je, ks, ke, iu_s, iu_e, iv_s, iv_e, iw_s, iw_e, &
                        min_i, min_j, min_k, np, myid, &
                        itag, dest, type, req(md_map%n), ierr
                integer status(MPI_STATUS_SIZE,md_map%n)
                integer, save :: ncalls = 0
                character(len=20) fname


! This local CFD domain is outside MD overlap zone 
                if ( md_map%n == 0 ) return 

               call mpi_comm_rank(CFD_COMM,myid,ierr)

                ncalls = ncalls + 1


!                write(0,*) "CFD, send_CFDvel: ", myid
                
                is = bbox_cfd%xbb(1,icoord(1,myid+1))
                ie = bbox_cfd%xbb(2,icoord(1,myid+1))
                
                iu_s = 2
                iv_s = 2
                iw_s = 1

                iu_e = iu_s + ie - is - 1 ! +1 - 1 !!
                iv_e = iv_s + ie - is - 1 
                iw_e = iw_s + ie - is - 1 
                
!                write(0,*) 'CFD send_CFDvel' , myid, ncalls, ie-is+1, nlx,nly,nlz

                min_i = minval(md_map%domains(1,:))
                min_j = minval(md_map%domains(3,:))
                min_k = minval(md_map%domains(5,:))

                do j = 1, 2

                        select case (j)
                        case (1)
                                vaux(1,:,:) = uc(iu_s:iu_e,2:jmax_overlap-1)
                        case (2)
                                vaux(1,:,:) = vc(iv_s:iv_e,2:jmax_overlap-1)
                        case (3) 
!                                vaux(:,:,:) = wc(1:1,iw_s:iw_e,1:jmax_overlap)
                        end select

                        do i = 1, md_map%n
                                dest = md_map%rank_list(i)
                                is = md_map%domains(1,i) - min_i + 1
                                ie = md_map%domains(2,i) - min_i + 1
                                js = md_map%domains(3,i) - min_j + 1
                                je = md_map%domains(4,i) - min_j + 1
                                ks = md_map%domains(5,i) - min_k + 1
                                ke = md_map%domains(6,i) - min_k + 1                                
                                np = (ie - is) * (je - js) * (ke - ks)

                                vbuf(1:np) = reshape(vaux(ks:ke-1,is:ie-1,js:je-1), (/ np /) )

	! Attention ncall could go over max tag value for long runs!!
                                itag = mod(j * ncalls, MPI_TAG_UB)
                                call mpi_send(vbuf, np, MPI_DOUBLE_PRECISION, dest, itag, CFD_MD_ICOMM, ierr)
!                                write(0,*) 'CFD sendCFD vel ', myid, j, i, itag,dest,np, is,ie,js,je,ks,ke,ierr        
                        enddo

!                        call mpi_waitall(md_map%n, req, MPI_STATUSES_IGNORE, ierr)

!                        write(fname,'(a,i0)') 'cfd_vel',3*(ncalls-1)+j

!                        call write_vector(vaux, fname, CFD_COMM_OVERLAP)

	!                                do i=1,novr
!                                 call mpi_get_count(status(1,i),mpi_double_precision,k,ierr)
	!                                 write(0,*) 'CFD send ', myid, ncalls, k, ' DP'
	!                                enddo

	!                                write(0,*) 'CFD send_CFDvel loop' , myid, j

!                                call mpi_barrier(CFD_BDY_MD_ICOMM,ierr)

                enddo
 

!                do k=1,nlz-1
!                        do j=2,jmax_overlap-1
!                                do i=iu_s,iu_e
!                                        write(5000+10*ncalls+myid,'(3E12.4)') uc(i,j), vc(i,j) !, wc(k,i,j)
!                                enddo
!                        enddo
!                enddo
                

!                write(0,*) 'CFD send_CFDvel finish' , myid, ncalls


	! This barrier does not work as this function si called inside an jblock == 1
	! condition. Is it really necessary?
!                call mpi_barrier(CFD_COMM, ierr)

	end subroutine send_CFDvel


	!
	!     debugging subroutines below
	!
!!$subroutine write_vector(q,fn, comm)
!!$	use mpi
!!$	use data_export, only : npx_cfd => npx,npz_cfd => npz
!!$	use messenger, only : myid, icoord
!!$	implicit none
!!$
!!$	real,dimension(:,:,:),intent(in) :: q
!!$	character(*),intent(in) :: fn
!!$	integer, intent(in) :: comm
!!$
!!$	integer  gsizes(3),lsizes(3),lstart(3),nt, &
!!$	ftype,fh,ierr, four_dp_type, ic, jc, kc, &
!!$	ibs, ibe, jbs,jbe, kbs, kbe, ibmin, &
!!$	ibmax, jbmin, jbmax, kbmin, kbmax, &
!!$	gs2(3), ls2(3), lb2(3), loc_suba
!!$	integer(kind=MPI_OFFSET_KIND) disp
!!$
!!$
!!$
!!$call mpi_comm_rank(comm, myid, ierr)
!!$
!!$	ic   = icoord(1,myid+1)
!!$	jc   = icoord(2,myid+1)
!!$kc   = icoord(3,myid+1)
!!$
!!$	ibs  = bbox_cfd%xbb(1,ic)
!!$	jbs  = bbox_cfd%ybb(1,jc)
!!$kbs  = bbox_cfd%zbb(1,kc)
!!$
!!$	ibe  = bbox_cfd%xbb(2,ic)
!!$	jbe  = bbox_cfd%ybb(2,ic) 
!!$kbe  = bbox_cfd%zbb(2,kc)
!!$
!!$	ibmin = bbox_cfd%xbb(1,1)
!!$	ibmax = bbox_cfd%xbb(2,npx_cfd)
!!$jbmin = bbox_cfd%ybb(1,1)
!!$	jbmax = jmax_overlap
!!$	kbmin = bbox_cfd%zbb(1,1)
!!$kbmax = bbox_cfd%zbb(2,npz_cfd)
!!$
!!$	gsizes(:) = (/  kbmax - kbmin+1, ibmax - ibmin+1, jbmax - jbmin + 1 /)
!!$!                psizes(:) = (/ npx_cfd, npz_cfd /)
!!$
!!$	lsizes(:) = (/ size(q,1) , size(q,2), size(q,3)  /)
!!$lstart(:) = (/  kbs - kbmin, (ic-1)*size(q,2), 0 /)
!!$
!!$	if (ic > 1) then
!!$	lsizes(2) = lsizes(2) -1
!!$	endif
!!$	if(ic > 2)then
!!$	lstart(2)=lstart(2)-ic+2
!!$	endif
!!$
!!$	disp = 0
!!$nt = lsizes(1) * lsizes(2) * lsizes(3)
!!$
!!$	write(0,*) 'CFD: write q', myid, gsizes, lsizes, lstart
!!$
!!$	gs2(:) = (/ size(q,1),size(q,2),size(q,3)  /)
!!$	ls2(:) = (/ size(q,1),size(q,2)-1,size(q,3) /)
!!$lb2(:) = (/ 0, 1, 0 /)
!!$
!!$	if (ic == 1)then
!!$ls2(2)=size(q,2)
!!$	lb2(2)=0
!!$	endif
!!$
!!$	!  three doubles to keep the speed value
!!$	!    call mpi_type_contiguous(4, MPI_DOUBLE_PRECISION, four_dp_type,ierr)
!!$!    call mpi_type_commit(four_dp_type, ierr)
!!$	four_dp_type = MPI_DOUBLE_PRECISION
!!$
!!$	call mpi_type_create_subarray(3, gsizes, lsizes, lstart, &
!!$			MPI_ORDER_FORTRAN, four_dp_type, ftype, ierr)
!!$call mpi_type_commit(ftype, ierr)
!!$	call mpi_type_create_subarray(3, gs2, ls2, lb2, &
!!$MPI_ORDER_FORTRAN, mpi_double_precision, loc_suba, ierr)
!!$call mpi_type_commit(loc_suba, ierr)
!!$
!!$                call mpi_file_open(comm, fn, &
!!$                        MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, fh, ierr)
!!$                call mpi_file_set_view(fh, disp, MPI_DOUBLE_PRECISION, ftype, &
!!$                        "native", MPI_INFO_NULL, IERR)
!!$                call mpi_file_write_all(fh, q, 1, loc_suba, &
!!$                        MPI_STATUS_IGNORE, ierr)
!!$                call mpi_file_close(fh, ierr)
!!$                call mpi_type_free(ftype, ierr)
!!$                call mpi_type_free(loc_suba,ierr)
!!$!                call mpi_type_free(four_dp_type,ierr)
!!$                
!!$                write(0,*) 'CFD: write q done', myid
!!$                
!!$        end subroutine write_vector
        
        

end module coupler_cfd_communication
