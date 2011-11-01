module coupler_md_communication
        use coupler_md_global_data
        implicit none
        save

 type cfd_box_sum
         integer np
         real(kind=kind(0.d0))  v(3)
         real(kind=kind(0.d0))  a(3)
 end type cfd_box_sum
!  real(kind=kind(0.d0)) sum(3)
!  integer np
! end type box_average

! type(box_average), allocatable ::  vel_bin(:,:)
! real(kind=kind(0.d0)), allocatable :: vel_average(:,:,:)

        real(kind(0.d0)), allocatable :: uc_bin(:,:,:), vc_bin(:,:,:,:), wc_bin(:,:,:) 

contains

        subroutine boundary_box_average(send_data)
!  computes MD average velocity in a box of size dx*dy*dz around a staggered FD grid point 
                use mpi
! MD modules
                use physical_constants_MD,      only : np
                use arrays_MD,                  only : r, v
                use messenger,                  only : myid, icomm_grid, icoord
                use coupler_md_setup,           only : dx, dz, y
                implicit none

                logical, intent(in) :: send_data

! tolerance for boundary overlap, 100th of sigma should do
                real(kind=kind(0.d0)), parameter :: epsilon=1.0d-2

                integer ixyz, displ, nbuff, nbuff_in, ib, kb, ic, kc, i, ip

!                write(0,*) 'md box-average: dx, dy, dz (sigma units), fsig', myid, dx,dz

! because of the staggered grid each a velocity componenet averahe
! must compute separately. Optimisation to be done later

                call compute_uc_average

                call compute_vc_average

!                call compute_wc_average

        contains 

                subroutine compute_uc_average
                        use coupler_md_setup, only : nlx, nlz, bbox, jmino, jmin => jmin_cfd
                        implicit none

                        integer ib, kb, ip, source, dest, ierr
                        real(kind=kind(0.d0)) rd(3)
                        logical, save :: first_time=.true.

                        if (first_time) then
                                first_time = .false.
                                allocate(uc_bin(2,nlz-1,nlx-1),stat=ierr)
                                uc_bin(:,:,:)  = 0.d0  
                        endif

                        do ip = 1, np
! using global particle coordinates
                                rd(:) = global_r(r(ip,:))

                                if ( rd(2) > y(jmin) .or. rd(2) < y(jmino) ) then
                                ! molecule outside the boundary layer
                                        cycle
                                endif

                                ib = ceiling((rd(1) - x(bbox%is)) / dx) + 0      ! staggered !!!
                                kb = ceiling((rd(3) - z(bbox%ks)) / dz)       ! the last z row unused    

                                if ( ib > 0 .and. ib <=  nlx  .and. &
                                     kb > 0 .and. kb <   nlz  ) then 
!  this particle are in this ranks domain
                                        uc_bin(1,kb,ib) = uc_bin(1,kb,ib) + v(ip,1)
                                        uc_bin(2,kb,ib) = uc_bin(2,kb,ib) + 1.d0 
                                else 
!                                       write(0,*) 'MD uc_average, outside domain rd', rd, ' bbox%bb ', bbox 
                                endif
                        end do

! debug   
!                         do i = 1, size(uc_bin,dim=2)
!                          write(0, '(a,I4,64F7.1)') 'MD myid uc_bin(2,..',myid,uc_bin(2,1,:)
!                         enddo



!                        write(0,*) 'MD uc sum in boxes', myid
!                        do i = 1, size(uc_bin,dim=2)
!                                write(0, '(a,I4,64E11.4)') 'MD myid uc_bin(1,..',myid, uc_bin(1,1,:)
!                        enddo
! send it to CFD        
                        if (send_data) then  
! testing                        uc_bin(1,:,:) = myid
!                                uc_bin(2,:,:) = 1.d0
                                call send_vel(uc_bin,nlz-1,nlx-1,1)   
                                uc_bin = 0.d0
                        endif

                end subroutine compute_uc_average


                subroutine compute_vc_average
                        use coupler_md_setup, only : nlx, nlz, bbox, jmino, jmin => jmin_cfd
                        implicit none
! this is the simplest one as there is no need the of data communication
                        integer  ib, jb, kb, ip, ierr
                        real(kind=kind(0.d0)) rd(3)
                        logical, save :: first_time = .true.

                        if( first_time ) then
                                first_time = .false.
                                allocate(vc_bin(2,nlz-1,nlx-1,1),stat=ierr)
                                vc_bin = 0.d0
                        endif

!                        write(2200+myid,*) bbox%bb

                        do ip = 1, np

! using global particle coordinates
                                rd(:) = global_r(r(ip,:))

!                                write(2200+myid,*) rd, y(jmino), y(jmin)

                                if (rd(2) >= y(jmino) .and. rd(2) <= y(jmin) ) then
                                        jb = 1
                                else
                                        cycle       
                                endif

! find the box indices

                                ib = ceiling((rd(1) - x(bbox%is)) / dx)
                                kb = ceiling((rd(3) - z(bbox%ks)) / dz)     

                                if ( ib > 0 .and. ib < nlx .and. &
                                     kb > 0 .and. kb < nlz ) then 
!  this particle are in this ranks domain
                                        vc_bin(1,kb,ib,jb) = vc_bin(1,kb,ib,jb) + v(ip,2)
                                        vc_bin(2,kb,ib,jb) = vc_bin(2,kb,ib,jb) + 1.d0
                                else 
!                                       write(0,*) 'MD vc_average, outside domain rd ', rd, ' bbox%bb ', bbox%bb 
                                endif
                        end do

! send to CFD
!                        write(0,*) 'MD, vc_average: got',  np, sum(vc_bin(2,:,:,:)), 'particles'

                        if (send_data ) then
! testing                        vc_bin(1,:,:,:) = 10*myid+1
!                                vc_bin(2,:,:,:) = 1.d0
                                call send_vel(vc_bin,nlz-1,nlx-1,1)    
                                vc_bin = 0.d0
                        endif

                end subroutine compute_vc_average



                subroutine compute_wc_average
                        use coupler_md_setup, only : nlx, nlz, bbox, bbox, jmin => jmin_cfd, jmino
                        implicit none

                        integer ib, kb, ip, source, dest, ierr
                        real(kind(0.d0)) rd(3)
                        logical, save :: first_time = .true.

                        if (first_time) then 
                                first_time = .false. 
                                allocate(wc_bin(2,nlz,nlx), stat = ierr)
                                wc_bin(:,:,:)  = 0.d0
                        endif

                        do ip = 1, np
! using global particle coordinates
                                rd(:) = global_r(r(ip,:))

                                if ( rd(2) > y(jmin) .or. rd(2) < y(jmino) ) then
                                ! molecule outside the boundary layer
                                        cycle
                                endif

                                ib = ceiling((rd(1) - x(bbox%is)) / dx)  
                                kb = nint((rd(3) - z(bbox%ks)) / dz) + 1 ! staggered    

                                if ( ib > 0 .and. ib <  nlx .and. &
                                     kb > 0 .and. kb <= nlz ) then 
!  this particle are in this ranks domain
                                        wc_bin(1,kb,ib) = wc_bin(1,kb,ib) + v(ip,3)
                                        wc_bin(2,kb,ib) = wc_bin(2,kb,ib) + 1.d0 
                                else 
!                                       write(0,*) 'MD wc_average, outside domain rd', rd, ' bbox%bb ', bbox%bb 
                                endif

                        enddo

! send to CFD
                        if (send_data) then
                                call send_vel(wc_bin,nlz,nlx,1)    
                                wc_bin(:,:,:) = 0.d0
                        endif

                end subroutine compute_wc_average

        end subroutine boundary_box_average


        subroutine send_vel(a,n1,n2,n3)
                use mpi
                use coupler_md_setup, only : cfd_map, bbox
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
                
                call mpi_comm_rank(MD_COMM,myid,ierr)
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


! get velocity fields from CFD for the force constrain needed in MD
        subroutine get_CFDvel
                use mpi
                use messenger, only : myid
                use coupler_md_setup, only : nlx, nly, nlz, bbox, cfd_map
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

                allocate(vbuf(np),stat=ierr)

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

!                         write(filename,'(a,i0,a,i0)') 'md_vel', 3*(ncalls-1),'_dim',id

!                         call write_vector(vaux(1,1,1),nlz,nlx,nly,filename, MD_COMM)
!                        call MPI_ERROR_STRING(status(MPI_ERROR,1),err_string,i,ierr)
!                        write(0,*) 'MD get_CFDvel error VC ', err_string(1:i)
!                         write(0,*) 'MD get_CFDvel loop' , myid, j

! two previous values of the CFD velocity field are needed for extrapolation
                         i    = itm1  ! swap the time indices in vel_fromCFD
                         itm1 = itm2
                         itm2   = i

                         do kb = 1, nlz - 1 
                                 do jb = 1, nly - 1 
                                         do ib = 1, nlx - 1 
                                                  vel_fromCFD(id,ib,jb,kb,itm1) = vaux(kb,ib,jb) !&
!                                                          0.5d0 * (vaux(kb,ib,jb) + vaux(kb + dk,ib+ di,jb + dj))
!                                                  write(6000+10*ncalls+myid,'(3E12.4)')  vel_fromCFD(id,ib,jb,kb,itm1)
                                         enddo
                                 enddo
                         enddo


                enddo ! id loop over dimensions


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


        subroutine simulation_apply_continuum_forces(iter)
                use arrays_MD,  only : r, v, a
                use physical_constants_MD, only : np
                use computational_constants_MD, only : halfdomain
                use coupler_md_global_data, only : x, y, z, dx, dz
                use coupler_md_setup, only : nlx, nly, nlz, dt_CFD,bbox
                use messenger, only : myid
                implicit none
                
                integer, intent(in) :: iter  ! iteration step, it assumes that each MD average
                                             ! start from iter = 1

                type(cfd_box_sum) :: box_average(bbox%ie - bbox%is, bbox%je - bbox%js, bbox%ke - bbox%ks)
                integer j, ib, jb, kb, nib, njb, nkb, ip, np_overlap
                integer list(4,np)
                real(kind=kind(0.d0)) inv_dtCFD, rd(3)
                integer :: ncalls = 0

!                write(0,*) 'MD: simulation_apply_continuum_forces', myid

! run through the particle, check if they are in the overlap region
! find the CFD box to which the particle belongs              
! attention to the particle that have left the domain boundaries 

! number of CFD cells in each direction

                nib = bbox%ie - bbox%is
                njb = bbox%je - bbox%js
                nkb = bbox%ke - bbox%ks

! What precisely is Y_boundary? It should be were the countinuum boundary is, the 
! averages for the boundary condition must be changed

                if (iter == 1) then
! get the previous value of CFD velocities
                    call  get_CFDvel
            endif

! at first CFD step we don't have two values to extrapolate CFD velocities, set inv_dtCFD=0
            if (ncalls == 0) then
                    inv_dtCFD = 0.0
            else
                    inv_dtCFD = 1.0/dt_CFD
            endif
            ncalls = ncalls + 1

            np_overlap = 0 ! number of particles in overlapping reg

            do kb = 1, nkb
                    do jb = 1, njb
                            do ib = 1, nib
                                    box_average(ib,jb,kb)%np   = 0
                                    box_average(ib,jb,kb)%v(:) = 0.0d0
                                    box_average(ib,jb,kb)%a(:) = 0.0d0
                            enddo
                    enddo
            enddo

                do ip = 1, np

! we need global MD coordinats to check if the particle is in the extended box.
! bbox%bb(:,:) are ok for they were used to build the MD domains
                 rd(:) = global_r(r(ip,:))

! struggling with the bottom boundary, bellow average boxes but with particles
!  for the moment let's work with 1 layer of MD blocks in 1 D
                        if ( rd(2) < y(bbox%js) .or.   rd(2) >= y(bbox%je) ) then
                                cycle 
                        else 
! non uniform grid in j direction                
                         jb=-666
                         do j =bbox%js+1, bbox%je
                                if( rd(2) <= y(j) ) then 
                                        !this is my cell, exit
                                        jb = j - bbox%js
                                        exit
                                endif
                          enddo

                        endif

! get the CFD cell coordinates   

                        
                        if (rd(1) < x(bbox%is) .or. rd(1) >= x(bbox%ie)) then
! this particle has left the domanin
!                                write(0,*) 'particle lost in x direction'
                                cycle
                        else 
                                ib = ceiling((rd(1) -  x(bbox%is))/ dx)
                        endif

                         
                        if (rd(3) < z(bbox%ks) .or. rd(3) >= z(bbox%ke) ) then
! this particle has left the domanin
!                                write(0,*) 'particle lost in z direction'
                                cycle
                        else
                                 kb = ceiling( (rd(3) - z(bbox%ks) ) / dz) 
                        endif

  

                        np_overlap = np_overlap + 1
                        list(1:4, np_overlap) = (/ ip, ib, jb, kb /)

                        box_average(ib,jb,kb)%np   =  box_average(ib,jb,kb)%np + 1
                        box_average(ib,jb,kb)%v(:) =  box_average(ib,jb,kb)%v(:) + v(ip,:)
                        box_average(ib,jb,kb)%a(:) =  box_average(ib,jb,kb)%a(:) + a(ip,:)

                enddo

! here we should have the cell coordinates for the particle ip which is 
! in the overlap region
! one has to tread separatley the particle that have left the domain
! compute the average force for each bin

!                write(0,*)'MD before average over bin'

!                call average_over_bin

!                write(0,*) 'MD: end simulation_apply_continuum_forces', myid

        contains

               subroutine average_over_bin
                        use computational_constants_MD, only : dt_MD => delta_t
                        use messenger, only : myid
                        implicit none

                        integer ib, jb, kb, i, ip, n
                        real(kind=kind(0.d0)) alpha(3), u_cfd_t_plus_dt(3), inv_dtMD


! set the continnum constrais for the particle in the bin

! speed extrapolation 

! add all up

                        inv_dtMD =1.d0/dt_MD

!                        write(400+10*ncalls+myid,'(a,I7,2E12.4)') "MD continuum np, vel_fromCFD1 : ", np_overlap, &
!                                                                  maxval(a(list(1,1:np_overlap),:)), &
!                                                                  minval(a(list(1,1:np_overlap),:))

                        do i = 1, np_overlap  
                                ip = list(1,i)
                                ib = list(2,i)
                                jb = list(3,i)
                                kb = list(4,i)

                                n = box_average(ib,jb,kb)%np

! use the following exptrapolation formula
! y = (y2-y1)/(x2-x1) * (x-x2) +y2

                                alpha(:) = inv_dtCFD*(vel_fromCFD(:,ib,jb,kb,itm1) - vel_fromCFD(:,ib,jb,kb,itm2))

                                u_cfd_t_plus_dt(:) = alpha(:) * (iter + 1)*dt_MD + vel_fromCFD(:,ib,jb,kb,itm1) 

                                a(ip,:) = a(ip,:) - box_average(ib,jb,kb)%a(:) / n - inv_dtMD * & 
                                          (box_average(ib,jb,kb)%v(:) / n - u_cfd_t_plus_dt(:))
                        enddo


!                        write(400+10*ncalls+myid,'(a,I7,2E12.4)') "MD continuum np, vel_fromCFD 2: ", np_overlap, &
!                                                                   maxval(a(list(1,1:np_overlap),:)), &
!                                                                   minval(a(list(1,1:np_overlap),:))
!                        write(400+10*ncalls+myid,'(a,2E12.4)')" inv_dtCFD, inv_dtMD ", inv_dtCFD, inv_dtMD
!                        do kb=1,nkb
!                        do jb=1,njb
!                        do ib=1,nib
!                        write(400+10*ncalls+myid,'(12E12.4,I7)') vel_fromCFD(:,ib,jb,kb,1), vel_fromCFD(:,ib,jb,kb,2),&
!                                                              box_average(ib,jb,kb)%v(:), box_average(ib,jb,kb)%a(:),&
!                                                              box_average(ib,jb,kb)%np
!                        enddo
!                        enddo
!                        enddo

end subroutine average_over_bin

end subroutine simulation_apply_continuum_forces


subroutine coupler_constrain_forces(np,pressure,r,a)
        use coupler_md_setup, only : bbox, jmax_overlap_cfd, halfdomain => half_domain_lengths
        implicit none 
        
        integer, intent(in)                  :: np
        real(kind=kind(0.d0)), intent(inout) :: a(:,:)
        real(kind=kind(0.d0)), intent(in)    :: pressure, r(:,:)

	! locals
        real(kind=kind(0.d0)), parameter :: eps = 0.d-2 ! avoid singular forces for molecules that 
                                                        ! happend to be very close to y(jmax_overlap)
        integer n
        real(kind=kind(0.d0)) p, yc, y2, y3
        
        
! Initial pressure is -ve - this line prevents problems        
        p = pressure
        if(pressure <= 0 ) then 
                p= 1.d0
        endif

        y2 = y(jmax_overlap_cfd -1)
        y3 = y(jmax_overlap_cfd) 

	do n = 1, np
! get the global value of y coordinate        
         yc  =  r(n,2) + halfdomain(2) + bbox%bb(1,2)

         if (yc  < y3 .and. yc >= y2  ) then
                 a(n,2)= a(n,2) - (yc-y2)/(1-(yc-y2)/(y3+eps-y2))*pressure
         endif
       enddo
        

end subroutine coupler_constrain_forces


function global_r(r) result(rg)
        use coupler_md_setup, only :  halfdomain => half_domain_lengths
        use coupler_md_setup, only : bbox
        implicit none
        
        real(kind(0.d0)),intent(in) :: r(3)
        real(kind(0.d0)) rg(3)
      
        rg(:) = r(:) + halfdomain(:) + bbox%bb(1,:)
        
end function global_r

subroutine exchange_vel_fromCFD_halo(ic)
!!$        USE mpi
!!$        USE messenger, only : icomm_grid, icoord
!!$        USE coupler_md_setup, only : nib,njb,nkb
!!$        implicit none
!!$
           integer, intent(in) :: ic
!!$
!!$        integer, parameter :: in=1, out=2, top=1, bottom=2, left=1, right=2
!!$
!!$        integer myid,sbx, sbz, idn, ide, ids, idw, irq(4), ierr
!!$        real(kind=kind(0.d0)) buff_x(3,njb,nkb,2,2), & ! exchange buffer along x
!!$                buff_z(3,0:nib+1,njb,2,2)              ! exchage buffer along z
!!$
!!$        sbx = 3 * njb * nkb
!!$        sbz = 3 * (nib + 2) * njb 
!!$
!!$        call mpi_comm_rank(icomm_grid,myid,ierr)
!!$
!!$! exchange along x
!!$
!!$        buff_x(:,:,:,left,out)  = vel_fromCFD(:,1,  1:njb,1:nkb,ic)
!!$        buff_x(:,:,:,right,out) = vel_fromCFD(:,nib,1:njb,1:nkb,ic)
!!$
!!$        call mpi_cart_rank(icomm_grid,icoord(1:3,myid+1)+(/-1-1,-1,-1/), idw, ierr)
!!$        call mpi_cart_rank(icomm_grid,icoord(1:3,myid+1)+(/-1+1,-1,-1/), ide, ierr)
!!$
!!$        call mpi_irecv(buff_x(1,1,1,left,in) , sbx,mpi_double_precision,idw,1,icomm_grid,irq(1),ierr)
!!$        call mpi_irecv(buff_x(1,1,1,right,in), sbx,mpi_double_precision,ide,1,icomm_grid,irq(2),ierr)
!!$
!!$        call mpi_isend(buff_x(1,1,1,left,out), sbx,mpi_double_precision,idw,1,icomm_grid,irq(3),ierr)
!!$        call mpi_isend(buff_x(1,1,1,right,out),sbx,mpi_double_precision,ide,1,icomm_grid,irq(4),ierr)
!!$
!!$       
!!$        
!!$        call mpi_waitall(4,irq,MPI_STATUSES_IGNORE,ierr)
!!$
!!$! exchage along z axis
!!$
!!$        buff_z(:,1:nib,:,top,   out) = vel_fromCFD(:,1:nib,1:njb,nkb,ic)
!!$        buff_z(:,1:nib,:,bottom,out) = vel_fromCFD(:,1:nib,1:njb,1  ,ic)
!!$
!!$         buff_z(:,0,    1:njb,top,out) = buff_x(:,1:njb,nkb,left, in)
!!$         buff_z(:,nib+1,1:njb,top,out) = buff_x(:,1:njb,nkb,right,in)
!!$
!!$         buff_z(:,0,    1:njb,bottom,out) = buff_x(:,1:njb,1,left,in)
!!$         buff_z(:,nib+1,1:njb,bottom,out) = buff_x(:,1:njb,1,right,in)         
!!$
!!$        call mpi_cart_rank(icomm_grid,icoord(1:3,myid+1)+(/-1,-1,-1+1/), idn, ierr)
!!$        call mpi_cart_rank(icomm_grid,icoord(1:3,myid+1)+(/-1,-1,-1-1/), ids, ierr)
!!$
!!$        call mpi_irecv(buff_z(1,0,1,top,in)   ,sbz,mpi_double_precision,idn,1,icomm_grid,irq(1),ierr)
!!$        call mpi_irecv(buff_z(1,0,1,bottom,in),sbz,mpi_double_precision,ids,1,icomm_grid,irq(2),ierr)
!!$ 
!!$        call mpi_isend(buff_z(1,0,1,top,out),   sbz,mpi_double_precision,idn,1,icomm_grid,irq(3),ierr)
!!$        call mpi_isend(buff_z(1,0,1,bottom,out),sbz,mpi_double_precision,ids,1,icomm_grid,irq(4),ierr)
!!$
!!$        call mpi_waitall(4,irq,MPI_STATUSES_IGNORE,ierr)
!!$
!!$! add the received velocities to vel_fromVFD array
!!$
!!$        vel_fromCFD(:, 0:nib+1, 1:njb, 0,    ic) = buff_z(:,0:nib+1,1:njb,bottom,in)
!!$        vel_fromCFD(:, 0:nib+1, 1:njb, nkb+1,ic) = buff_z(:,0:nib+1,1:njb,top   ,in)
!!$        vel_fromCFD(:, 0,       1:njb, 1:nkb,ic) = buff_x(:,1:njb,  1:nkb,left  ,in)
!!$        vel_fromCFD(:, nib+1,   1:njb, 1:nkb,ic) = buff_z(:,1:njb,  1:nkb,right ,in)
!!$

end subroutine exchange_vel_fromCFD_halo
! 
!     debugging subroutines below
! 
subroutine write_vector(q,n1,n2,n3,fn, comm)
!!$        use messenger, only : icoord_md => icoord
!!$        use mpi
!!$        implicit none
!!$        
        integer, intent(in) :: n1, n2, n3
        real(kind=kind(0.d0)),dimension(n1,n2,n3),intent(in) :: q
        character(*),intent(in) :: fn
        integer, intent(in) :: comm
!!$        
!!$        integer  gsizes(3),psizes(2),lsizes(3),lstart(3),nt, &
!!$                ftype,fh,ierr, three_dp_type, ic, kc, myid, &
!!$                gs2(3), ls2(3), lb2(3), loc_suba
!!$        integer(kind=MPI_OFFSET_KIND) disp
!!$        
!!$
!!$        
!!$call mpi_comm_rank(comm, myid, ierr)
!!$
!!$
!!$gsizes(:) = (/  bbox_md%zbb(2,npz_md) - bbox_md%zbb(1,1)+1, bbox_md%xbb(2,npx_md) - bbox_md%xbb(1,1)+1, n3  /)
!!$!psizes(:) = (/ npx_md, npz_md /)
!!$ic        = icoord_md(1,myid+1)
!!$kc        = icoord_md(3,myid+1)
!!$lsizes(:) = (/ n1, n2, n3 /)
!!$lstart(:) = (/ bbox_md%zbb(1,kc) - bbox_md%zbb(1,1), (ic-1)*n2, 0 /)
!!$
!!$if (ic > 1) then
!!$lsizes(2) = lsizes(2)-1
!!$endif
!!$
!!$if(ic > 2)then
!!$lstart(2)=lstart(2)-ic+2
!!$endif
!!$
!!$disp = 0
!!$nt = lsizes(1) * lsizes(2) * lsizes(3)
!!$
!!$write(0,*) 'MD write q', myid, gsizes, lsizes, lstart
!!$
!!$gs2(:) = (/ n1,n2,n3  /)
!!$ls2(:) = (/ n1,n2-1,n3 /)
!!$lb2(:) = (/ 0, 1, 0 /) 
!!$
!!$if (ic == 1)then
!!$ls2(2)=n2
!!$lb2(2)=0
!!$endif
!!$
!!$!  three doubles to keep the speed value 
!!$!                call mpi_type_contiguous(3, MPI_DOUBLE_PRECISION, three_dp_type,ierr)
!!$!                call mpi_type_commit(three_dp_type, ierr)
!!$three_dp_type = MPI_DOUBLE_PRECISION
!!$
!!$call mpi_type_create_subarray(3, gsizes, lsizes, lstart, &
!!$MPI_ORDER_FORTRAN, three_dp_type, ftype, ierr)
!!$call mpi_type_commit(ftype, ierr)
!!$
!!$call mpi_type_create_subarray(3, gs2, ls2, lb2, &
!!$MPI_ORDER_FORTRAN, mpi_double_precision, loc_suba, ierr)
!!$call mpi_type_commit(loc_suba, ierr)
!!$
!!$call mpi_file_open(comm, fn, &
!!$MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, fh, ierr)
!!$call mpi_file_set_view(fh, disp, MPI_DOUBLE_PRECISION, ftype, &
!!$"native", MPI_INFO_NULL, IERR)
!!$call mpi_file_write_all(fh, q, 1, loc_suba, &
!!$MPI_STATUS_IGNORE, ierr)
!!$call mpi_file_close(fh, ierr)
!!$call mpi_type_free(ftype, ierr)
!!$call mpi_type_free(loc_suba,ierr)
!!$!                call mpi_type_free(three_dp_type,ierr)

end subroutine write_vector


end module coupler_md_communication
