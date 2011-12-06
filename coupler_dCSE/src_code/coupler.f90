!=============================================================================
!
!                                  Coupler 
!
! Routines accessible from application ( molecular or continuum )
! after the name, in parenthesis, is the realm in which each routine must be called
!
! coupler_create_comm              (cfd+md) splits MPI_COMM_WORLD, create intercommunicator between CFD and MD
! coupler_get_cfd_info             (cfd)    brings CFD grid info and run parameters in the internal module
!					    coupler_internal_cfd                  
! coupler_get_md_info              (md)     does the same for MD in the module coupler_internal_md
! coupler_create_map               (cfd+md) creates correspondence maps between the CFD grid and MD domains
! coupler_constrain_forces         (md)     applies the the constrain force in the top Y domain in order to keep molecules inside 
! coupler_apply_continuum_forces   (md)     applies the continuum constrain force to molecules in a given sector of overlap 
!					    region (see Nie's et al article)
! coupler_uc_average_test          (md)     outputs uc component of Cartesian velocity in 4 for y layer centred around y=0 
! coupler_boundary_cell_average    (md)     computes averages of MD velocities and transfers them to CFD for the computation 
!					    of the boundary conditions 
! coupler_send_CFDvel              (cfd)    sends CFD velocities to MD for computation of the continuum constrains
! coupler_md_vel                   (cfd)    collects MD averages needed to compute boundary condition
! coupler_get                      (cfd+md) obtains coupler parameters from internal modules using a list of optional arguments
! coupler_get_save_period          (md)     returns save_period from coupler_internal_md
! coupler_get_average_period       (md)     returns average period from coupler_internal_md
! coupler_get_nsteps               (md)     returns nsteps (number of CFD simulation steps) from coupler_internal_md
! 
!  Lucian Anton, November 2011
!
!=============================================================================
module coupler
        use coupler_parameters
        implicit none
        save 

        ! dummy or active coupler ?
        logical, parameter :: COUPLER_IS_ACTIVE = .true.

contains

        subroutine coupler_create_comm(realm, realm_comm, ierror)
                use mpi
                use coupler_internal_common
                implicit none

                integer, intent(in) :: realm ! CFD or MD ?
                integer, intent(out):: realm_comm, ierror

                ! test if we have a CFD and a MD realm

                ierror=0

                call test_realms

                COUPLER_REALM = realm

                call create_comm

        contains

                subroutine create_comm
                        use mpi
                        use coupler_internal_common, only : COUPLER_COMM, CFD_MD_ICOMM
                        implicit none

                        integer ierr, myid, myid_comm, myid_comm_max, realm, &
                                iaux(2), jaux(2), remote_leader, comm, comm_size


                        realm = COUPLER_REALM

                        call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)

                        realm_comm   = MPI_COMM_NULL
                        COUPLER_COMM = MPI_COMM_NULL

                        ! get a communicator for each realm
                        call mpi_comm_split(MPI_COMM_WORLD,realm,myid,realm_comm,ierr)

                        ! get internal communicator for coupler operations
                        call mpi_comm_dup(realm_comm,COUPLER_COMM,ierr)

                        comm = COUPLER_COMM ! shorhand

                        ! create inter-communicators

                        ! Get the mpi_comm_world ranks that hold the largest ranks in cfd_comm and md_comm
                        call mpi_comm_rank(comm,myid_comm,ierr)
                        call mpi_comm_size(comm,comm_size,ierr)

                        iaux(:) = -1
                        jaux(:) = -1
                        if ( myid_comm == comm_size - 1) then
                                iaux(COUPLER_REALM) = myid
                        endif
                        call mpi_allreduce( iaux ,jaux, 2, MPI_INTEGER, MPI_MAX, &
                                MPI_COMM_WORLD, ierr)  

                        select case (COUPLER_REALM)
                        case (COUPLER_CFD)
                                remote_leader = jaux(COUPLER_MD)
                        case (COUPLER_MD)
                                remote_leader = jaux(COUPLER_CFD)
                        end select

                        !                        print*,color, jaux, remote_leader

                        call mpi_intercomm_create(comm, comm_size - 1, MPI_COMM_WORLD,&
                                remote_leader, 1, CFD_MD_ICOMM, ierr)

                        write(0,*) 'did (inter)communicators ', code_name(COUPLER_REALM), myid
                end subroutine create_comm


                subroutine test_realms
                        use mpi
                        implicit none

                        integer i, myid, nproc, nf, np, ierr
                        integer, allocatable :: ra(:)

                        call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)

                        if (myid == 0) then
                                call mpi_comm_size(mpi_comm_world, nproc, ierr)
                        	allocate(ra(nproc))
			else
				allocate(ra(1))	!Assign value arbitarily
                        endif

                        call mpi_gather(realm,1,MPI_INTEGER,ra,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

                        if (myid == 0) then

                                nf = 0
                                np = 0

                                do i =1, nproc
                                        if ( ra(i) == COUPLER_CFD ) then 
                                                nf = nf + 1
                                        else if ( ra(i) == COUPLER_MD ) then
                                                np = np +1
                                        else
                                                ierror = COUPLER_ERROR_REALM
                                                write(*,*) "wrong realm value in coupler_create_comm"
                                                call MPI_abort(MPI_COMM_WORLD,ierror,ierr)
                                        endif
                                enddo

                                if ( nf == 0 .or. np == 0) then 
                                        ierror = COUPLER_ERROR_ONE_REALM
                                        write(*,*) "CFD or MD  realm is missing in MPI_COMM_WORLD"
                                        call MPI_abort(MPI_COMM_WORLD,ierror,ierr)
                                endif



                        endif

                end subroutine test_realms

        end subroutine coupler_create_comm


        subroutine coupler_get_cfd_info(imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,&
                kmax,kmaxo,nsteps,x,y,z,dx,dz,npx,npy,npz,icoord,dt)
                use mpi
                use coupler_internal_common
                use coupler_internal_cfd, only : imino_ => imino, imin_ => imin, imax_ => imax, jmino_ => jmin, &
                        jmin_ => jmin, jmax_ => jmax, jmaxo_ => jmaxo, kmino_ => kmino, kmin_ => kmin, kmax_ => kmax, &
                        kmaxo_ => kmaxo, nsteps_ => nsteps, x_ => x, y_ => y, z_ => z, dx_ => dx, dz_ => dz, &
                        npx_ => npx, npy_ => npy, npz_ => npz, icoord_ => icoord, dt_ => dt, jmax_overlap, &
                        npx_md, npy_md, npz_md, nproc_md

                implicit none

                integer, intent(in) :: imino,imin,imax,imaxo,jmino,jmin,jmax,jmaxo,kmino,kmin,kmax,kmaxo,nsteps,&
                        npx,npy,npz,icoord(:,:)
                real(kind(0.d0)), intent(in)    :: x(:),y(:),z(:),dx,dz,dt

                integer i, ierr, myid, source, ngrid(3), ny_md, nblock(3), &
                        iaux(6)


                imino_ = imino; imin_ = imin; imax_ = imax; jmino_ = jmin
                        jmin_ = jmin; jmax_ = jmax; jmaxo_ = jmaxo; kmino_ = kmino; kmin_ = kmin; kmax_ = kmax
                        kmaxo_ = kmaxo; nsteps_ = nsteps;  dx_ = dx; dz_ = dz
                        npx_ = npx; npy_ = npy; npz_ = npz 

                allocate(x_(size(x)),stat=ierr); x_ = x
                allocate(y_(size(y)),stat=ierr); y_ = y
                allocate(z_(size(z)),stat=ierr); z_ = z
                allocate(icoord_(3,npx*npy*npz),stat=ierr)

                x_ = x; y_ = y; z_ = z;
                icoord_ = icoord

                !                 write(0,*) 'CFD exchange grid data'
                ! send CFD processor grid
                call mpi_comm_rank(COUPLER_COMM,myid,ierr)
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

        end subroutine coupler_get_cfd_info


        subroutine coupler_get_md_info(npxin,npyin,npzin,icoordin,dtin)
                use mpi
                use coupler_internal_common
                use coupler_internal_md
                implicit none

                integer, intent(in)          :: npxin, npyin, npzin, icoordin(:,:)
                real(kind(0.d0)), intent(in) :: dtin

                integer i, ierr, source, iaux(12)
                real(kind=kind(0.d0)) ra(2)


                call mpi_comm_rank(COUPLER_COMM,myid,ierr)

                npx = npxin; npy = npyin; npz = npzin
                nproc = npx * npy * npz

                dt_MD = dtin

                allocate(icoord(3,nproc),stat=ierr)

                icoord = icoordin

                ! write(0,*) 'MD exchange grid data'

                ! get CFD processor grid and the number of block in j direction
                call mpi_bcast( iaux, 4, MPI_INTEGER,&
                        0, CFD_MD_ICOMM,ierr)
                npx_cfd = iaux(1)
                npy_cfd = iaux(2)
                npz_cfd = iaux(3)
                jmax_overlap_cfd = iaux(4)
                nproc_cfd = npx_cfd * npy_cfd * npz_cfd

                ! send MD processor grid 
                if ( myid == 0 ) then
                        source=MPI_ROOT
                else
                        source=MPI_PROC_NULL
                endif

                call mpi_bcast((/ npx, npy, npz /), 3, MPI_INTEGER,&
                        source, CFD_MD_ICOMM,ierr)

                !! Test
                !	call mpi_comm_rank(MD_COMM,myid,ierr)!
                !
                !	write(0,*) 'exchange_grid_data: MD side', myid, bbox_cfd%xbb(1:2,1:npx_cfd),bbox_cfd%zbb(1:2,1:npz_cfd),&
                !	 bbox_md%xbb(1:2,1:npx_md),bbox_md%zbb(1:2,1:npz_md)    

                ! CFD mesh data 
                call mpi_bcast(iaux, 12, MPI_INTEGER, 0, CFD_MD_ICOMM,ierr) 
                imino = iaux(1); imin_cfd = iaux(2);  imax_cfd = iaux(3);  imaxo = iaux(4)
                jmino = iaux(5); jmin_cfd = iaux(6);  jmax_cfd = iaux(7);  jmaxo = iaux(8)
                kmino = iaux(9); kmin_cfd = iaux(10); kmax_cfd = iaux(11); kmaxo = iaux(12)
                allocate (x(imino:imaxo),y(jmino:jmaxo),z(kmino:kmaxo))


                call mpi_bcast(x,size(x),mpi_double_precision,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(y,size(y),mpi_double_precision,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(z,size(z),mpi_double_precision,0,CFD_MD_ICOMM,ierr)
                call mpi_bcast(ra,2,mpi_double_precision,0,CFD_MD_ICOMM,ierr)

                ! rescale all lengths to MD units to avoid disasters
                x = fsig * x; y = fsig * y; z = fsig * z 
                dx = fsig * ra(1); dz = fsig * ra(2)

                !		write(0,*) 'MD exchange grid data: imin0, imin ...', imino, imin_cfd, imax_cfd,imaxo, &
                !		 x(imino),x(imin_cfd), x(imax_cfd), x(imaxo)
                !		write(0,*) 'MD exchage grid data: recv dx, dz ', dx, dz, jmax_overlap_cfd

                ! get CFD nsteps
                call mpi_bcast(nsteps,1,mpi_integer,0,CFD_MD_ICOMM&
	 		&,ierr)

                ! get CFD dt
                call mpi_bcast(dt_CFD,1,mpi_double_precision,0&
	 		&,CFD_MD_ICOMM,ierr)
                ! should dt_CFD be scaled ?  see to it later
                dt_CFD = dt_CFD * FoP_time_ratio

                ! set the sizes of MD box

                DY_PURE_MD = 4.d0 * 5.12993d0 ! 4 nunits below CFD grid 
                !10.94378734741916534432d0 - (y(jmin_cfd) - y(jmino)) ! prototype
                !      
                xL_md = (x(imax_cfd) - x(imin_cfd))
                yL_md = (y(jmax_overlap_cfd) - y(jmino) +&
                        & DY_PURE_MD)		   
                zL_md = (z(kmax_cfd) - z(kmin_cfd))

                !	write(0,*) 'MD: exchange_grid... xL_md, yL_md, zL_md',&
                !	 & myid, xL_md, yL_md, zL_md
                !	write(0,*) 'MD: nsteps (CFD) ', nsteps



        end subroutine coupler_get_md_info


        subroutine coupler_create_map
                use mpi
                use coupler_internal_cfd
                use coupler_internal_md
                implicit none

                integer ierr

                if (COUPLER_REALM == COUPLER_CFD) then
                        call create_map_cfd
                else if (COUPLER_REALM .eq. COUPLER_MD) then
                        call create_map_md
                else
                        write(*,*) "Wrong COUPLER_REALM in coupler_create_map"
                        call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_REALM,ierr)
                end if

        end subroutine coupler_create_map


        subroutine coupler_constrain_forces(np,pressure,r,a)
                use coupler_internal_md, only : bbox, jmax_overlap_cfd, halfdomain => half_domain_lengths, y
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

                        if (yc  < y3 .and. yc >= y2 ) then
                                a(n,2)= a(n,2) - (yc-y2)/(1-(yc-y2)/(y3+eps-y2))*pressure
                        endif
                enddo


        end subroutine coupler_constrain_forces

        subroutine coupler_apply_continuum_forces(np,r,v,a,iter)
                use coupler_internal_common
                use coupler_internal_md, only : cfd_box_sum, halfdomain => half_domain_lengths, &
                                                x, y, z, dx, dz, global_r, jmin => jmin_cfd, &
                                                nlx, nly, nlz, dt_CFD,bbox, get_CFDvel, &
                                                jmax_overlap => jmax_overlap_cfd, myid
                implicit none

                real(kind=kind(0.d0)), dimension(:,:), intent(in) :: r,v
                real(kind=kind(0.d0)), dimension(:,:), intent(inout) :: a 
                integer, intent(in) :: np,iter  ! iteration step, it assumes that each MD average
                ! start from iter = 1

                type(cfd_box_sum) :: box_average(bbox%ie - bbox%is,1, bbox%ke - bbox%ks)
                integer j, ib, jb, kb, nib, njb, nkb, ip, np_overlap, jb_offset
                integer list(4,np)
                real(kind=kind(0.d0)) inv_dtCFD, rd(3)
                integer :: ncalls = 0

                ! run through the particle, check if they are in the overlap region
                ! find the CFD box to which the particle belongs              
                ! attention to the particle that have left the domain boundaries 

                ! number of CFD cells in each direction

                nib = bbox%ie - bbox%is
                njb = bbox%je - bbox%js
                nkb = bbox%ke - bbox%ks

                ! What precisely is Y_boundary? It should be were the countinuum boundary is, the 
                ! averages for the boundary condition must be changed

                if (iter .eq. 1) then
                        ! get the previous value of CFD velocities
                        call  get_CFDvel
                endif

                ! at first CFD step we don't have two values to extrapolate CFD velocities, set inv_dtCFD=0
                if (ncalls .eq. 0) then
                        inv_dtCFD = 0.0
                else
                        inv_dtCFD = 1.0/dt_CFD
                endif
                ncalls = ncalls + 1

                np_overlap = 0 ! number of particles in overlapping reg

                do kb = 1, ubound(box_average,dim=3)
                        do jb = 1, ubound(box_average,dim=2)
                                do ib = 1, ubound(box_average,dim=1)
                                        box_average(ib,jb,kb)%np   = 0
                                        box_average(ib,jb,kb)%v(:) = 0.0d0
                                        box_average(ib,jb,kb)%a(:) = 0.0d0
                                enddo
                        enddo
                enddo

                do ip = 1, np

                        ! we need global MD coordinats to check if the particle is in the extended box.
                        ! bbox%bb(:,:) are ok for they were used to build the MD domains
			rd(:) = r(ip,:)
                        rd(:) = global_r(rd)

                        ! struggling with the bottom boundary, below average boxes but with particles
                        !  for the moment let's work with 1 layer of MD blocks in 1 D
                        if ( rd(2) <= y(jmax_overlap-2) .or.   rd(2) >= y(jmax_overlap-1) ) then
                                cycle 
                        else 
                                ! non uniform grid in j direction                
                                jb=1
                                jb_offset=jmax_overlap-2-1
                                !                         do j =jmin+1, jmax_overlap
                                !                                if( rd(2) <= y(j) ) then 
                                !                                        !this is my cell index, exit
                                !                                        jb = j - jmin
                                !                                        exit
                                !                                endif
                                !                          enddo

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
                ! one has to treat separatley the particle that have left the domain
                ! compute the average force for each bin

                !                write(0,*)'MD before average over bin. np_overlap', np_overlap

                call average_over_bin

                !                write(0,*) 'MD: end simulation_apply_continuum_forces', myid

        contains

                subroutine average_over_bin
                        use coupler_internal_md, only : dt_MD, myid, itm1, itm2, vel_fromCFD 
                        implicit none

                        integer ib, jb, kb, i, ip, n
                        real(kind=kind(0.d0)) alpha(3), u_cfd_t_plus_dt(3), inv_dtMD, acfd


                        ! set the continnum constrais for the particle in the bin
                        ! speed extrapolation 
                        ! add all up

                        inv_dtMD =1.d0/dt_MD

                        !                        write(0,'(a,I7,2E12.4)') "MD continuum np, vel_fromCFD1 : ", np_overlap, &
                        !                                                                  maxval(a(list(1,1:np_overlap),:)), &
                        !                                                                  minval(a(list(1,1:np_overlap),:))

                        do i = 1, np_overlap  
                                ip = list(1,i)
                                ib = list(2,i)
                                jb = list(3,i)
                                kb = list(4,i)

                                n = box_average(ib,jb,kb)%np

                                !                                write(0,'(a,4I4,14E12.4)') "MD continuum force", ib,jb,kb,n,box_average(ib,jb,kb)%v(:), &
                                !                                        box_average(ib,jb,kb)%a(:),v(ip,:),a(ip,:),inv_dtMD,inv_dtCFD

                                if ( n .eq. 0 ) cycle
                                ! use the following exptrapolation formula
                                ! y = (y2-y1)/(x2-x1) * (x-x2) +y2

                                alpha(1) = inv_dtCFD*(vel_fromCFD(1,ib,jb+jb_offset,kb,itm1) - &
                                        vel_fromCFD(1,ib,jb+jb_offset,kb,itm2))

                                u_cfd_t_plus_dt(1) = alpha(1) * (iter + 1)*dt_MD + vel_fromCFD(1,ib,jb_offset,kb,itm1) 

                                acfd =  - box_average(ib,jb,kb)%a(1) / n - inv_dtMD * & 
                                        (box_average(ib,jb,kb)%v(1) / n - u_cfd_t_plus_dt(1))
                                a(ip,1) = a(ip,1) + acfd

                                !                                write(0,'(a,4I4,15E12.4)') "MD continuum force 2", ib,jb,kb,n, &
                                !                                 alpha(1),u_cfd_t_plus_dt(1),vel_fromCFD(1,ib,jb+jb_offset,kb,itm1),&
                                !                                 vel_fromCFD(1,ib,jb+jb_offset,kb,itm2),&
                                !                                 a(ip,1),acfd, r(ip,2) 

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

        end subroutine coupler_apply_continuum_forces



        subroutine coupler_uc_average_test(np,r,v,lwrite)
                use coupler_internal_common
                use coupler_internal_md, only : nlx, nlz, bbox, jmino, jmin => jmin_cfd,&
                                                global_r, x, dx, y, z, dz
                 
                implicit none

                integer, intent(in) :: np
                real(kind=kind(0.d0)), intent(in) :: r(:,:),v(:,:)
                logical, intent(in) :: lwrite

                integer ib, kb, jb, ip, myid, ierr
                real(kind=kind(0.d0)) rd(3), ymin, ymax, dy
                real(kind=kind(0.d0)),allocatable, save :: uc_bin(:,:,:,:)
                logical,save :: firsttime=.true.

                call mpi_comm_rank(COUPLER_COMM,myid,ierr)

                if(firsttime)then
                        firsttime = .false.
                        allocate(uc_bin(2,nlz-1,nlx-1,4))
                        uc_bin = 0.d0

                        if (myid .eq. 0) then 
                                open(45, file="md_vel.txt",position="rewind")
                                write(45,*)'# dx,dy,dz ', dx,y(jmin+1)-y(jmin),dz
                                close(45)
                        endif
                endif

                if (lwrite) then
                        call write_data
                        return
                endif


                dy = y(jmin+1) - y(jmin)
                ymin = y(jmin) - 2.d0 * dy
                ymax = y(jmin) + 2.d0 * dy 

                !        write(0,*)'MD uc test', np, dy, ymin,ymax

                do ip = 1, np
                        ! using global particle coordinates
			rd(:)=r(ip,:)
                        rd(:) = global_r(rd)

                        if ( rd(2) < ymin .or. rd(2) > ymax ) then
                                ! molecule outside the boundary layer
                                cycle
                        endif

                        ib = ceiling((rd(1) - x(bbox%is)) / dx) + 0      ! staggered !!!
                        kb = ceiling((rd(3) - z(bbox%ks)) / dz)       ! the last z row unused
                        jb = ceiling((rd(2) - ymin    )   / dy)

                        if ( ib > 0 .and. ib < nlx  .and. &
                                kb > 0 .and. kb < nlz  ) then 
                                !  this particle are in this ranks domain
                                uc_bin(1,kb,ib,jb) = uc_bin(1,kb,ib,jb) + v(ip,1)
                                uc_bin(2,kb,ib,jb) = uc_bin(2,kb,ib,jb) + 1.d0 
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

        contains 

                subroutine write_data
                        use mpi
                        use coupler_internal_md, only : nproc, imin_cfd, imax_cfd, kmin_cfd, kmax_cfd
                        implicit none

                        integer i, ibuff(2,2,0:nproc-1), ntot, nrecv, sa(nproc),req(nproc-1),  &
                                 ierr
                        real(kind(0.d0)),allocatable :: buff(:,:,:,:),buff_recv(:)

                        if(nproc > 1) then

                                ! works only for parallel decomposition in x and y direction
                                call mpi_gather((/bbox%is,bbox%ie,bbox%ks,bbox%ke/),4,MPI_INTEGER,&
                                        ibuff,4,MPI_INTEGER,0,COUPLER_COMM,ierr)

                                !               write(0,*) "MD write test data", myid, ibuff

                                if (myid .eq. 0) then

                                        ! the local bit first
                                        allocate(buff(2,kmax_cfd-kmin_cfd,imin_cfd:imax_cfd-1,4))

                                        buff = 0.d0
                                        buff(:,ibuff(1,2,0):ibuff(2,2,0)-1,ibuff(1,1,0):ibuff(2,1,0)-1,:) = &
                                                buff(:,ibuff(1,2,0):ibuff(2,2,0)-1,ibuff(1,1,0):ibuff(2,1,0)-1,:) + &
                                                uc_bin(:,1:nlz-1,1:nlx-1,:)


                                        ntot = 0
                                        do i=1,nproc-1
                                                ntot = ntot + 2*(ibuff(2,2,i)-ibuff(1,2,i))*(ibuff(2,1,i)-ibuff(1,1,i))*4
                                        enddo

                                        allocate(buff_recv(ntot))
                                        buff_recv(ntot) = 0.d0

                                        sa(1)=1

                                        do i=1,nproc-1

                                                nrecv = 2*(ibuff(2,2,i)-ibuff(1,2,i))*(ibuff(2,1,i)-ibuff(1,1,i))*4
                                                sa(i+1) = sa(i) + nrecv

                                                call mpi_irecv(buff_recv(sa(i)),nrecv,MPI_DOUBLE_PRECISION,&
                                                        i,1,COUPLER_COMM,req(i),ierr)

                                        enddo

                                        call mpi_waitall(nproc-1,req,MPI_STATUSES_IGNORE,ierr)

                                        do i =1, nproc-1
                                                buff(:,ibuff(1,2,i):ibuff(2,2,i)-1,ibuff(1,1,i):ibuff(2,1,i)-1,:) = &
                                                        buff(:,ibuff(1,2,i):ibuff(2,2,i)-1,ibuff(1,1,i):ibuff(2,1,i)-1,:) + &
                                                        reshape(buff_recv(sa(i):sa(i+1)-1), (/ 2,ibuff(2,2,i)-ibuff(1,2,i),ibuff(2,1,i)-ibuff(1,1,i),4 /))
                                        enddo
                                else

                                        call mpi_send(uc_bin,size(uc_bin),MPI_DOUBLE_PRECISION,0,1,COUPLER_COMM,ierr)

                                 endif
                        endif



                        if (nproc > 1) then
                                if (myid .eq. 0 ) then
                                        open(45,file="md_vel.txt",position="append")
                                        do i = 1,4
                                                write(45, '(100(E12.4,1x))') sum(buff(:,:,:,i),dim=2) 
                                        enddo
                                        write(45, '(1x/1x)')
                                        close(45)
                                endif
                        else
                                open(45,file="md_vel.txt",position="append")
                                do i = 1,4
                                        write(45, '(100(E12.4,1x))') sum(uc_bin(:,:,:,i),dim=2)
                                enddo
                                write(45, '(1x/1x)')
                                close(45)
                        endif



                        uc_bin = 0.d0



                end subroutine write_data

        end subroutine coupler_uc_average_test


        subroutine coupler_boundary_cell_average(np,r,v,send_data)
                !  computes MD average velocity in a box of size dx*dy*dz around a staggered FD grid point 
                use coupler_internal_md, only : myid, icoord, dx, dz, x, y, z, global_r, send_vel
                implicit none

                integer, intent(in) :: np
                real(kind=kind(0.d0)), dimension(:,:), intent(in) :: r,v
                logical, intent(in) :: send_data

                ! tolerance for boundary overlap, 100th of sigma should do
                real(kind=kind(0.d0)), parameter :: epsilon=1.0d-2

                integer ixyz, displ, nbuff, nbuff_in, ib, kb, ic, kc, i, ip

                !                write(0,*) 'md box-average: dx, dy, dz (sigma units), fsig', myid, dx,dz

                ! because of the staggered grid each a velocity component average
                ! must compute separately. Optimisation to be done later
                call compute_uc_average
                call compute_vc_average

                !                call compute_wc_average

        contains 

                subroutine compute_uc_average
                        use coupler_internal_md, only : nlx, nlz, bbox, jmino, jmin => jmin_cfd,&
                                                        uc_bin
                        implicit none

                        integer ib, kb, ip, source, dest, ierr
                        real(kind=kind(0.d0)) rd(3)
                        logical, save :: first_time=.true.

                        if (first_time) then
                                first_time = .false.
                                allocate(uc_bin(2,nlz-1,nlx-1),stat=ierr)
                                uc_bin(:,:,:)  = 0.d0  
                        endif

                        ! send it to CFD        
                        if (send_data) then  
                                ! testing                        uc_bin(1,:,:) = myid
                                !                                uc_bin(2,:,:) = 1.d0
                                call send_vel(uc_bin,nlz-1,nlx-1,1)   
                                uc_bin = 0.d0

                                return
                        endif


                        do ip = 1, np
                                ! using global particle coordinates
				rd(:) = r(ip,:)
                                rd(:) = global_r(rd)

                                if ( rd(2) > y(jmin) .or. rd(2) < y(jmino) ) then
                                        ! molecule outside the boundary layer
                                        cycle
                                endif

                                ib = ceiling((rd(1) - x(bbox%is)) / dx) + 0      ! staggered !!!
                                kb = ceiling((rd(3) - z(bbox%ks)) / dz)       ! the last z row unused    

                                if ( ib > 0 .and. ib <  nlx  .and. &
                                        kb > 0 .and. kb <  nlz  ) then 
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

                end subroutine compute_uc_average


                subroutine compute_vc_average
                        use coupler_internal_md, only : nlx, nlz, bbox, jmino, jmin => jmin_cfd, vc_bin
                        implicit none
                        ! this is the simplest one as there is no need for data communication
                        integer  ib, jb, kb, ip, ierr
                        real(kind=kind(0.d0)) rd(3)
                        logical, save :: first_time = .true.

                        if( first_time ) then
                                first_time = .false.
                                allocate(vc_bin(2,nlz-1,nlx-1,1),stat=ierr)
                                vc_bin = 0.d0
                        endif


                        if (send_data ) then
                                ! testing                        vc_bin(1,:,:,:) = 10*myid+1
                                !                                vc_bin(2,:,:,:) = 1.d0
                                call send_vel(vc_bin,nlz-1,nlx-1,1)    

                                vc_bin = 0.d0
                                return
                        endif


                        !                        write(2200+myid,*) bbox%bb

                        do ip = 1, np

                                ! using global particle coordinates
				rd(:) = r(ip,:)
                                rd(:) = global_r(rd)

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
                                        !write(0,*) 'MD vc_average, outside domain rd ', rd, ' bbox%bb ', bbox%bb 
                                endif
                        end do

                        ! send to CFD
                        !                        write(0,*) 'MD, vc_average: got',  np, sum(vc_bin(2,:,:,:)), 'particles'


                end subroutine compute_vc_average



                subroutine compute_wc_average(np)
                        use coupler_internal_md, only : nlx, nlz, bbox, bbox, jmin => jmin_cfd, jmino, wc_bin
                        implicit none

                        integer, intent(in) :: np  ! particle number, it might chage from step to step

                        integer ib, kb, ip, source, dest, ierr
                        real(kind(0.d0)) rd(3)
                        logical, save :: first_time = .true.

                        if (first_time) then 
                                first_time = .false. 
                                allocate(wc_bin(2,nlz,nlx), stat = ierr)
                                wc_bin(:,:,:)  = 0.d0
                        endif

                        do ip = 1, np
                                ! use global particle coordinates
				rd(:) = r(ip,:)
                                rd(:) = global_r(rd)

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
                                        !write(0,*) 'MD wc_average, outside domain rd', rd, ' bbox%bb ', bbox%bb 
                                endif

                        enddo

                        ! send to CFD
                        if (send_data) then
                                call send_vel(wc_bin,nlz,nlx,1)    
                                wc_bin(:,:,:) = 0.d0
                        endif

                end subroutine compute_wc_average

        end subroutine coupler_boundary_cell_average


        subroutine coupler_send_CFDvel(uc,vc)
                use mpi
                use coupler_internal_cfd, only : md_map, nlx, nly, nlz, CFD_COMM_OVERLAP, &
                                                 bbox_cfd, jmax_overlap, dx, dz, icoord
                use coupler_internal_common
                implicit none

                real(kind=kind(0.d0)) uc(:,:,:),vc(:,:,:)

                real(kind=kind(0.d0)) vaux(nlz-1,nlx-1,nly-1), vbuf((nlz-1)*(nlx-1)*(nly-1))
                integer i, j,k, is, ie, js, je, ks, ke, iu_s, iu_e, iv_s, iv_e, iw_s, iw_e, &
                        ku_s, ku_e, kv_s, kv_e, min_i, min_j, min_k, np, myid, &
                        itag, dest, type, req(md_map%n), ierr
                integer status(MPI_STATUS_SIZE,md_map%n)
                integer, save :: ncalls = 0
                character(len=20) fname


                ! This local CFD domain is outside MD overlap zone 
                if ( md_map%n .eq. 0 ) return 

                call mpi_comm_rank(COUPLER_COMM,myid,ierr)

                ncalls = ncalls + 1
                !                write(0,*) "CFD, send_CFDvel: ", myid
                is = bbox_cfd%xbb(1,icoord(1,myid+1))
                ie = bbox_cfd%xbb(2,icoord(1,myid+1))
                js = bbox_cfd%ybb(1,icoord(2,myid+1))
                je = bbox_cfd%ybb(2,icoord(2,myid+1))
                ks = bbox_cfd%zbb(1,icoord(3,myid+1))
                ke = bbox_cfd%zbb(2,icoord(3,myid+1))

                iu_s = 2
                ku_s = 1
                iv_s = 2
                kv_s = 1

                iu_e = iu_s + ie - is - 1 ! +1 - 1 !!
                iv_e = iv_s + ie - is - 1 
                ku_e = ku_s + ke - ks - 1
                kv_e = kv_s + ke - ks - 1 
                !                write(0,*) 'CFD send_CFDvel' , myid, ncalls, ie-is+1, nlx,nly,nlz
                min_i = minval(md_map%domains(1,:))
                min_j = minval(md_map%domains(3,:))
                min_k = minval(md_map%domains(5,:))

                do j = 1, 2
                        select case (j)
                        case (1)
                                vaux(:,:,:) = uc(ku_s:ku_e,iu_s:iu_e,2:jmax_overlap-1)
                        case (2)
                                vaux(:,:,:) = vc(kv_s:kv_e,iv_s:iv_e,2:jmax_overlap-1)
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


                ! This barrier does not work as this function si called inside an jblock .eq. 1
                ! condition. Is it really necessary?
                !                call mpi_barrier(CFD_COMM, ierr)

	end subroutine coupler_send_CFDvel


        subroutine coupler_md_vel(uc,vc,wc)
                use mpi
                use coupler_internal_common
                !                use data_export, only : ngz, i1_u, i2_u, j1_u, j2_u,&
                !                        i1_v, i2_v, j1_v, j2_v, i1_w, i2_w, j1_w, j2_w
                !                use messenger, only : myid, icoord
                use coupler_internal_cfd, only : nlx, nlz, recv_vel_MD, bbox_cfd, icoord
                implicit none

                real(kind=kind(0.d0)),dimension(:,:,:),intent(out) :: uc, vc, wc 

                integer is, ie, iu_s, iu_e, iv_s, iv_e, iw_s, iw_e, myid
                integer  ierr

                call mpi_comm_rank(COUPLER_COMM, myid, ierr)
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

                call  recv_vel_MD(uc, 1, nlz-1, 1, nlx-1, 1, 1, 0)
                call  recv_vel_MD(vc, 1, nlz-1, 1, nlx-1, 1, 1, 0)
		wc = 0.d0
                !                call  recv_vel_MD(wc, 1, ngz-1, i1_w, i2_w, 0, 0, 1)

        end subroutine coupler_md_vel

!-----------------------------------------------------------------------------
!
! function and subroutines that extract parameter forn internal modules 
!
!-----------------------------------------------------------------------------        
        subroutine coupler_get(xL_md,yL_md,zL_md)
                use coupler_internal_md, only : xL_md_ =>  xL_md, yL_md_ =>  yL_md, zL_md_ =>  zL_md
!                use coupler_internal_md
                implicit none

                real(kind(0.d0)), optional, intent(out) :: xL_md, yL_md, zL_md

                if (present(xL_md)) then 
                        xL_md = xL_md_
                endif

                if (present(yL_md)) then 
                        yL_md = yL_md_
                endif

                if (present(zL_md)) then 
                        zL_md = zL_md_
                endif
                
        end subroutine coupler_get


        function coupler_get_save_period() result(p)
                use coupler_internal_md, only : save_period
                implicit none

                integer p

                p = save_period
        end function coupler_get_save_period


        function coupler_get_average_period() result(p)
                use coupler_internal_md, only : average_period
                implicit none

                integer p

                p = average_period
        end function coupler_get_average_period


        function coupler_get_nsteps() result(n)
                 use coupler_internal_md, only : nsteps 

                 integer n

                 n = nsteps
         end function coupler_get_nsteps


        function coupler_md_get_dt_cfd() result(t)
                 use coupler_internal_md, only : dt_CFD  

                 real(kind=kind(0.d0)) t

                 t = dt_CFD
        end function coupler_md_get_dt_cfd
         

end module coupler
