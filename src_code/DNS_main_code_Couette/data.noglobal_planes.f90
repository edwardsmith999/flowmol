!=======================================================================
! Defining data index limits; handling data allocation, 
! data blocking, and data file input/output
!
! data_init()
! data_read()
! data_write()
! data_free()
! data_setLimits_1()*
! data_setLimits_2()*
! data_readData()
! data_writeData()
! data_writeZdata()
! data_isRoot(iflag)
! data_info(iunit)
!

module data
    use data_export
	use computation_parameters
	integer nfile			! file number for zdata files
end module

!=======================================================================
subroutine data_init()
	use data

	! Check dimensions in archive and export file
	! nx,ny,nz defined in export file
	call readInt("nx", ni)
	call readInt("ny", nj)
	call readInt("nz", nk)
	if (ni .ne. nx .or. nj .ne. ny .or. nk .ne. nz) &
		stop "Dimensions in archive and export file do not match"


	!  if (niz .ne. 1 .and. mod(niz,2) .ne. 0) &
	!  	stop "Incompatible dimension: niz must be even, unless niz=1"

	nfile = 0

	!-------------------------------------------------------------------
	! Initialize message passing
	!-------------------------------------------------------------------
	call messenger_init()

        ! Layout block ID
        iblock_1 = iblock
        jblock_1 = jblock
        kblock_1 = kblock
        iblock_2 = 1
        jblock_2 = 1
        kblock_2 = irank

	! Set block limits
	call data_setLimits_1()
	call data_setLimits_2()

	! Shorthand for primary layout

	nixb = nixb_1(iblock_1)
	niyb = niyb_1(jblock_1)
	nizb = nizb_1(kblock_1)

	nxb = nxb_1(iblock_1)
	nyb = nyb_1(jblock_1)
	nzb = nzb_1(kblock_1)

	ibmin = ibmin_1(iblock_1)
	jbmin = jbmin_1(jblock_1)
	kbmin = kbmin_1(kblock_1)

	ibmax = ibmax_1(iblock_1)
	jbmax = jbmax_1(jblock_1)
	kbmax = kbmax_1(kblock_1)

	ibmino = ibmino_1(iblock_1)
	jbmino = jbmino_1(jblock_1)
	kbmino = kbmino_1(kblock_1)

	ibmaxo = ibmaxo_1(iblock_1)
	jbmaxo = jbmaxo_1(jblock_1)
	kbmaxo = kbmaxo_1(kblock_1)

	ib1_ = imap_1(ibmin)
	jb1_ = jmap_1(jbmin)
	kb1_ = kmap_1(kbmin)
	ib2_ = imap_1(ibmax)
	jb2_ = jmap_1(jbmax)
	kb2_ = kmap_1(kbmax)


      !-------- Set Xhiaoua's local limits--------

      !---Primary Layout
        do ib=1,npx
          nlxb_(ib) = nxb_1(ib) - 2*nox
        end do

        do jb=1,npy
          nlyb_(jb) = nyb_1(jb) - 2*noy
        end do

        nlxb = nlxb_(iblock_1)
        nlyb = nlyb_(jblock_1)

      !---Secondary Layout (Note: 'noz2' because no halos in poisson)
        do kb=1,nproc
           nlzb_(kb) = nzb_2(kb) - 2*noz2
        end do
        nlzb = nlzb_(irank)

        !-----Advance{u} Limits----
         i1_u = 0 + nox1
         i2_u = nlxb - nox1
         if (iblock_1.eq.1)      i1_u = 2
         !if (iblock_1.eq.nbx_1)  i2_u = nlxb - 1	!TJ: For inflow-outflow code
         if (iblock_1.eq.nbx_1)  i2_u = nlxb 		!TJ: For x-periodic channel code
         niu = i2_u - i1_u + 1
         ! print*,'( RANK, i1_u,i2_u, niu ) = ', irank, i1_u,i2_u,niu

         j1_u = 0 + noy1
         j2_u   = nlyb - noy1
         if (jblock_1.eq.1)     j1_u = 1
         if (jblock_1.eq.nby_1) j2_u = nlyb -1
         nju = j2_u - j1_u + 1
         ! print*,'( RANK, j1_u,j2_u, nju ) = ', irank, j1_u,j2_u,nju

        !-----Advance{v} Limits----
         i1_v = 0 + nox1
         i2_v = nlxb - nox1
         if (iblock_1.eq.1)     i1_v = 1
         if (iblock_1.eq.nbx_1) i2_v = nlxb - 1
         niv = i2_v - i1_v + 1
         ! print*,'( RANK, i1_v,i2_v, niv ) = ', irank, i1_v,i2_v,niv

         j1_v = 0 + noy1
         j2_v = nlyb - noy1
         if (jblock_1.eq.1)     j1_v = 2
         if (jblock_1.eq.nby_1) j2_v = nlyb - 1
         njv = j2_v - j1_v + 1
         ! print*,'( RANK, j1_v,j2_v, njv ) = ', irank, j1_v,j2_v,njv

        !-----Advance{w} Limits----
         i1_w = 0 + nox1
         i2_w = nlxb - nox1
         if (iblock_1.eq.1)     i1_w = 1
         if (iblock_1.eq.nbx_1)  i2_w = nlxb - 1
         niw = i2_w - i1_w + 1
         ! print*,'( RANK, i1_w,i2_w, niw ) = ', irank, i1_w,i2_w,niw

         j1_w = 0 + noy1
         j2_w = nlyb - noy1
         if (jblock_1.eq.1)     j1_w = 1
         if (jblock_1.eq.nby_1) j2_w   = nlyb - 1
         njw = j2_w - j1_w + 1
         ! print*,'( RANK, j1_w,j2_w, njw ) = ', irank, j1_w,j2_w,njw

        !-----Poisson/fft qr Limits----
        !-- x-direction
         i1_T = nox1
         i2_T = nlxb-nox1
         if (iblock_1.eq.1  ) i1_T = 1
         if (iblock_1.eq.npx) i2_T = nlxb - 1

         iTmin_1(1) = 1
         iTmax_1(1) = iTmin_1(1) + (nlxb_(1)-nox1) -1
!#if nbx_1 > 1
         if (nbx_1.ne.1) then
            do ib=2,nbx_1-1
                    iTmin_1(ib) = iTmax_1(ib-1) + 1
                    iTmax_1(ib) = iTmin_1(ib) + (nlxb_(ib) - 2*nox1 +1) -1
            end do
            ib = nbx_1-1
            iTmin_1(nbx_1) = iTmax_1(ib) + 1
            iTmax_1(nbx_1) = iTmin_1(nbx_1) + (nlxb_(nbx_1)-1 - nox1 +1) -1
!#else
         else
            iTmax_1(1) = nlxb_(1) -1
         end if
!#endif

        !-- y-direction
         j1_T = noy1
         j2_T = nlyb-noy1
         if (jblock_1.eq.1  ) j1_T = 1
         if (jblock_1.eq.npy) j2_T = nlyb - 1

         jTmin_1(1) = 1
         jTmax_1(1) = jTmin_1(1) + (nlyb_(1) - noy1) -1
!#if nby_1 > 1
         if (nby_1.ne.1) then
            do jb=2,nby_1-1
                    jTmin_1(jb) = jTmax_1(jb-1) + 1
                    jTmax_1(jb) = jTmin_1(jb) + (nlyb_(jb) - 2*noy1 +1) -1
            end do
            jb = nby_1-1
            jTmin_1(nby_1) = jTmax_1(jb) + 1 ! avoids overzealos compiler
            jTmax_1(nby_1) = jTmin_1(nby_1) + (nlyb_(nby_1)-1 - noy1 +1) -1
!#else
         else
            jTmax_1(1) = nlyb_(1) - 1
         end if
!#endif

        !-- z-direction (transposed)
         k1_T = 1
         k2_T = nlzb-1

         kTmin_2(1) = 1
         kTmax_2(1) = kTmin_2(1) + (nlzb_(1)-1) -1
         do kb=2,nbz_2
                 kTmin_2(kb) = kTmax_2(kb-1) + 1
                 kTmax_2(kb) = kTmin_2(kb) + (nlzb_(kb)-1) - 1
         end do

	!--------------------------------
        !	 Print out some limits
	!--------------------------------
        !  if (irank.eq.1) then
        !     print*,'iTmin_1 = ', iTmin_1
        !     print*,'iTmax_1 = ', iTmax_1
	! 
        !     print*,'jTmin_1 = ', jTmin_1
        !     print*,'jTmax_1 = ', jTmax_1
	! 
        !     print*,'kTmin_2 = ', kTmin_2
        !     print*,'kTmax_2 = ', kTmax_2
	!  end if
	!--------------------------------

	!--------------------------------
        !  print*,'irank start'
        !  print*,'irank',irank
        !  print*,'nlyb',nlyb
        !  print*,'nyb_1',nyb_1
        !  print*,'niyb',niyb
        !  print*,'ny_1',ny_1
        !  print*,'niy_1',niy_1
        !  print*,'irank end '
	!--------------------------------

	return
end

subroutine data_read()
	use data
	call archive_isDefined("nfile", iflag)
	if (iflag == 1) call readInt("nfile", nfile)
	return
end

subroutine data_write()
  	use data
  	call writeInt("nx", nx)
  	call writeInt("ny", ny)
  	call writeInt("nz", nz)
	call writeInt("npx", npx)
	call writeInt("npy", npy)
	if (nfile .ne. 0) call writeInt("nfile", nfile)
	return
end

subroutine data_free()
	use data

	! Finialize message passing
	call messenger_free()

	return
end

!=======================================================================
subroutine data_setLimits_1()
	use data

!NEW WAY FOR MULTI-GRID BUT NOT WORKING ALL THE TIME
        nixb_1(:) = (nx-2*nox-1)/nbx_1 + 1
        niyb_1(:) = (ny-2*noy-1)/nby_1 + 1
        nizb_1(:) = (nz-2*noz-1)/nbz_1 + 1
        nrx = mod(nx-2*nox-1, nbx_1)
        nry = mod(ny-2*noy-1, nby_1)
        nrz = mod(nz-2*noz-1, nbz_1)

        nixb_1(1:nrx) = nixb_1(1:nrx) + 1
        niyb_1(1:nry) = niyb_1(1:nry) + 1
        nizb_1(1:nrz) = nizb_1(1:nrz) + 1

        nixb_1(1) = nixb_1(1)-(nox1-nox)
        niyb_1(1) = niyb_1(1)-(noy1-noy)

        nixb_1(nbx_1) = nixb_1(nbx_1)-(nox1-nox)
        niyb_1(nby_1) = niyb_1(nby_1)-(noy1-noy)

        nxb_1(:) = nixb_1(:) + 2*nox1
        nyb_1(:) = niyb_1(:) + 2*noy1
        nzb_1(:) = nizb_1(:) + 2*noz1
!END OF NEW WAY

!OLD WAY (WORKS)
!        nixb_1(:) = (nx-2*nox1-1)/nbx_1 + 1
!        niyb_1(:) = (ny-2*noy1-1)/nby_1 + 1
!        nizb_1(:) = (nz-2*noz1-1)/nbz_1 + 1
!        nrx = mod(nx-2*nox1-1, nbx_1)
!        nry = mod(ny-2*noy1-1, nby_1)
!        nrz = mod(nz-2*noz1-1, nbz_1)
!
!        nixb_1(1:nrx) = nixb_1(1:nrx) + 1
!        niyb_1(1:nry) = niyb_1(1:nry) + 1
!        nizb_1(1:nrz) = nizb_1(1:nrz) + 1
!
!        nxb_1(:) = nixb_1(:) + 2*nox1
!        nyb_1(:) = niyb_1(:) + 2*noy1
!        nzb_1(:) = nizb_1(:) + 2*noz1
!END OF OLD WAY

        ! Note the zero because I start from zero at the halo cells.
        ! ibmin/ibmax point to u-velocity location (needed for imap/ibmap,
        ! which we use all the time with (suxix, suxiy) at the u-location.

        ibmin_1(1) = 0 + nox1
        ibmax_1(1) = ibmin_1(1) + nixb_1(1) - 1
        do ib=2,nbx_1
                ibmin_1(ib) = ibmax_1(ib-1)
                ibmax_1(ib) = ibmin_1(ib) + nixb_1(ib) - 1
        end do

        jbmin_1(1) = 0 + noy1
        jbmax_1(1) = jbmin_1(1) + niyb_1(1) - 1
        do jb=2,nby_1
                jbmin_1(jb) = jbmax_1(jb-1)
                jbmax_1(jb) = jbmin_1(jb) + niyb_1(jb) - 1
        end do

        kbmin_1(1) = 0 + noz1
        kbmax_1(1) = kbmin_1(1) + nizb_1(1) - 1
        do kb=2,nbz_1
                kbmin_1(kb) = kbmax_1(kb-1)
                kbmax_1(kb) = kbmin_1(kb) + nizb_1(kb) - 1
        end do

        ibmino_1(:) = ibmin_1(:) - nox1
        jbmino_1(:) = jbmin_1(:) - noy1
        kbmino_1(:) = kbmin_1(:) - noz1

        ibmaxo_1(:) = ibmax_1(:) + nox1
        jbmaxo_1(:) = jbmax_1(:) + noy1
        kbmaxo_1(:) = kbmax_1(:) + noz1

        ! Index mapping
        ! Out-of-bounds maps to pad area
        imap_1 = -1
        jmap_1 = -1
        kmap_1 = -1

        ib = iblock_1
        do i=ibmino_1(ib),ibmaxo_1(ib)
                imap_1(i) = i-ibmino_1(ib)+0
        end do
        imap_1(i:nx-1) = nx

        do i=0,nx_1
                ibmap_1(i) = i+ibmino_1(ib)-0
        end do

        jb = jblock_1
        do j=jbmino_1(jb),jbmaxo_1(jb)
                jmap_1(j) = j-jbmino_1(jb)+0
        end do
        jmap_1(j:ny-1) = ny

        do j=0,ny_1
                jbmap_1(j) = j+jbmino_1(jb)-0
        end do

        kb = kblock_1
        do k=kbmino_1(kb),kbmaxo_1(kb)
                kmap_1(k) = k-kbmino_1(kb)+0
        end do
        kmap_1(k:nz-1) = nz

        do k=0,nz_1
                kbmap_1(k) = k+kbmino_1(kb)-0
        end do

	return
end

subroutine data_setLimits_2()
	use data

        nixb_2(:) = nixm / nbx_2 + 1
        niyb_2(:) = niym / nby_2 + 1
        !-- This (2*(ngzm/2)) is used to make sure we have
        !   an even set of wavenumbers
        nizb_2(:) = (2*(nizm/2)) / nbz_2 + 1

        nrx = mod(nixm, nbx_2)
        nry = mod(niym, nby_2)
        nrz = mod((2*(nizm/2)), nbz_2)
        nixb_2(1:nrx) = nixb_2(1:nrx) + 1
        niyb_2(1:nry) = niyb_2(1:nry) + 1
        nizb_2(1:nrz) = nizb_2(1:nrz) + 1

        nxb_2(:) = nixb_2(:) + 2*nox2
        nyb_2(:) = niyb_2(:) + 2*noy2
        nzb_2(:) = nizb_2(:) + 2*noz2

        ! Note the zero because I start from zero at the halo cells.
        ibmin_2(1) = 0 + nox2
        ibmax_2(1) = ibmin_2(1) + nixb_2(1) -1 
        do ib=2,nbx_2
                ibmin_2(ib) = ibmax_2(ib-1)
                ibmax_2(ib) = ibmin_2(ib) + nixb_2(ib) - 1
        end do

        jbmin_2(1) = 0 + noy2
        jbmax_2(1) = jbmin_2(1) + niyb_2(1) -1
        do jb=2,nby_2
                jbmin_2(jb) = jbmax_2(jb-1)
                jbmax_2(jb) = jbmin_2(jb) + niyb_2(jb) - 1
        end do

        ! Note the ONE because I start from ONE in poisson solver.
        kbmin_2(1) = 1 + noz2
        kbmax_2(1) = kbmin_2(1) + nizb_2(1) -1
        do kb=2,nbz_2
                kbmin_2(kb) = kbmax_2(kb-1)
                kbmax_2(kb) = kbmin_2(kb) + nizb_2(kb) - 1
        end do

        ibmino_2(:) = ibmin_2(:) - nox2
        jbmino_2(:) = jbmin_2(:) - noy2
        kbmino_2(:) = kbmin_2(:) - noz2

        ibmaxo_2(:) = ibmax_2(:) + nox2
        jbmaxo_2(:) = jbmax_2(:) + noy2
        kbmaxo_2(:) = kbmax_2(:) + noz2

        ! Index mapping
        ! Out-of-bounds maps to pad area
        imap_2 = nx_2
        jmap_2 = ny_2
        kmap_2 = nz_2

        ib = iblock_2
        do i=ibmino_2(ib),ibmaxo_2(ib)
                imap_2(i) = i-ibmino_2(ib)+0
        end do
        do i=0,nxm_2
                ibmap_2(i) = i+ibmino_2(ib)-0
        end do

        jb = jblock_2
        do j=jbmino_2(jb),jbmaxo_2(jb)
                jmap_2(j) = j-jbmino_2(jb)+0
        end do
        do j=0,nym_2
                jbmap_2(j) = j+jbmino_2(jb)-0
        end do

        ! Note the ONE because I start from ONE in poisson solver.
        kb = kblock_2
        do k=kbmino_2(kb),kbmaxo_2(kb)
                kmap_2(k) = k-kbmino_2(kb)+1
        end do
        do k=1,nzm_2
                kbmap_2(k) = k+kbmino_2(kb)-1
        end do

         return
end


!=======================================================================
subroutine data_readData()
	use data

        ! Use a record length of nz
        inquire(iolength=ilength) U(1,1,:)
        iunit = iopen()
        open (iunit, file=trim(prefix_dir)//"data", form="unformatted", &
                access="direct",recl=ilength, &
                status="old", iostat=ierr)
        if (ierr .ne. 0) stop "A data file is required"

        read (iunit, rec=1) ni, nj, nk, nvar
        ndata = 5
        if (nvar .ne. ndata) stop "Wrong number of variables in data file"
        if (ni .ne. nx .or. nj .ne. ny .or. nk .ne. nz) &
                stop "Wrong dimensions in data file"

        ! Global indices
        i1 = ibmino  ! ; if (ibmin == imin) i1 = imino
        i2 = ibmaxo  ! ; if (ibmax == imax) i2 = imaxo
        j1 = jbmino  ! ; if (jbmin == jmin) j1 = jmino
        j2 = jbmaxo  ! ; if (jbmax == jmax) j2 = jmaxo

        n = 1
        do j=j1,j2   
        do i=i1,i2
                read (iunit,rec=(n-1)*nx*ny+(j-1)*nx+i+1) &
                        (U(imap_1(i),jmap_1(j),k), k=1,nz_)
        enddo
        enddo   
        
        n = 2
        do j=j1,j2   
        do i=i1,i2
                read (iunit,rec=(n-1)*nx*ny+(j-1)*nx+i+1) &
                        (V(imap_1(i),jmap_1(j),k), k=1,nz_)
        enddo
        enddo   
        
        n = 3
        do j=j1,j2    
        do i=i1,i2
                read (iunit,rec=(n-1)*nx*ny+(j-1)*nx+i+1) &
                        (W(imap_1(i),jmap_1(j),k), k=1,nz_)
        enddo
        enddo   
        
        n = 4
        do j=j1,j2
        do i=i1,i2
                read (iunit,rec=(n-1)*nx*ny+(j-1)*nx+i+1) &
                        (P(imap_1(i),jmap_1(j),k), k=1,nz_)
        enddo
        enddo   

        n = 5 
        do j=j1,j2
        do i=i1,i2
        !T        read (iunit,rec=(n-1)*nx*ny+(j-1)*nx+i+1) &
        !T                (VT(imap_1(i),jmap_1(j),k), k=1,nz_)
        enddo   
        enddo

        close (iclose(iunit))
        
      if(stime.ne.0.0) then
      else
      end if

        return
end     


subroutine data_writeData()
	use data

        ! This routine writes simultaneously from all MPI processes to
        ! a single output file with no record written by more than one
        ! process; this works on the T3E with "assign -m on data" set;
        ! does not currently work on Origin 2000 system

        ! Use a record length of nz
        inquire(iolength=ilength) U(1,1,:)

        iunit = iopen()
        open (iunit, file=trim(prefix_dir)//"data", form="unformatted", &
                access="direct",recl=ilength,iostat=ierr)

        ni = nx
        nj = ny
        nk = nz
        nvar = 5
        ! Header is first record
        ! Record length should always exceed the length of the header
        ! On Origin2000 reclReal = 8, reclInt = 4
        ! Since nz >= 3 --> can fit at least 6 integers
        
        if (irank == iroot) write (iunit,rec=1) ni, nj, nk, nvar
        
        ! Global indices
        i1 = ibmin  ; if (ibmin == imin) i1 = imino
        i2 = ibmax  ; if (ibmax == imax) i2 = imaxo
        j1 = jbmin  ; if (jbmin == jmin) j1 = jmino
        j2 = jbmax  ; if (jbmax == jmax) j2 = jmaxo
 
        n = 1
        do j=j1,j2
        do i=i1,i2
                write (iunit,rec=(n-1)*nx*ny+(j-1)*nx+i+1) &
                        (U(imap_1(i),jmap_1(j),k), k=1,nz_)
        enddo   
        enddo
        
        n = 2
        do j=j1,j2
        do i=i1,i2
                write (iunit,rec=(n-1)*nx*ny+(j-1)*nx+i+1) &
                        (V(imap_1(i),jmap_1(j),k), k=1,nz_)
        enddo
        enddo


        n = 3
        do j=j1,j2
        do i=i1,i2
                write (iunit,rec=(n-1)*nx*ny+(j-1)*nx+i+1) &
                        (W(imap_1(i),jmap_1(j),k), k=1,nz_)
        enddo   
        enddo

        n = 4
        do j=j1,j2 
        do i=i1,i2 
                write (iunit,rec=(n-1)*nx*ny+(j-1)*nx+i+1) &
                        (P(imap_1(i),jmap_1(j),k), k=1,nz_)
        enddo
        enddo

        n = 5
        do j=j1,j2
        do i=i1,i2
        !T        write (iunit,rec=(n-1)*nx*ny+(j-1)*nx+i+1) &
        !T                (VT(imap_1(i),jmap_1(j),k), k=1,nz_)
        enddo  
        enddo   

!        call flush(iunit,ierr) 
        call flush(iunit)
        
        close (iclose(iunit))
        
        return
end     
        

subroutine data_writeYdata()
        use data
        character(LEN=10) :: fname
        data fname / "zdata.xxxx" /

        jloc = 13       ! y-plane to save
        if (jloc > ny) jloc = ny/2
        write(fname(7:10),"(i4.4)") nfile
        nfile = nfile + 1
        iunit = iopen()
        open (iunit, file=trim(prefix_dir)//fname, form="unformatted")
        write (iunit) nx, 1, nz, 4, jloc
        write (iunit) U(:,jloc,:) 
        write (iunit) V(:,jloc,:)
        write (iunit) W(:,jloc,:)
        write (iunit) P(:,jloc,:)
        close (iclose(iunit))

        return  
end  

subroutine data_writeZdata()
	use data
	character(LEN=10) :: fname
	data fname / "zdata.xxxx" /

!!	write(fname(7:10),"(i4.4)") nfile
!!	nfile = nfile + 1
!!	iunit = iopen()
!!	open (iunit, file=fname, form="formatted")
!!	write (iunit) nx, ny, 1, 5
!!	write (iunit) U(:,:,nz/2)
!!	write (iunit) V(:,:,nz/2)
!!	write (iunit) W(:,:,nz/2)
!!	write (iunit) P(:,:,nz/2)
!!	write (iunit) VT(:,:,nz/2)

!!            call gather(U, R)
!!            if (irank == iroot) then
!!            write(21,9991) stime, ((R(i,j,1),i=1,imax),j=2,2)
!!            write(21,9991) stime, ((R(i,j,1),i=1,imax),j=jmax/2,jmax/2)
!!            else
!!            end if
!!            call gather(P, R)
!!            if (irank == iroot) then
!!            write(21,9991) stime, ((R(i,j,1),i=1,imax),j=1,1)
!!            write(21,9991) stime, ((R(i,j,1),i=1,imax),j=jmax/2,jmax/2) 
!!           else
!!            end if
!! 
 9991       format(1x,6(1pe13.5))

!	close (iclose(iunit))

	return
end

!=======================================================================
subroutine data_isRoot(iflag)
	use data

	if (irank == iroot) then
		iflag = 1
	else
		iflag = 0
	end if

	return
end

!=======================================================================
subroutine data_info(iunit)
	use data

	return
end

