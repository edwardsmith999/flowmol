subroutine restrivp_l(pc_mt, pf_mt, ncx, ncy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! restriction :  determine coarse grid volume values using
!-----------------------------------------------------------------------
!                fine grid volume values
!                     pf = finer pressure cell
!                     pc = coarser pressure cell
!-----------------------------------------------------------------------
! 		 This subroutine is for local domains.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	use data_export
	implicit real (a-h,o-z)
	dimension pc_mt(0:ncx+1,0:ncy+1), pf_mt(0:2*ncx+1,0:2*ncy+1)

	do jc=1,ncy
	jf=2*jc-1
	do ic=1,ncx
		if=2*ic-1
		pc_mt(ic,jc)=pf_mt(if,jf)+pf_mt(if+1,jf)+pf_mt(if,jf+1)+pf_mt(if+1,jf+1)
	end do
	end do

	!cccccccccccccccccccccccccc
	! update border
	!cccccccccccccccccccccccccc
	call updateBorder_2D(pc_mt, ncx, ncy, id_x, 1, 1, ncx, 1)
	call updateBorder_2D(pc_mt, ncx, ncy, id_y, 2, 1, ncy, 1)

	!cccccccccccccccccccccccccc
	! boundary points
	!cccccccccccccccccccccccccc
	if(jblock==1  ) pc_mt(1:ncx,0    )=pc_mt(1:ncx,1  )
	if(jblock==npy) pc_mt(1:ncx,ncy+1)=pc_mt(1:ncx,ncy)

	if(iblock==1  ) pc_mt(0    ,:)=pc_mt(1  ,:)
	if(iblock==npx) pc_mt(ncx+1,:)=pc_mt(ncx,:)
	return
end

subroutine restrisu_l(suc_mt, suf_mt, ncx, ncy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! restriction :  determine coarse grid area vector values using
!-----------------------------------------------------------------------
!                fine grid area vector values
!                     pf = finer pressure cell
!                     pc = coarser pressure cell
!-----------------------------------------------------------------------
! 		 This subroutine is for local domains.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	use data_export
	implicit real (a-h,o-z)
	dimension suc_mt(0:ncx+1+1,0:ncy+1),suf_mt(0:2*ncx+1+1,0:2*ncy+1)

	do jc=1,ncy
	jf=2*jc-1
	do ic=1,ncx+1
		if=2*ic-1
		suc_mt(ic,jc)=suf_mt(if,jf)+suf_mt(if,jf+1)
	end do
	end do

	!cccccccccccccccccccccccccc
	! update border
	!cccccccccccccccccccccccccc
	!call updateBorder_2D_suv(suc_mt, ncx+1, ncy, id_x, 1, 2, ncx    )
	!----------------------------------------
	! If I don't use the above subroutine,
	! I need to call updateBorder_2D twice.
	!----------------------------------------
	call updateBorder_2D    (suc_mt, ncx+1, ncy, id_x, 1, 1, ncx  ,1)
	call updateBorder_2D    (suc_mt, ncx+1, ncy, id_x, 1, 2, ncx+1,1)
	call updateBorder_2D    (suc_mt, ncx+1, ncy, id_y, 2, 1, ncy  ,1)

	!cccccccccccccccccccccccccc
	! lower and upper b.c.
	!cccccccccccccccccccccccccc
	if(jblock==1  ) suc_mt(1:ncx+1,0    )=suc_mt(1:ncx+1,1  )
	if(jblock==npy) suc_mt(1:ncx+1,ncy+1)=suc_mt(1:ncx+1,ncy)
	!cccccccccccccccccccccccccc
	! left and right b.c.
	!cccccccccccccccccccccccccc
	if(iblock==1  ) suc_mt(0    ,:)=suc_mt(2  ,:)
	if(iblock==npx) suc_mt(ncx+2,:)=suc_mt(ncx,:)
	return
end

subroutine restrisv_l(svc_mt, svf_mt, ncx, ncy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! restriction :  determine coarse grid area vector values using
!-----------------------------------------------------------------------
!                fine grid area vector values
!                     pf = finer pressure cell
!                     pc = coarser pressure cell
!-----------------------------------------------------------------------
! 		 This subroutine is for local domains.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	use data_export
	implicit real (a-h,o-z)
	dimension svc_mt(0:ncx+1,0:ncy+1+1), svf_mt(0:2*ncx+1,0:2*ncy+1+1)

	do jc=1,ncy+1
	jf=2*jc-1
	do ic=1,ncx
		if=2*ic-1
		svc_mt(ic,jc)=svf_mt(if,jf)+svf_mt(if+1,jf)
	end do
	end do
	!cccccccccccccccccccccccccc
	! update border
	!cccccccccccccccccccccccccc
	call updateBorder_2D    (svc_mt, ncx, ncy+1, id_x, 1, 1, ncx  , 1)
	!call updateBorder_2D_suv(svc_mt, ncx, ncy+1, id_y, 2, 2, ncy     )
	!----------------------------------------
	! If I don't use the above subroutine,
	! I need to call updateBorder_2D twice.
	!----------------------------------------
	call updateBorder_2D    (svc_mt, ncx, ncy+1, id_y, 2, 1, ncy  , 1)
	call updateBorder_2D    (svc_mt, ncx, ncy+1, id_y, 2, 2, ncy+1, 1)

	!cccccccccccccccccccccccccc
	! left and right b.c.
	!cccccccccccccccccccccccccc
	if(iblock==1  ) svc_mt(0    ,1:ncy+1)=svc_mt(1  ,1:ncy+1)
	if(iblock==npx) svc_mt(ncx+1,1:ncy+1)=svc_mt(ncx,1:ncy+1)
	!cccccccccccccccccccccccccc
	! upper and lower b.c.
	!cccccccccccccccccccccccccc
	if(jblock==1  ) svc_mt(:,0    )=svc_mt(:,2  )
	if(jblock==npy) svc_mt(:,ncy+2)=svc_mt(:,ncy)
	return
end

subroutine restrisw_l(swc_mt, swf_mt, ncx, ncy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! restriction :  determine coarse grid area vector values using
!-----------------------------------------------------------------------
!                fine grid area vector values
!                     pf = finer pressure cell
!                     pc = coarser pressure cell
!-----------------------------------------------------------------------
! 		 This subroutine is for local domains.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	use data_export
	implicit real (a-h,o-z)
	dimension swc_mt(0:ncx+1,0:ncy+1), swf_mt(0:2*ncx+1,0:2*ncy+1)

	do jc=1,ncy
	jf=2*jc-1
	do ic=1,ncx
		if=2*ic-1
		swc_mt(ic,jc)=swf_mt(if,jf)+swf_mt(if+1,jf)+swf_mt(if,jf+1)+swf_mt(if+1,jf+1)
	end do
	end do
	!cccccccccccccccccccccccccc
	! update border
	!cccccccccccccccccccccccccc
	call updateBorder_2D( swc_mt, ncx, ncy, id_x, 1, 1, ncx, 1)
	call updateBorder_2D( swc_mt, ncx, ncy, id_y, 2, 1, ncy, 1)

	!cccccccccccccccccccccccccc
	! left and right b.c.
	!cccccccccccccccccccccccccc
	if(iblock==1  ) swc_mt(0    ,1:ncy)=swc_mt(1  ,1:ncy)
	if(iblock==npx) swc_mt(ncx+1,1:ncy)=swc_mt(ncx,1:ncy)
	!cccccccccccccccccccccccccc
	! upper and lower b.c.
	!cccccccccccccccccccccccccc
	if(jblock==1  ) swc_mt(:,0    )=swc_mt(:,  1)
	if(jblock==npy) swc_mt(:,ncy+1)=swc_mt(:,ncy)
	return
end

subroutine restrivp(pc_mt, pf_mt, ncx, ncy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! restriction :  determine coarse grid volume values using
!-----------------------------------------------------------------------
!                fine grid volume values
!                     pf = finer pressure cell
!                     pc = coarser pressure cell
!-----------------------------------------------------------------------
! 		 This subroutine is for global domains.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit real (a-h,o-z)
	dimension pc_mt(0:ncx+1,0:ncy+1), pf_mt(0:2*ncx+1,0:2*ncy+1)

	do jc=1,ncy
	jf=2*jc-1
	do ic=1,ncx
		if=2*ic-1
		pc_mt(ic,jc)=pf_mt(if,jf)+pf_mt(if+1,jf)+pf_mt(if,jf+1)+pf_mt(if+1,jf+1)
	end do
	end do
	!cccccccccccccccccccccccccc
	! boundary points
	!cccccccccccccccccccccccccc
	pc_mt(1:ncx,0    )=pc_mt(1:ncx,1  )
	pc_mt(1:ncx,ncy+1)=pc_mt(1:ncx,ncy)

	pc_mt(0    ,:)=pc_mt(1  ,:)
	pc_mt(ncx+1,:)=pc_mt(ncx,:)
	return
end


subroutine restrisu(suc_mt,suf_mt,ncx,ncy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! restriction :  determine coarse grid area vector values using
!-----------------------------------------------------------------------
!                fine grid area vector values
!                     pf = finer pressure cell
!                     pc = coarser pressure cell
!-----------------------------------------------------------------------
! 		 This subroutine is for global domains.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit real (a-h,o-z)
	dimension suc_mt(0:ncx+1+1,0:ncy+1), suf_mt(0:2*ncx+1+1,0:2*ncy+1)

	do jc=1,ncy
	jf=2*jc-1
	do ic=1,ncx+1
		if=2*ic-1
		suc_mt(ic,jc)=suf_mt(if,jf)+suf_mt(if,jf+1)
	end do
	end do
	!cccccccccccccccccccccccccc
	! lower and upper b.c.
	!cccccccccccccccccccccccccc
	suc_mt(1:ncx+1,0    )=suc_mt(1:ncx+1,1)
	suc_mt(1:ncx+1,ncy+1)=suc_mt(1:ncx+1,ncy)
	!cccccccccccccccccccccccccc
	! left and right b.c.
	!cccccccccccccccccccccccccc
	suc_mt(0    ,:)=suc_mt(2  ,:)
	suc_mt(ncx+2,:)=suc_mt(ncx,:)
	return
end


subroutine restrisv(svc_mt,svf_mt,ncx,ncy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! restriction :  determine coarse grid area vector values using
!-----------------------------------------------------------------------
!                fine grid area vector values
!                     pf = finer pressure cell
!                     pc = coarser pressure cell
!-----------------------------------------------------------------------
! 		 This subroutine is for global domains.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit real (a-h,o-z)
	dimension svc_mt(0:ncx+1,0:ncy+1+1), svf_mt(0:2*ncx+1,0:2*ncy+1+1)

	do jc=1,ncy+1
	jf=2*jc-1
	do ic=1,ncx
		if=2*ic-1
		svc_mt(ic,jc)=svf_mt(if,jf)+svf_mt(if+1,jf)
	end do
	end do
	!cccccccccccccccccccccccccc
	! left and right b.c.
	!cccccccccccccccccccccccccc
	svc_mt(0    ,1:ncy+1)=svc_mt(1  ,1:ncy+1)
	svc_mt(ncx+1,1:ncy+1)=svc_mt(ncx,1:ncy+1)
	!cccccccccccccccccccccccccc
	! upper and lower b.c.
	!cccccccccccccccccccccccccc
	svc_mt(:,0    )=svc_mt(:,2  )
	svc_mt(:,ncy+2)=svc_mt(:,ncy)
	return
end



subroutine preprocess_mt
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! allocate memory for Poisson coefficients at all levels
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	use data_export

	!TAZ  if (irank.eq.1)   &
		write(*,*) 'WARNING: mt int arrays: ngrid  ,nglevel_max  =',ngrid,nglevel_max

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! assign (restrict) area vector and volume at all levels
!-----------------------------------------------------------------------
	!ccccccccccccccccccccccccccccccccccccccccc
	! assign at the finest level
	!ccccccccccccccccccccccccccccccccccccccccc
	m = ngrid
	nxg = nxgm
	nyg = nygm
	nxl = mgl(m)%nx
	nyl = mgl(m)%ny
        ia = (iblock-1)*nxl + 1
        ib = ia + nxl - 1
        ja = (jblock-1)*nyl + 1
        jb = ja + nyl - 1
        mgl(m)%zvp (:,:)=vp    (ia-1:ib+1,ja-1:jb+1)
        mgl(m)%zsux(:,:)=suxix (ia-1:ib+2,ja-1:jb+1)
        mgl(m)%zsuy(:,:)=suxiy (ia-1:ib+2,ja-1:jb+1)
        mgl(m)%zsvx(:,:)=svetax(ia-1:ib+1,ja-1:jb+2)
        mgl(m)%zsvy(:,:)=svetay(ia-1:ib+1,ja-1:jb+2)

	!ccccccccccccccccccccccccccccccccccccccccc
	! restrict at lower levels
	!ccccccccccccccccccccccccccccccccccccccccc
	do m=ngrid-1,1,-1
	nxl = mgl(m)%nx
	nyl = mgl(m)%ny
	call restrivp_l(mgl(m)%zvp , mgl(m+1)%zvp , nxl, nyl)
	call restrisu_l(mgl(m)%zsux, mgl(m+1)%zsux, nxl, nyl)
	call restrisu_l(mgl(m)%zsuy, mgl(m+1)%zsuy, nxl, nyl)
	call restrisv_l(mgl(m)%zsvx, mgl(m+1)%zsvx, nxl, nyl)
	call restrisv_l(mgl(m)%zsvy, mgl(m+1)%zsvy, nxl, nyl)
	end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! evaluate Poisson coefficients at all levels
!---------------------------------------------------
! from fine to coarse
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	do m=ngrid,1,-1
	nxl = mgl(m)%nx
	nyl = mgl(m)%ny
		call pcoef_mgl( mgl(m)%zapw , mgl(m)%zape , &
				mgl(m)%zaps , mgl(m)%zapn , &
		 		mgl(m)%zapws, mgl(m)%zapwn, &
				mgl(m)%zapes, mgl(m)%zapen, &
				mgl(m)%zapp , mgl(m)%zvp  , &
				mgl(m)%zsux , mgl(m)%zsuy , &
				mgl(m)%zsvx , mgl(m)%zsvy , &
				mgl(m)%zedge, nxl, nyl)

	end do

	return
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!SY local only & 2D multigrid
subroutine pcoef_mgl(	apw_mt   , ape_mt   , aps_mt  , apn_mt  , &
			apws_mt  , apwn_mt  , apes_mt , apen_mt , &
			app_mt   , vp_mt    , suxix_mt, suxiy_mt, &
			svetax_mt, svetay_mt, aedgex  , nnx, nny)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! compute the Poisson equation coefficients at one
	! multigrid level with (nnx, nny) points
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	use data_export
	real	apw_mt (nnx,nny), ape_mt (nnx,nny), &
	    	aps_mt (nnx,nny), apn_mt (nnx,nny), &
		apws_mt(nnx,nny), apwn_mt(nnx,nny), &
		apes_mt(nnx,nny), apen_mt(nnx,nny), app_mt (nnx,nny), &
		suxix_mt (0:nnx+1+1,0:nny+1  ), suxiy_mt (0:nnx+1+1,0:nny+1  ), &
		svetax_mt(0:nnx+1  ,0:nny+1+1), svetay_mt(0:nnx+1  ,0:nny+1+1), &
	   	vp_mt    (0:nnx+1  ,0:nny+1  ), aedgex(2)

	integer	iedgex(2)

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! find the leading edge and trailing edge i index
	! pressure side and suction side have same index
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !-- Periodicity B.C.
        iedgex(1) = 1
        iedgex(2) = ngxm

        ia_block = iedgex(1)/nnx + 1
        ib_block = iedgex(2)/nnx + 1

        aedgex(1) = real(nnx+1)
        aedgex(2) = 0.0

        ! if (wall_flag.eq.1) then
        if (jblock.eq.1) then
         iedgex(1)  = max(  1 , iedgex(1)-(iblock-1)*nnx )
         iedgex(2)  = min( nnx, iedgex(2)-(iblock-1)*nnx )

         aedgex(1) = float(iedgex(1))
         aedgex(2) = float(iedgex(2))
        end if

	! if(iprob.eq.4) then
	! 	aedgex(1)=iedgex(1)
	! 	aedgex(2)=iedgex(2)
	! 	!TAZ  if (irank.eq.1) then
	! 		write(6,*) 'leading  edge located at i=    ',iedgex(1)
	! 		write(6,*) 'trailing edge located at i=    ',iedgex(2)
	! 	!TAZ  end if
	! end if

	!=================================================
	!     Initialize all the Poisson Coefficients
	!=================================================
	apws_mt = 0.0
	apw_mt  = 0.0
	apwn_mt = 0.0
	aps_mt  = 0.0
	app_mt  = 0.0
	apn_mt  = 0.0
	apes_mt = 0.0
	ape_mt  = 0.0
	apen_mt = 0.0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute a_ij associated with R_xi,e          (Page 37, 39-->40)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 31 j=1,nny
      do 31 i=1,nnx
!-----
        if(iblock.eq.npx.and.i.eq.nnx) then
          apws_mt(i,j)=apws_mt(i,j)
          apw_mt (i,j)=apw_mt (i,j)
          apwn_mt(i,j)=apwn_mt(i,j)

          aps_mt (i,j)=aps_mt (i,j)
          app_mt (i,j)=app_mt (i,j)
          apn_mt (i,j)=apn_mt (i,j)

          apes_mt(i,j)=apes_mt(i,j)
          ape_mt (i,j)=ape_mt (i,j)
          apen_mt(i,j)=apen_mt(i,j)

!-----
        else
!-----
          if((jblock.eq.npy.and.j.eq.nny).or. &
             (jblock.eq.1  .and.j.eq.1  )) then

            apws_mt(i,j)=apws_mt(i,j)
            apw_mt (i,j)=apw_mt (i,j)
            apwn_mt(i,j)=apwn_mt(i,j)

            aps_mt (i,j)=aps_mt (i,j)
            apn_mt (i,j)=apn_mt (i,j)

            apes_mt(i,j)=apes_mt(i,j)
            apen_mt(i,j)=apen_mt(i,j)

            vxie=(vp_mt(i,j)+vp_mt(i+1,j))/2.

            a1=suxix_mt(i+1,j)
            a2=suxiy_mt(i+1,j)

            b1=(suxix_mt(i+1,j)+suxix_mt(i+2,j))/2.
            b2=(suxiy_mt(i+1,j)+suxiy_mt(i+2,j))/2.

            ape_mt(i,j)=ape_mt(i,j)+(a1*b1+a2*b2)/vxie

            b1=(suxix_mt(i,j)+suxix_mt(i+1,j))/2.
            b2=(suxiy_mt(i,j)+suxiy_mt(i+1,j))/2.

            app_mt(i,j)=app_mt(i,j)-(a1*b1+a2*b2)/vxie
!-----------
          else

            apws_mt(i,j)=apws_mt(i,j)
            apw_mt (i,j)=apw_mt (i,j)
            apwn_mt(i,j)=apwn_mt(i,j)

            vxie=(vp_mt(i,j)+vp_mt(i+1,j))/2.

            a1=suxix_mt(i+1,j)
            a2=suxiy_mt(i+1,j)

            b1=(suxix_mt(i+1,j)+suxix_mt(i+2,j))/2.
            b2=(suxiy_mt(i+1,j)+suxiy_mt(i+2,j))/2.
            ape_mt(i,j)=ape_mt(i,j)+(a1*b1+a2*b2)/vxie

            b1=(suxix_mt(i,j)+suxix_mt(i+1,j))/2.
            b2=(suxiy_mt(i,j)+suxiy_mt(i+1,j))/2.
            app_mt(i,j)=app_mt(i,j)-(a1*b1+a2*b2)/vxie

            b1=(svetax_mt(i,j+1)+svetax_mt(i+1,j+1))/2./4.0
            b2=(svetay_mt(i,j+1)+svetay_mt(i+1,j+1))/2./4.0
            apn_mt (i,j)=apn_mt (i,j)+(a1*b1+a2*b2)/vxie
            app_mt (i,j)=app_mt (i,j)+(a1*b1+a2*b2)/vxie
            ape_mt (i,j)=ape_mt (i,j)+(a1*b1+a2*b2)/vxie
            apen_mt(i,j)=apen_mt(i,j)+(a1*b1+a2*b2)/vxie

            b1=(svetax_mt(i,j)+svetax_mt(i+1,j))/2./4.0
            b2=(svetay_mt(i,j)+svetay_mt(i+1,j))/2./4.0
            aps_mt (i,j)=aps_mt (i,j)-(a1*b1+a2*b2)/vxie
            app_mt (i,j)=app_mt (i,j)-(a1*b1+a2*b2)/vxie
            ape_mt (i,j)=ape_mt (i,j)-(a1*b1+a2*b2)/vxie
            apes_mt(i,j)=apes_mt(i,j)-(a1*b1+a2*b2)/vxie

          end if
!-----------
        end if
 31   continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute a_ij associated with R_xi,w      (page 37)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 32 j=1,nny
      do 32 i=1,nnx
!-----
        if(iblock.eq.1.and.i.eq.1) then
          apws_mt(i,j)=apws_mt(i,j)
          apw_mt (i,j)=apw_mt (i,j)
          apwn_mt(i,j)=apwn_mt(i,j)

          aps_mt(i,j)=aps_mt(i,j)
          app_mt(i,j)=app_mt(i,j)
          apn_mt(i,j)=apn_mt(i,j)

          apes_mt(i,j)=apes_mt(i,j)
          ape_mt (i,j)=ape_mt (i,j)
          apen_mt(i,j)=apen_mt(i,j)
!-----
        else
!-----
          if((jblock.eq.npy.and.j.eq.nny).or. &
             (jblock.eq.1  .and.j.eq.1  )) then

            apes_mt(i,j)=apes_mt(i,j)
            ape_mt (i,j)=ape_mt (i,j)
            apen_mt(i,j)=apen_mt(i,j)

            aps_mt (i,j)=aps_mt (i,j)
            apn_mt (i,j)=apn_mt (i,j)
	    apws_mt(i,j)=apws_mt(i,j) 
	    apwn_mt(i,j)=apwn_mt(i,j)

            vxiw=-(vp_mt(i,j)+vp_mt(i-1,j))/2.

            a1=suxix_mt(i,j)
            a2=suxiy_mt(i,j)

            b1=(suxix_mt(i,j)+suxix_mt(i+1,j))/2.
            b2=(suxiy_mt(i,j)+suxiy_mt(i+1,j))/2.

            app_mt(i,j)=app_mt(i,j)+(a1*b1+a2*b2)/vxiw

            b1=(suxix_mt(i,j)+suxix_mt(i-1,j))/2.
            b2=(suxiy_mt(i,j)+suxiy_mt(i-1,j))/2.

            apw_mt(i,j)=apw_mt(i,j)-(a1*b1+a2*b2)/vxiw
!-------------
          else

            apes_mt(i,j)=apes_mt(i,j)
            ape_mt (i,j)=ape_mt (i,j)
            apen_mt(i,j)=apen_mt(i,j)

            vxiw=-(vp_mt(i,j)+vp_mt(i-1,j))/2.

            a1=suxix_mt(i,j)
            a2=suxiy_mt(i,j)

            b1=(suxix_mt(i,j)+suxix_mt(i+1,j))/2.
            b2=(suxiy_mt(i,j)+suxiy_mt(i+1,j))/2.
            app_mt(i,j)=app_mt(i,j)+(a1*b1+a2*b2)/vxiw

            b1=(suxix_mt(i,j)+suxix_mt(i-1,j))/2.
            b2=(suxiy_mt(i,j)+suxiy_mt(i-1,j))/2.
            apw_mt(i,j)=apw_mt(i,j)-(a1*b1+a2*b2)/vxiw

            b1=(svetax_mt(i,j+1)+svetax_mt(i-1,j+1))/2./4.0
            b2=(svetay_mt(i,j+1)+svetay_mt(i-1,j+1))/2./4.0
            apn_mt (i,j)=apn_mt (i,j)+(a1*b1+a2*b2)/vxiw
            app_mt (i,j)=app_mt (i,j)+(a1*b1+a2*b2)/vxiw
            apw_mt (i,j)=apw_mt (i,j)+(a1*b1+a2*b2)/vxiw
            apwn_mt(i,j)=apwn_mt(i,j)+(a1*b1+a2*b2)/vxiw

            b1=(svetax_mt(i,j)+svetax_mt(i-1,j))/2./4.0
            b2=(svetay_mt(i,j)+svetay_mt(i-1,j))/2./4.0
            aps_mt (i,j)=aps_mt (i,j)-(a1*b1+a2*b2)/vxiw
            app_mt (i,j)=app_mt (i,j)-(a1*b1+a2*b2)/vxiw
            apw_mt (i,j)=apw_mt (i,j)-(a1*b1+a2*b2)/vxiw
            apws_mt(i,j)=apws_mt(i,j)-(a1*b1+a2*b2)/vxiw

          end if
!-------------
        end if
 32   continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute a_ij associated with R_eta,n    (page 37)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 33 j=1,nny
      do 33 i=1,nnx
!----
        if(jblock.eq.npy.and.j.eq.nny) then
          apws_mt(i,j)=apws_mt(i,j)
          apw_mt (i,j)=apw_mt (i,j)
          apwn_mt(i,j)=apwn_mt(i,j)

          aps_mt (i,j)=aps_mt (i,j)

          apes_mt(i,j)=apes_mt(i,j)
          ape_mt (i,j)=ape_mt (i,j)
          apen_mt(i,j)=apen_mt(i,j)

!----------
        else

          if((iblock.eq.npx.and.i.eq.nnx).or. &
             (iblock.eq.1  .and.i.eq.1  )) then

            apws_mt(i,j)=apws_mt(i,j)
            aps_mt (i,j)=aps_mt (i,j)
            apes_mt(i,j)=apes_mt(i,j)

            apw_mt (i,j)=apw_mt (i,j)
            ape_mt (i,j)=ape_mt (i,j)

            apwn_mt(i,j)=apwn_mt(i,j)
            apen_mt(i,j)=apen_mt(i,j)

            vetan=(vp_mt(i,j)+vp_mt(i,j+1))/2.

            a1=svetax_mt(i,j+1)
            a2=svetay_mt(i,j+1)

            b1=(svetax_mt(i,j+1)+svetax_mt(i,j+2))/2.
            b2=(svetay_mt(i,j+1)+svetay_mt(i,j+2))/2.
            apn_mt(i,j)=apn_mt(i,j)+(a1*b1+a2*b2)/vetan

            b1=(svetax_mt(i,j)+svetax_mt(i,j+1))/2.
            b2=(svetay_mt(i,j)+svetay_mt(i,j+1))/2.
            app_mt(i,j)=app_mt(i,j)-(a1*b1+a2*b2)/vetan
!--------
          else

            apws_mt(i,j)=apws_mt(i,j)
            aps_mt (i,j)=aps_mt (i,j)
            apes_mt(i,j)=apes_mt(i,j)

            vetan=(vp_mt(i,j)+vp_mt(i,j+1))/2.

            a1=svetax_mt(i,j+1)
            a2=svetay_mt(i,j+1)

            b1=(svetax_mt(i,j+1)+svetax_mt(i,j+2))/2.
            b2=(svetay_mt(i,j+1)+svetay_mt(i,j+2))/2.

            apn_mt(i,j)=apn_mt(i,j)+(a1*b1+a2*b2)/vetan

            b1=(svetax_mt(i,j)+svetax_mt(i,j+1))/2.
            b2=(svetay_mt(i,j)+svetay_mt(i,j+1))/2.
            app_mt(i,j)=app_mt(i,j)-(a1*b1+a2*b2)/vetan

            b1=(suxix_mt(i+1,j)+suxix_mt(i+1,j+1))/2./4.0
            b2=(suxiy_mt(i+1,j)+suxiy_mt(i+1,j+1))/2./4.0
            apn_mt (i,j)=apn_mt (i,j)+(a1*b1+a2*b2)/vetan
            app_mt (i,j)=app_mt (i,j)+(a1*b1+a2*b2)/vetan
            ape_mt (i,j)=ape_mt (i,j)+(a1*b1+a2*b2)/vetan
            apen_mt(i,j)=apen_mt(i,j)+(a1*b1+a2*b2)/vetan

            b1=(suxix_mt(i,j)+suxix_mt(i,j+1))/2./4.0
            b2=(suxiy_mt(i,j)+suxiy_mt(i,j+1))/2./4.0
            apn_mt (i,j)=apn_mt (i,j)-(a1*b1+a2*b2)/vetan
            app_mt (i,j)=app_mt (i,j)-(a1*b1+a2*b2)/vetan
            apw_mt (i,j)=apw_mt (i,j)-(a1*b1+a2*b2)/vetan
            apwn_mt(i,j)=apwn_mt(i,j)-(a1*b1+a2*b2)/vetan

          end if
!--------
        end if
 33   continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute a_ij associated with R_eta,s      (page 37)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 34 j=1,nny
      do 34 i=1,nnx
!----
        if(jblock.eq.1.and.j.eq.1) then

          apws_mt(i,j)=apws_mt(i,j)
          apw_mt (i,j)=apw_mt (i,j)
          apwn_mt(i,j)=apwn_mt(i,j)

          apn_mt (i,j)=apn_mt (i,j)

          apes_mt(i,j)=apes_mt(i,j)
          ape_mt (i,j)=ape_mt (i,j)
          apen_mt(i,j)=apen_mt(i,j)

!-----
        else

          if((iblock.eq.npx.and.i.eq.nnx).or. &
             (iblock.eq.1  .and.i.eq.1  )) then

            apwn_mt(i,j)=apwn_mt(i,j)
            apn_mt (i,j)=apn_mt (i,j)
            apen_mt(i,j)=apen_mt(i,j)

            apw_mt (i,j)=apw_mt (i,j)
            ape_mt (i,j)=ape_mt (i,j)

            apws_mt(i,j)=apws_mt(i,j)
            apes_mt(i,j)=apes_mt(i,j)

            vetas=-(vp_mt(i,j)+vp_mt(i,j-1))/2.

            a1=svetax_mt(i,j)
            a2=svetay_mt(i,j)

            b1=(svetax_mt(i,j)+svetax_mt(i,j+1))/2.
            b2=(svetay_mt(i,j)+svetay_mt(i,j+1))/2.
            app_mt(i,j)=app_mt(i,j)+(a1*b1+a2*b2)/vetas

            b1=(svetax_mt(i,j)+svetax_mt(i,j-1))/2.
            b2=(svetay_mt(i,j)+svetay_mt(i,j-1))/2.
            aps_mt(i,j)=aps_mt(i,j)-(a1*b1+a2*b2)/vetas
!-------------
          else

            apwn_mt(i,j)=apwn_mt(i,j)
            apn_mt (i,j)=apn_mt (i,j)
            apen_mt(i,j)=apen_mt(i,j)

            vetas=-(vp_mt(i,j)+vp_mt(i,j-1))/2.

            a1=svetax_mt(i,j)
            a2=svetay_mt(i,j)

            b1=(svetax_mt(i,j)+svetax_mt(i,j+1))/2.
            b2=(svetay_mt(i,j)+svetay_mt(i,j+1))/2.
            app_mt(i,j)=app_mt(i,j)+(a1*b1+a2*b2)/vetas

            b1=(svetax_mt(i,j)+svetax_mt(i,j-1))/2.
            b2=(svetay_mt(i,j)+svetay_mt(i,j-1))/2.
            aps_mt(i,j)=aps_mt(i,j)-(a1*b1+a2*b2)/vetas

            b1=(suxix_mt(i+1,j)+suxix_mt(i+1,j-1))/2./4.0
            b2=(suxiy_mt(i+1,j)+suxiy_mt(i+1,j-1))/2./4.0
            aps_mt (i,j)=aps_mt (i,j)+(a1*b1+a2*b2)/vetas
            app_mt (i,j)=app_mt (i,j)+(a1*b1+a2*b2)/vetas
            ape_mt (i,j)=ape_mt (i,j)+(a1*b1+a2*b2)/vetas
            apes_mt(i,j)=apes_mt(i,j)+(a1*b1+a2*b2)/vetas

            b1=(suxix_mt(i,j)+suxix_mt(i,j-1))/2./4.0
            b2=(suxiy_mt(i,j)+suxiy_mt(i,j-1))/2./4.0
            aps_mt (i,j)=aps_mt (i,j)-(a1*b1+a2*b2)/vetas
            app_mt (i,j)=app_mt (i,j)-(a1*b1+a2*b2)/vetas
            apw_mt (i,j)=apw_mt (i,j)-(a1*b1+a2*b2)/vetas
            apws_mt(i,j)=apws_mt(i,j)-(a1*b1+a2*b2)/vetas

          end if
!-------------
        end if
 34   continue

      return
      end
!===============================================================================
