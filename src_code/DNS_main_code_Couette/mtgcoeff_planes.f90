subroutine copy(aout,ain,nnx,nny)
	implicit real (a-h,o-z)
	dimension ain(0:nnx+1,0:nny+1),aout(0:nnx+1,0:nny+1)
	aout(:,:)=ain(:,:)
	return
end


subroutine copyap(aout,ain,nnx,nny)
	implicit real (a-h,o-z)
	dimension ain(nnx,nny),aout(nnx,nny)
	aout(:,:)=ain(:,:)
	return
end


subroutine copysu(aout,ain,nnx,nny)
	implicit real (a-h,o-z)
	dimension ain(0:nnx+2,0:nny+1),aout(0:nnx+2,0:nny+1)
	aout(:,:)=ain(:,:)
	return
end


subroutine copysv(aout,ain,nnx,nny)
	implicit real (a-h,o-z)
	dimension ain(0:nnx+1,0:nny+2),aout(0:nnx+1,0:nny+2)
	aout(:,:)=ain(:,:)
	return
end


subroutine restrivp(pc_mt,pf_mt,ncx,ncy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! restriction :  determine coarse grid volume values using
!-----------------------------------------------------------------------
!                fine grid volume values
!                     pf = finer pressure cell
!                     pc = coarser pressure cell
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit real (a-h,o-z)
	dimension pc_mt(0:ncx+1,0:ncy+1),pf_mt(0:2*ncx+1,0:2*ncy+1)

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
	! lower and upper b.c.
	!cccccccccccccccccccccccccc
	suc_mt(1:ncx+1,0    )=suc_mt(1:ncx+1,1)
	suc_mt(1:ncx+1,ncy+1)=suc_mt(1:ncx+1,ncy)
	!cccccccccccccccccccccccccc
	! left and right b.c.
	!cccccccccccccccccccccccccc
	suc_mt(0    ,:)=suc_mt(2,:)
	suc_mt(ncx+2,:)=suc_mt(ncx,:)
	return
end


subroutine restrisv(svc_mt,svf_mt,ncx,ncy)
	implicit real (a-h,o-z)
	dimension svc_mt(0:ncx+1,0:ncy+1+1),svf_mt(0:2*ncx+1,0:2*ncy+1+1)

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
	svc_mt(0    ,1:ncy+1)=svc_mt(1,1:ncy+1)
	svc_mt(ncx+1,1:ncy+1)=svc_mt(ncx,1:ncy+1)
	!cccccccccccccccccccccccccc
	! upper and lower b.c.
	!cccccccccccccccccccccccccc
	svc_mt(:,0    )=svc_mt(:,2)
	svc_mt(:,ncy+2)=svc_mt(:,ncy)
	return
end



subroutine preprocess_mt 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! allocate memory for Poisson coefficients at all levels
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	use data_export
	dimension zsuv(memlenap_zsuv)
	dimension isux_mt(ngrid),isuy_mt(ngrid),isvx_mt(ngrid),isvy_mt(ngrid)

	!TAZ  if (irank.eq.1)   &
		write(*,*) 'WARNING: mt int arrays: ngrid  ,nglevel_max  =',ngrid,nglevel_max 
	memap=0
	nnx=ngxm 
	nny=ngym 
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! pointer for coefficients in poisson equation at all levels
	!----  from fine to coarse ----
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do m = ngrid,1,-1
		iapw_mt (m)=malocap(nnx*nny        ,memap,memlenap_zap)
		iape_mt (m)=malocap(nnx*nny        ,memap,memlenap_zap)
		iaps_mt (m)=malocap(nnx*nny        ,memap,memlenap_zap)
		iapn_mt (m)=malocap(nnx*nny        ,memap,memlenap_zap)
		iapws_mt(m)=malocap(nnx*nny        ,memap,memlenap_zap)
		iapwn_mt(m)=malocap(nnx*nny        ,memap,memlenap_zap)
		iapes_mt(m)=malocap(nnx*nny        ,memap,memlenap_zap)
		iapen_mt(m)=malocap(nnx*nny        ,memap,memlenap_zap)
		iapp_mt (m)=malocap(nnx*nny        ,memap,memlenap_zap)
		ivp_mt  (m)=malocap((nnx+2)*(nny+2),memap,memlenap_zap)
		iedge_mt(m)=malocap(2              ,memap,memlenap_zap)
		nnx=nnx/2
		nny=nny/2
	end do
	!TAZ  if (irank.eq.1)   &
		write(*,*) 'WARNING: zap_mt       : memap+1,memlenap_zap =',memap+1,memlenap_zap 

	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	! temporary pointer for metric at all levels
	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	memap=0
	nnx=ngxm
	nny=ngym
	do m = ngrid,1,-1
		isux_mt(m)=malocap((nnx+3)*(nny+2),memap,memlenap_zsuv)
		isuy_mt(m)=malocap((nnx+3)*(nny+2),memap,memlenap_zsuv)
		isvx_mt(m)=malocap((nnx+2)*(nny+3),memap,memlenap_zsuv)
		isvy_mt(m)=malocap((nnx+2)*(nny+3),memap,memlenap_zsuv)
		nnx=nnx/2
		nny=nny/2
	end do
	!TAZ  if (irank.eq.1)  &
		write(*,*) 'WARNING: zsuv         : memap+1,memlenap_zsuv=',memap+1,memlenap_zsuv

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! pointer for solution, source term and residual at all levels
	! Xiahua defines these in 'poisson_newmtg.f90'
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	memap=0
	nnx=ngxm
	nny=ngym
	do m = ngrid,1,-1
		ip_mt  (m)=malocap((nnx+2)*(nny+2),memap,memlenap_z)
		irhs_mt(m)=malocap((nnx+2)*(nny+2),memap,memlenap_z)
		ires_mt(m)=malocap((nnx+2)*(nny+2),memap,memlenap_z)
		nnx=nnx/2
		nny=nny/2
	enddo
	!TAZ  if (irank.eq.1)  &
		write(*,*) 'WARNING: z_mt         : memap+1,memlenap_z   =',memap+1,memlenap_z


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! assign (restrict) area vector and volume at all levels
!-----------------------------------------------------------------------
	!ccccccccccccccccccccccccccccccccccccccccc
	! assign at the finest level
	!ccccccccccccccccccccccccccccccccccccccccc
	m = ngrid
	nnx = ngxm
	nny = ngym
	call copy  (zap_mt(ivp_mt (m)),vp    ,nnx,nny)
	call copysu(zsuv  (isux_mt(m)),suxix ,nnx,nny)
	call copysu(zsuv  (isuy_mt(m)),suxiy ,nnx,nny)
	call copysv(zsuv  (isvx_mt(m)),svetax,nnx,nny)
	call copysv(zsuv  (isvy_mt(m)),svetay,nnx,nny)
	!ccccccccccccccccccccccccccccccccccccccccc
	! restrict at lower levels
	!ccccccccccccccccccccccccccccccccccccccccc
	nnx=ngxm
	nny=ngym
	do m=ngrid-1,1,-1
  	nnx=nnx/2
  	nny=nny/2
  		call restrivp(zap_mt(ivp_mt (m)),zap_mt(ivp_mt (m+1)),nnx,nny)
  		call restrisu(zsuv  (isux_mt(m)),zsuv  (isux_mt(m+1)),nnx,nny)
  		call restrisu(zsuv  (isuy_mt(m)),zsuv  (isuy_mt(m+1)),nnx,nny)
  		call restrisv(zsuv  (isvx_mt(m)),zsuv  (isvx_mt(m+1)),nnx,nny)
  		call restrisv(zsuv  (isvy_mt(m)),zsuv  (isvy_mt(m+1)),nnx,nny)
	end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! evaluate Poisson coefficients at all levels
!---------------------------------------------------
! from fine to coarse
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	nnx=ngxm
	nny=ngym
	do m=ngrid,1,-1
		call pcoef_mt(	zap_mt(iapw_mt (m)),zap_mt(iape_mt (m)), &
				zap_mt(iaps_mt (m)),zap_mt(iapn_mt (m)), &
				zap_mt(iapws_mt(m)),zap_mt(iapwn_mt(m)), &
				zap_mt(iapes_mt(m)),zap_mt(iapen_mt(m)), &
				zap_mt(iapp_mt (m)),zap_mt(ivp_mt  (m)), &
				zsuv  (isux_mt (m)),zsuv  (isuy_mt (m)), &
				zsuv  (isvx_mt (m)),zsuv  (isvy_mt (m)), &
				zap_mt(iedge_mt(m)),nnx,nny) 
		nnx=nnx/2
		nny=nny/2
	end do

      return
end






function malocap(lenap,memap,memlenap)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! dynamical storage allocation, returns integer pointer to the
! starting position for lenap array elements in the array z.
! the malocap-1'th element in z contains the real value of lenap,
!-----------------------------------------------------------------------
!      lenap   : Size of array requested
!      memap   : Amount of memory already used/exhausted 
!                in the allocation stack
!      memlenap: Total memory of allocation stack
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	if(memap+lenap+1.gt.memlenap) stop 'inspf_mtficient memory in malocap'
	malocap=memap+1 
	memap=memap+lenap
	return
end




subroutine pcoef_mt(	apw_mt   ,ape_mt   ,aps_mt  ,apn_mt  , &
			apws_mt  ,apwn_mt  ,apes_mt ,apen_mt , &
			app_mt   ,vp_mt    ,suxix_mt,suxiy_mt, &
			svetax_mt,svetay_mt,aedgex  ,nnx,nny)
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! compute the Poisson equation coefficients at one
	! multigrid level with (nnx, nny) points
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real	apw_mt(nnx,nny),ape_mt (nnx,nny),aps_mt (nnx,nny), &
		apn_mt(nnx,nny),apws_mt(nnx,nny),apwn_mt(nnx,nny), &
		apes_mt(nnx,nny),apen_mt(nnx,nny),app_mt (nnx,nny), &
		suxix_mt(0:nnx+1+1,0:nny+1)  ,suxiy_mt (0:nnx+1+1,0:nny+1), &
		svetax_mt(0:nnx+1  ,0:nny+1+1),svetay_mt(0:nnx+1  ,0:nny+1+1), &
		aedgex(2),vp_mt(0:nnx+1,0:nny+1)
	integer iedgex(2)
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! find the leading edge and trailing edge i index
	! pressure side and suction side have same index
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	iedgex(1)=1
	iedgex(2)=nnx
	aedgex   = 1.0*iedgex

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
        if(i.eq.nnx) then
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
          if(j.eq.nny.or.j.eq.1) then

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
            apw_mt(i,j)=apw_mt(i,j)
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
            apn_mt(i,j)=apn_mt(i,j)+(a1*b1+a2*b2)/vxie
            app_mt(i,j)=app_mt(i,j)+(a1*b1+a2*b2)/vxie
            ape_mt(i,j)=ape_mt(i,j)+(a1*b1+a2*b2)/vxie
            apen_mt(i,j)=apen_mt(i,j)+(a1*b1+a2*b2)/vxie

            b1=(svetax_mt(i,j)+svetax_mt(i+1,j))/2./4.0
            b2=(svetay_mt(i,j)+svetay_mt(i+1,j))/2./4.0
            aps_mt(i,j)=aps_mt(i,j)-(a1*b1+a2*b2)/vxie
            app_mt(i,j)=app_mt(i,j)-(a1*b1+a2*b2)/vxie
            ape_mt(i,j)=ape_mt(i,j)-(a1*b1+a2*b2)/vxie
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
        if(i.eq.1) then
          apws_mt(i,j)=apws_mt(i,j)
          apw_mt(i,j) =apw_mt(i,j)
          apwn_mt(i,j)=apwn_mt(i,j)

          aps_mt(i,j)=aps_mt(i,j)
          app_mt(i,j)=app_mt(i,j)
          apn_mt(i,j)=apn_mt(i,j)

          apes_mt(i,j)=apes_mt(i,j)
          ape_mt(i,j) =ape_mt(i,j)
          apen_mt(i,j)=apen_mt(i,j)
!-----
        else 
!-----
          if(j.eq.nny.or.j.eq.1) then

            apes_mt(i,j)=apes_mt(i,j)
            ape_mt(i,j)=ape_mt(i,j)
            apen_mt(i,j)=apen_mt(i,j)

            aps_mt(i,j)=aps_mt(i,j)
            apn_mt(i,j)=apn_mt(i,j)

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
            ape_mt(i,j)=ape_mt(i,j)
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
            apn_mt(i,j)=apn_mt(i,j)+(a1*b1+a2*b2)/vxiw 
            app_mt(i,j)=app_mt(i,j)+(a1*b1+a2*b2)/vxiw 
            apw_mt(i,j)=apw_mt(i,j)+(a1*b1+a2*b2)/vxiw 
            apwn_mt(i,j)=apwn_mt(i,j)+(a1*b1+a2*b2)/vxiw 

            b1=(svetax_mt(i,j)+svetax_mt(i-1,j))/2./4.0
            b2=(svetay_mt(i,j)+svetay_mt(i-1,j))/2./4.0
            aps_mt(i,j)=aps_mt(i,j)-(a1*b1+a2*b2)/vxiw 
            app_mt(i,j)=app_mt(i,j)-(a1*b1+a2*b2)/vxiw 
            apw_mt(i,j)=apw_mt(i,j)-(a1*b1+a2*b2)/vxiw 
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
        if(j.eq.nny) then
          apws_mt(i,j)=apws_mt(i,j)
          apw_mt(i,j)=apw_mt(i,j)
          apwn_mt(i,j)=apwn_mt(i,j)

          aps_mt(i,j)=aps_mt(i,j)

          apes_mt(i,j)=apes_mt(i,j)
          ape_mt(i,j)=ape_mt(i,j)
          apen_mt(i,j)=apen_mt(i,j)

!----------
        else 

          if(i.eq.nnx.or.i.eq.1) then

            apws_mt(i,j)=apws_mt(i,j)
            aps_mt(i,j)=aps_mt(i,j)
            apes_mt(i,j)=apes_mt(i,j)

            apw_mt(i,j)=apw_mt(i,j)
            ape_mt(i,j)=ape_mt(i,j)

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
            aps_mt(i,j)=aps_mt(i,j)
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
            apn_mt(i,j)=apn_mt(i,j)+(a1*b1+a2*b2)/vetan
            app_mt(i,j)=app_mt(i,j)+(a1*b1+a2*b2)/vetan
            ape_mt(i,j)=ape_mt(i,j)+(a1*b1+a2*b2)/vetan
            apen_mt(i,j)=apen_mt(i,j)+(a1*b1+a2*b2)/vetan

            b1=(suxix_mt(i,j)+suxix_mt(i,j+1))/2./4.0
            b2=(suxiy_mt(i,j)+suxiy_mt(i,j+1))/2./4.0
            apn_mt(i,j)=apn_mt(i,j)-(a1*b1+a2*b2)/vetan
            app_mt(i,j)=app_mt(i,j)-(a1*b1+a2*b2)/vetan
            apw_mt(i,j)=apw_mt(i,j)-(a1*b1+a2*b2)/vetan
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
        if(j.eq.1) then
          apws_mt(i,j)=apws_mt(i,j)
          apw_mt(i,j)=apw_mt(i,j)
          apwn_mt(i,j)=apwn_mt(i,j)

          apn_mt(i,j)=apn_mt(i,j)

          apes_mt(i,j)=apes_mt(i,j)
          ape_mt(i,j)=ape_mt(i,j)
          apen_mt(i,j)=apen_mt(i,j)

!-----
        else 

          if(i.eq.nnx.or.i.eq.1) then

            apwn_mt(i,j)=apwn_mt(i,j)
            apn_mt(i,j)=apn_mt(i,j)
            apen_mt(i,j)=apen_mt(i,j)

            apw_mt(i,j)=apw_mt(i,j)
            ape_mt(i,j)=ape_mt(i,j)

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
            apn_mt(i,j)=apn_mt(i,j)
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
            aps_mt(i,j)=aps_mt(i,j)+(a1*b1+a2*b2)/vetas
            app_mt(i,j)=app_mt(i,j)+(a1*b1+a2*b2)/vetas
            ape_mt(i,j)=ape_mt(i,j)+(a1*b1+a2*b2)/vetas
            apes_mt(i,j)=apes_mt(i,j)+(a1*b1+a2*b2)/vetas

            b1=(suxix_mt(i,j)+suxix_mt(i,j-1))/2./4.0
            b2=(suxiy_mt(i,j)+suxiy_mt(i,j-1))/2./4.0
            aps_mt(i,j)=aps_mt(i,j)-(a1*b1+a2*b2)/vetas
            app_mt(i,j)=app_mt(i,j)-(a1*b1+a2*b2)/vetas
            apw_mt(i,j)=apw_mt(i,j)-(a1*b1+a2*b2)/vetas
            apws_mt(i,j)=apws_mt(i,j)-(a1*b1+a2*b2)/vetas

          end if
!-------------
        end if         
 34   continue

      return
      end
