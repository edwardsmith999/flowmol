subroutine copy_dct(aout,ain,nnx,nny)
	implicit real (a-h,o-z)
	dimension ain(0:nnx+1,0:nny+1),aout(0:nnx+1,0:nny+1)
	aout(:,:)=ain(:,:)
	return
end


subroutine copyap_dct(aout,ain,nnx,nny)
	implicit real (a-h,o-z)
	dimension ain(nnx,nny),aout(nnx,nny)
	aout(:,:)=ain(:,:)
	return
end


subroutine copysu_dct(aout,ain,nnx,nny)
	implicit real (a-h,o-z)
	dimension ain(0:nnx+2,0:nny+1),aout(0:nnx+2,0:nny+1)
	aout(:,:)=ain(:,:)
	return
end


subroutine copysv_dct(aout,ain,nnx,nny)
	implicit real (a-h,o-z)
	dimension ain(0:nnx+1,0:nny+2),aout(0:nnx+1,0:nny+2)
	aout(:,:)=ain(:,:)
	return
end


subroutine preprocess_dct 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! allocate memory for Poisson coefficients
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	use data_export
	dimension zsuv(memlenap_zsuv)
	dimension isux_mt(ngrid),isuy_mt(ngrid),isvx_mt(ngrid),isvy_mt(ngrid)

	if (irank.eq.1) then
		write(*,*) 'WARNING: dct int arrays: ngrid  ,nglevel_max  =',ngrid,nglevel_max 
	end if
	memap=0
	m    =ngrid
	nnx  =ngxm 
	nny  =ngym 
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! pointer for coefficients in poisson equation
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	iaps_mt (m)=malocap_dct(nnx*nny        ,memap,memlenap_zap)
	iapn_mt (m)=malocap_dct(nnx*nny        ,memap,memlenap_zap)
	iapp_mt (m)=malocap_dct(nnx*nny        ,memap,memlenap_zap)
	ivp_mt  (m)=malocap_dct((nnx+2)*(nny+2),memap,memlenap_zap)
	iedge_mt(m)=malocap_dct(2              ,memap,memlenap_zap)
	if (irank.eq.1) then
		write(*,*) 'WARNING: zap_dct       : memap+1,memlenap_zap =',memap+1,memlenap_zap 
	end if

	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	! temporary pointer for metric
	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	memap=0
	isux_mt(m)=malocap_dct((nnx+3)*(nny+2),memap,memlenap_zsuv)
	isuy_mt(m)=malocap_dct((nnx+3)*(nny+2),memap,memlenap_zsuv)
	isvx_mt(m)=malocap_dct((nnx+2)*(nny+3),memap,memlenap_zsuv)
	isvy_mt(m)=malocap_dct((nnx+2)*(nny+3),memap,memlenap_zsuv)
	if (irank.eq.1) then
		write(*,*) 'WARNING: zsuv_dct      : memap+1,memlenap_zsuv=',memap+1,memlenap_zsuv
	end if

	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	! pointer for solution, source term
	!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	memap=0
	ip_mt  (m)=malocap_dct((nnx+2)*(nny+2),memap,memlenap_z)
	irhs_mt(m)=malocap_dct((nnx+2)*(nny+2),memap,memlenap_z)
	if (irank.eq.1) then
		write(*,*) 'WARNING: z_dct         : memap+1,memlenap_z   =',memap+1,memlenap_z
	end if

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! assign area vector and volume
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
	call copy_dct  (zap_mt(ivp_mt (m)),vp    ,nnx,nny)
	call copysu_dct(zsuv  (isux_mt(m)),suxix ,nnx,nny)
	call copysu_dct(zsuv  (isuy_mt(m)),suxiy ,nnx,nny)
	call copysv_dct(zsuv  (isvx_mt(m)),svetax,nnx,nny)
	call copysv_dct(zsuv  (isvy_mt(m)),svetay,nnx,nny)

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! evaluate Poisson coefficients
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	call pcoef_dct(	zap_mt(iaps_mt (m)),zap_mt(iapn_mt (m)), &
			zap_mt(iapp_mt (m)),zap_mt(ivp_mt  (m)), &
			zsuv  (isux_mt (m)),zsuv  (isuy_mt (m)), &
			zsuv  (isvx_mt (m)),zsuv  (isvy_mt (m)), &
			zap_mt(iedge_mt(m)),nnx,nny) 

      return
end

function malocap_dct(lenap,memap,memlenap)
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
	malocap_dct=memap+1 
	memap=memap+lenap
	return
end

subroutine pcoef_dct(	aps_mt  ,apn_mt  , app_mt   ,vp_mt    , &
                        suxix_mt,suxiy_mt, svetax_mt,svetay_mt, &
                        aedgex  ,nnx,nny)
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! compute the Poisson equation coefficients
	! with (nnx, nny) points
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real	aps_mt (nnx,nny), apn_mt(nnx,nny), app_mt(nnx,nny),            &
		suxix_mt (0:nnx+1+1,0:nny+1  ),suxiy_mt (0:nnx+1+1,0:nny+1  ), &
		svetax_mt(0:nnx+1  ,0:nny+1+1),svetay_mt(0:nnx+1  ,0:nny+1+1), &
		aedgex(2),vp_mt(0:nnx+1,0:nny+1)
	integer iedgex(2)
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! find the leading edge and trailing edge i index
	! pressure side and suction side have same index
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	iedgex(1)=1
	iedgex(2)=nnx
	aedgex   =1.0*iedgex

	!=================================================
	!     Initialize all the Poisson Coefficients
	!=================================================
	aps_mt  = 0.0
	app_mt  = 0.0
	apn_mt  = 0.0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute a_ij associated with R_eta,n    (page 37)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 33 j=1,nny
      do 33 i=1,nnx
        if(j.eq.nny) then
          aps_mt(i,j)=aps_mt(i,j)
!----------
        else 
            aps_mt(i,j) =aps_mt(i,j)

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
        end if         
 33   continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute a_ij associated with R_eta,s      (page 37)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 34 j=1,nny
      do 34 i=1,nnx
        if(j.eq.1) then
          apn_mt(i,j)=apn_mt(i,j)
!-----
        else 
            apn_mt(i,j)=apn_mt(i,j)

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
        end if         
 34   continue

      return
      end
