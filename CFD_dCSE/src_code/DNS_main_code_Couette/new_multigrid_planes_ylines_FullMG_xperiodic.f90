

subroutine relax_b(	p_mt,rhs_mt,res_mt,nnx,nny,k_mt, &
			apw_mt,ape_mt,aps_mt,apn_mt,apws_mt,apwn_mt,apes_mt,apen_mt,app_mt, &
			vp_mt,aedgex,omega_mt,vkz_mt,p_mthat,a_mt,b_mt,c_mt,d_mt,f_mt,g_mt, &
			iprob_mt,resminb_mt,nitmaxb_mt)

	implicit real (a-h,o-z)
	dimension 	rhs_mt(0:nnx+1,0:nny+1),p_mt(0:nnx+1,0:nny+1),res_mt(0:nnx+1,0:nny+1)
	dimension	apw_mt(nnx,nny),ape_mt(nnx,nny),aps_mt(nnx,nny),   &
			apn_mt(nnx,nny),apws_mt(nnx,nny),apwn_mt(nnx,nny),   &
			apes_mt(nnx,nny),apen_mt(nnx,nny),app_mt(nnx,nny),   &
			vp_mt(0:nnx+1,0:nny+1),aedgex(2)
	dimension 	p_mthat(0:nnx+1,0:nny+1) 
	dimension 	a_mt(nny),b_mt(nny),c_mt(nny),d_mt(nny),f_mt(nny),g_mt(nny)

	!  tmp:  Adds the coefficients at every grid point (adds equations)
	dimension	tmp(0:nnx+1,0:nny+1)
	! tmp_l:  Adds coefficients of a single equation
	real		tmp_l, tmp_lmax, tmp_max

	p_mthat=0.0
	res_mt=0.0
 993    format(1x,10(1pe10.3))
 331    continue
	icounter = 0 
 9000   icounter = icounter+1
	resb_max=-1.e7

	tmp   = 0.;	tmp_max = 0.;	imax = 0;	jmax = 0;
	tmp_l = 0.;	tmp_lmax= 0.;

	do i=1,nnx
		do j=1,nny

			if (k_mt.eq.1) then
				tmp_l =	  ape_mt(i,j)+ apw_mt (i,j)+ apn_mt(i,j)+ aps_mt (i,j)+ apws_mt(i,j)  &
					+ apwn_mt(i,j)+ apes_mt(i,j)+ apen_mt(i,j) +  (app_mt(i,j)-vp_mt(i,j)*vkz_mt);
				tmp_lmax = amax1(tmp_l,tmp_lmax);
				if (abs(tmp_l).gt.1e-12) then
					print*, 'tmp_l = a+b+c = ',tmp_l
					STOP " PROBLEM IN NEW MULTIGRID "
				end if
			end if

			if (k_mt.eq.1) then
				jloc = (i-1)*nny + i

				tmp( i , j ) =  tmp( i ,j) + app_mt(i,j)-vp_mt(i,j)*vkz_mt
				tmp(i-1, j ) =  tmp(i-1,j) + apw_mt (i,j)
				tmp(i+1, j ) =  tmp(i+1,j) + ape_mt (i,j)

				tmp(i-1,j+1) =  tmp(i-1,j+1) + apwn_mt(i,j)
				tmp( i ,j+1) =  tmp( i ,j+1) + apn_mt (i,j)
				tmp(i+1,j+1) =  tmp(i+1,j+1) + apen_mt(i,j)

				tmp(i-1,j-1) =  tmp(i-1,j-1) + apws_mt(i,j)
				tmp( i ,j-1) =  tmp( i ,j-1) + aps_mt (i,j)
				tmp(i+1,j-1) =  tmp(i+1,j-1) + apes_mt(i,j)
			end if


			res_mt(i,j)=-rhs_mt(i,j)+(ape_mt(i,j)*p_mt(i+1,j  )+apw_mt (i,j)*p_mt(i-1,j  ) &
						+ apn_mt(i,j)*p_mt(i  ,j+1)+aps_mt (i,j)*p_mt(i  ,j-1) &
						+ apws_mt(i,j)*p_mt(i-1,j-1)+apwn_mt(i,j)*p_mt(i-1,j+1) &
						+ apes_mt(i,j)*p_mt(i+1,j-1)+apen_mt(i,j)*p_mt(i+1,j+1) &
						+(app_mt(i,j)-vp_mt(i,j)*vkz_mt)*p_mt(i,j))
			d_mt(j)=omega_mt*(res_mt(i,j)+apws_mt(i,j)*p_mthat(i-1,j-1)+ &
				apw_mt(i,j)*p_mthat(i-1,j)+apwn_mt(i,j)*p_mthat(i-1,j+1))
			c_mt(j)=-apn_mt(i,j)
			b_mt(j)=-app_mt(i,j)+vp_mt(i,j)*vkz_mt
			a_mt(j)=-aps_mt(i,j)
			resb_max=amax1(abs(res_mt(i,j)),resb_max)
		end do
		!cccccccccccccccccccccccccccccccccccc
		!      Thomas algorithm begin
		!cccccccccccccccccccccccccccccccccccc
		a_mt(1) = 0.
		c_mt(nny) = 0.
		f_mt(1) = c_mt(1) / b_mt(1)
		g_mt(1) = d_mt(1) / b_mt(1)
		!cccccccccccccccccccccccccccccccccccc
		!       forward elimination
		!cccccccccccccccccccccccccccccccccccc
		do j = 2, nny 
    			f_mt(j) = c_mt(j) /(b_mt(j) - a_mt(j)*f_mt(j-1)) 
			g_mt(j) =(d_mt(j)-a_mt(j)*g_mt(j-1))/(b_mt(j)-a_mt(j)*f_mt(j-1))
		end do
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!   back substitution (solution put into d(j))
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		d_mt(nny) = g_mt(nny)
		do j = nny-1, 1, -1
			d_mt(j) = g_mt(j) - f_mt(j)*d_mt(j+1)
		end do
		!cccccccccccccccccccccccccccccccccccc
		!      Thomas algorithm end
		!cccccccccccccccccccccccccccccccccccc
		do j=1,nny
			p_mthat(i,j)=d_mt(j)
		end do
	end do

	!ccccccccccccccccccccccccccccccccccccccccccccccc
	! TJ: enforce x-periodicty on p_mthat here
	!	- (nnx.eq.ngxm), (nny.eq.ngym)
	!	- global indices => don't need MPI
	!ccccccccccccccccccccccccccccccccccccccccccccccc
	p_mthat(  0,  :) = 0.5*(p_mthat(0,:)   + p_mthat(nnx,:))
	p_mthat(nnx+1,:) = 0.5*(p_mthat(nnx+1,:) + p_mthat(1,:))

	if (k_mt.eq.1) then	
		do ii=0,nnx+1
		do jj=0,nny+1
				if (amax1( tmp_max,tmp(ii,jj) ) .gt. tmp_max) then
					tmp_max = amax1( tmp_max,tmp(ii,jj) )
					imax = ii;
					jmax = jj;
				end if
				if (abs(tmp(ii,jj)).gt.1e-12) then
					print*, 'tmp = a+b+c = ',ii,jj,tmp(ii,jj)
					STOP " PROBLEM IN NEW MULTIGRID "
				end if
		end do
		end do
	end if

	!ccccccccccccccc
	!    update
	!ccccccccccccccc
	do j=1,nny
	do i=1,nnx
		p_mt(i,j)=p_mt(i,j)+p_mthat(i,j)
	end do
	end do

	!ccccccccccccccccccccccccccccccccccccccccccccccc
	! TJ: enforce x-periodicty on p_mt here
	!	- (nnx.eq.ngxm), (nny.eq.ngym)
	!ccccccccccccccccccccccccccccccccccccccccccccccc
	p_mt(  0,  :) = p_mt(0,:)     + p_mthat(0,:)
	p_mt(nnx+1,:) = p_mt(nnx+1,:) + p_mthat(nnx+1,:)

	if(resb_max.le.resminb_mt.or.icounter.ge.nitmaxb_mt) goto 9001
	goto 9000
 9001 	continue
	! if (k_mt.eq.1) then
	! 	print*,'(COARSE GRID) Max sum of any equation: tmp_lmax = ', tmp_lmax
	! 	print*,'(COARSE GRID) Max sum of ALL equation: (i,j, sum, app(i,j) )  =', &
	! 		imax, jmax, tmp_max, app_mt(imax,jmax)-vp_mt(imax,jmax)*vkz_mt
	! end if
	return
end






subroutine relax(p_mt,rhs_mt,res_mt,nnx,nny,k_mt,apw_mt,ape_mt,aps_mt, &
		apn_mt,apws_mt,apwn_mt,apes_mt,apen_mt,app_mt,vp_mt,aedgex, &
		omega_mt,vkz_mt,p_mthat,a_mt,b_mt,c_mt,d_mt,f_mt,g_mt,iprob_mt,npre_mt)  

	implicit real (a-h,o-z)
	dimension rhs_mt(0:nnx+1,0:nny+1),p_mt(0:nnx+1,0:nny+1),res_mt(0:nnx+1,0:nny+1)
			dimension apw_mt(nnx,nny),ape_mt(nnx,nny),aps_mt(nnx,nny), &
			apn_mt(nnx,nny),apws_mt(nnx,nny),apwn_mt(nnx,nny), &
			apes_mt(nnx,nny),apen_mt(nnx,nny),app_mt(nnx,nny), &
			vp_mt(0:nnx+1,0:nny+1),aedgex(2)

!	'p_mthat' is the delta value equal to 'p_mt(n+1) - p_mt(n)'
	dimension p_mthat(0:nnx+1,0:nny+1) 
	dimension a_mt(nny),b_mt(nny),c_mt(nny),d_mt(nny),f_mt(nny),g_mt(nny)

	! tmp :  Adds the coefficients at every grid point (adds equations)
	dimension	tmp(0:nnx+1,0:nny+1)
	! tmp_l: Adds coefficients of a single equation
	real		tmp_l, tmp_lmax, tmp_max

	p_mthat=0.0
	res_mt=0.0

do 333 jiter=1,npre_mt

	resb_max=-1.e7
	tmp   = 0.;	tmp_max = 0.;	imax = 0;	jmax=0;
	tmp_l = 0.;	tmp_lmax= 0.;

	do i=1,nnx
		do j=1,nny

			if (k_mt.eq.1) then
				tmp_l =	  ape_mt(i,j)+ apw_mt (i,j)+ apn_mt(i,j)+ aps_mt (i,j)+ apws_mt(i,j)  &
					+ apwn_mt(i,j)+ apes_mt(i,j)+ apen_mt(i,j) +  (app_mt(i,j)-vp_mt(i,j)*vkz_mt);
				tmp_lmax = amax1(tmp_l,tmp_lmax);
				if (abs(tmp_l).gt.1e-12) then
					print*, 'tmp_l = a+b+c = ',tmp_l
					STOP " PROBLEM IN NEW MULTIGRID "
				end if
			end if

			if (k_mt.eq.1) then
				jloc = (i-1)*nny + i

				tmp( i , j ) =  tmp( i ,j) + app_mt(i,j)-vp_mt(i,j)*vkz_mt
				tmp(i-1, j ) =  tmp(i-1,j) + apw_mt (i,j)
				tmp(i+1, j ) =  tmp(i+1,j) + ape_mt (i,j)

				tmp(i-1,j+1) =  tmp(i-1,j+1) + apwn_mt(i,j)
				tmp( i ,j+1) =  tmp( i ,j+1) + apn_mt (i,j)
				tmp(i+1,j+1) =  tmp(i+1,j+1) + apen_mt(i,j)

				tmp(i-1,j-1) =  tmp(i-1,j-1) + apws_mt(i,j)
				tmp( i ,j-1) =  tmp( i ,j-1) + aps_mt (i,j)
				tmp(i+1,j-1) =  tmp(i+1,j-1) + apes_mt(i,j)
			end if


			res_mt(i,j)=-rhs_mt(i,j)+(ape_mt (i,j)*p_mt(i+1,j  )+apw_mt (i,j)*p_mt(i-1,j  ) &
						+ apn_mt (i,j)*p_mt(i  ,j+1)+aps_mt (i,j)*p_mt(i  ,j-1) &
						+ apws_mt(i,j)*p_mt(i-1,j-1)+apwn_mt(i,j)*p_mt(i-1,j+1) &
						+ apes_mt(i,j)*p_mt(i+1,j-1)+apen_mt(i,j)*p_mt(i+1,j+1) &
						+(app_mt (i,j)-vp_mt(i,j)*vkz_mt)*p_mt(i,j))
			d_mt(j)=omega_mt*(res_mt(i,j)+apws_mt(i,j)*p_mthat(i-1,j-1)+   &
				apw_mt(i,j)*p_mthat(i-1,j)+apwn_mt(i,j)*p_mthat(i-1,j+1))
			c_mt(j) = -apn_mt(i,j)
			b_mt(j) = -app_mt(i,j)+vp_mt(i,j)*vkz_mt
			a_mt(j) = -aps_mt(i,j)
			resb_max=amax1(abs(res_mt(i,j)),resb_max)
		end do
		!ccccccccccccccccccccccccccccccccccccccccc
		!       Thomas algorithm begin
		!ccccccccccccccccccccccccccccccccccccccccc
		a_mt(1) = 0.
		c_mt(nny) = 0.
		f_mt(1) = c_mt(1) / b_mt(1)
		g_mt(1) = d_mt(1) / b_mt(1)
		!ccccccccccccccccccccccccccccccc
		!     forward elimination
		!ccccccccccccccccccccccccccccccc
		do j = 2, nny 
			f_mt(j) = c_mt(j) /(b_mt(j) - a_mt(j)*f_mt(j-1)) 
			g_mt(j) =(d_mt(j)-a_mt(j)*g_mt(j-1))/(b_mt(j)-a_mt(j)*f_mt(j-1))
		end do
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
		!    back substitution (solution put into d(j))
		!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
		d_mt(nny) = g_mt(nny)
		do j = nny-1, 1, -1
			d_mt(j) = g_mt(j) - f_mt(j)*d_mt(j+1)
		end do
		!ccccccccccccccccccccccccccccccc
		!    Thomas algorithm end
		!ccccccccccccccccccccccccccccccc
		do j=1,nny
			p_mthat(i,j)=d_mt(j)
		end do
	end do

	!ccccccccccccccccccccccccccccccccccccccccccccccc
	! TJ: enforce x-periodicty on p_mthat here
	!	- (nnx.eq.ngxm), (nny.eq.ngym)
	!	- global indices => don't need MPI
	!ccccccccccccccccccccccccccccccccccccccccccccccc
	p_mthat(  0,  :) = 0.5*(p_mthat(0,:)   + p_mthat(nnx,:))
	p_mthat(nnx+1,:) = 0.5*(p_mthat(nnx+1,:) + p_mthat(1,:))

	if (k_mt.eq.1) then
		do ii=0,nnx+1
		do jj=0,nny+1
				if ( amax1( tmp_max,tmp(ii,jj) ) .gt. tmp_max) then
					tmp_max = amax1( tmp_max,tmp(ii,jj) )
					imax = ii
					jmax = jj
				end if
				if (abs(tmp(ii,jj)).gt.1e-12) then
					print*, 'tmp = a+b+c = ',ii,jj,tmp(ii,jj)
					STOP " PROBLEM IN NEW MULTIGRID "
				end if
		end do
		end do
	end if

	if(iprob_mt.eq.1)  print*,'jiter =',jiter,'resb_max= ',resb_max

	!ccccccccccccccccccccc
	!      update
	!ccccccccccccccccccccc
	do j=1,nny
	do i=1,nnx
		p_mt(i,j)=p_mt(i,j)+p_mthat(i,j)
	end do
	end do

	!ccccccccccccccccccccccccccccccccccccccccccccccc
	! TJ: enforce x-periodicty on p_mt here
	!	- (nnx.eq.ngxm), (nny.eq.ngym)
	!ccccccccccccccccccccccccccccccccccccccccccccccc
	p_mt(  0,  :) = p_mt(0,:)     + p_mthat(0,:)
	p_mt(nnx+1,:) = p_mt(nnx+1,:) + p_mthat(nnx+1,:)

 333  continue
	! if (k_mt.eq.1) then
	! 	print*,'Max sum of any equation: tmp_lmax = ', tmp_lmax
	! 	print*,'Max sum of ALL equation: (i,j, sum, app(i,j)  )  =', & 
	! 		imax, jmax, tmp_max, app_mt(imax,jmax)-vp_mt(imax,jmax)*vkz_mt
	! end if
      return
end





subroutine resid_mt(	res_mt,p_mt,rhs_mt,nnx,nny,apw_mt,ape_mt,aps_mt,apn_mt, &
			apws_mt,apwn_mt,apes_mt,apen_mt,app_mt,vp_mt,vkz_mt)
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! returns minus residual for the model problem.
	! input quantities are p_mt(1:nnx,1:nny) and rhs_mt(1:nnx,1:nny),
	! while res_mt(1:nnx,1:nny) is returned
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit real (a-h,o-z)
	dimension 	res_mt(0:nnx+1,0:nny+1),rhs_mt(0:nnx+1,0:nny+1), &
			p_mt(0:nnx+1,0:nny+1)
	dimension 	apw_mt(nnx,nny),ape_mt(nnx,nny),aps_mt(nnx,nny), &
			apn_mt(nnx,nny),apws_mt(nnx,nny),apwn_mt(nnx,nny), &
			apes_mt(nnx,nny),apen_mt(nnx,nny),app_mt(nnx,nny), &
			vp_mt(0:nnx+1,0:nny+1)

	do j=1,nny
	do i=1,nnx
		res_mt(i,j)=rhs_mt(i,j)-(ape_mt (i,j)*p_mt(i+1,j  )+apw_mt (i,j)*p_mt(i-1,j  ) &
					+apn_mt (i,j)*p_mt(i  ,j+1)+aps_mt (i,j)*p_mt(i  ,j-1) &
					+apws_mt(i,j)*p_mt(i-1,j-1)+apwn_mt(i,j)*p_mt(i-1,j+1) &
					+apes_mt(i,j)*p_mt(i+1,j-1)+apen_mt(i,j)*p_mt(i+1,j+1) &
					+(app_mt (i,j)-vp_mt(i,j)*vkz_mt)*p_mt(i,j))
	end do
	end do

	!ccccccccccccccccccccccc
	!    boundary points
	!ccccccccccccccccccccccc
	res_mt(:,0    )=0.0
	res_mt(:,nny+1)=0.0
	res_mt(0    ,:)=0.0
	res_mt(nnx+1,:)=0.0

	return
end




subroutine restrict(pc_mt,pf_mt,ncx,ncy)
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! restriction :  determine coarse grid rhs using
	!                fine grid residual values
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

	!ccccccccccccccccccccccc
	!    boundary points
	!ccccccccccccccccccccccc
	pc_mt(:,0    )=0.0
	pc_mt(:,ncy+1)=0.0
	pc_mt(0    ,:)=0.0
	pc_mt(ncx+1,:)=0.0

	return
end



subroutine fill0(p_mt,nnx,nny)
	!ccccccccccccccccccccccccccccccc
	!   fill p_mt with zeros
	!ccccccccccccccccccccccccccccccc
	implicit real (a-h,o-z)
	dimension p_mt(0:nnx+1,0:nny+1)

	p_mt = 0.0
	return
end




subroutine addint(pf_mt,pc_mt,res_mt,nfx,nfy,aedgex,p_mthat,iprob_mt)
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! (*) coarse-to-fine interpolation, and adds res_mt to pf_mt
	! (*) nf is the fine-grid dimensions
	! (*) The coarse-grid solution is input as pc_mt(1:ncx,1:ncy),
	! where nc =  nf/2 + 1
	! (*) The fine grid solution is returned in pf_mt(1:nfx,1:nfy)
	! (*) res_mt(1:nfx,1:nfy) is used for temporary storage
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit real (a-h,o-z)
	dimension	res_mt(0:nfx+1,0:nfy+1),p_mthat(0:nfx+1,0:nfy+1),  & 
    			pc_mt(0:nfx/2+1,0:nfy/2+1),pf_mt(0:nfx+1,0:nfy+1),aedgex(2) 

	res_mt=0.0
	p_mthat=0.0
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! (*) coarse-to-fine prolongation by bilinear interpolation
	! (*) nf is the fine grid dimension, the coarse-grid solution
	! is input as pc_mt(0:ncx+1,0:ncy+1), where nc = nf/2 
	! (*) The fine-grid solution is returned in pf_mt(0:nfx+1,0:nfy+1)
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	ncx=nfx/2
	ncy=nfy/2
	!cccccccccccccccccccccccccccccccccccc
	!    interpolating veritcally
	!cccccccccccccccccccccccccccccccccccc
	do ic=0,ncx+1   
	do jc=1,ncy+1 
		jf1=2*jc-1 
		jf2=2*jc-2 
		p_mthat(ic,jf1)=0.25*pc_mt(ic,jc-1)+0.75*pc_mt(ic,jc)
		p_mthat(ic,jf2)=0.75*pc_mt(ic,jc-1)+0.25*pc_mt(ic,jc)
	end do
	end do
	!ccccccccccccccccccccccccccccccc
	! interpolating horizontally
	!ccccccccccccccccccccccccccccccc
	do jf=0,nfy+1
	do ic=1,ncx+1
		if1=2*ic-1
		if2=2*ic-2
		res_mt(if1,jf)=0.25*p_mthat(ic-1,jf)+0.75*p_mthat(ic,jf)
		res_mt(if2,jf)=0.75*p_mthat(ic-1,jf)+0.25*p_mthat(ic,jf)
	end do
	end do
	!cccccccccccccccccccccccccc
	!     boundary points
	!cccccccccccccccccccccccccc
	if(iprob_mt.ne.4) then
		do if=0,nfx+1
			!ccccccccccccccccccccccccccccccc
			!   ghost points are zero
			!ccccccccccccccccccccccccccccccc
			res_mt(if,0    )=0.0
			res_mt(if,nfy+1)=0.0
			!ccccccccccccccccccccccccccccccc
			!        neumann b.c.
			!ccccccccccccccccccccccccccccccc
  			res_mt(if,1  )=res_mt(if,2    )
  			res_mt(if,nfy)=res_mt(if,nfy-1)
		end do

		do jf=0,nfy+1
			!ccccccccccccccccccccccccccccccc
			!     ghost points are zero
			!ccccccccccccccccccccccccccccccc
  			res_mt(0    ,jf)=0.0
  			res_mt(nfx+1,jf)=0.0
			!ccccccccccccccccccccccccccccccc
			!        neumann b.c.
			!ccccccccccccccccccccccccccccccc
  			res_mt(1  ,jf)=res_mt(2    ,jf)
  			res_mt(nfx,jf)=res_mt(nfx-1,jf)
		end do
	else
		do if=0,nfx+1
			!ccccccccccccccccccccccccccccccc
			!    ghost points are zero
			!ccccccccccccccccccccccccccccccc
			res_mt(if,0    )=0.0
			res_mt(if,nfy+1)=0.0
			!ccccccccccccccccccccccccccccccc
			!        neumann b.c.
			!ccccccccccccccccccccccccccccccc
			res_mt(if,1  )=res_mt(if,2    )
			res_mt(if,nfy)=res_mt(if,nfy-1)
		end do

		do jf=0,nfy+1
			!ccccccccccccccccccccccccccccccc
			!    ghost points are zero
			!ccccccccccccccccccccccccccccccc
			res_mt(0    ,jf)=0.0
			res_mt(nfx+1,jf)=0.0
			!ccccccccccccccccccccccccccccccc
			!        neumann b.c.
			!ccccccccccccccccccccccccccccccc
			res_mt(1  ,jf)=res_mt(2    ,jf)
			res_mt(nfx,jf)=res_mt(nfx-1,jf)
		end do
	end if

	!ccccccccccccccccccccccccccccccccccccccccccccccc
	! TJ: enforce x-periodicty on res_mt 
	!	- (nnx.eq.ngxm), (nny.eq.ngym)
	!	- global indices => don't need MPI
	!ccccccccccccccccccccccccccccccccccccccccccccccc
	res_mt(0,:)     = res_mt(nfx,:)
	res_mt(nfx+1,:) = res_mt(1,:)

	do j=0,nfy+1 
	do i=0,nfx+1 
		pf_mt(i,j)=pf_mt(i,j)+res_mt(i,j)
	end do
	end do
	p_mthat=0.0
	res_mt=0.0
	return
end


