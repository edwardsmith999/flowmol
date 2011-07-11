
subroutine relax_b(    p_mt,     rhs_mt,     res_mt,    nnx,     nny,     k_mt,       &
	 	     apw_mt,     ape_mt,     aps_mt, apn_mt, apws_mt,  apwn_mt,       &
	            apes_mt,    apen_mt,     app_mt,  vp_mt,  aedgex, omega_mt,       &
   	             vkz_mt,       a_mt,       b_mt,   c_mt,    d_mt,     f_mt, g_mt, &
		   iprob_mt, resminb_mt, nitmaxb_mt)

	use data_export
	implicit real (a-h,o-z)

	dimension 	rhs_mt(0:nnx+1,0:nny+1), p_mt(0:nnx+1,0:nny+1), res_mt(0:nnx+1,0:nny+1)
	dimension	 vp_mt(0:nnx+1,0:nny+1)

	dimension	apw_mt(nnx,nny),  ape_mt(nnx,nny), &
			aps_mt(nnx,nny),  apn_mt(nnx,nny), &
		       apws_mt(nnx,nny), apwn_mt(nnx,nny), &
		       apes_mt(nnx,nny), apen_mt(nnx,nny), &
                        app_mt(nnx,nny)
	dimension 	aedgex(2)

	dimension 	a_mt(nny), b_mt(nny), c_mt(nny), d_mt(nny), f_mt(nny), g_mt(nny)
	dimension       p_mthat(0:nnx+1,0:nny+1)

	p_mthat=0.0
	res_mt=0.0
 993    format(1x,10(1pe10.3))
 331    continue

	icounter = 0
	resb_max=1.e7

	do while(resb_max.gt.resminb_mt.and.icounter.lt.nitmaxb_mt)
        icounter = icounter+1
	resb_max=-1.e7

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!SY  Two-color ordering (in the x direction) Gauss-Seidel method
	!    with line relaxation in the y direction
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	do ii=0,1
		ia = 2**ii
		ib = nnx-mod(ia,2)
		do i=ia,ib,2
			do j=1,nny
				res_mt(i,j) =-rhs_mt(i,j)+( ape_mt (i,j)*p_mt(i+1,j  )+apw_mt (i,j)*p_mt(i-1,j  ) &
							  + apn_mt (i,j)*p_mt(i  ,j+1)+aps_mt (i,j)*p_mt(i  ,j-1) &
							  + apws_mt(i,j)*p_mt(i-1,j-1)+apwn_mt(i,j)*p_mt(i-1,j+1) &
							  + apes_mt(i,j)*p_mt(i+1,j-1)+apen_mt(i,j)*p_mt(i+1,j+1) &
							  + (app_mt(i,j)-vp_mt(i,j)*vkz_mt        )*p_mt(i  ,j  ) )
				d_mt(j)=omega_mt*(res_mt(i,j)  &
					+ ii * (apws_mt(i,j)*p_mthat(i-1,j-1) + apw_mt (i,j)*p_mthat(i-1,j  ) &
					      + apwn_mt(i,j)*p_mthat(i-1,j+1) + apes_mt(i,j)*p_mthat(i+1,j-1) &
					      + ape_mt (i,j)*p_mthat(i+1,j  ) + apen_mt(i,j)*p_mthat(i+1,j+1) ) )
				c_mt(j)=-apn_mt(i,j)
				b_mt(j)=-app_mt(i,j)+vp_mt(i,j)*vkz_mt
				a_mt(j)=-aps_mt(i,j)
				resb_max=amax1(abs(res_mt(i,j)),resb_max)

			end do   ! j

			call triDiagonal_(id_y, a_mt, b_mt, c_mt, d_mt, nny)

			do j=1,nny
			   p_mthat(i,j)=d_mt(j)
			end do
		end do    ! i

		!-----------  Pass Zebra lines (left, right) ==> (0,1)  --------------
		if (ii.eq.0) then
	                call updateBorder_2D(p_mthat, nnx, nny , id_y, 2, 1, nny, 1)
			call poisson_zebra  (p_mthat( 1 ,0:nny+1), p_mthat(nnx+1,0:nny+1), nny+2, id_x, ii )
		end if
		if (ii.eq.1) then
	                call updateBorder_2D(p_mthat, nnx, nny , id_y, 2, 1, nny, 1)
			call poisson_zebra  (p_mthat(nnx,0:nny+1), p_mthat(  0  ,0:nny+1), nny+2, id_x, ii )
		end if

	end do  ! ii

	!ccccccccccccccc
	!    update
	!ccccccccccccccc
	ia = 0
	ib = nnx+1
	ja = 0
	jb = nny+1

	if(iblock==1  ) ia = 1
	if(iblock==npx) ib = nnx
	if(jblock==1  ) ja = 1
	if(jblock==npy) jb = nny

	do j=ja,jb
	do i=ia,ib
		p_mt(i,j)=p_mt(i,j)+p_mthat(i,j)
	end do
	end do

	call globalMax(resb_max, 1)

	end do ! do while

        if(irank.eq.1.and.k_mt.eq.1) write(6,991) icounter,resb_max

991    format(1x, 'The number of iterations at level   1 = ', 2x, I5, 2x,  &
              'and maximum residual at this level = ', 2x, 1pe12.5)

	return
end

subroutine relax(     p_mt,  rhs_mt,    res_mt,       nnx,     nny,     k_mt,       &
	       	    apw_mt,  ape_mt,    aps_mt,    apn_mt, apws_mt,  apwn_mt,       &
	           apes_mt, apen_mt,    app_mt,     vp_mt,  aedgex, omega_mt,       &
	            vkz_mt,    a_mt,      b_mt,      c_mt,    d_mt,     f_mt, g_mt, &
	          iprob_mt, npre_mt, resmin_mt, nitmax_mt,       m )

	use data_export
	implicit real (a-h,o-z)

	dimension 	rhs_mt(0:nnx+1,0:nny+1), p_mt(0:nnx+1,0:nny+1), res_mt(0:nnx+1,0:nny+1)
	dimension	 vp_mt(0:nnx+1,0:nny+1)

	dimension	apw_mt(nnx,nny),  ape_mt(nnx,nny), &
			aps_mt(nnx,nny),  apn_mt(nnx,nny), &
		       apws_mt(nnx,nny), apwn_mt(nnx,nny), &
		       apes_mt(nnx,nny), apen_mt(nnx,nny), &
                        app_mt(nnx,nny)
	dimension 	aedgex(2)
	dimension 	a_mt(nny), b_mt(nny), c_mt(nny), d_mt(nny), f_mt(nny), g_mt(nny)

!	'p_mthat' is the delta value equal to 'p_mt(n+1) - p_mt(n)'
	dimension      p_mthat(0:nnx+1,0:nny+1)

	p_mthat=0.0
	res_mt=0.0

	res_max=-1.e7

	do 333 jiter=1,npre_mt

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!SY  Two-color ordering (in the x direction) Gauss-Seidel method
	!    with line relaxation in the y direction
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		do ii=0,1
			ia = 2**ii
			ib = nnx-mod(ia,2)
			do i=ia,ib,2
				do j=1,nny
				res_mt(i,j) =-rhs_mt(i,j)+( ape_mt (i,j)*p_mt(i+1,j  )+apw_mt (i,j)*p_mt(i-1,j  ) &
							  + apn_mt (i,j)*p_mt(i  ,j+1)+aps_mt (i,j)*p_mt(i  ,j-1) &
							  + apws_mt(i,j)*p_mt(i-1,j-1)+apwn_mt(i,j)*p_mt(i-1,j+1) &
							  + apes_mt(i,j)*p_mt(i+1,j-1)+apen_mt(i,j)*p_mt(i+1,j+1) &
							  + (app_mt(i,j)-vp_mt(i,j)*vkz_mt        )*p_mt(i  ,j  ) )
				d_mt(j)=omega_mt*(res_mt(i,j)  &
					+ ii * (apws_mt(i,j)*p_mthat(i-1,j-1) + apw_mt (i,j)*p_mthat(i-1,j  ) &
					      + apwn_mt(i,j)*p_mthat(i-1,j+1) + apes_mt(i,j)*p_mthat(i+1,j-1) &
					      + ape_mt (i,j)*p_mthat(i+1,j  ) + apen_mt(i,j)*p_mthat(i+1,j+1) ) )
				c_mt(j)=-apn_mt(i,j)
				b_mt(j)=-app_mt(i,j)+vp_mt(i,j)*vkz_mt
				a_mt(j)=-aps_mt(i,j)
				res_max=amax1(abs(res_mt(i,j)),res_max)

				end do   ! j

				call triDiagonal_(id_y, a_mt, b_mt, c_mt, d_mt, nny)

            			do j=1,nny
					p_mthat(i,j)=d_mt(j)
	    			end do
			end do  ! i

			!-----------  Pass Zebra lines (left, right) ==> (0,1)  --------------
			if (ii.eq.0) then
	                        call updateBorder_2D(p_mthat, nnx, nny , id_y, 2, 1, nny, 1)
				call poisson_zebra  (p_mthat( 1 ,0:nny+1), p_mthat(nnx+1,0:nny+1), nny+2, id_x, ii )
			end if
			if (ii.eq.1) then	
	                        call updateBorder_2D(p_mthat, nnx, nny , id_y, 2, 1, nny, 1)
				call poisson_zebra  (p_mthat(nnx,0:nny+1), p_mthat(  0  ,0:nny+1), nny+2, id_x, ii )
			end if

		end do ! ii

	!ccccccccccccccccccccc
	!       update
	!ccccccccccccccccccccc
	ia = 0
	ib = nnx+1
	ja = 0
	jb = nny+1

	if(iblock==1  ) ia = 1
	if(iblock==npx) ib = nnx
	if(jblock==1  ) ja = 1
	if(jblock==npy) jb = nny

	do j=ja,jb
	do i=ia,ib
		p_mt(i,j)=p_mt(i,j)+p_mthat(i,j)
	end do
	end do

 333    continue  ! jiter

	call globalMax(res_max, 1)

        if(irank.eq.1.and.k_mt.eq.1) write(6,771) m,npre_mt,res_max

771    format(1x, 'The number of iterations at level', 1x, I3, 1x, '= ', 2x, I5, 2x,  &
              'and maximum residual at this level = ', 2x, 1pe12.5)

	return
end

subroutine resid_mt( res_mt,    p_mt,  rhs_mt,     nnx,    nny, apw_mt, ape_mt, aps_mt, apn_mt, &
 		    apws_mt, apwn_mt, apes_mt, apen_mt, app_mt,  vp_mt, vkz_mt)

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! returns minus residual for the model problem.
	! input quantities are p_mt(1:nnx,1:nny) and rhs_mt(1:nnx,1:nny),
	! while res_mt(1:nnx,1:nny) is returned
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit real (a-h,o-z)

	dimension 	res_mt(0:nnx+1,0:nny+1), rhs_mt(0:nnx+1,0:nny+1), &
			p_mt  (0:nnx+1,0:nny+1),  vp_mt(0:nnx+1,0:nny+1)

	dimension 	apw_mt (nnx,nny),ape_mt (nnx,nny),aps_mt (nnx,nny), &
			apn_mt (nnx,nny),apws_mt(nnx,nny),apwn_mt(nnx,nny), &
			apes_mt(nnx,nny),apen_mt(nnx,nny),app_mt (nnx,nny)

	do j=1,nny
	do i=1,nnx
		res_mt(i,j)=rhs_mt(i,j)- ( &
		                         ape_mt (i,j)*p_mt(i+1,j  )+apw_mt (i,j)*p_mt(i-1,j  ) &
					+apn_mt (i,j)*p_mt(i  ,j+1)+aps_mt (i,j)*p_mt(i  ,j-1) &
					+apws_mt(i,j)*p_mt(i-1,j-1)+apwn_mt(i,j)*p_mt(i-1,j+1) &
					+apes_mt(i,j)*p_mt(i+1,j-1)+apen_mt(i,j)*p_mt(i+1,j+1) &
					+(app_mt (i,j)-vp_mt(i,j)*vkz_mt       )*p_mt(i  ,j  ) )
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

subroutine restrict(pc_mt, pf_mt, ncx, ncy)

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! restriction :  determine coarse grid rhs using
	!                fine grid residual values
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

	!ccccccccccccccccccccccc
	!    boundary points
	!ccccccccccccccccccccccc
	pc_mt(:,0    )=0.0
	pc_mt(:,ncy+1)=0.0
	pc_mt(0    ,:)=0.0
	pc_mt(ncx+1,:)=0.0

	return
end

subroutine addint(pf_mt, pc_mt, res_mt, nfx, nfy, aedgex, iprob_mt)

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! (*) coarse-to-fine interpolation, and adds res_mt to pf_mt
	! (*) nf is the fine-grid dimensions
	! (*) The coarse-grid solution is input as pc_mt(1:ncx,1:ncy),
	! where nc =  nf/2 + 1
	! (*) The fine grid solution is returned in pf_mt(1:nfx,1:nfy)
	! (*) res_mt(1:nfx,1:nfy) is used for temporary storage
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	use data_export
	implicit real (a-h,o-z)

	dimension  res_mt(0:nfx+1  ,0:nfy+1  ),  &
		    pc_mt(0:nfx/2+1,0:nfy/2+1), pf_mt(0:nfx+1,0:nfy+1), aedgex(2)

	dimension  p_mthat(0:nfx/2+1,0:nfy+1)

	res_mt=0.0
	p_mthat=0.0
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! (*) coarse-to-fine prolongation by bilinear interpolation
	! (*) nf is the fine grid dimension, the coarse-grid solution
	! is input as pc_mt(0:ncx+1,0:ncy+1), where nc = nf/2
	! (*) The fine-grid solution is returned in pf_mt(0:nfx+1,0:nfy+1)
	!     NOTE: (p_mthat), (res_mt) are temporary storage only.
	!            They are set to zero at the end.
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	ncx=nfx/2
	ncy=nfy/2

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! interpolating veritcally
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
		!ccccccccccccccccccccccccccccccc
		!  in the y direction
		!ccccccccccccccccccccccccccccccc
		if (jblock.eq.1) then
			do ifa=0,nfx+1
				!ccccccccccccccccccccccccccccccc
				!   ghost points are zero
				!ccccccccccccccccccccccccccccccc
				res_mt(ifa,0  )=0.0
				!ccccccccccccccccccccccccccccccc
				!        neumann b.c.
				!ccccccccccccccccccccccccccccccc
				res_mt(ifa,1  )=res_mt(ifa,2    )
			end do
		end if
		if (jblock.eq.npy) then
			do ifa=0,nfx+1
				!ccccccccccccccccccccccccccccccc
				!   ghost points are zero
				!ccccccccccccccccccccccccccccccc
				res_mt(ifa,nfy+1)=0.0
				!ccccccccccccccccccccccccccccccc
				!        neumann b.c.
				!ccccccccccccccccccccccccccccccc
				res_mt(ifa,nfy  )=res_mt(ifa,nfy-1)
			end do
		end if
		!ccccccccccccccccccccccccccccccc
		!  in the x direction
		!ccccccccccccccccccccccccccccccc
		if (iblock.eq.1) then
			do jf=0,nfy+1
				!ccccccccccccccccccccccccccccccc
				!   ghost points are zero
				!ccccccccccccccccccccccccccccccc
				res_mt(0,jf   )=0.0
				!ccccccccccccccccccccccccccccccc
				!        neumann b.c.
				!ccccccccccccccccccccccccccccccc
				res_mt(1,jf   )=res_mt(2,jf     )
			end do
		end if
		if (iblock.eq.npx) then
			do jf=0,nfy+1
				!ccccccccccccccccccccccccccccccc
				!   ghost points are zero
				!ccccccccccccccccccccccccccccccc
				res_mt(nfx+1,jf)=0.0
				!ccccccccccccccccccccccccccccccc
				!        neumann b.c.
				!ccccccccccccccccccccccccccccccc
				res_mt(nfx  ,jf)=res_mt(nfx-1,jf)
			end do
		end if

	else

		!ccccccccccccccccccccccccccccccc
		!  in the y direction
		!ccccccccccccccccccccccccccccccc
		if (jblock.eq.1) then
			do ifa=0,nfx+1
				!ccccccccccccccccccccccccccccccc
				!   ghost points are zero
				!ccccccccccccccccccccccccccccccc
				res_mt(ifa,0  )=0.0
				!ccccccccccccccccccccccccccccccc
				!        neumann b.c.
				!ccccccccccccccccccccccccccccccc
				res_mt(ifa,1  )=res_mt(ifa,2    )
			end do
		end if
		if (jblock.eq.npy) then
			do ifa=0,nfx+1
				!ccccccccccccccccccccccccccccccc
				!   ghost points are zero
				!ccccccccccccccccccccccccccccccc
				res_mt(ifa,nfy+1)=0.0
				!ccccccccccccccccccccccccccccccc
				!        neumann b.c.
				!ccccccccccccccccccccccccccccccc
				res_mt(ifa,nfy  )=res_mt(ifa,nfy-1)
			end do
		end if
		!ccccccccccccccccccccccccccccccc
		!  in the x direction
		!ccccccccccccccccccccccccccccccc

		if (iblock.eq.1) then
			do jf=0,nfy+1
				!ccccccccccccccccccccccccccccccc
				!   ghost points are zero
				!ccccccccccccccccccccccccccccccc
				res_mt(0,jf   )=0.0
				!ccccccccccccccccccccccccccccccc
				!        neumann b.c.
				!ccccccccccccccccccccccccccccccc
				res_mt(1,jf   )=res_mt(2,jf     )
			end do
		end if
		if (iblock.eq.npx) then
			do jf=0,nfy+1
				!ccccccccccccccccccccccccccccccc
				!   ghost points are zero
				!ccccccccccccccccccccccccccccccc
				res_mt(nfx+1,jf)=0.0
				!ccccccccccccccccccccccccccccccc
				!        neumann b.c.
				!ccccccccccccccccccccccccccccccc
				res_mt(nfx  ,jf)=res_mt(nfx-1,jf)
			end do
		end if

	end if

	do j=0,nfy+1
	do i=0,nfx+1
		pf_mt(i,j)=pf_mt(i,j)+res_mt(i,j)
	end do
	end do

	res_mt=0.0

	return
end



subroutine find_max_res(res_mt, res_mtmax)
	!--------------------------------------------
	! Find the maximum residual
	!--------------------------------------------
	use data_export
	dimension res_mt(0:nixp+1,0:niyp+1)
	real      res_mtmax

	res_mtmax=-1.e7

	do i=1,nixp
	do j=1,niyp
		res_mtmax=max(abs(res_mt(i,j)), res_mtmax)
	end do
	end do

	call globalMax(res_mtmax, 1)
	return
end
