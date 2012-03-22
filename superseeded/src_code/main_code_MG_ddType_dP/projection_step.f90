
subroutine Update_Flux()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Update the velocity field, p is actually dt*p
! Equation (3.220) --------> look at note in margin next to (3.220)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	use data_export

	before = realClock()

	ka = 1
	kb = ngz-1

	do jj=j1_u, j2_u
	do ii=i1_u, i2_u
		i  = ibmap_1(ii)
		j  = jbmap_1(jj)

		!-----------------------------------------------------------------------
		! Construct flags to tell begining and end of local domain
		!-----------------------------------------------------------------------
		! jflag_start = [ 1,0,0,0,......,0,0,0 ]	( 1<= jj <= nlyb-1 )
		! jflag_end   = [ 0,0,0,0,......,0,0,1 ]	( 1<= jj <= nlyb-1 )
		!-----------------------------------------------------------------------
		jflag_start = 1/jj		! Only equal to one when (jj=1)
		jflag_end   = jj/(nlyb-1)	! Only equal to one when (jj=nlyb-1)


		!----------------------------Replace if statement with flags---------------------------
		! if(jj.ne.1.and.jj.ne.nlyb-1) then
		! 	u(ka:kb,ii,jj)= ust(ka:kb,ii,jj)-((suxix(i,j)*0.5*(suxix(i,j)+suxix(i+1,j)) &
                !      			+suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i+1,j)))*p(ka:kb,ii,jj) &
                !     			-(suxix(i,j)*0.5*(suxix(i,j)+suxix(i-1,j)) &
                !      			+suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i-1,j)))*p(ka:kb,ii-1,jj)) &
                !     			/((vp(i,j)+vp(i-1,j))/2.) &
                !    			-((suxix(i,j)*(svetax(i,j+1)+svetax(i-1,j+1))/2. &
                !      			+suxiy(i,j)*(svetay(i,j+1)+svetay(i-1,j+1))/2.) &
                !     			*(p(ka:kb,ii,jj)+p(ka:kb,ii-1,jj)          &
                !       			+p(ka:kb,ii,jj+1)+p(ka:kb,ii-1,jj+1))/4. &
                !     			-(suxix(i,j)*(svetax(i,j)+svetax(i-1,j))/2. &
                !      			+suxiy(i,j)*(svetay(i,j)+svetay(i-1,j))/2.) &
                !     			*(p(ka:kb,ii,jj)+p(ka:kb,ii-1,jj)  &
                !       			+p(ka:kb,ii,jj-1)+p(ka:kb,ii-1,jj-1))/4.) &
                !     			/((vp(i,j)+vp(i-1,j))/2.)
		! else if(jj.eq.nlyb-1) then    
		! 	u(ka:kb,ii,jj)=ust(ka:kb,ii,jj)-((suxix(i,j)*0.5*(suxix(i,j)+suxix(i+1,j)) &
                !      			+suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i+1,j)))*p(ka:kb,ii,jj) &
                !     			-(suxix(i,j)*0.5*(suxix(i,j)+suxix(i-1,j)) &
                !      			+suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i-1,j)))*p(ka:kb,ii-1,jj)) &
                !     			/((vp(i,j)+vp(i-1,j))/2.)
		! else if(jj.eq.1) then  
		! 	u(ka:kb,ii,jj)=ust(ka:kb,ii,jj)-((suxix(i,j)*0.5*(suxix(i,j)+suxix(i+1,j)) &
                !      			+suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i+1,j)))*p(ka:kb,ii,jj) &
                !     			-(suxix(i,j)*0.5*(suxix(i,j)+suxix(i-1,j)) &
                !      			+suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i-1,j)))*p(ka:kb,ii-1,jj)) &
                !     			/((vp(i,j)+vp(i-1,j))/2.)
		! end if
                u(ka:kb,ii,jj)= (1-jflag_start-jflag_end) * (    ust(ka:kb,ii,jj) -           &
                                                                ( ( suxix(i,j)*0.5*(suxix(i,j)+suxix(i+1,j))                      &
                                                                   +suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i+1,j)) )*p(ka:kb,ii  ,jj)   &
                                                                 -( suxix(i,j)*0.5*(suxix(i,j)+suxix(i-1,j))                      &
                                                                   +suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i-1,j)) )*p(ka:kb,ii-1,jj) ) &
                                                                /( (vp(i,j)+vp(i-1,j))/2. ) - &
                                                                ( ( suxix(i,j)*(svetax(i,j+1)+svetax(i-1,j+1))/2.                 &
                                                                   +suxiy(i,j)*(svetay(i,j+1)+svetay(i-1,j+1))/2. )               &
                                                                  *( p(ka:kb,ii,jj  )+p(ka:kb,ii-1,jj  )                          &
                                                                    +p(ka:kb,ii,jj+1)+p(ka:kb,ii-1,jj+1) )/4.                     &
                                                                 -( suxix(i,j)*(svetax(i,j  )+svetax(i-1,j  ))/2.                 &
                                                                   +suxiy(i,j)*(svetay(i,j  )+svetay(i-1,j  ))/2. )               &
                                                                  *( p(ka:kb,ii,jj  )+p(ka:kb,ii-1,jj  )                          &
                                                                    +p(ka:kb,ii,jj-1)+p(ka:kb,ii-1,jj-1) )/4.         )           &
                                                                /( (vp(i,j)+vp(i-1,j))/2. )                                       )

                u(ka:kb,ii,jj)=u(ka:kb,ii,jj) + jflag_end * (    ust(ka:kb,ii,jj) -           &
                                                                ( ( suxix(i,j)*0.5*(suxix(i,j)+suxix(i+1,j))                      &
                                                                   +suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i+1,j)) )*p(ka:kb,ii  ,jj)   &
                                                                 -( suxix(i,j)*0.5*(suxix(i,j)+suxix(i-1,j))                      &
                                                                   +suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i-1,j)) )*p(ka:kb,ii-1,jj) ) &
                                                                /( (vp(i,j)+vp(i-1,j))/2. )                                       )

                u(ka:kb,ii,jj)=u(ka:kb,ii,jj) + jflag_start * (  ust(ka:kb,ii,jj) -           &
                                                                ( ( suxix(i,j)*0.5*(suxix(i,j)+suxix(i+1,j))                      &
                                                                   +suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i+1,j)) )*p(ka:kb,ii  ,jj)   &
                                                                - ( suxix(i,j)*0.5*(suxix(i,j)+suxix(i-1,j))                      &
                                                                   +suxiy(i,j)*0.5*(suxiy(i,j)+suxiy(i-1,j)) )*p(ka:kb,ii-1,jj) ) &
                                                                /( (vp(i,j)+vp(i-1,j))/2. )                                       )

	end do
	end do


	after = realClock()
	cpu_updateu= (after - before) - t_over
	before = realClock()

	do jj=j1_v,j2_v
	do ii=i1_v,i2_v
		i = ibmap_1(ii)
		j = jbmap_1(jj)

		!-----------------------------------------------------------------------
		! Construct flags to tell begining and end of local domain
		!-----------------------------------------------------------------------
		! iflag_start = [ 1,0,0,0,......,0,0,0 ]	( 1<= ii <= nlxb-1 )
		! iflag_end   = [ 0,0,0,0,......,0,0,1 ]	( 1<= ii <= nlxb-1 )
		!-----------------------------------------------------------------------
		iflag_start = 1/ii		! Only equal to one when (ii=1)
		iflag_end   = ii/(nlxb-1)	! Only equal to one when (ii=nlxb-1)


		!----------------------------Replace if statement with flags---------------------------
		! if(ii.ne.1 .and. ii.ne.nlxb-1) then 
		! 	v(ka:kb,ii,jj)=vst(ka:kb,ii,jj)-((svetax(i,j)*(svetax(i,j)+svetax(i,j+1))/2. &
                !      			+svetay(i,j)*(svetay(i,j)+svetay(i,j+1))/2.)*p(ka:kb,ii,jj) &
                !     			-(svetax(i,j)*(svetax(i,j)+svetax(i,j-1))/2. &
                !      			+svetay(i,j)*(svetay(i,j)+svetay(i,j-1))/2.)*p(ka:kb,ii,jj-1)) &
                !     			/((vp(i,j)+vp(i,j-1))/2.) &
                !    			-((svetax(i,j)*(suxix(i+1,j)+suxix(i+1,j-1))/2. &
                !      			+svetay(i,j)*(suxiy(i+1,j)+suxiy(i+1,j-1))/2.) &
                !     			*(p(ka:kb,ii,jj)+p(ka:kb,ii+1,jj)  &
                !       			+p(ka:kb,ii,jj-1)+p(ka:kb,ii+1,jj-1))/4. &
                !     			-(svetax(i,j)*(suxix(i,j)+suxix(i,j-1))/2. &
                !      			+svetay(i,j)*(suxiy(i,j)+suxiy(i,j-1))/2.) &
                !     			*(p(ka:kb,ii,jj)+p(ka:kb,ii,jj-1)  &
                !       			+p(ka:kb,ii-1,jj)+p(ka:kb,ii-1,jj-1))/4.) &
                !     			/((vp(i,j)+vp(i,j-1))/2.)
		! else if (ii.eq.nlxb-1) then 
		! 	v(ka:kb,ii,jj)=vst(ka:kb,ii,jj)-((svetax(i,j)*(svetax(i,j)+svetax(i,j+1))/2. &
                !      			+svetay(i,j)*(svetay(i,j)+svetay(i,j+1))/2.)*p(ka:kb,ii,jj) &
                !     			-(svetax(i,j)*(svetax(i,j)+svetax(i,j-1))/2. &
                !      			+svetay(i,j)*(svetay(i,j)+svetay(i,j-1))/2.)*p(ka:kb,ii,jj-1)) &
                !     			/((vp(i,j)+vp(i,j-1))/2.)
		! else if(ii.eq.1) then 
		! 	v(ka:kb,ii,jj)=vst(ka:kb,ii,jj)-((svetax(i,j)*(svetax(i,j)+svetax(i,j+1))/2. & 
                !      			+svetay(i,j)*(svetay(i,j)+svetay(i,j+1))/2.)*p(ka:kb,ii,jj) & 
                !     			-(svetax(i,j)*(svetax(i,j)+svetax(i,j-1))/2. &
                !      			+svetay(i,j)*(svetay(i,j)+svetay(i,j-1))/2.)*p(ka:kb,ii,jj-1)) &
                !     			/((vp(i,j)+vp(i,j-1))/2.)
		! end if

                v(ka:kb,ii,jj)=(1-iflag_start-iflag_end) * (     vst(ka:kb,ii,jj) -           &
                                                                ( ( svetax(i,j)*(svetax(i,j)+svetax(i,j+1))/2.                      &
                                                                   +svetay(i,j)*(svetay(i,j)+svetay(i,j+1))/2. )*p(ka:kb,ii,jj  )   &
                                                                 -( svetax(i,j)*(svetax(i,j)+svetax(i,j-1))/2.                      &
                                                                   +svetay(i,j)*(svetay(i,j)+svetay(i,j-1))/2. )*p(ka:kb,ii,jj-1) ) &
                                                                /( (vp(i,j)+vp(i,j-1))/2. ) - &
                                                                ( ( svetax(i,j)*(suxix(i+1,j)+suxix(i+1,j-1))/2.                    &
                                                                   +svetay(i,j)*(suxiy(i+1,j)+suxiy(i+1,j-1))/2. )                  &
                                                                  *( p(ka:kb,ii,jj  )+p(ka:kb,ii+1,jj  )                            &
                                                                    +p(ka:kb,ii,jj-1)+p(ka:kb,ii+1,jj-1) )/4.                       &
                                                                 -( svetax(i,j)*(suxix(i  ,j)+suxix(i  ,j-1))/2.                    &
                                                                   +svetay(i,j)*(suxiy(i  ,j)+suxiy(i  ,j-1))/2. )                  &
                                                                  *( p(ka:kb,ii  ,jj)+p(ka:kb,ii  ,jj-1)                            &
                                                                    +p(ka:kb,ii-1,jj)+p(ka:kb,ii-1,jj-1) )/4.       )               &
                                                                /( (vp(i,j)+vp(i,j-1))/2. )                                         )

                v(ka:kb,ii,jj)=v(ka:kb,ii,jj) + iflag_end * (    vst(ka:kb,ii,jj) -            &
                                                                ( ( svetax(i,j)*(svetax(i,j)+svetax(i,j+1))/2.                      &
                                                                   +svetay(i,j)*(svetay(i,j)+svetay(i,j+1))/2. )*p(ka:kb,ii,jj  )   &
                                                                 -( svetax(i,j)*(svetax(i,j)+svetax(i,j-1))/2.                      &
                                                                   +svetay(i,j)*(svetay(i,j)+svetay(i,j-1))/2. )*p(ka:kb,ii,jj-1) ) &
                                                                /( (vp(i,j)+vp(i,j-1))/2. )                                         )

                v(ka:kb,ii,jj)=v(ka:kb,ii,jj) + iflag_start * (  vst(ka:kb,ii,jj) -            &
                                                                ( ( svetax(i,j)*(svetax(i,j)+svetax(i,j+1))/2.                      &
                                                                   +svetay(i,j)*(svetay(i,j)+svetay(i,j+1))/2. )*p(ka:kb,ii,jj  )   &
                                                                - ( svetax(i,j)*(svetax(i,j)+svetax(i,j-1))/2.                      &
                                                                   +svetay(i,j)*(svetay(i,j)+svetay(i,j-1))/2. )*p(ka:kb,ii,jj-1) ) &
                                                                /( (vp(i,j)+vp(i,j-1))/2. )                                         )

	end do
	end do

	after = realClock()
	cpu_updatev= (after - before) - t_over
	before = realClock()

	do jj=j1_w,j2_w
	do ii=i1_w,i2_w
		i = ibmap_1(ii)
		j = jbmap_1(jj)
                w(ka:kb,ii,jj)=wst(ka:kb,ii,jj)-( swz(i,j)*swz(i,j)*p(ka:kb    ,ii,jj) &
                                                 -swz(i,j)*swz(i,j)*p(ka-1:kb-1,ii,jj) )/vp(i,j)
	end do
	end do

	after = realClock()
	cpu_updatew= (after - before) - t_over
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! in z direction
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	u(0  ,:,:)=u(ngz-1,:,:)
	u(ngz,:,:)=u(1    ,:,:)

	v(0  ,:,:)=v(ngz-1,:,:)
	v(ngz,:,:)=v(1    ,:,:)

	w(0    ,:,:)=w(ngz-1,:,:)
	w(ngz  ,:,:)=w(1    ,:,:)
	w(ngz+1,:,:)=w(2    ,:,:)

	! if(ntime-ntime_.eq.nsteps) then			!(nsteps) is not defines and hence assumed (zero)
	! 	if (irank.eq.1) then
	! 		write(6,*) 'cpu_updateu = ', cpu_updateu
	! 		write(6,*) 'cpu_updatev = ', cpu_updatev 
	! 		write(6,*) 'cpu_updatew = ', cpu_updatew 
	! 	end if
	! end if

	!T  Exchange halos in x-direction , then in y-direction
	!T-- x-direction
	call updateBorder_lim(u, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 3)
	call updateBorder_lim(v, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 3)
	call updateBorder_lim(w, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 1)
	!T-- y-direction
	call updateBorder_lim(u, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 3)
	call updateBorder_lim(v, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 3)
	call updateBorder_lim(w, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 1)

	return
end

!========================================================================

subroutine Divergence_Check()
!ccccccccccccccccccccccccccccc
!  Compute the divergence
!ccccccccccccccccccccccccccccc
	use data_export

	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! find the maximum divergence
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	divmax=-1.e6
	do k=1,ngz-1
	do j=j1_T,j2_T
	do i=i1_T,i2_T
                ii = ibmap_1(i)
                jj = jbmap_1(j)

		! divmax=max(divmax, &
             	! 	abs(u(k,i+1,j)-u(k,i,j)+v(k,i,j+1)-v(k,i,j)+w(k+1,i,j)-w(k,i,j))  /  &
		! 	( abs(u(k,i+1,j))+abs(u(k,i,j))+abs(v(k,i,j+1))+abs(v(k,i,j))+abs(w(k+1,i,j))+abs(w(k,i,j)) )    )
		divmax=max(divmax, &
             		abs(u(k,i+1,j)-u(k,i,j)+v(k,i,j+1)-v(k,i,j)+w(k+1,i,j)-w(k,i,j))/vp(ii,jj)  )
		if(abs(w(k,i,j)).le.1.e-12) w(k,i,j)=0.0
	end do
	end do
	end do

	call globalMax(divmax, 1)

	if(mod(ntime,ipout).eq.0) then
	if (irank.eq.1)  write(6,3000) divmax
 3000     format(1x,' Maximum divergence is ', 1pe12.5)
	end if

	return
end

!========================================================================

subroutine Cart_to_Flux()
	use data_export
	real*8 :: vtemp(1:ngz-1), utemp(1:ngz-1)

	!cccccccccccccccccccccccccc
	!    Calculate U^{xi}
	!cccccccccccccccccccccccccc
	do j=0,nlyb
	do i=1,nlxb
		ii = ibmap_1(i)
		jj = jbmap_1(j)
		a1=suxix(ii,jj)
		a2=suxiy(ii,jj)
		vtemp(:)=( vc(1:ngz-1,i-1,j+1)+vc(1:ngz-1,i,j+1)   &
          		+vc(1:ngz-1,i-1,j  )+vc(1:ngz-1,i,j  ) )/4.
		u(1:ngz-1,i,j)=a1*uc(1:ngz-1,i,j)+a2*vtemp(:)
	end do
	end do

	if (iblock.eq.1) then
		u(1:ngz-1, 0, :)=2.0*u(1:ngz-1, 1,:)-u(1:ngz-1, 2 ,:)
	end if
	if (iblock.eq.npx) then
		u(1:ngz-1, nlxb+1,:)=2.0*u(1:ngz-1, nlxb,:)-u(1:ngz-1, nlxb-1,:)
	end if

	!cccccccccccccccccccccccccc
	!     Calculate U^{eta}
	!cccccccccccccccccccccccccc 
	do j=1,nlyb
	do i=0,nlxb
		ii = ibmap_1(i) 
		jj = jbmap_1(j)
		a1=svetax(ii,jj)
		a2=svetay(ii,jj)  
		utemp(:)=( uc(1:ngz-1,i+1,j-1)+uc(1:ngz-1,i+1,j)   &
          		+uc(1:ngz-1,i,j-1)  +uc(1:ngz-1,i,j)  )/4.
		v(1:ngz-1,i,j)=a1*utemp(:)+a2*vc(1:ngz-1,i,j)
	end do
	end do
      
	if (jblock.eq.1) then 
		v(1:ngz-1,:,0)=2.0*v(1:ngz-1,:,1)-v(1:ngz-1,:,2)
	end if
	if (jblock.eq.npy) then
		v(1:ngz-1,:,nlyb+1)=2.0*v(1:ngz-1,:,nlyb)-v(1:ngz-1,:,nlyb-1)
	end if    

	!ccccccccccccccccccccccccccccccc
	!	Calculate U^{z}        
	!ccccccccccccccccccccccccccccccc
	do j=0,nlyb   
	do i=0,nlxb
		ii = ibmap_1(i)
		jj = jbmap_1(j)
		!ccccccccccccccccccccccccccccccc
		! Compute S^{z} at w locations
		!ccccccccccccccccccccccccccccccc
		a3=swz(ii,jj)
		w(1:ngz-1,i,j)=a3*wc(1:ngz-1,i,j)
	end do
	end do

	!Use Periodicity Instead	w(0    ,:,:)=2.0*w(1  ,:,:)-w(2    ,:,:)
	!Use Periodicity Instead	w(ngz+1,:,:)=2.0*w(ngz,:,:)-w(ngz-1,:,:)	! w(ngz,:,:) not defined yet

	!ccccccccccccccccccccccccccccccc
	!	In z direction
	!ccccccccccccccccccccccccccccccc
	u(0  ,:,:)=u(ngz-1,:,:)
	u(ngz,:,:)=u(1    ,:,:)

	v(0  ,:,:)=v(ngz-1,:,:)
	v(ngz,:,:)=v(1    ,:,:)

	w(0    ,:,:)=w(ngz-1,:,:)
	w(ngz  ,:,:)=w(1    ,:,:)
	w(ngz+1,:,:)=w(2    ,:,:)
 
	return
end


!========================================================================

subroutine Flux_to_Cart
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Transform flux variables to Cartesian components
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	use data_export
	real*8 :: vtemp(1:ngz-1), utemp(1:ngz-1)

	!ccccccccccccccccccccccccccccccc
	! 	Calculate uc
	!ccccccccccccccccccccccccccccccc
	do ii=2,nlxb-1
	do jj=1,nlyb-1  
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		!--------------------------------------------------
		!    calculate S_{xi},S_{eta} at u locations
		!--------------------------------------------------
		a1=surxix(i,j)
		a2=surxiy(i,j)
		b1=(svretax(i-1,j)+svretax(i,j)+svretax(i,j+1)+svretax(i-1,j+1))/4.
		b2=(svretay(i-1,j)+svretay(i,j)+svretay(i,j+1)+svretay(i-1,j+1))/4.
		!-----------------------------------------------------------
		!    interpolate U^{eta} to get U^{eta} at u locations
		!-----------------------------------------------------------
		vtemp(1:ngz-1)=(v(1:ngz-1, ii-1, jj)+v(1:ngz-1, ii-1, jj+1)    &
               			+v(1:ngz-1, ii  , jj)+v(1:ngz-1, ii  , jj+1))/4.
		uc(1:ngz-1, ii,jj)=a1*u(1:ngz-1,ii,jj)+b1*vtemp(1:ngz-1)
	end do
	end do

	!ccccccccccccccccccccccccccccccc
	! 	Calculate vc
	!ccccccccccccccccccccccccccccccc
	do ii=1,nlxb-1
	do jj=1,nlyb
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		!-------------------------------------------------
		! calculate S_{xi},S_{eta} at v locations
		!-------------------------------------------------
		b1=svretax(i,j)
		b2=svretay(i,j)
		a1=(surxix(i,j-1)+surxix(i,j)+surxix(i+1,j)+surxix(i+1,j-1))/4.
		a2=(surxiy(i,j-1)+surxiy(i,j)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.
		!---------------------------------------------------------
		!   interpolate U^{xi} to get U^{xi} at v locations
		!---------------------------------------------------------
		utemp(1:ngz-1)=(u(1:ngz-1, ii  , jj-1)+u(1:ngz-1, ii  , jj)    &
               			+u(1:ngz-1, ii+1, jj-1)+u(1:ngz-1, ii+1, jj))/4.
		vc(1:ngz-1, ii, jj)=a2*utemp(1:ngz-1)+b2*v(1:ngz-1, ii, jj)
	end do
	end do
      
	!ccccccccccccccccccccccccccccccc
	!	Calculate wc
	!ccccccccccccccccccccccccccccccc
	do ii=1,nlxb-1
	do jj=1,nlyb-1
		i = ibmap_1(ii)
		j = jbmap_1(jj)
		!-------------------------
		!    calculate S_{z}
		!-------------------------
		a3=swrz(i,j)
		wc(1:ngz-1, ii, jj)=a3*w(1:ngz-1, ii, jj)
	end do
	end do

	!ccccccccccccccccccccccccccccccc
	!	In z direction
	!ccccccccccccccccccccccccccccccc
	uc(0  , :,:)=uc(ngz-1, :, :)
	uc(ngz, :,:)=uc(1    , :, :)

	vc(0  , :,:)=vc(ngz-1, :, :)
	vc(ngz, :,:)=vc(1    , :, :)

	wc(0    , :, :)=wc(ngz-1, :, :)
	wc(ngz  , :, :)=wc(1    , :, :)
	wc(ngz+1, :, :)=wc(2    , :, :)

	return
end



