!=======================================================================
! Calculating the outflow boundary condition
!
! outflow_convective(deltaT, Uout)
! outflow_massFlowConsistency(Uin, Uout, Vin, Vout)
! massFlowRate_X(rmdot, U1, i)*
! massFlowRate_Y(rmdot, V1, j)*
!

module outflow
	use data_export
	use mesh_export
end module

!=======================================================================
subroutine outflow_convective(deltaT, Uout)
	use outflow
	real Uout(0:nly+1, 0:ngz+1, 3)
	real dUout(0:nly+1, 0:ngz+1, 3)
	real ubulk				! bulk veclocity for laminar flow

	ia = nlxb
        !SY---need to be fixed for general case!!
        !deltaXi = 1./(x(imax+1) - x(imax))
        !deltaXMi = 1./(xm(imax) - xm(imax-1))
        deltaXMi = 1./(xpg(ngx,1) - xpg(ngx-1,1))

	! Use the so-called "convective condition"
	call readFloat("ubulk", ubulk)
	u_a = ubulk

	! Take u_a to be the maximum velocity at the exit plane
	u_a = 0.
	if (iblock .eq. npx) u_a = rmaxval(uc(:,ia-1:ia,:), 2*(nly+1)*(ngz+1))
	call globalMax(u_a, 1)
	!--- print*,'I am jblock =',jblock, 'and my convective u =', u_a

	! Make sure Courant number is less than 1 for stability
	if (u_a*deltaT*deltaXMi > 1.) u_a = 1./(deltaT*deltaXMi)

	! Convective condition: u,t = -c u,x
	if (iblock .eq. npx) then
                !------ uc ------
                do j=1,nlyb-1
                        jj = jbmap_1(j)
                        uc(1:ngz-1, nlxb, j)=(1.0-u_a*deltaT/(xpg(ngx,jj)-xpg(ngx-1,jj))) *uc(1:ngz-1, nlxb  ,j) &
                                                 +u_a*deltaT/(xpg(ngx,jj)-xpg(ngx-1,jj))  *uc(1:ngz-1, nlxb-1,j)
                end do
                uc(1:ngz-1,nlxb+1,0:nly)=2.0*uc(1:ngz-1,nlxb,0:nly)-uc(1:ngz-1,nlxb-1,0:nly)

                !------ vc ------
                do j=1,nlyb
                        jj = jbmap_1(j)
                        vc(1:ngz-1, nlxb, j)=(1.0-u_a*deltaT/(xpg(ngx,jj)-xpg(ngx-1,jj))) *vc(1:ngz-1, nlxb  ,j) &
                                                 +u_a*deltaT/(xpg(ngx,jj)-xpg(ngx-1,jj))  *vc(1:ngz-1, nlxb-1,j)
                end do

                !------ wc ------
                do j=1,nlyb-1
                        jj = jbmap_1(j)
                        wc(1:ngz-1, nlxb, j)=(1.0-u_a*deltaT/(xpg(ngx,jj)-xpg(ngx-1,jj))) *wc(1:ngz-1, nlxb  ,j) &
                                                 +u_a*deltaT/(xpg(ngx,jj)-xpg(ngx-1,jj))  *wc(1:ngz-1, nlxb-1,j)
                enddo
		if (jblock.eq.1)  &
			wc(1:ngz-1, nlxb, 0  )=2.0*wc(1:ngz-1, nlxb-1, 0  )-wc(1:ngz-1, nlxb-2, 0  )
		if (jblock.eq.npy)  &
			wc(1:ngz-1, nlxb,nlyb)=2.0*wc(1:ngz-1, nlxb-1,nlyb)-wc(1:ngz-1, nlxb-2,nlyb)

	end if

	return
end

!=======================================================================
subroutine outflow_massFlowConsistency(Uin, Uout, Vin, Vout)
	use outflow
	parameter (i1=i1_, i2=i2_, j1=j1_, j2=j2_, k1=k1_, k2=k2_)
	real Uin(0:nly+1, 0:ngz+1,3), Uout(0:nly+1, 0:ngz+1,3)
	real Vin(0:nlx+1, 0:ngz+1,3), Vout(0:nlx+1, 0:ngz+1,3)

	! Calculate mass flows entering and leaving the domain
	call massFlowRate_X(flowIn_X , Uin(0,0,1) ,  1 )
	call massFlowRate_X(flowOut_X, Uout(0,0,1), ngx)
	call massFlowRate_Y(flowIn_Y , Vin(0,0,2) ,  1 )
	call massFlowRate_Y(flowOut_Y, Vout(0,0,2), ngy)

	! Watch out for bad outflow condition
	! Uout, at least, must be an outflow
	if (flowOut_X <= spacing(1.)) then
		! Adjust Uout to fix things
		Uout(:,:,1) = 1.
		call massFlowRate_X(flowOut_X, Uout(0,0,1), ngx)
	end if

	! Uin is assumed to always be an inflow
	flowIn = flowIn_X &
	         +0.5*(flowIn_Y + abs(flowIn_Y)) &
	         -0.5*(flowOut_Y - abs(flowOut_Y))
	flowOut = flowOut_X &
	          -0.5*(flowIn_Y - abs(flowIn_Y)) &
	          +0.5*(flowOut_Y + abs(flowOut_Y))

	! Scale the outflow velocities to balance 
	! the inflow rate and net volume change rate
	ratio = flowIn / flowOut
	Uout(:,:,1) = ratio*Uout(:,:,1)
	if (flowIn_Y < 0.) Vin(:,:,2) = ratio*Vin(:,:,2)
	if (flowOut_Y > 0.) Vout(:,:,2) = ratio*Vout(:,:,2)

	return
end

!=======================================================================
! Integrate the velocity profile over a y-z plane
!
subroutine massFlowRate_X(rmdot, U1, i)
	use outflow
	real U1(0:nly+1, 0:ngz+1)

	ia = ibmin ; if (iblock ==  1 ) ia = imino
	ib = ibmax ; if (iblock == npx) ib = imaxo

	rmdot = 0.
	if (ia <= i .and. ib >= i) then
		do k=1,ngz-1
		do j=j1_T,j2_T
			dA = dy(j)*dz
			rmdot = rmdot + U1(j,k)*dA
		end do
		end do
	end if
	call globalSum(rmdot, 1)

	return
end

!=======================================================================
! Integrate the velocity profile over an x-z plane
!
subroutine massFlowRate_Y(rmdot, V1, j)
	use outflow
	real V1(0:nlx+1, 0:ngz+1)

	ja = jbmin ; if (jblock ==  1 ) ja = jmino
	jb = jbmax ; if (jblock == npy) jb = jmaxo

	rmdot = 0.
	if (ja <= j .and. jb >= j) then
		do k=1,ngz-1
		do i=i1_T,i2_T
			dA = dx*dz
			rmdot = rmdot + V1(i,k)*dA
		end do
		end do
	end if
	call globalSum(rmdot, 1)

	return
end

