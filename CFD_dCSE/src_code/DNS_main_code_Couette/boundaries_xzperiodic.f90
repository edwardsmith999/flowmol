!=======================================================================
! Boundary conditions for inflow-outflow boundary layer flow
!
! boundaries_init()
! boundaries_read()
! boundaries_write()
!
! CartesianBC(deltaT)
! FluxBC()
! poissonCoeffBC(ay, by, cy)
! poissonMatrixBC(A, B, C)
! poissonRhsBC()
!
! timeDependent_Inlet_BC(deltaT)
! boundaries_inflowOK(deltaT, iok)
! laminarInflowBC()
! mean_Time(deltaT)

module boundaries
	use data_export
	use mesh_export

	integer inflowType		! type of inflow condition
	integer inflowFluct		! type of inflow fluctuations
	real    ubulk			! bulk velocity for laminar inflow

        real Uin(0:nly+1,0:ngz+1,3), Uout(0:nly+1,0:ngz+1,3) !u,v,w left-right
        real Vin(0:nlx+1,0:ngz+1,3), Vout(0:nlx+1,0:ngz+1,3) !u,v,w top-bottom
        real Win(0:nlx+1,0:nly+1,3), Wout(0:nlx+1,0:nly+1,3) !u,v,w front-back

        real UatLeftBC(0:ngy+1,3)       !Left boundary inflow velocity profile
        real UppatLeftBC(0:ngy+1)       !Left boundary inflow velocity profile

	integer jref
	real    pref

end module
!========================================================================

subroutine boundaries_init()

!     The inlet flow is composed of a base flow and a perturbation.
!     Parameters are read from the archive

      use boundaries
      
	! Root process is assumed to be on inflow boundary
	call data_isRoot(iflag)
	if (iflag == 1) then
		if (iblock .ne. 1) &
			stop "Root process needs to be on inflow boundary"
	end if

	! Set up inflow and outflow BCs (initialize to zero)
        Uin  = 0.0
        Uout = 0.0
        Vin  = 0.0
        Vout = 0.0
        Win  = 0.0
        Wout = 0.0

        ! Read BCs for each face
        call readInt("BC_bottom", BC_bottom)
        call readInt("BC_top"   , BC_top   )
        call readInt("BC_left"  , BC_left  )
        call readInt("BC_right" , BC_right )
        call readInt("BC_front" , BC_front )
        call readInt("BC_back"  , BC_back  )

        ! Verify for periodic boundary conditions

        if ( (BC_bottom.eq.2 .and. BC_top.ne.2) .or. &
             (BC_top.eq.2 .and. BC_bottom.ne.2)) then
                stop "Both top and bottom boundaries must be periodic"
        end if

        ! Verify for inflow-outflow boundary conditions

        if ( (BC_left.eq.1 .and. BC_right.ne.1) .or. &
             (BC_right.eq.1 .and. BC_left.ne.1)) then
                stop "Left and right boundaries must be inflow/outlfow, respectively"
        end if

        ! Verify for periodic boundary conditions

        if ( (BC_left.eq.2 .and. BC_right.ne.2) .or. &
             (BC_right.eq.2 .and. BC_left.ne.2)) then
                stop "Both left and right boundaries must be periodic"
        end if

        ! Verify for periodic boundary conditions

        if ( (BC_front.eq.2 .and. BC_back.ne.2) .or. &
             (BC_back.eq.2 .and. BC_front.ne.2)) then
                stop "Both front and back boundaries must be periodic"
        end if


	!Determine fluctuations to superimpose onto base flow
      	!call readInt("inflowFluct", inflowFluct)
      	!call fluct_init(Uin, UinBlasius, UinppBlasius)

	! Setup the reference pressure
	! Should be located in a quiescent region of the flow
	! Reference pressure on upper boundary j=jmax, set to zero
	jref = jmax
	pref = 0.

	return
end


subroutine boundaries_read()
	use boundaries

	if (iblock.eq.1) then
			select case (inflowType)
			case (1)
				call inflow_read()
			end select
	end if

	return
end

subroutine boundaries_write()
	use boundaries

        call writeInt("BC_bottom", BC_bottom)
        call writeInt("BC_top"   , BC_top   )
        call writeInt("BC_left"  , BC_left  ) 
        call writeInt("BC_right" , BC_right )
        call writeInt("BC_front" , BC_front ) 
        call writeInt("BC_back"  , BC_back  )

	
	if (iblock == 1) then
			select case (inflowType)
			case (0)
				call writeFloat("ubulk", ubulk)
			case (1)
				call inflow_write()
			end select

			! Write information about inlet fluctuations
			call fluct_write()
	end if

	return
end
!===============================================================================
subroutine CartesianBC(deltaT) 
      	use boundaries
     
      	!===================================================================
      	! 	NO SLIP B.C. 
      	!===================================================================
	!--------------------------------------------------------------
	!	Bottom B.C. 
	!--------------------------------------------------------------
	if (jblock.eq.1) then
		
		uc(:, :, 0) = -uc(:, :, 1) 
		vc(:, :, 1) =  0.0            
		wc(:, :, 0) = -wc(:, :, 1)
	
		!Extend
		vc(:, :, 0) = vc(:, :, 2)        
        	
	end if

	!--------------------------------------------------------------
	!	 Top B.C.
	!-------------------------------------------------------------- 
	if (jblock.eq.npy) then
        	
        	uc(:, :, nlyb  ) =  1.0 + (0.5*(ypg(1,ngy)-ypg(1,ngy-1))) * ( (1.0 - 0.0) / (ypg(1,ngy)-ypg(1,1)) )
		vc(:, :, nlyb  ) =  0.0
		wc(:, :, nlyb  ) = -wc(:, :, nlyb-1)
         		
		!Extend
		vc(:, :, nlyb+1) = vc(:, :, nlyb-1) 
        	
	end if

	!-----------------------------------------------------------------------
        !       Periodic in streamwise (x) direction
        !-----------------------------------------------------------------------
	if (iblock.eq.1 .or. iblock.eq.npx) then
		if(npx.eq.1) then 
		
			uleftbc_temp(:,:,1) = uc(:,1,:)
                        uleftbc_temp(:,:,2) = uc(:,2,:)
		       	vleftbc_temp(:,:,1) = vc(:,1,:)
                        wleftbc_temp(:,:,1) = wc(:,1,:)
	
			urightbc_temp(:,:,1) = uc(:,nlxb-1,:)
			urightbc_temp(:,:,2) = uc(:, nlxb, :)
		       	vrightbc_temp(:,:,1) = vc(:,nlxb-1,:)
                        wrightbc_temp(:,:,1) = wc(:,nlxb-1,:)
	
		else
			! ==> (iblock.eq.1  ) receives ucbce
			! ==> (iblock.eq.npx) receives ucbcw
			call updateBorder_xperiodic(1)
		end if
	end if

        if (iblock.eq.1) then
		uc(:,0,:) = urightbc_temp(:,:,1)			   
		uc(:,1,:) = urightbc_temp(:,:,2)
		vc(:,0,:) = vrightbc_temp(:,:,1)			   
                wc(:,0,:) = wrightbc_temp(:,:,1)           	  
	end if

        if (iblock.eq.npx) then
		
		vc(:, nlxb, :) = vleftbc_temp(:,:,1)		   
                wc(:, nlxb, :) = wleftbc_temp(:,:,1)		 
 	
		!Extend
		uc(:,nlxb+1,:) = uleftbc_temp(:,:,2)		  
 
        end if

        !-----------------------------------------------------------------------
        !       Periodic in spanwise (z) direction
        !-----------------------------------------------------------------------
        uc(0    , :, :) = uc(ngz-1, :, :)
        uc(ngz  , :, :) = uc(1    , :, :)

        vc(0    , :, :) = vc(ngz-1, :, :)
        vc(ngz  , :, :) = vc(1    , :, :)

        wc(0    , :, :) = wc(ngz-1, :, :)
        wc(ngz  , :, :) = wc(1    , :, :)
        !DANGER ZONE: wc(  1  , :,:) = 0.5*(wc(1, :,:)+wc(ngz, :,:)) (wc only goes from 1:ngz-1 so wc(ngz,:,:) does not exist)
	!Extend
        wc(ngz+1, :, :) = wc(2    , :, :)

      	!----------------- Block boundaries -------------------
        !-- Exchange halos in x-direction , then in y-direction
        
	!-- x-direction
	call updateBorder_lim(uc, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 3)
        call updateBorder_lim(vc, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 3)
        call updateBorder_lim(wc, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 2)

        !-- y-direction
        call updateBorder_lim(uc, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 3)
        call updateBorder_lim(vc, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 3)
        call updateBorder_lim(wc, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 2)
      	
	return
end
!============================================================================================

subroutine FluxBC

!     Determine the BCs in terms of flux.     

      	use boundaries
      	!real*8  :: vtempb(ngz+1), vtempe(ngz+1)
      	!real*8  :: utempb(ngz+1), utempe(ngz+1)
      	real*8  :: vtempb(ngz-1), vtempe(ngz-1)
      	real*8  :: utempb(ngz-1), utempe(ngz-1)
      	real*8, dimension(2) :: Atmp

	!=======================================================================
        !       NO SLIP B.C.  
        !=======================================================================
	!-----------------------------------------------------------------------
      	!	Bottom B.C.     
      	!-----------------------------------------------------------------------
	if (jblock.eq.1) then
         	
       		u(:, :, 0) = -u(:, :, 1)
       		v(:, :, 1) =  0.0
       		w(:, :, 0) = -w(:, :, 1)

       	      	!Extend
		v(:, :, 0) = v(:, :, 2)	
         	
	end if

	! -----------------------------------------------------------------------
        !       Top B.C.     
        !-----------------------------------------------------------------------
	if (jblock.eq.npy) then
		!----------------------------------------
		! This section (Sec C) assumes knowing
		! Cartesian BC and calculating Flux BC.
		! (Use this section with Dirichlet top BC)
		!----------------------------------------
		!----------  U^{eta} (i.e. v) -----------
	      	v(:, :, nlyb  ) =  0.0
		!Extend
		v(:, :, nlyb+1) = v(:, :, nlyb-1)	

		!----------  U^{z} (i.e. w) -----------
	     	w(:, :, nlyb  ) = -w(:, :, nlyb-1)

		!----------  U^{xi} (i.e. u) -----------
		do i=1,nlxb
			ii = ibmap_1(i)
			vtempe(:)= ( vc(1:ngz-1, i-1,nlyb+1)+vc(1:ngz-1, i,nlyb+1)  &
					+vc(1:ngz-1, i-1,nlyb  )+vc(1:ngz-1, i,nlyb  ))/4.
			u(1:ngz-1,i,nlyb)=suxix(ii,ngy)*uc(1:ngz-1,i,nlyb)+suxiy(ii,ngy)*vtempe(:)
		end do
		!----------- End of Sec C ---------------
        end if

	!-----------------------------------------------------------------------
        !	Periodic in streamwise (x) direction     
        !-----------------------------------------------------------------------
        if (iblock.eq.1 .or. iblock.eq.npx) then
                if (npx.eq.1) then

                        uleftbc_temp(:,:,1) = u(:,1,:)
                        uleftbc_temp(:,:,2) = u(:,2,:)
                        vleftbc_temp(:,:,1) = v(:,1,:)
                        wleftbc_temp(:,:,1) = w(:,1,:)

                        urightbc_temp(:,:,1) = u(:,nlxb-1,:)
                        urightbc_temp(:,:,2) = u(:, nlxb, :)
                        vrightbc_temp(:,:,1) = v(:,nlxb-1,:)
                        wrightbc_temp(:,:,1) = w(:,nlxb-1,:)

                else
                        call updateBorder_xperiodic(2)
                end if

        end if

        if (iblock.eq.1) then

                u(:,0,:) = urightbc_temp(:,:,1)
                u(:,1,:) = urightbc_temp(:,:,2)
                v(:,0,:) = vrightbc_temp(:,:,1)
                w(:,0,:) = wrightbc_temp(:,:,1)

        end if

        if (iblock.eq.npx) then

                v(:, nlxb, :) = vleftbc_temp(:,:,1)
                w(:, nlxb, :) = wleftbc_temp(:,:,1)

                !Extend
                u(:,nlxb+1,:) = uleftbc_temp(:,:,2)

        end if

	!-----------------------------------------------------------------
        !	Periodic in spanwise (z) direction
        !-----------------------------------------------------------------
        u( 0 , :,:) = u(ngz-1, :,:)
        u(ngz, :,:) = u(  1  , :,:)
        
	v( 0 , :,:)=v(ngz-1, :,:)
        v(ngz, :,:)=v(  1  , :,:)

        w(  0  , :,:) = w(ngz-1, :,:)
        !TAZ DANGER ZONE: w(  1  , :,:) = 0.5*(w(1, :,:)+w(ngz, :,:))
        w(ngz  , :,:) = w(1, :,:)
	!Extend
        w(ngz+1, :,:) = w(2, :,:)

        !=======================================================================
        ! sum of mass flux along boundaries is zero
        !=======================================================================
        bmassin =0.0
        bmassout=0.0
        if (jblock.eq.1  ) then
        do i=i1_T,i2_T
           bmassin  = bmassin  + sum(v(1:ngzm, i, 1   )) !TJ-- equal to zero?
        end do
        end if

        if (jblock.eq.npy) then
        do i=i1_T,i2_T
           bmassin  = bmassin  - sum(v(1:ngzm, i, nlyb)) !TJ-- equal to zero?
        end do
        end if

        bmass = bmassin
        call globalSum(bmass,1)

        if (iblock.eq.1  ) then
        do j=j1_T,j2_T
           bmassin  = bmassin  + sum(u(1:ngzm, 1   , j))
        end do
        end if

        if (iblock.eq.npx) then
        do j=j1_T,j2_T
           bmassout = bmassout + sum(u(1:ngzm, nlxb, j))
        end do
        end if

        do 710 j=j1_T,j2_T
        do 710 i=i1_T,i2_T
           bmassin  = bmassin  + w(1,i,j)-w(ngz,i,j) !TJ-- equal to zero?
   710  continue 

        Atmp(1) = bmassin
        Atmp(2) = bmassout
        call globalSum(Atmp, 2)
        bmassin  = Atmp(1)
        bmassout = Atmp(2)

        if (irank.eq.1)   then
                write(6,*) 'Mass inflow  =  ', bmassin
                write(6,*) 'Mass outflow =  ', bmassout
        end if

        !if(abs(bmassout).ge.1.e-10) then 
        !   if (iblock.eq.npx) then
        !       u(1:ngzm,nlxb,1:nly-1)=u(1:ngzm,nlxb,1:nly-1)*bmassin/bmassout
        !   end if
        !end if

        bmassout=0.0
        if (iblock.eq.npx) then
           do j=j1_T,j2_T
              bmassout = bmassout + sum(u(1:ngzm, nlxb, j))
           end do
        end if

        call globalSum(bmassout, 1)
      
	if (irank.eq.1) then
                write(6,*) 'After adjustment'
                write(6,*) 'Mass inflow  =  ', bmassin
                write(6,*) 'Mass outflow =  ', bmassout
                write(6,*) 'mass corrected: inflow - outflow =  ',bmassout-bmassin
        end if

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! TJ: this is incorrect for streamwise periodicity
	!
	!if (iblock.eq.npx) then
        !	!Extend
        !	u(:, nlxb+1, 0:nly) = 2.0*u(:, nlxb, 0:nly)-u(:, nlxb-1, 0:nly)
        !end if
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! TJ: these boundary conditions have been enforced already
	!
        !TAZ DANGER ZONE:  Xiahua has u-periodicity here:
	!if (BC_front.eq.2) then
        !	u( 0 , :, :) = u(ngz-1, :, :)
        !	u(ngz, :, :) = u(1    , :, :)
        !end if
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	

!----------------- Block boundaries -------------------
      	!-- Exchange halos in x-direction , then in y-direction
      	!-- x-direction
	call updateBorder_lim(u, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 3)
      	call updateBorder_lim(v, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 3)
      	call updateBorder_lim(w, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 2)

      	!-- y-direction
      	call updateBorder_lim(u, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 3)
      	call updateBorder_lim(v, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 3)
      	call updateBorder_lim(w, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 2)

	 !print*, "u inside FluxBC"
         !do k=0,ngz
         !do j=0,nly
         	!print *, "u(", k, ",", nlxb+1, ",", j, ") :" ,u(k, nlxb+1, j)
         !end do
         !end do

return
end
!=======================================================================
subroutine poissonCoeffBC(ay, by, cy) 
        use boundaries
        real ay(ny), by(ny), cy(ny)

        ! Neumann condition on all boundaries

        ! Lower wall
        by(jmin) = by(jmin) + ay(jmin)
        ay(jmin) = 0.

        ! Upper wall
        by(jmax) = by(jmax) + cy(jmax)
        cy(jmax) = 0.

        return
end

subroutine poissonMatrixBC(A, B, C)
        use boundaries
        real A(nix+2,niz_2,niy_2), &
             B(nix+2,niz_2,niy_2), &
             C(nix+2,niz_2,niy_2)

        ! Fourier-space BC
        ! Special case of the zero wavenumber
        ka = max(kbmin_2(kblock_2), kmin) - kbmin_2(kblock_2) + 1
        kb = min(kbmax_2(kblock_2), kmin) - kbmin_2(kblock_2) + 1
        ja = max(jbmin_2(jblock_2), jref) - jbmin_2(jblock_2) + 1
        jb = min(jbmax_2(jblock_2), jref) - jbmin_2(jblock_2) + 1
        A(1, ka:kb, ja:jb) = 0.
        B(1, ka:kb, ja:jb) = 1.
        C(1, ka:kb, ja:jb) = 0.

        return
end

subroutine poissonRhsBC()
        use boundaries

	double precision	:: ust_left_halo(0:ngz,0:nly)
	double precision	:: ust_right_halo(0:ngz,0:nly)

        !T  Need a halo cell in {ust,vst,wst} before poisson solver
        !T-- x-direction
        call updateBorder_lim(ust, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 1)
        call updateBorder_lim(vst, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 1)
        call updateBorder_lim(wst, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 1)

        !T-- y-direction
        call updateBorder_lim(ust, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 1)
        call updateBorder_lim(vst, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 1)
        call updateBorder_lim(wst, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 1)

        !T-- Bottom Border
        if (jblock_1.eq.1) then
                ust(:,:,j1_u-1) = u(:,:,j1_u-1)
                vst(:,:,j1_v-1) = v(:,:,j1_v-1)
                wst(:,:,j1_u-1) = w(:,:,j1_u-1)
        end if

        !T-- Top Border
        if (jblock_1.eq.npy) then
                ust(:,:,j2_u+1) = u(:,:,j2_u+1)
                vst(:,:,j2_v+1) = v(:,:,j2_v+1)
                wst(:,:,j2_u+1) = w(:,:,j2_u+1)
        end if
	
	!ccccccccccccccccccccccccccccccccccccccccc
	! TJ: Periodicity in x
	!ccccccccccccccccccccccccccccccccccccccccc
	ust_left_halo  = 0.
	ust_right_halo = 0.

	if (iblock.eq.1 .or. iblock.eq.npx) then
		if (npx.eq.1) then
			ust_right_halo(:,:) = ust(:,  1 ,:)
			ust_left_halo(:,:)  = ust(:,nlxb,:)
		else
			if (iblock.eq.1) then
				ust_right_halo(:,:) = ust(:, 1, :)
				call poisson_xperiodic(ust_right_halo(0,0),ust_left_halo(0,0),((ngzm+2)*(nly+1)))
			else if (iblock.eq.npx) then
				ust_left_halo(:,:) = ust(:,nlxb,:)
				call poisson_xperiodic(ust_left_halo(0,0),ust_right_halo(0,0),((ngzm+2)*(nly+1)))
			end if
		end if
	end if

	if (iblock.eq.1) then
		ust(:,1,:) = ust_left_halo(:,:)
	end if

        !T-- Bottom Border
        !T-- Periodicity in Z-direction
        wst(ngz,:,:) = wst(1,:,:)

        return
end
!=======================================================================
subroutine timeDependent_Inlet_BC(deltaT)
	use boundaries
        
	!--------------------------------------------
	! 	Time-advance inflow BC
	!--------------------------------------------
	if (iblock == 1) then
		! Superimpose fluctuations on mean inlet
		!call mean_Time(deltaT)
		!call laminarInflowBC()
		!select case (inflowFluct)
		!case (0)
		!case (1:)
		!	call fluct_Time(deltaT)
		!end select

		select case (inflowFluct)
		case (0) ! No fluctuations
!		case (1) ! Random fluctuations
!			call fluct_random(Uin)
!		case (2) ! OS discrete mode fluctuations (in blayer)
!			call fluct_OSdiscrete(Uin)
!		case (3) ! OS continuous mode fluctuations (in freestream)
!			call fluct_OScontinuous(Uin)
!		case (3) ! OS continuous mode fluctuations (in freestream)
!			call fluct_OScontinuous(Uin)
!		case (4) ! OS discrete and continuous mode fluctuations
!			call fluct_freestream(Uin)
!		case (5:6) ! OS 2D mode + OS 3D oblique modes
!			call fluct_OS2D3D(Uin)
!		case (5:6) ! OS 2D & 3D modes
!			call fluct_OScontinuous(Uin)
		end select

	end if
	return
end

!=======================================================================
subroutine timeDependent_Wall_BC()
        use boundaries
        !--------------------------------------------
        !       Time-advance wall BC
        !--------------------------------------------
        if (jblock == 1) then
                Vin(:,:,1) = 0.0
                Vin(:,:,2) = 0.0
                Vin(:,:,3) = 0.0
!                Vin(:,:,3) = wbulk * cos(2.0*2.0*acos(0.0)*stime/wperiod)
        end if

        return
end

!==============================================================================
subroutine boundaries_inflowOK(deltaT, iok)
	use boundaries

	iok = 0

	if (ibmin == imin) then
		select case (inflowType)
		case (1)
			call inflow_ok(deltaT, iok)
		case default
			iok = 1
		end select
	end if

	rok = real(iok)
	call globalMax(rok, 1)
	iok = nint(rok)

	return
end
!=============================================================================
subroutine mean_Time(deltaT)
	use boundaries

	MeanTime = MeanTime + deltaT

	return
end
!==============================================================================================
