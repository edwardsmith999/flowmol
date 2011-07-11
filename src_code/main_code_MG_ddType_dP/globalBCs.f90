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
	!real    MeanOmega               ! frequency of oscillation of mean flow    (ubulk + MeanAmp*cos(MeanOmega*MeanTime))
	!real    MeanAmp                 ! amplitude of oscillation of mean flow
	!real    MeanTime                ! time for oscillations of mean flow

	real Uin(0:nly+1,0:ngz+1,3), Uout(0:nly+1,0:ngz+1,3) !u,v,w left-right
	real Vin(0:nlx+1,0:ngz+1,3), Vout(0:nlx+1,0:ngz+1,3) !u,v,w top-bottom
	real Win(0:nlx+1,0:nly+1,3), Wout(0:nlx+1,0:nly+1,3) !u,v,w front-back

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
	call readInt("BC_front" , BC_front )
	call readInt("BC_back"  , BC_back  )

	if ( (BC_bottom.eq.2 .and. BC_top.ne.2) .or. &
	     (BC_top.eq.2 .and. BC_bottom.ne.2)) then
		stop "Both top and bottom boundaries must be periodic"
	end if

	if ( (BC_front.eq.2 .and. BC_back.ne.2) .or. &
             (BC_back.eq.2 .and. BC_front.ne.2)) then
                stop "Both front and back boundaries must be periodic"
        end if

      	if (iblock.eq.1) then
		!Initialize mean time (same as simulation time)
		!call archive_isDefined("time", iflag)
		!if (iflag .eq. 1) then
		!	call readFloat("time", MeanTime)
		!else
		!	MeanTime = 0.0
		!end if

		!print *, "Meantime", MeanTime

      		! Determine base inflow
      		call readInt ("inflow", inflowType)
      		select case (inflowType)
      		case(0) !Laminar
           	call readFloat ("ubulk", ubulk)
		!call readFloat ("angle_attack", angle_attack)
		!call readFloat ("MeanOmega", MeanOmega)
		!call readFloat ("MeanAmp", MeanAmp)
		!call readFloat ("MeanTime", MeanTime)
		call laminarInflowBC()
      		case(1)  !Read inflow from file      
           	call inflow_init()
		call inflow_cache(ny, y, ym, nz, z, zm)
      		end select     
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
	call writeInt("BC_front" , BC_front )
	call writeInt("BC_back"  , BC_back  )


	if (iblock == 1) then
		select case (inflowType)
		case (0)
			call writeFloat("ubulk", ubulk)
			!call writeFloat ("angle_attack", angle_attack)
			!call writeFloat ("MeanOmega", MeanOmega)
			!call writeFloat ("MeanAmp", MeanAmp)
			!call writeFloat ("MeanTime", MeanTime)
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

!     Time advance the inflow BCs and 
!     determine the cartesian velocity BCs.     

      	use boundaries
      
	real*8  :: utemp(ngz+1)
     
	
	!print*, "nlx:", nlx, "nlxb:", nlxb, "nly:", nly	

      	!==================================================================
      	! 	INFLOW B.C.
      	!===================================================================
      	
	!--------------------------------------------------------------
	!	 Left B.C. 
	!-------------------------------------------------------------
	!Copy inflow BCs from Uin to uc,vc,wc
      	if (iblock ==  1 ) then
		do k=0,ngz
                	uc(k,  1  , 0:nly) =     Uin(0:nly, k, 1)
                	vc(k,  0  ,  :   ) = 2.0*Uin( :   , k, 2) - vc(k, 1,  :   )
                	wc(k,  0  , 0:nly) = 2.0*Uin(0:nly, k, 3) - wc(k, 1, 0:nly)
            	end do
                
		wc(ngz+1,  0, 0:nly) = 2.0*Uin(0:nly, ngz+1, 3) - wc(ngz+1, 1, 0:nly)
      		
		!Extend
		uc(:, 0, :) = 2.0*uc(:, 1, :) - uc(:, 2, :)
	end if

      	!==================================================================
      	! 	OUTFLOW B.C.
      	!==================================================================

	!--------------------------------------------------------------
	!	 Right B.C. 
      	!--------------------------------------------------------------
	 !print *, "uc(1," , nlxb  , ",", nlyb-1, "):", uc(1, nlxb  , nlyb-1)
         !print *, "uc(1," , nlxb-1, ",", nlyb-1, "):", uc(1, nlxb-1, nlyb-1)

	!Time-advance outflow BC
      	call outflow_convective(deltaT, Uout)
	
	!Extend
	if (iblock == npx) then
		uc(:, nlxb+1, 1:nlyb-1) = 2.0*uc(:, nlxb, 1:nlyb-1) - uc(:, nlxb-1, 1:nlyb-1)
		!print *, "uc(1," , nlxb+1, ",", nlyb-1, "):", uc(1, nlxb+1, nlyb-1)
		!print *, "uc(1," , nlxb  , ",", nlyb-1, "):", uc(1, nlxb  , nlyb-1)
		!print *, "uc(1," , nlxb-1, ",", nlyb-1, "):", uc(1, nlxb-1, nlyb-1) 
	end if
	
	!do k=0, ngz
	!do j=0, nly
	!	print*, "uc(", k, ",", nlxb+1, ",", j, "):" ,uc(k, nlxb+1, j)
	!end do
	!end do
	
	!=====================================================================
	!	FREE SLIP B.C.
	!=====================================================================
	!--------------------------------------------------------------
	!	 Bottom B.C. 
	!--------------------------------------------------------------
	SELECT CASE (BC_bottom)
	CASE(0)
                !Must set free-slip on fluxes and then compute corresponding Cartesian velocities
                !Setting free-slip BCs on cartesian velocities would be inaccurate as grid not Cartesian
                if (jblock.eq.1) then
                
                        u(:, :, 0) =  u(:, :, 1)
                        w(:, :, 0) =  w(:, :, 1)
                        v(:, :, 1) =  0.0
                        v(:, :, 0) = -v(:, :, 2)

                        !---Convert to Cartesian BC
                        jj=0 
                        j = jbmap_1(jj)
                        do ii=1,nlxb
                                i = ibmap_1(ii)

                                !------ Uc ------
                                a1 = surxix(i,j)
                                a2 = surxiy(i,j)
                                b1 = (svretax(i-1,j)+svretax(i,j)+svretax(i,j+1)+svretax(i-1,j+1))/4.
                                b2 = (svretay(i-1,j)+svretay(i,j)+svretay(i,j+1)+svretay(i-1,j+1))/4.

                                utemp(:) = ( v(:, ii-1, jj) + v(:, ii-1, jj+1)    &
                                           + v(:, ii  , jj) + v(:, ii  , jj+1) ) / 4.
                                uc(:, ii, jj) = a1*u(:, ii, jj) + b1*utemp(:)

                                !------ Vc ------
                                b1 = svretax(i,j)
                                b2 = svretay(i,j)
                                a1 = (surxix(i,j-1)+surxix(i,j)+surxix(i+1,j)+surxix(i+1,j-1))/4.
                                a2 = (surxiy(i,j-1)+surxiy(i,j)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.

                                utemp(:) = ( u(:, ii  , jj-1) + u(:, ii  , jj)    &
                                           + u(:, ii+1, jj-1) + u(:, ii+1, jj) ) / 4.
                                vc(:, ii, jj) = a2*utemp(:) + b2*v(:, ii, jj)

                                !------ Wc ------
                                i  = ibmap_1(ii)
                                a3 = swrz(i,j)
                                wc(:, ii, jj) = a3*w(:, ii, jj)
                        end do

			!Extend
			!extend vc to bottom halo
			vc(:, :, 0) = 2.0*vc(:, :, 1) - vc(:, :, 2)
			
			!extend uc to bottom left halo
			if (iblock == 1  ) then
                        	uc(1:ngz-1, 0     , jj)= 2.0*uc(1:ngz-1, 1   , jj) - uc(1:ngz-1, 2     , jj)
			end if
	
			!extend uc to bottom right halo
			if (iblock == npx) then
		                uc(1:ngz-1, nlxb+1, jj)= 2.0*uc(1:ngz-1, nlxb, jj) - uc(1:ngz-1, nlxb-1, jj)
			end if
               end if

	END SELECT

	!--------------------------------------------------------------
	!	Top B.C. 
	!--------------------------------------------------------------
	SELECT CASE (BC_top)
	CASE(0)
		!Must set free-slip on fluxes and then compute corresponding Cartesian velocities
		!Setting free-slip BCs on cartesian velocities would be inaccurate as grid not Cartesian
		if (jblock.eq.npy) then
		
			u(:, :, nlyb  ) =  u(:, :, nlyb-1)
			w(:, :, nlyb  ) =  w(:, :, nlyb-1)
			v(:, :, nlyb  ) =  0.
			v(:, :, nlyb+1) = -v(:, :, nlyb-1)

			!---Convert to Cartesian BC
			jj= nlyb
			j = jbmap_1(jj)
			do ii=1,nlxb
				i = ibmap_1(ii)

				!--- Uc ---
				a1 = surxix(i,j)
				a2 = surxiy(i,j)
				b1 = (svretax(i-1,j)+svretax(i,j)+svretax(i,j+1)+svretax(i-1,j+1))/4.
				b2 = (svretay(i-1,j)+svretay(i,j)+svretay(i,j+1)+svretay(i-1,j+1))/4.

				utemp(:) = ( v(:, ii-1, jj) + v(:, ii-1, jj+1)    &
				           + v(:, ii  , jj) + v(:, ii  , jj+1) ) / 4.
				uc(:, ii, jj) = a1*u(:, ii, jj) + b1*utemp(:)

				!--- Vc ---
				b1 = svretax(i,j)
				b2 = svretay(i,j)
				a1 = (surxix(i,j-1)+surxix(i,j)+surxix(i+1,j)+surxix(i+1,j-1))/4.
				a2 = (surxiy(i,j-1)+surxiy(i,j)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.

				utemp(:) = ( u(:, ii  , jj-1) + u(:, ii  , jj)    &
					   + u(:, ii+1, jj-1) + u(:, ii+1, jj) ) / 4.
				vc(:, ii, jj) = a2*utemp(:) + b2*v(:, ii, jj)

				!--- Wc ---
				i  = ibmap_1(ii)
				a3 = swrz(i,j)
				wc(:, ii, jj) = a3*w(:, ii, jj)
			end do

			!Extend
			!extend vc to top halo
               	 	vc(:, :, nlyb+1) = 2.0*vc(:, :, nlyb) - vc(:, :, nlyb-1)

			!extend uc to top left halo
			if (iblock == 1  ) then
                        	uc(1:ngz-1, 0     , jj)= 2.0*uc(1:ngz-1, 1   , jj) - uc(1:ngz-1, 2     , jj)
			end if

			!extend uc to top right halo
			if (iblock == npx) then
	                        uc(1:ngz-1, nlxb+1, jj)= 2.0*uc(1:ngz-1, nlxb, jj) - uc(1:ngz-1, nlxb-1, jj)
			end if

		end if
      
      	END SELECT

	!--------------------------------------------------------------
	!	Front B.C. 
	!--------------------------------------------------------------
	SELECT CASE (BC_front)
	CASE(0)
		!Can apply free-slip on Cartesian velocities as grid not curvilinear in z
		! ------ Uc ------
	  	uc(0, :, :) = uc(1, :, :)
         	
		! ------ Vc ------
		vc(0, :, :) = vc(1, :, :)
         	
		! ------ Wc ------
		wc(1, :, :) = 0.0

		!Extend
        	wc(0, :, :) = -wc(1, :, :)
         
		
	END SELECT

	!--------------------------------------------------------------
        !       Back B.C. 
        !--------------------------------------------------------------
	SELECT CASE (BC_back)
	CASE(0)
		!Can apply free-slip on Cartesian velocities as grid not curvilinear in z
		
		! ------ Uc ------
	  	uc(ngz  , :, :) = uc(ngz-1, :, :)
         	
		! ------ Vc ------
		vc(ngz  , :, :) = vc(ngz-1, :, :)
         
		! ------ Wc ------
		wc(ngz  , :, :) = 0.0 
         	
		!Extend
		wc(ngz+1, :, :) = -wc(ngz , :, :)
		
		
	END SELECT

      	!===================================================================
      	! 	NO SLIP B.C. 
      	!===================================================================
	
	!--------------------------------------------------------------
	!	Bottom B.C. 
	!--------------------------------------------------------------
	SELECT CASE (BC_bottom)
      	CASE(1)
		if (jblock.eq.1) then

			! ----- Uc ------
			uc(:, :, 0) = -uc(:, :, 1)
         	
			! ----- Vc ------
			vc(:, :, 1) = 0.0
		
			!Extend
			vc(:, :, 0) = 0.0
         
			! ----- Wc ------
			wc(:, :, 0) = -wc(:, :, 1)
         	
		end if
	END SELECT

	!--------------------------------------------------------------
	!	 Top B.C.
	!-------------------------------------------------------------- 
	SELECT CASE (BC_top)
        CASE(1)
		if (jblock.eq.npy) then
         	
			!------ Uc ------
         		uc(:, :, nlyb  ) = -uc(:, :, nlyb-1)
   
			!------ Vc ------
			vc(:, :, nlyb  ) = 0.0
         		
			!Extend
			vc(:, :, nlyb+1) = 0.0 
         	
			!------ Wc ------
			wc(:, :, nlyb  ) = -wc(:, :, nlyb-1)
         		
		end if
	END SELECT

	!--------------------------------------------------------------
	!	Front B.C. 
	!--------------------------------------------------------------
	SELECT CASE (BC_front)
        CASE(1)
		! ------ Uc ------
	  	uc(0, :, :) = -uc(1, :, :)
         	
		! ------ Vc ------
		vc(0, :, :) = -vc(1, :, :)
         	
		! ------ Wc ------
		wc(1, :, :) = 0.0
		
		!Extend
        	wc(0, :, :) = 0.0

         END SELECT

	!--------------------------------------------------------------
	!	Back B.C. 
	!--------------------------------------------------------------
	SELECT CASE (BC_back)
        CASE(1)
		! ------ Uc ------
	  	uc(ngz  , :, :) = -uc(ngz-1, :, :)
         	
		! ------ Vc ------
		vc(ngz  , :, :) = -vc(ngz-1, :, :)
         
		! ------ Wc ------
		wc(ngz  , :, :) = 0.0 

		!Extend
         	wc(ngz+1, :, :) = 0.0
	
	END SELECT

	!=======================================================================
	!	PERIODIC B.C.
	!=======================================================================
      	
	!-----------------------------------------------------------------------
      	!    Bottom and Top B.C.  !!!USE ONLY WITH PERIODIC POISSON ROUTINE!!!
      	!-----------------------------------------------------------------------
	SELECT CASE (BC_bottom)
      	CASE(2)
		print*, 'USE ONLY WITH PERIODIC POISSON ROUTINE!!!'

		if (jblock.eq.1 .or. jblock.eq.npy) then
                        if (npy.eq.1) then
                                ucbcs(:,:,1) = uc(:,:,1     )
                                wcbcs(:,:,1) = wc(:,:,1     )
                                vcbcs(:,:,1) = vc(:,:,1     )
                                vcbcs(:,:,2) = vc(:,:,2     )

                                ucbcn(:,:,1) = uc(:,:,nlyb-1)
                                wcbcn(:,:,1) = wc(:,:,nlyb-1)
                                vcbcn(:,:,1) = vc(:,:,nlyb-1)
                                vcbcn(:,:,2) = vc(:,:,nlyb  )
                        else
                                !-- Here, I use 'ucbcs/n' for saving the results form
                                !   boundary conditions on South/North faces respectively.
                                !   ==> (jblock=1  )  receives ucbcn
                                !   ==> (jblock=npy)  receivss ucbcs
                                call updateBorder_yperiodic(1)
                        end if
                end if
	
		if (jblock.eq.1) then

			! ----- Uc -----
             		uc(:, :, 0) = ucbcn(:, :, 1)
         	
			! ----- Vc ------
 			vc(:, :, 1) = 0.5*(vc(:, :, 1)+vcbcn(:, :, 2))
   	
			!Extend
			vc(:, :, 0) = vcbcn(:, :, 1)
         	
			! ----- Wc ------
             		wc(:, :, 0) = wcbcn(:, :, 1)
         	
		end if

		if (jblock.eq.npy) then

                        !------ Uc ------
                        uc(:, :, nlyb  ) = ucbcs(:, :, 1)

                        !------ Vc ------
                        vc(:, :, nlyb  ) = 0.5*(vc(:, :, nlyb) + vcbcs(:, :, 1))

                        !Extend
                        vc(:, :, nlyb+1) = vcbcs(:, :, 2)

                        !------ Wc ------
                        wc(:, :, nlyb  ) = wcbcs(:, :, 1)

                end if

	END SELECT

	!-----------------------------------------------------------------------
        !       Front and Back B.C.
        !-----------------------------------------------------------------------
	SELECT CASE (BC_front)
      	CASE(2)
		! ------ Uc ------
		uc(0    , :, :) = uc(ngz-1, :, :)
		uc(ngz  , :, :) = uc(1    , :, :)

		! ------ Vc ------
		vc(0    , :, :) = vc(ngz-1, :, :)
		vc(ngz  , :, :) = vc(1    , :, :)
		! ------ Wc ------
		wc(0    , :, :) = wc(ngz-1, :, :)
		wc(ngz  , :, :) = wc(1    , :, :)
		!DANGER ZONE: wc(  1  , :,:) = 0.5*(wc(1, :,:)+wc(ngz, :,:)) (wc only goes from 1:ngz-1 so wc(ngz,:,:) does not exist)
		!Extend
                wc(ngz+1, :, :) = wc(2    , :, :)

	END SELECT

	!do k=0, ngz
        !do j=0, nly
        !       print*, "uc(", k, ",", nlxb+1, ",", j, "):" ,uc(k, nlxb+1, j)
	!	print*, "uc(", k, ",", nlxb  , ",", j, "):" ,uc(k, nlxb  , j)
	!	print*, "uc(", k, ",", nlxb-1, ",", j, "):" ,uc(k, nlxb-1, j)
        !end do
        !end do

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
     

      	real*8  :: vtempb(ngz+1), vtempe(ngz+1)
      	real*8  :: utempb(ngz+1), utempe(ngz+1)
      	real*8, dimension(2) :: Atmp


      	!======================================================================
      	! 	INFLOW
      	!======================================================================

	!--------------------------------------------------------------
	!	Left B.C.
	!--------------------------------------------------------------
      	if (iblock.eq.1) then

		!----------  U^{xi} (i.e. u) -----------
      		do j=0,nlyb
           		jj = jbmap_1(j)
           		vtempb = ( vc(:, 0, j+1) + vc(:, 1, j+1)  &
                                 + vc(:, 0, j  ) + vc(:, 1, j  ) ) /4.
           		u(:, 1, j) = suxix(1, jj)*uc(:, 1, j) + suxiy(1, jj)*vtempb
      		end do
		
		!Extend
		u(:, 0 , 0:nly)=2.0*u(:, 1,  0:nly)-u(:, 2,  0:nly)
      		
		!----------  U^{eta} (i.e. v) -----------
      		do j=0,nlyb+1
            		jj = jbmap_1(j)
            		utempb = ( uc(:, 1, j-1) + uc(:, 1, j)  &
                                 + uc(:, 0, j-1) + uc(:, 0, j) ) /4.
            		v(:, 0, j) = svetax(0 ,jj)*utempb + svetay(0, jj)*vc(:, 0, j)
      		end do
      
      		!----------  U^{z} (i.e. w) -----------
      		do j=0,nlyb+1
            		jj = jbmap_1(j)
            		a3 = swz(0,jj)
            		w(:, 0, j) = a3*wc(:, 0, j)
      		end do
      		if (jblock.eq.1)  &
      		w(:, 0, 0   ) = 2.0*w(:, 1, 0   ) - w(:, 2, 0   )
      		if (jblock.eq.npy)  &
      		w(:, 0, nlyb) = 2.0*w(:, 1, nlyb) - w(:, 2, nlyb)

		
	end if
      
      	!========================================================================
      	!		OUTFLOW
      	!========================================================================
      	
	!--------------------------------------------------------------
        !       Right B.C.
        !--------------------------------------------------------------
	if (iblock.eq.npx) then

		!----------  U^{xi} (i.e. u) -----------
		do j=0,nlyb
			jj = jbmap_1(j)
			vtempe = ( vc(:, nlxb-1, j+1) + vc(:, nlxb, j+1)  &
				 + vc(:, nlxb-1, j  ) + vc(:, nlxb, j  ) ) /4.
			u(:, nlxb, j) = suxix(ngx, jj)*uc(:, nlxb, j) + suxiy(ngx, jj)*vtempe
		end do
		
		!Extend
		u(:, nlxb+1, 0:nly) = 2.0*u(:, nlxb, 0:nly)-u(:, nlxb-1, 0:nly)

		!print*, "u inside FluxBC"
		!do k=0,ngz
		!do j=0,nly
		!	print *, "u(", k, ",", nlxb+1, ",", j, ") :" ,u(k, nlxb+1, j) 	
		!end do
		!end do
	
		!----------  U^{eta} (i.e. v) -----------
		do j=0,nlyb+1
			jj = jbmap_1(j)
			utempe = ( uc(:, nlxb+1, j-1) + uc(:, nlxb+1, j)  &
			         + uc(:, nlxb  , j-1) + uc(:, nlxb  , j) ) /4.
			v(:, nlxb, j) = svetax(ngx, jj)*utempe + svetay(ngx, jj)*vc(:, nlxb, j)
		end do

		!----------  U^{z} (i.e. w) -----------
		do j=0,nlyb
			jj = jbmap_1(j)
			a3 = swz(ngx,jj)
			w(:, nlxb, j) = a3*wc(:, nlxb, j)
		end do
                if (jblock.eq.1)  &
                w(:, nlxb, 0   ) = 2.0*w(:, nlxb-1, 0   ) - w(:, nlxb-2, 0   )
                if (jblock.eq.npy)  &
                w(:, nlxb, nlyb) = 2.0*w(:, nlxb-1, nlyb) - w(:, nlxb-2, nlyb)

		
	end if

	!=======================================================================
      	!	FREE SLIP B.C. 	
      	!=======================================================================
	
	!---------------------------------------------------------------------
	!	Bottom B.C.
	!---------------------------------------------------------------------
	SELECT CASE (BC_bottom) 
	CASE(0)
		if (jblock.eq.1) then
         	
        	!----------  U^{xi} (i.e. u) -----------
		do i=0,nlxb+1
           		u(:, i, 0) = u(:, i, 1)
         	end do
         	
		!----------  U^{eta} (i.e. v) -----------
	 	do i=0,nlxb
              		v(:, i, 1) = 0.0
              		
			!Extend
			v(:, i, 0) = 2.0*v(:, i, 1) -v(:, i, 2)
         	end do
         	
		!----------  U^{z} (i.e. w) -----------	
		do i=0,nlxb
           		w(:, i, 0) = w(:, i, 1)
         	end do
         	
		end if
	END SELECT

	!---------------------------------------------------------------------
        !       Top B.C.
        !---------------------------------------------------------------------
        SELECT CASE (BC_top)
        CASE(0)
		if (jblock.eq.npy) then
         	
			!----------  U^{xi} (i.e. u) -----------
         		do i=0,nlxb+1
           			u(:, i, nlyb  ) = u(:, i, nlyb-1)
         		end do
         	
			!----------  U^{eta} (i.e. v) -----------   
         		do i=0,nlxb
              			v(:, i, nlyb  ) = 0.
              		
				!Extend
				v(:, i, nlyb+1) = 2.0*v(:, i, nlyb) - v(:, i, nlyb-1)
         		end do
         	
			!----------  U^{z} (i.e. w) -----------
      			do i=0,nlxb
             			w(:, i, nlyb  ) = w(:, i, nlyb-1)
      			end do
      		
		end if 
	END SELECT

	!------------------------------------------------------------------
	!	Front B.C.
	!------------------------------------------------------------------
	SELECT CASE (BC_front)
	CASE(0)	
		!-----------U^{xi} (i.e. u)-------------------
		do i=0,nlxb+1     
         		u(0, i, :) = u(1, i, :)
      		end do

		!------------ U^{eta} (i.e. v) ----------------
		do i=0,nlxb     
         		v(0, i, :) = v(1, i, :)
      		end do
	
		!------------ U^{z} (i.e. w) -----------------
		do i=0,nlxb
              		w(1, i, :) = 0.0
              	end do
		!Extend
		do i=0,nlxb
                        w(0, i, :) =  w(1, i, :)
                end do

	END SELECT
      	
	!------------------------------------------------------------------
        !       Back B.C.
        !------------------------------------------------------------------
	SELECT CASE (BC_back)
        CASE(0)
		!-----------U^{xi} (i.e. u)-------------------
		do i=0,nlxb+1     
         		u(ngz  , i, :) = u(ngz-1, i, :)
      		end do
      
		!------------ U^{eta} (i.e. v) ----------------
		do i=0,nlxb     
         		v(ngz  , i, :) = v(ngz-1, i, :)
      		end do
      
		!------------ U^{z} (i.e. w) -----------------
		do i=0,nlxb
              		w(ngz  , i, :) = 0.0
         	end do
		!Extend
		do i=0,nlxb
                        w(ngz+1, i, :) =  w(ngz  , i, :)
                end do

	END SELECT

	!=======================================================================
        !       NO SLIP B.C.  
        !=======================================================================

	!-----------------------------------------------------------------------
      	!	Bottom B.C.     
      	!-----------------------------------------------------------------------
	SELECT CASE(BC_bottom)
	CASE(1)
		if (jblock.eq.1) then
         	
		!----------  U^{xi} (i.e. u) -----------
		do i=0,nlxb+1
           		u(:, i, 0) = -u(:, i, 1)
         	end do
         	
		!----------  U^{eta} (i.e. v) -----------
         	do i=0,nlxb
              		v(:, i, 1) = 0.0
              		!Extend
			v(:, i, 0) = 0.0 
         	end do
         	
		!----------  U^{z} (i.e. w) -----------	
		do i=0,nlxb
           		w(:, i, 0) = -w(:, i, 1)
         	end do
         	
		end if
	END SELECT

	!-----------------------------------------------------------------------
        !       Top B.C.     
        !-----------------------------------------------------------------------
        SELECT CASE(BC_top)
        CASE(1)
		if (jblock.eq.npy) then
         	
			!----------  U^{xi} (i.e. u) -----------
         		do i=0,nlxb+1
           			u(:, i, nlyb  ) = -u(:, i, nlyb-1)
         		end do
         	
			!----------  U^{eta} (i.e. v) -----------   
         		do i=0,nlxb
              			v(:, i, nlyb  ) = 0.0
				!Extend
				v(:, i, nlyb+1) = 0.0
         		end do
         	
			!----------  U^{z} (i.e. w) -----------
      			do i=0,nlxb
             			w(:, i, nlyb  ) = -w(:, i, nlyb-1)
      			end do
      		
		end if 
	END SELECT
	
	!-----------------------------------------------------------------------
        !       Front B.C.     
        !-----------------------------------------------------------------------
        SELECT CASE(BC_front)
        CASE(1)
		!-----------U^{xi} (i.e. u)-------------------
		do i=0,nlxb+1
         		u(0, i, :) = -u(1, i, :)
      		end do
	
		!------------ U^{eta} (i.e. v) ---------------
		do i=0,nlxb
         		v(0, i, :) = -v(1, i, :)
      		end do
      
		!------------ U^{z} (i.e. w) -----------------
		do i=0,nlxb
              		w(1, i, :) = 0.0
			!Extend
			w(0, i, :) = 0.0
         	end do
	END SELECT

	!-----------------------------------------------------------------------
        !       Back B.C.     
        !-----------------------------------------------------------------------
        SELECT CASE(BC_back)
        CASE(1)
		
		!-----------U^{xi} (i.e. u)-------------------
		do i=0,nlxb+1
         		u(ngz  , i, :) = -u(ngz-1, i, :)
      		end do

		!------------ U^{eta} (i.e. v) ---------------
		do i=0,nlxb
         		v(ngz  , i, :) = -v(ngz-1, i, :)
      		end do
	
		!------------ U^{z} (i.e. w) -----------------
		do i=0,nlxb
              		w(ngz  , i, :) = 0.0
			!Extend
              		w(ngz+1, i, :) = 0.0
         	end do
	END SELECT

	!==========================================================================
	!	PERIODIC B.C.
	!==========================================================================

	!-----------------------------------------------------------------------
        !       Bottom and Top B.C.     
        !-----------------------------------------------------------------------
        SELECT CASE(BC_bottom)
        CASE(2)
      		!!!USE ONLY WITH PERIODIC POISSON ROUTINE!!!
		if (jblock.eq.1 .or. jblock.eq.npy) then
                if (npy.eq.1) then
                        vcbcs(:,:,1) = v (:,:,1)
                        vcbcs(:,:,2) = v (:,:,2)
			ucbcs(:,:,1) = u (:,:,1)
                  	!SY this is a bug. FIXED 
                  	!  wcbcs(:,:,2) = w (:,:,1)==> wcbcs(:,:,1) = w (:,:,1)
                        wcbcs(:,:,1) = w (:,:,1)

                        vcbcn(:,:,1) = v (:,:,nlyb-1)
                        vcbcn(:,:,2) = v (:,:,nlyb  )
			ucbcn(:,:,1) = u (:,:,nlyb-1)
                  	!SY this is a bug. FIXED 
                  	!  wcbcn(:,:,2) = w (:,:,nlyb-1)==> wcbcn(:,:,1) = w (:,:,nlyb-1)
                        wcbcn(:,:,1) = w (:,:,nlyb-1)
                else
                        !-- Here, I use 'ucbcs/n' for saving the results form
                        !   boundary conditions on South/North faces respectively.
                        !   ==> (jblock=1  )  receives ucbcn
                        !   ==> (jblock=npy)  receivss ucbcs
                        call updateBorder_yperiodic(2)
                end if
        	end if

        	if (jblock.eq.1) then
         	
		!----------  U^{xi} (i.e. u) -----------
		do i=0,nlxb+1
           		u(:, i, 0) = ucbcn(:, i, 1)
         	end do
         	
		!----------  U^{eta} (i.e. v) -----------
         	do i=0,nlxb
              		v (:, i, 0) = vcbcn(:, i, 1)
                  	!SY this is a bug. FIXED 
		  	! v (:, i, 1) = 0.5*(vc(:, i, 1) + vcbcn(:, i, 2)) ==>
	          	! v (:, i, 1) = 0.5*(v(:, i, 1) + vcbcn(:, i, 2))
			v (:, i, 1) = 0.5*(v(:, i, 1) + vcbcn(:, i, 2))
         	end do
         	
		!----------  U^{z} (i.e. w) -----------	
		do i=0,nlxb
           		w(:, i, 0) = wcbcn(:, i, 1)
         	end do
         	
		end if

		if (jblock.eq.npy) then

                        !----------  U^{xi} (i.e. u) -----------
                        do i=0,nlxb+1
                                u(:, i, nlyb  ) = ucbcs(:, i, 1)
                        end do

                        !----------  U^{eta} (i.e. v) -----------

                        do i=0,nlxb
                                v(:, i, nlyb  ) = 0.5*(v(:, i, nlyb) + vcbcs(:, i, 1))
                        end do
                        do i=0,nlxb
                                v(:, i, nlyb+1) = vcbcs(:, i, 2)
                        end do

                        !----------  U^{z} (i.e. w) -----------

                        do i=0,nlxb
                                w(:, i, nlyb  ) = wcbcs(:, i, 1)
                        end do

                end if

	END SELECT


	!-----------------------------------------------------------------------
        !       Front and Back B.C.     
        !-----------------------------------------------------------------------
        SELECT CASE(BC_front)
        CASE(2)
		!-----------U^{xi} (i.e. u)-------------------
      		do i=0,nlxb+1
			u(0    , i, :) = u(ngz-1, i, :)
			u(ngz  , i, :) = u(1    , i, :)
      		end do
		
		!------------ U^{eta} (i.e. v) ---------------
     		do i=0,nlxb   
        		v(0    , i, :) = v(ngz-1, i, :)
			v(ngz  , i, :) = v(1    , i, :)
      		end do
		
		!------------ U^{z} (i.e. w) -----------------	
		do i=0,nlxb
        		!DANGER ZONE: w(1    , i, :) = 0.5*(w(1, i, :)+ w(ngz, i, :))
                        w(ngz  , i, :) = w(1    , i, :)

      		end do
		!Extend
		 do i=0,nlxb
                        w(0    , i, :) = w(ngz-1, i, :)
                        w(ngz+1, i, :) = w(2    , i, :)

                end do

	END SELECT


        !=======================================================================
        ! sum of mass flux along boundaries is zero
        !=======================================================================
        bmassin =0.0
        bmassout=0.0
        if (jblock.eq.1  ) then
        do i=i1_T,i2_T
           bmassin  = bmassin  + sum(v(1:ngzm, i, 1   ))
        end do
        end if

        if (jblock.eq.npy) then
        do i=i1_T,i2_T
           bmassin  = bmassin  - sum(v(1:ngzm, i, nlyb))
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
           bmassin  = bmassin  + w(1,i,j)-w(ngz,i,j)
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

        if(abs(bmassout).ge.1.e-10) then 
           if (iblock.eq.npx) then
               u(1:ngzm,nlxb,1:nly-1)=u(1:ngzm,nlxb,1:nly-1)*bmassin/bmassout
           end if
        end if

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

        if (iblock.eq.npx) then
                !Extend
                u(:, nlxb+1, 0:nly) = 2.0*u(:, nlxb, 0:nly)-u(:, nlxb-1, 0:nly)
        end if

        !TAZ DANGER ZONE:  Xiahua has u-periodicity here:
	if (BC_front.eq.2) then
        	u( 0 , :, :) = u(ngz-1, :, :)
        	u(ngz, :, :) = u(1    , :, :)
        end if

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

        !T  Need a halo cell in {ust,vst,wst} before poisson solver
        !T-- x-direction
        call updateBorder_lim(ust, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 1)
        call updateBorder_lim(vst, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 1)
        call updateBorder_lim(wst, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 1)

        !T-- y-direction
        call updateBorder_lim(ust, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 1)
        call updateBorder_lim(vst, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 1)
        call updateBorder_lim(wst, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 1)

        !T Also need to copy {u,v,w} on the border to {ust,vst,wst} for Divergence
        !T-- Left Border
        if (iblock_1.eq.1) then
                ust(:,i1_u-1,:) = u(:,i1_u-1,:)
                vst(:,i1_v-1,:) = v(:,i1_v-1,:)
                wst(:,i1_v-1,:) = w(:,i1_v-1,:)
        end if

        !T-- Right Border
        if (iblock_1.eq.npx) then
                ust(:,i2_u+1,:) = u(:,i2_u+1,:)
                vst(:,i2_v+1,:) = v(:,i2_v+1,:)
                wst(:,i2_v+1,:) = w(:,i2_v+1,:)
        end if

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

        !SY-- front Border
	select case (BC_front)
	case(0:1) ! Slip or no-slip bc
                wst(1  ,:,:) = w(1  ,:,:)
	case(2)   ! Periodic bc
	end select

        !SY-- Back Border
	select case (BC_back)
	case(0:1) ! Slip or no-slip bc
                wst(ngz,:,:) = w(ngz,:,:)
	case(2)   ! Periodic bc
        	wst(ngz,:,:) = wst(1,:,:)
	end select

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
		call laminarInflowBC()
		select case (inflowFluct)
		case (0)
		case (1:)
			call fluct_Time(deltaT)
		end select

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

		! Boundary conditions for Uin
		
		if (jblock == 1) then
		select case (BC_bottom)
		case(0) ! slip    at lower wall
			Uin(0,:,1) =  Uin(1,:,1)
			Uin(0,:,3) =  Uin(1,:,3)
			Uin(1,:,2) =  Uin(2,:,2)
			Uin(0,:,2) =  Uin(1,:,2)
		case(1) ! No-slip at lower wall
			Uin(0,:,1) = -Uin(1,:,1)
			Uin(0,:,3) = -Uin(1,:,3)
			Uin(1,:,2) = 0.
			Uin(0,:,2) = 0.
		end select
		end if

		!------------------------------------------------
		!   V at j=jmaxp1 not used
		!  if (jbmax == jmax) Uin(jmap_1(jmaxo),:,2) = 0.
		!------------------------------------------------

		if (jblock == npy) then
		select case (BC_top)
		case(0) ! slip    at upper wall
			Uin(nlyb  ,:,1) =  Uin(nlyb-1,:,1)
			Uin(nlyb  ,:,3) =  Uin(nlyb-1,:,3)
			Uin(nlyb  ,:,2) =  Uin(nlyb-1,:,2)
			Uin(nlyb+1,:,2) =  Uin(nlyb  ,:,2)
		case(1) ! No-slip at upper wall
			Uin(nlyb  ,:,1) = -Uin(nlyb-1,:,1)
			Uin(nlyb  ,:,3) = -Uin(nlyb-1,:,3)
			Uin(nlyb  ,:,2) = 0.
			Uin(nlyb+1,:,2) = 0.
		end select
		end if

		if (npz == 1) then

		select case (BC_front)
		case(0) ! Slip    bc
			Uin(:,kmino,1) =  Uin(:,kmin  ,1)
			Uin(:,kmino,2) =  Uin(:,kmin  ,2)
			Uin(:,kmin ,3) =  Uin(:,kmin+1,3)
			Uin(:,kmino,3) =  Uin(:,kmin  ,3)
		case(1) ! No-slip bc
			Uin(:,kmino,1) = -Uin(:,kmin  ,1)
			Uin(:,kmino,2) = -Uin(:,kmin  ,2)
			Uin(:,kmin ,3) = 0.
			Uin(:,kmino,3) = 0.
		case(2) ! Periodicity in z
			Uin(:,kmino,:) =  Uin(:,kmax-1,:)
		end select
		select case (BC_back)
		case(0) ! Slip    bc
			Uin(:,kmax ,1) =  Uin(:,kmax-1,1)
			Uin(:,kmaxo,1) =  Uin(:,kmax  ,1)
			Uin(:,kmax ,2) =  Uin(:,kmax-1,2)
			Uin(:,kmaxo,2) =  Uin(:,kmax  ,2)
			Uin(:,kmax ,3) =  Uin(:,kmax-1,3)
			Uin(:,kmaxo,3) =  Uin(:,kmax  ,3)
		case(1) ! No-slip bc
			Uin(:,kmax ,1) = -Uin(:,kmax-1,1)
			Uin(:,kmaxo,1) =  Uin(:,kmax  ,1)
			Uin(:,kmax ,2) = -Uin(:,kmax-1,2)
			Uin(:,kmaxo,2) =  Uin(:,kmax  ,2)
			Uin(:,kmax ,3) = 0.
			Uin(:,kmaxo,3) = 0.
		case(2) ! Periodicity in z
			Uin(:,kmax ,:) =  Uin(:,kmin  ,:)
			Uin(:,kmaxo,:) =  Uin(:,kmin+1,:)
		end select

		endif
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
subroutine laminarInflowBC()
	use boundaries

        real yd, twopi, amp, wave_y, wave_z, omega_t, T_wave

        yd    = ypg(1,ngy)-ypg(1,1)
        twopi = 4.*acos(0.)
        amp = 0.01 * ubulk  ! amplitude of the wave
        T_wave = 1.0        ! time period of the wave
        omega_t = stime/T_wave

        ja = max(jbmino,jmino)
        jb = min(jbmaxo,jmaxo)
        if(jblock==1  ) ja = 1
        if(jblock==npy) jb = ngy-1
        do j=ja,jb
                jj = jmap_1(j)
                do k=1,ngzm
                wave_y=ypu(1,j)/yd
                wave_z=zpu(k)/zL
                !Uin(jj,k,1) = ubulk + amp * sin(twopi*wave_y)
                !Uin(jj,k,1) = ubulk + amp * sin(twopi*wave_z)
                !Uin(jj,k,1) = ubulk + amp * sin(twopi*omega_t)
                !Uin(jj,k,1) = ubulk + amp * sin(twopi*(wave_z-omega_t))
                enddo
        enddo

	!ja = max(jbmino,jmino)
	!jb = min(jbmaxo,jmaxo)
	!if (ja .lt.  1 )  ja=1			
	!if (jb .gt. ngy)  jb=ngy		
	!do j=ja,jb
	!	jj = jmap_1(j)
	!	do k=1,ngz
	!		Uin(jj,k,1) = (ubulk + MeanAmp*cos(MeanOmega*MeanTime))*cos(angle_attack)
	!		Uin(jj,k,2) = (ubulk + MeanAmp*cos(MeanOmega*MeanTime))*sin(angle_attack)
	!		Uin(jj,k,3) = 0
	!	enddo
	!enddo
	
	!Uin(:,:,1) = (ubulk + MeanAmp*cos(MeanOmega*MeanTime))*cos(angle_attack)
	!Uin(:,:,2) = (ubulk + MeanAmp*cos(MeanOmega*MeanTime))*sin(angle_attack)
	!Uin(:,:,3) = 0


	Uin(:,:,1) = ubulk  ! uniform flow
	Uin(:,:,2) = 0.
	Uin(:,:,3) = 0.

	return
end
!==========================================================================================
subroutine mean_Time(deltaT)
	use boundaries

	MeanTime = MeanTime + deltaT

	return
end
!==============================================================================================
