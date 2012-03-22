!=======================================================================
! Boundary conditions for inflow-outflow boundary layer flow
!
! boundaries_init()
! boundaries_read()
! boundaries_write()
!
! CartesianBC(deltaT)
! FluxBC()
! xzplaneBC(S, nj)
! poissonCoeffBC(ay, by, cy)
! poissonMatrixBC(A, B, C)
! poissonRhsBC()
! YsystemBC(A,B,C,R, iuvw)
!
! timeDependent_Inlet_BC(deltaT)
! boundaries_inflowOK(deltaT, iok)
! laminarInflowBC()*
!

module boundaries
        use data_export
        use mesh_export
        use coupler_cfd_global_data, only : use_coupling

        integer inflowType                ! type of inflow condition
        integer inflowFluct                ! type of inflow fluctuations
        real    ubulk                        ! bulk velocity for laminar inflow
        real, save :: Re                ! Re at inlet location

        real Uin(0:nly+1,0:ngz+1,3), Uout(0:nly+1,0:ngz+1,3)
        real Vin(0:nlx+1,0:ngz+1,3), Vout(0:nlx+1,0:ngz+1,3)

        integer jref
        real    pref

end module

!=======================================================================
subroutine boundaries_init()
        use boundaries

        ! Root process is assumed to be on inflow boundary
        call data_isRoot(iflag)
        if (iflag == 1) then
                if (iblock .ne. 1) &
                        stop "Root process needs to be on inflow boundary"
        end if

        ! Setup the inflow and outflow BC
        ! Initialize to zero
        Uin = 0.
        Uout = 0.
        Vin = 0.
        Vout = 0.

        ! Define Falkner-Skan parameters
        call readFloat("Re",Re)
        !-------------------------------------
        if (iblock.eq.1) then
                call readInt("inflow", inflowType)
                select case (inflowType)
                case (0) ! Laminar
                        call readFloat("ubulk", ubulk)
                        call laminarInflowBC()
                case (1) ! Read inflow from data file
                        call inflow_init()
                        call inflow_cache(ny, y, ym, nz, z, zm)
                end select

                ! Type of fluctuations to superimpose on Uin
                call readInt("inflowFluct", inflowFluct)
                call fluct_init(Uin)
        end if

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

        if (iblock == 1) then
                call writeInt("inflow", inflowType)
                select case (inflowType)
                case (0)
                        call writeFloat("Re",Re)
                        call writeFloat("ubulk", ubulk)
                case (1)
                        call inflow_write()
                end select

                ! Write information about inlet fluctuations
                call fluct_write()
        end if

        return
end

!=======================================================================
subroutine CartesianBC(deltaT)
        use boundaries
        real*8  :: utemp(ngz-1)

        write(0,*) 'CFD: CartesianBC: iblock, jblock, use_coupling ', iblock, jblock, use_coupling 

                !-------------------------------------------------
                !    Copy Inflow BC from (Uin) into (uc,vc,wc)
                !-------------------------------------------------
                if (iblock ==  1 ) then
                        do k=1,ngz-1
                                uc(k,  1  , 0:nly) = Uin (0:nly,k,1)
                                vc(k,  0  ,  :   ) = 2.0*Uin( :   ,k,2) - vc(k,1,  :  )
                                wc(k,  0  , 0:nly) = 2.0*Uin(0:nly,k,3) - wc(k,1,0:nly)
                        end do
                        wc(ngz,  0  , 0:nly) = 2.0*Uin(0:nly,ngz,3) - wc(ngz,1,0:nly)
                end if

                !--------------------------------------------
                !         Time-advance outflow BC
                !--------------------------------------------
                call outflow_convective(deltaT, Uout)

                if (jblock.eq.1) then

                        write(0,*) "CFD: use_coupling" , use_coupling
                        if ( use_coupling ) then

                                call md_vel(uc,vc,wc)
                                
                        else
                !--------------------------------------------
                !         No-slip at lower wall
                !---

                                uc(1:ngz-1, :, 0 )=-uc(1:ngz-1, :, 1)
                                vc(1:ngz-1, :, 1 )=0.
                                wc(1:ngz  , :, 0 )=-wc(1:ngz  , :, 1)
                         endif
                end if

                !--------------------------------------------
                !         Upper boundary condition
                !                (Neuman)
                !--------------------------------------------
                if (jblock == npy) then
                        !----------------------------------------
                        ! This section (Sec F) assumes knowing
                        ! Flux BC and calculating Cartesian BC.
                        ! (Use this section with No-Stress top BC)
                        !----------------------------------------
                        !--- NO Stress BC on Fluxes
                        u(1:ngz-1,:, nlyb ) =  u(1:ngz-1,:,nlyb-1)
                        w( 1:ngz ,:, nlyb ) =  w( 1:ngz ,:,nlyb-1)
                        v(1:ngz-1,:, nlyb ) =  0.
                        v(1:ngz-1,:,nlyb+1) = -v(1:ngz-1,:,nlyb-1)

                        !---Convert to Cartesian BC
                        jj= nlyb
                        j = jbmap_1(jj)
                        do ii=1,nlxb
                                i = ibmap_1(ii)

                                !--- uc ---
                                a1=surxix(i,j)
                                a2=surxiy(i,j)
                                b1=(svretax(i-1,j)+svretax(i,j)+svretax(i,j+1)+svretax(i-1,j+1))/4.
                                b2=(svretay(i-1,j)+svretay(i,j)+svretay(i,j+1)+svretay(i-1,j+1))/4.

                                utemp(1:ngz-1)=(v(1:ngz-1, ii-1, jj)+v(1:ngz-1, ii-1, jj+1)    &
                                                +v(1:ngz-1, ii  , jj)+v(1:ngz-1, ii  , jj+1))/4.
                                uc(1:ngz-1, ii,jj)=a1*u(1:ngz-1,ii,jj)+b1*utemp(1:ngz-1)

                                !--- vc ---
                                b1=svretax(i,j)
                                b2=svretay(i,j)
                                a1=(surxix(i,j-1)+surxix(i,j)+surxix(i+1,j)+surxix(i+1,j-1))/4.
                                a2=(surxiy(i,j-1)+surxiy(i,j)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.

                                utemp(1:ngz-1)=(u(1:ngz-1, ii  , jj-1)+u(1:ngz-1, ii  , jj)    &
                                                +u(1:ngz-1, ii+1, jj-1)+u(1:ngz-1, ii+1, jj))/4.
                                vc(1:ngz-1, ii, jj)=a2*utemp(1:ngz-1)+b2*v(1:ngz-1, ii, jj)

                                !--- wc ---
                                i = ibmap_1(ii)
                                a3=swrz(i,j)
                                wc(1:ngz-1, ii, jj)=a3*w(1:ngz-1, ii, jj)
                        end do

                        !----------- End of Sec F ---------------

                end if

                !--------------------------------------------
                !                EXTEND
                !--------------------------------------------
                if (iblock ==  1 ) then
                        uc(1:ngz-1,0,0:nly) = 2.0*uc (1:ngz-1,1,0:nly) - uc(1:ngz-1,2,0:nly)
                        if (jblock ==  1 )  then
                                uc(1:ngz-1, 0 ,   0 ) = -uc(1:ngz-1, 0 , 1 )
                                wc(1:ngz  , 0 ,   0 ) = -wc(1:ngz  , 0 , 1 )
                        end if
                        if (jblock == npy)  &
                                wc(1:ngz  , 0 , nlyb) = 2.0*wc (1:ngz  ,1, nlyb) - wc(1:ngz  ,2, nlyb)
                end if

                if (iblock == npx) then
                        uc(1:ngz-1, nlxb+1, 0:nly) = 2.0*uc (1:ngz-1, nlxb, 0:nly) - uc(1:ngz-1, nlxb-1, 0:nly)
                        if (jblock ==  1 )  then
                                uc(1:ngz-1, nlxb+1,  0  ) = -uc(1:ngz-1,nlxb+1, 1 )
                                wc(1:ngz  , nlxb  ,  0  ) = -wc(1:ngz  ,nlxb  , 1 )
                        end if
                        if (jblock == npy)  &
                                wc(1:ngz  , nlxb, nlyb) = 2.0*wc (1:ngz  ,nlxb-1, nlyb) - wc(1:ngz  ,nlxb-2, nlyb)
                end if

                if (jblock ==  1 ) &
                        vc(1:ngz-1, :, 0 )=0.

                if (jblock == npy) &
                        vc(1:ngz-1, :, nlyb+1)=2.0*vc(1:ngz-1, :,  nlyb)-vc(1:ngz-1, :, nlyb-1 )
                
                !--------------------------------------------
                !        Periodicity in Z-direction
                !--------------------------------------------
                uc(  0  , :,:) = uc(ngz-1, :, :)
                uc( ngz , :,:) = uc(  1  , :, :)

                vc(  0  , :,:) = vc(ngz-1, :,:)
                vc( ngz , :,:) = vc(  1  , :,:)

                wc(  0  , :,:) = wc(ngz-1, :, :)
                !TAZ DANGER ZONE: wc(  1  , :,:) = 0.5*(wc(1, :,:)+wc(ngz, :,:))
                wc( ngz , :,:) = wc(1, :, :)
                wc(ngz+1, :,:) = wc(2, :, :)

  !------------------------------------------------
  ! Ensure overall mass conservation is satisfied
  !------------------------------------------------
  !T  call outflow_massFlowConsistency(Uin, Uout, Vin, Vout)

                !----------------- Block boundaries -------------------
                !T-- HERE IS THE FIX: IMPORTANT 02/08/2002)
                !T-- Exchange halos in x-direction , then in y-direction
                !T-- x-direction
                call updateBorder_lim(uc, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 3)
                call updateBorder_lim(vc, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 3)
                call updateBorder_lim(wc, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 2)

                !T-- y-direction
                call updateBorder_lim(uc, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 3)
                call updateBorder_lim(vc, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 3)
                call updateBorder_lim(wc, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 2)

        return
end


!=======================================================================
subroutine FluxBC()
        use boundaries

        real*8  :: vtempb(ngz-1), vtempe(ngz-1)
        real*8  :: utempb(ngz-1), utempe(ngz-1)
        real*8, dimension(2) :: Atmp

        !-----------------------------------------------------------------
        !                Inflow condition (u,v,w) West
        !-----------------------------------------------------------------
        if (iblock.eq.1) then
                !----------  U^{eta} (i.e. v) -----------
                do j=1,nlyb
                        jj = jbmap_1(j)
                        utempb=(uc(1:ngz-1, 1 ,j-1)+uc(1:ngz-1, 1 ,j)  &
                                +uc(1:ngz-1, 0 ,j-1)+uc(1:ngz-1, 0 ,j)) /4.
                        v(1:ngz-1, 0 ,j)=svetax(0 ,jj)*utempb+svetay(0 ,jj)*vc(1:ngz-1, 0 ,j)
                end do

                !----------  U^{z} (i.e. w) -----------
                do j=1,nlyb-1
                        jj = jbmap_1(j)
                        a3=swz(0,jj)
                        w(1:ngz, 0,j)=a3*wc(1:ngz, 0,j)
                end do

                !----------  U^{xi} (i.e. u) -----------
                do j=1,nlyb-1
                        jj = jbmap_1(j)
                        vtempb = ( vc(1:ngz-1, 0   ,j+1)+vc(1:ngz-1, 1 ,j+1)   &
                                +vc(1:ngz-1, 0   , j )+vc(1:ngz-1, 1 , j ))/4.
                        u(1:ngz-1,  1 ,j)=suxix( 1 ,jj)*uc(1:ngz-1,  1 , j)+suxiy( 1 ,jj)*vtempb
                enddo
        end  if

        !-----------------------------------------------------------------
        !                Outflow condition (u,v,w) East
        !-----------------------------------------------------------------
        if (iblock.eq.npx) then
                !----------  U^{eta} (i.e. v) -----------
                do j=1,nlyb
                        jj = jbmap_1(j)
                        utempe=(uc(1:ngz-1, nlxb+1,j-1)+uc(1:ngz-1, nlxb+1,j)   &
                                +uc(1:ngz-1, nlxb  ,j-1)+uc(1:ngz-1, nlxb  ,j))  /4.
                        v(1:ngz-1,nlxb,j)=svetax(ngx,jj)*utempe+svetay(ngx,jj)*vc(1:ngz-1,nlxb,j)
                end do

                !----------  U^{z} (i.e. w) -----------
                do j=1,nlyb-1
                        jj = jbmap_1(j)
                        a3=swz(ngx,jj)
                        w(1:ngz, nlxb,j)=a3*wc(1:ngz, nlxb,j)
                enddo

                !----------  U^{xi} (i.e. u) -----------
                do j=1,nlyb-1
                        jj = jbmap_1(j)
                        vtempe = ( vc(1:ngz-1, nlxb-1, j+1)+vc(1:ngz-1, nlxb, j+1)  &
                                +vc(1:ngz-1, nlxb-1,  j )+vc(1:ngz-1, nlxb,  j ))/4.
                        u(1:ngz-1, nlxb,j)=suxix(ngx,jj)*uc(1:ngz-1, nlxb, j)+suxiy(ngx,jj)*vtempe
                enddo
        end if

        !-----------------------------------------------------------------
        !                        No-Slip at Lower Wall
        !-----------------------------------------------------------------
        if (jblock.eq.1) then
                v(1:ngz-1, 1:nlx-1, 1) = 0.
                w(1:ngz  , 0:nlx, 0 ) = -w(1:ngz  , 0:nlx,  1  )
                u(1:ngz-1, 0:nlx, 0 ) = -u(1:ngz-1, 0:nlx,  1  )
        end if

        !-----------------------------------------------------------------
        !                Blasius Velocity at Top of Domain
        !                        OR NO-STRESS
        !-----------------------------------------------------------------
        if (jblock.eq.npy) then
                !----------------------------------------
                ! This section (Sec C) assumes knowing
                ! Cartesian BC and calculating Flux BC.
                ! (Use this section with Dirichlet top BC)
                !----------------------------------------
                !----------  U^{eta} (i.e. v) -----------
                do i=1,nlxb
                        ii = ibmap_1(i)
                        utempe(:)= (  uc(1:ngz-1, i, nlyb  ) + uc(1:ngz-1, i+1, nlyb  )  &
                                     +uc(1:ngz-1, i, nlyb-1) + uc(1:ngz-1, i+1, nlyb-1) )/4.
                        v(1:ngz-1,i,nlyb)=svetay(ii,ngy)*vc(1:ngz-1,i,nlyb)+svetax(ii,ngy)*utempe(:)
                end do
                !----------  U^{z} (i.e. w) -----------
                do i=1,nlxb-1
                        ii = ibmap_1(i)
                        a3=swz(ii,ngy)
                        w(1:ngz, i,nlyb)=a3*wc(1:ngz, i,nlyb)
                end do
                !----------  U^{xi} (i.e. u) -----------
                do i=1,nlxb
                        ii = ibmap_1(i)
                        vtempe(:)= ( vc(1:ngz-1, i-1,nlyb+1)+vc(1:ngz-1, i,nlyb+1)  &
                                        +vc(1:ngz-1, i-1,nlyb  )+vc(1:ngz-1, i,nlyb  ))/4.
                        u(1:ngz-1,i,nlyb)=suxix(ii,ngy)*uc(1:ngz-1,i,nlyb)+suxiy(ii,ngy)*vtempe(:)
                end do
                !----------- End of Sec C ---------------

                !----------------------------------------
                ! This section (Sec F) assumes knowing
                ! FLUX BC DIRECTLY (Not from Cartesian BC)
                ! (Use this section with No-Stress top BC)
                !----------------------------------------
                !--- NO Stress BC
                u(1:ngz-1,:, nlyb ) =  u(1:ngz-1,:,nlyb-1)
                w( 1:ngz ,:, nlyb ) =  w( 1:ngz ,:,nlyb-1)
                v(1:ngz-1,:, nlyb ) =  0.
                v(1:ngz-1,:,nlyb+1) = -v(1:ngz-1,:,nlyb-1)
                !----------- End of Sec F ---------------
        end if


        !--------------------------------------------
        !                EXTEND
        !--------------------------------------------
        if (iblock ==  1 ) then
                u(1:ngz-1,0,0:nly) = 2.0*u (1:ngz-1,1,0:nly) - u(1:ngz-1,2,0:nly)
                if (jblock ==  1 )  then
                        u(1:ngz-1, 0 ,   0 ) = -u(1:ngz-1, 0 , 1 )
                        w(1:ngz  , 0 ,   0 ) = -w(1:ngz  , 0 , 1 )
                end if
                if (jblock == npy)  &
                        w(1:ngz  , 0 , nlyb) = 2.0*w (1:ngz  ,1, nlyb) - w(1:ngz  ,2, nlyb)
        end if

        if (iblock == npx) then
                u(1:ngz-1, nlxb+1, 0:nly) = 2.0*u (1:ngz-1, nlxb, 0:nly) - u(1:ngz-1, nlxb-1, 0:nly)
                if (jblock ==  1 )  then
                        u(1:ngz-1, nlxb+1,  0  ) = -u(1:ngz-1,nlxb+1, 1 )
                        w(1:ngz  , nlxb  ,  0  ) = -w(1:ngz  ,nlxb  , 1 )
                end if
                if (jblock == npy)  &
                        w(1:ngz  , nlxb, nlyb) = 2.0*w (1:ngz  ,nlxb-1, nlyb) - w(1:ngz  ,nlxb-2, nlyb)
        end if

        if (jblock ==  1 ) &
                v(1:ngz-1, :, 0 )=0.

        if (jblock == npy) &
                v(1:ngz-1, :, nlyb+1)=2.0*v(1:ngz-1, :,  nlyb)-v(1:ngz-1, :, nlyb-1 )
        
        !--------------------------------------------


        !-----------------------------------------------------------------
        !                Periodicity in Spanwise Direction
        !-----------------------------------------------------------------
        !----------  U^{eta} (i.e. v) -----------
        v( 0 , :,:)=v(ngz-1, :,:)
        v(ngz, :,:)=v(  1  , :,:)

        !----------  U^{z} (i.e. w) -----------
        w(  0  , :,:) = w(ngz-1, :,:)
        !TAZ DANGER ZONE: w(  1  , :,:) = 0.5*(w(1, :,:)+w(ngz, :,:))
        w(ngz  , :,:) = w(1, :,:)
        w(ngz+1, :,:) = w(2, :,:)

        !----------  U^{xi} (i.e. u) -----------
        u( 0 , :,:) = u(ngz-1, :,:)
        u(ngz, :,:) = u(  1  , :,:)

        !===========================================================================
        !------------------------------------------------
        ! Ensure overall mass conservation is satisfied
        !------------------------------------------------
        bmassin=0.0
        bmassout=0.0
        if (jblock.eq.1) then
                do i=i1_T,i2_T
                        bmassin = bmassin + sum(v(1:ngzm,i,1))
                end do
        end if

        if (jblock.eq.npy) then
                do i=i1_T,i2_T
                        bmassin = bmassin - sum(v(1:ngzm,i,nlyb))
                end do
        end if

        bmass = bmassin
        call globalSum(bmass,1)
        if (irank.eq.1) print*,'verticle mass flux = ',bmass

        if (iblock.eq.1) then
                do j=j1_T,j2_T
                        bmassin = bmassin + sum(u(1:ngzm,1,j))
                end do
        end if

        if (iblock.eq.npx) then
                do j=j1_T,j2_T
                        bmassout = bmassout + sum(u(1:ngzm, nlxb, j))
                end do
        end if

        do j=j1_T,j2_T
        do i=i1_T,i2_T
                bmassin = bmassin+w(1,i,j)-w(ngz,i,j)
        end do
        end do

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
                        bmassout = bmassout+ sum(u(1:ngzm, nlxb, j))
                end do
        end if

        call globalSum(bmassout, 1)

        if (irank.eq.1) then
                write(6,*) 'After adjustment'
                write(6,*) 'Mass inflow =  ', bmassin
                write(6,*) 'Mass outflow =  ', bmassout
                write(6,*) 'mass corrected: inflow - outflow =  ',bmassout-bmassin
        end if

        !ccccccccccccc
        !   extend
        !ccccccccccccc
        if (iblock.eq.1) then
                u(1:ngzm, 0 , 0:nly)=2.0*u(1:ngzm, 1,  0:nly)-u(1:ngzm, 2,  0:nly)
        end if

        if (iblock.eq.npx) then
                u(1:ngzm, nlxb+1, 0:nly)=2.0*u(1:ngzm, nlxb, 0:nly)-u(1:ngzm, nlxb-1, 0:nly)
        end if

        !TAZ DANGER ZONE:  Xiahua has u-periodicity here:
        u( 0 , :,:)=u(ngz-1, :,:)
        u(ngz, :,:)=u(1    , :,:)

        !----------------- Block boundaries -------------------
        !T-- HERE IS THE FIX: IMPORTANT 02/08/2002)
        !T-- Exchange halos in x-direction , then in y-direction
        !T-- x-direction
        call updateBorder_lim(u, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 3)
        call updateBorder_lim(v, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 3)
        call updateBorder_lim(w, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 2)

        !T-- y-direction
        call updateBorder_lim(u, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 3)
        call updateBorder_lim(v, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 3)
        call updateBorder_lim(w, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 2)

        return
end

!=======================================================================

subroutine xzplaneBC(S, nj)
        use boundaries
        real S(nx_1,nz_1,nj)

        ! Inflow and outflow: zero-order extrapolation
        if (ibmin == imin) S(imap_1(imino),:,:) = S(imap_1(imin),:,:)
        if (ibmax == imax) S(imap_1(imaxo),:,:) = S(imap_1(imax),:,:)

        ! Block boundaries, including periodicity
        call updateBorder(S, nx_1,nz_1,nj, id_x, 1)
        call updateBorder(S, nx_1,nz_1,nj, id_z, 2)

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

        !T-- Periodicity in Z-direction
        wst(ngz,:,:) = wst(1,:,:)

        return
end

!=======================================================================
subroutine YsystemBC(A,B,C,R, iuvw)
        use boundaries
        real A(nx_1,nz_1,ny_1), B(nx_1,nz_1,ny_1), C(nx_1,nz_1,ny_1), &
             R(nx_1,nz_1,ny_1)

        select case (iuvw)
        case (id_u)
                ! Lower boundary
                if (jbmin == jmin) then
                        A(:,:, jmap_1(jmino)) = 0.
                        B(:,:, jmap_1(jmino)) = 1.
                        C(:,:, jmap_1(jmino)) = 0.
                        R(:,:, jmap_1(jmino)) = Vin(:,:, 1)
                end if

                ! Upper boundary
                if (jbmax == jmax) then
                        A(:,:, jmap_1(jmaxo)) = 0.
                        B(:,:, jmap_1(jmaxo)) = 1.
                        C(:,:, jmap_1(jmaxo)) = 0.
                        R(:,:, jmap_1(jmaxo)) = Vout(:,:, 1)
                end if

        case (id_v)
                ! Lower boundary
                if (jbmin == jmin) then
                        A(:,:, jmap_1(jmino)) = 0.
                        B(:,:, jmap_1(jmino)) = 1.
                        C(:,:, jmap_1(jmino)) = 0.
                        R(:,:, jmap_1(jmino)) = Vin(:,:, 2)
                end if

                ! Upper boundary
                if (jbmax == jmax) then
                        A(:,:, jmap_1(jmax)) = 0.
                        B(:,:, jmap_1(jmax)) = 1.
                        C(:,:, jmap_1(jmax)) = 0.
                        R(:,:, jmap_1(jmax)) = Vout(:,:, 2)
                end if

                ! Values for V at j=jmaxo are not used
                if (jbmax == jmax) then
                        A(:,:, jmap_1(jmaxo)) = 0.
                        B(:,:, jmap_1(jmaxo)) = 1.
                        C(:,:, jmap_1(jmaxo)) = 0.
                        R(:,:, jmap_1(jmaxo)) = 0.
                end if

        case (id_w)
                ! Lower boundary
                if (jbmin == jmin) then
                        A(:,:, jmap_1(jmino)) = 0.
                        B(:,:, jmap_1(jmino)) = 1.
                        C(:,:, jmap_1(jmino)) = 0.
                        R(:,:, jmap_1(jmino)) = Vin(:,:, 3)
                end if

                ! Upper boundary
                if (jbmax == jmax) then
                        A(:,:, jmap_1(jmaxo)) = 0.
                        B(:,:, jmap_1(jmaxo)) = 1.
                        C(:,:, jmap_1(jmaxo)) = 0.
                        R(:,:, jmap_1(jmaxo)) = Vout(:,:, 3)
                end if

        end select

        return
end

!=======================================================================
subroutine timeDependent_Inlet_BC(deltaT)
        use boundaries
        !--------------------------------------------
        !         Time-advance inflow BC
        !--------------------------------------------
        if (iblock == 1) then
                ! Superimpose fluctuations on mean inlet

                call laminarInflowBC()

                select case (inflowFluct)
                case (0)
                case (1:)
                        call fluct_Time(deltaT)
                end select

                select case (inflowFluct)
                case (0) ! No fluctuations
!                case (1) ! Random fluctuations
!                        call fluct_random(Uin)
!                case (2) ! OS discrete mode fluctuations (in blayer)
!                        call fluct_OSdiscrete(Uin)
!                case (3) ! OS continuous mode fluctuations (in freestream)
!                        call fluct_OScontinuous(Uin)
!                case (3) ! OS continuous mode fluctuations (in freestream)
!                        call fluct_OScontinuous(Uin)
!                case (4) ! OS discrete and continuous mode fluctuations
!                        call fluct_freestream(Uin)
!                case (5:6) ! OS 2D mode + OS 3D oblique modes
!                        call fluct_OS2D3D(Uin)
!                case (5:6) ! OS 2D & 3D modes
!                        call fluct_OScontinuous(Uin)
                end select

                ! Boundary conditions for Uin
                
                ! No-slip at lower wall
                if (jblock == 1) then
                        !--- B.C for (u,w) at inlet
                        Uin(0,:,1) = -Uin(1,:,1)
                        Uin(0,:,3) = -Uin(1,:,3)
                        !--- B.C for (v) at inlet
                        Uin(1,:,2) = 0.
                        Uin(0,:,2) = 0.
                end if

                !T------------------------------------------------
                !T   V at j=jmaxp1 not used
                !T  if (jbmax == jmax) Uin(jmap_1(jmaxo),:,2) = 0.
                !T------------------------------------------------
                if (jblock == npy) Uin(nlyb+1,:,2) = Uin(nlyb,:,2)

                ! Periodicity in z
                if (npz.eq.1) then
                        Uin(:,kmino,:) = Uin(:,kmax,:)
                        Uin(:,kmaxo,:) = Uin(:,kmin,:)
                endif

        end if

        return
end


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



subroutine laminarInflowBC()
        use boundaries
                Uin(:,:,1) = ubulk
                Uin(:,:,2) = 0.
                Uin(:,:,3) = 0.
        return
end


