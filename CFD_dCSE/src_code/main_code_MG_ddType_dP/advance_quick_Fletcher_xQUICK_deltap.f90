
module advance_module
        use data_export
        integer	:: flag_newdt
        integer	:: flag_x, flag_y, flag_z
        real	:: x_QUICK
end module

subroutine advance_init()
        use advance_module
        integer	:: iflag
        !----------------------------------------
        ! (flag_newdt = 1) ==> changed dt
        ! (falg_newdt = 0) ==> same old dt
        ! Used in advance.f90 to figure out (r1,r2)
        !----------------------------------------
        call archive_isDefined("flag_newdt",iflag)
        if (iflag == 1) then
                call readInt("flag_newdt", flag_newdt)
        else
                flag_newdt = 1
        endif

        flag_x = flag_newdt
        flag_y = flag_newdt
        flag_z = flag_newdt

        !----------------------------------------
        ! (x_QUICK = 1.) ==> Fletcher's QUICK q = 1.*(3/8)/3
        ! (x_QUICK = x.) ==> Fletcher's QUICK q = x.*(3/8)/3
        !----------------------------------------
        call archive_isDefined("x_QUICK",iflag)
        if (iflag == 1) then
                call readFloat("x_QUICK", x_QUICK)
        else
                x_QUICK = 1.
        endif

        return
end

subroutine advance_write()
        use advance_module

        if (flag_x.ne.flag_y  .or.  flag_y.ne.flag_z)  &
                stop "flags in advance_xyz are not consistent"
        flag_newdt = flag_x * flag_y * flag_z
        call writeInt("flag_newdt", flag_newdt)

        call writeFloat("x_QUICK", x_QUICK)
        return
end


subroutine advancex
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! sgi scalar version
! calculate the convection, explicit diffusion and
! implicit diffusion terms
! Admas-Bashforth explicit scheme for convection term
! and explicit diffusion terms is applied here, 
! Crank-Nicolson implicit treatment of viscous term
! has already be incorporated when the equation was derived
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! use data_export
        use advance_module


!     if(ntime-ntime_.eq.0) then
!     if(ntime-ntime_.eq.1) then
        if(ntime.eq.0  .or.  flag_x.eq.1) then
                r1=1.
                r2=0.
                  conx = 0.0
                flag_x = 0
        else
                r1=1.5
                r2=0.5
        end if

        !ccccccccccccccccccccc
        !   for u equation
        !ccccccccccccccccccccc
        do 62 jj=j1_u,j2_u
        do 62 ii=i1_u,i2_u
        do 61 kk=1,ngz-1
                i = ibmap_1(ii)
                j = jbmap_1(jj)
                k = kk
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! convection terms
                !-----------------------------------------------------------------------
                ! U^{xi}_{e}U^{m}_{e}S_{m,e}
                !     (first terms in RHS brakets of eq 3.30)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                tta=(u(k,ii,jj)+u(k,ii+1,jj))/2.
                ttb=(v(k,ii,jj)+v(k,ii,jj+1))/2.
                a3=(surxix(i,j)+surxix(i+1,j))/2.
                b1=(surxiy(i,j)+surxiy(i+1,j))/2.
                b2=(svretax(i,j)+svretax(i,j+1))/2.
                b3=(svretay(i,j)+svretay(i,j+1))/2.
                
                a1=tta*tta*a3+tta*ttb*b2 
                a2=tta*tta*b1+tta*ttb*b3 
                !ccccccccccccccccccccccccccccccc
                !   add QUICK correction terms
                !ccccccccccccccccccccccccccccccc
                phi_m1=u(k,ii-1,jj)
                phi_  =u(k,ii  ,jj)
                phi_p1=u(k,ii+1,jj)
                phi_p2=u(k,ii+2,jj)
                
                a1=a1+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a3
                a2=a2+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                         -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b1 

                phi_m1=(v(k,ii-1,jj)+v(k,ii-2,jj)+v(k,ii-1,jj+1)+v(k,ii-2,jj+1))/4.0 
                phi_  =(v(k,ii  ,jj)+v(k,ii-1,jj)+v(k,ii  ,jj+1)+v(k,ii-1,jj+1))/4.0
                phi_p1=(v(k,ii  ,jj)+v(k,ii+1,jj)+v(k,ii  ,jj+1)+v(k,ii+1,jj+1))/4.0
                
                if(ii.ne.nlxb-1) then
                          phi_p2=(v(k,ii+1,jj)+v(k,ii+2,jj)+v(k,ii+1,jj+1)+v(k,ii+2,jj+1))/4.0
                else
                          phi_p2=(v(k,ii+1,jj)+v(k,ii+1,jj+1))/2.0
                end if 

                a1=a1+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b2
                a2=a2+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b3

                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! -U^{xi}_{w}U^{m}_{w}S_{m,w}
                !     (2nd terms in RHS brakets of eq 3.30)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                tta=(u(k,ii,jj)+u(k,ii-1,jj))/2.
                ttb=(v(k,ii-1,jj)+v(k,ii-1,jj+1))/2.
                a3=(surxix(i,j)+surxix(i-1,j))/2.
                b1=(surxiy(i,j)+surxiy(i-1,j))/2.
                b2=(svretax(i-1,j)+svretax(i-1,j+1))/2.
                b3=(svretay(i-1,j)+svretay(i-1,j+1))/2.

                a1=a1-(tta*tta*a3+tta*ttb*b2)
                a2=a2-(tta*tta*b1+tta*ttb*b3)
                !ccccccccccccccccccccccccccccccc
                !   add QUICK correction terms
                !ccccccccccccccccccccccccccccccc
                phi_m1=u(k,ii-2,jj)
                phi_  =u(k,ii-1,jj)
                phi_p1=u(k,ii  ,jj)
                phi_p2=u(k,ii+1,jj)

                a1=a1-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a3
                a2=a2-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b1

                if(ii.ne.2) then
                phi_m1=(v(k,ii-2,jj)+v(k,ii-3,jj)+v(k,ii-2,jj+1)+v(k,ii-3,jj+1))/4.0
                else
                phi_m1=(v(k,ii-2,jj)+v(k,ii-2,jj+1))/2.0
                end if    

                phi_  =(v(k,ii-1,jj)+v(k,ii-2,jj)+v(k,ii-1,jj+1)+v(k,ii-2,jj+1))/4.0
                phi_p1=(v(k,ii-1,jj)+v(k,ii  ,jj)+v(k,ii-1,jj+1)+v(k,ii  ,jj+1))/4.0
                phi_p2=(v(k,ii  ,jj)+v(k,ii+1,jj)+v(k,ii  ,jj+1)+v(k,ii+1,jj+1))/4.0

                a1=a1-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b2
                a2=a2-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b3

                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! U^{eta}_{n}U^{m}_{n}S_{m,n}
                !     (3rd terms in RHS brakets of eq 3.30)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                tta=(v(k,ii,jj+1)+v(k,ii-1,jj+1))/2.
                ttb=(u(k,ii,jj)+u(k,ii,jj+1))/2.
                a3=(surxix(i,j)+surxix(i,j+1))/2.
                b1=(surxiy(i,j)+surxiy(i,j+1))/2.
                b2=(svretax(i,j+1)+svretax(i-1,j+1))/2.
                b3=(svretay(i,j+1)+svretay(i-1,j+1))/2.

                a1=a1+(tta*ttb*a3+tta*tta*b2)
                a2=a2+(tta*ttb*b1+tta*tta*b3)
                !ccccccccccccccccccccccccccccccc
                !  add QUICK correction terms
                !ccccccccccccccccccccccccccccccc
                phi_m1=u(k,ii,jj-1) 
                phi_  =u(k,ii,jj  )  
                phi_p1=u(k,ii,jj+1)      
                if(jj.ne.nlyb-1) then
                        phi_p2=u(k,ii,jj+2)
                else
                        phi_p2=u(k,ii,jj+1)
                end if

                a1=a1+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a3
                a2=a2+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b1

                phi_m1=(v(k,ii-1,jj  )+v(k,ii,jj  )+v(k,ii-1,jj-1)+v(k,ii,jj-1))/4.0
                phi_  =(v(k,ii-1,jj+1)+v(k,ii,jj+1)+v(k,ii-1,jj  )+v(k,ii,jj  ))/4.0
                phi_p1=(v(k,ii-1,jj+2)+v(k,ii,jj+2)+v(k,ii-1,jj+1)+v(k,ii,jj+1))/4.0

                if(jj.ne.nlyb-1) then
                        phi_p2=(v(k,ii-1,jj+3)+v(k,ii,jj+3)+v(k,ii-1,jj+2)+v(k,ii,jj+2))/4.0
                else
                        phi_p2=(v(k,ii-1,jj+2)+v(k,ii,jj+2))/2.0
                end if

                a1=a1+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b2
                a2=a2+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b3

                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! -U^{eta}_{s}U^{m}_{s}S_{m,s}
                !     (4th terms in RHS brakets of eq 3.30)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                tta=(v(k,ii,jj)+v(k,ii-1,jj))/2.
                ttb=(u(k,ii,jj)+u(k,ii,jj-1))/2.
                a3=(surxix(i,j)+surxix(i,j-1))/2.
                b1=(surxiy(i,j)+surxiy(i,j-1))/2.
                b2=(svretax(i,j)+svretax(i-1,j))/2.
                b3=(svretay(i,j)+svretay(i-1,j))/2.

                a1=a1-(tta*ttb*a3+tta*tta*b2)
                a2=a2-(tta*ttb*b1+tta*tta*b3)
                !ccccccccccccccccccccccccccccccc
                !  add QUICK correction terms
                !ccccccccccccccccccccccccccccccc
                if(jj.ne.1) then
                        phi_m1=u(k,ii,jj-2)
                else
                        phi_m1=u(k,ii,jj-1)
                end if
                
                phi_  =u(k,ii,jj-1)
                phi_p1=u(k,ii,jj  )
                phi_p2=u(k,ii,jj+1)

                a1=a1-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a3
                a2=a2-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b1
                
                if(jj.ne.1) then
                        phi_m1=(v(k,ii-1,jj-1)+v(k,ii,jj-1)+v(k,ii-1,jj-2)+v(k,ii,jj-2))/4.0
                else
                        phi_m1=(v(k,ii-1,jj-1)+v(k,ii,jj-1))/2.0
                end if
                
                phi_  =(v(k,ii-1,jj  )+v(k,ii,jj  )+v(k,ii-1,jj-1)+v(k,ii,jj-1))/4.0
                phi_p1=(v(k,ii-1,jj+1)+v(k,ii,jj+1)+v(k,ii-1,jj  )+v(k,ii,jj  ))/4.0
                phi_p2=(v(k,ii-1,jj+2)+v(k,ii,jj+2)+v(k,ii-1,jj+1)+v(k,ii,jj+1))/4.0
                
                a1=a1-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b2
                a2=a2-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
                        -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b3
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! U^{z}_{f}U^{m}_{f}S_{m,f}
                !     (5th terms in RHS brakets of eq 3.30)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                tta=(u(k,ii,jj)+u(k+1,ii,jj))/2.
                ttb=((v(k  ,ii-1,jj)+v(k  ,ii-1,jj+1)+v(k  ,ii,jj)+v(k  ,ii,jj+1))/4. &
                        +(v(k+1,ii-1,jj)+v(k+1,ii-1,jj+1)+v(k+1,ii,jj)+v(k+1,ii,jj+1))/4.)/2. 
                a3=surxix(i,j)
                b1=surxiy(i,j)
                b2=(svretax(i,j)+svretax(i,j+1)+svretax(i-1,j)+svretax(i-1,j+1))/4.
                b3=(svretay(i,j)+svretay(i,j+1)+svretay(i-1,j)+svretay(i-1,j+1))/4.

                a1=a1+(w(k+1,ii,jj)+w(k+1,ii-1,jj))/2.*(tta*a3+ttb*b2)
                a2=a2+(w(k+1,ii,jj)+w(k+1,ii-1,jj))/2.*(tta*b1+ttb*b3)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! -U^{z}_{b}U^{m}_{b}S_{m,b}
                !     (6th/Last terms in RHS brakets of eq 3.30)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                tta=(u(k,ii,jj)+u(k-1,ii,jj))/2.
                ttb=((v(k  ,ii-1,jj)+v(k  ,ii-1,jj+1)+v(k  ,ii,jj)+v(k  ,ii,jj+1))/4.  &
                            +(v(k-1,ii-1,jj)+v(k-1,ii-1,jj+1)+v(k-1,ii,jj)+v(k-1,ii,jj+1))/4.)/2.
                a3=surxix(i,j)
                b1=surxiy(i,j)
                b2=(svretax(i,j)+svretax(i,j+1)+svretax(i-1,j)+svretax(i-1,j+1))/4.
                b3=(svretay(i,j)+svretay(i,j+1)+svretay(i-1,j)+svretay(i-1,j+1))/4.

                a1=a1-(w(k,ii,jj)+w(k,ii-1,jj))/2.*(tta*a3+ttb*b2)
                a2=a2-(w(k,ii,jj)+w(k,ii-1,jj))/2.*(tta*b1+ttb*b3)
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! compute the total convection terms here  
                ! i.e.  Multiply   (a1,a2).(Suxix,Suxiy) ------> eq(3.30)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                con=-a1*suxix(i,j)-a2*suxiy(i,j) 


                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! explicit diffusion terms 
                !-----------------------------------------------------------------------
                ! d_{xi,e,ex}            (eq 3.35)
                !-----------------------------------------------------------------------
                !       (S_E^xi S_eta,E + S_eta,E S_E^xi) U_E^eta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=suxix(i+1,j)
                a2=suxiy(i+1,j)
                b1=(svretax(i,j)+svretax(i+1,j)+svretax(i,j+1)+svretax(i+1,j+1))/4.
                b2=(svretay(i,j)+svretay(i+1,j)+svretay(i,j+1)+svretay(i+1,j+1))/4.
                tta=(v(k,ii,jj)+v(k,ii+1,jj)+v(k,ii,jj+1)+v(k,ii+1,jj+1))/4.
                c1=(a1*b1)*tta 
                c2=(a1*b2)*tta 
                c4=(a2*b1)*tta
                c5=(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S_P^xi S_eta,P + S_eta,P S_P^xi) U_P^eta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=suxix(i,j)
                a2=suxiy(i,j)
                b1=(svretax(i,j)+svretax(i-1,j)+svretax(i,j+1)+svretax(i-1,j+1))/4.
                b2=(svretay(i,j)+svretay(i-1,j)+svretay(i,j+1)+svretay(i-1,j+1))/4.
                tta=(v(k,ii,jj)+v(k,ii-1,jj)+v(k,ii,jj+1)+v(k,ii-1,jj+1))/4.
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S_ne^eta S_xi,ne + S_xi,ne S_ne^eta) U_ne^xi
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=svetax(i,j+1)
                a2=svetay(i,j+1)
                b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j+1)+surxix(i+1,j+1))/4.
                b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j+1)+surxiy(i+1,j+1))/4.
                tta=(u(k,ii,jj)+u(k,ii+1,jj)+u(k,ii,jj+1)+u(k,ii+1,jj+1))/4.
                c1=c1+(a1*b1)*tta
                c2=c2+(a1*b2)*tta
                c4=c4+(a2*b1)*tta
                c5=c5+(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S_ne^eta S_eta,ne + S_eta,ne S_ne^eta) U_ne^eta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=svetax(i,j+1)
                a2=svetay(i,j+1)
                b1=svretax(i,j+1)
                b2=svretay(i,j+1)
                tta=v(k,ii,jj+1)
                c1=c1+(a1*b1)*tta
                c2=c2+(a1*b2)*tta
                c4=c4+(a2*b1)*tta
                c5=c5+(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S_se^eta S_xi,se + S_xi,se S_se^eta) U_se^xi
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=svetax(i,j)
                a2=svetay(i,j)
                b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j-1)+surxix(i+1,j-1))/4.
                b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j-1)+surxiy(i+1,j-1))/4.
                tta=(u(k,ii,jj)+u(k,ii+1,jj)+u(k,ii,jj-1)+u(k,ii+1,jj-1))/4.
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S_se^eta S_eta,se + S_eta,se S_se^eta) U_se^eta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=svetax(i,j)
                a2=svetay(i,j)
                b1=svretax(i,j)
                b2=svretay(i,j)
                tta=v(k,ii,jj)
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !     calculate S^{xi}_{e}/(ReV_{e}) \bu d_{xi,e,ex}
                !     (1th term of RHS braket in 3.41)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i+1,j))/2.
                a2=(suxiy(i,j)+suxiy(i+1,j))/2.
                a3=vp(i,j) 
                e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
                e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3 
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! d_{xi,w,ex}                   (eq 3.36)
                !-----------------------------------------------------------------------
                !       (S_P^xi S_eta,P + S_eta,P S_P^xi) U_P^eta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=suxix(i,j)
                a2=suxiy(i,j)
                b1=(svretax(i,j)+svretax(i-1,j)+svretax(i,j+1)+svretax(i-1,j+1))/4.
                b2=(svretay(i,j)+svretay(i-1,j)+svretay(i,j+1)+svretay(i-1,j+1))/4.
                tta=(v(k,ii,jj)+v(k,ii-1,jj)+v(k,ii,jj+1)+v(k,ii-1,jj+1))/4.
                c1=(a1*b1)*tta
                c2=(a1*b2)*tta
                c4=(a2*b1)*tta
                c5=(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S_W^xi S_eta,W + S_eta,W S_W^xi) U_W^eta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=suxix(i-1,j)
                a2=suxiy(i-1,j)
                b1=(svretax(i-1,j)+svretax(i-2,j)+svretax(i-1,j+1)+svretax(i-2,j+1))/4.
                b2=(svretay(i-1,j)+svretay(i-2,j)+svretay(i-1,j+1)+svretay(i-2,j+1))/4.
                tta=(v(k,ii-1,jj)+v(k,ii-2,jj)+v(k,ii-1,jj+1)+v(k,ii-2,jj+1))/4.
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S^{eta}_{nw}S_{xi,nw}+S_{xi,nw}S^{eta}_{nw})U^{xi}_{nw}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=svetax(i-1,j+1)
                a2=svetay(i-1,j+1)
                b1=(surxix(i,j)+surxix(i-1,j)+surxix(i,j+1)+surxix(i-1,j+1))/4.
                b2=(surxiy(i,j)+surxiy(i-1,j)+surxiy(i,j+1)+surxiy(i-1,j+1))/4.
                tta=(u(k,ii,jj)+u(k,ii-1,jj)+u(k,ii,jj+1)+u(k,ii-1,jj+1))/4.
                c1=c1+(a1*b1)*tta
                c2=c2+(a1*b2)*tta
                c4=c4+(a2*b1)*tta
                c5=c5+(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S^{eta}_{nw}S_{eta,nw}+S_{eta,nw}S^{eta}_{nw})U^{eta}_{nw}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=svetax(i-1,j+1)
                a2=svetay(i-1,j+1)
                b1=svretax(i-1,j+1)
                b2=svretay(i-1,j+1)
                tta=v(k,ii-1,jj+1)
                c1=c1+(a1*b1)*tta
                c2=c2+(a1*b2)*tta
                c4=c4+(a2*b1)*tta
                c5=c5+(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{eta}_{sw}S_{xi,sw}+S_{xi,sw}S^{eta}_{sw})U^{xi}_{sw}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=svetax(i-1,j)
                a2=svetay(i-1,j)
                b1=(surxix(i,j)+surxix(i-1,j)+surxix(i,j-1)+surxix(i-1,j-1))/4.
                b2=(surxiy(i,j)+surxiy(i-1,j)+surxiy(i,j-1)+surxiy(i-1,j-1))/4.
                tta=(u(k,ii,jj)+u(k,ii-1,jj)+u(k,ii,jj-1)+u(k,ii-1,jj-1))/4.
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{eta}_{sw}S_{eta,sw}+S_{eta,sw}S^{eta}_{sw})U^{eta}_{sw}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=svetax(i-1,j)
                a2=svetay(i-1,j)
                b1=svretax(i-1,j)
                b2=svretay(i-1,j)
                tta=v(k,ii-1,jj)
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !    calculate -S^{xi}_{w}/(ReV_{w}) \bu d_{xi,w,ex}
                !    (2nd term of RHS braket in 3.41)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i-1,j))/2.
                a2=(suxiy(i,j)+suxiy(i-1,j))/2.
                a3=vp(i-1,j) 
                e1=e1-(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
                e2=e2-(a1*(c2+c4)+a2*(c5+c5))*visc/a3

                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! d_{xi,n,ex}                     (eq 3.37)
                !-----------------------------------------------------------------------
                !       (S^{eta}_{N}S_{eta,N}+S_{eta,N}S^{eta}_{N})U^{eta}_{N}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j+1)+svetax(i-1,j+1)+svetax(i,j+2)+svetax(i-1,j+2))/4.
                a2=(svetay(i,j+1)+svetay(i-1,j+1)+svetay(i,j+2)+svetay(i-1,j+2))/4.
                b1=(svretax(i,j+1)+svretax(i-1,j+1)+svretax(i,j+2)+svretax(i-1,j+2))/4.
                b2=(svretay(i,j+1)+svretay(i-1,j+1)+svretay(i,j+2)+svretay(i-1,j+2))/4.
                tta=(v(k,ii,jj+1)+v(k,ii-1,jj+1)+v(k,ii,jj+2)+v(k,ii-1,jj+2))/4.
                c1=(a1*b1)*tta
                c2=(a1*b2)*tta
                c4=(a2*b1)*tta
                c5=(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{eta}_{P}S_{eta,P}+S_{eta,P}S^{eta}_{P})U^{eta}_{P}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j)+svetax(i,j+1)+svetax(i-1,j)+svetax(i-1,j+1))/4.
                a2=(svetay(i,j)+svetay(i,j+1)+svetay(i-1,j)+svetay(i-1,j+1))/4.
                b1=(svretax(i,j)+svretax(i,j+1)+svretax(i-1,j)+svretax(i-1,j+1))/4.
                b2=(svretay(i,j)+svretay(i,j+1)+svretay(i-1,j)+svretay(i-1,j+1))/4.
                tta=(v(k,ii,jj)+v(k,ii,jj+1)+v(k,ii-1,jj)+v(k,ii-1,jj+1))/4.
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S^{xi}_{ne}S_{xi,ne}+S_{xi,ne}S^{xi}_{ne})U^{xi}_{ne}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i+1,j)+suxix(i,j+1)+suxix(i+1,j+1))/4.
                a2=(suxiy(i,j)+suxiy(i+1,j)+suxiy(i,j+1)+suxiy(i+1,j+1))/4.
                b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j+1)+surxix(i+1,j+1))/4.
                b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j+1)+surxiy(i+1,j+1))/4.
                tta=(u(k,ii,jj)+u(k,ii+1,jj)+u(k,ii,jj+1)+u(k,ii+1,jj+1))/4.
                c1=c1+(a1*b1)*tta
                c2=c2+(a1*b2)*tta
                c4=c4+(a2*b1)*tta
                c5=c5+(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S^{xi}_{ne}S_{eta,ne}+S_{eta,ne}S^{xi}_{ne})U^{eta}_{ne}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i+1,j)+suxix(i,j+1)+suxix(i+1,j+1))/4.
                a2=(suxiy(i,j)+suxiy(i+1,j)+suxiy(i,j+1)+suxiy(i+1,j+1))/4.
                b1=svretax(i,j+1)
                b2=svretay(i,j+1)
                tta=v(k,ii,jj+1)
                c1=c1+(a1*b1)*tta
                c2=c2+(a1*b2)*tta
                c4=c4+(a2*b1)*tta
                c5=c5+(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{xi}_{nw}S_{xi,nw}+S_{xi,nw}S^{xi}_{nw})U^{xi}_{nw}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i,j+1)+suxix(i-1,j)+suxix(i-1,j+1))/4.
                a2=(suxiy(i,j)+suxiy(i,j+1)+suxiy(i-1,j)+suxiy(i-1,j+1))/4.
                b1=(surxix(i,j)+surxix(i,j+1)+surxix(i-1,j)+surxix(i-1,j+1))/4.
                b2=(surxiy(i,j)+surxiy(i,j+1)+surxiy(i-1,j)+surxiy(i-1,j+1))/4.
                tta=(u(k,ii,jj)+u(k,ii,jj+1)+u(k,ii-1,jj)+u(k,ii-1,jj+1))/4.
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{xi}_{nw}S_{eta,nw}+S_{eta,nw}S^{xi}_{nw})U^{eta}_{nw}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i,j+1)+suxix(i-1,j)+suxix(i-1,j+1))/4.
                a2=(suxiy(i,j)+suxiy(i,j+1)+suxiy(i-1,j)+suxiy(i-1,j+1))/4.
                b1=svretax(i-1,j+1)
                b2=svretay(i-1,j+1)
                tta=v(k,ii-1,jj+1)
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !     calculate S^{eta}_{n}/(ReV_{n}) \bu d_{xi,n,ex}
                !     (3rd term of RHS braket in 3.41)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j+1)+svetax(i-1,j+1))/2.
                a2=(svetay(i,j+1)+svetay(i-1,j+1))/2.
                a3=(vp(i,j)+vp(i,j+1)+vp(i-1,j)+vp(i-1,j+1))/4.
                e1=e1+(a1*(c1+c1)+a2*(c4+c2))*visc/a3
                e2=e2+(a1*(c2+c4)+a2*(c5+c5))*visc/a3 
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! d_{xi,s,ex}                       (eq 3.38)
                !-----------------------------------------------------------------------
                !       (S^{eta}_{P}S_{eta,P}+S_{eta,P}S^{eta}_{P})U^{eta}_{P}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j)+svetax(i,j+1)+svetax(i-1,j)+svetax(i-1,j+1))/4.
                a2=(svetay(i,j)+svetay(i,j+1)+svetay(i-1,j)+svetay(i-1,j+1))/4.
                b1=(svretax(i,j)+svretax(i,j+1)+svretax(i-1,j)+svretax(i-1,j+1))/4.
                b2=(svretay(i,j)+svretay(i,j+1)+svretay(i-1,j)+svretay(i-1,j+1))/4.
                tta=(v(k,ii,jj)+v(k,ii,jj+1)+v(k,ii-1,jj)+v(k,ii-1,jj+1))/4.
                c1=(a1*b1)*tta
                c2=(a1*b2)*tta
                c4=(a2*b1)*tta
                c5=(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{eta}_{S}S_{eta,S}+S_{eta,S}S^{eta}_{S})U^{eta}_{S}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j)+svetax(i,j-1)+svetax(i-1,j)+svetax(i-1,j-1))/4.
                a2=(svetay(i,j)+svetay(i,j-1)+svetay(i-1,j)+svetay(i-1,j-1))/4.
                b1=(svretax(i,j)+svretax(i,j-1)+svretax(i-1,j)+svretax(i-1,j-1))/4.
                b2=(svretay(i,j)+svretay(i,j-1)+svretay(i-1,j)+svretay(i-1,j-1))/4.
                tta=(v(k,ii,jj)+v(k,ii,jj-1)+v(k,ii-1,jj)+v(k,ii-1,jj-1))/4.
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S^{xi}_{se}S_{xi,se}+S_{xi,se}S^{xi}_{se})U^{xi}_{se}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i+1,j)+suxix(i,j-1)+suxix(i+1,j-1))/4.
                a2=(suxiy(i,j)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i+1,j-1))/4.
                b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j-1)+surxix(i+1,j-1))/4.
                b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j-1)+surxiy(i+1,j-1))/4.
                tta=(u(k,ii,jj)+u(k,ii+1,jj)+u(k,ii,jj-1)+u(k,ii+1,jj-1))/4.
                c1=c1+(a1*b1)*tta
                c2=c2+(a1*b2)*tta
                c4=c4+(a2*b1)*tta
                c5=c5+(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S^{xi}_{se}S_{eta,se}+S_{eta,se}S^{xi}_{se})U^{eta}_{se}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i+1,j)+suxix(i,j-1)+suxix(i+1,j-1))/4.
                a2=(suxiy(i,j)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i+1,j-1))/4.
                b1=svretax(i,j)
                b2=svretay(i,j)
                tta=v(k,ii,jj)
                c1=c1+(a1*b1)*tta
                c2=c2+(a1*b2)*tta
                c4=c4+(a2*b1)*tta
                c5=c5+(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{xi}_{sw}S_{xi,sw}+S_{xi,sw}S^{xi}_{sw})U^{xi}_{sw}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i,j-1)+suxix(i-1,j)+suxix(i-1,j-1))/4.
                a2=(suxiy(i,j)+suxiy(i,j-1)+suxiy(i-1,j)+suxiy(i-1,j-1))/4.
                b1=(surxix(i,j)+surxix(i,j-1)+surxix(i-1,j)+surxix(i-1,j-1))/4.
                b2=(surxiy(i,j)+surxiy(i,j-1)+surxiy(i-1,j)+surxiy(i-1,j-1))/4.
                tta=(u(k,ii,jj)+u(k,ii,jj-1)+u(k,ii-1,jj)+u(k,ii-1,jj-1))/4.
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{xi}_{sw}S_{eta,sw}+S_{eta,sw}S^{xi}_{sw})U^{eta}_{sw}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i,j-1)+suxix(i-1,j)+suxix(i-1,j-1))/4.
                a2=(suxiy(i,j)+suxiy(i,j-1)+suxiy(i-1,j)+suxiy(i-1,j-1))/4.
                b1=svretax(i-1,j)
                b2=svretay(i-1,j)
                tta=v(k,ii-1,jj)
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !     calculate - S^{eta}_{s}/(ReV_{s}) \bu d_{xi,s,ex}
                !     (4th term of RHS braket in 3.41)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j)+svetax(i-1,j))/2.
                a2=(svetay(i,j)+svetay(i-1,j))/2.
                a3=(vp(i,j)+vp(i,j-1)+vp(i-1,j)+vp(i-1,j-1))/4.
                e1=e1-(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
                e2=e2-(a1*(c2+c4)+a2*(c5+c5))*visc/a3
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! d_{xi,f,ex}                (eq 3.39)
                !-----------------------------------------------------------------------
                !       (S^{z}_{F}S_{eta,F}+S_{eta,F}S^{z}_{F})U^{eta}_{F}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a3=(swz(i,j)+swz(i-1,j))/2.
                b1=(svretax(i,j)+svretax(i,j+1)+svretax(i-1,j)+svretax(i-1,j+1))/4.
                b2=(svretay(i,j)+svretay(i,j+1)+svretay(i-1,j)+svretay(i-1,j+1))/4.
                tta=(v(k+1,ii,jj)+v(k+1,ii,jj+1)+v(k+1,ii-1,jj)+v(k+1,ii-1,jj+1))/4.
                c7=(a3*b1)*tta
                c8=(a3*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{z}_{P}S_{eta,P}+S_{eta,P}S^{z}_{P})U^{eta}_{P}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a3=(swz(i,j)+swz(i-1,j))/2.
                b1=(svretax(i,j)+svretax(i,j+1)+svretax(i-1,j)+svretax(i-1,j+1))/4.
                b2=(svretay(i,j)+svretay(i,j+1)+svretay(i-1,j)+svretay(i-1,j+1))/4.
                tta=(v(k,ii,jj)+v(k,ii,jj+1)+v(k,ii-1,jj)+v(k,ii-1,jj+1))/4.
                c7=c7-(a3*b1)*tta
                c8=c8-(a3*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S^{xi}_{fe}S_{z,fe}+S_{z,fe}S^{xi}_{fe})U^{z}_{fe}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i+1,j))/2.0
                a2=(suxiy(i,j)+suxiy(i+1,j))/2.0
                b3=swrz(i,j)
                tta=w(k+1,ii,jj)
                c3=(a1*b3)*tta
                c6=(a2*b3)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{xi}_{fw}S_{z,fw}+S_{z,fw}S^{xi}_{fw})U^{z}_{fw}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i-1,j))/2.0
                a2=(suxiy(i,j)+suxiy(i-1,j))/2.0
                b3=swrz(i-1,j)
                tta=w(k+1,ii-1,jj)
                c3=c3-(a1*b3)*tta
                c6=c6-(a2*b3)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S^{eta}_{fn}S_{z,fn}+S_{z,fn}S^{eta}_{fn})U^{z}_{fn}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j+1)+svetax(i-1,j+1))/2.0
                a2=(svetay(i,j+1)+svetay(i-1,j+1))/2.0
                b3=(swrz(i,j)+swrz(i,j+1)+swrz(i-1,j)+swrz(i-1,j+1))/4.
                tta=(w(k+1,ii,jj)+w(k+1,ii-1,jj)+w(k+1,ii,jj+1)+w(k+1,ii-1,jj+1))/4.
                c3=c3+(a1*b3)*tta
                c6=c6+(a2*b3)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{eta}_{fs}S_{z,fs}+S_{z,fs}S^{eta}_{fs})U^{z}_{fs}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j)+svetax(i-1,j))/2.0
                a2=(svetay(i,j)+svetay(i-1,j))/2.0
                b3=(swrz(i,j)+swrz(i,j-1)+swrz(i-1,j)+swrz(i-1,j-1))/4.
                tta=(w(k+1,ii,jj)+w(k+1,ii-1,jj)+w(k+1,ii,jj-1)+w(k+1,ii-1,jj-1))/4.
                c3=c3-(a1*b3)*tta
                c6=c6-(a2*b3)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !     calculate  S^{z}_{f}/(ReV_{f}) \bu d_{xi,f,ex}
                !     (5th term of RHS braket in 3.41)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(swz(i,j)+swz(i-1,j))/2.
                a2=(vp(i,j)+vp(i-1,j))/2.
                e1=e1+(a1*(c7+c3))*visc/a2 
                e2=e2+(a1*(c8+c6))*visc/a2 
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! d_{xi,b,ex}                    (eq 3.40)
                !-----------------------------------------------------------------------
                !       (S^{z}_{P}S_{eta,P}+S_{eta,P}S^{z}_{P})U^{eta,P}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a3=(swz(i,j)+swz(i-1,j))/2.
                b1=(svretax(i,j)+svretax(i,j+1)+svretax(i-1,j)+svretax(i-1,j+1))/4.
                b2=(svretay(i,j)+svretay(i,j+1)+svretay(i-1,j)+svretay(i-1,j+1))/4.
                tta=(v(k,ii,jj)+v(k,ii,jj+1)+v(k,ii-1,jj)+v(k,ii-1,jj+1))/4.
                c7=(a3*b1)*tta
                c8=(a3*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{z}_{B}S_{eta,B}+S_{eta,B}S^{z}_{B})U^{eta,B}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a3=(swz(i,j)+swz(i-1,j))/2.
                b1=(svretax(i,j)+svretax(i,j+1)+svretax(i-1,j)+svretax(i-1,j+1))/4.
                b2=(svretay(i,j)+svretay(i,j+1)+svretay(i-1,j)+svretay(i-1,j+1))/4.
                tta=(v(k-1,ii,jj)+v(k-1,ii,jj+1)+v(k-1,ii-1,jj)+v(k-1,ii-1,jj+1))/4.
                                c7=c7-(a3*b1)*tta
                c8=c8-(a3*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S^{xi}_{be}S_{z,be}+S_{z,be}S^{xi}_{be})U^{z}_{be}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i+1,j))/2.0
                a2=(suxiy(i,j)+suxiy(i+1,j))/2.0
                a3=0.0
                b1=0.
                b2=0.
                b3=swrz(i,j)
                tta=w(k,ii,jj)
                c3=(a1*b3)*tta
                c6=(a2*b3)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{xi}_{bw}S_{z,bw}+S_{z,bw}S^{xi}_{bw})U^{z}_{bw}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i-1,j))/2.0
                a2=(suxiy(i,j)+suxiy(i-1,j))/2.0
                b3=swrz(i-1,j)
                tta=w(k,ii-1,jj)
                c3=c3-(a1*b3)*tta
                c6=c6-(a2*b3)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       (S^{eta}_{bn}S_{z,bn}+S_{z,bn}S^{eta}_{bn})U^{z}_{bn}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j+1)+svetax(i-1,j+1))/2.0
                a2=(svetay(i,j+1)+svetay(i-1,j+1))/2.0
                b3=(swrz(i,j)+swrz(i,j+1)+swrz(i-1,j)+swrz(i-1,j+1))/4.
                tta=(w(k,ii,jj)+w(k,ii-1,jj)+w(k,ii,jj+1)+w(k,ii-1,jj+1))/4.
                c3=c3+(a1*b3)*tta
                c6=c6+(a2*b3)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !       -(S^{eta}_{bs}S_{z,bs}+S_{z,bs}S^{eta}_{bs})U^{z}_{bs}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j)+svetax(i-1,j))/2.0
                a2=(svetay(i,j)+svetay(i-1,j))/2.0
                b3=(swrz(i,j)+swrz(i,j-1)+swrz(i-1,j)+swrz(i-1,j-1))/4.
                tta=(w(k,ii,jj)+w(k,ii-1,jj)+w(k,ii,jj-1)+w(k,ii-1,jj-1))/4.
                c3=c3-(a1*b3)*tta
                c6=c6-(a2*b3)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !     calculate  S^{z}_{b}/(ReV_{b}) \bu d_{xi,b,ex}
                !     (6st term of RHS braket in 3.41)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(swz(i,j)+swz(i-1,j))/2.
                a2=(vp(i,j)+vp(i-1,j))/2.
                e1=e1-(a1*(c7+c3))*visc/a2 
                e2=e2-(a1*(c8+c6))*visc/a2 
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! d_{xi,f,im}                
                !-----------------------------------------------------------------------
                ! NOTE: Front and Back seem to be treated EXPLICITLY (Fourier Transform?????)
                !      (1st element in 5th line  eq 3.34)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a3=(swz(i,j)+swz(i-1,j))/2.
                b1=surxix(i,j)
                b2=surxiy(i,j)
                tta=u(k+1,ii,jj)
                c7=(a3*b1)*tta
                c8=(a3*b2)*tta
                !---------------------------------------------------
                !       (2nd element in 5th line  eq 3.34)
                !---------------------------------------------------
                a3=(swz(i,j)+swz(i-1,j))/2.
                b1=surxix(i,j)
                b2=surxiy(i,j)
                tta=u(k,ii,jj)
                c7=c7-(a3*b1)*tta
                c8=c8-(a3*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !    calculate  S^{z}_{f}/(ReV_{f}) \bu d_{xi,f,im}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(swz(i,j)+swz(i-1,j))/2.
                a2=(vp(i,j)+vp(i-1,j))/2.
                e1=e1+(a1*c7)*visc/a2
                e2=e2+(a1*c8)*visc/a2

                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! d_{xi,b,im}
                !----------------------------------------------------------------------
                ! NOTE:  Front and Back seem to be treated EXPLICITLY
                !     (1st element in 6th line  eq 3.34)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a3=(swz(i,j)+swz(i-1,j))/2.
                b1=surxix(i,j)
                b2=surxiy(i,j)
                tta=u(k,ii,jj)
                c7=(a3*b1)*tta
                c8=(a3*b2)*tta
                !---------------------------------------------------
                !       (2nd element in 6th line  eq 3.34)
                !---------------------------------------------------
                a3=(swz(i,j)+swz(i-1,j))/2.
                b1=surxix(i,j)
                b2=surxiy(i,j)
                tta=u(k-1,ii,jj)
                c7=c7-(a3*b1)*tta
                c8=c8-(a3*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !    calculate  S^{z}_{b}/(ReV_{b}) \bu d_{xi,b,im}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(swz(i,j)+swz(i-1,j))/2.
                a2=(vp(i,j)+vp(i-1,j))/2.
                e1=e1-(a1*c7)*visc/a2
                e2=e2-(a1*c8)*visc/a2
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! calculate the total explicit diffusion terms here
                ! i.e.  Multiply   (e1,e2).(Suxix,Suxiy) ------> eq(3.41)
                ! We also add the total explicit convection term ('con')
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                con = suxix(i,j)*e1 + suxiy(i,j)*e2 + con
                
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! implicit diffusion terms
                !----------------------------------------------------------------------
                ! d_{xi,e,im}                
                !----------------------------------------------------------------------
                !       (1st elament in line 1  eq 3.34)
                !       (S^{xi}_{E}S_{xi,E}+S_{xi,E}S^{xi}_{E})U^{xi}_{E}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=suxix(i+1,j)
                a2=suxiy(i+1,j)
                b1=surxix(i+1,j)
                b2=surxiy(i+1,j)
                tta=u(k,ii+1,jj)
                c1=(a1*b1)*tta
                c2=(a1*b2)*tta
                c4=(a2*b1)*tta
                c5=(a2*b2)*tta
                !---------------------------------------------------
                !       (2nd elament in line 1  eq 3.34)
                !---------------------------------------------------
                a1=suxix(i,j)
                a2=suxiy(i,j)
                b1=surxix(i,j)
                b2=surxiy(i,j)
                tta=u(k,ii,jj)
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !    calculate S^{xi}_{e}/(ReV_{e}) \bu d_{xi,e,im}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i+1,j))/2.
                a2=(suxiy(i,j)+suxiy(i+1,j))/2.
                a3=vp(i,j)
                e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
                e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! d_{xi,w,im}
                !----------------------------------------------------------------------
                !     (1st element of 2nd line    eq 3.34)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=suxix(i,j)
                a2=suxiy(i,j)
                b1=surxix(i,j)
                b2=surxiy(i,j)
                tta=u(k,ii,jj)
                c1=(a1*b1)*tta
                c2=(a1*b2)*tta
                c4=(a2*b1)*tta
                c5=(a2*b2)*tta
                !---------------------------------------------------
                !       (2nd elament in 2nd line     eq 3.34)
                !---------------------------------------------------
                a1=suxix(i-1,j)
                a2=suxiy(i-1,j)
                b1=surxix(i-1,j)
                b2=surxiy(i-1,j)
                tta=u(k,ii-1,jj)
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !    calculate -S^{xi}_{w}/(ReV_{w}) \bu d_{xi,w,im}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(suxix(i,j)+suxix(i-1,j))/2.
                a2=(suxiy(i,j)+suxiy(i-1,j))/2.
                a3=vp(i-1,j) 
                e1=e1-(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
                e2=e2-(a1*(c2+c4)+a2*(c5+c5))*visc/a3 
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! d_{xi,n,im}
                !----------------------------------------------------------------------
                !     (1st element of 3rd line   eq 3.34)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j+1)+svetax(i-1,j+1)+svetax(i,j+2)+svetax(i-1,j+2))/4.
                a2=(svetay(i,j+1)+svetay(i-1,j+1)+svetay(i,j+2)+svetay(i-1,j+2))/4.
                b1=surxix(i,j+1)
                b2=surxiy(i,j+1)
                tta=u(k,ii,jj+1)
                c1=(a1*b1)*tta
                c2=(a1*b2)*tta
                c4=(a2*b1)*tta
                c5=(a2*b2)*tta
                !---------------------------------------------------
                !       (2nd elament in 3rd line     eq 3.34)
                !---------------------------------------------------
                a1=(svetax(i,j)+svetax(i,j+1)+svetax(i-1,j)+svetax(i-1,j+1))/4.
                a2=(svetay(i,j)+svetay(i,j+1)+svetay(i-1,j)+svetay(i-1,j+1))/4.
                b1=surxix(i,j)
                b2=surxiy(i,j)
                tta=u(k,ii,jj)
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !    calculate S^{eta}_{n}/(ReV_{n}) \bu d_{xi,n,im}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j+1)+svetax(i-1,j+1))/2.
                a2=(svetay(i,j+1)+svetay(i-1,j+1))/2.
                a3=(vp(i,j)+vp(i,j+1)+vp(i-1,j)+vp(i-1,j+1))/4.
                e1=e1+(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
                e2=e2+(a1*(c2+c4)+a2*(c5+c5))*visc/a3 

                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! d_{xi,s,im}
                !----------------------------------------------------------------------
                !     (1st element in 4th line    eq 3.34)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j)+svetax(i,j+1)+svetax(i-1,j)+svetax(i-1,j+1))/4.
                a2=(svetay(i,j)+svetay(i,j+1)+svetay(i-1,j)+svetay(i-1,j+1))/4.
                b1=surxix(i,j)
                b2=surxiy(i,j)
                tta=u(k,ii,jj)
                c1=(a1*b1)*tta
                c2=(a1*b2)*tta
                c4=(a2*b1)*tta
                c5=(a2*b2)*tta
                !---------------------------------------------------
                !       (2nd elament in 4th line     eq 3.34)
                !---------------------------------------------------
                a1=(svetax(i,j)+svetax(i,j-1)+svetax(i-1,j)+svetax(i-1,j-1))/4.
                a2=(svetay(i,j)+svetay(i,j-1)+svetay(i-1,j)+svetay(i-1,j-1))/4.
                b1=surxix(i,j-1)
                b2=surxiy(i,j-1)
                tta=u(k,ii,jj-1)
                c1=c1-(a1*b1)*tta
                c2=c2-(a1*b2)*tta
                c4=c4-(a2*b1)*tta
                c5=c5-(a2*b2)*tta
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                !    calculate -S^{eta}_{s}/(ReV_{s}) \bu d_{xi,s,im}
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(svetax(i,j)+svetax(i-1,j))/2.
                a2=(svetay(i,j)+svetay(i-1,j))/2.
                a3=(vp(i,j)+vp(i,j-1)+vp(i-1,j)+vp(i-1,j-1))/4.
                e1=e1-(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
                e2=e2-(a1*(c2+c4)+a2*(c5+c5))*visc/a3 
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! compute the total implicit diffusion terms     (eq 3.34)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                vis=suxix(i,j)*e1+suxiy(i,j)*e2

                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! Calculate pressure term (-dt*dp/dx)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                jflag_start = 1/jj              ! Only equal to one when (jj=1)
                jflag_end   = jj/(nlyb-1)       ! Only equal to one when (jj=nlyb-1)


                pres = (1-jflag_start-jflag_end) * &
                        ( -(  (suxix(i,j)*0.5*(suxix (i,j  )+suxix (i+1,j  ))                 &
                              +suxiy(i,j)*0.5*(suxiy (i,j  )+suxiy (i+1,j  )))*p(k,ii  ,jj)   &
                            - (suxix(i,j)*0.5*(suxix (i,j  )+suxix (i-1,j  ))                 &
                              +suxiy(i,j)*0.5*(suxiy (i,j  )+suxiy (i-1,j  )))*p(k,ii-1,jj) ) &
                           /( (vp(i,j)+vp(i-1,j))/2. )                                        &
                          -(  (suxix(i,j)*0.5*(svetax(i,j+1)+svetax(i-1,j+1))                 &
                              +suxiy(i,j)*0.5*(svetay(i,j+1)+svetay(i-1,j+1)))                &
                             *(p(k,ii,jj)+p(k,ii-1,jj)+p(k,ii,jj+1)+p(k,ii-1,jj+1))/4.        &
                            - (suxix(i,j)*0.5*(svetax(i,j  )+svetax(i-1,j  ))                 &
                              +suxiy(i,j)*0.5*(svetay(i,j  )+svetay(i-1,j  )))                &
                             *(p(k,ii,jj)+p(k,ii-1,jj)+p(k,ii,jj-1)+p(k,ii-1,jj-1))/4. )      &
                           /( (vp(i,j)+vp(i-1,j))/2. ) )

                pres = pres + jflag_end   * &
                        ( -(  (suxix(i,j)*0.5*(suxix (i,j  )+suxix (i+1,j  ))                 &
                              +suxiy(i,j)*0.5*(suxiy (i,j  )+suxiy (i+1,j  )))*p(k,ii  ,jj)   &
                            - (suxix(i,j)*0.5*(suxix (i,j  )+suxix (i-1,j  ))                 &
                              +suxiy(i,j)*0.5*(suxiy (i,j  )+suxiy (i-1,j  )))*p(k,ii-1,jj) ) &
                           /( (vp(i,j)+vp(i-1,j))/2. ) )

                pres = pres + jflag_start * &
                        ( -(  (suxix(i,j)*0.5*(suxix (i,j  )+suxix (i+1,j  ))                 &
                              +suxiy(i,j)*0.5*(suxiy (i,j  )+suxiy (i+1,j  )))*p(k,ii  ,jj)   &
                            - (suxix(i,j)*0.5*(suxix (i,j  )+suxix (i-1,j  ))                 &
                              +suxiy(i,j)*0.5*(suxiy (i,j  )+suxiy (i-1,j  )))*p(k,ii-1,jj) ) &
                           /( (vp(i,j)+vp(i-1,j))/2. ) )
                
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! calculate right-hand-side vector for u equation
                ! (x-component of RHS of equations 3.68, 3.80, 3.115)
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                a1=(vp(i,j)+vp(i-1,j))/2.
                rhs(k,ii,jj)=dt/a1*(r1*con-r2*conx(k,ii,jj)+vis)+pres
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                conx(k,ii,jj)=con
                
61           continue
62          continue
        
        return 
end






subroutine advancey
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the convection, explicit diffusion and
! implicit diffusion terms
! Admas-Bashforth explicit scheme for convection term
! and explicit diffusion terms is applied here, 
! Crank-Nicolson implicit treatment of viscous term
! has already be incorporated when the equation was
! derived
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! use data_export
        use advance_module


!     if(ntime-ntime_.eq.0) then
!     if(ntime-ntime_.eq.1) then
      if(ntime.eq.0  .or.  flag_y.eq.1) then
        r1=1.
        r2=0.
        cony = 0.0
        flag_y = 0
      else
        r1=1.5
        r2=0.5
      end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! for v equation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 60 jj=j1_v,j2_v
      do 60 ii=i1_v,i2_v  
      do 60 kk=1,ngz-1
        i = ibmap_1(ii)
        j = jbmap_1(jj)
        k = kk

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! convection terms 
!-----------------------------------------------------------------------
! U^{e}_{xi}U^{m}_{e}S_{m,e}
!     (first terms in RHS brakets of eq 3.44)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(u(k,ii+1,jj)+u(k,ii+1,jj-1))/2.
        a3=(surxix(i+1,j)+surxix(i+1,j-1))/2.
        b1=(surxiy(i+1,j)+surxiy(i+1,j-1))/2.
        ttb=(v(k,ii,jj)+v(k,ii+1,jj))/2.
        b2=(svretax(i,j)+svretax(i+1,j))/2.
        b3=(svretay(i,j)+svretay(i+1,j))/2.

        a1=tta*tta*a3+tta*ttb*b2
        a2=tta*tta*b1+tta*ttb*b3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! add QUICK correction terms
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        phi_m1=(u(k,ii-1,jj)+u(k,ii-1,jj-1)+u(k,ii  ,jj)+u(k,ii  ,jj-1))/4.0
        phi_  =(u(k,ii  ,jj)+u(k,ii  ,jj-1)+u(k,ii+1,jj)+u(k,ii+1,jj-1))/4.0
        phi_p1=(u(k,ii+1,jj)+u(k,ii+1,jj-1)+u(k,ii+2,jj)+u(k,ii+2,jj-1))/4.0

        if(ii.ne.nlxb-1) then
          phi_p2=(u(k,ii+2,jj)+u(k,ii+2,jj-1)+u(k,ii+3,jj)+u(k,ii+3,jj-1))/4.0
        else
          phi_p2=(u(k,ii+2,jj)+u(k,ii+2,jj-1))/2.0
        end if

        a1=a1+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a3
        a2=a2+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b1

        phi_m1=v(k,ii-1,jj)
        phi_  =v(k,ii  ,jj)
        phi_p1=v(k,ii+1,jj)

        if(ii.ne.nlxb-1) then
          phi_p2=v(k,ii+2,jj)
        else
          phi_p2=v(k,ii+1,jj)
        end if

        a1=a1+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b2 
        a2=a2+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b3 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! -U^{w}_{xi}U^{m}_{w}S_{m,w}
!     (2nd term in RHS brakets of eq 3.44)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(u(k,ii,jj)+u(k,ii,jj-1))/2.
        a3=(surxix(i,j)+surxix(i,j-1))/2.
        b1=(surxiy(i,j)+surxiy(i,j-1))/2.
        ttb=(v(k,ii,jj)+v(k,ii-1,jj))/2.
        b2=(svretax(i,j)+svretax(i-1,j))/2.
        b3=(svretay(i,j)+svretay(i-1,j))/2.

        a1=a1-(tta*tta*a3+tta*ttb*b2)
        a2=a2-(tta*tta*b1+tta*ttb*b3)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! add QUICK correction terms
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(ii.ne.1) then
          phi_m1=(u(k,ii-2,jj)+u(k,ii-2,jj-1)+u(k,ii-1,jj)+u(k,ii-1,jj-1))/4.0
        else
          phi_m1=(u(k,ii-1,jj)+u(k,ii-1,jj-1))/2.0
        end if

        phi_  =(u(k,ii-1,jj)+u(k,ii-1,jj-1)+u(k,ii  ,jj)+u(k,ii  ,jj-1))/4.0
        phi_p1=(u(k,ii  ,jj)+u(k,ii  ,jj-1)+u(k,ii+1,jj)+u(k,ii+1,jj-1))/4.0
        phi_p2=(u(k,ii+1,jj)+u(k,ii+1,jj-1)+u(k,ii+2,jj)+u(k,ii+2,jj-1))/4.0

        a1=a1-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a3
        a2=a2-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b1

        if(ii.ne.1) then
          phi_m1=v(k,ii-2,jj)
        else
          phi_m1=v(k,ii-1,jj)
        end if

        phi_  =v(k,ii-1,jj)
        phi_p1=v(k,ii  ,jj)
        phi_p2=v(k,ii+1,jj)

        a1=a1-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b2
        a2=a2-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! U^{eta}_{n}U^{m}_{n}S_{m,n}
!     (3rd term in RHS brakets of eq 3.44)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(v(k,ii,jj)+v(k,ii,jj+1))/2.
        ttb=(u(k,ii,jj)+u(k,ii+1,jj))/2.
        a3=(surxix(i,j)+surxix(i+1,j))/2.
        b1=(surxiy(i,j)+surxiy(i+1,j))/2.
        b2=(svretax(i,j)+svretax(i,j+1))/2.
        b3=(svretay(i,j)+svretay(i,j+1))/2.

        a1=a1+(tta*ttb*a3+tta*tta*b2) 
        a2=a2+(tta*ttb*b1+tta*tta*b3)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! add QUICK correction terms
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        phi_m1=(u(k,ii,jj-1)+u(k,ii+1,jj-1)+u(k,ii,jj-2)+u(k,ii+1,jj-2))/4. 
        phi_  =(u(k,ii,jj  )+u(k,ii+1,jj  )+u(k,ii,jj-1)+u(k,ii+1,jj-1))/4.
        phi_p1=(u(k,ii,jj+1)+u(k,ii+1,jj+1)+u(k,ii,jj  )+u(k,ii+1,jj  ))/4.

        if(jj.ne.nlyb-1) then
          phi_p2=(u(k,ii,jj+2)+u(k,ii+1,jj+2)+u(k,ii,jj+1)+u(k,ii+1,jj+1))/4.
        else
          phi_p2=(u(k,ii,jj+1)+u(k,ii+1,jj+1))/2.0
        end if

        a1=a1+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a3
        a2=a2+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b1

        phi_m1=v(k,ii,jj-1)
        phi_  =v(k,ii,jj  )
        phi_p1=v(k,ii,jj+1)

        if(jj.ne.nlyb-1) then
          phi_p2=v(k,ii,jj+2)
        else
          phi_p2=v(k,ii,jj+1)
        end if

        a1=a1+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b2
        a2=a2+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! -U^{eta}_{s}U^{m}_{s}S_{m,s}
!     (4th term in RHS brakets of eq 3.44)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(v(k,ii,jj)+v(k,ii,jj-1))/2.
        ttb=(u(k,ii,jj-1)+u(k,ii+1,jj-1))/2.
        a3=(surxix(i,j-1)+surxix(i+1,j-1))/2.
        b1=(surxiy(i,j-1)+surxiy(i+1,j-1))/2.
        b2=(svretax(i,j)+svretax(i,j-1))/2.
        b3=(svretay(i,j)+svretay(i,j-1))/2.

        a1=a1-(tta*ttb*a3+tta*tta*b2)
        a2=a2-(tta*ttb*b1+tta*tta*b3)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! add QUICK correction terms
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(jj.ne.2) then
          phi_m1=(u(k,ii,jj-2)+u(k,ii+1,jj-2)+u(k,ii,jj-3)+u(k,ii+1,jj-3))/4.
        else
          phi_m1=(u(k,ii,jj-2)+u(k,ii+1,jj-2))/2.0
        end if

        phi_  =(u(k,ii,jj-1)+u(k,ii+1,jj-1)+u(k,ii,jj-2)+u(k,ii+1,jj-2))/4.
        phi_p1=(u(k,ii,jj  )+u(k,ii+1,jj  )+u(k,ii,jj-1)+u(k,ii+1,jj-1))/4.
        phi_p2=(u(k,ii,jj+1)+u(k,ii+1,jj+1)+u(k,ii,jj  )+u(k,ii+1,jj  ))/4.

        a1=a1-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a3
        a2=a2-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b1

        if(jj.ne.2) then
          phi_m1=v(k,ii,jj-2)
        else
          phi_m1=v(k,ii,jj-1)
        end if

        phi_  =v(k,ii,jj-1)
        phi_p1=v(k,ii,jj  )
        phi_p2=v(k,ii,jj+1)

        a1=a1-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b2
        a2=a2-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*b3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! U^{z}_{f}U^{m}_{f}S_{m,f}
!     (5th term in RHS brakets of eq 3.44)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(v(k,ii,jj)+v(k+1,ii,jj))/2.
        ttb=((u(k  ,ii,jj)+u(k  ,ii,jj-1)+u(k  ,ii+1,jj)+u(k  ,ii+1,jj-1))/4. &
            +(u(k+1,ii,jj)+u(k+1,ii,jj-1)+u(k+1,ii+1,jj)+u(k+1,ii+1,jj-1))/4.)/2. 
        b1=svretax(i,j)
        b2=svretay(i,j)
        b3=(surxix(i,j)+surxix(i,j-1)+surxix(i+1,j)+surxix(i+1,j-1))/4.
        a3=(surxiy(i,j)+surxiy(i,j-1)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.

        a1=a1+(w(k+1,ii,jj)+w(k+1,ii,jj-1))/2.*(ttb*b3+tta*b1)
        a2=a2+(w(k+1,ii,jj)+w(k+1,ii,jj-1))/2.*(ttb*a3+tta*b2)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! -U^{z}_{b}U^{m}_{b}S_{m,b}
!     (Last term in RHS brakets of eq 3.44)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(v(k,ii,jj)+v(k-1,ii,jj))/2.
        ttb=((u(k  ,ii,jj)+u(k  ,ii,jj-1)+u(k  ,ii+1,jj)+u(k  ,ii+1,jj-1))/4. &
            +(u(k-1,ii,jj)+u(k-1,ii,jj-1)+u(k-1,ii+1,jj)+u(k-1,ii+1,jj-1))/4.)/2.
        a3=svretax(i,j)
        b1=svretay(i,j)
        b2=(surxix(i,j)+surxix(i,j-1)+surxix(i+1,j)+surxix(i+1,j-1))/4.
        b3=(surxiy(i,j)+surxiy(i,j-1)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.

        a1=a1-(w(k,ii,jj)+w(k,ii,jj-1))/2.*(ttb*b2+tta*a3)
        a2=a2-(w(k,ii,jj)+w(k,ii,jj-1))/2.*(ttb*b3+tta*b1)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute the total convection terms here
! i.e.  Multiply   (a1,a2).(Svetax,Svetay) ------> eq(3.44)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        con=-a1*svetax(i,j)-a2*svetay(i,j) 


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! explicit diffusion terms
!-----------------------------------------------------------------------
! d_{eta,e,ex}            (eq 3.48)
!-----------------------------------------------------------------------
!       (S_{E}^{xi}S_{xi,E}+S_{xi,E}S^{xi}_{E})U_{E}^{xi}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i+1,j)+suxix(i+1,j-1)+suxix(i+2,j)+suxix(i+2,j-1))/4.
        a2=(suxiy(i+1,j)+suxiy(i+1,j-1)+suxiy(i+2,j)+suxiy(i+2,j-1))/4.
        b1=(surxix(i+1,j)+surxix(i+1,j-1)+surxix(i+2,j)+surxix(i+2,j-1))/4.
        b2=(surxiy(i+1,j)+surxiy(i+1,j-1)+surxiy(i+2,j)+surxiy(i+2,j-1))/4.
        tta=(u(k,ii+1,jj)+u(k,ii+1,jj-1)+u(k,ii+2,jj)+u(k,ii+2,jj-1))/4.
        c1=(a1*b1)*tta
        c2=(a1*b2)*tta
        c4=(a2*b1)*tta
        c5=(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S_{P}^{xi}S_{xi,P}+S_{xi,P}S^{xi}_{P})U_{P}^{xi}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i,j-1)+suxix(i+1,j)+suxix(i+1,j-1))/4.
        a2=(suxiy(i,j)+suxiy(i,j-1)+suxiy(i+1,j)+suxiy(i+1,j-1))/4.
        b1=(surxix(i,j)+surxix(i,j-1)+surxix(i+1,j)+surxix(i+1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i,j-1)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.
        tta=(u(k,ii,jj)+u(k,ii,jj-1)+u(k,ii+1,jj)+u(k,ii+1,jj-1))/4.
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_{ne}^{eta}S_{eta,ne}+S_{eta,ne}S^{eta}_{ne})U_{eta}^{ne}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i+1,j)+svetax(i,j+1)+svetax(i+1,j+1))/4.
        a2=(svetay(i,j)+svetay(i+1,j)+svetay(i,j+1)+svetay(i+1,j+1))/4.
        b1=(svretax(i,j)+svretax(i+1,j)+svretax(i,j+1)+svretax(i+1,j+1))/4.
        b2=(svretay(i,j)+svretay(i+1,j)+svretay(i,j+1)+svretay(i+1,j+1))/4.
        tta=(v(k,ii,jj)+v(k,ii+1,jj)+v(k,ii,jj+1)+v(k,ii+1,jj+1))/4.
        c1=c1+(a1*b1)*tta
        c2=c2+(a1*b2)*tta
        c4=c4+(a2*b1)*tta
        c5=c5+(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_{ne}^{eta}S_{xi,ne}+S_{xi,ne}S^{eta}_{ne})U_{xi}^{ne}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i+1,j)+svetax(i,j+1)+svetax(i+1,j+1))/4.
        a2=(svetay(i,j)+svetay(i+1,j)+svetay(i,j+1)+svetay(i+1,j+1))/4.
        b1=surxix(i+1,j)
        b2=surxiy(i+1,j)
        tta=u(k,ii+1,jj)
        c1=c1+(a1*b1)*tta
        c2=c2+(a1*b2)*tta
        c4=c4+(a2*b1)*tta
        c5=c5+(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S_{se}^{eta}S_{eta,se}+S_{eta,se}S^{eta}_{se})U_{eta}^{se}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i+1,j)+svetax(i,j-1)+svetax(i+1,j-1))/4.
        a2=(svetay(i,j)+svetay(i+1,j)+svetay(i,j-1)+svetay(i+1,j-1))/4.
        b1=(svretax(i,j)+svretax(i+1,j)+svretax(i,j-1)+svretax(i+1,j-1))/4.
        b2=(svretay(i,j)+svretay(i+1,j)+svretay(i,j-1)+svretay(i+1,j-1))/4.
        tta=(v(k,ii,jj)+v(k,ii+1,jj)+v(k,ii,jj-1)+v(k,ii+1,jj-1))/4.
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_{ne}^{eta}S_{xi,ne}+S_{xi,ne}S^{eta}_{ne})U_{xi}^{ne}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i+1,j)+svetax(i,j-1)+svetax(i+1,j-1))/4.
        a2=(svetay(i,j)+svetay(i+1,j)+svetay(i,j-1)+svetay(i+1,j-1))/4.
        b1=surxix(i+1,j-1)
        b2=surxiy(i+1,j-1)
        tta=u(k,ii+1,jj-1)
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate S^{xi}_{e}/(ReV_{e}) \bu d_{eta,e,ex}
!     (1th term of RHS braket in 3.54)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i+1,j)+suxix(i+1,j-1))/2.
        a2=(suxiy(i+1,j)+suxiy(i+1,j-1))/2.
        a3=(vp(i,j)+vp(i,j-1)+vp(i+1,j)+vp(i+1,j-1))/4. 
        e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
        e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{eta,w,ex}        eq(3.49)
!-----------------------------------------------------------------------
!       (S_{P}^{xi}S_{xi,P}+S_{xi,P}S^{xi}_{P})U_{P}^{xi}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i,j-1)+suxix(i+1,j)+suxix(i+1,j-1))/4.
        a2=(suxiy(i,j)+suxiy(i,j-1)+suxiy(i+1,j)+suxiy(i+1,j-1))/4.
        b1=(surxix(i,j)+surxix(i,j-1)+surxix(i+1,j)+surxix(i+1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i,j-1)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.
        tta=(u(k,ii,jj)+u(k,ii,jj-1)+u(k,ii+1,jj)+u(k,ii+1,jj-1))/4.
        c1=(a1*b1)*tta
        c2=(a1*b2)*tta
        c4=(a2*b1)*tta
        c5=(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S_{W}^{xi}S_{xi,W}+S_{xi,W}S^{xi}_{W})U_{W}^{xi}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i,j-1)+suxix(i-1,j)+suxix(i-1,j-1))/4.
        a2=(suxiy(i,j)+suxiy(i,j-1)+suxiy(i-1,j)+suxiy(i-1,j-1))/4.
        b1=(surxix(i,j)+surxix(i,j-1)+surxix(i-1,j)+surxix(i-1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i,j-1)+surxiy(i-1,j)+surxiy(i-1,j-1))/4.
        tta=(u(k,ii,jj)+u(k,ii,jj-1)+u(k,ii-1,jj)+u(k,ii-1,jj-1))/4.
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_{nw}^{eta}S_{eta,nw}+S_{eta,nw}S^{eta}_{nw})U_{nw}^{eta}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i-1,j)+svetax(i,j+1)+svetax(i-1,j+1))/4.
        a2=(svetay(i,j)+svetay(i-1,j)+svetay(i,j+1)+svetay(i-1,j+1))/4.
        b1=(svretax(i,j)+svretax(i-1,j)+svretax(i,j+1)+svretax(i-1,j+1))/4.
        b2=(svretay(i,j)+svretay(i-1,j)+svretay(i,j+1)+svretay(i-1,j+1))/4.
        tta=(v(k,ii,jj)+v(k,ii-1,jj)+v(k,ii,jj+1)+v(k,ii-1,jj+1))/4.
        c1=c1+(a1*b1)*tta
        c2=c2+(a1*b2)*tta
        c4=c4+(a2*b1)*tta
        c5=c5+(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_{nw}^{eta}S_{xi,nw}+S_{xi,nw}S^{eta}_{nw})U_{nw}^{xi}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i-1,j)+svetax(i,j+1)+svetax(i-1,j+1))/4.
        a2=(svetay(i,j)+svetay(i-1,j)+svetay(i,j+1)+svetay(i-1,j+1))/4.
        b1=surxix(i,j)
        b2=surxiy(i,j)
        tta=u(k,ii,jj)
        c1=c1+(a1*b1)*tta
        c2=c2+(a1*b2)*tta
        c4=c4+(a2*b1)*tta
        c5=c5+(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S_{sw}^{eta}S_{eta,sw}+S_{eta,sw}S^{eta}_{sw})U_{sw}^{eta}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i-1,j)+svetax(i,j-1)+svetax(i-1,j-1))/4.
        a2=(svetay(i,j)+svetay(i-1,j)+svetay(i,j-1)+svetay(i-1,j-1))/4.
        b1=(svretax(i,j)+svretax(i-1,j)+svretax(i,j-1)+svretax(i-1,j-1))/4.
        b2=(svretay(i,j)+svretay(i-1,j)+svretay(i,j-1)+svretay(i-1,j-1))/4.
        tta=(v(k,ii,jj)+v(k,ii-1,jj)+v(k,ii,jj-1)+v(k,ii-1,jj-1))/4.
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S_{sw}^{eta}S_{xi,sw}+S_{xi,sw}S^{eta}_{sw})U_{sw}^{xi}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i-1,j)+svetax(i,j-1)+svetax(i-1,j-1))/4.
        a2=(svetay(i,j)+svetay(i-1,j)+svetay(i,j-1)+svetay(i-1,j-1))/4.
        b1=surxix(i,j-1)
        b2=surxiy(i,j-1)
        tta=u(k,ii,jj-1)
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate -S^{xi}_{w}/(ReV_{w}) \bu d_{eta,w,ex}
!     (2nd term of RHS braket in 3.54)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i,j-1))/2.
        a2=(suxiy(i,j)+suxiy(i,j-1))/2.
        a3=(vp(i,j)+vp(i-1,j)+vp(i-1,j-1)+vp(i,j-1))/4. 
        e1=e1-(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
        e2=e2-(a1*(c2+c4)+a2*(c5+c5))*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{eta,n,ex}            eq(3.50)
!-----------------------------------------------------------------------
!       (S_N^eta S_xi,N + S_xi,N S_N^eta) U_N^xi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j+1)
        a2=svetay(i,j+1)
        b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j+1)+surxix(i+1,j+1))/4.
        b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j+1)+surxiy(i+1,j+1))/4.
        tta=(u(k,ii,jj)+u(k,ii+1,jj)+u(k,ii,jj+1)+u(k,ii+1,jj+1))/4.
        c1=(a1*b1)*tta
        c2=(a1*b2)*tta
        c4=(a2*b1)*tta
        c5=(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_P^eta S_xi,P + S_xi,P S_P^eta) U_P^xi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j)
        a2=svetay(i,j)
        b1=(surxix(i,j)+surxix(i,j-1)+surxix(i+1,j)+surxix(i+1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i,j-1)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.
        tta=(u(k,ii,jj)+u(k,ii,jj-1)+u(k,ii+1,jj)+u(k,ii+1,jj-1))/4.
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_ne^xi S_eta,ne + S_eta,ne S_ne^xi) U_ne^eta 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i+1,j)
        a2=suxiy(i+1,j)
        b1=(svretax(i,j)+svretax(i+1,j)+svretax(i,j+1)+svretax(i+1,j+1))/4.
        b2=(svretay(i,j)+svretay(i+1,j)+svretay(i,j+1)+svretay(i+1,j+1))/4.
        tta=(v(k,ii,jj)+v(k,ii+1,jj)+v(k,ii,jj+1)+v(k,ii+1,jj+1))/4.
        c1=c1+(a1*b1)*tta
        c2=c2+(a1*b2)*tta
        c4=c4+(a2*b1)*tta
        c5=c5+(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_ne^xi S_xi,ne + S_xi,ne S_ne^xi) U_ne^xi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i+1,j)        
        a2=suxiy(i+1,j)
        b1=surxix(i+1,j)
        b2=surxiy(i+1,j)
        tta=u(k,ii+1,jj)
        c1=c1+(a1*b1)*tta
        c2=c2+(a1*b2)*tta
        c4=c4+(a2*b1)*tta
        c5=c5+(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_nw^xi S_eta,nw + S_eta,nw S_nw^xi) U_nw^eta 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i,j)
        a2=suxiy(i,j)
        b1=(svretax(i,j)+svretax(i,j+1)+svretax(i-1,j)+svretax(i-1,j+1))/4.
        b2=(svretay(i,j)+svretay(i,j+1)+svretay(i-1,j)+svretay(i-1,j+1))/4.
        tta=(v(k,ii,jj)+v(k,ii,jj+1)+v(k,ii-1,jj)+v(k,ii-1,jj+1))/4.
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_nw^xi S_xi,nw + S_xi,nw S_nw^xi) U_nw^xi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i,j)
        a2=suxiy(i,j)
        b1=surxix(i,j)
        b2=surxiy(i,j)
        tta=u(k,ii,jj)
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate S^{eta}_{n}/(ReV_{n}) \bu d_{eta,n,ex}
!     (3rd term of RHS braket in 3.54)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j+1))/2.
        a2=(svetay(i,j)+svetay(i,j+1))/2.
        a3=vp(i,j)
        e1=e1+(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
        e2=e2+(a1*(c2+c4)+a2*(c5+c5))*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{eta,s,ex}          eq(3.51)
!-----------------------------------------------------------------------
!       (S_{P}^{eta} S_{xi,P} + S_{xi,P} S_{P}^{eta}) U_{P}^{xi}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j)
        a2=svetay(i,j)
        b1=(surxix(i,j)+surxix(i,j-1)+surxix(i+1,j)+surxix(i+1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i,j-1)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.
        tta=(u(k,ii,jj)+u(k,ii,jj-1)+u(k,ii+1,jj)+u(k,ii+1,jj-1))/4.
        c1=(a1*b1)*tta
        c2=(a1*b2)*tta
        c4=(a2*b1)*tta
        c5=(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S_{S}^{eta} S_{xi,S} + S_{xi,S} S_{S}^{eta}) U_{S}^{xi}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j-1)
        a2=svetay(i,j-1)
        b1=(surxix(i,j-1)+surxix(i+1,j-1)+surxix(i,j-2)+surxix(i+1,j-2))/4.
        b2=(surxiy(i,j-1)+surxiy(i+1,j-1)+surxiy(i,j-2)+surxiy(i+1,j-2))/4.
        tta=(u(k,ii,jj-1)+u(k,ii+1,jj-1)+u(k,ii,jj-2)+u(k,ii+1,jj-2))/4.
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{se}S_{eta,se}+S_{eta,se}S^{xi}_{se})U^{eta}_{se}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i+1,j-1)
        a2=suxiy(i+1,j-1)
        b1=(svretax(i,j)+svretax(i+1,j)+svretax(i,j-1)+svretax(i+1,j-1))/4.
        b2=(svretay(i,j)+svretay(i+1,j)+svretay(i,j-1)+svretay(i+1,j-1))/4.
        tta=(v(k,ii,jj)+v(k,ii+1,jj)+v(k,ii,jj-1)+v(k,ii+1,jj-1))/4.
        c1=c1+(a1*b1)*tta
        c2=c2+(a1*b2)*tta
        c4=c4+(a2*b1)*tta
        c5=c5+(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{se}S_{xi,se}+S_{xi,se}S^{xi}_{se})U^{xi}_{se}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i+1,j-1)        
        a2=suxiy(i+1,j-1)
        b1=surxix(i+1,j-1)
        b2=surxiy(i+1,j-1)
        tta=u(k,ii+1,jj-1)
        c1=c1+(a1*b1)*tta
        c2=c2+(a1*b2)*tta
        c4=c4+(a2*b1)*tta
        c5=c5+(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{xi}_{sw}S_{eta,sw}+S_{eta,sw}S^{xi}_{sw})U^{eta}_{sw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i,j-1)
        a2=suxiy(i,j-1)
        b1=(svretax(i,j)+svretax(i,j-1)+svretax(i-1,j)+svretax(i-1,j-1))/4.
        b2=(svretay(i,j)+svretay(i,j-1)+svretay(i-1,j)+svretay(i-1,j-1))/4.
        tta=(v(k,ii,jj)+v(k,ii,jj-1)+v(k,ii-1,jj)+v(k,ii-1,jj-1))/4.
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{xi}_{sw}S_{xi,sw}+S_{xi,sw}S^{xi}_{sw})U^{xi}_{sw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i,j-1)
        a2=suxiy(i,j-1)
        b1=surxix(i,j-1)
        b2=surxiy(i,j-1)
        tta=u(k,ii,jj-1)
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate - S^{eta}_{s}/(ReV_{s}) \bu d_{eta,s,ex}
!     (4th term of RHS braket in 3.54)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j-1))/2.
        a2=(svetay(i,j)+svetay(i,j-1))/2.
        a3=vp(i,j-1)
        e1=e1-(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
        e2=e2-(a1*(c2+c4)+a2*(c5+c5))*visc/a3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{eta,f,ex}         eq(3.52)
!-----------------------------------------------------------------------
!       (S^{z}_{F}S_{xi,F}+S_{xi,F}S^{z}_{F})U^{xi}_{F}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j-1)+surxix(i+1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j-1)+surxiy(i+1,j-1))/4.
        tta=(u(k+1,ii,jj)+u(k+1,ii+1,jj)+u(k+1,ii,jj-1)+u(k+1,ii+1,jj-1))/4.
        c7=(a3*b1)*tta
        c8=(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{z}_{P}S_{xi,P}+S_{xi,P}S^{z}_{P})U^{xi}_{P}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j-1)+surxix(i+1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j-1)+surxiy(i+1,j-1))/4.
        tta=(u(k,ii,jj)+u(k,ii+1,jj)+u(k,ii,jj-1)+u(k,ii+1,jj-1))/4.
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!SY Error in the comment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{fe}S_{z,fe}+S_{z,fe}S^{xi}_{fe})U^{z}_{fe}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! The above comment should be replaced by the following one.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{fn}S_{z,fn}+S_{z,fn}S^{eta}_{fn})U^{z}_{fn}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j+1))/2.0
        a2=(svetay(i,j)+svetay(i,j+1))/2.0
        a3=0.0
        b1=0.0
        b2=0.0 
        b3=swrz(i,j)
        tta=w(k+1,ii,jj)
        c3=(a1*b3)*tta
        c6=(a2*b3)*tta
!SY Error in the comment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{xi}_{fw}S_{z,fw}+S_{z,fw}S^{xi}_{fw})U^{z}_{fw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! The above comment should be replaced by the following one.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{eta}_{fs}S_{z,fs}+S_{z,fs}S^{eta}_{fs})U^{z}_{fs}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j-1))/2.0
        a2=(svetay(i,j)+svetay(i,j-1))/2.0
        a3=0.0
        b1=0.0 
        b2=0.0 
        b3=swrz(i,j-1)
        tta=w(k+1,ii,jj-1)
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!SY Error in the comment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{fn}S_{z,fn}+S_{z,fn}S^{eta}_{fn})U^{z}_{fn}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! The above comment should be replaced by the following one.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{fe}S_{z,fe}+S_{z,fe}S^{xi}_{fe})U^{z}_{fe}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i+1,j)+suxix(i+1,j-1))/2.0
        a2=(suxiy(i+1,j)+suxiy(i+1,j-1))/2.0
        b3=(swrz(i,j)+swrz(i+1,j)+swrz(i,j-1)+swrz(i+1,j-1))/4.
        tta=(w(k+1,ii,jj)+w(k+1,ii,jj-1)+w(k+1,ii+1,jj)+w(k+1,ii+1,jj-1))/4.
        c3=c3+(a1*b3)*tta
        c6=c6+(a2*b3)*tta
!SY Error in the comment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{eta}_{fs}S_{z,fs}+S_{z,fs}S^{eta}_{fs})U^{z}_{fs}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! The above comment should be replaced by the following one.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{xi}_{fw}S_{z,fw}+S_{z,fw}S^{xi}_{fw})U^{z}_{fw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i,j-1))/2.0
        a2=(suxiy(i,j)+suxiy(i,j-1))/2.0
        b3=(swrz(i,j)+swrz(i-1,j)+swrz(i,j-1)+swrz(i-1,j-1))/4.
        tta=(w(k+1,ii,jj)+w(k+1,ii,jj-1)+w(k+1,ii-1,jj)+w(k+1,ii-1,jj-1))/4.
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate  S^{z}_{f}/(ReV_{f}) \bu d_{eta,f,ex}
!     (5th term of RHS braket in 3.54)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(swz(i,j)+swz(i,j-1))/2.
        a2=(vp(i,j)+vp(i,j-1))/2.
        e1=e1+(a1*(c7+c3))*visc/a2 
        e2=e2+(a1*(c8+c6))*visc/a2 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{eta,b,ex}              eq(3.53)
!-----------------------------------------------------------------------
!       (S^{z}_{P}S_{xi,P}+S_{xi,P}S^{z}_{P})U^{xi}_{P}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j-1)+surxix(i+1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j-1)+surxiy(i+1,j-1))/4.
        tta=(u(k,ii,jj)+u(k,ii+1,jj)+u(k,ii,jj-1)+u(k,ii+1,jj-1))/4.
        c7=(a3*b1)*tta
        c8=(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{z}_{B}S_{xi,B}+S_{xi,B}S^{z}_{B})U^{xi}_{B}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=(surxix(i,j)+surxix(i+1,j)+surxix(i,j-1)+surxix(i+1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i+1,j)+surxiy(i,j-1)+surxiy(i+1,j-1))/4.
        tta=(u(k-1,ii,jj)+u(k-1,ii+1,jj)+u(k-1,ii,jj-1)+u(k-1,ii+1,jj-1))/4.
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!SY Error in the comment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{be}S_{z,be}+S_{z,be}S^{xi}_{be})U^{z}_{be}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! The above comment should be replaced by the following one.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{bn}S_{z,bn}+S_{z,bn}S^{eta}_{bn})U^{z}_{bn}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j+1))/2.0
        a2=(svetay(i,j)+svetay(i,j+1))/2.0
        b3=swrz(i,j)
        tta=w(k,ii,jj)
        c3=(a1*b3)*tta
        c6=(a2*b3)*tta
!SY Error in the comment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{xi}_{bw}S_{z,bw}+S_{z,bw}S^{xi}_{bw})U^{z}_{bw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! The above comment should be replaced by the following one.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{eta}_{bs}S_{z,bs}+S_{z,bs}S^{eta}_{bs})U^{z}_{bs}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j-1))/2.0
        a2=(svetay(i,j)+svetay(i,j-1))/2.0
        b3=swrz(i,j-1)
        tta=w(k,ii,jj-1)
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!SY Error in the comment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{bn}S_{z,bn}+S_{z,bn}S^{eta}_{bn})U^{z}_{bn}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! The above comment should be replaced by the following one.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{be}S_{z,be}+S_{z,be}S^{xi}_{be})U^{z}_{be}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i+1,j)+suxix(i+1,j-1))/2.0
        a2=(suxiy(i+1,j)+suxiy(i+1,j-1))/2.0
        b3=(swrz(i,j)+swrz(i+1,j)+swrz(i,j-1)+swrz(i+1,j-1))/4.
        tta=(w(k,ii,jj)+w(k,ii,jj-1)+w(k,ii+1,jj)+w(k,ii+1,jj-1))/4.
        c3=c3+(a1*b3)*tta
        c6=c6+(a2*b3)*tta
!SY Error in the comment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{eta}_{bs}S_{z,bs}+S_{z,bs}S^{eta}_{bs})U^{z}_{bs}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! The above comment should be replaced by the following one.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{xi}_{bw}S_{z,bw}+S_{z,bw}S^{xi}_{bw})U^{z}_{bw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i,j-1))/2.0
        a2=(suxiy(i,j)+suxiy(i,j-1))/2.0
        b3=(swrz(i,j)+swrz(i-1,j)+swrz(i,j-1)+swrz(i-1,j-1))/4.
        tta=(w(k,ii,jj)+w(k,ii,jj-1)+w(k,ii-1,jj)+w(k,ii-1,jj-1))/4.
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate  S^{z}_{b}/(ReV_{b}) \bu d_{eta,b,ex}
!     (Last term of RHS braket in 3.54)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(swz(i,j)+swz(i,j-1))/2.
        a2=(vp(i,j)+vp(i,j-1))/2.
        e1=e1-(a1*(c7+c3))*visc/a2 
        e2=e2-(a1*(c8+c6))*visc/a2 


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{eta,f,im}
!-----------------------------------------------------------------------
! NOTE: Front and Back seem to be treated EXPLICITLY (Fourier Transform?????)
!      (1st element in 5th line  eq 3.47)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=svretax(i,j)
        b2=svretay(i,j)
        tta=v(k+1,ii,jj)
        c7=(a3*b1)*tta
        c8=(a3*b2)*tta
!---------------------------------------------------
!       (2nd element in 5th line  eq 3.47)
!---------------------------------------------------
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=svretax(i,j)
        b2=svretay(i,j)
        tta=v(k,ii,jj)
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate  S^{z}_{f}/(ReV_{f}) \bu d_{eta,f,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(swz(i,j)+swz(i,j-1))/2.
        a2=(vp(i,j)+vp(i,j-1))/2.
        e1=e1+(a1*c7)*visc/a2
        e2=e2+(a1*c8)*visc/a2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{eta,b,im}
!-----------------------------------------------------------------------
! NOTE: Front and Back seem to be treated EXPLICITLY (Fourier Transform?????)
!      (1st element in 6th line  eq 3.47)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=svretax(i,j)
        b2=svretay(i,j)
        tta=v(k,ii,jj)
        c7=(a3*b1)*tta
        c8=(a3*b2)*tta
!---------------------------------------------------
!       (2nd element in 6th line  eq 3.47)
!---------------------------------------------------
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=svretax(i,j)
        b2=svretay(i,j)
        tta=v(k-1,ii,jj)
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate  S^{z}_{b}/(ReV_{b}) \bu d_{eta,b,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(swz(i,j)+swz(i,j-1))/2.
        a2=(vp(i,j)+vp(i,j-1))/2.
        e1=e1-(a1*c7)*visc/a2
        e2=e2-(a1*c8)*visc/a2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the total explicit diffusion terms here
! i.e.  Multiply   (e1,e2).(Svetax,Svetay) ------> eq(3.54)
! We also add the total explicit convection term ('con')
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        con = svetax(i,j)*e1 + svetay(i,j)*e2 + con


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! implicit diffusion terms
!----------------------------------------------------------------------
! d_{xi,e,im}
!----------------------------------------------------------------------
!       (1st elament in line 1  eq 3.47)
!       (S^{xi}_{E}S_{eta,E}+S_{eta,E}S^{xi}_{E})U^{eta}_{E}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i+1,j)+suxix(i+2,j)+suxix(i+1,j-1)+suxix(i+2,j-1))/4.
        a2=(suxiy(i+1,j)+suxiy(i+2,j)+suxiy(i+1,j-1)+suxiy(i+2,j-1))/4.
        b1=svretax(i+1,j)
        b2=svretay(i+1,j)
        tta=v(k,ii+1,jj)
        c1=(a1*b1)*tta
        c2=(a1*b2)*tta
        c4=(a2*b1)*tta
        c5=(a2*b2)*tta
!------------------------------------------------------------
!       (2nd elament in line 1  eq 3.47)
!       -(S^{xi}_{P}S_{eta,P}+S_{eta,P}S^{xi}_{P})U^{eta}_{P}
!------------------------------------------------------------
        a1=(suxix(i,j)+suxix(i+1,j)+suxix(i,j-1)+suxix(i+1,j-1))/4.
        a2=(suxiy(i,j)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i+1,j-1))/4.
        b1=svretax(i,j)
        b2=svretay(i,j)
        tta=v(k,ii,jj)
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate S^{xi}_{e}/(ReV_{e}) \bu d_{eta,e,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i+1,j)+suxix(i+1,j-1))/2.
        a2=(suxiy(i+1,j)+suxiy(i+1,j-1))/2.
        a3=(vp(i,j)+vp(i+1,j)+vp(i,j-1)+vp(i+1,j-1))/4. 
        e1=(a1*(c1+c1)+a2*(c4+c2))*visc/a3
        e2=(a1*(c2+c4)+a2*(c5+c5))*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{eta,w,im}
!----------------------------------------------------------------------
!       (1st elament in 2nd line of eq 3.47)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i+1,j)+suxix(i,j-1)+suxix(i+1,j-1))/4.
        a2=(suxiy(i,j)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i+1,j-1))/4.
        b1=svretax(i,j)
        b2=svretay(i,j)
        tta=v(k,ii,jj)
        c1=(a1*b1)*tta
        c2=(a1*b2)*tta
        c4=(a2*b1)*tta
        c5=(a2*b2)*tta
!---------------------------------------------
!       (2nd element in 2nd line of eq 3.47)
!---------------------------------------------
        a1=(suxix(i,j)+suxix(i,j-1)+suxix(i-1,j)+suxix(i-1,j-1))/4.
        a2=(suxiy(i,j)+suxiy(i,j-1)+suxiy(i-1,j)+suxiy(i-1,j-1))/4.
        b1=svretax(i-1,j)
        b2=svretay(i-1,j)
        tta=v(k,ii-1,jj)
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate -S^{xi}_{w}/(ReV_{w}) \bu d_{eta,w,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i,j-1))/2.
        a2=(suxiy(i,j)+suxiy(i,j-1))/2.
        a3=(vp(i,j)+vp(i-1,j)+vp(i,j-1)+vp(i-1,j-1))/4.
        e1=e1-(a1*(c1+c1)+a2*(c4+c2))*visc/a3
        e2=e2-(a1*(c2+c4)+a2*(c5+c5))*visc/a3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{eta,n,im}
!----------------------------------------------------------------------
!       (1st elament in 3rd line of eq 3.47)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j+1)
        a2=svetay(i,j+1)
        b1=svretax(i,j+1)
        b2=svretay(i,j+1)
        tta=v(k,ii,jj+1)
        c1=(a1*b1)*tta
        c2=(a1*b2)*tta
        c4=(a2*b1)*tta
        c5=(a2*b2)*tta
!-----------------------------------------------
!       (2nd element in 3rd line of eq 3.47)
!-----------------------------------------------
        a1=svetax(i,j)
        a2=svetay(i,j)
        b1=svretax(i,j)
        b2=svretay(i,j)
        tta=v(k,ii,jj)
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate S^{eta}_{n}/(ReV_{n}) \bu d_{eta,n,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j+1))/2.
        a2=(svetay(i,j)+svetay(i,j+1))/2.
        a3=vp(i,j)
        e1=e1+(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
        e2=e2+(a1*(c2+c4)+a2*(c5+c5))*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{eta,s,im}
!----------------------------------------------------------------------
!       (1st elament in 4th line of eq 3.47)
!       (S^{eta}_{P}S_{eta,P}+S_{eta,P}S^{eta}_{P})U^{eta}_{P}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j)
        a2=svetay(i,j)
        b1=svretax(i,j)
        b2=svretay(i,j)
        tta=v(k,ii,jj)
        c1=(a1*b1)*tta
        c2=(a1*b2)*tta
        c4=(a2*b1)*tta
        c5=(a2*b2)*tta
!-----------------------------------------------
!       (2nd element in 4th line of eq 3.47)
!-----------------------------------------------
        a1=svetax(i,j-1) 
        a2=svetay(i,j-1)
        b1=svretax(i,j-1)
        b2=svretay(i,j-1)
        tta=v(k,ii,jj-1)
        c1=c1-(a1*b1)*tta
        c2=c2-(a1*b2)*tta
        c4=c4-(a2*b1)*tta
        c5=c5-(a2*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate -S^{eta}_{s}/(ReV_{s}) \bu d_{eta,s,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j-1))/2.
        a2=(svetay(i,j)+svetay(i,j-1))/2.
        a3=vp(i,j-1)
        e1=e1-(a1*(c1+c1)+a2*(c4+c2))*visc/a3 
        e2=e2-(a1*(c2+c4)+a2*(c5+c5))*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute the total implicit diffusion terms       (eq 3.47)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        vis=svetax(i,j)*e1+svetay(i,j)*e2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Calculate EXPLICIT pressure term (-dt*dp/dy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        iflag_start = 1/ii              ! Only equal to one when (ii=1)
        iflag_end   = ii/(nlxb-1)       ! Only equal to one when (ii=nlxb-1)

        pres=   (1-iflag_start-iflag_end) * &
                 ( - ( (svetax(i,j)*0.5*(svetax(i  ,j)+svetax(i  ,j+1))                 &
                       +svetay(i,j)*0.5*(svetay(i  ,j)+svetay(i  ,j+1)))*p(k,ii,jj  )   &
                      -(svetax(i,j)*0.5*(svetax(i  ,j)+svetax(i  ,j-1))                 &
                       +svetay(i,j)*0.5*(svetay(i  ,j)+svetay(i  ,j-1)))*p(k,ii,jj-1) ) &
                     /( (vp(i,j)+vp(i,j-1))/2. )                                        &
                   - ( (svetax(i,j)*0.5*(suxix (i+1,j)+suxix (i+1,j-1))                 &
                       +svetay(i,j)*0.5*(suxiy (i+1,j)+suxiy (i+1,j-1)))                &
                      *(p(k,ii,jj)+p(k,ii+1,jj)+p(k,ii,jj-1)+p(k,ii+1,jj-1))/4.         &
                      -(svetax(i,j)*0.5*(suxix (i  ,j)+suxix (i  ,j-1))                 &
                       +svetay(i,j)*0.5*(suxiy (i  ,j)+suxiy (i  ,j-1)))                &
                      *(p(k,ii,jj)+p(k,ii,jj-1)+p(k,ii-1,jj)+p(k,ii-1,jj-1))/4. )       &
                     /( (vp(i,j)+vp(i,j-1))/2. ) )

        pres=   pres + iflag_end *   &
                 ( - ( (svetax(i,j)*0.5*(svetax(i  ,j)+svetax(i  ,j+1))                 &
                       +svetay(i,j)*0.5*(svetay(i  ,j)+svetay(i  ,j+1)))*p(k,ii,jj  )   &
                      -(svetax(i,j)*0.5*(svetax(i  ,j)+svetax(i  ,j-1))                 &
                       +svetay(i,j)*0.5*(svetay(i  ,j)+svetay(i  ,j-1)))*p(k,ii,jj-1) ) &
                     /( (vp(i,j)+vp(i,j-1))/2. ) )

        pres=   pres + iflag_start * &
                 ( - ( (svetax(i,j)*0.5*(svetax(i  ,j)+svetax(i  ,j+1))                 &
                       +svetay(i,j)*0.5*(svetay(i  ,j)+svetay(i  ,j+1)))*p(k,ii,jj  )   &
                      -(svetax(i,j)*0.5*(svetax(i  ,j)+svetax(i  ,j-1))                 &
                       +svetay(i,j)*0.5*(svetay(i  ,j)+svetay(i  ,j-1)))*p(k,ii,jj-1) ) &
                     /( (vp(i,j)+vp(i,j-1))/2. ) )

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate right-hand-side vector for v equation
! (y-component of RHS of equations 3.68, 3.80, 3.115)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(vp(i,j)+vp(i,j-1))/2.
        rhs(k,ii,jj)=dt/a1*(r1*con-r2*cony(k,ii,jj)+vis)+pres
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        cony(k,ii,jj)=con

 60   continue

      return 
end 




subroutine advancez
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the convection, explicit diffusion and
! implicit diffusion terms
! Admas-Bashforth explicit scheme for convection term
! and explicit diffusion terms is applied here, 
! Crank-Nicolson implicit treatment of viscous term
! has already be incorporated when the equation was derived
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! use data_export
        use advance_module


!     if(ntime-ntime_.eq.0) then
!     if(ntime-ntime_.eq.1) then
      if(ntime.eq.0  .or.  flag_z.eq.1) then
        r1=1.
        r2=0.
        conz = 0.0
        flag_z = 0
      else
        r1=1.5
        r2=0.5
      end if
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! for w equation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 63 jj=j1_w,j2_w
      do 63 ii=i1_w,i2_w
      do 63 kk=1,ngz-1
        i = ibmap_1(ii)
        j = jbmap_1(jj)
        k = kk
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! convection terms
!-----------------------------------------------------------------------
! U^{xi}_{e}U^{m}_{e}S_{m,e}
!     (first terms in RHS brakets of eq 3.57)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(u(k,ii+1,jj)+u(k-1,ii+1,jj))/2.
        ttb=(w(k,ii,jj)+w(k,ii+1,jj))/2.

        a1=(swrz(i,j)+swrz(i+1,j))/2.
        a3=tta*ttb*a1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! add QUICK correction terms
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        phi_m1=w(k,ii-1,jj) 
        phi_  =w(k,ii  ,jj)
        phi_p1=w(k,ii+1,jj)

        if(ii.ne.nlxb-1) then
          phi_p2=w(k,ii+2,jj)
        else
          phi_p2=w(k,ii+1,jj)
        end if

        a3=a3+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a1 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! -U^{xi}_{w}U^{m}_{w}S_{m,w}
!     (2nd terms in RHS brakets of eq 3.57)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(u(k,ii,jj)+u(k-1,ii,jj))/2.
        ttb=(w(k,ii,jj)+w(k,ii-1,jj))/2.

        a1=(swrz(i,j)+swrz(i-1,j))/2.
        a3=a3-tta*ttb*a1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! add QUICK correction terms
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(ii.ne.1) then
          phi_m1=w(k,ii-2,jj)
        else
          phi_m1=w(k,ii-1,jj)
        end if

        phi_  =w(k,ii-1,jj)
        phi_p1=w(k,ii  ,jj)
        phi_p2=w(k,ii+1,jj)

        a3=a3-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a1 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! U^{eta}_{n}U^{m}_{n}S_{m,n}
!     (3rd terms in RHS brakets of eq 3.57)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(v(k,ii,jj+1)+v(k-1,ii,jj+1))/2.
        ttb=(w(k,ii,jj)+w(k,ii,jj+1))/2.

        a1=(swrz(i,j)+swrz(i,j+1))/2.
        a3=a3+tta*ttb*a1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! add QUICK correction terms
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        phi_m1=w(k,ii,jj-1)
        phi_  =w(k,ii,jj  )
        phi_p1=w(k,ii,jj+1)

        if(jj.ne.nlyb-1) then
          phi_p2=w(k,ii,jj+2)
        else
          phi_p2=w(k,ii,jj+1)
        end if

        a3=a3+(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! -U^{eta}_{s}U^{m}_{s}S_{m,s}
!     (4th terms in RHS brakets of eq 3.57)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(v(k,ii,jj)+v(k-1,ii,jj))/2.
        ttb=(w(k,ii,jj)+w(k,ii,jj-1))/2.

        a1=(swrz(i,j)+swrz(i,j-1))/2.
        a3=a3-tta*ttb*a1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! add QUICK correction terms
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(jj.ne.1) then
          phi_m1=w(k,ii,jj-2)
        else
          phi_m1=w(k,ii,jj-1)
        end if

        phi_  =w(k,ii,jj-1)
        phi_p1=w(k,ii,jj  )
        phi_p2=w(k,ii,jj+1)

        a3=a3-(-(1.0+sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_m1-2.0*phi_  +phi_p1) &
               -(1.0-sign(1.0,tta))/2.0*x_QUICK*1.0/8.0*tta*(phi_  -2.0*phi_p1+phi_p2))*a1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! U^{z}_{f}U^{m}_{f}S_{m,f}
!     (5th terms in RHS brakets of eq 3.57)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(w(k,ii,jj)+w(k+1,ii,jj))/2.
        a1=swrz(i,j)
        a3=a3+tta*tta*a1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! -U^{z}_{b}U^{m}_{b}S_{m,b}
!     (Last terms in RHS brakets of eq 3.57)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        tta=(w(k,ii,jj)+w(k-1,ii,jj))/2.
        a1=swrz(i,j) 
        a3=a3-tta*tta*a1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute the total convection terms here
! i.e.  Multiply   (a3).(Swz) ------> eq(3.57)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        con=-a3*swz(i,j)


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! explicit diffusion terms 
!-----------------------------------------------------------------------
! d_{z,e,ex}             (eq 3.61)
!-----------------------------------------------------------------------
!      (S^{eta}_{ne}S_{z,ne}+S_{z,ne}S^{eta}_{ne})U^{ne}_{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j+1)+svetax(i+1,j+1))/2.
        a2=(svetay(i,j+1)+svetay(i+1,j+1))/2.
        b3=(swrz(i,j)+swrz(i+1,j)+swrz(i,j+1)+swrz(i+1,j+1))/4.
        tta=(w(k,ii,jj)+w(k,ii+1,jj)+w(k,ii,jj+1)+w(k,ii+1,jj+1))/4.
        c3=(a1*b3)*tta
        c6=(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      -(S^{eta}_{se}S_{z,se}+S_{z,se}S^{eta}_{se})U^{se}_{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i+1,j))/2.
        a2=(svetay(i,j)+svetay(i+1,j))/2.
        b3=(swrz(i,j)+swrz(i+1,j)+swrz(i,j-1)+swrz(i+1,j-1))/4.
        tta=(w(k,ii,jj)+w(k,ii+1,jj)+w(k,ii,jj-1)+w(k,ii+1,jj-1))/4.
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      (S^{z}_{fe}S_{xi,fe}+S_{xi,fe}S^{z}_{fe})U^{fe}_{xi}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i+1,j))/2.
        b1=surxix(i+1,j)
        b2=surxiy(i+1,j)
        tta=u(k,ii+1,jj)
        c7=(a3*b1)*tta
        c8=(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      -(S^{z}_{be}S_{xi,be}+S_{xi,be}S^{z}_{be})U^{be}_{xi}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i+1,j))/2.
        b1=surxix(i+1,j)
        b2=surxiy(i+1,j)
        tta=u(k-1,ii+1,jj)
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      (S^{z}_{fe}S_{eta,fe}+S_{eta,fe}S^{z}_{fe})U^{eta}_{fe}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i+1,j))/2.
        b1=(svretax(i,j)+svretax(i+1,j)+svretax(i,j+1)+svretax(i+1,j+1))/4.
        b2=(svretay(i,j)+svretay(i+1,j)+svretay(i,j+1)+svretay(i+1,j+1))/4.
        tta=(v(k,ii,jj)+v(k,ii+1,jj)+v(k,ii,jj+1)+v(k,ii+1,jj+1))/4.
        c7=c7+(a3*b1)*tta
        c8=c8+(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      -(S^{z}_{be}S_{eta,be}+S_{eta,be}S^{z}_{be})U^{eta}_{be}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i+1,j))/2.
        b1=(svretax(i,j)+svretax(i+1,j)+svretax(i,j+1)+svretax(i+1,j+1))/4.
        b2=(svretay(i,j)+svretay(i+1,j)+svretay(i,j+1)+svretay(i+1,j+1))/4.
        tta=(v(k-1,ii,jj)+v(k-1,ii+1,jj)+v(k-1,ii,jj+1)+v(k-1,ii+1,jj+1))/4.
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate S^{xi}_{e}/(Re_eV_e) \bu d_{z,e,ex}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i+1,j)
        a2=suxiy(i+1,j)
        a3=(vp(i,j)+vp(i+1,j))/2. 
        e3=(a1*(c3+c7)+a2*(c6+c8))*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{z,w,ex}             eq(3.62)
!-----------------------------------------------------------------------
!       (S^{eta}_{nw}S_{z,nw}+S_{z,nw}S^{eta}_{nw})U^{nw}_{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j+1)+svetax(i-1,j+1))/2.
        a2=(svetay(i,j+1)+svetay(i-1,j+1))/2.
        b3=(swrz(i,j)+swrz(i-1,j)+swrz(i,j+1)+swrz(i-1,j+1))/4.
        tta=(w(k,ii,jj)+w(k,ii-1,jj)+w(k,ii,jj+1)+w(k,ii-1,jj+1))/4. 
        c3=(a1*b3)*tta
        c6=(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{eta}_{sw}S_{z,sw}+S_{z,sw}S^{eta}_{sw})U^{sw}_{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i-1,j))/2.
        a2=(svetay(i,j)+svetay(i-1,j))/2.
        b3=(swrz(i,j)+swrz(i-1,j)+swrz(i,j-1)+swrz(i-1,j-1))/4.
        tta=(w(k,ii,jj)+w(k,ii-1,jj)+w(k,ii,jj-1)+w(k,ii-1,jj-1))/4.
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{fw}S_{xi,fw}+S_{xi,fw}S^{z}_{fw})U^{xi}_{fw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i-1,j))/2.
        b1=surxix(i,j)
        b2=surxiy(i,j)
        tta=u(k,ii,jj)
        c7=(a3*b1)*tta
        c8=(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{z}_{bw}S_{xi,bw}+S_{xi,bw}S^{z}_{bw})U^{xi}_{bw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i-1,j))/2.
        b1=surxix(i,j)
        b2=surxiy(i,j)
        tta=u(k-1,ii,jj)
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{fw}S_{eta,fw}+S_{eta,fw}S^{z}_{fw})U^{eta}_{fw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i-1,j))/2.
        b1=(svretax(i,j)+svretax(i-1,j)+svretax(i,j+1)+svretax(i-1,j+1))/4.
        b2=(svretay(i,j)+svretay(i-1,j)+svretay(i,j+1)+svretay(i-1,j+1))/4.
        tta=(v(k,ii,jj)+v(k,ii-1,jj)+v(k,ii,jj+1)+v(k,ii-1,jj+1))/4.
        c7=c7+(a3*b1)*tta
        c8=c8+(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{z}_{bw}S_{eta,bw}+S_{eta,bw}S^{z}_{bw})U^{eta}_{bw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i-1,j))/2.
        b1=(svretax(i,j)+svretax(i-1,j)+svretax(i,j+1)+svretax(i-1,j+1))/4.
        b2=(svretay(i,j)+svretay(i-1,j)+svretay(i,j+1)+svretay(i-1,j+1))/4.
        tta=(v(k-1,ii,jj)+v(k-1,ii-1,jj)+v(k-1,ii,jj+1)+v(k-1,ii-1,jj+1))/4.
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate S^{xi}_{w}/(ReV_{w}) \bu d_{z,w,ex}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i,j)
        a2=suxiy(i,j)
        a3=(vp(i,j)+vp(i-1,j))/2. 
        e3=e3-(a1*(c3+c7)+a2*(c6+c8))*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{z,n,ex}            eq(3.63)
!-----------------------------------------------------------------------
!       (S^{xi}_{ne}S_{z,ne}+S_{z,ne}S^{xi}_{ne}) U^{z}_{ne}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i+1,j)+suxix(i+1,j+1))/2.
        a2=(suxiy(i+1,j)+suxiy(i+1,j+1))/2.
        b3=(swrz(i,j)+swrz(i+1,j)+swrz(i,j+1)+swrz(i+1,j+1))/4.
        tta=(w(k,ii,jj)+w(k,ii+1,jj)+w(k,ii,jj+1)+w(k,ii+1,jj+1))/4. 
        c3=(a1*b3)*tta
        c6=(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{xi}_{nw}S_{z,nw}+S_{z,nw}S^{xi}_{nw}) U^{z}_{nw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i,j+1))/2.
        a2=(suxiy(i,j)+suxiy(i,j+1))/2.
        b3=(swrz(i,j)+swrz(i-1,j)+swrz(i,j+1)+swrz(i-1,j+1))/4.
        tta=(w(k,ii,jj)+w(k,ii-1,jj)+w(k,ii,jj+1)+w(k,ii-1,jj+1))/4. 
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{nf}S_{xi,nf}+S_{xi,nf}S^{z}_{nf}) U^{xi}_{nf}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j+1))/2.
        b1=(surxix(i,j)+surxix(i,j+1)+surxix(i+1,j)+surxix(i+1,j+1))/4.
        b2=(surxiy(i,j)+surxiy(i,j+1)+surxiy(i+1,j)+surxiy(i+1,j+1))/4.
        tta=(u(k,ii,jj)+u(k,ii+1,jj)+u(k,ii,jj+1)+u(k,ii+1,jj+1))/4.
        c7=(a3*b1)*tta
        c8=(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{z}_{nb}S_{xi,nb}+S_{xi,nb}S^{z}_{nb}) U^{xi}_{nb}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j+1))/2.
        b1=(surxix(i,j)+surxix(i,j+1)+surxix(i+1,j)+surxix(i+1,j+1))/4.
        b2=(surxiy(i,j)+surxiy(i,j+1)+surxiy(i+1,j)+surxiy(i+1,j+1))/4.
        tta=(u(k-1,ii,jj)+u(k-1,ii+1,jj)+u(k-1,ii,jj+1)+u(k-1,ii+1,jj+1))/4.
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{nf}S_{eta,nf}+S_{eta,nf}S^{z}_{nf}) U^{eta}_{nf}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j+1))/2.
        b1=svretax(i,j+1)
        b2=svretay(i,j+1)
        tta=v(k,ii,jj+1)
        c7=c7+(a3*b1)*tta
        c8=c8+(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{z}_{nb}S_{eta,nb}+S_{eta,nb}S^{z}_{nb}) U^{eta}_{nb}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j+1))/2.
        b1=svretax(i,j+1)
        b2=svretay(i,j+1)
        tta=v(k-1,ii,jj+1)
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate S^{eta}_{n}/(Re_nV_{n}) \bu d_{z,n,ex}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j+1)
        a2=svetay(i,j+1)
        a3=(vp(i,j)+vp(i,j+1))/2.
        e3=e3+(a1*(c3+c7)+a2*(c6+c8))*visc/a3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{z,s,ex}            eq(3.64)
!-----------------------------------------------------------------------
!       (S^{xi}_{se}S_{z,se}+S_{z,se}S^{xi}_{se}) U^{z}_{se}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i+1,j)+suxix(i+1,j-1))/2.
        a2=(suxiy(i+1,j)+suxiy(i+1,j-1))/2.
        b3=(swrz(i,j)+swrz(i+1,j)+swrz(i,j-1)+swrz(i+1,j-1))/4.
        tta=(w(k,ii,jj)+w(k,ii+1,jj)+w(k,ii,jj-1)+w(k,ii+1,jj-1))/4.
        c3=(a1*b3)*tta
        c6=(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{xi}_{sw}S_{z,sw}+S_{z,sw}S^{xi}_{sw}) U^{z}_{sw}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i,j-1))/2.
        a2=(suxiy(i,j)+suxiy(i,j-1))/2.
        b3=(swrz(i,j)+swrz(i-1,j)+swrz(i,j-1)+swrz(i-1,j-1))/4.
        tta=(w(k,ii,jj)+w(k,ii-1,jj)+w(k,ii,jj-1)+w(k,ii-1,jj-1))/4.
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{sf}S_{xi,sf}+S_{xi,sf}S^{z}_{sf}) U^{xi}_{sf}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=(surxix(i,j)+surxix(i,j-1)+surxix(i+1,j)+surxix(i+1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i,j-1)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.
        tta=(u(k,ii,jj)+u(k,ii+1,jj)+u(k,ii,jj-1)+u(k,ii+1,jj-1))/4.
        c7=(a3*b1)*tta
        c8=(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{z}_{sb}S_{xi,sb}+S_{xi,sb}S^{z}_{sb}) U^{xi}_{sb}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=(surxix(i,j)+surxix(i,j-1)+surxix(i+1,j)+surxix(i+1,j-1))/4.
        b2=(surxiy(i,j)+surxiy(i,j-1)+surxiy(i+1,j)+surxiy(i+1,j-1))/4.
        tta=(u(k-1,ii,jj)+u(k-1,ii+1,jj)+u(k-1,ii,jj-1)+u(k-1,ii+1,jj-1))/4.
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{z}_{sf}S_{eta,sf}+S_{eta,sf}S^{z}_{sf})U^{eta}_{sf}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=svretax(i,j)
        b2=svretay(i,j)
        tta=v(k,ii,jj)
        c7=c7+(a3*b1)*tta
        c8=c8+(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       -(S^{z}_{sb}S_{eta,sb}+S_{eta,sb}S^{z}_{sb})U^{eta}_{sb}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=(swz(i,j)+swz(i,j-1))/2.
        b1=svretax(i,j)
        b2=svretay(i,j)
        tta=v(k-1,ii,jj)
        c7=c7-(a3*b1)*tta
        c8=c8-(a3*b2)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate - S^{eta}_{s}/(ReV_{s}) \bu d_{z,s,ex}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j)
        a2=svetay(i,j)
        a3=(vp(i,j)+vp(i,j-1))/2.
        e3=e3-(a1*(c3+c7)+a2*(c6+c8))*visc/a3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{z,f,im} 
!-----------------------------------------------------------------------
!       (1st element in 5th line  eq 3.60)
!       (S^{z}_{F}S_{z,F}+S_{z,F}S^{z}_{F}) U_{z}^{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=swz(i,j)
        b3=swrz(i,j)
        tta=w(k+1,ii,jj)
        c9=(a3*b3)*tta
!-----------------------------------------------
!       (2nd element in 5th line  eq 3.60)
!-----------------------------------------------
        a3=swz(i,j)
        b3=swrz(i,j)
        tta=w(k,ii,jj)
        c9=c9-(a3*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate  S^{z}_{f}/(ReV_{f}) \bu d_{z,f,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=swz(i,j)
        a3=vp(i,j)
        e3=e3+a1*(c9+c9)*visc/a3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{z,b,im}
!-----------------------------------------------------------------------
!       (1st element in 6th line of eq 3.60)
!       -(S^{z}_{B}S_{z,B}+S_{z,B}S^{z}_{B}) U_{B}^{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a3=swz(i,j)
        b3=swrz(i,j)
        tta=w(k,ii,jj)
        c9=(a3*b3)*tta
!-----------------------------------------------
!       (2nd element in 6th line  eq 3.60)
!-----------------------------------------------
        a3=swz(i,j)
        b3=swrz(i,j)
        tta=w(k-1,ii,jj)
        c9=c9-(a3*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate  -S^{z}_{b}/(ReV_{b}) \bu d_{z,b,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=swz(i,j)
        a3=vp(i,j)
        e3=e3-a1*(c9+c9)*visc/a3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the total explicit diffusion terms here
! i.e.  Multiply   (e3).(Swz) ------> eq(3.67)
! We also add the total explicit convection term ('con')
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        con = swz(i,j)*e3 + con


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! implicit diffusion terms
!----------------------------------------------------------------------
! d_{z,e,im}
!----------------------------------------------------------------------
!       (1st elament in line 1  eq 3.60)
!       (S^{xi}_{E}S_{z,E}+S_{z,E}S^{xi}_{E}) U_{E}^{z} 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i+1,j)+suxix(i+2,j))/2.
        a2=(suxiy(i+1,j)+suxiy(i+2,j))/2.
        b3=swrz(i+1,j)
        tta=w(k,ii+1,jj)
        c3=(a1*b3)*tta
        c6=(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (2nd elament in line 1  eq 3.60)
!       -(S^{xi}_{P}S_{z,P}+S_{z,P}S^{xi}_{P}) U_{P}^{z}      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i+1,j))/2.
        a2=(suxiy(i,j)+suxiy(i+1,j))/2.
        b3=swrz(i,j)
        tta=w(k,ii,jj)
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate S^{xi}_{e}/(ReV_{e}) \bu d_{z,e,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i+1,j)
        a2=suxiy(i+1,j)
        a3=(vp(i,j)+vp(i+1,j))/2. 
        e3=(a1*(c3)+a2*(c6))*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{z,w,im}
!----------------------------------------------------------------------
!       (1st elament in 2nd line of  eq 3.60)
!       (S^{xi}_{P}S_{z,P}+S_{z,P}S^{xi}_{P})U^{z}_{P}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i+1,j))/2.
        a2=(suxiy(i,j)+suxiy(i+1,j))/2.
        b3=swrz(i,j)
        tta=w(k,ii,jj)
        c3=(a1*b3)*tta
        c6=(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (2nd elament in 2nd line of  eq 3.60)
!       -(S^{xi}_{W}S_{z,W}+S_{z,W}S^{xi}_{W}) U_{W}^{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i-1,j))/2.
        a2=(suxiy(i,j)+suxiy(i-1,j))/2.
        b3=swrz(i-1,j)
        tta=w(k,ii-1,jj)
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate -S^{xi}_{w}/(ReV_{w}) \bu d_{z,w,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i,j)
        a2=suxiy(i,j)
        a3=(vp(i,j)+vp(i-1,j))/2.
        e3=e3-(a1*c3+a2*c6)*visc/a3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{z,n,im}
!----------------------------------------------------------------------
!       (1st elament in 3rd line of  eq 3.60)
!       (S^{eta}_{N}S_{z,N}+S_{z,N}S^{eta}_{N}) U_{N}^{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j+1)+svetax(i,j+2))/2.
        a2=(svetay(i,j+1)+svetay(i,j+2))/2.
        b3=swrz(i,j+1)
        tta=w(k,ii,jj+1)
        c3=(a1*b3)*tta
        c6=(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (2nd elament in 3rd line of  eq 3.60)
!       -(S^{eta}_{P}S_{z,P}+S_{z,P}S^{eta}_{P}) U_{P}^{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j+1))/2.
        a2=(svetay(i,j)+svetay(i,j+1))/2.
        b3=swrz(i,j)
        tta=w(k,ii,jj)
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate S^{eta}_{n}/(ReV_{n}) \bu d_{z,n,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j+1)
        a2=svetay(i,j+1)
        a3=(vp(i,j)+vp(i,j+1))/2. 
        e3=e3+(a1*c3+a2*c6)*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! d_{z,s,im}
!----------------------------------------------------------------------
!       (1st elament in 4th line of  eq 3.60)
!       (S^{eta}_{P}S_{z,P}+S_{z,P}S^{eta}_{P}) U_{P}^{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j+1))/2.
        a2=(svetay(i,j)+svetay(i,j+1))/2.
        b3=swrz(i,j)
        tta=w(k,ii,jj)
        c3=(a1*b3)*tta
        c6=(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (2nd elament in 4th line of  eq 3.60)
!       -(S^{eta}_{S}S_{z,S}+S_{z,S}S^{eta}_{S}) U_{S}^{z}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j-1))/2.
        a2=(svetay(i,j)+svetay(i,j-1))/2.
        b3=swrz(i,j-1)
        tta=w(k,ii,jj-1)
        c3=c3-(a1*b3)*tta
        c6=c6-(a2*b3)*tta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    calculate -S^{eta}_{s}/(ReV_{s}) \bu d_{z,s,im}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j)
        a2=svetay(i,j)
        a3=(vp(i,j)+vp(i,j-1))/2.
        e3=e3-(a1*c3+a2*c6)*visc/a3 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute the total implicit diffusion terms    (eq 3.60)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        vis=swz(i,j)*e3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Calculate EXPLICIT pressure term (-dt*dp/dz)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        pres= -(swz(i,j)*swz(i,j)*p(k  ,ii,jj) &
               -swz(i,j)*swz(i,j)*p(k-1,ii,jj))/vp(i,j)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate right-hand-side vector for w equation
! (z-component of RHS of equations 3.68, 3.80, 3.115)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        rhs(k,ii,jj)=dt/vp(i,j)*(r1*con-r2*conz(k,ii,jj)+vis)+pres
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        conz(k,ii,jj)=con

 63   continue

      return
end


