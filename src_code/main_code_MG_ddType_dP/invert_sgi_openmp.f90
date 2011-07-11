

      subroutine invertu
      use data_export

      dimension a(Nmax_1),b(Nmax_1),c(Nmax_1),d(Nmax_1)
      dimension    aa(i1_u:i2_u, j1_u:j2_u),  &
                   bb(i1_u:i2_u, j1_u:j2_u),  &
                   cc(i1_u:i2_u, j1_u:j2_u),  &
            rr(ngz-1, i1_u:i2_u, j1_u:j2_u)

      !T  dimension ust(0:ngz, 0:nlx+1, 0:nly)
      !T  equivalence(ust(0,0,0),uc(0,0,0))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (1-A2)(1-A3)(ustar-u(n))        pages 24,25,
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!$OMP parallel &
!$OMP default (shared) &
!$OMP private (a1,a2,a3,b1,b2,b3,d1,d2,d3) &
!$OMP private (c1,c2,c3,c4,c5,c6,c7,c8,c9) &
!$OMP private (aa,bb,cc) &
!$OMP private (e1,e2,e3,beta,tta) &
!$OMP private(a,b,c,d)

!$OMP do
!T    do 1000 kk=1,ngz-1
      do 75 jj=j1_u,j2_u
      do 75 ii=i1_u,i2_u
        i = ibmap_1(ii)
        j = jbmap_1(jj)
!T      k = kk

        beta=-dt/(2.*(vp(i,j)+vp(i-1,j))/2.)

!---------------------------------------
!-------------( eq 3.88 )---------------
!---------------------------------------
        a1=suxix(i+1,j)
        a2=suxiy(i+1,j)
        b1=surxix(i+1,j)
        b2=surxiy(i+1,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
        d1=(suxix(i,j)+suxix(i+1,j))/2.
        d2=(suxiy(i,j)+suxiy(i+1,j))/2.
        e1=(d1*c1+d2*c4)*visc/vp(i,j) 
        e2=(d1*c2+d2*c5)*visc/vp(i,j) 

        cc(ii,jj)=beta*(suxix(i,j)*e1+suxiy(i,j)*e2)

!---------------------------------------
!-------------( eq 3.89 )---------------
!---------------------------------------
        a1=suxix(i,j)
        a2=suxiy(i,j)
        b1=surxix(i,j)
        b2=surxiy(i,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!cccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{e}/(ReV_{e}) \bu d_{xi,e,im}       
!cccccccccccccccccccccccccccccccccccccccccccc
        d1=(suxix(i,j)+suxix(i+1,j))/2.
        d2=(suxiy(i,j)+suxiy(i+1,j))/2.
        e1=(d1*c1+d2*c4)*visc/vp(i,j) 
        e2=(d1*c2+d2*c5)*visc/vp(i,j)
!cccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{w}/(ReV_{w}) \bu d_{xi,w,im}  
!cccccccccccccccccccccccccccccccccccccccccccc
        d1=(suxix(i,j)+suxix(i-1,j))/2.
        d2=(suxiy(i,j)+suxiy(i-1,j))/2.
        e1=e1+(d1*c1+d2*c4)*visc/vp(i-1,j)
        e2=e2+(d1*c2+d2*c5)*visc/vp(i-1,j)

        bb(ii,jj)=1.- beta * ( suxix(i,j)*e1 + suxiy(i,j)*e2 )

!---------------------------------------
!-------------( eq 3.90 )---------------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{W}S_{xi,W}+S_{xi,W}S^{xi}_{W})  
!ccccccccccccccccccccccccccccccccccccccccccccccc
        a1=suxix(i-1,j)
        a2=suxiy(i-1,j)
        b1=surxix(i-1,j)
        b2=surxiy(i-1,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!cccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{w}/(ReV_{w}) \bu d_{xi,w,im}  
!cccccccccccccccccccccccccccccccccccccccccccc
        d1=(suxix(i,j)+suxix(i-1,j))/2.
        d2=(suxiy(i,j)+suxiy(i-1,j))/2.
        e1=(d1*c1+d2*c4)*visc/vp(i-1,j) 
        e2=(d1*c2+d2*c5)*visc/vp(i-1,j)

        aa(ii,jj)=beta*(suxix(i,j)*e1+suxiy(i,j)*e2)

 75   continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(idudx_1.eq.0) then
        do 1001 jj=j1_u,j2_u
          bb(i1_u,jj)=bb(i1_u,jj)+aa(i1_u,jj)
 1001   continue
      end if

      if(idudx_2.eq.0) then
        do 1002 jj=j1_u,j2_u
          bb(i2_u,jj)=bb(i2_u,jj)+cc(i2_u,jj)
 1002   continue
      end if

      if (iblock.eq.1) then
        do 79 jj=j1_u,j2_u
          aa(i1_u,jj)   = 0.
 79     continue
      end if
      
      if (iblock.eq.npx) then
        do 97 jj=j1_u,j2_u
          cc(i2_u,jj)=0.
 97     continue
      end if

      rr(:,:,:) = rhs(1:ngz-1, i1_u:i2_u, j1_u:j2_u)

!*****************************************************************************
!do jj=j1_u,j2_u
!    call triDiagonalM_mrhs(id_x, aa(i1_u,jj), bb(i1_u,jj),        &
!                 cc(i1_u,jj), rr(1,i1_u,jj), niu, 1, ngz-1)
! end do
!*****************************************************************************
!                              ORIGINAL
!*****************************************************************************
  do jj=j1_u,j2_u
    call triDiagonal_mrhs(id_x, aa(i1_u,jj), bb(i1_u,jj),                &
                          cc(i1_u,jj), rr(1,i1_u,jj), ngz-1, niu)
  end do 
!*****************************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! next solve for (1-A3)(ust-u(n)) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 83 jj=j1_u,j2_u
      do 83 ii=i1_u,i2_u
        i = ibmap_1(ii)
        j = jbmap_1(jj)

        beta=-dt/(2.*(vp(i,j)+vp(i-1,j))/2.)

!---------------------------------------
!-------------( eq 3.101 )--------------
!---------------------------------------
        a1=(svetax(i,j+1)+svetax(i-1,j+1)+svetax(i,j+2)+svetax(i-1,j+2))/4.
        a2=(svetay(i,j+1)+svetay(i-1,j+1)+svetay(i,j+2)+svetay(i-1,j+2))/4.
        b1=surxix(i,j+1)
        b2=surxiy(i,j+1)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{n}/(ReV_{n}) \bu d_{xi,n,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1=(svetax(i,j+1)+svetax(i-1,j+1))/2.
        d2=(svetay(i,j+1)+svetay(i-1,j+1))/2.
        tta=(vp(i,j)+vp(i,j+1)+vp(i-1,j)+vp(i-1,j+1))/4.
        e1=(d1*c1+d2*c4)*visc/tta
        e2=(d1*c2+d2*c5)*visc/tta

        cc(ii,jj)=beta*(suxix(i,j)*e1+suxiy(i,j)*e2) 

!---------------------------------------
!-------------( eq 3.102 )--------------
!---------------------------------------
        a1=(svetax(i,j)+svetax(i,j+1)+svetax(i-1,j)+svetax(i-1,j+1))/4.
        a2=(svetay(i,j)+svetay(i,j+1)+svetay(i-1,j)+svetay(i-1,j+1))/4.
        b1=surxix(i,j)
        b2=surxiy(i,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{n}/(ReV_{n}) \bu d_{xi,n,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1=(svetax(i,j+1)+svetax(i-1,j+1))/2.
        d2=(svetay(i,j+1)+svetay(i-1,j+1))/2.
        tta=(vp(i,j)+vp(i,j+1)+vp(i-1,j)+vp(i-1,j+1))/4.
        e1=(d1*c1+d2*c4)*visc/tta
        e2=(d1*c2+d2*c5)*visc/tta
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{xi,P}+S_{xi,P}S^{eta}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j+1)+svetax(i-1,j)+svetax(i-1,j+1))/4.
        a2=(svetay(i,j)+svetay(i,j+1)+svetay(i-1,j)+svetay(i-1,j+1))/4.
        b1=surxix(i,j)
        b2=surxiy(i,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{s}/(ReV_{s}) \bu d_{xi,s,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1=(svetax(i,j)+svetax(i-1,j))/2.
        d2=(svetay(i,j)+svetay(i-1,j))/2.
        tta=(vp(i,j)+vp(i,j-1)+vp(i-1,j)+vp(i-1,j-1))/4.
        e1=e1+(d1*c1+d2*c4)*visc/tta
        e2=e2+(d1*c2+d2*c5)*visc/tta

        bb(ii,jj)=1.-beta*(suxix(i,j)*e1+suxiy(i,j)*e2)

!---------------------------------------
!-------------( eq 3.103 )--------------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{S}S_{xi,S}+S_{xi,S}S^{eta}_{S})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j-1)+svetax(i-1,j)+svetax(i-1,j-1))/4.
        a2=(svetay(i,j)+svetay(i,j-1)+svetay(i-1,j)+svetay(i-1,j-1))/4.
        b1=surxix(i,j-1)
        b2=surxiy(i,j-1)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{s}/(ReV_{s}) \bu d_{xi,s,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1=(svetax(i,j)+svetax(i-1,j))/2.
        d2=(svetay(i,j)+svetay(i-1,j))/2.
        tta=(vp(i,j)+vp(i,j-1)+vp(i-1,j)+vp(i-1,j-1))/4.
        e1=(d1*c1+d2*c4)*visc/tta
        e2=(d1*c2+d2*c5)*visc/tta

        aa(ii,jj)=beta*(suxix(i,j)*e1+suxiy(i,j)*e2)

 83   continue


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(idudy_1.eq.0) then
        do 1003 ii=i1_u,i2_u
          bb(ii,j1_u)=bb(ii,j1_u)+aa(ii,j1_u)
 1003   continue
      end if

      if(idudy_2.eq.0) then
        do 1004 ii=i1_u,i2_u
          bb(ii,j2_u)=bb(ii,j2_u)+cc(ii,j2_u)
 1004   continue
      end if

      if (jblock.eq.1) then
         do 77 ii=i1_u,i2_u
           aa(ii,j1_u)=0.
 77      continue
      end if
      
      if (jblock.eq.npy) then
        do 91 ii=i1_u,i2_u
           cc(ii,j2_u)=0.
 91     continue
      end if

!T      do 100 ii=i1_u,i2_u
!T        do 701 jj=j1_u,i2_u
!T          a(jj)=aa(ii,jj)
!T          b(jj)=bb(ii,jj)
!T          c(jj)=cc(ii,jj) 
!T          d(jj)=rhs(ii,jj,k)
!T 701    continue
!T        call trida(ny-1,a(1),b(1),c(1),d(1))
!T        do 100 jj=j1_u,i2_u
!T          ust(ii,jj,k)=d(jj)+u(ii,jj,k) 
!T 100  continue


!*****************************************************************************
!   do i=i1_u,i2_u
!     call triDiagonal_mrhs(id_y, aa(i,j1_u), bb(i,j1_u),                &
!                           cc(i,j1_u), rr(1,i,j1_u), ngz-1, nju)
!   end do
!*****************************************************************************
!                              ORIGINAL
!*****************************************************************************
     call triDiagonalM_mrhs(id_y, aa(i1_u,j1_u), bb(i1_u,j1_u),        &
                  cc(i1_u,j1_u), rr(1,i1_u,j1_u), nju, niu, ngz-1)
!*****************************************************************************

         ust(1:ngz-1,i1_u:i2_u,j1_u:j2_u) = rr(:,:,:)  &
                                           + u(1:ngz-1,i1_u:i2_u,j1_u:j2_u)



!T 1000 continue
!$OMP end do
!$OMP end parallel
      return
      end







      subroutine invertv
      use data_export

      dimension ::       a(Nmax_1),b(Nmax_1),c(Nmax_1),d(Nmax_1)
      dimension ::       aa(i1_v:i2_v, j1_v:j2_v)
      dimension ::       bb(i1_v:i2_v, j1_v:j2_v)
      dimension ::       cc(i1_v:i2_v, j1_v:j2_v)
      dimension ::rr(ngz-1, i1_v:i2_v, j1_v:j2_v)

      !T dimension vst(0:ngz, 0:nlx, 0:nly+1 )
      !T equivalence(vst(0,0,0),vc(0,0,0))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (1-A2)(1-A3)(vstar-v(n)) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!$OMP parallel &
!$OMP default (shared) &
!$OMP private (a1,a2,a3,b1,b2,b3,d1,d2,d3) &
!$OMP private (c1,c2,c3,c4,c5,c6,c7,c8,c9) &
!$OMP private (aa,bb,cc) &
!$OMP private (e1,e2,e3,beta,tta) &
!$OMP private(a,b,c,d)

!$OMP do
!T    do 1000 kk=1,ngz-1
      do 75 jj=j1_v, j2_v
      do 75 ii=i1_v, i2_v 
        i = ibmap_1(ii)
        j = jbmap_1(jj)
!T      k = kk

        beta=-dt/(2.*(vp(i,j)+vp(i,j-1))/2.) 

!---------------------------------------
!-------------( eq 3.123 )--------------
!---------------------------------------
        a1=(suxix(i+1,j)+suxix(i+2,j)+suxix(i+1,j-1)+suxix(i+2,j-1))/4.
        a2=(suxiy(i+1,j)+suxiy(i+2,j)+suxiy(i+1,j-1)+suxiy(i+2,j-1))/4.
        b1=svretax(i+1,j)
        b2=svretay(i+1,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{e}/(ReV_{e}) \bu d_{eta,e,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1=(suxix(i+1,j)+suxix(i+1,j-1))/2.
        d2=(suxiy(i+1,j)+suxiy(i+1,j-1))/2.
        tta=(vp(i,j)+vp(i+1,j)+vp(i,j-1)+vp(i+1,j-1))/4. 
        e1=(d1*c1+d2*c4)*visc/tta
        e2=(d1*c2+d2*c5)*visc/tta 

        cc(ii,jj)=beta*(svetax(i,j)*e1+svetay(i,j)*e2)

!---------------------------------------
!-------------( eq 3.124 )--------------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{P}S_{eta,P}+S_{eta,P}S^{xi}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i+1,j)+suxix(i,j-1)+suxix(i+1,j-1))/4.
        a2=(suxiy(i,j)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i+1,j-1))/4.
        b1=svretax(i,j)
        b2=svretay(i,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{e}/(ReV_{e}) \bu d_{eta,e,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1=(suxix(i+1,j)+suxix(i+1,j-1))/2.
        d2=(suxiy(i+1,j)+suxiy(i+1,j-1))/2.
        tta=(vp(i,j)+vp(i+1,j)+vp(i,j-1)+vp(i+1,j-1))/4. 
        e1=(d1*c1+d2*c4)*visc/tta 
        e2=(d1*c2+d2*c5)*visc/tta 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      (S^{xi}_{P}S_{eta,P}+S_{eta,P}S^{xi}_{P})U^{eta}_{P}
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i+1,j)+suxix(i,j-1)+suxix(i+1,j-1))/4.
        a2=(suxiy(i,j)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i+1,j-1))/4.
        b1=svretax(i,j)
        b2=svretay(i,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{w}/(ReV_{w}) \bu d_{eta,w,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1=(suxix(i,j)+suxix(i,j-1))/2.
        d2=(suxiy(i,j)+suxiy(i,j-1))/2.
        tta=(vp(i,j)+vp(i-1,j)+vp(i,j-1)+vp(i-1,j-1))/4.
        e1=e1+(d1*c1+d2*c4)*visc/tta 
        e2=e2+(d1*c2+d2*c5)*visc/tta 

        bb(ii,jj)=1.-beta*(svetax(i,j)*e1+svetay(i,j)*e2)

!---------------------------------------
!-------------( eq 3.125 )--------------
!---------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccc
!        (S^{xi}_{W}S_{eta,W}+S_{eta,W}S^{xi}_{W})
!cccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i,j-1)+suxix(i-1,j)+suxix(i-1,j-1))/4.
        a2=(suxiy(i,j)+suxiy(i,j-1)+suxiy(i-1,j)+suxiy(i-1,j-1))/4.
        b1=svretax(i-1,j)
        b2=svretay(i-1,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{w}/(ReV_{w}) \bu d_{eta,w,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1=(suxix(i,j)+suxix(i,j-1))/2.
        d2=(suxiy(i,j)+suxiy(i,j-1))/2.
        tta=(vp(i,j)+vp(i-1,j)+vp(i,j-1)+vp(i-1,j-1))/4.
        e1=(d1*c1+d2*c4)*visc/tta 
        e2=(d1*c2+d2*c5)*visc/tta 

        aa(ii,jj)=beta*(svetax(i,j)*e1+svetay(i,j)*e2)

 75   continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(idvdx_1.eq.0) then
        do 1001 jj=j1_v,j2_v
          bb(i1_v,jj)=bb(i1_v,jj)+aa(i1_v,jj)
 1001   continue
      end if

      if(idvdx_2.eq.0) then
        do 1002 jj=j1_v,j2_v
          bb(i2_v,jj)=bb(i2_v,jj)+cc(i2_v,jj)
 1002   continue
      end if

      if (iblock.eq.1) then
        do 81 jj=j1_v,j2_v
          aa(i1_v ,jj)=0.
 81     continue
      end if

      if (iblock.eq.npx) then
        do 98 jj=j1_v,j2_v
          cc(i2_v ,jj)=0.
 98     continue
      end if


      rr(:,:,:) = rhs(1:ngz-1, i1_v:i2_v, j1_v:j2_v)

!*****************************************************************************
! do jj=j1_v,j2_v
!    call triDiagonalM_mrhs(id_x, aa(i1_v,jj), bb(i1_v,jj),        &
!                 cc(i1_v,jj), rr(1,i1_v,jj), niv, 1, ngz-1)
! end do
!*****************************************************************************
!                              ORIGINAL
!*****************************************************************************
   do jj=j1_v,j2_v
     call triDiagonal_mrhs(id_x, aa(i1_v,jj), bb(i1_v,jj),                &
                           cc(i1_v,jj), rr(1,i1_v,jj), ngz-1, niv)
   end do
!*****************************************************************************


 995  format(1x,5(1pe10.3)) 


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! next solve for (1-A3)(vst-v(n)) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do 93 jj=j1_v,j2_v
      do 93 ii=i1_v,i2_v
         i = ibmap_1(ii)
         j = jbmap_1(jj)

        beta=-dt/(2.*(vp(i,j)+vp(i,j-1))/2.)

!---------------------------------------
!-------------( eq 3.147 )--------------
!---------------------------------------
        a1=svetax(i,j+1)
        a2=svetay(i,j+1)
        b1=svretax(i,j+1)
        b2=svretay(i,j+1)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!cccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{n}/(ReV_{n}) \bu d_{eta,n,im}
!cccccccccccccccccccccccccccccccccccccccccccccc
        d1=(svetax(i,j)+svetax(i,j+1))/2.
        d2=(svetay(i,j)+svetay(i,j+1))/2.
        tta=vp(i,j)
        e1=(d1*c1+d2*c4)*visc/tta 
        e2=(d1*c2+d2*c5)*visc/tta 

        cc(ii,jj)=beta*(svetax(i,j)*e1+svetay(i,j)*e2)

!---------------------------------------
!-------------( eq 3.148 )--------------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{eta,P}+S_{eta,P}S^{eta}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j)
        a2=svetay(i,j)
        b1=svretax(i,j)
        b2=svretay(i,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!cccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{n}/(ReV_{n}) \bu d_{eta,n,im}
!cccccccccccccccccccccccccccccccccccccccccccccc
        d1=(svetax(i,j)+svetax(i,j+1))/2.
        d2=(svetay(i,j)+svetay(i,j+1))/2.
        tta=vp(i,j)
        e1=(d1*c1+d2*c4)*visc/tta
        e2=(d1*c2+d2*c5)*visc/tta 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{eta,P}+S_{eta,P}S^{eta}_{P})U^{eta}_{P}
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j)
        a2=svetay(i,j)
        b1=svretax(i,j)
        b2=svretay(i,j)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!cccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{s}/(ReV_{s}) \bu d_{eta,s,im}
!cccccccccccccccccccccccccccccccccccccccccccccc
        d1=(svetax(i,j)+svetax(i,j-1))/2.
        d2=(svetay(i,j)+svetay(i,j-1))/2.
        tta=vp(i,j-1)
        e1=e1+(d1*c1+d2*c4)*visc/tta
        e2=e2+(d1*c2+d2*c5)*visc/tta 

        bb(ii,jj)=1.-beta*(svetax(i,j)*e1+svetay(i,j)*e2)

!---------------------------------------
!-------------( eq 3.149 )--------------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{S}S_{eta,S}+S_{eta,S}S^{eta}_{S})
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
        a1=svetax(i,j-1) 
        a2=svetay(i,j-1)
        b1=svretax(i,j-1)
        b2=svretay(i,j-1)
        c1=(a1*b1+b1*a1)
        c2=(a1*b2+b1*a2)
        c4=(a2*b1+b2*a1)
        c5=(a2*b2+b2*a2)
!cccccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{s}/(ReV_{s}) \bu d_{eta,s,im}
!cccccccccccccccccccccccccccccccccccccccccccccc
        d1=(svetax(i,j)+svetax(i,j-1))/2.
        d2=(svetay(i,j)+svetay(i,j-1))/2.
        tta=vp(i,j-1)
        e1=(d1*c1+d2*c4)*visc/tta 
        e2=(d1*c2+d2*c5)*visc/tta 

        aa(ii,jj)=beta*(svetax(i,j)*e1+svetay(i,j)*e2)

 93   continue


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(idvdy_1.eq.0) then
        do 1003 ii=i1_v,i2_v
          bb(ii,j1_v)=bb(ii,j1_v)+aa(ii,j1_v)
 1003   continue
      end if

      if(idvdy_2.eq.0) then
        do 1004 ii=i1_v,i2_v
          bb(ii,j2_v)=bb(ii,j2_v)+cc(ii,j2_v)
 1004   continue
      end if

      if (jblock.eq.1) then
        do 89 ii=i1_v,i2_v
          aa(ii,j1_v)=0.
 89     continue
      end if

      if (jblock.eq.npy) then
        do 95 ii=i1_v,i2_v
          cc(ii,j2_v)=0.
 95     continue
      end if

!T      do 101 i=1,nx-1
!T        do 102 j=2,ny-1
!T          a(j)=aa(i,j)
!T          b(j)=bb(i,j)
!T          c(j)=cc(i,j)
!T          d(j)=rhs(i,j,k)
!T 102    continue
!T        call trida(ny-2,a(2),b(2),c(2),d(2))
!T        do 101 j=2,ny-1
!T          vst(i,j,k)=d(j)+v(i,j,k) 
!T 101  continue

!*****************************************************************************
!   do i=i1_v,i2_v
!     call triDiagonal_mrhs(id_y, aa(i,j1_v), bb(i,j1_v),                &
!                           cc(i,j1_v), rr(1,i,j1_v), ngz-1, njv)
!   end do
!*****************************************************************************
!                              ORIGINAL
!*****************************************************************************
     call triDiagonalM_mrhs(id_y, aa(i1_v,j1_v), bb(i1_v,j1_v),        &
                      cc(i1_v,j1_v), rr(1,i1_v,j1_v), njv, niv, ngz-1)
!*****************************************************************************

         vst(1:ngz-1,i1_v:i2_v,j1_v:j2_v) = rr(:,:,:)  &
                                            + v(1:ngz-1,i1_v:i2_v,j1_v:j2_v)


!T 1000  continue
!$OMP end do
!$OMP end parallel
      return
      end










      subroutine invertw
      use data_export

      dimension a(Nmax_1),b(Nmax_1),c(Nmax_1),d(Nmax_1)

      dimension     aa(i1_w:i2_w,j1_w:j2_w),  &
                    bb(i1_w:i2_w,j1_w:j2_w),  &
                    cc(i1_w:i2_w,j1_w:j2_w),  &
            rr(1:ngz-1,i1_w:i2_w,j1_w:j2_w)

      !T dimension wst(0:ngz+1, 0:nlx, 0:nly)
      !T equivalence(wst(0,0,0),wc(0,0,0))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (1-A2)(1-A3)(wstar-w(n)) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!$OMP parallel &
!$OMP default (shared) & 
!$OMP private (a1,a2,a3,b1,b2,b3,d1,d2,d3) &
!$OMP private (c1,c2,c3,c4,c5,c6,c7,c8,c9) &
!$OMP private (aa,bb,cc) &
!$OMP private (e1,e2,e3,beta,tta) &
!$OMP private(a,b,c,d)

!$OMP do
!T    do 1000 kk=1,ngz-1
      do 75 jj=j1_w,j2_w
      do 75 ii=i1_w,i2_w
        i = ibmap_1(ii)
        j = jbmap_1(jj)
!T      k = kk

        beta=-dt/(2.*vp(i,j))

!---------------------------------------
!-------------( eq 3.158 )--------------
!---------------------------------------
        a1=(suxix(i+1,j)+suxix(i+2,j))/2.
        a2=(suxiy(i+1,j)+suxiy(i+2,j))/2.
        b3=swrz(i+1,j)
        c3=(a1*b3)
        c6=(a2*b3)
!ccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{e}/(ReV_{e}) \bu d_{z,e,im}
!ccccccccccccccccccccccccccccccccccccccccccc
        d1=suxix(i+1,j)
        d2=suxiy(i+1,j)
        tta=(vp(i,j)+vp(i+1,j))/2. 
        e3=(d1*c3+d2*c6)*visc/tta 

        cc(ii,jj)=beta*(swz(i,j)*e3)

!---------------------------------------
!-------------( eq 3.159 )--------------
!---------------------------------------
        a1=(suxix(i,j)+suxix(i+1,j))/2.
        a2=(suxiy(i,j)+suxiy(i+1,j))/2.
        b3=swrz(i,j)
        c3=(a1*b3)
        c6=(a2*b3)
!ccccccccccccccccccccccccccccccccccccccccccc
!       S^{xi}_{e}/(ReV_{e}) \bu d_{z,e,im}
!ccccccccccccccccccccccccccccccccccccccccccc
        d1=suxix(i+1,j)
        d2=suxiy(i+1,j)
        tta=(vp(i,j)+vp(i+1,j))/2. 
        e3=(d1*c3+d2*c6)*visc/tta 
        a1=(suxix(i,j)+suxix(i+1,j))/2.
        a2=(suxiy(i,j)+suxiy(i+1,j))/2.
        b3=swrz(i,j)
        c3=(a1*b3)
        c6=(a2*b3)
!cccccccccccccccccccccccccccccccccccccccccccc
!       -S^{xi}_{w}/(ReV_{w}) \bu d_{z,w,im}
!cccccccccccccccccccccccccccccccccccccccccccc
        d1=suxix(i,j)
        d2=suxiy(i,j)
        tta=(vp(i,j)+vp(i-1,j))/2.
        e3=e3+(d1*c3+d2*c6)*visc/tta 

        bb(ii,jj)=1.-beta*(swz(i,j)*e3)

!---------------------------------------
!-------------( eq 3.160 )--------------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{xi}_{W}S_{z,W}+S_{z,W}S^{xi}_{W})
!ccccccccccccccccccccccccccccccccccccccccccccc
        a1=(suxix(i,j)+suxix(i-1,j))/2.
        a2=(suxiy(i,j)+suxiy(i-1,j))/2.
        b3=swrz(i-1,j)
        c3=(a1*b3)
        c6=(a2*b3)
!cccccccccccccccccccccccccccccccccccccccccccc
!       -S^{xi}_{w}/(ReV_{w}) \bu d_{z,w,im}
!cccccccccccccccccccccccccccccccccccccccccccc
        d1=suxix(i,j)
        d2=suxiy(i,j)
        tta=(vp(i,j)+vp(i-1,j))/2.
        e3=(d1*c3+d2*c6)*visc/tta 

        aa(ii,jj)=beta*(swz(i,j)*e3)

 75   continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(idwdx_1.eq.0) then
        do 1001 jj=j1_w,j2_w
          bb(i1_w,jj)=bb(i1_w,jj)+aa(i1_w,jj)
 1001   continue
      end if

      if(idwdx_2.eq.0) then
        do 1002 jj=j1_w,j2_w
          bb(i2_w,jj)=bb(i2_w,jj)+cc(i2_w,jj)
 1002   continue
      end if

      if (iblock.eq.1) then
        do 92 jj=j1_w,j2_w
          aa(i1_w, jj)=0.
 92     continue
      end if

      if (iblock.eq.npx) then
        do 94 jj=j1_w,j2_w
          cc(i2_w,jj)=0.
 94     continue
      end if


      rr(:,:,:) = rhs(1:ngz-1, i1_w:i2_w, j1_w:j2_w)

!*****************************************************************************
! do jj=j1_w,j2_w
!    call triDiagonalM_mrhs(id_x, aa(i1_w,jj), bb(i1_w,jj),        &
!                 cc(i1_w,jj), rr(1,i1_w,jj), niw, 1, ngz-1)
! end do
!*****************************************************************************
!                              ORIGINAL
!*****************************************************************************
     do jj=j1_w,j2_w
       call triDiagonal_mrhs(id_x, aa(i1_w,jj), bb(i1_w,jj),                &
                             cc(i1_w,jj), rr(1,i1_w,jj), ngz-1, niw)
     end do
!*****************************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! next solve for (1-A3)(wst-w(n)) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do 96 jj=j1_w,j2_w
      do 96 ii=i1_w,i2_w
        i = ibmap_1(ii)
        j = jbmap_1(jj)

        beta=-dt/(2.*vp(i,j))

!---------------------------------------
!-------------( eq 3.171 )--------------
!---------------------------------------
        a1=(svetax(i,j+1)+svetax(i,j+2))/2.
        a2=(svetay(i,j+1)+svetay(i,j+2))/2.
        b3=swrz(i,j+1)
        c3=(a1*b3)
        c6=(a2*b3)
!cccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{n}/(ReV_{n}) \bu d_{z,n,im}
!cccccccccccccccccccccccccccccccccccccccccccc
        d1=svetax(i,j+1)
        d2=svetay(i,j+1)
        tta=(vp(i,j)+vp(i,j+1))/2. 
        e3=(d1*c3+d2*c6)*visc/tta 

        cc(ii,jj)=beta*(swz(i,j)*e3) 

!---------------------------------------
!-------------( eq 3.172 )--------------
!---------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{z,P}+S_{z,P}S^{eta}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j+1))/2.
        a2=(svetay(i,j)+svetay(i,j+1))/2.
        b3=swrz(i,j)
        c3=(a1*b3)
        c6=(a2*b3)
!cccccccccccccccccccccccccccccccccccccccccccc
!       S^{eta}_{n}/(ReV_{n}) \bu d_{z,n,im}
!cccccccccccccccccccccccccccccccccccccccccccc
        d1=svetax(i,j+1)
        d2=svetay(i,j+1)
        tta=(vp(i,j)+vp(i,j+1))/2. 
        e3=(d1*c3+d2*c6)*visc/tta 
!ccccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{P}S_{z,P}+S_{z,P}S^{eta}_{P})
!ccccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j+1))/2.
        a2=(svetay(i,j)+svetay(i,j+1))/2.
        b3=swrz(i,j)
        c3=(a1*b3)
        c6=(a2*b3)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       -S^{eta}_{s}/(ReV_{s}) \bu d_{z,s,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1=svetax(i,j)
        d2=svetay(i,j)
        tta=(vp(i,j)+vp(i,j-1))/2.
        e3=e3+(d1*c3+d2*c6)*visc/tta 

        bb(ii,jj)=1.-beta*(swz(i,j)*e3)

!---------------------------------------
!-------------( eq 3.173 )--------------
!---------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccc
!       (S^{eta}_{S}S_{z,S}+S_{z,S}S^{xi}_{S})
!cccccccccccccccccccccccccccccccccccccccccccccc
        a1=(svetax(i,j)+svetax(i,j-1))/2.
        a2=(svetay(i,j)+svetay(i,j-1))/2.
        b3=swrz(i,j-1)
        c3=(a1*b3)
        c6=(a2*b3)
!ccccccccccccccccccccccccccccccccccccccccccccc
!       -S^{eta}_{s}/(ReV_{s}) \bu d_{z,s,im}
!ccccccccccccccccccccccccccccccccccccccccccccc
        d1=svetax(i,j)
        d2=svetay(i,j)
        tta=(vp(i,j)+vp(i,j-1))/2.
        e3=(d1*c3+d2*c6)*visc/tta 

        aa(ii,jj)=beta*swz(i,j)*e3

 96   continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Neumann b.c.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(idwdy_1.eq.0) then
        do 1003 ii=i1_w,i2_w
          bb(ii,j1_w)=bb(ii,j1_w)+aa(ii,j1_w)
 1003   continue
      end if

      if(idwdy_2.eq.0) then
        do 1004 ii=i1_w,i2_w
          bb(ii,j2_w)=bb(ii,j2_w)+cc(ii,j2_w)
 1004   continue
      end if

      if (jblock.eq.1) then
        do 112 ii=i1_w,i2_w
          aa(ii,j1_w)=0.
 112    continue
      end if

      if (jblock.eq.npy) then 
        do 121 ii=i1_w,i2_w
          cc(ii,j2_w)=0.
 121    continue
      end if

!T      do 115 i=1,nx-1
!T        do 116 j=1,ny-1
!T          a(j)=aa(i,j)
!T          b(j)=bb(i,j)
!T          c(j)=cc(i,j)
!T          d(j)=rhs(i,j,k)
!T 116    continue
!T        call trida(ny-1,a(1),b(1),c(1),d(1))
!T        do 115 j=1,ny-1
!T          wst(i,j,k)=d(j)+w(i,j,k)
!T 115  continue


!*****************************************************************************
!   do i=i1_w,i2_w
!     call triDiagonal_mrhs(id_y, aa(i,j1_w), bb(i,j1_w),                &
!                           cc(i,j1_w), rr(1,i,j1_w), ngz-1, njw)
!   end do
!*****************************************************************************
!                              ORIGINAL
!*****************************************************************************
         call triDiagonalM_mrhs(id_y, aa(i1_w,j1_w), bb(i1_w,j1_w),        &
                       cc(i1_w,j1_w), rr(1,i1_w,j1_w), njw, niw, ngz-1)
!*****************************************************************************


         wst(1:ngz-1,i1_w:i2_w,j1_w:j2_w) = rr(:,:,:)  &
                                           + w(1:ngz-1,i1_w:i2_w,j1_w:j2_w)


!T 1000  continue
!$OMP end do
!$OMP end parallel
      return
      end
