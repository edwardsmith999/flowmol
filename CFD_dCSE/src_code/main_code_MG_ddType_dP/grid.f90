!=======================================================================
! Mesh definitions for cartesian coordinates on staggered grid
!
! grid_init()
! anuint()
! grid_halo()
! grid_write()
!
!

module grid
        use data_export
end module

!=======================================================================
! Grid and Metric Definition
!

subroutine grid_init()
        use grid
        use data, only : file_dir
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! a. Read in "xpg", "ypg"
        ! 1. generate curvilinear grids X(xi,eta),Y(xi,eta),Z(z)
        ! 2. generate area vectors and volumes to define the mapping
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real, allocatable, dimension(:,:) ::  spxix ,spxiy
        real, allocatable, dimension(:,:) ::  spetax,spetay

        !cccccccccccccccccccc
        ! Read in the Grid
        !cccccccccccccccccccc
        !T  iunit = iopen()
        !T  open(iunit,file='grid.data',form='unformatted', &
        !T           status="old", iostat=ierr)
        !T  if (ierr .ne. 0) stop "grid.data file is required"
        !T  read(iunit) ni
        !T  read(iunit) xpg
        !T  read(iunit) nj
        !T  read(iunit) ypg
        !T  close(iclose(iunit))
        !T  if (ni .ne. nix .or. nj .ne. niy) &
        !T          stop "Dimensions in [grid.data] and [data.export.f90] files do not match"

        inquire(iolength=ilength) xpg
        iunit = iopen()
        open(iunit, file=trim(file_dir)//"grid.data", form="unformatted", access="direct", &
                recl=ilength, status="old", iostat=ierr)
        if (ierr .ne. 0) stop "grid.data file is required"
        read(iunit, rec=1) xpg
        read(iunit, rec=2) ypg
        close(iclose(iunit))

        allocate(spxix (ngx,ngy-1))
        allocate(spxiy (ngx,ngy-1))
        allocate(spetax(ngx-1,ngy))
        allocate(spetay(ngx-1,ngy))

        call readInt("iprob", iprob)
        twopi=2.*2.*acos(0.)
        pi=twopi/2.
        call readFloat("zL",alz)

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! calculate the location of p points in physical domain
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do k=1,ngz
                zeta=float(k-1)
                delta=2.0
                x1=delta/2.*(zeta/float(ngz-1)-1.)
                x2=delta/2.
                tanh1=(exp(x1)-exp(-x1))/(exp(x1)+exp(-x1))
                tanh2=(exp(x2)-exp(-x2))/(exp(x2)+exp(-x2))
        !        coef3=1.+tanh1/tanh2
                coef3=zeta/float(ngz-1)
                zpw(k)=(1.0-coef3)*0.0+coef3*alz
        end do

        do k=1,ngz
                zpg(k)=zpw(k)
        end do

        do k=1,ngz-1 
                zpp(k)=0.5*(zpw(k)+zpw(k+1))
                zpu(k)=0.5*(zpw(k)+zpw(k+1))
                zpv(k)=0.5*(zpw(k)+zpw(k+1))
        end do

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! calculate the location of u points in physical domain
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do j=1,ngy-1
        do i=1,ngx 
                xpu(i,j)=0.5*(xpg(i,j)+xpg(i,j+1))
                ypu(i,j)=0.5*(ypg(i,j)+ypg(i,j+1))
        end do
        end do

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! calculate the location of v points in physical domain
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do j=1,ngy
        do i=1,ngx-1
                xpv(i,j)=0.5*(xpg(i,j)+xpg(i+1,j))
                ypv(i,j)=0.5*(ypg(i,j)+ypg(i+1,j))
        end do
        end do

        do j=1,ngy-1
        do i=1,ngx-1
                xpp(i,j)=0.5*(xpu(i,j)+xpu(i+1,j))
                ypp(i,j)=0.5*(ypv(i,j)+ypv(i,j+1))
        end do
        end do

        do j=1,ngy-1
        do i=1,ngx-1
                xpw(i,j)=xpp(i,j)
                ypw(i,j)=ypp(i,j)
        end do
        end do

        pzk=alz/float(ngz-1)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the area vector for the p control volume
! faces (e,w,n,s,f,b): spxi(e,w),speta(n,s), spz(f,b)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !cccccccccccccccccccccccccc
        ! for point e,w
        !cccccccccccccccccccccccccc
        do j=1,ngy-1
        do i=1,ngx
                petax=xpg(i,j+1)-xpg(i,j)
                petay=ypg(i,j+1)-ypg(i,j)
                spxix(i,j)= petay*pzk 
                spxiy(i,j)=-petax*pzk
        end do
        end do
        !cccccccccccccccccccccccccc
        ! for point n,s  
        !cccccccccccccccccccccccccc
        do j=1,ngy
        do i=1,ngx-1
                pxix=xpg(i+1,j)-xpg(i,j)
                pxiy=ypg(i+1,j)-ypg(i,j)
                spetax(i,j)=-pzk*pxiy  
                spetay(i,j)= pzk*pxix
        end do
        end do
        !cccccccccccccccccccccccccc
        ! for point f,b
        !cccccccccccccccccccccccccc
        do j=1,ngy-1   
        do i=1,ngx-1
                pxix =(xpg(i+1,j)+xpg(i+1,j+1))/2.-(xpg(i,j)+xpg(i,j+1))/2.
                pxiy =(ypg(i+1,j)+ypg(i+1,j+1))/2.-(ypg(i,j)+ypg(i,j+1))/2.
                petax=(xpg(i,j+1)+xpg(i+1,j+1))/2.-(xpg(i,j)+xpg(i+1,j))/2.
                petay=(ypg(i,j+1)+ypg(i+1,j+1))/2.-(ypg(i,j)+ypg(i+1,j))/2.
                spz(i,j)=pxix*petay-pxiy*petax
        end do
        end do

        do j=1,ngy-1
                spz(0,j)  =spz(1    ,j)
                spz(ngx,j)=spz(ngx-1,j)
        end do
        do i=0,ngx
                spz(i,0)  =spz(i,1)
                spz(i,ngy)=spz(i,ngy-1)
        end do

        !ccccccccccccccccccccccccccccccccccccccccccccccccc
        ! calcuate the volume size for p control volume
        !ccccccccccccccccccccccccccccccccccccccccccccccccc
        do j=1,ngy-1
        do i=1,ngx-1
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! first find the positions of fne and bsw for p control volume
                ! calculate rfne-rbsw
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                rdifx=xpg(i+1,j+1)-xpg(i,j)
                rdify=ypg(i+1,j+1)-ypg(i,j)
                rdifz=pzk   
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! calculate spxiw+spseta+spbz
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                sumx=spxix(i,j)+spetax(i,j)
                sumy=spxiy(i,j)+spetay(i,j)
                sumz=spz(i,j)
                sumx=sumx/3.
                sumy=sumy/3.
                sumz=sumz/3.
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                ! perform inner product to obtain Vp
                !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                vp(i,j)=sumx*rdifx+sumy*rdify+sumz*rdifz
        end do
        end do

        do i=1,ngx-1
                vp(i,0  )=vp(i,1)
                vp(i,ngy)=vp(i,ngy-1)
        end do
        do j=0,ngy
                vp(0  ,j)=vp(1    ,j)
                vp(ngx,j)=vp(ngx-1,j)
        end do

        !ccccccccccccccccccccccccccccccc
        ! check that a cell is closed   
        !ccccccccccccccccccccccccccccccc
        spxmax=-1.e3
        spymax=-1.e3
        spzmax=-1.e3
        sptx=0.  
        spty=0.
        sptz=0.
        do j=1,ngy-1
        do i=1,ngx-1
                spsumx=spxix(i+1,j)-spxix(i,j)+spetax(i,j+1)-spetax(i,j)   
                spsumy=spxiy(i+1,j)-spxiy(i,j)+spetay(i,j+1)-spetay(i,j)   
                spsumz=spz(i,j)-spz(i,j)
                sptx=sptx+spsumx
                spty=spty+spsumy
                sptz=sptz+spsumz
                if(abs(spsumx).gt.spxmax) spxmax=abs(spsumx)
                if(abs(spsumy).gt.spymax) spymax=abs(spsumy)
                if(abs(spsumz).gt.spzmax) spzmax=abs(spsumz)
        end do
        end do
        if (irank.eq.1) then 
                write(6,999) spxmax,spymax,spzmax,sptx,spty,sptz
        end if

        !ccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! check that cell sizes sum to total domain volume
        !ccccccccccccccccccccccccccccccccccccccccccccccccccc
        vsum=0.
        do k=1,ngz-1
        do j=1,ngy-1
        do i=1,ngx-1
                vsum=vsum+vp(i,j)
        end do
        end do
        end do

        if (irank.eq.1) write(6,999) vsum
 999        format(1x,6(1pe10.3))

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute u location area vectors interior nodes
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,ngx
        do j=1,ngy-1
                suxix(i,j)=spxix(i,j)
                suxiy(i,j)=spxiy(i,j)
        end do
        end do
        !ccccccccccccccccccccccccccccccc
        ! lower and upper b.c.
        !ccccccccccccccccccccccccccccccc
        do i=1,ngx
                suxix(i,0)=suxix(i,1)
                suxiy(i,0)=suxiy(i,1)
                suxix(i,ngy)=suxix(i,ngy-1)
                suxiy(i,ngy)=suxiy(i,ngy-1)
        end do
        !ccccccccccccccccccccccccccccccc
        ! left and right b.c.
        !ccccccccccccccccccccccccccccccc
        do j=0,ngy
                suxix(0,j)=suxix(2,j)
                suxiy(0,j)=suxiy(2,j)
                suxix(ngx+1,j)=suxix(ngx-1,j)
                suxiy(ngx+1,j)=suxiy(ngx-1,j)
        end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! compute v location area vectors interior nodes
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
        do j=1,ngy
        do i=1,ngx-1
                svetax(i,j)=spetax(i,j)
                svetay(i,j)=spetay(i,j)
        end do
        end do
        !ccccccccccccccccccccccccccccccc
        ! left and right b.c.
        !ccccccccccccccccccccccccccccccc
        do j=1,ngy
                svetax(0,j)=svetax(1,j)
                svetay(0,j)=svetay(1,j)
                svetax(ngx,j)=svetax(ngx-1,j)
                svetay(ngx,j)=svetay(ngx-1,j)
        end do
        !ccccccccccccccccccccccccccccccc
        ! lower and upper b.c.
        !ccccccccccccccccccccccccccccccc
        do i=0,ngx
                svetax(i,0)=svetax(i,2)
                svetay(i,0)=svetay(i,2)
                svetax(i,ngy+1)=svetax(i,ngy-1)
                svetay(i,ngy+1)=svetay(i,ngy-1)
        end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the reciprocal area vectors for u locations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do j=1,ngy-1
        do i=1,ngx
                a1=(svetax(i-1,j+1)+svetax(i,j+1)+svetax(i-1,j)+svetax(i,j))/4. 
                a2=(svetay(i-1,j+1)+svetay(i,j+1)+svetay(i-1,j)+svetay(i,j))/4.
                b3=(spz(i-1,j)+spz(i,j))/2.
                c1=a2*b3
                c2=-a1*b3
                surxix(i,j)=c1/(suxix(i,j)*c1+suxiy(i,j)*c2)
                surxiy(i,j)=c2/(suxix(i,j)*c1+suxiy(i,j)*c2)
        end do
        end do
        !ccccccccccccccccccccccccccccccc
        ! lower and upper b.c.
        !ccccccccccccccccccccccccccccccc
        do i=1,ngx
                surxix(i,0)=surxix(i,1)
                surxiy(i,0)=surxiy(i,1)
                surxix(i,ngy)=surxix(i,ngy-1)
                surxiy(i,ngy)=surxiy(i,ngy-1)
        end do
        !ccccccccccccccccccccccccccccccc
        ! left and right b.c.
        !ccccccccccccccccccccccccccccccc
        do j=0,ngy
                surxix(0,j)=surxix(2,j)
                surxiy(0,j)=surxiy(2,j)
                surxix(ngx+1,j)=surxix(ngx-1,j)
                surxiy(ngx+1,j)=surxiy(ngx-1,j)
        end do

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the reciprocal area vectors for v locations
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do j=1,ngy 
        do i=1,ngx-1 
                b1=(suxix(i+1,j-1)+suxix(i+1,j)+suxix(i,j-1)+suxix(i,j))/4.
                b2=(suxiy(i+1,j-1)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i,j))/4.
                a3=(spz(i,j-1)+spz(i,j))/2.
                c1=-a3*b2
                c2=a3*b1
                svretax(i,j)=c1/(svetax(i,j)*c1+svetay(i,j)*c2)
                svretay(i,j)=c2/(svetax(i,j)*c1+svetay(i,j)*c2)
        end do
        end do
        !cccccccccccccccccccccccccc
        ! left and right b.c.
        !cccccccccccccccccccccccccc
        do j=1,ngy
                svretax(0 ,j)=svretax(1,j)
                svretay(0 ,j)=svretay(1,j)
                svretax(ngx,j)=svretax(ngx-1,j)
                svretay(ngx,j)=svretay(ngx-1,j)
        end do
        !cccccccccccccccccccccccccc
        ! lower and upper b.c.
        !cccccccccccccccccccccccccc
        do i=0,ngx
                svretax(i,0   )=svretax(i,2)
                svretay(i,0   )=svretay(i,2)
                svretax(i,ngy+1)=svretax(i,ngy-1)
                svretay(i,ngy+1)=svretay(i,ngy-1)
        end do

!ccccccccccccccccccccccccccccccccccccccccc
! area vectors for w locations
!ccccccccccccccccccccccccccccccccccccccccc
        do j=1,ngy-1
        do i=1,ngx-1
                swz(i,j)=spz(i,j)
        end do
        end do
        !cccccccccccccccccccccccccc
        ! left and right b.c. 
        !cccccccccccccccccccccccccc
        do j=1,ngy-1
                swz(0  ,j)=swz(1,j)
                swz(ngx,j)=swz(ngx-1,j)
        end do
        !cccccccccccccccccccccccccc
        ! lower and upper b.c.
        !cccccccccccccccccccccccccc
        do i=0,ngx
                swz(i,0  )=swz(i,1)
                swz(i,ngy)=swz(i,ngy-1)
        end do

        do i=0,ngx
        do j=0,ngy
                swrz(i,j)=1./swz(i,j)
        end do
        end do

        deallocate(spxix )
        deallocate(spxiy )
        deallocate(spetax)
        deallocate(spetay)

        return
end

subroutine anuint()
        use grid
        dimension counter_avp(nphase)
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! tell if Neumann B.cs are needed
        ! NOTE : Setting (idudx = 0) ==>  Newmann BC in effect (du/dx=0)
        !        Otherwise (idudx=1) ==>  No Newmann BC.
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! Initialize all Newmann BC to 1
        idudx_1=1
        idudx_2=1
        idvdx_1=1
        idvdx_2=1 
        idwdx_1=1
        idwdx_2=1
        
        idudy_1=1
        idudy_2=1
        idvdy_1=1
        idvdy_2=1
        idwdy_1=1 
        idwdy_2=1 

        idudz_1=1 
        idudz_2=1 
        idvdz_1=1   
        idvdz_2=1
        idwdz_1=1
        idwdz_2=1
        ! Set Newmann BC to 0/1 for border processors
        if (iblock.eq.1) then
                idudx_1=1 
                idvdx_1=1
                idwdx_1=1  
        end if
        if (iblock.eq.npx) then
                idudx_2=1
                idvdx_2=1
                idwdx_2=1
        end if
        if (jblock.eq.1) then
                idudy_1=1 
                idvdy_1=1 
                idwdy_1=1
        end if
        if (jblock.eq.npy) then
                idudy_2=1
                idvdy_2=1
                idwdy_2=1   
        end if
        if (kblock.eq.1) then
                idudz_1=1
                idvdz_1=1
                idwdz_1=1
        end if
        if (kblock.eq.npz) then
                idudz_2=1
                idvdz_2=1
                idwdz_2=1
        end if

        !if(iread_avt.eq.1) then
        !        open(29,file='average.input',form='unformatted')
        !        read(29) uavt,vavt,wavt,pavt,urmst,vrmst,wrmst,prmst,uvt,put,pvt,pwt,counter_avt
        !        close(29)
        !        ncounter_avt=int(counter_avt)
        !        write(6,*) 'average data read', ncounter_avt
        !else
        !        ncounter_avt=0
        !        uavt=0.0
        !        vavt=0.0
        !        wavt=0.0
        !        pavt=0.0
        !        urmst=0.0
        !        vrmst=0.0
        !        wrmst=0.0
        !        prmst=0.0
        !        uvt=0.0 
        !        put=0.0 
        !        pvt=0.0 
        !        pwt=0.0   
        !        ncounter_avp=0
        !end if   

        return
end



subroutine grid_halo()
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! correct upper and lower halo values
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use grid
        real, allocatable, dimension(:,:) :: spxix , spxiy
        real, allocatable, dimension(:,:) :: spetax,spetay

        allocate(spxix (ngx,ngy-1))
        allocate(spxiy (ngx,ngy-1))
        allocate(spetax(ngx-1,ngy))
        allocate(spetay(ngx-1,ngy))
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! perform inner product to obtain Vp
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,ngx-1
                vp(i,0 )=vp(i,1   )
                vp(i,ngy)=vp(i,ngy-1)
        end do

        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! compute u location area vectors interior nodes
        ! lower and upper b.c.
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,ngx
                suxix(i,0 ) = suxix(i,1)
                suxiy(i,0 ) = suxiy(i,1)
                suxix(i,ngy) = suxix(i,ngy-1)
                suxiy(i,ngy) = suxiy(i,ngy-1)
        end do
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! compute v location area vectors interior nodes
        ! lower and upper b.c.
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=0,ngx
                svetax(i,0   ) = svetax(i,2)
                svetay(i,0   ) = svetay(i,2)
                svetax(i,ngy+1) = svetax(i,ngy-1)    
                svetay(i,ngy+1) = svetay(i,ngy-1)    
        end do
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! calculate the reciprocal area vectors for u locations
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,ngx
        do j=1,ngy-1  
                a1=(svetax(i-1,j+1)+svetax(i,j+1)+svetax(i-1,j)+svetax(i,j))/4.
                a2=(svetay(i-1,j+1)+svetay(i,j+1)+svetay(i-1,j)+svetay(i,j))/4.
                b3=(spz(i-1,j)+spz(i,j))/2. 
                c1=a2*b3
                c2=-a1*b3
                surxix(i,j)=c1/(suxix(i,j)*c1+suxiy(i,j)*c2)
                surxiy(i,j)=c2/(suxix(i,j)*c1+suxiy(i,j)*c2)
        end do
        end do
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! lower and upper b.c.
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,ngx
                surxix(i,0 ) = surxix(i,1)
                surxiy(i,0 ) = surxiy(i,1)
                surxix(i,ngy) = surxix(i,ngy-1)
                surxiy(i,ngy) = surxiy(i,ngy-1)
        end do
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! calculate the reciprocal area vectors for v locations
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,ngx-1
        do j=1,ngy
                b1=(suxix(i+1,j-1)+suxix(i+1,j)+suxix(i,j-1)+suxix(i,j))/4.
                b2=(suxiy(i+1,j-1)+suxiy(i+1,j)+suxiy(i,j-1)+suxiy(i,j))/4.
                a3=(spz(i,j-1)+spz(i,j))/2.
                c1=-a3*b2
                c2=a3*b1
                svretax(i,j)=c1/(svetax(i,j)*c1+svetay(i,j)*c2)
                svretay(i,j)=c2/(svetax(i,j)*c1+svetay(i,j)*c2)
        end do
        end do
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! lower and upper b.c.
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=0,ngx
                svretax(i,0   ) = svretax(i,2)
                svretay(i,0   ) = svretay(i,2)
                svretax(i,ngy+1) = svretax(i,ngy-1)
                svretay(i,ngy+1) = svretay(i,ngy-1)
        end do
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! area vectors for w locations
        ! lower and upper b.c.
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=0,ngx
                swz(i,0 )=swz(i,1   )
                swz(i,ngy)=swz(i,ngy-1)
        end do
        do i=0,ngx
        do j=0,ngy
                swrz(i,j)=1./swz(i,j)
        end do
        end do

        deallocate(spxix )
        deallocate(spxiy )
        deallocate(spetax)
        deallocate(spetay)
        return
end


subroutine grid_write()
        use grid

        call writeFloat("xL", xL)
        call writeFloat("yL", yL)
        call writeFloat("zL", zL)
        call writeArray("y", y, ny)
        call writeArray("ym", ym, ny)

        ! Subclass responsibility
        call subMesh_write()

        return
end     


