
module data_export

!==============================================================
! Get user-specified parameters
!
      include "param.inc"

!==============================================================
! Global limits

      parameter (nox=1, noy=1, noz=1)
      parameter (nix =ngx , niy =ngy , niz =ngz )
      parameter (nixm=ngxm, niym=ngym, nizm=ngzm)
      parameter (nx=nix+2*nox, ny=niy+2*noy, nz=niz+2*noz  , nxyz=nx*ny*nz, Nmax=max(nx,ny,nz))
      parameter (nxm=nx-1    , nym=ny-1    , nzm=nz-1      )

      parameter (imin=0+nox, imax=imin+nix-1, imino=imin-nox, imaxo=imax+nox)
      parameter (jmin=0+noy, jmax=jmin+niy-1, jmino=jmin-noy, jmaxo=jmax+noy)
      parameter (kmin=0+noz, kmax=kmin+niz-1, kmino=kmin-noz, kmaxo=kmax+noz)


! Primary layout (QUICK Scheme)
! NOTE: here, i USED (NOX,NOY,NOZ) instead of (NOX1,NOY1,NOZ1) to have enough memory!!!!!
! VERY IMPORTANT!!!!
      parameter (nbx_1=npx, nby_1=npy, nbz_1=npz)
      parameter (nox1=3, noy1=3, noz1=1)
      parameter (  nix_1  = (nx-2*nox-1)/nbx_1 + 1 + min(1,mod(nx-2*nox-1,nbx_1))  )
      parameter (  niy_1  = (ny-2*noy-1)/nby_1 + 1 + min(1,mod(ny-2*noy-1,nby_1))  )
      parameter (  niz_1  = (nz-2*noz-1)/nbz_1 + 1 + min(1,mod(nz-2*noz-1,nbz_1))  )
      parameter (  nixyz_1= nix_1*niy_1*niz_1  )

! Transposed layout
      parameter (nbx_2=1, nby_2=1, nbz_2=nproc)
      parameter (nox2=1, noy2=1, noz2=0)
      parameter (  nix_2  = nix  )
      parameter (  niy_2  = niy  )
      parameter (  niz_2  = (niz-1)/nproc + 1 + min(1,mod(nizm,nproc))  )
      parameter (  nixyz_2= nix_2*niy_2*niz_2  )

! Poisson layout
      parameter (nixp = ngxm/npx, niyp=ngym/npy)

! Array dimensions (static)
!NOTE:  SHOULD USE (NOX1,NOY1) HERE.  Need to fix this to make it general!!!!
! FIX
      parameter (nx_1 = nix_1 + 2*nox1, nxm_1=nx_1-1)
      parameter (ny_1 = niy_1 + 2*noy1, nym_1=ny_1-1)
      parameter (nz_1 = niz_1 + 2*noz1, nzm_1=nz_1-1)
      parameter (nxyz_1=nx_1*ny_1*nz_1, Nmax_1=max(nx_1,ny_1,nz_1))
      parameter (i1_=0+nox1, i2_=i1_+nix_1-1)
      parameter (j1_=0+noy1, j2_=j1_+niy_1-1)
      parameter (k1_=0+noz1, k2_=k1_+niz_1-1)
      parameter (nx_=nx_1, ny_=ny_1, nz_=nz_1, nxyz_=nxyz_1)

      parameter (nx_2 = nix_2 + 2*nox2, nxm_2=nx_2-1)
      parameter (ny_2 = niy_2 + 2*noy2, nym_2=ny_2-1)
      parameter (nz_2 = niz_2 + 2*noz2, nzm_2=nz_2-1)
      parameter (nxyz_2=nx_2*ny_2*nz_2, Nmax_2=max(nx_1,ny_1,nz_1))

! Local Number of grid points for Xiaohua's code
      parameter (nlx =nx_1-2*nox, nly =ny_1-2*noy)
      parameter (nlxm=nlx-1     , nlym=nly-1     )
      parameter (nlz =nz_2 )
      parameter (nlzm=nlz-1)


! Block limits (dynamic)
      integer nxb_1 (nbx_1), nyb_1 (nby_1), nzb_1 (nbz_1)
      integer nixb_1(nbx_1), niyb_1(nby_1), nizb_1(nbz_1)
      integer ibmin_1(nbx_1), ibmax_1(nbx_1), ibmino_1(nbx_1), ibmaxo_1(nbx_1), &
              jbmin_1(nby_1), jbmax_1(nby_1), jbmino_1(nby_1), jbmaxo_1(nby_1), &
              kbmin_1(nbz_1), kbmax_1(nbz_1), kbmino_1(nbz_1), kbmaxo_1(nbz_1)

      integer nxb_2 (nbx_2), nyb_2 (nby_2), nzb_2 (nbz_2)
      integer nixb_2(nbx_2), niyb_2(nby_2), nizb_2(nbz_2)
      integer ibmin_2(nbx_2), ibmax_2(nbx_2), ibmino_2(nbx_2), ibmaxo_2(nbx_2), &
              jbmin_2(nby_2), jbmax_2(nby_2), jbmino_2(nby_2), jbmaxo_2(nby_2), &
              kbmin_2(nbz_2), kbmax_2(nbz_2), kbmino_2(nbz_2), kbmaxo_2(nbz_2)

      integer nlxb_(npx), nlyb_(npy), nlzb_(nproc)
      integer nlxb, nlyb, nlzb

! Shorthand for primary layout
      integer nxb, nyb, nzb
      integer nixb, niyb, nizb
      integer ibmin, ibmax, ibmino, ibmaxo
      integer jbmin, jbmax, jbmino, jbmaxo
      integer kbmin, kbmax, kbmino, kbmaxo
      integer ib1_, ib2_, jb1_, jb2_, kb1_, kb2_

! Indeces for advancing velocity
      integer i1_u, i2_u, niu
      integer i1_v, i2_v, niv
      integer i1_w, i2_w, niw

      integer j1_u, j2_u, nju
      integer j1_v, j2_v, njv
      integer j1_w, j2_w, njw

! Indeces for FFT and Transpose
      integer i1_T, j1_T, k1_T
      integer i2_T, j2_T, k2_T

      integer iTmin_1(nbx_1), iTmax_1(nbx_1),  &
              jTmin_1(nby_1), jTmax_1(nby_1),  &
              kTmin_2(nbz_2), kTmax_2(nbz_2)

! Indeces for Poisson Solver
      integer i1_P, i2_P
      integer j1_P, j2_P

! Index mapping
      integer imap_1 (0:nx-1)  , jmap_1 (0:ny-1)  , kmap_1 (0:nz-1)
      !TAZ--------------------------------------------------------------
      !TAZ  DANGER ZONE: (added some padding just because
      !TAZ               I need it in Jacobs' [boundaries.blayer.f90])
      !TAZ--------------------------------------------------------------
      !TAZ  integer ibmap_1(0:nx_1-1), jbmap_1(0:ny_1-1), kbmap_1(0:nz_1-1)
      !TAZ--------------------------------------------------------------
      integer ibmap_1(0:nx_1)  , jbmap_1(0:ny_1)  , kbmap_1(0:nz_1)
      integer imap_2 (0:nx-1)  , jmap_2 (0:ny-1)  , kmap_2 (0:nz-1)
      integer ibmap_2(0:nx_2-1), jbmap_2(0:ny_2-1), kbmap_2(0:nz_2-1)

! Block/Process ID
      integer irank, iroot, ierr
      integer iblock, jblock, kblock
      integer iblock_1, jblock_1, kblock_1
      integer iblock_2, jblock_2, kblock_2

! BC for each face of the domain
      integer BC_bottom, BC_top, BC_front, BC_back

!=======================================================================
! Grid-Metrics
	real     suxix  (0:ngx+1,0:ngy  ),suxiy  (0:ngx+1,0:ngy  )   &
		,svetax (0:ngx  ,0:ngy+1),svetay (0:ngx  ,0:ngy+1)   &
		,spz    (0:ngx  ,0:ngy  ),swz    (0:ngx  ,0:ngy  )   &
		,surxix (0:ngx+1,0:ngy  ),surxiy (0:ngx+1,0:ngy  )   &
		,svretax(0:ngx  ,0:ngy+1),svretay(0:ngx  ,0:ngy+1)   &
		,swrz   (0:ngx  ,0:ngy  ),vp     (0:ngx  ,0:ngy  )
	real    zpp(ngz-1),xpp(ngx-1,ngy-1),ypp(ngx-1,ngy-1),   &
		zpu(ngz-1),xpu(ngx  ,ngy-1),ypu(ngx  ,ngy-1),   &
		zpv(ngz-1),xpv(ngx-1,ngy  ),ypv(ngx-1,ngy  ),   &
		zpw(ngz  ),xpw(ngx-1,ngy-1),ypw(ngx-1,ngy-1),   &
		zpg(ngz  ),xpg(ngx  ,ngy  ),ypg(ngx  ,ngy  )

!=======================================================================
! Main data arrays

      ! real    p,u,v,w,uc,vc,wc

      real :: p(0:ngz,0:nlx  ,0:nly),dp (0:ngz,0:nlx,0:nly  ),   &
            u  (0:ngz,0:nlx+1,0:nly),v  (0:ngz,0:nlx,0:nly+1),w  (0:ngz+1,0:nlx,0:nly),   &
            uc (0:ngz,0:nlx+1,0:nly),vc (0:ngz,0:nlx,0:nly+1),wc (0:ngz+1,0:nlx,0:nly),   &
            ust(0:ngz,0:nlx+1,0:nly),vst(0:ngz,0:nlx,0:nly+1),wst(0:ngz+1,0:nlx,0:nly)

      equivalence(ust(0,0,0),uc(0,0,0))
      equivalence(vst(0,0,0),vc(0,0,0))
      equivalence(wst(0,0,0),wc(0,0,0))

      real          ucbcs(0:ngz  ,0:nlx+1,2),ucbcn(0:ngz  ,0:nlx+1,2),   &
                    ucbcw(0:ngz  ,0:nly  ,3),ucbce(0:ngz  ,0:nly  ,3),   &
                    ucbcb(0:nlx+1,0:nly  ,1),ucbcf(0:nlx+1,0:nly  ,1),   &
                    vcbcs(0:ngz  ,0:nlx  ,2),vcbcn(0:ngz  ,0:nlx  ,2),   &
                    vcbcw(0:ngz  ,0:nly+1,3),vcbce(0:ngz  ,0:nly+1,3),   &
                    vcbcb(0:nlx  ,0:nly+1,1),vcbcf(0:nlx  ,0:nly+1,1),   &
                    wcbcs(0:ngz+1,0:nlx  ,2),wcbcn(0:ngz+1,0:nlx  ,2),   &
                    wcbcw(0:ngz+1,0:nly  ,2),wcbce(0:ngz+1,0:nly  ,2),   &
                    wcbcb(0:nlx  ,0:nly  ,1),wcbcf(0:nlx  ,0:nly  ,1)

!------------------------------------------
! Some settings, and simulation parameters
!------------------------------------------
      integer       ipout,iprob,imodel, &
                    iavet,itime,iread,ipoisson,   &
                    idudx_1,idudx_2,idudy_1,idudy_2,idudz_1,idudz_2,   &
                    idvdx_1,idvdx_2,idvdy_1,idvdy_2,idvdz_1,idvdz_2,   &
                    idwdx_1,idwdx_2,idwdy_1,idwdy_2,idwdz_1,idwdz_2,   &
                    ncounter_avt,iread_avt,ncounter_avp(nphase)

      real    before,after,t_init,t_over
      real          alenx,aleny,alz,dt,visc, ufree,   &
                    angle_attack,counter_avt,cpu_wavet(ngzm)

      real          stime_, stime		! Start time, and current time of simulation
      integer       ntime_, ntime		! Starting step number, and current step

    !------------------------------------------
    ! Define id numbers for variables
      parameter (id_u=1, id_v=2, id_w=3, id_p=4, id_vt=5)

    ! Define id numbers for directions
      parameter (id_x=1, id_y=2, id_z=3)

    !---------------------------------------
    ! Need some vectors to hold the boundary
    ! since it is not all periodic in y

    !FIX      real  flxU(0:ngz,0:nlx+1,3), flxV(0:ngz,0:nlx,3), flxW(0:ngz+1,0:nlx,3)
    !FIX      real  flxP(0:ngz,0:nlx+1,1)

!-----------------------------------------------------------
! 'conx' ==> 'QUICK' correction terms (5-points Stencil)
!-----------------------------------------------------------
      ! real    conx,cony,conz,rhs,phatr
      real ::  conx(ngz-1, nlx-1, nly-1), &
               cony(ngz-1, nlx-1, nly-1), &
               conz(ngz-1, nlx-1, nly-1), &
               rhs(ngz+1  , 0:nlx ,0:nly  ), &
               qr (ngz+1  , 0:nlx, 0:nly  ), &
               phatr(ngzm , 0:nixp+1, 0:niyp+1 )
      equivalence(rhs(1,0,0),qr(1,0,0))

      real , allocatable :: qT(:,:,:), q2T(:,:,:)

	!==============================================================
	!     for multigrid Poisson solver
	!==============================================================

	!----- maximum grid levels -------------
	parameter (nglevel_max=5)
	!---------------------------------------

	!---------------------------------------------------------------
	!     I think 'vkz' should be 'vkz(2*(ngzm/2))'
	!     Note:  ngzm=niz-1
	!---------------------------------------------------------------
	real   vkz(niz-1),omegak(niz-1)

	real            resmin, resminb
	integer         nitmax, nitmaxb
	integer         ncycle,npre,npost,ngrid,i_allocate

	!--------------------------------------------------------------
	!SY     Derived Data Type for Local Domains (2D multigrid)
	!--------------------------------------------------------------
	type mglevel_local
	   integer nx, ny
	   real, dimension(:,:), allocatable :: zapw , zape , zaps , zapn ,         &
	                                        zapp , zapws, zapwn, zapes, zapen,  &
	       					zvp  , zsux , zsuy , zsvx , zsvy ,  &
	       					zp   , zrhs , zres
	   real, dimension(2) :: zedge
	end type

        type (mglevel_local), dimension(nglevel_max) :: mgl

!===============================================================
!     for grid turbulence
!--------------------------------------------------------------
!TAZ  take to (boundaries.f90) or some other routine.
!--------------------------------------------------------------
      parameter (nxalan=256/4 ,nyalan=640/4 ,nzalan=128/16+1)
      ! integer nxalan,nyalan,nzalan
    !===============================================================
    !     for grid turbulence
    !---------------------------------------------------------------
      real           ualan0(nxalan,nyalan,nzalan)   &
                    ,valan0(nxalan,nyalan,nzalan)   &
                    ,walan0(nxalan,nyalan,nzalan)   &
                    ,xalan(nxalan),yalan(nyalan),zalan(nzalan)


! Notation
! Jacob's Notation
!
! ngx        ! number of grid points in x-direction (from grid generator)
!            ! (equal to number of physical grid points)
! nix        ! number of interior points in x (Global Domain)
!            !
! nx         ! number of points in x (Global Domain including halos)
! nox        ! number of overlap in x
! nbx        ! number of blocks in x
! nxb        ! number of points in x in block
! nixb       ! number of interior points in x in block
! nxyz       ! total number of points
! nixyz      ! total number of interior points
! imin       ! minimum i index
! ibmin      ! minimum block i index
! imax       ! maximum i index
! ibmax      ! maximum block i index
! i          ! i point index
! ib         ! i block index
! imino      ! imin with overlap
! imaxo      ! imax with overlap

! Notation
! Tamer's Notation
!
! ngx        ! number of grid points in x-direction (from grid data)
! nix        ! number of interior grid points
! nx         ! Total number of points in x including all halo cells
! nox        ! number of overlap points in x (on each side)
! nbx        ! number of blocks in x
! ibmin      ! minimum block i index (Physical cell index)
! ibmax      ! maximum block i index (Physical cell index)
! i          ! i point index
! ib         ! i block index
! imino      ! imin with overlap
! imaxo      ! imax with overlap

! nlx        ! General number of local grid points (eq. to nx in Xiaohua)
! nlxb       ! Number of local grid points in block
! i1_u       ! First (i) index in advance(u)
! i1_T       ! First (i) index in FFT and Transpose
! iTmin_1    ! projection of (i1_T) to global coordinates

end module
