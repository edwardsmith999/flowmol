c-----------------------------------------------------------------------
c   DNS of a turbulent Couette flow with Re=400.
c   Reference: 
c
c
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      udiff =0.
      utrans=0.
      return

      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ffx = 0.0
      ffy = 0.0
      ffz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      qvol   = 0.0
      source = 0.0
      visc = param(2)

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      use CPL, only : CPL_create_comm, cfd_realm
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'  ! for nelx,nely,nelz
      common /nekmpi/ nekcomm

      integer :: ierr, CFD_COMM

      CFD_COMM = nekcomm

      call CPL_create_comm(cfd_realm,CFD_COMM,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      if (iside .eq. 3) then
         ux=1.0
      elseif (iside .eq. 1) then
         ux=-1.0
      endif
      uy=0.0
      uz=0.0

      temp=0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer idum
      save    idum 
      data    idum / 0 /

      if (idum.eq.0) idum = 99 + nid

c     Linear profile w/ random perturbations
      eps = .5
      ux  = y-1.  + eps*(ran1(idum)-.5)
      uy  =         eps*(ran1(idum)-.5)
      uz  =         eps*(ran1(idum)-.5)

      temp=0

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat   ! This routine to modify element vertices
      include 'SIZE'      ! _before_ mesh is generated, which 
      include 'TOTAL'     ! guarantees GLL mapping of mesh.

      common /cdsmag/ ediff(lx1,ly1,lz1,lelv)

!      one  = 1.
!      pi   = 4.*atan(one)
!      twopi= 2.*pi

!      xlen = twopi        !  Here, map domain onto pi x [-1,1] x 2pi
!      ylen = 2.0          !  Here, map domain onto pi x [-1,1] x 2pi
!      zlen = pi           !  Here, map domain onto pi x [-1,1] x 2pi

!      n=8*nelv
!      xmin = glmin(xc,n)
!      xmax = glmax(xc,n)
!      !ymin = glmin(yc,n)
!      !ymax = glmax(yc,n)
!      zmin = glmin(zc,n)
!      zmax = glmax(zc,n)

!      xscale = xlen/(xmax-xmin)
!      !yscale = ylen/(ymax-ymin)
!      zscale = zlen/(zmax-zmin)

!      do i=1,n
!         x       = xc(i,1)
!         xc(i,1) = xscale*x

!         !y       = yc(i,1)
!         !yc(i,1) = yscale*y-1.

!         z       = zc(i,1)
!         zc(i,1) = zscale*z
!      enddo

      param(59) = 0     ! all elements are _undeformed_

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2   ! This routine to modify mesh coordinates
      include 'SIZE'
      include 'TOTAL'

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'

      return
      end
c
c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)
      return
      end
