!==========================================================
! User-specified compile parameters
!
! Statically define problem size and number of processors
!==========================================================

! NOTE:  + (ngx) is the global number of grid points (cell faces) without halos
!        + Make sure (ngx-1) is a multiple of number of processors
!          in order to make sure the following works:
!             ==> QUICK SCHEME
!              -> Tridiagonalm.mpi.f90
!              -> Parallel I/O
!              -> A couple of more things
!        + In order to generalize ngx/ngxm, make sure to check entire code.
!        + Changing the tri-diagonal solver is nearly impossible? (may be easy)
!

!--Define Grid Dimension
   parameter (ngx = 8+1 , ngy = 8+1, ngz = 8+1, ngxyz =ngx *ngy *ngz  , nphase=26)
   parameter (ngxm= ngx-1 , ngym= ngy-1, ngzm= ngz-1, ngxyzm=ngxm*ngym*ngzm            )

   parameter (npx=1, npy=1, npz=1, nproc=npx*npy*npz)


