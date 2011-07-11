
! Exported mesh definitions

pointer x(:), xm(:)
pointer y(:), ym(:)
pointer z(:), zm(:)
pointer xpg(:,:), ypg(:,:)

common /mesh_/ xL, yL, zL, &
               x, xm, dx, &
               y, ym, dy, &
               z, zm, dz, &
               xpg, ypg

common /limits_/ nx , ny , nz , nxyz , &
                 nxm, nym, nzm,        &
                 nix, niy, niz, nixyz, &
                 imin, imax, iminm1, imaxp1, &
                 jmin, jmax, jminm1, jmaxp1, &
                 kmin, kmax, kminm1, kmaxp1

