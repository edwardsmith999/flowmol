module mesh_parameters
! Exported mesh definitions

pointer x(:), xm(:)
pointer y(:), ym(:)
pointer z(:), zm(:)
pointer xpg(:,:), ypg(:,:)

real            xL, yL, zL, &
               x, xm, dx, &
               y, ym, dy, &
               z, zm, dz, &
               xpg, ypg

integer          nx , ny , nz , nxyz , &
                 nxm, nym, nzm,        &
                 nix, niy, niz, nixyz, &
                 imin, imax, iminm1, iminp1, imaxm1, imaxp1, &
                 jmin, jmax, jminm1, jminp1, jmaxm1, jmaxp1, &
                 kmin, kmax, kminm1, kminp1, kmaxm1, kmaxp1

end module mesh_parameters
