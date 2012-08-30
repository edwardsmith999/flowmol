
module mesh_export

! Exported mesh definitions

real xL, yL, zL

allocatable x(:), xm(:)
real        dx, dxm, dxi, dxmi

allocatable y(:), ym(:), yi(:), ymi(:)
allocatable dy(:), dym(:), dyi(:), dymi(:)

allocatable z(:), zm(:)
real        dz, dzm, dzi, dzmi

! Weights and metrics for finite-volume calculations

allocatable cym1(:), cym2(:)

real        axu, axv, axw, axp
allocatable ayu(:), ayv(:), ayw(:), ayp(:)
real        azu, azv, azw, azp

end module

