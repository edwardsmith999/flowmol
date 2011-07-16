!=======================================================================
! Reading inflow data from a database file
! and interpolating in time and space to the mesh
!
! inflow_init()
! inflow_read()
! inflow_write()
! inflow_ok(deltaT, iok)
! inflow_getData(deltaT, Uin)
! inflow_cache(n2, y, ym, n3, z, zm)
!

module inflow
        integer iu_inflow              ! unit number for inflow data file
        real    fileStart              ! time at beginning of inflow file
        real    fileEnd                ! time at end of inflow file
        real    refTime                ! reference time for inflow generator
        real    time                   ! time since start of inflow
        real    time1, time2           ! times for inflow data
        real    timeMax                ! time available in inflow data file

        allocatable data1(:,:,:)       ! inflow data
        allocatable data2(:,:,:)
        allocatable data(:,:,:)

        integer IFny, IFnz             ! inflow mesh
        real, allocatable :: IFy(:), IFym(:)
        real    IFdz, IFdx

        allocatable jy(:), jym(:)      ! interpolation indices
        allocatable kz(:), kzm(:)
        allocatable wy1(:), wy2(:)     ! interpolation weights
        allocatable wym1(:), wym2(:)
        allocatable wz1(:), wz2(:)
        allocatable wzm1(:), wzm2(:)

        integer ny, nz                 ! size of destination mesh

end module

!=======================================================================
subroutine inflow_init()
        use inflow, data0 => data
        use data, only : file_dir

        ! inflow.dat file format:
        !   ny, nz
        !   y(ny), ym(ny)
        !   dz, dx
        !   time1
        !   data1
        !   time2
        !   data2
        !   ...

        ! Determine end time in inflow data file
        iu_inflow = iopen()
        open (iu_inflow, file=trim(file_dir)//"inflow.dat", &
              form="unformatted", status="old", position="append")
        backspace (iu_inflow)
        backspace (iu_inflow)
        read (iu_inflow) timeMax
        rewind (iu_inflow)
        fileEnd = timeMax

        ! Read first two data
        ! For the case that inflow_read() is not called
        read (iu_inflow) IFny, IFnz
        allocate(data1(IFny,IFnz,3), data2(IFny,IFnz,3))
        allocate(data0(IFny,IFnz,3))
        read (iu_inflow) d
        read (iu_inflow) d
        read (iu_inflow) time1
        read (iu_inflow) data1
        read (iu_inflow) time2
        read (iu_inflow) data2
        refTime = time1
        fileStart = time1
        time = time1
        rewind (iu_inflow)

        ! Read inflow generator mesh
        read (iu_inflow) IFny, IFnz
        allocate(IFy(IFny), IFym(IFny))
        read (iu_inflow) IFy, IFym
        read (iu_inflow) IFdz !#, IFdx

        return
end

subroutine inflow_read()
        use inflow

        call readFloat("IFstart", oldStart)
        call readFloat("IFend", oldEnd)
        if ((fileStart == oldStart .and. fileEnd == oldEnd) &
                .or. (fileStart > oldEnd) ) then
                ! Same or next inflow data file
                call readFloat("IFrefTime", refTime)
                call readFloat("IFtime", time)
                call readFloat("IFtime1", time1)
                call readFloat("IFtime2", time2)
                call readArray("IFdata1", data1, IFny*IFnz*3)
                call readArray("IFdata2", data2, IFny*IFnz*3)
                if (fileStart - oldEnd > 1.) then
                        ! Time gap larger than 1 time unit
                        stop "Inflow: database files are not contiguous"
                end if
        else
                ! Restart inflow database at beginning
        end if

        return
end

subroutine inflow_write()
        use inflow

        call writeFloat("IFstart", fileStart)
        call writeFloat("IFend", fileEnd)
        call writeFloat("IFrefTime", refTime)
        call writeFloat("IFtime", time)
        call writeFloat("IFtime1", time1)
        call writeFloat("IFtime2", time2)
        call writeArray("IFdata1", data1, IFny*IFnz*3)
        call writeArray("IFdata2", data2, IFny*IFnz*3)

        return
end

subroutine inflow_ok(deltaT, iok)
        use inflow

        ! Determine whether sufficient inflow data remains
        if (time + deltaT > timeMax) then
                iok = 0
        else
                iok = 1
        end if

        return
end

!=======================================================================
subroutine inflow_getData(deltaT, Uin)
        use inflow
        real Uin(ny,nz,3)

        time = time + deltaT

        do while (time > time2)
                time1 = time2
                data1 = data2
                read (iu_inflow, iostat=ierr) time2
                read (iu_inflow, iostat=ierr) data2
                if (ierr .ne. 0) stop "Out of inflow data"
        end do

        if (time < time1 .or. time > time2) &
                stop "Something wrong with inflow data"

        ! Interpolate in time between data1 and data2
        data = +(time2-time)/(time2-time1)*data1 &
               +(time-time1)/(time2-time1)*data2

        ! Interpolate in space
        do j=1,ny
        do k=1,nz
                j1 = jym(j)
                j2 = j1 + 1
                k1 = kzm(k)
                k2 = k1 + 1
                Uin(j,k,1) = wym1(j)*( wzm1(k)*data(j1,k1,1) &
                                      +wzm2(k)*data(j1,k2,1) ) &
                            +wym2(j)*( wzm1(k)*data(j2,k1,1) &
                                      +wzm2(k)*data(j2,k2,1) )
        end do
        end do

        do j=1,ny
        do k=1,nz
                j1 = jy(j)
                j2 = j1 + 1
                k1 = kzm(k)
                k2 = k1 + 1
                Uin(j,k,2) = wy1(j)*( wzm1(k)*data(j1,k1,2) &
                                     +wzm2(k)*data(j1,k2,2) ) &
                            +wy2(j)*( wzm1(k)*data(j2,k1,2) &
                                     +wzm2(k)*data(j2,k2,2) )
        end do
        end do

        do j=1,ny
        do k=1,nz
                j1 = jym(j)
                j2 = j1 + 1
                k1 = kz(k)
                k2 = k1 + 1
                Uin(j,k,3) = wym1(j)*( wz1(k)*data(j1,k1,3) &
                                      +wz2(k)*data(j1,k2,3) ) &
                            +wym2(j)*( wz1(k)*data(j2,k1,3) &
                                      +wz2(k)*data(j2,k2,3) )
        end do
        end do

        return
end

subroutine inflow_cache(n2, y, ym, n3, z, zm)
        use inflow
        real y(n2), ym(n2), z(n3), zm(n3)
        real IFz(IFnz), IFzm(IFnz)

        ! Pre-compute interpolation weights

        ny = n2
        nz = n3

        allocate(jy(ny), jym(ny))
        allocate(kz(nz), kzm(nz))

        allocate(wy1(ny), wy2(ny))
        allocate(wym1(ny), wym2(ny))
        allocate(wz1(nz), wz2(nz))
        allocate(wzm1(nz), wzm2(nz))

        ! Subtle issue involving grid alignment and interpolation:
        !    All interior points (k=2...nz-1) must lie between alignment 
        !    points (z=0, z=L).  Otherwise, end points at k=2 or k=nz-1 
        !    could end up outside the valid interpolation region.
        !    This is particularly subtle for staggered grids.

        ! Generate the uniform IFz-mesh
        do k=1,IFnz
                IFz(k) = (k-1)*IFdz
                IFzm(k) = (k-1-0.5)*IFdz
        end do

        do j=1,ny
                jy(j) = lowerIndex(IFy, y(j), IFny-1)
                j1 = jy(j)
                j2 = j1 + 1
                wy1(j) = (IFy(j2)-y(j))/(IFy(j2)-IFy(j1))
                wy2(j) = (y(j)-IFy(j1))/(IFy(j2)-IFy(j1))
        end do

        do j=1,ny
                jym(j) = lowerIndex(IFym, ym(j), IFny)
                j1 = jym(j)
                j2 = j1 + 1
                wym1(j) = (IFym(j2)-ym(j))/(IFym(j2)-IFym(j1))
                wym2(j) = (ym(j)-IFym(j1))/(IFym(j2)-IFym(j1))
        end do

        do k=1,nz
                kz(k) = lowerIndex(IFz, z(k), IFnz)
                k1 = kz(k)
                k2 = k1 + 1
                wz1(k) = (IFz(k2)-z(k))/(IFz(k2)-IFz(k1))
                wz2(k) = (z(k)-IFz(k1))/(IFz(k2)-IFz(k1))

                kzm(k) = lowerIndex(IFzm, zm(k), IFnz)
                k1 = kzm(k)
                k2 = k1 + 1
                wzm1(k) = (IFzm(k2)-zm(k))/(IFzm(k2)-IFzm(k1))
                wzm2(k) = (zm(k)-IFzm(k1))/(IFzm(k2)-IFzm(k1))
        end do

        ! Do not allow extrapolation
        where (wy1 < -0.1 .or. wy1 > 1.1) wy1 = 0.
        where (wy2 < -0.1 .or. wy2 > 1.1) wy2 = 0.
        where (wym1 < -0.1 .or. wym1 > 1.1) wym1 = 0.
        where (wym2 < -0.1 .or. wym2 > 1.1) wym2 = 0.
        where (wz1 < -0.1 .or. wz1 > 1.1) wz1 = 0.
        where (wz2 < -0.1 .or. wz2 > 1.1) wz2 = 0.
        where (wzm1 < -0.1 .or. wzm1 > 1.1) wzm1 = 0.
        where (wzm2 < -0.1 .or. wzm2 > 1.1) wzm2 = 0.

        return
end

