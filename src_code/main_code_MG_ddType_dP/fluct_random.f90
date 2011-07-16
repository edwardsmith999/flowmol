subroutine fluct_random(Uin)
        use fluct

        real Uin (ny, nz, 3)
        real Uinp(ny, nz, 3)
        real flowIn1,  flowIn2
        real amplitude
        real umean
        integer :: j, k
        integer :: i1, i2, j1, j2, k1, k2

        ! Add random fluctuations

        ! Initialize random generator
        iseed = 7
        call randomSeed(iseed)

        ! Fluctuation amplitude is typically 30%
        amplitude = 0.3

        do k= kmin, kmax
        do j=jmin, jmax
                umean = max(Uin(j,k,1), Uin(j,k,2), Uin(j,k,3))  
                Uinp(j,k,1) = umean*amplitude*2.*(random()-0.5)
                Uinp(j,k,2) = umean*amplitude*2.*(random()-0.5)
                Uinp(j,k,3) = umean*amplitude*2.*(random()-0.5)
        end do
        end do

        !Adjust fluctuations such that mass flux doesnt change with time
        i1=ibmin; i2=ibmax; j1=jbmin; j2=jbmax; k1=kbmin; k2=kbmax;
        call massFlowRate_X(flowIn1, Uinp(1,1,1), 1)
        call massFlowRate_X(flowIn2, abs(Uinp(1:ny,1:nz,1)), 1)
        if (flowIn2 > 0.) Uinp(:,:,1) = Uinp(:,:,1) -abs(Uinp(:,:,1))* flowIn1/flowIn2

        Uin = Uin + Uinp

        return


        return
end
