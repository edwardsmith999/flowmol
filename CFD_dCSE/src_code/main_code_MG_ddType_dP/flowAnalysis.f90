!=======================================================================
! Perform flow analysis somewhat diagnostic in nature.
! These can either be after each timestep or at the
! end of the run.
!
! flowAnalysis_init()
! flowAnalysis_read()
! flowAnalysis_write()
! flowAnalysis_blayerinfo()
!

module flowAnalysis
	use data_export
	use mesh_export
end module

!=======================================================================
subroutine flowAnalysis_init()
	use flowAnalysis
	return
end

subroutine flowAnalysis_read()
	use flowAnalysis
	return
end

subroutine flowAnalysis_write()
	use flowAnalysis
	return
end

subroutine flowAnalysis_blayerinfo()
! Compute boundary layer integral quantities at each i location

	use flowAnalysis
	parameter (nflowInfo = 1)
	real Re				! Suffix zero implies at the inlet of the domain
	real flowInfo(nflowInfo,ngx)
	integer ii, jj, Jinfty

	! Define Falkner-Skan parameters
	call readFloat("Re",Re)
	call readFloat("ubulk", ubulk)

	!-------------------------------------
	i1 = iTmin_1(iblock) ; if (iblock ==  1 ) i1 =  1  ; i1 = imap_1(i1)
	i2 = iTmax_1(iblock) ; if (iblock == npx) i2 = ngx ; i2 = imap_1(i2)
	j1 = jTmin_1(jblock) ; if (jblock ==  1 ) j1 =  1  ; j1 = jmap_1(j1)
	j2 = jTmax_1(jblock) ; if (jblock == npy) j2 = ngy ; j2 = jmap_1(j2)

	flowInfo = 0.

	do i=i1,i2
		ii = ibmap_1(i)

		!--- Cf = 2 nu dudy/Uinf^2
		if (jblock == 1) then
			flowInfo(1,ii) = 2./ (Re * ubulk**2) * uc(ngzm/2,i,1) / ( dym(1)/2. )
			! print*, 'a',flowInfo(1,ii) 
			flowInfo(1,ii) = 2./ (Re * ubulk**2) * uc(ngzm/2,i,1) / ( 0.5*(ypg(ii,2)-ypg(ii,1)) )
			! print*, 'b',flowInfo(1,ii) 
		end if
	enddo

	call globalSum(flowInfo, nflowInfo*ngx)

	! Write file header
	if (irank == iroot) then
		iunit = iopen()
		open (iunit, file="flowAnalysis", form="formatted")
		write(iunit,'(a)') '% ------------------------------------------------------ '
		write(iunit,'(a)') '%       Re        ubulk '
		write(iunit,'(a,f15.8,f15.8/)') '%',Re, ubulk
		write(iunit,'(a)') '% ------------------------------------------------------ '
		write(iunit,'(1a,2(5x,a) )') &
			'%','Re_x', 'Cf'

		! Write solutions to file
		do i=1,ngx	! Uc is defined from 1:ngx
			xtemp = x(i)
			write(iunit,'(es10.4,1(f11.6))') Re*ubulk*xtemp, flowInfo(1,i)
		enddo
		close (iclose(iunit))
	endif

	return
end

