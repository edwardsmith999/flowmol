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
    use computation_parameters
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

	!----- Calculate Re_{tau}=SQRT(Re*(<dudy>_xz))
	!	- calculate Re_{tau} on top and bottom walls
	!	- average Re_{tau}
	!

	!-----------------------------------------------------------------------
	! TJ: This is a quick way to validate the navier-slip BC is being
	!     enforced correctly. If it is enforced correctly, then
	!     Lslip = u|wall / dudy|wall will be satisfied over the top
	!     and bottom walls
	!    
	! TJ: What about w|wall? There is no mean-flow in W at t=0, only a
	!     perturbation.
	!
	!
	!-----------------------------------------------------------------------

	!--- Initialize wall velocities and x-z planes ---
	uc_bot    = 0.
	uc_top    = 0.
	uc_xz_top = 0.
	uc_xz_bot = 0.
	du_bot    = 0.
	du_top    = 0.

	wc_bot    = 0.
	wc_top    = 0.
	wc_xz_top = 0.
	wc_xz_bot = 0.
	dw_bot    = 0.
	dw_top    = 0.
	

	!--- Initialize gradient x-z planes ---
	dudy_xz_bot = 0.
	dudy_xz_top = 0.
	dwdy_xz_bot = 0.
	dwdy_xz_top = 0.

	!--- Initialize friction Reynolds number ----
	Re_tau_bot  = 0.
	Re_tau_top  = 0.
	Re_tau_avg  = 0.

	!--- Define dS ----
	dS          = 1./DBLE(ngxm*ngzm)

	do i=i1,i2
	do k=1,ngzm
		kk = k
		ii = ibmap_1(i)
	
	!--- Loop over bottom wall ---

		if (jblock.eq.1) then
		
			uc_bot = 	  0.5*(uc(k,i,1)+uc(k,i,0))
			wc_bot =          0.5*(wc(k,i,1)+wc(k,i,0))

			du_bot = uc(k,i,1) - uc_bot
			dw_bot = wc(k,i,1) - wc_bot

			dy_bot =          0.5*(ypg(ii,2)-ypg(ii,1))

			dudy_xz_bot = dudy_xz_bot + ABS(du_bot/dy_bot)
			dwdy_xz_bot = dwdy_xz_bot + ABS(dw_bot/dy_bot)

			uc_xz_bot   = uc_xz_bot + uc_bot
			wc_xz_bot   = wc_xz_bot + wc_bot

		end if
	
	!--- Loop over top wall ---

		if (jblock.eq.npy) then
			
			uc_top =	  0.5*(uc(k,i,nlyb-1)+uc(k,i,nlyb))
			wc_top =	  0.5*(wc(k,i,nlyb-1)+wc(k,i,nlyb))

			du_top = uc(k,i,nlyb-1) - uc_top
			dw_top = wc(k,i,nlyb-1) - wc_top
		
			dy_top =          0.5*(ypg(ii,ngy-1)-ypg(ii,ngy))
			
			dudy_xz_top = dudy_xz_top + ABS(du_top/dy_top)
			dwdy_xz_top = dwdy_xz_top + ABS(dw_top/dy_top)

			uc_xz_top   = uc_xz_top + uc_top 
			wc_xz_top   = wc_xz_top + wc_top 

		end if

	end do
	end do

	!--- Sum over processes ---
	call globalSum(uc_xz_bot,   1)
	call globalSum(uc_xz_top,   1)
	call globalSum(wc_xz_bot,   1)
	call globalSum(wc_xz_top,   1)

	call globalSum(dudy_xz_bot, 1)
	call globalSum(dudy_xz_top, 1)
	call globalSum(dwdy_xz_bot, 1)
	call globalSum(dwdy_xz_top, 1)

	!--- Compute the XZ-planar averages ---
	uc_xz_bot   =   uc_xz_bot*dS
	uc_xz_top   =   uc_xz_top*dS
	wc_xz_bot   =   wc_xz_bot*dS
	wc_xz_top   =   wc_xz_top*dS
	dudy_xz_bot = dudy_xz_bot*dS
	dudy_xz_top = dudy_xz_top*dS
	dwdy_xz_bot = dwdy_xz_bot*dS
	dwdy_xz_top = dwdy_xz_top*dS

	!--- Compute friction Reynolds number ---
	Re_tau_bot = SQRT(Re*dudy_xz_bot)
	Re_tau_top = SQRT(Re*dudy_xz_top)

	!--- Compute average friction Reynolds number
	Re_tau_avg = SQRT(0.5*Re*(dudy_xz_bot+dudy_xz_top))

	!--- Print results to run_history ---
	if (irank.eq.1) write(6,*)'Re_tau_bot   =',Re_tau_bot
	if (irank.eq.1) write(6,*)'Re_tau_top   =',Re_tau_top
	if (irank.eq.1) write(6,*)'Re_tau_avg   =',Re_tau_avg

	if (irank.eq.1) write(6,*)'uc_xz_bot    =',uc_xz_bot
	if (irank.eq.1) write(6,*)'uc_xz_top    =',uc_xz_top
	if (irank.eq.1) write(6,*)'dudy_xz_bot  =',dudy_xz_bot
	if (irank.eq.1) write(6,*)'dudy_xz_top  =',dudy_xz_top

	if (irank.eq.1) write(6,*)'wc_xz_bot    =',wc_xz_bot
	if (irank.eq.1) write(6,*)'wc_xz_top    =',wc_xz_top
	if (irank.eq.1) write(6,*)'dwdy_xz_bot  =',dwdy_xz_bot
	if (irank.eq.1) write(6,*)'dwdy_xz_top  =',dwdy_xz_top
	
	if (irank.eq.1) then 
		if (abs(dudy_xz_bot) .gt. 1.d-12) then 
		write(6,*)'Lslip_u_bot  =',uc_xz_bot/dudy_xz_bot
                else
                write(6,*)'Lslip_u_bot  diverges,,uc_xz_bot,dudy_xz_bot   ',uc_xz_bot, dudy_xz_bot
                endif
   
                if (abs(dudy_xz_top) .gt. 1.d-12) then 
                write(6,*)'Lslip_u_top  =',uc_xz_bot/dudy_xz_top
                else
                write(6,*)'Lslip_u_top  diverges,,uc_xz_bot,dudy_xz_top   ',uc_xz_bot,dudy_xz_top
                endif

                if(abs(dwdy_xz_bot) .gt. 1.d-12) then 
        	write(6,*)'Lslip_w_bot  =',wc_xz_bot/dwdy_xz_bot
		else
		write(6,*)'Lslip_w_bot  diverges,,wc_xz_bot,dwdy_xz_bot   ',wc_xz_bot, dwdy_xz_bot
		endif

    		if (abs(dwdy_xz_top) .gt. 1.d-12) then
		write(6,*)'Lslip_w_top  =',wc_xz_bot/dwdy_xz_top
                else
                write(6,*)'Lslip_w_top  diverges,,wc_xz_bot,dwdy_xz_top   ',wc_xz_bot,dwdy_xz_top
                endif
	endif
	!---- Calculate skin friction, Cf
	!
	!
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
		open (iunit, file=trim(prefix_dir)//"flowAnalysis", form="formatted")
		write(iunit,'(a)') '% ------------------------------------------------------ '
		write(iunit,'(a)') '%       Re        ubulk '
		write(iunit,'(a,f15.8,f15.8/)') '%',Re, ubulk
		write(iunit,'(a)') '% ------------------------------------------------------ '
		write(iunit,'(1a,6(5x,a) )') &
			'%','Re_x', 'Cf'

		! Write solutions to file
		!do i=1,ngx	! Uc is defined from 1:ngx
		!	xtemp = x(i)
		!	write(iunit,'(es10.3,1(f11.6))') Re*ubulk*xtemp, flowInfo(1,i) 
		!
		!enddo
		
		close (iclose(iunit))
	endif

	return
end

