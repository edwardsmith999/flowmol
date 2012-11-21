!=======================================================================
! Highest algorithm level; controls overall simulation and runs the time loop
!
! simulation_init()
! simulation_read()
! simulation_write()
! simulation_run()
! simulation_info(iunit)
!

module simulation

        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	use data_export
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	implicit none
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	INTERFACE
		real*8 function realClock()
		end function realClock
	END INTERFACE

	logical iexist
	CHARACTER(LEN=47) ::   FN2
	CHARACTER(LEN=41) ::   NAME

	real cpu_cflcal , cpu_bc   , cpu_advx , cpu_invx    , cpu_advy
	real cpu_invy   , cpu_advz , cpu_invz , cpu_poisson , cpu_update
	real cpu_divchk , cpu_toxyz, cpu_init , cpu_output  , cpu_total_i , cpu_total
	real Re         , restart

        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	integer nsteps       ! number of steps to run
	real    timeLimit    ! stopping time

	integer ianim        ! whether to save data for animation
	real    saveInterval ! time interval to save data
	integer nsave        ! interval number for saving data

	integer idct         ! option for discrete cosine transform ((0) -> multigrid, (1) -> dct)
	integer iKEt         ! option for compute KE ((0) -> no, (1) -> yes)
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	real	:: TimeLoop_time

end module

!=======================================================================
subroutine simulation_init()
	use simulation

	cpu_cflcal =0.0
	cpu_bc     =0.0
	cpu_advx   =0.0
	cpu_invx   =0.0
	cpu_advy   =0.0
	cpu_invy   =0.0
	cpu_advz   =0.0
	cpu_invz   =0.0
	cpu_poisson=0.0
	cpu_update =0.0
	cpu_divchk =0.0
	cpu_toxyz  =0.0
	cpu_init   =0.0
	cpu_output =0.0
	cpu_total_i=0.0
	cpu_total  =0.0
	t_init     = realClock()
	t_over     = realClock() - t_init
	before     = realClock()

	call arrayzero
	call readFloat("Re",Re)
	visc = 1.0/Re

	call readInt("nsteps", nsteps)
	call readInt("ipout", ipout)
	call readFloat("timeLimit", timeLimit)

	dt    = 0.
	stime = 0.
	ntime = 0

	call readInt("ianim", ianim)
	if (ianim > 0) call readFloat("save_dt", saveInterval)
	nsave = 0

	call readInt("idct", idct)
	if (idct == 1) then
		call poisson_init_dct()  ! Discrete cosine transform
		write(6,*) 'FFT in Z and X directions!!!'
	else
		call poisson_init()      ! Multigrid
	end if
	
	call readInt("iKEt", iKEt)

	call CFLcontrol_init()
	call advance_init()

	return
end

subroutine simulation_read()
	use simulation
	call readFloat("time", stime)
	call readInt("ntime", ntime)

	!T  Gets (dt)
	call CFLcontrol_read()

	return
end

subroutine simulation_write()
	use simulation
	call writeInt("nsteps", nsteps)
	call writeFloat("time", stime)
	call writeInt("ntime", ntime)
	call writeFloat("timeLimit", timeLimit)
	call writeInt("ianim", ianim)

	call writeFloat("Re",Re)

	call advance_write()
	call CFLcontrol_write()

	return
end

!=======================================================================
subroutine simulation_run()
	use simulation
#if USE_COUPLER
    use continuum_coupler_socket
#endif

	stime_ = stime
	ntime_ = ntime

	TimeLoop_time = realClock()

	do while (ntime < ntime_ + nsteps)
		!--------------------
		!    Compute CFL
		!--------------------
		before = realClock()
		if(mod(ntime-ntime_,ipout).eq.0) &
			call CFLcontrol_maxCFL()
		after = realClock()
		cpu_cflcal= (after - before) - t_over  

#if USE_COUPLER
		call socket_coupler_send_velocity ! send velocities needed to set the continuum forces in MD
#endif

		!----------------------------------------
		!    Set Boundary Conditions
		!----------------------------------------
		before = realClock()
			call timeDependent_Inlet_BC(dt)		!--- Inlet FST + Convective outflow
			call timeDependent_Wall_BC()		!--- Wall oscillation
			call CartesianBC(dt)			!--- No Slip   + Blasius upper surface
			call FluxBC()				!--- From cartesian
			call CopyBC_out()			!--- Copy BC to free (uc,vc,wc)
		after = realClock()
		cpu_bc = (after - before) - t_over

		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! I. Advance the convection, explicit diffusion and
		!    implicit diffusion terms to get intermediate velocity
		! Note: 'invertu' saves the results of (ust) in (uc)
		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		!------------- u-equation -------------
		before = realClock()
			call advancex
		after = realClock()
		cpu_advx = (after - before) - t_over

		before = realClock()  
			call invertu
		after = realClock()
		cpu_invx = (after - before) - t_over

		!------------- v-equation -------------
		before = realClock()
			call advancey
		after = realClock()
		cpu_advy = (after - before) - t_over
        
		before = realClock()
			call invertv
		after = realClock()
		cpu_invy = (after - before) - t_over
      
		!------------- w-equation -------------

		before = realClock()
			call advancez
		after = realClock()
		cpu_advz = (after - before) - t_over

		before = realClock()
			call invertw
		after = realClock() 
		cpu_invz = (after - before) - t_over

		! do j=j1_w,j2_w
		! do i=i1_w,i2_w
		! do k=0,ngz+1
		!	if (abs(w  (k,i,j)).gt.0.0) print*, '(i,j,k, w  ) = ', i,j,k, w  (k,i,j)
		!	if (abs(wst(k,i,j)).gt.0.0) print*, '(i,j,k, wst) = ', i,j,k, wst(k,i,j)
		! enddo
		! enddo
		! enddo

		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	II.	POISSON	 SOLVER
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!------------ B.C. on (ust)  -------------
		call poissonRhsBC()

		!--------------- Solver ------------------
		before = realClock()
		if (idct == 1) then
			!-------------------------------------------
			! FFT in z-dir. + Cosine transform in x-dir.
			!-------------------------------------------
			call poisson_fftz_dctx(cpu_poisson)
		else
			!------------------------------------------
			! FFT in z-dir. + Multigrid in x & y dir.
			!------------------------------------------
			!call poisson_fftz_mtgxy
			call poisson_fftz_mtgxy(cpu_poisson)
		end if
		after = realClock()
		cpu_poisson = (after - before) - t_over
        
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	III.	Projection Step
		!ccccccccccccccccccccccccccccccccccccccccccccccccc
		!	Update (u,v,w) from (ust,vst,wst) & P
		!-------------------------------------------------
		before = realClock()
			call Update_Flux()
		after = realClock()
		cpu_update = (after - before) - t_over
        
		before = realClock()
			call Divergence_Check()
		after = realClock()
		cpu_divchk = (after - before) - t_over

		!------- Compute Cartesian Vel. --------
		before = realClock()   
			call CopyBC_in()		!--- Copy BC to free (uc,vc,wc)
			call Flux_to_Cart()
		after = realClock()
		cpu_toxyz = (after - before) - t_over

		!------- Compute flowAnalysis -----
	        ! Write flow analysis file
	        call flowAnalysis_blayerinfo()

		!------- Compute Statistics --------
		call statistics_sample

		!---------------------------
		!  INCREMENT TIME COUNTERS
		!---------------------------
		stime = stime + dt
		ntime = ntime + 1

		!cccccccccccccccccccccccccccccccccccccccccc
		!		ANIMATION
		!cccccccccccccccccccccccccccccccccccccccccc
		!------ Write out data for animation-------
		if (ianim == 1) then
			if (nint(stime/saveInterval) >= nsave) then
				nsave = nint(stime/saveInterval) + 1
				!JACOBS	call data_writeZdata()
 				!	call WriteZplane()
			end if
		end if
		if (ianim == 2) then
			if (nint(stime/saveInterval) >= nsave) then
				nsave = nint(stime/saveInterval) + 1
				!JACOBS	call data_writeYdata()
 				!	call WriteYplane()
			end if
		end if
		if (ianim == 12) then
			if (nint(stime/saveInterval) >= nsave) then
				nsave = nint(stime/saveInterval) + 1
				!JACOBS	call data_writeYdata()
 					! call WriteZplane()
 					! call WriteYplane()
					! call WriteSubDomain('SubDomYdble.', 1,ngx,1, 3, 21,2, 1,ngz,1)
					! call WriteSubDomain('SubDomZdble.', 1,ngx,1, 1,ngy,1, 1, 7 ,2)
					! call WriteSubDomain('SubDom_dble.', 1,ngx,2, 1,ngy,2, 1,ngz,2)
					!call WriteSubDomain('SubDom_dble.', 1,ngx,1, 1,ngy,1, 1,ngz,1)
					call WriteSubDomain_cc('SubDom_dble.', 1,ngx,1, 0,ngy,1, 1,ngz,1) !Cell centered version
			end if
		end if
		!	if (ianim == 3) then
		!		if (nint(stime/saveInterval) >= nsave) then
		!			nsave = nint(stime/saveInterval) + 1
		!			!JACOBS	call data_writeXdata()
		!		end if
		!	end if

		!cccccccccccccccccccccccccccccccccccccccccc
		!	Checks for Stopping the code
		!cccccccccccccccccccccccccccccccccccccccccc
		!--------- Check for stopping time --------
		if (timeLimit .ne. 0. .and. stime >= timeLimit) then
			call signalShell("done")
			exit
		end if

		!------- Check for user intervention ------
		if(mod(ntime-ntime_,20).eq.0) then
			inquire(file="resubmit", exist=iexist)
			if (iexist) then
				call signalShell("resubmit")
				! Remove file (only one process actually does this)
				! call system('rm resubmit')
				exit
			end if
			inquire(file="done", exist=iexist)
			if (iexist) then
				call signalShell("done")
				! Remove file (only one process actually does this)
				! call system('rm done')
				exit
			end if
		end if

		!--------- Check for CPU time --------
		tremain = timeRemaining()
		call globalMax(tremain, 1)
		if (tremain < 0.) then
			call signalShell("resubmit")
			exit
		end if

		cpu_total_i = cpu_cflcal+cpu_bc+cpu_advx+cpu_invx+cpu_advy+cpu_invy+ &
		cpu_advz+cpu_invz+cpu_update+cpu_divchk+cpu_toxyz+cpu_poisson

		if (irank.eq.1) then
			write(6,'(A70,A70)') & 
				'CPU cflcal    bc       advx      invx      advy      invy      advz   ,', &
				'   invz     update    divchk     toxyz  tot-poison poisson   total_i  '
			write(6,'(14(1PE10.3))') &
				cpu_cflcal,cpu_bc,cpu_advx,cpu_invx,cpu_advy,cpu_invy,cpu_advz, &
				cpu_invz,cpu_update,cpu_divchk,cpu_toxyz,cpu_total_i-cpu_poisson,cpu_poisson,cpu_total_i
		end if

		cpu_total = cpu_total + cpu_total_i

	end do
           
	print*,'Time loop timing = ', TimeLoop_time - realClock()

	if (ntime >= ntime_ + nsteps) then
		if (timeLimit .ne. 0. .and. stime + 0.5*dt < timeLimit) then
			call signalShell("continue")
		else
			call signalShell("done")
		end if
	end if

	if (irank.eq.iroot) then 
		write(6,'(A13,1PE10.3)') 'CPU total  = ',cpu_total
	end if

	return
end


!=======================================================================
subroutine simulation_info(iunit)
	use simulation

	call sectionTitle(iunit, "Simulation info")
	write (iunit, 11) ntime-ntime_, stime_, stime, stime-stime_

11  format (4x, "Number of time steps this run : ", I13 / &
	        4x, "Time at beginning of run      : ", f13.6 / &
	        4x, "Time at end of run            : ", f13.6 / &
	        4x, "Time advancement              : ", f13.6)

	!T	call equations_info(iunit)

	return
end

