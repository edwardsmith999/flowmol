!-----------------------------------------------------------------------------
!
!                 Check_Rebuild and write rescue files                               
! Check to see if neighbourlist needs to be rebuit and returns flag true
! Also check if rescue file needs to be written
!
!-----------------------------------------------------------------------------
module module_checkrebuild

	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD

end module module_checkrebuild
!----------------------------------------------------------------------------------

subroutine simulation_checkrebuild(rebuild)
	use module_checkrebuild
	use interfaces, only: error_abort
	implicit none
	
	logical,save		   :: just_written_snapshot=.true.
    logical                :: abort
	integer                :: n, proc_abort
	integer, save		   :: rb_count, total_rb=0
	integer, intent(out)   :: rebuild
	double precision       :: vmax, t2, dt
	double precision, save :: rmax = 0.d0, t1, average_rb_count
	double precision,dimension(:),allocatable :: vmagnitude
	double precision,dimension(:,:),allocatable,save :: rdisp  !Displacement from initial state used in checkrebuild
	rebuild = 0

	!Trigger rebuild if record to be taken on next timestep
	if (mod(iter+1,tplot) .eq. 0) rebuild = 1

	!Check if backup file should be saved based on elapsed time
	!since last backup has been saved
	if (just_written_snapshot .eqv. .true.) then
		just_written_snapshot = .false.
		call cpu_time(t1)  !Get start cpu time
	endif
	call cpu_time(t2)  !Get current cpu time
	dt = t2-t1		   !Difference in times
	call globalMax(dt) !Ensure that all processors use the same time
	if (dt .gt. rescue_snapshot_freq) then
		if (irank.eq.iroot) then
			print('(2(a,f10.5))'), ' It`s been ', dt, & 
				' seconds since last rescue and requested freq is ', rescue_snapshot_freq
			print('(a,i8)'), ' Rescue microstate written to ./results/final_state at iter ', iter
		end if
		just_written_snapshot = .true.
		call messenger_syncall
		call parallel_io_final_state
	endif

    !Abort can be forced by creating file -- global check required to ensure all abort together
	!NOTE - there may still be a bug in this -- is an ABORTABORT file is created but any of the
	!processes has passed beyond this check then code could hang...

	inquire(file=trim(prefix_dir)//'ABORTABORT',exist=abort)
    if (abort) then
        proc_abort = 1
    else
        proc_abort = 0
    endif
	call globalMaxInt(proc_abort)
    if (proc_abort .eq. 1) abort = .true.
	if (abort) then
		print*, 'File ABORTABORT detected. Writing final_state file...'
		call messenger_syncall
		call parallel_io_final_state	
		call error_abort('Restart file written. Simulation aborted.')
	end if

    select case(rebuild_criteria)
    case(0) !Rapaport vmax criterion
	    !Evaluate velocity magnitude
	    allocate(vmagnitude(np))
	    do n = 1, np    ! Loop over all molecules
		    vmagnitude(n) = dot_product(v(:,n),v(:,n)) 
	    enddo
	    !Obtain maximum velocity
	    vmax = maxval(vmagnitude)
	    vmax = sqrt(vmax)
	    deallocate(vmagnitude)

	    !Calculate maximum possible displacement based on maximum velocity
	    rmax = rmax + vmax * delta_t

	    !Check if maximum displacment has exceeded extra distance
	    if (rmax .gt. 0.5d0*delta_rneighbr) rebuild = 1

    case(1) !fixed

		if (mod(iter,fixed_rebuild) .eq. 0) rebuild = 1

    case(2) !displacement
        if (allocated(rdisp) .and. size(rdisp,2) .ne. np) then
            deallocate(rdisp)
        endif
        if (.not. allocated(rdisp)) then
            allocate(rdisp(3,np)); rdisp = 0.d0
        endif
        ! Loop over all molecules and add displacement
	    do n = 1, np    
            rdisp(:,n) = rdisp(:,n) + v(:,n) * delta_t
            !rmax = max(rmax,rdisp(:,n))
            if (any(abs(rdisp(:,n)) .gt. 0.5d0*delta_rneighbr)) then
                !print'(2i6,4f10.5)', iter, n, abs(rdisp(:,n)), delta_rneighbr
                rebuild = 1;
                exit
            endif
        enddo
    end select

!	!Trigger rebuild if FIXED REBUILD specified in input
!	if (fixed_rebuild_flag .eq. 1) then
!		if (mod(iter,fixed_rebuild) .eq. 0) rebuild = 1
!		return
!	endif

!	!Evaluate velocity magnitude
!	allocate(vmagnitude(np))
!	do n = 1, np    ! Loop over all molecules
!		vmagnitude(n) = dot_product(v(:,n),v(:,n)) 
!	enddo
!	!Obtain maximum velocity
!	vmax = maxval(vmagnitude)
!	vmax = sqrt(vmax)
!	deallocate(vmagnitude)

!	!Calculate maximum possible displacement based on maximum velocity
!	rmax = rmax + vmax * delta_t

!	!Check if maximum displacment has exceeded extra distance
!	if (rmax .gt. 0.5d0*delta_rneighbr) rebuild = 1

	!Ensure all processor flags are the same - if any rebuild
	!then all should rebuild
	call globalMaxInt(rebuild)

	!If rebuild is happening ...
	if(rebuild .eq. 1) then
		!Set local processor rmax to zero as rebuild is global
		rmax = 0.d0
        if (rebuild_criteria .eq. 2) then
            deallocate(rdisp) 
        endif
		!Reset count and store average rebuild frequency
		average_rb_count = average_rb_count + dble(rb_count)
		rb_count = 1; total_rb = total_rb + 1
		if (irank .eq. iroot) then
			if (mod(total_rb,100) .eq. 0) then	!Print rebuild stats every 1000 rebuilds
		   		!open(30,file=trim(prefix_dir)//'/results/rebuild_stats',access='append')
				!write(30,'(a,i7,a,f10.5,a)'), 'Up to iteration ', iter, & 
				!	' On average, rebuild happened every ', average_rb_count/total_rb, ' steps'
				!close(30,status='keep')
				print('(a,i7,a,f10.5,a)'), 'Up to iteration ', iter, & 
					' On average, rebuild happened every ', average_rb_count/dble(total_rb), ' steps'
			endif
		endif
	else
		rb_count = rb_count + 1
	endif

end subroutine simulation_checkrebuild
