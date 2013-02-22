!-----------------------------------------------------------------------------
!
!                              Record Final Properties
! Output properties from simualtion output to file
!
!-----------------------------------------------------------------------------

module module_final_record

	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use calculated_properties_MD
	use polymer_info_MD

end module module_final_record
!----------------------------------------------------------------------------------

subroutine finish_final_record
use module_final_record
use interfaces, only: error_abort
implicit none

	call simulation_compute_forces
	call evaluate_macroscopic_properties
	if (irank .eq. iroot) then
		select case(integration_algorithm)
		case(leap_frog_verlet)
			print('(a,i8,a)'), ' Results from final state of simulation at iter ',iter-1,':'
		case(velocity_verlet)
			print('(a,i8,a)'), ' Results from final state of simulation at iter ',iter,':'
		end select
	end if

   select case(integration_algorithm)
   case(leap_frog_verlet)
       call print_macroscopic_properties(Nsteps)
   case(velocity_verlet)
       call print_macroscopic_properties(Nsteps)
   end select
	
	!Write values of distribution functions
	!write(12,'(a)') 'Velocity frequency distribution'
	!do n=1,nbins(1)
	!	write(12,'(2(f10.5))') (n-0.5d0)*binsize, normalisedbin(n) 
	!enddo

	call messenger_syncall() !Make sure all processes have finished writing

	!Write simualtion properties and final position & velocity to 
	!unformatted output file to allow restart
	call parallel_io_final_state

	!Reformat positions recorded into the correct form for VMD to use
	if(irank .eq. iroot) then
		select case (vmd_outflag)
		case(0)
		case(1)
			call reformat_dcd
		case(2)
			call reformat_dcd_sl
		case(3)
			call reformat_dcd
			call reformat_dcd_halo
		case(4)
			call reformat_dcd_true
		case default
		end select 
	endif

	!Close all output files
	if (irank .eq. iroot) then   !Close Pressure tensor and viscosity output file
		if (pressure_outflag .ne. 0) close(2,status='keep')   !Close Pressure tensor and viscosity output file
	endif
	close(4,status='keep') !Close Control volume output file

	!close(12,status='keep') !Close statistics output file
	!close(14,status='keep') !Close velocity contour output file
	if (nproc .eq. 1) close(18,status='keep') !Close velocity slice output

	!Final timing printout
	!if (irank .eq. iroot) then
	!	print*, 'Elapsed time so far is', elapsedtime, 'timesteps'
	!	print*, '                  END OF SIMULATION                         '
	!endif

	select case (macro_outflag)
	case(2,4)
		close(10,status='keep') !Keep macroscopic_properties
	case default
	end select

	!Write config to be used in DL_POLY input
	!call write_DL_POLY_config
	
end subroutine finish_final_record

!--------------------------------------------------------------------------------------------
!Re-format the file into binary with header used by .dcd file format

!Header added by David Trevelyan, DCD file format from:
!http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
!Information at the above url is a little misleading in places, comments below try to clarify. 

subroutine reformat_dcd
	use module_final_record
	use computational_constants_MD
	implicit none

	integer							:: i			!--Dummy variables
	integer                         :: plot_mod
	integer							:: NSET,vmd_sets!--Number of frames
	integer							:: ISTRT		!--Starting frame
	integer							:: NSAVC		!--Number of frames per coordinate save
	integer							:: NATOMNFREAT		!--Number of fixed atoms
	integer							:: NTITLE		!--Number of 80-character strings in title (set as 2)
	integer							:: NATOM		!--Number of atoms
	integer, parameter 				:: LongInt = selected_int_kind (8)
	integer(kind=LongInt)			:: bufsize, starti, endi
	integer, dimension (5)			:: FIVEZ		!--According to documentation, five zeros, but first is actually NSET
	integer, dimension (9)			:: NINEZ		!--Nine zeros
	character(len=4)				:: HDR			!--Header string, value 'CORD'
	character(len=80), dimension(2)	:: TITLE	        !--Title(s)
	real 							:: time_start, time_end
	real,allocatable,dimension(:)   :: Xbuf, Ybuf, Zbuf	!--Buffers used to copy from direct access to binary
	double precision				:: DELTA		!--Timestep between frames

	print*, 'Initialising trajectory file reformat to *.dcd. For large'
	print*, 'systems or long runs this may take some time...' 
	
	call cpu_time(time_start)

	!Determine size of file datatype
	!inquire(file='testfile.dcd', recl=datasize)
	!print*, 'datasize', datasize
 	!call MPI_type_size(MPI_real,datasize,ierr)
	!print*, 'datasize', datasize

	if (Nvmd_intervals.eq.0) then
		vmd_sets = (Nsteps-initialstep+1)/tplot
	else
		vmd_sets = vmd_count-1
	endif

	!Set header information	
	HDR			=	'CORD'				!header text
	NSET		=	vmd_sets			!number of recorded frames
	ISTRT		=	initialstep			!the starting timestep
	NSAVC		=	tplot				!number of timesteps between dcd frame saves
	FIVEZ(1)	=	NSET				!not sure why
	FIVEZ(2:5)	=	0					!buffer zeros
	NATOMNFREAT	=	0					!number of fixed atoms?
	DELTA		=	delta_t				!delta_t (x-plor is double, charmm is real)
	NINEZ(:)	=	0					!buffer zeros
	NTITLE		=	2					!number of 80-character strings in title
	TITLE(1)	=	'  Simulation record file '	!
	TITLE(2)    =	'   Written in serial or parallel   '	!
	NATOM		=	globalnp			!number of particles

	!Get size of arrays required to store data
	bufsize = NSET*globalnp

	!Avoid reformat if filesize is too large
	if (bufsize .gt. 2147483647) then
		print*, 'File size too large for VMD reformat - reformat vmd_temp.dcd locally"
		return
	else
		
		allocate(Xbuf(bufsize))
		allocate(Ybuf(bufsize))
		allocate(Zbuf(bufsize))

		!Read position information from file
		!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
		open (unit=17, file=trim(prefix_dir)//"results/vmd_temp.dcd",access='stream')
	    
		!Open unit 6 (stdout) with fortran carriage control 
		open (unit=6, carriagecontrol='fortran')  
		plot_mod = max(1,NSET/100)
		write(*,'(a)') ' VMD reformat read completion:  '

		!Read temp trajectory file
		do i=1,NSET
			starti = globalnp*(i-1)+1
			endi   = globalnp*i
			read(18) Xbuf(starti:endi)
			read(18) Ybuf(starti:endi)
			read(18) Zbuf(starti:endi)
			!read(18) Xbuf(globalnp*(i-1)+1:globalnp*i)
			!read(18) Ybuf(globalnp*(i-1)+1:globalnp*i)
			!read(18) Zbuf(globalnp*(i-1)+1:globalnp*i)
			if (mod(i,plot_mod) .eq. 0) then
				call progress(100*i/NSET)
			end if
		enddo
		close(17,status='delete')

		!Open binary .dcd file and write header information	
		open(unit=3, file=trim(prefix_dir)//"results/vmd_out.dcd",status='replace', form="unformatted")
		
		write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
		write(3) NTITLE, TITLE(1), TITLE(2)
		write(3) NATOM

		write(*,'(a)') ' VMD reformat write completion: '
		do i=1,NSET
			starti = (i-1)*globalnp+1
			endi   = i*globalnp
			write(4) Xbuf(starti:endi)
			write(4) Ybuf(starti:endi)
			write(4) Zbuf(starti:endi)
			!write(4) Xbuf((i-1)*globalnp+1:i*globalnp)
			!write(4) Ybuf((i-1)*globalnp+1:i*globalnp)
			!write(4) Zbuf((i-1)*globalnp+1:i*globalnp)
			if (mod(i,plot_mod) .eq. 0) then
				call progress(100*i/NSET)
			end if
		enddo

		close(3,status='keep')

		deallocate(Xbuf)
		deallocate(Ybuf)
		deallocate(Zbuf)

		call cpu_time(time_end)

	 	print '(a,g10.2,a)', ' Reformatted to *.dcd in', time_end - time_start, ' seconds.'
	endif

end subroutine reformat_dcd

subroutine progress(j)  
implicit none  

	integer(kind=4)   :: j,k,bark
	character(len=78) :: bar

	bar(1:7) = " ???% |"
	do k=8,77
		bar(k:k) = " "
	end do
	bar(78:78) = "|"

	write(unit=bar(2:4),fmt="(i3)") j

	bark = int(dble(j)*70.d0/100.d0)

	do k=1,bark
		!bar(7+k:7+k)=char(2) 
		bar(7+k:7+k)="="
	enddo  

	! print the progress bar.  
	write(unit=6,fmt="(a1,a1,a78)") '+',char(13),bar  

	return  

end subroutine progress 

!--------------------------------------------------------------------------------------------
!Re-format the file into binary with header used by .dcd file format

!Header added by David Trevelyan, DCD file format from:
!http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
!Information at the above url is a little misleading in places, comments below try to clarify. 

subroutine reformat_dcd_true
	use module_final_record
	use computational_constants_MD
	implicit none

	integer							:: i			!--Dummy variables
	integer                         :: plot_mod
	integer							:: NSET,vmd_sets!--Number of frames
	integer							:: ISTRT		!--Starting frame
	integer							:: NSAVC		!--Number of frames per coordinate save
	integer							:: NATOMNFREAT	!--Number of fixed atoms
	integer							:: NTITLE		!--Number of 80-character strings in title (set as 2)
	integer							:: NATOM		!--Number of atoms
	integer, parameter 				:: LongInt = selected_int_kind (8)
	integer(kind=LongInt)			:: bufsize, starti, endi
	integer, dimension (5)			:: FIVEZ		!--According to documentation, five zeros, but first is actually NSET
	integer, dimension (9)			:: NINEZ		!--Nine zeros
	character(len=4)				:: HDR			!--Header string, value 'CORD'
	character(len=80), dimension(2)	:: TITLE	        !--Title(s)
	real 							:: time_start, time_end
	real,allocatable,dimension(:)   :: Xbuf, Ybuf, Zbuf	!--Buffers used to copy from direct access to binary
	double precision				:: DELTA		!--Timestep between frames


	print*, 'Initialising trajectory file reformat to *.dcd. For large'
	print*, 'systems or long runs this may take some time...' 
	
	call cpu_time(time_start)

	!Determine size of file datatype
	!inquire(file='testfile.dcd', recl=datasize)
	!print*, 'datasize', datasize
 	!call MPI_type_size(MPI_real,datasize,ierr)
	!print*, 'datasize', datasize

	if (Nvmd_intervals.eq.0) then
		vmd_sets = (Nsteps-initialstep+1)/tplot
	else
		vmd_sets = vmd_count-1
		!do i = 1,size(vmd_intervals,2)
		!	vmd_sets = vmd_sets + (vmd_intervals(2,i) - vmd_intervals(1,i))
		!enddo
	endif

	!Set header information	
	HDR			=	'CORD'				!header text
	NSET		=	vmd_sets			!number of recorded frames
	ISTRT		=	initialstep			!the starting timestep
	NSAVC		=	tplot				!number of timesteps between dcd frame saves
	FIVEZ(1)	=	NSET				!not sure why
	FIVEZ(2:5)	=	0					!buffer zeros
	NATOMNFREAT	=	0					!number of fixed atoms?
	DELTA		=	delta_t				!delta_t (x-plor is double, charmm is real)
	NINEZ(:)	=	0					!buffer zeros
	NTITLE		=	2					!number of 80-character strings in title
	TITLE(1)	=	'  Simulation record file '	!
	TITLE(2)    =	'   Written in serial or parallel   '	!
	NATOM		=	globalnp			!number of particles

	bufsize = NSET*globalnp
	allocate(Xbuf(bufsize))
	allocate(Ybuf(bufsize))
	allocate(Zbuf(bufsize))

	!Read position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=18, file=trim(prefix_dir)//"results/vmd_temp_true.dcd",access='stream')

	!Open unit 6 (stdout) with fortran carriage control 
	open (unit=6, carriagecontrol='fortran')  
	plot_mod = max(1,NSET/100)
	write(*,'(a)') ' VMD reformat read completion:  '

	!Read temp trajectory file
	do i=1,NSET
		starti = globalnp*(i-1)+1
		endi   = globalnp*i
		read(18) Xbuf(starti:endi)
		read(18) Ybuf(starti:endi)
		read(18) Zbuf(starti:endi)
		!read(18) Xbuf(globalnp*(i-1)+1:globalnp*i)
		!read(18) Ybuf(globalnp*(i-1)+1:globalnp*i)
		!read(18) Zbuf(globalnp*(i-1)+1:globalnp*i)
		if (mod(i,plot_mod) .eq. 0) then
			call progress(100*i/NSET)
		end if
	enddo

	close(18,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=4, file=trim(prefix_dir)//"results/vmd_out_true.dcd",status='replace', form="unformatted")
	
	write(4) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(4) NTITLE, TITLE(1), TITLE(2)
	write(4) NATOM

	write(*,'(a)') ' VMD reformat write completion: '
	do i=1,NSET
		starti = (i-1)*globalnp+1
		endi   = i*globalnp
		write(4) Xbuf(starti:endi)
		write(4) Ybuf(starti:endi)
		write(4) Zbuf(starti:endi)
		!write(4) Xbuf((i-1)*globalnp+1:i*globalnp)
		!write(4) Ybuf((i-1)*globalnp+1:i*globalnp)
		!write(4) Zbuf((i-1)*globalnp+1:i*globalnp)
		if (mod(i,plot_mod) .eq. 0) then
			call progress(100*i/NSET)
		end if
	enddo

	close(4,status='keep')

	deallocate(Xbuf)
	deallocate(Ybuf)
	deallocate(Zbuf)

	call cpu_time(time_end)

 	print '(a,g10.2,a)', ' Reformatted to *.dcd in', time_end - time_start, ' seconds.'

end subroutine reformat_dcd_true

!--------------------------------------------------------------------------------------------
! SPLIT SOLID AND LIQUID REGIONS INTO SEPERATE MOLECULE FILES
!Re-format the file into binary with header used by .dcd file format

!Header added by David Trevelyan, DCD file format from:
!http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
!Information at the above url is a little misleading in places, comments below try to clarify. 


subroutine reformat_dcd_sl
	use module_final_record
	use computational_constants_MD
	implicit none

	integer							:: n, i			!--Dummy variables
	integer                         :: plot_mod
	integer							:: NSET,vmd_sets!--Number of frames
	integer							:: ISTRT		!--Starting frame
	integer							:: NSAVC		!--Number of frames per coordinate save
	integer							:: NATOMNFREAT		!--Number of fixed atoms
	integer							:: NTITLE		!--Number of 80-character strings in title (set as 2)
	integer							:: NATOM		!--Number of atoms
	integer, dimension (5)			:: FIVEZ		!--According to documentation, five zeros, but first is actually NSET
	integer, dimension (9)			:: NINEZ		!--Nine zeros
	character(len=4)				:: HDR			!--Header string, value 'CORD'
	character(len=80), dimension(2)	:: TITLE	        !--Title(s)
	real 							:: time_start, time_end
	real,allocatable,dimension(:)   :: Xbuf, Ybuf, Zbuf	!--Buffers used to copy from direct access to binary
	double precision				:: DELTA		!--Timestep between frames

	print*, 'Initialising trajectory file reformat to *.dcd. For large'
	print*, 'systems or long runs this may take some time...' 
	
	call cpu_time(time_start)

	if (Nvmd_intervals.eq.0) then
		vmd_sets = (Nsteps-initialstep+1)/tplot
	else
		vmd_sets = vmd_count-1
		!do i = 1,size(vmd_intervals,2)
		!	vmd_sets = vmd_sets + (vmd_intervals(2,i) - vmd_intervals(1,i))
		!enddo
	endif

	HDR			=	'CORD'				!header text
	NSET		=	vmd_sets			!number of recorded frames
	ISTRT		=	initialstep			!the starting timestep
	NSAVC		=	tplot				!number of timesteps between dcd frame saves
	FIVEZ(1)	=	NSET				!not sure why
	FIVEZ(2:5)	=	0					!buffer zeros
	NATOMNFREAT	=	0					!number of fixed atoms?
	DELTA		=	delta_t				!delta_t (x-plor is double, charmm is real)
	NINEZ(:)	=	0					!buffer zeros
	NTITLE		=	2					!number of 80-character strings in title
	TITLE(1)	=	'  Simulation record file '	!
	TITLE(2)    =	'   Written in serial/parallel   '	!
	NATOM		=	globalnp			!number of particles

	allocate(Xbuf(NSET*globalnp))
	allocate(Ybuf(NSET*globalnp))
	allocate(Zbuf(NSET*globalnp))

	!Read Solid molecule position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file=trim(prefix_dir)//"results/vmd_solid_temp.dcd",access='direct',recl=1)
	
	!Open unit 6 (stdout) with fortran carriage control 
	open (unit=6, carriagecontrol='fortran')  
	plot_mod = max(1,NSET/100)
	write(*,'(a)') ' VMD reformat read completion (solid):  '

	do i=1,NSET
		do n=1,globalnp
			read(17,rec=(i-1)*nd*globalnp+n) Xbuf(n+globalnp*(i-1))
			read(17,rec=(i-1)*nd*globalnp+n+globalnp) Ybuf(n+globalnp*(i-1))
			read(17,rec=(i-1)*nd*globalnp+n+2*globalnp) Zbuf(n+globalnp*(i-1))
		enddo
		if (mod(i,plot_mod) .eq. 0) then
			call progress(100*i/NSET)
		end if
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file=trim(prefix_dir)//"results/vmd_solid_out.dcd",status='replace', form="unformatted")
	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	write(*,'(a)') ' VMD reformat write completion (solid):  '
	do i=1,NSET
		write(3) Xbuf((i-1)*globalnp+1:i*globalnp)
		write(3) Ybuf((i-1)*globalnp+1:i*globalnp)
		write(3) Zbuf((i-1)*globalnp+1:i*globalnp)
		if (mod(i,plot_mod) .eq. 0) then
			call progress(100*i/NSET)
		end if
	enddo

	close(3,status='keep')

	!Read liquid molecule position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file=trim(prefix_dir)//"results/vmd_liquid_temp.dcd",access='direct',recl=1)

	write(*,'(a)') ' VMD reformat read completion (liquid):  '
	do i=1,NSET
		do n=1,globalnp
			read(17,rec=(i-1)*nd*globalnp+n) 			Xbuf(n+globalnp*(i-1))
			read(17,rec=(i-1)*nd*globalnp+n+globalnp) 	Ybuf(n+globalnp*(i-1))
			read(17,rec=(i-1)*nd*globalnp+n+2*globalnp) Zbuf(n+globalnp*(i-1))
		enddo
		if (mod(i,plot_mod) .eq. 0) then
			call progress(100*i/NSET)
		end if
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file=trim(prefix_dir)//"results/vmd_liquid_out.dcd",status='replace', form="unformatted")
	write(*,'(a)') ' VMD reformat write completion (liquid):  '
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	do i=1,NSET
		write(3) Xbuf((i-1)*globalnp+1:i*globalnp)
		write(3) Ybuf((i-1)*globalnp+1:i*globalnp)
		write(3) Zbuf((i-1)*globalnp+1:i*globalnp)
		if (mod(i,plot_mod) .eq. 0) then
			call progress(100*i/NSET)
		end if
	enddo

	close(3,status='keep')
	
	deallocate(Xbuf)
	deallocate(Ybuf)
	deallocate(Zbuf)

	call cpu_time(time_end)


 	print '(a,g10.2,a)', ' Reformatted to *.dcd in', time_end - time_start, ' seconds.'

end subroutine reformat_dcd_sl

!--------------------------------------------------------------------------------------------
!Re-format halo cell output for vmd ~ useful for debuggin purposes

subroutine reformat_dcd_halo
	use module_final_record
	implicit none

	integer							:: n, i			!--Dummy variables
	integer                         :: plot_mod
	integer							:: NSET,vmd_sets!--Number of frames
	integer							:: ISTRT		!--Starting frame
	integer							:: NSAVC		!--Number of frames per coordinate save
	integer							:: NATOMNFREAT		!--Number of fixed atoms
	integer							:: NTITLE		!--Number of 80-character strings in title (set as 2)
	integer							:: NATOM		!--Number of atoms
	integer, dimension (5)			:: FIVEZ		!--According to documentation, five zeros, but first is actually NSET
	integer, dimension (9)			:: NINEZ		!--Nine zeros
	character(len=4)				:: HDR			!--Header string, value 'CORD'
	character(len=80), dimension(2)	:: TITLE	        !--Title(s)
	real 							:: time_start, time_end
	real,allocatable,dimension(:)   :: Xbuf, Ybuf, Zbuf	!--Buffers used to copy from direct access to binary
	double precision				:: DELTA		!--Timestep between frames

	print*, 'Initialising trajectory file reformat to *.dcd. For large'
	print*, 'systems or long runs this may take some time...' 
	
	call cpu_time(time_start)

	!Determine size of file datatype
	!inquire(file='testfile.dcd', recl=datasize)
	!print*, 'datasize', datasize
 	!call MPI_type_size(MPI_real,datasize,ierr)
	!print*, 'datasize', datasize

	if (Nvmd_intervals.eq.0) then
		vmd_sets = (Nsteps-initialstep+1)/tplot
	else
		vmd_sets = vmd_count-1
		!do i = 1,size(vmd_intervals,2)
		!	vmd_sets = vmd_sets + (vmd_intervals(2,i) - vmd_intervals(1,i))
		!enddo
	endif

	!Set header information	
	HDR			=	'CORD'				!header text
	NSET		=	vmd_sets			!number of recorded frames
	ISTRT		=	0					!the starting timestep
	NSAVC		=	1					!number of timesteps between dcd frame saves
	FIVEZ(1)	=	NSET				!not sure why
	FIVEZ(2:5)	=	0					!buffer zeros
	NATOMNFREAT	=	0					!number of fixed atoms?
	DELTA		=	delta_t				!delta_t (x-plor is double, charmm is real)
	NINEZ(:)	=	0					!buffer zeros
	NTITLE		=	2					!number of 80-character strings in title
	TITLE(1)	=	'  Simulation record file '	!
	TITLE(2)    =	'   Written in parallel   '	!
	NATOM		=	extralloc			!number of particles

	allocate(Xbuf(NSET*extralloc))
	allocate(Ybuf(NSET*extralloc))
	allocate(Zbuf(NSET*extralloc))

	!Read Solid molecule position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file=trim(prefix_dir)//"results/vmd_halo_temp.dcd",access='direct',recl=1)
	
	!Open unit 6 (stdout) with fortran carriage control 
	open (unit=6, carriagecontrol='fortran')  
	plot_mod = max(1,NSET/100)
	write(*,'(a)') ' VMD reformat read completion (halo):  '

	do i=1,NSET
		do n=1,extralloc
			read(17,rec=(i-1)*nd*extralloc+n) Xbuf(n+extralloc*(i-1))
			read(17,rec=(i-1)*nd*extralloc+n+extralloc) Ybuf(n+extralloc*(i-1))
			read(17,rec=(i-1)*nd*extralloc+n+2*extralloc) Zbuf(n+extralloc*(i-1))
		enddo
		if (mod(i,plot_mod) .eq. 0) then
			call progress(100*i/NSET)
		end if
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file=trim(prefix_dir)//"results/vmd_halo_out.dcd",status='replace', form="unformatted")
	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	write(*,'(a)') ' VMD reformat write completion (halo):  '
	do i=1,NSET
		write(3) Xbuf((i-1)*extralloc+1:i*extralloc)
		write(3) Ybuf((i-1)*extralloc+1:i*extralloc)
		write(3) Zbuf((i-1)*extralloc+1:i*extralloc)
		if (mod(i,plot_mod) .eq. 0) then
			call progress(100*i/NSET)
		end if
	enddo

	close(3,status='keep')

	deallocate(Xbuf)
	deallocate(Ybuf)
	deallocate(Zbuf)

	call cpu_time(time_end)

 	print '(a,g10.2,a)', ' Reformatted to *.dcd in', time_end - time_start, ' seconds'

end subroutine reformat_dcd_halo

!==============================================================================================================

!----------------------------------------------------------------------------------
!Write DL_POLY style CONFIG file of molecules in simulation
subroutine write_DL_POLY_config
use module_final_record
implicit none

	integer			i

	open(1000,file=trim(prefix_dir)//'CONFIG')
	write(1000,*) 'Argon'					!For Argon LJ units
	write(1000,*) 2, 1, np					!2 ~ use position, velocity & acceleration input, 1 3D cube periodic BC, np
	write(1000,*) 3.4*domain(1),	0.d0,		0.d0		!Domain size in Angstroms
	write(1000,*) 0.d0,			3.4*domain(2),	0.d0
	write(1000,*) 0.d0,			0.d0,		3.4*domain(3)
	do i=1,np
		write(1000,*) 'Ar', i
		write(1000,*) 3.4*r(:,i) !Convert to Angstroms
		write(1000,*) (3.4/2.16) *v(:,i)	!Convert to Angstroms/ps
		write(1000,*) a(:,i)
	enddo

	close(1000,status='keep')

end subroutine write_DL_POLY_config


