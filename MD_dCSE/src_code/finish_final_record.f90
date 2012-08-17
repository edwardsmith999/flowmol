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

	!Print out final results from the simulation
	if (irank .eq. iroot) then
		select case(potential_flag)
		case(0)
			select case(macro_outflag)
			case(1:2)
				print*, 'Results from last iteration of simulation'
				print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f10.4)', &
				Nsteps,';',vsum,';', v2sum,';', temperature,';', &
				kinenergy,';',potenergy,';',totenergy,';',pressure
			case(3:4)
				print*, 'Results from last iteration of simulation'
				print '(1x,i8,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4)', &
				Nsteps,';',vsum,';',temperature,';',&
				kinenergy,';',potenergy,';',totenergy,';',pressure
			case default
			end select
		case(1)
			select case(macro_outflag)
			case(1:2)
				print*, 'Results from last iteration of simulation'
				print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f10.4,a,f10.4,a,f10.4)', &
				Nsteps,';',vsum,';', v2sum,';', temperature,';', &
				kinenergy,';',potenergy,';',totenergy,';',pressure,';',etevtcf,';',R_g
			case(3:4)
				print*, 'Results from last iteration of simulation'
				print '(1x,i8,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f6.3,a,f6.2)', &
				Nsteps,';',vsum,';',temperature,';', &
				kinenergy,';',potenergy,';',totenergy,';',pressure,';',etevtcf,';',R_g
			case default
			end select
		case default
		end select
	end if
	
	!Set up final print out table
	!print '(8a)', &
	!'Final Iteration; Simulation time;  Mean Temp;   Mean KE; SD KE; Mean PE; SD PE; Mean TE; SD TE; Mean Pressure;'

	!Write values of distribution functions
	!write(12,'(a)') 'Velocity frequency distribution'
	!do n=1,nbins(1)
	!	write(12,'(2(f10.5))') (n-0.5d0)*binsize, normalisedbin(n) 
	!enddo

	call messenger_syncall() !Make sure all processes have finished writing

	!Write simualtion properties and final position & velocity to 
	!unformatted output file to allow restart
	call parallel_io_final_state

	!Reformat positions recorded into the correct from for VMD to use
	if(irank .eq. iroot) then
		select case (vmd_outflag)
		case(0)
		case(1)
			select case(potential_flag)
			case(0)
				call reformat_dcd
			case(1)
				call reformat_dcd_true
			case default
			end select
		case(2)
			select case(potential_flag)
			case(0)
				call reformat_dcd_sl
			case default
				call error_abort('Reformat_dcd_sl not available for polymers')
			end select
		case(3)
			select case(potential_flag)
			case(0)
				call reformat_dcd
				call reformat_dcd_halo
			case default
				call error_abort('Reformat_dcd_halo not available for polymers')
			end select
		case default
		end select 
	endif

	if (vmd_outflag.ne.0 .and. potential_flag.eq.1) call build_psf

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

	integer							:: n, i			!--Dummy variables
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

	print*, 'Generating final VMD.dcd ouput file - for large systems or'
	print*, 'long runs this may take some time'

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

	allocate(Xbuf(NSET*globalnp))
	allocate(Ybuf(NSET*globalnp))
	allocate(Zbuf(NSET*globalnp))

	!Read position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file=trim(prefix_dir)//"results/vmd_temp.dcd",access='stream')

	do i=1,NSET
		read(17) Xbuf(globalnp*(i-1)+1:globalnp*i)
		read(17) Ybuf(globalnp*(i-1)+1:globalnp*i)
		read(17) Zbuf(globalnp*(i-1)+1:globalnp*i)
		if (mod(i,100) .eq. 0) print*, 'Reading % complete =', (100.d0*i/NSET)
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file=trim(prefix_dir)//"results/vmd_out.dcd",status='replace', form="unformatted")
	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	do i=1,NSET
		write(3) Xbuf((i-1)*globalnp+1:i*globalnp)
		write(3) Ybuf((i-1)*globalnp+1:i*globalnp)
		write(3) Zbuf((i-1)*globalnp+1:i*globalnp)
		if (mod(i,100) .eq. 0) print*, 'Writing % complete =', (100.d0*i/NSET)
	enddo

	close(3,status='keep')

	deallocate(Xbuf)
	deallocate(Ybuf)
	deallocate(Zbuf)

	call cpu_time(time_end)

 	print '(a,g10.2,a)', 'Generated final VMD.dcd ouput file in', time_end - time_start, ' seconds'

end subroutine reformat_dcd

!--------------------------------------------------------------------------------------------
!Re-format the file into binary with header used by .dcd file format

!Header added by David Trevelyan, DCD file format from:
!http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
!Information at the above url is a little misleading in places, comments below try to clarify. 

subroutine reformat_dcd_true
	use module_final_record
	use computational_constants_MD
	implicit none

	integer							:: n, i			!--Dummy variables
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

	print*, 'Generating final VMD.dcd ouput file - for large systems or'
	print*, 'long runs this may take some time'

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

	allocate(Xbuf(NSET*globalnp))
	allocate(Ybuf(NSET*globalnp))
	allocate(Zbuf(NSET*globalnp))

	!Read position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=18, file=trim(prefix_dir)//"results/vmd_temp_true.dcd",access='stream')

	do i=1,NSET
		read(18) Xbuf(globalnp*(i-1)+1:globalnp*i)
		read(18) Ybuf(globalnp*(i-1)+1:globalnp*i)
		read(18) Zbuf(globalnp*(i-1)+1:globalnp*i)
		if (mod(i,100) .eq. 0) print*, 'Reading % complete =', (100.d0*i/NSET)
	enddo

	close(18,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=4, file=trim(prefix_dir)//"results/vmd_out_true.dcd",status='replace', form="unformatted")
	
	write(4) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(4) NTITLE, TITLE(1), TITLE(2)
	write(4) NATOM

	do i=1,NSET
		write(4) Xbuf((i-1)*globalnp+1:i*globalnp)
		write(4) Ybuf((i-1)*globalnp+1:i*globalnp)
		write(4) Zbuf((i-1)*globalnp+1:i*globalnp)
		if (mod(i,100) .eq. 0) print*, 'Writing % complete =', (100.d0*i/NSET)
	enddo

	close(4,status='keep')

	deallocate(Xbuf)
	deallocate(Ybuf)
	deallocate(Zbuf)

	call cpu_time(time_end)

 	print '(a,g10.2,a)', 'Generated final VMD.dcd ouput file in', time_end - time_start, ' seconds'

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

	print*, 'Generating final VMD.dcd ouput file - for large systems or'
	print*, 'long runs this may take some time'

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
	TITLE(2)    =	'   Written in serial/parallel   '	!
	NATOM		=	globalnp			!number of particles

	allocate(Xbuf(NSET*globalnp))
	allocate(Ybuf(NSET*globalnp))
	allocate(Zbuf(NSET*globalnp))

	!Read Solid molecule position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file=trim(prefix_dir)//"results/vmd_solid_temp.dcd",access='direct',recl=1)

	do i=1,NSET
	do n=1,globalnp
		read(17,rec=(i-1)*nd*globalnp+n) Xbuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+globalnp) Ybuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+2*globalnp) Zbuf(n+globalnp*(i-1))
		if (mod(i,100) .eq. 0) print*, 'Reading solid % complete =', (100.d0*i/NSET)
	enddo
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file=trim(prefix_dir)//"results/vmd_solid_out.dcd",status='replace', form="unformatted")
	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	do i=1,NSET
		write(3) Xbuf((i-1)*globalnp+1:i*globalnp)
		write(3) Ybuf((i-1)*globalnp+1:i*globalnp)
		write(3) Zbuf((i-1)*globalnp+1:i*globalnp)
		if (mod(i,100) .eq. 0) print*, 'Writing solid % complete =', (100.d0*i/NSET)
	enddo

	close(3,status='keep')

	!Read liquid molecule position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file=trim(prefix_dir)//"results/vmd_liquid_temp.dcd",access='direct',recl=1)

	do i=1,NSET
	do n=1,globalnp
		read(17,rec=(i-1)*nd*globalnp+n) 			Xbuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+globalnp) 	Ybuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+2*globalnp) Zbuf(n+globalnp*(i-1))
		if (mod(i,100) .eq. 0) print*, 'Reading liquid % complete =', (100.d0*i/NSET)
	enddo
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file=trim(prefix_dir)//"results/vmd_liquid_out.dcd",status='replace', form="unformatted")
	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	do i=1,NSET
		write(3) Xbuf((i-1)*globalnp+1:i*globalnp)
		write(3) Ybuf((i-1)*globalnp+1:i*globalnp)
		write(3) Zbuf((i-1)*globalnp+1:i*globalnp)
		if (mod(i,100) .eq. 0) print*, 'Writing liquid % complete =', (100.d0*i/NSET)
	enddo

	close(3,status='keep')
	
	deallocate(Xbuf)
	deallocate(Ybuf)
	deallocate(Zbuf)

	call cpu_time(time_end)

 	print '(a,g10.2,a)', 'Generated final VMD.dcd ouput file in', time_end - time_start, ' seconds'

end subroutine reformat_dcd_sl

!--------------------------------------------------------------------------------------------
!Re-format halo cell output for vmd ~ useful for debuggin purposes

subroutine reformat_dcd_halo
	use module_final_record
	implicit none

	integer							:: n, i			!--Dummy variables
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


	print*, 'Generating final VMD.dcd ouput file - for large systems or'
	print*, 'long runs this may take some time'

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

	do i=1,NSET
	do n=1,extralloc
		read(17,rec=(i-1)*nd*extralloc+n) Xbuf(n+extralloc*(i-1))
		read(17,rec=(i-1)*nd*extralloc+n+extralloc) Ybuf(n+extralloc*(i-1))
		read(17,rec=(i-1)*nd*extralloc+n+2*extralloc) Zbuf(n+extralloc*(i-1))
		if (mod(i,100) .eq. 0) print*, 'Reading % complete =', (100.d0*i/NSET)
	enddo
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file=trim(prefix_dir)//"results/vmd_halo_out.dcd",status='replace', form="unformatted")
	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	do i=1,NSET
		write(3) Xbuf((i-1)*extralloc+1:i*extralloc)
		write(3) Ybuf((i-1)*extralloc+1:i*extralloc)
		write(3) Zbuf((i-1)*extralloc+1:i*extralloc)
		if (mod(i,100) .eq. 0) print*, 'Writing % complete =', (100.d0*i/NSET)
	enddo

	close(3,status='keep')

	deallocate(Xbuf)
	deallocate(Ybuf)
	deallocate(Zbuf)

	call cpu_time(time_end)

 	print '(a,g10.2,a)', 'Generated final VMD.dcd ouput file in', time_end - time_start, ' seconds'

end subroutine reformat_dcd_halo

!==============================================================================================================

subroutine build_psf
	use module_final_record
	use polymer_info_MD
	implicit none

	integer :: i,j,n,item,molno,sc
	integer :: NTITLE, NATOM, NBONDS
	integer	:: write_items
	integer :: group, bin_expo
	integer, allocatable, dimension(:,:) :: bonds	
	integer, allocatable, dimension(:)	 ::	res_ID
	integer, allocatable, dimension(:)	 ::	glob_sc
	integer, allocatable, dimension(:,:) ::	glob_bf
	character(len=4), allocatable, dimension(:) :: seg_name, res_name, atom_name, atom_type
	real, allocatable, dimension(:)	:: charge, mass
	
	NTITLE = 1												! How many 'REMARKS' lines you want
	NATOM  = globalnp										! Determine total number of atoms
	NBONDS = 0
	do n=1,np
		do sc=1,nmonomers
            j=0
			group = ceiling(real(sc)/real(intbits))
			bin_expo = mod(sc,intbits)-1
			if (bin_expo.eq.-1) bin_expo = intbits - 1
            if ( btest(monomer(n)%bin_bflag(group),bin_expo) ) j = 1
            NBONDS = NBONDS + j
		end do
	end do
	call globalSumInt(NBONDS)
	NBONDS = int(NBONDS/2)
	
	print*, 'Generating polymer topology file polymer_topol.psf'

	allocate(seg_name(NATOM))								! Determine segment names for each atom
	allocate(res_ID(NATOM))									! Determine molecule ID for each atom
	allocate(glob_sc(NATOM))								! Determine molecule ID for each atom
	allocate(glob_bf(NATOM,4))                              ! Determine molecule ID for each atom
	allocate(res_name(NATOM))								! Determine name for each molecule
	allocate(atom_name(NATOM))								! Determine name for each atom
	allocate(atom_type(NATOM))								! Determine type for each atom
	allocate(charge(NATOM))									! Determine charge for each atom
	allocate(mass(NATOM))									! Determine mass of each atom

	res_ID(:) = 0
	glob_sc(:) = 0
	glob_bf(:,:) = 0
	do n=1,np
		molno            = monomer(n)%glob_no
		res_ID(molno)    = monomer(n)%chainID 
		glob_sc(molno)   = monomer(n)%subchainID
		glob_bf(molno,:) = monomer(n)%bin_bflag(:)
	end do

	call globalSumIntVect(res_ID,globalnp)
	call globalSumIntVect(glob_sc,globalnp)
	call globalSumIntTwoDim(glob_bf,globalnp,4)

	do n=1,globalnp
		select case (res_ID(n))
		case(0)
			seg_name(n)  = 'N'
			res_name(n)  = 'SOL'
			atom_name(n) = 'N'
			atom_type(n) = 'N'
			mass(n)      = 1.00794
		case(1:)
			seg_name(n)  = 'C'
			res_name(n)  = 'POL'
			atom_name(n) = 'C'
			atom_type(n) = 'C'
			mass(n)      = 1.00794
		case default
		end select
		charge(n)    = 0.00000
	end do

	if (irank.eq.iroot) then
		
		open(unit=1, file='results/polymer_topol.psf', status='replace', form='formatted')

		! Header
		write(1,'(a3)') 'PSF'
		write(1,'(a1)')
		write(1,'(i8,a)') NTITLE, ' !NTITLE'
		write(1,'(a9,a)') "REMARKS ","FENE polymer protein structure file, written by MD_dCSE"
			
		! Atoms
		write(1,'(a1)')
		write(1,'(i8,a)') NATOM, ' !NATOM'

		do i=1,NATOM
			write(1,'(i8,3a,i4,6a,2f10.5,i1)') i,' ',seg_name(i),' ',res_ID(i),'&
					 & ',res_name(i),' ',atom_name(i),' ',atom_type(i),charge(i),mass(i),0 
		end do

		! Bonds
		write(1,'(a1)')
		write(1,'(i8,a)') NBONDS, ' !NBONDS'
	
		write_items = 4*ceiling(NBONDS/4.0)!todo change for branched polymers
		allocate(bonds(write_items,2))
		bonds(:,:)=0

		item=1                                                   !Initialise bond item number
		do i=1,globalnp                                          !Loop through all molecules
			do j=1,i-1                                           !Avoid double counting
				if (res_ID(i).eq.res_ID(j)) then                 !If same global chainID
					
					!If j is bonded to i then add pair to items	
					do n=1,nmonomers
						group    = ceiling(real(n)/real(intbits))
						bin_expo = mod(n,intbits)-1
						if (bin_expo.eq.-1) bin_expo = intbits - 1

                        if(btest(glob_bf(i,group),bin_expo) .and. glob_sc(j) .eq. n) then
								bonds(item,1) = i
								bonds(item,2) = j
								item=item+1
						end if

					end do
	
				end if
			end do
		end do

		!Write all bonds to .psf file
		do i=1,write_items,4
			write(1,'(8i8)')bonds(i,  1), bonds(i,  2),&
							bonds(i+1,1), bonds(i+1,2),&
							bonds(i+2,1), bonds(i+2,2),&
							bonds(i+3,1), bonds(i+3,2)
		end do

		! Angles
		! Dihedrals
		! Impropers
		! Donors
		! Acceptors
		! NNB
		! NGRP

		close(1, status='keep')	

		deallocate(bonds)

	end if
	
	deallocate(seg_name)								! Determine segment names for each atom
	deallocate(res_ID)									! Determine molecule ID for each atom
	deallocate(res_name)								! Determine name for each molecule
	deallocate(atom_name)								! Determine name for each atom
	deallocate(atom_type)								! Determine type for each atom
	deallocate(charge)									! Determine charge for each atom
	deallocate(mass)									! Determine mass of each atom
	if(allocated(bonds)) deallocate(bonds)
	deallocate(glob_sc)
	deallocate(glob_bf)

end subroutine build_psf


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
		write(1000,*) 3.4*r(i,:) !Convert to Angstroms
		write(1000,*) (3.4/2.16) *v(i,:)	!Convert to Angstroms/ps
		write(1000,*) a(i,:)
	enddo

	close(1000,status='keep')

end subroutine write_DL_POLY_config


