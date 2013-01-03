!======================================================================
!	        	Parallel Read Write Subroutines               

! Parallel i/o

! --- INPUT ROUTINES ---
! setup_inputs				Read input file
! setup_restart_inputs		Read restart files and input file
! setup_restart_microstate	Read Initial configuration from restart file

! --- OUTPUT ROUTINES ---
! parallel_io_final_state 	Write final configuration and simulation details for restart
! parallel_io_vmd			Output for VMD
! parallel_io_vmd_sl		Output for VMD with seperate solid/lquid regions
! parallel_io_vmd_halo		Output Halos
! CV AVERAGING
! mass_averaging			slice/CV
! cumulative_mass			slice/CV
! mass_snapshot				CV
! momentum_averaging		slice/CV
! cumulative_momentum		slice/CV
! momentum_snapshot			CV
! pressure_averaging		domain/CV
! cumulative_pressure		domain/CV
! FLUX AVERAGING
! mass_flux_averaging		CV_surface
! cumulative_mass_flux		CV_surface
! momentum_flux_averaging	plane(MOP)/CV_surface
! cumulative_momentum_flux	plane(MOP)/CV_surface
! surface_pressure			plane(MOP)/CV_surface
! cumulative_pressure		plane(MOP)/CV_surface
! 
!======================================================================

module module_parallel_io
	use mpi
	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use polymer_info_MD
	use shear_info_MD
	use calculated_properties_MD
	use messenger, only : MD_COMM
	use interfaces
	implicit none

	integer		:: restartfileid, fileid, fileidtrue !File name used for parallel i/o

	!Generic interface so write arrays can be used with both integers and reals
	interface write_arrays
		module procedure iwrite_arrays_1, rwrite_arrays_1, iwrite_arrays, rwrite_arrays 
	end interface write_arrays

	private  iwrite_arrays_1, rwrite_arrays_1, iwrite_arrays, rwrite_arrays

contains


!====================================================================
!			Array writing subroutines
!--------------------------------------------------------------------


! --- Write integer arrays ---

!	1D array wrapper
subroutine iwrite_arrays_1(temp,nresults,outfile,outstep)

	integer, intent(in)						:: nresults,outstep
	integer, dimension(:,:,:),intent(in)	:: temp
	character(*),intent(in) 				:: outfile

	integer, dimension(:,:,:,:),allocatable	:: some_array

	allocate(some_array(size(temp,1),size(temp,2),size(temp,3),1))
	some_array(:,:,:,1) = temp(:,:,:)
	call iwrite_arrays(some_array,nresults,outfile,outstep)
	deallocate(some_array)

end subroutine iwrite_arrays_1

!	Main iwrite Routine
subroutine iwrite_arrays(some_array,nresults,outfile,outstep)
	use messenger, only : icomm_grid
	implicit none

	integer, intent(in)						:: nresults,outstep
	integer, dimension(:,:,:,:),intent(in)	:: some_array!(nx,ny,nz,nresults)
	character(*),intent(in) 				:: outfile

	integer									:: n, fh
	integer 								:: MEM_FLAG = 0
	integer 								:: FILE_FLAG = 0
	integer									:: global_cnt,int_size,datatype
	integer									:: status(mpi_status_size)
	integer			   						:: filetype, memtype
	integer (kind=MPI_offset_kind)			:: offset
	integer, dimension(3)					:: gsizes, lsizes, memsizes
	integer, dimension(3)					:: global_indices, local_indices
	integer, dimension(:,:),allocatable 	:: proc_lsizes 

	integer, allocatable,dimension(:,:,:)	:: OutBuffer

	datatype = MPI_INTEGER
	call MPI_TYPE_SIZE(datatype, int_size, ierr)

	call MPI_file_open(icomm_grid, outfile, MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   					MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	!  Note:  MPI assumes here that numbering starts from zero
	!  Since numbering starts from (one), subtract (one) from every index
	!-------------------------------------------------------
	global_cnt 	= (outstep-1)*gnbins(1)*gnbins(2)*gnbins(3)*nresults
	offset		= global_cnt * int_size
	gsizes 		= gnbins
	lsizes		= nbins
	local_indices(:) = (/  0  , 0 , 0 /)
	!Calculate global_indices
	global_indices(:)= (/  0  , 0 , 0 /)
	allocate(proc_lsizes(3,nproc))
	call globalGather(lsizes,proc_lsizes,3)
	global_indices(1) = sum(proc_lsizes(1,1:iblock-1))
	global_indices(2) = sum(proc_lsizes(2,1:jblock-1))
	global_indices(3) = sum(proc_lsizes(3,1:kblock-1))
	deallocate(proc_lsizes)

	!Allocate ouput buffer
	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))
	memsizes = lsizes

	do n =1,nresults
		!Copy to outbuffer
		!print*, n,1+nhb(1),nbins(1)+nhb(1),1+nhb(2),nbins(2)+nhb(2),1+nhb(3),nbins(3)+nhb(3)
		OutBuffer =  some_array(1+nhb(1):nbins(1)+nhb(1),1+nhb(2):nbins(2)+nhb(2),1+nhb(3):nbins(3)+nhb(3),n)
		!Update array datatype and reset fileview to correct section of array
		CALL Create_commit_fileview(gsizes,lsizes,global_indices,offset,datatype,FILE_FLAG,filetype,fh)
		!Update local array datatype to ignore halo cells
		CALL Create_commit_subarray(memsizes,lsizes,local_indices,datatype,MEM_FLAG,memtype)
		!Write to file
		CALL MPI_file_write_all(fh, OutBuffer, 1, memtype, status, ierr)

		!Calculate global cnt offset
		global_cnt = global_cnt + gnbins(1)*gnbins(2)*gnbins(3)
		offset = global_cnt * int_size
	
	enddo

	deallocate(OutBuffer)
	CALL MPI_FILE_CLOSE(fh, ierr)

	!==========================================================
	!     FREE DATA TYPES
	!----------------------------------------------------------
	CALL MPI_BARRIER(icomm_grid,IERR)
	call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
	call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

end subroutine iwrite_arrays


! --- Double precision arrays ---

!	1D array wrapper
subroutine rwrite_arrays_1(temp,nresults,outfile,outstep)

	integer, intent(in)								:: nresults,outstep
	double precision, dimension(:,:,:),intent(in)	:: temp
	character(*),intent(in) 						:: outfile

	double precision, dimension(:,:,:,:),allocatable	:: some_array

	allocate(some_array(size(temp,1),size(temp,2),size(temp,3),1))
	some_array(:,:,:,1) = temp(:,:,:)
	call rwrite_arrays(some_array,nresults,outfile,outstep)
	deallocate(some_array)

end subroutine rwrite_arrays_1

!	Main rwrite Routine
subroutine rwrite_arrays(some_array,nresults,outfile,outstep)
	use messenger, only : icomm_grid
	implicit none

	integer, intent(in)								:: nresults,outstep
	double precision, dimension(:,:,:,:),intent(in)	:: some_array!(nx,ny,nz,nresults)
	character(*),intent(in) 						:: outfile

	integer								:: n, fh
	integer 							:: MEM_FLAG = 0
	integer 							:: FILE_FLAG = 0
	integer								:: global_cnt,dp_size,datatype
	integer								:: status(mpi_status_size)
	integer			   					:: filetype, memtype
	integer (kind=MPI_offset_kind)		:: offset
	integer, dimension(3)				:: gsizes, lsizes, memsizes
	integer, dimension(3)				:: global_indices, local_indices
	integer, dimension(:,:),allocatable :: proc_lsizes 
	double precision, allocatable,dimension(:,:,:) 		:: OutBuffer
	!print'(a,3i8,9f10.5)','somearray', iter,irank,size(some_array),some_array(6,6,6,:)

	datatype = MPI_DOUBLE_PRECISION
	call MPI_TYPE_SIZE(datatype, dp_size, ierr)

	call MPI_file_open(icomm_grid, outfile, MPI_MODE_WRONLY+MPI_MODE_CREATE, &
			   					MPI_INFO_NULL, fh, ierr)

	!--------- DEFINE LIMITS (FILE & LOCAL SUBARRAY) -------
	!  Note:  MPI assumes here that numbering starts from zero
	!  Since numbering starts from (one), subtract (one) from every index
	!-------------------------------------------------------
	global_cnt 	= (outstep-1)*gnbins(1)*gnbins(2)*gnbins(3)*nresults
	offset		= global_cnt * dp_size
	gsizes 		= gnbins
	lsizes		= nbins
	local_indices(:) = (/  0  , 0 , 0 /)
	!Calculate global_indices
	global_indices(:)= (/  0  , 0 , 0 /)
	allocate(proc_lsizes(3,nproc))
	call globalGather(lsizes,proc_lsizes,3)
	global_indices(1) = sum(proc_lsizes(1,1:iblock-1))
	global_indices(2) = sum(proc_lsizes(2,1:jblock-1))
	global_indices(3) = sum(proc_lsizes(3,1:kblock-1))
	deallocate(proc_lsizes)

	!Allocate ouput buffer
	allocate(OutBuffer(lsizes(1),lsizes(2),lsizes(3)))
	memsizes = lsizes

	do n =1,nresults
		!Copy to outbuffer
		OutBuffer =  some_array(1+nhb(1):nbins(1)+nhb(1),1+nhb(2):nbins(2)+nhb(2),1+nhb(3):nbins(3)+nhb(3),n)
		!print*,irank, n, OutBuffer(5,5,5),global_cnt,offset
		!Update array datatype and reset fileview to correct section of array
		CALL Create_commit_fileview(gsizes,lsizes,global_indices,offset,datatype,FILE_FLAG,filetype,fh)
		!Update local array datatype to ignore halo cells
		CALL Create_commit_subarray(memsizes,lsizes,local_indices,datatype,MEM_FLAG,memtype)
		!Write to file
		CALL MPI_file_write_all(fh, OutBuffer, 1, memtype, status, ierr)

		!Calculate global cnt offset
		global_cnt = global_cnt + gnbins(1)*gnbins(2)*gnbins(3)
		offset = global_cnt * dp_size
	
	enddo

	deallocate(OutBuffer)
	CALL MPI_FILE_CLOSE(fh, ierr)

	!==========================================================
	!     FREE DATA TYPES
	!----------------------------------------------------------
	CALL MPI_BARRIER(icomm_grid,IERR)
	call MPI_TYPE_FREE(memtype,ierr) ; MEM_FLAG = 0
	call MPI_TYPE_FREE(filetype,ierr); FILE_FLAG = 0

end subroutine rwrite_arrays

subroutine Create_commit_fileview(gsizes,lsizes,global_indices,offset,datatype,FILE_FLAG,filetype,fh)
	implicit none

	integer,intent(inout)						:: filetype
	integer,intent(in)							:: datatype,fh
	integer (kind=MPI_offset_kind),intent(in)	:: offset
	integer, dimension(3),intent(in)			:: gsizes, lsizes, global_indices
	integer,intent(inout)						:: FILE_FLAG

	if (FILE_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(filetype,ierr)
		FILE_FLAG = 0
	end if
	CALL MPI_TYPE_CREATE_SUBARRAY(nd, gsizes, lsizes, global_indices, &
									MPI_ORDER_FORTRAN, datatype,  filetype, ierr)
	CALL MPI_TYPE_COMMIT(filetype, ierr)
	CALL MPI_FILE_SET_VIEW(fh, offset, datatype, filetype, &
			       			'native', MPI_INFO_NULL, ierr)
	FILE_FLAG = 1
	
end subroutine


subroutine Create_commit_subarray(memsizes,lsizes,local_indices,datatype,MEM_FLAG,memtype)
	implicit none

	integer,intent(inout)			:: memtype
	integer,intent(in)				:: datatype
	integer, dimension(3),intent(in):: memsizes, lsizes, local_indices
	integer,intent(inout)			:: MEM_FLAG


	if (MEM_FLAG.eq.1) then
		CALL MPI_TYPE_FREE(memtype,ierr)
		MEM_FLAG = 0
	end if
	CALL MPI_TYPE_CREATE_SUBARRAY(nd, memsizes, lsizes, local_indices, &
			MPI_ORDER_FORTRAN, datatype,  memtype, ierr)
	CALL MPI_TYPE_COMMIT(memtype, ierr)
	MEM_FLAG = 1
	
end subroutine

end module module_parallel_io

!======================================================================
!	        		INPUTS			              =
!======================================================================


!=============================================================================
! setup_command_arguments
! Checks for command-line arguments passed to the program and assigns the 
! relevant values.
!-----------------------------------------------------------------------------
subroutine setup_command_arguments
use module_parallel_io
implicit none
	
	integer	:: i,argcount
	logical :: restart_file_exists, input_file_exists
	character(len=200) :: arg,nextarg
	character(len=270) :: errorstr

	!Set default values in case they aren't specified by user
	restart = .false.			
	input_file_exists = .false.
	restart_file_exists = .false.
	input_file = trim(prefix_dir)//'MD.in'
	initial_microstate_file = trim(prefix_dir)//'final_state'

	argcount = command_argument_count()

	if (argcount.gt.0) then			!If more than 0 arguments	
				
		do i=1,argcount			!Loop through all arguments...
				
			call get_command_argument(i,arg)	!Reading two at once
			if ( i < argcount) then 
				call get_command_argument(i+1,nextarg)
			else
				nextarg=''
			endif
				
			if (trim(arg).eq.'-r' .and. nextarg(1:1).ne.'-') then
				initial_microstate_file = trim(nextarg)
				inquire(file=initial_microstate_file, exist=restart) !Check file exists
				if (restart .eqv. .false.) then
					write(errorstr,'(3a)') "Restart file not found, please check file ", & 
										trim(initial_microstate_file)," exists or remove -r flag"
					call error_abort(errorstr)
				endif
			end if

			if (trim(arg).eq.'-i' .and. nextarg(1:1).ne.'-') then
				input_file = trim(nextarg)
			end if

		end do

	end if

	inquire(file=input_file, exist=input_file_exists)	!Check file exists
	if (input_file.eq.'./MD.in'.and..not. input_file_exists) input_file = './default.in'
	inquire(file=input_file, exist=input_file_exists)
	
	if(.not. input_file_exists) then
		print*, 'Input file ', trim(input_file), ' not found. Stopping simulation.'
		call error_abort
	end if

end subroutine setup_command_arguments

!-----------------------------------------------------------------------------
! Subroutine:	setup_inputs
! Author(s):	David Trevelyan & Edward Smith
! Description:
!		The input file MD.in contains capitalised keywords followed by
!		numerical values. The "locate" subroutine rewinds to the beginning
!		of MD.in and scans each line until the keyword is matched.
!
!		Consequently, the file position is set for the next statement to read
!		the line underneath the previously "located" keyword. 
!-----------------------------------------------------------------------------

subroutine setup_inputs
	use module_parallel_io
	use librarymod, only : locate
	implicit none
	
	integer 			:: n, tvalue(8)

	call random_seed(size=n)
	allocate(seed(n))

	!Read input file
	call setup_read_input

	rcutoff2= rcutoff**2		!Useful definition to save computational time
	initialstep = 0				!Set initial step to one to start
	
	if (seed(1)==seed(2)) then
		! Randomisations 
		call random_seed(get=seed(1:n))
		call date_and_time(values=tvalue)
		seed=IEOR(tvalue(8)+irank,seed)
	else 
		!Assign different random number seed to each processor
		seed =  irank
	endif

	!Assign seed to random number generator
	call random_seed(put=seed(1:n))

	elapsedtime = delta_t*Nsteps !Set elapsed time to end of simualtion

end subroutine setup_inputs

!-------------------------------------------------------------------------------------
!                    Restart Simulation inputs
! Set up inputs on every processor, based on the final state of a previous simulation
!-------------------------------------------------------------------------------------

subroutine setup_restart_inputs
	use module_parallel_io
	use librarymod, only : locate
	implicit none

	integer							:: n
	integer							:: prev_nproc
	integer 						:: extrasteps
	integer 						:: checkint
    integer(MPI_OFFSET_KIND)        :: ofs, header_ofs
    integer(selected_int_kind(18))  :: header_pos
	double precision 				:: checkdp
	character(400)					:: error_message

	!Allocate random number seed
	call random_seed(size=n)
	allocate(seed(n))

	!Read input file
	call setup_read_input
	extrasteps = Nsteps

	!=====================================================================================================!
	!========================   R E A D    R E S T A R T    H E A D E R   ================================!
	!Check if values from input file are different and alert user - all processors have
	!read the same file so only need to check on one processor

	!Open on a single process and broadcast if different
	if (irank .eq. iroot) then
            
	    call MPI_File_open(MPI_COMM_SELF, initial_microstate_file, & 
		MPI_MODE_RDONLY , MPI_INFO_NULL, restartfileid, ierr)

	    ! read the size of offset, just in case we hit a system that still uses 32 bits addresses
	    ! read in a 8 bit integer
	    ofs = -8
	    call mpi_file_seek(restartfileid,ofs,mpi_seek_end,ierr)
	    call mpi_file_read(restartfileid,header_pos	   ,1,mpi_integer8,MPI_STATUS_IGNORE,ierr)
	    header_ofs = header_pos
	    call MPI_FILE_SEEK(restartfileid,header_ofs,MPI_SEEK_SET,ierr)

	    call MPI_File_read(restartfileid,globalnp        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,checkint   	 ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		if (checkint .ne. initialnunits(1)) then
			print*, 'Discrepancy between x domain size', &
					'in input & restart file - restart file will be used', checkint, initialnunits(1)
			initialnunits(1) = checkint
		endif
	    call MPI_File_read(restartfileid,checkint   	 ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		if (checkint .ne. initialnunits(2)) then
			print*, 'Discrepancy between y domain size', &
					'in input & restart file - restart file will be used', checkint, initialnunits(2)
			initialnunits(2) = checkint
		endif
	    call MPI_File_read(restartfileid,checkint   	 ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		if (checkint .ne. initialnunits(3)) then
			print*, 'Discrepancy between z domain size', &
					'in input & restart file - restart file will be used', checkint, initialnunits(3)
			initialnunits(3) = checkint
		endif
	    call MPI_File_read(restartfileid,Nsteps          ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !tplot
	    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !seed
	    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !seed
	    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !periodic
	    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !periodic
	    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !periodic
	    call MPI_File_read(restartfileid,checkint  		 ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !potential_flag
		if (checkint .ne. potential_flag) then
		    print*, 'Discrepancy between potential_flag', &
					  'in input & restart file - restart file will be used'
			potential_flag = checkint
		endif
	    call MPI_File_read(restartfileid,prev_rtrue_flag ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !rtrue_flag
		if (prev_rtrue_flag .ne. rtrue_flag) then
		    print*, 'Discrepancy between rtrue_flag', &
					'in input & restart file - current file will be used,', &
			        ' but rtrue will still be read from restart file.'  
			print*, 'prev_rtrue_flag:', prev_rtrue_flag
			print*, 'rtrue_flag:', rtrue_flag
		endif
		call MPI_File_read(restartfileid,checkint  		 ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !solvent_flag
		if (checkint .ne. solvent_flag) then
		    print*, 'Discrepancy between solvent_flag', &
					  'in input & restart file - restart file will be used'
			solvent_flag = checkint
		endif

		call MPI_File_read(restartfileid,checkint       ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr) !nmonomers
		if (checkint.ne.nmonomers) then
			print*, 'Discrepancy between nmonomers', &
				  'in input & restart file - restart file will be used'
			nmonomers = checkint
		endif

		procnp = 0; proc_reorder = 0; 	prev_nproc = 1
		!Small debugging run (nproc<27) - if proc mismatch use serial reordering (all read everything and discard)
		if(npx .le. 3 .and. npy .le. 3 .and. npz .le. 3) then
			error_message = 'Small debug run (less than 3 x 3 x 3 processors). &
							&Molecules will be assigned to correct processors - all read everything and discard'
		    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
			if (checkint .ne. npx) 	proc_reorder = 1; prev_nproc = prev_nproc*checkint
		    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
			if (checkint .ne. npy) 	proc_reorder = 1; prev_nproc = prev_nproc*checkint
		    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
			if (checkint .ne. npz) 	proc_reorder = 1; prev_nproc = prev_nproc*checkint
			if (proc_reorder.eq.0) then
				do n=1,prev_nproc			!Loop through all processors and store for restart
					call MPI_File_read(restartfileid,procnp(n),1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
				enddo
				do n=1,prev_nproc			!Loop through all processors and store for restart
					call MPI_File_read(restartfileid,proctethernp(n),1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
				enddo
			else
				do n=1,prev_nproc			!Loop through all processors and discard
					call MPI_File_read(restartfileid,checkint,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
				enddo
				do n=1,prev_nproc			!Loop through all processors and discard
					call MPI_File_read(restartfileid,checkint,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
				enddo
			endif
		!Large run - if proc mismatch, stop code and suggest reordering in serial
		else
			error_message = 'Number of processors in input file does not match the restart file.            &
							&Options:                                                                        & 
							&1) Set processors in input file to zeros to use restart proc topology.         & 
							&2) Reorder restart file in serial for current proc topology                    & 
							&3) Run in serial or with less than 3x3x3 processors                        '
		    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
			if (checkint .ne. npx) call error_abort(error_message)
		    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
			if (checkint .ne. npy) call error_abort(error_message)
		    call MPI_File_read(restartfileid,checkint        ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
			if (checkint .ne. npz) call error_abort(error_message)
			prev_nproc = npx*npy*npz	!If parallel run, number of molecules per processor written
			do n=1,prev_nproc			!Loop through all processors and store for restart
				call MPI_File_read(restartfileid,procnp(n),1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
			enddo
			do n=1,prev_nproc			!Loop through all processors and store for restart
				call MPI_File_read(restartfileid,proctethernp(n),1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
			enddo
		endif
	    call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
		if (checkdp .ne. density) then
			print*, 'Discrepancy between system density', &
					'in input & restart file - restart file will be used'
			density = checkdp
		endif
	    call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
		if (checkdp .ne. rcutoff) then
			print*, 'Discrepancy between cut off radius', &
					'in input & restart file - restart file will be used'
			rcutoff = checkdp
		endif
	    call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)	!delta_t 
	    call MPI_File_read(restartfileid,elapsedtime     ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,simtime         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
	    call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)	!k_c 
		if (checkdp.ne.k_c) then
			print*, 'Discrepancy between k_c', &
				'in input & restart file - restart file will be used'
			k_c = checkdp
		endif
	    call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)	!R_0
		if (checkdp.ne.R_0) then
			print*, 'Discrepancy between R_0', &
					  'in input & restart file - restart file will be used'
			R_0 = checkdp
		endif
	    call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)	!eps_pp
		if (checkdp.ne.eps_pp) then
			print*, 'Discrepancy between eps_pp', &
					  'in input & restart file - restart file will be used'
			eps_pp = checkdp
		endif
	    call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)	!eps_ps
		if (checkdp.ne.eps_ps) then
			print*, 'Discrepancy between eps_ps', &
					  'in input & restart file - restart file will be used'
			eps_ps = checkdp
		endif
	    call MPI_File_read(restartfileid,checkdp         ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)	!eps_ss
		if (checkdp.ne.eps_ss) then
			print*, 'Discrepancy between eps_ss', &
					  'in input & restart file - restart file will be used'
			eps_ss = checkdp
		endif
	    call MPI_File_read(restartfileid,delta_rneighbr  ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)	!delta_rneighbr

		call MPI_File_close(restartfileid,ierr)
	
	endif

	!---------------Broadcast data read by root to all other processors-------------------------!
	! temporary np, exact value will be fixed after reading r and v arrays
	! np is needed in set_parameters_outputs that is called before microstate initialisation
	! when restart is true
	call MPI_BCAST(globalnp,          1,MPI_integer,iroot-1,MD_COMM,ierr)
	np = globalnp/nproc
	call MPI_BCAST(rcutoff,           1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
	rcutoff2 = rcutoff**2             !Useful definition to save computational time
	call MPI_BCAST(initialnunits,     3,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(Nsteps,            1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(potential_flag,    1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(prev_rtrue_flag,   1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(solvent_flag,      1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(proc_reorder, 	  1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(procnp, size(procnp),MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(proctethernp, size(proctethernp),MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(nmonomers,         1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(density,           1,MPI_double_precision,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(inputtemperature,  1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
	call MPI_BCAST(elapsedtime,       1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
	call MPI_BCAST(simtime,           1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
	call MPI_BCAST(k_c,               1,MPI_double_precision,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(R_0,               1,MPI_double_precision,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(eps_pp,            1,MPI_double_precision,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(eps_ps,            1,MPI_double_precision,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(eps_ss,            1,MPI_double_precision,iroot-1,MD_COMM,ierr)

	elapsedtime = elapsedtime + delta_t*extrasteps 	!Set elapsed time to end of simualtion
	initialstep = Nsteps         					!Set plot count to final plot of last
	Nsteps = Nsteps + extrasteps 					!Establish final iteration step based on previous

end subroutine setup_restart_inputs

!------------------------------------------------------------------------------

subroutine setup_restart_microstate
	use module_parallel_io
	implicit none

	integer                                      :: i,n,nl,procassign
	integer                                      :: pos
	integer                                      :: dp_datasize
	integer(kind=MPI_OFFSET_KIND)                :: disp, procdisp
	integer, dimension(:), allocatable           :: bufsize
	double precision                             :: tagtemp
	double precision, dimension (nd)             :: rtemp,vtemp,rtruetemp,rtethertemp
	double precision, dimension (nsdmi)          :: monomertemp
	double precision, dimension (:), allocatable :: buf !Temporary variable

	allocate(bufsize(nproc))
	bufsize = 0
	
	!Determine size of datatypes
  	call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

	!Open restart file on all processor
	call MPI_FILE_OPEN(MD_COMM, initial_microstate_file, & 
		MPI_MODE_RDONLY , MPI_INFO_NULL, restartfileid, ierr)

	select case(proc_reorder)
	case(0)

		bufsize(irank) = 2*nd*procnp(irank) + 1*procnp(irank)
		if (prev_rtrue_flag.eq.1) then
			bufsize(irank) = bufsize(irank) + nd*procnp(irank)
		end if
		bufsize(irank) = bufsize(irank) + nd*proctethernp(irank)
		if (potential_flag .eq. 1) then
			bufsize(irank) = bufsize(irank) + nsdmi*procnp(irank)
		end if
		call globalSumIntVect(bufsize,nproc)

		!Obtain displacement of each processor using procs' np from restart file
		!with 6 position (3 wrapped, 3 unwrapped) and 3 velocity components
		!for each molecule
		procdisp = 0
		do i=1,irank-1
			procdisp = procdisp + bufsize(i)*dp_datasize
		end do
		!Obtain location to write in file
		disp =  procdisp

		allocate(buf(bufsize(irank)))
		!Set each processor to that location and write particlewise
		call MPI_FILE_SET_VIEW(restartfileid, disp, MPI_double_precision, & 
							   MPI_double_precision, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_READ_ALL(restartfileid, buf, bufsize(irank), MPI_double_precision, & 
							   MPI_STATUS_IGNORE, ierr) !Read position from file
	
		pos = 1
		do nl=1,procnp(irank)

			tag(nl) = nint(buf(pos))
			pos = pos + 1
			!Correct to local coordinates
			r(1,nl) = buf(pos)  -domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
			r(2,nl) = buf(pos+1)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
			r(3,nl) = buf(pos+2)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)
			pos = pos + 3
			!Read velocities
			v(1,nl) = buf(pos)
			v(2,nl) = buf(pos+1)
			v(3,nl) = buf(pos+2)
			pos = pos + 3
			if (prev_rtrue_flag.eq.1) then
				!Read true positions
				rtrue(1,nl) = buf(pos)
				rtrue(2,nl) = buf(pos+1)
				rtrue(3,nl) = buf(pos+2)
				pos = pos + 3
			end if
			if (any(tag(nl).eq.tether_tags)) then
				!Read tether position, corrected to local coords
				rtether(1,nl) = buf(pos)  -domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
				rtether(2,nl) = buf(pos+1)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
				rtether(3,nl) = buf(pos+2)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)
				pos = pos + 3
			end if
			if (potential_flag.eq.1) then
				!Read monomer data
				monomer(nl)%chainID        = nint(buf(pos))
				monomer(nl)%subchainID     = nint(buf(pos+1))
				monomer(nl)%funcy          = nint(buf(pos+2))
				monomer(nl)%glob_no        = nint(buf(pos+3))
				monomer(nl)%bin_bflag(1:4) = nint(buf(pos+4:pos+7))
				pos = pos + 8
			end if

		enddo

		np = procnp(irank)
	
		deallocate(buf)

	case(1)	!Reorder flag triggered

		nl = 0		!Reset local molecules count nl

		!---------- For all molecule positions ------------
		!Move through location of position co-ordinates
		do n=1,globalnp

			!---------------  READ ONE MOL -------------------------!
			!Read tag
			call MPI_FILE_READ_ALL(restartfileid, tagtemp, 1, MPI_DOUBLE_PRECISION, &
								   MPI_STATUS_IGNORE, ierr)
			!Read position	
			call MPI_FILE_READ_ALL(restartfileid, rtemp, 3, MPI_DOUBLE_PRECISION, &
								   MPI_STATUS_IGNORE, ierr)
			!Read velocity
			call MPI_FILE_READ_ALL(restartfileid, vtemp, 3, MPI_DOUBLE_PRECISION, &
								   MPI_STATUS_IGNORE, ierr)
			if (prev_rtrue_flag .eq. 1) then
				call MPI_FILE_READ_ALL(restartfileid, rtruetemp, 3, MPI_DOUBLE_PRECISION, &
									   MPI_STATUS_IGNORE, ierr)
			end if
			if (any(nint(tagtemp).eq.tether_tags)) then
				call MPI_FILE_READ_ALL(restartfileid, rtethertemp, 3, MPI_DOUBLE_PRECISION, &
									   MPI_STATUS_IGNORE, ierr)
			end if
			if (potential_flag.eq.1) then
				call MPI_FILE_READ_ALL(restartfileid, monomertemp, nsdmi, MPI_DOUBLE_PRECISION, &
									   MPI_STATUS_IGNORE, ierr)
			end if
			!------------ END READ ONE MOL -------------------------!

			!Use integer division to determine which processor to assign molecule to
			procassign = ceiling((rtemp(1)+globaldomain(1)/2.d0)/domain(1))
			if (procassign .ne. iblock) cycle
			procassign = ceiling((rtemp(2)+globaldomain(2)/2.d0)/domain(2))
			if (procassign .ne. jblock) cycle
			procassign = ceiling((rtemp(3)+globaldomain(3)/2.d0)/domain(3))
			if (procassign .ne. kblock) cycle

			!If molecules is in the domain then add to processor's total
			nl = nl + 1 !Local molecule count

			tag(nl) = nint(tagtemp)
			!Correct to local coordinates
			r(1,nl) = rtemp(1)-domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
			r(2,nl) = rtemp(2)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
			r(3,nl) = rtemp(3)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)

			v(:,nl) = vtemp(:)

			if (prev_rtrue_flag .eq. 1) then
				!Read true unwrapped positions
				rtrue(1,nl) = rtruetemp(1)
				rtrue(2,nl) = rtruetemp(2)
				rtrue(3,nl) = rtruetemp(3)
			end if	

			if (any(tag(nl).eq.tether_tags)) then
				!Read tethered positions
				rtether(1,nl) = rtethertemp(1)-domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
				rtether(2,nl) = rtethertemp(2)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
				rtether(3,nl) = rtethertemp(3)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)
			end if
			if (potential_flag.eq.1) then
				monomer(nl)%chainID        = nint(monomertemp(1))
				monomer(nl)%subchainID     = nint(monomertemp(2))
				monomer(nl)%funcy          = nint(monomertemp(3))
				monomer(nl)%glob_no        = nint(monomertemp(4))
				monomer(nl)%bin_bflag(1:4) = nint(monomertemp(5:8))
			end if
			if (mod(n,1000) .eq. 0) print'(a,f10.2)', & 
				' Redistributing molecules to input processor topology - % complete =', (100.d0*n/globalnp)
		enddo
	
		np = nl	!Correct local number of particles on processor

	case default

		call error_abort('processor re-ordering flag incorrect in restart microstate')

	end select

	! Determine number of chains by global maximum of chainID
	if (potential_flag .eq. 1) then
		nchains = maxval(monomer(:)%chainID)
		call globalMaxInt(nchains)
	end if

	!Close file used to load initial state and remove if called "final_state" 
	!to prevent confusion with final state of current run
	call MPI_FILE_CLOSE(restartfileid, ierr)
	if (initial_microstate_file .eq. './results/final_state') then
		call MPI_FILE_DELETE(initial_microstate_file, MPI_INFO_NULL, ierr)
	endif

	if (irank.eq.iroot) print*, 'Molecular tags have been read from restart file.'
	do n = 1,np
		call read_tag(n)		!Read tag and assign properties
	enddo

	!call setup_initialise_velocities_TG_parallel
	deallocate(bufsize)

end subroutine setup_restart_microstate

!======================================================================
!=	        							OUTPUTS			              =
!======================================================================

!------------------------------------------------------------------------
!Write positions and velocities of molecules to a file for restart

subroutine parallel_io_final_state
	use module_parallel_io
	use polymer_info_MD
	implicit none

	integer				   							:: n, i
	integer                                         :: pos
	integer 			   							:: dp_datasize
	integer, dimension(:), allocatable              :: bufsize
	integer(kind=MPI_OFFSET_KIND)      				:: disp, procdisp, filesize
    integer(kind=selected_int_kind(18))     		:: header_pos
	double precision, dimension(:,:), allocatable 	:: rglobal,rtetherglobal
	double precision, dimension(:)  , allocatable   :: buf

	allocate(bufsize(nproc))
	bufsize = 0
	
	!Rebuild simulation before recording final state
	call linklist_deallocateall	   		!Deallocate all linklist components
	call sendmols			   			!Exchange particles between processors
	call sort_mols						!Improved spatial locality helps restart on different proc topology
	call assign_to_cell	  	   			!Re-build linklist every timestep
	call messenger_updateborders(1)	   	!Update borders between processors
	call assign_to_neighbourlist	   	!Setup neighbourlist

	!Build array of number of particles on neighbouring
	!process' subdomains on current proccess
	call globalGathernp

	!Determine size of datatypes
  	call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

	!Adjust r according to actual location for storage according
	!to processor topology with r = 0 at centre
	allocate(rglobal(3,np))
	rglobal(1,:) = r(1,1:np)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
	rglobal(2,:) = r(2,1:np)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
	rglobal(3,:) = r(3,1:np)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

	allocate(rtetherglobal(3,np))
	rtetherglobal(1,:) = rtether(1,1:np)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
	rtetherglobal(2,:) = rtether(2,1:np)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
	rtetherglobal(3,:) = rtether(3,1:np)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

	!Remove previous final state file
	call MPI_FILE_DELETE(trim(prefix_dir)//'results/final_state', MPI_INFO_NULL, ierr)

	!Open file on all processors
	call MPI_FILE_OPEN(	MD_COMM,trim(prefix_dir)//'results/final_state', & 
						MPI_MODE_RDWR + MPI_MODE_CREATE, & 
						MPI_INFO_NULL, restartfileid, ierr)

	! Calculate buffer size	---------------------------------------!
	! Attention: np is changed inside reorder_data%sendmols call
	! Add 2*nd for r and v, 1 for tag
	bufsize(irank) = 2*procnp(irank)*nd + 1*procnp(irank)
	! If rtrue on, add nd
	if (rtrue_flag .eq. 1) then
		bufsize(irank) = bufsize(irank) + procnp(irank)*nd
	end if
	! For any tethered, add nd for rtether and count proctethernp
	proctethernp = 0
	do n=1,np
		if (any(tag(n).eq.tether_tags)) then
			proctethernp(irank) = proctethernp(irank) + 1
			bufsize(irank) = bufsize(irank) + nd
		end if
	end do
	! If polymer sim, add space for polymer info
	if (potential_flag .eq. 1) then
		bufsize(irank) = bufsize(irank) + nsdmi*procnp(irank)
	end if
	!--------------------------------------------------------------!

	!Collect all bufsizes and proctethernps
	call globalSumIntVect(bufsize,nproc)
	call globalSumIntVect(proctethernp,nproc)

	! Allocate buffer to be written to final_state file
	allocate(buf(bufsize(irank)))
	!Populate buffer -----------------------------------------------------!
	pos = 1
	do n = 1,np
		buf(pos) = real(tag(n),kind(0.d0)); pos = pos + 1
		buf(pos:pos+2) = rglobal(:,n);      pos = pos + 3
		buf(pos:pos+2) = v(:,n);            pos = pos + 3
		if (rtrue_flag .eq. 1) then
			buf(pos:pos+2) = rtrue(:,n);    pos = pos + 3
		end if
		if (any(tag(n) .eq. tether_tags)) then
			buf(pos:pos+2) = rtetherglobal(:,n);  pos = pos + 3
		end if
		if (potential_flag .eq. 1) then
			buf(pos)     = real(monomer(n)%chainID,kind(0.d0))
			buf(pos+1)   = real(monomer(n)%subchainID,kind(0.d0))
			buf(pos+2)   = real(monomer(n)%funcy,kind(0.d0))
			buf(pos+3)   = real(monomer(n)%glob_no,kind(0.d0))
			buf(pos+4:pos+7) = real(monomer(n)%bin_bflag(1:4),kind(0.d0))
			pos = pos + 8
		end if
	end do
	!---------------------------------------------------------------------!

	! Obtain file displacement for each process		
	procdisp = 0
	do i=1,irank-1
		procdisp = procdisp + bufsize(i)*dp_datasize
	enddo
	disp = procdisp
	!Set each processor to that location
	call MPI_FILE_SET_VIEW(restartfileid, disp, MPI_double_precision, & 
						   MPI_double_precision, 'native', MPI_INFO_NULL, ierr)
	! Write buffer to own space in file
	call MPI_FILE_WRITE(restartfileid, buf,bufsize(irank), & 
							MPI_double_precision, MPI_STATUS_IGNORE, ierr)
	!Close file on all processors
	call MPI_FILE_CLOSE(restartfileid, ierr)
	!This barrier is needed in order to get the correct file size in the next write
	call MPI_BARRIER(MD_COMM, ierr)	
 

	!----------------Write header------------------------
	!Written at the end for performance and simplicity reasons 
	!(See Gropp, lusk & Thakur Using MPI-2)

	!Write the header with one processor only
	if (irank .eq. iroot) then

        call MPI_file_open(MPI_COMM_SELF,trim(prefix_dir)//'results/final_state', & 
			MPI_MODE_WRONLY, MPI_INFO_NULL, restartfileid, ierr)
                
        if (ierr .ne. 0) then 
			write(0,*) "MD parallel_io: error in MPI_File open"
        endif

        call MPI_File_get_size(restartfileid,filesize,ierr)
        disp = filesize

        call MPI_file_set_view(restartfileid,disp,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierr)

        call MPI_File_write(restartfileid,sum(procnp)   ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,initialnunits ,3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		if (iter .lt. Nsteps ) then
	        call MPI_File_write(restartfileid,iter      ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		else
	        call MPI_File_write(restartfileid,Nsteps    ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		endif
        call MPI_File_write(restartfileid,tplot         ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,seed          ,2,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,periodic      ,3,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,potential_flag,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,rtrue_flag    ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,solvent_flag  ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,nmonomers     ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		call MPI_File_write(restartfileid,npx           ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		call MPI_File_write(restartfileid,npy           ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		call MPI_File_write(restartfileid,npz           ,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		call MPI_File_write(restartfileid,procnp,size(procnp),MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		call MPI_File_write(restartfileid,proctethernp,size(proctethernp),MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
		call MPI_File_write(restartfileid,density       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,rcutoff       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,delta_t       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,elapsedtime   ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,simtime       ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,k_c           ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,R_0           ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,eps_pp        ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,eps_ps        ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,eps_ss        ,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
        call MPI_File_write(restartfileid,delta_rneighbr,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

        header_pos = filesize ! just in case offset kind is 32 bit, rather improbable these days  !!!
        call MPI_File_write(restartfileid,header_pos,1,MPI_INTEGER8,MPI_STATUS_IGNORE,ierr)
        call MPI_file_close(restartfileid, ierr)

	endif

	deallocate(rglobal)	
	deallocate(buf)

end subroutine parallel_io_final_state

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer							:: procdisp
	integer							:: i, datasize
	real,dimension(np)				:: Xbuf, Ybuf, Zbuf
	real,dimension(globalnp)        :: Xbufglob
	real,dimension(globalnp)        :: Ybufglob
	real,dimension(globalnp)        :: Zbufglob
	integer                         :: n,globmolno
	integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

	Xbufglob = 0.0           !Initialise to zero so that global array...
	Ybufglob = 0.0           !...may be found by summation in parallel
	Zbufglob = 0.0           !------------------------------------------

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

	!Determine size of real datatype
 	call MPI_type_size(MPI_real,datasize,ierr)

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre
	select case(potential_flag)
	case(0)
		Xbuf(:) = r(1,1:np)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
		Ybuf(:) = r(2,1:np)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
		Zbuf(:) = r(3,1:np)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

		procdisp = 0
		!Obtain displacement of each processor using all other procs' np
		do i = 1, irank -1
			procdisp = procdisp + procnp(i)*datasize
		enddo

		!Open file on all processors
		call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/vmd_temp.dcd',      & 
		                   MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		                   MPI_INFO_NULL, fileid,     ierr)

		!-------------Write X coordinates--------------------
		!Obtain location to write in file

		!If intervals set to zero then full simulation recorded
		if (Nvmd_intervals.eq.0) then
			disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
				+ procdisp	
		else
			!Otherwise, calculate number of previous intervals
			disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp	
		endif

		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
		                       MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
		!Write information to file
		call MPI_FILE_WRITE_ALL(fileid,     Xbuf,     np, MPI_REAL, & 
		                        MPI_STATUS_IGNORE, ierr) 

		!-------------Write Y coordinates--------------------
		!Obtain location to write in file
		!If intervals set to zero then full simulation recorded
		if (Nvmd_intervals.eq.0) then
			disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
				+ procdisp	&
				+ globalnp * datasize			  	!Y Coordinate location
		else
			!Otherwise, calculate number of previous intervals
			disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp		&
				+ globalnp * datasize			  	!Y Coordinate location
		endif

		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
		                       MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
		!Write information to file
		call MPI_FILE_WRITE_ALL(fileid,     Ybuf,     np, MPI_REAL, & 
		                        MPI_STATUS_IGNORE, ierr) 

		!-------------Write Z coordinates--------------------

		!If intervals set to zero then full simulation recorded
		if (Nvmd_intervals.eq.0) then
			disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
				+ procdisp	&
				+ 2* globalnp * datasize			  	!Z Coordinate location
		else
			!Otherwise, calculate number of previous intervals
			disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp		&
				+ 2* globalnp * datasize			  	!Z Coordinate location
		endif

		!Obtain location to write in file
		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
		                       MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
		!Write information to file
		call MPI_FILE_WRITE_ALL(fileid,     Zbuf,     np, MPI_REAL, & 
		                        MPI_STATUS_IGNORE, ierr) 

		!-------------- CLOSE -------------------------------	
		call MPI_FILE_CLOSE(fileid, ierr) 
		call MPI_FILE_CLOSE(fileidtrue, ierr) 
	
	case(1)

		!Build sparse individual "global" buffers according to global molecular ID of each monomer
		do n=1,np
			globmolno           = monomer(n)%glob_no
			Xbufglob(globmolno) = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
			Ybufglob(globmolno) = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
			Zbufglob(globmolno) = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
		end do

		call globalSumVectReal(Xbufglob,globalnp)  !Global summation to complete global buffer
		call globalSumVectReal(Ybufglob,globalnp)
		call globalSumVectReal(Zbufglob,globalnp)
	
		call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/vmd_temp.dcd', & 
		                   MPI_MODE_RDWR + MPI_MODE_CREATE,       & 
		                   MPI_INFO_NULL, fileid, ierr)

		disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize   !Current iteration

		!Write X positions---------------------------------------------
		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
		                       'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(fileid,     Xbufglob,     globalnp, MPI_REAL, &      !Write buffer
		                        MPI_STATUS_IGNORE, ierr)

		disp = disp + globalnp*datasize                                          !Update file disp

		!Write Y positions---------------------------------------------
		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
		                       'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(fileid,     Ybufglob,     globalnp, MPI_REAL, &      !Write buffer
		                        MPI_STATUS_IGNORE, ierr)

		disp = disp + globalnp*datasize

		!Write Z positions---------------------------------------------
		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
		                       'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(fileid,     Zbufglob,     globalnp, MPI_REAL, &      !Write buffer
		                        MPI_STATUS_IGNORE, ierr)
	
		call MPI_FILE_CLOSE(fileid, ierr) 
	
	case default
	end select


end subroutine parallel_io_vmd

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd_true
	use module_parallel_io
	implicit none

	integer							:: procdisp
	integer							:: i, datasize
	real,dimension(np)				:: Xbuf, Ybuf, Zbuf
	real,dimension(globalnp)        :: Xbufglob
	real,dimension(globalnp)        :: Ybufglob
	real,dimension(globalnp)        :: Zbufglob
	integer                         :: n,globmolno
	integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

	Xbufglob = 0.0           !Initialise to zero so that global array...
	Ybufglob = 0.0           !...may be found by summation in parallel
	Zbufglob = 0.0           !------------------------------------------

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

	!Determine size of real datatype
 	call MPI_type_size(MPI_real,datasize,ierr)

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre
	select case(potential_flag)
	case(0)
		Xbuf(:) = rtrue(1,1:np)
		Ybuf(:) = rtrue(2,1:np)
		Zbuf(:) = rtrue(3,1:np)

		procdisp = 0
		!Obtain displacement of each processor using all other procs' np
		do i = 1, irank -1
			procdisp = procdisp + procnp(i)*datasize
		enddo

		!Open file on all processors
		call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/vmd_temp_true.dcd', & 
		                   MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		                   MPI_INFO_NULL, fileid, ierr)

		!-------------Write X coordinates--------------------
		!Obtain location to write in file

		!If intervals set to zero then full simulation recorded
		if (Nvmd_intervals.eq.0) then
			disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
				+ procdisp	
		else
			!Otherwise, calculate number of previous intervals
			disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp	
		endif

		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
		                       MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
		!Write information to file
		call MPI_FILE_WRITE_ALL(fileid,     Xbuf,     np, MPI_REAL, & 
		                        MPI_STATUS_IGNORE, ierr) 

		!-------------Write Y coordinates--------------------
		!Obtain location to write in file
		!If intervals set to zero then full simulation recorded
		if (Nvmd_intervals.eq.0) then
			disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
				+ procdisp	&
				+ globalnp * datasize			  	!Y Coordinate location
		else
			!Otherwise, calculate number of previous intervals
			disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp		&
				+ globalnp * datasize			  	!Y Coordinate location
		endif

		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
		                       MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		!Write information to file
		call MPI_FILE_WRITE_ALL(fileid,     Ybuf,     np, MPI_REAL, & 
		                        MPI_STATUS_IGNORE, ierr) 

		!-------------Write Z coordinates--------------------

		!If intervals set to zero then full simulation recorded
		if (Nvmd_intervals.eq.0) then
			disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
				+ procdisp	&
				+ 2* globalnp * datasize			  	!Z Coordinate location
		else
			!Otherwise, calculate number of previous intervals
			disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp		&
				+ 2* globalnp * datasize			  	!Z Coordinate location
		endif

		!Obtain location to write in file
		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL,          & 
		                       MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
		!Write information to file
		call MPI_FILE_WRITE_ALL(fileid,     Zbuf,     np, MPI_REAL, & 
		                        MPI_STATUS_IGNORE, ierr) 

		!-------------- CLOSE -------------------------------	
		call MPI_FILE_CLOSE(fileid, ierr) 
	
	case(1)

		!Build sparse individual "global" buffers according to global molecular ID of each monomer
		do n=1,np
			globmolno           = monomer(n)%glob_no
			Xbufglob(globmolno) = rtrue(1,n)
			Ybufglob(globmolno) = rtrue(2,n)
			Zbufglob(globmolno) = rtrue(3,n)
		end do

		call globalSumVectReal(Xbufglob,globalnp)  !Global summation to complete global buffer
		call globalSumVectReal(Ybufglob,globalnp)
		call globalSumVectReal(Zbufglob,globalnp)
	
		call MPI_FILE_OPEN(MD_COMM,trim(prefix_dir)//'results/vmd_temp_true.dcd', & 
		                   MPI_MODE_RDWR + MPI_MODE_CREATE,       & 
		                   MPI_INFO_NULL, fileid, ierr)

		disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize   !Current iteration

		!Write X positions---------------------------------------------
		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
		                       'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(fileid,     Xbufglob,     globalnp, MPI_REAL, &      !Write buffer
		                        MPI_STATUS_IGNORE, ierr)

		disp = disp + globalnp*datasize                                          !Update file disp

		!Write Y positions---------------------------------------------
		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
		                       'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(fileid,     Ybufglob,     globalnp, MPI_REAL, &      !Write buffer
		                        MPI_STATUS_IGNORE, ierr)

		disp = disp + globalnp*datasize

		!Write Z positions---------------------------------------------
		call MPI_FILE_SET_VIEW(fileid,     disp, MPI_REAL, MPI_REAL,          &      !Find position
		                       'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(fileid,     Zbufglob,     globalnp, MPI_REAL, &      !Write buffer
		                        MPI_STATUS_IGNORE, ierr)
			
		call MPI_FILE_CLOSE(fileid, ierr) 
		
	case default
	end select

end subroutine parallel_io_vmd_true

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd_sl
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer							:: procdisp
	integer							:: i,n, datasize
	real,dimension(np)				:: Xbuf, Ybuf, Zbuf
	integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

	!Determine size of real datatype
 	call MPI_type_size(MPI_real,datasize,ierr)

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre
	do n = 1, np
		select case(tag(n))
		case(0,4) 	!Liquid Molecules
			Xbuf(n) = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
			Ybuf(n) = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
			Zbuf(n) = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
		case default	!Solid molecules
			Xbuf(n) = -halfdomain(1)
			Ybuf(n) = -halfdomain(2)
			Zbuf(n) = -halfdomain(3)
		end select
	enddo

	procdisp = 0
	!Obtain displacement of each processor using all other procs' np
	do i = 1, irank -1
		procdisp = procdisp + procnp(i)*datasize
	enddo

	!================================
	!    Write liquid Molecules	=
	!================================
	!Open file on all processors
	call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/vmd_liquid_temp.dcd', & 
		MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		MPI_INFO_NULL, fileid, ierr)

	!-------------Write X coordinates--------------------

	!If intervals set to zero then full simulation recorded
	if (Nvmd_intervals.eq.0) then
		disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
			+ procdisp	
	else
		!Otherwise, calculate number of previous intervals
		disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp	
	endif


	!print*, irank, 'x disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Xbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Y coordinates--------------------

	!If intervals set to zero then full simulation recorded
	if (Nvmd_intervals.eq.0) then
		disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
			+ procdisp	&
			+ globalnp * datasize			  	!Y Coordinate location
	else
		!Otherwise, calculate number of previous intervals
		disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp	&
				+ globalnp * datasize			  	!Y Coordinate location
	endif

	!print*, irank, 'y disp', disp
	
	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Ybuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Z coordinates--------------------

	!If intervals set to zero then full simulation recorded
	if (Nvmd_intervals.eq.0) then
		disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
				+ procdisp	&
				+ 2 * globalnp * datasize			  	!Z Coordinate location
	else
		!Otherwise, calculate number of previous intervals
		disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp	&
				+ 2 * globalnp * datasize			  	!Z Coordinate location
	endif

	!print*, irank, 'z disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Zbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!Close file on all processors
	call MPI_FILE_CLOSE(fileid, ierr) 

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre
	do n = 1, np
		select case(tag(n))
		case(0, 4) 	!Liquid Molecules
			Xbuf(n) = -halfdomain(1)
			Ybuf(n) = -halfdomain(2)
			Zbuf(n) = -halfdomain(3)
		case default	!Solid molecules
			Xbuf(n) = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
			Ybuf(n) = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
			Zbuf(n) = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
		end select
	enddo

	!================================
	!    Write Solid Molecules	=
	!================================

	!Open file on all processors
	call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/vmd_solid_temp.dcd', & 
		MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		MPI_INFO_NULL, fileid, ierr)

	!-------------Write X coordinates--------------------
		
	!If intervals set to zero then full simulation recorded
	if (Nvmd_intervals.eq.0) then
		!Obtain location to write in file
		disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
			+ procdisp	
	else
		!Otherwise, calculate number of previous intervals
		disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp	
	endif

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Xbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Y coordinates--------------------

	!If intervals set to zero then full simulation recorded
	if (Nvmd_intervals.eq.0) then
		!Obtain location to write in file
		disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
			+ procdisp	&
			+ globalnp * datasize			  	!Y Coordinate location
	else
		!Otherwise, calculate number of previous intervals
		disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp	&
				+ globalnp * datasize			  	!Y Coordinate location
	endif

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Ybuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Z coordinates--------------------

	!If intervals set to zero then full simulation recorded
	if (Nvmd_intervals.eq.0) then
		!Obtain location to write in file
		disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
				+ procdisp	&
				+ 2 * globalnp * datasize			  	!Z Coordinate location
	else
		!Otherwise, calculate number of previous intervals
		disp = vmd_count * nd * globalnp * datasize & !Current iteration
				+ procdisp	&
				+ 2 * globalnp * datasize			  	!Z Coordinate location
	endif


	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Zbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!Close file on all processors
	call MPI_FILE_CLOSE(fileid, ierr) 

end subroutine parallel_io_vmd_sl

!------------------------------------------------------------------------
!Write positions of molecules to a file using one single write instead of 3

subroutine parallel_io_vmd_optimised
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer							:: procdisp
	integer							:: i, datasize
	real,dimension(np*3)			:: buf
	integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

	!Determine size of real datatype
 	call MPI_type_size(MPI_real,datasize,ierr)

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre
	do i=1,np
		buf(3*(i-1)+1) = r(1,i)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
		buf(3*(i-1)+2) = r(2,i)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
		buf(3*(i-1)+3) = r(3,i)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
	enddo

	procdisp = 0
	!Obtain displacement of each processor using all other procs' np
	do i = 1, irank -1
		procdisp = procdisp + procnp(i)*datasize
	enddo

	!Open file on all processors
	call MPI_FILE_OPEN(MD_COMM, trim(prefix_dir)//'results/vmd_temp.dcd', & 
						MPI_MODE_RDWR + MPI_MODE_CREATE, & 
						MPI_INFO_NULL, fileid, ierr)

	!-------------Write XYZ coordinates--------------------

	!If intervals set to zero then full simulation recorded
	if (Nvmd_intervals.eq.0) then
		!Obtain location to write in file
		disp =((iter-initialstep+1)/real((tplot),kind(0.d0))-1)*nd*globalnp*datasize & !Current iteration
			+ procdisp
	else
		!Otherwise, calculate number of previous intervals
		disp = vmd_count * nd * globalnp * datasize & !Current iteration
			+ procdisp
	endif

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, buf, np*3, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!Close file on all processors
	call MPI_FILE_CLOSE(fileid, ierr) 

end subroutine parallel_io_vmd_optimised

!------------------------------------------------------------------------
!Write positions of molecules in halo to a file

subroutine parallel_io_vmd_halo
	use module_parallel_io
	implicit none

	call error_abort("Cannot print vmd halos in parallel simulation - run in serial")

end subroutine parallel_io_vmd_halo

!---------------------------------------------------------------------------------
! Write value of last output iteration

subroutine update_simulation_progress_file
	use module_parallel_io
	implicit none

	if (irank .eq. iroot) then
		open (unit=99999, file=trim(prefix_dir)//"results/simulation_progress")
		write(99999,*) iter
		close(99999,status='keep')
	endif

end subroutine update_simulation_progress_file


!---------------------------------------------------------------------------------
! Record mass in a slice through the domain

subroutine mass_slice_io(ixyz)
	use module_parallel_io
	use calculated_properties_MD
	use messenger
	use interfaces
	implicit none

	integer							:: ixyz,jxyz,kxyz
	integer                         :: ijkblock(3)
	integer							:: slicefileid, int_datasize
	integer(kind=MPI_OFFSET_KIND)   :: disp

	!Get two directions orthogonal to slice direction
	kxyz = mod(ixyz,3)+1
	jxyz = mod(ixyz+1,3)+1

	ijkblock = (/iblock,jblock,kblock/)

	!Sum over all bins using directional sub communicators and gather on {ijk}block=1
	call SubcommSum(slice_mass, nbins(ixyz), jxyz)
	call SubcommSum(slice_mass, nbins(ixyz), kxyz)

	!Only root processor in each directional subcomm writes data
	if (icoord(jxyz,irank) .eq. 1 .and. icoord(kxyz,irank) .eq. 1) then

		!Determine size of datatypes
		call MPI_type_size(MPI_Integer,int_datasize,ierr)

		!Only processors on directional subcomm write
		call MPI_FILE_OPEN(icomm_xyz(ixyz), trim(prefix_dir)//'./results/mslice', & 
				   			MPI_MODE_WRONLY+ MPI_MODE_CREATE , & 
				   			MPI_INFO_NULL, slicefileid, ierr)

		!Obtain displacement of current record
		disp =   ((iter-initialstep+1)/(tplot*Nmass_ave) - 1) 	  	&		!Current iteration
		       * gnbins(ixyz)*int_datasize 	&		!Global record size
		       + nbins(ixyz)*int_datasize*(ijkblock(ixyz)-1)		!Processor location

		call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_INTEGER, & 
	 				MPI_INTEGER, 'native', MPI_INFO_NULL, ierr)

		call MPI_FILE_WRITE_ALL(slicefileid,slice_mass,nbins(ixyz), MPI_INTEGER, & 
						MPI_STATUS_IGNORE, ierr)

		call MPI_FILE_CLOSE(slicefileid, ierr)

	endif 

end subroutine mass_slice_io

!---------------------------------------------------------------------------------
! Record mass in 3D bins throughout domain

subroutine mass_bin_io(CV_mass_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer,intent(in)	:: CV_mass_out(nbinso(1),nbinso(2),nbinso(3))
	character(4),intent(in)	:: io_type

	integer	:: CVmasscopy(nbinso(1),nbinso(2),nbinso(3))
	integer :: m,nresults
	character(30) :: filename, outfile

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) 'results/m', io_type
	outfile = trim(prefix_dir)//filename

	!Swap halo surface fluxes to get correct values for all cells
	nresults = 1

	!Copy CV_mass_out so it is not changed
	CVmasscopy = CV_mass_out
	call iswaphalos(CVmasscopy,nbinso(1),nbinso(2),nbinso(3),nresults)

	!Calculate record number timestep
	if (io_type .eq. 'snap') then
		select case(CV_conserve)
		case(0)
			m = (iter-initialstep+1)/(tplot*Nmflux_ave) + 1 !Initial snapshot taken
		case(1)
			m = (iter-initialstep+1)/(Nmflux_ave) + 1 !Initial snapshot taken
		case default
			call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')	
		end select
	else
		m = (iter-initialstep+1)/(tplot*Nmass_ave)
	endif

	!Write mass to file
	call write_arrays(CVmasscopy,nresults,outfile,m)

end subroutine mass_bin_io

!---------------------------------------------------------------------------------
! Record velocity in a slice through the domain

subroutine velocity_slice_io(ixyz)
	use module_parallel_io
	use calculated_properties_MD
	use messenger
	implicit none

	integer							:: ixyz,jxyz,kxyz,ijkblock(3)
	integer							:: slicefileid, dp_datasize
	integer,dimension(3)			:: idims
	integer(kind=MPI_OFFSET_KIND)   :: disp

	!Write mass
	call mass_slice_io(ixyz)

	!Get two directions orthogonal to slice direction
	kxyz = mod(ixyz,3)+1
	jxyz = mod(ixyz+1,3)+1
	idims(1) = npx; idims(2) = npy; idims(3) = npz

	ijkblock = (/iblock,jblock,kblock/)

	!Sum over all bins using directional sub communicators and gather on root
	call SubcommSum(slice_momentum(:,1), nbins(ixyz), jxyz)
	call SubcommSum(slice_momentum(:,1), nbins(ixyz), kxyz)
	call SubcommSum(slice_momentum(:,2), nbins(ixyz), jxyz)
	call SubcommSum(slice_momentum(:,2), nbins(ixyz), kxyz)
	call SubcommSum(slice_momentum(:,3), nbins(ixyz), jxyz)
	call SubcommSum(slice_momentum(:,3), nbins(ixyz), kxyz)

	!Only root processor in each directional subcomm writes data
	if (icoord(jxyz,irank) .eq. 1 .and. icoord(kxyz,irank) .eq. 1) then

		call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

		!Only processors on directional subcomm write
		call MPI_FILE_OPEN(icomm_xyz(ixyz), trim(prefix_dir)//'./results/vslice', & 
				   MPI_MODE_WRONLY + MPI_MODE_CREATE , & 
				   MPI_INFO_NULL, slicefileid, ierr)

		!Obtain displacement of x record
		disp =   ((iter-initialstep+1)/(tplot*Nmass_ave) - 1) 	  	&	!Current iteration
		       * nd*gnbins(ixyz)*dp_datasize 	&	!times record size
		       + nbins(ixyz)*dp_datasize*(ijkblock(ixyz)-1)		!Processor location

		call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
	 				MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,1),nbins(ixyz), & 
					MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

		!Obtain displacement of y record
		disp =   ((iter-initialstep+1)/(tplot*Nmass_ave) - 1) 	  	&	!Current iteration
		       * nd*gnbins(ixyz)*dp_datasize 	&	!Record size
		       + nbins(ixyz)*dp_datasize*(ijkblock(ixyz)-1)	&	!Processor location
		       + nbins(ixyz)*dp_datasize*idims(ixyz)		!after x data 

		call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
	 				MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,2),nbins(ixyz), & 
					MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

		!Obtain displacement of z record
		disp =   ((iter-initialstep+1)/(tplot*Nmass_ave) - 1) 	  	&	!Current iteration
		       * nd*gnbins(ixyz)*dp_datasize 	&	!Record size
		       + nbins(ixyz)*dp_datasize*(ijkblock(ixyz)-1)	&	!Processor location
		       + 2*nbins(ixyz)*dp_datasize*idims(ixyz)		!after x & y data 

		call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
	 				MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,3),nbins(ixyz), & 
					MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

		call MPI_FILE_CLOSE(slicefileid, ierr)

	endif 

end subroutine velocity_slice_io

!------------------------------------------------------------------------
!A routine with each proc writing its own bins in binary

subroutine velocity_bin_io(CV_mass_out,CV_momentum_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer					:: m,nresults
	integer					:: CV_mass_out(nbinso(1),nbinso(2),nbinso(3))
	double precision		:: CV_momentum_out(nbinso(1),nbinso(2),nbinso(3),nd)
	character(4)			:: io_type
	character(30)			:: filename,outfile

	!Write mass bins
	call mass_bin_io(CV_mass_out,io_type)

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) 'results/v', io_type
	outfile = trim(prefix_dir)//filename

	! Swap Halos
	nresults = nd
	call rswaphalos(CV_momentum_out,nbinso(1),nbinso(2),nbinso(3),nresults)

	!Setup arrays
	if (io_type .eq. 'snap') then
		!CV_momentum_out = CV_momentum_out / (tplot*Nvflux_ave)
		select case(CV_conserve)
		case(0)
			m = (iter-initialstep+1)/(tplot*Nvflux_ave) + 1 !Initial snapshot taken
		case(1)
			m = (iter-initialstep+1)/(Nvflux_ave) + 1 !Initial snapshot taken
		case default
			call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')	
		end select
	else
		!CV_momentum_out = CV_momentum_out / (tplot*Nvel_ave)
		m = (iter-initialstep+1)/(tplot*Nvel_ave)
	endif

	!Write out arrays
	call write_arrays(CV_momentum_out,nresults,outfile,m)

end subroutine velocity_bin_io

!---------------------------------------------------------------------------------
! Record temperature in a slice through the domain

subroutine temperature_slice_io(ixyz)
use module_parallel_io
use calculated_properties_MD
implicit none

	integer		:: ixyz
	
	! Remove the following lines when parallelising, they are only here
	! to satisfy the compiler in pedantic mode
	integer :: dummy
	dummy = ixyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine temperature_slice_io

!---------------------------------------------------------------------------------
! Record temperature in 3D bins throughout domain

subroutine temperature_bin_io(CV_mass_out,CV_temperature_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: m,nresults
	integer				:: CV_mass_out(nbinso(1),nbinso(2),nbinso(3))
	double precision	:: CV_temperature_out(nbinso(1),nbinso(2),nbinso(3))
	character(4)		:: io_type
	character(30)		:: filename, outfile

	!Write mass bins
	call mass_bin_io(CV_mass_out,io_type)

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) 'results/T', io_type
	outfile = trim(prefix_dir)//filename

	nresults = 1

	! Swap Halos
	call rswaphalos(CV_temperature_out,nbinso(1),nbinso(2),nbinso(3),nresults)

	if (io_type .eq. 'snap') then
		!CV_temperature_out = CV_temperature_out / (tplot*Nvflux_ave)
		stop "Error in temp_io - Temperature SNAP called and NTflux_ave not defined"
		!m = iter/(NTflux_ave) + 1 !Initial snapshot taken
	else
		!CV_temperature_out = CV_temperature_out / (tplot*Nvel_ave)
		m = (iter-initialstep+1)/(tplot*NTemp_ave)
	endif
	!Write temperature to file
	call write_arrays(CV_temperature_out,nresults,outfile,m)

end subroutine temperature_bin_io

!---------------------------------------------------------------------------------
! Record velocity in 3D bins throughout domain

subroutine energy_bin_io(CV_energy_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: m,nresults
	double precision	:: CV_energy_out(nbinso(1),nbinso(2),nbinso(3))
	character(4)		:: io_type
	character(30)		:: filename, outfile

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) 'results/e', io_type
	outfile = trim(prefix_dir)//filename

	nresults = 1
	!---------------Correct for surface fluxes on halo cells---------------
	! Swap Halos
	call rswaphalos(CV_energy_out,nbinso(1),nbinso(2),nbinso(3),nresults)

	if (io_type .eq. 'snap') then
		!CV_energy_out = CV_energy_out / (tplot*Nvflux_ave)
		select case(CV_conserve)
		case(0)
			m = (iter-initialstep+1)/(tplot*Neflux_ave) + 1 !Initial snapshot taken
		case(1)
			m = (iter-initialstep+1)/(Neflux_ave) + 1 !Initial snapshot taken
		case default
			call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')	
		end select
	else
		!CV_energy_out = CV_energy_out / (tplot*Nvel_ave)
		m = (iter-initialstep+1)/(tplot*Neflux_ave)
	endif

	!Write Energy to file
	call write_arrays(CV_energy_out,nresults,outfile,m)

end subroutine energy_bin_io

!---------------------------------------------------------------------------------

subroutine virial_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none
	integer		:: m, length

	!call globalAverage(Pxy, 9)

	!Write virial pressure to file
	m = (iter-initialstep+1)/(tplot*Nstress_ave)
	if (irank .eq. iroot) then
		inquire(iolength=length) Pxy
		open (unit=7, file=trim(prefix_dir)//'results/pvirial',form='unformatted',access='direct',recl=length)
		write(7,rec=m) Pxy
		close(7,status='keep')
	endif

end subroutine virial_stress_io

!---------------------------------------------------------------------------------

subroutine VA_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer											:: ixyz, jxyz, m, nresults
	double precision								:: binvolume
	double precision,dimension(:,:,:,:),allocatable	:: buf

	!Add kinetic and configurational to Pxybin total
	Pxybin(:,:,:,:,:) = 	vvbin(:,:,:,:,:) 		& 
					      + rfbin(  2:nbins(1)+1, 	& 
						      		2:nbins(2)+1, 	& 
						      		2:nbins(3)+1,:,:)/2.d0

	!Average over samples
	Pxybin = Pxybin / Nstress_ave
	vvbin  = vvbin  / Nstress_ave 
	rfbin  = rfbin  / Nstress_ave 

	!Calculate Virial from sum of Volume Averaged Stresses
	do ixyz = 1,3
	do jxyz = 1,3
		Pxy(ixyz,jxyz) = sum(Pxybin(:,:,:,ixyz,jxyz))
	enddo
	enddo
	Pxy = Pxy / volume

	!Write out Virial stress
	call virial_stress_io

	!VA pressure per bin
	binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
	Pxybin = Pxybin / binvolume
	vvbin  = vvbin  / binvolume
	rfbin  = rfbin  / (2.d0*binvolume)

	!Write VA pressure to file
	nresults = 9
	m = (iter-initialstep+1)/(tplot*Nstress_ave)

	select case (split_kin_config)
	case(0)
		!print'(a,3i8,9f10.5)','Pxybin', iter,irank,size(Pxybin), & 
		!											Pxybin(5,5,5,1,1),Pxybin(5,5,5,2,1),Pxybin(5,5,5,3,1), &
		!											Pxybin(5,5,5,1,2),Pxybin(5,5,5,2,2),Pxybin(5,5,5,3,2), &
		!											Pxybin(5,5,5,1,3),Pxybin(5,5,5,2,3),Pxybin(5,5,5,3,3)

		!Write sum of kinetic and configurational
		!Allocate buf with halo padding and 3x3 stresses reordered as 9 vector.
		allocate(buf(nbins(1)+2,nbins(2)+2,nbins(3)+2,nresults))
		buf = 0.d0; 
		buf(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,1:9) = &
			reshape(Pxybin,(/nbins(1),nbins(2),nbins(3),nresults/))
		call write_arrays(buf,nresults,trim(prefix_dir)//'results/pVA',m)
	case(1)
		!Kinetic
		!Allocate buf with halo padding and 3x3 stresses reordered as 9 vector.
		allocate(buf(nbins(1)+2,nbins(2)+2,nbins(3)+2,nresults))
		buf = 0.d0; 
		buf(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,1:9) = &
			reshape(vvbin,(/nbins(1),nbins(2),nbins(3),nresults/))
		call write_arrays(buf,nresults,trim(prefix_dir)//'results/pVA_k',m)
		!Configurational
		!Allocate buf with halo padding and 3x3 stresses reordered as 9 vector.
		allocate(buf(nbins(1)+2,nbins(2)+2,nbins(3)+2,nresults))
		buf = 0.d0; 
		buf(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,1:9) = &
			reshape(rfbin,(/nbins(1),nbins(2),nbins(3),nresults/))
		call write_arrays(buf,nresults,trim(prefix_dir)//'results/pVA_c',m)
	case default
		stop 'Error in VA/virial extra flag to split_kinetic_& configuartional parts'
	end select

	!deallocate(buf)

end subroutine VA_stress_io

!===================================================================================
!Integrate virial pressure to get autocorrelations (Green Kubo) viscosity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine viscosity_io
	use module_parallel_io
	use physical_constants_MD
	use calculated_properties_MD
	use librarymod, only : integrate_trap
	implicit none

	integer				:: m, length
	double precision	:: viscosity

	!call globalAverage(Pxycorrel, Nvisc_ave)

	if (irank .eq. iroot) then
		call integrate_trap(Pxycorrel,tplot*delta_t,Nstress_ave,viscosity)

		viscosity = (viscosity*volume)/(3.0*Nstress_ave*Nvisc_ave*inputtemperature)

		!Write viscosity to file
		m = (iter-initialstep+1)/(tplot*Nstress_ave*Nvisc_ave)
		inquire(iolength=length) viscosity
		open (unit=7, file=trim(prefix_dir)//'results/visc',form='unformatted',access='direct',recl=length)
		write(7,rec=m) viscosity
		close(7,status='keep')
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine viscosity_io

subroutine viscometrics_io(eta,N1,N2,P)
implicit none

	real(kind(0.d0)), intent(in) :: eta, N1, N2, P

	!Remove all these lines, they are only here to satisfy pedantic compiler	
	real(kind(0.d0)) :: dummy
	dummy = eta + N1 + N2 + P

	print*, 'Viscometrics_io not developed in parallel'

end subroutine viscometrics_io

!=================================================================================
! CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV
! CV 			Record Fluxes accross surfaces of Control Volumes				CV
! CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV
!=================================================================================

!---------------------------------------------------------------------------------
! Record mass fluxes accross surfaces of Control Volumes

subroutine mass_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: m,nresults
	character(30)		:: filename, outfile

	!Work out correct filename for i/o type
	filename='results/mflux'
	outfile = trim(prefix_dir)//filename

	!Include halo surface fluxes to get correct values for all cells
	nresults = 6
	call iswaphalos(mass_flux,nbinso(1),nbinso(2),nbinso(3),nresults)

	!Calculate record number timestep
	select case(CV_conserve)
	case(0)
		m = (iter-initialstep+1)/(tplot*Nmflux_ave)
	case(1)
		m = (iter-initialstep+1)/(Nmflux_ave)
	case default
		call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
	end select

	!Write mass to file
	call write_arrays(mass_flux,nresults,outfile,m)
	!call write_arrays(mass_flux,nbinso(1),nbinso(2),nbinso(3),nresults,trim(prefix_dir)//'results/mflux',m)

end subroutine mass_flux_io

!---------------------------------------------------------------------------------
! Record momentum fluxes accross surfaces of Control Volumes

subroutine momentum_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer											:: ixyz,m,nresults
	double precision								:: binface
	double precision,allocatable,dimension(:,:,:,:)	:: temp

	! Swap Halos
	nresults = 18
	call rswaphalos(momentum_flux,nbinso(1),nbinso(2),nbinso(3),nresults)

	do ixyz = 1,3
		binface	      = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
			     		(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
		momentum_flux(:,:,:,:,ixyz  )=momentum_flux(:,:,:,:,ixyz  )/(binface) !Bottom
		momentum_flux(:,:,:,:,ixyz+3)=momentum_flux(:,:,:,:,ixyz+3)/(binface) !Top
	enddo

	!Divide momentum flux by averaing period tau=delta_t*Nvflux_ave if CV_conserve=1
	!or Divide momentum flux by sum of the Nvflux_ave times delta_t averaging periods 
	!as sample is taken every tplot steps. The output is then a representaive momentum flux.
	momentum_flux = momentum_flux/(delta_t*Nvflux_ave)

	!Write momentum to file
	select case(CV_conserve)
	case(0)
		m = (iter-initialstep+1)/(tplot*Nvflux_ave)
	case(1)
		m = (iter-initialstep+1)/(Nvflux_ave)
	case default
		call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
	end select

	allocate(temp(size(momentum_flux,1),size(momentum_flux,2),size(momentum_flux,3),18))
	temp = reshape(momentum_flux,(/ size(momentum_flux,1),size(momentum_flux,2),size(momentum_flux,3),18 /))
	call write_arrays(temp,nresults,trim(prefix_dir)//'results/vflux',m)

end subroutine momentum_flux_io

!---------------------------------------------------------------------------------
! Record stress accross plane

subroutine MOP_stress_io(ixyz_in)
	use module_parallel_io
	use calculated_properties_MD
	use messenger
	implicit none

	integer, intent(in) :: ixyz_in

	integer							:: m, ixyz, jxyz, kxyz
	integer							:: slicefileid, dp_datasize
	integer,dimension(3)			:: idims
	integer(kind=MPI_OFFSET_KIND)   :: disp

	!This line added to satisfy pedantic compiler. Remove if you want.
	ixyz = ixyz_in

	ixyz = vflux_outflag

	!Divide by number of samples taken
	Pxy_plane = Pxy_plane/(Nstress_ave)

	!Divide by area of domain and factor of 4 for interactions
	Pxy_plane = Pxy_plane/(4*domain(1)*domain(3))

	!Write plane pressures to file
	m = (iter-initialstep+1)/(tplot*Nvflux_ave)
!************WRITTEN ROUGHLY AND NOT TESTED************************
	!Write mass
	call mass_slice_io(ixyz)

	!Get two directions orthogonal to slice direction
	kxyz = mod(ixyz,3)+1
	jxyz = mod(ixyz+1,3)+1
	idims(1) = npx; idims(2) = npy; idims(3) = npz

	!Sum over all bins using directional sub communicators and gather on root
	call SubcommSum(Pxy_plane(1,:), nplanes, jxyz)
	call SubcommSum(Pxy_plane(1,:), nplanes, kxyz)
	call SubcommSum(Pxy_plane(2,:), nplanes, jxyz)
	call SubcommSum(Pxy_plane(2,:), nplanes, kxyz)
	call SubcommSum(Pxy_plane(3,:), nplanes, jxyz)
	call SubcommSum(Pxy_plane(3,:), nplanes, kxyz)

	!Only root processor in each directional subcomm writes data
	if (icoord(jxyz,irank) .eq. 1 .and. icoord(kxyz,irank) .eq. 1) then

		call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

		!Only processors on directional subcomm write
		call MPI_FILE_OPEN(icomm_xyz(ixyz), trim(prefix_dir)//'/results/pplane', & 
				   MPI_MODE_WRONLY + MPI_MODE_CREATE , & 
				   MPI_INFO_NULL, slicefileid, ierr)

		select case(vflux_outflag)
			case(1)
				!Obtain displacement of x record
				disp =   (iter/(tplot*Nmass_ave) - 1) 	  	&	!Current iteration
				       * nd*gnplanes*dp_datasize 			&	!times record size
				       + nplanes*dp_datasize*(iblock-1)			!Processor location
			case(2)
					!Obtain displacement of x record
				disp =   (iter/(tplot*Nmass_ave) - 1) 	  	&	!Current iteration
				       * nd*gnplanes*dp_datasize 			&	!times record size
				       + nplanes*dp_datasize*(jblock-1)			!Processor location
			case(3)
				!Obtain displacement of x record
				disp =   (iter/(tplot*Nmass_ave) - 1) 	  	&	!Current iteration
				       * nd*gnplanes*dp_datasize 			&	!times record size
				       + nplanes*dp_datasize*(kblock-1)			!Processor location
			case default
				call error_abort('MOP average output called with incorrect vflux_outflag')
		end select

		call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_DOUBLE_PRECISION, & 
	 				MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_WRITE_ALL(slicefileid,Pxy_plane,nplanes, & 
					MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierr)

		call MPI_FILE_CLOSE(slicefileid, ierr)

	endif 
!************WRITTEN ROUGHLY AND NOT TESTED************************

end subroutine MOP_stress_io

!---------------------------------------------------------------------------------
! Record stress accross surfaces of Control Volumes

subroutine surface_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer												:: ixyz,m,nresults
	double precision									:: binface
	double precision,dimension(:,:,:,:),allocatable		:: temp

	! Swap Halos
	nresults = 18
	call rswaphalos(Pxyface,nbinso(1),nbinso(2),nbinso(3),nresults)

	do ixyz = 1,3
		binface		  = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
			     		(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
		Pxyface(:,:,:,:,ixyz  ) = 0.25d0 * Pxyface(:,:,:,:,ixyz  )/binface !Bottom
		Pxyface(:,:,:,:,ixyz+3) = 0.25d0 * Pxyface(:,:,:,:,ixyz+3)/binface !Top
	enddo

	!Integration of stress using trapizium rule requires multiplication by timestep
	!so delta_t cancels upon division by tau=delta_t*Nvflux_ave resulting in division by Nvflux_ave
	Pxyface = Pxyface/Nvflux_ave

	!Write surface pressures to file
	select case(CV_conserve)
	case(0)
		m = (iter-initialstep+1)/(Nvflux_ave*tplot)
	case(1)
		m = (iter-initialstep+1)/(Nvflux_ave)
	case default
		call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
	end select

	!Write surface pressures to file
	allocate(temp(size(momentum_flux,1),size(momentum_flux,2),size(momentum_flux,3),18))
	temp = reshape(Pxyface,(/ size(momentum_flux,1),size(momentum_flux,2),size(momentum_flux,3),18 /))
	call write_arrays(temp,nresults,trim(prefix_dir)//'results/psurface',m)

end subroutine surface_stress_io

!---------------------------------------------------------------------------------
! Record energy fluxes accross surfaces of Control Volumes

subroutine energy_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer					:: ixyz,m,nresults
	double precision		:: binface

	!Include halo surface fluxes to get correct values for all cells
	nresults = 6
	call rswaphalos(energy_flux,nbinso(1),nbinso(2),nbinso(3),nresults)

	do ixyz = 1,3
		binface	      = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
			     		(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
		energy_flux(:,:,:,ixyz  )=energy_flux(:,:,:,ixyz  )/binface !Bottom
		energy_flux(:,:,:,ixyz+3)=energy_flux(:,:,:,ixyz+3)/binface !Top
	enddo

	!Divide energy flux by averaing period tau=delta_t*Neflux_ave if CV_conserve=1
	!or Divide sample energy flux by equivalent averaging period delta_t*Neflux_ave
	energy_flux = energy_flux/(delta_t*Neflux_ave)

	!Write energy flux to file
	select case(CV_conserve)
	case(0)
		m = (iter-initialstep+1)/(Neflux_ave*tplot)
	case(1)
		m = (iter-initialstep+1)/(Neflux_ave)
	case default
		call error_abort('CV_conserve value used for flux averages is incorrectly defined - should be 0=off or 1=on')
	end select

	call write_arrays(energy_flux,nresults,trim(prefix_dir)//'results/eflux',m)

end subroutine energy_flux_io

!---------------------------------------------------------------------------------
! Record stress times velocity (power) accross surfaces of Control Volumes

subroutine surface_power_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer					:: ixyz,m,nresults
	double precision		:: binface

	!Include halo surface stresses to get correct values for all cells
	nresults = 6
	call rswaphalos(Pxyvface,nbinso(1),nbinso(2),nbinso(3),nresults)

	do ixyz = 1,3
		binface		  = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
			     		(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
		Pxyvface(:,:,:,ixyz  ) = 0.25d0 * Pxyvface(:,:,:,ixyz  )/binface !Bottom
		Pxyvface(:,:,:,ixyz+3) = 0.25d0 * Pxyvface(:,:,:,ixyz+3)/binface !Top
	enddo

	!Integration of stress using trapizium rule requires multiplication by timestep
	Pxyvface = Pxyvface/Neflux_ave

	!Write surface pressures * velocity to file
	select case(CV_conserve)
	case(0)
		m = (iter-initialstep+1)/(Neflux_ave*tplot)
	case(1)
		m = (iter-initialstep+1)/(Neflux_ave)
	case default
		call error_abort('CV_conserve value used forsurface power is incorrectly defined - should be 0=off or 1=on')
	end select
	call write_arrays(Pxyvface,nresults,trim(prefix_dir)//'results/esurface',m)

end subroutine surface_power_io

!---------------------------------------------------------------------------------
! Record  energy accross plane

subroutine MOP_energy_io(ixyz)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer		:: ixyz

	! Remove when parallelising, only here to satisfy pedantic compiler
	integer :: dummy
	dummy = ixyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine MOP_energy_io


!-----------------------------------------------------------------------------
! Write macroscopic properties to file
!-----------------------------------------------------------------------------
subroutine macroscopic_properties_header
use module_parallel_io
use calculated_properties_MD
implicit none

	if (irank .eq. iroot) then	
		open(unit=10,file=trim(prefix_dir)//'results/macroscopic_properties',status='replace')
		if (potential_flag.eq.0) then
			write(10,'(2a)') &
			' iter; simtime; VSum; V^2Sum; Temp;', &
			' KE; PE; TE; Pressure;'
		else if (potential_flag.eq.1) then
			write(10,'(2a)') &
			' iter; simtime; VSum; V^2Sum; Temp;', &
			' KE; PE (LJ); PE (FENE); PE (Tot); TE; Pressure; Etevtcf; R_g '
		end if
	endif
	call macroscopic_properties_record

end subroutine macroscopic_properties_header


subroutine macroscopic_properties_record
use module_parallel_io
use calculated_properties_MD
implicit none

	if (irank .eq. iroot) then
		if (potential_flag.eq.0) then	
			write(10,'(1x,i8,a,f15.4,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f10.4)'), &
			iter,';',simtime,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy,';',totenergy,';',pressure
		else if (potential_flag.eq.1) then
			write(10, '(1x,i8,a,f15.4,a,f15.4,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f19.15,a,f19.15,a,f10.4,a,f10.4,a,f10.4)') &
			iter,';',simtime,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy_LJ,';',potenergy_FENE,';',potenergy,';',totenergy,';',pressure,';',etevtcf,';',R_g
		end if
	endif

end subroutine macroscopic_properties_record

!-----------------------------------------------------------------------------
! Write end-to-end vector time correlation function
!-----------------------------------------------------------------------------
subroutine etev_io
	use module_parallel_io
	implicit none

	integer :: etev_fid

	if (irank .eq. iroot) then
	
		if (iter.eq.etevtcf_iter0) then
			call MPI_FILE_DELETE(trim(prefix_dir)//'results/etev', MPI_INFO_NULL, ierr)
			call MPI_FILE_OPEN(MPI_COMM_SELF, trim(prefix_dir)//'results/etev', &
			                   MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, etev_fid, ierr)
			call MPI_FILE_WRITE(etev_fid,etev,nchains*nd,MPI_DOUBLE_PRECISION, &
			                    MPI_STATUS_IGNORE, ierr)
		else if (iter.gt.etevtcf_iter0) then
			call MPI_FILE_OPEN(MPI_COMM_SELF, trim(prefix_dir)//'results/etev', &
			                   MPI_MODE_WRONLY + MPI_MODE_APPEND, MPI_INFO_NULL, etev_fid, ierr)
			call MPI_FILE_WRITE(etev_fid,etev,nchains*nd,MPI_DOUBLE_PRECISION, &
			                    MPI_STATUS_IGNORE, ierr)
		end if
		
		call MPI_FILE_CLOSE(etev_fid,ierr)

	endif

end subroutine etev_io

!Write radius of gyration to output file
subroutine r_gyration_io
use module_parallel_io
implicit none

	integer :: Rg_fid

	if (irank .eq. iroot) then

		if (iter .eq. r_gyration_iter0) then

			call MPI_FILE_DELETE    (trim(prefix_dir)//'results/r_gyration', MPI_INFO_NULL, ierr)
			call MPI_FILE_OPEN      (MPI_COMM_SELF, trim(prefix_dir)//'results/r_gyration',       &
			                         MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,            &
			                         Rg_fid, ierr)
			call MPI_FILE_WRITE_ALL (Rg_fid, R_g, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr) 

		else if (iter .gt. r_gyration_iter0) then

			call MPI_FILE_OPEN      (MPI_COMM_SELF, trim(prefix_dir)//'results/r_gyration',       &
			                         MPI_MODE_WRONLY + MPI_MODE_APPEND, MPI_INFO_NULL,            &
			                         Rg_fid, ierr)
			call MPI_FILE_WRITE_ALL (Rg_fid, R_g, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr) 

		end if
		
		call MPI_FILE_CLOSE(Rg_fid, ierr)

	end if
	
end subroutine r_gyration_io 

subroutine rdf3d_io
	use module_parallel_io
	implicit none

	print*, 'Need to parallelise RDF i/o'

end subroutine rdf3d_io

subroutine rdf_io
	use module_parallel_io
	implicit none

	print*, 'Need to parallelise RDF i/o'

end subroutine rdf_io

subroutine ssf_io 
	use module_parallel_io
	implicit none

	print*, 'Need to parallelise SSF i/o'

end subroutine ssf_io

