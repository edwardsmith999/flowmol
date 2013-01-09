!=============================================================================
!				Ensemble average - the double precision edition
! Fortran code to read simulation output and add to cumulative ensemble 
! average file
!-----------------------------------------------------------------------------

program ensemble_average_dp
implicit none

	integer										:: i,length,ensemble_count,int_or_dp
	integer										:: nrecords,datasize,file_size,file_size_ensemble,sizeofdata
	integer,dimension(3)						:: gnbins
	integer,dimension(:),allocatable			:: ibuf1,ibuf2
	double precision,dimension(:),allocatable	:: buf1,buf2
	logical										:: ex,finalise
	character(30)								:: current_filename,ensemble_filename,averagetype
	character(30),dimension(8)					:: typearray

	integer,parameter							:: isint=0,isdp=1

	ex = .false.
	gnbins = (/ 43, 8, 8 /)

	!Take in command line arguments
	call setup_command_arguments(averagetype,finalise)
	current_filename = './'//averagetype
	ensemble_filename = './ensemble/'//averagetype

	!Determine size of datatype for averagetype specified on commandline
	select case(averagetype)
	!Integers
	case('mbins')
		datasize =    gnbins(1)*gnbins(2)*gnbins(3)
		int_or_dp = isint
	case('msnap')
		datasize =    gnbins(1)*gnbins(2)*gnbins(3)
		int_or_dp = isint
	case('mflux')
		datasize = 6* gnbins(1)*gnbins(2)*gnbins(3)
		int_or_dp = isint
	!Double precision
	case('Tbins')
		datasize =    gnbins(1)*gnbins(2)*gnbins(3)
		int_or_dp = isdp
	case('vbins')
		datasize = 3 *gnbins(1)*gnbins(2)*gnbins(3)
		int_or_dp = isdp
	case('vsnap')
		datasize = 3 *gnbins(1)*gnbins(2)*gnbins(3)
		int_or_dp = isdp
	case('vflux')
		datasize = 18*gnbins(1)*gnbins(2)*gnbins(3)
		int_or_dp = isdp
	case('psurface')
		datasize = 18*gnbins(1)*gnbins(2)*gnbins(3)
		int_or_dp = isdp
	!Clean up
	case('clean_all')
		typearray = (/'mbins','msnap','mflux','Tbins','vbins','vsnap','vflux','psurface'/)
		do i=1,8
			averagetype = typearray(i)
			open (unit=1, file='./ensemble/'//averagetype)
			close (unit=1,status='delete')
			open (unit=1, file='./ensemble/ensemble_count'//averagetype)
			close (unit=1,status='delete')
		enddo
		stop "Deleted all ensemble averaged files"
	end select

	!Allocate buffers for integer or double precision
	select case(int_or_dp)
	case(isint)
		allocate(ibuf1(datasize))
		allocate(ibuf2(datasize))
		inquire(iolength=length) ibuf1
		sizeofdata = 4
	case(isdp)
		allocate(buf1(datasize))
		allocate(buf2(datasize))
		inquire(iolength=length) buf1
		sizeofdata = 8
	end select
	
	!Establish filesize and number of records
	inquire(file=current_filename,size=file_size)
	nrecords = file_size/(datasize*sizeofdata)

	!print*, averagetype,current_filename ,ensemble_filename,datasize,file_size,sizeofdata,nrecords,length,size(buf1),size(ibuf1)

	!Check file exists 
	inquire(file=ensemble_filename,exist=ex)
	if (ex .eq. .false.) then

		!If file doesn't exist, create a file full of zeros
		open (unit=2, file=ensemble_filename,access='direct',recl=length,action='write',status='new')
		select case(int_or_dp)
		case(isint)
			ibuf2 = 0
			do i = 1 , nrecords
				write(2, rec=i) ibuf2(:)
			enddo
		case(isdp)
			buf2 = 0.d0
			do i = 1 , nrecords
				write(2, rec=i) buf2(:)
			enddo
		end select
		close(2,status='keep')

		!Reset count of number of averages
		open (unit=3, file='./ensemble/ensemble_count'//averagetype)
		ensemble_count = 1
		print'(3a,i8)', 'number of ensembles averaged over for ',averagetype,'=', ensemble_count
		write(3,*) ensemble_count
		close (unit=3,status='keep')
	else
		!Keep count of number of averages
		open (unit=3, file='./ensemble/ensemble_count'//averagetype)
		read(3,*) ensemble_count
		ensemble_count = ensemble_count + 1
		print'(3a,i8)', 'number of ensembles averaged over for ',averagetype,' = ', ensemble_count
		backspace(3)	!Overwrite previous value
		write(3,*) ensemble_count
		close (unit=3,status='keep')
	endif

	!Check filesizes agree for ensemble average an latest contribution
	inquire(file=ensemble_filename,size=file_size_ensemble)
	if (file_size_ensemble .ne. file_size) stop "Ensemble and averaged filesize do not agree"

	!Open files and average
	open (unit=1, file= current_filename,access='direct',recl=length,action='read')
	open (unit=2, file=ensemble_filename,access='direct',recl=length,action='readwrite')

	!Ensemble average for integer or double precision
	select case(int_or_dp)
	case(isint)
		do i = 1 , nrecords
			read(1, rec=i) ibuf1(:)
			read(2, rec=i) ibuf2(:)
			ibuf2(:) = ibuf2(:) + ibuf1(:)
			if (finalise) ibuf2 = nint(dble(ibuf2)/ensemble_count) !Take average if final record
			write(2, rec=i) ibuf2(:)
		enddo
	case(isdp)
		do i = 1 , nrecords
			read(1, rec=i) buf1(:)
			read(2, rec=i) buf2(:)
			buf2(:) = buf2(:) + buf1(:)
			if (finalise) buf2 = buf2/ensemble_count !Take average if final record
			write(2, rec=i) buf2(:)
		enddo
	end select

	close(1,status='keep')
	close(2,status='keep')

	select case(int_or_dp)
	case(isint)
		deallocate(ibuf1)
		deallocate(ibuf2)
	case(isdp)
		deallocate(buf1)
		deallocate(buf2)
	end select


end program ensemble_average_dp


!=============================================================================
!							command_arguments
! Checks for command-line arguments passed to the program and assigns the 
! relevant values.
!-----------------------------------------------------------------------------
subroutine setup_command_arguments(arg,finalise)
	implicit none
	
	integer						:: i,argcount
	logical,intent(out)			:: finalise
	character(*),intent(out) 	:: arg
	character(200) 				:: argtemp

	finalise = .false.
	argcount = command_argument_count()
	select case(argcount)
	case (1) 
		call get_command_argument(1,argtemp)			
		arg = trim(argtemp)
	case (2) 
		call get_command_argument(1,argtemp)			
		arg = trim(argtemp)
		finalise = .true.
	case (3:) 
		stop "Too many arguments specified - arg or arg + finalise required"
	case default
		stop "please specify input arguments averagetype filename"
	end select

end subroutine setup_command_arguments
