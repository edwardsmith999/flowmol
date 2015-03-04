! A program, run in serial which takes the current restartfile and
! writes it in a form which can be read by VMD

program finalstate2vmdout
	implicit none

	integer :: i,j,k,m,n,pos,endpos
	integer :: np,globalnp,Nsteps,tplot,potential_flag,nmonomers
	integer :: rtrue_flag, solvent_flag
	integer :: tag
	integer							:: NSET			!--Number of frames
	integer							:: ISTRT		!--Starting frame
	integer							:: NSAVC		!--Number of frames per coordinate save
	integer							:: NATOMNFREAT		!--Number of fixed atoms
	integer							:: NTITLE		!--Number of 80-character strings in title (set as 2)
	integer							:: NATOM		!--Number of atoms
	integer, parameter 				:: LongInt = selected_int_kind (18)
	integer(kind=LongInt)			:: bufsize, starti, endi
	integer, dimension (5)			:: FIVEZ		!--According to documentation, five zeros, but first is actually NSET
	integer, dimension (9)			:: NINEZ		!--Nine zeros
	character(len=4)				:: HDR			!--Header string, value 'CORD'
	character(len=80), dimension(2)	:: TITLE	        !--Title(s)
	double precision				:: DELTA		!--Timestep between frames
	real,allocatable,dimension(:) :: Xbuf, Ybuf, Zbuf, buf_reduced

	character(40)					:: filename_out
	integer,parameter :: nd = 3
    integer(selected_int_kind(18)) :: header_pos,end_pos ! 8 byte integer for header address
	integer,dimension(2) :: seed
	integer,dimension(3) :: initialnunits, pa,periodic
	integer,dimension(:),allocatable :: procnp, proctethernp
	integer,dimension(:,:,:),allocatable :: npcount, tethernpcount,fileunit
	double precision :: density,rcutoff,delta_t,elapsedtime,k_c,R_0,delta_rneighbr
	double precision :: tagbuf, dpbuf, large, xmin, xmax, ymin, ymax, zmin, zmax
	double precision,dimension(3) :: domain,halfdomain,globaldomain
	double precision,dimension(3) :: dp3buf

	character(80) :: input_file,initial_microstate_file
	integer, dimension(5), parameter :: tether_tags = (/3,5,6,7,10/)

	call setup_command_arguments

	!=======================================================================
	!
	!				Read finalstate file and store positions
	!
	!=======================================================================


	!Read Header first to get globalnp
	open(2,file=initial_microstate_file,form='unformatted',action='read',access='stream',position='append')
    inquire(2,POS=end_pos) 				! go the end of file
    read(2,pos=end_pos-8) header_pos 	! header start is in the last 8 bytes
    header_pos = header_pos +1 			! for compatibility with MPI IO we store header_pos - 1 in final state 
	read(2,pos=header_pos) globalnp
	read(2) initialnunits(1)                   	!x dimension box split into number of cells
	read(2) initialnunits(2)                   	!y dimension box split into number of cells
	read(2) initialnunits(3)					!z dimension box split into number of cells

	!Discard input file values
	read(2) Nsteps  	    
	read(2) tplot
	read(2) seed(1)
	read(2) seed(2)
	read(2) periodic(1)
	read(2) periodic(2)
	read(2) periodic(3)
	read(2) potential_flag
	read(2) rtrue_flag 

	close(2,status='keep')

	!Get size of arrays required to store data
	allocate(Xbuf(globalnp))
	allocate(Ybuf(globalnp))
	allocate(Zbuf(globalnp))

	large = 100000.d0
	!xmax =  large;	xmin = -large
	xmax =  large;	xmin = -large
	ymax =  1.d0;	ymin = -1.d0
	zmax =  large;	zmin = -large
	!zmax =  large;	zmin = -large

	!---------- For all molecule positions ------------
	!Move through location of position co-ordinates
	open(2,file=initial_microstate_file, form='unformatted',action='read', access='stream',position='rewind')
	np = 0; m = 1
	do n=1,globalnp
		!Read molecular tag and convert in integer
		read(2) tagbuf; tag = nint(tagbuf)	
		!Read molecular position
		read(2) dp3buf

		if ((dp3buf(1) .lt. xmax .and. dp3buf(1) .gt. xmin) .and. & 
			(dp3buf(2) .lt. ymax .and. dp3buf(2) .gt. ymin) .and. & 
			(dp3buf(3) .lt. zmax .and. dp3buf(3) .gt. zmin)  ) then
				Xbuf(m) = real(dp3buf(1))
				Ybuf(m) = real(dp3buf(2))
				Zbuf(m) = real(dp3buf(3))
				np = np + 1
				m = m + 1
		endif

		read(2) dp3buf !Read 3dp to temp buffer and discard
		if (rtrue_flag.eq.1) then
			read(2) dp3buf  !Read 3dp to temp buffer and discard
		end if
		if (any(tag.eq.tether_tags)) then 
			read(2) dp3buf   !Read 3dp to temp buffer and discard
		end if
		if (potential_flag.eq.1) stop 'Not yet developed for polymers'
	enddo
	close(2,status='keep')	

	if (np .ne. globalnp) then
		allocate(buf_reduced(np)); buf_reduced = Xbuf(1:np)
		deallocate(Xbuf); allocate(Xbuf(np)); 
		Xbuf = buf_reduced; deallocate(buf_reduced)

		allocate(buf_reduced(np)); buf_reduced = Ybuf(1:np)
		deallocate(Ybuf); allocate(Ybuf(np)); 
		Ybuf = buf_reduced; deallocate(buf_reduced)

		allocate(buf_reduced(np)); buf_reduced = Zbuf(1:np)
		deallocate(Zbuf); allocate(Zbuf(np)); 
		Zbuf = buf_reduced; deallocate(buf_reduced)
	else

	endif

	!=======================================================================
	!
	!				Write header and positions to VMD file
	!
	!=======================================================================

	!Set header information	
	HDR			=	'CORD'				!header text
	NSET		=	1					!number of recorded frames
	ISTRT		=	0					!the starting timestep
	NSAVC		=	1					!number of timesteps between dcd frame saves
	FIVEZ(1)	=	NSET				!not sure why
	FIVEZ(2:5)	=	0					!buffer zeros
	NATOMNFREAT	=	0					!number of fixed atoms?
	DELTA		=	delta_t				!delta_t (x-plor is double, charmm is real)
	NINEZ(:)	=	0					!buffer zeros
	NTITLE		=	2					!number of 80-character strings in title
	TITLE(1)	=	'  Simulation record file '	!
	TITLE(2)    =	'   Written in serial or parallel   '	!
	NATOM		=	np			!number of particles

	!Open binary .dcd file 
	write (filename_out,'(2a)') "vmd_out.dcd"	
	open(unit=3, file=filename_out,status='replace', form="unformatted")
	!write header information	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM
	!write positions
	write(3) Xbuf
	write(3) Ybuf
	write(3) Zbuf
	close(3,status='keep')

	deallocate(Xbuf)
	deallocate(Ybuf)
	deallocate(Zbuf)

contains

!=============================================================================
! setup_command_arguments
! Checks for command-line arguments passed to the program and assigns the 
! relevant values.
!-----------------------------------------------------------------------------
subroutine setup_command_arguments
	implicit none
	
	integer				:: i,argcount
	logical 			:: restart, initial_microstate_file_exists
	character(len=200) 	:: arg,nextarg


	!Set default values in case they aren't specified by user
	initial_microstate_file_exists = .false.
	initial_microstate_file = 'final_state'

	argcount = command_argument_count()

	if (argcount.gt.0) then			!If more than 0 arguments	
		do i=1,argcount 		!Loop through all arguments...
			call get_command_argument(i,arg)	!Reading two at once
                            if ( i < argcount) then 
                                    call get_command_argument(i+1,nextarg)
                            else
                                    nextarg=''
                            endif
			
			if (trim(arg).eq.'-r' .and. nextarg(1:1).ne.'-') then
				initial_microstate_file = trim(nextarg)
				inquire(file=initial_microstate_file, exist=restart) !Check file exists
			end if			
		end do
	end if

	inquire(file=initial_microstate_file, exist=initial_microstate_file_exists)	!Check file exists

	if(.not. initial_microstate_file_exists) then
		print*, 'Restart file ', trim(initial_microstate_file), ' not found or not specified.'
		print*, 'Input should be of the form ./finalstate2vmdout.exe -r final_state '
		stop
	end if

end subroutine setup_command_arguments



!-----------------------------------------------------------------------------
! Subroutine:	locate(keyword)
! Author(s):	David Trevelyan
! Description:
!		The file opened with 'fileid' is scanned sequentially by line until the 
!		character string 'keyword' matches the beginning of what is read.
!		The file position is thus set to the line after the matched keyword.
!
!		If a keyword is not found, it is assumed that the file is
!		incomplete and the program is stopped.
!-----------------------------------------------------------------------------

subroutine locate(fileid,keyword,required,input_present)
	implicit none
	
	character*(*),intent(in)		:: keyword			! Input keyword	
	integer,intent(in)				:: fileid			! File unit number
	logical,intent(in)				:: required			! Flag to check if input is required
	logical,intent(out),optional	:: input_present	! Optional flag passed back if present in intput

	character*(100)			:: linestring		! First 100 characters in a line
	integer					:: keyword_length	! Length of input keyword
	integer					:: io				! File status flag

	!Check if required, if not then check if input check variable is include
	!and terminate if missing
	if (.not. required) then
		if(present(input_present)) then
			input_present = .true.	!Assume true until not found
		else
			stop  "ERROR IN LOCATE - If input not required, extra logical argument&
				& must be included to check if variable is specified in input file"
			
		endif
	endif

	keyword_length = len(keyword)
	rewind(fileid)	! Go to beginning of input file
	do
		read (fileid,*,iostat=io) linestring			! Read first 100 characters
		if (linestring(1:keyword_length).eq.keyword) exit	! If the first characters match keyword then exit loop

		if (io.ne.0) then	! If end of file is reached
			if (.not. required) then
				!print*, keyword, ' - Not specified in input file - default value taken'
				input_present = .false.
				return
			else
				print*, "ERROR IN LOCATE -  Required input "//trim(keyword)//"& 
					& not found in input file. Stopping simulation."
				stop
			endif
		end if

	end do	

end subroutine locate



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

end program finalstate2vmdout
