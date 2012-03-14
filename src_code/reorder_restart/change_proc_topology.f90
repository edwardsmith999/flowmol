! A program, run in serial which takes the current restart file and
! reorders it to allow restart with a different number of processors than the previous
! run. This may need updating if changes to the rest of the code have adjusted the 
! order of the header or the number of values stored...

program change_proc_topology
	implicit none

	logical										:: found_in_input
	integer										:: i,j,k,l,m,n,nl,writefreq
	integer										:: globalnp,Nsteps,tplot,potential_flag,nmonomers
	integer 									:: checkint, file_size
	integer 									:: prev_nproc,npx,npy,npz,prev_npx,prev_npy,prev_npz
	integer,parameter							:: nd = 3
    integer(selected_int_kind(18))  			:: ofs,header_ofs,header_pos,end_pos ! 8 byte integer for header address
	integer,dimension(2)						:: seed
	integer,dimension(3)						:: initialnunits, pa,periodic
	integer,dimension(:)	,allocatable		:: procnp
	integer,dimension(:,:,:),allocatable		:: npcount
	double precision 									:: density,checkdp,rcutoff,delta_t,elapsedtime,k_c,R_0,delta_rneighbr
	double precision,dimension(3)						:: domain,halfdomain,globaldomain
	double precision,dimension(6)						:: buf
	double precision,dimension(:),allocatable			:: rv
	double precision,dimension(:,:,:,:,:),allocatable	:: rv_buf
	character(30)										:: readin_format,filename,input_file,initial_microstate_file
	character(400)										:: error_message

	call setup_command_arguments

	!Freqency at which to write out buffer to files
	writefreq = 100
	
    ! Get the periodic constrains form MD.in
    ! and the processor topology description
    open(1,file=trim(input_file))
    call locate(1,'PROCESSORS',.true.)
    read(1,*) npx
    read(1,*) npy
    read(1,*) npz
	close(1,status='keep')      !Close input file

	allocate(procnp(npx*npy*npz))
	allocate(npcount(npx,npy,npz))
	allocate(rv_buf(npx,npy,npz,writefreq,2*nd))
	
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
	read(2) nmonomers
	read(2) prev_npx !npx
	read(2) prev_npy !npy
	read(2) prev_npz !npz
	print'(2(a,3i8))', ' Previous Topology -- ',prev_npx,prev_npy,prev_npz ,' New Topology -- ',npx,npy,npz
	if (prev_npx .eq. npx .and. prev_npy.eq.npy.and.prev_npz .eq. npz) 	stop "No re-order needed"
	prev_nproc = prev_npx*prev_npy*prev_npz	
	do n=1,prev_nproc			!Loop through all processors and discard information
		read(2) checkint
	enddo
	read(2) density 
	read(2) rcutoff
	read(2) delta_t
	read(2) elapsedtime	
	read(2) k_c
	read(2) R_0
	read(2) delta_rneighbr
	close(2,status='keep')

	globaldomain(:) = initialnunits(:) &       !Size domain based on required density
						/((density/4)**(1.d0/nd)) 
	domain(1) = 	globaldomain(1)/npx
	domain(2) = 	globaldomain(2)/npy
	domain(3) = 	globaldomain(3)/npz
	halfdomain(:) = 0.5d0*domain(:)      !Useful defintion used later

	nl = 0; npcount = 0; procnp=0	!Reset local molecules count nl
	open(2,file=initial_microstate_file, form='unformatted',action='read', access='stream',position='rewind')
	!---------- For all molecule positions ------------
	!Move through location of position co-ordinates
	do n=1,globalnp

		read(2) buf		!Read position to buffer

		!Use integer division to determine which processor to assign molecule to
		pa(1) = ceiling((buf(1)+globaldomain(1)/2.d0)/domain(1))
		pa(2) = ceiling((buf(2)+globaldomain(2)/2.d0)/domain(2))
		pa(3) = ceiling((buf(3)+globaldomain(3)/2.d0)/domain(3))

		!If molecules is in the domain then add to processor's total
		npcount(pa(1),pa(2),pa(3)) = npcount(pa(1),pa(2),pa(3)) + 1 !Local molecule count
		nl = npcount(pa(1),pa(2),pa(3))
		!Read position and velocity
		rv_buf(pa(1),pa(2),pa(3),nl,1:6) = buf(1:6)

		!print'(7i8,6f10.5)', n, nl,pa,npcount, rv_buf(pa(1),pa(2),pa(3),nl,1:6)

		if (nl .eq. writefreq .or. n .eq. globalnp) then
			m = 0
			do i=1,npx
			do j=1,npy
			do k=1,npz
				write(filename,'(i3,a,i3,a,i3,a)'),i,'x',j,'x',k,'.temp'
				open(3,file=trim(filename), form='unformatted', action='write',access='stream',position='append')
				do l = 1,npcount(i,j,k)
					buf(1:3) = rv_buf(i,j,k,l,1:3); buf(4:6)= rv_buf(i,j,k,l,4:6)
					write(3) buf
				enddo
				close(3,status='keep')
				m = m + 1
				procnp(m) = procnp(m) + npcount(i,j,k)
			enddo
			enddo
			enddo
			!Reset write frequency counter
			nl = 0; npcount=0; rv_buf = 0.d0
			print'(a,f10.2)', & 
				'Redistributing molecules to new processor topology - % complete =', (100.d0*n/globalnp)
		endif
	enddo

	close(2,status='keep')
	deallocate(rv_buf)
	
	!Write to single file
	m = 0
	open(2,file='./final_state2', form='unformatted',action='write', access='stream',position='rewind',status='replace')
	do i=1,npx
	do j=1,npy
	do k=1,npz
		write(filename,'(i3,a,i3,a,i3,a)'),i,'x',j,'x',k,'.temp'
		m = m + 1
		! simpler version using inquire
		!inquire(file=filename,size=file_size)
		!PRINT*, PROCNP,2*ND*PROCNP(M)*8, file_size
		allocate(rv(2*nd*procnp(m)))
		open(3,file=trim(filename), form='unformatted',action='read', access='stream',position='rewind')
		read(3) rv(:)
		close(3,status='delete')
		!open(2,file='./final_state2', form='unformatted',action='write', access='stream',position='append')
		write(2) rv(:)
		deallocate(rv)
	enddo
	enddo
	enddo
	close(2,status='keep')

	!Write integer data at end of file	
	open(2,file='./final_state2', form='unformatted',access='stream',position='append')

	! header address
	inquire(2,POS=header_pos)
        
	write(2) globalnp           !Number of particles
	write(2) initialnunits(1) 	!x dimension split into number of cells
	write(2) initialnunits(2) 	!y dimension box split into number of cells
	write(2) initialnunits(3) 	!z dimension box split into number of cells
	write(2) Nsteps           	!Number of computational steps
	write(2) tplot            	!Frequency at which to record results
	write(2) seed(1)          	!Random number seed value 1
	write(2) seed(2)          	!Random number seed value 2
	write(2) periodic(1)  	 	!Boundary condition flags
	write(2) periodic(2)   		!Boundary condition flags
	write(2) periodic(3)	   	!Boundary condition flags
	write(2) potential_flag   	!Polymer/LJ potential flag
	write(2) nmonomers		   	!Polymer chain length
	write(2) npx                !Processors (npx) for new topology
	write(2) npy                !Processors (npy) for new topology
	write(2) npz                !Processors (npz) for new topology
	write(2) procnp				!Number of molecules per processors

	write(2) density          !Density of system
	write(2) rcutoff          !Cut off distance for particle interaction
	write(2) delta_t          !Size of time step
	write(2) elapsedtime      !Total elapsed time of all restarted simulations
	write(2) k_c		   	  !FENE spring constant
	write(2) R_0		      !FENE spring max elongation
	write(2) delta_rneighbr	  !Extra distance used for neighbour list cell size
    write(2) header_pos-1     ! -1 for MPI IO compatibility
	close(2,status='keep') 	  !Close final_state file

	close(2,status='keep')

	deallocate(procnp)
	deallocate(npcount)

contains

!=============================================================================
! setup_command_arguments
! Checks for command-line arguments passed to the program and assigns the 
! relevant values.
!-----------------------------------------------------------------------------
subroutine setup_command_arguments
	implicit none
	
	integer				:: i,argcount
	logical 			:: restart, restart_file_exists, input_file_exists
	character(len=200) 	:: arg,nextarg


	!Set default values in case they aren't specified by user
	restart = .false.			
	input_file_exists = .false.
	restart_file_exists = .false.
	input_file = 'MD.in'
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


end program change_proc_topology
