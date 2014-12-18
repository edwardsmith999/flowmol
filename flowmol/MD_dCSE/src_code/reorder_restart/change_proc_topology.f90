! A program, run in serial which takes the current restart file and
! reorders it to allow restart with a different number of processors than the previous
! run. This may need updating if changes to the rest of the code have adjusted the 
! order of the header or the number of values stored...

program change_proc_topology
	implicit none

	logical :: molecules_outside=.false.
	integer :: i,j,k,m,n,pos,endpos
	integer :: globalnp,Nsteps,tplot,potential_flag,nmonomers
	integer :: checkint 
	integer :: prev_nproc,npx,npy,npz,prev_npx,prev_npy,prev_npz
	integer :: npxyz(3)
	integer :: rtrue_flag, solvent_flag
	integer :: tag
	integer,parameter :: nd = 3
    integer(selected_int_kind(18)) :: header_pos,end_pos ! 8 byte integer for header address
	integer,dimension(2) :: seed
	integer,dimension(3) :: initialnunits, pa,periodic
	integer,dimension(:),allocatable :: procnp, proctethernp
	integer,dimension(:,:,:),allocatable :: npcount, tethernpcount,fileunit
	double precision :: density,rcutoff,delta_t,elapsedtime,k_c,R_0,delta_rneighbr
	double precision :: simtime, eps_pp, eps_ps, eps_ss, dpbuf
	double precision,dimension(3) :: domain,halfdomain,globaldomain
	double precision,dimension(3) :: rbuf, vbuf, rtruebuf, rtetherbuf
	double precision,dimension(13) :: rv_buf
	character(30) :: input_file,initial_microstate_file
	integer, dimension(5), parameter :: tether_tags = (/3,5,6,7,10/)
	type monomer_info
		SEQUENCE !For MPI convenience
		integer :: chainID !Integer value: chain number 
		integer :: subchainID !Integer value: bead number
		integer :: funcy !Functionality of the monomer
		integer :: glob_no !Global molecule number
		integer :: bin_bflag(4) !Integer for bit manipulation to find 
                                !bond flags. For more info see function 
                                !get_bondflag
		! THE TOTAL NUMBER OF ITEMS IN THIS DATA TYPE MUST ALSO BE STORED 
        ! IN THE VARIABLE nsdmi
	end type monomer_info
    type(monomer_info) :: monomerbuf
    double precision, dimension(8) :: monomer_dpbuf

	call setup_command_arguments

    ! Get the periodic constrains form MD.in
    ! and the processor topology description
    open(1,file=trim(input_file))
    call locate(1,'PROCESSORS',.true.)
    read(1,*) npx
    read(1,*) npy
    read(1,*) npz
	close(1,status='keep')      !Close input file
	npxyz = (/npx,npy,npz/)

	allocate(procnp(npx*npy*npz))
	allocate(proctethernp(npx*npy*npz))
	allocate(npcount(npx,npy,npz))
	allocate(tethernpcount(npx,npy,npz))
	allocate(fileunit(npx,npy,npz))
	
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
	read(2) solvent_flag 
	read(2) nmonomers
	read(2) prev_npx !npx
	read(2) prev_npy !npy
	read(2) prev_npz !npz
	print'(2(a,3i8))', ' Previous Topology -- ',prev_npx,prev_npy,prev_npz ,' New Topology -- ',npx,npy,npz
	if (prev_npx .eq. npx .and. prev_npy.eq.npy.and.prev_npz .eq. npz) 	stop "No re-order needed"
	prev_nproc = prev_npx*prev_npy*prev_npz	
	do n=1,prev_nproc			!Loop through all processors and discard information
		read(2) checkint        !procnp
	enddo
	do n=1,prev_nproc			!Loop through all processors and discard information
		read(2) checkint        !proctethernp
	enddo
    read(2) globaldomain(1)
    if (globaldomain(1) .lt. 3.d0) then
        print*,  'Warning -- Since revision 673, globaldomain included in finalstate ', &
                 ' this appears to be an older restart file.', & 
                  '(or globaldomain in x = ',globaldomain(1),'). ', &
                  'It is assumed that global domain is not specified and density = ', globaldomain(1)
        density = globaldomain(1)
        globaldomain(:) = initialnunits(:) &       !Size domain based on required density
        					/((density/4)**(1.d0/nd)) 
    else
        read(2) globaldomain(2)
        read(2) globaldomain(3)
	    read(2) density 
    endif
	read(2) rcutoff
	read(2) delta_t
	read(2) elapsedtime	
	read(2) simtime 
	read(2) k_c
	read(2) R_0
	read(2) eps_pp 
	read(2) eps_ps 
	read(2) eps_ss 
	read(2) delta_rneighbr
	close(2,status='keep')

	domain(1) = 	globaldomain(1)/npx
	domain(2) = 	globaldomain(2)/npy
	domain(3) = 	globaldomain(3)/npz
	halfdomain(:) = 0.5d0*domain(:)      !Useful defintion used later

	npcount = 0; tethernpcount = 0; procnp=0; proctethernp=0	!Reset local molecules count 
	open(2,file=initial_microstate_file, form='unformatted',action='read', access='stream',position='rewind')

	m = 99
	do i = 1,npx
	do j = 1,npy
	do k = 1,npz
		m = m+1
		fileunit(i,j,k) = m
		open(m,form='unformatted',action='write',access='stream',status='replace',position='append')
	end do
	end do
	end do

	open (unit=6,form='formatted',carriagecontrol='fortran')

	print*, 'Reformatting restart file:'
	!---------- For all molecule positions ------------
	!Move through location of position co-ordinates
	do n=1,globalnp

		read(2) dpbuf; tag = nint(dpbuf)		!Read position to buffer
		read(2) rbuf
		read(2) vbuf
		if (rtrue_flag.eq.1) then
			read(2) rtruebuf
		end if
		if (any(tag.eq.tether_tags)) then
			read(2) rtetherbuf  
		end if
		if (potential_flag.eq.1) then
            !Read monomer data
            read(2) monomer_dpbuf
        end if


		!Use integer division to determine which processor to assign molecule to
		pa(1) = ceiling((rbuf(1)+globaldomain(1)/2.d0)/domain(1))
		pa(2) = ceiling((rbuf(2)+globaldomain(2)/2.d0)/domain(2))
		pa(3) = ceiling((rbuf(3)+globaldomain(3)/2.d0)/domain(3))


		!Capture molecules that are slightly outside the domain and print a warning
		if (any(pa.eq.0) .or. any(pa.gt.npxyz)) then
			where (pa .eq. 0) pa = 1
			where (pa .gt. npxyz) pa = npxyz
			molecules_outside = .true.
		end if

		!Read position and velocity
		pos = 1	
		rv_buf(pos) = dble(tag); pos = pos + 1
		rv_buf(pos:pos+2) = rbuf; pos = pos + 3
		rv_buf(pos:pos+2) = vbuf; pos = pos + 3
		if (rtrue_flag.eq.1) then
			rv_buf(pos:pos+2) = rtruebuf; pos = pos + 3
		end if
		if (any(tag.eq.tether_tags)) then
			rv_buf(pos:pos+2) = rtetherbuf; pos = pos + 3
		end if
		endpos = pos - 1

		write(fileunit(pa(1),pa(2),pa(3))) rv_buf(1:endpos)
        if (potential_flag .eq. 1) then
            write(fileunit(pa(1),pa(2),pa(3))) monomer_dpbuf
        end if

		npcount(pa(1),pa(2),pa(3)) = npcount(pa(1),pa(2),pa(3)) + 1
		if (any(tag.eq.tether_tags)) then
			tethernpcount(pa(1),pa(2),pa(3)) = tethernpcount(pa(1),pa(2),pa(3)) + 1
		end if

        if (mod(100*(n/globalnp),10)) call progress(nint(100.d0*n/globalnp))

	end do

	m = 1
	do i = 1,npx
	do j = 1,npy
	do k = 1,npz

		procnp(m) = npcount(i,j,k)
		proctethernp(m) = tethernpcount(i,j,k)
		m = m + 1

		close(fileunit(i,j,k),status='keep')

	end do
	end do
	end do

	if (molecules_outside .eq. .true.) then
		print*, ''
		print*, ''
		print*, ''
		print*, '       *****************************     '
		print*, 'WARNING: At least one molecule was found outside the domain, '
		print*, 'so was/were assigned to the nearest processor(s).'
		print*, '       *****************************     '
		print*, ''
		print*, ''
		print*, ''
	end if

	call system('cat fort.* > final_state2')
	call system('rm fort.*')

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
	write(2) rtrue_flag         !
	write(2) solvent_flag       !Solvent on/off flag
	write(2) nmonomers		   	!Polymer chain length
	write(2) npx                !Processors (npx) for new topology
	write(2) npy                !Processors (npy) for new topology
	write(2) npz                !Processors (npz) for new topology
	write(2) procnp				!Number of molecules per processors
	write(2) proctethernp				!Number of molecules per processors

    write(2) globaldomain(1)
    write(2) globaldomain(2)
    write(2) globaldomain(3)
	write(2) density          !Density of system
	write(2) rcutoff          !Cut off distance for particle interaction
	write(2) delta_t          !Size of time step
	write(2) elapsedtime      !Total elapsed time of all restarted simulations
	write(2) simtime          !Total elapsed time of all restarted simulations
	write(2) k_c		   	  !FENE spring constant
	write(2) R_0		      !FENE spring max elongation
	write(2) eps_pp           !Soddemann potential parameter
	write(2) eps_ps           !Soddemann potential parameter
	write(2) eps_ss           !Soddemann potential parameter
	write(2) delta_rneighbr	  !Extra distance used for neighbour list cell size
    write(2) header_pos-1     ! -1 for MPI IO compatibility
	close(2,status='keep') 	  !Close final_state file

	close(2,status='keep')

	deallocate(procnp)
	deallocate(proctethernp)
	deallocate(npcount)
	deallocate(tethernpcount)

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

end program change_proc_topology
