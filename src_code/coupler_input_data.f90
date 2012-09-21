!=============================================================================
!
!						 Coupler input
!
! Modules containing initialisation data for CFD and MD which 
! can be read from coupler input file.
!
! This offers a consistent set of parameters for both CFD and MD when 
! running coupled.
!
! Some flexibility is provided. If some keywords from coupler.in or the whole file is missing
! the couper will try to use CFD and MD parammeters from their respective input files with 
! some adjustments 
! 
!
module coupler_input_data
	implicit none

	type cfd_domain_sizes
		SEQUENCE 							! useful for MPI
		integer tag 						! tells from which input file domain values are taken
		integer ndim
		character(len=8) :: units, cell_type
		real(kind=kind(0.d0)) x,y,z
	end type cfd_domain_sizes

	type cfd_ncell_counts
		SEQUENCE 				! useful for MPI
		integer tag 			! tells from which input file domain values are taken
		integer x,y,z			! number of cells in CFD physical domain
	end type cfd_ncell_counts

	type cfd_md_overlap
		SEQUENCE 				! useful for MPI
		integer tag 			! tells from which input file domain values are taken
		integer y_overlap 		! number of cells that overlaps CFD physical domain with MD domain
	end type cfd_md_overlap

	type cfd_parameters
		SEQUENCE
		type(cfd_domain_sizes) domain 
		type(cfd_ncell_counts) ncells
		type(cfd_md_overlap)   overlap
	end type cfd_parameters

	type section_type
		integer, pointer :: tag => null()
		character(len=64) str
	end type section_type

	! data that can be provided in coupler.in
	! physical
	real(kind=kind(0.d0)) density
	integer,target :: density_tag
	!CFD
	type(cfd_parameters), target :: cfd_coupler_input 
	! MD
	real(kind=kind(0.d0)) 	:: md_ly_extension 			! sigma units so far
	integer, target 		:: md_ly_extension_tag		! MD extension below CFD grid in y direction		  
	integer 				:: md_average_period		! collect data for velocity average every ... MD step
	integer, target 		:: md_average_period_tag
	integer 				:: md_save_period			! save data for velocity profile every ... CFD step
	integer, target 		:: md_save_period_tag	

	!stop_tag used in request_abort - declared in coupler_module

	! staggered_averages - declared in coupler_module
	
	integer				 :: md_steps_per_dt_cfd	  ! number of MD steps per CFD step
	integer, target		 :: md_steps_per_dt_cfd_tag

	integer 			 :: md_cfd_match_cellsize		!Force MD cells to be an integer multiple of CFD cell size
	integer, target 	 :: md_cfd_match_cellsize_tag

	character(len=64)	:: cfd_code_name ! used to set cfd_code_id in MD
	integer		   		:: cfd_code_id
	integer, target		:: cfd_code_id_tag

	! auxiliary list for easy manipulation
	integer,parameter 		:: nsections=12 ! total number of section in input files
	type(section_type) section(nsections)

contains

subroutine read_coupler_input
	use mpi
	use coupler_parameters
	use coupler_module
	implicit none 

	integer ndim, myid, myid_cfd, iroot_global, ierr
	logical have_input

	call setup_inputs(have_input)

	if (have_input) then		   
		if (myid_cfd == 0) then
			! Read input data on rank 0 and broadcast
			call read_input_file
			call pack_bcast_input			  
		else
			! receive input data from rank 0
			call bcast_unpack_input
		endif

		call check_inputs
	else
		if (myid_cfd == 0 ) then
			write(0,*) "WARNING: No coupler input file, will use CFD domain data and adjust accordingly"
		endif
	endif

contains

	subroutine setup_inputs(have_input)
		implicit none

		logical,intent(out)		:: have_input

		! find rank 0 in CFD_COMM and open COUPLER.in file on it
		call mpi_comm_rank(COUPLER_GLOBAL_COMM, myid, ierr)
		myid_cfd = -1; iroot_global = -1
		if (COUPLER_REALM == COUPLER_CFD) then
			call mpi_comm_rank(COUPLER_REALM_COMM, myid_cfd,ierr)
		endif

		have_input = .false.
		if (myid_cfd == 0) then 
			inquire(file="COUPLER.in",exist=have_input)
			iroot_global = myid
		endif

		! I cannot bcast directly have_input because it is not guaranteed that
		! myid == myid_cfd == 0
		call mpi_allreduce(MPI_IN_PLACE,have_input,1,MPI_LOGICAL,MPI_LOR,COUPLER_GLOBAL_COMM,ierr)
		! set the root for the global communicator
		call mpi_allreduce(MPI_IN_PLACE,iroot_global,1,MPI_INTEGER,MPI_MAX,COUPLER_GLOBAL_COMM,ierr)

		! set default values for coupler input data
		section(1)%str = "DENSITY"
		section(1)%tag => density_tag
		density_tag  = VOID
		density	  = 0.d0

		section(2)%str = "CFD_DOMAIN"
		section(2)%tag => cfd_coupler_input%domain%tag
		cfd_coupler_input%domain%tag   		= VOID
		cfd_coupler_input%domain%ndim  		= 0
		cfd_coupler_input%domain%units 		= "sigma"
		cfd_coupler_input%domain%cell_type 	= "fcc"
		cfd_coupler_input%domain%x	 		= 0.d0
		cfd_coupler_input%domain%y	 		= 0.d0
		cfd_coupler_input%domain%z	 		= 0.d0

		section(3)%str = "CFD_NCELLS"
		section(3)%tag => cfd_coupler_input%ncells%tag
		cfd_coupler_input%ncells%tag   = VOID
		cfd_coupler_input%ncells%x	 = 0
		cfd_coupler_input%ncells%y	 = 0
		cfd_coupler_input%ncells%z	 = 0

		section(4)%str = "CFD_OVERLAP"
		section(4)%tag => cfd_coupler_input%overlap%tag
		cfd_coupler_input%overlap%y_overlap = 0

		! No initialisation of the MD values associated with the tags because 
		!defaults are set in coupler_internal_md module		 
		section(5)%str = "MD_LY_EXTENSION"
		section(5)%tag => md_ly_extension_tag 
		md_ly_extension_tag   = VOID

		section(6)%str = "MD_AVERAGE_PERIOD"
		section(6)%tag => md_average_period_tag
		md_average_period_tag = VOID

		section(7)%str = "MD_SAVE_PERIOD"
		section(7)%tag => md_save_period_tag
		md_save_period_tag	= VOID

		section(8)%str = "STOP_TAG"
		section(8)%tag => stop_request_tag
		stop_request_tag = VOID

		section(9)%str = "STAGGERED_AVERAGES"
		section(9)%tag => staggered_averages_tag
		staggered_averages_tag = VOID

		section(10)%str = "CFD_CODE_ID"
		section(10)%tag => cfd_code_id_tag
		cfd_code_id_tag = VOID

		section(11)%str = "MD_STEPS_PER_DT_CFD"
		section(11)%tag => md_steps_per_dt_cfd_tag
		section(11)%tag = VOID
		!Force MD cells to be an integer multiple of CFD cell size
		section(12)%str = "MD_CFD_MATCH_CELLSIZE"		
		section(12)%tag => md_cfd_match_cellsize_tag
		section(12)%tag = VOID
		md_cfd_match_cellsize = 0

	end subroutine setup_inputs

	subroutine read_input_file
		use coupler_parameters
		implicit none

		integer io,i
		logical have_data
		character(len=100) linein,line

		open(34,file="COUPLER.in", status="old", action="read")

		do 
			read(34,*,iostat=io) linein
			if (io /= 0) exit

			line=adjustl(linein)
			if(line(1:1) == "#") cycle

			select case (line)
			case ("DENSITY","density","Density")
				density_tag = CPL
				read(34,*) density
			case ("CFD_DOMAIN","Cfd_Domain","Cfd_domain","cfd_domain")
				cfd_coupler_input%domain%tag = CPL
				read(34,*) ndim 
				cfd_coupler_input%domain%ndim = ndim
				if ( ndim < 1 .or. ndim > 3) then 
					write(*,*) "Wrong value for CFD number of dimensions", ndim
					call MPI_Abort(COUPLER_GLOBAL_COMM,COUPLER_ERROR_INPUT_FILE, ierr)
				endif
				read(34,*) cfd_coupler_input%domain%units
					! Test if units are known ? 
				read(34,*) cfd_coupler_input%domain%cell_type
				read(34,*) cfd_coupler_input%domain%x
				if (ndim > 1) then
					read(34,*) cfd_coupler_input%domain%y
				endif
				if (ndim > 2) then
					read(34,*) cfd_coupler_input%domain%z
				endif
			case("CFD_NCELLS","Cfd_Ncells","Cfd_ncells","cfd_ncells")
				cfd_coupler_input%ncells%tag = CPL
				read(34,*) cfd_coupler_input%ncells%x
				if (ndim > 1) then 
					read(34,*) cfd_coupler_input%ncells%y
				endif
				if (ndim > 2) then 
					read(34,*) cfd_coupler_input%ncells%z
				endif
			case("OVERLAP_CELLS","Overlap_Cells","Overlap_cells","overlap_cells")
				cfd_coupler_input%overlap%tag = CPL
				!if (ndim > 1) then 
				read(34,*) cfd_coupler_input%overlap%y_overlap
				!endif
			case("MD_LY_EXTENSION","md_ly_extension","Md_ly_extension")
				md_ly_extension_tag = CPL
				read(34,*) md_ly_extension
			case("MD_AVERAGE_PERIOD","md_average_period")
				md_average_period_tag = CPL
				read(34,*) md_average_period
			case("MD_SAVE_PERIOD","md_save_period")
				md_save_period_tag = CPL
				read(34,*) md_save_period
			case("STOP_REQUEST","Stop_request","stop_request")
				stop_request_tag = CPL
				read(34,*) stop_request_name
				case("STAGGERED_AVERAGES","Staggered_averages", "staggered_averages")
					staggered_averages_tag = CPL
					read(34,*) staggered_averages
				case("CFD_CODE_ID","cfd_code_id","Cfd_code_id")
					cfd_code_id_tag = CPL
					read(34,*) cfd_code_name
					select case(cfd_code_name)
					case ("COUETTE_SERIAL","couette_serial")
						cfd_code_id = couette_serial
					case("COUETTE_PARALLEL","couette_parallel")
						cfd_code_id = couette_parallel
					case default
						write(0,*) "unknow CFD_CODE_NAME, check COUPLER.in file"
						call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_INPUT_FILE,ierr)
					end select
				case("MD_STEPS_PER_DT_CFD","md_steps_per_dt_cfd")
					md_steps_per_dt_cfd_tag = CPL
					read(34,*) md_steps_per_dt_cfd
				case("MD_CFD_MATCH_CELLSIZE","md_cfd_match_cellsize")
					md_cfd_match_cellsize_tag = CPL
					read(34,*) md_cfd_match_cellsize
				case default
				! unkown input section abort
				write(0,*) "unknow coupler input section,",line,", check COUPLER.in file"
				call MPI_Abort(MPI_COMM_WORLD,COUPLER_ERROR_INPUT_FILE,ierr)
			end select
		enddo
		close(34)

		!report missing sections from input
		do i=1,nsections
			if (section(i)%tag == VOID ) then 
				write(0,*) " WARNING: no ", trim(section(i)%str), "  section found in coupler input file"
			endif
		end do

		write(0,*) "read COUPLER.in done"

	end subroutine read_input_file

	subroutine pack_bcast_input
		implicit none

		integer, parameter :: sbuff = 1024,scaux=64
		integer position, ierr
		character(len=sbuff) buffer
		character(len=scaux) caux

		position = 0
		call mpi_pack((/ cfd_coupler_input%domain%tag, cfd_coupler_input%domain%ndim /),2,MPI_INTEGER, &
						buffer,sbuff,position,COUPLER_REALM_COMM,ierr)
		call mpi_pack(cfd_coupler_input%domain%units,8,MPI_CHARACTER,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)
		call mpi_pack(cfd_coupler_input%domain%cell_type,8,MPI_CHARACTER,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)
		call mpi_pack(cfd_coupler_input%domain%x,1,MPI_DOUBLE_PRECISION,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)
		call mpi_pack(cfd_coupler_input%domain%y,1,MPI_DOUBLE_PRECISION,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)
		call mpi_pack(cfd_coupler_input%domain%z,1,MPI_DOUBLE_PRECISION,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)

		call mpi_pack(cfd_coupler_input%ncells%tag,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)
		call mpi_pack(cfd_coupler_input%ncells%x,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)
		call mpi_pack(cfd_coupler_input%ncells%y,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)
		call mpi_pack(cfd_coupler_input%ncells%z,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)
		call mpi_pack(cfd_coupler_input%overlap%y_overlap,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)

		call mpi_pack(density,1,MPI_DOUBLE_PRECISION,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)
		call mpi_pack(density_tag,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_REALM_COMM,ierr)

		call mpi_bcast(position,1,MPI_INTEGER,0,COUPLER_REALM_COMM,ierr)
		call mpi_bcast(buffer,position,MPI_PACKED,0,COUPLER_REALM_COMM,ierr)

		! pack MD data. This can be optimized, for a start we send variables with a non-void tag 
		position = 0
		if (md_ly_extension_tag == CPL) then
			write(caux,'(a)')"MD_LY_EXTENSION"
			call mpi_pack(caux,scaux,MPI_CHARACTER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(md_ly_extension_tag,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(md_ly_extension	,1,MPI_DOUBLE_PRECISION,buffer,sbuff,position,COUPLER_ICOMM,ierr)
		endif

		if (md_average_period_tag == CPL) then
			write(caux,'(a)')"MD_AVERAGE_PERIOD"
			call mpi_pack(caux,scaux,MPI_CHARACTER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(md_average_period_tag,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(md_average_period	,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
		endif

		if (md_save_period_tag == CPL) then
			write(caux,'(a)')"MD_SAVE_PERIOD"
			call mpi_pack(caux,scaux,MPI_CHARACTER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(md_save_period_tag,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(md_save_period	,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
		endif

		if (cfd_code_id_tag == CPL) then
			write(caux,'(a)')"CFD_CODE_ID"
			call mpi_pack(caux,scaux,MPI_CHARACTER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(cfd_code_id_tag,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(cfd_code_id	,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
		endif

		if (md_steps_per_dt_cfd_tag == CPL) then
			write(caux,'(a)')"MD_STEPS_PER_DT_CFD"
			call mpi_pack(caux,scaux,MPI_CHARACTER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(md_steps_per_dt_cfd_tag,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(md_steps_per_dt_cfd	,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
		endif

		if (md_cfd_match_cellsize_tag == CPL) then
			write(caux,'(a)')"MD_CFD_MATCH_CELLSIZE_TAG"
			call mpi_pack(caux,scaux,MPI_CHARACTER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(md_cfd_match_cellsize_tag,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
			call mpi_pack(md_cfd_match_cellsize	,1,MPI_INTEGER,buffer,sbuff,position,COUPLER_ICOMM,ierr)
		endif

		! broadcast the packed data to MD realm
		call mpi_bcast(position,1,MPI_INTEGER,0,COUPLER_REALM_COMM,ierr)
		call mpi_bcast(position,1,MPI_INTEGER,MPI_ROOT,COUPLER_ICOMM,ierr)
		if ( position > 0) then
			call mpi_bcast(buffer,position,MPI_PACKED,MPI_ROOT,COUPLER_ICOMM,ierr)
		endif

	end subroutine pack_bcast_input

	subroutine bcast_unpack_input
 		implicit none

		integer, parameter :: sbuff = 1024, scaux=64
		integer position, count, ierr
		character(len=sbuff) buffer
		character(len=scaux) caux

		if (COUPLER_REALM == COUPLER_CFD) then 
			
			call mpi_bcast(count,1,MPI_INTEGER,0,COUPLER_REALM_COMM,ierr)
			call mpi_bcast(buffer,count,MPI_PACKED,0,COUPLER_REALM_COMM,ierr)
			position=0
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%domain%tag,1,MPI_INTEGER,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%domain%ndim,1,MPI_INTEGER,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%domain%units,8,MPI_CHARACTER,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%domain%cell_type,8,MPI_CHARACTER,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%domain%x,1,MPI_DOUBLE_PRECISION,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%domain%y,1,MPI_DOUBLE_PRECISION,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%domain%z,1,MPI_DOUBLE_PRECISION,COUPLER_REALM_COMM,ierr)

			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%ncells%tag,1,MPI_INTEGER,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%ncells%x,1,MPI_INTEGER,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%ncells%y,1,MPI_INTEGER,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%ncells%z,1,MPI_INTEGER,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,cfd_coupler_input%overlap%y_overlap,1,MPI_INTEGER,COUPLER_REALM_COMM,ierr)

			call mpi_unpack(buffer,sbuff,position,density,1,MPI_DOUBLE_PRECISION,COUPLER_REALM_COMM,ierr)
			call mpi_unpack(buffer,sbuff,position,density_tag,1,MPI_INTEGER,COUPLER_REALM_COMM,ierr)

			! intercommunicator null sends
			call mpi_bcast(count,1,MPI_INTEGER,0,COUPLER_REALM_COMM,ierr)
			call mpi_bcast(count,1,MPI_INTEGER,MPI_PROC_NULL,COUPLER_ICOMM,ierr)
			if (count > 0) then
				call mpi_bcast(buffer,count,MPI_PACKED,MPI_PROC_NULL,COUPLER_ICOMM,ierr)
			endif
		else 
			! MD side receives
			call mpi_bcast(count,1,MPI_INTEGER,0,COUPLER_ICOMM,ierr)
			if (count > 0) then
				call mpi_bcast(buffer,count,MPI_PACKED,0,COUPLER_ICOMM,ierr)
			endif
			if ( count > 0 ) then 
				position = 0
				call mpi_unpack(buffer,sbuff,position,caux,scaux,MPI_CHARACTER,COUPLER_ICOMM,ierr)
				do 
					select case (caux)
					case("MD_LY_EXTENSION")
						call mpi_unpack(buffer,sbuff,position,md_ly_extension_tag,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
						call mpi_unpack(buffer,sbuff,position,md_ly_extension	,1,MPI_DOUBLE_PRECISION,COUPLER_ICOMM,ierr)
					case("MD_AVERAGE_PERIOD")
						call mpi_unpack(buffer,sbuff,position,md_average_period_tag,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
						call mpi_unpack(buffer,sbuff,position,md_average_period   ,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
					case("MD_SAVE_PERIOD")
						call mpi_unpack(buffer,sbuff,position,md_save_period_tag,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
						call mpi_unpack(buffer,sbuff,position,md_save_period	,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
					case("CFD_CODE_ID")
						call mpi_unpack(buffer,sbuff,position,cfd_code_id_tag,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
						call mpi_unpack(buffer,sbuff,position,cfd_code_id	,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
					case("MD_STEPS_PER_DT_CFD")
						call mpi_unpack(buffer,sbuff,position,md_steps_per_dt_cfd_tag,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
						call mpi_unpack(buffer,sbuff,position,md_steps_per_dt_cfd	,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
					case("MD_CFD_MATCH_CELLSIZE_TAG")
						call mpi_unpack(buffer,sbuff,position,md_cfd_match_cellsize_tag,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
						call mpi_unpack(buffer,sbuff,position,md_cfd_match_cellsize	,1,MPI_INTEGER,COUPLER_ICOMM,ierr)
					case default
						call MPI_Abort(COUPLER_GLOBAL_COMM,COUPLER_ERROR_READ_INPUT,ierr)
					end select
					if (position >= count) exit 
					call mpi_unpack(buffer,sbuff,position,caux,scaux,MPI_CHARACTER,COUPLER_ICOMM,ierr)
				enddo

			endif
		endif

	end subroutine bcast_unpack_input


	subroutine check_inputs
		implicit none

		! broadcast stop_request across communicator if section id present in COUPLER.in
		call mpi_bcast(stop_request_tag,1,MPI_INTEGER,iroot_global,COUPLER_GLOBAL_COMM,ierr)
		if (stop_request_tag == CPL) then 
			call mpi_bcast(stop_request_name,len(stop_request_name),MPI_CHARACTER,iroot_global,COUPLER_GLOBAL_COMM,ierr)
			stop_request_activated = .true.
		endif

		! the save for staggered averages
		call mpi_bcast(staggered_averages_tag,1,MPI_INTEGER,iroot_global,COUPLER_GLOBAL_COMM,ierr)
		if (staggered_averages_tag == CPL) then 
			call mpi_bcast(staggered_averages,size(staggered_averages),MPI_LOGICAL,iroot_global,COUPLER_GLOBAL_COMM,ierr)
		endif

	end subroutine check_inputs

end subroutine read_coupler_input


subroutine locate(fileid,keyword,have_data)
	implicit none
	
	integer,intent(in)			:: fileid			! File unit number
	character(len=*),intent(in)	:: keyword			! Input keyword	
	logical,intent(out)			:: have_data		! Flag to check if input is required

	character*(100)				:: linestring		! First 100 characters in a line
	integer						:: keyword_length	! Length of input keyword
	integer						:: io				! File status flag

	have_data = .false.
	keyword_length = len(keyword)
	rewind(fileid)	! Go to beginning of input file
	do
	   read (fileid,'(a)',iostat=io) linestring				! Read first 100 characters
	   if (io.ne.0) exit									! If end of file is reached
	   if (linestring(1:keyword_length).eq.keyword) then	! If the first characters match keyword then exit loop
		  have_data = .true.
		  exit
	   endif
	end do

end subroutine locate


end module coupler_input_data
