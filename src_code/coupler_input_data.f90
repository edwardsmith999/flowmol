module coupler_input_data

contains

subroutine read_coupler_input
	use coupler_module
	implicit none

	integer :: infileid=873457263
	logical :: found

	open(infileid,file='COUPLER.in',status="old",action="read", &
				  form="formatted")

	call locate(infileid,'DENSITY_CFD',found)
	if (found) then
		read(infileid,*) density_cfd
	else
		call error_abort("Density not specified in coupler input file.")
	end if
	
	call locate(infileid,'OVERLAP_EXTENTS',found)
	if (found) then
		read(infileid,*) icmin_olap
		read(infileid,*) icmax_olap
		read(infileid,*) jcmin_olap
		read(infileid,*) jcmax_olap
		read(infileid,*) kcmin_olap
		read(infileid,*) kcmax_olap
	else
		call error_abort("Ovelap extents unspecified in coupler input file.")
	end if

	call locate(infileid,'TIMESTEP_RATIO',found)
	if (found) then
		read(infileid,*) timestep_ratio !TODO name change
	else
		timestep_ratio = VOID
	end if
	
	call locate(infileid,'MATCH_CELLSIZE',found)
	if (found) then
		read(infileid,*) md_cfd_match_cellsize
	else
		md_cfd_match_cellsize = 0
	end if

	close(infileid,status="keep")

end subroutine read_coupler_input

subroutine locate(fileid,keyword,have_data)
	implicit none
	
	integer,intent(in)			:: fileid               ! File unit number
	character(len=*),intent(in)	:: keyword              ! Input keyword	
	logical,intent(out)			:: have_data            ! Flag: input found

	character*(100)				:: linestring           ! First 100 chars
	integer						:: keyword_length       ! Length of keyword
	integer						:: io                   ! File status flag

	keyword_length = len(keyword)
	rewind(fileid)
	
	! Loop until end of file or keyword found
	do
		! Read first 100 characters of line
		read (fileid,'(a)',iostat=io) linestring

		! If end of file is reached, exit
		if (io.ne.0) then 
			have_data = .false.
			exit
		end if
		
		! If the first characters match keyword, exit
		if (linestring(1:keyword_length).eq.keyword) then
			have_data = .true.
			exit
		endif

	end do

end subroutine locate

end module coupler_input_data
