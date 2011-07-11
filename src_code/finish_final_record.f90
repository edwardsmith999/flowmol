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

end module module_final_record
!----------------------------------------------------------------------------------

subroutine finish_final_record
use module_final_record
implicit none

	!Print out final results from the simulation
	if (irank .eq. iroot) then
		print*, 'Results from last iteration of simulation'
		print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f19.15,a,f19.15,a,f19.15,a,f10.4)', &
			Nsteps,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy,';',totenergy,';',pressure
	endif

	!Set up final print out table
	!print '(8a)', &
	!'Final Iteration; Simulation time;  Mean Temp;   Mean KE; SD KE; Mean PE; SD PE; Mean TE; SD TE; Mean Pressure;'

	!Write values of distribution functions
	!write(12,'(a)') 'Velocity frequency distribution'
	!do n=1,nbins(1)
	!	write(12,'(2(f10.5))') (n-0.5d0)*binsize, normalisedbin(n) 
	!enddo

	!write(12,'(a)') 'Radial distribution function'
	!do n = 1, nshells
	!	write(12,'(2(f10.5))') (n-0.5d0)*delta_r, RDF(n)
	!enddo

	call messenger_syncall() !Make sure all processes have finished writing

	!Write simualtion properties and final position & velocity to 
	!unformatted output file to allow restart
	call parallel_io_final_state

	!Reformat positions recorded into the correct from for VMD to use
	if(irank .eq. iroot) then
		if(vmd_outflag .eq. 1) call reformat_dcd
		if(vmd_outflag .eq. 2) call reformat_dcd_sl
		if(vmd_outflag .eq. 3) call reformat_dcd_halo
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

end subroutine finish_final_record

!--------------------------------------------------------------------------------------------
!Re-format the file into binary with header used by .dcd file format

!Header added by David Trevelyan, DCD file format from:
!http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
!Information at the above url is a little misleading in places, comments below try to clarify. 


subroutine reformat_dcd
	use module_final_record
	implicit none

	integer				:: n, i			!--Dummy variables
	integer				:: NSET			!--Number of frames
	integer				:: ISTRT		!--Starting frame
	integer				:: NSAVC		!--Number of frames per coordinate save
	integer				:: NATOMNFREAT		!--Number of fixed atoms
	integer				:: NTITLE		!--Number of 80-character strings in title (set as 2)
	integer				:: NATOM		!--Number of atoms
	integer, dimension (5)		:: FIVEZ		!--According to documentation, five zeros, but first is actually NSET
	integer, dimension (9)		:: NINEZ		!--Nine zeros
	character(len=4)		:: HDR			!--Header string, value 'CORD'
	character(len=80), dimension(2)	:: TITLE	        !--Title(s)
	real,allocatable,dimension(:)   :: Xbuf, Ybuf, Zbuf	!--Buffers used to copy from direct access to binary
	double precision		:: DELTA		!--Timestep between frames


	print*, 'Generating final VMD.dcd ouput file - for large systems this may take some time'

	!Determine size of file datatype
	!inquire(file='testfile.dcd', recl=datasize)
	!print*, 'datasize', datasize
 	!call MPI_type_size(MPI_real,datasize,ierr)
	!print*, 'datasize', datasize

	!Set header information	
	HDR		=	'CORD'				!header text
	NSET		=	(Nsteps/tplot)			!number of recorded frames
	ISTRT		=	0				!the starting timestep
	NSAVC		=	1				!number of timesteps between dcd frame saves
	FIVEZ(1)	=	NSET				!not sure why
	FIVEZ(2:5)	=	0				!buffer zeros
	NATOMNFREAT	=	0				!number of fixed atoms?
	DELTA		=	delta_t				!delta_t (x-plor is double, charmm is real)
	NINEZ(:)	=	0				!buffer zeros
	NTITLE		=	2				!number of 80-character strings in title
	TITLE(1)	=	'  Simulation record file '	!
	TITLE(2)        =	'   Written in parallel   '	!
	NATOM		=	globalnp			!number of particles

	allocate(Xbuf(NSET*globalnp))
	allocate(Ybuf(NSET*globalnp))
	allocate(Zbuf(NSET*globalnp))

	!Read position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file="results/vmd_temp.dcd",access='direct',recl=1)

	do i=1,NSET
	do n=1,globalnp
		read(17,rec=(i-1)*nd*globalnp+n) Xbuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+globalnp) Ybuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+2*globalnp) Zbuf(n+globalnp*(i-1))
	enddo
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file="results/vmd_out.dcd",status='replace', form="unformatted")
	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	do i=1,NSET
		write(3) Xbuf((i-1)*globalnp+1:i*globalnp)
		write(3) Ybuf((i-1)*globalnp+1:i*globalnp)
		write(3) Zbuf((i-1)*globalnp+1:i*globalnp)
	enddo

	close(3,status='keep')

end subroutine reformat_dcd


!--------------------------------------------------------------------------------------------
! SPLIT SOLID AND LIQUID REGIONS INTO SEPERATE MOLECULE FILES
!Re-format the file into binary with header used by .dcd file format

!Header added by David Trevelyan, DCD file format from:
!http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
!Information at the above url is a little misleading in places, comments below try to clarify. 


subroutine reformat_dcd_sl
	use module_final_record
	implicit none

	integer				:: n, i			!--Dummy variables
	integer				:: NSET			!--Number of frames
	integer				:: ISTRT		!--Starting frame
	integer				:: NSAVC		!--Number of frames per coordinate save
	integer				:: NATOMNFREAT		!--Number of fixed atoms
	integer				:: NTITLE		!--Number of 80-character strings in title (set as 2)
	integer				:: NATOM		!--Number of atoms
	integer, dimension (5)		:: FIVEZ		!--According to documentation, five zeros, but first is actually NSET
	integer, dimension (9)		:: NINEZ		!--Nine zeros
	character(len=4)		:: HDR			!--Header string, value 'CORD'
	character(len=80), dimension(2)	:: TITLE	        !--Title(s)
	real,allocatable,dimension(:)   :: Xbuf, Ybuf, Zbuf	!--Buffers used to copy from direct access to binary
	double precision		:: DELTA		!--Timestep between frames


	print*, 'Generating final VMD.dcd ouput file - for large systems this may take some time'

	!Determine size of file datatype
	!inquire(file='testfile.dcd', recl=datasize)
	!print*, 'datasize', datasize
 	!call MPI_type_size(MPI_real,datasize,ierr)
	!print*, 'datasize', datasize

	!Set header information	
	HDR		=	'CORD'				!header text
	NSET		=	(Nsteps/tplot)			!number of recorded frames
	ISTRT		=	0				!the starting timestep
	NSAVC		=	1				!number of timesteps between dcd frame saves
	FIVEZ(1)	=	NSET				!not sure why
	FIVEZ(2:5)	=	0				!buffer zeros
	NATOMNFREAT	=	0				!number of fixed atoms?
	DELTA		=	delta_t				!delta_t (x-plor is double, charmm is real)
	NINEZ(:)	=	0				!buffer zeros
	NTITLE		=	2				!number of 80-character strings in title
	TITLE(1)	=	'  Simulation record file '	!
	TITLE(2)        =	'   Written in parallel   '	!
	NATOM		=	globalnp			!number of particles

	allocate(Xbuf(NSET*globalnp))
	allocate(Ybuf(NSET*globalnp))
	allocate(Zbuf(NSET*globalnp))

	!Read Solid molecule position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file="results/vmd_solid_temp.dcd",access='direct',recl=1)

	do i=1,NSET
	do n=1,globalnp
		read(17,rec=(i-1)*nd*globalnp+n) Xbuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+globalnp) Ybuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+2*globalnp) Zbuf(n+globalnp*(i-1))
	enddo
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file="results/vmd_solid_out.dcd",status='replace', form="unformatted")
	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	do i=1,NSET
		write(3) Xbuf((i-1)*globalnp+1:i*globalnp)
		write(3) Ybuf((i-1)*globalnp+1:i*globalnp)
		write(3) Zbuf((i-1)*globalnp+1:i*globalnp)
	enddo

	close(3,status='keep')

	!Read liquid molecule position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file="results/vmd_liquid_temp.dcd",access='direct',recl=1)

	do i=1,NSET
	do n=1,globalnp
		read(17,rec=(i-1)*nd*globalnp+n) Xbuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+globalnp) Ybuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+2*globalnp) Zbuf(n+globalnp*(i-1))
	enddo
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file="results/vmd_liquid_out.dcd",status='replace', form="unformatted")
	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	do i=1,NSET
		write(3) Xbuf((i-1)*globalnp+1:i*globalnp)
		write(3) Ybuf((i-1)*globalnp+1:i*globalnp)
		write(3) Zbuf((i-1)*globalnp+1:i*globalnp)
	enddo

	close(3,status='keep')

end subroutine reformat_dcd_sl

!--------------------------------------------------------------------------------------------
!Re-format halo cell output for vmd ~ useful for debuggin purposes

subroutine reformat_dcd_halo
	use module_final_record
	implicit none

	integer				:: n, i			!--Dummy variables
	integer				:: NSET			!--Number of frames
	integer				:: ISTRT		!--Starting frame
	integer				:: NSAVC		!--Number of frames per coordinate save
	integer				:: NATOMNFREAT		!--Number of fixed atoms
	integer				:: NTITLE		!--Number of 80-character strings in title (set as 2)
	integer				:: NATOM		!--Number of atoms
	integer, dimension (5)		:: FIVEZ		!--According to documentation, five zeros, but first is actually NSET
	integer, dimension (9)		:: NINEZ		!--Nine zeros
	character(len=4)		:: HDR			!--Header string, value 'CORD'
	character(len=80), dimension(2)	:: TITLE	        !--Title(s)
	real,allocatable,dimension(:)   :: Xbuf, Ybuf, Zbuf	!--Buffers used to copy from direct access to binary
	double precision		:: DELTA		!--Timestep between frames


	print*, 'Generating final halo VMD.dcd ouput file - for large systems this may take some time'

	!Determine size of file datatype
	!inquire(file='testfile.dcd', recl=datasize)
	!print*, 'datasize', datasize
 	!call MPI_type_size(MPI_real,datasize,ierr)
	!print*, 'datasize', datasize

	!Set header information	
	HDR		=	'CORD'				!header text
	NSET		=	(Nsteps/tplot)			!number of recorded frames
	ISTRT		=	0				!the starting timestep
	NSAVC		=	1				!number of timesteps between dcd frame saves
	FIVEZ(1)	=	NSET				!not sure why
	FIVEZ(2:5)	=	0				!buffer zeros
	NATOMNFREAT	=	0				!number of fixed atoms?
	DELTA		=	delta_t				!delta_t (x-plor is double, charmm is real)
	NINEZ(:)	=	0				!buffer zeros
	NTITLE		=	2				!number of 80-character strings in title
	TITLE(1)	=	'  Simulation record file '	!
	TITLE(2)        =	'   Written in parallel   '	!
	NATOM		=	extralloc			!number of particles

	allocate(Xbuf(NSET*extralloc))
	allocate(Ybuf(NSET*extralloc))
	allocate(Zbuf(NSET*extralloc))

	!Read Solid molecule position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file="results/vmd_halo_temp.dcd",access='direct',recl=1)

	do i=1,NSET
	do n=1,extralloc
		read(17,rec=(i-1)*nd*extralloc+n) Xbuf(n+extralloc*(i-1))
		read(17,rec=(i-1)*nd*extralloc+n+extralloc) Ybuf(n+extralloc*(i-1))
		read(17,rec=(i-1)*nd*extralloc+n+2*extralloc) Zbuf(n+extralloc*(i-1))
	enddo
	enddo

	close(17,status='delete')

	!Open binary .dcd file and write header information	
	open(unit=3, file="results/vmd_halo_out.dcd",status='replace', form="unformatted")
	
	write(3) HDR, NSET, ISTRT, NSAVC, FIVEZ, NATOMNFREAT, DELTA, NINEZ
	write(3) NTITLE, TITLE(1), TITLE(2)
	write(3) NATOM

	do i=1,NSET
		write(3) Xbuf((i-1)*extralloc+1:i*extralloc)
		write(3) Ybuf((i-1)*extralloc+1:i*extralloc)
		write(3) Zbuf((i-1)*extralloc+1:i*extralloc)
	enddo

	close(3,status='keep')

end subroutine reformat_dcd_halo
