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
		if (potential_flag.eq.0) then	
			print*, 'Results from last iteration of simulation'
			print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.4)', &
			Nsteps,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy,';',totenergy,';',pressure
		else if (potential_flag.eq.1) then
			print*, 'Results from last iteration of simulation'
			print '(1x,i8,a,f15.4,a,f15.4,a,f10.4,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.5,a,f10.4)', &
			Nsteps,';',vsum,';', v2sum,';', temperature,';', &
			kinenergy,';',potenergy_LJ,';',potenergy_FENE,';',potenergy,';',totenergy,';',pressure
		end if
	end if
	
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
		if(vmd_outflag .eq. 3) then
			call reformat_dcd
			call reformat_dcd_halo
		end if
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

	if (macro_outflag.eq.2) close(10,status='keep') !Keep macroscopic_properties
		
	if (potential_flag.eq.1) call build_psf

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

        real time_start, time_end
 
	print*, 'Generating final VMD.dcd ouput file - for large systems this may take some time ...'

        call cpu_time(time_start)

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

        write(0,*)'globalnp ', globalnp

	!Read position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file=trim(prefix_dir)//"results/vmd_temp.dcd",access='stream')

	do i=1,NSET
		read(17) Xbuf(globalnp*(i-1)+1:globalnp*i)
		read(17) Ybuf(globalnp*(i-1)+1:globalnp*i)
		read(17) Zbuf(globalnp*(i-1)+1:globalnp*i)
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
	enddo

	close(3,status='keep')

        call cpu_time(time_end)

        print '(a,g7.0,a)', 'Generated final VMD.dcd ouput file in', time_end - time_start, ' seconds'



end subroutine reformat_dcd


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
	open (unit=17, file=trim(prefix_dir)//"results/vmd_solid_temp.dcd",access='direct',recl=1)

	do i=1,NSET
	do n=1,globalnp
		read(17,rec=(i-1)*nd*globalnp+n) Xbuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+globalnp) Ybuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+2*globalnp) Zbuf(n+globalnp*(i-1))
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
	enddo

	close(3,status='keep')

	!Read liquid molecule position information from file
	!RECORD LENGTH IS 1 WHICH IN FORTRAN IS A 4 BYTE BLOCKS (REAL, INT BUT NOT DP) 	
	open (unit=17, file=trim(prefix_dir)//"results/vmd_liquid_temp.dcd",access='direct',recl=1)

	do i=1,NSET
	do n=1,globalnp
		read(17,rec=(i-1)*nd*globalnp+n) Xbuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+globalnp) Ybuf(n+globalnp*(i-1))
		read(17,rec=(i-1)*nd*globalnp+n+2*globalnp) Zbuf(n+globalnp*(i-1))
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
	open (unit=17, file=trim(prefix_dir)//"results/vmd_halo_temp.dcd",access='direct',recl=1)

	do i=1,NSET
	do n=1,extralloc
		read(17,rec=(i-1)*nd*extralloc+n) Xbuf(n+extralloc*(i-1))
		read(17,rec=(i-1)*nd*extralloc+n+extralloc) Ybuf(n+extralloc*(i-1))
		read(17,rec=(i-1)*nd*extralloc+n+2*extralloc) Zbuf(n+extralloc*(i-1))
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
	enddo

	close(3,status='keep')

end subroutine reformat_dcd_halo

!==============================================================================================================

subroutine build_psf
	use module_final_record
	use polymer_info_MD
	implicit none

	integer :: i,j,item,nchains
	integer :: NTITLE, NATOM, NBONDS
	integer	:: write_items
	integer, allocatable, dimension(:,:) :: bonds	
	integer, allocatable, dimension(:)	::	res_ID
	character(len=4), allocatable, dimension(:) :: seg_name, res_name, atom_name, atom_type
	real, allocatable, dimension(:)	:: charge, mass
	
	nchains = ceiling(np/dble(chain_length))		
	
	NTITLE = 1												! How many 'REMARKS' lines you want
	NATOM  = globalnp										! Determine total number of atoms
	NBONDS = (chain_length-1)*nchains						! Determine total number of bonds

	print*, 'Generating polymer topology file polymer_topol.psf'
	open(unit=1, file='results/polymer_topol.psf', status='replace', form='formatted')
	
	! Header
	write(1,'(a3)') 'PSF'
	write(1,'(a1)')
	write(1,'(i8,a)') NTITLE, ' !NTITLE'
	write(1,'(a9,a)') "REMARKS ","David's first attempt at writing a *.psf file in FORTRAN"
		
	! Atoms
	write(1,'(a1)')
	write(1,'(i8,a)') NATOM, ' !NATOM'

	allocate(seg_name(NATOM))								! Determine segment names for each atom
	allocate(res_ID(NATOM))									! Determine molecule ID for each atom
	allocate(res_name(NATOM))								! Determine name for each molecule
	allocate(atom_name(NATOM))								! Determine name for each atom
	allocate(atom_type(NATOM))								! Determine type for each atom
	allocate(charge(NATOM))									! Determine charge for each atom
	allocate(mass(NATOM))									! Determine mass of each atom
	seg_name(:)='C'
	res_ID(:) =	polyinfo_mol(:)%chainID 
	res_name(:) = 'CBN'
	atom_name(:)='C'
	atom_type(:)='C'
	charge(:)=0.00000
	mass(:)=1.00794

	do i=1,NATOM
		write(1,'(i8,3a,i4,6a,2f10.5,i1)') i,' ',seg_name(i),' ',res_ID(i),'&
                 ',res_name(i),' ',atom_name(i),' ',atom_type(i),charge(i),mass(i),0 
	end do

	! Bonds
	write(1,'(a1)')
	write(1,'(i8,a)') NBONDS, ' !NBONDS'
	
	write_items = 4*ceiling(NBONDS/4.0)
	allocate(bonds(write_items,2))
	bonds(:,:)=0

	item=1
	do i=1,globalnp
		do j=1,i-1
			if (polyinfo_mol(i)%chainID.eq.polyinfo_mol(j)%chainID) then
				if (abs(polyinfo_mol(i)%subchainID-polyinfo_mol(j)%subchainID).eq.1) then
					bonds(item,1) = i
					bonds(item,2) = j
					item=item+1
				end if
			end if
		end do
	end do

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

end subroutine build_psf

