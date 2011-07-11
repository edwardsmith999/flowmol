!======================================================================
!	        	Parallel Read Write Subroutines               =
!======================================================================

module module_parallel_io
	use mpi
	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use calculated_properties_MD

	integer		:: restartfileid, fileid !File name used for parallel i/o

end module

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer				:: procdisp
	integer				:: i, datasize
	real,dimension(np)		:: Xbuf, Ybuf, Zbuf
	integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

	!Determine size of real datatype
 	call MPI_type_size(MPI_real,datasize,ierr)

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre
	Xbuf(:) = r(:,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
	Ybuf(:) = r(:,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
	Zbuf(:) = r(:,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

	procdisp = 0
	!Obtain displacement of each processor using all other procs' np
	do i = 1, irank -1
		procdisp = procdisp + procnp(i)*datasize
	enddo

	!Open file on all processors
	call MPI_FILE_OPEN(MPI_COMM_WORLD, 'results/vmd_temp.dcd', & 
		MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		MPI_INFO_NULL, fileid, ierr)

	!-------------Write X coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp				  	!Processor location

	!print*, irank, 'x disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Xbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Y coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp &				  	!Processor location
		+ globalnp * datasize			  	!Y Coordinate location

	!print*, irank, 'y disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Ybuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Z coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp &				  	!Processor location
		+ 2 * globalnp * datasize		  	!Z Coordinate location

	!print*, irank, 'z disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Zbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!call MPI_FILE_GET_SIZE(fileid, resultsize, ierr)
	!print*, 'Vmd result filesize', resultsize

	!Close file on all processors
	call MPI_FILE_CLOSE(fileid, ierr) 

end subroutine parallel_io_vmd

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd_sl
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer				:: procdisp
	integer				:: i,n, datasize
	real,dimension(np)		:: Xbuf, Ybuf, Zbuf
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
		case(0)
			Xbuf(:) = r(:,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
			Ybuf(:) = r(:,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
			Zbuf(:) = r(:,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
		case(1:)
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
	call MPI_FILE_OPEN(MPI_COMM_WORLD, 'results/vmd_liquid_temp.dcd', & 
		MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		MPI_INFO_NULL, fileid, ierr)

	!-------------Write X coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp				  	!Processor location

	!print*, irank, 'x disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Xbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Y coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp &				  	!Processor location
		+ globalnp * datasize			  	!Y Coordinate location

	!print*, irank, 'y disp', disp
	
	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
		
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Ybuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Z coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp &				  	!Processor location
		+ 2 * globalnp * datasize		  	!Z Coordinate location

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
		case(0)
			Xbuf(n) = -halfdomain(1)
			Ybuf(n) = -halfdomain(2)
			Zbuf(n) = -halfdomain(3)
		case(1:)
			Xbuf(n) = r(n,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
			Ybuf(n) = r(n,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
			Zbuf(n) = r(n,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
		end select
	enddo

	!================================
	!    Write Solid Molecules	=
	!================================

	!Open file on all processors
	call MPI_FILE_OPEN(MPI_COMM_WORLD, 'results/vmd_solid_temp.dcd', & 
		MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		MPI_INFO_NULL, fileid, ierr)

	!-------------Write X coordinates--------------------
		
	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp				  	!Processor location

	!print*, irank, 'x disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Xbuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Y coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp &				  	!Processor location
		+ globalnp * datasize			  	!Y Coordinate location

	!print*, irank, 'y disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, Ybuf, np, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!-------------Write Z coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp &				  	!Processor location
		+ 2 * globalnp * datasize		  	!Z Coordinate location

	!print*, irank, 'z disp', disp

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

	integer				:: procdisp
	integer				:: i, datasize
	real,dimension(np*3)		:: buf
	integer(kind=MPI_OFFSET_KIND)   :: disp!, resultsize

	!Build array of number of particles on neighbouring
	!processe's subdomain on current proccess
	call globalGathernp

	!Determine size of real datatype
 	call MPI_type_size(MPI_real,datasize,ierr)

	!Load buffers with single precision r and adjust according
	!to processor topology with r = 0 at centre

	!buf(1:np) = r(:,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
	!buf((np+1):(2*np)) = r(:,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
	!buf((2*np+1):(3*np)) = r(:,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

	do i=1,np
		buf(3*(i-1)+1) = r(i,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
		buf(3*(i-1)+2) = r(i,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
		buf(3*(i-1)+3) = r(i,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
	enddo

	procdisp = 0
	!Obtain displacement of each processor using all other procs' np
	do i = 1, irank -1
		procdisp = procdisp + procnp(i)*datasize
	enddo

	!Open file on all processors
	call MPI_FILE_OPEN(MPI_COMM_WORLD, 'results/vmd_temp.dcd', & 
		MPI_MODE_RDWR + MPI_MODE_CREATE, & 
		MPI_INFO_NULL, fileid, ierr)

	!-------------Write X coordinates--------------------

	!Obtain location to write in file
	disp =(iter/real((tplot),kind(0.d0))-1) * nd * globalnp * datasize & !Current iteration
		+ procdisp				  	!Processor location

	!print*, irank, 'x disp', disp

	call MPI_FILE_SET_VIEW(fileid, disp, MPI_REAL, & 
 		MPI_REAL, 'native', MPI_INFO_NULL, ierr)
	
	!Write information to file
	call MPI_FILE_WRITE_ALL(fileid, buf, np*3, MPI_REAL, & 
		MPI_STATUS_IGNORE, ierr) 

	!Close file on all processors
	call MPI_FILE_CLOSE(fileid, ierr) 

end subroutine parallel_io_vmd_optimised


!------------------------------------------------------------------------

subroutine parallel_io_vmd_halo

	print*, 'VMD halo output not possible for parallel code'

end subroutine parallel_io_vmd_halo

!------------------------------------------------------------------------
!Write positions and velocities of molecules to a file for restart

subroutine parallel_io_final_state
	use module_parallel_io
	implicit none
	!include 'mpif.h'

	integer				:: n, i
	integer				:: procdisp	!Processor's displacement
	integer				:: dp_datasize, filesize
	integer(kind=MPI_OFFSET_KIND)   :: disp
	double precision, dimension(nd)	:: Xwrite	!Temporary variable used in write

	!Rebuild simulation before recording final state
	call linklist_deallocateall	   !Deallocate all linklist components
	call sendmols			   !Exchange particles between processors
	call assign_to_cell	  	   !Re-build linklist every timestep
	call messenger_updateborders	   !Update borders between processors
	call assign_to_halocell		   !Re-build linklist
	call assign_to_neighbourlist	   !Setup neighbourlist

	!Build array of number of particles on neighbouring
	!process' subdomains on current proccess
	call globalGathernp

	!Determine size of datatypes
  	call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

	!Adjust r according to actual location for storage according
	!to processor topology with r = 0 at centre
	r(:,1) = r(:,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
	r(:,2) = r(:,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
	r(:,3) = r(:,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)

	!Obtain displacement of each processor using all other procs' np
	!with 3 position and 3 velocity components for each molecule
	procdisp = 0
	do i = 1, irank -1
		procdisp = procdisp + 2*nd*procnp(i)*dp_datasize
	enddo

	!Remove previous final state file
	call MPI_FILE_DELETE('results/final_state', MPI_INFO_NULL, ierr)

	!Open file on all processors
	call MPI_FILE_OPEN(MPI_COMM_WORLD, 'results/final_state', & 
			MPI_MODE_RDWR + MPI_MODE_CREATE, & 
			MPI_INFO_NULL, restartfileid, ierr)

	!-------------Write coordinates--------------------

	!Obtain location to write in file
	disp =  procdisp

	!Set each processor to that location and write particlewise
	call MPI_FILE_SET_VIEW(restartfileid, disp, MPI_double_precision, & 
 		MPI_double_precision, 'native', MPI_INFO_NULL, ierr)

	do n = 1, np
		Xwrite = r(n,:) !Load into temp in case r dimensions are non contiguous
		call MPI_FILE_WRITE(restartfileid, Xwrite, nd, MPI_double_precision, & 
					MPI_STATUS_IGNORE, ierr) 
		Xwrite = v(n,:) !Load into temp in case v dimensions are non contiguous
		call MPI_FILE_WRITE(restartfileid, Xwrite, nd, MPI_double_precision, & 
					MPI_STATUS_IGNORE, ierr) 
	enddo

	!----------------Write header------------------------
	!Written at the end for performance and simplicity reasons 
	!(See Gropp, lusk & Thakur Using MPI-2)

	!Obtain location to write in file
	disp =  2*nd*globalnp*dp_datasize

	!Obtain location of end of file
	call MPI_File_get_size(restartfileid,filesize,ierr)

	!File point incremented by write calls of both int and DP 
	!so etype and fileview set in bytes
	call MPI_FILE_SET_VIEW(restartfileid, disp, MPI_byte, & 
 		MPI_byte, 'native', MPI_INFO_NULL, ierr)

	!Write with one processor only
	if (irank .eq. iroot) then

		!Double precision data
		call MPI_FILE_WRITE(restartfileid,   density,      1, MPI_double_precision, & 
					MPI_STATUS_IGNORE, ierr) !Density of system
		call MPI_FILE_WRITE(restartfileid,   rcutoff,      1, MPI_double_precision, & 
					MPI_STATUS_IGNORE, ierr) !Cut off distance for particle interaction
		call MPI_FILE_WRITE(restartfileid,inputtemperature,1,MPI_double_precision, & 
					MPI_STATUS_IGNORE, ierr) !Define initial temperature
		call MPI_FILE_WRITE(restartfileid, 	delta_t,    1, MPI_double_precision, & 
					MPI_STATUS_IGNORE, ierr) !Size of time step
		call MPI_FILE_WRITE(restartfileid,   elapsedtime,   1, MPI_double_precision, & 
					MPI_STATUS_IGNORE, ierr) !Total elapsed time of all restarted simulations

		!Integer Data
		call MPI_FILE_WRITE(restartfileid,globalnp, 	   1, MPI_integer, & 
					MPI_STATUS_IGNORE, ierr) !Number of particles
		call MPI_FILE_WRITE(restartfileid,initialnunits(:),3,MPI_integer, & 
					MPI_STATUS_IGNORE, ierr) !x, y & z dimension split into number of cells
		call MPI_FILE_WRITE(restartfileid, 	Nsteps,     1, MPI_integer, & 
					MPI_STATUS_IGNORE, ierr) !Number of computational steps
		call MPI_FILE_WRITE(restartfileid, 	tplot,      1, MPI_integer, & 
					MPI_STATUS_IGNORE, ierr) !Frequency at which to record results
		call MPI_FILE_WRITE(restartfileid, 	seed(1),    1, MPI_integer, & 
					MPI_STATUS_IGNORE, ierr) !Random number seed value 1
		call MPI_FILE_WRITE(restartfileid, 	seed(2),    1, MPI_integer, & 
					MPI_STATUS_IGNORE, ierr) !Random number seed value 2

	endif

	!Close file on all processors
	call MPI_FILE_CLOSE(restartfileid, ierr) 

end subroutine parallel_io_final_state

!------------------------------------------------------------------------------
! Set up inputs on every processor, based on the final state of a previous simulation

subroutine setup_restart_inputs
	use module_parallel_io
	implicit none
	!include 'mpif.h'

	integer				:: n, k
	integer 			:: extrasteps
	integer				:: filesize,int_filesize,dp_filesize
	integer 			:: checkint
	integer(kind=MPI_OFFSET_KIND)   :: disp
	double precision 		:: checkdp

	!Allocate random number seed
	call random_seed
	call random_seed(size=n)
	allocate(seed(n))

	!Open on a single process and broadcast
	if (irank .eq. iroot) then

		!Call function library to get file size
		call get_file_size('final_state',filesize)

		!File size is in bytes and integer fortran records are blocks of 4 bytes
		int_filesize = filesize/4

		!Open file to read integers
		open(13,file='final_state', form="unformatted", access="direct",recl=1)

		read(13,rec=int_filesize-7) globalnp	    !Number of particles
		read(13,rec=int_filesize-6) initialnunits(1)  !x dimension split into number of cells
		read(13,rec=int_filesize-5) initialnunits(2)  !y dimension box split into number of cells
		read(13,rec=int_filesize-4) initialnunits(3)  !z dimension box split into number of cells
		read(13,rec=int_filesize-3) Nsteps  	    !Number of elapsed computational steps

		close(13,status='keep')

		!File size is in bytes and fortran double precision records are blocks of 8 bytes
		dp_filesize = filesize/8

		!Reopen file to read doubles
		open(13,file='final_state', form="unformatted", access="direct",recl=2)

		read(13,rec=dp_filesize-8) density		!Density of system
		read(13,rec=dp_filesize-7) rcutoff		!Cut off distance for particle interaction
		read(13,rec=dp_filesize-6) inputtemperature	!Define initial temperature
		read(13,rec=dp_filesize-4) elapsedtime   	!Elapsed simulation time to date

		close(13,status='keep')

		!Check if values from input file are different and alert user - all processors have
		!the same values so only need to check on one processor
		open(11,file='input')
	
		read(11,* ) checkdp          !Density of system
		if (checkdp .ne. density) print*, 'Discrepancy between system density', &
		'in input & restart file - restart file will be used'
		read(11,* ) checkdp          !Cut off distance for particle interaction
		if (checkdp .ne. rcutoff) print*, 'Discrepancy between cut off radius', &
		'in input & restart file - restart file will be used'
		read(11,* ) checkdp	     !Define initial temperature
		if (checkdp .ne. inputtemperature) print*, 'Discrepancy between initial temperature', &
		'in input & restart file - restart file will be used'
		read(11,* ) checkint	     !x dimension split into number of cells
		if (checkint .ne. initialnunits(1)) print*, 'Discrepancy between x domain size', &
		'in input & restart file - restart file will be used'
		read(11,* ) checkint         !y dimension box split into number of cells
		if (checkint .ne. initialnunits(2)) print*, 'Discrepancy between y domain size', &
		'in input & restart file - restart file will be used'
		read(11,* ) checkint	     !z dimension box split into number of cells
		if (nd == 3) then
		if (checkint .ne. initialnunits(3)) print*, 'Discrepancy between z domain size', &
		'in input & restart file - restart file will be used'
		endif

		!Get number of extra steps, timestep and plot frequency from input file
	
		read(11,* ) extrasteps       !Number of computational steps
		read(11,* ) delta_t          !Size of time step
		read(11,* ) tplot            !Frequency at which to record results
		read(11,* ) delta_rneighbr   !Extra distance used for neighbour cell
		read(11,* ) seed(1)	     !Random number seed value 1
		read(11,* ) seed(2)	     !Random number seed value 2

		!Flag to determine if output is switched on
		read(11,* ) vmd_outflag
		read(11,* ) macro_outflag	
		read(11,* ) slice_outflag
		read(11,* ) Nslice_ave
		read(11,* ) pressure_outflag
		read(11,* ) viscsample
		read(11,* ) Nvisc_ave

		close(11,status='keep')      !Close input file

	endif

	!Broadcast data read by root to all other processors
	call MPI_BCAST(globalnp,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(density,1,MPI_double_precision,iroot-1,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(rcutoff,1,MPI_double_precision,iroot-1,MPI_COMM_WORLD,ierr) 
	rcutoff2= rcutoff**2         !Useful definition to save computational time
	call MPI_BCAST(inputtemperature,1,MPI_double_precision,iroot-1,MPI_COMM_WORLD,ierr) 
	call MPI_BCAST(initialnunits(1),1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(initialnunits(2),1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(initialnunits(3),1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(Nsteps,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(extrasteps,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(delta_t,1,MPI_double_precision,iroot-1,MPI_COMM_WORLD,ierr) 
	call MPI_BCAST(tplot,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(delta_rneighbr,1,MPI_double_precision,iroot-1,MPI_COMM_WORLD,ierr) 
	call MPI_BCAST(vmd_outflag,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(macro_outflag,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr) 
	call MPI_BCAST(slice_outflag,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr) 
	call MPI_BCAST(Nslice_ave,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr) 
	call MPI_BCAST(pressure_outflag,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(viscsample,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)  
	call MPI_BCAST(Nvisc_ave,1,MPI_integer,iroot-1,MPI_COMM_WORLD,ierr)

	elapsedtime = elapsedtime + delta_t*extrasteps !Set elapsed time to end of simualtion
	initialstep = Nsteps         !Set plot count to final plot of last
	Nsteps = Nsteps + extrasteps !Establish final iteration step based on previous

	!print*,irank,globalnp,density,rcutoff,Nsteps,elapsedtime,inputtemperature,initialnunits, &
	!	elapsedtime,extrasteps,delta_t,tplot,delta_rneighbr
	!call MPI_BARRIER(MPI_COMM_WORLD, ierr)

end subroutine setup_restart_inputs

!------------------------------------------------------------------------------

subroutine setup_restart_microstate
	use module_parallel_io
	implicit none
	!include 'mpif.h'

	integer 			:: n, nl, ixyz, procassign
	integer				:: dp_datasize
	integer(kind=MPI_OFFSET_KIND)   :: disp
	double precision, dimension (nd):: rc !Temporary variable

	!Determine size of datatypes
  	call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

	!Open restart file on all processor
	call MPI_FILE_OPEN(MPI_COMM_WORLD, 'final_state', & 
		MPI_MODE_RDONLY , MPI_INFO_NULL, restartfileid, ierr)

	nl = 0		!Reset local molecules count nl

	!read positions
	do n=1,globalnp

		!---------- For all molecule positions ------------

		!Move through location of position co-ordinates
		disp =  2 * (n-1) * nd * dp_datasize				 

		!Set each processor to that location and write particlewise
		call MPI_FILE_SET_VIEW(restartfileid, disp, MPI_double_precision, & 
 					MPI_double_precision, 'native', MPI_INFO_NULL, ierr)
		call MPI_FILE_READ_ALL(restartfileid, rc(:), nd, MPI_double_precision, & 
					MPI_STATUS_IGNORE, ierr) !Read position from file

		!Use integer division to determine which processor to assign molecule to
		procassign = ceiling((rc(1)+globaldomain(1)/2.d0)/domain(1))
		if (procassign .ne. iblock) cycle
		procassign = ceiling((rc(2)+globaldomain(2)/2.d0)/domain(2))
		if (procassign .ne. jblock) cycle
		procassign = ceiling((rc(3)+globaldomain(3)/2.d0)/domain(3))
		if (procassign .ne. kblock) cycle

		!print*, irank, n

		!If molecules is in the domain then add to processor's total
		nl = nl + 1 !Local molecule count

		!Correct to local coordinates
		r(nl,1) = rc(1)-domain(1)*(iblock-1)+halfdomain(1)*(npx-1)
		r(nl,2) = rc(2)-domain(2)*(jblock-1)+halfdomain(2)*(npy-1)
		r(nl,3) = rc(3)-domain(3)*(kblock-1)+halfdomain(3)*(npz-1)

		!print'(2(i8,a),3(f10.5,a))', iblock, ';', n, ';', r(nl,1), ';', r(nl,2), ';', r(nl,3), ';'

		!Read corresponding velocities
		call MPI_FILE_READ(restartfileid, v(nl,:), nd, MPI_double_precision, & 
					MPI_STATUS_IGNORE, ierr) !Read position from file

	enddo

	np = nl	!Correct local number of particles on processor

	!Close and remove final state file from directory
	call MPI_FILE_CLOSE(restartfileid, ierr) 
	call MPI_FILE_DELETE('final_state', MPI_INFO_NULL, ierr)

end subroutine setup_restart_microstate


!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!SERIAL CODE NOT PARALLELISED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine velocity_average_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	open(10,file='results/v_average') !Open velocity slice file

end subroutine velocity_average_header

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!SERIAL CODE NOT PARALLELISED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine velocity_average_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer	:: ibin, jbin, kbin

	write(10,'(4(a,i8))')	& 
	 	'; ', iter , 	&
		'; ', nbins(1), 	&
		'; ', nbins(2) , 	&
		'; ', nbins(3)

	do kbin = 1,nbins(3)
	do ibin = 1,nbins(1)
	do jbin = 1,nbins(2)
		write(10,'( 3(a, f10.5), a , i18)')		& 
		 	'; ', meanvelbin(ibin,jbin,kbin,1), 	&
			'; ', meanvelbin(ibin,jbin,kbin,2), 	&
			'; ', meanvelbin(ibin,jbin,kbin,3), 	&
			'; ',countvelbin(ibin,jbin,kbin  )
	enddo
	enddo
	enddo

end subroutine velocity_average_io

!------------------------------------------------------------------------
!Write velocity slice information to a file

subroutine parallel_slice_header
	use module_parallel_io
	implicit none

	if (irank .eq. iroot) then
		open(18,file='results/vslice') !Open velocity slice file
		write(18,'(i10,a,f10.5)') globalnp,';', globaldomain(2)
		write(18,'(i10,a,i10)') globalnbins(1),';', Nsteps
		write(18,'(i10,a,i10)') tplot,';', Nslice_ave
	endif

end subroutine parallel_slice_header


subroutine parallel_slice_io
	use messenger
	use module_parallel_io
	implicit none

	integer						:: ibin, jbin, kbin
	integer, dimension(globalnbins(1))			:: globalcountvel
	double precision, dimension(globalnbins(1))	:: globalmeanvel

	call MPI_gather(meanvel,nbins(1),MPI_double_precision,  & 
			globalmeanvel,nbins(1),MPI_double_precision, & 
			iroot-1,icomm_xyz(2),ierr)

	call MPI_gather(countvel,nbins(1),MPI_integer, & 
			globalcountvel,nbins(1),MPI_integer, &
			iroot-1,icomm_xyz(2),ierr)

	if (irank .eq. iroot) then

		ibin = 1	
		kbin = 1
		do jbin =1 , globalnbins(1)
			write(18,'(f20.5,a,i10)') globalmeanvel(jbin) & 
					    ,';', globalcountvel(jbin)
		enddo
	endif

end subroutine parallel_slice_io


!------------------------------------------------------------------------
!A large scale routine with each proc writing its own bins in binary
!Write velocity slice information to a file

subroutine parallel_slice_header_large_scale
	use module_parallel_io
	implicit none

	integer		:: slicefileid

	!Open file on all processors
	call MPI_FILE_OPEN(MPI_COMM_WORLD, 'results/vslice', & 
			MPI_MODE_RDWR + MPI_MODE_CREATE, & 
			MPI_INFO_NULL, slicefileid, ierr)

	!Root processes writes file information
	if (irank .eq. iroot) then
		call MPI_FILE_WRITE(slicefileid,globalnp, 1, MPI_integer, & 
				MPI_STATUS_IGNORE, ierr) !Number of particles
		call MPI_FILE_WRITE(slicefileid,globaldomain(2), 1, MPI_double_precision, & 
				MPI_STATUS_IGNORE, ierr) !Domain size
		call MPI_FILE_WRITE(slicefileid,globalnbins(1), 1, MPI_integer, & 
				MPI_STATUS_IGNORE, ierr) !Number of bins
		call MPI_FILE_WRITE(slicefileid,Nsteps, 1, MPI_integer, & 
				MPI_STATUS_IGNORE, ierr) !Number of simulation steps
		call MPI_FILE_WRITE(slicefileid,tplot, 1, MPI_integer, & 
				MPI_STATUS_IGNORE, ierr) !Frequency of plots
		call MPI_FILE_WRITE(slicefileid,Nslice_ave, 1, MPI_integer, & 
				MPI_STATUS_IGNORE, ierr) !Number of results used in slice average
	endif

	call MPI_FILE_CLOSE(slicefileid, ierr) 

end subroutine parallel_slice_header_large_scale



subroutine parallel_slice_io_large_scale
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer					:: slicefileid, int_datasize, dp_datasize
	integer(kind=MPI_OFFSET_KIND)   	:: disp!, resultsize

	!Determine size of datatypes
	call MPI_type_size(MPI_Integer,int_datasize,ierr)
	call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

	call MPI_FILE_OPEN(MPI_COMM_WORLD, 'results/vslice', & 
			MPI_MODE_WRONLY , & 
			MPI_INFO_NULL, slicefileid, ierr)

	!WRITE BINNED VELOCITY SUMMATION

	!Obtain displacement of current record
	disp =(iter/real((tplot*Nslice_ave),kind(0.d0))-1) *       &	!Current iteration
	      globalnbins(1) *(int_datasize+dp_datasize) & !Current iteration	
	      + nbins(1)*(dp_datasize)*(jblock-1)		!Processor location
	!Set current processes view
	call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_double_precision, & 
 			MPI_double_precision, 'native', MPI_INFO_NULL, ierr)
	!Write data
	call MPI_FILE_WRITE_ALL(slicefileid,meanvel(:),nbins(1), MPI_double_precision, & 
			MPI_STATUS_IGNORE, ierr) !Velocity bins

	!WRITE BINNED MOLECULAR COUNT

	!Obtain displacement of current record
	disp =(iter/real((tplot*Nslice_ave),kind(0.d0))-1) *		& !Current iteration
	      globalnbins(1) * (int_datasize+dp_datasize)  & !Current iteration	
		+ nbins(1)*(int_datasize)*(jblock-1) 	& !Processor location
		+ globalnbins(1)* dp_datasize 		  !After binned velocity

	!Set current processes view
	call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_Integer, & 
				MPI_Integer, 'native', MPI_INFO_NULL, ierr)

	!Write data
	call MPI_FILE_WRITE_ALL(slicefileid,countvel(:),nbins(1), MPI_Integer, & 
			MPI_STATUS_IGNORE, ierr) !molecular tally bins

	call MPI_FILE_CLOSE(slicefileid, ierr) 

end subroutine parallel_slice_io_large_scale

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!SERIAL CODE NOT PARALLELISED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine stress_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	open(2,file='results/pressure_visc') !Open pressure field file

end subroutine stress_header

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!SERIAL CODE NOT PARALLELISED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine virial_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	write(2,'(a,i8)') 'Record of stress tensor at iter ', iter
	write(2,'(3(a, f18.15))')	   '; ', Pxy(1,1), &
					   '; ', Pxy(1,2), &
					   '; ', Pxy(1,3)
	write(2,'(3(a, f18.15))') 	   '; ', Pxy(2,1), &
					   '; ', Pxy(2,2), &
					   '; ', Pxy(2,3)
	write(2,'(3(a, f18.15))') 	   '; ', Pxy(3,1), &
					   '; ', Pxy(3,2), &
					   '; ', Pxy(3,3)

end subroutine virial_stress_io

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!SERIAL CODE NOT PARALLELISED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine VA_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer	:: ibin, jbin, kbin

	write(2,'(3(a,i8),6(a,f10.5))')	'; ', nbins(1),'; ',nbins(2),'; ',nbins(3),	&
					'; ', Pxy(1,1),'; ',Pxy(2,2),'; ',Pxy(3,3),	& 
					'; ', Pxy(1,2),'; ',Pxy(1,3),'; ',Pxy(2,3)
	do kbin = 1,nbins(3)
	do ibin = 1,nbins(1)
	do jbin = 1,nbins(2)
		write(2,'(9(a, f10.5))')   '; ', Pxybin(ibin,jbin,kbin,1,1), &
					   '; ', Pxybin(ibin,jbin,kbin,1,2), &
					   '; ', Pxybin(ibin,jbin,kbin,1,3), &
					   '; ', Pxybin(ibin,jbin,kbin,2,1), &
					   '; ', Pxybin(ibin,jbin,kbin,2,2), &
					   '; ', Pxybin(ibin,jbin,kbin,2,3), &
					   '; ', Pxybin(ibin,jbin,kbin,3,1), &
					   '; ', Pxybin(ibin,jbin,kbin,3,2), &
					   '; ', Pxybin(ibin,jbin,kbin,3,3)
	enddo
	enddo
	enddo

end subroutine VA_stress_io

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!SERIAL CODE NOT PARALLELISED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MOP_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer	:: ibin, jbin, kbin, n

	write(2,'(i8)') iter
	do n = 1, nplanes
		write(2,'(i8,3f18.8)'),n, Pxy_plane(:,n)
	enddo

end subroutine MOP_stress_io
!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!SERIAL CODE NOT PARALLELISED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine control_volume_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	open(4,file='results/control_volume') !Open control volume file

end subroutine control_volume_header

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!SERIAL CODE NOT PARALLELISED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Control_Volume_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: ibin, jbin, kbin
	double precision,dimension(3,6)	:: Pxyplane

	write(4,'(3(a,i8),6(a,f10.5))')	'; ', nbins(1),'; ',nbins(2),'; ',nbins(3),	&
					'; ', Pxy(1,1),'; ',Pxy(2,2),'; ',Pxy(3,3),	& 
					'; ', Pxy(1,2),'; ',Pxy(1,3),'; ',Pxy(2,3)

end subroutine Control_Volume_io
!---------------------------------------------------------------------------------
