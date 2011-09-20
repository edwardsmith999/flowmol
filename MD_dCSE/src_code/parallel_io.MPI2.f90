!======================================================================
!	        	Parallel Read Write Subroutines               

! Parallel i/o

! --- INPUT ROUTINES ---
! setup_inputs			Read input file
! setup_restart_inputs		Read restart files and input file
! setup_restart_microstate	Read Initial configuration from restart file

! --- OUTPUT ROUTINES ---
! parallel_io_final_state 	Write final configuration and simulation details for restart
! parallel_io_vmd		Output for VMD
! parallel_io_vmd_sl		Output for VMD with seperate solid/lquid regions
! parallel_io_vmd_halo		Output Halos
! CV AVERAGING
! mass_averaging		slice/CV
! cumulative_mass		slice/CV
! mass_snapshot			CV
! momentum_averaging		slice/CV
! cumulative_momentum		slice/CV
! momentum_snapshot		CV
! pressure_averaging		domain/CV
! cumulative_pressure		domain/CV
! FLUX AVERAGING
! mass_flux_averaging		CV_surface
! cumulative_mass_flux		CV_surface
! momentum_flux_averaging	plane(MOP)/CV_surface
! cumulative_momentum_flux	plane(MOP)/CV_surface
! surface_pressure		plane(MOP)/CV_surface
! cumulative_pressure		plane(MOP)/CV_surface
! 
!======================================================================

module module_parallel_io
	use mpi
	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use calculated_properties_MD
	use messenger, only : MD_COMM

	integer		:: restartfileid, fileid !File name used for parallel i/o

end module

!======================================================================
!	        		INPUTS			              =
!======================================================================

!=============================================================================
! Input values used to set up the simulation such as number of dimensions and
! number of particles
!
!-----------------------------------------------------------------------------

subroutine setup_inputs
use module_parallel_io
implicit none
  
	integer :: k, n

	call random_seed
	call random_seed(size=n)
	allocate(seed(n))

	open(1,file=trim(file_dir)//'input')

	!Input physical co-efficients
	read(1,* ) density          !Density of system
	read(1,* ) rcutoff          !Cut off distance for particle interaction
	rcutoff2= rcutoff**2         !Useful definition to save computational time
	read(1,* ) inputtemperature !Define initial temperature
	read(1,* ) initialnunits(1) !x dimension split into number of cells
	read(1,* ) initialnunits(2) !y dimension box split into number of cells

	if (nd .eq. 3) then
	read(1,* ) initialnunits(3) !z dimension box split into number of cells
	else
	read(1,* ) k		     !Read into dummy variable as value not used
	endif
	!Input computational co-efficients
	read(1,* ) Nsteps           !Number of computational steps
	read(1,* ) delta_t          !Size of time step
	read(1,* ) tplot            !Frequency at which to record results
	initialstep = 0   	     !Set initial step to one to start
	read(1,* ) delta_rneighbr   !Extra distance used for neighbour cell
	read(1,* ) seed(1)	     !Random number seed value 1
	read(1,* ) seed(2)	     !Random number seed value 2

	!Flag to determine if output is switched on
	read(1,* ) vmd_outflag
	read(1,* ) macro_outflag
	read(1,* ) mass_outflag
	read(1,* ) Nmass_ave
	read(1,* ) velocity_outflag
	read(1,* ) Nvel_ave
	read(1,* ) pressure_outflag
	read(1,* ) Nstress_ave
	read(1,* ) viscosity_outflag
	read(1,* ) Nvisc_ave
	read(1,* ) mflux_outflag
	read(1,* ) Nmflux_ave
	read(1,* ) vflux_outflag
	read(1,* ) Nvflux_ave

	!Flags to determine if periodic boundaries required
	read(1,* ) periodic(1)
	read(1,* ) periodic(2)
	read(1,* ) periodic(3)

	close(1,status='keep')      !Close input file

	if (seed(1)==seed(2)) then
		call random_seed
		call random_seed(get=seed(1:n))
	endif
	
	!Assign different random number seed to each processor
	seed = seed * irank
	!Assign seed to random number generator
	call random_seed(put=seed(1:n))

	elapsedtime = 1.d0*delta_t*Nsteps !Set elapsed time to end of simualtion

end subroutine setup_inputs

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
		read(11,* ) velocity_outflag
		read(11,* ) Nvel_ave
		read(11,* ) pressure_outflag
		read(11,* ) Nstress_ave
		read(11,* ) Nvisc_ave

		close(11,status='keep')      !Close input file

	endif

	!Broadcast data read by root to all other processors
	call MPI_BCAST(globalnp,1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(density,1,MPI_double_precision,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(rcutoff,1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
	rcutoff2= rcutoff**2         !Useful definition to save computational time
	call MPI_BCAST(inputtemperature,1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
	call MPI_BCAST(initialnunits(1),1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(initialnunits(2),1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(initialnunits(3),1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(Nsteps,1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(extrasteps,1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(delta_t,1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
	call MPI_BCAST(tplot,1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(delta_rneighbr,1,MPI_double_precision,iroot-1,MD_COMM,ierr) 
	call MPI_BCAST(vmd_outflag,1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(macro_outflag,1,MPI_integer,iroot-1,MD_COMM,ierr) 
	call MPI_BCAST(velocity_outflag,1,MPI_integer,iroot-1,MD_COMM,ierr) 
	call MPI_BCAST(Nvel_ave,1,MPI_integer,iroot-1,MD_COMM,ierr) 
	call MPI_BCAST(pressure_outflag,1,MPI_integer,iroot-1,MD_COMM,ierr)
	call MPI_BCAST(Nstress_ave,1,MPI_integer,iroot-1,MD_COMM,ierr)  
	call MPI_BCAST(Nvisc_ave,1,MPI_integer,iroot-1,MD_COMM,ierr)

	elapsedtime = elapsedtime + delta_t*extrasteps !Set elapsed time to end of simualtion
	initialstep = Nsteps         !Set plot count to final plot of last
	Nsteps = Nsteps + extrasteps !Establish final iteration step based on previous

	!print*,irank,globalnp,density,rcutoff,Nsteps,elapsedtime,inputtemperature,initialnunits, &
	!	elapsedtime,extrasteps,delta_t,tplot,delta_rneighbr
	!call MPI_BARRIER(MD_COMM, ierr)

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
	call MPI_FILE_OPEN(MD_COMM, 'final_state', & 
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

!======================================================================
!	        		OUTPUTS			              =
!======================================================================

!---------------------------------------------------------------------------------
! Write Simulation Parameter File contain all data required to completely recreate
! simulation and to be used for post processing
!
subroutine simulation_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	open(3,file=trim(file_dir)//'results/simulation_header')


	write(3,*) 'Number of Dimensions ;  nd ;', nd
	write(3,*) 'Number of Particles ;  globalnp ;', globalnp
	write(3,*) 'Time Step - delta t ;   delta_t ;',  delta_t
	write(3,*) 'Total number of steps ; Nsteps;',  Nsteps - initialstep
	write(3,*) 'Starting step of simulation ;  initialstep ;', initialstep
	write(3,*) 'Generate output file every steps ;   tplot ;',  tplot
	write(3,*) 'Density ; density ;',density
	write(3,*) 'Initial Temperature ;   inputtemperature ;',  inputtemperature
	write(3,*) 'Cut off distance ;  rcutoff ;', rcutoff
	write(3,*) 'Neighbour List Delta r ;  delta_rneighbr ;', delta_rneighbr
	write(3,*) 'Initial FCC unit size in x ;  initialunitsize_x ;', initialunitsize(1)
	write(3,*) 'Initial FCC unit size in y ;  initialunitsize_y ;', initialunitsize(2)
	write(3,*) 'Initial FCC unit size in z ;  initialunitsize_z ;', initialunitsize(3)
	write(3,*) 'Domain in x ;  globaldomain_x  ;', globaldomain(1) 
	write(3,*) 'Domain in y ;  globaldomain_y  ;', globaldomain(2) 
	write(3,*) 'Domain in z ;  globaldomain_z  ;', globaldomain(3) 
	write(3,*) 'Domain volume ;  volume ;', volume
	write(3,*) 'Periodic Boundary Conditions in x ;  periodic_x ;', periodic(1)
	write(3,*) 'Periodic Boundary Conditions in y ;  periodic_y ;', periodic(2)
	write(3,*) 'Periodic Boundary Conditions in z ;  periodic_z ;', periodic(3)
	write(3,*) 'Computational cells in x ;  globalncells_x ;',  ncells(1)*npx
	write(3,*) 'Computational cells in y ;  globalncells_y  ;', ncells(2)*npy 
	write(3,*) 'Computational cells in z ;  globalncells_z  ;', ncells(3)*npz 
	write(3,*) 'Of size in x ;  cellsidelength_x ;', cellsidelength(1)
	write(3,*) 'Of size in y ;  cellsidelength_y ;', cellsidelength(2)
	write(3,*) 'Of size in z ;  cellsidelength_z ;', cellsidelength(3)
	write(3,*) 'Number of processors in x ;  npx ;', npx
	write(3,*) 'Number of processors in y ;  npy ;', npy
	write(3,*) 'Number of processors in z ;  npz ;', npz
	write(3,*) 'Cells per Processor in x ;  nicellxl ;', nicellxl
	write(3,*) 'Cells per Processor in y ;  nicellyl ;', nicellyl
	write(3,*) 'Cells per Processor in z ;  nicellzl ;', nicellzl
	write(3,*) 'Cells per Processor including Halos in x ;  ncellxl ;', ncellxl
	write(3,*) 'Cells per Processor including Halos in y ;  ncellyl ;', ncellyl
	write(3,*) 'Cells per Processor including Halos in z ;  ncellzl ;', ncellzl
	write(3,*) '1st Random seed ;  seed_1 ;', seed(1)
	write(3,*) '2nd Random seed ;  seed_2 ;', seed(2)
	write(3,*)  'VMD flag ;  vmd_outflag ;', vmd_outflag
	write(3,*)  'macro flag ;  macro_outflag	 ;', macro_outflag	
	write(3,*)  'velocity flag ;  velocity_outflag ;', velocity_outflag
	write(3,*)  'Pressure flag ;  pressure_outflag ;', pressure_outflag
	write(3,*)  'velocity average steps ;  Nvel_ave ;', Nvel_ave
	write(3,*)  'pressure average steps ;  Nstress_ave ;', Nstress_ave
	write(3,*)  'viscosity average samples ;  Nvisc_ave ;', Nvisc_ave
	write(3,*)  'Velocity/stress Averaging Bins in x ;  globalnbins_x ;', globalnbins(1)
	write(3,*)  'Velocity/stress Averaging Bins in y ;  globalnbins_y ;', globalnbins(2)
	write(3,*)  'Velocity/stress Averaging Bins in z ;  globalnbins_z ;', globalnbins(3)
	write(3,*)  'Of size in x ;  binsize_x  ;', globaldomain(1)/globalnbins(1) 
	write(3,*)  'Of size in y ;  binsize_y  ;', globaldomain(2)/globalnbins(2) 
	write(3,*)  'Of size in z ;  binsize_z  ;', globaldomain(3)/globalnbins(3) 
	write(3,*)  'Bins per Processor in x ;  nbins_x ;', nbins(1)
	write(3,*)  'Bins per Processor in y ;  nbins_y ;', nbins(2)
	write(3,*)  'Bins per Processor in z ;  nbins_z ;', nbins(3)
	write(3,*)  'Number of Bins on outer Surface of each processor ;  nsurfacebins ;', nsurfacebins
	write(3,*)  'Domain split into Planes for Pressure Averaging ; nplanes  ;',nplanes 
	write(3,*)  'Separated by distance ;  planespacing  ;', planespacing 
	write(3,*)  'with first plane at ;  planes ;', planes(1)

	close(3,status='keep')

end subroutine simulation_header


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
	call MPI_FILE_DELETE(trim(file_dir)//'results/final_state', MPI_INFO_NULL, ierr)

	!Open file on all processors
	call MPI_FILE_OPEN(MD_COMM,trim(file_dir)//'results/final_state', & 
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
	call MPI_FILE_OPEN(MD_COMM,trim(file_dir)//'results/vmd_temp.dcd', & 
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
		case(0, 4) 	!Liquid Molecules
			Xbuf(:) = r(:,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
			Ybuf(:) = r(:,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
			Zbuf(:) = r(:,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
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
	call MPI_FILE_OPEN(MD_COMM, trim(file_dir)//'results/vmd_liquid_temp.dcd', & 
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
		case(0, 4) 	!Liquid Molecules
			Xbuf(n) = -halfdomain(1)
			Ybuf(n) = -halfdomain(2)
			Zbuf(n) = -halfdomain(3)
		case default	!Solid molecules
			Xbuf(n) = r(n,1)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
			Ybuf(n) = r(n,2)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
			Zbuf(n) = r(n,3)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
		end select
	enddo

	!================================
	!    Write Solid Molecules	=
	!================================

	!Open file on all processors
	call MPI_FILE_OPEN(MD_COMM, trim(file_dir)//'results/vmd_solid_temp.dcd', & 
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
	call MPI_FILE_OPEN(MD_COMM, trim(file_dir)//'results/vmd_temp.dcd', & 
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
!Write positions of molecules in halo to a file

subroutine parallel_io_vmd_halo
	use module_parallel_io
	implicit none

	stop "Cannot print vmd halos in parallel simulation - run in serial"

end subroutine parallel_io_vmd_halo

!---------------------------------------------------------------------------------
! Write value of last output iteration

subroutine update_simulation_progress_file
	use module_parallel_io
	implicit none

	if (irank .eq. iroot) then
		open (unit=99999, file="results/simulation_progress")
		write(99999,*) iter
		close(99999,status='keep')
	endif

end subroutine update_simulation_progress_file


!---------------------------------------------------------------------------------
! Record mass in a slice through the domain

subroutine mass_slice_io(ixyz)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer		:: ixyz,jxyz,kxyz, m, length

	!Sum over all bins using directional sub communicators
	!call SubcommSumIntVect(slice_mass, nbins(jxyz), jxyz)
	!call SubcommSumIntVect(slice_mass, nbins(kxyz), kxyz)

	!Gather on a single processor to write out
	!call MPI_gather(slice_mass,nbins(1),MPI_integer, & 
	!		globalslice_mass,nbins(1),MPI_integer, &
	!		iroot-1,icomm_xyz(2),ierr)

	m = iter/(tplot*Nmass_ave)

	!Write mass slice to file
	inquire(iolength=length) slice_mass(1:nbins(ixyz))
	open (unit=5, file="results/mslice",form="unformatted",access='direct',recl=length)
	write(5,rec=m) slice_mass(1:nbins(ixyz))
	close(5,status='keep')

end subroutine mass_slice_io

!---------------------------------------------------------------------------------
! Record mass in 3D bins throughout domain

subroutine mass_bin_io(CV_mass_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: i,j,k,n
	integer				:: m,length
	integer				:: CV_mass_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	character(4)			:: io_type
	character(13)			:: filename

	!Work out correct filename for i/o type
	write(filename, '(a9,a4)' ) trim(file_dir)//'results/m', io_type

	!Include halo surface fluxes to get correct values for all cells
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!do n = 1, nsurfacebins
	!	i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
		!Change in number of Molecules in halo cells
	!	CV_mass_out(modulo((i-2),nbins(1))+2, & 
	!		    modulo((j-2),nbins(2))+2, & 
	!		    modulo((k-2),nbins(3))+2) = & 
	!		    CV_mass_out(modulo((i-2),nbins(1))+2, & 
	!				modulo((j-2),nbins(2))+2, & 
	!				modulo((k-2),nbins(3))+2) &
	!					+ CV_mass_out(i,j,k)
	!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (io_type .eq. 'snap') then
		m = iter/(Nmflux_ave) + 1 !Initial snapshot taken
	else
		m = iter/(tplot*Nmass_ave)
	endif

	!Write mass to file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!inquire(iolength=length) CV_mass_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1)
	!open (unit=5,file=filename,form="unformatted",access='direct',recl=length)
	!write(5,rec=m) CV_mass_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1)
	!close(5,status='keep')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine mass_bin_io

!---------------------------------------------------------------------------------
! Record velocity in a slice through the domain

subroutine velocity_slice_io(ixyz)
use module_parallel_io
use calculated_properties_MD
implicit none

	integer		:: ixyz,jxyz,kxyz, m, length

	!Write mass
	call mass_slice_io(ixyz)

	!Get two directions orthogonal to slice direction
	kxyz = mod(ixyz,3)+1
	jxyz = mod(ixyz+1,3)+1

	!Sum over all bins using directional sub communicators
	!call SubcommSumVect(slice_momentum, nbins(jxyz), jxyz)
	!call SubcommSumVect(slice_momentum, nbins(kxyz), kxyz)

	!Gather on a single processor to write out
	!call MPI_gather(slice_momentum,nbins(1),MPI_double_precision,  & 
	!		globalslice_momentum,nbins(1),MPI_double_precision, & 
	!		iroot-1,icomm_xyz(2),ierr)

	!Write velocity to file
	if (irank .eq. iroot) then

		m = iter/(tplot*Nvel_ave)
		inquire(iolength=length) slice_momentum(1:nbins(ixyz),:)
		open (unit=6, file="results/vslice",form="unformatted",access='direct',recl=length)
		write(6,rec=m) slice_momentum(1:nbins(ixyz),:)
		close(6,status='keep')
	endif

end subroutine velocity_slice_io


!------------------------------------------------------------------------
!A large scale routine with each proc writing its own bins in binary
!Write velocity slice information to a file


subroutine parallel_slice_io_large_scale
	use module_parallel_io
	implicit none
	!include 'mpif.h' 

	integer					:: slicefileid, int_datasize, dp_datasize
	integer(kind=MPI_OFFSET_KIND)   	:: disp!, resultsize

	!Determine size of datatypes
	call MPI_type_size(MPI_Integer,int_datasize,ierr)
	call MPI_type_size(MPI_double_precision,dp_datasize,ierr)

	call MPI_FILE_OPEN(MD_COMM, trim(file_dir)//'results/vslice', & 
			MPI_MODE_WRONLY , & 
			MPI_INFO_NULL, slicefileid, ierr)

	!WRITE BINNED VELOCITY SUMMATION

	!Obtain displacement of current record
	disp =(iter/real((tplot*Nvel_ave),kind(0.d0))-1) *       &	!Current iteration
	      globalnbins(1) *(int_datasize+3*dp_datasize) & 		!Current iteration	
	      + nbins(1)*(3*dp_datasize)*(jblock-1)			!Processor location
	!Set current processes view
	call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_double_precision, & 
 			MPI_double_precision, 'native', MPI_INFO_NULL, ierr)
	!Write data
	call MPI_FILE_WRITE_ALL(slicefileid,slice_momentum(:,:),nd*nbins(1), MPI_double_precision, & 
			MPI_STATUS_IGNORE, ierr) !Velocity bins

	!WRITE BINNED MOLECULAR COUNT

	!Obtain displacement of current record
	disp =(iter/real((tplot*Nvel_ave),kind(0.d0))-1) *		& !Current iteration
	      globalnbins(1) * (int_datasize+dp_datasize)  & !Current iteration	
		+ nbins(1)*(int_datasize)*(jblock-1) 	& !Processor location
		+ globalnbins(1)* dp_datasize 		  !After binned velocity

	!Set current processes view
	call MPI_FILE_SET_VIEW(slicefileid, disp, MPI_Integer, & 
				MPI_Integer, 'native', MPI_INFO_NULL, ierr)

	!Write data
	call MPI_FILE_WRITE_ALL(slicefileid,slice_mass(:),nbins(1), MPI_Integer, & 
			MPI_STATUS_IGNORE, ierr) !molecular tally bins

	call MPI_FILE_CLOSE(slicefileid, ierr) 

end subroutine parallel_slice_io_large_scale

subroutine velocity_bin_io(CV_mass_out,CV_momentum_out,io_type)
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: n,m,i,j,k
	integer				:: length
	integer				:: CV_mass_out(nbins(1)+2,nbins(2)+2,nbins(3)+2)
	double precision		:: CV_momentum_out(nbins(1)+2,nbins(2)+2,nbins(3)+2,3)
	character(4)			:: io_type
	character(13)			:: filename

	!Write mass bins
	call mass_bin_io(CV_mass_out,io_type)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Work out correct filename for i/o type
	!write(filename, '(a9,a4)' ) 'results/v', io_type

	!---------------Correct for surface fluxes on halo cells---------------
	!Include halo surface fluxes to get correct values for all cells
	!do n = 1, nsurfacebins
	!	i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  

	!	!Change in Momentum in halo cells
	!	CV_momentum_out(modulo((i-2),nbins(1))+2, & 
	!		      	modulo((j-2),nbins(2))+2, & 
	!		      	modulo((k-2),nbins(3))+2,:) = & 
	!			CV_momentum_out(modulo((i-2),nbins(1))+2,& 
	!					modulo((j-2),nbins(2))+2,&
	!					modulo((k-2),nbins(3))+2,:) & 
	!						+ CV_momentum_out(i,j,k,:)
	!enddo

	if (io_type .eq. 'snap') then
		!CV_momentum_out = CV_momentum_out / (tplot*Nvflux_ave)
		m = iter/(Nvflux_ave) + 1 !Initial snapshot taken
	else
		!CV_momentum_out = CV_momentum_out / (tplot*Nvel_ave)
		m = iter/(tplot*Nvel_ave)
	endif
	!Write velocity to file
	!inquire(iolength=length) CV_momentum_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	!open (unit=6, file=filename,form="unformatted",access='direct',recl=length)
	!write(6,rec=m) CV_momentum_out(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	!close(6,status='keep')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine velocity_bin_io

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine virial_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

end subroutine virial_stress_io

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine VA_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none


end subroutine VA_stress_io

!===================================================================================
!Integrate virial pressure to get autocorrelations (Green Kubo) viscosity

subroutine viscosity_io
use module_parallel_io
use physical_constants_MD
use calculated_properties_MD
implicit none

	integer			:: m, length
	double precision	:: viscosity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!call intergrate_trap(Pxycorrel,tplot*delta_t,Nstress_ave,viscosity)

	!viscosity = (viscosity*volume)/(3.0*Nstress_ave*Nvisc_ave*inputtemperature)

	!Write viscosity to file
	!m = iter/(tplot*Nstress_ave*Nvisc_ave)
	!inquire(iolength=length) viscosity
	!open (unit=7, file="results/visc",form="unformatted",access='direct',recl=length)
	!write(7,rec=m) viscosity
	!close(7,status='keep')

	!Pxycorrel = 0.d0	!Reset Pxycorrel to zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine viscosity_io

!=================================================================================
! Record Fluxes accross surfaces of Control Volumes


!---------------------------------------------------------------------------------
! Record mass fluxes accross surfaces of Control Volumes

subroutine mass_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer		:: i,j,k,n,m,length

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Include halo surface fluxes to get correct values for all cells
	!do n = 1, nsurfacebins
	!	i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
		!Flux over halo cells
	!	mass_flux(modulo((i-2),nbins(1))+2,modulo((j-2),nbins(2))+2,modulo((k-2),nbins(3))+2,:) = & 
	!	mass_flux(modulo((i-2),nbins(1))+2,modulo((j-2),nbins(2))+2,modulo((k-2),nbins(3))+2,:) + mass_flux(i,j,k,:)
	!enddo

	!Write six CV surface mass fluxes to file
	!m = iter/(Nmflux_ave)
	!inquire(iolength=length) mass_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	!open (unit=8, file="results/mflux",form="unformatted",access='direct',recl=length)
	!write(8,rec=m) mass_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:)
	!close(8,status='keep')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine mass_flux_io

!---------------------------------------------------------------------------------
! Record momentum fluxes accross surfaces of Control Volumes

subroutine momentum_flux_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: ixyz,i,j,k,n,m,length
	double precision		:: binface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Include halo surface fluxes to get correct values for all cells
	!do n = 1, nsurfacebins
	!	i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
		!Flux over halo cells
	!	momentum_flux(	modulo((i-2),nbins(1))+2, & 
	!		      	modulo((j-2),nbins(2))+2, & 
	!		      	modulo((k-2),nbins(3))+2,:,:) = & 
	!			momentum_flux(	modulo((i-2),nbins(1))+2,& 
	!					modulo((j-2),nbins(2))+2,&
	!					modulo((k-2),nbins(3))+2,:,:) & 
	!						+ momentum_flux(i,j,k,:,:)
	!enddo

	!do ixyz = 1,3
	!	binface	      = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
	!		     	(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
	!	momentum_flux(:,:,:,:,ixyz  )=momentum_flux(:,:,:,:,ixyz  )/(binface) !Bottom
	!	momentum_flux(:,:,:,:,ixyz+3)=momentum_flux(:,:,:,:,ixyz+3)/(binface) !Top
	!enddo

	!momentum_flux = momentum_flux/(delta_t*Nvflux_ave)
	!Write momnetum flux to file
	!m = iter/(Nvflux_ave)
	!inquire(iolength=length) momentum_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	!open (unit=9, file="results/vflux",form="unformatted",access='direct',recl=length)
	!write(9,rec=m) momentum_flux(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	!close(9,status='keep')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine momentum_flux_io

!---------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MOP_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none


end subroutine MOP_stress_io

!---------------------------------------------------------------------------------
! Record stress accross surfaces of Control Volumes

subroutine surface_stress_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer				:: ixyz,i,j,k,n,m,length
	double precision,dimension(3)	:: binface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Include halo surface stresses to get correct values for all cells
	!do n = 1, nsurfacebins
	!	i = surfacebins(n,1); j = surfacebins(n,2); k = surfacebins(n,3)  
		!Set Stresses to value of halo cells
	!	Pxyface(	modulo((i-2),nbins(1))+2, & 
	!		      	modulo((j-2),nbins(2))+2, & 
	!		      	modulo((k-2),nbins(3))+2,:,:) = & 
	!			Pxyface(	modulo((i-2),nbins(1))+2,& 
	!					modulo((j-2),nbins(2))+2,&
	!					modulo((k-2),nbins(3))+2,:,:) & 
	!						     + Pxyface(i,j,k,:,:)

	!enddo


	!do ixyz = 1,3
	!	binface(ixyz) = (domain(modulo(ixyz  ,3)+1)/nbins(modulo(ixyz  ,3)+1))* & 
	!		     	(domain(modulo(ixyz+1,3)+1)/nbins(modulo(ixyz+1,3)+1))
	!	Pxyface(:,:,:,:,ixyz  ) = 0.25d0 * Pxyface(:,:,:,:,ixyz  )/binface(ixyz) !Bottom
	!	Pxyface(:,:,:,:,ixyz+3) = 0.25d0 * Pxyface(:,:,:,:,ixyz+3)/binface(ixyz) !Top
	!enddo

	!Integration of stress using trapizium rule requires multiplication by timestep
	!Pxyface = delta_t*Pxyface!/Nvflux_ave

	!Write surface pressures to file
	!m = iter/(Nvflux_ave)
	!inquire(iolength=length) Pxyface(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	!open (unit=9, file="results/psurface",form="unformatted",access='direct',recl=length)
	!write(9,rec=m) Pxyface(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,:,:)
	!close(9,status='keep')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NEED TO PARALLELISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine surface_stress_io


