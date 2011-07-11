!======================================================================
!	        	Parallel Read Write Subroutines               =
!======================================================================
!Serial emulation of parallel i/o

module module_parallel_io
	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD

end module


!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd
	use module_parallel_io
	implicit none

	integer				:: i, n
	real,dimension(np)		:: Xbuf, Ybuf, Zbuf

	Xbuf(:) = r(:,1)
	Ybuf(:) = r(:,2)
	Zbuf(:) = r(:,3)

	i = iter / tplot

	open (unit=17, file="results/vmd_temp.dcd",access='direct',recl=1)
	
	do n = 1,np
		write(17,rec=(i-1)*nd*np+n) Xbuf(n)
		write(17,rec=(i-1)*nd*np+n+np) Ybuf(n)
		write(17,rec=(i-1)*nd*np+n+2*np) Zbuf(n)
	enddo

	close(17,status='keep')

end subroutine parallel_io_vmd

!------------------------------------------------------------------------
!Write positions of molecules to a file

subroutine parallel_io_vmd_sl
	use module_parallel_io
	implicit none

	integer				:: i, n
	real,dimension(np)		:: Xbuf, Ybuf, Zbuf
	real,dimension(3)		:: rhalfdomain

	Xbuf(:) = r(:,1)
	Ybuf(:) = r(:,2)
	Zbuf(:) = r(:,3)

	rhalfdomain(1) = -halfdomain(1)
	rhalfdomain(2) = -halfdomain(2)
	rhalfdomain(3) = -halfdomain(3)

	i = iter / tplot


	!---Write liquid molecules---
	open (unit=17, file="results/vmd_liquid_temp.dcd",access='direct',recl=1)
	
	do n = 1,np
		select case (tag(n))
		case(0)
			write(17,rec=(i-1)*nd*np+n) Xbuf(n)
			write(17,rec=(i-1)*nd*np+n+np) Ybuf(n)
			write(17,rec=(i-1)*nd*np+n+2*np) Zbuf(n)
		case(1:)
			write(17,rec=(i-1)*nd*np+n) rhalfdomain(1)
			write(17,rec=(i-1)*nd*np+n+np) rhalfdomain(2)
			write(17,rec=(i-1)*nd*np+n+2*np) rhalfdomain(3)
		end select
	enddo

	close(17,status='keep')

	!---Write solid molecules---
	open (unit=17, file="results/vmd_solid_temp.dcd",access='direct',recl=1)

	do n = 1,np
		select case (tag(n))
		case(0)
			write(17,rec=(i-1)*nd*np+n) rhalfdomain(1)
			write(17,rec=(i-1)*nd*np+n+np) rhalfdomain(2)
			write(17,rec=(i-1)*nd*np+n+2*np) rhalfdomain(3)
		case(1:)
			write(17,rec=(i-1)*nd*np+n) Xbuf(n)
			write(17,rec=(i-1)*nd*np+n+np) Ybuf(n)
			write(17,rec=(i-1)*nd*np+n+2*np) Zbuf(n)
		end select
	enddo

	close(17,status='keep')

end subroutine parallel_io_vmd_sl

!------------------------------------------------------------------------
!Write positions of molecules in halo to a file

subroutine parallel_io_vmd_halo
	use module_parallel_io
	implicit none

	integer				:: i, n
	real,dimension(halo_np)		:: Xbuf, Ybuf, Zbuf
	real,dimension(3)		:: rhalfdomain

	Xbuf(:) = r(np+1:np+halo_np,1)
	Ybuf(:) = r(np+1:np+halo_np,2)
	Zbuf(:) = r(np+1:np+halo_np,3)

	rhalfdomain(1) = -halfdomain(1)
	rhalfdomain(2) = -halfdomain(2)
	rhalfdomain(3) = -halfdomain(3)

	i = iter / tplot

	open (unit=17, file="results/vmd_halo_temp.dcd",access='direct',recl=1)
	
	do n = 1,halo_np
		write(17,rec=(i-1)*nd*extralloc+n) Xbuf(n)
		write(17,rec=(i-1)*nd*extralloc+n+extralloc) Ybuf(n)
		write(17,rec=(i-1)*nd*extralloc+n+2*extralloc) Zbuf(n)
	enddo

	!Write up to extralloc to give constant number of molecule in vmd file
	do n = halo_np,halo_np+extralloc
		write(17,rec=(i-1)*nd*extralloc+n+3*extralloc) rhalfdomain(1)
		write(17,rec=(i-1)*nd*extralloc+n+4*extralloc) rhalfdomain(2)
		write(17,rec=(i-1)*nd*extralloc+n+5*extralloc) rhalfdomain(3)
	enddo

	close(17,status='keep')

end subroutine parallel_io_vmd_halo

!------------------------------------------------------------------------------
!Serial version of parallel code to print final_state for restart

subroutine parallel_io_final_state
	use module_parallel_io
	implicit none

	integer :: ixyz,n

	!Rebuild simulation before recording final state
	call sendmols			   !Exchange particles between processors
	call linklist_deallocateall	   !Deallocate all linklist components
	call assign_to_cell	  	   !Re-build linklist every timestep
	call messenger_updateborders	   !Update borders between processors
	call assign_to_halocell		   !Re-build linklist
	call assign_to_neighbourlist	   !Setup neighbourlist

	!Remove previous final state file
	open(13,file='results/final_state')
	close(13,status='delete')

	open(13,file='results/final_state', form='unformatted',access='direct',recl=2)
	!Written in this form so each molecule's information is together to allow 
	!re-allocation to seperate processors
	do n=1,np
		!Write particle n's positions
		do ixyz=1,nd
			write(13,rec=6*(n-1)+(ixyz-1)+1) r(n,ixyz) 
		enddo
		!Write particle n's velocities
		do ixyz=1,nd
			write(13,rec=6*(n-1)+(ixyz-1)+4) v(n,ixyz)   
		enddo

		!print'(a,i8,6f10.5)','written', n, r(n,:), v(n,:)
	enddo

	write(13,rec=6*np+1) density           !Density of system
	write(13,rec=6*np+2) rcutoff           !Cut off distance for particle interaction
	write(13,rec=6*np+3) inputtemperature  !Define initial temperature
	write(13,rec=6*np+4) delta_t           !Size of time step
	write(13,rec=6*np+5) elapsedtime       !Total elapsed time of all restarted simulations

	close(13,status='keep') !Close final_state file

	!Re-open file with different record length to write Integer Data
	open(13,file='results/final_state', form='unformatted',access='direct',recl=1)

	write(13,rec=12*np+11) np               !Number of particles
	write(13,rec=12*np+12) initialnunits(1) !x dimension split into number of cells
	write(13,rec=12*np+13) initialnunits(2) !y dimension box split into number of cells
	write(13,rec=12*np+14) initialnunits(3) !z dimension box split into number of cells
	write(13,rec=12*np+15) Nsteps           !Number of computational steps
	write(13,rec=12*np+16) tplot            !Frequency at which to record results
	write(13,rec=12*np+17) seed(1)          !Random number seed value 1
	write(13,rec=12*np+18) seed(2)          !Random number seed value 2

	close(13,status='keep') !Close final_state file

end subroutine parallel_io_final_state


!----------------------------------------------------------------------------------
!
!                                Restart Simulation inputs
! Set up inputs based on the final state of a previous simulation
!
!---------------------------------------------------------------------------------

subroutine setup_restart_inputs
	use module_parallel_io
	implicit none

	integer				:: n, k
	integer 			:: extrasteps
	integer				:: filesize,int_filesize,dp_filesize
	integer 			:: checkint
	double precision 		:: checkdp

	!Allocate random number seed
	call random_seed
	call random_seed(size=n)
	allocate(seed(n))

	!Call function library to get file size
	call get_file_size('final_state',filesize)

	!File size is in bytes and integer fortran records are blocks of 4 bytes
	int_filesize = filesize/4

	!Open file to read integers
	open(13,file='final_state', form="unformatted", access="direct",recl=1)

	read(13,rec=int_filesize-7) np		    !Number of particles
	globalnp = np				    !Global np and local np same in serial
	read(13,rec=int_filesize-6) initialnunits(1)!x dimension split into number of cells
	read(13,rec=int_filesize-5) initialnunits(2)!y dimension box split into number of cells
	read(13,rec=int_filesize-4) initialnunits(3)!z dimension box split into number of cells
	read(13,rec=int_filesize-3) Nsteps  	    !Number of elapsed computational steps

	close(13,status='keep')

	!File size is in bytes and fortran double precision records are blocks of 8 bytes
	dp_filesize = filesize/8

	!Reopen file to read doubles
	open(13,file='final_state', form="unformatted", access="direct",recl=2)

	read(13,rec=dp_filesize-8) density		!Density of system
	read(13,rec=dp_filesize-7) rcutoff		!Cut off distance for particle interaction
	rcutoff2= rcutoff**2         !Useful definition to save computational time
	read(13,rec=dp_filesize-6) inputtemperature	!Define initial temperature
	read(13,rec=dp_filesize-4) elapsedtime   	!Elapsed simulation time to date

	close(13,status='keep')

	!Check if values from input file are different and alert user - all processors have
	!read the same file so only need to check on one processor

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

	elapsedtime = elapsedtime + delta_t*extrasteps !Set elapsed time to end of simualtion
	initialstep = Nsteps         !Set plot count to final plot of last
	Nsteps = Nsteps + extrasteps !Establish final iteration step based on previous

	print*,irank,globalnp,density,rcutoff,Nsteps,elapsedtime,inputtemperature,initialnunits, &
		elapsedtime,extrasteps,delta_t,tplot,delta_rneighbr

end subroutine setup_restart_inputs
!----------------------------------------------------------------------------------
!
!                                Restart Microstate
! Set up position and velocity of each molecule (the microstate) based on the final
! state of a previous simulation
!
!---------------------------------------------------------------------------------

subroutine setup_restart_microstate
	use module_parallel_io
	implicit none

	integer :: ixyz, n

	!Open file at first recorded value
	open(13,file='final_state', form='unformatted', access='direct',recl=2)

	!Read positions
	do n=1,globalnp
		!Read position from file
		do ixyz=1,nd
			read(13,rec=6*(n-1)+(ixyz-1)+1) r(n,ixyz)
		enddo
	enddo

	call setup_tag				!Setup location of fixed molecules
	do n = 1,np
		call read_tag(n)		!Read tag and assign properties
	enddo

	do n=1,globalnp
		!Read corresponding velocities
		do ixyz=1,nd
			read(13,rec=6*(n-1)+(ixyz-1)+4) v(n,ixyz)
		enddo

		!print'(a,i8,6f10.5)','read', n, r(n,:), v(n,:)
	enddo

	close(13,status='delete') !Close final_state file

end subroutine setup_restart_microstate

!---------------------------------------------------------------------------------

subroutine local_temperature_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	open(3,file='results/local_temperature') !Open local temperature file
	write(3,'(i8)') nbins(1)

end subroutine local_temperature_header


!---------------------------------------------------------------------------------

subroutine velocity_average_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	open(10,file='results/v_average') !Open velocity slice file

end subroutine velocity_average_header

!---------------------------------------------------------------------------------

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


!---------------------------------------------------------------------------------

subroutine parallel_slice_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	open(18,file='results/vslice') !Open velocity slice file
	write(18,'(i10,a,f10.5)') globalnp,';', globaldomain(2)
	write(18,'(i10,a,i10)') globalnbins(1),';', Nsteps
	write(18,'(i10,a,i10)') tplot,';', Nslice_ave

end subroutine parallel_slice_header

!---------------------------------------------------------------------------------

subroutine parallel_slice_io
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	integer	:: ibin, jbin, kbin

	ibin = 1
	kbin = 1
	do jbin =1 , nbins(2)
		write(18,'(f20.5,a,i10)') meanvel(jbin),';', countvel(jbin)
	enddo

end subroutine parallel_slice_io

!---------------------------------------------------------------------------------

subroutine stress_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	open(2,file='results/pressure_visc') !Open pressure field file

end subroutine stress_header

!---------------------------------------------------------------------------------

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

subroutine control_volume_header
	use module_parallel_io
	use calculated_properties_MD
	implicit none

	open(4,file='results/control_volume') !Open control volume file

end subroutine control_volume_header

!---------------------------------------------------------------------------------

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
