!==========================================================================
!                     CALCULATE AND RECORD PROPERTIES
!==========================================================================
! Evaluate properties of interest from the simulation. 
! This includes macroscopic properties including sum of velocity/velocity^2, 
! kinetic energy, potential energy, temperature and pressure.
! Also calculated is the diffusion; velocity distributions and Boltzmanns
! H function; the radial distribution function; the velocity magnitude by cells 
! to be used for contours; plane based uni-directional velocity slices; 
! the stress tensor calculated from a virial (bulk or local for a homogenous fluid),
! Volume averaged (local assignment based on proportion of force/momentum is in a 
! given cell) or Method Of Planes (force/momentum flux over a plane).
! Many of these outputs are interpreted using MATLAB scripts included
!
! simulation_record  			   			Top level function choosing logging based on inputs
! evaluate_macroscopic_properties           Macroscopic properties
! evaluate_properties_vdistribution			Maxwell Boltzmann distribution
! evaluate_properties_radialdist			Radial distribution
! evaluate_properties_diffusion				Diffusion

! mass_averaging							slice/CV
! cumulative_mass							slice/CV
! velocity_averaging						slice/CV
! cumulative_velocity						slice/CV
! pressure_averaging						virial/CV
! cumulative_pressure						virial/CV
! simulation_compute_kinetic_VA		
! simulation_compute_kinetic_VA_cells
! mass_flux_averaging						CV_surface
! cumulative_mass_flux						CV_surface
! mass_snapshot								CV
! momentum_flux_averaging					MOP_plane/CV_surface
! cumulative_momentum_flux					MOP_plane/CV_surface
! momentum_snapshot							CV
! pressure_tensor_forces					virial
! pressure_tensor_forces_VA					CV
! control_volume_forces						CV_surface
! control_volume_stresses					CV_surface
! pressure_tensor_forces_MOP				MOP_plane
! evaluate_properties_cellradialdist

!===================================================================================
! 			FULL DOMAIN VIRIAL STRESS CALCULATION
! Calculate pressure_tensor for entire domain- average over all xy, yz and xz 
! components performed when calculating the interactions of particles
!===================================================================================
!===================================================================================
! 			VOLUME AVERAGE STRESS TENSOR
! Stress calculated for each cell by assinging to cells based on location of 
! molecules (kinetic) and fraction of pairwise interaction in cell (potential)
! First used in J. F. Lutsko, J. Appl. Phys 64(3), 1152 (1988) and corrected in
! J. Cormier, J.M. Rickman and T. J. Delph (2001) Journal of Applied Physics,
! Volume 89, No. 1 "Stress calculations in atomistic simulations of perfect 
! and imperfect solids" 
!===================================================================================
!===================================================================================
!	PRESSURE TENSOR CALCULATED USING METHOD OF PLANES (MOP)
! See B.D.Todd, D.J.Evans and P.J.Daivis (1995) Phys Review E, Vol 52,2
!===================================================================================

module module_record

	use interfaces
	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use calculated_properties_MD
	use polymer_info_MD
	use librarymod

	double precision :: vel

end module module_record
!==========================================================================
! Top level routine which calls all recording statements throughout the 
! code.

subroutine simulation_record
	use module_record
	use CV_objects, only : CVcheck_mass, CV_debug
	implicit none

	integer			:: vmd_iter
	integer,save	:: i = 1

	if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
		call mass_flux_averaging(mflux_outflag)				!Average mass flux before movement of particles
		call momentum_flux_averaging(vflux_outflag)         !Average momnetum flux after movement of particles
		call energy_flux_averaging(eflux_outflag)			!Average energy flux after movement of particles
	endif

	!---------------Only record every tplot iterations------------------------
	!--------------Only evaluate every teval iterations-----------------------
	if (mod(iter,tplot) .eq. 0) then
		!Do nothing	
	else if (mod(iter,teval) .eq. 0) then
		call evaluate_microstate_pressure
		return
	else
		return	
	end if 
	!--------------Only evaluate every teval iterations-----------------------
	!---------------Only record every tplot iterations------------------------

	!Parallel output for molecular positions
	if (vmd_outflag.ne.0 .and. size(vmd_intervals,2).ge.i) then
		vmd_iter = iter-initialstep+1
		if (vmd_iter.ge.vmd_intervals(1,i).and.vmd_iter.lt.vmd_intervals(2,i)) then
			select case(vmd_outflag)
			case(1)
				call parallel_io_vmd
			case(2)
				call parallel_io_vmd_sl
			case(3)
				call parallel_io_vmd
				call parallel_io_vmd_halo
			case(4)
				call parallel_io_vmd_true
			case default
				call error_abort('Unrecognised vmd_outflag in simulation_record')
			end select
			vmd_count = vmd_count + 1
		else if (vmd_iter.ge.vmd_intervals(2,i)) then
			i = i + 1			
		endif
	endif
	
	if (potential_flag.eq.1) then
		select case(etevtcf_outflag)
		case (1)
			call etevtcf_calculate_parallel
		case (2)
			call etevtcf_calculate_parallel
			call etev_io
		case default
		end select
		
		select case(r_gyration_outflag)
		case (1)
			call r_gyration_calculate_parallel
		case (2)
			call r_gyration_calculate_parallel
			call r_gyration_io
		case default
		end select
	end if

	!Obtain each processe's subdomain's macroscopic 
	!properties; gather on root process and record
	if (macro_outflag .ne. 0) then

		call evaluate_macroscopic_properties

		select case(integration_algorithm)
		case(leap_frog_verlet)
		   call print_macroscopic_properties(iter-1)   
		case(velocity_verlet)
		   call print_macroscopic_properties(iter) 
		end select

		select case(macro_outflag)
		case(2,4)
			call macroscopic_properties_record
		case default
		end select

	endif

	!Obtain and record velocity distributions
	if (vdist_flag .eq. 1) then
		call evaluate_properties_vdistribution
	endif

	!Obtain and record radial distributions
	select case (rdf_outflag)
	case(1)
		call evaluate_properties_rdf
		call rdf_io
	case(2)
		call evaluate_properties_rdf3d
		call rdf3d_io
	case default
	end select

	!Obtain and record static structure factor
	select case (ssf_outflag)
	case(1)
		call evaluate_properties_ssf
		call ssf_io
	case default
	end select

	!Obtain and record velocity and mass
	if (velocity_outflag .ne. 0) call velocity_averaging(velocity_outflag)

	!Obtain and record mass only
	if (velocity_outflag .eq. 0 .and. &
		mass_outflag .ne. 0) call mass_averaging(mass_outflag)

	!Obtain and record temperature
	if (temperature_outflag .ne. 0)	call temperature_averaging(temperature_outflag)

	!Obtain and record molecular diffusion
	!call evaluate_properties_diffusion

	!Calculate pressure tensor
	if (pressure_outflag .ne. 0) then
		call pressure_averaging(pressure_outflag)
	end if

	call update_simulation_progress_file

end subroutine simulation_record

!==========================================================================
!Calculate kinetic and potential energy as well as temperature and pressure

subroutine evaluate_macroscopic_properties
	use module_record
	implicit none

	integer :: n,ixyz

	vsum  = 0.d0                                                ! Reset all sums
	v2sum = 0.d0                                                ! Reset all sums

	! Potential Component
	select case(potential_flag)
	case(0)
		potenergysum = sum(potenergymol(1:np))
		call globalSum(potenergysum)
	case(1)
		potenergysum_LJ = sum(potenergymol_LJ(1:np))
		potenergysum_FENE = sum(potenergymol_FENE(1:np))
		potenergysum = sum(potenergymol_LJ(1:np) + potenergymol_FENE(1:np))
		call globalSum(potenergysum_LJ)
		call globalSum(potenergysum_FENE)
		call globalSum(potenergysum)
	case default
		call error_abort("Unrecognised potential flag in simulation_record")
	end select

	! Kinetic Component
	select case(integration_algorithm)
	case(leap_frog_verlet)
		do n = 1, np 									! Loop over all particles
		do ixyz = 1, nd									! Loop over all dimensions
			vel   = v(ixyz,n) + 0.5d0*a(ixyz,n)*delta_t	! Velocity must shifted half a timestep
			vsum  = vsum + vel							! Add up all molecules' velocity components
			v2sum = v2sum + vel**2						! Add up all molecules' velocity squared components  
		enddo
		enddo
	case(velocity_verlet) 								! If velocity Verlet algorithm
		do n = 1, np
		do ixyz = 1, nd
			vel   = v(ixyz,n)
			vsum  = vsum+vel
			v2sum = v2sum + vel**2          			! Sum all velocity squared components
		enddo
		enddo
	end select
        
	!Obtain global sums for all parameters
	call globalSum(vsum)
	call globalSum(v2sum)
	virial = sum(virialmol(1:np))
	call globalSum(virial)

	kinenergy   = (0.5d0 * v2sum) / real(globalnp,kind(0.d0))
	potenergy   = potenergysum /(2.d0*real(globalnp,kind(0.d0))) + Potential_sLRC !N.B. extra 1/2 as all interactions calculated
	if (potential_flag.eq.1) then
		potenergy_LJ= potenergysum_LJ/(2.d0*real(globalnp,kind(0.d0))) + Potential_sLRC
		potenergy_FENE= potenergysum_FENE/(2.d0*real(globalnp,kind(0.d0)))
	end if
	!print'(4(a,f18.8))', ' <PE>= ',potenergy, & 
	!						 ' std(PE) = ',sqrt(sum((potenergymol(1:np)-potenergy)**2)/(2.d0*real(globalnp,kind(0.d0)))), & 
	!						 ' max= ',maxval(potenergymol(1:np)),' min= ',minval(potenergymol(1:np))
	if (maxval(potenergymol(1:np)) .gt. 10000) then
		print*, 'np = ', np
		print*, 'max(potenergymol) = ', maxval(potenergymol(1:np))
		do n=1,np
			write(3000+irank,'(i10,f28.4,3f10.4)'), n , potenergymol(n), r(:,n)
		enddo
		print*, 'Simulation aborted because max PE has reached an unreasonably high value.'
		print*, 'Inspect fort.(3000+irank) for n, potenergymol, r.'
		stop
	endif
	totenergy   = kinenergy + potenergy
	temperature = v2sum / real(nd*globalnp,kind(0.d0))
	if (any(periodic.gt.1)) temperature = get_temperature_PUT()
	pressure    = (density/(globalnp*nd))*(v2sum+virial/2.d0) + Pressure_sLRC !N.B. virial/2 as all interactions calculated

end subroutine evaluate_macroscopic_properties

subroutine evaluate_microstate_pressure
	use module_record
	implicit none

	integer :: n
	real(kind(0.d0)) :: vtemp(nd)

	v2sum = 0.d0

	! Kinetic part of Virial
	select case(integration_algorithm)
	case(leap_frog_verlet)

		do n = 1, np
			vtemp(:) = v(:,n) + 0.5d0*a(:,n)*delta_t ! Velocity must shifted half a timestep
			v2sum = v2sum + dot_product(vtemp,vtemp)
		enddo

	case(velocity_verlet)

		do n = 1, np
			v2sum = v2sum + dot_product(v(:,n),v(:,n))
		enddo

	end select

	! Configurational part of Virial
	virial = sum(virialmol(1:np))
	call globalSum(virial)

	! Instantaneous pressure 
	pressure = (density/(globalnp*nd))*(v2sum+virial/2.d0) + Pressure_sLRC !N.B. virial/2 as all interactions calculated

end subroutine evaluate_microstate_pressure


subroutine print_macroscopic_properties(it)
use module_record
implicit none

	integer, intent(in) :: it

	if (irank .eq. iroot) then
		select case(potential_flag)
		case(0)
			select case(macro_outflag)
			case(1:2)
				print '(1x,i8,a,f10.3,a,e12.2,a,f10.2,a,f7.3,a,f19.15,a,f19.15,a,f19.15,a,f10.4)', &
				it,';', simtime,';',vsum,';', v2sum,';', temperature,';', &
				kinenergy,';',potenergy,';',totenergy,';',pressure
			case(3:4)
				print '(1x,i7,a,f9.3,a,e9.2,a,f7.3,a,f8.4,a,f8.4,a,f8.4,a,f8.4)', &
				it,';', simtime,';',vsum,';', temperature,';', &
				kinenergy,';',potenergy,';',totenergy,';',pressure
			case default
			end select
		case(1)
			select case(macro_outflag)
			case(1:2)
				print '(1x,i8,a,f10.3,a,e10.3,a,f10.2,a,f7.3,a,f15.11,a,f15.11,a,f15.11,a,f10.4,a,f7.4,a,f9.3)', &
				it,';',simtime,';',vsum,';', v2sum,';', temperature,';', &
				kinenergy,';',potenergy,';',totenergy,';',pressure,';',etevtcf,';',R_g
			case(3:4)
				print '(1x,i7,a,f8.3,a,e8.1,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f6.3,a,f5.1)', &
				it,';', simtime,';',vsum,';', temperature,';', &
				kinenergy,';',potenergy,';',totenergy,';',pressure,';',etevtcf,';',R_g
			case default
			end select
		case default
			call error_abort("Invalid potential flag in input file")
		end select
	end if

end subroutine print_macroscopic_properties

subroutine print_mol_escape_error(n)
	use arrays_MD, only: r,v,a,tag
	use computational_constants_MD, only: irank, iblock, jblock, kblock, iter, ensemble, tag_move
	use messenger, only: globalise
	use librarymod, only : get_new_fileunit
	implicit none

	integer, intent(in) :: n
	real(kind(0.d0)) :: rglob(3)
	character(len=19)   :: filename
	character(len=4)    :: ranknum
	integer             :: fileunit
	logical             :: op


	rglob = globalise(r(:,n))

	! Store processor-specific filename
	write(ranknum, '(i4)'), irank	
	write(filename,'(a15,a)'), "mol_escape_err.", trim(adjustl(ranknum))

	! Warn user of escape
	print('(a,i8,a,i4,3a)'),' Molecule ',n,' on process ', &
		  irank, ' is outside the domain and halo cells. Inspect ', &
		  filename, 'for more information.'

	! Find unused file unit number
	fileunit = get_new_fileunit()

	! Open file and write escape info
	open(unit=fileunit,file=filename,access='append')

		write(fileunit,'(a,i6,a,i4,a)'),' Molecule ',n,' on process ', &
			  irank, ' is outside the domain and halo cells.'
		write(fileunit,'(a,i8,a)'),' At iteration ',iter,' it is located at: '
		write(fileunit,'(a,e20.5)'),   '    rx: ', r(1,n)
		write(fileunit,'(a,e20.5)'),   '    ry: ', r(2,n)
		write(fileunit,'(a,e20.5,a)'), '    rz: ', r(3,n), ','
		write(fileunit,'(a,e20.5)'),   '    globalrx: ', rglob(1) 
		write(fileunit,'(a,e20.5)'),   '    globalry: ', rglob(2) 
		write(fileunit,'(a,e20.5,a)'), '    globalrz: ', rglob(3), ','
		write(fileunit,'(a)'),         ' with velocity: '
		write(fileunit,'(a,e20.5)'),   '    vx: ', v(1,n)
		write(fileunit,'(a,e20.5)'),   '    vy: ', v(2,n)
		write(fileunit,'(a,e20.5)'),   '    vz: ', v(3,n)
		write(fileunit,'(a)'),         ' and acceleration: '
		write(fileunit,'(a,e20.5)'),   '    ax: ', a(1,n)
		write(fileunit,'(a,e20.5)'),   '    ay: ', a(2,n)
		write(fileunit,'(a,e20.5)'),   '    az: ', a(3,n)
		if (ensemble.eq.tag_move) then
			write(fileunit,'(a)'),         ' Molecular tag: '
			write(fileunit,'(a,i4)'),   '    tag: ', tag(n)
		end if
		write(fileunit,'(a,3i4)'),' Processor coords: ',iblock,jblock,kblock

	close(fileunit,status='keep')

end subroutine print_mol_escape_error


!===================================================================================
!Molecules grouped into velocity ranges to give vel frequency distribution graph
!and calculate Boltzmann H function
subroutine evaluate_properties_vdistribution
	use module_record
	implicit none

	integer          :: n
	integer          :: cbin
	double precision :: Hfunction
	double precision,dimension(:),allocatable :: vmagnitude

	!Calculate matrix of velocity magnitudes and bin in histogram
	allocate(vmagnitude(np))
	do n = 1, np    ! Loop over all particles
		vmagnitude(n) = dot_product(v(:,n),v(:,n)) 
		!Assign to bins using integer division
 		cbin = ceiling(vmagnitude(n)/binsize)   !Establish current bin
		if (cbin > nbins(1)) cbin = nbins(1) 		!Prevents out of range values
		if (cbin < 1 ) cbin = 1        		!Prevents out of range values
		vfd_bin(cbin) = vfd_bin(cbin)+1		!Add one to current bin
	enddo
	deallocate(vmagnitude)

	!Normalise bins to use for output and H function - so all bins add up to one
	!Total stored molecules must be equal to number of molecules (np) times the
	!number of times bins have been assigned (iter/tplot) with the 1.d0 to avoid integer
	!division
	normalisedvfd_bin=0
	do n=1,nbins(1)
		normalisedvfd_bin(n) = vfd_bin(n)/((np*iter)*1.d0/tplot) 
	enddo

	!Calculate Boltzmann H function using discrete defintion as in
	!Rapaport p37. N.B. velocity at middle or range is used for vn
	Hfunction = 0.d0 !Reset H function before re-calculating
	do n=1,nbins(1)
		if (normalisedvfd_bin(n) .ne. 0) then
			Hfunction=Hfunction+normalisedvfd_bin(n)*log(normalisedvfd_bin(n)/(((n-0.5d0)*binsize)**(nd-1)))
		endif
	enddo
	write(12,'(a,i5, a, f20.10)') 'Boltzmann H function at iteration ', iter , ' is ', Hfunction

	!Write values of bin to file to follow evolution of distribution function
	write(12,'(a)') 'Velocity frequency distribution'
	do n=1,nbins(1) 
		write(12,'(2(f10.5))') (n-0.5d0)*binsize, normalisedvfd_bin(n) 
	enddo

end subroutine evaluate_properties_vdistribution

!==============================================================================
!Calculate Radial distribution function (rdf) using all pairs

subroutine evaluate_properties_rdf
use module_record
implicit none

	integer                                    :: i,j,bin
	integer, save                              :: hist_count=1
	double precision                           :: rmag,dr,Nideal,Nbin,dV
	double precision, dimension(nd)            :: rij

	!Shell "width"
	dr = rdf_rmax/real(rdf_nbins,kind(0.d0))

	!Add to histogram of radial separations
	do i = 1,np-1
	do j = i+1,np

		!Find distance between i and j (rmag)
		rij(:) = r(:,i) - r(:,j)
		rij(:) = rij(:) - domain(:)*anint(rij(:)/domain(:))
		rmag   = sqrt(dot_product(rij,rij))
	
		!Assign to bin, one entry for each molecule
		bin    = int(rmag/dr) + 1
		if (bin.le.rdf_nbins) rdf_hist(bin) = rdf_hist(bin) + 2

	end do
	end do

	!Normalise histogram to find g(r) (rdf)
	do bin = 1,rdf_nbins
		rmag     = real(bin - 1.d0)*dr        !rmag now means r in g(r)
		select case(nd)
		case(2)
			dV   = pi*((rmag+dr)**2.d0 - rmag**2.d0)
		case(3)
			dV   = (4.d0*pi/3.d0)*((rmag+dr)**3.d0 - rmag**3.d0)
		end select
		Nideal   = density*dV
		Nbin     = real(rdf_hist(bin))/real(np*hist_count)
		rdf(bin) = Nbin/Nideal
	end do

	hist_count = hist_count + 1

end subroutine evaluate_properties_rdf

subroutine evaluate_properties_rdf3d
use module_record
implicit none

	integer          :: i,j,bin,ixyz
	integer, save    :: hist_count=1
	double precision :: x,dx,Nideal,Nbin,dV
	double precision :: xij

	!Shell "width"
	dx = rdf_rmax/real(rdf_nbins,kind(0.d0))

	!Add to histogram of radial separations
	do i = 1,np-1
	do j = i+1,np

		do ixyz = 1,nd

			!Find distance between i and j (rmag)
			xij = r(ixyz,i) - r(ixyz,j)
			xij = xij - domain(ixyz)*anint(xij/domain(ixyz))
			x   = abs(xij)			

			!Assign to bin, one entry for each molecule
			bin    = int(x/dx) + 1
			if (bin.le.rdf_nbins) then
				rdf3d_hist(bin,ixyz) = rdf3d_hist(bin,ixyz) + 2
			end if

		end do

	end do
	end do

	!Normalise histogram to find g(r) (rdf)
	do bin = 1,rdf_nbins
		do ixyz = 1,nd

			x = real(bin - 1.d0)*dx
			select case(nd)
			case(2)
				dV = 2.d0*dx*domain(1)*domain(2)/domain(ixyz) 
			case(3)
				!Both sides
				dV = 2.d0*dx*domain(1)*domain(2)*domain(3)/domain(ixyz)
			end select
			Nideal = density*dV
			Nbin = real(rdf3d_hist(bin,ixyz))/real(np*hist_count)
			rdf3d(bin,ixyz) = Nbin/Nideal

		end do
	end do

	hist_count = hist_count + 1

end subroutine evaluate_properties_rdf3d

!Use a test particle to plot the potential field
subroutine simulation_write_potential_field(xmin,xmax,ymin,ymax,z,res,filenum)
	implicit none

	integer, intent(in) :: res, filenum	
	real(kind(0.d0)), intent(in) :: xmin, xmax, ymin, ymax, z

	integer :: i,j
	real(kind(0.d0)) :: x, y, dx, dy
	real(kind(0.d0)) :: U
	real(kind(0.d0)) :: dummy(3)

	dx = (xmax - xmin)/real(res)
	dy = (ymax - ymin)/real(res)

	do i = 0,res

		do j = 0,res

				x  = xmin + i*dx
				y  = ymin + j*dy 

				call compute_force_and_potential_at((/x,y,z/),U,dummy)

				write(filenum,*) x, y, U

		end do

		write(filenum,*) ' ' 

	end do

end subroutine simulation_write_potential_field

!============================================================================!
! Evaluate static structure factor S(k) in three dimensions
subroutine evaluate_properties_ssf
use module_record
implicit none

	integer            :: i,j,posi,posj
	integer            :: n
	integer, save      :: hist_count=1
	double precision   :: c,s
	double precision, dimension(nd)  :: k,dk

	!Wavevectors restricted to fit inside periodic box lengths
	dk(:) = 2.d0*pi/domain(:)

	k = 0.d0
	do i=-ssf_nmax,ssf_nmax
		do j=-ssf_nmax,ssf_nmax
	
		k(ssf_ax1) = i*dk(ssf_ax1)
		k(ssf_ax2) = j*dk(ssf_ax2)

		c = 0.d0
		s = 0.d0	
		do n=1,np
			c = c + cos(dot_product(k,r(:,n)))
			s = s + sin(dot_product(k,r(:,n)))
		end do

		if (i.eq.0.and.j.eq.0) then
			c = 0.d0
			s = 0.d0
		end if

		posi = i + ssf_nmax + 1		
		posj = j + ssf_nmax + 1		
		ssf_hist(posi,posj) = ssf_hist(posi,posj) + c*c + s*s

		end do
	end do
	
	ssf(:,:) = ssf_hist(:,:)/(real(np)*hist_count)
	hist_count = hist_count + 1

end subroutine evaluate_properties_ssf

!===================================================================================
!Diffusion function calculated

!subroutine evaluate_properties_diffusion
!	use module_record
!	implicit none
!
!	integer          :: n
!
!	diffusion = 0
!	do n=1,np
!		diffusion(:)=diffusion(:)+(rtrue(:,n) - rinitial(:,n))**2
!	enddo
!	meandiffusion(:) = diffusion(:)/(2.d0*nd*np*(delta_t*iter/tplot))
!	!print*, 'Instantanous diffusion', diffusion
!	print*, 'time average of diffusion', meandiffusion
!
!	!Using the auto-correlation function
!	do n=1,np
!		diffusion(:) = v(:,n) - 0.5 * a(:,n) * delta_t
!	enddo
!
!end subroutine evaluate_properties_diffusion

!=============================================================================
!Evaluate viscometric data
subroutine evaluate_viscometrics
use calculated_properties_MD, only: Pxy
use shear_info_MD,            only: le_sd, le_sp, le_rp, le_sr
implicit none

	double precision :: eta              ! Shear viscosity
	double precision :: N1               ! First normal stress difference
	double precision :: N2               ! Second normal stress difference
	double precision :: P                ! Hydrostatic pressure
	integer          :: x,y,z            ! Notation (see below)

	x = le_sd                            ! Shear direction
	y = le_sp                            ! Shear plane
	z = le_rp                            ! Remaining plane

	eta = Pxy(x,y)/le_sr                 ! Shear stress / shear rate
	N1  = Pxy(x,x) - Pxy(y,y)            ! Normal stresses
	N2  = Pxy(y,y) - Pxy(z,z)            ! Normal stresses
	P   = -(1.0/3.0)*(Pxy(x,x) + &
	                  Pxy(y,y) + &
	                  Pxy(z,z)     )     ! Hydrostatic pressure

	eta = -1.d0*eta
	N1  = -1.d0*N1
	N2  = -1.d0*N2
	P   = -1.d0*P

	call viscometrics_io(eta,N1,N2,P)

end subroutine evaluate_viscometrics

!=============================================================================
!Calculate end-to-end time correlation function of FENE chain
subroutine etevtcf_calculate_parallel
	use module_record
	implicit none
	
	integer :: nbond,i,j
	integer :: chain, i_sub, j_sub, funcy
	double precision :: etev_prod, etev_prod_sum
	double precision :: etev2, etev2_sum
	double precision, dimension(nd) :: rij

	if (iter.eq.etevtcf_iter0) then	!Initialise end-to-end vectors at t_0

		allocate(etev(nchains,nd))
		allocate(etev_0(nchains,nd))
		etev_0 = 0.d0
		do i=1,np
			chain = monomer(i)%chainID
			i_sub = monomer(i)%subchainID
			funcy = monomer(i)%funcy
			do nbond=1,funcy
				j      = bond(nbond,i)
				j_sub  = monomer(j)%subchainID
				if (j_sub.lt.i_sub) cycle  !Avoid counting backwards
				rij(:) = r(:,j) - r(:,i)
				rij(:) = rij(:) - domain(:)*anint(rij(:)/domain(:))
				etev_0(chain,:) = etev_0(chain,:) + rij(:)
			end do
		end do
		call globalSumTwoDim(etev_0,nchains,nd)

	end if

	!--------------------------------------------------------------------
	!Calculate all end-to-end vectors for file output
	etev = 0.d0
	do i=1,np
		chain = monomer(i)%chainID
		i_sub = monomer(i)%subchainID
		funcy = monomer(i)%funcy
		do nbond=1,funcy
			j             = bond(nbond,i)     ! Molecule number j is nth bond to i
			j_sub         = monomer(j)%subchainID  ! Find subchain ID of mol j
			if (j_sub.lt.i_sub) cycle              ! Avoid counting backwards
			rij(:)        = r(:,j) - r(:,i)                     
			rij(:)        = rij(:) - domain(:)*anint(rij(:)/domain(:))
			etev(chain,:) = etev(chain,:) + rij(:)
		end do
	end do
	call globalSumTwoDim(etev,nchains,nd)
	!--------------------------------------------------------------------

	!--------------------------------------------------------------------
	!Running calculation for stdout...	
	etev_prod_sum	= 0.d0
	etev2_sum 		= 0.d0
	do chain=1,nchains
		etev_prod		  = dot_product(etev(chain,:),etev_0(chain,:))
		etev2			  = dot_product(etev(chain,:),etev(chain,:))			
		etev_prod_sum	  = etev_prod_sum + etev_prod
		etev2_sum		  = etev2_sum + etev2		
	end do
	etevtcf = etev_prod_sum/etev2_sum              !Sample counts cancel
	!-------------------------------------------------------------------

end subroutine etevtcf_calculate_parallel

!=============================================================================
!Calculate radius of gyration
subroutine r_gyration_calculate_parallel
use module_record
use linked_list
implicit none

	integer :: i,chainID
	double precision :: R_g2
	double precision, dimension(nd) :: rij
	double precision, dimension(:,:), allocatable :: r_cm

	allocate(r_cm(nchains,nd))

	!Calculate center of masses
	r_cm = 0.d0

	do i=1,np
		chainID = monomer(i)%chainID
		if (chainID .eq. 0) cycle
		r_cm(chainID,:) = r_cm(chainID,:) + rtrue(:,i) 
	end do

	call globalSumTwoDim(r_cm,nchains,nd)

	do chainID = 1,nchains
		r_cm(chainID,:) = r_cm(chainID,:)/nmonomers
	end do

	!Calculate R_g2
	R_g2 = 0.d0
	do i=1,np
		chainID = monomer(i)%chainID
		if (chainID.eq.0) cycle
		rij = rtrue(:,i) - r_cm(chainID,:)
		R_g2 = R_g2 + dot_product(rij,rij)
	end do
	call globalSum(R_g2)
	R_g2 = R_g2/(nmonomers*nchains)

	R_g  = R_g2**0.5d0

	deallocate(r_cm)

end subroutine r_gyration_calculate_parallel

!=============================================================================
!		RECORD MASS AT LOCATION IN SPACE
! Either by binning or taking slices on a plane, the mass in the molecular
! system is recorded and output
!===================================================================================

!Calculate averaged velocity components of each bin or slice with 2D slice in 
!ixyz = 1,2,3 or in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine mass_averaging(ixyz)
	use module_record
	use field_io, only : mass_slice_io,mass_bin_io,mass_bin_cpol_io
	use concentric_cylinders, only: cyl_mass
	implicit none

	integer							:: ixyz
	integer, save					:: average_count=-1

	average_count = average_count + 1
	call cumulative_mass(ixyz)
	if (average_count .eq. Nmass_ave) then
		average_count = 0
		!Determine bin size

		select case(ixyz)
		case(1:3)
			call mass_slice_io(ixyz)
			!Reset mass slice
			slice_mass = 0
		case(4)
			call mass_bin_io(volume_mass,'bins')
			!Reset mass slice
			volume_mass = 0
		case(5)
			call mass_bin_cpol_io(cyl_mass)
			cyl_mass = 0
		case default
			call error_abort("Error input for velocity averaging incorrect")
		end select

		!Collect mass for next step
		!call cumulative_mass(ixyz)

	endif


end subroutine mass_averaging

!-----------------------------------------------------------------------------------
!Add mass to running total, with 2D slice in ixyz = 1,2,3 or
!in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine cumulative_mass(ixyz)
	use module_record
	use linked_list
	use messenger, only: globalise
	use concentric_cylinders
	implicit none

	integer         				:: n, ixyz
	integer         				:: cbin
	integer,dimension(3)			:: ibin
	double precision				:: slicebinsize
	double precision,dimension(3) 	:: mbinsize 

	integer :: br, bt, bz
	real(kind(0.d0)) :: fluiddomain_cyl(3), rglob(3), rpol(3)

	select case(ixyz)
	!Mass measurement is a number of 2D slices through the domain
	case(1:3)
		slicebinsize = domain(ixyz) / nbins(ixyz)
		do n = 1, np    ! Loop over all particles
			!Assign to bins using integer division
			cbin = ceiling((r(ixyz,n)+halfdomain(ixyz))/slicebinsize)!Establish current bin
			if (cbin > nbins(ixyz)) cbin = nbins(ixyz) 		 !Prevents out of range values
			if (cbin < 1 ) cbin = 1        				 !Prevents out of range values
			slice_mass(cbin)= slice_mass(cbin)+1      			 !Add one to current bin
		enddo
	!Mass measurement for 3D bins throughout the domain
	case(4)
		mbinsize(:) = domain(:) / nbins(:) 
		do n = 1,np
			!Add up current volume mass densities
			ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb
			volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + 1
		enddo
	case(5)

		! Cylindrical polar coordinates                                       -
		fluiddomain_cyl(1) = r_io - r_oi
		fluiddomain_cyl(2) = 2.d0 * pi 
		fluiddomain_cyl(3) = domain(3) 
		mbinsize(:) = fluiddomain_cyl(:)/cpol_bins(:)				

		do n = 1,np

			!Find cyl pol coords, rr=0 at inner cylinder
			rglob    = globalise(r(:,n))
			rpol     = cpolariser(rglob)

			!Shift z component to be between 0 < r_z < domain(3),
			!      theta to 0 < theta < 2*pi (i.e. not -pi to pi),
			!      r to r_oi < r (< r_io) 
 			rpol(1)  = rpol(1) - r_oi
			rpol(2)  = modulo(rpol(2),2.d0*pi)
			rpol(3)  = r(3,n) + halfdomain(3) 

			!Add to cylindrical bins
			br = ceiling(rpol(1)/mbinsize(1)) 
			bt = ceiling(rpol(2)/mbinsize(2)) 
			bz = ceiling(rpol(3)/mbinsize(3)) + cpol_nhbz

			!Ignore cylinder molecules and stray liquid mols
			if ( br .gt. cpol_bins(1) ) cycle
 			if ( br .lt. 1 ) cycle
			if ( tag(n) .eq. cyl_teth_thermo_rotate ) cycle

			cyl_mass(br,bt,bz)  = cyl_mass(br,bt,bz)  + 1

		end do

	case default 
		call error_abort("Mass Binning Error")
	end select
 
end subroutine cumulative_mass

!===================================================================================
!		RECORD VELOCITY AT LOCATION IN SPACE
! Either by binning or taking slices on a plane, the velocity field in the molecular
! system is recorded and output
!===================================================================================

!Calculate averaged velocity components of each bin or slice with 2D slice in 
!ixyz = 1,2,3 or in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine velocity_averaging(ixyz)
	use module_record
	use field_io , only : velocity_slice_io,velocity_bin_io,velocity_bin_cpol_io
	use concentric_cylinders, only: cyl_mass, cyl_mom
	use linked_list
	implicit none

	integer				:: ixyz, n
	integer,dimension(3):: ib
	integer, save		:: average_count=-1
	double precision,dimension(3) 	:: Vbinsize, temp

	average_count = average_count + 1
	call cumulative_velocity(ixyz)

	!Save streaming velocity for temperature averages
	if (temperature_outflag .ne. 0 .and. peculiar_flag .ne. 0) then

		! NOTE THE peculiar_flag CAN BE SET TO 0 AND
		! THE UNBIAS TEMPERATURE USUALLY CALCULATED FROM
		! T_{unbias} = (1/3N) * \sum_i^N m_i * (vi - u)*(vi - u)
		! CAN INSTEAD BE CALCULATED FROM:
		! T_{unbias} = (1/3N) * \sum_i^N m_i * vi*vi - u^2/3
		stop "Peculiar momentum functionality removed -- please calculate using T_{unbias} = (1/3N) * \sum_i^N m_i*vi*vi - u^2/3"

		!Determine bin size
		!Vbinsize(:) = domain(:) / nbins(:)

		!Get instantanous temperature field included swapped halos
		!sm = volume_mass
		!sv = volume_momentum
		!call rswaphalos(sm,nbinso(1),nbinso(2),nbinso(3),1)
		!call rswaphalos(sv,nbinso(1),nbinso(2),nbinso(3),3)

		do n=1,np
			!Save streaming velocity per molecule
			ib(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
			!U(:,n) =  sv(ib(1),ib(2),ib(3),:) / sm(ib(1),ib(2),ib(3))
			U(:,n) =  volume_momentum(ib(1),ib(2),ib(3),:) / volume_mass(ib(1),ib(2),ib(3))
		enddo

	endif

	if (average_count .eq. Nvel_ave) then
		average_count = 0

		select case(ixyz)
			case(1:3)
				call velocity_slice_io(ixyz)
				!Reset velocity slice
				slice_mass = 0
				slice_momentum  = 0.d0
			case(4)
				call velocity_bin_io(volume_mass,volume_momentum,'bins')
				!Reset velocity bins
				volume_mass = 0
				volume_momentum = 0.d0
			case(5)
				call velocity_bin_cpol_io(cyl_mass,cyl_mom)
				cyl_mass = 0
				cyl_mom = 0.d0
			case default
				call error_abort("Error input for velocity averaging incorrect")
			end select

			!Collect velocities for next step
			!call cumulative_velocity(ixyz)

	endif

end subroutine velocity_averaging

subroutine evaluate_U
	use module_record
	use linked_list
	implicit none

	integer				:: n
	integer,dimension(3):: ib
	double precision,dimension(3) 	:: Vbinsize

	integer, dimension(:,:,:), allocatable :: mbin
	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: vbin

	integer :: x,y,z,c

	Vbinsize(:) = domain(:) / nbins(:)
	
	x = nbins(1) + 2*nhb(1)
	y = nbins(2) + 2*nhb(2)
	z = nbins(3) + 2*nhb(3)
	c = 3
	
	allocate(vbin(x,y,z,c))
	allocate(mbin(x,y,z))
	
	vbin = 0.d0
	mbin = 0
		
	do n = 1, np

		! Get bin
		ib(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
		! Add v, m to bin
		vbin(ib(1),ib(2),ib(3),:) = vbin(ib(1),ib(2),ib(3),:) + v(:,n)
		mbin(ib(1),ib(2),ib(3)) = mbin(ib(1),ib(2),ib(3)) + 1 

	end do

	do n = 1, np

		! Get bin
		ib(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
		! Get U from vbin / mbin
		U(:,n) = vbin(ib(1),ib(2),ib(3),:)/real(mbin(ib(1),ib(2),ib(3)),kind(0.d0))
	
	end do
	
	deallocate(vbin)
	deallocate(mbin)

end subroutine evaluate_U

!-----------------------------------------------------------------------------------
!Add velocities to running total, with 2D slice in ixyz = 1,2,3 or
!in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine cumulative_velocity(ixyz)
	use module_record
	use concentric_cylinders
	use messenger, only: globalise, localise
	use linked_list
	implicit none

	integer							:: n,ixyz
	integer         				:: cbin
	integer		,dimension(3)		:: ibin
	double precision				:: slicebinsize
	double precision,dimension(3) 	:: Vbinsize 
	
	integer :: br, bt, bz
	real(kind(0.d0)) :: fluiddomain_cyl(3), rglob(3), rpol(3), vpol(3)

	!In case someone wants to record velocity in a simulation without sliding walls!?!?
	if (ensemble .ne. tag_move) then
		allocate(slidev(3,np)); slidev = 0.d0
	endif

	select case(ixyz)
	!Velocity measurement is a number of 2D slices through the domain
	case(1:3)

		slicebinsize = domain(ixyz) / nbins(ixyz)

		do n = 1, np    ! Loop over all particles
			!Assign to bins using integer division
			cbin = ceiling((r(ixyz,n)+halfdomain(ixyz))/slicebinsize)!Establish current bin
			if (cbin > nbins(ixyz)) cbin = nbins(ixyz) 		 !Prevents out of range values
			if (cbin < 1 ) cbin = 1        				 !Prevents out of range values
			slice_mass(cbin)= slice_mass(cbin)+1      			 !Add one to current bin
			slice_momentum(cbin,:) = slice_momentum(cbin,:)+v(:,n) + slidev(:,n) 	 !Add streamwise velocity to current bin
		enddo

	!Velocity measurement for 3D bins throughout the domain
	case(4)

		!Determine bin size
		Vbinsize(:) = domain(:) / nbins(:)

		!Reset Control Volume momentum 
		do n = 1,np
			!Add up current volume mass and momentum densities
			ibin(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
			volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + 1
			volume_momentum(ibin(1),ibin(2),ibin(3),:) = volume_momentum(ibin(1),ibin(2),ibin(3),:) & 
														+ v(:,n) + slidev(:,n)
		enddo

	case(5)

		! Cylindrical polar coordinates                                       -
		fluiddomain_cyl(1) = r_io - r_oi
		fluiddomain_cyl(2) = 2.d0 * pi 
		fluiddomain_cyl(3) = domain(3) 
		Vbinsize(:) = fluiddomain_cyl(:)/cpol_bins(:)				

		do n = 1,np

			!Find cyl pol coords, rr=0 at inner cylinder
			rglob    = globalise(r(:,n))
			rpol     = cpolariser(rglob)
			vpol     = cpolarisev(v(:,n),rpol(2))

			!Shift z component to be between 0 < r_z < domain(3),
			!      theta to 0 < theta < 2*pi (i.e. not -pi to pi),
			!      r to r_oi < r (< r_io) 
 			rpol(1)  = rpol(1) - r_oi
			rpol(2)  = modulo(rpol(2),2.d0*pi)
			rpol(3)  = r(3,n) + halfdomain(3) 

			!Add to cylindrical bins
			br = ceiling(rpol(1)/Vbinsize(1)) 
			bt = ceiling(rpol(2)/Vbinsize(2)) 
			bz = ceiling(rpol(3)/Vbinsize(3)) + cpol_nhbz

			!Ignore cylinder molecules and stray liquid mols
			if ( br .gt. cpol_bins(1) ) cycle
 			if ( br .lt. 1 ) cycle
			if ( tag(n) .eq. cyl_teth_thermo_rotate ) cycle

			cyl_mass(br,bt,bz)  = cyl_mass(br,bt,bz)  + 1
			cyl_mom(br,bt,bz,:) = cyl_mom(br,bt,bz,:) + vpol

		enddo
	
	case default 
		call error_abort("Velocity Binning Error")
	end select

	if (ensemble .ne. tag_move) then
		deallocate(slidev)
	endif

	 
end subroutine cumulative_velocity

!===================================================================================
!		RECORD TEMPERATURE AT LOCATION IN SPACE
! Either by binning or taking slices on a plane, the temperature field in the molecular
! system is recorded and output
!===================================================================================

!Calculate averaged temperature components of each bin or slice with 2D slice in 
!ixyz = 1,2,3 or in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine temperature_averaging(ixyz)
	use module_record
	use field_io, only : temperature_slice_io,temperature_bin_io,temperature_bin_cpol_io
	use linked_list
	use concentric_cylinders
	implicit none

	integer				:: ixyz
	integer, save		:: average_count=-1
	
	average_count = average_count + 1
	call cumulative_temperature(ixyz)
	if (average_count .eq. NTemp_ave) then
		average_count = 0

		select case(ixyz)
		case(1:3)
			call temperature_slice_io(ixyz)
			!Reset temperature slice
			if (velocity_outflag .ne. ixyz) slice_mass = 0
			slice_temperature  = 0.d0
		case(4)
			call temperature_bin_io(volume_mass,volume_temperature,'bins')
			!Reset temperature bins
			if (velocity_outflag .ne. 4) volume_mass = 0
			volume_temperature = 0.d0
		case(5)
			call temperature_bin_cpol_io(cyl_mass,cyl_KE)
			if (velocity_outflag .ne. 5) cyl_mass = 0
			cyl_KE = 0.d0
		case default
			stop "Error input for temperature averaging incorrect"
		end select

		!Collect temperature for next step
		!call cumulative_temperature(ixyz)

	endif

end subroutine temperature_averaging

!-----------------------------------------------------------------------------------
!Add velocities to running total, with 2D slice in ixyz = 1,2,3 or
!in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine cumulative_temperature(ixyz)
	use module_record
	use concentric_cylinders
	use messenger, only: globalise, localise
	use linked_list
	implicit none

	integer							:: n,ixyz
	integer         				:: cbin
	integer		,dimension(3)		:: ibin
	double precision				:: slicebinsize
	double precision,dimension(3) 	:: Tbinsize 

	integer :: br, bt, bz
	real(kind(0.d0)) :: fluiddomain_cyl(3), rglob(3), rpol(3), vpol(3)

	!In case someone wants to record velocity in a simulation without sliding walls!?!?
	if (ensemble .ne. tag_move) then
		allocate(slidev(3,np)); slidev = 0.d0
	endif

	select case(ixyz)
	!temperature measurement is a number of 2D slices through the domain
	case(1:3)
		slicebinsize = domain(ixyz) / nbins(ixyz)
		do n = 1, np    ! Loop over all particles
			!Assign to bins using integer division
			cbin = ceiling((r(ixyz,n)+halfdomain(ixyz))/slicebinsize)!Establish current bin
			if (cbin > nbins(ixyz)) cbin = nbins(ixyz) 		 !Prevents out of range values
			if (cbin < 1 ) cbin = 1        				 !Prevents out of range values
			if (velocity_outflag .ne. ixyz) & 
			slice_mass(cbin)= slice_mass(cbin)+1      			 !Add one to current bin
			slice_temperature(cbin) = slice_temperature(cbin) & 
					+ dot_product((v(:,n)+slidev(:,n)),(v(:,n)+slidev(:,n))) 	 !Add streamwise temperature to current bin
		enddo

	!Temperature measurement for 3D bins throughout the domain
	case(4)

		!Determine bin size
		Tbinsize(:) = domain(:) / nbins(:)
		
		!call evaluate_U 

		!Reset Control Volume momentum 
		do n = 1,np
			!Add up current volume mass and temperature densities
			ibin(:) = ceiling((r(:,n)+halfdomain(:))/Tbinsize(:)) + nhb
			if (velocity_outflag .ne. 4) & 
			volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + 1
			!Note - the streaming term is removed but includes sliding so this must be added back on
			if (peculiar_flag .eq. 0) then
				!if (mod(iter,1000) .eq. 0 	.and. &
				!	ibin(1) .eq. 4 			.and. &
				!	ibin(2) .eq. 2 			.and. &
				!	ibin(3) .eq. 4				) then
				!	print*, iter, ibin, dot_product(v(:,n),v(:,n)), v(:,n)
				!endif
				volume_temperature(ibin(1),ibin(2),ibin(3)) = volume_temperature(ibin(1),ibin(2),ibin(3)) & 
										+ dot_product((v(:,n)+slidev(:,n)),(v(:,n)+slidev(:,n)))
			else
				volume_temperature(ibin(1),ibin(2),ibin(3)) = volume_temperature(ibin(1),ibin(2),ibin(3)) & 
										+ dot_product((v(:,n)-U(:,n)+slidev(:,n)), & 
													  (v(:,n)-U(:,n)+slidev(:,n)))
				!write(958,'(2i8,5f10.5)'),iter,n,r(2,n),U(1,n),v(1,n),dot_product(v(:,n),v(:,n)),dot_product((v(:,n)-U(:,n)+slidev(:,n)), & 
				!									  (v(:,n)-U(:,n)+slidev(:,n)))
			endif

		enddo

	case(5)
	
		! Cylindrical polar coordinates                                       -
		fluiddomain_cyl(1) = r_io - r_oi
		fluiddomain_cyl(2) = 2.d0 * pi 
		fluiddomain_cyl(3) = domain(3) 
		Tbinsize(:) = fluiddomain_cyl(:)/cpol_bins(:)				

		do n = 1,np

			!Find cyl pol coords, rr=0 at inner cylinder
			rglob    = globalise(r(:,n))
			rpol     = cpolariser(rglob)
			!v^2 is independent of coordinate system
			!vpol     = cpolarisev(v(:,n),rpol(2))

			!Shift z component to be between 0 < r_z < domain(3),
			!      theta to 0 < theta < 2*pi (i.e. not -pi to pi),
			!      r to r_oi < r (< r_io) 
 			rpol(1)  = rpol(1) - r_oi
			rpol(2)  = modulo(rpol(2),2.d0*pi)
			rpol(3)  = r(3,n) + halfdomain(3) 

			!Add to cylindrical bins
			br = ceiling(rpol(1)/Tbinsize(1)) 
			bt = ceiling(rpol(2)/Tbinsize(2)) 
			bz = ceiling(rpol(3)/Tbinsize(3)) + cpol_nhbz

			!Ignore cylinder molecules and stray liquid mols
			if ( br .gt. cpol_bins(1) ) cycle
 			if ( br .lt. 1 ) cycle
			if ( tag(n) .eq. cyl_teth_thermo_rotate ) cycle

			cyl_mass(br,bt,bz) = cyl_mass(br,bt,bz) + 1
			cyl_KE(br,bt,bz) = cyl_KE(br,bt,bz) + dot_product(v(:,n),v(:,n))

		enddo

	case default 
		stop "Temperature Binning Error"
	end select

	if (ensemble .ne. tag_move) then
		deallocate(slidev)
	endif
	 
end subroutine cumulative_temperature

!===================================================================================
!		RECORD PRESSURE AT LOCATION IN SPACE
! Either by Volume Average style binning or virial expression for the whole domain
!===================================================================================

subroutine pressure_averaging(ixyz)
	use field_io, only : virial_stress_io,VA_stress_io,VA_stress_cpol_io
	use module_record
	implicit none

	integer			:: ixyz
	integer, save	:: sample_count, average_count

	sample_count = sample_count + 1
	call cumulative_pressure(ixyz,sample_count)
	if (sample_count .eq. Nstress_ave) then
		sample_count = 0
		Pxy = Pxy/dble(Nstress_ave)
		Pxyzero = Pxy		!Update Pxy(0) value

		select case(ixyz)
		case(1)
		!FULL DOMAIN VIRIAL STRESS CALCULATION
			!print'(a,10f12.5)', 'cumulative stress', Pxy,(Pxy(1,1)+Pxy(2,2)+Pxy(3,3))/3.d0
			call virial_stress_io
			Pxy = 0.d0
		case(2)
		!VA STRESS CALCULATION
			call VA_stress_io
			Pxybin = 0.d0
			vvbin  = 0.d0
			rfbin  = 0.d0
		case(3)
			call VA_stress_cpol_io
			Pxybin = 0.d0
			vvbin  = 0.d0
			rfbin  = 0.d0
		case default 
			call error_abort("Average Pressure Binning Error")
		end select
		if(viscosity_outflag .eq. 1) then
			average_count = average_count+1
			if (average_count .eq. Nvisc_ave) then
				call viscosity_io
				average_count = 0
			endif
		endif
	endif

end subroutine pressure_averaging

!===================================================================================
!Add pressure_tensor to running total

subroutine cumulative_pressure(ixyz,sample_count)
	use module_record
	use shear_info_MD, only: le_sp, le_sr, le_sd
	implicit none

	integer								:: sample_count,n,ixyz,jxyz,kxyz
	double precision, dimension(3)		:: velvect
	double precision, dimension(3)      :: rglob
	double precision, dimension(3,3)	:: Pxytemp


	Pxytemp = 0.d0
	Pxymol = 0.d0

	!Factor of 2 as every interaction calculated
	rfmol = rfmol / 2.d0

	select case(ixyz)
	case(1)
		!FULL DOMAIN VIRIAL STRESS CALCULATION
		do n = 1, np    ! Calculate pressure tensor for all molecules

			select case(integration_algorithm)
			case(leap_frog_verlet)
				!Calculate velocity at time t (v is at t+0.5delta_t due to use of verlet algorithm)
				velvect(:) = v(:,n) + 0.5d0 * a(:,n) * delta_t
			case(velocity_verlet)     
				!Velocity is already at time t for Velocity Verlet algorithm                              
				velvect(:) = v(:,n)
			end select

			!Adjustment for Lees Edwards sliding boundries
			if (any(periodic.gt.1)) then
				rglob(1) = r(1,n)-(halfdomain(1)*(npx-1))+domain(1)*(iblock-1)
				rglob(2) = r(2,n)-(halfdomain(2)*(npy-1))+domain(2)*(jblock-1)
				rglob(3) = r(3,n)-(halfdomain(3)*(npz-1))+domain(3)*(kblock-1)
				velvect(le_sd) = velvect(le_sd) - rglob(le_sp)*le_sr 
				!velvect(:) = velvect(:) - U(:,n)
			end if

			do jxyz = 1,3
			do kxyz = 1,3
				Pxymol(n,kxyz,jxyz) = Pxymol(n,kxyz,jxyz)    &
						      		+ velvect(kxyz)*velvect(jxyz) &
						      		+ rfmol(n,kxyz,jxyz)
			enddo
			enddo

			!Sum of molecules to obtain pressure tensor for entire domain
			Pxytemp(:,:) = Pxytemp(:,:) + Pxymol(n,:,:)

		enddo

		!Sum pressure tensor over all processors -- Should this be Pxytemp???
		call globalSumVect(Pxytemp(:,1), nd)	
		call globalSumVect(Pxytemp(:,2), nd) 	
		call globalSumVect(Pxytemp(:,3), nd)

		!Divide sum of stress by volume -- Should this include sample_count - divide by Nstress above??
		Pxytemp = Pxytemp / volume

		!Add current pressure to cumulative average
		Pxy = Pxy + Pxytemp

		!Reset position force tensor before calculation
		rfmol = 0.d0
	case(2)
		!VOLUME AVERAGE STRESS CALCULATION
		!Calculate Position (x) force for configurational part of stress tensor
		call  simulation_compute_rfbins(1,nbins(1)+2,1,nbins(2)+2,1,nbins(3)+2)
		!Calculate velocity (x) velocity for kinetic part of stress tensor
		if (nbins(1) .eq. ncells(1) .and. & 
		    nbins(2) .eq. ncells(2) .and. & 
		    nbins(3) .eq. ncells(3)) then
			call  simulation_compute_kinetic_VA_cells(2,nbins(1)+1,2,nbins(2)+1,2,nbins(3)+1)
		else
			call  simulation_compute_kinetic_VA(2,nbins(1)+1,2,nbins(2)+1,2,nbins(3)+1)
		endif

		!Add results to cumulative total
		Pxybin(:,:,:,:,:) =     vvbin(:,:,:,:,:) 				& 
						      + rfbin(  2:nbins(1)+1, 			& 
							      		2:nbins(2)+1, 			& 
							      		2:nbins(3)+1,:,:)/2.d0

		!NOTE: REBUILD AT (mod(iter+1,tplot) .eq. 0) WHEN RECORD AFTER FORCE CALCULATION
		!Reset bin force tensor before next calculation
	  	!rfbin = 0.d0; vvbin = 0.d0

		!Calculate Virial from sum of Volume Averaged Stresses
		do kxyz = 1,3
		do jxyz = 1,3
			Pxytemp(kxyz,jxyz) = sum(Pxybin(:,:,:,kxyz,jxyz))
		enddo
		enddo

		!Divide sum of stress by volume and number of samples for each bin
		Pxytemp = Pxytemp / (volume*sample_count)

		!Add current pressure to cumulative average
		Pxy = Pxy + Pxytemp
		!print'(a,i8,3f18.4)','Sum of all VA ',iter,(sum(rfbin(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,1,1)) + & 
		!										    sum(rfbin(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,2,2)) + & 
		!											sum(rfbin(2:nbins(1)+1,2:nbins(2)+1,2:nbins(3)+1,3,3))) & 
		!												/ (6.d0 * volume*Nstress_ave), & 
		!										   (sum(vvbin(:,:,:,1,1)) + & 
		!										    sum(vvbin(:,:,:,2,2)) + & 
		!											sum(vvbin(:,:,:,3,3)))/ (3.d0 * volume*Nstress_ave), & 
		!											(Pxy(1,1) + Pxy(2,2) + Pxy(3,3))/3.d0
	case(3)

		call simulation_compute_rfbins_cpol(1,nbins(1)+2,1,nbins(2)+2,1,nbins(3)+2)
		call simulation_compute_kinetic_VA_cpol(2,nbins(1)+1,2,nbins(2)+1,2,nbins(3)+1)
		Pxybin = vvbin + rfbin/2.d0

	case default 
		call error_abort("Cumulative Pressure Averaging Error")
	end select

	!Calculate correrlation between current value of Pxy and intial value of sample
	!Average 3 components of Pxycorrel to improve statistics
	if(viscosity_outflag .eq. 1) then
		if(sample_count .ne. 0) then
			Pxycorrel(sample_count) = Pxycorrel(sample_count) + Pxytemp(2,3)*Pxyzero(2,3)
			Pxycorrel(sample_count) = Pxycorrel(sample_count) + Pxytemp(1,3)*Pxyzero(1,3)
			Pxycorrel(sample_count) = Pxycorrel(sample_count) + Pxytemp(1,2)*Pxyzero(1,2)
		endif
	endif


end subroutine cumulative_pressure

!----------------------------------------------------------------------------------
!Compute kinetic part of stress tensor

subroutine simulation_compute_kinetic_VA(imin,imax,jmin,jmax,kmin,kmax)
	use module_record
	use physical_constants_MD
	implicit none

	integer								:: imin, jmin, kmin, imax, jmax, kmax
	integer         					:: n, ixyz,jxyz
	integer 							:: ibin, jbin, kbin
	double precision, dimension(3)		:: VAbinsize, velvect

	!vvbin = 0.d0

	!Determine bin size
	VAbinsize(:) = domain(:) / nbins(:)

	! Add kinetic part of pressure tensor for all molecules
	do n = 1, np

		!Assign to bins using integer division
		ibin = ceiling((r(1,n)+halfdomain(1))/VAbinsize(1))	!Establish current bin
		if (ibin .gt. imax) cycle ! ibin = maxbin		!Prevents out of range values
		if (ibin .lt. imin) cycle ! ibin = minbin		!Prevents out of range values
		jbin = ceiling((r(2,n)+halfdomain(2))/VAbinsize(2)) 	!Establish current bin
		if (jbin .gt. jmax) cycle ! jbin = maxbin 		!Prevents out of range values
		if (jbin .lt. jmin) cycle ! jbin = minbin		!Prevents out of range values
		kbin = ceiling((r(3,n)+halfdomain(3))/VAbinsize(3)) 	!Establish current bin
		if (kbin .gt. kmax) cycle ! kbin = maxbin		!Prevents out of range values
		if (kbin .lt. kmin) cycle ! kbin = minbin		!Prevents out of range values

		select case(integration_algorithm)
		case(leap_frog_verlet)
			!Calculate velocity at time t (v is at t+0.5delta_t due to use of verlet algorithm)
			velvect(:) = v(:,n) + 0.5d0 * a(:,n) * delta_t
			!Velocity is already at time t for Velocity Verlet algorithm
		case(velocity_verlet)                                   
			velvect(:) = v(:,n)
		end select

		do ixyz = 1,3
		do jxyz = 1,3
			vvbin(ibin,jbin,kbin,ixyz,jxyz) = vvbin(ibin,jbin,kbin,ixyz,jxyz)	&
					       		  + velvect(ixyz) * velvect(jxyz)
		enddo
		enddo

	enddo

end subroutine simulation_compute_kinetic_VA

subroutine simulation_compute_kinetic_VA_cpol(imin,imax,jmin,jmax,kmin,kmax)
	use module_record
	use physical_constants_MD
	use concentric_cylinders
	use messenger, only: globalise
	use librarymod, only: cpolariser, cpolariseT, outerprod
	implicit none

	integer, intent(in) :: imin, jmin, kmin, imax, jmax, kmax
	integer :: n, ixyz,jxyz
	integer :: br, bt, bz
	real(kind(0.d0)) :: VAbinsize(3), velvect(3)
	real(kind(0.d0)) :: ripol(3), vvpol(3,3)

	!Determine bin sizes
	VAbinsize(1) = (r_io - r_oi) / cpol_bins(1)
	VAbinsize(2) = 2.d0*pi       / cpol_bins(2)
	VAbinsize(3) = domain(3)     / cpol_bins(3)

	! Add kinetic part of pressure tensor for all molecules
	do n = 1, np

		!Assign to bins using integer division
		br = ceiling((r(1,n)+halfdomain(1))/VAbinsize(1))	!Establish current bin
		bt = ceiling((r(2,n)+halfdomain(2))/VAbinsize(2)) 	!Establish current bin
		bz = ceiling((r(3,n)+halfdomain(3))/VAbinsize(3)) 	!Establish current bin

		select case(integration_algorithm)
		case(leap_frog_verlet)
			!Calculate velocity at time t (v is at t+0.5delta_t due to use of verlet algorithm)
			velvect(:) = v(:,n) + 0.5d0 * a(:,n) * delta_t
			!Velocity is already at time t for Velocity Verlet algorithm
		case(velocity_verlet)                                   
			velvect(:) = v(:,n)
		end select

		ripol = cpolariser(globalise(r(:,n)))
		vvpol = cpolariseT(outerprod(velvect,velvect),ripol(2))

		! Binning conventions 
		ripol(1) = ripol(1) - r_oi
		ripol(2) = modulo(ripol(2),2.d0*pi)
		ripol(3) = r(3,n) + halfdomain(3) 

		! Cylindrical bins
		br = ceiling(ripol(1)/VAbinsize(1)) 
		bt = ceiling(ripol(2)/VAbinsize(2)) 
		bz = ceiling(ripol(3)/VAbinsize(3)) + cpol_nhbz

		!Ignore molecules not in fluid region
		if (br .gt. cpol_bins(1)) cycle
		if (br .lt. 1)            cycle
		if (bt .gt. cpol_bins(2)) cycle
		if (bt .lt. 1)            cycle
		if (bz .gt. cpol_bins(3)) cycle
		if (bz .lt. 1)            cycle
		
		vvbin(br,bt,bz,:,:) = vvbin(br,bt,bz,:,:) + vvpol(:,:)

	enddo

end subroutine simulation_compute_kinetic_VA_cpol

!----------------------------------------------------------------------------------
!Compute kinetic part of stress tensor ONLY IF BINSIZE = CELLSIZE

subroutine simulation_compute_kinetic_VA_cells(imin,imax,jmin,jmax,kmin,kmax)
	use module_record
	use physical_constants_MD
	use linked_list
	implicit none

	integer,intent(in)				:: imin, jmin, kmin, imax, jmax, kmax
	integer                         :: i, ixyz, jxyz   !Define dummy index
	integer 						:: ibin, jbin, kbin,icell, jcell, kcell
	integer                         :: cellnp, molnoi
	double precision, dimension(3)	:: velvect
	type(node), pointer 	        :: oldi, currenti

	!vvbin = 0.d0

	! Add kinetic part of pressure tensor for all cells
	do kcell=kmin, kmax
	do jcell=jmin, jmax
	do icell=imin, imax

		!ASSUME Cell same size as bins
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list
		ibin = icell-1; jbin = jcell-1; kbin = kcell-1	

		do i = 1,cellnp					!Step through each particle in list 
			molnoi = oldi%molno 	 	!Number of molecule

			select case(integration_algorithm)
			case(leap_frog_verlet)
				!Calculate velocity at time t (v is at t+0.5delta_t due to use of verlet algorithm)
				velvect(:) = v(:,molnoi) + 0.5d0 * a(:,molnoi) * delta_t
				!Velocity is already at time t for Velocity Verlet algorithm
			case(velocity_verlet)                                   
				velvect(:) = v(:,molnoi)
			end select

			do ixyz = 1,3
			do jxyz = 1,3
				vvbin(ibin,jbin,kbin,ixyz,jxyz) = vvbin(ibin,jbin,kbin,ixyz,jxyz) &
								  + velvect(ixyz) * velvect(jxyz)
			enddo
			enddo

			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo
	enddo
	enddo

	nullify(oldi)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required

end subroutine simulation_compute_kinetic_VA_cells

!=================================================================================
! CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV
! CV 																			CV
! CV   _____             _             _   _   _       _                      	CV
! CV  /  __ \           | |           | | | | | |     | |                     	CV
! CV  | /  \/ ___  _ __ | |_ _ __ ___ | | | | | | ___ | |_   _ _ __ ___   ___ 	CV
! CV  | |    / _ \| '_ \| __| '__/ _ \| | | | | |/ _ \| | | | | '_ ` _ \ / _ \	CV
! CV  | \__/\ (_) | | | | |_| | | (_) | | \ \_/ / (_) | | |_| | | | | | |  __/	CV
! CV  \_____/\___/|_| |_|\__|_|  \___/|_|  \___/ \___/|_|\__,_|_| |_| |_|\___|	CV
! CV 																		    CV
! CV 			Record Fluxes accross surfaces of Control Volumes				CV
! CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV
!=================================================================================
!				CONTROL VOLUME FORMULATION OF MASS AND STRESS
! Based on the control volume formulation of mass and momentum, it is possible to
! define change in properties in terms of surface fluxes over a volume in space
! Fluxes need to be called every timestep by setting CV_CONSERVE to 1 in input
!=================================================================================

!-----------------------------------------------------------------------------------
! Control Volume mass continuity
!-----------------------------------------------------------------------------------

subroutine mass_flux_averaging(ixyz)
	!use field_io, only : mass_flux_io
	use module_record
	use CV_objects, only : CVcheck_mass, CV_debug
	implicit none

	integer			:: ixyz
	integer, save	:: sample_count

	!Only average if mass averaging turned on
	if (ixyz .eq. 0) return

	call cumulative_mass_flux
	sample_count = sample_count + 1
	if (sample_count .eq. Nmflux_ave) then
		if (CV_debug) call CVcheck_mass%check_error(2,nbins(1)+1,2,nbins(2)+1,2,nbins(3)+1,iter,irank)
		call mass_flux_io
		sample_count = 0
		mass_flux = 0
		call mass_snapshot
	endif

end subroutine mass_flux_averaging


!===================================================================================
! Mass Flux over a surface of a bin -- only single crossings assumed

!subroutine cumulative_mass_flux
!	use module_record
!	implicit none

!	integer							:: ixyz, n
!	integer		,dimension(3)		:: ibin1,ibin2,crossplane
!	double precision,dimension(3)	:: mbinsize, ri1, ri2

	!Determine bin size
!	mbinsize(:) = domain(:) / nbins(:)

!	do n = 1,np

!		ri1(:) = r(:,n) 							!Molecule i at time t
!		ri2(:) = r(:,n)	-delta_t*v(:,n)				!Molecule i at time t-dt

		!Assign to bins before and after using integer division
!		ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb
!		ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb

		!Replace Signum function with this functions which gives a
		!check for plane crossing and the correct sign 
!		crossplane(:) =  ibin1(:) - ibin2(:)

!		if (sum(abs(crossplane(:))) .ne. 0) then

			!Find which direction the surface is crossed
			!For simplicity, if more than one surface has been crossed surface fluxes of intermediate cells
			!are not included. This assumption => more reasonable as Delta_t => 0 or Delta_r => ∞
			!imaxloc = maxloc(abs(crossplane))
!			ixyz = imaxloc(abs(crossplane))

			!Add mass flux to the new bin surface count and take from the old
!			mass_flux(ibin1(1),ibin1(2),ibin1(3),ixyz+3*heaviside(-dble(crossplane(ixyz)))) = & 
!				mass_flux(ibin1(1),ibin1(2),ibin1(3),ixyz+3*heaviside(-dble(crossplane(ixyz)))) & 
!					+ abs(crossplane(ixyz))
!			mass_flux(ibin2(1),ibin2(2),ibin2(3),ixyz+3*heaviside(dble(crossplane(ixyz)))) = & 
!				mass_flux(ibin2(1),ibin2(2),ibin2(3),ixyz+3*heaviside(dble(crossplane(ixyz)))) &
!					- abs(crossplane(ixyz))
!		endif

!	enddo

!end subroutine cumulative_mass_flux

!===================================================================================
! Mass Flux over a surface of a bin
! Includes all intermediate bins

subroutine cumulative_mass_flux
	use module_record
	implicit none

	integer							:: ixyz,jxyz,i,j,k,n
	integer							:: planeno
	integer							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	double precision				:: crossplane,rplane,shift
	double precision,dimension(3)	:: mbinsize,velvect,crossface
	double precision,dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	do n = 1,np

		ri1(:) = r(:,n) 							!Molecule i at time t
		ri2(:) = r(:,n)	- delta_t*v(:,n)			!Molecule i at time t-dt
		ri12   = ri1 - ri2							!Molecule i trajectory between t-dt and t
		where (ri12 .eq. 0.d0) ri12 = 0.000001d0

		!Assign to bins before and after using integer division
		ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb(:)
		ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb(:)

		!Replace Signum function with this functions which gives a
		!check for plane crossing and the correct sign 
		crossface(:) =  ibin1(:) - ibin2(:)

		if (sum(abs(crossface(:))) .ne. 0) then

			do i = ibin1(1),ibin2(1),sign(1,ibin2(1)-ibin1(1))
			do j = ibin1(2),ibin2(2),sign(1,ibin2(2)-ibin1(2))
			do k = ibin1(3),ibin2(3),sign(1,ibin2(3)-ibin1(3))

				cbin(1) = i; cbin(2) = j; cbin(3) = k

				bintop(:) = (cbin(:)-1*nhb(:)  )*mbinsize(:)-halfdomain(:)
				binbot(:) = (cbin(:)-1*nhb(:)-1)*mbinsize(:)-halfdomain(:)

				!Calculate the plane intersect of trajectory with surfaces of the cube
				Pxt=(/ 			bintop(1), 		     & 
						ri1(2)+(ri12(2)/ri12(1))*(bintop(1)-ri1(1)), & 
						ri1(3)+(ri12(3)/ri12(1))*(bintop(1)-ri1(1))  	/)
				Pxb=(/ 			binbot(1), 		     & 
						ri1(2)+(ri12(2)/ri12(1))*(binbot(1)-ri1(1)), & 
						ri1(3)+(ri12(3)/ri12(1))*(binbot(1)-ri1(1))  	/)
				Pyt=(/	ri1(1)+(ri12(1)/ri12(2))*(bintop(2)-ri1(2)), & 
							bintop(2), 		     & 
						ri1(3)+(ri12(3)/ri12(2))*(bintop(2)-ri1(2))  	/)
				Pyb=(/	ri1(1)+(ri12(1)/ri12(2))*(binbot(2)-ri1(2)), &
							binbot(2), 		     & 
						ri1(3)+(ri12(3)/ri12(2))*(binbot(2)-ri1(2))  	/)
				Pzt=(/	ri1(1)+(ri12(1)/ri12(3))*(bintop(3)-ri1(3)), & 
						ri1(2)+(ri12(2)/ri12(3))*(bintop(3)-ri1(3)), &
							bintop(3) 			/)
				Pzb=(/	ri1(1)+(ri12(1)/ri12(3))*(binbot(3)-ri1(3)), &
						ri1(2)+(ri12(2)/ri12(3))*(binbot(3)-ri1(3)), & 
							binbot(3) 			/)

				onfacexb =0.5d0*(sign(1.d0,binbot(1) - ri2(1)) 	 & 
						       - sign(1.d0,binbot(1) - ri1(1)))* &
								(heaviside(bintop(2) - Pxb(2)) 	 &
						       - heaviside(binbot(2) - Pxb(2)))* &
								(heaviside(bintop(3) - Pxb(3)) 	 &
						       - heaviside(binbot(3) - Pxb(3)))
				onfaceyb =0.5d0*(sign(1.d0,binbot(2) - ri2(2))   &
						       - sign(1.d0,binbot(2) - ri1(2)))* &
								(heaviside(bintop(1) - Pyb(1))   &
						       - heaviside(binbot(1) - Pyb(1)))* &
								(heaviside(bintop(3) - Pyb(3))   &
						       - heaviside(binbot(3) - Pyb(3)))
				onfacezb =0.5d0*(sign(1.d0,binbot(3) - ri2(3))   &
						       - sign(1.d0,binbot(3) - ri1(3)))* &
								(heaviside(bintop(1) - Pzb(1))   &
						       - heaviside(binbot(1) - Pzb(1)))* &
								(heaviside(bintop(2) - Pzb(2))   &
						       - heaviside(binbot(2) - Pzb(2)))

				onfacext =0.5d0*(sign(1.d0,bintop(1) - ri2(1))   &
						       - sign(1.d0,bintop(1) - ri1(1)))* &
								(heaviside(bintop(2) - Pxt(2))   &
						       - heaviside(binbot(2) - Pxt(2)))* &
								(heaviside(bintop(3) - Pxt(3))   &
						       - heaviside(binbot(3) - Pxt(3)))
				onfaceyt =0.5d0*(sign(1.d0,bintop(2) - ri2(2))   &
						       - sign(1.d0,bintop(2) - ri1(2)))* &
								(heaviside(bintop(1) - Pyt(1))   &
						       - heaviside(binbot(1) - Pyt(1)))* &
								(heaviside(bintop(3) - Pyt(3))   &
						       - heaviside(binbot(3) - Pyt(3)))
				onfacezt =0.5d0*(sign(1.d0,bintop(3) - ri2(3))   &
						       - sign(1.d0,bintop(3) - ri1(3)))* &
								(heaviside(bintop(1) - Pzt(1))   &
							   - heaviside(binbot(1) - Pzt(1)))* &
								(heaviside(bintop(2) - Pzt(2))   &
						       - heaviside(binbot(2) - Pzt(2)))

				jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer

				!Add Mass flux over face
				mass_flux(cbin(1),cbin(2),cbin(3),1) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),1) & 
				      + nint(dble(onfacexb)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),2) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),2) & 
				      + nint(dble(onfaceyb)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),3) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),3) &
				      + nint(dble(onfacezb)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),4) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),4) &
				      - nint(dble(onfacext)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),5) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),5) &
				      - nint(dble(onfaceyt)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),6) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),6) &
				      - nint(dble(onfacezt)*abs(crossface(jxyz)))

				!if (onfacexb .ne. 0) print*, n, i,j,k,ibin1,ibin2,bintop,halfdomain

				!if (cbin(1) .ge. nbins(1)+1 .or. cbin(1) .le. 2) then
				!	print'(4i8,6f10.5)',iter, cbin, momentum_flux(cbin(1),cbin(2),cbin(3),:,1),momentum_flux(cbin(1),cbin(2),cbin(3),:,4)
				!endif

			enddo
			enddo
			enddo

		endif

	enddo

end subroutine cumulative_mass_flux

!===================================================================================
! Control Volume snapshot of the mass in a given bin

subroutine mass_snapshot
	use module_record
	use field_io, only : mass_bin_io
	use CV_objects, only : CVcheck_mass, CV_debug
	implicit none

	integer										:: n
	integer		,dimension(3)					:: ibin
	integer		,dimension(:,:,:)  ,allocatable	:: volume_mass_temp
	double precision,dimension(3)				:: mbinsize

	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	!Allocate temporary array for mass and momentum in volume
	allocate(volume_mass_temp(nbinso(1),nbinso(2),nbinso(3)))

	!Reset Control Volume momentum 
	volume_mass_temp = 0
	do n = 1,np
		!Add up current volume momentum densities
		ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb
		volume_mass_temp(ibin(1),ibin(2),ibin(3)) = volume_mass_temp(ibin(1),ibin(2),ibin(3)) + 1
	enddo

	!Output Control Volume momentum change and fluxes
	call mass_bin_io(volume_mass_temp,'snap')
	!Create copy of previous timestep Control Volume mass and calculate time evolution
	if (CV_debug) then
		call CVcheck_mass%update_dXdt(volume_mass_temp(:,:,:))
	endif

	deallocate(volume_mass_temp)

end subroutine mass_snapshot


!===================================================================================
! Control Volume Momentum continuity
!===================================================================================

subroutine momentum_flux_averaging(ixyz)
	!use field_io, only :  momentum_flux_io,surface_stress_io, & 
	!					  external_force_io,MOP_stress_io
	use module_record
	use control_volume, only : check_CV_conservation
	use CV_objects, only : CV_debug, CVcheck_momentum
	implicit none

	integer				:: ixyz,icell,jcell,kcell
	integer, save		:: sample_count

	if (vflux_outflag .eq. 0) return

	call cumulative_momentum_flux(ixyz)
	sample_count = sample_count + 1
	if (sample_count .eq. Nvflux_ave) then

		select case(ixyz)
		case(1:3)
			!MOP momentum flux and stresses
			call MOP_stress_io(ixyz)
			Pxy_plane = 0.d0
		case(4)
			!CV momentum flux and stress
			call momentum_flux_io
			momentum_flux = 0.d0
			call momentum_snapshot
			if (external_force_flag .ne. 0 .or. ensemble .eq. tag_move) then
				call external_force_io
				F_ext_bin = 0.d0
			endif
		case default 
			call error_abort("Momentum flux and pressure averaging Error")
		end select

		sample_count = 0

	endif

	!Write forces out at time t before snapshot/final fluxes
	!as both use velocity at v(t-dt/2)
	if (sample_count .eq. Nvflux_ave-1) then
		call surface_stress_io
		Pxyface = 0.d0
		!Debug flag to check CV conservation in code
		if (CV_debug) call CVcheck_momentum%check_error(2,nbins(1)+1,2,nbins(2)+1,2,nbins(3)+1,iter,irank)
	endif

end subroutine momentum_flux_averaging

!===================================================================================
! Momentum Flux over a surface of a bin including all intermediate bins

subroutine cumulative_momentum_flux(ixyz)
	use module_record
	use CV_objects, only : CV_debug, CVcheck_momentum2
	implicit none

	integer							:: ixyz,jxyz,i,j,k,n
	integer							:: planeno
	integer							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	double precision				:: crossplane,rplane,shift
	double precision,dimension(3)	:: mbinsize,velvect,crossface
	double precision,dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

	select case(ixyz)
	case(1:3)
		!MOP momentum flux
		!Shift by half difference between value rounded down and actual value
		!to ensure same distance to top and bottom plane from domain edge
		shift = 0.5d0*(domain(ixyz) - planespacing * (nplanes-1))

		!Add molecular velocities to the configuration stresses
		do n=1,np

			!Replace Signum function with this functions which gives a
			!check for plane crossing and the correct sign 
			crossplane = ceiling((r(ixyz,n)+halfdomain(ixyz)-shift)/planespacing) & 
				    	-ceiling((r(ixyz,n)-delta_t*v(ixyz,n)+halfdomain(ixyz)-shift)/planespacing)

			if (crossplane .ne. 0) then

				!Obtain nearest plane number by integer division 
				!and retrieve location of plane from array
				planeno = ceiling((r(ixyz,n)+halfdomain(ixyz)-shift) 	& 
					  /planespacing)-heaviside(dble(crossplane))+1
				if (planeno .lt. 1) planeno = 1
				if (planeno .gt. nplanes) planeno = nplanes
				rplane = planes(planeno)

				!Calculate velocity at time of intersection
				!crosstime = (r(ixyz,n) - rplane)/v(ixyz,n)
				velvect(:) = v(:,n) !- a(:,n) * crosstime

				!if (crosstime/delta_t .gt. 1.d0)&
                !                    call error_abort("error in kinetic MOP")

				!Obtain stress for three components on y plane
				Pxy_plane(:,planeno) = Pxy_plane(:,planeno) + velvect(:)!*crossplane

			endif
		enddo

	case(4)
		!CV momentum flux
		!Determine bin size
		mbinsize(:) = domain(:) / nbins(:)

		do n = 1,np	

			!Get velocity at v(t+dt/2) from v(t-dt/2)
			velvect(:) = v(:,n)
			ri1(:) = r(:,n) 							!Molecule i at time t
			ri2(:) = r(:,n)	- delta_t*velvect			!Molecule i at time t-dt
			ri12   = ri1 - ri2							!Molecule i trajectory between t-dt and t
			where (ri12 .eq. 0.d0) ri12 = 0.000001d0

			!Assign to bins before and after using integer division
			ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb(:)
			ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb(:)

			!Replace Signum function with this functions which gives a
			!check for plane crossing and the correct sign 
			crossface(:) =  ibin1(:) - ibin2(:)

			if (sum(abs(crossface(:))) .ne. 0) then

				do i = ibin1(1),ibin2(1),sign(1,ibin2(1)-ibin1(1))
				do j = ibin1(2),ibin2(2),sign(1,ibin2(2)-ibin1(2))
				do k = ibin1(3),ibin2(3),sign(1,ibin2(3)-ibin1(3))

					cbin(1) = i; cbin(2) = j; cbin(3) = k

					bintop(:) = (cbin(:)-1*nhb(:)  )*mbinsize(:)-halfdomain(:)
					binbot(:) = (cbin(:)-1*nhb(:)-1)*mbinsize(:)-halfdomain(:)

					!Calculate the plane intersect of trajectory with surfaces of the cube
					Pxt=(/ 			bintop(1), 		     & 
							ri1(2)+(ri12(2)/ri12(1))*(bintop(1)-ri1(1)), & 
							ri1(3)+(ri12(3)/ri12(1))*(bintop(1)-ri1(1))  	/)
					Pxb=(/ 			binbot(1), 		     & 
							ri1(2)+(ri12(2)/ri12(1))*(binbot(1)-ri1(1)), & 
							ri1(3)+(ri12(3)/ri12(1))*(binbot(1)-ri1(1))  	/)
					Pyt=(/	ri1(1)+(ri12(1)/ri12(2))*(bintop(2)-ri1(2)), & 
								bintop(2), 		     & 
							ri1(3)+(ri12(3)/ri12(2))*(bintop(2)-ri1(2))  	/)
					Pyb=(/	ri1(1)+(ri12(1)/ri12(2))*(binbot(2)-ri1(2)), &
								binbot(2), 		     & 
							ri1(3)+(ri12(3)/ri12(2))*(binbot(2)-ri1(2))  	/)
					Pzt=(/	ri1(1)+(ri12(1)/ri12(3))*(bintop(3)-ri1(3)), & 
							ri1(2)+(ri12(2)/ri12(3))*(bintop(3)-ri1(3)), &
								bintop(3) 			/)
					Pzb=(/	ri1(1)+(ri12(1)/ri12(3))*(binbot(3)-ri1(3)), &
							ri1(2)+(ri12(2)/ri12(3))*(binbot(3)-ri1(3)), & 
								binbot(3) 			/)

					onfacexb =0.5d0*(sign(1.d0,binbot(1) - ri2(1)) 	 & 
							       - sign(1.d0,binbot(1) - ri1(1)))* &
									(heaviside(bintop(2) - Pxb(2)) 	 &
							       - heaviside(binbot(2) - Pxb(2)))* &
									(heaviside(bintop(3) - Pxb(3)) 	 &
							       - heaviside(binbot(3) - Pxb(3)))
					onfaceyb =0.5d0*(sign(1.d0,binbot(2) - ri2(2))   &
							       - sign(1.d0,binbot(2) - ri1(2)))* &
									(heaviside(bintop(1) - Pyb(1))   &
							       - heaviside(binbot(1) - Pyb(1)))* &
									(heaviside(bintop(3) - Pyb(3))   &
							       - heaviside(binbot(3) - Pyb(3)))
					onfacezb =0.5d0*(sign(1.d0,binbot(3) - ri2(3))   &
							       - sign(1.d0,binbot(3) - ri1(3)))* &
									(heaviside(bintop(1) - Pzb(1))   &
							       - heaviside(binbot(1) - Pzb(1)))* &
									(heaviside(bintop(2) - Pzb(2))   &
							       - heaviside(binbot(2) - Pzb(2)))

					onfacext =0.5d0*(sign(1.d0,bintop(1) - ri2(1))   &
							       - sign(1.d0,bintop(1) - ri1(1)))* &
									(heaviside(bintop(2) - Pxt(2))   &
							       - heaviside(binbot(2) - Pxt(2)))* &
									(heaviside(bintop(3) - Pxt(3))   &
							       - heaviside(binbot(3) - Pxt(3)))
					onfaceyt =0.5d0*(sign(1.d0,bintop(2) - ri2(2))   &
							       - sign(1.d0,bintop(2) - ri1(2)))* &
									(heaviside(bintop(1) - Pyt(1))   &
							       - heaviside(binbot(1) - Pyt(1)))* &
									(heaviside(bintop(3) - Pyt(3))   &
							       - heaviside(binbot(3) - Pyt(3)))
					onfacezt =0.5d0*(sign(1.d0,bintop(3) - ri2(3))   &
							       - sign(1.d0,bintop(3) - ri1(3)))* &
									(heaviside(bintop(1) - Pzt(1))   &
								   - heaviside(binbot(1) - Pzt(1)))* &
									(heaviside(bintop(2) - Pzt(2))   &
							       - heaviside(binbot(2) - Pzt(2)))

					jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer

					!Calculate velocity at time of intersection
					!crosstime = (r(jxyz,n) - rplane)/v(jxyz,n)
					!velvect(:) = v(:,n) !- a(:,n) * crosstime
					!Change in velocity at time of crossing is not needed as velocity assumed constant 
					!for timestep and changes when forces are applied.

					!Add Momentum flux over face
					momentum_flux(cbin(1),cbin(2),cbin(3),:,1) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,1) & 
					      - velvect(:)*dble(onfacexb)*abs(crossface(jxyz))
					momentum_flux(cbin(1),cbin(2),cbin(3),:,2) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,2) & 
					      - velvect(:)*dble(onfaceyb)*abs(crossface(jxyz))
					momentum_flux(cbin(1),cbin(2),cbin(3),:,3) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,3) &
					      - velvect(:)*dble(onfacezb)*abs(crossface(jxyz))
					momentum_flux(cbin(1),cbin(2),cbin(3),:,4) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,4) &
					      + velvect(:)*dble(onfacext)*abs(crossface(jxyz))
					momentum_flux(cbin(1),cbin(2),cbin(3),:,5) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,5) &
					      + velvect(:)*dble(onfaceyt)*abs(crossface(jxyz))
					momentum_flux(cbin(1),cbin(2),cbin(3),:,6) = & 
						momentum_flux(cbin(1),cbin(2),cbin(3),:,6) &
					      + velvect(:)*dble(onfacezt)*abs(crossface(jxyz))

					!Add instantanous Momentum flux to CV record
					if (CV_debug) then
    					CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,1) = & 
    						CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,1) &
							   - velvect(:)*dble(onfacexb)*abs(crossface(jxyz))
    					CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,2) = & 
    						CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,2)  & 
    					      - velvect(:)*dble(onfaceyb)*abs(crossface(jxyz))
    					CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,3) = & 
    						CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,3)  & 
    					      - velvect(:)*dble(onfacezb)*abs(crossface(jxyz))
    					CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,4) = & 
    						CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,4)  & 
    					      + velvect(:)*dble(onfacext)*abs(crossface(jxyz))
    					CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,5) = & 
    						CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,5)  & 
    					      + velvect(:)*dble(onfaceyt)*abs(crossface(jxyz))
    					CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,6) = & 
    						CVcheck_momentum2%flux(cbin(1),cbin(2),cbin(3),:,6)  & 
    					      + velvect(:)*dble(onfacezt)*abs(crossface(jxyz))			
					endif

				enddo
				enddo
				enddo

			endif

		enddo
	case default 
		call error_abort("Cumulative Momentum flux Error")
	end select

end subroutine cumulative_momentum_flux

!===================================================================================
! Control Volume snapshot of momentum in a given bin

subroutine momentum_snapshot
	use field_io, only : velocity_bin_io
	use module_record
	implicit none

	integer											:: n
	integer		,dimension(3)						:: ibin
	integer		,dimension(:,:,:)  ,allocatable		:: volume_mass_temp
	double precision								:: binvolume
	double precision,dimension(3)					:: mbinsize
	double precision,dimension(:,:,:,:),allocatable :: volume_momentum_temp

	mbinsize(:) = domain(:) / nbins(:)

	!Allocate temporary array for mass and momentum in volume
	allocate(volume_mass_temp(nbinso(1),nbinso(2),nbinso(3)))
	allocate(volume_momentum_temp(nbinso(1),nbinso(2),nbinso(3),3  ))

	!Reset Control Volume momentum 
	volume_mass_temp = 0
	volume_momentum_temp = 0.d0
	do n = 1,np
		!Add up current volume momentum densities
		ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb
		volume_mass_temp(ibin(1),ibin(2),ibin(3)) = volume_mass_temp(ibin(1),ibin(2),ibin(3)) + 1
		volume_momentum_temp(ibin(1),ibin(2),ibin(3),:) = volume_momentum_temp(ibin(1),ibin(2),ibin(3),:) + v(:,n)
	enddo
	binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
	volume_momentum_temp = volume_momentum_temp/binvolume

	!Output Control Volume momentum change and fluxes
	call velocity_bin_io(volume_mass_temp,volume_momentum_temp,'snap')

	deallocate(volume_mass_temp)
	deallocate(volume_momentum_temp)

end subroutine momentum_snapshot


!===================================================================================
! Control Volume Energy continuity
!===================================================================================

subroutine energy_flux_averaging(ixyz)
	!use field_io, only : energy_flux_io,surface_power_io,MOP_energy_io
	use module_record
	implicit none

	integer				:: ixyz
	integer, save		:: sample_count

	if (eflux_outflag .eq. 0) return

	call cumulative_energy_flux(ixyz)
	sample_count = sample_count + 1
	if (sample_count .eq. Neflux_ave) then

		select case(ixyz)
		case(1:3)
			!MOP energy flux and Power (stresses*velocity)
			call MOP_energy_io(ixyz)
			Pxy_plane = 0.d0
		case(4)
			!CV energy flux and Power (stresses*velocity)
			call energy_flux_io
			call surface_power_io
			energy_flux = 0.d0
			Pxyvface = 0.d0
			call energy_snapshot
		case default 
			call error_abort("Energy flux averaging Error")
		end select

		sample_count = 0

	endif

end subroutine energy_flux_averaging

!===================================================================================
! Energy Flux over a surface of a bin including all intermediate bins

subroutine cumulative_energy_flux(ixyz)
	use module_record
	implicit none

	integer							:: ixyz,jxyz,i,j,k,n
	integer							:: planeno,onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	double precision				:: crosstime,crossplane,rplane,shift,energy
	double precision,dimension(3)	:: mbinsize,velvect,crossface
	double precision,dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

	select case(ixyz)
	case(1:3)
		!MOP energy flux
		!Shift by half difference between value rounded down and actual value
		!to ensure same distance to top and bottom plane from domain edge
		shift = 0.5d0*(domain(ixyz) - planespacing * (nplanes-1))

		!Add molecular velocities to the configuration stresses
		do n=1,np

			!Replace Signum function with this functions which gives a
			!check for plane crossing and the correct sign 
			crossplane = ceiling((r(ixyz,n)+halfdomain(ixyz)-shift)/planespacing) & 
				    	-ceiling((r(ixyz,n)-delta_t*v(ixyz,n)+halfdomain(ixyz)-shift)/planespacing)

			if (crossplane .ne. 0) then

				!Obtain nearest plane number by integer division 
				!and retrieve location of plane from array
				planeno = ceiling((r(ixyz,n)+halfdomain(ixyz)-shift) 	& 
					  /planespacing)-heaviside(dble(crossplane))+1
				if (planeno .lt. 1) planeno = 1
				if (planeno .gt. nplanes) planeno = nplanes
				rplane = planes(planeno)

				!Calculate energy at intersection
				velvect(:) = v(:,n) + 0.5d0*a(:,n)*delta_t
				energy = 0.5d0 * (dot_product(velvect,velvect) + potenergymol(n))

				if (crosstime/delta_t .gt. 1.d0) call error_abort("error in kinetic MOP")

				!Obtain stress for three components on y plane
				Pxyv_plane(planeno) = Pxyv_plane(planeno) + energy!*crossplane

			endif
		enddo

	case(4)
		!CV momentum flux
		!Determine bin size
		mbinsize(:) = domain(:) / nbins(:)

		do n = 1,np

			ri1(:) = r(:,n)					!Molecule i at time  t
			ri2(:) = r(:,n)-delta_t*v(:,n)	!Molecule i at time t-dt
			ri12   = ri1 - ri2				!Molecule i trajectory between t-dt and t
			where (ri12 .eq. 0.d0) ri12 = 0.000001d0

			!Assign to bins before and after using integer division
			ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb(:)
			ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb(:)

			!Replace Signum function with this functions which gives a
			!check for plane crossing and the correct sign 
			crossface(:) =  ibin1(:) - ibin2(:)

			if (sum(abs(crossface(:))) .ne. 0) then

				do i = ibin1(1),ibin2(1),sign(1,ibin2(1)-ibin1(1))
				do j = ibin1(2),ibin2(2),sign(1,ibin2(2)-ibin1(2))
				do k = ibin1(3),ibin2(3),sign(1,ibin2(3)-ibin1(3))

					cbin(1) = i; cbin(2) = j; cbin(3) = k

					bintop(:) = (cbin(:)-1*nhb(:)  )*mbinsize(:)-halfdomain(:)
					binbot(:) = (cbin(:)-1*nhb(:)-1)*mbinsize(:)-halfdomain(:)

					!Calculate the plane intersect of trajectory with surfaces of the cube
					Pxt=(/ 			bintop(1), 		     & 
							ri1(2)+(ri12(2)/ri12(1))*(bintop(1)-ri1(1)), & 
							ri1(3)+(ri12(3)/ri12(1))*(bintop(1)-ri1(1))  	/)
					Pxb=(/ 			binbot(1), 		     & 
							ri1(2)+(ri12(2)/ri12(1))*(binbot(1)-ri1(1)), & 
							ri1(3)+(ri12(3)/ri12(1))*(binbot(1)-ri1(1))  	/)
					Pyt=(/	ri1(1)+(ri12(1)/ri12(2))*(bintop(2)-ri1(2)), & 
									bintop(2), 		     & 
							ri1(3)+(ri12(3)/ri12(2))*(bintop(2)-ri1(2))  	/)
					Pyb=(/	ri1(1)+(ri12(1)/ri12(2))*(binbot(2)-ri1(2)), &
									binbot(2), 		     & 
							ri1(3)+(ri12(3)/ri12(2))*(binbot(2)-ri1(2))  	/)
					Pzt=(/	ri1(1)+(ri12(1)/ri12(3))*(bintop(3)-ri1(3)), & 
							ri1(2)+(ri12(2)/ri12(3))*(bintop(3)-ri1(3)), &
									bintop(3) 			/)
					Pzb=(/	ri1(1)+(ri12(1)/ri12(3))*(binbot(3)-ri1(3)), &
							ri1(2)+(ri12(2)/ri12(3))*(binbot(3)-ri1(3)), & 
									binbot(3) 			/)

					onfacexb =0.5d0*(sign(1.d0,binbot(1) - ri2(1)) 	 & 
						      	   - sign(1.d0,binbot(1) - ri1(1)))* &
								(heaviside(bintop(2) - Pxb(2)) 	 &
						       - heaviside(binbot(2) - Pxb(2)))* &
								(heaviside(bintop(3) - Pxb(3)) 	 &
						       - heaviside(binbot(3) - Pxb(3)))
					onfaceyb =0.5d0*(sign(1.d0,binbot(2) - ri2(2))   &
						       	   - sign(1.d0,binbot(2) - ri1(2)))* &
								(heaviside(bintop(1) - Pyb(1))   &
						       - heaviside(binbot(1) - Pyb(1)))* &
								(heaviside(bintop(3) - Pyb(3))   &
						       - heaviside(binbot(3) - Pyb(3)))
					onfacezb =0.5d0*(sign(1.d0,binbot(3) - ri2(3))   &
						       	   - sign(1.d0,binbot(3) - ri1(3)))* &
								(heaviside(bintop(1) - Pzb(1))   &
						       - heaviside(binbot(1) - Pzb(1)))* &
								(heaviside(bintop(2) - Pzb(2))   &
						       - heaviside(binbot(2) - Pzb(2)))

					onfacext =0.5d0*(sign(1.d0,bintop(1) - ri2(1))   &
						       	   - sign(1.d0,bintop(1) - ri1(1)))* &
								(heaviside(bintop(2) - Pxt(2))   &
						       - heaviside(binbot(2) - Pxt(2)))* &
				            	(heaviside(bintop(3) - Pxt(3))   &
						       - heaviside(binbot(3) - Pxt(3)))
					onfaceyt =0.5d0*(sign(1.d0,bintop(2) - ri2(2))   &
						       	   - sign(1.d0,bintop(2) - ri1(2)))* &
								(heaviside(bintop(1) - Pyt(1))   &
						       - heaviside(binbot(1) - Pyt(1)))* &
								(heaviside(bintop(3) - Pyt(3))   &
						       - heaviside(binbot(3) - Pyt(3)))
					onfacezt =0.5d0*(sign(1.d0,bintop(3) - ri2(3))   &
						       	   - sign(1.d0,bintop(3) - ri1(3)))* &
								(heaviside(bintop(1) - Pzt(1))   &
    						   - heaviside(binbot(1) - Pzt(1)))* &
								(heaviside(bintop(2) - Pzt(2))   &
						       - heaviside(binbot(2) - Pzt(2)))

					jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer

					!Calculate velocity at time of intersection
					!crosstime = (r(jxyz,n) - rplane)/v(jxyz,n)
					velvect(:) = v(:,n) + 0.5d0*a(:,n)*delta_t
					energy = 0.5d0 * (dot_product(velvect,velvect) + potenergymol(n))
					!Change in velocity at time of crossing is not needed as velocity assumed constant 
					!for timestep and changes when forces are applied.

					!Add Energy flux over face
					energy_flux(cbin(1),cbin(2),cbin(3),1) = & 
						energy_flux(cbin(1),cbin(2),cbin(3),1) & 
					      - energy*dble(onfacexb)*abs(crossface(jxyz))
					energy_flux(cbin(1),cbin(2),cbin(3),2) = & 
						energy_flux(cbin(1),cbin(2),cbin(3),2) & 
					      - energy*dble(onfaceyb)*abs(crossface(jxyz))
					energy_flux(cbin(1),cbin(2),cbin(3),3) = & 
						energy_flux(cbin(1),cbin(2),cbin(3),3) &
					      - energy*dble(onfacezb)*abs(crossface(jxyz))
					energy_flux(cbin(1),cbin(2),cbin(3),4) = & 
						energy_flux(cbin(1),cbin(2),cbin(3),4) &
					      + energy*dble(onfacext)*abs(crossface(jxyz))
					energy_flux(cbin(1),cbin(2),cbin(3),5) = & 
						energy_flux(cbin(1),cbin(2),cbin(3),5) &
					      + energy*dble(onfaceyt)*abs(crossface(jxyz))
					energy_flux(cbin(1),cbin(2),cbin(3),6) = & 
						energy_flux(cbin(1),cbin(2),cbin(3),6) &
					      + energy*dble(onfacezt)*abs(crossface(jxyz))

				enddo
				enddo
				enddo

			endif

		enddo
	case default 
		call error_abort("Cumulative Energy flux Error")
	end select

end subroutine cumulative_energy_flux

!===================================================================================
! Control Volume snapshot of momentum in a given bin

subroutine energy_snapshot
	use field_io, only : energy_bin_io
	use librarymod
	use module_record
	implicit none

	integer											:: n
	integer		,dimension(3)						:: ibin
	double precision								:: binvolume, energy
	double precision,dimension(3)					:: mbinsize,velvect
	double precision,dimension(:,:,:),allocatable 	:: volume_energy_temp

	mbinsize(:) = domain(:) / nbins(:)

	!Allocate temporary array for energy in volume
	allocate(volume_energy_temp(nbinso(1),nbinso(2),nbinso(3)))

	!Reset Control Volume momentum 
	volume_energy_temp = 0.d0
	do n = 1,np
		!Add up current volume momentum densities
		ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb(:)
		velvect(:) = v(:,n) + 0.5d0*a(:,n)*delta_t
		energy = 0.5d0 * (dot_product(velvect,velvect) + potenergymol(n))
		volume_energy_temp(ibin(1),ibin(2),ibin(3)) = volume_energy_temp(ibin(1),ibin(2),ibin(3)) + energy
	enddo

	binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
	volume_energy_temp = volume_energy_temp/(binvolume)

	!Output Control Volume momentum change and fluxes
	call energy_bin_io(volume_energy_temp,'snap')

	deallocate(volume_energy_temp)

end subroutine energy_snapshot

!====================================================================================
!
!	Force based statistics called from computation of forces routine
!
!====================================================================================
! Called from force routine - adds up all molecular interactions

subroutine pressure_tensor_forces(molno, rij, accijmag)
	use module_record
	implicit none

	integer										:: ixyz, jxyz
	integer,intent(in)							:: molno
	double precision,intent(in)     			:: accijmag    !Non directional component of acceleration
	double precision,dimension(3),intent(in)   	:: rij         !vector between particles i and j

	do ixyz = 1,3
	do jxyz = 1,3
		rfmol(molno,ixyz,jxyz) = rfmol(molno,ixyz,jxyz) + accijmag*rij(ixyz)*rij(jxyz)
	enddo
	enddo

end subroutine pressure_tensor_forces

!====================================================================================
! CONFIGURATIONAL INTERACTIONS ASSIGNED TO CONTAINING BIN - HALF PER MOLECULE 
! Use a configurational expression, similar to the virial with interaction partitioned 
! between bins with HALF given per bin and nothing in intermediate bins
! This is also called the Harasima contour in NAMD (MD package) and the paper it references
! Jacob Sonne,a Flemming Y. Hansen, and Günther H. Peters J.CHEM.PHYS. 122, 124903 (2005)
! The original paper by A. Harasima, Advances in Chemical Physics, Volume 1 P201 is unavailable

subroutine pressure_tensor_forces_H(ri,rj,rij,accijmag)
	use module_record
	implicit none

	integer											:: ixyz, jxyz
	integer											:: diff
	integer,dimension(3)							:: ibin, jbin, bindiff 
	double precision,intent(in)            			:: accijmag    !Non directional component of acceleration
	double precision,dimension(3), intent(in)		:: rij, ri, rj
	double precision,dimension(3)					:: VAbinsize

	!================================================================
	!= Establish bins for molecules & check number of required bins	=
	!================================================================

	do ixyz = 1,nd

		VAbinsize(ixyz) = domain(ixyz) / nbins(ixyz)
		if (VAbinsize(ixyz) .lt. cellsidelength(ixyz)) stop "Binsize bigger than cellsize ~ Not ok for volume averaging"

		!Determine current bins using integer division
		ibin(ixyz) = ceiling((ri(ixyz)+halfdomain(ixyz))/VAbinsize(ixyz)) + 1 !Establish current i bin
		jbin(ixyz) = ceiling((rj(ixyz)+halfdomain(ixyz))/VAbinsize(ixyz)) + 1 !Establish current j bin

		!Check number of bins between molecules
		bindiff(ixyz) = abs(ibin(ixyz) - jbin(ixyz)) + 1

	enddo

	!================================================================
	!=			Assign to bins				=
	!================================================================

	!Check difference between bins i and j
	diff = bindiff(1)+bindiff(2)+bindiff(3)

	!Ignore values outside of domain resulting from shifted bin 
	if (minval(ibin) .lt. 1) diff = 0
	if (minval(jbin) .lt. 1) diff = 0
	if (maxval(ibin) .gt. maxval(nbins)+2) diff = 0
	if (maxval(jbin) .gt. maxval(nbins)+2) diff = 0

	select case(diff)
	!================Skip Force addition===========================
	case(0)
	!Do Nothing

	!================Molecules in same bin===========================
	case(3)

		!------Add molecules to bin-----
		do ixyz = 1,nd
		do jxyz = 1,nd
			!Factor of two as molecule i and molecule j are both counted in bin i
			rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) = rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) &
				  + accijmag*rij(ixyz)*rij(jxyz)
		enddo
		enddo

	!===================Interactions split over 2 cells only==========
	case(4:6)

		!------Add molecules to bin-----
		!Molecule i and j contribution split between bins
		do ixyz = 1,nd
		do jxyz = 1,nd
			!-----------Molecule i bin-----------
			rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) = rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) &
				  + 0.5d0*accijmag*rij(ixyz)*rij(jxyz)

			!-----------Molecule j bin-----------
			rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) = rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) &
				  + 0.5d0*accijmag*rij(ixyz)*rij(jxyz)

		enddo
		enddo
			
	case default

		stop "Harasima Stress AVERAGING ERROR"

	end select

end subroutine pressure_tensor_forces_H

!====================================================================================
! VOLUME AVERAGE CONFIGURATIONAL EXPRESSION
! Author: David Trevelyan, July 2013
! Linear trajectory path sampled to find approximate values of l_ij (less accurate, 
! but much easier to understand, and no problems with polar coordinates)
subroutine pressure_tensor_forces_VA_trap(ri,rj,accijmag)
	use computational_constants_MD, only: domain, halfdomain, VA_line_samples
	use calculated_properties_MD, only: nbins, rfbin
	use librarymod, only: outerprod
	
	implicit none

	real(kind(0.d0)), intent(in) :: accijmag
	real(kind(0.d0)), dimension(3), intent(in) :: ri, rj

	integer :: ss
	integer :: Ns
	real(kind(0.d0)) :: s, ds, rs(3), rij(3), Fij(3), VAbinsize(3), bin(3), rF(3,3)

	VAbinsize(:) = domain(:) / nbins(:)
	rij = rj - ri
	rF = outerprod(rij, accijmag*rij)

	! Split line l_ij into segments of size ds
	Ns = VA_line_samples
	ds = 1.d0 / real(Ns, kind(0.d0))
	! First sample at midpoint of first segment 
	s = 0.5d0*ds 

	! Loop over all samples, s varies from 0 to 1
	do ss = 1, Ns

		! Position of sample on line
		rs(:) = ri(:) + s*rij(:)	

		! Don't count if sample is outside the domain (will be picked up
		! by neighbouring processors)
		if ( .not. any( abs(rs(:)) .gt. halfdomain(:) ) ) then

			bin(:) = ceiling((rs(:)+halfdomain(:))/VAbinsize(:)) + 1
			rfbin(bin(1),bin(2),bin(3),:,:) =  &
			rfbin(bin(1),bin(2),bin(3),:,:) + rF(:,:)/real(Ns,kind(0.d0))

		end if

		s = s + ds	

	end do	
	
end subroutine pressure_tensor_forces_VA_trap

subroutine pressure_tensor_forces_VA_trap_cpol(ri,rj,accijmag)
	use concentric_cylinders
	use computational_constants_MD, only: domain, halfdomain, VA_line_samples
	use physical_constants_MD, only: pi
	use calculated_properties_MD, only: rfbin
	use librarymod, only: outerprod, cartesianiser, cpolariser, cpolariseT
	use messenger, only: localise, globalise
	implicit none

	real(kind(0.d0)), intent(in) :: accijmag
	real(kind(0.d0)), intent(in) :: ri(3), rj(3)

	integer :: ss
	integer :: br, bt, bz
	integer :: Ns
	real(kind(0.d0)) :: s, ds, rs(3), rij(3), Fij(3), rs_cart(3), VAbinsize(3), rF(3,3)
	real(kind(0.d0)) :: ripol(3), rjpol(3), rijpol(3)

	! Calculate relevant polar quantities
	ripol = cpolariser(globalise(ri)) 
	rjpol = cpolariser(globalise(rj))
	rijpol = rjpol - ripol 
	! Periodic in theta so take minimum image
	rijpol(2) = rijpol(2) - nint(rijpol(2)/(2.d0*pi))*2.d0*pi

	! Store rij * Fij outer product tensor (cartesian)
	rij = rj - ri
	rF = outerprod(rij, accijmag*rij)
	! Transform to polar coordinates
	rF = cpolariseT(rF,ripol(2)) 

	! Bin sizes
	VAbinsize(1) = (r_io - r_oi) / cpol_bins(1)
	VAbinsize(2) = 2.d0*pi       / cpol_bins(2)
	VAbinsize(3) = domain(3)     / cpol_bins(3)

	! First sample at midpoint of first segment 
	Ns = VA_line_samples
	ds = 1.d0 / real(Ns, kind(0.d0))
	s = 0.5d0*ds 

	! Loop over all samples, s varies from 0 to 1
	do ss = 1, Ns

		! Position of sample on line
		rs(:) = ripol(:) + s*rijpol(:)	
		rs_cart = localise(cartesianiser(rs))

		! Don't count if sample is outside the domain (will be picked up
		! by neighbouring processors)
		if (  .not. any(abs(rs_cart(:)).gt.halfdomain(:))  ) then

			! Binning conventions 
 			rs(1)  = rs(1) - r_oi
			rs(2)  = modulo(rs(2),2.d0*pi)
			rs(3)  = rs_cart(3) + halfdomain(3) 

			!Add to cylindrical bins
			br = ceiling(rs(1)/VAbinsize(1)) 
			bt = ceiling(rs(2)/VAbinsize(2)) 
			bz = ceiling(rs(3)/VAbinsize(3)) + cpol_nhbz

			!Ignore molecules in cylinder region
			if ( br .ge. 1 .and. br .le. cpol_bins(1) ) then
				rfbin(br,bt,bz,:,:) =  &
				rfbin(br,bt,bz,:,:) + rF(:,:)/real(Ns,kind(0.d0))
			end if

		end if

		s = s + ds	

	end do	
	
end subroutine pressure_tensor_forces_VA_trap_cpol


!===============================================================================
! VOLUME AVERAGE CONFIGURATIONAL EXPRESSION
! Use a configurational expression, similar to the virial with interaction partitioned 
! between bins based on the share of the interaction between two molecules in a given bin
! N.B. Assume no more than 1 bin between bins containing the molecules
! This is the Volume Average stress of Lutsko, although this is also called the
! Irving Kirkwood contour in NAMD (MD package) and the paper it references
! Jacob Sonne,a Flemming Y. Hansen, and Günther H. Peters J.CHEM.PHYS. 122, 124903 (2005)
! 								╭∩╮（︶︿︶）╭∩╮﻿
subroutine pressure_tensor_forces_VA(ri,rj,rij,accijmag)
	use module_record
	implicit none

	integer											:: ixyz, jxyz,i,j,k,l,n
	integer											:: diff
	integer,dimension(3)							:: ibin, jbin, bindiff 
	integer,dimension(:,:)		   ,allocatable		:: interbin, interbindiff
	double precision,intent(in)            			:: accijmag    !Non directional component of acceleration
	double precision,dimension(3), intent(in)		:: rij, ri, rj
	double precision,dimension(3)					:: VAbinsize, normal, p,temp1,temp2,temp3
	double precision,dimension(:,:)	   ,allocatable	:: intersection
	double precision,dimension(:,:,:)  ,allocatable	:: MLfrac !Magnitude of fraction of stress
	double precision,dimension(:,:,:,:),allocatable	:: Lfrac  !Fraction of stress in a given cell


	!================================================================
	!= Establish bins for molecules & check number of required bins	=
	!================================================================

	do ixyz = 1,nd

		VAbinsize(ixyz) = domain(ixyz) / nbins(ixyz)
		if (VAbinsize(ixyz) .lt. cellsidelength(ixyz)) &
                    call error_abort("Binsize bigger than cellsize ~ Not ok for volume averaging")

		!Determine current bins using integer division
		ibin(ixyz) = ceiling((ri(ixyz)+halfdomain(ixyz))/VAbinsize(ixyz)) + 1 !Establish current i bin
		jbin(ixyz) = ceiling((rj(ixyz)+halfdomain(ixyz))/VAbinsize(ixyz)) + 1 !Establish current j bin

		!Check number of bins between molecules
		bindiff(ixyz) = abs(ibin(ixyz) - jbin(ixyz)) + 1

	enddo

	!================================================================
	!=			Assign to bins				=
	!================================================================

	!Check difference between bins i and j
	diff = bindiff(1)+bindiff(2)+bindiff(3)


	!Ignore values outside of domain resulting from shifted bin 
	if (minval(ibin) .lt. 1) diff = 0
	if (minval(jbin) .lt. 1) diff = 0
	if (maxval(ibin) .gt. maxval(nbins)+2) diff = 0
	if (maxval(jbin) .gt. maxval(nbins)+2) diff = 0

	select case(diff)
	!================Skip Force addition===========================
	case(0)
	!Do Nothing
	!print*, maxval(ri(:)+halfdomain(:)), 'is outside of domain', domain(1)
	!================Molecules in same bin===========================
	case(3)

		!------Add molecules to bin-----
		do ixyz = 1,nd
		do jxyz = 1,nd
			!Factor of two as molecule i and molecule j are both counted in bin i
			rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) = rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) &
				  + accijmag*rij(ixyz)*rij(jxyz)
		enddo
		enddo

	!===================Interactions split over 2 cells only==========
	case(4)

		!Allocate array for vectors from molecules to bin walls
		allocate(Lfrac(bindiff(1),bindiff(2),bindiff(3),nd))
		!Allocate arrays for magnitude of vector Lfrac
		allocate(MLfrac(bindiff(1),bindiff(2),bindiff(3)))
		!Allocate array for one bin intersection point
		allocate(intersection(3,1))
		Lfrac = 0.d0
		MLfrac = 0.d0

		do ixyz = 1,nd
			!Test to see over which coordinate molecules are in different bins
			if (bindiff(ixyz) .ne. 1) then

				!Set normal to plane for ixyz
				normal = 0.d0
				normal(ixyz) = 1.d0

				!Establish location of plane between ri and rj
				if (ri(ixyz) .lt. rj(ixyz)) then
					p(ixyz) = (ibin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				elseif (ri(ixyz) .gt. rj(ixyz)) then
					p(ixyz) = (jbin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				endif

				!Calculate location of intersection of rij and plane
				call plane_line_intersect(intersection(:,1),normal,p,ri,rj)

				!Calculate vectors from ri & rj to intersect
				Lfrac(1,1,1,:) = ri(:)-intersection(:,1)
				Lfrac(bindiff(1),bindiff(2),bindiff(3),:) = rj(:)-intersection(:,1)
				
			endif
		enddo

		!Calculate magnitude of 3 vector components
		do ixyz = 1,3
			MLfrac(:,:,:) = MLfrac(:,:,:) + Lfrac(:,:,:,ixyz)**2
		enddo
		MLfrac(:,:,:) = MLfrac(:,:,:)**0.5d0
		!Normalise to one
		MLfrac(:,:,:) = MLfrac(:,:,:)/magnitude(rij(:))

		!------Add molecules to bin-----
		!Molecule i and j contribution split between bins
		do ixyz = 1,nd
		do jxyz = 1,nd
			!-----------Molecule i bin-----------
			rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) = rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) &
				  + accijmag*rij(ixyz)*rij(jxyz)*MLfrac(1,1,1)

			!-----------Molecule j bin-----------
			rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) = rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) &
				  + accijmag*rij(ixyz)*rij(jxyz)*MLfrac(bindiff(1),bindiff(2),bindiff(3))

		enddo
		enddo

		deallocate(intersection)
		deallocate(Lfrac)
		deallocate(MLfrac)
		
	!==============Interactions split over intermediate cell===========
	case(5)

		!Allocate array for vectors from molecules to bin walls
		allocate(Lfrac(bindiff(1),bindiff(2),bindiff(3),nd))
		!Allocate arrays for magnitude of vector Lfrac
		allocate(MLfrac(bindiff(1),bindiff(2),bindiff(3)))
		!Allocate array for bin intersection points
		allocate(intersection(3,2))
		!Allocate array for intersection points
		allocate(interbin(3,1))
		allocate(interbindiff(3,1))

		!Set intersection location array to zero
		n = 1
		intersection = 0.d0
		Lfrac = 0.d0
		MLfrac = 0.d0

		do ixyz = 1,nd

			!Test to see over which coordinate molecules are in different bins
			if (bindiff(ixyz) .ne. 1) then

				!Set normal to plane for ixyz
				normal = 0.d0
				normal(ixyz) = 1.d0

				!Establish location of plane between ri and rj
				if (ri(ixyz) .lt. rj(ixyz)) then
					p(ixyz) = (ibin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				elseif (ri(ixyz) .ge. rj(ixyz)) then
					p(ixyz) = (jbin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				endif

				!Calculate location of intersection of rij and plane
				call plane_line_intersect(intersection(:,n),normal,p,ri,rj)

				n = n + 1

			endif
		enddo

		!Take average of 2 cell side intersections to determine intermediate cell
		do ixyz=1,3
			interbin(ixyz,1) = ceiling((0.5d0*(intersection(ixyz,1) &
					+intersection(ixyz,2))+halfdomain(ixyz))/VAbinsize(ixyz))+1
			interbindiff(ixyz,1) = abs(interbin(ixyz,1)-ibin(ixyz)) + 1
		enddo



		!Check which plane is closest to i and which corresponds to j
		temp1 = ri(:)-intersection(:,1)
		temp2 = ri(:)-intersection(:,2)
		!if (magnitude(ri(:)-intersection(:,1)) .le.  & 
  	    !   	    magnitude(ri(:)-intersection(:,2))) then
		if (magnitude(temp1) .le.magnitude(temp1)) then
			i = 1
			j = 2
		else
			i = 2
			j = 1
		endif

		!Fix for vectors going directly through vertex of bin
		if (all(interbindiff(:,1) .eq. 1)) then
			!Ensure not in same bin as 1
			if (bindiff(1)+bindiff(2) .eq. 3) then
				interbindiff(3,1) = 2 
			else 
				interbindiff(1,1) = 2
			endif
		endif
		if (all(interbindiff(:,1) .eq. bindiff(:))) then
			!Ensure not in same bin as bindiff
			if (bindiff(1)+bindiff(2) .eq. 3) then
				interbindiff(3,1) = 1 
			else 
				interbindiff(1,1) = 1
			endif
		endif

		!Calculate vectors from ri to intersect and rj to intersect
		Lfrac(1,1,1,:) = ri(:)-intersection(:,i)
		Lfrac(bindiff(1),bindiff(2),bindiff(3),:) = rj(:)-intersection(:,j)
		Lfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1),:) = intersection(:,i)-intersection(:,j)

		!Calculate magnitude of 3 vector components
		do ixyz = 1,3
			MLfrac(:,:,:) = MLfrac(:,:,:) + Lfrac(:,:,:,ixyz)**2
		enddo
		MLfrac(:,:,:) = MLfrac(:,:,:)**0.5d0

		!Normalise to one
		MLfrac(:,:,:) = MLfrac(:,:,:)/magnitude(rij(:))

		!------Add stress component to bin weighted by line length-----
		!Intermediate bin is either in domain or in halo. 
		!For halo bins the stress is added to the cell the halo represents. 
		!For domain cells it is added directly to that cell
		do ixyz = 1,nd
		do jxyz = 1,nd

			!-----------Molecule i bin-----------
			rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) = rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(1,1,1)

			!-----------Intermediate Bin-----------
			!If both bins are in the domain then contribution is added for both molecules
			rfbin(interbin(1,1),interbin(2,1),interbin(3,1),ixyz,jxyz) = &
				rfbin(interbin(1,1),interbin(2,1),interbin(3,1),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1))

			!-----------Molecule j bin-----------
			rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) = rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(bindiff(1),bindiff(2),bindiff(3))

		enddo
		enddo

		deallocate(Lfrac)
		deallocate(MLfrac)
		deallocate(interbin)
		deallocate(intersection)
		deallocate(interbindiff)

	!===========Interactions split over 2 intermediate cell===============
	case (6)

		!Allocate array for vectors from molecules to bin walls
		allocate(Lfrac(bindiff(1),bindiff(2),bindiff(3),nd))
		!Allocate arrays for magnitude of vector Lfrac
		allocate(MLfrac(bindiff(1),bindiff(2),bindiff(3)))
		!Allocate array for bin intersection points
		allocate(intersection(3,3))
		!Allocate array for intersection points
		allocate(interbin(3,2))
		allocate(interbindiff(3,2))

		!Set intersection location array to zero
		n = 1
		intersection = 0.d0

		Lfrac = 0.d0
		MLfrac = 0.d0

		do ixyz = 1,nd

			!Test to see over which coordinate molecules are in different bins
			if (bindiff(ixyz) .ne. 1) then

				!Set normal to plane for ixyz
				normal = 0.d0
				normal(ixyz) = 1.d0

				!Establish location of plane between ri and rj
				if (ri(ixyz) .lt. rj(ixyz)) then
					p(ixyz) = (ibin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				elseif (ri(ixyz) .gt. rj(ixyz)) then
					p(ixyz) = (jbin(ixyz)-1)*VAbinsize(ixyz)-halfdomain(ixyz)
				endif

				!Calculate location of intersection of rij and plane
				call plane_line_intersect(intersection(:,n),normal,p,ri,rj)

				n = n+1

			endif
		enddo

		!Determine which intersection covers both intermediate cells
		! |(1)-(2)| > |(1)-(3)| > |(2)-(3)| then 1,2 must cover
		! both cells while 1,3 and 2,3 are the cell intercepts
		temp1 = intersection(:,1)-intersection(:,2)
		temp2 = intersection(:,3)-intersection(:,2)
		temp3 = intersection(:,1)-intersection(:,3)
		!if (magnitude(intersection(:,1)-intersection(:,2)) .gt. & 
		!    magnitude(intersection(:,3)-intersection(:,2))) then
		!	if (magnitude(intersection(:,1)-intersection(:,2)) .gt. & 
		!	    magnitude(intersection(:,1)-intersection(:,3))) then
		if (magnitude(temp1).gt.magnitude(temp2)) then
			if (magnitude(temp1).gt.magnitude(temp3)) then
				k = 1
				l = 2
			else
				k = 1
				l = 3
			endif
		else
			!if (magnitude(intersection(:,3)-intersection(:,2)) .gt. & 
			!    magnitude(intersection(:,1)-intersection(:,3))) then
			if (magnitude(temp2).gt.magnitude(temp3)) then
				k = 2
				l = 3
			else
				k = 1
				l = 3
			endif
		endif

		!Take average of cell side intersections to determine intermediate bins
		!k and l are used to determine intermediate bins
		do ixyz=1,3
		
			interbin(ixyz,1) = ceiling((0.5d0*(intersection(ixyz,(6-k-l)) &
					    +intersection(ixyz,l))+halfdomain(ixyz))/VAbinsize(ixyz))+1
			interbindiff(ixyz,1) = abs(interbin(ixyz,1)-ibin(ixyz)) + 1

			interbin(ixyz,2) = ceiling((0.5d0*(intersection(ixyz,(6-k-l)) &
					    +intersection(ixyz,k))+halfdomain(ixyz))/VAbinsize(ixyz))+1
			interbindiff(ixyz,2) = abs(interbin(ixyz,2)-ibin(ixyz)) + 1

		enddo

		!Check which plane is closest to i
		temp1 = ri(:)-intersection(:,1)
		temp2 = ri(:)-intersection(:,2)
		temp3 = ri(:)-intersection(:,3)
		if (magnitude(temp1).lt.magnitude(temp2)) then
			if (magnitude(temp1).lt.magnitude(temp3)) then
		!if (magnitude(ri(:)-intersection(:,1)) .lt.  & 
  	    !   	    magnitude(ri(:)-intersection(:,2))) then
		!	if (magnitude(ri(:)-intersection(:,1)) .lt.  & 
  		 !      	    magnitude(ri(:)-intersection(:,3))) then
				i = 1
			else
				i = 3
			endif
		else
			if (magnitude(temp2) .lt. magnitude(temp3)) then
			!if (magnitude(ri(:)-intersection(:,2)) .lt.  & 
  	       	!  	    magnitude(ri(:)-intersection(:,3))) then
				i = 2
			else
				i = 3
			endif
		endif

		!Check which plane is closest to j 
		temp1 = rj(:)-intersection(:,1)
		temp2 = rj(:)-intersection(:,2)
		temp3 = rj(:)-intersection(:,3)
		!if (magnitude(rj(:)-intersection(:,1)) .lt.  & 
  	    !   	    magnitude(rj(:)-intersection(:,2))) then
		!	if (magnitude(rj(:)-intersection(:,1)) .lt.  & 
  		!      	    magnitude(rj(:)-intersection(:,3))) then
		if (magnitude(temp1).lt.magnitude(temp2))then
			if (magnitude(temp1).lt.magnitude(temp3)) then
				j = 1
			else
				j = 3
			endif
		else
		!	if (magnitude(rj(:)-intersection(:,2)) .lt.  & 
  	    !  	  	    magnitude(rj(:)-intersection(:,3))) then
			if (magnitude(temp2).lt. magnitude(temp3)) then
				j = 2
			else
				j = 3
			endif
		endif

		!Calculate vectors from ri to intersect & rj to intersect
		Lfrac(1,1,1,:) = ri(:)-intersection(:,i)
		Lfrac(bindiff(1),bindiff(2),bindiff(3),:) = rj(:)-intersection(:,j)

		!Calculate vectors in two intermediate cells
		Lfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1),:) = intersection(:,l)-intersection(:,(6-k-l))
		Lfrac(interbindiff(1,2),interbindiff(2,2),interbindiff(3,2),:) = intersection(:,k)-intersection(:,(6-k-l))

		!Calculate magnitude of 3 vector components
		do ixyz = 1,3
			MLfrac(:,:,:) = MLfrac(:,:,:) + Lfrac(:,:,:,ixyz)**2
		enddo
		MLfrac(:,:,:) = MLfrac(:,:,:)**0.5d0

		!Normalise to one
		MLfrac(:,:,:) = MLfrac(:,:,:)/magnitude(rij(:))

		!------Add stress component to bin weighted by line length-----
		!Intermediate bin is either in domain or in halo. 
		!For halo bins the stress is added to the cell the halo represents. 
		!For domain cells it is added directly to that cell
		do ixyz = 1,nd
		do jxyz = 1,nd
			!-----------Molecule i bin-----------
			rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) = rfbin(ibin(1),ibin(2),ibin(3),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(1,1,1)

			!-----------1st Intermediate Bin-----------
			!Intermediate and i bin in domain, j in halo - add intermediate for molecule i
			rfbin(interbin(1,1),interbin(2,1),interbin(3,1),ixyz,jxyz) = &
				rfbin(interbin(1,1),interbin(2,1),interbin(3,1),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1))

			!-----------2nd Intermediate Bin-----------
			!Intermediate and i bin in domain, j in halo - add intermediate for molecule i
			rfbin(interbin(1,2),interbin(2,2),interbin(3,2),ixyz,jxyz) = &
				rfbin(interbin(1,2),interbin(2,2),interbin(3,2),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(interbindiff(1,2),interbindiff(2,2),interbindiff(3,2))

			!-----------Molecule j bin-----------
			rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) = rfbin(jbin(1),jbin(2),jbin(3),ixyz,jxyz) &
				+ accijmag*rij(ixyz)*rij(jxyz)*MLfrac(bindiff(1),bindiff(2),bindiff(3))
		enddo
		enddo

		deallocate(intersection)
		deallocate(Lfrac)
		deallocate(MLfrac)
		deallocate(interbin)
		deallocate(interbindiff)
		
	case default

	call error_abort("VOLUME AVERAGING ERROR")

	end select

end subroutine pressure_tensor_forces_VA

!===================================================================================
!Forces over the surface of a Volume

subroutine control_volume_forces(fij,ri,rj,molnoi,molnoj)
use module_record
implicit none

	integer							:: molnoi, molnoj
	integer,dimension(3)			:: ibin, jbin
	double precision,dimension(3)	:: ri, rj, fij,crossplane,fsurface
	double precision,dimension(3)	:: Fbinsize, bintopi, binboti, bintopj, binbotj

	!Determine bin size
	Fbinsize(:) = domain(:) / nbins(:)

	!Assign to bins using integer division
	ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:)) + nhb	!Establish current bin
	jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:)) + nhb 	!Establish current bin

	crossplane(:) =  dble(ibin(:)-jbin(:))

	bintopi(:) = (ibin(:)-1)*Fbinsize(:)-halfdomain(:)
	binboti(:) = (ibin(:)-2)*Fbinsize(:)-halfdomain(:)
	bintopj(:) = (jbin(:)-1)*Fbinsize(:)-halfdomain(:)
	binbotj(:) = (jbin(:)-2)*Fbinsize(:)-halfdomain(:)

	!Add for molecule i
	if(molnoi .le. np) then
		fsurface = fij(:)* dble((heaviside(bintopi(1)-ri(1))-heaviside(binboti(1)-ri(1)))* & 
			  		(heaviside(bintopi(2)-ri(2))-heaviside(binboti(2)-ri(2)))* & 
			  		(heaviside(bintopi(3)-ri(3))-heaviside(binboti(3)-ri(3)))- & 
			  		(heaviside(bintopi(1)-rj(1))-heaviside(binboti(1)-rj(1)))* & 
			  		(heaviside(bintopi(2)-rj(2))-heaviside(binboti(2)-rj(2)))* & 
			  		(heaviside(bintopi(3)-rj(3))-heaviside(binboti(3)-rj(3))))
		volume_force(ibin(1),ibin(2),ibin(3),:,1) = volume_force(ibin(1),ibin(2),ibin(3),:,1) + fsurface*delta_t
	endif

	!Add for molecule j
	if(molnoj .le. np) then
		fsurface = fij(:)* dble((heaviside(bintopj(1)-ri(1))-heaviside(binbotj(1)-ri(1)))* & 
			  		(heaviside(bintopj(2)-ri(2))-heaviside(binbotj(2)-ri(2)))* & 
			  		(heaviside(bintopj(3)-ri(3))-heaviside(binbotj(3)-ri(3)))- & 
			  		(heaviside(bintopj(1)-rj(1))-heaviside(binbotj(1)-rj(1)))* & 
			  		(heaviside(bintopj(2)-rj(2))-heaviside(binbotj(2)-rj(2)))* & 
			  		(heaviside(bintopj(3)-rj(3))-heaviside(binbotj(3)-rj(3))))
		volume_force(jbin(1),jbin(2),jbin(3),:,1) = volume_force(jbin(1),jbin(2),jbin(3),:,1) + fsurface*delta_t
	endif

end subroutine control_volume_forces

!===================================================================================
!Forces over the surface of a Volume

subroutine control_volume_stresses(fij,ri,rj,molnoi)
    use module_record
	use CV_objects, only : CV_debug,CVcheck_momentum2
    implicit none

	integer							:: i,j,k,ixyz,molnoi
	integer							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	integer,dimension(3)			:: cbin, ibin, jbin
	double precision,dimension(3)	:: ri,rj,rij,fij,fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb,velvect
	double precision,dimension(3)	:: Fbinsize, bintop, binbot

	!Calculate rij
	rij = ri - rj
	!Prevent Division by zero
	do ixyz = 1,3
		if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
	enddo

	!Determine bin size
	Fbinsize(:) = domain(:) / nbins(:)

	!Assign to bins using integer division
	ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
	jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin

	if (ibin(1) .eq. jbin(1) .and. ibin(2) .eq. jbin(2) .and. ibin(3) .eq. jbin(3)) return
		
	do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
	do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
	do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

		cbin(1) = i; cbin(2) = j; cbin(3) = k

		bintop(:) = (cbin(:)-1*nhb(:)  )*Fbinsize(:)-halfdomain(:)
		binbot(:) = (cbin(:)-1*nhb(:)-1)*Fbinsize(:)-halfdomain(:)

		!Calculate the plane intersect of line with surfaces of the cube
		Pxt=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
		Pxb=(/ binbot(1),ri(2)+(rij(2)/rij(1))*(binbot(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(binbot(1)-ri(1))  /)
		Pyt=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
		Pyb=(/ri(1)+(rij(1)/rij(2))*(binbot(2)-ri(2)), binbot(2),ri(3)+(rij(3)/rij(2))*(binbot(2)-ri(2))  /)
		Pzt=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)
		Pzb=(/ri(1)+(rij(1)/rij(3))*(binbot(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(binbot(3)-ri(3)), binbot(3) /)

		onfacexb =  	(sign(1.d0,binbot(1)- rj(1)) - sign(1.d0,binbot(1)- ri(1)))* &
						(heaviside(bintop(2)-Pxb(2)) - heaviside(binbot(2)-Pxb(2)))* &
						(heaviside(bintop(3)-Pxb(3)) - heaviside(binbot(3)-Pxb(3)))
		onfaceyb =  	(sign(1.d0,binbot(2)- rj(2)) - sign(1.d0,binbot(2)- ri(2)))* &
						(heaviside(bintop(1)-Pyb(1)) - heaviside(binbot(1)-Pyb(1)))* &
						(heaviside(bintop(3)-Pyb(3)) - heaviside(binbot(3)-Pyb(3)))
		onfacezb =  	(sign(1.d0,binbot(3)- rj(3)) - sign(1.d0,binbot(3)- ri(3)))* &
						(heaviside(bintop(1)-Pzb(1)) - heaviside(binbot(1)-Pzb(1)))* &
						(heaviside(bintop(2)-Pzb(2)) - heaviside(binbot(2)-Pzb(2)))

		onfacext =  	(sign(1.d0,bintop(1)- rj(1)) - sign(1.d0,bintop(1)- ri(1)))* &
						(heaviside(bintop(2)-Pxt(2)) - heaviside(binbot(2)-Pxt(2)))* &
	            		(heaviside(bintop(3)-Pxt(3)) - heaviside(binbot(3)-Pxt(3)))
		onfaceyt = 		(sign(1.d0,bintop(2)- rj(2)) - sign(1.d0,bintop(2)- ri(2)))* &
						(heaviside(bintop(1)-Pyt(1)) - heaviside(binbot(1)-Pyt(1)))* &
						(heaviside(bintop(3)-Pyt(3)) - heaviside(binbot(3)-Pyt(3)))
		onfacezt =  	(sign(1.d0,bintop(3)- rj(3)) - sign(1.d0,bintop(3)- ri(3)))* &
						(heaviside(bintop(1)-Pzt(1)) - heaviside(binbot(1)-Pzt(1)))* &
						(heaviside(bintop(2)-Pzt(2)) - heaviside(binbot(2)-Pzt(2)))

		!Stress acting on face over volume
		Pxyface(cbin(1),cbin(2),cbin(3),:,1) = Pxyface(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*dble(onfacexb)
		Pxyface(cbin(1),cbin(2),cbin(3),:,2) = Pxyface(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*dble(onfaceyb)
		Pxyface(cbin(1),cbin(2),cbin(3),:,3) = Pxyface(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*dble(onfacezb)
		Pxyface(cbin(1),cbin(2),cbin(3),:,4) = Pxyface(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*dble(onfacext)
		Pxyface(cbin(1),cbin(2),cbin(3),:,5) = Pxyface(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*dble(onfaceyt)
		Pxyface(cbin(1),cbin(2),cbin(3),:,6) = Pxyface(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*dble(onfacezt)

		!Add instantanous stress to CV record
		if (CV_debug) then
    		CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,1) = & 
				CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*dble(onfacexb)
    		CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,2) = & 
				CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*dble(onfaceyb)
    		CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,3) = & 
				CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*dble(onfacezb)
    		CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,4) = & 
				CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*dble(onfacext)
    		CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,5) = & 
				CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*dble(onfaceyt)
    		CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,6) = & 
				CVcheck_momentum2%Pxy(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*dble(onfacezt)
		endif

		!Stress acting on face over volume
		if (eflux_outflag .ne. 0) then
			velvect(:) = v(:,molnoi) 
			!velvect(:) = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)
			Pxyvface(cbin(1),cbin(2),cbin(3),1) = Pxyvface(cbin(1),cbin(2),cbin(3),1) + dot_product(fij,velvect)*dble(onfacexb)
			Pxyvface(cbin(1),cbin(2),cbin(3),2) = Pxyvface(cbin(1),cbin(2),cbin(3),2) + dot_product(fij,velvect)*dble(onfaceyb)
			Pxyvface(cbin(1),cbin(2),cbin(3),3) = Pxyvface(cbin(1),cbin(2),cbin(3),3) + dot_product(fij,velvect)*dble(onfacezb)
			Pxyvface(cbin(1),cbin(2),cbin(3),4) = Pxyvface(cbin(1),cbin(2),cbin(3),4) + dot_product(fij,velvect)*dble(onfacext)
			Pxyvface(cbin(1),cbin(2),cbin(3),5) = Pxyvface(cbin(1),cbin(2),cbin(3),5) + dot_product(fij,velvect)*dble(onfaceyt)
			Pxyvface(cbin(1),cbin(2),cbin(3),6) = Pxyvface(cbin(1),cbin(2),cbin(3),6) + dot_product(fij,velvect)*dble(onfacezt)
		endif

		!Force applied to volume
		fsurface(:) = 0.d0
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfacexb - onfacext)
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfaceyb - onfaceyt)
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfacezb - onfacezt)
		volume_force(cbin(1),cbin(2),cbin(3),:,1) = volume_force(cbin(1),cbin(2),cbin(3),:,1) + fsurface*delta_t

	enddo
	enddo
	enddo

end subroutine control_volume_stresses


!===================================================================================
!Forces over the surface of a Volume optmised for computational efficiency

subroutine control_volume_stresses_opt(fij,ri,rj,molnoi)
	use module_record
	implicit none

	integer							:: i,j,k,ixyz,molnoi,molnoj
	integer,dimension(3)			:: cbin, ibin, jbin, Si
	double precision,dimension(3)	:: ri,rj,rij,fij,fsurface,Px,Py,Pz,sgnjit,sgnjib,onfaceb,onfacet,velvect
	double precision,dimension(3)	:: Fbinsize, bintop, binbot

	!Calculate rij
	rij = ri - rj
	!Prevent Division by zero
	do ixyz = 1,3
		if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
	enddo

	!Determine bin size
	Fbinsize(:) = domain(:) / nbins(:)

	!Assign to bins using integer division
	ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
	jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin

	if (ibin(1) .eq. jbin(1) .and. ibin(2) .eq. jbin(2) .and. ibin(3) .eq. jbin(3)) return
		
	do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
	do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
	do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

		cbin(1) = i; cbin(2) = j; cbin(3) = k

		bintop(:) = (cbin(:)-1*nhb(:)  )*Fbinsize(:)-halfdomain(:)
		binbot(:) = (cbin(:)-1*nhb(:)-1)*Fbinsize(:)-halfdomain(:)

		!Calculate the plane intersect of line with surfaces of the cube
		Px=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
		Py=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
		Pz=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)

		sgnjit(:)= sign(1.d0,bintop(:)- rj(:)) - sign(1.d0,bintop(:)- ri(:))
		sgnjib(:)= sign(1.d0,binbot(:)- rj(:)) - sign(1.d0,binbot(:)- ri(:))

		Si(1) =	(heaviside(bintop(2)-Px(2)) - heaviside(binbot(2)-Px(2)))* &
				(heaviside(bintop(3)-Px(3)) - heaviside(binbot(3)-Px(3)))
		Si(2) =	(heaviside(bintop(1)-Py(1)) - heaviside(binbot(1)-Py(1)))* &
				(heaviside(bintop(3)-Py(3)) - heaviside(binbot(3)-Py(3)))
		Si(3) =	(heaviside(bintop(1)-Pz(1)) - heaviside(binbot(1)-Pz(1)))* &
				(heaviside(bintop(2)-Pz(2)) - heaviside(binbot(2)-Pz(2)))

		onfaceb(:) = sgnjib(:)*dble(Si(:))
		onfacet(:) = sgnjit(:)*dble(Si(:))

		!Stress acting on face over volume
		Pxyface(cbin(1),cbin(2),cbin(3),:,1) = Pxyface(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*onfaceb(1)
		Pxyface(cbin(1),cbin(2),cbin(3),:,2) = Pxyface(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*onfaceb(2)
		Pxyface(cbin(1),cbin(2),cbin(3),:,3) = Pxyface(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*onfaceb(3)
		Pxyface(cbin(1),cbin(2),cbin(3),:,4) = Pxyface(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*onfacet(1)
		Pxyface(cbin(1),cbin(2),cbin(3),:,5) = Pxyface(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*onfacet(2)
		Pxyface(cbin(1),cbin(2),cbin(3),:,6) = Pxyface(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*onfacet(3)

		!Stress acting on face over volume
		if (eflux_outflag .ne. 0) then
			velvect(:) = v(:,molnoi) 
			!if (molnoi .gt. np) print*, velvect(1)
			!velvect(:) = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)
			Pxyvface(cbin(1),cbin(2),cbin(3),1) = Pxyvface(cbin(1),cbin(2),cbin(3),1) + dot_product(fij,velvect)*onfaceb(1)
			Pxyvface(cbin(1),cbin(2),cbin(3),2) = Pxyvface(cbin(1),cbin(2),cbin(3),2) + dot_product(fij,velvect)*onfaceb(2)
			Pxyvface(cbin(1),cbin(2),cbin(3),3) = Pxyvface(cbin(1),cbin(2),cbin(3),3) + dot_product(fij,velvect)*onfaceb(3)
			Pxyvface(cbin(1),cbin(2),cbin(3),4) = Pxyvface(cbin(1),cbin(2),cbin(3),4) + dot_product(fij,velvect)*onfacet(1)
			Pxyvface(cbin(1),cbin(2),cbin(3),5) = Pxyvface(cbin(1),cbin(2),cbin(3),5) + dot_product(fij,velvect)*onfacet(2)
			Pxyvface(cbin(1),cbin(2),cbin(3),6) = Pxyvface(cbin(1),cbin(2),cbin(3),6) + dot_product(fij,velvect)*onfacet(3)
		endif

		!Force applied to volume
		fsurface(:) = 0.d0
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*(onfaceb(1) - onfacet(1))
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*(onfaceb(2) - onfacet(2))
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*(onfaceb(3) - onfacet(3))
		volume_force(cbin(1),cbin(2),cbin(3),:,1) = volume_force(cbin(1),cbin(2),cbin(3),:,1) + fsurface*delta_t

	enddo
	enddo
	enddo

end subroutine control_volume_stresses_opt


!===================================================================================
!Forces over the surface of a Volume further optmised for computational efficiency

!subroutine control_volume_stresses_opt_2(fij,ri,rj,molnoi)
!use module_record
!implicit none

!	integer							:: i,j,k,ixyz,molnoi
	!integer							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
!	integer,dimension(3)			:: cbin, ibin, jbin
!	integer,dimension(18)			:: hfacelimits
!	double precision,dimension(3)	:: ri,rj,rij,fij,fsurface,Px,Py,Pz,Si,sgnjit,sgnjib,onfaceb,onfacet,velvect
!	double precision,dimension(3)	:: Fbinsize, bintop, binbot
!	double precision,dimension(18)	:: facelimits

	!Calculate rij
!	rij = ri - rj
	!Prevent Division by zero
!	do ixyz = 1,3
!		if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
!	enddo

	!Determine bin size
!	Fbinsize(:) = domain(:) / nbins(:)

	!Assign to bins using integer division
!	ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!	jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!
!	if (ibin(1) .eq. jbin(1) .and. ibin(2) .eq. jbin(2) .and. ibin(3) .eq. jbin(3)) return
		
!	do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
!	do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
!	do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

!		cbin(1) = i; cbin(2) = j; cbin(3) = k

!		bintop(:) = (cbin(:)-1*nhb(:)  )*Fbinsize(:)-halfdomain(:)
!		binbot(:) = (cbin(:)-1*nhb(:)-1)*Fbinsize(:)-halfdomain(:)

		!Calculate the plane intersect of line with surfaces of the cube
!		Px=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
!		Py=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
!		Pz=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)

!		facelimits(1 :3 ) = bintop-Px	
!		facelimits(4 :6 ) = bintop-Py 		
!		facelimits(7 :9 ) = bintop-Pz
!		facelimits(10:12) = binbot-Px	
!		facelimits(13:15) = binbot-Py 		
!		facelimits(16:18) = binbot-Pz

!		hfacelimits = heaviside(facelimits)

!		Si(1) =	(hfacelimits(2) - hfacelimits(11))* &
!				(hfacelimits(3) - hfacelimits(12))
!		Si(2) =	(hfacelimits(4) - hfacelimits(13))* &
!				(hfacelimits(6) - hfacelimits(15))
!		Si(3) =	(hfacelimits(7) - hfacelimits(16))* &
!				(hfacelimits(8) - hfacelimits(17))

!		sgnjit(:)= sign(1.d0,bintop(:)- rj(:)) - sign(1.d0,bintop(:)- ri(:))
!		sgnjib(:)= sign(1.d0,binbot(:)- rj(:)) - sign(1.d0,binbot(:)- ri(:))

!		onfaceb =  	sgnjib*Si
!		onfacet =  	sgnjit*Si

		!Stress acting on face over volume
!		Pxyface(cbin(1),cbin(2),cbin(3),:,1) = Pxyface(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*dble(onfaceb(1))
!		Pxyface(cbin(1),cbin(2),cbin(3),:,2) = Pxyface(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*dble(onfaceb(2))
!		Pxyface(cbin(1),cbin(2),cbin(3),:,3) = Pxyface(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*dble(onfaceb(3))
!		Pxyface(cbin(1),cbin(2),cbin(3),:,4) = Pxyface(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*dble(onfacet(1))
!		Pxyface(cbin(1),cbin(2),cbin(3),:,5) = Pxyface(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*dble(onfacet(2))
!		Pxyface(cbin(1),cbin(2),cbin(3),:,6) = Pxyface(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*dble(onfacet(3))

!		!Stress acting on face over volume
!		if (eflux_outflag .ne. 0) then
!			velvect(:) = v(:,molnoi) 
			!if (molnoi .gt. np) print*, velvect(1)
			!velvect(:) = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)
!			Pxyvface(cbin(1),cbin(2),cbin(3),1) = Pxyvface(cbin(1),cbin(2),cbin(3),1) + dot_product(fij,velvect)*dble(onfaceb(1))
!			Pxyvface(cbin(1),cbin(2),cbin(3),2) = Pxyvface(cbin(1),cbin(2),cbin(3),2) + dot_product(fij,velvect)*dble(onfaceb(2))
!			Pxyvface(cbin(1),cbin(2),cbin(3),3) = Pxyvface(cbin(1),cbin(2),cbin(3),3) + dot_product(fij,velvect)*dble(onfaceb(3))
!			Pxyvface(cbin(1),cbin(2),cbin(3),4) = Pxyvface(cbin(1),cbin(2),cbin(3),4) + dot_product(fij,velvect)*dble(onfacet(1))
!			Pxyvface(cbin(1),cbin(2),cbin(3),5) = Pxyvface(cbin(1),cbin(2),cbin(3),5) + dot_product(fij,velvect)*dble(onfacet(2))
!			Pxyvface(cbin(1),cbin(2),cbin(3),6) = Pxyvface(cbin(1),cbin(2),cbin(3),6) + dot_product(fij,velvect)*dble(onfacet(3))
!		endif

		!Force applied to volume
!		fsurface(:) = 0.d0
!		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfaceb(1) - onfacet(1))
!		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfaceb(2) - onfacet(2))
!		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfaceb(3) - onfacet(3))
!		volume_force(cbin(1),cbin(2),cbin(3),:,1) = volume_force(cbin(1),cbin(2),cbin(3),:,1) + fsurface*delta_t

!	enddo
!	enddo
!	enddo

!end subroutine control_volume_stresses_opt_2

!====================================================================================

subroutine pressure_tensor_forces_MOP(pnxyz,ri,rj,rij,accijmag)
	use module_record
	implicit none

	integer							:: n
	integer							:: pnxyz	 !Plane normal direction
	integer							:: planenoi,planenoj
	double precision                :: shift, plane !Plane normal components i and j
	double precision                :: accijmag      !Non directional component of acceleration
	double precision,dimension(3)   :: ri, rj, rij   !Vector between particles i and j
	double precision,dimension(3)   :: Pyb           !Location of intercept with plane

	!Shift by half difference between value rounded down and actual value
	!to ensure same distance to top and bottom plane from domain edge
	shift = 0.5d0*(domain(pnxyz)-planespacing*(nplanes-1))

	!Obtain nearest plane number by integer division
	planenoi = nint((ri(2)+halfdomain(pnxyz)-shift)/planespacing)+1
	if (planenoi .lt. 1)	   planenoi = 1
	if (planenoi .gt. nplanes) planenoi = nplanes
	planenoj = nint((rj(2)+halfdomain(pnxyz)-shift)/planespacing)+1
	if (planenoj .lt. 1)	   planenoj = 1
	if (planenoj .gt. nplanes) planenoj = nplanes

	!Use calculated plane numbers and check all planes between
	do n = planenoi,planenoj,sign(1,planenoj-planenoi)
		plane = planes(n)

		!Calculate intersection with plane to check if interaction is within domain
		Pyb=(/ri(1)+(rij(1)/rij(2))*(plane-ri(2)), plane,ri(3)+(rij(3)/rij(2))*(plane-ri(2)) /)

		!Using sign function (Heaviside function is less efficient as no inbuilt Fortran Heaviside)
		Pxy_plane(:,n) =  Pxy_plane(:,n) + accijmag * rij(:) * & 
				( sign(1.d0,ri(2)   -    plane ) -      sign(1.d0,rj(2)  -  plane) )* &
				(heaviside(halfdomain(1)-Pyb(1)) - heaviside(-halfdomain(1)-Pyb(1)))* &
				(heaviside(halfdomain(3)-Pyb(3)) - heaviside(-halfdomain(3)-Pyb(3)))

	enddo

end subroutine pressure_tensor_forces_MOP

!===================================================================================
! Record external forces applied to molecules inside a volume

subroutine record_external_forces(F,ri)
	use module_record, only : domain,halfdomain, nbins, nhb, F_ext_bin
	implicit none

	double precision,dimension(3),intent(in):: F,ri

	integer	,dimension(3)					:: ibin
	double precision,dimension(3)			:: mbinsize

	!Determine bin size and bin
	mbinsize(:) = domain(:) / nbins(:)
	ibin(:) = ceiling((ri(:)+halfdomain(:))/mbinsize(:)) + nhb(:)

	!Add external force to bin
	print'(7i9,6f10.5)', shape(F_ext_bin), ibin, F, ri
	F_ext_bin(ibin(1),ibin(2),ibin(3),:) = & 
		F_ext_bin(ibin(1),ibin(2),ibin(3),:) + F(:)

end subroutine record_external_forces

!===================================================================================
!Calculate Radial distribution function (RDF) using cell method
!subroutine evaluate_properties_cellradialdist(molnoi)
!use module_record
!implicit none
!
!	integer                         :: j, ixyz   !Define dummy index
!	integer                         :: icell, jcell
!	integer                         :: icellshift, jcellshift
!	integer                         :: adjacentcellnp
!	integer							:: molnoi, molnoj
!	type(node), pointer 	        :: oldj, currentj

	!Find cell location of specified molecules
!	icell = ceiling((r(1,molnoi)+halfdomain(1))/cellsidelength(1))+1 !Add 1 due to halo
!	jcell = ceiling((r(2,molnoi)+halfdomain(2))/cellsidelength(2))+1 !Add 1 due to halo
!	ri = r(:,molnoi)

!	do jcellshift = -1,1
!	do icellshift = -1,1
!		oldj => cell%head(icell+icellshift,jcell+jcellshift)%point
!		adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift)
!
!		do j = 1,adjacentcellnp          !Step through all j for each i
!			rij2=0                   !Set rij^2 to zero
!			molnoj = oldj%molno !Number of molecule
!			rj = r(:,molnoj)         !Retrieve rj
!					
!			currentj => oldj
!			oldj => currentj%next    !Use pointer in datatype to obtain next item in list

!			if(molnoi==molnoj) cycle !Check to prevent interaction with self

!			do ixyz=1,nd
!				rij(ixyz) = ri(ixyz) - rj(ixyz)                 !Evaluate distance between particle i and j
!				if (abs(rij(ixyz)) > halfdomain(ixyz)) then  !N.B array operator-corresponding element compared 
!					rij(ixyz) = rij(ixyz) - sign(domain(ixyz),rij(ixyz))  !Sign of rij applied to domain
!				endif
!				rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
!			enddo
!		enddo
!	enddo
!	enddo

!	nullify(currentj)        !Nullify current as no longer required
!	nullify(oldj)            !Nullify old as no longer required

!end subroutine evaluate_properties_cellradialdist
