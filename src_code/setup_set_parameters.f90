!-----------------------------------------------------------------------------
!
!                               Set Parameters
! Set up domain size and cut off regions; allocate space for arrays
! and calculate initial velocity magnitudes based on temperature
!
! setup_set_parameters ! Main routine calling all setups
! set_parameters_allocate(1 or 2) !Allocate arrays in 2 parts given by intput
! set_parameters_global_domain
! set_parameters_cells
! set_parameters_setlimits
! set_parameters_outputs
! setup_linklist

!
!-----------------------------------------------------------------------------

module module_set_parameters 

	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use calculated_properties_MD
	use linked_list

end module module_set_parameters 
!----------------------------------------------------------------------------------

subroutine setup_set_parameters
	use module_set_parameters
	implicit none

	integer		:: i

	!Calculate shift in lennard-Jones potential based on cutoff
	potshift = 4*(rcutoff**-12 - rcutoff**-6)

	!Allocate arrays based on number of dimensions
	call set_parameters_allocate(1)

	!call set_parameters_domain
	call set_parameters_global_domain
	call set_parameters_cells
	call set_parameters_setlimits

	!Allocate array sizes for position, velocity and acceleration
	call set_parameters_allocate(2)

	!Zero quantities and arrays
	r = 0.d0
	v = 0.d0
	a = 0.d0
	fix = 1		!Set fix to one (unfixed)
	thermostat = 0	!Set all molecules to unthermostatted
	zeta= 0.d0	!Set Nose Hoover thermostat scaling property to zero
	rfmol = 0.d0
	halo_np = 0
	Pxyzero = 0.d0
	vmagnitude = 0.d0

	call set_parameters_outputs
	call setup_linklist
	call establish_surface_bins
	!call establish_surface_cells(-1)
	!call establish_surface_cells2(-1)
	!call establish_gpusurface_cells2(0)
	!call CUDA_setup

end subroutine setup_set_parameters


!===========================================================================================
!Allocate arrays first based on number of dimensions (n=1) then using extra allocation (n=2)
subroutine set_parameters_allocate(n)
	use module_set_parameters
	implicit none

	integer :: ixyz, n

	select case(n)
	case(1)
		!Allocate arrays based on number of dimensions
		allocate(initialunitsize(nd))
		allocate(domain(nd)) 
		allocate(halfdomain(nd))
		allocate(ncells(nd))
		allocate(cellsidelength(nd))
		allocate(halfcellsidelength(nd))
		!Pressure tensor
		allocate(Pxy(nd,nd))
		allocate(Pxyzero(nd,nd))
	case(2)
		!Calculate required extra allocation of molecules to allow copied Halos
		!using ratio of halo to domain volume with safety factor
		extralloc = 0
		do ixyz =1,nd
			extralloc = extralloc + &
			ceiling(((2.d0*(3*ncells(ixyz)**2+ &
			6*ncells(ixyz)+4))/ncells(ixyz)**3))*np
		enddo
		extralloc = extralloc/nd  + 300  !Average of all 3 dimensions inc safety factor

		!Allocate array sizes for position, velocity and acceleration
		allocate(r(np+extralloc,nd))
		allocate(rtrue(np+extralloc,nd)) !Used to establish diffusion - r with no periodic BC
		allocate(rinitial(np+extralloc,nd))
		allocate(rijsum(np+extralloc,nd)) !Sum of rij for each i, used for SLLOD algorithm
		allocate(v(np+extralloc,nd))
		allocate(vmagnitude(np+extralloc))
		allocate(a(np+extralloc,nd))

		!Allocate arrays use to fix molecules and allow sliding
		allocate(tag(np+extralloc))
		allocate(fix(np+extralloc,nd))
		allocate(slidev(np+extralloc,nd))
		allocate(thermostat(np+extralloc,nd))

		!Allocate potential energy and virial per molecule array
		allocate(potenergymol(np+extralloc))
		allocate(virialmol(np+extralloc))

		!Allocate pressure tensors
		allocate(rfmol(np+extralloc,nd,nd))
		allocate(Pxymol(np+extralloc,nd,nd))

	end select

end subroutine set_parameters_allocate

!-----------------------------------------------------------------------------

subroutine set_parameters_domain
	use module_set_parameters
	implicit none

	integer                :: ixyz

	np=1      !Set number of particles to unity for loop below
	volume=1  !Set domain size to unity for loop below
	do ixyz=1,nd
		domain(ixyz) = initialnunits(ixyz) &       !Size domain based on required density
		/((density/4)**(1.d0/nd)) 
		halfdomain(ixyz) = 0.5d0*domain(ixyz)      !Useful defintion used later
		np = np*initialnunits(ixyz)                !One particle per unit cell
		volume = volume*domain(ixyz)		   !Volume based on size of domain
	enddo
	if (nd .eq. 3) np=4*np   !FCC structure in 3D had 4 molecules per unit cell

	!Global number of particles
	globalnp = np
	call globalSumInt(globalnp)

end subroutine set_parameters_domain

!-----------------------------------------------------------------------------

subroutine set_parameters_global_domain
	use module_set_parameters
	implicit none

	integer                :: ixyz

	globalnp=1      !Set number of particles to unity for loop below
	volume=1	!Set domain size to unity for loop below
	do ixyz=1,nd
		globaldomain(ixyz) = initialnunits(ixyz) & 	!Size domain based on required density
		/((density/4.d0)**(1.d0/nd))
		globalnp = globalnp*initialnunits(ixyz)		!One particle per unit cell
		volume = volume*globaldomain(ixyz)		!Volume based on size of domain
	enddo

	globalnp=4*globalnp   !FCC structure in 3D had 4 molecules per unit cell

	!Initially assume molecules per processor are evenly split  - corrected after position setup
	np = globalnp / nproc					

	domain(1) = globaldomain(1) / npx			!determine domain size per processor
	domain(2) = globaldomain(2) / npy			!determine domain size per processor
	domain(3) = globaldomain(3) / npz			!determine domain size per processor

	do ixyz=1,nd
		halfdomain(ixyz) = 0.5d0*domain(ixyz)			!Useful definition
	enddo

	!Establish initial size of single unit to initialise microstate
	do ixyz=1,nd
		initialunitsize(ixyz) = globaldomain(ixyz) / initialnunits(ixyz)
	enddo

end subroutine set_parameters_global_domain

!-----------------------------------------------------------------------------------------

subroutine set_parameters_cells
	use module_set_parameters
	implicit none

	integer :: ixyz

	!Calculate size of neighbour list region
	rneighbr2 = (rcutoff + delta_rneighbr)**2

	!Calculate number of cells based on domain size and rcutoff rounding
	!down to give fewer cells but to ensure cells are all at least rcutoff
	do ixyz=1,nd
		ncells(ixyz)=floor(domain(ixyz)/(rcutoff+delta_rneighbr))
	enddo

	if (ncells(1)<3 .or. ncells(2)<3 .or. ncells(3)<3) then
		print*, ncells(1),'    in x and ', ncells(2), '    in y' , ncells(3), '    in z' 
		stop "WARNING - DOMAIN SHOULD HAVE AT LEAST 3 CELLS', &
		'IN X, Y AND Z - INCREASE NUMBER OF UNITS IN INPUT"
	endif

	!Determine side length of cells after rounding
	do ixyz=1,nd
		cellsidelength(ixyz) = domain(ixyz)/ncells(ixyz)
	enddo

end subroutine set_parameters_cells

!-----------------------------------------------------------------------------------------

subroutine set_parameters_setlimits
	use module_set_parameters
	implicit none

	nicellxl  = ncells(1)
	nicellyl  = ncells(2)
	nicellzl  = ncells(3)

	ncellxl  = nicellxl + 2 * nh
	ncellyl  = nicellyl + 2 * nh
	ncellzl  = nicellzl + 2 * nh

end subroutine set_parameters_setlimits

!----------------------------------------------------------------------------
!Nullify linklist head pointers and zero contents counts

subroutine setup_linklist
use module_set_parameters 
implicit none

	integer :: i
	integer :: icell, jcell, kcell
	integer :: ibin, jbin, kbin

	!Cell lists
	allocate(cell%head(ncells(1)+2,ncells(2)+2,ncells(3)+2))
	allocate(cell%cellnp(ncells(1)+2,ncells(2)+2,ncells(3)+2))
	do icell =1,ncells(1)+2
	do jcell =1,ncells(2)+2
	do kcell =1,ncells(3)+2
		cell%cellnp(icell,jcell,kcell)=0 !Zero number of molecules per cell
		nullify(cell%head(icell,jcell,kcell)%point) !Nullify cell head pointers
	enddo
	enddo
	enddo

	!Bin lists
	allocate(bin%head(1,nbins(1),1))
	allocate(bin%cellnp(1,nbins(1),1))
	do ibin =1,1
	do jbin =1,nbins(1)
	do kbin =1,1
		bin%cellnp(ibin,jbin,kbin)=0 !Zero number of molecules per cell
		nullify(bin%head(ibin,jbin,kbin)%point) !Nullify cell head pointers
	enddo
	enddo
	enddo

	!Passed lists
	pass%sendnp = 0			!Zero number of molecules in neighbour list
	nullify(pass%head)		!Nullify neighbour list head pointer 

end subroutine setup_linklist

!-----------------------------------------------------------------------------------------

subroutine set_parameters_outputs
	use module_set_parameters
	implicit none

	integer				:: n
	double precision		:: maxv !Maximum possible velocity of all molecules
	double precision		:: shift

	!Use definition of temperature and re-arrange to define an average velocity minus 3 degrees of
	!freedom - this is to fix the momentum of the domain boundaries 
	initialvel = sqrt(nd * (1.d0 - 1.d0/np)*inputtemperature)

	!Allocate bins used for calculating simulation properties
	globalnbins(1) = ncells(1) !ncells(1)*4*floor(rcutoff)!Total number of domain bins
	globalnbins(2) = ncells(2) !ncells(1)*4*floor(rcutoff)!Total number of domain bins
	globalnbins(3) = ncells(3) !ncells(1)*4*floor(rcutoff)!Total number of domain bins
	nbins(1) = nint(globalnbins(1)/dble(npy))	!Share global evenly between processes
	nbins(2) = nint(globalnbins(2)/dble(npy))	!Share global evenly between processes
	nbins(3) = nint(globalnbins(3)/dble(npy))	!Share global evenly between processes
	!Obtain global number of bins after rounding to given same number per process
	globalnbins(1) = nbins(1)
	globalnbins(2) = nbins(2)
	globalnbins(3) = nbins(3)
	if (slice_outflag .eq. 1) call SubcommSumInt(globalnbins(1),2)

	!print'(a,3i8,2(a,3f10.5))','nbins ',nbins,' binsize ', domain/nbins, ' Domain size ', domain

	!nbins(1) = ncells(2) 
	!nbins(1) = ceiling(np/10.d0)    	!Set number of equal sized velocity ranges based on 1/10 number of molecules
	allocate(vfd_bin(nbins(1)))           	!Allocate storage space for frequency tally over time
	allocate(normalisedvfd_bin(nbins(1))) 	!Allocate storage space for normalised frequency tally over time
	vfd_bin = 0 		       		!Set initial molecular frequency count to zero

	!Define maximum possible velocity of a given molecule to determine size of each bin
	maxv=initialvel*3.0            		!Assume molecule will not have more than 3 time its initial velocity 
	binsize = maxv/nbins(1)

	!Allocate and define number of shells used for Radial distribution function (RDF)
	nshells = ceiling(np / 2.d0)		!Define number of shells
	allocate(shell(nshells))		!Allocate array to tally positions
	allocate(RDF(nshells))		   	!Allocate array for radial distribution function
	shell = 0 				!Set initial molecular frequency count to zero

	!Define distance over which to evaluate RDF to determine delta_r	
	rd = 7.d0				   
	delta_r = rd/nshells			!Width of each shell

	!Allocate array for diffusion to number of dimensions
	allocate(diffusion(nd))
	allocate(meandiffusion(nd))
	diffusion = 0 !diffusion set to zero before sum over all molecules

	!Allocate pressure tensor correlation record length
	allocate(Pxycorrel(viscsample))

	!Allocated arrays for velocity slice
	allocate(meanvel(nbins(1)))
	allocate(countvel(nbins(1)))
	meanvel = 0.d0
	countvel = 0

	!Allocated bins for velocity averaging
	allocate(meanvelbin(nbins(1),nbins(2),nbins(3),3))
	allocate(countvelbin(nbins(1),nbins(2),nbins(3)))
	meanvelbin  = 0.d0
	countvelbin = 0

	!Allocated Nose Hoover local PUT thermstat bins
	allocate(zeta_array(nbins(1),nbins(2),nbins(3)))
	zeta_array = 0.d0
	!call local_temperature_header

	!Allocate pressure bin for Stress volume averaging
	allocate(rfbin(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,3))
	allocate(vvbin(nbins(1),nbins(2),nbins(3),3,3))
	allocate(Pxybin(nbins(1),nbins(2),nbins(3),3,3))
	allocate(Pxyface(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,6))
	rfbin  = 0.d0
	Pxybin = 0.d0

	!Allocate bins for control volume mass fluxes
	allocate(volume_mass(nbins(1)+2,nbins(2)+2,nbins(3)+2))
	allocate(  mass_flux(nbins(1)+2,nbins(2)+2,nbins(3)+2,6))
	volume_mass = 0
	mass_flux   = 0

	!Allocate bins for control volume momentum fluxes and forces
	allocate(volume_momentum(nbins(1)+2,nbins(2)+2,nbins(3)+2,3))
	allocate(momentum_flux(nbins(1)+2,nbins(2)+2,nbins(3)+2,6,3))
	allocate(volume_force(nbins(1)+2,nbins(2)+2,nbins(3)+2,3,2))
	momentum_flux = 0.d0
	volume_momentum = 0.d0
	volume_force = 0.d0

	!Allocated Bins for Nose Hoover Stress Control
	allocate(Gxybins(nbins(1),nbins(2),nbins(3),3,3))
	Gxybins = 0.d0

	!Allocate array for Stress Method of Planes
	planespacing = 0.1*cellsidelength(2)
	nplanes = floor(domain(2)/planespacing)

	!Shift by half difference between value rounded down and actual value
	!to ensure same distance to top and bottom plane from domain edge
	shift = 0.5d0*(domain(2) - planespacing * (nplanes-1))
	allocate(planes(nplanes))
	allocate(Pxy_plane(3,nplanes))

	!print*, 'number of planes',nplanes
	!print*, 'planespacing', planespacing, 'times nplanes =', & 
	!	planespacing*(nplanes),'shift',shift, 'domain',domain(2)

	!Setup location of planes
	do n = 1, nplanes
		planes(n) = planespacing*(n-1) + shift - halfdomain(2)
	enddo

	!Write header for velocity slice record
	if (slice_outflag .eq. 1) call parallel_slice_header
	if (slice_outflag .eq. 2) call velocity_average_header
	if (slice_outflag .eq. 3) then
		call parallel_slice_header
		call velocity_average_header
	endif

	!Open pressure tensor and viscosity record file 
	if (pressure_outflag .ne. 0) call stress_header
	call control_volume_header

end subroutine set_parameters_outputs

!-------------------------------------------------------------------
!Establish and store indices of cells which are on the outer domain surface specified 
!by size of buf
! buf = -1 Halo cells
! buf =  0 Domain outer surface cell 
! buf =  N Outer Surface-N

subroutine establish_surface_cells(buf)
	use module_set_parameters
	implicit none

	integer		:: n
	integer		:: icell, jcell, kcell
	integer		:: buf

	nsurfacecells =2*((ncells(1)-2*buf)*(ncells(2)-2*buf) &
			+  (ncells(1)-2-2*buf)*(ncells(3)-2*buf) &
		        +  (ncells(2)-2-2*buf)*(ncells(3)-2-2*buf))

	allocate(surfacecell(nsurfacecells,3))
	
	n = 1
	do kcell=1, ncells(1)+2
	do jcell=1, ncells(2)+2
	do icell=1, ncells(3)+2
		!Remove inner part of domain
		if((icell .gt. (2+buf) .and. icell .lt. (ncells(1)+1-buf)) .and. &
		   (jcell .gt. (2+buf) .and. jcell .lt. (ncells(2)+1-buf)) .and. &
		   (kcell .gt. (2+buf) .and. kcell .lt. (ncells(3)+1-buf))) cycle
		!Remove outer cells leaving only 1 layer of surface cells
		if((icell .lt. (2+buf) .or. icell .gt. (ncells(1)+1-buf)) .or. &
		   (jcell .lt. (2+buf) .or. jcell .gt. (ncells(2)+1-buf)) .or. &
		   (kcell .lt. (2+buf) .or. kcell .gt. (ncells(3)+1-buf))) cycle

		surfacecell(n,1)=icell
		surfacecell(n,2)=jcell
		surfacecell(n,3)=kcell
		n = n + 1
	enddo
	enddo
	enddo

	print*, 'Number of Surface Cells', nsurfacecells
	!do n = 1,nsurfacecells
	!	print*, 'surface cells',surfacecell(n,:)
	!enddo

end subroutine establish_surface_cells

!-------------------------------------------------------------------
!Establish and store indices of bins which are on the outer domain

subroutine establish_surface_bins
	use module_set_parameters
	implicit none

	integer		:: n
	integer		:: icell, jcell, kcell

	nsurfacecells =2*((nbins(1)+2)*(nbins(2)+2) &
			+  (nbins(1))*(nbins(3)+2) &
		        +  (nbins(2))*(nbins(3)))

	allocate(surfacecell(nsurfacecells,3))
	
	n = 1
	do kcell=1, nbins(1)+2
	do jcell=1, nbins(2)+2
	do icell=1, nbins(3)+2
		!Remove inner part of domain
		if((icell .gt. (1) .and. icell .lt. (nbins(1)+2)) .and. &
		   (jcell .gt. (1) .and. jcell .lt. (nbins(2)+2)) .and. &
		   (kcell .gt. (1) .and. kcell .lt. (nbins(3)+2))) cycle
		!Remove outer cells leaving only 1 layer of surface cells
		if((icell .lt. (1) .or. icell .gt. (nbins(1)+2)) .or. &
		   (jcell .lt. (1) .or. jcell .gt. (nbins(2)+2)) .or. &
		   (kcell .lt. (1) .or. kcell .gt. (nbins(3)+2))) cycle

		surfacecell(n,1)=icell
		surfacecell(n,2)=jcell
		surfacecell(n,3)=kcell
		n = n + 1
	enddo
	enddo
	enddo

	print*, 'Number of Surface Cells', nsurfacecells
	!do n = 1,nsurfacecells
	!	print*, 'surface cells',surfacecell(n,:)
	!enddo

end subroutine establish_surface_bins
