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

! *********CUDA ROUTINES - WILL BE IN SEPERATE FILE SOON************
! establish_surface_cells2(CPU_buffer_region_size)a
! establish_gpusurface_cells2(CPU_buffer_region_size-1)
! CUDA_setup
!
!
!-----------------------------------------------------------------------------

module module_set_parameters 

	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use calculated_properties_MD
	use linked_list
	use polymer_info_MD
	use concentric_cylinders

end module module_set_parameters 
!------------------------------------------------------------------------------

subroutine setup_set_parameters
	use module_set_parameters
	use interfaces, only : error_abort
	use librarymod, only : build_hilbert
	implicit none

	integer							:: i,iblk,jblk,kblk
	integer,dimension(3)			:: nblocks
	double precision,dimension(3)	:: blocksidelength

	!This has alreay been done in initialise for a coupled run
#if (USE_COUPLER == 0)	
   	call set_parameters_global_domain
	call set_parameters_cells
#endif
	!call set_parameters_setlimits

	!Calculate shift in lennard-Jones potential based on cutoff
	potshift = 4.d0*(1.d0/rcutoff**12 - 1.d0/rcutoff**6)

	!Calculate correction to lennard-Jones potential/pressure based on cutoff
	if (sLRC_flag .ne. 0) then
		potential_sLRC = 8.d0*pi*density      *(1.d0/(9.d0*rcutoff**9) - 1.d0/(3.d0*rcutoff**3))
		Pressure_sLRC  = 8.d0*pi*density**2.d0*(4.d0/(9.d0*rcutoff**9) - 2.d0/(3.d0*rcutoff**3))
	else
		potential_sLRC = 0.d0; Pressure_sLRC = 0.d0;
	endif

	!Allocate array sizes for position, velocity and acceleration
	call set_parameters_allocate

	!Zero quantities and arrays
	r = 0.d0
	v = 0.d0
	a = 0.d0

	zeta= 0.d0	!Set Nose Hoover thermostat scaling property to zero
	rfmol = 0.d0
	halo_np = 0

	call set_parameters_outputs
	call setup_linklist
	call establish_surface_bins
	call establish_halo_bins
	call establish_surface_cells
	call establish_halo_cells

	!Setup arrays for sorting algorithms
	select case(sort_flag)
	case(0)
		!No sort - do nothing
	case(1)
		!Build array of ordered cells
		blocksidelength = sortblocksize*cellsidelength
		nblocks = ceiling(domain/blocksidelength)
		allocate(Hcurve(nblocks(1),nblocks(2),nblocks(3)))
		do iblk=1,nblocks(1)
		do jblk=1,nblocks(2)
		do kblk=1,nblocks(3)
			Hcurve(iblk,jblk,kblk) = iblk + nblocks(1)*(jblk-1) & 
									      + nblocks(1)*nblocks(2)*(kblk-1)
		enddo
		enddo
		enddo
	case(2)
		!Hilbert curve between blocks of cells
		blocksidelength = sortblocksize*cellsidelength
		nblocks = ceiling(domain/blocksidelength)
		call build_hilbert(nblocks,Hcurve)
	case default
		call error_abort('Incorrect value of sort_flag')
	end select

	!call establish_surface_cells(-1)
	!call establish_surface_cells2(-1)
	!call establish_gpusurface_cells2(0)
	!call CUDA_setup

	!Setup shear info
	call setup_shear_parameters

	!Setup polymer info
	select case(potential_flag)
	case(0)
	case(1)
		call setup_polymer_info
	end select

	!Setup boundary forces to prevent molecular escape
	if (any(bforce_flag.ne.0)) then
		teval = 1
	else
		teval = tplot
	end if
#if USE_COUPLER
	teval = 1
#endif

	!Setup external forces applied to regions of the domain
	do i = 1,6
		!If F_ext_limits > domain extents set to domain extents
		if(F_ext_limits(i) .gt. 0.5d0*globaldomain(ceiling(real(i)/2.d0))) then
			!print*, 'b4', i,F_ext_limits(i)
			F_ext_limits(i) =  0.5d0*globaldomain(ceiling(real(i)/2.d0))
			!print*, 'at', i,F_ext_limits(i)
		elseif (F_ext_limits(i) .lt. -0.5d0*globaldomain(ceiling(real(i)/2.d0))) then
			!print*, 'b4', i,F_ext_limits(i)
			F_ext_limits(i) = -0.5d0*globaldomain(ceiling(real(i)/2.d0))
			!print*, 'at', i,F_ext_limits(i)
		endif

		!if(abs(F_ext_limits(i)) .gt. 0.5d0*globaldomain(ceiling(real(i)/2.d0))) then
		!	print*, 'b4', i,F_ext_limits(i)
		!	F_ext_limits(i) = ((-1)**i)*0.5d0*globaldomain(ceiling(real(i)/2.d0))
		!	print*, 'at', i,F_ext_limits(i)
		!endif
	enddo

end subroutine setup_set_parameters


!===========================================================================================
!Allocate arrays based on number of particles, np, using extra allocation

subroutine set_parameters_allocate
	use module_set_parameters
	use shear_info_MD
	implicit none

	integer :: ixyz

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
	allocate(r(nd,np+extralloc))
	allocate(v(nd,np+extralloc))
	allocate(a(nd,np+extralloc))

	!Allocate potential energy and virial per molecule array
	allocate(potenergymol(np+extralloc))
	allocate(potenergymol_LJ(np+extralloc))
	allocate(virialmol(np+extralloc))

	!Check if rtrue required
	if (r_gyration_outflag .eq. 1 .or. vmd_outflag       .eq. 4) rtrue_flag = 1
	if (rtrue_flag.eq.1) then
		allocate(rtrue(nd,np+extralloc)) !Used to establish diffusion - r with no periodic BC
		allocate(vtrue(nd,np+extralloc)) !Used to establish diffusion - r with no periodic BC
	endif

	!allocate(rijsum(nd,np+extralloc)) !Sum of rij for each i, used for SLLOD algorithm
	!allocate(vmagnitude(np+extralloc))

	!Arrays used for DPD thermostat
	if (ensemble .eq. nvt_DPD) then
		allocate(theta(nd,np+extralloc))
		allocate(aD(nd,np+extralloc))
		allocate(aR(nd,np+extralloc))
		call random_number(theta)
	endif

	!Allocate arrays use to fix molecules and allow sliding
	allocate(tag(np+extralloc)); tag = free
	if (ensemble .eq. tag_move) then
		allocate(fix(nd,np+extralloc)); fix = 1	!default - set fix to one (unfixed)
		allocate(rtether(nd,np+extralloc))
		allocate(slidev(nd,np+extralloc))
	endif

	!Allocate pressure tensors
	if (pressure_outflag .eq. 1) then
		allocate(rfmol(np+extralloc,nd,nd))
		allocate(Pxymol(np+extralloc,nd,nd))
	endif

	!Allocate bulk shear array
	if (any(periodic.eq.2)) then
		allocate(mol_wrap_integer(np))
	endif

end subroutine set_parameters_allocate

!-----------------------------------------------------------------------------
subroutine setup_polymer_info
	use module_set_parameters
	use polymer_info_MD
    use interfaces
	implicit none
	
	!Allocate polymer arrays
	allocate(bond(max_funcy,np+extralloc))
	bond = 0
	allocate(bondcount(np+extralloc))
	bondcount = 0
	allocate(monomer(np+extralloc))
	allocate(potenergymol_FENE(np+extralloc))

	etevtcf = 0.d0
	R_g     = 0.d0
	
	if (iter .gt. etevtcf_iter0)    etevtcf_iter0    = iter
	if (iter .gt. r_gyration_iter0) r_gyration_iter0 = iter

	intbits = bit_size(monomer(1)%bin_bflag(1))

end subroutine setup_polymer_info

!-----------------------------------------------------------------------------
subroutine setup_shear_parameters
use interfaces
use module_set_parameters
use shear_info_MD
implicit none

	integer :: i

	! DEFAULTS	
	le_sp = 2
	le_sd = 1
	le_rp = 3
	le_st = 0.d0
	le_sx = 0.d0	

	! CHANGE FROM DEFAULT VALUES
	do i=1,nd
		if (periodic(i).eq.2) le_sp = i
	end do

	if (any(periodic.gt.1)) then	
		select case(le_sp + le_sd)
		case(3)
			le_rp = 3
		case(4)
			le_rp = 2
		case(5)
			le_rp = 1
		case default
			call error_abort('Shear plane and shear direction must be different and 1,2 or 3')
		end select
	else
		le_sv = 0.d0
		le_sr = 0.d0
	end if

	if (define_shear_as.eq.0) le_sr = le_sv/domain(le_sd)
	if (define_shear_as.eq.1) le_sv = le_sr*domain(le_sd)
	
	if (iter .lt. le_i0) then
		le_st = 0.d0
		le_sx = 0.d0
	else
	!	le_i0 = iter
	end if

end subroutine setup_shear_parameters

!-----------------------------------------------------------------------------
!Setup domain based on density and number of initial units specified

subroutine set_parameters_global_domain
	use module_set_parameters
	use interfaces, only: error_abort
	implicit none

	integer                :: ixyz, extranp
	real(kind(0.d0))       :: lat

	select case (initial_config_flag)
	case (0)

		volume=1	!Set domain size to unity for loop below
		do ixyz=1,nd
			globaldomain(ixyz) = initialnunits(ixyz) & 	!Size domain based on required density
			/((density/4.d0)**(1.d0/nd))
			volume = volume*globaldomain(ixyz)		!Volume based on size of domain
		enddo

		! no need to fix globalnp if we have it already
		if(.not. restart) then
			globalnp=1      !Set number of particles to unity for loop below
			do ixyz=1,nd
				globalnp = globalnp*initialnunits(ixyz)		!One particle per unit cell
			enddo
			globalnp=4*globalnp   !FCC structure in 3D had 4 molecules per unit cell
		endif

		!Initially assume molecules per processor are evenly split  - corrected after position setup
		np = globalnp / nproc					

		domain(1) = globaldomain(1) / real(npx, kind(0.d0))			!determine domain size per processor
		domain(2) = globaldomain(2) / real(npy, kind(0.d0))			!determine domain size per processor
		domain(3) = globaldomain(3) / real(npz, kind(0.d0))			!determine domain size per processor

		do ixyz=1,nd
			halfdomain(ixyz) = 0.5d0*domain(ixyz)			!Useful definition
		enddo

		!Establish initial size of single unit to initialise microstate
		do ixyz=1,nd
			initialunitsize(ixyz) = globaldomain(ixyz) / initialnunits(ixyz)
		enddo

	case (1)

		select case (config_special_case)
		case ('sparse_fene')

			volume = product(globaldomain(1:3))
			if (.not.restart) then
				globalnp = nmonomers*nchains
			end if
			domain(1) = globaldomain(1) / real(npx, kind(0.d0))
			domain(2) = globaldomain(2) / real(npy, kind(0.d0))
			domain(3) = globaldomain(3) / real(npz, kind(0.d0))
			halfdomain(:) = 0.5d0*domain(:)

			!Initially assume molecules per processor are evenly split 
			! - corrected after position setup
			np = globalnp / nproc
		
		case ('dense_fene')
			
			globaldomain(:) = initialnunits(:)/((density/4.d0)**(1.d0/nd))
			initialunitsize(:) = globaldomain(:) / initialnunits(:)
			volume = product(globaldomain(1:3))

			! no need to fix globalnp if we have it already
			if(.not. restart) then
				globalnp = 4*product(initialnunits(1:3))
			endif

			domain(1) = globaldomain(1) / real(npx, kind(0.d0))
			domain(2) = globaldomain(2) / real(npy, kind(0.d0))
			domain(3) = globaldomain(3) / real(npz, kind(0.d0))
			halfdomain(:) = 0.5d0*domain(:)

			!Initially assume molecules per processor are evenly split 
			! - corrected after position setup
			np = globalnp / nproc
			
		case ('concentric_cylinders')

			globaldomain(:) = initialnunits(:)/((cyl_density/4.d0)**(1.d0/nd))
			volume = product(globaldomain(1:3))

			!Initially assume molecules per processor are evenly split 
			! - corrected after position setup
			globalnp = 4*product(initialnunits(1:3))
			np = globalnp / nproc					

			domain(1) = globaldomain(1) / real(npx, kind(0.d0))	
			domain(2) = globaldomain(2) / real(npy, kind(0.d0))
			domain(3) = globaldomain(3) / real(npz, kind(0.d0))
			halfdomain(:) = 0.5d0*domain(:)

			!Establish initial size of single unit to initialise microstate
			initialunitsize(:) = globaldomain(:) / initialnunits(:)
	
			specular_flag = specular_radial	

			r_oo = 0.5*initialunitsize(1)*cyl_units_oo
			r_io = 0.5*initialunitsize(1)*cyl_units_io
			r_oi = 0.5*initialunitsize(1)*cyl_units_oi
			r_ii = 0.5*initialunitsize(1)*cyl_units_ii

			! turn off periodicity in x and y direction
			periodic(1) = 0
			periodic(2) = 0
			periodic(3) = 1
		
		case ('fill_cylinders')

			!Read globalnp, globaldomain and r_oo etc
			call parallel_io_cyl_footer(cyl_file)


			volume = product(globaldomain(1:3))
			! Get closest number of FCC units in global domain from input
			initialnunits(:) = floor(globaldomain(:)*((density/4.d0)**(1.d0/nd)))
			! Get actual density with this integer number of FCC units
			density = product(initialnunits(1:3))*4 / volume
			! Guess how many extra molecules will fill the cylinder volume
			! so that arrays are allocated with enough space to add fluid
			extranp = int(ceiling(density * &
			                      globaldomain(3)*pi*(r_io**2.0-r_oi**2.0)))
			initialunitsize(:) = globaldomain(:) / initialnunits(:)

			domain(1) = globaldomain(1) / real(npx, kind(0.d0))	
			domain(2) = globaldomain(2) / real(npy, kind(0.d0))
			domain(3) = globaldomain(3) / real(npz, kind(0.d0))
			halfdomain(:) = 0.5d0*domain(:)

			!Initially assume molecules per processor are evenly split 
			! - corrected after position setup
			globalnp = cyl_np + extranp
			np = globalnp / nproc

			! turn off periodicity in x and y direction
			periodic(1) = 0
			periodic(2) = 0
			periodic(3) = 1

			! Ensure no rotation on setup	
			omega = 0.d0
		
			!print('(a,3f12.3,a,3f12.3,a,i10,a,i4,a,3l4,a,f12.3,a,f12.3)'), &
			!'globaldomain',globaldomain,' domain',domain,'np', np, 'nproc',&
			! nproc, 'periodic', periodic, 'density', density, 'volume', volume

		case ('rotate_cylinders')

			if (.not.restart) then
				stop 'Rotate cylinders requires a restart file.'
			end if

			!Read globalnp, globaldomain and r_oo etc
			call parallel_io_cyl_footer(cyl_file)

			volume = product(globaldomain(1:3))
			initialnunits(:) = globaldomain(:)*((density/4.d0)**(1.d0/nd))
			initialunitsize(:) = globaldomain(:) / initialnunits(:)

			domain(1) = globaldomain(1) / real(npx, kind(0.d0))
			domain(2) = globaldomain(2) / real(npy, kind(0.d0))
			domain(3) = globaldomain(3) / real(npz, kind(0.d0))
			halfdomain(:) = 0.5d0*domain(:)

			!Initially assume molecules per processor are evenly split 
			! - corrected after position setup
			np = globalnp / nproc

			! turn off periodicity in x and y direction
			periodic(1) = 0
			periodic(2) = 0
			periodic(3) = 1

			!print('(a,3f12.3,a,3f12.3,a,i10,a,i4,a,3l4,a,f12.3,a,f12.3)'), &
			!'globaldomain',globaldomain,' domain',domain,'np', np, 'nproc',&
			! nproc, 'periodic', periodic, 'density', density, 'volume', volume

		case default

			call error_abort("set_parameters_global_domain must be corrected for this special case")	

		end select

	case default

		call error_abort("Unrecognised initial_config_flag in set_parameters_global_domain")

	end select

end subroutine set_parameters_global_domain

!-----------------------------------------------------------------------------------------

subroutine set_parameters_cells
	use interfaces
	use module_set_parameters
	use polymer_info_MD
	implicit none

	integer :: ixyz
	double precision :: rneighbr

	!Calculate size of neighbour list region
	rneighbr = rcutoff + delta_rneighbr
	rneighbr2 = (rcutoff + delta_rneighbr)**2

	select case(potential_flag)
	case(1)
		select case(solvent_flag)
		case(0:1)
			if (rneighbr < R_0) then
				rneighbr = R_0 
				rneighbr2 = R_0**2
				if(irank.eq.iroot) print*, 'Neighbour list distance rneighbr set to &
						& maximum elongation of polymer spring, ',R_0
			end if
		case(2)
			if (rneighbr < sod_cut) then
				rcutoff   = sod_cut
				rcutoff2  = sod_cut2
				rneighbr  = rcutoff + delta_rneighbr
				rneighbr2 = rneighbr**2.d0
			end if
		case default
			call error_abort('ERROR - Unrecognised solvent_flag in set_parameters_cells')
		end select
	case default
	end select

	!Calculate number of cells based on domain size and rcutoff rounding
	!down to give fewer cells but to ensure cells are all at least rcutoff

	do ixyz=1,nd
		ncells(ixyz)=floor(domain(ixyz)/(rcutoff+delta_rneighbr))
	enddo

	if (ncells(1)<3 .or. ncells(2)<3 .or. ncells(3)<3) then
		print*, 'NCELLS:'
		print*, ncells(1),'    in x and ', ncells(2), '    in y' , ncells(3), '    in z' 
		call  error_abort( "ERROR - DOMAIN SHOULD HAVE AT LEAST 3 CELLS, &
		 					& IN X, Y AND Z - INCREASE NUMBER OF UNITS IN INPUT")
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
	use interfaces
	implicit none

	integer					:: n
	double precision		:: maxv !Maximum possible velocity of all molecules
	double precision		:: shift

	!Use definition of temperature and re-arrange to define an average velocity minus 3 degrees of
	!freedom - this is to fix the momentum of the domain boundaries 
	initialvel = sqrt(nd * (1.d0 - 1.d0/np)*inputtemperature)

	!Allocate bins used for calculating simulation properties
	gnbins(1) = npx*ncells(1) !Total number of domain bins
 	gnbins(2) = npy*ncells(2) !Total number of domain bins
	gnbins(3) = npz*ncells(3) !Total number of domain bins

	nbins(1) = nint(gnbins(1)/dble(npx))	!Share global evenly between processes
	nbins(2) = nint(gnbins(2)/dble(npy))	!Share global evenly between processes
	nbins(3) = nint(gnbins(3)/dble(npz))	!Share global evenly between processes

	!Obtain global number of bins after rounding to given same number per process
	gnbins(1) = nbins(1)
	gnbins(2) = nbins(2)
	gnbins(3) = nbins(3)
	call SubcommSum(gnbins(1),1)	!Sum up over all x processes
	call SubcommSum(gnbins(2),2)	!Sum up over all y processes
	call SubcommSum(gnbins(3),3)	!Sum up over all z processes

	!Calculate number of halo bins from ratio of cells to bins
	nhb = ceiling(dble(nbins)/dble(ncells))
	nbinso = nbins+2*nhb

	!Velocity binning routines
	!nbins(1) = ceiling(np/10.d0)    	!Set number of equal sized velocity ranges based on 1/10 number of molecules
	if (vdist_flag .eq. 1) then
		allocate(vfd_bin(nbins(1)))           	!Allocate storage space for frequency tally over time
		allocate(normalisedvfd_bin(nbins(1))) 	!Allocate storage space for normalised frequency tally over time
		vfd_bin = 0 		       		!Set initial molecular frequency count to zero

		!Define maximum possible velocity of a given molecule to determine size of each bin
		maxv=initialvel*3.0            		!Assume molecule will not have more than 3 time its initial velocity 
		binsize = maxv/nbins(1)
	endif

	!Allocate and define number of shells used for Radial distribution function (rdf)
	if (rdf_outflag .eq. 1) then
		allocate(rdf(rdf_nbins))                        !Allocate array for radial distribution function
		allocate(rdf_hist(rdf_nbins))                   !Allocate array to tally positions
		rdf= 0.d0
		rdf_hist= 0
	elseif(rdf_outflag .eq. 2) then
		allocate(rdf3d(rdf_nbins,nd))                   !Allocate array for radial distribution function
		allocate(rdf3d_hist(rdf_nbins,nd))              !Allocate array to tally positions
		rdf3d_hist= 0
		rdf3d= 0.d0
	endif
	!Allocate and define arrays for static structure factor
	if (ssf_outflag .eq. 1) then
		allocate(ssf_hist(2*ssf_nmax+1,2*ssf_nmax+1))   !Allocate array to tally positions
		allocate(ssf(2*ssf_nmax+1,2*ssf_nmax+1))        !Allocate array for radial distribution function
		ssf= 0.d0
		ssf_hist= 0.d0
	endif

	!Allocate array for diffusion to number of dimensions
	!allocate(diffusion(nd))
	!allocate(meandiffusion(nd))
	!diffusion = 0 !diffusion set to zero before sum over all molecules

	!Allocate pressure tensor correlation record length
	if(viscosity_outflag .eq. 1) then
		allocate(Pxycorrel(Nstress_ave))
		Pxycorrel = 0.d0
	endif

	!Allocated arrays for velocity slice
	if (velocity_outflag.ne.0 .and. velocity_outflag.lt.4) then
		allocate(slice_momentum(nbins(velocity_outflag),3))
		allocate(slice_mass(nbins(velocity_outflag)))
		slice_momentum = 0.d0
		slice_mass = 0
	else
		if (mass_outflag.ne.0 .and. mass_outflag.lt.4) then
			allocate(slice_mass(nbins(mass_outflag)))
			slice_mass = 0
		endif
	endif

	!Allocated bins for velocity averaging
	if (velocity_outflag.eq.4) then
		allocate(volume_momentum(nbinso(1),nbinso(2),nbinso(3),3  ))
		volume_momentum = 0.d0
		mass_outflag = 4	!Mass binning required too
	endif

	!Allocated bins for temperature averaging
	if (temperature_outflag.eq.4) then
		allocate(volume_temperature(nbinso(1),nbinso(2),nbinso(3)))
		volume_temperature = 0.d0
		mass_outflag = 4	!Mass binning required too
		!Allocate and zero peculiar momentum binning array
		if (peculiar_flag .ne. 0) then
			allocate(u(nd,np+extralloc)); u = 0.d0
			if (velocity_outflag.ne.4) then
				call error_abort("set_parameters_outputs Error -- Temperature outflag on with &
								 &perculiar momentum but velocity binning is off. Please switch &
								 &VELOCITY_OUTFLAG 1st option to 4 or TEMPERATURE_OUTFLAG 3rd option to 0")
			endif
		endif
	endif

	!Allocate mass bins if they haven't been already allocated (and they are needed)
	if (mass_outflag.eq.4) then
		if (.not. allocated(volume_mass)) then
			allocate(volume_mass(nbinso(1),nbinso(2),nbinso(3)))
			volume_mass = 0
		endif
	endif

	! Allocate cylindrical polar bins

	if (config_special_case.eq.'rotate_cylinders') then

		! No halos in r or theta
		cpol_bins(1) = gcpol_bins(1) ! r and theta bins are stored globally...
		cpol_bins(2) = gcpol_bins(2) ! on all processors in a z-plane
		cpol_bins(3) = nint(gcpol_bins(3)/dble(npz))

		!Obtain global number of bins after rounding to given same number per process
		gcpol_bins(3) = cpol_bins(3)
		call SubcommSum(gcpol_bins(3),3) !Sum up over all z processes

		! Only add halos to allocation in z direction	
		cpol_binso = cpol_bins
		cpol_nhbz = ceiling(dble(cpol_bins(3))/dble(ncells(3)))
		cpol_binso(3) = cpol_bins(3) + 2*cpol_nhbz

		if (velocity_outflag.eq.5) then
			allocate(cyl_mom(cpol_binso(1),cpol_binso(2),cpol_binso(3),3))
			allocate(cyl_mass(cpol_binso(1),cpol_binso(2),cpol_binso(3)))
			cyl_mom  = 0.d0
			cyl_mass = 0
		else if (mass_outflag.eq.5) then
			allocate(cyl_mass(cpol_binso(1),cpol_binso(2),cpol_binso(3)))
			cyl_mass = 0
		end if
	
	end if

	!Allocated Nose Hoover local PUT thermstat bins
	!allocate(zeta_array(nbins(1),nbins(2),nbins(3)))
	!zeta_array = 0.d0
	!call local_temperature_header

	!Pressure tensor
	allocate(Pxy(nd,nd))
	allocate(Pxyzero(nd,nd))
	Pxy = 0.d0
	Pxyzero = 0.d0
	if (pressure_outflag .eq. 2) then
		!Allocate pressure bin for Stress volume averaging
		allocate( rfbin( nbinso(1), nbinso(2),nbinso(3),3,3))
		allocate( vvbin( nbins (1),  nbins(2),  nbins(3),3,3  ))
		allocate( Pxybin(nbins (1),  nbins(2),  nbins(3),3,3  ))
		rfbin  = 0.d0
		Pxybin = 0.d0
	endif

	!Allocated Bins for Nose Hoover Stress Control
	!allocate(Gxybins(nbins(1),nbins(2),nbins(3),3,3))
	!Gxybins = 0.d0

	!Allocate array for Stress Method of Planes and/or 
	!allocate bins for control volume momentum fluxes and forces
	planespacing = cellsidelength(2)
	select case(vflux_outflag)
		case(1)
			gnplanes = npx*floor(domain(1)/planespacing)
			nplanes  = nint(gnplanes/dble(npx))
			gnplanes = nplanes
			call SubcommSum(gnplanes,1)	!Sum up over all x processes
			!Shift by half difference between value rounded down and actual value
			!to ensure same distance to top and bottom plane from domain edge
			shift = 0.5d0*(domain(1) - planespacing * (nplanes-1))
			allocate(planes(nplanes))
			allocate(Pxy_plane(3,nplanes))
			!Setup location of planes
			do n = 1, nplanes
				planes(n) = planespacing*(n-1) + shift - halfdomain(1)
			enddo
		case(2)
			gnplanes = npy*floor(domain(2)/planespacing)
			nplanes  = nint(gnplanes/dble(npy))
			gnplanes = nplanes
			call SubcommSum(gnplanes,2)	!Sum up over all y processes
			!Shift by half difference between value rounded down and actual value
			!to ensure same distance to top and bottom plane from domain edge
			shift = 0.5d0*(domain(2) - planespacing * (nplanes-1))
			allocate(planes(nplanes))
			allocate(Pxy_plane(3,nplanes))
			!Setup location of planes
			do n = 1, nplanes
				planes(n) = planespacing*(n-1) + shift - halfdomain(2)
			enddo
		case(3)
			gnplanes = npz*floor(domain(3)/planespacing)
			nplanes  = nint(gnplanes/dble(npz))
			gnplanes = nplanes
			call SubcommSum(gnplanes,3)	!Sum up over all z processes
			!Shift by half difference between value rounded down and actual value
			!to ensure same distance to top and bottom plane from domain edge
			shift = 0.5d0*(domain(3) - planespacing * (nplanes-1))
			allocate(planes(nplanes))
			allocate(Pxy_plane(3,nplanes))
			!Setup location of planes
			do n = 1, nplanes
				planes(n) = planespacing*(n-1) + shift - halfdomain(3)
			enddo
		case(4)
			if (.not.(allocated(volume_momentum))) 	allocate(volume_momentum(nbinso(1),nbinso(2),nbinso(3),3  ))
			allocate( Pxyface(nbinso(1),nbinso(2),nbinso(3),3,6))
			allocate(  momentum_flux(nbinso(1),nbinso(2),nbinso(3),3,6))
			allocate(   volume_force(nbinso(1),nbinso(2),nbinso(3),3,2))
			if (external_force_flag .ne. 0 .or. ensemble .eq. tag_move) then
				allocate(F_ext_bin(nbinso(1),nbinso(2),nbinso(3),3))
				F_ext_bin = 0.d0
			endif
			momentum_flux 	= 0.d0
			volume_momentum = 0.d0
			volume_force 	= 0.d0
			!Allocate bins for control volume mass fluxes
			if (.not.(allocated(volume_mass)))  allocate(volume_mass(nbinso(1),nbinso(2),nbinso(3)  ))
			allocate(  mass_flux(nbinso(1),nbinso(2),nbinso(3),6))
			volume_mass = 0
			mass_flux   = 0
		case default
			!Allocate bins for control volume mass fluxes
			if (mflux_outflag .eq. 1) then
				if (.not. allocated(volume_mass)) &
				allocate(volume_mass(nbinso(1),nbinso(2),nbinso(3)  ))
				allocate(  mass_flux(nbinso(1),nbinso(2),nbinso(3),6))
				volume_mass = 0
				mass_flux   = 0
			endif
	end select

	!Allocate bins for control volume energy fluxes and forces*velocity
	if (eflux_outflag .eq. 4) then
		allocate(  energy_flux(nbinso(1),nbinso(2),nbinso(3),6))
		allocate( Pxyvface(nbinso(1),nbinso(2),nbinso(3),6))
		energy_flux 	= 0.d0; Pxyvface = 0.d0
	endif

#if USE_COUPLER
	! Check end of maximum VMD intervals is not greater than the number of steps Nsteps
	! which has been changed by the coupler.
	!Nsteps = initialstep + coupler_md_get_nsteps() * coupler_md_get_md_steps_per_dt_cfd()
	if (vmd_outflag .ne. 0 .and. Nvmd_intervals .ne. 0) then
		if (maxval(vmd_intervals) .gt. Nsteps) then
			print'(2(a,i8))', 'Value specified for end of final vmd_interval = ' & 
							, maxval(vmd_intervals), ' but Nsteps = ', Nsteps 
			call error_abort("Specified VMD interval greater than Nsteps")
		endif
	endif
#endif

end subroutine set_parameters_outputs

!-------------------------------------------------------------------
!Establish and store indices of cells which are on the outer domain

subroutine establish_surface_cells
	use module_set_parameters
	implicit none

	integer		:: n
	integer		:: icell, jcell, kcell

	nsurfacecells=	2*( ncells(1)   * ncells(2) &
					+  (ncells(3)-2)* ncells(2) &
		        	+  (ncells(3)-2)*(ncells(1)-2))

	allocate(surfacecells(nsurfacecells,3))

	n = 1
	do kcell=1, ncells(3)+2
	do jcell=1, ncells(2)+2
	do icell=1, ncells(1)+2

		!Remove inner part of domain
		if((icell .gt. (2) .and. icell .lt. (ncells(1)+1)) .and. &
		   (jcell .gt. (2) .and. jcell .lt. (ncells(2)+1)) .and. &
		   (kcell .gt. (2) .and. kcell .lt. (ncells(3)+1))) cycle
		!Remove outer cells leaving only 1 layer of surface cells
		if((icell .lt. (2) .or. icell .gt. (ncells(1)+1)) .or. &
		   (jcell .lt. (2) .or. jcell .gt. (ncells(2)+1)) .or. &
		   (kcell .lt. (2) .or. kcell .gt. (ncells(3)+1))) cycle

		surfacecells(n,1)=icell
		surfacecells(n,2)=jcell
		surfacecells(n,3)=kcell
		n = n + 1

	enddo
	enddo
	enddo

end subroutine establish_surface_cells

!-------------------------------------------------------------------
!Establish and store indices of cells which are in the halo

subroutine establish_halo_cells
	use module_set_parameters
	implicit none

	integer		:: n
	integer		:: icell, jcell, kcell

	nhalocells  =	2*((ncells(1)+2)*(ncells(2)+2) &
					+  (ncells(3)  )*(ncells(2)+2) &
					+  (ncells(3)  )*(ncells(1)  ))

	allocate(halocells(nhalocells,3))

	n = 1
	do kcell=1, ncells(3)+2
	do jcell=1, ncells(2)+2
	do icell=1, ncells(1)+2

		!Remove inner part of domain
		if((icell .gt. (1) .and. icell .lt. (ncells(1)+2)) .and. &
		   (jcell .gt. (1) .and. jcell .lt. (ncells(2)+2)) .and. &
		   (kcell .gt. (1) .and. kcell .lt. (ncells(3)+2))) cycle

		halocells(n,1)=icell
		halocells(n,2)=jcell
		halocells(n,3)=kcell
		n = n + 1

	enddo
	enddo
	enddo

end subroutine establish_halo_cells

!-------------------------------------------------------------------
!Establish and store indices of bins which are on the outer domain

subroutine establish_surface_bins
	use module_set_parameters
	use interfaces, only : error_abort
	implicit none

	integer		:: n
	integer		:: ibin, jbin, kbin

	if (any(nbins(:) .eq. 1)) then
		call error_abort("Error in surface bins -- nbins must be greater than 1")
	else
		nsurfacebins=	2*( nbins(1)   * nbins(2) &
						+  (nbins(3)-2)* nbins(2) &
			        	+  (nbins(3)-2)*(nbins(1)-2))
	endif

	allocate(surfacebins(nsurfacebins,3))

	n = 1
	do kbin=1, nbins(3)+2
	do jbin=1, nbins(2)+2
	do ibin=1, nbins(1)+2

		!Remove inner part of domain
		if((ibin .gt. (2) .and. ibin .lt. (nbins(1)+1)) .and. &
		   (jbin .gt. (2) .and. jbin .lt. (nbins(2)+1)) .and. &
		   (kbin .gt. (2) .and. kbin .lt. (nbins(3)+1))) cycle
		!Remove outer bins leaving only 1 layer of surface bins
		if((ibin .lt. (2) .or. ibin .gt. (nbins(1)+1)) .or. &
		   (jbin .lt. (2) .or. jbin .gt. (nbins(2)+1)) .or. &
		   (kbin .lt. (2) .or. kbin .gt. (nbins(3)+1))) cycle

		surfacebins(n,1)=ibin
		surfacebins(n,2)=jbin
		surfacebins(n,3)=kbin
		n = n + 1

	enddo
	enddo
	enddo

end subroutine establish_surface_bins

!-------------------------------------------------------------------
!Establish and store indices of bins which are in the halo

subroutine establish_halo_bins
	use module_set_parameters
	implicit none

	integer		:: n
	integer		:: ibin, jbin, kbin

	nhalobins  =	2*((nbins(1)+2)*(nbins(2)+2) &
					+  (nbins(3)  )*(nbins(2)+2) &
					+  (nbins(3)  )*(nbins(1)  ))

	allocate(halobins(nhalobins,3))

	n = 1
	do kbin=1, nbins(3)+2
	do jbin=1, nbins(2)+2
	do ibin=1, nbins(1)+2

		!Remove inner part of domain
		if((ibin .gt. (1) .and. ibin .lt. (nbins(1)+2)) .and. &
		   (jbin .gt. (1) .and. jbin .lt. (nbins(2)+2)) .and. &
		   (kbin .gt. (1) .and. kbin .lt. (nbins(3)+2))) cycle

		halobins(n,1)=ibin
		halobins(n,2)=jbin
		halobins(n,3)=kbin
		n = n + 1

	enddo
	enddo
	enddo

end subroutine establish_halo_bins

!-------------------------------------------------------------------
!Establish and store indices of cells which are on the outer domain surface specified 
!by size of buf
! buf = -1 Halo cells
! buf =  0 Domain outer surface cell 
! buf =  N Outer Surface-N

!subroutine establish_surface_cells(buf)
!	use module_set_parameters
!	implicit none

!	integer		:: n
!	integer		:: icell, jcell, kcell
!	integer		:: buf

!	nsurfacecells =	2*((ncells(1)-2*buf)*(ncells(2)-2*buf) &
!					+  (ncells(1)-2-2*buf)*(ncells(3)-2*buf) &
!		        	+  (ncells(2)-2-2*buf)*(ncells(3)-2-2*buf))

!	allocate(surfacecell(nsurfacecells,3))
	
!	n = 1
!	do kcell=1, ncells(1)+2
!	do jcell=1, ncells(2)+2
!	do icell=1, ncells(3)+2
		!Remove inner part of domain
!		if((icell .gt. (2+buf) .and. icell .lt. (ncells(1)+1-buf)) .and. &
!		   (jcell .gt. (2+buf) .and. jcell .lt. (ncells(2)+1-buf)) .and. &
!		   (kcell .gt. (2+buf) .and. kcell .lt. (ncells(3)+1-buf))) cycle
		!Remove outer cells leaving only 1 layer of surface cells
!		if((icell .lt. (2+buf) .or. icell .gt. (ncells(1)+1-buf)) .or. &
!		   (jcell .lt. (2+buf) .or. jcell .gt. (ncells(2)+1-buf)) .or. &
!		   (kcell .lt. (2+buf) .or. kcell .gt. (ncells(3)+1-buf))) cycle

!		surfacecell(n,1)=icell
!		surfacecell(n,2)=jcell
!		surfacecell(n,3)=kcell
!		n = n + 1
!	enddo
!	enddo
!	enddo

!end subroutine establish_surface_cells
