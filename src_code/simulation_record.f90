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
    use boundary_MD, only: bforce_pdf_measure
	!use module_set_parameters, only : velPDF, velPDFMB, velPDF_array

	double precision :: vel


contains

    function get_bin(r) result(bin)
        use computational_constants_MD, only : halfdomain, nhb
        use calculated_properties_MD, only : binsize
        implicit none

        double precision,intent(in),dimension(3) :: r
	    integer,dimension(3) 					 :: bin

        bin = ceiling((r+halfdomain)/binsize)+nhb

    end function get_bin

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
    integer, save :: vmd_skip_count=0

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
        
        vmd_skip_count = vmd_skip_count + 1 
        if (vmd_skip_count .eq. vmd_skip) then

            vmd_skip_count = 0

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

        else
            
            continue
            
        end if

	endif
	
    !Get polymer statistics
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

	!Obtain and record velocity distributions
	if (vPDF_flag .ne. 0) call velocity_PDF_averaging(vPDF_flag)

	!Record boundary force PDFs
	if (bforce_pdf_measure .ne. 0) call bforce_pdf_stats

	!Obtain and record temperature
	if (temperature_outflag .ne. 0)	call temperature_averaging(temperature_outflag)

	!Obtain and record temperature
	if (energy_outflag .ne. 0)	call energy_averaging(energy_outflag)

	!Obtain and record molecular diffusion
	!call evaluate_properties_diffusion

	!Calculate pressure tensor
	if (pressure_outflag .ne. 0) then
		call pressure_averaging(pressure_outflag)
	end if

	if (heatflux_outflag .ne. 0) then
		call heatflux_averaging(heatflux_outflag)
	end if

	call update_simulation_progress_file

end subroutine simulation_record

!==========================================================================
!Calculate kinetic and potential energy as well as temperature and pressure

subroutine evaluate_macroscopic_properties
	use module_record
	use messenger, only : globalise
    use messenger_data_exchange, only : globalSum
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
			write(3000+irank,'(i10,f28.4,6f10.4)'), n , potenergymol(n), r(:,n), globalise(r(:,n))
		enddo
		print*, 'Simulation aborted because max PE has reached an unreasonably high value.'
		print*, 'Inspect fort.(3000+irank) for n, potenergymol, r, r_global'
		call error_abort("STOPPING CODE")
	endif
	totenergy   = kinenergy + potenergy
	temperature = v2sum / real(nd*globalnp,kind(0.d0))
	if (any(periodic.gt.1)) temperature = get_temperature_PUT()
	pressure    = (density/(globalnp*nd))*(v2sum+virial/2.d0) + Pressure_sLRC !N.B. virial/2 as all interactions calculated

end subroutine evaluate_macroscopic_properties

subroutine evaluate_microstate_pressure
	use module_record
    use messenger_data_exchange, only : globalSum
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
				print '(1x,i8,a,f10.3,a,e12.2,a,e10.2,a,f7.3,a,f19.15,a,f19.15,a,f19.15,a,f10.4)', &
				it,';', simtime,';',vsum,';', v2sum,';', temperature,';', &
				kinenergy,';',potenergy,';',totenergy,';',pressure
			case(3:4)
				print '(1x,i7,a,f9.3,a,e9.2,a,e8.1,a,f8.4,a,f8.4,a,f8.4,a,f8.4)', &
				it,';', simtime,';',vsum,';', temperature,';', &
				kinenergy,';',potenergy,';',totenergy,';',pressure
			case default
			end select
		case(1)
			select case(macro_outflag)
			case(1:2)
				print '(1x,i8,a,f10.3,a,e10.3,a,e10.2,a,f7.3,a,f15.11,a,f15.11,a,f15.11,a,f10.4,a,f7.4,a,f9.3)', &
				it,';',simtime,';',vsum,';', v2sum,';', temperature,';', &
				kinenergy,';',potenergy,';',totenergy,';',pressure,';',etevtcf,';',R_g
			case(3:4)
				print '(1x,i7,a,f8.3,a,e8.1,a,e8.1,a,f7.3,a,f7.3,a,f7.3,a,f7.3,a,f6.3,a,f5.1)', &
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
!and calculate Boltzmann H function on a bin by bin basis


subroutine velocity_PDF_averaging(ixyz)
	use module_record
	use librarymod, only : Maxwell_Boltzmann_vel,Maxwell_Boltzmann_speed
    use field_io, only : velocity_PDF_slice_io
	use module_set_parameters, only : velPDF, velPDFMB, velPDF_array
	implicit none

    integer,intent(in)      :: ixyz

	integer                 :: n, i,j,k, jxyz,kxyz
	integer, save		    :: average_count=-1

	integer,dimension(3)	:: cbin
	integer,dimension(:,:),allocatable          :: pdfx,pdfy,pdfz

	double precision 	                        :: Hfunction,vfactor,const
	double precision,save 	                    :: meanstream, meannp
	double precision,dimension(3)	            :: peculiarv,binsize_
	double precision,dimension(:),allocatable 	:: vmagnitude,binloc

	average_count = average_count + 1
	call cumulative_velocity_PDF

	!Write and reset PDF FOR EACH Y SLAB
	if (average_count .eq. Nvpdf_ave) then
        average_count = 0

		select case(ixyz)
		case(1:3)

           	const = sqrt(temperature)
		    !Allocate arrays based on cell 1,1,1 (assuming all identical)
		    allocate(pdfx(nbins(ixyz)+nhb(ixyz),velPDF_array(nhb(1)+1,nhb(2)+1,nhb(3)+1,1)%nbins)); pdfx =0.d0
		    allocate(pdfy(nbins(ixyz)+nhb(ixyz),velPDF_array(nhb(1)+1,nhb(2)+1,nhb(3)+1,2)%nbins)); pdfy =0.d0
		    allocate(pdfz(nbins(ixyz)+nhb(ixyz),velPDF_array(nhb(1)+1,nhb(2)+1,nhb(3)+1,3)%nbins)); pdfz =0.d0
		    allocate(binloc,source=velPDF_array(nhb(1)+1,nhb(2)+1,nhb(3)+1,1)%binvalues())

	        kxyz = mod(ixyz,3)+1
	        jxyz = mod(ixyz+1,3)+1
		    do j = nhb(ixyz)+1,nbins(ixyz)+nhb(ixyz)
        		do i = nhb(jxyz)+1,nbins(jxyz)+nhb(jxyz)
        		do k = nhb(kxyz)+1,nbins(kxyz)+nhb(kxyz)
                    !Collect values per slice for PDF
                    pdfx(j-1,:) = pdfx(j-1,:) + velPDF_array(i,j,k,1)%hist
                    pdfy(j-1,:) = pdfy(j-1,:) + velPDF_array(i,j,k,2)%hist
                    pdfz(j-1,:) = pdfz(j-1,:) + velPDF_array(i,j,k,3)%hist
	            enddo
	            enddo

		    enddo
            call velocity_PDF_slice_io(ixyz,pdfx(:,:),pdfy(:,:),pdfz(:,:))

            deallocate(pdfx,pdfy,pdfz)
            deallocate(binloc)

        	!Write values of bin to file to follow evolution of moments
!            write(14,'(12f12.6)') velPDF_array(5,3 ,5,1)%moments(0),  velPDF_array(5,3 ,5,1)%moments(1),  & 
!                                  velPDF_array(5,3 ,5,1)%moments(2),  velPDF_array(5,3 ,5,1)%moments(3),  &
!                                  velPDF_array(5,10,5,1)%moments(0),  velPDF_array(5,10,5,1)%moments(1),  & 
!                                  velPDF_array(5,10,5,1)%moments(2),  velPDF_array(5,10,5,1)%moments(3),  &
!                                  velPDF_array(5,18,5,1)%moments(0),  velPDF_array(5,18,5,1)%moments(1),  & 
!                                  velPDF_array(5,18,5,1)%moments(2),  velPDF_array(5,18,5,1)%moments(3)

        case(4)

            !OUTPUT A PDF FOR EVERY CELL
!		    do i = 2,nbins(1)+1
!		    do j = 2,nbins(2)+1
!		    do k = 2,nbins(3)+1

!			    !Normalise bins to use for output and H function - so all bins add up to one
!			    allocate(pdf,source=velPDF_array(i,j,k)%normalise())
!			    allocate(binloc,source=velPDF_array(i,j,k)%binvalues())
!			    do n=1,size(pdf,1) 
!				    write(10000+iter/100,'(3i4,2(f10.5))') i,j,k, binloc(n), pdf(n)
!			    enddo

!			    deallocate(pdf)
!			    deallocate(binloc)
!		    enddo
!		    enddo
!		    enddo

        end select

		do j = nhb(2)+1,nbins(2)+nhb(2)
    	do i = nhb(1)+1,nbins(1)+nhb(1)
    	do k = nhb(3)+1,nbins(3)+nhb(3)
        do n = 1,nd
		    velPDF_array(i,j,k,n)%hist = 0
        enddo
        enddo
        enddo
        enddo
		velPDFMB%hist = 0
		meanstream = 0.d0
		meannp = 0.d0


	endif


contains

subroutine cumulative_velocity_PDF
    use module_set_parameters, only : velPDF_array
    implicit none

    integer      :: ixyz

	!Calculate streaming velocity
	binsize_ = domain/nbins

	!Add vmagnitude to PDF histogram for each bin
	allocate(vmagnitude(1))
	do n = 1, np    ! Loop over all particles
		!cbin(:) = ceiling((r(:,n)+0.5d0*domain(:))/binsize_(:))+1
        cbin = get_bin(r(:,n))
        do ixyz = 1,nd
    		vmagnitude(1)=v(ixyz,n)
            call velPDF_array(cbin(1),cbin(2),cbin(3),ixyz)%update(vmagnitude)
        enddo
	enddo
	deallocate(vmagnitude)

end subroutine

end subroutine velocity_PDF_averaging





!!===================================================================================
!!Molecules grouped into velocity ranges to give vel frequency distribution graph
!!and calculate Boltzmann H function
!subroutine evaluate_properties_vdistribution
!	use module_record
!	use librarymod, only : Maxwell_Boltzmann_vel,Maxwell_Boltzmann_speed
!	implicit none

!	integer          :: n, cbin, nvbins
!	double precision :: Hfunction,vfactor,const,streamvel
!	double precision,save :: meanstream, meannp
!	double precision,dimension(3)		:: peculiarv
!	double precision,dimension(:),allocatable :: vmagnitude,normalisedvfd_bin,normalisedvfdMB_bin,binloc

!	!Calculate matrix of velocity magnitudes and bin in histogram
!	allocate(vmagnitude(np))

!	!Calculate streaming velocity
!	streamvel = 0.d0
!	do n = 1, np    ! Loop over all particles
!		streamvel = streamvel + v(1,n)
!	enddo

!	!Add vmagnitude to PDF histogram
!	do n = 1, np    ! Loop over all particles
!		peculiarv = (/ v(1,n)-streamvel/np, v(2,n), v(3,n) /)
!		!vmagnitude(n) = sqrt(dot_product(peculiarv,peculiarv))
!		vmagnitude(n) = v(1,n)
!	enddo
!	call velPDF%update(vmagnitude)

!	!Keep cumulative velocity over ae=veraging period (for analytical comparison)
!	meanstream = meanstream + streamvel
!	meannp     = meannp     + np

!	!Calculate maxwell boltzmann distribution for comparison
!	do n = 1, np    ! Loop over all particles
!		!vmagnitude(n) = Maxwell_Boltzmann_speed(T=temperature,u=streamvel/np)
!		vmagnitude(n) = Maxwell_Boltzmann_vel(T=temperature,u=streamvel/np)
!	enddo
!	call velPDFMB%update(vmagnitude)
!	deallocate(vmagnitude)

!	!Calculate Boltzmann H function using discrete defintion as in
!	!Rapaport p37. N.B. velocity at middle or range is used for vn
!	Hfunction = velPDF%Hfunction()

!	!Write and reset RDF after a number of stored records
!	const = sqrt(temperature)
!	if (mod(iter,1000) .eq. 0) then
!		!Normalise bins to use for output and H function - so all bins add up to one
!		allocate(normalisedvfd_bin,source=velPDF%normalise())
!		allocate(normalisedvfdMB_bin,source=velPDFMB%normalise())
!		allocate(binloc,source=velPDF%binvalues())
!		do n=1,size(normalisedvfd_bin,1) 
!			write(12,'(5(f10.5))') binloc(n), normalisedvfd_bin(n),normalisedvfdMB_bin(n), & 
!									sqrt(2/pi)*((binloc(n)**2)*exp((-binloc(n)**2)/(2*const**2))/(const**3)), & 
!								    (1.d0/(const*sqrt(2.d0*pi)))*exp( -((binloc(n)-meanstream/meannp)**2.d0)/(2.d0*const**2.d0) ) 
!		enddo

!    	!Write values of bin to file to follow evolution of moments
!    	write(14,'(8f17.10)') velPDF%moments(0),  velPDF%moments(1),  velPDF%moments(2),  velPDF%moments(3), & 
!    					      velPDFMB%moments(0),velPDFMB%moments(1),velPDFMB%moments(2),  velPDFMB%moments(3)

!    	!Write values of bin to file to follow evolution of H-function
!    	write(13,'(a,i5, a, f20.10)') 'Boltzmann H function at iteration ', iter , ' is ', Hfunction

!		velPDF%hist = 0
!		velPDFMB%hist = 0
!		meanstream = 0.d0
!		meannp = 0.d0
!	endif

!end subroutine evaluate_properties_vdistribution



!!===================================================================================
!!Molecules grouped into velocity ranges to give vel frequency distribution graph
!!and calculate Boltzmann H function on a bin by bin basis

!subroutine evaluate_properties_vdistribution_perbin
!	use module_record
!	use librarymod, only : Maxwell_Boltzmann_vel,Maxwell_Boltzmann_speed
!	implicit none

!	integer  :: n, i,j,k, nvbins, outinterval!, sumhist
!	integer,dimension(3)	:: cbin
!	double precision 	:: Hfunction,vfactor,const
!	double precision,save 	:: meanstream, meannp
!	double precision,dimension(3)	:: peculiarv,binsize_
!	double precision,dimension(:),allocatable 	:: vmagnitude,normalisedvfd_bin,normalisedvfdMB_bin,binloc
!	double precision,dimension(:,:,:),allocatable 	:: streamvel

!    outinterval = Nvel_ave*tplot

!	!Calculate streaming velocity
!	binsize_ = domain/nbins
!	allocate(streamvel(nbins(1)+2,nbins(2)+2,nbins(3)+2))
!	streamvel = 0.d0
!	do n = 1, np    ! Loop over all particles
!		cbin(:) = ceiling((r(:,n)+0.5d0*domain(:))/binsize_(:))+1
!		streamvel(cbin(1),cbin(2),cbin(3)) = streamvel(cbin(1),cbin(2),cbin(3)) + v(1,n)
!	enddo

!	!Add vmagnitude to PDF histogram for each bin
!	allocate(vmagnitude(1))
!	do n = 1, np    ! Loop over all particles
!		cbin(:) = ceiling((r(:,n)+0.5d0*domain(:))/binsize_(:))+1
!		!peculiarv = (/ v(1,n)-streamvel/np, v(2,n), v(3,n) /)
!		!vmagnitude = sqrt(dot_product(peculiarv,peculiarv))
!		vmagnitude = v(1,n)
!		call velPDF_array(cbin(1),cbin(2),cbin(3))%update(vmagnitude)
!	enddo
!	deallocate(vmagnitude)

!	!Keep cumulative velocity over averaging period (for analytical comparison)
!	meanstream = meanstream + sum(streamvel)
!	meannp     = meannp     + np

!	!Calculate maxwell boltzmann distribution for comparison
!	allocate(vmagnitude(np))
!	do n = 1, np    ! Loop over all particles
!		!vmagnitude(n) = Maxwell_Boltzmann_speed(T=temperature,u=streamvel/np)
!		vmagnitude(n) = Maxwell_Boltzmann_vel(T=temperature,u=meanstream/meannp)
!	enddo
!	call velPDFMB%update(vmagnitude)

!	!Calculate Boltzmann H function using discrete defintion as in
!	!Rapaport p37. N.B. velocity at middle or range is used for vn
!	!Hfunction = velPDF%Hfunction()

!	!Write and reset PDF FOR EACH Y SLAB
!	const = sqrt(temperature)

!	if (mod(iter,outinterval) .eq. 0) then

!        !print*, "WARNING -- WRITING OUT PDF INFORMATION"
!        

!		!Allocate arrays based on cell 1,1,1 (assuming all identical)
!		allocate(normalisedvfd_bin(size(velPDF_array(2,2,2)%normalise(),1))); normalisedvfd_bin =0.d0
!		allocate(binloc,source=velPDF_array(2,2,2)%binvalues())
!		allocate(normalisedvfdMB_bin,source=velPDFMB%normalise())
!        !sumhist = 0
!		do j = 2,nbins(2)+1
!    		do i = 2,nbins(1)+1
!    		do k = 2,nbins(3)+1
!                !Collect values per slice for PDF
!                normalisedvfd_bin(:) = normalisedvfd_bin(:) + velPDF_array(i,j,k)%normalise()
!                !sumhist = sumhist + sum(velPDF_array(i,j,k)%hist)
!	        enddo
!	        enddo
!            normalisedvfd_bin = normalisedvfd_bin/(nbins(1)*nbins(3))
!            !Write values per slice to file
!	        do n=1,size(normalisedvfd_bin,1) 
!		        write(10000+iter/outinterval,'(i5,5(f10.5))') j, binloc(n), normalisedvfd_bin(n),normalisedvfdMB_bin(n), & 
!								        sqrt(2/pi)*((binloc(n)**2)*exp((-binloc(n)**2)/(2*const**2))/(const**3)), & 
!								        (1.d0/(const*sqrt(2.d0*pi)))*exp( -((binloc(n)-meanstream/meannp)**2.d0)/(2.d0*const**2.d0) ) 
!	        enddo
!		enddo
!        !print*, 'SUM of HIST', sumhist,normalisedvfd_bin
!        deallocate(normalisedvfd_bin)
!        deallocate(normalisedvfdMB_bin)
!        deallocate(binloc)

!		do j = 2,nbins(2)+1
!    	do i = 2,nbins(1)+1
!    	do k = 2,nbins(3)+1
!		    velPDF_array(i,j,k)%hist = 0
!        enddo
!        enddo
!        enddo
!		velPDFMB%hist = 0
!		meanstream = 0.d0
!		meannp = 0.d0

!        !OUTPUT A PDF FOR EVERY CELL
!!		do i = 2,nbins(1)+1
!!		do j = 2,nbins(2)+1
!!		do k = 2,nbins(3)+1

!!			!Normalise bins to use for output and H function - so all bins add up to one
!!			allocate(normalisedvfd_bin,source=velPDF_array(i,j,k)%normalise())
!!			allocate(normalisedvfdMB_bin,source=velPDFMB%normalise())
!!			allocate(binloc,source=velPDF_array(i,j,k)%binvalues())
!!			do n=1,size(normalisedvfd_bin,1) 
!!				write(10000+iter/100,'(3i4,5(f10.5))') i,j,k, binloc(n), normalisedvfd_bin(n),normalisedvfdMB_bin(n), & 
!!										sqrt(2/pi)*((binloc(n)**2)*exp((-binloc(n)**2)/(2*const**2))/(const**3)), & 
!!										(1.d0/(const*sqrt(2.d0*pi)))*exp( -((binloc(n)-meanstream/meannp)**2.d0)/(2.d0*const**2.d0) ) 
!!			enddo

!!!			!Write values of bin to file to follow evolution of moments
!!!			write(20000+i+(j-1)*nbins(1)+(k-1)*nbins(1)*nbins(2),'(8f17.10)') velPDF%moments(0),  velPDF%moments(1),  velPDF%moments(2),  velPDF%moments(3), & 
!!!							      velPDFMB%moments(0),velPDFMB%moments(1),velPDFMB%moments(2),  velPDFMB%moments(3)

!!!			!Write values of bin to file to follow evolution of H-function
!!!			write(30000+i+(j-1)*nbins(1)+(k-1)*nbins(1)*nbins(2),'(a,i5, a, f20.10)') 'Boltzmann H function at iteration ', iter , ' is ', Hfunction

!!			deallocate(normalisedvfd_bin)
!!			deallocate(normalisedvfdMB_bin)
!!			deallocate(binloc)
!!		enddo
!!		enddo
!!		enddo

!!		velPDF%hist = 0
!!		velPDFMB%hist = 0
!!		meanstream = 0.d0
!!		meannp = 0.d0
!	
!		!stop "Stop after one write of evaluate_properties_vdistribution_perbin"
!	endif

!end subroutine evaluate_properties_vdistribution_perbin


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
	use module_compute_forces, only : compute_force_and_potential_at
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
    use messenger_data_exchange, only : globalSum
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
                if (j .eq. 0) then
                    call missing_bond_error(i,nbond)
                end if
				j_sub  = monomer(j)%subchainID
				if (j_sub.lt.i_sub) cycle  !Avoid counting backwards
				rij(:) = r(:,j) - r(:,i)
                rij(:)        = rij(:) - &
                                globaldomain(:)*anint(rij(:)/globaldomain(:))
				etev_0(chain,:) = etev_0(chain,:) + rij(:)
			end do
		end do
		call globalSum(etev_0,nchains,nd)

	end if

	!--------------------------------------------------------------------
	!Calculate all end-to-end vectors for file output
	etev = 0.d0
	do i=1,np
		chain = monomer(i)%chainID
		i_sub = monomer(i)%subchainID
		funcy = monomer(i)%funcy
		do nbond=1,funcy
			j             = bond(nbond,i)          ! Molecule number j is nth bond to i
			j_sub         = monomer(j)%subchainID  ! Find subchain ID of mol j
			if (j_sub.lt.i_sub) cycle              ! Avoid counting backwards
			rij(:)        = r(:,j) - r(:,i)                     
			rij(:)        = rij(:) - &
                            globaldomain(:)*anint(rij(:)/globaldomain(:))
			etev(chain,:) = etev(chain,:) + rij(:)
		end do
	end do
	call globalSum(etev,nchains,nd)
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
    use messenger_data_exchange, only : globalSum
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

	call globalSum(r_cm,nchains,nd)

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

            select case(split_pol_sol_stats)
            case(0)
                call mass_bin_io(volume_mass,'bins')
                volume_mass = 0
            case(1)
                call mass_bin_io(volume_mass_s,'solv')
                call mass_bin_io(volume_mass_p,'poly')
                call mass_bin_io(volume_mass,'bins')
                volume_mass_s = 0
                volume_mass_p = 0
                volume_mass = 0
            case default
                call error_abort('Invalid split_pol_sol_stats in m averaging')
            end select
                
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
	use librarymod, only : cpolariser
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
        select case (split_pol_sol_stats)
        case(0)
            do n = 1,np
                !Add up current volume mass densities
                !ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb
                ibin(:) = get_bin(r(:,n))
                volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + 1
            enddo
        case(1)
            do n = 1,np
                !Add up current volume mass densities
                !ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb
                ibin = get_bin(r(:,n))
                if (monomer(n)%chainID .eq. 0) then
                    !Skip wall molecules
                    if (any(tag(n).eq.tether_tags)) cycle
                    volume_mass_s(ibin(1),ibin(2),ibin(3)) = &
                    volume_mass_s(ibin(1),ibin(2),ibin(3)) + 1
                else
                    volume_mass_p(ibin(1),ibin(2),ibin(3)) = &
                    volume_mass_p(ibin(1),ibin(2),ibin(3)) + 1
                end if
            enddo
            volume_mass = volume_mass_s + volume_mass_p
        case default 
            call error_abort('Mass output not implemented for this potential flag')
        end select 

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
	double precision,dimension(3) 	:: Vbinsize

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

		do n=1,np
			!Save streaming velocity per molecule
			!ib(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
            ib(:) = get_bin(r(:,n))
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

                select case(split_pol_sol_stats)
                case(0)
                    call velocity_bin_io(volume_mass, volume_momentum,'bins')
                    volume_mass = 0
                    volume_momentum = 0.d0
                case(1)
                    call velocity_bin_io(volume_mass_s, volume_momentum_s, 'solv')
                    call velocity_bin_io(volume_mass_p, volume_momentum_p, 'poly')
                    call velocity_bin_io(volume_mass, volume_momentum, 'bins')
                    volume_mass_s = 0
                    volume_mass_p = 0
                    volume_mass = 0
                    volume_momentum_s = 0.d0
                    volume_momentum_p = 0.d0
                    volume_momentum = 0.d0 
                case default
                    call error_abort('Invalid split_pol_sol_stats in m averaging')
                end select

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

!-----------------------------------------------------------------------------------
!Add velocities to running total, with 2D slice in ixyz = 1,2,3 or
!in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine cumulative_velocity(ixyz)
	use module_record
	use concentric_cylinders
	use messenger, only: globalise, localise
	use linked_list
    use librarymod, only : cpolariser, cpolarisev
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

        Vbinsize(:) = domain(:) / nbins(:) 
        select case (split_pol_sol_stats)
        case(0)
            !Reset Control Volume momentum 
            do n = 1,np
                !Add up current volume mass and momentum densities
                !ibin(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
                ibin(:) = get_bin(r(:,n))
                volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + 1
                volume_momentum(ibin(1),ibin(2),ibin(3),:) = volume_momentum(ibin(1),ibin(2),ibin(3),:) & 
                                                            + v(:,n) + slidev(:,n)
            enddo
        case(1)
            !Reset Control Volume momentum 
            do n = 1,np
                !Add up current volume mass and momentum densities
                !ibin(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
                ibin(:) = get_bin(r(:,n))
                if (monomer(n)%chainID .eq. 0) then
                    if (any(tag(n).eq.tether_tags)) cycle
                    volume_mass_s(ibin(1),ibin(2),ibin(3)) = volume_mass_s(ibin(1),ibin(2),ibin(3)) + 1
                    volume_momentum_s(ibin(1),ibin(2),ibin(3),:) = volume_momentum_s(ibin(1),ibin(2),ibin(3),:) & 
                                                                + v(:,n) + slidev(:,n)
                else
                    volume_mass_p(ibin(1),ibin(2),ibin(3)) = volume_mass_p(ibin(1),ibin(2),ibin(3)) + 1
                    volume_momentum_p(ibin(1),ibin(2),ibin(3),:) = volume_momentum_p(ibin(1),ibin(2),ibin(3),:) & 
                                                                + v(:,n) + slidev(:,n)
                end if
            enddo
            volume_mass = volume_mass_s + volume_mass_p
            volume_momentum = volume_momentum_s + volume_momentum_p
        case default 
            call error_abort('Momentum output not implemented for this potential flag')
        end select 

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
    use librarymod, only : cpolariser
	implicit none

	integer							:: n,ixyz
	integer         				:: cbin
	integer		,dimension(3)		:: ibin
	double precision				:: slicebinsize
	double precision,dimension(3) 	:: Tbinsize 

	integer :: br, bt, bz
	real(kind(0.d0)) :: fluiddomain_cyl(3), rglob(3), rpol(3)

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
			!ibin(:) = ceiling((r(:,n)+halfdomain(:))/Tbinsize(:)) + nhb
            ibin(:) = get_bin(r(:,n))
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

			if (velocity_outflag .ne. 5) then
                cyl_mass(br,bt,bz) = cyl_mass(br,bt,bz) + 1
            end if

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
!		RECORD ENERGY AT LOCATION IN SPACE
! Either by binning or taking slices on a plane, the energy field in the molecular
! system is recorded and output
!===================================================================================

!Calculate averaged energy components of each bin or slice with 2D slice in 
!ixyz = 1,2,3 or in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine energy_averaging(ixyz)
	use module_record
	use field_io, only : energy_bin_io
	use linked_list
	use concentric_cylinders
	implicit none

	integer				:: ixyz
	integer, save		:: average_count=-1
	
	average_count = average_count + 1
	call cumulative_energy(ixyz)
    !print*, 'Total energy = ',average_count, sum(volume_energy)/(average_count*np), 0.5*sum(potenergymol(1:np)/np)
	if (average_count .eq. Nenergy_ave) then
		average_count = 0

		select case(ixyz)
		case(1:3)
            stop "No slices for energy binning"
		case(4)
            call energy_bin_io(volume_energy,'bins')
			!Reset energy bins
			volume_energy = 0.d0
		case(5)
		    stop "No_cpol for energy binning"
		case default
			stop "Error input for energy averaging incorrect"
		end select

	endif

end subroutine energy_averaging

!-----------------------------------------------------------------------------------
!Add energy to running total, with 2D slice in ixyz = 1,2,3 or
!in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine cumulative_energy(ixyz)
	use module_record
	use concentric_cylinders
	use messenger, only: globalise, localise
	use linked_list
    use librarymod, only : cpolariser
	implicit none

	integer							:: n,ixyz
	integer         				:: cbin
	integer		,dimension(3)		:: ibin
	double precision				:: slicebinsize, energy
	double precision,dimension(3) 	:: Tbinsize, velvect

	select case(ixyz)
	!energy measurement is a number of 2D slices through the domain
	case(1:3)
        stop "No slices for energy binning"
	!energy measurement for 3D bins throughout the domain
	case(4)
 
		do n = 1,np
			!Add up current volume mass and energy densities
            ibin(:) = get_bin(r(:,n))
			if (peculiar_flag .eq. 0) then
		        velvect(:) = v(:,n) + 0.5d0*a(:,n)*delta_t! + slidev(:,n)
		        energy = 0.5d0*(dot_product(velvect,velvect)+potenergymol(n))

				volume_energy(ibin(1),ibin(2),ibin(3)) = & 
                    volume_energy(ibin(1),ibin(2),ibin(3)) + energy
            else
                stop "No peculiar removal in energy - do it in post processing"
			endif

		enddo

	case(5)
	    stop "No_cpol for energy binning"
	case default 
		stop "energy Binning Error"
	end select
	 
end subroutine cumulative_energy


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
    use messenger_data_exchange, only : globalSum
	implicit none

	integer								:: sample_count,n,ixyz,jxyz,kxyz
	double precision, dimension(3)		:: velvect
	double precision, dimension(3)      :: rglob
	double precision, dimension(3,3)	:: Pxytemp


	Pxytemp = 0.d0

	select case(ixyz)
	case(1)

		Pxymol = 0.d0

		!Factor of 2 as every interaction calculated
		rfmol = rfmol / 2.d0

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
		call globalSum(Pxytemp(:,1), nd)	
		call globalSum(Pxytemp(:,2), nd) 	
		call globalSum(Pxytemp(:,3), nd)

		!Divide sum of stress by volume -- Should this include sample_count - divide by Nstress above??
		Pxytemp = Pxytemp / volume

		!Add current pressure to cumulative average
		Pxy = Pxy + Pxytemp

		!Reset position force tensor before calculation
		rfmol = 0.d0
	case(2)
		!VOLUME AVERAGE STRESS CALCULATION
		!Calculate Position (x) force for configurational part of stress tensor
		call  simulation_compute_rfbins!(1,nbins(1)+2*nhb(1), & 
                                       ! 1,nbins(2)+2*nhb(2), & 
                                       ! 1,nbins(3)+2*nhb(3))
		!Calculate velocity (x) velocity for kinetic part of stress tensor
		if (nbins(1) .eq. ncells(1) .and. & 
		    nbins(2) .eq. ncells(2) .and. & 
		    nbins(3) .eq. ncells(3)) then
			call  simulation_compute_kinetic_VA_cells(2,ncells(1)+1,2,ncells(2)+1,2,ncells(3)+1)
		else
			call  simulation_compute_kinetic_VA(1,nbins(1), & 
                                                1,nbins(2), & 
                                                1,nbins(3))
		endif

		!Add results to cumulative total
		Pxybin(:,:,:,:,:) =     vvbin(:,:,:,:,:) 				& 
						      + rfbin(  1+nhb(1):nbins(1)+nhb(1), 			& 
							      		1+nhb(2):nbins(2)+nhb(2), 			& 
							      		1+nhb(3):nbins(3)+nhb(3),:,:)/2.d0

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
		call simulation_compute_kinetic_VA_cpol()
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
	integer 							:: bin(3)
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

subroutine simulation_compute_kinetic_VA_cpol()
	use module_record
	use physical_constants_MD
	use concentric_cylinders
	use messenger, only: globalise
	use librarymod, only: cpolariser, cpolariseT, outerprod
	implicit none

	integer :: n
	integer :: br, bt, bz
	real(kind(0.d0)) :: VAbinsize(3), velvect(3)
	real(kind(0.d0)) :: ripol(3), vvpol(3,3)

	!Determine bin sizes
	VAbinsize(1) = (r_io - r_oi) / cpol_bins(1)
	VAbinsize(2) = 2.d0*pi       / cpol_bins(2)
	VAbinsize(3) = domain(3)     / cpol_bins(3)

	! Add kinetic part of pressure tensor for all molecules
	do n = 1, np

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
        if (tag(n) .eq. cyl_teth_thermo_rotate) cycle
	
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
		oldi  =>  cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list
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

			currenti  =>  oldi
			oldi  =>  currenti%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo
	enddo
	enddo

	nullify(oldi)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required

end subroutine simulation_compute_kinetic_VA_cells




!---------------------------------------------------------------
! Contains all routines for Volume Average Stress Calculation

module Volume_average_pressure


contains

!---------------------------------------------------------------
! Top level for stress calculation
! Author: Edward Smith (edward.smith05@imperial.ac.uk)
!Inputs: 
!       ri       - Position of molecules i 
!       rj       - Position of molecules j
!       accijmag - Magnitude of force
!       domain   - Total domain size
!       nbins    - Number of averaging bins in domain
!       nhb      - Number of halos bins (for periodic boundaries/parallel simulations)
!  VA_calcmethod - Method of VA calculation
!       rfbin    - Returned array with length of interaction
!     
! N.B. molecular positon from -halfdomain to +halfdomain       


subroutine pressure_tensor_forces_VA(ri, rj, rF, domain,  & 
                                     nbins, nhb, rfbin,   & 
                                     VA_calcmethod,       &
                                     VA_line_samples)
    use computational_constants_MD, only : cellsidelength
    use module_record, only : get_bin

    integer,intent(in)                                      :: VA_calcmethod
    integer,intent(in),optional                             :: VA_line_samples
    integer,dimension(3),intent(in)                         :: nbins,nhb
	real(kind(0.d0)), dimension(:,:), intent(in)            :: rF
	double precision,dimension(3), intent(in)		        :: ri, rj, domain
	double precision,dimension(:,:,:,:,:), intent(inout)	:: rfbin

    integer                                                 :: VA_line_samples_
	double precision,dimension(3)                           :: halfdomain

    halfdomain = 0.5d0 * domain

	!Select requested configurational line partition methodology
	select case(VA_calcmethod)
	case(0)
		call pressure_tensor_forces_H(ri,rj,rF)
	case(1)
        if (present(VA_line_samples)) then
            VA_line_samples_ = VA_line_samples
        else
            VA_line_samples_ = 20
        endif
		call pressure_tensor_forces_VA_trap(ri,rj,rF,VA_line_samples_)
	case(2)
		call pressure_tensor_forces_VA_exact(ri,rj,rF)
	case default
		stop "Error - VA_calcmethod incorrect"
	end select 

contains

!====================================================================================
! CONFIGURATIONAL INTERACTIONS ASSIGNED TO CONTAINING BIN - HALF PER MOLECULE 
! Author: Edward Smith (edward.smith05@imperial.ac.uk)
! Use a configurational expression, similar to the virial with interaction partitioned 
! between bins with HALF given per bin and nothing in intermediate bins
! This is also called the Harasima contour in NAMD (MD package) and the paper it references
! Jacob Sonne,a Flemming Y. Hansen, and Günther H. Peters J.CHEM.PHYS. 122, 124903 (2005)
! The original paper by A. Harasima, Advances in Chemical Physics, Volume 1 P201 is unavailable

subroutine pressure_tensor_forces_H(ri,rj,rF)
	implicit none

	real(kind(0.d0)), dimension(:,:), intent(in)    :: rF
	double precision,dimension(3), intent(in)		:: ri, rj

	integer											:: ixyz, jxyz
	integer											:: diff
	integer,dimension(3)							:: ibin, jbin, bindiff 
	double precision,dimension(3)					:: VAbinsize, rij

    rij = rj - ri

	!================================================================
	!= Establish bins for molecules & check number of required bins	=
	!================================================================

    ibin(:) = get_bin(ri)
    jbin(:) = get_bin(rj)
    bindiff(:) = abs(ibin(:) - jbin(:)) + 1

	!================================================================
	!=			Assign to bins				=
	!================================================================

	!Check difference between bins i and j
	diff = bindiff(1)+bindiff(2)+bindiff(3)

	!Ignore values outside of domain resulting from shifted bin 
	if (minval(ibin) .lt. 1) diff = 0
	if (minval(jbin) .lt. 1) diff = 0
	if (maxval(ibin) .gt. maxval(nbins+2*nhb)) diff = 0
	if (maxval(jbin) .gt. maxval(nbins+2*nhb)) diff = 0

	select case(diff)
	!================Skip Force addition===========================
	case(0)
	!Do Nothing

	!================Molecules in same bin===========================
	case(3)

		!------Add molecules to bin-----

		!Factor of two as molecule i and molecule j are both counted in bin i
		rfbin(ibin(1),ibin(2),ibin(3),:,:) = & 
            rfbin(ibin(1),ibin(2),ibin(3),:,:) + rF(:,:)


	!===================Interactions split over 2 or more cells ==========
	case(4:)

		!------Add molecules to bin-----
		!Molecule i and j contribution split between bins

		!-----------Molecule i bin-----------
		rfbin(ibin(1),ibin(2),ibin(3),:,:) = & 
            rfbin(ibin(1),ibin(2),ibin(3),:,:) + 0.5d0*rF(:,:)

		!-----------Molecule j bin-----------
		rfbin(jbin(1),jbin(2),jbin(3),:,:) = & 
            rfbin(jbin(1),jbin(2),jbin(3),:,:) + 0.5d0*rF(:,:)
		
	case default

        print*, 'bindiff = ', ibin, jbin, bindiff
		stop "Harasima Stress AVERAGING ERROR"

	end select

end subroutine pressure_tensor_forces_H


!====================================================================================
! VOLUME AVERAGE CONFIGURATIONAL EXPRESSION
! Author: David Trevelyan, July 2013 (david.trevelyan06@imperial.ac.uk)
! Linear trajectory path sampled to find approximate values of l_ij 
!(less accurate, but much easier to understand)

subroutine pressure_tensor_forces_VA_trap(ri,rj,rF,VA_line_samples)
	use librarymod, only: outerprod
	implicit none

    integer,intent(in)                           :: VA_line_samples
	real(kind(0.d0)), dimension(:,:), intent(in) :: rF
	real(kind(0.d0)), dimension(3),   intent(in) :: ri, rj

	integer :: ss
	integer :: Ns, bin(3)
	real(kind(0.d0)) :: s, ds, rs(3), rij(3), Fij(3), VAbinsize(3)

	VAbinsize(:) = domain(:) / nbins(:)
	rij = rj - ri
	!rF = outerprod(rij, accijmag*rij)

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

			!bin(:) = ceiling((rs(:)+halfdomain(:))/VAbinsize(:)) + 1
            bin(:) = get_bin(rs(:))
			rfbin(bin(1),bin(2),bin(3),:,:) =  &
			rfbin(bin(1),bin(2),bin(3),:,:) + rF(:,:)/real(Ns,kind(0.d0))

		end if

		s = s + ds	

	end do	
	
end subroutine pressure_tensor_forces_VA_trap

!===============================================================================
! VOLUME AVERAGE CONFIGURATIONAL EXPRESSION
! Author: Edward Smith (edward.smith05@imperial.ac.uk)
! Use a configurational expression, similar to the virial with interaction partitioned 
! between bins based on the share of the interaction between two molecules in a given bin
! N.B. Assume no more than 1 bin between bins containing the molecules
! This is the Volume Average stress of Lutsko, although this is also called the
! Irving Kirkwood contour in NAMD (MD package) and the paper it references
! Jacob Sonne,a Flemming Y. Hansen, and Günther H. Peters J.CHEM.PHYS. 122, 124903 (2005)
! 								╭∩╮（︶︿︶）╭∩╮﻿

subroutine pressure_tensor_forces_VA_exact(ri,rj,rF)
	use librarymod, only : magnitude, plane_line_intersect
	implicit none

	real(kind(0.d0)), dimension(:,:), intent(in)    :: rF
	double precision,dimension(3), intent(in)		:: ri, rj


	integer											:: ixyz, i,j,k,l,n
	integer											:: diff
	integer,dimension(3)							:: ibin, jbin, bindiff 
	integer,dimension(:,:)		   ,allocatable		:: interbin, interbindiff
	double precision,dimension(3)					:: VAbinsize, normal, p,temp1,temp2,temp3
	double precision,dimension(3)                   :: rij
	double precision,dimension(:,:)	   ,allocatable	:: intersection
	double precision,dimension(:,:,:)  ,allocatable	:: MLfrac !Magnitude of fraction of stress
	double precision,dimension(:,:,:,:),allocatable	:: Lfrac  !Fraction of stress in a given cell

    !Define rij
    rij = rj - ri

	!================================================================
	!= Establish bins for molecules & check number of required bins	=
	!================================================================

	do ixyz = 1,3

		VAbinsize(ixyz) = domain(ixyz) / nbins(ixyz)
		if (VAbinsize(ixyz) .lt. cellsidelength(ixyz)) &
                     stop "Binsize bigger than cellsize ~ Not ok for volume averaging"

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
	if (maxval(ibin) .gt. maxval(nbins+2*nhb)) diff = 0
	if (maxval(jbin) .gt. maxval(nbins+2*nhb)) diff = 0

	select case(diff)
	!================Skip Force addition===========================
	case(0)
	!Do Nothing

	!================Molecules in same bin===========================
	case(3)

		!------Add molecules to bin-----

		!Factor of two as molecule i and molecule j are both counted in bin i
		rfbin(ibin(1),ibin(2),ibin(3),:,:) =& 
            rfbin(ibin(1),ibin(2),ibin(3),:,:) + rF(:,:)

	!===================Interactions split over 2 cells only==========
	case(4)

		!Allocate array for vectors from molecules to bin walls
		allocate(Lfrac(bindiff(1),bindiff(2),bindiff(3),3))
		!Allocate arrays for magnitude of vector Lfrac
		allocate(MLfrac(bindiff(1),bindiff(2),bindiff(3)))
		!Allocate array for one bin intersection point
		allocate(intersection(3,1))
		Lfrac = 0.d0
		MLfrac = 0.d0

		do ixyz = 1,3
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

		!-----------Molecule i bin-----------
		rfbin(ibin(1),ibin(2),ibin(3),:,:) = & 
            rfbin(ibin(1),ibin(2),ibin(3),:,:) &
			  + rF(:,:)*MLfrac(1,1,1)

		!-----------Molecule j bin-----------
		rfbin(jbin(1),jbin(2),jbin(3),:,:) = & 
            rfbin(jbin(1),jbin(2),jbin(3),:,:) &
			  + rF(:,:)*MLfrac(bindiff(1),bindiff(2),bindiff(3))

		deallocate(intersection)
		deallocate(Lfrac)
		deallocate(MLfrac)
		
	!==============Interactions split over intermediate cell===========
	case(5)

		!Allocate array for vectors from molecules to bin walls
		allocate(Lfrac(bindiff(1),bindiff(2),bindiff(3),3))
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

		do ixyz = 1,3

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

		!-----------Molecule i bin-----------
		rfbin(ibin(1),ibin(2),ibin(3),:,:) = & 
            rfbin(ibin(1),ibin(2),ibin(3),:,:) + rF(:,:)*MLfrac(1,1,1)

		!-----------Intermediate Bin-----------
		!If both bins are in the domain then contribution is added for both molecules
		rfbin(interbin(1,1),interbin(2,1),interbin(3,1),:,:) = &
			rfbin(interbin(1,1),interbin(2,1),interbin(3,1),:,:) &
			+ rF(:,:)*MLfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1))

		!-----------Molecule j bin-----------
		rfbin(jbin(1),jbin(2),jbin(3),:,:) = & 
            rfbin(jbin(1),jbin(2),jbin(3),:,:) &
			+ rF(:,:)*MLfrac(bindiff(1),bindiff(2),bindiff(3))

		deallocate(Lfrac)
		deallocate(MLfrac)
		deallocate(interbin)
		deallocate(intersection)
		deallocate(interbindiff)

	!===========Interactions split over 2 intermediate cell===============
	case (6)

		!Allocate array for vectors from molecules to bin walls
		allocate(Lfrac(bindiff(1),bindiff(2),bindiff(3),3))
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

		do ixyz = 1,3

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
		if (magnitude(temp1).gt.magnitude(temp2)) then
			if (magnitude(temp1).gt.magnitude(temp3)) then
				k = 1
				l = 2
			else
				k = 1
				l = 3
			endif
		else
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
				i = 1
			else
				i = 3
			endif
		else
			if (magnitude(temp2) .lt. magnitude(temp3)) then
				i = 2
			else
				i = 3
			endif
		endif

		!Check which plane is closest to j 
		temp1 = rj(:)-intersection(:,1)
		temp2 = rj(:)-intersection(:,2)
		temp3 = rj(:)-intersection(:,3)
		if (magnitude(temp1).lt.magnitude(temp2))then
			if (magnitude(temp1).lt.magnitude(temp3)) then
				j = 1
			else
				j = 3
			endif
		else
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
		Lfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1),:) = & 
                intersection(:,l)-intersection(:,(6-k-l))
		Lfrac(interbindiff(1,2),interbindiff(2,2),interbindiff(3,2),:) = & 
                intersection(:,k)-intersection(:,(6-k-l))

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

		!-----------Molecule i bin-----------
		rfbin(ibin(1),ibin(2),ibin(3),:,:) = & 
            rfbin(ibin(1),ibin(2),ibin(3),:,:) + rF(:,:)*MLfrac(1,1,1)

		!-----------1st Intermediate Bin-----------
		!Intermediate and i bin in domain, j in halo - add intermediate for molecule i
		rfbin(interbin(1,1),interbin(2,1),interbin(3,1),:,:) = &
			rfbin(interbin(1,1),interbin(2,1),interbin(3,1),:,:) &
			+ rF(:,:)*MLfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1))

		!-----------2nd Intermediate Bin-----------
		!Intermediate and i bin in domain, j in halo - add intermediate for molecule i
		rfbin(interbin(1,2),interbin(2,2),interbin(3,2),:,:) = &
			rfbin(interbin(1,2),interbin(2,2),interbin(3,2),:,:) &
			+ rF(:,:)*MLfrac(interbindiff(1,2),interbindiff(2,2),interbindiff(3,2))

		!-----------Molecule j bin-----------
		rfbin(jbin(1),jbin(2),jbin(3),:,:) = & 
            rfbin(jbin(1),jbin(2),jbin(3),:,:) &
			+ rF(:,:)*MLfrac(bindiff(1),bindiff(2),bindiff(3))

		deallocate(intersection)
		deallocate(Lfrac)
		deallocate(MLfrac)
		deallocate(interbin)
		deallocate(interbindiff)
		
	case default

	    stop "VOLUME AVERAGING ERROR"

	end select

end subroutine pressure_tensor_forces_VA_exact

end subroutine pressure_tensor_forces_VA

end module Volume_average_pressure



!========================================================================
!Compute Volume Averaged stress using all cells including halos

subroutine simulation_compute_rfbins!(imin, imax, jmin, jmax, kmin, kmax)
    use Volume_average_pressure
	use module_compute_forces
	use librarymod, only: outerprod
	implicit none

	!integer,intent(in)  			:: imin, jmin, kmin, imax, jmax, kmax

	integer                         :: i, j, ixyz !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp 
	integer							:: molnoi, molnoj
	integer							:: icellmin,jcellmin,kcellmin,icellmax,jcellmax,kcellmax
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	double precision,dimension(3)	:: vi_t, cellsperbin
	real(kind(0.d0)), dimension(3,3):: rF
	real(kind(0.d0)), dimension(3,1):: rFv
	!rfbin = 0.d0
	!allocate(rijsum(nd,np+extralloc)) !Sum of rij for each i, used for SLLOD algorithm
	!rijsum = 0.d0

	!Calculate bin to cell ratio
	cellsperbin = 1.d0/binspercell !ceiling(ncells(1)/dble(nbins(1)))

    ! Still need to loop over every cell (i.e. get all interactions) if
    ! bins are bigger than cells
	where (cellsperbin .ge. 1.d0) cellsperbin = 1.d0

	!Get cell number from bin numbers
!	icellmin = (imin-1)*cellsperbin(1)+1+(1-cellsperbin(1))
!	icellmax =  imax   *cellsperbin(1)  +(1-cellsperbin(1))
!	jcellmin = (jmin-1)*cellsperbin(2)+1+(1-cellsperbin(2))
!	jcellmax =  jmax   *cellsperbin(2)  +(1-cellsperbin(2))
!	kcellmin = (kmin-1)*cellsperbin(3)+1+(1-cellsperbin(3))
!	kcellmax =  kmax   *cellsperbin(3)  +(1-cellsperbin(3))

    icellmin = 1; icellmax = ncells(1) + 2
    jcellmin = 1; jcellmax = ncells(2) + 2
    kcellmin = 1; kcellmax = ncells(3) + 2

	do kcell=kcellmin, kcellmax
	do jcell=jcellmin, jcellmax 
	do icell=icellmin, icellmax 

	
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		do i = 1,cellnp					!Step through each particle in list 
			molnoi = oldi%molno 	 	!Number of molecule
			ri = r(:,molnoi)         	!Retrieve ri

			do kcellshift = -1,1
			do jcellshift = -1,1
			do icellshift = -1,1

				!Prevents out of range values in i
				if (icell+icellshift .lt. icellmin) cycle
				if (icell+icellshift .gt. icellmax) cycle
				!Prevents out of range values in j
				if (jcell+jcellshift .lt. jcellmin) cycle
				if (jcell+jcellshift .gt. jcellmax) cycle
				!Prevents out of range values in k
				if (kcell+kcellshift .lt. kcellmin) cycle
				if (kcell+kcellshift .gt. kcellmax) cycle

				oldj => cell%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
				adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

				do j = 1,adjacentcellnp          !Step through all j for each i

					molnoj = oldj%molno 	 !Number of molecule
					rj = r(:,molnoj)         !Retrieve rj

					currentj => oldj
					oldj => currentj%next    !Use pointer in datatype to obtain next item in list

					if(molnoi==molnoj) cycle !Check to prevent interaction with self

					rij2=0                   !Set rij^2 to zero
					rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j
					rij2 = dot_product(rij,rij)	!Square of vector calculated

					if (rij2 < rcutoff2) then

                        !---------------------------------------
                        ! - Get volume average pressure tensor -

				        if (pressure_outflag .eq. 2) then
						    !Linear magnitude of acceleration for each molecule
						    invrij2 = 1.d0/rij2                 !Invert value
						    accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)
                            rf = outerprod(rij, accijmag*rij)

						    !Select requested configurational line partition methodology
                            call pressure_tensor_forces_VA(ri, rj, rF, domain,  & 
                                                           nbins, nhb,rfbin,    & 
                                                           VA_calcmethod,       &
                                                           VA_line_samples)
                        endif
                        !----------------------------------------
                        ! - Get volume average pressure heating -
				        if (heatflux_outflag .eq. 2) then
                            !Get the velocity, v, at time t 
                            ! ( This is the reason we need to do this after
                            !   the force calculation so we know a(t) ) 
                            vi_t(:) = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)

						    !Select requested configurational line partition methodology
                            rf = outerprod(rij, accijmag*rij)
                            do ixyz = 1,3
                                rfv(ixyz,1) = dot_product(rf(ixyz,:),vi_t(:))
                            enddo
                            call pressure_tensor_forces_VA(ri, rj, rFv, domain, & 
                                                           nbins, nhb, rfvbin,  & 
                                                           VA_heatflux_calcmethod, &
                                                           VA_heatflux_line_samples)
                        endif

					endif
				enddo
			enddo
			enddo
			enddo
			currenti => oldi
			oldi => currenti%next !Use pointer in datatype to obtain next item in list
		enddo
	enddo
	enddo
	enddo

    ! Add FENE contribution if it's there
    if (potential_flag .eq. 1) then
        call add_FENE_contribution
    end if

	nullify(oldi)      	!Nullify as no longer required
	nullify(oldj)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

contains

    subroutine add_FENE_contribution
        use polymer_info_MD
        use Volume_average_pressure, only : pressure_tensor_forces_VA
	    use librarymod, only: outerprod
        implicit none

        integer :: b

        do molnoi=1,np+halo_np

            ri(:) = r(:,molnoi) !Retrieve ri(:)
            do b=1,monomer(molnoi)%funcy

                molnoj = bond(b,molnoi)
                if (molnoj.eq.0) cycle

                rj(:)  = r(:,molnoj)
                rij(:) = ri(:) - rj(:)
                rij2   = dot_product(rij,rij)
                accijmag = -k_c/(1-(rij2/(R_0**2)))	!(-dU/dr)*(1/|r|)
                rf = outerprod(rij,rij*accijmag)
                call pressure_tensor_forces_VA(ri, rj, rf, domain,  & 
                                                nbins, nhb, rfbin,  &
                                                VA_calcmethod=1,    & 
                                                VA_line_samples=VA_line_samples)

            end do	

        end do

    end subroutine

end subroutine simulation_compute_rfbins


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

        ! Straight line calculation (DT thinks this is not the correct
        ! way to do it)
!       rs_cart = ri + s*rij 
!
!		if (  .not. any(abs(rs_cart(:)).gt.halfdomain(:))  ) then
!
!			! Binning conventions 
!           rs = cpolariser(globalise(rs_cart))
!			rs(1)  = rs(1) - r_oi
!			rs(2)  = modulo(rs(2),2.d0*pi)
!			rs(3)  = rs_cart(3) + halfdomain(3) 
!
!			!Add to cylindrical bins
!			br = ceiling(rs(1)/VAbinsize(1)) 
!			bt = ceiling(rs(2)/VAbinsize(2)) 
!			bz = ceiling(rs(3)/VAbinsize(3)) + cpol_nhbz
!
!			!Ignore molecules in cylinder region
!			if ( br .ge. 1 .and. br .le. cpol_bins(1) ) then
!				rfbin(br,bt,bz,:,:) =  &
!				rfbin(br,bt,bz,:,:) + rF(:,:)/real(Ns,kind(0.d0))
!			end if
!
!       end if
!		s = s + ds	


	end do	
	
end subroutine pressure_tensor_forces_VA_trap_cpol


!Cylindrical polar version
subroutine simulation_compute_rfbins_cpol(imin, imax, jmin, jmax, kmin, kmax)
	use module_compute_forces
	use librarymod, only: cpolariser
	use messenger, only: globalise
	implicit none

	integer, intent(in) :: imin, jmin, kmin, imax, jmax, kmax

	integer :: i, j
	integer	:: icell, jcell, kcell
	integer :: icellshift, jcellshift, kcellshift
	integer :: cellnp, adjacentcellnp
	integer	:: molnoi, molnoj
	type(node), pointer :: oldi, currenti, oldj, currentj

	double precision,dimension(3)	:: cellsperbin

	!Calculate bin to cell ratio
	cellsperbin = 1.d0/binspercell !ceiling(ncells(1)/dble(nbins(1)))
	where (cellsperbin .gt. 1.d0) cellsperbin = 1.d0

	do kcell=(kmin-1)*cellsperbin(3)+1, kmax*cellsperbin(3)
	do jcell=(jmin-1)*cellsperbin(2)+1, jmax*cellsperbin(2)
	do icell=(imin-1)*cellsperbin(1)+1, imax*cellsperbin(1)
	
		cellnp = cell%cellnp(icell,jcell,kcell)
		oldi => cell%head(icell,jcell,kcell)%point

		do i = 1,cellnp

			molnoi = oldi%molno
			ri = r(:,molnoi)

			do kcellshift = -1,1
			do jcellshift = -1,1
			do icellshift = -1,1

				!Prevents out of range values in i
				if (icell+icellshift .lt. imin) cycle
				if (icell+icellshift .gt. imax) cycle
				!Prevents out of range values in j
				if (jcell+jcellshift .lt. jmin) cycle
				if (jcell+jcellshift .gt. jmax) cycle
				!Prevents out of range values in k
				if (kcell+kcellshift .lt. kmin) cycle
				if (kcell+kcellshift .gt. kmax) cycle

				oldj => cell%head(icell+icellshift, &
				                  jcell+jcellshift, &
				                  kcell+kcellshift) % point

				adjacentcellnp = cell%cellnp(icell+icellshift, &
				                             jcell+jcellshift, &
				                             kcell+kcellshift)

				do j = 1,adjacentcellnp

					molnoj = oldj%molno
					rj = r(:,molnoj)

					currentj => oldj
					oldj => currentj%next 

					! Prevent interaction with self
					if ( molnoi == molnoj ) cycle 

					rij(:) = ri(:) - rj(:)
					rij2 = dot_product(rij,rij)

					if (rij2 < rcutoff2) then

						invrij2 = 1.d0/rij2
						accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)

						call pressure_tensor_forces_VA_trap_cpol(ri,rj,accijmag)

					endif

				enddo

			enddo
			enddo
			enddo

			currenti => oldi
			oldi => currenti%next

		enddo

	enddo
	enddo
	enddo

    ! Add FENE contribution if it's there
    if (potential_flag .eq. 1) then
        call add_FENE_contribution
    end if

	nullify(oldi)
	nullify(oldj)
	nullify(currenti)
	nullify(currentj)

contains

    subroutine add_FENE_contribution
        use polymer_info_MD
        implicit none

        integer :: b

        do molnoi=1,np

            ri(:) = r(:,molnoi) !Retrieve ri(:)
            do b=1,monomer(molnoi)%funcy

                molnoj = bond(b,molnoi)
                if (molnoj.eq.0) cycle

                rj(:)  = r(:,molnoj)
                rij(:) = ri(:) - rj(:)
                rij2   = dot_product(rij,rij)
                accijmag = -k_c/(1-(rij2/(R_0**2)))	!(-dU/dr)*(1/|r|)

                call pressure_tensor_forces_VA_trap_cpol(ri,rj,accijmag)

            end do	

        end do

    end subroutine

end subroutine simulation_compute_rfbins_cpol




!===================================================================================
!		RECORD HEAT FLUX AT LOCATION IN SPACE
! Either by Volume Average style binning or virial expression for the whole domain
!===================================================================================

subroutine heatflux_averaging(ixyz)
	use field_io, only : VA_heatflux_io
	use module_record
	implicit none

	integer			:: ixyz
	integer, save	:: sample_count = 0

	sample_count = sample_count + 1
	call cumulative_heatflux(ixyz,sample_count)
	if (sample_count .eq. Nheatflux_ave) then
		sample_count = 0

		select case(ixyz)
		case(1)
		!FULL DOMAIN VIRIAL heatflux CALCULATION
            stop "Error - no global heatflux calculation"
		case(2)
		!VA heatflux CALCULATION
			call VA_heatflux_io
			evbin  = 0.d0
			rfvbin  = 0.d0
            heatfluxbin = 0.d0
		case(3)
            stop "Error - no cpol bins heatflux calculation"
		case default 
			call error_abort("Average heatflux Binning Error")
		end select
	endif

end subroutine heatflux_averaging

!===================================================================================
!Add heatflux_tensor to running total

subroutine cumulative_heatflux(ixyz,sample_count)
	use module_record
    use messenger_data_exchange, only : globalSum
	implicit none

	integer								:: sample_count,n,ixyz,jxyz,kxyz
	double precision, dimension(3)		:: velvect
	double precision, dimension(3)      :: rglob
	double precision, dimension(3,3)	:: Pxytemp

	Pxytemp = 0.d0

	select case(ixyz)
	case(1)
        stop "Error - Cumulative heatflux  "
	case(2)
		!VOLUME AVERAGE heatflux CALCULATION
		!Don't calculate Position (x) force dot v for configurational part of heatflux tensor
        !if it has already been calculated in pressure calculation
        if (pressure_outflag .ne. 2 .or. Nheatflux_ave .ne. Nstress_ave) then
            call simulation_compute_rfbins
        endif

		!Calculate velocity (x) velocity for kinetic part of heatflux tensor
        call  simulation_compute_energy_VA(1,nbins(1), 1,nbins(2), 1,nbins(3))

		!Add results to cumulative total
		heatfluxbin(:,:,:,:)  =     evbin(:,:,:,:) 				                & 
						         + rfvbin( 1+nhb(1):nbins(1)+nhb(1), 			& 
							      	   	   1+nhb(2):nbins(2)+nhb(2), 			& 
							      		   1+nhb(3):nbins(3)+nhb(3),:,1)/2.d0

	case default 
		call error_abort("Cumulative heatflux Averaging Error")
	end select

end subroutine cumulative_heatflux

!----------------------------------------------------------------------------------
!Compute kinetic part of heatflux tensor

subroutine simulation_compute_energy_VA(imin,imax,jmin,jmax,kmin,kmax)
	use module_record
	use physical_constants_MD
	implicit none

	integer								:: imin, jmin, kmin, imax, jmax, kmax
	integer         					:: n, ixyz,jxyz
	integer 							:: ibin, jbin, kbin
	integer 							:: bin(3)
    double precision                    :: energy
	double precision, dimension(3)		:: VAbinsize, velvect

	!Determine bin size
	VAbinsize(:) = domain(:) / nbins(:)

	! Add kinetic part of heatflux tensor for all molecules
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

        energy = 0.5d0 * (dot_product(velvect,velvect) + potenergymol(n))
		evbin(ibin,jbin,kbin,:) = evbin(ibin,jbin,kbin,:) & 
                                  + velvect(:) * energy

	enddo

end subroutine simulation_compute_energy_VA


subroutine bforce_pdf_stats
    use boundary_MD
    use statistics_io, only: bforce_pdf_write
    implicit none

	integer, save		:: average_count=-1

    average_count = average_count + 1
    call collect_bforce_pdf_data

    if (average_count .eq. bforce_pdf_Nave) then
        average_count = 0
        call bforce_pdf_write
    end if 

end subroutine bforce_pdf_stats

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

subroutine mass_flux_averaging(flag)
	!use field_io, only : mass_flux_io
	use module_record
	use CV_objects, only : CVcheck_mass, CV_debug!,CV_sphere_mass
	implicit none

	integer			                :: flag
    integer,dimension(3)            :: thermbinstop,thermbinsbot
    double precision,dimension(3)   :: mbinsize
	integer, save	                :: sample_count

	!Only average if mass averaging turned on
	if (flag .eq. 0) return

	call cumulative_mass_flux
	sample_count = sample_count + 1
	if (sample_count .eq. Nmflux_ave) then
		if (CV_debug .ne. 0) then
    		mbinsize(:) = domain(:) / nbins(:)
            thermbinstop = ceiling(thermstattop/mbinsize)
            thermbinsbot = ceiling(thermstatbottom/mbinsize)
            
		    call CVcheck_mass%check_error(1+nhb(1)+thermbinsbot(1),nbins(1)+nhb(1)-thermbinstop(1), & 
										  1+nhb(2)+thermbinsbot(2),nbins(2)+nhb(2)-thermbinstop(2), & 
										  1+nhb(3)+thermbinsbot(3),nbins(3)+nhb(3)-thermbinstop(3),iter,irank)
			!call CV_sphere_mass%check_error(1,1,1,1,1,1,iter,irank)
	    endif
		call mass_flux_io
		sample_count = 0
		mass_flux = 0
		call mass_snapshot
	endif

end subroutine mass_flux_averaging

!===================================================================================
! Mass Flux over a surface of a bin
! Includes all intermediate bins

subroutine cumulative_mass_flux
	use module_record
    use librarymod, only : imaxloc, heaviside  =>  heaviside_a1
    !use CV_objects, only : CV_sphere_mass
    implicit none

	integer							:: jxyz,i,j,k,n
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	double precision				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	double precision,dimension(3)	:: mbinsize,crossface
	double precision,dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	do n = 1,np

		ri1(:) = r(:,n) 							!Molecule i at time t
		ri2(:) = r(:,n)	- delta_t*v(:,n)			!Molecule i at time t-dt
		ri12   = ri1 - ri2							!Molecule i trajectory between t-dt and t
		where (ri12 .eq. 0.d0) ri12 = 0.000001d0
		
		!call CV_sphere_mass%Add_spherical_CV_fluxes(ri2,ri1)

		!Assign to bins before and after using integer division
        ibin1(:) =  get_bin(ri1)
        ibin2(:) =  get_bin(ri2)

!		ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb(:)
!		ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb(:)

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
				      + nint(dble(onfacext)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),5) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),5) &
				      + nint(dble(onfaceyt)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),6) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),6) &
				      + nint(dble(onfacezt)*abs(crossface(jxyz)))
				      

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
	use CV_objects, only : CVcheck_mass, CV_debug!, CV_sphere_mass
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
		ibin(:) = get_bin(r(:,n)) 
		!ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb
		volume_mass_temp(ibin(1),ibin(2),ibin(3)) = volume_mass_temp(ibin(1),ibin(2),ibin(3)) + 1
		!call  CV_sphere_mass%Add_spherical_CV_mass(r(:,n))
	enddo

	!Output Control Volume momentum change and fluxes
	call mass_bin_io(volume_mass_temp,'snap')
	!Create copy of previous timestep Control Volume mass and calculate time evolution
	if (CV_debug .ne. 0) then
		call CVcheck_mass%update_dXdt(volume_mass_temp)
		!call CV_sphere_mass%update_dXdt(CV_sphere_mass%Xtemp)
	endif

	deallocate(volume_mass_temp)

end subroutine mass_snapshot


!===================================================================================
! Control Volume Momentum continuity
!===================================================================================


!===================================================================================
! Momentum Flux over a surface of a bin including all intermediate bins

module cumulative_momentum_flux_mod

contains

subroutine cumulative_momentum_flux(r_,v_,momentum_flux_,notcrossing)
	use module_record, only : vflux_outflag, domain, halfdomain, planespacing, CV_debug, & 
							  delta_t, planes, Pxy_plane, nplanes, np, nbins, nhb, iter
	use CV_objects, only : CV_constraint!, CV_sphere_momentum
    use librarymod, only : imaxloc, heaviside  =>  heaviside_a1
	use interfaces, only : error_abort
    use module_record, only : get_bin
	implicit none

	double precision,dimension(:,:),allocatable,intent(in) 			:: r_,v_
	double precision,dimension(:,:,:,:,:),allocatable,intent(inout) :: momentum_flux_
	integer,dimension(:),allocatable,intent(out),optional			:: notcrossing


	integer							:: ixyz,jxyz,i,j,k,n
	integer							:: planeno
	integer		,dimension(3)		:: ibin1,ibin2,cbin
    double precision				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	double precision				:: crossplane,rplane,shift
	double precision,dimension(3)	:: mbinsize,velvect,crossface
	double precision,dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

	!Allocate array if required and assume all are not crossing
	if (present(notcrossing)) then
		if (allocated(notcrossing)) deallocate(notcrossing)
		allocate(notcrossing(np)); notcrossing = 1
	endif

	ixyz = vflux_outflag

	select case(vflux_outflag)
	case(1:3)
		!MOP momentum flux
		!Shift by half difference between value rounded down and actual value
		!to ensure same distance to top and bottom plane from domain edge
		shift = 0.5d0*(domain(ixyz) - planespacing * (nplanes-1))

		!Add molecular velocities to the configuration stresses
		do n=1,np

			!Replace Signum function with this functions which gives a
			!check for plane crossing and the correct sign 
			crossplane = ceiling((r_(ixyz,n)+halfdomain(ixyz)-shift)/planespacing) & 
				    	-ceiling((r_(ixyz,n)-delta_t*v_(ixyz,n)+halfdomain(ixyz)-shift)/planespacing)

			if (crossplane .ne. 0) then

				!Obtain nearest plane number by integer division 
				!and retrieve location of plane from array
				planeno = ceiling((r_(ixyz,n)+halfdomain(ixyz)-shift) 	& 
					  /planespacing)-heaviside(dble(crossplane))+1
				if (planeno .lt. 1) planeno = 1
				if (planeno .gt. nplanes) planeno = nplanes
				rplane = planes(planeno)

				!Calculate velocity at time of intersection
				!crosstime = (r_(ixyz,n) - rplane)/v_(ixyz,n)
				velvect(:) = v_(:,n) !- a(:,n) * crosstime

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

			!Get velocity at v_(t+dt/2) from v_(t-dt/2)
			velvect(:) = v_(:,n)
			ri1(:) = r_(:,n) 							!Molecule i at time t
			ri2(:) = r_(:,n)	- delta_t*velvect			!Molecule i at time t-dt
			ri12   = ri1 - ri2							!Molecule i trajectory between t-dt and t
			where (ri12 .eq. 0.d0) ri12 = 0.000001d0

			!call CV_sphere_momentum%Add_spherical_CV_fluxes(velvect,ri2,ri1)

			!Assign to bins before and after using integer division
            ibin1(:) =  get_bin(ri1)
            ibin2(:) =  get_bin(ri2)

!			ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb(:)
!			ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb(:)

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
					!crosstime = (r_(jxyz,n) - rplane)/v_(jxyz,n)
					velvect(:) = v_(:,n) !- a(:,n) * crosstime
					!Change in velocity at time of crossing is not needed as velocity assumed constant 
					!for timestep and changes when forces are applied.
!					if (iter .eq. 1737) then
!					if (cbin(1) .eq. 2 .and. cbin(2) .eq. 7 .and. cbin(3) .eq. 2 .or. &
!						cbin(1) .eq. 2 .and. cbin(2) .eq. 7 .and. cbin(3) .eq. 3 .or. &
!						cbin(1) .eq. 3 .and. cbin(2) .eq. 7 .and. cbin(3) .eq. 3 ) then
!						if (abs(onfacexb) .gt. 0.000001) print'(a,i8,3i4,f16.10,a,f16.10,a,f16.10,a,2f16.10)', & 
!							'Flux botx', n,cbin, ri2(1), ' => ', binbot(1), ' => ', ri1(1), ' v= ', velvect(1),sign(1.d0,onfacexb)
!						if (abs(onfaceyb) .gt. 0.000001) print'(a,i8,3i4,f16.10,a,f16.10,a,f16.10,a,2f16.10)', & 
!							'Flux boty', n,cbin, ri2(2), ' => ', binbot(2), ' => ', ri1(2), ' v= ', velvect(2),sign(1.d0,onfaceyb)
!						if (abs(onfacezb) .gt. 0.000001) print'(a,i8,3i4,f16.10,a,f16.10,a,f16.10,a,2f16.10)', & 
!							'Flux botz', n,cbin, ri2(3), ' => ', binbot(3), ' => ', ri1(3), ' v= ', velvect(3),sign(1.d0,onfacezb)
!						if (abs(onfacext) .gt. 0.000001) print'(a,i8,3i4,f16.10,a,f16.10,a,f16.10,a,2f16.10)', & 
!							'Flux topx', n,cbin, ri2(1), ' => ', bintop(1), ' => ', ri1(1), ' v= ', velvect(1),sign(1.d0,onfacext)
!						if (abs(onfaceyt) .gt. 0.000001) print'(a,i8,3i4,f16.10,a,f16.10,a,f16.10,a,2f16.10)', & 
!							'Flux topy', n,cbin, ri2(2), ' => ', bintop(2), ' => ', ri1(2), ' v= ', velvect(2),sign(1.d0,onfaceyt)
!						if (abs(onfacezt) .gt. 0.000001) print'(a,i8,3i4,f16.10,a,f16.10,a,f16.10,a,2f16.10)', & 
!							'Flux topz', n,cbin, ri2(3), ' => ', bintop(3), ' => ', ri1(3), ' v= ', velvect(3),sign(1.d0,onfacezt)
!					endif
!					endif

					!Add Momentum flux over face
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,1) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,1) & 
					      + velvect(:)*dble(onfacexb)*abs(crossface(jxyz))
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,2) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,2) & 
					      + velvect(:)*dble(onfaceyb)*abs(crossface(jxyz))
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,3) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,3) &
					      + velvect(:)*dble(onfacezb)*abs(crossface(jxyz))
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,4) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,4) &
					      + velvect(:)*dble(onfacext)*abs(crossface(jxyz))
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,5) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,5) &
					      + velvect(:)*dble(onfaceyt)*abs(crossface(jxyz))
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,6) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,6) &
					      + velvect(:)*dble(onfacezt)*abs(crossface(jxyz))

				enddo
				enddo
				enddo

				!Record mask of molecules which are currently crossing
 				if (present(notcrossing)) notcrossing(n) = 0

			endif

		enddo
	case default 
		call error_abort("Cumulative Momentum flux Error")
	end select

end subroutine cumulative_momentum_flux

end module cumulative_momentum_flux_mod

subroutine momentum_flux_averaging(flag)
	!use field_io, only :  momentum_flux_io,surface_stress_io, & 
	!					  external_force_io,MOP_stress_io
	use module_record
	use cumulative_momentum_flux_mod, only : cumulative_momentum_flux
	use CV_objects, only : CV_debug, CVcheck_momentum!, CV_constraint!, CV_sphere_momentum
	implicit none

	integer,intent(in)	:: flag
	integer				::icell,jcell,kcell,n
    integer,dimension(3):: thermbinstop,thermbinsbot
	double precision,dimension(3)	:: mbinsize
	integer, save		:: sample_count

	if (flag .eq. 0) return

	call cumulative_momentum_flux(r,v,momentum_flux)
	sample_count = sample_count + 1
	if (sample_count .eq. Nvflux_ave) then

		select case(flag)
		case(1:3)
			!MOP momentum flux and stresses
			call MOP_stress_io(flag)
			Pxy_plane = 0.d0
		case(4)
			!CV momentum flux and stress
			call momentum_flux_io
			momentum_flux = 0.d0
			call momentum_snapshot
			if (external_force_flag .ne. 0 .or. & 
				ensemble .eq. tag_move     .or. & 
				CVforce_flag .ne. VOID) then
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
        if (CV_debug .ne. 0) then

            !Currently, CV check will not work for thermostatted molecules
            !so exclude these from checked region
    		mbinsize(:) = domain(:) / nbins(:)
            if (ensemble .eq. 6) then
                thermbinstop = ceiling(thermstattop/mbinsize)
                thermbinsbot = ceiling(thermstatbottom/mbinsize)
            else
                thermbinstop = 0
                thermbinsbot = 0
            endif
	    	    call CVcheck_momentum%check_error(1+nhb(1)+thermbinsbot(1),nbins(1)+nhb(1)-thermbinstop(1), & 
	    										  1+nhb(2)+thermbinsbot(2),nbins(2)+nhb(2)-thermbinstop(2), & 
	    										  1+nhb(3)+thermbinsbot(3),nbins(3)+nhb(3)-thermbinstop(3),iter,irank)
        endif
	endif

end subroutine momentum_flux_averaging

!===================================================================================
! Control Volume snapshot of momentum in a given bin

subroutine momentum_snapshot
	use field_io, only : velocity_bin_io
	use module_record
	!use CV_objects, only : CV_sphere_momentum
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
        ibin(:) =  get_bin(r(:,n))
		!ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb

		volume_mass_temp(ibin(1),ibin(2),ibin(3)) = volume_mass_temp(ibin(1),ibin(2),ibin(3)) + 1
		volume_momentum_temp(ibin(1),ibin(2),ibin(3),:) = volume_momentum_temp(ibin(1),ibin(2),ibin(3),:) + v(:,n)
		!call CV_sphere_momentum%Add_spherical_CV_velocity(r(:,n),v(:,n))
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


!===================================================================================
! Energy Flux over a surface of a bin including all intermediate bins

module cumulative_energy_flux_mod

contains

subroutine cumulative_energy_flux(r_,v_,energy_flux_)
	use module_record, only : eflux_outflag, domain, halfdomain, planespacing, CV_debug, & 
							  delta_t, planes, Pxyv_plane, nplanes, np, nbins, nhb, potenergymol,potenergymol_mdt,a
	use CV_objects, only : CVcheck_energy !, CV_sphere_momentum
    use librarymod, only : imaxloc, heaviside  =>  heaviside_a1
	use interfaces, only : error_abort
    use module_record, only : get_bin
	implicit none

	double precision,dimension(:,:),allocatable,intent(in) 			:: r_,v_
	double precision,dimension(:,:,:,:),allocatable,intent(inout) :: energy_flux_

	integer							:: ixyz,jxyz,i,j,k,n,planeno
	double precision				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	double precision				:: crosstime,crossplane,rplane,shift,energy
	double precision				:: crosstimetop,crosstimebot,frac,delta_t_cross,potenergy_cross
	double precision,dimension(3)	:: mbinsize,velvect,crossface
	double precision,dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

	ixyz = eflux_outflag

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
			crossplane = ceiling((r_(ixyz,n)+halfdomain(ixyz)-shift)/planespacing) & 
				    	-ceiling((r_(ixyz,n)-delta_t*v_(ixyz,n)+halfdomain(ixyz)-shift)/planespacing)

			if (crossplane .ne. 0) then

				!Obtain nearest plane number by integer division 
				!and retrieve location of plane from array
				planeno = ceiling((r_(ixyz,n)+halfdomain(ixyz)-shift) 	& 
					  /planespacing)-heaviside(dble(crossplane))+1
				if (planeno .lt. 1) planeno = 1
				if (planeno .gt. nplanes) planeno = nplanes
				rplane = planes(planeno)

				!Calculate energy at intersection
				velvect(:) = v_(:,n) + 0.5d0*a(:,n)*delta_t
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

			ri1(:) = r_(:,n)					!Molecule i at time  t
			ri2(:) = r_(:,n)-delta_t*v_(:,n)	!Molecule i at time t-dt
			ri12   = ri1 - ri2				!Molecule i trajectory between t-dt and t
			where (ri12 .eq. 0.d0) ri12 = 0.000001d0

			!Assign to bins before and after using integer division
            ibin1(:) =  get_bin(ri1)
            ibin2(:) =  get_bin(ri2)

			!ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb(:)
			!ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb(:)

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
!  					crosstimetop = (bintop(jxyz) - ri2(jxyz))/ v_(jxyz,n)
!  					crosstimebot = (binbot(jxyz) - ri2(jxyz))/ v_(jxyz,n)
!  					if (crosstimetop .gt. 0.d0 .and. crosstimetop .lt. delta_t) then
!  						delta_t_cross = crosstimetop
! 	 					frac = delta_t_cross/delta_t
!  					elseif (crosstimebot .gt. 0.d0 .and. crosstimebot .lt. delta_t) then
!  						delta_t_cross = crosstimebot
!  					else
!  						stop "Error -- crossing doesn't happen in energy flux"
!  					endif
!  					velvect(:) = v_(:,n) + 0.5d0*a(:,n)*delta_t_cross
!  					potenergy_cross =   potenergymol(n)
! 					frac = delta_t_cross/delta_t
!  					potenergy_cross =   potenergymol(n)*frac  & 
!  									  + potenergymol_mdt(n)*(1-frac)
!  					energy = 0.5d0 * (dot_product(velvect,velvect) + potenergy_cross)

					velvect(:) = v_(:,n) + 0.5d0*a(:,n)*delta_t
					energy = 0.5d0 * (dot_product(velvect,velvect) + potenergymol(n))


! 					print'(a,8f10.5)', 'Energy flux', delta_t_cross,potenergymol(n),potenergymol_mdt(n), & 
! 													  frac, dot_product(velvect,velvect), & 
! 													dot_product(v_(:,n) + 0.5d0*a(:,n)*delta_t,v_(:,n) + 0.5d0*a(:,n)*delta_t), & 
! 													energy,0.5d0 * (dot_product(velvect,velvect) + potenergymol(n))

					!print'(2(a,f16.12))', ' V^2 ', 0.5d0 * dot_product(velvect,velvect) ,' potential_energy ',potenergymol(n)

					!Change in velocity at time of crossing is not needed as velocity assumed constant 
					!for timestep and changes when forces are applied.

					!Add Energy flux over face
					energy_flux_(cbin(1),cbin(2),cbin(3),1) = & 
						energy_flux_(cbin(1),cbin(2),cbin(3),1) & 
					      + energy*dble(onfacexb)*abs(crossface(jxyz))
					energy_flux_(cbin(1),cbin(2),cbin(3),2) = & 
						energy_flux_(cbin(1),cbin(2),cbin(3),2) & 
					      + energy*dble(onfaceyb)*abs(crossface(jxyz))
					energy_flux_(cbin(1),cbin(2),cbin(3),3) = & 
						energy_flux_(cbin(1),cbin(2),cbin(3),3) &
					      + energy*dble(onfacezb)*abs(crossface(jxyz))
					energy_flux_(cbin(1),cbin(2),cbin(3),4) = & 
						energy_flux_(cbin(1),cbin(2),cbin(3),4) &
					      + energy*dble(onfacext)*abs(crossface(jxyz))
					energy_flux_(cbin(1),cbin(2),cbin(3),5) = & 
						energy_flux_(cbin(1),cbin(2),cbin(3),5) &
					      + energy*dble(onfaceyt)*abs(crossface(jxyz))
					energy_flux_(cbin(1),cbin(2),cbin(3),6) = & 
						energy_flux_(cbin(1),cbin(2),cbin(3),6) &
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

end module cumulative_energy_flux_mod


!===================================================================================
! Collect Energy Flux over a surface and write data after specified time period

subroutine energy_flux_averaging(flag)
	!use field_io, only : energy_flux_io,surface_power_io,MOP_energy_io
	use module_record
	use CV_objects, only :  CVcheck_energy
	use cumulative_energy_flux_mod, only : cumulative_energy_flux
	implicit none

	integer,intent(in)	:: flag
	integer, save		:: sample_count

    integer,dimension(3):: thermbinstop,thermbinsbot
	double precision,dimension(3)	:: ebinsize

	if (flag .eq. 0) return

	!Integration of surface power using trapizium rule, get current value and
    !add to segment to running total using previous value
	!call simulation_compute_power(1, nbins(1), 1, nbins(2), 1, nbins(3))

	!call simulation_compute_power(1,nbins(1)+2*nhb(1), & 
    !                              1,nbins(2)+2*nhb(2), & 

	call simulation_compute_power!(1+nhb(1), nbins(1)+nhb(1), & 
                                 ! 1+nhb(2), nbins(2)+nhb(2), & 
                                 ! 1+nhb(3), nbins(3)+nhb(3))

    !Integration of stress using trapizium rule requires multiplication by timestep
	Pxyvface_integrated = Pxyvface_integrated + 0.5d0 * (Pxyvface_mdt + Pxyvface) 
	Pxyvface_mdt = Pxyvface
	Pxyvface = 0.d0

	call cumulative_energy_flux(r,v,energy_flux)
	sample_count = sample_count + 1
	if (sample_count .eq. Neflux_ave) then

		select case(flag)
		case(1:3)
			!MOP energy flux and Power (stresses*velocity)
			call MOP_energy_io(flag)
			Pxy_plane = 0.d0
		case(4)
			!CV energy flux and Power (stresses*velocity)
			call energy_flux_io
			energy_flux = 0.d0
			call energy_snapshot
			if (external_force_flag .ne. 0 .or. & 
				ensemble .eq. tag_move     .or. & 
				CVforce_flag .ne. VOID) then
				call external_forcev_io
				Fv_ext_bin = 0.d0
			endif
		case default 
			call error_abort("Energy flux averaging Error")
		end select

		sample_count = 0

	endif

	!Write forces out at time t before snapshot/final fluxes
	!as both use velocity at v(t-dt/2)
	if (sample_count .eq. Neflux_ave-1) then
		call surface_power_io
        Pxyvface_integrated = 0.d0

		!Debug flag to check CV conservation in code
		if (CV_debug .ne. 0) then
    		ebinsize(:) = domain(:) / nbins(:)
            if (ensemble .eq. 6) then
                thermbinstop = ceiling(thermstattop/ebinsize)
                thermbinsbot = ceiling(thermstatbottom/ebinsize)
            else
                thermbinstop = 0
                thermbinsbot = 0
            endif
		    call CVcheck_energy%check_error(1+nhb(1)+thermbinsbot(1),nbins(1)+nhb(1)-thermbinstop(1), & 
											1+nhb(2)+thermbinsbot(2),nbins(2)+nhb(2)-thermbinstop(2), & 
											1+nhb(3)+thermbinsbot(3),nbins(3)+nhb(3)-thermbinstop(3),iter,irank)
	   endif
	endif

end subroutine energy_flux_averaging

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
        ibin(:) = get_bin(r(:,n))
		!ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb(:)
		velvect(:) = v(:,n) + 0.5d0*a(:,n)*delta_t
		energy = 0.5d0 * (dot_product(velvect,velvect) + potenergymol(n))

		!if (all(ibin .eq. 3)) then
		!	print'(a,2i5,6f12.7)','E__ in bin 3', iter, n, velvect,0.5d0*(dot_product(velvect,velvect)),0.5d0*potenergymol(n), energy
		!endif
		!if (abs(sum(velvect)) .gt. 0.000001d0 .or. abs(potenergymol(n)) .gt. 0.000001d0) then
		!	print'(5i6,a,3f16.12,a,f16.12)', iter, n, ibin, ' velvect ', velvect ,' potential_energy ',potenergymol(n)
		!endif
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







!===================================================================================
!Forces over the surface of a Volume

subroutine control_volume_forces(fij,ri,rj,molnoi,molnoj)
    use module_record
    use librarymod, only : heaviside  =>  heaviside_a1
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
! Stresses over each of the six surfaces of the cuboid

subroutine control_volume_stresses(fij,ri,rj)
    use module_record
	use CV_objects, only : CV_debug,CV_constraint
    use librarymod, only : heaviside  =>  heaviside_a1
    implicit none


	double precision,intent(in),dimension(3)	:: ri,rj,fij

	integer							:: i,j,k,ixyz,face
	integer,dimension(3)			:: cbin, ibin, jbin
    double precision				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb,fijvi,fijvj
	double precision,dimension(3)	:: rij,fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb,velvect
	double precision,dimension(3)	:: Fbinsize, bintop, binbot, vi_t,vj_t,vi_tmdt,vj_tmdt

	!double precision,allocatable,dimension(:,:,:),save 	:: fij_dmt
	!double precision,allocatable,dimension(:,:,:),save 	:: fij_dmt

	!if (.not. allocated(fij_dmt)) then
	!	allocate(fij_dmt(3,np+extralloc,np+extralloc))
	!	fij_dmt = 0.d0
	!endif

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
		if (CVforce_flag .ne. VOID .and. iter-initialstep+1 .ge. CVforce_starttime) then
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,1) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*dble(onfacexb)
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,2) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*dble(onfaceyb)
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,3) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*dble(onfacezb)
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,4) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*dble(onfacext)
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,5) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*dble(onfaceyt)
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,6) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*dble(onfacezt)
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
! Stresses times velocity over each of the six surfaces of the cuboid

subroutine control_volume_power(fij,ri,rj,vi_t)
    use module_record
	use CV_objects, only : CV_debug,CV_constraint
    use librarymod, only : heaviside  =>  heaviside_a1
    implicit none


	double precision,intent(in),dimension(3)	:: ri,rj,fij
	double precision,dimension(3),intent(in)	:: vi_t


	integer							:: i,j,k,ixyz,face
	integer,dimension(3)			:: cbin, ibin, jbin
    double precision				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb,fijvi,fijvj
	double precision,dimension(3)	:: rij,fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb,velvect
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
	ibin(:) = get_bin(ri) !ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
	jbin(:) = get_bin(rj) !ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin

	if (ibin(1) .eq. jbin(1) .and. ibin(2) .eq. jbin(2) .and. ibin(3) .eq. jbin(3)) return
		
	do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
	do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
	do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

		cbin(1) = i; cbin(2) = j; cbin(3) = k

		bintop(:) = (cbin(:)-nhb(:)  )*Fbinsize(:)-halfdomain(:)
		binbot(:) = (cbin(:)-nhb(:)-1)*Fbinsize(:)-halfdomain(:)

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

		!Stress work acting on face over volume - add current value 
        ! to use in trapizium rule later 
		fijvi = dot_product(fij,vi_t)

		Pxyvface(cbin(1),cbin(2),cbin(3),1) = & 
			Pxyvface(cbin(1),cbin(2),cbin(3),1) + fijvi*onfacexb
		Pxyvface(cbin(1),cbin(2),cbin(3),2) = & 
			Pxyvface(cbin(1),cbin(2),cbin(3),2) + fijvi*onfaceyb
		Pxyvface(cbin(1),cbin(2),cbin(3),3) = &
			Pxyvface(cbin(1),cbin(2),cbin(3),3) + fijvi*onfacezb
		Pxyvface(cbin(1),cbin(2),cbin(3),4) = &
			Pxyvface(cbin(1),cbin(2),cbin(3),4) + fijvi*onfacext
		Pxyvface(cbin(1),cbin(2),cbin(3),5) = &
			Pxyvface(cbin(1),cbin(2),cbin(3),5) + fijvi*onfaceyt
		Pxyvface(cbin(1),cbin(2),cbin(3),6) = &
			Pxyvface(cbin(1),cbin(2),cbin(3),6) + fijvi*onfacezt

	enddo
	enddo
	enddo

end subroutine control_volume_power



!===================================================================================
!Forces times velocity over the surface of a Volume

module get_timesteps_module

contains

subroutine get_timesteps(ncrossings,bin,bin_mdt,Fbinsize,rc,vc,delta_t,delta_t_list,count_t)
	use computational_constants_MD, only : halfdomain, nhb, iter
	implicit none

	integer,intent(in)								:: ncrossings, bin(3),bin_mdt(3)
	integer,intent(inout)							:: count_t
	double precision,intent(in)						:: delta_t
	double precision,dimension(3),intent(in)		:: rc,vc,Fbinsize
	double precision,dimension(:),allocatable,intent(inout)	:: delta_t_list

	integer											:: i,j,k,m
	double precision,dimension(3)					:: bintop,binbot
	double precision,dimension(6)					:: crosstime


	if (ncrossings .eq. 0) return

       !Otherwise, get time spent in each cell
	do i = bin(1),bin_mdt(1),sign(1,bin_mdt(1)-bin(1))
    do j = bin(2),bin_mdt(2),sign(1,bin_mdt(2)-bin(2))
    do k = bin(3),bin_mdt(3),sign(1,bin_mdt(3)-bin(3))

		!Get bin top and bottom
	    bintop(1) = (i-1*nhb(1)  )*Fbinsize(1)-halfdomain(1)
	    bintop(2) = (j-1*nhb(2)  )*Fbinsize(2)-halfdomain(2)
	    bintop(3) = (k-1*nhb(3)  )*Fbinsize(3)-halfdomain(3)
	    binbot(1) = (i-1*nhb(1)-1)*Fbinsize(1)-halfdomain(1)
	    binbot(2) = (j-1*nhb(2)-1)*Fbinsize(2)-halfdomain(2)
	    binbot(3) = (k-1*nhb(3)-1)*Fbinsize(3)-halfdomain(3)

	    !Calculate the plane intersect of line with surfaces of the cube
		crosstime(1) = (rc(1) - bintop(1))/ vc(1)
		crosstime(2) = (rc(1) - binbot(1))/ vc(1)
		crosstime(3) = (rc(2) - bintop(2))/ vc(2)
		crosstime(4) = (rc(2) - binbot(2))/ vc(2)
		crosstime(5) = (rc(3) - bintop(3))/ vc(3)
		crosstime(6) = (rc(3) - binbot(3))/ vc(3)

! 		print'(a,i6,3i3,6f12.7)','Crossingtimes',iter,i,j,k,crosstime
!  		print'(a,i6,3i3,9f12.7)','Surfaces',iter,i,j,k,	bintop(1),rc(1),binbot(1), & 
!  														bintop(2),rc(2),binbot(2), & 
!  														bintop(3),rc(3),binbot(3) 

		!Add any crossings within the time period to the list of crossings
		do m = 1,6
			if (crosstime(m) .gt. 0.d0 .and. crosstime(m) .lt. delta_t) then
				count_t = count_t + 1
				delta_t_list(count_t) = crosstime(m)
			endif
		enddo

	enddo
	enddo
	enddo

end subroutine get_timesteps


subroutine get_CV_surface_contributions(ibin,jbin,ri,rj,Fbinsize,value,CV_Face_value,molnoi,molnoj)
	use computational_constants_MD, only : halfdomain, nhb, iter
    use librarymod, only : heaviside => heaviside_a1
	implicit none

	integer, dimension(3),intent(in)								:: ibin,jbin
	integer,intent(in)												:: molnoi,molnoj !TEMPTEMP
	double precision, dimension(3),intent(in)						:: ri, rj, Fbinsize
	double precision, intent(in)									:: value
	double precision,dimension(:,:,:,:),allocatable, intent(inout)	:: CV_Face_value

	integer										:: i,j,k,ixyz, face
    double precision							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	double precision,dimension(3)   			:: Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
	double precision,dimension(3)				:: rij,  bintop, binbot


	!If same bin, nothing to do here
	if (ibin(1) .eq. jbin(1) .and. ibin(2) .eq. jbin(2) .and. ibin(3) .eq. jbin(3)) return

	!Get interaction line
	rij = ri - rj

	!Prevent Division by zero
	do ixyz = 1,3
		if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
	enddo
		
	!Loop through all intermediate bins, check surface to add to and then add
	do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
	do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
	do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

		!Get bin top and bottom
	    bintop(1) = (i-1*nhb(1)  )*Fbinsize(1)-halfdomain(1)
	    bintop(2) = (j-1*nhb(2)  )*Fbinsize(2)-halfdomain(2)
	    bintop(3) = (k-1*nhb(3)  )*Fbinsize(3)-halfdomain(3)
	    binbot(1) = (i-1*nhb(1)-1)*Fbinsize(1)-halfdomain(1)
	    binbot(2) = (j-1*nhb(2)-1)*Fbinsize(2)-halfdomain(2)
	    binbot(3) = (k-1*nhb(3)-1)*Fbinsize(3)-halfdomain(3)

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

		!Value acting on face
		CV_Face_value(i,j,k,1) = CV_Face_value(i,j,k,1) + value*onfacexb
		CV_Face_value(i,j,k,2) = CV_Face_value(i,j,k,2) + value*onfaceyb
		CV_Face_value(i,j,k,3) = CV_Face_value(i,j,k,3) + value*onfacezb
		CV_Face_value(i,j,k,4) = CV_Face_value(i,j,k,4) + value*onfacext
		CV_Face_value(i,j,k,5) = CV_Face_value(i,j,k,5) + value*onfaceyt
		CV_Face_value(i,j,k,6) = CV_Face_value(i,j,k,6) + value*onfacezt

  	if (i .eq. 3 .and. j .eq. 3 .and. k .eq. 3) then
   	if (any(abs((/onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb/)) .gt. 0.00001)) then
   	!if (abs(onfacext+onfacexb+onfaceyt+onfaceyb+onfacezt+onfacezb) .gt. 0.00001) then
   			if (abs(onfacexb) .gt. 0.000001) face = 1 
   			if (abs(onfaceyb) .gt. 0.000001) face = 2 
   			if (abs(onfacezb) .gt. 0.000001) face = 3 
   			if (abs(onfacext) .gt. 0.000001) face = 4 
   			if (abs(onfaceyt) .gt. 0.000001) face = 5 
   			if (abs(onfacezt) .gt. 0.000001) face = 6 
   				print'(a,6i5,2f12.4,6i3,f6.0,5f4.0)','Int pp box 3', iter,molnoi,molnoj,i,j,k,value, & 
   																	CV_Face_value(i,j,k,face), & 
   																	ibin,jbin,onfacext,onfacexb,onfaceyt, & 
   																	onfaceyb,onfacezt,onfacezb
   	!endif
   	endif
   	endif

		!print'(a,3i4,7f10.5)', 'IN surface cv routine ',i,j,k, CV_Face_value(i,j,k,:), value

	enddo
	enddo
	enddo

end subroutine get_CV_surface_contributions

end module get_timesteps_module

!subroutine control_volume_power_partialint(fij,ri,rj,vi_mhdt,vj_mhdt,ai_mdt,aj_mdt,ai,aj,molnoi,molnoj)
!    use module_record
!	use CV_objects, only : CV_debug,CV_constraint
!    use librarymod, only : heaviside => heaviside_a1, bubble_sort
!	use get_timesteps_module
!    implicit none

!	integer,intent(in)						 	:: molnoi, molnoj
!	double precision,dimension(3),intent(in) 	:: fij,ri,rj,vi_mhdt,vj_mhdt,ai_mdt,aj_mdt,ai,aj


!	integer										:: i,j,k,m,n,ixyz,ncrossings,ncrossingsi,ncrossingsj,count_t
!	integer,dimension(3)						:: cbin, ibin, jbin, ibin_mdt, jbin_mdt
!    double precision							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
!    double precision							:: invrij2,rij2_mdt,fijvi,fijvi_mdt,fijvi_trapz,delta_t_portion,eps = 0.0001d0
!	double precision,dimension(3)   			:: fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb,velvect
!	double precision,dimension(3)   			:: vi,vj,vi_mdt,vj_mdt
!	double precision,dimension(3)   			:: ri_mdt,rj_mdt,rij_mdt,fij_mdt,ri_p,rj_p
!	double precision,dimension(3)				:: Fbinsize, bintop, binbot
!	double precision,dimension(6)				:: crosstime
!	double precision,dimension(:),allocatable	:: delta_t_list,delta_t_array

!	character(30)	:: frmtstr

!	!Determine bin size
!	Fbinsize(:) = domain(:) / nbins(:)

!! 	!Old velocities
!! 	vi_mhdt = vi_hdt - ai_mdt*delta_t
!! 	vj_mhdt = vj_hdt - aj_mdt*delta_t
!! 	
!! 	!Old positions
!! 	ri_mdt = ri - vi_mhdt*delta_t
!! 	rj_mdt = rj - vj_mhdt*delta_t

!! 	!Get fij, vi and vj
!! 	vi = vi_hdt + 0.5*ai_mdt*delta_t		!This is not okay -- ai_pdt should be used...
!! 	vj = vi_hdt + 0.5*ai_mdt*delta_t		!This is not okay -- ai_pdt should be used...

!! 	!Get fij_mdt, vi_mdt and vj_mdt
!! 	rij_mdt(:)= ri_mdt(:) - rj_mdt(:)   					!Evaluate distance between particle i and j
!! 	invrij2  = 1.d0/(dot_product(rij_mdt,rij_mdt))					!Invert value
!! 	fij_mdt = 48.d0*(invrij2**7-0.5d0*invrij2**4)*rij_mdt
!! 	vi_mdt = vi_hdt - 0.5*ai_mdt*delta_t
!! 	vj_mdt = vi_hdt - 0.5*ai_mdt*delta_t


!	!Assuming acceleration passed in are at time t-dt/2

!	!Velocities for dot product @ t and t-dt
!	vi(:) 	  = vi_mhdt + 0.5d0*delta_t*ai(:)
!	vj(:) 	  = vj_mhdt + 0.5d0*delta_t*aj(:)

!	vi_mdt(:) = vi_mhdt - 0.5d0*delta_t*ai_mdt(:)
!	vj_mdt(:) = vj_mhdt - 0.5d0*delta_t*aj_mdt(:)

!	!Velocity=>Positions=>fij at t-dt
!	!vi_mhdt(:) = vi_hdt - delta_t*ai_mdt(:)
!	!vj_mhdt(:) = vj_hdt - delta_t*aj_mdt(:)

! 	ri_mdt(:) = ri - delta_t*vi_mhdt(:)
! 	rj_mdt(:) = rj - delta_t*vj_mhdt(:)

! 	rij_mdt(:)= ri_mdt(:) - rj_mdt(:)   					!Evaluate distance between particle i and j
!	rij2_mdt = dot_product(rij_mdt,rij_mdt)
!	if (rij2_mdt < rcutoff2) then
!	 	invrij2  = 1.d0/rij2_mdt
! 		fij_mdt = 48.d0*(invrij2**7-0.5d0*invrij2**4)*rij_mdt
!	else
!		fij_mdt = 0.d0
!	endif

!	!Trapizium rule calculation of power to add during portion of timestep
!	fijvi 	  = dot_product(fij    ,  vi  )
!	fijvi_mdt = dot_product(fij_mdt,vi_mdt)


!! 	vi_phdt(:) = vi_hdt + delta_t*ai(:)
!! 	vj_phdt(:) = vj_hdt + delta_t*aj(:)

!! 	ri_pdt(:) = ri + delta_t*vi_phdt(:)
!! 	rj_pdt(:) = rj + delta_t*vj_phdt(:)

! 	!Get fij_mdt, vi_mdt and vj_mdt
!!  	rij_pdt(:)= ri_pdt(:) - rj_pdt(:)   					!Evaluate distance between particle i and j
!!  	invrij2  = 1.d0/(dot_product(rij_pdt,rij_pdt))					!Invert value
!!  	fij_pdt = 48.d0*(invrij2**7-0.5d0*invrij2**4)*rij_pdt



!	!Assign to bins using integer division
!	ibin(:) 	= ceiling((ri(:)	+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!	jbin(:) 	= ceiling((rj(:)	+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!	ibin_mdt(:) = ceiling((ri_mdt(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish previous bin
!	jbin_mdt(:) = ceiling((rj_mdt(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish previous bin

!	!Get number of crossings
!	ncrossingsi =  abs(ibin(1) - ibin_mdt(1)) + abs(ibin(2) - ibin_mdt(2)) + abs(ibin(3) - ibin_mdt(3))
!	ncrossingsj =  abs(jbin(1) - jbin_mdt(1)) + abs(jbin(2) - jbin_mdt(2)) + abs(jbin(3) - jbin_mdt(3)) 
!	ncrossings = ncrossingsi+ncrossingsj

!	!Get portion of timestep in each cell
!	count_t = 0
!	if (ncrossings .eq. 0) then
!		allocate(delta_t_array(2))
!		delta_t_array(1) = 0.d0
!		delta_t_array(2) = delta_t
!	else
!		! Each molecular crossing time is counted twice -- once for the cell it's 
!		! leaving and once for the cell it's moving into
!		allocate(delta_t_list(2*ncrossings))
!		if (ncrossingsi .ne. 0) print'(a,2i4,6f12.8,6i4)', 'Molecule i',iter, molnoi, ri,ri_mdt,ibin,ibin_mdt
!		if (ncrossingsj .ne. 0) print'(a,2i4,6f12.8,6i4)', 'Molecule j',iter, molnoj, rj,rj_mdt,jbin,jbin_mdt

!    	!Check if molecule i has changed bin
!		call  get_timesteps(ncrossingsi,ibin,ibin_mdt,Fbinsize,ri,vi_mhdt,delta_t,delta_t_list,count_t)

!    	!Check if molecule j has changed bin
!		call  get_timesteps(ncrossingsj,jbin,jbin_mdt,Fbinsize,rj,vj_mhdt,delta_t,delta_t_list,count_t)

!		!Get crossing times in chronological order
!		call bubble_sort(delta_t_list)

!		!Sanity checks -- count_t == ncrossings & sum(delta_t_list) == delta_t
!		if (ncrossings .ne. 0.5*count_t) then
!			print'(a,2i12)', ' count_t == ncrossings ',count_t/2,ncrossings
!			stop "Error - crossing values not equal"
!		endif
!		if (any(delta_t_list .gt. delta_t)) then
!			stop "Error - delta t not ok in energy CV"
!		endif

!		!print'(i4,a,2i12,a,2f10.5)',iter,' count_t == ncrossings ',ncrossings,count_t/2, & 
!		!						 ' sum(delta_t_list) == delta_t ',delta_t, sum(delta_t_list)

!		!Set initial and final time and copy timestep arrays 
!		allocate(delta_t_array(ncrossings+2))
!		delta_t_array(1) = 0.d0
!		delta_t_array(ncrossings+2) = delta_t
!		n = 2
!		do i = 1,size(delta_t_list),2
!			if (delta_t_list(i) .ne. delta_t_list(i+1)) then
!				stop "Error - Successive timesteps not same in delta_t array"
!			endif 
!			delta_t_array(n) = delta_t_list(i)
!			n = n + 1
!		enddo

!	endif

!	if (ncrossings .ge. 1) then
! 	if (all(ibin .eq. 3)  .or. &
! 		all(jbin .eq. 3) 	) then
!!  	if (all(ibin .eq. 3) .and. any(jbin .ne. 3) .or. &
!!  		all(jbin .eq. 3) .and. any(ibin .ne. 3)	) then
!		print'(a,12i4)','Predicted cells', ibin,ibin_mdt,jbin,jbin_mdt
!	endif
!	endif

!	!Calculate surface contributions
!	do n=1,ncrossings+1
!	
!		delta_t_portion = delta_t_array(n+1)-delta_t_array(n)

!		!Calculate intermediate positions
!		ri_p = ri - vi_mhdt*delta_t_portion + eps
!		rj_p = rj - vj_mhdt*delta_t_portion + eps

!		ibin(:) 	= ceiling((ri_p(:)	+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!		jbin(:) 	= ceiling((rj_p(:)	+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin

!		if (ncrossings .ge. 1) then
! 		if (all(ibin .eq. 3) .and. any(jbin .ne. 3) .or. &
! 			all(jbin .eq. 3) .and. any(ibin .ne. 3)	) then
!				print'(a,12i4,3f14.9)','Delta_t',iter,n,molnoi,molnoj,ibin,jbin, & 
!												 ncrossingsi,ncrossingsj, delta_t_array(n),& 
!												 delta_t_portion,delta_t_array(n+1)
!		endif
!		endif

!		!Calculate and add fraction of interation
!		fijvi_trapz = 0.5d0*delta_t_portion * (fijvi + fijvi_mdt)/delta_t
!		call  get_CV_surface_contributions(ibin,jbin,ri_p,rj_p,Fbinsize,fijvi_trapz,Pxyvface2,molnoi,molnoj)

!!  		if (all(ibin .eq. 3) .and. any(jbin .ne. 3)) then
!! 			print'(a,3i4,10f11.6)', 'Saving ints',iter,molnoi,molnoj,delta_t_portion , & 
!! 													0.25*Pxyvface2(ibin(1),ibin(2),ibin(3),:)/5.0, &  !binsize
!! 													fijvi,fijvi_mdt,fijvi_trapz
!! 		endif
!! 		if (all(jbin .eq. 3) .and. any(ibin .ne. 3)) then
!! 			print'(a,3i4,10f11.6)', 'Saving ints',iter,molnoi,molnoj,delta_t_portion , &
!! 													0.25*Pxyvface2(jbin(1),jbin(2),jbin(3),:)/5.0, & 
!! 													fijvi,fijvi_mdt,fijvi_trapz
!! 		endif

!! 		if (all(ibin .eq. 3) .and. any(jbin .ne. 3) .or. &
!! 			all(jbin .eq. 3) .and. any(ibin .ne. 3)	) then
!! 			print'(a,3i4,10f11.7)', 'Current',iter,molnoi,molnoj,fij_mdt,vi_mdt,fijvi_mdt,0.25*fijvi_trapz/(5.0*delta_t)
!! 			print'(a,3i4,10f11.7)', 'Future ',iter,molnoi,molnoj,fij,vi,fijvi,0.25*fijvi_trapz/(5.0*delta_t)
!! 		endif
!		!print'(a,12f10.5)', 'Power after', Pxyvface(ibin(1),ibin(2),ibin(3),:), Pxyvface(jbin(1),jbin(2),jbin(3),:)

!	enddo

!end subroutine control_volume_power_partialint




! !===================================================================================
! !Forces over the surface of a Volume optmised for computational efficiency

! subroutine control_volume_stresses_opt(fij,ri,rj,molnoi)
! 	use module_record
!     use librarymod, only : heaviside  =>  heaviside_a1
! 	implicit none

! 	integer							:: i,j,k,ixyz,molnoi,molnoj
! 	integer,dimension(3)			:: cbin, ibin, jbin, Si
! 	double precision,dimension(3)	:: ri,rj,rij,fij,fsurface,Px,Py,Pz,sgnjit,sgnjib,onfaceb,onfacet,velvect
! 	double precision,dimension(3)	:: Fbinsize, bintop, binbot

! 	!Calculate rij
! 	rij = ri - rj
! 	!Prevent Division by zero
! 	do ixyz = 1,3
! 		if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
! 	enddo

! 	!Determine bin size
! 	Fbinsize(:) = domain(:) / nbins(:)

! 	!Assign to bins using integer division
! 	ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
! 	jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin

! 	if (ibin(1) .eq. jbin(1) .and. ibin(2) .eq. jbin(2) .and. ibin(3) .eq. jbin(3)) return
! 		
! 	do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
! 	do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
! 	do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

! 		cbin(1) = i; cbin(2) = j; cbin(3) = k

! 		bintop(:) = (cbin(:)-1*nhb(:)  )*Fbinsize(:)-halfdomain(:)
! 		binbot(:) = (cbin(:)-1*nhb(:)-1)*Fbinsize(:)-halfdomain(:)

! 		!Calculate the plane intersect of line with surfaces of the cube
! 		Px=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
! 		Py=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
! 		Pz=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)

! 		sgnjit(:)= sign(1.d0,bintop(:)- rj(:)) - sign(1.d0,bintop(:)- ri(:))
! 		sgnjib(:)= sign(1.d0,binbot(:)- rj(:)) - sign(1.d0,binbot(:)- ri(:))

! 		Si(1) =	(heaviside(bintop(2)-Px(2)) - heaviside(binbot(2)-Px(2)))* &
! 				(heaviside(bintop(3)-Px(3)) - heaviside(binbot(3)-Px(3)))
! 		Si(2) =	(heaviside(bintop(1)-Py(1)) - heaviside(binbot(1)-Py(1)))* &
! 				(heaviside(bintop(3)-Py(3)) - heaviside(binbot(3)-Py(3)))
! 		Si(3) =	(heaviside(bintop(1)-Pz(1)) - heaviside(binbot(1)-Pz(1)))* &
! 				(heaviside(bintop(2)-Pz(2)) - heaviside(binbot(2)-Pz(2)))

! 		onfaceb(:) = sgnjib(:)*dble(Si(:))
! 		onfacet(:) = sgnjit(:)*dble(Si(:))

! 		!Stress acting on face over volume
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,1) = Pxyface(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*onfaceb(1)
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,2) = Pxyface(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*onfaceb(2)
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,3) = Pxyface(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*onfaceb(3)
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,4) = Pxyface(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*onfacet(1)
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,5) = Pxyface(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*onfacet(2)
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,6) = Pxyface(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*onfacet(3)

! 		!Stress acting on face over volume
! 		if (eflux_outflag .ne. 0) then
! 			velvect(:) = v(:,molnoi) 
! 			!if (molnoi .gt. np) print*, velvect(1)
! 			!velvect(:) = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),1) = Pxyvface(cbin(1),cbin(2),cbin(3),1) + dot_product(fij,velvect)*onfaceb(1)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),2) = Pxyvface(cbin(1),cbin(2),cbin(3),2) + dot_product(fij,velvect)*onfaceb(2)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),3) = Pxyvface(cbin(1),cbin(2),cbin(3),3) + dot_product(fij,velvect)*onfaceb(3)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),4) = Pxyvface(cbin(1),cbin(2),cbin(3),4) + dot_product(fij,velvect)*onfacet(1)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),5) = Pxyvface(cbin(1),cbin(2),cbin(3),5) + dot_product(fij,velvect)*onfacet(2)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),6) = Pxyvface(cbin(1),cbin(2),cbin(3),6) + dot_product(fij,velvect)*onfacet(3)
! 		endif

! 		!Force applied to volume
! 		fsurface(:) = 0.d0
! 		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*(onfaceb(1) - onfacet(1))
! 		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*(onfaceb(2) - onfacet(2))
! 		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*(onfaceb(3) - onfacet(3))
! 		volume_force(cbin(1),cbin(2),cbin(3),:,1) = volume_force(cbin(1),cbin(2),cbin(3),:,1) + fsurface*delta_t

! 	enddo
! 	enddo
! 	enddo

! end subroutine control_volume_stresses_opt


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
    use librarymod, only : heaviside  =>  heaviside_a1
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

module module_record_external_forces

contains

subroutine record_external_forces(F,ri,vi)
	use module_record, only : domain,halfdomain, nbins, nhb
	use calculated_properties_MD, only :  F_ext_bin,Fv_ext_bin
	use computational_constants_MD, only : eflux_outflag, delta_t,iter
	use librarymod, only : imaxloc
	implicit none

	double precision,dimension(3),intent(in):: F,ri
	double precision,dimension(3),intent(in),optional :: vi

	integer									:: jxyz
	integer	,dimension(3)					:: ibin, ibin_pt, crossface
	double precision						:: crosstimetop,crosstimebot,delta_t_cross, Fiextvi
	double precision,dimension(3)			:: mbinsize, ri_pt, bintop, binbot

	!Determine bin size and bin
	mbinsize(:) = domain(:) / nbins(:)
	ibin(:) = ceiling((ri(:)+halfdomain(:))/mbinsize(:)) + nhb(:)

	!Add external force to bin
	F_ext_bin(ibin(1),ibin(2),ibin(3),:) = & 
		F_ext_bin(ibin(1),ibin(2),ibin(3),:) + F(:)

	if (eflux_outflag .eq. 4) then
		if ( present(vi)) then
			Fiextvi = dot_product(F(:),vi(:))
			Fv_ext_bin(ibin(1),ibin(2),ibin(3)) = & 
				Fv_ext_bin(ibin(1),ibin(2),ibin(3)) + Fiextvi
		else
			!Velocity assumed to be zero if not supplied
		endif
	endif

end subroutine record_external_forces

end module module_record_external_forces


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
!		oldj  =>  cell%head(icell+icellshift,jcell+jcellshift)%point
!		adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift)
!
!		do j = 1,adjacentcellnp          !Step through all j for each i
!			rij2=0                   !Set rij^2 to zero
!			molnoj = oldj%molno !Number of molecule
!			rj = r(:,molnoj)         !Retrieve rj
!					
!			currentj  =>  oldj
!			oldj  =>  currentj%next    !Use pointer in datatype to obtain next item in list

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
