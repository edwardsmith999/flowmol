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

	real(kind(0.d0)) :: vel

contains

    function get_bin(r) result(bin)
        use computational_constants_MD, only : halfdomain, nhb
        use calculated_properties_MD, only : binsize
        implicit none

        real(kind(0.d0)),intent(in),dimension(3) :: r
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

    integer, save   :: vmd_skip_count=0
    integer, save   :: vmdintervalno = 1

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
	if (vmd_outflag.ne.0 .and. size(vmd_intervals,2).ge.vmdintervalno) then
        vmd_skip_count = vmd_skip_count + 1 
        if (vmd_skip_count .eq. vmd_skip) then
            vmd_skip_count = 0
            call parallel_io_write_vmd(vmdintervalno,vmd_count)
        endif
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
	if (momentum_outflag .ne. 0) call momentum_averaging(momentum_outflag)

	!Obtain and record mass only
	if (momentum_outflag .eq. 0 .and. &
		mass_outflag .ne. 0) call mass_averaging(mass_outflag)

    !Cluster analysis or average bin based tracking of liquid vapour interfaces
    if (cluster_analysis_outflag .eq. 1) then
        call get_interface_from_clusters()
    elseif (cluster_analysis_outflag .eq. 2) then
        call sl_interface_from_binaverage()
    endif

	!Obtain and record velocity distributions
	if (vPDF_flag .ne. 0) call velocity_PDF_averaging(vPDF_flag)

	!Record boundary force PDFs
	if (bforce_pdf_measure .ne. 0) call bforce_pdf_stats

	!Obtain and record temperature
	if (temperature_outflag .ne. 0)	call temperature_averaging(temperature_outflag)

	!Obtain and record energy
	if (energy_outflag .ne. 0)	call energy_averaging(energy_outflag)

	!Obtain and record density on a surface
    if (msurf_outflag .ne. 0) call surface_density_averaging(msurf_outflag)

	!Obtain and record molecular diffusion
	!call evaluate_properties_diffusion

	!Calculate pressure tensor
	if (pressure_outflag .ne. 0) then
		call pressure_averaging(pressure_outflag)
	end if

    !Calculate heat flux information
	if (heatflux_outflag .ne. 0) then
		call heatflux_averaging(heatflux_outflag)
	end if

	!Write current timestep to a progress file
	call update_simulation_progress_file()

end subroutine simulation_record

!==========================================================================
!Calculate kinetic and potential energy as well as temperature and pressure

subroutine evaluate_macroscopic_properties
	use module_record
	use messenger, only : globalise
    use messenger_data_exchange, only : globalSum
    use module_set_parameters, only : mass
	implicit none

	integer :: n,ixyz
    real(kind(0.d0))    :: msum

    msum  = 0.d0
	vsum  = 0.d0                                                ! Reset all sums
	mv2sum = 0.d0                                                ! Reset all sums

	! Potential Component
	select case(potential_flag)
	case(0)
		potenergysum = sum(potenergymol(1:np))
		call globalSum(potenergysum)
	case(1)
		potenergysum_LJ = sum(potenergymol_LJ(1:np))
		potenergysum_POLY = sum(potenergymol_POLY(1:np))
		potenergysum = sum(potenergymol_LJ(1:np) + potenergymol_POLY(1:np))
		call globalSum(potenergysum_LJ)
		call globalSum(potenergysum_POLY)
		call globalSum(potenergysum)
	case default
		call error_abort("Unrecognised potential flag in simulation_record")
	end select

	! Kinetic Component
	select case(integration_algorithm)
	case(leap_frog_verlet)
		do n = 1, np 									! Loop over all particles
        msum = msum + mass(n)
		do ixyz = 1, nd									! Loop over all dimensions
			vel   = v(ixyz,n) + 0.5d0*a(ixyz,n)*delta_t	! Velocity must shifted half a timestep
			vsum  = vsum + mass(n)*vel							! Add up all molecules' velocity components
			mv2sum = mv2sum + mass(n)*vel**2			! Add up all molecules' velocity squared components  
		enddo
		enddo
	case(velocity_verlet) 								! If velocity Verlet algorithm
		do n = 1, np
        msum = msum + mass(n)
		do ixyz = 1, nd
			vel   = v(ixyz,n)
			vsum  = vsum + mass(n)*vel
			mv2sum = mv2sum + mass(n)*vel**2          	! Sum all velocity squared components
		enddo
		enddo
	end select
        
	!Obtain global sums for all parameters
	call globalSum(msum)
	call globalSum(vsum)
	call globalSum(mv2sum)
	virial = sum(virialmol(1:np))
	call globalSum(virial)

	kinenergy   = (0.5d0 * mv2sum) / real(globalnp,kind(0.d0)) 
	potenergy   = (0.5d0 * potenergysum) / real(globalnp,kind(0.d0)) + Potential_sLRC !N.B. extra 1/2 as all interactions calculated
	if (potential_flag.eq.1) then
		potenergy_LJ= (0.5d0 * potenergysum_LJ)/real(globalnp,kind(0.d0)) + Potential_sLRC
		potenergy_POLY= (0.5d0 * potenergysum_POLY)/real(globalnp,kind(0.d0))
	end if

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
	temperature = mv2sum / real(nd*globalnp,kind(0.d0))
	if (any(periodic.gt.1)) temperature = get_temperature_PUT()
	pressure    = (density/real(globalnp*nd,kind(0.d0)))*(mv2sum+virial/2.d0) + Pressure_sLRC !N.B. virial/2 as all interactions calculated

!	kinenergy   = (0.5d0 * mv2sum) / msum
!	potenergy   = potenergysum /(2.d0*msum) + Potential_sLRC !N.B. extra 1/2 as all interactions calculated
!	if (potential_flag.eq.1) then
!		potenergy_LJ= potenergysum_LJ/(2.d0*msum) + Potential_sLRC
!		potenergy_POLY= potenergysum_POLY/(2.d0*msum)
!	end if
!	!print'(4(a,f18.8))', ' <PE>= ',potenergy, & 
!	!						 ' std(PE) = ',sqrt(sum((potenergymol(1:np)-potenergy)**2)/(2.d0*real(msum,kind(0.d0)))), & 
!	!						 ' max= ',maxval(potenergymol(1:np)),' min= ',minval(potenergymol(1:np))
!	if (maxval(potenergymol(1:np)) .gt. 10000) then
!		print*, 'np = ', np
!		print*, 'max(potenergymol) = ', maxval(potenergymol(1:np))
!		do n=1,np
!			write(3000+irank,'(i10,f28.4,6f10.4)'), n , potenergymol(n), r(:,n), globalise(r(:,n))
!		enddo
!		print*, 'Simulation aborted because max PE has reached an unreasonably high value.'
!		print*, 'Inspect fort.(3000+irank) for n, potenergymol, r, r_global'
!		call error_abort("STOPPING CODE")
!	endif
!	totenergy   = kinenergy + potenergy
!	temperature = mv2sum / (real(nd,kind(0.d0))*msum)
!	if (any(periodic.gt.1)) temperature = get_temperature_PUT()
!	pressure    = (density/(real(nd,kind(0.d0))*msum))*(mv2sum+virial/2.d0) + Pressure_sLRC !N.B. virial/2 as all interactions calculated

    !print'(a,i8,3f20.10)', 'pressure   ', iter, mv2sum+virial/2.d0, mv2sum , virial

end subroutine evaluate_macroscopic_properties

subroutine evaluate_microstate_pressure
	use module_record
    use messenger_data_exchange, only : globalSum
    use module_set_parameters, only : mass
	implicit none

	integer :: n
	real(kind(0.d0)) :: vtemp(nd)
    real(kind(0.d0)) :: msum

    msum  = 0.d0
	mv2sum = 0.d0

	! Kinetic part of Virial
	select case(integration_algorithm)
	case(leap_frog_verlet)

		do n = 1, np
            msum = msum + mass(n)
			vtemp(:) = v(:,n) + 0.5d0*a(:,n)*delta_t ! Velocity must shifted half a timestep
			mv2sum = mv2sum + mass(n) * dot_product(vtemp,vtemp)
		enddo

	case(velocity_verlet)

		do n = 1, np
            msum = msum + mass(n)
			mv2sum = mv2sum + mass(n) * dot_product(v(:,n),v(:,n))
		enddo

	end select

	! Configurational part of Virial
	virial = sum(virialmol(1:np))
	call globalSum(virial)

	! Instantaneous pressure 
	pressure = (density/(globalnp*nd))*(mv2sum+virial/2.d0) + Pressure_sLRC !N.B. virial/2 as all interactions calculated

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
				it,';', simtime,';',vsum,';', mv2sum,';', temperature,';', &
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
				it,';',simtime,';',vsum,';', mv2sum,';', temperature,';', &
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

	real(kind(0.d0)) 	                        :: Hfunction,vfactor,const
	real(kind(0.d0)),save 	                    :: meanstream, meannp
	real(kind(0.d0)),dimension(3)	            :: peculiarv,binsize_
	real(kind(0.d0)),dimension(:),allocatable 	:: vmagnitude,binloc

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
		    !allocate(binloc,source=velPDF_array(nhb(1)+1,nhb(2)+1,nhb(3)+1,1)%binvalues())
		    allocate(binloc(velPDF_array(nhb(1)+1,nhb(2)+1,nhb(3)+1,1)%nbins))
            binloc = velPDF_array(nhb(1)+1,nhb(2)+1,nhb(3)+1,1)%binvalues()

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
!	real(kind(0.d0)) :: Hfunction,vfactor,const,streamvel
!	real(kind(0.d0)),save :: meanstream, meannp
!	real(kind(0.d0)),dimension(3)		:: peculiarv
!	real(kind(0.d0)),dimension(:),allocatable :: vmagnitude,normalisedvfd_bin,normalisedvfdMB_bin,binloc

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
!	real(kind(0.d0)) 	:: Hfunction,vfactor,const
!	real(kind(0.d0)),save 	:: meanstream, meannp
!	real(kind(0.d0)),dimension(3)	:: peculiarv,binsize_
!	real(kind(0.d0)),dimension(:),allocatable 	:: vmagnitude,normalisedvfd_bin,normalisedvfdMB_bin,binloc
!	real(kind(0.d0)),dimension(:,:,:),allocatable 	:: streamvel

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
	real(kind(0.d0))                           :: rmag,dr,Nideal,Nbin,dV
	real(kind(0.d0)), dimension(nd)            :: rij

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
	real(kind(0.d0)) :: x,dx,Nideal,Nbin,dV
	real(kind(0.d0)) :: xij

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
	real(kind(0.d0))   :: c,s
	real(kind(0.d0)), dimension(nd)  :: k,dk

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

	real(kind(0.d0)) :: eta              ! Shear viscosity
	real(kind(0.d0)) :: N1               ! First normal stress difference
	real(kind(0.d0)) :: N2               ! Second normal stress difference
	real(kind(0.d0)) :: P                ! Hydrostatic pressure
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
	integer :: chainID, i_sub, j_sub, funcy
	real(kind(0.d0)) :: etev_prod, etev_prod_sum
	real(kind(0.d0)) :: etev2, etev2_sum
	real(kind(0.d0)), dimension(nd) :: rij

	if (iter.eq.etevtcf_iter0) then	!Initialise end-to-end vectors at t_0

		allocate(etev(nchains,nd))
		allocate(etev_0(nchains,nd))
		etev_0 = 0.d0
		do i=1,np
			chainID = monomer(i)%chainID
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
				etev_0(chainID,:) = etev_0(chainID,:) + rij(:)
			end do
		end do
		call globalSum(etev_0,nchains,nd)

	end if

	!--------------------------------------------------------------------
	!Calculate all end-to-end vectors for file output
	etev = 0.d0
	do i=1,np
		chainID = monomer(i)%chainID
		i_sub = monomer(i)%subchainID
		funcy = monomer(i)%funcy
		do nbond=1,funcy
			j             = bond(nbond,i)          ! Molecule number j is nth bond to i
			j_sub         = monomer(j)%subchainID  ! Find subchain ID of mol j
			if (j_sub.lt.i_sub) cycle              ! Avoid counting backwards
			rij(:)        = r(:,j) - r(:,i)                     
			rij(:)        = rij(:) - &
                            globaldomain(:)*anint(rij(:)/globaldomain(:))
			etev(chainID,:) = etev(chainID,:) + rij(:)
		end do
	end do
	call globalSum(etev,nchains,nd)
	!--------------------------------------------------------------------

	!--------------------------------------------------------------------
	!Running calculation for stdout...	
	etev_prod_sum	= 0.d0
	etev2_sum 		= 0.d0
	do chainID=1,nchains
		etev_prod		  = dot_product(etev(chainID,:),etev_0(chainID,:))
		etev2			  = dot_product(etev(chainID,:),etev(chainID,:))			
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
	real(kind(0.d0)) :: R_g2
	real(kind(0.d0)), dimension(nd) :: rij
	real(kind(0.d0)), dimension(:,:), allocatable :: r_cm

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
    use module_set_parameters, only : mass
	implicit none

	integer         				:: n, ixyz
	integer         				:: cbin
	integer,dimension(3)			:: ibin
	real(kind(0.d0))				:: slicebinsize
	real(kind(0.d0)),dimension(3) 	:: mbinsize 

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
			if (cbin < 1 ) cbin = 1        				                !Prevents out of range values
			slice_mass(cbin)= slice_mass(cbin) + mass(n)      			 !Add mass to current bin
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
                volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + mass(n)
            enddo
        case(1:2)
            !Seperate logging for polymer and non-polymer 
            do n = 1,np
                !Add up current volume mass densities
                !ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb
                ibin = get_bin(r(:,n))
                if (monomer(n)%chainID .eq. 0) then
                    !Skip wall molecules
                    if (split_pol_sol_stats .eq. 1 .and. &
                        any(tag(n).eq.tether_tags)) cycle
                    volume_mass_s(ibin(1),ibin(2),ibin(3)) = &
                    volume_mass_s(ibin(1),ibin(2),ibin(3)) + mass(n)
                else
                    volume_mass_p(ibin(1),ibin(2),ibin(3)) = &
                    volume_mass_p(ibin(1),ibin(2),ibin(3)) + mass(n)
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
!		RECORD MOMENTUM AT LOCATION IN SPACE
! Either by binning or taking slices on a plane, the momentum field in the molecular
! system is recorded and output
!===================================================================================

!Calculate averaged momentum components of each bin or slice with 2D slice in 
!ixyz = 1,2,3 or in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine momentum_averaging(ixyz)
	use module_record
	use field_io , only : momentum_slice_io,momentum_bin_io,momentum_bin_cpol_io
	use concentric_cylinders, only: cyl_mass, cyl_mom
	use linked_list
	implicit none

	integer				:: ixyz, n
	integer,dimension(3):: ib
	integer, save		:: average_count=-1
	real(kind(0.d0)),dimension(3) 	:: Vbinsize

	average_count = average_count + 1
	call cumulative_momentum(ixyz)

	!Save streaming momentum for temperature averages
	if (temperature_outflag .ne. 0 .and. peculiar_flag .ne. 0) then

		! NOTE THE peculiar_flag CAN BE SET TO 0 AND
		! THE UNBIAS TEMPERATURE USUALLY CALCULATED FROM
		! T_{unbias} = (1/3N) * \sum_i^N m_i * (vi - u)*(vi - u)
		! CAN INSTEAD BE CALCULATED FROM:
		! T_{unbias} = (1/3N) * \sum_i^N m_i * vi*vi - u^2/3
		call error_abort("Peculiar momentum functionality removed -- "//&
                         "please calculate using T_{unbias} = (1/3N) "//&
                         "* \sum_i^N m_i*vi*vi - u^2/3")

		do n=1,np
			!Save streaming momentum per molecule
			!ib(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
            ib(:) = get_bin(r(:,n))
			U(:,n) =  volume_momentum(ib(1),ib(2),ib(3),:) / volume_mass(ib(1),ib(2),ib(3))
		enddo

	endif

	if (average_count .eq. Nvel_ave) then
		average_count = 0

		select case(ixyz)
			case(1:3)
				call momentum_slice_io(ixyz)
				!Reset momentum slice
				slice_mass = 0
				slice_momentum  = 0.d0
			case(4)

                select case(split_pol_sol_stats)
                case(0)
                    call momentum_bin_io(volume_mass, volume_momentum,'bins')
                    volume_mass = 0
                    volume_momentum = 0.d0
                case(1:2)
                    call momentum_bin_io(volume_mass_s, volume_momentum_s, 'solv')
                    call momentum_bin_io(volume_mass_p, volume_momentum_p, 'poly')
                    call momentum_bin_io(volume_mass, volume_momentum, 'bins')
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
				call momentum_bin_cpol_io(cyl_mass,cyl_mom)
				cyl_mass = 0
				cyl_mom = 0.d0
			case default
				call error_abort("Error input for momentum averaging incorrect")
			end select

			!Collect velocities for next step
			!call cumulative_momentum(ixyz)

	endif

end subroutine momentum_averaging

!-----------------------------------------------------------------------------------
!Add velocities to running total, with 2D slice in ixyz = 1,2,3 or
!in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine cumulative_momentum(ixyz)
	use module_record
	use concentric_cylinders
	use messenger, only: globalise, localise
	use linked_list
    use librarymod, only : cpolariser, cpolarisev
    use module_set_parameters, only : mass
	implicit none

	integer							:: n,ixyz
	integer         				:: cbin
	integer		,dimension(3)		:: ibin
	real(kind(0.d0))				:: slicebinsize
	real(kind(0.d0)),dimension(3) 	:: Vbinsize 
	
	integer :: br, bt, bz
	real(kind(0.d0)) :: fluiddomain_cyl(3), rglob(3), rpol(3), vpol(3)

	!In case someone wants to record momentum in a simulation without sliding walls!?!?
	if (ensemble .ne. tag_move) then
		allocate(slidev(3,np)); slidev = 0.d0
	endif

	select case(ixyz)
	!momentum measurement is a number of 2D slices through the domain
	case(1:3)

		slicebinsize = domain(ixyz) / nbins(ixyz)

		do n = 1, np    ! Loop over all particles
			!Assign to bins using integer division
			cbin = ceiling((r(ixyz,n)+halfdomain(ixyz))/slicebinsize)!Establish current bin
			if (cbin > nbins(ixyz)) cbin = nbins(ixyz) 		 !Prevents out of range values
			if (cbin < 1 ) cbin = 1        				 !Prevents out of range values
			slice_mass(cbin)= slice_mass(cbin) + mass(n)      			 !Add one to current bin
			slice_momentum(cbin,:) = slice_momentum(cbin,:) + mass(n) * (v(:,n) + slidev(:,n)) 	 !Add streamwise momentum to current bin
		enddo

	!momentum measurement for 3D bins throughout the domain
	case(4)

        Vbinsize(:) = domain(:) / nbins(:) 
        select case (split_pol_sol_stats)
        case(0)
            !Reset Control Volume momentum 
            do n = 1,np
                !Add up current volume mass and momentum densities
                !ibin(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
                ibin(:) = get_bin(r(:,n))
                volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + mass(n)
                volume_momentum(ibin(1),ibin(2),ibin(3),:) = volume_momentum(ibin(1),ibin(2),ibin(3),:) & 
                                                            + mass(n)*(v(:,n) + slidev(:,n))
            enddo
        case(1)
            !Reset Control Volume momentum 
            do n = 1,np
                !Add up current volume mass and momentum densities
                !ibin(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
                ibin(:) = get_bin(r(:,n))
                if (monomer(n)%chainID .eq. 0) then
                    if (any(tag(n).eq.tether_tags)) cycle
                    volume_mass_s(ibin(1),ibin(2),ibin(3)) = volume_mass_s(ibin(1),ibin(2),ibin(3)) + mass(n)
                    volume_momentum_s(ibin(1),ibin(2),ibin(3),:) = volume_momentum_s(ibin(1),ibin(2),ibin(3),:) & 
                                                                + mass(n)*(v(:,n) + slidev(:,n))
                else
                    volume_mass_p(ibin(1),ibin(2),ibin(3)) = volume_mass_p(ibin(1),ibin(2),ibin(3)) + mass(n)
                    volume_momentum_p(ibin(1),ibin(2),ibin(3),:) = volume_momentum_p(ibin(1),ibin(2),ibin(3),:) & 
                                                                + mass(n)*(v(:,n) + slidev(:,n))
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
		call error_abort("momentum Binning Error")
	end select

	if (ensemble .ne. tag_move) then
		deallocate(slidev)
	endif

	 
end subroutine cumulative_momentum

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
			if (momentum_outflag .ne. ixyz) slice_mass = 0
			slice_temperature  = 0.d0
		case(4)
			call temperature_bin_io(volume_mass,volume_temperature,'bins')
			!Reset temperature bins
			if (momentum_outflag .ne. 4) volume_mass = 0
			volume_temperature = 0.d0
		case(5)
			call temperature_bin_cpol_io(cyl_mass,cyl_KE)
			if (momentum_outflag .ne. 5) cyl_mass = 0
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
    use module_set_parameters, only : mass
	implicit none

	integer							:: n,ixyz
	integer         				:: cbin
	integer		,dimension(3)		:: ibin
	real(kind(0.d0))				:: slicebinsize
	real(kind(0.d0)),dimension(3) 	:: Tbinsize 

	integer :: br, bt, bz
	real(kind(0.d0)) :: fluiddomain_cyl(3), rglob(3), rpol(3)

	!In case someone wants to record momentum in a simulation without sliding walls!?!?
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
			if (momentum_outflag .ne. ixyz) & 
			slice_mass(cbin)= slice_mass(cbin) + mass(n)     			 !Add mass to current bin
			slice_temperature(cbin) = slice_temperature(cbin) & 
					+ mass(n)*dot_product((v(:,n)+slidev(:,n)),(v(:,n)+slidev(:,n))) 	 !Add streamwise temperature to current bin
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
			if (momentum_outflag .ne. 4) & 
			volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + mass(n)
			!Note - the streaming term is removed but includes sliding so this must be added back on
			if (peculiar_flag .eq. 0) then
				!if (mod(iter,1000) .eq. 0 	.and. &
				!	ibin(1) .eq. 4 			.and. &
				!	ibin(2) .eq. 2 			.and. &
				!	ibin(3) .eq. 4				) then
				!	print*, iter, ibin, dot_product(v(:,n),v(:,n)), v(:,n)
				!endif
				volume_temperature(ibin(1),ibin(2),ibin(3)) = volume_temperature(ibin(1),ibin(2),ibin(3)) & 
										+ mass(n)*dot_product((v(:,n)+slidev(:,n)),(v(:,n)+slidev(:,n)))
			else
				volume_temperature(ibin(1),ibin(2),ibin(3)) = volume_temperature(ibin(1),ibin(2),ibin(3)) & 
										+ mass(n)*dot_product((v(:,n)-U(:,n)+slidev(:,n)), & 
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

			if (momentum_outflag .ne. 5) then
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
    use module_set_parameters, only : mass
	implicit none

	integer							:: n,ixyz
	integer         				:: cbin
	integer		,dimension(3)		:: ibin
	real(kind(0.d0))				:: slicebinsize, energy
	real(kind(0.d0)),dimension(3) 	:: Tbinsize, velvect

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
		        energy = 0.5d0*(mass(n)*dot_product(velvect,velvect)+potenergymol(n))

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
!            print'(a,i8,3f20.10)', 'Pressure VA', iter, &
!                                    sum(Pxybin(:,:,:,1,1)+Pxybin(:,:,:,2,2)+Pxybin(:,:,:,3,3)), &
!                                    sum(vvbin(:,:,:,1,1)+vvbin(:,:,:,2,2)+vvbin(:,:,:,3,3)), &
!                                    sum(rfbin(1+nhb(1):nbins(1)+nhb(1),   & 
!                                              1+nhb(2):nbins(2)+nhb(2),   & 
!                                              1+nhb(3):nbins(3)+nhb(3),1,1) & 
!                                       +rfbin(1+nhb(1):nbins(1)+nhb(1),   & 
!                                              1+nhb(2):nbins(2)+nhb(2),   & 
!                                              1+nhb(3):nbins(3)+nhb(3),2,2) & 
!                                       +rfbin(1+nhb(1):nbins(1)+nhb(1),   & 
!                                              1+nhb(2):nbins(2)+nhb(2),   & 
!                                              1+nhb(3):nbins(3)+nhb(3),3,3))
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
	real(kind(0.d0)), dimension(3)		:: velvect
	real(kind(0.d0)), dimension(3)      :: rglob
	real(kind(0.d0)), dimension(3,3)	:: Pxytemp


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
		!Calculate mass [velocity (x) velocity] for kinetic part of stress tensor
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
    use module_set_parameters, only : mass
	implicit none

	integer								:: imin, jmin, kmin, imax, jmax, kmax
	integer         					:: n, ixyz,jxyz
	integer 							:: ibin, jbin, kbin
	integer 							:: bin(3)
	real(kind(0.d0)), dimension(3)		:: VAbinsize, velvect

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
					       		  + mass(n) * velvect(ixyz) * velvect(jxyz)
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
    use module_set_parameters, only : mass
	implicit none

	integer,intent(in)				:: imin, jmin, kmin, imax, jmax, kmax
	integer                         :: i, ixyz, jxyz   !Define dummy index
	integer 						:: ibin, jbin, kbin,icell, jcell, kcell
	integer                         :: cellnp, molnoi
	real(kind(0.d0)), dimension(3)	:: velvect
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
								  + mass(i) * velvect(ixyz) * velvect(jxyz)
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
	real(kind(0.d0)),dimension(3), intent(in)		        :: ri, rj, domain
	real(kind(0.d0)),dimension(:,:,:,:,:), intent(inout)	:: rfbin

    integer                                                 :: VA_line_samples_
	real(kind(0.d0)),dimension(3)                           :: halfdomain

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
	real(kind(0.d0)),dimension(3), intent(in)		:: ri, rj

	integer											:: ixyz, jxyz
	integer											:: diff
	integer,dimension(3)							:: ibin, jbin, bindiff 
	real(kind(0.d0)),dimension(3)					:: VAbinsize, rij

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
	real(kind(0.d0)),dimension(3), intent(in)		:: ri, rj


	integer											:: ixyz, i,j,k,l,n
	integer											:: diff
	integer,dimension(3)							:: ibin, jbin, bindiff 
	integer,dimension(:,:)		   ,allocatable		:: interbin, interbindiff
	real(kind(0.d0)),dimension(3)					:: VAbinsize, normal, p,temp1,temp2,temp3
	real(kind(0.d0)),dimension(3)                   :: rij
	real(kind(0.d0)),dimension(:,:)	   ,allocatable	:: intersection
	real(kind(0.d0)),dimension(:,:,:)  ,allocatable	:: MLfrac !Magnitude of fraction of stress
	real(kind(0.d0)),dimension(:,:,:,:),allocatable	:: Lfrac  !Fraction of stress in a given cell

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
	use librarymod, only: get_outerprod
	implicit none

	!integer,intent(in)  			:: imin, jmin, kmin, imax, jmax, kmax

	integer                         :: i, j, ixyz !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp 
	integer							:: molnoi, molnoj
	integer							:: icellmin,jcellmin,kcellmin,icellmax,jcellmax,kcellmax
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	real(kind(0.d0)),dimension(3)	:: vi_t, cellsperbin
	!real(kind(0.d0)), dimension(3,3):: rF
	real(kind(0.d0)), dimension(:,:), allocatable :: rF
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
						    accijmag = get_accijmag(invrij2, molnoi, molnoj) !48.d0*(invrij2**7-0.5d0*invrij2**4)
                            call get_outerprod(rij,rij*accijmag,rf)
                            !rf = outerprod(rij, accijmag*rij)

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
                            !rf = outerprod(rij, accijmag*rij)
                            call get_outerprod(rij,rij*accijmag,rf)
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
        call add_POLY_contribution
    end if

	nullify(oldi)      	!Nullify as no longer required
	nullify(oldj)      	!Nullify as no longer required
	nullify(currenti)      	!Nullify as no longer required
	nullify(currentj)      	!Nullify as no longer required

contains

    subroutine add_POLY_contribution
        use polymer_info_MD
        use Volume_average_pressure, only : pressure_tensor_forces_VA
	    !use librarymod, only: outerprod
	    use librarymod, only: get_outerprod
        use module_set_parameters, only : get_poly_accijmag
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
                !accijmag = -k_c/(1-(rij2/(R_0**2)))	!(-dU/dr)*(1/|r|)
                accijmag = -get_poly_accijmag(rij2, molnoi, molnoj)
                call get_outerprod(rij,rij*accijmag,rf)
                !rf = outerprod(rij,rij*accijmag)
                call pressure_tensor_forces_VA(ri, rj, rf, domain,  & 
                                                nbins, nhb, rfbin,  &
                                                VA_calcmethod=1,    & 
                                                VA_line_samples=VA_line_samples)

            end do	

        end do

    end subroutine add_POLY_contribution

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

	real(kind(0.d0)),dimension(3)	:: cellsperbin

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
						accijmag = get_accijmag(invrij2, molnoi, molnoj) !48.d0*(invrij2**7-0.5d0*invrij2**4)

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
        use module_set_parameters, only : get_poly_accijmag
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
                accijmag =  -get_poly_accijmag(rij2, molnoi, molnoj) !-k_c/(1-(rij2/(R_0**2)))	!(-dU/dr)*(1/|r|)

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
	real(kind(0.d0)), dimension(3)		:: velvect
	real(kind(0.d0)), dimension(3)      :: rglob
	real(kind(0.d0)), dimension(3,3)	:: Pxytemp

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

		!Calculate mass velocity dot [velocity (x) velocity] for kinetic part of heatflux tensor
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
    use module_set_parameters, only : mass
	implicit none

	integer								:: imin, jmin, kmin, imax, jmax, kmax
	integer         					:: n, ixyz,jxyz
	integer 							:: ibin, jbin, kbin
	integer 							:: bin(3)
    real(kind(0.d0))                    :: energy
	real(kind(0.d0)), dimension(3)		:: VAbinsize, velvect

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

        energy = 0.5d0 * (mass(n)*dot_product(velvect,velvect) + potenergymol(n))
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
	use module_record, only : Nmflux_ave, domain, nbins, nhb, & 
                              thermstattop, thermstatbottom, &
                              iter, irank, mass_flux
	use CV_objects, only : CV_debug, CVcheck_mass, check_error_mass !,CV_sphere_mass
	implicit none

	integer			                :: flag
    integer,dimension(3)            :: thermbinstop,thermbinsbot
    real(kind(0.d0)),dimension(3)   :: mbinsize
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
            
            !E.S. this causes a compiler seg fault for 
            !     ifort version 13.0.1 which is fixed by 
            !     replacing 
            !     "call CVcheck_mass%check_error( ... "   
            !     with 
            !     "check_error_mass(CVcheck_mass, ... "
            call check_error_mass(CVcheck_mass, &

             1+nhb(1)+thermbinsbot(1),nbins(1)+nhb(1)-thermbinstop(1), & 
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

!!===================================================================================
!! Mass Flux over a surface of a bin
!! Includes all intermediate bins

subroutine cumulative_mass_flux
	use module_record
    use librarymod, only : imaxloc, heaviside  =>  heaviside_a1
    use module_set_parameters, only : mass
    !use CV_objects, only : CV_sphere_mass
    implicit none

	integer							:: jxyz,i,j,k,n
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	real(kind(0.d0)),dimension(3)	:: mbinsize,crossface
	real(kind(0.d0)),dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

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
				      + mass(n)*nint(dble(onfacexb)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),2) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),2) & 
				      + mass(n)*nint(dble(onfaceyb)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),3) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),3) &
				      + mass(n)*nint(dble(onfacezb)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),4) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),4) &
				      + mass(n)*nint(dble(onfacext)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),5) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),5) &
				      + mass(n)*nint(dble(onfaceyt)*abs(crossface(jxyz)))
				mass_flux(cbin(1),cbin(2),cbin(3),6) = & 
					mass_flux(cbin(1),cbin(2),cbin(3),6) &
				      + mass(n)*nint(dble(onfacezt)*abs(crossface(jxyz)))
				      

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

!!===================================================================================
!! Control Volume snapshot of the mass in a given bin

subroutine mass_snapshot
	use module_record
	use field_io, only : mass_bin_io
	use CV_objects, only : CVcheck_mass, CV_debug!, CV_sphere_mass
    use module_set_parameters, only : mass
	implicit none

	integer										:: n
	integer		,dimension(3)					:: ibin
	real(kind(0.d0)),dimension(:,:,:)  ,allocatable	:: volume_mass_temp
	real(kind(0.d0)),dimension(3)				:: mbinsize

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
		volume_mass_temp(ibin(1),ibin(2),ibin(3)) = volume_mass_temp(ibin(1),ibin(2),ibin(3)) + mass(n)
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


!!===================================================================================
!! Control Volume Momentum continuity
!!===================================================================================


!!===================================================================================
!! Momentum Flux over a surface of a bin including all intermediate bins

module cumulative_momentum_flux_mod

contains

subroutine cumulative_momentum_flux(r_,v_,momentum_flux_,notcrossing)
	use module_record, only : vflux_outflag, domain, halfdomain, planespacing, CV_debug, & 
							  delta_t, planes, Pxy_plane, nplanes, np, nbins, nhb, iter
	use CV_objects, only : CV_constraint!, CV_sphere_momentum
    use librarymod, only : imaxloc, heaviside  =>  heaviside_a1
	use interfaces, only : error_abort
    use module_record, only : get_bin
    use module_set_parameters, only : mass
	implicit none

	real(kind(0.d0)),dimension(:,:),allocatable,intent(in) 			:: r_,v_
	real(kind(0.d0)),dimension(:,:,:,:,:),allocatable,intent(inout) :: momentum_flux_
	integer,dimension(:),allocatable,intent(out),optional			:: notcrossing


	integer							:: ixyz,jxyz,i,j,k,n
	integer							:: planeno
	integer		,dimension(3)		:: ibin1,ibin2,cbin
    real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	real(kind(0.d0))				:: crossplane,rplane,shift
	real(kind(0.d0)),dimension(3)	:: mbinsize,velvect,crossface
	real(kind(0.d0)),dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

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
				Pxy_plane(:,planeno) = Pxy_plane(:,planeno) + mass(n)*velvect(:)!*crossplane

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
					      + mass(n)*velvect(:)*dble(onfacexb)*abs(crossface(jxyz))
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,2) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,2) & 
					      + mass(n)*velvect(:)*dble(onfaceyb)*abs(crossface(jxyz))
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,3) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,3) &
					      + mass(n)*velvect(:)*dble(onfacezb)*abs(crossface(jxyz))
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,4) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,4) &
					      + mass(n)*velvect(:)*dble(onfacext)*abs(crossface(jxyz))
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,5) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,5) &
					      + mass(n)*velvect(:)*dble(onfaceyt)*abs(crossface(jxyz))
					momentum_flux_(cbin(1),cbin(2),cbin(3),:,6) = & 
						momentum_flux_(cbin(1),cbin(2),cbin(3),:,6) &
					      + mass(n)*velvect(:)*dble(onfacezt)*abs(crossface(jxyz))

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
	real(kind(0.d0)),dimension(3)	:: mbinsize
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
	use field_io, only : momentum_bin_io
	use module_record
    use module_set_parameters, only : mass
	!use CV_objects, only : CV_sphere_momentum
	implicit none

	integer											:: n
	integer		,dimension(3)						:: ibin
	real(kind(0.d0)) ,dimension(:,:,:)  ,allocatable		:: volume_mass_temp
	real(kind(0.d0))								:: binvolume
	real(kind(0.d0)),dimension(3)					:: mbinsize
	real(kind(0.d0)),dimension(:,:,:,:),allocatable :: volume_momentum_temp

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

		volume_mass_temp(ibin(1),ibin(2),ibin(3)) = volume_mass_temp(ibin(1),ibin(2),ibin(3)) + mass(n)
		volume_momentum_temp(ibin(1),ibin(2),ibin(3),:) = volume_momentum_temp(ibin(1),ibin(2),ibin(3),:) + mass(n)*v(:,n)
	enddo
	binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
	volume_momentum_temp = volume_momentum_temp/binvolume

	!Output Control Volume momentum change and fluxes
	call momentum_bin_io(volume_mass_temp,volume_momentum_temp,'snap')

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
							  delta_t, planes, Pxyv_plane, nplanes, np, nbins, nhb, & 
                              potenergymol, get_bin, a
	use CV_objects, only : CVcheck_energy !, CV_sphere_momentum
    use librarymod, only : imaxloc, heaviside  =>  heaviside_a1
	use interfaces, only : error_abort
    use module_set_parameters, only : mass
	implicit none

	real(kind(0.d0)),dimension(:,:),allocatable,intent(in) 			:: r_,v_
	real(kind(0.d0)),dimension(:,:,:,:),allocatable,intent(inout) :: energy_flux_

	integer							:: ixyz,jxyz,i,j,k,n,planeno
	real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	real(kind(0.d0))				:: crosstime,crossplane,rplane,shift,energy
	real(kind(0.d0))				:: crosstimetop,crosstimebot,frac,delta_t_cross,potenergy_cross
	real(kind(0.d0)),dimension(3)	:: mbinsize,velvect,crossface
	real(kind(0.d0)),dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

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
				energy = 0.5d0 *  (mass(n)*dot_product(velvect,velvect) + potenergymol(n))

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
					velvect(:) = v_(:,n) + 0.5d0*a(:,n)*delta_t
					energy = 0.5d0 * ( mass(n)*dot_product(velvect,velvect) + potenergymol(n))
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
	use module_record
	use CV_objects, only :  CVcheck_energy
	use cumulative_energy_flux_mod, only : cumulative_energy_flux
	implicit none

	integer,intent(in)	:: flag
	integer, save		:: sample_count

    integer,dimension(3):: thermbinstop,thermbinsbot
	real(kind(0.d0)),dimension(3)	:: ebinsize

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
    use module_set_parameters, only : mass
	implicit none

	integer											:: n
	integer		,dimension(3)						:: ibin
	real(kind(0.d0))								:: binvolume, energy
	real(kind(0.d0)),dimension(3)					:: mbinsize,velvect
	real(kind(0.d0)),dimension(:,:,:),allocatable 	:: volume_energy_temp

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
		energy = 0.5d0 * ( mass(n)*dot_product(velvect,velvect) + potenergymol(n))

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






subroutine surface_density_averaging(flag)
	!use field_io, only : mass_flux_io
	use module_record
	implicit none

	integer			                :: flag
    real(kind(0.d0)),dimension(3)   :: mbinsize
	integer, save	                :: sample_count

	!Only average if mass averaging turned on
	if (flag .eq. 0) return

	call cumulative_surface_density
	sample_count = sample_count + 1
	if (sample_count .eq. Nsurfm_ave) then
		call surface_density_io
		sample_count = 0
		surface_density = 0
	endif

end subroutine surface_density_averaging

!===================================================================================
! Density of molecules found on the surface of a bin 
! Includes all intermediate bins, methodology from 
! " A technique for the calculation of mass, energy, and momentum densities
!   at planes in molecular dynamics simulations"
!  By Peter J. Daivis,a) Karl P. Travis, and B. D. Todd
! 

subroutine cumulative_surface_density
	use module_record
    use librarymod, only : imaxloc, heaviside  =>  heaviside_a1
    use module_set_parameters, only : mass
    !use CV_objects, only : CV_sphere_mass
    implicit none

	integer							:: jxyz,i,j,k,n
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	real(kind(0.d0)),dimension(3)	:: mbinsize,crossface,velvect
	real(kind(0.d0)),dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

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
				velvect(:) = v(:,n) !- a(:,n) * crosstime

				!Add Density on surface for molecules crossing face
				surface_density(cbin(1),cbin(2),cbin(3),1) = & 
					surface_density(cbin(1),cbin(2),cbin(3),1) & 
				      + mass(n)*nint(dble(onfacexb)*crossface(jxyz))/abs(velvect(1))
				surface_density(cbin(1),cbin(2),cbin(3),2) = & 
					surface_density(cbin(1),cbin(2),cbin(3),2) & 
				      + mass(n)*nint(dble(onfaceyb)*crossface(jxyz))/abs(velvect(2))
				surface_density(cbin(1),cbin(2),cbin(3),3) = & 
					surface_density(cbin(1),cbin(2),cbin(3),3) &
				      + mass(n)*nint(dble(onfacezb)*crossface(jxyz))/abs(velvect(3))
				surface_density(cbin(1),cbin(2),cbin(3),4) = & 
					surface_density(cbin(1),cbin(2),cbin(3),4) &
				      + mass(n)*nint(dble(onfacext)*crossface(jxyz))/abs(velvect(1))
				surface_density(cbin(1),cbin(2),cbin(3),5) = & 
					surface_density(cbin(1),cbin(2),cbin(3),5) &
				      + mass(n)*nint(dble(onfaceyt)*crossface(jxyz))/abs(velvect(2))
				surface_density(cbin(1),cbin(2),cbin(3),6) = & 
					surface_density(cbin(1),cbin(2),cbin(3),6) &
				      + mass(n)*nint(dble(onfacezt)*crossface(jxyz))/abs(velvect(3))
				      

				!if (onfacexb .ne. 0) print*, n, i,j,k,ibin1,ibin2,bintop,halfdomain

				!if (cbin(1) .ge. nbins(1)+1 .or. cbin(1) .le. 2) then
				!	print'(4i8,6f10.5)',iter, cbin, momentum_flux(cbin(1),cbin(2),cbin(3),:,1),momentum_flux(cbin(1),cbin(2),cbin(3),:,4)
				!endif

			enddo
			enddo
			enddo

		endif

	enddo

end subroutine cumulative_surface_density

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
	real(kind(0.d0)),intent(in)     			:: accijmag    !Non directional component of acceleration
	real(kind(0.d0)),dimension(3),intent(in)   	:: rij         !vector between particles i and j

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
	real(kind(0.d0)),dimension(3)	:: ri, rj, fij,crossplane,fsurface
	real(kind(0.d0)),dimension(3)	:: Fbinsize, bintopi, binboti, bintopj, binbotj

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


	real(kind(0.d0)),intent(in),dimension(3)	:: ri,rj,fij

	integer							:: i,j,k,ixyz,face
	integer,dimension(3)			:: cbin, ibin, jbin
    real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb,fijvi,fijvj
	real(kind(0.d0)),dimension(3)	:: rij,fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb,velvect
	real(kind(0.d0)),dimension(3)	:: Fbinsize, bintop, binbot, vi_t,vj_t,vi_tmdt,vj_tmdt

	!real(kind(0.d0)),allocatable,dimension(:,:,:),save 	:: fij_dmt
	!real(kind(0.d0)),allocatable,dimension(:,:,:),save 	:: fij_dmt

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


	real(kind(0.d0)),intent(in),dimension(3)	:: ri,rj,fij
	real(kind(0.d0)),dimension(3),intent(in)	:: vi_t


	integer							:: i,j,k,ixyz,face
	integer,dimension(3)			:: cbin, ibin, jbin
    real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb,fijvi,fijvj
	real(kind(0.d0)),dimension(3)	:: rij,fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb,velvect
	real(kind(0.d0)),dimension(3)	:: Fbinsize, bintop, binbot

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













!!===================================================================================
!!Forces times velocity over the surface of a Volume

!module get_timesteps_module

!contains

!subroutine get_timesteps(ncrossings,bin,bin_mdt,Fbinsize,rc,vc,delta_t,delta_t_list,count_t)
!	use computational_constants_MD, only : halfdomain, nhb, iter
!	implicit none

!	integer,intent(in)								:: ncrossings, bin(3),bin_mdt(3)
!	integer,intent(inout)							:: count_t
!	real(kind(0.d0)),intent(in)						:: delta_t
!	real(kind(0.d0)),dimension(3),intent(in)		:: rc,vc,Fbinsize
!	real(kind(0.d0)),dimension(:),allocatable,intent(inout)	:: delta_t_list

!	integer											:: i,j,k,m
!	real(kind(0.d0)),dimension(3)					:: bintop,binbot
!	real(kind(0.d0)),dimension(6)					:: crosstime


!	if (ncrossings .eq. 0) return

!       !Otherwise, get time spent in each cell
!	do i = bin(1),bin_mdt(1),sign(1,bin_mdt(1)-bin(1))
!    do j = bin(2),bin_mdt(2),sign(1,bin_mdt(2)-bin(2))
!    do k = bin(3),bin_mdt(3),sign(1,bin_mdt(3)-bin(3))

!		!Get bin top and bottom
!	    bintop(1) = (i-1*nhb(1)  )*Fbinsize(1)-halfdomain(1)
!	    bintop(2) = (j-1*nhb(2)  )*Fbinsize(2)-halfdomain(2)
!	    bintop(3) = (k-1*nhb(3)  )*Fbinsize(3)-halfdomain(3)
!	    binbot(1) = (i-1*nhb(1)-1)*Fbinsize(1)-halfdomain(1)
!	    binbot(2) = (j-1*nhb(2)-1)*Fbinsize(2)-halfdomain(2)
!	    binbot(3) = (k-1*nhb(3)-1)*Fbinsize(3)-halfdomain(3)

!	    !Calculate the plane intersect of line with surfaces of the cube
!		crosstime(1) = (rc(1) - bintop(1))/ vc(1)
!		crosstime(2) = (rc(1) - binbot(1))/ vc(1)
!		crosstime(3) = (rc(2) - bintop(2))/ vc(2)
!		crosstime(4) = (rc(2) - binbot(2))/ vc(2)
!		crosstime(5) = (rc(3) - bintop(3))/ vc(3)
!		crosstime(6) = (rc(3) - binbot(3))/ vc(3)

!! 		print'(a,i6,3i3,6f12.7)','Crossingtimes',iter,i,j,k,crosstime
!!  		print'(a,i6,3i3,9f12.7)','Surfaces',iter,i,j,k,	bintop(1),rc(1),binbot(1), & 
!!  														bintop(2),rc(2),binbot(2), & 
!!  														bintop(3),rc(3),binbot(3) 

!		!Add any crossings within the time period to the list of crossings
!		do m = 1,6
!			if (crosstime(m) .gt. 0.d0 .and. crosstime(m) .lt. delta_t) then
!				count_t = count_t + 1
!				delta_t_list(count_t) = crosstime(m)
!			endif
!		enddo

!	enddo
!	enddo
!	enddo

!end subroutine get_timesteps


!subroutine get_CV_surface_contributions(ibin,jbin,ri,rj,Fbinsize,value,CV_Face_value,molnoi,molnoj)
!	use computational_constants_MD, only : halfdomain, nhb, iter
!    use librarymod, only : heaviside => heaviside_a1
!	implicit none

!	integer, dimension(3),intent(in)								:: ibin,jbin
!	integer,intent(in)												:: molnoi,molnoj !TEMPTEMP
!	real(kind(0.d0)), dimension(3),intent(in)						:: ri, rj, Fbinsize
!	real(kind(0.d0)), intent(in)									:: value
!	real(kind(0.d0)),dimension(:,:,:,:),allocatable, intent(inout)	:: CV_Face_value

!	integer										:: i,j,k,ixyz, face
!    real(kind(0.d0))							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
!	real(kind(0.d0)),dimension(3)   			:: Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
!	real(kind(0.d0)),dimension(3)				:: rij,  bintop, binbot


!	!If same bin, nothing to do here
!	if (ibin(1) .eq. jbin(1) .and. ibin(2) .eq. jbin(2) .and. ibin(3) .eq. jbin(3)) return

!	!Get interaction line
!	rij = ri - rj

!	!Prevent Division by zero
!	do ixyz = 1,3
!		if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
!	enddo
!		
!	!Loop through all intermediate bins, check surface to add to and then add
!	do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
!	do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
!	do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

!		!Get bin top and bottom
!	    bintop(1) = (i-1*nhb(1)  )*Fbinsize(1)-halfdomain(1)
!	    bintop(2) = (j-1*nhb(2)  )*Fbinsize(2)-halfdomain(2)
!	    bintop(3) = (k-1*nhb(3)  )*Fbinsize(3)-halfdomain(3)
!	    binbot(1) = (i-1*nhb(1)-1)*Fbinsize(1)-halfdomain(1)
!	    binbot(2) = (j-1*nhb(2)-1)*Fbinsize(2)-halfdomain(2)
!	    binbot(3) = (k-1*nhb(3)-1)*Fbinsize(3)-halfdomain(3)

!		!Calculate the plane intersect of line with surfaces of the cube
!		Pxt=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
!		Pxb=(/ binbot(1),ri(2)+(rij(2)/rij(1))*(binbot(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(binbot(1)-ri(1))  /)
!		Pyt=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
!		Pyb=(/ri(1)+(rij(1)/rij(2))*(binbot(2)-ri(2)), binbot(2),ri(3)+(rij(3)/rij(2))*(binbot(2)-ri(2))  /)
!		Pzt=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)
!		Pzb=(/ri(1)+(rij(1)/rij(3))*(binbot(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(binbot(3)-ri(3)), binbot(3) /)

!		onfacexb =  	(sign(1.d0,binbot(1)- rj(1)) - sign(1.d0,binbot(1)- ri(1)))* &
!						(heaviside(bintop(2)-Pxb(2)) - heaviside(binbot(2)-Pxb(2)))* &
!						(heaviside(bintop(3)-Pxb(3)) - heaviside(binbot(3)-Pxb(3)))
!		onfaceyb =  	(sign(1.d0,binbot(2)- rj(2)) - sign(1.d0,binbot(2)- ri(2)))* &
!						(heaviside(bintop(1)-Pyb(1)) - heaviside(binbot(1)-Pyb(1)))* &
!						(heaviside(bintop(3)-Pyb(3)) - heaviside(binbot(3)-Pyb(3)))
!		onfacezb =  	(sign(1.d0,binbot(3)- rj(3)) - sign(1.d0,binbot(3)- ri(3)))* &
!						(heaviside(bintop(1)-Pzb(1)) - heaviside(binbot(1)-Pzb(1)))* &
!						(heaviside(bintop(2)-Pzb(2)) - heaviside(binbot(2)-Pzb(2)))

!		onfacext =  	(sign(1.d0,bintop(1)- rj(1)) - sign(1.d0,bintop(1)- ri(1)))* &
!						(heaviside(bintop(2)-Pxt(2)) - heaviside(binbot(2)-Pxt(2)))* &
!	            		(heaviside(bintop(3)-Pxt(3)) - heaviside(binbot(3)-Pxt(3)))
!		onfaceyt = 		(sign(1.d0,bintop(2)- rj(2)) - sign(1.d0,bintop(2)- ri(2)))* &
!						(heaviside(bintop(1)-Pyt(1)) - heaviside(binbot(1)-Pyt(1)))* &
!						(heaviside(bintop(3)-Pyt(3)) - heaviside(binbot(3)-Pyt(3)))
!		onfacezt =  	(sign(1.d0,bintop(3)- rj(3)) - sign(1.d0,bintop(3)- ri(3)))* &
!						(heaviside(bintop(1)-Pzt(1)) - heaviside(binbot(1)-Pzt(1)))* &
!						(heaviside(bintop(2)-Pzt(2)) - heaviside(binbot(2)-Pzt(2)))

!		!Value acting on face
!		CV_Face_value(i,j,k,1) = CV_Face_value(i,j,k,1) + value*onfacexb
!		CV_Face_value(i,j,k,2) = CV_Face_value(i,j,k,2) + value*onfaceyb
!		CV_Face_value(i,j,k,3) = CV_Face_value(i,j,k,3) + value*onfacezb
!		CV_Face_value(i,j,k,4) = CV_Face_value(i,j,k,4) + value*onfacext
!		CV_Face_value(i,j,k,5) = CV_Face_value(i,j,k,5) + value*onfaceyt
!		CV_Face_value(i,j,k,6) = CV_Face_value(i,j,k,6) + value*onfacezt

!  	if (i .eq. 3 .and. j .eq. 3 .and. k .eq. 3) then
!   	if (any(abs((/onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb/)) .gt. 0.00001)) then
!   	!if (abs(onfacext+onfacexb+onfaceyt+onfaceyb+onfacezt+onfacezb) .gt. 0.00001) then
!   			if (abs(onfacexb) .gt. 0.000001) face = 1 
!   			if (abs(onfaceyb) .gt. 0.000001) face = 2 
!   			if (abs(onfacezb) .gt. 0.000001) face = 3 
!   			if (abs(onfacext) .gt. 0.000001) face = 4 
!   			if (abs(onfaceyt) .gt. 0.000001) face = 5 
!   			if (abs(onfacezt) .gt. 0.000001) face = 6 
!   				print'(a,6i5,2f12.4,6i3,f6.0,5f4.0)','Int pp box 3', iter,molnoi,molnoj,i,j,k,value, & 
!   																	CV_Face_value(i,j,k,face), & 
!   																	ibin,jbin,onfacext,onfacexb,onfaceyt, & 
!   																	onfaceyb,onfacezt,onfacezb
!   	!endif
!   	endif
!   	endif

!		!print'(a,3i4,7f10.5)', 'IN surface cv routine ',i,j,k, CV_Face_value(i,j,k,:), value

!	enddo
!	enddo
!	enddo

!end subroutine get_CV_surface_contributions

!end module get_timesteps_module

!!subroutine control_volume_power_partialint(fij,ri,rj,vi_mhdt,vj_mhdt,ai_mdt,aj_mdt,ai,aj,molnoi,molnoj)
!!    use module_record
!!	use CV_objects, only : CV_debug,CV_constraint
!!    use librarymod, only : heaviside => heaviside_a1, bubble_sort
!!	use get_timesteps_module
!!    implicit none

!!	integer,intent(in)						 	:: molnoi, molnoj
!!	real(kind(0.d0)),dimension(3),intent(in) 	:: fij,ri,rj,vi_mhdt,vj_mhdt,ai_mdt,aj_mdt,ai,aj


!!	integer										:: i,j,k,m,n,ixyz,ncrossings,ncrossingsi,ncrossingsj,count_t
!!	integer,dimension(3)						:: cbin, ibin, jbin, ibin_mdt, jbin_mdt
!!    real(kind(0.d0))							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
!!    real(kind(0.d0))							:: invrij2,rij2_mdt,fijvi,fijvi_mdt,fijvi_trapz,delta_t_portion,eps = 0.0001d0
!!	real(kind(0.d0)),dimension(3)   			:: fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb,velvect
!!	real(kind(0.d0)),dimension(3)   			:: vi,vj,vi_mdt,vj_mdt
!!	real(kind(0.d0)),dimension(3)   			:: ri_mdt,rj_mdt,rij_mdt,fij_mdt,ri_p,rj_p
!!	real(kind(0.d0)),dimension(3)				:: Fbinsize, bintop, binbot
!!	real(kind(0.d0)),dimension(6)				:: crosstime
!!	real(kind(0.d0)),dimension(:),allocatable	:: delta_t_list,delta_t_array

!!	character(30)	:: frmtstr

!!	!Determine bin size
!!	Fbinsize(:) = domain(:) / nbins(:)

!!! 	!Old velocities
!!! 	vi_mhdt = vi_hdt - ai_mdt*delta_t
!!! 	vj_mhdt = vj_hdt - aj_mdt*delta_t
!!! 	
!!! 	!Old positions
!!! 	ri_mdt = ri - vi_mhdt*delta_t
!!! 	rj_mdt = rj - vj_mhdt*delta_t

!!! 	!Get fij, vi and vj
!!! 	vi = vi_hdt + 0.5*ai_mdt*delta_t		!This is not okay -- ai_pdt should be used...
!!! 	vj = vi_hdt + 0.5*ai_mdt*delta_t		!This is not okay -- ai_pdt should be used...

!!! 	!Get fij_mdt, vi_mdt and vj_mdt
!!! 	rij_mdt(:)= ri_mdt(:) - rj_mdt(:)   					!Evaluate distance between particle i and j
!!! 	invrij2  = 1.d0/(dot_product(rij_mdt,rij_mdt))					!Invert value
!!! 	fij_mdt = 48.d0*(invrij2**7-0.5d0*invrij2**4)*rij_mdt
!!! 	vi_mdt = vi_hdt - 0.5*ai_mdt*delta_t
!!! 	vj_mdt = vi_hdt - 0.5*ai_mdt*delta_t


!!	!Assuming acceleration passed in are at time t-dt/2

!!	!Velocities for dot product @ t and t-dt
!!	vi(:) 	  = vi_mhdt + 0.5d0*delta_t*ai(:)
!!	vj(:) 	  = vj_mhdt + 0.5d0*delta_t*aj(:)

!!	vi_mdt(:) = vi_mhdt - 0.5d0*delta_t*ai_mdt(:)
!!	vj_mdt(:) = vj_mhdt - 0.5d0*delta_t*aj_mdt(:)

!!	!Velocity=>Positions=>fij at t-dt
!!	!vi_mhdt(:) = vi_hdt - delta_t*ai_mdt(:)
!!	!vj_mhdt(:) = vj_hdt - delta_t*aj_mdt(:)

!! 	ri_mdt(:) = ri - delta_t*vi_mhdt(:)
!! 	rj_mdt(:) = rj - delta_t*vj_mhdt(:)

!! 	rij_mdt(:)= ri_mdt(:) - rj_mdt(:)   					!Evaluate distance between particle i and j
!!	rij2_mdt = dot_product(rij_mdt,rij_mdt)
!!	if (rij2_mdt < rcutoff2) then
!!	 	invrij2  = 1.d0/rij2_mdt
!! 		fij_mdt = 48.d0*(invrij2**7-0.5d0*invrij2**4)*rij_mdt
!!	else
!!		fij_mdt = 0.d0
!!	endif

!!	!Trapizium rule calculation of power to add during portion of timestep
!!	fijvi 	  = dot_product(fij    ,  vi  )
!!	fijvi_mdt = dot_product(fij_mdt,vi_mdt)


!!! 	vi_phdt(:) = vi_hdt + delta_t*ai(:)
!!! 	vj_phdt(:) = vj_hdt + delta_t*aj(:)

!!! 	ri_pdt(:) = ri + delta_t*vi_phdt(:)
!!! 	rj_pdt(:) = rj + delta_t*vj_phdt(:)

!! 	!Get fij_mdt, vi_mdt and vj_mdt
!!!  	rij_pdt(:)= ri_pdt(:) - rj_pdt(:)   					!Evaluate distance between particle i and j
!!!  	invrij2  = 1.d0/(dot_product(rij_pdt,rij_pdt))					!Invert value
!!!  	fij_pdt = 48.d0*(invrij2**7-0.5d0*invrij2**4)*rij_pdt



!!	!Assign to bins using integer division
!!	ibin(:) 	= ceiling((ri(:)	+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!!	jbin(:) 	= ceiling((rj(:)	+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!!	ibin_mdt(:) = ceiling((ri_mdt(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish previous bin
!!	jbin_mdt(:) = ceiling((rj_mdt(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish previous bin

!!	!Get number of crossings
!!	ncrossingsi =  abs(ibin(1) - ibin_mdt(1)) + abs(ibin(2) - ibin_mdt(2)) + abs(ibin(3) - ibin_mdt(3))
!!	ncrossingsj =  abs(jbin(1) - jbin_mdt(1)) + abs(jbin(2) - jbin_mdt(2)) + abs(jbin(3) - jbin_mdt(3)) 
!!	ncrossings = ncrossingsi+ncrossingsj

!!	!Get portion of timestep in each cell
!!	count_t = 0
!!	if (ncrossings .eq. 0) then
!!		allocate(delta_t_array(2))
!!		delta_t_array(1) = 0.d0
!!		delta_t_array(2) = delta_t
!!	else
!!		! Each molecular crossing time is counted twice -- once for the cell it's 
!!		! leaving and once for the cell it's moving into
!!		allocate(delta_t_list(2*ncrossings))
!!		if (ncrossingsi .ne. 0) print'(a,2i4,6f12.8,6i4)', 'Molecule i',iter, molnoi, ri,ri_mdt,ibin,ibin_mdt
!!		if (ncrossingsj .ne. 0) print'(a,2i4,6f12.8,6i4)', 'Molecule j',iter, molnoj, rj,rj_mdt,jbin,jbin_mdt

!!    	!Check if molecule i has changed bin
!!		call  get_timesteps(ncrossingsi,ibin,ibin_mdt,Fbinsize,ri,vi_mhdt,delta_t,delta_t_list,count_t)

!!    	!Check if molecule j has changed bin
!!		call  get_timesteps(ncrossingsj,jbin,jbin_mdt,Fbinsize,rj,vj_mhdt,delta_t,delta_t_list,count_t)

!!		!Get crossing times in chronological order
!!		call bubble_sort(delta_t_list)

!!		!Sanity checks -- count_t == ncrossings & sum(delta_t_list) == delta_t
!!		if (ncrossings .ne. 0.5*count_t) then
!!			print'(a,2i12)', ' count_t == ncrossings ',count_t/2,ncrossings
!!			stop "Error - crossing values not equal"
!!		endif
!!		if (any(delta_t_list .gt. delta_t)) then
!!			stop "Error - delta t not ok in energy CV"
!!		endif

!!		!print'(i4,a,2i12,a,2f10.5)',iter,' count_t == ncrossings ',ncrossings,count_t/2, & 
!!		!						 ' sum(delta_t_list) == delta_t ',delta_t, sum(delta_t_list)

!!		!Set initial and final time and copy timestep arrays 
!!		allocate(delta_t_array(ncrossings+2))
!!		delta_t_array(1) = 0.d0
!!		delta_t_array(ncrossings+2) = delta_t
!!		n = 2
!!		do i = 1,size(delta_t_list),2
!!			if (delta_t_list(i) .ne. delta_t_list(i+1)) then
!!				stop "Error - Successive timesteps not same in delta_t array"
!!			endif 
!!			delta_t_array(n) = delta_t_list(i)
!!			n = n + 1
!!		enddo

!!	endif

!!	if (ncrossings .ge. 1) then
!! 	if (all(ibin .eq. 3)  .or. &
!! 		all(jbin .eq. 3) 	) then
!!!  	if (all(ibin .eq. 3) .and. any(jbin .ne. 3) .or. &
!!!  		all(jbin .eq. 3) .and. any(ibin .ne. 3)	) then
!!		print'(a,12i4)','Predicted cells', ibin,ibin_mdt,jbin,jbin_mdt
!!	endif
!!	endif

!!	!Calculate surface contributions
!!	do n=1,ncrossings+1
!!	
!!		delta_t_portion = delta_t_array(n+1)-delta_t_array(n)

!!		!Calculate intermediate positions
!!		ri_p = ri - vi_mhdt*delta_t_portion + eps
!!		rj_p = rj - vj_mhdt*delta_t_portion + eps

!!		ibin(:) 	= ceiling((ri_p(:)	+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!!		jbin(:) 	= ceiling((rj_p(:)	+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin

!!		if (ncrossings .ge. 1) then
!! 		if (all(ibin .eq. 3) .and. any(jbin .ne. 3) .or. &
!! 			all(jbin .eq. 3) .and. any(ibin .ne. 3)	) then
!!				print'(a,12i4,3f14.9)','Delta_t',iter,n,molnoi,molnoj,ibin,jbin, & 
!!												 ncrossingsi,ncrossingsj, delta_t_array(n),& 
!!												 delta_t_portion,delta_t_array(n+1)
!!		endif
!!		endif

!!		!Calculate and add fraction of interation
!!		fijvi_trapz = 0.5d0*delta_t_portion * (fijvi + fijvi_mdt)/delta_t
!!		call  get_CV_surface_contributions(ibin,jbin,ri_p,rj_p,Fbinsize,fijvi_trapz,Pxyvface2,molnoi,molnoj)

!!!  		if (all(ibin .eq. 3) .and. any(jbin .ne. 3)) then
!!! 			print'(a,3i4,10f11.6)', 'Saving ints',iter,molnoi,molnoj,delta_t_portion , & 
!!! 													0.25*Pxyvface2(ibin(1),ibin(2),ibin(3),:)/5.0, &  !binsize
!!! 													fijvi,fijvi_mdt,fijvi_trapz
!!! 		endif
!!! 		if (all(jbin .eq. 3) .and. any(ibin .ne. 3)) then
!!! 			print'(a,3i4,10f11.6)', 'Saving ints',iter,molnoi,molnoj,delta_t_portion , &
!!! 													0.25*Pxyvface2(jbin(1),jbin(2),jbin(3),:)/5.0, & 
!!! 													fijvi,fijvi_mdt,fijvi_trapz
!!! 		endif

!!! 		if (all(ibin .eq. 3) .and. any(jbin .ne. 3) .or. &
!!! 			all(jbin .eq. 3) .and. any(ibin .ne. 3)	) then
!!! 			print'(a,3i4,10f11.7)', 'Current',iter,molnoi,molnoj,fij_mdt,vi_mdt,fijvi_mdt,0.25*fijvi_trapz/(5.0*delta_t)
!!! 			print'(a,3i4,10f11.7)', 'Future ',iter,molnoi,molnoj,fij,vi,fijvi,0.25*fijvi_trapz/(5.0*delta_t)
!!! 		endif
!!		!print'(a,12f10.5)', 'Power after', Pxyvface(ibin(1),ibin(2),ibin(3),:), Pxyvface(jbin(1),jbin(2),jbin(3),:)

!!	enddo

!!end subroutine control_volume_power_partialint




!! !===================================================================================
!! !Forces over the surface of a Volume optmised for computational efficiency

!! subroutine control_volume_stresses_opt(fij,ri,rj,molnoi)
!! 	use module_record
!!     use librarymod, only : heaviside  =>  heaviside_a1
!! 	implicit none

!! 	integer							:: i,j,k,ixyz,molnoi,molnoj
!! 	integer,dimension(3)			:: cbin, ibin, jbin, Si
!! 	real(kind(0.d0)),dimension(3)	:: ri,rj,rij,fij,fsurface,Px,Py,Pz,sgnjit,sgnjib,onfaceb,onfacet,velvect
!! 	real(kind(0.d0)),dimension(3)	:: Fbinsize, bintop, binbot

!! 	!Calculate rij
!! 	rij = ri - rj
!! 	!Prevent Division by zero
!! 	do ixyz = 1,3
!! 		if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
!! 	enddo

!! 	!Determine bin size
!! 	Fbinsize(:) = domain(:) / nbins(:)

!! 	!Assign to bins using integer division
!! 	ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!! 	jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin

!! 	if (ibin(1) .eq. jbin(1) .and. ibin(2) .eq. jbin(2) .and. ibin(3) .eq. jbin(3)) return
!! 		
!! 	do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
!! 	do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
!! 	do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

!! 		cbin(1) = i; cbin(2) = j; cbin(3) = k

!! 		bintop(:) = (cbin(:)-1*nhb(:)  )*Fbinsize(:)-halfdomain(:)
!! 		binbot(:) = (cbin(:)-1*nhb(:)-1)*Fbinsize(:)-halfdomain(:)

!! 		!Calculate the plane intersect of line with surfaces of the cube
!! 		Px=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
!! 		Py=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
!! 		Pz=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)

!! 		sgnjit(:)= sign(1.d0,bintop(:)- rj(:)) - sign(1.d0,bintop(:)- ri(:))
!! 		sgnjib(:)= sign(1.d0,binbot(:)- rj(:)) - sign(1.d0,binbot(:)- ri(:))

!! 		Si(1) =	(heaviside(bintop(2)-Px(2)) - heaviside(binbot(2)-Px(2)))* &
!! 				(heaviside(bintop(3)-Px(3)) - heaviside(binbot(3)-Px(3)))
!! 		Si(2) =	(heaviside(bintop(1)-Py(1)) - heaviside(binbot(1)-Py(1)))* &
!! 				(heaviside(bintop(3)-Py(3)) - heaviside(binbot(3)-Py(3)))
!! 		Si(3) =	(heaviside(bintop(1)-Pz(1)) - heaviside(binbot(1)-Pz(1)))* &
!! 				(heaviside(bintop(2)-Pz(2)) - heaviside(binbot(2)-Pz(2)))

!! 		onfaceb(:) = sgnjib(:)*dble(Si(:))
!! 		onfacet(:) = sgnjit(:)*dble(Si(:))

!! 		!Stress acting on face over volume
!! 		Pxyface(cbin(1),cbin(2),cbin(3),:,1) = Pxyface(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*onfaceb(1)
!! 		Pxyface(cbin(1),cbin(2),cbin(3),:,2) = Pxyface(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*onfaceb(2)
!! 		Pxyface(cbin(1),cbin(2),cbin(3),:,3) = Pxyface(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*onfaceb(3)
!! 		Pxyface(cbin(1),cbin(2),cbin(3),:,4) = Pxyface(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*onfacet(1)
!! 		Pxyface(cbin(1),cbin(2),cbin(3),:,5) = Pxyface(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*onfacet(2)
!! 		Pxyface(cbin(1),cbin(2),cbin(3),:,6) = Pxyface(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*onfacet(3)

!! 		!Stress acting on face over volume
!! 		if (eflux_outflag .ne. 0) then
!! 			velvect(:) = v(:,molnoi) 
!! 			!if (molnoi .gt. np) print*, velvect(1)
!! 			!velvect(:) = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)
!! 			Pxyvface(cbin(1),cbin(2),cbin(3),1) = Pxyvface(cbin(1),cbin(2),cbin(3),1) + dot_product(fij,velvect)*onfaceb(1)
!! 			Pxyvface(cbin(1),cbin(2),cbin(3),2) = Pxyvface(cbin(1),cbin(2),cbin(3),2) + dot_product(fij,velvect)*onfaceb(2)
!! 			Pxyvface(cbin(1),cbin(2),cbin(3),3) = Pxyvface(cbin(1),cbin(2),cbin(3),3) + dot_product(fij,velvect)*onfaceb(3)
!! 			Pxyvface(cbin(1),cbin(2),cbin(3),4) = Pxyvface(cbin(1),cbin(2),cbin(3),4) + dot_product(fij,velvect)*onfacet(1)
!! 			Pxyvface(cbin(1),cbin(2),cbin(3),5) = Pxyvface(cbin(1),cbin(2),cbin(3),5) + dot_product(fij,velvect)*onfacet(2)
!! 			Pxyvface(cbin(1),cbin(2),cbin(3),6) = Pxyvface(cbin(1),cbin(2),cbin(3),6) + dot_product(fij,velvect)*onfacet(3)
!! 		endif

!! 		!Force applied to volume
!! 		fsurface(:) = 0.d0
!! 		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*(onfaceb(1) - onfacet(1))
!! 		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*(onfaceb(2) - onfacet(2))
!! 		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*(onfaceb(3) - onfacet(3))
!! 		volume_force(cbin(1),cbin(2),cbin(3),:,1) = volume_force(cbin(1),cbin(2),cbin(3),:,1) + fsurface*delta_t

!! 	enddo
!! 	enddo
!! 	enddo

!! end subroutine control_volume_stresses_opt


!!===================================================================================
!!Forces over the surface of a Volume further optmised for computational efficiency

!!subroutine control_volume_stresses_opt_2(fij,ri,rj,molnoi)
!!use module_record
!!implicit none

!!	integer							:: i,j,k,ixyz,molnoi
!	!integer							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
!!	integer,dimension(3)			:: cbin, ibin, jbin
!!	integer,dimension(18)			:: hfacelimits
!!	real(kind(0.d0)),dimension(3)	:: ri,rj,rij,fij,fsurface,Px,Py,Pz,Si,sgnjit,sgnjib,onfaceb,onfacet,velvect
!!	real(kind(0.d0)),dimension(3)	:: Fbinsize, bintop, binbot
!!	real(kind(0.d0)),dimension(18)	:: facelimits

!	!Calculate rij
!!	rij = ri - rj
!	!Prevent Division by zero
!!	do ixyz = 1,3
!!		if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
!!	enddo

!	!Determine bin size
!!	Fbinsize(:) = domain(:) / nbins(:)

!	!Assign to bins using integer division
!!	ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!!	jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+nhb(:)	!Establish current bin
!!
!!	if (ibin(1) .eq. jbin(1) .and. ibin(2) .eq. jbin(2) .and. ibin(3) .eq. jbin(3)) return
!		
!!	do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
!!	do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
!!	do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

!!		cbin(1) = i; cbin(2) = j; cbin(3) = k

!!		bintop(:) = (cbin(:)-1*nhb(:)  )*Fbinsize(:)-halfdomain(:)
!!		binbot(:) = (cbin(:)-1*nhb(:)-1)*Fbinsize(:)-halfdomain(:)

!		!Calculate the plane intersect of line with surfaces of the cube
!!		Px=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
!!		Py=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
!!		Pz=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)

!!		facelimits(1 :3 ) = bintop-Px	
!!		facelimits(4 :6 ) = bintop-Py 		
!!		facelimits(7 :9 ) = bintop-Pz
!!		facelimits(10:12) = binbot-Px	
!!		facelimits(13:15) = binbot-Py 		
!!		facelimits(16:18) = binbot-Pz

!!		hfacelimits = heaviside(facelimits)

!!		Si(1) =	(hfacelimits(2) - hfacelimits(11))* &
!!				(hfacelimits(3) - hfacelimits(12))
!!		Si(2) =	(hfacelimits(4) - hfacelimits(13))* &
!!				(hfacelimits(6) - hfacelimits(15))
!!		Si(3) =	(hfacelimits(7) - hfacelimits(16))* &
!!				(hfacelimits(8) - hfacelimits(17))

!!		sgnjit(:)= sign(1.d0,bintop(:)- rj(:)) - sign(1.d0,bintop(:)- ri(:))
!!		sgnjib(:)= sign(1.d0,binbot(:)- rj(:)) - sign(1.d0,binbot(:)- ri(:))

!!		onfaceb =  	sgnjib*Si
!!		onfacet =  	sgnjit*Si

!		!Stress acting on face over volume
!!		Pxyface(cbin(1),cbin(2),cbin(3),:,1) = Pxyface(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*dble(onfaceb(1))
!!		Pxyface(cbin(1),cbin(2),cbin(3),:,2) = Pxyface(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*dble(onfaceb(2))
!!		Pxyface(cbin(1),cbin(2),cbin(3),:,3) = Pxyface(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*dble(onfaceb(3))
!!		Pxyface(cbin(1),cbin(2),cbin(3),:,4) = Pxyface(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*dble(onfacet(1))
!!		Pxyface(cbin(1),cbin(2),cbin(3),:,5) = Pxyface(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*dble(onfacet(2))
!!		Pxyface(cbin(1),cbin(2),cbin(3),:,6) = Pxyface(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*dble(onfacet(3))

!!		!Stress acting on face over volume
!!		if (eflux_outflag .ne. 0) then
!!			velvect(:) = v(:,molnoi) 
!			!if (molnoi .gt. np) print*, velvect(1)
!			!velvect(:) = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)
!!			Pxyvface(cbin(1),cbin(2),cbin(3),1) = Pxyvface(cbin(1),cbin(2),cbin(3),1) + dot_product(fij,velvect)*dble(onfaceb(1))
!!			Pxyvface(cbin(1),cbin(2),cbin(3),2) = Pxyvface(cbin(1),cbin(2),cbin(3),2) + dot_product(fij,velvect)*dble(onfaceb(2))
!!			Pxyvface(cbin(1),cbin(2),cbin(3),3) = Pxyvface(cbin(1),cbin(2),cbin(3),3) + dot_product(fij,velvect)*dble(onfaceb(3))
!!			Pxyvface(cbin(1),cbin(2),cbin(3),4) = Pxyvface(cbin(1),cbin(2),cbin(3),4) + dot_product(fij,velvect)*dble(onfacet(1))
!!			Pxyvface(cbin(1),cbin(2),cbin(3),5) = Pxyvface(cbin(1),cbin(2),cbin(3),5) + dot_product(fij,velvect)*dble(onfacet(2))
!!			Pxyvface(cbin(1),cbin(2),cbin(3),6) = Pxyvface(cbin(1),cbin(2),cbin(3),6) + dot_product(fij,velvect)*dble(onfacet(3))
!!		endif

!		!Force applied to volume
!!		fsurface(:) = 0.d0
!!		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfaceb(1) - onfacet(1))
!!		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfaceb(2) - onfacet(2))
!!		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfaceb(3) - onfacet(3))
!!		volume_force(cbin(1),cbin(2),cbin(3),:,1) = volume_force(cbin(1),cbin(2),cbin(3),:,1) + fsurface*delta_t

!!	enddo
!!	enddo
!!	enddo

!!end subroutine control_volume_stresses_opt_2

!!====================================================================================

subroutine pressure_tensor_forces_MOP(pnxyz,ri,rj,rij,accijmag)
	use module_record
    use librarymod, only : heaviside  =>  heaviside_a1
	implicit none

	integer							:: n
	integer							:: pnxyz	 !Plane normal direction
	integer							:: planenoi,planenoj
	real(kind(0.d0))                :: shift, plane !Plane normal components i and j
	real(kind(0.d0))                :: accijmag      !Non directional component of acceleration
	real(kind(0.d0)),dimension(3)   :: ri, rj, rij   !Vector between particles i and j
	real(kind(0.d0)),dimension(3)   :: Pyb           !Location of intercept with plane

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

!!===================================================================================
!! Record external forces applied to molecules inside a volume

module module_record_external_forces

contains

subroutine record_external_forces(F,ri,vi)
	use module_record, only : domain,halfdomain, nbins, nhb
	use calculated_properties_MD, only :  F_ext_bin,Fv_ext_bin
	use computational_constants_MD, only : eflux_outflag, delta_t,iter
	use librarymod, only : imaxloc
	implicit none

	real(kind(0.d0)),dimension(3),intent(in):: F,ri
	real(kind(0.d0)),dimension(3),intent(in),optional :: vi

	integer									:: jxyz
	integer	,dimension(3)					:: ibin, ibin_pt, crossface
	real(kind(0.d0))						:: crosstimetop,crosstimebot,delta_t_cross, Fiextvi
	real(kind(0.d0)),dimension(3)			:: mbinsize, ri_pt, bintop, binbot

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


!subroutine evaluate_U
!	use module_record
!	use linked_list
!	implicit none

!	integer				:: n
!	integer,dimension(3):: ib
!	real(kind(0.d0)),dimension(3) 	:: Vbinsize

!	real(kind(0.d0)), dimension(:,:,:), allocatable :: mbin
!	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: vbin

!	integer :: x,y,z,c

!	Vbinsize(:) = domain(:) / nbins(:)
!	
!	x = nbins(1) + 2*nhb(1)
!	y = nbins(2) + 2*nhb(2)
!	z = nbins(3) + 2*nhb(3)
!	c = 3
!	
!	allocate(vbin(x,y,z,c))
!	allocate(mbin(x,y,z))
!	
!	vbin = 0.d0
!	mbin = 0
!		
!	do n = 1, np

!		! Get bin
!		ib(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
!		! Add v, m to bin
!		vbin(ib(1),ib(2),ib(3),:) = vbin(ib(1),ib(2),ib(3),:) + v(:,n)
!		mbin(ib(1),ib(2),ib(3)) = mbin(ib(1),ib(2),ib(3)) + 1 

!	end do

!	do n = 1, np

!		! Get bin
!		ib(:) = ceiling((r(:,n)+halfdomain(:))/Vbinsize(:)) + nhb
!		! Get U from vbin / mbin
!		U(:,n) = vbin(ib(1),ib(2),ib(3),:)/real(mbin(ib(1),ib(2),ib(3)),kind(0.d0))
!	
!	end do
!	
!	deallocate(vbin)
!	deallocate(mbin)

!end subroutine evaluate_U


module cluster_analysis
    use linked_list, only : clusterinfo, node
	implicit none

contains

    subroutine build_clusters(self, rd)
	    use module_compute_forces
        use computational_constants_MD, only : iter
	    use interfaces, only : error_abort
	    implicit none

        type(clusterinfo),intent(inout)    :: self
        double precision, intent(in)        :: rd

        integer                         :: i, sumi
        !Initialised cluster list
        if (.not. allocated(self%Nlist)) then
            allocate(self%Nlist(np+extralloc))
            allocate(self%inclust(np+extralloc))
            allocate(self%head(np+extralloc))
        endif

        self%Nclust = 0
	    do i = 1,np+extralloc
            self%inclust(i) = 0
		    self%Nlist(i) = 0	        !Zero number of molecules in cluster list
		    nullify(self%head(i)%point) !Nullify cluster list head pointer 
	    enddo

        !Call to build clusters from neighbour and cell lists
        call build_from_cellandneighbour_lists(self, cell, neighbour, rd, r, np, skipwalls_=.true.)

        !Remove all empty cluster references
        call CompressClusters(self)

    end subroutine build_clusters

    subroutine get_cluster_properties(self, rd)
        use physical_constants_MD, only : pi, np
        use physical_constants_MD, only : tethereddisttop, tethereddistbottom
        use computational_constants_MD, only : iter, thermo_tags, thermo, free, globaldomain 
        use librarymod, only : imaxloc, get_Timestep_FileName, least_squares, get_new_fileunit
        use minpack_fit_funcs_mod, only : fn, cubic_fn, curve_fit
        use arrays_MD, only : tag, r
        use librarymod, only : heaviside  =>  heaviside_a1
        implicit none

        type(clusterinfo),intent(inout)    :: self
        double precision, intent(in)       :: rd

        logical                         :: first_time=.true., print_debug
        character(32)                   :: filename, debug_outfile
        integer                         :: n,pid,i,j,resolution,fileunit
        double precision                :: tolerence, m, c, cl_angle, theta_i, yi
        double precision, dimension(3)  :: bintopi, binboti, ri 
        double precision, dimension(4)  :: pt = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
        double precision, dimension(4)  :: pb = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
        double precision,dimension(6)   :: extents
        double precision,dimension(:),allocatable :: x,y,f
        double precision,dimension(:,:),allocatable :: rnp, extents_grid


        resolution = 10; tolerence = rd
        call cluster_global_extents(self, imaxloc(self%Nlist), extents)
        call cluster_extents_grid(self, imaxloc(self%Nlist), 1, resolution, & 
                                  extents_grid)!, debug_outfile='./results/maxcell_top')
        call cluster_outer_mols(self, imaxloc(self%Nlist), tolerence=tolerence, dir=1, & 
                                rmols=rnp, extents=extents_grid)!, debug_outfile='./results/clust_edge_top')

        !Curve fits to clusers
        allocate(x(size(rnp,2)),y(size(rnp,2)))
        x = rnp(2,:); y = rnp(1,:)
        !Linear
        call least_squares(x, y, m, c)
        !Cubic using minpack
        fn => cubic_fn
        call curve_fit(fn, x, y, pt, f)
        deallocate(x,y)
        cl_angle = 90.d0+atan(m)*180.d0/pi
    	fileunit = get_new_fileunit()
        if (first_time) then
            open(unit=fileunit,file='./results/linecoeff_top',status='replace')
        else
            open(unit=fileunit,file='./results/linecoeff_top',access='append')
        endif
        write(fileunit,'(i12, 7f15.8)'), iter, m, c, cl_angle, pt
        !write(fileunit,'(i12, 3(a,f10.5))'), iter, ' Top line    y = ', m, ' x + ',c , ' angle = ', cl_angle
        close(fileunit,status='keep')

        call cluster_extents_grid(self, imaxloc(self%Nlist), 4, resolution, &
                                  extents_grid )!, debug_outfile='./results/maxcell_bot')
        call cluster_outer_mols(self, imaxloc(self%Nlist), tolerence=tolerence, dir=4, & 
                                rmols=rnp, extents=extents_grid)!, debug_outfile='./results/clust_edge_bot')

        !Curve fits to clusers
        allocate(x(size(rnp,2)),y(size(rnp,2)))
        x = rnp(2,:); y = rnp(1,:)
        !Linear
        call least_squares(x, y, m, c)
        !Cubic using minpack
        fn => cubic_fn
        call curve_fit(fn, x, y, pb, f)
        deallocate(x,y)
        cl_angle = 90.d0+atan(m)*180.d0/pi
    	fileunit = get_new_fileunit()
        if (first_time) then
            open(unit=fileunit,file='./results/linecoeff_bot',status='replace')
            first_time = .false.
        else
            open(unit=fileunit,file='./results/linecoeff_bot',access='append')
        endif
        write(fileunit,'(i12, 7f15.8)'), iter, m, c, cl_angle, pb
        !write(fileunit,'(i12, 3(a,f10.5))'), iter, ' Bottom line y = ', m, ' x + ',c  , ' angle = ', cl_angle
        close(fileunit,status='keep')

        !If cluster has broken up, stop simulation
        call check_for_cluster_breakup(self)

        !Front/back surfaces in y
        bintopi(2) =  0.5d0*globaldomain(2) - tethereddisttop(2)
        binboti(2) = -0.5d0*globaldomain(2) + tethereddistbottom(2)
        
        !Left/Right cluser based surfaces in z
        bintopi(3) =  0.5d0*globaldomain(3) 
        binboti(3) = -0.5d0*globaldomain(3) 


        print_debug = .false.
        debug_outfile = './results/CV_mols'
        if (print_debug) then
            pid = get_new_fileunit()
            call get_Timestep_FileName(iter,debug_outfile,filename)
            open(unit=pid,file=trim(filename),status='replace')
            do n =1,np

                ri(:) = r(:,n)
                yi = ri(2)
                !Top/bottom surfaces in x
                bintopi(1) = pt(4)*yi**3.d0 + pt(3)*yi**2.d0 + pt(2)*yi + pt(1)
                binboti(1) = pb(4)*yi**3.d0 + pb(3)*yi**2.d0 + pb(2)*yi + pb(1)

                !Use CV function
		        theta_i = dble((heaviside(bintopi(1)-ri(1))-heaviside(binboti(1)-ri(1)))* & 
			               	   (heaviside(bintopi(2)-ri(2))-heaviside(binboti(2)-ri(2)))* & 
			              	   (heaviside(bintopi(3)-ri(3))-heaviside(binboti(3)-ri(3))))

                if (theta_i .eq. 1.d0) then
                    write(pid,'(i10,6f18.9)') n, ri, 0.d0, 0.d0, 0.d0
                else
                    write(pid,'(i10,6f18.9)') n, 0.d0, 0.d0, 0.d0, ri
                endif

            enddo
            close(pid,status='keep')
        endif

        ! - - -Set cluster molecules to be thermostatted - - -

        !First set all thermostatted liquid molecules to unthermostatted
        !do n = 1,np
        !   ! if (tag(n) .eq. thermo) 
        !    tag(n) = free
        !enddo
        !Then set molecules in biggest cluster to thermostatted
        !call thermostat_cluster(self, imaxloc(self%Nlist))

        !Print Biggest cluster
!        call get_Timestep_FileName(iter,'./results/Big_clust',filename)
!        open(unit=1042,file=trim(filename),access='append')
!        do i = 1, self%Nclust
!            if (self%Nlist(i) .ne. maxval(self%Nlist)) then
!                print*, 'skipping', i, 'with ', self%Nlist(i)
!                cycle
!            endif
!            call linklist_printneighbourlist(self, i, 1042)
!        enddo
!        close(1042,status='keep')


    end subroutine get_cluster_properties

    subroutine check_for_cluster_breakup(self)
        use librarymod, only : bubble_sort_r
        use interfaces, only : error_abort
        implicit none

        type(clusterinfo),intent(inout)    :: self

        double precision,dimension(:),allocatable :: cluster_sizes

        ! Sort clusters by size and check if more than one big one!
        allocate(cluster_sizes(self%Nclust))
        cluster_sizes = dble(self%Nlist)
        call bubble_sort_r(cluster_sizes)

        if ((cluster_sizes(1) - cluster_sizes(2))/cluster_sizes(1) .gt. 0.4d0) then

        else
            print*, 'It appears clusters have broken up'
            print'(a,8f10.1,e18.8)', 'CLUSTER DETAILS ', cluster_sizes(1:8), (cluster_sizes(1) - cluster_sizes(2))/cluster_sizes(1)

            !print*, 'It appears clusters have broken up -- writing final state and exiting'
            !Exit Gracefully
	        !call messenger_syncall
	        !call parallel_io_final_state	
	        !call error_abort('Restart file written. Simulation aborted.')
        endif


    end subroutine check_for_cluster_breakup


    !This routine is highly experimental and has not been fully tested

    subroutine thermostat_cluster(self, clustNo)
        use arrays_MD, only : tag
        use computational_constants_MD, only : thermo_tags, thermo
        implicit none

        type(clusterinfo),intent(in)    :: self
        integer, intent(in)             :: clustNo

        integer                         :: n, molno, Nmols
	    type(node), pointer 	        :: old, current

        Nmols = self%Nlist(clustNo)
        if (Nmols .gt. 0) then
            current => self%head(clustNo)%point
            !Loop through all cluster molecules
            do n = 1,Nmols

                !Get molecule number and step to the next link list item
                molno = current%molno
                old => current%next      
                current => old

                !if (any(tag(n) .ne. thermo_tags)) then
                    !print*, n, tag(n)
                    tag(n) = thermo
                !endif

            enddo
        endif

    end subroutine thermostat_cluster


    subroutine build_from_cellandneighbour_lists(self, cell, neighbour, rd, rmols, nmols, skipwalls_)
	    use module_compute_forces, only: cellinfo, neighbrinfo, rj, rij, ri,&
                                         delta_rneighbr, rcutoff, rij2, &
                                         moltype
	    use interfaces, only : error_abort
        implicit none

        type(clusterinfo),intent(inout) :: self
        type(cellinfo),intent(in)       :: cell
        type(neighbrinfo),intent(in)    :: neighbour

        integer,intent(in)              :: nmols
        logical, intent(in),optional    :: skipwalls_
        double precision, intent(in)    :: rd
        double precision, intent(in), dimension(:,:) :: rmols

        logical                         :: skipwalls
	    integer                         :: i, j !Define dummy index
	    integer							:: icell, jcell, kcell, Nchecked
	    integer                         :: cellnp 
	    integer							:: molnoi, molnoj, noneighbrs
        double precision                :: rd2
	    type(node), pointer 	        :: oldi, currenti, noldj,ncurrentj

        if (present(skipwalls_)) then
            skipwalls = skipwalls_
        else
            skipwalls = .false.
        endif

        rd2 = rd**2.d0

        do molnoi = 1, nmols

	        ri = rmols(:,molnoi)         	!Retrieve ri
            if (skipwalls .and. moltype(molnoi) .eq. 2) cycle !Don't include wall molecules

            ! If interface cutoff is less that interaction rcutoff
            ! then we can use the neighbourlist to get molecules in 
            ! interface region (N.B. need to use all interations)
            if (rd .le. rcutoff + delta_rneighbr) then

                noneighbrs = neighbour%Nlist(molnoi)	!Determine number of elements in neighbourlist
	            noldj => neighbour%head(molnoi)%point		!Set old to head of neighbour list
                Nchecked = 0
	            do j = 1,noneighbrs							!Step through all pairs of neighbours i and j

		            molnoj = noldj%molno			        !Number of molecule j

                    !if (molnoj .gt. np) cycle               !Ignore halo values

                    !if (moltype(molnoj) .eq. 2 .and. & 
                    !    .not.( any(tag(molnoj).eq.tether_tags))) stop "ERROR -- moltype not same as tethered"

		            rj(:) = rmols(:,molnoj)			            !Retrieve rj
		            rij(:)= ri(:) - rj(:)   	            !Evaluate distance between particle i and j
		            rij2  = dot_product(rij,rij)            !Square of vector calculated

		            if (rij2 .lt. rd2) then
                        if (skipwalls .and. moltype(molnoj) .eq. 2) then
                            call AddBondedPair(self, molnoi, molnoi)
                        else
                            call AddBondedPair(self, molnoi, molnoj)
                        endif
                        Nchecked = Nchecked + 1
                    endif

		            ncurrentj => noldj
		            noldj => ncurrentj%next !Use pointer in datatype to obtain next item in list
                enddo
                !If no neighbours, add molecule to its own cluster list
                if (Nchecked .eq. 0) then
                    call AddBondedPair(self, molnoi, molnoi)
                endif
            else
                call error_abort("Error in build cluster -- rd must be less than neighbourlist cutoff")
            endif

		    !currenti => oldi
		    !oldi => currenti%next !Use pointer in datatype to obtain next item in list
        enddo

    end subroutine


    subroutine AddBondedPair(self, molnoi, molnoj)
        use linked_list, only : linklist_checkpushneighbr, linklist_merge
        implicit none

        type(clusterinfo),intent(inout)    :: self
        integer, intent(in)                :: molnoi, molnoj

        integer :: nc, nci, ncj, cbig, csmall, m
	    type(node), pointer 	        :: old, current

        !Special case adds one molecule only
        if (molnoi .eq. molnoj) then
            if (self%inclust(molnoi) .eq. 0) then
                self%Nclust = self%Nclust + 1
                nc = self%Nclust
                call linklist_checkpushneighbr(self, nc, molnoi)
                self%inclust(molnoi) = nc
            endif
            return
        endif

        !If molecule i is NOT already in a cluster
        if (self%inclust(molnoi) .eq. 0) then
            !and molecule j is also NOT in a cluster
            if (self%inclust(molnoj) .eq. 0) then
                !Create a new cluster
                self%Nclust = self%Nclust + 1
                nc = self%Nclust
                !Add both molecules to it
                self%inclust(molnoi) = nc
                self%inclust(molnoj) = nc
                !Add to cluster linked lists
                call linklist_checkpushneighbr(self, nc, molnoi)
                call linklist_checkpushneighbr(self, nc, molnoj)

            !But molecule j is in a cluster
            else
                !Get cluster number and add one more to it
                nc = self%inclust(molnoj)
                !Add molecule i to same cluster as j
                self%inclust(molnoi) = nc
                !Add molecule i to cluster linked lists
                call linklist_checkpushneighbr(self, nc, molnoi)

            endif
        !If molecule i is in a cluster
        else
            !But molecule j is NOT in a cluster
            if (self%inclust(molnoj) .eq. 0) then
                !Get cluster number and add one more to it
                nc = self%inclust(molnoi)
                !Add molecule i to same cluster as j
                self%inclust(molnoj) = nc
                !Add molecule i to cluster linked lists
                call linklist_checkpushneighbr(self, nc, molnoj)

            !Molecule j is also in a cluster
            else
                !Load cluster numbers and check if they are the same
                nci = self%inclust(molnoi); ncj = self%inclust(molnoj)
                if (nci .ne. ncj) then
                    !Get biggest and smallest cluster
                    if (self%Nlist(nci) .ge. self%Nlist(ncj)) then
                        cbig = nci; csmall = ncj
                    else
                        cbig = ncj; csmall = nci
                    endif

                    !Change all small cluster references to big
                    current => self%head(csmall)%point
                    do m = 1,self%Nlist(csmall)
                        self%inclust(current%molno) = cbig
                        old => current%next      
                        current => old
                    enddo

                    !Add smaller cluster linked lists to bigger one
                    call linklist_merge(self, keep=cbig, delete=csmall)

                else
                    !If already in the same cluster, nothing to do
                endif
            endif
        endif

    end subroutine AddBondedPair

    !Remove any gaps in the list of clusters
    subroutine CompressClusters(self)
        implicit none

        type(clusterinfo),intent(inout) :: self

        integer                         :: m, j, nc
	    type(node), pointer 	        :: old, current

        !Loop though all clusters
        nc = 0
        do j = 1,self%Nclust
            !For all clusters which are not empty
            if (self%Nlist(j) .gt. 0) then
                !Copy to next sequential array position
                nc = nc + 1
                self%Nlist(nc) = self%Nlist(j)
                self%head(nc)%point => self%head(j)%point
                current => self%head(nc)%point
                !Redefine cluster molecules is included in
                do m = 1,self%Nlist(nc)
                    self%inclust(current%molno) = nc
                    old => current%next      
                    current => old
                enddo
            endif

            self%Nclust = nc
        enddo

    end subroutine CompressClusters


    subroutine cluster_centre_of_mass(self, clustNo, COM)
        use arrays_MD, only : r
        use module_set_parameters, only : mass
        implicit none

        type(clusterinfo),intent(in)    :: self

        integer, intent(in)             :: clustNo
        double precision,dimension(3),intent(out)   :: COM

        integer                         :: m,n,Nmols
        double precision                :: msum
        double precision,dimension(3)   :: MOI
	    type(node), pointer 	        :: old, current

        !For all clusters which are not empty
        msum = 0.d0; MOI = 0.d0
        Nmols = self%Nlist(clustNo)
        if (Nmols .gt. 0) then
            current => self%head(clustNo)%point
            !Loop through all cluster molecules
            do n = 1,Nmols
                m = mass(n)
                MOI = MOI + m*r(:,current%molno)
                msum = msum + m
                old => current%next      
                current => old
            enddo
        endif
        COM = MOI/msum
    end subroutine cluster_centre_of_mass

    subroutine cluster_global_extents(self, clustNo, extents)
        use arrays_MD, only : r
        use module_set_parameters, only : mass
        use computational_constants_MD, only : halfdomain
        implicit none

        type(clusterinfo),intent(in)    :: self

        integer, intent(in)             :: clustNo
        double precision,dimension(6),intent(out)   :: extents

        integer                         :: ixyz,n,Nmols
	    type(node), pointer 	        :: old, current

        !For all clusters which are not empty
        extents = 0.d0
        Nmols = self%Nlist(clustNo)
        if (Nmols .gt. 0) then
            current => self%head(clustNo)%point
            !Loop through all cluster molecules         
            do n = 1,Nmols
                do ixyz = 1,3
                    extents(ixyz)   = max(r(ixyz,current%molno),extents(ixyz))
                    extents(ixyz+3) = min(r(ixyz,current%molno),extents(ixyz+3))
                enddo
                old => current%next      
                current => old
            enddo
        endif
    
        do ixyz = 1,3
            !print*, ixyz, extents(ixyz), extents(ixyz+3), halfdomain(ixyz)
            if (extents(ixyz) .gt. halfdomain(ixyz)) extents(ixyz) = halfdomain(ixyz)
            if (extents(ixyz+3) .lt. -halfdomain(ixyz) ) extents(ixyz+3) = -halfdomain(ixyz)
        enddo

    end subroutine cluster_global_extents

    !Generate extents for the 2D surface on a cell by cell basis
    subroutine cluster_extents_grid(self, clustNo, dir, resolution, & 
                                    extents, maxmols, debug_outfile)
        use arrays_MD, only : r
        use module_set_parameters, only : mass
        use computational_constants_MD, only : halfdomain, iter
        use librarymod, only : get_new_fileunit,get_Timestep_FileName
        implicit none

        type(clusterinfo),intent(in)    :: self

        integer, intent(in)             :: clustNo, dir, resolution
        double precision,dimension(:,:),allocatable,intent(out)   :: extents
        double precision,dimension(:,:),allocatable,intent(out),optional   :: maxmols
        character(*),intent(in),optional :: debug_outfile

        character(32)                   :: filename
        logical                         :: print_debug = .false.
        integer                         :: i,j,ixyz,jxyz,kxyz,n,molno,Nmols
        integer                         :: jcell, kcell, surftodim
        double precision,dimension(6)   :: global_extents
        double precision,dimension(2)   :: cellsidelength, clusterwidth
	    type(node), pointer 	        :: old, current

        if (present(debug_outfile)) print_debug = .true.

        !Allocate array of extents
        allocate(extents(resolution,resolution))
        if (present(maxmols)) allocate(maxmols(resolution,resolution))
        if (dir .le. 3) then
            surftodim = 0
        elseif (dir .gt. 3) then
            surftodim = 3
        endif
        !Get directional index and orthogonal directions
        ixyz = dir-surftodim
        kxyz = mod(ixyz+1,3)+1 
        jxyz = 6 - ixyz - kxyz

        !Get global limits of cluster and cell spacing in orthogonal directions
        call cluster_global_extents(self, clustNo, global_extents)
        clusterwidth(1) = (global_extents(jxyz) - global_extents(jxyz+3))
        clusterwidth(2) = (global_extents(kxyz) - global_extents(kxyz+3))
        cellsidelength(1) = clusterwidth(1)/dble(resolution)
        cellsidelength(2) = clusterwidth(2)/dble(resolution)

        !For all clusters which are not empty
        if (dir .le. 3) then
            extents(:,:) = global_extents(ixyz+3)
        elseif (dir .gt. 3) then
            extents(:,:) = global_extents(ixyz)
        endif

        Nmols = self%Nlist(clustNo)
        if (Nmols .gt. 0) then
            current => self%head(clustNo)%point
            !Loop through all cluster molecules
            do n = 1,Nmols
                !Get molecule number and step to the next link list item
                molno = current%molno
                old => current%next      
                current => old

                !Get cell indices in orthogonal direction
	            jcell = ceiling((r(jxyz,molno)-global_extents(jxyz+3))/cellsidelength(1))
                kcell = ceiling((r(kxyz,molno)-global_extents(kxyz+3))/cellsidelength(2))

                !Ignore halo molecules
                if (jcell .lt. 1) cycle
                if (kcell .lt. 1) cycle
                if (jcell .gt. resolution) cycle
                if (kcell .gt. resolution) cycle

                !Get extents in requested direction for current cell
                if (dir .le. 3) then
                    extents(jcell,kcell) = max(r(ixyz,molno),extents(jcell,kcell))
                elseif (dir .gt. 3) then  
                    extents(jcell,kcell) = min(r(ixyz,molno),extents(jcell,kcell))
                endif
                if (present(maxmols) .and. & 
                    r(ixyz,molno) .eq. extents(jcell,kcell)) then
                    maxmols(jcell,kcell) = molno            
                endif

            enddo
        endif

        if (print_debug) then
            call get_Timestep_FileName(iter,debug_outfile,filename)
            n = get_new_fileunit()
            open(unit=n,file=trim(filename),status='replace')
            do i = 1,resolution
            do j = 1,resolution
                write(n,'(2i6,3f10.5)') i,j,extents(i,j),& 
                                 (i-0.5d0)*(global_extents(2)-global_extents(5)) & 
                                  /dble(resolution)+global_extents(5), &
                                 (j-0.5d0)*(global_extents(3)-global_extents(6)) & 
                                  /dble(resolution)+global_extents(6)
            enddo
            enddo
            close(n,status='keep')
        endif
    
    end subroutine cluster_extents_grid


    subroutine cluster_outer_mols(self, clustNo, tolerence, dir, rmols, extents, debug_outfile)
        use arrays_MD, only : r
        use module_set_parameters, only : mass
        use computational_constants_MD, only : halfdomain, iter
	    use interfaces, only : error_abort
        use librarymod, only : get_new_fileunit, get_Timestep_FileName
        implicit none

        type(clusterinfo),intent(in)    :: self

        integer, intent(in)             :: clustNo, dir
        double precision,intent(in)     :: tolerence
        double precision,dimension(:,:),intent(in),optional  :: extents
        double precision,dimension(:,:),allocatable,intent(out) :: rmols
        character(*),intent(in),optional :: debug_outfile

        character(32)                   :: filename
        logical                         :: print_debug = .false.

        integer                         :: ixyz,jxyz,kxyz,i,j,n,m, Nmols, molno
        integer                         :: jcell, kcell, surftodim
        double precision,dimension(6)               :: global_extents
        double precision,dimension(:,:),allocatable :: rtemp, extents_, molband
        double precision,dimension(2)   :: cellsidelength, clusterwidth
	    type(node), pointer 	        :: old, current

        if (present(debug_outfile)) print_debug = .true.

        !Get directional index and orthogonal directions
        if (dir .le. 3) then
            surftodim = 0
        elseif (dir .gt. 3 .and. dir .le. 6) then
            surftodim = 3
        else
            call error_abort("Error in cluster_outer_mols -- dir should be between 1 and 6")
        endif
        ixyz = dir-surftodim
        kxyz = mod(ixyz+1,3)+1 
        jxyz = 6 - ixyz - kxyz

        !If no extents supplies, use global cluster values
        call cluster_global_extents(self, clustNo, global_extents)
        clusterwidth(1) = (global_extents(jxyz) - global_extents(jxyz+3))
        clusterwidth(2) = (global_extents(kxyz) - global_extents(kxyz+3))
        if (present(extents)) then
            allocate(extents_(size(extents,1),size(extents,2)))
            allocate(molband(size(extents,1),size(extents,2)))
            extents_ = extents
            cellsidelength(1) = clusterwidth(1)/dble(size(extents_,1))
            cellsidelength(2) = clusterwidth(2)/dble(size(extents_,2))
        else
            allocate(extents_(1,1))
            allocate(molband(1,1))
            extents_(1,1) = global_extents(dir)
            cellsidelength(:) = clusterwidth(:)
        endif

        !Get band of molecules within tolerence 
        if (dir .le. 3) then
            molband = extents_ - tolerence
        elseif (dir .gt. 3) then
            molband = extents_ + tolerence
        endif

        Nmols = self%Nlist(clustNo)
        allocate(rtemp(3,Nmols))
        m = 0
        !For all clusters which are not empty
        if (Nmols .gt. 0) then
            current => self%head(clustNo)%point
            !Loop through all cluster molecules
            do n = 1,Nmols
                !Get molecule number and step to the next link list item
                molno = current%molno
                old => current%next      
                current => old

                !Get cell indices in orthogonal direction
	            jcell = ceiling((r(jxyz,molno)-global_extents(jxyz+3))/cellsidelength(1))
                kcell = ceiling((r(kxyz,molno)-global_extents(kxyz+3))/cellsidelength(2))

                !print'(2i6,6f10.5)', jcell, kcell, r(jxyz,molno), r(kxyz,molno), cellsidelength, clusterwidth

                !Ignore halo molecules
                if (jcell .lt. 1) cycle
                if (kcell .lt. 1) cycle
                if (jcell .gt. size(extents_,1)) cycle
                if (kcell .gt. size(extents_,2)) cycle

                !Check if in top band
                if (dir .le. 3) then
                    if (r(ixyz,molno) .lt. extents_(jcell,kcell) .and. & 
                        r(ixyz,molno) .gt. molband(jcell,kcell)) then
                        m = m + 1
                        rtemp(:,m) = r(:,molno)
                    endif
                elseif (dir .gt. 3) then  
                !Check if in bottom band
                    if (r(ixyz,molno) .gt. extents_(jcell,kcell) .and. & 
                        r(ixyz,molno) .lt. molband(jcell,kcell)) then
                        m = m + 1
                        rtemp(:,m) = r(:,molno)
                    endif
                endif

            enddo
        endif

        !Copy total array to array size of outer molecules
        allocate(rmols(3,m))
        rmols = rtemp(:,1:m)

        if (print_debug) then
            call get_Timestep_FileName(iter,debug_outfile,filename)
            n = get_new_fileunit()
            open(unit=n,file=trim(filename),status='replace')
            do i =1,size(rmols,2)
                write(n,'(i6,3f10.5)') i, rmols(:,i)
            enddo
            close(n,status='keep')
        endif

    end subroutine cluster_outer_mols


!    subroutine build_debug_clusters(self, rd)
!        use librarymod, only : get_Timestep_FileName, get_new_fileunit
!	    use module_compute_forces, only : iter, np, r, halfdomain, cellsidelength, nh, ncells, rneighbr2
!        use linked_list, only : build_cell_and_neighbourlist_using_debug_positions, cellinfo,neighbrinfo
!	    implicit none

!        type(clusterinfo),intent(inout) :: self
!        double precision, intent(in)    :: rd

!	    integer		:: i, cellnp, n, m, testmols, molcount, clustno, molnoi, molnoj, adjacentcellnp, noneighbrs, j
!	    integer		:: icell, jcell, kcell, icellshift, kcellshift, jcellshift, fileunit
!        double precision                :: rd2, rij2
!        double precision,dimension(3)   :: ri, rj, rij
!        double precision,dimension(4)   :: rand
!        double precision, dimension(:,:),allocatable   :: rdebug
!        character(20)                    :: filename
!        type(cellinfo)                  :: celldebug
!        type(neighbrinfo)               :: neighbourdebug
!	    type(node), pointer	            :: old, current



!        rd2 = rd**2.d0
!        testmols = 200

!        allocate(rdebug(3,testmols))
!        rdebug = 0.d0

!        do n = 1,testmols
!            call random_number(rand)
!            if (rand(1) .le. 1.d0/4.d0) then
!                rdebug(:,n) = (/ 0.d0, 5.d0, 0.d0 /) + (rand(2:4)-0.5d0) * 2.d0
!            elseif (rand(1) .gt. 1.d0/4.d0 .and. rand(1) .le. 2.d0/4.d0) then
!                rdebug(:,n) = (/  6.d0, -6.d0, 0.d0 /) + (rand(2:4)-0.5d0) * 2.d0
!            elseif (rand(1) .gt. 2.d0/4.d0 .and. rand(1) .le. 3.d0/4.d0) then
!                rdebug(:,n) = (/ -7.d0, -7.d0, 0.d0 /) + (rand(2:4)-0.5d0) * 2.d0
!            else
!                rdebug(:,n) = (/ -7.d0, -7.d0, 0.d0 /) + (/ (rand(2)) * 10.d0, (rand(3)) * 14.d0 , 0.d0 /)
!            endif
!            !print'(i6,3f10.5)', n, rdebug(:,n)
!        enddo


!        !Try an all pairs solution
!!        do n = 1,testmols
!!            ri(:) = rdebug(:,n)
!!            do m = n,testmols
!!                if (n .eq. m) cycle
!!                rj(:) = rdebug(:,m)			            !Retrieve rj
!!			    rij(:)= ri(:) - rj(:)   	            !Evaluate distance between particle i and j
!!			    rij2  = dot_product(rij,rij)            !Square of vector calculated
!!	            if (rij2 .lt. rd2) then
!!                    call AddBondedPair(self, n, m)
!!                endif
!!            enddo
!!        enddo

!        !Or build all cell and neighbour lists
!        call build_cell_and_neighbourlist_using_debug_positions(rdebug, celldebug, neighbourdebug)
!        !Then call original routine to build
!        call build_from_cellandneighbour_lists(self, celldebug, neighbourdebug, rd, rdebug, testmols, skipwalls_=.false.)

!        !And print test of neighlists
!!        do molno = 1, testmols

!!	        noneighbrs = neighbourdebug%Nlist(molno)  	!Determine number of elements in neighbourlist
!!	        old => neighbourdebug%head(molno)%point		!Set old to head of neighbour list

!!	        if(noneighbrs == 0) print'(i6, a,3f10.5)', molno, 'has 0 neighbours', rdebug(:,molno)

!!	        current => old ! make current point to head of list
!!	        do j=1,noneighbrs
!!                
!!		        !print*, 'more items in linked list?: ', associated(old%next)
!!          		print'(i6, 2(a,i8),3f10.5)', j, ' Linklist print called for molecule ', molno,' j = ', current%molno, rdebug(:,current%molno)
!!		        if (associated(old%next) .eqv. .true. ) then !Exit if null
!!			        old => current%next ! Use pointer in datatype to obtain next item in list
!!			        current => old      ! make current point to old - move alone one
!!		        endif
!!	        enddo

!!	        nullify(current)                    !Nullify current as no longer required
!!	        nullify(old)                        !Nullify old as no longer required

!!        enddo

!        !PRINT CLUSTER DEBUGGING INFO
!        call CompressClusters(self)

!    	fileunit = get_new_fileunit()
!        call get_Timestep_FileName(iter,'./Big_clust',filename)
!        open(unit=fileunit,file=trim(filename),access='append')

!        molcount = 0
!        do clustno = 1, self%Nclust

!	        noneighbrs = self%Nlist(clustno)  	!Determine number of elements in neighbourlist
!	        old => self%head(clustno)%point		!Set old to head of neighbour list

!	        if(noneighbrs == 0) print*, clustno, 'cluster has no elements'

!	        current => old ! make current point to head of list
!	        do j=1,noneighbrs
!                molcount =  molcount + 1
!          		write(fileunit,'(2i6, 3(a,i8),3f10.5)'), molcount, j, ' of ', self%Nlist(clustno), & 
!                                             ' Linklist for cluster ', clustno,' j = ', & 
!                                             current%molno, rdebug(:,current%molno)
!		        if (associated(old%next) .eqv. .true. ) then !Exit if null
!			        old => current%next ! Use pointer in datatype to obtain next item in list
!			        current => old      ! make current point to old - move alone one
!		        endif
!	        enddo

!	        nullify(current)                    !Nullify current as no longer required
!	        nullify(old)                        !Nullify old as no longer required

!        enddo

!        close(fileunit,status='keep')

!     

!    end subroutine build_debug_clusters

    subroutine destroy_clusters(self)
        use linked_list, only : linklist_deallocate_cluster
        implicit none

        type(clusterinfo),intent(inout)    :: self

        call linklist_deallocate_cluster(self)

    end subroutine destroy_clusters


end module cluster_analysis

subroutine get_interface_from_clusters()
    use cluster_analysis, only : build_clusters, & 
                                 get_cluster_properties, & 
                                 destroy_clusters
    use linked_list, only : cluster
    implicit none

    double precision    :: rd

    rd = 1.5d0
    call build_clusters(cluster, rd)
    call get_cluster_properties(cluster, rd)
    call destroy_clusters(cluster)

end subroutine get_interface_from_clusters


!This routine uses the mass averaging to obtain the gas/liquid interface

module sl_interface_mod

contains

    !This routine uses the mass averaging to obtain the gas/liquid interface

    !Requires volume_mass to be up to date and divided by (volume*Nmass_ave)
    ! (OR use mass bins and multiply liquiddensity/massdensity by (volume*Nmass_ave)

    subroutine get_fluid_liquid_interface_cells(density_bins, & 
                                                soliddensity, &
                                                liquiddensity, & 
                                                gasdensity, &
                                                interfacecells, &
                                                liquidcells)
        use computational_constants_MD, only : iter, binspercell, Nmass_ave
        use physical_constants_MD, only : globalnp
        implicit none


        double precision, intent(in)                               :: liquiddensity,gasdensity,soliddensity
        double precision, dimension(:,:,:),allocatable,intent(in)  :: density_bins

        integer, dimension(:,:),allocatable,optional,intent(out)   :: interfacecells
        integer, dimension(:,:),allocatable,optional,intent(out)   :: liquidcells

	    integer                         :: ixyz, i, j, k, is, js, ks, ncount
        double precision                :: averagegldensity,averagesldensity
        double precision,dimension(3)   :: cellsperbin
        double precision, dimension(:,:),allocatable               :: interfacelist_temp, liquidcells_temp
        logical, dimension(:,:,:),allocatable                      :: liquidbins, gasbins, interfacebin


        !Sanity check
        if (sum(density_bins(2:size(density_bins,1)-1, & 
                             2:size(density_bins,2)-1, &
                             2:size(density_bins,3)-1)) .ne. &
                            (Nmass_ave-1)*globalnp) then
            print*, 'Warning in interface check - number of molecules in mass averaged bins is greater than total'
            print*, 'e.g.', globalnp, &
                    sum(    density_bins                   &
                            (                              &
                               2:size(density_bins,1)-1,   &
                               2:size(density_bins,2)-1,   &
                               2:size(density_bins,3)-1    &
                            )                              &
                       )/(Nmass_ave-1)
        endif
        print*, 'Densities', &
                sum&
                (&
                    density_bins &
                    (&
                        2:size(density_bins,1)-1,&
                        2:size(density_bins,2)-1,&
                        2:size(density_bins,3)-1 &
                    )&
                ),&
                product&
                (&
                    shape&
                    (&
                        density_bins&
                        (&
                            2:size(density_bins,1)-1,&
                            2:size(density_bins,2)-1,&
                            2:size(density_bins,3)-1&
                        )&
                    )&
                ),&
                liquiddensity, gasdensity

	    cellsperbin = 1.d0/binspercell !ceiling(ncells(1)/dble(nbins(1)))

        averagesldensity = 0.5d0 * (liquiddensity + soliddensity)
        averagegldensity = 0.5d0 * (liquiddensity + gasdensity)
        allocate(liquidbins(size(density_bins,1), &
                            size(density_bins,2), &
                            size(density_bins,3)))
        allocate(gasbins   (size(density_bins,1), &
                            size(density_bins,2), &
                            size(density_bins,3)))
        liquidbins = .false.; gasbins = .false.
        do i = 1,size(density_bins,1)
        do j = 1,size(density_bins,2)
        do k = 1,size(density_bins,3)
            !print'(4i7,f10.5,2l)', i,j,k,density_bins(i,j,k), averagegldensity,density_bins(i,j,k) .gt. averagegldensity,density_bins(i,j,k) .le. averagegldensity

            if (density_bins(i,j,k) .gt. averagesldensity) then
                open(unit=2015,file='./solidcells',access='append')
                write(2015,'(5i5,i8,f10.5)') iter,ncount,i,j,k,density_bins(i,j,k),averagesldensity
                close(2015,status='keep')         
            elseif (density_bins(i,j,k) .gt. averagegldensity) then
                liquidbins(i,j,k) = .true.
                open(unit=2010,file='./liquidcells',access='append')
                write(2010,'(5i5,i8,f10.5)') iter,ncount,i,j,k,density_bins(i,j,k),averagegldensity
                close(2010,status='keep')
            elseif (density_bins(i,j,k) .le. averagegldensity) then
                gasbins(i,j,k) = .true.
                open(unit=2020,file='./gascells',access='append')
                write(2020,'(5i5,i8,f10.5)') iter,ncount,i,j,k,density_bins(i,j,k),averagegldensity
                close(2020,status='keep')
            endif
        enddo
        enddo
        enddo

        !If list of liquid cells is requested, get any cell which is true 
        if (present(liquidcells)) then
            allocate(liquidcells_temp(3,product(shape(density_bins))))
            do i = 1,size(density_bins,1)
            do j = 1,size(density_bins,2)
            do k = 1,size(density_bins,3)
                if (liquidbins(i,j,k)) then
                    ncount = ncount + 1
                    liquidcells_temp(:,ncount) = (/i,j,k/)
                endif
            enddo
            enddo
            enddo
            !Reduce array to size of count
            allocate(liquidcells(3,ncount))
            liquidcells = liquidcells_temp(:,1:ncount)
            deallocate(liquidcells_temp)
        endif

        if (present(interfacecells)) then

            !If liquid bin is next to a gas bin then must be an interface position
            allocate(interfacebin(size(density_bins,1), &
                                  size(density_bins,2), &
                                  size(density_bins,3)))
            interfacebin = .false.; 
            ncount = 0
            allocate(interfacelist_temp(3,product(shape(density_bins))))
            do i = 2,size(density_bins,1)-1
            do j = 2,size(density_bins,2)-1
            do k = 2,size(density_bins,3)-1
	            do ks = -1,1
	            do js = -1,1
	            do is = -1,1

				    !Prevents out of range values in i
				    if (i+is .lt.           2           ) cycle
				    if (i+is .gt. size(density_bins,1)-1) cycle
				    !Prevents out of range values in j
				    if (j+js .lt.           2           ) cycle
				    if (j+js .gt. size(density_bins,2)-1) cycle
				    !Prevents out of range values in k
				    if (k+ks .lt.           2           ) cycle
				    if (k+ks .gt. size(density_bins,3)-1) cycle

                    if (gasbins(i,j,k) .and. liquidbins(i+is,j+js,k+ks) .and. (.not. interfacebin(i,j,k))) then
                        interfacebin(i,j,k) = .true.
                        ncount = ncount + 1
                        interfacelist_temp(:,ncount) = (/i,j,k/)
                        open(unit=90210,file='./interfacecells',access='append')
                        write(90210,'(6i12,f10.5)') iter,ncount,i,j,k,density_bins(i,j,k),averagegldensity
                        close(90210,status='keep')
                    endif

                    if (gasbins(i+is,j+js,k+ks) .and. liquidbins(i,j,k)  .and. (.not. interfacebin(i,j,k))) then
                        interfacebin(i,j,k) = .true.
                        ncount = ncount + 1
                        interfacelist_temp(:,ncount) = (/i,j,k/)
                        open(unit=90210,file='./interfacecells',access='append')
                        write(90210,'(6i12,f10.5)') iter,ncount,i,j,k,density_bins(i,j,k),averagegldensity
                        close(90210,status='keep')
                    endif
                enddo
                enddo
                enddo
            enddo
            enddo
            enddo

            !Reduce array to size of count
            allocate(interfacecells(3,ncount))
            interfacecells = interfacelist_temp(:,1:ncount)
            deallocate(interfacelist_temp)
        endif

    end subroutine get_fluid_liquid_interface_cells

    !===============================
    !Inputs:
    !       A -- 3 x N array of values
    !Outputs:
    !       Omega -- 3 x 3 covariance matrix
    !

    subroutine get_covariance(A,omega,normalise)
        use librarymod, only : outerprod
        implicit none

        logical,intent(in),optional :: normalise
        real(kind(0.d0)),dimension(:,:),allocatable,intent(in) :: A

        real(kind(0.d0)),dimension(3,3),intent(out) :: omega

        integer                         :: j
        real(kind(0.d0)),dimension(3) :: rj

        omega = 0.d0
        do j = 1,size(A,2)
            rj = A(:,j)
            omega = omega + outerprod(rj,rj)
        enddo

        if (present(normalise)) then
            if (normalise) omega = omega/size(A,2)
        endif

    end subroutine get_covariance

    !Input -- molecules within rd of molecule i

    subroutine get_surface_at_ri(ri, rarray, surface)
        use librarymod, only : get_eigenvec3x3
        !TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
        use computational_constants_MD, only : iter
        !TEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMPTEMP
        implicit none

        real(kind(0.d0)),dimension(3),intent(in) ::  ri
        real(kind(0.d0)),dimension(:,:),allocatable,intent(in) ::  rarray
        real(kind(0.d0)),dimension(3,3),intent(out) :: surface 

        integer :: i, mineig
        real(kind(0.d0)),dimension(3)   :: rave,eigval
        real(kind(0.d0)),dimension(3,3) :: omega, eigvec

        real(kind(0.d0)),dimension(:,:),allocatable :: rcent

        if (size(rarray,2) .eq. 0) then
            print'(a,3f10.5,a)', 'Molecule at ', ri ,' has no neighbours!'
            return
        endif

        !Get average position of all molecules
        rave = sum(rarray,2)/size(rarray,2)

        !Subtract average from all molecules
        allocate(rcent(size(rarray,1),size(rarray,2)))
        do i = 1,size(rarray,2)
            rcent(:,i) = rarray(:,i) - rave(:)
        enddo

        !Get covariance (omega) of rcert
        call get_covariance(rcent,omega,normalise=.true.)

        !Get eigenvalues and eigenvectors of covariance (omega)
        call  get_eigenvec3x3(omega,eigvec,eigval)

        mineig = minloc(eigval,1)
        surface(:,1) = eigvec(:,mineig)
        surface(:,2) = eigvec(:,mod(mineig,3)+1)
        surface(:,3) = eigvec(:,mod(mineig+1,3)+1)

        open(unit=1984,file='./eigenvalues',access='append')
        write(1984,'(i8,6f10.4,3f20.7)'), iter, ri, rave(:), eigval
        close(1984,status='keep')

 !dot_product(surface(:,1),surface(:,2)), dot_product(surface(:,1),surface(:,3)),surface

        !do i = 1,3
        !    print'(a,2i5,4f15.4,6f10.4)','Eigenstuff', i,minloc(eigval,1), omega(i,:),eigval(i),eigvec(i,:),surface(i,:)
        !enddo

    end subroutine get_surface_at_ri


    subroutine get_molecules_within_rd(celllist, rd, rarray)
	    use module_compute_forces
	    use librarymod, only: outerprod
        use computational_constants_MD, only : iter
	    use interfaces, only : error_abort
	    implicit none

        real(kind(0.d0)), intent(in)                     :: rd
        integer, dimension(:,:),allocatable,intent(in)   :: celllist

        real(kind(0.d0)), dimension(:,:),allocatable  :: rarray

	    integer                         :: i, j, n, ixyz !Define dummy index
	    integer							:: icell, jcell, kcell, ncount
	    integer                         :: icellshift, jcellshift, kcellshift
	    integer                         :: cellnp, adjacentcellnp 
	    integer							:: molnoi, molnoj, noneighbrs, cellshifts
	    integer							:: icellmin,jcellmin,kcellmin,icellmax,jcellmax,kcellmax


        real(kind(0.d0))                :: rd2
        real(kind(0.d0)), dimension(:,:),allocatable :: rneigh
        real(kind(0.d0)), dimension(:,:,:),allocatable :: cellsurface, surfacei

	    type(node), pointer 	        :: oldi, currenti, oldj, currentj, noldj,ncurrentj

        rd2 = rd**2.d0
        if (force_list .ne. 2) then
            call error_abort("Error in get_molecules_within_rc -- full "//&
                             "neightbour list should be used with interface tracking")
        end if

        allocate(cellsurface(3,3,size(celllist,2)))
        do n = 1,size(celllist,2)

            icell = celllist(1,n)
            jcell = celllist(2,n)
            kcell = celllist(3,n)

		    if ( icell .lt. 2 .or. icell .gt. ncells(1)+1 .or. &
			     jcell .lt. 2 .or. jcell .gt. ncells(2)+1 .or. &
			     kcell .lt. 2 .or. kcell .gt. ncells(3)+1      ) then
			     print*, 'Warning - celllist is in halo'
			    return
		    end if

		    cellnp = cell%cellnp(icell,jcell,kcell)
		    oldi => cell%head(icell,jcell,kcell)%point !Set oldi to first molecule in list

            !print'(5i7,2f10.5)', n, icell, jcell, kcell, cellnp, rd, rcutoff + delta_rneighbr

            allocate(surfacei(3,3,cellnp))
		    do i = 1,cellnp					!Step through each particle in list 
			    molnoi = oldi%molno 	 	!Number of molecule
			    ri = r(:,molnoi)         	!Retrieve ri

                ! If interface cutoff is less that interaction rcutoff
                ! then we can use the neighbourlist to get molecules in 
                ! interface region (N.B. need to use all interations)
                if (rd .le. rcutoff + delta_rneighbr) then

	                noneighbrs = neighbour%Nlist(molnoi)	!Determine number of elements in neighbourlist
		            noldj => neighbour%head(molnoi)%point		!Set old to head of neighbour list
                    allocate(rneigh(3,noneighbrs))
                    ncount = 0
		            do j = 1,noneighbrs							!Step through all pairs of neighbours i and j

			            molnoj = noldj%molno			        !Number of molecule j
			            rj(:) = r(:,molnoj)			            !Retrieve rj
			            rij(:)= ri(:) - rj(:)   	            !Evaluate distance between particle i and j
			            rij2  = dot_product(rij,rij)            !Square of vector calculated

			            if (rij2 < rd2) then
                            write(451,'(4i8,6f11.5)') iter,icell,jcell,kcell,ri(:),rj(:)
                            ncount = ncount + 1
                            rneigh(:,ncount) = rj(:)
                        endif

			            ncurrentj => noldj
			            noldj => ncurrentj%next !Use pointer in datatype to obtain next item in list
                    enddo
                    allocate(rarray(3,ncount))
                    rarray = rneigh(:,1:ncount)
                    deallocate(rneigh)
                    call get_surface_at_ri(ri, rarray, surfacei(:,:,i))
                    write(452,'(i5,12f10.4)'), i, ri(:), surfacei(:,:,i)
                    deallocate(rarray)

                ! If the interface cutoff is greater than rcutoff
                ! then we can't use neighbour list and need to loop 
                ! over a greater number of adjacent cells
                elseif (rd .gt. rcutoff + delta_rneighbr) then

                    stop "May work, probably won't so check get_molecules_within_rd for rd > rcutoff + delta_rneighbr"

!                    cellshifts = ceiling(rd/(rcutoff + delta_rneighbr))

!			        do kcellshift = -cellshifts,cellshifts
!			        do jcellshift = -cellshifts,cellshifts
!			        do icellshift = -cellshifts,cellshifts

!				        !Prevents out of range values in i
!				        if (icell+icellshift .lt. icellmin) cycle
!				        if (icell+icellshift .gt. icellmax) cycle
!				        !Prevents out of range values in j
!				        if (jcell+jcellshift .lt. jcellmin) cycle
!				        if (jcell+jcellshift .gt. jcellmax) cycle
!				        !Prevents out of range values in k
!				        if (kcell+kcellshift .lt. kcellmin) cycle
!				        if (kcell+kcellshift .gt. kcellmax) cycle

!				        oldj => cell%head(icell+icellshift,jcell+jcellshift,kcell+kcellshift)%point
!				        adjacentcellnp = cell%cellnp(icell+icellshift,jcell+jcellshift,kcell+kcellshift)

!				        do j = 1,adjacentcellnp      !Step through all j for each i

!					        molnoj = oldj%molno 	 !Number of molecule
!					        rj = r(:,molnoj)         !Retrieve rj

!					        if(molnoi==molnoj) cycle !Check to prevent interaction with self

!					        rij2=0                   !Set rij^2 to zero
!					        rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j
!					        rij2 = dot_product(rij,rij)	!Square of vector calculated

!					        if (rij2 < rd2) then
!                                write(451,'(7i5,6f10.5)') iter, icell,jcell,kcell, & 
!                                                         icell+icellshift,jcell+jcellshift, & 
!                                                         kcell+kcellshift, ri(:),rj(:)  
!                            endif

!					        currentj => oldj
!					        oldj => currentj%next    !Use pointer in datatype to obtain next item in list

!				        enddo
!			        enddo
!			        enddo
!			        enddo
                endif

		        currenti => oldi
		        oldi => currenti%next !Use pointer in datatype to obtain next item in list

            enddo
            !Sum up cellnp molecules in cell
            if (cellnp .ne. 0) then
                cellsurface(:,:,n) = sum(surfacei,3)/dble(cellnp)
            endif

            deallocate(surfacei)

            !print'(a,3i8,9f10.5)','cell sum', icell,jcell,kcell,cellsurface
	    enddo

	    nullify(oldi)      	!Nullify as no longer required
	    nullify(oldj)      	!Nullify as no longer required
	    nullify(currenti)      	!Nullify as no longer required
	    nullify(currentj)      	!Nullify as no longer required

    end subroutine get_molecules_within_rd



end module sl_interface_mod

!Calculate interface location
subroutine sl_interface_from_binaverage()
    use sl_interface_mod
    use interfaces, only : error_abort
    use physical_constants_MD, only : density
    use calculated_properties_MD, only : nbins, volume_mass
	use computational_constants_MD, only: domain, iter, Nmass_ave,& 
                                          tplot, & 
                                          liquid_density, gas_density
    implicit none

    real(kind(0.d0))   :: binvolume, input_soliddensity, input_liquiddensity, input_gasdensity, rd
    integer,dimension(:,:),allocatable    :: interfacecells
    real(kind(0.d0)),dimension(:,:),allocatable    :: rarray

    if (Nmass_ave .eq. 0) call error_abort("Error in sl_interface_from_binaverage -- Interface checking requires mass binning")
    if (mod(iter/tplot+1,(Nmass_ave)) .eq. 0) then
        !allocate(density_bins(size(volume_mass,1),size(volume_mass,2),size(volume_mass,3)))
    	binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
        !density_bins = volume_mass/(binvolume*Nmass_ave)
        input_soliddensity = density*binvolume*(Nmass_ave-1)
        input_liquiddensity = liquid_density*binvolume*(Nmass_ave-1)
        input_gasdensity = gas_density*binvolume*(Nmass_ave-1)
        !print*, 'interface output', mod(iter/tplot,Nmass_ave),iter/tplot,iter,tplot, & 
        !                            Nmass_ave,input_gasdensity,gas_density, & 
        !                            input_liquiddensity,liquid_density,binvolume
        call get_fluid_liquid_interface_cells(volume_mass, &
                                              input_soliddensity, &  
                                              input_liquiddensity, & 
                                              input_gasdensity, &
                                              interfacecells)
        rd = 1.5d0
        call get_molecules_within_rd(interfacecells,rd,rarray)
    else
        !print*, 'NO output', mod(iter/tplot,Nmass_ave), iter,tplot,Nmass_ave
    endif

end subroutine sl_interface_from_binaverage



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
