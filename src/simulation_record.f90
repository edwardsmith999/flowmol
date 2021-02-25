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
!   Written by Edward ( ͡° ͜ʖ ͡°)﻿ Smith unless otherwise stated
!===================================================================================
module module_record

	use interfaces
	use physical_constants_MD
	use computational_constants_MD
	use arrays_MD
	use calculated_properties_MD
	use polymer_info_MD
    use intrinsic_module, only : intrinsic_surface_real, intrinsic_surface_bilinear

#if __INTEL_COMPILER > 1200
    use boundary_MD, only: bforce_pdf_measure
	!use module_set_parameters, only : velPDF, velPDFMB, velPDF_array

#endif
	real(kind(0.d0)) :: vel
    integer :: qm, qu
	class(intrinsic_surface_real), pointer :: ISR, ISR_mdt
	type(intrinsic_surface_real), target		:: ISR_r, ISR_mdt_r	! declare an instance
	type(intrinsic_surface_bilinear), target		:: ISR_b, ISR_mdt_b	! declare an instance
	!real(kind(0.d0)), dimension(:,:,:), allocatable :: q_vectors
	real(kind(0.d0)), dimension(:,:,:,:), allocatable :: Abilinear


	abstract interface
		function fn_get_bin(r)
	        integer,dimension(3) 					 :: fn_get_bin
            real(kind(0.d0)),intent(in),dimension(3) :: r
		end function fn_get_bin
	end interface

	abstract interface
		function fn_get_bin_molno(n)
	        integer,dimension(3) :: fn_get_bin_molno
            integer,intent(in)   :: n
		end function fn_get_bin_molno
	end interface

   procedure (fn_get_bin),   pointer :: get_bin   => null ()
   procedure (fn_get_bin_molno),   pointer :: get_bin_molno   => null ()

contains

    function bin_from_integer_division(r) result(bin)
        use computational_constants_MD, only : halfdomain, nhb
        use calculated_properties_MD, only : binsize
        implicit none

        real(kind(0.d0)),intent(in),dimension(3) :: r
	    integer,dimension(3) 					 :: bin

        bin = ceiling((r+halfdomain)/binsize)+nhb

    end function bin_from_integer_division

    function bin_molno_from_integer_division(n) result(bin)
        use computational_constants_MD, only : halfdomain, nhb
        use calculated_properties_MD, only : binsize
        use arrays_MD, only : r
        implicit none

        integer,intent(in)       :: n
	    integer,dimension(3) 	 :: bin

        bin = ceiling((r(:,n)+halfdomain)/binsize)+nhb

    end function bin_molno_from_integer_division

    function zeta(y, z)
        use computational_constants_MD, only : globaldomain
        use messenger, only : globalise
        implicit none

        real(kind(0.d0)),intent(in) :: y, z
	    real(kind(0.d0))            :: zeta

	    real(kind(0.d0))            :: yg, zg

        yg = globalise(y,2)
        zg = globalise(z,3)
        zeta = 5.d0*sin(2.d0*3.141592*yg/globaldomain(2))

    end function zeta

    subroutine debug_sine_to_bilinear_surface_coeff()!Abilinear)
        use calculated_properties_MD, only : binsize, nbins
        use computational_constants_MD, only : globaldomain
        use messenger, only : globalise_bin
        implicit none

        !real(kind(0.d0)),intent(in),dimension(:,:),allocatable ::Abilinear

        logical, save    :: first_time = .true.
        integer          :: j, k, gbin(3)
        real(kind(0.d0)) :: ysb, yst, zsb, zst, y(2,2), z(2,2), P(2,2), A(2,2)

        if (first_time) then
            allocate(Abilinear(2,2,nbins(2)+2, nbins(3)+2))
            first_time = .false.
        endif

	    do j = 1,nbins(2)+2
	    do k = 1,nbins(3)+2
            gbin = globalise_bin((/ 1, j, k /))
            ysb = (gbin(2)-3)*binsize(2)
            yst = (gbin(2)-2)*binsize(2)
            zsb = (gbin(3)-3)*binsize(3)
            zst = (gbin(3)-2)*binsize(3)
            P(1,1) = zeta(ysb, zsb); y(1,1) = ysb; z(1,1) = zsb 
            P(2,1) = zeta(yst, zsb); y(2,1) = yst; z(2,1) = zsb 
            P(1,2) = zeta(ysb, zst); y(1,2) = ysb; z(1,2) = zst    
            P(2,2) = zeta(yst, zst); y(2,2) = yst; z(2,2) = zst    
!            P(1,1) = 5.d0*sin(2.d0*3.141592*ysb/globaldomain(2))
!            P(2,1) = 5.d0*sin(2.d0*3.141592*yst/globaldomain(2))
!            P(1,2) = 5.d0*sin(2.d0*3.141592*ysb/globaldomain(2))
!            P(2,2) = 5.d0*sin(2.d0*3.141592*yst/globaldomain(2))
            !print*, j, k, ysb, yst, P(1,1), P(2,1)
            call save_bin_bilinear_surface_coeff_unitsquare(y, z, P, A)
            Abilinear(:,:,j,k) = A(:,:)
        enddo
        enddo

    end subroutine debug_sine_to_bilinear_surface_coeff

    !Save coefficiients to matrix for bilinear 
    !This works assuming x and y between zero and one
    subroutine save_bin_bilinear_surface_coeff_unitsquare(x, y, P, A)

        real(kind(0.d0)),intent(in)  :: x(2,2), y(2,2)
        real(kind(0.d0)),intent(in)  :: P(2,2)
	    real(kind(0.d0)),intent(out) :: A(2,2)

        A(1,1) = P(1,1)
        A(2,1) = P(2,1) - P(1,1)
        A(1,2) = P(1,2) - P(1,1)
        A(2,2) = P(2,2) + P(1,1) - (P(2,1) + P(1,2))

    end subroutine save_bin_bilinear_surface_coeff_unitsquare


    !A version getting the explicit location from the intrinsic surface
    function bin_from_full_intrinsic(r) result(bin)
        implicit none

        real(kind(0.d0)),intent(in),dimension(3) :: r
	    integer,dimension(3) 					 :: bin

		bin = ISR_mdt%get_bin(r)

    end function bin_from_full_intrinsic


    !A version getting the explicit location from the intrinsic surface
    function bin_molno_from_full_intrinsic(n) result(bin)
        use computational_constants_MD, only : halfdomain, nhb
        use calculated_properties_MD, only : binsize, nbins
        use arrays_MD, only : r, intnscshift
        implicit none

        integer,intent(in)          :: n
	    integer,dimension(3) 		:: bin

        integer :: ixyz,i,j,k

        i=ISR%normal
        j=ISR%ixyz
        k=ISR%jxyz

        bin(i) = ceiling((r(i,n)+halfdomain(i)-intnscshift(n)+0.5d0*binsize(i))/binsize(i))+nhb(i)
        bin(j) = ceiling((r(j,n)+halfdomain(j))/binsize(j))+nhb(j)
        bin(k) = ceiling((r(k,n)+halfdomain(k))/binsize(k))+nhb(k)

        !Prevents out of range values
        ixyz = 1
        if (bin(ixyz) > nbins(ixyz)+nhb(ixyz)) then
            bin(ixyz) = nbins(ixyz)+nhb(ixyz)
        elseif (bin(ixyz) < 1 ) then
            bin(ixyz) = 1   
        endif

    end function bin_molno_from_full_intrinsic

    !Get surface position from positions for a known bin
    ! e.g. use with a bilinear surface as follows
    ! i = ceiling((r(2)+halfdomain(2))/binsize(2))+nhb(2)
    ! j = ceiling((r(3)+halfdomain(3))/binsize(3))+nhb(3)
    ! P(:,:) = Abilinear(:,:,i,j)
    ! x = xsurfpos(y, z, P)
    function xsurfpos(y, z, Abilinear)
        implicit none

        real(kind(0.d0)),intent(in) :: y, z, Abilinear(2,2)

        integer                     :: jxyz, kxyz
	    real(kind(0.d0))            :: xsurfpos

        xsurfpos = 0.d0
        do jxyz=1,2
        do kxyz=1,2
            xsurfpos = xsurfpos + Abilinear(jxyz, kxyz)*(y**(jxyz-1))*(z**(kxyz-1))
        enddo
        enddo

    end function xsurfpos

    function bin_from_bilinear(r) result(bin)
        use computational_constants_MD, only : halfdomain, nhb
        use calculated_properties_MD, only : binsize, nbins
        implicit none

        real(kind(0.d0)),intent(in),dimension(3) :: r
	    integer,dimension(3) 					 :: bin
        real(kind(0.d0))                         :: elevation

        integer :: ixyz, binshift

        bin = bin_from_integer_division(r)

        !Use convention that no halos in Abilinear values (as no idea how to parallise all this)
        elevation = xsurfpos(r(2), r(3), Abilinear(:,:,bin(2)-nhb(2), bin(3)-nhb(3)))
        binshift = ceiling(elevation/binsize(1))
        bin(1) = bin(1) - binshift


        !print*, "Bin before = ", ceiling((r(1)+halfdomain(1))/binsize(1))+nhb(1), "bin after = ", bin(1)

        !Prevents out of range values
        do ixyz=1,3
    		if (bin(ixyz) > nbins(ixyz)+nhb(ixyz)) then
                bin(ixyz) = nbins(ixyz)+nhb(ixyz)
            elseif (bin(ixyz) < 1 ) then
                bin(ixyz) = 1   
            endif
        enddo

    end function bin_from_bilinear



    subroutine get_crossings(r1, r2, bin1, bin2, n, rc, crossings)
        use computational_constants_MD, only : halfdomain, nhb
        use calculated_properties_MD, only : binsize
        implicit none

        integer, intent(in)                      :: n   !normal
        integer, intent(in), dimension(3) 	     :: bin1, bin2
        real(kind(0.d0)),intent(in),dimension(3) :: r1, r2
        
        logical, intent(out)                    :: crossings
        real(kind(0.d0)),intent(out),dimension(:,:), allocatable :: rc

        integer                 :: t1, t2, i, j, maxbin, minbin
        double precision        :: pt, r12(3)

        if (bin1(n) .eq. bin2(n)) then
            crossings = .false.
        else
            crossings = .true.

            !Tangents
            t1 = mod(n,3)+1
            t2 = mod(n+1,3)+1
            minbin = min(bin1(n), bin2(n))
            maxbin = max(bin1(n), bin2(n))

	        r12   = r1 - r2							!Molecule i trajectory between t-dt and t
	        where (abs(r12) .lt. 0.000001d0) r12 = 0.000001d0

            !print*, "rc size", maxbin-minbin, bin1, bin2
            allocate(rc(3, maxbin-minbin))
            j = 1
            do i = minbin, maxbin-1
		        pt = (i-1*nhb(n))*binsize(n)-halfdomain(n)
                rc(n,j) = pt
                rc(t1,j) = r1(t1)+(r12(t1)/r12(n))*(pt-r1(n))            
                rc(t2,j) = r1(t2)+(r12(t2)/r12(n))*(pt-r1(n))
                !if (i .eq. 120) then
                !    print*, "crossings=", iter, i, rc(:,j)
                !endif
                j = j + 1
            enddo

        endif

    end subroutine get_crossings


	!Identical to get crossings but with a shift of 0.5
    subroutine get_crossings_test(r1, r2, bin1, bin2, n, rc, crossings)
        use computational_constants_MD, only : halfdomain, nhb
        use calculated_properties_MD, only : binsize
        implicit none

        integer, intent(in)                      :: n   !normal
        integer, intent(in), dimension(3) 	     :: bin1, bin2
        real(kind(0.d0)),intent(in),dimension(3) :: r1, r2
        
        logical, intent(out)                    :: crossings
        real(kind(0.d0)),intent(out),dimension(:,:), allocatable :: rc

        integer                 :: t1, t2, i, j, maxbin, minbin
        double precision        :: pt, r12(3)

        if (bin1(n) .eq. bin2(n)) then
            crossings = .false.
        else
            crossings = .true.

            !Tangents
            t1 = mod(n,3)+1
            t2 = mod(n+1,3)+1
            minbin = min(bin1(n), bin2(n))
            maxbin = max(bin1(n), bin2(n))

	        r12   = r1 - r2							!Molecule i trajectory between t-dt and t
	        where (abs(r12) .lt. 0.000001d0) r12 = 0.000001d0

            allocate(rc(3, maxbin-minbin))
            j = 1
            do i = minbin, maxbin-1
		        pt = (i-1*nhb(n)-0.5d0)*binsize(n)-halfdomain(n)
                rc(n,j) = pt
                rc(t1,j) = r1(t1)+(r12(t1)/r12(n))*(pt-r1(n))            
                rc(t2,j) = r1(t2)+(r12(t2)/r12(n))*(pt-r1(n))
                j = j + 1
            enddo

        endif

    end subroutine get_crossings_test


    ! subroutine get_crossings_bilinear(r1, r2, bin1, bin2, n, rc, crossings, cbins)
        ! use computational_constants_MD, only : halfdomain, nhb
        ! use calculated_properties_MD, only : binsize
        ! use bilnear_intersect, only : line_plane_intersect
        ! implicit none

        ! integer, intent(in)                      :: n   !normal
        ! integer, intent(in), dimension(3) 	     :: bin1, bin2
        ! real(kind(0.d0)),intent(in),dimension(3) :: r1, r2
        
        ! logical, intent(out)                    :: crossings
        ! integer,intent(out), optional, dimension(:), allocatable :: cbins
        ! real(kind(0.d0)),intent(out),dimension(:,:), allocatable :: rc

        ! integer                 :: t1, t2, i, j, k, jb, kb, ixyz
        ! integer                 :: maxbin, minbin, minTbin, maxTbin
        ! real(kind(0.d0))        :: pt, r12(3)

        ! integer :: flag, ss, Ns, cross, bc, tempsize
        ! integer, dimension(3) :: bin, bin_mdt
        ! integer, dimension(:), allocatable :: cbinstemp

        ! real(kind(0.d0)) :: yrb, yrt, zrb, zrt, s, ds
        ! real(kind(0.d0)) :: y(2,2), z(2,2), P(2,2,3), A(2,2)
        ! real(kind(0.d0)), dimension(:), allocatable :: elevation, yt, yb, zt, zb
        ! real(kind(0.d0)),dimension(3,2)  :: intersect, normal
        ! real(kind(0.d0)),dimension(:,:), allocatable  :: points, temp

        ! if (bin1(n) .eq. bin2(n)) then
            ! crossings = .false.
        ! else
            ! crossings = .true.

            ! !if (all(ISR_mdt%get_bin(r1, nbins, nhb) .ne. bin1)) stop "get_crossings_bilinear - Bin check failure"
            ! !if (all(ISR_mdt%get_bin(r2, nbins, nhb) .ne. bin2)) stop "get_crossings_bilinear - Bin check failure"

            ! !Tangents
            ! if (ISR_mdt%normal .ne. n) stop "Error in get_crossings_bilinear, normal not consistent with surface"
            ! if (ISR_mdt%ixyz .ne. mod(n,3)+1) stop "Error in get_crossings_bilinear, normal not consistent with surface"
            ! if (ISR_mdt%jxyz .ne. mod(n+1,3)+1) stop "Error in get_crossings_bilinear, normal not consistent with surface"
            ! t1 = ISR_mdt%ixyz !mod(n,3)+1
            ! t2 = ISR_mdt%jxyz !mod(n+1,3)+1

            ! minbin = min(bin1(n), bin2(n))
            ! maxbin = max(bin1(n), bin2(n))

	        ! r12   = r2 - r1  !Molecule i trajectory between t-dt and t
	        ! where (abs(r12) .lt. 0.000001d0) r12 = 0.000001d0

            ! !Get positions of molecules which
            ! !describe a bounding box
            ! yrb = min(r1(t1), r2(t1))
            ! yrt = max(r1(t1), r2(t1))
            ! zrb = min(r1(t2), r2(t2))
            ! zrt = max(r1(t2), r2(t2))

            ! allocate(points(4,3))
            ! !Factor of two as can be maximum of two per surface
			! tempsize = 1
			! do i=1,3
				! minTbin = min(bin1(i), bin2(i))
				! maxTbin = max(bin1(i), bin2(i))
				! tempsize = tempsize * (maxTbin-minTbin+1)
			! enddo
            ! allocate(temp(3, 2*tempsize))
            ! if (present(cbins)) allocate(cbinstemp(2*tempsize))
			
            ! !loop over bin range
            ! bc = 1
            ! do i = minbin, maxbin-1

                ! !Use either distance between molecules or binside
                ! !to create a bilinear patch
                ! if (binsize(t1) .gt. yrt-yrb) then 
                    ! allocate(yb(1), yt(1))
                    ! yb(1) = yrb
                    ! yt(1) = yrt
                ! else
                    ! !Check all intermediate bins in y and z
                    ! minTbin = min(bin1(t1), bin2(t1))
                    ! maxTbin = max(bin1(t1), bin2(t1))
                    ! allocate(yb(maxTbin-minTbin+1))
                    ! allocate(yt(maxTbin-minTbin+1))
                    ! j = 1
                    ! do jb = minTbin, maxTbin
                        ! yb(j) = (jb-nhb(t1)-1)*binsize(t1)-0.5d0*ISR_mdt%box(t1)
                        ! yt(j) = (jb-nhb(t1)  )*binsize(t1)-0.5d0*ISR_mdt%box(t1)
                        ! j = j + 1
                    ! enddo
                ! endif

                ! !Use either distance between molecules or binside
                ! !to create a bilinear patch
                ! if (binsize(t2) .gt. zrt-zrb) then 
                    ! allocate(zb(1), zt(1))
                    ! zb(1) = zrb
                    ! zt(1) = zrt
                ! else
                    ! minTbin = min(bin1(t2), bin2(t2))
                    ! maxTbin = max(bin1(t2), bin2(t2))
                    ! allocate(zb(maxTbin-minTbin+1))
                    ! allocate(zt(maxTbin-minTbin+1))
                    ! j = 1
                    ! do kb = minTbin, maxTbin
                        ! zb(j) = (kb-nhb(t2)-1)*binsize(t2)-0.5d0*ISR_mdt%box(t2)
                        ! zt(j) = (kb-nhb(t2)  )*binsize(t2)-0.5d0*ISR_mdt%box(t2)
                        ! j = j + 1
                    ! enddo
                ! endif

                ! j = 1
                ! do jb = 1, size(yb,1)
                ! do kb = 1, size(zb,1)
				
                    ! points(:,n) = 0.5d0*(r1(ISR_mdt%normal)+r2(ISR_mdt%normal)) !Not used
                    ! points(1,t1) = yb(jb); points(1,t2) = zb(kb)
                    ! points(2,t1) = yb(jb); points(2,t2) = zt(kb) 
                    ! points(3,t1) = yt(jb); points(3,t2) = zb(kb) 
                    ! points(4,t1) = yt(jb); points(4,t2) = zt(kb) 

                    ! !Get surface in x as a function of y and z
                    ! call ISR_mdt%get_surface(points, elevation, include_zeromode=.true.)

                    ! !Shift by current bin
! !                    elevation(:) = elevation(:) & 
! !                                  + (i-1*nhb(n)+0.5d0)*binsize(n) &
! !                                  -0.5d0*ISR_mdt%box(ISR_mdt%normal)
                    ! elevation(:) = elevation(:) & 
                                  ! + (i-1*nhb(n)-0.5d0)*binsize(n) &
                                  ! -0.5d0*ISR_mdt%box(ISR_mdt%normal)  !HALF SHIFT

                    ! points(:,n) = elevation(:)  !Normal not necessarily used but for completeness

                    ! !Get a patch of P values to use in bilinear patch - line calculation
                    ! P(:,:,1) = reshape(points(:,n ), (/2,2/))
                    ! P(:,:,2) = reshape(points(:,t1), (/2,2/))
                    ! P(:,:,3) = reshape(points(:,t2), (/2,2/))

                    ! !print'(a,2i5,18f10.5)', "x_bilinear", i,j, p(1,1,:), P(1,2,:), P(2,1,:), p(2,2,:), r1(:), r12(:)

                    ! !Special case of flat surface causes problems so need to handle separatly
                    ! if (P(1,1,1) .eq. P(2,1,1) .and. &
                        ! P(2,1,1) .eq. P(1,2,1) .and. &
                        ! P(1,2,1) .eq. P(2,2,1)) then
                        ! intersect = -666
	                    ! intersect(1,1) = P(1,1,1)
                        ! intersect(2,1) = r1(t1)+(r12(t1)/r12(n))*(intersect(1,1)-r1(n))            
                        ! intersect(3,1) = r1(t2)+(r12(t2)/r12(n))*(intersect(1,1)-r1(n))
                        ! flag = 1
                    ! else
                        ! call line_plane_intersect(r1, r12, P, intersect, normal, flag)
                    ! endif

                    ! print*, 'line_plane_intersect output', i, r1(1), elevation(1), r2(1), bin1(1)-nhb(1), yb, yt, zb, zt, intersect
                    ! !Loop over intersects and add to temp
                    ! do ixyz=1,size(intersect,2)
                        ! if (all(abs(intersect(:,ixyz)+666) .gt. 1e-7)) then
                            ! if (size(intersect,2) .ne. 1) print'(a,6i5,9f10.5)', "intrsect", ixyz, i, jb, kb, j,& 
                                                                                 ! bc, r1, intersect(:,ixyz), r2
                            ! temp(:,bc) = intersect(:,ixyz)
                            ! if (present(cbins)) cbinstemp(bc) = i
                            ! j = j + 1
                            ! bc = bc + 1
                        ! endif
                    ! enddo

! !                    if (i .eq. 120) then
 ! !                       print*, "crossings=", iter, i, temp(:,j-1) !yb(jb), zb(kb), yt(jb), zt(kb),
                        ! !stop
 ! !                   endif

                ! enddo
                ! enddo

                ! deallocate(yb, yt, zb, zt)

            ! enddo

            ! print'(a,4i6,6f10.5)', "Crossings ", maxbin-minbin, size(yb,1), size(zb,1), bc, r1, r2

            ! !Copy array of intersects to rc
            ! if (bc .gt. 1) then
                ! allocate(rc(3, bc-1))
                ! rc = temp(:,1:bc-1)
                ! if (present(cbins)) then
                    ! allocate(cbins(bc-1))
                    ! cbins = cbinstemp(1:bc-1)
                ! endif
            ! else
                ! print*, "Warning in get_crossings_bilinear - no crossings found ", minbin, maxbin, r1, r2
                ! ! Crossing of bilinear differs from Fourier surface
                ! ! Walk line between two points and try to get location of crossing
                ! Ns = 100
                ! if (allocated(points)) deallocate(points) 
                ! allocate(points(Ns+2,3))
            	! ds = 1.d0 / real(Ns, kind(0.d0))
            	! ! First sample at r1 
	            ! s = -ds
	            ! ! Loop over all samples, s varies from 0 to 1
                ! bin_mdt = bin1
	            ! do ss = 1, Ns+2
		            ! ! Position of sample on line
		            ! points(ss,:) = r1(:) + s*r12(:)
                    ! bin = ISR_mdt%get_bin(points(ss,:), nbins, nhb)
                    ! !print*, "Points on line", s, points(ss,:), bin, bin_mdt, &
                    ! !bin(ISR_mdt%normal) .ne. bin_mdt(ISR_mdt%normal), ISR_mdt%normal
                    ! !If bin changes then must be a crossing
                    ! cross = bin(ISR_mdt%normal) - bin_mdt(ISR_mdt%normal) 
                    ! if (cross .ne. 0) then
                        ! allocate(rc(3, 1))
                        ! rc(:, 1) = r1(:) + (s-0.5d0*ds)*r12(:)
                        ! !print*, "Crossing found", rc, bin, s-0.5d0*ds
                        ! if (present(cbins)) then
                            ! allocate(cbins(1))
                            ! if (cross .gt. 0) then
                                ! cbins(1) = bin_mdt(ISR_mdt%normal)
                            ! else 
                                ! cbins(1) = bin(ISR_mdt%normal)
                            ! endif
                        ! endif
                        ! exit
                    ! endif
                    ! bin_mdt = bin
		            ! s = s + ds
	            ! end do	
                ! !call ISR_mdt%get_surface(points, elevation)
                ! !print*, "Elevations", elevation+ (bin1(n)-1*nhb(n))*binsize(n)-0.5d0*ISR_mdt%box(ISR_mdt%normal)
                ! !stop "Error get_crossings_bilinear - interactions must be missed"
            ! endif
        ! endif

    ! end subroutine get_crossings_bilinear


    ! subroutine get_crossings_bilinear_opt(r1, r2, bin1, bin2, n, rc, crossings, cbins)
        ! use computational_constants_MD, only : halfdomain, nhb
        ! use calculated_properties_MD, only : binsize
        ! use bilnear_intersect, only : line_plane_intersect
        ! implicit none

        ! integer, intent(in)                      :: n   !normal
        ! integer, intent(in), dimension(3) 	     :: bin1, bin2
        ! real(kind(0.d0)),intent(in),dimension(3) :: r1, r2
        
        ! logical, intent(out)                    :: crossings
        ! integer,intent(out), optional, dimension(:), allocatable :: cbins
        ! real(kind(0.d0)),intent(out),dimension(:,:), allocatable :: rc

        ! integer                 :: t1, t2, i, j, k, jb, kb, ixyz
        ! integer                 :: maxbin, minbin, minTbin, maxTbin
        ! real(kind(0.d0))        :: pt, r12(3)

        ! integer :: flag, ss, Ns, cross, bc, tempsize
        ! integer, dimension(3) :: bin, bin_mdt
        ! integer, dimension(:), allocatable :: cbinstemp

        ! real(kind(0.d0)) :: yrb, yrt, zrb, zrt, s, ds
        ! real(kind(0.d0)),save :: maxerror
        ! real(kind(0.d0)) :: y(2,2), z(2,2), P(2,2,3), A(2,2), rb(3)
        ! real(kind(0.d0)), dimension(:), allocatable :: elevation
        ! real(kind(0.d0)),dimension(3,2)  :: intersect, normal, intersect2
        ! real(kind(0.d0)),dimension(:,:), allocatable  :: points, points2, temp

        ! if (bin1(n) .eq. bin2(n)) then
            ! crossings = .false.
        ! else
            ! crossings = .true.

            ! !Tangents
            ! if (ISR_mdt%normal .ne. n) stop "Error in get_crossings_bilinear, normal not consistent with surface"
            ! if (ISR_mdt%ixyz .ne. mod(n,3)+1) stop "Error in get_crossings_bilinear, normal not consistent with surface"
            ! if (ISR_mdt%jxyz .ne. mod(n+1,3)+1) stop "Error in get_crossings_bilinear, normal not consistent with surface"
            ! t1 = ISR_mdt%ixyz !mod(n,3)+1
            ! t2 = ISR_mdt%jxyz !mod(n+1,3)+1

            ! minbin = min(bin1(n), bin2(n))
            ! maxbin = max(bin1(n), bin2(n))

	        ! r12   = r2 - r1  !Molecule i trajectory between t-dt and t
	        ! where (abs(r12) .lt. 0.000001d0) r12 = 0.000001d0

            ! !Get positions of molecules which
            ! !describe a bounding box
            ! yrb = min(r1(t1), r2(t1))
            ! yrt = max(r1(t1), r2(t1))
            ! zrb = min(r1(t2), r2(t2))
            ! zrt = max(r1(t2), r2(t2))

            ! allocate(points(4,3))
            ! allocate(points2(4,3))
            ! !Factor of two as can be maximum of two per surface
			! tempsize = 1
			! do i=1,3
				! minTbin = min(bin1(i), bin2(i))
				! maxTbin = max(bin1(i), bin2(i))
				! tempsize = tempsize * (maxTbin-minTbin+1)
			! enddo
            ! allocate(temp(3, 2*tempsize))
            ! if (present(cbins)) allocate(cbinstemp(2*tempsize))

        	! ! First sample at r1 
            ! Ns = maxbin-minbin
        	! ds = 1.d0 / real(Ns, kind(0.d0))
            ! s = ds
			
            ! !loop over bin range
            ! bc = 1
            ! do i = minbin, maxbin-1

	            ! ! Loop over all samples, s varies from 0 to 1
	            ! ! Position of sample on line
	            ! rb(:) = r1(:) + s*r12(:)
                ! yrb = min(r1(t1), rb(t1))
                ! yrt = max(r1(t1), rb(t1))
                ! zrb = min(r1(t2), rb(t2))
                ! zrt = max(r1(t2), rb(t2))
                ! !print'(a,3i5,10f10.5)', "patch", i, bc, Ns, s, r1(:), rb(:), r2(:)
                ! !print'(a,3i5,8f10.5)', "patch", i, bc, Ns, s, rb(:), yrb, yrt, zrb, zrt
	            ! s = s + ds

                ! !Use distance between molecules
                ! !to create a bilinear patch
                ! points(:,n) = rb(n) !Not used
                ! points(1,t1) = yrb; points(1,t2) = zrb
                ! points(2,t1) = yrb; points(2,t2) = zrt 
                ! points(3,t1) = yrt; points(3,t2) = zrb 
                ! points(4,t1) = yrt; points(4,t2) = zrt

                ! !Get surface in x as a function of y and z
                ! call ISR_mdt%get_surface(points, elevation, include_zeromode=.true.)

                ! !Normal direction used in line plane intersect
                ! !Shift by current bin
                ! points(:,n) = elevation(:) & 
                              ! + (i-1*nhb(n)-0.5d0)*binsize(n) &
                              ! -0.5d0*ISR_mdt%box(n) 

                ! !Get a patch of P values to use in ` patch & line calculation
                ! P(:,:,1) = reshape(points(:,n ), (/2,2/))
                ! P(:,:,2) = reshape(points(:,t1), (/2,2/))
                ! P(:,:,3) = reshape(points(:,t2), (/2,2/))

                ! !print'(a,2i5,18f10.5)', "x_bilinear", i,j, p(1,1,:), P(1,2,:), P(2,1,:), p(2,2,:), r1(:), r12(:)

                ! !Special case of flat surface causes problems so need to handle separatly
                ! if (P(1,1,1) .eq. P(2,1,1) .and. &
                    ! P(2,1,1) .eq. P(1,2,1) .and. &
                    ! P(1,2,1) .eq. P(2,2,1)) then
                    ! intersect = -666
                    ! intersect(1,1) = P(1,1,1)
                    ! intersect(2,1) = r1(t1)+(r12(t1)/r12(n))*(intersect(1,1)-r1(n))            
                    ! intersect(3,1) = r1(t2)+(r12(t2)/r12(n))*(intersect(1,1)-r1(n))
                    ! flag = 1
                ! else
                    ! call line_plane_intersect(r1, r12, P, intersect, normal, flag)
                ! endif
				
				
				! !Get surface sampled from bins
				! ! call ISR_mdt%get_sampled_surface(rb, points2)
				! ! points2(:,n) = points2(:,n) & 
							  ! ! + (i-1*nhb(n)-0.5d0)*binsize(n) &
                              ! ! -0.5d0*ISR_mdt%box(n)
				! ! P(:,:,1) = reshape(points2(:,n ), (/2,2/))
                ! ! P(:,:,2) = reshape(points2(:,t1), (/2,2/))
                ! ! P(:,:,3) = reshape(points2(:,t2), (/2,2/))
				! ! call line_plane_intersect(r1, r12, P, intersect2, normal, flag)
				! ! if (abs(intersect2(1,1)+666) .lt. 1e-7) then
					! ! intersect2(:,1) = r1(:) + (s-0.5d0*ds)*r12(:)
				! ! endif
				! ! !if (all(abs(intersect2(:,1)+666) .gt. 1e-7)) then
					! ! print'(a,13f10.5)', "Elevation check ", & 
						! ! r1, rb, points2(:,n), intersect(:,1)-intersect2(:,1)
				! !endif

				
                ! !print*, 'line_plane_intersect output', i, r1(1), elevation(1), r2(1), bin1(1)-nhb(1), yrb, yrt, zrb, zrt, intersect
                ! !Loop over intersects and add to temp
                ! do ixyz=1,size(intersect,2)
                    ! if (all(abs(intersect(:,ixyz)+666) .gt. 1e-7)) then
                        ! !if (size(intersect,2) .ne. 1) print'(a,3i5,9f10.5)', "intrsect", ixyz, i, bc, r1, intersect(:,ixyz), r2
                        ! temp(:,bc) = intersect(:,ixyz)
                        ! if (present(cbins)) cbinstemp(bc) = i
                        ! bc = bc + 1
                    ! endif
                ! enddo

            ! enddo

            ! !print'(a,4i6,6f10.5)', "Crossings ", maxbin-minbin, size(yb,1), size(zb,1), bc, r1, r2

            ! !Copy array of intersects to rc
            ! if (bc .gt. 1) then
                ! allocate(rc(3, bc-1))
                ! rc = temp(:,1:bc-1)
                ! if (present(cbins)) then
                    ! allocate(cbins(bc-1))
                    ! cbins = cbinstemp(1:bc-1)
                ! endif
            ! else
                ! print*, "Warning in get_crossings_bilinear - no crossings found ", minbin, maxbin, r1, r2
                ! ! Crossing of bilinear differs from Fourier surface
                ! ! Walk line between two points and try to get location of crossing
                ! Ns = maxbin-minbin+1
                ! if (allocated(points)) deallocate(points)
                ! allocate(points(Ns+2,3))
            	! ds = 1.d0 / real(Ns, kind(0.d0))
            	! ! First sample at r1 
	            ! s = -ds
	            ! ! Loop over all samples, s varies from 0 to 1
                ! bin_mdt = bin1
	            ! do ss = 1, Ns+2
		            ! ! Position of sample on line
		            ! points(ss,:) = r1(:) + s*r12(:)
                    ! bin = ISR_mdt%get_bin(points(ss,:), nbins, nhb)
                    ! !print*, "Points on line", s, points(ss,:), bin, bin_mdt, &
                    ! !bin(ISR_mdt%normal) .ne. bin_mdt(ISR_mdt%normal), ISR_mdt%normal
                    ! !If bin changes then must be a crossing
                    ! cross = bin(ISR_mdt%normal) - bin_mdt(ISR_mdt%normal) 
                    ! if (cross .eq. 1) then
                        ! allocate(rc(3, 1))
                        ! rc(:, 1) = r1(:) + (s-0.5d0*ds)*r12(:)
                        ! !print*, "Crossing found", rc, bin, s-0.5d0*ds
                        ! if (present(cbins)) then
                            ! allocate(cbins(1))
                            ! if (cross .gt. 0) then
                                ! cbins(1) = bin_mdt(ISR_mdt%normal)
                            ! else 
                                ! cbins(1) = bin(ISR_mdt%normal)
                            ! endif
                        ! endif
                        ! exit
                    ! endif
                    ! bin_mdt = bin
		            ! s = s + ds
	            ! end do	
                ! !call ISR_mdt%get_surface(points, elevation)
                ! !print*, "Elevations", elevation+ (bin1(n)-1*nhb(n))*binsize(n)-0.5d0*ISR_mdt%box(ISR_mdt%normal)
                ! !stop "Error get_crossings_bilinear - interactions must be missed"
            ! endif
        ! endif

    ! end subroutine get_crossings_bilinear_opt



    ! subroutine get_crossings_bilinear_lookup(r1, r2, bin1, bin2, n, rc, crossings, cbins)
      ! use computational_constants_MD, only : halfdomain, nhb
        ! use calculated_properties_MD, only : binsize
        ! use bilnear_intersect, only : line_plane_intersect
        ! implicit none

        ! integer, intent(in)                      :: n   !normal
        ! integer, intent(in), dimension(3) 	     :: bin1, bin2
        ! real(kind(0.d0)),intent(in),dimension(3) :: r1, r2
        
        ! logical, intent(out)                    :: crossings
        ! integer,intent(out), optional, dimension(:), allocatable :: cbins
        ! real(kind(0.d0)),intent(out),dimension(:,:), allocatable :: rc

        ! integer                 :: t1, t2, i, j, k, jb, kb, ixyz
        ! integer                 :: maxbin, minbin, minTbin, maxTbin
        ! real(kind(0.d0))        :: pt, r12(3)

        ! integer :: flag, ss, Ns, cross, bc, tempsize
        ! integer, dimension(3) :: bin, bin_mdt
        ! integer, dimension(:), allocatable :: cbinstemp

        ! real(kind(0.d0)) :: yrb, yrt, zrb, zrt, s, ds
        ! real(kind(0.d0)),save :: maxerror
        ! real(kind(0.d0)) :: y(2,2), z(2,2), P(2,2,3), A(2,2), rb(3)
        ! real(kind(0.d0)), dimension(:), allocatable :: elevation
        ! real(kind(0.d0)),dimension(3,2)  :: intersect, normal
        ! real(kind(0.d0)),dimension(:,:), allocatable  :: points, temp

        ! if (bin1(n) .eq. bin2(n)) then
            ! crossings = .false.
        ! else
            ! crossings = .true.

            ! !Tangents
            ! if (ISR_mdt%normal .ne. n) stop "Error in get_crossings_bilinear, normal not consistent with surface"
            ! if (ISR_mdt%ixyz .ne. mod(n,3)+1) stop "Error in get_crossings_bilinear, normal not consistent with surface"
            ! if (ISR_mdt%jxyz .ne. mod(n+1,3)+1) stop "Error in get_crossings_bilinear, normal not consistent with surface"
            ! t1 = ISR_mdt%ixyz !mod(n,3)+1
            ! t2 = ISR_mdt%jxyz !mod(n+1,3)+1

            ! minbin = min(bin1(n), bin2(n))
            ! maxbin = max(bin1(n), bin2(n))

	        ! r12   = r2 - r1  !Molecule i trajectory between t-dt and t
	        ! where (abs(r12) .lt. 0.000001d0) r12 = 0.000001d0

            ! !Get positions of molecules which
            ! !describe a bounding box
            ! yrb = min(r1(t1), r2(t1))
            ! yrt = max(r1(t1), r2(t1))
            ! zrb = min(r1(t2), r2(t2))
            ! zrt = max(r1(t2), r2(t2))

            ! allocate(points(4,3))
            ! !Factor of two as can be maximum of two per surface
			! tempsize = 1
			! do i=1,3
				! minTbin = min(bin1(i), bin2(i))
				! maxTbin = max(bin1(i), bin2(i))
				! tempsize = tempsize * (maxTbin-minTbin+1)
			! enddo
            ! allocate(temp(3, 2*tempsize))
            ! if (present(cbins)) allocate(cbinstemp(2*tempsize))

        	! ! First sample at r1 
            ! Ns = maxbin-minbin
        	! ds = 1.d0 / real(Ns, kind(0.d0))
            ! s = ds
			
            ! !loop over bin range
            ! bc = 1
            ! do i = minbin, maxbin-1

	            ! ! Loop over all samples, s varies from 0 to 1
	            ! ! Position of sample on line
	            ! rb(:) = r1(:) + s*r12(:)
	            ! s = s + ds

				! !Get surface sampled from bins
				! call ISR_mdt%get_sampled_surface(rb, points)
				! points(:,n) = points(:,n) & 
							  ! + (i-1*nhb(n)-0.5d0)*binsize(n) &
                              ! -0.5d0*ISR_mdt%box(n)
				! P(:,:,1) = reshape(points(:,n ), (/2,2/))
                ! P(:,:,2) = reshape(points(:,t1), (/2,2/))
                ! P(:,:,3) = reshape(points(:,t2), (/2,2/))
				! call line_plane_intersect(r1, r12, P, intersect, normal, flag)
				! !If intersect not found, just assume it's half way
				! if (abs(intersect(1,1)+666) .lt. 1e-7) then
					! intersect(:,1) = r1(:) + (s-0.5d0*ds)*r12(:)
				! endif
						
				! !print'(a,2i5,18f10.5)', "x_bilinear", i,j, p(1,1,:), P(1,2,:), P(2,1,:), p(2,2,:), r1(:), r12(:)

                ! !Special case of flat surface causes problems so need to handle separatly
                ! if (P(1,1,1) .eq. P(2,1,1) .and. &
                    ! P(2,1,1) .eq. P(1,2,1) .and. &
                    ! P(1,2,1) .eq. P(2,2,1)) then
                    ! intersect = -666
                    ! intersect(1,1) = P(1,1,1)
                    ! intersect(2,1) = r1(t1)+(r12(t1)/r12(n))*(intersect(1,1)-r1(n))            
                    ! intersect(3,1) = r1(t2)+(r12(t2)/r12(n))*(intersect(1,1)-r1(n))
                    ! flag = 1
                ! else
                    ! call line_plane_intersect(r1, r12, P, intersect, normal, flag)
                ! endif

                ! !print*, 'line_plane_intersect output', i, r1(1), elevation(1), r2(1), bin1(1)-nhb(1), yrb, yrt, zrb, zrt, intersect
                ! !Loop over intersects and add to temp
                ! do ixyz=1,size(intersect,2)
                    ! if (all(abs(intersect(:,ixyz)+666) .gt. 1e-7)) then
                        ! !if (size(intersect,2) .ne. 1) print'(a,3i5,9f10.5)', "intrsect", ixyz, i, bc, r1, intersect(:,ixyz), r2
                        ! temp(:,bc) = intersect(:,ixyz)
                        ! if (present(cbins)) cbinstemp(bc) = i
                        ! bc = bc + 1
                    ! endif
                ! enddo

            ! enddo

            ! !print'(a,4i6,6f10.5)', "Crossings ", maxbin-minbin, size(yb,1), size(zb,1), bc, r1, r2

            ! !Copy array of intersects to rc
            ! if (bc .gt. 1) then
                ! allocate(rc(3, bc-1))
                ! rc = temp(:,1:bc-1)
                ! if (present(cbins)) then
                    ! allocate(cbins(bc-1))
                    ! cbins = cbinstemp(1:bc-1)
                ! endif
            ! else
                ! print*, "Warning in get_crossings_bilinear - no crossings found ", minbin, maxbin, r1, r2
                ! ! Crossing of bilinear differs from Fourier surface
                ! ! Walk line between two points and try to get location of crossing
                ! Ns = maxbin-minbin+1
                ! if (allocated(points)) deallocate(points)
                ! allocate(points(Ns+2,3))
            	! ds = 1.d0 / real(Ns, kind(0.d0))
            	! ! First sample at r1 
	            ! s = -ds
	            ! ! Loop over all samples, s varies from 0 to 1
                ! bin_mdt = bin1
	            ! do ss = 1, Ns+2
		            ! ! Position of sample on line
		            ! points(ss,:) = r1(:) + s*r12(:)
                    ! bin = ISR_mdt%get_bin(points(ss,:), nbins, nhb)
                    ! !print*, "Points on line", s, points(ss,:), bin, bin_mdt, &
                    ! !bin(ISR_mdt%normal) .ne. bin_mdt(ISR_mdt%normal), ISR_mdt%normal
                    ! !If bin changes then must be a crossing
                    ! cross = bin(ISR_mdt%normal) - bin_mdt(ISR_mdt%normal) 
                    ! if (cross .eq. 1) then
                        ! allocate(rc(3, 1))
                        ! rc(:, 1) = r1(:) + (s-0.5d0*ds)*r12(:)
                        ! !print*, "Crossing found", rc, bin, s-0.5d0*ds
                        ! if (present(cbins)) then
                            ! allocate(cbins(1))
                            ! if (cross .gt. 0) then
                                ! cbins(1) = bin_mdt(ISR_mdt%normal)
                            ! else 
                                ! cbins(1) = bin(ISR_mdt%normal)
                            ! endif
                        ! endif
                        ! exit
                    ! endif
                    ! bin_mdt = bin
		            ! s = s + ds
	            ! end do	
                ! !call ISR_mdt%get_surface(points, elevation)
                ! !print*, "Elevations", elevation+ (bin1(n)-1*nhb(n))*binsize(n)-0.5d0*ISR_mdt%box(ISR_mdt%normal)
                ! !stop "Error get_crossings_bilinear - interactions must be missed"
            ! endif
        ! endif


    ! end subroutine get_crossings_bilinear_lookup



end module module_record
!==========================================================================
! Top level routine which calls all recording statements throughout the 
! code.

subroutine setup_assign_get_bins_fn()
	use module_record
	implicit none

    !Cluster analysis or average bin based tracking of liquid vapour interfaces
    if (cluster_analysis_outflag .eq. 1 .and.  & 
        any(intrinsic_interface_outflag .eq. (/1,2/))) then
        call get_interface_from_clusters()
		get_bin => bin_from_full_intrinsic
		get_bin_molno => bin_molno_from_full_intrinsic
    elseif (cluster_analysis_outflag .eq. 2) then
        call sl_interface_from_binaverage()
    else
        get_bin => bin_from_integer_division
        get_bin_molno => bin_molno_from_integer_division
    endif

end subroutine setup_assign_get_bins_fn

subroutine simulation_record
	use module_record
	use CV_objects, only : CVcheck_mass, CV_debug
	implicit none

    integer, save   :: vmd_skip_count=0
    integer, save   :: vmdintervalno = 1

    !Cluster analysis or average bin based tracking of liquid vapour interfaces
    if (cluster_analysis_outflag .eq. 1) then
        call get_interface_from_clusters()
    elseif (cluster_analysis_outflag .eq. 2) then
        call sl_interface_from_binaverage()
    endif

	if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then
		call mass_flux_averaging(mflux_outflag)				!Average mass flux before movement of particles
		call momentum_flux_averaging(vflux_outflag)         !Average momentum flux after movement of particles
		call energy_flux_averaging(eflux_outflag)			!Average energy flux after movement of particles
        call surface_density_averaging(msurf_outflag)	    !Obtain and record density on a surface
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

#if __INTEL_COMPILER > 1200
	!Obtain and record velocity distributions
	if (vPDF_flag .ne. 0) call velocity_PDF_averaging(vPDF_flag)

	!Record boundary force PDFs
	if (bforce_pdf_measure .ne. 0) call bforce_pdf_stats
#endif

	!Obtain and record temperature
	if (temperature_outflag .ne. 0)	call temperature_averaging(temperature_outflag)

    !Obtain and record energy
	if (energy_outflag .ne. 0)	call energy_averaging(energy_outflag)

    !Obtain and record centre of mass
    if (centre_of_mass_outflag .ne. 0)	call centre_of_mass_averaging(centre_of_mass_outflag)

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
			write(3000+irank,'(i10,f28.4,6f10.4)') n , potenergymol(n), r(:,n), globalise(r(:,n))
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
	write(ranknum, '(i4)') irank	
	write(filename,'(a15,a)') "mol_escape_err.", trim(adjustl(ranknum))

	! Warn user of escape
	print('(a,i8,a,i4,3a)'),' Molecule ',n,' on process ', &
		  irank, ' is outside the domain and halo cells. Inspect ', &
		  filename, 'for more information.'

	! Find unused file unit number
	fileunit = get_new_fileunit()

	! Open file and write escape info
	open(unit=fileunit,file=filename,access='append')

		write(fileunit,'(a,i6,a,i4,a)')' Molecule ',n,' on process ', &
			  irank, ' is outside the domain and halo cells.'
		write(fileunit,'(a,i8,a)')' At iteration ',iter,' it is located at: '
		write(fileunit,'(a,e20.5)')   '    rx: ', r(1,n)
		write(fileunit,'(a,e20.5)')   '    ry: ', r(2,n)
		write(fileunit,'(a,e20.5,a)') '    rz: ', r(3,n), ','
		write(fileunit,'(a,e20.5)')   '    globalrx: ', rglob(1) 
		write(fileunit,'(a,e20.5)')   '    globalry: ', rglob(2) 
		write(fileunit,'(a,e20.5,a)') '    globalrz: ', rglob(3), ','
		write(fileunit,'(a)')         ' with velocity: '
		write(fileunit,'(a,e20.5)')   '    vx: ', v(1,n)
		write(fileunit,'(a,e20.5)')   '    vy: ', v(2,n)
		write(fileunit,'(a,e20.5)')   '    vz: ', v(3,n)
		write(fileunit,'(a)')         ' and acceleration: '
		write(fileunit,'(a,e20.5)')   '    ax: ', a(1,n)
		write(fileunit,'(a,e20.5)')   '    ay: ', a(2,n)
		write(fileunit,'(a,e20.5)')   '    az: ', a(3,n)
		if (ensemble.eq.tag_move) then
			write(fileunit,'(a)')         ' Molecular tag: '
			write(fileunit,'(a,i4)')   '    tag: ', tag(n)
		end if
		write(fileunit,'(a,3i4)')' Processor coords: ',iblock,jblock,kblock

	close(fileunit,status='keep')

end subroutine print_mol_escape_error





!===================================================================================
!Molecules grouped into velocity ranges to give vel frequency distribution graph
!and calculate Boltzmann H function on a bin by bin basis

#if __INTEL_COMPILER > 1200
subroutine velocity_PDF_averaging(ixyz)
	use module_record
	use librarymod, only : Maxwell_Boltzmann_vel,Maxwell_Boltzmann_speed
    use field_io, only : velocity_PDF_slice_io
	use module_set_parameters, only : velPDF, velPDFMB, velPDF_array
	implicit none

    integer,intent(in)      :: ixyz

	integer                 :: n, i,j,k, jxyz,kxyz
	integer, save		    :: sample_count=-1

	integer,dimension(3)	:: cbin
	integer,dimension(:,:),allocatable          :: pdfx,pdfy,pdfz

	real(kind(0.d0)) 	                        :: Hfunction,vfactor,const
	real(kind(0.d0)),save 	                    :: meanstream, meannp
	real(kind(0.d0)),dimension(3)	            :: peculiarv,binsize_
	real(kind(0.d0)),dimension(:),allocatable 	:: vmagnitude,binloc

	sample_count = sample_count + 1
	call cumulative_velocity_PDF

	!Write and reset PDF FOR EACH Y SLAB
	if (sample_count .eq. Nvpdf_ave) then
        sample_count = 0

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
        cbin = get_bin_molno(n)
        do ixyz = 1,nd
    		vmagnitude(1)=v(ixyz,n)
            call velPDF_array(cbin(1),cbin(2),cbin(3),ixyz)%update(vmagnitude)
        enddo
	enddo
	deallocate(vmagnitude)

end subroutine

end subroutine velocity_PDF_averaging
#endif




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
			dV   = pi*((rmag+dr)**2 - rmag**2)
		case(3)
			dV   = (4.d0*pi/3.d0)*((rmag+dr)**3 - rmag**3)
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
subroutine simulation_write_potential_field(xmin,xmax,ymin,ymax,zmin, & 
                                            zmax,xres,yres,zres, &
                                            filenum,casetype)
	use module_compute_forces, only : compute_force_and_potential_at
	implicit none

	integer, intent(in) :: xres,yres,zres, filenum, casetype
	real(kind(0.d0)), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax

	integer :: i,j,k
	real(kind(0.d0)) :: x, y,z, dx, dy, dz
	real(kind(0.d0)) :: U
	real(kind(0.d0)) :: rijave(3), f(3), rf(3,3)

	dx = (xmax - xmin)/real(xres)
	dy = (ymax - ymin)/real(yres)
	dz = (zmax - zmin)/real(zres)

	do i = 0,xres
    do j = 0,yres
	do k = 0,zres

		x  = xmin + i*dx
		y  = ymin + j*dy 
		z  = zmin + k*dz 

        select case(casetype)
        case(0)
			call compute_force_and_potential_at((/x,y,z/),U,f)
		    write(filenum,'(4e12.4)') x, y, z, U
        case(1)
		    call compute_force_and_potential_at((/x,y,z/),U,f,rf=rf)
		    write(filenum,'(12e12.4)') x, y, z, rf
        case(2)
		    call compute_force_and_potential_at((/x,y,z/),U,f,rijave=rijave)
		    write(filenum,'(12f12.4)') x, y, z, rijave
        case(3)
            !Write the lot
		    call compute_force_and_potential_at((/x,y,z/),U,f,rmin=2.d0,rf=rf,rijave=rijave)
		    write(filenum,'(4e12.4)') x, y, z, U
		    write(filenum+10000,'(12e12.4)') x, y, z, rf
		    write(filenum+20000,'(12f12.4)') x, y, z, rijave
        end select

    enddo
    enddo
    enddo

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

	R_g  = sqrt(R_g2)

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
	integer, save					:: sample_count = -1

	sample_count = sample_count + 1
	call cumulative_mass(ixyz)
	if (sample_count .eq. Nmass_ave) then
		sample_count = 0

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

        select case (split_pol_sol_stats)
        case(0)
            do n = 1,np
                !Add up current volume mass densities
                ibin(:) = get_bin_molno(n)
                volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + mass(n)
            enddo
        case(1:2)
            !Seperate logging for polymer and non-polymer 
            do n = 1,np
                !Add up current volume mass densities
                ibin = get_bin_molno(n)
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
	integer, save		:: sample_count = -1
	real(kind(0.d0)),dimension(3) 	:: Vbinsize

	sample_count = sample_count + 1
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
            ib(:) = get_bin_molno(n)
			U(:,n) = volume_momentum(ib(1),ib(2),ib(3),:) / volume_mass(ib(1),ib(2),ib(3))
		enddo

	endif

	if (sample_count .eq. Nvel_ave) then
		sample_count = 0

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

        select case (split_pol_sol_stats)
        case(0)
            do n = 1,np
                !Add up current volume mass and momentum densities
                ibin(:) = get_bin_molno(n)
                volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + mass(n)
                volume_momentum(ibin(1),ibin(2),ibin(3),:) = volume_momentum(ibin(1),ibin(2),ibin(3),:) & 
                                                            + mass(n)*(v(:,n) + slidev(:,n))
            enddo
        case(1)
            do n = 1,np
                !Add up current volume mass and momentum densities
                ibin(:) = get_bin_molno(n)
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
	integer, save		:: sample_count = -1
	
	sample_count = sample_count + 1
	call cumulative_temperature(ixyz)
	if (sample_count .eq. NTemp_ave) then
		sample_count = 0

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

		!Reset Control Volume momentum 
		do n = 1,np

			!Add up current volume mass and temperature densities
            ibin(:) = get_bin_molno(n)
			if (momentum_outflag .ne. 4) & 
			volume_mass(ibin(1),ibin(2),ibin(3)) = volume_mass(ibin(1),ibin(2),ibin(3)) + mass(n)
			!Note - the streaming term is removed but includes sliding so this must be added back on
			if (peculiar_flag .eq. 0) then
				volume_temperature(ibin(1),ibin(2),ibin(3)) = volume_temperature(ibin(1),ibin(2),ibin(3)) & 
										+ mass(n)*dot_product((v(:,n)+slidev(:,n)),(v(:,n)+slidev(:,n)))
			else
				volume_temperature(ibin(1),ibin(2),ibin(3)) = volume_temperature(ibin(1),ibin(2),ibin(3)) & 
										+ mass(n)*dot_product((v(:,n)-U(:,n)+slidev(:,n)), & 
													          (v(:,n)-U(:,n)+slidev(:,n)))
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
	integer, save		:: sample_count = -1
	
	sample_count = sample_count + 1
	call cumulative_energy(ixyz)
    !print*, 'Total energy = ',sample_count, sum(volume_energy)/(sample_count*np), 0.5*sum(potenergymol(1:np)/np)
	if (sample_count .eq. Nenergy_ave) then
		sample_count = 0

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
			ibin(:) = get_bin_molno(n)
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
!		RECORD CENTRE OF MASS AT LOCATION IN SPACE
! By binning, the centre of mass field in the molecular
! system is recorded and output
!===================================================================================

!Calculate averaged centre of mass in 3D bins
!-----------------------------------------------------------------------------------

subroutine centre_of_mass_averaging(ixyz)
	use module_record, only : centre_of_mass, Ncom_ave, error_abort
	use field_io , only : centre_of_mass_bin_io
	implicit none

	integer				            :: ixyz
	integer,dimension(3)            :: ib
	integer, save		            :: sample_count = -1
	real(kind(0.d0)),dimension(3) 	:: Vbinsize

	sample_count = sample_count + 1
	call cumulative_centre_of_mass(ixyz)

	if (sample_count .eq. Ncom_ave) then
		sample_count = 0

		select case(ixyz)
		case(1:3)
            stop "Centre of Mass is not possible on slice"
		case(4)
            call centre_of_mass_bin_io(centre_of_mass)
            centre_of_mass = 0.d0
		case default
			call error_abort("Error input for Centre of Mass incorrect")
		end select
	endif

end subroutine centre_of_mass_averaging

!-----------------------------------------------------------------------------------
!Add velocities to running total, with 2D slice in ixyz = 1,2,3 or
!in 3D bins when ixyz =4
!-----------------------------------------------------------------------------------

subroutine cumulative_centre_of_mass(ixyz)
	use module_record, only : tag, tether_tags, r, nbins, domain, error_abort, &
                              halfdomain, centre_of_mass, np, get_bin, nhb
    use module_set_parameters, only : mass
    use module_record, only : get_bin_molno
	implicit none

	integer							:: n,ixyz
	integer		,dimension(3)		:: ibin
	real(kind(0.d0)),dimension(3) 	:: COM, mbinsize, bin_centre

	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)
	
	select case(ixyz)
	!COM measurement is a number of 2D slices through the domain
	case(1:3)
        stop "Centre of Mass is not possible on slice"
	!COM measurement for 3D bins throughout the domain
	case(4)
        do n = 1,np
            !Add up current centre_of_mass
            if (any(tag(n).eq.tether_tags)) cycle
            ibin(:) = get_bin_molno(n)
            bin_centre(:) = (ibin(:)-1*nhb(:)-0.5d0)*mbinsize(:)-halfdomain(:)
            COM(:) = mass(n) * (r(:,n) - bin_centre(:))
            !COM(:) = mass(n) * r(:,n)
            centre_of_mass(ibin(1),ibin(2),ibin(3),:) = & 
                centre_of_mass(ibin(1),ibin(2),ibin(3),:) + COM(:)
        enddo
	case default 
		call error_abort("Centre of Mass Binning Error")
	end select
	 
end subroutine cumulative_centre_of_mass




!===================================================================================
!		RECORD PRESSURE AT LOCATION IN SPACE
! Either by Volume Average style binning or virial expression for the whole domain
!===================================================================================

subroutine pressure_averaging(ixyz)
	use field_io, only : virial_stress_io,VA_stress_io,VA_stress_cpol_io
	use module_record
	implicit none

	integer			:: ixyz
	integer, save	:: sample_count = 0

	sample_count = sample_count + 1
	call cumulative_pressure(ixyz,sample_count)
	if (sample_count .eq. Nstress_ave) then
		sample_count = 0
		Pxy = Pxy/dble(Nstress_ave)
		Pxyzero = Pxy		!Update Pxy(0) value

		select case(ixyz)
		case(1)
		    !FULL DOMAIN VIRIAL STRESS CALCULATION
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
			sample_count = sample_count + 1
			if (sample_count .eq. Nvisc_ave) then
				call viscosity_io
				sample_count = 0
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

	integer								:: sample_count
    integer                             :: n,ixyz,jxyz,kxyz
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

	! Add kinetic part of pressure tensor for all molecules
	do n = 1, np

        !Determine bin size (removing halo as vvbins arrays doesn't include halo)
        bin(:) = get_bin_molno(n)-nhb

		!Assign to bins using integer division
		ibin = bin(1) !ceiling((r(1,n)+halfdomain(1))/VAbinsize(1))	!Establish current bin
		if (ibin .gt. imax) cycle ! ibin = maxbin		!Prevents out of range values
		if (ibin .lt. imin) cycle ! ibin = minbin		!Prevents out of range values
		jbin = bin(2) !ceiling((r(2,n)+halfdomain(2))/VAbinsize(2)) 	!Establish current bin
		if (jbin .gt. jmax) cycle ! jbin = maxbin 		!Prevents out of range values
		if (jbin .lt. jmin) cycle ! jbin = minbin		!Prevents out of range values
		kbin = bin(3) !ceiling((r(3,n)+halfdomain(3))/VAbinsize(3)) 	!Establish current bin
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
!       array    - Returned array with length of interaction
!     
! N.B. molecular positon from -halfdomain to +halfdomain       


subroutine pressure_tensor_forces_VA(ri, rj, rF, domain,  & 
                                     nbins, nhb, array,   & 
                                     VA_calcmethod,       &
                                     VA_line_samples)
    use computational_constants_MD, only : cellsidelength
    use module_record, only : get_bin
    implicit none

    integer,intent(in)                                      :: VA_calcmethod
    integer,intent(in),optional                             :: VA_line_samples
    integer,dimension(3),intent(in)                         :: nbins,nhb
	real(kind(0.d0)), dimension(:,:), intent(in)            :: rF
	real(kind(0.d0)),dimension(3), intent(in)		        :: ri, rj, domain
	real(kind(0.d0)),dimension(:,:,:,:,:), intent(inout)	:: array

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
            VA_line_samples_ = 0 !Auto select if zero
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
		array(ibin(1),ibin(2),ibin(3),:,:) = & 
            array(ibin(1),ibin(2),ibin(3),:,:) + rF(:,:)


	!===================Interactions split over 2 or more cells ==========
	case(4:)

		!------Add molecules to bin-----
		!Molecule i and j contribution split between bins

		!-----------Molecule i bin-----------
		array(ibin(1),ibin(2),ibin(3),:,:) = & 
            array(ibin(1),ibin(2),ibin(3),:,:) + 0.5d0*rF(:,:)

		!-----------Molecule j bin-----------
		array(jbin(1),jbin(2),jbin(3),:,:) = & 
            array(jbin(1),jbin(2),jbin(3),:,:) + 0.5d0*rF(:,:)
		
	case default

        print*, 'bindiff = ', ibin, jbin, bindiff
		stop "Harasima Stress AVERAGING ERROR"

	end select

end subroutine pressure_tensor_forces_H


!====================================================================================
! VOLUME AVERAGE CONFIGURATIONAL EXPRESSION
! Author: Edward Smith, David Trevelyan, July 2019
! Linear trajectory path sampled to find approximate values of l_ij 
!(less accurate, but much easier to understand)

subroutine pressure_tensor_forces_VA_trap(ri, rj, rF, VA_line_samples)
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
	!Auto select line segments so minimum of one segment per bin
	if (VA_line_samples .eq. 0) then
		Ns = ceiling(maxval(abs(rij)/VAbinsize))+1
		!print*, "Trap bins", Ns !, rij, rij/VAbinsize, ceiling(maxval(rij/VAbinsize))+1, ceiling(maxval(abs(rij)/VAbinsize))+1
	else
		Ns = VA_line_samples
	endif
	
	! Split line l_ij into segments of size ds
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


            bin(:) = get_bin(rs(:))
			array(bin(1),bin(2),bin(3),:,:) =  &
			array(bin(1),bin(2),bin(3),:,:) + rF(:,:)/real(Ns,kind(0.d0))

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
		array(ibin(1),ibin(2),ibin(3),:,:) =& 
            array(ibin(1),ibin(2),ibin(3),:,:) + rF(:,:)

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
		MLfrac(:,:,:) = sqrt(MLfrac(:,:,:))
		!Normalise to one
		MLfrac(:,:,:) = MLfrac(:,:,:)/magnitude(rij(:))

		!------Add molecules to bin-----
		!Molecule i and j contribution split between bins

		!-----------Molecule i bin-----------
		array(ibin(1),ibin(2),ibin(3),:,:) = & 
            array(ibin(1),ibin(2),ibin(3),:,:) &
			  + rF(:,:)*MLfrac(1,1,1)

		!-----------Molecule j bin-----------
		array(jbin(1),jbin(2),jbin(3),:,:) = & 
            array(jbin(1),jbin(2),jbin(3),:,:) &
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
		MLfrac(:,:,:) = sqrt(MLfrac(:,:,:))

		!Normalise to one
		MLfrac(:,:,:) = MLfrac(:,:,:)/magnitude(rij(:))

		!------Add stress component to bin weighted by line length-----
		!Intermediate bin is either in domain or in halo. 
		!For halo bins the stress is added to the cell the halo represents. 
		!For domain cells it is added directly to that cell

		!-----------Molecule i bin-----------
		array(ibin(1),ibin(2),ibin(3),:,:) = & 
            array(ibin(1),ibin(2),ibin(3),:,:) + rF(:,:)*MLfrac(1,1,1)

		!-----------Intermediate Bin-----------
		!If both bins are in the domain then contribution is added for both molecules
		array(interbin(1,1),interbin(2,1),interbin(3,1),:,:) = &
			array(interbin(1,1),interbin(2,1),interbin(3,1),:,:) &
			+ rF(:,:)*MLfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1))

		!-----------Molecule j bin-----------
		array(jbin(1),jbin(2),jbin(3),:,:) = & 
            array(jbin(1),jbin(2),jbin(3),:,:) &
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
		MLfrac(:,:,:) = sqrt(MLfrac(:,:,:))

		!Normalise to one
		MLfrac(:,:,:) = MLfrac(:,:,:)/magnitude(rij(:))

		!------Add stress component to bin weighted by line length-----
		!Intermediate bin is either in domain or in halo. 
		!For halo bins the stress is added to the cell the halo represents. 
		!For domain cells it is added directly to that cell

		!-----------Molecule i bin-----------
		array(ibin(1),ibin(2),ibin(3),:,:) = & 
            array(ibin(1),ibin(2),ibin(3),:,:) + rF(:,:)*MLfrac(1,1,1)

		!-----------1st Intermediate Bin-----------
		!Intermediate and i bin in domain, j in halo - add intermediate for molecule i
		array(interbin(1,1),interbin(2,1),interbin(3,1),:,:) = &
			array(interbin(1,1),interbin(2,1),interbin(3,1),:,:) &
			+ rF(:,:)*MLfrac(interbindiff(1,1),interbindiff(2,1),interbindiff(3,1))

		!-----------2nd Intermediate Bin-----------
		!Intermediate and i bin in domain, j in halo - add intermediate for molecule i
		array(interbin(1,2),interbin(2,2),interbin(3,2),:,:) = &
			array(interbin(1,2),interbin(2,2),interbin(3,2),:,:) &
			+ rF(:,:)*MLfrac(interbindiff(1,2),interbindiff(2,2),interbindiff(3,2))

		!-----------Molecule j bin-----------
		array(jbin(1),jbin(2),jbin(3),:,:) = & 
            array(jbin(1),jbin(2),jbin(3),:,:) &
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

	integer                         :: i, j, ixyz, jxyz !Define dummy index
	integer							:: icell, jcell, kcell
	integer                         :: icellshift, jcellshift, kcellshift
	integer                         :: cellnp, adjacentcellnp 
	integer							:: molnoi, molnoj
	integer							:: icellmin,jcellmin,kcellmin,icellmax,jcellmax,kcellmax
	type(node), pointer 	        :: oldi, currenti, oldj, currentj

	real(kind(0.d0)),dimension(3)	:: vi_t, cellsperbin
	!real(kind(0.d0)), dimension(3,3):: rF
	real(kind(0.d0)), dimension(1,1):: one
	real(kind(0.d0)), dimension(:,:), allocatable :: rF
	real(kind(0.d0)), dimension(3,1):: rFv

	real(kind(0.d0)), dimension(:,:,:,:,:), allocatable :: zeros
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

	do kcell = kcellmin, kcellmax
	do jcell = jcellmin, jcellmax 
	do icell = icellmin, icellmax 

	
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

				        !Linear magnitude of acceleration for each molecule
				        invrij2 = 1.d0/rij2                 !Invert value
				        accijmag = get_accijmag(invrij2, molnoi, molnoj) !48.d0*(invrij2**7-0.5d0*invrij2**4)
                        call get_outerprod(rij,rij*accijmag,rf)

                        !------------------------------------------------------------
                        ! - Get volume average pressure tensor and pressure heating -
				        if (pressure_outflag .eq. 2 .and. heatflux_outflag .eq. 2) then

                            !Get the velocity, v, at time t 
                            ! ( This is the reason we need to do this after
                            !   the force calculation so we know a(t) ) 
                            vi_t(:) = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)

						    !Select requested configurational line partition methodology
                            do ixyz = 1,3
                                rfv(ixyz,1) = dot_product(rf(ixyz,:),vi_t(:))
                            enddo

                            !Merge both together as line calculation is expensive
                            !especially if non-uniform grid
!                            if (VA_calcmethod .eq. VA_heatflux_calcmethod .and. &
!                                VA_line_samples .eq. VA_heatflux_line_samples) then
!                                if (.not. allocated(zeros)) then
!                                    allocate(zeros(size(rfbin,1), & 
!                                                   size(rfbin,2), & 
!                                                   size(rfbin,3), 1, 1))
!                                    zeros = 0.d0
!                                else
!                                    zeros = 0.d0
!                                endif
                            !    one = 1.d0
!                                call pressure_tensor_forces_VA(ri, rj, one, domain, & 
!                                                               nbins, nhb, zeros,   & 
!                                                               VA_calcmethod,       &
!                                                               VA_line_samples)
!                                do ixyz = 1,3
!                                    rFvbin(:,:,:,ixyz,1) = rFvbin(:,:,:,ixyz,1) + zeros(:,:,:,1,1)*rFv(ixyz,1)
!                                    do jxyz = 1,3
!                                        rFbin(:,:,:,ixyz,jxyz) = rFbin(:,:,:,ixyz,jxyz) + zeros(:,:,:,1,1)*rF(ixyz,jxyz)
!                                    enddo
!                                enddo
!                                
!    
!                            else

						        !Select requested configurational line partition methodology
                                call pressure_tensor_forces_VA(ri, rj, rF, domain,  & 
                                                               nbins, nhb,rfbin,    & 
                                                               VA_calcmethod,       &
                                                               VA_line_samples)

                                call pressure_tensor_forces_VA(ri, rj, rFv, domain, & 
                                                               nbins, nhb, rfvbin,  & 
                                                               VA_heatflux_calcmethod, &
                                                               VA_heatflux_line_samples)
!                            endif
                        else
                            !---------------------------------------
                            ! - Get volume average pressure tensor -
				            if (pressure_outflag .eq. 2) then
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

                                do ixyz = 1,3
                                    rfv(ixyz,1) = dot_product(rf(ixyz,:),vi_t(:))
                                enddo
                                call pressure_tensor_forces_VA(ri, rj, rFv, domain, & 
                                                               nbins, nhb, rfvbin,  & 
                                                               VA_heatflux_calcmethod, &
                                                               VA_heatflux_line_samples)
                            endif
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

	do kcell=int((kmin-1)*cellsperbin(3))+1, int(kmax*cellsperbin(3))
	do jcell=int((jmin-1)*cellsperbin(2))+1, int(jmax*cellsperbin(2))
	do icell=int((imin-1)*cellsperbin(1))+1, int(imax*cellsperbin(1))
	
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
	integer, save	:: sample_count = -1

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

	integer								:: sample_count
    integer                             :: n,ixyz,jxyz,kxyz
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
        call  simulation_compute_energy_VA(1, nbins(1), 1, nbins(2), 1, nbins(3))

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

	! Add kinetic part of heatflux tensor for all molecules
	do n = 1, np

        !Determine bin size (removing halo as evbins arrays doesn't include halo)
        bin(:) = get_bin_molno(n)-nhb

		!Assign to bins using integer division
		ibin = bin(1)	!Establish current bin
		if (ibin .gt. imax) cycle ! ibin = maxbin		!Prevents out of range values
		if (ibin .lt. imin) cycle ! ibin = minbin		!Prevents out of range values
		jbin = bin(2) 	!Establish current bin
		if (jbin .gt. jmax) cycle ! jbin = maxbin 		!Prevents out of range values
		if (jbin .lt. jmin) cycle ! jbin = minbin		!Prevents out of range values
		kbin = bin(3)	!Establish current bin
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

#if __INTEL_COMPILER > 1200
subroutine bforce_pdf_stats
    use boundary_MD
    use statistics_io, only: bforce_pdf_write
    implicit none

	integer, save		:: sample_count=-1

    sample_count = sample_count + 1
    call collect_bforce_pdf_data

    if (sample_count .eq. bforce_pdf_Nave) then
        sample_count = 0
        call bforce_pdf_write
    end if 

end subroutine bforce_pdf_stats
#endif

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


module flux_opt

contains

subroutine cumulative_flux_opt(ri1, ri2, fluxes, quantity, ISR)
	use module_record, only : get_bin, get_crossings, nd, intrinsic_surface_real, &
							ISR_b, iter
    use librarymod, only : CV_surface_flux, imaxloc
    use module_set_parameters, only : mass
	use CV_objects, only : CVcheck_momentum
    implicit none

    real(kind(0.d0)),dimension(3), intent(in)	:: ri1,ri2
	real(kind(0.d0)),dimension(:), allocatable, intent(in) :: quantity
	real(kind(0.d0)),dimension(:,:,:,:,:), allocatable, intent(inout)	:: fluxes

	class(intrinsic_surface_real), intent(in), optional	:: ISR

    logical, save                   :: first_time=.true.

    logical                         :: crossings, crossings_test, changed
	integer							:: jxyz,i,j,k,n,normal,n1,t1,t2,ixyz, Xcount
	integer		,dimension(3)		:: bin1,bin2,cbin,bs
	real(kind(0.d0)),parameter		:: eps = 1e-12
	real(kind(0.d0))        		:: crossdir
    real(kind(0.d0)),dimension(3)	:: ri12, rci

    integer, dimension(:), allocatable :: cbins
    real(kind(0.d0)),dimension(:,:),allocatable :: points
    real(kind(0.d0)),dimension(:,:), allocatable :: rc, rcx, rcy, rcz, dSdr, rc_test

	changed = .false.

    if (present(ISR)) then
        bin1(:) = ISR%get_bin(ri1)
        bin2(:) = ISR%get_bin(ri2)
		n1=ISR%normal
		t1=ISR%ixyz
		t2=ISR%jxyz
    else
        bin1 = get_bin(ri1)
        bin2 = get_bin(ri2)
    endif

	Xcount = 0
    do normal=1,3

        if (present(ISR) .and. (normal .eq. n1)) then

            !Redefine bins here
            call ISR%get_crossings(ri1, ri2, bin1, bin2, normal, rc, crossings, cbins)

			!DEBUG Match to flat surface case DEBUG
			! call get_crossings_test(ri1, ri2, bin1, bin2, normal, rc_test, crossings_test)
			! if (crossings .and. crossings_test) then
				! if (size(rc) .ne. size(rc_test)) then
					! print*, "Sizes not equal", size(rc,2), size(rc_test,2)
					! do i=1,size(rc,2)
						! print*, i, rc(:,i)
					! enddo
					! do i=1,size(rc_test,2)
						! print*, i, rc_test(:,i)
					! enddo
				! else
					! do i=1,size(rc_test,2)
						! if (abs(sum(rc_test(:,i)-rc(:,i))) .gt. 1e-7) then	
							! print'(a,i4, 7f10.5)', "TEST", i, abs(sum(rc_test(:,i)-rc(:,i))), rc_test(:,i), rc(:,i)
						! endif
					! enddo
				! endif
			! else if (crossings) then
				! print*, "No flat crossing", size(rc)
			! else if (crossings_test) then
				! print*, "No intrinsic crossin", size(rc_test)
			! endif
			!DEBUG  Match to flat surface case DEBUG

        else
            call get_crossings(ri1, ri2, bin1, bin2, normal, rc, crossings)
        endif
	
        if (crossings) then

		    !print'(a,i9,8i5,12f10.5)', "Ncrossings ", iter, normal, bin1, bin2, size(rc,2), ri1, rc(:,1), rc(:,2), ri2
            bs = 0
            bs(normal) = 1
            ri12   = ri1 - ri2		        !Molecule i trajectory between t-dt and t
            where (abs(ri12) .lt. 0.000001d0) ri12 = 0.000001d0

            if (present(ISR) .and. (normal .eq. n1)) then
                if (allocated(points)) deallocate(points)
                allocate(points(size(rc,2), nd))
                points(:,1) = rc(1,:)
                points(:,2) = rc(2,:)
                points(:,3) = rc(3,:)
                call ISR%get_surface_derivative(points, dSdr)
            endif

            do i=1,size(rc,2)
                rci = rc(:,i)
                rci(normal) = rci(normal) + eps
				if (present(ISR)) then
					!Whole point of cbins return was to avoid another expensive get_bin call
					!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
					cbin(:) = ISR%get_bin(rci)
				else
					cbin(:) = get_bin(rci)
				endif

                if (present(ISR) .and. (normal .eq. n1)) then
                    !if (cbins(i)+1 .ne. cbin(normal)) print'(a,i9,6i5,9f10.5)', "cross no ", iter, i, & 
                    !                                size(rc,2), cbin, cbins(i)+1, ri1, rci, ri2
                    !Set normal bin based on solution from bilinear crossinfs
                    cbin(normal) = cbins(i)+1
					!if (cbins(i)+1 .ne. cbin(normal)) print*, cbin(normal), cbins(i)+1 
                    crossdir = sign(1.d0,(ri12(n1) - ri12(t1)*dSdr(i,1) - ri12(t2)*dSdr(i,2)))
				else
                    crossdir  = sign(1.d0, ri12(normal))
                endif

    		    !if (size(rc,2) .gt. 1) print'(a,i9,6i5,9f10.5)', "cross no ", iter, i, size(rc,2), & 
                !                       ISR%get_bin(rci, nbins, nhb), cbins(i), ri1, rci, ri2
                !if (cbin(normal) .ne. cbins(normal)) stop "Error"
				!print*, "CBINS", rci, cbin, bin1, bin2
                fluxes(cbin(1),cbin(2),cbin(3),:,normal) = & 
                    fluxes(cbin(1),cbin(2),cbin(3),:,normal) + crossdir*quantity(:)
                fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),:,normal+3) = & 
                    fluxes(cbin(1),cbin(2),cbin(3),:,normal)

				! if (present(ISR)) then
				! if ((cbin(1) .ge. CVcheck_momentum%debug_CV(1) .and. & 
					 ! cbin(1) .lt. CVcheck_momentum%debug_CV(1)+CVcheck_momentum%debug_CV_range(1)) .and. & 
				    ! (cbin(2) .ge. CVcheck_momentum%debug_CV(2) .and. &
					 ! cbin(2) .lt. CVcheck_momentum%debug_CV(2)+CVcheck_momentum%debug_CV_range(2)) .and. & 
				    ! (cbin(3) .ge. CVcheck_momentum%debug_CV(3) .and. & 
					 ! cbin(3) .lt. CVcheck_momentum%debug_CV(3)+CVcheck_momentum%debug_CV_range(3))) then
					! print'(a,i8,2i3, 3i5,6f10.5,7f10.4)', "t crossing ", iter, i, normal, cbin, rci, & 
								! 0.25*quantity(1)/(ISR%binsize(2)*ISR%binsize(3)), &
								! 0.25*quantity(1)/(ISR%binsize(1)*ISR%binsize(3)), &
								! 0.25*quantity(1)/(ISR%binsize(1)*ISR%binsize(2)),  &
								! 0.25d0*fluxes(cbin(1),cbin(2),cbin(3),1,1)/(ISR%binsize(2)*ISR%binsize(3)), & 
								! 0.25d0*fluxes(cbin(1),cbin(2),cbin(3),1,2)/(ISR%binsize(1)*ISR%binsize(3)), & 
								! 0.25d0*fluxes(cbin(1),cbin(2),cbin(3),1,3)/(ISR%binsize(1)*ISR%binsize(2)), &
								! 0.25d0*fluxes(cbin(1),cbin(2),cbin(3),1,4)/(ISR%binsize(2)*ISR%binsize(3)), & 
								! 0.25d0*fluxes(cbin(1),cbin(2),cbin(3),1,5)/(ISR%binsize(1)*ISR%binsize(3)), & 
								! 0.25d0*fluxes(cbin(1),cbin(2),cbin(3),1,6)/(ISR%binsize(1)*ISR%binsize(2)), &
						   ! (fluxes(cbin(1),cbin(2),cbin(3),1,1)-fluxes(cbin(1),cbin(2),cbin(3),1,4))/ISR%binsize(1) &
						  ! +(fluxes(cbin(1),cbin(2),cbin(3),1,2)-fluxes(cbin(1),cbin(2),cbin(3),1,5))/ISR%binsize(2) &
						  ! +(fluxes(cbin(1),cbin(2),cbin(3),1,3)-fluxes(cbin(1),cbin(2),cbin(3),1,6))/ISR%binsize(3)
					! write(586410,*) iter, normal, ri1, rci, ri2, quantity
					! Xcount = Xcount + 1
					! write(292847,*) iter, cbin(1),cbin(2),cbin(3), ISR_b%binsize, ISR_b%indices_to_points(cbin(1),cbin(2),cbin(3))
					! changed = .true.
				! endif
				! if ((cbin(1)-bs(1) .ge. CVcheck_momentum%debug_CV(1) .and. & 
					 ! cbin(1)-bs(1) .lt. CVcheck_momentum%debug_CV(1)+CVcheck_momentum%debug_CV_range(1)) .and. & 
				    ! (cbin(2)-bs(2) .ge. CVcheck_momentum%debug_CV(2) .and. & 
					 ! cbin(2)-bs(2) .lt. CVcheck_momentum%debug_CV(2)+CVcheck_momentum%debug_CV_range(2)) .and. & 
				    ! (cbin(3)-bs(3) .ge. CVcheck_momentum%debug_CV(3) .and. & 
					 ! cbin(3)-bs(3) .lt. CVcheck_momentum%debug_CV(3)+CVcheck_momentum%debug_CV_range(3))) then
					! print'(a,i8,2i3, 3i5,6f10.5,7f10.4)', "b crossing ", iter, i, normal, cbin-bs, rci, & 
								! 0.25*quantity(1)/(ISR%binsize(2)*ISR%binsize(3)), &
								! 0.25*quantity(1)/(ISR%binsize(1)*ISR%binsize(3)), & 
								! 0.25*quantity(1)/(ISR%binsize(1)*ISR%binsize(2)),  &
								! 0.25d0*fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,1)/(ISR%binsize(2)*ISR%binsize(3)), &
								! 0.25d0*fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,2)/(ISR%binsize(1)*ISR%binsize(3)), &
								! 0.25d0*fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,3)/(ISR%binsize(1)*ISR%binsize(2)), &
								! 0.25d0*fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,4)/(ISR%binsize(2)*ISR%binsize(3)), &
								! 0.25d0*fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,5)/(ISR%binsize(1)*ISR%binsize(3)), &
								! 0.25d0*fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,6)/(ISR%binsize(1)*ISR%binsize(2)), &
						   ! (fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,1)& 
						   ! -fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,4))/ISR%binsize(1) &
						  ! +(fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,2)& 
						   ! -fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,5))/ISR%binsize(2) &
						  ! +(fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,3)& 
						   ! -fluxes(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),1,6))/ISR%binsize(3)
					! write(586410,*) iter, normal, ri1, rci, ri2, quantity
					! Xcount = Xcount + 1
					! write(292847,*) iter, cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3), ISR_b%binsize, & 
						! ISR_b%indices_to_points(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3))
					! changed = .true.
				! endif
				! endif

				! if (changed) then
					! print*, ri1, ri2, bin1, bin2, rci
				! endif


            enddo
            deallocate(rc)
        endif

	enddo

	!if (changed .and. Xcount .ne. 2) then
		!print*, "Crossing count", Xcount, ri1, ri2
		! do i =136,137
		! do j =18,19
		! do k =18,19
			! print'(a,i8,3i5, 6f11.4)', "fluxes ", iter, i, j, k, &
				! 0.25d0*fluxes(i,j,k,1,1)/(ISR%binsize(2)*ISR%binsize(3)), &
				! 0.25d0*fluxes(i,j,k,1,2)/0.069480940251555504, &
				! 0.25d0*fluxes(i,j,k,1,3)/0.069480940251555504, &
				! 0.25d0*fluxes(i,j,k,1,4)/(ISR%binsize(2)*ISR%binsize(3)), &
				! 0.25d0*fluxes(i,j,k,1,5)/0.069480940251555504, &
				! 0.25d0*fluxes(i,j,k,1,6)/0.069480940251555504	
		! enddo
		! enddo
		! enddo
	!endif

end subroutine cumulative_flux_opt


subroutine cumulative_mass_flux_opt
	use module_record, only : cluster_analysis_outflag, np, r, v, & 
                              delta_t, intrinsic_interface_outflag, & 
							  mass_flux, ISR_mdt
    use module_set_parameters, only : mass
    implicit none

    logical                         :: use_bilinear

    integer :: n
    real(kind(0.d0)),dimension(3)	:: ri1,ri2
	real(kind(0.d0)),dimension(:), allocatable :: quantity
	real(kind(0.d0)),dimension(:,:,:,:,:), allocatable	:: fluxes

    if (cluster_analysis_outflag .eq. 1 .and.  & 
        any(intrinsic_interface_outflag .eq. (/1,2/))) then
        use_bilinear = .true.
    else
        use_bilinear = .false.
    endif

    allocate(quantity(1))
    allocate(fluxes(size(mass_flux,1), size(mass_flux,2), & 
                    size(mass_flux,3), 1, size(mass_flux,4)))
    fluxes(:,:,:,1,:) = mass_flux(:,:,:,:)

	do n = 1,np
		ri1(:) = r(:,n) 							!Molecule i at time t
		ri2(:) = r(:,n)	- delta_t*v(:,n)			!Molecule i at time t-dt
        quantity(1) = mass(n)
		if (use_bilinear) then
			call cumulative_flux_opt(ri1, ri2, fluxes, quantity, ISR_mdt)
		else
			call cumulative_flux_opt(ri1, ri2, fluxes, quantity)
		endif
    enddo
    mass_flux(:,:,:,:) = fluxes(:,:,:,1,:)

end subroutine cumulative_mass_flux_opt

subroutine cumulative_momentum_flux_opt(r_,v_,momentum_flux_,notcrossing)
	use module_record, only : cluster_analysis_outflag, np, & 
                             delta_t, intrinsic_interface_outflag, &
							 ISR_mdt
   use module_set_parameters, only : mass
   implicit none

	real(kind(0.d0)),dimension(:,:),allocatable,intent(in) 			:: r_,v_
	real(kind(0.d0)),dimension(:,:,:,:,:),allocatable,intent(inout) :: momentum_flux_
	integer,dimension(:),allocatable,intent(out),optional			:: notcrossing

   logical                         :: use_bilinear

   integer :: n
   real(kind(0.d0)),dimension(3)	:: ri1,ri2
	real(kind(0.d0)),dimension(:), allocatable :: quantity
	real(kind(0.d0)),dimension(:,:,:,:,:), allocatable	:: fluxes

	!Allocate array if required and assume all are not crossing
	if (present(notcrossing)) then
		if (allocated(notcrossing)) deallocate(notcrossing)
		allocate(notcrossing(np)); notcrossing = 1
	endif

   if (cluster_analysis_outflag .eq. 1 .and.  & 
       any(intrinsic_interface_outflag .eq. (/1,2/))) then
       use_bilinear = .true.
   else
       use_bilinear = .false.
   endif

   allocate(quantity(3))

	do n = 1,np
		ri1(:) = r_(:,n) 							!Molecule i at time t
		ri2(:) = r_(:,n)	- delta_t*v_(:,n)			!Molecule i at time t-dt
       quantity(:) = v_(:,n)
		if (use_bilinear) then
			call cumulative_flux_opt(ri1, ri2, momentum_flux_, quantity, ISR_mdt)
		else
			call cumulative_flux_opt(ri1, ri2, momentum_flux_, quantity)
		endif
		!Record mask of molecules which are currently crossing
		if (present(notcrossing)) notcrossing(n) = 0
   enddo

end subroutine cumulative_momentum_flux_opt


end module flux_opt

!===================================================================================
! Stresses over each of the six surfaces of the cuboid

subroutine control_volume_stresses_opt(fij, ri, rj)
    use module_record, only : ISR, Pxyface, void, initialstep, &
							  CVforce_flag, CVforce_starttime, &
							  cluster_analysis_outflag, iter, &
							  intrinsic_interface_outflag
	use CV_objects, only : CV_constraint
	use flux_opt, only : cumulative_flux_opt
    implicit none

	real(kind(0.d0)),intent(in),dimension(3)	:: ri,rj,fij

    !integer, save :: Ncount
    integer :: ixyz, bini(3), binj(3)
    real(kind(0.d0)) :: zeromode, window, zl(3)

	real(kind(0.d0)),dimension(:), allocatable :: quantity

    !print*, "ISR_mdt%coeff", maxval(ISR_mdt%coeff), sum(ISR_mdt%coeff)

    !Limit calc to part of domain for debugging
    !window = 3.d0
    !call ISR_mdt%get_zero_mode(zeromode)
    !zl = (/ zeromode, 0.d0, 0.d0 /)
    !do ixyz=1,3
        !if (bini(ixyz) .ge. CV_constraint%debug_CV(ixyz) .and. & 
        !    bini(ixyz) .lt. CV_constraint%debug_CV(ixyz)+CV_constraint%debug_CV_range(ixyz)) then
        !    cycle
        !else
        !    return
        !endif
     !   if (ri(ixyz) .lt. -window+zl(ixyz) .or. & 
     !       rj(ixyz) .lt. -window+zl(ixyz) .or. &
     !       ri(ixyz) .gt.  window+zl(ixyz) .or. &
     !       rj(ixyz) .gt.  window+zl(ixyz)) return
    !enddo
    !print*,  get_bin(ri), get_bin(rj)

	allocate(quantity(3))
	quantity(:) = 2.d0*fij(:)
	if (cluster_analysis_outflag .eq. 1 .and.  & 
	   any(intrinsic_interface_outflag .eq. (/1,2/))) then! .and. &
        !maxval(ISR_mdt%coeff) .gt. 1e-4) then
		call cumulative_flux_opt(ri, rj, Pxyface, quantity, ISR)
	else
		call cumulative_flux_opt(ri, rj, Pxyface, quantity)
	endif

	!Add instantanous stress to CV record
	if (CVforce_flag .ne. VOID .and. iter-initialstep+1 .ge. CVforce_starttime) then
		CV_constraint%Pxy = Pxyface
	endif

    !Debug count interactions
    !Ncount = Ncount + 1
    !if (mod(Ncount,10000) .eq. 0) print*, "cumulative_flux_opt", Ncount, fij, ri, rj

end subroutine control_volume_stresses_opt

!-----------------------------------------------------------------------------------
! Control Volume mass continuity
!-----------------------------------------------------------------------------------

subroutine mass_flux_averaging(flag)
	!use field_io, only : mass_flux_io
	use module_record, only : Nmflux_ave, domain, nbins, nhb, & 
                              thermstattop, thermstatbottom, &
                              iter, irank, mass_flux 
	use CV_objects, only : CV_debug, CVcheck_mass, check_error_mass !,CV_sphere_mass
    use boundary_MD, only : specular_wall
	use flux_opt, only : cumulative_mass_flux_opt
	implicit none

	integer			                :: flag
    integer,dimension(3)            :: thermbinstop,thermbinsbot
    real(kind(0.d0)),dimension(3)   :: mbinsize
	integer, save	                :: sample_count = 0
    integer,dimension(3):: skipbinstop,skipbinsbot

    double precision :: t1, t2
    double precision, save :: tsum=0.d0

	!Only average if mass averaging turned on
	if (flag .eq. 0) return

    !call cpu_time(t1)
	call cumulative_mass_flux_opt()

    !call cpu_time(t2)
    !tsum = tsum + (t2-t1)
    !print*, "Tsum = ", iter, tsum
	sample_count = sample_count + 1
	if (sample_count .eq. Nmflux_ave) then
		if (CV_debug .ne. 0) then
			call mass_flux_io()
			sample_count = 0
			mass_flux = 0
			call mass_snapshot()
    		mbinsize(:) = domain(:) / nbins(:)
            skipbinstop = ceiling((thermstattop + specular_wall)/mbinsize)
            skipbinsbot = ceiling((thermstatbottom + specular_wall)/mbinsize)
            
            !E.S. this causes a compiler seg fault for 
            !     ifort version 13.0.1 which is fixed by 
            !     replacing 
            !     "call CVcheck_mass%check_error( ... "   
            !     with 
            !     "check_error_mass(CVcheck_mass, ... "
            call check_error_mass(CVcheck_mass, &
                                  1+nhb(1)+skipbinsbot(1),nbins(1)+nhb(1)-skipbinstop(1), & 
								  1+nhb(2)+skipbinsbot(2),nbins(2)+nhb(2)-skipbinstop(2), & 
								  1+nhb(3)+skipbinsbot(3),nbins(3)+nhb(3)-skipbinstop(3),iter,irank)
			!call CV_sphere_mass%check_error(1,1,1,1,1,1,iter,irank)
	    endif
	endif

end subroutine mass_flux_averaging

!!===================================================================================
!! Mass Flux over a surface of a bin
!! Includes all intermediate bins

subroutine cumulative_mass_flux
	use module_record
    use librarymod, only : CV_surface_flux, imaxloc
    use module_set_parameters, only : mass
    !use CV_objects, only : CV_sphere_mass
    implicit none

	integer							:: jxyz,i,j,k,n
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	real(kind(0.d0)),dimension(3)	:: mbinsize,crossface
	real(kind(0.d0)),dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
	real(kind(0.d0)),dimension(1)	:: quantity
	real(kind(0.d0)),dimension(6)	:: CV
	real(kind(0.d0)),dimension(1,6)	:: fluxes

	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	do n = 1,np

		ri1(:) = r(:,n) 							!Molecule i at time t
		ri2(:) = r(:,n)	- delta_t*v(:,n)			!Molecule i at time t-dt
		ri12   = ri1 - ri2							!Molecule i trajectory between t-dt and t
		where (abs(ri12) .lt. 0.000001d0) ri12 = 0.000001d0
		
		!call CV_sphere_mass%Add_spherical_CV_fluxes(ri2,ri1)

		!Assign to bins before and after using integer division
        ibin1(:) =  get_bin(ri1)
        ibin2(:) =  get_bin(ri2)

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

                CV(1:3) = binbot(:); CV(4:6) = bintop(:)
                quantity(1) = mass(n)
                call CV_surface_flux(ri1, ri2, CV, quantity, 1, fluxes)

				jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer

				!Add Mass flux over face
                mass_flux(cbin(1),cbin(2),cbin(3),:) = & 
                    mass_flux(cbin(1),cbin(2),cbin(3),:) + fluxes(1,:) * abs(crossface(jxyz))

			enddo
			enddo
			enddo

		endif

	enddo

end subroutine cumulative_mass_flux

!!===================================================================================
!! Control Volume snapshot of the mass in a given bin

subroutine mass_snapshot()
	use module_record
	use field_io, only : mass_bin_io
	use CV_objects, only : CVcheck_mass, CV_debug!, CV_sphere_mass
    use module_set_parameters, only : mass
	implicit none

	integer										    :: n
	integer		,dimension(3)					    :: ibin
	real(kind(0.d0)),dimension(:,:,:)  ,allocatable	:: volume_mass_temp

	!Allocate temporary array for mass and momentum in volume
	allocate(volume_mass_temp(nbinso(1),nbinso(2),nbinso(3)))

	!Reset Control Volume momentum 
	volume_mass_temp = 0
	do n = 1,np

		!Add up current volume momentum densities
		ibin(:) = get_bin_molno(n)
		volume_mass_temp(ibin(1),ibin(2),ibin(3)) = & 
            volume_mass_temp(ibin(1),ibin(2),ibin(3)) + mass(n)
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
    use librarymod, only : imaxloc, CV_surface_flux
#if ASSMBLY_HEAVISIDES
    use librarymod, only : heaviside  =>  heaviside_a1
#else
    use librarymod, only : heaviside
#endif
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
	real(kind(0.d0)),dimension(6)	:: CV
	real(kind(0.d0)),dimension(3,6)	:: fluxes

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
			ri2(:) = r_(:,n) - delta_t*velvect			!Molecule i at time t-dt
			ri12   = ri1 - ri2							!Molecule i trajectory between t-dt and t
			where (abs(ri12) .lt. 0.000001d0) ri12 = 0.000001d0

			!call CV_sphere_momentum%Add_spherical_CV_fluxes(velvect,ri2,ri1)

			!Assign to bins before and after using integer division
            ibin1(:) =  get_bin(ri1)
            ibin2(:) =  get_bin(ri2)

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

                    CV(1:3) = binbot(:); CV(4:6) = bintop(:)

					jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer

					!Calculate velocity at time of intersection
					!crosstime = (r_(jxyz,n) - rplane)/v_(jxyz,n)
					velvect(:) = v_(:,n) !- a(:,n) * crosstime

                    call CV_surface_flux(ri1, ri2, CV, velvect, 3, fluxes)

					!Add Momentum flux over face
                    momentum_flux_(cbin(1),cbin(2),cbin(3),:,:) = & 
                        momentum_flux_(cbin(1),cbin(2),cbin(3),:,:) & 
                            + fluxes(:,:) * abs(crossface(jxyz))

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
	use flux_opt, only : cumulative_momentum_flux_opt
	use cumulative_momentum_flux_mod, only : cumulative_momentum_flux
	use CV_objects, only : CV_debug, CVcheck_momentum!, CV_constraint!, CV_sphere_momentum
    use boundary_MD, only : specular_wall
	implicit none

	integer,intent(in)	:: flag
	integer				::icell,jcell,kcell,n
    integer,dimension(3):: skipbinstop,skipbinsbot
	real(kind(0.d0)),dimension(3)	:: mbinsize
	integer, save		:: sample_count = 0

	if (flag .eq. 0) return

	call cumulative_momentum_flux_opt(r,v,momentum_flux)
	sample_count = sample_count + 1
	if (sample_count .eq. Nvflux_ave) then

		select case(flag)
		case(1:3)
			!MOP momentum flux and stresses
			call MOP_stress_io(flag)
			Pxy_plane = 0.d0
		case(4)
			!CV momentum flux and stress
			call momentum_flux_io()
			momentum_flux = 0.d0
			call momentum_snapshot()
			if (external_force_flag .ne. 0 .or. & 
				ensemble .eq. tag_move     .or. & 
				CVforce_flag .ne. VOID) then
				call external_force_io()
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
		call surface_stress_io()
		Pxyface = 0.d0
		!Debug flag to check CV conservation in code
        if (CV_debug .ne. 0) then

            !Currently, CV check will not work for thermostatted molecules
            !so exclude these from checked region
    		mbinsize(:) = domain(:) / nbins(:)
            if (ensemble .eq. 6) then
                skipbinstop = ceiling((thermstattop + specular_wall)/mbinsize)
                skipbinsbot = ceiling((thermstatbottom + specular_wall)/mbinsize)
            else
                skipbinstop = 0
                skipbinsbot = 0
            endif
    	    call CVcheck_momentum%check_error(1+nhb(1)+skipbinsbot(1), & 
                                              nbins(1)+nhb(1)-skipbinstop(1), & 
    										  1+nhb(2)+skipbinsbot(2), & 
                                              nbins(2)+nhb(2)-skipbinstop(2), & 
    										  1+nhb(3)+skipbinsbot(3), & 
                                              nbins(3)+nhb(3)-skipbinstop(3),iter,irank)
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
        ibin(:) =  get_bin_molno(n)
		!ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb

		volume_mass_temp(ibin(1),ibin(2),ibin(3)) = & 
            volume_mass_temp(ibin(1),ibin(2),ibin(3)) + mass(n)
		volume_momentum_temp(ibin(1),ibin(2),ibin(3),:) = & 
            volume_momentum_temp(ibin(1),ibin(2),ibin(3),:) + mass(n)*v(:,n)
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
    use librarymod, only : imaxloc, CV_surface_flux
#if ASSMBLY_HEAVISIDES
    use librarymod, only : heaviside  =>  heaviside_a1
#else
    use librarymod, only : heaviside
#endif
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
	real(kind(0.d0)),dimension(1)	:: quantity
	real(kind(0.d0)),dimension(6)	:: CV
	real(kind(0.d0)),dimension(1,6)	:: fluxes

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
			where (abs(ri12) .lt. 0.000001d0) ri12 = 0.000001d0

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
                    CV(1:3) = binbot(:); CV(4:6) = bintop(:)

					jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer

					!Calculate velocity at time of intersection
					velvect(:) = v_(:,n) + 0.5d0*a(:,n)*delta_t
					energy = 0.5d0 * ( mass(n)*dot_product(velvect,velvect) + potenergymol(n))

					!Change in velocity at time of crossing is not needed as velocity assumed constant 
					!for timestep and changes when forces are applied.
                    quantity(1) = energy
                    call CV_surface_flux(ri1, ri2, CV, quantity, 1, fluxes)
    
					!Add Energy flux over face
                    energy_flux_(cbin(1),cbin(2),cbin(3),:) = &
                        energy_flux_(cbin(1),cbin(2),cbin(3),:) + fluxes(1,:)*abs(crossface(jxyz))

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
    use boundary_MD, only : specular_wall
	implicit none

	integer,intent(in)	:: flag
	integer, save		:: sample_count = -1

    integer,dimension(3):: skipbinstop,skipbinsbot
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
                skipbinstop = ceiling((thermstattop + specular_wall)/ebinsize)
                skipbinsbot = ceiling((thermstatbottom + specular_wall)/ebinsize)
            else
                skipbinstop = 0
                skipbinsbot = 0
            endif
		    call CVcheck_energy%check_error(1+nhb(1)+skipbinstop(1), & 
                                            nbins(1)+nhb(1)-skipbinstop(1), & 
											1+nhb(2)+skipbinstop(2), & 
                                            nbins(2)+nhb(2)-skipbinstop(2), & 
											1+nhb(3)+skipbinstop(3), & 
                                            nbins(3)+nhb(3)-skipbinstop(3),iter,irank)
	   endif
	endif

end subroutine energy_flux_averaging

!===================================================================================
! Control Volume snapshot of momentum in a given bin

subroutine energy_snapshot
	use field_io, only : energy_bin_io
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
        ibin(:) = get_bin_molno(n)
		!ibin(:) = ceiling((r(:,n)+halfdomain(:))/mbinsize(:)) + nhb(:)
		velvect(:) = v(:,n) + 0.5d0*a(:,n)*delta_t
		energy = 0.5d0 * ( mass(n)*dot_product(velvect,velvect) + potenergymol(n))

		volume_energy_temp(ibin(1),ibin(2),ibin(3)) = volume_energy_temp(ibin(1),ibin(2),ibin(3)) + energy
	enddo

	binvolume = (domain(1)/nbins(1))*(domain(2)/nbins(2))*(domain(3)/nbins(3))
	volume_energy_temp = volume_energy_temp/(binvolume)

	!Output Control Volume momentum change and fluxes
	call energy_bin_io(volume_energy_temp,'snap')

	deallocate(volume_energy_temp)

end subroutine energy_snapshot



!===================================================================================
! Density of molecules found on the surface of a bin 
! Includes all intermediate bins, methodology from 
! " A technique for the calculation of mass, energy, and momentum densities
!   at planes in molecular dynamics simulations"
!  By Peter J. Daivis, Karl P. Travis, and B. D. Todd


subroutine surface_density_averaging(flag)
	!use field_io, only : mass_flux_io
	use module_record
	implicit none

	integer			                :: flag
    real(kind(0.d0)),dimension(3)   :: mbinsize
	integer, save	                :: sample_count = 0

	!Only average if mass averaging turned on
	if (flag .eq. 0) return

	call cumulative_surface_density()
	sample_count = sample_count + 1
	if (sample_count .eq. Nsurfm_ave) then
		call surface_density_io()
		sample_count = 0
		surface_density = 0
	endif

end subroutine surface_density_averaging


! ----------------------------------------------------------------------------------

subroutine cumulative_surface_density
	use module_record
    use librarymod, only : imaxloc,  CV_surface_crossing
    use module_set_parameters, only : mass
    !use CV_objects, only : CV_sphere_mass
    implicit none

	integer							:: jxyz,i,j,k,n
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	real(kind(0.d0)),dimension(3)	:: mbinsize,crossface,velvect
	real(kind(0.d0)),dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
	real(kind(0.d0)),dimension(6)	:: CV, onface
	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	do n = 1,np

		ri1(:) = r(:,n) 							!Molecule i at time t
		ri2(:) = r(:,n)	- delta_t*v(:,n)			!Molecule i at time t-dt
		ri12   = ri1 - ri2							!Molecule i trajectory between t-dt and t
		where (abs(ri12) .lt. 0.000001d0) ri12 = 0.000001d0
		
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
                CV(1:3) = binbot(:); CV(4:6) = bintop(:)
                call CV_surface_crossing(ri1, ri2, CV, onface)

				jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer
				velvect(:) = v(:,n) !- a(:,n) * crosstime

				!Add Density on surface for molecules crossing face
				surface_density(cbin(1),cbin(2),cbin(3),1) = & 
					surface_density(cbin(1),cbin(2),cbin(3),1) & 
				      + mass(n)*nint(dble(onface(1))*crossface(jxyz))/abs(velvect(1))
				surface_density(cbin(1),cbin(2),cbin(3),2) = & 
					surface_density(cbin(1),cbin(2),cbin(3),2) & 
				      + mass(n)*nint(dble(onface(2))*crossface(jxyz))/abs(velvect(2))
				surface_density(cbin(1),cbin(2),cbin(3),3) = & 
					surface_density(cbin(1),cbin(2),cbin(3),3) &
				      + mass(n)*nint(dble(onface(3))*crossface(jxyz))/abs(velvect(3))
				surface_density(cbin(1),cbin(2),cbin(3),4) = & 
					surface_density(cbin(1),cbin(2),cbin(3),4) &
				      + mass(n)*nint(dble(onface(4))*crossface(jxyz))/abs(velvect(1))
				surface_density(cbin(1),cbin(2),cbin(3),5) = & 
					surface_density(cbin(1),cbin(2),cbin(3),5) &
				      + mass(n)*nint(dble(onface(5))*crossface(jxyz))/abs(velvect(2))
				surface_density(cbin(1),cbin(2),cbin(3),6) = & 
					surface_density(cbin(1),cbin(2),cbin(3),6) &
				      + mass(n)*nint(dble(onface(6))*crossface(jxyz))/abs(velvect(3))
				      

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
#if ASSMBLY_HEAVISIDES
    use librarymod, only : heaviside  =>  heaviside_a1
#else
    use librarymod, only : heaviside
#endif
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

subroutine control_volume_stresses(fij, ri, rj)
    use module_record
	use CV_objects, only : CV_debug,CV_constraint
    use librarymod, only :  CV_surface_crossing
    implicit none


	real(kind(0.d0)),intent(in),dimension(3)	:: ri,rj,fij

	integer							:: i,j,k,ixyz,n,face
	integer,dimension(3)			:: cbin, ibin, jbin
    real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb,fijvi,fijvj
	real(kind(0.d0)),dimension(3)	:: rij,fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb,velvect
	real(kind(0.d0)),dimension(3)	:: Fbinsize, bintop, binbot, vi_t,vj_t,vi_tmdt,vj_tmdt
	real(kind(0.d0)),dimension(6)	:: CV, onface

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

		!Calculate the plane intersect of trajectory with surfaces of the cube
        CV(1:3) = binbot(:); CV(4:6) = bintop(:)
        call CV_surface_crossing(ri, rj, CV, onface)
        onface = 2.d0 * onface

		!Stress acting on face over volume
!        do n =1,6
!    		Pxyface(cbin(1),cbin(2),cbin(3),:,n) = Pxyface(cbin(1),cbin(2),cbin(3),:,n) + fij(:)*dble(onface(n))
!		    if (CVforce_flag .ne. VOID .and. iter-initialstep+1 .ge. CVforce_starttime) then
!        		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,n) = & 
!				    CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,n) + fij(:)*dble(onface(n))
!            endif
!        enddo

		Pxyface(cbin(1),cbin(2),cbin(3),:,1) = Pxyface(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*dble(onface(1))
		Pxyface(cbin(1),cbin(2),cbin(3),:,2) = Pxyface(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*dble(onface(2))
		Pxyface(cbin(1),cbin(2),cbin(3),:,3) = Pxyface(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*dble(onface(3))
		Pxyface(cbin(1),cbin(2),cbin(3),:,4) = Pxyface(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*dble(onface(4))
		Pxyface(cbin(1),cbin(2),cbin(3),:,5) = Pxyface(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*dble(onface(5))
		Pxyface(cbin(1),cbin(2),cbin(3),:,6) = Pxyface(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*dble(onface(6))

		!Add instantanous stress to CV record
		if (CVforce_flag .ne. VOID .and. iter-initialstep+1 .ge. CVforce_starttime) then
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,1) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*dble(onface(1))
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,2) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*dble(onface(2))
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,3) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*dble(onface(3))
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,4) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*dble(onface(4))
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,5) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*dble(onface(5))
    		CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,6) = & 
				CV_constraint%Pxy(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*dble(onface(6))
		endif

		!Force applied to volume
		fsurface(:) = 0.d0
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onface(1) - onface(4))
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onface(2) - onface(5))
		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onface(3) - onface(6))
		volume_force(cbin(1),cbin(2),cbin(3),:,1) = volume_force(cbin(1),cbin(2),cbin(3),:,1) + fsurface*delta_t

	enddo
	enddo
	enddo


end subroutine control_volume_stresses




!===================================================================================
! Stresses times velocity over each of the six surfaces of the cuboid

subroutine control_volume_power(fij, ri, rj, vi_t)
    use module_record
	use CV_objects, only : CV_debug,CV_constraint
    use librarymod, only : CV_surface_crossing
    implicit none


	real(kind(0.d0)),intent(in),dimension(3)	:: ri,rj,fij
	real(kind(0.d0)),dimension(3),intent(in)	:: vi_t


	integer							:: i,j,k,ixyz,face
	integer,dimension(3)			:: cbin, ibin, jbin
    real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb,fijvi,fijvj
	real(kind(0.d0)),dimension(3)	:: rij,fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb,velvect
	real(kind(0.d0)),dimension(3)	:: Fbinsize, bintop, binbot
	real(kind(0.d0)),dimension(6)	:: CV, onface

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

		!Calculate the plane intersect of trajectory with surfaces of the cube
        CV(1:3) = binbot(:); CV(4:6) = bintop(:)
        call CV_surface_crossing(ri, rj, CV, onface)
        onface = 2.d0 * onface

		!Stress work acting on face over volume - add current value 
        ! to use in trapizium rule later 
		fijvi = dot_product(fij,vi_t)

		Pxyvface(cbin(1),cbin(2),cbin(3),1) = & 
			Pxyvface(cbin(1),cbin(2),cbin(3),1) + fijvi*onface(1)
		Pxyvface(cbin(1),cbin(2),cbin(3),2) = & 
			Pxyvface(cbin(1),cbin(2),cbin(3),2) + fijvi*onface(2)
		Pxyvface(cbin(1),cbin(2),cbin(3),3) = &
			Pxyvface(cbin(1),cbin(2),cbin(3),3) + fijvi*onface(3)
		Pxyvface(cbin(1),cbin(2),cbin(3),4) = &
			Pxyvface(cbin(1),cbin(2),cbin(3),4) + fijvi*onface(4)
		Pxyvface(cbin(1),cbin(2),cbin(3),5) = &
			Pxyvface(cbin(1),cbin(2),cbin(3),5) + fijvi*onface(5)
		Pxyvface(cbin(1),cbin(2),cbin(3),6) = &
			Pxyvface(cbin(1),cbin(2),cbin(3),6) + fijvi*onface(6)

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
#if ASSMBLY_HEAVISIDES
    use librarymod, only : heaviside  =>  heaviside_a1
#else
    use librarymod, only : heaviside
#endif
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




module cubic_surface_CV

    contains


!========================================================================================================
!    !CV objects
!========================================================================================================
!
!    type :: CV_object
!	    integer						 			    :: N_ave, Nvals
!        integer, parameter                          :: ct_mass=1, ct_momentum=2, ct_energy=3
!	    real(kind(0.d0))						    :: delta_t
!        real(kind(0.d0)), dimension(4)              :: pb, pt, pb_mdt, pt_mdt
!	    real(kind(0.d0)),dimension(:),allocatable	:: dXdt, X, X_mdt, X_oldvol
!	    real(kind(0.d0)),dimension(:),allocatable	:: dsurf_top, dsurf_bot
!	    real(kind(0.d0)),dimension(:,:),allocatable :: flux, stress
!	    contains
!		    procedure :: initialise    => initialise_CV
!		    procedure :: update_X      => update_X_CV
!		    procedure :: update_flux   => update_flux_CV
!		    procedure :: update_stress => update_stress_CV
!		    procedure :: update_dXdt   => update_dXdt_CV
!		    procedure :: update_pbpt   => update_pbpt_CV
!		    procedure :: conservation  => conservation_CV
!    end type CV_object


!    subroutine initialise_CV(self, pt, pb, bintop, binbot, cnsvtype, write_debug=.false.)
!        use computational_constants_MD, only : iter
!        implicit none

!		class(CV_object) :: self

!        integer, intent(in)                           :: cnsvtype
!        logical, intent(in), optional                 :: write_debug
!        double precision, dimension(4), intent(in)    :: pt, pb
!        double precision, dimension(3), intent(in)    :: bintop, binbot

!        integer                                       :: i

!        self%pt = pt; self%pb = pb
!        self%pt_mdt = pt; self%pb_mdt = pb
!        self%bintop = bintop
!        self%binbot = binbot
!        self%cnsvtype = cnsvtype
!        self%write_debug = write_debug

!        !Select type of averaging
!        if (cnsvtype .eq. self%ct_mass) then
!            self%Nvals = 1
!        elseif (cnsvtype .eq. self%ct_momentum) then
!            self%Nvals = 3
!        elseif (cnsvtype .eq. self%ct_energy) then
!	        self%Nvals = 1
!        endif

!        !Select type of averaging
!        if (cnsvtype .eq. ct_mass) then
!            Nvals = 1
!        elseif (cnsvtype .eq. ct_momentum) then
!            Nvals = 3
!        elseif (cnsvtype .eq. ct_energy) then
!	        Nvals = 1
!        endif
!        allocate(self%X_oldvol(self%Nvals))
!        allocate(self%X_mdt(self%Nvals))
!        allocate(self%dX_dt(self%Nvals))
!        allocate(self%dsurf_top(self%Nvals))
!        allocate(self%dsurf_bot(self%Nvals))
!        allocate(self%X_stress(self%Nvals,6))
!        allocate(self%X_flux(self%Nvals,6))
!        allocate(self%X(self%Nvals))

!    end subroutine initialise_CV

!    subroutine update_X_CV()
!        implicit none

!        ! Get all molecules in single CV for whole cluster
!        self%X_mdt = self%X
!        call get_clustCV_out(self%pt, self%pb, self%X, self%nvals, self%cnsvtype, self%write_debug)

!    end subroutine update_X_CV


!    subroutine update_ptpb_CV(pt, pb)
!        implicit none

!        self%pb_mdt = self%pb
!        self%pt_mdt = self%pt
!        self%pb = pb; self%pt = pt
!        
!        ! Assume molecules are fixed and surface is evolving to work out which molecules 
!        ! are now outside/inside per surface. This is done here by defining a 
!        ! CV function based on the volume described as the surface moves from 
!        ! old position to its new position (Negative values are possible here
!        ! and count as molecules which have left in the book keeping).
!        call get_clustCV_out(self%pt, self%pt_mdt, self%dsurf_top, self%nvals, self%cnsvtype, self%write_debug)
!        call get_clustCV_out(self%pb, self%pb_mdt, self%dsurf_bot, self%nvals, self%cnsvtype, self%write_debug)

!    end subroutine update_ptpb_CV

!    subroutine update_flux()
!        implicit none

!        ! Get all surface crossings over all CV surfaces
!        call get_all_surface_crossings(self%pt_mdt, self%pb_mdt,   self%X_cross, & 
!                                       self%nvals,  self%cnsvtype, self%write_debug)

!    end subroutine update_flux

!    subroutine conservation_CV()
!        implicit none

!        ! CV time evolution is equal to molecules crossing surface and molecule change
!        ! due to movement of surface
!        if (abs(sum(self%dX_dt - sum(self%X_cross,2) - (self%dsurf_top-self%dsurf_bot))) .gt. 1e-8) then
!            select case (self%cnsvtype)
!            case(self%ct_mass)
!                print*, "ERROR IN CV FUNCTION FOR MASS CLUSTER"
!            case(self%ct_momentum)
!                print*, "ERROR IN CV FUNCTION FOR MOMENTUM CLUSTER"
!            case(self%ct_energy)
!                print*, "ERROR IN CV FUNCTION FOR ENERGY CLUSTER"
!            case default
!                stop "ERROR in CLUSTER CV, UNKNOWN cnsvtype TYPE"
!            end select
!            !call get_clustCV_out(pt_mdt, pb_mdt, X, nvals, cnsvtype, write_debug=.true.)
!            !call get_all_surface_crossings(pt_mdt, pb_mdt, X_cross, nvals, cnsvtype, write_debug=.true.)
!            do i =1,self%nvals
!                print'(a,2i8,8f13.5)', 'cluster_CV =', iter, i, self%X(i), self%X_cross(i,:), & 
!                                                       self%X(i)-self%X_mdt(i)
!                print'(a,2i8,4f10.5)', 'dsurf', iter, i, self%X_oldvol(i)-self%X(i), & 
!                                    self%dsurf_top(i)-self%dsurf_bot(i), & 
!                                    self%dsurf_top(i), self%dsurf_bot(i)
!            enddo
!            !Calculate new number of molecules in CV for new surfaces
!            !X_oldvol = X
!            !call get_clustCV_out(pt, pb, X, nvals, cnsvtype, write_debug=.false.)
!        else
!            do i =1,self%nvals
!                write(586482,'(i12,i4,10f18.8)'), iter,i, self%X(i), self%X(i)-self%X_mdt(i), & 
!                                                self%X_cross(i,:), -self%dsurf_top(i), self%dsurf_bot(i)
!            enddo
!        endif

!    end subroutine conservation_CV

!========================================================================================================

    ! Calculates change of mass in a control volume and flux over surfaces, then checks the balance is exact.
    ! This is calculated in seperate parts instead of working out the time  
    ! of crossing of each molecule on a moving surface (which is mathematically required/correct).
    ! In practice, determining the moving surface/molecule crossing time is complex 
    ! and would require iteration to find the time of crossing over up to 3 changing surfaces (cubic sides, 
    ! flat faces and the varying intersection between these two).
    ! Instead, the assumption here is:
    ! 1) The surface is fixed and molecule move. 
    ! 2) The molecules are fixed and the surface moves.
    ! Time change in mass is then the sum of both of these seperate contributions
    ! This is assumed to be reasonable as it is not explicitly determined when the surface changes
    ! or if this time evolution is physically meaningful.

    subroutine cluster_CV_fn(pt, pb, bintop, binbot, cnsvtype)
        use computational_constants_MD, only : iter, irank, delta_t
        use librarymod, only : get_new_fileunit, get_Timestep_FileName
        implicit none

        integer, intent(in)                            :: cnsvtype
        double precision, dimension(3), intent(in)     :: bintop, binbot
        double precision, dimension(4), intent(in)     :: pt, pb

        integer                                        :: i, nvals, pid
        integer, parameter                             :: ct_mass=1, ct_momentum=2, ct_energy=3
        character(33)                                  :: filename, debug_outfile
        double precision                               :: conserved
        double precision,dimension(:),allocatable      :: X_mdt, dX_dt
        double precision,dimension(:),allocatable      :: X_oldvol, dsurf_top, dsurf_bot 
        double precision,dimension(:,:),allocatable    :: X_cross, X_stress
        double precision,dimension(:),allocatable,save :: X, CVcount
        logical                                        :: first_time=.true.
        double precision, dimension(4), save           :: pt_mdt, pb_mdt

        !Select type of averaging
        select case (cnsvtype)
        case(ct_mass)
            Nvals = 1
        case(ct_momentum)
            Nvals = 3
        case(ct_energy)
            Nvals = 1
        case default
            stop "ERROR in CLUSTER CV, UNKNOWN cnsvtype TYPE"
        end select

        allocate(X_oldvol(Nvals))
        allocate(X_mdt(Nvals))
        allocate(dX_dt(Nvals))
        allocate(dsurf_top(Nvals))
        allocate(dsurf_bot(Nvals))
        allocate(X_stress(Nvals,6)); X_stress = 0.d0
        allocate(X_cross(Nvals,6))

        if (first_time .eqv. .true.) then
            pt_mdt = pt; pb_mdt = pb
            allocate(X(Nvals))
            first_time = .false.
        endif

        ! Get all molecules in single CV for whole cluster
        X_mdt = X
        call get_clustCV_out(pt_mdt, pb_mdt, bintop, binbot, &
                             X, nvals, cnsvtype, write_debug=.false.)

        ! Get all surface crossings over all CV surfaces
        call get_all_surface_crossings(pt_mdt, pb_mdt, bintop, binbot, &
                                       X_cross, nvals, cnsvtype, write_debug=.false.)

        ! Get all stresses acting over CV surface
        if (cnsvtype .eq. ct_momentum) then
            call get_all_surface_stresses_cells(pt_mdt, pb_mdt, bintop, binbot, &
                                                X_stress, 3, cnsvtype, write_debug=.false.)
            X_stress = X_stress * delta_t
        endif

        !Calculate new number of molecules in CV for new surfaces
        X_oldvol = X
        call get_clustCV_out(pt, pb, bintop, binbot, & 
                             X, nvals, cnsvtype, write_debug=.false.)

        ! Assume molecules are fixed and surface is evolving to work out which molecules 
        ! are now outside/inside per surface. This is done here by defining a 
        ! CV function based on the volume described as the surface moves from 
        ! old position to its new position (Negative values are possible here
        ! and count as molecules which have left in the book keeping).
        call get_clustCV_out(pt, pt_mdt, bintop, binbot, &
                             dsurf_top, nvals, cnsvtype, write_debug=.false.)
        call get_clustCV_out(pb, pb_mdt, bintop, binbot, &
                             dsurf_bot, nvals, cnsvtype, write_debug=.false.)

        !Check for error and write out detail
        dX_dt = X-X_mdt

        ! CV time evolution is equal to molecules crossing surface and molecule change
        ! due to movement of surface
        conserved = abs(sum(dX_dt - sum(X_cross,2) - sum(X_stress,2) - (dsurf_top-dsurf_bot)))
        if (conserved .gt. 1e-8) then
            select case (cnsvtype)
            case(ct_mass)
                print*, "ERROR IN CV FUNCTION FOR MASS CLUSTER"
                do i =1,nvals
				    print'(a,i8,2i4,7f11.5)','Error_clustCV_mass', iter,irank, i, & 
					    conserved, 0.d0, sum(X_cross(i,:)), dX_dt(i), 0.d0, X_mdt(i), X(i)                
                    print'(a,2i8,4f10.5)', 'dsurf', iter, i, X_oldvol(i)-X(i), & 
                                        dsurf_top(i)-dsurf_bot(i), dsurf_top(i), dsurf_bot(i)
                enddo
            case(ct_momentum)
                print*, "ERROR IN CV FUNCTION FOR MOMENTUM CLUSTER"
                do i =1,nvals
				    print'(a,i8,2i4,8f11.5)','Error_clustCV_mom', iter,irank, i, & 
					    conserved, sum(X_stress(i,:)), sum(X_cross(i,:)), dX_dt(i), 0.d0, X_mdt(i), X(i)
				    print'(a,3i8,6f11.5)','clustCV_stresses',iter,irank, i,X_stress(i,:) 
                    print'(a,2i8,4f10.5)', 'dsurf', iter, i, X_oldvol(i)-X(i), &
                                        dsurf_top(i)-dsurf_bot(i), dsurf_top(i), dsurf_bot(i)
                enddo
            case(ct_energy)
                print*, "ERROR IN CV FUNCTION FOR ENERGY CLUSTER"
            case default
                stop "ERROR in CLUSTER CV, UNKNOWN cnsvtype TYPE"
            end select
            !call get_clustCV_out(pt_mdt, pb_mdt, X, nvals, cnsvtype, write_debug=.true.)
            call get_all_surface_crossings(pt_mdt, pb_mdt, bintop, binbot, &
                                           X_cross, nvals, cnsvtype, write_debug=.true.)

            call get_all_surface_stresses_cells(pt_mdt, pb_mdt, bintop, binbot, &
                                                X_stress, 3, cnsvtype, write_debug=.true.)

            !Calculate new number of molecules in CV for new surfaces
            !X_oldvol = X
            !call get_clustCV_out(pt, pb, X, nvals, cnsvtype, write_debug=.false.)
        else
            do i =1,nvals
				!print'(a,i8,i4,7f11.5)','clustCV_mass', iter,irank, & 
				!	conserved, 0.d0, sum(X_cross,2), dX_dt, 0.d0, X_mdt, X
                write(586482,'(i12,i4,11f18.8)') iter,i, conserved, X(i), dX_dt(i), X_cross(i,:), -dsurf_top(i), dsurf_bot(i)
            enddo
        endif

        !===================================
        !TEMP outputting routine
        !===================================

        debug_outfile = './results/clust_CV_dXdt'
        pid = get_new_fileunit()
        call get_Timestep_FileName(iter,debug_outfile,filename)
        open(unit=pid,file=trim(filename),status='replace')
        write(pid,'(i12,6f18.8)') iter, X,X_mdt
        close(pid,status='keep')

        debug_outfile = './results/clust_CV_stress'
        pid = get_new_fileunit()
        call get_Timestep_FileName(iter,debug_outfile,filename)
        open(unit=pid,file=trim(filename),status='replace')
        write(pid,'(i12,18f18.8)') iter, X_stress
        close(pid,status='keep')

        debug_outfile = './results/clust_CV_flux'
        pid = get_new_fileunit()
        call get_Timestep_FileName(iter,debug_outfile,filename)
        open(unit=pid,file=trim(filename),status='replace')
        write(pid,'(i12,18f18.8)') iter, X_cross
        close(pid,status='keep')

        debug_outfile = './results/clust_CV_surf'
        pid = get_new_fileunit()
        call get_Timestep_FileName(iter,debug_outfile,filename)
        open(unit=pid,file=trim(filename),status='replace')
        write(pid,'(i12,6f18.8)') iter, dsurf_top, dsurf_bot
        close(pid,status='keep')

        !===================================
        !TEMP outputting routine
        !===================================

        !Save Control Volume's previous surfaces for next timestep
        pt_mdt = pt; pb_mdt = pb

    end subroutine cluster_CV_fn

! To get the volume of the CV we need integration under curve
!    function int_surface_fn(p0, y)
!        implicit none

!        double precision :: y, int_surface_fn
!        double precision,dimension(4) :: p0

!        int_surface_fn =  (1.d0/4.d0)*p0(4)*y**4 + (1.d0/3.d0)*p0(3)*y**3 &
!                        + (1.d0/2.d0)*p0(2)*y**2 + p0(1)*y

!    end function int_surface_fn

!    function surface_fn_volume(pt, pb, yt, yb)
!        implicit none

!        double precision :: yt, yb
!        double precision,dimension(4) :: pt, pb

!        volume =  int_surface_fn(pt,yt)-int_surface_fn(pt,yb) &
!                 +int_surface_fn(pb,yt)-int_surface_fn(pb,yb) 

!    end function surface_fn_volume

    function surface_fn(p0, yi)
        implicit none

        double precision :: yi, surface_fn
        double precision,dimension(4) :: p0

        surface_fn = p0(4)*yi**3 + p0(3)*yi**2 + p0(2)*yi + p0(1)

    end function surface_fn

    function dsurface_fndyi(p0, yi)
        implicit none

        double precision :: yi, dsurface_fndyi
        double precision,dimension(4) :: p0

        dsurface_fndyi = 3.d0*p0(4)*yi**2 + 2.d0*p0(3)*yi + p0(2)

    end function dsurface_fndyi


    subroutine get_clustCV_out(pt, pb, bintop, binbot, &
                               clustCV_out, nvals, &
                               cnsvtype, write_debug)
        use physical_constants_MD, only : np, halo_np
        use computational_constants_MD, only : iter, delta_t
        use librarymod, only : get_Timestep_FileName, get_new_fileunit
        use arrays_MD, only : r, v, a
        use module_set_parameters, only : mass, potenergymol
        implicit none

        logical, intent(in)                        :: write_debug
        integer, intent(in)                        :: nvals, cnsvtype
        double precision, dimension(3)             :: bintop, binbot
        double precision, dimension(4), intent(in) :: pt, pb
        double precision, dimension(nvals), intent(out)  :: clustCV_out

        integer                            :: n, pid
        integer, parameter                 :: ct_mass=1, ct_momentum=2, ct_energy=3
        character(33)                      :: filename, debug_outfile
        double precision                   :: energy, theta_i
        double precision, dimension(3)     :: ri, velvect
        double precision, dimension(nvals) :: clustCV, qnty

        !Save previous timestep clustCV and reset
        clustCV = 0.d0

        !Plot all molecules inside the liquid cluster control volume
        if (write_debug) then
            debug_outfile = './results/CV_mols'
            pid = get_new_fileunit()
            call get_Timestep_FileName(iter,debug_outfile,filename)
            open(unit=pid,file=trim(filename),status='replace')
        endif

        !Loop over all molecules and halos
        do n =1,np+halo_np

            ri(:) = r(:,n)
            if (cnsvtype .eq. ct_mass) then
                qnty = mass(n)
            elseif (cnsvtype .eq. ct_momentum) then
                qnty = v(:,n)
            elseif (cnsvtype .eq. ct_energy) then
		        velvect(:) = v(:,n) + 0.5d0*a(:,n)*delta_t
		        energy = 0.5d0 * ( mass(n)*dot_product(velvect,velvect) + potenergymol(n))
                qnty = energy
            endif

            call CV_cluster(pt, pb, bintop, binbot, ri, theta_i)
            clustCV = clustCV + qnty * theta_i

            if (write_debug) then
                if (theta_i .eq. 1.d0) then
                    write(pid,'(i10,6f18.9)') n, ri, 0.d0, 0.d0, 0.d0
                else
                    write(pid,'(i10,6f18.9)') n, 0.d0, 0.d0, 0.d0, ri
                endif
            endif

        enddo
        if (write_debug) then
            close(pid,status='keep')
        endif

        clustCV_out = clustCV

    end subroutine get_clustCV_out

    ! A control volume with two intrinsic surfaces in the x directions
    ! and flat surfaces in the y and z directions
    subroutine CV_cluster(pt, pb, bintop, binbot, ri, theta_i)
#if ASSMBLY_HEAVISIDES
		use librarymod, only : heaviside  =>  heaviside_a1
#else
		use librarymod, only : heaviside
#endif
        implicit none

        double precision, dimension(3), intent(in)  :: ri
        double precision, dimension(4), intent(in)  :: pt, pb
        double precision, dimension(3), intent(in)  :: bintop, binbot
        double precision, intent(out)               :: theta_i

        double precision, dimension(3)              :: top, bot        

        !Left/Right cluster based surfaces in x {xsurf = f(yi)}
        top = bintop; bot = binbot
        top(1) = surface_fn(pt, ri(2))
        bot(1) = surface_fn(pb, ri(2))

        !Use CV function
        theta_i = dble(( heaviside(top(1)-ri(1))   & 
                        -heaviside(bot(1)-ri(1)))* & 
	               	   ( heaviside(top(2)-ri(2))   &
                        -heaviside(bot(2)-ri(2)))* & 
	              	   ( heaviside(top(3)-ri(3))   & 
                        -heaviside(bot(3)-ri(3))))

    end subroutine CV_cluster

    ! Add up surface crossings for cubic surfaces and flat surfaces
    subroutine get_all_surface_crossings(pt, pb, bintop, binbot, & 
                                         surfacecross_out, nvals, & 
                                         cnsvtype, write_debug)
        use computational_constants_MD, only : iter, delta_t
        use physical_constants_MD, only : np, halo_np
        use arrays_MD, only : r, v, a
        use module_set_parameters, only : mass, potenergymol
        use librarymod, only : get_Timestep_FileName, get_new_fileunit
        implicit none

        logical, intent(in)                         :: write_debug
        integer, intent(in)                         :: nvals, cnsvtype
        double precision, dimension(3),intent(in)   :: bintop, binbot
        double precision, dimension(4), intent(in)  :: pt, pb
        double precision, dimension(nvals,6), intent(out) :: surfacecross_out

        integer                             :: n, pid
        integer, parameter                  :: ct_mass=1, ct_momentum=2, ct_energy=3
        character(33)                       :: filename, debug_outfile
        double precision                    :: energy
        double precision,dimension(nvals)   :: Ncross, qnty
        double precision,dimension(3)       :: ri1, ri2, vi, velvect
        double precision,dimension(nvals,4) :: Ncross4
        double precision,dimension(nvals,6) :: sc, sc_

        if (write_debug) then
            debug_outfile = './results/CV_surface_mols'
            pid = get_new_fileunit()
            call get_Timestep_FileName(iter,debug_outfile,filename)
            open(unit=pid,file=trim(filename),status='replace')
        endif

        !Loop over all molecules and halos
        sc = 0.d0  
        do n =1,np+halo_np
            !Get position and velocity for molecule
            ri1 = r(:,n);  vi = v(:,n); 
    	    ri2(:) = ri1(:) - delta_t*vi(:)	!Molecule i at time t-dt
            if (cnsvtype .eq. ct_mass) then
                qnty = mass(n)
            elseif (cnsvtype .eq. ct_momentum) then
                qnty = v(:,n)
            elseif (cnsvtype .eq. ct_energy) then
		        velvect(:) = v(:,n) + 0.5d0*a(:,n)*delta_t
		        energy = 0.5d0 * ( mass(n)*dot_product(velvect,velvect) + potenergymol(n))
                qnty = energy
            endif
            ! A function with uses position 1, 2 and qnty to update the array sc which
            ! contains all surface crossings
            call get_all_surface_crossings_r12(pt, pb, bintop, binbot, &
                                               ri1, ri2, qnty, nvals, & 
                                               sc_, write_debug)
            sc = sc + sc_
    
        enddo

        if (write_debug) then
            if (sum(sc(:,2)+sc(:,3)+sc(:,5)+sc(:,6)) .ne. 0) then
                print*, 'Total crossing in y and z is not zero = ',  sc
            endif
            close(pid,status='keep')
        endif

        surfacecross_out = sc

    end subroutine get_all_surface_crossings

    ! Add up stresses of all crossings for cubic surfaces and flat surfaces

    subroutine get_all_surface_stresses(pt, pb, bintop, binbot, & 
                                        surfacecross_out, nvals, & 
                                        cnsvtype, write_debug)
        use computational_constants_MD, only : iter, delta_t, pass_vhalo
        use physical_constants_MD, only : np, halo_np, rcutoff2
        use arrays_MD, only : r, v, a
        use module_set_parameters, only : mass, potenergymol, get_accijmag
        use librarymod, only : get_Timestep_FileName, get_new_fileunit
        use linked_list, only : neighbour, node
        implicit none

        logical, intent(in)                         :: write_debug
        integer, intent(in)                         :: nvals, cnsvtype
        double precision, dimension(3),intent(in)   :: bintop, binbot
        double precision, dimension(4), intent(in)  :: pt, pb
        double precision, dimension(nvals,6), intent(out) :: surfacecross_out

        integer                             :: n, pid
        integer							    :: molnoi, molnoj
        integer							    :: noneighbrs
        type(node), pointer		            :: old, current
        character(33)                       :: filename, debug_outfile
        double precision                    :: rij2, invrij2, accijmag 
        double precision,dimension(3)       :: ri, rj, rij, fij
        double precision,dimension(nvals,6) :: sc_, sc

        !Loop over all molecules and halos
        if (pass_vhalo .eq. 0) stop "ERROR in get_all_surface_stresses -- pass_vhalo should be on for CV cluster"
        sc = 0.d0  
        do molnoi =1,np

            noneighbrs = neighbour%Nlist(molnoi)	!Determine number of elements in neighbourlist
		    old => neighbour%head(molnoi)%point			!Set old to head of neighbour list
		    ri(:) = r(:,molnoi) - delta_t*v(:,molnoi)	!Retrieve ri

            do n = 1,noneighbrs							!Step through all pairs of neighbours molnoi and j

			    molnoj = old%molno						!Number of molecule j
    		    rj(:) = r(:,molnoj) - delta_t*v(:,molnoj)	!Retrieve rj (note VHALOS MUST BE PASSED!!!)
    			!rj(:) = r(:,molnoj)						!Retrieve rj
			    rij(:) = ri(:) - rj(:)   				!Evaluate distance between particle i and j
			    rij2 = dot_product(rij,rij)				!Square of vector calculated

			    if (rij2 < rcutoff2) then
				    invrij2  = 1.d0/rij2                !Invert value
                    accijmag = get_accijmag(invrij2, molnoi, molnoj)

                    ! A function with uses position 1, 2 and qnty to update the array sc which
                    ! contains all surface crossings
					if (molnoj .gt. np) then
    				    fij(:) = accijmag*rij(:)
                        PRINT*, "THIS USED TO INCLUDE A HALF IN THE NON HALO VALUES BUT THAT WAS ALWAYS WRONG."
                        PRINT*, "THE DXDT = STRESS OVER CUBIC SURFACES BUT THE OTHER SURFACES DO NOT BALANCE"
                        PRINT*, "THIS APPEARS TO BE AN ISSUE WITH MOLECULAR INTERACTIONS WHICH CROSS FROM HALOS"
                        PRINT*, "PERHAPS WE NEED TO LOOP THROUGH HALO MOLECULES TO GET INTERACTIONS WITH HERE"
                        stop "Error in get_all_surface_stresses -- Change this routine for cell list with all halo molecules"
                    else
						fij(:) = accijmag*rij(:)
                    endif
                    call get_all_surface_crossings_r12(pt, pb, bintop, binbot, &
                                                       ri, rj, fij, nvals, sc_, write_debug)

                    sc = sc + sc_
					!if ((molnoj .gt. np) .and. (any(abs(sc_) .gt. 1e-8))) then
                    !    print'(l, 2i8,18f10.5)', any(abs(sc) .gt. 1e-8), molnoi, molnoj, sc
                    !endif

                endif
			    current => old
			    old => current%next !Use pointer in datatype to obtain next item in list
            enddo
        enddo

    	nullify(old)      		!Nullify as no longer required
    	nullify(current)      	!Nullify as no longer required

!        if (write_debug) then
!            if (sum(sc(:,2)+sc(:,3)+sc(:,5)+sc(:,6)) .ne. 0) then
!                print*, 'Total crossing in y and z is not zero = ',  sc
!            endif
!            close(pid,status='keep')
!        endif

        surfacecross_out = sc

    end subroutine get_all_surface_stresses




    subroutine get_all_surface_stresses_cells(pt, pb, bintop, binbot, & 
                                              surfacecross_out, nvals, & 
                                              cnsvtype, write_debug)
        use computational_constants_MD, only : iter, delta_t, pass_vhalo, & 
                                               binspercell, ncells
        use physical_constants_MD, only : np, halo_np, rcutoff2
        use arrays_MD, only : r, v, a
        use module_set_parameters, only : get_accijmag
        use librarymod, only : get_Timestep_FileName, get_new_fileunit
        use linked_list, only : cell
        use linked_list, only : node
        implicit none

        logical, intent(in)                         :: write_debug
        integer, intent(in)                         :: nvals, cnsvtype
        double precision, dimension(3),intent(in)   :: bintop, binbot
        double precision, dimension(4), intent(in)  :: pt, pb
        double precision, dimension(nvals,6), intent(out) :: surfacecross_out

        integer                             :: n, pid
        integer							    :: molnoi, molnoj
	    integer                             :: i, j, ixyz !Define dummy index
	    integer							    :: icell, jcell, kcell
	    integer                             :: icellshift, jcellshift, kcellshift
	    integer                             :: cellnp, adjacentcellnp 
	    integer 							:: icellmin,jcellmin,kcellmin
	    integer 							:: icellmax,jcellmax,kcellmax
        integer, parameter                  :: ct_mass=1, ct_momentum=2, ct_energy=3
	    type(node), pointer 	            :: oldi, currenti, oldj, currentj
        character(33)                       :: filename, debug_outfile
        double precision                    :: rij2, invrij2, accijmag
        double precision,dimension(nvals)   :: qnty
        double precision,dimension(3)       :: ri, rj, rij, fij
        double precision,dimension(nvals,6) :: sc_, sc

	    real(kind(0.d0)),dimension(3)	:: vi_t, cellsperbin
	    real(kind(0.d0)), dimension(:,:), allocatable :: rF
	    real(kind(0.d0)), dimension(3,1):: rFv

        if (pass_vhalo .eq. 0) stop "ERROR in get_all_surface_stresses -- pass_vhalo should be on for CV cluster"
        sc = 0.d0  

	    !Calculate bin to cell ratio
	    cellsperbin = 1.d0/binspercell !ceiling(ncells(1)/dble(nbins(1)))

        ! Still need to loop over every cell (i.e. get all interactions) if
        ! bins are bigger than cells
	    where (cellsperbin .ge. 1.d0) cellsperbin = 1.d0

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
	   		    ri(:) = r(:,molnoi) - delta_t*v(:,molnoi)	!Retrieve ri
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
            		    rj(:) = r(:,molnoj) - delta_t*v(:,molnoj)	!Retrieve rj (note VHALOS MUST BE PASSED!!!)

					    currentj => oldj
					    oldj => currentj%next    !Use pointer in datatype to obtain next item in list

					    if(molnoi==molnoj) cycle !Check to prevent interaction with self

					    rij2=0                   !Set rij^2 to zero
					    rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j
					    rij2 = dot_product(rij,rij)	!Square of vector calculated

					    if (rij2 < rcutoff2) then
				            invrij2  = 1.d0/rij2                !Invert value
                            accijmag = get_accijmag(invrij2, molnoi, molnoj)

                            ! A function with uses position 1, 2 and qnty to update the array sc which
                            ! contains all surface crossings
            			    fij(:) = 0.5d0*accijmag*rij(:)
                            if (cnsvtype .eq. ct_momentum) then
                                qnty(:) = fij(:)
                            elseif (cnsvtype .eq. ct_energy) then
                                vi_t = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)
                                qnty = dot_product(fij, vi_t)
                            endif
                            call get_all_surface_crossings_r12(pt, pb, bintop, binbot, &
                                                               ri, rj, qnty, nvals, sc_, write_debug)

                            sc = sc + sc_
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

        surfacecross_out = sc

	    nullify(oldi)      	!Nullify as no longer required
	    nullify(oldj)      	!Nullify as no longer required
	    nullify(currenti)      	!Nullify as no longer required
	    nullify(currentj)      	!Nullify as no longer required

    end subroutine get_all_surface_stresses_cells

    subroutine get_all_surface_crossings_r12(pt, pb, bintop, binbot, & 
                                             ri1, ri2, qnty, nvals, & 
                                             sc, write_debug)
        implicit none

        logical, intent(in)                        :: write_debug
        integer, intent(in)                        :: nvals
        double precision, dimension(3),intent(in)  :: ri1, ri2
        double precision, dimension(3),intent(in)  :: bintop, binbot
        double precision, dimension(4),intent(in)  :: pt, pb
        double precision, dimension(nvals),intent(in)  :: qnty
        double precision,dimension(nvals,6),intent(out) :: sc

        double precision, dimension(nvals)   :: Ncross
        double precision, dimension(nvals,4) :: Ncross4

        !Top cubic surface crossing
        call get_cubic_surface_crossing(pt, bintop, binbot, &
                                        ri1, ri2, nvals, qnty, &
                                        Ncross, write_debug)
        sc(:,1) = + Ncross(:)
        !Bottom cubic surface crossing
        call get_cubic_surface_crossing(pb, bintop, binbot, &
                                        ri1, ri2, nvals, qnty, &
                                        Ncross, write_debug)
        sc(:,4) = - Ncross(:)
        !periodic boundries and tethered molecule surface crossing
        !Note sign convention seems wrong here, opposite of curved surfaces...
        call get_plane_surface_crossing(pt, pb, bintop, binbot, & 
                                        ri1, ri2, nvals, qnty, &
                                        Ncross4, write_debug)
        sc(:,2) = - Ncross4(:,1)
        sc(:,3) = - Ncross4(:,2)
        sc(:,5) = + Ncross4(:,3)
        sc(:,6) = + Ncross4(:,4)

    end subroutine get_all_surface_crossings_r12

    ! Crossings over cubic surface

    subroutine get_cubic_surface_crossing(p0, bintop, binbot, & 
                                          ri1, ri2, N, qnty, &
                                          Ncross, write_debug)
        use computational_constants_MD, only : iter
#if ASSMBLY_HEAVISIDES
    use librarymod, only : heaviside  =>  heaviside_a1
#else
    use librarymod, only : heaviside
#endif
        use PolynomialRoots, only : QuadraticRoots, LinearRoot, cubicroots, SolvePolynomial
        implicit none

        integer, intent(in)                        :: N
        logical, intent(in)                        :: write_debug
        double precision, dimension(N),intent(in)  :: qnty
        double precision, dimension(3),intent(in)  :: ri1, ri2
        double precision, dimension(3),intent(in)  :: bintop, binbot
        double precision, dimension(4),intent(in)  :: p0
        double precision, dimension(N),intent(out) :: Ncross

        integer                           :: i, pid, code, crossings
        double precision, parameter       :: tol=1e-14
        double precision                  :: dS_i, tcross, m,c,crosssign, dsdy
        double precision, dimension(4)    :: p0l
        double precision, dimension(3)    :: vi, ri12, rcross
        complex(KIND(1.0D0))              :: temp
        complex(KIND(1.0D0)),dimension(3) :: z

        !reset surface crossing
        Ncross = 0

        !Get intersection of moving molecule and current surface
	    ri12   = ri1 - ri2		        !Molecule i trajectory between t-dt and t
        if (abs(ri12(2)) .gt. tol) then
            m = ri12(1)/ri12(2)
        else
            m = ri12(1)/tol
        endif
        c = ri1(1)-ri1(2)*m

        ! Solution of cubic and line:
        ! p0(4)*x**3 + p0(3)*x**2 + p0(2)*x + p0(1) = (mx + c) 
        ! p0(4)*x**3 + p0(3)*x**2 + (p0(2)-m)*x + p0(1)-c = 0
        ! then find roots to determine points of crossing
        p0l(:) = p0(:)
        p0l(1) = p0(1) - c
        p0l(2) = p0(2) - m

        code = 0
!        call SolvePolynomial(0.d0, p0l(4), p0l(3), p0l(2),p0l(1), &
!                             code, z(1),z(2),z(3),temp)
        if (abs(p0l(2)) .lt. tol .and. abs(p0l(3)) .lt. tol .and.  & 
            abs(p0l(4)) .lt. tol) then
            !flat (and line is exactly parallel?)
            z(1) = cmplx(0.d0, 1.d0, kind(1.d0))
            z(2) = cmplx(0.d0, 1.d0, kind(1.d0))
            z(3) = cmplx(0.d0, 1.d0, kind(1.d0))
        elseif (abs(p0l(3)) .lt. tol .and. abs(p0l(4)) .lt. tol) then
            !Linear
            z(1) = cmplx(-p0l(1)/p0l(2), 0.d0, kind(1.d0))
            z(2) = cmplx(0.d0, 1.d0, kind(1.d0))
            z(3) = cmplx(0.d0, 1.d0, kind(1.d0))
        elseif (abs(p0l(4)) .lt. tol) then
            !Quadratic
            call QuadraticRoots(p0l(1:3), z)
            z(3) = cmplx(0.d0, 1.d0, kind(1.d0))
        else
            !Cubic
            call CubicRoots(p0l, z)
        endif

        !Check if any roots are real
        crossings = 0
        do i = 1, size(z)
            if (abs(imag(z(i))) .lt. tol) then

                !Get time of crossing
                if (abs(ri12(2)) .gt. tol) then
                    tcross = (ri1(2) - dble(z(i))) / ri12(2)
                else
                    tcross = (ri1(2) - dble(z(i))) / tol
                endif
                rcross = (/ ri1(1)-ri12(1)*tcross, & 
                            ri1(2)-ri12(2)*tcross, & 
                            ri1(3)-ri12(3)*tcross /)

                !Get surface crossing function
                dS_i = dble((heaviside( tcross )           -heaviside(tcross - 1.d0))* & 
                            (heaviside(bintop(2)-rcross(2))-heaviside(binbot(2)-rcross(2)))* & 
                      	    (heaviside(bintop(3)-rcross(3))-heaviside(binbot(3)-rcross(3))))

!                if (iter .eq. 62 .and. i .ge. 2 .and. abs(real(z(i))) .lt. 1.d0) then
!                    print'(i4, 13f10.5)', iter, z(i), ri1, ri2, dS_i, tcross, rcross
!                endif

                if (abs(dS_i) .gt. tol) then

                    !Get surface crossing direction
                    dsdy = dsurface_fndyi(p0, rcross(2))
                    crosssign = sign(1.d0,(-ri12(1) + ri12(2)*dsdy))
                    !Get normal and tangent of surface
                    !norm => x - surface_fn(rcross(2)) = -(1.d0/dsdy)*(y - rcross(2))
                    !tang => x - surface_fn(rcross(2)) = dsdy * (y - rcross(2))
                    ! and project qnty (if mom vector only) along normal and tangent 
                    ! using qnty = (/ dot_product(qnty, (/norm, tang, z/), & 
                    !                 dot_product(qnty, (/tang, norm, z/))
                    Ncross(:) = Ncross(:) + qnty(:) * crosssign * dS_i
                    if (write_debug) then
                        crossings = crossings + 1
                        print('(a,3i6,2f6.2,10f10.5)'), 'xing', iter, i, crossings, & 
                                                         dS_i, crosssign, tcross, & 
                                                         ri1(:), ri2(:), rcross(:)
                        !It seems to be possible to get real roots which 
                        !don't actually give a zero value in the original function!?
                        if (surface_fn(p0, rcross(2)) .gt. 1e-8) then
                            print'(a,f15.10,7f12.3)', "root isn't a zero of function", & 
                                                      surface_fn(p0, rcross(2)), z, tcross
                            !stop 'ERROR -- root not zero'
                        endif
                    endif

                    if (crossings .eq. 3) print'(a,2i4,6f10.5)', 'three roots', iter, i,  ri1(:), ri2(:)

                endif

            endif
        enddo

    end subroutine get_cubic_surface_crossing

    !crossings for all other flat surfaces

    subroutine get_plane_surface_crossing(pt, pb, bintop, binbot, & 
                                          ri1, ri2, nvals, qnty, &
                                          Ncross, write_debug)
        use computational_constants_MD, only : iter, delta_t
#if ASSMBLY_HEAVISIDES
		use librarymod, only : heaviside  =>  heaviside_a1
#else
		use librarymod, only : heaviside
#endif
        implicit none

        integer, intent(in)                           :: nvals
        logical, intent(in)                           :: write_debug
        double precision, dimension(nvals), intent(in):: qnty
        double precision, dimension(3),intent(in)     :: ri1, ri2
        double precision, dimension(3),intent(in)     :: bintop, binbot
        double precision, dimension(4),intent(in)     :: pt, pb
        double precision, dimension(nvals,4), intent(out) :: Ncross

        integer                         :: i
        double precision, dimension(3)  :: ri12, rcross, top, bot
        double precision, dimension(3)  :: Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
        real(kind(0.d0))				:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
        complex(kind(0.d0))             :: z
        double precision, parameter     :: tol=1e-8

        !Copy to variable to define polynomial surfaces as bintop/binbot intent in 
        top = bintop; bot = binbot

        !reset surface crossing
        Ncross = 0.d0

        !Get intersection of moving molecule and plane surfaces
	    ri12   = ri1 - ri2		        !Molecule i trajectory between t-dt and t
	    where (abs(ri12) .lt. 0.000001d0) ri12 = 0.000001d0

		Pyt=(/	ri1(1)+(ri12(1)/ri12(2))*(top(2)-ri1(2)), & 
					top(2), 		     & 
				ri1(3)+(ri12(3)/ri12(2))*(top(2)-ri1(2))  	/)
		Pyb=(/	ri1(1)+(ri12(1)/ri12(2))*(bot(2)-ri1(2)), &
					bot(2), 		     & 
				ri1(3)+(ri12(3)/ri12(2))*(bot(2)-ri1(2))  	/)
		Pzt=(/	ri1(1)+(ri12(1)/ri12(3))*(top(3)-ri1(3)), & 
				ri1(2)+(ri12(2)/ri12(3))*(top(3)-ri1(3)), &
					top(3) 			/)
		Pzb=(/	ri1(1)+(ri12(1)/ri12(3))*(bot(3)-ri1(3)), &
				ri1(2)+(ri12(2)/ri12(3))*(bot(3)-ri1(3)), & 
					bot(3) 			/)

        !Get x value of top and bottom surfaces

        !Y TOP SURFACE
        top(1) = surface_fn(pt, Pyt(2))
        bot(1) = surface_fn(pb, Pyt(2))
		onfaceyt =0.5d0*(sign(1.d0,top(2) - ri2(2))   &
				       - sign(1.d0,top(2) - ri1(2)))* &
						(heaviside(top(1) - Pyt(1))   &
				       - heaviside(bot(1) - Pyt(1)))* &
						(heaviside(top(3) - Pyt(3))   &
				       - heaviside(bot(3) - Pyt(3)))

        !Y BOTTOM SURFACE
        top(1) = surface_fn(pt, Pyb(2))
        bot(1) = surface_fn(pb, Pyb(2))
		onfaceyb =0.5d0*(sign(1.d0,bot(2) - ri2(2))   &
				       - sign(1.d0,bot(2) - ri1(2)))* &
						(heaviside(top(1) - Pyb(1))   &
				       - heaviside(bot(1) - Pyb(1)))* &
						(heaviside(top(3) - Pyb(3))   &
				       - heaviside(bot(3) - Pyb(3)))

        !Z TOP SURFACE
        top(1) = surface_fn(pt, Pzt(2))
        bot(1) = surface_fn(pb, Pzt(2))
		onfacezt =0.5d0*(sign(1.d0,top(3) - ri2(3))   &
				       - sign(1.d0,top(3) - ri1(3)))* &
						(heaviside(top(1) - Pzt(1))   &
					   - heaviside(bot(1) - Pzt(1)))* &
						(heaviside(top(2) - Pzt(2))   &
				       - heaviside(bot(2) - Pzt(2)))

        !Z BOTTOM SURFACE
        top(1) = surface_fn(pt, Pzb(2))
        bot(1) = surface_fn(pb, Pzb(2))
		onfacezb =0.5d0*(sign(1.d0,bot(3) - ri2(3))   &
				       - sign(1.d0,bot(3) - ri1(3)))* &
						(heaviside(top(1) - Pzb(1))   &
				       - heaviside(bot(1) - Pzb(1)))* &
						(heaviside(top(2) - Pzb(2))   &
				       - heaviside(bot(2) - Pzb(2)))

        !Ncross(1) = Ncross(1) + int(onfaceyt - onfaceyb + onfacezt - onfacezb)
        Ncross(:,1) = Ncross(:,1) + qnty(:)*onfaceyt
        Ncross(:,2) = Ncross(:,2) + qnty(:)*onfacezt
        Ncross(:,3) = Ncross(:,3) + qnty(:)*onfaceyb
        Ncross(:,4) = Ncross(:,4) + qnty(:)*onfacezb

        if (write_debug) then
!            if (onfaceyt .ne. 0.d0) then
!                print('(a,4f6.2,i4,9f8.3)'), 'yplanet xing',  &
!                    onfaceyt, onfaceyb, onfacezt, onfacezb, &
!                    int(onfaceyt - onfaceyb + onfacezt - onfacezb), &
!                    ri1, ri2, Pyt
!            endif
!            if (onfaceyb .ne. 0.d0) then
!                print('(a,4f6.2,i4,9f8.3)'), 'yplaneb xing',  &
!                    onfaceyt, onfaceyb, onfacezt, onfacezb, &
!                    int(onfaceyt - onfaceyb + onfacezt - onfacezb), &
!                    ri1, ri2, Pyb
!            endif
!            if (onfacezt .ne. 0.d0) then
!                print('(a,4f6.2,i4,9f8.3)'), 'zplanet xing',  &
!                    onfaceyt, onfaceyb, onfacezt, onfacezb, &
!                    int(onfaceyt - onfaceyb + onfacezt - onfacezb), &
!                    ri1, ri2, Pzt
!            endif
!            if (onfacezb .ne. 0.d0) then
!                print('(a,4f6.2,i4,9f8.3)'), 'zplaneb xing',  & 
!                     onfaceyt, onfaceyb, onfacezt, onfacezb, & 
!                    int(onfaceyt - onfaceyb + onfacezt - onfacezb), &
!                    ri1, ri2, Pzb
!                !write(pid,'(2i6,f6.2,7f10.5,3f18.12)'), n, i, dS_i, tcross, ri(:), ri(:)-vi(:)*delta_t, rcross(:)
!            endif

        endif

    end subroutine get_plane_surface_crossing

end module cubic_surface_CV




!=======================================================================
!
!   C L U S T E R   A N A L Y S I S
!
!-----------------------------------------------------------------------

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

        double precision, intent(in)       :: rd

        integer                            :: i, allocsize

        !Initialised cluster list
        allocsize = np+extralloc
        if (.not. allocated(self%Nlist)) then
            allocate(self%Nlist(allocsize))
            allocate(self%head(allocsize))
        endif
        if (.not. allocated(self%inclust)) then
            allocate(self%inclust(allocsize))
        endif
        if (.not. allocated(self%clusterngbrs)) then
            allocate(self%clusterngbrs(allocsize))
        endif

        self%Nclust = 0
	    do i = 1,allocsize
            self%inclust(i) = 0
            !self%clusterngbrs(i) = 0
		    self%Nlist(i) = 0	!Zero number of molecules in cluster list
		    nullify(self%head(i)%point)!Nullify cluster list head pointer 
	    enddo

        !Call to build clusters from neighbour and cell lists
        if (force_list .eq. 1) then
            call build_from_cell_list(self, cell, rd,  & 
                                      r, np, skipwalls_=.true.)
        else if (force_list .eq. 2) then
            call build_from_neighbour_list(self, neighbour, rd,  & 
                                           r, np, skipwalls_=.true.)
        else 
            call error_abort("Error in build_clusters -- full cell or "//&
                             "neightbour list should be used with interface tracking")
        endif

        !Remove all empty cluster references
        call CompressClusters(self)

    end subroutine build_clusters

    subroutine fit_surface(x, y, z, fittype, p0, f, debug_outfile)
        use computational_constants_MD, only : iter 
        use physical_constants_MD, only : pi
        use minpack_fit_funcs_mod, only : fn, cubic_fn, cubic_fn2D, & 
                                          cubic_fn2D_16coeff, curve_fit
        use librarymod, only : least_squares, get_new_fileunit
        implicit none

        integer, intent(in) :: fittype
        character(23), intent(in), optional :: debug_outfile
        double precision,dimension(:), &
            allocatable, intent(in) :: x, y, z
        double precision,dimension(:), &
            allocatable, intent(out) :: p0, f 

        logical :: first_time=.true.
        integer :: fileunit
        integer, parameter :: DEBUG=0, ONEDIM=1, TWODIM=2, TWODIM16=3
        double precision   :: m, c, cl_angle
        !Save previous surface as initial guess
        double precision,dimension(:), &
            allocatable, save :: p0_ 
        if (present(debug_outfile)) then
            fileunit = get_new_fileunit()
            if (first_time .eqv. .true.) then
                open(unit=fileunit,file=debug_outfile,status='replace')
                first_time = .false.
                if ((fittype .eq. DEBUG) .or. (fittype .eq. ONEDIM)) then
                    allocate(p0_(4)); p0_ = 0.d0 
                else if (fittype .eq. TWODIM) then
                   allocate(p0_(8)); p0_ = 0.d0
                else if (fittype .eq. TWODIM16) then
                    allocate(p0_(16)); p0_ = 0.d0
                endif
            else
                open(unit=fileunit,file=debug_outfile,access='append')
            endif
        endif

        !Linear with angle
        call least_squares(x, y, m, c)
        cl_angle = 90.d0+atan(m)*180.d0/pi

        if (fittype .eq. DEBUG) then
            allocate(p0(4)); p0 = p0_
            allocate(f(size(x,1)))
            f = 0.d0
            !Set dummy values of CV surfaces
            if (minval(y) .lt. 0.d0) then
                !p0 = (/-3.d0+cos(2*3.14159*iter/1000), 0.1d0, 0.02d0, -0.001d0  /)
                p0 = (/-3.d0, -0.1d0, -0.01d0, -0.001d0  /)
            else            
                !p0 = (/ 3.d0+sin(2*3.14159*iter/1000), 0.1d0, 0.02d0, -0.001d0  /)
                p0 = (/ 3.d0, 0.1d0,  0.01d0, -0.001d0  /)
            endif
            !p0 = (/ 1.d0, 0.5d0,  0.1d0, -0.80d0  /)
            !p0 = (/-0.2d0, 0.5d0,  0.2d0, -0.00d0  /)
            if (present(debug_outfile)) then
                write(fileunit,'(i12, 7f15.8)') iter, m, c, cl_angle, p0
            endif
        else if (fittype .eq. ONEDIM) then
            allocate(p0(4)); p0 = p0_
            fn => cubic_fn
            call curve_fit(fn, x, y, p0, f)
            if (present(debug_outfile)) then
                write(fileunit,'(i12, 7f15.8)') iter, m, c, cl_angle, p0
            endif
        else if (fittype .eq. TWODIM) then
            allocate(p0(8)); p0 = p0_
            fn => cubic_fn2D
            call curve_fit(fn, x, y, p0, f, z)
            if (present(debug_outfile)) then
                write(fileunit,'(i12, 11f15.8)') iter, m, c, cl_angle, p0
            endif
       else if (fittype .eq. TWODIM16) then
            allocate(p0(16)); p0 = p0_
            fn => cubic_fn2D_16coeff
            call curve_fit(fn, x, y, p0, f, z)
            if (present(debug_outfile)) then
                write(fileunit,'(i12, 19f15.8)') iter, m, c, cl_angle, p0
            endif
        endif
        !Save soln  as initial guess for next time
        p0_ = p0

        if (present(debug_outfile)) close(fileunit,status='keep')

    end subroutine fit_surface


    subroutine get_cluster_properties(self, rd, min_ngbr)
        use physical_constants_MD, only : np, nd, tethereddisttop, tethereddistbottom
        use computational_constants_MD, only : iter, tplot, thermo_tags, thermo, &
                                               free, globaldomain, intrinsic_interface_outflag, &
                                               II_normal, II_alpha, II_tau, II_eps, II_ns, II_topbot, &
                                               mflux_outflag, Nsurfevo_outflag, nhb, & 
                                               CA_generate_xyz, CA_generate_xyz_res, &
											   mflux_outflag, vflux_outflag, CV_conserve
        use librarymod, only : imaxloc, get_Timestep_FileName, least_squares, get_new_fileunit, & 
								write_wave_xyz, write_waveobj
        use minpack_fit_funcs_mod, only : fn, cubic_fn, curve_fit
        use arrays_MD, only : tag, r, intnscshift, glob_no
#if ASSMBLY_HEAVISIDES
    use librarymod, only : heaviside  =>  heaviside_a1
#else
    use librarymod, only : heaviside
#endif
        use intrinsic_module, only : fit_intrinsic_surface_bilinear, fit_intrinsic_surface_modes
        use calculated_properties_MD, only : nbins, nbinso, binsize, mass_surface_flux
        use module_record, only : Abilinear, ISR, ISR_mdt, ISR_r, ISR_mdt_r, ISR_b, ISR_mdt_b
        use interfaces, only : error_abort
        use cubic_surface_CV, only : cluster_CV_fn
        use module_record, only : get_bin, get_bin_molno

        implicit none

        type(clusterinfo),intent(inout)    :: self
        integer, intent(in)                :: min_ngbr
        double precision, intent(in)       :: rd

        logical, save :: write_cluster_header=.true., first_sample=.true.
        logical,save                    :: first_time=.true., first_time_coeff=.true.
        character(32)                   :: filename, debug_outfile
        integer                         :: n,i,j,ixyz, jxyz,resolution,fittype,normal,topbot, clustNo, bins(3)
        integer, dimension(:), allocatable :: pivots
        double precision                :: zeromode, maprange
        double precision                :: tolerance, alpha, tau, eps, ns, area
        double precision, dimension(3)  :: bintop, binbot, box !Used in CV
        double precision,dimension(6)   :: extents
        double precision,dimension(4)   :: ptin, pbin
        double precision,dimension(:),allocatable :: x,y,z,f,pt,pb, elevation
        double precision,dimension(:),allocatable,save :: coeffmdt
        double precision,dimension(:,:),allocatable :: rnp, extents_grid
        double precision,dimension(:,:),allocatable :: molnos, clustmolnos
        double precision,dimension(:,:),allocatable, save :: points
        double precision,dimension(:,:,:),allocatable :: vertices

        double precision,dimension(:,:),allocatable,save :: intrnsc_smplemdt
        double precision,dimension(:,:,:,:),allocatable,save :: Abilinearmdt

        integer, save :: writeiter=0

        clustNo = imaxloc(self%Nlist)

        !Different intrinsic surface types
        if (intrinsic_interface_outflag .eq. 0) then
            !Do nothing, no intrinsic interface
            return
        else if (any(intrinsic_interface_outflag .eq. (/ 1, 2 /))) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Fit intrinsic (sine/cosine) surface !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Intrinsic surface coefficients
            normal = II_normal ! 3
            alpha =  II_alpha ! 0.5d0
            tau =  II_tau ! 1.d0
            eps =  II_eps !0.00000001d0
            ns =  II_ns !0.8d0
            topbot = II_topbot

            if (first_time) then
				if (intrinsic_interface_outflag .eq. 1) then
					ISR => ISR_r
					ISR_mdt => ISR_mdt_r
				elseif (intrinsic_interface_outflag .eq. 2) then
					ISR => ISR_b
					ISR_mdt => ISR_mdt_b
				endif
                call ISR%initialise(globaldomain, normal, alpha, eps, nbins, nhb, topbot)   ! initialise
                call ISR_mdt%initialise(globaldomain, normal, alpha, eps, nbins, nhb, topbot)   ! initialise
                first_time = .false.
			!else
				!do i=1,np
					!print*, i, r(:,i), get_bin(r(:,i)), get_bin_molno(i)
					!if (any(get_bin(r(:,i)) .ne. get_bin_molno(i))) then 
					!	print*, "Error in ISR%get_bin",i,r(:,i), get_bin(r(:,i)), get_bin_molno(i)
					!endif
				!enddo

            endif


            !Only recheck external molecules of cluster every tplot timesteps
            if (CV_conserve .eq. 1 .or. mod(iter,tplot) .eq. 0) then

                !Get cluster data into array
                call cluster_to_array(self, clustNo, r, min_ngbr, rnp)
				
                !x normal is tested appears to work
                if (allocated(points)) deallocate(points)
                allocate(points(size(rnp,2), size(rnp,1)))
                points = 0.d0
                do i =1, size(rnp,2)
                    points(i,1) = rnp(1,i)
                    points(i,2) = rnp(2,i)
                    points(i,3) = rnp(3,i)
                enddo

				!Copy previous object before updating
				!if (.not. first_time_coeff) then
				!ISR%coeff = 0.d0
				ISR_mdt%coeff = ISR%coeff
				if (intrinsic_interface_outflag .eq. 2) then
					ISR_mdt%Abilinear = ISR%Abilinear
					ISR_mdt%intrnsc_smple = ISR%intrnsc_smple
				endif
				!endif

                !Get surface in terms of modes
				call ISR%fit_intrinsic_surface(points, tau, ns, pivots)
 				
				!Write out surface modes to file
				!call ISR_mdt%write_modes(iter)

				!Save initial surface for debugging
				!if (first_time_coeff) then
					! allocate(coeffmdt(size(ISR%coeff,1)))
					! coeffmdt(:) = ISR%coeff(:)
					! if (intrinsic_interface_outflag .eq. 2) then
						! allocate(Abilinearmdt(2,2,size(ISR%Abilinear,3),size(ISR%Abilinear,4)))
						! Abilinearmdt = ISR%Abilinear
						! allocate(intrnsc_smplemdt(size(ISR%intrnsc_smple,1),size(ISR%intrnsc_smple,2)))
						! intrnsc_smplemdt = ISR%intrnsc_smple
					! endif
					! first_time_coeff = .false.
				!else
					!print*, "DEBUG in get_cluster_properties, setting coeff to zero"
					!ISR%coeff = 0.d0
					! print*, "DEBUG in get_cluster_properties, setting coeff to intial"
					! ISR%coeff(:)=coeffmdt(:)
					! !ISR%coeff(313)=ISR%coeff(313)+ 1.0*sin((iter-100000)/100.d0) !shift from sin(0) mode

					! !ISR_mdt%coeff = 0.d0 
					! !ISR%coeff(313)= 2.d0*sin((iter-100000-1)/1000.d0) !shift from sin(0) mode
					!ISR%coeff(312)=0.5d0 !sin (2*pi/Lx)
					! !ISR%coeff(314)=0.5d0
					! if (intrinsic_interface_outflag .eq. 2) then
					!ISR%intrnsc_smple = 0.d0
						! ISR%intrnsc_smple = intrnsc_smplemdt
						! !ISR%intrnsc_smple = ISR%intrnsc_smple + 1.0*sin((iter-100000)/100.d0)
					!ISR%Abilinear = 0.d0
						! ISR%Abilinear = Abilinearmdt
						! !ISR%Abilinear(1,1,:,:) = ISR%Abilinear(1,1,:,:) + 1.0*sin((iter-100000)/100.d0)
					! endif

				!endif

				!Get surface crossings due to surface's evolution
				if (Nsurfevo_outflag .ne. 0) then
					call surface_evolution(ISR, ISR_mdt, .false.)
				endif

				!DEBUG - write surface out
				if (CA_generate_xyz .eq. 1) then
                    if (CA_generate_xyz_res .gt. 0) then
                        print*, "RESOLUTION NOT ZERO", CA_generate_xyz_res
    					call ISR%sample_surface(vertices, nbins=(/1, CA_generate_xyz_res, CA_generate_xyz_res/), &
                                                writeiter=writeiter)
                        !Default size writes bilinear as well
    					call ISR%sample_surface(vertices, writeiter=writeiter)
                        writeiter = writeiter + 1
                    else
    					call ISR%sample_surface(vertices, writeiter=writeiter)
                        writeiter = writeiter + 1
                    endif
					call write_wave_xyz(vertices)
					!Store pivots in intrinsic surface to plot
					if (allocated(ISR%pivots)) deallocate(ISR%pivots)
					allocate(ISR%pivots(size(pivots,1)))
					ISR%pivots = pivots
				elseif (CA_generate_xyz .eq. 2) then
                    if (CA_generate_xyz_res .gt. 0) then
    					call ISR%sample_surface(vertices, (/1, CA_generate_xyz_res, CA_generate_xyz_res/))
                    else
    					call ISR%sample_surface(vertices)
                    endif
					call write_waveobj(vertices, iter)
				endif

				!Get shift for intrinsic surface for each molecule
				deallocate(points)
				allocate(points(np, nd))
				points(:,1) = r(1,1:np)
				points(:,2) = r(2,1:np)
				points(:,3) = r(3,1:np)
				call ISR_mdt%get_surface(points, elevation, include_zeromode=.true.)
				intnscshift(1:np) = elevation(1:np)

            endif

        elseif (intrinsic_interface_outflag .eq. 3) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !     Fit linear and cubic surface    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            resolution = 10; tolerance = rd
            fittype = 3   

            call cluster_global_extents(self, clustNo, extents)

            !Get molecules on top surface
            call cluster_extents_grid(self, clustNo, 1, resolution, & 
                                      extents_grid)!, debug_outfile='./results/maxcell_top')
            call cluster_outer_mols(self, clustNo, tolerance=tolerance, dir=1, & 
                                    rmols=rnp, extents=extents_grid, debug_outfile='./results/clust_edge_top')

            !Curve fits to clusters
            allocate(x(size(rnp,2)), y(size(rnp,2)), z(size(rnp,2)))
            x = rnp(2,:); y = rnp(1,:); z = rnp(3,:)
            call fit_surface(x, y, z, fittype, pt, f, debug_outfile='./results/linecoeff_top')
            deallocate(x,y,z)

            !Get molecules on bottom surface
            call cluster_extents_grid(self, imaxloc(self%Nlist), 4, resolution, &
                                      extents_grid )!, debug_outfile='./results/maxcell_bot')
            call cluster_outer_mols(self, imaxloc(self%Nlist), tolerance=tolerance, dir=4, & 
                                    rmols=rnp, extents=extents_grid)!, debug_outfile='./results/clust_edge_bot')

            !Curve fits to clusters
            allocate(x(size(rnp,2)), y(size(rnp,2)), z(size(rnp,2)))
            x = rnp(2,:); y = rnp(1,:); z = rnp(3,:)
            call fit_surface(x, y, z, fittype, pb, f, debug_outfile='./results/linecoeff_bot')
            deallocate(x,y,z)

            !pb = (/ c, m, 0.0d0, 0.0d0  /)

            !No surface in x needed
            bintop(1) = 1e18
            binbot(1) = -1e18

            !Set a small CV in x, y and z at the surface 
            bintop(2) = 2.d0; binbot(2) = -2.d0
            bintop(3) = 2.d0; binbot(3) = -2.d0
            !pb = (/pt(1)-4.d0, 0.0d0, 0.00d0, 0.000d0  /)

            !Top/bottom surfaces in y
            !bintop(2) = 0.5d0*globaldomain(2) - tethereddisttop(2)
            !binbot(2) = -0.5d0*globaldomain(2) + tethereddistbottom(2)
            
            !Front/back surfaces in z
            !bintop(3) = 0.5d0*globaldomain(3) - tethereddisttop(3)
            !binbot(3) = -0.5d0*globaldomain(3) + tethereddistbottom(3)

            !Apply CV analysis to control volume with moving interface
            ptin = (/ 3.d0, -0.1d0, -0.01d0, -0.001d0  /) !pt(1:4)
            pbin = (/-3.d0, -0.1d0, -0.01d0, -0.001d0  /) !pb(1:4)
            !call cluster_CV_fn(pt, pb, bintop, binbot, 1)
            call cluster_CV_fn(ptin, pbin, bintop, binbot, 2)
            !call cluster_CV_fn(pt, pb, bintop, binbot, 3)

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


            !If cluster has broken up, stop simulation
            call check_for_cluster_breakup(self)


        else
            call error_abort('Error - Unknown intrinsic_interface_outflag')
        endif


        !if (iter .gt. 100200) stop


    end subroutine get_cluster_properties


    subroutine write_cluster_xyz(self, min_ngbr)
        use computational_constants_MD, only : extralloc, globaldomain
        use physical_constants_MD, only : np, nd
        use librarymod, only : get_Timestep_FileName, get_new_fileunit, imaxloc
        use interfaces, only : error_abort
        use arrays_MD, only : r, glob_no
        use module_record, only : ISR
        implicit none

        type(clusterinfo),intent(inout)    :: self
        integer, intent(in)                :: min_ngbr

        integer                         :: i, N, clustno, mainclusterNo, fileunit, Nrecords
        integer                         :: countwritten
        logical,save                    :: first_time=.true.

        integer, dimension(:), allocatable :: clust
        double precision,dimension(:,:),allocatable :: rnp, rnp_other, rnp_ex

        Nrecords = np! + extralloc
        mainclusterNo = imaxloc(self%Nlist)

        !Write vmd xyz file for debugging
        if (first_time) then
            fileunit = get_new_fileunit()
            open(fileunit, file="./all_clusters.xyz",status='replace')
            write(fileunit,*) Nrecords
            write(fileunit,*) ""
            first_time = .false.

        else
            fileunit = get_new_fileunit()
            open(fileunit, file="./all_clusters.xyz", access='append')
            write(fileunit,*) ""
            write(fileunit,*) Nrecords
        endif

        countwritten = 0
        do clustno =1,self%Nclust
            N = self%Nlist(clustno)
            if (N .ne. 0) then
                !print*, clustno, N, mainclusterNo
                if (clustno .eq. mainclusterNo) then
                    if (allocated(rnp)) deallocate(rnp)
                    if (allocated(rnp_ex)) deallocate(rnp_ex)
                    call cluster_to_array(self, clustno, r, min_ngbr, rnp, rnp_ex)
                    !Write interface first
                    if (allocated(ISR%pivots)) then
                        do i =1, size(ISR%pivots,1)
                            print*, "writing pivots", i, ISR%pivots(i), rnp(:,ISR%pivots(i))
                            write(fileunit,'(a,3f18.8)') "C", rnp(:,ISR%pivots(i))
                            countwritten = countwritten + 1
                        enddo
                    endif
                    !Then rest of cluster without interface
                    do i=1,size(rnp,2) !Cluster size can be larger than rnp as min_ngbr excluded
                        if (any(ISR%pivots .eq. i)) cycle
                        write(fileunit,'(a1,3f18.8)') "O", rnp(:,i)
                        countwritten = countwritten + 1
                    enddo
                    !countwritten = countwritten + size(rnp,2)
                    !Write excluded molecules by minimum neightbour cutoff
                    do i=1,size(rnp_ex,2) 
                        write(fileunit,'(a1,3f18.8)') "H", rnp_ex(:,i)
                    enddo
                    countwritten = countwritten + size(rnp_ex,2)

                else
                    if (allocated(rnp_other)) deallocate(rnp_other)
                    call cluster_to_array(self, clustno, r, 0, rnp_other)
                    do i=1,size(rnp_other,2)
                        write(fileunit,'(a1,3f18.8)') "B", rnp_other(:,i)
                    enddo
                    countwritten = countwritten + size(rnp_other,2)
                endif
            endif
        enddo

        !Sanity check, have all molecules been written either to main cluster or other
        if (countwritten .ne. Nrecords) then 
            print*, size(ISR%pivots,1), size(rnp,2), size(rnp_ex,2), size(rnp_other,2), countwritten,  Nrecords
            stop "Error in write_cluster_xyz - clust_main+clust_others .ne. total_particles"
        endif
        close(fileunit)

    end subroutine write_cluster_xyz



    subroutine write_cluster_xyz_multifile(self, min_ngbr)
        use computational_constants_MD, only : extralloc, globaldomain
        use physical_constants_MD, only : np, nd
        use librarymod, only : get_Timestep_FileName, get_new_fileunit, imaxloc
        use interfaces, only : error_abort
        use arrays_MD, only : r, glob_no
        use module_record, only : ISR
        implicit none

        type(clusterinfo),intent(inout)    :: self
        integer, intent(in)                :: min_ngbr

        integer                         :: i, N, clustno, mainclusterNo, fileunit, Nrecords
        integer                         :: fileunit1, fileunit2, fileunit3
        integer                         :: countwritten1, countwritten2, countwritten3
        logical,save                    :: first_time=.true.

        integer, dimension(:), allocatable :: clust
        double precision,dimension(:,:),allocatable :: rnp, rnp_ex

        Nrecords = np! + extralloc
        mainclusterNo = imaxloc(self%Nlist)

        !Write vmd xyz file for debugging
        if (first_time) then
            fileunit1 = get_new_fileunit()
            open(fileunit1, file="./cluster_main.xyz",status='replace')
            write(fileunit1,*) Nrecords
            write(fileunit1,*) ""
            fileunit2 = get_new_fileunit()
            open(fileunit2, file="./cluster_others.xyz",status='replace')
            write(fileunit2,*) Nrecords
            write(fileunit2,*) ""
            fileunit3 = get_new_fileunit()
            open(fileunit3, file="./cluster_interface.xyz",status='replace')
            write(fileunit3,*) Nrecords
            write(fileunit3,*) ""
            first_time = .false.

        else
            fileunit1 = get_new_fileunit()
            open(fileunit1, file="./cluster_main.xyz", access='append')
            write(fileunit1,*) Nrecords
            write(fileunit1,*) ""
            fileunit2 = get_new_fileunit()
            open(fileunit2, file="./cluster_others.xyz", access='append')
            write(fileunit2,*) Nrecords
            write(fileunit2,*) ""
            fileunit3 = get_new_fileunit()
            open(fileunit3, file="./cluster_interface.xyz", access='append')
            write(fileunit3,*) Nrecords
            write(fileunit3,*) ""
        endif

        countwritten1 = 0; countwritten2 = 0; countwritten3 = 0
        do clustno =1,self%Nclust
            N = self%Nlist(clustno)
            if (N .ne. 0) then
                !print*, clustno, N, mainclusterNo
                if (clustno .eq. mainclusterNo) then
                    call cluster_to_array(self, clustno, r, min_ngbr, rnp, rnp_ex)
                    do i=1,size(rnp,2) !Cluster size can be larger than rnp as min_ngbr excluded
                        write(fileunit1,'(a1,3f18.8)') "O", rnp(:,i)
!                        if (any(rnp(:,i)+0.5d0*globaldomain .gt. globaldomain)) then 
!                            print*, "write_cluster_xyz outside domain", i, rnp(:,i), clustno, np
!                            stop "Error in write_cluster_xyz"
!                        endif
                    enddo
                    countwritten1 = countwritten1 + size(rnp,2)
                    !Write excluded molecules by minimum neightbour cutoff to cluster_others.xyz file
                    do i=1,size(rnp_ex,2) 
                        write(fileunit2,'(a1,3f18.8)') "H", rnp_ex(:,i)
                    enddo
                    countwritten2 = countwritten2 + size(rnp_ex,2)
                    !print*, "CLUSTER SIZE ADDED = ", size(rnp,2) + size(rnp_ex,2), N
                    if (allocated(ISR%pivots)) then
                        do i =1, size(ISR%pivots,1)
                            !print'(a,i6,3f10.5,3i6)', "Cluster interface mols and bins", i, points(ISR%pivots(i),:),& 
                            !         ISR%get_bin(points(ISR%pivots(i),:), nbins, nhb)!, & 
                                 !ISR_mdt%get_bin(points(ISR%pivots(i),:), nbins, nhb)
                            write(fileunit3,'(a,3f18.8)') "C", rnp(:,ISR%pivots(i))
                        enddo
                        countwritten3 = countwritten3 + size(ISR%pivots,1)
                    endif
                else
                    call cluster_to_array(self, clustno, r, 0, rnp)
                    do i=1,size(rnp,2)
                        write(fileunit2,'(a1,3f18.8)') "H", rnp(:,i)
                    enddo
                    countwritten2 = countwritten2 + size(rnp,2)
                endif
                deallocate(rnp)
                if (allocated(rnp_ex)) deallocate(rnp_ex)
            endif
        enddo

        !Sanity check, have all molecules been written either to main cluster or other
        if (countwritten1 + countwritten2 .ne. Nrecords) then 
            print*, countwritten1, countwritten2, Nrecords
            stop "Error in write_cluster_xyz - clust_main+clust_others .ne. total_particles"
        endif

        do i=countwritten1+1, Nrecords
            write(fileunit1,'(a1,3f18.8)') "c", 0.d0, 0.d0, 0.d0
        enddo
        close(fileunit1)
        do i=countwritten2+1, Nrecords
            write(fileunit2,'(a1,3f18.8)') "c", 0.d0, 0.d0, 0.d0
        enddo
        close(fileunit2)
        do i=countwritten3+1, Nrecords
            write(fileunit3,'(a,3f18.8)') "c", 0.d0, 0.d0, 0.d0
        enddo
        close(fileunit3)

    end subroutine write_cluster_xyz_multifile



    subroutine surface_evolution(ISR, ISR_mdt, write_debug)
        use physical_constants_MD, only : np, halo_np
        use computational_constants_MD, only : iter, delta_t, nhb, halfdomain, &
											mflux_outflag, vflux_outflag, eflux_outflag, &
                                            Nmflux_ave, Nvflux_ave, Neflux_ave
        use librarymod, only : get_Timestep_FileName, get_new_fileunit
        use arrays_MD, only : r, v, a
        use module_set_parameters, only : mass, potenergymol
        use intrinsic_module, only : intrinsic_surface_real
        use calculated_properties_MD, only : binsize, nbins, &
					mass_surface_flux, momentum_surface_flux, energy_surface_flux
        use module_record, only : get_bin_molno
    	use CV_objects, only : CVcheck_mass, CV_debug
        implicit none

        logical, intent(in)                              :: write_debug
    	class(intrinsic_surface_real), intent(in)	     :: ISR, ISR_mdt

        integer                            :: n, i, b, pid, minbin, maxbin, ib, jb, kb
        integer                            :: n1, t1, t2
    	integer, save		               :: countmass=0, countmomentum=0, countenergy=0
        integer, parameter                 :: ct_mass=1, ct_momentum=2, ct_energy=3
        integer,dimension(3)               :: temp, bin, bin_mdt
        character(33)                      :: filename, debug_outfile
        double precision                   :: energy, crossdir, bintop, binbot
        double precision, dimension(2)     :: cross
        double precision, dimension(3)     :: ri, vi, velvect, rc, bs
        double precision, dimension(:), allocatable :: clustCV, s, smdt

        double precision, dimension(:,:), allocatable ::  points

        !Plot all molecules inside the liquid cluster control volume
        if (write_debug) then
            debug_outfile = './results/CV_mols'
            pid = get_new_fileunit()
            call get_Timestep_FileName(iter,debug_outfile,filename)
            open(unit=pid,file=trim(filename),status='replace')
        endif

        if (mflux_outflag .ne. 0) mass_surface_flux = 0.d0
        if (vflux_outflag .ne. 0) momentum_surface_flux = 0.d0

        n1=ISR_mdt%normal
        t1=ISR_mdt%ixyz
        t2=ISR_mdt%jxyz

        !Loop over all molecules and halos
        do n =1,np!+halo_np

			!Project velocity of particle forward to include 
			!forces applied to i before volume evolved
			ri(:) = r(:,n)
			vi(:) = v(:,n) + delta_t*a(:,n)

            !Get bins, only normal part can be different
            bin(:) = ISR%get_bin(ri)
            bin_mdt(:) = ISR_mdt%get_bin(ri)

            !If bin has changed, moving surface must have crossed them
            if (bin(n1) .ne. bin_mdt(n1)) then

                minbin = min(bin(n1), bin_mdt(n1))
                maxbin = max(bin(n1), bin_mdt(n1))
                crossdir  = sign(1.d0, dble(bin(n1)-bin_mdt(n1)))

                if (write_debug) then
                    write(pid,'(2i10,4f18.9)') minbin, maxbin, crossdir, rc
                endif
                if ((minbin .gt. 1) .and. (maxbin .lt. nbins(n1)+nhb(n1))) then
                    !More generally for non x evolutions
                    do i = minbin, maxbin-1
                        !Ensure within limits
                        if (i .gt. nbins(n1)+nhb(n1)) then
                            b = nbins(n1)+nhb(n1)
                        elseif (i .lt. 1 ) then
                            b = 2
                        else
                            b = i+1
                        endif					
							
						if (mflux_outflag .ne. 0) then
							mass_surface_flux(b, bin(2), bin(3), n1) = &
								mass_surface_flux(b, bin(2), bin(3), n1) + crossdir*mass(n)
							mass_surface_flux(b-1,bin(2),bin(3),n1+3) = & 
								mass_surface_flux(b,bin(2),bin(3),n1)
						endif
						if (vflux_outflag .ne. 0) then		
							momentum_surface_flux(b, bin(2), bin(3), :, n1) = &
								momentum_surface_flux(b, bin(2), bin(3), :, n1) + crossdir*vi
							momentum_surface_flux(b-1,bin(2),bin(3), :, n1+3) = & 
								momentum_surface_flux(b,bin(2),bin(3), :, n1)
						endif
						if (eflux_outflag .ne. 0) then
							velvect(:) = v(:,n) + 0.5d0*a(:,n)*delta_t
							energy = 0.5d0 * ( mass(n)*dot_product(velvect,velvect) + potenergymol(n))
							energy_surface_flux(b, bin(2), bin(3), n1) = &
								energy_surface_flux(b, bin(2), bin(3), n1) + crossdir*energy
							energy_surface_flux(b-1,bin(2),bin(3), n1+3) = & 
								energy_surface_flux(b,bin(2),bin(3), n1)
						endif
                    enddo
                endif
            endif

        enddo
        if (write_debug) then
            close(pid,status='keep')
        endif

        !Increment counters
	    countmass = countmass + 1
	    countmomentum = countmomentum + 1
	    countenergy = countenergy + 1
		if (mflux_outflag .ne. 0 .and. countmass .eq. Nmflux_ave) then
			call surface_evolution_mass_flux_io()
            countmass = 0
            mass_surface_flux = 0.d0
		endif
		if (vflux_outflag .ne. 0 .and. countmomentum .eq. Nvflux_ave) then		
			call surface_evolution_momentum_flux_io()
            countmomentum = 0
            momentum_surface_flux = 0.d0
		endif
		if (eflux_outflag .ne. 0 .and. countenergy .eq. Neflux_ave) then
			stop "surface_evolution - energy not developed"
			!call surface_evolution_energy_flux_io()
            countenergy = 0
            energy_surface_flux = 0.d0
		endif
		
    end subroutine surface_evolution

    ! A control volume with two intrinsic surfaces in the x directions
    ! and flat surfaces in the y and z directions
    subroutine CV_cluster_time(ISR, ISR_mdt, bintop, binbot, ri, cross)
#if ASSMBLY_HEAVISIDES
		use librarymod, only : heaviside  =>  heaviside_a1
#else
		use librarymod, only : heaviside
#endif     
        use intrinsic_module, only : intrinsic_surface_real
        implicit none

        double precision, dimension(3), intent(in)  :: ri
        double precision, intent(in)  :: bintop, binbot
    	class(intrinsic_surface_real), intent(in)	:: ISR, ISR_mdt
        double precision, dimension(2), intent(out) :: cross

        double precision     :: top, topmdt, bot, botmdt
        double precision, dimension(:), allocatable :: s, smdt      
        double precision, dimension(:,:), allocatable :: points
        
        !Left/Right cluster based surfaces in x {xsurf = f(yi)}
        top = bintop; bot = binbot

        allocate(points(1,3))
        points(1,1) = ri(1)
        points(1,2) = ri(2)
        points(1,3) = ri(3)
        call ISR%get_surface(points, s)
        call ISR_mdt%get_surface(points, smdt)

        !Use CV function to get top surface change
        top = bintop + s(1)
        topmdt = bintop + smdt(1)
        cross(1) = dble(( heaviside(top   -ri(1))   & 
                         -heaviside(topmdt-ri(1))))

        !if (cross(1) .ne. 0) print*, "surface moving cross = ", top, ri(1), topmdt, s(1), smdt(1)

        !Use CV function to get bottom surface change
        bot = binbot + s(1)
        botmdt = binbot + smdt(1)
        cross(2) = dble(( heaviside(bot   -ri(1))   & 
                         -heaviside(botmdt-ri(1))))

        !if (cross(2) .ne. 0) print*, "surface moving cross = ", bot, ri(1), botmdt, s(1), smdt(1)

    end subroutine CV_cluster_time




    !Plot density in a number of fluid bins aligned with the CV surface

!    subroutine CV_density_binning(p0)
!        use physical_constants_MD, only : tethereddisttop, tethereddistbottom
!        use computational_constants_MD, only : iter, globaldomain, delta_t
!        use physical_constants_MD, only : np, halo_np
!        use arrays_MD, only : r, v
!        use librarymod, only : get_Timestep_FileName, get_new_fileunit
!        use librarymod, only : heaviside  =>  heaviside_a1
!        implicit none

!        logical :: first_time = .true.
!        integer :: nbins, bin, n, pid
!        double precision,dimension(4),intent(in) :: p0
!        double precision                             :: surface, dx, width, yi, theta_i
!        double precision, dimension(3)  :: bintopi, binboti, ri
!        double precision, dimension(:), allocatable  :: surface_fitted_density
!        double precision, parameter     :: tol=1e-8

!        !Front/back surfaces in y
!        bintopi(2) =  0.5d0*globaldomain(2) - tethereddisttop(2)
!        binboti(2) = -0.5d0*globaldomain(2) + tethereddistbottom(2)
!        
!        !Left/Right cluser based surfaces in z
!        bintopi(3) =  0.5d0*globaldomain(3) - tethereddisttop(3)
!        binboti(3) = -0.5d0*globaldomain(3) + tethereddistbottom(3)

!        !Plot all molecules inside the liquid cluster control volume
!        nbins = 100
!        allocate(surface_fitted_density(nbins))
!        surface_fitted_density = 0.d0
!        width = 20.d0
!        dx = width/dble(nbins)    !Should use cluser extents
!        do n =1,np

!            ri(:) = r(:,n)
!            yi = ri(2)
!            surface = surface_fn(p0, yi)
!            do bin = 1,nbins
!                !Top/bottom surfaces in x
!                bintopi(1) = surface - dx*(bin-1) + 0.5*width
!                binboti(1) = surface - dx*bin     + 0.5*width

!                !Use CV function
!    	        theta_i = dble((heaviside(bintopi(1)-ri(1))-heaviside(binboti(1)-ri(1)))* & 
!    		               	   (heaviside(bintopi(2)-ri(2))-heaviside(binboti(2)-ri(2)))* & 
!    		              	   (heaviside(bintopi(3)-ri(3))-heaviside(binboti(3)-ri(3))))

!                if (abs(theta_i - 1.d0) .lt. tol) then
!                    !write(666600+iter,'(i8,5f10.5)'), bin, bintopi(1), binboti(1), ri
!                    surface_fitted_density(bin) = surface_fitted_density(bin) + 1.d0
!                    exit !Found the bin for this molecule, skip to the next molecule 
!            endif
!        enddo

!        enddo
!        !Write molecular density in slice to file
!        pid = get_new_fileunit()
!        if (first_time) then
!            open(unit=pid,file='./results/CV_binning',status='replace')
!            first_time = .false.
!        else
!            open(unit=pid,file='./results/CV_binning',access='append')
!        endif
!        write(pid,'(100f10.5)') surface_fitted_density
!        close(pid,status='keep')

!    end subroutine CV_density_binning



    subroutine check_for_cluster_breakup(self)
        use librarymod, only : bubble_sort_r
        use interfaces, only : error_abort
        implicit none

        type(clusterinfo),intent(inout)    :: self

        double precision,dimension(:),allocatable :: cluster_sizes

        ! Sort clusters by size and check if more than one big one!
        allocate(cluster_sizes(size(self%Nlist,1)))
        cluster_sizes = dble(self%Nlist)
        call bubble_sort_r(cluster_sizes)

        if ((cluster_sizes(1) - cluster_sizes(2))/cluster_sizes(1) .gt. 0.4d0) then

        else
            print*, 'It appears clusters have broken up'
            print'(a,8f10.1,e18.8)', 'CLUSTER DETAILS ', cluster_sizes(1:8), & 
                   (cluster_sizes(1) - cluster_sizes(2))/cluster_sizes(1)

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


    subroutine build_from_cell_list(self, cell, rd, & 
                                    rmols, nmols, skipwalls_)
	    use module_compute_forces, only: cellinfo, rj, rij, ri,&
                                         delta_rneighbr, rcutoff, rij2, &
                                         moltype
        use computational_constants_MD, only : Mie_potential, ncells
	    use interfaces, only : error_abort
        implicit none

        type(clusterinfo),intent(inout) :: self
        type(cellinfo),intent(in)       :: cell

        integer,intent(in)              :: nmols
        logical, intent(in),optional    :: skipwalls_
        double precision, intent(in)    :: rd
        double precision, intent(in), dimension(:,:) :: rmols

        logical                              :: skipwalls
	    integer                              :: i, j, m, n !Define dummy index
	    integer							     :: icell, jcell, kcell, Nchecked, Nnghbrs
        integer                              :: adjacentcellnp, icellshift, jcellshift, kcellshift
	    integer                              :: cellnp, molnoi, molnoj, noneighbrs
        double precision                     :: rd2

	    type(node), pointer 	             :: oldi, currenti, oldj, currentj

        if (present(skipwalls_)) then
            skipwalls = skipwalls_
        else
            skipwalls = .false.
        endif

        rd2 = rd**2

        self%clusterngbrs = 0
	    do kcell=2, ncells(3)+1
	    do jcell=2, ncells(2)+1
	    do icell=2, ncells(1)+1

		    cellnp = cell%cellnp(icell,jcell,kcell)
		    oldi => cell%head(icell,jcell,kcell)%point !Set old to first molecule in list

		    do i = 1,cellnp					!Step through each particle in list 

			    molnoi = oldi%molno 	 	!Number of molecule
			    ri = rmols(:,molnoi)         	!Retrieve ri
                
                Nchecked = 0
			    do kcellshift = -1,1
			    do jcellshift = -1,1
			    do icellshift = -1,1

				    oldj => cell%head(icell+icellshift, & 
                                      jcell+jcellshift, &
                                      kcell+kcellshift)%point
				    adjacentcellnp = cell%cellnp(icell+icellshift, & 
                                                 jcell+jcellshift, &
                                                 kcell+kcellshift)

				    do j = 1,adjacentcellnp			!Step through all j for each i

					    molnoj = oldj%molno			!Number of molecule
                        if (molnoj .gt. nmols) cycle
					    rj = rmols(:,molnoj)			!Retrieve rj

					    currentj => oldj
					    oldj => currentj%next		!Use pointer in datatype to obtain next item in list

					    if(molnoi==molnoj) cycle	!Check to prevent interaction with self

					    rij2=0						!Set rij^2 to zero
					    rij(:) = ri(:) - rj(:)		!Evaluate distance between particle i and j
					    rij2 = dot_product(rij,rij)	!Square of vector calculated

                        if (skipwalls .and. (Mie_potential .ne. 0)) then
                            if (any(moltype(molnoi) .eq. (/ 2, 9 /) )) cycle !Don't include wall molecules
                        endif

	                    if (rij2 .lt. rd2) then
                            !print*, molnoi, molnoj, rij2, self%Nclust
                            if (skipwalls .and. (Mie_potential .ne. 0)) then
                                if (any(moltype(molnoj) .eq. (/ 2, 9 /))) then
                                    call AddBondedPair(self, molnoi, molnoi)
                                else
                                    call AddBondedPair(self, molnoi, molnoj)
                                endif
                            else
                                call AddBondedPair(self, molnoi, molnoj)
                            endif
                            Nchecked = Nchecked + 1
                        endif
                    enddo
			    enddo
			    enddo
			    enddo
			    currenti => oldi
			    oldi => currenti%next !Use pointer in datatype to obtain next item in list

               !If no neighbours, add molecule to its own cluster list
                if (Nchecked .eq. 0) then
                    call AddBondedPair(self, molnoi, molnoi)
                    Nchecked = 1
                endif

                !Store number of clusterable molecules close to each molecule
                self%clusterngbrs(molnoi) = Nchecked
		    enddo
	    enddo
	    enddo
	    enddo


    end subroutine build_from_cell_list


    subroutine build_from_neighbour_list(self, neighbour, rd, & 
                                         rmols, nmols, skipwalls_)
	    use module_compute_forces, only: neighbrinfo, rj, rij, ri,&
                                         delta_rneighbr, rcutoff, rij2, &
                                         moltype
        use computational_constants_MD, only : Mie_potential
	    use interfaces, only : error_abort
        implicit none

        type(clusterinfo),intent(inout) :: self
        type(neighbrinfo),intent(in)    :: neighbour

        integer,intent(in)              :: nmols
        logical, intent(in),optional    :: skipwalls_
        double precision, intent(in)    :: rd
        double precision, intent(in), dimension(:,:) :: rmols

        logical                              :: skipwalls
	    integer                              :: i, j, m, n !Define dummy index
	    integer							     :: Nchecked, Nnghbrs
	    integer                              :: molnoi, molnoj, noneighbrs
        double precision                     :: rd2
	    type(node), pointer 	             :: old, current

        if (present(skipwalls_)) then
            skipwalls = skipwalls_
        else
            skipwalls = .false.
        endif

        rd2 = rd**2

        self%clusterngbrs = 0
        do molnoi = 1, nmols

	        ri = rmols(:,molnoi)         	!Retrieve ri
            if (skipwalls .and. (Mie_potential .ne. 0)) then
                if (any(moltype(molnoi) .eq. (/ 2, 9 /) )) cycle !Don't include wall molecules
            endif

            noneighbrs = neighbour%Nlist(molnoi)	    !Determine number of elements in neighbourlist
            old => neighbour%head(molnoi)%point		!Set old to head of neighbour list

            !A list of bonded pairs which we add to cluster linked list
            !only if they have more than min_ngbrs
            Nchecked = 0
            !Step through all neighbours i 
            do j = 1, noneighbrs

	            molnoj = old%molno			        !Number of molecule j
                if (molnoj .gt. nmols) cycle
	            rj(:) = rmols(:,molnoj)			            !Retrieve rj
	            rij(:)= ri(:) - rj(:)   	            !Evaluate distance between particle i and j
	            rij2  = dot_product(rij,rij)            !Square of vector calculated

	            if (rij2 .lt. rd2) then
                    if (skipwalls .and. (Mie_potential .ne. 0)) then
                        if (any(moltype(molnoj) .eq. (/ 2, 9 /))) then
                            call AddBondedPair(self, molnoi, molnoi)
                        else
                            call AddBondedPair(self, molnoi, molnoj)
                        endif
                    else
                        call AddBondedPair(self, molnoi, molnoj)
                    endif
                    Nchecked = Nchecked + 1
                endif
	            current => old
	            old => current%next !Use pointer in datatype to obtain next item in list
            enddo

            !If no neighbours, add molecule to its own cluster list
            if (Nchecked .eq. 0) then
                call AddBondedPair(self, molnoi, molnoi)
                Nchecked = 1
            endif

            !Store number of clusterable molecules close to each molecule
            self%clusterngbrs(molnoi) = Nchecked

        enddo

    end subroutine build_from_neighbour_list


    subroutine AddBondedPair(self, molnoi, molnoj)
        use linked_list, only : linklist_checkpushneighbr, linklist_merge
	    use interfaces, only : error_abort
        use arrays_MD, only : r
        implicit none

        type(clusterinfo),intent(inout)    :: self
        integer, intent(in)                :: molnoi, molnoj

        integer :: nc, nci, ncj, cbig, csmall, m
	    type(node), pointer 	        :: old, current

        !Special case adds one molecule only
        if (molnoi .eq. molnoj) then
            if (self%inclust(molnoi) .eq. 0) then
                self%Nclust = self%Nclust + 1
                !print'(a,3i5,3f10.5)', "NEw mol in add bonded pair", molnoi, self%Nclust, self%inclust(molnoi), r(:,molnoi)
                if (self%Nclust .gt. self%maxclusts) then
            		call error_abort("Increase maxcluster in clusterinfo linklist")
                endif
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
                if (self%Nclust .gt. self%maxclusts) then
            		call error_abort("Increase maxcluster in clusterinfo linklist")
                endif
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
                    !Now we've merge two linklists, it seems we should be able
                    !to reduce count by one BUT CANNOT as leaves a gap which
                    !we later remove with compress
                    !self%Nclust = self%Nclust - 1

                else
                    !If already in the same cluster, nothing to do
                endif
            endif
        endif

    end subroutine AddBondedPair

    !Remove any molecules from clusters with fewer than target 
    !numbers of neighbours
    !THIS DOES NOT WORK, NOT SIMPLE TO STRIP OUT MOLECULES
!    subroutine StripClusterSubNgbr(self, min_ngbr)
!        implicit none

!        type(clusterinfo),intent(inout) :: self
!        integer, intent(in)             :: min_ngbr

!        integer                         :: Nmols, molno, j, n
!	    type(node), pointer 	        :: old, current, pop, pold, pcurrent

!        !Loop though all clusters
!        do j = 1,self%Nclust
!            !For all clusters which are not empty
!            Nmols = self%Nlist(j)
!            if (Nmols .gt. 0) then
!                current => self%head(j)%point
!                do n = 1,Nmols
!                    molno = current%molno
!                    if (self%clusterngbrs(molno) .le. min_ngbr) then
!                        !Remove molno from cluster
!                        !Check if popped molecule is last in list
!                        pop => current
!                        if (associated(pop%next)) then
!                            !If just an element in list, remove it
!                            old     => pop%next        !Set old to next item in list
!                            current => pop%previous    !Set current to previous item in list
!                            deallocate(pop)				!Destroy pop

!                            !print*, n, Nmols, self%clusterngbrs(molno)

!                            old%previous => current   !Previous pointer connects to old (pop missed out)
!                            current%next => old     	!Next pointer connects to current (pop missed out)
!                        endif
!                    endif
!                    old => current%next
!                    current => old
!                enddo
!            endif
!        enddo

!    end subroutine StripClusterSubNgbr


    !Remove any gaps in the list of clusters
    subroutine CompressClusters(self)
        implicit none

        type(clusterinfo),intent(inout) :: self

        integer                         :: m, j, onc, nc
	    type(node), pointer 	        :: old, current

        !Loop though all clusters
        onc = self%Nclust !old number of clusters
        nc = 0
        !Loop though all clusters
        do j = 1,onc
            !For all clusters which are not empty
            if (self%Nlist(j) .gt. 0) then
                !print*, "CompressClusters", j, nc, self%Nclust, onc, self%Nlist(nc+1), self%Nlist(j)
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

        enddo
        self%Nclust = nc

!        do j = 1,self%Nclust
!            print*, "CompressClusters", j, self%Nclust, self%Nlist(j)
!        enddo

    end subroutine CompressClusters


    subroutine cluster_to_array(self, clustNo, qnty_array, min_ngbr, array, array_exlude)
        use computational_constants_MD, only : globaldomain
        use physical_constants_MD, only : np
        implicit none

        type(clusterinfo),intent(in)    :: self
        integer, intent(in)             :: clustNo, min_ngbr
		!qnty_array is the vector of molecular properties to get (r, v, a, molno)
        double precision, dimension(:,:), allocatable, intent(in) :: qnty_array
        double precision, dimension(:,:), allocatable, intent(out) :: array
        double precision, dimension(:,:), allocatable, intent(out), optional :: array_exlude

        integer                         :: n, m, o, Nmols
	    type(node), pointer 	        :: old, current
        double precision, dimension(:,:), allocatable :: temp, temp_ex


        !For all clusters which are not empty
        Nmols = self%Nlist(clustNo)
        if (Nmols .gt. 0) then
            allocate(temp(size(qnty_array,1),Nmols)); m = 0
            if (present(array_exlude)) then
                allocate(temp_ex(size(qnty_array,1),Nmols)); o = 0
            endif
            current => self%head(clustNo)%point
            !Loop through all cluster molecules         

            do n = 1,Nmols
                !Check not halo and minimum number of neighbours
                if (current%molno .le. np) then
                    if (self%clusterngbrs(current%molno) .ge. min_ngbr) then
                        m = m + 1
                        temp(:,m) = qnty_array(:,current%molno)
                        !if (any(r(:,current%molno)+0.5d0*globaldomain .gt. globaldomain)) then 
                        !    print*, "Cluster to array - outside domain", m, r(:,current%molno), current%molno
                            !stop "Error in write_cluster_xyz"
                        !endif
                    else
                        if (present(array_exlude)) then
                            o = o + 1
                            temp_ex(:,o) = qnty_array(:, current%molno)
                        endif
                        !print'(a,2i5,3f10.5,2(a,i3),a)', "Cluster to array ", clustNo, current%molno, r(:,current%molno), & 
                        !        " has ", self%clusterngbrs(current%molno) ," which is fewer than ", min_ngbr, " neighbours"
                    endif
                    !print*, n, m, current%molno, temp(:,m)
                endif
                old => current%next
                current => old
            enddo
        endif
        allocate(array(size(qnty_array,1),m)); array=0.d0
        array(:,:) = temp(:,1:m)
        if (present(array_exlude)) then
            allocate(array_exlude(size(qnty_array,1),o)); array_exlude=0.d0
            array_exlude(:,:) = temp_ex(:,1:o)
        endif

    end subroutine cluster_to_array

    subroutine print_cluster(self, clustNo)
        use arrays_MD, only : r
        use physical_constants_MD, only : np
        implicit none

        type(clusterinfo),intent(in)    :: self

        integer, intent(in)             :: clustNo
        integer                         :: n,Nmols
	    type(node), pointer 	        :: old, current

        !For all clusters which are not empty
        Nmols = self%Nlist(clustNo)
        if (Nmols .gt. 0) then
            current => self%head(clustNo)%point
            !Loop through all cluster molecules         
            do n = 1,Nmols
                if (n .le. np) then
                    print*, "molno = ", current%molno, " position = ", r(:,current%molno), np, Nmols
                endif
                old => current%next
                current => old
            enddo
        endif

    end subroutine print_cluster

    subroutine destroy_clusters(self)
        use linked_list, only : linklist_deallocate_cluster
        implicit none

        type(clusterinfo),intent(inout)    :: self

        call linklist_deallocate_cluster(self)

    end subroutine destroy_clusters



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


    subroutine cluster_outer_mols(self, clustNo, tolerance, dir, rmols, extents, debug_outfile)
        use arrays_MD, only : r
        use module_set_parameters, only : mass
        use computational_constants_MD, only : halfdomain, iter
	    use interfaces, only : error_abort
        use librarymod, only : get_new_fileunit, get_Timestep_FileName
        implicit none

        type(clusterinfo),intent(in)    :: self

        integer, intent(in)             :: clustNo, dir
        double precision,intent(in)     :: tolerance
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

        !If no extents supplied, use global cluster values
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

        !Get band of molecules within tolerance 
        if (dir .le. 3) then
            molband = extents_ - tolerance
        elseif (dir .gt. 3) then
            molband = extents_ + tolerance
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

                   ! print'(a, 2i6,6f10.5)', "Cluster", jcell, kcell, r(jxyz,molno), r(kxyz,molno), cellsidelength, clusterwidth

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
                        !print'(a, 2i6,5f10.5)', "Top Cluster", jcell, kcell, _extents(jcell,kcell),molband(jcell,kcell), r(:,molno)
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

        !Create dummy data
        !call build_debug_surface(rmols, global_extents)

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

    !Build a surface with the form x = a_{ij} * y^i * z^j 
    subroutine build_debug_surface(rmols, global_extents)
        use computational_constants_MD, only : halfdomain, iter
        use librarymod, only : get_new_fileunit, get_Timestep_FileName
        implicit none

        double precision,dimension(6),intent(in)         :: global_extents
        double precision,dimension(:,:),allocatable,intent(inout) :: rmols

        integer                         :: i,j,m,n
        double precision                :: rand
        double precision,dimension(4,4) :: a
        logical :: first_time = .true.

        rmols = 0.d0
        call random_number(a)
        a = a*0.00001d0
        a(1,1) = 14.d0

        n = get_new_fileunit()
        if (first_time) then
            open(unit=n,file="./results/debug_surface",status='replace')
            first_time = .false.
        else
            open(unit=n,file="./results/debug_surface",access='append')
        endif
        write(n,'(i8, 16f18.12)') iter, a 
        close(n,status='keep')
        
        do m=1,size(rmols,2)
            call random_number(rand)
            rmols(2,m) = 2.d0*rand*halfdomain(2)-halfdomain(2)
            call random_number(rand)
            rmols(3,m) = 2.d0*rand*halfdomain(3)-halfdomain(3)
            call random_number(rand)
            rmols(1,m) = 2.d0*rand !Add a little initial noise to surface
            do i=1,4
            do j=1,4
                rmols(1,m) = rmols(1,m) + a(i,j)*rmols(2,m)**(i-1) * rmols(3,m)**(j-1)
            enddo
            enddo
        enddo


    end subroutine build_debug_surface


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


end module cluster_analysis

subroutine get_interface_from_clusters()
    use cluster_analysis, only : build_clusters, & 
                                 get_cluster_properties, & 
                                 write_cluster_xyz, &
                                 destroy_clusters
    use linked_list, only : cluster
    use computational_constants_MD, only : CA_rd, CA_min_nghbr, CA_generate_xyz, & 
                                           iter, tplot, intrinsic_interface_outflag, &
                                           CV_conserve
    implicit none

    integer             :: min_ngbr
    double precision    :: rd

    min_ngbr = CA_min_nghbr !1000000
    rd = CA_rd !1.5d0
    if (mod(iter,tplot) .eq. 0 .or. CV_conserve .eq. 1) then
        call build_clusters(cluster, rd)
        if (any(intrinsic_interface_outflag .eq. (/1,2/))) then
            call get_cluster_properties(cluster, rd, min_ngbr)
        endif
        if (CA_generate_xyz .eq. 1) call write_cluster_xyz(cluster, min_ngbr)
        call destroy_clusters(cluster)
    endif

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
                    sum(density_bins( 2:size(density_bins,1)-1,   &
                                      2:size(density_bins,2)-1,   &
                                      2:size(density_bins,3)-1 ))/(Nmass_ave-1)
        endif
        print*, 'Densities', sum(density_bins(2:size(density_bins,1)-1,&
                                              2:size(density_bins,2)-1,&
                                              2:size(density_bins,3)-1)), &
                   product(shape(density_bins(2:size(density_bins,1)-1,&
                                              2:size(density_bins,2)-1,&
                                              2:size(density_bins,3)-1))),liquiddensity, gasdensity

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
        write(1984,'(i8,6f10.4,3f20.7)') iter, ri, rave(:), eigval
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

        rd2 = rd**2
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
                if (rd .le. rcutoff + minval(delta_rneighbr)) then

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
                    write(452,'(i5,12f10.4)') i, ri(:), surfacei(:,:,i)
                    deallocate(rarray)

                ! If the interface cutoff is greater than rcutoff
                ! then we can't use neighbour list and need to loop 
                ! over a greater number of adjacent cells
                elseif (rd .gt. rcutoff + minval(delta_rneighbr)) then

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


!===================================================================================
! Mass Flux over a surface of a bin
! Includes all intermediate bins

!subroutine cumulative_mass_flux_opt!(quantity, fluxes)
!	use module_record
!    use librarymod, only : CV_surface_flux, imaxloc!, heaviside  =>  heaviside_a1
!    use module_set_parameters, only : mass
!    !use CV_objects, only : CV_sphere_mass
!    implicit none

!!	real(kind(0.d0)),dimension(1), intent(in)	    :: quantity
!!	real(kind(0.d0)),dimension(1,6), intent(inout)	:: fluxes

!    logical, save                   :: first_time=.true.


!    logical                         :: crossings, use_bilinear
!	integer							:: jxyz,i,j,k,n,normal
!	integer		,dimension(3)		:: bin1,bin2,cbin,bs
!	real(kind(0.d0)),parameter		:: eps = 1e-12
!	real(kind(0.d0))        		:: crossdir
!	real(kind(0.d0)),dimension(3)	:: ri1,ri2,ri12, rci

!    integer, dimension(:), allocatable :: cbins
!    real(kind(0.d0)),dimension(:,:),allocatable :: points
!    real(kind(0.d0)),dimension(:,:), allocatable :: rc, rcx, rcy, rcz, dSdr

!    if (cluster_analysis_outflag .eq. 1 .and.  & 
!        any(intrinsic_interface_outflag .eq. (/1,2/))) then
!        use_bilinear = .true.
!    else
!        use_bilinear = .false.
!    endif

!	do n = 1,np!+halo_np

!		ri1(:) = r(:,n) 							!Molecule i at time t
!		ri2(:) = r(:,n)	- delta_t*v(:,n)			!Molecule i at time t-dt

!        if (use_bilinear) then
!            bin1(:) = ISR_mdt%get_bin(ri1, nbins, nhb)
!            bin2(:) = ISR_mdt%get_bin(ri2, nbins, nhb)
!        else
!            bin1 = get_bin(ri1)
!            bin2 = get_bin(ri2)
!        endif

!        do normal=1,3
!            if (use_bilinear .and. normal .eq. 1) then
!                !Redefine bins here
!                call get_crossings_bilinear(ri1, ri2, bin1, bin2, normal, rc, crossings, cbins)
!            else
!                call get_crossings(ri1, ri2, bin1, bin2, normal, rc, crossings)
!            endif
!            if (crossings) then
!                bs = 0
!                bs(normal) = 1
!	            ri12   = ri1 - ri2		        !Molecule i trajectory between t-dt and t
!	            where (abs(ri12) .lt. 0.000001d0) ri12 = 0.000001d0

!                if (normal .eq. 1 .and. use_bilinear) then
!                    if (allocated(points)) deallocate(points)
!                    allocate(points(size(rc,2), nd))
!                    points(:,1) = rc(1,:)
!                    points(:,2) = rc(2,:)
!                    points(:,3) = rc(3,:)
!                    call ISR_mdt%get_surface_derivative(points, dSdr)
!                endif

!                !print'(a,i9,8i5,12f10.5)', "Ncrossings ", iter, normal, bin1, bin2, size(rc,2), ri1, rc(:,1), rc(:,2), ri2
!                do i =1,size(rc,2)
!                    rci = rc(:,i)
!                    rci(normal) = rci(normal) + eps
!                    cbin(:) = ISR_mdt%get_bin(rci, nbins, nhb)
!                    if (normal .eq. 1 .and. use_bilinear) then
!                        cbin(normal) = cbins(normal)+1
!                        crossdir = sign(1.d0,(ri12(1) - ri12(2)*dSdr(i,1) - ri12(3)*dSdr(i,2)))
!                    else
!!                        cbin(:) =  get_bin(rci)
!                        !print'(a,i9,12i5,9f10.5)', "Ncrossings ", iter, normal, i, bin1, bin2, cbin, size(rc,2), ri1, rci, ri2
!                        crossdir  = sign(1.d0, ri12(normal))
!                    endif

!                    mass_flux(cbin(1),cbin(2),cbin(3),normal) = & 
!                        mass_flux(cbin(1),cbin(2),cbin(3),normal) + crossdir*mass(n)
!                    mass_flux(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),normal+3) = & 
!                        mass_flux(cbin(1),cbin(2),cbin(3),normal)

!                enddo
!                deallocate(rc)
!            endif
!        enddo

!	enddo

!end subroutine cumulative_mass_flux_opt




!subroutine cumulative_mass_flux_bilinear_only!(quantity, fluxes)
!	use module_record
!    use librarymod, only : CV_surface_flux, imaxloc!, heaviside  =>  heaviside_a1
!    use module_set_parameters, only : mass
!    !use CV_objects, only : CV_sphere_mass
!    implicit none

!!	real(kind(0.d0)),dimension(1), intent(in)	    :: quantity
!!	real(kind(0.d0)),dimension(1,6), intent(inout)	:: fluxes

!    logical, save                   :: first_time=.true.


!    logical                         :: crossings, use_bilinear
!	integer							:: jxyz,i,j,k,n,normal
!	integer		,dimension(3)		:: bin1,bin2,cbin,bs
!	real(kind(0.d0)),parameter		:: eps = 1e-12
!	real(kind(0.d0))        		:: crossdir
!	real(kind(0.d0)),dimension(3)	:: ri1,ri2,ri12, rci

!    integer, dimension(:), allocatable :: cbins
!    real(kind(0.d0)),dimension(:,:),allocatable :: points
!    real(kind(0.d0)),dimension(:,:), allocatable :: rc, rcx, rcy, rcz, dSdr

!    if (cluster_analysis_outflag .eq. 1 .and.  & 
!        any(intrinsic_interface_outflag .eq. (/1,2/))) then

!	    do n = 1,np

!		    ri1(:) = r(:,n) 							!Molecule i at time t
!		    ri2(:) = r(:,n)	- delta_t*v(:,n)			!Molecule i at time t-dt

!            normal=ISR_mdt%normal
!            bin1(:) = ISR_mdt%get_bin(ri1, nbins, nhb)
!            bin2(:) = ISR_mdt%get_bin(ri2, nbins, nhb)
!            call get_crossings_bilinear(ri1, ri2, bin1, bin2, normal, rc, crossings, cbins)

!            if (crossings) then
!                bs = 0
!                bs(normal) = 1
!                ri12   = ri1 - ri2		        !Molecule i trajectory between t-dt and t
!                where (abs(ri12) .lt. 0.000001d0) ri12 = 0.000001d0

!                if (allocated(points)) deallocate(points)
!                allocate(points(size(rc,2), nd))
!                points(:,1) = rc(1,:)
!                points(:,2) = rc(2,:)
!                points(:,3) = rc(3,:)
!                call ISR_mdt%get_surface_derivative(points, dSdr)

!                do i =1,size(rc,2)
!                    rci = rc(:,i)
!                    rci(normal) = rci(normal) + eps 

!                    cbin(:) = ISR_mdt%get_bin(rci, nbins, nhb)
!                    cbin(normal) = cbins(normal)+1
!                    crossdir = sign(1.d0,(ri12(1) - ri12(2)*dSdr(i,1) - ri12(3)*dSdr(i,2)))

!                    mass_flux(cbin(1),cbin(2),cbin(3),normal) = & 
!                        mass_flux(cbin(1),cbin(2),cbin(3),normal) + crossdir*mass(n)
!                    mass_flux(cbin(1)-bs(1),cbin(2)-bs(2),cbin(3)-bs(3),normal+3) = & 
!                        mass_flux(cbin(1),cbin(2),cbin(3),normal)

!                enddo
!                deallocate(rc)
!            endif
!        enddo
!    endif


!end subroutine cumulative_mass_flux_bilinear_only








