!-----------------------------------------------------------------------------
!
!							   Set Parameters
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
	use polymer_info_MD
	use concentric_cylinders
	use librarymod, only : PDF
    implicit none

	type(PDF) 									:: velPDF, velPDFMB
	type(PDF),allocatable,dimension(:,:,:,:) 	:: velPDF_array

    !LJ parameters
    double precision           :: potshift !Shift in Lennard Jones potential due to cutoff

    !Generalised Mie-potential parameters
    integer,parameter :: ntypes = 9
    character(30),dimension(ntypes)             :: moltype_names
    double precision,dimension(ntypes)          :: mass_lookup
    double precision,dimension(ntypes,ntypes)   :: epsilon_lookup, sigma_lookup, &    
                                                   lambdar_lookup, lambdaa_lookup, &
                                                   C_lookup, potshift_lookup, &
                                                   k_lookup, r0_lookup, equil_sep_lookup, &
                                                   alpha_lookup

    double precision,dimension(ntypes,ntypes,ntypes) :: angular_k_lookup, angular_r0_lookup

	abstract interface
		function fn_mass(i)
			double precision              :: fn_mass
            integer, intent(in)           :: i
		end function fn_mass
	end interface

	abstract interface
		function fn_accijmag(invrij2, i, j)
			double precision              :: fn_accijmag
            integer, intent(in)           :: i, j
			double precision, intent (in) :: invrij2
		end function fn_accijmag
	end interface

	abstract interface
		function fn_force(invrij2, rij, i, j)
			double precision,dimension(3) :: fn_force
            integer, intent(in)           :: i, j
            double precision, intent(in),dimension(3)   :: rij
			double precision, intent (in) :: invrij2
		end function fn_force
	end interface

	abstract interface
		function fn_energy(invrij2, i, j)
			double precision              :: fn_energy
            integer, intent(in)           :: i, j
			double precision, intent (in) :: invrij2
		end function fn_energy
	end interface

   procedure (fn_mass),     pointer :: mass         => null ()
   procedure (fn_accijmag), pointer :: get_accijmag => null ()
   procedure (fn_force),    pointer :: get_force    => null ()
   procedure (fn_energy),   pointer :: get_energy   => null ()

   procedure (fn_accijmag), pointer :: get_poly_accijmag => null ()
   procedure (fn_force),    pointer :: get_poly_force    => null ()
   procedure (fn_energy),   pointer :: get_poly_energy   => null ()

contains

    !LJ or Mie force calculation functions
    function LJ_mass(i)

        integer, intent(in)             :: i
        double precision                :: LJ_mass

        LJ_mass = 1.d0

    end function LJ_mass

    function LJ_accijmag(invrij2, i, j)

        integer, intent(in)             :: i, j
        double precision, intent(in)    :: invrij2
        double precision                :: LJ_accijmag

        LJ_accijmag = 48.d0 * ( invrij2**7 - 0.5d0*invrij2**4 )

    end function LJ_accijmag

    function LJ_force(invrij2, rij, i, j)

        integer, intent(in)                         :: i, j
        double precision, intent(in)                :: invrij2
        double precision,dimension(3), intent(in)   :: rij
        double precision,dimension(3)               :: LJ_force

        LJ_force = LJ_accijmag(invrij2, i, j)*rij

    end function LJ_force

    function LJ_energy(invrij2,  i, j)

        integer, intent(in)             :: i, j
        double precision, intent(in)    :: invrij2
        double precision                :: LJ_energy

        LJ_energy = 4.d0*( invrij2**6 - invrij2**3 )-potshift

    end function LJ_energy




    function Mie_mass(i)
        use arrays_MD, only : moltype

        integer, intent(in)             :: i
        double precision                :: Mie_mass

        Mie_mass = mass_lookup(moltype(i))

    end function Mie_mass

    function Mie_accijmag(invrij2, i, j)
        use arrays_MD, only : moltype

        integer, intent(in)             :: i, j
        double precision, intent(in)    :: invrij2
        double precision                :: Mie_accijmag

        double precision                :: C, sigmaij, epsilonij, &
                                              lambdar, lambdaa, alpha

        epsilonij = epsilon_lookup(moltype(i),moltype(j))
        sigmaij   = sigma_lookup(moltype(i),moltype(j))
        lambdar   = lambdar_lookup(moltype(i),moltype(j))
        lambdaa   = lambdaa_lookup(moltype(i),moltype(j))
        C         = C_lookup(moltype(i),moltype(j))
        alpha     = alpha_lookup(moltype(i),moltype(j))

        Mie_accijmag = C*epsilonij*(       lambdar*invrij2**(0.5d0*lambdar+1) & 
                                    -alpha*lambdaa*invrij2**(0.5d0*lambdaa+1) )
    end function Mie_accijmag

    function Mie_force(invrij2, rij, i, j)

        integer, intent(in)                         :: i, j
        double precision, intent(in)                :: invrij2
        double precision, intent(in),dimension(3)   :: rij
        double precision,dimension(3)               :: Mie_force

        Mie_force = rij*Mie_accijmag(invrij2, i, j)

    end function Mie_force


    function Mie_energy(invrij2, i, j)
        use arrays_MD, only : moltype

        integer, intent(in)             :: i, j
        double precision, intent(in)    :: invrij2
        double precision                :: Mie_energy

        double precision                :: C, sigmaij, epsilonij, &
                                           lambdar, lambdaa, alpha,&
                                           potshift

        epsilonij = epsilon_lookup(moltype(i),moltype(j))
        sigmaij   = sigma_lookup(moltype(i),moltype(j))
        lambdar   = lambdar_lookup(moltype(i),moltype(j))
        lambdaa   = lambdaa_lookup(moltype(i),moltype(j))
        C         = C_lookup(moltype(i),moltype(j))
        potshift  = potshift_lookup(moltype(i),moltype(j))
        alpha     = alpha_lookup(moltype(i),moltype(j))

        Mie_energy = C*epsilonij*(      invrij2**(0.5d0*lambdar) & 
                                 -alpha*invrij2**(0.5d0*lambdaa) ) - potshift

    end function Mie_energy


    !Functions for FENE
    function FENE_accijmag(rij2, i, j)
        use polymer_info_MD, only : k_c, R_0

        integer, intent(in)             :: i, j
        double precision, intent(in)    :: rij2
        double precision                :: FENE_accijmag

		if(rij2.ge.R_0**2)	call polymer_bond_error(i,j)
        FENE_accijmag =  -k_c/(1-(rij2/(R_0**2)))

    end function FENE_accijmag

    function FENE_force(rij2, rij, i, j)

        integer, intent(in)                         :: i, j
        double precision, intent(in)                :: rij2
        double precision, intent(in),dimension(3)   :: rij
        double precision,dimension(3)               :: FENE_force

        FENE_force = rij*FENE_accijmag(rij2, i, j)

    end function FENE_force

    function FENE_energy(rij2, i, j)
        use polymer_info_MD, only : k_c, R_0

        integer, intent(in)             :: i, j
        double precision, intent(in)    :: rij2
        double precision                :: FENE_energy

        FENE_energy = 0.5d0*k_c*R_0*R_0*dlog(1.d0-(rij2/(R_0**2)))

    end function FENE_energy

    !Functions for harmonic potential
    function harmonic_accijmag(rij2, i, j)
        use arrays_MD, only : moltype

        integer, intent(in)             :: i, j
        double precision, intent(in)    :: rij2
        double precision                :: k_harmonic, r0
        double precision                :: rij_mag, harmonic_accijmag

        k_harmonic = k_lookup(moltype(i),moltype(j))
        r0 = r0_lookup(moltype(i),moltype(j))

        rij_mag = sqrt(rij2)
        harmonic_accijmag = -k_harmonic * (rij_mag - r0)/rij_mag

    end function


    function harmonic_force(rij2, rij, i, j)
        use arrays_MD, only : moltype

        integer, intent(in)             :: i, j
        double precision, intent(in)    :: rij2
        double precision,dimension(3), intent(in)    :: rij
        double precision                :: k_harmonic, r0
        double precision,dimension(3)   :: harmonic_force

        k_harmonic = k_lookup(moltype(i),moltype(j))
        r0 = r0_lookup(moltype(i),moltype(j))
        harmonic_force = rij * harmonic_accijmag(rij2, i, j) 

       ! print'(4i6,3f10.5,3f16.2)', i,j,moltype(i),moltype(j),k_harmonic,r0,sqrt(rij2),harmonic_force

    end function harmonic_force

    function harmonic_energy(rij2, i, j)
        use arrays_MD, only : moltype

        integer, intent(in)             :: i, j
        double precision, intent(in)    :: rij2
        double precision                :: harmonic_energy, rij_mag
        double precision                :: k_harmonic, r0

        k_harmonic = k_lookup(moltype(i),moltype(j))
        r0 = r0_lookup(moltype(i),moltype(j))

        rij_mag = sqrt(rij2)
        harmonic_energy = 0.5d0*k_harmonic*(rij_mag-r0)**2.d0
        !harmonic_energy = 0.5d0*k_harmonic*(dot_product(rij-r0,rij-r0))

    end function harmonic_energy


    ! Add function for harmonic angular potential
    ! for polymer of the form   o-o-o-o-o
    !                             ^ ^ ^
    !                             | | |
    ! we will check               i j k 

    function angular_harmonic_force(rij, rjk, i, j, k)

        integer, intent(in)                         :: i, j, k
        double precision, intent(in),dimension(3)   :: rij, rjk
        double precision                            :: k_ijk, theta_0
        double precision                            :: mag_rij, mag_rjk
        double precision                            :: theta_ijk, dot_rijrjk, accijmag
        double precision,dimension(3,3)             :: angular_harmonic_force

        k_ijk = angular_k_lookup(moltype(i),moltype(j),moltype(k))
        theta_0 = angular_r0_lookup(moltype(i),moltype(j),moltype(k))

        !Check if no angular interaction present between molecules
        if (k_ijk .lt. 1e-5) return

        mag_rij = sqrt(dot_product(rij,rij))
        mag_rjk = sqrt(dot_product(rjk,rjk))
        theta_ijk = acos(dot_product(rij,rjk)/(mag_rij*mag_rjk))

        if (theta_ijk .lt. 1e-5) then
            angular_harmonic_force = 0.d0
            return
        endif



        !Force on molecule i
        angular_harmonic_force(1,:) =-(k_ijk * (theta_ijk - theta_0) & 
                                      / (mag_rij * sin(theta_ijk))) &
                                       * ((rjk/mag_rjk) - (rij/mag_rij)& 
                                       *cos(theta_ijk))
        !Force on molecule k
        angular_harmonic_force(3,:) =-(k_ijk * (theta_ijk - theta_0) & 
                                      / (mag_rjk * sin(theta_ijk))) & 
                                       * ((rij/mag_rij) - (rjk/mag_rjk)& 
                                       *cos(theta_ijk))
        !Force on central molecule j
        angular_harmonic_force(2,:) = - angular_harmonic_force(1,:) & 
                                      - angular_harmonic_force(3,:)

!        accijmag = -(k_ijk * (theta_ijk - theta_0)/(sin(theta_ijk)))
!        angular_harmonic_force(1,:)=accijmag* ((rjk/mag_rjk) & 
!                                              -(rij/mag_rij)*cos(theta_ijk)) & 
!                                                /mag_rij
!        angular_harmonic_force(3,:)=accijmag* ((rij/mag_rij) & 
!                                              -(rjk/mag_rjk)*cos(theta_ijk)) & 
!                                                /mag_rjk
!        angular_harmonic_force(2,:) = - angular_harmonic_force(1,:) & 
!                                      - angular_harmonic_force(3,:)

    end function angular_harmonic_force



    function angular_harmonic_energy(rij, rjk, i,j,k)

        integer, intent(in)                         :: i, j, k
        double precision, intent(in),dimension(3)   :: rij, rjk
        double precision                            :: k_ijk, theta_0
        double precision                            :: mag_rij, mag_rjk
        double precision                            :: theta_ijk, dot_rijrjk
        double precision                            :: angular_harmonic_energy

        k_ijk = angular_k_lookup(moltype(i),moltype(j),moltype(k))
        theta_0 = angular_r0_lookup(moltype(i),moltype(j),moltype(k))

        mag_rij = sqrt(dot_product(rij,rij))
        mag_rjk = sqrt(dot_product(rjk,rjk))
        theta_ijk = acos(dot_product(rij,rjk)/(mag_rij*mag_rjk))
        angular_harmonic_energy = 0.5d0 * k_ijk * (theta_ijk - theta_0)**2

    end function angular_harmonic_energy	


    subroutine polymer_bond_error(molnoi, molnoX)
		use interfaces
		implicit none

		integer, intent(in) :: molnoi, molnoX
		real(kind(0.d0)) :: rglobi(3),rglobX(3), rij2

		rglobi(1) = r(1,molnoi) - halfdomain(1)*(npx-1) + domain(1)*(iblock - 1)   
		rglobi(2) = r(2,molnoi) - halfdomain(2)*(npy-1) + domain(2)*(jblock - 1)   
		rglobi(3) = r(3,molnoi) - halfdomain(3)*(npz-1) + domain(3)*(kblock - 1)   

		rglobX(1) = r(1,molnoX) - halfdomain(1)*(npx-1) + domain(1)*(iblock - 1)   
		rglobX(2) = r(2,molnoX) - halfdomain(2)*(npy-1) + domain(2)*(jblock - 1)   
		rglobX(3) = r(3,molnoX) - halfdomain(3)*(npz-1) + domain(3)*(kblock - 1)   

        rij2 = dot_product(r(:,molnoi),r(:,molnoX))

		print*, 'irank: ', irank
		print '(a,i6,a,i8,a,i8,a,f8.5,a,f8.5,a,f12.5)', & 
					'Bond broken at iter ',iter,': atoms ',molnoi,' and ',molnoX,' are separated by ', &
					rij2**0.5,', which is greater than the allowed limit of ', R_0, &
					'. Stopping simulation, total time elapsed = ', iter*delta_t
		print '(a)', 'Atomic positions:'
		print '(a,i8,a,f10.5,a,f10.5,a,f10.5)', 'Atom ',molnoi,' is located at global position', &
		                                        rglobi(1),' ',rglobi(2),' ',rglobi(3) 
		print '(a,i8,a,f10.5,a,f10.5,a,f10.5)', 'Atom ',molnoX,' is located at global position', &
		                                        rglobX(1),' ',rglobX(2),' ',rglobX(3) 

		print '(a,i8,a)', 'Monomer information for atom ', molnoi,':'
		print '(a,i8)', 'ChainID: '   , monomer(molnoi)%chainID
		print '(a,i8)', 'SubchainID: ', monomer(molnoi)%subchainID
		print '(a,i8)', 'Funcy: '     , monomer(molnoi)%funcy
		print '(a,i8)', 'Glob_no: '   , monomer(molnoi)%glob_no
		print '(a,i8)', 'Bin_bflag: ' , monomer(molnoi)%bin_bflag
		print '(a,i8)', 'Tag: '       , tag(molnoi) 

		print '(a,i8,a)', 'Monomer information for atom ', molnoX,':'
		print '(a,i8)', 'ChainID: '   , monomer(molnoX)%chainID
		print '(a,i8)', 'SubchainID: ', monomer(molnoX)%subchainID
		print '(a,i8)', 'Funcy: '     , monomer(molnoX)%funcy
		print '(a,i8)', 'Glob_no: '   , monomer(molnoX)%glob_no
		print '(a,i8)', 'Bin_bflag: ' , monomer(molnoX)%bin_bflag
		print '(a,i8)', 'Tag: '       , tag(molnoX) 

		if (molnoX.gt.np) then
			call error_abort('Halo!')
		else
			call error_abort('')
		end if

	end subroutine polymer_bond_error


subroutine setup_mie_potential
    use interfaces, only : error_abort
    implicit none

    integer :: i, j, ids(4)


    ! ------------Mie Potential--------------
    ! phi = A(lambda_r,lambda_a) * epsilon * [ (sigma/r)^lambda_r - (sigma/r)^lambda_a)]
    ! A(lambda_r,lambda_a) = (lambda_r/(lambda_r-lambda_a))*(lambda_r/lambda_a)^(lambda_a/(lambda_r-lambda_a))
    ! 
    ! ------------Including two phase repulsion between species --------------
    ! phi = A(lambda_r,lambda_a) * epsilon * [ (sigma/r)^lambda_r - alpha * (sigma/r)^lambda_a)]
    ! Where alpha is the relative magnitude of attrative and repulsive terms (allows opposite signs)

    !Initialise all as zero
    mass_lookup = 0.d0
    epsilon_lookup = 0.d0
    sigma_lookup   = 0.d0
    lambdar_lookup = 0.d0
    lambdaa_lookup = 0.d0
    alpha_lookup = 1.d0

    !Setup table of cross potentials 
    !1 == Argon;
    !moltype_names(1) = '           Ar           '
    moltype_names(1)    = 'Ar' !' Ar '
    mass_lookup(1)      = 1.0d0
    epsilon_lookup(1,1) = 1.d0
    sigma_lookup(1,1)   = 1.d0
    lambdar_lookup(1,1) = 12.d0
    lambdaa_lookup(1,1) = 6.d0

    !2 == Wall; 
    !moltype_names(2) = '          Wall          '
    moltype_names(2)    = 'S' !' Wall '
    mass_lookup(2)      = 1.d0
    epsilon_lookup(2,2) = 1.d0
    sigma_lookup(2,2)   = 1.d0
    lambdar_lookup(2,2) = 12.d0
    lambdaa_lookup(2,2) = 6.d0


    !Liquid Argon and wall
    epsilon_lookup(2,1) = eij_wall(1)
    epsilon_lookup(8,2) = eij_wall(1)
    epsilon_lookup(9,1) = eij_wall(2)
    epsilon_lookup(9,8) = eij_wall(2)

    !1-2 == Wall/{D,M,CM} hydrophobic/strong wall interaction
    ids = (/ 3,4,5,7 /)
    do i =1,size(ids)
        epsilon_lookup(ids(i),2) = eij_wall(1)
        epsilon_lookup(9,ids(i)) = eij_wall(2)
    enddo

    !select case (wall_liquid)
    !case(No_wetting) 
    !    epsilon_lookup(2,1) = 1.0d0
    !    do i =1,size(ids)
    !        epsilon_lookup(ids(i),2) = 1.0d0
    !    enddo
    !case(Paraffin_Water) 
    !    epsilon_lookup(2,1) = 0.5d0
    !    do i =1,size(ids)
    !        epsilon_lookup(ids(i),2) = 0.5d0
    !    enddo
    !case(Superhydrophobic)
    !    do i =1,size(ids)
    !        epsilon_lookup(ids(i),2) = 0.01d0
    !    enddo
    !case(Superspreading)
    !    !Superspreading requires this as 1.4 according to Panos
         !for water and hydrophobic parts CM, M and D

    !1-2 == Wall/Water [hydrophilic or weak wall interaction]
    ! set to same as surfactant
    !epsilon_lookup(3,2) = epsilon_lookup(ids(1),2)
    !epsilon_lookup(9,3) = epsilon_lookup(9,ids(1))

    !Wall/EO [hydrophillic or weak wall interaction]
    !epsilon_lookup(6,2)  = 0.5d0     
    !case(Cross_rules)
    !end select

    !end select


    ! SAFT gamma Mie values ₀₁₂₃₄₅₆₇₈₉	
    !(from Theodorakis, Muller, Craster & Matar (2014))
    !3 == Water "W" 2{H₂O} molecules per bead
    !moltype_names(3) = '         2H₂O          '
    moltype_names(3)    = 'W'!' 2H2O '
    mass_lookup(3)      = 0.8179d0
    epsilon_lookup(3,3) = 0.8129d0
    sigma_lookup(3,3)   = 0.8584d0
    lambdar_lookup(3,3) = 8.d0
    lambdaa_lookup(3,3) = 6.d0

    !4 == SAFT "M" {(CH₃)₃--Si--O½} molecules per bead
    !moltype_names(4) = '     (CH₃)₃--Si--O½     '
    moltype_names(4)    = 'M'!' 3CH3SiOh '
    mass_lookup(4)      = 1.8588d0
    epsilon_lookup(4,4) = 0.8998d0
    sigma_lookup(4,4)   = 1.2398d0
    lambdar_lookup(4,4) = 26.d0
    lambdaa_lookup(4,4) = 6.d0

    !5 == SAFT "D" {O½--(CH₃)₂--Si--O½} molecules per bead
    !moltype_names(5) = '   O½--(CH₃)₂--Si--O½   '
    moltype_names(5)    = 'D'!' Oh2CH3SiOh '
    mass_lookup(5)      = 1.6833d0
    epsilon_lookup(5,5) = 0.5081d0
    sigma_lookup(5,5)   = 1.0702d0
    lambdar_lookup(5,5) = 13.90d0
    lambdaa_lookup(5,5) = 6.d0

    !6 == SAFT ether "EO" {--CH₂--O--CH₂--} molecules per bead
    !moltype_names(6)   = '    --CH₂--O--CH₂--     '
    moltype_names(6)    = 'EO'!' CH2OCH2 '
    mass_lookup(6)      = 1.0000d0
    epsilon_lookup(6,6) = 0.8067d0
    sigma_lookup(6,6)   = 0.9307d0
    lambdar_lookup(6,6) = 19.d0
    lambdaa_lookup(6,6) = 6.d0

    !7 == SAFT alkane "CM" {--CH₂--CH₂--CH₂--} molecules per bead
    !moltype_names(7) = '   --CH₂--CH₂--CH₂--    '
    moltype_names(7)    = 'CM'!' 3CH2 '
    mass_lookup(7)      = 0.9552d0
    epsilon_lookup(7,7) = 0.7000d0
    sigma_lookup(7,7)   = 1.0000d0
    lambdar_lookup(7,7) = 15.d0
    lambdaa_lookup(7,7) = 6.d0

    !8 == Second phase of Argon;
    !moltype_names(1) = '           Ar           '
    moltype_names(8)    = 'rA' !' Ar '
    mass_lookup(8)      = 1.0d0
    epsilon_lookup(8,8) = 1.d0
    sigma_lookup(8,8)   = 1.d0
    lambdar_lookup(8,8) = 12.d0
    lambdaa_lookup(8,8) = 6.d0

    !The two phases of argon attract each other less strongly
    alpha_lookup(1,8) = -1.d0
    alpha_lookup(8,1) = -1.d0

    !9 == Wall type 2 ; 
    !moltype_names(9) = '          Wall          '
    moltype_names(9)    = 'S2' !' Wall '
    mass_lookup(9)      = 1.d0
    epsilon_lookup(9,9) = 1.d0
    sigma_lookup(9,9)   = 1.d0
    lambdar_lookup(9,9) = 12.d0
    lambdaa_lookup(9,9) = 6.d0

    !Define adjusted cross potential interactions (tuned by prior simulation)
    !ether and Water
    epsilon_lookup(6,3) = 0.9756d0
    !SAFT adjusted water--CH3 interaction from prior studies
    epsilon_lookup(7,3) = 0.5081d0
    !SAFT adjusted D--M interaction from prior studies
    epsilon_lookup(5,4) = 0.7114d0
    !SAFT adjusted alkane--ether interaction from prior studies
    epsilon_lookup(7,6) = 0.7154d0

    !Define chain interactions

    !Default to zero for anything which shouldn't interact!
    k_lookup = 0.d0 ;         r0_lookup = 0.d0
    k_lookup(5,4) = 295.3322; r0_lookup(5,4) = 1.1550
    k_lookup(6,5) = 295.3322; r0_lookup(6,5) = 1.0004
    k_lookup(6,6) = 295.3322; r0_lookup(6,6) = 0.9307
    k_lookup(7,6) = 295.3322; r0_lookup(7,6) = 0.9653
    k_lookup(7,7) = 295.3322; r0_lookup(7,7) = 1.0000
    !This is assumed -- not in paper
    k_lookup(6,4) = 295.3322; r0_lookup(6,4) = 1.0000

    !Angular interactions
    angular_k_lookup = 0.d0;          angular_r0_lookup = 0.d0
    angular_k_lookup(6,6,6) = 4.3196; angular_r0_lookup(6,6,6) = 2.75064
    angular_k_lookup(6,6,7) = 4.3196; angular_r0_lookup(6,6,7) = 2.75064
    angular_k_lookup(6,7,7) = 4.3196; angular_r0_lookup(6,7,7) = 2.75064
    angular_k_lookup(6,7,6) = 4.3196; angular_r0_lookup(6,7,6) = 2.75064
    angular_k_lookup(7,6,6) = 4.3196; angular_r0_lookup(7,6,6) = 2.75064
    angular_k_lookup(7,7,6) = 4.3196; angular_r0_lookup(7,7,6) = 2.75064
    angular_k_lookup(7,7,7) = 4.3196; angular_r0_lookup(7,7,7) = 2.75064
    angular_k_lookup(7,6,7) = 4.3196; angular_r0_lookup(7,6,7) = 2.75064

    ! Define anything that isn't already defined
    ! Epsilon and lambda cross rules from:
    ! Thomas Lafitte, Anastasia Apostolakou, Carlos Avendaño, 
    ! Amparo Galindo, Claire S. Adjiman, Erich A. Müller
    ! and George Jackson (2013) "Accurate statistical 
    ! associating fluid theory for chain molecules formed from Mie segments"
    ! The Journal of Chemical Physics 139, 154504 ; doi: 10.1063/1.4819786
    do i = 1,ntypes
    do j = 1,ntypes

        if (lambdar_lookup(i,j) .lt. 1e-5) then
            lambdar_lookup(i,j) = sqrt( (lambdar_lookup(i,i)-3.d0) &
                                       *(lambdar_lookup(j,j)-3.d0)) + 3.d0
        endif

        if (lambdaa_lookup(i,j) .lt. 1e-5) then
            lambdaa_lookup(i,j) = sqrt( (lambdaa_lookup(i,i)-3.d0) &
                                       *(lambdaa_lookup(j,j)-3.d0)) + 3.d0
        endif

        ! Replace undefined sigma from cross rules using Lorentz-Berthelot
        if (sigma_lookup(i,j) .lt. 1e-5) then
            sigma_lookup(i,j) = 0.5d0*(sigma_lookup(i,i)+sigma_lookup(j,j))
        endif
        if (epsilon_lookup(i,j) .lt. 1e-5) then
            epsilon_lookup(i,j) = (sqrt((sigma_lookup(i,i)**3)   & 
                                       *(sigma_lookup(j,j)**3))  &
                                       /(sigma_lookup(i,j)**3))  & 
                                  *sqrt(epsilon_lookup(i,i)      &
                                       *epsilon_lookup(j,j))
        endif

        !Get C_array
        C_lookup(i,j) = (lambdar_lookup(i,j)/(lambdar_lookup(i,j)-lambdaa_lookup(i,j))) & 
                        *(lambdar_lookup(i,j)/lambdaa_lookup(i,j))**(lambdaa_lookup(i,j) & 
                        /(lambdar_lookup(i,j)-lambdaa_lookup(i,j)))

        potshift_lookup(i,j) = get_energy(1.d0/rcutoff**2,i,j)

    enddo
    enddo

    ! Enforce Symmetry of lookup tables
    do i = 1,ntypes
    do j = i,ntypes
        !if (lambdar_lookup(i,j) .lt. 1e-5) then
            lambdar_lookup(i,j) = lambdar_lookup(j,i)
        !endif

        !if (lambdaa_lookup(i,j) .lt. 1e-5) then
            lambdaa_lookup(i,j) = lambdaa_lookup(j,i)
        !endif
        !if (sigma_lookup(i,j) .lt. 1e-5) then
            sigma_lookup(i,j) = sigma_lookup(j,i)
        !endif
        !if (epsilon_lookup(i,j) .lt. 1e-5) then
            epsilon_lookup(i,j) = epsilon_lookup(j,i)
        !endif
        !if (C_lookup(i,j) .lt. 1e-5) then
            C_lookup(i,j) = C_lookup(j,i)
        !endif
        !if (potshift_lookup(i,j) .lt. 1e-5) then
            potshift_lookup(i,j) = potshift_lookup(j,i)
        !endif

            k_lookup(i,j)  = k_lookup(j,i)
            r0_lookup(i,j) = r0_lookup(j,i)

        !How do we enforce symmetry on a rank 3 tensor?!?
        !Done manually above as only 2 molecules have interactions
!        do k = i,ntypes
!            angular_k_lookup(i,j,k)  = angular_k_lookup(k,i,j)
!            angular_r0_lookup(i,j,k) = angular_k_lookup(k,i,j)
!        enddo

    enddo
    enddo

    call get_equilibrium_seperations()

    !Print all properties
    if (irank .eq. iroot) then
        do i = 1,ntypes
        do j = 1,ntypes
            print'(a,2i6,2a5,f10.5,3a4,6f10.5)', 'SAFT PARAMETERS',i,j, & 
                                        moltype_names(i), ' m= ', mass_lookup(i), ' to ', & 
                                        moltype_names(j), ' m= ', mass_lookup(j), &
                                        sigma_lookup(i,j), & 
                                        epsilon_lookup(i,j), &
                                        lambdar_lookup(i,j), &
                                        lambdaa_lookup(i,j), &
                                        equil_sep_lookup(i,j)
        enddo
        enddo
    endif

    contains

        !Calculate the equilbrium seperations -- assumes moltype has not been allocated yet

        subroutine get_equilibrium_seperations
            implicit none

            integer          :: i,j,n
            real(kind(0.d0)) :: rand, accijmag,error, step, val, rij2, invrij2
            real(kind(0.d0)) :: temp(ntypes)

            ! Get equilibrium separation of Mie and harmonic potential
            ! Could probably solve this exactly or using a crompulent
            ! numerical technique but I can't be bothered so I'm 
            ! using a Monte Carlo/steepest descent type method...


            ! Save first ntypes moltypes in case they've already been defined.
            do n =1,ntypes
                temp(n) = moltype(n)
                moltype(n) = n
            enddo

            !Loop through all possible interactions
            do i = 1,ntypes
            do j = 1,ntypes
                !Get initial guess
                n = 0
                val = 100d0
                error = val
                equil_sep_lookup(i,j) = 1.d0
                accijmag = 1.d0
                call random_number(rand)
                rij2 = 10.d0*(rand-0.5d0)
                !Steepest descent type thing until error in force is small
                do while (error .gt. 1e-12) 
                    call random_number(rand)
                    step = min(0.1d0,error)
                    rij2 = equil_sep_lookup(i,j) + sign(step*rand,accijmag)
                    invrij2 = 1.d0 / rij2
                    accijmag = harmonic_accijmag(rij2, i, j) + Mie_accijmag(invrij2, i, j)
                    if ( abs(accijmag) .lt. abs(val)) then
                        equil_sep_lookup(i,j) = rij2
                        error = abs(val)-abs(accijmag)
                        !print'(i7,2i3,4(a,f20.10))',n, i,j, ' Error = ',  error , ' Previous =', val, ' new = ',accijmag , ' Seperation = ', sqrt(rij2)
                        val = min(abs(val),abs(accijmag))
                    endif

                    !Exit conditions if failure or non-forces
                    if (abs(accijmag) .lt. 1e-5) exit ! print'(3i6,2f10.5)', n,i,j, accijmag, val
                    n = n + 1
                    if (n .gt. 1000000) then
                        print*, 'failed to find equilibrium distance'
                        print'(i7,2i3,4(a,f20.10))',n, i,j, ' Error = ',  &
                        error , ' Previous =', val, ' new = ',accijmag ,  &
                        ' Separation = ', sqrt(rij2)
                        exit                        
                    endif

                enddo
                equil_sep_lookup(i,j) = sqrt(equil_sep_lookup(i,j))
            enddo
            enddo

            !Sanity check
            do i = 1,ntypes
            do j = 1,ntypes
                rij2 = equil_sep_lookup(i,j)**2
                invrij2 = 1.d0 / rij2
                accijmag = Mie_accijmag(invrij2, i, j) + harmonic_accijmag(rij2, i, j)
                if (accijmag .gt. 1e-5) print*, 'WARNING -- equilibrium force not zero for', i,j, accijmag
            enddo
            enddo
            
            !Copy first ntypes moltypes back to array.
            do n =1,ntypes
                moltype(n) = temp(n)
            enddo


        end subroutine get_equilibrium_seperations



end subroutine setup_mie_potential



end module module_set_parameters 
!------------------------------------------------------------------------------

subroutine setup_set_parameters
	use module_set_parameters
	use boundary_MD, only: bforce_flag
	use interfaces, only : error_abort
	use librarymod, only : build_hilbert
	implicit none

	integer							:: i,iblk,jblk,kblk
	integer,dimension(3)			:: nblocks,nvalues
	real(kind(0.d0)),dimension(3)	:: blocksidelength

	!This has alreay been done in initialise for a coupled run
#if (USE_COUPLER == 0)	
   	call set_parameters_global_domain
	call set_parameters_cells
#endif
	!call set_parameters_setlimits

	!Allocate array sizes for position, velocity and acceleration
	call set_parameters_allocate

	!Zero quantities and arrays
	r = 0.d0
	v = 0.d0
	a = 0.d0

	zeta= 0.d0	!Set Nose Hoover thermostat scaling property to zero
	halo_np = 0

	call set_parameters_outputs
	call setup_linklist
	call establish_surface_bins
	call establish_halo_bins
	call establish_surface_cells
	call establish_halo_cells

	nvalues = (/ min(ncells(1),nbins(1)), & 
				 min(ncells(2),nbins(2)), & 
				 min(ncells(3),nbins(3)) /)
	call establish_halo_cellbins(nvalues)

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
			F_ext_limits(i) =  0.5d0*globaldomain(ceiling(real(i)/2.d0))
		elseif (F_ext_limits(i) .lt. -0.5d0*globaldomain(ceiling(real(i)/2.d0))) then
			F_ext_limits(i) = -0.5d0*globaldomain(ceiling(real(i)/2.d0))
		endif
	enddo

    !Choose potenital and force calculation method to use
    if ( Mie_potential  .eq. 0) then
        mass => LJ_mass
        get_accijmag => LJ_accijmag
        get_force => LJ_force
        get_energy => LJ_energy
        get_poly_accijmag => FENE_accijmag
        get_poly_force => FENE_force
        get_poly_energy => FENE_energy
    elseif (Mie_potential  .eq. 1) then
        mass => Mie_mass
        get_accijmag => Mie_accijmag
        get_force => Mie_force
        get_energy => Mie_energy
        get_poly_accijmag => harmonic_accijmag
        get_poly_force => harmonic_force
        get_poly_energy => harmonic_energy
    else
        call error_abort("Error in simulation_compute_forces -- Mie potential flag is incorrectly specified")
    end if

    !Setup Mie potential if LJ not used
    if (mie_potential .eq. 0) then
    	!Calculate shift in lennard-Jones potential based on cutoff
        potshift = 4.d0*(1.d0/rcutoff**12 - 1.d0/rcutoff**6)
    elseif (mie_potential .eq. 1) then
        call setup_mie_potential
    endif

	!Calculate correction to lennard-Jones potential/pressure based on cutoff
	if (sLRC_flag .ne. 0) then
		potential_sLRC = 8.d0*pi*density	  *(1.d0/(9.d0*rcutoff**9) - 1.d0/(3.d0*rcutoff**3))
		Pressure_sLRC  = 8.d0*pi*density**2.d0*(4.d0/(9.d0*rcutoff**9) - 2.d0/(3.d0*rcutoff**3))
	else
		potential_sLRC = 0.d0; Pressure_sLRC = 0.d0;
	endif

    !If Mie potential used, standard LRC are no longer valid
    if (Mie_potential .ne. 0) then
        potential_sLRC = 0.d0; Pressure_sLRC = 0.d0;
    endif

end subroutine setup_set_parameters


!===========================================================================================
!Allocate arrays based on number of particles, np, using extra allocation

subroutine set_parameters_allocate
	use module_set_parameters
	use shear_info_MD
	implicit none

	integer :: ixyz, n
    double precision    :: temp

	!Calculate required extra allocation of molecules to allow copied Halos
	!using ratio of halo to domain volume (with safety factor)
	extralloc = 0
	do ixyz =1,nd
		extralloc = extralloc + &
		ceiling(((2.d0*(3*ncells(ixyz)**2+ &
		                6*ncells(ixyz)+4)) &
                         /ncells(ixyz)**3))*np
	enddo
	!extralloc = extralloc/nd  + 300  !Average of all 3 dimensions inc safety factor
    !Set to 2000 here as start case was 2 phase with uneven distribution
	extralloc = extralloc/nd  + 2000 

	!Allocate array sizes for position, velocity and acceleration
	allocate(r(nd,np+extralloc))
	allocate(v(nd,np+extralloc))
	allocate(a(nd,np+extralloc))
!	if (eflux_outflag .eq. 4) then
!		allocate(a_old(nd,np+extralloc))
!	endif

	!Allocate potential energy and virial per molecule array
	allocate(potenergymol(np+extralloc))
	!allocate(potenergymol_mdt(np+extralloc))
	allocate(potenergymol_LJ(np+extralloc))
	allocate(virialmol(np+extralloc))

	!Check if rtrue required
	if (r_gyration_outflag .eq. 1 .or. vmd_outflag	   .eq. 4) rtrue_flag = 1
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

    !If molecular types are used
    if (Mie_potential .eq. 1) then
        allocate(moltype(np+extralloc)) 
        !Default value is 2 (models 2 x water per bead with Mie)
        moltype(1:np) = default_moltype

    endif

	!Allocate arrays use to fix molecules and allow sliding
	allocate(tag(np+extralloc)); tag = free
	if (ensemble .eq. tag_move) then
		allocate(fix(nd,np+extralloc)); fix = 1	!default - set fix to one (unfixed)
		allocate(rtether(nd,np+extralloc))
		allocate(slidev(nd,np+extralloc))
	endif

    !If necessary, allocate global molecular number
    if (global_numbering .ne. 0) then
        allocate(glob_no(np+extralloc))
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
	allocate(potenergymol_POLY(np+extralloc))

	etevtcf = 0.d0
	R_g	 = 0.d0
	
	if (iter .gt. etevtcf_iter0)	etevtcf_iter0	= iter
	if (iter .gt. r_gyration_iter0) r_gyration_iter0 = iter

	!intbits = bit_size(monomer(1)%bin_bflag(1))

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
	use boundary_MD, only: specular_flag, specular_flat, specular_wall, &
						   specular_radial
	use interfaces, only: error_abort
	implicit none

	integer				:: ixyz, extranp

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
			globalnp=1	  !Set number of particles to unity for loop below
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
		
		case ('dense_fene','fene_solution','single_fene')
			
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

		case('solid_liquid','polymer_brush', & 
             'droplet2D','droplet3D','2phase', & 
             '2phase_surfactant_solution', & 
             '2phase_surfactant_atsurface', &
              '2phase_LJ')

			volume=1	!Set domain size to unity for loop below
			do ixyz=1,nd
				globaldomain(ixyz) = initialnunits(ixyz) & 	!Size domain based on required density
				/((density/4.d0)**(1.d0/nd))
				volume = volume*globaldomain(ixyz)		!Volume based on size of domain
			enddo

			! no need to fix globalnp if we have it already
			if(.not. restart) then
				globalnp=1	  !Set number of particles to unity for loop below
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
		
		case ('fill_cylinders','fill_cylinders_fene_solution')

			!Read globalnp, globaldomain and r_oo etc
			call parallel_io_cyl_footer(cyl_file)


			volume = product(globaldomain(1:3))
			! Get closest number of FCC units in global domain from input
			initialnunits(:) = ceiling(globaldomain(:)*((density/4.d0)**(1.d0/nd)))
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

			call guess_block_cylinder_np()
			!np = globalnp / nproc

			! turn off periodicity in x and y direction
			periodic(1) = 0
			periodic(2) = 0
			if (specular_flag .eq. specular_flat .and. &
				specular_wall(3) .ne. 0.0) then
				periodic(3) = 0
			else
				periodic(3) = 1
			end if

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
			call guess_block_cylinder_np()
			!np = globalnp / nproc

			! turn off periodicity in x and y direction
			periodic(1) = 0
			periodic(2) = 0
			if (specular_flag .eq. specular_flat .and. &
				specular_wall(3) .ne. 0.0) then
				periodic(3) = 0
			else
				periodic(3) = 1
			end if

			!print('(a,3f12.3,a,3f12.3,a,i10,a,i4,a,3l4,a,f12.3,a,f12.3)'), &
			!'globaldomain',globaldomain,' domain',domain,'np', np, 'nproc',&
			! nproc, 'periodic', periodic, 'density', density, 'volume', volume

			! Calculate omega ramping parameters
			omega = omega_i
			omega_rampiters = nint(omega_ramplength/delta_t)
			omega_ramplength = real(omega_rampiters,kind(0.d0))*delta_t

		case default

			call error_abort("Error in set_parameters_global_domain -- must be corrected for this special case")	

		end select

	case default

		call error_abort("Unrecognised initial_config_flag in set_parameters_global_domain")

	end select

contains

	subroutine guess_block_cylinder_np()
		use circle_rectangle_intersection
		use messenger_data_exchange, only : PlaneSum
		implicit none

		real(kind(0.d0)) :: A, A_o, A_i, Atotal, ratio
		real(kind(0.d0)) :: vx(4), vy(4) 

		call block_vertices(iblock,jblock,vx,vy) 
		call circle_rectangle_intersection_area(vx, vy, r_oo, npx, npy, A_o)
		call circle_rectangle_intersection_area(vx, vy, r_ii, npx, npy, A_i)
		A = A_o - A_i
	   
		Atotal = A 
		call PlaneSum(Atotal, 3)

		ratio = A / Atotal
		np = globalnp * ratio 

		!print*, irank, np
		
		!call messenger_syncall
		!stop

	end subroutine guess_block_cylinder_np

	subroutine block_vertices(iblock, jblock, vx, vy)
		implicit none

		integer, intent(in) :: iblock, jblock
		real(kind(0.d0)), intent(out) :: vx(4), vy(4)

		real(kind(0.d0)) :: ivec(2), jvec(2), vertices(4,2)

		ivec = (/domain(1), 0.d0/)
		jvec = (/0.d0, domain(2)/)

		!	2		   
		!	x -------- x 3
		!	|		   |
		!	|		   |
		!	|		   |
		!	|		   |
		! 1 x -------- x
		!			   4  

		vertices(1,:) = (iblock-1)*ivec + (jblock-1)*jvec
		vertices(2,:) = vertices(1,:) + jvec 
		vertices(3,:) = vertices(1,:) + ivec + jvec
		vertices(4,:) = vertices(1,:) + ivec

		vx(:) = vertices(:,1) - globaldomain(1)/2.0
		vy(:) = vertices(:,2) - globaldomain(2)/2.0

	end subroutine block_vertices


end subroutine set_parameters_global_domain


!-----------------------------------------------------------------------------------------

subroutine set_parameters_cells
	use interfaces
	use module_set_parameters
	use polymer_info_MD
	implicit none

	integer :: ixyz
	real(kind(0.d0)) :: rneighbr

	!Calculate size of neighbour list region
	rneighbr  = rcutoff + delta_rneighbr
	rneighbr2 = rneighbr**2

	select case(potential_flag)
	case(1)
		select case(solvent_flag)
		case(0)
			if (rneighbr < R_0) then
				rneighbr = R_0 
				rneighbr2 = R_0**2
				if(irank.eq.iroot) print*, 'Neighbour list distance rneighbr set to &
						& maximum elongation of polymer spring, ',R_0
			end if
		case(1)
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
		print*, ncells(1),'	in x and ', ncells(2), '	in y' , ncells(3), '	in z' 
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
!Setup all paramters for all output parameters

subroutine set_parameters_outputs
	use module_set_parameters
	use interfaces
	use boundary_MD, only: bforce_pdf_measure, bforce_pdf, bforce_pdf_nsubcells, &
						   bforce_pdf_nbins, bforce_pdf_min, bforce_pdf_max
	use CV_objects, only : CVcheck_mass,CVcheck_momentum, & 
						   CV_constraint,CVcheck_energy, CV_debug!,CV_sphere_momentum,CV_sphere_mass
    use messenger, only : localise_bin
	implicit none

	integer					:: n,i,j,k, ixyz
	real(kind(0.d0))		:: shift


	!Use definition of temperature and re-arrange to define an average velocity minus 3 degrees of
	!freedom - this is to fix the momentum of the domain boundaries 
	initialvel = sqrt(nd * (1.d0 - 1.d0/np)*inputtemperature)

	!Allocate bins used for calculating simulation properties
	gnbins(1) = binspercell(1)*npx*ncells(1) !Total number of domain bins
 	gnbins(2) = binspercell(2)*npy*ncells(2) !Total number of domain bins
	gnbins(3) = binspercell(3)*npz*ncells(3) !Total number of domain bins

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
    binsize = domain/nbins

	!Velocity PDF binning routines
	select case(vPDF_flag)
	case(1:4)
		velPDFMB = PDF(NPDFbins,-PDFvlims,PDFvlims)
		allocate(velPDF_array(nbinso(1),nbinso(2),nbinso(3),nd))
		do i = 1,nbinso(1)
		do j = 1,nbinso(2)
		do k = 1,nbinso(3)
		do ixyz = 1,nd
			velPDF_array(i,j,k,ixyz) = PDF(NPDFbins,-PDFvlims,PDFvlims)
		enddo
		enddo
		enddo
		enddo
	case(5)
		!Instantiate whole domain PDF object
		velPDF   = PDF(NPDFbins,-PDFvlims,PDFvlims)
		velPDFMB = PDF(NPDFbins,-PDFvlims,PDFvlims)
	end select

	if (bforce_pdf_measure.ne.0) then
		allocate(bforce_pdf(nd,bforce_pdf_nsubcells))
		do j=1,bforce_pdf_nsubcells
			do i = 1,nd
				bforce_pdf(i,j) = PDF(bforce_pdf_nbins, bforce_pdf_min, bforce_pdf_max)
			end do
		end do
	end if

	!Allocate and define number of shells used for Radial distribution function (rdf)
	if (rdf_outflag .eq. 1) then
		allocate(rdf(rdf_nbins))						!Allocate array for radial distribution function
		allocate(rdf_hist(rdf_nbins))				   !Allocate array to tally positions
		rdf= 0.d0
		rdf_hist= 0
	elseif(rdf_outflag .eq. 2) then
		allocate(rdf3d(rdf_nbins,nd))				   !Allocate array for radial distribution function
		allocate(rdf3d_hist(rdf_nbins,nd))			  !Allocate array to tally positions
		rdf3d_hist= 0
		rdf3d= 0.d0
	endif
	!Allocate and define arrays for static structure factor
	if (ssf_outflag .eq. 1) then
		allocate(ssf_hist(2*ssf_nmax+1,2*ssf_nmax+1))   !Allocate array to tally positions
		allocate(ssf(2*ssf_nmax+1,2*ssf_nmax+1))		!Allocate array for radial distribution function
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

	!Allocated arrays for momentum slice
	if (momentum_outflag.ne.0 .and. momentum_outflag.lt.4) then
		allocate(slice_momentum(nbins(momentum_outflag),3))
		allocate(slice_mass(nbins(momentum_outflag)))
		slice_momentum = 0.d0
		slice_mass = 0
	else
		if (mass_outflag.ne.0 .and. mass_outflag.lt.4) then
			allocate(slice_mass(nbins(mass_outflag)))
			slice_mass = 0
		endif
	endif

	!Allocated bins for momentum averaging
	if (momentum_outflag.eq.4) then
		if (.not. allocated(volume_momentum)) then
			allocate(volume_momentum(nbinso(1),nbinso(2),nbinso(3),3))
		endif
        volume_momentum = 0.d0
        mass_outflag = 4	!Mass binning required too
        select case (split_pol_sol_stats)
        case(0)
            !Do nothing extra
        case(1)
		    if (.not. allocated(volume_momentum_s)) then
                allocate(volume_momentum_s(nbinso(1),nbinso(2),nbinso(3),3))
            end if
		    if (.not. allocated(volume_momentum_p)) then
                allocate(volume_momentum_p(nbinso(1),nbinso(2),nbinso(3),3))
            end if
			volume_momentum_s = 0.d0; volume_momentum_p = 0.d0
        case default
            call error_abort('Invalid potential flag in volume_momentum allocate')
        end select
	endif

	!Allocated bins for temperature averaging
	if (temperature_outflag.eq.4) then
		allocate(volume_temperature(nbinso(1),nbinso(2),nbinso(3)))
		volume_temperature = 0.d0
		mass_outflag = 4	!Mass binning required too
		!Allocate and zero peculiar momentum binning array
		if (peculiar_flag .ne. 0) then
			allocate(u(nd,np+extralloc)); u = 0.d0
			if (momentum_outflag.ne.4) then
				call error_abort("set_parameters_outputs Error -- Temperature outflag on with &
								 &perculiar momentum but momentum binning is off. Please switch &
								 &MOMENTUM_OUTFLAG 1st option to 4 or TEMPERATURE_OUTFLAG 3rd option to 0")
			endif
		endif
	endif

	!Allocated bins for energy averaging
	if (energy_outflag.eq.4) then
		allocate(volume_energy(nbinso(1),nbinso(2),nbinso(3)))
		volume_energy = 0.d0
	endif


	!Allocate mass bins if they haven't been already allocated (and they are needed)
	if (mass_outflag.eq.4) then
		if (.not. allocated(volume_mass)) then
			allocate(volume_mass(nbinso(1),nbinso(2),nbinso(3)))
			volume_mass = 0
		endif
        select case (split_pol_sol_stats)
        case(0)
            !Do nothing extra
        case(1)
		    if (.not. allocated(volume_mass_s)) then
                allocate(volume_mass_s(nbinso(1),nbinso(2),nbinso(3)))
            end if
		    if (.not. allocated(volume_mass_p)) then
                allocate(volume_mass_p(nbinso(1),nbinso(2),nbinso(3)))
            end if
			volume_mass_s = 0; volume_mass_p = 0
        case default
            call error_abort('Invalid potential flag in volume_mass allocate')
        end select
	endif

	! Allocate cylindrical polar bins

	if (config_special_case .eq. 'rotate_cylinders' .or. &
		config_special_case .eq. 'fill_cylinders' .or. &
		config_special_case .eq. 'fill_cylinders_fene_solution' &
		) then

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

		! Allocate bins for field io
		if (temperature_outflag .eq. 5) then

			allocate(cyl_KE(cpol_binso(1),cpol_binso(2),cpol_binso(3)))
			cyl_KE = 0.d0

			if ( .not. allocated(cyl_mass) ) then

				allocate(cyl_mass(cpol_binso(1),cpol_binso(2),cpol_binso(3)))
				cyl_mass = 0

			end if

		endif

		if ( momentum_outflag .eq. 5 ) then

			allocate(cyl_mom(cpol_binso(1),cpol_binso(2),cpol_binso(3),3))
			cyl_mom  = 0.d0

			if ( .not. allocated(cyl_mass) ) then
				allocate(cyl_mass(cpol_binso(1),cpol_binso(2),cpol_binso(3)))
				cyl_mass = 0
			end if

		end if

		if ( mass_outflag .eq. 5 ) then

			if ( .not. allocated(cyl_mass) ) then

				allocate(cyl_mass(cpol_binso(1),cpol_binso(2),cpol_binso(3)))
				cyl_mass = 0

			end if

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
		allocate( rfbin( nbinso(1), nbinso(2), nbinso(3),3,3  ))
		allocate( vvbin( nbins (1),  nbins(2),  nbins(3),3,3  ))
		allocate( Pxybin(nbins (1),  nbins(2),  nbins(3),3,3  ))
		rfbin  = 0.d0
		vvbin = 0.d0
		Pxybin = 0.d0
	else if (pressure_outflag .eq. 3) then
		!Allocate pressure bin for Stress volume averaging in polar coords
		allocate( rfbin(  cpol_binso(1), cpol_binso(2), cpol_binso(3) ,3,3 ))
		allocate( vvbin(  cpol_binso(1), cpol_binso(2), cpol_binso(3) ,3,3 ))
		allocate( Pxybin( cpol_binso(1), cpol_binso(2), cpol_binso(3) ,3,3 ))
		rfbin  = 0.d0
		vvbin = 0.d0
		Pxybin = 0.d0
	endif

    if (heatflux_outflag .eq. 2) then
		!Allocate pressure bin for Stress volume averaging
		allocate( rfvbin( nbinso(1), nbinso(2), nbinso(3), 3, 1 ))
		allocate( evbin( nbins (1),  nbins(2),  nbins(3),3  ))
		rfvbin  = 0.d0
		evbin = 0.d0
    elseif (heatflux_outflag .eq. 0) then
        !pass
    else
        call error_abort("Heatflux_outflag only developed for case 2")
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
			if (external_force_flag .ne. 0 .or. &
				ensemble .eq. tag_move .or.	 &
				CVforce_flag .ne. VOID) then
				allocate(F_ext_bin(nbinso(1),nbinso(2),nbinso(3),3))
				F_ext_bin = 0.d0
			endif
			momentum_flux 	= 0.d0
			volume_momentum = 0.d0
			volume_force 	= 0.d0
			if (CV_debug .eq. 1) then
				call CVcheck_mass%initialise(nbins,nhb,domain,delta_t,Nmflux_ave)   ! initialize CVcheck
				call CVcheck_momentum%initialise(nbins,nhb,domain,delta_t,Nvflux_ave)   ! initialize CVcheck
				call CV_constraint%initialise(nbins,nhb,domain,delta_t,Nvflux_ave)   ! initialize CV constraint object
				call CVcheck_energy%initialise(nbins,nhb,domain,delta_t,Neflux_ave)   ! initialize CVcheck
            elseif (CV_debug .eq. 2) then
                print*, 'bin', debug_CV,localise_bin(debug_CV)
				call CVcheck_mass%initialise(nbins,nhb,domain, & 
                                             delta_t,Nmflux_ave, & 
                                             localise_bin(debug_CV))     ! initialize CVcheck
				call CVcheck_momentum%initialise(nbins,nhb,domain, & 
                                                 delta_t,Nvflux_ave, & 
                                                 localise_bin(debug_CV)) ! initialize CVcheck
				call CV_constraint%initialise(nbins,nhb,domain,delta_t,Nvflux_ave)   ! initialize CV constraint object
				call CVcheck_energy%initialise(nbins,nhb,domain, & 
                                               delta_t,Neflux_ave, & 
                                               localise_bin(debug_CV))   ! initialize CVcheck
            
				!call CV_sphere_mass%initialise((/1,1,1/))	
				!call CV_sphere_momentum%initialise_sphere((/1,1,1/),collect_spherical=.false.)	
			endif
			!Allocate bins for control volume mass fluxes
			if (.not.(allocated(volume_mass)))  allocate(volume_mass(nbinso(1),nbinso(2),nbinso(3)))
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
				if (CV_debug .eq. 1) then
					call CVcheck_mass%initialise(nbins,nhb,domain,delta_t,Nmflux_ave)   ! initialize CVcheck
                elseif (CV_debug .eq. 2) then
					call CVcheck_mass%initialise(nbins,nhb,domain, & 
                                                delta_t,Nmflux_ave, & 
                                                localise_bin(debug_CV))   ! initialize CVcheck
				endif
			endif
	end select

	!Allocate bins for control volume energy fluxes and forces*velocity
	if (eflux_outflag .eq. 4) then
		allocate(  energy_flux(nbinso(1),nbinso(2),nbinso(3),6))
		allocate( Pxyvface(nbinso(1),nbinso(2),nbinso(3),6))
		allocate( Pxyvface_mdt(nbinso(1),nbinso(2),nbinso(3),6))
		allocate( Pxyvface_integrated(nbinso(1),nbinso(2),nbinso(3),6))
		energy_flux 	= 0.d0; Pxyvface = 0.d0; 
		Pxyvface_mdt=0.d0; Pxyvface_integrated = 0.d0
		if (external_force_flag .ne. 0 .or. & 
			ensemble .eq. tag_move .or.	 & 
			CVforce_flag .ne. VOID) then
			allocate(Fv_ext_bin(nbinso(1),nbinso(2),nbinso(3)))
			Fv_ext_bin = 0.d0
		endif
	elseif (eflux_outflag .eq. 1 .or. eflux_outflag .eq. 2 .or. eflux_outflag .eq. 3) then
		stop "Error - eflux MOP is not coded!"
	endif


    if (msurf_outflag .eq. 1) then
        allocate(surface_density(nbinso(1),nbinso(2),nbinso(3),6))
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
! Establish and store indices of cells/bins which are in the halo
! for when the bin2cell ratio is smaller/larger than unity

subroutine establish_halo_cellbins(nvalues)
	use module_set_parameters
	implicit none

	integer, intent(in)	:: nvalues(3)

	integer		:: n
	integer		:: i, j, k

	nhalocellbins  =	2*((nvalues(1)+2)*(nvalues(2)+2) &
						+  (nvalues(3)  )*(nvalues(2)+2) &
						+  (nvalues(3)  )*(nvalues(1)  ))

	allocate(halocellbins(nhalocellbins,3))

	n = 1
	do k=1, nvalues(3)+2
	do j=1, nvalues(2)+2
	do i=1, nvalues(1)+2

		!Remove inner part of domain
		if((i .gt. (1) .and. i .lt. (nvalues(1)+2)) .and. &
		   (j .gt. (1) .and. j .lt. (nvalues(2)+2)) .and. &
		   (k .gt. (1) .and. k .lt. (nvalues(3)+2))) cycle

		halocellbins(n,1)=i
		halocellbins(n,2)=j
		halocellbins(n,3)=k
		n = n + 1

	enddo
	enddo
	enddo

end subroutine establish_halo_cellbins

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
!					+  (ncells(2)-2-2*buf)*(ncells(3)-2-2*buf))

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
