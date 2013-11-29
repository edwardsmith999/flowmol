
!===================================================================================
!
!	Routines add forces and fluxes per molecule to running totals
!
!===================================================================================



!Fortran object oriented solution to CV conservation checking
module CV_objects
	implicit none

type :: check_CV_mass
	integer,dimension(:,:,:,:),allocatable 	:: flux
	integer,dimension(:,:,:),allocatable	:: dXdt, X, X_minus_t, X_minus_2t
	contains
		procedure :: initialise  => initialise_mass
		procedure :: update_dXdt => update_dXdt_mass
		procedure :: check_error => check_error_mass
		procedure :: swap_halos  => swap_halos_mass
end type check_CV_mass

type :: check_CV_momentum
	double precision,dimension(:,:,:,:,:),allocatable 	:: flux, Pxy,  Pxy_minus_t
	double precision,dimension(:,:,:,:),allocatable		:: dXdt, X, X_minus_t, X_minus_2t, & 
														   F_ext, totalflux, totalpressure
	contains
		procedure :: initialise  => initialise_momentum
		procedure :: update_dXdt => update_dXdt_momentum
		procedure :: update_flux => update_flux_momentum
		procedure :: update_Pxy  => update_Pxy
		procedure :: update_F_ext=> update_F_ext
		procedure :: check_error => check_error_momentum
		procedure :: swap_halos  => swap_halos_momentum
end type check_CV_momentum

type :: check_CV_energy
	double precision,dimension(:,:,:,:),allocatable 	:: flux, Pxyv,  Pxyv_minus_t
	double precision,dimension(:,:,:),allocatable		:: dXdt, X, X_minus_t, X_minus_2t, & 
														   Fv_ext, totalflux, totalpower
	contains
		procedure :: initialise  => initialise_energy
		procedure :: update_dXdt => update_dXdt_energy
		procedure :: update_flux => update_flux_energy
		procedure :: update_Pxy  => update_Pxyv
		procedure :: update_F_ext=> update_Fv_ext
		procedure :: check_error => check_error_energy
		procedure :: swap_halos  => swap_halos_energy

end type check_CV_energy

type, extends(check_CV_mass) :: sphereObj_mass

    double precision                        :: radius = 2.d0
    integer,dimension(:,:,:),allocatable	:: Xtemp
    integer,dimension(:,:,:,:),allocatable	:: fluxTemp
contains
	procedure :: Add_spherical_CV_fluxes => Add_spherical_CV_mass_fluxes
    procedure :: Add_spherical_CV_mass   => Add_spherical_CV_mass
    procedure :: check_error             => check_error_sphere_mass
	procedure :: update_dXdt             => update_dXdt_sphere_mass
				
end type sphereObj_mass

type, extends(check_CV_momentum) :: sphereObj_mom

	logical											:: collect_spherical
    double precision    							:: radius = 2.d0
    double precision,dimension(:,:,:,:),allocatable	:: Xtemp,FsurfaceTemp, Fsurface
contains
	procedure :: initialise_sphere		   => initialise_sphere_momentum
	procedure :: Add_spherical_CV_forces   => Add_spherical_CV_forces
	procedure :: Add_spherical_CV_fluxes   => Add_spherical_CV_fluxes
    procedure :: Add_spherical_CV_velocity => Add_spherical_CV_velocity
    procedure :: check_error               => check_error_sphere_momentum
	!procedure :: update_dXdt =>    update_dXdt_sphere_momentum
				
end type sphereObj_mom

	!Check CV conservation
	type(check_CV_mass)		:: CVcheck_mass							! declare an instance of CV checker
	type(check_CV_momentum)	:: CVcheck_momentum, CVcheck_momentum2	! declare an instance of CV checker
	type(check_CV_energy)	:: CVcheck_energy						! declare an instance of CV checker

    !CV spherical object
    type(sphereObj_mass) :: CV_sphere_mass
    type(sphereObj_mom) :: CV_sphere_momentum


contains

!===================================================
!	   M a s s   C o n t r o l   V o l u m e
!===================================================

	!Constructor for object
	subroutine initialise_mass(self, nb)
		implicit none

		! initialize shape objects
		class(check_CV_mass) 		:: self

		integer, dimension(3),intent(in) :: nb

		allocate(self%flux(nb(1),nb(2),nb(3),6))
		allocate(self%dXdt(nb(1),nb(2),nb(3)))
		allocate(self%X(nb(1),nb(2),nb(3)))
		allocate(self%X_minus_t(nb(1),nb(2),nb(3)))
		allocate(self%X_minus_2t(nb(1),nb(2),nb(3)))

		self%flux 		= 0
		self%dXdt 		= 0
		self%X 			= 0
		self%X_minus_t 	= 0
		self%X_minus_2t = 0

	end subroutine initialise_mass

	!Update time evolution and store previous two values
	subroutine update_dXdt_mass(self, X)
		implicit none

		! initialize shape objects
		class(check_CV_mass) :: self

		integer,dimension(:,:,:),allocatable,intent(in) :: X

		self%X_minus_t  = self%X
		self%X 		    = X

		self%dXdt = self%X - self%X_minus_t

	end subroutine update_dXdt_mass


	!Swap halos on edges of processor boundaries
	subroutine swap_halos_mass(self,nb)
		use messenger_bin_handler, only : swaphalos
		implicit none

		integer							 :: nresults
		integer, dimension(3),intent(in) :: nb

		! initialize shape objects
		class(check_CV_mass) :: self

    	!Include halo surface fluxes to get correct values for all cells
    	nresults = 6
    	call swaphalos(self%flux,nb(1),nb(2),nb(3),nresults)

	end subroutine swap_halos_mass


	!Check error for specified range of bins
	subroutine check_error_mass(self,imin,imax,jmin,jmax,kmin,kmax,iter,irank)
		implicit none

		! initialize shape objects
		class(check_CV_mass) :: self

		integer,intent(in) :: iter,irank,imin,imax,jmin,jmax,kmin,kmax

		integer 		:: i,j,k
		integer,save 	:: first_time = 0
		logical		 	:: check_ok

		!First call doesn't have difference in time yet so skip
		if (first_time .lt. 2) then
			first_time = first_time + 1
			return
		endif

		check_ok = .true.

		do i = imin,imax
		do j = jmin,jmax
		do k = kmin,kmax

			if(sum(self%flux(i,j,k,:))-self%dXdt(i,j,k) .ne. 0) then
				print'(a,i8,4i4,4i8)','Error_cubeCV_mass', iter,irank,i,j,k, & 
					sum(self%flux(i,j,k,:)),self%dXdt(i,j,k),self%X_minus_t(i,j,k),self%X(i,j,k)
				check_ok = .false.
			endif

		enddo
		enddo
		enddo

	end subroutine check_error_mass

!===================================================
!	M o m e n t u m   C o n t r o l   V o l u m e
!===================================================

	!Constructor for object
	subroutine initialise_momentum(self, nb)
		implicit none

		! initialize shape objects
		class(check_CV_momentum) 		:: self

		integer, dimension(3),intent(in) :: nb

		allocate(self%flux(nb(1),nb(2),nb(3),3,6))
		allocate(self%Pxy(nb(1),nb(2),nb(3),3,6))
		allocate(self%Pxy_minus_t(nb(1),nb(2),nb(3),3,6))
		allocate(self%dXdt(nb(1),nb(2),nb(3),3))
		allocate(self%X(nb(1),nb(2),nb(3),3))
		allocate(self%X_minus_t(nb(1),nb(2),nb(3),3))
		allocate(self%F_ext(nb(1),nb(2),nb(3),3))
		allocate(self%totalflux(nb(1),nb(2),nb(3),3))
		allocate(self%totalpressure(nb(1),nb(2),nb(3),3))
		!allocate(self%X_minus_2t(nb(1),nb(2),nb(3),3))

		self%flux 		= 0.d0
		self%Pxy  		= 0.d0
		self%Pxy_minus_t= 0.d0
		self%dXdt 		= 0.d0
		self%X 			= 0.d0
		self%X_minus_t 	= 0.d0
		self%F_ext		= 0.d0
		self%totalflux  = 0.d0
		self%totalpressure = 0.d0
		!self%X_minus_2t = 0.d0

	end subroutine initialise_momentum

	!Update time evolution and store previous two values
	subroutine update_dXdt_momentum(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_momentum) :: self

		double precision,dimension(:,:,:,:),intent(in) :: X

		self%X_minus_t  = self%X
		self%X 		  = X

		self%dXdt = self%X - self%X_minus_t

	end subroutine update_dXdt_momentum

	!Update time evolution and store previous two values
	subroutine update_F_ext(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_momentum) :: self

		double precision,dimension(:,:,:,:),allocatable,intent(in) :: X

		self%F_ext = X

	end subroutine update_F_ext

	!Update time evolution and store previous two values
	subroutine update_flux_momentum(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_momentum) :: self

		double precision,dimension(:,:,:,:,:),allocatable,intent(in) :: X

		self%flux = X

	end subroutine update_flux_momentum

	!Update time evolution and store previous two values
	subroutine update_Pxy(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_momentum) :: self

		double precision,dimension(:,:,:,:,:),allocatable,intent(in) :: X

		self%Pxy_minus_t = self%Pxy
		self%Pxy = X

	end subroutine update_Pxy

	!Swap halos on edges of processor boundaries
	subroutine swap_halos_momentum(self,nb)
		use messenger_bin_handler, only : swaphalos
		implicit none

		integer, dimension(3),intent(in) :: nb

		integer							 :: nresults
		double precision,dimension(:,:,:,:),allocatable :: temp

		! initialize shape objects
		class(check_CV_momentum) :: self

    	! Include halo surface fluxes, stress and external forces 
		! to get correct values for all cells
    	nresults = 18 + 18 + 3
		allocate(temp(nb(1),nb(2),nb(3),nresults))
		temp(:,:,:,1 :18) = reshape(self%flux,(/ nb(1),nb(2),nb(3),18 /))
		temp(:,:,:,19:36) = reshape(self%Pxy ,(/ nb(1),nb(2),nb(3),18 /))
		temp(:,:,:,37:39) = self%F_ext(:,:,:,:)
    	call swaphalos(temp,nb(1),nb(2),nb(3),nresults)
		self%flux = reshape(temp(:,:,:,1 :18),(/ nb(1),nb(2),nb(3),3,6 /))
		self%Pxy  = reshape(temp(:,:,:,19:36),(/ nb(1),nb(2),nb(3),3,6 /))
		self%F_ext = temp(:,:,:,37:39)
		deallocate(temp)

	end subroutine swap_halos_momentum


	!Return the total flux/stresses totals for a given CV
	subroutine return_CV_totals(self)
		use computational_constants_MD, only :delta_t
		implicit none

		! initialize shape objects
		class(check_CV_momentum) :: self

	    !Calculate total CV flux and change in mass
	    self%totalflux(:,:,:,:) =((self%flux(:,:,:,:,1)+self%flux(:,:,:,:,4)) &
	              				 +(self%flux(:,:,:,:,2)+self%flux(:,:,:,:,5)) &
	             				 +(self%flux(:,:,:,:,3)+self%flux(:,:,:,:,6)))/delta_t

	    !Totalpressure = totalpressure*delta_t
	    self%totalpressure(:,:,:,:)=0.25d0 * ( self%Pxy_minus_t(:,:,:,:,1)  & 
											  -self%Pxy_minus_t(:,:,:,:,4)) &
	                  					    +( self%Pxy_minus_t(:,:,:,:,2)  &
											  -self%Pxy_minus_t(:,:,:,:,5)) &
	                  						+( self%Pxy_minus_t(:,:,:,:,3)  &
											  -self%Pxy_minus_t(:,:,:,:,6))

	end subroutine return_CV_totals

	!Check error for specified range of bins
	subroutine check_error_momentum(self,imin,imax,jmin,jmax,kmin,kmax,iter,irank)
		! initialize shape objects
		use computational_constants_MD, only : domain,delta_t,Nvflux_ave
		use calculated_properties_MD, only : nbins
		implicit none

		class(check_CV_momentum) :: self

		integer,intent(in) :: iter,irank,imin,imax,jmin,jmax,kmin,kmax

		logical							:: check_ok
		integer 						:: i,j,k
		integer,save 					:: first_time = 0
		double precision				:: conserved
		double precision,dimension(3)	:: binsize,totalpressure,totalflux,F_ext,dvelocitydt

		!First call doesn't have difference in time yet so skip
		if (first_time .lt. 2) then
			first_time = first_time + 1
			return
		endif

		binsize = domain/nbins

		!print'(2i4,f13.5,3i8)', iter, irank, maxval(self%F_ext), maxloc(self%F_ext)
		!print'(a,4i8,4f13.8)', 'Inside  object', irank,6,6,6, self%F_ext(6,6,6,:),sum(self%F_ext(6,6,6,:))

		check_ok = .true.

		do i = imin,imax
		do j = jmin,jmax
		do k = kmin,kmax


		    !Calculate total CV flux and change in mass
		    totalflux =(self%flux(i,j,k,:,1)+self%flux(i,j,k,:,4))/binsize(1) &
		              +(self%flux(i,j,k,:,2)+self%flux(i,j,k,:,5))/binsize(2) &
		              +(self%flux(i,j,k,:,3)+self%flux(i,j,k,:,6))/binsize(3)

		    !Totalpressure = totalpressure*delta_t
		    totalpressure =(self%Pxy_minus_t(i,j,k,:,1)-self%Pxy_minus_t(i,j,k,:,4))/binsize(1) &
		                  +(self%Pxy_minus_t(i,j,k,:,2)-self%Pxy_minus_t(i,j,k,:,5))/binsize(2) &
		                  +(self%Pxy_minus_t(i,j,k,:,3)-self%Pxy_minus_t(i,j,k,:,6))/binsize(3)
			F_ext = self%F_ext(i,j,k,:)/product(binsize)

			!drhou/dt
		    dvelocitydt =  self%dxdt(i,j,k,:)/(delta_t*Nvflux_ave)

		    !Verify that CV momentum is exactly conservative
		    conserved = sum(totalpressure-totalflux-dvelocitydt-F_ext)
			if(abs(conserved) .gt. 0.000000001d0) then
			!if (i .eq. 3 .and. j .eq. 3 .and. k .eq. 3) then
				print'(a,i8,4i4,7f10.5)','Error_in_momentum_flux', iter,irank,i,j,k, & 
					 conserved, sum(totalpressure),-sum(totalflux),sum(dvelocitydt), & 
					+sum(F_ext), sum(self%X(i,j,k,:)),   & 
					 sum(self%X_minus_t(i,j,k,:))
				check_ok = .false.
			endif

		enddo
		enddo
		enddo

		!if (check_ok .eq. .false.) then
		!	stop "Error in momentum flux"
		!endif

	end subroutine check_error_momentum



!===================================================
!	E n e r g y   C o n t r o l   V o l u m e
!===================================================

	!Constructor for object
	subroutine initialise_energy(self, nb)
		implicit none

		! initialize shape objects
		class(check_CV_energy) 		:: self

		integer, dimension(3),intent(in) :: nb

		allocate(self%flux(nb(1),nb(2),nb(3),6))
		allocate(self%Pxyv(nb(1),nb(2),nb(3),6))
		allocate(self%Pxyv_minus_t(nb(1),nb(2),nb(3),6))
		allocate(self%dXdt(nb(1),nb(2),nb(3)))
		allocate(self%X(nb(1),nb(2),nb(3)))
		allocate(self%X_minus_t(nb(1),nb(2),nb(3)))
		allocate(self%Fv_ext(nb(1),nb(2),nb(3)))
		allocate(self%totalflux(nb(1),nb(2),nb(3)))
		allocate(self%totalpower(nb(1),nb(2),nb(3)))
		!allocate(self%X_minus_2t(nb(1),nb(2),nb(3),3))

		self%flux 		= 0.d0
		self%Pxyv  		= 0.d0
		self%Pxyv_minus_t= 0.d0
		self%dXdt 		= 0.d0
		self%X 			= 0.d0
		self%X_minus_t 	= 0.d0
		self%Fv_ext		= 0.d0
		self%totalflux  = 0.d0
		self%totalpower = 0.d0
		!self%X_minus_2t = 0.d0

	end subroutine initialise_energy

	!Update time evolution and store previous two values
	subroutine update_dXdt_energy(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_energy) :: self

		double precision,dimension(:,:,:),intent(in) :: X

		self%X_minus_t  = self%X
		self%X 		  = X

		self%dXdt = self%X - self%X_minus_t

	end subroutine update_dXdt_energy

	!Update time evolution and store previous two values
	subroutine update_Fv_ext(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_energy) :: self

		double precision,dimension(:,:,:),allocatable,intent(in) :: X

		self%Fv_ext = X

	end subroutine update_Fv_ext

	!Update time evolution and store previous two values
	subroutine update_flux_energy(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_energy) :: self

		double precision,dimension(:,:,:,:),allocatable,intent(in) :: X

		self%flux = X

	end subroutine update_flux_energy

	!Update time evolution and store previous two values
	subroutine update_Pxyv(self, X)
		implicit none
		! initialize shape objects
		class(check_CV_energy) :: self

		double precision,dimension(:,:,:,:),allocatable,intent(in) :: X

		self%Pxyv_minus_t = self%Pxyv
		self%Pxyv = X

	end subroutine update_Pxyv

	!Swap halos on edges of processor boundaries
	subroutine swap_halos_energy(self,nb)
		use messenger_bin_handler, only : swaphalos
		implicit none

		integer, dimension(3),intent(in) 				:: nb

		integer							 				:: nresults
		double precision,dimension(:,:,:,:),allocatable :: temp

		! initialize shape objects
		class(check_CV_energy) :: self

    	! Include halo surface fluxes, stress and external forces 
		! to get correct values for all cells
    	nresults = 6 + 6 + 1
		allocate(temp(nb(1),nb(2),nb(3),nresults))
		temp(:,:,:,1 :6) = self%flux
		temp(:,:,:,7:12) = self%Pxyv
		temp(:,:,:,13  ) = self%Fv_ext
    	call swaphalos(temp,nb(1),nb(2),nb(3),nresults)
		self%flux  = temp(:,:,:,1 :6)
		self%Pxyv   = temp(:,:,:,7:12)
		self%Fv_ext = temp(:,:,:,13)
		deallocate(temp)

	end subroutine swap_halos_energy

	!Check error for specified range of bins
	subroutine check_error_energy(self,imin,imax,jmin,jmax,kmin,kmax,iter,irank)
		! initialize shape objects
		use computational_constants_MD, only : domain,delta_t,Neflux_ave
		use calculated_properties_MD, only : nbins
		implicit none

		class(check_CV_energy) :: self

		integer,intent(in) :: iter,irank,imin,imax,jmin,jmax,kmin,kmax

		logical							:: check_ok
		integer 						:: i,j,k
		integer,save 					:: first_time = 0
		double precision				:: conserved,totalpower,totalflux,Fv_ext,denergydt
		double precision,dimension(3)	:: binsize

		!First call doesn't have difference in time yet so skip
		if (first_time .lt. 2) then
			first_time = first_time + 1
			return
		endif

		binsize = domain/nbins

		!print'(2i4,f13.5,3i8)', iter, irank, maxval(self%Fv_ext), maxloc(self%Fv_ext)
		!print'(a,4i8,4f13.8)', 'Inside  object', irank,6,6,6, self%Fv_ext(6,6,6,:),sum(self%Fv_ext(6,6,6,:))

		check_ok = .true.

		do i = imin,imax
		do j = jmin,jmax
		do k = kmin,kmax

		    !Calculate total CV flux and change in mass
		    totalflux =(self%flux(i,j,k,1)+self%flux(i,j,k,4))/binsize(1) &
		              +(self%flux(i,j,k,2)+self%flux(i,j,k,5))/binsize(2) &
		              +(self%flux(i,j,k,3)+self%flux(i,j,k,6))/binsize(3)

		    !totalpower = totalpower*delta_t
		    totalpower = (self%Pxyv_minus_t(i,j,k,1)-self%Pxyv_minus_t(i,j,k,4))/binsize(1) &
		                +(self%Pxyv_minus_t(i,j,k,2)-self%Pxyv_minus_t(i,j,k,5))/binsize(2) &
		                +(self%Pxyv_minus_t(i,j,k,3)-self%Pxyv_minus_t(i,j,k,6))/binsize(3)
			Fv_ext = self%Fv_ext(i,j,k)/product(binsize)

			!drhou/dt
		    denergydt =  self%dxdt(i,j,k)/(delta_t*Neflux_ave)

		    !Verify that CV momentum is exactly conservative
		    conserved = totalpower-totalflux-denergydt-Fv_ext
			!if(conserved .gt. 0.000000001d0) then
			!if (abs(Fv_ext) .gt. 0.000001) then
			if (i .eq. 3 .and. j .eq. 3 .and. k .eq. 3) then
				print'(a,i8,4i4,9f10.5)','Error_in_energy_flux', iter,irank,i,j,k, & 
					 conserved, totalpower,-totalflux,denergydt, & 
					+Fv_ext, self%X(i,j,k),self%X_minus_t(i,j,k)
				check_ok = .false.
			endif

		enddo
		enddo
		enddo

		!if (check_ok .eq. .false.) then
		!	stop "Error in energy flux"
		!endif

	end subroutine check_error_energy



!===================================================
!	S p h e r i c a l   C o n t r o l   V o l u m e
!===================================================


! - - - - MASS - - - -


    !Check if molecule is inside volume
    subroutine Add_spherical_CV_mass(self,ri)
        use librarymod, only : sphereCV
    	implicit none
    	
        class(sphereObj_mass) :: self

    	real(kind(0.d0)),dimension(3),intent(in)	:: ri

    	!Add for molecule i
    	if (.not. allocated(self%Xtemp)) then
    	    allocate(self%Xtemp(1,1,1))
    	    self%Xtemp = 0
        endif
        
        !if ( sphereCV(ri,self%radius) .gt. 0.00000001) print'(3f10.5,i5,f10.5)', ri, self%Xtemp(1,1,1), sphereCV(ri,self%radius)
    	self%Xtemp(1,1,1) = self%Xtemp(1,1,1) + sphereCV(ri,self%radius)

    end subroutine Add_spherical_CV_mass

    !Update time evolution and store previous two values
    subroutine update_dXdt_sphere_mass(self, X)
    	implicit none

    	! initialize shape objects
    	class(sphereObj_mass) :: self

    	integer,dimension(:,:,:),allocatable,intent(in) :: X

    	self%X_minus_t = self%X
    	self%X 		   = X

    	self%dXdt = self%X - self%X_minus_t

    end subroutine update_dXdt_sphere_mass


    !Check if molecule has crossed surface of CV
    subroutine Add_spherical_CV_mass_fluxes(self,ri,ri_tp1)
        use librarymod, only : sphereCV
    	implicit none
    	
        class(sphereObj_mass) :: self

    	real(kind(0.d0)),dimension(3),intent(in)	:: ri,ri_tp1

    	!Add for molecule i
    	if (.not. allocated(self%fluxTemp)) then
    	    allocate(self%fluxTemp(1,1,1,6))
    	    self%fluxTemp = 0
        endif

    	!Add for molecule i
    	self%fluxTemp(1,1,1,1) = self%fluxTemp(1,1,1,1) + (sphereCV(ri_tp1,self%radius)-sphereCV(ri,self%radius))

    end subroutine Add_spherical_CV_mass_fluxes



    !Check error for specified range of bins
    subroutine check_error_sphere_mass(self,imin,imax,jmin,jmax,kmin,kmax,iter,irank)
    	! initialize shape objects
    	use computational_constants_MD, only : domain,delta_t,Nmflux_ave
    	use calculated_properties_MD, only : nbins
    	implicit none

    	class(sphereObj_mass) :: self

    	integer,intent(in) :: iter,irank,imin,imax,jmin,jmax,kmin,kmax

    	integer 						:: i,j,k
    	integer,save 					:: first_time = 0
    	integer         				:: conserved, dmdt, totalflux
    	double precision,dimension(3)	:: binsize
    	
    	!First call doesn't have difference in time yet so skip
    	if (first_time .lt. 2) then
    		first_time = first_time + 1
    		print'(a)', '                      iter    irank     i       j       k    cnsvd   fluxes   dmdt     m      mt-dt'
    		return
    	endif

    	!Calculate total CV flux and change in mass
    	totalflux =sum(self%flux(1,1,1,:))
        
    	!drhou/dt
        dmdt =  self%dxdt(1,1,1) !(delta_t*Nvflux_ave)

        !Verify that CV momentum is exactly conservative
        conserved = totalflux-dmdt
        
    	if(sum(self%flux(1,1,1,:))-self%dXdt(1,1,1) .ne. 0) then
    		print'(a,11i8)','Error_sphere_mass', iter,irank,1,1,1, & 
    				 conserved, totalflux,dmdt,self%X(1,1,1),   & 
    				 self%X_minus_t(1,1,1)
    	endif

    	self%flux = self%fluxTemp
        self%fluxTemp = 0
        self%Xtemp = 0
        
    end subroutine check_error_sphere_mass




    ! - - - - MOMENTUM - - - -

    !Constructor for sphere object extends momentum base object
    subroutine initialise_sphere_momentum(self,nb,collect_spherical)
    	implicit none

    	! initialize shape objects
    	class(sphereObj_mom) 		:: self

		logical, optional, intent(in)	 :: collect_spherical
    	integer, dimension(3),intent(in) :: nb

    	!Call constructor of underlying momentum CV object
    	call self%initialise(nb)

    	!Extend to include new variables
    	allocate(self%Fsurface(1,1,1,3))
    	self%Fsurface = 0.d0
        allocate(self%FsurfaceTemp(1,1,1,3))
        self%FsurfaceTemp = 0.d0

		!Set flag to convert collected values to spherical coordinates
		if (present(collect_spherical)) then
			self%collect_spherical = collect_spherical
		else
			self%collect_spherical = .false.
		endif

    end subroutine initialise_sphere_momentum


    !Check if interaction between molecules crosses CV
    subroutine Add_spherical_CV_forces(self,fij,ri,rj)
        use librarymod, only : sphereCV, sphereiser, sphereisev
    	implicit none
    	
        class(sphereObj_mom) :: self

    	real(kind(0.d0)),dimension(3),intent(in)	:: ri,rj,fij
    	real(kind(0.d0)),dimension(3)				:: f_vect, rs

		!Convert collected value to spherical coordinates
		if (self%collect_spherical) then
			rs	= sphereiser(ri)
			f_vect = sphereisev(fij,rs(2),rs(3))
		else
			f_vect = fij
		endif

    	!Add for molecule i
    	self%FsurfaceTemp(1,1,1,:) = self%FsurfaceTemp(1,1,1,:) & 
    			- f_vect*(sphereCV(rj,self%radius)-sphereCV(ri,self%radius))

    end subroutine Add_spherical_CV_forces

    !Check if molecule is inside volume
    subroutine Add_spherical_CV_velocity(self,ri,vi)
        use librarymod, only : sphereCV, sphereiser, sphereisev
    	implicit none
    	
        class(sphereObj_mom) :: self

    	real(kind(0.d0)),dimension(3),intent(in)	:: ri,vi
    	real(kind(0.d0)),dimension(3)				:: vel_vect,rs

		!Convert collected value to spherical coordinates
		if (self%collect_spherical) then
			rs	= sphereiser(ri)
			vel_vect = sphereisev(vi,rs(2),rs(3))
		else
			vel_vect = vi
		endif

    	!Add for molecule i
    	if (.not. allocated(self%Xtemp)) then
    	    allocate(self%Xtemp(1,1,1,3))
    	    self%Xtemp = 0.d0
        endif
    	self%Xtemp(1,1,1,:) = self%Xtemp(1,1,1,:) + vel_vect*sphereCV(ri,self%radius)

    end subroutine Add_spherical_CV_velocity


    !Check if molecule has crossed surface of CV
    subroutine Add_spherical_CV_fluxes(self,vi,ri,ri_tp1)
        use librarymod, only : sphereCV, sphereiser, sphereisev
    	implicit none
    	
        class(sphereObj_mom) :: self

    	real(kind(0.d0)),dimension(3),intent(in)	:: vi,ri,ri_tp1

    	real(kind(0.d0)),dimension(3)				:: vel_vect,rs

		!Convert collected value to spherical coordinates
		if (self%collect_spherical) then
			rs	= sphereiser(ri)
			vel_vect = sphereisev(vi,rs(2),rs(3))
		else
			vel_vect = vi
		endif

    	!Add for molecule i
    	self%flux(1,1,1,:,1) = self%flux(1,1,1,:,1) - vel_vect* & 
				(sphereCV(ri_tp1,self%radius)-sphereCV(ri,self%radius))

    end subroutine Add_spherical_CV_fluxes


    !Update time evolution and store previous two values
    subroutine update_dXdt_sphere_momentum(self, X)
    	implicit none

    	! initialize objects
    	class(sphereObj_mom) :: self

    	double precision,dimension(:,:,:,:),allocatable,intent(in) :: X

    	self%X_minus_t  = self%X
    	self%X 		  = X
    		
    	self%dXdt = self%X - self%X_minus_t

    end subroutine update_dXdt_sphere_momentum



    !Check error for specified range of bins
    subroutine check_error_sphere_momentum(self,imin,imax,jmin,jmax,kmin,kmax,iter,irank)
    	! initialize shape objects
    	use computational_constants_MD, only : domain,delta_t,Nvflux_ave
    	use calculated_properties_MD, only : nbins
		use librarymod, only : sphereiser
    	implicit none

    	! initialize objects
    	class(sphereObj_mom) :: self

    	integer,intent(in) :: iter,irank,imin,imax,jmin,jmax,kmin,kmax

    	logical							:: check_ok
    	integer 						:: i,j,k
    	integer,save 					:: first_time = 0
    	double precision				:: conserved
    	double precision,dimension(3)	:: binsize,totalpressure,totalflux,F_ext,dvelocitydt

    	!First call doesn't have difference in time yet so skip
    	if (first_time .lt. 2) then
    		first_time = first_time + 1
    		print'(a)', '                iter irank  i  j  k   conserved    Forces     fluxes      dvdt      Fext       v        vt-dt'
    		return
    	endif

        !Calculate total CV flux and change in mass
        totalflux(:) = self%flux(1,1,1,:,1)

        !Totalpressure = totalpressure*delta_t
        totalpressure(:) = 0.5d0*delta_t*self%Fsurface(1,1,1,:)
    	    
        call self%update_dXdt(self%Xtemp)
        self%Xtemp = 0.d0

    	!drhou/dt
        dvelocitydt =  self%dxdt(1,1,1,:)/Nvflux_ave !(delta_t*Nvflux_ave)


		if (self%collect_spherical) then
			conserved = dot_product(totalpressure,totalpressure) & 
						-dot_product(totalflux,totalflux)        &
					    -dot_product(dvelocitydt,dvelocitydt)

        	if(abs(conserved) .gt. 0.000000001d0) then
				print'(a,i8,3f18.12)','Error_sphere_mom', iter, & 
					 dot_product(totalpressure,totalpressure), & 
					-dot_product(totalflux,totalflux), & 
					 dot_product(dvelocitydt,dvelocitydt)


			endif

		else
            !Verify that CV momentum is exactly conservative
            conserved = sum(totalpressure-totalflux-dvelocitydt-F_ext)

        	if(abs(conserved) .gt. 0.000000001d0) then
        		print'(a,i8,4i4,7f11.5)','Error_sphere_mom', iter,irank,1,1,1, & 
        				 conserved, sum(totalpressure),-sum(totalflux),sum(dvelocitydt), & 
        			+sum(F_ext), sum(self%X(1,1,1,:)),   & 
        			 sum(self%X_minus_t(1,1,1,:))
        	endif
		endif

! 		if (any(sphereiser(dvelocitydt) .lt. 0.000000000001)) then
! 			!Do nothing
! 		else
! 			print'(a,4f11.5)','qqqq', -sum(sphereiser(totalflux)),sum(sphereiser(dvelocitydt)),-sum(sphereiser(totalflux))/sum(sphereiser(dvelocitydt)),-sum(sphereiser(totalflux)+sphereiser(dvelocitydt))
! 		endif

		if (iter .gt. 5) then
		!	write(57,'(a,i8,3f18.12)'),'Error_sphere_mom', iter, & 
		!			 dot_product(totalpressure,totalpressure),-dot_product(totalflux,totalflux),dot_product(dvelocitydt,dvelocitydt)

 			write(57,'(a,i8,9f14.9)'),'Error_sphere_mom', iter, & 
 				 totalpressure,-totalflux,dvelocitydt

			write(58,'(a,i8,9f11.5)'),'Error_sphere_mom', iter, & 
					 sphereiser(totalpressure),-sphereiser(totalflux),sphereiser(dvelocitydt)
		endif


    	CV_sphere_momentum%Fsurface = CV_sphere_momentum%FsurfaceTemp
    	CV_sphere_momentum%flux = 0.d0; CV_sphere_momentum%FsurfaceTemp = 0.d0;

    end subroutine check_error_sphere_momentum


    end module CV_objects






! !===================================================================================
! ! Momentum Flux over a surface of a bin including all intermediate bins

! subroutine get_molecule_CV_momentum_flux(molno)
! 	use computational_constants_MD, only : domain, delta_t, halfdomain, nhb
! 	use calculated_properties_MD, only : nbins, momentum_flux
! 	use arrays_MD, only : r, v
! 	use librarymod, only : heaviside, imaxloc
! 	implicit none

! 	integer,intent(in)				:: molno

! 	integer							:: i,j,k,jxyz
! 	integer							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
! 	integer		,dimension(3)		:: ibin1,ibin2,cbin
! 	double precision,dimension(3)	:: mbinsize,velvect,crossface
! 	double precision,dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

! 	!CV momentum flux
! 	!Determine bin size
! 	mbinsize(:) = domain(:) / nbins(:)

! 	!Get velocity at v(t+dt/2) from v(t-dt/2)
! 	velvect(:) = v(:,molno)
! 	ri1(:) = r(:,molno) 						!Molecule i at time t
! 	ri2(:) = r(:,molno)	- delta_t*velvect(:)	!Molecule i at time t-dt
! 	ri12   = ri1 - ri2							!Molecule i trajectory between t-dt and t
! 	where (ri12 .eq. 0.d0) ri12 = 0.000001d0

! 	!Assign to bins before and after using integer division
! 	ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb(:)
! 	ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb(:)

! 	!Replace Signum function with this functions which gives a
! 	!check for plane crossing and the correct sign 
! 	crossface(:) =  ibin1(:) - ibin2(:)

! 	if (sum(abs(crossface(:))) .ne. 0) then

! 		do i = ibin1(1),ibin2(1),sign(1,ibin2(1)-ibin1(1))
! 		do j = ibin1(2),ibin2(2),sign(1,ibin2(2)-ibin1(2))
! 		do k = ibin1(3),ibin2(3),sign(1,ibin2(3)-ibin1(3))

! 			cbin(1) = i; cbin(2) = j; cbin(3) = k

! 			bintop(:) = (cbin(:)-1*nhb(:)  )*mbinsize(:)-halfdomain(:)
! 			binbot(:) = (cbin(:)-1*nhb(:)-1)*mbinsize(:)-halfdomain(:)

! 			!Calculate the plane intersect of trajectory with surfaces of the cube
! 			Pxt=(/ 			bintop(1), 		     & 
! 					ri1(2)+(ri12(2)/ri12(1))*(bintop(1)-ri1(1)), & 
! 					ri1(3)+(ri12(3)/ri12(1))*(bintop(1)-ri1(1))  	/)
! 			Pxb=(/ 			binbot(1), 		     & 
! 					ri1(2)+(ri12(2)/ri12(1))*(binbot(1)-ri1(1)), & 
! 					ri1(3)+(ri12(3)/ri12(1))*(binbot(1)-ri1(1))  	/)
! 			Pyt=(/	ri1(1)+(ri12(1)/ri12(2))*(bintop(2)-ri1(2)), & 
! 						bintop(2), 		     & 
! 					ri1(3)+(ri12(3)/ri12(2))*(bintop(2)-ri1(2))  	/)
! 			Pyb=(/	ri1(1)+(ri12(1)/ri12(2))*(binbot(2)-ri1(2)), &
! 						binbot(2), 		     & 
! 					ri1(3)+(ri12(3)/ri12(2))*(binbot(2)-ri1(2))  	/)
! 			Pzt=(/	ri1(1)+(ri12(1)/ri12(3))*(bintop(3)-ri1(3)), & 
! 					ri1(2)+(ri12(2)/ri12(3))*(bintop(3)-ri1(3)), &
! 						bintop(3) 			/)
! 			Pzb=(/	ri1(1)+(ri12(1)/ri12(3))*(binbot(3)-ri1(3)), &
! 					ri1(2)+(ri12(2)/ri12(3))*(binbot(3)-ri1(3)), & 
! 						binbot(3) 			/)

! 			onfacexb =0.5d0*(sign(1.d0,binbot(1) - ri2(1)) 	 & 
! 					       - sign(1.d0,binbot(1) - ri1(1)))* &
! 							(heaviside(bintop(2) - Pxb(2)) 	 &
! 					       - heaviside(binbot(2) - Pxb(2)))* &
! 							(heaviside(bintop(3) - Pxb(3)) 	 &
! 					       - heaviside(binbot(3) - Pxb(3)))
! 			onfaceyb =0.5d0*(sign(1.d0,binbot(2) - ri2(2))   &
! 					       - sign(1.d0,binbot(2) - ri1(2)))* &
! 							(heaviside(bintop(1) - Pyb(1))   &
! 					       - heaviside(binbot(1) - Pyb(1)))* &
! 							(heaviside(bintop(3) - Pyb(3))   &
! 					       - heaviside(binbot(3) - Pyb(3)))
! 			onfacezb =0.5d0*(sign(1.d0,binbot(3) - ri2(3))   &
! 					       - sign(1.d0,binbot(3) - ri1(3)))* &
! 							(heaviside(bintop(1) - Pzb(1))   &
! 					       - heaviside(binbot(1) - Pzb(1)))* &
! 							(heaviside(bintop(2) - Pzb(2))   &
! 					       - heaviside(binbot(2) - Pzb(2)))

! 			onfacext =0.5d0*(sign(1.d0,bintop(1) - ri2(1))   &
! 					       - sign(1.d0,bintop(1) - ri1(1)))* &
! 							(heaviside(bintop(2) - Pxt(2))   &
! 					       - heaviside(binbot(2) - Pxt(2)))* &
! 							(heaviside(bintop(3) - Pxt(3))   &
! 					       - heaviside(binbot(3) - Pxt(3)))
! 			onfaceyt =0.5d0*(sign(1.d0,bintop(2) - ri2(2))   &
! 					       - sign(1.d0,bintop(2) - ri1(2)))* &
! 							(heaviside(bintop(1) - Pyt(1))   &
! 					       - heaviside(binbot(1) - Pyt(1)))* &
! 							(heaviside(bintop(3) - Pyt(3))   &
! 					       - heaviside(binbot(3) - Pyt(3)))
! 			onfacezt =0.5d0*(sign(1.d0,bintop(3) - ri2(3))   &
! 					       - sign(1.d0,bintop(3) - ri1(3)))* &
! 							(heaviside(bintop(1) - Pzt(1))   &
! 						   - heaviside(binbot(1) - Pzt(1)))* &
! 							(heaviside(bintop(2) - Pzt(2))   &
! 					       - heaviside(binbot(2) - Pzt(2)))

! 			jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer

! 			!Calculate velocity at time of intersection
! 			!crosstime = (r(jxyz,n) - rplane)/v(jxyz,n)
! 			!velvect(:) = v(:,n) !- a(:,n) * crosstime
! 			!Change in velocity at time of crossing is not needed as velocity assumed constant 
! 			!for timestep and changes when forces are applied.

! 			!Add Momentum flux over face
! 			momentum_flux(cbin(1),cbin(2),cbin(3),:,1) = & 
! 				momentum_flux(cbin(1),cbin(2),cbin(3),:,1) & 
! 			      - velvect(:)*dble(onfacexb)*abs(crossface(jxyz))
! 			momentum_flux(cbin(1),cbin(2),cbin(3),:,2) = & 
! 				momentum_flux(cbin(1),cbin(2),cbin(3),:,2) & 
! 			      - velvect(:)*dble(onfaceyb)*abs(crossface(jxyz))
! 			momentum_flux(cbin(1),cbin(2),cbin(3),:,3) = & 
! 				momentum_flux(cbin(1),cbin(2),cbin(3),:,3) &
! 			      - velvect(:)*dble(onfacezb)*abs(crossface(jxyz))
! 			momentum_flux(cbin(1),cbin(2),cbin(3),:,4) = & 
! 				momentum_flux(cbin(1),cbin(2),cbin(3),:,4) &
! 			      + velvect(:)*dble(onfacext)*abs(crossface(jxyz))
! 			momentum_flux(cbin(1),cbin(2),cbin(3),:,5) = & 
! 				momentum_flux(cbin(1),cbin(2),cbin(3),:,5) &
! 			      + velvect(:)*dble(onfaceyt)*abs(crossface(jxyz))
! 			momentum_flux(cbin(1),cbin(2),cbin(3),:,6) = & 
! 				momentum_flux(cbin(1),cbin(2),cbin(3),:,6) &
! 			      + velvect(:)*dble(onfacezt)*abs(crossface(jxyz))

! 		enddo
! 		enddo
! 		enddo

! 		!if ((ibin1(1) .eq. 4 .and. ibin1(2) .eq. 4 .and. ibin1(3) .eq. 4) .or. &
! 		!	(ibin2(1) .eq. 4 .and. ibin2(2) .eq. 4 .and. ibin2(3) .eq. 4)) then
! 		!	print'(a,8i4,18f7.3)', 'mol',iter,ibin1,ibin2, molno,momentum_flux(4,4,4,:,:)
! 	!endif
! 			

! 	endif

! end subroutine get_molecule_CV_momentum_flux

! !===================================================================================
! !Forces over the surface of a Volume

! subroutine get_molecule_CV_forces(fij,ri,rj,molnoi,molnoj)
! 	use computational_constants_MD, only : domain, delta_t, halfdomain, nhb
! 	use calculated_properties_MD, only : nbins, momentum_flux, volume_force
! 	use physical_constants_MD, only: np,nd
! 	!use arrays_MD, only : r, v
! 	use librarymod, only : heaviside
! 	implicit none

! 	integer,intent(in)							:: molnoi,molnoj
! 	double precision,dimension(3),intent(in)	:: ri,rj,fij

! 	integer,dimension(3)			:: ibin, jbin
! 	double precision,dimension(3)	:: crossplane,fsurface
! 	double precision,dimension(3)	:: Fbinsize, bintopi, binboti, bintopj, binbotj

! 	!Determine bin size
! 	Fbinsize(:) = domain(:) / nbins(:)

! 	!Assign to bins using integer division
! 	ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:)) + nhb	!Establish current bin
! 	jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:)) + nhb 	!Establish current bin

! 	crossplane(:) =  dble(ibin(:)-jbin(:))

! 	bintopi(:) = (ibin(:)-1)*Fbinsize(:)-halfdomain(:)
! 	binboti(:) = (ibin(:)-2)*Fbinsize(:)-halfdomain(:)
! 	bintopj(:) = (jbin(:)-1)*Fbinsize(:)-halfdomain(:)
! 	binbotj(:) = (jbin(:)-2)*Fbinsize(:)-halfdomain(:)

! 	!Add for molecule i
! 	if(molnoi .le. np) then
! 		fsurface = fij(:)* dble((heaviside(bintopi(1)-ri(1))-heaviside(binboti(1)-ri(1)))* & 
! 			  		(heaviside(bintopi(2)-ri(2))-heaviside(binboti(2)-ri(2)))* & 
! 			  		(heaviside(bintopi(3)-ri(3))-heaviside(binboti(3)-ri(3)))- & 
! 			  		(heaviside(bintopi(1)-rj(1))-heaviside(binboti(1)-rj(1)))* & 
! 			  		(heaviside(bintopi(2)-rj(2))-heaviside(binboti(2)-rj(2)))* & 
! 			  		(heaviside(bintopi(3)-rj(3))-heaviside(binboti(3)-rj(3))))
! 		volume_force(ibin(1),ibin(2),ibin(3),:,1) = volume_force(ibin(1),ibin(2),ibin(3),:,1) + fsurface*delta_t
! 	endif

! 	!Add for molecule j
! 	if(molnoj .le. np) then
! 		fsurface = fij(:)* dble((heaviside(bintopj(1)-ri(1))-heaviside(binbotj(1)-ri(1)))* & 
! 			  		(heaviside(bintopj(2)-ri(2))-heaviside(binbotj(2)-ri(2)))* & 
! 			  		(heaviside(bintopj(3)-ri(3))-heaviside(binbotj(3)-ri(3)))- & 
! 			  		(heaviside(bintopj(1)-rj(1))-heaviside(binbotj(1)-rj(1)))* & 
! 			  		(heaviside(bintopj(2)-rj(2))-heaviside(binbotj(2)-rj(2)))* & 
! 			  		(heaviside(bintopj(3)-rj(3))-heaviside(binbotj(3)-rj(3))))
! 		volume_force(jbin(1),jbin(2),jbin(3),:,1) = volume_force(jbin(1),jbin(2),jbin(3),:,1) + fsurface*delta_t
! 	endif

! end subroutine get_molecule_CV_forces

! !===================================================================================
! !Forces over the surface of a Volume

! subroutine get_molecule_CV_stresses(fij,ri,rj,molnoi)
! 	use computational_constants_MD, only : domain, delta_t, halfdomain, nhb, eflux_outflag
! 	use calculated_properties_MD, only : nbins, momentum_flux, volume_force, Pxyface, Pxyvface
! 	use physical_constants_MD, only: np,nd
! 	use arrays_MD, only : r, v
! 	use librarymod, only : heaviside
! 	implicit none


! 	integer,intent(in)							:: molnoi
! 	double precision,dimension(3),intent(in)	:: ri,rj,fij

! 	integer							:: i,j,k,ixyz
! 	integer							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
! 	integer,dimension(3)			:: cbin, ibin, jbin
! 	double precision,dimension(3)	:: rij,fsurface,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb,velvect
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
! 		Pxt=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
! 		Pxb=(/ binbot(1),ri(2)+(rij(2)/rij(1))*(binbot(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(binbot(1)-ri(1))  /)
! 		Pyt=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
! 		Pyb=(/ri(1)+(rij(1)/rij(2))*(binbot(2)-ri(2)), binbot(2),ri(3)+(rij(3)/rij(2))*(binbot(2)-ri(2))  /)
! 		Pzt=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)
! 		Pzb=(/ri(1)+(rij(1)/rij(3))*(binbot(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(binbot(3)-ri(3)), binbot(3) /)

! 		onfacexb =  	(sign(1.d0,binbot(1)- rj(1)) - sign(1.d0,binbot(1)- ri(1)))* &
! 						(heaviside(bintop(2)-Pxb(2)) - heaviside(binbot(2)-Pxb(2)))* &
! 						(heaviside(bintop(3)-Pxb(3)) - heaviside(binbot(3)-Pxb(3)))
! 		onfaceyb =  	(sign(1.d0,binbot(2)- rj(2)) - sign(1.d0,binbot(2)- ri(2)))* &
! 						(heaviside(bintop(1)-Pyb(1)) - heaviside(binbot(1)-Pyb(1)))* &
! 						(heaviside(bintop(3)-Pyb(3)) - heaviside(binbot(3)-Pyb(3)))
! 		onfacezb =  	(sign(1.d0,binbot(3)- rj(3)) - sign(1.d0,binbot(3)- ri(3)))* &
! 						(heaviside(bintop(1)-Pzb(1)) - heaviside(binbot(1)-Pzb(1)))* &
! 						(heaviside(bintop(2)-Pzb(2)) - heaviside(binbot(2)-Pzb(2)))

! 		onfacext =  	(sign(1.d0,bintop(1)- rj(1)) - sign(1.d0,bintop(1)- ri(1)))* &
! 						(heaviside(bintop(2)-Pxt(2)) - heaviside(binbot(2)-Pxt(2)))* &
! 	            		(heaviside(bintop(3)-Pxt(3)) - heaviside(binbot(3)-Pxt(3)))
! 		onfaceyt = 		(sign(1.d0,bintop(2)- rj(2)) - sign(1.d0,bintop(2)- ri(2)))* &
! 						(heaviside(bintop(1)-Pyt(1)) - heaviside(binbot(1)-Pyt(1)))* &
! 						(heaviside(bintop(3)-Pyt(3)) - heaviside(binbot(3)-Pyt(3)))
! 		onfacezt =  	(sign(1.d0,bintop(3)- rj(3)) - sign(1.d0,bintop(3)- ri(3)))* &
! 						(heaviside(bintop(1)-Pzt(1)) - heaviside(binbot(1)-Pzt(1)))* &
! 						(heaviside(bintop(2)-Pzt(2)) - heaviside(binbot(2)-Pzt(2)))

! 		!Stress acting on face over volume
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,1) = Pxyface(cbin(1),cbin(2),cbin(3),:,1) + fij(:)*dble(onfacexb)
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,2) = Pxyface(cbin(1),cbin(2),cbin(3),:,2) + fij(:)*dble(onfaceyb)
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,3) = Pxyface(cbin(1),cbin(2),cbin(3),:,3) + fij(:)*dble(onfacezb)
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,4) = Pxyface(cbin(1),cbin(2),cbin(3),:,4) + fij(:)*dble(onfacext)
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,5) = Pxyface(cbin(1),cbin(2),cbin(3),:,5) + fij(:)*dble(onfaceyt)
! 		Pxyface(cbin(1),cbin(2),cbin(3),:,6) = Pxyface(cbin(1),cbin(2),cbin(3),:,6) + fij(:)*dble(onfacezt)

! 		!Stress acting on face over volume
! 		if (eflux_outflag .ne. 0) then
! 			velvect(:) = v(:,molnoi) 
! 			!velvect(:) = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),1) = Pxyvface(cbin(1),cbin(2),cbin(3),1) + dot_product(fij,velvect)*dble(onfacexb)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),2) = Pxyvface(cbin(1),cbin(2),cbin(3),2) + dot_product(fij,velvect)*dble(onfaceyb)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),3) = Pxyvface(cbin(1),cbin(2),cbin(3),3) + dot_product(fij,velvect)*dble(onfacezb)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),4) = Pxyvface(cbin(1),cbin(2),cbin(3),4) + dot_product(fij,velvect)*dble(onfacext)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),5) = Pxyvface(cbin(1),cbin(2),cbin(3),5) + dot_product(fij,velvect)*dble(onfaceyt)
! 			Pxyvface(cbin(1),cbin(2),cbin(3),6) = Pxyvface(cbin(1),cbin(2),cbin(3),6) + dot_product(fij,velvect)*dble(onfacezt)
! 		endif

! 		!Force applied to volume
! 		fsurface(:) = 0.d0
! 		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfacexb - onfacext)
! 		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfaceyb - onfaceyt)
! 		fsurface(:) = fsurface(:) + 0.25d0*fij(:)*dble(onfacezb - onfacezt)
! 		volume_force(cbin(1),cbin(2),cbin(3),:,1) = volume_force(cbin(1),cbin(2),cbin(3),:,1) + fsurface*delta_t

! 	enddo
! 	enddo
! 	enddo

! end subroutine get_molecule_CV_stresses


! !===================================================================================
! !
! !	Routines to obtain forces and fluxes for entire CV
! !
! !===================================================================================

! !-----------------------------------------------------------------------------------
! ! Flux of molecules over bin surface requires all molecules in current cell and 
! ! surrounding cells to be checked

! subroutine get_CV_momentum_flux(icell,jcell,kcell,isumflux,Flux)
! 	use module_linklist, only : node, cell
! 	use arrays_MD, only : r,v,a
! 	use computational_constants_MD, only : delta_t, domain,halfdomain, nhb,iter
! 	use calculated_properties_MD, only	 : nbins
! 	use librarymod, only : imaxloc
! 	implicit none

! 	integer,intent(in)										:: icell,jcell,kcell
! 	double precision,intent(inout)							:: isumflux
! 	double precision,optional,dimension(3,6),intent(out)	:: Flux

! 	integer								:: n,molno
! 	integer								:: icellshift,jcellshift,kcellshift,adjacentcellnp
! 	double precision,dimension(3)		:: mbinsize,ri1,ri2,Fsurface,velvect
! 	type(node), pointer		 	        :: old, current


! 	!Determine bin size
! 	mbinsize(:) = domain(:) / nbins(:)

! 	!Calculate bin surface fluxes
! 	do kcellshift = -1,1
! 	do jcellshift = -1,1
! 	do icellshift = -1,1
! 		old =>  cell%head(icell+icellshift, & 
! 				  		  jcell+jcellshift, & 
! 				  		  kcell+kcellshift)%point
! 		adjacentcellnp = cell%cellnp(icell+icellshift, & 
! 					     			 jcell+jcellshift, & 
! 					     			 kcell+kcellshift)

! 		do n = 1,adjacentcellnp          !Step through all adjacent cells' molecules

! 			molno = old%molno			!Number of molecule

! 			velvect(:) = v(:,molno)
! 			ri1(:) = r(:,molno) 							!Molecule i at time t
! 			ri2(:) = r(:,molno)	- delta_t*velvect			!Molecule i at time t-dt

! 			Fsurface = 0.d0
! 			! *********************************************************************************
! 			!Calculate flux over surface only if molecule is entering/leaving bin of interest
! 			if (present(Flux)) then
! 				call get_Flux(icell,jcell,kcell,velvect,ri1,ri2,Fsurface,Flux)	
! 			else
! 				call get_Flux(icell,jcell,kcell,velvect,ri1,ri2,Fsurface)	
! 			endif

! 			!print'(a,i4,6i2,i4,18f6.3)', '!CV',iter,icell,jcell,kcell,icell+icellshift,jcell+jcellshift,kcell+kcellshift,molno,Flux

! 			isumflux = isumflux + Fsurface(1)

! 			! *********************************************************************************
! 			current => old
! 			old => current%next    !Use pointer in datatype to obtain next item in list

! 		enddo
! 	enddo
! 	enddo
! 	enddo

! contains

! 	!-----------------------------------------------------------------------------------
! 	!Flux over the surface of a bin

! 	subroutine get_Flux(icell,jcell,kcell,velvect,ri1,ri2,Fsurface,flux)
! 		use computational_constants_MD, only : domain,halfdomain, iter
! 		use calculated_properties_MD, only	 : nbins
! 		use librarymod, only : heaviside, imaxloc
! 		implicit none

! 		integer,intent(in)										:: icell,jcell,kcell
! 		double precision,dimension(3),intent(in)				:: ri1, ri2, velvect
! 		double precision,dimension(3),intent(out)				:: Fsurface
! 		double precision,dimension(3,6),intent(inout),optional	:: Flux

! 		integer									:: ixyz,jxyz,kxyz,i,j,k,ii,jj,kk,n
! 		integer									:: planeno
! 		integer									:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
! 		integer		,dimension(3)				:: ibin1,ibin2,cbin
! 		double precision						:: crosstime,crossplane,rplane,shift
! 		double precision,dimension(3)			:: mbinsize,crossface
! 		double precision,dimension(3)			:: ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
! 		double precision,dimension(3,3,3,3)		:: Fsurfacebins
! 		double precision,dimension(3,3,3,3,6)	:: Fluxbins

! 		!CV momentum flux
! 		!Determine bin size
! 		mbinsize(:) = domain(:) / nbins(:)

! 		ri12   = ri1 - ri2		!Molecule i trajectory between t-dt and t
! 		where (ri12 .eq. 0.d0) ri12 = 0.000001d0

! 		!Assign to bins before and after using integer division
! 		ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb(:)
! 		ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb(:)

! 		!Replace Signum function with this functions which gives a
! 		!check for plane crossing and the correct sign 
! 		crossface(:) =  ibin1(:) - ibin2(:)

! 		if (sum(abs(crossface(:))) .ne. 0) then

! 			Fluxbins = 0.d0

! 			do i = ibin1(1),ibin2(1),sign(1,ibin2(1)-ibin1(1))
! 			do j = ibin1(2),ibin2(2),sign(1,ibin2(2)-ibin1(2))
! 			do k = ibin1(3),ibin2(3),sign(1,ibin2(3)-ibin1(3))

! 				cbin(1) = i; cbin(2) = j; cbin(3) = k

! 				bintop(:) = (cbin(:)-1*nhb(:)  )*mbinsize(:)-halfdomain(:)
! 				binbot(:) = (cbin(:)-1*nhb(:)-1)*mbinsize(:)-halfdomain(:)

! 				!Calculate the plane intersect of trajectory with surfaces of the cube
! 				Pxt=(/ 			bintop(1), 		     & 
! 						ri1(2)+(ri12(2)/ri12(1))*(bintop(1)-ri1(1)), & 
! 						ri1(3)+(ri12(3)/ri12(1))*(bintop(1)-ri1(1))  	/)
! 				Pxb=(/ 			binbot(1), 		     & 
! 						ri1(2)+(ri12(2)/ri12(1))*(binbot(1)-ri1(1)), & 
! 						ri1(3)+(ri12(3)/ri12(1))*(binbot(1)-ri1(1))  	/)
! 				Pyt=(/	ri1(1)+(ri12(1)/ri12(2))*(bintop(2)-ri1(2)), & 
! 								bintop(2), 		     & 
! 						ri1(3)+(ri12(3)/ri12(2))*(bintop(2)-ri1(2))  	/)
! 				Pyb=(/	ri1(1)+(ri12(1)/ri12(2))*(binbot(2)-ri1(2)), &
! 								binbot(2), 		     & 
! 						ri1(3)+(ri12(3)/ri12(2))*(binbot(2)-ri1(2))  	/)
! 				Pzt=(/	ri1(1)+(ri12(1)/ri12(3))*(bintop(3)-ri1(3)), & 
! 						ri1(2)+(ri12(2)/ri12(3))*(bintop(3)-ri1(3)), &
! 								bintop(3) 			/)
! 				Pzb=(/	ri1(1)+(ri12(1)/ri12(3))*(binbot(3)-ri1(3)), &
! 						ri1(2)+(ri12(2)/ri12(3))*(binbot(3)-ri1(3)), & 
! 								binbot(3) 			/)

! 				onfacexb =0.5d0*(sign(1.d0,binbot(1) - ri2(1)) 	 & 
! 						       - sign(1.d0,binbot(1) - ri1(1)))* &
! 								(heaviside(bintop(2) - Pxb(2)) 	 &
! 						       - heaviside(binbot(2) - Pxb(2)))* &
! 								(heaviside(bintop(3) - Pxb(3)) 	 &
! 						       - heaviside(binbot(3) - Pxb(3)))
! 				onfaceyb =0.5d0*(sign(1.d0,binbot(2) - ri2(2))   &
! 						       - sign(1.d0,binbot(2) - ri1(2)))* &
! 								(heaviside(bintop(1) - Pyb(1))   &
! 						       - heaviside(binbot(1) - Pyb(1)))* &
! 								(heaviside(bintop(3) - Pyb(3))   &
! 						       - heaviside(binbot(3) - Pyb(3)))
! 				onfacezb =0.5d0*(sign(1.d0,binbot(3) - ri2(3))   &
! 						       - sign(1.d0,binbot(3) - ri1(3)))* &
! 								(heaviside(bintop(1) - Pzb(1))   &
! 						       - heaviside(binbot(1) - Pzb(1)))* &
! 								(heaviside(bintop(2) - Pzb(2))   &
! 						       - heaviside(binbot(2) - Pzb(2)))

! 				onfacext =0.5d0*(sign(1.d0,bintop(1) - ri2(1))   &
! 						       - sign(1.d0,bintop(1) - ri1(1)))* &
! 								(heaviside(bintop(2) - Pxt(2))   &
! 						       - heaviside(binbot(2) - Pxt(2)))* &
! 			            		(heaviside(bintop(3) - Pxt(3))   &
! 						       - heaviside(binbot(3) - Pxt(3)))
! 				onfaceyt =0.5d0*(sign(1.d0,bintop(2) - ri2(2))   &
! 						       - sign(1.d0,bintop(2) - ri1(2)))* &
! 								(heaviside(bintop(1) - Pyt(1))   &
! 						       - heaviside(binbot(1) - Pyt(1)))* &
! 								(heaviside(bintop(3) - Pyt(3))   &
! 					    	   - heaviside(binbot(3) - Pyt(3)))
! 				onfacezt =0.5d0*(sign(1.d0,bintop(3) - ri2(3))   &
! 						       - sign(1.d0,bintop(3) - ri1(3)))* &
! 								(heaviside(bintop(1) - Pzt(1))   &
! 						       - heaviside(binbot(1) - Pzt(1)))* &
! 								(heaviside(bintop(2) - Pzt(2))   &
! 						       - heaviside(binbot(2) - Pzt(2)))

! 				fsurface(:) = 0.d0
! 				fsurface(:) = fsurface(:) - 0.5d0*velvect(:)*dble(onfacexb - onfacext)
! 				fsurface(:) = fsurface(:) - 0.5d0*velvect(:)*dble(onfaceyb - onfaceyt)
! 				fsurface(:) = fsurface(:) - 0.5d0*velvect(:)*dble(onfacezb - onfacezt)

! 				jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer

! 				!Add Momentum flux over face
! 				ii = modulo(i,3)+1; jj = modulo(j,3)+1; kk = modulo(k,3)+1
! 				Fluxbins(ii,jj,kk,:,1) = Fluxbins(ii,jj,kk,:,1) - velvect(:)*dble(onfacexb)*abs(crossface(jxyz))
! 				Fluxbins(ii,jj,kk,:,2) = Fluxbins(ii,jj,kk,:,2) - velvect(:)*dble(onfaceyb)*abs(crossface(jxyz))
! 				Fluxbins(ii,jj,kk,:,3) = Fluxbins(ii,jj,kk,:,3) - velvect(:)*dble(onfacezb)*abs(crossface(jxyz))
! 				Fluxbins(ii,jj,kk,:,4) = Fluxbins(ii,jj,kk,:,4) + velvect(:)*dble(onfacext)*abs(crossface(jxyz))
! 				Fluxbins(ii,jj,kk,:,5) = Fluxbins(ii,jj,kk,:,5) + velvect(:)*dble(onfaceyt)*abs(crossface(jxyz))
! 				Fluxbins(ii,jj,kk,:,6) = Fluxbins(ii,jj,kk,:,6) + velvect(:)*dble(onfacezt)*abs(crossface(jxyz))


! 			enddo
! 			enddo
! 			enddo
! 			
! 			!Add surface force to current bin
! 			!Fsurface(:) = sum(Fluxbins(modulo(icell,3)+1, 	& 
! 			!	     		   modulo(jcell,3)+1, 	& 
! 			!	    		   modulo(kcell,3)+1,:,:),2)

! 			!Take flux from central bin only
! 			if (present(Flux))then
! 				Flux(:,:) = Flux(:,:) +  Fluxbins(modulo(icell,3)+1, 	& 
! 				     		     	 			  modulo(jcell,3)+1, 	& 
! 				    		     	 	          modulo(kcell,3)+1,:,:)
! 			endif


! 			!if ((ibin1(1) .eq. 4 .and. ibin1(2) .eq. 4 .and. ibin1(3) .eq. 4) .or. &
! 			!	(ibin2(1) .eq. 4 .and. ibin2(2) .eq. 4 .and. ibin2(3) .eq. 4)) then
! 			!	print'(a,i4,6i2,i4,18f6.3)', 'CV ',iter,ibin1,ibin2,molno,Flux!Fluxbins(modulo(icell,3)+1,modulo(jcell,3)+1,modulo(kcell,3)+1,:,:)!
! 			!endif
! 					
! 		endif

! 	end subroutine get_Flux

! end subroutine get_CV_momentum_flux

! !-----------------------------------------------------------------------------------
! ! Forces between current bin and surrounding bins

! subroutine get_CV_force(icell,jcell,kcell,isumforce,Traction)
! 	use module_linklist, only : node, cell
! 	use arrays_MD, only : r
! 	use computational_constants_MD, only : domain,halfdomain
! 	use calculated_properties_MD, only	 : nbins, Pxyface
! 	use physical_constants_MD, only 	 : rcutoff2,np,nd
! 	implicit none

! 	integer,intent(in)										:: icell,jcell,kcell
! 	double precision,intent(inout)							:: isumforce
! 	double precision,dimension(3,6),optional,intent(inout)	:: Traction

! 	integer							:: n,j,ixyz,molnoi,molnoj
! 	integer							:: icellshift,jcellshift,kcellshift,cellnp,adjacentcellnp
! 	double precision				:: rij2, invrij2, accijmag
! 	double precision,dimension(3)	:: ri,rj,rij,fij,Fsurface
! 	type(node), pointer		 	 	:: oldi, currenti, oldj, currentj

! 	cellnp = cell%cellnp(icell,jcell,kcell)
! 	oldi => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list

! 	!Calculate averages for bin
! 	do n = 1, cellnp    ! Loop over all particles

! 		molnoi = oldi%molno	!Number of molecule
! 		ri = r(:,molnoi)	!Retrieve ri

! 		!Calculate bin surface Forces
! 		do kcellshift = -1,1
! 		do jcellshift = -1,1
! 		do icellshift = -1,1
! 			oldj => cell%head(icell+icellshift, & 
! 					  jcell+jcellshift, & 
! 					  kcell+kcellshift)%point
! 			adjacentcellnp = cell%cellnp(icell+icellshift, & 
! 						     jcell+jcellshift, & 
! 						     kcell+kcellshift)

! 			do j = 1,adjacentcellnp          !Step through all j for each i

! 				molnoj = oldj%molno 	 !Number of molecule
! 				rj = r(molnoj,:)         !Retrieve rj

! 				currentj => oldj
! 				oldj => currentj%next    !Use pointer in datatype to obtain next item in list

! 				if(molnoi==molnoj) cycle !Check to prevent interaction with self

! 				rij2=0                   !Set rij^2 to zero
! 				rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j

! 				!rij2 = dot_product(rij)
! 				do ixyz=1,nd
! 					rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
! 				enddo

! 				if (rij2 < rcutoff2) then

! 					invrij2 = 1.d0/rij2                 !Invert value

! 					!Linear magnitude of acceleration for each molecule
! 					accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)

! 					!Get force and add to bin total
! 					fij = accijmag*rij(:)
! 					Fsurface = 0.d0

! 					!print*, 'FORCE', fij
! 					if (present(Traction)) then
! 						call get_Traction(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface,Traction)
! 					else
! 						call get_Fsurface(molnoi,molnoj,fij,ri,rj,Fsurface)
! 					endif
! 					isumforce = isumforce +  Fsurface(1)

! 				endif
! 			enddo
! 		enddo
! 		enddo
! 		enddo

! 		currenti => oldi
! 		oldi => currenti%next !Use pointer in datatype to obtain next item in list
! 	enddo

! contains


! 	!-----------------------------------------------------------------------------------
! 	!Forces over all of the surface of a bin

! 	subroutine get_Fsurface(molnoi,molnoj,fij,ri,rj,Fsurface)
! 		use librarymod, only : heaviside
! 		implicit none

! 		integer,intent(in)							:: molnoi,molnoj
! 		double precision,dimension(3),intent(in)	:: ri, rj, fij
! 		double precision,dimension(3),intent(out)	:: Fsurface

! 		integer							:: ixyz
! 		integer,dimension(3)			:: ibin, jbin
! 		double precision,dimension(3)	:: Fbinsize, bintop, binbot, crossplane

! 		!Determine bin size
! 		Fbinsize(:) = domain(:) / nbins(:)

! 		!Assign to bins using integer division
! 		ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+1  !Establish current bin

! 		bintop(:) = (ibin(:)-1)*Fbinsize(:)-halfdomain(:)
! 		binbot(:) = (ibin(:)-2)*Fbinsize(:)-halfdomain(:)

! 		!Add for molecule i
! 		if(molnoi .le. np) then
! 			Fsurface = fij(:)* dble((heaviside(bintop(1)-ri(1))-heaviside(binbot(1)-ri(1)))* & 
! 				  					(heaviside(bintop(2)-ri(2))-heaviside(binbot(2)-ri(2)))* & 
! 				  					(heaviside(bintop(3)-ri(3))-heaviside(binbot(3)-ri(3)))- & 
! 				  					(heaviside(bintop(1)-rj(1))-heaviside(binbot(1)-rj(1)))* & 
! 				  					(heaviside(bintop(2)-rj(2))-heaviside(binbot(2)-rj(2)))* & 
! 				  					(heaviside(bintop(3)-rj(3))-heaviside(binbot(3)-rj(3))))
! 		endif

! 	end subroutine get_Fsurface

! 	!-----------------------------------------------------------------------------------
! 	! Tractions on one surface of a bin

! 	subroutine get_Traction(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface,Traction)
! 		use librarymod, only : heaviside
! 		implicit none

! 		integer,intent(in)										:: molnoi,molnoj
! 		double precision,dimension(3),intent(in)				:: ri,rj,fij
! 		double precision,dimension(3),intent(out)				:: Fsurface
! 		double precision,dimension(3,6),intent(inout),optional	:: Traction

! 		integer									:: i,j,k,ixyz,n,tempi
! 		integer									:: icell,jcell,kcell
! 		integer									:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
! 		integer,dimension(3)					:: cbin, ibin, jbin
! 		double precision						:: binforce
! 		double precision,dimension(3)			:: rij,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
! 		double precision,dimension(3)			:: Fbinsize, bintop, binbot
! 		double precision,dimension(3,3,3,3,6)	:: Tractionbins
! 		!Calculate rij
! 		rij = ri - rj
! 		!Prevent Division by zero
! 		do ixyz = 1,3
! 			if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
! 		enddo

! 		!Determine bin size
! 		Fbinsize(:) = domain(:) / nbins(:)

! 		!Assign to bins using integer division
! 		ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+1	!Establish current bin
! 		jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+1 	!Establish current bin

! 		do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
! 		do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
! 		do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

! 			cbin(1) = i; cbin(2) = j; cbin(3) = k

! 			bintop(:) = (ibin(:)-1)*Fbinsize(:)-halfdomain(:)
! 			binbot(:) = (ibin(:)-2)*Fbinsize(:)-halfdomain(:)

! 			if(present(Traction)) then

! 				!Calculate the plane intersect of line with surfaces of the cube
! 				Pxt=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
! 				Pxb=(/ binbot(1),ri(2)+(rij(2)/rij(1))*(binbot(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(binbot(1)-ri(1))  /)
! 				Pyt=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
! 				Pyb=(/ri(1)+(rij(1)/rij(2))*(binbot(2)-ri(2)), binbot(2),ri(3)+(rij(3)/rij(2))*(binbot(2)-ri(2))  /)
! 				Pzt=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)
! 				Pzb=(/ri(1)+(rij(1)/rij(3))*(binbot(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(binbot(3)-ri(3)), binbot(3) /)

! 				onfacexb =   (sign(1.d0,binbot(1)- rj(1))   &
! 							- sign(1.d0,binbot(1)- ri(1)))* &
! 							 (heaviside(bintop(2)-Pxb(2))   &
! 							- heaviside(binbot(2)-Pxb(2)))* &
! 							 (heaviside(bintop(3)-Pxb(3))   & 
! 							- heaviside(binbot(3)-Pxb(3)))
! 				onfaceyb =	 (sign(1.d0,binbot(2)- rj(2))	&
! 							- sign(1.d0,binbot(2)- ri(2)))* &
! 							 (heaviside(bintop(1)-Pyb(1))	& 
! 							- heaviside(binbot(1)-Pyb(1)))* &
! 							 (heaviside(bintop(3)-Pyb(3))	&
! 							- heaviside(binbot(3)-Pyb(3)))
! 				onfacezb =	 (sign(1.d0,binbot(3)- rj(3)) 	&
! 							- sign(1.d0,binbot(3)- ri(3)))* &
! 							 (heaviside(bintop(1)-Pzb(1))	&
! 							- heaviside(binbot(1)-Pzb(1)))* &
! 							 (heaviside(bintop(2)-Pzb(2))	&
! 							- heaviside(binbot(2)-Pzb(2)))

! 				onfacext =	 (sign(1.d0,bintop(1)- rj(1))	&
! 							- sign(1.d0,bintop(1)- ri(1)))* &
! 							 (heaviside(bintop(2)-Pxt(2))	&
! 							- heaviside(binbot(2)-Pxt(2)))* &
! 			            	 (heaviside(bintop(3)-Pxt(3))	&
! 							- heaviside(binbot(3)-Pxt(3)))
! 				onfaceyt =	 (sign(1.d0,bintop(2)- rj(2))	&
! 							- sign(1.d0,bintop(2)- ri(2)))* &
! 							 (heaviside(bintop(1)-Pyt(1))	&
! 							- heaviside(binbot(1)-Pyt(1)))* &
! 							 (heaviside(bintop(3)-Pyt(3))	&
! 							- heaviside(binbot(3)-Pyt(3)))
! 				onfacezt =	 (sign(1.d0,bintop(3)- rj(3))	&
! 							- sign(1.d0,bintop(3)- ri(3)))* &
! 							 (heaviside(bintop(1)-Pzt(1))	&
! 							- heaviside(binbot(1)-Pzt(1)))* &
! 							 (heaviside(bintop(2)-Pzt(2))	&
! 							- heaviside(binbot(2)-Pzt(2)))

! 				!Prevent halo molecules from being included but include molecule which have left domain 
! 				!before rebuild has been triggered.
! 				if (molnoi .gt. np .or. molnoj .gt. np) then
! 					if (cbin(1) .lt. 2 .or. cbin(1) .gt. nbins(1)+1) cycle
! 					if (cbin(2) .lt. 2 .or. cbin(2) .gt. nbins(2)+1) cycle
! 					if (cbin(3) .lt. 2 .or. cbin(3) .gt. nbins(3)+1) cycle
! 				endif

! 				!Stress acting on face over volume
! 				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,1) = & 
! 					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,1) + fij(:)*dble(onfacexb)
! 				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,2) = & 
! 					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,2) + fij(:)*dble(onfaceyb)
! 				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,3) = & 
! 					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,3) + fij(:)*dble(onfacezb)
! 				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,4) = & 
! 					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,4) + fij(:)*dble(onfacext)
! 				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,5) = & 
! 					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,5) + fij(:)*dble(onfaceyt)
! 				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,6) = & 
! 					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,6) + fij(:)*dble(onfacezt)	


! 				!Force applied to volume
! 				!fsurface(:) = 0.25d0*fij(:)*dble(onfacexb - onfacext)
! 				!fsurface(:) = 0.25d0*fij(:)*dble(onfaceyb - onfaceyt)
! 				!fsurface(:) = 0.25d0*fij(:)*dble(onfacezb - onfacezt)

! 				!Add surface force to current bin
! 				!Traction(:,1) = Traction(:,1) + 0.25d0*fij(:)*dble(onfacexb)
! 				!Traction(:,2) = Traction(:,2) + 0.25d0*fij(:)*dble(onfaceyb)
! 				!Traction(:,3) = Traction(:,3) + 0.25d0*fij(:)*dble(onfacezb)
! 				!Traction(:,4) = Traction(:,4) + 0.25d0*fij(:)*dble(onfacext)
! 				!Traction(:,5) = Traction(:,5) + 0.25d0*fij(:)*dble(onfaceyt)
! 				!Traction(:,6) = Traction(:,6) + 0.25d0*fij(:)*dble(onfacezt)

! 				!Add for molecule i
! 				Fsurface = fij(:)* dble((heaviside(bintop(1)-ri(1))-heaviside(binbot(1)-ri(1)))* & 
! 					  					(heaviside(bintop(2)-ri(2))-heaviside(binbot(2)-ri(2)))* & 
! 					  					(heaviside(bintop(3)-ri(3))-heaviside(binbot(3)-ri(3)))- & 
! 					  					(heaviside(bintop(1)-rj(1))-heaviside(binbot(1)-rj(1)))* & 
! 					  					(heaviside(bintop(2)-rj(2))-heaviside(binbot(2)-rj(2)))* & 
! 					  					(heaviside(bintop(3)-rj(3))-heaviside(binbot(3)-rj(3))))

! 				!if (onfaceyb.ne.0.or.onfaceyt.ne.0) print'(9i8)', iter, molnoi, molnoj,np, ibin,onfaceyb,onfaceyt

! 			else
! 				!Add for molecule i
! 				Fsurface = fij(:)* dble((heaviside(bintop(1)-ri(1))-heaviside(binbot(1)-ri(1)))* & 
! 					  					(heaviside(bintop(2)-ri(2))-heaviside(binbot(2)-ri(2)))* & 
! 					  					(heaviside(bintop(3)-ri(3))-heaviside(binbot(3)-ri(3)))- & 
! 					  					(heaviside(bintop(1)-rj(1))-heaviside(binbot(1)-rj(1)))* & 
! 					  					(heaviside(bintop(2)-rj(2))-heaviside(binbot(2)-rj(2)))* & 
! 					  					(heaviside(bintop(3)-rj(3))-heaviside(binbot(3)-rj(3))))
! 			endif

! 		enddo
! 		enddo
! 		enddo

! 		!Take flux from central bin only
! 		if (present(Traction))then
! 			Traction(:,:) = Traction(:,:) +  Tractionbins(modulo(icell,3)+1, 	& 
! 							     		     			  modulo(jcell,3)+1, 	& 
! 			    						     			  modulo(kcell,3)+1,:,:)


! 			if (icell .eq. 5 .and. kcell .eq. 3) then
! 				print'(3i8,4f10.5)',icell,jcell,kcell, Tractionbins(modulo(icell,3)+1, 	& 
! 							     		     			  modulo(jcell,3)+1, 	& 
! 			    						     			  modulo(kcell,3)+1,1,2),& 
! 											 Tractionbins(modulo(icell,3)+1, 	& 
! 							     		     			  modulo(jcell,3)+1, 	& 
! 			    						     			  modulo(kcell,3)+1,1,5) & 
! 						, Pxyface(icell,jcell,kcell,1,2),Pxyface(icell,jcell,kcell,1,5)
! 			endif
! 		endif
! 			
! 	end subroutine get_Traction


! end subroutine get_CV_force

