!=============================================================================
!				   
! Routine which interface with the coupler to the CFD code
! Corresponding empty shell dummy routine used for uncoupled calculation
!
!  Lucian Anton, November 2011
!
! Gave up dummy subroutines for preprocessor flags
!
! LA, February 2012
!
!-----------------------------------------------------------------------------


module md_coupler_socket
	implicit none

#if USE_COUPLER 

    ! CFD id
    integer cfd_code_id 

! type used in applying continuum force
	type cfd_box_sum
		integer np
		real(kind=kind(0.d0))  v(3)
		real(kind=kind(0.d0))  a(3)
	end type cfd_box_sum

 real(kind(0.d0)), allocatable :: vel_cfd(:,:,:,:,:), vbuff(:,:,:,:)
 type(cfd_box_sum),allocatable :: box_average(:,:,:)

contains

!=============================================================================
! Setup initial times based on coupled calculation
!-----------------------------------------------------------------------------
subroutine socket_coupler_init
    use interfaces
	use coupler
	use computational_constants_MD, only : npx,npy,npz,delta_t,elapsedtime,Nsteps,initialstep, delta_t
	use messenger, only	      :  myid, icoord, icomm_grid
	implicit none

 	integer  naverage, nsteps_cfd
	real(kind(0.d0)) :: delta_t_CFD

	! if coupled calculation prepare exchange layout
	call coupler_md_init(npx,npy,npz,icoord,icomm_grid,delta_t)

	! fix NSTEPS for the coupled case
    nsteps_cfd = coupler_md_get_nsteps()
    naverage   = coupler_md_get_md_steps_per_dt_cfd()
    
	Nsteps = initialstep + nsteps_cfd * naverage
	elapsedtime = elapsedtime + nsteps_cfd * naverage * delta_t

	if (myid .eq. 0) then 
		write(*,'(4(a,/),I7,/a,/a,E10.4,a/,a,/a)') &
				"*********************************************************************",	&
 				"WARNING - this is a coupled run which resets the number	      ", 		&
				" of extrasteps to:						   ", 		&
				"								     ", 		&
				nsteps_cfd*naverage,  	&
				"								     ", 		& 
				" The elapsed time was changed accordingly to: ", elapsedtime, " s    ", 		&
				" The value of NSTEPS parameter form input file was discarded.	", 		&
				"*********************************************************************"   
	endif 

end subroutine socket_coupler_init

!=============================================================================
!Apply coupling forces so MD => CFD
!-----------------------------------------------------------------------------
subroutine socket_coupler_apply_continuum_forces(iter)
	use coupler
	use computational_constants_MD, only : delta_t
	use physical_constants_MD, only : np
	use arrays_MD, only :r,v,a
	implicit none
	
	integer, intent(in) :: iter

	integer :: iter_average, Naverage
	real(kind(0.d0)) :: delta_t_CFD
	logical, save :: first_time=.true.
	save Naverage

	if (first_time) then
		first_time = .false.
		Naverage = coupler_md_get_md_steps_per_dt_cfd()
	endif
	iter_average = mod(iter-1, Naverage)+1

!	coupler library version ( deprecated )
!	call coupler_md_apply_continuum_forces(np,r,v,a,iter_average)	!Apply force based on Nie,Chen an Robbins coupling

!  use this module version
	call apply_continuum_forces(iter_average)

	
end subroutine socket_coupler_apply_continuum_forces

!=============================================================================
! Calculate averages of MD to pass to CFD
!-----------------------------------------------------------------------------
subroutine socket_coupler_average(iter)
	use computational_constants_MD, only : initialstep, delta_t
	use physical_constants_MD, only : np
	use arrays_MD, only :r,v
	use coupler
	implicit none

	integer, intent(in) :: iter
	
	integer :: iter_cfd, iter_average, Naverage, save_period, average_period
	logical, save :: first_time=.true.
	save save_period, average_period, Naverage

    if (first_time) then 
	    first_time     = .false.
	    save_period    = coupler_md_get_save_period()    	! period to save uc data (for debugging, testing)
	    average_period = coupler_md_get_average_period() 	! collection interval in the average cycle
	    Naverage = coupler_md_get_md_steps_per_dt_cfd()  	! number of steps in MD average cycle
    endif

    iter_average = mod(iter-1, Naverage)+1	! current step
    
    iter_cfd     = (iter-initialstep)/Naverage +1 ! CFD corresponding step

    if ( mod(iter_average,average_period) .eq. 0 ) then
	    call coupler_md_boundary_cell_average(np,r,v,send_data=.false.) ! accumlate velocities

        ! collect uc data every save_period cfd iteration but discard the first one which cfd uses for initialisation
	    !if ( mod(iter_cfd-1,save_period) .eq. 0 .and. iter_cfd > 1) then
		!    call coupler_uc_average_test(np,r,v,lwrite=.false.)
	    !endif
    endif

	! Send accumulated results to CFD at the end of average cycle 
	if (iter_average .eq. Naverage) then 
	    call coupler_md_boundary_cell_average(np,r,v,send_data=.true.)
	    !if (mod(iter_cfd-1,save_period) .eq. 0 .and. iter_cfd > 1) then
		!    call coupler_uc_average_test(np,r,v,lwrite=.true.)
	    !endif
	endif
	    	
end subroutine socket_coupler_average


!=============================================================================
! Apply force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
!-----------------------------------------------------------------------------
subroutine apply_continuum_forces(iter)
	use coupler
    use physical_constants_MD, only : np
	use computational_constants_MD, only : delta_t, nh, ncells, cellsidelength, halfdomain, delta_rneighbr
	use linked_list
	use arrays_MD, only : r, v, a
    use messenger, only : myid
	implicit none

	integer, intent(in) :: iter ! iteration step, it assumes that each MD average
								! starts from iter = 1

    type(cfd_grid_info) cfd_box
	integer i, j, k, js, je, ib, jb, kb, nib, njb, nkb, ip,m, np_overlap
	integer, allocatable :: list(:,:)
	real(kind=kind(0.d0))  inv_dtCFD,ymin,ymax,xmin,xmax,zmin,zmax,dx_cfd,dz_cfd
	type(node), pointer :: current => null()
	integer, save :: ncalls = 0, itm1=1, itm2=2
    logical, save :: firsttime=.true., overlap
    save nib,njb,nkb,xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD

	! run through the particle, check if they are in the overlap region
	! find the CFD box to which the particle belongs		  
	! attention to the particle that have left the domain boundaries 

	! check if this MD domain overlap with the region in which continuum
	! force is applied
    if (firsttime)then
        firsttime=.false.
        call coupler_md_get(overlap_with_continuum_force=overlap)

        if (overlap) then 
		! get cfd cell sizes
!		call coupler_md_get(ymin_continuum_force=ymin, ymax_continuum_force=ymax, xmin_cfd_grid=xmin, &
!			xmax_cfd_grid=xmax, zmin_cfd_grid=zmin, zmax_cfd_grid=zmax, dx_cfd=dx_cfd, dz_cfd=dz_cfd)
            call coupler_md_get(cfd_box=cfd_box)

            cfd_code_id = coupler_md_get_cfd_id()

            ! number of CFD cells in each direction
            nib = cfd_box%imax - cfd_box%imin 
            njb = 1
            nkb = cfd_box%kmax - cfd_box%kmin 

        ! layer extend in local coordinates, i.e. centered on MD box
            xmin = cfd_box%xmin
            xmax = cfd_box%xmax
            ymin = cfd_box%y(cfd_box%jmax-2)
            ymax = cfd_box%y(cfd_box%jmax-1)
            zmin = cfd_box%zmin
            zmax = cfd_box%zmax
        
            dx_cfd = cfd_box%dx
            dz_cfd = cfd_box%dz
            
            allocate(vel_cfd(3,nib,njb,nkb,2),vbuff(1,nkb,nib+1,njb))
            vel_cfd = 0.d0
            allocate(box_average(nib,njb,nkb))

		! vel_fromCFD cell index from which continum velocity is collected
		!jb_constrain =	  njb - 1 ! the second row of cells from the top
		
        endif
    endif
	
	if (iter .eq. 1) then
        ! get the previous value of CFD velocities
		! this call must be global at the momment
		! because the whole CFD grid is transferred
		! it will be optimised later
		!call coupler_md_get_CFDvel(overlap, vel_cfd)

        call coupler_recv_grid_data(vbuff,index_transpose=(/2,3,1/))

        if ( .not. overlap) return


        !swap the time indices and put the cfd velocities work array
        itm1 = mod(itm1,2)+1
        itm2 = mod(itm2,2)+1
        
        select case (cfd_code_id)
        case (couette_parallel)
            do k=1,nkb
                do j=1,njb
                    do i=1,nib
                        vel_cfd(1,i,j,k,itm1)= 0.5d0*(vbuff(1,k,i,j)+vbuff(1,k,i+1,j))
                    enddo
                enddo
            enddo
        case (couette_serial)
            do k=1,nkb
                do j=1,njb
                    do i=1,nib
                        vel_cfd(1,i,j,k,itm1)= vbuff(1,k,i,j)
                    enddo
                enddo
            enddo
        end select

        !write(800+myid, *) itm1, itm2        
        !do i = 1, nib
        !    write(800+myid,'(1000(E11.3,1x))') vel_cfd(1,i,njb,1:nkb,itm1), vel_cfd(1,i,njb,1:nkb,itm2)
        !enddo
        !write(800+myid,*)
        !call flush(800+myid)

       ! at first CFD step we don't have two values to extrapolate CFD velocities, set inv_dtCFD=0
        if (ncalls .eq. 0) then
            inv_dtCFD = 0.0
        else
            inv_dtCFD = 1.0/coupler_md_get_dt_cfd()
        endif
        ncalls = ncalls + 1	

	endif

	  

	do kb = 1, ubound(box_average,dim=3)
		do jb = 1, ubound(box_average,dim=2)
			do ib = 1, ubound(box_average,dim=1)
				box_average(ib,jb,kb)%np   = 0
				box_average(ib,jb,kb)%v(:) = 0.0d0
				box_average(ib,jb,kb)%a(:) = 0.0d0
			enddo
		enddo
	enddo

	! get the range of j cell index in y direction

	js = min(ncells(2),ceiling((ymin+halfdomain(2))/cellsidelength(2))) + nh
	je = min(ncells(2),ceiling((ymax+halfdomain(2))/cellsidelength(2))) + nh

    !write(0,*) 'md socket ncells(2) js je', ncells(2), js, je

!	find the maximum number of molecules and allocate a list array	   
	np_overlap = 0 ! number of particles in overlapping reg

    allocate(list(4,np))

	!do k = nh + 1, nh + 1 + ncells(3)
	!	do j = js, je
	!		do i = nh + 1, nh + 1 + ncells(1) 
	!			np_overlap =np_overlap + cell%cellnp(i,j,k)
	!		enddo
	!	enddo
	!enddo

	!allocate(list(4,np_overlap))

    !write(0,*)" md socket, np_overlap", np_overlap

	np_overlap = 0

    do m = 1,np
       if ( r(m,2) >  ymin            .and. r(m,2) < ymax          .and. &
           r(m,1) >= -halfdomain(1)  .and. r(m,1) < halfdomain(1) .and. &
           r(m,3) >= -halfdomain(3)  .and. r(m,3) < halfdomain(3) ) then
           ib = ceiling( (r(m,1) -	 xmin) / dx_cfd)
           jb = 1
           kb = ceiling( (r(m,3) - zmin ) / dz_cfd)
           
           np_overlap = np_overlap + 1
           list(1:4, np_overlap) = (/ m, ib, jb, kb /)
           
           box_average(ib,jb,kb)%np   =  box_average(ib,jb,kb)%np + 1
           box_average(ib,jb,kb)%v(:) =  box_average(ib,jb,kb)%v(:) + v(m,:)
           box_average(ib,jb,kb)%a(:) =  box_average(ib,jb,kb)%a(:) + a(m,:)
       endif
   enddo

!      The cell base version does not work 
!      Revise !!
!
!!$    do k = nh + 1, nh + 1 + ncells(3)
!!$		do j = js, je
!!$			do i = nh + 1, nh + 1 + ncells(1) 
!!$
!!$				if ( cell%cellnp(i,j,k) == 0) cycle
!!$			 
!!$                current => cell%head(i,j,k)%point
!!$				
!!$				do ip = 1, cell%cellnp(i,j,k) 
!!$
!!$					m = current%molno
!!$					
!!$					! get the CFD cell coordinates	
!!$                    ! apply the force only for the particle which are in domain
!!$                    ! exchange of particles outside of domain to be added later 
!!$					if ( r(m,2) >  ymin            .and. r(m,2) < ymax          .and. &
!!$                         r(m,1) >= -halfdomain(1)  .and. r(m,1) < halfdomain(1) .and. &
!!$                         r(m,3) >= -halfdomain(3)  .and. r(m,3) < halfdomain(3) ) then
!!$                        ib = ceiling( (r(m,1) -	 xmin) / dx_cfd)
!!$                        jb = 1
!!$                        kb = ceiling( (r(m,3) - zmin ) / dz_cfd)
!!$
!!$                        np_overlap = np_overlap + 1
!!$                        list(1:4, np_overlap) = (/ m, ib, jb, kb /)
!!$                        
!!$                        box_average(ib,jb,kb)%np   =  box_average(ib,jb,kb)%np + 1
!!$                        box_average(ib,jb,kb)%v(:) =  box_average(ib,jb,kb)%v(:) + v(ip,:)
!!$                        box_average(ib,jb,kb)%a(:) =  box_average(ib,jb,kb)%a(:) + a(ip,:)
!!$                    endif
!!$				 
!!$					current => current%next
!!$				enddo
!!$			enddo
!!$		enddo
!!$	enddo
	
	! here we should have the cell coordinates for the particle ip which is 
	! in the overlap region
	! one has to treat separatley the particle that have left the domain
	! compute the average force for each bin

	!write(0,*)'MD before average over bin. np_overlap', np_overlap

	call average_over_bin

	!write(0,*) 'MD: end simulation_apply_continuum_forces', myid

contains

!=============================================================================
! Get velocity from CFD and apply force to molecule
!-----------------------------------------------------------------------------

	subroutine average_over_bin
		implicit none

		integer ib, jb, kb, i, ip, n
		real(kind=kind(0.d0)) alpha(3), u_cfd_t_plus_dt(3), inv_dtMD, acfd


		! set the continnum constraints for the particle in the bin
		! speed extrapolation 
		! add all up
		inv_dtMD =1.d0/delta_t

		!write(0,'(a,I7,2E12.4)') "MD continuum np, vel_fromCFD1 : ", np_overlap, &
		!						  maxval(a(list(1,1:np_overlap),:)), &
		!						  minval(a(list(1,1:np_overlap),:))

		do i = 1, np_overlap
			ip = list(1,i)
			ib = list(2,i)
			jb = list(3,i)
			kb = list(4,i)

			n = box_average(ib,jb,kb)%np

			!write(900+myid,'(a,5I4,14E12.4)') "MD continuum force", iter,ib,jb,kb,n,box_average(ib,jb,kb)%v(1), &
			!	box_average(ib,jb,kb)%a(1),v(ip,1),a(ip,1),inv_dtMD,inv_dtCFD

			if ( n .eq. 0 ) cycle

			! using the following exptrapolation formula for continuum velocity
			! y = (y2-y1)/(x2-x1) * (x-x2) +y2
            alpha(1) = inv_dtCFD*(vel_cfd(1,ib,1,kb,itm1) - &
                    vel_cfd(1,ib,1,kb,itm2))

			u_cfd_t_plus_dt(1) = alpha(1) * (iter + 1)*delta_t + vel_cfd(1,ib,1,kb,itm1) 

			acfd =	- box_average(ib,jb,kb)%a(1) / n - inv_dtMD * & 
				( box_average(ib,jb,kb)%v(1) / n - u_cfd_t_plus_dt(1) )
			a(ip,1) = a(ip,1) + acfd

			!	write(900+myid,'(a,4I4,15E12.4)') "MD continuum force 2", ib,jb,kb,n, &
			!	 alpha(1),u_cfd_t_plus_dt(1),vel_cfd(1,ib,1,kb,itm1),&
			!	 vel_cfd(1,ib,1,kb,itm2), a(ip,1),acfd, r(ip,2) 

		enddo


		!	write(400+10*ncalls+myid,'(a,I7,2E12.4)') "MD continuum np, vel_fromCFD 2: ", np_overlap, &
		!						   maxval(a(list(1,1:np_overlap),:)), &
		!						   minval(a(list(1,1:np_overlap),:))
		!	write(400+10*ncalls+myid,'(a,2E12.4)')" inv_dtCFD, inv_dtMD ", inv_dtCFD, inv_dtMD
		!	do kb=1,nkb
		!	do jb=1,njb
		!	do ib=1,nib
		!		write(400+10*ncalls+myid,'(12E12.4,I7)') vel_fromCFD(:,ib,jb,kb,1), vel_fromCFD(:,ib,jb,kb,2),&
		!						  box_average(ib,jb,kb)%v(:), box_average(ib,jb,kb)%a(:),&
		!						  box_average(ib,jb,kb)%np
		!	enddo
		!	enddo
		!	enddo

	end subroutine average_over_bin

end subroutine apply_continuum_forces

!=============================================================================
! Apply force from Nie et al (2004) paper to fix molecular velocity to
! continuum value inside the overlap region. 
!-----------------------------------------------------------------------------

subroutine apply_continuum_forces_ES(iter)
	use coupler
	use computational_constants_MD, only : delta_t,nh,halfdomain,ncells,cellsidelength,initialstep,Nsteps
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	implicit none

	integer, intent(in) 				:: iter ! iteration step, it assumes that each MD average
												! starts from iter = 1

    type(cfd_grid_info)                 :: cfd_box
	integer								:: i,j,k,js,je,ib,jb,kb,nib,njb,nkb,ip,m,np_overlap
	integer         					:: n, molno, ixyz, cellnp
	integer         					:: cbin, averagecount
	integer								:: icell, jcell, kcell
	integer								:: isummol
	double precision					::  inv_dtCFD,ymin,ymax,xmin,xmax,zmin,zmax,dx_cfd,dy_cfd,dz_cfd
	double precision					:: isumvel, isumacc
	double precision					:: t_fract, average
	double precision, dimension(4)		:: continuum_u
    integer, save                       :: ncalls = 0, itm1=1, itm2=2
	logical, save						:: firsttime, overlap
	type(node), pointer 	        	:: old, current
    save nib,njb,nkb,xmin,xmax,ymin,ymax,zmin,zmax,dx_cfd,dz_cfd,inv_dtCFD

	! run through the particle, check if they are in the overlap region
	! find the CFD box to which the particle belongs		  
	! attention to the particle that have left the domain boundaries 

	! check if this MD domain overlap with the region in which continuum
	! force is applied
    if (firsttime)then
        firsttime=.false.
        call coupler_md_get(overlap_with_continuum_force=overlap)

        if (overlap) then 
		! get cfd cell sizes
		!call coupler_md_get(ymin_continuum_force=ymin, ymax_continuum_force=ymax, xmin_cfd_grid=xmin, &
		!	xmax_cfd_grid=xmax, zmin_cfd_grid=zmin, zmax_cfd_grid=zmax, dx_cfd=dx_cfd, dz_cfd=dz_cfd)


            call coupler_md_get(cfd_box=cfd_box)
            cfd_code_id = coupler_md_get_cfd_id()

		! number of CFD cells in each direction
            nib = cfd_box%imax - cfd_box%imin
            njb = 1
            nkb = cfd_box%kmax - cfd_box%kmin
		
            ! layer extend in local coordinates, i.e. centered on MD box
            xmin = cfd_box%xmin
            xmax = cfd_box%xmax
            ymin = cfd_box%y(cfd_box%jmax-2)
            ymax = cfd_box%y(cfd_box%jmax-1)
            zmin = cfd_box%zmin
            zmax = cfd_box%zmax
        
            dx_cfd = cfd_box%dx
            dz_cfd = cfd_box%dz

            allocate(vel_cfd(1,nib,njb,nkb,2),vbuff(1,nkb,nib+1,njb))
            vel_cfd = 0.d0
            ! vel_fromCFD cell index from which continum velocity is collected
            !jb_constrain =	  njb - 1 ! the second row of cells from the top
		
        endif
    endif
	
	if (iter .eq. 1) then
		! get the previous value of CFD velocities
		! this call must be global at the moment
		! because the whole CFD grid is transferred
		! it will be optimised later
        call coupler_recv_grid_data(vbuff,index_transpose=(/2,3,1/))
	
        if ( .not. overlap) return

        select case (cfd_code_id)
        case (couette_parallel)
            do k=1,nkb
                do j=1,njb
                    do i=1,nib
                        vel_cfd(1,i,j,k,itm1)= 0.5d0*(vbuff(1,k,i,j)+vbuff(1,k,i+1,j))
                    enddo
                enddo
            enddo
        case (couette_serial)
            do k=1,nkb
                do j=1,njb
                    do i=1,nib
                        vel_cfd(1,i,j,k,itm1)= vbuff(1,k,i,j)
                    enddo
                enddo
            enddo
        end select

        ! at first CFD step we don't have two values to extrapolate CFD velocities, set inv_dtCFD=0
        if (ncalls .eq. 0) then
            inv_dtCFD = 0.0
        else
            inv_dtCFD = 1.0/coupler_md_get_dt_cfd()
        endif
        ncalls = ncalls + 1

    endif

	! get the range of j cell index in y direction
	js = min(ncells(2),ceiling(ymin/cellsidelength(2))) + nh
	je = min(ncells(2),ceiling(ymax/cellsidelength(2))) + nh

	!Fraction of continuum timestep which has passed
	t_fract = dble((iter - initialstep)) / dble((Nsteps-initialstep))

	!Linear extrapolation between velocity at t and t+1 and save in continuum_u
	!Taken for cell at nx/2 and top domain cell (ny + 1) to (ny + 1) - overlap
	!N.B. domain runs from 2 to ny + 1 due to halos
	!continuum_u(1) = uc_t_minus_1(nint(nx/2.d0),3)*(1.d0-t_fract)  &
	!	 	         + uc(nint(nx/2.d0),3)*      t_fract
	!continuum_u(2) = uc_t_minus_1(nint(nx/2.d0),3)*(1.d0-t_fract)  &
	!	 	         + uc(nint(nx/2.d0),3)*      t_fract
	!continuum_u(3) = uc_t_minus_1(nint(nx/2.d0),4)*(1.d0-t_fract)  &
	!	 	         + uc(nint(nx/2.d0),4)*      t_fract
	!continuum_u(4) = uc_t_minus_1(nint(nx/2.d0),4)*(1.d0-t_fract)  &
	!	 	         + uc(nint(nx/2.d0),4)*      t_fract

	average = 0.d0
	averagecount = 0

	!Apply force to top three bins in y
	!ASSUME Cell same size as bins and one continuum cell is two MD cells
	!do jcell= (ncells(2)+1)-3,(ncells(2)+1)			!Loop through 4 y cells in controlled region
	do jcell= js,je
 		cbin = jcell - (ncells(2)+1)+4			!Local no. of overlap cell from 1 to overlap

		!print'(3i8,4f20.7)', iter, jcell, cbin, continuum_u

		!Reset acceleration and velocity sums
		isummol = 0
		isumvel = 0.d0
		isumacc = 0.d0

		!do icell=2 , ncells(1)+1	!Loop through all x cells
		!do kcell=2 , ncells(3)+1	!Loop through all z cells
    	do kcell = nh+1,ncells(3)+1+nh
		do icell = nh+1,ncells(1)+1+nh

			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list

			!Calculate averages for bin
			do n = 1, cellnp    ! Loop over all particles
				molno = old%molno 	 !Number of molecule
				!Assign to bins using integer division
				isumvel = isumvel + v(molno,1) 	!Add streamwise velocity to current bin
				isumacc = isumacc + a(molno,1) 	!Add acceleration to current bin
				isummol = isummol + 1
				current => old
				old => current%next 
			enddo

		enddo
		enddo

		!Get average velocity and acceleration in bin
		if (isummol .ne. 0) then
			isumacc = isumacc/real(isummol,kind(0.d0))
		 	isumvel = isumvel/real(isummol,kind(0.d0))
		endif

		!print*, isumacc, isumvel, isummol

		!do icell=2 , ncells(1)+1	!Loop through all x cells
		!do kcell=2 , ncells(3)+1	!Loop through all z cells
    	do kcell = nh+1,ncells(3)+1+nh
		do icell = nh+1,ncells(1)+1+nh

			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list
			
			!Apply coupling force as Nie, Chen and Robbins (2004), using
			!Linear extrapolation of velocity
			do n = 1, cellnp    ! Loop over all particles
				molno = old%molno !Number of molecule

				a(molno,1)= a(molno,1) - isumacc   &
					    -(isumvel-vel_cfd(1,icell,1,kcell,1))/delta_t

				current => old
				old => current%next 

				!if (molno .gt. np) stop "Force applied to halo molecules"

				!average = average - isumacc   &
				!	    -(isumvel-continuum_u(cbin))/delta_t
				!averagecount = averagecount + 1

			enddo

		enddo
		enddo
	enddo

	!print'(a,f10.5,a,f18.9)', 'MD_velocity ',sum(continuum_u(:))/3 , & 
	!			' average force applied to MD molecules ', average/(averagecount)

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine apply_continuum_forces_ES

!--------------------------------------------------------------------------------------
! Apply force to match velocity to continuum using CV formulation of Nie et al 2004
! which results in matching of force and flux on a CV.

subroutine simulation_apply_continuum_forces_CV(iter)
	use computational_constants_MD, only : delta_t,nh,domain,halfdomain,ncells,cellsidelength,initialstep,Nsteps
	use calculated_properties_MD, only : nbins
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	implicit none

	integer, intent(in) 				:: iter ! iteration step, it assumes that each MD average
												! starts from iter = 1

	integer								:: i,j,k,js,je,ib,jb,kb,nib,njb,nkb,ip,m,np_overlap
	integer         					:: n, molno, ixyz, cellnp
	integer         					:: cbin, averagecount
	integer								:: icell, jcell, kcell
	integer								:: isummol
	integer, save 						:: ncalls = 0
	double precision					:: inv_dtCFD,ymin,ymax,xmin,xmax,zmin,zmax,dx_cfd,dy_cfd,dz_cfd
	double precision					:: isumflux, isumforce, isumvel
	double precision					:: t_fract, average
	double precision, dimension(3)		:: ri
	double precision, dimension(4)		:: continuum_res, continuum_Fs
	double precision, dimension(4)		:: continuum_u
	double precision, dimension(:,:,:,:,:),allocatable 	:: vel_cfd

	logical								:: overlap
	type(node), pointer 	        	:: old, current

	!Linear extrapolation between force at t and t+1 and save in continuum_F
	!Taken for cell at nx/2 and top domain cell (ny + 1) to (ny + 1)-overlap
	!N.B. domain runs from 2 to ny + 1 due to halos
	!continuum_res(1) = xresidual_t_minus_1(nint(nx/2.d0),3)*(1.d0-t_fract)  &
	!	 	           + xresidual(nint(nx/2.d0),3)*      t_fract
	!continuum_res(2) = xresidual_t_minus_1(nint(nx/2.d0),3)*(1.d0-t_fract)  &
	!	 	           + xresidual(nint(nx/2.d0),3)*      t_fract
	!continuum_res(3) = xresidual_t_minus_1(nint(nx/2.d0),4)*(1.d0-t_fract)  &
	!	 	           + xresidual(nint(nx/2.d0),4)*      t_fract
	!continuum_res(4) = xresidual_t_minus_1(nint(nx/2.d0),4)*(1.d0-t_fract)  &
	!	 	           + xresidual(nint(nx/2.d0),4)*      t_fract

	!Fraction of continuum timestep which has passed
	t_fract = dble((iter - initialstep)) / dble((Nsteps-initialstep))

	!Apply force to top three bins in y
	!ASSUME Cell same size as bins and one continuum cell is two MD cells
	do jcell= (ncells(2)+1)-3,(ncells(2)+1)	!Loop through 4 y cells in controlled region
 	cbin = jcell - (ncells(2)+1)+4		!Local no. of overlap cell from 1 to overlap
	do icell=2 , ncells(1)+1		!Loop through all x cells
	do kcell=2 , ncells(3)+1		!Loop through all z cells

		!Reset flux and force sums
		isummol   = 0

		!Calculate flux and force averages for bin
		call compute_bin_surface_flux(icell,jcell,kcell,isumflux)
		call compute_force_surrounding_bins(icell,jcell,kcell,isumforce)

		cellnp = cell%cellnp(icell,jcell,kcell)
		old => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list

		isumvel = 0.d0

		!Calculate velocity and molecular averages for bin
		do n = 1, cellnp    ! Loop over all particles

			molno = old%molno	!Number of molecule
			isumvel = isumvel + v(molno,1) 	!Add streamwise velocity to current bin
			isummol = isummol + 1

			current => old
			old => current%next !Use pointer in datatype to obtain next item in list
		enddo

		!All constraint terms are multiplied by { (m_j \theta_j)/(\sum_i^N m_i^2 \theta_i^2) }
		!Which for unit mass reduces to { 1/(sum inside volume) }
		if (isummol .ne. 0) then
			isumflux     = isumflux     / real(isummol,kind(0.d0))
		 	isumforce    = isumforce    / real(isummol,kind(0.d0))
			continuum_Fs = continuum_res/ real(isummol,kind(0.d0))
			if(icell .eq. 5 .and. kcell .eq. 5 .and. jcell .eq. 11 ) then 
				!print'(a,i8,4f10.5)','bin and forces',cbin,continuum_Fs
				write(99,'(3i8,4f18.9)')iter,jcell,isummol, isumflux, isumforce, isumvel,continuum_Fs(cbin)
			endif
		endif

		old => cell%head(icell,jcell,kcell)%point 	!Reset old to first molecule in list
		
		!Apply coupling force for CV form of Nie, Chen and Robbins (2004)
		do n = 1, cellnp    ! Loop over all particles
			molno = old%molno !Number of molecule

			a(molno,1) = a(molno,1) + ((isumflux - isumforce) + continuum_Fs(cbin))

			current => old
			old => current%next 

			!if (molno .gt. np) stop "Force applied to halo molecules"
		enddo

	enddo
	enddo
	enddo

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

end subroutine simulation_apply_continuum_forces_CV

!===================================================================================
! Flux of molecules over bin surface requires all molecules in current cell and 
! surrounding cells to be checked

subroutine compute_bin_surface_flux(icell,jcell,kcell,isumflux)
	use computational_constants_MD, only : delta_t,nh,domain,halfdomain,ncells,cellsidelength,initialstep,Nsteps
	use calculated_properties_MD, only : nbins
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	implicit none

	integer,intent(in)					:: icell,jcell,kcell
	double precision,intent(out)		:: isumflux

	integer								:: n,molno
	integer								:: icellshift,jcellshift,kcellshift,adjacentcellnp
	integer,dimension(3)				:: ibin1, ibin2
	double precision,dimension(3)		:: Fbinsize,ri1,ri2,ri12,Fsurface,velvec
	type(node), pointer		 	        :: old, current

	isumflux = 0.d0

	!Determine bin size
	Fbinsize(:) = domain(:) / nbins(:)

	!Calculate bin surface Forces
	do kcellshift = -1,1
	do jcellshift = -1,1
	do icellshift = -1,1
		old =>  cell%head(icell+icellshift, & 
				  jcell+jcellshift, & 
				  kcell+kcellshift)%point
		adjacentcellnp = cell%cellnp(icell+icellshift, & 
					     jcell+jcellshift, & 
					     kcell+kcellshift)

		do n = 1,adjacentcellnp          !Step through all adjacent cells' molecules

			molno = old%molno			!Number of molecule

			velvec(:) = v(molno,:) + delta_t *a(molno,:) 	!Velocity at t calculated from acceleration
			ri1(:)    = r(molno,:) + delta_t * velvec(:)	!Position at t calculated from velocity
			ri2(:)    = r(molno,:)				!Molecule i at time t-dt

			current => old
			old => current%next    !Use pointer in datatype to obtain next item in list

			!Check if molecule crosses surface, if it has, 
			!add to surface flux bin total

			Fsurface = 0.d0
			! *********************************************************************************
			!Calculate flux over surface only if molecule is entering/leaving bin of interest
			!ibin1(:) = ceiling((ri1(:)+halfdomain(:))/Fbinsize(:))+1	!Establish current bin
			!ibin2(:) = ceiling((ri2(:)+halfdomain(:))/Fbinsize(:))+1 !Establish current bin

			!if (icell.eq.ibin1(1).and.jcell.eq.ibin1(2).and.kcell.eq.ibin1(3).or.&
			!    icell.eq.ibin2(1).and.jcell.eq.ibin2(2).and.kcell.eq.ibin2(3)) then
			!	call get_Fsurface(icell,jcell,kcell,molno,molno,velvec,ri1,ri2,Fsurface)
			!	if(icell .eq. 5 .and. kcell .eq. 5 .and. jcell .eq. 11 ) then 
			!!		if (any(Fsurface .ne. 0.d0)) print'(2i8,6f10.5)', iter,molno, & 
			!!			dble(ibin1(:)-ibin2(:)),velvec
			!!	endif
			!endif
			! *********************************************************************************
			call get_Flux(icell,jcell,kcell,molno,molno,velvec,ri1,ri2,Fsurface)

			!if(icell .eq. 5 .and. kcell .eq. 5 ) then 
			!	if (any(Fsurface .ne. 0.d0)) print'(5i8,6f10.5)', iter,icell,jcell,kcell,molno,velvec,Fsurface
			!endif
			isumflux = isumflux + Fsurface(1)

		enddo
	enddo
	enddo
	enddo

contains

	!===================================================================================
	!Flux over the surface of a bin

	subroutine get_Flux(icell,jcell,kcell,molnoi,molnoj,velvect,ri1,ri2,Fsurface)
		use coupler
		use librarymod, only : heaviside
		implicit none

		integer,intent(in)							:: icell,jcell,kcell,molnoi,molnoj
		double precision,dimension(3),intent(in)	:: ri1, ri2, velvect
		double precision,dimension(3),intent(out)	:: Fsurface

		integer										:: ixyz,jxyz,kxyz,i,j,k,n
		integer										:: planeno
		integer										:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
		integer		,dimension(1)					:: imaxloc
		integer		,dimension(3)					:: ibin1,ibin2,cbin
		double precision							:: crosstime,crossplane,rplane,shift
		double precision,dimension(3)				:: mbinsize,crossface
		double precision,dimension(3)				:: ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
		double precision,dimension(3,3,3,3)			:: Fsurfacebins

		!CV momentum flux
		!Determine bin size
		mbinsize(:) = domain(:) / nbins(:)

		ri12   = ri1 - ri2		!Molecule i trajectory between t-dt and t
		where (ri12 .eq. 0.d0) ri12 = 0.000001d0

		!Assign to bins before and after using integer division
		ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:))+1
		ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:))+1

		!Replace Signum function with this functions which gives a
		!check for plane crossing and the correct sign 
		crossface(:) =  ibin1(:) - ibin2(:)

		!if(icell .eq. 5 .and. kcell .eq. 5 ) then 
		!	if (any(crossface .eq. 1)) print*, ibin1(:),ibin2(:)
		!endif

		if (sum(abs(crossface(:))) .ne. 0) then

			Fsurfacebins = 0.d0

			do i = ibin1(1),ibin2(1),sign(1,ibin2(1)-ibin1(1))
			do j = ibin1(2),ibin2(2),sign(1,ibin2(2)-ibin1(2))
			do k = ibin1(3),ibin2(3),sign(1,ibin2(3)-ibin1(3))

				cbin(1) = i; cbin(2) = j; cbin(3) = k

				bintop(:) = (cbin(:)-1)*mbinsize(:)-halfdomain(:)
				binbot(:) = (cbin(:)-2)*mbinsize(:)-halfdomain(:)

				!Calculate the plane intersect of trajectory with surfaces of the cube
				Pxt=(/ 							  bintop(1),		 & 
						ri1(2)+(ri12(2)/ri12(1))*(bintop(1)-ri1(1)), & 
						ri1(3)+(ri12(3)/ri12(1))*(bintop(1)-ri1(1))  	/)
				Pxb=(/ 							  binbot(1),		 & 
						ri1(2)+(ri12(2)/ri12(1))*(binbot(1)-ri1(1)), & 
						ri1(3)+(ri12(3)/ri12(1))*(binbot(1)-ri1(1))  	/)
				Pyt=(/	ri1(1)+(ri12(1)/ri12(2))*(bintop(2)-ri1(2)), & 
												  bintop(2),		 & 
						ri1(3)+(ri12(3)/ri12(2))*(bintop(2)-ri1(2))  	/)
				Pyb=(/	ri1(1)+(ri12(1)/ri12(2))*(binbot(2)-ri1(2)), &
												  binbot(2),		 & 
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


				!Add Momentum flux over face
				Fsurfacebins(modulo(cbin(1),3)+1, 	& 
					     modulo(cbin(2),3)+1, 	& 
					     modulo(cbin(3),3)+1,:) = 	& 
						 Fsurfacebins(modulo(cbin(1),3)+1, 	& 
							      modulo(cbin(2),3)+1, 	& 
							      modulo(cbin(3),3)+1,:) 	& 
							 		- velvect(:)*dble(onfacexb - onfacext) & 
							 		- velvect(:)*dble(onfaceyb - onfaceyt) & 
							 		- velvect(:)*dble(onfacezb - onfacezt)

				!if(icell .eq. 5 .and. kcell .eq. 5 ) then 
				!	print'(12i8)', cbin, ibin1,ibin2,icell,jcell,kcell
				!	print'(12i8)', modulo(cbin(:),3)+1,modulo(ibin1(:),3)+1,modulo(ibin2(:),3)+1 & 
				!	      	      ,modulo(icell,3)+1,modulo(jcell,3)+1,modulo(kcell,3)+1
				!	print'(6f10.5,6i8)', velvect(:),Fsurface,onfacexb,onfacext,onfaceyb,onfaceyt,onfacezb,onfacezt
				!endif

			enddo
			enddo
			enddo

			!if(icell .eq. 5 .and. kcell .eq. 5 ) then 
			!	print'(6f10.5,6i8)', velvect(:),Fsurface,onfacexb,onfacext,onfaceyb,onfaceyt,onfacezb,onfacezt
			!endif
			

			Fsurface(:) = Fsurfacebins(modulo(icell,3)+1, 	& 
				     		   modulo(jcell,3)+1, 	& 
				    		   modulo(kcell,3)+1,:)

		endif

	end subroutine get_Flux

end subroutine compute_bin_surface_flux

!===================================================================================
! Forces between current bin and surrounding bins

subroutine compute_force_surrounding_bins(icell,jcell,kcell,isumforce)
	use physical_constants_MD, only : nd, np, rcutoff2
	use computational_constants_MD, only : domain,halfdomain
	use calculated_properties_MD, only : nbins
	use arrays_MD, only : r, v, a
	use linked_list, only : node, cell
	use coupler
	implicit none

	integer,intent(in)					:: icell,jcell,kcell
	double precision,intent(out)		:: isumforce

	integer								:: n,j,ixyz,molnoi,molnoj
	integer								:: icellshift,jcellshift,kcellshift,cellnp,adjacentcellnp
	integer,dimension(3)				:: ibin, jbin
	double precision					:: rij2, invrij2, accijmag
	double precision,dimension(3)		:: ri,rj,rij,fij,Fsurface
	type(node), pointer		 	        :: oldi, currenti, oldj, currentj

	!print'(a,4i8,4f10.5)', 'Before input', icell,jcell,kcell,molnoi,ri(:),isumforce

	isumforce = 0.d0

	cellnp = cell%cellnp(icell,jcell,kcell)
	oldi => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list

	!Calculate averages for bin
	do n = 1, cellnp    ! Loop over all particles

		molnoi = oldi%molno	!Number of molecule
		ri = r(molnoi,:)	!Retrieve ri

		!Calculate bin surface Forces
		do kcellshift = -1,1
		do jcellshift = -1,1
		do icellshift = -1,1
			oldj => cell%head(icell+icellshift, & 
					  jcell+jcellshift, & 
					  kcell+kcellshift)%point
			adjacentcellnp = cell%cellnp(icell+icellshift, & 
						     jcell+jcellshift, & 
						     kcell+kcellshift)

			do j = 1,adjacentcellnp          !Step through all j for each i

				molnoj = oldj%molno 	 !Number of molecule
				rj = r(molnoj,:)         !Retrieve rj

				currentj => oldj
				oldj => currentj%next    !Use pointer in datatype to obtain next item in list

				if(molnoi==molnoj) cycle !Check to prevent interaction with self

				rij2=0                   !Set rij^2 to zero
				rij(:) = ri(:) - rj(:)   !Evaluate distance between particle i and j

				!rij2 = dot_product(rij)
				do ixyz=1,nd
					rij2 = rij2+rij(ixyz)*rij(ixyz) !Square of vector calculated
				enddo

				if (rij2 < rcutoff2) then

					invrij2 = 1.d0/rij2                 !Invert value

					!Linear magnitude of acceleration for each molecule
					accijmag = 48.d0*(invrij2**7-0.5d0*invrij2**4)

					!Get force and add to bin total
					fij = accijmag*rij(:)
					Fsurface = 0.d0
					!print*, 'FORCE', fij
					call get_Fsurface(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface)
					isumforce = isumforce +  Fsurface(1)

				endif
			enddo
		enddo
		enddo
		enddo

		currenti => oldi
		oldi => currenti%next !Use pointer in datatype to obtain next item in list
	enddo

contains

	!===================================================================================
	!Forces over the surface of a bin

	subroutine get_Fsurface(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface)
		use coupler
		use librarymod, only : heaviside
		implicit none

		integer,intent(in)				:: icell,jcell,kcell,molnoi,molnoj
		double precision,dimension(3),intent(in)	:: ri, rj, fij
		double precision,dimension(3),intent(out)	:: Fsurface

		integer									:: ixyz
		integer,dimension(3)					:: ibin, jbin
		double precision,dimension(3)			:: Fbinsize, bintopi, binboti, bintopj, binbotj, crossplane

		!Determine bin size
		Fbinsize(:) = domain(:) / nbins(:)

		!Assign to bins using integer division
		ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+1	!Establish current bin
		jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+1 	!Establish current bin

		bintopi(:) = (ibin(:)-1)*Fbinsize(:)-halfdomain(:)
		binboti(:) = (ibin(:)-2)*Fbinsize(:)-halfdomain(:)
		bintopj(:) = (jbin(:)-1)*Fbinsize(:)-halfdomain(:)
		binbotj(:) = (jbin(:)-2)*Fbinsize(:)-halfdomain(:)

		!Add for molecule i
		if(molnoi .le. np) then
			Fsurface = fij(:)* dble((heaviside(bintopi(1)-ri(1))-heaviside(binboti(1)-ri(1)))* & 
							  		(heaviside(bintopi(2)-ri(2))-heaviside(binboti(2)-ri(2)))* & 
							  		(heaviside(bintopi(3)-ri(3))-heaviside(binboti(3)-ri(3)))- & 
							  		(heaviside(bintopi(1)-rj(1))-heaviside(binboti(1)-rj(1)))* & 
							  		(heaviside(bintopi(2)-rj(2))-heaviside(binboti(2)-rj(2)))* & 
							  		(heaviside(bintopi(3)-rj(3))-heaviside(binboti(3)-rj(3))))
		endif

		!Add for molecule j
		if(molnoj .le. np) then
			Fsurface = fij(:)* dble((heaviside(bintopj(1)-ri(1))-heaviside(binbotj(1)-ri(1)))* & 
							  		(heaviside(bintopj(2)-ri(2))-heaviside(binbotj(2)-ri(2)))* & 
							  		(heaviside(bintopj(3)-ri(3))-heaviside(binbotj(3)-ri(3)))- & 
							  		(heaviside(bintopj(1)-rj(1))-heaviside(binbotj(1)-rj(1)))* & 
							  		(heaviside(bintopj(2)-rj(2))-heaviside(binbotj(2)-rj(2)))* & 
							  		(heaviside(bintopj(3)-rj(3))-heaviside(binbotj(3)-rj(3))))
		endif

	end subroutine get_Fsurface


end subroutine compute_force_surrounding_bins

#endif

end module md_coupler_socket
