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

 	integer nsteps_cfd, naverage
	real(kind(0.d0)) :: delta_t_CFD

	! if coupled calculation prepare exchange layout
	call coupler_md_init(npx,npy,npz,icoord,icomm_grid,delta_t)

	! fix NSTEPS for the coupled case
	delta_t_CFD = coupler_md_get_dt_cfd()
	nsteps_cfd  = coupler_md_get_nsteps()

    naverage = int(delta_t_cfd/delta_t)

    if (naverage == 0) then 
        call error_abort("ratio dt_cfd/dt_md < 1, coupler will not work !!!")
    endif

	Nsteps = initialstep + nsteps_cfd * naverage
	elapsedtime = elapsedtime + nsteps_cfd * delta_t_CFD

	if (myid .eq. 0) then 
		write(*,'(4(a,/),a,I7,/a,/a,E10.4,a,/,a,/a)') &
				"*********************************************************************", 		&
 				"WARNING - this is a coupled run which resets the number	      ", 		&
				" of extrasteps to:						   ", 		&
				"								     ", 		&
				" nstep_cfd*(delta_t_cdf/delta_t_md) = ",nsteps_cfd*naverage            , 	&
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
		Naverage = int(coupler_md_get_dt_cfd()/delta_t)
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
	    save_period    = coupler_md_get_save_period()    ! period to save uc data (for debugging, testing)
	    average_period = coupler_md_get_average_period() ! collection interval in the average cycle
	    Naverage = int(coupler_md_get_dt_cfd() /delta_t)	   ! number of steps in MD average cycle
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


#endif

end module md_coupler_socket
