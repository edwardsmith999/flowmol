!=======================================================================
! Serial Message-passing routines to emulate parallel code
!
!===============================
!  Key Messenger Subroutines   =
!===============================
! messenger_invoke()
! messenger_init()
! messenger_syncall()
! messenger_free()
!===============================
!  Border Update Subroutines   =
!===============================
! messenger_updateborders()
! updateface(ixyz)
! updateedge(face1,face2)
! updatecorners()
!===============================
! Molecule Sending Subroutines =messenger_updateborders
!===============================
! sendmols()
! sendface(ixyz)
! sendedge(face1,face2)
! sendcorners()
!!===============================
! Send-Probe-Receive Subroutines=
!================================
! sendproberecv(recvsize,sendsize,sendbuffer,pos,isource,idest)
! pairedsendproberecv(recvsize,sendsize,sendbuffer,pos,isource,idest,pairsplit)
! NBsendproberecv(recvsize,sendsize,sendbuffer,pos,length,isource,idest)
!===============================
! Data gathering subroutines   =
!===============================
! globalsyncreduce(A, na, meanA, maxA, minA)
! globalSum(A, na)
! globalMax(A, na)
! globalMin(A, na)
! globalAverage(A, na)
!

module messenger

	use physical_constants_MD
	use computational_constants_MD
	use linked_list
	use polymer_info_MD	
	use shear_info_MD

        integer :: MD_COMM                     	! global communicator
	integer :: myid                         ! my process rank
	integer :: idroot                       ! rank of root process

	! Grid topology
	integer :: icomm_grid                   ! comm for grid topology
	integer, allocatable :: icoord(:,:)     ! proc grid coordinates
	integer	:: icomm_xyz(3) 		! Directional subcomms

end module

!======================================================================
!			Key Messenger Subroutines                     =
!======================================================================
subroutine messenger_invoke()
	use messenger

        npx = 1; npy = 1;  npz = 1;  nproc = 1
	print*, 'Serial version of parallel code'

	return
end


subroutine messenger_init()
	use messenger
	implicit none

	if (nproc .ne. 1) stop "Serial code - Param.inc should be 1x1x1"

	irank  = 1
	iblock = 1
	jblock = 1
	kblock = 1
	iroot  = 1

        ! allocate arrays that depend on topology parameters

        allocate(ibmin(npx), ibmax(npx), ibmino(npx), ibmaxo(npx), &
		jbmin(npy), jbmax(npy), jbmino(npy), jbmaxo(npy), &
		kbmin(npz), kbmax(npz), kbmino(npz), kbmaxo(npz), stat=ierr)
        if (ierr /= 0) then 
                write(*,*) "Error allocating topology arrays in messenger_init"
                stop
        endif

        allocate(icoord(3,nproc),stat=ierr)
        if (ierr /= 0) then 
                write(*,*) "Error allocating icoord in messenger_init"
                stop
        endif

	return
end

subroutine messenger_proc_topology
	use messenger

	return
end


subroutine messenger_syncall()
	use messenger

	return
end

subroutine messenger_free()
	use messenger

	print*, 'End run - Serial version of parallel code'

	return
end

!======================================================================
!			Border Update Subroutines                     =
!======================================================================
subroutine messenger_updateborders(rebuild)
use messenger
implicit none

	integer, intent(in) :: rebuild
	
!todo other cases 
	if (all(periodic.lt.2)) then
	 	call messenger_updateborders_quiescent(rebuild)
	else
		call messenger_updateborders_leesedwards(rebuild)
	end if

end subroutine messenger_updateborders

subroutine messenger_updateborders_quiescent(rebuild)
	use messenger
	implicit none

	integer, intent(in) :: rebuild
	halo_np = 0 !fixedhalo_np

	!Update faces of domain
	if (periodic(1) .eq. 1) then
		call updatefacedown(1)
		call updatefaceup(1)
	endif
	if (periodic(2) .eq. 1) then
		call updatefacedown(2)
		call updatefaceup(2)
	endif
	if (periodic(3) .eq. 1) then
		call updatefacedown(3)
		call updatefaceup(3)
	endif

	!Update edges of domain
	if (periodic(1)+periodic(2) .eq. 2) then
		call updateedge(1,2)
	endif
	if (periodic(2)+periodic(3) .eq. 2) then
		call updateedge(2,3)
	endif
	if (periodic(1)+periodic(3) .eq. 2) then
		call updateedge(1,3)
	endif

	!Check if more array space 
	!should be allocated
	!call allocatecheck()

	!Update Corners of domain
	if (periodic(1)+periodic(2)+periodic(3) .eq. 3) then
		call updatecorners
	endif

	if (rebuild.eq.1) call assign_to_halocell(np+1,np+halo_np)

	return

end

!-------------------------------------------------------------------------------------------
! subroutine: messenger_updateborders_leesedwards
! author: David Trevelyan
! 
! input:	rebuild		- integer, 1 if rebuilding and 0 otherwise
! output:	rebuild		
!
! description:
!
!	messenger_updateborders_leesedwards is the master subroutine that may be used to
!	apply the Lees-Edwards (1972) sliding boundary conditions (note that it must be used
!	in conjunction with sendmols_leesedwards) in any direction. It is possible to correctly
!	update the positions of molecules in the halo cells with only three calls to the
!	subroutine update_plane. The first call updates only the halo cells that border the
!	domain "faces". The second call to update_plane updates the halo cells in a second
!	cartesian direction but also includes the halo cells that were previously established
!	in the first call. Finally, the third call to update_plane includes the halo cells
!	updated in the first two calls - this completes the entire halo correctly.
! 
!-------------------------------------------------------------------------------------------
subroutine messenger_updateborders_leesedwards(rebuild)
use messenger
implicit none

	integer :: rebuild,n

	halo_np = 0
	
	if (iter .ge. shear_iter0) then
		shear_time = dble(iter - shear_iter0)*delta_t
		shear_distance = shear_time*shear_velocity	
	end if
	
	call update_plane(shear_plane,shear_direction,.false.,shear_remainingplane,.false.,rebuild)
	call update_plane(shear_direction,shear_plane,.true.,shear_remainingplane,.false.,rebuild)
	call update_plane(shear_remainingplane,shear_direction,.true.,shear_plane,.true.,rebuild)

	return	

end subroutine messenger_updateborders_leesedwards

!===========================================================================================
!			Periodic Boundaries
!===========================================================================================

!-------------------------------------------------------------------------------------------
! subroutine: update_plane
! author: David Trevelyan
! basis: subroutine updatefacedown(ixyz), author: Ed Smith
!
! inputs:	copyplane	-	plane in which face cells are copied to halo cells
!			loop1plane	-	loop direction 1
!			loop1halos	-	logical argument for inclusion of halo cells in copy routine
!			loop2plane	-	loop direction 2
!			loop2halos	-	logical argument for inclusion of halo cells in copy routine
!			rebuild		-	1 if rebuilding, 0 otherwise
!
! description:
!	
!	update_plane extends the functionality of updatefacedown(ixyz) in order to enable
!	Lees-Edwards sliding boundary conditions. Images of molecules that are assigned to
!	the face cells in the input plane "copyplane" are copied to the corresponding halo 
!	cells. The inputs "loop1plane" and "loop2plane" denote the remaining cartesian 
!	directions, over which looping considers all cells in the "copyplane" surface. 
!	
!	The logical input parameters "loop1halos" and "loop2halos" specify whether the 
!	halo cells in the "loop1plane" and "loop2plane" are included in the image-making
!	routine. If "loop1halos" is .true., the loop in the "loop1plane" direction is taken
!	over cells 1 (bottom halo) to ncells()+2 (top halo), rather than from 2 to ncells()+1.
!
!	At the end of the routine, assign_to_halocell is called so that halo molecules may
!	be correctly copied in subsequent calls to update-plane in the remaining cartesian
!	directions.
!
!-------------------------------------------------------------------------------------------
subroutine update_plane(copyplane,loop1plane,loop1halos,loop2plane,loop2halos,rebuild)
use messenger
use arrays_MD
implicit none

	integer :: p,q,n,m
	integer	:: rebuild
	integer :: copyplane, loop1plane, loop2plane
	integer	:: loop1plane_lower, loop1plane_upper, loop2plane_lower, loop2plane_upper
	integer	:: cellnp, molno, molno_start, molno_finish, startnp 
	integer,dimension(3) :: xyzcell

	logical :: loop1halos, loop2halos
	type(node), pointer	 :: old,current

	startnp = halo_np													!Starting number of halo particles
	m = halo_np
	molno_start = np+m+1

	if (loop1halos) then
		loop1plane_lower = 1											!Set loop limits to include halo cells
		loop1plane_upper = ncells(loop1plane)+2
	else
		loop1plane_lower = 2											!Set loop limits to not include halo cells
		loop1plane_upper = ncells(loop1plane)+1
	end if
	
	if (loop2halos) then
		loop2plane_lower = 1
		loop2plane_upper = ncells(loop2plane)+2
	else
		loop2plane_lower = 2
		loop2plane_upper = ncells(loop2plane)+1	
	end if
	
	!Make images of bottom face cells 
	xyzcell(copyplane) = 2								!Loop over all cells in bottom face of domain
	do p = loop1plane_lower,loop1plane_upper							
	do q = loop2plane_lower,loop2plane_upper
		xyzcell(loop1plane) = p							!Can't iterate an array element
		xyzcell(loop2plane) = q
		cellnp = cell%cellnp(xyzcell(1),xyzcell(2),xyzcell(3))			!Number of particles in the cell
		old => cell%head(xyzcell(1),xyzcell(2),xyzcell(3))%point 		!Set old to top of link list
		do n=1,cellnp
			m = m + 1							!Count one molecule
			molno = old%molno						!Obtain molecule number
			r(np+m,:) = r(molno,:) 		    				!Copy molecule
			v(np+m,:) = v(molno,:)                          !copy velocity
			r(np+m,copyplane) = r(np+m,copyplane) + domain(copyplane)   	!Move to other side of domain
			
			if (potential_flag.eq.1) then					!Polymer IDs copied too
				polyinfo_mol(np+m)%chainID    = polyinfo_mol(molno)%chainID
				polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
			end if
			
			if (copyplane.eq.shear_plane) then
				if (rebuild.eq.1) then				!If rebuilding...
					mol_wrap_integer(molno) = &		!Molecular wrap integer kept the same until next rebuild
					floor((r(np+m,shear_direction)+halfdomain(shear_direction)+shear_distance)/(domain(shear_direction)))
				end if
				r(np+m,shear_direction) = 	r(np+m,shear_direction) + &	!Slide and wrap
									  		(shear_distance-mol_wrap_integer(molno)*domain(shear_direction))
				v(np+m,shear_direction) =   v(np+m,shear_direction) + shear_velocity
			end if

			current => old			!Use current to move to next
			old => current%next 		!Use pointer to obtain next item in list
		enddo
	enddo
	enddo
	
	!Make images of top face cells 
	xyzcell(copyplane) = ncells(copyplane)+1
	do p = loop1plane_lower,loop1plane_upper
	do q = loop2plane_lower,loop2plane_upper
		xyzcell(loop1plane) = p
		xyzcell(loop2plane) = q
		cellnp = cell%cellnp(xyzcell(1),xyzcell(2),xyzcell(3))			!Number of particles in the cell
		old => cell%head(xyzcell(1),xyzcell(2),xyzcell(3))%point 		!Set old to top of link list
		do n=1,cellnp
			m = m + 1													!Count one molecule
			molno = old%molno		    								!Obtain molecule number
			r(np+m,:) = r(molno,:) 		    							!Copy molecule
			v(np+m,:) = v(molno,:)                          !copy velocity
			r(np+m,copyplane) = r(np+m,copyplane) - domain(copyplane)   !Move to other side of domain
			
			if (potential_flag.eq.1) then								!Polymer IDs copied too
				polyinfo_mol(np+m)%chainID    = polyinfo_mol(molno)%chainID
				polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
			end if
		
			if (copyplane.eq.shear_plane) then
				if (rebuild.eq.1) then									!If rebuilding...
					mol_wrap_integer(molno) = &							!Molecular wrap integer kept the same until next rebuild
					-floor((r(np+m,shear_direction)+halfdomain(shear_direction)-shear_distance)/(domain(shear_direction)))
				end if
				r(np+m,shear_direction) = 	r(np+m,shear_direction) - & !Slide and wrap
									  		(shear_distance-mol_wrap_integer(molno)*domain(shear_direction))
				v(np+m,shear_direction) =   v(np+m,shear_direction) - shear_velocity
			end if
		
			current => old			    								!Use current to move to next
			old => current%next 		    							!Use pointer to obtain next item in list
		enddo
	enddo
	enddo

	halo_np = halo_np + m - startnp	
	molno_finish = np+m	

	if (rebuild.eq.1) call assign_to_halocell(molno_start,molno_finish)!Assign all new copies to halo cells
	
	nullify(current)
	nullify(old)

	return

end subroutine update_plane

subroutine updatefacedown(ixyz)
	use messenger
	use arrays_MD
	implicit none

	integer			:: icell, jcell, kcell, ixyz, n, m
	integer			:: cellnp, molno, startnp
	type(node), pointer  	:: old, current

	startnp = halo_np
	m = halo_np + 1 !Set added molecule counter to one more than last added

	select case (ixyz)
	case (1)
		icell = 2
		do jcell =2,ncells(2)+1
		do kcell =2,ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point     !Set old to top of link list
			do n=1,cellnp
				molno = old%molno		    !Obtain molecule number
				r(np+m,:) = r(molno,:) 		    !Copy molecule
				v(np+m,:) = v(molno,:)                          !copy velocity
				r(np+m,1) = r(np+m,1) + domain(1)   !Move to other side of domain

				if (potential_flag.eq.1) then
					polyinfo_mol(np+m)%chainID 	  = polyinfo_mol(molno)%chainID
					polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
				end if

				current => old			    !Use current to move to next
				old => current%next 		    !Use pointer to obtain next item in list
				!print*, r(np+m,1), r(np+m,2), r(np+m,3)
				m = m + 1			    !Update counter of new molecules
			enddo
		enddo
		enddo
	case (2)
		jcell = 2
		do icell =2,ncells(1)+1
		do kcell =2,ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point 	    !Set old to top of link list
			do n=1,cellnp
				molno = old%molno		    !Obtain molecule number
				r(np+m,:) = r(molno,:) 		    !Copy molecule
				v(np+m,:) = v(molno,:)                          !copy velocity
				r(np+m,2) = r(np+m,2) + domain(2)   !Move to other side of domain

				if (potential_flag.eq.1) then
					polyinfo_mol(np+m)%chainID 	  = polyinfo_mol(molno)%chainID
					polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
				end if

				current => old			    !Use current to move to next
				old => current%next 		    !Use pointer to obtain next item in list
				!print*, r(np+m,1), r(np+m,2), r(np+m,3)
				m = m + 1			    !Update counter of new molecules
			enddo
		enddo
		enddo
	case (3)
		kcell = 2
		do icell =2,ncells(1)+1
		do jcell =2,ncells(2)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point !Set old to top of link list
			do n=1,cellnp
				molno = old%molno		    !Obtain molecule number
				r(np+m,:) = r(molno,:) 		    !Copy molecule
				v(np+m,:) = v(molno,:)                          !copy velocity
				r(np+m,3) = r(np+m,3) + domain(3)   !Move to other side of domain

				if (potential_flag.eq.1) then
					polyinfo_mol(np+m)%chainID 	  = polyinfo_mol(molno)%chainID
					polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
				end if

				current => old			    !Use current to move to next
				old => current%next 		    !Use pointer to obtain next item in list
				!print*, r(np+m,1), r(np+m,2), r(np+m,3)
				m = m + 1			    !Update counter of new molecules
			enddo
		enddo
		enddo
	case default
		stop "updateBorder: invalid value for ixyz"
	end select
	
	!Update global number of particles
	halo_np = halo_np + m-1 - startnp

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return

end

!------------------------------------------------------------------------------

subroutine updatefaceup(ixyz)
	use messenger
	use arrays_MD
	implicit none

	integer			:: icell, jcell, kcell, ixyz, n, m 
	integer			:: cellnp, molno, startnp
	type(node), pointer 	:: old, current

	startnp = halo_np
	m = halo_np + 1 !Set added molecule counter to one more than last added

	select case (ixyz)
	case (1)
		icell = ncells(1)+1
		do jcell =2,ncells(2)+1
		do kcell =2,ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point !Set old to top of link list
			do n=1,cellnp
				molno = old%molno		    !Obtain molecule number
				r(np+m,:) = r(molno,:) 		    !Copy molecule
				v(np+m,:) = v(molno,:)                          !copy velocity
				r(np+m,1) = r(np+m,1) - domain(1)   !Move to other side of domain

				if (potential_flag.eq.1) then
					polyinfo_mol(np+m)%chainID 	  = polyinfo_mol(molno)%chainID
					polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
				end if

				current => old			    !Use current to move to next
				old => current%next 		    !Use pointer to obtain next item in list
				!print*, r(np+m,1), r(np+m,2), r(np+m,3)
				m = m + 1			    !Update counter of new molecules
			enddo
		enddo
		enddo
	case (2)
		jcell = ncells(2)+1
		do icell =2,ncells(1)+1
		do kcell =2,ncells(3)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point !Set old to top of link list
			do n=1,cellnp
				molno = old%molno		    !Obtain molecule number
				r(np+m,:) = r(molno,:) 		    !Copy molecule
				v(np+m,:) = v(molno,:)                          !copy velocity
				r(np+m,2) = r(np+m,2) - domain(2)   !Move to other side of domain

				if (potential_flag.eq.1) then
					polyinfo_mol(np+m)%chainID 	  = polyinfo_mol(molno)%chainID
					polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
				end if

				current => old			    !Use current to move to next
				old => current%next 		    !Use pointer to obtain next item in list
				!print*, r(np+m,1), r(np+m,2), r(np+m,3)
				m = m + 1			    !Update counter of new molecules
			enddo
		enddo
		enddo
	case (3)
		kcell = ncells(3)+1
		do icell =2,ncells(1)+1
		do jcell =2,ncells(2)+1
			cellnp = cell%cellnp(icell,jcell,kcell)
			old => cell%head(icell,jcell,kcell)%point !Set old to top of link list
			do n=1,cellnp
				molno = old%molno		    !Obtain molecule number
				r(np+m,:) = r(molno,:) 		    !Copy molecule
				v(np+m,:) = v(molno,:)                          !copy velocity
				r(np+m,3) = r(np+m,3) - domain(3)   !Move to other side of domain

				if (potential_flag.eq.1) then
					polyinfo_mol(np+m)%chainID 	  = polyinfo_mol(molno)%chainID
					polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
				end if

				current => old			    !Use current to move to next
				old => current%next 		    !Use pointer to obtain next item in list
				!print*, r(np+m,1), r(np+m,2), r(np+m,3)
				m = m + 1			    !Update counter of new molecules
			enddo
		enddo
		enddo
	case default
		stop "updateBorder: invalid value for ixyz"
	end select

	!Update global number of particles
	halo_np = halo_np + m-1 - startnp

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return

end


!Halo Edge Cells

subroutine updateedge(face1, face2)
	use messenger
	use arrays_MD
	implicit none

	integer			:: icell, jcell, kcell, ixyz, face1, face2, i, n, m
	integer			:: cellnp, molno, startnp
	integer, dimension(3,4) :: edge1, edge2
	type(node), pointer  	:: old, current

	!Set up all 12 edges
	edge1(1,:) = (/2, 2, ncells(2)+1, ncells(2)+1/)
	edge2(1,:) = (/2, ncells(3)+1, 2, ncells(3)+1/)
	edge1(2,:) = (/2, 2, ncells(1)+1, ncells(1)+1/)
	edge2(2,:) = (/2, ncells(3)+1, 2, ncells(3)+1/)
	edge1(3,:) = (/2, 2, ncells(1)+1, ncells(1)+1/)
	edge2(3,:) = (/2, ncells(2)+1, 2, ncells(2)+1/)

	startnp = halo_np
	m = halo_np + 1 !Set added molecule counter to one more than last added
	ixyz = 6 - face1 - face2 !Determine coordinate along edge

	select case (ixyz)
        case (1)
		do i = 1,4 !Counter for each edge along the x-axis 
		do icell = 2, ncells(1)+1 !Move along x-axis
			cellnp = cell%cellnp(icell,edge1(1,i),edge2(1,i))
			old => cell%head(icell,edge1(1,i),edge2(1,i))%point !Set old to top of link list
			do n=1,cellnp
				molno = old%molno		    !Obtain molecule number
				r(np+m,:) = r(molno,:) 	 	    !Copy molecule
				v(np+m,:) = v(molno,:)                          !copy velocity
				r(np+m,2) = r(np+m,2) &  !Move to other side of domain
				+ sign(1,ncells(2)-edge1(1,i))*domain(2)
				r(np+m,3) = r(np+m,3) &  !Move to other side of domain
				+ sign(1,ncells(3)-edge2(1,i))*domain(3)

				if (potential_flag.eq.1) then
					polyinfo_mol(np+m)%chainID 	  = polyinfo_mol(molno)%chainID
					polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
				end if

				current => old			    !Use current to move to next
				old => current%next 		    !Use pointer to obtain next item in list
				!print*, r(np+m,1), r(np+m,2), r(np+m,3)
				m = m + 1                !Update counter of new molecules
			enddo
		enddo
		enddo
        case (2)
		do i = 1,4 !Counter for each edge along the x-axis 
		do jcell = 2, ncells(2)+1 !Move along x-axis
			cellnp = cell%cellnp(edge1(2,i),jcell,edge2(2,i))
			old => cell%head(edge1(2,i),jcell,edge2(2,i))%point !Set old to top of link list
			do n=1,cellnp
				molno = old%molno		    !Obtain molecule number
				r(np+m,:) = r(molno,:) 	 	    !Copy molecule
				v(np+m,:) = v(molno,:)                          !copy velocity
				r(np+m,1) = r(np+m,1) &  !Move to other side of domain
				+ sign(1,ncells(1)-edge1(2,i))*domain(1)
				r(np+m,3) = r(np+m,3) &  !Move to other side of domain
				+ sign(1,ncells(3)-edge2(2,i))*domain(3)

				if (potential_flag.eq.1) then
					polyinfo_mol(np+m)%chainID 	  = polyinfo_mol(molno)%chainID
					polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
				end if

				current => old			    !Use current to move to next
				old => current%next 		    !Use pointer to obtain next item in list
				!print*, r(np+m,1), r(np+m,2), r(np+m,3)
				m = m + 1                !Update counter of new molecules
			enddo
		enddo
		enddo

        case (3)
		do i = 1,4 !Counter for each edge along the x-axis 
		do kcell = 2, ncells(3)+1 !Move along x-axis
			cellnp = cell%cellnp(edge1(3,i),edge2(3,i),kcell)
			old => cell%head(edge1(3,i),edge2(3,i),kcell)%point !Set old to top of link list
			do n=1,cellnp
				molno = old%molno		    !Obtain molecule number
				r(np+m,:) = r(molno,:) 	 	    !Copy molecule
				v(np+m,:) = v(molno,:)                          !copy velocity
				r(np+m,1) = r(np+m,1) &  !Move to other side of domain
				+ sign(1,ncells(1)-edge1(3,i))*domain(1)
				r(np+m,2) = r(np+m,2) &  !Move to other side of domain
				+ sign(1,ncells(2)-edge2(3,i))*domain(2)

				if (potential_flag.eq.1) then
					polyinfo_mol(np+m)%chainID 	  = polyinfo_mol(molno)%chainID
					polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
				end if

				current => old			    !Use current to move to next
				old => current%next 		    !Use pointer to obtain next item in list
				!print*, r(np+m,1), r(np+m,2), r(np+m,3)
				m = m + 1                !Update counter of new molecules
			enddo
		enddo
		enddo
	case default
		stop "updateBorder: invalid value for ixyz"
	end select

	!Update global number of particles
	halo_np = halo_np + m-1 - startnp

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return

end


!Halo Corner Cells

subroutine updatecorners()
	use messenger
	use arrays_MD
	implicit none

	integer			:: i, n, m
	integer			:: cellnp, molno, startnp
	integer, dimension(8)   :: icornercell, jcornercell, kcornercell
	type(node), pointer  	:: old, current

	!Set up all 8 corners
 	icornercell = (/ 2, 2, 2, 2, ncells(1)+1, ncells(1)+1, ncells(1)+1, ncells(1)+1/)
 	jcornercell = (/ 2, 2, ncells(2)+1, ncells(2)+1, 2, 2, ncells(2)+1, ncells(2)+1/) 
 	kcornercell = (/ 2, ncells(3)+1, 2, ncells(3)+1, 2, ncells(3)+1, 2, ncells(3)+1/)

	startnp = halo_np
	m = halo_np + 1 !Set added molecule counter to one more than last added

	do i=1,8 !Counter for each corner cell
		cellnp = cell%cellnp(icornercell(i),jcornercell(i),kcornercell(i))
		old => cell%head(icornercell(i),jcornercell(i),kcornercell(i))%point !Set old to top of link list
		do n=1,cellnp
			molno = old%molno	 !Obtain molecule number
			r(np+m,:) = r(molno,:) 	 !Copy molecule
			v(np+m,:) = v(molno,:)                          !copy velocity
			r(np+m,1) = r(np+m,1) &  !Move to other side of domain
			+ sign(1,ncells(1)-icornercell(i))*domain(1)
			r(np+m,2) = r(np+m,2) &  !Move to other side of domain
			+ sign(1,ncells(2)-jcornercell(i))*domain(2)
			r(np+m,3) = r(np+m,3) &  !Move to other side of domain
			+ sign(1,ncells(3)-kcornercell(i))*domain(3)

			if (potential_flag.eq.1) then
				polyinfo_mol(np+m)%chainID 	  = polyinfo_mol(molno)%chainID
				polyinfo_mol(np+m)%subchainID = polyinfo_mol(molno)%subchainID
			end if

			current => old			    !Use current to move to next
			old => current%next 		    !Use pointer to obtain next item in list
			!print*, r(np+m,1), r(np+m,2), r(np+m,3)
			m = m + 1                !Update counter of new molecules
		enddo
	enddo

	!Update global number of particles
	halo_np = halo_np + m-1 - startnp
	!print*, 'halonp =', halo_np

	nullify(current)        !Nullify current as no longer required
	nullify(old)            !Nullify old as no longer required

	return

end

!---------------------------------------------------------------------------------
!Check if array allocation is big enough for halo molecules

subroutine allocatecheck()
	use messenger
	use arrays_MD
	implicit none

	double precision, dimension(:,:), allocatable :: tempr

	!Check allocated space is less than half filled at half way
	if (nint((size(r,1)-np)/2.) < halo_np) then
		print*, 'Increasing halo array allocation space'
		allocate(tempr(size(r,1),nd))
		tempr = r
		deallocate(r)
		allocate(r((np+halo_np+ceiling(np*0.001)+20),nd)) !Increase in size
		r = tempr
		deallocate(tempr)
	endif

end

!======================================================================
!			Molecule Transfer Subroutines                 =
!======================================================================
subroutine sendmols()
use messenger
implicit none
	
	if (all(periodic.lt.2)) then
		call sendmols_quiescent
	else
		call sendmols_leesedwards
	end if

end subroutine sendmols

subroutine sendmols_quiescent()
	use messenger
	use arrays_MD
	implicit none

	integer :: ixyz, n

	do n=1,np

		!Domain goes from -halfdomain to +halfdomain
		if (r(n,1) >= halfdomain(1)) then   !Above +halfdomain
			r(n,1) = r(n,1) - domain(1) !Move to other side of domain
		end if
		if (r(n,1) < -halfdomain(1)) then   !Below -halfdomain
			r(n,1) = r(n,1) + domain(1) !Move to other side of domain
		endif

		!Domain goes from -halfdomain to +halfdomain
		if (r(n,2) >= halfdomain(2)) then   !Above +halfdomain
			r(n,2) = r(n,2) - domain(2) !Move to other side of domain
		end if
		if (r(n,2) < -halfdomain(2)) then   !Below -halfdomain
			r(n,2) = r(n,2) + domain(2) !Move to other side of domain
		endif

		!Domain goes from -halfdomain to +halfdomain
		if (r(n,3) >= halfdomain(3)) then   !Above +halfdomain
			r(n,3) = r(n,3) - domain(3) !Move to other side of domain
		end if
		if (r(n,3) < -halfdomain(3)) then   !Below -halfdomain
			r(n,3) = r(n,3) + domain(3) !Move to other side of domain
		endif

	enddo

	return
end

subroutine sendmols_leesedwards()
use messenger
use arrays_MD
implicit none
	
	integer :: ixyz,n

	if (iter.lt.shear_iter0) then
		call sendmols_quiescent
		return
	end if
	
	shear_time = dble(iter - shear_iter0)*delta_t
	shear_distance = shear_time*shear_velocity
	wrap_integer = floor(shear_time*shear_velocity/domain(shear_direction))
	
	do n=1,np

		!---- Slide and wrap in shearing plane first --------------------!
		if (r(n,shear_plane) .ge. halfdomain(shear_plane)) then   									!Above +halfdomain
			r(n,shear_plane) = r(n,shear_plane) - domain(shear_plane) 								!Move to other side of domain
			r(n,shear_direction) = r(n,shear_direction) - (shear_distance - wrap_integer*domain(shear_direction))
			v(n,shear_direction) = v(n,shear_direction) - shear_velocity
		end if
		if (r(n,shear_plane) .lt. -halfdomain(shear_plane)) then   									!Below -halfdomain
			r(n,shear_plane) = r(n,shear_plane) + domain(shear_plane) 								!Move to other side of domain
			r(n,shear_direction) = r(n,shear_direction) + (shear_distance - wrap_integer*domain(shear_direction))
			v(n,shear_direction) = v(n,shear_direction) + shear_velocity
		endif
		!----------------------------------------------------------------!

		if (r(n,shear_direction) >= halfdomain(shear_direction)) then   							!Above +halfdomain
			r(n,shear_direction) = r(n,shear_direction) - domain(shear_direction) 					!Move to other side of domain
		end if			
		if (r(n,shear_direction) < -halfdomain(shear_direction)) then   							!Below -halfdomain
			r(n,shear_direction) = r(n,shear_direction) + domain(shear_direction) 					!Move to other side of domain
		endif

		if (r(n,shear_remainingplane) >= halfdomain(shear_remainingplane)) then   					!Above +halfdomain
			r(n,shear_remainingplane) = r(n,shear_remainingplane) - domain(shear_remainingplane) 	!Move to other side of domain
		end if
		if (r(n,shear_remainingplane) < -halfdomain(shear_remainingplane)) then   					!Below -halfdomain
			r(n,shear_remainingplane) = r(n,shear_remainingplane) + domain(shear_remainingplane) 	!Move to other side of domain
		endif
		
	enddo
	
	return

end subroutine sendmols_leesedwards
!======================================================================
!			Data Transfer Subroutines                     =
!======================================================================

subroutine send_VA_interaction(Rfbin_halo,intercbin)
	use messenger
	use calculated_properties_MD
	implicit none

	integer				:: ibin, jbin, kbin
	integer,dimension(3)		:: intercbin
	double precision,dimension(3,3)	:: Rfbin_halo

	ibin = modulo((intercbin(1)-1),nbins(1))+1
	jbin = modulo((intercbin(2)-1),nbins(1))+1
	kbin = modulo((intercbin(3)-1),nbins(1))+1

	!Map from halo to bin it represents
	rfbin(ibin,jbin,kbin,:,:) = rfbin(ibin,jbin,kbin,:,:) + Rfbin_halo(:,:)

	return
end
!======================================================================
!			Data gathering subroutines                    =
!======================================================================

subroutine globalbroadcast(A,na,broadprocid)
	use messenger
	implicit none

	integer			:: na, broadprocid
	double precision	:: A

	A = A

	return
end

subroutine globalsyncreduce(A, na, meanA, maxA, minA)
	use messenger
	implicit none
	
        integer na
	double precision, intent(in)  :: A(na)
	double precision, intent(out) :: meanA(na), maxA(na), minA(na)

	meanA = A
	maxA = A
	minA = A

	return
end

subroutine globalSum(A)
	use messenger
	
	double precision A

	A = A

	return
end

subroutine globalSumInt(A)
	use messenger
	
	integer A

	A = A

	return
end


subroutine globalMaxInt(A)
	use messenger
	
	integer A

	A = A

	return
end

subroutine globalSumVect(A, na)
	use messenger
	
        integer na
	double precision A(na)

	A = A

	return
end

subroutine globalSumIntVect(A, na)
	use messenger
	
        integer na
	integer A(na)

	A = A

	return
end


subroutine globalMaxIntVect(A, na)
	use messenger
	
        integer na
	integer A(na)

	A = A

	return
end

subroutine globalMaxVect(A, na)
	use messenger
	
        integer na
	double precision A(na)

	A = A

	return
end

subroutine globalMinVect(A, na)
	use messenger
	
        integer na
	double precision A(na)

	A = A

	return
end

subroutine globalAverage(A, na)
	use messenger
	
        integer na
	double precision A(na)

	A = A

	return
end

subroutine globalGathernp()
	use physical_constants_MD
	use messenger
	
	globalnp = np

	return
end

!----Sum routines over global sub communitcators

subroutine SubcommSumInt(A, ixyz)
	use messenger
	!include "mpif.h"

        integer, intent(in) :: ixyz !Direction of sub-comm
	integer	A
	integer buf

	A = A

	return
end

subroutine SubcommSumVect(A, na, ixyz)
	use messenger
	!include "mpif.h"

        integer, intent(in) :: na, ixyz !Direction of sub-comm
	double precision A(na)
	double precision buf(na)

	A = A

	return
end

subroutine SubcommSumIntVect(A, na, ixyz)
	use messenger
	!include "mpif.h"

        integer, intent(in) :: na, ixyz !Direction of sub-comm
	integer	A(na)
	integer buf(na)

	A = A 

	return
end

subroutine error_abort(msg)
        use mpi
        implicit none
        
        character(len=*), intent(in), optional :: msg

        if(present(msg)) then
                write(*,*) msg
        endif

        stop

end subroutine error_abort
