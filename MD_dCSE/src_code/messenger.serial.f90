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
! Molecule Sending Subroutines =
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

end module

!======================================================================
!			Key Messenger Subroutines                     =
!======================================================================
subroutine messenger_invoke()
	use messenger

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

subroutine messenger_updateborders()
	use messenger
	implicit none

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

	return

end

!-------------------------------------------------------------------
!			Periodic Boundries
!-------------------------------------------------------------------

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
				r(np+m,1) = r(np+m,1) + domain(1)   !Move to other side of domain
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
				r(np+m,2) = r(np+m,2) + domain(2)   !Move to other side of domain
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
				r(np+m,3) = r(np+m,3) + domain(3)   !Move to other side of domain
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
				r(np+m,1) = r(np+m,1) - domain(1)   !Move to other side of domain
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
				r(np+m,2) = r(np+m,2) - domain(2)   !Move to other side of domain
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
				r(np+m,3) = r(np+m,3) - domain(3)   !Move to other side of domain
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
				r(np+m,2) = r(np+m,2) &  !Move to other side of domain
				+ sign(1,ncells(2)-edge1(1,i))*domain(2)
				r(np+m,3) = r(np+m,3) &  !Move to other side of domain
				+ sign(1,ncells(3)-edge2(1,i))*domain(3)
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
				r(np+m,1) = r(np+m,1) &  !Move to other side of domain
				+ sign(1,ncells(1)-edge1(2,i))*domain(1)
				r(np+m,3) = r(np+m,3) &  !Move to other side of domain
				+ sign(1,ncells(3)-edge2(2,i))*domain(3)
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
				r(np+m,1) = r(np+m,1) &  !Move to other side of domain
				+ sign(1,ncells(1)-edge1(3,i))*domain(1)
				r(np+m,2) = r(np+m,2) &  !Move to other side of domain
				+ sign(1,ncells(2)-edge2(3,i))*domain(2)
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
			r(np+m,1) = r(np+m,1) &  !Move to other side of domain
			+ sign(1,ncells(1)-icornercell(i))*domain(1)
			r(np+m,2) = r(np+m,2) &  !Move to other side of domain
			+ sign(1,ncells(2)-jcornercell(i))*domain(2)
			r(np+m,3) = r(np+m,3) &  !Move to other side of domain
			+ sign(1,ncells(3)-kcornercell(i))*domain(3)
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

