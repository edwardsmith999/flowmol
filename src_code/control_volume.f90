
!===================================================================================
!
!	Routines add forces and fluxes per molecule to running totals
!
!===================================================================================

module control_volume
	implicit none

contains

subroutine check_CV_conservation(mflux_outflag,vflux_outflag,eflux_outflag)
	use calculated_properties_MD, only : dmdt,mass_flux,momentum_flux,Pxyface,dmvdt,nbins, & 
										 volume_mass_pdt, volume_momentum_pdt,binsize
	use computational_constants_MD, only : iter, irank, tplot,delta_t,Nvflux_ave
	use CV_objects, only : CVcheck_mass, CVcheck_momentum
	implicit none

	integer,intent(in)	:: mflux_outflag,vflux_outflag,eflux_outflag

	integer				:: i,j,k

	call CVcheck_mass%check_error(2,nbins(1)+1,2,nbins(2)+1,2,nbins(3)+1,iter,irank)
	call CVcheck_momentum%check_error(2,nbins(1)+1,2,nbins(2)+1,2,nbins(3)+1,iter,irank)

	!call CVcheck_mass%check_error(8,8,11,11,8,8,iter,irank)


	!CV snapshots taken every tplot 
	!if (mod(iter,tplot) .ne. 0) return



!	do i = 2,nbins(1)+1
!	do j = 2,nbins(2)+1
!	do k = 2,nbins(3)+1

		!if (mflux_outflag .eq. 1) then
			!if(sum(mass_flux(i,j,k,:))-dmdt(i,j,k) .ne. 0) then
		!		 print'(a,i8,4i4,3i8)', & 
		!			'Error in mass flux', iter,irank,i,j,k,sum(mass_flux(i,j,k,:)),volume_mass_pdt(i,j,k),dmdt(i,j,k)
			!endif
		!endif

		!if (vflux_outflag .eq. 4) then
		!	if (abs(sum(sum(momentum_flux(i,j,k,:,:))+sum(Pxyface(i,j,k,:,:))-dmvdt(i,j,k,:))) .gt. 0.0001d0) then
		!		print'(a,5i5,9f10.5)','Error in momentum flux ', &
		!			    iter,i,j,k,irank,sum(momentum_flux(i,j,k,1,:)),sum(Pxyface(i,j,k,1,:)),volume_momentum_pdt(i,j,k,1) & 
		!				       	   		,sum(momentum_flux(i,j,k,2,:)),sum(Pxyface(i,j,k,2,:)),volume_momentum_pdt(i,j,k,2) & 
		!				     			,sum(momentum_flux(i,j,k,3,:)),sum(Pxyface(i,j,k,3,:)),volume_momentum_pdt(i,j,k,3)
		!	endif
		!endif

		!if (eflux_outflag .eq. 4) then
		!	stop "CV Energy debug is not developed"
		!endif

!	enddo
!	enddo
!	enddo

end subroutine check_CV_conservation


!===================================================================================
! Momentum Flux over a surface of a bin including all intermediate bins

subroutine get_molecule_CV_momentum_flux(molno)
	use computational_constants_MD, only : domain, delta_t, halfdomain, nhb
	use calculated_properties_MD, only : nbins, momentum_flux
	use arrays_MD, only : r, v
	use librarymod, only : heaviside, imaxloc
	implicit none

	integer,intent(in)				:: molno

	integer							:: jxyz,i,j,k
	integer							:: planeno
	integer							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	integer		,dimension(3)		:: ibin1,ibin2,cbin
	double precision				:: crossplane,rplane,shift
	double precision,dimension(3)	:: mbinsize,velvect,crossface
	double precision,dimension(3)	:: ri1,ri2,ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb

	!CV momentum flux
	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	!Get velocity at v(t+dt/2) from v(t-dt/2)
	velvect(:) = v(:,molno)
	ri1(:) = r(:,molno) 						!Molecule i at time t
	ri2(:) = r(:,molno)	- delta_t*velvect(:)	!Molecule i at time t-dt
	ri12   = ri1 - ri2							!Molecule i trajectory between t-dt and t
	where (ri12 .eq. 0.d0) ri12 = 0.000001d0

	!Assign to bins before and after using integer division
	ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb(:)
	ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb(:)

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
			!crosstime = (r(jxyz,n) - rplane)/v(jxyz,n)
			!velvect(:) = v(:,n) !- a(:,n) * crosstime
			!Change in velocity at time of crossing is not needed as velocity assumed constant 
			!for timestep and changes when forces are applied.

			!Add Momentum flux over face
			momentum_flux(cbin(1),cbin(2),cbin(3),:,1) = & 
				momentum_flux(cbin(1),cbin(2),cbin(3),:,1) & 
			      - velvect(:)*dble(onfacexb)*abs(crossface(jxyz))
			momentum_flux(cbin(1),cbin(2),cbin(3),:,2) = & 
				momentum_flux(cbin(1),cbin(2),cbin(3),:,2) & 
			      - velvect(:)*dble(onfaceyb)*abs(crossface(jxyz))
			momentum_flux(cbin(1),cbin(2),cbin(3),:,3) = & 
				momentum_flux(cbin(1),cbin(2),cbin(3),:,3) &
			      - velvect(:)*dble(onfacezb)*abs(crossface(jxyz))
			momentum_flux(cbin(1),cbin(2),cbin(3),:,4) = & 
				momentum_flux(cbin(1),cbin(2),cbin(3),:,4) &
			      + velvect(:)*dble(onfacext)*abs(crossface(jxyz))
			momentum_flux(cbin(1),cbin(2),cbin(3),:,5) = & 
				momentum_flux(cbin(1),cbin(2),cbin(3),:,5) &
			      + velvect(:)*dble(onfaceyt)*abs(crossface(jxyz))
			momentum_flux(cbin(1),cbin(2),cbin(3),:,6) = & 
				momentum_flux(cbin(1),cbin(2),cbin(3),:,6) &
			      + velvect(:)*dble(onfacezt)*abs(crossface(jxyz))

		enddo
		enddo
		enddo

		!if ((ibin1(1) .eq. 4 .and. ibin1(2) .eq. 4 .and. ibin1(3) .eq. 4) .or. &
		!	(ibin2(1) .eq. 4 .and. ibin2(2) .eq. 4 .and. ibin2(3) .eq. 4)) then
		!	print'(a,8i4,18f7.3)', 'mol',iter,ibin1,ibin2, molno,momentum_flux(4,4,4,:,:)
	!endif
			

	endif

end subroutine get_molecule_CV_momentum_flux

!===================================================================================
!Forces over the surface of a Volume

subroutine get_molecule_CV_forces(fij,ri,rj,molnoi,molnoj)
	use computational_constants_MD, only : domain, delta_t, halfdomain, nhb
	use calculated_properties_MD, only : nbins, momentum_flux, volume_force
	use physical_constants_MD, only: np,nd
	!use arrays_MD, only : r, v
	use librarymod, only : heaviside
	implicit none

	integer,intent(in)							:: molnoi,molnoj
	double precision,dimension(3),intent(in)	:: ri,rj,fij

	integer,dimension(3)			:: ibin, jbin
	double precision,dimension(3)	:: crossplane,fsurface
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

end subroutine get_molecule_CV_forces

!===================================================================================
!Forces over the surface of a Volume

subroutine get_molecule_CV_stresses(fij,ri,rj,molnoi)
	use computational_constants_MD, only : domain, delta_t, halfdomain, nhb, eflux_outflag
	use calculated_properties_MD, only : nbins, momentum_flux, volume_force, Pxyface, Pxyvface
	use physical_constants_MD, only: np,nd
	use arrays_MD, only : r, v
	use librarymod, only : heaviside
	implicit none


	integer,intent(in)							:: molnoi
	double precision,dimension(3),intent(in)	:: ri,rj,fij

	integer							:: i,j,k,ixyz
	integer							:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
	integer,dimension(3)			:: cbin, ibin, jbin
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

		!Stress acting on face over volume
		if (eflux_outflag .ne. 0) then
			velvect(:) = v(:,molnoi) 
			!velvect(:) = v(:,molnoi) + 0.5d0*delta_t*a(:,molnoi)
			Pxyvface(cbin(1),cbin(2),cbin(3),1) = Pxyvface(cbin(1),cbin(2),cbin(3),1) + dot_product(fij,velvect)*dble(onfacexb)
			Pxyvface(cbin(1),cbin(2),cbin(3),2) = Pxyvface(cbin(1),cbin(2),cbin(3),2) + dot_product(fij,velvect)*dble(onfaceyb)
			Pxyvface(cbin(1),cbin(2),cbin(3),3) = Pxyvface(cbin(1),cbin(2),cbin(3),3) + dot_product(fij,velvect)*dble(onfacezb)
			Pxyvface(cbin(1),cbin(2),cbin(3),4) = Pxyvface(cbin(1),cbin(2),cbin(3),4) + dot_product(fij,velvect)*dble(onfacext)
			Pxyvface(cbin(1),cbin(2),cbin(3),5) = Pxyvface(cbin(1),cbin(2),cbin(3),5) + dot_product(fij,velvect)*dble(onfaceyt)
			Pxyvface(cbin(1),cbin(2),cbin(3),6) = Pxyvface(cbin(1),cbin(2),cbin(3),6) + dot_product(fij,velvect)*dble(onfacezt)
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

end subroutine get_molecule_CV_stresses


!===================================================================================
!
!	Routines to obtain forces and fluxes for entire CV
!
!===================================================================================

!-----------------------------------------------------------------------------------
! Flux of molecules over bin surface requires all molecules in current cell and 
! surrounding cells to be checked

subroutine get_CV_momentum_flux(icell,jcell,kcell,isumflux,Flux)
	use module_linklist, only : node, cell
	use arrays_MD, only : r,v,a
	use computational_constants_MD, only : delta_t, domain,halfdomain, nhb,iter
	use calculated_properties_MD, only	 : nbins
	use librarymod, only : imaxloc
	implicit none

	integer,intent(in)										:: icell,jcell,kcell
	double precision,intent(inout)							:: isumflux
	double precision,optional,dimension(3,6),intent(out)	:: Flux

	integer								:: n,molno, jxyz
	integer								:: icellshift,jcellshift,kcellshift,adjacentcellnp
	integer,dimension(3)				:: ibin1, ibin2, crossface
	double precision,dimension(3)		:: mbinsize,ri1,ri2,ri12,Fsurface,velvect
	type(node), pointer		 	        :: old, current


	!Determine bin size
	mbinsize(:) = domain(:) / nbins(:)

	!Calculate bin surface fluxes
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

			velvect(:) = v(:,molno)
			ri1(:) = r(:,molno) 							!Molecule i at time t
			ri2(:) = r(:,molno)	- delta_t*velvect			!Molecule i at time t-dt

			Fsurface = 0.d0
			! *********************************************************************************
			!Calculate flux over surface only if molecule is entering/leaving bin of interest
			if (present(Flux)) then
				call get_Flux(icell,jcell,kcell,velvect,ri1,ri2,Fsurface,Flux)	
			else
				call get_Flux(icell,jcell,kcell,velvect,ri1,ri2,Fsurface)	
			endif

			!print'(a,i4,6i2,i4,18f6.3)', '!CV',iter,icell,jcell,kcell,icell+icellshift,jcell+jcellshift,kcell+kcellshift,molno,Flux

			isumflux = isumflux + Fsurface(1)

			! *********************************************************************************
			current => old
			old => current%next    !Use pointer in datatype to obtain next item in list

		enddo
	enddo
	enddo
	enddo

contains

	!-----------------------------------------------------------------------------------
	!Flux over the surface of a bin

	subroutine get_Flux(icell,jcell,kcell,velvect,ri1,ri2,Fsurface,flux)
		use computational_constants_MD, only : domain,halfdomain, iter
		use calculated_properties_MD, only	 : nbins
		use librarymod, only : heaviside, imaxloc
		implicit none

		integer,intent(in)										:: icell,jcell,kcell
		double precision,dimension(3),intent(in)				:: ri1, ri2, velvect
		double precision,dimension(3),intent(out)				:: Fsurface
		double precision,dimension(3,6),intent(inout),optional	:: Flux

		integer									:: ixyz,jxyz,kxyz,i,j,k,ii,jj,kk,n
		integer									:: planeno
		integer									:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
		integer		,dimension(3)				:: ibin1,ibin2,cbin
		double precision						:: crosstime,crossplane,rplane,shift
		double precision,dimension(3)			:: mbinsize,crossface
		double precision,dimension(3)			:: ri12,bintop,binbot,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
		double precision,dimension(3,3,3,3)		:: Fsurfacebins
		double precision,dimension(3,3,3,3,6)	:: Fluxbins

		!CV momentum flux
		!Determine bin size
		mbinsize(:) = domain(:) / nbins(:)

		ri12   = ri1 - ri2		!Molecule i trajectory between t-dt and t
		where (ri12 .eq. 0.d0) ri12 = 0.000001d0

		!Assign to bins before and after using integer division
		ibin1(:) = ceiling((ri1+halfdomain(:))/mbinsize(:)) + nhb(:)
		ibin2(:) = ceiling((ri2+halfdomain(:))/mbinsize(:)) + nhb(:)

		!Replace Signum function with this functions which gives a
		!check for plane crossing and the correct sign 
		crossface(:) =  ibin1(:) - ibin2(:)

		if (sum(abs(crossface(:))) .ne. 0) then

			Fluxbins = 0.d0

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

				fsurface(:) = 0.d0
				fsurface(:) = fsurface(:) - 0.5d0*velvect(:)*dble(onfacexb - onfacext)
				fsurface(:) = fsurface(:) - 0.5d0*velvect(:)*dble(onfaceyb - onfaceyt)
				fsurface(:) = fsurface(:) - 0.5d0*velvect(:)*dble(onfacezb - onfacezt)

				jxyz = imaxloc(abs(crossface))	!Integer array of size 1 copied to integer

				!Add Momentum flux over face
				ii = modulo(i,3)+1; jj = modulo(j,3)+1; kk = modulo(k,3)+1
				Fluxbins(ii,jj,kk,:,1) = Fluxbins(ii,jj,kk,:,1) - velvect(:)*dble(onfacexb)*abs(crossface(jxyz))
				Fluxbins(ii,jj,kk,:,2) = Fluxbins(ii,jj,kk,:,2) - velvect(:)*dble(onfaceyb)*abs(crossface(jxyz))
				Fluxbins(ii,jj,kk,:,3) = Fluxbins(ii,jj,kk,:,3) - velvect(:)*dble(onfacezb)*abs(crossface(jxyz))
				Fluxbins(ii,jj,kk,:,4) = Fluxbins(ii,jj,kk,:,4) + velvect(:)*dble(onfacext)*abs(crossface(jxyz))
				Fluxbins(ii,jj,kk,:,5) = Fluxbins(ii,jj,kk,:,5) + velvect(:)*dble(onfaceyt)*abs(crossface(jxyz))
				Fluxbins(ii,jj,kk,:,6) = Fluxbins(ii,jj,kk,:,6) + velvect(:)*dble(onfacezt)*abs(crossface(jxyz))


			enddo
			enddo
			enddo
			
			!Add surface force to current bin
			!Fsurface(:) = sum(Fluxbins(modulo(icell,3)+1, 	& 
			!	     		   modulo(jcell,3)+1, 	& 
			!	    		   modulo(kcell,3)+1,:,:),2)

			!Take flux from central bin only
			if (present(Flux))then
				Flux(:,:) = Flux(:,:) +  Fluxbins(modulo(icell,3)+1, 	& 
				     		     	 			  modulo(jcell,3)+1, 	& 
				    		     	 	          modulo(kcell,3)+1,:,:)
			endif


			!if ((ibin1(1) .eq. 4 .and. ibin1(2) .eq. 4 .and. ibin1(3) .eq. 4) .or. &
			!	(ibin2(1) .eq. 4 .and. ibin2(2) .eq. 4 .and. ibin2(3) .eq. 4)) then
			!	print'(a,i4,6i2,i4,18f6.3)', 'CV ',iter,ibin1,ibin2,molno,Flux!Fluxbins(modulo(icell,3)+1,modulo(jcell,3)+1,modulo(kcell,3)+1,:,:)!
			!endif
					
		endif

	end subroutine get_Flux

end subroutine get_CV_momentum_flux

!-----------------------------------------------------------------------------------
! Forces between current bin and surrounding bins

subroutine get_CV_force(icell,jcell,kcell,isumforce,Traction)
	use module_linklist, only : node, cell
	use arrays_MD, only : r
	use computational_constants_MD, only : domain,halfdomain
	use calculated_properties_MD, only	 : nbins, Pxyface
	use physical_constants_MD, only 	 : rcutoff2,np,nd
	implicit none

	integer,intent(in)										:: icell,jcell,kcell
	double precision,intent(inout)							:: isumforce
	double precision,dimension(3,6),optional,intent(inout)	:: Traction

	integer							:: n,j,ixyz,molnoi,molnoj
	integer							:: icellshift,jcellshift,kcellshift,cellnp,adjacentcellnp
	integer,dimension(3)			:: ibin, jbin
	double precision				:: rij2, invrij2, accijmag
	double precision,dimension(3)	:: ri,rj,rij,fij,Fsurface
	type(node), pointer		 	 	:: oldi, currenti, oldj, currentj

	cellnp = cell%cellnp(icell,jcell,kcell)
	oldi => cell%head(icell,jcell,kcell)%point 	!Set old to first molecule in list

	!Calculate averages for bin
	do n = 1, cellnp    ! Loop over all particles

		molnoi = oldi%molno	!Number of molecule
		ri = r(:,molnoi)	!Retrieve ri

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
					if (present(Traction)) then
						call get_Traction(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface,Traction)
					else
						call get_Fsurface(molnoi,molnoj,fij,ri,rj,Fsurface)
					endif
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


	!-----------------------------------------------------------------------------------
	!Forces over all of the surface of a bin

	subroutine get_Fsurface(molnoi,molnoj,fij,ri,rj,Fsurface)
		use librarymod, only : heaviside
		implicit none

		integer,intent(in)							:: molnoi,molnoj
		double precision,dimension(3),intent(in)	:: ri, rj, fij
		double precision,dimension(3),intent(out)	:: Fsurface

		integer							:: ixyz
		integer,dimension(3)			:: ibin, jbin
		double precision,dimension(3)	:: Fbinsize, bintop, binbot, crossplane

		!Determine bin size
		Fbinsize(:) = domain(:) / nbins(:)

		!Assign to bins using integer division
		ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+1  !Establish current bin

		bintop(:) = (ibin(:)-1)*Fbinsize(:)-halfdomain(:)
		binbot(:) = (ibin(:)-2)*Fbinsize(:)-halfdomain(:)

		!Add for molecule i
		if(molnoi .le. np) then
			Fsurface = fij(:)* dble((heaviside(bintop(1)-ri(1))-heaviside(binbot(1)-ri(1)))* & 
				  					(heaviside(bintop(2)-ri(2))-heaviside(binbot(2)-ri(2)))* & 
				  					(heaviside(bintop(3)-ri(3))-heaviside(binbot(3)-ri(3)))- & 
				  					(heaviside(bintop(1)-rj(1))-heaviside(binbot(1)-rj(1)))* & 
				  					(heaviside(bintop(2)-rj(2))-heaviside(binbot(2)-rj(2)))* & 
				  					(heaviside(bintop(3)-rj(3))-heaviside(binbot(3)-rj(3))))
		endif

	end subroutine get_Fsurface

	!-----------------------------------------------------------------------------------
	! Tractions on one surface of a bin

	subroutine get_Traction(icell,jcell,kcell,molnoi,molnoj,fij,ri,rj,Fsurface,Traction)
		use librarymod, only : heaviside
		implicit none

		integer,intent(in)										:: molnoi,molnoj
		double precision,dimension(3),intent(in)				:: ri,rj,fij
		double precision,dimension(3),intent(out)				:: Fsurface
		double precision,dimension(3,6),intent(inout),optional	:: Traction

		integer									:: i,j,k,ixyz,n,tempi
		integer									:: icell,jcell,kcell
		integer									:: onfacext,onfacexb,onfaceyt,onfaceyb,onfacezt,onfacezb
		integer,dimension(3)					:: cbin, ibin, jbin
		double precision						:: binforce
		double precision,dimension(3)			:: rij,Pxt,Pxb,Pyt,Pyb,Pzt,Pzb
		double precision,dimension(3)			:: Fbinsize, bintop, binbot
		double precision,dimension(3,3,3,3,6)	:: Tractionbins
		!Calculate rij
		rij = ri - rj
		!Prevent Division by zero
		do ixyz = 1,3
			if (abs(rij(ixyz)) .lt. 0.000001d0) rij(ixyz) = sign(0.000001d0,rij(ixyz))
		enddo

		!Determine bin size
		Fbinsize(:) = domain(:) / nbins(:)

		!Assign to bins using integer division
		ibin(:) = ceiling((ri(:)+halfdomain(:))/Fbinsize(:))+1	!Establish current bin
		jbin(:) = ceiling((rj(:)+halfdomain(:))/Fbinsize(:))+1 	!Establish current bin

		do i = ibin(1),jbin(1),sign(1,jbin(1)-ibin(1))
		do j = ibin(2),jbin(2),sign(1,jbin(2)-ibin(2))
		do k = ibin(3),jbin(3),sign(1,jbin(3)-ibin(3))

			cbin(1) = i; cbin(2) = j; cbin(3) = k

			bintop(:) = (ibin(:)-1)*Fbinsize(:)-halfdomain(:)
			binbot(:) = (ibin(:)-2)*Fbinsize(:)-halfdomain(:)

			if(present(Traction)) then

				!Calculate the plane intersect of line with surfaces of the cube
				Pxt=(/ bintop(1),ri(2)+(rij(2)/rij(1))*(bintop(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(bintop(1)-ri(1))  /)
				Pxb=(/ binbot(1),ri(2)+(rij(2)/rij(1))*(binbot(1)-ri(1)),ri(3)+(rij(3)/rij(1))*(binbot(1)-ri(1))  /)
				Pyt=(/ri(1)+(rij(1)/rij(2))*(bintop(2)-ri(2)), bintop(2),ri(3)+(rij(3)/rij(2))*(bintop(2)-ri(2))  /)
				Pyb=(/ri(1)+(rij(1)/rij(2))*(binbot(2)-ri(2)), binbot(2),ri(3)+(rij(3)/rij(2))*(binbot(2)-ri(2))  /)
				Pzt=(/ri(1)+(rij(1)/rij(3))*(bintop(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(bintop(3)-ri(3)), bintop(3) /)
				Pzb=(/ri(1)+(rij(1)/rij(3))*(binbot(3)-ri(3)), ri(2)+(rij(2)/rij(3))*(binbot(3)-ri(3)), binbot(3) /)

				onfacexb =   (sign(1.d0,binbot(1)- rj(1))   &
							- sign(1.d0,binbot(1)- ri(1)))* &
							 (heaviside(bintop(2)-Pxb(2))   &
							- heaviside(binbot(2)-Pxb(2)))* &
							 (heaviside(bintop(3)-Pxb(3))   & 
							- heaviside(binbot(3)-Pxb(3)))
				onfaceyb =	 (sign(1.d0,binbot(2)- rj(2))	&
							- sign(1.d0,binbot(2)- ri(2)))* &
							 (heaviside(bintop(1)-Pyb(1))	& 
							- heaviside(binbot(1)-Pyb(1)))* &
							 (heaviside(bintop(3)-Pyb(3))	&
							- heaviside(binbot(3)-Pyb(3)))
				onfacezb =	 (sign(1.d0,binbot(3)- rj(3)) 	&
							- sign(1.d0,binbot(3)- ri(3)))* &
							 (heaviside(bintop(1)-Pzb(1))	&
							- heaviside(binbot(1)-Pzb(1)))* &
							 (heaviside(bintop(2)-Pzb(2))	&
							- heaviside(binbot(2)-Pzb(2)))

				onfacext =	 (sign(1.d0,bintop(1)- rj(1))	&
							- sign(1.d0,bintop(1)- ri(1)))* &
							 (heaviside(bintop(2)-Pxt(2))	&
							- heaviside(binbot(2)-Pxt(2)))* &
			            	 (heaviside(bintop(3)-Pxt(3))	&
							- heaviside(binbot(3)-Pxt(3)))
				onfaceyt =	 (sign(1.d0,bintop(2)- rj(2))	&
							- sign(1.d0,bintop(2)- ri(2)))* &
							 (heaviside(bintop(1)-Pyt(1))	&
							- heaviside(binbot(1)-Pyt(1)))* &
							 (heaviside(bintop(3)-Pyt(3))	&
							- heaviside(binbot(3)-Pyt(3)))
				onfacezt =	 (sign(1.d0,bintop(3)- rj(3))	&
							- sign(1.d0,bintop(3)- ri(3)))* &
							 (heaviside(bintop(1)-Pzt(1))	&
							- heaviside(binbot(1)-Pzt(1)))* &
							 (heaviside(bintop(2)-Pzt(2))	&
							- heaviside(binbot(2)-Pzt(2)))

				!Prevent halo molecules from being included but include molecule which have left domain 
				!before rebuild has been triggered.
				if (molnoi .gt. np .or. molnoj .gt. np) then
					if (cbin(1) .lt. 2 .or. cbin(1) .gt. nbins(1)+1) cycle
					if (cbin(2) .lt. 2 .or. cbin(2) .gt. nbins(2)+1) cycle
					if (cbin(3) .lt. 2 .or. cbin(3) .gt. nbins(3)+1) cycle
				endif

				!Stress acting on face over volume
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,1) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,1) + fij(:)*dble(onfacexb)
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,2) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,2) + fij(:)*dble(onfaceyb)
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,3) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,3) + fij(:)*dble(onfacezb)
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,4) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,4) + fij(:)*dble(onfacext)
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,5) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,5) + fij(:)*dble(onfaceyt)
				Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,6) = & 
					Tractionbins(modulo(i,3)+1,modulo(j,3)+1,modulo(k,3)+1,:,6) + fij(:)*dble(onfacezt)	


				!Force applied to volume
				!fsurface(:) = 0.25d0*fij(:)*dble(onfacexb - onfacext)
				!fsurface(:) = 0.25d0*fij(:)*dble(onfaceyb - onfaceyt)
				!fsurface(:) = 0.25d0*fij(:)*dble(onfacezb - onfacezt)

				!Add surface force to current bin
				!Traction(:,1) = Traction(:,1) + 0.25d0*fij(:)*dble(onfacexb)
				!Traction(:,2) = Traction(:,2) + 0.25d0*fij(:)*dble(onfaceyb)
				!Traction(:,3) = Traction(:,3) + 0.25d0*fij(:)*dble(onfacezb)
				!Traction(:,4) = Traction(:,4) + 0.25d0*fij(:)*dble(onfacext)
				!Traction(:,5) = Traction(:,5) + 0.25d0*fij(:)*dble(onfaceyt)
				!Traction(:,6) = Traction(:,6) + 0.25d0*fij(:)*dble(onfacezt)

				!Add for molecule i
				Fsurface = fij(:)* dble((heaviside(bintop(1)-ri(1))-heaviside(binbot(1)-ri(1)))* & 
					  					(heaviside(bintop(2)-ri(2))-heaviside(binbot(2)-ri(2)))* & 
					  					(heaviside(bintop(3)-ri(3))-heaviside(binbot(3)-ri(3)))- & 
					  					(heaviside(bintop(1)-rj(1))-heaviside(binbot(1)-rj(1)))* & 
					  					(heaviside(bintop(2)-rj(2))-heaviside(binbot(2)-rj(2)))* & 
					  					(heaviside(bintop(3)-rj(3))-heaviside(binbot(3)-rj(3))))

				!if (onfaceyb.ne.0.or.onfaceyt.ne.0) print'(9i8)', iter, molnoi, molnoj,np, ibin,onfaceyb,onfaceyt

			else
				!Add for molecule i
				Fsurface = fij(:)* dble((heaviside(bintop(1)-ri(1))-heaviside(binbot(1)-ri(1)))* & 
					  					(heaviside(bintop(2)-ri(2))-heaviside(binbot(2)-ri(2)))* & 
					  					(heaviside(bintop(3)-ri(3))-heaviside(binbot(3)-ri(3)))- & 
					  					(heaviside(bintop(1)-rj(1))-heaviside(binbot(1)-rj(1)))* & 
					  					(heaviside(bintop(2)-rj(2))-heaviside(binbot(2)-rj(2)))* & 
					  					(heaviside(bintop(3)-rj(3))-heaviside(binbot(3)-rj(3))))
			endif

		enddo
		enddo
		enddo

		!Take flux from central bin only
		if (present(Traction))then
			Traction(:,:) = Traction(:,:) +  Tractionbins(modulo(icell,3)+1, 	& 
							     		     			  modulo(jcell,3)+1, 	& 
			    						     			  modulo(kcell,3)+1,:,:)


			if (icell .eq. 5 .and. kcell .eq. 3) then
				print'(3i8,4f10.5)',icell,jcell,kcell, Tractionbins(modulo(icell,3)+1, 	& 
							     		     			  modulo(jcell,3)+1, 	& 
			    						     			  modulo(kcell,3)+1,1,2),& 
											 Tractionbins(modulo(icell,3)+1, 	& 
							     		     			  modulo(jcell,3)+1, 	& 
			    						     			  modulo(kcell,3)+1,1,5) & 
						, Pxyface(icell,jcell,kcell,1,2),Pxyface(icell,jcell,kcell,1,5)
			endif
		endif
			
	end subroutine get_Traction


end subroutine get_CV_force


end module control_volume
