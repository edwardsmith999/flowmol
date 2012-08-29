!=================================================================
! Copy BC from (uc,vc,wc) to border arrays to save space
!
! CopyBC_out()
! CopyBC_in()
!

module CopyBC
	use data_export
	real dummy		! A module must at least have one variable
end module

!=======================================================================
subroutine CopyBC_out()
	use CopyBC
		!------
		! Move boundary data from (uc,vc,wc) to (--bc) for north and south boundaries.
		!------
		if (jblock.eq.1) then
			ucbcs(:,:,1)=uc(:,:,0    )
			ucbcs(:,:,2)=uc(:,:,1    )

			vcbcs(:,:,1)=vc(:,:,0    )
			vcbcs(:,:,2)=vc(:,:,1    )

			wcbcs(:,:,1)=wc(:,:,0    )
			wcbcs(:,:,2)=wc(:,:,1    )
		end if

		if (jblock.eq.npy) then
			ucbcn(:,:,1)=uc(:,:,nlyb-1)
			ucbcn(:,:,2)=uc(:,:,nlyb  )

			vcbcn(:,:,1)=vc(:,:,nlyb  )
			vcbcn(:,:,2)=vc(:,:,nlyb+1)

			wcbcn(:,:,1)=wc(:,:,nlyb-1)
			wcbcn(:,:,2)=wc(:,:,nlyb  )
		end if

		!--------
		! Move boundary data from (uc,vc,wc) to (--bc) for east and west boundaries.
		!--------
		!if (iblock.eq.1) then
		!	ucbcw(:,:,1)=uc(:,0    ,:)
		!	ucbcw(:,:,2)=uc(:,1    ,:) 
		!	ucbcw(:,:,3)=uc(:,2    ,:)
		!
		!	vcbcw(:,:,1)=vc(:,0    ,:)
		!	vcbcw(:,:,2)=vc(:,1    ,:)
		!	vcbcw(:,:,3)=vc(:,2    ,:)
		!
		!	wcbcw(:,:,1)=wc(:,0    ,:)
		!	wcbcw(:,:,2)=wc(:,1    ,:)
		!end if
       		
		!if (iblock.eq.npx) then
		!	ucbce(:,:,1)=uc(:,nlxb-1,:)
		!	ucbce(:,:,2)=uc(:,nlxb  ,:)
		!	ucbce(:,:,3)=uc(:,nlxb+1,:)
		!
		!	vcbce(:,:,1)=vc(:,nlxb-2,:)
		!	vcbce(:,:,2)=vc(:,nlxb-1,:)
		!	vcbce(:,:,3)=vc(:,nlxb  ,:)
		!
		!	wcbce(:,:,1)=wc(:,nlxb-1,:)
		!	wcbce(:,:,2)=wc(:,nlxb  ,:)
		!end if

		!------
		!Move boundary data from (uc,vc,wc) to (--bc) for front and back boundaries.
		!------
		ucbcb(:,:,1)=uc(1    ,:,:)
		ucbcf(:,:,1)=uc(ngz-1,:,:)

		vcbcb(:,:,1)=vc(1    ,:,:)
		vcbcf(:,:,1)=vc(ngz-1,:,:)

		wcbcb(:,:,1)=wc(1  ,:,:)
		wcbcf(:,:,1)=wc(ngz,:,:)

	return
end


subroutine CopyBC_in()
	use CopyBC
		!------
		! Move boundary data from (--bc) to (uc,vc,wc) for north and south boundaries.
		!------
		if (jblock.eq.1) then
			uc(:,:,0    ) = ucbcs(:,:,1)
			uc(:,:,1    ) = ucbcs(:,:,2)

			vc(:,:,0    ) = vcbcs(:,:,1)
			vc(:,:,1    ) = vcbcs(:,:,2)

			wc(:,:,0    ) = wcbcs(:,:,1)
			wc(:,:,1    ) = wcbcs(:,:,2)
		end if

		if (jblock.eq.npy) then
			uc(:,:,nlyb-1) = ucbcn(:,:,1)
			uc(:,:,nlyb  ) = ucbcn(:,:,2)

			vc(:,:,nlyb  ) = vcbcn(:,:,1)
			vc(:,:,nlyb+1) = vcbcn(:,:,2)

			wc(:,:,nlyb-1) = wcbcn(:,:,1)
			wc(:,:,nlyb  ) = wcbcn(:,:,2)
		end if

		!--------
		! Move boundary data from (--bc) to (uc,vc,wc) for east and west boundaries.
		!--------
		!if (iblock.eq.1) then
		!	uc(:,0    ,:) = ucbcw(:,:,1)
		!	uc(:,1    ,:) = ucbcw(:,:,2)
		!	uc(:,2    ,:) = ucbcw(:,:,3)
		!
		!	vc(:,0    ,:) = vcbcw(:,:,1)
		!	vc(:,1    ,:) = vcbcw(:,:,2)
		!	vc(:,2    ,:) = vcbcw(:,:,3)
		!
		!	wc(:,0    ,:) = wcbcw(:,:,1)
		!	wc(:,1    ,:) = wcbcw(:,:,2)
		!end if
		
		!if (iblock.eq.npx) then
		!	uc(:,nlxb-1,:) = ucbce(:,:,1)
		!	uc(:,nlxb  ,:) = ucbce(:,:,2)
		!	uc(:,nlxb+1,:) = ucbce(:,:,3)
		!
		!	vc(:,nlxb-2,:) = vcbce(:,:,1)
		!	vc(:,nlxb-1,:) = vcbce(:,:,2)
		!	vc(:,nlxb  ,:) = vcbce(:,:,3)
		!
		!	wc(:,nlxb-1,:) = wcbce(:,:,1)
		!	wc(:,nlxb  ,:) = wcbce(:,:,2)
		!end if

		!------
		!Move boundary data from (--bc) to (uc,vc,wc) for front and back boundaries.
		!------
		uc(1    ,:,:) = ucbcb(:,:,1)
		uc(ngz-1,:,:) = ucbcf(:,:,1)

		vc(1    ,:,:) = vcbcb(:,:,1)
		vc(ngz-1,:,:) = vcbcf(:,:,1)

		wc(1  ,:,:) = wcbcb(:,:,1)
		wc(ngz,:,:) = wcbcf(:,:,1)

	return
end
