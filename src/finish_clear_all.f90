!-----------------------------------------------------------------------------
!
!                               Clear all
! Deallocates all allocated arrays and destroys all linklists
!
!-----------------------------------------------------------------------------

module module_clear_all

	use computational_constants_MD
	use physical_constants_MD
	use arrays_MD
	use linked_list
	use calculated_properties_MD
	use shear_info_MD
	use polymer_info_MD
	use module_set_parameters, only : velPDF,velPDF_array

end module module_clear_all
!----------------------------------------------------------------------------------

subroutine finish_clear_all
    use module_clear_all
    implicit none

    integer :: i,j,k,n

	call linklist_deallocateall() !Final deallocation of all linked lists
	
	!Deallocate all allocated arrays
	deallocate(r)
	deallocate(v)
	deallocate(a)

	deallocate(potenergymol)
	deallocate(potenergymol_LJ)
	deallocate(virialmol)

	deallocate(seed)
	deallocate(cell%head)
	deallocate(cell%cellnp)
	deallocate(bin%head)
	deallocate(bin%cellnp)

	if (rtrue_flag.eq.1) then
		deallocate(rtrue)
		deallocate(vtrue)
	endif

	!deallocate(rijsum)
	!deallocate(vmagnitude)

	if (ensemble .eq. nvt_DPD) then
		deallocate(theta)
		deallocate(aD)
		deallocate(aR)
	endif

	if (ensemble .eq. tag_move) then
		deallocate(tag)
		deallocate(fix)
		deallocate(rtether)
		deallocate(slidev)
	endif

	if (vPDF_flag .eq. 5) then
		call velPDF%destroy
    elseif (vPDF_flag .ne. 0) then

        do i =1,nbins(1)+2
        do j =1,nbins(2)+2
        do k =1,nbins(3)+2
        do n =1,nd
            call velPDF_array(i,j,k,n)%destroy
        enddo
        enddo
        enddo
        enddo
	endif

	!deallocate(diffusion)
	!deallocate(meandiffusion)

	deallocate(Pxy)
	deallocate(Pxyzero)
	if (pressure_outflag .eq. 1) then
		deallocate(rfmol)
		deallocate(Pxymol)
	endif
	if (pressure_outflag .eq. 2) then
		deallocate(rfbin)
		deallocate(vvbin)
		deallocate(Pxybin)
	endif
	if(viscosity_outflag .eq. 1) then
		deallocate(Pxycorrel)
	endif

	if (allocated(rdf))             deallocate(rdf)	
	if (allocated(rdf3d))           deallocate(rdf3d)	
	if (allocated(ssf))             deallocate(ssf)	
	if (allocated(rdf_hist))        deallocate(rdf_hist)	
	if (allocated(rdf3d_hist))      deallocate(rdf3d_hist)	
	if (allocated(ssf_hist))        deallocate(ssf_hist)	
	if (allocated(mass_flux))		deallocate(mass_flux)
	if (allocated(volume_mass))		deallocate(volume_mass)

	if (allocated(momentum_flux))	deallocate(momentum_flux)
	if (allocated(Pxybin))			deallocate(Pxybin)
	if (allocated(volume_force))	deallocate(volume_force)
	if (allocated(volume_momentum))	deallocate(volume_momentum)
	if (allocated(Pxyface))			deallocate(Pxyface)

	if (allocated(energy_flux))     deallocate(energy_flux)
	if (allocated(Pxyvface))     	deallocate(Pxyvface)


	if (any(periodic.eq.2)) then
		deallocate(mol_wrap_integer)
	endif

	select case (potential_flag)
	case(0)
	case(1)
!		deallocate(monomer)
!		deallocate(bond)
!		deallocate(bondcount)
!		deallocate(potenergymol_FENE)
		if(allocated(etev))   deallocate(etev)
		if(allocated(etev_0)) deallocate(etev_0)
	end select
	
	if (allocated(vmd_intervals))      deallocate(vmd_intervals)
	if (allocated(volume_temperature)) deallocate(volume_temperature)
	if (allocated(slice_momentum))     deallocate(slice_momentum)
	if (allocated(slice_mass))         deallocate(slice_mass)
	if (allocated(planes))			   deallocate(planes)
	if (allocated(Pxy_plane))	       deallocate(Pxy_plane)
	if (nproc .ne. 1) deallocate(procnp)	!Only allocated in parallel code

end subroutine finish_clear_all
