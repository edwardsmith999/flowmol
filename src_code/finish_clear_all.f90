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

end module module_clear_all
!----------------------------------------------------------------------------------

subroutine finish_clear_all
use module_clear_all
implicit none

	call linklist_deallocateall !Final deallocation of all linked lists
	call linklist_deallocate_bins
	
	!Deallocate all allocated arrays
	deallocate(initialunitsize)
	deallocate(domain) 
	deallocate(halfdomain)
	deallocate(ncells)
	deallocate(cellsidelength)
	deallocate(halfcellsidelength)
	deallocate(r)
	deallocate(rtrue)
	deallocate(rinitial)
	deallocate(rijsum)
	deallocate(v)
	deallocate(vmagnitude)
	deallocate(a)
	deallocate(theta)
	deallocate(aD)
	deallocate(aR)
	deallocate(vfd_bin)
	deallocate(shell)
	deallocate(RDF)	
	deallocate(normalisedvfd_bin)
	deallocate(seed)
	deallocate(cell%head)
	deallocate(cell%cellnp)
	deallocate(bin%head)
	deallocate(bin%cellnp)
	deallocate(diffusion)
	deallocate(meandiffusion)
	deallocate(tag)
	deallocate(fix)
	deallocate(slidev)
	deallocate(thermostat)
	deallocate(rfmol)
	deallocate(rfbin)
	deallocate(Pxy)
	deallocate(Pxymol)
	deallocate(Pxybin)
	deallocate(Pxyface)
	deallocate(Pxyzero)
	deallocate(Pxycorrel)


	if (allocated(volume_momentum))	deallocate(volume_momentum)
	if (allocated(momentum_flux))	deallocate(momentum_flux)
	if (allocated(volume_force))	deallocate(volume_force)
	if (allocated(volume_mass))		deallocate(volume_mass)
	if (allocated(mass_flux))		deallocate(mass_flux)
	deallocate(potenergymol)
	deallocate(potenergymol_LJ)
	deallocate(virialmol)
	deallocate(mol_wrap_integer)
	
	select case (potential_flag)
	case(0)
	case(1)
!		deallocate(monomer)
!		deallocate(bond)
!		deallocate(bondcount)
!		deallocate(potenergymol_FENE)
!		deallocate(etev)
!		deallocate(etev_0)
	end select
	
	if (allocated(vmd_intervals))      deallocate(vmd_intervals)
	if (allocated(volume_temperature)) deallocate(volume_temperature)
	if (allocated(slice_momentum))     deallocate(slice_momentum)
	if (allocated(slice_momentumbin))  deallocate(slice_momentumbin)
	if (allocated(slice_mass))         deallocate(slice_mass)
	if (allocated(slice_massbin))      deallocate(slice_massbin)
	if (allocated(planes))			   deallocate(planes)
	if (allocated(planes))		       deallocate(Pxy_plane)
	if (nproc .ne. 1) deallocate(procnp)	!Only allocated in parallel code

end subroutine finish_clear_all
