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
	deallocate(Pxy_plane)
	deallocate(Pxycorrel)
	deallocate(volume_mass)
	deallocate(volume_momentum)
	deallocate(momentum_flux)
	deallocate(volume_force)
	deallocate(mass_flux)
	deallocate(potenergymol)
	deallocate(virialmol)
	deallocate(planes)

	if (nproc .ne. 1) deallocate(procnp)	!Only allocated in parallel code

end subroutine finish_clear_all
