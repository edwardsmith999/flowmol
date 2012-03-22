!=======================================================================
! Routines for reading and writing to temporary storage arrays.
! The SSD (solid-state-device) is a form of storage commonly used
! on the Cray that is cheaper than memory but faster than disk.
!
! Charles Pierce, July 1995
!
! SSDfile_init()
! SSDfile_free()
! SSDRead(ctag, ipart, data, isize)
! SSDReadTag(ctag, data, isize)
! SSDReadIndex(index, data, isize)
! SSDWrite(ctag, ipart, data, isize)
! SSDWriteTag(ctag, data, isize)
! SSDWriteIndex(index, data, isize)
! SSDfile_dataSize(ctag, ipart, isize)
! SSDfile_isDefined(ctag, ipart, iflag)
! SSDfile_indexForTag(ctag, ipart, index, iwrite)*
!

module SSDfile
	parameter (maxTag=16, maxPart=32, maxEntry=64)

	character*16 tag(maxTag)         ! name tags for variables
	integer itag(maxTag, maxPart)    ! pointer to directory entry
	integer ntag                     ! number of tags
	integer nentry                   ! number of entries
	integer isizes(maxEntry)         ! sizes of entries

	! Definitions for memory file
	type rpointer
		real, pointer, dimension(:) :: p
	end type
	type (rpointer) memory(maxEntry)

end module

subroutine SSDfile_init()
	use SSDfile

	tag = ""
	itag = 0
	ntag = 0
	nentry = 0
	isizes = 0

	do i=1,maxEntry
		if (associated(memory(i)%p)) &
			deallocate(memory(i)%p)
	end do

	return
end

subroutine SSDfile_free()
	use SSDfile

	do i=1,maxEntry
		if (associated(memory(i)%p)) &
			deallocate(memory(i)%p)
	end do

	return
end

!=======================================================================
subroutine SSDRead(ctag, ipart, data, isize)
	use SSDfile
	character*(*) ctag

	call SSDfile_indexForTag(ctag, ipart, index, 0)
	if (index == 0) then
		print *, "WARNING: SSDfile: name not found: ", ctag, " ", ipart
		stop
	end if
	call SSDReadIndex(index, data, isize)

	return
end

subroutine SSDReadTag(ctag, data, isize)
	use SSDfile
	character*(*) ctag
	call SSDRead(ctag, 1, data, isize)
	return
end

subroutine SSDReadIndex(index, data, isize)
	use SSDfile
	real data(isize)

	if (index < 1 .or. index > maxEntry) stop "SSDfile: invalid index"
	if (isize .ne. isizes(index)) stop "SSDfile: wrong isize"

	data = memory(index)%p

	return
end

!=======================================================================
subroutine SSDWrite(ctag, ipart, data, isize)
	use SSDfile
	character*(*) ctag
	call SSDfile_indexForTag(ctag, ipart, index, 1)
	call SSDWriteIndex(index, data, isize)
	return
end

subroutine SSDWriteTag(ctag, data, isize)
	use SSDfile
	character*(*) ctag
	call SSDWrite(ctag, 1, data, isize)
	return
end

subroutine SSDWriteIndex(index, data, isize)
	use SSDfile
	real data(isize)

	if (index <=0 .or. index > maxEntry) stop "SSDfile: invalid index"
	if (isize <= 0) stop "SSDfile: isize must be > 0"

	if (isize .ne. isizes(index)) then
		if (associated(memory(index)%p)) deallocate(memory(index)%p)
		allocate(memory(index)%p(isize))
		isizes(index) = isize
	end if

	memory(index)%p = data

	return
end

!=======================================================================
subroutine SSDfile_dataSize(ctag, ipart, isize)
	use SSDfile
	character*(*) ctag

	call SSDfile_indexForTag(ctag, ipart, index, 0)
	if (index .ne. 0) then
		isize = isizes(index)
	else
		print *, "WARNING: SSDfile: name not found: ", ctag, " ", ipart
		stop
	end if

	return
end

subroutine SSDfile_isDefined(ctag, ipart, iflag)
	use SSDfile
	character*(*) ctag

	call SSDfile_indexForTag(ctag, ipart, index, 0)
	if (index .ne. 0) then
		iflag = 1
	else
		iflag = 0
	end if

	return
end

subroutine SSDfile_indexForTag(ctag, ipart, index, iwrite)
	use SSDfile
	character*(*) ctag

	if (ipart <= 0 .or. ipart > maxPart) &
		stop "SSDfile: invalid ipart reference"

	do i=1,ntag
		if (tag(i) == ctag) then
			index = itag(i,ipart)
			if (index == 0 .and. iwrite == 1) then
				nentry = nentry + 1
				index = nentry
				itag(i,ipart) = index
			end if
			return
		end if
	end do

	if (iwrite == 1) then
		! Add new field
		if (ntag < maxTag) then
			ntag = ntag + 1
			tag(ntag) = ctag
			nentry = nentry + 1
			index = nentry
			itag(ntag,ipart) = index
		else
			stop "SSDfile: maxTag insufficient"
		end if
	else
		index = 0
	end if

	return
end

