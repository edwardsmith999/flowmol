!=======================================================================
! General facility for archiving variables
! Variables are referenced by character name tag
!
! Charles Pierce, July 1995
!
! archive_init()
! readInt(ctag, ivalue)
! readFloat(ctag, rvalue)
! readIntArray(ctag, iarray, n)
! readArray(ctag, array, n)
! writeInt(ctag, ivalue)
! writeFloat(ctag, rvalue)
! writeIntArray(ctag, iarray, n)
! writeArray(ctag, array, n)
! archive_arraySize(ctag, isize)
! archive_isDefined(ctag, iflag)
! archive_rename(ctag_old, ctag_new)
! archive_delete(ctag)
! archive_fieldForTag(ctag, ifield, iwrite)*
! archive_clear()
! archive_read(iunit, iformat)
! archive_write(iunit, iformat)
! archive_info(iunit)
!

module archive
	! Max number of fields = 128
	! Field length = 16
	parameter (maxFields=128, fieldLen=16)

	type rpointer
		real, pointer, dimension(:) :: p
	end type

	integer         nfields              ! number of fields
	character*16    tag(maxFields)       ! name identifying field
	character       type(maxFields)      ! data type, 'f' or 'i'
	type (rpointer) value(maxFields)     ! data value, stored as real
end module

!=======================================================================
subroutine archive_init()
	use archive
	call archive_clear()
	return
end

!=======================================================================
subroutine readInt(ctag, ivalue)
	use archive
	character*(*) ctag
	integer       ivalue
	call readIntArray(ctag, ivalue, 1)
	return
end

subroutine readFloat(ctag, rvalue)
	use archive
	character*(*) ctag
	real          rvalue
	call readArray(ctag, rvalue, 1)
	return
end

subroutine readIntArray(ctag, iarray, n)
	use archive
	character*(*) ctag
	integer       iarray(n)

	call archive_fieldForTag(ctag, ifield, 0)
	if (ifield .ne. 0) then
		if (type(ifield) .ne. "i") then
			print *, "archive: ", tag(ifield), " is not integer"
			stop
		end if
		isize = size(value(ifield)%p)
		if (n .ne. isize) then
			print *, "archive: wrong size for ", tag(ifield)
			stop
		end if
		iarray = nint(value(ifield)%p)
	else
		iarray = 0
		print *, "WARNING: archive: name not found: ", ctag
	end if

	return
end

subroutine readArray(ctag, array, n)
	use archive
	character*(*) ctag
	real          array(n)

	call archive_fieldForTag(ctag, ifield, 0)
	if (ifield .ne. 0) then
		if (type(ifield) .ne. "f") then
			print *, "archive: ", tag(ifield), " is not float"
			stop
		end if
		isize = size(value(ifield)%p)
		if (n .ne. isize) then
			print *, "archive: wrong size for ", tag(ifield)
			stop
		end if
		array = value(ifield)%p
	else
		array = 0.
		print *, "WARNING: archive: name not found: ", ctag
	end if

	return
end

!=======================================================================
subroutine writeInt(ctag, ivalue)
	use archive
	character*(*) ctag
	integer       ivalue
	call writeIntArray(ctag, ivalue, 1)
	return
end

subroutine writeFloat(ctag, rvalue)
	use archive
	character*(*) ctag
	real          rvalue
	call writeArray(ctag, rvalue, 1)
	return
end

subroutine writeIntArray(ctag, iarray, n)
	use archive
	character*(*) ctag
	integer       iarray(n)

	call archive_fieldForTag(ctag, ifield, 1)
	if (associated(value(ifield)%p)) &
		deallocate(value(ifield)%p)
	allocate(value(ifield)%p(n))
	value(ifield)%p = float(iarray)
	type(ifield) = "i"

	return
end

subroutine writeArray(ctag, array, n)
	use archive
	character*(*) ctag
	real          array(n)

	call archive_fieldForTag(ctag, ifield, 1)
	if (associated(value(ifield)%p)) &
		deallocate(value(ifield)%p)
	allocate(value(ifield)%p(n))
	value(ifield)%p = array
	type(ifield) = "f"

	return
end

!=======================================================================
subroutine archive_arraySize(ctag, isize)
	use archive
	character*(*) ctag

	call archive_fieldForTag(ctag, ifield, 0)
	if (ifield .ne. 0) then
		isize = size(value(ifield)%p)
	else
		isize = 0
	end if

	return
end

subroutine archive_isDefined(ctag, iflag)
	use archive
	character*(*) ctag

	call archive_fieldForTag(ctag, ifield, 0)
	if (ifield .ne. 0) then
		iflag = 1
	else
		iflag = 0
	end if

	return
end

subroutine archive_rename(ctag_old, ctag_new)
	use archive
	character*(*) ctag_old, ctag_new

	call archive_fieldForTag(ctag_old, ifield, 0)
	if (ifield .ne. 0) then
		tag(ifield) = ctag_new
		ctag_old = "NULL"
	end if

	return
end

subroutine archive_delete(ctag)
	use archive
	character*(*) ctag

	call archive_fieldForTag(ctag, ifield, 0)
	if (ifield .ne. 0) then
		tag(ifield) = "NULL"
		deallocate(value(ifield)%p)
		allocate(value(ifield)%p(1))
		value(ifield)%p = 0.
		ctag = "NULL"
	end if

	return
end

!=======================================================================
! Internal (private) routine
subroutine archive_fieldForTag(ctag, ifield, iwrite)
	use archive
	character*(*) ctag

	do i=1,nfields
		if (tag(i) == ctag) then
			ifield = i
			return
		end if
	end do

	if (iwrite == 1) then
		! Add new field
		if (nfields < maxFields) then
			nfields = nfields + 1
			tag(nfields) = ctag
			ifield = nfields
		else
			stop "archive: maxFields insufficient"
		end if
	else
		ifield = 0
	end if

	return
end

!=======================================================================
subroutine archive_clear()
	use archive

	nfields = 0
	tag = ""
	type = ""
	do i=1,maxFields
		if (associated(value(i)%p)) &
			deallocate(value(i)%p)
	end do

	return
end

!=======================================================================
subroutine archive_read(iunit, iformat)
	use archive
	character*80 buffer
	character*16 ctag
	character*16 charValue
	real         floatValue
	integer      intValue
	character    ctype
	allocatable  array(:), iarray(:)

	select case (iformat)
	case (0) ! Formatted
		ierr = 0
		do while (ierr == 0)
			read (iunit, "(a)", iostat=ierr) buffer
			if (buffer(1:1) .ne. '!' .and. &
			    buffer(1:1) .ne. '#' .and. &
			    trim(buffer) .ne. "" .and. &
			    ierr == 0) then
				read (buffer, *) charValue, ctag
				if (index(charValue, ".") .ne. 0) then
					! Floating-point value
					read (charValue, *) floatValue
					call writeFloat(ctag, floatValue)
				else
					! Integer value
					read (charValue, *) intValue
					call writeInt(ctag, intValue)
				end if
			end if
		end do
	case (1) ! Unformatted
		read (iunit) n
		do i=1,n
			read (iunit) ctag
			read (iunit) ctype
			read (iunit) len
			allocate(array(len))
			read (iunit) array
			if (ctype == "f") then
				call writeArray(ctag, array, len)
			else if (ctype == "i") then
				allocate(iarray(len))
				iarray = nint(array)
				call writeIntArray(ctag, iarray, len)
				deallocate(iarray)
			end if
			deallocate(array)
		end do
	case (2) ! Binary
		! call bread(iunit, data, isize)

	end select

	return
end

!=======================================================================
subroutine archive_write(iunit, iformat)
	use archive

	select case (iformat)
	case (0) ! Formatted
		do i=1,nfields
			isize = size(value(i)%p)
			if (isize == 1) then
				select case (type(i))
				case ("f")
					write (iunit, *) tag(i), value(i)%p
				case ("i")
					write (iunit, *) tag(i), nint(value(i)%p)
				case default
					stop "invalid type specifier"
				end select
			else
				write (iunit, *) tag(i)
				select case (type(i))
				case ("f")
					do k=1,isize
						write (iunit, *) k, value(i)%p(k)
					end do
				case ("i")
					do k=1,isize
						write (iunit, *) k, nint(value(i)%p(k))
					end do
				case default
					stop "invalid type specifier"
				end select
			end if
		end do
	case (1) ! Unformatted
		write (iunit) nfields
		do i=1,nfields
			write (iunit) tag(i)
			write (iunit) type(i)
			write (iunit) size(value(i)%p)
			write (iunit) value(i)%p
		end do
	case (2) ! Binary
		! call bwrite(iunit, data, isize)

	end select

	return
end

!=======================================================================
subroutine archive_info(iunit)
	use archive

	call sectionTitle(iunit, "Variables in Archive")

	do i=1,nfields
		! Only first element of arrays is printed
		n = size(value(i)%p)
		if (n > 1) then
			write (iunit, "($,a)") " *"
		else
			write (iunit, "($,a)") "  "
		end if

		select case (type(i))
		case ("f")
			write (iunit, "($,a10, f10.4)") tag(i), value(i)%p(1)
		case ("i")
			write (iunit, "($,a10, i10)") &
			tag(i), nint(value(i)%p(1))
		case default
			stop "invalid type specifier"
		end select

		if (mod(i,3) == 0 .or. i == nfields) then
			write (iunit, *)
		else
			write (iunit, "($,a)") "  |"
		end if
	end do

	write (iunit, 11) nfields

11  format (/ 2x, "Number of variables: ", i4)

	return
end

