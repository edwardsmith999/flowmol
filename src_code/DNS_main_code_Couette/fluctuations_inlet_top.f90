!=======================================================================
! Reading inflow data from a database file
! and interpolating in time and space to the mesh
!
! fluct_init(Uin)
! fluct_read()
! fluct_write()
! fluct_Time(dt)
! fluct_free()

! fluct_random(Uin)
! fluct_OSdiscrete(Uin)
! fluct_OScontinuous(Uin)
! fluct_OS2D3D(Uin)
! fluct_freestream(Uin)

! fluct_getomega(omega)
!

module fluct
	use data_export
	use mesh_export

	integer	inflowFluct					! Type of inflow fluctuation
	real	ftime						! Time since start of fluctuations
	real 	dftime						! Current timestep
end module

!=======================================================================
subroutine fluct_init(Uin)
	use fluct
        real Uin(0:nly+1, 0:ngz+1, 3)

	! Type of fluctuations to superimpose on Uin
	call readInt("inflowFluct", inflowFluct)

	select case (inflowFluct)
	case (0) ! No fluctuations
!	case (1) ! Random fluctuations
!		! Initialize random generator
!		iseed = 7
!		call randomSeed(iseed)
!	case (2) ! OS discrete mode fluctuations (in blayer)
!	case (3) ! OS continuous mode fluctuations (in freestream)
!	case (4) ! freestream turbulence of intensity Tu
!	case (5) ! 2D + 3D disturbances
!	case (6) ! Two continuous_mode disturbances
!	case (6) ! Two 2D disturbances
	end select

	! Some fluctutation types need to know the time
	select case (inflowFluct)
	case (1:) ! No fluctuations
		! Send current time to freestream.f90
		call archive_isDefined("time",iflag)
		if (iflag == 1) then
			call readFloat("time", ftime)
		else
			ftime = 0.
		endif
		! Send current timestep to freestream.f90
		call archive_isDefined("dt",iflag)
		if (iflag == 1) then
			call readFloat("dt", dftime)
		else
			dftime = 0.1
			print*,'No dt written to database yet'
		endif
	end select


	return
end



subroutine fluct_read()
	use fluct

	return
end



subroutine fluct_write()
	use fluct

	call writeInt("inflowFluct", inflowFluct)
	select case (inflowFluct)
	case (0:3)
	case (4)
	end select

	return
end

subroutine fluct_Time(deltaT)
	use fluct

	ftime = ftime + deltaT

	return
end


subroutine fluct_free()
	use fluct

	select case (inflowFluct)
	case (4) ! freestream turbulence of intensity Tu
	end select

	return
end
