
!=======================================================================
! High-level routines for starting and stopping the code
!
! main_init()
! main_free()
! main_restore()
! main_save()
!
!=======================================================================
module cfd_control
        use data_export
        save
 
! implicit none

contains


subroutine main_init()
        use data, only : file_dir

        !---- Make sure MPI Starts----
        call messenger_invoke()

        !--- SYSTEM CALL --- call system("ps -l -u tzaki | grep 'a.out'")
        ! Read archive file

! perhaps a call to an initialisation subroutine for 
! couper data would fit here as an extension to the if ...
!        if ( use_coupling ) then 
!                file_dir="cfd_data/"
!        endif

        call archive_init()
        iunit = iopen()
        open (iunit, file=trim(file_dir)//"archive", form="unformatted", &
              status="old", iostat=ierr)
        if (ierr .ne. 0) stop "An archive file is required"
        call archive_read(iunit, 1)
        close (iclose(iunit))

        ! Read input file
        iunit = iopen()
        open (iunit, file=trim(file_dir)//"input", status="old", iostat=ierr)
        if (ierr == 0) then
                ! Input file is present
                call archive_read(iunit, 0)
                close (iclose(iunit))
        end if

        ! Initialize SSD, data, and mesh first
        ! Then bounaries, and finally, simulation
        !T  call SSDfile_init()
        call data_init()
        call mesh_init()
                call grid_init()
                call anuint()
                call grid_halo()
        call boundaries_init()
        call simulation_init()
        call statistics_init()
        


        return
end subroutine main_init

subroutine main_free()
        use data, only : file_dir

        ! Output summary report
        call data_isRoot(iflag)
        if (iflag == 1) then
                ireport = iopen()
                open (ireport, file=trim(file_dir)//"report")
                call archive_info(ireport)
                call simulation_info(ireport)
                call data_info(ireport)
                close (iclose(ireport))
        end if

        ! Free simulation
        call simulation_free()

        ! Free objects
        call data_free()

        !T Free File system stuff
        !T call SSDfile_free()

        return
end subroutine main_free

subroutine main_restore()

        ! Determine restart condition
        call readInt("irestart", irestart)

        if (irestart .ne. 0) then
                call data_read()
                !T--- Get (ntime, stime)
                call simulation_read()
                call boundaries_read()
                !T--- Read data file
                call ReadRstrt()
        else
                !T--- Read data file
                call ReadRstrt()

                !----------------- Block boundaries -------------------
                !T-- HERE IS THE FIX: IMPORTANT 01/04/2003)
                !T-- Exchange halos in x-direction , then in y-direction
                !T-- x-direction
                call updateBorder_lim(uc, ngz  ,nlx+1,nly  , id_x, 2, i1_u, i2_u, 3)
                call updateBorder_lim(vc, ngz  ,nlx  ,nly+1, id_x, 2, i1_v, i2_v, 3)
                call updateBorder_lim(wc, ngz+1,nlx  ,nly  , id_x, 2, i1_v, i2_v, 2)

                !T-- y-direction
                call updateBorder_lim(uc, ngz  ,nlx+1,nly  , id_y, 3, j1_u, j2_u, 3)
                call updateBorder_lim(vc, ngz  ,nlx  ,nly+1, id_y, 3, j1_v, j2_v, 3)
                call updateBorder_lim(wc, ngz+1,nlx  ,nly  , id_y, 3, j1_u, j2_u, 2)

                call Cart_to_Flux()
        end if

        call statistics_read()

        return
end subroutine main_restore

subroutine main_save()
        use data, only : file_dir
  
        integer iflag,irestart,iunit,ireport      

        ! Free uPrime, vPrime, wPrime
        call fluct_free()

        ! Write data file
        call WriteRstrt()

        ! Write flow analysis file
        call flowAnalysis_blayerinfo()

        ! Clear archive
        call archive_clear()

        ! Archive objects
        !TAZ DANGER ZONE (May be I want "irestart" to be specified in input file)
        !TAZ ANSWER:  input file overwrites archive at read-in
        irestart = 1
        call writeInt("irestart", irestart)
        call data_write()
        call mesh_write()
        call simulation_write()
        call boundaries_write()
        call statistics_write()

        ! Write archive file
        call data_isRoot(iflag)
        if (iflag == 1) then
                iunit = iopen()
                open (iunit, file=trim(file_dir)//"archive", form="unformatted")
                call archive_write(iunit, 1)
                close (iclose(iunit))
        end if

        return
end subroutine main_save

end module cfd_control
