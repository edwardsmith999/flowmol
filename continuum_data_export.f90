module continuum_data_export
        implicit none

        integer, parameter :: npx = 1, npy = 1, npz = 1, nproc = npx * npy * npz
	integer icoord(3,nproc)           ! proc grid coordinates

        character(len=64) file_dir        ! output subdirectory

end module continuum_data_export
