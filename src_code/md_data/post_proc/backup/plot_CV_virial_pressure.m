%Read Header file
read_header

for i=1:Nsteps
    [velocity_snapshot_t(i),velocity_flux_t(i),pressure_surface_t(i)] = read_vflux(i,resultfile_dir,globalnbins,nd)
end