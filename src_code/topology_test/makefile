SOURCES= coupler_parameters.f90 coupler_module.f90 coupler_input_data.f90 coupler_internal_cfd.f90 coupler_internal_md.f90 coupler.f90 test_coupler.f90

a.out: *.f90
	mpif90 -O3 $(SOURCES)

debug: *.f90
	mpif90 -O0 -g -traceback -check all $(SOURCES)

mpidebug: *.f90
	mpif90 -O0 -g -traceback -check all -mpe=mpilog mpe=mpicheck -mpe=mpitrace $(SOURCES)

clean:
	rm a.out *.mod info_*
