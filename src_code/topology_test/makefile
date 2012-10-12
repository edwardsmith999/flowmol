SOURCES= coupler_parameters.f90 coupler_module.f90 tocoupler.f90 tests.f90  

a.out: *.f90
	mpif90 -O3 -mpe=mpilog $(SOURCES)

debug: *.f90
	mpif90 -O0 -g -traceback -check all $(SOURCES)

mpidebug: *.f90
	mpif90 -O0 -g -traceback -check all -mpe=mpicheck -mpe=mpitrace $(SOURCES)

clean:
	rm a.out *.mod info_*
