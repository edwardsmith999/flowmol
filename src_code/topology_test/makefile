SOURCES=  ../coupler_module.f90 ../coupler.f90 ./test_coupler.f90

a.out: *.f90
	mpif90 -O3 $(SOURCES)

profile: *.f90
	mpif90 -O3 -pg $(SOURCES)

debug: *.f90
	mpif90 -O0 -g -traceback -check all $(SOURCES)

mpidebug: *.f90
	mpif90 -O0 -g -traceback -check all -mpe=mpilog -mpe=mpicheck  $(SOURCES)

clean:
	rm a.out *.mod info_*
