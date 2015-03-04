cleanables = a.out obj/ *.mod *.pyc fort.* *.a info_*
f90_files = coupler.f90 coupler_module.f90 CPL.f90 
test_src = test_coupler.f90
flags = -O0 -traceback -g -check all
objs = $(addprefix obj/, $(f90_files:.f90=.o))
test_obj = $(addprefix obj/, $(test_src:.f90=.o))

## Compilation rules
a.out: obj $(objs) $(test_obj)
	ar cr libcoupler.a obj/*.o
	mpif90 -o a.out $(test_obj) -L./ -lcoupler
obj:
	[ -d obj ] || mkdir obj
obj/%.o: %.f90
	mpif90 $(flags) -c $< -o $@
clean:
	rm -rf $(cleanables)

## Dependencies
$(test_obj): $(test_src)
	mpif90 $(flags) -c $(test_src) -o $(test_obj)	
obj/coupler.o: obj/coupler_module.o 
obj/CPL.o: obj/coupler_module.o obj/coupler.o
