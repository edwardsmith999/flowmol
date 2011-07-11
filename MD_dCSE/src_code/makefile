#=======================================================================
# File definitions
#=======================================================================

F90_MD_FILES = modules.f90 functions_lib.f90 linklist.f90 molecular_properties.f90 messenger.serial.f90 parallel_io.serial.f90 setup_inputs.f90 setup_set_parameters.f90 setup_initialise_microstate.f90 setup_initial_record.f90 simulation_compute_forces.f90 simulation_move_particles.f90 simulation_checkrebuild.f90 simulation_record.f90 finish_final_record.f90 finish_clear_all.f90 md.f90 md_main.f90 external_forces.f90

F90_CUDA_MD_FILES = modules.f90 functions_lib.f90 linklist.f90 molecular_properties.f90 messenger.serial.f90 parallel_io.serial.f90 setup_inputs.f90 setup_set_parameters.f90 setup_initialise_microstate.f90 setup_initial_record.f90 simulation_compute_forces_CUDA.f90 simulation_move_particles.f90 simulation_checkrebuild.f90 simulation_record.f90 finish_final_record.f90 finish_clear_all.f90 md.f90 md_main.f90 

CU_CUDA_MD_FILES = inter.cu inter2.cu

F90_CONTINUUM = continuum_modules.f90 continuum_advance_time.f90 continuum_calculate_flux.f90 continuum_CFL.f90 continuum_finish.f90 continuum_initial_record.f90 continuum_mesh_gen.f90 continuum_read_inputs.f90 continuum_record.f90 continuum_set_BC.f90 continuum_set_parameters.f90 continuum_setup_macrostate.f90 continuum.f90 continuum_main.f90

F90_MD_FILES_P = modules.f90 functions_lib.f90 linklist.f90 molecular_properties.f90 external_forces.f90 messenger.MPI.f90 parallel_io.MPI2.f90 setup_inputs.f90 setup_set_parameters.f90 setup_initialise_microstate.f90 setup_initial_record.f90 simulation_compute_forces.f90 simulation_move_particles.f90 simulation_checkrebuild.f90 simulation_record.f90 finish_final_record.f90 finish_clear_all.f90 md.f90 md_main.f90

LIBS = -L /usr/local/cuda/lib64/ -lcudart
LIBSEMU	= -L /usr/local/cuda/lib64/ -lcudartemu
PARAM = param.inc
FLAGS =  -O3 -funroll-loops 
FLAGS_CUDA = -O3 --ptxas-options=-v -arch=sm_11 -pg -g  
MD_EXE = md.exe
CONTINUUM_EXE = continuum.exe

#=======================================================================
# Compiler
#=======================================================================

#Re-use f90 file names with .o instead of .f90
O_MD = $(F90_MD_FILES:.f90=.o)
O_CUDA = $(CU_CUDA_MD_FILES:.cu=.o) $(F90_CUDA_MD_FILES:.f90=.o)
O_CONTINUUM = $(F90_CONTINUUM:.f90=.o)
O_MD_P = $(F90_MD_FILES_P:.f90=.o)

#Make everything depend on parameter files
$(O_MD): $(PARAM)
$(O_CUDA): $(PARAM)
$(O_CONTINUUM): $(PARAM)
$(O_MD_P): $(PARAM)

F90 =  mpif90
CU = nvcc
VPATH = obj
.SUFFIXES: .exe .o .f90 .cu .inc

#=======================================================================
# Commands
#=======================================================================
default: 
	@echo "Please add flag serial (s), parallel (p) or type help for options"
s:
	@make md.exe
s_continuum:
	@make continuum.exe
p:
	@make parallel_md.exe
s_cuda:
	@make md_CUDA.exe
s_cuda_emulate:
	@make md_CUDA_emu.exe "FLAGS = -O0 -g -traceback -CB -fpe0 -g" "FLAGS_CUDA = -deviceemu -G -g"
debug_s_cuda:
	@make md_CUDA.exe "FLAGS = -O0 -g -traceback -CB -fpe0" "FLAGS_CUDA = -G -g -O0"
optimised_s:
	@make md.exe "FLAGS = -fast -funroll-loops -msse4.1 -ip"
optimised_p:
	@make parallel_md.exe "FLAGS = -fast -funroll-loops -msse4.1"
debug_s:
	@make md.exe "FLAGS = -O0 -g -traceback -CB -fpe0"
debug_p:
	@make parallel_md.exe "FLAGS = -O0 -g -traceback -CB"
full_debug_s:
	@make md.exe "FLAGS = -O0 -traceback -CB -check all -warn all,nodec,interfaces -gen_interfaces"
full_debug_p:
	@make parallel_md.exe "FLAGS = -O0 -traceback -CB -check all -warn all,nodec,interfaces -gen_interfaces"
profile_s:
	@make md.exe "FLAGS = -O0 -pg"
profile_p:
	@make parallel_md.exe "FLAGS = -mpe=mpilog "
md.exe: obj $(O_MD)
	@cd obj; $(F90) -o $(MD_EXE) $(O_MD) ; mv md.exe ..
md_CUDA.exe: obj $(O_CUDA)
	@cd obj; $(F90) -o $(MD_EXE) $(O_CUDA) $(LIBS) ; mv md.exe ..
md_CUDA_emu.exe: obj $(O_CUDA)
	@cd obj; $(F90) -o $(MD_EXE) $(O_CUDA) $(LIBSEMU) ; mv md.exe ..
continuum.exe: obj $(O_CONTINUUM)
	@cd obj; $(F90) -o $(CONTINUUM_EXE) $(O_CONTINUUM) ; mv continuum.exe ..
parallel_md.exe: obj $(O_MD_P)
	@cd obj; $(F90) -o $(MD_EXE) $(O_MD_P) ; mv md.exe ..
obj:
	mkdir obj
help:
	@echo "Options"
	@echo "s			Optimised serial build"
	@echo "p			Optimised parallel build"
	@echo "s_cuda			Optimised serial build with CUDA GPGPU force optimisation"
	@echo "s_cuda_emulate		Emulate CUDA GPU code on the CPU for debugging"
	@echo "optimised_s		Serial build optimised for intel proccessors"
	@echo "optimised_p		Parallel build optimised for intel proccessors"
	@echo "debug_s			Serial build with no optimisation and basic debug flags"
	@echo "debug_p			Parallel build with no optimisation and basic debug flags"
	@echo "full_debug_s		Serial build with no optimisation and extended debug flags"
	@echo "full_debug_p		Parallel build with no optimisation and extended debug flags"
	@echo "profile_s		Serial build with gprof or other profiler data generated"
	@echo "clean			Deletes all .mod, .obj and other temporary files"
clean:
	rm -rf obj *.exe *.mod *.f90~ *__genmod.f90 *__genmod.mod *~ 

#=======================================================================
# Compilation rules
#=======================================================================
.f90.o:
	$(F90) $(FLAGS) -c $*.f90 -o obj/$*.o
.cu.o:
	$(CU) $(FLAGS_CUDA) -c $*.cu -o obj/$*.o
