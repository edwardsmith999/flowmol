#=======================================================================
# application name used in coupler build
#=======================================================================

override NAME := stream_md

#=======================================================================
# File definitions
#=======================================================================

F90_MD_FILES =  interfaces.f90 modules.f90 functions_lib.f90 messenger.serial.f90 md_coupler_socket.f90 setup_read_input.f90  parallel_io.serial.f90 linklist.f90 molecular_properties.f90 setup_set_parameters.f90 setup_initialise_microstate.f90 setup_initial_record.f90  simulation_compute_forces.f90 simulation_move_particles.f90 simulation_checkrebuild.f90 simulation_record.f90 finish_final_record.f90 finish_clear_all.f90 md.f90 md_main.f90 external_forces.f90 simulation_move_particles_vv.f90 simulation_move_particles_lfv.f90

F90_MD_FILES_P = $(patsubst parallel_io.serial.f90,parallel_io.MPI2.f90,$(patsubst messenger.serial.f90,messenger.MPI.f90,$(F90_MD_FILES)))

F90_CUDA_MD_FILES = $(F90_MD_FILES:simulation_compute_forces.f90=simulation_compute_forces_CUDA.f90) 
CU_CUDA_MD_FILES = inter.cu inter2.cu

# object files directory
OBJ_DIR := obj

# Use coupler?
USE_COUPLER := no
#
# path to .inc file the choosen platform. 
# This is needed for build started in this directory
#
MAKE_INC_PATH :=  ../../platforms
#
# get platform to build, compiler names and flags, check the make.inc directory
#
#Check for default file
ifndef PLATFORM
    PLATFORM := $(shell cat $(MAKE_INC_PATH)/PLATFORM_default)
endif
 
ifndef PLATFORM
  platform_list := $(shell ls $(MAKE_INC_PATH) | grep  '.inc' | sed 's/.inc//' )
  $(error The PLATFORM variable must be specified. Try one of the following: $(platform_list) )
else
  include $(MAKE_INC_PATH)/$(PLATFORM).inc
  #Create default file for current platform
  $(shell echo $(PLATFORM) > $(MAKE_INC_PATH)/PLATFORM_default)
endif

 LDFLAGS := $(FFLAGS)

ifeq ($(strip $(USE_COUPLER)),yes)
  COUPLER_LIB = ../../lib/coupler/$(PLATFORM)/$(BUILD)
  FFLAGS  +=  -I$(COUPLER_LIB)
  LDFLAGS +=  -L$(COUPLER_LIB) -lcoupler
endif

MD_EXE = md.exe

#Re-use f90 file names with .o instead of .f90
O_MD = $(addprefix obj/, $(F90_MD_FILES:.f90=.o))
O_CUDA = $(CU_CUDA_MD_FILES:.cu=.o) $(F90_CUDA_MD_FILES:.f90=.o)
O_CONTINUUM = $(F90_CONTINUUM:.f90=.o)
O_MD_P = $(addprefix obj/,$(F90_MD_FILES_P:.f90=.o))

CU = nvcc
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
	@make md.exe "FLAGS = -O3"
optimised_p:
	@make parallel_md.exe "BUILD=opt"
debug_s:
	@make md.exe "BUILD=debug"
debug_p:
	@make parallel_md.exe "BUILD=debug"
full_debug_s:
	@make md.exe "FLAGS = -O0 -traceback -CB -check all -warn all,nodec,interfaces -gen_interfaces"
full_debug_p:
	@make parallel_md.exe "FLAGS = -O0 -traceback -CB -check all -warn all,nodec,interfaces -gen_interfaces"
profile_s:
	@make md.exe "FLAGS = -O0 -pg"
profile_p:
	@make parallel_md.exe "FFLAGS = -mpe=mpilog "
md.exe: obj $(O_MD)
	@$(F90) -o $(MD_EXE) $(O_MD) $(LDFLAGS)  
md_CUDA.exe: obj $(O_CUDA)
	@cd obj; $(F90) -o $(MD_EXE) $(O_CUDA) $(LIBS) ; mv md.exe ..
md_CUDA_emu.exe: obj $(O_CUDA)
	@cd obj; $(F90) -o $(MD_EXE) $(O_CUDA) $(LIBSEMU) ; mv md.exe ..
continuum.exe: obj $(O_CONTINUUM)
	@cd obj; $(F90) -o $(CONTINUUM_EXE) $(O_CONTINUUM) ; mv continuum.exe ..
parallel_md.exe: obj $(O_MD_P)
	@$(F90) -o $(MD_EXE) $(O_MD_P)  $(LDFLAGS)
obj:
	[ -d obj ] || mkdir obj
help:
	@echo "======================================================================================"
	@echo "MOLECULAR DYNAMICS Options"
	@echo "s				Optimised serial build"
	@echo "p				Optimised parallel build"
	@echo "s_cuda			Optimised serial build with CUDA GPGPU force optimisation"
	@echo "s_cuda_emulate	Emulate CUDA GPU code on the CPU for debugging"
	@echo "optimised_s		Serial build optimised for intel proccessors"
	@echo "optimised_p		Parallel build optimised for intel proccessors"
	@echo "debug_s			Serial build with no optimisation and basic debug flags"
	@echo "debug_p			Parallel build with no optimisation and basic debug flags"
	@echo "full_debug_s		Serial build with no optimisation and extended debug flags"
	@echo "full_debug_p		Parallel build with no optimisation and extended debug flags"
	@echo "profile_s		Serial build with gprof or other profiler data generated"
	@echo "======================================================================================"
	@echo "GENERAL Options"
	@echo "clean			Deletes all .mod, .obj and other temporary files"
	@echo "======================================================================================"
clean:
	rm -rf obj *.exe *.mod *.f90~ *__genmod.f90 *__genmod.mod *~ 

#=======================================================================
# Compilation rules
#=======================================================================
obj/%.o: %.f90
	$(F90) $(FPP_FLAGS) $(FFLAGS) -c $< -o $@
.cu.o:
	$(CU) $(FLAGS_CUDA) -c $*.cu -o obj/$*.o
#
# Dependecies
#

#obj/parallel_io.MPI2.o : obj/messenger.MPI.o
#obj/messenger.MPI.o : obj/setup_inputs.o obj/coupler_md_global_data.o obj/modules.o
#obj/coupler_md_global_data.o : obj/modules.o
#obj/finish_final_record.o : obj/setup_inputs.o 
#obj/md.o : obj/coupler_md_communication.o obj/coupler_md_setup.o
#obj/setup_inputs.o : obj/coupler_md_global_data.o
#obj/setup_set_parameters.o : obj_null/coupler_null.o
#obj/coupler_md_setup.o : obj/coupler_md_global_data.o obj/messenger.MPI.o
#obj/coupler_md_communication.o : obj/coupler_md_setup.o obj/modules.o

