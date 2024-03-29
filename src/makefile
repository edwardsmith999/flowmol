############################################################
# __  __ _____                    _         __ _ _         #
# |  \/  |  __ \                  | |       / _(_) |       #
# | \  / | |  | |  _ __ ___   __ _| | _____| |_ _| | ___   #
# | |\/| | |  | | | '_ ` _ \ / _` | |/ / _ \  _| | |/ _ \  #
# | |  | | |__| | | | | | | | (_| |   <  __/ | | | |  __/  #
# |_|  |_|_____/  |_| |_| |_|\__,_|_|\_\___|_| |_|_|\___|  #
# 														   #
############################################################

# application name used in coupler build
override NAME := stream_md

GIT_VERSION := $(shell git rev-parse HEAD)

#=======================================================================
# File definitions
#=======================================================================
F90_MD_FILES = version.f90 minpack.f90 interfaces.f90 bilinear_intersect.f90 functions_lib.f90 intrinsic_interface.f90 modules.f90 linklist.f90 messenger.serial.f90 setup_read_input.f90 control_volume.f90 setup_set_parameters.f90 parallel_io.serial.f90  simulation_compute_forces.f90 simulation_record.f90  molecular_properties.f90 md_coupler_socket.f90 setup_initialise_microstate.f90 setup_initial_record.f90 simulation_checkrebuild.f90 finish_final_record.f90 finish_clear_all.f90 md.f90 md_main.f90 external_forces.f90 simulation_move_particles_vv.f90 simulation_move_particles_lfv.f90

F90_MD_FILES_P = $(patsubst parallel_io.serial.f90,parallel_io.MPI2.f90,$(patsubst messenger.serial.f90,messenger.MPI.f90,$(F90_MD_FILES)))

F90_CUDA_MD_FILES = $(F90_MD_FILES:simulation_compute_forces.f90=simulation_compute_forces_CUDA.f90) 
CU_CUDA_MD_FILES = inter.cu inter2.cu

ASSEMBLY_FILES = Assembly_Heaviside_fn.s

LAPACK_FILES = psudo_inv_lapack.f90

# object files directory
OBJ_DIR := obj

# Use coupler default value
#USE_COUPLER=no 

# path to .inc file the choosen platform. 
# This is needed for build started in this directory
#
MAKE_INC_PATH := ../platforms

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

#
# Get platform to build, compiler names and flags, check the make.inc directory
#
#Check for default file
#ifndef USE_COUPLER
#  $(warning USE_COUPLER not define for code $(NAME), setting uncoupled run)
#  FPP_FLAGS     := -cpp -DUSE_COUPLER=0
#else 
#  FPP_FLAGS_no  := -cpp -DUSE_COUPLER=0 
#  FPP_FLAGS_yes := -cpp -DUSE_COUPLER=1
#  FPP_FLAGS     := $(FPP_FLAGS_$(USE_COUPLER))
#endif

# Linking flags to use
LDFLAGS := $(FFLAGS)
LAPACKFLAGS:=$(filter-out -fno-underscoring,$(FFLAGS))

#Executable name
MD_EXE = parallel_md.exe

#Re-use .f90/.cu/.s file names with .o instead of .f90
O_MD = $(addprefix obj/, $(F90_MD_FILES:.f90=.o))
O_CUDA = $(CU_CUDA_MD_FILES:.cu=.o) $(F90_CUDA_MD_FILES:.f90=.o)
O_CONTINUUM = $(F90_CONTINUUM:.f90=.o)
O_MD_P = $(addprefix obj/,$(F90_MD_FILES_P:.f90=.o))
ifndef NO_ASSMBLY_HEAVISIDES
  O_ASSEMBLY_FILES = $(addprefix obj/,$(ASSEMBLY_FILES:.s=.o))
endif

ifdef USE_LAPACK
  O_LAPACK_FILES = $(addprefix obj/,$(LAPACK_FILES:lapack.f90=.o))
  LDFLAGS=-llapack -lblas  #system libraries
  #LDFLAGS=lapack-3.8.0/lib* #local build
endif

ifdef PROFILE
  FFLAGS += -pg
  LDFLAGS += -pg
endif

# If coupled, link to static library coupler.a
ifeq ($(strip $(USE_COUPLER)),yes)
  #CPL_LIBRARY_PATH = /home/es205/scratch/coupled_turbulent_case/cpl-library/cpl
  CPL_LIBRARY_PATH = /home/es205/codes/cpl/cpl-library/
  CPL_LIB_DIR = $(CPL_LIBRARY_PATH)/lib
  CPL_INC_DIR = $(CPL_LIBRARY_PATH)/include/cpl
  FFLAGS  += -I$(CPL_INC_DIR) # Object generation flags
  LDFLAGS += -L$(CPL_LIB_DIR) -Wl,-rpath=$(CPL_LIB_DIR) -lcpl
endif


CU = nvcc
.SUFFIXES: .exe .o .f90 .cu .inc

#=======================================================================
# Commands
#=======================================================================
default: 
	@echo "Please add flag serial (s), parallel (p) or type help for options"
s:
	@make serial_md.exe
p: version
	@make $(MD_EXE)
p_lapack:
	@make MD_EXE_LAPACK "FPP_FLAGS=-cpp -DUSE_LAPACK=1 -fPIC" "LDFLAGS=lapack-3.8.0/lib*"
p_sys_lapack:
	@make MD_EXE_LAPACK "FPP_FLAGS=-cpp -DUSE_LAPACK=1 -fPIC" "LDFLAGS=-llapack -lblas"
s_cuda:
	@make md_CUDA.exe
s_cuda_emulate:
	@make md_CUDA_emu.exe "FFLAGS = -O0 -g -traceback -CB -fpe0 -g" "FLAGS_CUDA = -deviceemu -G -g"
debug_s_cuda:
	@make md_CUDA.exe "FFLAGS = -O0 -g -traceback -CB -fpe0" "FLAGS_CUDA = -G -g -O0"
optimised_s:
	@make serial_md.exe "FFLAGS = -O3"
optimised_p:
	@make $(MD_EXE) "BUILD=opt"
debug_s:
	@make serial_md.exe "BUILD=debug"
debug_p:
	@make $(MD_EXE) "BUILD=debug"
fastdebug_p:
	@make $(MD_EXE) "FFLAGS = -O3 -g -traceback -CB -fpe0"   # "BUILD=debug"
full_debug_s:
	@make serial_md.exe "FFLAGS = -g -O0 -traceback -CB -check all -warn all,nodec,interfaces -gen_interfaces"
full_debug_p:
	@make $(MD_EXE) "FFLAGS = -g -O0 -traceback -CB -check all -warn all,nodec,interfaces -gen_interfaces"
profile_s:
	@make serial_md.exe "FFLAGS = -O3 -funroll-loops -pg -g -fpp"
memprofile_s:
	@make serial_md.exe "FFLAGS = -O3 -g"; valgrind --tool=cachegrind ./md.exe
profile_p:
	@make $(MD_EXE) "FFLAGS = -O3 -pg -funroll-loops -Jobj -fno-underscoring -cpp -DUSE_LAPACK=1 -fPIC" "LDFLAGS=-llapack -lblas"
version:
	@echo "function version()" > version.f90
	@echo "   character(100) :: version" >> version.f90
	@echo "   version='$(GIT_VERSION)'" >> version.f90
	@echo "end function" >> version.f90

#Make executables
md.exe: $(MD_EXE)
$(MD_EXE): obj $(O_ASSEMBLY_FILES) $(O_LAPACK_FILES) $(O_MD_P)
	@echo "MD executable is being relinked"
	@echo "LDFLAGS = $(LDFLAGS)"
	$(F90) -o $@ $(O_ASSEMBLY_FILES) $(O_LAPACK_FILES) $(O_MD_P) $(LDFLAGS)
MD_EXE_LAPACK: obj $(O_ASSEMBLY_FILES) $(O_LAPACK_FILES) $(O_MD_P)
	@echo "MD executable is being relinked"
	@echo "LDFLAGS = $(LDFLAGS)"
	$(F90) -o $@ $(O_ASSEMBLY_FILES) $(O_LAPACK_FILES) $(O_MD_P) $(LDFLAGS)
serial_md.exe: obj $(O_MD) 
	@$(F90) -o $(MD_EXE) $(O_MD) $(LDFLAGS)  
md_CUDA.exe: obj $(O_CUDA)
	@cd obj; $(F90) -o $(MD_EXE) $(O_CUDA) $(LIBS) ; mv md.exe ..
md_CUDA_emu.exe: obj $(O_CUDA)
	@cd obj; $(F90) -o $(MD_EXE) $(O_CUDA) $(LIBSEMU) ; mv md.exe ..
continuum.exe: obj $(O_CONTINUUM)
	@cd obj; $(F90) -o $(CONTINUUM_EXE) $(O_CONTINUUM) ; mv continuum.exe ..
obj:
	[ -d obj ] || mkdir obj
help:
	@echo "======================================================================================"
	@echo "MOLECULAR DYNAMICS Options"
	@echo "s				Optimised serial build"
	@echo "p				Optimised parallel build"
	@echo "p_lapack			Optimised parallel build with LAPACK (assumes compiled ./lapack-3.8.0 exists)"
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
	$(F90) $(FFLAGS) -c $< -o $@
.cu.o:
	$(CU) $(FLAGS_CUDA) -c $*.cu -o obj/$*.o
obj/%.o: %.s
	$(F90) $(FFLAGS) -c $< -o $@

obj/%.o: %lapack.f90
	$(F90) $(LAPACKFLAGS) -c $< -o $@
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

