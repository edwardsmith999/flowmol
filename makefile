#
# Name ( useful in coupler builds)
#
override NAME := COUTTE

#=======================================================================
# File definitions
#=======================================================================

F90_MD_FILES = modules.f90 functions_lib.f90 linklist.f90 molecular_properties.f90 messenger.serial.f90 parallel_io.serial.f90 setup_set_parameters.f90 setup_initialise_microstate.f90 setup_initial_record.f90 simulation_compute_forces.f90 simulation_move_particles.f90 simulation_checkrebuild.f90 simulation_record.f90 finish_final_record.f90 finish_clear_all.f90 md.f90 md_main.f90 external_forces.f90

F90_CUDA_MD_FILES = $(F90_MD_FILES:simulation_compute_forces.f90=simulation_compute_forces_CUDA.f90) 

CU_CUDA_MD_FILES = inter.cu inter2.cu

F90_CONTINUUM = data_export.f90  continuum_modules.f90 continuum_coupler_socket_init$(COUPLER_SOCKET).f90 messenger.MPI.noglobal_blocks.f90 continuum_advance_time.f90 continuum_calculate_flux.f90 continuum_coupler_socket$(COUPLER_SOCKET).f90 continuum_CFL.f90 continuum_finish.f90 continuum_mesh_gen.f90 continuum_read_inputs.f90 continuum_record.f90 continuum_set_BC.f90 continuum_set_parameters.f90 continuum_setup_macrostate.f90 continuum.f90 continuum_main.f90

F90_MD_FILES_P = $(patsubst parallel_io.serial.f90,parallel_io.MPI2.f90,$(patsubst messenger.serial.f90,messenger.MPI.f90,$(F90_MD_FILES)))

F90_COUPLED = $(F90_CONTINUUM:continuum_main.f90=) $(F90_MD_FILES:md_main.f90=MD_continuum_main.f90) setup_coupling.f90 coupler.f90

#
# get platform to build, compiler names and flags, check the make.inc directory
#

# use coupler ? 
USE_COUPLER      := no

# coupler path
COUPLER_PATH_yes := ../coupler_dCSE/src_code/obj
COUPLER_PATH_no  := ./obj
COUPLER_PATH     := $(COUPLER_PATH_$(USE_COUPLER))

# coupler object files in coupler director
O_COUPLER_no     :=
O_COUPLER_yes    := $(COUPLER_PATH)/coupler_cfd_global_data.o $(COUPLER_PATH)/coupler_cfd_setup.o $(COUPLER_PATH)/coupler_cfd_communication.o
O_COUPLER        := $(O_COUPLER_$(USE_COUPLER))   

# sockek files
COUPLER_SOCKET_yes :=
COUPLER_SOCKET_no  :=_dummy
COUPLER_SOCKET     := $(COUPLER_SOCKET_$(USE_COUPLER))#

#Check for default file
ifndef PLATFORM
    PLATFORM := $(shell cat ../make.inc/PLATFORM_default)
endif
 
ifndef PLATFORM
  platform_list := $(shell ls ../make.inc | grep  '.inc' | sed 's/.inc//' )
  $(error The PLATFORM variable must be specified. Try one of the following: $(platform_list) )
else
# 
#   get the branch ( cfd or md)
#
  branch = cfd# $(if $(findstring CFD_dCSE,$(PWD)),cfd,$(if $(findstring MD_dCSE,$(PWD)),md,$(error "error in find branch")))
  include ../make.inc/$(PLATFORM).inc
  #Create default file for current platform
  $(shell echo $(PLATFORM) > ../make.inc/PLATFORM_default)
endif

# computed flags with values from PLATFORM files
FFLAGS   := $(FLAGS_$(COMP)_$(VERSION)) -I$(COUPLER_PATH)
LDFLAGS := $(FLAGS)

CONTINUUM_EXE = continuum.exe

#=======================================================================
# Compiler
#=======================================================================

#Re-use f90 file names with .o instead of .f90
O_MD = $(F90_MD_FILES:.f90=.o)
O_CUDA = $(CU_CUDA_MD_FILES:.cu=.o) $(F90_CUDA_MD_FILES:.f90=.o)
O_CONTINUUM = $(addprefix obj/,$(F90_CONTINUUM:.f90=.o))
O_COUPLED = $(F90_COUPLED:.f90=.o)
O_MD_P = $(F90_MD_FILES_P:.f90=.o)

#Make everything depend on parameter files
$(O_MD): $(PARAM)
$(O_CUDA): $(PARAM)
$(O_CONTINUUM): $(PARAM)
$(O_COUPLED) : $(PARAM)
$(O_MD_P): $(PARAM)

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
s_coupled:
	@make coupled.exe 
debug_s_coupled:
	@make coupled.exe 
p:
	@make parallel_md.exe
s_cuda:
	@make md_CUDA.exe
s_cuda_emulate:
	@make md_CUDA_emu.exe 
debug_s_cuda:
	@make md_CUDA.exe 
optimised_s:
	@make md.exe 
optimised_p:
	@make parallel_md.exe 
debug_s:
	@make md.exe 
debug_p:
	@make parallel_md.exe 
full_debug_s:
	@make md.exe 
full_debug_p:
	@make parallel_md.exe 
profile_s:
	@make md.exe 
profile_p:
	@make parallel_md.exe 
md.exe: obj $(O_MD)
	@cd obj; $(F90) -o $(MD_EXE) $(O_MD) ; mv md.exe ..
md_CUDA.exe: obj $(O_CUDA)
	@cd obj; $(F90) -o $(MD_EXE) $(O_CUDA) $(LIBS) ; mv md.exe ..
md_CUDA_emu.exe: obj $(O_CUDA)
	@cd obj; $(F90) -o $(MD_EXE) $(O_CUDA) $(LIBSEMU) ; mv md.exe ..
continuum.exe: obj $(O_CONTINUUM)
	$(F90) -o $(CONTINUUM_EXE) $(O_CONTINUUM) $(O_COUPLER)
coupled.exe: obj $(O_COUPLED)
	$(F90) -o $(COUPLED_EXE) $(O_COUPLED) 
parallel_md.exe: obj $(O_MD_P)
	@cd obj; $(F90) -o $(MD_EXE) $(O_MD_P) ; mv md.exe ..
obj:
	[ -d obj ] || mkdir obj
help:
	@echo "======================================================================================"
	@echo "MOLECULAR DYNAMICS Options"
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
	@echo "======================================================================================"
	@echo "CONTINUUM Options"
	@echo "s_continuum		Optimised serial build"
	@echo "======================================================================================"
	@echo "COUPLED Options"
	@echo "s_coupled		Optimised serial build"
	@echo "debug_s_coupled		Serial build with no optimisation and basic debug flags"
	@echo "======================================================================================"
	@echo "GENERAL Options"
	@echo "clean			Deletes all .mod, .obj and other temporary files"
	@echo "======================================================================================"
clean:
	rm -rf obj *.exe *.mod *.f90~ *__genmod.f90 *__genmod.mod *~ 

#=======================================================================
# Compilation rules
#=======================================================================
obj/%.o : %.f90
	$(F90) -c $(FFLAGS) $< -o $@
.cu.o:
	$(CU) $(FLAGS_CUDA) -c $*.cu -o obj/$*.o

#
# Dependecies
#
#obj/continuum_coupler_socket.o : $(COUPLER_PATH)/obj/%.o