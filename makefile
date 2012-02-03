#
# Name ( useful in coupler builds)
#
override NAME := COUETTE

#=======================================================================
# File definitions
#=======================================================================

F90_CONTINUUM = continuum_data_export.f90  continuum_modules.f90 continuum_coupler_socket_init$(COUPLER_SOCKET).f90 continuum_messenger.f90 continuum_advance_time.f90 continuum_calculate_flux.f90 continuum_coupler_socket$(COUPLER_SOCKET).f90 continuum_CFL.f90 continuum_finish.f90 continuum_mesh_gen.f90 continuum_read_inputs.f90 continuum_record.f90 continuum_set_BC.f90 continuum_set_parameters.f90 continuum_setup_macrostate.f90 continuum.f90 continuum_main.f90

F90_COUPLED = $(F90_CONTINUUM:continuum_main.f90=) $(F90_MD_FILES:md_main.f90=MD_continuum_main.f90) setup_coupling.f90 coupler.f90

#
# get platform to build, compiler names and flags, check the make.inc directory
#

#
# objec file dir
#

OBJ_DIR := obj

# use coupler ? 
USE_COUPLER      := no

# coupler path
COUPLER_PATH_yes := ../coupler_dCSE/src_code/obj
COUPLER_PATH_no  := ../coupler_dCSE/src_code/obj_null
COUPLER_PATH     := $(COUPLER_PATH_$(USE_COUPLER))

# coupler object files in coupler director
O_COUPLER_no     := $(COUPLER_PATH)/coupler_null.o
O_COUPLER_yes    := $(COUPLER_PATH)/coupler.o $(COUPLER_PATH)/coupler_parameters.o $(COUPLER_PATH)/coupler_internal_common.o $(COUPLER_PATH)/coupler_internal_md.o $(COUPLER_PATH)/coupler_internal_cfd.o $(COUPLER_PATH)/coupler_input_data.o
O_COUPLER        := $(O_COUPLER_$(USE_COUPLER))   

# socket files
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
O_CONTINUUM = $(addprefix obj/,$(F90_CONTINUUM:.f90=.o))
O_COUPLED = $(F90_COUPLED:.f90=.o)

#Make everything depend on parameter files
$(O_CONTINUUM): $(PARAM)
$(O_COUPLED) : $(PARAM)

CU = nvcc
.SUFFIXES: .exe .o .f90 .cu .inc


#=======================================================================
# Commands
#=======================================================================
default: 
	@echo "Please add flag serial (s), parallel (p) or type help for options"
s_continuum:
	@make continuum.exe
continuum.exe: obj $(O_CONTINUUM)
	$(F90) -o $(CONTINUUM_EXE) $(O_CONTINUUM) $(O_COUPLER)
obj:
	[ -d obj ] || mkdir obj
help:
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

