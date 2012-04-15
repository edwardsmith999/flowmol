#
# Name ( useful in coupler builds)
#
override NAME := serial_couette

#=======================================================================
# File definitions
#=======================================================================

F90_CONTINUUM = continuum_data_export.f90  continuum_modules.f90 continuum_messenger.f90  continuum_coupler_socket.f90  continuum_advance_time.f90 continuum_calculate_flux.f90 continuum_CFL.f90 continuum_finish.f90 continuum_mesh_gen.f90 continuum_read_inputs.f90 continuum_record.f90 continuum_set_BC.f90 continuum_set_parameters.f90 continuum_setup_macrostate.f90 continuum.f90 continuum_main.f90

F90_COUPLED = $(F90_CONTINUUM:continuum_main.f90=) $(F90_MD_FILES:md_main.f90=MD_continuum_main.f90) setup_coupling.f90 coupler.f90

# needed in platform
OBJ_DIR := obj
USE_COUPLER := no
#
# make.inc path, needed because some .inc files contain include statements
#
MAKE_INC_PATH :=  ../platforms
#
# get te platform to build: compiler names and flags, see the make.inc directory
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

LDFLAGES := $(FFLAGS)

ifeq ($(strip $(USE_COUPLER)),yes)
  COUPLER_LIB = ../lib/coupler/$(PLATFORM)/$(BUILD)
  FFLAGS  +=  -I$(COUPLER_LIB)
  LDFLAGS +=  -L$(COUPLER_LIB) -lcoupler
endif


CONTINUUM_EXE = continuum.exe

#=======================================================================
# Compiler
#=======================================================================

#Re-use f90 file names with .o instead of .f90
O_CONTINUUM = $(addprefix $(OBJ_DIR)/,$(F90_CONTINUUM:.f90=.o))
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
continuum.exe: $(O_CONTINUUM)
	$(F90) -o $(CONTINUUM_EXE)  $(O_CONTINUUM)  $(LDFLAGS)
$(O_CONTINUUM): | $(OBJ_DIR)
$(OBJ_DIR) :
	mkdir -p $@
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
	rm -rf $(OBJ_DIR) *.exe *.mod *.f90~ *__genmod.f90 *__genmod.mod *~ 

#=======================================================================
# Compilation rules
#=======================================================================
$(OBJ_DIR)/%.o : %.f90
	$(F90) -c $(FPP_FLAGS) $(FFLAGS) $< -o $@
.cu.o:
	$(CU) $(FLAGS_CUDA) -c $*.cu -o obj/$*.o

#
# Dependecies
#

