#
# useful ID for the top make file of coupler system
#
override NAME := coupler
#
# Coupler files
#

F90_FILES      := coupler_parameters.f90 coupler_module.f90 coupler.f90 coupler_internal_cfd.f90 coupler_internal_md.f90  coupler_input_data.f90
#
# object file dir, is needed by platform
#
OBJ_DIR := ./obj
#
# make.inc path, needed because some .inc files contain include statements
#
MAKE_INC_PATH :=  ../../platforms
#
# pick the build platform from platform directory
#
ifndef PLATFORM
    PLATFORM := $(shell cat $(MAKE_INC_PATH)/PLATFORM_default)
endif

ifndef PLATFORM
  platform_list := $(shell ls $(MAKE_INC_PATH)  | grep  '.inc' | sed 's/.inc//' )
  $(error The PLATFORM variable must be specified. Try one of the following "PLATFORM= $(platform_list)" )
else
  include $(MAKE_INC_PATH)/$(PLATFORM).inc
  #Create default file for current platform
  $(shell echo $(PLATFORM) > $(MAKE_INC_PATH)/PLATFORM_default)
endif

#
# object files
#
O_FILES := $(addprefix $(OBJ_DIR)/, $(F90_FILES:.f90=.o))
#
# library directory
#
LIB_DIR := ../../lib/coupler/$(PLATFORM)/$(BUILD)

#
# list of interface modules to be moved in library directory
#
INTERFACE_MODULES :=  ./obj/*.mod # avoids the trouble with cray compiler capitalised names $(O_FILES:.o=.mod)
#INTERFACE_MODULES :=  *.mod 

$(LIB_DIR)/libcoupler.a : $(O_FILES) 
	[ -d $(LIB_DIR) ] || mkdir -p $(LIB_DIR) && ar cr $@ $(O_FILES) && cp $(INTERFACE_MODULES) $(LIB_DIR) 

$(O_FILES) : | $(OBJ_DIR) 

$(OBJ_DIR) :
	mkdir $(OBJ_DIR)

.PHONY: clean 
clean :
	rm -rf $(O_FILES) *.mod *.a  $(LIB_DIR)/*
help:
	@echo "======================================================================================"
	@echo "Coupler library build "
	@echo "Options are as follows:"
	@echo ""
	@echo "clean - deletes all build files and library "
	@echo "======================================================================================"


#=======================================================================
# Compilation rules
#=======================================================================
$(O_FILES) : $(OBJ_DIR)/%.o : %.f90
	$(F90) $(FFLAGS) -c $< -o $@


#=======================================================================
# Dependencies
#=======================================================================
$(OBJ_DIR)/coupler.o : $(addprefix $(OBJ_DIR)/, coupler_parameters.o coupler_input_data.o)
$(OBJ_DIR)/coupler_internal_cfd.o $(OBJ_DIR)/coupler_internal_md.o :  $(addprefix $(OBJ_DIR)/, coupler_parameters.o coupler_input_data.o) 
$(OBJ_DIR)/coupler_input_data.o :  $(addprefix $(OBJ_DIR)/, coupler_parameters.o )


