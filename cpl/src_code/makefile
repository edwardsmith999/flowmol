###################################################################
#   _____ _____  _                        _         __ _ _        #
#  / ____|  __ \| |                      | |       / _(_) |       #
# | |    | |__) | |       _ __ ___   __ _| | _____| |_ _| | ___   #
# | |    |  ___/| |      | '_ ` _ \ / _` | |/ / _ \  _| | |/ _ \  #
# | |____| |    | |____  | | | | | | (_| |   <  __/ | | | |  __/  #
#  \_____|_|    |______| |_| |_| |_|\__,_|_|\_\___|_| |_|_|\___|  #
#                                                                 #
###################################################################
# Compiles the coupler routines into a library which can be linked
# to the MD and CFD codes during compilation

#
# useful ID for the top make file of coupler system
#
override NAME := coupler

#
# Coupler files
#
F90_FILES      := coupler_module.f90 coupler.f90 CPL.f90
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
CPL_LIB =  $(LIB_DIR)/libcoupler.a

#
# list of interface modules to be moved in library directory
#
INTERFACE_MODULES :=  ./obj/*.mod # avoids the trouble with cray compiler capitalised names $(O_FILES:.o=.mod)

#
# Make library file and copies to directory
#
$(CPL_LIB): $(O_FILES) 
	@echo "Coupler has changed -- libcoupler is being rebuilt"
	@[ -d $(LIB_DIR) ] || mkdir -p $(LIB_DIR) && ar cr $@ $(O_FILES) && ar s $@ $(O_FILES) && cp $(INTERFACE_MODULES) $(LIB_DIR) 

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
$(OBJ_DIR)/coupler.o : $(addprefix $(OBJ_DIR)/, coupler_module.o)

