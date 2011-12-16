#
# targets for CFD and MD
#
MD_TARGET  := p
CFD_TARGET := opt

F90 := ftn
#
#  object dir
#
OBJ_DIR = obj
#
# Coupler files
#
F90_COUPLER_FILES      := coupler.f90 coupler_parameters.f90 coupler_internal_cfd.f90 coupler_internal_md.f90 coupler_internal_common.f90 
F90_NULL_COUPLER_FILES := coupler_null.f90
O_COUPLER_FILES        := $(addprefix obj/,$(F90_COUPLER_FILES:.f90=.o))
O_NULL_COUPLER_FILES   := $(addprefix obj_null/,$(F90_NULL_COUPLER_FILES:.f90=.o))

#
#  CFD sector
#
CFD_SRC_PATH := ../../CFD_dCSE/src_code/main_code_MG_ddType_dP
#  Redefine CFD path to couette solver if couette built
clean_couette couette : CFD_SRC_PATH := ../../Couette_serial 
couette : CFD_TARGET := continuum.exe

#
# MD sector
#
MD_SRC_PATH := ../../MD_dCSE/src_code

# get platform to build, compiler names and flags, check the make.inc directory
#
#Check for default file

# hack for coupler flags 
branch := md

ifndef PLATFORM
    PLATFORM := $(shell cat ../../make.inc/PLATFORM_default)
endif
 
ifndef PLATFORM
  platform_list := $(shell ls ../../make.inc | grep  '.inc' | sed 's/.inc//' )
  $(error The PLATFORM variable must be specified. Try one of the following "PLATFORM= $(platform_list)" )
else
  include ../../make.inc/$(PLATFORM).inc
  #Create default file for current platform
  $(shell echo $(PLATFORM) > ../../make.inc/PLATFORM_default)
endif

FLAGS := $(FLAGS_$(COMP)_$(VERSION))

.PHONY. : clean_all clean_cfd clean_md cfd md create_obj_dir cfd_coupler couette couette_md

all :
	@echo "Top level coupled build - Type help for options"

couette_md : couette md 

cfd : 
	cd $(CFD_SRC_PATH) && $(MAKE) USE_COUPLER=yes $(CFD_TARGET)
md  : coupler
	cd $(MD_SRC_PATH)  && $(MAKE) USE_COUPLER=yes $(MD_TARGET)
couette : coupler 
	cd $(CFD_SRC_PATH) && $(MAKE) USE_COUPLER=yes $(CFD_TARGET)

coupler : create_obj_dir $(O_COUPLER_FILES) 

null_coupler :
	$(MAKE) OBJ_DIR=obj_null create_obj_dir $(O_NULL_COUPLER_FILES)

create_obj_dir :
	[ -d $(OBJ_DIR) ] || mkdir $(OBJ_DIR)

clean_all : clean_null_coupler clean_coupler clean_cfd clean_md clean_couette

clean_couette_md : clean_couette clean_md clean_coupler

clean_cfd :
	cd  $(CFD_SRC_PATH) && $(MAKE) clean 
clean_md :
	cd  $(MD_SRC_PATH)  && $(MAKE) clean  
clean_couette :
	cd  $(CFD_SRC_PATH) && $(MAKE) clean
clean_null_coupler :
	$(MAKE) OBJ_DIR=obj_null clean_coupler
clean_coupler :
	rm -rf $(OBJ_DIR) *.mod *.f90~ *__genmod.f90 *__genmod.mod *~ 
help:
	@echo "======================================================================================"
	@echo "TOP LEVEL COUPLED BUILD - individual MD/CFD can be built using make in MD or CFD "
	@echo "directories directly as well as from this folder. Options are as follows:"
	@echo ""
	@echo "couette_md - build simple diffusive couette flow solver coupled with MD code"
	@echo "cfd - Build full DNS CFD code only "
	@echo "md - Build MD code only"
	@echo "null_coupler - build coupler null files to allow MD code to be built in file"
	@echo "GENERAL Options"
	@echo "clean_all - Deletes all .mod, .obj and other temporary files from everywhere"
	@echo "clean_couette_md - Deletes all .mod, .obj and other temporary files from couette & MD"
	@echo "clean_cfd - Deletes all .mod, .obj and other temporary files from cfd"
	@echo "clean_md - Deletes all .mod, .obj and other temporary files from md"
	@echo "clean_couette - Deletes all .mod, .obj and other temporary files from couette"
	@echo "======================================================================================"


#=======================================================================
# Compilation rules
#=======================================================================
$(O_COUPLER_FILES) : obj/%.o : %.f90
	$(F90) $(FLAGS) -c $< -o $@
$(O_NULL_COUPLER_FILES) : obj_null/%.o : %.f90
	$(F90) $(FLAGS) -c $< -o $@

#=======================================================================
# Dependencies
#=======================================================================
obj/coupler.o : obj/coupler_internal_cfd.o obj/coupler_internal_md.o obj/coupler_parameters.o obj/coupler_internal_common.o
obj/coupler_internal_cfd.o obj/coupler_internal_md.o : obj/coupler_parameters.o obj/coupler_internal_common.o

