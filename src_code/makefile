#
# targets for CFD and MD
#
MD_TARGET  := p
CFD_TARGET := opt

F90 := ftn

#
# Define CFD coupler files
#
F90_CFD_FILES := coupler_cfd_global_data.f90 coupler_cfd_setup.f90 coupler_cfd_communication.f90
O_CFD         := $(addprefix obj/,$(F90_CFD_FILES:.f90=.o))

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
md  :
	cd $(MD_SRC_PATH)  && $(MAKE) USE_COUPLER=yes $(MD_TARGET)
couette : cfd_coupler 
	cd $(CFD_SRC_PATH) && $(MAKE) USE_COUPLER=yes $(CFD_TARGET)

cfd_coupler : create_obj_dir $(O_CFD) 

create_obj_dir :
	[ -d obj ] || mkdir obj

clean_all : clean_cfd clean_md clean_couette clean_cfd_coupler

clean_couette_md : clean_couette clean_md clean_cfd_coupler

clean_cfd :
	cd  $(CFD_SRC_PATH) && $(MAKE) clean 
clean_md :
	cd  $(MD_SRC_PATH)  && $(MAKE) clean  
clean_couette :
	cd  $(CFD_SRC_PATH) && $(MAKE) clean
clean_cfd_coupler :
	rm -rf obj *.mod *.f90~ *__genmod.f90 *__genmod.mod *~ 
help:
	@echo "======================================================================================"
	@echo "TOP LEVEL COUPLED BUILD - individual MD/CFD can be built using make in MD or CFD "
	@echo "directories directly as well as from this folder. Options are as follows:"
	@echo ""
	@echo "couette_md - build simple diffusive couette flow solver coupled with MD code"
	@echo "cfd - Build full DNS CFD code only "
	@echo "md - Build MD code only"
	@echo "GENERAL Options"
	@echo "clean_all - Deletes all .mod, .obj and other temporary files from everywhere"
	@echo "clean_couette_md - Deletes all .mod, .obj and other temporary files from couette & MD"
	@echo "clean_cfd - Deletes all .mod, .obj and other temporary files from cfd"
	@echo "clean_md - Deletes all .mod, .obj and other temporary files from md"
	@echo "clean_couette - Deletes all .mod, .obj and other temporary files from couette"
	@echo "clean_cfd_coupler - Deletes all .mod, .obj and temporary files from cfd and coupler"
	@echo "======================================================================================"


#=======================================================================
# Compilation rules
#=======================================================================
$(O_CFD): obj/%.o : %.f90
	$(F90) $(FLAGS) -c $< -o $@
