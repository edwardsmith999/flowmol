#
# targets for CFD and MD
#
MD_TARGET  := p
CFD_TARGET := opt

F90 := ftn

#
# CFD coupler files
#
F90_CFD_FILES := coupler_cfd_global_data.f90 coupler_cfd_setup.f90 coupler_cfd_communication.f90
O_CFD         := $(addprefix obj/,$(F90_CFD_FILES:.f90=.o))
#
#  CFD sector
#
CFD_SRC_PATH := ../../CFD_dCSE/src_code/main_code_MG_ddType_dP
#
# MD sector
#
MD_SRC_PATH := ../../MD_dCSE/src_code

branch := cfd
include ../../make.inc/$(PLATFORM).inc
FLAGS := $(FLAGS_$(COMP)_$(VERSION))

.PHONY. : clean_all clean_cfd clean_md cfd md create_obj_dir cfd_coupler couette couette_md

all :
	@echo "pick a task"

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

clean_couette couette : CFD_SRC_PATH := ../../Couette_serial 
couette : CFD_TARGET := continuum.exe


clean_all : clean_cfd clean_md clean_couette clean_cfd_coupler

clean_couette_md : clean_couete clean_md clean_cfd_coupler

clean_cfd :
	cd  $(CFD_SRC_PATH) && $(MAKE) clean 
clean_md :
	cd $(MD_SRC_PATH)  && $(MAKE) clean  
clean_couette :
	cd  $(CFD_SRC_PATH) && $(MAKE) clean
clean_cfd_coupler :
	rm -fr obj

#=======================================================================
# Compilation rules
#=======================================================================
$(O_CFD): obj/%.o : %.f90
	$(F90) $(FLAGS) -c $< -o $@
