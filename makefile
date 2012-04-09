# Top makefile for the coupler
#
# The organisation schema is as follow:
#
# three top sectors: CFD MD COUPLER
#
# CFD and MD sectors contain sets of codes for cfd and md computations
# COUPLER is a library

#
#
# CFD sector
#
# Couette solver
#
# serial
clean_serial_couette serial_couette serial_couette_solo : CFD_SRC_PATH := ./Couette_serial 
serial_couette serial_couette_solo : CFD_TARGET := continuum.exe
#
#parallel
#
couette clean_couette couette_solo clean_couette_solo : CFD_SRC_PATH := ./CFD_dCSE/src_code/DNS_main_code_Couette 
couette clean_couette couette_solo clean_couette_solo : MAKEFILE_NAME := -f makefile.planes_fftz_fftx 
couette clean_couette couette_solo clean_couette_solo : CFD_TARGET := parallel_couette.exe 

#
# MD sector 
#
# so far we have only onne MD code, no need for target specific values as for CFD (see above)
#
MD_SRC_PATH := ./MD_dCSE/src_code
MD_TARGET  := p
#
# Coupler sector
#
COUPLER_SRC_PATH := coupler_dCSE/src_code

#
# variables for submakes
#
export BUILD := opt


# Tasks

.PHONY. :serial_couette_md serial_couette md coupler serial_couette_solo md_solo clean_all  clean_couette_md clean_couette clean_md clean_coupler clean_couette_solo clean_md_solo  

all :
	@echo "Top level coupled build - Type help for options"

couette_md : couette md 
serial_couette_md :serial_couette md 

md : coupler
	cd $(MD_SRC_PATH)  && $(MAKE) USE_COUPLER=yes $(MD_TARGET)
couette : coupler
	cd $(CFD_SRC_PATH) && $(MAKE) $(MAKEFILE_NAME) USE_COUPLER=yes $(CFD_TARGET)
serial_couette :  coupler
	cd $(CFD_SRC_PATH) && $(MAKE) USE_COUPLER=yes $(CFD_TARGET)
coupler : 
	cd $(COUPLER_SRC_PATH) && $(MAKE) 
md_solo :
	cd $(MD_SRC_PATH)  && $(MAKE) USE_COUPLER=no $(MD_TARGET)
serial_couette_solo : 
	cd $(CFD_SRC_PATH) && $(MAKE) USE_COUPLER=no $(CFD_TARGET)
couette_solo :
	cd $(CFD_SRC_PATH) && $(MAKE) $(MAKEFILE_NAME) USE_COUPLER=no $(CFD_TARGET)


clean_all : clean_coupler clean_cfd clean_md 

clean_couette_md : clean_couette_solo clean_md

clean_serial_couette_md : clean_serial_couette clean_md 

clean_md :
	cd  $(MD_SRC_PATH)  && $(MAKE) clean  
clean_cfd :
	cd  $(CFD_SRC_PATH) && $(MAKE) $(MAKEFILE_NAME) clean
clean_coupler :
	cd $(COUPLER_SRC_PATH) && $(MAKE) clean
clean_couette clean_couette_solo : clean_cfd
clean_serial_couette_solo clean_serial_couette : clean_cfd
clean_md_solo : clean_md

help:
	@echo "======================================================================================"
	@echo "TOP LEVEL COUPLED BUILD - individual MD/CFD can be built using make in MD or CFD "
	@echo "directories directly as well as from this folder. Tasks are as follows:"
	@echo ""
	@echo "couette_md - build simple diffusive couette flow solver coupled with MD code"
	@echo "cfd - Builds  DNS CFD code only but liked witj coupler library "
	@echo "md - Builds MD code only  but liked witj coupler library"
	@echo "cfd_solo - Builds standalone  DNS CFD code "
	@echo "md_solo  - Builds standalone MD executable "
	@echo "GENERAL Options"
	@echo "clean_all - Deletes all .mod, .obj, coupler library and other temporary files from everywhere"
	@echo "clean_couette_md - Deletes all .mod, .obj and other temporary files from couette & MD"
	@echo "clean_md - Deletes all .mod, .obj and other temporary files from md"
	@echo "clean_couette - Deletes all .mod, .obj and other temporary files from couette"
	@echo "clean_md_solo - The same as clean_md"
	@echo "clean_couette_solo - The same as clean_couette"
	@echo "======================================================================================"




