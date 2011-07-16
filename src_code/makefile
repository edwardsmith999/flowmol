#
#  CFD sector
#
CFD_SRC_PATH := ../../CFD_dCSE/src_code/main_code_MG_ddType_dP
#
# MD sector
#
MD_SRC_PATH := ../../MD_dCSE/src_code

#
# tasks
#

all : cfd md 

cfd : 
	cd $(CFD_SRC_PATH) && make all "VERSION=debug" "USE_COUPLER=yes" 
md  :
	cd $(MD_SRC_PATH)  && make debug_p "USE_COUPLER=yes"

.PHONY. : clean
clean :
	cd  $(CFD_SRC_PATH) && make clean && cd - && cd $(MD_SRC_PATH)  && make clean  