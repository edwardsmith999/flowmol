#
# targets for CFD and MD
#
MD_TARGET  := p
CFD_TARGET := opt
#
#  CFD sector
#
CFD_SRC_PATH := ../../CFD_dCSE/src_code/main_code_MG_ddType_dP
#
# MD sector
#
MD_SRC_PATH := ../../MD_dCSE/src_code

.PHONY. : clean_all clean_cfd clean_md cfd md

all : cfd md 

cfd : 
	cd $(CFD_SRC_PATH) && $(MAKE) USE_COUPLER=yes $(CFD_TARGET)
md  :
	cd $(MD_SRC_PATH)  && $(MAKE) USE_COUPLER=yes $(MD_TARGET)

clean_all : clean_cfd clean_md

clean_cfd :
	cd  $(CFD_SRC_PATH) && $(MAKE) clean 

clean_md :
	cd $(MD_SRC_PATH)  && $(MAKE) clean  

