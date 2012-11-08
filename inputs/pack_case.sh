#!/bin/sh
#
# Packs input folders and adds to subversion server to share with other users
#
if [ $# -lt 1 ]; then
	echo "Input of the form pack_case.sh [Case Name]"
else
	mkdir $1
	cp ./../CFD_dCSE/src_code/grid_generation/input.file ./$1
	cp ./../CFD_dCSE/src_code/setup/input.setup ./$1
	cp ./../CFD_dCSE/src_code/setup/uy_input ./$1
	cp ./../CFD_dCSE/src_code/main_code/param.inc ./$1
	cp ./../CFD_dCSE/src_code/main_code/input ./$1
	cp ./../coupler_dCSE/src_code/COUPLER.in ./$1
	cp ./../coupler_dCSE/src_code/md_data/MD.in ./$1
	svn add $1
	echo "Case" $1 "added to input folder - use upack_case.sh [" $1 "] in ./input folder to use" > ./$1/log 
	svn ci $1 -F ./$1/log
fi
