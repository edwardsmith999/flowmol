#!/bin/sh
#
# Unpacks input folders generated by the pack_case command
#

if [ $# -lt 1 ]; then
	echo "Input of the form unpack_case.sh [Case Name]"
else
	echo "WARNING - this will replace all input files in the CFD and coupler directories"
	read -p "Are you sure? (y/n) "
	if [ "$REPLY" == "y" ]; then
		cp ./$1/input.file  ./../CFD_dCSE/src_code/grid_generation/
		cp ./$1/input.setup ./../CFD_dCSE/src_code/setup/
		cp ./$1/param.inc  ./../CFD_dCSE/src_code/main_code/
		cp ./$1/input  ./../CFD_dCSE/src_code/main_code/
		cp ./$1/COUPLER.in  ./../coupler_dCSE/src_code/
		cp ./$1/MD.in  ./../coupler_dCSE/src_code/md_data/
		echo "Case" $1 "unpacked in input folder"
	fi
fi
