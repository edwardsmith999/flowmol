#!/bin/sh
#
if [ $# -lt 1 ]; then
	echo "Input of the form cgrep 'search string'"
else
	echo "==============COUPLER FILES====================="
	grep -n -i $1  *.f90
	echo "==============CFD FILES====================="
	grep -n -i $1  ./../../CFD_dCSE/src_code/main_code/*.f90
	echo "==============MD FILES====================="
	grep -n -i $1  ./../../MD_dCSE/src_code/*.f90
fi
