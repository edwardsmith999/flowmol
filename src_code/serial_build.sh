#!/bin/sh
#
BUILDtype=debug
if [ $# -lt 1 ]; then
	cd ./../../;  make BUILD=$BUILDtype serial_couette md  ; cd ./coupler_dCSE/src_code/;
else
	cd ./../../; make clean_all;  make BUILD=$BUILDtype serial_couette md ; cd ./coupler_dCSE/src_code/;
fi
