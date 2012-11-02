#!/bin/sh
#
BUILDtype=debug
if [ $# -lt 1 ]; then
	cd ./../../;  make BUILD=$BUILDtype couette_md ;cd ./coupler_dCSE/src_code/;
else
	cd ./../../; make clean_all;  make BUILD=$BUILDtype couette_md;cd ./coupler_dCSE/src_code/;
fi
