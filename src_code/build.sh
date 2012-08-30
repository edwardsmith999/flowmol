#!/bin/sh
#
if [ $# -lt 1 ]; then
	cd ./../../;  make BUILD=debug couette_md;cd ./coupler_dCSE/src_code/;
else
	cd ./../../; make clean_all;  make BUILD=debug couette_md;cd ./coupler_dCSE/src_code/;
fi
