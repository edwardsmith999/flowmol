#!/bin/sh
set -e

#Check paramters agree
checkparam() {
	# Checks if any params.
	if [ -z $1 ]; then
		echo "No parameters passed to function - checkparam."
		exit 0
	else
		if [ $1 -ne $2 ]; then
			echo ""
			echo "Grid generation file and setup gridsize do not agree"
			echo "Grid generation = $1"
			echo "Setup = $2"
			echo ""
	 		exit 0
		fi
	fi
} 

# Change the grid based on the input parameters in the grid_generation input file

#Check input files agree

#Domain sizes
Lx=$(awk '/Lx/  {print $1}' ./DNS_grid_generation_Couette/input.file)
xL=$(awk '/xL/  {print $1}' ./DNS_setup_Couette/input.setup)
px=$(awk '/ngx/ {print $3}' ./DNS_main_code_Couette/param.inc)

#checkparam $Lx $xL

awk '/Ly/  { Ly = $1}' ./DNS_grid_generation_Couette/input.file
awk '/yL/  { yL = $1}' ./DNS_setup_Couette/input.setup
awk '/ngy/  {print $1}' ./DNS_main_code_Couette/param.inc

#checkparam $Ly $yL

#awk '/ngx/ { ngx = $1}' ./DNS_grid_generation_Couette/input.file
#awk '/ngy/ { ngy = $1}' ./DNS_grid_generation_Couette/input.file

cd ./DNS_grid_generation_Couette
rm ./Gen_grid.data.exe
ifort        -r8 -o Gen_grid.data.exe  main.f90 mesh_tanh_stretch.f90
./Gen_grid.data.exe
if [ $# -eq 1 ]; then
	rm 	./Gen_sym_grid.exe
	ifort        -r8 -o Gen_sym_grid.exe main_sym_grid.f90
	./Gen_sym_grid.exe
	mv grid.data.2 ./grid.data
fi
cp grid.data ./../DNS_setup_Couette/
cd ./../DNS_setup_Couette
./a.out
cd ./../
