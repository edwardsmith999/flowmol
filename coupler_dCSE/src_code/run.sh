#!/bin/sh
#
if [ $# -lt 1 ]; then
	echo ""
	echo "Number of processors (#) for BOTH MD and CFD required to run coupled code"
	echo "Parameter files under MD and CFD files define this, details below:"
	echo ""
	echo "========================================================================="
	echo "MD input file:"
	grep -A3 "PROCESSORS" ./md_data/MD_coupled.in
	echo "========================================================================="
	echo "CFD param.inc files:"
	grep "parameter (npx" ./../../CFD_dCSE/src_code/main_code/param.inc
	echo "========================================================================="
	echo "Input should be of form: ./coupled.exe #"
 	exit 1
fi

#Make results directory if not present
mkdir -p ./couette_data/results

#Clean it up if it is present
rm ./couette_data/data ./couette_data/ucvcwc.dble.* ./couette_data/uuvvww.dble.* ./couette_data/conold.dble.* ./couette_data/pressure_ph.* ./couette_md/pres_p* ./couette_data/archive* ./couette_data/report ./couette_data/SubDom_dble*

#Check version of setup/grid/etc is current
echo " "
echo " ***************************** grid_generation *******************************"
grep "Lx"  ./../../CFD_dCSE/src_code/grid_generation/input.file
grep "Ly"  ./../../CFD_dCSE/src_code/grid_generation/input.file
grep "ngx"  ./../../CFD_dCSE/src_code/grid_generation/input.file
grep "ngy"  ./../../CFD_dCSE/src_code/grid_generation/input.file
echo " "
echo " ********************************** setup ************************************"
grep "xL"  ./../../CFD_dCSE/src_code/setup/input.setup
grep "yL"  ./../../CFD_dCSE/src_code/setup/input.setup
grep "nix"  ./../../CFD_dCSE/src_code/setup/input.setup
grep "niy"  ./../../CFD_dCSE/src_code/setup/input.setup
grep "niz"  ./../../CFD_dCSE/src_code/setup/input.setup
echo " "
echo " ******** Checking difference in setup and grid folder - if not text below then ok *******"
echo " "
diff ./../../CFD_dCSE/src_code/setup/grid.data ./../../CFD_dCSE/src_code/grid_generation/grid.data
echo " "
echo " ***************************** main_code *******************************"
grep "er (ngx"  ./../../CFD_dCSE/src_code/main_code/param.inc
echo " "
echo " ***********************************************************************************"
echo " "
echo " "
echo " ** SIMULATION STARTS * SIMULATION STARTS * SIMULATION STARTS * SIMULATION STARTS **"
echo " "
echo " "

#Move required files to directory
echo "  0" > data
mv ./data ./couette_data/
cp ./../../CFD_dCSE/src_code/main_code/input ./couette_data/
cp ./../../CFD_dCSE/src_code/main_code/parallel_couette.exe ./couette_data/
cp ./../../CFD_dCSE/src_code/setup/ucvcwc.dble.000000 ./couette_data/
cp ./../../CFD_dCSE/src_code/setup/uuvvww.dble.000000 ./couette_data/
cp ./../../CFD_dCSE/src_code/setup/conold.dble.000000 ./couette_data/
cp ./../../CFD_dCSE/src_code/setup/pressure_ph.000000 ./couette_data/
cp ./../../CFD_dCSE/src_code/setup/grid.data ./couette_data/
cp ./../../CFD_dCSE/src_code/setup/archive ./couette_data/
cp ./../../CFD_dCSE/src_code/setup/report ./couette_data/
cp ./couette_data/archive ./couette_data/archive.000000

#Run Coupled code
mpirun -n $1 ./../../MD_dCSE/src_code/md.exe -i ./md_data/MD.in $2 : -n 1 ./couette_data/parallel_couette.exe 
#mpiexec -n $1 xterm -geometry 150x20+1000+0 -hold -e  ./../../MD_dCSE/src_code/md.exe -i ./md_data/MD_coupled.in $2 : -n 4 xterm -geometry 150x20+1000+220 -hold -e ./couette_data/parallel_couette.exe 
#mpiexec -n 1 ./couette_data/parallel_couette.exe  : -n $1 ./../../MD_dCSE/src_code/md.exe -i ./md_data/MD_coupled.in $2 # -r ./md_data/results/r0
#mpiexec -n 1 xterm -geometry 150x20+1000+0 -hold -e gdb ./couette_data/parallel_couette.exe :  -n $1 xterm -geometry 150x20+1000+220 -hold -e gdb ./../../MD_dCSE/src_code/md.exe -i /md_data/MD_coupled.in : -n 1 xterm -geometry 150x20+100+240 -hold -e gdb ./../../MD_dCSE/src_code/md.exe -i /md_data/MD_coupled.in

