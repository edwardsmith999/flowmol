#!/bin/sh
#
if [ $# -lt 1 ]; then
	echo ""
	echo "Number of processors (#) for BOTH MD and CFD required to run coupled code"
	echo "Parameter files under MD and CFD files define this, details below:"
	echo ""
	echo "========================================================================="
	echo "MD input file:"
	grep -A3 "PROCESSORS" ./../../MD_dCSE/src_code/MD.in
	echo "========================================================================="
	echo "CFD param.inc files:"
	grep "parameter (npx" ./../../CFD_dCSE/src_code/DNS_main_code_Couette/param.inc
	echo "========================================================================="
	echo "Input should be of form: ./coupled.exe #"
 	exit 1
fi

#Make results directory if not present
mkdir -p ./couette_data/results

#Clean it up if it is present
rm ./couette_data/data ./couette_data/ucvcwc.dble.* ./couette_data/uuvvww.dble.* ./couette_data/conold.dble.* ./couette_data/pressure_ph.* ./couette_data/archive* ./couette_data/report ./couette_data/SubDom_dble*

#Check version of setup/grid/etc is current
echo " "
echo " ***************************** DNS_grid_generation_Couette *******************************"
grep "Lx"  ./../../CFD_dCSE/src_code/DNS_grid_generation_Couette/input.file
grep "Ly"  ./../../CFD_dCSE/src_code/DNS_grid_generation_Couette/input.file
grep "ngx"  ./../../CFD_dCSE/src_code/DNS_grid_generation_Couette/input.file
grep "ngy"  ./../../CFD_dCSE/src_code/DNS_grid_generation_Couette/input.file
echo " "
echo " ********************************** DNS_setup_Couette ************************************"
grep "xL"  ./../../CFD_dCSE/src_code/DNS_setup_Couette/input.setup
grep "yL"  ./../../CFD_dCSE/src_code/DNS_setup_Couette/input.setup
grep "nix"  ./../../CFD_dCSE/src_code/DNS_setup_Couette/input.setup
grep "niy"  ./../../CFD_dCSE/src_code/DNS_setup_Couette/input.setup
grep "niz"  ./../../CFD_dCSE/src_code/DNS_setup_Couette/input.setup
echo " "
echo " ******** Checking difference in setup and grid folder - if not text below then ok *******"
echo " "
diff ./../../CFD_dCSE/src_code/DNS_setup_Couette/grid.data ./../../CFD_dCSE/src_code/DNS_grid_generation_Couette/grid.data
echo " "
echo " ***************************** DNS_main_code_Couette *******************************"
grep "er (ngx"  ./../../CFD_dCSE/src_code/DNS_main_code_Couette/param.inc
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
cp ./../../CFD_dCSE/src_code/DNS_main_code_Couette/input ./couette_data/
cp ./../../CFD_dCSE/src_code/DNS_main_code_Couette/parallel_couette.exe ./couette_data/
cp ./../../CFD_dCSE/src_code/DNS_setup_Couette/ucvcwc.dble.000000 ./couette_data/
cp ./../../CFD_dCSE/src_code/DNS_setup_Couette/uuvvww.dble.000000 ./couette_data/
cp ./../../CFD_dCSE/src_code/DNS_setup_Couette/conold.dble.000000 ./couette_data/
cp ./../../CFD_dCSE/src_code/DNS_setup_Couette/pressure_ph.000000 ./couette_data/
cp ./../../CFD_dCSE/src_code/DNS_setup_Couette/grid.data ./couette_data/
cp ./../../CFD_dCSE/src_code/DNS_setup_Couette/archive ./couette_data/
cp ./../../CFD_dCSE/src_code/DNS_setup_Couette/report ./couette_data/
cp ./couette_data/archive ./couette_data/archive.000000

#Run Coupled code
#mpiexec -n 1  ./../../Couette_serial/continuum.exe : -n $1 ./../../MD_dCSE/src_code/md.exe 
mpiexec -n 1 ./couette_data/parallel_couette.exe  : -n $1 ./../../MD_dCSE/src_code/md.exe -i ./md_data/MD_coupled.in 
#mpiexec -n 1 xterm -geometry 150x20+1000+0 -hold -e ./couette_data/parallel_couette.exe : -n $1 xterm -geometry 150x20+1000+220 -hold -e ./../../MD_dCSE/src_code/md.exe -i ./../../MD_dCSE/src_code/MD_coupled.in

