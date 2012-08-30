#!/bin/sh
#
if [ $# -lt 1 ]; then
	echo ""
	echo "Number of processors (#) for CFD required to run"
	echo "Parameter files under CFD files define this, details below:"
	echo "========================================================================="
	echo "CFD param.inc files:"
	grep "parameter (npx" ./DNS_main_code_Couette/param.inc
	echo "========================================================================="
	echo "Input should be of form: ./run.sh #"
 	exit 1
fi

#Make results directory if not present
mkdir -p ./results

#Clean it up if it is present
rm ./results/data ./results/ucvcwc.dble.* ./results/uuvvww.dble.* ./results/conold.dble.* ./results/pressure_ph.* ./results/archive* ./results/report ./results/SubDom_dble.*

#Check version of setup/grid/etc is current
echo " "
echo " ***************************** DNS_grid_generation_Couette *******************************"
grep "Lx"  ./DNS_grid_generation_Couette/input.file
grep "Ly"  ./DNS_grid_generation_Couette/input.file
grep "ngx"  ./DNS_grid_generation_Couette/input.file
grep "ngy"  ./DNS_grid_generation_Couette/input.file
echo " "
echo " ********************************** DNS_setup_Couette ************************************"
grep "xL"  ./DNS_setup_Couette/input.setup
grep "yL"  ./DNS_setup_Couette/input.setup
grep "nix"  ./DNS_setup_Couette/input.setup
grep "niy"  ./DNS_setup_Couette/input.setup
grep "niz"  ./DNS_setup_Couette/input.setup
echo " "
echo " ***********************Check difference in setup and grid folder *******************************"
echo " "
diff ./DNS_setup_Couette/grid.data ./DNS_grid_generation_Couette/grid.data
echo " "
echo " ***************************** DNS_main_code_Couette *******************************"
grep "er (ngx"  ./DNS_main_code_Couette/param.inc
echo " "
echo " ***********************************************************************************"
echo " "
echo " "
echo " ** SIMULATION STARTS * SIMULATION STARTS * SIMULATION STARTS * SIMULATION STARTS **"
echo " "
echo " "

#Move required files to directory
echo "  0" > data
mv ./data ./results/
cp ./DNS_main_code_Couette/input ./results/
cp ./DNS_main_code_Couette/parallel_couette.exe ./results/
cp ./DNS_setup_Couette/ucvcwc.dble.000000 ./results/
cp ./DNS_setup_Couette/uuvvww.dble.000000 ./results/
cp ./DNS_setup_Couette/conold.dble.000000 ./results/
cp ./DNS_setup_Couette/pressure_ph.000000 ./results/
cp ./DNS_setup_Couette/grid.data ./results/
cp ./DNS_setup_Couette/archive ./results/
cp ./DNS_setup_Couette/report ./results/
cp ./results/archive ./results/archive.000000

#Run Coupled code
#mpiexec -n 1  ./../../Couette_serial/continuum.exe : -n $1 ./../../MD_dCSE/src_code/md.exe
cd ./results
mpiexec -n $1 ./parallel_couette.exe
#mpiexec -n 1 xterm -geometry 150x20+1000+0 -hold -e ./results/parallel_couette.exe : -n $1 xterm -geometry 150x20+1000+220 -hold -e ./../../MD_dCSE/src_code/md.exe -i ./../../MD_dCSE/src_code/MD_coupled.in

