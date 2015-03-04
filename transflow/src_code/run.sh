#!/bin/sh
#
if [ $# -lt 1 ]; then
	echo ""
	echo "Number of processors (#) for CFD required to run"
	echo "Parameter files under CFD files define this, details below:"
	echo "========================================================================="
	echo "CFD param.inc files:"
	grep "parameter (npx" ./main_code/param.inc
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
echo " ***************************** grid_generation *******************************"
grep "Lx"  ./grid_generation/input.file
grep "Ly"  ./grid_generation/input.file
grep "ngx"  ./grid_generation/input.file
grep "ngy"  ./grid_generation/input.file
echo " "
echo " ********************************** setup ************************************"
grep "xL"  ./setup/input.setup
grep "yL"  ./setup/input.setup
grep "nix"  ./setup/input.setup
grep "niy"  ./setup/input.setup
grep "niz"  ./setup/input.setup
echo " "
echo " ***********************Check difference in setup and grid folder *******************************"
echo " "
diff ./setup/grid.data ./grid_generation/grid.data
echo " "
echo " ***************************** main_code *******************************"
grep "er (ngx"  ./main_code/param.inc
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
cp ./main_code/input ./results/
cp ./main_code/parallel_couette.exe ./results/
cp ./setup/ucvcwc.dble.000000 ./results/
cp ./setup/uuvvww.dble.000000 ./results/
cp ./setup/conold.dble.000000 ./results/
cp ./setup/pressure_ph.000000 ./results/
cp ./setup/grid.data ./results/
cp ./setup/archive ./results/
cp ./setup/report ./results/
cp ./results/archive ./results/archive.000000

#Run Coupled code
#mpiexec -n 1  ./../../Couette_serial/continuum.exe : -n $1 ./../../MD_dCSE/src_code/md.exe
cd ./results
mpiexec -n $1 ./parallel_couette.exe
#mpiexec -n 1 xterm -geometry 150x20+1000+0 -hold -e ./results/parallel.exe : -n $1 xterm -geometry 150x20+1000+220 -hold -e ./../../MD_dCSE/src_code/md.exe -i ./../../MD_dCSE/src_code/MD_coupled.in

