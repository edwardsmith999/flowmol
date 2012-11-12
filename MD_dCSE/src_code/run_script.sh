#! /bin/bash
#PBS -o /work/es205/codes/coupled/coupler_dCSE/src_code/mdresult.out
#PBS -S /bin/bash
#PBS -l select=2:ncpus=8:icib=true:mem=16000mb
#PBS -l walltime=24:00:00 
#PBS -N COUETTE_TEST
#PBS -q pqtzaki
module load intel-suite
module load fftw/2.1.5-double
module load mpi
cd /work/es205/codes/coupled/coupler_dCSE/src_code
##############################################################
# 		Prepare the simulations
##############################################################
#Make results directory if not present
mkdir -p ./couette_data/results
#Clean it up if it is present
rm ./couette_data/data ./couette_data/ucvcwc.dble.* ./couette_data/uuvvww.dble.* ./couette_data/conold.dble.* ./couette_data/pressure_ph.* ./couette_data/pres_p* ./couette_data/archive* ./couette_data/report ./couette_data/SubDom_dble*
#Check version of setup/grid/etc is current
echo " "
echo " ***************************** DNS_grid_generation_Couette *******************************"
grep "Lx"  ./../../CFD_dCSE/src_code/grid_generation/input.file
grep "Ly"  ./../../CFD_dCSE/src_code/grid_generation/input.file
grep "ngx"  ./../../CFD_dCSE/src_code/grid_generation/input.file
grep "ngy"  ./../../CFD_dCSE/src_code/grid_generation/input.file
echo " "
echo " ********************************** DNS_setup_Couette ************************************"
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
echo " ***************************** DNS_main_code_Couette *******************************"
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
cp ./couette_data/report ./couette_data/report.000000
date
mpiexec heterostart 8 ./couette_data/parallel_couette.exe 8 ./../../MD_dCSE/src_code/md.exe -i ./md_data/MD.in   >>testy
date
