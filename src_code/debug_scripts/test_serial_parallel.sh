#!/bin/sh
#
plat=intel

echo "Running build and run test for serial and parallel code. Check debug_log for full outputs"
echo "See Serial.out and parallel.out to check results - same restart file so should give same result"
echo "Difference is printed out at end of script for comparison"

#Clean previous tests
rm debug_log 
rm serial.out parallel.out temp_restart.out diff_serialparallel.out final_state > ./debug_log
cd ./../
make PLATFORM=$plat clean >> ./debug_scripts/debug_log

#Build serial code to generate initial configuration
make debug_s >> ./debug_scripts/debug_log
echo "Serial Build OK"

#Setup initial state
./md.exe -i ./debug_scripts/input_setup_restart.in > ./debug_scripts/temp_restart.out 
cat ./debug_scripts/temp_restart.out  >> ./debug_scripts/debug_log
cp ./results/final_state ./debug_scripts/
echo "Restart file generated OK"

#Run serial code to generate output
mv ./results/final_state ./debug_scripts/
mv ./debug_scripts/final_state ./debug_scripts/debug_initial
./md.exe -i ./debug_scripts/input_2048periodic.in -r ./debug_scripts/debug_initial > ./debug_scripts/serial.out
mv ./results/macroscopic_properties ./debug_scripts/serial.macro
echo "Serial run OK"
cat ./debug_scripts/serial.out  >> ./debug_scripts/debug_log


#Run parallel code to generate output
make clean >> ./debug_scripts/debug_log
make debug_p >> ./debug_scripts/debug_log
echo "Parallel Build OK"
mpiexec -n 8 ./md.exe -i ./debug_scripts/input_2048periodic.in -r ./debug_scripts/debug_initial > ./debug_scripts/parallel.out 
cat ./debug_scripts/parallel.out >> ./debug_scripts/debug_log
mv ./results/macroscopic_properties ./debug_scripts/parallel.macro
echo "Parallel run OK"
make clean

#Check serial and parallel outputs are the same
#awk -F "\"*;\"*" '{print $2}' ./debug_scripts/serial.out   > ./debug_scripts/s_$2
#awk -F "\"*;\"*" '{print $2}' ./debug_scripts/parallel.out > ./debug_scripts/p_$2
#diff -B --suppress-common-lines ./debug_scripts/serial.out ./debug_scripts/parallel.out >./debug_scripts/diff_serialparallel.out

./debug_scripts/test_rundiff.py
echo "Difference between serial and parallel macroscopic properties"
cat ./debug_scripts/diff.macro

echo "set title \"Serial - Parallel Total Energies\"
set xlabel \"iter\"
set ylabel \"Energy (LJU)\"
plot \"./debug_scripts/diff.macro\" u 1:7 w p lc 1 t \"TE\"" | gnuplot -persist
