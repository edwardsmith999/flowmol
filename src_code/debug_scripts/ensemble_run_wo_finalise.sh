#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Runs an ensemble of parallel_md.exe simulations"
	echo "Input of form ensemble_run.sh [no of runs]" #[no simultaneous procs]"
	echo "NOTE - Assumes code is setup correctly and random initial state is on"
    exit 1
fi

#debug flag so everything is printed to screen
#set -x

#Rebuild code so nbins is correct (MUST BE CHANGED IN ensemble_average.f90 TO AGREE WITH SIMULATIONS)
ifort -o ./results/ensemble_average.exe ./results/ensemble_average.f90
mkdir -p ./results/ensemble

#clean all previous ensembles
cd ./../results
./ensemble_average.exe clean_all
echo "clean previous ensemble averages"
cd ./../

#Assuming code is setup correctly and random initial state is on
#if [ $2 -gt 1 ]; then
#	echo "Runs ensemble on multiple processors simultaneously"
#	cd ./results/
#	for (( i=1; i<$2; i++ ))
#	do
#		dirname = ensemble_proc//i
#		mkdir dirname
#		cp ./../MD.in ./dirname
#		cp ./../parallel_md.exe ./dirname
#		mkdir dirname/results
#	done
#fi

for (( i=1; i<$1; i++ ))
do
	#run code
	mpiexec -n 1 ./parallel_md.exe 

	#ensemble average outputs
	cd ./results
	./ensemble_average.exe mbins
	./ensemble_average.exe mflux
	./ensemble_average.exe msnap
	./ensemble_average.exe vbins
	./ensemble_average.exe vflux
	./ensemble_average.exe psurface
	./ensemble_average.exe vsnap

	#output run number here
	echo "run number = ", $i

	#move back to code directory for next run
	cd ./../

done

#Run final case and store sum of all ensemble of runs 

#run code
mpiexec -n 1 ./parallel_md.exe

#ensemble average outputs
cd ./results

# Note finalise can be uncommented to divide sum of all ensembles by number of ensembles
./ensemble_average.exe mbins #finalise
./ensemble_average.exe mflux #finalise
./ensemble_average.exe msnap #finalise
./ensemble_average.exe vbins #finalise
./ensemble_average.exe vflux #finalise
./ensemble_average.exe psurface #finalise
./ensemble_average.exe vsnap #finalise

#output run number here
echo "Final number of ensemble averages = ", $i
echo "No averaging performed - divide results by number of ensembles to get average values"

#Copy header file to ensemble folder
cp simulation_header ./ensemble/
cp simulation_progress ./ensemble/
