#!/bin/sh

# Command line args
nproc_md=$1
nproc_cfd=$2
md_extra_args="$3 $4" 

# Useful variables
MD_DIR=../../MD_dCSE/src_code
#MD_EXE=md.exe
MD_EXE=colonel_mustard
CFD_DIR=../../CFD_dCSE/src_code
#CFD_EXE=parallel_couette.exe
CFD_EXE=professor_plum

#Make results directory if not present
mkdir -p ./couette_data/results

# Clean up old data files
old_files= ./couette_data/data ./couette_data/ucvcwc.dble.* ./couette_data/uuvvww.dble.* ./couette_data/conold.dble.* ./couette_data/pressure_ph.* ./couette_md/pres_p* ./couette_data/archive* ./couette_data/report ./couette_data/SubDom_dble* ./couette_data/parallel_couette.exe
rm $old_files

# If sanity check succeeds and returns 0, run simulation
if ./checksims.py $nproc_md $nproc_cfd; then

	# Move required files to directory
	echo "  0" > data
	mv ./data ./couette_data/
	cp $CFD_DIR/main_code/input ./couette_data/
	cp $CFD_DIR/main_code/parallel_couette.exe ./couette_data/$CFD_EXE
	cp $CFD_DIR/setup/ucvcwc.dble.000000 ./couette_data/
	cp $CFD_DIR/setup/uuvvww.dble.000000 ./couette_data/
	cp $CFD_DIR/setup/conold.dble.000000 ./couette_data/
	cp $CFD_DIR/setup/pressure_ph.000000 ./couette_data/
	cp $CFD_DIR/setup/grid.data ./couette_data/
	cp $CFD_DIR/setup/archive ./couette_data/
	cp $CFD_DIR/setup/report ./couette_data/
	cp ./couette_data/archive ./couette_data/archive.000000

    cp $MD_DIR/parallel_md.exe ./md_data/$MD_EXE

	#Run Coupled code
	mpirun -n $nproc_md ./md_data/$MD_EXE -i ./md_data/MD.in $md_extra_args : -n $nproc_cfd ./couette_data/$CFD_EXE

fi
