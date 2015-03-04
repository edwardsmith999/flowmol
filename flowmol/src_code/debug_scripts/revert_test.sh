#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Checks out a previous subversion revision to tmp directory, copies "
	echo "current input file and attempts to make debug_p and run. " 
	echo "Input of form revert_test.sh [revision No.]"
    exit 1
fi

DIR=$( pwd )
input_file="/MD.in"
echo $DIR$input_file
cd /tmp/
rm -rf ./coupled_code/
svn co -r $1  http://svn.ma.ic.ac.uk/subversion/edward/MDNS_repo/branch ./coupled_code/
cd ./coupled_code/MD_dCSE/src_code/
cp $DIR$input_file ./
make PLATFORM=intel debug_p
mpiexec -n 2 ./parallel_md.exe

