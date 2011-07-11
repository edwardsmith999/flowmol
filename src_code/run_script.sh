#!/bin/bash
#PBS -o $WORK/my_first_dns/data/myOutput.out
#PBS -S /bin/bash
#PBS -l select=1:ncpus=4
#PBS -l walltime=01:00:00 
#PBS -N my_first_dns
module load intel-suite
module load fftw/2.1.5-double
module load mpi
cd $WORK/my_first_dns/data/
date
mpiexec ./a.out >> run_history
date

