#!/bin/bash
#PBS -o /work/tjelly/RESEARCH/DNS_COUETTE
#PBS -S /bin/bash
#PBS -l select=1:ncpus=4:icib=true
#PBS -l walltime=48:00:00 
#PBS -N Couette
#PBS -q pqtzaki
module load intel-suite
module load fftw/2.1.5-double
module load mpi
cd /work/tjelly/RESEARCH/DNS_COUETTE
date
mpiexec ./a.out >> run_history
date


