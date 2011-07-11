 #!/bin/bash
#PBS -o /work/USERNAME/110405_MD_dCSE/results/mdresult.out
#PBS -S /bin/bash
#PBS -l select=1:ncpus=4
#PBS -l walltime=00:30:00 
#PBS -N MD
module load intel-suite
module load mpi
cd /work/USERNAME/110405_MD_dCSE/
date
mpiexec ./md.exe >> history
date

