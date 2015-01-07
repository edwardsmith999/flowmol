#! /usr/bin/env python2.7
import numpy as np
import string
import itertools
import sys
import multiprocessing

sys.path.append('../../../')
import simwraplib as swl
from misclib import round_to_n

# Number of threads and runs per thread
ncpus = 12
maxlicenses = ncpus

# Inputs that are the same for every thread
srcdir =  '../../../../MD_dCSE/src_code/'
basedir = srcdir
executables = './parallel_md.exe'

# Specify information needed for each run
inputfile = 'MD_2Dinterface.in'
outputfile = 'MD.out'
finish = [{'final_state':'final_state'}]

# Specify input file changes for each thread 
inputs1 = swl.InputDict({'WALLSLIDEV':[round_to_n(i,3) for i in np.arange(0.0,1.0,0.05)]})
changes = inputs1.expand()
filenames = changes.filenames()

threadlist =[]
for thread in range(0,len(changes)):
     rundir = srcdir + '../runs/' + filenames[thread]

     run = swl.MDRun(
                  srcdir,
                  basedir,
                  rundir,
                  executables,
                  inputfile,
                  outputfile,
                  inputchanges=changes[thread],
                  finishargs = {},
                  dryrun=False
                  )

     runlist = [run]
     threadlist.append(runlist)
     print('Run in directory '  + rundir + ' and dryrun is '  + str(run.dryrun))

# Run the study
study = swl.Study(threadlist,ncpus)

