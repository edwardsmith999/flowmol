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
srcdir =  '../../../../../MD_dCSE/src_code/'
basedir = srcdir
executables = './parallel_md.exe'

# Specify information needed for each run
inputfile = 'MD_2Dinterface.in'
outputfile = 'MD.out'
finish = [{'final_state':'final_state'}]

# Specify input file changes for each thread 
inputs1 = swl.InputDict({'EIJ_WALL':[round_to_n(i,3) for i in np.arange(0.4,3.0,0.1)]})
inputs2 = swl.InputDict({'INITIALNUNITS':[]})
inputs3 = swl.InputDict({'REPEAT':[0,1,2]})
#Get array for domain sizes based on expected droplet behaviour
spread = [400, 30, 8]; hydrophobic = [200,60,8]
#At 40 degrees the angle is:
spreadthresh = 1.8
for i in inputs1['EIJ_WALL']:
    if float(i) > spreadthresh:
        inputs2['INITIALNUNITS'].append(spread)
    elif float(i) <= spreadthresh:
        inputs2['INITIALNUNITS'].append(hydrophobic)

#inputs2 = swl.InputDict({'LIQUIDDENSITY':[0.6]})
#inputs2 = swl.InputDict({'LIQUIDDENSITY':[round_to_n(i,3) for i in np.arange(0.3,0.9,0.3)]})
#inputs3 = swl.InputDict({'GASDENSITY':[0.001,0.01,0.1]})


changes = (inputs1 + inputs2) * inputs3
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

