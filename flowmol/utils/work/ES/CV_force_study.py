#! /usr/bin/env python2.7
import numpy as np
import string
import itertools
import sys
import multiprocessing

sys.path.insert(0, './../../')
from MDRun import MDRun
from MDStudy import MDStudy
from InputClass import InputList, InputDict

# Number of threads and runs per thread
ncpus = 21

# Inputs that are the same for every thread
srcdir =  './../MD_dCSE/src_code/'
basedir = srcdir
executables = './parallel_md.exe'

# Specify information needed for each run
inputfile = 'MD.in'
outputfile = 'MD.out'
finish = [{'final_state':'final_state'}]

# Specify input file changes for each thread 
inputs1 = InputDict({'INITIALNUNITS':[[5,5,5],[10,10,10],[20,20,20],[40,40,40]]})
inputs2 = InputDict({'BIN2CELLRATIO':[[0.5,0.5,0.5],[0.25,0.25,0.25],[0.125,0.125,0.125],[0.0625,0.0625,0.0625]]})
inputs3 = InputDict({'CV_FORCES':[['.true.',0],['.true.',3]]})

threadlist =[]
changes = (inputs1 + inputs2) * inputs3

print(changes)

for thread in range(0,len(changes)):
     filename = str(changes[thread].items())
     rundir = './../MD_dCSE/runs/' + changes.filenames()[thread]

     run = MDRun(
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
study = MDStudy(threadlist,ncpus)

finishargs = {'gnuplot_script':'./gnuplot_script.plt'}
