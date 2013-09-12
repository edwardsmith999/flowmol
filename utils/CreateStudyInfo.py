#! /usr/bin/env python
import numpy as np
from MDRun import MDRun
from MDStudy import MDStudy

# Number of threads and runs per thread
ncpus = 8
ncpusperrun = 2
nthreads = 15
nrunsperthread = 3
maxlicenses = ncpus / ncpusperrun

# Inputs that are the same for every thread
srcdir =  './../MD_dCSE/src_code/'
basedir = 'base'
executable = './creamsoda'
cylinderfile = None

#inputfiles = ['']*nrunsperthread
#outputfiles = ['']*nrunsperthread
#finish = [{}]*nrunsperthread
#initstates = [None]*nrunsperthread
#restartfiles = [None]*nrunsperthread

# Specify information needed for each run
inputfiles = ['ramp.in',
              'equil.in',
              'calcs.in']

outputfiles = ['ramp.out',
               'equil.out',
               'calcs.out']

finish = [{'final_state':'ramp.state'},
          {'final_state':'equil.state'},
          {}]

initstates = ['initial.state',
              None,
              None]

restartfiles = [None,
                'ramp.state',
                'equil.state']

runchanges = [{'NSTEPS':100},
              {'NSTEPS':200},
              {'NSTEPS':400}]

# Specify variables for each thread 
temperatures = np.linspace(0.5,2.0,num=nthreads)
temperatures = np.around(temperatures,decimals=2)

threadchanges = []
rundirs = []

for iThread in range(nthreads):

    T = str(temperatures[iThread])
    threadchanges.append({'INPUTTEMPERATURE':T})
    rundirs.append('TEMP_'+T)

#
#
#
#

# Build up lists of MDRun objects
threadlist = []
for iThread in range(nthreads):

    change_by_thread = threadchanges[iThread]
    rundir = rundirs[iThread] 

    runlist = []
    for iRun in range(nrunsperthread):

        change_by_run = runchanges[iRun]

        combined_changes = dict(change_by_thread.items() + 
                                change_by_run.items())

        run = MDRun(
                    srcdir,
                    basedir,
                    rundir,
                    executable,
                    inputfiles[iRun],
                    outputfiles[iRun],
                    inputchanges=combined_changes,
                    initstate=initstates[iRun],
                    restartfile=restartfiles[iRun],
                    cylinderfile=cylinderfile,
                    finishargs=finish[iRun]
                   )

        runlist.append(run)
    
    threadlist.append(runlist)

#
#
#
#
#
# Run the study

study = MDStudy(threadlist,maxlicenses)
