#! /usr/bin/env python
import numpy as np
from MDRun import MDRun
from MDStudy import MDStudy

# Number of threads and runs per thread
ncpus = 8
nproc = 2
nthreads = 15
nrunsperthread = 3
maxlicenses = ncpus / nproc

# Inputs that are the same for every thread
srcdir =  './../MD_dCSE/src_code/'
basedir = 'base'
executable = './creamsoda'
cylinderfile = None

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

# Specify file changes for the run
nsteps = [100,200,3000]
# Specify input file changes for each thread 
temperatures = np.linspace(0.5,2.0,num=nthreads)

# Build up the lists of changes

threadlist = []

for iThread in range(nthreads):

    change_by_thread = {}
    change_by_thread['INPUTTEMPERATURE'] = str(temperatures[iThread])

    rundir = 'run' + str(iThread)

    runlist = []
    for iRun in range(nrunsperthread):

        change_by_run = {}
        change_by_run['NSTEPS'] = str(nsteps[iRun]) 
        combined_changes = dict(change_by_thread.items() + change_by_run.items())

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


# Run the study
study = MDStudy(threadlist,maxlicenses)
