#! /usr/bin/env python
import numpy as np
from MDRun import MDRun
from MDStudy import MDStudy

# Inputs that are the same for every thread
srcdir =  '../branch/MD_dCSE/src_code/'
basedir = '../base'
executable = './professor_plum'
cylinderfile = 'cyl_0.6' 

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

initstates = ['fill.state_228',
              None,
              None]

restartfiles = [None,
                'ramp.state',
                'equil.state']

walltimes = ['04:00:00',
             '04:00:00',
             '30:00:00']

# Specify variables for each thread 
# speeds = np.arange(0.04,0.11,0.005)
speeds = [0.06]

# Build up lists of MDRun objects
threadlist = []
for speed in speeds:

    rundir = '../study/'+str(speed) 

    runlist = []
    for iRun in range(3):

        jobname = str(speed) + '_' + inputfiles[iRun].replace('.in','.job')

        if ( inputfiles[iRun] == 'ramp.in' ):

            changes = {'ROTATE_CYLINDERS':[0.0,speed,200]}

        else:
            
            changes = {'ROTATE_CYLINDERS':[speed,speed,0.0]}

        run = MDRun(
                    srcdir,
                    basedir,
                    rundir,
                    executable,
                    inputfiles[iRun],
                    outputfiles[iRun],
                    inputchanges=changes,
                    initstate=initstates[iRun],
                    restartfile=restartfiles[iRun],
                    cylinderfile=cylinderfile,
                    finishargs=finish[iRun],
                    walltime=walltimes[iRun],
                    jobname=jobname
                   )

        runlist.append(run)
    
    threadlist.append(runlist)

#
#
#
#
#
# Run the study
dummy = 0
study = MDStudy(threadlist,dummy)
