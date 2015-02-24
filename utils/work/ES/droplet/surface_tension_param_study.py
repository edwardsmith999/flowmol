#! /usr/bin/env python2.7
import numpy as np
import matplotlib.pyplot as plt
import string
import itertools
import sys
import multiprocessing

sys.path.append('../../../')
import simwraplib as swl
import postproclib as ppl
from misclib import round_to_n
from plot_pressures import get_surface_tension

# Number of threads and runs per thread
ncpus = 12
maxlicenses = ncpus

# Inputs that are the same for every thread
srcdir =  '../../../../MD_dCSE/src_code/'
basedir = srcdir
executables = './parallel_md.exe'

# Specify information needed for each run
inputfile = 'MD_2phase_LJ.in'
outputfile = 'MD.out'
finish = [{'final_state':'final_state'}]

# Specify input file changes for each thread 
inputs1 = swl.InputDict({'INPUTTEMPERATURE': [round_to_n(i,3) for i in np.arange(1.0,2.5,0.1)]})
inputs2 = swl.InputDict({'LIQUIDDENSITY': [round_to_n(i,3) for i in np.arange(0.85,0.55,-0.02)]})
inputs3 = swl.InputDict({'RCUTOFF': [4.0]*20})
changes = inputs1 + inputs2 + inputs3
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
#study = swl.Study(threadlist,ncpus)

#Post process
#style = ['k','k--','b','r','g']; i = 0
#ax = plt.subplot(111)
#for thread in range(0,len(changes)):
#    rundir = srcdir + '../runs/' + filenames[thread] + '/results/'

#    #Get y MD profile
#    PPObj = ppl.MD_PostProc(rundir)
#    endrec = PPObj.plotlist['psurface'].maxrec-3; startrec = endrec-30

#    x, rho = PPObj.plotlist['rho'].profile(axis=0,startrec=startrec, endrec = endrec) 
#    x, T = PPObj.plotlist['T'].profile(axis=0,startrec=startrec, endrec = endrec)
#    T_estimate = np.mean(T[rho > 0.4])

#    if (thread in [2,4,5,6,7]):
#        label = 'T ' + str(round_to_n(T_estimate,3))+ " rcutoff " + round_to_n(float(PPObj.plotlist['psurface'].Raw.header.rcutoff),3)
#        ax.plot(x, rho[:,0], style[i], label=label  )
#        i += 1
#    
#    print(rundir, 'T,VA,CV', get_surface_tension(rundir, plotstuff=False))

#ax.legend()
#ax.set_ylim((0.0,1.0))
#plt.show()

