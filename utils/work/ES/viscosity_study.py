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
ncpus = 16
maxlicenses = ncpus

# Inputs that are the same for every thread
srcdir =  '../../../MD_dCSE/src_code/'
basedir = srcdir
executables = './parallel_md.exe'
initstate = './initial_state'


# Specify information needed for each run
inputfile = 'MD_visc.in'
outputfile = 'MD.out'
finish = [{'final_state':'final_state'}]

# Specify input file changes for each thread 
#inputs1 = InputDict({'TEMPERATURE':[np.arange(0.3,1.6,0.1).tolist()]})
inputs1 = InputDict({'INPUTTEMPERATURE':[0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,2.0]})
inputs2 = InputDict({'LIQUIDDENSITY':[0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]})
inputs3 = InputDict({'WALLSLIDEV':[1,2]})

pythonscriptdir = ['../../../MD_dCSE/src_code/post_proc/python/',
                   '../../../MD_dCSE/src_code/post_proc/python/work/ES/']

pyscriptname = 'Get_Visc'
args = ['thread,threadlist[0][0].rundir','thread,threadlist[thread][0].rundir']

#Generate all permutations of inputs
threadlist =[]
changes = inputs1 * inputs2 * inputs3

for thread in range(0,len(changes)):
    filename = str(changes[thread].items())
    rundir = '../../../MD_dCSE/runs/VISCSTUDY_' + changes.filenames()[thread]
    
    run = MDRun(
                 srcdir,
                 basedir,
                 rundir,
                 executables,
                 inputfile,
                 outputfile,
                 inputchanges=changes[thread],
                 finishargs = {},
                 dryrun=True
                 )

    runlist = [run]
    threadlist.append(runlist)
    print('Run in directory '  + rundir + ' and dryrun is '  + str(run.dryrun))

# Run the study
#Go to directory, import function and call
#study = MDStudy(threadlist,maxlicenses)


#Post processing for study
studyfinishargs = {'python_script':[pythonscriptdir,pyscriptname,args]}

for key,value in studyfinishargs.iteritems():
    if key == 'python_script':
        for dirs in value[0]:
            sys.path.append(dirs)
        try:
            pp_mod = __import__(value[1])
        except:
            raise
        pp_obj = pp_mod.Get_Visc()
        for thread in range(0,len(changes)):
            pp_obj.run(threadlist[thread][0])
