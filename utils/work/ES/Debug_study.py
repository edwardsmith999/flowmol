#! /usr/bin/env python2.7
import numpy as np
import string
import itertools
import sys
import multiprocessing

from MDRun import MDRun
from MDStudy import MDStudy

# Number of threads and runs per thread
ncpus = 8
nthreads = 12
nrunsperthread = 1
maxlicenses = ncpus

# Inputs that are the same for every thread
srcdir =  './../MD_dCSE/src_code/'
basedir = srcdir
executables = './parallel_md.exe'
initstate = './initial_state'


# Specify information needed for each run
inputfile = 'MD.in'
outputfile = 'MD.out'
finish = [{'final_state':'final_state'}]

# Specify input file changes for each thread 
inputs = {  'BIN2CELLRATIO':[[1.0,1.0,1.0],[0.5,0.5,0.5],[2.0,2.0,2.0]],
             'PROCESSORS':[[1,1,1],[2,2,2]]}
#             'MASS_OUTFLAG':[0,4],
#             'VELOCITY_OUTFLAG':[0,4],
#             'TEMPERATURE_OUTFLAG':[0,4],
#             'PRESSURE_OUTFLAG':[0,2],
#             'MFLUX_OUTFLAG':[0,1],
#             'VFLUX_OUTFLAG':[0,4]}

#inputs = { 'PROCESSORS':[[1,1,1],[2,1,1],[1,2,1],[1,1,2],[2,2,1],[1,2,2],[2,1,2],[2,2,2]]}

pythonscriptdir = ['./../MD_dCSE/src_code/post_proc/python/',
                   './../MD_dCSE/src_code/post_proc/python/work/ES/']
pyscriptname = 'compare_results'
args = ['thread,threadlist[0][0].rundir','thread,threadlist[thread][0].rundir']

#finishargs = {'python_script':[pythonscriptdir,pyscriptname,[args]]}

#Generate all permutations of inputs
product = [x for x in apply(itertools.product, inputs.values())]
changes = [dict(zip(inputs.keys(), p)) for p in product]
valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
threadlist =[]

for thread in range(0,len(changes)):
    filename = str(changes[thread].items())
    rundir = './../MD_dCSE/runs/'
    rundir += ''.join(c for c in filename if c in valid_chars[6:])
    
    run = MDRun(
                 srcdir,
                 basedir,
                 rundir,
                 executables,
                 inputfile,
                 outputfile,
                 initstate=initstate,
                 inputchanges=changes[thread],
                 finishargs = {},
                 dryrun=False
                 )

    runlist = [run]
    threadlist.append(runlist)
    print('Run in directory '  + rundir + ' and dryrun is '  + str(run.dryrun))

# Run the study
#Go to directory, import function and call
study = MDStudy(threadlist,maxlicenses)


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
        #print('Calling post processing function ' + value[1] + 
        #    'in directory '+ ' '.join(value[0])+'with arguments '+', '.join(value[2]))
        pp_obj = pp_mod.CompareResults()
        for thread in range(0,len(changes)):
            #print(thread,threadlist[0][0].rundir,thread,threadlist[thread][0].rundir)
            pp_obj.run(threadlist[0][0],threadlist[thread][0])
            #pp_obj.run(eval(value[2][0]),eval(value[2][1]))
