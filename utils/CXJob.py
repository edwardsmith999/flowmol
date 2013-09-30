#! /usr/bin/env python
import time
import os
import subprocess as sp

class CXJob:

    def __init__(self,
                 platform,
                 rundir,
                 jobname,
                 nproc,
                 walltime,
                 exec_cmd,
                 queue='pqtzaki',
                 icib='true'):

        # Store absolute run directory so we can move around
        # after submission
        absrundir = os.path.abspath(rundir)

        # Calculate number of nodes required
        if (nproc%8 != 0):
            quit('nprocs not a factor of 8. Aborting.')
        else:
            ncpus = 8
            select = nproc / ncpus
            
        # Create script
        if (platform == 'cx1'):

            script = ("#!/bin/bash \n"+
            "#PBS -N %s \n"+
            "#PBS -l walltime=%s \n"+
            "#PBS -l select=%s:ncpus=%s:icib=%s \n"+
            "#PBS -q %s \n"
            ) % (jobname, walltime, select, ncpus, icib, queue)

        elif (platform == 'cx2'):

            script = ("#!/bin/bash \n"+
            "#PBS -N %s \n"+
            "#PBS -l walltime=%s \n"+
            "#PBS -l select=%s:ncpus=%s:mpiprocs=8:ompthreads=1:mem=23500mb \n"
            ) % (jobname, walltime, select, ncpus)

        else:

            quit('Unrecognised platform in CXJob')

        script += 'module load intel-suite\n'
        script += 'module load mpi\n\n\n'
        script += 'cd ' + absrundir + '\n\n'
        script += 'date\n\n'
        script += exec_cmd + '\n\n'
        script += 'date\n'

        # Store variables
        self.script = script
        self.jobname = jobname
        self.rundir = rundir

        return

    def submit(self,blocking=False):

        # Open job submission file and write script to it
        job_file = self.jobname
        fobj = open(self.rundir+job_file,'w')
        fobj.write(self.script)
        fobj.close()

        # Submit script and store the submission ID
        sub_id = sp.Popen(["qsub", job_file], cwd=self.rundir,
                          stdout=sp.PIPE).communicate()[0]

        sub_id = sub_id.strip()

        # Alert user of job submission
        print('Submitted '+self.rundir+job_file+': '+sub_id)
      
        # If blocking, check qstat until sub_id no longer appears
        if (blocking):

            while True:

                qstat = sp.Popen(["qstat"],stdout=sp.PIPE).communicate()[0]
                if (sub_id not in qstat):
                    break 

                # Don't flood cx1 with qstat requests
                time.sleep(10)

        return

