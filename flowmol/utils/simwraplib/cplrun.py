import shlex
import os
import shutil as sh
import subprocess as sp

from simwraplib.run import Run
from simwraplib.inpututils import InputMod
from simwraplib.platform import get_platform
from simwraplib.hpc import CXJob

class CPLRun(Run):
    
    def __init__(self, 
                 mdrun, 
                 cfdrun, 
                 srcdir='../coupler_dCSE/src_code/',
                 basedir=None,
                 rundir='../coupler_dCSE/src_code/',
                 inputfile='COUPLER.in',
                 outputfile='COUPLER.out',
                 inputchanges={},
                 jobname='default_jobname',
                 walltime='24:00:00',
                 icib='true',
                 queue='pqtzaki',
                 dryrun=False
                ):

        self.mdrun = mdrun
        self.cfdrun = cfdrun
        self.srcdir = srcdir
        self.basedir = basedir
        self.rundir = rundir
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.inputchanges = inputchanges
        self.dryrun = dryrun

        if (self.basedir == None):
            if (self.mdrun.basedir == self.cfdrun.basedir):
                self.basedir = self.mdrun.basedir
            else:
                quit('Unable to obtain base directory for coupler inputs')

        if (self.srcdir[-1] != '/'): self.srcdir += '/'
        if (self.rundir[-1] != '/'): self.rundir += '/'
        if (self.basedir[-1] != '/'): self.basedir += '/'

        self.platform = get_platform()
        # Store more values if cx1
        if (self.platform == 'cx1'):
            self.jobname = jobname
            self.walltime = walltime 
            self.icib = icib 
            self.queue = queue 
        elif (self.platform == 'cx2'):
            self.jobname = jobname
            self.walltime = walltime

        # Set input modifier to be normal kind
        self.inputmod = InputMod

    def setup(self):

        # Create dir and copy coupler input file
        self.create_rundir()
        sh.copy(self.basedir+self.inputfile,self.rundir+self.inputfile)

        # Enforce directory structure for now
        self.mdrunsubdir = './md_data/'
        self.cfdrunsubdir = './couette_data/'
        print('Resetting MDRun rundir to ' + self.rundir + self.mdrunsubdir)
        print('Resetting CFDRun rundir to ' + self.rundir + self.mdrunsubdir)
        self.mdrun.rundir = self.rundir + self.mdrunsubdir 
        self.cfdrun.rundir = self.rundir + self.cfdrunsubdir 

        self.prepare_inputs()
        self.mdrun.setup()
        self.cfdrun.setup()
       
    def get_nprocs(self):
        self.mdprocs = self.mdrun.get_nprocs() 
        self.cfdprocs = self.cfdrun.get_nprocs() 
        return self.mdprocs + self.cfdprocs

    def execute(self, blocking=False, **kwargs):
        
        if (self.dryrun):
            print('DRYRUN -- no execution in ' + self.rundir)
            return
        
        if (self.platform == 'cx1' or
            self.platform == 'cx2' ):
            self.execute_cx(blocking=blocking)
        else:
            self.execute_local(blocking=blocking)

    def execute_cx(self, blocking=False):

        
        # MD part of command string
        cmd =  'mpiexec heterostart '
        cmd += '{0:d} {1:s}{2:s} '.format(
              self.cfdprocs, self.cfdrunsubdir, self.cfdrun.executable)
        cmd += '{0:d} {1:s}{2:s} -i {1:s}{3:s} '.format(
              self.mdprocs, self.mdrunsubdir, self.mdrun.executable, 
              self.mdrun.inputfile)

        if (self.mdrun.initstate != None):
            cmd += ' -r ' + self.mdrunsubdir + self.mdrun.initstate
        elif (self.mdrun.restartfile != None):
            cmd += ' -r ' + self.mdrun.restartfile 
       
        nprocs = self.get_nprocs()
 
        # Create CXJob object based on the platform
        if (self.platform == 'cx1'):

            job = CXJob(self.rundir, self.jobname, nprocs, self.walltime, 
                        cmd, queue=self.queue, icib=self.icib)

        elif (self.platform == 'cx2'):

            job = CXJob(self.rundir, self.jobname, nprocs, self.walltime, cmd) 

        else:

            quit('Unrecognised platform in execute_cx')

        # Submit the job
        job.submit(blocking=blocking)

    def execute_local(self, blocking=False):

        # MD part of command string
        cmd = 'mpirun -n {0:d} {1:s}{2:s} -i {1:s}{3:s}'.format(
              self.mdprocs, self.mdrunsubdir, self.mdrun.executable, 
              self.mdrun.inputfile)

        if (self.mdrun.initstate != None):
            cmd += ' -r ' + self.mdrunsubdir + self.mdrun.initstate
        elif (self.mdrun.restartfile != None):
            cmd += ' -r ' + self.mdrun.restartfile 

        # CFD part of command string
        cmd += ' : -n {0:d} {1:s}{2:s} '.format(
               self.cfdprocs,self.cfdrunsubdir,self.cfdrun.executable)

        stdoutfile = self.rundir+self.outputfile
        stderrfile = self.rundir+self.outputfile+'_err'
        fstout = open(stdoutfile,'w')
        fsterr = open(stderrfile,'w')
        split_cmdstg = shlex.split(cmd)

        # Run the executable
        self.proc = sp.Popen(split_cmdstg, cwd=self.rundir, stdin=None, 
                             stdout=fstout, stderr=fsterr)

        #If blocking, wait here
        if blocking:
            self.proc.wait()

        return
        
    def finish(self):
        self.mdrun.finish()
        self.cfdrun.finish()
