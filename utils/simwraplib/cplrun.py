import shlex
import os
import shutil as sh
import subprocess as sp

from simwraplib.run import Run
from simwraplib.platform import get_platform

class CPLRun(Run):
    
    def __init__(self, 
                 mdrun, 
                 cfdrun, 
                 srcdir='../coupler_dCSE/src_code/',
                 basedir=None,
                 rundir='../coupler_dCSE/src_code/',
                 inputfile='COUPLER.in',
                 outputfile='COUPLER.out',
                 dryrun=False
                ):

        self.mdrun = mdrun
        self.cfdrun = cfdrun
        self.srcdir = srcdir
        self.basedir = basedir
        self.rundir = rundir
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.dryrun = dryrun

        if (self.basedir == None):
            quit('You must specify a base directory which contains'+
                 ' all input files.') 

        if (self.srcdir[-1] != '/'): self.srcdir += '/'
        if (self.rundir[-1] != '/'): self.rundir += '/'
        if (self.basedir[-1] != '/'): self.basedir += '/'

        self.platform = get_platform()


    def setup(self):

        # Create dir and copy coupler input file
        self.create_rundir()
        sh.copy(self.basedir+self.inputfile,self.rundir+self.inputfile)

        # Enforce directory structure for now
        self.mdrunsubdir = 'md_data/'
        self.cfdrunsubdir = 'couette_data/'
        print('Resetting MDRun rundir to ' + self.rundir + self.mdrunsubdir)
        print('Resetting CFDRun rundir to ' + self.rundir + self.mdrunsubdir)
        self.mdrun.rundir = self.rundir + self.mdrunsubdir 
        self.cfdrun.rundir = self.rundir + self.cfdrunsubdir 

        self.mdrun.setup()
        self.cfdrun.setup()

    def prepare_inputs(self):
        print("CPLRun can't change inputs yet")
       
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
            quit("Can't run on cx1 or cx2 with the CFDRun wrapper yet.")
        else:
            self.execute_local(blocking=blocking)

    def execute_local(self, blocking=False):

        # MD part of command string
        cmd = 'mpirun -n {0:d} ./md_data/{1:s} -i ./md_data/{2:s}'.format(
              self.mdprocs, self.mdrun.executable, self.mdrun.inputfile)

        if (self.mdrun.initstate != None):
            cmd += ' -r ' + self.mdrunsubdir + self.mdrun.initstate
        elif (self.mdrun.restartfile != None):
            cmd += ' -r ' + self.mdrun.restartfile 

        # CFD part of command string
        cmd += ' : -n {0:d} ./couette_data/{1:s}'.format(
               self.cfdprocs, self.cfdrun.executable)

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
