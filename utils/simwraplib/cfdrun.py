import os
import shutil as sh
import subprocess as sp
import shlex

from simwraplib.platform import get_platform
from simwraplib.run import Run

class CFDRun(Run):

    def __init__(self,
                 srcdir='../CFD_dCSE/src_code/',
                 basedir=None,
                 rundir='../CFD_dCSE/src_code/main_code/',
                 executable='./parallel_couette.exe',
                 gridinputfile='input.file',
                 setupinputfile='input.setup',
                 inputfile='input',
                 outputfile='CFD.out',
                 dryrun=False
                ):

        self.srcdir = srcdir
        self.rundir = rundir
        self.basedir = basedir
        self.executable = executable
        self.setupinputfile = setupinputfile
        self.gridinputfile = gridinputfile
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.dryrun = dryrun

        if (self.basedir == None):
            quit('You must specify a base directory which contains'+
                 ' all input files.') 

        if (self.srcdir[-1] != '/'): self.srcdir += '/'
        if (self.rundir[-1] != '/'): self.rundir += '/'
        if (self.basedir[-1] != '/'): self.basedir += '/'
        if (self.executable[0:2] != './'): 
            self.executable = './' + self.executable

        # Work out what machine we're on
        self.platform = get_platform()

    def setup(self, existscheck=False): 
        
        self.create_rundir(existscheck=existscheck)       

        # Snapshot of the source directory
        cmd = 'tar -cPf ' + self.rundir + 'src.tar ' + self.srcdir
        sp.Popen(cmd,shell=True)

        # Store useful directory paths
        griddir = self.rundir + 'grid_generation/'
        setupdir = self.rundir + 'setup/'
        libdir = self.rundir + 'lib/'

        #########   # Perform setup  ###########
        # cp grid gen folder from source
        try:
            sh.copytree(self.srcdir+'grid_generation/',griddir)
        except OSError as e:
            pass

        # rm grid.data and gen_grid_data.exe
        for delfile in ['grid.data', 'Gen_grid.data.exe']:
            try:
                os.remove(griddir+delfile)
            except OSError as e:
                pass

        # put input.file into folder
        sh.copy(self.basedir+self.gridinputfile,griddir+self.gridinputfile)
        # compile & run grid gen
        for cmd in [
            'ifort -r8 -o Gen_grid.data.exe main.f90 mesh_tanh_stretch.f90',
            './Gen_grid.data.exe']:
            comp = sp.Popen(cmd, shell=True, cwd=griddir)
            comp.wait() 

        # cp setup and lib folders from source
        try:
            sh.copytree(self.srcdir+'setup/',setupdir)
        except OSError as e:
            pass
        try:
            sh.copytree(self.srcdir+'lib/',libdir)
        except OSError as e:
            pass
        # copy grid.data and setup input file to setup folder
        sh.copy(griddir+'grid.data',setupdir+'grid.data')
        sh.copy(self.basedir+self.setupinputfile,setupdir+self.setupinputfile)
        # compile and run setup
        for cmd in ['make simple', './a.out']:
            comp = sp.Popen(cmd, shell=True, cwd=setupdir)
            comp.wait()

        # Copy required run files into it from base
        basefiles = [ self.inputfile, self.executable ]
        for f in basefiles:
            sh.copy(self.basedir+f, self.rundir+f)
        setupfiles = [ 
                      'ucvcwc.dble.000000',
                      'uuvvww.dble.000000',
                      'conold.dble.000000',
                      'pressure_ph.000000',
                      'grid.data',
                      'archive',
                      'report',
                     ]
        for f in setupfiles:
            sh.copy(setupdir+f, self.rundir+f)
        sh.copy(self.rundir+'archive',self.rundir+'archive.000000')
        
        with open(self.rundir+'data','w') as f:
            f.write('0')
            f.close()
        
        self.prepare_inputs()

    def prepare_inputs(self, *args, **kwargs):
        print('CFD cannot currently change inputs')
    
    def get_nprocs(self, *args, **kwargs):  
       
        with open(self.srcdir+'main_code/param.inc') as f:
            for line in f:
                if ('npx' in line):
                    npx = int(line.split('npx=')[1].split(',')[0])
                    npy = int(line.split('npy=')[1].split(',')[0])
                    npz = int(line.split('npz=')[1].split(',')[0])
        nprocs = npx*npy*npz
        return nprocs 

    def execute(self, blocking=False, **kwargs):
        
        if (self.dryrun):
            print('DRYRUN -- no execution in ' + self.rundir)
            return
        
        if (self.platform == 'cx1' or
            self.platform == 'cx2' ):
            quit("Can't run on cx1 or cx2 with the CFDRun wrapper yet.")
        else:
            self.execute_local(blocking=blocking)
   
    def execute_local(self, blocking=False, **kwargs):

        # Store the number of processors required
        nprocs = self.get_nprocs()

        #Build runstring
        cmdstg = 'mpiexec -n ' + str(nprocs) + ' ' + self.executable 
        print(self.rundir + ':\t' + cmdstg)

        #Setup standard out and standard error files
        stdoutfile = self.rundir+self.outputfile
        stderrfile = self.rundir+self.outputfile+'_err'
        fstout = open(stdoutfile,'w')
        fsterr = open(stderrfile,'w')
        split_cmdstg = shlex.split(cmdstg)

        #Execute subprocess and create subprocess object
        self.proc = sp.Popen(split_cmdstg, cwd=self.rundir, stdin=None, 
                             stdout=fstout, stderr=fsterr)

        #If blocking, wait here
        if blocking:
            self.proc.wait()

        return

    def execute_cx(self, blocking=False):
        
        nprocs = self.get_nprocs()
        cmd = 'mpiexec ' + self.executable 
        cmd += ' > ' + self.outputfile
        cmd += ' 2> ' + self.outputfile + '_err'

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
 
    def finish(self, *args, **kwargs):
        print('No finish implementation yet for CFDRun')

