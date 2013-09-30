#! /usr/bin/env python
import os
import errno
import shlex
import subprocess as sp
import shutil as sh
import string

import userconfirm as uc
from Platform import get_platform
from InputUtils import InputMod
from CXJob import CXJob
from GnuplotUtils import GnuplotUtils

class MDRun:
        
    """ 
        When instatiated, MDRun will create a new folder as a "run" 
        directory (if it doesn't already exist) and copy the necessary files 
        from a "base" directory into it. All the information necessary for
        execution of the run is stored from the following inputs to the
        constructor: 

            Directories:

                srcdir  - path to source code folder
                basedir - path from which input/restart/etc files are copied
                rundir  - path in which the run will take place, files are 
                          copied from basedir to here.

            Copied files:

                executable   - name of executable (e.g. parallel_md.exe)
                inputfile    - input file name
                extrafiles   - a List of additional files copied to rundir
                initstate    - initial state file that is TO BE COPIED 
                               FROM THE BASE DIRECTORY
                restartfile  - initial state "restart" file, assumed to be
                               already located at the given path that is
                               RELATIVE TO THE RUN DIRECTORY
                cylinderfile - cylinder restart info file
   
            New files: 

                outputfile - output file name

            Other:

                inputchanges - dictionary of changes to make to the inputfile
                               once copied from the base directory
                finishargs   - dictionary of keywords and associated objects
                               that specify a range of actions to perform
                               when the execution of the run has finished. See
                               the comments at the top of finish() for more
                               info. 
                


        Example usage from a higher level:

            run = MDRun('../MD_dCSE/src_code/', etc.)
            run.setup_directory()
            run.execute()
            run.finish()

    """

    def __init__(self, 
                 srcdir,
                 basedir,
                 rundir,
                 executable,
                 inputfile,
                 outputfile,
                 inputchanges={},
                 initstate=None,
                 restartfile=None,
                 cylinderfile=None,
                 jobname='default_jobname',
                 walltime='24:00:00',
                 icib='true',
                 queue='pqtzaki',
                 extrafiles=None, 
                 finishargs={},
                 dryrun=False):

        self.srcdir = srcdir
        self.basedir = basedir
        self.rundir = rundir  
        self.executable = executable
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.inputchanges = inputchanges
        self.initstate = initstate
        self.restartfile = restartfile
        self.cylinderfile = cylinderfile 
        self.extrafiles = extrafiles
        self.finishargs = finishargs
        self.dryrun = dryrun

        # Add slashes to end of folders if they aren't already there
        if (self.srcdir[-1] != '/'): self.srcdir += '/'
        if (self.rundir[-1] != '/'): self.rundir += '/'
        if (self.basedir[-1] != '/'): self.basedir += '/'

        # Check initstate and restartfile are not both specified
        if (initstate != None and restartfile != None):
            quit('Error: both initial state and restart files are not None')

        # Keep a list of files to iterate over later
        self.copyfiles = [executable, inputfile]
        if (initstate): self.copyfiles.append(initstate)
        if (cylinderfile): self.copyfiles.append(cylinderfile)
        if (extrafiles): 
            for f in extrafiles:
                self.copyfiles.append(f)

        # Work out what machine we're on
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

    def change_inputs(self,extrachanges=None):

        """
            Make alterations to the base input file (specified on 
            construction) that will be copied into the run directory.
        
            The input "changes" should be a dictionary of the form:
            
                changes = { 'DENSITY': 0.8, 'INPUTTEMPERATURE': 1.0 , ...}    
                
        """

        mod = InputMod(self.rundir+self.inputfile)

        #If additional changes, add these to the input changes
        if (extrachanges):
            self.inputchanges.update(extrachanges)

        for key in self.inputchanges:
            values = self.inputchanges[key]
            mod.replace_input(key,values)    
        
        return

    def setup_directory(self, existscheck=True):

        """
            Create a new run directory and copy the relevant files into
            it. 

        """

        # Create run directory (and results dir inside it). If it already
        # exists, ask the user if they want to carry on.
        try:

            os.makedirs(self.rundir+'/results')

            # Copy post_proc folder from base directory (if not there already)
            sh.copytree(self.basedir+'post_proc/',self.rundir+'post_proc/')
            # Make a snapshot of the source code and store in a tarball
            cmd = 'tar -cPf ' + self.rundir+'src.tar ' + self.srcdir+'*.f90'
            sp.Popen(cmd,shell=True)

            print('Created '+self.rundir)

        except OSError as e:

            if existscheck:

                if (e.errno == errno.EEXIST):
                    message = ('Directory ' + self.rundir + ' already exists.\n' +
                               'Continue anyway? (files could be overwritten)')
                    print(message)
                    go = uc.confirm(prompt=message,resp='y')
                    if (not go):
                        quit('Stopping.')
                else:
                    quit('Error creating directory.')
            else:

                pass

        # Copy files and save new locations to instance variables
        for f in self.copyfiles:

            # Do nothing if the files are the same
            if (self.basedir+f == self.rundir+f):
                pass
            else:
                sh.copy(self.basedir+f, self.rundir+f)

        # Make changes to the input file once it has been copied
        self.change_inputs()

        return

    def get_nprocs(self):

        with open(self.rundir+self.inputfile,'r') as f:

            for line in f:
                if ('PROCESSORS' in line):
                    npx = int(f.next()) 
                    npy = int(f.next()) 
                    npz = int(f.next()) 
                    break
        
        return npx*npy*npz

    def execute(self, blocking=False, nprocs=0):

        """
            Wrapper for execute_cx1, execute_local, and eventually others.

        """

        if (self.dryrun):
            print('DRYRUN -- no execution in ' + self.rundir)
            return
        
        if (self.platform == 'cx1' or
            self.platform == 'cx2' ):
            self.execute_cx(blocking=blocking)
        else:
            self.execute_local(blocking=blocking, nprocs=nprocs)

    def execute_cx(self, blocking=False):

        """
            Runs an executable from the directory specified  
            during instatiation of the object. 

        """ 

        cpuspernode = 8
        nprocs = self.get_nprocs()
        
        cmd = 'mpiexec ' + self.executable 
        cmd += ' -i ' + self.inputfile 

        if self.initstate != None:
            cmd += ' -r ' + self.initstate
        elif self.restartfile != None:
            cmd += ' -r ' + self.restartfile 

        if self.cylinderfile != None:
            cmd += ' -c ' + self.cylinderfile 

        cmd += ' > ' + self.outputfile

        # Create CXJob object based on the platform
        if (self.platform == 'cx1'):

            job = CXJob(self.platform, self.rundir, self.jobname, nprocs, 
                        self.walltime, cmd, queue=self.queue, icib=self.icib)

        elif (self.platform == 'cx2'):

            job = CXJob(self.platform, self.rundir, self.jobname, nprocs, 
                        self.walltime, cmd) 

        else:

            quit('Unrecognised platform in execute_cx')

        # Submit the job
        job.submit(blocking=blocking)
    
        
    def execute_local(self, blocking=False, nprocs=0):

        """
            Runs an executable from the directory specified  
            during instatiation of the object. 

        """ 

        # Store the number of processors required
        if (nprocs == 0):
            nprocs = self.get_nprocs()

        #Build runstring
        cmdstg = 'mpiexec -n ' + str(nprocs) + ' ' + self.executable 
        cmdstg += ' -i ' + self.inputfile 

        if self.initstate != None:
            cmdstg += ' -r ' + self.initstate
        elif self.restartfile != None:
            cmdstg += ' -r ' + self.restartfile 

        if self.cylinderfile != None:
            cmdstg += ' -c ' + self.cylinderfile 

        print(self.rundir + ':\t' + cmdstg)

        #Setup standard out and standard error files
        fstout = open(self.rundir+self.outputfile,'w')
        fsterr = open(self.rundir+self.outputfile+'_err','w')
        split_cmdstg = shlex.split(cmdstg)

        #Execute subprocess and create subprocess object
        self.proc = sp.Popen(split_cmdstg, cwd=self.rundir, stdin=None, 
                             stdout=fstout, stderr=fsterr)

        #If blocking, wait here
        if blocking:
            self.proc.wait()

        return

    def finish(self):
        
        """
            Perform a selection of actions once the simulation has finished.
   
            self.finishargs must be a dictionary of the form:
            
                {'keyword': object}
 
            Keyword options:
        
                final_state - when not None, move the results/final_state
                              file to a specified string (object) location.

        """

        print('Simulation in directory ' + self.rundir + 'has finished')

        for key,value in self.finishargs.iteritems():

            if key == 'final_state':
                
                src = self.rundir + 'results/final_state'
                dst = self.rundir + value
                print('Moving ' + src + ' to ' + dst)
                sh.move(src,dst)

            if key == 'gnuplot_script':

                outfile = 'tempout'
                gnufile = GnuplotUtils(value)
                gnufile.specify_outputfile(rundir=self.rundir,
                                           outfilename=outfile)
   
                #Run gnuplot and generate outputs
                cmdstr = ' gnuplot ' + value
                gnurun = sp.Popen(shlex.split(cmdstr),cwd=self.rundir)
                print('Running gnuplot script ' + value + ' with output ' + outfile)
                gnurun.wait()

                #Copy results back to calling directory and name based on rundir
                sh.copy(self.rundir+outfile, outfile)
                valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
                newname = ''.join(c for c in outfile+self.rundir if c in valid_chars[6:])
                os.rename(outfile, newname)

        return

