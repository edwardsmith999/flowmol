#! /usr/bin/env python
import os
import errno
import tarfile
import subprocess as sp
import shutil as sh
import shlex

import userconfirm as uc
from InputUtils import InputMod

svn_post_proc = ('http://svn.ma.ic.ac.uk/subversion/edward/MDNS_repo/' +
                 'branch/MD_dCSE/src_code/post_proc/')

class MDRun:
        
    """ 
        When instatiated, MDRun will create a new folder as a run 
        directory (if it doesn't already exist) and copy the necessary files 
        into it. 

        Example usage from a higher level:

            run = MDRun('../MD_dCSE/src_code/',...)
            run.change_inputs(changes)
            run.setup_directory()
            run.execute()

    """
    

    def __init__(self, 
                 srcdir,
                 basedir,
                 rundir,
                 executable,
                 inputfile,
                 outputfile, 
                 restartfile=None,
                 cylinderfile=None):

        self.srcdir = srcdir    
        self.basedir = basedir    
        self.rundir = rundir  
        self.executable = executable
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.restartfile = restartfile
        self.cylinderfile = cylinderfile 

        # Add slashes to end of folders if they aren't already there
        if (self.srcdir[-1] != '/'): self.srcdir += '/'
        if (self.rundir[-1] != '/'): self.rundir += '/'
        if (self.basedir[-1] != '/'): self.basedir += '/'

        # Keep a list of files to iterate over later
        self.files = [executable, inputfile]
        if (restartfile): self.files.append(restartfile)
        if (cylinderfile): self.files.append(cylinderfile)

    def change_inputs(self,changes):

        """
            Make alterations to the base input file (specified on 
            construction) that will be copied into the run directory.
        
            The input "changes" should be a dictionary of the form:
            
                changes = { 'DENSITY': 0.8, 'INPUTTEMPERATURE': 1.0 , ...}    
                
        """

        mod = InputMod(self.rundir+self.inputfile)
        
        for key in changes:
            values = changes[key]
            mod.replace_input(key,values)    

    def setup_directory(self,existscheck=True):

        """
            Create a new run directory and copy the relevant files into
            it. 

        """

        # Create run directory (and results dir inside it). If it already
        # exists, ask the user if they want to carry on.
        try:

            os.makedirs(self.rundir+'/results')

            # Check out post_proc folder (if not there already)
            cmd = ('svn co ' + svn_post_proc + ' ' + self.rundir +'post_proc')
            log = sp.check_output(cmd.split())    
            # Make a snapshot of the source code and store in a tarball
            cmd = 'tar -cPf ' +  self.rundir+'src.tar '  +  self.srcdir+'*.f90'
            sp.Popen(cmd,shell=True)

            print('Created '+self.rundir)

        except OSError as e:

            if existscheck:
                if (e.errno == errno.EEXIST):
                    message = ('Directory ' + self.rundir + ' already exists.\n' +
                               'Continue anyway? (files could be overwritten)')
                    go = uc.confirm(prompt=message,resp='y')
                    if (not go):
                        quit('Stopping.')
                else:
                    quit('Error creating directory.')
            else:
                pass

        # Copy files and save new locations to instance variables
        for f in self.files:

            # Do nothing if the files are the same
            if (self.basedir+f == self.rundir+f):
                pass
            else:
                sh.copy(self.basedir+f, self.rundir+f)


    def execute(self,blocking=False):

        """
            Runs an executable from the directory specified  
            during instatiation of the object. 

        """    

        #Build runstring
        cmdstg = self.executable + ' -i ' + self.inputfile 
        if self.restartfile != None:
            cmdstg += ' -r ' + self.restartfile 
        if self.cylinderfile != None:
            cmdstg += ' -c ' + self.cylinderfile 

        #Setup standard out and standard error files
        fstout = open(self.rundir+self.outputfile,'w')
        fsterr = open(self.rundir+self.outputfile+'_err','w')
        split_cmdstg = shlex.split(cmdstg)

        #Execute subprocess and create subprocess object
        self.proc = sp.Popen(split_cmdstg, cwd=self.rundir, stdin=None, stdout=fstout, stderr=fsterr)

        #If blocking, wait here
        if blocking:
            self.proc.wait()
