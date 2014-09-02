import os
import shutil as sh
import subprocess as sp

import simwraplib.userconfirm as uc

class Run:
    
    def __init__(self, *args, **kwargs):
        raise NotImplementedError
    
    def setup(self, *args, **kwargs):
        raise NotImplementedError

    def prepare_inputs(self, extrachanges=None, **kwargs):

        """
            Make alterations to the base input file (specified on 
            construction) that will be copied into the run directory.
        
            The input "changes" should be a dictionary of the form:
            
                changes = { 'DENSITY': 0.8, 'INPUTTEMPERATURE': 1.0 , ...}    
                
        """

        mod = self.inputmod(self.rundir+self.inputfile)

        #If additional changes, add these to the input changes
        if (extrachanges):
            self.inputchanges.update(extrachanges)

        for key in self.inputchanges:
            values = self.inputchanges[key]
            mod.replace_input(key,values)    
        
        return

    def get_nprocs(self, *args, **kwargs):  
        raise NotImplementedError

    def execute(self, *args, **kwargs):
        raise NotImplementedError 
    
    def finish(self, *args, **kwargs):
        raise NotImplementedError

    def create_rundir(self, existscheck=False, **kwargs):

        # Create run directory (and results dir inside it). If it already
        # exists, ask the user if they want to carry on.
        try:

            os.makedirs(self.rundir+'/results')
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
