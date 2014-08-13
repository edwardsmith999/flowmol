import os
import shutil as sh
import subprocess as sp

import simwraplib.userconfirm as uc

class Run:
    
    def __init__(self, *args, **kwargs):
        raise NotImplementedError
    
    def setup(self, *args, **kwargs):
        raise NotImplementedError

    def prepare_inputs(self, *args, **kwargs):
        raise NotImplementedError

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
