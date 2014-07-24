import os
from cplfields import *
from postproc import PostProc, NoResultsInDir

class CPL_PostProc(PostProc):

    """ 
        Post processing class for CFD runs
    """

    def __init__(self,resultsdir,**kwargs):
        self.resultsdir = resultsdir
        self.plotlist = {} 

        # Check directory exists before instantiating object and check 
        # which files associated with plots are in directory
        if (not os.path.isdir(self.resultsdir)):
            print("Directory " +  self.resultsdir + " not found")
            raise IOError

        try:
            fobj = open(self.resultsdir + 'results/coupler_header','r')
        except IOError:
            raise NoResultsInDir

        possibles = {'CPL Velocity': CPL_vField,
                     'CPL Stress': CPL_PField}

        self.plotlist = {}
        for key, field in possibles.items(): 
            try:
                self.plotlist[key] = field(self.resultsdir)
            except AssertionError:
                pass 
