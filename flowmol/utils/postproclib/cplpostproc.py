import os
from cplfields import *
from postproc import PostProc
from pplexceptions import NoResultsInDir, DataNotAvailable
from mdpostproc import MD_PostProc
from cfdpostproc import CFD_PostProc
from serial_cfdpostproc import Serial_CFD_PostProc

class CPL_PostProc(PostProc):

    """ 
        Post processing class for Coupled runs
    """

    def __init__(self,resultsdir,**kwargs):
        self.resultsdir = resultsdir
        # Check directory exists before instantiating object and check 
        # which files associated with plots are in directory
        if (not os.path.isdir(self.resultsdir)):
            print("Directory " +  self.resultsdir + " not found")
            raise IOError

        self.plotlist = {}
        try:
            fobj = open(self.resultsdir + 'results/coupler_header','r')
        except IOError:
            raise NoResultsInDir

        #print(self.md_pp)
        #print(self.cfd_pp)

        #For MD and Transflow DNS code
#        if ('Transflow_u' in self.md_field.plotlist.keys() and
#            'CFD Velocity' in self.cfd_field.plotlist.keys()):
#            u = CPL_vField(self.resultsdir)
#            self.plotlist.update({'CPL Velocity':u})

#        if (any([('pVA','pVA_k','pVA_c') in self.md_field.plotlist.keys()]) and
#            'CPL Stress' in self.cfd_field.plotlist.keys()):
#            u = CPL_stressField(self.resultsdir)
#            self.plotlist.update({'CPL Velocity':u})

#        #For MD and Serial Couette solver
#        if ('u' in self.md_field.plotlist.keys() and
#            'u' in self.cfd_field.plotlist.keys()):
#            u = CPL_Serial_CFD_vField(self.resultsdir)
#            self.plotlist.update({'CPL Velocity':u})

#        if ('rho u' in self.md_field.plotlist.keys() and
#            'u' in self.cfd_field.plotlist.keys()):
#            u = CPL_Serial_CFD_momField(self.resultsdir)
#            self.plotlist.update({'CPL Momentum':u})

#        if ('psurface' in self.md_field.plotlist.keys() and
#            'CFD surface Tau_xx' in self.cfd_field.plotlist.keys()):
#            psurface = CPL_Serial_CFD_stressField(self.resultsdir)
#            self.plotlist.update({'CPL surface':psurface})

        possibles = {'CPL Velocity': [CPL_vField, CPL_Serial_CFD_vField],
                     'CPL Momentum': [CPL_Serial_CFD_momField],
                     'CPL Stress': [CPL_stressField,CPL_Serial_CFD_stressField]}

        for key, fields in possibles.items():
            for field in fields:
                print(field)
                try:
                    self.plotlist[key] = field(self.resultsdir)
                except AssertionError:
                    pass
                except DataNotAvailable:
                    pass
        
