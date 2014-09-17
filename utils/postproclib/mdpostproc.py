import os
import numpy as np
import sys
import math as maths
import glob
#import collections

from mdfields import *
from headerdata import *
from postproc import PostProc
from pplexceptions import NoResultsInDir

class MD_PostProc(PostProc):

    """ 
        Post processing class for MD runs
    """

    def __init__(self,resultsdir,**kwargs):
        self.resultsdir = resultsdir
        self.plotlist = {} #collections.OrderedDict
        self.error = {}
        self.name = self.resultsdir.split('/')[-2]

        # Check directory exists before instantiating object and check 
        # which files associated with plots are in directory
        self.potentialfiles = ( "mslice", "mbins", "msnap","vslice", "vbins", 
                                "vsnap","pvirial", "pVA", "pVA_k","pVA_c", 
                                "visc", "mflux","vflux", "pplane", "psurface",
                                "esnap", "eflux", "eplane","esurface", "Fvext", 
                                "viscometrics", "rdf", "rdf3d", "ssf", "Fext",
                                "Tbins", "vPDF", "msolv", "mpoly", "vsolv",
                                "vpoly")        

        if (not os.path.isdir(self.resultsdir)):
            print("Directory " +  self.resultsdir + " not found")
            raise IOError
            
        self.fields_present = []
        for fname in self.potentialfiles:
            if (glob.glob(self.resultsdir+fname)):
                self.fields_present.append(fname)
            if (glob.glob(self.resultsdir+fname+'.*')):
                self.fields_present.append(fname.strip().split('.')[0])

        self.fieldfiles1 = list(set(self.fields_present) & set(self.potentialfiles)) 

        try:
            Header1 = MDHeaderData(self.resultsdir)
        except IOError:
            raise NoResultsInDir

        #Mass
        if 'mbins' in (self.fieldfiles1):
            m1 = MD_mField(self.resultsdir, **kwargs)
            self.plotlist.update({'mbins':m1})
            d1 = MD_dField(self.resultsdir, **kwargs)
            self.plotlist.update({'rho':d1})

        if 'msolv' in (self.fieldfiles1):
            m1 = MD_mField(self.resultsdir, fname='msolv', **kwargs)
            self.plotlist.update({'msolv':m1})
            d1 = MD_dField(self.resultsdir, fname='msolv', **kwargs)
            self.plotlist.update({'rho_solv':d1})

        if 'mpoly' in (self.fieldfiles1):
            m1 = MD_mField(self.resultsdir, fname='mpoly', **kwargs)
            self.plotlist.update({'mpoly':m1})
            d1 = MD_dField(self.resultsdir, fname='mpoly', **kwargs)
            self.plotlist.update({'rho_poly':d1})

        #Momentum
        if 'vbins' in (self.fieldfiles1):
            M1 = MD_pField(self.resultsdir, **kwargs)
            self.plotlist.update({'vbins':M1})
            M1 = MD_momField(self.resultsdir, **kwargs)
            self.plotlist.update({'rho u':M1})

        if 'vsolv' in (self.fieldfiles1):
            M1 = MD_pField(self.resultsdir,fname='vsolv',**kwargs)
            self.plotlist.update({'vsolv':M1})
            M1 = MD_momField(self.resultsdir,fname='vsolv',**kwargs)
            self.plotlist.update({'rho u_solv':M1})

        if 'vpoly' in (self.fieldfiles1):
            M1 = MD_pField(self.resultsdir,fname='vpoly',**kwargs)
            self.plotlist.update({'vpoly':M1})
            M1 = MD_momField(self.resultsdir,fname='vpoly',**kwargs)
            self.plotlist.update({'rho u_poly':M1})

        #Kinetic energy
        if 'Tbins' in (self.fieldfiles1):
            KE1 = MD_EField(self.resultsdir, **kwargs)
            self.plotlist.update({'Tbins':KE1})

        #Mass snapshots
        if 'msnap' in (self.fieldfiles1):
            m1 = MD_mField(self.resultsdir,fname='msnap', **kwargs)
            self.plotlist.update({'msnap':m1})
            m1 = MD_dField(self.resultsdir,fname='msnap', **kwargs)
            self.plotlist.update({'rho_snap':m1})

        #Velocity snapshots
        if 'vsnap' in (self.fieldfiles1):
            v1 = MD_pField(self.resultsdir,fname='vsnap', **kwargs)
            self.plotlist.update({'vsnap':v1})
            v1 = MD_momField(self.resultsdir,fname='vsnap', **kwargs)
            self.plotlist.update({'rhou_snap':v1})

        #VA stress
        if 'pVA' in (self.fieldfiles1):
            P1 = MD_pVAField(self.resultsdir,fname='pVA', **kwargs)
            self.plotlist.update({'pVA':P1})
        if 'pVA_k' in (self.fieldfiles1):
            P1 = MD_pVAField(self.resultsdir,fname='pVA_k', **kwargs)
            self.plotlist.update({'pVA_k':P1})
        if 'pVA_c' in (self.fieldfiles1):
            P1 = MD_pVAField(self.resultsdir,fname='pVA_c', **kwargs)
            self.plotlist.update({'pVA_c':P1})

        #CV fluxes
        if 'mflux' in (self.fieldfiles1):
            flux1 = MD_mfluxField(self.resultsdir,'mflux', **kwargs)
            self.plotlist.update({'mflux':flux1})

        if 'vflux' in (self.fieldfiles1):
            flux1 = MD_pCVField(self.resultsdir,'vflux', **kwargs)
            self.plotlist.update({'vflux':flux1})

        #External forces
        if 'Fext' in (self.fieldfiles1):
            Fext1 = MD_FField(self.resultsdir,'Fext', **kwargs)
            self.plotlist.update({'Fext':Fext1})

        #CV energy snapshot
        if 'esnap' in (self.fieldfiles1):
            esnap1 = MD_EField(self.resultsdir,'esnap', **kwargs)
            self.plotlist.update({'esnap':esnap1})


        #CV stresses
        if 'psurface' in (self.fieldfiles1):
            stress1 = MD_pCVField(self.resultsdir,'psurface', **kwargs)
            self.plotlist.update({'psurface':stress1})

        #Stress Heating
        if ('vflux'   in (self.fieldfiles1) and
           'psurface' in (self.fieldfiles1)):
            stress1 = MD_CVStressheat_Field(self.resultsdir, **kwargs)
            self.plotlist.update({'stress_heating':stress1})

        #CV energy fluxes
        if 'eflux' in (self.fieldfiles1):
            eflux1 = MD_efluxField(self.resultsdir,'eflux', **kwargs)
            self.plotlist.update({'eflux':eflux1})

        #CV surface power
        if 'esurface' in (self.fieldfiles1):
            energy1 = MD_efluxField(self.resultsdir,'esurface', **kwargs)
            self.plotlist.update({'esurface':energy1})

        #CV Energy due to external body 
        if 'Fvext' in (self.fieldfiles1):
            Fvext1 = MD_EField(self.resultsdir,'Fvext', **kwargs)
            self.plotlist.update({'Fvext':Fvext1})

        #Velocity
        if ('mbins' in (self.fieldfiles1) and 'vbins' in (self.fieldfiles1)):
            try:
                v1 = MD_vField(self.resultsdir, **kwargs)
                self.plotlist.update({'u':v1})
            except DataMismatch:
                pass

            try:
                v1 = MD_rhouuField(self.resultsdir, **kwargs)
                self.plotlist.update({'rhouu':v1})
            except DataMismatch:
                pass

        #strain
        if ('mbins' in (self.fieldfiles1) and 'vbins' in (self.fieldfiles1)):
            try:
                v1 = MD_strainField(self.resultsdir, **kwargs)
                self.plotlist.update({'Strain':v1})
            except DataMismatch:
                pass


        #Vorticity
        if ('mbins' in (self.fieldfiles1) and 'vbins' in (self.fieldfiles1)):
            try:
                v1 = MD_vortField(self.resultsdir, **kwargs)
                self.plotlist.update({'Vorticity':v1})
            except DataMismatch:
                pass

        #Dissipation
        if ('mbins' in (self.fieldfiles1) and 'vbins' in (self.fieldfiles1)):
            try:
                v1 = MD_dissipField(self.resultsdir, **kwargs)
                self.plotlist.update({'Dissipation':v1})
            except DataMismatch:
                pass

        #Velocity snapshot
        if ('msnap' in (self.fieldfiles1) and 'vsnap' in (self.fieldfiles1)):
            try:
                v1 = MD_vField(self.resultsdir,rectype='snap', **kwargs)
                self.plotlist.update({'u_snap':v1})
            except DataMismatch:
                pass
        #Temperature
        if ('mbins' in (self.fieldfiles1) and 
            'vbins' in (self.fieldfiles1) and 
            'Tbins' in (self.fieldfiles1)):
            try:
                T1 = MD_TField(self.resultsdir, **kwargs)
                self.plotlist.update({'T':T1})
            except DataMismatch:
                pass

            try:
                T1 = MD_dTdrField(self.resultsdir, **kwargs)
                self.plotlist.update({'dTdr':T1})
            except DataMismatch:
                pass

        if (len(self.plotlist) == 0):
            raise NoResultsInDir

