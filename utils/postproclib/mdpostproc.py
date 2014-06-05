import os
import numpy as np
import sys
import math as maths
import glob
#import collections

from mdfields import *
from headerdata import *

class NoResultsInDir(Exception):
    pass

class MD_PostProc:

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
                                "Tbins", "vPDF" )        

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
            self.plotlist.update({'\\rho':d1})

        #Momentum
        if 'vbins' in (self.fieldfiles1):
            M1 = MD_pField(self.resultsdir, **kwargs)
            self.plotlist.update({'vbins':M1})
            M1 = MD_momField(self.resultsdir, **kwargs)
            self.plotlist.update({'\\rho u':M1})

        #Kinetic energy
        if 'Tbins' in (self.fieldfiles1):
            KE1 = MD_EField(self.resultsdir, **kwargs)
            self.plotlist.update({'Tbins':KE1})

        #Mass snapshots
        if 'msnap' in (self.fieldfiles1):
            m1 = MD_mField(self.resultsdir,fname='msnap', **kwargs)
            self.plotlist.update({'msnap':m1})
            m1 = MD_dField(self.resultsdir,fname='msnap', **kwargs)
            self.plotlist.update({'\\rho_{snap}':m1})

        #Velocity snapshots
        if 'vsnap' in (self.fieldfiles1):
            v1 = MD_pField(self.resultsdir,fname='vsnap', **kwargs)
            self.plotlist.update({'vsnap':v1})
            v1 = MD_momField(self.resultsdir,fname='vsnap', **kwargs)
            self.plotlist.update({'\\rhou_{snap}':v1})

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
        if 'vflux' in (self.fieldfiles1):
            flux1 = MD_pfluxField(self.resultsdir,'vflux', **kwargs)
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
            stress1 = MD_pfluxField(self.resultsdir,'psurface', **kwargs)
            self.plotlist.update({'psurface':stress1})

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
            v1 = MD_vField(self.resultsdir, **kwargs)
            self.plotlist.update({'u':v1})

        #Velocity snapshot
        if ('msnap' in (self.fieldfiles1) and 'vsnap' in (self.fieldfiles1)):
            v1 = MD_vField(self.resultsdir,rectype='snap', **kwargs)
            self.plotlist.update({'u_{snap}':v1})

        #Temperature
        if ('mbins' in (self.fieldfiles1) and 
            'vbins' in (self.fieldfiles1) and 
            'Tbins' in (self.fieldfiles1)):
            T1 = MD_TField(self.resultsdir, **kwargs)
            self.plotlist.update({'T':T1})

        if (len(self.plotlist) == 0):
            raise NoResultsInDir 

    def available_output_string(self):
        print('\nAvailable outputs in ' + self.resultsdir + ' include:\n')
        print('\t{0:^24s}\t{1:>10s}'.format('field', 'records'))
        print('\t{0:^24s}\t{1:^10s}'.format('-'*24, '-'*10))
        for key,field in self.plotlist.items():
            line = '\t{0:<24s}\t{1:>10d}'.format(key, field.maxrec)
            print(line)
            #print(key + ' with ', str(field[0].maxrec), ' records')

    def __repr__(self):
        return "MD_PostProc()"
    def __str__(self):
        self.available_output_string()
        return "\nInstance of MD_PostProc()"

    def read_output_file(self,keyword=None):

        """
            Reads from output file and strips out keywords using grep

        """

        tempfile = 'temp'
        self.extract_from_stdout(keyword=keyword,tempfile=tempfile)

        try:
            array = np.genfromtxt(tempfile,invalid_raise=False)
        except:
            print("ERROR in read_output_file")
            raise
            array = []

        return array

    def extract_from_stdout(self,keyword=None,tempfile='temp'):

        """
            Reads from stdout output file and writes to temp

        """
        cmd = "cat " + self.rundir + "/" + self.outputfile
        with open(tempfile,'w+') as f:
            if (keyword == None):
                args = shlex.split(cmd)
                out=subprocess.Popen(args, stdout=f)
            else:
                proc1 = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE)
                proc2 = subprocess.Popen(shlex.split('grep ' + keyword),stdin=proc1.stdout,
                             stdout=f,stderr=subprocess.PIPE)

                proc1.stdout.close() # Allow proc1 to receive a SIGPIPE if proc2 exits.
                err=proc2.communicate()

        return err

    def file_len(self,fname):
        try:
            with open(fname) as f:
                for i, l in enumerate(f):
                    pass
            return i + 1
        except:
            return 0
