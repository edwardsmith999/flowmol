import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import math as maths

from MDFields import *
from HeaderData import *


class MD_PostProc:

    """ 
        Post processing class for MD runs
    """

    def __init__(self,resultsdir):
        self.resultsdir = resultsdir
        self.plotlist = {}
        self.error = {}
        self.name = self.resultsdir.split('/')[-2]

        self.cmap = plt.cm.RdYlBu_r


        # Check directory exists before instantiating object and check 
        # which files associated with plots are in directory
        self.potentialfiles = ( "mslice", "mbins", "msnap","vslice", "vbins", 
                                "vsnap","pvirial", "pVA", "pVA_k","pVA_c", 
                                "visc", "mflux","vflux", "pplane", "psurface",
                                "esnap", "eflux", "eplane","esurface", 
                                "viscometrics", "rdf", "rdf3d", "ssf", "Fext",
                                "Tbins", "vPDF" )        

        if os.path.isdir(self.resultsdir):
            self.fieldfiles1 = list(set(os.listdir(self.resultsdir)) & set(self.potentialfiles))
            Header1 = HeaderData(open(self.resultsdir + 'simulation_header','r'))
        else:
            print("Directory " +  self.resultsdir + " not found")
            return

        nbintop = round(float(Header1.tethdisttop2)/float(Header1.binsize2))
        nbinbot = round(float(Header1.tethdistbot2)/float(Header1.binsize2))
        binlimits = [None,(int(nbinbot),int(Header1.gnbins2)-int(nbintop)),None]

        #Mass
        if 'mbins' in (self.fieldfiles1):
            m1 = MD_mField(self.resultsdir)
            self.plotlist.update({'mbins':m1})

        #Momentum
        if 'vbins' in (self.fieldfiles1):
            M1 = MD_pField(self.resultsdir)
            self.plotlist.update({'vbins':M1})

        #Kinetic energy
        if 'Tbins' in (self.fieldfiles1):
            KE1 = MD_KEField(self.resultsdir)
            self.plotlist.update({'Tbins':KE1})

        #Mass snapshots
        if 'msnap' in (self.fieldfiles1):
            m1 = MD_mField(self.resultsdir,fname='msnap')
            self.plotlist.update({'msnap':m1})

        #Velocity snapshots
        if 'vsnap' in (self.fieldfiles1):
            v1 = MD_pField(self.resultsdir,fname='vsnap')
            self.plotlist.update({'vsnap':v1})

        #VA stress
        if 'pVA' in (self.fieldfiles1):
            P1 = MD_pVAField(self.resultsdir)
            self.plotlist.update({'pVA':P1})
        elif 'pVA_k' in (self.fieldfiles1):
            P1 = MD_pVAField(self.resultsdir,fname='pVA_k')
            self.plotlist.update({'pVA_k':P1})

            P1 = MD_pVAField(self.resultsdir,fname='pVA_c')
            self.plotlist.update({'pVA_c':P1})

        #CV fluxes
        if 'vflux' in (self.fieldfiles1):
            flux1 = MD_pfluxField(self.resultsdir,'vflux')
            self.plotlist.update({'vflux':flux1})

        #CV stresses
        if 'psurface' in (self.fieldfiles1):
            stress1 = MD_pfluxField(self.resultsdir,'psurface')
            self.plotlist.update({'psurface':stress1})

        #Velocity
        if ('mbins' in (self.fieldfiles1) and 'vbins' in (self.fieldfiles1)):
            v1 = MD_vField(self.resultsdir)
            self.plotlist.update({'vel':v1})

        #Velocity snapshot
        if ('msnap' in (self.fieldfiles1) and 'vsnap' in (self.fieldfiles1)):
            v1 = MD_vField(self.resultsdir,rectype='snap')
            self.plotlist.update({'vel_snap':v1})

        #Temperature
        if ('mbins' in (self.fieldfiles1) and 
            'vbins' in (self.fieldfiles1) and 
            'Tbins' in (self.fieldfiles1)):
            T1 = MD_TField(self.resultsdir)
            self.plotlist.update({'Temprt':T1})

        # ============================================================================
        # Useful Parameters
        self.nbins = int(Header1.gnbins1)*int(Header1.gnbins2)*int(Header1.gnbins3)
        self.binsize = float(Header1.binsize1)*float(Header1.binsize2)*float(Header1.binsize3)
        inspectfile = 'mbins'
        figname = inspectfile
        filebytes = os.path.getsize(self.resultsdir+inspectfile)
        inspectbytesperbin = 4
        self.maxrec = filebytes / (inspectbytesperbin*self.nbins) 

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


if __name__ == "__main__":
    import matplotlib

    fdir = '/home/es205/scratch/Re400/iter390870_to_667054/'
    #fdir = '../MD_dCSE/src_code/results/'
    PP_fielddict = MD_PostProc(fdir)
    print(PP_fielddict)

    #plot data
    #var = raw_input("Please enter choice of field: ")
    var = 'velocity'
    field = PP_fielddict.plotlist[var]
    startrec = 5

    #nrecs = 1
    #spectrum = field.power_spectrum(startrec=startrec,endrec=startrec+nrecs-1,preavgaxes=(3),fftaxes=(0),postavgaxes=(1,2))
    #plt.plot(spectrum[:,0],'-o',label=str(nrecs))

    var = 'velocity'
    field = PP_fielddict.plotlist[var]

    x, z, v = field.contour(axes=(0,2),startrec=startrec,endrec=startrec+80,binlimits=[None,(128,128),None])
    cmap = plt.cm.RdYlBu_r
    contour = plt.contourf(x,z,v[:,:,0],40,cmap=cmap)
    print(v.shape)
    plt.colorbar(contour)
    plt.show()
    #quit()

    #setup figure
    #fig = plt.figure()
    #ax = fig.add_subplot(111)


    #for nrecs in [20,100]:
    spectrum = field.power_spectrum(startrec=startrec,endrec=startrec+1-1,preavgaxes=(1,3),fftaxes=(0),postavgaxes=(0),binlimits=[None,(128,128),None])
    #plt.plot(spectrum[:,0],'-o',label=str(nrecs*1600))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.plot(spectrum[:,0])
    spectrum = field.power_spectrum(startrec=startrec,endrec=startrec+1-1,preavgaxes=(1,3),fftaxes=(0),postavgaxes=(2),binlimits=[None,(128,128),None])
    plt.plot(spectrum[:,0])
    #contour = ax.pcolormesh(spectrum[:,:,0],cmap=cmap,norm=matplotlib.colors.LogNorm())

#    spectrum = field.power_spectrum(startrec=startrec,endrec=startrec+100-1,preavgaxes=(1,3),fftaxes=(0,2),postavgaxes=(),binlimits=[None,(128,128),None])
#    #plt.plot(spectrum[:,0],'-o',label=str(nrecs*1600))
#    #fig = plt.figure()
#    #ax = fig.add_subplot(111)
#    contour = ax.contour(spectrum[:,:,0],color='k',norm=matplotlib.colors.LogNorm())

#    plt.colorbar(contour)
    plt.show()
    plt.clf()

    #plt.legend()
    #plt.show()


