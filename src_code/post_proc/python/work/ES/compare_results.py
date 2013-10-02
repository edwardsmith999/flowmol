import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#sys.path.insert(0, './../../')
sys.path.insert(0,'/home/es205/codes/coupled/MD_dCSE/src_code/post_proc/python')

from MDFields import *
from MDPlotData import MD_PlotData
from HeaderData import *

class CompareResults():

    def __init__(self):
        self.plotlist = {}
        self.error = {}

    def plot_and_save(  self,
                        rmin,
                        rmax,
                        cnt,
                        outputaxis_handle=None,
                        filename='figname'):

        # New axes on cleared figure
        if outputaxis_handle == None:
            fig = plt.figure()
            ax = fig.add_subplot(221, frame_on=False)

        #Plot all possible comparisons
        error_found = False
        for key,value in self.plotlist.items():
            if (value[0].get_bins(rmin,rmax)[0].shape == value[1].get_bins(rmin,rmax)[0].shape):
                error = value[0].get_bins(rmin,rmax)[0] - value[1].get_bins(rmin,rmax)[0]
                totalerror = error.sum()+error.max()+error.min()
                if (totalerror > 1e-5):
                    print("At record number " + str(rmin) + " Error in " + key + " = " + str(totalerror))
                    error_found = True

                    #Attempt to Plot slices of the errors
                    #CONTOUR
                    #x = np.arange(-2.0, 2.0, 4.0/error.shape[0])
                    #y = np.arange(-2.0, 2.0, 4.0/error.shape[1])
                    #X, Y = np.meshgrid(x, y)
                    #CS = ax.contourf(X, Y, Z/Z.max())
                    figname = self.name + key

                    #Attempt to Plot slices of the errors in 3D
                    ax = fig.add_subplot(221, frame_on=False)
                    Z = error[:,:,error.shape[0]/4.0,1]
                    CS = ax.imshow(Z, extent=[0, 1, 0, 1])
                    plt.colorbar(CS); CS.set_clim(-1.0,1.0)
                    plt.title('xy slice')

                    ax = fig.add_subplot(222, frame_on=False)
                    Z = error[:,error.shape[0]/4.0,:,1]
                    CS = ax.imshow(Z, extent=[0, 1, 0, 1])
                    plt.colorbar(CS); CS.set_clim(-1.0,1.0)
                    plt.title('xz slice')

                    ax = fig.add_subplot(223, frame_on=False)
                    Z = error[error.shape[0]/4.0,:,:,1]
                    CS = ax.imshow(Z, extent=[0, 1, 0, 1])
                    plt.colorbar(CS); CS.set_clim(-1.0,1.0)
                    plt.title('yz slice')

                    ax = fig.add_subplot(224, frame_on=False)
                    Z = error.mean((1,2))[:,:]
                    CS = ax.imshow(Z, extent=[0, 1, 0, 1])
                    plt.colorbar(CS); CS.set_clim(-1.0,1.0)
                    plt.title('sum over x and y with z vs 4th index')

                    # Save figure and clear
                    plt.savefig(figname+'.'+"%05d"%rmin+'.png')
                    plt.clf()
            else:
                #print("At record number " + str(rmin) + " Array sizes differ for " + key)
                error_found = True

        return error_found




    def run(self,run1,run2):

        self.files1 = run1.rundir
        self.files2 = run2.rundir
        self.name = self.files2.split('/')[-2]
        print("Comparing results in directory " +  self.files1 + " \n"
              " with results in directory " +  self.files2 )

        self.files1 = self.files1 + '/results/'
        self.files2 = self.files2 + '/results/'

        # Check directory exists before instantiating object and check 
        # which files associated with plots are in directory
        self.potentialfiles = ( "mslice", "mbins", "msnap","vslice", "vbins", 
                                "vsnap","pvirial", "pVA", "pVA_k","pVA_c", 
                                "visc", "mflux","vflux", "pplane", "psurface",
                                "esnap", "eflux", "eplane","esurface", 
                                "viscometrics", "rdf", "rdf3d", "ssf", "Fext ",
                                "Tbins" )        
        if os.path.isdir(self.files1):
            self.fieldfiles1 = list(set(os.listdir(self.files1)) & set(self.potentialfiles))
        if os.path.isdir(self.files2):
            self.fieldfiles2 = list(set(os.listdir(self.files2)) & set(self.potentialfiles))

        # Global Objects
        Header1 = HeaderData(open(self.files1 + 'simulation_header','r'))
        Header2 = HeaderData(open(self.files2 + 'simulation_header','r'))

        #Mass
        if 'mbins' in (self.fieldfiles1 and self.fieldfiles2):
            m1 = MassBins(self.files1)
            m2 = MassBins(self.files2)
            self.plotlist.update({'mbins':[m1,m2]})

        #Momentum
        if 'vbins' in (self.fieldfiles1 and self.fieldfiles2):
            M1 = MomBins(self.files1)
            M2 = MomBins(self.files2)
            self.plotlist.update({'vbins':[M1,M2]})

        #Temperature
        if 'Tbins' in (self.fieldfiles1 and self.fieldfiles2):
            T1 = KEBins(self.files1)
            T2 = KEBins(self.files2)
            self.plotlist.update({'Tbins':[T1,T2]})

        #VA stress
        if 'pVA' in (self.fieldfiles1 and self.fieldfiles2):
            P1 = PBins(self.files1)
            P2 = PBins(self.files2)
            self.plotlist.update({'pVA':[P1,P2]})

        #CV fluxes
        if 'vflux' in (self.fieldfiles1 and self.fieldfiles2):
            flux1 = CV(self.files1,'vflux')
            flux2 = CV(self.files2,'vflux')
            self.plotlist.update({'CV_fluxes':[flux1,flux2]})

        #CV stresses
        if 'psurface' in (self.fieldfiles1 and self.fieldfiles2):
            stress1 = CV(self.files1,'psurface')
            stress2 = CV(self.files2,'psurface')
            self.plotlist.update({'CV_stresses':[stress1,stress2]})

        # ============================================================================
        # Useful Parameters
        nbins = int(Header1.gnbins1)*int(Header1.gnbins2)*int(Header1.gnbins3)
        inspectfile = 'mbins'
        figname = inspectfile
        filebytes = os.path.getsize(self.files1+inspectfile)
        inspectbytesperbin = 4
        maxrec = filebytes / (inspectbytesperbin*nbins) 

        # ============================================================================
        # SWEEP THROUGH ALL VALUES
        drec = 1
        cnt = 0
        error_found = False
        for rec in range(0,maxrec-1,drec):
            rmin = rec
            rmax = rec + drec
            error_found = self.plot_and_save(rmin,rmax,cnt)
            cnt += 1
            if error_found == False:
                print("No Error in " + ' '.join(self.plotlist.keys()) + 
                      " between " + str(rmin) + ' and ' +  str(rmax))
            #else:
                #print("********* Error in " + ' '.join(self.plotlist.keys()) + 
                #      " between " + str(rmin) + ' and ' +  str(rmax)+ " **********")
               # print("This may be a result of different numbers of bins between cases!")
