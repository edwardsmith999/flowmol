import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle
import MDAnalysis

sys.path.append('../../../')
import postproclib as ppl

#Get MD profile data
def get_surface_tension(fdir='./results/', plotstuff=False):

    def __get_stuff(axis, startrec, endrec):

        x, rho = PPObj.plotlist['rho'].profile(axis=0,startrec=startrec, endrec = endrec)
        x, T = PPObj.plotlist['T'].profile(axis=0,startrec=startrec, endrec = endrec)
        x, pVA_c = PPObj.plotlist['pVA_c'].profile(axis=0,startrec=startrec, endrec = endrec)
        x, pVA_k = PPObj.plotlist['pVA_k'].profile(axis=0,startrec=startrec, endrec = endrec)
        x, pVA = PPObj.plotlist['pVA'].profile(axis=0,startrec=startrec, endrec = endrec)

        x, psurface = PPObj.plotlist['psurface'].profile(axis=0,startrec=startrec, endrec = endrec)
        x, vflux = PPObj.plotlist['vflux'].profile(axis=0,startrec=startrec, endrec = endrec)
        CV_pressure = psurface + vflux

        density = float(PPObj.plotlist['psurface'].Raw.header.liquid_density)
        rcutoff = float(PPObj.plotlist['psurface'].Raw.header.rcutoff)
        sLRC = 8.*np.pi*density**2.*(4./(9.*np.power(rcutoff,9.)) - 2./(3.*np.power(rcutoff,3.)))

        return x, rho, T, pVA_c, pVA_k, pVA, psurface, vflux, CV_pressure, sLRC

    #Get y MD profile
    PPObj = ppl.MD_PostProc(fdir)
    startrec = 20; endrec = PPObj.plotlist['psurface'].maxrec-3
    nrecs = (endrec-startrec)+1

    #Get average temperature in run
    mp = np.genfromtxt(fdir + 'macroscopic_properties', delimiter=';', names=True)
    Nave = int(PPObj.plotlist['psurface'].Raw.header.Nmass_ave)
    temperature = np.mean(mp['Temp'][startrec*Nave:])

    #Load all required items
    x, rho, T, pVA_c, pVA_k, pVA, psurface, vflux, CV_pressure, sLRC = __get_stuff(axis=0, startrec=startrec, endrec = endrec)

    #get surface tension
    T_estimate = np.mean(T[rho > 0.4])
    integrand_VA = pVA[:,0] - 0.5*(pVA[:,4] + pVA[:,8])
    integrand_CV = 0.5*(  CV_pressure[:,0] - 0.5*(CV_pressure[:,4] + CV_pressure[:,8]) 
                        + CV_pressure[:,9] - 0.5*(CV_pressure[:,13] + CV_pressure[:,17]))

    #PLOT pressures as a funciton of x
    if plotstuff:

        #Not we have a psf file, better to use this (access to bond, etc)
        try:
            mols = MDAnalysis.Universe(fdir+'polymer_topol.psf', fdir+'vmd_out.dcd')
            DCDObj = mols.trajectory
        except:
            DCDObj = MDAnalysis.coordinates.DCD.DCDReader(fdir+"vmd_out.dcd")

        try:
            maxvmdrec = len(DCDObj)-1
            vmd = DCDObj[maxvmdrec][:] 
            vmd[:,0] = vmd[:,0] + 0.5*float(PPObj.plotlist['psurface'].Raw.header.globaldomain1)
            vmd[:,1] = vmd[:,1]/10
            plt.plot( vmd[:,0],vmd[:,1], 'k.', ms=0.3)
        except IndexError:
            print(endrec*Nave, 'is outside of vmd range')
        

        comp = [0, 4, 8]
        plt.plot(x, rho[:,0], 'k-')
        plt.plot(x, T[:,0], 'c-')
        plt.plot(x, pVA_c[:,comp], 'r')
        plt.plot(x, pVA_k[:,comp], 'b')
        plt.plot(x, pVA[:,comp], 'y')

        comp = [0, 4, 8]
        plt.plot(x, psurface[:,comp], 'ro')
        plt.plot(x, vflux[:,comp], 'bo')
        plt.plot(x, CV_pressure[:,comp], 'yo')

        comp = [9, 13, 17]
        plt.plot(x, psurface[:,comp], 'rx')
        plt.plot(x, vflux[:,comp], 'bx')
        plt.plot(x, CV_pressure[:,comp], 'yx')
        plt.show()

        plt.plot(integrand_VA)
        plt.plot(integrand_CV)
        plt.show()

    dx = float(PPObj.plotlist['pVA'].Raw.header.binsize1)
    gamma_VA = 0.5 * np.trapz(integrand_VA, x=None, dx=dx)
    gamma_CV = 0.5 * np.trapz(integrand_CV, x=None, dx=dx)

    #Try getting per value
#    gamma_VA_hist = []; gamma_CV_hist = []
#    for rec in range(startrec, endrec):
#        x, rho, T, pVA_c, pVA_k, pVA, psurface, vflux, CV_pressure, sLRC = __get_stuff(axis=0,startrec=rec, endrec = rec)

#        integrand_VA = pVA[:,0] - 0.5*(pVA[:,4] + pVA[:,8])
#        integrand_CV = 0.5*(  CV_pressure[:,0] - 0.5*(CV_pressure[:,4] + CV_pressure[:,8]) 
#                            + CV_pressure[:,9] - 0.5*(CV_pressure[:,13] + CV_pressure[:,17]))

#        dx = float(PPObj.plotlist['pVA'].Raw.header.binsize1)
#        gamma_VA_hist.append(0.5 * np.trapz(integrand_VA, x=None, dx=dx))
#        gamma_CV_hist.append(0.5 * np.trapz(integrand_CV, x=None, dx=dx))

#    print(np.mean(np.array(gamma_VA_hist)), np.mean(np.array(gamma_CV_hist)))

    return temperature, T_estimate, gamma_VA, gamma_CV, sLRC


if __name__ == "__main__":
    print(get_surface_tension(fdir = '../../../../MD_dCSE/src_code/results/', plotstuff=True))



#    print(np.mean(np.array(gamma_VA_hist)), np.mean(np.array(gamma_CV_hist)))

#comp = [1, 2, 3]

#x, pVA_c = PPObj.plotlist['pVA_c'].profile(axis=0,startrec=startrec, endrec = endrec)
#plt.plot(x, pVA_c[:,comp], 'g')

#x, pVA_k = PPObj.plotlist['pVA_k'].profile(axis=0,startrec=startrec, endrec = endrec)
#plt.plot(x, pVA_k[:,comp], 'k')

#x, psurface = PPObj.plotlist['psurface'].profile(axis=0,startrec=startrec)
#plt.plot(x, psurface[:,comp], 'go')

#x, vflux = PPObj.plotlist['vflux'].profile(axis=0,startrec=startrec)
#plt.plot(x, vflux[:,comp], 'ko')

#CV_pressure = psurface + vflux
#plt.plot(x, CV_pressure[:,comp], 'yo')
#plt.show()
#comp = range(9)
#x, pVA = PPObj.plotlist['pVA'].profile(axis=0,startrec=startrec, endrec = endrec)
#plt.plot(x, pVA[:,comp], 'r')
#plt.show()





#CV_pressure = psurface + vflux
#plt.plot(x, CV_pressure[:,comp], 'ro')


