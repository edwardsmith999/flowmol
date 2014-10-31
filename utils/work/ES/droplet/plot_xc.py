import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle
from scipy.optimize import curve_fit

sys.path.append('../../../')
import postproclib as ppl

def update_plot(fdir,Lmin,Lmax):

    def flip(a,origin=0.):
        return -a[::-1]

    def linear_fn(x, b, c):
        return b * x + c

    def linear_fn_zero_origin(x, b):
        return b * x

    def quadratic_fn(x, a, b, c):
        return a*np.power(x,2.) + b * x + c

    def quadratic_fn_zero_origin(x, a, b):
        return a*np.power(x,2.) + b * x


    #Get data
    PPObj = ppl.MD_PostProc(fdir)
    rhoObj = PPObj.plotlist['rho']
    startrec = 0
    endrec = rhoObj.maxrec
    rho = rhoObj.read(startrec,endrec)
    x = rhoObj.grid[0]
    plotfreq = rhoObj.plotfreq
    tplot=int(rhoObj.Raw.header.tplot)
    delta_t = float(rhoObj.Raw.header.delta_t)
    trec = plotfreq*tplot*delta_t
    t = np.arange(startrec*trec,(endrec+1)*trec,trec)
    X,T = np.meshgrid(x,t,indexing='ij')

    #Plot contourmap of edge position against time
    f = plt.figure()
    ax1 = plt.subplot2grid((2,2), (0,0), rowspan=2)
    ax2 = plt.subplot2grid((2,2), (0,1))
    ax3 = plt.subplot2grid((2,2), (1,1))

    #Take gradient of density 
    #drhodx, drhodz = np.gradient(np.mean(rho[:,0:1,:,:,0],(1,2)))
    #cm=ax1.pcolormesh(X,T,np.sqrt(drhodx*drhodx),cmap=plt.cm.RdYlBu_r,shading='gouraud')

    cm=ax1.pcolormesh(X,T,np.mean(rho[:,0:1,:,:,0],(1,2)),cmap=plt.cm.RdYlBu_r,shading='gouraud')
    ax1.contour(X,T,np.mean(rho[:,0:1,:,:,0],(1,2))>Lmin,color='k')
    #ax.pcolormesh(np.mean(rho[:,0:1,:,:,0],(1,2))>0.2,cmap=plt.cm.bone_r,alpha=0.2)
    ax1.axis('tight')

    #Plot locations of edges only using Lmin and Lmax
    edgedensity = ((Lmax>np.mean(rho[:,0:1,:,:,0],(1,2))) & 
                   (Lmin<np.mean(rho[:,0:1,:,:,0],(1,2)))    )
    x = X[edgedensity] - np.max(x[:])/2.
    t = T[edgedensity]

    #Plot all data
    ax2.plot(x,t,'bo')

    #Mirror data, filter zeros and plot
    indrop = 1e-5
    x[x<0.] = -x[x<0.]
    x = x[t>indrop]  
    t = t[t>indrop]
    ax2.plot(x,t,'ro',ms=2.5)

    #Fit line
    #Zero minimum value and 
    x = x - x.min() 
    popxt, pcovx = curve_fit(linear_fn_zero_origin, x, t, (1.))
    t_linfit = linear_fn_zero_origin(x, *popxt)
    popxt, pcovx = curve_fit(quadratic_fn_zero_origin, x, t, (1.,1.))
    t_quadfit = quadratic_fn_zero_origin(x, *popxt)

    ax3.plot(x,t,'bo')
    ax3.plot(x,t_linfit,'r-')
    ax3.plot(x,t_quadfit,'y-')

    f.subplots_adjust(left=0.2)
    cbar_ax = f.add_axes([0.05, 0.1, 0.025, 0.8])
    f.colorbar(cm, cax=cbar_ax)
    plt.show()


#fdir = '/home/es205/scratch/droplet/2D_e1p4/'
#fdir = '/home/es205/scratch/droplet/2D_e0p6/'
fdir = '/home/es205/codes/superspreading/coupled_code/MD_dCSE/src_code/results/'


#Minimin and maximum values for liquid
Lmin = 0.38
Lmax = 0.39
update_plot(fdir,Lmin,Lmax)

