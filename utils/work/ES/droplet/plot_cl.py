import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle
from scipy.optimize import curve_fit
import scipy
from skimage.morphology import skeletonize, closing
from skimage import feature

sys.path.append('../../../')
import postproclib as ppl

def update_plot(fdir,Lmin,Lmax,startrec=0,endrec=None,showplot=True):

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
    temp = fdir.split('WALLSLIDEV')[1]
    wallspeed = float(temp.split('/')[0])
    PPObj = ppl.MD_PostProc(fdir)
    rhoObj = PPObj.plotlist['rho']
    Lx = float(rhoObj.Raw.header.globaldomain1)
    Ly = float(rhoObj.Raw.header.globaldomain2)
    binsize = [float(rhoObj.Raw.header.binsize1),
               float(rhoObj.Raw.header.binsize2),
               float(rhoObj.Raw.header.binsize3)]
    halfbinsize = [0.5*i for i in binsize]
    botwall = [float(rhoObj.Raw.header.tethdistbot1),
               float(rhoObj.Raw.header.tethdistbot2),
               float(rhoObj.Raw.header.tethdistbot3)]
    topwall = [float(rhoObj.Raw.header.tethdisttop1),
               float(rhoObj.Raw.header.tethdisttop2),
               float(rhoObj.Raw.header.tethdisttop3)]
    #Add an extra 1 here as only inside domain of interest
    botwallbins = [int(botwall[i]/binsize[i])+1 for i in range(3)] 
    topwallbins = [int(topwall[i]/binsize[i])+1 for i in range(3)] 

    if endrec == None:
        endrec = rhoObj.maxrec
    rho = rhoObj.read(startrec,endrec)
    x = rhoObj.grid[0]
    y = rhoObj.grid[1][botwallbins[1]:-topwallbins[1]]
    X,Y = np.meshgrid(x,y,indexing='ij')

    #Plot contourmap of edge position against time
    if showplot:
        f = plt.figure()
        ax1 = plt.subplot2grid((2,2), (0,0), rowspan=2)
        ax2 = plt.subplot2grid((2,2), (0,1))
        ax3 = plt.subplot2grid((2,2), (1,1))

    #Take gradient of density 
    #drhodx, drhodz = np.gradient(np.mean(rho[:,0:1,:,:,0],(1,2)))
    #cm=ax1.pcolormesh(X,Y,np.sqrt(drhodx*drhodx),cmap=plt.cm.RdYlBu_r,shading='gouraud')
    meanrho = np.mean(rho[:,botwallbins[1]:-topwallbins[1],:,:,0],(2,3))

    #Plot locations of edges only using Lmin and Lmax
    fluid =  (Lmax<meanrho)
    #vapour = (Lmin<meanrho)
    fluid = scipy.ndimage.morphology.binary_closing(fluid)
    fluid = scipy.ndimage.morphology.binary_fill_holes(fluid)
    vapour = ~ fluid

    fluid.astype('bool')
    vapour.astype('bool')

    edgedensity = ~(  (fluid ^ np.roll(vapour,1,0)))# | (fluid ^ np.roll(vapour,-1,0)))
    #               | (fluid ^ np.roll(vapour,1,1)) | (fluid ^ np.roll(vapour,-1,1)) )
    #edgedensity = fluid #skeletonize(fluid)
    #edgedensity = feature.canny(fluid, sigma=3)

    x = X[edgedensity]
    y = Y[edgedensity]
    #Filter zeros, zero minimum value and plot
    indrop = 1e-5
    x = x[y>indrop]  
    y = y[y>indrop]
    x = x - x.min()
    y = y - y.min()
    if showplot:
        ax2.plot(x, y, 'bo')
    #Mirror data and plot
    x[x>x.max()/2.] = x[x<x.max()/2.]

    #Plot contour of values
    if showplot:
        try:
            cm=ax1.pcolormesh(X,Y,meanrho,cmap=plt.cm.RdYlBu_r,shading='gouraud')
            ax1.contour(X,Y,edgedensity,color='k')
        except ValueError:  #raised if `y` is empty.
            print('No interface found for Lmin and Lmax values',Lmin,Lmax)
            f.subplots_adjust(left=0.2)
            cbar_ax = f.add_axes([0.05, 0.1, 0.025, 0.8])
            f.colorbar(cm, cax=cbar_ax)
            plt.show()

        #ax.pcolormesh(meanrho>0.2,cmap=plt.cm.bone_r,alpha=0.2)
        ax1.axis('tight')
        ax1.set_xlabel('$x$')
        ax1.set_ylabel(r'$y$')

    #plot
    if showplot:
        ax2.plot(x, y, 'ro', ms=2.5)
        ax2.set_xlabel(r'$x$')
        ax2.set_ylabel(r'$y$')

    #Fit line      
    popxt, pcovx = curve_fit(linear_fn_zero_origin, x, y, (1.))
    y_linfit = linear_fn_zero_origin(x, *popxt)
    #popxt, pcovx = curve_fit(quadratic_fn_zero_origin, x, y, (1.,1.))
    #y_quadfit = quadratic_fn_zero_origin(x, *popxt)

    if showplot:
        ax3.plot(y,x,'bo')
        ax3.plot(y_linfit,x,'r-')
        #ax3.plot(y_quadfit,x,'y-')
        ax3.set_ylabel(r'$x$')
        ax3.set_xlabel(r'$y$')

        f.subplots_adjust(left=0.2)
        cbar_ax = f.add_axes([0.05, 0.1, 0.025, 0.8])
        f.colorbar(cm, cax=cbar_ax)
        plt.show()

    angle = 180-np.arctan(popxt)*180./np.pi
    print('Wall speed =',wallspeed, 'Line gradient = ', popxt, 'Angle = ', angle)

    return wallspeed, x, y_linfit, angle

#Minimin and maximum values for liquid
Lmin = 0.4
Lmax = 0.55
fbase = '/home/es205/results/droplet/tether_walls/Thompson_Robbins_1989/ewall_alkene_water/sliding/SLIDING_study/coupled_code/MD_dCSE/runs/'
fruns = ['WALLSLIDEV0000/results/','WALLSLIDEV0050/results/','WALLSLIDEV0100/results/','WALLSLIDEV0150/results/','WALLSLIDEV0200/results/']#,'WALLSLIDEV0200_highres/results/']
wlist = []; xlist = []; ylist = []; alist=[]
for f in fruns:
    fdir = fbase + f
    w, x, y, a = update_plot(fdir,Lmin,Lmax,startrec = 10, showplot=True)
    wlist.append(w)
    xlist.append(x)
    ylist.append(y)

for i in range(len(fruns)):
    plt.plot(xlist[i],ylist[i])
plt.show()

