import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle
from scipy.optimize import curve_fit
import scipy
from skimage.morphology import skeletonize, closing
from skimage import feature
import traceback
import re
from operator import itemgetter, attrgetter, methodcaller

sys.path.append('../../../')
import postproclib as ppl
from postproclib.pplexceptions import NoResultsInDir

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

def update_plot(fdir,Lmax,startrec=0,endrec=None,showplot=True):

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
    temp = temp.split('/')[0]
    wallspeed = float(temp[0] + '.' + temp[1:])
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
    meanrho = np.mean(rho[:,botwallbins[1]:-topwallbins[1],:,:,0],(2,3))

    #Plot locations of edges only using Lmax
    fluid =  (Lmax<meanrho)
    fluid = scipy.ndimage.morphology.binary_closing(fluid)
    fluid = scipy.ndimage.morphology.binary_fill_holes(fluid)
    vapour = ~ fluid

    fluid.astype('bool')
    vapour.astype('bool')

    edgedensity = ~(  (fluid ^ np.roll(vapour,1,0)))# | (fluid ^ np.roll(vapour,-1,0)))
    #               | (fluid ^ np.roll(vapour,1,1)) | (fluid ^ np.roll(vapour,-1,1)) )
    #edgedensity = fluid #skeletonize(fluid)
    #edgedensity = feature.canny(fluid, sigma=3)

    #Plot contour of values
    if showplot:
        try:
            cm=ax1.pcolormesh(X,Y,meanrho,cmap=plt.cm.RdYlBu_r,shading='gouraud')
            ax1.contour(X,Y,edgedensity,color='k')
        except ValueError:  #raised if `y` is empty.
            print('No interface found for and Lmax values',Lmax)
            f.subplots_adjust(left=0.2)
            cbar_ax = f.add_axes([0.05, 0.1, 0.025, 0.8])
            f.colorbar(cm, cax=cbar_ax)
            plt.show()

        #ax.pcolormesh(meanrho>0.2,cmap=plt.cm.bone_r,alpha=0.2)
        ax1.axis('tight')
        ax1.set_xlabel('$x$')
        ax1.set_ylabel(r'$y$')

    x = X[edgedensity]
    y = Y[edgedensity]

    #Filter zeros, zero minimum value and plot
    indrop = 1e-5
    x = x[y>indrop]  
    y = y[y>indrop]
    x = x - x.min()
    y = y - y.min()
    if showplot:
        ax2.plot(x, y, 's',alpha=0.3)
        #ax2.plot(x[x>x.max()/2.], y[x>x.max()/2.], 'rx')
        #ax2.plot(x[x<x.max()/2.], y[x<x.max()/2.], 'rs',alpha=0.4)

    #Split data by surface
    x1 = []; x2 = []
    y1 = []; y2 = []
    mid = x.mean()
    for i in range(x.shape[0]):
        if x[i] > mid:
            x2.append(x[i])
            y2.append(y[i])
            #shift = 2.*(x[i]-mid)
            #x1.append(-x[i] + 2.*mid + 2.*(x[i]-mid))
        else:
            x1.append(x[i])
            y1.append(y[i])
    
    x1 = np.array(x1)
    x2 = np.array(x2)
    y1 = np.array(y1)
    y2 = np.array(y2)

    #x[x>x.max()/2.] = x[x<x.max()/2.]

    #plot
    if showplot:
        ax2.plot(x1, y1, 'bo', ms=2.5)
        ax2.plot(x2, y2, 'ro', ms=2.5)
        ax2.set_xlabel(r'$x$')
        ax2.set_ylabel(r'$y$')

    #Fit line      
    popxt1, pcovx1 = curve_fit(linear_fn, x1, y1, (1., 0.))
    y1_linfit = linear_fn(x1, *popxt1)
    popxt2, pcovx2 = curve_fit(linear_fn, x2, y2, (1., 0.))
    y2_linfit = linear_fn(x2, *popxt2)

#    popxtq, pcovxq = curve_fit(quadratic_fn_zero_origin, x1, y1, (1.,1.))
#    y1_quadfit = quadratic_fn_zero_origin(x1, *popxtq)

    if showplot:
        ax3.plot(y1,x1,  'bo')
        ax3.plot(y1_linfit,x1,'b-')
        #ax3.plot(y1_quadfit,x1,'b--')

        ax3.plot(y2,x2, 'ro')
        ax3.plot(y2_linfit,x2,'r-')

        ax3.set_ylabel(r'$x$')
        ax3.set_xlabel(r'$y$')

        f.subplots_adjust(left=0.2)
        cbar_ax = f.add_axes([0.05, 0.1, 0.025, 0.8])
        f.colorbar(cm, cax=cbar_ax)
        plt.show()

    popxt = 0.5*(popxt1[0] + popxt2[0])
    angle = 180-np.arctan(popxt)*180./np.pi
    print('Wall speed =',wallspeed, 'Line gradient = ', popxt, 'Angle = ', angle)
    y_linfit = 0.5*(y1_linfit[0] + y2_linfit[0])
    return wallspeed, x, y_linfit, angle

#Minimin and maximum values for liquid
Lmax = 0.4
fbase = '/home/es205/results/droplet/tether_walls/Thompson_Robbins_1989/LJ_wall/coupled_code/MD_dCSE/runs/'

#Get directories and order by speed
fruns = [ v for v in os.listdir(fbase)  if v.endswith('0') ]
temps = [fdir.split('WALLSLIDEV')[1] for fdir in fruns]
wallspeeds = [temp.split('/')[0] for temp in temps]
wallspeeds = [float(wallspeed[0] + '.' + wallspeed[1:]) for wallspeed in wallspeeds]
out = sorted(zip(fruns, wallspeeds), key=itemgetter(1))
fruns = [x[0]  + '/results/' for x in out]
wallspeeds = sorted(wallspeeds)
print(fruns)

#Loop through directories
wlist = []; xlist = []; ylist = []; alist=[]
for i, f in enumerate(fruns):
    fdir = fbase + f 
    print('Trying fdir = ', fdir, wallspeeds[i], wallspeeds[i]< 0.1)
    if wallspeeds[i] < 0.1:
        try:
            w, x, y, a = update_plot(fdir,Lmax,startrec = 10, showplot=True)
        except NoResultsInDir:
            print(traceback.format_exc())
            w = 0; x = 0; y = 0; a = 0
        except:
            raise
            
        wlist.append(w)
        xlist.append(x)
        ylist.append(y)
        alist.append(a)
wlist = np.array(wlist)
alist = np.array(alist)
alist[0] = 90; alist[1] = 90
plt.plot(wlist,alist,'o')
plt.show()

