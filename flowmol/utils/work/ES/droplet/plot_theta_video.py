import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle
from scipy.optimize import curve_fit

sys.path.append('../../../')
import postproclib as ppl

#fdir = '/home/es205/scratch/droplet/2D_e1p6/'
fdir = '/home/es205/scratch/droplet/2D_e0p6/'

#Get data
PPObj = ppl.MD_PostProc(fdir)
rhoObj = PPObj.plotlist['rho']
startrec = 0
endrec = rhoObj.maxrec
x = rhoObj.grid[0]
z = rhoObj.grid[1]
X,Z = np.meshgrid(x,z,indexing='ij')

xwindow = 1.0
minxi = int(x.shape[0]*(5.-xwindow)/10.)
maxxi = int(x.shape[0]*(5.+xwindow)/10.)
minx = x[minxi]; maxx = x[maxxi]; 
xlims = range(minxi,maxxi)
Lmin = 0.15; Lmax = 0.5

f = plt.figure()
outname = "Video"
for rec in range(startrec,endrec):
    print(rec)
    rho = rhoObj.read(rec,rec)
    edgedensity = ((Lmax>np.mean(rho[:,:,:,:,0],3)) & 
                   (Lmin<np.mean(rho[:,:,:,:,0],3))    )

    #Plot contourmap of edge position in x against z
    ax1 = plt.subplot2grid((1,1), (0,0))


    cm=ax1.pcolormesh(X[xlims,:],Z[xlims,:],np.mean(rho[xlims,:,:,:,0],(2,3)),cmap=plt.cm.RdYlBu_r,shading='gouraud')

    ax1.set_aspect('equal')
    ax1.autoscale(tight=True)
    #ax1.set_ylim((0,z.max()))
    #ax1.set_xlim((minxi,maxxi))
    #plt.show()
    f.savefig(outname + '{:05d}'.format(rec)+'.png')


cmd = "ffmpeg -f image2 -r 10 -i " + outname + "%05d.png -sameq anim.mov -pass 2"

os.system(cmd)
