import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle

sys.path.append('../../../')
import postproclib as ppl



#fdir = '/home/es205/scratch/droplet/2D_e1p4/'
fdir = '/home/es205/results/droplet/tether_walls/Thompson_Robbins_1989/ewall_alkene_water/sliding/SLIDING_study/coupled_code/MD_dCSE/runs/WALLSLIDEV0200_highres/results/'
xwindow = 5.
#startrecs = [6,10,50,80]; endrecs = [7,12,55,90]

#fdir = '../../../../MD_dCSE/src_code/results/'
#xwindow = 1.0


PPObj = ppl.MD_PostProc(fdir)
rhoObj = PPObj.plotlist['rho']
maxrec = rhoObj.maxrec

startrecs = [10,20,maxrec-5]; endrecs = [11,21,maxrec]

f, axs = plt.subplots(nrows=len(startrecs))

for i,ax in enumerate(axs):
    startrec = startrecs[i]; endrec = endrecs[i]

    PPObj = ppl.MD_PostProc(fdir)
    rhoObj = PPObj.plotlist['rho']
    rho = rhoObj.read(startrec,endrec)

    uvwObj = PPObj.plotlist['rho u']
    uvw = uvwObj.read(startrec,endrec)

    X,Y = np.meshgrid(uvwObj.grid[0],uvwObj.grid[1],index='ij')
    Rho = np.mean(rho[:,:,:,:,0],(2,3))
    U = np.mean(uvw[:,:,:,:,0],(2,3)); V = np.mean(uvw[:,:,:,:,1],(2,3))

    xlims = range(int(uvw.shape[0]*(5.-xwindow)/10.),int(uvw.shape[0]*(5.+xwindow)/10.))
    speed = np.sqrt(U*U + V*V)
    lw = 2.5*speed/speed.max()
    im=ax.streamplot(X[:,xlims], Y[:,xlims], U[xlims,:].T, V[xlims,:].T, 
                     color=speed[xlims,:].T, linewidth=lw[xlims,:].T, 
                     cmap=plt.cm.RdYlBu_r, density=[3., 3])
    ax.axis('scaled')

f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
f.colorbar(im.lines, cax=cbar_ax)
plt.show()
