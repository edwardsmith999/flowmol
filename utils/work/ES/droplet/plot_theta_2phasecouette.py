import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle
import scipy
from scipy.optimize import curve_fit
import MDAnalysis

sys.path.append('../../../')
import postproclib as ppl
import gen_surface

#fdir = '/home/es205/scratch/droplet/2D_e1p4/'
#xwindow = 3.0

#fdir = '/home/es205/scratch/droplet/2D_e0p6/'
#xwindow = 1.0

fdir = '/home/es205/results/droplet/tether_walls/Thompson_Robbins_1989/ewall_alkene_water/sliding/SLIDING_study/coupled_code/MD_dCSE/runs/WALLSLIDEV0200/results/'
xwindow = 5.


def get_fluid_vapour(rho,Lmax):

    fluid =  (Lmax<rho)
    fluid = scipy.ndimage.morphology.binary_closing(fluid)
    fluid = scipy.ndimage.morphology.binary_fill_holes(fluid)
    vapour = ~ fluid

    return fluid, vapour

def get_edge_cells(density,Lmin,Lmax):
    #Return a mask array of edgecells
    edgemask = ((Lmax>density) &(Lmin<density))
    #return arrays with cells on edge
    return edgemask, np.where(edgemask)

        
def get_edge_molecules(mols,cellloc,halfcellsize=[0.5,0.5]):
    Lmax = cellloc + halfcellsize
    Lmin = cellloc - halfcellsize
    return ((Lmax>mols) &(Lmin<mols))

#Get data
PPObj = ppl.MD_PostProc(fdir)
rhoObj = PPObj.plotlist['rho']
DCDObj = MDAnalysis.coordinates.DCD.DCDReader(fdir+"vmd_out.dcd")
startrec = 0
endrec = rhoObj.maxrec
x = rhoObj.grid[0]
z = rhoObj.grid[1]
X,Z = np.meshgrid(x,z,indexing='ij')

intialstep = int(rhoObj.Raw.header.initialstep)
plotfreq = rhoObj.plotfreq
tplot=int(rhoObj.Raw.header.tplot)
delta_t = float(rhoObj.Raw.header.delta_t)

minxi = int(x.shape[0]*(5.-xwindow)/10.)
maxxi = int(x.shape[0]*(5.+xwindow)/10.)-1
minx = x[minxi]; maxx = x[maxxi]; 
Lx = float(rhoObj.Raw.header.globaldomain1)
Ly = float(rhoObj.Raw.header.globaldomain2)
binsize = [float(rhoObj.Raw.header.binsize1),
           float(rhoObj.Raw.header.binsize2),
           float(rhoObj.Raw.header.binsize3)]
halfbinsize = [0.5*i for i in binsize]
xlims = range(minxi,maxxi)
Lmin = 0.05; Lmax = 0.7

botwall = [float(rhoObj.Raw.header.tethdistbot1),
           float(rhoObj.Raw.header.tethdistbot2),
           float(rhoObj.Raw.header.tethdistbot3)]
topwall = [float(rhoObj.Raw.header.tethdisttop1),
           float(rhoObj.Raw.header.tethdisttop2),
           float(rhoObj.Raw.header.tethdisttop3)]
botwallbins = [int(botwall[i]/binsize[i]) for i in range(3)] 
topwallbins = [int(topwall[i]/binsize[i]) for i in range(3)] 

print(topwallbins,botwallbins)

f = plt.figure()
#plt.ion()
#plt.show()
step = 10
for rec in range(startrec,endrec,step):
    print('Record ' + str(rec) + ' of ' + str(endrec))
    rho = rhoObj.read(rec,rec+step)
    vmd = DCDObj[rec*100][:]

    x = X[xlims,botwallbins[1]:-topwallbins[1]]
    z = Z[xlims,botwallbins[1]:-topwallbins[1]]
    density = np.mean(rho[xlims,botwallbins[1]:-topwallbins[1],:,:,0],(2,3))
    edge, index = get_edge_cells(density,Lmin,Lmax)

    fluid,vapour = get_fluid_vapour(density,Lmax)
    fluid.astype('bool')
    vapour.astype('bool')
    edgedensity = ~(  (fluid ^ np.roll(vapour,1,0)))

    #Plot contourmap of edge position in x against z
    ax1 = plt.subplot2grid((3,1), (0,0))#, rowspan=2)
    ax2 = plt.subplot2grid((3,1), (1,0))
    ax3 = plt.subplot2grid((3,1), (2,0))
    #ax4 = plt.subplot2grid((2,2), (1,1))

    cm = ax1.pcolormesh(x,z,density,
                        cmap=plt.cm.RdYlBu_r,shading='gouraud')
    #ax1.contour(x,z,edgedensity,color='k')
    #ax1.contour(x,z,edge,color='k')
    ax2.plot(X[edge],Z[edge],'o')
    molwindow = (((maxx-0.5*Lx)>vmd[:,0]) & ((minx-0.5*Lx)<vmd[:,0]))

    ax3.plot(vmd[molwindow,0],vmd[molwindow,1],'ko',ms=0.1)

    #surfacemol = np.genfromtxt('../../../../MD_dCSE/src_code/fort.451')

    timestep = intialstep + (rec*plotfreq-1)*tplot
    print(timestep)
    #rectime = surfacemol[:,0] == timestep
    #ax3.plot((surfacemol[rectime,1]-1)*binsize[0]-0.5*Lx,
    #         (surfacemol[rectime,2]-1)*binsize[1]-0.5*Ly,'bs')
    #ax3.plot(surfacemol[rectime,7],surfacemol[rectime,8],'rx',ms=0.1)

    #surfacenormal = np.genfromtxt('../../../../MD_dCSE/src_code/fort.452')
    #ax3.quiver(surfacenormal[:,1],surfacenormal[:,2],surfacenormal[:,4],surfacenormal[:,5])

    #surfacevar = np.genfromtxt('../../../../MD_dCSE/src_code/fort.1984')
    #rectime = surfacevar[:,0] == timestep
#    print(rectime)
#    #surfacevar = surfacevar[rectime,:]
#    varmax = 10.
#    varmin = 0.
#    vartest = surfacevar[:,9]# np.mean(d[:,6:8],1)
#    varfilter = ((varmax>vartest) & ((varmin<vartest)))
#    try:
#        cm=ax3.scatter(surfacevar[varfilter,1],surfacevar[varfilter,2],c=vartest[varfilter],cmap=plt.cm.RdYlBu_r)
#    except:
#        pass
##    for i in index[0]:
#        for j in index[1]:
#            get_edge_molecules(vmd[molwindow,0:1],[x[i],z[j]],halfcellsize=halfbinsize)
    

    #ax3.plot(vmd[:,0],vmd[:,1],'o',ms=0.1)
    #img = gen_surface.gen_surface(X[edge],Z[edge],np.zeros(Z[edge[:,:,0]].shape))
    #ax3.imshow(img)

    ax1.set_aspect('equal')
    ax2.set_aspect('equal')
    ax3.set_aspect('equal')
    ax1.set_ylim((0,z.max()))
    ax2.set_ylim((0,z.max()))
    ax3.set_ylim((-0.5*z.max(),0.5*z.max()))
    ax1.set_xlim((minxi,maxxi))
    ax2.set_xlim((minx, maxx))
    ax3.set_xlim((minx-0.5*Lx, maxx-0.5*Lx))

    ax1.autoscale(tight=True)
    f.subplots_adjust(left=0.2)
    cbar_ax = f.add_axes([0.05, 0.1, 0.025, 0.8])
    f.colorbar(cm, cax=cbar_ax)

    plt.show()
    #plt.pause(0.1)
    ax1.cla()
    ax2.cla()
    ax3.cla()

