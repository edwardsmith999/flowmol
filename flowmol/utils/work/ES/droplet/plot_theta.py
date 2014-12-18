import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle
from scipy.optimize import curve_fit
import MDAnalysis

sys.path.append('../../../')
import postproclib as ppl
import gen_surface

#fdir = '/home/es205/scratch/droplet/2D_e1p4/'
#xwindow = 3.0

#fdir = '/home/es205/scratch/droplet/2D_e0p6/'
#xwindow = 1.0

fdir = '/home/es205/codes/superspreading/coupled_code/MD_dCSE/src_code/results/'
xwindow = 3.0


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
maxxi = int(x.shape[0]*(5.+xwindow)/10.)
minx = x[minxi]; maxx = x[maxxi]; 
Lx = float(rhoObj.Raw.header.globaldomain1)
Ly = float(rhoObj.Raw.header.globaldomain2)
binsize = [float(rhoObj.Raw.header.binsize1),
           float(rhoObj.Raw.header.binsize2)]
halfbinsize = [0.5*i for i in binsize]
xlims = range(minxi,maxxi)
Lmin = 0.15; Lmax = 0.5



f = plt.figure()
#plt.ion()
#plt.show()
for rec in range(startrec,endrec,1):
    print('Record ' + str(rec) + ' of ' + str(endrec))
    rho = rhoObj.read(rec,rec)
    vmd = DCDObj[rec][:]

    density = np.mean(rho[:,:,:,:,0],(2,3))
    edge, index = get_edge_cells(density,Lmin,Lmax)


    #Plot contourmap of edge position in x against z
    ax1 = plt.subplot2grid((3,1), (0,0))#, rowspan=2)
    ax2 = plt.subplot2grid((3,1), (1,0))
    ax3 = plt.subplot2grid((3,1), (2,0))
    #ax4 = plt.subplot2grid((2,2), (1,1))

    cm = ax1.pcolormesh(X[xlims,:],Z[xlims,:],np.mean(rho[xlims,:,:,:,0],(2,3)),cmap=plt.cm.RdYlBu_r,shading='gouraud')
    ax2.plot(X[edge],Z[edge],'o')
    molwindow = (((maxx-0.5*Lx)>vmd[:,0]) & ((minx-0.5*Lx)<vmd[:,0]))

    ax3.plot(vmd[molwindow,0],vmd[molwindow,1],'ko',ms=0.1)

    surfacemol = np.genfromtxt('../../../../MD_dCSE/src_code/fort.451')

    timestep = intialstep + (rec*plotfreq-1)*tplot
    print(timestep)
    rectime = surfacemol[:,0] == timestep
    #ax3.plot((surfacemol[rectime,1]-1)*binsize[0]-0.5*Lx,
    #         (surfacemol[rectime,2]-1)*binsize[1]-0.5*Ly,'bs')
    ax3.plot(surfacemol[rectime,7],surfacemol[rectime,8],'rx',ms=0.1)

    surfacenormal = np.genfromtxt('../../../../MD_dCSE/src_code/fort.452')
    ax3.quiver(surfacenormal[:,1],surfacenormal[:,2],surfacenormal[:,4],surfacenormal[:,5])

    surfacevar = np.genfromtxt('../../../../MD_dCSE/src_code/fort.1984')
    rectime = surfacevar[:,0] == timestep
    print(rectime)
    surfacevar = surfacevar[rectime,:]
    varmax = 10.
    varmin = 0.
    vartest = surfacevar[:,9]# np.mean(d[:,6:8],1)
    varfilter = ((varmax>vartest) & ((varmin<vartest)))
    try:
        cm=ax3.scatter(surfacevar[varfilter,1],surfacevar[varfilter,2],c=vartest[varfilter],cmap=plt.cm.RdYlBu_r)
    except:
        pass
#    for i in index[0]:
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
    #cbar_ax = f.add_axes([0.05, 0.1, 0.025, 0.8])
    #f.colorbar(cm, cax=cbar_ax)

    plt.show()
    #plt.pause(0.1)
    ax1.cla()
    ax2.cla()
    ax3.cla()

