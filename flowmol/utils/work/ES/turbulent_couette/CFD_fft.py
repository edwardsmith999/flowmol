import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../../')
import postproclib as ppl

#CFD laminar data
fdir = '/home/es205/codes/coupled/CFD_dCSE/src_code/results/'
fdir = './'
startrec = 3; endrec = 3
recs = [1,2,3]
styles = ['k-','k--','k:']
labels = ['$Re=4000$','$Re=1000$','$Re=400$']

#CFD turbulent data
#fdir = '/home/es205/scratch/data_DNS/onehundredflowthrough_vs_CFD/'
#startrec = 100; endrec = 150

CFD_pp = ppl.CFD_PostProc(fdir)	

vfield = CFD_pp.plotlist['CFD Velocity']

domain = [vfield.grid[0][-1],vfield.grid[2][-1]]
domainsize = [[domain[0],    domain[1]],
              [domain[0]/4., domain[1]/4.],
              [domain[0]/10.,domain[1]/10.]]
volume = np.empty(3)
#volume[0]=vfield.grid[0][-1]*vfield.grid[1][-1]*vfield.grid[2][-1]
volume[0]=1.0
volume[1]=volume[0]/64.
volume[2]=volume[0]/1000.
f, ax = plt.subplots(1, 2)
for i, rec in enumerate(recs):
    vcont = vfield.read(startrec=rec,endrec=rec)
    x,y,z=vfield.grid
    X,Z = np.meshgrid(x,z)
    print(rec)

#    figure = plt.figure()
#    axtemp = figure.add_subplot(111)
#    cmap = matplotlib.cm.RdYlBu_r
#    colormesh = axtemp.pcolormesh(X, Z, vcont[:,128,:,0,0], cmap=cmap)
#    cbar = figure.colorbar(colormesh)
#    figure.show()

    vspectras = []
    vspectras.append(vfield.power_spectrum(startrec=rec,endrec=rec,
                                           preavgaxes=(3),fftaxes=(0),
                                           postavgaxes=(2)))

    vspectras.append(vfield.power_spectrum(startrec=rec,endrec=rec,
                                           preavgaxes=(3),fftaxes=(2),
                                           postavgaxes=(0)).transpose(1,0,2))

    grid = vfield.grid[0::2]
    for num, vspectra in enumerate(vspectras):
        locs = [5,int(vspectra.shape[1]/4.),int(vspectra.shape[1]/2.)]
        #k = 2.*np.pi*np.arange(vspectra.shape[0])/domainsize[i][num]
        k = 2.*np.pi*np.arange(vspectra.shape[0])
    #    ax[num].plot(vspectra[:,5,0],'o-',alpha=0.4)
    #    ax[num].plot(vspectra[:,45,0],'o-',alpha=0.4)
        ax[num].plot(k[:-1],vspectra[:-1,locs[2],0],styles[i],alpha=0.8,label=labels[i])
#        ax[num].plot(k[:-1],vspectra[:-1,locs[2],0]*volume[i],styles[i],alpha=0.8,label=labels[i])
        
        #ax[num].set_xlim([0.9, 50])
        #ax[num].set_ylim([1e-16, 1e2])
        ax[num].set_xscale('log')
        ax[num].set_yscale('log')

        ax[num].set_xlabel('$k$')
        ax[num].set_ylabel('$E$')
        plt.legend(loc='best')

plt.show()

