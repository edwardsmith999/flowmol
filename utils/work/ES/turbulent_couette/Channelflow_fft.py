import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import sys

sys.path.append('../../../')
import postproclib as ppl

#Load the colormap -- full range of colors can be got with
# colourscheme(0) to colourscheme(256)
# Good values for lines:
#  Dark blue = 0
#  Light blue = 64
#  Yellow = 150
#  Light Red = 175
#  Red = 256
colourscheme=plt.cm.RdYlBu_r
dashes = [(None,None),[5,5],[1,3],[5,3,1,3],[5,3,1,3,1,3]]
colors = [colourscheme(i) for i in range(0,257,64)]
colors.insert(0,(0.,0.,0.,0.))


#CFD laminar and turbulent data
fdirs = ['/media/My Passport/Work/CFD_laminar/', '/media/My Passport/Work/CFD_minimal_channel/']
#fdirs = ['/home/es205/scratch/data_DNS/data_DNS_longrun/laminar/', '/home/es205/scratch/data_DNS/onehundredflowthrough_vs_MD/']
styles = [colors[2], colors[4]]
startrec = 100; endrec = 100
f, ax = plt.subplots(3, 1)
rho = 0.3
component = 0

for i, fdir in enumerate(fdirs):

    vfield = ppl.Channelflow_vField(fdir)
    vdata = vfield.read(startrec=startrec,endrec=endrec)
    ax[2].plot(vfield.grid[1],np.mean(vdata[:,:,:,:,component],(0,2,3)),'o-')

    #vfield = ppl.Channelflow_uuField(fdir)
    #vdata = rho*vdata

    vspectras = []
    vspectras.append(vfield.power_spectrum(data=vdata,
                                           preavgaxes=(3),fftaxes=(0),
                                           postavgaxes=(2)))

    vspectras.append(vfield.power_spectrum(data=vdata,
                                           preavgaxes=(3),fftaxes=(2),
                                           postavgaxes=(0)).transpose(1,0,2))
    machineeps = 1.11e-16
    for num, vspectra in enumerate(vspectras):
        aliasloc = vspectra.shape[0]-2 #int((2./3.)*vspectra.shape[0])
        locs = [5,int(vspectra.shape[1]/4.),int(vspectra.shape[1]/2.)]
#        ax[num].plot(range(1,vspectra[:aliasloc,:,:].shape[0]+1),
#                             vspectra[:aliasloc,locs[0],0]+machineeps,
#                             styles[i],lw=2.,alpha=0.8)
#        ax[num].plot(range(1,vspectra[:aliasloc,:,:].shape[0]+1),
#                             vspectra[:aliasloc,locs[1],0]+machineeps,
#                             styles[i]+':',lw=2.,alpha=0.8)
        ax[num].plot(range(1,vspectra[:aliasloc,:,:].shape[0]+1),
                             vspectra[:aliasloc,locs[2],0]+machineeps,'-o',
                             color=styles[i],lw=2.,alpha=0.8)


        #ax[num].set_xlim([0.9, 50])
        #ax[num].set_ylim([1e-17, 10])
        ax[num].set_xscale('log')
        ax[num].set_yscale('log')
        ax[num].set_xlabel('$k$')
        ax[num].set_ylabel('$E$')

plt.show()
#plt.savefig('Channelflow_fft_' + str(vfield.labels[component]) + '.png')

