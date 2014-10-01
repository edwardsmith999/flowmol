import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../../')
import postproclib as ppl


#CFD laminar and turbulent data
fdirs = ['/media/My Passport/Work/CFD_laminar/', '/media/My Passport/Work/CFD_minimal_channel/']
styles = ['r-', 'k-']
startrec = 300; endrec = 300
f, ax = plt.subplots(2, 1)

for i, fdir in enumerate(fdirs):

    CFD_pp = ppl.channelflow_PostProc(fdir)	
    vfield = CFD_pp.plotlist['channelflow Velocity']

    vspectras = []
    vspectras.append(vfield.power_spectrum(startrec=startrec,endrec=endrec,
                                           preavgaxes=(3),fftaxes=(0),
                                           postavgaxes=(2)))

    vspectras.append(vfield.power_spectrum(startrec=startrec,endrec=endrec,
                                           preavgaxes=(3),fftaxes=(2),
                                           postavgaxes=(0)).transpose(1,0,2))
    machineeps = 1.11e-16
    for num, vspectra in enumerate(vspectras):
        aliasloc = int((2./3.)*vspectra.shape[0])
        locs = [5,int(vspectra.shape[1]/4.),int(vspectra.shape[1]/2.)]
        #ax[num].plot(vspectra[:,locs[0],0],styles[i],alpha=0.8)
        #ax[num].plot(vspectra[:,locs[1],0],styles[i]+'-',alpha=0.8)
        ax[num].plot(vspectra[:aliasloc,locs[2],0]+machineeps,styles[i],alpha=0.8)

        ax[num].set_xlim([0.9, 50])
        ax[num].set_ylim([1e-17, 10])
        ax[num].set_xscale('log')
        ax[num].set_yscale('log')
        ax[num].set_xlabel('$k$')
        ax[num].set_ylabel('$E$')

plt.savefig('Channelflow_fft.png')

