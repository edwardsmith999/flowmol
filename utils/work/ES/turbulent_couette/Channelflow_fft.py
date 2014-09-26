import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../../')
import postproclib as ppl


#CFD laminar and turbulent data
fdirs = ['/media/My Passport/Work/CFD_laminar/', '/media/My Passport/Work/CFD_minimal_channel/']
styles = ['-o', '-']
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

    for num, vspectra in enumerate(vspectras):
        locs = [5,vspectra.shape[1]/4.,vspectra.shape[1]/2.]
        ax[num].plot(vspectra[:,5,0],styles[i]+'k',alpha=0.4)
        ax[num].plot(vspectra[:,45,0],styles[i]+'r',alpha=0.4)
        ax[num].plot(vspectra[:,99,0],styles[i]+'b',alpha=0.4)

        #ax[num].set_ylim([1e-16, 1e2])
        ax[num].set_xscale('log')
        ax[num].set_yscale('log')

plt.savefig('CFD_fft.png')

