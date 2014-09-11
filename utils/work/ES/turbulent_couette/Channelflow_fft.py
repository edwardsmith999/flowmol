import postproclib as ppl
from postproclib.channelflowpostproc import channelflow_PostProc
import matplotlib.pyplot as plt
import itertools

#CFD laminar data
fdir = '/media/My Passport/Work/CFD_minimal_channel/'
fdir = '/media/My Passport/Work/CFD_minimal_channel/'

startrec = 1050; endrec = 1051
#CFD turbulent data
#fdir = '/home/es205/scratch/data_DNS/onehundredflowthrough_vs_CFD/'
#startrec = 100; endrec = 150

CFD_pp = channelflow_PostProc(fdir)	

vfield = CFD_pp.plotlist['channelflow Velocity']

vspectras = []
vspectras.append(vfield.power_spectrum(startrec=startrec,endrec=endrec,
                                       preavgaxes=(3),fftaxes=(0),
                                       postavgaxes=(2)))

vspectras.append(vfield.power_spectrum(startrec=startrec,endrec=endrec,
                                       preavgaxes=(3),fftaxes=(2),
                                       postavgaxes=(0)).transpose(1,0,2))

f, ax = plt.subplots(2, 1)
for num, vspectra in enumerate(vspectras):
    ax[num].plot(vspectra[:,5,0],'o-',alpha=0.4)
    ax[num].plot(vspectra[:,45,0],'o-',alpha=0.4)
    ax[num].plot(vspectra[:,99,0],'o-',alpha=0.4)

    #ax[num].set_ylim([1e-5, 1e2])
    ax[num].set_xscale('log')
    ax[num].set_yscale('log')

plt.show()

