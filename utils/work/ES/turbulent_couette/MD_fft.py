import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../../')
import postproclib as ppl



#MD laminar and turbulent data
fdirs = [ '/media/My Passport/Work/MD_laminar/Re400_transition/iter350000_to_1050000/results/',
          '/media/My Passport/Work/MD_laminar/Re400_transition/iter350000_to_1050000/results/',
          '/media/My Passport/Work/MD_laminar/Re400_transition/iter350000_to_1050000/results/',
#          '/media/My Passport/Work/MD_turbulence/iter4384274_to_5000000/',
          '/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/',
          '/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/',
          '/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/']#,
#          '/media/My Passport/Work/MD_turbulence/iter5000000_to_5600000/']#,
#          ]
styles = ['r-','r--','r:', 'k-','k--','k:']#,'k-.']
startrecs = [50,50,50,2500,2500,2500]#,100]
endrecs = [50,60,100,2500,2510,2550]#,100]
datatype = ['u','u','u','u','u','u']#,'u_snap']
f, ax = plt.subplots(2, 1)

for i, fdir in enumerate(fdirs):
    print(fdir)
    try:
        MD_pp = ppl.MD_PostProc(fdir)	
        vfield = MD_pp.plotlist[datatype[i]]
        startrec = startrecs[i]; endrec=endrecs[i]
        vspectras = []
        vspectras.append(vfield.power_spectrum(startrec=startrec,endrec=endrec,
                                               preavgaxes=(3),fftaxes=(0),
                                               postavgaxes=(2)))

        vspectras.append(vfield.power_spectrum(startrec=startrec,endrec=endrec,
                                               preavgaxes=(3),fftaxes=(2),
                                               postavgaxes=(0)).transpose(1,0,2))

        for num, vspectra in enumerate(vspectras):
            locs = [5,vspectra.shape[1]/4.,vspectra.shape[1]/2.]
            #ax[num].plot(vspectra[:,locs[0],0],styles[i],alpha=0.8)
            #ax[num].plot(vspectra[:,locs[1],0],styles[i]+'-',alpha=0.8)
            #ax[num].plot(vspectra[:,locs[2],0],styles[i].replace('-',':'),alpha=0.8)
            if i == 7:
                vspectra = vspectra*73000
            ax[num].plot(vspectra[:,locs[2],0],styles[i],alpha=0.8)

            ax[num].set_xlim([0.9, 50])
            ax[num].set_ylim([1e-5, 1e0])
            #ax[num].set_ylim([1e-16, 1e2])
            ax[num].set_xscale('log')
            ax[num].set_yscale('log')

            ax[num].set_xlabel('$k$')
            ax[num].set_ylabel('$E$')

    except:
        pass    

plt.savefig('MD_fft.png')

