import matplotlib.pyplot as plt
import numpy as np
import sys
import cPickle as pickle

sys.path.append('../../../')
import postproclib as ppl



#MD laminar and turbulent data
fdir = '/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/'
styles = ['k-', 'r-','b-']
startrec = 1000; endrec = 2500
recs = endrec - startrec + 1
f, ax = plt.subplots(2, 1)


MD_pp = ppl.MD_PostProc(fdir)	
vfield = MD_pp.plotlist['u']
vspectras = []

xfft = (vfield.power_spectrum(startrec=startrec,endrec=startrec,
                                       preavgaxes=(3),fftaxes=(0),
                                       postavgaxes=(2)))

zfft = (vfield.power_spectrum(startrec=startrec,endrec=startrec,
                                       preavgaxes=(3),fftaxes=(2),
                                       postavgaxes=(0)).transpose(1,0,2))
saved_xfft = np.empty((xfft.shape[0],3,recs))
saved_zfft = np.empty((zfft.shape[0],3,recs))
for rec in range(startrec,endrec):
    print(rec)
    try:
        xfft = (vfield.power_spectrum(startrec=rec,endrec=rec,
                                               preavgaxes=(3),fftaxes=(0),
                                               postavgaxes=(2)))

        zfft = (vfield.power_spectrum(startrec=rec,endrec=rec,
                                               preavgaxes=(3),fftaxes=(2),
                                               postavgaxes=(0)).transpose(1,0,2))
    except:
        pass

    locs = [5,xfft.shape[1]/4.,xfft.shape[1]/2.]
    for i, loc in enumerate(locs):
        saved_xfft[:,i,rec] = xfft[:,loc,0]
        saved_zfft[:,i,rec] = zfft[:,loc,0]

    if np.mod(rec,30) == 0:
        print('writing time history')
        pickle.dump( saved_xfft, open( "save_MDxfft.p", "wb" ) )
        pickle.dump( saved_zfft, open( "save_MDzfft.p", "wb" ) )

plt.savefig('MD_fft.png')

