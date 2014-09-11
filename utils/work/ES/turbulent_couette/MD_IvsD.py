import matplotlib.pyplot as plt
import itertools
import numpy as np
import sys

sys.path.append('../../../')
import postproclib as ppl
from postproclib.mdpostproc import MD_PostProc
from misclib import tvregdiff

fdirs = ['/media/My Passport/Work/MD_laminar/Re400_transition/iter350000_to_1050000/results', '/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins64x256x64/','/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/']
startrecs = [0,0,960]
endrecs   = [185,980,3500]
skip = 300

I = []; D = []
for i, fdir in enumerate(fdirs):
    startrec = startrecs[i]
    endrec   = endrecs[i]
    recrange = endrec - startrec + 1

    dudrField   = ppl.MD_strainField(fdir)
    dissipField = ppl.MD_dissipField(fdir)
    volumes     = dudrField.Raw.get_binvolumes()

    xyz = dudrField.grid
    [X,Y,Z] = np.meshgrid(xyz[0],xyz[1],xyz[2],index='ij')

    for rec in range(startrec,endrec,skip):
        try:
            dudr = dudrField.read(startrec=rec,endrec=rec)
            dissip = dissipField.read(startrec=rec,endrec=rec)
        except:
            continue

        Irec = ( (np.sum(   dudr[:, 1,:,0,1],(0,1))
                 /np.sum(volumes[:, 1,:,0],(0,1)))
                +(np.sum(   dudr[:,-1,:,0,1],(0,1))
                 /np.sum(volumes[:,-1,:,0],(0,1))) )
        Drec = (  np.sum( dissip[:,2:-2,:,0,0],(0,1,2))
                 /np.sum(volumes[:,2:-2,:,0],(0,1,2)))
        I.append(Irec)
        D.append(Drec)
        print('Rec= ', rec, 'D= ', Drec , 'I= ', Irec)

plt.plot(D,I,'-o',alpha=0.7)
plt.show()

#Function differentiate noisy data (didn't seem to work...)
#dudy_noiseless = tvregdiff.TVRegDiff(u[:,0], 10, 5, plotflag=0)
#print(dudy_noiseless)
#plt.plot(y,dudy_noiseless[1:])
