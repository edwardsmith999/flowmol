import matplotlib.pyplot as plt
import itertools
import numpy as np
import sys

sys.path.append('../../../')
import postproclib as ppl
from misclib import chaos_test

fdirs = ['/media/My Passport/Work/CFD_laminar/','/media/My Passport/Work/CFD_minimal_channel/']
startrecs = [0,0]
endrecs   = [360,3500]
skip = 10

I = []; D = []
for i, fdir in enumerate(fdirs):
    startrec = startrecs[i]
    endrec   = endrecs[i]
    recrange = endrec - startrec + 1

    dudrField   = ppl.Channelflow_strainField(fdir)
    dissipField = ppl.Channelflow_dissipField(fdir)
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
        print('Rec= ', rec, 'D= ', Drec/8.4595006495855801e-07 , 'I= ', Irec/-1.8381574557365301e-06)

plt.plot(D,I,'-o',alpha=0.7)
plt.show()

chaos_test(I)
chaos_test(D)

#Function differentiate noisy data (didn't seem to work...)
#dudy_noiseless = tvregdiff.TVRegDiff(u[:,0], 10, 5, plotflag=0)
#print(dudy_noiseless)
#plt.plot(y,dudy_noiseless[1:])
