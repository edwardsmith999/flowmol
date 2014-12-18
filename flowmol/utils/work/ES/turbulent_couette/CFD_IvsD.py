import matplotlib.pyplot as plt
import itertools
import numpy as np
import sys
import cPickle as pickle

sys.path.append('../../../')
import postproclib as ppl
from misclib import chaos_test

#fdirs = ['/home/es205/scratch/data_DNS/data_DNS_longrun/laminar/', '/home/es205/scratch/data_DNS/data_DNS_longrun/run0/','/home/es205/scratch/data_DNS/onehundredflowthrough_vs_MD/']
#startrecs = [0,0,0]
#endrecs   = [143,5000,3500]
#skip = 1

#Ihist = []; Dhist = []
#for i, fdir in enumerate(fdirs):
#    startrec = startrecs[i]
#    endrec   = endrecs[i]
#    recrange = endrec - startrec + 1

#    vField = ppl.Channelflow_vField(fdir)
#    dudrField   = ppl.Channelflow_strainField(fdir)
#    dissipField = ppl.Channelflow_dissipField(fdir)
#    volumes     = dudrField.Raw.get_binvolumes()

#    xyz = dudrField.grid
#    [X,Y,Z] = np.meshgrid(xyz[0],xyz[1],xyz[2],index='ij')

#    for rec in range(startrec,endrec,skip):
#        try:
#            u = vField.read(startrec=rec,endrec=rec)
#            dudr = dudrField.read(startrec=rec,endrec=rec)
#            dissip = dissipField.read(startrec=rec,endrec=rec)
#        except:
#            continue

#        Irec = ( (np.sum(   dudr[:, 1,:,0,1],(0,1))
#                 /np.sum(volumes[:, 1,:,0],(0,1)))
#                +(np.sum(   dudr[:,-2,:,0,1],(0,1))
#                 /np.sum(volumes[:,-2,:,0],(0,1))) )
#        Drec = (  np.sum( dissip[:,2:-3,:,0,0],(0,1,2))
#                 /np.sum(volumes[:,2:-3,:,0],(0,1,2)))

#        #plt.plot(xyz[1], np.mean( u[:,:,:,0,0],(0,2)))
#        #plt.plot(xyz[1], np.sum( dissip[:,:,:,0,0],(0,2))
#        #           /np.sum(volumes[:,:,:,0],(0,2)))
#        #plt.show()

#        Ihist.append(Irec)
#        Dhist.append(Drec)
#        print(fdir, 'Rec= ', rec, 'D= ', Drec , 'I= ', Irec)
#        #print(fdir, 'Rec= ', rec, 'D= ', Drec/8.4595006495855801e-07 , 'I= ', Irec/-1.8381574557365301e-06)

#pickle.dump( Ihist, open( "Ihist_CFD.p", "wb" ) )
#pickle.dump( Dhist, open( "Dhist_CFD.p", "wb" ) )

Ihist = pickle.load(open( "Ihist_CFD.p", "rb" ) )
Dhist = pickle.load(open( "Dhist_CFD.p", "rb" ) )

Ihist[5143] = Ihist[5142]
laminarend = 143
oneflowthroughend = 4999
plt.plot(Ihist[0:laminarend]/np.mean(Ihist[0:laminarend]),'bo',ms=1.,alpha=0.4)
plt.plot(Ihist[laminarend:oneflowthroughend]/np.mean(Ihist[0:laminarend]), 'r-',ms=1.,alpha=0.4)
plt.plot(Ihist[oneflowthroughend:]/np.mean(Ihist[0:laminarend]), 'k-',lw=3.,alpha=0.4)
plt.show()



plt.plot(Ihist[0:laminarend]/np.mean(Ihist[0:laminarend]), Dhist[0:laminarend]/np.mean(Dhist[0:laminarend]),'bo',ms=1.,alpha=0.4)
plt.plot(Ihist[laminarend:oneflowthroughend]/np.mean(Ihist[0:laminarend]), Dhist[laminarend:oneflowthroughend]/np.mean(Dhist[0:laminarend]),'r-',ms=1.,alpha=0.4)
plt.plot(Ihist[oneflowthroughend:]/np.mean(Ihist[0:laminarend]), Dhist[oneflowthroughend:]/np.mean(Dhist[0:laminarend]),'k-',lw=3.,alpha=0.4)
plt.show()

#print('I chaotic = ', chaos_test.zero_one_test(Ihist))
#print('D chaotic = ', chaos_test.zero_one_test(Dhist))

#Function differentiate noisy data (didn't seem to work...)
#dudy_noiseless = tvregdiff.TVRegDiff(u[:,0], 10, 5, plotflag=0)
#print(dudy_noiseless)
#plt.plot(y,dudy_noiseless[1:])
