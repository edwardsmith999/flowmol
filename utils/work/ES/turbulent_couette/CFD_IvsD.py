import matplotlib.pyplot as plt
import itertools
import numpy as np
import sys
import cPickle as pickle

sys.path.append('../../../')
import postproclib as ppl
from misclib import chaos_test

fdirs = ['/home/es205/scratch/data_DNS/data_DNS_longrun/laminar/', 
         '/home/es205/scratch/data_DNS/data_DNS_longrun/run0/',
         '/home/es205/scratch/data_DNS/onehundredflowthrough_vs_MD/']
startrecs = [0,0,0]
endrecs   = [143,5000,3500]
flowtype = ['laminar','turbulent','turbulent']
skip = 100


CFout_laminar = np.genfromtxt('./Channelflowout_IvsD_laminar')
CFout_longrun = np.genfromtxt('./Channelflowout_IvsD_run0')
CFout_oneflowthu = np.genfromtxt('./Channelflowout_IvsD_onehundredflowthu')

#plt.plot(CFout_longrun[:,2],CFout_longrun[:,5]/np.mean(CFout_laminar[:,5]))
#plt.plot(CFout_oneflowthu[:,2],CFout_oneflowthu[:,5]/np.mean(CFout_laminar[:,5]))
#plt.show()

Ihist = []; Dhist = []
wallcell = 1
for i, fdir in enumerate(fdirs):
    startrec = startrecs[i]
    endrec   = endrecs[i]
    recrange = endrec - startrec + 1

    vField = ppl.Channelflow_vField(fdir)
    dudrField   = ppl.Channelflow_strainField(fdir)
    dissipField = ppl.Channelflow_dissipField(fdir)
    volumes     = dudrField.Raw.get_binvolumes()

    xyz = dudrField.grid
    [X,Y,Z] = np.meshgrid(xyz[0],xyz[1],xyz[2],index='ij')

    for rec in range(startrec,endrec,skip):
        try:
            u = vField.read(startrec=rec,endrec=rec)
            dudr = dudrField.read(startrec=rec,endrec=rec)
            dissip = dissipField.read(startrec=rec,endrec=rec)
        except:
            raise
            continue

        y = dudrField.grid[1]

#        print(wallcell,np.mean(   dudr[:,wallcell,:,0,1],(0,1)),np.mean(   dudr[:,-wallcell,:,0,1],(0,1)))

#        plt.plot(y, np.mean(   dudr[:, :,:,0,1],(0,2)),'b-')
#        plt.plot(y[wallcell],np.mean(   dudr[:, wallcell,:,0,1],(0,1)),'bo')
#        plt.plot(y[-wallcell],np.mean(   dudr[:,-(wallcell+1),:,0,1],(0,1)),'bo')
#        plt.plot(y,np.mean( dissip[:, :,:,0,0],(0,2)),'r-')
#        plt.plot(y[(wallcell+1):-(wallcell+2)],np.mean( dissip[:,(wallcell+1):-(wallcell+2),:,0,0],(0,2)),'ro')
#        plt.show()

#        Irec = ( (np.sum(   dudr[:, wallcell,:,0,1],(0,1))
#                 /np.sum(volumes[:, wallcell,:,0],(0,1)))
#                +(np.sum(   dudr[:,-(wallcell+1),:,0,1],(0,1))
#                 /np.sum(volumes[:,-(wallcell+1),:,0],(0,1))) )
#        Drec = (  np.sum( dissip[:,(wallcell+1):-(wallcell+2),:,0,0],(0,1,2))
#                 /np.sum(volumes[:,(wallcell+1):-(wallcell+2),:,0],(0,1,2)))

        Irec =0.5*(  np.mean(dudr[:,  wallcell,:,0,1],(0,1))
                   + np.mean(dudr[:,-(wallcell+1),:,0,1],(0,1)))
        Drec =       np.mean(dissip[:,(wallcell+1):-(wallcell+2),:,0,0],(0,1,2))

        #plt.plot(xyz[1], np.mean( u[:,:,:,0,0],(0,2)))
        #plt.plot(xyz[1], np.sum( dissip[:,:,:,0,0],(0,2))
        #           /np.sum(volumes[:,:,:,0],(0,2)))
        #plt.show()

        Ihist.append(Irec)
        Dhist.append(Drec)
        print(fdir, 'Rec= ', rec, 'D= ', Drec , 'I= ', Irec)

    if flowtype[i] is 'laminar':
       Dlam = np.mean(Dhist) 
       Ilam = np.mean(Ihist) 
       print('LAMINAR D = ', Dlam, 'I', Ilam)
        #print(fdir, 'Rec= ', rec, 'D= ', Drec/8.4595006495855801e-07 , 'I= ', Irec/-1.8381574557365301e-06)

#pickle.dump( Ihist, open( "Ihist_CFD.p", "wb" ) )
#pickle.dump( Dhist, open( "Dhist_CFD.p", "wb" ) )

#Ihist = pickle.load(open( "Ihist_CFD.p", "rb" ) )
#Dhist = pickle.load(open( "Dhist_CFD.p", "rb" ) )

Ihist[5143] = Ihist[5142]
laminarend = 143
oneflowthroughend = 5150# 4999
#plt.plot(Ihist[0:laminarend]/np.mean(Ihist[0:laminarend]),'bo',ms=1.,alpha=0.4)
#plt.plot(Ihist[laminarend:oneflowthroughend]/np.mean(Ihist[0:laminarend]), 'r-',ms=1.,alpha=0.4)
#plt.plot(Ihist[oneflowthroughend:]/np.mean(Ihist[0:laminarend]), 'k-',lw=3.,alpha=0.4)
#plt.show()

print(len(Ihist),CFout_oneflowthu.shape,CFout_laminar.shape,CFout_longrun.shape)

#plt.plot(Ihist[oneflowthroughend:]/np.mean(Ihist[0:laminarend]),'b')
#plt.plot(3.*(1.-CFout_oneflowthu[:,2])/np.mean(1.-CFout_oneflowthu[:,2]),'r')
#plt.show()
fact = 0.5
plt.plot(Ihist[0:laminarend]/np.mean(Ihist[0:laminarend]), fact*np.array(Dhist[0:laminarend]/np.mean(Dhist[0:laminarend])),'bo',ms=1.,alpha=0.4)
plt.plot(Ihist[laminarend:oneflowthroughend]/np.mean(Ihist[0:laminarend]), fact*np.array(Dhist[laminarend:oneflowthroughend]/np.mean(Dhist[0:laminarend])),'r-',ms=1.,alpha=0.4)
plt.plot(Ihist[oneflowthroughend:]/np.mean(Ihist[0:laminarend]), fact*np.array(Dhist[oneflowthroughend:]/np.mean(Dhist[0:laminarend])),'k-',lw=3.,alpha=0.4)
start = 10
plt.plot(3.*(1.-CFout_oneflowthu[start:,2])/np.mean(1.-CFout_oneflowthu[start:,2]),CFout_oneflowthu[start:,5]/np.mean(CFout_laminar[:,5]),'b')
plt.plot(3.*(1.-CFout_longrun[start:,2])/np.mean(1.-CFout_longrun[start:,2]),CFout_longrun[start:,5]/np.mean(CFout_laminar[:,5]),'k')
plt.plot(np.linspace(0.0,5.0,100),np.linspace(0.0,5.0,100),'--',color=[.7,.7,.7])
plt.show()

#print('I chaotic = ', chaos_test.zero_one_test(Ihist))
#print('D chaotic = ', chaos_test.zero_one_test(Dhist))

#Function differentiate noisy data (didn't seem to work...)
#dudy_noiseless = tvregdiff.TVRegDiff(u[:,0], 10, 5, plotflag=0)
#print(dudy_noiseless)
#plt.plot(y,dudy_noiseless[1:])
