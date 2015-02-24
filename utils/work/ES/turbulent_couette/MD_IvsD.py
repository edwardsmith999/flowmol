import matplotlib.pyplot as plt
import itertools
import numpy as np
import sys
import cPickle as pickle

sys.path.append('/home/es205/codes/coupled/utils/')
import postproclib as ppl
from postproclib.mdpostproc import MD_PostProc
from postproclib.pplexceptions import DataNotAvailable
#from misclib import tvregdiff

fdirs = ['/home/es205/scratch/Re400_transition/iter350000_to_1050000/results', 
         '/home/es205/scratch/Re43p5/',
         '/home/es205/scratch/Re43p5/iter_101000_to_201000/results/',
         '/home/es205/scratch/Re43p5/iter_201000_to_301000/results/']#,
         #'/home/es205/scratch/Re96/up2iter100000/',
         #'/home/es205/scratch/Re96/iter100000_2_300000/',
         #'/home/es205/scratch/Re400/COMBINED_iter0_to_5600000/bins64x256x64/']
         #'/home/es205/scratch/Re400/COMBINED_iter0_to_5600000/bins84x198x50/']
startrecs = [0,  
             0,
             0,
             0,
             0,
             80,
             980 ]
endrecs   = [185,
             390,
             390,
             390,
             1190,
             970,
             3500]
sizeratio = [1., 1., 1., 1., 1., 1.]
flowtype = ['laminar','turbulent','turbulent','turbulent','turbulent','turbulent']
skip = 1; avetime = 100

startrecs = [i + int(0.5*avetime) for i in startrecs]
endrecs   = [i - int(0.5*avetime) for i in endrecs]

I = []; D = []
for i, fdir in enumerate(fdirs):
    print(i,fdir)
    startrec = startrecs[i]
    endrec   = endrecs[i]
    recrange = endrec - startrec + 1

    dudrField   = ppl.MD_strainField(fdir)
    dissipField = ppl.MD_dissipField(fdir)
    volumes     = dudrField.Raw.get_gridvolumes()
    y = dudrField.grid[1]
    #TEST OF HIGH PASS FILTER CODE
    #dissip = dissipField.read(startrec=int(0.5*recrange),endrec=int(0.5*recrange)+1,highpassfilter=0.9)

   # xyz = dudrField.grid
    #[X,Y,Z] = np.meshgrid(xyz[0],xyz[1],xyz[2],index='ij')
    wallcell = 4
    for rec in range(startrec,endrec,skip):
        try:
            # AS OF 27/12/14 DUDR AND CONSEQUENTLY DISSIP HAS A BUG WHERE
            # TENSOR IS NOT SYMMETRICAL (PROBLEM WITH FIELD GRAD ROUTINE I THINK!!)
            dudr   = dudrField.read(  startrec=rec-int(0.5*avetime),
                                        endrec=rec+int(0.5*avetime),
                                     missingrec='skip')
            
            dissip = dissipField.read(startrec=rec-int(0.5*avetime),
                                        endrec=rec+int(0.5*avetime),
                                      missingrec='skip')
        except DataNotAvailable:
            continue

#        plt.plot(y, np.mean(   dudr[:, :,:,0,:],(0,2)),'b-')
#        plt.show()

#        print(wallcell,np.mean(   dudr[:,wallcell,:,0,1],(0,1)),np.mean(   dudr[:,-wallcell,:,0,1],(0,1)))
#        plt.plot(y, np.mean(   dudr[:, :,:,0,1],(0,2)),'b-')
#        plt.plot(y, np.mean(   dudr[:, :,:,0,3],(0,2)),'k--')
#        plt.plot(y[wallcell],np.mean(   dudr[:, wallcell,:,0,1],(0,1)),'bo')
#        plt.plot(y[-wallcell],np.mean(   dudr[:,-wallcell,:,0,1],(0,1)),'bo')
#        plt.plot(y,np.mean( dissip[:, :,:,0,0],(0,2)),'r-')
#        plt.plot(y[(wallcell+1):-(wallcell+1)],np.mean( dissip[:,(wallcell+1):-(wallcell+1),:,0,0],(0,2)),'ro')
#        plt.show()

#        Irec =0.5*( (np.sum(   dudr[:, wallcell,:,0,1],(0,1))
#                    /np.sum(volumes[:, wallcell,:,0],(0,1)))
#                   +(np.sum(   dudr[:,-wallcell,:,0,1],(0,1))
#                    /np.sum(volumes[:,-wallcell,:,0],(0,1))) )
#        Drec = (  np.sum( dissip[:,wallcell+1:-(wallcell+1),:,0,0],(0,1,2))
#                 /np.sum(volumes[:,wallcell+1:-(wallcell+1),:,0],(0,1,2)))

        Irec =0.5*(  np.mean(dudr[:, wallcell,:,0,1],(0,1))
                   + np.mean(dudr[:,-wallcell,:,0,1],(0,1)))/sizeratio[i]
        #Drec = np.mean(np.power(dudr[:,(wallcell+1):-(wallcell+1),:,0,1],2.))
        Drec =       np.mean(dissip[:,(wallcell+1):-(wallcell+1),:,0,0],(0,1,2))/np.power(sizeratio[i],2)
        I.append(Irec)
        D.append(Drec)
        print('Rec= ', rec, 'D= ', Drec , 'I= ', Irec)

        if( np.mod(rec,100) == 0):
            print('100 iter dump at ', rec)
            pickle.dump( I, open( "Ihist_MD.p", "wb" ) )
            pickle.dump( D, open( "Dhist_MD.p", "wb" ) )


    print(flowtype[i],flowtype[i] is 'laminar',np.mean(D))
    if flowtype[i] is 'laminar':
       Dlam = np.mean(D) 
       Ilam = np.mean(I) 
       print('LAMINAR D = ', Dlam, 'I', Ilam)

pickle.dump( I, open( "Ihist_MD.p", "wb" ) )
pickle.dump( D, open( "Dhist_MD.p", "wb" ) )

plt.plot(I/Ilam,D/Dlam,'-o',alpha=0.7)
plt.show()

#Function differentiate noisy data (didn't seem to work...)
#dudy_noiseless = tvregdiff.TVRegDiff(u[:,0], 10, 5, plotflag=0)
#print(dudy_noiseless)
#plt.plot(y,dudy_noiseless[1:])
