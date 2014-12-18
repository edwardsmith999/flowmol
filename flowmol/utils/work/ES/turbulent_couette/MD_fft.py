import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import sys
import os

sys.path.append('../../../')
import postproclib as ppl

#Load the colormap -- full range of colors can be got with
# colourscheme(0) to colourscheme(256)
# Good values for lines:
#  Dark blue = 0
#  Light blue = 64
#  Yellow = 150
#  Light Red = 175
#  Red = 256
colourscheme=plt.cm.RdYlBu_r
dashes = [(None,None),[5,5],[1,3],[5,3,1,3],[5,3,1,3,1,3]]
colors = [colourscheme(i) for i in range(0,257,64)]
colors.insert(0,(0.,0.,0.,0.))

#MD laminar and turbulent data
#fdirs = [ '/media/My Passport/Work/MD_laminar/Re400_transition/iter350000_to_1050000/results/',
#          '/media/My Passport/Work/MD_laminar/Re400_transition/iter350000_to_1050000/results/',
#          '/media/My Passport/Work/MD_laminar/Re400_transition/iter350000_to_1050000/results/',
#          '/media/My Passport/Work/MD_laminar/Re400_transition/iter350000_to_1050000/results/',
#          '/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/',
#          '/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/',
#          '/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/',
#          '/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/',
#          '/media/My Passport/Work/MD_turbulence/COMBINED_iter0_to_5600000/bins84x198x50/',
#          '/media/My Passport/Work/MD_turbulence/iter5000000_to_5600000/']
fdirs = ['/home/es205/scratch/Re400_transition/iter350000_to_1050000/results/',
          '/home/es205/scratch/Re400_transition/iter350000_to_1050000/results/',
          '/home/es205/scratch/Re400_transition/iter350000_to_1050000/results/',
          '/home/es205/scratch/Re400_transition/iter350000_to_1050000/results/',
          '/home/es205/scratch/Re400/iter3734274_to_4384274/',
          '/home/es205/scratch/Re400/COMBINED_iter0_to_5600000/bins84x198x50/',
          '/home/es205/scratch/Re400/COMBINED_iter0_to_5600000/bins84x198x50/',
          '/home/es205/scratch/Re400/COMBINED_iter0_to_5600000/bins84x198x50/',
          '/home/es205/scratch/Re400/COMBINED_iter0_to_5600000/bins84x198x50/',
          '/home/es205/scratch/Re400/COMBINED_iter0_to_5600000/bins84x198x50/']
styles = [{'color':colors[0],'dashes':dashes[0]},
          {'color':colors[0],'dashes':dashes[1]},
          {'color':colors[0],'dashes':dashes[2]},
          {'color':colors[0],'dashes':dashes[3]},
          {'color':colors[1],'dashes':dashes[0]},
          {'color':colors[2],'dashes':dashes[0]},
          {'color':colors[3],'dashes':dashes[0]},
          {'color':colors[4],'dashes':dashes[0]},
          {'color':colors[5],'dashes':dashes[0]},
          {'color':colors[1],'dashes':dashes[1]}]
for i in styles:
    i['alpha'] = 0.8

startrecs = [20,20,20,20,166,2500,2500,2500,2500,2500]
endrecs = [20,30,70,180,166,2500,2510,2550,2660,3499]
datatype = ['bins','bins','bins','bins','snap',
            'bins','bins','bins','bins','bins']
labels = []
for i in range(len(startrecs)):
    labels.append('')

    if fdirs[i] is '/home/es205/scratch/Re400_transition/iter350000_to_1050000/results/':
        labels[i] += 'Laminar '
    elif fdirs[i] is '/home/es205/scratch/Re400/COMBINED_iter0_to_5600000/bins84x198x50/':
        labels[i] += 'Turbulent '
    elif fdirs[i] is '/home/es205/scratch/Re400/iter3734274_to_4384274/':
        labels[i] += 'Turbulent '
    else:
        print(fdirs[i])
        quit('Error -- unknown if Laminar or Turbulent')

    if datatype[i] is 'bins':
        labels[i] += "avrgd timesteps = " + str((endrecs[i]-startrecs[i]+1)*64) 
    elif datatype[i] is 'snap':
        labels[i] += "avrgd timesteps = " + str((endrecs[i]-startrecs[i]+1)) 
    else:
        quit('Error -- unknown bins or snap datatype')

f, ax = plt.subplots(3, 1)
component = 0

lines = []
for i, fdir in enumerate(fdirs):
    print(fdir)
    try:
        vfield = ppl.MD_vField(fdir,rectype=datatype[i])
        #vfield = ppl.MD_rhouuField(fdir)
        startrec = startrecs[i]; endrec=endrecs[i]
        vdata = vfield.averaged_data(startrec=startrec,endrec=endrec,
                                     avgaxes=(3),missingrec='skip')

        print(fdir,vdata.shape)

#            for rec in range(startrec,endrec):
#                try:
#                    vdata = vfield.read(startrec=rec,endrec=rec)
#                except

        if datatype[i] is 'snap':
            vdata = vdata*1217

        #ax[2].pcolormesh(vfield.grid[2], vfield.grid[0], 
        #                 np.mean(vdata[:,vdata.shape[1]/2.,:,:,component],(2)), 
        #                 cmap=colourscheme)
        ax[2].plot(vfield.grid[1],np.mean(vdata[:,:,:,component],(0,2)))

        vspectras = []
        vspectras.append(vfield.power_spectrum(data=vdata,
                                               preavgaxes=(),fftaxes=(0),
                                               postavgaxes=(2)))

        vspectras.append(vfield.power_spectrum(data=vdata,
                                               preavgaxes=(),fftaxes=(2),
                                               postavgaxes=(0)).transpose(1,0,2))

        for num, vspectra in enumerate(vspectras):
            locs = [5,vspectra.shape[1]/4.,vspectra.shape[1]/2.]
            #ax[num].plot(vspectra[:,locs[0],0],styles[i],alpha=0.8)
            #ax[num].plot(vspectra[:,locs[1],0],styles[i]+'-',alpha=0.8)
            #ax[num].plot(vspectra[:,locs[2],0],styles[i].replace('-',':'),alpha=0.8)

            #line = ax[num].plot(range(1,vspectra.shape[0]+1),vspectra[:,locs[2],1],**style)

            line = ax[num].plot(range(1,vspectra.shape[0]+1),vspectra[:,locs[2],1],alpha=0.8,lw=2.,label=labels[i])
            line[0].set_dashes(styles[i]['dashes'])
            line[0].set_color(styles[i]['color'])
            lines.append(line)
#            ax[num].set_xlim([0.9, 50])
#            ax[num].set_ylim([1e-5, 1e0])
            #ax[num].set_ylim([1e-16, 1e2])
            ax[num].set_xscale('log')
            ax[num].set_yscale('log')

            ax[num].set_xlabel('$k$')
            ax[num].set_ylabel('$E$')

    except:
        raise    


#Set default figure properties
ax[0].legend(loc='best')
ax[1].legend(loc='best')
#plt.rc('font', size=24)
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)

figname = 'MD_fft_' + str(vfield.labels[component])
plt.savefig(figname+'pdf', bbox_inches='tight')
os.system('convert ' + figname + '.pdf ' +  figname + '.png')

