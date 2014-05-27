import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

from MDFields import *
from HeaderData import *
from MD_PostProc import MD_PostProc


#def init(field):
#    initpos = int(len(field.grid[1])/2.)
#    initrec = int(field.maxrec/2.)
#    ax1, ax2, data = field.contour(axes=naxes,startrec=initrec,endrec=initrec,binlimits=[None,(initpos,initpos),None])
#    colormesh = ax.pcolormesh(ax1, ax2, data[:,1:,component],cmap=cmap,shading='gourand')

#    return initpos, initrec, ax1, ax2, data, colormesh

def update(val):
    pos = int(np.rint(sbin.val))
    rec = int(np.rint(srec.val))
    ax1, ax2, data = field.contour(axes=naxes,startrec=rec,endrec=rec,binlimits=[None,(pos,pos),None])
    colormesh.set_array(data[:,1:,component].ravel())
    ax.figure.canvas.draw_idle()
    print(pos,rec,field,colormesh)

def plotfield(label):
    print(label)
    field = fielddict.plotlist[label]
    initpos = int(len(field.grid[1])/2.)
    initrec = int(field.maxrec/2.)
    ax1, ax2, data = field.contour(axes=naxes,startrec=initrec,endrec=initrec,binlimits=[None,(initpos,initpos),None])
    #colormesh = ax.pcolormesh(ax1, ax2, data[:,1:,component],cmap=cmap)
    colormesh.set_array(data[:,1:,component].ravel())

    #initpos, initrec = init(field)[0:2]
    ax.figure.canvas.draw_idle()
    print(field,label)


component = 0
cmap = plt.cm.RdYlBu_r
naxes = (0,2)

#Setup array of values to plot
#fdir = '../MD_dCSE/src_code/results/'
fdir = '/home/es205/scratch/Re400/iter1918000_to_2233899/'
fielddict = MD_PostProc(fdir)
print(fielddict)

#Select initial values
field = fielddict.plotlist.values()[8]
print("Starting with plot of " + fielddict.plotlist.keys()[8])
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
initpos = int(len(field.grid[1])/2.)
initrec = int(field.maxrec/2.)
ax1, ax2, data = field.contour(axes=naxes,startrec=initrec,endrec=initrec,binlimits=[None,(initpos,initpos),None])
colormesh = ax.pcolormesh(ax1, ax2, data[:,1:,component],cmap=cmap,shading='gourand')

#Setup sliders 
axcolor = 'snow'
axbin = plt.axes([0.25, 0.1 , 0.65, 0.03], axisbg=axcolor)
axrec = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
sbin = Slider(axbin, 'Bin', 0, len(field.grid[1]), valinit=initpos)
srec = Slider(axrec, 'Record', 0, field.maxrec, valinit=initrec)

#Setup  buttons
rax = plt.axes([0.025, 0.3, 0.15, 0.5], axisbg=axcolor)
radio = RadioButtons(rax, fielddict.plotlist.keys(), active=0, activecolor='k')

#Setup main plot axis
ax.axis([ax1.min(), ax1.max(), ax2.min(), ax2.max()])

#Functionality
srec.on_changed(update)
sbin.on_changed(update)
radio.on_clicked(plotfield)
plt.show()
