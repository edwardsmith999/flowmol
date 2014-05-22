import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

from MDFields import *
from HeaderData import *

component = 0
cmap = plt.cm.RdYlBu_r
naxes = (0,2)

#fdir = '../MD_dCSE/src_code/results/'
fdir = '/home/es205/scratch/Re400/iter1918000_to_2233899/'
field = MD_vField(fdir)

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
initpos = int(len(field.grid[1])/2.)
initrec = int(field.maxrec/2.)
ax1, ax2, data = field.contour(axes=naxes,startrec=initrec,endrec=initrec,binlimits=[None,(initpos,initpos),None])
#contour = plt.pcolormesh(ax1,ax2,data[:,:,component],cmap=cmap)
#plt.colorbar(contour)
#plt.show()
#quit()
colormesh = plt.pcolormesh(ax1, ax2, data[:,:,component],cmap=cmap)
#print(dir(colormesh))
#quit()

plt.axis([ax1.min(), ax1.max(), ax2.min(), ax2.max()])

axcolor = 'lightgoldenrodyellow'
axbin = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
axrec  = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)

sbin = Slider(axbin, 'Bin', 0, len(field.grid[1]), valinit=initpos)
srec = Slider(axrec, 'Record', 0, field.maxrec, valinit=initrec)

def update(val):
    pos = int(np.rint(sbin.val))
    rec = int(np.rint(srec.val))
    ax1, ax2, data = field.contour(axes=naxes,startrec=rec,endrec=rec,binlimits=[None,(pos,pos),None])
    colormesh.set_array(data[:,:,component].ravel())
    #colormesh.set_array(data[:,:,component])
    fig.canvas.draw_idle()

srec.on_changed(update)
sbin.on_changed(update)

#resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
#button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
#def reset(event):
#    srec.reset()
#    sbin.reset()
#button.on_clicked(reset)

plt.show()
