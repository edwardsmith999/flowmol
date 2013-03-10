#! /usr/bin/env python
from CPLPlot import *

#wdir = '/home/djt06/Documents/Academia/PhD/Code/Data/2013/02/19/NCER/branch/coupler_dCSE/src_code/'
wdir = './../../'

CPLfig = PlotCoupled(wdir)
CPLfig.start_animation()
#CPLfig.anim.save('movie',fps=20)
plt.show()
