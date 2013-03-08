#! /usr/bin/env python
from PlotCoupled import *

#wdir = '/home/djt06/Documents/Academia/PhD/Code/Data/2013/02/19/NCER/branch/coupler_dCSE/src_code/'
#wdir = '/home/djt06/Desktop/LJ2/branch/coupler_dCSE/src_code/'
wdir = './../../'

CPLfig = PlotCoupled(wdir)
for rec in [7,15,31,63,127]:
	CPLfig.plot_singlerec(rec)

plt.show()
