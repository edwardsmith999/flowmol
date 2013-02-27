#! /usr/bin/env python
from PlotCoupled import *

#wdir = '/home/djt06/Documents/Academia/PhD/Code/Data/2013/02/19/NCER/branch/coupler_dCSE/src_code/'

CPLfig = PlotCoupled(wdir)
for rec in [1,3,7,15,31]:
	CPLfig.plot_singlerec(rec)

plt.show()
