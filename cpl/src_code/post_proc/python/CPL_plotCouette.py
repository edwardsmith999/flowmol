#! /usr/bin/env python
from CPLPlot import *

wdir = './../../'

CPLfig = CPL_Plot(wdir,v=True,P=False)
for rec in [7,15,31,63,127]:
	CPLfig.plot_singlerec(rec)

plt.show()
