#! /usr/bin/env python
from CPLPlot import *

wdir = './../../'

CPLfig = CPL_Plot(wdir,v=True,P=True)
CPLfig.start_animation()
#CPLfig.anim.save('movie',fps=20)
plt.show()
