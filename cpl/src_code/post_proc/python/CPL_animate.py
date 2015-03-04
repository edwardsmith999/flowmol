#! /usr/bin/env python
from CPLPlot import *
import sys 

""" 
	Last command line argument can be "save", in which case
	the animation is saved and not shown.
"""

wdir = './../../'
CPLfig = CPL_Plot(wdir,v=True,P=True)

if (sys.argv[-1] == 'save'):
	CPLfig.start_animation(frames=1000)
	CPLfig.anim.save('./movie.mp4',fps=30,bitrate=1800,writer='mencoder')
else:
	CPLfig.start_animation()
	plt.show()
