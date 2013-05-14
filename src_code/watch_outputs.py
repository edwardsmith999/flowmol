#! /usr/bin/env python
import sys
import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class DataLine:

	colors = itertools.cycle(['b','g','r','c','m','y','k'])

	def __init__(self, linedef):

		self.filename = linedef[0] 
		self.haxis = linedef[1] 
		self.vaxis = linedef[2]

		label = self.filename + " " + self.vaxis
		color = next(DataLine.colors)

		self.line2D = matplotlib.lines.Line2D([],[],color=color,label=label)
	
	def update_linedata(self):

		filedata = np.genfromtxt(self.filename,delimiter=';',names=True)
		self.line2D.set_xdata(filedata[self.haxis])
		self.line2D.set_ydata(filedata[self.vaxis])
	
class OutputFigure:

	def __init__(self,linedefs):

		# Own figure and axes
		self.fobj = plt.figure()
		self.ax = self.fobj.add_subplot(111)

		# Create DataLine instances 
		self.datalines = []
		for linedef in linedefs:
			self.datalines.append(DataLine(linedef))

		# Add line2D objects to axis
		for dataline in self.datalines:
			self.ax.add_line(dataline.line2D)
			self.ax.set_xlabel(dataline.haxis)

		# Setup animation
		self.anim = animation.FuncAnimation(
		                                    self.fobj,
		                                    self.update_axis,
		                                    interval=50
	   	                                  )
		self.ax.legend()

	def update_axis(self,dummy):

		for dataline in self.datalines:
			dataline.update_linedata()

		self.ax.relim()
		self.ax.autoscale_view(True,True,True)
	

class WatchFigures:

	def __init__(self,figuredefs):
	
		figures = []	
		for figuredef in figuredefs:
			figures.append(OutputFigure(figuredef))
			

def ParseArguments():

	arguments = sys.argv[1:]

	if (sys.argv[2] == 'help'):

		message = (
					". Watch_outputs.py requires arguments that \n" +
					"specify which figures you would like to watch during \n"+
					"the simulation. They take the form: \n\n"+
					"  ./watch_outputs.py filename1 horizontal_axis1 \n" +
					"                     {vertical_axes1} + filename2 \n "+
					"                     horizontal_axis2 {vertical_axes2}\n" +
					"                     ++ new figure arguments... \n" +
					"\n\n " +
					"Any number of extra figures may be generated simply by \n"+
					"separating arguments with TWO + signs (++). The arguments \n" +
					"horizontal_axis and vertical_axes must be EXACT \n" +
					"matches to the first-line header of the columns in \n" +
					"the output file you are trying to follow. Let's say, \n"+
					"for example, you'd like to follow two figures. One for \n"+
					"plotting the energies, and the other for plotting the \n" +
					"pressure. The appropriate call would look like: \n\n" +
					"  ./watch_outputs.py results/macroscopic_properties \n" +
					"                     simtime KE + \n" +
					"                     results/macroscipic_properties \n" +
					"                     simtime PE + \n" +
					"                     results/macroscopic_properties \n" +
					"                     simtime TE     ++ \n" +
					"                     results/macroscopic_properties \n" +
					"                     simtime Pressure \n"
				  )
		print(message)
		quit()

	linedef = []
	linedefs = []
	figuredefs = []
	for arg in arguments:

		linedef.append(arg)

		if ( '+' in arg ):
			linedef.pop()

		if ( arg == '+' ):

			linedefs.append(linedef)	
			linedef = []

		elif ( arg == '++' ):

			linedefs.append(linedef)
			figuredefs.append(linedefs)
			linedef = []
			linedefs = []

	linedefs.append(linedef)
	figuredefs.append(linedefs)

	return figuredefs


##############################################################################
# Parse the command line arguments and try plotting!

figure_definitions = ParseArguments()
Watch = WatchFigures(figure_definitions)
plt.show()
