#! /usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class WatchOutput:

	def __init__(self,filename,figureinfo):	

		self.filename = filename
		self.haxis = figureinfo[0]
		self.vaxis = figureinfo[1]
	
		self.fig = plt.figure()
		self.ax  = self.fig.add_subplot(111)

		self.lines = []
		for i in range(len(self.vaxis)):
			self.lines.append(self.ax.plot([])[0]) # First item of ax.plot
			self.lines[i].set_data([],[])
	
		self.setup_anim()

	def setup_anim(self):
		self.anim = animation.FuncAnimation(
		                                    self.fig,
		                                    self.update,
		                                    self.data_gen,
		                                    interval=50
		                                  )
	def update(self,data):
		self.ax.relim()
		self.ax.autoscale_view(True,True,True)
		for i in range(len(self.vaxis)):
			haxis = self.haxis
			vaxis = self.vaxis[i]
			self.lines[i].set_xdata(data[haxis])
			self.lines[i].set_ydata(data[vaxis])

	def data_gen(self):
		while True:
			filedata = np.genfromtxt(self.filename,delimiter=';',names=True)
			yield filedata

class Figure:

	def __init__(self,args):

		fname = args[0]
		haxis = args[1]
		vaxis = args[2:]

		self.filename = fname 
		self.figdef = [haxis,vaxis] 
	
# -------------------------------------------

def SetupFigures(args):

	try:
		figargs = []
		figures = []
		for arg in args:
			if (arg != '+'):
				figargs.append(arg)
			else:
				figures.append(Figure(figargs))
				figargs = []		
		figures.append(Figure(figargs))
		return figures 
	except:
		message = (
		            "It looks like you did something wrong, so here's \n" +
		            "some help. Watch_outputs.py requires arguments that \n" +
		            "specify which figures you would like to watch during \n"+
		            "the simulation. They take the form: \n\n"+
		            "  ./watch_outputs.py filename1 horizontal_axis1 \n" +
		            "                     {vertical_axes1} + filename2 \n "+
		            "                     horizontal_axis2 {vertical_axes2}\n" +
		            "                     + ... \n" +
		            "\n\n " +
		            "Any number of extra figures may be generated simply by \n"+
		            "separating arguments with + signs. The arguments \n" +
		            "horizontal_axis and vertical_axes must be EXACT \n" +
		            "matches to the first-line header of the columns in \n" +
		            "the output file you are trying to follow. Let's say, \n"+
		            "for example, you'd like to follow two figures. One for \n"+
		            "plotting the energies, and the other for plotting the \n" +
		            "pressure. The appropriate call would look like: \n\n" +
		            "  ./watch_outputs.py results/macroscopic_properties \n" +
		            "                     simtime KE PE TE + \n" +
		            "                     results/macroscopic_properties \n" +
		            "                     simtime Pressure \n"
		          )
		print(message)
		quit()

# --------------------------------------------

figures = SetupFigures(sys.argv[1:])
	
for f in figures:
	fig = WatchOutput(f.filename,f.figdef)

plt.show()
