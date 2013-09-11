#! /usr/bin/env python
import sys
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

processname = sys.argv[1]
outfile = sys.argv[2]

class WatchFigure:

	def __init__(self):

		# Total memory history, iteration history
		self.tmem_history = []
		self.iter_history = []
		
		# Own figure and axes
		self.fobj = plt.figure()
		self.ax = self.fobj.add_subplot(111)
		# Empty line data
		self.line2D = matplotlib.lines.Line2D([],[])
		# Add line to axis
		self.ax.add_line(self.line2D)

		# Setup animation
		self.anim = animation.FuncAnimation(
											self.fobj,
											self.update_axis,
											interval=500
		                                   )


	def update_axis(self,dummy):

		try:

			# Get total memory usage tmem
			cmd = ['tail','-1',outfile]
			iter = int(sp.check_output(cmd).split(';')[0])
			cmd = ('ps -eo comm,vsz').split()
			processes = [p for p in sp.check_output(cmd).split('\n') if processname in p]
			tmem = 0
			for p in processes:
				mem = int(p.split()[1])	
				tmem += mem
	
			# Update data	
			self.tmem_history.append(tmem)
			self.iter_history.append(iter)

		except:
			
			print('Warning: failed to read data.')
			

		# Plot, relimit, autoscale	
		self.line2D.set_xdata(xrange(len(self.iter_history)))
		self.line2D.set_ydata(self.tmem_history)
		self.ax.relim()
		self.ax.autoscale_view(True,True,True)

		#print('Process search string: \t\t' + processname)
		#print('Current iteration:     \t\t' + str(iter))
		#print('Total memory usage:    \t\t' + str(tmem))

Object = WatchFigure()
plt.show()
