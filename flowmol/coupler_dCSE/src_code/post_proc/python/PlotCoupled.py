#! /usr/bin/env python2.7
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from CPLData import *

"""
	class PlotCoupled():
		
		def __init__(self,homefilepath):	
		def start_animation(self,lastrec=False):
			def update_animation(self,dummy,lastrec):
				def update_MD():
				def update_an_sol():
				def update_CFD(): 
		def plot_rec(self,rec):


	Instances of PlotCoupled must be created with specification of the
	absolute file path of the coupled home directory. At the time of writing
	this is located at _____/branch/coupler_dCSE/src_code/.
	
	The file path must be of string format and must include the trailing
	forward-slash character.

"""

# Class to loop through data and plot
class PlotCoupled():
	
	def __init__(self,homefilepath):	

		# Keep track of records
		self.rec = 0
		# Create coupled results object
		self.CPL = CPLData(homefilepath)

		# Initialise figure axes for velocity profiles
		self.fig = plt.figure()
		self.ax  = self.fig.add_subplot(111)
		self.ax.set_xlim([-0.1,1.1])
		self.ax.set_ylim([-0.1,1.1])

		# Initialise MD points object (no data yet)
		self.MDpts = self.ax.plot(
		                          [], 
		                          'ro', 
		                          mec='none', 
		                          label='MD'
		                          )[0] 
		self.MDcon = self.ax.plot(
		                          [],
		                          'mo', 
		                          mec='none',
		                          label='Constrained MD'
		                          )[0] 
		self.MDpts.set_data([],[])
		self.MDcon.set_data([],[])

		# Initialise empty line object for analytical soln
		self.an_sol = self.ax.plot(
		                           [],
		                           'b--', 
		                           label='Analytical'
		                          )[0]
		self.an_sol.set_data([],[])

		# Initialise CFD points objects (no data yet)
		self.CFDpts = self.ax.plot(
		                           [],
		                           'gx',
		                           mew=2.0,
		                           label='CFD'
		                          )[0]
		self.CFDhal = self.ax.plot(
		                           [],
		                           'o',
		                           mew=2.0, 
		                           mec='g',
		                           ms=10.0,
		                           mfc='none',
		                           label='CFD Halo'
		                          )[0]
		self.CFDpts.set_data([],[])
		self.CFDhal.set_data([],[])

		# Set legend and axes labels	
		self.ax.legend(loc='lower right')
		self.ax.set_xlabel('y/L_y')	
		self.ax.set_ylabel('u/U_wall')	

	def start_animation(self,lastrec=False):
		# Create animated figure when called. The optional input lastrec,
		# when True, will watch only the very last records of each sim. This
		# is useful when attempting to "watch" the progress of a coupled run.

		def update_animation(dummy,lastrec):
			# Update data function required by matplotlib's FuncAnimation 

			def update_MD():
				# Get velocity profile for MD region. If the call fails, 
				# try to get the final record by setting last=True (see 
				# get_vprofile routine in MDData.py)
				try:
					bcenters,vprofile = self.CPL.MD.get_vprofile(last=lastrec)
				except:
					bcenters,vprofile = self.CPL.MD.get_vprofile(last=True)
				yspace = bcenters/self.CPL.yL
				self.MDpts.set_xdata(vprofile)
				self.MDpts.set_ydata(yspace)

				# Update constrained cell highlighters
				botcon = self.CPL.MDcnstbins[0]
				topcon = self.CPL.MDcnstbins[-1]+1 # To be inclusive
				self.MDcon.set_xdata(vprofile[botcon:topcon])
				self.MDcon.set_ydata(yspace[botcon:topcon])

			def update_an_sol():
				# Update the analytical solution: find time, then get start-up
				# Couette flow solution based on that time and a previously 
				# defined Reynolds number and geometry (see 
				# CouetteAnalytical.py).	
				t = self.CPL.get_rectime(self.rec)
				yspace, vprofile = self.CPL.analytical.get_vprofile(t)
			
				# Map normalised Couette solution to fluid region (i.e.
				# considering the MD wall)
				self.an_sol.set_xdata(vprofile)
				self.an_sol.set_ydata(yspace)

			def update_CFD():
				# Get velocity profile from CFD data. 
				try:
					vprofile = self.CPL.CFD.get_vprofile(last=lastrec)
				except:
					vprofile = self.CPL.CFD.get_vprofile(last=True)
				yspace = self.CPL.CFDyspace / self.CPL.yL 
				self.CFDpts.set_xdata(vprofile)
				self.CFDpts.set_ydata(yspace)

				halos = np.array( 
				                  [ [vprofile[0], vprofile[-1] ],
				                    [ yspace[0] ,  yspace[-1]  ] ]
				                ) 
				self.CFDhal.set_xdata(halos[0])
				self.CFDhal.set_ydata(halos[1])

			# Perform update
			if (lastrec == False): self.rec += 1
			update_MD()
			update_an_sol()
			update_CFD()

		# Setup the animation based on functions above
		self.anim = ani.FuncAnimation(
		                               self.fig,
		                               update_animation,
		                               fargs=[lastrec],
		                               interval=40
		                             )

	def plot_singlerec(self,rec):
		# Plot a single record chosen by the user

		# MD	
		bincenters,vprofile = self.CPL.MD.get_vprofile(rec=rec)
		yspace = bincenters/self.CPL.yL
		self.ax.plot(vprofile,yspace,'ro', mec='none',label='MD') 

		botcon = self.CPL.MDcnstbins[0]
		topcon = self.CPL.MDcnstbins[-1]+1 # To be inclusive
		self.ax.plot(vprofile[botcon:topcon],yspace[botcon:topcon],
		             'mo',mec='none', label='Constrained MD')

		# Analytical
		t = self.CPL.get_rectime(rec)
		yspace, vprofile = self.CPL.analytical.get_vprofile(t)
		self.ax.plot(vprofile,yspace)

		# CFD
		vprofile = self.CPL.CFD.get_vprofile(rec=rec)
		yspace = self.CPL.CFDyspace / self.CPL.yL
		self.ax.plot(vprofile,yspace) 

