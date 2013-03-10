#! /usr/bin/env python2.7
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from CPLData import *

"""

	class CPL_Subplot()              - Class containing axes, data, update()
	class Plot_CPL(CPL_profile)      - Create one or more CPL_Subplot objects
	                                   and call their update routines 

	Instances of Plot_CPL must be created with specification of the
	file path of the coupled home directory. At the time of writing
	this is located at _____/branch/coupler_dCSE/src_code/. It is also
	possible to set the velocity and pressure profile plot on and off
	by initialising with v=bool and P=bool, where bool is True or False.

	Example:
		wdir = './../../'
		CPLPlot = Plot_CPL(wdir,v=True,P=False)

	The file path must be of string format and must include the trailing
	forward-slash character.

"""

class CPL_Subplot():
	
	def __init__(self,CPLobj,ptype,figobj,rows,cols,num):

		# Save CPL data object and profile type
		self.CPLobj = CPLobj
		self.ptype = ptype

		# Set variables based on profile type
		if (ptype=='velocity'):
			xlab = 'u/U_wall'
			ylab = 'y/y_L'
			xlim = [-0.1,1.1]
			ylim = [-0.1,1.1]
		elif (ptype=='pressure'):
			xlab = 'P_xy/P_0'
			ylab = 'y/y_L'
			xlim = [-0.1,0.1]
			ylim = [-0.1,1.1]

		# Initialise axes
		self.ax = figobj.add_subplot(rows,cols,num)
		self.ax.set_xlim(xlim)
		self.ax.set_ylim(ylim)
		# Initialise MD points object (no data yet)
		self.MDpts  = self.ax.plot([],'ro',mec='none',label='MD')[0] 
		self.MDcon  = self.ax.plot([],'mo',mec='none',label='Constr. MD')[0] 
		# Initialise empty line object for analytical soln
		self.an_sol = self.ax.plot([],'b--',label='Analytical')[0]
		# Initialise CFD points objects (no data yet)
		self.CFDpts = self.ax.plot([],'gx',mew=2.0,label='CFD')[0]
		self.CFDhal = self.ax.plot([],'o',mew=2.0,mec='g',ms=10.0,mfc='none',
								   label='CFD Halo')[0]
		# Set legend and axes labels	
		self.ax.set_xlabel(xlab)	
		self.ax.set_ylabel(ylab)	
		plt.figlegend(
		               ( self.MDpts, self.MDcon, 
		                 self.CFDpts, self.CFDhal,self.an_sol ),
		               ( "MD","C->P","CFD","CFD Halo","Analytical"),
		               'upper center',ncol=5
		             )
		# Initialise empty data sets
		self.MDpts.set_data([],[])
		self.MDcon.set_data([],[])
		self.an_sol.set_data([],[])
		self.CFDpts.set_data([],[])
		self.CFDhal.set_data([],[])

	def get_MD_profile(self,rec,last=False):
		# Return the velocity or pressure profile from MD Data object
		if (self.ptype=='velocity'):
			return self.CPLobj.MD.get_vprofile(rec,last=last)	
		elif (self.ptype=='pressure'):
			return self.CPLobj.MD.get_Pprofile(rec,last=last)	

	def get_analytical(self,rec):
		# Return the analytical velocity or pressure profile 
		t = self.CPLobj.get_rectime(rec)
		if (self.ptype=='velocity'):
			return self.CPLobj.analytical.get_vprofile(t)
		elif (self.ptype=='pressure'):
			return self.CPLobj.analytical.get_Pprofile(t)
	
	def get_CFD_profile(self,rec,last=False):
		# Return the CFD velocity or pressure profile 
		if (self.ptype=='velocity'):
			return self.CPLobj.CFD.get_vprofile(rec,last=last)	
		elif (self.ptype=='pressure'):
			return self.CPLobj.CFD.get_Pprofile(rec,last=last)	
	
	def update(self,rec,lastrec=False):
		# Update plot data to the next, last or specified record

		def update_MD():
			try:
				bcenters,profile = self.get_MD_profile(rec,last=lastrec)
			except:
				bcenters,profile = self.get_MD_profile(rec,last=True)
			yspace = bcenters/self.CPLobj.yL
			self.MDpts.set_xdata(profile)
			self.MDpts.set_ydata(yspace)
			# Update constrained cell highlighters
			bot = self.CPLobj.MDcnstbins[0]
			top = self.CPLobj.MDcnstbins[-1]+1 # To be inclusive
			self.MDcon.set_xdata(profile[bot:top])
			self.MDcon.set_ydata(yspace[bot:top])

		def update_analytical():
			yspace, profile = self.get_analytical(rec)
			self.an_sol.set_xdata(profile)
			self.an_sol.set_ydata(yspace)

		def update_CFD():
			try:
				vprofile = self.get_CFD_profile(rec,last=lastrec)
			except:
				vprofile = self.get_CFD_profile(rec,last=True)
			yspace = self.CPLobj.CFDyspace / self.CPLobj.yL 
			self.CFDpts.set_xdata(vprofile)
			self.CFDpts.set_ydata(yspace)
			halos = np.array( 
							  [ [vprofile[0], vprofile[-1] ],
								[ yspace[0] ,  yspace[-1]  ] ]
							) 
			self.CFDhal.set_xdata(halos[0])
			self.CFDhal.set_ydata(halos[1])

		# Perform update
		update_MD()
		update_analytical()
		update_CFD()

# Class to loop through data and plot
class CPL_Plot(CPL_Subplot):

	def __init__(self,homefilepath,v=True,P=True):	

		# Store on/off switches for v and P plots
		self.v = v
		self.P = P

		# Create coupled results object
		self.CPL = CPLData(homefilepath)

		# Calculate number of rows and columns
		rows = 1
		cols = int(v) + int(P)

		# Initialise figure axes for velocity profiles
		self.fig = plt.figure(figsize=[cols*7.,6.])

		# Initialise CPL_Subplot objects
		num = 0
		if (v==True):
			num += 1
			self.vplot = CPL_Subplot(self.CPL,'velocity',self.fig,
			                         rows,cols,num)
		if (P==True):
			num += 1
			self.Pplot = CPL_Subplot(self.CPL,'pressure',self.fig,
			                         rows,cols,num)
		# Keep track of records
		self.rec = 0

	def update_animation(self,dummy,lastrec):
		# Update data function required by matplotlib's FuncAnimation 
		if (lastrec == False): self.rec += 1
		if (self.v==True): self.vplot.update(self.rec,lastrec)
		if (self.P==True): self.Pplot.update(self.rec,lastrec)

	def start_animation(self,lastrec=False):
		# Create animated figure when called. The optional input lastrec,
		# when True, will watch only the very last records of each sim. This
		# is useful when attempting to "watch" the progress of a coupled run.
		self.anim = ani.FuncAnimation(
		                               self.fig,
		                               self.update_animation,
		                               fargs=[lastrec],
		                               interval=100
		                             )
	
	def plot_singlerec(self,rec):

		if (self.v==True): 
			# MD	
			bincenters,vprofile = self.vplot.get_MD_profile(rec)
			yspace = bincenters/self.CPL.yL
			self.vplot.ax.plot(vprofile,yspace,'ro', mec='none') 
			bot = self.CPL.MDcnstbins[0]
			top = self.CPL.MDcnstbins[-1]+1 # To be inclusive
			self.vplot.ax.plot(vprofile[bot:top],yspace[bot:top],'mo',mec='none')

			# Analytical
			yspace, vprofile = self.vplot.get_analytical(rec)
			self.vplot.ax.plot(vprofile,yspace,'b--')

			# CFD
			vprofile = self.vplot.get_CFD_profile(rec)
			yspace   = self.CPL.CFDyspace / self.CPL.yL
			self.vplot.ax.plot(vprofile,yspace,'gx',mew=2.0) 
			halos = np.array( 
							  [ [vprofile[0], vprofile[-1] ],
								[ yspace[0] ,  yspace[-1]  ] ]
							)
			self.vplot.ax.plot(halos[0],halos[1],'o',mew=2.0,mec='g',ms=10.0,mfc='none')

		if (self.P==True):
			pass
		
