from MDFields import *

class MD_PlotData():

	"""
		MD_PlotData: Return the arguments required by matplotlib's plotting
		             functions for a specific figure type.  
		Author: David Trevelyan, April 2013

		Each class method should return arguments required by a single 
		matplotlib plotting function for a single figure. Each time you need 
		a new figure, please save a minimal method in this class so that 
		similar figures can be (hopefully!) trivially plotted in future.

		IMPORTANT: It should be possible to get all the information from the
		classes in MDFields.py. If the functionality you need is not there, 
		please add it to those classes rather than changing the MD_RawData 
		class or cluttering this MD_PlotData class. 

		This class will likely get quite long, so it may be sensible to 
		break it up into categories. It's not necessary just yet, though,
		as of April 2013.

	"""

	def __init__(self,fdir,cpol_bins=False):
		self.fdir = fdir
		self.cpol_bins = cpol_bins

	def get_vplane_quiver_args(self,plane,haxis,vaxis,minrec,maxrec):

		# Instantiate velocity data object
		vData = VBins(self.fdir,cpol_bins=self.cpol_bins)
		# Extract 3D velocity field averaged over 1D of bins
		vplane, binspaces = vData.get_field(minrec,maxrec,sumaxes=(plane))
		# Get bin center positions on both axes for every field point
		X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')
		# Extract components of velocity for each bin
		U = vplane[:,:,haxis]
		V = vplane[:,:,vaxis]
		return X, Y, U, V

	def get_vplane_streamplot_args(self,*args,**kw):
	
		# Same arguments as quiver, just transposed	
		X, Y, U, V = self.get_vplane_quiver_args(*args,**kw)
		X = X.transpose()
		Y = Y.transpose()
		U = U.transpose()
		V = V.transpose()
		return X, Y, U, V

	def get_vslice_plot_args(self,axis,minrec,maxrec):

		# Get which axes to average over
		sumaxes = []	
		for ax in range(3):
			if (ax != axis): sumaxes.append(ax)
		sumaxes = tuple(sumaxes)

		# Instantiate velocity data object
		vData = VBins(self.fdir,cpol_bins=self.cpol_bins)
		# Get 3D velocity field in 3D bins
		vslice, binspaces = vData.get_field(minrec,maxrec,sumaxes=sumaxes)
	
		return binspaces[axis], vslice

	def get_vplane_splot_args(self,plane,haxis,vaxis,component,minrec,maxrec):

		# Instantiate velocity data object
		vData = VBins(self.fdir,cpol_bins=self.cpol_bins)
		# Extract 3D velocity field averaged over 1D of bins
		vplane, binspaces = vData.get_field(minrec,maxrec,sumaxes=(plane))
		# Get bin center positions on both axes for every field point
		X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')
		# Extract components of velocity for each bin
		U = vplane[:,:,component]

		return X, Y, U		

	def get_mslice_plot_args(self,axis,minrec,maxrec):

		# Get which axes to average over
		avgaxes = []	
		for ax in range(3):
			if (ax != axis): avgaxes.append(ax)
		avgaxes = tuple(avgaxes)

		mData = MassBins(self.fdir,cpol_bins=self.cpol_bins)
		mslice, binspaces = mData.get_bins(minrec,maxrec,sumaxes=(avgaxes))
		#Take zeroth component to 'Squeeze' array to lower dimensionality
		mslice = mslice[:,0] 

		return binspaces[axis], mslice

	def get_mplane_splot_args(self,plane,haxis,vaxis,minrec,maxrec):

		# Instantiate mass data object
		mData = MassBins(self.fdir,cpol_bins=self.cpol_bins)
		# Extract 3D velocity field averaged over 1D of bins
		mplane, binspaces = mData.get_field(minrec,maxrec,sumaxes=(plane))
		# Get bin center positions on both axes for every field point
		X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')
		#Take zeroth component to 'Squeeze' array to lower dimensionality
		m = mplane[:,:,0] 

		return X, Y, m		

	def get_pVA_prof_args(self,axis,component,minrec,maxrec):
		
		# Get which axes to average over
		avgaxes = []	
		for ax in range(3):
			if (ax != axis): avgaxes.append(ax)
		avgaxes = tuple(avgaxes)
		
		pVA_obj = PBins(self.fdir,cpol_bins=self.cpol_bins)
		pslice, binspaces = pVA_obj.get_field(minrec,maxrec,meanaxes=(avgaxes))
	
		return binspaces[axis], pslice[:,component]

	def get_Tplane_splot_args(self,plane,haxis,vaxis,minrec,maxrec):

		# Instantiate temperature data object
		TData = TBins(self.fdir,cpol_bins=self.cpol_bins)
		# Extract 3D velocity field averaged over 1D of bins
		Tplane, binspaces = TData.get_field(minrec,maxrec,sumaxes=(plane))
		# Get bin center positions on both axes for every field point
		X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')
		#Take zeroth component to 'Squeeze' array to lower dimensionality
		T = Tplane[:,:,0] 

		return X, Y, T		

	def get_vfield_energy_spectra(self,plane,component,minrec,maxrec,tavg_rec,
	                              fftaxis=None,ffttime=False):

		class vfield_spectra:

			def __init__(self,DataObj):
				# Initialise velocity time series to nothing
				self.v_ts = []
				self.vDataObj = DataObj
				# Keep track of which axes have been averaged, summed, etc.
				# This is kept so the user can just specify which axis they
				# want to fft, without having to worry about how the dimensions
				# of the array have changed on averaging/summing etc.
				# THE AXES CORRESPOND TO: [X,Y,Z,TIME,COMPONENT]
				self.axesreduced = [False]*5
				# Local copy of which axes represent real space (each entry is
				# "popped" out when reduced (see self.axesreduced)
				self.realspacepopped = [True]*5
			
			def populate_timeseries(self,minrec,maxrec,drec,sumplane=None):

				# Read time series of velocity field records, averaging over
				# length "drec" before storing in array and summing mass/mom
				# bins in direction "sumplane" so velocity field is averaged
				# in that direction.
				for rec in range(minrec,maxrec,drec):
					temp_vfield, binspaces = self.vDataObj.get_field(rec,
					                         rec+drec,sumaxes=(sumplane))
					self.v_ts.append(temp_vfield)
				self.v_ts = np.array(self.v_ts)

				# Mark "reduction" of sumplane axis
				self.axesreduced[sumplane] = True
				# Remove sumplane from realspace list
				self.realspacepopped.pop(sumplane)

				# Put time axis in third position, component in 4th, remaining
				# spatial coordinates in 1st and 2nd (i.e. 0th and 1st),
				# so now we have u(axis1,axis2,time,component), reduced over 
				# the sumplane axis
				self.v_ts = np.transpose(self.v_ts,axes=(1,2,0,3))

			def extract_component(self,component):
				self.v_ts = self.v_ts[:,:,:,component]
				# 4 for pre sumaxes in populate_timeseries 	
				self.axesreduced[4] = True 
				# 3 for post sumaxes in populate_timeseries
				self.realspacepopped.pop(3)

			def fft(self,fftaxis,window=False):

				def apply_window(a):
					a = a * self.window
					return a

				if (self.axesreduced[fftaxis]):
					print('You\'re trying to FFT over an axis that is ' +
					      'already reduced')
					quit()

				newaxis = fftaxis - np.sum(self.axesreduced[:fftaxis])

				if (window == True):
					N = self.v_ts.shape[newaxis]
					self.window = np.hanning(N)
					self.wss = np.sum(self.window**2.0)/N
					np.apply_along_axis(apply_window,newaxis,self.v_ts)		

				# Perform FFT
				self.v_ts = np.fft.fft(self.v_ts,axis=newaxis)
				# Mark axis as no longer representing real space
				self.realspacepopped[newaxis] = False

			def set_energyfield(self,window=False):

				# Work out which axes we can average over (i.e. any that we
				# haven't Fourier transformed)
				axes = np.where(np.array(self.realspacepopped)==True)[0]

				# Energy field
				self.E = np.abs(self.v_ts)**2.0
		
				if (window == True):
					self.E = self.E / self.wss 

				# Average remaining dimensions
				self.E = np.mean(self.E,axis=tuple(axes))

				# Delete duplicate parts of the energy spectrum and double
				# energy contributions of middle wavenumbers
				n = len(self.E)/2 + 1
				self.E = self.E[:n]
				self.E[1:-1] = self.E[1:-1] * 2.0

		# Create field reading object
		self.vDataObj = VBins(self.fdir,cpol_bins=self.cpol_bins)	
		# Create spectra calculating object
		VField = vfield_spectra(self.vDataObj)

		# Read the time series of the velocity field desired
		VField.populate_timeseries(minrec,maxrec,tavg_rec,sumplane=plane)
		# Extract the cartesian component we are interested in
		VField.extract_component(component)

		# FFT in space or time depending on input to function, and 
		# calculate power spectra.
		if (fftaxis != None):
			VField.fft(fftaxis)
			VField.set_energyfield()

		if (ffttime != False):
			VField.fft(3,window=True)
			VField.set_energyfield(window=True)
	
		# Return energy normalised by N	
		return VField.E/len(VField.E)
