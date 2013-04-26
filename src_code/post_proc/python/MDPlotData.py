from MDFields import *

class MD_PlotData():

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
		avgaxes = []	
		for ax in range(3):
			if (ax != axis): avgaxes.append(ax)
		avgaxes = tuple(avgaxes)

		# Instantiate velocity data object
		vData = VBins(self.fdir,cpol_bins=self.cpol_bins)
		# Get 3D velocity field in 3D bins
		vslice, binspaces = vData.get_field(minrec,maxrec,meanaxes=avgaxes)
	
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
		mslice, binspaces = mData.get_bins(minrec,maxrec,meanaxes=(avgaxes))

		return binspaces[axis], mslice
