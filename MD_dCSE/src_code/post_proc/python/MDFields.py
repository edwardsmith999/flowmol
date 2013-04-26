import numpy as np
from MDData import MD_RawData

# Abstract Field Class
class Field():

	def __init__(self,fdir,cpol_bins=False): 

		self.fdir = fdir
		self.cpol_bins = cpol_bins	

	def get_bins(self,minrec,maxrec,sumaxes=(),meanaxes=()):

		Raw = MD_RawData(self.fdir,self.fname,self.dtype,self.nperbin,self.cpol_bins)
		nbins, binspaces = Raw.get_bintopology()

		bins = Raw.get_bindata(seekrec=minrec)
		for rec in range(minrec+1,maxrec):
			bins += Raw.get_bindata(seekrec=rec)
		bins = bins/float(maxrec - minrec)

		if (sumaxes != ()):
			bins = bins.sum(axis=sumaxes)

		if (meanaxes != ()):
			bins = bins.mean(axis=meanaxes)
		
		return bins, binspaces


# Mass field	
class MassBins(Field):

	fname = 'mbins'
	dtype = 'i'
	nperbin = 1
	get_field = Field.get_bins

# Momentum field	
class MomBins(Field):

	fname = 'vbins'
	dtype = 'd'
	nperbin = 3
	get_field = Field.get_bins

# Pressure field
class PBins(Field):

	fname = 'pVA'
	dtype = 'd'
	nperbin = 9
	get_field = Field.get_bins

# Mass slice?
#class MassSlice(Field):
#	
#	fname = 'mslice'
#	dtype = 'i'
#	nperbin = 1
#	get_field = Field.get_slice

# Velocity field
class VBins():

	def __init__(self,fdir,cpol_bins):
		self.mdata = MassBins(fdir,cpol_bins)
		self.pdata = MomBins(fdir,cpol_bins)

	def get_field(self,minrec,maxrec,sumaxes=()):	

		"""
		    Get the velocity field from file vbins averaged over
		    time records minrec->maxrec, AND averaged over 
		    spatial directions specified by sumaxes.
			
			sumaxes   - *tuple* of *int* between 0-2.
		    minrec    - *int*
			maxrec    - *int*

		"""
	
		msum, binspaces = self.mdata.get_bins(minrec,maxrec,sumaxes=sumaxes)
		psum, binspaces = self.pdata.get_bins(minrec,maxrec,sumaxes=sumaxes)

		# Divide and patch any NaNs
		vfield = np.divide(psum,msum) 
		vfield[np.isnan(vfield)] = 0.0
		
		return vfield, binspaces


