import numpy as np
from MDRawData import MD_RawData

# Abstract Field Class
class Field():

	"""
		Field:   An abstract data averaging base class to be inherited by 
		         derived bins/slices classes.
		Authors: David Trevelyan & Ed Smith, April 2013

		Field contains the information necessary to read a binary data file
		using the MD_RawData class, and average the information spatially
		and temporally.

		The following attributes are required by MD_RawData and should be 
		specified by the derived classes:
			
			fname   - file name of binary data file of interest, 
			dtype   - binary data type string, 'i'/'d' for int/float, and
			nperbin - number of values per bin. 
		
		The important functions are get_bins and get_slices. For ease of use,
		they SHOULD BE APPROPRIATELY ALIASED AS get_field IN THE DERIVED 
		CLASSES (so the user only has to remember one function name, regardless
		of whether the data output is in slices or bins). 
		
		Functionality:

			minrec, maxrec - specify a minimum and maximum record in time
			                 over which to average the data.

			sumaxes        - specify axes over which to sum the data in the
			                 4D array.
			
			meanaxes       - specify axes over which to average the data in 
			                 the 4D array.

		EXAMPLE USE OF GET_BINS:
				
		Read data and average over records 100 to 199:
	
			get_bins(100,200) 
				
		Read data, average over recs 100 to 199, and sum values in x-z plane
		to get the TOTAL in a y-slice:
	
			get_bins(100,200,sumaxes=(0,2)) 

		Read data, average over recs 100 to 199, and average values over the 
		z-direction to get an x/y field:

			get_bins(100,200,meanaxes=(2)) 

	"""

	def __init__(self,fdir,cpol_bins=False): 

		self.fdir = fdir
		self.cpol_bins = cpol_bins	

	def get_bins(self,minrec,maxrec,sumaxes=(),meanaxes=()):
		
		"""
			Read data file using MD_RawData class, average over time record
			range specified by minrec to maxrec. Sum over sumaxes (tuple), or
			average in directions specified by meanaxes (tuple).

			CURRENTLY CANNOT DO BOTH SUM AND MEAN BECAUSE THE SHAPE OF THE 
			ARRAY CHANGES AFTER YOU DO ONE OR THE OTHER.

		"""

		if ( sumaxes != () ) and ( meanaxes != () ):
			print("""Currently cannot specify both sumaxes and meanaxes
			       because the shape of the array of interest changes once
			       it is summed/averaged for the first time. It is possible
			       to determine how meanaxes should be changed based on the
			       specification of sumaxes, but this just hasn\'t been
			       done yet.""")
			quit()

		# Create raw data reading object and get the topology of the bins
		Raw = MD_RawData(self.fdir,self.fname,self.dtype,self.nperbin,self.cpol_bins)
		nbins, binspaces = Raw.get_bintopology()

		# Get bin data from file and get the mean over time records (axis 4)
		bins = Raw.get_bindata(minrec,nrecs=maxrec-minrec)
		bins = bins.mean(axis=4) 

		# Sum or mean as appropriate
		if (sumaxes != ()):
			bins = bins.sum(axis=sumaxes)
			#meanaxes = meanaxes - len(sumaxes)
		if (meanaxes != ()):
			bins = bins.mean(axis=meanaxes)
		
		return bins, binspaces

	def get_slices(self,minrec,maxrec,sumaxes=(),meanaxes=()):
	
		print("""get_slices has not yet been developed in the Field class.
		       Aborting post-processing.""") 	
		quit()

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
	
		print('Getting velocity field from records ' + str(minrec) + ' to ' 
		      + str(maxrec) + ', averaging over axes ' + str(sumaxes) + '.')
	
		msum, binspaces = self.mdata.get_bins(minrec,maxrec,sumaxes=sumaxes)
		psum, binspaces = self.pdata.get_bins(minrec,maxrec,sumaxes=sumaxes)

		# Divide and patch any NaNs
		vfield = np.divide(psum,msum) 
		vfield[np.isnan(vfield)] = 0.0
		
		return vfield, binspaces


