#! /usr/bin/env python
import numpy as np 
from HeaderData import *

"""

	MD_RawData Class
	Author: David Trevelyan, April 2013

	The MD_RawData class is associated with both a data file (e.g.
	mbins, pVA, etc.) and a results directory in which the 
	simulation_header is located. This header contains necessary 
	information that is used to reshape the 1D array of data read 
	from the file into a format that is easy to process later on.
	
	MD_RawData can read any binary data file from the MD code, but
	requires knowledge of the data type (i.e. integer or real) and 
	the number of values per averaging "bin" to read (i.e. velocity
	would require 3 values per bin, one for each cartesian direction).
	
	The optional argument cpol_bins is true if the data has been
	averaged in the MD code using cylindrical polar bins: in this 
	case the only difference lies in creating the mesh of the 
	bin topology, because the domain size in cylindrical polar
	coordinates is different to the cartesian domain size that is
	written to the simulation_header.

	This class is designed for use in two ways:
		
		1) As a contained object within another that packages
	 	   and returns the data in a 3D field format (from which 
		   it may be averaged/summed/sliced into a new format that 
		   you may want to plot), or

		2) For inspection of the numerical values in any binary
		   file.

"""

class MD_RawData:
	
	def __init__(self,fdir,fname,dtype,nperbin,cpol_bins=False):

		"""
			fdir       -  file directory containing results, string
			fname      -  file path from which to read raw data, string
			dtype      -  datatype string, 'i' for integer, 'd' for float
			nperbin    -  number of items to read per bin, integer
			cpol_bins  -  boolean flag indicating cylindrical polar bins
		"""

		self.fdir = fdir
		self.cpol_bins = cpol_bins
		self.header = HeaderData(open(fdir+'simulation_header','r'))

		try: 
			self.fobj = open(fdir+fname,'rb')
			print('Detected file ' + fname )
		except:
			print('Unable to find file ' + fname )
			quit()
		
		self.dtype = dtype
		self.nperbin = nperbin
		self.nbins, self.binspaces = self.get_bintopology()

	def get_bintopology(self):

		"""
			Returns:
			
				gnbins   - A length-3 list of the number of bins in each
						   direction, and
				binspace - A length-3 list of numpy linspaces specifying
						   the locations of the center of each bin in a
						   uniform grid (one linspace for each direction)

		"""
		
		gnbins  = ([ int(self.header.gnbins1), 
		             int(self.header.gnbins2),
		             int(self.header.gnbins3) ])

		if (self.cpol_bins == True):
			domain = ([ float(self.header.r_io) - float(self.header.r_oi), 
			            2.0*np.pi,
			            float(self.header.globaldomain3) ])
		else:
			domain = ([ float(self.header.globaldomain1),
			            float(self.header.globaldomain2),
			            float(self.header.globaldomain3) ])

		binspaces = [] 
		for ixyz in range(3):
			binsize = np.divide(domain[ixyz],gnbins[ixyz])
			botbincenter = binsize/2.0 
			topbincenter = gnbins[ixyz]*binsize - binsize/2.0
			binspaces.append(np.linspace(botbincenter,
			                           topbincenter,
			                           num=gnbins[ixyz]))

		return gnbins, binspaces


	def get_bindata(self, seekrec=0, whence=0):

		"""
			Inputs:

				seekrec - seek a specific record with this integer
				whence  - specify where to start seeking and which 
			              direction

			Return:
				
				bindata - 4D array of data in one record that was
				          read from the binary data file. The size
				          is (nbinsx, nbinsy, nbinsz, nperbin) or
				          the equivalent in cylindrical polar.
				
		"""

		recitems = np.product(self.nbins)*self.nperbin 

		if (self.dtype == 'i'):
			recbytes = 4*recitems
		elif (self.dtype == 'd'):
			recbytes = 8*recitems
		else:
			quit('Unrecognised data type in read_bins')

		# Seek to correct point in the file
		# seekrec=0 and whence=0 (i.e. start) by default
		seekbyte = seekrec*recbytes
		self.fobj.seek(seekbyte,whence)

		# Get data and reshape with fortran array ordering
		bindata = np.fromfile(self.fobj,dtype=self.dtype,count=recitems)	
		bindata = np.reshape( bindata,
		                      [ self.nbins[0],
		                        self.nbins[1],
		                        self.nbins[2],
		                        self.nperbin ],
		                      order='F' )
		return bindata
