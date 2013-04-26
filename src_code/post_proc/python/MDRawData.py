#! /usr/bin/env python
import numpy as np 
from HeaderData import *

class MD_RawData:
	
	def __init__(self,fdir,fname,dtype,nperbin,cpol_bins=False):

		# 	fdir       -  file directory containing results 
		# 	fname      -  file path from which to read raw data
		#	dtype      -  datatype string
		# 	nperbin    -  number of items to read per bin
		#   cpol_bins  -  boolean flag indicating cylindrical polar bins	

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

		# Inputs:
		# 	seekrec    -  seek a specific record with this integer
		# 	whence     -  specify where to start seeking and which direction

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
