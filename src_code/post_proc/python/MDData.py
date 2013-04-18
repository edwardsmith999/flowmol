#! /usr/bin/env python
import numpy as np 
from HeaderData import *

class MD_RawData:
	
	def __init__(self,fdir):

		self.fdir = fdir
		headerobj   = open(fdir+'simulation_header','r')
		self.header = HeaderData(headerobj)


	def get_bintopology(self):
		
		gnbins  = ([ int(self.header.gnbins1), 
		             int(self.header.gnbins2),
		             int(self.header.gnbins3) ])

		if (self.cpol_bins == True):
			domain = ([ float(self.header.r_oo), 2.0*np.pi,
			            float(self.header.globaldomain3) ])
		else:
			domain = ([ float(self.header.globaldomain1),
			            float(self.header.globaldomain2),
			            float(self.header.globaldomain3) ])

		binmesh = [] 
		for ixyz in range(3):
			binsize = np.divide(domain[ixyz],gnbins[ixyz])
			botbincenter = binsize/2.0 
			topbincenter = gnbins[ixyz]*binsize - binsize/2.0
			binmesh.append(np.linspace(botbincenter,
			                           topbincenter,
			                           num=gnbins[ixyz]))

		return gnbins, binmesh


	def get_bindata(self, fobj, dtype, nbins, nperbin, lastrec=False,
	                seekrec=0, whence=1):
		# Inputs:
		# 	fobj       -  file object to be read
		#	dtype      -  datatype string
		# 	nbins      -  list/array of number of bins (length 3)
		# 	nperbin    -  number of items to read per bin
		# 	lastrec    -  bool set to True if only reading final record
		# 	seekrec    -  seek a specific record with this integer
		# 	whence     -  specify where to start seeking and which direction

		recitems = np.product(nbins)*nperbin 

		if (dtype == 'i'):
			recbytes = 4*recitems
		elif (dtype == 'd'):
			recbytes = 8*recitems
		else:
			quit('Unrecognised data type in read_bins')

		# Seek to correct point in the file
		# seekrec=0 and whence=1 (i.e. continue) by default
		if (lastrec == True): 
			seekbyte = -1.0*recbytes
			fobj.seek(seekbyte,2)
		else:
			seekbyte = seekrec*recbytes
			fobj.seek(seekbyte,whence)

		# Get data and reshape with fortran array ordering
		bindata = np.fromfile(fobj,dtype=dtype,count=recitems)	
		bindata = np.reshape(bindata,[nbins[0],nbins[1],nbins[2],nperbin],
		                     order='F')
		return bindata


	def get_field(self, fobj, dtype, nperbin, rec, last=False):

		# Get number of bins and their locations
		nbins, binmesh = self.get_bintopology()

		# If user wants "final" profile
		if (last == True):
			bindata = self.get_bindata( fobj, dtype, nbins, nperbin,
			                            lastrec=True )
		# Else if rec is specified
		else:
			bindata = self.get_bindata( fobj, dtype, nbins, nperbin,
			                            seekrec=rec, whence=0 )

		return bindata, binmesh

class MD_PlotData:

	def __init__(self,fdir,cpol_bins=False):
		self.Raw = MD_RawData(fdir)
		self.Raw.cpol_bins = cpol_bins

	def get_vplane(self,plane,haxis,vaxis,rec,last=False):

		# Inputs:
		# 	plane     - plane over which to avg field data (eg. 0=x,1=y,2=z)
		# 	haxis     - horizontal axis (eg. 0=x,1=y,2=z)
		# 	vaxis     - vertical axis
		# 	rec       - record to retrieve from raw data
		# 	last      - flag to override rec and get final record if True
	
		# File objects to be passed to self.get_field
		mobj  = open(self.Raw.fdir+'mbins','rb')
		vobj  = open(self.Raw.fdir+'vbins','rb')

		# Get 3D momentum and mass fields
		mbins, binmesh = self.Raw.get_field(mobj,'i',1,rec,last=last)
		vbins, binmesh = self.Raw.get_field(vobj,'d',3,rec,last=last)

		# Sum over axis specified in input for averaging 
		msum = mbins.sum(axis=plane)
		vsum = vbins.sum(axis=plane)

		# Calculate velocity field and patch any NaNs
		vplane = np.divide(vsum,msum)
		vplane[np.isnan(vplane)] = 0.0

		# Get bin center positions on both axes for every field point
		ax1, ax2 = np.meshgrid(binmesh[haxis],binmesh[vaxis],indexing='ij')

		# Extract desired velocity components
		v1 = vplane[:,:,haxis]
		v2 = vplane[:,:,vaxis]

		return ax1, ax2, v1, v2

	def get_avgd_vplane(self,plane,haxis,vaxis,minrec,maxrec):

		if ( minrec >= maxrec ):
			print('Min/Max records incorrect in get_avgd_vplane')
			quit()

		ax1, ax2, v1, v2  = self.get_vplane(plane,haxis,vaxis,0)
		ax1 *= 0.0
		ax2 *= 0.0
		v1  *= 0.0
		v2  *= 0.0

		for i in range(minrec,maxrec):
			ax1, ax2, _v1, _v2 = self.get_vplane(plane,haxis,vaxis,i)
			v1 += _v1
			v2 += _v2
		v1 = v1 / (maxrec - minrec) 
		v2 = v2 / (maxrec - minrec) 

		return ax1, ax2, v1, v2

	def get_vprofile(self,axis,rec,last=False):

		sumaxes = []	
		for ax in range(3):
			if (ax != axis): sumaxes.append(ax)
		sumaxes = tuple(sumaxes)

		# File objects to be passed to self.get_field
		mobj  = open(self.Raw.fdir+'mbins','rb')
		vobj  = open(self.Raw.fdir+'vbins','rb')

		# Get 3D momentum and mass fields
		mbins, binmesh = self.Raw.get_field(mobj,'i',1,rec,last=last)
		vbins, binmesh = self.Raw.get_field(vobj,'d',3,rec,last=last)

		# Sum over axis specified in input for averaging 
		try:
			msum = mbins.sum(axis=sumaxes)
			vsum = vbins.sum(axis=sumaxes)
		except:
			print('The sum of mbins and vbins over certain axes  ' +
				  'failed in MD_PlotData.get_vprofile(). This is ' +
				  'most often caused when your version of numpy (' +
				  str(np.__version__) + ') is earlier than the '   +
				  'required version (1.7.1)')
			quit()

		# Calculate velocity field and patch any NaNs
		vprof = np.divide(vsum,msum)
		vprof[np.isnan(vprof)] = 0.0
	
		return binmesh[axis], vprof

	def get_avgd_vprofile(self,axis,minrec,maxrec):

		if ( minrec >= maxrec ):
			print('Min/Max records incorrect in get_avgd_vprofile')
			quit()

		ax, vprof = self.get_vprofile(axis=axis,rec=0)
		ax *= 0.0
		vprof *= 0.0
		for i in range(minrec,maxrec):
			ax, v = self.get_vprofile(axis=axis,rec=i)
			vprof += v
		vprof = vprof / (maxrec - minrec) 
	
		return ax, vprof
