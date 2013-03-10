#! /usr/bin/env python
import numpy as np 
from HeaderData import *

class MDData:
	
	def __init__(self,fdir):
		fobj        = open(fdir+'simulation_header','r') # Object
		self.header = HeaderData(fobj)
		self.mbins  = open(fdir+'mbins','rb')
		self.vbins  = open(fdir+'vbins','rb')
		self.Pbins  = open(fdir+'pVA','rb')
		self.bincenters = self.get_bincenters()

	def read_bins(self,fobj,dtype,nperbin,lastrec=False,seekrec=0,whence=1):
	# Read next record in object (unless call asks for final or a seek pos)

		nbins   = np.array([int(self.header.gnbins1),
		                    int(self.header.gnbins2),
		                    int(self.header.gnbins3)])
		recitems = np.product(nbins)*nperbin 

		if (dtype == 'i'):
			recbytes = 4*recitems
		elif (dtype == 'd'):
			recbytes = 8*recitems
		else:
			quit('Unrecognised data type in read_bins')

		if (lastrec == True): # Seek from end of file
			seekbyte = -1.0*recbytes
			fobj.seek(seekbyte,2)
		else: # seekrec=0 and whence=1 (i.e. continue) by default
			seekbyte = seekrec*recbytes
			fobj.seek(seekbyte,whence)

		# Fortran array ordering
		bindata  = np.fromfile(fobj,dtype=dtype,count=recitems)	
		bindata  = np.reshape(bindata,[nbins[0],nbins[1],nbins[2],nperbin],order='F')
		return bindata

	def get_vprofile(self,last=False,rec=-1):

		# If user wants "final" profile
		if (last == True):
			massbins = self.read_bins(self.mbins,'i',1,lastrec=True)
			velobins = self.read_bins(self.vbins,'d',3,lastrec=True)
		# Else if rec is (not stupidly) specified
		elif (rec > -1):
			massbins = self.read_bins(self.mbins,'i',1,seekrec=rec,whence=0)
			velobins = self.read_bins(self.vbins,'d',3,seekrec=rec,whence=0)
		# Otherwise just work through file continuing from current pos
		else:
			massbins = self.read_bins(self.mbins,'i',1)
			velobins = self.read_bins(self.vbins,'d',3)

		# Old version of np.sum (axis doesn't accept tuple)
		summ  = massbins.sum(axis=0).sum(axis=1)
		sumv  = velobins.sum(axis=0).sum(axis=1)

		vprofile = []
		for ybin in range(int(self.header.gnbins2)):
			vxbin = sumv[ybin][0]/float(summ[ybin])
			vprofile.append(vxbin)
		yspace = self.bincenters

		return yspace, vprofile	
	
	def get_Pprofile(self,last=False,rec=-1):	
		yspace = np.linspace(0.0,1.0,num=100)
		Pprofile = np.power(yspace,2.0)

		# If user wants "final" profile
		if (last == True):
			Pbins = self.read_bins(self.Pbins,'d',9,lastrec=True)
		# Else if rec is (not stupidly) specified
		elif (rec > -1):
			Pbins = self.read_bins(self.Pbins,'d',9,seekrec=rec,whence=0)
		# Otherwise just work through file continuing from current pos
		else:
			Pbins = self.read_bins(self.Pbins,'d',9)

		Pslice = Pbins.mean(axis=0).mean(axis=1)

		Pprofile = []
		for ybin in range(int(self.header.gnbins2)):
			Pxybin = Pslice[ybin][1]
			Pprofile.append(Pxybin)
		yspace = self.bincenters
		
		return yspace,Pprofile
		
	def get_bincenters(self):
		nbins   =   int(self.header.gnbins2)
		binsize = float(self.header.binsize2)
		botbincenter = binsize/2.0 
		topbincenter = nbins*binsize - binsize/2.0
		centers = np.linspace(botbincenter,topbincenter,num=nbins)
		return centers 
