#! /usr/bin/env python
import numpy as np 

class CFDData:
	
	def __init__(self,fdir):
		self.fdir = fdir
		self.get_reportdata()
		self.subdoms = self.get_subdomlist()

	def get_reportdata(self):
		# Get nx,ny,nz,xL,yL,zL,dx,dy,dz
		fobj = open(self.fdir+'report','r')
		report = fobj.readlines()[3:6] # Lines with info in
		for line in report:
			linepairs = line.split('|')
			for pair in linepairs:
				varname = pair.split()[0]
				varval  = pair.split()[1]
				vars(self)[varname] = varval	
		self.nrx = int(self.nx)-2 # number of subdom records in x
		self.nry = int(self.ny)-1 # no idea why only -1
		self.nrz = int(self.nz)-2
		self.xL = float(self.xL)
		self.yL = float(self.yL)
		self.zL = float(self.zL)
		self.dx = self.xL/float(self.nrx+2-3)
		self.dy = self.yL/float(self.nry+1-3)
		self.dz = self.zL/float(self.nrz+2-3)

	def get_subdomlist(self):
		import os

		def get_int(name):
			string, integer = name.split('.')
			return int(integer)

		subdoms = []
		for filename in os.listdir(self.fdir):
			if (filename.find('SubDom') != -1):
				subdoms.append(filename)

		subdoms = sorted(subdoms,key=get_int)
		return subdoms

	def get_vprofile(self,last=False,rec=-1):

		if (last==True):
			# Get SubDom list again and pop last value
			fpath = self.fdir + self.get_subdomlist().pop()
		elif (rec > -1):
			# If specific record is specified by the user	
			fpath = self.fdir + self.get_subdomlist().pop(rec)
		else:
			# Pop first value from current SubDom list
			fpath = self.fdir + self.subdoms.pop(0)
	
		fobj = open(fpath,'rb')
		subdata = self.read_subfile(fobj,self.nrx,self.nry,self.nrz)
		vx = subdata.mean(axis=0).mean(axis=0) # Avg over x and z dirs

		vprofile = []
		for y in range(self.nry):
			vprofile.append(vx[y][0])

		return vprofile

	def read_subfile(self,fobj,nrx,nry,nrz):
		subdata = np.fromfile(fobj,dtype='d')
		nvar = len(subdata)/(nrx*nry*nrz)
		# Ordered z,x,y for some reason
		subdata = np.reshape(subdata,[nrz,nrx,nry,nvar],order='F')
		return subdata

