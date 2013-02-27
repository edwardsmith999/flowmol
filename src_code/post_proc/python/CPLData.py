#! /usr/bin/env python
import numpy as np 
from HeaderData import *
from MDData import *
from CFDData import *
from CouetteAnalytical import *

class CPLData:

	viscosity = 1.6	
	U_wall    = 1.0

	def __init__(self,fdir):
		fobj           = open(fdir+'/results/coupler_header','r')
		self.header    = HeaderData(fobj)
		self.MD        = MDData(fdir+'/md_data/results/')
		self.CFD       = CFDData(fdir+'/couette_data/')
		self.yL        = self.get_coupleddomain()
		self.CFDyspace = self.get_CFDyspace()
		Re, U, L       = self.get_Couetteparams()
		self.analytical = CouetteAnalytical(Re,U,L)

	def get_coupleddomain(self):
		yL_md = float(self.MD.header.globaldomain2)
		yL_cfd = float(self.CFD.yL)
		dy_cfd = float(self.CFD.dy)
		ncy_olap = int(self.header.jcmax_olap)-int(self.header.jcmin_olap) + 1
		yL_olap = dy_cfd*ncy_olap
		yL = yL_md + yL_cfd - yL_olap
		return yL

	def get_CFDyspace(self):
		CFDtop = self.yL + 0.5*float(self.CFD.dy)
		CFDbot = self.yL - float(self.CFD.yL) - 0.5*float(self.CFD.dy)
		yspace = np.linspace(CFDbot,CFDtop,num=self.CFD.nry)/self.yL
		return yspace	

	def get_Couetteparams(self):
		U = self.U_wall
		L = self.yL 
		Re = float(self.header.density_cfd) * U * 1.0 / self.viscosity
		return Re, U, L

	def get_rectime(self,rec):
		initial_offset = 0 
		tplot = float(self.MD.header.tplot)
		delta_t = float(self.MD.header.delta_t)
		Nmass_ave = float(self.MD.header.Nmass_ave)
		t = (initial_offset+rec-0.5)*delta_t*Nmass_ave*tplot
		return t
