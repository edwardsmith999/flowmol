#! /usr/bin/env python
import numpy as np 
from HeaderData import *
from MDData import *
from CFDData import *
from CouetteAnalytical import *

class CPLData:

	def __init__(self,fdir):
		fobj            = open(fdir+'/results/coupler_header','r')
		self.header     = HeaderData(fobj)
		self.MD         = MDData(fdir+'/md_data/results/')
		self.CFD        = CFDData(fdir+'/couette_data/')
		self.yL         = self.get_coupleddomain()
		self.CFDyspace  = self.get_CFDyspace()
		self.MDcnstbins = self.get_constrainedMDbins()

		Re,U,L,wall     = self.get_Couetteparams()
		self.analytical = CouetteAnalytical(Re,U,L,wall)

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
		yspace = np.linspace(CFDbot,CFDtop,num=self.CFD.nry)
		return yspace	

	def get_Couetteparams(self):
		viscosity = 1.6	
		U = 1.0 
		wall = float(self.MD.header.tethdistbot2)
		L = self.yL - wall
		Ln = (self.yL - wall) / self.yL
		Re = float(self.header.density_cfd) * U * Ln / viscosity
		return Re, U, L, wall

	def get_rectime(self,rec):
		initial_offset = -1 
		tplot = float(self.MD.header.tplot)
		delta_t = float(self.MD.header.delta_t)
		Nmass_ave = float(self.MD.header.Nmass_ave)
		t = (initial_offset+rec-0.5)*delta_t*Nmass_ave*tplot
		return t
		
	def get_constrainedMDbins(self):
		cmin = int(self.header.jcmin_cnst)-1+1 # -1 for Fortran index, 
		cmax = int(self.header.jcmax_cnst)-1+1 # +1 for halo
		dy   = self.CFD.dy
		ymin = self.CFDyspace[cmin] - 0.5*dy 
		ymax = self.CFDyspace[cmax] + 0.5*dy 

		MDbinsize = float(self.MD.header.binsize2)
		const_bins = []
		MDbins = self.MD.bincenters  
		for binID in range(len(MDbins)):
			bincen = MDbins[binID]
			binbot = bincen - 0.5*MDbinsize
			bintop = bincen + 0.5*MDbinsize
			if (bintop <= ymax and bintop > ymin or
			    binbot <= ymax and binbot > ymin):
				const_bins.append(binID)		 

		return const_bins 
