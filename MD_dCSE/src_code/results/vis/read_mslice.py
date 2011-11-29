#! /usr/bin/env python
from numpy import *
from array import array
import read_header
import math
import os

Nsteps 	  		= int(read_header.Nsteps)
tplot  	  		= int(read_header.tplot)
Nmass_ave 		= int(read_header.Nmass_ave)
m_outflag		= int(read_header.mass_outflag)
globalnbins 	= [int(read_header.globalnbins1),int(read_header.globalnbins2),int(read_header.globalnbins3)]
Nmass_records 	= int(math.floor(Nsteps/(tplot*Nmass_ave)))
nbins			= int(globalnbins[m_outflag-1])

f = open('mslice','rb')
mslice = array('i')
mslice.fromfile(f,nbins*Nmass_records)
mslice = reshape(mslice,(Nmass_records,nbins))
f.close()
