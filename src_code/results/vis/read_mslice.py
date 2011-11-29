#! /usr/bin/env python
######################################################################################
# read_mslice.py
# author: David Trevelyan
#
# read_mslice.py reads the binary data from the file 'mslice' and stores it in 
# a numpy array. This can then be used as part of read_vslice to plot the average
# velocity profile from a simulation.
# 
# ------------------------------------------------------------------------------------
from numpy import *
from array import array
import read_header
import math
import os

# Retrieve important variables from simulation_header with read_header.py
Nsteps 	  		= int(read_header.Nsteps)
tplot  	  		= int(read_header.tplot)
Nmass_ave 		= int(read_header.Nmass_ave)
m_outflag		= int(read_header.mass_outflag)
globalnbins 	= [int(read_header.globalnbins1),int(read_header.globalnbins2),int(read_header.globalnbins3)]
Nmass_records 	= int(math.floor(Nsteps/(tplot*Nmass_ave)))
nbins			= int(globalnbins[m_outflag-1])

f = open('mslice','rb')								# Create binary file object
mslice = array('i')									# Initialise array of integers
mslice.fromfile(f,nbins*Nmass_records)				# Read entire binary file into array
mslice = reshape(mslice,(Nmass_records,nbins))		# Reshape array into 3 dimensions
f.close()											# Close file object
