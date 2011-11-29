#! /usr/bin/env python
from numpy import *
from array import array
import math
import os
import read_header
import read_mslice

def mean(values):
	return float(sum(values)/len(values))

nd				= int(read_header.nd)
Nsteps 			= int(read_header.Nsteps)
tplot			= int(read_header.tplot)
Nvel_ave		= int(read_header.Nvel_ave)   
Nvel_records	= int(math.floor(Nsteps/(tplot*Nvel_ave)))
v_outflag		= int(read_header.velocity_outflag)
globalnbins 	= [int(read_header.globalnbins1),int(read_header.globalnbins2),int(read_header.globalnbins3)]
shear_direction = int(read_header.shear_direction)-1
nbins			= globalnbins[v_outflag-1]			# Fortran starts counting at 1, python counts from 0

mslice = read_mslice.mslice

f = open('vslice','rb')								# Create file object f to read from vslice
#dble_filesize = int(os.path.getsize('vslice')/8)	# Calculate number of double precision records
vslice = array('d')									# Create vslice array in which to store data
vslice.fromfile(f,Nvel_records*nd*nbins)					# Store vslice data
vslice = reshape(vslice,(Nvel_records,nd,nbins)) 	# Reshape array into number of records, dimensions and bins
f.close()

for cell in range(nbins):
	outstring = str(float(cell+0.5)/nbins).rjust(16) + str(mean(vslice[Nvel_records/2:,shear_direction,cell]/mslice[Nvel_records/2:,cell])).rjust(32)
	print(outstring)
