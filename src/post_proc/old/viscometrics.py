#! /usr/bin/env python
# ----------------------------------------------------------------------------
# viscometrics.py
# author: David Trevelyan 
# 
# viscometrics.py results/pvirial and calculates the (time-averaged) shear 
# viscosity, first and second normal stress differences and the hydrostatic 
# pressure. A histogram of the shear viscosity records is plotted by 
# uncommenting the line "show()". Regardless, the data for plotting the
# histogram is written to /post_proc/viscosity_histogram.
# ============================================================================
import numpy as np
from pylab import *
from array import array
import os
import read_header

x     = int(read_header.le_sd)-1                  # Python counts from zero
y     = int(read_header.le_sp)-1                  # Python counts from zero
z     = int(read_header.le_rp)-1                  # Python counts from zero
le_sr = float(read_header.le_sr)
 
stressfile = '../results/pvirial'

f = open(stressfile,'rb')                         # Binary file object
Ndoubles = os.path.getsize(stressfile)/8          # Number of doubles in file
Nrecords = Ndoubles/9                             # Number of Pxy tensors

data = array('d')                                 
data.fromfile(f,Ndoubles)                         # Read data from file
f.close()                                         # Close file

Pxy = zeros((Nrecords,3,3))                       # Create array of stress
rec = 0                                           # tensors
for i in range(0,Ndoubles,9):                     # Extract correct data
	Pxy[rec,x,x] = data[i+0] 
	Pxy[rec,y,x] = data[i+1] 
	Pxy[rec,z,x] = data[i+2] 
	Pxy[rec,x,y] = data[i+3] 
	Pxy[rec,y,y] = data[i+4] 
	Pxy[rec,z,y] = data[i+5] 
	Pxy[rec,x,z] = data[i+6] 
	Pxy[rec,y,z] = data[i+7] 
	Pxy[rec,z,z] = data[i+8] 
	rec += 1
                                                  # ___ Create arrays ________
eta = zeros(Nrecords)                             # Array of viscosities
N1  = zeros(Nrecords)                             # Arrays of first and second
N2  = zeros(Nrecords)                             # stress differences
P   = zeros(Nrecords)                             # Array of pressures
                                                  # --------------------------
for rec in range(Nrecords):
	eta[rec] = -(Pxy[rec,x,y]/le_sr)              # Calculate quantities
	N1[rec]  = -(Pxy[rec,x,x] - Pxy[rec,y,y])
	N2[rec]  = -(Pxy[rec,y,y] - Pxy[rec,z,z])
	P[rec]   = np.mean([Pxy[rec,x,x],Pxy[rec,y,y],Pxy[rec,z,z]])

outstring = ''                                    # Prepare outstrings
outstring += 'Average viscometric data: \n\n'
outstring += '\tShear viscosity: ' + str(np.mean(eta)) + '\n'
outstring += '\tFirst normal stress difference: ' + str(np.mean(N1)) + '\n'
outstring += '\tSecond normal stress difference: ' + str(np.mean(N2)) + '\n'
outstring += '\tHydrostatic pressure: ' + str(np.mean(P)) + '\n'
print(outstring)

mu = np.mean(eta)                                 # Mean viscosity
sig = np.std(eta)                                 # Std dev viscosity
histogram = hist(eta,bins=50,range=(mu-4*sig,mu+4*sig))
#show()                                           # Uncomment to plot hist

etafile = './viscosity_histogram'
f = open(etafile,'w')
for i in range(len(histogram[0])):                # Write histogram data
	x = str(round(histogram[1][i],3)).rjust(16)
	y = str(histogram[0][i]).rjust(16)
	f.write(x+y+'\n')
f.close()

