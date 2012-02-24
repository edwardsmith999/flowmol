#! /usr/bin/env python
###############################################################
# calc_etevtcf_fft.py
# author: David Trevelyan
#
# calc_etevtcf_fft.py calculates the end-to-end vector time
# correlation function for the end-to-end vectors recorded 
# throughout a polymer simulation. The file etev contains
# a time history of the end-to-end vectors of all chains in 
# the domain every TPLOT timesteps. 
# 
# The end-to-end vectors are read from the binary file into
# a large numpy array using numpy.fromfile(). These data are
# then rearranged into three arrays (one for each x,y and z
# direction): the time history of the x-componen of the end-
# to-end vector for chain i is stored in etev_x[i] (note we
# count from 0 in python, so chainIDs in the Fortran code are
# 1 greater than here). 
# 
# We seek the time-correlation of the inner product of end-to-
# end vectors, hence we simply add the auto-correlations of
# each cartesian component. The auto-correlation functions for
# every chain i are calculated, using the fast Fourier trans-
# form method bundled in scipy.signal.fftconvolve(). Note that
# we compute the correlation by convolving in reverse. 
#
# The ensemble average of the inner product auto-correlations,
# C(t), is computed and normalised to the range 0<=C(t)<=1.
# 
# USAGE:
#     - This program must be called from the /results folder.
#     - It outputs firstly (perhaps in bad practice), using
#       the stderr channel, a dialogue of calculation progress.
#       This is followed by a stdout output of C(t). So, for
#       example, to write C(t) to output file "end-to-end":
#
#            :$ ./vis/calc_etevtcf_fft.py > end-to-end
#       
#       The above command would allows the user to write C(t)
#       to a file whilst at the same time being able to 
#       track the progress of the program. To surpress the 
#       progress information, use:
#
#       :$ ./vis/calc_etevtcf_fft.py 1> end-to-end 2>null
#
# -------------------------------------------------------------
import os
import sys
import read_header
import numpy as np
import scipy.signal as sp

# Get important information from simulation_header with read_header.py
delta_t		= float(read_header.delta_t)
tplot		= int(read_header.tplot)
Nsteps		= int(read_header.Nsteps)
nchains     = int(read_header.nchains)


etev = np.fromfile('etev',dtype=float)                  # read from binary file
Ndt = len(etev)/(3*nchains)                             # find number of timesteps
etev_x = np.empty((nchains,Ndt),float)         
etev_y = np.empty((nchains,Ndt),float)
etev_z = np.empty((nchains,Ndt),float)

# begin rearrange of etev
for i in range(nchains):
	etev_x[i]  = etev[0+i::3*nchains]
	etev_y[i]  = etev[nchains+i::3*nchains]
	etev_z[i]  = etev[2*nchains+i::3*nchains]
del(etev)
# end rearrange

# calculate auto-correlation functions
C = np.zeros(Ndt)
for i in range(nchains):
	progress = 'Calculating auto-correlation function ' + str(i+1) + ' of ' + str(nchains) + '...\n'
	sys.stderr.write(progress)
	auto_cx = sp.fftconvolve(etev_x[i],etev_x[i][::-1]) # correlation is convolution in reverse
	auto_cy = sp.fftconvolve(etev_y[i],etev_y[i][::-1])	
	auto_cz = sp.fftconvolve(etev_z[i],etev_z[i][::-1])	
	l = len(auto_cx)	                                # only take positive half
	C += auto_cx[l/2:l]+auto_cy[l/2:l]+auto_cz[l/2:l]

C = C/float(nchains)                                    # ensemble average
C = C/C[0]                                              # normalise

for i in range(len(C)):					
	outstring = str(i*tplot*delta_t).rjust(32) + str(C[i]).rjust(32)
	print(outstring)
