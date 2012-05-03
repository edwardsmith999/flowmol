#! /usr/bin/env python
###############################################################################
# read_etevtcf.py
# author: David Trevelyan
# 
# read_etevtcf reads the end-to-end vector time correlation function data
# from the binary file 'etevtcf'. This routine is designed to be called
# from the results directory, and the output can be easily piped to a file
# of the user's choice. 
#
# EXAMPLE:
# 
#		./vis/read_etevtcf.py > etevtcf_data
#		gnuplot
#	  > set log x
#	  >	plot [0.001:10000][0:1] "./etevtcf_data" u 1:2 w p
#
# -----------------------------------------------------------------------------
from numpy import *
from array import array
import os
import read_header

# Get important information from simulation_header with read_header.py
delta_t = float(read_header.delta_t)
tplot   = int(read_header.tplot)
Nsteps  = int(read_header.Nsteps)

fileloc = '../results/etevtcf'

f = open(fileloc,'rb')                          # Create binary file object
dble_filesize = int(os.path.getsize(fileloc)/8) # Number of doubles
etevtcf = array('d')                            # Init. array of type double
etevtcf.fromfile(f,dble_filesize)               # Read entirety of f into array
f.close()                                       # Close file object

# Print simulation time and corresponding etevtcf columns
for i in range(len(etevtcf)):					
	outstring = str(i*tplot*delta_t).rjust(32) + str(etevtcf[i]).rjust(32)
	print(outstring)
