#! /usr/bin/env python
###############################################################################
# read_vslice.py - PRINTS THE VELOCITY PROFILE ACROSS A NORMALISED DOMAIN
# author: David Trevelyan
#
# read_vslice.py reads all the 'snapshot' velocity profiles from the file
# 'vslice' and calculates a set of (Nvel_records/nt) mean profiles that show
# the time evolution of the velocity field. 
# 
# Time resolution can be easily adjusted by changing the numerical value of 
# "nt", i.e. adjusting how many values in the vslice array (created in this
# script) are passed to the function 'mean(values)'. Setting "nt" to 1 will
# result in a velocity profile that has been averaged over ALL the "snapshots"
# in results/vslice. Setting "nt" to 2 will result in two profiles: an average
# of the first half of the simulation is printed before the average over the
# second half (and so on...). If you do choose to adjust "nt", don't forget
# to change the loop that prepares the output string:
# 
# 	for t in ________________ :           (for example, in range(1:nt))
#		for cell in range(nbins):
#			...
#			...
#
# to something appropriate.
#
# The output print string is designed to be easily piped and is formatted as 
# a normalised position across the domain followed by the average velocity 
# for the bin centered at the corresponding 'unnormalised' position. 
# read_vslice is designed to be called from the post_proc folder.
#
# EXAMPLE (if v_outflag=2 and le_sd=1):
#
#		./python/read_vslice.py > vprofile
#		gnuplot
#	  > set xlabel "u_x"
#	  > set ylabel "y/Y"
#	  > plot [][0:1] "vprofile" index 0 using 2:1 with points,
#                    "vprofile" index 1 using 2:1 with points... etc.
#                                                                             
# -----------------------------------------------------------------------------
from numpy import *
from pylab import *
from array import array
import math
import os
import read_header
import read_mslice

def mean(values):                                 # Averages a list of values
	return float(sum(values)/len(values))

nd            = int(read_header.nd)               # Number of dimensions
Nsteps        = int(read_header.Nsteps)           # Number of timesteps
tplot         = int(read_header.tplot)            # Plotting frequency
Nvel_ave      = int(read_header.Nvel_ave)         # Averaging frequency
le_sd         = int(read_header.le_sd)-1          # L-E shear direction 
                                                  # (python counts from 0)
initialstep   = int(read_header.initialstep)      # Initial timestep
v_outflag     = int(read_header.velocity_outflag) # Slice direction
gnbins        =[int(read_header.gnbins1),         # Number of bins in each
                int(read_header.gnbins2),         # direction
                int(read_header.gnbins3)]
Nvel_records  = int(math.floor(
                   (Nsteps-initialstep)/
                   (tplot*Nvel_ave)))
nbins         = gnbins[v_outflag-1]    			  # Python counts from 0

mslice = read_mslice.mslice

f = open('../results/vslice','rb')                # File object f
vslice = array('d')                               # Array to store data
vslice.fromfile(f,Nvel_records*nd*nbins)          # Get vslice data
vslice = reshape(vslice,(Nvel_records,nd,nbins))  # Reshape array
f.close()

nt        = 16                  # Split vslice into nt time regions
vprofile  = [[0.0]*nbins]*nt    # Initialise array of average velocity profiles
sp_coords = [0.0]*nbins         # Initialise shear plane coordinates

for t in [0,1,3,7,15]:          # Plot profiles for time in powers of 2
	for cell in range(nbins):                     # Loop over nbins cells
		sp_coords[cell] = float(cell+0.5)/nbins   # Center of cell coordinates
		varray = vslice[t*(Nvel_records/nt)+1:(t+1)*(Nvel_records/nt),
		                le_sd,cell]               # Array of velocities for
		                                          # this cell over Nvel_records/
		                                          # nt timesteps
		marray = mslice[t*(Nvel_records/nt)+1:(t+1)*(Nvel_records/nt),
		                cell]                     # Total mass over Nvel_records
		                                          # /nt timsteps
		vprofile[t][cell] = mean(varray/marray)   # Average velocity profile
		outstring = (str(sp_coords[cell]).rjust(16)
		             + str(vprofile[t][cell]).rjust(32))
		print(outstring)                          # Write for gnuplot too
	print('\n')                                   # Gnuplot index separator
	plot(vprofile[t][:],sp_coords,'x')            # Matplotlib plot

ylabel('Shear plane coordinate / shear plane domain length')
xlabel('Velocity (LJU)')
show()                                            # Show plot
