#! /usr/bin/env python
###############################################################################
# read_header.py
# author: David Trevelyan
#
# This script is designed for importing into other python routines that require
# information from the file 'simulation_header'
#
# -----------------------------------------------------------------------------

f = open('./../../results/simulation_header','r')    # Create ascii file object
headerdata = f.readlines()                      # Read all lines into list
for i in range(len(headerdata)):                # Loop all lines 
	headerdata[i] = headerdata[i].split(';')    # Split list elements into
	                                            # sub-lists by semicolon 
	for j in range(len(headerdata[i])):         # Loop through sublist
		headerdata[i][j] = (headerdata[i][j].
		                    strip().replace(
		                    '(','').replace(
		                    ')',''))            # Strip white space and ()
	vars()[headerdata[i][1]] = headerdata[i][2] # Create variables with name 
	                                            # headerdata[i][1] and value
	                                            # headerdata[i][2]
f.close()

f = open('./../../results/simulation_progress','r')  # Create ascii file object
Nsteps = f.read().strip()                       # Read number of timesteps 
f.close()
