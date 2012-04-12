#! /usr/bin/env python
########################################################################
# read_r_gyration.py
# author: David Trevelyan
#
# read_r_gyration.py reads binary data from the file 'r_gyration' and
# sequentially prints each value on a new line. Each value is the
# ensemble average radius of gyration calculated by MD_dCSE. The time
# at which the radius of gyration is recorded is printed in column 1,  
# and the radius of gyration itself in column 2. 
#
# ----------------------------------------------------------------------
import os
import struct
import read_header

# Get important information from simulation_header with read_header.py
delta_t		= float(read_header.delta_t)
tplot		= int(read_header.tplot)


f = open('r_gyration','rb')                       # Create file object
num_recs = int(os.path.getsize('r_gyration')/8)   # Find number of records

for i in range(num_recs):
	R_g       = f.read(8)                         # Read binary data
	R         = struct.unpack('@d',R_g)           # Unpack as double precision number
	outstring = str(i*tplot*delta_t).rjust(32) + str(R[0]).rjust(30)
	print(outstring)

f.close()	
