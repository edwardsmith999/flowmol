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
from array import array
import numpy as np
import os
import struct
import read_header

# Get important information from simulation_header with read_header.py
delta_t		= float(read_header.delta_t)
tplot		= int(read_header.tplot)

f = open('r_gyration','rb')                       # Create file object
num_recs = int(os.path.getsize('r_gyration'))     # Find number of records
num_recs = (num_recs-4)/16 + 1                    # Bizarre fortran binary format
R_g = f.read(4)                                   # Padding at beginning

for i in range(num_recs):
	R_g = f.read(8)
	R = struct.unpack('@d',R_g)
	outstring = str(i*tplot*delta_t).rjust(32) + str(R[0]).rjust(30)
	print(outstring)
	R_g = f.read(8)                               # Padding between records

f.close()	
