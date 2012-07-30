#! /usr/bin/env python
##############################################################################
# read_rdf.py
# author: David Trevelyan
#
# desription:
#     read_rdf.py reads the discretised radial distribution function
#     from the binary file "radialdist" in the results folder.
#     "radialdist" consists of double precision values for g(r) over the
#     distance 0 < r < rmax, split into nbins regions of width dr. 
#_____________________________________________________________________________
import os
import struct
import read_header

nbins = int(read_header.rdf_nbins)
rmax  = float(read_header.rdf_rmax)
dr    = rmax/float(nbins)

fname = '../results/rdf'
f = open(fname,'rb')

for i in range(nbins):
	data = f.read(8)
	g    = struct.unpack('@d',data)
	outstring = str((i+0.5)*dr).rjust(32) + str(g[0]).rjust(32)
	print(outstring)

f.close()
