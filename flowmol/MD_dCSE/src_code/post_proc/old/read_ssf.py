#! /usr/bin/env python
##############################################################################
# read_ssf.py
# author: David Trevelyan
#
# desription:
#     read_struct.py reads the discretised static structure factor
#     from the binary file "stat_struct" in the results folder.
#     "stat_struct" consists of double precision values for S(k) over 
#     k-space k=i*dk (i=1,2...nmax), where dk = 2*pi/domain. 
#_____________________________________________________________________________
import os
import math
import struct
import read_header

domain    = [0.0]*3
domain[0] = float(read_header.globaldomain1)
domain[1] = float(read_header.globaldomain2)
domain[2] = float(read_header.globaldomain3)

dk        = [0.0]*3
dk[0]     = 2.0*math.pi/domain[0] 
dk[1]     = 2.0*math.pi/domain[1] 
dk[2]     = 2.0*math.pi/domain[2] 
nmax      = int(read_header.ssf_nmax)
nk        = nmax*2 + 1
axis1     = int(read_header.ssf_ax1) - 1 # fortran counts from 1
axis2     = int(read_header.ssf_ax2) - 1 # fortran counts from 1

fname = '../results/ssf'
f = open(fname,'rb')

print('ky,kx,S(kx,ky)')

for pos2 in range(nk): 
	n2 = pos2 - nmax
	k2 = n2*dk[axis2]
	for pos1 in range(nk):
		n1 = pos1 - nmax
		k1 = n1*dk[axis1]
		data = f.read(8)
		S    = struct.unpack('@d',data)
		outstring = str(k2).rjust(32) + str(k1).rjust(32) + str(S[0]).rjust(32)
		print(outstring)
	print('')

f.close()
