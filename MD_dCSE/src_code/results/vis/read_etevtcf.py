#! /usr/bin/env python
from numpy import *
from array import array
import os
import read_header

delta_t		= float(read_header.delta_t)
tplot		= int(read_header.tplot)
Nsteps		= int(read_header.Nsteps)

f = open('etevtcf','rb')
dble_filesize = int(os.path.getsize('etevtcf')/8)
etevtcf = array('d')
etevtcf.fromfile(f,dble_filesize)
f.close()

for i in range(len(etevtcf)):
	outstring = str(i*tplot*delta_t).rjust(32) + str(etevtcf[i]).rjust(32)
	print(outstring)
