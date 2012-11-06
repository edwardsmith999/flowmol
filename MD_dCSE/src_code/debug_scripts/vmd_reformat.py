#! /usr/bin/env python

# Gets number of molecules from header file and 
# calls fortran routine vmd_reformat.exe to convert 
# vmd.temp into readbale files
# Written by David Trevelyan

import os

headerfile = './../results/simulation_header'
# Extract np from header
fobj = open(headerfile,'r')
np = 0
while np==0:
	line = fobj.readline().split(';')
	if (line[1].strip() == 'globalnp'):
		np   = int(line[2].strip())
	
# Either call VMD_reformat with np
cmd = './vmd_reformat.exe ' + str(np)
os.system(cmd)
