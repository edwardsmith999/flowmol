#! /usr/bin/env python
from testmodule import *
import os

# ------------------------------------------------------------------------- #

def set_defaults():

	global npx_md
	global npy_md
	global npz_md
	global npx_cfd
	global npy_cfd
	global npz_cfd

	npx_md = 8
	npy_md = 4
	npz_md = 1

	npx_cfd = 2
	npy_cfd = 2
	npz_cfd = 1

# ------------------------------------------------------------------------- #

# Clean and/or setup logs directory
if (os.path.exists('./logs')):
	os.system('rm -r ./logs')
os.mkdir('./logs')

# Set default values of inputs
set_defaults()

# Parameter study
npxmdlist = [2,4,8,16]
for NPXMDTEST in npxmdlist:
	job = RunClass( NPXMDTEST, npy_md, npz_md, npx_cfd, npy_cfd, npz_cfd )
	print(job)
	job.execute()
	job.concatenate()
	job.checkvalues()
