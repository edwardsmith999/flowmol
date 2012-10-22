#! /usr/bin/env python2.7
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
npxmdlist  = [2,4,6,8,10,12]
npycfdlist = [2,3,4,5,6,7,8]
for npx_md in npxmdlist:
	for npy_cfd in npycfdlist:
		job = RunClass( npx_md, npy_md, npz_md, npx_cfd, npy_cfd, npz_cfd )
		print(job)
		job.execute()
		job.concatenate()
		job.checkvalues()
