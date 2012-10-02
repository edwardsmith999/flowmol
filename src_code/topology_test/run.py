#! /usr/bin/env python
from os import system 

fobj    = open('procs','r')

npx_md  = int(fobj.readline()[0])
npy_md  = int(fobj.readline()[0])
npz_md  = int(fobj.readline()[0])
npx_cfd = int(fobj.readline()[0])
npy_cfd = int(fobj.readline()[0])
npz_cfd = int(fobj.readline()[0])

fobj.close()	

nproc = npx_md*npy_md*npz_md + npx_cfd*npy_cfd*npz_cfd

cmd = 'mpiexec -n ' + str(nproc) + ' ./a.out'
system(cmd)
