#! /usr/bin/env python
from os import system 

fobj    = open('procs','r')

npx_md  = int(fobj.readline()[0:4])
npy_md  = int(fobj.readline()[0:4])
npz_md  = int(fobj.readline()[0:4])
npx_cfd = int(fobj.readline()[0:4])
npy_cfd = int(fobj.readline()[0:4])
npz_cfd = int(fobj.readline()[0:4])

fobj.close()	

nproc = npx_md*npy_md*npz_md + npx_cfd*npy_cfd*npz_cfd
cmd = 'mpiexec -n ' + str(nproc) + ' ./a.out'

print(cmd)
system(cmd)

cmd = 'cat fort.10* > info_realms && rm fort.10*'
system(cmd)
cmd = 'cat fort.20* > info_olap_md && rm fort.20*'
system(cmd)
cmd = 'cat fort.3* > info_graph && rm fort.3*'
system(cmd)
cmd = 'cat fort.4* > info_MD_send && rm fort.4*'
system(cmd)
cmd = 'cat fort.5* > info_CFD_recv && rm fort.5*'
system(cmd)
cmd = 'cat fort.6* > info_MDCFD && rm fort.6*'
system(cmd)
