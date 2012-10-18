#! /usr/bin/env python
from os import system 
import check

fobj    = open('input','r')

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

cmd = 'cat fort.11* > info_MD_recv && rm fort.11*'
system(cmd)
cmd = 'cat fort.1* > info_realms && rm fort.1*'
system(cmd)
cmd = 'cat fort.2* > info_olap_md && rm fort.2*'
system(cmd)
cmd = 'cat fort.3* > info_graph && rm fort.3*'
system(cmd)
cmd = 'cat fort.4* > info_MD_send && rm fort.4*'
system(cmd)
cmd = 'cat fort.5* > info_CFD_recv && rm fort.5*'
system(cmd)
cmd = 'cat fort.6* > info_maps && rm fort.6*'
system(cmd)
cmd = 'cat fort.7* > info_scatter_md && rm fort.7*'
system(cmd)
cmd = 'cat fort.8* > info_gather_cfd && rm fort.8*'
system(cmd)
cmd = 'cat fort.9* > info_CFD_send && rm fort.9*'
system(cmd)


check.gatherscattervals('info_scatter_md')
check.gatherscattervals('info_gather_cfd')

check.gatherscattervals('info_CFD_send')
check.gatherscattervals('info_MD_send')

check.gatherscattervals('info_MD_recv')
check.gatherscattervals('info_CFD_recv')
