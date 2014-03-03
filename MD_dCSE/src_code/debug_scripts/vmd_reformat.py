#! /usr/bin/env python2.7

# Gets number of molecules from header file and 
# calls fortran routine vmd_reformat.exe to convert 
# vmd.temp into readable files
# Written by David Trevelyan

import os

class VmdReformat:

    def __init__(self,fdir):

        headerfile = fdir + 'simulation_header'
        # Extract np from header
        fobj = open(headerfile,'r')
        np = 0
        while np==0:
            line = fobj.readline().split(';')
            if (line[1].strip() == 'globalnp'):
                np   = int(line[2].strip())
            
        # Build and call VMD_reformat with np from header
        os.system('ifort -O3 -o vmd_reformat.exe vmd_reformat.f90')
        cmd = './vmd_reformat.exe ' + str(np)
        os.system(cmd)

if __name__ == "__main__": 

    filepath = './../results/'

    VMDobj = VmdReformat(filepath)

    try:
        with open('./' + filepath + '/vmd_out.dcd'): pass
        os.system('vmd ' + filepath + '/vmd_out.dcd')
    except IOError:
        print 'vmd_out.dcd file missing, reorder unsuccessful'

