#! /usr/bin/env python2.7

#Opens VMD with overlayed fields
import os
import numpy as np
import csv
import sys

from misc_lib import Chdir
from MDFields import MD_mField, MD_vField, MD_TField, MD_momField
from HeaderData import HeaderData

class VMDFields:

    def __init__(self,fieldobj,fdir='../MD_dCSE/src_code/results/'):
        self.fdir = fdir
        self.fieldobj = fieldobj

        #Read Header
        with open(fdir + 'simulation_header','r') as f:
            self.header = HeaderData(f)

        #Get averaging time per record 
        if (isinstance(fieldobj, MD_mField)):
            self.Nave = self.header.Nmass_ave
        elif (isinstance(fieldobj, MD_vField)):
            self.Nave = self.header.Nvel_ave
        elif (isinstance(fieldobj, MD_TField)):
            self.Nave = self.header.NTemp_ave
        elif (isinstance(fieldobj, MD_momField)):
            self.Nave = self.header.Nvel_ave

        #Create VMD vol_data folder
        self.vol_dir = fdir + './vmd/vol_data/'
        if not os.path.exists(self.vol_dir):
            os.makedirs(self.vol_dir)

    def write_vmd_header(self):

        #Write VMD intervals
        print('Writing VMD Header')
        with open(self.vol_dir + '/vmd_header','w+') as f:
            f.write(self.header.tplot + '\n')
            f.write(self.header.delta_t + '\n')
            f.write(self.Nave + '\n')

    def write_vmd_intervals(self):

        #Get VMD intervals from header
        starts = []; ends = []
        for i in dir(self.header):
            if (i.find('vmd_start') == 0):
                starts.append(int(vars(self.header)[i]))
            if (i.find('vmd_end') == 0):
                ends.append(int(vars(self.header)[i]))

        #Write VMD intervals
        print('Writing VMD intervals data')
        starts.sort(); ends.sort()
        self.vmdintervals = zip(starts,ends)
        with open(self.vol_dir + '/vmd_intervals','w+') as f:
            for i in self.vmdintervals:
                f.write(str(i[0]) + '\n' + str(i[1]) + '\n')

    #Write range of dx files based on VMD intervals
    def write_dx_range(self):
        print('Writing dx files intervals data')
        for i in self.vmdintervals:
            fieldrecstart = i[0]/(int(self.header.tplot)*int(self.Nave))
            fieldrecend   = i[1]/(int(self.header.tplot)*int(self.Nave))
            self.fieldobj.write_dx_file(fieldrecstart,fieldrecend)

            #for rec in range(fieldrecstart,fieldrecend):

#                 field = self.fieldobj.read(rec,rec)
#                 self.write_dx_file(field[:,:,:,:,0],rec)

#     def write_dx_file(self,field,rec):

#         """
#            Write MD field to dx file format which is primarily
#            useful for importing into VMD, see 
#            http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dxplugin.html
#            and format website http://www.opendx.org/index2.php
#         """

#         #Save file
#         dxFileName = self.vol_dir + '/MD' + str(rec) + '.dx'
#         eps = 0.5 #VMD seems to leave gaps round the edges
#         with open(dxFileName,'w+') as f:
#             Nx = float(self.header.gnbins1)
#             Ny = float(self.header.gnbins2)
#             Nz = float(self.header.gnbins3)
#             dx = (float(self.header.globaldomain1)+2.0*eps)/(float(self.header.gnbins1)-1.0)
#             dy = (float(self.header.globaldomain2)+2.0*eps)/(float(self.header.gnbins2)-1.0)
#             dz = (float(self.header.globaldomain3)+2.0*eps)/(float(self.header.gnbins3)-1.0)
#             originx = -float(self.header.globaldomain1)/2.0-eps
#             originy = -float(self.header.globaldomain2)/2.0-eps
#             originz = -float(self.header.globaldomain3)/2.0-eps

#             # - - Write Header - -
#             f.write("object 1 class gridpositions counts%8.0f%8.0f%8.0f\n" % (Nx,Ny,Nz))
#             f.write("origin%16g%16g%16g\n" % (originx,originy,originz))
#             f.write("delta %16g 0 0\n" % dx)
#             f.write("delta 0 %16g 0\n" % dy)
#             f.write("delta 0 0 %16g\n" % dz)
#             f.write("object 2 class gridconnections counts%8.0f%8.0f%8.0f\n" % (Nx,Ny,Nz))
#             f.write("object 3 class array type double rank 0 items%8.0f follows\n" % (Nx*Ny*Nz))

#             # - - Write Data - -
#             import matplotlib.pyplot as plt
#             plt.contourf(field[:,:,3,0])
#             plt.show()
#             col=1
#             for j in range(0,int(Ny)):
#                 for i in range(0,int(Nx)):
#                     for k in range(0,int(Nz)):
#                         f.write("%16E" % field[i,j,k,0])
#                         col=col+1
#                         if (col>3):
#                             f.write(' \n')
#                             col=1

#             # - - Write Footer - - 
#             if (col != 1):
#                 f.write('           \n')
#             f.write('object "Untitled" call field \n')

if __name__ == "__main__":

    fdir='../MD_dCSE/src_code/results/'

    fieldtypes = {'mbins','vbins','Tbins'}

    if(len(sys.argv) == 1):
        print("No field type specified, options include: " 
              + str(fieldtypes) + " Setting default vbins")
        objtype = 'vbins'
    else:
        objtype = sys.argv[1]

    if (objtype == 'mbins'):
        fobj = MD_mField(fdir)
    elif (objtype == 'vbins'):
        fobj = MD_vField(fdir)
    elif (objtype == 'Tbins'):
        fobj = MD_TField(fdir)
    elif (objtype == 'momentum'):
		fobj = MD_momField(fdir)
    else:
        quit("Argument " + sys.argv[1] + " is not a valid field type")

    vmdobj = VMDFields(fobj,fdir)
    vmdobj.write_vmd_header()
    vmdobj.write_vmd_intervals()
    vmdobj.write_dx_range()
    with Chdir(fdir + './vmd/'):
        command = "vmd -e " + "./plot_MD_field.vmd"
        os.system(command)
