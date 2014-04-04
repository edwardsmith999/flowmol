#! /usr/bin/env python2.7

#Opens VMD with overlayed fields
import os
import numpy as np
import csv
import sys

from misc_lib import Chdir
from MDFields import (MD_mField, MD_vField, MD_TField, 
                      MD_momField, MD_dField, MD_CVField)
from HeaderData import HeaderData
from WriteColorMap import WriteColorMap

class VMDFields:
    
    """
        Class to create reformat mass, momentun, temperature, etc 
        fields and run VMD with molecules coloured by these fields
    """

    def __init__(self,fieldobj,fdir='../MD_dCSE/src_code/results/'):
        self.fdir = fdir
        self.fieldobj = fieldobj

        #Read Header
        with open(fdir + 'simulation_header','r') as f:
            self.header = HeaderData(f)

        #Check simulation progress file to check if run 
        #is still going or has finished prematurely (crashed)
        with open(self.fdir+'simulation_progress', 'r') as f:
            simulation_time = f.readline()
        for i in dir(self.header):
            if (i.find('Nsteps') == 0):
                self.Nsteps = int(vars(self.header)[i])
                if (int(simulation_time) < self.Nsteps):
                    self.Nsteps = int(simulation_time)
                    self.finished = False
                    self.reformat_vmdtemp()
                else:
                    self.finished = True

        #Get VMD intervals from header
        self.starts = []; self.ends = []
        for i in dir(self.header):
            if (i.find('vmd_start') == 0):
                start = int(vars(self.header)[i])
                #print(start,int(simulation_time),start < int(simulation_time))
                if (start < int(simulation_time)):
                    self.starts.append(start)
            if (i.find('vmd_end') == 0):
                end = int(vars(self.header)[i])
                #print(end,int(simulation_time),end < int(simulation_time))
                if (end < float(simulation_time)):
                    self.ends.append(end)
        #If part way though an interval, set maximum to last iteration run
        if (len(self.starts) > len(self.ends)):
            self.ends.append(int(simulation_time))

        #Get averaging time per record 
        if (isinstance(fieldobj, MD_mField)):
            self.Nave = self.header.Nmass_ave
        elif (isinstance(fieldobj, MD_vField)):
            self.Nave = self.header.Nvel_ave
        elif (isinstance(fieldobj, MD_TField)):
            self.Nave = self.header.NTemp_ave
        elif (isinstance(fieldobj, MD_dField)):
            self.Nave = self.header.Nmass_ave
        elif (isinstance(fieldobj, MD_momField)):
            self.Nave = self.header.Nvel_ave
        elif (isinstance(fieldobj, MD_CVField)):
            self.Nave = self.header.Nvflux_ave

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

        #Write VMD intervals
        print('Writing VMD intervals data')
        self.starts.sort(); self.ends.sort()
        self.vmdintervals = zip(self.starts,self.ends)
        with open(self.vol_dir + '/vmd_intervals','w+') as f:
            for i in self.vmdintervals:
                f.write(str(i[0]) + '\n' + str(i[1]) + '\n')

    #Write range of dx files based on VMD intervals
    def write_dx_range(self,component=0,clims=None):
        #Clean previous files
        outdir = self.fdir+"./vmd/vol_data/"
        filelist = [ f for f in os.listdir(outdir + ".") if f.endswith(".dx") ]
        for f in filelist:
            os.remove(outdir+f)

        #Write range of files in all intervals
        print('Writing dx files intervals data',self.vmdintervals)
        clims_array = []
        for i in self.vmdintervals:
            fieldrecstart = i[0]/(int(self.header.tplot)*int(self.Nave))
            fieldrecend   = i[1]/(int(self.header.tplot)*int(self.Nave))
            #If limits are not specified, store time history for all intervals and average
            if (clims == None):
                clims_array.append(self.fieldobj.write_dx_file(fieldrecstart,fieldrecend,component=component))
            elif (len(clims) != 2):
                quit("Error in write_dx_range - clims should be tuple length 2 of form (cmin,cmax)")
            else:
                dummy = self.fieldobj.write_dx_file(fieldrecstart,fieldrecend,component=component)

        #Write maximum and minimum values for colourbar
        if (clims == None):
            clims = np.max(clims_array,axis=0)
        with open(self.vol_dir + '/colour_range','w+') as f:
            f.write(str(clims[0]) + '\n' + str(clims[1]) + '\n')

    def writecolormap(self,cmap='RdYlBu_r'):
        cmap_writer = WriteColorMap(cmap,1024)
        cmap_writer.write(self.vol_dir)

    def reformat_vmdtemp(self):

        """ 
            If run has not finished, attempt to build and
            run fortran code to reorder temp files to a
            useful form
        """

        #Get number of molecules from header file
        for i in dir(self.header):
            if (i.find('globalnp') == 0):
                np = int(vars(self.header)[i])

        # Build and call VMD_reformat with np from header
        with Chdir(fdir + '../debug_scripts/'):
            cmd = "ifort -O3 -o vmd_reformat.exe vmd_reformat.f90"
            os.system(cmd)
            cmd = './vmd_reformat.exe ' + str(np)
            os.system(cmd)

if __name__ == "__main__":

    fdir='../MD_dCSE/src_code/results/'

    fieldtypes = {'mbins','vbins','Tbins',
                  'density','momentum','CV_config',
                  'CV_kinetic','CV_total'}

    if(len(sys.argv) == 1):
        print("No field type specified, options include: " 
              + str(fieldtypes) + " Setting default vbins")
        objtype = 'vbins'
        component = 0
    else:
        objtype = sys.argv[1]
        if(len(sys.argv) == 2):
            print("No components direction specified, setting default = 0")
            component = 0
        else:
            component = sys.argv[2]

    if (objtype == 'mbins'):
        fobj = MD_mField(fdir)
    elif (objtype == 'vbins'):
        fobj = MD_vField(fdir)
    elif (objtype == 'Tbins'):
        fobj = MD_TField(fdir)
    elif (objtype == 'density'):
        fobj = MD_dField(fdir)
    elif (objtype == 'momentum'):
        fobj = MD_momField(fdir)
    elif (objtype == 'CV_config'):
        fobj = MD_CVField(fdir,fname='psurface')
    elif (objtype == 'CV_kinetic'):
        fobj = MD_CVField(fdir,fname='vflux')
    elif (objtype == 'CV_total'):
        fobj = MD_CVField(fdir,fname='total')
    else:
        quit("Argument " + sys.argv[1] + " is not a valid field type")

    vmdobj = VMDFields(fobj,fdir)
    vmdobj.write_vmd_header()
    vmdobj.write_vmd_intervals()
    vmdobj.write_dx_range(component=component)
    vmdobj.writecolormap('RdYlBu')
    with Chdir(fdir + './vmd/'):
        command = "vmd -e " + "./plot_MD_field.vmd"
        os.system(command)
