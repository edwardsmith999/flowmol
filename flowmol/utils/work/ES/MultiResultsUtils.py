#! /usr/bin/env/ python
import os
import re
import sys
import shutil as sh
import numpy as np

from ResultsUtils import *
from MDFields import *
from HeaderData import *
from MD_PostProc import MD_PostProc

fdir_list = []
results_dir_base = '/home/es205/scratch/Re400/'

# BINS 64 x 256 x 64
#fdir_list.append(results_dir_base+'iter0_to_90870/')
fdir_list.append(results_dir_base+'iter90870_to_390870/')
fdir_list.append(results_dir_base+'iter390870_to_667054/')
fdir_list.append(results_dir_base+'iter667054_to_967054/')
fdir_list.append(results_dir_base+'iter967054_to_1267054/')
fdir_list.append(results_dir_base+'iter1267054_to_1567054/')

# BINS 84 x 198 x 50 with CV stats included (at eratic intervals in time!?)
fdir_list.append(results_dir_base+'iter1567054_to_1918000/')
fdir_list.append(results_dir_base+'iter1918000_to_2233899/')
fdir_list.append(results_dir_base+'iter2233899_to_2500000/')
fdir_list.append(results_dir_base+'iter2500000_to_2800000/')
fdir_list.append(results_dir_base+'iter2800000_to_3100005/')
fdir_list.append(results_dir_base+'iter3100005_to_3734274/')
fdir_list.append(results_dir_base+'iter3734274_to_4384274/')
fdir_list.append(results_dir_base+'iter4384274_to_5000000/')
fdir_list.append(results_dir_base+'iter5000000_to_5600000/')

outdir = '/home/es205/scratch/Re400/COMBINED_iter0_to_5600000/'
RU = ResultsUtils(outdir)
finalrec = 56
for fdir in fdir_list:
    fielddict = MD_PostProc(fdir)
    #get initialrec for each folder (With first runs at high resolution)
    #initialrec = int(finalrec)
    for keys, field in fielddict.plotlist.items():
        if keys in ['mbins','vbins','Tbins']:
            #Get filepath and recordsize
            filepath = fdir + keys
            if field.Raw.dtype == 'i':
                typesize = 4
            elif field.Raw.dtype == 'd':
                typesize = 8
            else:
                quit("Error - datatype not recognised")
            recbytes = (typesize*field.Raw.nperbin*np.product(field.Raw.nbins))
            if ((field.Raw.get_maxrec()+1)*recbytes != os.path.getsize(filepath)):
                print(filepath,(field.Raw.get_maxrec()+1)*recbytes != os.path.getsize(filepath))
                quit("File sizes do not match")

            #Get starting and final record of data in folder
            initialstep = int(field.Raw.header.initialstep)
            with open(fdir+'simulation_progress', 'r') as f:
                finalrec_progress = int(f.readline())-initialstep
                if ( float(field.Raw.header.Nsteps)
                    -float(finalrec_progress)
                    >float(field.Raw.header.tplot)):
                    finalstep = initialstep + int(finalrec_progress)
                    finished = False
                    print("Simulation not finished! = ", finalstep)
                elif ((int(finalrec_progress) == int(field.Raw.header.Nsteps)) |
                     ( float(field.Raw.header.Nsteps)
                      -float(finalrec_progress)
                      <float(field.Raw.header.tplot))):
                    finalstep = initialstep+ int(field.Raw.header.Nsteps)
                    finished = True
                    #print("Simulation finished = ", finalstep)
                else:
                    quit("Error -- simulation progress thinks run is longer than specified number of Nsteps!!!")

            plotfreq = float(field.plotfreq) * float(field.Raw.header.tplot)

            initialrec = np.floor( float(initialstep)/float(plotfreq) )
            finalrec = np.floor( float(finalstep)/float(plotfreq) )
            print(finished,"Writing files",keys," in ", fdir, "with starting recno = ",int(initialrec), "To final recno = ",int(initialrec)+os.path.getsize(filepath)/recbytes, "into folder = ", RU.outdir)

#            print("dir iter = ", initialstep, " to ",finalstep, 
#                  "dir start rec=", int(initialrec),"predicted finish = ", int(finalrec)-1,
#                  "dir no. records = ",os.path.getsize(filepath)/float(recbytes), 
#                  "final record no. = ", int(initialrec)+os.path.getsize(filepath)/float(recbytes))

           

            #RU.DismemberResults(fdir,keys, recbytes,initialrec=int(initialrec))
    



