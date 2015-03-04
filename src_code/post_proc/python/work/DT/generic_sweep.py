import os
import numpy as np
import matplotlib.pyplot as plt
from MDPlotData import MD_PlotData
from HeaderData import *

cpol_bins = True
drec = 20
figname = 'generic'
loc = '../../results/'
inspectfile = 'vbins'
inspectbytesperbin = 3*8 # 8 is the size of 1 double precision 

# Global Objects
fig = plt.figure()
Header = HeaderData(open(loc+'simulation_header','r'))
Getter = MD_PlotData(loc,cpol_bins=cpol_bins)

def plot_and_save(rmin,rmax,cnt):

    # New axes on cleared figure
    ax = fig.add_subplot(111)

    # YOUR PLOT CODE GOES HERE 

    # Save and clear
    plt.savefig(figname+'.'+"%05d"%rmax+'.png')
    plt.clf()

# ============================================================================
# Useful Parameters
nbins = int(Header.gnbins1)*int(Header.gnbins2)*int(Header.gnbins3)
filebytes = os.path.getsize(loc+inspectfile)
maxrec = filebytes / (inspectbytesperbin*nbins) 
if (cpol_bins):
    r_oi = float(Header.r_oi)

# ============================================================================
# SWEEP THROUGH PLOTS
maxrec = maxrec - maxrec%drec
cnt = 0
for rec in range(drec/2,maxrec,drec):
    rmin = rec - drec/2
    rmax = rec + drec/2
    plot_and_save(rmin,rmax,cnt)
    cnt += 1
