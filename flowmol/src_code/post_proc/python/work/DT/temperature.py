import os
import numpy as np
import matplotlib.pyplot as plt
from MDPlotData import MD_PlotData
from HeaderData import *

# Objects
Header = HeaderData(open('../../results/simulation_header','r'))
Getter = MD_PlotData('../../results/',cpol_bins=False)

# Parameters
#Nmass_ave = int(Header.Nmass_ave) 
#Lz = float(Header.globaldomain3) 
#r_ii = float(Header.r_ii)
#r_oi = float(Header.r_oi)
#r_io = float(Header.r_io)
#r_oo = float(Header.r_oo)
nbins = int(Header.gnbins1)*int(Header.gnbins2)*int(Header.gnbins3)
filebytes = os.path.getsize('../../results/vbins')
maxrec = filebytes / (3*8*nbins) # 8 is the size of 1 double precision 

X, Y, T = Getter.get_Tplane_splot_args(2,0,1,0,maxrec)

# PLOTS
fig = plt.figure()
ax = fig.add_subplot(111)
ax.contourf(X,Y,T)
plt.show()
# ...

# SHOW/SAVE
#ax.savefig('generic.png')
#plt.show()
