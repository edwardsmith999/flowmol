import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from MDPlotData import MD_PlotData
from HeaderData import HeaderData

# Objects
Header = HeaderData(open('../../results/simulation_header','r'))
Getter = MD_PlotData('../../results/',cpol_bins=True)

# Parameters
Nmass_ave = int(Header.Nmass_ave) 
Lz = float(Header.globaldomain3) 
r_ii = float(Header.r_ii)
r_oi = float(Header.r_oi)
r_io = float(Header.r_io)
r_oo = float(Header.r_oo)
nbins = int(Header.gnbins1)*int(Header.gnbins2)*int(Header.gnbins3)
filebytes = os.path.getsize('../../results/mbins')
maxrec = filebytes / (4*nbins) # 4 is the size of a single integer

# Read profile
r, mprofile = Getter.get_mslice_plot_args(0,0,maxrec)
r = r + r_oi # Correct to r_oi < r < r_io
mprofile = mprofile / Nmass_ave # Average per snapshot

dr = r[1] - r[0] # Assuming uniform
binvolumes = 2.0 * np.pi * r * dr * Lz # Cylindrical element

density = mprofile / binvolumes

# Data for gnuplot
gpdata = zip(density,r)
for point in gpdata:
	line = str(point[0]) + '\t' + str(point[1])
	print(line)


# Plot
mpl.rc('text',usetex=True)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('$\\rho$')
ax.set_ylabel('$r$')
ax.axhspan(r_ii,r_oi,color='gray')
ax.axhspan(r_io,r_oo,color='gray')
ax.plot(density,r,'k')
ax.plot(density,r,'ko')

# Savefig
plt.savefig('density.png',bbox_inches='tight')
# Show
plt.show()
