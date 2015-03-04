import matplotlib.pyplot as plt
from MDPlotData import *

Getter = MD_PlotData('../../results/',cpol_bins=True)
r, vprofile = Getter.get_vslice_plot_args(0,0,160)
plt.plot(vprofile,r,'o')
plt.show()
