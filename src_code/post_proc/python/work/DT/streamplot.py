import matplotlib.pyplot as plt
from MDPlotData import *

Getter = MD_PlotData('../../results/',cpol_bins=True)

r, z, vr, vz = Getter.get_vplane_streamplot_args(1,0,2,0,160)
plt.streamplot(r,z,vr,vz)
plt.show()
