import numpy as np
import matplotlib.pyplot as plt
from MDPlotData import MD_PlotData

Getter = MD_PlotData('../../results/')
T, Z, VR = Getter.get_vplane_splot_args(0,1,2,0,0,200)
plt.contourf(T,Z,VR)
plt.show()
