
from matplotlib.mlab import griddata
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


mp = np.genfromtxt('./results/macroscopic_properties',delimiter=';',names=True)

fig, axa = plt.subplots(2,2,figsize=(10,8))
plt.ion(); n=0
first_time = 0
xres = 51
yres = 51
zres = 51
axs = axa.ravel()

xmin = -7.; xmax = 7.
ymin = -7.; ymax = 7.
zmin = -7.; zmax = 7.

#for i in range(20,1000,10):
i = 40

filenames = []
#Plots with 30000 are potential energy
#filenames.append("./fort.30{:03d}".format(i))
#Plots with 40000 are virial
#filenames.append("./fort.40{:03d}".format(i))
#filenames.append("./fort.40{:03d}".format(i))
#Plots with 50000 are rijrij
filenames.append("./fort.50{:03d}".format(i))

for j, filename in enumerate(filenames):
    print(filename)
    a = np.genfromtxt(filename)
    x=a[:,0]; y=a[:,1]; z=a[:,2]

    if j is 2:
        u=a[:,4]
    else:
        u=a[:,3]

    #The reshape solution
    xi = x.reshape((xres,yres,zres))
    yi = y.reshape((xres,yres,zres))
    zi = z.reshape((xres,yres,zres))
    ui = u.reshape((xres,yres,zres))

    from mayavi import mlab
    mlab.contour3d(xi, yi, zi, ui)
    mlab.show()

