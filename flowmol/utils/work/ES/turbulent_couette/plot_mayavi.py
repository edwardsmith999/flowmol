from mayavi.mlab import quiver3d, figure, colorbar
import numpy as np

b = np.genfromtxt('./fort.6000')
f = figure(bgcolor=(1.,1.,1.),fgcolor=(0.,0.,0.))
quiver3d(b[:,2],b[:,3],b[:,4],b[:,5],b[:,6],b[:,7], line_width=1,mode='arrow',scale_mode='vector')
colorbar()
#quiver3d(b[:,2],b[:,3],b[:,4],b[:,8],b[:,9],b[:,10], line_width=1,mode='cone',scale_mode='vector')

