import matplotlib.pyplot as plt
from MDPlotData import *

Getter = MD_PlotData('../../results/',cpol_bins=True)
fig = plt.figure(figsize=(3,12))
ax = fig.add_subplot(111)

drec = 250
cnt = 0
for rec in range(0,250,drec):

	r, z, vr, vz = Getter.get_vplane_streamplot_args(1,0,2,rec,rec+drec)
	ax.streamplot(r,z,vr,vz)

	plt.savefig('streamplot.'+'%05d'%cnt+'.jpg',bbox_inches='tight')
	ax.cla()

	cnt += 1
