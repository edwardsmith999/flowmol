from MDPlotData import MD_PlotData
import matplotlib.pyplot as plt
import numpy as np

PlotObj = MD_PlotData('../../../../coupler_dCSE/src_code/md_data/results/')

di = 5
cnt = 0
for i in range(di/2,1000,1):
	minrec = i - di/2 
	maxrec = i + di/2 
	print(minrec,maxrec)
	X,Y,U,V = PlotObj.get_vplane_streamplot_args(2,0,1,minrec,maxrec)
	#plt.pcolormesh(X,Y,U)
	speed = np.sqrt(U*U + V*V)
	lw = 5*speed/speed.max()
	plt.streamplot(X,Y,U,V,color=U,linewidth=lw,cmap=plt.cm.autumn)
	#plt.quiver(X,Y,U,V,color='k')

	fname = 'streamplot.' + "%05d"%cnt + '.png'
	plt.savefig(fname)
	plt.cla()
	cnt += 1
