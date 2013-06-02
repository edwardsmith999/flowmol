#! /usr/bin/env python 


from MDPlotData import MD_PlotData
import matplotlib.pyplot as plt
import numpy as np

PlotObj = MD_PlotData('./../../results/')

di = 3
cnt = 0


#X,Y,U,V = PlotObj.get_vplane_quiver_args(2,0,1,5,5)

#print(U)

#plt.contourf(X, Y, U)
#plt.show()

for i in range(di/2,1000,1):
	minrec = i - di/2 
	maxrec = i + di/2 
	print(minrec,maxrec)
	#X,Y,U,V = PlotObj.get_vplane_streamplot_args(2,0,1,minrec,maxrec)
	X,Y,U,V = PlotObj.get_vplane_quiver_args(2,0,1,minrec,maxrec)
	plt.contourf(X, Y, U)
	#plt.show()
	#plt.pcolormesh(X,Y,U)
	#speed = np.sqrt(U*U + V*V)
	#lw = 5*speed/speed.max()
	#plt.streamplot(X,Y,U,V,color=U,density=2.0,linewidth=lw,cmap=plt.cm.autumn)
	#plt.quiver(X,Y,U,V,color='k')

	fname = 'streamplot.' + "%05d"%cnt + '.png'
	plt.savefig(fname)
	plt.cla()
	cnt += 1
