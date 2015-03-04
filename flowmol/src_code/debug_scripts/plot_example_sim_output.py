#!/usr/bin/env python

import pylab as pl
import matplotlib.cm as cm
from numpy import *

a=open('./../results/try_sim_output', 'r').readlines()
a=a[0]
a=a.replace(')','')
a=a.replace('\n','')
a=a.replace(' ','')
a=a.replace('[','')
a=a.replace(']','')
a=a.split('(')
density = zeros(len(a))
temperature = zeros(len(a))
pressure = zeros(len(a))
stdpressure = zeros(len(a))
c = zeros((len(a),4))
n = 0
for b in a:
	b=b.replace('\n','')
	b = b.split(',')
	if (b[0] == ''):
		print b[0]	
	else:
		density[n] = float(b[0])
		temperature[n] = float(b[1])
		pressure[n] = float(b[2])
		stdpressure[n] = float(b[3])
		n = n + 1
		c[n,0:4] = [float(b[0]),float(b[1]),float(b[2]),float(b[3])]
		print str(c[n,0]) + " " + str(c[n,1]) + " " + str(c[n,2]) + " " + str(c[n,3])

#Plot
Ntemp = 15
fig = pl.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel('Temperature')
ax1.set_ylabel('Pressure')
for i in range(0,5):
	n = 1+Ntemp + 2*Ntemp*i
	m = n+Ntemp
	print c[n:m,0]
	colors = cm.rainbow(linspace(0, 1, 5))
	pl.scatter(c[n:m,1],c[n:m,2],c=colors[i])
ax1.set_xlim(0.0, 1.7)
#ax1.set_ylim(-2.0, 3.0)
pl.show()


#Plot
#fig = pl.figure()
#ax1 = fig.add_subplot(111)
#ax1.set_xlabel('Density')
#ax1.set_ylabel('Pressure')
#ax1.set_xlim(-0.1, 1.1)
#ax1.set_ylim(-3.0, 4.0)

#Plot points coloured by temperature with errorbars moved behind (can't colour errorbars apparently)
#scatter_kwargs = {"zorder":100}
#error_kwargs = {"lw":.5, "zorder":0}
#pl.scatter(density,pressure,c=temperature,**scatter_kwargs)
#pl.errorbar(density,pressure,yerr=stdpressure,fmt=None, marker=None, mew=0,**error_kwargs )
#pl.show()
