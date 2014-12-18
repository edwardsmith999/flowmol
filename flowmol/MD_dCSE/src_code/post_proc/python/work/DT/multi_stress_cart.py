#! /usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../../')
from MDPlotData import MD_PlotData
from HeaderData import *

# Objects
Header = HeaderData(open('../../../../results/simulation_header','r'))
Getter = MD_PlotData('../../../../results/',cpol_bins=False)

# Parameters
nbins = int(Header.gnbins1)*int(Header.gnbins2)*int(Header.gnbins3)
filebytes = os.path.getsize('../../../../results/pVA')
maxrec = filebytes / (9*8*nbins) # 8 is the size of 1 double precision 
print('Maxrec = '+str(maxrec))

# PLOTS
fig = plt.figure(figsize=(18.,15.))

drec = 2
maxrec = maxrec - maxrec%drec
cnt = 0

subrows = 3
subcols = 3

for rec in range(drec/2,maxrec,drec):

	minr = rec - drec/2
	maxr = rec + drec/2
	print(str(rec)+' of '+str(maxrec))

	def plot_pressure_tensor():
		
	#	y, P_k = Getter.get_pVA_prof_args('pVA_k',1,minr,maxr,peculiar=True)
	#	y, P_c = Getter.get_pVA_prof_args('pVA_c',1,minr,maxr)
	#	P = P_k + P_c
		y, P = Getter.get_pVA_prof_args('pVA',1,minr,maxr,peculiar=False)

		names = ['$P_{xx}$','$P_{xy}$','$P_{xz}$'       ,
		         '$P_{yx}$','$P_{yy}$','$P_{yz}$',
		         '$P_{zx}$','$P_{zy}$','$P_{zz}$'       ]

		comps = [  0  ,  1  ,  2  ,
		           3  ,  4  ,  5  , 
		           6  ,  7  ,  8  ]

		defs  = zip(names,comps)

	
		for nc in defs:

			name = nc[0]
			component = nc[1]
			# Add one for position of subplot in figure
			subcnt = component + 1
			if (component in [0,4,8]): 
				limits = (0.0,8.0)
				color = 'c' 
			elif (component in [1,3]): 
				limits = (-0.2,0.)
				color = 'm'
			else:
				limits = (-0.2,0.2)
				color = 'y'

			ax = fig.add_subplot(subrows,subcols,subcnt)
			ax.locator_params(tight=True, nbins=4)	
			ax.set_xlim(limits)
			ax.set_xlabel(name)
			ax.set_ylabel('$y$')
#			ax.plot(P_k[:,component],y,'o',mfc='none',mec=color,label='k')
#			ax.plot(P_c[:,component],y,'x',color=color,label='c')
			ax.plot(  P[:,component],y,'o',color=color,label='t')

#		y, v = Getter.get_vslice_plot_args(1,minr,maxr)
#		dy = y[1] - y[0] 
#		dvdy = 0.0*v[:,1]
#		for s in np.arange(1,len(dvdy)):
#			dvdy[s] = (v[s,0] - v[s-1,0])/dy
#		dvdy[0] = dvdy[1]
#		visco = -P[:,1] / dvdy
#		ax = fig.add_subplot(subrows,subcols,12)
#		ax.set_ylabel('$y$')
#		ax.set_ylim((min(y),max(y)))
#		ax.plot(visco,y,'yo',label='$\\mu$')
#		ax.legend(loc=0)

	def plot_other():
		
		y, v = Getter.get_vslice_plot_args(1,minr,maxr)
		ax = fig.add_subplot(subrows,subcols,10)
		ax.set_ylabel('$y$')
		ax.set_xlabel('$v$')
		ax.set_ylim((min(y),max(y)))
		ax.plot(v[:,0],y,'ro',label='$v_x$')
		ax.plot(v[:,1],y,'go',label='$v_y$')
		ax.plot(v[:,2],y,'bo',label='$v_z$')
		ax.legend(loc=0)	

		y, density = Getter.get_density_prof_args(1,minr,maxr)
		ax = fig.add_subplot(subrows,subcols,11)
		ax.set_ylabel('$r$')
		ax.set_xlabel('$\\rho$')
		ax.set_xlim((0.0,1.0))
		ax.set_ylim((min(y),max(y)))
		ax.plot(density,y,'ko')

	def save_and_clear():
		plt.savefig('multiplot.'+"%05d"%cnt+'.png')
		plt.clf()

	plot_pressure_tensor()
#	plot_other()
	save_and_clear()
	cnt += 1
