import os
import numpy as np
import matplotlib.pyplot as plt
from MDPlotData import MD_PlotData
from HeaderData import *

# Objects
Header = HeaderData(open('../../results/simulation_header','r'))
Getter = MD_PlotData('../../results/')

# Parameters
nbins = int(Header.gnbins1)*int(Header.gnbins2)*int(Header.gnbins3)
filebytes = os.path.getsize('../../results/pVA_k')
maxrec = filebytes / (9*8*nbins) # 8 is the size of 1 double precision 

# PLOTS
fig = plt.figure(figsize=(18.,15.))

drec = 200
maxrec = maxrec - maxrec%drec
cnt = 0

subrows = 4
subcols = 3

for rec in range(drec/2,maxrec,drec):

	minr = rec - drec/2
	maxr = rec + drec/2
	print(str(rec)+' of '+str(maxrec))

	def plot_pressure_tensor():
		
		y, P_k = Getter.get_pVA_prof_args('pVA_k',1,minr,maxr,peculiar=True)
		#y, P_kt = Getter.get_pVA_prof_args('pVA_k',1,minr,maxr,peculiar=False)
		y, P_c = Getter.get_pVA_prof_args('pVA_c',1,minr,maxr)
		P = P_k + P_c

		names = ['$P_{xx}$','$P_{xy}$','$P_{xz}$',
		         '$P_{yx}$','$P_{yy}$','$P_{yz}$',
		         '$P_{zx}$','$P_{zy}$','$P_{zz}$']

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
				limits = (0.0,2.5)
				color = 'c' 
			elif (component in [1,3]): 
				limits = (-0.04,0.01)
				color = 'm'
			else:
				limits = (-0.1,0.1)
				color = 'y'

			ax = fig.add_subplot(subrows,subcols,subcnt)
			ax.locator_params(tight=True, nbins=4)	
			ax.set_xlim(limits)
			ax.set_xlabel(name)
			ax.set_ylabel('$y$')
			ax.plot(P_k[:,component],y,'o',mfc='none',mec=color,label='k')
			#ax.plot(P_kt[:,component],y,'1',mfc='none',mec=color,label='kt')
			ax.plot(P_c[:,component],y,'x',color=color,label='c')
			ax.plot(  P[:,component],y,'o',color=color,label='t')
			
			if (component == 1):
				ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
				          ncol=4, mode="expand", borderaxespad=0.)

	def plot_velocity_fields():
		
		y, x, v = Getter.get_vplane_splot_args(2,0,1,minr,maxr)
		vx = v[:,:,0] 
		ax = fig.add_subplot(subrows,subcols,10)
		ax.set_xlabel('$v_x(x,y)$')
		ax.set_ylabel('$y$')
		contour = ax.contourf(y,x,vx,40)
		plt.colorbar(contour)

		y, v = Getter.get_vslice_plot_args(1,minr,maxr)
		ax = fig.add_subplot(subrows,subcols,11)
		ax.set_ylabel('$y$')
		ax.set_xlabel('$v$')
		#ax.set_xlim((-0.5,2.5))
		#ax.set_ylim((min(z),max(z)))
		ax.plot(v[:,0],y,'ro',label='$v_x$')
		ax.plot(v[:,1],y,'go',label='$v_y$')
		ax.plot(v[:,2],y,'bo',label='$v_z$')
		ax.legend(loc=0)	

		#z, v = Getter.get_vslice_plot_args(2,minr,maxr)
		#ax = fig.add_subplot(subrows,subcols,12)
		#ax.set_ylabel('$z$')
		#ax.set_xlabel('$v$')
		##ax.set_xlim((-0.5,7.5))
		##ax.set_ylim((min(r),max(r)))
		#ax.plot(v[:,0],z,'r',label='$v_x$')
		#ax.plot(v[:,1],z,'g',label='$v_y$')
		#ax.plot(v[:,2],z,'b',label='$v_z$')
		#ax.legend(loc=0)	

	def plot_other():

		#theta, z, P = Getter.get_pVA_splot_args(0,1,2,minr,maxr)
		#Prt = P[:,:,1] 
		#ax = fig.add_subplot(subrows,subcols,12)
		#ax.set_xlabel('$P_{r\\theta}(\\theta,z)$')
		#ax.set_ylabel('$z$')
		#contour = ax.contourf(theta, z, Prt , levels=np.arange(-0.1,0.4,0.01))
		#plt.colorbar(contour)

		y, P = Getter.get_pVA_prof_args(1,minr,maxr)
		ax = fig.add_subplot(subrows,subcols,12)
		ax.set_ylabel('$y$')
		ax.set_xlabel('$P$')
		ax.set_xlim((-0.1,0.3))
		ax.set_ylim((min(z),max(z)))
		ax.plot(P[:,1],y,'m',label='$P_{rt}$')
		ax.plot(P[:,2],y,'y',label='$P_{rz}$')
		ax.legend(loc=0)	

	def save_and_clear():
		plt.savefig('multiplot.'+"%05d"%cnt+'.png')
		plt.clf()

	plot_pressure_tensor()
	plot_velocity_fields()
	#plot_other()
	save_and_clear()
	cnt += 1
