import os
import numpy as np
import matplotlib.pyplot as plt
from MDPlotData import MD_PlotData
from HeaderData import *

# Objects
Header = HeaderData(open('../../results/simulation_header','r'))
Getter = MD_PlotData('../../results/',cpol_bins=True)

# Parameters
nbins = int(Header.gnbins1)*int(Header.gnbins2)*int(Header.gnbins3)
r_ii = float(Header.r_ii)
r_oi = float(Header.r_oi)
r_io = float(Header.r_io)
r_oo = float(Header.r_oo)
filebytes = os.path.getsize('../../results/pVA_k')
maxrec = filebytes / (9*8*nbins) # 8 is the size of 1 double precision 
print('Maxrec = '+str(maxrec))

# PLOTS
fig = plt.figure(figsize=(18.,15.))

drec = 40
maxrec = maxrec - maxrec%drec
cnt = 0

subrows = 4
subcols = 3

for rec in range(drec/2,maxrec,drec):

	minr = rec - drec/2
	maxr = rec + drec/2
	print(str(rec)+' of '+str(maxrec))

	def plot_pressure_tensor():
		
		r, P_k = Getter.get_pVA_prof_args('pVA_k',0,minr,maxr,peculiar=True)
		r, P_kt = Getter.get_pVA_prof_args('pVA_k',0,minr,maxr,peculiar=False)
		r, P_c = Getter.get_pVA_prof_args('pVA_c',0,minr,maxr)
		P = P_k + P_c

		names = ['$P_{rr}$'       ,'$P_{r\\theta}$'      ,'$P_{rz}$'       ,
		         '$P_{\\theta r}$','$P_{\\theta\\theta}$','$P_{\\theta z}$',
		         '$P_{zr}$'       ,'$P_{z\\theta}$'      ,'$P_{zz}$'       ]

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
				limits = (0.0,6.0)
				color = 'c' 
			elif (component in [1,3]): 
				limits = (0.0,0.2)
				color = 'm'
			else:
				limits = (-0.2,0.2)
				color = 'y'

			ax = fig.add_subplot(subrows,subcols,subcnt)
			ax.locator_params(tight=True, nbins=4)	
			ax.set_xlim(limits)
			ax.set_xlabel(name)
			ax.set_ylabel('$r$')
			ax.plot(P_k[:,component],r,'o',mfc='none',mec=color)
			ax.plot(P_kt[:,component],r,'1',mfc='none',mec=color)
			ax.plot(P_c[:,component],r,'x',color=color)
			ax.plot(  P[:,component],r,'o',color=color)

	def plot_velocity_fields():
		
		theta, z, v = Getter.get_vplane_splot_args(0,1,2,minr,maxr)
		vtheta = v[:,:,1] 
		ax = fig.add_subplot(subrows,subcols,10)
		ax.set_xlabel('$v_\\theta(\\theta,z)$')
		ax.set_ylabel('$z$')
		contour = ax.contourf(theta,z,vtheta,40)
		plt.colorbar(contour)

		z, v = Getter.get_vslice_plot_args(2,minr,maxr)
		ax = fig.add_subplot(subrows,subcols,11)
		ax.set_ylabel('$z$')
		ax.set_xlabel('$v$')
		ax.set_xlim((-0.5,2.5))
		ax.set_ylim((min(z),max(z)))
		ax.plot(v[:,0],z,'r',label='$v_r$')
		ax.plot(v[:,1],z,'g',label='$v_\\theta$')
		ax.plot(v[:,2],z,'b',label='$v_z$')
		ax.legend(loc=0)	

		r, v = Getter.get_vslice_plot_args(0,minr,maxr)
		ax = fig.add_subplot(subrows,subcols,12)
		ax.set_ylabel('$r$')
		ax.set_xlabel('$v$')
		ax.set_xlim((-0.5,7.5))
		ax.set_ylim((min(r),max(r)))
		ax.plot(v[:,0],r,'ro',label='$v_r$')
		ax.plot(v[:,1],r,'go',label='$v_\\theta$')
		ax.plot(v[:,2],r,'bo',label='$v_z$')
		ax.legend(loc=0)	

	def plot_other():

		#theta, z, P = Getter.get_pVA_splot_args(0,1,2,minr,maxr)
		#Prt = P[:,:,1] 
		#ax = fig.add_subplot(subrows,subcols,12)
		#ax.set_xlabel('$P_{r\\theta}(\\theta,z)$')
		#ax.set_ylabel('$z$')
		#contour = ax.contourf(theta, z, Prt , levels=np.arange(-0.1,0.4,0.01))
		#plt.colorbar(contour)

		z, P = Getter.get_pVA_prof_args(2,minr,maxr)
		ax = fig.add_subplot(subrows,subcols,12)
		ax.set_ylabel('$z$')
		ax.set_xlabel('$P$')
		ax.set_xlim((-0.1,0.3))
		ax.set_ylim((min(z),max(z)))
		ax.plot(P[:,1],z,'m',label='$P_{rt}$')
		ax.plot(P[:,2],z,'y',label='$P_{rz}$')
		ax.legend(loc=0)	

	def save_and_clear():
		plt.savefig('multiplot.'+"%05d"%cnt+'.png')
		plt.clf()

	plot_pressure_tensor()
	plot_velocity_fields()
	#plot_other()
	save_and_clear()
	cnt += 1
