import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.optimize import curve_fit
import scipy
from matplotlib import rc

sys.path.append('../../../')
import postproclib as ppl

def flip(a):
    return -a[::-1]

#Get MD profile data
fdir = '/home/es205/scratch/Re400/Highres_iter_2800000_16ybinpercell/results'

densityObj = ppl.MD_dField(fdir)
momentumObj = ppl.MD_momField(fdir)
velocityObj = ppl.MD_vField(fdir)
startrec = 12; endrec = 14; nrecs = (endrec-startrec)+1
y_MD,rho_MD = densityObj.profile(axis=1,startrec=startrec,endrec=endrec)
y_MD,mom_MD = momentumObj.profile(axis=1,startrec=startrec,endrec=endrec)
y_MD,u_MD   = velocityObj.profile(axis=1,startrec=startrec,endrec=endrec)

#Plot limits
mid =  len(mom_MD)/2.

#Divide by number of records
mom_MD = mom_MD/nrecs

#Take off laminar solution
#mom_MD[:,0] = mom_MD[:,0] - np.linspace(-1.0,1.0,mom_MD.shape[0])

#Take of one from every cell
#mom_MD[:,0] = np.ones(mom_MD[:,0].shape)

#fig, axs = plt.subplots(2)
#fluiddensity = 0.3
#for ax in axs:
#    ax.plot(rho_MD[:,0],'-s')
#    ax.plot(mom_MD[:,0]/fluiddensity,'-o')
#    ax.plot(mom_MD[:,0]/rho_MD[:,0],'-^')
#    ax.plot(u_MD[:,0],'-x')
#axs[1].set_xscale('log')
#plt.show()


#Start of fluid part of channel
fluid_strt = 38
#Start of MD stacking sublayer
MDlayer_strt = fluid_strt
MDlayer_end = 85
MDlayer_mid = int(MDlayer_strt + 0.5*(MDlayer_end - MDlayer_strt))
#Start of viscous sublayer
sublayer_strt = MDlayer_end
sublayer_end = 400
sublayer_mid = int(sublayer_strt + 0.5*(sublayer_end - sublayer_strt))
#Buffer region
bufferlayer_strt = sublayer_end 
bufferlayer_end = 1400
bufferlayer_mid = int(bufferlayer_strt + 0.5*(bufferlayer_end - bufferlayer_strt))
#Loglaw start
loglaw_strt = bufferlayer_end
loglaw_end = mid
loglaw_mid = int(loglaw_strt + 0.5*(loglaw_end - loglaw_strt))
#Location to use for wall stress
wall_stress_loc = 60


#Scaling parameters
visc = 1.8; fluiddensity = 0.3
dudy = np.diff(u_MD[:,0])/np.diff(y_MD)

#def func(x, a, b, c, d, e, f, g, h, i):
#  return g*np.cos(h*x) + i*x**6. + f*x**5. + e*x**4. + d*x**3. + c*x**2. + b*x + a

def func(x, A):
    return A*np.log(np.tan(x*0.25*np.pi/(0.5*x.max())))

#p0 = scipy.array([1])
#coeffs, matcov = curve_fit(func, y_MD, u_MD[:,0], p0)
#yaj = func(y_MD, *coeffs)
yaj = func(y_MD, 2.5)
plt.plot(y_MD, u_MD[:,0],'x',y_MD,yaj,'r-')
plt.show()




dudy_top =  np.max(dudy[wall_stress_loc:-wall_stress_loc])
dudy_bot = -np.min(dudy[wall_stress_loc:-wall_stress_loc])
dudy_ave = 0.5*(dudy_top + dudy_bot)
u_tau = np.power(dudy_ave*visc,0.5)
delta_tau = np.power(visc/dudy_ave,0.5)

#Get plus units
y_plus_wall = y_MD[0:fluid_strt]/delta_tau;
y_plus_MD  = (y_MD[fluid_strt:mid]-y_MD[fluid_strt])/delta_tau;
y_plus_CFD = (y_MD[sublayer_strt:mid]-y_MD[sublayer_strt])/delta_tau;
u_plus_MD = mom_MD/(fluiddensity*u_tau);
u_plus_MD_loc = u_MD/u_tau

#Setup analytical solutions
#kappa = 0.41; alpha = 5.5
kappa = 1.0; alpha = -3.15
#y_plus_wall = (y_MD[sublayer_strt:mid])/delta_tau
strt = sublayer_strt-MDlayer_strt
#Shift velocity to zero at start and add current MD velocity

visc_sublayer = y_plus_MD + u_plus_MD_loc[fluid_strt,0]  #-y_plus_MD[strt] + u_plus_MD[sublayer_strt,0]
loglaw = (1./kappa) * np.log(y_plus_MD) + alpha #+ u_plus_MD_loc[fluid_strt,0] 

fsize = 16
fig = plt.figure()
#ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(111)

#ax2.plot(y_plus_wall,0.5*(u_plus_MD[0:fluid_strt,0]+flip(u_plus_MD[-fluid_strt:,0])),'b-',alpha=0.8)
ax2.plot(y_plus_MD,0.5*(u_plus_MD[fluid_strt:mid,0]+flip(u_plus_MD[mid:-fluid_strt,0])),'r-',alpha=0.8)
ax2.plot(y_plus_MD,0.5*(u_plus_MD_loc[fluid_strt:mid,0]+flip(u_plus_MD_loc[mid:-fluid_strt,0])),'b--')



plt.rc('font', size=fsize)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

labley = -14
ax2.text(y_plus_MD[1],labley,r'$\underbrace{\phantom{a + b + c + d + e + f + g + h }}_{\mathcal{MD} \; Stick-Slip}$',fontsize=fsize)
ax2.text(y_plus_MD[sublayer_strt-MDlayer_strt],labley,r'$\underbrace{\phantom{a + b + c + d \;\;\; }}_{Viscous}$',fontsize=fsize)
ax2.text(y_plus_MD[bufferlayer_strt-MDlayer_strt],labley,r'$\underbrace{\phantom{a + b + c  \; }}_{Buffer}$',fontsize=fsize)
ax2.text(y_plus_MD[loglaw_strt-MDlayer_strt],labley,r'$\underbrace{\phantom{a + b \; }}_{Loglaw}$',fontsize=fsize)

ax2.plot(y_plus_MD,visc_sublayer,'k-')
ax2.plot(y_plus_MD[bufferlayer_strt:],loglaw[bufferlayer_strt:],'k--')
#ax1.set_xscale('log')
ax2.set_xscale('log')
ax2.set_xlim([6e-3,31])
xticks=ax2.get_xticks().tolist()

ax2.set_ylim([-15,1.01])

ax2.set_xlabel('$y^+$')
ax2.set_ylabel('$u^+$')

f.subplots_adjust(wspace=0)
plt.rc('font', size=24)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

#plt.show()
plt.savefig('./law_of_the_wall.pdf')
plt.savefig('/home/es205/Documents/Turbulent_Couette/Resubmission_butwhere/law_of_the_wall.pdf')

