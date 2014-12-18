import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('../../../')
import postproclib as ppl

def flip(a):
    return -a[::-1]

#Get MD profile data
fdir = '/media/My Passport/Work/MD_turbulence/Highres_iter_2800000_16ybinpercell/results/'

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
MDlayer_end = 70
MDlayer_mid = int(MDlayer_strt + 0.5*(MDlayer_end - MDlayer_strt))
#Start of viscous sublayer
sublayer_strt = MDlayer_end
sublayer_end = 500
sublayer_mid = int(sublayer_strt + 0.5*(sublayer_end - sublayer_strt))
#Buffer region
bufferlayer_strt = sublayer_end 
bufferlayer_end = 1500
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
dudy_top =  np.max(dudy[wall_stress_loc:-wall_stress_loc])
dudy_bot = -np.min(dudy[wall_stress_loc:-wall_stress_loc])
dudy_ave = 0.5*(dudy_top + dudy_bot)
u_tau = np.power(dudy_ave*visc,0.5)
delta_tau = np.power(visc/dudy_ave,0.5)

#Get plus units
y_plus_MD  = (y_MD[fluid_strt:mid]-y_MD[fluid_strt])/delta_tau;
y_plus_CFD = (y_MD[sublayer_strt:mid]-y_MD[sublayer_strt])/delta_tau;
u_plus_MD = mom_MD/(fluiddensity*u_tau);
u_plus_MD_loc = u_MD/u_tau

#Setup analytical solutions
#kappa = 0.41; alpha = 5.1
kappa = 1.0; alpha = -3.3
#y_plus_wall = (y_MD[sublayer_strt:mid])/delta_tau
strt = sublayer_strt-MDlayer_strt
#Shift velocity to zero at start and add current MD velocity

visc_sublayer = y_plus_MD + u_plus_MD_loc[fluid_strt,0]  #-y_plus_MD[strt] + u_plus_MD[sublayer_strt,0]
loglaw = (1./kappa) * np.log(y_plus_MD) + alpha #+ u_plus_MD_loc[fluid_strt,0] 

print(visc_sublayer,loglaw)

fsize = 16
fig = plt.figure()
#ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(111)

print(y_plus_MD.shape,rho_MD.shape,fluid_strt,mid)

#ax1.plot(y_plus_MD,rho_MD[fluid_strt:mid],'r-')
#ax1.plot(y_plus_MD,flip(rho_MD[mid:-fluid_strt]),'b-')
ax2.plot(y_plus_MD,u_plus_MD[fluid_strt:mid,0],'r-')
ax2.plot(y_plus_MD,flip(u_plus_MD[mid:-fluid_strt,0]),'b-')

ax2.plot(y_plus_MD,u_plus_MD_loc[fluid_strt:mid,0],'r--')
ax2.plot(y_plus_MD,flip(u_plus_MD_loc[mid:-fluid_strt,0]),'b--')


plt.rc('font', size=fsize)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

ax2.text(y_plus_MD[1],-8,r'$\underbrace{\phantom{a + b + c + d + e + f}}_{Subviscous}$',fontsize=fsize)
ax2.text(y_plus_MD[sublayer_strt-MDlayer_strt],-8,r'$\underbrace{\phantom{a + b + c + d + e + f   \; }}_{Viscous}$',fontsize=fsize)
ax2.text(y_plus_MD[bufferlayer_strt-MDlayer_strt],-8,r'$\underbrace{\phantom{a + b + \; \; }}_{Buffer}$',fontsize=fsize)
ax2.text(y_plus_MD[loglaw_strt-MDlayer_strt],-8,r'$\underbrace{\phantom{a + b \; }}_{Loglaw}$',fontsize=fsize)

ax2.plot(y_plus_MD,visc_sublayer,'k-')
ax2.plot(y_plus_MD,loglaw,'k--')
#ax1.set_xscale('log')
ax2.set_xscale('log')
#ax1.set_ylim([0.2,0.55])
ax2.set_ylim([-9,1.01])

plt.show()

