import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle

sys.path.append('../../../')
import postproclib as ppl


#def plotm(x,y,**kwargs):

#    


#Get MD profile data
#fdir = '/home/es205/results/heat_flux/initial_study/'
fdir = '/home/es205/results/heat_flux/run_2/'
#fdir = '/home/es205/results/heat_flux/initial_VA_study/iter10000_20000/'
#fdir = '../../../../MD_dCSE/src_code/results/'
#startrec = 30; endrec = 995
startrec = 1; endrec = 15
nrecs = (endrec-startrec)+1

#Get y MD profile
PPObj = ppl.MD_PostProc(fdir)
varname = PPObj.plotlist.keys()[0].replace(' ','_')
vars()[varname] = PPObj.plotlist[PPObj.plotlist.keys()[0]]
print(PPObj,varname)
temp = eval(varname + '.profile(axis=1,startrec=' + str(startrec) 
                                     + ',endrec=' + str(endrec) + ')')
y_MD = temp[0]

vardict = {}; objdict ={}
for plotObj in PPObj.plotlist.keys():
    #Get Field Object Name 
    varname = plotObj.replace(' ','_')
    print(varname)
    vars()[varname] = PPObj.plotlist[plotObj]
    objdict[varname] = PPObj.plotlist[plotObj]
    
    #Read Field Object Profile Data
    temp = eval(varname + '.profile(axis=1,startrec=' + str(startrec) 
                                         + ',endrec=' + str(endrec) + ')')
    vars()[varname + '_p'] = temp[1] 
    vardict[varname + '_p'] = temp[1]


# Save a dictionary into a pickle file.
pickle.dump( objdict, open( "heatflux_obj.p", "wb" ) )
pickle.dump( vardict, open( "heatflux.p", "wb" ) )

grid = mbins.grid
bsx = np.mean(np.diff(grid[0]))
bsy = np.mean(np.diff(grid[1]))
bsz = np.mean(np.diff(grid[2]))
V = bsx*bsy*bsz
Ax = bsy*bsz
Ay = bsx*bsz
Az = bsx*bsy
Lx = float(mbins.Raw.header.globaldomain1)
Ly = float(mbins.Raw.header.globaldomain2)
Lz = float(mbins.Raw.header.globaldomain3)
hbsy = bsy/2.
y_MD_surf = y_MD #+ hbsy

#Analytical
utop = 1.; ubot = -1.; U = utop - ubot
Ttop = 1.05  ; Tbot = 1.05   
#dTdy_top =  ; dTdy_bot = 
fluiddensity = 0.8 #mbins.Raw.header.liquiddensity
walldensity = mbins.Raw.header.density
visc = 1.8; condct = 0.5
dt = float(mbins.Raw.header.delta_t)

Lyliquid = Ly - float(mbins.Raw.header.tethdistbot2) - float(mbins.Raw.header.tethdisttop2)
wallbinsbot = int( np.ceil(float(mbins.Raw.header.tethdistbot2)/float(mbins.Raw.header.binsize2)))
wallbinstop = int( np.ceil(float(mbins.Raw.header.tethdisttop2)/float(mbins.Raw.header.binsize2)))
liquidbins = int(Lyliquid/float(mbins.Raw.header.binsize2))
liquidstart = wallbinsbot
liquidend   = wallbinsbot+liquidbins
if (int(mbins.Raw.header.nbins2) != wallbinsbot+wallbinstop+liquidbins):
    print('Error = ', int(mbins.Raw.header.nbins2), wallbinsbot, wallbinstop, liquidbins)

#Analytical solutions
u_analy = np.empty(y_MD.shape)
u_analy[0:wallbinsbot] = ubot
u_analy[liquidstart:liquidend] = np.linspace(ubot,utop,liquidbins)
u_analy[liquidend:] = utop

tau_analy = np.zeros(y_MD.shape)
tau_analy[wallbinsbot:(wallbinsbot+liquidbins)] = visc * U/Lyliquid

lin_heatflux_analy = np.zeros(y_MD.shape)
#lin_heatflux_analy[wallbinsbot:(wallbinsbot+liquidbins)] = condct * np.linspace(dTdy_bot,dTdy_top,liquidbins)

rows = 4
columns = 2 
fig, axs = plt.subplots(rows,columns)
fig.set_figwidth(15*1.9411764705882353)
fig.set_figheight(15)

titles=[]
for row in range(rows):
	titles.append([])

styles = {'VAc':'r-', 'CVc':'rx','VAk':'b-','CVk':'bs','VA':'g-','CV':'go'}
ms = 5; lt = 1

#Mass/energy plot
titles[0].append('Mass Energy')
axs[0,0].plot(y_MD,rho_p[:,0],'-s',label=r'$\rho$', markersize=ms)
#Divide bins by no. of flux records 625 means snaps and bins have the same value. Doesn't mean this value is correct/meaningful!
#    axs[0,0].plot(y_MD,esnap_p[:,0],'-^',label='$e_{snap}$', markersize=ms)
#    axs[0,0].plot(y_MD,ebins_p[:,0]/625.,'-^',label='$e_{bins}$', markersize=ms) 
axs[0,0].plot(y_MD,Energy_p[:,0],'-^',label='$Energy$', markersize=ms)
axs[0,0].plot(y_MD,rhoEnergy_p[:,0]/(nrecs*fluiddensity),'o-',label=r'$\rho Energy/\rho$', markersize=ms)

#Velocity/momentum plot
titles[1].append('Velocity Momentum')
axs[1,0].plot(y_MD,u_analy[:],'-k',label=r'$u_{analy} $', lw=4)

axs[1,0].plot(y_MD,rho_u_p[:,0]/(nrecs*fluiddensity),'-r',label=r'$\rho u^{VA}/\rho$', markersize=ms)
axs[1,0].plot(y_MD,u_p[:,0],'--b',label=r'$u^{VA}$', markersize=ms)

axs[1,0].plot(y_MD,rho_uCV_p[:,0]/(fluiddensity),'or',label=r'$ \rho u^{CV}/\rho$', markersize=ms)
axs[1,0].plot(y_MD,uCV_p[:,0],'^b',label=r'$u_{CV}$', markersize=ms)

axs[1,0].set_ylim((-2.0,2.0))

#Shear Pressure plot
titles[2].append('Shear Pressure')

#MOP/CV method
ixyz = 1; bs = np.mean(np.diff(y_MD))
PiMOP_p = psurface_p + vflux_p
axs[2,0].plot(y_MD,-tau_analy[:],'-k',label=r'$\tau_{analy} $', lw=4)
axs[2,0].plot(y_MD_surf,psurface_p[:,ixyz],styles['CVc'],label=r'$ \sigma_{yx}^{MOP} $', markersize=ms)
axs[2,0].plot(y_MD_surf,vflux_p[:,ixyz],styles['CVk'],label=r'$ \kappa_{yx}^{MOP} $', markersize=ms)
axs[2,0].plot(y_MD_surf,PiMOP_p[:,ixyz],styles['CV'],label=r'$ \Pi_{yx}^{MOP} $', markersize=ms)
axs[2,0].plot(y_MD_surf,rho_uuCV_p[:,ixyz],'c^',label=r'$ \rho u v^{MOP} $', markersize=ms)

#VA method
ixyz = 1
PiVA_p = pVA_c_p + pVA_k_p
axs[2,0].plot(y_MD,pVA_c_p[:,ixyz],styles['VAc'],label=r'$ \sigma_{xy}^{\mathcal{VA}} $', lw=lt)
axs[2,0].plot(y_MD,pVA_k_p[:,ixyz],styles['VAk'],label=r'$ \kappa_{xy}^{\mathcal{VA}} $', lw=lt)
axs[2,0].plot(y_MD,PiVA_p[:,ixyz],styles['VA'],label=r'$ \Pi_{xy}^{\mathcal{VA}} $', lw=lt)
axs[2,0].plot(y_MD,rhouu_p[:,ixyz],'c-',label=r'$ \rho u v $', lw=lt)

#Plot Temperature
titles[0].append('Temperature')
axs[0,1].plot(y_MD,T_p[:,0],'-^',label='$T$', markersize=ms)
#axs[0,1].plot(y_MD,dTdr_p[:,1]/nrecs,'-^',label=r'$\frac{dT}{dy}$', markersize=ms)
axs[0,1].set_ylim((0.8,1.3))

#Plot energy qx component
titles[1].append('Energy Equation Terms in x')
axs[1,1].plot(y_MD,esurface_p[:,0],styles['CVc'],label='$fijvi^{MOP}$', markersize=ms)
axs[1,1].plot(y_MD,eflux_p[:,0],styles['CVk'],label='$evi^{MOP}$', markersize=ms)
axs[1,1].plot(y_MD,esurface_p[:,0]+eflux_p[:,0],styles['CV'],label='$[fijvidS+ evi]^{MOP}$', markersize=ms)
axs[1,1].plot(y_MD,CV_stressheat_p[:,0],'c^',label='$Stress Heating MOP$', markersize=ms)
axs[1,1].plot(y_MD,CVrhouE_p[:,0],'md',label=r'$\rho u E ^{MOP}$', markersize=ms)


axs[1,1].plot(y_MD,hfVA_c_p[:,0],styles['VAc'],label='$fijvi^{VA}$', markersize=ms)
axs[1,1].plot(y_MD,hfVA_k_p[:,0],styles['VAk'],label='$evi^{VA}$', markersize=ms)
axs[1,1].plot(y_MD,hfVA_c_p[:,0]+hfVA_k_p[:,0],styles['VA'],label='$[fijvidS+ evi]^{VA}$', markersize=ms)
axs[1,1].plot(y_MD,pVA_stressheat_p[:,0],'c-',label='$Stress Heating VA$', markersize=ms)
axs[1,1].plot(y_MD,rhouE_p[:,0],'m-',label=r'$\rho u E^{VA}$', markersize=ms)


#Plot energy qy component
titles[2].append('Energy Equation Terms in y')
axs[2,1].plot(y_MD,esurface_p[:,1],styles['CVc'],label='$fijvi^{MOP}$', markersize=ms)
axs[2,1].plot(y_MD,eflux_p[:,1],styles['CVk'],label='$evi^{MOP}$', markersize=ms)
axs[2,1].plot(y_MD,esurface_p[:,1]+eflux_p[:,1],styles['CV'],label='$[fijvidS+ evi] ^{MOP}$', markersize=ms)
axs[2,1].plot(y_MD,CV_stressheat_p[:,1],'c^',label='$Stress Heating^{MOP}$', markersize=ms)
axs[2,1].plot(y_MD,CVrhouE_p[:,1],'md',label=r'$\rho u E^{MOP}$', markersize=ms)

axs[2,1].plot(y_MD,hfVA_c_p[:,1],styles['VAc'],label='$fijvi^{VA}$', markersize=ms)
axs[2,1].plot(y_MD,hfVA_k_p[:,1],styles['VAk'],label='$evi^{VA}$', markersize=ms)
axs[2,1].plot(y_MD,hfVA_c_p[:,1]+hfVA_k_p[:,1],styles['VA'],label='$[fijvidS+ evi]^{VA}$', markersize=ms)
axs[2,1].plot(y_MD,pVA_stressheat_p[:,1],'c-',label='$Stress Heating^{VA}$', markersize=ms)
axs[2,1].plot(y_MD,rhouE_p[:,1],'m-',label=r'$\rho u E^{VA}$', markersize=ms)


#Direct Pressure plot
titles[3].append('Direct Pressure')

#MOP/CV method
PiMOP_p = psurface_p + vflux_p
axs[3,0].plot(y_MD,PiMOP_p[:,0],styles['VAc'], label=r'$ \Pi_{xx} ^{MOP} $', markersize=ms)
axs[3,0].plot(y_MD,PiMOP_p[:,4],styles['VAk'],label=r'$ \Pi_{yy} ^{MOP} $', markersize=ms)
axs[3,0].plot(y_MD,PiMOP_p[:,8],styles['VA'], label=r'$ \Pi_{zz} ^{MOP} $', markersize=ms)

#for ixyz in [0,4,8]:
    #axs[3,0].plot(y_MD_surf,psurface_p[:,ixyz],'-o',label=r'$ \sigma_{xx} dS_x^- $', markersize=ms)
    #axs[3,0].plot(y_MD_surf,vflux_p[:,ixyz],'-b',label=r'$ \kappa_{xx} dS_x^- $', markersize=ms)
    #axs[3,0].plot(y_MD_surf,PiMOP_p[:,ixyz],'-x',label=r'$ \Pi_{xx} dS_x^- $', markersize=ms)


#VA method
PiVA_p = pVA_c_p + pVA_k_p
axs[3,0].plot(y_MD,PiVA_p[:,0], 'ro',label=r'$ \Pi_{xx}^{\mathcal{VA}} $', markersize=ms)
axs[3,0].plot(y_MD,PiVA_p[:,4], 'bx',label=r'$ \Pi_{yy}^{\mathcal{VA}} $', markersize=ms)
axs[3,0].plot(y_MD,PiVA_p[:,8], 'g^',label=r'$ \Pi_{zz}^{\mathcal{VA}} $', markersize=ms)

#for ixyz in [0,4,8]:
    #axs[3,0].plot(y_MD,pVA_c_p[:,ixyz],'-x',label=r'$ \sigma_{xx}^{\mathcal{VA}} $', markersize=ms,alpha=0.8)
    #axs[3,0].plot(y_MD,pVA_k_p[:,ixyz],'-r',label=r'$ \kappa_{xx}^{\mathcal{VA}} $', markersize=ms,alpha=0.8)
    #axs[3,0].plot(y_MD,PiVA_p[:,ixyz], 'o',label=r'$ \Pi_{xx}^{\mathcal{VA}} $', markersize=ms,alpha=0.4)
    #axs[3,0].plot(y_MD,rhouu_p[:,ixyz],'-y',label=r'$ \rho u u $', markersize=ms,alpha=0.8)


#Plot energy q components
#a,q_approx_p = q_approx.profile(1,startrec,endrec,k=-5.,c=10.,f=-0.5)
dy = np.mean(np.diff(T.grid[1]))
dTdy = np.gradient(T_p[:,0],dy)
dudy = np.gradient(u_p[:,0],dy)
k=-5.; c=10.; f=-0.5
q_approx_p = np.empty((dudy.shape[0],3))
q_approx_p[:,0] = (k + 3. * f * np.power(dudy[:],2) ) * dTdy[:] 
q_approx_p[:,1] = c * dudy[:] * dTdy[:]
q_approx_p[:,2] = 0.

titles[3].append(r'Heat Flux')
axs[3,1].plot(y_MD,q_p[:,0],'-r',label='$q_x^{VA}$', markersize=ms)
axs[3,1].plot(y_MD,q_p[:,1],'--b',label='$q_y^{VA}$', markersize=ms)
axs[3,1].plot(y_MD,q_CV_p[:,0],'or',label='$q_x ^{MOP}$', markersize=ms)
axs[3,1].plot(y_MD,q_CV_p[:,1],'sb',label='$q_y ^{MOP}$', markersize=ms)
axs[3,1].plot(y_MD,q_approx_p[:,0],'k-',label='$q_x approx$', lw=3)
axs[3,1].plot(y_MD,q_approx_p[:,1],'k--',label='$q_y approx$', lw=3)


axs[3,1].set_ylim((-0.3,0.3))
#axs[3,1].plot(y_MD_surf,esurface_p[:,1],'-rs',label='$fijvidS_y^+$', markersize=ms)
#axs[3,1].plot(y_MD_surf,eflux_p[:,1],'-bx',label='$evidS_y^+$', markersize=ms)
#axs[3,1].plot(y_MD_surf,esurface_p[:,1]+eflux_p[:,1],'-g^',label='$[fijvidS+ evi] dS_y^+$', markersize=ms)
#axs[3,1].plot(y_MD_surf,stress_heating_p[:,1],'-k.',label='$Stress Heating MOP$', markersize=ms)
#axs[3,1].plot(y_MD_surf,pVAheat_p[:,1],'--y',label='$Stress Heating VA$', markersize=ms)



fsize = 8
plt.rc('font', size=fsize)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

axs[0,0].legend(loc='best')
axs[0,1].legend(loc='best')
axs[1,0].legend(loc='best',ncol=2)
axs[1,1].legend(loc='best',ncol=5)
axs[2,0].legend(loc='best',ncol=9)
axs[2,1].legend(loc='best',ncol=5)
axs[3,0].legend(loc='best',ncol=3)
axs[3,1].legend(loc='best',ncol=3)

# Set Titles
for j in range(columns):
    for i in range(rows):
	    axs[i,j].set_title(titles[i][j])


# Save the full figure...
#plt.show()
figname = 'Heat_flux.pdf'
plt.savefig(figname)
os.system('evince ' + figname)


# Save just the portion _inside_ the second axis's boundaries
#for j in range(columns):
#	for i in range(rows):
#		extent = axs[i,j].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#		# Pad the saved area by 10% in the x-direction and 20% in the y-direction
#		fig.savefig(titles[i][j].replace(' ','_') + '.png', bbox_inches=extent.expanded(1.1, 1.2))






#Top surfaces for verification
#ixyz = ixyz+9
#axs[2,0].plot(y_MD_surf,psurface_p[:,ixyz],'-gs',label=r'$ \sigma_{yx} dS_y^+ $',alpha=0.8, markersize=ms)
#axs[2,0].plot(y_MD_surf,vflux_p[:,ixyz],'-gx',label=r'$ \kappa_{yx} dS_y^+ $',alpha=0.8, markersize=ms)
#axs[2,0].plot(y_MD_surf,psurface_p[:,ixyz]+vflux_p[:,ixyz],'-go',label=r'$ \Pi_{yx} dS_y^+ $',alpha=0.8, markersize=ms)

#axs[1,1].plot(y_MD_surf,esurface_p[:,3],'-gs',label='$fijvidS_x^-$',alpha=0.8, markersize=ms)
#axs[1,1].plot(y_MD_surf,eflux_p[:,3],'-gx',label='$evidS_x^-$',alpha=0.8, markersize=ms)
#axs[1,1].plot(y_MD_surf,esurface_p[:,3]+eflux_p[:,3],'-go',label='$evidS_x^-$',alpha=0.8, markersize=ms)

#axs[2,1].plot(y_MD_surf,esurface_p[:,4],'-gs',label='$fijvidS_y^-$',alpha=0.8, markersize=ms)
#axs[2,1].plot(y_MD_surf,eflux_p[:,4],'-gx',label='$evidS_x^-$',alpha=0.8, markersize=ms)
#axs[2,1].plot(y_MD_surf,esurface_p[:,4]-eflux_p[:,4],'-go',label='$evidS_x^-$',alpha=0.8, markersize=ms)






