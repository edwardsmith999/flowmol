import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle
from scipy.optimize import curve_fit

sys.path.append('../../../')
import postproclib as ppl

# Load a dictionary into a pickle file.
varsDict = pickle.load( open( "heatflux.p", "rb" ) )
varsObjs = pickle.load( open( "heatflux_obj.p", "rb" ) )
nrecs = varsDict['nrecs']
print(varsDict.keys())
print(varsObjs.keys())
for plotObj in varsDict.keys():
    varname = plotObj.replace(' ','_')
    print(varname)
    vars()[varname] = varsDict[varname]

for plotObj in varsObjs.keys():
    varname = plotObj.replace(' ','_')
    print(varname)
    vars()[varname] = varsObjs[varname]

y_MD = mbins.grid[1]

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
utop = 1.4; ubot = -1.4; U = utop - ubot
Ttop = 1.05  ; Tbot = 1.05   
#dTdy_top =  ; dTdy_bot = 
fluiddensity = 0.01273239545 #mbins.Raw.header.liquiddensity
walldensity = mbins.Raw.header.density
visc = 0.14; condct = 0.5
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

rows = 2
columns = 2 
fig, axs = plt.subplots(rows,columns)
fig.set_figwidth(15*1.9411764705882353)
fig.set_figheight(15)

titles=[]
for row in range(rows):
	titles.append([])

styles = {'VAc':'r-', 'CVc':'rx','VAk':'b-','CVk':'bs','VA':'g-','CV':'go'}
ms = 5; lt = 1

def linear_fn(x, b, c):
    return  b * x + c

def quadratic_fn(x, a, b, c):
    return a*np.power(x,2.) + b * x + c

def qx_approx_fn(x, c):
    return c * dTdy[x] * dudy[x]

def qy_approx_fn(x, k, f):
    return (k + 3. * f * np.power(dudy[x],2) ) * dTdy[x] 

#Plot energy q components
liquidstart = liquidstart +2
liquidend = liquidend -3
#a,q_approx_p = q_approx.profile(1,startrec,endrec,k=-5.,c=10.,f=-0.5)
dy = np.mean(np.diff(T.grid[1]))
dTdy = np.gradient(T_p[:,0],dy)
dudy = np.gradient(u_p[:,0],dy)

k=-5.; f=3.; c=5.
liquidbins = range(liquidstart,liquidend)
q_approx_p = np.empty((dudy[liquidbins].shape[0],3))
q_approx_p[:,0] = qx_approx_fn(liquidbins,c)  #(k + 3. * f * np.power(dudy[:],2) ) * dTdy[:] 
q_approx_p[:,1] = qy_approx_fn(liquidbins,k,f) #c * dudy[:] * dTdy[:]
q_approx_p[:,2] = 0.

popTt, pcovT = curve_fit(quadratic_fn, y_MD[liquidbins], T_p[liquidbins,0], (1., 1., 1.))
print('Equation for T(y) = ', popTt[0], r'$y^2 + $ ', popTt[1],  r'$y +$ ',popTt[2],'with errors ' , pcovT)
T_approx_ls = quadratic_fn(y_MD[liquidbins], *popTt)
axs[0,0].plot(y_MD[liquidbins],T_p[liquidbins],label=r'$T$')
axs[0,0].plot(y_MD[liquidbins],T_approx_ls,label=r'$T^{LS}$')
axs[0,0].legend()
dTLSdy = 2.* popTt[0] * y_MD + popTt[1]

efl = 5    # exclude_for_linfit
popxt, pcovx = curve_fit(qx_approx_fn, liquidbins, q_p[liquidbins,0], (1.))
qx_approx_ls = qx_approx_fn(liquidbins, *popxt)
poptlinqx, pcovlinqx = curve_fit(linear_fn, y_MD[liquidbins[efl:-efl]], q_p[liquidbins[efl:-efl],0], (1., 0.))
qx_linfit_ls = linear_fn(y_MD[liquidbins[efl:-efl]], *poptlinqx)

popty, pcovy = curve_fit(qy_approx_fn, liquidbins, q_p[liquidbins,1], (1., 1.))
qy_approx_ls = qy_approx_fn(liquidbins, *popty)
poptlinqy, pcovlinqy = curve_fit(linear_fn, y_MD[liquidbins], q_p[liquidbins,1], (1., 1.))
qy_linfit_ls = linear_fn(y_MD[liquidbins], *poptlinqy)



print(popxt,popty)


print(qy_approx_ls.shape, y_MD[liquidbins].shape,qx_approx_ls.shape)

axs[0,1].plot(y_MD[liquidbins],dudy[liquidbins],'r',label=r'$dudy$')
axs[0,1].plot(y_MD[liquidbins],np.power(dudy[liquidbins],2.),'r--',label=r'$(dudy)^2$')
axs[0,1].plot(y_MD[liquidbins],dTdy[liquidbins],'b',label=r'$dTdy$')
axs[0,1].plot(y_MD[liquidbins],dTLSdy[liquidbins],'k',label=r'$dTdy_{linfit}$')
axs[0,1].plot(y_MD[liquidbins],dTdy[liquidbins]*dudy[liquidbins],'g',label=r'$dTdy*dudy$')
axs[0,1].legend()

#axs[1,0].plot(y_MD[liquidbins],qx_approx_ls,'b-',label='$q_x^{ls}$', markersize=ms)
axs[1,0].plot(y_MD[liquidbins[efl:-efl]],q_approx_p[efl:-efl,0],'r-',label=r'$NLC \;  q_x$')
axs[1,0].plot(y_MD[liquidbins[efl:-efl]],q_p[liquidbins[efl:-efl],0],'b-',label='$q_x^{VA}$')
axs[1,0].plot(y_MD[liquidbins[efl:-efl]],qx_linfit_ls,'k-',label='$q_x^{linfit}$')
axs[1,0].legend()

#axs[1,1].plot(y_MD[liquidbins],qy_approx_ls,'b-',label='$q_y^{ls}$', markersize=ms)
axs[1,1].plot(y_MD[liquidbins],q_approx_p[:,1],'r-',label=r'$NLC \; q_y$')
axs[1,1].plot(y_MD[liquidbins],q_p[liquidbins,1],'b-',label='$q_y^{VA}$')
axs[1,1].plot(y_MD[liquidbins],qy_linfit_ls,'k-',label='$q_y^{linfit}$')
axs[1,1].legend()
plt.show()


fsize = 8
plt.rc('font', size=fsize)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

axs[0].legend(loc='best')


# Save the full figure...
plt.show()
#figname = 'Heat_flux.pdf'
#plt.savefig(figname)
#os.system('evince ' + figname)


# Save just the portion _inside_ the second axis's boundaries
#for j in range(columns):
#	for i in range(rows):
#		extent = axs[0][i,j].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#		# Pad the saved area by 10% in the x-direction and 20% in the y-direction
#		fig.savefig(titles[i][j].replace(' ','_') + '.png', bbox_inches=extent.expanded(1.1, 1.2))






#Top surfaces for verification
#ixyz = ixyz+9
#axs[0][2,0].plot(y_MD_surf,psurface_p[:,ixyz],'-gs',label=r'$ \sigma_{yx} dS_y^+ $',alpha=0.8, markersize=ms)
#axs[0][2,0].plot(y_MD_surf,vflux_p[:,ixyz],'-gx',label=r'$ \kappa_{yx} dS_y^+ $',alpha=0.8, markersize=ms)
#axs[0][2,0].plot(y_MD_surf,psurface_p[:,ixyz]+vflux_p[:,ixyz],'-go',label=r'$ \Pi_{yx} dS_y^+ $',alpha=0.8, markersize=ms)

#axs[0][1,1].plot(y_MD_surf,esurface_p[:,3],'-gs',label='$fijvidS_x^-$',alpha=0.8, markersize=ms)
#axs[0][1,1].plot(y_MD_surf,eflux_p[:,3],'-gx',label='$evidS_x^-$',alpha=0.8, markersize=ms)
#axs[0][1,1].plot(y_MD_surf,esurface_p[:,3]+eflux_p[:,3],'-go',label='$evidS_x^-$',alpha=0.8, markersize=ms)

#axs[0][2,1].plot(y_MD_surf,esurface_p[:,4],'-gs',label='$fijvidS_y^-$',alpha=0.8, markersize=ms)
#axs[0][2,1].plot(y_MD_surf,eflux_p[:,4],'-gx',label='$evidS_x^-$',alpha=0.8, markersize=ms)
#axs[0][2,1].plot(y_MD_surf,esurface_p[:,4]-eflux_p[:,4],'-go',label='$evidS_x^-$',alpha=0.8, markersize=ms)






