import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import cPickle as pickle
from scipy.optimize import curve_fit

sys.path.append('../../../')
import postproclib as ppl

# Load a dictionary into a pickle file.
nrecs = 1596
varsDict = pickle.load( open( "heatflux.p", "rb" ) )
varsObjs = pickle.load( open( "heatflux_obj.p", "rb" ) )
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

rows = 2
columns = 4
fig, axs = plt.subplots(rows,columns)
fig.set_figwidth(15*1.9411764705882353)
fig.set_figheight(15)

titles=[]
for row in range(rows):
	titles.append([])

styles = {'VAc':'r-', 'CVc':'kx','VAk':'b-','CVk':'ks','VA':'g-','CV':'ko'}
ms = 2; lt = 1

def linear_fn(x, b, c):
    return  b * x + c

def quadratic_fn(x, a, b, c):
    return a*np.power(x,2.) + b * x + c

titles[1].append('Energy Equation Terms in x')



liquidbins = range(liquidstart,liquidend)

ixyz = 0
axs[ixyz,0].plot(y_MD[liquidbins],pVA_stressheat_p[liquidbins,ixyz],styles['VAc'],label='$\Pi \cdot u$', markersize=ms)
axs[ixyz,0].plot(y_MD[liquidbins],hfVA_c_p[liquidbins,ixyz],'g',label='$fijvi^{VA}$', markersize=ms)
axs[ixyz,0].plot(y_MD[liquidbins],esurface_p[liquidbins,ixyz]-CV_stressheat_p[liquidbins,ixyz],styles['CVc'],label='$fijvi^{MOP}-\Pi \cdot u$', markersize=ms)

axs[ixyz,1].plot(y_MD[liquidbins],hfVA_k_p[liquidbins,ixyz],styles['VAk'],label=r'$evi^{VA}$', markersize=ms)
axs[ixyz,1].plot(y_MD[liquidbins],rhouE_p[liquidbins,ixyz],'g',label=r'$\rho uE$', markersize=ms)
axs[ixyz,1].plot(y_MD[liquidbins],eflux_p[liquidbins,ixyz]-CVrhouE_p[liquidbins,ixyz],styles['CVk'],label='$evi^{MOP} - \rho u E$', markersize=ms)


axs[ixyz,2].plot(y_MD[liquidbins],esurface_p[liquidbins,ixyz]+eflux_p[liquidbins,ixyz],styles['CV'],label='$[fijvidS+ evi]^{MOP}$', markersize=ms)
axs[ixyz,2].plot(y_MD[liquidbins],hfVA_c_p[liquidbins,ixyz]+hfVA_k_p[liquidbins,ixyz],styles['VA'],label='$[fijvidS+ evi]^{VA}$', markersize=ms)


q_VA = ( hfVA_c_p[liquidbins,ixyz]-pVA_stressheat_p[liquidbins,ixyz]
        +hfVA_k_p[liquidbins,ixyz]-rhouE_p[liquidbins,ixyz]         )
q_CV= ( esurface_p[liquidbins,ixyz]-CV_stressheat_p[liquidbins,ixyz]
        +  eflux_p[liquidbins,ixyz]-CVrhouE_p[liquidbins,ixyz]      )
axs[ixyz,3].plot(y_MD[liquidbins],q_CV,styles['CV'],label='$q^{MOP}$', markersize=ms)
axs[ixyz,3].plot(y_MD[liquidbins],q_VA,styles['VA'],label='$q^{VA}$', markersize=ms)

ixyz = 1
axs[ixyz,0].plot(y_MD[liquidbins],pVA_stressheat_p[liquidbins,ixyz],styles['VAc'],label='$\Pi \cdot u$', markersize=ms)
axs[ixyz,0].plot(y_MD[liquidbins],hfVA_c_p[liquidbins,ixyz],'g',label='$fijvi^{VA}$', markersize=ms)
axs[ixyz,0].plot(y_MD[liquidbins],esurface_p[liquidbins,ixyz]-CV_stressheat_p[liquidbins,ixyz],styles['CVc'],label='$fijvi^{MOP}-\Pi \cdot u$', markersize=ms)


axs[ixyz,1].plot(y_MD[liquidbins],hfVA_k_p[liquidbins,ixyz],'g',label='$evi^{VA}$', markersize=ms)
axs[ixyz,1].plot(y_MD[liquidbins],rhouE_p[liquidbins,ixyz],styles['VAk'],label='$rhouE_p$', markersize=ms)
axs[ixyz,1].plot(y_MD[liquidbins],eflux_p[liquidbins,ixyz]-CVrhouE_p[liquidbins,ixyz],styles['CVk'],label=r'$evi^{MOP}- \rho u E$', markersize=ms)


axs[ixyz,2].plot(y_MD[liquidbins],esurface_p[liquidbins,ixyz]+eflux_p[liquidbins,ixyz],styles['CV'],label='$[fijvidS+ evi]^{MOP}$', markersize=ms)
axs[ixyz,2].plot(y_MD[liquidbins],hfVA_c_p[liquidbins,ixyz]+hfVA_k_p[liquidbins,ixyz],styles['VA'],label='$[fijvidS+ evi]^{VA}$', markersize=ms)

q_VA = ( hfVA_c_p[liquidbins,ixyz]-pVA_stressheat_p[liquidbins,ixyz]
        +hfVA_k_p[liquidbins,ixyz]-rhouE_p[liquidbins,ixyz]         )
q_CV= ( esurface_p[liquidbins,ixyz]-CV_stressheat_p[liquidbins,ixyz]
        +  eflux_p[liquidbins,ixyz]-CVrhouE_p[liquidbins,ixyz]      )
axs[ixyz,3].plot(y_MD[liquidbins],q_CV,styles['CV'],label='$q^{MOP}$', markersize=ms)
axs[ixyz,3].plot(y_MD[liquidbins],q_VA,styles['VA'],label='$q^{VA}$', markersize=ms)


for ixyz in range(0,2):
    axs[ixyz,0].legend(loc='best')
    axs[ixyz,1].legend(loc='best')
    axs[ixyz,2].legend(loc='best')
    axs[ixyz,3].legend(loc='best')

plt.savefig('./energy_terms.pdf')
#plt.show()


