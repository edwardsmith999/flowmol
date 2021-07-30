showplots = list(range(14))+[1.5,2.5,2.75,2.9,9.75]

import matplotlib.pyplot as plt
import numpy as np
import sys

#sys.call("git clone https://github.com/montefra/mimic_alpha.git")
import mimic_alpha as ma

sys.path.insert(0, '/home/es205/codes/flowmol_et_al/utils/misclib/')
#from misc import pastel

ppdir = '/home/es205/codes/python/pyDataView'
sys.path.append(ppdir)
import postproclib as ppl

import pickle
try:
    with open("summary.p", "rb") as f:
        data = pickle.load(f, encoding="latin1")
except FileNotFoundError:
    raise FileNotFoundError("No summary.p, run Flowmol simulation and" + 
                            " then use Analyse_data_and_Pickle.py")

def interp(a, indx, pad=True):
    l = [a[i+1,indx]+a[i,indx] for i in range(a.shape[0]-1)]
    if pad:
        l.append(a[-1,0])
    return 0.5*np.array(l)

def integrate(integrand, mn, mx, dx, inttype="trapz"):
    if inttype is "simp":
        from scipy.integrate import simps
        return np.array([0] + [simps(integrand[mn:n], dx=dx, even="first") 
                       for n in range(mn+1,mx)])
    elif inttype is "trapz":
        return np.array([np.trapz(integrand[mn:n], dx=dx) 
                       for n in range(mn,mx)])
    else:
        raise IOError("inttype should be simp or trapz")

def get_integrand(P):

    PN = P[:,0]
    PT =  0.5*(P[:,4]+P[:,8])
    integrand = PN - PT

    return integrand

fsize = 14
plt.rc('font', size=fsize)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

#Load the colormap -- full range of colors can be got with
# colourscheme(0) to colourscheme(256)
# Good values for lines:
#  Dark blue = 0
#  Light blue = 64
#  Yellow = 150
#  Light Red = 175
#  Red = 256
#colourscheme=plt.cm.RdYlBu_r
#dashes = [(None,None),[5,5],[1,3],[5,3,1,3],[5,3,1,3,1,3]]
#colors = [colourscheme(i) for i in range(0,257,32)]
#c1 = colors[2]
#c2 = colors[-1]
#c3 = colors[-2]
#c4 = colors[1]

colours = ['DarkBlue','SkyBlue','LightGreen','Orange','Red', "Black"]
colours = [plt.cm.tab10(i) for i in range(6)]
c1, c2, c3, c4, c5, c6 = [c for c in colours]



Obj = ppl.MD_dummyField("./results")
dx = float(Obj.header.binsize1)
Lx = float(Obj.header.globaldomain1)

#Get profile
x = data["x"]
rho = data["rho"]
dsurf = data["dsurf_vflux"]
psurf =data["psurface"]
vflux = data["vflux"]
pVAc = data["pVA_c"]
pVAk = data["pVA_k"]

pIK1c = data["pIK1_c"]
pIK1k = data["pIK1_k"]

#Shift so surface at zero
xp = x - dx/2.


#Patch molecules in surface
xmn = -3.0
xmx = 2.0
sb = rho.argmax()
smn = np.array([x>xmn-2*dx]).argmax()
smx = np.array([x>xmx+2*dx]).argmax()
indx = range(smn,smx)
#indx.pop(sb)



#####################################
#
#  Plot Surface of tension without Normal  
#
#####################################
if 1 in showplots:
    fig, ax = plt.subplots(1,2, sharey=True)
    indx = range(0,psurf.shape[0])

    PCV = psurf[indx,:]+vflux[indx,:]
    #PCV[:,0] = interp(psurf,0)[indx] + interp(vflux,0)[indx] + dsurf[indx,0]
    #PCV[:,0] = psurf[indx,0]+vflux[indx,0]+dsurf[indx,0]
    PCV_c = psurf[indx,:]
    PCV_c[:,0] = interp(psurf,0)[indx]
    pVA = pVAc[indx,:]+pVAk[indx,:]
    pIK1= pIK1c[indx,:]+pIK1k[indx,:]

    intmin = -3
    intmax = 2.5
    xi = xp[indx]
    intrange = (xi>intmin) & (xi<intmax)
    imin = np.argmax(xi>intmin)
    imax = np.argmin(xi<intmax)

    integrandCV = 0.5*(PCV[:,4]+PCV[:,8])  #get_integrand(PCV)
    STCV = integrate(integrandCV*xi, imin, imax, dx)
    ax[1].plot(xi[intrange], STCV, c=c2)
    ax[0].plot(xi, integrandCV*xi, c=c2,label=r"$\int x \left[P_N-P_T\right]_{CV} dx$", zorder=9)

    integrandIK1 = 0.5*(pIK1[:,4]+pIK1[:,8]) #get_integrand(pIK1)
    STIK1 = integrate(integrandIK1*xi, imin, imax, dx)
    ax[1].plot(xi[intrange], STIK1, c=c1)
    ax[0].plot(xi, integrandIK1*xi, c=c1, label=r"$\int x \left[P_N-P_T\right]_{VA} dx$", zorder=9)

    integrandVA =  0.5*(pVA[:,4]+pVA[:,8]) # get_integrand(pVA)
    STVA = integrate(integrandVA*xi, imin, imax, dx)
    ax[1].plot(xi[intrange], STVA, c=c3)
    ax[0].plot(xi, integrandVA*xi, c=c3, label=r"$\int x \left[P_N-P_T\right]_{VA} dx$", zorder=9)

    for i in range(ax.shape[0]):
        ax[i].plot(xi, np.zeros(xi.shape[0]), "-k")

        ax[i].set_xlim([-3,2.])
        ax[i].set_ylim([-1.5,1.5])
        ax[i].set_xlabel("$z$")

    #ax[0].set_ylabel("PT*z and Surface of Tension")
    #plt.legend()

    #plt.savefig("PT_surface_of_tension.pdf", bbox_inches="tight")
    #plt.savefig("PT_surface_of_tension.eps", bbox_inches="tight")
    plt.show()

#########################################################################
#
#  Plot of AVERAGE IK1 only 
#
#########################################################################
if 1.5 in showplots:
    fig, ax = plt.subplots(1,2, sharey=True)
    a = 1.0

    #IK1 Pressures
    ax[0].plot(x[indx], pIK1c[indx,0],'--', c=c1, label=r"IK1 $\Pi^{c}_{N}$", alpha=a)
    ax[1].plot(x[indx], pIK1c[indx,4],'--', c=c1, label=r"IK1 $\Pi^{c}_{T}$", alpha=a) 
    ax[0].plot(x[indx], pIK1k[indx,0],'--', c=c2, label=r"IK1 $\Pi^{k}_{N}$", alpha=a) 
    ax[1].plot(x[indx], pIK1k[indx,4],'--', c=c2, label=r"IK1 $\Pi^{k}_{T}$", alpha=a) 
    ax[0].plot(x[indx], pIK1k[indx,0]+pIK1c[indx,0],'--', c=c3, label=r"IK1 $\Pi_{N}$") 
    ax[1].plot(x[indx], pIK1k[indx,4]+pIK1c[indx,4],'--', c=c3, label=r"IK1 $\Pi_{T}$") 


    for i in range(ax.shape[0]):
        ax[i].set_xlim([xmn,xmx])
        ax[i].set_xlabel("$z$")

    #ax[0].set_ylim([-1.5,1.2])
    ax[0].set_ylabel("Pressure")

    plt.savefig("IK1.pdf", bbox_inches="tight")
    #Needed as we lose transparency in eps otherwise
    import os
    os.system("pdftops -eps IK1.pdf IK1.eps")
    plt.show()


#########################################################################
#
#  Plot of AVERAGE VA, IK1 and CV stress FAINT COMPONENTS
#
#########################################################################



if 2 in showplots:
    fig, ax = plt.subplots(1,2, sharey=True)
    a = 0.05
    # c1bak = c1
    # c2bak = c2
    # c1 = [c1[i] for i in range(3)]+[a] #[ma.colorAlpha_to_rgb(c1, a)[0][i] for i in range(3)] + [1.]
    # c2 = [c2[i] for i in range(3)]+[a] #[ma.colorAlpha_to_rgb(c2, a)[0][i] for i in range(3)] + [1.]

    #CV pressures
    ax[0].plot(xp[indx], psurf[indx,0],'-', c=c1, label=r"CV $\Pi^{c}_{N}$", alpha=a)
    ax[1].plot(xp[indx], psurf[indx,4],'-', c=c1, label=r"{CV} $\Pi^{c}_{T}$", alpha=a) 
    ax[0].plot(xp[indx], vflux[indx,0],'-', c=c2, label=r"{CV} $\Pi^{k}_{N}$", alpha=a)
    ax[1].plot(xp[indx], vflux[indx,4],'-', c=c2, label=r"{CV} $\Pi^{k}_{T}$", alpha=a)
    #ax[0].plot(xp[indx], psurf[indx,0]+vflux[indx,0],'-', c=c3, label=r"{CV} $\Pi_{N}$") 
    ax[0].plot(xp[indx], psurf[indx,0]+vflux[indx,0]+dsurf[indx,0],'-', c=c3, label=r"{CV} $\Pi_{N}$") 

    ax[1].plot(xp[indx], psurf[indx,4]+vflux[indx,4],'-', c=c3, label=r"{CV} $\Pi_{T}$")

    #IK1 Pressures
    ax[0].plot(x[indx], pIK1c[indx,0],'--', c=c1, label=r"IK1 $\Pi^{c}_{N}$", alpha=a)
    ax[1].plot(x[indx], pIK1c[indx,4],'--', c=c1, label=r"IK1 $\Pi^{c}_{T}$", alpha=a) 
    ax[0].plot(x[indx], pIK1k[indx,0],'--', c=c2, label=r"IK1 $\Pi^{k}_{N}$", alpha=a) 
    ax[1].plot(x[indx], pIK1k[indx,4],'--', c=c2, label=r"IK1 $\Pi^{k}_{T}$", alpha=a) 
    ax[0].plot(x[indx], pIK1k[indx,0]+pIK1c[indx,0],'--', c=c3, label=r"IK1 $\Pi_{N}$") 
    ax[1].plot(x[indx], pIK1k[indx,4]+pIK1c[indx,4],'--', c=c3, label=r"IK1 $\Pi_{T}$") 
    #VA pressures
    ax[0].scatter(xp[indx], pVAc[indx,0], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5, alpha=a)
    ax[1].scatter(xp[indx], pVAc[indx,4], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{T}$",zorder=5, alpha=a)
    ax[0].scatter(xp[indx], pVAk[indx,0], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5, alpha=a)
    ax[1].scatter(xp[indx], pVAk[indx,4], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{T}$",zorder=5, alpha=a)
    ax[0].scatter(xp[indx], pVAk[indx,0]+pVAc[indx,0], c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{N}$",zorder=5)
    ax[1].scatter(xp[indx], pVAk[indx,4]+pVAc[indx,4], c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{T}$",zorder=5)

    for i in range(ax.shape[0]):
        ax[i].plot(xp[indx], np.zeros(xp[indx].shape[0]), "-k", lw=2.,zorder=0)
        ax[i].set_xlim([xmn,xmx])
        ax[i].set_xlabel("$z$")

    ax[0].set_ylim([-1.5,1.2])
    ax[0].set_ylabel("Pressure")

    plt.savefig("Intrinsic.pdf", bbox_inches="tight")
    #Needed as we lose transparency in eps otherwise
    import os
    os.system("pdftops -eps Intrinsic.pdf Intrinsic.eps")
    plt.show()


#########################################################################
#
#  Plot of VA, IK1 and CV stress divided by density
#
#########################################################################
if 2.5 in showplots:
    fig, ax = plt.subplots(1,2, sharey=True)
    a = 1

    #CV pressures
    ax[0].plot(xp[indx], interp(psurf,0)[indx]/rho[indx,0],'-', c=c1, label=r"CV $\Pi^{c}_{N}$", alpha=a)
    ax[1].plot(xp[indx], interp(psurf,4)[indx]/rho[indx,0],'-', c=c1, label=r"{CV} $\Pi^{c}_{T}$", alpha=a) 
    ax[0].plot(xp[indx], interp(vflux,0)[indx]/rho[indx,0],'-', c=c2, label=r"{CV} $\Pi^{k}_{N}$", alpha=a)
    ax[1].plot(xp[indx], interp(vflux,4)[indx]/rho[indx,0],'-', c=c2, label=r"{CV} $\Pi^{k}_{T}$", alpha=a)
    ax[0].plot(xp[indx], ( interp(psurf,0)[indx]
                          +interp(vflux,0)[indx]
                          +interp(dsurf,0)[indx])/rho[indx,0],
                          '-', c=c3, label=r"{CV} $\Pi_{N}$") 

    ax[1].plot(xp[indx], (interp(psurf,4)[indx]+interp(vflux,4)[indx])/rho[indx,0],'-', c=c3, label=r"{CV} $\Pi_{T}$")

    #VA pressures
    ax[0].scatter(xp[indx], pVAc[indx,0]/rho[indx,0], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5, alpha=a)
    ax[1].scatter(xp[indx], pVAc[indx,4]/rho[indx,0], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{T}$",zorder=5, alpha=a)
    ax[0].scatter(xp[indx], pVAk[indx,0]/rho[indx,0], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5, alpha=a)
    ax[1].scatter(xp[indx], pVAk[indx,4]/rho[indx,0], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{T}$",zorder=5, alpha=a)
    ax[0].scatter(xp[indx], (pVAk[indx,0]+pVAc[indx,0])/rho[indx,0], c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{N}$",zorder=5)
    ax[1].scatter(xp[indx], (pVAk[indx,4]+pVAc[indx,4])/rho[indx,0], c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{T}$",zorder=5)

    for i in range(ax.shape[0]):
        ax[i].set_xlim([xmn,xmx])
        ax[i].set_xlabel("$z$")

    #ax[0].set_ylim([-1.5,1.2])
    ax[0].set_ylabel("Pressure")

    plt.savefig("Intrinsic_density_normalised.pdf", bbox_inches="tight")
    #Needed as we lose transparency in eps otherwise
    import os
    os.system("pdftops -eps Intrinsic_density_normalised.pdf Intrinsic_density_normalised.eps")
    plt.show()


if 2.75 in showplots:

    fig, ax = plt.subplots(1,2)
    a = 1

    #CV pressures
    P_rho_CV = (interp(psurf,4)[indx]+interp(vflux,4)[indx])/rho[indx,0]
    ax[0].plot(xp[indx], P_rho_CV,'-', c=c3, label=r"{CV} $\Pi_{T}$") 
    ax[1].plot(xp[indx], P_rho_CV,'-', c=c3, label=r"{CV} $\Pi_{T}$")

    #VA pressures
    P_rho_VA = (pVAk[indx,4]+pVAc[indx,4])/rho[indx,0]
    ax[0].scatter(xp[indx], P_rho_VA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{T}$",zorder=5)
    ax[1].scatter(xp[indx], P_rho_VA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{T}$",zorder=5)

    ymn = -1.
    ymx = 0.5
    ax[0].fill_between([xmn,xmx],ymn*np.ones(2), ymx*np.ones(2), facecolor=[.5,.5,.5], alpha=0.3, zorder=10)
    #ax[0].hlines(xmn,xmx,ymn,'k',lw=2)
    #ax[0].hlines(xmn,xmx,ymx,'k',lw=2)

    for i in range(ax.shape[0]):
        ax[i].set_xlim([xmn,xmx])
        ax[i].set_xlabel("$z$")

    ax[1].set_ylim([ymn, ymx])
    ax[0].set_ylabel(r"$P/\rho$")

    plt.savefig("Intrinsic_tangent_density_normalised.pdf", bbox_inches="tight")
    #Needed as we lose transparency in eps otherwise
    import os
    os.system("pdftops -eps Intrinsic_tangent_density_normalised.pdf Intrinsic_tangent_density_normalised.eps")
    plt.show()


if 2.9 in showplots:

    from scipy.optimize import curve_fit
    def func(x, a,b):
        return np.exp(-x)*(1.+a*np.sin(2*x*np.pi + b))

    fig, ax = plt.subplots(1,1)
    a = 0.1

    #Attempt to fit decaying exponent
    xmn = -5.
    xmx = 3.0
    sb = rho.argmax()
    smn = np.array([x>xmn-2*dx]).argmax()
    smx = np.array([x>xmx+2*dx]).argmax()
    indx = range(smn,smx)


    ax.plot(xp[indx], interp(psurf,4)[indx]/rho[indx,0],'-', c=c1, alpha=a)
    ax.scatter(xp[indx], pVAc[indx,4]/rho[indx,0], c=c1, edgecolors='k', s=20.,zorder=5, alpha=a)

    ax.plot(xp[indx], interp(vflux,4)[indx]/rho[indx,0],'-', c=c2, alpha=a)
    ax.scatter(xp[indx], pVAk[indx,4]/rho[indx,0], c=c2, edgecolors='k', s=20. ,zorder=5, alpha=a)


    P_rho_CV = (interp(psurf,4)[indx]+interp(vflux,4)[indx])/rho[indx,0]
    P_rho_VA = (pVAk[indx,4]+pVAc[indx,4])/rho[indx,0]

    loc = np.array([P_rho_VA>-1.]).argmin()
    xdata = -xp[indx][loc:0:-1]
    ydata = -P_rho_VA[loc:0:-1]
    #ax.plot(xdata, ydata, 'o')


    ax.plot(xp[indx], P_rho_CV,'-', c=c3, label=r"{CV} $\Pi_{T}$")
    ax.scatter(xp[indx], P_rho_VA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{T}$",zorder=5)

    popt, pcov = curve_fit(func, xdata, ydata)
    print(popt)
    #Add plot   
    xlin = np.linspace(xdata.min(), xdata.max(), 100)
    #ax.plot(-xlin, -func(xlin, *popt), 'k--', lw=2)
    ax.plot(-xlin, -func(xlin, -2.5, 0.5), 'k-', lw=2)

    ax.set_xlim([-5, 3])

    ax.set_ylim([-1.5, 1])
    ax.set_ylabel(r"$P_T/\rho$")
    ax.set_xlabel(r"$z$")

    ax.plot(xp[indx], np.zeros(xp[indx].shape[0]), "-k", lw=2.,zorder=0)

    ax.text(-4.9, 0.3, r"$ e^{-z}\left[ 1 - 2.5\sin(2\pi z + 0.5) \right]$")

    plt.savefig("P_over_rho_fitted.pdf", bbox_inches="tight")
    #Needed as we lose transparency in eps otherwise
    import os
    os.system("pdftops -eps P_over_rho_fitted.pdf P_over_rho_fitted.eps")
    plt.show()


#####################################
#
#  Plot of VA, IK1 and CV stress NO AVERAGE
#
#####################################

if 3 in showplots:
    fig, ax = plt.subplots(1,2, sharey=True)

    ax[0].plot(xp[indx], rho[indx,0],'-', c=[0.3,0.3,0.3], label=r"$\rho$")
    ax[1].plot(xp[indx], rho[indx,0],'-', c=[0.3,0.3,0.3], label=r"$\rho$")

    #CV pressures
    ax[0].plot(xp[indx], psurf[indx,0],'-', c=c1, label=r"CV $\Pi^{c}_{N}$")
    ax[1].plot(xp[indx], psurf[indx,4],'-', c=c1, label=r"{CV} $\Pi^{c}_{T}$") 
    ax[0].plot(xp[indx], vflux[indx,0]+dsurf[indx,0],'-', c=c2, label=r"{CV} $\Pi_{N}$") 
    ax[1].plot(xp[indx], vflux[indx,4],'-', c=c2, label=r"{CV} $\Pi^{k}_{T}$")

    #IK1 Pressures
    ax[0].plot(x[indx], pIK1c[indx,0],'--', c=c1, label=r"IK1 $\Pi^{c}_{N}$")
    ax[1].plot(x[indx], pIK1c[indx,4],'--', c=c1, label=r"IK1 $\Pi^{c}_{T}$") 
    ax[0].plot(x[indx], pIK1k[indx,0],'--', c=c2, label=r"IK1 $\Pi^{k}_{N}$") 
    ax[1].plot(x[indx], pIK1k[indx,4],'--', c=c2, label=r"IK1 $\Pi^{k}_{T}$") 
    #VA pressures
    ax[0].scatter(xp[indx], pVAc[indx,0], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax[1].scatter(xp[indx], pVAc[indx,4], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{T}$",zorder=5)
    ax[0].scatter(xp[indx], pVAk[indx,0], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax[1].scatter(xp[indx], pVAk[indx,4], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{T}$",zorder=5)

    for i in range(ax.shape[0]):
        ax[i].plot(xp[indx], np.zeros(xp[indx].shape[0]), "-k", lw=2.,zorder=0)
        ax[i].set_xlim([xmn,xmx])
        ax[i].set_xlabel("$z$")

    ax[0].set_ylim([-1.5,1.2])
    ax[0].set_ylabel("Pressure")

    plt.savefig("Intrinsic_nosum.pdf", bbox_inches="tight")
    plt.savefig("Intrinsic_nosum.eps", bbox_inches="tight")
    plt.show()



#####################################
#
#  Plot of VA, IK1 and CV stress JUST AVERAGE
#
#####################################

if 4 in showplots:
    fig, ax = plt.subplots(1,2, sharey=True)

    #CV pressures
    ax[0].plot(xp[indx], psurf[indx,0]+vflux[indx,0]+dsurf[indx,0],'-', c=c3, label=r"{CV} $\Pi_{N}$") 
    ax[1].plot(xp[indx], psurf[indx,4]+vflux[indx,4],'-', c=c3, label=r"{CV} $\Pi_{T}$")

    #IK1 Pressures
    ax[0].plot(x[indx], pIK1k[indx,0]+pIK1c[indx,0],'--', c=c3, label=r"IK1 $\Pi_{N}$") 
    ax[1].plot(x[indx], pIK1k[indx,4]+pIK1c[indx,4],'--', c=c3, label=r"IK1 $\Pi_{T}$") 
    #VA pressures
    ax[0].scatter(xp[indx], pVAk[indx,0]+pVAc[indx,0], c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{N}$",zorder=5)
    ax[1].scatter(xp[indx], pVAk[indx,4]+pVAc[indx,4], c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{T}$",zorder=5)

    for i in range(ax.shape[0]):
        ax[i].plot(xp[indx], np.zeros(xp[indx].shape[0]), "-k", lw=2.,zorder=0)
        ax[i].set_xlim([xmn,xmx])
        ax[i].set_xlabel("$z$")

    ax[0].set_ylim([-1.5,1.2])
    ax[0].set_ylabel("Pressure")

    plt.savefig("Intrinsic_Ave.pdf", bbox_inches="tight")
    plt.savefig("Intrinsic_Ave.eps", bbox_inches="tight")
    plt.show()









#####################################
#
#  Plot of VA, IK1 and CV stress NO AVERAGE
#
#####################################

if 5 in showplots:

    fig, ax = plt.subplots(3,2, sharex=True, sharey='row', figsize=(8,10))
    fig.subplots_adjust(wspace=0, hspace=0)
    #plt.tight_layout()

    #CV pressures
    ax[0,0].plot(xp[indx], psurf[indx,0],'-', c=c1, label=r"CV $\Pi^{c}_{N}$")
    ax[0,1].plot(xp[indx], psurf[indx,4],'-', c=c1, label=r"{CV} $\Pi^{c}_{T}$") 
    ax[0,0].plot(xp[indx], vflux[indx,0]+dsurf[indx,0],'-', c=c2, label=r"{CV} $\Pi_{N}$") 
    ax[0,1].plot(xp[indx], vflux[indx,4],'-', c=c2, label=r"{CV} $\Pi^{k}_{T}$")

    #IK1 Pressures
    ax[0,0].plot(x[indx], pIK1c[indx,0],'--', c=c1, label=r"IK1 $\Pi^{c}_{N}$")
    ax[0,1].plot(x[indx], pIK1c[indx,4],'--', c=c1, label=r"IK1 $\Pi^{c}_{T}$") 
    ax[0,0].plot(x[indx], pIK1k[indx,0],'--', c=c2, label=r"IK1 $\Pi^{k}_{N}$") 
    ax[0,1].plot(x[indx], pIK1k[indx,4],'--', c=c2, label=r"IK1 $\Pi^{k}_{T}$") 
    #VA pressures
    ax[0,0].scatter(xp[indx], pVAc[indx,0], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax[0,1].scatter(xp[indx], pVAc[indx,4], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{T}$",zorder=5)
    ax[0,0].scatter(xp[indx], pVAk[indx,0], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax[0,1].scatter(xp[indx], pVAk[indx,4], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{T}$",zorder=5)

    for i in range(ax.shape[0]):
        ax[i,0].set_xlim([xmn,xmx])

    ax[0,0].set_ylim([-1.5,1.2])
    ax[0,0].set_ylabel("Pressure")

    #CV pressures
    ax[1,0].plot(xp[indx], psurf[indx,0]+vflux[indx,0]+dsurf[indx,0],'-', c=c3, label=r"{CV} $\Pi_{N}$") 
    ax[1,1].plot(xp[indx], psurf[indx,4]+vflux[indx,4],'-', c=c3, label=r"{CV} $\Pi_{T}$")

    #IK1 Pressures
    ax[1,0].plot(x[indx], pIK1k[indx,0]+pIK1c[indx,0],'--', c=c3, label=r"IK1 $\Pi_{N}$") 
    ax[1,1].plot(x[indx], pIK1k[indx,4]+pIK1c[indx,4],'--', c=c3, label=r"IK1 $\Pi_{T}$") 
    #VA pressures
    ax[1,0].scatter(xp[indx], pVAk[indx,0]+pVAc[indx,0], c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{N}$",zorder=5)
    ax[1,1].scatter(xp[indx], pVAk[indx,4]+pVAc[indx,4], c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi_{T}$",zorder=5)

    for i in range(ax.shape[0]):
        ax[i,1].set_xlim([xmn,xmx])

    ax[1,0].set_ylim([-1.5,1.2])
    ax[1,0].set_ylabel("Pressure")



    indx = range(0,psurf.shape[0])
    PCV = psurf[indx,:]+vflux[indx,:]
    PCV[:,0] = psurf[indx,0]+vflux[indx,0]+dsurf[indx,0]
    pVA = pVAc[indx,:]+pVAk[indx,:]
    pIK1= pIK1c[indx,:]+pIK1k[indx,:]

    intmin = -5
    intmax = 2
    xi = xp[indx]
    intrange = (xi>intmin) & (xi<intmax)
    imin = np.argmax(xi>intmin)
    imax = np.argmin(xi<intmax)


    integrandCV = PCV[:,0]
    STCV = integrate(integrandCV, imin, imax, dx)
    ax[2,0].plot(xi[intrange], STCV, c=c3)

    integrandIK1 = pIK1[:,0]
    STIK1 = integrate(integrandIK1, imin, imax, dx)
    ax[2,0].plot(xi[intrange], STIK1, '--', c=c3)

    integrandVA = pVA[:,0]
    STVA = integrate(integrandVA, imin, imax, dx)
    #ax[0].plot(xp[intrange], STVA, c=c3)
    ax[2,0].scatter(xp[intrange], STVA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)


    integrandCV = 0.5*(PCV[:,4]+PCV[:,8])
    STCV = integrate(integrandCV, imin, imax, dx)
    ax[2,1].plot(xi[intrange], -STCV, c=c3)

    integrandIK1 = 0.5*(pIK1[:,4]+pIK1[:,8])
    STIK1 = integrate(integrandIK1, imin, imax, dx)
    ax[2,1].plot(xi[intrange], -STIK1, '--', c=c3)

    integrandVA = 0.5*(pVA[:,4]+pVA[:,8])
    STVA = integrate(integrandVA, imin, imax, dx)
    #ax[1].plot(xp[intrange], -STVA, c=c3)
    ax[2,1].scatter(xp[intrange], -STVA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)

    for i in range(2):
        ax[2,i].plot(xi, np.zeros(xi.shape[0]), "-k")
        ax[2,i].set_xlim([-3,2.])
        ax[2,i].set_xlabel("$z$")

    #ax[0].set_ylim([-1.,1.])
    ax[2,0].set_ylim([-0.5,1.0])

    ax[2,0].set_ylabel("Surface Tension ")
    plt.tight_layout()
    plt.savefig("Intrinsic_all.pdf", bbox_inches="tight")
    plt.savefig("Intrinsic_all.eps", bbox_inches="tight")
    plt.show()





if 6 in showplots:

    fig, ax = plt.subplots(1,2, sharey=True)


    indx = range(0,psurf.shape[0])
    PCV = psurf[indx,:]+vflux[indx,:]
    PCV[:,0] = psurf[indx,0]+vflux[indx,0]+dsurf[indx,0]
    pVA = pVAc[indx,:]+pVAk[indx,:]
    pIK1= pIK1c[indx,:]+pIK1k[indx,:]

    intmin = -5
    intmax = 2
    xi = xp[indx]
    intrange = (xi>intmin) & (xi<intmax)
    imin = np.argmax(xi>intmin)
    imax = np.argmin(xi<intmax)

    integrandCV = PCV[:,0]
    STCV = integrate(integrandCV, imin, imax, dx)
    ax[0].plot(xi[intrange], STCV, c=c3)

    integrandIK1 = pIK1[:,0]
    STIK1 = integrate(integrandIK1, imin, imax, dx)
    ax[0].plot(xi[intrange], STIK1, '--', c=c3)

    integrandVA = pVA[:,0]
    STVA = integrate(integrandVA, imin, imax, dx)
    #ax[0].plot(xp[intrange], STVA, c=c3)
    ax[0].scatter(xp[intrange], STVA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)

    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
    axin = zoomed_inset_axes(ax[0], zoom=5, loc=1)
    axin.plot(xi[intrange], STCV, c=c3)
    axin.plot(xi[intrange], STIK1, '--', c=c3)
    axin.scatter(xp[intrange], STVA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    axin.set_xlim([1.35,1.93])
    axin.set_ylim([0.0,0.05])
    mark_inset(ax[0], axin, loc1=3, loc2=4, fc="none", ec="0.5")

    integrandCV = 0.5*(PCV[:,4]+PCV[:,8])
    STCV = integrate(integrandCV, imin, imax, dx)
    ax[1].plot(xi[intrange], -STCV, c=c3)

    integrandIK1 = 0.5*(pIK1[:,4]+pIK1[:,8])
    STIK1 = integrate(integrandIK1, imin, imax, dx)
    ax[1].plot(xi[intrange], -STIK1, '--', c=c3)

    integrandVA = 0.5*(pVA[:,4]+pVA[:,8])
    STVA = integrate(integrandVA, imin, imax, dx)
    #ax[1].plot(xp[intrange], -STVA, c=c3)
    ax[1].scatter(xp[intrange], -STVA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)


    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
    axin1 = zoomed_inset_axes(ax[1], zoom=5, loc=4,bbox_to_anchor=(0, .05, 1., 1.),
                   bbox_transform=ax[1].transAxes,)

    #pos1 = axin1.get_position() # get the original position 
    #pos2 = [pos1.x0, pos1.y0 + 20,  pos1.width, pos1.height] 
    #axin1.set_position(pos2) 

    axin1.plot(xi[intrange], -STCV, c=c3)
    axin1.plot(xi[intrange], -STIK1, '--', c=c3)
    axin1.scatter(xp[intrange], -STVA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    axin1.set_xlim([1.35,1.93])
    axin1.set_ylim([0.525,0.575])
    p = mark_inset(ax[1], axin1, loc1=1, loc2=2, fc="none", ec="0.5")
    
    for i in range(2):
        ax[i].plot(xi, np.zeros(xi.shape[0]), "-k", lw=2.,zorder=0)
        ax[i].set_xlim([-3,2.])
        ax[i].set_xlabel("$z$")
    axin.plot(xi, np.zeros(xi.shape[0]), "-k", lw=2.,zorder=0)

    #ax[0].set_ylim([-1.,1.])
    ax[0].set_ylim([-0.6,0.9])

    ax[0].set_ylabel(r"$\int_{-\infty}^{z} P(z') dz'$")



    plt.savefig("Intrinsic_ST.pdf", bbox_inches="tight")
    plt.savefig("Intrinsic_ST.eps", bbox_inches="tight")
    plt.show()



#fig, ax = plt.subplots(1,2, sharey=True)

##IK1 Pressures
#ax[0].plot(x[indx], pIK1c[indx,0],'--', c=c1, label=r"IK1 $\Pi^{c}_{N}$")
#ax[1].plot(x[indx], pIK1c[indx,4],'--', c=c1, label=r"IK1 $\Pi^{c}_{T}$") 
#ax[0].plot(x[indx], pIK1k[indx,0],'--', c=c2, label=r"IK1 $\Pi^{k}_{N}$") 
#ax[1].plot(x[indx], pIK1k[indx,4],'--', c=c2, label=r"IK1 $\Pi^{k}_{T}$") 
#ax[0].plot(x[indx], pIK1k[indx,0]+pIK1c[indx,0],'--', c=c3, label=r"IK1 $\Pi_{N}$") 
#ax[1].plot(x[indx], pIK1k[indx,4]+pIK1c[indx,4],'--', c=c3, label=r"IK1 $\Pi_{T}$") 
##VA pressures
#ax[0].plot(xp[indx], pVAc[indx,0], '-', c=c1, label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
#ax[1].plot(xp[indx], pVAc[indx,4], '-', c=c1, label=r"{VA} $\Pi^{c}_{T}$",zorder=5)
#ax[0].plot(xp[indx], pVAk[indx,0], '-', c=c2, label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
#ax[1].plot(xp[indx], pVAk[indx,4], '-', c=c2, label=r"{VA} $\Pi^{c}_{T}$",zorder=5)
#ax[0].plot(xp[indx], pVAk[indx,0]+pVAc[indx,0], '-', c=c3, label=r"{VA} $\Pi_{N}$",zorder=5)
#ax[1].plot(xp[indx], pVAk[indx,4]+pVAc[indx,4], '-', c=c3,  label=r"{VA} $\Pi_{T}$",zorder=5)

#for i in range(ax.shape[0]):
#    ax[i].set_xlim([xmn,xmx])
#    ax[i].set_xlabel("$z$")

#ax[0].set_ylim([-1.5,1.2])
#ax[0].set_ylabel("Pressure")

#plt.savefig("Intrinsic_noCV.pdf", bbox_inches="tight")
#plt.savefig("Intrinsic_noCV.eps", bbox_inches="tight")
#plt.show()









if 7 in showplots:

    fig, ax = plt.subplots(1,2, sharey=True)

    #CV pressures
    ax[0].plot(xp[indx], psurf[indx,0]+vflux[indx,0], c=c2, label=r"{CV} $\Pi_{N}$") 
    ax[1].plot(xp[indx], psurf[indx,4]+vflux[indx,4], c=c2, label=r"{CV} $\Pi_{T}$")

    #IK1 Pressures
    ax[0].plot(x[indx], pIK1k[indx,0]+pIK1c[indx,0], c=c1, label=r"IK1 $\Pi_{N}$") 
    ax[1].plot(x[indx], pIK1k[indx,4]+pIK1c[indx,4], c=c1, label=r"IK1 $\Pi_{T}$") 

    #VA pressures
    ax[0].plot(xp[indx], pVAk[indx,0]+pVAc[indx,0], c=c3, label=r"{VA} $\Pi_{N}$",zorder=5)
    ax[1].plot(xp[indx], pVAk[indx,4]+pVAc[indx,4], c=c3, label=r"{VA} $\Pi_{T}$",zorder=5)

    for i in range(ax.shape[0]):
        ax[i].set_xlim([xmn,xmx])
        ax[i].set_xlabel("$z$")

    ax[0].set_ylim([-1.5,1.2])
    ax[0].set_ylabel("Pressure")

    plt.savefig("Intrinsic_sums.pdf", bbox_inches="tight")
    plt.savefig("Intrinsic_sums.eps", bbox_inches="tight")
    plt.show()



#fig, ax = plt.subplots(1,1)

##CV pressures
#ax.plot(xp[indx], psurf[indx,0]+vflux[indx,0]-(psurf[indx,4]+vflux[indx,4]),'-', c=c1, label=r"{CV} $\Pi_{N}$") 

##IK1 Pressures
#ax.plot(x[indx], pIK1k[indx,0]+pIK1c[indx,0]-(pIK1k[indx,4]+pIK1c[indx,4]),'--', c=c2, label=r"IK1 $\Pi_{N}$") 

##VA pressures
#ax.scatter(xp[indx], pVAk[indx,0]+pVAc[indx,0]-(pVAk[indx,4]+pVAc[indx,4]), c=c4, edgecolors='k', s=20., label=r"{VA} $\Pi_{N}$",zorder=5)

#ax.set_xlim([xmn,xmx])
#ax.set_xlabel("$z$")

#ax.set_ylim([-1.5,1.2])
#ax.set_ylabel("Pressure")

#plt.savefig("Intrinsic_PNmPT.pdf", bbox_inches="tight")
#plt.savefig("Intrinsic_PNmPT.eps", bbox_inches="tight")
#plt.show()



#quit()

#####################################
#
#  Plot of VA vs CV stress
#
#####################################




##Plot everything
#fig, ax = plt.subplots(1,1)

##Plot line and hidden region for intrinsic interface
#for i in range(10):
#    ax.vlines(xp[sb],-2,2.,color='w', lw=2*i, zorder=5+i, alpha=1-0.1*i)
#ax.vlines(xp[sb],-2,2.,color='k', zorder=30)

##Density
#ax.plot(xp[indx], rho[indx,0],'-', c=[0.7,0.7,0.7], label=r"$\rho$")

##CV pressures
#ax.plot(xp[indx], psurf[indx,0],'-', c=c1, label=r"CV $\Pi^{c}_{N}$")
#ax.plot(xp[indx], psurf[indx,4],'-', c=c4, label=r"{CV} $\Pi^{c}_{T}$") 
#ax.plot(xp[indx], vflux[indx,4],'-', c=c5, label=r"{CV} $\Pi^{k}$")

##IK1 Pressures
#ax.plot(xp[indx], pIK1c[indx,0],'-', c=c2, label=r"IK1 $\Pi^{c}_{N}$")
#ax.plot(xp[indx], pIK1c[indx,4],'-', c=c3, label=r"IK1 $\Pi^{c}_{T}$") 
#ax.plot(xp[indx], pIK1k[indx,0],'-', c=c3, label=r"IK1 $\Pi^{k}$") 

##VA pressures
#ax.scatter(xp[indx], pVAc[indx,0], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
#ax.scatter(xp[indx], pVAc[indx,4], c=c4, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{T}$",zorder=5)
#ax.scatter(xp[indx], pVAk[indx,4], c=c5, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}$",zorder=5)

#ax.set_xlim([xmn,xmx])
#ax.set_ylim([-1.5,1.2])
#ax.set_xlabel("$z$")
#ax.set_ylabel("Pressure and Density")

#plt.savefig("VA_vs_CV.pdf", bbox_inches="tight")
#plt.savefig("VA_vs_CV.eps", bbox_inches="tight")
#plt.show()



#####################################
#
#  Plot CV balance
#
#####################################

if 8 in showplots:
    fig, ax = plt.subplots(1,1)
    indx = range(0,psurf.shape[0])

    #Plot line and hidden region for intrinsic interface
    #for i in range(10):
    #    ax.vlines(xp[sb],-2,2.,color='w', lw=2*i, zorder=5+i, alpha=1-0.1*i)
    #ax.vlines(xp[sb],-2,2.,color='k', zorder=30)

    dxidx = dsurf[indx,12]
    dxidy = dsurf[indx,15]
    ax.plot(xp[indx], psurf[indx,0]-0.5*(dxidx+dxidy),'-', c=c1,label=r"{CV} $\Pi^{c}_{N}$")
    ax.plot(xp[indx], 0.5*(dxidx+dxidy),'-', c=c6,label=r"x_{ij} \partial \xi / \partial x$")
    ax.plot(xp[indx], psurf[indx,0],'-', c=c2,label=r"x_{ij} \partial \xi / \partial x$")
    ax.plot(xp[indx], vflux[indx,0],'-', c=c5,label=r"{CV} $\Pi^{k}_{N}$")
    ax.plot(xp[indx], dsurf[indx,0],'-', c=c3,label=r"$\partial \xi / \partial t$")
    ax.plot(xp[indx], vflux[indx,0]+dsurf[indx,0],'-', c=c4,label=r"{CV} $\Pi^{k}_{N}+\partial \xi / \partial t$")
    ax.plot(xp[indx], psurf[indx,0]+vflux[indx,0]+dsurf[indx,0],'-', c='k',label="total conserved")
    #ax.legend()

    ax.set_xlim([-3,2.])
    ax.set_ylim([-1.1,1.5])
    ax.set_xlabel("$z$")
    ax.set_ylabel("Surface Pressure")

    plt.savefig("CV_terms.pdf", bbox_inches="tight")
    plt.savefig("CV_terms.eps", bbox_inches="tight")
    plt.show()





#####################################
#
#  Plot CV balance inc VA
#
#####################################

if 9 in showplots:

    fig, ax = plt.subplots(1,1)
    indx = range(0,psurf.shape[0])

    dxidt = interp(dsurf,0)  #0.5*(dsurf[1:,0]+dsurf[:-1,0])
    dxidx = interp(dsurf,12)  #0.5*(dsurf[1:,12]+dsurf[:-1,12])
    dxidy = interp(dsurf,15)  #0.5*(dsurf[1:,15]+dsurf[:-1,15])
    dxidr = 0.5*(dxidx+dxidy)
    ax.plot(xp[indx], psurf[indx,0]-dxidr,'-', c=c4,label=r"{CV} $\Pi^{c}_{N}$")
    ax.plot(xp[indx], dxidr,'-', c=c6,label=r"x_{ij} \partial \xi / \partial x$")
    ax.plot(xp[indx], psurf[indx,0],'-', c=c1,label=r"x_{ij} \partial \xi / \partial x$")
    ax.scatter(x[indx], pVAc[indx,0]+dxidr[indx], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax.plot(xp[indx], vflux[indx,0],'-', c=c5,label=r"{CV} $\Pi^{k}_{N}$")
    ax.plot(xp[indx], dsurf[indx,0],'-', c='k',label=r"$\partial \xi / \partial t$")
    ax.plot(xp[indx], vflux[indx,0]+dsurf[indx,0],'-', c=c2,label=r"{CV} $\Pi^{k}_{N}+\partial \xi / \partial t$")
    ax.scatter(x[indx], pVAk[indx,0]+dxidt[indx], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax.plot(xp[indx], psurf[indx,0]+vflux[indx,0]+dsurf[indx,0],'-', c=c3,label="total conserved")
    ax.scatter(x[indx], pVAc[indx,0]+dxidr[indx]+pVAk[indx,0]+dxidt[indx], c=c3, edgecolors='k', s=20.,label="total conserved")

    ax.set_xlim([-3,2.])
    ax.set_ylim([-1.1,1.5])
    ax.set_xlabel("$z$")
    ax.set_ylabel("Surface Pressure")

    plt.savefig("CV_terms_incVA.pdf", bbox_inches="tight")
    plt.savefig("CV_terms_incVA.eps", bbox_inches="tight")
    plt.show()



#####################################
#
#  Plot VA with correction only
#
#####################################

if 9.5 in showplots:

    fig, ax = plt.subplots(1,1)
    indx = range(0,psurf.shape[0])

    dxidt = interp(dsurf,0)  #0.5*(dsurf[1:,0]+dsurf[:-1,0])
    dxidx = interp(dsurf,12)  #0.5*(dsurf[1:,12]+dsurf[:-1,12])
    dxidy = interp(dsurf,15)  #0.5*(dsurf[1:,15]+dsurf[:-1,15])
    dxidr = 0.5*(dxidx+dxidy)

    ax.plot(xp[indx], dxidr,'-', c=c6,label=r"x_{ij} \partial \xi / \partial x$")

    ax.scatter(x[indx], pVAc[indx,0], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax.plot(xp[indx], dxidt,'-', c='k',label=r"$\partial \xi / \partial t$")

    ax.scatter(x[indx], pVAk[indx,0], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)

    ax.set_xlim([-3,2.])
    ax.set_ylim([-1.2,1.2])
    ax.set_xlabel("$z$")
    ax.set_ylabel("Surface Pressure")

    plt.savefig("VA_correction_breakdown.pdf", bbox_inches="tight")
    plt.savefig("VA_correction_breakdown.eps", bbox_inches="tight")
    plt.show()


if 9.75 in showplots:

    fig, ax = plt.subplots(1,2, sharey=True)
    indx = range(0,psurf.shape[0])

    dxidt = interp(dsurf,0)  #0.5*(dsurf[1:,0]+dsurf[:-1,0])
    dxidx_k = interp(dsurf,3)
    dxidy_k = interp(dsurf,6)
    dxidr_k = 0.5*(dxidx_k+dxidy_k)
    dxidx = interp(dsurf,12)
    dxidy = interp(dsurf,15)
    dxidr = 0.5*(dxidx+dxidy)


    ax[0].plot(xp[indx], -dxidr_k,'-', c=c4)

    ax[0].plot(xp[indx], -dxidr,'-', c=c6,label=r"x_{ij} \partial \xi / \partial x$")
    ax[0].scatter(x[indx], pVAc[indx,0], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax[0].plot(xp[indx], dxidt,'-', c=c5,label=r"$\partial \xi / \partial t$")
    ax[0].scatter(x[indx], pVAk[indx,0], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)

    ax[1].plot(xp[indx], interp(psurf,0)[indx],'-', c=c1,label=r"x_{ij} \partial \xi / \partial x$")
    ax[1].plot(xp[indx], interp(vflux,0)[indx]+dxidt[indx],'-', c=c2,label=r"{CV} $\Pi^{k}_{N}+\partial \xi / \partial t$")
    ax[1].plot(xp[indx], interp(psurf,0)[indx]+interp(vflux,0)[indx]+dxidt[indx],'-', c=c3,label="total conserved")


    ax[1].scatter(xp[indx], pVAc[indx,0]+dxidr[indx], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax[1].scatter(xp[indx], pVAk[indx,0]+dxidt[indx], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax[1].scatter(xp[indx], pVAc[indx,0]+dxidr[indx]+pVAk[indx,0]+dxidt[indx], c=c3, edgecolors='k', s=20.,label="total conserved")
    
    ax[0].set_ylim([-1.2,1.2])
    ax[0].set_ylabel(r"$P_N$")

    for a in ax:
        a.set_xlim([-3,2.])
        a.set_xlabel("$z$")
        a.plot(xp[indx], np.zeros(xp[indx].shape[0]), "-k", lw=2.,zorder=0)


    plt.savefig("VA_breakdown_CVvsVA.pdf", bbox_inches="tight")
    plt.savefig("VA_breakdown_CVvsVA.eps", bbox_inches="tight")
    plt.show()


#####################################
#
#  Plot VA vs CV
#
#####################################

if 10 in showplots:
    fig, ax = plt.subplots(1,1)
    indx = range(1,psurf.shape[0]-1)

    dxidt = interp(dsurf,0)  #0.5*(dsurf[1:,0]+dsurf[:-1,0])
    dxidx = interp(dsurf,12)  #0.5*(dsurf[1:,12]+dsurf[:-1,12])
    dxidy = interp(dsurf,15)  #0.5*(dsurf[1:,15]+dsurf[:-1,15])
    dxidr = 0.5*(dxidx+dxidy)
    ax.plot(xp[indx], interp(psurf,0)[indx],'-', c=c1,label=r"x_{ij} \partial \xi / \partial x$")
    ax.plot(xp[indx], interp(vflux,0)[indx]+dxidt[indx],'-', c=c2,label=r"{CV} $\Pi^{k}_{N}+\partial \xi / \partial t$")
    ax.plot(xp[indx], interp(psurf,0)[indx]+interp(vflux,0)[indx]+dxidt[indx],'-', c=c3,label="total conserved")


    ax.scatter(xp[indx], pVAc[indx,0]+dxidr[indx], c=c1, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax.scatter(xp[indx], pVAk[indx,0]+dxidt[indx], c=c2, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)
    ax.scatter(xp[indx], pVAc[indx,0]+dxidr[indx]+pVAk[indx,0]+dxidt[indx], c=c3, edgecolors='k', s=20.,label="total conserved")
    #ax.legend()

    ax.set_xlim([-3,2.])
    ax.set_ylim([-1.2,1.2])
    ax.set_xlabel("$z$")
    ax.set_ylabel("Surface Pressure")

    plt.savefig("CV_VA_corrected_terms.pdf", bbox_inches="tight")
    plt.savefig("CV_VA_corrected_terms.eps", bbox_inches="tight")
    plt.show()




#####################################
#
#  Plot Surface Tension 
#
#####################################

if 11 in showplots:



    fig, ax = plt.subplots(1,2, sharey=True)
    indx = range(0,psurf.shape[0])

    PCV = psurf[indx,:]+vflux[indx,:]
    #PCV[:,0] = interp(psurf,0)[indx] + interp(vflux,0)[indx] + dsurf[indx,0]
    PCV[:,0] = psurf[indx,0]+vflux[indx,0]+dsurf[indx,0]
    PCV_c = psurf[indx,:]
    PCV_c[:,0] = interp(psurf,0)[indx]
    pVA = pVAc[indx,:]+pVAk[indx,:]
    pIK1= pIK1c[indx,:]+pIK1k[indx,:]

    intmin = -6
    intmax = 3
    xi = xp[indx]
    intrange = (xi>intmin) & (xi<intmax)
    imin = np.argmax(xi>intmin)
    imax = np.argmin(xi<intmax)

    integrandCV = get_integrand(PCV)
    STCV = integrate(integrandCV, imin, imax, dx)
    ax[1].plot(xi[intrange], STCV, c=c2)
    ax[0].plot(xi, integrandCV, c=c2,label=r"$\int \left[P_N-P_T\right]_{CV} dx$", zorder=9)

    integrandIK1 = get_integrand(pIK1)
    STIK1 = integrate(integrandIK1, imin, imax, dx)
    ax[1].plot(xi[intrange], STIK1, c=c1)
    ax[0].plot(xi, integrandIK1, c=c1, label=r"$\int \left[P_N-P_T\right]_{VA} dx$", zorder=9)

    integrandVA = get_integrand(pVA)
    STVA = integrate(integrandVA, imin, imax, dx)
    ax[1].plot(xi[intrange], STVA, c=c3)
    ax[0].plot(xi, integrandVA, c=c3, label=r"$\int \left[P_N-P_T\right]_{VA} dx$", zorder=9)

    for i in range(ax.shape[0]):
        ax[i].plot(xi, np.zeros(xi.shape[0]), "-k")

        ax[i].set_xlim([-3,2.])
        ax[i].set_ylim([-0.6,1.6])
        ax[i].set_xlabel("$z$")

    ax[0].set_ylabel("Pressure and Surface Tension")
    #plt.legend()


    plt.savefig("surface_tension.pdf", bbox_inches="tight")
    plt.savefig("surface_tension.eps", bbox_inches="tight")
    plt.show()





#####################################
#
#  Plot Surface of Tension 
#
#####################################
if 12 in showplots:
    fig, ax = plt.subplots(1,2, sharey=True)
    indx = range(0,psurf.shape[0])

    PCV = psurf[indx,:]+vflux[indx,:]
    #PCV[:,0] = interp(psurf,0)[indx] + interp(vflux,0)[indx] + dsurf[indx,0]
    #PCV[:,0] = psurf[indx,0]+vflux[indx,0]+dsurf[indx,0]
    PCV_c = psurf[indx,:]
    PCV_c[:,0] = interp(psurf,0)[indx]
    pVA = pVAc[indx,:]+pVAk[indx,:]
    pIK1= pIK1c[indx,:]+pIK1k[indx,:]

    intmin = -3
    intmax = 2.5
    xi = xp[indx]
    intrange = (xi>intmin) & (xi<intmax)
    imin = np.argmax(xi>intmin)
    imax = np.argmin(xi<intmax)

    integrandCV = get_integrand(PCV)
    STCV = integrate(integrandCV*xi, imin, imax, dx)
    ax[1].plot(xi[intrange], STCV, c=c2)
    ax[0].plot(xi, integrandCV*xi, c=c2,label=r"$\int x \left[P_N-P_T\right]_{CV} dx$", zorder=9)

    integrandIK1 = get_integrand(pIK1)
    STIK1 = integrate(integrandIK1*xi, imin, imax, dx)
    ax[1].plot(xi[intrange], STIK1, c=c1)
    ax[0].plot(xi, integrandIK1*xi, c=c1, label=r"$\int x \left[P_N-P_T\right]_{VA} dx$", zorder=9)

    integrandVA = get_integrand(pVA)
    STVA = integrate(integrandVA*xi, imin, imax, dx)
    ax[1].plot(xi[intrange], STVA, c=c3)
    ax[0].plot(xi, integrandVA*xi, c=c3, label=r"$\int x \left[P_N-P_T\right]_{VA} dx$", zorder=9)

    for i in range(ax.shape[0]):
        ax[i].plot(xi, np.zeros(xi.shape[0]), "-k")

        ax[i].set_xlim([-3,2.])
        ax[i].set_ylim([-1.5,1.5])
        ax[i].set_xlabel("$z$")

    ax[0].set_ylabel("Pressure*z and Surface of Tension")
    plt.legend()

    plt.savefig("surface_of_tension.pdf", bbox_inches="tight")
    plt.savefig("surface_of_tension.eps", bbox_inches="tight")
    plt.show()


#####################################
#
#  Plot integral  of normal
#
#####################################


#fig, ax = plt.subplots(1,2)

# integrandCV = PCV[:,0]
# STCV = integrate(integrandCV, imin, imax, dx)
# ax[1].plot(xi[intrange], STCV, c=c2)
# ax[0].plot(xi, integrandCV, c=c2,label=r"$\int \left[P_N\right]_{CV} dx$", zorder=9)

# integrandIK1 = pIK1[:,0]
# STIK1 = integrate(integrandIK1, imin, imax, dx)
# ax[1].plot(xi[intrange], STIK1, c=c1)
# ax[0].plot(xi, integrandIK1, c=c1, label=r"$\int \left[P_N\right]_{VA} dx$", zorder=9)

# integrandVA = pVA[:,0]
# STVA = integrate(integrandVA, imin, imax, dx)
# ax[1].plot(xi[intrange], STVA, c=c3)
# ax[0].plot(xi, integrandVA, c=c3, label=r"$\int \left[P_N\right]_{VA} dx$", zorder=9)

# for i in range(ax.shape[0]):
    # ax[i].plot(xi, np.zeros(xi.shape[0]), "-k")
    # ax[i].set_xlim([-3,2.])
    # ax[i].set_xlabel("$z$")

# ax[0].set_ylim([-1.,1.])
# ax[1].set_ylim([-0.2,0.2])
# ax[0].set_ylabel("Pressure")

# ax[1].set_ylabel("Surface Tension ")

# plt.savefig("Normal_integral.pdf", bbox_inches="tight")
# plt.savefig("Normal_integral.eps", bbox_inches="tight")
# plt.show()


# fig, ax = plt.subplots(1,2)

# integrandCV = PCV[:,4]
# STCV = integrate(integrandCV, imin, imax, dx)
# ax[1].plot(xi[intrange], STCV, c=c2)
# ax[0].plot(xi, integrandCV, c=c2,label=r"$\int \left[P_T\right]_{CV} dx$", zorder=9)

# integrandIK1 = pIK1[:,4]
# STIK1 = integrate(integrandIK1, imin, imax, dx)
# ax[1].plot(xi[intrange], STIK1, c=c1)
# ax[0].plot(xi, integrandIK1, c=c1, label=r"$\int \left[P_T\right]_{VA} dx$", zorder=9)

# integrandVA = pVA[:,4]
# STVA = integrate(integrandVA, imin, imax, dx)
# ax[1].plot(xi[intrange], STVA, c=c3)
# ax[0].plot(xi, integrandVA, c=c3, label=r"$\int \left[P_T\right]_{VA} dx$", zorder=9)

# for i in range(ax.shape[0]):
    # ax[i].plot(xi, np.zeros(xi.shape[0]), "-k")
    # ax[i].set_xlim([-3,2.])
    # ax[i].set_xlabel("$z$")

# #ax[0].set_ylim([-1.,1.])
# #ax[1].set_ylim([-0.2,0.2])
# ax[0].set_ylabel("Pressure")

# ax[1].set_ylabel("Surface Tension ")

# plt.savefig("Tangent_integral.pdf", bbox_inches="tight")
# plt.savefig("Tangent_integral.eps", bbox_inches="tight")
# plt.show()






#####################################
#
#  Split surface tension integral into parts
#
#####################################


# for i, c in enumerate([0,4]):
    # integrandCV = PCV[:,c]
    # STCV = integrate(integrandCV, imin, imax, dx)
    # ax[i].plot(xi[intrange], STCV, c=c2)

    # integrandIK1 = pIK1[:,c]
    # STIK1 = integrate(integrandIK1, imin, imax, dx)
    # ax[i].plot(xi[intrange], STIK1, c=c1)

    # integrandVA = pVA[:,c]
    # STVA = integrate(integrandVA, imin, imax, dx)
    # ax[i].plot(xi[intrange], STVA, c=c3)

    # ax[i].plot(xi, np.zeros(xi.shape[0]), "-k")
    # ax[i].set_xlim([-3,2.])
    # ax[i].set_xlabel("$z$")

# #ax[0].set_ylim([-1.,1.])
# #ax[1].set_ylim([-0.2,0.2])

# ax[0].set_ylabel("Surface Tension ")

# plt.savefig("Surface_integral_components.pdf", bbox_inches="tight")
# plt.savefig("Surface_integral_components.eps", bbox_inches="tight")
# plt.show()

if 13 in showplots:

    indx = range(0,psurf.shape[0])
    PCV = psurf[indx,:]+vflux[indx,:]
    PCV[:,0] = psurf[indx,0]+vflux[indx,0]+dsurf[indx,0]
    pVA = pVAc[indx,:]+pVAk[indx,:]
    pIK1= pIK1c[indx,:]+pIK1k[indx,:]

    intmin = -5
    intmax = 2
    xi = xp[indx]
    intrange = (xi>intmin) & (xi<intmax)
    imin = np.argmax(xi>intmin)
    imax = np.argmin(xi<intmax)


    fig, ax = plt.subplots(1,2, sharey=True)

    integrandCV = PCV[:,0]
    STCV = integrate(integrandCV, imin, imax, dx)
    ax[0].plot(xi[intrange], STCV, c=c3)

    integrandIK1 = pIK1[:,0]
    STIK1 = integrate(integrandIK1, imin, imax, dx)
    ax[0].plot(xi[intrange], STIK1, '--', c=c3)

    integrandVA = pVA[:,0]
    STVA = integrate(integrandVA, imin, imax, dx)
    #ax[0].plot(xp[intrange], STVA, c=c3)
    ax[0].scatter(xp[intrange], STVA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)


    integrandCV = 0.5*(PCV[:,4]+PCV[:,8])
    STCV = integrate(integrandCV, imin, imax, dx)
    ax[1].plot(xi[intrange], -STCV, c=c3)

    integrandIK1 = 0.5*(pIK1[:,4]+pIK1[:,8])
    STIK1 = integrate(integrandIK1, imin, imax, dx)
    ax[1].plot(xi[intrange], -STIK1, '--', c=c3)

    integrandVA = 0.5*(pVA[:,4]+pVA[:,8])
    STVA = integrate(integrandVA, imin, imax, dx)
    #ax[1].plot(xp[intrange], -STVA, c=c3)
    ax[1].scatter(xp[intrange], -STVA, c=c3, edgecolors='k', s=20., label=r"{VA} $\Pi^{c}_{N}$",zorder=5)

    for i in range(2):
        ax[i].plot(xi, np.zeros(xi.shape[0]), "-k")
        ax[i].set_xlim([-3,2.])
        ax[i].set_xlabel("$z$")

    #ax[0].set_ylim([-1.,1.])
    #ax[1].set_ylim([-0.2,0.2])

    ax[0].set_ylabel("$\int_{-\infty}^{z} P(z^{\prime}) dz^{\prime}")

    plt.savefig("Surface_integral_components_symbols.pdf", bbox_inches="tight")
    plt.savefig("Surface_integral_components_symbols.eps", bbox_inches="tight")
    plt.show()





