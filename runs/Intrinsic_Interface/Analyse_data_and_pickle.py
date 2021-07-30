import matplotlib.pyplot as plt
import numpy as np
import sys
import cPickle as pickle

def load_dtype(PPObj, dtype, startrec=1):
    print("loading "+ dtype)
    Obj = PPObj.plotlist[dtype]
    endrec=Obj.maxrec-1
    x, d = Obj.profile(axis=0, 
               startrec=startrec, 
               endrec=endrec)
    return d

ppdir = '/home/es205/pyDataView'
sys.path.append(ppdir)
import postproclib as ppl

fdir = './results/'

try:
    data = pickle.load(open("summary.p","r"))
except FileNotFoundError:
    data = {}

#Density
if "rho" not in data:
    PPObj = ppl.All_PostProc(fdir)
    data["rho"] = load_dtype(PPObj, "rho")
    pickle.dump(data, open("summary.p","w+"))

#VA Data
if "pVA_c" not in data:
    PPObj = ppl.All_PostProc(fdir)
    data["pVA_c"] = load_dtype(PPObj, "pVA_c")
    data["pVA_k"] = load_dtype(PPObj, "pVA_k")
    pickle.dump(data, open("summary.p","w+"))

#CV data
if "psurface" not in data:
    PPObj = ppl.All_PostProc(fdir)
    for dtype in ["dsurf_vflux", "psurface", "vflux"]: 
        if dtype not in data:
            data[dtype] = load_dtype(PPObj, dtype)
            pickle.dump(data, open("summary.p","w+"))

fig, ax = plt.subplots(1,1)
dx = data["x"][2] - data["x"][1]
ax.plot(data["x"], data["pVA_c"][:,0], label="Volume Average Configurational")
ax.plot(data["x"], data["pVA_k"][:,0], label="Volume Average Kinetic")
ax.plot(data["x"]-dx/2., data["psurface"][:,0], label="Surface Flux Configurational")
ax.plot(data["x"]-dx/2., data["vflux"][:,0], label="Surface Flux Kinetic")

ax.plot(data["x"]-dx/2., data["dsurf_vflux"][:,0], label="dxi/dt")
ax.plot(data["x"]-dx/2., data["dsurf_vflux"][:,15], label="dxi/dy")

plt.legend()
plt.show()
