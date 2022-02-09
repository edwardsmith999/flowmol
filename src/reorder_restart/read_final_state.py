
import numpy as np
import os
import struct

class final_state:

    def __init__(self, fname= "./final_state", tether_tags = [3,5,6,7,10], verbose=False):
        self.fname = fname
        self.tether_tags = tether_tags

        #Get filesize and read headersize
        self.size = os.path.getsize(fname)
        self.headersize = np.fromfile(fname, dtype=np.int64, offset=self.size-8)
        with open(fname, "rb") as f:
            f.seek(self.headersize[0])
            self.binaryheader = f.read()

        self.read_header(verbose=verbose)

    def read_header(self, verbose=False):

        #Assume 14 doubles and work out integers
        ndbls = 14; ntrlint=4
        noints = int((len(self.binaryheader) - ndbls*8)/4)-ntrlint
        self.fmtstr = str(noints) + "i"+str(ndbls) +"d"+ str(ntrlint) + "i"
        self.hdata = list(struct.unpack(self.fmtstr, self.binaryheader))
        self.htypes = ["globalnp", "initialunits1", 
                      "initialunits2", "initialunits3", 
                      "Nsteps", "tplot", "seed1", "seed2",
                      "periodic1", "periodic2", "periodic3",
                      "potential_flag","rtrue_flag","solvent_flag",
                      "nmonomers","npx","xpy","npz"]

        nproc = int(self.hdata[15])*int(self.hdata[16])*int(self.hdata[17])
        self.nproc = nproc
        [self.htypes.append("procnp"+str(p)) for p in range(nproc)]
        [self.htypes.append("proctethernp"+str(p)) for p in range(nproc)]
        [self.htypes.append(i) for i in 
                        ["globaldomain1", "globaldomain2",
                        "globaldomain3", "density", "rcutoff",
                        "delta_t", "elapsedtime", "simtime", 
                        "k_c","R_0", "eps_pp", "eps_ps", 
                        "eps_ss", "delta_rneighbr",
                        "mie_potential","global_numbering",
                        "headerstart","fileend"]]

        self.headerDict = {}
        for i in range(len(self.hdata)):
            self.headerDict[self.htypes[i]]=self.hdata[i]

        if verbose:
            for k, i in self.headerDict.items():
                print(k,i)

        #Read integers in header
#        self.dtype = np.int32
#        self.dtypesize = np.dtype(self.dtype).itemsize
#        headerdata = np.fromfile(self.fname, dtype=self.dtype, offset=self.headersize)
#        self.initialunits = np.zeros(3); self.periodic = np.zeros(3); self.seed = np.zeros(2)
#        self.N = headerdata[0]
#        self.initialunits[0] = headerdata[1]
#        self.initialunits[1] = headerdata[2]
#        self.initialunits[2] = headerdata[3]
#        self.Nsteps = headerdata[4]              #Number of computational steps
#        self.tplot  = headerdata[5]              #Frequency at which to record results
#        self.seed[0]  = headerdata[6]            #Random number seed value 1
#        self.seed[1]  = headerdata[7]            #Random number seed value 2
#        self.periodic[0] = headerdata[8]          #Boundary condition flags
#        self.periodic[1] = headerdata[9]          #Boundary condition flags
#        self.periodic[2]    = headerdata[10]       #Boundary condition flags
#        self.potential_flag = headerdata[11]       #Polymer/LJ potential flag
#        self.rtrue_flag  = headerdata[12]       #
#        self.solvent_flag  = headerdata[13]     #Solvent on/off flag
#        self.nmonomers= headerdata[14]               #Polymer chain length
#        self.npx  = headerdata[15]              #Processors (npx) for new topology
#        self.npy  = headerdata[16]              #Processors (npy) for new topology
#        self.npz   = headerdata[17]             #Processors (npz) for new topology
#        nproc = self.npx*self.npy*self.npz
#        self.nproc = nproc
#        self.procnp = np.zeros(nproc); self.proctethernp = np.zeros(nproc)
#        for p in range(nproc):
#            self.procnp[p] = headerdata[17+p+1] #Number of molecules per processors
#        for p in range(nproc):
#            self.proctethernp[p] = headerdata[17+nproc+p+1]    #Number of tethered molecules per processors
#        self.mie_potential = headerdata[-4]        
#        self.global_numbering = headerdata[-3]  

#        if verbose:
#            print(headerdata)  

#        #Read the double precision part
#        dpheaderdata = np.fromfile(self.fname, dtype=np.double, 
#                                   offset=self.headersize+(18+2*nproc)*self.dtypesize)

#        if verbose:
#            print(dpheaderdata) 
#        self.globaldomain = np.zeros(3)
#        self.globaldomain[0] = dpheaderdata[0]
#        self.globaldomain[1] = dpheaderdata[1]
#        self.globaldomain[2] = dpheaderdata[2]
#        self.density = dpheaderdata[3]
#        self.rcutoff = dpheaderdata[4]
#        self.delta_t = dpheaderdata[5]
#        self.elapsedtime    = dpheaderdata[6]
#        self.simtime = dpheaderdata[7]
#        self.k_c = dpheaderdata[8]
#        self.R_0 = dpheaderdata[9]
#        self.eps_pp  = dpheaderdata[10]
#        self.eps_ps  = dpheaderdata[11]
#        self.eps_ss  = dpheaderdata[12]
#        self.delta_rneighbr = dpheaderdata[13]

    def read_moldata(self):

        #Read the rest of the data
        data = np.fromfile(self.fname, dtype=np.double, count=int(self.headersize/8))

        #Allocate arrays
        h = self.headerDict
        N = h["globalnp"]#self.N
        self.tag = np.zeros(N)
        self.r = np.zeros([N,3])
        self.v = np.zeros([N,3])
        self.rtether = np.zeros([N,3])
        self.Ntethered = 0

        #Create arrays for molecular removal
        self.Nnew = N
        self.delmol = np.zeros(N)
        self.molecules_deleted=False

        if (h["rtrue_flag"]):
            self.rtrue = np.zeros([N,3])
        if (h["mie_potential"]):
            self.moltype = np.zeros(N)
        if (h["global_numbering"]):
            self.globnum = np.zeros(N)
        if (h["potential_flag"]):
            self.potdata = np.zeros([N,8])

        i = 0
        for n in range(N):
            self.tag[n] = data[i]; i += 1
            self.r[n,:] = data[i:i+3]; i += 3
            self.v[n,:] = data[i:i+3]; i += 3

            if (h["rtrue_flag"]):
                self.rtrue[n,:] = data[i:i+3]; i += 3
            if (self.tag[n] in self.tether_tags):
                self.rtether[n,:] = data[i:i+3]; i += 3
                self.Ntethered += 1
            if (h["mie_potential"]):
                self.moltype[n] = data[i]; i += 1
            if (h["global_numbering"]):
                self.globnum[n] = data[i]; i += 1
            if (h["potential_flag"]):
                self.potdata[n,:] = data[i:i+8]; i += 8

        return self.tag, self.r, self.v

    def plot_molecules(self, ax=None):

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
        ax.scatter(self.r[:,0], self.r[:,1], self.r[:,2], c=self.tag[:])

    def remove_molecules(self, rpos, radius, rdim=0):

        h = self.headerDict
        N = h["globalnp"]

        rmapped = np.zeros(3)
        for n in range(N):
            rmapped[:] = self.r[n,:] - rpos[:]
            #Set zero along direction
            if (rdim != 3):
                rmapped[rdim] = 0.
            #Spherical or cylindrical radius
            rspherical2 = np.dot(rmapped,rmapped)    #Get position in spherical coordinates
            rspherical = np.sqrt(rspherical2)
            #theta = np.acos(rmapped[2]/rspherical)
            #phi = np.atan(rmapped[1]/rmapped[0])
            if (rspherical < radius and self.delmol[n] != 1):
                print(n, self.Nnew, rspherical, radius)
                self.delmol[n] = 1           
                self.Nnew -= 1                   
                self.molecules_deleted = True

    def write_moldata(self, outfile=None, verbose=False):

        #Default to same filename with a 2
        if (outfile is None):
            outfile = self.fname + "2"

        h = self.headerDict
        N = h["globalnp"]

        #Values are the number of values per molecule including all 
        vals = (7 + 3*h["rtrue_flag"] + h["mie_potential"]
                + h["global_numbering"] + 8*h["potential_flag"])
        data = np.zeros(N*vals+ 3*self.Ntethered)

        #Start a new global numbering if any molecules have been deleted
        if (self.molecules_deleted):
            newglob = 1

        #Loop and write all data
        i = 0
        for n in range(N):

            if self.delmol[n] == 1:
                continue

            data[i] = self.tag[n]; i += 1
            data[i:i+3] = self.r[n,:]; i += 3
            data[i:i+3] = self.v[n,:]; i += 3
            #print(n, i, data[i-7:i])

            if (h["rtrue_flag"]):
                data[i:i+3] = self.rtrue[n,:]; i += 3
            if (tag[n] in self.tether_tags):
                data[i:i+3] = self.rtether[n,:]; i += 3
            if (h["mie_potential"]):
                data[i] = self.moltype[n]; i += 1
            if (h["global_numbering"]):
                if (self.molecules_deleted):
                    data[i] = newglob; newglob += 1; i += 1
                else:
                    data[i] = self.globnum[n]; i += 1
            if (h["potential_flag"]):
                data[i:i+8] = self.potdata[n,:]; i += 8

        #Write data to file
        data.tofile(open(outfile, "w+"))

        #If number of molecules has changed, reset to 1x1x1 processors
        if (self.Nnew != h["globalnp"]):
            print("N=", N, "Nnew=", self.Nnew)
            h["globalnp"] = self.Nnew
            h["npx"] = 1; h["npy"] = 1; h["npz"] = 1
            h["procnp0"] = self.Nnew
            proctethernp = 0
            for p in range(self.nproc):
                proctethernp += h["proctethernp"+str(p)]
            h["proctethernp0"] = proctethernp
        #Update hdata
        for i in range(len(self.hdata)):
            if (verbose or self.hdata[i] != self.headerDict[self.htypes[i]]):
                    print("UPDATE", i, self.htypes[i], "before=", self.hdata[i], 
                          "after=", self.headerDict[self.htypes[i]])
            self.hdata[i] = self.headerDict[self.htypes[i]]

        #Update binaryheader
        binaryheader = struct.pack(self.fmtstr, *self.hdata)

        #Write header at end of file
        #self.size = os.path.getsize(outfile)
        with open(outfile, "ab") as f:
            #f.seek(self.headersize[0])
            f.write(binaryheader)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d

    fig = plt.figure(); ax = []
    ax.append(fig.add_subplot(1,2,1,projection='3d'))
    ax.append(fig.add_subplot(1,2,2,projection='3d'))

    fs = final_state("./final_state", verbose=False)
    tag, r, v = fs.read_moldata()
    fs.plot_molecules(ax[0])
    
    fs.remove_molecules([0.,0.,0.],4,0)
    fs.write_moldata("./final_state_hole", verbose=False)

    fst = final_state("./final_state_hole", verbose=True)
    tag, r, v = fst.read_moldata()

    fst.plot_molecules(ax[1])
    plt.show()
