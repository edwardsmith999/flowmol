import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import glob
import re
import sys
import os

from streamlines import streamplot
sys.path.append('../../../utils/')
from postproclib.rawdata import RawData



class FEA_RawData(RawData):
    
    def __init__(self,fdir,fname,dtype,npercell):

        """
            fdir       -  file directory containing results, string
            fname      -  file path from which to read raw data, string
            dtype      -  datatype string, 'i' for integer, 'd' for float
            npercell    -  number of items to read per bin, integer
        """

        if (fdir[-1] != '/'): fdir += '/' 
        self.fdir = fdir
        self.fname = fname
        self.dtype = dtype
        self.npercell = npercell
        
        self.header = self.read_header(fdir)
        self.grid = self.get_gridtopology(0,0)
        self.maxrec = self.get_maxrec()

    def read_header(self,fdir):
        
        self.nx = 601
        self.ny = 1
        self.nz = 100       

    def get_maxrec(self):
        return len(self.get_subdomlist(self.fname))

    def get_subdomlist(self,fname):

        def get_int(name,fname):
            name = name.replace(fname+'.','')
            try:
                integer, extension = name.split('.')
                integer = int(integer)
            except ValueError:
                print("Unrecognised fileno:", name.split('.'))
                integer = 10000000000

            return integer

        subdoms = []
        for filename in os.listdir(self.fdir):
            if (filename.find(fname) != -1):
                subdoms.append(filename)

        subdoms = sorted(subdoms,key=lambda x: get_int(x,fname))

        return subdoms
      

    def get_gridtopology(self, startrec, endrec, **kwargs):

        xgrid = self.read_fname(startrec, endrec, 'xgrid', nz = 1, npercell=1, **kwargs) 
        zgrid = self.read_fname(startrec, endrec, 'zgrid', npercell=1, **kwargs)
        
        return [xgrid,zgrid]

    def read(self, startrec, endrec, **kwargs):

        uw = self.read_fname(startrec, endrec, self.fname, **kwargs)
        
        return uw

    def read_fname(self, startrec, endrec, fname, nx=None, ny=None, nz=None,
                   npercell=None,binlimits=None, verbose=False, 
                        quit_on_error=True,wallnormaldir=1):

        """
            Get details for specified file extension of the form
                            fname.0000000.DAT
            where nx,ny,nz or npercell can be overidden for 1D case
        """

        if (nx == None):
            nx = self.nx
        if (ny == None):
            ny = self.ny
        if (nz == None):
            nz = self.nz
        if (npercell == None):
            npercell = self.npercell
        
        nrecs = endrec - startrec + 1
        nxyz = [nx, ny, nz]
        lower = np.empty(3); upper = np.empty(3) 

        self.fname_records  = self.get_subdomlist(fname)

        if nrecs > len(self.fname_records):
            print('Number of records ', nrecs , ' greater than ', len(self.fname_records) ,
                 ' available:', self.fname_records)
        elif startrec + nrecs > len(self.fname_records):
            print('Range of records ', startrec, ' to ', endrec , ' outside ', 
                  startrec, ' to ', len(self.fname_records)+startrec,
                  ' available:', self.fname_records)

        # If bin limits are specified, return only those within range
        for axis in range(3):
            if (binlimits):
                if (binlimits[axis] == None):
                    lower[axis] = 0
                    upper[axis] = int(nxyz[axis])
                else:
                    lower[axis] = int(binlimits[axis][0]) 
                    upper[axis] = int(binlimits[axis][1])
            else:
                lower[axis] = 0
                upper[axis] = int(nxyz[axis])

        # Efficient memory allocation
        grid = np.empty((upper[0]-lower[0],
                         upper[1]-lower[1],
                         upper[2]-lower[2],nrecs,self.npercell))

        data =  np.empty((nx*ny*nz*nrecs*npercell))

        # Loop through files and insert data
        for plusrec in range(0,nrecs):

            try:
                fpath = self.fdir + self.fname_records[startrec+plusrec]
                if verbose:
                    print('Filename = ', fpath)
            except IndexError:
                raise DataNotAvailable
            try:
                with open(fpath,'rb') as fobj:
                    data = np.fromfile(fobj,dtype=self.dtype)
            except:
                raise

            bindata = np.reshape( data,[nx,ny,nz,npercell], order='F')

            # insert into array
            grid[:,:,:,plusrec,:] = bindata[lower[0]:upper[0],
                                            lower[1]:upper[1],
                                            lower[2]:upper[2], :]

         
        return grid


if __name__ == "__main__":

    fdir = '/home/es205/codes/superspreading/coupled_code/FEA_superspread/src_code/results/'
    fname = 'uwgrid'
    fObj  = FEA_RawData(fdir,fname,'d',2)

    print('Number of records = ', fObj.maxrec)
    for rec in range(10,fObj.maxrec,10):

        grid = fObj.get_gridtopology(rec,rec)
        uw = fObj.read(rec,rec,verbose=True)
        print(rec)

        x = grid[0][:,0,:,0,0]
        z = grid[1][:,0,:,0,0]
        X,Z = np.meshgrid(x,z[0,:],indexing='ij')
        u = uw[:,0,:,0,0]
        w = uw[:,0,:,0,1]
        speed = np.sqrt(u*u + w*w)
        size = [10.*z[i,:].max() for i in range(z.shape[0])]
        #print([z[i,:].max() for i in range(z.shape[0])])
        #scaled_u = (u - u.min()) / u.ptp()
        #colors = plt.cm.RdYlBu_r(scaled_u)

        plt.plot(x,np.zeros(x.shape[0]),'k-')
        #for i in range(0,z.shape[1]):
        #    cm = plt.scatter(x[:,0], z[:,i],  s=size, c=u[:,i], marker='o', lw = 0,  alpha=0.7, cmap=plt.cm.RdYlBu_r)

        #Check mass continity is satisfied
#        dx = np.gradient(x[:,0])
#        dz = np.gradient(z)
#        dudx, dudz = np.gradient(u.T,dx)
#        dudx = dudx.T
#        dwdx, dwdz = np.gradient(w,dz)

#        f, axs = plt.subplots(nrows=2)
#        cm = axs[0].pcolormesh(X, z, -dwdz[0], cmap=plt.cm.RdYlBu_r)
#        f.colorbar(cm)
#        cm = axs[1].pcolormesh(X, z, dudx, cmap=plt.cm.RdYlBu_r)
#        f.colorbar(cm)
#        plt.show()

        #NOTE THE -w here is not in line with expected!
        skipx = 5; skipz = 5
#        cm = plt.pcolormesh(X[::skipx,::skipz], z[::skipx,::skipz], speed[::skipx,::skipz], cmap=plt.cm.RdYlBu_r, vmin=0.0, vmax=0.001)
#        plt.colorbar(cm)
#        plt.quiver(X[::skipx,::skipz], z[::skipx,::skipz], u[::skipx,::skipz], -w[::skipx,::skipz], scale=0.1,minlength=0.1)


        #Streamplot only works on a uniform grid!!
        f, ax = plt.subplots(nrows=1)
        lw = 2.5*speed/speed.max()
        sp = streamplot(ax, Z[::skipx,::skipz], X[::skipx,::skipz],u[::skipx,::skipz], -w[::skipx,::skipz], 
                            linewidth=lw[::skipx,::skipz])
        #im=plt.streamplot(Z[::skipx,::skipz], X[::skipx,::skipz], u[::skipx,::skipz], w[::skipx,::skipz], 
        #                 linewidth=lw[::skipx,::skipz], density=[3., 3])
        #im=plt.streamplot(Z[::skipx,::skipz], X[::skipx,::skipz], u[::skipx,::skipz], w[::skipx,::skipz])
        #plt.colorbar(im.lines)

        plt.show()

