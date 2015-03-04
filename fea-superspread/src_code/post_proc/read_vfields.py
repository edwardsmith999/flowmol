import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import griddata
import glob
import re
import sys
import os

from streamlines import nonuniform_streamplot
sys.path.append('../../../utils/')
from postproclib.rawdata import RawData

#def streams(ax, xx, yy, u, v, 
#            density = None, color=None, lw=None, 
#            base_map=False, plotgrid=False, **kwargs):

#    if density == None:
#        density = 1.

#    #Streamplot defaults to 30 x 30 seeds for density of 1
#    if isinstance(density,(int,float)):
#        res = [30*density, 30*density]
#    elif isinstance(density,list) and len(density) == 2:
#        res = [30*i for i in density]
#    else:
#        raise ValueError("density should be float or 2-tuple")

#    x = np.linspace(xx.min(), xx.max(), res[0])
#    y = np.linspace(yy.min(), yy.max(), res[1])
#    xi, yi = np.meshgrid(x,y)

#    #then, interpolate your data onto this grid:
#    px = xx.flatten()
#    py = yy.flatten()
#    pu = u.flatten()
#    pv = v.flatten()

#    gu = griddata(zip(px,py), pu, (xi,yi))
#    gv = griddata(zip(px,py), pv, (xi,yi))

#    #Get speed for color and linewidth
#    if lw == None or color == None or plotgrid:
#        speed = np.sqrt(u*u + w*w)
#        if lw == None:
#            pspeed = speed.flatten()
#            gspeed = griddata(zip(px,py), pspeed, (xi,yi))
#            linewidth = 6*gspeed/np.nanmax(gspeed)
#            lw = linewidth
#        if color == None:
#            color = gspeed

#    #now, you can use x, y, gu, gv and gspeed in streamplot:
#    if base_map:
#        xx,yy = ax(xx,yy)
#        xi,yi = ax(xi,yi)

#    #Add grid to plot
#    if plotgrid:
#        ax.contour(xx,yy,speed, colors='k', alpha=0.4)
#        ax.plot(xx,yy,'-k',alpha=0.3)
#        ax.plot(xx.T,yy.T,'-k',alpha=0.3)
#        ax.plot(xi,yi,'-b',alpha=0.1)
#        ax.plot(xi.T,yi.T,'-b',alpha=0.1)

#    #Call streamplot with mapped grid
#    c = ax.streamplot(x, y, gu, gv, color=color, 
#                      linewidth=lw, density=density, **kwargs) 

#    return c


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

        class HeaderData:
            def __init__(self,fdir):
                headerfile = fdir + './SPARAMETERS'
                with open(headerfile) as f:
                    dump = f.read()
                    data = dump.split('\n')
                    for i in data:
                        if (i.find('=')>0):
                            nocomment = i.split('!')[0]
                            nocomment = nocomment.replace(" ","")
                            varname, varval = nocomment.split('=')
                            vars(self)[varname] = varval.replace('d','e')

                    self.nx = int(self.NXELa) + int(self.NXELb) + 1
                    self.ny = 1
                    self.nz = int(self.NZELa)

        return HeaderData(fdir)

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
            nx = int(self.header.nx)
        if (ny == None):
            ny = int(self.header.ny)
        if (nz == None):
            nz = int(self.header.nz)
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
                    print('Filename = ', fpath, ' is record ', plusrec, ' of ', nrecs)
            except IndexError:
                raise DataNotAvailable
            try:
                with open(fpath,'rb') as fobj:
                    data = np.fromfile(fobj,dtype=self.dtype)
            except:
                raise

            bindata = np.reshape( data,[nx,ny,nz,npercell], order='F')

            # insert into array (N.B. calling this grid for the data 
            #                    on a grid is a terrible name...)
            grid[:,:,:,plusrec,:] = bindata[lower[0]:upper[0],
                                            lower[1]:upper[1],
                                            lower[2]:upper[2], :]

         
        return grid


if __name__ == "__main__":

    fdir = '/home/es205/codes/superspreading/coupled_code/FEA_superspread/src_code/results/HLratio_study/eps2_0p18/'
    #fdir = '/home/es205/codes/superspreading/coupled_code/FEA_superspread/src_code/results/HLratio_study/eps2_0p1'
    #fdir = '/home/es205/codes/superspreading/coupled_code/FEA_superspread/src_code/results/HLratio_study/eps2_0p02/'
    fname = 'uwgrid'
    truescale = True
    fObj  = FEA_RawData(fdir,fname,'d',2)

    print('Number of records = ', fObj.maxrec)
    n = 0
    startrec = 1; endrec = fObj.maxrec; skiprec = 10
    plots = int(np.ceil((endrec-startrec)/skiprec))+1
    ratio = 2.
    f,ax = plt.subplots(int(plots/ratio),int(ratio))
    axs = ax.reshape(plots) 
    for rec in range(startrec,endrec,skiprec):
        print(rec,plots)

        grid = fObj.get_gridtopology(rec,rec)
        uw = fObj.read(rec,rec,verbose=True)
        print(rec)

        x = grid[0][:,0,:,0,0]
        z = grid[1][:,0,:,0,0]
        eps = np.sqrt(float(fObj.header.eps2))
        if truescale:
            H = 1.0
            L = H/eps
        else:
            L = 1.0
            H = 1.0
        x = x * L
        z = z * H
        X,Z = np.meshgrid(x,z[0,:],indexing='ij')
        u = uw[:,0,:,0,0]
        w = uw[:,0,:,0,1]
        speed = np.sqrt(u*u + w*w)
        size = [10.*z[i,:].max() for i in range(z.shape[0])]
        #print([z[i,:].max() for i in range(z.shape[0])])
        #scaled_u = (u - u.min()) / u.ptp()
        #colors = plt.cm.RdYlBu_r(scaled_u)


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
        skipx = 50; skipz = 16

        #axs[n].quiver(X[::skipx,::skipz], z[::skipx,::skipz], u[::skipx,::skipz],
        #             -w[::skipx,::skipz])#,angles='xy',scale_units='xy',scale=0.1)
        cs = nonuniform_streamplot(axs[n], X, z, u, -w, density=1., lw = 1., color='k')# cmap=plt.cm.RdYlBu_r)
        cm = axs[n].pcolormesh(X, z, speed, cmap=plt.cm.RdYlBu_r,vmin=0.000001,vmax = 0.1,norm=matplotlib.colors.LogNorm(),alpha=0.9)
        axs[n].plot(x,np.zeros(x.shape[0]),'k-')
        axs[n].axis('scaled')
        axs[n].set_xlim((0.0,1.8*L))
        axs[n].set_ylim((0.0,1.0))


        #Streamplot only works on a uniform grid!!
        #f, ax = plt.subplots(nrows=1)
        #lw = 5.*speed/speed.max()
        #sp = streamplot(ax, Z[::skipx,::skipz], X[::skipx,::skipz],u[::skipx,::skipz], -w[::skipx,::skipz], 
        #                    linewidth=lw[::skipx,::skipz])
        #im=axs[n].streamplot(Z[::skipx,::skipz], X[::skipx,::skipz], u[::skipx,::skipz], -w[::skipx,::skipz])#, 
                             #linewidth=lw[::skipx,::skipz], density=[2., 2])
        #cm = im.lines
        #im=plt.streamplot(Z[::skipx,::skipz], X[::skipx,::skipz], u[::skipx,::skipz], w[::skipx,::skipz])
        #plt.colorbar(im.lines)

        #Update axis number
        n += 1

    plt.suptitle("eps2 = " + fObj.header.eps2)
    f.subplots_adjust(left=0.2)
    cbar_ax = f.add_axes([0.05, 0.1, 0.025, 0.8])
    try:
        f.colorbar(cs.lines, cax=cbar_ax)
    except:
        pass
    try:
        f.colorbar(cm, cax=cbar_ax)
    except:
        pass
    plt.show()

