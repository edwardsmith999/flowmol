#! /usr/bin/env python
import numpy as np
import subprocess as sp
import os

class Channelflow_RawData:
    
    def __init__(self,fdir):
        self.fdir = fdir
        self.get_channelflow_utils()
        self.subdomlist = self.get_subdomlist()
        self.grid = self.get_grid()
        self.maxrec = len(self.subdomlist)-1 # count from 0

    def get_channelflow_utils(self):

        """
            Find dirctory of compiled code and
            save a range of routine used for postprocessing
        """

        # Routine for converting spectral flowfied (.h5 format) to spectral flowfield (.ff format)
        self.fieldconvert = self.__syslocate('fieldconvert')
        # Routine for converting spectral flowfied (.h5 format) to realspace ascii and geom (.asc and .geom format)
        self.ascii2field = self.__syslocate('ascii2field')
        # Routine to get field geometry from spectral flowfied (.h5 format)
        self.fieldprops = self.__syslocate('fieldprops')

    def __syslocate(self,name):
        locout = sp.check_output(['locate',str(name)])
        return locout.split('\n')[0]

    def get_grid(self):
        """
            Get details for CFD grid from 1st subdomain
        """

        self.filename = self.fdir + self.subdomlist[0] 
        try:
            rawgeomdata = sp.check_output([self.fieldprops,'-g', self.filename])
        except:
            raise

        # Split into list and keep only variables
        variables = [x for x in rawgeomdata.split('\n') if x.find('==') != -1]

        #Evaluate and store as class variables
        for x in variables:
            try:
                exec('self.'+ x.replace('==','=').strip(' '))
            except:
                print("Can't save" + ' self.'+ x.replace('==','=').strip(' '))

        # Number of grid points in main code
        self.nx = int(self.Nx)
        self.ny = int(self.Ny)
        self.nz = int(self.Nz)
        # Number of cell-centered values written to files
        # -1 for cell centers rather than grid points
        # -2 for not writing halos (except in y-direction)
        # Therefore -3 in x and z, -1 in y
        self.nrx = int(self.Nx) # number of subdom grid records in x
        self.nry = int(self.Ny) # +2 halos
        self.nrz = int(self.Nz)
        # Domain lengths
        self.xL = float(self.Lx)
        self.yL = float(self.Ly)
        self.zL = float(self.Lz)
        # Grid spacing
        self.dx = self.xL/float(self.nx)
        self.dy = self.yL/float(self.ny)
        self.dz = self.zL/float(self.nz)
        # Linspaces of cell centers, accounting for halos written in y
        gridx = np.linspace( self.dx/2., self.xL -self.dx/2., num=self.nrx)
        gridy_linear = np.linspace(-self.dy/2., self.yL +self.dy/2., num=self.nry)
        gridy = self.linear2cosinegrid(gridy_linear)
        gridz = np.linspace( self.dz/2., self.zL -self.dz/2., num=self.nrz)
        grid = [gridx,gridy,gridz]

        return grid

    def get_subdomlist(self):

        def get_int(name):
            name = name.replace('u','')
            integer, extension = name.split('.')
            return int(integer)

        subdoms = []
        for filename in os.listdir(self.fdir):
            if (filename.find('.h5') != -1):
                subdoms.append(filename)

        subdoms = sorted(subdoms,key=get_int)
        return subdoms

    def convert(self,filename):
        """
            Convert *.h5 format to ascii *.ff format
        """
        filebase = filename.replace('.h5','')
        if (os.path.isfile(filebase + '.asc')):
            pass
        else: 
            syscall = self.ascii2field + ' -p ' + filebase + '.h5 ' + filebase + '.asc'
            print("Output file in ascii format not founts, " + 
                  "attempting to convert using field2ascii routine")
            try:
                os.system(syscall)
            except:
                raise

        return filebase + '.asc'

    def linear2cosinegrid(self,lingrid):

        """
           Convert a linear grid (from linspace) to a cosine grid
        """
        cosgrid = np.empty(len(lingrid))
        for i in range(0,len(lingrid)):
            cosgrid[i] = -np.cos((lingrid[i]*np.pi)/(self.Ny-1))

        return cosgrid
 

    def read(self,startrec,endrec, verbose=False, quit_on_error=True):

        nrecs = endrec - startrec + 1
        # Efficient memory allocation
        subdata = np.empty((self.nrx,self.nry,self.nrz,nrecs,3))

        # Loop through files and insert data
        for plusrec in range(0,nrecs):

            fpath = self.fdir + self.get_subdomlist().pop(startrec+plusrec)
            data = self.read_field(fpath)

            # insert into array
            subdata[:,:,:,plusrec,:] = data 
         
        return subdata

    # Read channelflow field
    def read_field(self,fpath):

        # Efficient memory allocation
        data = np.empty((self.nrz,self.nrx,self.nry,3))

        fileasc = self.convert(fpath)
        with open(fileasc,'r') as f:
            for nx in range(self.nx):
                for ny in range(self.ny):
                    for nz in range(self.nz):
                        try:
                            data[nz,nx,ny,0] = float(f.readline())
                            data[nz,nx,ny,1] = float(f.readline())
                            data[nz,nx,ny,2] = float(f.readline())
                        except:
                            data[nz,nx,ny,0] = 0.0
                            data[nz,nx,ny,1] = 0.0
                            data[nz,nx,ny,2] = 0.0
        return np.transpose(data,(1,2,0,3))


