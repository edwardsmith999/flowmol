#! /usr/bin/env python
import numpy as np
import os

class CFD_RawData:
    
    def __init__(self,fdir):
        self.fdir = fdir
        self.grid = self.get_grid()
        self.subdomlist = self.get_subdomlist()
        self.npercell = self.get_npercell()
        self.maxrec = len(self.subdomlist)-1 # count from 0

    def get_grid(self):
        fobj = open(self.fdir+'report','r')
        report = fobj.readlines()[3:6] # Lines with info in
        for line in report:
            linepairs = line.split('|')
            for pair in linepairs:
                varname = pair.split()[0]
                varval  = pair.split()[1]
                vars(self)[varname] = varval    
        # Number of grid points in main code
        self.nx = int(self.nx)
        self.ny = int(self.ny)
        self.nz = int(self.nz)
        # Number of cell-centered values written to files
        # -1 for cell centers rather than grid points
        # -2 for not writing halos (except in y-direction)
        # Therefore -3 in x and z, -1 in y
        self.nrx = int(self.nx)-3 # number of subdom grid records in x
        self.nry = int(self.ny)-3+2 # +2 halos
        self.nrz = int(self.nz)-3
        # Domain lengths
        self.xL = float(self.xL)
        self.yL = float(self.yL)
        self.zL = float(self.zL)
        # Grid spacing
        self.dx = self.xL/float(self.nx-3)
        self.dy = self.yL/float(self.ny-3)
        self.dz = self.zL/float(self.nz-3)
        # Linspaces of cell centers, accounting for halos written in y
        gridx = np.linspace( self.dx/2., self.xL -self.dx/2., num=self.nrx)
        gridy = np.linspace(-self.dy/2., self.yL +self.dy/2., num=self.nry)
        gridz = np.linspace( self.dz/2., self.zL -self.dz/2., num=self.nrz)
        grid = [gridx,gridy,gridz]
        # Debug info
        #print('nxyz = ', [self.nx, self.ny, self.nz])
        #print('nrxyz = ', [self.nrx, self.nry, self.nrz])
        #print('xyzL = ',[self.xL, self.yL, self.zL])
        #print('dxyz = ',[self.dx, self.dy, self.dz])
        #print(grid)
        return grid 

    def get_subdomlist(self):

        def get_int(name):
            string, integer = name.split('.')
            return int(integer)

        subdoms = []
        for filename in os.listdir(self.fdir):
            if (filename.find('SubDom') != -1):
                subdoms.append(filename)

        if (len(subdoms) == 0):
            raise IOError 

        subdoms = sorted(subdoms,key=get_int)
        return subdoms

    def get_npercell(self):
        dprealbytes = 8 # 8 for dp float
        ngridpoints = self.nrx * self.nry * self.nrz
        filepath = self.fdir + self.subdomlist[0]
        filesize = os.path.getsize(filepath)
        npercell = filesize / (dprealbytes*ngridpoints) 
        return npercell

    def read(self,startrec,endrec):

        nrecs = endrec - startrec + 1
        # Efficient memory allocation
        subdata = np.empty((self.nrx,self.nry,self.nrz,nrecs,self.npercell))

        # Loop through files and insert data
        for plusrec in range(0,nrecs):

            fpath = self.fdir + self.get_subdomlist().pop(startrec+plusrec)
            with open(fpath,'rb') as fobj:
                data = np.fromfile(fobj,dtype='d')
                # zxy ordered in file
                data = np.reshape(data,[self.nrz,self.nrx,self.nry,self.npercell],
                                  order='F')
                # change to xyz ordering
                data = np.transpose(data,(1,2,0,3))
                # insert into array
                subdata[:,:,:,plusrec,:] = data 
         
        return subdata
