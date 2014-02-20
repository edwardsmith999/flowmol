#! /usr/bin/env python
import numpy as np

class Field():

    """
        Base class to be inherited by MDField, CFDField and CPLField.

        Authors: David Trevelyan & Ed Smith 2014

        Abstract base class template that generally specifies how data
        should be processed and returned in postprocessing routines.
        Field should be inherited by MDField (which will then, in turn,
        be inherited by MD_mField for mass, etc.). 

        Fields are, essentially, data reformatters. They must be instantiated
        with a RawData object, which is used to do the actual reading, and
        the methods in Field re-package that data in an easily plottable
        format.
    
        TOP-LEVEL EXAMPLE
        
            v = MD_vField(fdir)
            y, v3 = v.profile(axis=1)
            plt.plot(v3[:,0], y)
        
        This will instantiate a velocity field (see MDFields for details of 
        the complex inheritance and containment), read and store the data
        from multiple files to construct a velocity profile, and finally
        will plot the x-component of velocity against the y-axis.

    """
    
    def __init__(self, RawDataObj):
        self.Raw = RawDataObj
        self.fdir = self.Raw.fdir
        self.grid = self.Raw.grid
        self.maxrec = self.Raw.maxrec
        self.header = self.Raw.header

    def read(self,startrec,endrec):

        """
            TO BE OVERRIDDEN IN COMPLEX FIELDS.
            Method that returns grid data that is read by Raw.
            
        """
        if (endrec > self.maxrec):
            quit('Record ' + str(endrec) + ' is greater than the maximum '
                 'available (' + str(self.maxrec) + '). Aborting.')

        grid_data = self._read(startrec,endrec)
        return grid_data
    
    def _read(self,startrec,endrec):
        grid_data = self.Raw.read(startrec,endrec)
        return grid_data
    
    def averaged_data(self,startrec,endrec,avgaxes=(),avgtime=True):

        """
            TO BE OVERRIDDEN IN COMPLEX FIELDS.
            Average the data in the user-specified way.
        """ 
        # Read 4D time series from startrec to endrec
        grid_data = self.read(startrec,endrec)

        # Time axis is the final (fourth) axis...
        if (avgtime):
            grid_data = np.mean(grid_data,axis=4)
           
        # ... so spatial axes can be averaged here without changing their
        # indices.
        if (avgaxes != ()):
            grid_data = np.mean(grid_data,axis=avgaxes) 

        #return avg_data
        return grid_data 

    def contour(self,axes,startrec=0,endrec=None):

        """
            NOT TO BE OVERRIDDEN UNLESS ABSOLUTELY NECESSARY
            Wrapper for averaged_data, returns easily plottable data 
        """

        avgaxes = [0,1,2]
        avgaxes.remove(axes[0])
        avgaxes.remove(axes[1])
        avgaxes = tuple(avgaxes)

        if (endrec==None): 
            endrec = self.maxrec

        data = self.averaged_data(startrec,endrec,avgaxes=avgaxes)
        X, Y = np.meshgrid(self.grid[axes[0]],self.grid[axes[1]],indexing='ij')
        return X, Y, data 
 
    def contour_timeseries(self,axes):
        pass
 
    def profile(self,axis,startrec=0,endrec=None):

        """
            NOT TO BE OVERRIDDEN UNLESS ABSOLUTELY NECESSARY
            Wrapper for averaged_data, returns easily plottable data 
        """

        avgaxes = [0,1,2]
        avgaxes.remove(axis)
        avgaxes = tuple(avgaxes)

        if (endrec==None): 
            endrec = self.maxrec

        data = self.averaged_data(startrec,endrec,avgaxes=avgaxes)
        return self.grid[axis], data
    
    def profile_timeseries(self,axis):
        pass


    def write_dx_file(self,startrec,endrec,writedir=None):

        """
           Write MD field to dx file format which is primarily
           useful for importing into VMD, see 
           http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dxplugin.html
           and format website http://www.opendx.org/index2.php
        """

        #Get field data
        for rec in range(startrec,endrec):
            data = self.read(startrec=rec,endrec=rec)

            Nx, Ny, Nz = [len(gridxyz) for gridxyz in self.grid]
            dx,dy,dz = [(self.grid[i][1] - self.grid[i][0]) for i in range(3)]
            Lx = float(Nx) * dx; Ly = float(Ny) * dy; Lz = float(Nz) * dz
            originx = -Lx/2.0
            originy = -Ly/2.0
            originz = -Lz/2.0

            #Get file name
            if (writedir == None):
                writedir = self.fdir + '/vmd/vol_data/'

            dxFileName = writedir + 'DATA' + str(rec) + '.dx'

            #Write data
            with open(dxFileName,'w+') as f:

                # - - Write Header - -
                f.write("object 1 class gridpositions counts%8.0f%8.0f%8.0f\n" % (Nx,Ny,Nz))
                f.write("origin%16g%16g%16g\n" % (originx,originy,originz))
                f.write("delta %16g 0 0\n" % dx)
                f.write("delta 0 %16g 0\n" % dy)
                f.write("delta 0 0 %16g\n" % dz)
                f.write("object 2 class gridconnections counts%8.0f%8.0f%8.0f\n" % (Nx,Ny,Nz))
                f.write("object 3 class array type double rank 0 items%8.0f follows\n" % (Nx*Ny*Nz))

                # - - Write Data - -
                col=1
                for i in range(Nx):
                    for j in range(Ny):
                        for k in range(Nz):
                            f.write("%16E" % (1.0*j + 2.0)) # data[i,j,k,0,0])
                            col=col+1
                            if (col>3):
                                f.write(' \n')
                                col=1

                # - - Write Footer - - 
                if (col != 1):
                    f.write('           \n')
                f.write('object "'+str(self).split('.')[1].split(' ')[0]+'" class field \n')

