#! /usr/bin/env python
import numpy as np
import scipy.ndimage

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

    def read(self,startrec,endrec,**kwargs):

        """
            TO BE OVERRIDDEN IN COMPLICATED FIELDS.
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
    
    def averaged_data(self,startrec,endrec,avgaxes=(),avgtime=True,**kwargs):

        """
            TO BE OVERRIDDEN IN COMPLICATED FIELDS.
            Average the data in the user-specified way.
        """
        
        # Read 4D time series from startrec to endrec
        grid_data = self.read(startrec,endrec,**kwargs)

        # Time axis is the final (fourth) axis...
        if (avgtime):
            grid_data = np.mean(grid_data,axis=4)
           
        # ... so spatial axes can be averaged here without changing their
        # indices.
        if (avgaxes != ()):
            grid_data = np.mean(grid_data,axis=avgaxes) 

        #return avg_data
        return grid_data 

    def contour(self,axes,startrec=0,endrec=None,**kwargs):

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

        data = self.averaged_data(startrec,endrec,avgaxes=avgaxes,**kwargs)
        # Need version 1.7.1 of numpy or higher
        X, Y = np.meshgrid(self.grid[axes[0]],self.grid[axes[1]],indexing='ij')
        return X, Y, data 
 
    def contour_timeseries(self,axes):
        pass
 
    def profile(self,axis,startrec=0,endrec=None,**kwargs):

        """
            NOT TO BE OVERRIDDEN UNLESS ABSOLUTELY NECESSARY
            Wrapper for averaged_data, returns easily plottable data 
        """

        avgaxes = [0,1,2]
        avgaxes.remove(axis)
        avgaxes = tuple(avgaxes)

        if (endrec==None): 
            endrec = self.maxrec

        data = self.averaged_data(startrec,endrec,avgaxes=avgaxes,**kwargs)
        return self.grid[axis], data
    
    def profile_timeseries(self,axis):
        pass


    def write_dx_file(self,startrec,endrec,writedir=None,component=0):

        """
           Write MD field to dx file format which is primarily
           useful for importing into VMD, see 
           http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dxplugin.html
           and format website http://www.opendx.org/index2.php
           NOTE -- VMD dx format assumes data points are located at the
                   cell vertices while Field class and it's children contain
                   cell centred data
        """

        #Get field data
        datamin = []; datamax = []
        for rec in range(startrec,endrec):

            data = self.read(startrec=rec,endrec=rec)

            #Return minimum and maximum values
            datamin.append(np.min(data[:,:,:,component,:]))
            datamax.append(np.max(data[:,:,:,component,:]))

            Nx, Ny, Nz = data.shape[0], data.shape[1], data.shape[2]
            dx,dy,dz = [(self.grid[i][1] - self.grid[i][0]) for i in range(3)]
            Lx = float(Nx) * dx; Ly = float(Ny) * dy; Lz = float(Nz) * dz
            originx = -Lx/2.0
            originy = -Ly/2.0
            originz = -Lz/2.0

            data = self.cellcentre2vertex(data[:,:,:,component,0])
            Nx_v, Ny_v, Nz_v = data.shape[0], data.shape[1], data.shape[2]

            #Get file name
            if (writedir == None):
                writedir = self.fdir + '/vmd/vol_data/'

            dxFileName = writedir + 'DATA' + str(rec) + '.dx'

            #Write data
            with open(dxFileName,'w+') as f:

                # - - Write Header - -
                #dx_ = dx*1.25; dy_ = dy*1.25; dz_ = dz*1.25
                f.write("object 1 class gridpositions counts%8.0f%8.0f%8.0f\n" % (Nx_v,Ny_v,Nz_v))
                f.write("origin%16g%16g%16g\n" % (originx,originy,originz))
                f.write("delta %16g 0 0\n" % dx)
                f.write("delta 0 %16g 0\n" % dy)
                f.write("delta 0 0 %16g\n" % dz)
                f.write("object 2 class gridconnections counts%8.0f%8.0f%8.0f\n" % (Nx_v,Ny_v,Nz_v))
                f.write("object 3 class array type double rank 0 items%8.0f follows\n" % (Nx_v*Ny_v*Nz_v))

                # - - Write Data - -
                col=1
                for i in range(Nx_v):
                    for j in range(Ny_v):
                        for k in range(Nz_v):
                            f.write("%16E" %  data[i,j,k])
                            col=col+1
                            if (col>3):
                                f.write(' \n')
                                col=1

                # - - Write Footer - - 
                if (col != 1):
                    f.write('           \n')
                f.write('object "'+str(self).split('.')[1].split(' ')[0]+'" class field \n')

        return np.mean(datamin), np.mean(datamax)


    def cellcentre2vertex(self,celldata):

        """
           Routine to return grid data on an array one larger than the existing
           cell centred data - currently uses zoom for simplicity

        """
        Nx, Ny, Nz = celldata.shape[0], celldata.shape[1], celldata.shape[2]
        vertexdata = scipy.ndimage.zoom(celldata,((Nx+1)/float(Nx),
                                                  (Ny+1)/float(Ny),
                                                  (Nz+1)/float(Nz)))
        return vertexdata

        # Innards
#        for i in range(1,Nx):
#            for j in range(1,Ny):
#                for k in range(1,Nz):
#                    vertexdata[i,j,k] = np.mean(celldata[i-1:i+2,j-1:j+2,k-1:k+2]) 

