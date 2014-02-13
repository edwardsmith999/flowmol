#! /usr/bin/env python
import numpy as np
import copy

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
        self.grid = self.Raw.grid
        self.maxrec = self.Raw.maxrec

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

