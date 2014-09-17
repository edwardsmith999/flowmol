#! /usr/bin/env python
import numpy as np
import scipy.ndimage
import scipy.interpolate as interp
import scipy.ndimage.interpolation as interp2
import matplotlib.pyplot as plt

class OutsideRecRange(Exception):
    pass

class Field():

    """
        Abstract base class to be inherited by MDField, CFDField and CPLField.

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

        In line with RawData readers, axes correspond to:

            [ spatial ax1, spatial ax2, spatial ax3, record, component ] 
            [     0      ,      1     ,      2     ,    3  , 4 (&-1)   ] 

    """
    
    def __init__(self, RawDataObj):
        self.Raw = RawDataObj
        self.fdir = self.Raw.fdir
        self.grid = self.Raw.grid
        self.maxrec = self.Raw.maxrec

    def read(self,startrec,endrec,**kwargs):

        """
            TO BE OVERRIDDEN IN COMPLICATED FIELDS.
            Method that returns grid data that is read by Raw.
            
        """
        if (endrec > self.maxrec):
            print('Record ' + str(endrec) + ' is greater than the maximum '
                  'available (' + str(self.maxrec) + ').')
            raise OutsideRecRange

        grid_data = self.Raw.read(startrec,endrec,**kwargs)
        return grid_data
    
    def averaged_data(self,startrec,endrec,avgaxes=(),**kwargs):

        """
            TO BE OVERRIDDEN IN COMPLICATED FIELDS.
            Average the data in the user-specified way.
        """
        
        # Read 4D time series from startrec to endrec
        grid_data = self.read(startrec,endrec,**kwargs)
           
        # Average over axes
        if (avgaxes != ()):
            grid_data = np.mean(grid_data,axis=avgaxes) 

        #return avg_data
        return grid_data 

    def contour(self,axes,startrec=0,endrec=None,**kwargs):

        """
            NOT TO BE OVERRIDDEN UNLESS ABSOLUTELY NECESSARY
            Wrapper for averaged_data, returns easily plottable data 
        """

        avgaxes = [0,1,2,3]
        avgaxes.remove(axes[0])
        avgaxes.remove(axes[1])
        avgaxes = tuple(avgaxes)

        if (endrec==None): 
            endrec = self.maxrec

        data = self.averaged_data(startrec,endrec,avgaxes=avgaxes,**kwargs)
        # Need version 1.7.1 of numpy or higher
        X, Y = np.meshgrid(self.grid[axes[0]],self.grid[axes[1]],indexing='ij')
        return X, Y, data 
 
    def profile(self,axis,startrec=0,endrec=None,**kwargs):

        """
            NOT TO BE OVERRIDDEN UNLESS ABSOLUTELY NECESSARY
            Wrapper for averaged_data, returns easily plottable data 
        """

        avgaxes = [0,1,2,3]
        avgaxes.remove(axis)
        avgaxes = tuple(avgaxes)

        if (endrec==None): 
            endrec = self.maxrec

        data = self.averaged_data(startrec,endrec,avgaxes=avgaxes,**kwargs)
        return self.grid[axis], data

    def power_spectrum(self,startrec=0,endrec=None,preavgaxes=(),fftaxes=(),
                       postavgaxes=(), windowaxis=None, verify_Parseval=True,
                       savefile=None,**kwargs):

        """
            Calculates power spectrum
        """

        class AxisManager():

            def __init__(self):
                self.nreduced = 0
                self.axisactive = [True]*5

            def reduce_axes(self,axes):
                for axis in axes:
                    self.nreduced += 1
                    self.axisactive[axis] = False 

            def current_axis_number(self,axis_request):
                if (not self.axisactive[axis_request]):
                    return None
                n = 0
                for otheraxis in range(axis_request):
                    if (not self.axisactive[otheraxis]):
                        n += 1
                return int(axis_request - n)
            
            def current_axes_numbers(self,axes_request):
                newaxes = [] 
                for axis in list(axes_request):
                    newaxes.append(self.current_axis_number(axis))
                return tuple(newaxes)

        def managed_mean(axismanager, data, avgaxes):
            newaxes = axismanager.current_axes_numbers(avgaxes)
            if (None in newaxes):
                quit("Can't average over an axis that has been reduced")
            avgdata = np.mean(data, axis=newaxes) 
            axismanager.reduce_axes(avgaxes)
            return avgdata

        def managed_fft(axismanager, data, fftaxes):
            newaxes = axismanager.current_axes_numbers(fftaxes)
            if (None in newaxes):
                quit("Can't fft over an axis that has been reduced")
            fftdata = np.fft.fftn(data,axes=newaxes)
            return fftdata

        def managed_energyfield(axismanager, data, fftaxes):
            
            def add_negatives(a):
                b = a
                for i in range(1,len(b)/2):
                    b[i] += a[-i]
                return b

            newfftaxes = axismanager.current_axes_numbers(fftaxes)

            # Count N
            N = 1
            for axis in newfftaxes:
                N = N * data.shape[axis]

            # Initial energy in every wavenumber (pos and neg)
            E = np.abs(data)**2.0
             
            #Add negative contributions to positive wavenumbers
            nd = len(E.shape)
            slicer = [slice(None)]*nd
            for axis in newfftaxes:

                E = np.apply_along_axis(add_negatives,axis,E)
                # Discard negative parts after adding to positive
                k = E.shape[axis]
                mid = int(np.ceil(float(k)/2.0))
                cutout = np.s_[0:mid+1:1]
                slicer[axis] = cutout
                E = E[slicer]

                # Refresh slice for new axis calculation
                slicer[axis] = slice(None)

            return E/N

        def managed_window(axismanager, data, windowaxis):

            def window_axis_function(a, window):
                a = a * window
                return a

            newaxis = axismanager.current_axis_number(windowaxis)

            N = data.shape[newaxis]
            window = np.hanning(N)

            # Save "window summed and squared" (see Numerical Recipes)
            wss = np.sum(window**2.0)/float(N)

            # Apply window
            windoweddata = np.apply_along_axis(window_axis_function, 
                                               newaxis, data, window)

            return windoweddata, wss 

        # ---------------------------------------------------------------- 
        # Checks
        if (not isinstance(preavgaxes,tuple)):
            try:
                preavgaxes = tuple([preavgaxes])
            except:
                print('Failed to make preavgaxes and fftaxes a tuple')
        if (not isinstance(fftaxes,tuple)):
            try:
                fftaxes = tuple([fftaxes])
            except:
                print('Failed to make preavgaxes and fftaxes a tuple')
        if (not isinstance(postavgaxes,tuple)):
            try:
                postavgaxes = tuple([postavgaxes])
            except:
                print('Failed to make preavgaxes and fftaxes a tuple')

        if (4 in postavgaxes or 4 in fftaxes or 4 in preavgaxes):
            message = "WARNING: you're asking me to average or fft over "
            message += "each component of the field. I don't know how to "
            message += "deal with this right now. Aborting."
            quit(message)

        if (windowaxis):
            if (windowaxis in preavgaxes or windowaxis in postavgaxes):
                message = "Warning: you're asking me to window over an axis "
                message += "that you're eventually going to average over. "
                message += "Aborting."
                quit(message)
            if (windowaxis not in fftaxes):
                message = "Warning: you're asking me to window over an axis "
                message += "that won't be Fourier transformed. This makes no "
                message += "sense. Aborting."
        
        # ---------------------------------------------------------------- 
        # Do the process 

        if (endrec==None):
            endrec = self.maxrec
         
        axisman = AxisManager() 
        data = self.read(startrec, endrec,**kwargs)
        data = managed_mean(axisman, data, preavgaxes)

        if (windowaxis):
            data, wss = managed_window(axisman, data, windowaxis)

        if (verify_Parseval):
            Esumreal = np.sum(np.abs(data)**2.0)

        fftdata = managed_fft(axisman, data, fftaxes)
        energy = managed_energyfield(axisman, fftdata, fftaxes)

        if (verify_Parseval):
            Esumfft = np.sum(energy)
            ratio = abs(Esumreal - Esumfft)/Esumreal 
            perc = (1. - ratio)*100.
            print('Parseval thm (discounting window): ' + "%9.6f"%perc + '%')

        if (windowaxis):
            energy = energy / wss

        energy = managed_mean(axisman, energy, postavgaxes)

        if (savefile):
            with open(savefile,'w') as f:
                f.write(energy)
        
        return energy


    def grad(self,vdata, dx=None, 
                         dy=None, 
                         dz=None):

        """
            Return the gradient of a vector field
        """

        #if (dxyz is None):
        #    dxyz=[self.Raw.dx,self.Raw.dy,self.Raw.dz]

        if (dx is None):
            dx=self.Raw.dx
        if (dy is None):
            dy=self.Raw.dy
        if (dz is None):
            dz=self.Raw.dz

        grad_temp = np.empty((vdata.shape[0],
                              vdata.shape[1],
                              vdata.shape[2],
                              vdata.shape[4]))
        ndims = vdata.shape[4]
        gradv = np.empty(list(vdata.shape[:-1]) + [3*ndims])
        for rec in range(gradv.shape[-2]):
            for ixyz in range(ndims):
                grad_temp = np.gradient(vdata[:,:,:,rec,ixyz], 
                                        dx, dy, dz)

                for jxyz in range(3):
                    c = 3*ixyz + jxyz
                    gradv[:,:,:,rec,c] = grad_temp[jxyz]

        return gradv

    def write_dx_file(self,startrec,endrec,writedir=None,component=0,**kwargs):

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

            data = self.read(startrec=rec,endrec=rec,**kwargs)

            #Return minimum and maximum values
            datamin.append(np.min(data[:,:,:,:,component]))
            datamax.append(np.max(data[:,:,:,:,component]))

            Nx, Ny, Nz = data.shape[0], data.shape[1], data.shape[2]
            dx,dy,dz = [(self.grid[i][1] - self.grid[i][0]) for i in range(3)]
            Lx = float(Nx) * dx; Ly = float(Ny) * dy; Lz = float(Nz) * dz
            originx = -Lx/2.0
            originy = -Ly/2.0
            originz = -Lz/2.0

            data = self.cellcentre2vertex(data[:,:,:,0,component])
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

    # Write ascii type field
    def write_asciifield(self,startrec,endrec,
                        writedir=None,maptocosine=True,
                        flipdir=[False,True,False],**kwargs):

        #Get file name
        if (writedir == None):
            writedir = self.fdir

        Nxyz = [Nx,Ny,Nz]
        for n,flip in enumerate(flipdir):
            if flip:
                mindir[n] = Nxyz[n]
                maxdir[n] = 0
            else:
                mindir[n] = 0
                maxdir[n] = Nxyz[n] 

        data = self.read(startrec=startrec,endrec=endrec,**kwargs)

        outfiles =[]
        for rec in range(startrec,endrec):

            FileName = writedir + 'u' + str(rec) + '.asc'
            outfiles.append(FileName)
            with open(FileName,'w+') as f:
                for nx in range(mindir[0],maxdir[0]):
                    for ny in range(mindir[1],maxdir[1]):
                        for nz in range(mindir[2],maxdir[2]):
                            for ixyz in data.shape[4]:
                                if maptocosine:
                                    f.write(str(
                                                self.map_data_lineartocosine(
                                                data[nx,ny,nz,rec,ixyz],
                                                self.Raw.Ny,
                                                self.Raw.a,
                                                self.Raw.b))
                                            +"\n")
                                else:
                                    f.write(str(
                                                data[nx,ny,nz,rec,ixyz])
                                            +"\n")

        return outfiles

    def map_data_lineartocosine(self,values_on_linear_grid,Ny,a,b):
        """
            Map data on a linear grid to a cosine grid 
        """
        ycells = np.linspace(0, Ny, Ny)
        ylin = np.linspace(a, b, Ny)
        ycos = 0.5*(b+a) - 0.5*(b-a)*np.cos((ycells*np.pi)/(Ny-1))
        plt.plot(ylin,values_on_linear_grid,'o-',alpha=0.4,label='lineartocosine Before')
        values_on_cosine_grid = interp.griddata(ylin, values_on_linear_grid, 
                                                ycos, method='cubic',
                                                fill_value=values_on_linear_grid[-1])
        plt.plot(ycos,values_on_cosine_grid,'x-',label='lineartocosine After')
        plt.legend()
        plt.show()
        return values_on_cosine_grid

    def map_data_cosinetolinear(self,values_on_cosine_grid,Ny,a,b):
            """
                Map data on a cosine grid to a linear grid 
            """
            ycells = np.linspace(0, Ny, Ny)
            ylin = np.linspace(a, b, Ny)
            ycos = 0.5*(b+a) - 0.5*(b-a)*np.cos((ycells*np.pi)/(Ny-1))
            plt.plot(ycos,values_on_cosine_grid,'x-',label='cosinetolinear Before')
            values_on_linear_grid = interp.griddata(ycos, values_on_cosine_grid, 
                                                    ylin, method='cubic',
                                                    fill_value=values_on_cosine_grid[-1])
            #values_on_linear_grid = interp2.map_coordinates(values_on_cosine_grid,ycos,output=ylin)
            plt.plot(ylin,values_on_linear_grid,'o-',alpha=0.4,label='cosinetolinear After')
            plt.legend()
            plt.show()
            return values_on_linear_grid


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

