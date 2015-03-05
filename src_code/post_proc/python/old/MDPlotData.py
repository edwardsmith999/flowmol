from MDFields import *

class MD_PlotData():

    """
        MD_PlotData: Return the arguments required by matplotlib's plotting
                     functions for a specific figure type.  
        Author: David Trevelyan, April 2013 
        Each class method should return arguments required by a single 
        matplotlib plotting function for a single figure. Each time you need 
        a new figure, please save a minimal method in this class so that 
        similar figures can be (hopefully!) trivially plotted in future.

        IMPORTANT: It should be possible to get all the information from the
        classes in MDFields.py. If the functionality you need is not there, 
        please add it to those classes rather than changing the MD_RawData 
        class or cluttering this MD_PlotData class. 

        This class will likely get quite long, so it may be sensible to 
        break it up into categories. It's not necessary just yet, though,
        as of April 2013.

    """

    def __init__(self,fdir,cpol_bins=False):
        self.fdir = fdir
        self.cpol_bins = cpol_bins

    def get_vplane_quiver_args(self,plane,haxis,vaxis,minrec,maxrec):

        # Instantiate velocity data object
        vData = VBins(self.fdir,cpol_bins=self.cpol_bins)
        # Extract 3D velocity field averaged over 1D of bins
        vplane, binspaces = vData.get_field(minrec,maxrec,sumaxes=(plane))
        # Get bin center positions on both axes for every field point
        X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')
        # Extract components of velocity for each bin
        U = vplane[:,:,haxis]
        V = vplane[:,:,vaxis]
        return X, Y, U, V

    def get_vplane_streamplot_args(self,*args,**kw):
    
        # Same arguments as quiver, just transposed    
        X, Y, U, V = self.get_vplane_quiver_args(*args,**kw)
        X = X.transpose()
        Y = Y.transpose()
        U = U.transpose()
        V = V.transpose()
        return X, Y, U, V

    def get_vslice_plot_args(self,axis,minrec,maxrec,binlimits=None):

        # Get which axes to average over
        sumaxes = []    
        for ax in range(3):
            if (ax != axis): sumaxes.append(ax)
        sumaxes = tuple(sumaxes)

        # Instantiate velocity data object
        vData = VBins(self.fdir,cpol_bins=self.cpol_bins)
        # Get 3D velocity field in 3D bins
        vslice, binspaces = vData.get_field(minrec,maxrec,sumaxes=sumaxes,
                                            binlimits=binlimits)
    
        return binspaces[axis], vslice

    def get_vplane_splot_args(self,plane,haxis,vaxis,minrec,maxrec,
                              binlimits=None):

        # Instantiate velocity data object
        vData = VBins(self.fdir,cpol_bins=self.cpol_bins)
        # Extract 3D velocity field averaged over 1D of bins
        vplane, binspaces = vData.get_field(minrec,maxrec,sumaxes=(plane),
                                            binlimits=binlimits)
        # Get bin center positions on both axes for every field point
        X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')
        # Return all components of velocity for each bin
        U = vplane

        return X, Y, U      

    def get_mslice_plot_args(self,axis,minrec,maxrec):

        # Get which axes to average over
        avgaxes = []    
        for ax in range(3):
            if (ax != axis): avgaxes.append(ax)
        avgaxes = tuple(avgaxes)

        mData = MassBins(self.fdir,cpol_bins=self.cpol_bins)
        mslice, binspaces = mData.get_field(minrec,maxrec,sumaxes=(avgaxes))
        #Take zeroth component to 'Squeeze' array to lower dimensionality
        mslice = mslice[:,0] 

        return binspaces[axis], mslice

    def get_mplane_splot_args(self,plane,haxis,vaxis,minrec,maxrec):

        # Instantiate mass data object
        mData = MassBins(self.fdir,cpol_bins=self.cpol_bins)
        # Extract 3D velocity field averaged over 1D of bins
        mplane, binspaces = mData.get_field(minrec,maxrec,sumaxes=(plane))
        # Get bin center positions on both axes for every field point
        X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')
        #Take zeroth component to 'Squeeze' array to lower dimensionality
        m = mplane[:,:,0] 

        return X, Y, m        

    def get_density_prof_args(self,axis,minrec,maxrec,binlimits=None):

        # Get which axes to average over
        avgaxes = []    
        for ax in range(3):
            if (ax != axis): avgaxes.append(ax)
        avgaxes = tuple(avgaxes)
    
        ddata = DensityBins(self.fdir,cpol_bins=self.cpol_bins)
        density, binspaces = ddata.get_field(minrec,maxrec,meanaxes=(avgaxes),
                                             binlimits=binlimits)
        density = density[:,0]

        return binspaces[axis], density

    def get_density_splot_args(self,plane,haxis,vaxis,minrec,maxrec,
                               binlimits=None):
    
        ddata = DensityBins(self.fdir,cpol_bins=self.cpol_bins)
        density, binspaces = ddata.get_field(minrec,maxrec,meanaxes=(plane),
                                             binlimits=binlimits)
        X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')
        density = density[:,:,0]

        return X, Y, density

    def get_pVA_prof_args(self,filename,axis,minrec,maxrec,peculiar=False):
        
        # Get which axes to average over
        avgaxes = []    
        for ax in range(3):
            if (ax != axis): avgaxes.append(ax)
        avgaxes = tuple(avgaxes)
    
        # Instantiate pVA Field object    
        pVA_obj = pVABins(self.fdir,filename,cpol_bins=self.cpol_bins)

        # Extract 3D pressure tensor field averaged over 1D of bins
        Pslices, binspaces = pVA_obj.get_field(minrec,maxrec,
                             meanaxes=(avgaxes),peculiar=peculiar)
    
        return binspaces[axis], Pslices

    def get_pVA_splot_args(self,filename,plane,haxis,vaxis,minrec,maxrec,
                           peculiar=False):
        
        # Instantiate VA data object
        pVA_obj = pVABins(self.fdir,filename,cpol_bins=self.cpol_bins)

        # Extract 3D pressure tensor field 
        Pplane, binspaces = pVA_obj.get_field(minrec,maxrec,meanaxes=(plane),
                                              peculiar=peculiar)
        X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')
        return X, Y, Pplane

    def get_CV_splot_args(self,filename,plane,haxis,vaxis,minrec,maxrec):

        # Instantiate CV data object
        CV_obj = CV_data(self.fdir,filename,cpol_bins=self.cpol_bins)

        # Extract 3D velocity field averaged over 1D of bins
        Pplane, binspaces = CV_obj.get_field(minrec,maxrec,meanaxes=(plane))
        X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')
        return X, Y, Pplane

    def get_Tplane_splot_args(self,plane,haxis,vaxis,minrec,maxrec,
                              peculiar=True):

        # Instantiate temperature data object
        TData = TBins(self.fdir,cpol_bins=self.cpol_bins)

        # Extract 3D velocity field averaged over 1D of bins
        Tplane, binspaces = TData.get_field(minrec,maxrec,sumaxes=(plane),
                                            peculiar=peculiar)

        # Get bin center positions on both axes for every field point
        X, Y = np.meshgrid(binspaces[haxis],binspaces[vaxis],indexing='ij')

        return X, Y, Tplane

    def get_Tslice_plot_args(self,axis,minrec,maxrec,peculiar=True):

        # Get which axes to average over
        avgaxes = []    
        for ax in range(3):
            if (ax != axis): avgaxes.append(ax)
        avgaxes = tuple(avgaxes)

        TData = TBins(self.fdir,cpol_bins=self.cpol_bins)
        Tslice, binspaces = TData.get_field(minrec,maxrec,sumaxes=(avgaxes),
                                            peculiar=peculiar)

        return binspaces[axis], Tslice


    def get_vfield_energy_spectra(self, meanaxes, component, minrec, maxrec,
                                  tavg_rec, binlimits=None, fftaxis=None,
                                  ffttime=False, savefile=None, 
                                  verify_Parseval=False):

        """
            meanaxes - tuple (pre-fft avg axes)
            component - integer (velocity component/direction)
            minrec - integer
            maxrec - integer
            tavg_rec - integer (running average width in records)
            binlimits - [(min,max),(),()] (only pull data from certain bins
            fftaxis - integer
            ffttime - bool
            savefile - string
            verify_Parseval - bool

        """

        class vfield_spectra:

            def __init__(self,DataObj):
                # Initialise velocity time series to nothing
                self.v_ts = []
                self.vDataObj = DataObj
                # Keep track of which axes have been averaged, summed, etc.
                # This is kept so the user can just specify which axis they
                # want to fft, without having to worry about how the dimensions
                # of the array have changed on averaging/summing etc.
                # THE AXES CORRESPOND TO: [X,Y,Z,TIME,COMPONENT]
                self.axesreduced = [False]*5
                # Local copy of which axes represent real space (each entry is
                # "popped" out when reduced (see self.axesreduced)
                self.realspacepopped = [True]*5
                # Energy normalisation factor, initialise to one
                self.N = 1.0
                # Initialise windowed to be False
                self.windowed = False

            def markaxisfourierspace(self,axis):
                newaxis = int(axis - np.sum(self.axesreduced[:axis]))
                self.realspacepopped[newaxis] = False

            def markaxisreduced(self,axis):
                newaxis = int(axis - np.sum(self.axesreduced[:axis]))
                self.axesreduced[axis] = True
                self.realspacepopped.pop(newaxis)

            def populate_timeseries(self,minrec,maxrec,drec,meanaxes):

                # Override meanaxes to be a tuple if possible
                if (not isinstance(meanaxes,tuple)):
                    try:
                        meanaxes = tuple([meanaxes,]) 
                    except:
                        print('Failed to make meanaxes a tuple')

                # Read time series of velocity field records, averaging over
                # length "drec" before storing in array.
                # maxrec+1 in range to catch actual maxrec, drec-1 in call to 
                # get the right number of records (consider drec=1) 

                # Loop over records and append to time series 
                for rec in range(minrec,maxrec+1,drec):
                    v, binspaces = self.vDataObj.get_field(rec,rec+drec-1,
                                                 binlimits=binlimits,
                                                 sumaxes=(meanaxes))
                    self.v_ts.append(v)

                # Turn list into array
                self.v_ts = np.array(self.v_ts)

                # Put time axis in penultimate position, component in last, 
                # remaining spatial coordinates at start, so e.g. for single 
                # meanaxis we get u(axis1,axis2,time,component)

                if (len(meanaxes) == 0):
                    self.v_ts = np.transpose(self.v_ts,axes=(1,2,3,0,4))
                elif (len(meanaxes) == 1):
                    self.v_ts = np.transpose(self.v_ts,axes=(1,2,0,3))
                elif (len(meanaxes) == 2):
                    self.v_ts = np.transpose(self.v_ts,axes=(1,0,2))
                elif (len(meanaxes) == 3):
                    # Array will already be in the correct order
                    pass

                for axis in meanaxes:
                    self.markaxisreduced(axis)

            def extract_component(self,component):

                # Get index array for all elements
                slicer = [np.arange(i) for i in self.v_ts.shape]
                # Set only desired component (final position)
                slicer[-1] = np.array([component])
                slicer = np.ix_(*slicer)
                
                # Component extraction
                self.v_ts = self.v_ts[slicer]

                # Collapse final dimension (which is only of length 1 anyway)
                self.v_ts = np.reshape(self.v_ts,self.v_ts.shape[:-1])
                self.markaxisreduced(4)

            def apply_time_window(self):

                def window_axis_function(a):
                    a = a * self.window
                    return a

                timeaxis = 3
                newaxis = timeaxis - np.sum(self.axesreduced[:timeaxis])

                N = self.v_ts.shape[newaxis]
                #self.window = np.hanning(N)
                self.window = np.ones(N)
                self.window[-1] = 0.0
                self.window[0] = 0.0

                # Save "window summed and squared" (see Numerical Recipes)
                self.wss = np.sum(self.window**2.0)/N

                # Apply window
                self.v_ts = np.apply_along_axis(window_axis_function,newaxis,self.v_ts)        
                self.windowed = True

            def fft(self,fftaxis,window=False):

                if (self.axesreduced[fftaxis]):
                    print('You\'re trying to FFT over an axis that is ' +
                          'already reduced')
                    quit()

                newaxis = fftaxis - np.sum(self.axesreduced[:fftaxis])

                # Perform FFT
                self.v_ts = np.fft.fft(self.v_ts,axis=newaxis)
                self.N = self.N * self.v_ts.shape[newaxis]

                # Mark axis as no longer representing real space
                self.markaxisfourierspace(fftaxis) 

            def set_energyfield(self):

                # Work out which axes we can average over (i.e. any that we
                # haven't Fourier transformed)
                axes = np.where(np.array(self.realspacepopped)==True)[0]

                # Energy field
                self.E = np.abs(self.v_ts)**2.0

                if (self.windowed == True):
                    self.E = self.E / self.wss 

                # Average remaining dimensions
                self.E = np.mean(self.E,axis=tuple(axes))
                #self.E = np.sum(self.E,axis=tuple(axes))

                # Delete duplicate parts of the energy spectrum and double
                # energy contributions of middle wavenumbers
                cutout = [np.s_[0:i/2+1:1] for i in self.E.shape] 
                double = [np.s_[1:i/2  :1] for i in self.E.shape] 
                self.E = self.E[cutout]
                self.E[double] = self.E[double] * 2.0

        # Create field reading object
        self.vDataObj = VBins(self.fdir,cpol_bins=self.cpol_bins)    
        # Create spectra calculating object
        VField = vfield_spectra(self.vDataObj)
        # Read the time series of the velocity field desired
        VField.populate_timeseries(minrec,maxrec,tavg_rec,meanaxes)
        # Extract the cartesian component we are interested in
        VField.extract_component(component)
        # Apply window if we're going to fft over time axis
        if (ffttime):
            VField.apply_time_window()

        # If we want to verify P's thm, store real space energy
        if (verify_Parseval):
            Esumreal = np.sum(np.abs(VField.v_ts)**2.0)

        # FFT in space or time depending on input to function, and 
        # calculate power spectra.
        if (fftaxis):
            # Override fftaxis to be a tuple if possible
            if (not isinstance(fftaxis,tuple)):
                try:
                    fftaxis = tuple([fftaxis,]) 
                except:
                    print('Failed to make fftaxis a tuple')
            for axis in fftaxis:
                VField.fft(axis)

        if (ffttime):
            timeaxis = 3
            VField.fft(timeaxis)

        # If we want to verify P's thm, store Fourier space energy
        # (before we undo the window)
        if (verify_Parseval):
            Esumfft = np.sum(np.abs(VField.v_ts)**2.0)/VField.N
            ratio = abs(Esumreal - Esumfft)/Esumreal 
            perc = (1. - ratio)*100.
            print('Parseval thm (discounting window): ' + "%5.2f"%perc + '%')

        # Set the energy field
        VField.set_energyfield()

        if (savefile):
           with open(savefile,'w') as f:
                f.write(VField.E/VField.N)
    
        # Return energy field (don't normalise yet)
        return VField.E/VField.N
