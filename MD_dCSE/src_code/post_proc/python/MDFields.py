
import os
import glob
import numpy as np

from MDRawData import MD_RawData
from HeaderData import HeaderData

# Abstract Field Class
class Field():

    """
        Field:   An abstract data averaging base class to be inherited by 
                 derived field classes.
        Authors: David Trevelyan & Ed Smith, April 2013

        Field contains the information necessary to read a binary data file
        using the MD_RawData class, and average the information spatially
        and temporally.

        The following attributes are required by MD_RawData and should be 
        specified by the derived classes:
            
            fname   - file name of binary data file of interest, 
            dtype   - binary data type string, 'i'/'d' for int/float, and
            nperbin - number of values per bin. 
        
        Functionality:

            minrec, maxrec - specify a minimum and maximum record in time
                             over which to average the data.

            sumaxes        - specify axes over which to sum the data in the
                             4D array.
            
            meanaxes       - specify axes over which to average the data in 
                             the 4D array.

        EXAMPLE USE OF GET_BINS:
                
        Read data and average over records 100 to 199:
    
            get_bins(100,200) 
                
        Read data, average over recs 100 to 199, and sum values in x-z plane
        to get the TOTAL in a y-slice:
    
            get_bins(100,200,sumaxes=(0,2)) 

        Read data, average over recs 100 to 199, and average values over the 
        z-direction to get an x/y field:

            get_bins(100,200,meanaxes=(2)) 

    """

    def __init__(self,fdir,cpol_bins=False): 

        self.fdir = fdir
        self.cpol_bins = cpol_bins

        # Create raw data reading object and get the topology of the bins
        self.Raw = MD_RawData(self.fdir,self.fname,self.dtype,self.nperbin,
                         self.cpol_bins)
        self.nbins, self.binspaces = self.Raw.get_bintopology()

        # Store maximum possible record index
        self.maxrecindex = self.get_maxrec() 
       
    def get_maxrec(self):

        if (glob.glob(self.fdir+self.fname)):

            filesize = os.path.getsize(self.fdir+self.fname)
            if (self.dtype == 'i'):
                maxrec = filesize/(4*self.nperbin*np.prod(self.nbins)) - 1
            elif (self.dtype == 'd'):
                maxrec = filesize/(8*self.nperbin*np.prod(self.nbins)) - 1
            else:
                quit('Unrecognised dtype in Field.get_nrecs')

        elif (glob.glob(self.fdir+self.fname+'.*')):

            filelist = glob.glob(self.fdir+self.fname+'.*')
            sortedlist = sorted(filelist)
            maxrec = int(sortedlist[-1].split('.')[-1])
            
        else:
            print('Neither ' + fname + ' nor ' + fname + '.* exist.')
            quit()

        return maxrec 
         

    def get_bins(self, minrec, maxrec, sumaxes=(), meanaxes=(), meantime=True,
                 sumtime=False, binlimits=None):
        
        """
            Read data file using MD_RawData class, average over time record
            range specified by minrec to maxrec. Sum over sumaxes (tuple), or
            average in directions specified by meanaxes (tuple).
            CURRENTLY CANNOT DO BOTH SUM AND MEAN BECAUSE THE SHAPE OF THE 
            ARRAY CHANGES AFTER YOU DO ONE OR THE OTHER.

        """

        if ( sumaxes != () ) and ( meanaxes != () ):
            print("""Currently cannot specify both sumaxes and meanaxes
                   because the shape of the array of interest changes once
                   it is summed/averaged for the first time. It is possible
                   to determine how meanaxes should be changed based on the
                   specification of sumaxes, but this just hasn\'t been
                   done yet.""")
            quit()

        if (maxrec > self.maxrecindex):
            quit('Record ' + str(maxrec) + ' is greater than the maximum '
                 'available (' + str(self.maxrecindex) + '). Aborting.')

        # Get bin data from file and get the mean/sum over time records (axis 4)
        binspaces = self.binspaces
        bins = self.Raw.get_bindata(minrec,nrecs=maxrec-minrec+1)

        # If bin limits are specified, return only those within range
        if (binlimits):
            # Initialise slice list as every index in bins
            s = [np.arange(i) for i in bins.shape]
            # Loop over axes and change slicer limits (and binspaces)
            for axis in np.arange(3):
                if (binlimits[axis]):
                    # Prepare slice array (+1 for python slicing convention)
                    s[axis] = np.arange(binlimits[axis][0],
                                        binlimits[axis][1]+1) 
                    # Chop out unwanted binspace entries (+1 for python slice)
                    binspaces[axis] = binspaces[axis][binlimits[axis][0]:
                                                      binlimits[axis][1]+1] 
            # Convert slice list to proper shape for numpy fancy indexing
            slicer = np.ix_(*s) 
            # Delete entries not in slicer
            bins = bins[slicer]

        if (meantime==True):
            bins = bins.mean(axis=4) 
        elif (sumtime==True):
            bins = bins.sum(axis=4) 

        # Sum or mean as appropriate
        if (sumaxes != ()):
            bins = bins.sum(axis=sumaxes)
        elif (meanaxes != ()):
            bins = bins.mean(axis=meanaxes)
       
        return bins, binspaces

    def get_binvolumes(self,binlimits=None):

        binspaces = self.binspaces
    
        if (self.cpol_bins == True):    

            r_oi = float(self.Raw.header.r_oi)
            r, theta, z = np.meshgrid((binspaces[0]+r_oi),
                                       binspaces[1],
                                       binspaces[2],
                                       indexing='ij')


            dr     = binspaces[0][1] - binspaces[0][0]
            dtheta = binspaces[1][1] - binspaces[1][0]
            dz     = binspaces[2][1] - binspaces[2][0]

            r = r + 0.5*dr
            binvolumes = r*dr*dtheta*dz

        else:

            x, y, z = np.meshgrid(binspaces[0],binspaces[1],binspaces[2],
                                  indexing='ij')


            dx = binspaces[0][1] - binspaces[0][0]
            dy = binspaces[1][1] - binspaces[1][0]
            dz = binspaces[2][1] - binspaces[2][0]

            binvolumes = np.ones(x.shape)*dx*dy*dz

        # If bin limits are specified, return only those within range
        if (binlimits):

            # Initialise slice list as every index in bins
            s = [np.arange(i) for i in binvolumes.shape]
            # Loop over axes and change slicer limits
            for axis in np.arange(len(binvolumes.shape)):
                if (binlimits[axis]):
                    # (+1 for python slicing convention) 
                    s[axis] = np.arange(binlimits[axis][0],
                                        binlimits[axis][1]+1) 
            # Convert slice list to proper shape for numpy fancy indexing
            slicer = np.ix_(*s) 
            # Delete entries not in slicer
            binvolumes = binvolumes[slicer]

        # Ensure binvolumes is the right shape for subsequent
        # broadcasting with other fields
        binvolumes = np.expand_dims(binvolumes,-1)
        return binvolumes

# Mass field    
class MassBins(Field):

    fname = 'mbins'
    dtype = 'i'
    nperbin = 1
    get_field = Field.get_bins

# Momentum field    
class MomBins(Field):

    fname = 'vbins'
    dtype = 'd'
    nperbin = 3
    get_field = Field.get_bins

# Pressure fields
class PBins(Field):

    dtype = 'd'
    nperbin = 9
    get_field = Field.get_bins

    def __init__(self,fdir,fname='pVA',cpol_bins=False):

        """
            PBins requires the specification of a filename by the
            user, allowing any of pVA, virial, or separate kinetic
            and configurational parts to be plotted with the same
            field class functionality.

        """

        Field.__init__(self,fdir,cpol_bins=cpol_bins)
        self.fname = fname  


# CV fields
class CV(Field):

    dtype = 'd'
    nperbin = 18
    get_field = Field.get_bins

    def __init__(self,fdir,fname,cpol_bins=False):

        """
            CV requires the specification of a filename by the
            user, allowing any of vfluxes, psurface, vsnap or F_ext
            to be plotted correctly

        """

        Field.__init__(self,fdir,cpol_bins=cpol_bins)
        self.fname = fname    


# Kinetic energy field
class KEBins(Field):
    
    fname = 'Tbins'
    dtype = 'd'
    nperbin = 1
    get_field = Field.get_bins

# Velocity field
class VBins():

    def __init__(self,fdir,cpol_bins=False):
        self.mdata = MassBins(fdir,cpol_bins)
        self.pdata = MomBins(fdir,cpol_bins)

    def get_field(self, minrec, maxrec, sumaxes=(), sumtime=True,
                  binlimits=None):

        """
            Get the velocity field from files vbins/mbins averaged over
            time records minrec->maxrec, AND averaged over 
            spatial directions specified by sumaxes.
            
            sumaxes   - *tuple* of *int* between 0-2.
            minrec    - *int*
            maxrec    - *int*

        """
    
        message = ('V field, recs ' + str(minrec) + ':' 
                    + str(maxrec) + ', avg ax ' + str(sumaxes))
        if (binlimits):
            message += ', bins ' + str(binlimits)
        print(message) 
    
        msum, binspaces = self.mdata.get_bins(minrec,maxrec,sumaxes=sumaxes,
                                              meantime=False,sumtime=sumtime,
                                              binlimits=binlimits)
        psum, binspaces = self.pdata.get_bins(minrec,maxrec,sumaxes=sumaxes,
                                              meantime=False,sumtime=sumtime,
                                              binlimits=binlimits)

        # Divide and patch any NaNs
        vfield = np.divide(psum,msum) 
        vfield[np.isnan(vfield)] = 0.0
        
        return vfield, binspaces

# Temperature field
class TBins():

    def __init__(self,fdir,cpol_bins):
        self.mdata = MassBins(fdir,cpol_bins)
        self.pdata = MomBins(fdir,cpol_bins)
        self.KEdata = KEBins(fdir,cpol_bins)

    def get_field(self, minrec, maxrec, sumaxes=(), peculiar=True,
                  binlimits=None): 


        """
            Get the temperature field from files Tbins, mbins and vbins 
            averaged over time records minrec->maxrec, AND averaged over 
            spatial directions specified by sumaxes.
            
            sumaxes   - *tuple* of *int* between 0-2.
            minrec    - *int*
            maxrec    - *int*
            peculiar  - take into account streaming velocity

        """
    
        print('Getting temperature field from records ' + str(minrec) + ' to ' 
              + str(maxrec) + ', averaging over axes ' + str(sumaxes) + '.')
    
        mfield, binspaces = self.mdata.get_bins(minrec, maxrec, meantime=False,
                                                sumtime=True,
                                                binlimits=binlimits)
        pfield, binspaces = self.pdata.get_bins(minrec, maxrec, meantime=False,
                                                sumtime=True,
                                                binlimits=binlimits)
        KEfield, binspaces = self.KEdata.get_bins(minrec, maxrec,
                                                  meantime=False, sumtime=True,
                                                  binlimits=binlimits)


        # Temperature (no streaming consideration)
        Tfield = np.divide(KEfield,(3.0*mfield))
        Tfield[np.isnan(Tfield)] = 0.0
        #Tfield = Tfield[:,:,:,0,:]
        Tfield = Tfield[:,:,:,0]

        # Remove average of streaming component
        if (peculiar==True):

            print('Average samples for streaming velocity: ' 
                   + str(np.mean(mfield)) )

            vfield = np.divide(pfield,mfield)
            vfield[np.isnan(vfield)] = 0.0
            v2field = np.sum((vfield**2.0),3)
            Tfield = Tfield - (1./3.)*v2field

        # Avg time
        #Tfield = np.mean(Tfield,3)

        # Avg space
        Tfield = np.mean(Tfield,sumaxes)
        
        return Tfield, binspaces

# Pressure fields
class pVABins():

    def __init__(self, fdir, fname, cpol_bins=False):
        self.fdir = fdir
        self.fname = fname
        self.cpol_bins = cpol_bins
        self.Pobj = PBins(fdir,fname,cpol_bins)

    def get_field(self, minrec, maxrec, meanaxes=(), peculiar=False,
                  binlimits=None):

        print('Getting '+self.fname+' field from recs ' + str(minrec) + ' to ' 
              + str(maxrec) + ', meanaxes = ' + str(meanaxes) + ', peculiar = ' 
              + str(peculiar) )

        # Read raw data file    
        Pfield, binspaces = self.Pobj.get_field(minrec, maxrec,
                                                binlimits=binlimits)  

        # Take off square of peculiar momenta if specified
        if (peculiar==True):


            if (self.fname=='pVA_c'):
                message = ('\n *** \n Removing the peculiar velocity from '
                +' the configurational part \n of the stress tensor is '
                +' entirely nonsensical! I will ignore this instruction.\n'
                +' ***\n')
                print(message)
        
            else:   

                # Get mean velocity and density field
                vData = VBins(self.fdir,cpol_bins=self.cpol_bins)
                dData = DensityBins(self.fdir,cpol_bins=self.cpol_bins)
                vfield, binspaces = vData.get_field(minrec,maxrec,
                                                    binlimits=binlimits)
                dfield, binspaces = dData.get_field(minrec,maxrec,
                                                    binlimits=binlimits)

                # Find outer product of v*v and reshape to 1x9 rather than 3x3
                vvfield = np.einsum('...j,...k->...jk',vfield,vfield)
                vvshapelist = list(vvfield.shape)
                newshape = tuple(vvshapelist[0:-2]+[9])
                vvfield  = np.reshape(vvfield,newshape)
    
                # Remove square of streaming velocity
                Pfield = Pfield - dfield*vvfield

        # Find the mean over the axes specified by the user
        Pfield = np.mean(Pfield,axis=meanaxes)

        return Pfield, binspaces



# Pressure fields
class CV_data():

    def __init__(self,fdir,cpol_bins=False):
        self.fdir = fdir
        self.cpol_bins = cpol_bins
        self.fluxobj   = CV(fdir,'vflux',cpol_bins)
        self.stressobj = CV(fdir,'psurface',cpol_bins)

    def get_field(self,minrec,maxrec,meanaxes=(),peculiar=False,    
                  binlimits=None):

        print('Getting '+'vflux & stress'+' fields from recs ' + str(minrec) +
              ' to ' + str(maxrec) + ', meanaxes = ' + str(meanaxes))

        # Read raw data file    
        flux, binspaces = self.fluxobj.get_field(minrec,maxrec,
                                                 binlimits=binlimits)  
        stress, binspaces = self.stressobj.get_field(minrec,maxrec,
                                                     binlimits=binlimits)     

        return flux, stress, binspaces


class DensityBins():

    def __init__(self,fdir,cpol_bins):
        self.mdata = MassBins(fdir,cpol_bins)
        self.header = HeaderData(open(fdir+'simulation_header','r'))
        self.Nmass_ave = int(self.header.Nmass_ave)

    def get_field(self, minrec, maxrec, meanaxes=(), binlimits=None):

        print('Getting density field from recs ' + str(minrec) + ' to ' 
              + str(maxrec) + ', meanaxes = ' + str(meanaxes))

        msum, binspaces = self.mdata.get_bins(minrec, maxrec, meantime=True,
                                              sumtime=False, 
                                              binlimits=binlimits)
        mfield  = np.divide(msum,self.Nmass_ave)

        binvolumes = self.mdata.get_binvolumes(binlimits=binlimits)
        density = np.divide(mfield,binvolumes)

        # Take an average over the meanaxes if specified
        density = np.mean(density,axis=meanaxes)
        
        return density, binspaces
