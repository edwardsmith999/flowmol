#! /usr/bin/env python
import numpy as np
from Field import Field
from MDRawData import MD_RawData

# ============================================================================
# MDField base class

class MDField(Field):
       
    def __init__(self,fdir,cpol_bins=False):
        Raw = MD_RawData(fdir, self.fname, self.dtype, 
                         self.nperbin,cpol_bins)
        Field.__init__(self,Raw)

# ============================================================================
# MDField derived classes, but calculated by the main code
class MD_mField(MDField):

    fname = 'mbins'
    dtype = 'i'
    nperbin = 1

class MD_pField(MDField):

    fname = 'vbins'
    dtype = 'd'
    nperbin = 3

class MD_KEField(MDField):
    
    fname = 'Tbins'
    dtype = 'd'
    nperbin = 1

class MD_PField(MDField):

    dtype = 'd'
    nperbin = 9

    def __init__(self,fdir,fname='pVA',cpol_bins=False):

        """
            MD_PField requires the specification of a filename by the
            user, allowing any of pVA, virial, or separate kinetic
            and configurational parts to be plotted with the same
            MDField class functionality.

        """

        self.fname = fname  
        MDField.__init__(self,fdir,cpol_bins=cpol_bins)

class MD_CVField(MDField):

    dtype = 'd'
    nperbin = 18

    def __init__(self,fdir,fname,cpol_bins=False):

        """
            CV requires the specification of a filename by the
            user, allowing any of vfluxes, psurface, vsnap or F_ext
            to be plotted correctly

        """

        if (fname in ("psurface","vflux")):
            self.fname = fname    
            MDField.__init__(self,fdir,cpol_bins=cpol_bins)
        elif (fname is "total"):
            self.CVflux = MD_CVField(fdir,fname="vflux",cpol_bins=cpol_bins)
            self.CVsurface = MD_CVField(fdir,fname="psurface",cpol_bins=cpol_bins)
            self.read = self.read_both
            Field.__init__(self,self.CVsurface.Raw)
        else:
            quit("Output type not recognised, should be psurface, vflux or total")

    def read_both(self,startrec,endrec):

        CVfluxdata = self.CVflux.read(startrec,endrec)
        CVsurfacedata = self.CVsurface.read(startrec,endrec)

        # Add together
        vdata = CVsurfacedata + CVfluxdata

        return vdata 
            
# ============================================================================
# Complex fields that inherit MDField AND contain MDField objects, require 
# extra calculations. "Read" and "average_data" routines are commonly 
# overridden.
class MD_vField(MDField):

    def __init__(self,fdir,cpol_bins=False):
        self.mField = MD_mField(fdir,cpol_bins)
        self.pField = MD_pField(fdir,cpol_bins)

        Field.__init__(self,self.mField.Raw)

    def read(self,startrec,endrec):

        mdata = self.mField.read(startrec,endrec)
        pdata = self.pField.read(startrec,endrec)

        # Divide and patch any NaNs
        vdata = np.divide(pdata,mdata) 
        vdata[np.isnan(vdata)] = 0.0

        return vdata 

    def averaged_data(self,startrec,endrec,avgaxes=(),avgtime=True):
        
        # Read 4D time series from startrec to endrec
        mdata = self.mField.read(startrec,endrec)
        pdata = self.pField.read(startrec,endrec)

        # Time axis is the final (fourth) axis...
        if (avgtime):
            mdata = np.sum(mdata,axis=4)
            pdata = np.sum(pdata,axis=4)
           
        # ... so spatial axes can be averaged here without changing their
        # indices.
        if (avgaxes != ()):
            mdata = np.sum(mdata,axis=avgaxes) 
            pdata = np.sum(pdata,axis=avgaxes) 

        # Divide and patch any NaNs
        vdata = np.divide(pdata,mdata) 
        vdata[np.isnan(vdata)] = 0.0

        return vdata 

class MD_pVAField(MDField):

    def __init__(self, fdir, fname, cpol_bins=False):
        self.fname = fname
        self.cpol_bins = cpol_bins
        self.PField = MD_PField(fdir,fname,cpol_bins)
        Field.__init__(self,self.PField.Raw)

    def read(self,startrec,endrec,peculiar=True):

        # Read 4D time series from startrec to endrec
        Pdata = self.PField.read(startrec,endrec)  

        # Take off square of peculiar momenta if specified
        if (peculiar==True):

            if (self.fname=='pVA_c'):
                #message = ('\n *** \n Removing the peculiar velocity from '
                #+' the configurational part \n of the stress tensor is '
                #+' entirely nonsensical! I will ignore this instruction.\n'
                #+' ***\n')
                #print(message)
                pass

            else:   

                # Get mean velocity and density field
                vField = MD_vField(self.fdir,cpol_bins=self.cpol_bins)
                dField = MD_dField(self.fdir,cpol_bins=self.cpol_bins)
                vdata = vField.read(startrec,endrec)
                ddata = dField.read(startrec,endrec)

                # Find outer product of v*v and reshape to 1x9 rather than 3x3
                vvdata = np.einsum('...j,...k->...jk',vdata,vdata)
                vvshapelist = list(vvfield.shape)
                newshape = tuple(vvshapelist[0:-2]+[9])
                vvdata = np.reshape(vvdata,newshape)
    
                # Remove square of streaming velocity
                Pdata = Pdata - ddata*vvdata

        return Pdata 

class MD_TField(MDField):

    def __init__(self,fdir,cpol_bins=False):
        self.mField = MD_mField(fdir,cpol_bins)
        self.pField = MD_pField(fdir,cpol_bins)
        self.KEField = MD_KEField(fdir,cpol_bins)
        Field.__init__(self,self.KEField.Raw)

    def read(self,startrec,endrec,peculiar=True):

        mdata = self.mField.read(startrec,endrec)
        KEdata = self.KEField.read(startrec,endrec)

        # Temperature (no streaming consideration)
        Tdata = np.divide(KEdata,(3.0*mdata))
        Tdata[np.isnan(Tdata)] = 0.0

        # Remove average of streaming component
        if (peculiar==True):
            #print('Average samples for streaming velocity: ' 
            #       + str(np.mean(mfield)) )
            pdata = self.pField.read(startrec,endrec)
            vdata = np.divide(pdata,mdata)
            vdata[np.isnan(vdata)] = 0.0
            v2data= np.sum((vdata**2.0),3)
            Tdata = Tdata - (1./3.)*v2data[:,:,:,np.newaxis]

        return Tdata 

    def averaged_data(self,startrec,endrec,peculiar=True,avgaxes=(),
                      avgtime=True):
        
        # Read 4D time series from startrec to endrec
        mdata = self.mField.read(startrec,endrec)
        KEdata = self.KEField.read(startrec,endrec)

        # Consider streaming velocity
        if (peculiar):
            pdata = self.pField.read(startrec,endrec)
            vdata = np.divide(pdata,mdata)
            vdata[np.isnan(vdata)] = 0.0
            v2data = np.sum((vdata**2.0),3,keepdims=True)

        # Time axis is the final (fourth) axis...
        if (avgtime):
            mdata = np.sum(mdata,axis=4)
            KEdata = np.sum(KEdata,axis=4)
           
        # ... so spatial axes can be averaged here without changing their
        # indices.
        if (avgaxes != ()):
            mdata = np.sum(mdata,axis=avgaxes) 
            KEdata = np.sum(KEdata,axis=avgaxes) 

        # Temperature (no streaming consideration)
        Tdata = np.divide(KEdata,(3.0*mdata))
        Tdata[np.isnan(Tdata)] = 0.0

        # Remove streaming velocity
        if (peculiar):
            if (avgtime):
                v2data = np.mean(v2data,axis=4)
            if (avgaxes != ()):
                v2data = np.mean(v2data,axis=avgaxes) 
            Tdata = Tdata - (1./3.)*v2data

        return Tdata 

# Density field
class MD_dField(MDField):
    
    def __init__(self,fdir,cpol_bins=False):
        self.mField = MD_mField(fdir,cpol_bins)
        Field.__init__(self,self.mField.Raw)

    def read(self, startrec, endrec):

        binvolumes = self.mField.Raw.get_binvolumes()
        binvolumes = np.expand_dims(binvolumes,axis=-1)
        Nmass_ave = self.mField.Raw.header.Nmass_ave

        # Read 4D time series from startrec to endrec
        mdata = self.mField.read(startrec, endrec)
        mdata = np.divide(mdata,float(Nmass_ave))

        density = np.divide(mdata,binvolumes)
        
        return density

    def averaged_data(self,startrec,endrec,avgaxes=(),avgtime=True):

        binvolumes = self.mField.Raw.get_binvolumes()
        Nmass_ave = self.mField.Raw.header.Nmass_ave

        # Read 4D time series from startrec to endrec
        mdata = self.mField.read(startrec, endrec)
        mdata = np.divide(mdata,float(Nmass_ave))

        # Time axis is the final (fourth) axis...
        if (avgtime):
            mdata = np.mean(mdata,axis=4)
        # ... so spatial axes can be averaged here without changing their
        # indices.
        if (avgaxes != ()):
            mdata = np.sum(mdata,axis=avgaxes) 
            binvolumes = np.sum(binvolumes,axis=avgaxes) 
        
        density = np.divide(mdata,binvolumes)

        return density 

# Momentum density field
class MD_momField(MDField):
    
    def __init__(self,fdir,cpol_bins=False):
        self.pField = MD_pField(fdir,cpol_bins)
        Field.__init__(self,self.pField.Raw)

    def read(self, startrec, endrec):

        binvolumes = self.pField.Raw.get_binvolumes()
        binvolumes = np.expand_dims(binvolumes,axis=-1)
        N_ave = self.pField.Raw.header.Nvel_ave

        # Read 4D time series from startrec to endrec
        pdata = self.pField.read(startrec, endrec)
        pdata = np.divide(pdata,float(N_ave))

        momdensity = np.divide(pdata,binvolumes)
        
        return momdensity

    def averaged_data(self,startrec,endrec,avgaxes=(),avgtime=True):

        binvolumes = self.pField.Raw.get_binvolumes()
        N_ave = self.pField.Raw.header.Nvel_ave

        # Read 4D time series from startrec to endrec
        pdata = self.pField.read(startrec, endrec)
        pdata = np.divide(mdata,float(N_ave))

        # Time axis is the final (fourth) axis...
        if (avgtime):
            pdata = np.mean(pdata,axis=4)
        # ... so spatial axes can be averaged here without changing their
        # indices.
        if (avgaxes != ()):
            pdata = np.sum(pdata,axis=avgaxes) 
            binvolumes = np.sum(binvolumes,axis=avgaxes) 
        
        momdensity = np.divide(pdata,binvolumes)

        return momdensity 


## THIS SEEMS TO BE JUST A CONTAINER TO DO TWO AT ONCE: THIS IS MORE
## SUITED TO MD_PLOTDATA, I THINK. THE OTHER CONTAINER CLASSES ARE THERE
## FOR WHEN YOU NEED TO CALCULATE A FIELD WE DON'T CALCULATE DIRECTLY IN
## THE CODE, E.G. V = MOM/MASS
#class CV_data():
#
#    def __init__(self,fdir,cpol_bins=False):
#        self.fdir = fdir
#        self.cpol_bins = cpol_bins
#        self.fluxobj   = CV(fdir,'vflux',cpol_bins)
#        self.stressobj = CV(fdir,'psurface',cpol_bins)
#
#    def get_field(self,minrec,maxrec,meanaxes=(),peculiar=False,    
#                  binlimits=None):
#
#        print('Getting '+'vflux & stress'+' fields from recs ' + str(minrec) +
#              ' to ' + str(maxrec) + ', meanaxes = ' + str(meanaxes))
#
#        # Read raw data file    
#        flux, binspaces = self.fluxobj.get_field(minrec,maxrec,
#                                                 binlimits=binlimits)  
#        stress, binspaces = self.stressobj.get_field(minrec,maxrec,
#                                                     binlimits=binlimits)     
#
#        return flux, stress, binspaces
