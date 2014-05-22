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

    """
        MD_mField manages mass field data in the form of
        molecular counts with 1 integer data type per bin 
        e.g. fnames = [mbins], msnap  (default in [])
    """

    dtype = 'i'
    nperbin = 1

    def __init__(self,fdir,fname='mbins',cpol_bins=False):

        self.fname = fname  
        MDField.__init__(self,fdir,cpol_bins=cpol_bins)

class MD_pField(MDField):

    """
        MD_pField manages velocity field data in the form of
        molecular velocity summed with 3 double
        precision real data type per bin
        e.g. fnames = [vbins], vsnap (default in [])
    """

    dtype = 'd'
    nperbin = 3

    def __init__(self,fdir,fname='vbins',cpol_bins=False):

        self.fname = fname  
        MDField.__init__(self,fdir,cpol_bins=cpol_bins)

class MD_KEField(MDField):

    """
        MD_TField manages temperature field data in the form of
        molecular velocity squared with 1 double precision 
        real data type per bin            
        e.g. fnames = [Tbins] (default in [])
    """
    
    dtype = 'd'
    nperbin = 1

    def __init__(self,fdir,fname='Tbins',cpol_bins=False):

        self.fname = fname  
        MDField.__init__(self,fdir,cpol_bins=cpol_bins)

class MD_mfluxField(MDField):

    """
        MD_mfluxField manages mass flux field data in the form of
        molecular count over 6 cubic bin surfaces with 6 integer 
        data types per bin            
        e.g. fnames = [mflux] (default in [])
    """
    
    dtype = 'i'
    nperbin = 6

    def __init__(self,fdir,fname='mflux',cpol_bins=False):

        self.fname = fname  
        MDField.__init__(self,fdir,cpol_bins=cpol_bins)

class MD_PField(MDField):

    """
        MD_PField requires the specification of a filename by the
        user, allowing any of pVA or separate kinetic pVA_k
        and configurational parts pVA_c to be plotted with the same
        MDField class functionality.
        e.g. fnames = [pVA], pVA_k, pVA_c (default in [])
    """

    dtype = 'd'
    nperbin = 9

    def __init__(self,fdir,fname='pVA',cpol_bins=False):
        self.fname = fname
        if (fname in ("pVA","pVA_k","pVA_c")):
            MDField.__init__(self,fdir,cpol_bins=cpol_bins)
        else:
            quit("Output type not recognised, should be pVA, pVA_k or pVA_c")

class MD_pfluxField(MDField):

    """
        MD_vfluxField manages velcoity flux field data in the form of
        velocity/stress sum over 6 cubic bin surfaces with 18 double 
        precision real data types per bin
        e.g. fnames = totalflux, vflux, psurface (no default)
    """

    dtype = 'd'
    nperbin = 18

    def __init__(self,fdir,fname,cpol_bins=False):

        if (fname in ("psurface","vflux")):
            self.fname = fname    
            MDField.__init__(self,fdir,cpol_bins=cpol_bins)
        else:
            quit("Output type not recognised, should be psurface, vflux or total")

# ============================================================================
# Complex fields that inherit MDField AND contain MDField objects, require 
# extra calculations. "Read" and "average_data" routines are commonly 
# overridden.
class MD_vField(MDField):

    def __init__(self,fdir,rectype='bins',cpol_bins=False):
        if (rectype == 'bins'):
            self.mField = MD_mField(fdir,fname='mbins',cpol_bins=cpol_bins)
            self.pField = MD_pField(fdir,fname='vbins',cpol_bins=cpol_bins)
        elif (rectype == 'snap'):
            self.mField = MD_mField(fdir,fname='msnap',cpol_bins=cpol_bins)
            self.pField = MD_pField(fdir,fname='vsnap',cpol_bins=cpol_bins)

        Field.__init__(self,self.mField.Raw)

    def read(self,startrec,endrec,**kwargs):

        mdata = self.mField.read(startrec,endrec,**kwargs)
        pdata = self.pField.read(startrec,endrec,**kwargs)

        # Divide and patch any NaNs
        vdata = np.divide(pdata,mdata) 
        vdata[np.isnan(vdata)] = 0.0

        return vdata 

    def averaged_data(self,startrec,endrec,avgaxes=(),**kwargs):
        
        mdata = self.mField.read(startrec,endrec,**kwargs)
        pdata = self.pField.read(startrec,endrec,**kwargs)

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
        self.PField = MD_PField(fdir,fname,cpol_bins=cpol_bins)
        Field.__init__(self,self.PField.Raw)

    def read(self,startrec,endrec,peculiar=True,**kwargs):

        # Read 4D time series from startrec to endrec
        Pdata = self.PField.read(startrec,endrec,**kwargs)  

        # Take off square of peculiar momenta if specified
        if (peculiar==True):

            if (self.fname=='pVA_c'):
                message = ('\n *** \n Removing the peculiar velocity from '
                +' the configurational part \n of the stress tensor is '
                +' entirely nonsensical! I will ignore this instruction.\n'
                +' ***\n')
                print(message)
                pass

            else:   

                # Get mean velocity and density field
                vField = MD_vField(self.fdir,cpol_bins=self.cpol_bins)
                dField = MD_dField(self.fdir,cpol_bins=self.cpol_bins)
                vdata = vField.read(startrec,endrec,**kwargs)
                ddata = dField.read(startrec,endrec,**kwargs)

                # Find outer product of v*v and reshape to 1x9 rather than 3x3
                nrecs = endrec-startrec+1
                vvdata = np.einsum('abcdj,abcdk->abcdjk',vdata,vdata)
                vvshapelist = list(vvdata.shape)
                newshape = tuple(vvshapelist[0:4]+[9])
                vvdata = np.reshape(vvdata,newshape)
    
                # Remove square of streaming velocity
                Pdata = Pdata - ddata*vvdata

        return Pdata 

class MD_TField(MDField):

    def __init__(self,fdir,cpol_bins=False):
        self.mField = MD_mField(fdir,cpol_bins=cpol_bins)
        self.pField = MD_pField(fdir,cpol_bins=cpol_bins)
        self.KEField = MD_KEField(fdir,cpol_bins=cpol_bins)
        Field.__init__(self,self.KEField.Raw)

    def read(self,startrec,endrec,peculiar=True,**kwargs):

        mdata = self.mField.read(startrec,endrec,**kwargs)
        KEdata = self.KEField.read(startrec,endrec,**kwargs)

        # Temperature (no streaming consideration)
        Tdata = np.divide(KEdata,(3.0*mdata))
        Tdata[np.isnan(Tdata)] = 0.0

        # Remove average of streaming component
        if (peculiar==True):
            #print('Average samples for streaming velocity: ' 
            #       + str(np.mean(mfield)) )
            pdata = self.pField.read(startrec,endrec,**kwargs)
            vdata = np.divide(pdata,mdata)
            vdata[np.isnan(vdata)] = 0.0
            v2data = np.sum((vdata**2.0),axis=4,keepdims=True)
            Tdata = Tdata - (1./3.)*v2data

        return Tdata 

    def averaged_data(self,startrec,endrec,peculiar=True,avgaxes=()):
        
        # Read 4D time series from startrec to endrec
        mdata = self.mField.read(startrec,endrec,**kwargs)
        KEdata = self.KEField.read(startrec,endrec,**kwargs)

        # Consider streaming velocity
        if (peculiar):
            pdata = self.pField.read(startrec,endrec,**kwargs)
            vdata = np.divide(pdata,mdata)
            vdata[np.isnan(vdata)] = 0.0
            v2data = np.sum((vdata**2.0),axis=4,keepdims=True)

        if (avgaxes != ()):
            mdata = np.sum(mdata,axis=avgaxes) 
            KEdata = np.sum(KEdata,axis=avgaxes) 

        # Temperature (no streaming consideration)
        Tdata = np.divide(KEdata,(3.0*mdata))
        Tdata[np.isnan(Tdata)] = 0.0

        # Remove streaming velocity
        if (peculiar):
            if (avgaxes != ()):
                v2data = np.mean(v2data,axis=avgaxes) 
            Tdata = Tdata - (1./3.)*v2data

        return Tdata 

# Density field
class MD_dField(MDField):
    
    def __init__(self,fdir,cpol_bins=False):
        self.mField = MD_mField(fdir,cpol_bins=cpol_bins)
        Field.__init__(self,self.mField.Raw)

    def read(self, startrec, endrec,**kwargs):

        binvolumes = self.mField.Raw.get_binvolumes()
        binvolumes = np.expand_dims(binvolumes,axis=-1)
        Nmass_ave = self.mField.Raw.header.Nmass_ave

        # Read 4D time series from startrec to endrec
        mdata = self.mField.read(startrec, endrec,**kwargs)
        mdata = np.divide(mdata,float(Nmass_ave))

        density = np.divide(mdata,binvolumes)
        
        return density

    def averaged_data(self,startrec,endrec,avgaxes=()):

        nrecs = endrec - startrec + 1
        binvolumes = self.mField.Raw.get_binvolumes()
        binvolumes = np.expand_dims(binvolumes,axis=-1)
        Nmass_ave = self.mField.Raw.header.Nmass_ave

        # Read 4D time series from startrec to endrec
        mdata = self.mField.read(startrec, endrec,**kwargs)
        mdata = np.divide(mdata,float(Nmass_ave))

        if (avgaxes != ()):
            mdata = np.sum(mdata,axis=avgaxes) 
            # binvolumes should only be length=1 in time & component axis
            binvolumes = np.sum(binvolumes,axis=avgaxes) 
        
        density = np.divide(mdata,binvolumes*nrecs)

        return density 

# Momentum density field
class MD_momField(MDField):
    
    def __init__(self,fdir,cpol_bins=False):
        self.pField = MD_pField(fdir,cpol_bins=cpol_bins)
        Field.__init__(self,self.pField.Raw)

    def read(self, startrec, endrec,**kwargs):

        binvolumes = self.pField.Raw.get_binvolumes()
        binvolumes = np.expand_dims(binvolumes,axis=-1)
        N_ave = self.pField.Raw.header.Nvel_ave

        # Read 4D time series from startrec to endrec
        pdata = self.pField.read(startrec, endrec,**kwargs)
        pdata = np.divide(pdata,float(N_ave))

        momdensity = np.divide(pdata,binvolumes)
        
        return momdensity

    def averaged_data(self,startrec,endrec,avgaxes=(),**kwargs):

        binvolumes = self.pField.Raw.get_binvolumes()
        binvolumes = np.expand_dims(binvolumes,axis=-1)
        N_ave = self.pField.Raw.header.Nvel_ave

        # Read 4D time series from startrec to endrec
        pdata = self.pField.read(startrec, endrec,**kwargs)
        pdata = np.divide(pdata,float(N_ave))

        if (avgaxes != ()):
            pdata = np.sum(pdata,axis=avgaxes) 
            # binvolumes should only be length=1 in time & component axis 
            binvolumes = np.sum(binvolumes,axis=avgaxes) 
        
        momdensity = np.divide(pdata,binvolumes)

        return momdensity 

class MD_CVField(MDField):



    def read_both(self,startrec,endrec,**kwargs):

        CVfluxdata = self.CVflux.read(startrec,endrec,**kwargs)
        CVsurfacedata = self.CVsurface.read(startrec,endrec,**kwargs)

        # Add together
        vdata = CVsurfacedata + CVfluxdata

        return vdata 
