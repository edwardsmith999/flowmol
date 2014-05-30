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
        self.labels = ["mag"]
        MDField.__init__(self,fdir,cpol_bins=cpol_bins)
        self.nperbin = self.Raw.nperbin
        self.plotfreq = self.Raw.header.Nmass_ave

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
        self.labels = ["x","y","z"]
        MDField.__init__(self,fdir,cpol_bins=cpol_bins)
        self.nperbin = self.Raw.nperbin
        self.plotfreq = self.Raw.header.Nvel_ave

class MD_FField(MDField):

    """
        MD_FField manages Body Force field data in the form of
        applied Force per molecule summed with 3 double
        precision real data type per bin
        e.g. fnames = [Fext] (default in [])
    """

    dtype = 'd'
    nperbin = 3

    def __init__(self,fdir,fname='Fext',cpol_bins=False):

        self.fname = fname
        self.labels = ["x","y","z"]
        MDField.__init__(self,fdir,cpol_bins=cpol_bins)
        self.nperbin = self.Raw.nperbin
        self.plotfreq = self.Raw.header.Nvel_ave


class MD_EField(MDField):

    """
        MD_TField manages energy field data in the form of
        molecular velocity squared and potential energy with 1 
        double precision real data type per bin            
        e.g. fnames = [Tbins], esnap, Fvext (default in [])
    """
    
    dtype = 'd'
    nperbin = 1

    def __init__(self,fdir,fname='Tbins',cpol_bins=False):

        self.fname = fname
        self.labels = ["mag"]
        MDField.__init__(self,fdir,cpol_bins=cpol_bins)
        self.nperbin = self.Raw.nperbin
        self.plotfreq = self.Raw.header.NTemp_ave

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
        self.labels = ["xtop","ytop","ztop","xbottom","ybottom","zbottom"]
        MDField.__init__(self,fdir,cpol_bins=cpol_bins)
        self.nperbin = self.Raw.nperbin
        self.plotfreq = self.Raw.header.Nmflux_ave

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
        self.labels = ["xx","xy","xz","yx","yy","yz","zx","zy","zz"]
        if (fname in ("pVA","pVA_k","pVA_c")):
            MDField.__init__(self,fdir,cpol_bins=cpol_bins)
        else:
            quit("Output type not recognised, should be pVA, pVA_k or pVA_c")
        self.nperbin = self.Raw.nperbin
        self.plotfreq = self.Raw.header.Nstress_ave

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
            self.labels = ["xxtop","yxtop","zxtop",
                           "xytop","yytop","zytop",
                           "xztop","yztop","zztop",
                           "xxbottom","yxbottom","zxbottom",
                           "xybottom","yybottom","zybottom",
                           "xybottom","yybottom","zzbottom"]
            MDField.__init__(self,fdir,cpol_bins=cpol_bins)
            self.nperbin = self.Raw.nperbin
            self.plotfreq = self.Raw.header.Nvflux_ave
        else:
            quit("Output type not recognised, should be psurface, vflux or total")


class MD_efluxField(MDField):

    """
        MD_efluxField manages energy flux field data in the form of
        energy/esurface sum over 6 cubic bin surfaces with 6 double 
        precision real data types per bin
        e.g. fnames = totalflux, vflux, psurface (no default)
    """

    dtype = 'd'
    nperbin = 6

    def __init__(self,fdir,fname,cpol_bins=False):

        if (fname in ("esurface","eflux")):
            self.fname = fname
            self.labels = ["xtop","ytop","ztop",
                           "xbottom","ybottom","zbottom"]
            MDField.__init__(self,fdir,cpol_bins=cpol_bins)
            self.nperbin = self.Raw.nperbin
            self.plotfreq = self.Raw.header.Neflux_ave
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
        self.nperbin = self.pField.nperbin
        self.labels = self.pField.labels
        if (self.mField.plotfreq == self.pField.plotfreq):
            self.plotfreq = self.pField.plotfreq
        else:
            quit("Error in MD_vfield -- Nmass_ave differs from Nvel_ave")

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
        self.nperbin = self.PField.nperbin
        self.labels = self.PField.labels
        self.plotfreq = self.PField.plotfreq

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
        self.KEField = MD_EField(fdir,cpol_bins=cpol_bins)
        Field.__init__(self,self.KEField.Raw)
        self.nperbin = self.KEField.nperbin
        if ((self.mField.plotfreq == self.pField.plotfreq) &
            (self.mField.plotfreq == self.KEField.plotfreq)):
            self.plotfreq = self.KEField.plotfreq
            self.labels = self.KEField.labels
        else:
            quit("Error in MD_Tfield -- Nmass_ave differs from Nvel_ave and/or NTemp_ave")

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

    def averaged_data(self,startrec,endrec,peculiar=True,avgaxes=(),**kwargs):
        
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
        self.nperbin = self.mField.nperbin
        self.plotfreq = self.mField.plotfreq
        self.labels = self.mField.labels
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
        self.nperbin = self.pField.nperbin
        self.plotfreq = self.pField.plotfreq
        self.labels = self.pField.labels
    def read(self, startrec, endrec,**kwargs):

        binvolumes = self.pField.Raw.get_binvolumes()
        binvolumes = np.expand_dims(binvolumes,axis=-1)

        # Read 4D time series from startrec to endrec
        pdata = self.pField.read(startrec, endrec,**kwargs)
        pdata = np.divide(pdata,float(self.plotfreq))

        momdensity = np.divide(pdata,binvolumes)
        
        return momdensity

    def averaged_data(self,startrec,endrec,avgaxes=(),**kwargs):

        binvolumes = self.pField.Raw.get_binvolumes()
        binvolumes = np.expand_dims(binvolumes,axis=-1)

        # Read 4D time series from startrec to endrec
        pdata = self.pField.read(startrec, endrec,**kwargs)
        pdata = np.divide(pdata,float(self.plotfreq))

        if (avgaxes != ()):
            pdata = np.sum(pdata,axis=avgaxes) 
            # binvolumes should only be length=1 in time & component axis 
            binvolumes = np.sum(binvolumes,axis=avgaxes) 
        
        momdensity = np.divide(pdata,binvolumes)

        return momdensity 

class MD_CVvField(MDField):

    def __init__(self,fdir,cpol_bins=False):
        self.CVvsnapField = MD_pField(fdir='vsnap',cpol_bins=cpol_bins)
        self.CVvfluxField = MD_pfluxField(fdir='vflux',cpol_bins=cpol_bins)
        self.CVpsurfField = MD_pfluxField(fdir='psurface',cpol_bins=cpol_bins)

    def check_conservation(self,startrec,endrec,**kwargs):

        CVvsnap = self.CVvsnapField(startrec,endrec+1,**kwargs)
        CVfluxdata = self.CVvfluxField.read(startrec,endrec+1,**kwargs)
        CVsurfacedata = self.CVpsurfField.read(startrec,endrec+1,**kwargs)
        Fext = self.MD_FField.read(startrec,endrec+1,**kwargs)

        CVfluxdata = np.reshape(CVfluxdata,[ self.nbins[0],
                                             self.nbins[1],
                                             self.nbins[2],
                                            3 , 6, self.nrecs ], order='F')

        CVsurfacedata = np.reshape(CVsurfacedata,[ self.nbins[0],
                                                   self.nbins[1],
                                                   self.nbins[2],
                                                  3 , 6, self.nrecs ], order='F')

        dmvdt = np.zeros(self.nbins[0],self.nbins[1],self.nbins[2],3,nrecs)
        for n in range(startrec,enrec+1):
            dmvdt[:,:,:,n] = CVvsnap[:,:,:,:,n+1] - CVvsnap[:,:,:,:,n]

        totalflux[:,:,:,:] = ( (CVfluxdata[:,:,:,:,4,:] - CVfluxdata[:,:,:,:,1,:]) 
                              +(CVfluxdata[:,:,:,:,5,:] - CVfluxdata[:,:,:,:,2,:])
                              +(CVfluxdata[:,:,:,:,6,:] - CVfluxdata[:,:,:,:,3,:]))

        totalforce[:,:,:,:] =( (CVsurfacedata[:,:,:,:,4,:] - CVsurfacedata[:,:,:,:,1,:]) 
                              +(CVsurfacedata[:,:,:,:,5,:] - CVsurfacedata[:,:,:,:,2,:])
                              +(CVsurfacedata[:,:,:,:,6,:] - CVsurfacedata[:,:,:,:,3,:]))

        conserved =  dmvdt + totalflux + totalforce + Fext

        return conserved

