#! /usr/bin/env python
import numpy as np
from field import Field
from mdrawdata import MD_RawData
from pplexceptions import DataMismatch, DataNotAvailable

# ============================================================================
# MDField base class

class MDField(Field):
       
    def __init__(self,fdir):
        Raw = MD_RawData(fdir, self.fname, self.dtype, 
                         self.nperbin)
        Field.__init__(self,Raw)
        self.header = self.Raw.header
        self.cpol_bins = bool(int(self.header.cpol_bins))
        if (self.cpol_bins):
            self.axislabels = ['r','theta','z']
        else:
            self.axislabels = ['x','y','z']

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

    def __init__(self,fdir,fname='mbins'):
        
        self.fname = fname
        MDField.__init__(self,fdir)
        self.labels = ["mag"]
        self.nperbin = self.Raw.nperbin
        if fname in ['mbins','mpoly','msolv']:
            self.plotfreq = self.Raw.header.Nmass_ave
        elif fname == 'msnap':
            self.plotfreq = self.Raw.header.Nmflux_ave


class MD_pField(MDField):

    """
        MD_pField manages velocity field data in the form of
        molecular velocity summed with 3 double
        precision real data type per bin
        e.g. fnames = [vbins], vsnap (default in [])
    """

    dtype = 'd'
    nperbin = 3

    def __init__(self,fdir,fname='vbins'):

        self.fname = fname
        MDField.__init__(self,fdir)
        self.labels = self.axislabels
        self.nperbin = self.Raw.nperbin
        if fname in ['vbins','vpoly','vsolv']:
            self.plotfreq = self.Raw.header.Nvel_ave
        elif fname == 'vsnap':
            self.plotfreq = self.Raw.header.Nvflux_ave


class MD_FField(MDField):

    """
        MD_FField manages Body Force field data in the form of
        applied Force per molecule summed with 3 double
        precision real data type per bin
        e.g. fnames = [Fext] (default in [])
    """

    dtype = 'd'
    nperbin = 3

    def __init__(self,fdir,fname='Fext'):

        self.fname = fname
        MDField.__init__(self,fdir)
        self.labels = self.axislabels 
        self.nperbin = self.Raw.nperbin
        self.plotfreq = self.Raw.header.Nvel_ave


class MD_EField(MDField):

    """
        MD_EField manages energy field data in the form of
        molecular velocity squared and potential energy with 1 
        double precision real data type per bin            
        e.g. fnames = [Tbins], esnap, Fvext (default in [])
    """
    
    dtype = 'd'
    nperbin = 1

    def __init__(self,fdir,fname='Tbins'):

        self.fname = fname
        MDField.__init__(self,fdir)
        self.labels = ["mag"]
        self.nperbin = self.Raw.nperbin
        if fname == 'Tbins':
            self.plotfreq = self.Raw.header.NTemp_ave
        elif fname == 'esnap':
            self.plotfreq = self.Raw.header.Neflux_ave


class MD_mfluxField(MDField):

    """
        MD_mfluxField manages mass flux field data in the form of
        molecular count over 6 cubic bin surfaces with 6 integer 
        data types per bin            
        e.g. fnames = [mflux] (default in [])
    """
    
    dtype = 'i'
    nperbin = 6

    def __init__(self,fdir,fname='mflux'):

        self.fname = fname
        MDField.__init__(self,fdir)
        self.labels = ["xtop","ytop","ztop","xbottom","ybottom","zbottom"]
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

    def __init__(self,fdir,fname='pVA'):
        self.fname = fname
        if (fname in ("pVA","pVA_k","pVA_c")):
            MDField.__init__(self,fdir)
        else:
            print("Output type not recognised, should be pVA, pVA_k or pVA_c")
            raise DataMismatch

        x = self.axislabels[0]; y = self.axislabels[1]; z = self.axislabels[2]
        self.labels = [x+x,x+y,x+z,
                       y+x,y+y,y+z,
                       z+x,z+y,z+z]
        self.nperbin = self.Raw.nperbin
        self.plotfreq = self.Raw.header.Nstress_ave

class MD_stressField(MD_PField):
   
    """
        PField multiplied by -1, useful for coupled output

    """ 

    def read(self,startrec,endrec,**kwargs):

        if (endrec > self.maxrec):
            quit('Record ' + str(endrec) + ' is greater than the maximum '
                 'available (' + str(self.maxrec) + '). Aborting.')
        
        grid_data = -1.0 * self.Raw.read(startrec,endrec,**kwargs) 
        return grid_data

class MD_pfluxField(MDField):

    """
        MD_vfluxField manages velocity flux field data in the form of
        velocity/stress sum over 6 cubic bin surfaces with 18 double 
        precision real data types per bin
        e.g. fnames = totalflux, vflux, psurface (no default)
    """

    dtype = 'd'
    nperbin = 18

    def __init__(self,fdir,fname):

        if (fname in ("psurface","vflux")):
            self.fname = fname
            self.labels = ["xxtop","yxtop","zxtop",
                           "xytop","yytop","zytop",
                           "xztop","yztop","zztop",
                           "xxbottom","yxbottom","zxbottom",
                           "xybottom","yybottom","zybottom",
                           "xybottom","yybottom","zzbottom"]
            MDField.__init__(self,fdir)
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

    def __init__(self,fdir,fname):

        if (fname in ("esurface","eflux")):
            self.fname = fname
            self.labels = ["xtop","ytop","ztop",
                           "xbottom","ybottom","zbottom"]
            MDField.__init__(self,fdir)
            self.nperbin = self.Raw.nperbin
            self.plotfreq = self.Raw.header.Neflux_ave
        else:
            quit("Output type not recognised, should be psurface, vflux or total")

# ============================================================================
# Complex fields that inherit MDField AND contain MDField objects, require 
# extra calculations. "Read" and "average_data" routines are commonly 
# overridden. Parameters for the complex field are usually inherited from
# one of the sub-fields.

class MD_complexField(MDField):
   
    def inherit_parameters(self, subfieldobj):
        self.header = subfieldobj.Raw.header
        self.nperbin = subfieldobj.nperbin
        self.cpol_bins = subfieldobj.cpol_bins
        self.plotfreq = subfieldobj.plotfreq
        self.axislabels = subfieldobj.axislabels
        self.labels = subfieldobj.labels


class MD_vField(MD_complexField):

    def __init__(self,fdir,rectype='bins'):
        if (rectype == 'bins'):
            self.mField = MD_mField(fdir,fname='mbins')
            self.pField = MD_pField(fdir,fname='vbins')
        elif (rectype == 'snap'):
            self.mField = MD_mField(fdir,fname='msnap')
            self.pField = MD_pField(fdir,fname='vsnap')

        Field.__init__(self,self.mField.Raw)
        self.inherit_parameters(self.pField)

        if (self.mField.plotfreq == self.pField.plotfreq):
            self.plotfreq = self.pField.plotfreq
        else:
            print("Error in MD_vfield -- Nmass_ave differs from Nvel_ave")
            raise DataMismatch


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

class MD_strainField(MD_vField):

    def __init__(self,fdir,rectype='bins'):
        self.vField = MD_vField(fdir)

        Field.__init__(self,self.vField.Raw)
        self.inherit_parameters(self.vField)
        self.labels = ["dudx","dudy","dudz",
                       "dvdx","dvdy","dvdz",
                       "dwdx","dwdy","dwdz"]
        self.nperbin = 9

    def read(self,startrec,endrec, binlimits=None,**kwargs):
        vdata = self.vField.read(startrec, endrec, 
                                 binlimits=None)

        straindata = self.grad(vdata)

        if (binlimits):

            # Defaults
            lower = [0]*3
            upper = [i for i in straindata.shape] 
    
            for axis in range(3):
                if (binlimits[axis] == None):
                    continue
                else:
                    lower[axis] = binlimits[axis][0] 
                    upper[axis] = binlimits[axis][1] 

            straindata = straindata[lower[0]:upper[0],
                                    lower[1]:upper[1],
                                    lower[2]:upper[2], :, :]

        return straindata

    def averaged_data(self,startrec,endrec,avgaxes=(),**kwargs):

        straindata = self.read(startrec, endrec,**kwargs)
        
        if (avgaxes != ()):
            return np.sum(straindata,axis=avgaxes) 

class MD_vortField(MD_vField):

    def __init__(self,fdir,rectype='bins'):
        self.vField = MD_vField(fdir)
        self.strainField = MD_strainField(fdir)

        Field.__init__(self,self.vField.Raw)
        self.inherit_parameters(self.strainField)
        self.labels = ["x","y","z"]
        self.nperbin = 3

    def read(self,startrec,endrec, binlimits=None,**kwargs):
        dudr = self.strainField.read(startrec, endrec, 
                                      binlimits=None)

        vortdata = np.empty([dudr.shape[0],dudr.shape[1],
                             dudr.shape[2],dudr.shape[3],self.nperbin])
        vortdata[:,:,:,:,0] = ( dudr[:,:,:,:,7]
                               -dudr[:,:,:,:,5])
        vortdata[:,:,:,:,1] = ( dudr[:,:,:,:,2]
                               -dudr[:,:,:,:,6])
        vortdata[:,:,:,:,2] = ( dudr[:,:,:,:,3]
                               -dudr[:,:,:,:,1])

        if (binlimits):

            # Defaults
            lower = [0]*3
            upper = [i for i in vortdata.shape] 
    
            for axis in range(3):
                if (binlimits[axis] == None):
                    continue
                else:
                    lower[axis] = binlimits[axis][0] 
                    upper[axis] = binlimits[axis][1] 

            vortdata = vortdata[lower[0]:upper[0],
                                lower[1]:upper[1],
                                lower[2]:upper[2], :, :]

        return  vortdata

    def averaged_data(self,startrec,endrec,avgaxes=(),**kwargs):

        vdata = self.read(startrec, endrec, 
                           **kwargs)
        
        if (avgaxes != ()):
            return np.sum(vdata,axis=avgaxes)




class MD_dissipField(MD_vField):

    def __init__(self,fdir,rectype='bins'):
        self.vField = MD_vField(fdir)
        self.strainField = MD_strainField(fdir)

        Field.__init__(self,self.vField.Raw)
        self.inherit_parameters(self.strainField)
        self.labels = ["mag"]
        self.nperbin = 1

    def read(self,startrec,endrec, binlimits=None,**kwargs):
        dudr = self.strainField.read(startrec, endrec, 
                                      binlimits=None)

        vortdata = np.empty([dudr.shape[0],dudr.shape[1],
                             dudr.shape[2],dudr.shape[3],self.nperbin])
        vortdata[:,:,:,:,0] = ( np.power(dudr[:,:,:,:,0],2.)
                               +np.power(dudr[:,:,:,:,1],2.)
                               +np.power(dudr[:,:,:,:,2],2.))


        if (binlimits):

            # Defaults
            lower = [0]*3
            upper = [i for i in vortdata.shape] 
    
            for axis in range(3):
                if (binlimits[axis] == None):
                    continue
                else:
                    lower[axis] = binlimits[axis][0] 
                    upper[axis] = binlimits[axis][1] 

            vortdata = vortdata[lower[0]:upper[0],
                                lower[1]:upper[1],
                                lower[2]:upper[2], :, :]

        return  vortdata

    def averaged_data(self,startrec,endrec,avgaxes=(),**kwargs):

        vdata = self.read(startrec, endrec, 
                           **kwargs)
        
        if (avgaxes != ()):
            return np.sum(vdata,axis=avgaxes) 


class MD_rhouuField(MD_complexField):

    def __init__(self, fdir):

        # Get mean velocity and density field
        self.fdir = fdir
        self.vField = MD_vField(self.fdir)
        self.dField = MD_dField(self.fdir)
        Field.__init__(self,self.vField.Raw)
        self.inherit_parameters(self.vField)
        self.labels = ['uu','uv','uw',
                       'vu','vv','vw',
                       'wu','wv','ww']
        self.nperbin = 9

    def read(self,startrec,endrec,peculiar=True,**kwargs):
        vdata = self.vField.read(startrec,endrec,**kwargs)
        ddata = self.dField.read(startrec,endrec,**kwargs)

        # Find outer product of v*v and reshape to 1x9 rather than 3x3
        nrecs = endrec-startrec+1
        vvdata = np.einsum('abcdj,abcdk->abcdjk',vdata,vdata)
        vvshapelist = list(vvdata.shape)
        newshape = tuple(vvshapelist[0:4]+[9])
        vvdata = np.reshape(vvdata,newshape)

        return vvdata 

class MD_pVAField(MD_complexField):

    def __init__(self, fdir, fname):

        self.fname = fname
        self.PField = MD_PField(fdir,fname)
        Field.__init__(self,self.PField.Raw)
        self.inherit_parameters(self.PField)

    def read(self,startrec,endrec,peculiar=True,verbose=False,**kwargs):

        # Read 4D time series from startrec to endrec
        Pdata = self.PField.read(startrec,endrec,**kwargs)  

        # Take off square of peculiar momenta if specified
        if (peculiar==True):

            if (self.fname=='pVA_c'):
                if (verbose == True):
                    message = ('\n *** \n Removing the peculiar velocity from '
                    +' the configurational part \n of the stress tensor is '
                    +' entirely nonsensical! I will ignore this instruction.\n'
                    +' ***\n')
                    print(message)
                pass

            else:   

                rhovvField = MD_rhouuField(self.fdir)
                rhovvdata =  rhovvField.read(startrec,endrec,**kwargs)

                # Remove square of streaming velocity
                Pdata = Pdata - rhovvdata

        return Pdata 

class MD_TField(MD_complexField):

    def __init__(self,fdir):
        self.mField = MD_mField(fdir)
        self.pField = MD_pField(fdir)
        self.KEField = MD_EField(fdir)
        Field.__init__(self,self.KEField.Raw)
        self.inherit_parameters(self.KEField)

        if ((self.mField.plotfreq == self.pField.plotfreq) &
            (self.mField.plotfreq == self.KEField.plotfreq)):
            self.plotfreq = self.KEField.plotfreq
            self.axislabels = self.KEField.axislabels
            self.labels = self.KEField.labels
        else:
            print("Error in MD_Tfield -- Nmass_ave differs from Nvel_ave and/or NTemp_ave")
            raise DataMismatch

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


class MD_dTdrField(MD_TField):

    def __init__(self,fdir,rectype='bins'):
        self.TField = MD_TField(fdir)

        Field.__init__(self,self.TField.Raw)
        self.inherit_parameters(self.TField)
        self.labels = ["dTdx","dTdy","dTdz"]
        self.nperbin = 3

    def read(self,startrec,endrec, binlimits=None,**kwargs):
        Tdata = self.TField.read(startrec, endrec, 
                                 binlimits=None)

        dTdr = self.grad(Tdata)

        if (binlimits):

            # Defaults
            lower = [0]*3
            upper = [i for i in dTdr.shape] 
    
            for axis in range(3):
                if (binlimits[axis] == None):
                    continue
                else:
                    lower[axis] = binlimits[axis][0] 
                    upper[axis] = binlimits[axis][1] 

            dTdr = dTdr[lower[0]:upper[0],
                        lower[1]:upper[1],
                        lower[2]:upper[2], :, :]

        return dTdr

    def averaged_data(self,startrec,endrec,avgaxes=(),**kwargs):

        dTdr = self.read(startrec, endrec,**kwargs)
        
        if (avgaxes != ()):
            return np.sum(dTdr,axis=avgaxes) 

# Density field
class MD_dField(MD_complexField):
    
    def __init__(self,fdir,fname='mbins'):
        self.mField = MD_mField(fdir,fname)
        Field.__init__(self,self.mField.Raw)
        self.inherit_parameters(self.mField)

    def read(self, startrec, endrec, binlimits=None, **kwargs):

        binvolumes = self.mField.Raw.get_binvolumes(binlimits=binlimits)
        binvolumes = np.expand_dims(binvolumes,axis=-1)

        # Read 4D time series from startrec to endrec
        mdata = self.mField.read(startrec, endrec, binlimits=binlimits)
        mdata = np.divide(mdata,float(self.plotfreq))

        density = np.divide(mdata,binvolumes)
        
        return density

    def averaged_data(self,startrec,endrec,avgaxes=(),binlimits=None, **kwargs):

        nrecs = endrec - startrec + 1
        binvolumes = self.mField.Raw.get_binvolumes(binlimits=binlimits)
        binvolumes = np.expand_dims(binvolumes,axis=-1)

        # Read 4D time series from startrec to endrec
        mdata = self.mField.read(startrec, endrec, binlimits=binlimits)
        mdata = np.divide(mdata,float(self.plotfreq))

        if (avgaxes != ()):
            mdata = np.sum(mdata,axis=avgaxes) 
            # binvolumes should only be length=1 in time & component axis
            binvolumes = np.sum(binvolumes,axis=avgaxes) 

        density = np.divide(mdata,binvolumes*nrecs)

        return density 

# Momentum density field
class MD_momField(MD_complexField):
    
    def __init__(self,fdir,fname='vbins'):
        self.pField = MD_pField(fdir,fname)
        Field.__init__(self,self.pField.Raw)
        self.inherit_parameters(self.pField)

    def read(self, startrec, endrec, binlimits=None, **kwargs):

        binvolumes = self.pField.Raw.get_binvolumes(binlimits=binlimits)
        binvolumes = np.expand_dims(binvolumes,axis=-1)

        # Read 4D time series from startrec to endrec
        pdata = self.pField.read(startrec, endrec,binlimits=binlimits,**kwargs)
        pdata = np.divide(pdata,float(self.plotfreq))

        momdensity = np.divide(pdata,binvolumes)
        
        return momdensity

    def averaged_data(self,startrec,endrec,binlimits=None,avgaxes=(),**kwargs):

        binvolumes = self.pField.Raw.get_binvolumes(binlimits=binlimits)
        binvolumes = np.expand_dims(binvolumes,axis=-1)

        # Read 4D time series from startrec to endrec
        pdata = self.pField.read(startrec,endrec,binlimits=binlimits,**kwargs)
        pdata = np.divide(pdata,float(self.plotfreq))

        if (avgaxes != ()):
            pdata = np.sum(pdata,axis=avgaxes) 
            # binvolumes should only be length=1 in time & component axis 
            binvolumes = np.sum(binvolumes,axis=avgaxes) 
        
        momdensity = np.divide(pdata,binvolumes)

        return momdensity 



# MD CV Fields...

class MD_pCVField(MD_complexField):

    def __init__(self, fdir, fname):

        self.fname = fname
        self.pfluxField = MD_pfluxField(fdir,fname)
        Field.__init__(self,self.pfluxField.Raw)
        self.inherit_parameters(self.pfluxField)
        self.warnagain = True

    def read(self,startrec,endrec,peculiar=True,verbose=False,**kwargs):

        # Read 4D time series from startrec to endrec
        pflux = self.pfluxField.read(startrec,endrec,**kwargs)  

        # Take off square of peculiar momenta if specified
        if (peculiar==True):

            if (self.fname=='psurface'):
                if (verbose == True):
                    message = ('\n *** \n Removing the peculiar velocity from '
                    +' the configurational part \n of the stress tensor is '
                    +' entirely nonsensical! I will ignore this instruction.\n'
                    +' ***\n')
                    print(message)
                pass

            else:   

                rhovvField = MD_rhouuField(self.fdir)
                rhovvdata =  rhovvField.read(startrec,endrec,**kwargs)

                # Remove square of streaming velocity
                for ixyz in range(0,pflux.shape[4]):
                    pflux[:,:,:,:,ixyz] = ( pflux[:,:,:,:,ixyz] 
                                           - rhovvdata[:,:,:,:,np.mod(
                                            ixyz,rhovvdata.shape[4])])

        return pflux 



class MD_CVStressheat_Field(MD_complexField):

    def __init__(self,fdir):
        self.fdir = fdir
        self.vField = MD_vField(fdir)
        self.vfluxField = MD_pCVField(fdir,fname='vflux')
        self.psurfaceField = MD_pCVField(fdir,fname='psurface')
        Field.__init__(self,self.vField.Raw)
        self.inherit_parameters(self.vField)
        self.labels = ["Pi_dot_u","Pi_dot_v","Pi_dot_w"]
        self.nperbin = 3

    def read(self,startrec,endrec,peculiar=True,**kwargs):

        u = self.vField.read(startrec,endrec,**kwargs)
        vflux = self.vfluxField.read(startrec,endrec,peculiar=True,**kwargs)
        psurface = self.psurfaceField.read(startrec,endrec,peculiar=False,**kwargs)

        vflux = np.reshape(vflux,[ vflux.shape[0],
                                   vflux.shape[1],
                                   vflux.shape[2],
                                   vflux.shape[3], 3, 6],                                            
                                   order='F')

        psurface = np.reshape(psurface,[ vflux.shape[0],
                                         vflux.shape[1],
                                         vflux.shape[2],
                                         vflux.shape[3], 3, 6],                                            
                                         order='F')

        Pi = psurface + vflux
        Pidotu = np.einsum('abcdij,abcdi->abcdj',Pi,u)

        return Pidotu


## MD CV Fields...
#class MD_CVvField(MDField):

#    def __init__(self,fdir):
#        self.CVvsnapField = MD_pField(fdir='vsnap')
#        self.CVvfluxField = MD_pfluxField(fdir='vflux')
#        self.CVpsurfField = MD_pfluxField(fdir='psurface')

#    def check_conservation(self,startrec,endrec,**kwargs):

#        CVvsnap = self.CVvsnapField(startrec,endrec+1,**kwargs)
#        CVfluxdata = self.CVvfluxField.read(startrec,endrec+1,**kwargs)
#        CVsurfacedata = self.CVpsurfField.read(startrec,endrec+1,**kwargs)
#        Fext = self.MD_FField.read(startrec,endrec+1,**kwargs)

#        CVfluxdata = np.reshape(CVfluxdata,[ self.nbins[0],
#                                             self.nbins[1],
#                                             self.nbins[2],
#                                            3 , 6, self.nrecs ], order='F')

#        CVsurfacedata = np.reshape(CVsurfacedata,[ self.nbins[0],
#                                                   self.nbins[1],
#                                                   self.nbins[2],
#                                                  3 , 6, self.nrecs ], order='F')

#        dmvdt = np.zeros(self.nbins[0],self.nbins[1],self.nbins[2],3,nrecs)
#        for n in range(startrec,enrec+1):
#            dmvdt[:,:,:,n] = CVvsnap[:,:,:,:,n+1] - CVvsnap[:,:,:,:,n]

#        totalflux[:,:,:,:] = ( (CVfluxdata[:,:,:,:,4,:] - CVfluxdata[:,:,:,:,1,:]) 
#                              +(CVfluxdata[:,:,:,:,5,:] - CVfluxdata[:,:,:,:,2,:])
#                              +(CVfluxdata[:,:,:,:,6,:] - CVfluxdata[:,:,:,:,3,:]))

#        totalforce[:,:,:,:] =( (CVsurfacedata[:,:,:,:,4,:] - CVsurfacedata[:,:,:,:,1,:]) 
#                              +(CVsurfacedata[:,:,:,:,5,:] - CVsurfacedata[:,:,:,:,2,:])
#                              +(CVsurfacedata[:,:,:,:,6,:] - CVsurfacedata[:,:,:,:,3,:]))

#        conserved =  dmvdt + totalflux + totalforce + Fext

#        return conserved

