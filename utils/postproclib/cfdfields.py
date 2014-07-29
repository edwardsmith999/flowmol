#! /usr/bin/env python
import numpy as np

from field import Field
from cfdrawdata import CFD_RawData

# ============================================================================
# CFDField base class

class CFDField(Field):

    def __init__(self,fdir):
        Raw = CFD_RawData(fdir)
        Field.__init__(self,Raw)
        self.axislabels = ['x','y','z']

# ============================================================================
# CFDField derived classes, but calculated by the main code
class CFD_vField(CFDField):

    nperbin = 3 

    def __init__(self,fdir):
        CFDField.__init__(self,fdir)
        assert self.Raw.npercell > 0
        self.labels = ['u','v','w']
    
    def read(self,startrec,endrec,binlimits=None,**kwargs):
        subdata = CFDField.read(self,startrec,endrec,binlimits=binlimits,
                                **kwargs) 
        v = subdata[:,:,:,:,0:3]
        return v 
        
class CFD_PField(CFDField):

    nperbin = 1

    def __init__(self,fdir):
        CFDField.__init__(self,fdir)
        assert self.Raw.npercell > 3
        self.labels = ['p']

    def read(self,startrec,endrec,binlimits=None,**kwargs):
        subdata = CFDField.read(self,startrec,endrec,binlimits=binlimits,
                                **kwargs) 
        P = subdata[:,:,:,:,3:4]
        return P 


class CFD_StressField(CFDField):

    nperbin = 9    
    def __init__(self,fdir):
        CFDField.__init__(self,fdir)
        assert self.Raw.npercell > 4
        x = self.axislabels[0]; y = self.axislabels[1]; z = self.axislabels[2]
        self.labels = [x+x,x+y,x+z,
                       y+x,y+y,y+z,
                       z+x,z+y,z+z]

    def read(self,startrec,endrec,binlimits=None,**kwargs):
        subdata = CFDField.read(self,startrec,endrec,binlimits=binlimits,
                                **kwargs) 
        P = subdata[:,:,:,:,4:]
        return P 
# =============================================================================
# Complex fields that require extra calculations. 
class CFD_complexField(CFDField):
    pass

class CFD_mugradvField(CFD_complexField):
  
    nperbin = 9
 
    def __init__(self, fdir):
        self.vField = CFD_vField(fdir)
        CFD_complexField.__init__(self, fdir)
        x = self.axislabels[0]; y = self.axislabels[1]; z = self.axislabels[2]
        self.labels = [x+x,x+y,x+z,
                       y+x,y+y,y+z,
                       z+x,z+y,z+z]
        self.rho = None

    def set_rho(self, rho):
        self.rho = rho
        
    def read(self, startrec, endrec, binlimits=None, **kwargs):

        if (self.rho == None):
            print('CFD_mugradvField requires rho, set by ' +
                  'CFD_mugradvField.set_rho(rho).')
 
        vdata = self.vField.read(startrec, endrec, binlimits=binlimits, 
                                 **kwargs)   

        dx = self.vField.Raw.dx
        dy = self.vField.Raw.dy
        dz = self.vField.Raw.dz
        gradv = np.empty(list(vdata.shape[:-1]) + [9])
        for rec in range(gradv.shape[-2]):
            for ixyz in range(3):
                for jxyz in range(3):
                    c = 3*ixyz + jxyz
                    gradv[:,:,:,rec,c] = (
                        np.gradient(vdata[:,:,:,rec,ixyz], dx, dy, dz)[jxyz]
                    )

        nugradv = self.vField.Raw.nu*gradv
        try:
            mugradv = np.multiply(nugradv, self.rho)
            return mugradv
        except TypeError:
            print('Rho not set, returning nugradv')
            return nugradv

