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
# ============================================================================
# Complex fields that require extra calculations. 
