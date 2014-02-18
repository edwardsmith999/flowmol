#! /usr/bin/env python
import numpy as np
from Field import Field
from CFDRawData import CFD_RawData

# ============================================================================
# CFDField base class

class CFDField(Field):

    def __init__(self,fdir):
        self.Raw = CFD_RawData(fdir)
        self.grid = self.Raw.grid
        self.maxrec = self.Raw.maxrec

# ============================================================================
# CFDField derived classes, but calculated by the main code
class CFD_vField(CFDField):
    
    def read(self,startrec,endrec):
        subdata = CFDField.read(self,startrec,endrec) 
        v = subdata[:,:,:,0:3,:]
        return v 

class CFD_PField(CFDField):
    
    def read(self,startrec,endrec):
        subdata = CFDField.read(self,startrec,endrec) 
        P = subdata[:,:,:,3:,:]
        return P 

# ============================================================================
# Complex fields that require extra calculations. 
