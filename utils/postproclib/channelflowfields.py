#! /usr/bin/env python
import numpy as np
from field import Field
from channelflowrawdata import Channelflow_RawData

# ============================================================================
# CFDField base class

class ChannelflowField(Field):

    def __init__(self,fdir):
        Raw = Channelflow_RawData(fdir)
        Field.__init__(self,Raw)
        self.axislabels = ['x','y','z']

# ============================================================================
# CFDField derived classes, but calculated by the main code
class Channelflow_vField(ChannelflowField):

    nperbin = 3 
    labels = ['u','v','w']

    def read(self,startrec,endrec,**kwargs):
        subdata = ChannelflowField.read(self,startrec,endrec)
        v = subdata[:,:,:,:,0:3]
        return v 

# ============================================================================
# Complex fields that require extra calculations. 
