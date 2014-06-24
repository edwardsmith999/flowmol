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
        subdata = ChannelflowField.read(self,startrec,endrec,**kwargs)

        def add_laminar(vin):
            #ny = np.linspace(1,self.Raw.Ny,self.Raw.Ny)
            #v_laminar = 2.0*ny/self.Raw.Ny-1.0
            v_laminar = self.Raw.cosinegrid(a=-1.0, b=1.0, Npoints=self.Raw.nry)
            return vin - v_laminar

        v = np.apply_along_axis(add_laminar,1,subdata)
        return v 

# ============================================================================
# Complex fields that require extra calculations. 
