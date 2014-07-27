#! /usr/bin/env python
import numpy as np

from field import Field
from serial_cfdrawdata import Serial_CFD_RawData

# ============================================================================
# CFDField base class

class Serial_CFDField(Field):

    def __init__(self,fdir):
        Raw = Serial_CFD_RawData(fdir, self.fname, self.dtype, 
                         self.nperbin)
        Field.__init__(self,Raw)
        self.header = self.Raw.header
        self.axislabels = ['x','y','z']

# ============================================================================
# CFDField derived classes, but calculated by the main code
class Serial_CFD_vField(Serial_CFDField):

    dtype = 'd'
    nperbin = 3

    def __init__(self,fdir,fname='continuum_vbins'):

        self.fname = fname
        Serial_CFDField.__init__(self,fdir)
        self.labels = self.axislabels
        self.nperbin = self.Raw.nperbin
        self.plotfreq = self.Raw.header.continuum_tplot
        assert self.Raw.nperbin > 0
        self.labels = ['u','v','w']


class Serial_CFD_StressField(Serial_CFDField):

    dtype = 'd'
    nperbin = 18

    def __init__(self,fdir,fname='continuum_vbins'):

        if (fname in ("continuum_tau_xx")):
            self.fname = fname
            Serial_CFDField.__init__(self,fdir)
            self.labels = ["xxtop","yxtop","zxtop",
                           "xytop","yytop","zytop",
                           "xztop","yztop","zztop",
                           "xxbottom","yxbottom","zxbottom",
                           "xybottom","yybottom","zybottom",
                           "xybottom","yybottom","zzbottom"]
            Serial_CFDField.__init__(self,fdir)
            self.nperbin = self.Raw.nperbin
            self.plotfreq = self.Raw.header.continuum_tplot

