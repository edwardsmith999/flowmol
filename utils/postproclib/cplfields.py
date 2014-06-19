#! /usr/bin/env python
import numpy as np
try:
    import skimage.transform as skit
    skit_imported=True
except:
    skit_imported=False
    skit_imported_fail_message = (
              'You have not been able to import the skimage.transform '
            + '\npackage, which is required for MDField resampling. Try '
            + '\nCPLField.?_both instead, where ? is profile/contour/etc, '
            + '\nin which no resampling is performed.') 

import cfdfields
import mdfields
from headerdata import HeaderData
from field import Field

class CPLField(Field):
    
    def __init__(self,fdir):
        self.fdir = fdir
        self.header = HeaderData(open(fdir+'results/coupler_header')) 
        self.md_field = self.MDFieldType(fdir+'md_data/results/') 
        self.cfd_field = self.CFDFieldType(fdir+'couette_data/')
        self.olap_cells = np.array(
            [int(self.header.icmax_olap) - int(self.header.icmin_olap) + 1,
             int(self.header.jcmax_olap) - int(self.header.jcmin_olap) + 1,
             int(self.header.kcmax_olap) - int(self.header.kcmin_olap) + 1 ])
        self.cfd_dxyz = np.array([(self.cfd_field.grid[i][1] - 
                                   self.cfd_field.grid[i][0])
                                  for i in range(3)])
        self.cpl_dxyz = np.copy(self.cfd_dxyz)
        self.md_xyzL = np.array([float(self.md_field.header.globaldomain1),
                                 float(self.md_field.header.globaldomain2),
                                 float(self.md_field.header.globaldomain3)])
        self.cfd_cells = np.array([int(len(self.cfd_field.grid[i])) 
                                  for i in range(3)])
        self.md_cells = np.array([int(len(self.md_field.grid[i])) 
                                 for i in range(3)])
        self.md_cfdcells = np.rint(
                                   np.divide(self.md_xyzL,self.cfd_dxyz)
                                  ).astype(int)
        self.cfd_halos = [0,1,0]
        self.cpl_cells = (self.md_cfdcells + self.cfd_cells - self.olap_cells - 
                         self.cfd_halos)
        self.md_grid, self.cfd_grid = self.get_grids()
        self.maxrec = self.md_field.maxrec 
        self.grid = self.get_cpl_grid()

        self.labels = self.md_field.labels
        self.axislabels = self.md_field.axislabels
        self.Raw = self.md_field.Raw

    def get_grids(self):

        md_grid = self.md_field.grid
        cfd_grid = self.cfd_field.grid

        # Add MD domain (minus olap size) to CFD y grid
        md_yL = float(self.md_field.header.globaldomain2)
        olap_yL = self.olap_cells[1] * self.cfd_dxyz[1]
        cfd_grid[1] += md_yL - olap_yL

        return md_grid, cfd_grid  

    def get_cpl_grid(self):
        gridx = (np.arange(0,self.cpl_cells[0])+0.5)*self.cpl_dxyz[0]
        gridy = (np.arange(0,self.cpl_cells[1])+0.5)*self.cpl_dxyz[1]
        gridz = (np.arange(0,self.cpl_cells[2])+0.5)*self.cpl_dxyz[2]
        grid = [gridx,gridy,gridz]
        return grid
        
    def read(self,startrec,endrec,**kwargs):
        
        if (not skit_imported):
            quit(skit_imported_fail_message)

        md_data = self.md_field.read(startrec,endrec)
        cfd_data = self.cfd_field.read(startrec,endrec)
        ndims =  md_data.shape[-1]
        nrecs = endrec - startrec + 1

        md_coarse = np.empty([self.md_cfdcells[0],self.md_cfdcells[1],
                              self.md_cfdcells[2],nrecs,ndims])
        zoom = np.divide(self.md_cfdcells.astype(float),
                         self.md_cells.astype(float))

        #import matplotlib.pyplot as plt
        for rec in range(nrecs):
            for comp in range(ndims):
                md_coarse[:,:,:,rec,comp]=skit.resize(md_data[:,:,:,rec,comp],
                                                             self.md_cfdcells)

        jstart = self.olap_cells[1] - 1 + self.cfd_halos[1]
        cfd_data = cfd_data[:,jstart:,:,:,:]
        md_coarse = md_coarse[:,:-1,:,:,:]
        #jstart = self.olap_cells[1] + self.cfd_halos[1]
        #cfd_data = cfd_data[:,jstart:,:,:,:]

        cpl_data = np.concatenate((md_coarse,cfd_data),axis=1)
        return cpl_data
    
    def profile_both(self,axis,startrec=0,endrec=None,**kwargs):
         
        avgaxes = [0,1,2,3]
        avgaxes.remove(axis)
        avgaxes = tuple(avgaxes)

        if (endrec==None): 
            endrec = self.maxrec

        md_data = self.md_field.averaged_data(startrec,endrec,avgaxes=avgaxes,
                                              **kwargs)
        cfd_data = self.cfd_field.averaged_data(startrec,endrec,
                                                avgaxes=avgaxes, **kwargs)
        return self.md_grid[axis], md_data, self.cfd_grid[axis], cfd_data 

class CPL_vField(CPLField):
    nperbin = 3
    MDFieldType = mdfields.MD_vField 
    CFDFieldType = cfdfields.CFD_vField

class CPL_PField(CPLField):
    nperbin = 9
    MDFieldType = mdfields.MD_PField 
    CFDFieldType = cfdfields.CFD_StressField
