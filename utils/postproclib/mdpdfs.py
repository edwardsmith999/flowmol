# /usr/bin/env python
import glob
import numpy as np
import struct

from headerdata import HeaderData
from mdrawdata import MD_RawData
    
class MD_BForcePDFs(MD_RawData):

    def __init__(self,fdir):
       
        if (fdir[-1] != '/'): fdir += '/' 
        fname = 'bforce_pdf' 
        MD_RawData.__init__(self,fdir,fname,'i',3) 

        self.nsubcells = int(self.header.bforce_pdf_nsubcells)
        self.histbins = int(self.header.bforce_pdf_nbins)
        self.minval = float(self.header.bforce_pdf_min)
        self.maxval = float(self.header.bforce_pdf_max)
        self.recitems = self.nsubcells * self.histbins * self.nperbin

        self.histbinsize = (self.maxval - self.minval)/float(self.histbins)
        self.F = []
        for ixyz in range(self.nperbin):
            minbincenter = self.minval + self.histbinsize/2.0
            maxbincenter = self.maxval - self.histbinsize/2.0
            self.F = np.linspace(minbincenter,maxbincenter,num=self.histbins)
        self.binlefts = self.F - self.histbinsize/2.0
        self.binrights = self.F + self.histbinsize/2.0
    
    def get_bintopology(self):
        return None, None
    
    def get_maxrec(self):

        if (not self.separate_outfiles): 
            print('Warning, BForcePDFs miscalculates self.maxrec when '+
                  'separate_outfiles is off, returning 0.')
            return 0 
        else:
            return MD_RawData.get_maxrec(self)

    def read(self,startrec,endrec):

        # Store how many records are to be read
        nrecs = endrec - startrec + 1 
        data  = np.empty(nrecs*self.recitems)
    
        # Check whether the records are written separately
        if (self.separate_outfiles):

            # Loop through files and append data
            for plusrec in range(0,nrecs):

                filepath = self.fdir+self.fname+'.'+"%07d"%(startrec+plusrec)
                with open(filepath,'rb') as fobj:
                    istart = plusrec * self.recitems
                    iend = istart + self.recitems
                    data[istart:iend] = np.fromfile(fobj,dtype='i')

        else:

            with open(self.fdir+self.fname,'rb') as fobj:
                recbytes = 4 * self.recitems # 4 for integers
                seekbyte = startrec*recbytes
                fobj.seek(seekbyte)
                data = np.fromfile(fobj,dtype='i',count=nrecs*self.recitems)  

        data = np.reshape(data, [self.nsubcells, self.histbins, self.nperbin, nrecs], 
                          order='F')
        return data

    def mean_pdfs(self,startrec,endrec):
        data = self.read(startrec, endrec)
        data = np.sum(data, axis=-1) 
        itgls = np.sum(data, axis=1, keepdims=True)*self.histbinsize
        pdfs = np.divide(data,itgls)
        return pdfs 

    def prepare_inputfile(self, savefile, startrec, endrec): 
        data = self.mean_pdfs(startrec, endrec) 
        with open(savefile, 'wb') as fobj:
            fobj.write(struct.pack('i',self.nsubcells)) 
            fobj.write(struct.pack('i',self.histbins)) 
            fobj.write(struct.pack('d',self.minval)) 
            fobj.write(struct.pack('d',self.maxval)) 
            fobj.write(struct.pack('i',self.nperbin)) 
            np.ravel(data,order='F').tofile(fobj)
