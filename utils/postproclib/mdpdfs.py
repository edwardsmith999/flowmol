# /usr/bin/env python
import glob
import numpy as np
import struct

from headerdata import MDHeaderData
from mdrawdata import MD_RawData

    
class MD_boundaryPDFs():

    """
        Abstract, to be inherited.
    """

    def __init__(self):
        self.recitems = self.nsubcells * self.histbins * self.nperbin
        self.maxrec = self.get_maxrec()
        self.set_bintopology()

    def set_bintopology(self):

        self.histbinsize = (self.maxval - self.minval)/float(self.histbins)
        minbincenter = self.minval + self.histbinsize/2.0
        maxbincenter = self.maxval - self.histbinsize/2.0
        self.F = np.linspace(minbincenter,maxbincenter,num=self.histbins)
        self.binlefts = self.F - self.histbinsize/2.0
        self.binrights = self.F + self.histbinsize/2.0

        self.subcellsize = (self.cutoff)/float(self.nsubcells)
        minbincenter = self.histbinsize/2.0
        maxbincenter = self.cutoff - self.histbinsize/2.0
        self.dy = np.linspace(minbincenter,maxbincenter,num=self.nsubcells)
    
    def get_maxrec(self):

        if (glob.glob(self.fdir+self.fname)):
            self.separate_outfiles = False
        elif (glob.glob(self.fdir+self.fname+'.*')):
            self.separate_outfiles = True 
        else:
            print('Neither ' + self.fname + ' nor ' + self.fname + '.* exist.')
            quit()

        if (not self.separate_outfiles): 
            filesize = os.path.getsize(self.fdir+self.fname)
            # 4 for integers
            maxrec = filesize/(4*self.recitems) - 1
        else:
            filelist = glob.glob(self.fdir+self.fname+'.*')
            sortedlist = sorted(filelist)
            maxrec = int(sortedlist[-1].split('.')[-1])

        return maxrec 

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
        return self.dy, self.F, pdfs 

    def prepare_inputfile(self, savefile, startrec, endrec): 
        data = self.mean_pdfs(startrec, endrec) 
        with open(savefile, 'wb') as fobj:
            fobj.write(struct.pack('i',self.nsubcells)) 
            fobj.write(struct.pack('i',self.histbins)) 
            fobj.write(struct.pack('d',self.minval)) 
            fobj.write(struct.pack('d',self.maxval)) 
            fobj.write(struct.pack('d',self.cutoff)) 
            fobj.write(struct.pack('i',self.nperbin)) 
            np.ravel(data,order='F').tofile(fobj)

class MD_bforcePDFs(MD_boundaryPDFs):

    def __init__(self,fdir):

        if (fdir[-1] != '/'): fdir += '/' 
        self.fdir = fdir
        self.fname = 'bforce_pdf'
        self.nperbin = 3
        self.header = MDHeaderData(fdir)
        self.nsubcells = int(self.header.bforce_pdf_nsubcells)
        self.histbins = int(self.header.bforce_pdf_nbins)
        self.minval = float(self.header.bforce_pdf_min)
        self.maxval = float(self.header.bforce_pdf_max)
        self.cutoff = float(self.header.bforce_pdf_cutoff)
        MD_boundaryPDFs.__init__(self)

class MD_bUPDFs(MD_boundaryPDFs):

    def __init__(self,fdir):
        if (fdir[-1] != '/'): fdir += '/' 
        self.fdir = fdir
        self.fname = 'bU_pdf'
        self.nperbin = 1
        self.header = MDHeaderData(fdir)
        self.nsubcells = int(self.header.bforce_pdf_nsubcells)
        self.histbins = int(self.header.bU_pdf_nbins)
        self.minval = float(self.header.bU_pdf_min)
        self.maxval = float(self.header.bU_pdf_max)
        self.cutoff = float(self.header.bforce_pdf_cutoff)
        MD_boundaryPDFs.__init__(self)
