# /usr/bin/env python
import glob
import numpy as np
import struct
import os
from headerdata import MDHeaderData
from mdrawdata import MD_RawData

    
class MD_PDFs():

    """
        Abstract, to be inherited.
    """

    def __init__(self):
        self.recitems = self.nPDFs * self.histbins * self.nperbin
        self.maxrec = self.get_maxrec()
    
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

        data = np.reshape(data, [self.nPDFs, self.histbins, 
                                 self.nperbin, nrecs], order='F')
        return data

    def mean_pdfs(self,startrec,endrec):
        data = self.read(startrec, endrec)
        data = np.sum(data, axis=-1) #Sum over all time records
        #itgls= np.sum(data, axis=1, keepdims=True)*self.histbinsize
        itgls= np.trapz(data,dx=self.histbinsize,axis=1)
        pdfs = np.empty(data.shape)
        for j in range(data.shape[0]):
            for ixyz in range(2):
                pdfs[j,:,ixyz] = np.divide(data[j,:,ixyz],itgls[j,ixyz])
        return self.dy, self.histbinloc, pdfs 

    def prepare_inputfile(self, savefile, startrec, endrec): 
        data = self.mean_pdfs(startrec, endrec) 
        with open(savefile, 'wb') as fobj:
            fobj.write(struct.pack('i',self.nPDFs)) 
            fobj.write(struct.pack('i',self.histbins)) 
            fobj.write(struct.pack('d',self.minval)) 
            fobj.write(struct.pack('d',self.maxval)) 
            fobj.write(struct.pack('d',self.cutoff)) 
            fobj.write(struct.pack('i',self.nperbin)) 
            np.ravel(data,order='F').tofile(fobj)

class MD_vPDFs(MD_PDFs):

    def __init__(self,fdir):

        if (fdir[-1] != '/'): fdir += '/' 
        self.fdir = fdir
        self.fname = 'vPDF'
        self.nperbin = 3
        self.header = MDHeaderData(fdir)
        self.PDF_flag = int(self.header.vPDF_flag)-1
        self.nPDFs = int(eval("self.header.gnbins"+str(self.PDF_flag+1)))
        self.NPDF_ave = int(self.header.NvPDF_ave)
        self.histbins = int(self.header.NPDFbins)
        self.minval = -float(self.header.PDFvlims)
        self.maxval =  float(self.header.PDFvlims)
        MD_PDFs.__init__(self)
        self.dy, self.histbinloc, self.histbinsize = self.get_bintopology()


    def get_bintopology(self):

        """
            Returns:
            
                binspaces - A length-3 list of numpy linspaces specifying
                            the locations of the center of each bin in a
                            uniform grid (one linspace for each direction)

        """
        
        #Get data values of histogram bins
        histbinsize = (self.maxval - self.minval)/float(self.histbins)
        minbincenter = self.minval + histbinsize/2.0
        maxbincenter = self.maxval - histbinsize/2.0
        histbinloc = np.linspace(minbincenter,maxbincenter,num=self.histbins)
 

        #Get spatial location of bins
        gnbins  = ([ int(self.header.gnbins1), 
                     int(self.header.gnbins2),
                     int(self.header.gnbins3) ])

        domain = ([ float(self.header.globaldomain1),
                    float(self.header.globaldomain2),
                    float(self.header.globaldomain3) ])

        ixyz = self.PDF_flag
        binsize = np.divide(domain[ixyz],gnbins[ixyz])
        botbincenter = binsize/2.0 
        topbincenter = gnbins[ixyz]*binsize - binsize/2.0
        dy=np.linspace(botbincenter,
                       topbincenter,
                       num=gnbins[ixyz])
        

        return dy, histbinloc, histbinsize


class MD_boundaryPDFs(MD_PDFs):

    def __init__(self,fdir):

        self.header = MDHeaderData(fdir)
        self.nPDFs = int(self.header.bforce_pdf_nsubcells)
        self.histbins = int(self.header.bforce_pdf_nbins)
        self.minval = float(self.header.bforce_pdf_min)
        self.maxval = float(self.header.bforce_pdf_max)
        self.cutoff = float(self.header.rcutoff)
        MD_PDFs.__init__(self)
        self.dy, self.histbinloc, self.histbinsize = self.get_bintopology()
    
    def get_bintopology(self):

        #Get data values of histogram bins
        histbinsize = (self.maxval - self.minval)/float(self.histbins)
        minbincenter = self.minval + histbinsize/2.0
        maxbincenter = self.maxval - histbinsize/2.0
        histbinloc = np.linspace(minbincenter,maxbincenter,num=self.histbins)

        #Get spatial location of bins
        minbincenter = histbinsize/2.0
        maxbincenter = self.cutoff - histbinsize/2.0
        dy = np.linspace(minbincenter,maxbincenter,num=self.nPDFs)

        return dy, histbinloc, histbinsize


class MD_bforcePDFs(MD_boundaryPDFs):

    def __init__(self,fdir):

        if (fdir[-1] != '/'): fdir += '/' 
        self.fdir = fdir
        self.fname = 'bforce_pdf'
        self.nperbin = 3
        MD_boundaryPDFs.__init__(self,fdir)


class MD_bUPDFs(MD_boundaryPDFs):

    def __init__(self,fdir):
        if (fdir[-1] != '/'): fdir += '/' 
        self.fdir = fdir
        self.fname = 'bU_pdf'
        self.nperbin = 1
        MD_boundaryPDFs.__init__(self,fdir)


