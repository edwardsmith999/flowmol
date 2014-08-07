#! /usr/bin/env python
import numpy as np 
import glob
import os
import sys

from headerdata import MDHeaderData
from pplexceptions import DataNotAvailable

"""

    MD_RawData Class
    Author: David Trevelyan, April 2013
    Updated: February 2014

    The MD_RawData class is associated with both a data file (e.g.
    mbins, pVA, etc.) and a results directory in which the 
    simulation_header is located. This header contains necessary 
    information that is used to reshape the 1D array of data read 
    from the file into a format that is easy to process later on.
    
    MD_RawData can read any binary data file from the MD code, but
    requires knowledge of the data type (i.e. integer or real) and 
    the number of values per averaging "bin" to read (i.e. velocity
    would require 3 values per bin, one for each cartesian direction).
    
    The header variable cpol_bins is true if the data has been
    averaged in the MD code using cylindrical polar bins: in this 
    case the only difference lies in creating the mesh of the 
    bin topology, because the domain size in cylindrical polar
    coordinates is different to the cartesian domain size that is
    written to the simulation_header.

    This class is designed for use in two ways:
        
        1) As a contained object within another that packages
           and returns the data in a 3D field format (from which 
           it may be averaged/summed/sliced into a new format that 
           you may want to plot), i.e. an MDField object, or

        2) For inspection of the numerical values in any binary
           data file in a results folder from the MD code.

"""

class MD_RawData:
    
    def __init__(self,fdir,fname,dtype,nperbin):

        """
            fdir       -  file directory containing results, string
            fname      -  file path from which to read raw data, string
            dtype      -  datatype string, 'i' for integer, 'd' for float
            nperbin    -  number of items to read per bin, integer
        """

        if (fdir[-1] != '/'): fdir += '/' 
        self.fdir = fdir
        self.fname = fname

        if (glob.glob(fdir+fname)):
            self.separate_outfiles = False
        elif (glob.glob(fdir+fname+'.*')):
            self.separate_outfiles = True 
        else:
            print('Neither ' + fname + ' nor ' + fname + '.* exist.')
            raise DataNotAvailable

        self.header = MDHeaderData(fdir)
        try:
            self.cpol_bins = bool(int(self.header.cpol_bins))
        except AttributeError:
            self.cpol_bins = False
            self.header.cpol_bins = False
        self.dtype = dtype
        self.nperbin = nperbin
        self.nbins, self.grid = self.get_bintopology()
        self.maxrec = self.get_maxrec()

    def get_bintopology(self):

        """
            Returns:
            
                gnbins    - A length-3 list of the number of bins in each
                            direction, and
                binspaces - A length-3 list of numpy linspaces specifying
                            the locations of the center of each bin in a
                            uniform grid (one linspace for each direction)

        """
        
        gnbins  = ([ int(self.header.gnbins1), 
                     int(self.header.gnbins2),
                     int(self.header.gnbins3) ])

        if (self.cpol_bins == True):
            domain = ([ float(self.header.r_io) - float(self.header.r_oi), 
                        2.0*np.pi,
                        float(self.header.globaldomain3) ])
        else:
            domain = ([ float(self.header.globaldomain1),
                        float(self.header.globaldomain2),
                        float(self.header.globaldomain3) ])

        binspaces = [] 
        for ixyz in range(3):
            binsize = np.divide(domain[ixyz],gnbins[ixyz])
            botbincenter = binsize/2.0 
            topbincenter = gnbins[ixyz]*binsize - binsize/2.0
            binspaces.append(np.linspace(botbincenter,
                                       topbincenter,
                                       num=gnbins[ixyz]))

        return gnbins, binspaces

    def get_binvolumes(self,binlimits=None):

        binspaces = self.grid
    
        if (self.cpol_bins == True):    

            r_oi = float(self.header.r_oi)
            r, theta, z = np.meshgrid((binspaces[0]+r_oi),
                                       binspaces[1],
                                       binspaces[2],
                                       indexing='ij')


            dr     = binspaces[0][1] - binspaces[0][0]
            dtheta = binspaces[1][1] - binspaces[1][0]
            dz     = binspaces[2][1] - binspaces[2][0]

            r = r + 0.5*dr
            binvolumes = r*dr*dtheta*dz

        else:

            x, y, z = np.meshgrid(binspaces[0],binspaces[1],binspaces[2],
                                  indexing='ij')


            dx = binspaces[0][1] - binspaces[0][0]
            dy = binspaces[1][1] - binspaces[1][0]
            dz = binspaces[2][1] - binspaces[2][0]

            binvolumes = np.ones(x.shape)*dx*dy*dz

        # If bin limits are specified, return only those within range
        if (binlimits):

            # Defaults
            lower = [0]*3
            upper = [i for i in binvolumes.shape] 
    
            for axis in range(3):
                if (binlimits[axis] == None):
                    continue
                else:
                    lower[axis] = binlimits[axis][0] 
                    upper[axis] = binlimits[axis][1] 

            binvolumes = binvolumes[lower[0]:upper[0],
                                    lower[1]:upper[1],
                                    lower[2]:upper[2]]
                
        # Ensure binvolumes is the right shape for subsequent
        # broadcasting with other fields
        binvolumes = np.expand_dims(binvolumes,-1)
        return binvolumes

    def get_maxrec(self):

        if (glob.glob(self.fdir+self.fname)):

            filesize = os.path.getsize(self.fdir+self.fname)
            if (self.dtype == 'i'):
                maxrec = filesize/(4*self.nperbin*np.prod(self.nbins)) - 1
            elif (self.dtype == 'd'):
                maxrec = filesize/(8*self.nperbin*np.prod(self.nbins)) - 1
            else:
                sys.exit('Unrecognised dtype in MD_RawData.get_maxrec')

        elif (glob.glob(self.fdir+self.fname+'.*')):

            filelist = glob.glob(self.fdir+self.fname+'.*')
            sortedlist = sorted(filelist)
            maxrec = int(sortedlist[-1].split('.')[-1])
            
        else:
            print('Neither ' + self.fname + ' nor ' + self.fname + '.* exist.')
            sys.exit

        return maxrec


    def read(self, startrec, endrec, binlimits=None, verbose=False, quit_on_error=True):

        """
            Required inputs:

                startrec - seek a specific record with this integer, count
                           from 0.
                endrec   - record at which to finish (integer)

            Return:
                
                bindata - 4D array of data in one record that was
                          read from the binary data file. The size
                          is (nbinsx, nbinsy, nbinsz, nperbin) or
                          the equivalent in cylindrical polar.
                
        """

        #return_zeros if data cannot be obtained
        return_zeros = False
      
        # Store how many records are to be read
        nrecs = endrec - startrec + 1 
        # Allocate enough memory in the C library to efficiently insert
        # into bindata
        recitems = np.product(self.nbins)*self.nperbin
        bindata  = np.empty(nrecs*recitems)

        # Check whether the records are written separately
        # If so
        if (self.separate_outfiles):

            # Loop through files and append data
            for plusrec in range(0,nrecs):

                filepath = self.fdir+self.fname+'.'+"%07d"%(startrec+plusrec)
                try: 
                    fobj = open(filepath,'rb')
                except:
                    if quit_on_error:
                        sys.exit('Unable to find file ' + filepath)    
                    else:
                        print('Unable to find file ' + filepath)
                        return_zeros = True

                istart = plusrec*recitems
                iend = istart + recitems
                if (verbose):
                    print('Reading {0:s} rec {1:5d}'.format(
                          self.fname,startrec+plusrec))
                if return_zeros:
                    bindata = np.zeros([ self.nbins[0],self.nbins[1],
                                         self.nbins[2],self.nperbin ,nrecs ])
                else:
                    bindata[istart:iend] = np.fromfile(fobj,dtype=self.dtype)
                    fobj.close()

       # Else
        else:

            try: 
                fobj = open(self.fdir+self.fname,'rb')
            except:
                if quit_on_error:
                    sys.exit('Unable to find file ' + self.fname)    
                else:
                    print('Unable to find file ' + self.fname)
                    return_zeros = True

            # Seek to correct point in the file
            if (self.dtype == 'i'):
                recbytes = 4*recitems
            elif (self.dtype == 'd'):
                recbytes = 8*recitems
            else:
                if quit_on_error:
                    sys.exit('Unrecognised data type in read_bins')
                else:
                    print('Unrecognised data type in read_bins')   
            seekbyte = startrec*recbytes
            fobj.seek(seekbyte)

            if (verbose):
                print('Reading {0:s} recs {1:5d} to {2:5d}'.format(
                      self.fname,startrec,endrec))

            # Get data and reshape with fortran array ordering
            if return_zeros:
                bindata = np.zeros([ self.nbins[0],self.nbins[1],
                                     self.nbins[2],self.nperbin ,nrecs ])
            else:
                bindata = np.fromfile(fobj, dtype=self.dtype,
                                      count=nrecs*recitems)  

            fobj.close()

        if (verbose):
            print('Reshaping and transposing {0:s} '.format(self.fname))

        # Reshape bindata
        bindata = np.reshape( bindata,
                             [ self.nbins[0],
                               self.nbins[1],
                               self.nbins[2],
                               self.nperbin ,
                               nrecs ],
                              order='F')
        bindata = np.transpose(bindata, (0,1,2,4,3))

        # If bin limits are specified, return only those within range
        if (binlimits):

            if (verbose):
                print('bindata.shape = {0:s}'.format(str(bindata.shape)))
                print('Extracting bins {0:s} from {1:s} '.format(
                      str(binlimits),self.fname))
            # Defaults
            lower = [0]*3
            upper = [i for i in bindata.shape] 
    
            for axis in range(3):
                if (binlimits[axis] == None):
                    continue
                else:
                    lower[axis] = binlimits[axis][0] 
                    upper[axis] = binlimits[axis][1] 

            bindata = bindata[lower[0]:upper[0],
                              lower[1]:upper[1],
                              lower[2]:upper[2], :, :]


            if (verbose):
                print('new bindata.shape = {0:s}'.format(str(bindata.shape)))

        return bindata
        
