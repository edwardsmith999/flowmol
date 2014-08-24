import glob
import numpy as np
import scipy.signal as spsi
from postproclib.pplexceptions import DataNotAvailable
from postproclib.headerdata import MDHeaderData 

class EndToEnd():

    fname = 'etev'
    
    def __init__(self, fdir):

        if (fdir[-1] != '/'): fdir += '/' 
        self.fdir = fdir
        self.fname = 'etev'

        if (glob.glob(fdir+self.fname)):
            self.separate_outfiles = False
        elif (glob.glob(fdir+self.fname+'.*')):
            self.separate_outfiles = True 
        else:
            print('Neither '+self.fname+' nor '+self.fname+'.* exist.')
            raise DataNotAvailable

        self.header = MDHeaderData(fdir)
        self.nchains = int(self.header.nchains)

        if (not self.separate_outfiles):

            filesize = os.path.getsize(self.fdir+self.fname)
            maxrec = filesize/(8*self.nchains*3) - 1

        else:

            filelist = glob.glob(self.fdir+self.fname+'.*')
            sortedlist = sorted(filelist)
            self.maxrec = int(sortedlist[-1].split('.')[-1])

    def read(self, startrec, endrec, quit_on_error=True):
       
        nrecs = endrec - startrec + 1 
        # Allocate enough memory in the C library to efficiently insert
        # into bindata
        recitems = self.nchains * 3 
        data = np.empty(nrecs*recitems)

        if (not self.separate_outfiles):

            try: 
                fobj = open(self.fdir+self.fname,'rb')
            except:
                if quit_on_error:
                    sys.exit('Unable to find file ' + self.fname)    
                else:
                    print('Unable to find file ' + self.fname)

            # Seek to correct point in the file
            recitems = self.nchains*3 # 3 for 3D
            seekbyte = startrec*8*recitems # 8 for dp
            fobj.seek(seekbyte)

            # Get data and reshape with fortran array ordering
            data = np.fromfile(fobj, dtype='d', count=nrecs*recitems)
    
        else:
            
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

                istart = plusrec*recitems
                iend = istart + recitems
                data[istart:iend] = np.fromfile(fobj,dtype='d')
                fobj.close()

        data = np.reshape(data, (self.nchains, 3, nrecs), order='F')
        return data

    def inclinations(self, axis, startrec=0, endrec=None):

        """
            Axis defines a plane, this routine calculates the mean 
            inclination of all polymer end-to-end vectors to it.
        """

        if (endrec==None):
            endrec = self.maxrec

        R = self.read(startrec, endrec)
        magR = np.sqrt(R[:,0,:]**2 + R[:,1,:]**2 + R[:,2,:]**2) 
        #print(magR[0])
        #magR = np.sqrt(np.einsum('ijk,ijk->ik',R,R))
        #print(magR[0])
        Rhat = np.divide(R, magR[:,np.newaxis,:])

        Rhatdotaxishat = Rhat[:,axis,:]
        Rhatdotaxishat = np.clip(Rhatdotaxishat,-1.0,1.0)
        theta = np.pi/2.0 - np.arccos(Rhatdotaxishat)

        return theta

    def inclinations_distribution(self, axis, bins=25, startrec=0, endrec=None):

        if (endrec == None):
            endrec = self.maxrec

        theta  = self.inclinations(axis, startrec, endrec)
        hist, edges = np.histogram(theta, bins=bins, density=True)
        # NB len(edges) = len(hist) + 1
        return hist, edges 

    def time_selfcorrelation(self, startrec=0, endrec=None, verbose=False):

        tplot = int(self.header.tplot)
        delta_t = float(self.header.delta_t)

        if (endrec==None): 
            endrec=self.maxrec

        nrecs = endrec - startrec + 1 

        R = self.read(startrec, endrec)
        Rx = R[:,0,:]
        Ry = R[:,1,:]
        Rz = R[:,2,:]

        C = np.zeros(nrecs)
        for i in range(self.nchains):

            if (verbose):
                progress = ('Calculating auto-correlation function '
                            + str(i+1) + ' of ' + str(self.nchains) 
                            + '...')                       
                print(progress)

            # Correlation is # convolution in reverse
            auto_cx = spsi.fftconvolve(Rx[i,:],Rx[i,:][::-1]) 
            auto_cy = spsi.fftconvolve(Ry[i,:],Ry[i,:][::-1]) 
            auto_cz = spsi.fftconvolve(Rz[i,:],Rz[i,:][::-1]) 
            # Only positive half
            l = len(auto_cx) 
            C += auto_cx[l/2:]+auto_cy[l/2:]+auto_cz[l/2:]
        
        C = C/float(self.nchains)
        C = C/C[0]
    
        tmin = 0.0
        tmax = tplot*delta_t*(len(C)-1)
        t = np.linspace(tmin, tmax, num=len(C))
        return t, C
        
