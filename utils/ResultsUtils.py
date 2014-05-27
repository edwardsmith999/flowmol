#! /usr/bin/env/ python
import os
import re
import sys
import shutil as sh

class ResultsUtils:

    def __init__(self,outdir):
        self.outdir = outdir
        pass

    def ConcatenateResults(self,fdir,cleanup=False):

        # Compile regular expression for 7-digit integer record number
        recint = re.compile('\d{7}')

        # Define sorting key by integer in file name to be confident python 
        # will always sort by number order
        def get_int(name):
            string, integer = name.split('.')
            return int(integer)

        # Get list of files that need catting
        filedict = {}
        for filename in os.listdir(fdir):

            if (recint.search(filename)):
                string, integer = filename.split('.') 
                if (string not in filedict.keys()):
                    filedict[string] = [filename]
                else:
                    filedict[string].append(filename)

        # Cat the files
        for catfile, filelist in filedict.iteritems():

            print('Concatenating to '+catfile+'...')

            # Sort file names by integer suffix
            filelist = sorted(filelist,key=get_int)

            # Open file object for appending        
            catf = open(fdir+catfile,'ab')

            # Loop over file names and append bytes to catfile
            for filename in filelist: 

                # Open record file obj and append bytes to catf
                recf = open(fdir+filename, 'rb')
                sh.copyfileobj(recf, catf)
                recf.close()

                # Delete record file
                if (cleanup):
                    os.remove(fdir+filename)

            catf.close()

        return

def DismemberResults(self,filename,recbytes,initialrec=0,outdir='./',cleanup=False): 
   
    victim = open(filename, 'rb')
    bodybags = os.path.getsize(filename)/recbytes
    bagfilenamelist = []

    for bagnumber in range(bodybags):
        fileindex = bagnumber + initialrec
        bagfilename = filename + '.' + "%07d"%bagnumber
        tissue = victim.read(recbytes)
        bag = open(self.outdir+bagfilename, 'wb')
        bag.write(tissue)
        bagfilenamelist.append(bagfilename) 
  
    victim.close() 
    if (cleanup):
        os.remove(filename)

    return bagfilenamelist


    ##################

if __name__ == "__main__":

    message =('\nResultUtils is designed to combine or split result files \n'
                + '  Inputs should be in the form: \n'
                + '    "ResultsUtils.py -c filename" to combine files of form filename.0000001 into filename \n'
                + '    "ResultsUtils.py -d filename" to split files of form filename into filename.0000001 \n')
    print(message)
    RU = ResultsUtils()

    for arg in sys.argv[1:]:

        if (arg == '-c'):
            ix = sys.argv.index(arg)
            fdir = sys.argv[ix+1]
            if ('--cleanup' in sys.argv
                or '-cleanup' in sys.argv):
                cleanup = True
            else:
                cleanup = False
            RU.ConcatenateResults(fdir,cleanup=cleanup)

        elif (arg == '-d'):
            ix = sys.argv.index(arg)
            fdir = sys.argv[ix+1]
            recbytes = sys.argv[ix+2]
            RU.DismemberResults(fdir,recbytes)
