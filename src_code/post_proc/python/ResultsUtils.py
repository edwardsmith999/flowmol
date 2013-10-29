#! /usr/bin/env/ python
import os
import re
import shutil as sh

def ConcatenateResults(fdir,cleanup=False):

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
            if (cleanup) os.remove(fdir+filename)

        catf.close()

    return

def DismemberResults(filepath, recbytes): 
   
    victim = open(filepath, 'rb')
    bodybags = os.path.getsize(filepath)/recbytes

    for bagnumber in range(bodybags):
        bagfilename = filepath + '.' + "%07d"%bagnumber
        tissue = victim.read(recbytes)
        bag = open(bagfilename, 'wb')
        bag.write(tissue)
  
    victim.close() 
    os.remove(filepath)
