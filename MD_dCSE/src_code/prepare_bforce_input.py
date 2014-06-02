#! /usr/bin/env python
import sys
sys.path.insert(0,'../../utils')

import postproclib as ppl

fdirdefault = './results/'
fdirmessage = ('Folder of bforce_pdf records to be averaged ' +
               '[{0:s}]: '.format(fdirdefault))
fdir = raw_input(fdirmessage) or fdirdefault

pdfs_obj = ppl.MD_BForcePDFs(fdir)
startdefault = 0
enddefault = pdfs_obj.maxrec

startmessage = ('Starting record [{0:d}]: ').format(startdefault)
startrec = int(raw_input(startmessage) or startdefault)

endmessage = ('Ending record [{0:d}]: ').format(enddefault)
endrec = int(raw_input(endmessage) or enddefault)

pdfs_obj.prepare_inputfile('bforce.input', startrec, endrec)
print('Bforce PDF input file written to ./bforce.input, averaged over records'
      + ' {0:d} to {1:d}.'.format(startrec, endrec))
