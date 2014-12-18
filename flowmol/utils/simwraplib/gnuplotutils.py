#!/usr/bin/env python

class GnuplotUtils:

    def __init__(self,gnuplotscript):
        self.gnuplotscript = gnuplotscript

    def specify_outputfile(self,rundir,outfilename,extracommands=None):

        """ 
            Specify the name of the gnuplot output file 
            in the gnuplot script
        """
        try:
            #Attempt to get name of gnuplot output
            inf = open(self.gnuplotscript,'r')
            outf = open(rundir+self.gnuplotscript,'w')

            #Rewrite input file in each run directory with required modifications
            for line in inf:
                #Replace output name with specified name
                if ('set output' in line):
                    outf.write('set output ' + '"' + outfilename + '"\n')
                #Plot is last line in file
                elif('plot' in line):
                    #Functionality to add extra command lines
                    if extracommands:
                        for line in extracommands:
                            outf.write(line)
                    # Write set output file before plot
                    outf.write('set output ' + '"' + outfilename + '"\n')
                    outf.write(line)
                else:
                    #Write line
                    outf.write(line)
            inf.close()
            outf.close()
        except IOError as e:
            print "Gnuplot file I/O error({0}): {1}".format(e.errno, e.strerror)
