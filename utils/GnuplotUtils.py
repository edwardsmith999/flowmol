#!/usr/bin/env python

class GnuplotUtils:

    def __init__(self,gnuplotscript):
        self.gnuplotscript = gnuplotscript

    def specify_outputfile(self,rundir,outfilename):    

        """ 
            Specify the name of the gnuplot output file 
            in the gnuplot script
        """
        try:
            #Attempt to get name of gnuplot output
            inf = open(self.gnuplotscript,'r')
            outf = open(rundir+self.gnuplotscript,'w')

            for line in inf:
                if ('set output' in line):
                    outf.write('set output ' + '"' + outfilename + '"\n')
                elif('plot' in line):
                    outf.write('set output ' + '"' + outfilename + '"\n')
                    outf.write(line)
                else:
                    outf.write(line)
            inf.close()
            outf.close()
        except IOError as e:
            print "Gnuplot file I/O error({0}): {1}".format(e.errno, e.strerror)


            #    found = False
            #    gnuplotfile = open(value)
            #    for line in gnuplotfile:
            #        outputline = line.find('set output')
            #        if(outputline != -1):
            #            outfile=line[outputline+10:].split('"')[1]
            #            found = True
            #    if (found != True):
            #        quit('Gnuplot file does not contain line: output file "NAME.png"')
            #except IOError as e:
            #    print "Gnuplot file I/O error({0}): {1}".format(e.errno, e.strerror)
