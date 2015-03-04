#!/usr/bin/env python
import sys
import argparse
import os

class Chdir:          
    """
       Wrapper to move from current directory to new directory
       and return when using with

       Example usage:

       with Chdir('./../'):
           os.system('./a.out')
    """
    def __init__( self, newPath ):  
        self.savedPath = os.getcwd()
        self.newPath = newPath

    def __enter__( self ):
        os.chdir('./'+self.newPath)

    def __exit__( self, etype, value, traceback):
        os.chdir(self.savedPath )

#Keyword arguments
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('pptype',help='Postprocessing to run, either vis or vmd', default=None)
argns, unknown = parser.parse_known_args()
args = vars(argns)

scriptdir = '../../utils/'
sys.path.insert(0,scriptdir)
import run_visualiser as vis
import run_vmd as vmd

if args['pptype'] == 'vis':

    with Chdir(scriptdir):
        vis.run_visualiser(parser)

elif args['pptype'] == 'vmd':
    parser.set_defaults(fdir=os.path.realpath('./results')+'/')
    with Chdir(scriptdir):
        vmd.run_vmd(parser)

else:
    sys.exit("Error, specify vis or vmd post processing")

