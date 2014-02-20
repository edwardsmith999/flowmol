#! /usr/bin/env python
# Routines which are not specific to MD/CFD or CPL code
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
        os.chdir(self.newPath)

    def __exit__( self ):
        os.chdir( self.savedPath )
