import shutil as sh

from simwraplib.inpututils import InputMod

class CFDInputMod(InputMod):
    
    def __init__(self,filename):
        InputMod.__init__(self,filename) 
  
    def replace_input(self, keyword, keyval):
         
        """ 
            Replace value before all appearances of
            keyword with keyval 

        """
        f = open(self.filename,'r')
        fout = open(self.filename+'.tmp','w')
        
        lines = f.readlines()
        keywordfound = False
        for line in lines:

            try:
                name = line.split()[1]
            except IndexError:
                name = None

            if (name == keyword):
                keywordfound = True
                value = line.split()[0]
                fout.write(line.replace(value,keyval))
            else:
                fout.write(line)

        if (not keywordfound):
            print('Input string '+keyword+' not found.')
        sh.move(self.filename+'.tmp',self.filename) 
