#!/usr/bin/env python
import os

class InputMod:

    def __init__(self,filename):
        self.filename = filename
    
    def replace_input(self,keyword,keyvals):    

        """ 
            Replace N values underneath the first appearance of
            keyword with keyvals (length N)

        """

        found = False
        key_lno = 0 # Keyword linenumber

        for line in open(self.filename):
        
            key_lno += 1

            # Take into account keyword might be "turned off" by a
            # comment character before it
            if (line[0:len(keyword)]   == keyword or 
                line[1:len(keyword)+1] == keyword ): 

                # Mark the keyword as found 
                found = True

                # Ensure keyword is activated (i.e. not commented out)
                sedstr = ( "sed -i '" + str(key_lno) + "s/.*/" + keyword + 
                           "/' " + self.filename ) 
                os.system(sedstr)

                # Values start on next line
                val_lno = key_lno + 1

                if type(keyvals) is list:

                    for val in keyvals:

                        if (val != None):

                            sedstr = ( "sed -i '" + str(val_lno) + "s/.*/" + 
                                        str(val) + "/' " + self.filename ) 
                            os.system(sedstr)
                    
                        val_lno += 1

                else:

                    sedstr = ( "sed -i '" + str(val_lno) + "s/.*/" + 
                                str(keyvals) + "/' " + self.filename ) 
                    os.system(sedstr)

                # Stop looping through the file
                break
        
        if ( found == False ):

            with open(self.filename,'a') as f:
                f.write(keyword+'\n')
                for keyval in keyvals:
                    f.write(keyval+'\n')
    
            quit('Input string ' + keyword + 
                 ' not found, appended to file instead.')
