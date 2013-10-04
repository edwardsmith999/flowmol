import string

#List of Dictonary classes with added routines to add and multiple inputs
class InputList(list):
    def __init__(self,*arg,**kw):
        super(InputList, self).__init__(*arg, **kw)

        self.valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)

    #Define Addition operation as elementwise addition
    def __add__(self, x):

        if (type(x) == InputDict):
            Listx = x.expand()
        elif (type(x) == InputList):
            Listx = x
        else:
            raise TypeError("Unsupported type " + str(type(x))  
                            + " for input addition" + 
                            " -- must be InputList or InputDict type")

        returnlist = InputList()
        ziplist = zip(self,Listx)
        for entry in ziplist:
            tempdict = {}
            for dic in entry:
                tempdict.update(dic)
            returnlist.append(tempdict)

        return returnlist

    __radd__=__add__

    #Define multiplication operation as all permutations of inputs
    def __mul__(self,x):

        if (type(x) == InputDict):
            Listx = x.expand()  
        elif (type(x) == InputList):
            Listx = x
        else:
            raise TypeError("Unsupported type " + str(type(x))  
                            + " for input multiplication" + 
                            " -- must be InputList or InputDict type")

        returnlist = InputList()
        for entry1 in self:
            for entry2 in Listx:
                newdict = {}
                newdict.update(entry1)
                newdict.update(entry2)
                returnlist.append(newdict)

        return returnlist

    __rmul__=__mul__

    #Define Addition operation as elementwise addition
    def zip_inputs(self, x):

        returnlist = self + x

        return returnlist

    #Generate all permutations of inputs and return filenames
    def outer_product_inputs(self,x):

        returnlist = self*x

        return returnlist

    def filenames(self):

        #Generate list containing filenames
        filenames=[(''.join(c for c in str(name.items()) 
                          if c in self.valid_chars[6:]))
                          for name in self]
         
        return filenames

#Dictonary class with added routines to add and multiple inputs
class InputDict(dict):
    def __init__(self,*arg,**kw):
        super(InputDict, self).__init__(*arg, **kw)

    #Expand InputDict with multiple values per entry into InputList
    def expand(self):

        expansion = self.values()[0]
        returnlist = InputList({self.keys()[0]:e} for e in expansion)
        
        return returnlist       

    #Wrapper to convert InputDict to InputList then add
    def __add__(self, x):

        # Convert to lists
        templist = self.expand()
        returnlist = templist + x

        return returnlist

    __radd__=__add__

    #Wrapper to convert InputDict to InputList then multiply
    def __mul__(self,x):

        # Convert to lists
        templist = self.expand()
        returnlist = templist * x

        return returnlist

    __rmul__=__mul__
