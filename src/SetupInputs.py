
import pprint
import re

class SetupInputs():

    """
        Create a dictionary like object
        from flowmol setup_read_input
        which has all keywords and inputs which 
        can be read

        Format is as follows

            InputsDict = {}
            for key, item in FlowmolInputDict.items():
                try:
                    InputsDict[key] = {}
                    helpstr = FlowmolInputs.get_helpstring(key)
                    InputsDict[key]["HELP"] = helpstr
                    InputsDict[key]["vars"] = {i:"0" for i in item}
                    #InputsDict[key] = {"vars":{i:"0" for i in item}}
                    print(key, InputsDict[key])
                except KeyError:
                    print("key ", key, " not found")

        Examples include:

            InputsDict = {"INPUT":{"HELP":"THis is text to describe variable", "vars":{"name":"2"}}, 
                         "THING":{"HELP":"different help text", "vars":{"xcells":"1","ycells":"2","zcells":"3"}},
                          "STR":{"HELP":"Example of a string with a really long help example", 
                           "vars":{"string":"Hello","logical":".true.","int":"2","float":"3.14159",
                           "List":{"names":['NVE', 'NVT', 'Tag Move system'], "numbers":[0,1,6]}}}}

    """

    def __init__(self, fsetup='./setup_read_input.f90'):

        with open(fsetup,'r') as f:
            out = f.read()

        self.listout = out.split('\n')
        self.dictinput = {}

    def __str__(self):
        pprint.pprint(self.dictinput)

    def remove_dup(self,seq):
        """
            Remove duplicates in list
        """
        seen = set()
        seen_add = seen.add
        return [ x for x in seq if x not in seen and not seen_add(x)]

    def get_values(self):

        self.values = []
        for item in self.listout:
            if (item.find('read(1') != -1):
                #print(re.split("read\(.*?\)",item)[-1].strip(' ').split('\t')[0].split("!")[0])
                self.values.append(re.split('\*\\)|ios\\)',
                                 item)[-1].strip(' ').split('\t')[0])

        #pprint.pprint(self.values)
        return self.remove_dup(self.values)

    def get_keys(self):

        self.keys = []
        for item in self.listout:
            if (item.find(' locate') != -1):
                self.keys.append(re.split(',',item)[1].strip("'"))

        #pprint.pprint(self.keys)
        return self.remove_dup(self.keys)

    def get_items(self):

        keyword = 'default'; located = False
        for item in self.listout:
            if (item.find(' locate') != -1):
                keyword = re.split(',',item)[1].strip("'")
                variables = []; located = True
            if (item.find('read(1') != -1 and located):
                variables.append(re.split('\*\\)|ios\\)',
                                 item)[-1].strip(' ').split('\t')[0].split("!")[0].strip(" "))
                #print(item,keyword,variables)
                self.dictinput[keyword] = variables

        #pprint.pprint(self.dictinput)
        return self.dictinput

    def get_optional(self, key):
        for item in self.listout:
            if (item.find(' locate') != -1):
                if (item.find(key) != -1):
                    if (".true." in item):
                        optional = True
                    elif (".false." in item):
                        optional = False
                    else:
                        raise IOError("Argument should be optional or not", item)
                    break
        return optional

    def get_helpstring(self, key):
        for n in range(len(self.listout)):
            if (self.listout[n].find(' locate') != -1):
                if (self.listout[n].find(key) != -1):
                    break
        misscount = 5; txtfound = False
        for r in range(n, 0, -1):
            if (self.listout[r].find('!#') != -1 or 
                self.listout[r].find('! #') != -1):
                txtfound = True
                if ("###" in self.listout[r]):
                    break
            else:
                misscount -= 1
            if (misscount < 0):
                break

        helpstr = ""
        if txtfound:
            for i in range(r, n):
                helpstr += self.listout[i] + "\n"

        return helpstr

    def variables_from_string(self, string):
        out = []
        for s in string.split("\n"):
            if "[" in s and "]" in s:
                varnum = int(s.split("[")[1].split("]")[0])
                splt = s.split("]")[1].split("-")
                try:
                    #Integer in numbered input
                    val = int(splt[0])
                except ValueError:
                    #Or store as string
                    val = splt[0]
                name = "-".join(splt[1:])
                out.append([varnum, val, name])
        return out

if __name__ == "__main__":
    a = SetupInputs()
    b = a.get_keys()
    c = a.get_items()
    d = a.get_values()

    #print(d)

    for key in b:
        try:
            print(key,c[key],  a.get_optional(key))
            helpstr = a.get_helpstring(key)
            vars = a.variables_from_string(helpstr)
            print(helpstr, vars)
        except KeyError:
            print("key ", key, " not found")

        
           



   

