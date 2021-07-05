
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


    def get_conditional(self, var, maxsteps=1000):

        conditions = []; convariables = []
        for l, item in enumerate(self.listout):
            condition = None; convariable = []
            if (item.find('if ') != -1 and item.lower().find(var.lower()) != -1):

                #Look for comments and ignore
                if (item.lstrip()[0] == "!"):
                    continue

                #print(var, item, item.find(var) != -1,  item.find('if') != -1)
                if (item.find('ios') != -1):
                    continue
                # First get what condition is extracting
                # between brackets "()" handling newline "&"
                # assuming more than 10 is unlikely
                condition = item
                #print("condition = ", condition, var, item.find('if') != -1)
                for i in range(1,10):
                    if (self.listout[l+i].find('&') == -1):
                        condition += self.listout[l+i]
                    else:
                        break
                #Single space on "if (" followed by a ") then" at the line end 
                #assumed, currently true and prevents
                #issues with other brackets in conditionals
                condition = condition.replace("&","").split("if (")[1]
                if (item.find(') then') != -1):
                    condition = condition.split(") then")[0]
                else:
                    condition = condition.split(")")[0]

                if (condition == "found_in_input"):
                    continue
                #print("Conditional = ", condition)
                nestif = 1
                for i in range(1, maxsteps):
                    try:
                        line = self.listout[l+i]
                    except IndexError:
                        print("End of file with conditional ", condition, " not found, passing")
                        pass
                    #print(condition, nestif,  line)
                    #Then look for a read statement
                    if (line.find('read(') != -1):
                        convariable.append(re.split('\*\\)|ios\\)',
                                            line)[-1].strip(' ').split('\t')[0].split("!")[0].strip(" "))
                    elif (line.find(' locate') != -1):                
                        convariable.append(re.split(',',line)[1].strip("'"))

                    #Handle if statements before end of main block
                    if (line.find('if ') != -1 and line.find('then') != -1):
                        nestif += 1
                    elif (line.find('if ') != -1 and line.find('&') != -1):
                        for n in range(1,10):
                            fline = self.listout[l+i+n]
                            print(fline,fline.find('&'),  fline.find('then'))
                            if (fline.find('&') != -1):
                                pass
                            elif (fline.find('then') != -1):
                                nestif += 1
                                break
                            #else:
                            #    raise IOError("Runaway argument in SetupInputs get_conditional due to ampersand and if statements ")

                    #Go until endif
                    if (line.find('endif') != -1 or line.find('end if') != -1):
                        nestif -= 1
                        if (nestif == 0):
                            break
                    if (line.find('elseif') != -1 or line.find('else if') != -1):
                        if (nestif == 1):
                            break

                #Only add them if they are non zero
                if (condition != None and convariable != []):
                    conditions.append(condition)
                    convariables.append(convariable)

        return conditions, convariables

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

    def fortran_ifstatement(self, string, varcheck=None, varDict=None):
        logical = None
        if (".and." in string):
            for s in string.split(".and."):
                logical = logical and self.fortran_logical(s, varcheck, varDict)
        elif (".or." in string):
            for s in string.split(".or."):
                logical = logical or self.fortran_logical(s, varcheck, varDict)
        else:
            logical = self.fortran_logical(string, varcheck, varDict)
        return logical

    def fortran_logical(self, string, varcheck=None, varDict=None, debug=False):

        import operator

        if (varcheck != None and varDict != None):
            raise IOError("Either varcheck or varDict needs to be supplied")

        if (".eq." in string):
            var, con = string.split(".eq.")
            logicalop = operator.eq
        elif (".ne." in string):
            var, con = string.split(".ne.")
            logicalop = operator.ne
        elif (".gt." in string):
            var, con = string.split(".gt.")
            logicalop = operator.gt
        elif (".ge." in string):
            var, con = string.split(".ge.")
            logicalop = operator.ge
        elif (".lt." in string):
            var, con = string.split(".lt.")
            logicalop = operator.lt
        elif (".le." in string):
            var, con = string.split(".le.")
            logicalop = operator.le
        else:
            logical = None
            return logical

        if (varcheck != None):
            #Compare either integers or strings
            try:
                logical = logicalop(int(varcheck), int(con))
            except ValueError:
                logical = logicalop(varcheck, con)
            if debug:
                print("fortran logical", var.strip(), logicalop, con, varcheck,  logical)
        elif (varDict != None):
            #Compare either integers or strings
            try:
                logical = logicalop(int(varDict[var.strip()]), int(con))
            except ValueError:
                logical = logicalop(str(varDict[var.strip()]).strip().lower(), 
                                    con.strip().lower())
            if debug:
                print("fortran logical", var.strip(), con, varDict[var.strip()], logical)
        else:
            raise IOError("varcheck or varDict needs to be supplied")

        return logical

if __name__ == "__main__":
    a = SetupInputs()
    b = a.get_keys()
    c = a.get_items()
    d = a.get_values()

    #print(d)

    for key in b:
        try:
            helpstr = a.get_helpstring(key)
            vars = a.variables_from_string(helpstr)
            print("Key = ", key, c[key],  a.get_optional(key))

            #print(helpstr, vars)
        except KeyError:
            print("key ", key, " not found")

    #Next plot the conditional dependence
    print("\n\n\n\n\n")
    for var in a.get_values():
        b = a.get_conditional(var)
        #if b[0]:
        #print("Options = ", var)
        for i, z in enumerate(zip(b[0], b[1])):
                #if (z[1]  != []):
            #print(i, z)

            testDict = {}
            s = z[0]
            testDict[var] = 6
            try:
                print("Logical test", var, s, a.fortran_ifstatement(s, varDict=testDict), z[1])
            except KeyError:
                print("logical and key not matching", var, s)
            except TypeError:
                print("Not a flag", var, s)



   

