
import pprint
import re

class SetupInputs():

    def __init__(self,fsetup='../MD_dCSE/src_code/setup_read_input.f90'):

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
            if (item.find('read') != -1):
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
            if (item.find('read') != -1 and located):
                variables.append(re.split('\*\\)|ios\\)',
                                 item)[-1].strip(' ').split('\t')[0])
                #print(item,keyword,variables)
                self.dictinput[keyword] = variables

        #pprint.pprint(self.dictinput)
        return self.dictinput


if __name__ == "__main__":
    a = SetupInputs()
    b = a.get_keys()
    c = a.get_items()

    for i in b:
        print(i,c[i])



        
           



   

