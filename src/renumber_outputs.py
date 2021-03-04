import os
import glob

fdir = "./interface/results/"
files = glob.glob(fdir + "*.000*")
files.sort()
no = 0; name =''
for f in files:
    namep = name
    name = f.split(".")[1].split("/")[-1]
    nop = no
    no= int(f.split(".")[-1])
    if (no-nop != 1):
        if (name == namep):
            shift = no-nop+1
            print(f, nop-no)
        else:
            shift = 0
    newf = name+".{:007d}".format(no-shift)
    #print("MUST TURN ON os.replace", f, fdir+newf)
    os.replace(f, fdir+newf)


