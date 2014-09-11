import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import glob
import re

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

#Get Headers from first file
fpath = '../results/'
filenames = ['VPROFILE.0000001.DAT','CPROFILE.0000001.DAT']
for filename in filenames:
    print('Available fields in ' + filename + ':')
    data = np.genfromtxt(fpath+filename,skip_header=3,names=True)
    for name in data.dtype.names:
        print(name)
grid_x, grid_y = np.mgrid[0:1:100j, 0:1:100j]

#Get time evolving properties
filenames = ['TIMEDATA1.DAT','TIMEDATA2.DAT','ERROR.DAT']
tdata = {}
for filename in filenames:
    print('Available fields in ' + filename + ':')
    data = np.genfromtxt(fpath+filename,names=True)
    #tdata[filename.replace('.DAT','')] = data
    for name in data.dtype.names:
        print(name)
        tdata[name] = data[name]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(tdata['Time'],tdata['Xc']-1.,'-',alpha=0.7)
plt.xlabel('Time')
plt.ylabel('$x_c - 1$')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim(1e-2,11)
ax.set_xlim(1e-1,1e5)
plt.show()

#for name in tdata.keys():
#    print(name)
#    plt.plot(tdata['LogTime'],tdata[name])
#    plt.xlabel('LogTime')
#    plt.ylabel(name)
#    plt.show()

#All files
files = glob.glob('./*.DAT')
vfiles = [s for s in files if "VPROFILE" in s]
cfiles = [s for s in files if "CPROFILE" in s]
sort_nicely(vfiles)
sort_nicely(cfiles)

for filename in vfiles:
    data = np.genfromtxt(filename,skip_header=3,names=True)
    plt.plot(data['X'],data['Cax'],'-',alpha=0.3)
    plt.plot(data['X'],data['Cs'],'-',alpha=0.3)
    plt.xlim((0.0,10.0))
    plt.xlabel('X')
    plt.ylabel('Cax')
    plt.draw()
    plt.pause(0.1)
    plt.cla()

#Attempt to iterpolate from known points to contour plot
#for name in data.dtype.names:
    #grid_z0 = griddata((data['X'],data['Z']), data[name], (grid_x, grid_y), method='nearest')

#        plt.imshow(grid_z0)
#        plt.colorbar()
#        plt.title(name)
#        plt.show()
#X, Z = np.meshgrid(data['X'],data['Z'])
#plt.contour(X,Z,data['C'])

#plt.quiver(data['X'],data['Z'],data['C'],data['Ca'])

#plt.plot(data['x'],data['Z'])

