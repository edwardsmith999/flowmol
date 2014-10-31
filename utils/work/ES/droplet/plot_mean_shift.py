"""
=============================================
Apply mean-shift clustering algorithm to MD
=============================================

"""

print(__doc__)

import numpy as np
import scipy
import matplotlib.pyplot as plt
from itertools import cycle

import MDAnalysis
from sklearn.cluster import MeanShift, estimate_bandwidth

def plot_DCD(DCDObj):
    f = plt.figure()
    ax = f.add_subplot(111)
    X = DCDObj[1][:]
    line, = ax.plot(X[:,0],X[:,1],'o',ms=0.1)
    ax.set_aspect('equal')
    ax.autoscale(tight=True)
    plt.ion()
    plt.show()
    for frame in DCDObj:
        X = frame[:]
        line.set_xdata(X[:,0])
        line.set_ydata(X[:,1])
        plt.pause(0.0000001)
        plt.draw()

###############################################################################
# Load MD data
fdir = "/home/es205/scratch/droplet/2D_e1p4/"
frame = 100
DCDObj = MDAnalysis.coordinates.DCD.DCDReader(fdir+"vmd_out.dcd")
XYZ = DCDObj[frame][:]
XY = XYZ[:,0:2]
#plot_DCD(DCDObj)


###############################################################################
# Compute clustering with MeanShift

# The following bandwidth can be automatically detected using
bandwidth = estimate_bandwidth(XY, quantile=0.15, n_samples=500)
ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=False)

#ms = MeanShift(bin_seeding=True)
ms.fit(XY)
labels = ms.labels_
cluster_centers = ms.cluster_centers_

labels_unique = np.unique(labels)
labels_positive = labels_unique[labels_unique>=0]
n_clusters_ = len(labels_positive)

print("number of estimated clusters : %d" % n_clusters_)

###############################################################################
# Plot result

plt.figure(1)
plt.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(XY[my_members, 0], XY[my_members, 1], col + '.')
    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
plt.title('Estimated number of clusters: %d' % n_clusters_)
#plt.imshow(XY)
plt.show()
