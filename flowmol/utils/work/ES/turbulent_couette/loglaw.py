# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np

A = 1/0.41; B = 5.1
y = np.arange(0.0001,1000,0.001)
u = A*np.log(y) + B
fig = plt.figure()
ax2 = fig.add_subplot(111)
plt.plot(y,u)
ax2.set_xscale('log')
plt.show()
