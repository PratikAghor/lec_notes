"""
This code plots the streamlines of the stagnation point flow
u = \alpha x, v = -\alpha y

\psi = c => xy = c
Author: Pratik Aghor
"""
###########################################
from numpy import pi, cos, sin, exp, sqrt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
###########################################
"""
params
"""
alpha = 1
c = np.array([1, 2, 3, -1, -2, -3])
###########################################
"""
define x = [0.1, l],
"""
Nx = 100
l = 2
x = np.linspace(-l, l, num=Nx)
###########################################

###########################################
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()  # Get current axes
for i in range(0, len(c)):
    c_ = c[i]
    ax.plot(x, c_/(alpha*x), label = r"$c = $ " + str(c_))  # Plot U
    ax.set_xlabel(r'$x$', fontsize=20)  # Set x label
    ax.set_ylabel(r'$y$', fontsize=20)  # Set y label
    ax.legend(loc=1)

plt.ylim(-l, l)
fig.savefig('xy_c_streamlines.png')
###########################################
