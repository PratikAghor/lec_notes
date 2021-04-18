"""
Plotting the gauge pressure p - p0
in the thin film beneath a knife
Author: Pratik Aghor
"""
#####################
import numpy as np
from numpy import sin, cos, tan, pi, arange, zeros, sqrt, exp
import matplotlib
import matplotlib.pyplot as plt
#####################
dx = 0.01
x = arange(0, 1 + dx, dx)
eps = 1e-1
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()  # Get current axes
ax.plot(x, (1/eps**3)*(x**2-x), linewidth = 2, label=r'$p-p_{0}$')  # Plot sigma vs k
ax.set_xlabel(r'$x$', fontsize=20)  # Set x label
ax.set_ylabel(r'$p-p_{0}$', fontsize=20)  # Set y label

plt.tight_layout()

fig.savefig('knife_gauge_p.png')

#####################
#####################
