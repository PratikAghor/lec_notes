"""
Plot the laminar jet profile
Author: Pratik Aghor
"""
#####################
import numpy as np
from numpy import sin, cos, tan, pi, arange, zeros, sqrt, exp, cosh
import matplotlib
import matplotlib.pyplot as plt
#####################
ymin = -5
ymax = 5
dy = 0.1
y = arange(ymin, ymax+dy, dy)

u = 1.0/(cosh(y)*cosh(y))

fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()  # Get current axes
ax.plot(u, y, linewidth = 2, color = 'r', label=r'$u$')  # Plot sigma vs k
ax.set_xlabel(r'$u$', fontsize=20)  # Set x label
ax.set_ylabel(r'$y$', fontsize=20)  # Set y label

plt.tight_layout()

fig.savefig('laminar_jet_profile.png')

#####################
#####################
