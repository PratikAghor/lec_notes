"""
Plot steady state film thickness
for a pendant droplet on an inclined substrate


Author: Pratik Aghor
"""
#####################
import numpy as np
from numpy import sin, cos, tan, pi, arange, zeros, sqrt, exp
import matplotlib
import matplotlib.pyplot as plt
#####################
# steday state h as a function of the bond number
# and other parameters
def get_h(B, alpha, x, V0):
    c1 = (-24.*V0 + B*sin(alpha))/2.
    c2 = (72.*V0 - B*sin(alpha))/12.
    h = (-B*sin(alpha))*((x**3)/6.) + c1* ((x**2)/2.) + c2*x
    return h
#####################
dx = 0.01
x = arange(0, 1 + dx, dx)
alpha = pi/6.
V0Array = np.array([1., 50, 100])
BArray = np.array([1, 50, 100])
nB = len(BArray)
nV0 = len(V0Array)
#####################
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()  # Get current axes
for i in range(0, nB):
    B = BArray[i]
    V0 = V0Array[0] # fix V0
    h = get_h(B, alpha, x, V0)
    ax.plot(x, h, linewidth = 2, label=r'$B = $'+str(B))  # Plot h vs x
    ax.set_xlabel(r'$x$', fontsize=20)  # Set x label
    ax.set_ylabel(r'$h$', fontsize=20)  # Set y label

ax.plot(x, zeros(len(x)), linewidth = 2, linestyle = '--', color = 'k', label=r'$x = 0$')  # Plot h vs x

ax.legend()
plt.tight_layout()

fig.savefig('pendant_droplet_plot_h_fix_V0_vary_B.png')

#####################
fig = plt.figure(2)  # Create a figure instance
ax = fig.gca()  # Get current axes
for i in range(0, nB):
    B = BArray[0] # fix B
    V0 = V0Array[i]
    h = get_h(B, alpha, x, V0)
    ax.plot(x, h, linewidth = 2, label=r'$V0 = $'+str(V0))  # Plot h vs x
    ax.set_xlabel(r'$x$', fontsize=20)  # Set x label
    ax.set_ylabel(r'$h$', fontsize=20)  # Set y label


ax.plot(x, zeros(len(x)), linewidth = 2, linestyle = '--', color = 'k', label=r'$x = 0$')  # Plot h vs x
ax.legend()
plt.tight_layout()

fig.savefig('pendant_droplet_plot_h_fix_B_vary_V0.png')


#####################
