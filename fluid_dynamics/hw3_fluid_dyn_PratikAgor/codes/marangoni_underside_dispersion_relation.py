"""
Plotting the dispersion relation for
Linear stability of a uniformly thick film
lining the underside of a rigid flat horizontal substrate
sigma = (B/3 - C\Lambda/2)k^{2} - (1/3)k^{4}
Author: Pratik Aghor
"""
#####################
import numpy as np
from numpy import sin, cos, tan, pi, arange, zeros, sqrt, exp
import matplotlib
import matplotlib.pyplot as plt
#####################
kmin = 0
kmax = 0.4
dk = 0.01

k = arange(kmin, kmax + dk, dk)

aArray = [-0.1, 0., 0.1] # alias for (2B - 3C\Lambda)

fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()  # Get current axes

for i in range(0, len(aArray)):
    a = aArray[i]
    sigma = a*k**2 - k**4
    ax.plot(k, sigma, linewidth = 2, label=r'$a = $'+str(a))  # Plot sigma vs k

ax.plot(k, zeros(len(k)), color = 'k', label=r'$\sigma = 0$')
ax.set_xlabel(r'$k$', fontsize=20)  # Set x label
ax.set_ylabel(r'$\sigma$', fontsize=20)  # Set y label
ax.legend(loc=2)
ax.set_ylim([-3e-3, 3e-3])
plt.tight_layout()

fig.savefig('marangoni_dispersion_reln.png')

#####################
#####################
