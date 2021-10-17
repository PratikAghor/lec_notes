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
kmax = 2
dk = 0.01
A = 0.5
B = 0.5

k = arange(kmin, kmax + dk, dk)

LArray = [0.3, 0.5, 0.7, 0.9] # Lambda array

#####################
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()  # Get current axes

for i in range(0, len(LArray)):
    Lambda = LArray[i]

    sigma=-(k**4/3)*(1 + 2/((A*(1-Lambda)-2)-B*(k**2)))
    ax.plot(k, sigma, linewidth = 2, label=r'$\Lambda = $'+str(Lambda))  # Plot sigma vs k

ax.plot(k, zeros(len(k)), color = 'k', label=r'$\sigma = 0$')
ax.set_xlabel(r'$k$', fontsize=20)  # Set x label
ax.set_ylabel(r'$\sigma$', fontsize=20)  # Set y label
ax.legend(loc=2)
ax.set_xlim([0, 1])
ax.set_ylim([-5e-3, 5e-3])
plt.tight_layout()

fig.savefig('thin_liq_films_deformable_substrate_dispersion_reln.png')

#####################
"""
Marginal stability curve:
"""
fig = plt.figure(2)  # Create a figure instance
ax = fig.gca()  # Get current axes

ax.plot(k, 1- (B/A)*k**2, color = 'k', label=r'$\Lambda_{N}$')
ax.plot(k, zeros(len(k)), color = 'k', linestyle='--')

ax.set_xlabel(r'$k$', fontsize=20)  # Set x label
ax.set_ylabel(r'$\Lambda$', fontsize=20)  # Set y label
ax.legend(loc=2)
# ax.set_ylim([-5e-3, 5e-3])
plt.tight_layout()

fig.savefig('thin_liq_films_deformable_substrate_marginal_stab.png')

#####################

#####################
