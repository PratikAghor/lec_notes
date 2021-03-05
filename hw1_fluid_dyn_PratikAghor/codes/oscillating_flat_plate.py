"""
This code plots the solution of the oscillating flat plate problem
at 8 different times during one cycle.

Also known as the Stokes' second problem.
u(y, t) = U_{0} \exp{\left(- \sqrt\frac{\Omega}{2\nu}  y \right)} \sin{\left(- \sqrt\frac{\Omega}{2\nu}  y  + \Omega t\right)

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
Omega = 1.
nu = 1.
U0 = 0.1

freq = sqrt(nu/(2.0*Omega))
###########################################
"""
define y = [0, 7], t = [0, 2pi/Omega]
"""
Ny = 100
Nt = 9
y = np.linspace(0.0, 7.0, num=Ny)
t = np.linspace(0.0, 2.0*pi/Omega, num=Nt)
print("y = ", y)
print("t = ", t)
###########################################
"""
define u(y, t)
"""
u = np.zeros((Nt, Ny))

for i in range(0, Nt):
    for j in range(0, Ny):
        u[i, j] = exp(-freq * y[j])* sin(-freq * y[j] + Omega * t[i])

u *= U0
###########################################
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()  # Get current axes
for i in range(0, Nt-1):
    ax.plot(u[i, :], y, label = r"$\frac{t}{T} = $ " + str(i))  # Plot U
    ax.set_xlabel(r'$u$', fontsize=20)  # Set x label
    ax.set_ylabel(r'$y$', fontsize=20)  # Set y label
    ax.legend(loc=1)

fig.savefig('u_vs_y_diff_t.png')
###########################################
