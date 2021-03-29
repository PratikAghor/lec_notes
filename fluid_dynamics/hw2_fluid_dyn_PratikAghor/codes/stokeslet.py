"""
Velocity field due to a single Stokeslet at r=0
psi = -(3/4)r \sin^{2}{\theta}

u ~ [2 cos(theta)/r]e_{r} - [sin(theta)/r]e_{theta} in polar
u ~ [2x/(x^2+y^2)]e_{r} - [y/(x^2+y^2)]e_{theta}
u ~ [(2x^2 + y^2)/(x^2+y^2)^(3/2)]e_{x} + [xy/(x^2+y^2)^(3/2)]e_{theta}

The last expression comes from ux = ur*cos(theta) - utheta*sin(theta)
uy = ur*sin(theta) + utheta*cos(theta),
with tan(theta) = y/x

Author: Pratik Aghor
"""
#####################
import numpy as np
from numpy import sin, cos, tan, arccos, cosh, tanh, arccosh, pi, arange, zeros, sqrt
import matplotlib
import matplotlib.pyplot as plt
#####################
x0 = -10
x1 = 10
y0 = -10
y1 = 10
C= 0.1
nx = 64; ny = 64;

x = np.linspace(x0, x1, nx)
y = np.linspace(y0, y1, ny)

Y, X = np.meshgrid(y, x)

ux = zeros((ny, nx))
uy = zeros((ny, nx))
psi = zeros((ny, nx))
for i in range(0, ny):
    for j in range(0, nx):
        if(x[j]**2 + y[i]**2 > 0.1):
            C = 0.1
        else:
            C = 0.01
        ux[i, j] = C*(2.0*x[j]**2 + y[i]**2)/(x[j]**2 + y[i]**2)**(3/2)
        uy[i, j] = C*x[j]*y[i]/(x[j]**2 + y[i]**2)**(3/2)
        psi[i, j] = -(3.0/4.0)*y[i]**2/(x[j]**2 + y[i]**2)**(1/2)


#####################
skip = (slice(None, None, 8), slice(None, None, 8))
fig = plt.figure()
ax = fig.add_subplot(111, polar=False)
ax.quiver(X[skip], Y[skip], ux[skip], uy[skip])
ax.contour(Y, X, psi, colors='blue')
ax.set_xlabel(r'$x$')  # Set x label
ax.set_ylabel(r'$y$')  # Set y label
plt.savefig('stokeslet.png')
#####################
