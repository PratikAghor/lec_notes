"""
Oseen flow past a sphere

\psi = (r^2/2 + 1/4r)sin^2(theta) - (3/(2Re))(1+c)(1-exp(-Re r (1-c)/2 ) )
=>
\psi = ((x^2+y^2)/2 + 1/4((x^2+y^2)^1/2)))*(y/(x^2+y^2)) - (3/(2Re))(1+c)(1-exp(-Re (x^2+y^2)^1/2 (1-c)/2 ) )
with tan(theta) = y/x

Author: Pratik Aghor
"""
#####################
import numpy as np
from numpy import sin, cos, tan, pi, arange, zeros, sqrt, exp
import matplotlib
import matplotlib.pyplot as plt
#####################
x0 = -5
x1 = 5
y0 = -5
y1 = 5

nx = 128; ny = 128;
Re = 1
x = np.linspace(x0, x1, nx)
y = np.linspace(y0, y1, ny)

Y, X = np.meshgrid(y, x)

psi = zeros((ny, nx))
for i in range(0, ny):
    for j in range(0, nx):
        r = sqrt(x[j]**2 + y[i]**2)
        if (r>=1):
            sintheta = y[i]/sqrt(x[j]**2 + y[i]**2)
            c = x[j]/sqrt(x[j]**2 + y[i]**2)
            psi[i, j] = (r**2/2.0 + 1.0/(4.0*r))*sintheta**2 - (3.0/(2.0*Re))*(1.0+c)*(1.0- exp(-Re*r*(1.0-c)/2.0 ))


#####################
levels = np.arange(0.1, 5, 0.5)
fig = plt.figure()
ax = fig.add_subplot(111)
circle1 = plt.Circle(( 0 , 0 ), 1, color='k')
ax.contour(Y, X, psi, colors='blue', levels= levels)
ax.add_patch(circle1)
plt.axis('equal')
plt.ylim(-4., 4)
plt.xlim(-1, 1)
ax.set_xlabel(r'$x$')  # Set x label
ax.set_ylabel(r'$y$')  # Set y label
plt.savefig('oseen_Re_'+str(Re)+'.png')
#####################
