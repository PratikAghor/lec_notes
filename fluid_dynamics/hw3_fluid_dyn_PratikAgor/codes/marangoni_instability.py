"""
Solving the  Marangoni linear stability
problem, convection driven by gradients in
the surface tension, not bouyancy.

temperature perturbation: theta = thetaHat * exp(i k x + sigma t )
streamfunction perturbation: psi = psiHat * exp(i k x + sigma t )

In the current system: we write thetaHat equations first

The eigenvalue problem is L x = sigma M x,
where x = [thetaHat, psiHat]^{T}

Also, BCs for Marangoni convection:
at z = 0:
    theta = 0  => thetaHat = 0
    d(psi)/dx = 0 => psiHat = 0
    d(psi)/dz = 0 => D (psiHat) = 0

at z = 1:
    d(theta)/dz = 0 => D(thetaHat) = 0
    d(psi)/dx = 0 => psiHat = 0
    d^2(psi)/dz^2 = -const*d(theta)/dx => D2(psiHat) = -LambdaTilde* ik * thetaHat


Author: Pratik Aghor
"""
#####################
import numpy as np
from numpy import sin, cos, tan, pi, arange, zeros, sqrt, exp, eye
import matplotlib
import matplotlib.pyplot as plt
from scipy.linalg import eig

from cheb import *

#####################
"""
Define inputs
"""
Lambda = 1
Q0 = 1
H = 1
kappa = 1e-4
mu = 1e-2

LambdaTilde = 1e-2 # Lambda*Q0*H**2/(kappa*mu)

Ny = 33
n = Ny+1
k = 2.0*pi*arange(1, 5, 1)
Idn = eye(n)

D, zc = cheb(Ny)
k_sigma = np.zeros((len(k), 2))
#####################
"""
convert cheb grid to physical grid
zc = a*z + b
zc = [-1, 1], z= [0, 1]
zc = 2*z - 1

d()/dz = [d()/dzc]*(dzc/dz) ...chain rule

d()/dz = a*D

"""
a = 2
b = -1

D = a*D     # redefine D
D2 = np.matmul(D, D)
D4 = np.matmul(D2, D2)
#####################

for i in range(0, len(k)):
    kx = k[i]
    ksq = kx**2

    L = np.block([ \
    [(D2 - ksq*Idn ), -1j*kx*Idn], \
    [np.zeros((n, n)), (D4 - 2*ksq*D2 + ksq*ksq*Idn) ] \
    ])

    M = np.block([ \
    [Idn, np.zeros((n, n))], \
    [np.zeros((n, n)), np.zeros((n, n)) ]\
    ])

    """
    BCs for Marangoni convection
    first apply BCs at z = 1
    since Cheb grid goes from 1 to -1
    """
    # D(thetaHat) = 0 at z = 1 or zc = 1
    L[0, :] = 0.
    L[0, 0:n] = D[0, :]
    M[0, :] = 0.

    # psiHat = 0 at z = 1 or zc = 1
    L[n, :] = 0.
    L[n, n] = 1.
    M[n, :] = 0.

    # LambdaTilde* ik * thetaHat + D2(psiHat) = 0
    # at z = 1 or zc = 1
    L[n+1, :] = 0.
    L[n+1, 0:n] = LambdaTilde * 1j * kx
    L[n+1, n:2*n] = D2[0, :]
    M[n+1, :] = 0.

    # thetaHat = 0. at at z = 0 or zc = -1
    L[n-1, :] = 0.
    L[n-1, 0] = 1.
    M[n-1, :] = 0.

    # psiHat = 0 at at z = 0 or zc = -1
    L[2*n-1, :] = 0.
    L[2*n-1, n] = 1.
    M[2*n-1, :] = 0

    #  D (psiHat) = 0 at z = 0 or zc = -1
    L[2*n-2, :] = 0.
    L[2*n-1, n:2*n] = D[n-1, :]
    M[2*n-1, :] = 0.
    #############################
    sigma, v = eig(L, b=M, check_finite=True)
    sigma = np.sort(sigma)
    # print("Eigenvalues = ", w)
    k_sigma[i, 0] = kx
    k_sigma[i, 1] = np.real(sigma[0]) # store the largest eigenvalue

#############################
# plot the dispersion relation
# sigma vs k
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()  # Get current axes

ax.plot(k_sigma[:, 0], k_sigma[:, 1], color = 'k', linewidth = 2)  # Plot sigma vs k
ax.set_xlabel(r'$k$', fontsize=20)  # Set x label
ax.set_ylabel(r'$\sigma$', fontsize=20)  # Set y label
plt.tight_layout()

fig.savefig('marangoni_instability_dispersion_reln.png')
#####################
#####################
