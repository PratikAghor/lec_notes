from numpy import pi, cos, arange, ones
import numpy as np

def cheb(N):
    '''Chebyshev polynomial differentiation matrix.
       Ref.: Trefethen's 'Spectral Methods in MATLAB' book.
       author: Pratik Aghor
    '''
    xc  = cos(pi*arange(0,N+1)/N) # xc = Chebyshev-Gauss-Lobatto grid


    D   = np.zeros((N+1, N+1))
    c   = ones(N+1); c[0] = 2.0; c[N] = 2.0

    D[0, 0] = (1.0 + 2.0*N*N)/6.0

    for j in range(0, N+1):
        for k in range(0, N+1):
            if(j == k and j > 0 and j < N):
                D[k, k] = -xc[k]/(2.0*(1.0-xc[k]**2))
            elif(j!=k):
                D[j, k] = (c[j]/c[k])*( ((-1)**(j + k))/(xc[j] - xc[k]) )


    D[N, N] = -(1.0 + 2.0*N*N)/6.0

    return D,xc
