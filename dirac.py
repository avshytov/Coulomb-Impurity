import numpy as np
import scipy
from scipy import linalg
import cmath
import math
from math import *
from pylab import *

N = 3
H = np.zeros((2*N,2*N), dtype=complex)
r = np.zeros((N))
j = 1j
rmin = 1.0
rmax = 10.0


for i in range (0,N):
    r[i] = rmin +  i*(rmax-rmin) / N

for y in range (0,N):
    if y == 0:
        a = r[y+1] - r[y]
    else:
        a = r[y] - r[y-1]

    if y==0:
        H[0,1] = -1.0 * j / a
        H[1,0] = 1.0 * j /a
    else:
        H[2*y+1,2*y]= 1.0 * j /a
        H[2*y,2*y+1]= -1.0 * j / a
        H[2*y,2*y-1]= 1.0 * j /a
        H[2*y-1,2*y]= -1.0 * j / a
print H
#w, vr =  scipy.linalg.eigh(H)
#for i in range (N-5, N+5):
#    u = vr[i,:]
#    c = norm((np.dot(H,u) - w[i]*u))
#    print i, c
#    plot (r, u, label='i=%d'%i)







