import numpy as np
import scipy
from scipy import linalg
import cmath
import math

N = 10
H = np.zeros((2*N,2*N), dtype=complex)
rho = 3.0
m = 1
r = np.zeros((N))
k = 1
j = 1j
rmin = 1
rmax = 10


for i in range (0,N):
    r[i] = rmin * math.exp(math.log(rmax/rmin)/(N - 1.0) * i) 

for x in range (0,N):
    for y in range (0,N):
        if y == 0:
            h = r[y+1] - r[y]
        else:
            h = r[y] - r[y-1]

        a = math.sin(k*h) / h
        b = 2.0 / h * (math.sin(k * h / 2.0)**2)

        H[2*x,2*y]=0
        H[(2*x)+1,(2*y)]= a - j*b
        H[(2*x),(2*y)+1]= a + j*b
          
#print H
print scipy.linalg.eigh(H)









