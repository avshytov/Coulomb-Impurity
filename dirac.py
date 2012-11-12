import numpy as np
import scipy
from scipy import linalg
import cmath
import math

N = 1
i = cmath.sqrt(-1.0)
H = np.zeros((2*N,2*N), dtype=complex)
rho = 3.0
m = 1

for y in range (0,2*N):
    for x in range (0,y):
        if (x+y)%2 != 0:
            H[x,y] = -i*((m + 0.5)/rho)
            H[y,x] = scipy.conj(H[x,y])
print H
print scipy.linalg.eigh(H)

#test = np.array([(0, i),(-i, 0)])


#print test
#print scipy.linalg.eigh(test)



