import numpy as np
import scipy
from scipy import linalg
import cmath
import math

N = 3
H = np.zeros((2*N,2*N), dtype=complex)
rho = 3.0
m = 1
h = 1
k = 1
j = 1j

a = math.sin(k*h) / h
b =((math.sin(((h * k)/2))))
print a
print b

for x in range (0,N):
    for y in range (0,3):
        H[2*x,2*y]=0
        H[(2*x)+1,(2*y)]= a - j*b
        H[(2*x),(2*y)+1]= a + j*b
  
  
  
print H
#print scipy.linalg.eigh(H)









