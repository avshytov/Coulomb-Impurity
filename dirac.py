import numpy as np
import scipy
from scipy import linalg
import cmath
import math
from math import *
from pylab import *
import time

from scipy import special

N = 1000
H = np.zeros((2*N,2*N), dtype=complex)
P = np.zeros((2*N,2*N), dtype=complex)
M = np.zeros((2*N,2*N), dtype=complex)
U = np.zeros((2*N,2*N), dtype=complex)
r = np.zeros((N))
pot = np.zeros((N))
psi_up = np.zeros((N))
psi_down = np.zeros((N))
modpsi = np.zeros((N))
j = 1j
m = -1.0
#m = 1.5
rmin = 0.01
rmax = 10.0
t_start = time.time()

for i in range (0,N):
    r[i] = rmin +  i*(rmax-rmin) / N
    pot[i] = 0 #-1/r[i]

for y in range (0,N):
    if y == 0:
        a = r[y+1] - r[y]
        P[0,1] = -1.0 * j / a
        P[1,0] = 1.0 * j /a
    else:
        a = r[y] - r[y-1]
        P[2*y+1,2*y]= 1.0 * j /a
        P[2*y,2*y+1]= -1.0 * j / a
        P[2*y,2*y-1]= 1.0 * j /a
        P[2*y-1,2*y]= -1.0 * j / a
    M[2*y,2*y+1]= -1.0 * j * (m + 0.5) / r[y]
    M[2*y+1,2*y]= 1.0 * j * (m + 0.5) / r[y]
    U[2*y, 2*y] = pot[y]
    U[2*y +1, 2*y +1] = pot[y] 

H = P + M + U
t_fill = time.time()
print "fill matrix: ", t_fill - t_start
print "diagonalising ... "
w, vr =  scipy.linalg.eigh(H)
t_diag = time.time()
print "diag: ", t_diag - t_fill 
print w
hist(w, bins=40)
for i in range (N-5, N+5):
    u = vr[:,i]
    u_up = u[::2]
    u_down = u[1::2]
    u_up_real = u_up.real 
    u_up_imag = u_up.imag
    u_down_real = u_down.real
    u_down_imag = u_down.imag 
    c = norm((np.dot(H,u) - w[i]*u))
    psi_up = u_up_real / sqrt(r)
    psi_down = u_down_imag / sqrt(r)
    modpsi = (psi_up**2 + psi_down**2)

    print 'c = %d' %c 
    figure()
    kr = abs(w[i]) * r
    jm = special.jn(abs(m), kr)
    jm1 = special.jn(abs(m + 1), kr)
          
    C1 = max(abs(u_up)  /sqrt(r)) / max(jm)
    C2 = max(abs(u_down)/sqrt(r)) / max(jm1)
    if (w[i] <  0):
       C1 *= -1
    
    plot(r, (u_up_real / sqrt(r)), label='up r i=%d' %i)
    plot(r, (u_up_imag / sqrt(r)), label='up i i=%d' %i)
    plot(r, (u_down_real / sqrt(r)), label='down r i=%d' %i)
    plot(r, (u_down_imag / sqrt(r)), label='down i i=%d' %i)
    title('Energy state %f' %i)
#     figtext(0.2, 0.85, 'Energy state %d' %n)
#    figure()
    plot(r, modpsi, label='charge dens.')
    plot(r, -C1 * jm,  '--', label='J_m')
    plot(r, -C2 * jm1, '--', label='J_{m + 1}')    
    legend()
show()

