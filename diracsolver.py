import numpy as np
import scipy
from scipy import linalg
import cmath
import math
from math import *
from pylab import *
import time
from scipy import special

def diracham(r,pot):
    N = len(r)
    H = np.zeros((2*N,2*N), dtype=complex)
    P = np.zeros((2*N,2*N), dtype=complex)
    M = np.zeros((2*N,2*N), dtype=complex)
    U = np.zeros((2*N,2*N), dtype=complex)
    psi_up = np.zeros((N))
    psi_down = np.zeros((N))
    modpsi = np.zeros((N))
    totmodpsi = np.zeros((N))
    j = 1j
    m = 1.5
    t_start = time.time()
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
    t_end = time.time()
    print "fill matrix: ", t_end - t_start
    print "diagonalising ... "
    w, vr =  scipy.linalg.eigh(H)
    t_diag = time.time()
    np.savetxt("eigvecs.dat", vr)
    print "diag: ", t_diag - t_end
    for i in range (0, N):
        u = vr[:,i]
        u_up = u[::2]
        u_down = u[1::2]
        u_up_real = u_up.real
        u_up_imag = u_up.imag
        u_down_real = u_down.real
        u_down_imag = u_down.imag
        c = norm((np.dot(H,u) - w[i]*u))
        psi_up = u_up_real# / sqrt(r)
        psi_down = u_down_imag# / sqrt(r)
        modpsi = (psi_up**2 + psi_down**2)
        totmodpsi += modpsi
    plot(r, totmodpsi, label='total charge density?')    
    legend()
    show()
    return w

def DOS(eigvals):
    N = 2000
    N0 = len(eigvals)
    N0min = N0/2 - 100
    N0max = N0/2 + 100
    Emin = eigvals[N0min]
    Emax = eigvals[N0max]
    dos = np.zeros((N))
    E = np.zeros((N))
    gam = np.pi / 1000
    print eigvals
    print Emin, Emax, N0, E[N0/2]
    for i in range (0,N):
        E[i] = Emin + i * (Emax - Emin)/N
    for n in range (N0min, N0max):
        for i in range (0,N):
            dos[i] += (1/np.pi) * gam /(gam**2 + (E[i]-eigvals[n])**2)
    plot(E, dos)
    show()
    return 0

if __name__ == '__main__':
   N = 1000
   rmin = 0.01
   rmax = 25.0
   r = zeros((N))
   pot = zeros((N))
   for i in range (0,N):
       r[i] = rmin +  i*(rmax-rmin) / N
       pot[i] =0.0
   eigvals = diracham(r, pot)
   
   print eigvals
   DOS(eigvals)
