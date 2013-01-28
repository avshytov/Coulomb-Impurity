import numpy as np
import scipy
from scipy import linalg
import cmath
import math
from math import *
from pylab import *
import time
from scipy import special

def diracham(r,pot,mlist):
    N = len(r)
    b = len(mlist)
    H = np.zeros((2*N,2*N), dtype=complex)
    P = np.zeros((2*N,2*N), dtype=complex)
    M = np.zeros((2*N,2*N), dtype=complex)
    U = np.zeros((2*N,2*N), dtype=complex)
    Emat = np.zeros((2*N,b))
    cdtens = np.zeros((b,N,2*N))
    psi_up = np.zeros((N))
    psi_down = np.zeros((N))
    modpsi = np.zeros((N))
    totmodpsi = np.zeros((N))
    j = 1j
   # t_start = time.time()
    for m in range (0,b): 
        print "Calculating Momentum Channel:", mlist[m]
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
        #for m in range (0,b):
         #   print "Calculating Momentum Channel:", 
            M[2*y,2*y+1]= -1.0 * j * (mlist[m] + 0.5) / r[y]
            M[2*y+1,2*y]= 1.0 * j * (mlist[m] + 0.5) / r[y]
            U[2*y, 2*y] = pot[y]
            U[2*y +1, 2*y +1] = pot[y]
            H = P + M + U
             #t_end = time.time()
            # print "fill matrix: ", t_end - t_start
        print "diagonalising... "
        w, vr =  scipy.linalg.eigh(H)
        Emat[:,m] = w[:]
        print "Resolving wavefunctions"
        for i in range (0,2*N):
            u = vr[:,i]
            u_up = u[::2]
            u_down = u[1::2]
            u_up_real = u_up.real
            u_up_imag = u_up.imag
            u_down_real = u_down.real
            u_down_imag = u_down.imag
            #c = norm((np.dot(H,u) - w[i]*u))
            psi_up = u_up_real
            psi_down = u_down_imag
            modpsi = (psi_up**2 + psi_down**2)
            cdtens[m,:,i] = modpsi[:]
            totmodpsi += modpsi
#            if i > (N - 5) and i < (N + 5):
#                print "Plotting state %d" %i
#                figure()
#                plot(r,u_up_real, label='u up real')
#                plot(r,u_up_imag, label='u up imag')
#                plot(r,u_down_real, label='u down real')
#                plot(r,u_down_imag, label='u down imag')
#                legend()
#                figure()
#                plot(r,totmodpsi, label='charge density m %f' %mlist[m])
#                legend()
#                show()
    np.save("cdtens",cdtens)
#    np.save("emat",Emat)
    return Emat, cdtens

def DOS(Emat, mlist ,r):
    N0 = len(Emat)
    c = len(Emat[0])
    N = 1 * N0  #PUT BACK TO 5000
    E = np.zeros((N))
    N0min = 0
    N0max = N0-1
    dos = np.zeros((N))
    dosmat = np.zeros((c,N))
    doschan = np.zeros((N))
    rmax = r[N0/2 -1.0]
    gam = np.pi * 0.7 / rmax    #edit constant
    Emax = 0
    Emin = 0
    wf = np.load("cdtens.npy")
    for k in range (0,c):
        if Emax < Emat[N0max, k]:
            Emax = Emat[N0max, k]
        if Emin > Emat[N0min, k]:
            Emin = Emat[N0min, k]
    for i in range (0,N):
        E[i] = Emin + i * (Emax - Emin)/N
    for m in range (0,c):
        mlab = mlist[m]
        for n in range (N0min, N0max):
            if Emat[n,m] < 0.001 and Emat[n,m] > -0.001:          
                print "Zero energy at n element =  %d" %n
                print "m = %d"  %mlist[m]
                rhoo = cdtens[m,:,n]
                figure()
                plot(r, rhoo, label='charge density energy=%f' %Emat[n,m])
                title('Momentum Channel %d' %mlist[m])
                figtext(0.2, 0.85, 'Energy state %d' %n)
                legend()
               # show()
            for i in range (0,N):
                dosmat[m,i] += (1/np.pi)*gam/(gam**2+(E[i]-Emat[n,m])**2) 
     #plot(E, doschan, label='Density of states for momentum channel %f' %mlab)
    for j in range (0,c):
#        doschan[:] = dosmat[j,:]
#        plot(E,doschan, label='DOS m %d' %mlist[j])
#        legend()
        dos[:] += dosmat[j,:]
    show()
#    dosmat = 2.0 * dosmat  ####### !!!! 
    figure()
    plot(E, (2* dos), label='%d positive momentum channels (doubled), all n' %c) ###!!!! 
    title('Global Density of States')
    legend()
    show()
    np.save("dosmat", dosmat)
    return E, dosmat

if __name__ == '__main__':
   N = 500
   rmin = 0.01
   rmax = 25.0
   r = zeros((N))
   pot = zeros((N))
   a = 10
   mlist = zeros((2*a + 1))
#   mlist = [-1, 0, 1, 2,3,4,5,6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
#   mlist[0] = -2
   for i in range (0,N):
       r[i] = rmin +  i*(rmax-rmin) / N
       pot[i] = 0.0
   for k in range (0, (2*a + 1)):
       mlist[k] = k - a - 1
   print "Momentum Channels:",  mlist
   Emat, cdtens = diracham(r, pot, mlist)
   E, dosmat = DOS(Emat, mlist ,r)
#   LDOS = np.dot(cdtens, dosmat)
#   np.save("LDOS", LDOS)
#   pcolor(r,E,ldos)
#   colorbar()
#   show()





