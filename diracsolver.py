import numpy as np
import scipy
from scipy import linalg
import cmath
import math
from math import *
from pylab import *
import time
from scipy import special
from ldos import *

def diracham(r,pot,mlist):
    N = len(r)
    b = len(mlist)
    H = np.zeros((2*N,2*N), dtype=complex)
    P = np.zeros((2*N,2*N), dtype=complex)
    M = np.zeros((2*N,2*N), dtype=complex)
    U = np.zeros((2*N,2*N), dtype=complex)
    A = np.zeros((2*N,2*N), dtype=complex)
    Emat = np.zeros((2*N,b))
    cdtens = np.zeros((b,N,2*N))
    psi_up = np.zeros((N))
    psi_down = np.zeros((N))
    modpsi = np.zeros((N))
    totmodpsi = np.zeros((N))
    dr = np.zeros((N))
    dr[1:] = r[1:] - r[0:-1]
    dr[0] = r[1] - r[0]
    np.save("drvec", dr)
    j = 1j
    B = 1
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
            M[2*y,2*y+1]= -1.0 * j * (mlist[m] + 0.5) / r[y]
            M[2*y+1,2*y]= 1.0 * j * (mlist[m] + 0.5) / r[y]
            U[2*y, 2*y] = pot[y]
            U[2*y +1, 2*y +1] = pot[y]
            A[2*y,2*y+1]= -1.0 * j / 2 * B * r[y]#math.sqrt(r[y])
            A[2*y+1,2*y]= 1.0 * j / 2 *B * r[y]#math.sqrt(r[y])
            H = P + M + U - A
        print "diagonalising... "
        w, vr =  scipy.linalg.eigh(H)
        Emat[:,m] = w[:]
        
        if False:
            ea = -2
            eb = 2
            ivals = [t[0] for t in enumerate(w) if t[1] > ea and t[1] < eb]
            evals = [w[t] for t in ivals]
            if mlist[m] == 0:
                ens = list(evals)
            else:
                ens.extend(evals)
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
            modpsi = (abs(u_up)**2 + abs(u_down)**2) / (2 * np.pi * r * dr)
            cdtens[m,:,i] = modpsi[:]
            totmodpsi += modpsi
            if False: # i >= (N - 2) and i <= (N + 1):
                if mlist[m] < 0:
                    print "Plotting state %d" %i
                    figure()
                    plot(r,u_up_real, label='u up real')
                    plot(r,u_up_imag, label='u up imag')
                    plot(r,u_down_real, label='u down real')
                    plot(r,u_down_imag, label='u down imag')
                    legend()
                    title('Momentum Channel %d' %mlist[m])
                    figtext(0.2, 0.85, 'Energy state %d' %i)
                    figtext(0.2, 0.80, 'Energy = %f' %w[i])
                figure()
                plot(r,totmodpsi, label='charge density m %f' %mlist[m])
                legend()
#    hist(ens, bins=40)
    show()
    np.save("cdtens",cdtens)
    np.save("emat",Emat)
#    plot(r,totmodpsi*r, label='charge density m %f' %mlist[m])
#    show()
    return Emat, cdtens

def DOS(Emat, mlist ,r):
    N0 = len(Emat)
    c = len(Emat[0])
    N = 4 * N0
    E = np.zeros((N))
    N0min = 0
    N0max = N0-1
    dos = np.zeros((N))
    dostens = np.zeros((c,N0,N))
    doschan = np.zeros((N))
    rmax = r[N0/2 -1.0]
    gam = np.pi * 0.4 / rmax  
    Emax = 10
    Emin = -Emax
    wf = np.load("cdtens.npy")
    for i in range (0,N):
        E[i] = Emin + float(i) * (Emax - Emin)/N
    for m in range (0,c):
        mlab = mlist[m]
        for n in range (N0min, N0max):
            if Emat[n,m] < 0.001 and Emat[n,m] > -0.001:          
                print "Zero energy at n element =  %d" %n
                print "m = %d"  %mlist[m]
                print "Eigenvalue:", Emat[n,m]
                rhoo = cdtens[m,:,n]
                figure()
                plot(r, rhoo, label='charge density energy=%f' %Emat[n,m])
                title('Momentum Channel %d' %mlist[m])
                figtext(0.2, 0.85, 'Energy state %d' %n)
                legend()
               # show()
            for i in range (0,N): 
                dostens[m,n,i]+=2.0/np.pi*gam**3/(gam**2+(E[i]-Emat[n,m])**2)**2 
    for m in range (0,c):
        for n in range (N0min, N0max):
            dos[:] += dostens[m,n,:]
    dos = 2.0 * dos
    show()
    dostens = 2.0 * dostens  ####### !!!! 
    
    figure()
    plot(E, dos, label='%d momentum channels, 50 n' %c) ###!!!! 
    title('Global Density of States')
    legend()
    show()
    np.save("dostens", dostens)
    np.save("globdos(pos)", dos)
    return E, dostens

if __name__ == '__main__':
   N = 300
   rmin = 0.01
   rmax = 25.0
   r = zeros((N))
   pot = zeros((N))
   a = 3
#   mlist = zeros((2*a + 1))
   mlist = np.array(range(0,a))
#   mlist[0] = 0
   for i in range (0,N):
       r[i] = rmin +  i*(rmax-rmin) / N
#       pot[i] = -1.0 / 1.0 / r[i]
   print "Momentum Channels:",  mlist
   np.save("rvec",r)
   Emat, cdtens = diracham(r, pot, mlist)
   E, dostens = DOS(Emat, mlist ,r)
   np.save("mlist",mlist)
   np.save("Evec",E)
   np.save("potvec",pot)




