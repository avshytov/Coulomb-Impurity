import numpy as np
import scipy
from scipy import linalg
import cmath
import math
from math import *
from pylab import *
from scipy import special
from ldos import *
import time


def diracham(r,pot,mlist,B0):
    N = len(r)
    b = len(mlist)
    H = np.zeros((2*N,2*N), dtype=complex)
    P = np.zeros((2*N,2*N), dtype=complex)
    M = np.zeros((2*N,2*N), dtype=complex)
    U = np.zeros((2*N,2*N), dtype=complex)
    Ematp = np.zeros((2*N,b))
    Ematn = Ematp
    cdtens = np.zeros((b,N,2*N))
    cdtensp = cdtens
    cdtensn = cdtens
    psi_up = np.zeros((N))
    psi_down = np.zeros((N))
    modpsi = np.zeros((N))
    totmodpsi = np.zeros((N))
    dr = np.zeros((N))
    dr[1:] = r[1:] - r[0:-1]
    dr[0] = r[1] - r[0]
    np.save("drvec", dr)
    j = 1j
    timestart1 = time.time() ####

    for B in [B0, -B0]:
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
                M[2*y,2*y+1]= -1.0*j*(mlist[m]+(B*r[y]**2/2)+0.5)/r[y]
                M[2*y+1,2*y]= 1.0* j*(mlist[m]+(B*r[y]**2/2)+0.5)/r[y]
                U[2*y, 2*y] = pot[y]
                U[2*y +1, 2*y +1] = pot[y]
                H = P + M + U
            print "diagonalising... "
            w, vr =  scipy.linalg.eigh(H)
            if B == B0:
                Ematp[:,m] = w[:]
            elif B == (-1.0 * B0) and B!= 0.0:
                Ematn[:,m] = w[:]
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
            iplot = []
            for Eplot in []:#[-0.24, -0.04]:
                Ei = list(enumerate (w))
                Ei.sort (lambda x, y: cmp(abs(x[1]- Eplot), abs(y[1]- Eplot)))
                iplot.append(Ei[0][0])
            print "iplot:", iplot
            for i in range (0,2*N):
                u = vr[:,i]
                u_up = u[::2]
                u_down = u[1::2]
                u_up_real = u_up.real
                u_up_imag = u_up.imag
                u_down_real = u_down.real
                u_down_imag = u_down.imag
                modpsi =(abs(u_up)**2+abs(u_down)**2)/(2*np.pi*r*dr)
                if B == B0:
                    cdtensp[m,:,i] = modpsi[:]
                elif B == (-1.0 * B0) and B!=0.0:
                    cdtensn[m,:,i] = modpsi[:]
                totmodpsi += modpsi
                if i in iplot:
                #if  i >= (N - 2) and i <= (N + 1):
                    if False: #mlist[m] < 0:
                        print "Plotting state %d" %i
                        figure()
                        plot(r,u_up_real, label='u up real')
                        plot(r,u_up_imag, label='u up imag')
                        plot(r,u_down_real, label='u down real')
                        plot(r,u_down_imag, label='u down imag')
                        legend()
                        title('Momentum Channel %d' %mlist[m])
                    figure()
                    plot(r,totmodpsi,label='charge density m %f' %mlist[m])
                    legend()
        timeend1 = time.time()
        print "Time taken:", timeend1 - timestart1
    cdtens = cdtensp + cdtensn
    np.save("cdtens",cdtens)
    return Ematp, Ematn, cdtens

def DOS(Ematp, Ematn, mlist ,r):
    N0 = len(Ematp)
    c = len(mlist)
    N = 4 * N0
    E = np.zeros((N))
    N0min = 0
    N0max = N0-1
    dos = np.zeros((N))
    dostens = np.zeros((c,N0,N))
    doschan = np.zeros((N))
    rmax = r[N0/2 -1.0]
    gam = np.pi * 0.8 / rmax  
    Emax = 24.0
    Emin = -Emax
    wf = np.load("cdtens.npy")
    timestart2 = time.time()
    A = 2.0/np.pi*gam**3 
    for i in range (0,N):
        E[i] = Emin + float(i) * (Emax - Emin)/N
    for m in range (0,c):
        mlab = mlist[m]
        for n in range (N0min, N0max):
            if False: # Emat[n,m] < 0.001 and Emat[n,m] > -0.001:          
                print "Zero energy at n element =  %d" %n
                print "m = %d"  %mlist[m]
                print "Eigenvalue:", Emat[n,m]
                rhoo = cdtens[m,:,n]
                figure()
                plot(r, rhoo, label='charge density energy=%f' %Emat[n,m])
                title('Momentum Channel %d' %mlist[m])
                figtext(0.2, 0.85, 'Energy state %d' %n)
                legend()
                show()
            for i in range (0,N): 
                dostens[m,n,i]+=A/(gam**2+(E[i]-Ematp[n,m])**2)**2
                dostens[m,n,i]+=A/(gam**2+(E[i]-Ematp[n,m])**2)**2
 
    for m in range (0,c):
        for n in range (N0min, N0max):
            dos[:] += dostens[m,n,:]
    np.save("dostens", dostens)
    np.save("globdos(pos)", dos)
    timeend2 = time.time()
    print "Total time:", timeend2 - timestart2 
    return E, dostens

if __name__ == '__main__':
   N = 200
   rmin = 0.01
   rmax = 25.0
   B0 = 0.0
   r = zeros((N))
   pot = zeros((N))
   a = 1
   Ustr = 0.9
   info = np.zeros((2))
   info[0] = Ustr
   info[1] = B0
#   mlist = zeros((2*a + 1))
   mlist = np.array(range(0,a))
#   mlist[0] = 0
   for i in range (0,N):
       r[i] = rmin +  i*(rmax-rmin) / N
       pot[i] = -Ustr / r[i]
   print "Momentum Channels:",  mlist
   np.save("rvec",r)
   Ematp, Ematn, cdtens = diracham(r, pot, mlist,B0)
   E, dostens = DOS(Ematp, Ematn, mlist ,r)
   np.save("mlist",mlist)
   np.save("Evec",E)
   np.save("potvec",pot)
   np.save("EMinfo", info)



