import math
import cmath
import numpy as np
import scipy
from scipy import special
import matplotlib
from matplotlib import pylab
from pylab import *

Us = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
Emin = -1.0
Emax = -0.0
mlist = np.arange(15)
r0 = 1.0
r = np.linspace(0.01,25,1000)
comprhos = np.zeros((len(r),len(Us)))
theorhos = np.zeros((len(r),len(Us)))

def highm(Ef, r, Mmax, U0):
    Jsum = np.zeros((len(r)))
    #Efr = Ef + Uvals(U0,r)/2.0                                                           
    for i in range (-Mmax, Mmax+1):
        Jsum += 2.0 * (special.jn(i,(Ef*r)))**2
    drhohm = 2.0 - Jsum
    drhohm += (special.jn(-Mmax,(Ef*r)))**2
    drhohm -= (special.jn(1+Mmax,(Ef*r)))**2
    drhohm *= - abs(Ef) / 4.0 / np.pi * Uvals(U0,r, r0)
    return drhohm

def Uvals(Beta, rvals, r0):
    return -Beta /np.sqrt(rvals**2 + r0**2)

def RPAvals(beta, rvals, r0):
    return beta / 16.0 * (-Uvals(1.0,rvals, r0))**3

for i in range (0, len(Us)):
    U0 = Us[i]
    print 'Calculating U0 = ', U0
    A = np.load('allrhos-U0=%g-N=1000-Mmax=%d-B=0-Emin=%g-Emax=%g.npy' 
                %(U0,np.max(mlist),Emin,Emax))
    if norm(r-A[0,:]) > 1e-6:
        print 'ERROR r GRID MISMATCH!!! U0 =', U0
    rhotot = A[1,:]+A[8,:]-Uvals(U0,r,r0)**2/4.0/np.pi
    rhotot += highm(Emax,A[0,:],np.max(mlist),U0)
    rhotot -= highm(Emin,A[0,:],np.max(mlist),U0)
    comprhos[:,i] = rhotot[:]
    if abs(U0) >= 0.5:
        gamma = np.sqrt(U0**2-0.25)
        sgn = U0 / abs(U0)
        theorhos[:,i] = gamma*sgn/2.0/(np.pi**2)*Uvals(1.0,r,r0)**2
        theorhos[:,i] += RPAvals(sgn*0.5,r,r0)

if True:
    for i in range (0, len(Us)):
        U0 = Us[i]
        figure()
        title('Total drho, U0 = %g' %U0)
        plot(r,comprhos[:,i], label='sim')
        plot(r,-comprhos[:,(len(Us)-1-i)], '--', label='Opposite U0')
        plot(r,theorhos[:,i], label='Theory')
        legend()

if True:
    figure()
    for i in range (0, len(Us)):
      #  figure()
        title('loglog Total drho')
        loglog(r, abs(comprhos[:,i]), label='%g' %Us[i])
    loglog(r, -1.0/2.0/np.pi**2*Uvals(1.0,r,r0),'c--', label='1/r')
    loglog(r, 1.0/2.0/np.pi**2*Uvals(1.0,r,r0)**2,'k--', label='1/r^2')
    loglog(r, -1.0/2.0/np.pi**2*Uvals(1.0,r,r0)**3, 'r--', label='1/r^3')
       # loglog(r,abs(theorhos[:,i]), label='Theory')
    legend()

if False:
    for i in range (0, len(Us)):
        figure()
        title('Difference of Sim - Theory U0 = %g' %Us[i])
        diff = comprhos[:,i]-theorhos[:,i]
        plot(r, diff)
#        figure()
#        title('Difference of Sim - Theory U0 = %g loglog' %Us[i])
#        loglog(r,abs(diff))
show()
