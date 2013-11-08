import math
import cmath
import numpy as np
import scipy
from scipy import special
import matplotlib
from matplotlib import pylab
from pylab import *
import diracsolver
import denscheck
from denscheck import polaris_generic

Us = [-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, -0.05,
       0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]

N = 1000
r0 = 1.0
F = 1.01340014 + 0.05991426 / np.sqrt(N) + 7.12091516 / N #### Correction  
Mmax = 14
Emax0 = 0.0
Emaxs = [1.0]#, 1.3, 1.7, 2.0]
Emin0 = 0.0
Emins = [-1.0]#, -1.3, -1.7, -2.0]
ldosfile0 = np.load('ldos-U0=0-N=%d-Mmax=%d-B=0.npz' %(N,Mmax))
LDOS0 = ldosfile0['ldos']
rbg = ldosfile0['r']
gambg = ldosfile0['gam']
Ev0 = ldosfile0['Ev']
rho0s = np.zeros((len(rbg),len(Emins),2))
RPA = np.zeros((len(rbg),len(Emaxs)))
totalrho = np.zeros((len(rbg),len(Emins),2))
rhoband = np.zeros((len(rbg),len(Emins),2))
rhosim = np.zeros((len(rbg),len(Emins),2))
collection = []
Q = np.zeros((len(Us),2))

def Uq_Coulomb(q):
    return 2.0 * np.pi / q * np.exp(-q*r0)

def Uvals(Beta, rvals):
    return -Beta / np.sqrt(rvals**2 + r0**2)

def highm(Ef, r, Mmax, U0):
    Jsum = np.zeros((len(r)))
    Efr = Ef + Uvals(U0,r)/2.0
    for i in range (-Mmax, Mmax+1):
        Jsum += 2.0 * (special.jn(i,(Ef*r)))**2
    drhohm = 2.0 - Jsum
    drhohm += (special.jn(-Mmax,(Ef*r)))**2
    drhohm -= (special.jn(1+Mmax,(Ef*r)))**2
    drhohm *= -abs(Ef) / 4.0 / np.pi * Uvals(U0,r)
    return drhohm

for i in range (0,len(Emins)):
    rho0s[:,i,0] = diracsolver.find_rho(Ev0,rbg,LDOS0,Emins[i],Emax0)
    rho0s[:,i,1] = diracsolver.find_rho(Ev0,rbg,LDOS0,Emin0,Emaxs[i])
    RPA0 = - 1.0 / 16.0 * Uvals(1.0,rbg)**3 
    RPA[:,i] = polaris_generic(rbg,Emins[i],Uq_Coulomb)


for t in range (0,len(Us)):
    x=Us[t]
    figure()
    title('drho = %g' %x)
    ldosfile = np.load('ldos-U0=%g-N=%d-Mmax=%d-B=0.npz' %(x,N,Mmax))
    mlist = ldosfile['mlist']
    gam = ldosfile['gam']
    r = ldosfile['r']
    Ev = ldosfile['Ev']
    LDOS = ldosfile['ldos']
    rhoTFe = Uvals(x,r)**2 / 4.0 / np.pi 
    if abs(len(rbg)-len(r)) != 0 or abs(gam-gambg) > 1e-6:
        print 'MISMATCHING LDOS'
    test = np.load('allrhos-U0=%g-N=%d-Mmax=%d-B=0-Emin=-1-Emax=-0.npy' %(x,N,Mmax))
    for i in range (0,len(Emins)):
        rhosim[:,i,0] = diracsolver.find_rho(Ev,r,LDOS,Emins[i],Emax0)
        rhosim[:,i,1] = diracsolver.find_rho(Ev,r,LDOS,Emin0,Emaxs[i])
        rhoband[:,i,0] = F*(rhosim[:,i,0] - rho0s[:,i,0])
        rhoband[:,i,1] = F*(rhosim[:,i,1] - rho0s[:,i,1])
        totalrho[:,i,0] = rhoband[:,i,0] + RPA[:,i]*x - rhoTFe
        totalrho[:,i,1] = rhoband[:,i,1] - RPA[:,i]*x - rhoTFe
        totalrho[:,i,0] += highm(Emax0, r, Mmax, x)
        totalrho[:,i,0] -= highm(Emins[i], r, Mmax, x)
        totalrho[:,i,1] += highm(Emaxs[i], r, Mmax, x)
        totalrho[:,i,1] -= highm(Emin0, r, Mmax, x)
        plot(r,totalrho[:,i,0], label='%g-%g' %(Emins[i], Emax0))                        
        plot(r,-totalrho[:,i,1], label='%g-%g' %(Emin0, Emaxs[i]))                 
    plot(r,x*RPA0,'k--', label='RPA')         
    legend()

    rcut = 20
    imax = 400
    print "Ef step at r=", r[rcut]


 ## Low r values
    Qh = 0.5 * totalrho[0,0,0] * r[0]
    Qe = 0.5 * totalrho[0,0,1] * r[0]
    for k in range (1,rcut):
        Qh += totalrho[k,0,0]*r[k]
        Qe += totalrho[k,0,1]*r[k]
 ## High r values
 ## Change window if we want to cut
    Qh += 0.5 * totalrho[-1,0,0] * r[imax-1]
    Qe += 0.5 * totalrho[-1,0,1] * r[imax-1]
    for k in range (rcut,imax-1):
        Qh += totalrho[k,0,0]*r[k]
        Qe += totalrho[k,0,1]*r[k]

    Q[t,0] = Qh*(r[1]-r[0])*2.0*np.pi
    Q[t,1] = Qe*(r[1]-r[0])*2.0*np.pi
    Us = np.array(Us)
#Qt = (( 1.0 / np.sqrt(r0**2+r[0]**2))
#      - (1.0/ np.sqrt(r0**2 + r[N-1]**2)))
#Qt *= np.pi / 8.0 * r0
Qt = np.pi / 8.0

alp = Us #/ 4.0 / np.pi

figure()
title('Total Charge vs Potential Strength rmax=%g'%r[imax-1])
plot(Us,Q[:,0], label='Holes')
plot(Us,Q[:,1], label='Electrons')
plot(-Us,Q[:,1], '--', label='Elec vs -U0')
plot(-Us,-Q[:,0], '--', label='-Holes vs -U0')
plot(Us, (np.pi/8.0*alp + 0.19*alp**3), label='B-S')
plot(Us, Us*Qt, label='RPA Theory')
legend()

show()
