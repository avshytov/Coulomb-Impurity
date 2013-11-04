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
def Uq_Coulomb(q):
    return 2.0 * np.pi / q * np.exp(-q*r0)

Us = [-0.5, -0.4, -0.3, -0.2, -0.1, -0.05, 
       0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
N = 1000
r0 = 1.0
F = 1.01340014 + 0.05991426 / np.sqrt(N) + 7.12091516 / N #### Correction  
Mmax = 14
Emax0 = 0.0
Emaxs = [1.0, 1.4, 1.6, 1.8]
Emin0 = 0.0
Emins = -np.array(Emaxs)
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
Q = np.zeros((len(Us),len(Emaxs),2))

def highm(Ef, r, Mmax, U0):
    ###### r0 is set as 1.0 ###### 
    Jsum = np.zeros((len(r)))
    #Efr = Ef + Uvals(U0,r)/2.0
    for i in range (-Mmax, Mmax+1):
        Jsum += 2.0 * (special.jn(i,(Ef*r)))**2
    drhohm = 2.0 - Jsum
    drhohm += (special.jn(-Mmax,(Ef*r)))**2
    drhohm -= (special.jn(1+Mmax,(Ef*r)))**2
    drhohm *= abs(Ef) / 4.0 / np.pi * U0 / np.sqrt(r**2 + 1.0)  
    return drhohm

for i in range (0,len(Emins)):
    rho0s[:,i,0] = diracsolver.find_rho(Ev0,rbg,LDOS0,Emins[i],Emax0)
    rho0s[:,i,1] = diracsolver.find_rho(Ev0,rbg,LDOS0,Emin0,Emaxs[i])
    RPA0 = 1.0 / 16.0 / (np.sqrt(rbg**2 + r0**2))**3 
    RPA[:,i] = polaris_generic(rbg,Emaxs[i],Uq_Coulomb)


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
    rhoTFe = x**2 / 4.0 / np.pi / (r**2 + r0**2)
    if abs(len(rbg)-len(r)) != 0 or abs(gam-gambg) > 1e-6:
        print 'MISMATCHING LDOS'
    for i in range (0,len(Emins)):
        #figure()
        #title('drho = %g' %x)
        
        rhosim[:,i,0] = diracsolver.find_rho(Ev,r,LDOS,Emins[i],Emax0)
        rhosim[:,i,1] = diracsolver.find_rho(Ev,r,LDOS,Emin0,Emaxs[i])
        rhoband[:,i,0] = rhosim[:,i,0] - rho0s[:,i,0]
        rhoband[:,i,1] = rhosim[:,i,1] - rho0s[:,i,1]
        #collection.append(rhos)
        totalrho[:,i,0] = F*rhoband[:,i,0] + RPA[:,i]*x - rhoTFe
        totalrho[:,i,1] = F*rhoband[:,i,1] - RPA[:,i]*x - rhoTFe
        totalrho[:,i,0] += highm(Emax0, r, Mmax, x)
        totalrho[:,i,0] -= highm(Emins[i], r, Mmax, x)
        totalrho[:,i,1] += highm(Emaxs[i], r, Mmax, x)
        totalrho[:,i,1] += highm(Emin0, r, Mmax, x)
        plot(r,totalrho[:,i,0], label='%g-%g' %(Emins[i], Emax0))                        
        plot(r,-totalrho[:,i,1], '--', label='%g-%g' %(Emin0, Emaxs[i]))                 
    plot(r,x*RPA0,'k*', label='RPA')         

#        kmax = len(r)/2
#        Qh = 0.5*(totalrho[0,i,0]*r[0]+totalrho[kmax,i,0]*r[kmax])
#        Qe = 0.5*(totalrho[0,i,1]*r[0]+totalrho[kmax,i,1]*r[kmax])
#        for k in range (1,kmax):
#            Qh += totalrho[k,i,0]*r[k]
#            Qe += totalrho[k,i,1]*r[k]
#        Q[t,i,0] = Qh*(r[1]-r[0])*2.0*np.pi
#        Q[t,i,1] = Qe*(r[1]-r[0])*2.0*np.pi
#        Us = np.array(Us)
#    title('Total Charge against Potential Strength')
#    plot(Us,Q[:,i,0], label='%g, %g' %(Emins[i],Emax0))
#    plot(Us,Q[:,i,1], label='%g, %g' %(Emin0,Emaxs[i]))
#    plot(Us, (np.pi/8.0*Us + 0.19*Us**3), label='B-S')


    legend()

show()
