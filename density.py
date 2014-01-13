import numpy as np
import util
from numpy import linalg
import math
import coulomb
import scipy
from scipy import special
from scipy.interpolate import splev, splrep
from diracsolver import (makeDirac, solveDirac, dos_bg, diracLDOS, find_rho, getDensity,
                         prepareLDOS, bandDensity)
import mkrpa2
import rpakernel
from odedirac import (_sgn, odedos_m, doscalc, rhocalc) 

#def backgroundDensity(r, mlist, Temp):
#    rho_0 = bandDensity(r, 0.0 * r, mlist, BTemp)
#    return rho_0
    
def RPA_kernel(r, kF):
    Q_rpa = mkrpa2.RPA_inter(r)  
    if (kF * r.max() > 0.01):
        Q_rpa += mkrpa2.RPA_intra(r, kF)
    return Q_rpa

class GrapheneResponse:
    def __init__ (self, r, Ef, **kwargs):
        params = {
           'Temp'   : 0.01, 
            'Ecut'  : -1.5, 
            'B'     : 0.0, 
            'Mmax'  : 16
        }
        params.update(**kwargs)
        self.Emin = min(Ef, params['Ecut']) # Ecut could be both positive
        self.Emax = max(Ef, params['Ecut']) # and negative
        self.Temp = params['Temp']
        self.r    = r
        self.rexp = util.make_exp_grid(r.min(), r.max(), len(r))
        self.B = params['B']
        self.mlist =  np.array(range(0, params['Mmax']))
        
        if False:  # Use this to skip kernel calculation; helpful in debugging
           self.Q_Emin = np.zeros((len (self.rexp), len(self.rexp)))#RPA_kernel(self.rexp, abs(self.Emin))
           self.Q_Emax = np.zeros(np.shape(self.Q_Emin)) #RPA_kernel(self.rexp, abs(self.Emax))
           self.Qm_Emax = np.zeros(np.shape(self.Q_Emin))#rpakernel.kernel_m(self.rexp, self.mlist, abs(self.Emax))
           self.Qm_Emin = np.zeros(np.shape(self.Q_Emax))#rpakernel.kernel_m(self.rexp, self.mlist, abs(self.Emin))
        else:
           # Full RPA kernel
           self.Q_Emin  = RPA_kernel(self.rexp, abs(self.Emin))
           self.Q_Emax  = RPA_kernel(self.rexp, abs(self.Emax))
           # m-resolved RPA kernel
           self.Qm_Emax = rpakernel.kernel_m(self.rexp, self.mlist, abs(self.Emax))
           self.Qm_Emin = rpakernel.kernel_m(self.rexp, self.mlist, abs(self.Emin))
        self.rho_0 = self.diracDensity(np.zeros(np.shape(r)))
        N = len(self.r)
          
    if False:
        def diracDensity(self, U):
            return bandDensity(self.r, U, self.mlist, self.B, 
                               self.Emin, self.Emax, self.Temp)
    
    if True:
        def diracDensity(self, U):
            return rhocalc(self.Emin, self.Emax, self.r, 
                           U, self.mlist)
                       
    def bandResponse(self, U):
        rho_b = self.diracDensity(U)
        return (rho_b - self.rho_0) * self.F
        
    def apply_kernel(self, Q, U):
        U_exp = util.gridswap(self.r, self.rexp, U)
        rho_exp = np.dot(Q, U_exp)
        return util.gridswap(self.rexp, self.r, rho_exp)
    
    def rho_RPA(self, U):
        return self.apply_kernel(self.Q_Emax, U)
    
    def rho_RPA_m(self, U):
        return self.apply_kernel(self.Qm_Emax - self.Qm_Emin, U)
    
    def highm(self, U): # Calculation of high-m correction within RPA
        Qmax =  self.Q_Emax - self.Qm_Emax
        Qmin =  self.Q_Emin - self.Qm_Emin
        return -self.apply_kernel(Qmax - Qmin, U)
        
    def highm_old(self, E, U): # old way to calculate high-m correction
        Mmax = self.mlist[-1]
        Jsum = np.zeros((len(self.r)))
        Efr = E + U/2.0
        for i in range (-Mmax, Mmax+1):
          Jsum += special.jn(i, E*self.r)**2 + special.jn(i + 1, E*self.r)**2
          #Jsum += 2.0 * (special.jn(i,(Ef*r)))**2
        drhohm = 2.0 - Jsum
        #drhohm += (special.jn(-Mmax,(Ef*r)))**2
        #drhohm -= (special.jn(1+Mmax,(Ef*r)))**2
        drhohm *= -abs(E) / 4.0 / np.pi * U
        return drhohm

    def highm2(self, U): # calculate quadratic part of high-m correction
        eps = 1e-5 
        highm_max1 = self.highm_old(self.Emax + eps, U)
        highm_max2 = self.highm_old(self.Emax - eps, U)
        highm_min1 = self.highm_old(self.Emin + eps, U)
        highm_min2 = self.highm_old(self.Emin - eps, U)
        dnu_max = (highm_max1 - highm_max2) / 2.0 / eps * _sgn(self.Emax)
        dnu_min = (highm_min1 - highm_min2) / 2.0 / eps * _sgn(self.Emin)
        return (dnu_max - dnu_min) * U / 2.0;
    
    
    def seaContribution(self, U):
        sgnE = 1.0
        if (self.Emin < 0): sgnE = -1.0
        rho = np.zeros(np.shape(U))
        rho = self.highm(U) * sgnE + self.highm2(U)
        #rho  = self.highm( self.Emax,  U)
        #rho -= self.highm( self.Emin,  U)
        
        rho += U**2 / 4.0 / np.pi * sgnE # Quadratic contribution
        
        rho += self.apply_kernel(self.Q_Emin, U) 
        return rho
        
    def rho_U(self, U):    
        rho  = self.bandResponse(U)
        rho += self.seaContribution(U) 
        Rmax = 30.0
        imax = np.abs(self.r - Rmax).argmin()
        rho_rpa = self.apply_kernel(self.Q_Emax, U)
        rho[imax:] = rho_rpa[imax:]
        return rho 


if __name__ == '__main__':
   rmin = 0.01
   rmax = 50.0
   N = 500
   r = util.make_lin_grid(rmin, rmax, N) 
   Ef = -0.2
   Ecut = -3.0
   graphene = GrapheneResponse(r, -1e-4, Ecut=Ecut, Mmax=31)
   
   Z = 0.25
   r_0 = 1.0
   U = Z / np.sqrt(r**2 + r_0**2)
   rho_th = -Z * r_0 / 16.0 / np.sqrt(r**2 + r_0**2)**3
   
   rho_RPA = graphene.rho_RPA(U)
   rho_U   = graphene.rho_U(U)
   rho_Ub  = graphene.bandResponse(U) 
   rho_RPAb = graphene.apply_kernel(graphene.Q_Emax - graphene.Q_Emin, U)
   rho_0 = (graphene.Emax**2 - graphene.Emin**2) / 4.0 / math.pi
   rho_1 = graphene.diracDensity(U)   

   if True:
       imin = 0
       imax = 199
       rmin = r[0]
       rmax = r[-1]
       dr = r[imin+1]-r[imin]
       Ugrid = []
       Qtots = []
       Qths = []
       print 'Total Charge Calculation: rmin=', rmin, 'rmax=', rmax
       for Z0 in [0.1, 0.2, 0.3, 0.4, 0.5]:
           print 'Calculating for Z0=',Z0 
           Ugrid.append(Z0)
           Qtheory = (( 1.0 / np.sqrt(r_0**2+rmin**2))
                      - (1.0/ np.sqrt(r_0**2 + rmax**2)))
           Qtheory *= Z0 * np.pi / 8.0
           Us = Z0 / np.sqrt(r**2 + r_0**2)
           rhotot = graphene.rho_U(Us)
           
           import pylab as pl
           pl.title('Total density for Z0=%g'%Z0)
           pl.plot(r,rhotot, label='total charge density')
           pl.legend()
           pl.show()
           pl.figure()

           Qsim = 0.5 * dr * rhotot[imin] * r[imin]
           Qsim += 0.5 * dr * rhotot[imax] * r[imax]
           for i in range (imin+1, imax):
               Qsim += dr * rhotot[i] * r[i]
           Qsim *= -2.0 * np.pi
           Qtots.append(Qsim)
           Qths.append(Qtheory)
       Ugrid = np.array(Ugrid)
       Q_sim = np.array(Qtots)
       Q_linear = np.array(Qths)
       Q_bs = (np.pi/8.0*Ugrid + 0.19*Ugrid**3)

   import pylab as pl
   pl.plot(graphene.r, graphene.rho_0, label='rho_0')
   pl.plot(graphene.r,          rho_1, label='rho_0 + U')
   pl.plot(graphene.r, graphene.rho_0 + rho_U, label='sum')
   pl.legend()
   pl.figure()
   pl.plot(r, rho_RPA, label='RPA response (full)')
   pl.plot(r, rho_RPAb, label='RPA response (band)')
   pl.plot(r, rho_U,   label='response from sim')
   pl.plot(r, rho_Ub,   label='response from sim (band)')
   pl.plot(r, rho_th,  label='expected')
   pl.legend()
   pl.figure()
   pl.loglog(r, np.abs(rho_RPA), label='RPA response (full)')
   pl.loglog(r, np.abs(rho_RPAb), label='RPA response (band)')
   pl.loglog(r, np.abs(rho_U),   label='response from sim')
   pl.loglog(r, np.abs(rho_Ub),   label='response from sim (band)')
   pl.loglog(r, np.abs(rho_th),  label='expected')
   pl.legend()
   pl.figure()
   pl.plot(Ugrid, Q_sim, label='sim')
   pl.plot(Ugrid, Q_linear, label='Linear')
   pl.plot(Ugrid, Q_bs, label='B-S')
   pl.legend()
   pl.show()
