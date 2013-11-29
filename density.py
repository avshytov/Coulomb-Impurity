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


#def backgroundDensity(r, mlist, Temp):
#    rho_0 = bandDensity(r, 0.0 * r, mlist, BTemp)
#    return rho_0
    
def RPA_kernel(r, kF):
    Q = mkrpa2.RPA_inter(r)  
    if (kF * r.max() > 0.01):
        Q += mkrpa2.RPA_intra(r, kF)
    return Q

class GrapheneResponse:
    def __init__ (self, r, Ef, **kwargs):
        params = {
           'Temp'   : 0.01, 
            'Ecut'  : -1.5, 
            'B'     : 0.0, 
            'Mmax'  : 10
        }
        params.update(**kwargs)
        self.Emin = params['Ecut']
        self.Emax = Ef
        self.Temp = params['Temp']
        self.r    = r
        self.rexp = util.make_exp_grid(r.min(), r.max(), len(r))
        self.B = params['B']
        self.mlist =  np.array(range(0, params['Mmax']))
        
        self.Q_Emin = RPA_kernel(self.rexp, abs(self.Emin))
        self.Q_Emax = RPA_kernel(self.rexp, abs(self.Emax))
        self.rho_0 = self.diracDensity(np.zeros(np.shape(r)))
        #### Correction
        N = len(self.r)
        self.F = 1.01340014 + 0.05991426 / np.sqrt(N) + 7.12091516 / N 
        
    def diracDensity(self, U):
        return bandDensity(self.r, U, self.mlist, self.B, 
                           self.Emin, self.Emax, self.Temp)
                           
    def bandResponse(self, U):
        rho_b = self.diracDensity(U)
        return (rho_b - self.rho_0) * self.F
        
    def apply_kernel(self, Q, U):
        Uexp = util.gridswap(self.r, self.rexp, U)
        rho_exp = np.dot(Q, Uexp)
        return util.gridswap(self.rexp, self.r, rho_exp)
    
    def rho_RPA(self, U):
        return self.apply_kernel(self.Q_Emax, U)
    
    def highm(self, E, U):
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

    
    def seaContribution(self, U):
        rho = np.zeros(np.shape(U))
        #rho  = self.highm( self.Emax,  U)
        #rho -= self.highm( self.Emin,  U)
        
        #rho += -U**2 / 4.0 / np.pi # Quadratic contribution
        
        rho += self.apply_kernel(self.Q_Emin, U) 
        return rho
        
    def rho_U(self, U):    
        rho  = self.bandResponse(U)
        rho += self.seaContribution(U) 
        #Rmax = 20.0
        #imax = np.abs(r - Rmax).argmin()
        #rho_rpa = self.apply_kernel(self.Q_Emax, U)
        #rho[imax:] = rho_rpa[imax:]
        return rho 


if __name__ == '__main__':
   rmin = 0.02
   rmax = 40.0
   N = 200
   r = util.make_lin_grid(rmin, rmax, N) 
   graphene = GrapheneResponse(r, 0.0, Ecut=-1.5)
   
   Z = 0.25
   r_0 = 1.0
   U = Z / np.sqrt(r**2 + r_0**2)
   rho_th = -Z * r_0 / 16.0 / np.sqrt(r**2 + r_0**2)**3
   
   rho_RPA = graphene.rho_RPA(U)
   rho_U   = graphene.rho_U(U)
   rho_Ub  = graphene.bandResponse(U) 
   rho_RPAb = graphene.apply_kernel(graphene.Q_Emax - graphene.Q_Emin, U)
   
   import pylab as pl
   rho_0 = (graphene.Emax**2 - graphene.Emin**2) / 4.0 / math.pi
   rho_1 = graphene.diracDensity(U)
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
   pl.show()
