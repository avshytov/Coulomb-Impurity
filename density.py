import numpy as np
import util
from numpy import linalg
import math
import coulomb
import scipy
from scipy import special
from scipy import interpolate
import mkrpa2
import rpakernel
import rpam
import rpacorr
import deltafunc
import odedirac
from odedirac import (_sgn) 

class Kernel:
    def __init__ (self, r, Q):
        self.r = r
        self.Q = Q
        
    def __call__ (self, r, U):
        Us = interpolate.splrep(r, U)
        Un = interpolate.splev(self.r, Us)
        Xn = np.dot(self.Q, Un)
        Xs = interpolate.splrep(self.r, Xn)
        return interpolate.splev(r, Xs, der=0)

class RPA_tot: 
    def __init__ (self, r_inter, r_intra, kF):
       Q_inter = mkrpa2.RPA_inter(r_inter.r, r_inter.label)
       self.Q_inter = Kernel(r_inter.r, Q_inter)
       Q_intra = np.zeros((len(r_inter.r), len(r_inter.r)))
       if (kF * r_intra.r.max() > 0.01):
           Q_intra = mkrpa2.RPA_intra(r_intra.r, kF, r_intra.label)
       self.Q_intra = Kernel(r_intra.r, Q_intra)
       
    def __call__ (self, r, U): 
        return self.Q_inter(r, U) + self.Q_intra(r, U)
    
class RPA_m:
    def __init__ (self, r_inter, r_intra, kF, mlist):
       Qm_inter = rpakernel.kernel_m_inter(r_inter.r, mlist, r_inter.label)
       self.Qm_inter = Kernel(r_inter.r, Qm_inter)
       Qm_intra = np.zeros((len(r_intra.r), len(r_intra.r)))
       if (kF * r_intra.r.max() > 0.01):
           Qm_intra = rpam.kernel_m_intra(r_intra.r, mlist, kF, '3' + r_intra.label) 
           #rpakernel.kernel_m_intra(r, mlist, kF)
       self.Qm_intra = Kernel(r_intra.r, Qm_intra)
    def __call__ (self, r, U): 
        return self.Qm_inter(r, U) + self.Qm_intra(r, U)
        
class Correction:
    def __init__ (self, r, Q):
        self.r = r
        self.Q = Q
    def __call__ (self, r, U):
        Us = interpolate.splrep(r, U)
        Un = interpolate.splev(self.r, Us)
        # Assume U(r) = U[-1] * r[-1]/r for now
        # One can improve this approximation by including 
        # further powers, (r[-1]/r)^2, etc., as well as 
        # a constant term. 
        Xn = self.Q[:, 1] * Un[-1]       
        Xs = interpolate.splrep(self.r, Xn)
        return interpolate.splev(r, Xs, der=0)
        
class RPA_corr: 
    def __init__ (self, r_inter, r_intra, kF, mlist):
        Q_inter = rpacorr.RPA_corr_inter(r_inter.r, r_inter.label)
        if (kF * max(r_intra.r) > 0.01):
            Q_intra = rpacorr.RPA_corr_intra(r_intra.r, kF, r_intra.label)
        else:
            Q_intra = np.zeros((len(r), len(Q_inter[0, :])))
        self.C_inter = Correction(r_inter.r, Q_inter)
        self.C_intra = Correction(r_intra.r, Q_intra)
        
    def __call__ (self, r, U):
        return self.C_intra(r, U) + self.C_inter(r, U)
    
    
class Grid:
    def __init__ (self, r, label):
        self.r = r
        self.label = label
    
class GrapheneResponse:
    def __init__ (self, r, Ef, **kwargs):
        params = {
           'Temp'   : 0.01, 
            'Ecut'  : -1.5, 
            'B'     : 0.0, 
            'Mmax'  : 16, 
            'grid_intra' : 'lin', 
            'grid_inter' : 'exp'
        }
        params.update(**kwargs)
        self.Emin = min(Ef, params['Ecut']) # Ecut could be both positive
        self.Emax = max(Ef, params['Ecut']) # and negative
        self.Temp = params['Temp']
        self.r    = r
        self.rexp = util.make_exp_grid(r.min(), r.max(), len(r))
        self.B = params['B']
        self.mlist =  np.array(range(0, params['Mmax']))
        
        grid_dict = {
            'lin' : self.r, 
            'exp' : self.rexp 
        }
        
        r_inter = Grid(grid_dict[params['grid_inter']], params['grid_inter'])
        r_intra = Grid(grid_dict[params['grid_intra']], params['grid_intra'])
        
        if False:  # Use this to skip kernel calculation; helpful in debugging
           self.Q_Emin = np.zeros((len (self.rexp), len(self.rexp)))#RPA_kernel(self.rexp, abs(self.Emin))
           self.Q_Emax = np.zeros(np.shape(self.Q_Emin)) #RPA_kernel(self.rexp, abs(self.Emax))
        else:
           # Full RPA kernel
           self.Q_Emin  = RPA_tot(r_inter, r_intra, abs(self.Emin)) 
           self.Q_Emax  = RPA_tot(r_inter, r_intra, abs(self.Emax)) 
           #Kernel(self.rexp, RPA_kernel(self.rexp, abs(self.Emin)))
           #self.Q_Emax  = Kernel(self.rexp, RPA_kernel(self.rexp, abs(self.Emax)))
           # m-resolved RPA kernel
        if False:
           self.Qm_Emax = self.Q_Emax; #np.zeros(np.shape(self.Q_Emin))#rpakernel.kernel_m(self.rexp, self.mlist, abs(self.Emax))
           self.Qm_Emin = self.Q_Emin; #np.zeros(np.shape(self.Q_Emax))#rpakernel.kernel_m(self.rexp, self.mlist, abs(self.Emin))
        else:
           self.Qm_Emin = RPA_m(r_inter, r_intra, abs(self.Emin), self.mlist) 
           self.Qm_Emax = RPA_m(r_inter, r_intra, abs(self.Emax), self.mlist) 
           #rpakernel.kernel_m(self.rexp, self.mlist, abs(self.Emax))
           #self.Qm_Emin = 
           #rpakernel.kernel_m(self.rexp, self.mlist, abs(self.Emin))
        self.rho_0 = self.diracDensity(np.zeros(np.shape(r)))
        self.Q_corr = RPA_corr(r_inter, r_intra, abs(self.Emax), self.mlist)
        N = len(self.r)
          
    def diracDensity(self, U):
        return odedirac.rhocalc(self.Emin, self.Emax, self.r, U, self.mlist)
                       
    def bandResponse(self, U):
        rho_b = self.diracDensity(U)
        return (rho_b - self.rho_0) 
        
    def apply_kernel(self, Q, U):
        U_exp = util.gridswap(self.r, self.rexp, U)
        rho_exp = np.dot(Q, U_exp)
        return util.gridswap(self.rexp, self.r, rho_exp)
    
    def rho_RPA(self, U):
        return self.Q_Emax(self.r, U) 
    
    def rho_RPA_m(self, U):
        return self.Qm_Emax(self.r, U) - self.Qm_Emin(self.r, U)
    
    def highm(self, U): # Calculation of high-m correction within RPA
        r = self.r
        highm_max = - (self.Q_Emax(r, U) - self.Qm_Emax(r, U))
        highm_min = - (self.Q_Emin(r, U) - self.Qm_Emin(r, U))
        return highm_max - highm_min
        
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
    
    def U2(self, E, U):
        s = np.zeros(np.shape(self.r))
        #
        # The derivative of epsilon * sum(J_m^2 (epsilon*r) + J_{m+1}^2)
        # is sum( (2m + 1) J_m^2 + (2m + 3) J_{m + 1}^2). 
        # This could be shown with the help of recursion relations 
        #
        for m in self.mlist: 
            b  = special.jn(m,     E * self.r)
            b1 = special.jn(m + 1, E * self.r)
            s += b**2  + b1**2 
        m0 = self.mlist[0]
        m1 = self.mlist[-1] + 1
        s += 2.0 * m0 * special.jn( m0, E * self.r)**2
        s -= 2.0 * m1 * special.jn( m1, E * self.r)**2
        s *= 2.0
        if (E < 0): s *= -1.0
        #import pylab as pl
        #pl.plot(r, s)
        #pl.show()
        return U * U / 8.0 / np.pi * s 
        
    
    def seaContribution(self, U):
        sgnE = 1.0
        if (self.Emin < 0): sgnE = -1.0
        rho = np.zeros(np.shape(U))
        #rho = self.highm(U) * sgnE + self.highm2(U) 
        #rho  = self.highm( self.Emax,  U)
        #rho -= self.highm( self.Emin,  U)
        
        #rho += 1.0*U**2 / 4.0 / np.pi * sgnE # Quadratic contribution
        rho += self.U2(self.Emin, U)
        
        rho += self.Qm_Emin(self.r, U)
        rho += self.Q_Emax(self.r, U) - self.Qm_Emax(self.r, U)
        return rho
        
    def rho_U(self, U):    
        rho  = self.bandResponse(U)
        rho += self.seaContribution(U) 
        rho += self.Q_corr(self.r, U) 
        #Rmax = 45.0
        #imax = np.abs(self.r - Rmax).argmin()
        #rho_rpa = self.Q_Emax(self.r, U) #self.apply_kernel(self.Q_Emax, U)
        #rho[imax:] = rho_rpa[imax:]
        return rho 

            
if __name__ == '__main__':
   rmin = 0.01
   rmax = 50.0
   N = 500
   r = util.make_lin_grid(rmin, rmax, N) 
   Ef = -0.2
   Ecut = -3.0
   Mmax = 31
   #graphene = GrapheneResponse(r, -1e-4, Ecut=Ecut, Mmax=31)
   graphene = GrapheneResponse(r, -1e-4, Ecut=Ecut, Mmax=Mmax, grid_intra='lin')
   
   Z = 0.02
   r_0 = 1.0
   #U = Z / np.sqrt(r**2 + r_0**2)
   #def Ur(rx):
   #    if rx > r_0: 
   #       return -Z/rx 
   #    return -Z / r_0
   #U = np.vectorize(Ur)(r)
   #def rho_b(rx):
   #     #if rx > r_0: return 0.0; 
   #     #return Z / math.pi / r_0**2
   #     return r_0**2 / math.pi / (rx**2 + r_0**2)**2 * Z
   delta_func = deltafunc.DeltaGauss(r_0)
   #delta_func = deltafunc.DeltaCubic(r_0)
   rho_bare = Z * delta_func.rho(r)
   U = Z * delta_func.U(r)
   
   rexp = util.make_exp_grid(0.001, 50.0, 1000)
   Qcoul = Kernel(rexp, coulomb.kernel(rexp))
   #rho_bare = np.vectorize(rho_b)(r)
   #U = Qcoul(r, rho_bare)
   #U[0:3] = U[3]
   if True:
      import pylab as pl
      pl.plot (r, U)
      pl.plot (r, Z / np.sqrt(r**2 + r_0**2))
      pl.figure()
      pl.loglog(r, np.abs(U))
      pl.loglog(r, Z / np.sqrt(r**2 + r_0**2))
      pl.show ()
   rho_th = - math.pi / 8.0 * rho_bare
   #rho_th = -Z * r_0 / 16.0 / np.sqrt(r**2 + r_0**2)**3
   
   rho_RPA = graphene.rho_RPA(U)
   rho_U   = graphene.rho_U(U)
   rho_Ub  = graphene.bandResponse(U) 
   rho_c = graphene.Q_corr(r, U)
   rho_RPAb = graphene.Q_Emax(r, U) - graphene.Q_Emin(r, U)
   rho_RPAbm = graphene.Qm_Emax(r, U) - graphene.Qm_Emin(r, U)
   #rho_RPAb = graphene.apply_kernel(graphene.Q_Emax - graphene.Q_Emin, U)
   #rho_RPAbm = graphene.apply_kernel(graphene.Qm_Emax - graphene.Qm_Emin, U)
   rho_0 = (graphene.Emax**2 - graphene.Emin**2) / 4.0 / math.pi
   rho_1 = graphene.diracDensity(U)   

   if True:
       imin = 0
       r_stop = 50.0; 
       imax = np.abs(r - r_stop).argmin()
       rmin = r[0]
       rmax = r[-1]
       dr = r[imin+1]-r[imin]
       Zvals = []
       Qtots = []
       Qths = []
       print 'Total Charge Calculation: rmin=', rmin, 'rmax=', rmax
       #for Z0 in [0.1]: 
       #for Z0 in [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
       for Z0 in np.linspace(-0.5, 0.5, 41):
       #for Z0 in [0.6, 0.8, 1.0, 1.2, 1.4]:
           print 'Calculating for Z0=',Z0 
           Zvals.append(Z0)
           Qtheory = (( r_0 / np.sqrt(r_0**2+rmin**2))
                      - (r_0/ np.sqrt(r_0**2 + rmax**2)))
           Qtheory *= Z0 * np.pi / 8.0
           #Us = Z0 / np.sqrt(r**2 + r_0**2)
           Us = Z0 / Z * U
           rhotot = graphene.rho_U(Us)
           
           rho_rpa = rho_bare * (-math.pi/8.0)
           #-Z0 / 16.0 * r_0 / np.sqrt(r**2 + r_0**2)**3
           
           
           Qsim = 0.5 * dr * rhotot[imin] * r[imin]
           Qsim += 0.5 * dr * rhotot[imax] * r[imax]
           Qrpa = 0.5 * dr * rho_rpa[imin] * r[imin]
           Qrpa += 0.5 * dr * rhotot[imax] * r[imax]
           for i in range (imin+1, imax):
               Qsim += dr * rhotot[i] * r[i]
               Qrpa += dr * rho_rpa[i] * r[i]
           Qsim *= -2.0 * np.pi
           Qrpa *= -2.0 * np.pi
           print "Z0 = ", Z0, "Qsim = ", Qsim, "Qrpa = ", Qrpa, "th: ", Qtheory
           import pylab as pl
           if False:
              pl.figure()
              pl.title('Total density for Z0=%g'%Z0)
              pl.plot(r,rhotot, label='total charge density')
              pl.plot(r,rho_rpa, label='RPA charge density')
              pl.plot(r, rho_rpa - rhotot, label='diff')
              pl.legend()
              pl.figure()
              pl.title('Total density for Z0=%g'%Z0)
              pl.plot(r,r*rhotot, label='r*total charge density')
              pl.plot(r,r*rho_rpa, label='r*RPA charge density')
              pl.plot(r, r*(rho_rpa - rhotot), label='r*diff')
              pl.legend()
              pl.show()
           f = open('rho-Z=%g.dat' % (Z0), 'w')
           for i in range(len(r)):
               f.write('%g\t%g\n' % (r[i], rhotot[i]))
           f.close()
           Qtots.append(Qsim)
           Qths.append(Qtheory)
       f = open("Q-data.dat", 'w')
       for i in range(len(Zvals)):
           f.write("%g\t%g\t%g\n" % (Zvals[i], Qtots[i], Qths[i]))
       f.close()
       Zvals = np.array(Zvals)
       Q_sim = np.array(Qtots)
       Q_linear = np.array(Qths)
       Q_bs = (np.pi/8.0*Zvals + 0.19*Zvals**3)
       pl.figure()
       pl.plot(Zvals, Q_sim, label='sim')
       pl.plot(Zvals, Q_linear, label='Linear')
       pl.plot(Zvals, Q_bs, label='B-S')
       pl.legend()

   import pylab as pl
   pl.figure()
   pl.plot(graphene.r, graphene.rho_0, label='rho_0')
   pl.plot(graphene.r,          rho_1, label='rho_0 + U')
   pl.plot(graphene.r, graphene.rho_0 + rho_U, label='sum')
   pl.title('Ecut = %g' % graphene.Emin)
   pl.legend()
   pl.figure()
   pl.plot(r, rho_RPA, label='RPA response (full)')
   pl.plot(r, rho_RPAb, label='RPA response (band)')
   pl.plot(r, rho_RPAbm, label='RPA response (band-m)')
   pl.plot(r, rho_U,   label='response from sim')
   pl.plot(r, rho_c,   label='correction from r->inf')
   pl.plot(r, rho_Ub,   label='response from sim (band)')
   pl.plot(r, rho_th,  label='expected')
   pl.plot(r, rho_Ub - rho_RPAbm, label='rho_Ub - rho_RPA_b(m)')
   pl.plot(r, rho_Ub - rho_RPAb, label='rho_Ub - rho_RPA_b')
   pl.title('Ecut = %g' % graphene.Emin)
   pl.legend()
   pl.figure()
   pl.loglog(r, np.abs(rho_RPA), label='RPA response (full)')
   pl.loglog(r, np.abs(rho_RPAb), label='RPA response (band)')
   pl.loglog(r, np.abs(rho_U),   label='response from sim')
   pl.loglog(r, np.abs(rho_c),   label='correction from r->inf')
   pl.loglog(r, np.abs(rho_Ub),   label='response from sim (band)')
   pl.loglog(r, np.abs(rho_th),  label='expected')
   pl.legend()
   pl.title('Ecut = %g' % graphene.Emin)
   pl.figure()
   pl.plot(r, rho_RPAb - rho_RPAbm, label='rho_RPA_U - rho_RPA_bm')
   pl.title("diff rpa resp Ecut=%g" % graphene.Emin)
   pl.legend()
   f = open('rho-u-Z=%g-Ecut=%g-Mmax=%d.dat' % (Z, Ecut, Mmax), 'w')
   for i in range(len(r)):
       f.write('%g\t%g\n' % (r[i], rho_U[i]))
   f.close()
   pl.show()
