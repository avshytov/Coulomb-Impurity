import numpy as np
from numpy import linalg
import math
import coulomb
import scipy
from scipy.interpolate import splev, splrep
from pylab import *
import mkrpa2
import util

def cheb_perm(r):
    """
        Construct the permutation which determines the 
        ordering of relaxation times tau.
        This is done recursively
    """
    if (r <= 0): return [ 0 ]
    if (r == 1):
        return [0, 1]
    prev = cheb_perm(r - 1)
    out = []
    nr = 2 * len(prev)
    for j in range(len(prev)):
        out.append(prev[j])
        out.append(nr - prev[j] - 1)
    return out;     

def cheb_tau(r, lam_max, lam_min):
    """ 
        Chebyshev set of relaxation parameters
        for a spectrum limited by eigenvalues lam_max and lam_min
    """
    i = 2**r; 
    tau_ch = np.zeros((i))
    av = (lam_max + lam_min)/2.0
    dif = (lam_max - lam_min) / 2.0
    for j in range(0, i):
        phi = (2. * j + 1.) * math.pi / 2.0 / float(i)
        tau_ch[j] = 1.0/(av + dif * math.cos(phi))
    return tau_ch; 

def make_tau_set(r, tau_max, tau_min):
    """
       Generate the set of relaxation parameters
    """
    #perm = [1, 8, 4, 5, 2, 7, 3, 6]
    #tau_max = 0.3
    #tau_min = 0.003; 
    taus = np.zeros((2**r))
    for i in range(len(taus)):
        taus[i] = tau_max * math.exp(math.log(tau_min/tau_max) * float(i)/(len(taus) - 1)) 
    #taus = [0.3, 0.2, 0.1, 0.05, 0.03, 0.01, 0.005, 0.003]
    #taus = cheb(3, lam_max, lam_min)
    perm = cheb_perm(r)
    print "perm: ", perm
    out = []
    for i in range(len(taus)):
        out.append(taus[perm[i] - 1])
    return out; 
    
def empty_callback(it, r, U, rho, U1, rho1):
    pass

def pylab_callback(it, r, U, rho, U1, rho1):
    import pylab as pl
    j = 0
    if it % 5 == 0:
                j+=1
                if j < 8:
                    pl.figure(0)
                    pl.title('Potential')
                    pl.plot(r,U, label='it=%d' %it)
                    pl.legend()
                    pl.figure(1)
                    pl.title('Rho')
                    pl.plot(r,rho, label='it=%d' %it)
                    pl.legend()
                else:
                    pl.figure(0)
                    pl.title('Potential')
                    pl.plot(r,U,'--', label='it=%d' %it)
                    pl.legend()
                    pl.figure(1)
                    pl.title('Rho')
                    pl.plot(r,rho,'--', label='it=%d' %it)
                    pl.legend()
                if j == 13:
                    j = 0
                    pl.figure(0)
                    pl.plot(r,U0, label='U0')
                    pl.legend()
                    #pl.figure(2)
                    #pl.plot(Uerror, label='U error')
                    #pl.plot(rhoerror, label='rho error')
                    #pl.legend()
                    pl.show()
                    
def solve_coulomb(rho_U, Uext, r, tau_u_set, tau_rho_set, **kwarg):
           
    params = {
       'it'               : 0,
       'zero_endpoint'    : True,
       'display_callback' : empty_callback, 
       'fname_template'   : "data/coulombsim-it=%d.npz"   
    } 
    
    params.update(kwarg)
    print params
    
    zero_endpoint = params['zero_endpoint']
    it            = params['it']
    display_callback = params['display_callback']
    fname_template = params['fname_template']
    
    rexp = util.make_exp_grid(r.min(), r.max(), len(r))
    C = coulomb.kernel( rexp )
    U = np.array( Uext )
    rho = np.zeros( (len(r)) )
    it = 0
    j = 0
    zero_U = True
    
    def mkFilename(it):
       fname = fname_template % (it)
       print "fname: ", fname
       return fname
    
    if it > 0: # Load previous solution, if possible
       data = np.load( mkFilename( it ) )
       rho = data['rho']
       U = data['U']
       
    while True:
        
        tau_U   = tau_u_set  [it % len(tau_u_set)]
        tau_rho = tau_rho_set[it % len(tau_rho_set)]
        
        it += 1
        
        rho = util.gridswap(r, rexp, rho)
        U1 = np.dot( C, rho )
        U1 = util.gridswap(rexp,  r, U1)
        rho = util.gridswap(rexp, r, rho)
        U1 += Uext
        
        if zero_endpoint: U1 -= U1[-1];
        
        err_U = linalg.norm(U - U1)
        U += tau_U * (U1 - U)
        if zero_endpoint: U -= U[-1];
        
        rho1 = rho_U(U, r)
        err_rho = linalg.norm(rho1 - rho)
        rho += tau_rho * (rho1 - rho)
        
        print "U error", err_U, "rho error", err_rho, "it", it
        np.savez(mkFilename( it ), rho=rho, U=U, rho1=rho1, U1=U1, r=r)
        
        
        display_callback(it, r, U, rho, U1, rho1)
        
        if (err_U < (1e-7)) and (err_rho < (1e-7)):
            break

    return U, rho


if __name__ == '__main__':
   Nf = 4.0 ###
   alpha = 2.5
   rmin = 0.01
   rmax = 100.0
   N = 1000
   Z = 1.0
   r_0 = 1.0
   
   import util
   r = util.make_exp_grid(rmin, rmax, N)
   #r = np.array(rexp)
   #print rexp
   print r
   #import sys
   #sys.exit(1)
   
   tau_u_set = make_tau_set(8, 0.1, 0.01)
   tau_rho_set = list(tau_u_set)
   
   def rho_minus12(U, r):
       return -U / r / 2.0 / np.pi**2 * Nf * alpha

   Uext = Z / np.sqrt (r**2 + r_0**2)

   U, rho = solve_coulomb(rho_minus12, Uext, r, tau_u_set, tau_rho_set) 
   
   ra = 5.0
   rb = 15.0

   ivals = [t[0] for t in enumerate(r) if t[1] < rb and t[1] > ra]
   xvals = [r[t] for t in ivals]
   yvals = [(abs(U[t])) for t in ivals]

   grad = np.polyfit(log(xvals), log(yvals), 1)
   print grad
   fitvals = exp(grad[1] + grad[0]*log(xvals))

#figure()
#plot (r, U)
#figure()
#loglog (r, abs(U), label='U')
#loglog (r, 1.0/r, label='1/r')
#loglog(r, 1.0/r/r, label='1/r^2')
#loglog(xvals, abs(fitvals), label='Fit: alpha = %g' % (-grad[0]))
#legend()
#show()
