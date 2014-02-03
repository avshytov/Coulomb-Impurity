import numpy as np
import math
from scipy import special
from scipy import linalg
from scipy import integrate
import RPA

#
# This module is supposed to calculate Coulomb kernel, which 
# allows one to find the potential U(r) created by a centrally
# symmetric charge density rho(r) as
# 
#    U = np.dot(Q, rho)
#

def charge(rho, r):
    s = 0.0
    for i in range (len(r) - 1):
        dr = r[i + 1] - r[i]
        s += (rho[i] * r[i] + rho[i + 1] * r[i + 1])/2.0 * dr   
    return s * 2.0 * math.pi


#
#  The main routine. 
#    Input:
#       r -- the grid on which rho(r) and U(r) are specified
#   
#       Returns: kernel matrix
#
def do_kernel(r):  
    N = len(r)
    K = np.zeros ((N, N))
    
    for i in range (0, N):
        ri = r[i]
        if (i > 0):
            dri = r[i] - r[i - 1]
        else:
            dri = r[i + 1] - r[i]
        for j in range (i + 1, N):
            rj  = r[j]
            drj = r[j] - r[j - 1]
            mij = 4.0 * ri * rj / (ri + rj)**2
            kval = special.ellipk(mij)  # elliptic integral
            K[i, j] = kval / (ri + rj) * 4.0 * drj 
            K[j, i] = kval / (ri + rj) * 4.0 * dri 
        if i != 0:
           # contributions of left and right neighbours
           K[i, i] = 2.0/ri * (math.log(8.0*ri/dri) + 1.0) * dri
        else:
           K[i, i]  = 1.0 / ri * (math.log(8.0*ri/dri) + 1.0) * dri
           K[i, i] += math.log(8.0 ) + 1.0 # interval [0, r[0]] is 
                                           # to be treated separately
        K[i, -1] *= 0.5 # 0.5 factor from standard integration rule
	#K[i, 0]  *= 0.5
	def f(r):
	    m = 4 * ri * r / (r + ri)**2
	    return special.ellipk(m) / (r + ri) / r
	I, eps = integrate.quad(f, r[-1], np.inf, limit=100)
	K[i, -1] += I * 4 * r[-1]  # Right endpoint correction
    
    return K  
        
def kernel( r ):
    fname = "data/coulomb-kernel-Rmin=%g-Rmax=%g-N=%g.npz" % (r.min(), r.max(), len(r))
    try:
        data = np.load(fname)
        print "Loading Coulomb kernel from", fname
        print "Check r"
        assert linalg.norm(r - data['r'])<1e-8, "r vectors match"
        return data['C']
    except:
        import traceback
        traceback.print_exc()
        print "cannot load", fname, ": recalculating"
        C = do_kernel(r) * r
        np.savez(fname, C=C, r=r)
        return C


