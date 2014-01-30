import numpy as np
import pylab as pl
from scipy import special
import math

def psi1(z): 
    h = 1e-5
    z1 = z + h; 
    z2 = z - h
    return (special.psi(z1) - special.psi(z2)) / (z1 - z2); 

def Lambda_m(m, alpha):
    kappa = m + 0.5
    gam = math.sqrt(kappa**2 - alpha**2)
    z = gam - 1j * alpha
    gl = special.gammaln(z)
    c1 = gl.imag
    c2 = - 0.5 * math.atan(alpha / gam) 
    p = special.psi(z) * z
    c3 = - p.imag
    c4 = alpha / kappa / 2; 
    c5 = - alpha * kappa * psi1(kappa)
    return (c1 + c2 + c3 + c4 + c5) * 2.0 / math.pi

def Q(alpha):
    print alpha
    s = 0.0
    for m in range(100):
        s += Lambda_m(m, alpha)
    return math.pi / 8.0 * alpha + s

alphas = np.arange(-0.499, 0.50, 0.001)
qs = np.vectorize(Q)(alphas)
q_rpa = math.pi / 8.0 * alphas
q_bs = q_rpa + 0.19 * alphas**3 
pl.plot (alphas, qs, label='M-T')
pl.plot (alphas, q_rpa, label='RPA')
pl.plot (alphas, q_bs, label='B-S')
#alphapos = 
pl.legend()
f = open("milstein.dat", 'w')
for i in range(len(alphas)):
    f.write("%g\t%g\t%g\t%g\n" % (alphas[i], qs[i], q_rpa[i], q_bs[i]))
pl.show()
