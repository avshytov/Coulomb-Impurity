import numpy as np
import math
from scipy import interpolate

def make_exp_grid(rmin, rmax, N):
    xi = np.linspace(0.0, math.log(rmax/rmin), N)
    return rmin * np.exp(xi)

def make_lin_grid(rmin, rmax, N):
    return np.linspace(rmin, rmax, N)

def gridswap(r1,r2,f1):
    spl = interpolate.splrep(r1,f1)
    f2  = interpolate.splev(r2,spl)
    return f2
            
def makeFileTemplate(s):
    return s + "-it=%d.npz"
