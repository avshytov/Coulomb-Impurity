import sys
import numpy as np
from pylab import *
import odedirac

#
#   A script to analyse DOS peak
#   python dos-peak.py file.npz Emin Emax
#   

if __name__ == '__main__':
    filename = sys.argv[1]
    Emin = float(sys.argv[2])
    Emax = float(sys.argv[3])
    data=np.load(filename)
    r = data['r']
    U = data['U']
    
    Mmax = 10
    Espace = np.linspace(Emin,Emax,200)
    mlist = np.arange(0.0, Mmax+0.5, 1.0)
    dos = odedirac.doscalc(Espace, r, U, mlist)

    figure()
    for ri in [0.1, 1.0, 2.0, 3.0, 5.0, 10.0]:
        i = np.abs(r - ri).argmin()
        plot(Espace, dos[i, :], label='r = %g' % r[i])
    legend()
    
    show()
