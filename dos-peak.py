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
    Emin0 = float(sys.argv[2])
    Emax0 = float(sys.argv[3])
    data=np.load(filename)
    N0 = 200
    rtest = 1.0
    threshold = 1e-7
    r = data['r']
    U = data['U']
    itest = (r-rtest).argmin()
    Mmax = 10
    mlist = np.arange(0.0, Mmax+0.5, 1.0)

    Emax = Emax0
    Emin = Emin0

    while (Emax - Emin) > threshold:
        Espace = np.linspace(Emin,Emax,N0)
        dos = odedirac.doscalc(Espace, r, U, mlist)
        ipeak = dos[itest,:].argmax()
        Emin = Espace[ipeak - 1]
        Emax = Espace[ipeak + 1]

    Epeak = Espace[ipeak]
    print "Peak position =", Epeak

    figure()
    title('LDOS at Emax=%g' %Epeak)
    plot(r,dos[:,ipeak])
    show()

    if False:
        figure()
        plot(Espace, dos[itest, :], label='peak')
        legend()
        show()
