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
    N0 = 1000
    rtest = 1.0
    threshold = 1e-7
    r = data['r']
    U = data['U']
    itest = (r-rtest).argmin()
    Mmax = 10
    mlist = np.arange(0.0, Mmax+0.5, 1.0)

    Emax = Emax0
    Emin = Emin0

    Espace0 = np.linspace(Emin,Emax,N0)
    Espace = np.array(Espace0)
    dos0 = odedirac.doscalc(Espace0, r, U, mlist)
    dos = np.array(dos0)
    ipeak = dos0[itest,:].argmax()

    while (Emax - Emin)/N0 > threshold:
        if ipeak > 0:
            Emin = Espace[ipeak - 1]
        else: Emin = Espace[0]
        if ipeak < len(Espace) - 1:
            Emax = Espace[ipeak + 1]
        else: Emax = Espace[-1]
        Espace = np.linspace(Emin,Emax,N0)
        dos = odedirac.doscalc(Espace, r, U, mlist)
        ipeak = dos[itest,:].argmax()

    Epeak = Espace[ipeak]
    print "Peak position =", Epeak

    chi_u, chi_d = odedirac.odepsi_m(Espace,r,U,0)

    def regions(E):
        return (E - U)**2 - (0.0 + 0.5)**2 / r**2

    figure()
    title('Wavefunctions for state E=%g' %Epeak)
    plot(r,chi_u[:,ipeak], label='chi up')
    plot(r, chi_d[:,ipeak], label='chi down')
    legend()

    figure()
    plot(r, U, label='Potential')
    plot(r,dos[:,ipeak], label='peak DOS(E=%g)' %Epeak)
    plot(r,dos0[:,0], label='DOS(E=%g)' %Espace0[0])
    plot(r,dos0[:,-1], label='DOS(E=%g)' %Espace0[-1])
    plot(r,1e3*regions(Epeak), '--', label='SCALED Pr^2')
    axis([r[0], r[-1], 1.1*U[0], 1.1*dos[0,ipeak]])
    legend()

    if True:
        figure()
        plot(Espace0, dos0[itest, :], label='peak')
        legend()
    
    show()
