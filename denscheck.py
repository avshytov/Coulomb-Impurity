import math
import numpy as np
import scipy
from scipy import special
from scipy import integrate
import matplotlib
from matplotlib import pylab
from pylab import *

def polaris(r, Elims,Ustr,r0):
    nr = len(r)
    func = np.zeros(nr)

    def f(q):
        kf1 = Elims[1]
        kf0 = Elims[0]
        pimin_t = np.pi / 8 * q
        pimin_b = np.pi / 8 * q
        if q <=  2 * kf1:
            piplus_t = kf1 - pimin_t
        elif q > 2 * kf1:
            piplus_t  = kf1 - kf1/2*np.sqrt(1-(4*kf1**2 / q**2))
            piplus_t += -q/4 * np.arcsin(2*kf1 / q)
        if q <=  2 * kf0:
            piplus_b = kf0 - pimin_b
        elif q > 2 * kf0:
            piplus_b  = kf0 - kf0/2*np.sqrt(1-(4*kf0**2 / q**2))
            piplus_b += -q/4 *np.arcsin(2*kf0 / q)
        dPi = (piplus_t + pimin_t) - (pimin_b + piplus_b)        
        U = (2*np.pi) * np.exp(-q*r0) / q * Ustr
        J = special.jn(0,q*ri)
        f = q / (2*np.pi)**2 * J * U * dPi
        return f
    for a in range(0,nr):
        ri = r[a]
        func[a], eps = integrate.quad(f,0,Inf)
    plot(r, func, label='polarisation operator')
    return func


def drho(r,Elims,Ustr,r0):
    N = len(r)
    rho0 = np.load('charge-density-U=0-m=10-N=%d.npy' %N)
    Emax = Elims[1]
    Emin = Elims[0]
    U0 = -1 /np.sqrt(r**2 + r0**2)
    U = Ustr * U0
    rho = np.load('charge-density-U=%g-m=10-N=%d.npy' %(Ustr,N))
    difft = - (Emax - Emin) * U / (2.0 * np.pi)
    diffs = rho - rho0
#    figure()
    plot(r, diffs, label='delta rho - sim')
    plot(r, difft, label='delta rho - theory')
    title('U =- %g/r' %Ustr)
    ylim(0,0.1)
    legend()
#    figure()
#    loglog(r, diffs, label='delta rho - sim')
#    loglog(r, difft, label='delta rho - theory')
#    title('U =- %g/r Log-Log' %a)
#    legend()
    return diffs



if __name__ == '__main__':
    r = np.load('rvec.npy')
    r0 = r[10]
    Elims = np.load('Elims.npy')
    print 'Elims =', Elims
    Ustr = [0.001, 0.01, 0.025, 0.05, 0.075, 0.10, 0.2]
    ratios = np.zeros((len(r),len(Ustr)))
    for a in range (0,len(Ustr)):
        U = Ustr[a]
        figure()
        diffp = polaris(r, Elims, U, r0)
        diffs = drho(r, Elims, U, r0)
        ratios[:,a] = diffp / diffs
    figure()
    for b in range (0,len(Ustr)):
        plot(r, ratios[:,b], label='U = %g' %Ustr[b])
        title('Ratio of Polarisation to Simulation')
        legend()
    show()
