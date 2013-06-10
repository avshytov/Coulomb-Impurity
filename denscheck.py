import math
import numpy as np
import scipy
from scipy import special
from scipy import integrate
import matplotlib
from matplotlib import pylab
from pylab import *


def Pi_mono(kf,q):
    #if abs(kf) < 1e-2: return q / 16.0; 
    pimin = np.pi * q / 8.0
    kfabs = abs(kf)
    if abs(q) <= 2.0 * kfabs:
            piplus = kfabs - pimin
    else:
            piplus  = kfabs - kfabs / 2.0 * np.sqrt( 1 - ((2.0*kfabs / q)**2))
            piplus += -q / 4.0 * np.arcsin( 2.0 * kfabs / q)
    return (pimin + piplus) / 2.0 / np.pi

if False: 
   qvals   = np.linspace(0.0, 10.0, 1000)
   figure()
   for kf in [1.0, 2.0]: 
       pivals = np.vectorize(lambda q: Pi_mono(kf, q)) (qvals)
       plot (qvals, pivals, label='kF = %g' % kf)
   plot (qvals, 1.0 /16.0 * qvals, 'k--', label='kf = 0 limit')
   legend()
   show() 


def polaris(r, Elims,Ustr,r0):
    nr = len(r)
    func = np.zeros(nr)

    for a in range(0,nr):
        def f(q):
            Pi_1 = Pi_mono( Elims[1], q )
            Pi_0 = Pi_mono( Elims[0], q )        
            dPi = Pi_1 - Pi_0
            U = (2 * np.pi) * np.exp( - q * r0 ) / q * Ustr
            J = special.jn( 0, q * r[a] )
            return  q * J * U * dPi
        func[a], eps = integrate.quad(f,0,Inf)
    func = func / (2*np.pi)
    plot(r, func, label='polarisation operator')
    return func

def drho(r,Elims,Ustr,r0):
    N = len(r)    
    check = np.zeros((len(r)))
    Emax = Elims[1]
    Emin = Elims[0]
    drho = np.zeros((N))
    nm = 10
    U0 = -1 /np.sqrt(r**2 + r0**2)
    U = Ustr * U0
    mdrho = np.load("mdrhos-cou-U=%g-B=0-m=%d-grid=%d-E=%g-%g.npy" 
                   %(Ustr, nm, N, Emin, Emax))
    for m in range (0,nm):
         drho[:] += mdrho[:,m,3] 
    difft = - (abs(Emax) - abs(Emin)) * U / (2.0 * np.pi)
    plot(r, drho, label='delta rho - sim')
    plot(r, difft, 'r--', label='delta rho - theory')
    title('Change in charge density; U = %g, Energy window: %g to %g'  %(-Ustr, Emin, Emax))
#    ylim(0,0.1)
    legend()
    show()
    figure()
    loglog(r, r**(-2.0), 'k--', label='r^-2')
    loglog(r, r**(-1.0), 'k--', label='r^-1')
    loglog(r, drho, label='delta rho - sim')
    loglog(r, difft, label='delta rho - theory')
    title('U =- %g/r Log-Log' %Ustr)
#    legend()
#    show()
    return drho


if __name__ == '__main__':
    r = np.load('rvec.npy')
    r0 = r[16]
    Elims = np.load('Elims.npy')
    print 'Elims =', Elims
    Ustr = [0.001, -0.001, 1.3, -1.3, 0.8, -0.8] 
    ratios = np.zeros((len(r),len(Ustr)))
    for a in range (0,len(Ustr)):
        U =  Ustr[a]
        print "Checking U = ", -U, "/r"
        figure()
        diffp = polaris(r, Elims, U, r0)
        diffs = drho(r, Elims, U, r0)
        ratios[:,a] = diffs / diffp
        loglog(r, diffp, label='Polarisation Operator')
        legend()
        show()
    figure()
   # np.save('ratios-N=%d-m=1-E=%g-%g' %(len(r), Elims[0], Elims[1]), ratios)
    plot([r[0],r[-1]], [1,1], 'k--')
    for b in range (0,len(Ustr)):
        plot(r, ratios[:,b], label='U = %g' %Ustr[b])
        title('Ratio of Simulation / Polarisation')
        legend()
    show()
