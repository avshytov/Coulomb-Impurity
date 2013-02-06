import numpy as np
import math
import matplotlib
from matplotlib import pylab
from pylab import *


def ldoscalc():
    dostens = np.load("dostens.npy")
    cdtens = np.load("cdtens.npy")
    r = np.load("rvec.npy")
    m0 = len(cdtens)
    r0 = len(cdtens[0])
    i0 = len(dostens[0,0])
    n0 = len(dostens[0])
    cdmat = np.zeros((r0,n0))
    dosmat = np.zeros((n0,i0))
    ldosmat = np.zeros((r0,i0))
    for m in range (0, m0):
        print "Calculating ldos for:", (m+1), "/", m0
        cdmat = cdtens[m,:,:]
        dosmat = dostens[m,:,:]
        ldosmat += np.dot(cdmat, dosmat)
    ldosmat = np.transpose(ldosmat)
    ldosmat = ldosmat[:,:] / r
    np.save("ldosmat",ldosmat)
    return ldosmat

def ldosplot(ldosmat):
    r = np.load("rvec.npy")
    nldosmat = 2 * np.pi * r * ldosmat
    E = np.load("Evec.npy")
    ldostot = np.zeros(len(E))
    print "Producing colour plot..."
    if False:
        pcolor(r,E,ldosmat, vmin=0.0, vmax=0.15)
        colorbar()
        ylim(-10.0, 10.0)
        xlim(0.0, 12.5)
        show()
        figure()

    for r1 in range (0, len(r)):
        ldostot[:] += nldosmat[:, r1]
    figure()
    plot(E, ldostot, label='Global Density of States')
    legend()
    show()
    rq = 0.2
    rp = 0.4
    ivals = [t[0] for t in enumerate(E) if t[1] > rq and t[1] < rp]
    xvals = [E[t] for t in ivals]
    yvals = [ldostot[t] for t in ivals]
    grad1 = np.polyfit(xvals, yvals, 1)
    print "gradient of ldost", grad1[0]




    for rs in [1.0, 2.0, 3.0, 4.0]:
        si = int (rs / 25.0 * 500)
        ri = r[si]
        print rs, si, ri
        
        plot (E, ldosmat[:, si], label='r = %g' % ri) 
        ra = 0.5
        rb = 1.5
        ivals = [t[0] for t in enumerate(E) if t[1] < rb and t[1] > ra]
        xvals = [E[t] for t in ivals]
        yvals = [ldosmat[t, si] for t in ivals]
        plot (xvals, yvals, label='slice %f' %si)
        grad = np.polyfit(xvals, yvals, 1)
        print grad[0], "for r =", ri 
    xlim(-2.5, 2.5)
    legend()
    show()


    return 0

if __name__ == '__main__':
    ldosmat =  ldoscalc()
    ldosplot(ldosmat)
    
