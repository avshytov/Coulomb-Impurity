import numpy as np
import math
import scipy
from scipy import special
import matplotlib
from matplotlib import pylab
from pylab import *


def ldoscalc():
#    dostens = np.load("dostens.npy")
    r = np.load("rvec.npy")
    info = np.load("EMinfo.npy")
    mlist = np.load("mlist.npy")
    nm = len(mlist)
    nr = len(r)
    dostens = np.load("dostens.npy" %(info[0], info[1], nm, nr))
    cdtens = np.load("cdtens-U=%g-B=%g-ms=%d-N=%d.npy" %(info[0], info[1],nm,nr))
    ni = len(dostens[0,0])
    nn = len(dostens[0])
    cdmat = np.zeros((nr,nn))
    dosmat = np.zeros((nn,ni))
    ldosmat = np.zeros((nr,ni))
    for m in range (0, nm):
        print "Calculating ldos for:", (m+1), "/", nm
        cdmat = cdtens[m,:,:]
        dosmat = dostens[m,:,:]
        ldosmat += np.dot(cdmat, dosmat)
    ldosmat = np.transpose(ldosmat)
    ldosmat = ldosmat[:,:] 
    fname = "ldosmat-U=%g-m=%d-B=%g-grid=%d" %(info[0], nm, info[1], nr)
    print "Saving the data into", fname
    np.save(fname, ldosmat) 
    return ldosmat

def ldosplot(ldosmat):
    info = np.load("EMinfo.npy")
    mlist = np.load("mlist.npy")
    cdtens = np.load("cdtens.npy")
    r = np.load("rvec.npy")
    dr = np.load("drvec.npy")
    pot = np.load("potvec.npy")
    nm = len(mlist)
    ldos0 = np.load("ldosmat-U=0-m=%d-B=0-grid=200.npy" %nm)  
    nldosmat = 2 * np.pi * dr * r * ldosmat
    E = np.load("Evec.npy")
    ldostot = np.zeros((len(E)))
    nu = np.zeros((len(E)))
    nu_m = np.zeros((len(E)))
    rmax = r[len(r)-1]
    samps = np.arange(rmax -1) + 1
    cd = np.zeros((len(r)))
        
    if False:
        for r1 in range (0, len(r)):
        #Plots GDOS
            ldostot[:] += nldosmat[:, r1]
        figure()
        plot(E, ldostot, label='Global Density of States')
        legend()
        rq = 0.1
        rp = 0.5
        ivals = [t[0] for t in enumerate(E) if t[1] > rq and t[1] < rp]
        xvals = [E[t] for t in ivals]
        yvals = [ldostot[t] for t in ivals]
        grad1 = np.polyfit(xvals, yvals, 1)
        print "gradient of ldost", grad1[0]
        show()

    if True:
        ### Plots charge density
        Elims = np.zeros((2))
        Emin = 0.0
        Emax = 1.0
        Elims[0] = Emin
        Elims[1] = Emax
#        test = 1.0 / 4.0 / np.pi * (Emax**2 - Emin**2)
#        Etest = vectorize(lambda x: test)(r)
        e1min = (np.abs(E-Emin)).argmin()
        e1max = (np.abs(E-Emax)).argmin()
        #for e1 in range (e1min, e1max):
        #    cd[:] += ldosmat[e1,:] * (E[e1 + 1] - E[e1])
        figure()
        plot (E[e1min-40:e1max+40], ldosmat[e1min-40:e1max+40, 0], label='r=r[0]')
        plot (E[e1min-40:e1max+40], ldosmat[e1min-40:e1max+40, 1], label='r=r[1]')
        plot (E[e1min-40:e1max+40], ldosmat[e1min-40:e1max+40, 2], 'o', label='r=r[2]')
        plot (E[e1min-40:e1max+40], ldosmat[e1min-40:e1max+40, 3], label='r=r[3]')
        plot (E[e1min-40:e1max+40], ldosmat[e1min-40:e1max+40, 5], label='r=r[5]')
        plot (E[e1min-40:e1max+40], ldosmat[e1min-40:e1max+40, 10], label='r=r[10]')
        plot (E[e1min-40:e1max+40], ldosmat[e1min-40:e1max+40, 20], label='r=r[20]')
        title( "m=%d" %mlist[0])
        legend()
        show()
        for e1 in range (e1min+1, e1max-1):
            cd[:] += ldosmat[e1,:] * (E[e1 + 1] - E[e1])
        cd[:] += ldosmat[e1min,:] * (E[e1min + 1] - E[e1min]) * 0.5
        cd[:] += ldosmat[e1max - 1,:] * (E[e1max-1] - E[e1max -2]) * 0.5
            
        np.save("charge-density-U=%g-B=%g-m=%d-N=%d-E=%g-%g" %(info[0], info[1],nm, len(r), Emin, Emax), cd)
        figure()
        plot(r,  cd, label='Charge Density')
#        plot(r, Etest, '--', label='test')
        np.save('Elims', Elims)
        title('Charge Density: U=%g, B=%g, Emin=%g, Emax=%g, ms=%d' %(info[0], info[1], Emin, Emax, nm))
        legend()
        show()

    if False:
        print "Plotting slices..."
        for rs in [0.2, 0.4, 0.6, 0.8, 1.0]: #1.0, 2.0, 3.0, 4.0]: #samps:
            si = int (rs / r[(len(r)-1)] * len(r))
            ri = r[si]
           # shift = pot[si]
#            figure()
            plot (E, ldosmat[:, si], label='LDOS r = %g' % ri)
#            plot (E, ldos0[:, si], label='Free LDOS')
#            plot ((E + shift), ldos0[:,si], label='High E shift')
            if False:
                ra = 0.5
                rb = 1.5
                ivals = [t[0] for t in enumerate(E) if t[1] < rb and t[1] > ra]
                xvals = [E[t] for t in ivals]
                yvals = [ldosmat[t, si] for t in ivals]
                plot (xvals, yvals, label='slice %f' %si)
                grad = np.polyfit(xvals, yvals, 1)
                print grad[0], "for r =", ri 
            xlim(-6,6)
            legend()
            title('LDOS, r=%g, U=%g, B=%g' %(ri, info[0], info[1]))
          #  figure() ##!!! Only required for Magcheck loop 10/3/13 
        show()

    if False:
        ### Plots slices and calculates gradient
        print "Calculating analytical LDOS"
        for rs in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.2]: #samps:
            ldosmat = np.load('ldosmat-U=0-m=1-B=0-grid=%d.npy' %len(r))
            #figure()
            si = int (rs / r[(len(r)-1)] * len(r))
            ri = r[si]
            kr = abs(E) * ri
            nu = zeros ((shape(E)))
            for m in range (0,nm):
                jm = special.jn(abs(mlist[m]),kr)
                jm1 = special.jn(abs(mlist[m]+1),kr)
                nu_m = abs(E) / 2.0 / np.pi * (jm**2 + jm1**2)
                nu += nu_m
            plot (E, ldosmat[:, si], label='Simulated r = %g' % ri)
            plot(E, nu) 
           # ra = 0.1
           # rb = 0.5
           # ivals = [t[0] for t in enumerate(E) if t[1] < rb and t[1] > ra]
           # xvals = [E[t] for t in ivals]
           # yvals = [nu[t] for t in ivals]
           # grad = np.polyfit(xvals, yvals, 1)
           # print grad[0], "for theory r =", ri
            xlim(-2.5, 2.5)
            legend()
        show()

    if False:
        ###Colour plot
        print "Producing colour plot..."
        pcolor(r,E,ldosmat, vmin=0.0, vmax=0.5)
        colorbar()
        ylim(-10.0, 10.0)
        xlim(0.0, 20.0)# 12.5)
        title('LDOS: U=%g, B=%g' %(info[0], info[1]))
        show()
        figure()

    return 0

if __name__ == '__main__':
    ldosmat =  ldoscalc()
    ldosplot(ldosmat)
    
